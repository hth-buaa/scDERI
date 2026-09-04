#' @title Adaptive Cell Grouping Along Pseudotime Trajectories
#'
#' @description Performs adaptive cell grouping along pseudotime trajectories for each cell type.
#' It first fits a logarithmic model to determine the minimum cells per group based on total cell count,
#' then iteratively adjusts group sizes to ensure all genes have sufficient non-zero expression/activity values in each group.
#' Returns aggregated gene profiles and pseudotime endpoints for downstream dynamic modeling.
#'
#' @param interest_cell_type_genes_pseudotime_info List. Output from \code{\link{get_genes_pseudotime_info}}.
#' @param points_x_for_fitting_nc_and_nmin Integer vector of length 3, strictly increasing positive values.
#'   Total cell counts (x-axis) for fitting the minimum-cells-per-group logarithmic model.
#' @param points_y_for_fitting_nc_and_nmin Integer vector of length 3, strictly increasing positive values.
#'   Corresponding minimum cells per group (y-axis) for the logarithmic model. Defaults to \code{c(5, 10, 15)}.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#'
#' @return Nested list named by cell type, each entry contains:
#' \itemize{
#'   \item \code{Branches_n_group}: Integer vector. Number of valid groups per branch; returns 0 if grouping fails.
#'   \item \code{Branches_n_k}: List of integer vectors. Exact cell count of each group per branch.
#'   \item \code{n_min}: Numeric vector. Fitted minimum cell threshold per branch.
#'   \item \code{Branches_T_k}: List of numeric vectors. Pseudotime value at the endpoint of each group.
#'   \item \code{Branches_Tslot_k}: List of numeric vectors. Length of pseudotime interval for each group.
#'   \item \code{Branches_TFE_T}: List of data frames. Aggregated mean of non-zero expression values for TFs (genes as rows, groups as columns).
#'   \item \code{Branches_TGA_T}: List of data frames. Aggregated mean of non-zero activity values for TGs.
#'   \item \code{Branches_TGE_T}: List of data frames. Aggregated mean of non-zero expression values for TGs.
#' }
#'
#' @details For each cell type and branch, the workflow executes as follows:
#' \enumerate{
#'   \item Validates input fitting points to ensure they are strictly increasing positive integers.
#'   \item Fits a logarithmic model \eqn{y = c \cdot \ln(ax + b)} via the three given points, where x = total cells in the branch and y = minimum cells per group. The fitted value is rounded up to obtain \code{n_min}.
#'   \item Runs adaptive grouping on pseudotime-ordered cells: each group starts with \code{n_min} cells; if any gene in TFE, TGA, or TGE matrix has fewer than \code{n_min} non-zero values, one cell is added incrementally until the constraint is met or remaining cells are exhausted.
#'   \item Handles edge cases: remaining cells fewer than \code{n_min} are merged into the previous group; if the first group cannot satisfy the constraint, grouping fails and returns 0 groups.
#'   \item Calculates endpoint pseudotime and interval length for each group.
#'   \item Computes mean of non-zero expression/activity values per gene per group, and returns aggregated matrices.
#' }
#'
#' @note
#' \itemize{
#'   \item The three fitting points directly determine \code{n_min} and thus the number of groups. Lower \code{n_min} produces more groups but increases sparsity risk.
#'   \item If grouping returns 0 groups or fewer groups than expected, you can lower the y values of fitting points to reduce \code{n_min}, or tune upstream parameters in \code{\link{get_genes_pseudotime_info}} to reduce gene sparsity.
#'   \item For reliable downstream modeling, 8–15 groups per branch is recommended. Too few groups reduce temporal resolution; too many groups weaken statistical power per group.
#'   \item The function creates a \code{2.2 Data Processing - Cell Grouping} directory, temporarily switches the working directory, and saves logarithmic model fitting plots per cell type.
#'   \item Parallel computation relies on \code{\link[parallel]{mclapply}} and is not supported on Windows.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_genes_pseudotime_info = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#' points_x_for_fitting_nc_and_nmin = c(280, 450, 650)
#' points_y_for_fitting_nc_and_nmin = c(5, 10, 15)
#' ncores = parallel::detectCores() - 1 # In Linux.
#' # ncores = 1 # In windows.
#' interest_cell_type_group = cell_grouping(
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   points_x_for_fitting_nc_and_nmin = points_x_for_fitting_nc_and_nmin,
#'   points_y_for_fitting_nc_and_nmin = points_y_for_fitting_nc_and_nmin,
#'   ncores = ncores
#' )
#' base::saveRDS(
#'   interest_cell_type_group,
#'   file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' }
cell_grouping = function(
    interest_cell_type_genes_pseudotime_info,
    points_x_for_fitting_nc_and_nmin,
    points_y_for_fitting_nc_and_nmin = c(5, 10, 15),
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Grouping cells based on pseudotime points and valid (non-zero) counts of expression and activity values. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.2 Data Processing - Cell Grouping"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Windows compatibility for parallel computing
  if (.Platform$OS.type == "windows" && ncores > 1) {
    message("Parallel computing is not supported on Windows. ncores is forced to 1.")
    ncores = 1
  }

  ### Validate input parameters
  validate_points = function(x, name = "x") {
    if (!base::is.numeric(x)) {
      base::stop("Parameter '", name, "' must be a numeric vector.")
    }
    if (base::length(x) != 3) {
      base::stop("Parameter '", name, "' must have exactly 3 elements.")
    }
    if (base::any(x <= 0)) {
      base::stop("Parameter '", name, "' must contain only positive values.")
    }
    if (base::any(x %% 1 != 0)) {
      base::stop("Parameter '", name, "' must contain only integers.")
    }
    if (!base::all(base::diff(x) > 0)) {
      base::stop("Parameter '", name, "' must be strictly increasing.")
    }
    base::message("Validation passed for parameter '", name, "'.")
  }
  message("Validating input parameter points_x_for_fitting_nc_and_nmin.")
  validate_points(points_x_for_fitting_nc_and_nmin, name = "points_x_for_fitting_nc_and_nmin")
  message("Validating input parameter points_y_for_fitting_nc_and_nmin.")
  validate_points(points_y_for_fitting_nc_and_nmin, name = "points_y_for_fitting_nc_and_nmin")

  ### Defining a Function for fitting the relationship between total cell count x and minimum cells per group y: y = c*ln(ax+b).
  ### Parameters a,b,c are solved via three fitting points (x-axis = total cells per branch, y-axis = minimum cells per group).
  calculate_n_min = function(Branches_n_cell,
                             x = c(200, 1200, 2000),
                             y = c(5, 10, 15),
                             xmax = 5000,
                             cell_type)
  {
    # Calculating the Minimum Cell Count Per Group Based on Branch Cell Counts
    base::sapply(
      base::seq_along(Branches_n_cell), function(n) {
        message("Organizing Data Point Coordinates.")
        points = base::list()
        for (i in 1:3) {points[[i]] = base::data.frame(x = x, y = y)}
        data = base::do.call(rbind, points)
        # data = base::data.frame(x = x, y = y)
        x1 = data$x[1]
        x2 = data$x[2]
        y1 = data$y[1]
        y2 = data$y[2]
        if (x1 == x2) stop("Error: The first two points cannot have identical x-values.")

        a_init = 1
        b_init = 1
        c_init = stats::median(data$y)
        message("Initializing parameters a = 1, b = 1, c = ", c_init, " for y=c*ln(ax+b).")

        message("Ensuring parameter validity.")
        if (base::any(a_init * data$x + b_init <= 0)) {
          message("Adding minor offset to ensure logarithmic domain validity.")
          b_init = b_init + 1e-6
          if (base::any(a_init * data$x + b_init <= 0)) {
            stop("Unable to find valid initial parameters. Please verify the input data.")
          }
        }

        message("Fitting using the nls.lm algorithm.")
        fit = minpack.lm::nlsLM(
          y ~ c * base::log(a * x + b),
          data = data,
          start = base::list(a = a_init, b = b_init, c = c_init),
          # start = base::list(a = 0.001, b = 1, c = 20),
          control = minpack.lm::nls.lm.control(maxiter = 1000)
        )
        a = stats::coef(fit)["a"]
        b = stats::coef(fit)["b"]
        c = stats::coef(fit)["c"]
        message("Fitting results: a = ", a, ", b = ", b, ", c = ", c, ".")

        # Plotting function: y = ceiling(c*ln(ax+b))
        message("Plotting function: y = ceiling(c*ln(ax+b)).")
        x = base::seq(100, xmax, by = 1)
        y = base::ceiling(c * base::log(a * x + b))
        df = base::data.frame(x = x, y = y)

        x_annotate = 200
        step = 300
        while((next_val = x_annotate[base::length(x_annotate)] + step) <= xmax) {
          x_annotate = c(x_annotate, next_val)
          step = step + 100
        }
        y_annotate = base::ceiling(c * base::log(a * x_annotate + b))
        annotate_df = base::data.frame(
          x = x_annotate, y = y_annotate,
          label = base::sprintf("(%.0f, %.0f)", base::ceiling(x_annotate), base::ceiling(y_annotate))
        )

        p = ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_line(color = "blue", linewidth = 1) +
          ggplot2::geom_hline(yintercept = 0, color = "gray50") +
          ggplot2::geom_vline(xintercept = 0, color = "gray50") +

          ggplot2::geom_point(data = annotate_df, ggplot2::aes(x = x, y = y),
                              color = "red", fill = NA, shape = 21, size = 3, stroke = 1.5) +

          ggrepel::geom_text_repel(data = annotate_df, ggplot2::aes(label = label),
                                   nudge_x = 50, nudge_y = 0.2 * base::max(df$y),
                                   size = 3.5, segment.color = "gray60",
                                   box.padding = 0.5, max.overlaps = Inf) +

          ggplot2::labs(
            title = "y = base::ceiling(c*log(a*x+b))",
            subtitle = base::paste0("a = ", base::format(base::round(a, 4), nsmall = 4),
                                    ", b = ", base::format(base::round(b, 4), nsmall = 4),
                                    ", c = ", base::format(base::round(c, 4), nsmall = 4)),
            x = "Number of cells in the branch",
            y = "Minimum number of cells per group in the branch") +

          ggplot2::scale_x_continuous(labels = function(x) base::format(x, big.mark = ",")) +
          ggplot2::scale_y_continuous(labels = function(y) base::format(base::round(y, 1), nsmall = 1)) +

          ggplot2::theme_minimal(base_family = "sans") +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40"),
            axis.text = ggplot2::element_text(size = 10),
            panel.grid.minor = ggplot2::element_blank()
          )

        ggplot2::ggsave(base::paste0(cell_type, " --- y = ceiling(c*ln(ax+b)).png"), width = 6, height = 5, plot = p)

        return(base::ceiling(c * base::log(a * Branches_n_cell[[n]] + b) + 0.01))
      }
    )
  }

  ### Defining a function for individual cell type processing.
  cell_grouping_process = function(cell_type,
                                   Branches_n_cell,
                                   Branches_TFE, Branches_TGA, Branches_TGE,
                                   Branches,
                                   calculate_n_min,
                                   points_x_for_fitting_nc_and_nmin, points_y_for_fitting_nc_and_nmin)
  {
    ## Partitioning Branches_n_cell[[i]] cells into Branches_n_group[[i]] groups based on k~ordering.
    ## Group k contains Branches_n_k[[i]][k] cells, corresponding to cells from position sum(Branches_n_k[[i]][1:k]) - Branches_n_k[[i]][k] + 1 to sum(Branches_n_k[[i]][1:k]) in the k~ ordered sequence.
    Branches_n_k = base::list()
    for (i in base::seq_along(Branches)) {
      message("Calculating the minimum cell count per group based on branch cell counts for branch ", i, ".")
      n_min = calculate_n_min(
        Branches_n_cell = Branches_n_cell,
        x = points_x_for_fitting_nc_and_nmin,
        y = points_y_for_fitting_nc_and_nmin,
        xmax = base::max(c(5000, base::unlist(Branches_n_cell))),
        cell_type = cell_type
      )
      message("Calculating the cell group counts per group for branch ", i, ".")
      Branches_n_k[[i]] = numeric(base::ceiling(Branches_n_cell[[i]] / n_min[i]))
    }
    message("Initializing minimum cell count thresholds for all branches: ", paste(n_min, collapse = ","), ".")
    message("Initializing cell group counts for all branches.")

    message("Preparing to execute cell grouping algorithm.")
    Branches_n_group = base::list() # Storing number of groups.
    Branches_k = base::list() # Group index k (for loop iteration)
    Branches_n_r = base::list() # Ungrouped cell count (for loop traversal)
    for (i in base::seq_along(Branches)) {
      message("Commencing cell grouping for branch ", i, ".")
      if (n_min[i] > Branches_n_cell[[i]]) {
        # If the minimum cell count threshold exceeds the total cell count, it indicates an excessively high threshold setting that prevents grouping.
        message("Minimum cell count threshold exceeds total cell count, indicating an excessively high threshold that prevents grouping.")
        Branches_n_group[[i]] = 0
        Branches_n_r[[i]] = Branches_n_cell[[i]]
        Branches_k[[i]] = 0
        Branches_n_k[[i]] = 0
        next
      }
      current_step = 1 # Algorithm execution initiated, proceeding to Step 1.
      repeat {
        if (current_step == 1) {
          # Step 1: Preparing to partition group k=1; initially there are n_r = n_cell ungrouped cells.
          Branches_k[[i]] = 1 # Preparing to partition group k=1.
          Branches_n_r[[i]] = Branches_n_cell[[i]] # Currently there are n_r = n_cell ungrouped cells (all cells remain ungrouped).
          current_step = 2 # Proceeding to Step 2.

        } else if (current_step == 2) {
          # Step 2: Assigning n_k = n_min cells to group k
          Branches_n_k[[i]][Branches_k[[i]]] = n_min[i] # Assigning n_k = n_min cells to group k
          Branches_n_r[[i]] = Branches_n_r[[i]] - Branches_n_k[[i]][Branches_k[[i]]]  # Currently there are n_r = n_r - n_k ungrouped cells remaining.
          current_step = 3 # Proceeding to Step 3.

        } else if (current_step == 3) {
          # Step 3: Obtaining counts of non-missing values for:
            # Expression vectors TF^E_{i, k~} for each TFG in the n_k cells,
            # Activity vectors TG^A_{j, k~} for each TG,
            # Expression vectors TG^E_{j, k~} for each TG
          # If any gene's non-missing values count is less than n_min, the current pre-assigned cell count for this group needs to be increased; otherwise, viable grouping cannot be achieved.
            # If Branches_n_r[[i]] == 0, it indicates no ungrouped cells remain; cell allocation cannot be increased, making this group unviable for assignment.
              # Termination condition: If k == 1 (indicating an attempt to partition the first cell group), grouping becomes impossible. The group count remains 0 and the grouping process terminates.
              # Otherwise, merge these cells into the previous group, reducing the group count k to k-1, and increase the cell count of the previous group by Branches_n_k[[i]][Branches_k[[i]]], then terminate the grouping process.
            # If Branches_n_r[[i]] > 0, increment current group cell count by 1 and return to this step for re-evaluation.
          # Otherwise, grouping is successful: k = k + 1, proceed to the next step.
          if (!all(
            all(base::rowSums(Branches_TFE[[i]][, (base::sum(Branches_n_k[[i]]) - Branches_n_k[[i]][Branches_k[[i]]] + 1) : base::sum(Branches_n_k[[i]])] != 0) >= n_min[i]),
            all(base::rowSums(Branches_TGA[[i]][, (base::sum(Branches_n_k[[i]]) - Branches_n_k[[i]][Branches_k[[i]]] + 1) : base::sum(Branches_n_k[[i]])] != 0) >= n_min[i]),
            all(base::rowSums(Branches_TGE[[i]][, (base::sum(Branches_n_k[[i]]) - Branches_n_k[[i]][Branches_k[[i]]] + 1) : base::sum(Branches_n_k[[i]])] != 0) >= n_min[i])
          )) { # If the count of non-missing values for any gene is less than n_min, it indicates that the number of cells currently pre-assigned to this group needs to be increased; otherwise, valid grouping cannot be achieved.
            if (Branches_n_r[[i]] == 0) { # No ungrouped cells remain at this point; cell allocation cannot be increased, making this group unviable for assignment.
              if (Branches_k[[i]] == 1) { # If the first cell group (Group 1)  is currently attempted to partition, grouping becomes impossible. The group count defaults to 0 and the process terminates.
                Branches_n_group[[i]] = 0
                Branches_n_r[[i]] = Branches_n_cell[[i]]
                Branches_k[[i]] = 0
                Branches_n_k[[i]] = 0
                break
              } else { # If the k-th cell group is currently attempted to partition, merge these cells into the (k-1)-th group, resulting in a total of k-1 groups, and terminate the grouping process.
                Branches_n_k[[i]][Branches_k[[i]] - 1] = Branches_n_k[[i]][Branches_k[[i]] - 1] + Branches_n_k[[i]][Branches_k[[i]]]
                Branches_k[[i]] = Branches_k[[i]] - 1
                Branches_n_r[[i]] = 0
                while (base::length(Branches_n_k[[i]]) > 0 && Branches_n_k[[i]][base::length(Branches_n_k[[i]])] == 0) {
                  Branches_n_k[[i]] = Branches_n_k[[i]][-base::length(Branches_n_k[[i]])] # Remove the trailing 0 in Branches_n_k[[i]] caused by initialization.
                }
                Branches_n_group[[i]] = Branches_k[[i]]
                Branches_n_k[[i]] = Branches_n_k[[i]][-base::length(Branches_n_k[[i]])]
                break
              }
            }
            if (Branches_n_r[[i]] > 0) { # If ungrouped cells still remain, allocate one additional cell and re-evaluate whether the non-zero count of gene expression or activity values in the current partition meets the threshold n_min.
              Branches_n_k[[i]][Branches_k[[i]]] = Branches_n_k[[i]][Branches_k[[i]]] + 1
              Branches_n_r[[i]] = Branches_n_r[[i]] - 1
              current_step = 3 # Proceeding to Step 3.
            }
            if (Branches_n_r[[i]] < 0) {
              # Edge case: This scenario should not occur during normal operation. It may result from either:
                # Omission of the initial validation check (where minimum cell threshold > total cellsindicates excessive threshold setting preventing grouping), or
                # Failure to perform the Step 5 prerequisite check (where attempting to partition a new group fails due to insufficient remaining cells).
              break
            }
          } else { # If the non-zero counts of gene expression or activity values in the currently partitioned cells all meet the threshold n_min, then the grouping of the k-th group is successful, and the process proceeds to partition the (k+1)-th group.
            if (Branches_n_r[[i]] == 0) { # Previous group allocation concurrently exhausts all available cells.
              while (base::length(Branches_n_k[[i]]) > 0 && Branches_n_k[[i]][base::length(Branches_n_k[[i]])] == 0) {
                Branches_n_k[[i]] = Branches_n_k[[i]][-base::length(Branches_n_k[[i]])] # Remove the trailing 0 in Branches_n_k[[i]] caused by initialization.
              }
              Branches_n_group[[i]] = Branches_k[[i]]
              Branches_n_r[[i]] = 0
              break
            } else {
              Branches_k[[i]] = Branches_k[[i]] + 1
              current_step = 4 # Proceeding to Step 4.
            }
          }

        } else if (current_step == 4) {
          # Step 4: Preparing to partition a new cell group (Group k).
          # If n_r < n_min, it indicates insufficient remaining cells to form a new group.
          # Merge these cells into Group k-1, resulting in k-1 total groups, and terminate grouping. Otherwise, return to Step 2.
          if (Branches_n_r[[i]] < n_min[i]) { # If n_r < n_min, it indicates insufficient remaining cells to form a new group. Merge these cells into group k-1, resulting in k-1 total groups, and terminate grouping.
            Branches_n_k[[i]][Branches_k[[i]]] = Branches_n_r[[i]]
            Branches_n_k[[i]][Branches_k[[i]]-1] = Branches_n_k[[i]][Branches_k[[i]]-1] + Branches_n_k[[i]][Branches_k[[i]]]
            Branches_k[[i]] = Branches_k[[i]] - 1
            Branches_n_r[[i]] = 0
            while (base::length(Branches_n_k[[i]]) > 0 && Branches_n_k[[i]][base::length(Branches_n_k[[i]])] == 0) {
              Branches_n_k[[i]] = Branches_n_k[[i]][-base::length(Branches_n_k[[i]])] # Remove the trailing 0 in Branches_n_k[[i]] caused by initialization.
            }
            Branches_n_k[[i]] = Branches_n_k[[i]][-base::length(Branches_n_k[[i]])]
            Branches_n_group[[i]] = Branches_k[[i]] # Recording final group count.
            break
          } else { # Sufficient remaining cells are available to form a new group.
            current_step = 2 # Proceeding to Step 2.
          }
        }

      }
      message("Finishing cell grouping for branch ", i, ".")
    }
    for (i in base::seq_along(Branches)) {
      if (base::length(Branches_n_k[[i]]) != Branches_n_group[[i]]) {
        stop("Grouping process failed for branch ", i, " with invalid cell group count records.")
      }
      if (base::sum(Branches_n_k[[i]]) != Branches_n_cell[[i]]) {
        stop("Grouping process failed for branch ", i, " because ungrouped cells remain, and the sum of grouped cells is less than the total cell count.")
      }
      if (Branches_n_r[[i]] != 0) {
        stop("Grouping process failed for branch ",i , " with ungrouped cells remaining.")
      }
    }

    ## Setting pseudotime t_k and interval length for each cell group.
    message("Setting pseudotime t_k and interval length for each cell group.")
    Branches_t_k = base::list()
    Branches_Tslot_k = base::list()
    for (i in base::seq_along(Branches)) {
      Branches_t_k[[i]] = Branches[[i]]$Pseudotime[
        base::match(base::colnames(Branches_TFE[[i]])[cumsum(Branches_n_k[[i]])], Branches[[i]]$`k~`)
      ]
      Branches_Tslot_k[[i]] = base::diff(c(0, base::as.numeric(Branches_t_k[[i]])))
    }

    ## Calculate mean of non-zero expression/activity values for each gene within each cell group.
    message("Setting aggregate expression values TFG^E_k and activity values TG^A_k, TG^E_k for each gene in every cell group.")
    # Constructing a function: The values of vector n_k_vector represent the number of elements in each segment of vector gene_value_vector, and compute the mean of each segment.
    calculate_segment_means = function(gene_value_vector, n_k_vector) {
      gene_value_vector = base::as.numeric(gene_value_vector)
      n_k_vector = base::as.numeric(n_k_vector)
      m = base::length(n_k_vector)
      starts = base::cumsum(c(1, n_k_vector[-m]))
      ends = base::cumsum(n_k_vector)
      means = base::vapply(1:m, function(k) {
        segment = gene_value_vector[starts[k]:ends[k]]
        segment = segment[segment != 0]
        if(base::length(segment) == 0) {return(NA)}
        base::mean(segment)
      }, numeric(1))
      return(means)
    }
    # Constructing a function to compute segmented means for each row of dataframe X (segmentation based on n_k_vector).
    apply_segment_means = function(X, n_k_vector) {
      X = base::as.data.frame(X)
      n_k_vector = base::as.numeric(n_k_vector)
      Y = base::as.data.frame(base::t(base::apply(X, 1, function(row) {
        calculate_segment_means(base::as.numeric(row), n_k_vector)
      })))
      base::colnames(Y) = 1:base::length(base::colnames(Y))
      return(Y)
    }
    # Setting aggregate expression values TFG^E_k and activity values TG^A_k, TG^E_k for each gene in every cell group (to be used as observational values for evolutionary equations).
    Branches_TFE_T = base::list()
    Branches_TGA_T = base::list()
    Branches_TGE_T = base::list()
    for (i in base::seq_along(Branches)) {
      Branches_TFE_T[[i]] = apply_segment_means(Branches_TFE[[i]], Branches_n_k[[i]])
      Branches_TGA_T[[i]] = apply_segment_means(Branches_TGA[[i]], Branches_n_k[[i]])
      Branches_TGE_T[[i]] = apply_segment_means(Branches_TGE[[i]], Branches_n_k[[i]])
    }

    return(base::list("Branches_n_group" = Branches_n_group,
                      "Branches_n_k" = Branches_n_k,
                      "n_min" = n_min,
                      "Branches_T_k" = Branches_t_k,
                      "Branches_Tslot_k" = Branches_Tslot_k,
                      "Branches_TFE_T" = Branches_TFE_T,
                      "Branches_TGA_T" = Branches_TGA_T,
                      "Branches_TGE_T" = Branches_TGE_T))
  }

  ### Invoking cell grouping function in parallel.
  interest_cell_type_group = parallel::mclapply(
    X = base::names(interest_cell_type_genes_pseudotime_info),
    FUN = function(cell_type) {
      message("Preparing for cell grouping across all branches in cell type ", cell_type, ".")
      cell_grouping_process(
        cell_type = cell_type,
        Branches_n_cell = interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]],
        Branches_TFE = interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]],
        Branches_TGA = interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGA"]],
        Branches_TGE = interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGE"]],
        Branches = interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches"]],
        calculate_n_min = calculate_n_min,
        points_x_for_fitting_nc_and_nmin = points_x_for_fitting_nc_and_nmin,
        points_y_for_fitting_nc_and_nmin = points_y_for_fitting_nc_and_nmin
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_group) = base::names(interest_cell_type_genes_pseudotime_info)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Grouping cells based on pseudotime points and valid (non-zero) counts of expression and activity values. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_group)
}

#' @title Pseudotime Ordering and Interactive Branch Partitioning for scRNA-seq Data
#'
#' @description
#' Builds developmental trajectories and performs pseudotime ordering for each cell type using Monocle2 (DDRTree algorithm).
#' Supports both automatic and manual parameter tuning for trajectory construction, and divides cells into branches via interactive user input based on visualization results.
#'
#' @param interest_cell_type_data List. The output of function \code{\link{get_interest_cell_type_data}},
#'   containing filtered scRNA-seq objects and gene sets for each cell type.
#' @param auto_param_selection Logical. Whether to enable automatic parameter selection for DDRTree. Defaults to \code{TRUE}.
#'   When \code{TRUE}, trajectory is constructed with default DDRTree parameters.
#'   When \code{FALSE}, you can manually adjust trajectory topology via parameters below.
#' @param maxIter Integer. Maximum iterations for DDRTree. Only used when \code{auto_param_selection = FALSE}.
#' @param sigma Numeric. Bandwidth parameter for DDRTree. Only used when \code{auto_param_selection = FALSE}.
#' @param lambda Numeric. Regularization parameter for DDRTree. Only used when \code{auto_param_selection = FALSE}.
#' @param ncenter Integer. Number of centroids for DDRTree. Only used when \code{auto_param_selection = FALSE}.
#' @param param.gamma Numeric. Gamma parameter for DDRTree. Only used when \code{auto_param_selection = FALSE}.
#' @param tol Numeric. Convergence tolerance for DDRTree. Only used when \code{auto_param_selection = FALSE}.
#'
#' @return
#' A nested list structure:
#' \itemize{
#'   \item First level: Named list with cell types as names
#'   \item Second level: List of branches for each cell type, ordered by user input sequence
#'   \item Each branch is a data frame with rows as cells, sorted by pseudotime in ascending order, containing two columns:
#'     \itemize{
#'       \item \code{Pseudotime}: Pseudotime value of each cell
#'       \item \code{State}: Assigned cell state from trajectory inference
#'     }
#' }
#'
#' @details
#' For each cell type, the function executes the following workflow:
#' \enumerate{
#'   \item Converts the scRNA-seq count matrix into a Monocle CellDataSet object
#'   \item Estimates size factors and dispersions for expression normalization
#'   \item Sets the input gene set (\code{interest_gene_for_GRN}) as features for trajectory inference
#'   \item Performs dimensionality reduction and trajectory construction via DDRTree (auto or manual parameter mode)
#'   \item Orders cells along the trajectory to assign pseudotime values
#'   \item Exports trajectory visualization plots and saves the full CDS object
#'   \item Accepts interactive user input: number of branches and cell state sequences for each branch
#'   \item Extracts cells for each defined branch, sorts by pseudotime, and returns the result
#' }
#'
#' @section Interactive Input:
#' The function pauses for user input per cell type. You need to:
#' \itemize{
#'   \item Enter the total number of branches (positive integer) based on the generated trajectory plots
#'   \item Enter cell state IDs (space-separated integers) for each branch in sequence, matching the State labels in the visualization
#' }
#' Please check the output images in the cell type subdirectory before entering values.
#'
#' @note
#' \itemize{
#'   \item Depends on \code{monocle} and \code{DDRTree} packages. It is recommended to use \code{igraph} version 2.0.3, as higher versions may cause compatibility issues with Monocle.
#'   \item Creates a top-level directory \code{2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment}, with subdirectories named by cell type.
#'     Visualization plots and the CDS object are saved in the corresponding subdirectory. Working directory is restored on exit.
#'   \item If the cell count is too small or the default trajectory topology is not satisfactory, set \code{auto_param_selection = FALSE} and tune DDRTree parameters manually.
#'   \item This is an interactive function. It will pause for input once per cell type during execution.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform pseudotime analysis and branch partitioning
#' library(DDRTree)
#' packageVersion("igraph")
#' remove.packages("igraph")
#' install_version("igraph", version = "2.0.3")
#' interest_cell_type_data = base::readRDS(paste0(
#'   "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'   "/interest_cell_type_data.rds"
#' ))
#' interest_cell_type_Branches = order_pseudotime_and_divide_branches(interest_cell_type_data)
#' # ...
#' # Please enter the number of branches based on the pseudotime 2D visualization:
#' #   6
#' # ...
#' # Please enter the integer vector of cell state order for branch 1 (space-separated):
#' #   1 11
#' # Please enter the integer vector of cell state order for branch 2 (space-separated):
#' #   1 2 10
#' # Please enter the integer vector of cell state order for branch 3 (space-separated):
#' #   1 2 3 9
#' # Please enter the integer vector of cell state order for branch 4 (space-separated):
#' #   1 2 3 4 8
#' # Please enter the integer vector of cell state order for branch 5 (space-separated):
#' #   1 2 3 4 5 6
#' # Please enter the integer vector of cell state order for branch 6 (space-separated):
#' #   1 2 3 4 5 7
#' # ...
#'
#' interest_cell_type_Branches = order_pseudotime_and_divide_branches(
#'    interest_cell_type_data = interest_cell_type_data,
#'    auto_param_selection = FALSE,
#'    ncenter = 3
#' )
#' # ...
#' # Please enter the number of branches based on the pseudotime 2D visualization:
#' #    6
#' # ...
#' # Please enter the integer vector of cell state order for branch 1 (space-separated):
#' #    1
#' # ...
#'
#' base::saveRDS(
#'   interest_cell_type_Branches,
#'   file = paste0(
#'     "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'     "/interest_cell_type_Branches.rds"
#'   )
#' )
#' }
order_pseudotime_and_divide_branches = function(
    interest_cell_type_data,
    auto_param_selection = TRUE,
    maxIter = 20,
    sigma = 0.001,
    lambda = NULL,
    ncenter = NULL,
    param.gamma = 10,
    tol = 0.001
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Performing pseudotime ordering and branch partitioning for cell types of interest. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Defining a function for pseudotime ordering and branch partitioning of individual cell types.
  order_pseudotime_and_divide_branches_process = function (interest_scRNAseq_for_GRN, interest_gene_for_GRN) {
    # Converting counts expression matrix to sparse matrix format and creating a CellDataSet object.
    message("Converting counts expression matrix to sparse matrix format and creating a CellDataSet object.")
    interest_RNAsparseMatrix_for_GRN =
      methods::as(base::as.matrix(SeuratObject::JoinLayers(interest_scRNAseq_for_GRN)@assays$RNA$counts), 'sparseMatrix')
    interest_cds = monocle::newCellDataSet(
      cellData = interest_RNAsparseMatrix_for_GRN,
      phenoData =
        methods::new('AnnotatedDataFrame',
                     data = SeuratObject::JoinLayers(interest_scRNAseq_for_GRN)@meta.data),
      featureData =
        methods::new('AnnotatedDataFrame',
                     data = base::data.frame(gene_short_name = base::row.names(interest_RNAsparseMatrix_for_GRN),
                                             row.names = base::row.names(interest_RNAsparseMatrix_for_GRN))),
      lowerDetectionLimit = 0.1,
      expressionFamily = VGAM::negbinomial.size()
    )

    # Estimating normalization factors and dispersion.
    message("Estimating normalization factors and dispersion.")
    interest_cds = BiocGenerics::estimateSizeFactors(interest_cds)
    interest_cds = BiocGenerics::estimateDispersions(interest_cds)

    # Specifying genes for subsequent trajectory inference analysis.
    message("Specifying genes for subsequent trajectory inference analysis.")
    interest_cds = monocle::setOrderingFilter(interest_cds, interest_gene_for_GRN)

    # Dimensionality reduction.
    message("Dimensionality reduction.")
    # library(DDRTree)
    if (auto_param_selection) {
      interest_cds = monocle::reduceDimension(
        interest_cds, reduction_method = 'DDRTree', max_components = 2
      )
    } else{
      interest_cds = monocle::reduceDimension(
        cds = interest_cds, max_components = 2, reduction_method = 'DDRTree',
        auto_param_selection = auto_param_selection,
        maxIter = maxIter, sigma = sigma, lambda = lambda,
        ncenter = ncenter, param.gamma = param.gamma, tol = tol
      )
    }

    # Ordering cells by trajectory (pseudotime).
    message("Ordering cells by trajectory (pseudotime).")
    message("It is recommended to use igraph R package version 2.0.3, as later versions are not compatible with monocle.")
    interest_cds = monocle::orderCells(interest_cds)

    message("Saving the CDS variable 'interest_cds' for monocle pseudotime analysis.")
    base::saveRDS(interest_cds, file = "interest_cds.rds")

    # Generating 2D visualizations by pseudotime and cell state saved in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment.
    message("Generating 2D visualizations by pseudotime and cell states saved in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment.")
    ggplot2::ggsave("interest_AllCellState_in_trajectory.png", width = 5, height = 5,
                    plot = monocle::plot_cell_trajectory(interest_cds, color_by = "State"))
    ggplot2::ggsave(
      "interest_EveryCellState_in_trajectory.png",
      width = 3 * base::ceiling(
        base::length(base::unique(Biobase::pData(interest_cds)$State)) /
          base::ceiling(base::sqrt(base::length(base::unique(Biobase::pData(interest_cds)$State))))
      ),
      height = 3 * base::ceiling(
        base::sqrt(base::length(base::unique(Biobase::pData(interest_cds)$State)))
      ),
      plot = monocle::plot_cell_trajectory(interest_cds, color_by = "State") +
        ggplot2::facet_wrap(
          ~State,
          nrow = base::ceiling(base::sqrt(base::length(base::unique(Biobase::pData(interest_cds)$State)))),
          axes = "all"
        )
    )
    ggplot2::ggsave("interest_Pseudotime_in_trajectory.png", width = 5, height = 5,
                    plot = monocle::plot_cell_trajectory(interest_cds, color_by = "Pseudotime"))

    # Interactively input the number of branches.
    message("Interactively input the number of branches (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).")
    n_branch = base::as.integer(base::readline(prompt = "Please enter the number of branches based on the pseudotime 2D visualization (e.g., 3) (e.g., 5): "))
    if (base::is.na(n_branch) || n_branch <= 0.9999999) {
      stop("Invalid input. Please enter a positive integer for the number of branches.")
    }

    # Interactively input integer vector of cell state order for each branch.
    message("Interactively input the integer vector of cell state order for each branch (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).")
    branches = base::list()
    for (i in 1:n_branch) {
      base::cat(base::sprintf("Please enter the integer vector of cell state order for branch %d (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5): ", i))
      integer_vector = base::as.integer(base::strsplit(base::readline(), " ")[[1]])
      if (base::any(base::is.na(integer_vector))) {
        stop(base::sprintf("Invalid input. Please ensure the input for branch %d is an integer vector.", i))
      }
      branches[[i]] = integer_vector
    }

    # Extracting pseudotime and cell states for each branch.
    message("Extracting pseudotime and cell states for each branch.")
    Branches = base::list()
    for (i in base::seq_along(branches)) {
      Branches[[i]] = Biobase::pData(interest_cds)[
        Biobase::pData(interest_cds)$State %in% branches[[i]],
      ][, c("Pseudotime", "State")]
      Branches[[i]] = Branches[[i]][base::order(Branches[[i]][, 1]), ]
    }

    return(Branches)
  }

  ### Executing the defined function for each cell type.
  interest_cell_type_Branches = stats::setNames(
    base::vector("list", base::length(interest_cell_type_data)), base::names(interest_cell_type_data)
  )
  for (name in base::names(interest_cell_type_Branches)) {
    # Start.
    message("Performing Pseudotime Ordering and Branch Partitioning for ", name, ".")
    original_dir1 = base::getwd()
    new_folder1 = name
    if (!base::dir.exists(new_folder1)) {
      base::dir.create(new_folder1, recursive = TRUE)
      message("Folder created: ", new_folder1, ".")
    } else {message("Folder already exists: ", new_folder, ".")}
    base::setwd(new_folder1)
    message("The current working directory has been switched to: ", base::getwd(), ".")

    # Do.
    interest_cell_type_Branches[[name]] = order_pseudotime_and_divide_branches_process(
      interest_scRNAseq_for_GRN = interest_cell_type_data[[name]]$interest_scRNAseq_for_GRN,
      interest_gene_for_GRN = interest_cell_type_data[[name]]$interest_gene_for_GRN
    )

    # End.
    base::setwd(original_dir1)
    message("The current working directory has been switched to: ", base::getwd(), ".")
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Performing pseudotime ordering and branch partitioning for cell types of interest. Time: ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_Branches)
}

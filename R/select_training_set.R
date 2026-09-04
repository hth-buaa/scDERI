#' @title Select Training Sets for Gene Regulatory Network Modeling
#'
#' @description
#' Extracts high-confidence TF-TG interactions from branch-level initial GRNs to
#' construct training sets for downstream regulatory strength prediction.
#' It implements multi-layer threshold filtering, optional refinement with TF coverage guarantee,
#' and appends pseudotime-resolved expression/activity features for each branch.
#'
#' @param interest_cell_type_branch_iGRN Named list. Output from \code{\link{get_branch_iGRN}}, containing branch-specific iGRN matrices for each cell type.
#' @param interest_cell_type_tau Named numeric vector. Output from \code{\link{get_cell_type_tau}}, with cell-type-specific binarization thresholds.
#' @param interest_cell_type_group Named list. Output from \code{\link{cell_grouping}}, containing cell grouping results and aggregated gene expression/activity matrices per branch.
#' @param top_n_for_max_theta_in_each_regulon Positive integer.
#' Number of top interactions retained per TF and per TG separately;
#' their union forms the initial training set. Defaults to 1.
#' @param max_theta_retention_rate_in_tops Numeric between 0 and 1.
#' Quantile retention ratio of regulatory strength within the initial training set. Defaults to 1 (no filtering).
#' @param max_theta_retention_rate_in_all Numeric between 0 and 1.
#' Quantile retention ratio of regulatory strength across the global candidate set. Defaults to 1 (no filtering).
#' @param retention_number_of_tops Positive integer.
#' Absolute number of top-ranked interactions by regulatory strength. Defaults to 1000.
#' @param refine Logical.
#' Whether to perform additional mutual-top refinement on the training set. Defaults to \code{TRUE}.
#' @param keep_all_TF Logical. Whether to ensure every TF has at least one interaction retained after refinement.
#' Only works when \code{refine = TRUE}. Defaults to \code{TRUE}.
#' @param top_ratio_in_refine Numeric between 0 and 1.
#' Quantile retention ratio of regulatory strength after the refinement step.
#' Only works when \code{refine = TRUE}. Defaults to 1 (no filtering).
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#'
#' @return Nested list with the structure \code{list[[cell_type]][[branch]]}:
#' \itemize{
#'   \item Each branch entry is a data frame with TF-TG pairs as rows.
#'   Columns include \code{TF}, \code{TG}, \code{theta_i},
#'   plus time-point-specific features \code{TFE_K}, \code{TGA_K}, \code{TGE_K}.
#'   \item An \code{all} entry per cell type stores the global merged training set without time-point features.
#' }
#'
#' @details The workflow is executed independently for each cell type:
#' \enumerate{
#'   \item Generates all TF-TG pairs for each branch, attaches regulatory strength \code{theta_i}, and filters out pairs below the cell-type-specific \code{tau} threshold. Merges and deduplicates pairs across all branches to form the global candidate set.
#'   \item Builds the initial training set as the union of top-N interactions ranked by \code{theta_i} for each TF and for each TG separately.
#'   \item Marks three filter criteria on the initial training set (filtering is not applied at this stage):
#'     \itemize{
#'       \item Meets the within-training-set quantile threshold
#'       \item Meets the global candidate-set quantile threshold
#'       \item Ranks within the absolute top-N by regulatory strength
#'     }
#'   \item Optional refinement step (when \code{refine = TRUE}):
#'     \itemize{
#'       \item Retains pairs where the TF is in the top regulator set of the TG AND the TG is in the top target set of the TF
#'       \item Applies an additional quantile filter by \code{theta_i}
#'       \item If \code{keep_all_TF = TRUE}, adds back the strongest interaction for each missing TF to ensure full TF coverage
#'     }
#'   \item Applies all three pre-marked filter criteria; only pairs satisfying all three are retained.
#'   \item Maps the final training set back to each branch: removes pairs absent from the current branch, and appends corresponding time-resolved TF expression, TG activity, and TG expression features.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates a \code{4.1 Build Prediction Model - Select Training Set} directory and temporarily switches the working directory into it; the original path is restored on exit.
#'   \item Numeric filter parameters are automatically clipped to valid ranges, and integer parameters are rounded internally.
#'   \item Parallel computation relies on \code{\link[parallel]{mclapply}} and only works on Unix-like systems.
#'   \item The refinement step can substantially reduce the number of interactions. Adjust \code{top_ratio_in_refine} to balance coverage and confidence.
#'   \item Training set selection is fully based on regulatory strength from the iGRN. Ensure consistent upstream data processing.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_branch_iGRN.rds")
#' interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
#' interest_cell_type_group = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' top_n_for_max_theta_in_each_regulon = 1
#' max_theta_retention_rate_in_tops = 0.8
#' max_theta_retention_rate_in_all = 0.8
#' retention_number_of_tops = 1000
#' refine = TRUE
#' keep_all_TF = TRUE
#' top_ratio_in_refine = 0.8
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_branch_training_set = select_training_set(
#'   interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
#'   interest_cell_type_tau  = interest_cell_type_tau ,
#'   interest_cell_type_group = interest_cell_type_group,
#'   top_n_for_max_theta_in_each_regulon = top_n_for_max_theta_in_each_regulon,
#'   max_theta_retention_rate_in_tops = max_theta_retention_rate_in_tops,
#'   max_theta_retention_rate_in_all = max_theta_retention_rate_in_all,
#'   retention_number_of_tops = retention_number_of_tops,
#'   refine = refine,
#'   keep_all_TF = keep_all_TF,
#'   top_ratio_in_refine = top_ratio_in_refine,
#'   ncores = ncores
#' )
#' }
select_training_set = function(interest_cell_type_branch_iGRN,
                               interest_cell_type_tau,
                               interest_cell_type_group,
                               top_n_for_max_theta_in_each_regulon = 1,
                               max_theta_retention_rate_in_tops = 1,
                               max_theta_retention_rate_in_all = 1,
                               retention_number_of_tops = 1000,
                               refine = TRUE,
                               keep_all_TF = TRUE,
                               top_ratio_in_refine = 1,
                               ncores = 1)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Partitioning Training Sets. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "4.1 Build Prediction Model - Select Training Set"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Input validation.
  if (!is.numeric(top_n_for_max_theta_in_each_regulon) ||
      top_n_for_max_theta_in_each_regulon < 1 ||
      top_n_for_max_theta_in_each_regulon != round(top_n_for_max_theta_in_each_regulon)) {
    stop("'top_n_for_max_theta_in_each_regulon' must be a positive integer.")
  }
  if (!is.numeric(max_theta_retention_rate_in_tops) ||
      max_theta_retention_rate_in_tops < 0 || max_theta_retention_rate_in_tops > 1) {
    stop("'max_theta_retention_rate_in_tops' must be between 0 and 1.")
  }
  if (!is.logical(refine) || length(refine) != 1) {
    stop("'refine' must be TRUE or FALSE.")
  }

  # Core processing function: select training set for one cell type (all branches combined)
  select_training_set_process = function(branch_GRN,
                                         branch_tau,
                                         top_n_for_max_theta_in_each_regulon = 1,
                                         max_theta_retention_rate_in_tops = 1,
                                         max_theta_retention_rate_in_all = 1,
                                         retention_number_of_tops = 1000,
                                         refine = TRUE,
                                         keep_all_TF = TRUE,
                                         top_ratio_in_refine = 1,
                                         TFE_T,
                                         TGA_T,
                                         TGE_T) {
    message("Generating all row-column (TG-TF) combinations, and adding regulatory strength theta for each TG_TF pair.")
    n_branch = base::length(branch_GRN)
    all_combinations = base::list()

    # Generate all TF-TG pairs per branch, attach theta_i, filter by tau threshold
    for (n in 1:n_branch) {
      all_combinations[[base::paste0("branch_", n, "_GRN")]] = base::expand.grid(
        TG = base::rownames(branch_GRN[[n]]), TF = base::colnames(branch_GRN[[n]]), stringsAsFactors = FALSE
      )
      all_combinations[[base::paste0("branch_", n, "_GRN")]]$theta_i = base::apply(
        all_combinations[[base::paste0("branch_", n, "_GRN")]], 1, function(x) branch_GRN[[n]][x[1], x[2]]
      )
      all_combinations[[base::paste0("branch_", n, "_GRN")]] = base::subset(
        all_combinations[[base::paste0("branch_", n, "_GRN")]],
        theta_i >= branch_tau
      )
    }

    # Combine pairs from all branches, remove duplicates to get global candidate set
    all_combinations[["GRN"]] = base::do.call(base::rbind, all_combinations)
    all_combinations[["GRN"]] = all_combinations[["GRN"]][!base::duplicated(all_combinations[["GRN"]]), ]

    # Stratified filtering: build training set from global candidate pairs with multiple thresholds
    stratified_slice = function(data,
                                top_n_for_max_theta_in_each_regulon = 1,
                                max_theta_retention_rate_in_tops = 1,
                                max_theta_retention_rate_in_all = 1,
                                retention_number_of_tops = 1000,
                                refine = TRUE,
                                keep_all_TF = TRUE,
                                top_ratio_in_refine = 1) {
      # Auto-clip parameters to valid ranges and coerce to integer where needed
      message("Normalizing and clipping filter parameters to valid ranges.")
      top_n_for_max_theta_in_each_regulon = base::ifelse(
        top_n_for_max_theta_in_each_regulon < 1.5,
        1,
        base::round(top_n_for_max_theta_in_each_regulon)
      )
      max_theta_retention_rate_in_tops = base::pmin(
        base::pmax(max_theta_retention_rate_in_tops, 0), 1
      )
      max_theta_retention_rate_in_all = base::pmin(
        base::pmax(max_theta_retention_rate_in_all, 0), 1
      )

      message("Training set strategy: take union of top ", top_n_for_max_theta_in_each_regulon,
              " interactions per TF and per TG; 3 filters are applied after refinement.")

      message("Select the top ", top_n_for_max_theta_in_each_regulon, " thetas with the largest values for each TF.")
      TFregulon_top_n = base::list()
      for (TF in unique(data$TF)) {
        TFregulon = data[data$TF == TF, ]
        TFregulon_top_n[[TF]] = dplyr::slice_max(
          .data = TFregulon,
          order_by = theta_i,
          n = top_n_for_max_theta_in_each_regulon,
          with_ties = FALSE
        )
      }
      TFregulon_top_n = data.table::rbindlist(TFregulon_top_n)

      message("Select the top ", top_n_for_max_theta_in_each_regulon, " thetas with the largest values for each TG.")
      TGregulon_top_n = base::list()
      for (TG in unique(data$TG)) {
        TGregulon = data[data$TG == TG, ]
        TGregulon_top_n[[TG]] = dplyr::slice_max(
          .data = TGregulon,
          order_by = theta_i,
          n = top_n_for_max_theta_in_each_regulon,
          with_ties = FALSE
        )
      }
      TGregulon_top_n = data.table::rbindlist(TGregulon_top_n)

      # Build initial training set: union of top-N interactions per TF and per TG
      message("Building initial training set as union of top-N per-TF and per-TG pairs.")
      train_set = base::unique(base::rbind(TFregulon_top_n, TGregulon_top_n))

      # Mark three filter criteria on the initial training set (applied after refinement)
      message("Marking three filter criteria on the initial training set.")
      train_set$max_theta_retention_rate_in_tops = base::ifelse(
        train_set$theta_i >= stats::quantile(
          train_set$theta_i,
          probs = 1 - max_theta_retention_rate_in_tops,
          na.rm = TRUE
        ), TRUE, FALSE
      )
      train_set$max_theta_retention_rate_in_all = base::ifelse(
        train_set$theta_i >= stats::quantile(
          data$theta_i,
          probs = 1 - max_theta_retention_rate_in_all,
          na.rm = TRUE
        ), TRUE, FALSE
      )
      train_set = dplyr::mutate(
        train_set,
        retention_number_of_tops = dplyr::row_number(dplyr::desc(theta_i)) <= retention_number_of_tops
      )
      message("There are ", base::nrow(train_set), " regulons in training set.")
      message("There are ", base::length(base::unique(train_set$TF)), " TFs in training set.")
      message("TFs: ", base::paste(base::sort(base::unique(train_set$TF)), collapse = ", "), ".")
      message("There are ", base::length(base::unique(train_set$TG)), " TGs in training set.")
      message("TGs: ", base::paste(base::sort(base::unique(train_set$TG)), collapse = ", "), ".")

      # Refinement step: filter by global top regulon membership and additional theta quantile
      if (refine == TRUE) {
        message("Running refinement step: global regulon membership filtering + theta quantile cutoff.")
        train_set0 = train_set
        train_set = train_set[
          (train_set$TF %in% base::unique(TGregulon_top_n$TF)) &
            (train_set$TG %in% base::unique(TFregulon_top_n$TG)),
        ]
        message("Retaining top ", top_ratio_in_refine * 100, "% of interactions by theta_i after refinement.")
        train_set = train_set[train_set$theta_i >= stats::quantile(
          train_set$theta_i, probs = 1 - top_ratio_in_refine , na.rm = TRUE
        ), ]

        message("There are ", base::nrow(train_set), " regulons in training set after refining.")
        message("There are ", base::length(base::unique(train_set$TF)), " TFs in training set after refining.")
        message("TFs: ", base::paste(base::sort(base::unique(train_set$TF)), collapse = ", "), ".")
        message("There are ", base::length(base::unique(train_set$TG)), " TGs in training set after refining.")
        message("TGs: ", base::paste(base::sort(base::unique(train_set$TG)), collapse = ", "), ".")

        # Ensure full TF coverage: add back missing TFs with their strongest interaction
        if (keep_all_TF == TRUE) {
          TFs_filter = base::setdiff(base::unique(train_set0$TF), base::unique(train_set$TF))

          message("There are ", base::length(TFs_filter), " regulons that are needed to add.")
          message("Added TFs: ", base::paste(base::sort(TFs_filter), collapse = ", "), ".")

          # message("Whether add: ", base::length(TFs_filter) > 0)
          if (base::length(TFs_filter) > 0) {
            for (TF_filter in TFs_filter) {
              # message("Added TF ", TF_filter, ".")
              TF_filter_train_set = train_set0[train_set0$TF == TF_filter, , drop = FALSE]
              TF_filter_train_set = TF_filter_train_set[base::which.max(TF_filter_train_set$theta_i), , drop = FALSE]
              # message("Added sample: ")
              # base::print(TF_filter_train_set)
              train_set = base::rbind(train_set, TF_filter_train_set)
            }

            message("There are ", base::nrow(train_set), " regulons in training set after refining.")
            message("There are ", base::length(base::unique(train_set$TF)), " TFs in training set after refining.")
            message("TFs: ", base::paste(base::sort(base::unique(train_set$TF)), collapse = ", "), ".")
            message("There are ", base::length(base::unique(train_set$TG)), " TGs in training set after refining.")
            message("TGs: ", base::paste(base::sort(base::unique(train_set$TG)), collapse = ", "), ".")
          }
        }
      }

      # Apply all three pre-marked filter criteria; retain only pairs passing all three
      message("Applying final triple filters and formatting output training set.")
      train_set = train_set[rowSums(train_set[, c(
        "max_theta_retention_rate_in_tops",
        "max_theta_retention_rate_in_all",
        "retention_number_of_tops")]) == 3, ]
      train_set = base::as.data.frame(train_set[, c("TF", "TG", "theta_i")])
      # message("The training set is:")
      # base::print(train_set)
      train_set = train_set[!base::duplicated(train_set), ]
      rownames(train_set) = base::paste0(train_set$TF, "_to_", train_set$TG)
      # message("After adding rownames, the training set is:")
      # message("The training set is:")
      # base::print(train_set)
      return(train_set)
    }
    slices_GRN = stratified_slice(
      data = all_combinations[["GRN"]],
      top_n_for_max_theta_in_each_regulon = top_n_for_max_theta_in_each_regulon,
      max_theta_retention_rate_in_tops = max_theta_retention_rate_in_tops,
      max_theta_retention_rate_in_all = max_theta_retention_rate_in_all,
      retention_number_of_tops = retention_number_of_tops,
      refine = refine,
      keep_all_TF = keep_all_TF,
      top_ratio_in_refine = top_ratio_in_refine
    )

    message("Adding cell group states (TFE, TGA, TGE).")
    # Attach branch-specific cell state features; remove pairs absent from the branch
    slices_branch_GRN = base::list()
    for (n in 1:n_branch) {
      slices_branch_GRN[[base::paste0("branch", n)]] = slices_GRN
      for (TG_TF in base::rownames(slices_branch_GRN[[base::paste0("branch", n)]])) {
        TF = base::strsplit(TG_TF, "_to_")[[1]][1]
        TG = base::strsplit(TG_TF, "_to_")[[1]][2]
        if (TF %in% base::colnames(branch_GRN[[n]]) && TG %in% base::rownames(branch_GRN[[n]])) {
          for (K in 1:ncol(TFE_T[[n]])) {
            slices_branch_GRN[[base::paste0("branch", n)]][TG_TF, base::paste0("TFE_", K)] = TFE_T[[n]][TF, K]
            slices_branch_GRN[[base::paste0("branch", n)]][TG_TF, base::paste0("TGA_", K)] = TGA_T[[n]][TG, K]
            slices_branch_GRN[[base::paste0("branch", n)]][TG_TF, base::paste0("TGE_", K)] = TGE_T[[n]][TG, K]
          }
        } else {
          slices_branch_GRN[[base::paste0("branch", n)]] = slices_branch_GRN[[base::paste0("branch", n)]][
            !rownames(slices_branch_GRN[[base::paste0("branch", n)]]) %in% TG_TF,
          ]
        }
      }
    }
    slices_branch_GRN[["all"]] = slices_GRN
    return(slices_branch_GRN)
  }
  interest_cell_type_branch_training_set = parallel::mclapply(
    X = base::names(interest_cell_type_branch_iGRN),
    FUN = function(cell_type) {
      message("Partitioning Training Sets from the data of ", cell_type, ".")
      select_training_set_process(
        branch_GRN = interest_cell_type_branch_iGRN[[cell_type]],
        branch_tau = interest_cell_type_tau[cell_type],
        top_n_for_max_theta_in_each_regulon = top_n_for_max_theta_in_each_regulon,
        max_theta_retention_rate_in_tops = max_theta_retention_rate_in_tops,
        max_theta_retention_rate_in_all = max_theta_retention_rate_in_all,
        retention_number_of_tops = retention_number_of_tops,
        refine = refine,
        keep_all_TF = keep_all_TF,
        top_ratio_in_refine = top_ratio_in_refine,
        TFE_T = interest_cell_type_group[[cell_type]][["Branches_TFE_T"]],
        TGA_T = interest_cell_type_group[[cell_type]][["Branches_TGA_T"]],
        TGE_T = interest_cell_type_group[[cell_type]][["Branches_TGE_T"]]
      )
    },
    mc.cores = ncores
  )

  ### End.
  base::names(interest_cell_type_branch_training_set) = base::names(interest_cell_type_branch_iGRN)
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Partitioning Training Sets. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_branch_training_set)
}

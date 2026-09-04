#' @title Prepare Data for GRN Inference
#'
#' @description
#' Generates inference-ready datasets by combining branch-level iGRN TF-TG pairs
#' with time-series cell state features from cell grouping results. Supports partition-based
#' parallel processing for large datasets, and outputs both branch-level full feature sets
#' and a global de-duplicated regulatory pair set.
#'
#' @param interest_cell_type_branch_iGRN List. Output from \code{\link{get_branch_iGRN}},
#'   containing branch-specific iGRN matrices (rows: target genes, columns: TFs) for each cell type.
#' @param interest_cell_type_group List. Output from \code{\link{cell_grouping}},
#'   containing aggregated TFE/TGA/TGE matrices per branch.
#' @param n_part Positive integer. Number of data partitions for parallel feature binding. Defaults to \code{1}.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#'
#' @return
#' Nested list organized by cell type. Each cell type contains:
#' \itemize{
#'   \item \code{branch1}, \code{branch2}, ...: Branch-specific data frames, each row corresponds to one TF-TG pair.
#'     Columns include \code{TG}, \code{TF}, \code{theta_i}, plus time-point features
#'     \code{TFE_K}, \code{TGA_K}, \code{TGE_K} for each time point K. Pairs with \code{theta_i = 0} are filtered out.
#'   \item \code{GRN}: Global de-duplicated TF-TG pair data frame across all branches,
#'     containing only \code{TG}, \code{TF} and \code{theta_i} columns (no time-series features).
#'     Pairs with \code{theta_i = 0} are filtered out.
#' }
#'
#' @details
#' The function runs independently for each cell type with the following workflow:
#' \enumerate{
#'   \item Generates all TG-TF combinations from each branch's iGRN matrix, and extracts
#'     the corresponding regulatory strength \code{theta_i} for each pair.
#'   \item Merges TF-TG pairs across all branches and removes duplicates to form the global GRN set.
#'   \item Splits each branch's pair set into \code{n_part} partitions, and appends
#'     time-series cell state features (TFE, TGA, TGE) to each partition in parallel.
#'   \item Combines processed partitions back into a complete branch-level dataset.
#'   \item Filters out TF-TG pairs with \code{theta_i} equal to 0 for both branch-level and global datasets.
#'   \item Returns the nested list with branch-level full feature sets and the global pair set.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates a \code{5 Infer GRN} directory and temporarily switches
#'     the working directory, restoring the original path on exit.
#'   \item Partition-based parallelism accelerates feature binding for large pair sets;
#'     adjust \code{n_part} according to dataset size and available memory.
#'   \item Output automatically excludes pairs with zero regulatory strength, retaining
#'     only predicted regulatory relationships.
#'   \item The \code{GRN} element provides a global, de-duplicated view of all regulatory
#'     pairs across branches, without time-series features.
#'   \item Parallel computation relies on \code{parallel::mclapply}; set \code{ncores = 1}
#'     on Windows systems.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_branch_iGRN.rds")
#' interest_cell_type_group = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' n_part = ncores
#' interest_cell_type_branch_data_for_GRN_infer = get_data_for_inferring_GRN(
#'   interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
#'   interest_cell_type_group = interest_cell_type_group,
#'   n_part = n_part,
#'   ncores = ncores
#' )
#' base::saveRDS(
#'   interest_cell_type_branch_data_for_GRN_infer,
#'   file = "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds"
#' )
#' }
get_data_for_inferring_GRN = function(
    interest_cell_type_branch_iGRN,
    interest_cell_type_group,
    n_part = 1,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Getting data for inferring GRN for each branch. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "5 Infer GRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Process.
  interest_cell_type_branch_data_for_GRN_infer = base::list()
  for (cell_type in names(interest_cell_type_branch_iGRN)) {
    message("Process cell type ", cell_type, ".")
    interest_cell_type_branch_data_for_GRN_infer[[cell_type]] = base::list()
    branch_GRN = interest_cell_type_branch_iGRN[[cell_type]]

    # Generate all TF-TG pairs for each branch and extract regulatory strength theta_i.
    message("Generating all row-column (TG-TF) combinations, and adding regulatory strength theta for each TF-TG pair.")
    n_branch = base::length(branch_GRN)
    all_combinations = base::list()
    for (n in 1:n_branch) {
      all_combinations[[base::paste0("branch_", n, "_GRN")]] = base::expand.grid(
        TG = base::rownames(branch_GRN[[n]]), TF = base::colnames(branch_GRN[[n]]), stringsAsFactors = FALSE
      )
      all_combinations[[base::paste0("branch_", n, "_GRN")]]$theta_i = base::apply(
        all_combinations[[base::paste0("branch_", n, "_GRN")]], 1, function(x) branch_GRN[[n]][x[1], x[2]]
      )
      rownames(all_combinations[[base::paste0("branch_", n, "_GRN")]]) = paste0(
        all_combinations[[base::paste0("branch_", n, "_GRN")]]$TF,
        "_to_",
        all_combinations[[base::paste0("branch_", n, "_GRN")]]$TG
      )
    }
    # Merge all branch pairs and remove duplicates to get global TF-TG set.
    all_combinations[["GRN"]] = base::do.call(rbind, all_combinations)
    all_combinations[["GRN"]] = all_combinations[["GRN"]][!base::duplicated(all_combinations[["GRN"]]), ]
    rownames(all_combinations[["GRN"]]) = paste0(
      all_combinations[["GRN"]]$TF,
      "_to_",
      all_combinations[["GRN"]]$TG
    )

    # Add time-series cell state features (TFE, TGA, TGE) for each branch.
    message("Adding cell group states (TFE, TGA, TGE).")
    for (n in 1:n_branch) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      GRN = all_combinations[[base::paste0("branch_", n, "_GRN")]]
      N_re_all = base::nrow(GRN)
      message("In branch ", n, ", there are ", N_re_all, " potential regulatory relationships before filtering.")

      # Partition data for parallel processing.
      message("Partitioning data.")
      each = base::floor(N_re_all / n_part)
      # each_last = N_re_all - each * (n_part - 1)
      GRN_list = base::list()
      for (s in 1:n_part) {
        if (s != n_part) {
          GRN_list[[s]] = GRN[(each * (s - 1) + 1):(each * s), ]
        } else (
          GRN_list[[s]] = GRN[(each * (s - 1) + 1):N_re_all, ]
        )
      }

      # Parallel processing: add cell state features to each partition.
      message("Adding cell group states (TFE, TGA, TGE) in parallel.")
      process_GRN_slice = function(s,
                                   GRN_list = GRN_list,
                                   TFs = base::colnames(branch_GRN[[n]]),
                                   TGs = base::rownames(branch_GRN[[n]]),
                                   TFE_T = interest_cell_type_group[[cell_type]][["Branches_TFE_T"]],
                                   TGA_T = interest_cell_type_group[[cell_type]][["Branches_TGA_T"]],
                                   TGE_T = interest_cell_type_group[[cell_type]][["Branches_TGE_T"]]) {
        current_branch = GRN_list[[s]]
        N_re = base::nrow(current_branch)
        res = base::rownames(current_branch)
        for (TF_TG in res) {
          if (current_branch[TF_TG, "TF"] %in% TFs && current_branch[TF_TG, "TG"] %in% TGs) {
            # Fill in TFE, TGA, TGE values for all time points.
            for (K in 1:ncol(TFE_T[[n]])) {
              current_branch[TF_TG, c(
                base::paste0("TFE_", K), base::paste0("TGA_", K), base::paste0("TGE_", K)
              )] = c(TFE_T[[n]][current_branch[TF_TG, "TF"], K],
                     TGA_T[[n]][current_branch[TF_TG, "TG"], K],
                     TGE_T[[n]][current_branch[TF_TG, "TG"], K])
            }
          } else {
            # Remove pairs that are not in the current branch's gene set.
            current_branch = current_branch[rownames(current_branch) != TF_TG, ]
          }
        }
        return(current_branch)
      }
      data_for_prediction_list = parallel::mclapply(
        X = 1:n_part,
        FUN = function(s) {
          process_GRN_slice(
            s,
            GRN_list = GRN_list,
            TFs = base::colnames(branch_GRN[[n]]),
            TGs = base::rownames(branch_GRN[[n]]),
            TFE_T = interest_cell_type_group[[cell_type]][["Branches_TFE_T"]],
            TGA_T = interest_cell_type_group[[cell_type]][["Branches_TGA_T"]],
            TGE_T = interest_cell_type_group[[cell_type]][["Branches_TGE_T"]]
          )
        },
        mc.cores = ncores
      )
      df = base::do.call("rbind", data_for_prediction_list)

      # Filter out pairs with zero regulatory strength.
      message("Eliminating potential TF-TG pairs with a theta_i of zero.")
      df = df[df$theta_i != 0, , drop = FALSE]
      interest_cell_type_branch_data_for_GRN_infer[[cell_type]][[base::paste0("branch", n)]] = df
      message("In branch ", n, ", there are ", base::nrow(df), " potential regulatory relationships after filtering.")
    }
    # Process global GRN set (filter zero theta pairs).
    df_all = all_combinations[["GRN"]]
    df_all = df_all[df_all$theta_i != 0, , drop = FALSE]
    interest_cell_type_branch_data_for_GRN_infer[[cell_type]][["GRN"]] = df_all
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Getting data for inferring GRN for each branch. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_branch_data_for_GRN_infer)
}

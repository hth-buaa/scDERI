#' @title Extract Branch-Specific Initial Gene Regulatory Networks (iGRNs)
#'
#' @description
#' This function extracts branch-specific initial gene regulatory networks (iGRNs) from the complete cell-type-specific iGRNs.
#' It creates subnetworks for each branch within each cell type by filtering the global TF-TG association matrix to include only the TFs and TGs present in each specific branch.
#'
#' @param interest_cell_type_iGRN The output of function \code{\link{get_iGRN_by_TFBS_pwm_by_JASPAR2024}} or \code{\link{filter_iGRN_by_ChIP_Atlas_priors}}.
#' @param interest_cell_type_group The output of function \code{\link{cell_grouping}}.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#'
#' @return
#' A nested list structure:
#' \itemize{
#'   \item First level: named by cell types
#'   \item Second level: list of branch-specific iGRNs for each cell type
#' }
#' Each branch iGRN is a data frame with target genes as rows and transcription factors as columns.
#' Values are normalized regulatory strengths inherited from the complete network.
#' All-zero rows (target genes) and all-zero columns (TFs) are automatically removed from the output.
#'
#' @details
#' For each cell type and each branch, the function executes the following steps:
#' \enumerate{
#'   \item Extracts the TF set and TG set of the current branch from the cell grouping result.
#'   \item Initializes an all-zero subnetwork matrix matching the branch's gene dimensions.
#'   \item Fills the subnetwork with regulatory strength values from the complete cell-type iGRN for matched TF-TG pairs.
#'   \item Removes all-zero rows and all-zero columns from the subnetwork.
#'   \item Returns organized branch-specific networks for downstream branch-level analysis.
#' }
#'
#' @note
#' \itemize{
#'   \item This function does not recalculate regulatory strengths; it only extracts and filters existing network values.
#'   \item The function creates (if missing) and temporarily switches to the \code{3 get iGRN/} directory during execution,
#'     and restores the original working directory on completion.
#'   \item Parallel processing is performed at the cell type level to improve efficiency for multi-cell-type inputs.
#'   \item The sparsity of branch networks depends on the overlap between branch-specific genes and the complete network.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
#' interest_cell_type_group = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_branch_iGRN = get_branch_iGRN(
#'   interest_cell_type_iGRN = interest_cell_type_iGRN,
#'   interest_cell_type_group = interest_cell_type_group,
#'   ncores = ncores
#' )
#' base::saveRDS(
#'   interest_cell_type_branch_iGRN,
#'   file = "./3 get iGRN/interest_cell_type_branch_iGRN.rds"
#' )
#' }
get_branch_iGRN = function(
    interest_cell_type_iGRN,
    interest_cell_type_group,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Extracting iGRN for each branch from every cell type's iGRN. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Defining inner function.
  get_branch_iGRN_process = function(cell_type, interest_cell_type_iGRN, interest_cell_type_group)
  {
    branch_iGRN = base::list()
    for (n in 1:base::length(interest_cell_type_group[[cell_type]][["n_min"]])) {
      branch_TGs = base::rownames(interest_cell_type_group[[cell_type]][["Branches_TGE_T"]][[n]])
      branch_TFs = base::rownames(interest_cell_type_group[[cell_type]][["Branches_TFE_T"]][[n]])
      branch_iGRN[[n]] = base::as.data.frame(base::matrix(0, nrow = base::length(branch_TGs), ncol = base::length(branch_TFs)))
      base::rownames(branch_iGRN[[n]]) = branch_TGs
      base::colnames(branch_iGRN[[n]]) = branch_TFs
      for (TG in branch_TGs) {
        for (TF in branch_TFs) {
          if (TG %in% base::rownames(interest_cell_type_iGRN[[cell_type]]) &&
              TF %in% base::colnames(interest_cell_type_iGRN[[cell_type]])) {
            branch_iGRN[[n]][TG, TF] = interest_cell_type_iGRN[[cell_type]][TG, TF]
          }
        }
      }
      branch_iGRN[[n]] = branch_iGRN[[n]][
        base::rowSums(branch_iGRN[[n]]) > 0,
        base::colSums(branch_iGRN[[n]]) > 0,
        drop = FALSE
      ]
    }
    return(branch_iGRN)
  }

  ### Invoking functions in parallel.
  interest_cell_type_branch_iGRN = parallel::mclapply(
    X = base::names(interest_cell_type_iGRN),
    FUN = function(cell_type) {
      message("Preparing to extract iGRN for each branch in cell type ", cell_type, ".")
      get_branch_iGRN_process(
        cell_type = cell_type,
        interest_cell_type_iGRN = interest_cell_type_iGRN[cell_type],
        interest_cell_type_group = interest_cell_type_group[cell_type]
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_branch_iGRN) = base::names(interest_cell_type_iGRN)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Extracting iGRN for each branch from every cell type's iGRN. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_branch_iGRN)
}

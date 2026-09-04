#' @title Consolidate Cell Type-Specific Initial Gene Regulatory Networks (iGRNs)
#'
#' @description
#' This function consolidates branch-specific initial gene regulatory networks (iGRNs) into comprehensive cell type-specific iGRNs.
#' It aggregates information from all branches within each cell type to create unified TF-TG association matrices that represent the complete regulatory landscape for each cell type.
#'
#' @param interest_cell_type_iGRN_all_TGTF_pairs The output of function \code{get_iGRN_by_TFBS_pwm_by_JASPAR2024}.
#' @param interest_cell_type_branch_iGRN The output of function \code{get_branch_iGRN}.
#'
#' @returns
#' A named list where each element corresponds to a cell type and contains:
#' \itemize{
#' \item A data frame representing the consolidated iGRN for the cell type
#' \item Rows: All target genes (TGs) found across all branches of the cell type
#' \item Columns: All transcription factors (TFs) found across all branches of the cell type
#' \item Values: Regulatory strength scores from the complete iGRN (values range 0-1)
#' \item Matrix is sparse with many zeros indicating no regulatory relationship
#' }
#' The consolidation process preserves the original regulatory strength values from the complete iGRN while ensuring all TFs and TGs from all branches are represented.
#'
#' @details
#' The function performs the following operations for each cell type:
#' \enumerate{
#' \item Identifies all unique TGs and TFs across all branches of the cell type
#' \item Creates a unified matrix with dimensions encompassing all identified TGs and TFs
#' \item Populates the matrix with regulatory strength values from the complete iGRN
#' \item Sets values to 0 where no regulatory relationship exists in the complete iGRN
#' \item Returns organized cell type-specific networks for comprehensive analysis
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item This function serves as an aggregation step, creating comprehensive cell type-level networks
#' \item The output matrices include all regulatory relationships from the complete iGRN
#' \item Missing values in the complete iGRN are represented as zeros in the output
#' \item The function does not perform any additional network inference or filtering
#' \item Output directory structure is maintained for consistency with other functions
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_iGRN_all_TGTF_pairs = base::readRDS(
#'   "./3 get iGRN/interest_cell_type_iGRN_all_TGTF_pairs.rds"
#' )
#' interest_cell_type_branch_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_branch_iGRN.rds")
#' interest_cell_type_iGRN = get_iGRN(
#'   interest_cell_type_iGRN_all_TGTF_pairs = interest_cell_type_iGRN_all_TGTF_pairs,
#'   interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN
#' )
#' base::saveRDS(interest_cell_type_iGRN, file = "./3 get iGRN/interest_cell_type_iGRN.rds")
#' }
get_iGRN = function(
    interest_cell_type_iGRN_all_TGTF_pairs = interest_cell_type_iGRN_all_TGTF_pairs,
    interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Extracting branch-specific iGRN from each cell type's iGRN ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Do.
  interest_cell_type_iGRN = base::list()
  for (cell_type in base::names(interest_cell_type_branch_iGRN)) {
    message("Extracting branch-specific iGRN from iGRN of ", cell_type, ".")
    TG_cell_type = base::sort(base::unique(base::unlist(base::lapply(
      interest_cell_type_branch_iGRN[[cell_type]],
      base::rownames)
    )))
    TF_cell_type = base::sort(base::unique(base::unlist(base::lapply(
      interest_cell_type_branch_iGRN[[cell_type]],
      base::colnames)
    )))
    interest_cell_type_iGRN[[cell_type]] = base::as.data.frame(base::matrix(
      0, nrow = base::length(TG_cell_type), ncol = base::length(TF_cell_type),
      dimnames = base::list(TG_cell_type, TF_cell_type)
    ))
    for (j in TG_cell_type) {
      for (i in TF_cell_type) {
        if (base::is.null(interest_cell_type_iGRN_all_TGTF_pairs[[cell_type]][j, i])) {
          interest_cell_type_iGRN[[cell_type]][j, i] = 0
        } else {
          interest_cell_type_iGRN[[cell_type]][j, i] = interest_cell_type_iGRN_all_TGTF_pairs[[cell_type]][j, i]
        }
      }
    }
    interest_cell_type_iGRN[[cell_type]] = interest_cell_type_iGRN[[cell_type]][
      base::rowSums(interest_cell_type_iGRN[[cell_type]] != 0) > 0,
      base::colSums(interest_cell_type_iGRN[[cell_type]] != 0) > 0,
      drop = FALSE
    ]
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Extracting branch-specific iGRN from each cell type's iGRN ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_iGRN)
}

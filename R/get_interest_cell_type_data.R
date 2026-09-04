#' @title Process Cell-Type-Specific Single-Cell Multi-Omics Data
#'
#' @description
#' Extracts scRNA-seq and scATAC-seq data for specified cell types,
#' performs quality control by filtering genes and cells based on zero-expression rates,
#' and aligns the two omics datasets to ensure full consistency.
#' The output is formatted for downstream GRN construction and pseudotime analysis.
#'
#' @param alpha_gene Numeric value between 0 and 1. Zero-expression rate threshold for gene filtering.
#'   Genes with a zero-expression proportion higher than this threshold will be removed.
#'   Defaults to \code{1} (no extra filtering beyond all-zero gene removal).
#' @param alpha_cell Numeric value between 0 and 1. Zero-expression rate threshold for cell filtering.
#'   Cells with a zero-expression proportion higher than this threshold will be removed.
#'   Defaults to \code{1} (no extra filtering beyond all-zero cell removal).
#' @param interest_cell_type Character vector. Target cell type names to extract.
#'   Must match values in the \code{seurat_annotations} metadata column of the input Seurat objects.
#' @param Basic_Info List. Output from the \code{\link{get_expression_and_activity_matrix}} function,
#'   which must contain \code{scRNAseq_for_GRN}, \code{scATACseq_for_GRN}, \code{gene_for_GRN},
#'   \code{TGs} and \code{TFs} elements.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#'
#' @return A named list with one entry per input cell type (named by cell type).
#'   Each entry contains:
#'   \itemize{
#'     \item \code{interest_scRNAseq_for_GRN}: Filtered scRNA-seq Seurat object
#'     \item \code{interest_scATACseq_for_GRN}: Filtered scATAC-seq Seurat object, with genes and cells fully aligned to the scRNA-seq data
#'     \item \code{interest_gene_for_GRN}: Gene set for GRN construction (union of filtered transcription factors and target genes)
#'     \item \code{interest_TGs}: Filtered target genes retained in the current cell type
#'     \item \code{interest_TFs}: Filtered transcription factors retained in the current cell type
#'   }
#'
#' @details
#' For each specified cell type, the function executes the following workflow:
#' \enumerate{
#'   \item Subsets the scRNA-seq object by target cell type, and removes genes and cells with all-zero counts
#'   \item Further filters genes and cells using the \code{alpha_gene} and \code{alpha_cell} zero-expression rate thresholds
#'   \item Subsets the scATAC-seq object by the same cell type, and aligns its gene and cell dimensions to the filtered scRNA-seq data
#'   \item Intersects the input gene sets (full gene list, TFs, TGs) with retained genes in the scRNA-seq data, and defines the final gene set as the union of filtered TFs and TGs
#'   \item Trims both omics datasets to the final gene set, verifies complete consistency of row and column names, and returns the processed data
#' }
#'
#' @note
#' \itemize{
#'   \item Both \code{scRNAseq_for_GRN} and \code{scATACseq_for_GRN} must contain a \code{seurat_annotations} metadata column for cell type annotation.
#'   \item This function creates a subdirectory named \code{2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment} in the current working directory, temporarily switches the working directory into it, and restores the original directory on exit.
#'   \item Parallel computation relies on \code{\link[parallel]{mclapply}}, which does not support multicore on Windows.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type = c("CD4 TEM", "cDc")
#' # interest_cell_type = "Treg"
#' Basic_Info = base::readRDS(
#'   "./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.rds"
#' )
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_data = get_interest_cell_type_data(
#'   interest_cell_type = interest_cell_type,
#'   Basic_Info = Basic_Info,
#'   ncores = ncores
#' )
#' base::saveRDS(
#'   interest_cell_type_data,
#'   file = paste0(
#'     "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'     "/interest_cell_type_data.rds"
#'   )
#' )
#' }
get_interest_cell_type_data = function(
    alpha_gene = 1,
    alpha_cell = 1,
    interest_cell_type = interest_cell_type,
    Basic_Info = Basic_Info,
    ncores = 1
)
{
  ### Start.
  if (!all(interest_cell_type %in% Basic_Info$scRNAseq_for_GRN$seurat_annotations)) {
    base::stop("Some cell types in interest_cell_type are not found in seurat_annotations.")
  }
  if (alpha_gene < 0 || alpha_gene > 1 || alpha_cell < 0 || alpha_cell > 1) {
    base::stop("alpha_gene and alpha_cell must be between 0 and 1.")
  }

  t_start = base::Sys.time()
  message("Run: Extracting Data for Cell Types of Interest. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Inner function.
  get_interest_cell_type_data_process = function(alpha_gene = 1,
                                                 alpha_cell = 1,
                                                 cell_type, scRNAseq_for_GRN, scATACseq_for_GRN, gene_for_GRN, TGs, TFs)
  {
    message("The cell label is ", cell_type, ".")

    # Extracting expression matrix for cell types of interest (retaining genes/cells with non-zero expression values).
    message("Extracting expression matrix for cell types of interest (retaining genes/cells with non-zero expression values).")
    interest_scRNAseq_for_GRN = base::subset(scRNAseq_for_GRN, subset = seurat_annotations == cell_type)
    counts_mat <- SeuratObject::GetAssayData(interest_scRNAseq_for_GRN, layer = "counts")
    interest_scRNAseq_for_GRN = interest_scRNAseq_for_GRN[Matrix::rowSums(counts_mat) > 0, ]
    counts_mat <- SeuratObject::GetAssayData(interest_scRNAseq_for_GRN, layer = "counts")
    interest_scRNAseq_for_GRN = interest_scRNAseq_for_GRN[, Matrix::colSums(counts_mat) > 0]

    # Retaining genes with zero-expression rate <= alpha_gene and cells with zero-expression rate alpha_cell.
    message("Retaining genes with zero-expression rate <= ", alpha_gene,
            " and cells with zero-expression rate <= ", alpha_cell, ".")
    zero_frequency =
      base::rowMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(interest_scRNAseq_for_GRN, layer = "counts"))) == 0)
    interest_scRNAseq_for_GRN =
      base::subset(interest_scRNAseq_for_GRN, features = base::names(zero_frequency[zero_frequency <= alpha_gene]))

    zero_frequency =
      base::colMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(interest_scRNAseq_for_GRN, layer = "counts"))) == 0)
    interest_scRNAseq_for_GRN =
      base::subset(interest_scRNAseq_for_GRN, cells = base::names(zero_frequency[zero_frequency <= alpha_cell]))

    # Extracting and aligning activity matrix to match the filtered expression matrix dimensions.
    message("Extracting and aligning activity matrix to match the filtered expression matrix dimensions.")
    interest_scATACseq_for_GRN = base::subset(scATACseq_for_GRN, subset = seurat_annotations == cell_type)
    interest_scATACseq_for_GRN =
      interest_scATACseq_for_GRN[base::rownames(interest_scRNAseq_for_GRN), base::colnames(interest_scRNAseq_for_GRN)]

    # Finalizing gene sets and ensuring consistency between RNA and ATAC data.
    message("Finalizing gene sets and ensuring consistency between RNA and ATAC data.")
    interest_gene_for_GRN = base::intersect(gene_for_GRN, base::rownames(interest_scRNAseq_for_GRN))
    interest_TGs = base::sort(base::intersect(TGs, interest_gene_for_GRN))
    interest_TFs = base::sort(base::intersect(TFs, interest_gene_for_GRN))
    interest_gene_for_GRN = base::sort(base::union(interest_TGs, interest_TFs))
    interest_scRNAseq_for_GRN =
      interest_scRNAseq_for_GRN[base::sort(base::intersect(base::rownames(interest_scRNAseq_for_GRN), interest_gene_for_GRN)), ]
    interest_scATACseq_for_GRN = interest_scATACseq_for_GRN[base::rownames(interest_scRNAseq_for_GRN), base::colnames(interest_scRNAseq_for_GRN)]

    base::stopifnot(base::identical(
      base::rownames(interest_scRNAseq_for_GRN),
      base::rownames(interest_scATACseq_for_GRN)
    ))
    base::stopifnot(
      base::identical(base::sort(base::rownames(interest_scRNAseq_for_GRN)), base::sort(base::rownames(interest_scATACseq_for_GRN))),
      base::identical(base::sort(base::colnames(interest_scRNAseq_for_GRN)), base::sort(base::colnames(interest_scATACseq_for_GRN)))
    )

    # End.
    return(base::list("interest_scRNAseq_for_GRN" = interest_scRNAseq_for_GRN,
                      "interest_scATACseq_for_GRN" = interest_scATACseq_for_GRN,
                      "interest_gene_for_GRN" = interest_gene_for_GRN,
                      "interest_TGs" = interest_TGs,
                      "interest_TFs" = interest_TFs))
  }

  ### Iterating inner function.
  interest_cell_type_data = parallel::mclapply(
    X = interest_cell_type,
    FUN = function(cell_type) {get_interest_cell_type_data_process(
      alpha_gene = alpha_gene,
      alpha_cell = alpha_cell,
      cell_type = cell_type,
      scRNAseq_for_GRN = Basic_Info[["scRNAseq_for_GRN"]],
      scATACseq_for_GRN = Basic_Info[["scATACseq_for_GRN"]],
      gene_for_GRN = Basic_Info[["gene_for_GRN"]],
      TGs = Basic_Info[["TGs"]],
      TFs = Basic_Info[["TFs"]]
    )},
    mc.cores = ncores)
  base::names(interest_cell_type_data) = interest_cell_type

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd())
  t_end = base::Sys.time()
  message("Finish: Extracting Data for Cell Types of Interest. Time: ", t_end)
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_data)
}

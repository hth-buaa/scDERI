#' @title Obtain Expression and Gene Activity Matrices for GRN Construction
#'
#' @description
#' Processes single-cell multi-omics data to extract expression matrices of target genes (TGs) and transcription factors (TFs),
#' and computes a gene activity matrix from scATAC-seq data for downstream Gene Regulatory Network (GRN) construction.
#'
#' @param TGs Character vector of target gene symbols. Defaults to \code{results_identify_TGs$TGs}.
#' \code{results_identify_TGs} is the output of function \code{\link{identify_TGs}} or \code{\link{filter_peak_annotations}}.
#' @param TFs Character vector of transcription factor symbols, output from \code{\link{get_TFs_from_priors}}.
#' @param scRNAseq Seurat object of single-cell RNA-seq data. See \code{\link{identify_TGs}} for details.
#' @param scATACseq Seurat object of single-cell ATAC-seq data. See \code{\link{identify_TGs}} for details.
#' @param Ensdb \code{EnsDb} genomic annotation object, passed to \code{\link[Signac]{GetGRangesFromEnsDb}}.
#' Must match the species and genome build of input data. Defaults to \code{EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86}.
#' @param flanking_region Positive numeric.
#' Base pairs to extend for gene activity calculation: applied both upstream of the transcription start site (TSS) and downstream of the transcription termination site (TTS).
#' Passed to both \code{extend.upstream} and \code{extend.downstream} in \code{\link[Signac]{GeneActivity}}.
#' Defaults to 10000.
#' @param path Absolute path to the scATAC-seq fragment file (\code{.tsv.gz}).
#' A matching \code{.tbi} index file must exist in the same directory.
#' For example, \code{"/mnt/d/pbmc-equation/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"}.
#'
#' @returns
#' A list with 5 elements:
#' \itemize{
#'   \item \code{scRNAseq_for_GRN}: Seurat object of scRNA-seq data, subsetted to genes shared with the scATAC-seq dataset
#'   \item \code{scATACseq_for_GRN}: Seurat object of scATAC-seq data with a \code{ACTIVITY} assay, subsetted to shared genes;
#'     \code{ACTIVITY} is set as the default assay
#'   \item \code{gene_for_GRN}: Character vector of genes present in both scRNA-seq and scATAC-seq datasets
#'   \item \code{TGs}: Filtered target genes retained in both omics datasets
#'   \item \code{TFs}: Filtered transcription factors retained in both omics datasets
#' }
#'
#' @details
#' The function executes the following main steps:
#' \enumerate{
#'   \item Creates a dedicated output directory and temporarily switches the working directory (restored automatically on exit)
#'   \item Combines TGs and TFs into a unified gene set, and extracts the corresponding expression subset from scRNA-seq data
#'   \item Imports genomic annotations from EnsDb, links the fragment file to the scATAC-seq object, and computes the gene activity matrix
#'   \item Identifies genes shared by both omics datasets, and filters all matrices and gene lists to maintain consistency
#'   \item Returns standardized multi-omics data ready for GRN inference
#' }
#'
#' @note
#' \itemize{
#'   \item The input fragment file must be paired with a \code{.tbi} index file in the same directory
#'   \item Gene activity scores are calculated from ATAC fragment counts across gene bodies plus upstream TSS and downstream TTS flanking regions
#'   \item Stable operation is only guaranteed on Linux. Genome coordinate style conversion via \code{GenomeInfoDb} has high failure rate on Windows
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run in linux.
#' library(EnsDb.Hsapiens.v86)
#' results_identify_TGs = base::readRDS(
#'   "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds"
#' )
#' TFs = readRDS("./1.1 Data Preprocessing - Get TFs From Database/TFs.rds")
#' scRNAseq = base::readRDS("./scRNAseq.rds")
#' scATACseq = base::readRDS("./scATACseq.rds")
#' flanking_region = 10000
#' Basic_Info = get_expression_and_activity_matrix(
#'   TGs = results_identify_TGs$TGs,
#'   TFs = TFs,
#'   scRNAseq = scRNAseq,
#'   scATACseq = scATACseq,
#'   Ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
#'   flanking_region = flanking_region,
#'   path = "/mnt/d/pbmc-equation/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
#' )
#' base::saveRDS(
#'   Basic_Info,
#'   file = "./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.rds"
#' )
#' }
get_expression_and_activity_matrix = function(
    TGs = results_identify_TGs$TGs,
    TFs,
    scRNAseq,
    scATACseq,
    Ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
    flanking_region = 10000,
    path
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Obtaining Counts Matrices for Gene Expression and Activity. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "1.3 Data Preprocessing - Get Expression And Activity Matrix"
  if (!base::dir.exists(new_folder)) {base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Extracting gene expression matrix.
  gene_for_GRN = base::union(TGs, TFs)
  message("Extracting gene expression matrix.")
  scRNAseq_for_GRN = base::subset(scRNAseq, features = gene_for_GRN)

  ### Generating gene activity matrix (based on ATAC-seq fragment counts in gene bodies and flanking regions of TSS per cell)
  message("Annotating chromatin fragments in scATAC-seq data with genomic annotation.")
  annotations = Signac::GetGRangesFromEnsDb(ensdb = Ensdb)
  GenomeInfoDb::seqlevelsStyle(annotations) = 'UCSC'
  Signac::Annotation(scATACseq) = annotations

  message("Adding fragment file. Both the fragment file and its .tbi index must exist in the directory specified by 'path'.")
  fragment = Signac::CreateFragmentObject(path = path)
  fragment@cells = SeuratObject::Cells(scATACseq)
  Signac::Fragments(scATACseq) = NULL
  Signac::Fragments(scATACseq) = fragment

  message("Adding an 'ACTIVITY' matrix to scATAC-seq data in Seurat object.")
  gene.activities = Signac::GeneActivity(scATACseq,
                                         features = base::rownames(scRNAseq),
                                         extend.upstream = flanking_region,
                                         extend.downstream = flanking_region)
  scATACseq[["ACTIVITY"]] = SeuratObject::CreateAssayObject(counts = gene.activities)
  SeuratObject::DefaultAssay(scATACseq) = "ACTIVITY"
  scATACseq_for_GRN = base::subset(scATACseq, features = gene_for_GRN)

  ### Update data.
  message("Update data.")
  gene_for_GRN = base::intersect(base::rownames(scRNAseq_for_GRN), base::rownames(scATACseq_for_GRN))
  TGs = base::intersect(TGs, gene_for_GRN)
  TFs = base::intersect(TFs, gene_for_GRN)
  scRNAseq_for_GRN = base::subset(scRNAseq, features = gene_for_GRN)
  scATACseq_for_GRN = base::subset(scATACseq, features = gene_for_GRN)
  SeuratObject::DefaultAssay(scATACseq_for_GRN) = "ACTIVITY"

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Obtaining Counts Matrices for Gene Expression and Activity. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(base::list("scRNAseq_for_GRN" = scRNAseq_for_GRN,
                    "scATACseq_for_GRN" = scATACseq_for_GRN,
                    "gene_for_GRN" = gene_for_GRN,
                    "TGs" = TGs,
                    "TFs" = TFs))
}

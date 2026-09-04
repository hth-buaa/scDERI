#' @title Annotate ATAC Peaks and Identify Promoter-associated Target Genes
#'
#' @description Performs genomic annotation on scATAC-seq chromatin accessibility peaks,
#' identifies target genes associated with open promoter regions, and retains only genes
#' detected in the matched scRNA-seq data. Generates annotation tables and visualization files.
#'
#' @param genome_for_anno A \code{TxDb} object for genome annotation.
#' Defaults to \code{TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene}.
#' @param fragment_annotation_range Positive numeric. Annotation range around TSS in base pairs,
#' used to define promoter-proximal regions. Defaults to 500000.
#' Smaller values reduce runtime but may exclude distal regulatory elements.
#' @param scRNAseq
#' A Seurat object containing scRNA-seq data (the initially loaded data) and gene names as rownames.
#' The \code{scRNAseq@meta.data} contains the \code{"seurat_annotation"} column, which describes cell labels.
#' Used to filter target genes detectable in transcriptome data.
#' @param scATACseq
#' A Seurat object containing scATAC-seq data (the initially loaded data) with files \code{_fragments.tsv.gz.tbi} and \code{_fragments.tsv.gz} (in \code{Signac} framework),
#' whose peak ranges are used for annotation.
#' The \code{scATACseq@meta.data} contains the seurat_annotation column, which describes cell labels.
#' @param annoDb An \code{OrgDb} object for organism-level gene ID conversion.
#' Defaults to \code{org.Hs.eg.db::org.Hs.eg.db}.
#' @param flankDistance Non-negative numeric. Distance (bp) extended upstream and downstream
#' of each peak to search for flanking genes. Defaults to 10000.
#' See \code{\link[ChIPseeker]{annotatePeak}} for implementation details.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{peak_anno_promoter}{Data frame of annotated peaks located in promoter regions,
#'     with added flanking gene symbols, and all associated genes present in scRNA-seq.}
#'   \item{TGs}{Character vector of unique target gene symbols.}
#' }
#'
#' @details The function follows a 4-step workflow:
#' \enumerate{
#'   \item \strong{Peak annotation}: Annotate all ATAC peaks via \code{\link[ChIPseeker]{annotatePeak}}
#'     with the specified TSS range, and append flanking gene information.
#'   \item \strong{Quality filtering}: Remove peaks with missing ENSEMBL ID, gene symbol or annotation.
#'     Save the full annotation table and a genomic feature distribution bar plot.
#'   \item \strong{Promoter peak extraction}: Subset peaks annotated as promoter regions.
#'     Parse flanking transcript IDs, remove version suffixes, and convert to official gene symbols.
#'   \item \strong{Target gene filtering}: Intersect promoter-associated genes with genes in scRNA-seq
#'     as the final target gene set, and subset the promoter peak table accordingly.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates an output directory and temporarily switches the working directory
#'     during execution, and restores the original path upon completion.
#'   \item Output files (\code{peak_anno.png}, \code{peak_anno.csv}, \code{TGs.csv}) are saved
#'     under the \code{1.2 Data Preprocessing - Identify TGs By Annotation} directory.
#'   \item The term "promoter" here is defined broadly by the TSS range set via \code{fragment_annotation_range}.
#'   \item The scRNA-seq data and scATAC-seq data need to be loaded as objects named scRNAseq and scATACseq,
#'    respectively, and they must come from the same cells with identical cell labels.
#'    This is a strict input requirement for scDERI, although the function will still run even if it is not met.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' # Load scRNA-seq data and scATAC-seq data
#' # as object scRNAseq and scATACseq with the same cells and cell labels.
#' # scRNAseq = SeuratData::LoadData("pbmcMultiome", "pbmc.rna")
#' # scATACseq = SeuratData::LoadData("pbmcMultiome", "pbmc.atac")
#' scRNAseq = base::readRDS("./scRNAseq.rds")
#' scATACseq = base::readRDS("./scATACseq.rds")
#' addFlankGeneInfo = base::min(10000, fragment_annotation_range)
#' results_identify_TGs = identify_TGs(
#'   genome_for_anno = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   fragment_annotation_range = 500000,
#'   scRNAseq = scRNAseq,
#'   scATACseq = scATACseq,
#'   annoDb = org.Hs.eg.db::org.Hs.eg.db,
#'   flankDistance = 10000
#' )
#' base::saveRDS(
#'   results_identify_TGs,
#'   file = "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds"
#' )
#' }
identify_TGs = function(genome_for_anno = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                        fragment_annotation_range = 500000,
                        scRNAseq,
                        scATACseq,
                        annoDb = org.Hs.eg.db::org.Hs.eg.db,
                        flankDistance = 10000)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Annotating Chromatin Fragments and Identifying TGs. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "1.2 Data Preprocessing - Identify TGs By Annotation"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  message("Annotate chromatin fragments using the annotatePeak function from the ChIPseeker package.")
  message("The genome annotation database used is ", GenomeInfoDb::genome(genome_for_anno)[1], ".")

  ### Peak annotation.
  peak_anno = ChIPseeker::annotatePeak(
    peak = scATACseq@assays[["ATAC"]]@ranges,
    tssRegion = c(-fragment_annotation_range, fragment_annotation_range),
    TxDb = genome_for_anno,
    annoDb = annoDb,
    assignGenomicAnnotation = TRUE,
    addFlankGeneInfo = TRUE,
    flankDistance = flankDistance,
    sameStrand = FALSE,
    ignoreOverlap = TRUE
  )

  ### Savethe annotation results.
  message("Save the annotation results.")
  ggplot2::ggsave("peak_anno.png", plot = ChIPseeker::plotAnnoBar(peak_anno), width = 10, height = 6)
  peak_anno = base::as.data.frame(peak_anno) %>% dplyr::filter(!is.na(ENSEMBL)) %>% dplyr::filter(!is.na(SYMBOL)) %>% dplyr::filter(!is.na(annotation))
  utils::write.csv(peak_anno, file = "peak_anno.csv", row.names = FALSE)
  base::rownames(peak_anno) = IRanges::paste(peak_anno$seqnames, peak_anno$start, peak_anno$end, sep = "-")

  ### Obtain TGs.
  message("Identify target genes and process flank gene annotations.")
  peak_anno_promoter = peak_anno[base::grepl("promoter", peak_anno$annotation, ignore.case = TRUE), ]
  all_tx_ids = base::unique(base::unlist(base::strsplit(peak_anno_promoter$flank_txIds, ";")))
  all_tx_ids_clean = base::sub("\\.[0-9]+$", "", all_tx_ids)
  tx_to_symbol = clusterProfiler::bitr(
    all_tx_ids_clean,
    fromType = "ENSEMBLTRANS",
    toType = "SYMBOL",
    OrgDb = annoDb
  )
  peak_anno_promoter$flank_geneSymbols = base::sapply(
    peak_anno_promoter$flank_txIds,
    function(tx_string) {
      if (base::is.na(tx_string)) return(NA)
      tx_list = base::unlist(base::strsplit(tx_string, ";"))
      tx_list_clean = base::sub("\\.[0-9]+$", "", tx_list)
      symbols = base::unique(tx_to_symbol$SYMBOL[base::match(tx_list_clean, tx_to_symbol$ENSEMBLTRANS)])
      symbols = symbols[!base::is.na(symbols)]
      base::paste(symbols, collapse = ";")
    }
  )
  TGs = base::unique(base::intersect(peak_anno_promoter$SYMBOL, base::rownames(scRNAseq)))
  message(
    "Total promoter-associated genes: ", base::length(base::unique(peak_anno_promoter$SYMBOL)),
    "; retained in scRNA-seq: ", base::length(TGs)
  )
  # tryCatch(
  #   utils::write.csv(TGs, file = "TGs.csv", row.names = FALSE), error = function(e) {
  #     stop("The file could not be saved: ", e$message)
  #   }
  # )
  tryCatch({
    utils::write.csv(data.frame(target_gene = TGs), file = "TGs.csv", row.names = FALSE)
  }, error = function(e) {
    warning("Failed to save TGs.csv: ", e$message)
  })

  ### Filter promoter annotation to retain genes expressed in scRNA-seq
  peak_anno_promoter = peak_anno_promoter %>% dplyr::filter(SYMBOL %in% TGs)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Annotating Chromatin Fragments and Identifying TGs. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(base::list("peak_anno_promoter" = peak_anno_promoter, "TGs" = TGs))
}

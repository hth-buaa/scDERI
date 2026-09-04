#' @title Infer Predictive Enhancer-driven Gene Regulatory Networks (peGRNs)
#'
#' @description
#' Integrates initial gene regulatory networks (iGRNs), predictive gene regulatory networks (pGRNs),
#' and single-cell multi-omics features to construct peak-level enhancer-driven gene regulatory networks (peGRNs).
#' Returns results at two confidence levels (loose and strict) for each cell type.
#'
#' @param scRNAseq Seurat object of scRNA-seq gene expression data, used to compute gene-level expression features. See also \code{\link{identify_TGs}}.
#' @param scATACseq Seurat object of scATAC-seq chromatin accessibility data, used to compute peak-level accessibility features.
#' @param results_identify_TGs List output from \code{\link{identify_TGs}} or \code{\link{filter_peak_annotations}}.
#'   Must contain the \code{peak_anno_promoter} element with peak-to-gene promoter mapping (columns \code{peak} and \code{SYMBOL}).
#' @param interest_cell_type_iGRN Output from \code{\link{get_iGRN_by_TFBS_pwm_by_JASPAR2024}} or \code{\link{filter_iGRN_by_ChIP_Atlas_priors}}.
#'   A nested list with cell types as top-level names; each element is a TF-TG regulatory strength matrix (TGs as rows, TFs as columns).
#' @param interest_cell_type_pGRN Output from \code{\link{filter_GRN}}.
#'   A nested list with cell types as top-level names; each cell type contains \code{pred_summary_loose} and
#'   \code{pred_summary_strict} data frames with predicted TF-TG regulatory strengths and fitting losses.
#' @param FPS_factor Numeric between 0 and 1. Expansion factor for per-TF min-max scaling of PWM match scores.
#'   Widens the scaling range based on global score boundaries to reduce boundary saturation. Defaults to \code{1.0}.
#'
#' @returns
#' A nested list organized by cell type. Each cell type contains two data frames of TF-peak-TG regulatory edges:
#' \itemize{
#'   \item \code{loose}: Edges at loose confidence level
#'   \item \code{strict}: Edges at strict confidence level
#' }
#'
#' Each row represents one TF-peak-TG edge, with columns grouped as follows:
#' \enumerate{
#'   \item Core identifiers: \code{TF}, \code{peak}, \code{TG}
#'   \item TF-peak binding features: Per-TF min-max scaled PWM match score (\code{FP-S}) and its normalized rank (\code{FP-S-R})
#'   \item TF-TG regulatory features:
#'     Predicted strength from pGRN (\code{FG-theta_p}),
#'     raw fitting loss (\code{FG-L}) and edge-set normalized fitting loss (\code{FG-NL}),
#'     initial strength from iGRN (\code{FG-theta_i}) and its normalized rank (\code{FG-theta_i-R})
#'   \item Multi-omics features: Non-zero ratio, non-zero mean, overall mean and their normalized ranks,
#'     for TF expression (prefixed with \code{F-}), peak accessibility (prefixed with \code{P-}) and TG expression (prefixed with \code{G-})
#' }
#'
#' @details
#' The function executes in two main stages:
#'
#' **1. Global preprocessing (shared across all cell types)**
#' - Construct peak-to-TG promoter mapping from input annotations, and filter to retain only peaks and genes
#'   present in the corresponding scATAC-seq and scRNA-seq matrices.
#' - Compute peak-level accessibility features (after TF-IDF normalization) and gene-level expression features
#'   (after log-normalization) using all cells. For each feature, both raw values and normalized ranks are calculated.
#'
#' **2. Per-cell-type network construction**
#' - Load precomputed TF-peak PWM matching results for the cell type, and aggregate multiple scores per
#'   TF-peak-TG pair by mean.
#' - Perform per-TF min-max scaling on PWM match scores. The scaling range is expanded by \code{FPS_factor}
#'   using the global maximum/minimum scores to alleviate boundary saturation. Global normalized ranks of
#'   scaled scores are also computed.
#' - For each confidence level:
#'   - Filter TF-peak-TG triples to keep only pairs present in both the corresponding pGRN confidence edge set
#'     and the valid peak-TG promoter set.
#'   - Append TF-TG regulatory metrics: initial strength from iGRN (ranked against all positive iGRN values),
#'     predicted strength and fitting loss from pGRN (loss normalized within the current edge set).
#'   - Attach precomputed omics features for the TF, peak and TG of each edge.
#'   - Sort and reorder columns to produce the final edge table.
#'
#' @note
#' \itemize{
#'   \item Requires precomputed files named \code{<cell_type>_TF_match_peak_results.csv} (TF-peak PWM matching results)
#'     under the \code{3 get iGRN/} working subdirectory (one file per cell type).
#'   \item Internally creates a \code{5 Infer GRN} directory and switches the working directory into it during execution;
#'     the original working directory is restored before function exit.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' setwd(Working_directory)
#' scRNAseq_filtered = readRDS(
#'   "./1.2 Data Preprocessing - Identify TGs By Annotation/scRNAseq_filtered.rds"
#' )
#' scATACseq = readRDS("./scATACseq.rds")
#' results_identify_TGs = readRDS(
#'   "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds"
#' )
#' interest_cell_type_iGRN = readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
#' interest_cell_type_pGRN = readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' interest_cell_type_peGRN = infer_eGRN(
#'   scRNAseq = scRNAseq_filtered,
#'   scATACseq = scATACseq,
#'   results_identify_TGs = results_identify_TGs,
#'   interest_cell_type_iGRN = interest_cell_type_iGRN,
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   FPS_factor = 1.0
#' )
#' saveRDS(interest_cell_type_peGRN, file = "./5 Infer GRN/interest_cell_type_peGRN.rds")
#' }
infer_eGRN = function(
    scRNAseq,
    scATACseq,
    results_identify_TGs,
    interest_cell_type_iGRN,
    interest_cell_type_pGRN,
    FPS_factor = 1.0
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Inferring Predictive Enhancer-driven Gene Regulatory Networks (peGRNs). Time: ", t_start, ".")
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

  ### Pre-check: verify all required TF-peak result files exist.
  for (cell_type in base::names(interest_cell_type_pGRN)) {
    file_path = base::paste0(
      original_dir, "/3 get iGRN/", cell_type, "_TF_match_peak_results.csv"
    )
    if (!base::file.exists(file_path)) {
      base::stop(base::paste0("File not found: ", file_path))
    }
  }

  ### Output list: peGRN results per cell type
  interest_cell_type_peGRN = base::list()

  ### Step 1: Build global peak-TG (promoter) mapping and compute shared omics features.
  # Features are computed using all cells, shared across all cell types and confidence levels.
  message("Building peak-TG promoter mapping.")
  peak_TG = results_identify_TGs[["peak_anno_promoter"]][, c("peak", "SYMBOL"), drop = FALSE]
  base::colnames(peak_TG)[base::colnames(peak_TG) == "SYMBOL"] = "TG"
  # Keep only peaks and genes present in the omics matrices.
  peak_TG = peak_TG[peak_TG$peak %in% base::rownames(scATACseq), ]
  peak_TG = peak_TG[peak_TG$TG %in% base::rownames(scRNAseq), ]
  peak_TG = peak_TG[base::order(peak_TG$peak, peak_TG$TG), ]
  base::rownames(peak_TG) = NULL
  peak_TG[, "peak_TG"] = base::paste0(peak_TG$peak, "_to_", peak_TG$TG)

  # Peak-level features from scATAC-seq.
  message("Obtain peak-level accessibility features.")
  scATACseq = Signac::RunTFIDF(scATACseq)
  atac_data = SeuratObject::GetAssayData(scATACseq, layer = "data")
  nonzero_n_peak = Matrix::rowSums(atac_data > 0)
  peak_features = base::data.frame(
    nonzero_freq = nonzero_n_peak / base::ncol(atac_data),
    nonzero_mean = Matrix::rowSums(atac_data) / nonzero_n_peak,
    overall_mean = Matrix::rowMeans(atac_data)
  )
  # Normalized rank versions of each peak feature.
  peak_features[base::paste0(base::colnames(peak_features), "-norm_rank")] =
    base::lapply(peak_features, function(x) {
      r = base::rank(x, ties.method = "min", na.last = "keep")
      r / base::max(r, na.rm = TRUE)
    })

  # Attach peak features to peak-TG table.
  message("[Peak Features] Obtain Non-Zero Accessibility Ratio (P-NZAR).")
  peak_TG[, "P-NZAR"] = peak_features[peak_TG$peak, "nonzero_freq"]
  message("[Peak Features] Obtain normalized Non-Zero Accessibility Ratio Ranks (P-NZAR-R).")
  peak_TG[, "P-NZAR-R"] = peak_features[peak_TG$peak, "nonzero_freq-norm_rank"]

  message("[Peak Features] Obtain Non-Zero Accessibility Average (P-NZAA).")
  peak_TG[, "P-NZAA"] = peak_features[peak_TG$peak, "nonzero_mean"]
  message("[Peak Features] Obtain normalized Non-Zero Accessibility Average Ranks (P-NZAA-R).")
  peak_TG[, "P-NZAA-R"] = peak_features[peak_TG$peak, "nonzero_mean-norm_rank"]

  message("[Peak Features] Obtain Accessibility Average (P-AA).")
  peak_TG[, "P-AA"] = peak_features[peak_TG$peak, "overall_mean"]
  message("[Peak Features] Obtain normalized Accessibility Average Ranks (P-AA-R).")
  peak_TG[, "P-AA-R"] = peak_features[peak_TG$peak, "overall_mean-norm_rank"]

  # Gene-level features from scRNA-seq (shared by both TFs and TGs).
  message("Obtain gene-level expression features (shared by TFs and TGs).")
  scRNAseq = Seurat::NormalizeData(scRNAseq)
  rna_data = SeuratObject::GetAssayData(scRNAseq, layer = "data")
  nonzero_n_TG = Matrix::rowSums(rna_data > 0)
  gene_features = base::data.frame(
    nonzero_freq = nonzero_n_TG / base::ncol(rna_data),
    nonzero_mean = Matrix::rowSums(rna_data) / nonzero_n_TG,
    overall_mean = Matrix::rowMeans(rna_data)
  )
  # Normalized rank versions of each gene feature.
  gene_features[base::paste0(base::colnames(gene_features), "-norm_rank")] =
    base::lapply(gene_features, function(x) {
      r = base::rank(x, ties.method = "min", na.last = "keep")
      r / base::max(r, na.rm = TRUE)
    })

  # Attach TG features to peak-TG table.
  message("[TG Features] Obtain Non-Zero Expression Ratio (G-NZER).")
  peak_TG[, "G-NZER"] = gene_features[peak_TG$TG, "nonzero_freq"]
  message("[TG Features] Obtain normalized Non-Zero Expression Ratio Ranks (G-NZER-R).")
  peak_TG[, "G-NZER-R"] = gene_features[peak_TG$TG, "nonzero_freq-norm_rank"]

  message("[TG Features] Obtain Non-Zero Expression Average (G-NZEA).")
  peak_TG[, "G-NZEA"] = gene_features[peak_TG$TG, "nonzero_mean"]
  message("[TG Features] Obtain normalized Non-Zero Expression Average Ranks (G-NZEA-R).")
  peak_TG[, "G-NZEA-R"] = gene_features[peak_TG$TG, "nonzero_mean-norm_rank"]

  message("[TG Features] Obtain Expression Average (G-EA).")
  peak_TG[, "G-EA"] = gene_features[peak_TG$TG, "overall_mean"]
  message("[TG Features] Obtain normalized Expression Average Ranks (G-EA-R).")
  peak_TG[, "G-EA-R"] = gene_features[peak_TG$TG, "overall_mean-norm_rank"]

  ### Step 2: Build peGRN for each cell type at loose / strict confidence levels.
  for (cell_type in base::names(interest_cell_type_pGRN)) {
    interest_cell_type_peGRN[[cell_type]] = base::list()

    for (level in c("loose", "strict")) {
      message("Obtain ", cell_type, " ", level, " peGRN with features of regulatory components.")

      # Load TF-peak-TG mapping from precomputed PWM-TFBS results
      message("[", cell_type, " ", level, " peGRN] Obtain TF-peak-TG relationships.")
      TF_peak_TG = base::as.data.frame(data.table::fread(base::paste0(
        original_dir, "/3 get iGRN/", cell_type, "_TF_match_peak_results.csv"
      ))[, c("Potential_TFs", "peaks" ,"TGs", "match_score")])
      base::colnames(TF_peak_TG) = c("TF", "peak" ,"TG", "match_score")
      # Aggregate multiple PWM match scores for the same TF-peak-TG pair by taking the mean.
      TF_peak_TG = stats::aggregate(. ~ TF + peak + TG, data = TF_peak_TG, FUN = mean)

      message("[", cell_type, " ", level, " peGRN] Obtain TF-related features.")
      TFs = base::sort(base::unique(TF_peak_TG$TF))
      TGs = base::sort(base::unique(TF_peak_TG$TG))

      # TF-peak PWM match score (FP-S) with per-TF min-max scaling.
      message("[", cell_type, " ", level, " peGRN] [TF-Peak Features] Obtain normalized matching Score between TFBS PWMs and peak bases (FP-S).")
      TF_peak_TG[, "FP-S"] = TF_peak_TG[, "match_score"]
      TF_peak_TG[, "match_score"] = NULL
      TF_S_max = base::max(TF_peak_TG[, "FP-S"])
      TF_S_min = base::min(TF_peak_TG[, "FP-S"])
      # Per-TF min-max scaling of PWM match scores
      # Scaling range is expanded by FPS_factor based on global score range to avoid boundary saturation
      for (tf in TFs) {
        tf_S = TF_peak_TG[TF_peak_TG$TF == tf, "FP-S", drop = FALSE]
        tf_S_max = base::max(tf_S)
        tf_S_min = base::min(tf_S)
        # Expand scaling range slightly to avoid boundary saturation
        scale_max = tf_S_max + (TF_S_max - tf_S_max) * FPS_factor
        scale_min = tf_S_min - (tf_S_min - TF_S_min) * FPS_factor
        tf_S = (tf_S - scale_min) / (scale_max - scale_min)
        TF_peak_TG[TF_peak_TG$TF == tf, "FP-S"] = tf_S
      }
      # Rank-normalized version of FP-S
      message("[", cell_type, " ", level, " peGRN] [TF-Peak Features] Obtain normalized Ranks of normalized matching Score between TF PWMs and peak bases (FP-S-R).")
      r = base::rank(TF_peak_TG[, "FP-S", drop = TRUE], ties.method = "min", na.last = "keep")
      TF_peak_TG[, "FP-S-R"] = r / base::max(r, na.rm = TRUE)

      # Filter edges: retain only TF-TG pairs present in pGRN and peak-TG pairs with valid promoter annotation
      message("[", cell_type, " ", level, " peGRN] Refining TF-peak-TG relationships.")
      TF_peak_TG[, "TF_TG"] = base::paste0(TF_peak_TG$TF, "_to_", TF_peak_TG$TG)
      TF_TG = base::rownames(interest_cell_type_pGRN[[cell_type]][[base::paste0("pred_summary_", level)]])
      TF_peak_TG = TF_peak_TG[TF_peak_TG$TF_TG %in% TF_TG, ]

      TF_peak_TG[, "peak_TG"] = base::paste0(TF_peak_TG$peak, "_to_", TF_peak_TG$TG)
      TF_peak_TG = TF_peak_TG[TF_peak_TG$peak_TG %in% peak_TG$peak_TG, ]

      # TF expression features (same metric set as TG features).
      message("[", cell_type, " ", level, " peGRN] [TF Features] Obtain Non-Zero Expression Ratio (F-NZER).")
      TF_peak_TG[, "F-NZER"] = gene_features[TF_peak_TG$TF, "nonzero_freq"]
      message("[", cell_type, " ", level, " peGRN] [TF Features] Obtain normalized Non-Zero Expression Ratio Ranks (F-NZER-R).")
      TF_peak_TG[, "F-NZER-R"] = gene_features[TF_peak_TG$TF, "nonzero_freq-norm_rank"]

      message("[", cell_type, " ", level, " peGRN] [TF Features] Obtain Non-Zero Expression Average (F-NZEA).")
      TF_peak_TG[, "F-NZEA"] = gene_features[TF_peak_TG$TF, "nonzero_mean"]
      message("[", cell_type, " ", level, " peGRN] [TF Features] Obtain normalized Non-Zero Expression Average Ranks (F-NZEA-R).")
      TF_peak_TG[, "F-NZEA-R"] = gene_features[TF_peak_TG$TF, "nonzero_mean-norm_rank"]

      message("[", cell_type, " ", level, " peGRN] [TF Features] Obtain Expression Average (F-EA).")
      TF_peak_TG[, "F-EA"] = gene_features[TF_peak_TG$TF, "overall_mean"]
      message("[", cell_type, " ", level, " peGRN] [TF Features] Obtain normalized Expression Average Ranks (F-EA-R).")
      TF_peak_TG[, "F-EA-R"] = gene_features[TF_peak_TG$TF, "overall_mean-norm_rank"]

      # TF-TG regulatory features from iGRN (theta_i).
      message("[", cell_type, " ", level, " peGRN] [TF-TG Features] Obtain theta_i (FG-theta_i).")
      TF_peak_TG[, "FG-theta_i"] = 0
      tftg_all = base::unique(TF_peak_TG[, "TF_TG"])
      for (tftg in tftg_all) {
        tf = base::sub("_to_.*", "", tftg)
        tg = base::sub(".*_to_", "", tftg)
        TF_peak_TG[TF_peak_TG$TF_TG == tftg, "FG-theta_i" ] =
          interest_cell_type_iGRN[[cell_type]][tg, tf]
      }
      # Rank normalization only applies to positive theta_i; pairs with theta_i = 0 keep rank 0.
      message("[", cell_type, " ", level, " peGRN] [TF-TG Features] Obtain theta_i Ranks (FG-theta_i-R).")
      theta_i_all = base::sort(base::unlist(interest_cell_type_iGRN[[cell_type]]))
      theta_i_all = theta_i_all[theta_i_all > 0]
      theta_i_rank_all = base::unique(base::data.frame(
        value = theta_i_all,
        norm_rank = base::rank(theta_i_all, ties.method = "min") / base::length(theta_i_all)
      ))
      TF_peak_TG[, "FG-theta_i-R"] = 0
      for (theta_i in TF_peak_TG[, "FG-theta_i"]) {
        TF_peak_TG[TF_peak_TG$`FG-theta_i` == theta_i, "FG-theta_i-R"] =
          theta_i_rank_all[theta_i_rank_all$value == theta_i, "norm_rank"]
      }

      # TF-TG regulatory features from pGRN (theta_p, fit_loss).
      message("[", cell_type, " ", level, " peGRN] [TF-TG Features] Obtain theta_p (FG-theta_p).")
      message("[", cell_type, " ", level, " peGRN] [TF-TG Features] Obtain fit Loss (FG-L).")
      TF_peak_TG[, "FG-theta_p"] = 0
      TF_peak_TG[, "FG-L"] = 0
      for (tftg in tftg_all) {
        TF_peak_TG[TF_peak_TG$TF_TG == tftg, "FG-theta_p"] =
          interest_cell_type_pGRN[[cell_type]][[base::paste0("pred_summary_", level)]][tftg, "theta_p"]
        TF_peak_TG[TF_peak_TG$TF_TG == tftg, "FG-L"] =
          interest_cell_type_pGRN[[cell_type]][[base::paste0("pred_summary_", level)]][tftg, "fit_loss"]
      }
      # Normalized fit loss (scaled to 0-1 range, lower is better)
      message("[", cell_type, " ", level, " peGRN] [TF-TG Features] Obtain Normalized fit Loss (FG-NL).")
      max_loss = base::max(TF_peak_TG[, "FG-L"])
      min_loss = base::min(TF_peak_TG[, "FG-L"])
      TF_peak_TG[, "FG-NL"] = (TF_peak_TG[, "FG-L"] - min_loss) / (max_loss - min_loss)

      # Attach all peak-TG level features to the TF-peak-TG table.
      message("[", cell_type, " ", level, " peGRN] [Peak-TG Features] Add all peak-TG features.")
      PG_feature_names = c(
        "P-NZAR", "P-NZAR-R", "P-NZAA", "P-NZAA-R", "P-AA", "P-AA-R",
        "G-NZER", "G-NZER-R", "G-NZEA", "G-NZEA-R", "G-EA", "G-EA-R"
      )
      TF_peak_TG[, PG_feature_names] = 0
      peaktg_all = base::unique(TF_peak_TG$peak_TG)
      for (peaktg in peaktg_all) {
        TF_peak_TG[TF_peak_TG$peak_TG == peaktg, PG_feature_names] =
          peak_TG[peak_TG$peak_TG == peaktg, PG_feature_names]
      }

      message("[", cell_type, " ", level, " peGRN] Sort out the results.")
      TF_peak_TG = TF_peak_TG[, c(
        "TF", "peak", "TG",
        "FP-S", "FP-S-R",
        "FG-theta_p",
        "FG-L", "FG-NL",
        "FG-theta_i", "FG-theta_i-R",
        "F-NZER", "F-NZER-R", "F-NZEA", "F-NZEA-R", "F-EA", "F-EA-R",
        "P-NZAR", "P-NZAR-R", "P-NZAA", "P-NZAA-R", "P-AA", "P-AA-R",
        "G-NZER", "G-NZER-R", "G-NZEA", "G-NZEA-R", "G-EA", "G-EA-R"
      )]
      TF_peak_TG = TF_peak_TG[base::order(TF_peak_TG$TF, TF_peak_TG$peak, TF_peak_TG$TG), ]
      base::rownames(TF_peak_TG) = NULL

      interest_cell_type_peGRN[[cell_type]][[level]] = TF_peak_TG
    }
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Inferring Predictive Enhancer-driven Gene Regulatory Networks (peGRNs). Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_peGRN)
}

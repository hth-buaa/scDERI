#' @title Extract Pseudotemporal Gene Expression and Activity Profiles along Pseudotime
#'
#' @description
#' For each cell type and pseudotime branch, this function calculates gene importance scores,
#' filters genes and cells, then extracts pseudotime-series matrices of transcription factor (TF) expression,
#' target gene (TG) expression and chromatin activity, ready for downstream GRN inference.
#'
#' @param interest_cell_type_Branches List. Output from \code{\link{order_pseudotime_and_divide_branches}}.
#' @param interest_cell_type_data List. Output from \code{\link{get_interest_cell_type_data}}.
#' @param alpha_TF Numeric (0–1), default 0.63. Quantile threshold for TF filtering.
#'   TFs with importance score above this quantile are retained.
#' @param alpha_TG Numeric (0–1), default 0.69. Quantile threshold for TG filtering.
#'   TGs with importance score above this quantile are retained.
#' @param alpha_cell Numeric (0–1), default 1 (no filtering). Zero-rate threshold for cell filtering.
#'   Cells with zero-expression rate **in both RNA and ATAC** below or equal to this value are retained.
#' @param w_TFNER Numeric (-0.2 to 1). Weight for TF Non-Zero Expression Ratio in TFIS (TF importance score) calculation.
#' @param w_TFNERR Numeric (-0.2 to 1). Weight for TF Non-zero Expression Ratio normalized Rank (max corresponds to 1, min to 0) in TFIS calculation.
#' @param w_TFNEA Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Expression Average in TFIS calculation.
#' @param w_TFNEAR Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Expression Average Rank (max corresponds to 1, min to 0) in TFIS calculation.
#' @param w_TFNAR Numeric (-0.2 to 1). Weight for TF Non-Zero Activity Ratio in TFIS calculation.
#' @param w_TFNARR Numeric (-0.2 to 1). Weight for TF Non-zero Activity Ratio normalized Rank (max corresponds to 1, min to 0) in TFIS calculation.
#' @param w_TFNAA Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Activity Average in TFIS calculation.
#' @param w_TFNAAR Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Activity Average Rank (max corresponds to 1, min to 0) in TFIS calculation.
#' @param w_TGNER Numeric (-0.2 to 1). Weight for TG Non-Zero Expression Ratio in TGIS (TG importance score) calculation.
#' @param w_TGNERR Numeric (-0.2 to 1). Weight for TG Non-zero Expression Ratio normalized Rank (max corresponds to 1, min to 0) in TGIS calculation.
#' @param w_TGNEA Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Expression Average in TGIS calculation.
#' @param w_TGNEAR Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Expression Average Rank (max corresponds to 1, min to 0) in TGIS calculation.
#' @param w_TGNAR Numeric (-0.2 to 1). Weight for TG Non-Zero Activity Ratio in TGIS calculation.
#' @param w_TGNARR Numeric (-0.2 to 1). Weight for TG Non-zero Activity Ratio normalized Rank (max corresponds to 1, min to 0) in TGIS calculation.
#' @param w_TGNAA Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Activity Average in TGIS calculation.
#' @param w_TGNAAR Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Activity Average Rank (max corresponds to 1, min to 0) in TGIS calculation.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#' @param sparsity_cut Nested list or \code{NULL} (default). Optional post-filtering by non-zero rate.
#'   Top-level names must match cell types; each element is a data frame with rows \code{TFE, TGA, TGE}
#'   and columns \code{branch1, branch2, ...}. Values must be between 0 and 1.
#'   Invalid format is silently reset to \code{NULL} with a message.
#'   Re-filters genes based on sparsity thresholds after the initial filtering step using a series of alpha cutoffs.
#'
#' @returns
#' A nested list named by cell type. Each entry contains:
#' \itemize{
#'   \item \code{Branches_TFE}: List of data frames. Rows = TFs, columns = cells ordered by pseudotime, values = normalized expression.
#'   \item \code{Branches_TGA}: List of data frames. Rows = TGs, columns = cells ordered by pseudotime, values = normalized chromatin activity.
#'   \item \code{Branches_TGE}: List of data frames. Rows = TGs, columns = cells ordered by pseudotime, values = normalized expression.
#'   \item \code{Branches_n_cell}: List of integers. Number of retained cells per branch.
#'   \item \code{Branches}: List of data frames. Original branch annotation added with cell-pseudotime index mapping.
#'   \item \code{Branches_TFs}: List of character vectors. Retained TF names per branch.
#'   \item \code{Branches_TGs}: List of character vectors. Retained TG names per branch.
#' }
#'
#' @details
#' For each cell type and branch, the workflow proceeds as follows:
#' \enumerate{
#'   \item Extract scRNA-seq and scATAC-seq count matrices for cells in the current branch.
#'   \item Compute 8 features for TFs and TGs separately: non-zero ratio and max-normalized mean of non-zero values,
#'         for both expression and activity, plus rank-normalized versions of all four metrics.
#'   \item Calculate weighted importance scores TFIS and TGIS using the input weight parameters.
#'   \item Filter genes by score quantiles (\code{alpha_TF}, \code{alpha_TG});
#'         filter cells by zero-expression rate in both omics layers (\code{alpha_cell}).
#'   \item Map cells to integer pseudotime indices and sort cells in ascending pseudotime order.
#'   \item Normalize RNA and ATAC matrices with \code{\link[Seurat]{NormalizeData}}.
#'   \item Extract and organize TFE, TGA, TGE temporal matrices.
#'   \item Optionally apply \code{sparsity_cut} for per-matrix non-zero rate post-filtering,
#'         and sync gene lists across matrices.
#'   \item Return the structured nested list.
#' }
#'
#' @note
#' \itemize{
#'   \item Alpha thresholds and feature weights are core hyperparameters that directly affect the number of retained genes and cells.
#'   \item Parallel computing via \code{\link[parallel]{mclapply}} only works on Unix-like systems; ncores is forced to 1 on Windows.
#'   \item The function creates a \code{2.2 Data Processing - Cell Grouping} directory and temporarily switches the working directory,
#'         which is restored on exit.
#'   \item Check gene non-zero rates after running: excessively sparse matrices may degrade downstream cell grouping performance.
#'         To reduce sparsity, increase weights for non-zero ratio features, or raise alpha thresholds to keep fewer high-quality genes.
#'   \item When \code{sparsity_cut} is applied, TGA and TGE matrices are intersected to keep consistent TG sets.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_Branches = base::readRDS(base::paste0(
#'   "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'   "/interest_cell_type_Branches.rds"
#' ))
#' interest_cell_type_data = base::readRDS(base::paste0(
#'   "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'   "/interest_cell_type_data.rds"
#' ))
#' alpha_TF = 0.7; alpha_TG = 0.7; alpha_cell = 0.9
#' w_TFNER = 1; w_TFNERR = 0.6; w_TFNEA = 0.8; w_TFNEAR = 0.4
#' w_TFNAR = 0.9; w_TFNARR = 0.5; w_TFNAA = 0.7; w_TFNAAR = 0.3
#' w_TGNER = 0.3; w_TGNERR = 0.1; w_TGNEA = 0.3; w_TGNEAR = 0.1
#' w_TGNAR = 1; w_TGNARR = 1; w_TGNAA = 0.3; w_TGNAAR = 0.1
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#'
#' # Example 1: without sparsity_cut
#' interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
#'   interest_cell_type_Branches = interest_cell_type_Branches,
#'   interest_cell_type_data = interest_cell_type_data,
#'   alpha_TF = alpha_TF, alpha_TG = alpha_TG, alpha_cell = alpha_cell,
#'   w_TFNER = w_TFNER, w_TFNERR = w_TFNERR, w_TFNEA = w_TFNEA, w_TFNEAR = w_TFNEAR,
#'   w_TFNAR = w_TFNAR, w_TFNARR = w_TFNARR, w_TFNAA = w_TFNAA, w_TFNAAR = w_TFNAAR,
#'   w_TGNER = w_TGNER, w_TGNERR = w_TGNERR, w_TGNEA = w_TGNEA, w_TGNEAR = w_TGNEAR,
#'   w_TGNAR = w_TGNAR, w_TGNARR = w_TGNARR, w_TGNAA = w_TGNAA, w_TGNAAR = w_TGNAAR,
#'   ncores = ncores
#' )
#'
#' # Example 2: with sparsity_cut (one cell type, one branch)
#' sparsity_cut = list(
#'   cell_name1 = data.frame(
#'     branch1 = c(TFE = 0.00, TGA = 0.00, TGE = 0.03),
#'     branch2 = c(TFE = 0.00, TGA = 0.00, TGE = 0.03),
#'     branch3 = c(TFE = 0.00, TGA = 0.00, TGE = 0.03),
#'     row.names = c("TFE", "TGA", "TGE")
#'   ),
#'   cell_name2 = data.frame(
#'     branch1 = c(TFE = 0.00, TGA = 0.00, TGE = 0.03),
#'     branch2 = c(TFE = 0.00, TGA = 0.00, TGE = 0.03),
#'     branch3 = c(TFE = 0.00, TGA = 0.00, TGE = 0.03),
#'     row.names = c("TFE", "TGA", "TGE")
#'   ),
#' )
#' interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
#'   interest_cell_type_Branches = interest_cell_type_Branches,
#'   interest_cell_type_data = interest_cell_type_data,
#'   alpha_TF = alpha_TF, alpha_TG = alpha_TG, alpha_cell = alpha_cell,
#'   w_TFNER = w_TFNER, w_TFNERR = w_TFNERR, w_TFNEA = w_TFNEA, w_TFNEAR = w_TFNEAR,
#'   w_TFNAR = w_TFNAR, w_TFNARR = w_TFNARR, w_TFNAA = w_TFNAA, w_TFNAAR = w_TFNAAR,
#'   w_TGNER = w_TGNER, w_TGNERR = w_TGNERR, w_TGNEA = w_TGNEA, w_TGNEAR = w_TGNEAR,
#'   w_TGNAR = w_TGNAR, w_TGNARR = w_TGNARR, w_TGNAA = w_TGNAA, w_TGNAAR = w_TGNAAR,
#'   ncores = ncores,
#'   sparsity_cut = sparsity_cut
#' )
#' # or rely solely on sparsity_cut:
#' interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
#'   interest_cell_type_Branches = interest_cell_type_Branches,
#'   interest_cell_type_data = interest_cell_type_data,
#'   alpha_TF = , alpha_TG = 0, alpha_cell = 1,
#'   ncores = ncores,
#'   sparsity_cut = sparsity_cut
#' )
#'
#' # Check minimum nonzero frequency
#' min_nonzero_freq = 1.0
#' for (i in seq_along(interest_cell_type_genes_pseudotime_info)) {
#'   for (j in 1:3) {
#'     for (k in seq_along(interest_cell_type_genes_pseudotime_info[[i]][[j]])) {
#'       df = interest_cell_type_genes_pseudotime_info[[i]][[j]][[k]]
#'       nonzero_freq = rowSums(df != 0) / ncol(df)
#'       min_nonzero_freq = min(min_nonzero_freq, nonzero_freq)
#'     }
#'   }
#' }
#' cat("min_nonzero_freq: ", min_nonzero_freq, "\n"); rm(min_nonzero_freq, df, nonzero_freq)
#' base::saveRDS(
#'   interest_cell_type_genes_pseudotime_info,
#'   file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#' }
#'
get_genes_pseudotime_info = function(
    interest_cell_type_Branches,
    interest_cell_type_data,
    alpha_TF = 0.63, alpha_TG = 0.69, alpha_cell = 0.99,
    w_TFNER = 0.69, w_TFNERR = 1.00, w_TFNEA = 0.67, w_TFNEAR = 0.00,
    w_TFNAR = 0.11, w_TFNARR = 0.07, w_TFNAA = 0.56, w_TFNAAR = 0.12,
    w_TGNER = 0.30, w_TGNERR = 0.88, w_TGNEA = 0.23, w_TGNEAR = 0.12,
    w_TGNAR = 0.48, w_TGNARR = 0.83, w_TGNAA = 0.36, w_TGNAAR = 0.07,
    ncores = 1,
    sparsity_cut = NULL
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Obtaining Gene Expression/Activity Vectors along Pseudotime. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.2 Data Processing - Cell Grouping"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("Working directory set to: ", base::getwd(), ".")

  ### Validate parameters.
  message("Validating sparsity_cut.")
  if (!base::is.null(sparsity_cut)) {
    valid = TRUE

    # Must be a list
    if (!base::is.list(sparsity_cut)) {
      valid = FALSE
    }
    # Top-level names must match cell types
    else if (!base::setequal(base::names(sparsity_cut), base::names(interest_cell_type_Branches))) {
      valid = FALSE
    }
    else {
      # Validate each cell_type entry
      for (cell_type in base::names(interest_cell_type_Branches)) {
        x = sparsity_cut[[cell_type]]
        expected_cols = base::paste0("branch", base::seq_along(interest_cell_type_Branches[[cell_type]]))

        # Must be a data.frame
        if (!base::is.data.frame(x)) { valid = FALSE; break }
        # Rownames must be TFE, TGA, TGE
        if (!base::setequal(base::rownames(x), c("TFE", "TGA", "TGE"))) { valid = FALSE; break }
        # Colnames must match branch count
        if (!base::setequal(base::colnames(x), expected_cols)) { valid = FALSE; break }
        # All values must be numeric and non-negative
        x_mat = base::as.matrix(x)
        if (!base::is.numeric(x_mat)) { valid = FALSE; break }
        if (base::anyNA(x_mat) || base::any(x_mat < 0) || base::any(x_mat > 1)) { valid = FALSE; break }
      }
    }

    if (!valid) {
      message("Invalid sparsity_cut format; resetting to NULL.")
      sparsity_cut = NULL
    }
  }

  message("Validating alpha params.")
  if (alpha_TF < 0 || alpha_TF > 1 || alpha_TG < 0 || alpha_TG > 1 || alpha_cell < 0 || alpha_cell > 1) {
    base::stop("alpha_TF, alpha_TG, alpha_cell must be between 0 and 1.")
  }

  message("Validating input object structure.")
  required_elements = c("interest_scRNAseq_for_GRN", "interest_scATACseq_for_GRN", "interest_gene_for_GRN", "interest_TGs", "interest_TFs")
  if (!all(required_elements %in% base::names(interest_cell_type_data[[1]]))) {
    base::stop("interest_cell_type_data format is incorrect. Please use the output of get_interest_cell_type_data.")
  }

  if (.Platform$OS.type == "windows") {
    ncores = 1
    message("Parallel computing is not supported on Windows. ncores is forced to 1.")
  }

  ### Define per-cell-type processing function.
  get_genes_pseudotime_info_process = function(
    Branches,
    interest_scRNAseq_for_GRN, interest_scATACseq_for_GRN,
    interest_gene_for_GRN, interest_TGs, interest_TFs,
    alpha_TF, w_TFNER, w_TFNERR, w_TFNEA, w_TFNEAR, w_TFNAR, w_TFNARR, w_TFNAA, w_TFNAAR,
    alpha_TG, w_TGNER, w_TGNERR, w_TGNEA, w_TGNEAR, w_TGNAR, w_TGNARR, w_TGNAA, w_TGNAAR,
    alpha_cell)
  {
    # Initialize per-branch storage lists
    message("Initializing per-branch storage lists.")
    Branches_RNA = base::list()
    Branches_ATAC = base::list()
    Branches_gene_for_GRN = base::list()
    Branches_TGs = base::list()
    Branches_TFs = base::list()
    Branches_n_cell = base::list()

    # Feature calculation helper: computes NER/NEA/NAR/NAA + ranks + weighted IS score
    calculate_features = function(exp_mat, # gene expression matrix (genes x cells)
                                  act_mat, # gene activity matrix (genes x cells)
                                  genes,   # gene name vector (TF or TG)
                                  weights, # 8 feature weights (NER, NERR, NEA, NEAR, NAR, NARR, NAA, NAAR)
                                  prefix)  # column name prefix ("TF" or "TG")
    {
      # Non-zero rates (NER / NAR)
      ner = base::rowMeans(exp_mat != 0)
      nar = base::rowMeans(act_mat != 0)

      # Mean of non-zero values (NEA / NAA)
      exp_na = exp_mat; exp_na[exp_na == 0] = NA
      nea = base::rowMeans(exp_na, na.rm = TRUE); rm(exp_na)

      act_na = act_mat; act_na[act_na == 0] = NA
      naa = base::rowMeans(act_na, na.rm = TRUE); rm(act_na)

      # Fill NA with 0
      ner[base::is.na(ner)] = 0
      nea[base::is.na(nea)] = 0
      nar[base::is.na(nar)] = 0
      naa[base::is.na(naa)] = 0

      # Normalize by max
      max_ner = base::max(ner)
      if (max_ner > 0) {
        ner = ner / max_ner
      } else {
        stop("The maximum of the NER feature is zero, causing normalization to fail.")
      }

      max_nea = base::max(nea)
      if (max_nea > 0) {
        nea = nea / max_nea
      } else {
        stop("The maximum of the NEA feature is zero, causing normalization to fail.")
      }

      max_nar = base::max(nar)
      if (max_nar > 0) {
        nar = nar / max_nar
      } else {
        stop("The maximum of the NAR feature is zero, causing normalization to fail.")
      }

      max_naa = base::max(naa)
      if (max_naa > 0) {
        naa = naa / max_naa
      } else {
        stop("The maximum of the NAA feature is zero, causing normalization to fail.")
      }

      # Rank-normalize
      n = length(genes)
      nerr = base::rank(ner) / n
      near = base::rank(nea) / n
      narr = base::rank(nar) / n
      naar = base::rank(naa) / n

      # Weighted importance score
      feat_mat = base::cbind(ner, nerr, nea, near, nar, narr, naa, naar)
      is_score = base::as.vector(feat_mat %*% weights)

      # Output feature data frame
      df = data.frame(
        feat_mat,
        IS = is_score,
        row.names = genes,
        check.names = FALSE
      )
      base::colnames(df) = paste0(prefix, c("NER", "NERR", "NEA", "NEAR", "NAR", "NARR", "NAA", "NAAR", "IS"))

      return(df)
    }

    for (i in base::seq_along(Branches)) {
      message("Extracting RNA count matrix for branch ", i, ".")
      Branches_RNA[[i]] = interest_scRNAseq_for_GRN[, base::rownames(Branches[[i]])]
      rna_assay = SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "counts")

      message("Extracting ATAC count matrix for branch ", i, ".")
      Branches_ATAC[[i]] = interest_scATACseq_for_GRN[base::rownames(Branches_RNA[[i]]), base::colnames(Branches_RNA[[i]])]
      atac_assay = SeuratObject::GetAssayData(Branches_ATAC[[i]], layer = "counts")

      message("Computing TF feature matrix for branch ", i, ".")
      TF_features = calculate_features(
        exp_mat = base::as.matrix(rna_assay[interest_TFs, ]),
        act_mat = base::as.matrix(atac_assay[interest_TFs, ]),
        genes = interest_TFs,
        weights = c(w_TFNER, w_TFNERR, w_TFNEA, w_TFNEAR, w_TFNAR, w_TFNARR, w_TFNAA, w_TFNAAR),
        prefix = "TF"
      )

      message("Computing TG feature matrix for branch ", i, ".")
      TG_features = calculate_features(
        exp_mat = base::as.matrix(rna_assay[interest_TGs, ]),
        act_mat = base::as.matrix(atac_assay[interest_TGs, ]),
        genes = interest_TGs,
        weights = c(w_TGNER, w_TGNERR, w_TGNEA, w_TGNEAR, w_TGNAR, w_TGNARR, w_TGNAA, w_TGNAAR),
        prefix = "TG"
      )
      rm(rna_assay, atac_assay)

      message("Filtering genes for branch ", i, ".")
      TFs_filter = base::rownames(TF_features)[TF_features$TFIS > stats::quantile(TF_features$TFIS, alpha_TF)]
      TGs_filter = base::rownames(TG_features)[TG_features$TGIS > stats::quantile(TG_features$TGIS, alpha_TG)]
      genes_filter = base::sort(base::unique(c(TFs_filter, TGs_filter)))
      Branches_RNA[[i]] = subset(Branches_RNA[[i]], features = genes_filter)
      Branches_ATAC[[i]] = subset(Branches_ATAC[[i]], features = genes_filter)
      rm(genes_filter)

      message("Filtering cells for branch ", i, ".")
      zero_freq_RNA = base::colMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "counts")) == 0))
      zero_freq_ATAC = base::colMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_ATAC[[i]], layer = "counts")) == 0))
      cells_filter = base::intersect(
        base::names(zero_freq_RNA[zero_freq_RNA <= alpha_cell]),
        base::names(zero_freq_ATAC[zero_freq_ATAC <= alpha_cell])
      )
      rm(zero_freq_RNA, zero_freq_ATAC)
      Branches_RNA[[i]] = subset(Branches_RNA[[i]], cells = cells_filter)
      Branches_ATAC[[i]] = subset(Branches_ATAC[[i]], cells = cells_filter)
      rm(cells_filter)

      message("Collecting TGs and TFs for branch ", i, ".")
      Branches_gene_for_GRN[[i]] = base::intersect(interest_gene_for_GRN, base::rownames(Branches_RNA[[i]]))
      Branches_TGs[[i]] = base::intersect(base::intersect(interest_TGs, Branches_gene_for_GRN[[i]]), TGs_filter); rm(TGs_filter)
      Branches_TFs[[i]] = base::intersect(base::intersect(interest_TFs, Branches_gene_for_GRN[[i]]), TFs_filter); rm(TFs_filter)

      message("Counting retained cells for branch ", i, ".")
      Branches_gene_for_GRN[[i]] = base::union(Branches_TGs[[i]], Branches_TFs[[i]])
      Branches_RNA[[i]] = Branches_RNA[[i]][base::intersect(base::rownames(Branches_RNA[[i]]), interest_gene_for_GRN), ]
      Branches_ATAC[[i]] = Branches_ATAC[[i]][base::rownames(Branches_RNA[[i]]), ]
      Branches_n_cell[[i]] = base::ncol(Branches_RNA[[i]])
    }
    rm(calculate_features)

    # Map cells to pseudotime indices k~
    message("Mapping cells to pseudotime indices per branch.")
    for (i in base::seq_along(Branches)) {
      message("Building cell-to-pseudotime map for branch ", i, ".")
      Branches[[i]]["cell"] = base::rownames(Branches[[i]])
      Branches[[i]]["k~"] = 1:base::nrow(Branches[[i]])
      name_mapping = stats::setNames(Branches[[i]]$`k~`, Branches[[i]]$cell)

      message("Renaming cells by pseudotime index for branch ", i, ".")
      base::colnames(Branches_RNA[[i]]) = name_mapping[base::colnames(Branches_RNA[[i]])]
      base::colnames(Branches_ATAC[[i]]) = name_mapping[base::colnames(Branches_ATAC[[i]])]
    }

    # Normalize matrices
    message("Normalizing expression and activity matrices per branch.")
    for (i in base::seq_along(Branches)) {
      message("Normalizing scRNA-seq for branch ", i, ".")
      Branches_RNA[[i]] = Seurat::NormalizeData(Branches_RNA[[i]])

      message("Normalizing scATAC-seq for branch ", i, ".")
      Branches_ATAC[[i]] = Seurat::NormalizeData(Branches_ATAC[[i]])
    }

    # Extract temporal vectors for TFE / TGA / TGE
    message("Initializing TFE / TGA / TGE temporal vector lists.")
    Branches_TFE = base::list()
    Branches_TGA = base::list()
    Branches_TGE = base::list()
    for (i in base::seq_along(Branches)) {
      message("Extracting TF temporal expression for branch ", i, ".")
      Branches_TFE[[i]] =
        base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "data")))[Branches_TFs[[i]], ]
      Branches_TFE[[i]] = Branches_TFE[[i]][, base::order(base::as.integer(base::colnames(Branches_TFE[[i]])))]

      message("Extracting TG temporal activity for branch ", i, ".")
      Branches_TGA[[i]] =
        base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_ATAC[[i]], layer = "data")))[Branches_TGs[[i]], ]
      Branches_TGA[[i]] = Branches_TGA[[i]][, base::order(base::as.integer(base::colnames(Branches_TGA[[i]])))]

      message("Extracting TG temporal expression for branch ", i, ".")
      Branches_TGE[[i]] =
        base::as.data.frame(SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "data"))[Branches_TGs[[i]], ]
      Branches_TGE[[i]] = Branches_TGE[[i]][, base::order(base::as.integer(base::colnames(Branches_TGE[[i]])))]
    }

    return(base::list("Branches_TFE" = Branches_TFE,
                      "Branches_TGA" = Branches_TGA,
                      "Branches_TGE" = Branches_TGE,
                      "Branches_n_cell" = Branches_n_cell,
                      "Branches" = Branches,
                      "Branches_TFs" = Branches_TFs,
                      "Branches_TGs" = Branches_TGs))
  }

  ### Run per-cell-type processing in parallel.
  interest_cell_type_genes_pseudotime_info = parallel::mclapply(
    X = base::names(interest_cell_type_Branches),
    FUN = function(cell_type) {
      message("Extracting TFE/TGA/TGE pseudotime vectors for cell type ", cell_type, ".")
      get_genes_pseudotime_info_process(
        Branches = interest_cell_type_Branches[[cell_type]],
        interest_scRNAseq_for_GRN = interest_cell_type_data[[cell_type]][["interest_scRNAseq_for_GRN"]],
        interest_scATACseq_for_GRN = interest_cell_type_data[[cell_type]][["interest_scATACseq_for_GRN"]],
        interest_gene_for_GRN = interest_cell_type_data[[cell_type]][["interest_gene_for_GRN"]],
        interest_TGs = interest_cell_type_data[[cell_type]][["interest_TGs"]],
        interest_TFs = interest_cell_type_data[[cell_type]][["interest_TFs"]],
        alpha_TF = alpha_TF,
        w_TFNER = w_TFNER, w_TFNERR = w_TFNERR, w_TFNEA = w_TFNEA, w_TFNEAR = w_TFNEAR,
        w_TFNAR = w_TFNAR, w_TFNARR = w_TFNARR, w_TFNAA = w_TFNAA, w_TFNAAR = w_TFNAAR,
        alpha_TG = alpha_TG,
        w_TGNER = w_TGNER, w_TGNERR = w_TGNERR, w_TGNEA = w_TGNEA, w_TGNEAR = w_TGNEAR,
        w_TGNAR = w_TGNAR, w_TGNARR = w_TGNARR, w_TGNAA = w_TGNAA, w_TGNAAR = w_TGNAAR,
        alpha_cell = alpha_cell
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_genes_pseudotime_info) = base::names(interest_cell_type_Branches)

  ### Apply sparsity-based gene filtering
  if (!base::is.null(sparsity_cut)) {
    message("Applying sparsity cutoff filtering.")
    for (cell_type in base::names(interest_cell_type_Branches)) {
      # Filter TFE / TGA / TGE matrices by non-zero rate threshold
      for (feat in base::rownames(sparsity_cut[[cell_type]])) {
        interest_cell_type_genes_pseudotime_info[[cell_type]][[paste0("Branches_", feat)]] = lapply(
          base::seq_along(interest_cell_type_genes_pseudotime_info[[cell_type]][[paste0("Branches_", feat)]]),
          function(n) {
            df = interest_cell_type_genes_pseudotime_info[[cell_type]][[paste0("Branches_", feat)]][[n]]
            df[base::rowMeans(df != 0) >= sparsity_cut[[cell_type]][feat, paste0("branch", n)], ]
          }
        )
      }
      # Sync TFs / TGs lists with filtered matrices
      for (n in base::seq_along(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]])) {
        interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFs"]][[n]] =
          base::rownames(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]][[n]])
        TGs = base::intersect(
          base::rownames(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGA"]][[n]]),
          base::rownames(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGE"]][[n]])
        )
        interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGA"]][[n]] =
          interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGA"]][[n]][TGs, ]
        interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGE"]][[n]] =
          interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGE"]][[n]][TGs, ]
        interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGs"]][[n]] = TGs
      }; rm(TGs)
    }
  }

  ### End.
  base::setwd(original_dir)
  message("Working directory restored to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Obtaining Gene Expression/Activity Vectors along Pseudotime. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_genes_pseudotime_info)
}

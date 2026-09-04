#' @title Construct Cell-Type-Specific Initial GRNs (iGRNs) via JASPAR PWM Matching
#'
#' @description For each cell type,
#' this function scans sequences of promoter-associated accessible peaks with transcription factor binding site position weight matrices from the JASPAR2024 database,
#' and constructs an initial gene regulatory network matrix.
#' Rows represent target genes, columns represent transcription factors, and values are normalized PWM matching scores.
#'
#' @param interest_cell_type_genes_pseudotime_info List. Output from \code{\link{get_genes_pseudotime_info}}, providing the target TF and TG sets.
#' @param results_identify_TGs List. Output from \code{\link{identify_TGs}} or \code{\link{filter_peak_annotations}}.
#'   Must contain a \code{peak_anno_promoter} element with annotation of promoter-associated ATAC peaks.
#' @param genome \code{BSgenome} object. Reference genome for sequence extraction.
#'   Defaults to \code{BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38}.
#'   \code{BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19} and other species-specific BSgenome objects are also supported.
#' @param min_score_for_matchPWM Character. Minimum score threshold for PWM matching, passed to \code{\link[Biostrings]{matchPWM}}.
#'   Accepts percentage strings such as \code{"80\%"}. Higher thresholds produce sparser networks with higher specificity.
#'   Defaults to \code{"80\%"}.
#' @param species Character. Species name used to query the JASPAR2024 database. Defaults to \code{"Homo sapiens"}. See in \code{\link{get_TFs_from_priors}}.
#' @param collection Character vector. JASPAR matrix collections to use. Defaults to CORE and multiple additional commonly used collections. See in \code{\link{get_TFs_from_priors}}.
#' @param output_predicted_TFBS Logical. Whether to include predicted TFBS sequences in the output CSV.
#'   Defaults to \code{FALSE}, which exports only gene, TF, max score and peak name, and reduces runtime and output file size.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#'
#' @return
#' A named list with one element per cell type. Each element is a data frame used as an association matrix:
#' \itemize{
#'   \item Rows: Target gene symbols
#'   \item Columns: Transcription factor symbols
#'   \item Values: Normalized PWM matching scores, ranging from 0 to approximately 0.95
#' }
#' A value of 0 indicates no significant TF-TG regulatory pair was predicted.
#'
#' @details
#' The function constructs initial GRNs independently for each cell type with the following workflow:
#' \enumerate{
#'   \item \strong{Gene set restriction}: Extracts the union of TFs and TGs across all branches from the pseudotime result to ensure consistency with upstream analysis.
#'   \item \strong{Promoter peak mapping}: Retrieves genomic coordinates of promoter-associated accessible peaks from the input annotation, and retains only peaks linked to target TGs.
#'   \item \strong{PWM retrieval and conversion}: Fetches TF PWMs from JASPAR2024 for the specified species and collections, and filters to retain only target TFs.
#'     Log-odds PWMs are converted to base frequency matrices using the formula \eqn{2^{PWM} / 4}.
#'   \item \strong{Parallel sequence matching}: Extracts nucleotide sequences of all peaks from the reference genome, and scans each PWM against all peak sequences in parallel.
#'     The maximum matching score per peak is recorded.
#'   \item \strong{Score aggregation}: Aggregates peak-level scores to TF-TG pairs by calculating the mean score across all peaks associated with the same gene.
#'   \item \strong{Score normalization}: Divides all scores by 1.05 times the global maximum score, capping values at approximately 0.95 to avoid numerical saturation in downstream modeling.
#'   \item \strong{Result export}: Saves per-cell-type matching details (CSV) and intermediate results (RDS) in the \code{3 get iGRN} directory.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates a \code{3 get iGRN} directory and temporarily switches the working directory into it; the original directory is restored on exit.
#'   \item \code{min_score_for_matchPWM} is the core parameter controlling network sparsity. Adjust it according to the desired balance between specificity and recall.
#'   \item Only TFs present in both the input gene set and the JASPAR database are retained in the final matrix.
#'   \item For TF-TG pairs linked to multiple promoter peaks, the final score is the mean of the maximum matching scores of all corresponding peaks.
#'   \item Parallel computation relies on \code{\link[parallel]{mclapply}}, which is only supported on Linux and macOS.
#'   \item Required dependencies include \code{JASPAR2024}, \code{TFBSTools}, \code{Biostrings}, a matching \code{BSgenome} package,
#'     as well as \code{data.table}, \code{dplyr}, \code{tidyr} and \code{tibble}.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_genes_pseudotime_info = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#' results_identify_TGs = base::readRDS(
#'   "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds"
#' )
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' min_score_for_matchPWM = "80%"
#' ncores = ceiling(parallel::detectCores() / 2) # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_iGRN = get_iGRN_by_TFBS_pwm_by_JASPAR2024(
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   results_identify_TGs = results_identify_TGs,
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   min_score_for_matchPWM = min_score_for_matchPWM,
#'   species = "Homo sapiens",
#'   collection = c(
#'     "CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM",
#'     "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"
#'   ),
#'   output_predicted_TFBS = FALSE,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_iGRN, file = "./3 get iGRN/interest_cell_type_iGRN.rds")
#' }
get_iGRN_by_TFBS_pwm_by_JASPAR2024 = function(
    interest_cell_type_genes_pseudotime_info,
    results_identify_TGs,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    min_score_for_matchPWM = "80%",
    species = "Homo sapiens",
    collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII",
                   "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"),
    output_predicted_TFBS = FALSE,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Constructing Initial Gene Regulatory Networks (GRNs) from Transcription Factor Binding Site (TFBS) Position Weight Matrices (PWMs). Time: ", t_start, ".")
  message("Only considering TFs and TGs that appear in pseudotime analysis.")
  message("Generating a TG-TF association matrix with known TFs filtered by TF-TFBS match strength (rows: TGs, columns: TFs, regulatory strength based on TF-TFBS match scores).")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Defining inner function.
  get_iGRN_by_TFBS_pwm_by_JASPAR2024_process = function(
    cell_type = cell_type,
    TGs,
    TFs,
    peak_anno_promoter,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    min_score_for_matchPWM = "80%",
    ncores = 1,
    species = "Homo sapiens",
    collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII",
                   "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"),
    output_predicted_TFBS = FALSE)
  {
    # Constructing TG-Promoter mapping.
    message("Constructing TG-Promoter mapping.")
    map_from_TGs_to_promoter = function(peak_anno_promoter) {
      TGs_and_promoter = base::data.frame(
        chromosome_name = peak_anno_promoter$seqnames,
        start_position = peak_anno_promoter$start,
        end_position = peak_anno_promoter$end,
        TGs = peak_anno_promoter$SYMBOL,
        peak = peak_anno_promoter$peak,
        stringsAsFactors = FALSE
      )
      TGs_and_promoter$chromosome_name = base::gsub("^chr", "", TGs_and_promoter$chromosome_name)
      TGs_and_promoter$chromosome_name = base::as.character(TGs_and_promoter$chromosome_name)
      return(TGs_and_promoter)
    }
    TGs_and_promoter = map_from_TGs_to_promoter(peak_anno_promoter)
    TGs_and_promoter = TGs_and_promoter[TGs_and_promoter$TGs %in% TGs, , drop = FALSE]

    # Defining a function to extract promoter sequences.
    get_peak_sequence = function(chromosome, start, end, genome = genome) {
      Biostrings::getSeq(genome, base::paste0("chr", chromosome), start, end)
    }

    # Retrieving PWMs of TFBS from JASPAR database (typical TF binding length: 5-20 bp)
    message("Retrieving PWMs of TFBS from JASPAR database (typical TF binding length: 5-20 bp).")
    pwm = TFBSTools::getMatrixSet(
      x = JASPAR2024::db(JASPAR2024::JASPAR2024()),
      opts = base::list(all = FALSE, species = species, collection = collection,
                        all_versions = FALSE, matrixtype = "PWM")
    )

    # Retaining only TFs present in pseudotime analysis data.
    message("Retaining only TFs present in pseudotime analysis data.")
    pwm@listData = base::Filter(function(x) x@name %in% TFs, pwm@listData)

    # Convert log-odds PWM to base frequency matrix (formula: 2^PWM / 4, assuming uniform 0.25 background).
    # Rows: nucleotides (A/C/G/T), columns: binding positions, values: occurrence frequency per position.
    message("Normalize each column of the position weight matrix to sum to 1 (rows: nucleotide bases, columns: binding site positions, values: base occurrence probabilities, column sum constraint: 1).")
    for (i in base::seq_along(pwm@listData)) {
      pwm@listData[[i]]@profileMatrix = 2 ^ (pwm@listData[[i]]@profileMatrix) / 4
    }

    # Calculating the sequence of each promoter.
    # Using parallel::mclapply for parallel computing (Only effective on Linux systems).
    message("Calculating the sequence of each promoter.")
    message("The number of CPU cores for parallel::mclapply parallel computing is: ", ncores)
    promoter_seq = base::list()
    promoter_seq = parallel::mclapply(1:base::dim(TGs_and_promoter)[1], function(i) {
      get_peak_sequence(chromosome = TGs_and_promoter[i, "chromosome_name"],
                        start = TGs_and_promoter[i, "start_position"],
                        end = TGs_and_promoter[i, "end_position"],
                        genome = genome)
    }, mc.cores = ncores)

    # Processing each PWM and matching to TFGs based on PWM patterns.
    # Using parallel::mclapply for parallel computing (Only effective on Linux systems).
    message("Processing each PWM and matching against peak sequences.")
    message("The minimum score threshold for PWM matching to TFGs is: ", min_score_for_matchPWM)
    message("The number of CPU cores for parallel::mclapply parallel computing is: ", ncores)
    TF_match_peak_results = base::list()
    length_pwm = base::length(base::names(pwm))
    for (name in base::names(pwm)) {
      pos = base::match(name, base::names(pwm))
      if (pos %% 10 == 1 | pos == length_pwm) {
        message(pos, "th in ", length_pwm)
      }
      results = parallel::mclapply(1:base::nrow(TGs_and_promoter), function(i) {
        predicted_TFBS = Biostrings::matchPWM(pwm = pwm[[name]]@profileMatrix,
                                              subject = promoter_seq[[i]],
                                              min.score = min_score_for_matchPWM,
                                              with.score = TRUE)
        gene_name = TGs_and_promoter$TGs[i]
        peak_name = TGs_and_promoter$peak[i]
        if (
          !base::is.null(predicted_TFBS) &&
          base::length(gene_name) != 0 &&
          base::length(peak_name) != 0 &&
          !base::is.na(base::mean(predicted_TFBS@elementMetadata@listData[["score"]]))
        ) {
          if (output_predicted_TFBS == TRUE) {
            result_line = base::paste(gene_name,
                                      pwm[[name]]@name,
                                      predicted_TFBS,
                                      base::max(predicted_TFBS@elementMetadata@listData[["score"]]),
                                      peak_name,
                                      sep = ",")
          } else {
            result_line = base::paste(gene_name,
                                      pwm[[name]]@name,
                                      "-",
                                      base::max(predicted_TFBS@elementMetadata@listData[["score"]]),
                                      peak_name,
                                      sep = ",")
          }
          return(result_line)
        } else {
          return(NULL)
        }
      }, mc.cores = ncores)
      results = base::unlist(base::Filter(base::Negate(is.null), results))

      TF_match_peak_results[[name]] = results
    }

    message("Saving parallel computation results to './3 get iGRN/", cell_type, "_TF_match_peak_results.rds'")
    base::saveRDS(TF_match_peak_results, file = base::paste0(cell_type, "_TF_match_peak_results.rds"))

    # Writing results to CSV file.
    message("Saving parallel computation results to './3 get iGRN/", cell_type, "_TF_match_peak_results.csv'")
    file_conn = base::file(base::paste0(cell_type, "_TF_match_peak_results.csv"), "w")
    base::writeLines("TGs, Potential_TFs, Predicted_TFBS, match_score, peaks", file_conn)
    for (result_list in TF_match_peak_results) {
      base::writeLines(as.character(result_list), file_conn)
    }
    base::close(file_conn)

    # Generating TF-TG association matrix (rows: TGs, columns: TFs, values: match scores normalized by 1.05 times max match score).
    message("Generating TF-TG association matrix (rows: TGs, columns: TFs, values: match scores normalized by 1.05 times max match score).")
    potential_regulation = data.table::fread(base::paste0(cell_type, "_TF_match_peak_results.csv"))
    potential_regulation$peaks = NULL
    potential_regulation = potential_regulation[potential_regulation$Predicted_TFBS != "" & !base::is.na(potential_regulation$Predicted_TFBS), ]
    potential_regulation$Predicted_TFBS = NULL
    potential_regulation = potential_regulation[potential_regulation$Potential_TFs %in% TFs, ]
    # library(dplyr)
    potential_regulation = potential_regulation %>%
      dplyr::group_by(TGs, Potential_TFs) %>%
      dplyr::summarise(match_score = base::mean(match_score)) %>%
      dplyr::ungroup()
    M_TGs_TFs = tidyr::pivot_wider(potential_regulation, names_from = Potential_TFs, values_from = match_score, values_fill = 0)
    M_TGs_TFs = tibble::column_to_rownames(M_TGs_TFs, "TGs")
    M_TGs_TFs = M_TGs_TFs / (base::max(M_TGs_TFs) * 1.05)

    message("Generating an association matrix containing ", nrow(M_TGs_TFs), " TGs and ", ncol(M_TGs_TFs), " TFs with a sparsity of ", base::round(base::sum(M_TGs_TFs == 0) / (base::nrow(M_TGs_TFs) * base::ncol(M_TGs_TFs)), 3))
    return(M_TGs_TFs)
  }

  ### Running inner function.
  interest_cell_type_iGRN_all_TGTF_pairs = parallel::mclapply(
    X = base::names(interest_cell_type_genes_pseudotime_info),
    FUN = function(cell_type) {
      message("Preparing to construct iGRN for cell type ", cell_type, ".")
      get_iGRN_by_TFBS_pwm_by_JASPAR2024_process(
        cell_type = cell_type,
        TGs = base::sort(base::unique(base::Reduce(
          base::union,
          interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGs"]]
        ))),
        TFs = base::sort(base::unique(base::Reduce(
          base::union,
          interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFs"]]
        ))),
        peak_anno_promoter = results_identify_TGs[["peak_anno_promoter"]],
        genome = genome,
        min_score_for_matchPWM = min_score_for_matchPWM,
        ncores = ncores,
        species = species,
        collection = collection,
        output_predicted_TFBS = output_predicted_TFBS
      )},
    mc.cores = 1)
  base::names(interest_cell_type_iGRN_all_TGTF_pairs) = base::names(interest_cell_type_genes_pseudotime_info)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Constructing Initial Gene Regulatory Networks (GRNs) from Transcription Factor Binding Site (TFBS) Position Weight Matrices (PWMs). Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_iGRN_all_TGTF_pairs)
}

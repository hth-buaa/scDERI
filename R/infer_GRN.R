#' @title Infer Gene Regulatory Networks with Trained Models
#'
#' @description
#' Predicts optimal regulatory strengths and corresponding prediction losses for all candidate TF-TG pairs using pre-trained scDERI models.
#' Supports subset processing by cell type and branch, chunk-based parallel computation, and automatic saving of branch-level results.
#'
#' @param interest_cell_type_branch_data_for_GRN_infer List. Output from \code{\link{get_data_for_inferring_GRN}}, containing branch-level TF-TG pairs with time-series expression features.
#' @param interest_cell_type_tau Named numeric vector. Output from \code{\link{get_cell_type_tau}}, cell-type-specific binarization thresholds.
#' @param interest_cell_type_branch_model_train Nested list. Output from \code{\link{model_train}}, containing trained model parameters for each branch.
#' @param interest_cell_type_group List. Output from \code{\link{cell_grouping}}, containing pseudotime interval lengths for each branch.
#' @param select_cell_type Character vector or \code{NULL}.
#' Cell types to process. \code{NULL} means processing all cell types. Defaults to \code{NULL}.
#' @param select_branch Numeric vector or \code{NULL}. Branch indices to process. \code{NULL} means processing all branches. Defaults to \code{NULL}.
#' @param each Positive integer. Maximum number of regulatory pairs per parallel chunk. Adjust to control peak memory usage. Defaults to 1000000.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#'
#' @return
#' A nested list named \code{interest_cell_type_pGRN}, organized by cell type and branch. Each branch contains two data frames:
#' \itemize{
#'   \item \code{theta_p}: Rows are candidate TF-TG pairs to predict, columns are trained models. Values are optimized regulatory strengths.
#'   \item \code{loss}: Rows and columns correspond to \code{theta_p}. Values are the minimum prediction loss corresponding to the optimal theta.
#' }
#' Each branch's result is also saved as a separate RDS file in the \code{5 Infer GRN} directory.
#'
#' @details
#' The inference workflow for each selected branch is as follows:
#' \enumerate{
#'   \item Splits all candidate TF-TG pairs into multiple chunks according to the \code{each} parameter.
#'   \item Processes chunks in parallel: for each candidate pair, iterates over all trained models and uses Brent's method to find the theta value (within 0.01–0.99) that minimizes prediction loss.
#'   \item Merges results from all chunks to form branch-level theta_p and loss matrices.
#'   \item Automatically saves the branch-level result as an RDS file to the \code{5 Infer GRN} directory.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and temporarily switches to the \code{5 Infer GRN} directory, and restores the original working directory on exit.
#'   \item The \code{each} parameter controls peak memory usage. Reduce its value if you encounter out-of-memory errors.
#'   \item Use \code{select_cell_type} and \code{select_branch} to process data in batches, which is suitable for large datasets and supports interrupted recovery.
#'   \item Parallel computation relies on \code{parallel::mclapply}. Set \code{ncores = 1} on Windows systems.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_data_for_GRN_infer = base::readRDS(
#'   "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds"
#' )
#' interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
#' interest_cell_type_branch_model_train = base::readRDS(
#'   "./4.2 Build Prediction Model/interest_cell_type_branch_model_train.rds"
#' )
#' interest_cell_type_group = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' each = 1000
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#'
#' # way 1: one step (If the running memory is sufficient)
#' interest_cell_type_pGRN = infer_GRN(
#'   interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
#'   interest_cell_type_tau = interest_cell_type_tau,
#'   interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
#'   interest_cell_type_group = interest_cell_type_group,
#'   select_cell_type = NULL,
#'   select_branch = NULL,
#'   each = each,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#'
#' # way 2: more steps but stable
#' # If the process is terminated due to memory issues,
#' # you can manually continue running it branch by branch
#' # (saving the results to interest_cell_type_pGRN).
#' interest_cell_type_pGRN = base::list()
#' for (cell_type in base::names(interest_cell_type_branch_data_for_GRN_infer)) {
#'   message("Infer cell type ", cell_type, ".")
#'   interest_cell_type_pGRN[[cell_type]] = base::list()
#'   for (n in 1:(base::length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)) {
#'     message("Infer cell type ", cell_type, " branch ", n, ".")
#'     interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]] = infer_GRN(
#'       interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
#'       interest_cell_type_tau = interest_cell_type_tau,
#'       interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
#'       interest_cell_type_group = interest_cell_type_group,
#'       select_cell_type = cell_type,
#'       select_branch = n,
#'       each = each,
#'       ncores = ncores
#'     )[[cell_type]][[base::paste0("branch", n)]]
#'   }
#' }
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
infer_GRN = function(
    interest_cell_type_branch_data_for_GRN_infer,
    interest_cell_type_tau,
    interest_cell_type_branch_model_train,
    interest_cell_type_group,
    select_cell_type = NULL,
    select_branch = NULL,
    each = 1000000,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Inferring GRN. Time: ", t_start, ".")
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

  ### Initialize output result list.
  interest_cell_type_pGRN = base::list()

  ### Determine cell types to process (all by default, or user-specified subset).
  if (base::is.null(select_cell_type)) {
    range_cell_type = base::names(interest_cell_type_branch_data_for_GRN_infer)
  } else {
    range_cell_type = select_cell_type
  }

  for (cell_type in range_cell_type) {
    interest_cell_type_pGRN[[cell_type]] = base::list()
    tau = as.numeric(interest_cell_type_tau[cell_type])

    # Determine branches to process (all by default, or user-specified subset).
    if (base::is.null(select_branch)) {
      range_branch = 1:(base::length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)
    } else {
      range_branch = select_branch
    }

    for (n in range_branch) {
      message("Process branch ", n, " in cell type ", cell_type, ".")

      # Get all TF-TG pairs to be predicted.
      re_pred = base::rownames(interest_cell_type_branch_data_for_GRN_infer[[cell_type]][[base::paste0("branch", n)]])
      N_pred_all = base::length(re_pred)
      message(N_pred_all, " potential regulatory relationships will be predicted.")

      # Get names of all trained models (each trained TF-TG pair corresponds to one model).
      re_train = base::sort(base::names(interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]]))
      message(base::length(re_train), " trained regulatory models will be applied for prediction.")

      # Partition prediction pairs into small chunks for parallel processing.
      message("Partitioning data.")
      n_part = base::ceiling(N_pred_all / each)
      each_last = N_pred_all - each * (n_part - 1)

      # Control progress output frequency: print progress every idx chunks.
      idx = max(ceiling(n_part / 50), 2)
      re_pred_list = base::list()
      for (s in 1:n_part) {
        if (s %% idx ==1) {
          message(paste0("Partitioning data progress: ", s, " / ", n_part, ", ", round(s/n_part * 100), "%"))
        }
        if (s != n_part) {
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):(each * s)]
        } else (
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):N_pred_all]
        )
      }

      # Initialize result containers for each chunk.
      pGRN_list = base::list()
      for (s in 1:n_part) {
        if (s %% idx ==1) {
          message(base::paste0("Initializing result matrices progress: ", s, " / ", n_part, ", ", base::round(s/n_part * 100), "%"))
        }
        pGRN_list[[s]] = base::list()
        pGRN_list[[s]][["theta_p"]] = base::data.frame(base::matrix(
          0, nrow = base::length(re_pred_list[[s]]), ncol = base::length(re_train),
          dimnames = base::list(re_pred_list[[s]], re_train))
        )
        pGRN_list[[s]][["loss"]] = base::data.frame(base::matrix(
          Inf, nrow = base::length(re_pred_list[[s]]), ncol = base::length(re_train),
          dimnames = base::list(re_pred_list[[s]], re_train))
        )
      }; rm(idx)

      # Extract common feature columns and parameters shared by all pairs in the branch.
      message("Inferring GRN of branch ", n, " of cell type ", cell_type, ".")
      state_data = interest_cell_type_branch_data_for_GRN_infer[[cell_type]][[base::paste0("branch", n)]]
      TFE_col = base::grep("TFE_", base::colnames(state_data), value = TRUE)
      TGA_col = base::grep("TGA_", base::colnames(state_data), value = TRUE)
      TGE_col = base::grep("TGE_", base::colnames(state_data), value = TRUE)
      Length = base::length(TFE_col)
      seqLength = base::seq_len(base::length(TFE_col) - 1)
      train_params = lapply(interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]], \(sublist) sublist[["best_params"]])
      Tslot_K = interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]]

      # Core worker function: process one chunk of TF-TG pairs.
      process_data_slice = function(pGRN = base::list(theta_p = NA, loss = NA),
                                    RE_PRED = re_pred_list[[s]],
                                    ncols = base::length(re_train),
                                    STATE_DATA = state_data[re_pred_list[[s]], -(1:3)],
                                    TFE_col = TFE_col,
                                    TGA_col = TGA_col,
                                    TGE_col = TGE_col,
                                    Length = Length,
                                    seqLength = seqLength,
                                    re_train = re_train,
                                    train_params = train_params,
                                    tau = tau,
                                    Tslot_K = Tslot_K) {
        for (re in RE_PRED) {
          # Temporary containers for current TF-TG pair.
          theta_p = base::data.frame(base::matrix(0, nrow = 1, ncol = ncols, dimnames = base::list(re, re_train)))
          loss = base::data.frame(base::matrix(Inf, nrow = 1, ncol = ncols, dimnames = base::list(re, re_train)))

          # Extract time-series expression features for current pair.
          TFE = base::as.numeric(STATE_DATA[re, TFE_col])
          TGA = base::as.numeric(STATE_DATA[re, TGA_col])
          TGE = base::as.numeric(STATE_DATA[re, TGE_col])

          # Predict with each trained model.
          for (re_tr in re_train) {
            pa_best = train_params[[re_tr]]

            # Optimize theta to minimize prediction loss via Brent's method.
            predict_result = stats::optim(par = tau, fn = function(theta) {
              R = scDERI::R_cal(
                theta_TF_TG = theta, tau = tau, r1 = 200,
                r2 = pa_best["r.r2"],
                r3 = 1, r4 = 1000, r5 = 1
              )

              a1 = pa_best["alpha.alpha1"]
              a2 = pa_best["alpha.alpha2"]

              b1 = pa_best["beta.beta1"]
              b2 = pa_best["beta.beta2"]
              b3 = pa_best["beta.beta3"]

              s1 = pa_best["s.s1"]
              s2 = pa_best["s.s2"]

              u11 = pa_best["u.u11"]
              u21 = pa_best["u.u21"]
              u311 = pa_best["u.u311"]
              u321 = pa_best["u.u321"]
              u331 = pa_best["u.u331"]

              u12 = pa_best["u.u12"]
              u22 = pa_best["u.u22"]
              u312 = pa_best["u.u312"]
              u322 = pa_best["u.u322"]
              u332 = pa_best["u.u332"]

              TGE_tilde = pa_best["E~.E~_TG"]
              TFE_tilde = pa_best["E~.E~_TF"]

              U1_tilde = scDERI::Hill_cal(x = TGE_tilde, Dissociation_Constant = u11, Hill_Coefficient = u12)
              U2_tilde = scDERI::Hill_cal(x = TFE_tilde, Dissociation_Constant = u21, Hill_Coefficient = u22)
              U1 = scDERI::Hill_cal(x = TGE[-Length], Dissociation_Constant = u11, Hill_Coefficient = u12)
              U2 = scDERI::Hill_cal(x = TFE[-Length], Dissociation_Constant = u21, Hill_Coefficient = u22)
              U31 = scDERI::Hill_cal(x = TGA[-Length], Dissociation_Constant = u311, Hill_Coefficient = u312)
              U32 = scDERI::Hill_cal(x = TGE[-Length], Dissociation_Constant = u321, Hill_Coefficient = u322)
              U33 = scDERI::Hill_cal(x = TGE_tilde, Dissociation_Constant = u331, Hill_Coefficient = u332)

              v1 = pa_best["v.v1"]
              v2 = pa_best["v.v2"]
              v3 = pa_best["v.v3"]

              base::mean(base::sapply(2:Length, function(K) {
                prev_idx = K - 1
                2 * (TFE[prev_idx] + a1 * R * scDERI::S_cal(s = s1, U = U1[prev_idx], U_tilde = U1_tilde) + b1 - TFE[K])^2 +
                  2 * (TGA[prev_idx] + a2 * R * scDERI::S_cal(s = s2, U = U2[prev_idx], U_tilde = U2_tilde) + b2 - TGA[K])^2 +
                  (TGE[prev_idx] + R * (v1 * U31[prev_idx] - v2 * U32[prev_idx] - v3 * U33) * Tslot_K[K] + b3 - TGE[K])^2
              }))

            }, lower = 0.01, upper = 0.99, method = "Brent")
            theta_p[re, re_tr] = predict_result$par
            loss[re, re_tr] = predict_result$value
          }

          # Filter out columns containing "." to retain valid model results.
          loss = loss[!base::grepl(".", base::names(loss), fixed = TRUE)]
          theta_p = theta_p[!base::grepl(".", base::names(theta_p), fixed = TRUE)]

          pGRN[[re]] = base::list(theta_p = theta_p, loss = loss)
        }

        # Merge pair-level results into two unified data frames for the chunk.
        pGRN[["theta_p"]] = base::as.data.frame(data.table::rbindlist(base::lapply(
          pGRN[3:base::length(pGRN)], function(x) x[["theta_p"]]
        )))
        base::rownames(pGRN[["theta_p"]]) = base::names(pGRN[3:base::length(pGRN)])

        pGRN[["loss"]] = base::as.data.frame(data.table::rbindlist(base::lapply(
          pGRN[3:base::length(pGRN)], function(x) x[["loss"]]
        )))
        base::rownames(pGRN[["loss"]]) = base::names(pGRN[3:base::length(pGRN)])

        return(pGRN)
      }

      # Parallel computation over all data chunks.
      prediction_list = parallel::mclapply(
        X = 1:n_part,
        FUN = function(s) {process_data_slice(
          pGRN = base::list(theta_p = NA, loss = NA),
          RE_PRED = re_pred_list[[s]],
          ncols = base::length(re_train),
          STATE_DATA = state_data[re_pred_list[[s]], -(1:3)],
          TFE_col = TFE_col,
          TGA_col = TGA_col,
          TGE_col = TGE_col,
          Length = Length,
          seqLength = seqLength,
          re_train = re_train,
          train_params = train_params,
          tau = tau,
          Tslot_K = Tslot_K
        )},
        mc.cores = ncores
      ); gc()

      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]] = base::list()

      # Merge theta_p results across all chunks.
      message("Summarizing theta_p of this branch.")
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["theta_p"]] =
        base::do.call("rbind", base::lapply(prediction_list, function(x) x[["theta_p"]]))

      # Merge loss results across all chunks.
      message("Summarizing loss of this branch.")
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["loss"]] =
        base::do.call("rbind", base::lapply(prediction_list, function(x) x[["loss"]]))

      # Save branch-level result as RDS file.
      message("Saving the predictive results of this branch.")
      base::saveRDS(
        interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]],
        file = base::paste0("./interest_cell_type_pGRN.", base::gsub(" ", "_", cell_type), ".branch", n, ".rds")
      )
    }
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Inferring GRN. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_pGRN)
}

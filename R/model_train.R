#' @title Train scDERI Gene Regulatory Network Models
#'
#' @description Trains gene regulatory network models for each TF-TG pair across cell types and branches
#' using a hybrid optimization strategy combining genetic algorithm global search and gradient ascent local refinement.
#' Supports multiple early stopping criteria and parallel computing at the regulon level.
#'
#' @param interest_cell_type_branch_training_set List. Output from \code{\link{select_training_set}},
#'   containing branch-level training data with expression and activity features.
#' @param interest_cell_type_branch_init_params List. Output from \code{\link{set_init_params}},
#'   containing initialized parameters and bounds for each TF-TG pair.
#' @param interest_cell_type_group List. Output from \code{\link{cell_grouping}},
#'   containing pseudotime interval lengths for each branch.
#' @param interest_cell_type_tau Named numeric vector. Output from \code{\link{get_cell_type_tau}},
#'   cell-type-specific binarization thresholds.
#' @param interest_cell_type_genes_pseudotime_info List. Output from \code{\link{get_genes_pseudotime_info}},
#'   containing cell counts per branch for weighted integration of results.
#' @param max_epochs Positive integer. Maximum number of training epochs. Defaults to 1000.
#' @param early_stop Positive integer. Integer. Number of consecutive epochs without loss improvement to trigger early stopping. Defaults to 10.
#' @param eps_theta Tiny positive numeric. Convergence threshold for predicted regulatory strength versus observed strength. Defaults to 1e-3.
#' @param eps_loss Tiny positive numeric. Threshold for negligible loss improvement. Defaults to 1e-5.
#' @param popSize_ga Integer. Population size for genetic algorithm. Defaults to 50. See \code{popSize} in \code{\link[GA]{ga}}.
#' @param maxiter_ga Positive integer. Maximum generations per genetic algorithm run. Defaults to 150. See \code{maxiter} in \code{\link[GA]{ga}}.
#' @param pcrossover_ga Numeric between 0 and 1. Crossover probability for genetic algorithm. Defaults to 0.9. See \code{pcrossover} in \code{\link[GA]{ga}}.
#' @param pmutation_ga Numeric between 0 and 1. Mutation probability for genetic algorithm, between 0 and 1. Defaults to 0.5. See \code{pmutation} in \code{\link[GA]{ga}}.
#' @param parallel_ga Logical. Whether to enable parallel computing within genetic algorithm. Defaults to \code{FALSE}. See \code{parallel} in \code{\link[GA]{ga}}.
#' @param seed_ga Integer. Random seed for reproducible genetic algorithm results. Defaults to 123. See \code{seed} in \code{\link[GA]{ga}}.
#' @param iterations_grad Positive integer. Number of gradient ascent steps per epoch. Defaults to 30.
#' @param ncores_grad Positive integer. Number of cores for parallel gradient computation. Defaults to 1.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#'
#' @return Nested list organized by cell type, each entry contains:
#' \itemize{
#'   \item \code{best_pred_result}: Integrated data frame across all branches, including predicted regulatory strength,
#'     binarized results, fitting loss, and training metadata.
#'   \item \code{branch*}: List of detailed training results for each regulon in the corresponding branch,
#'     including training data, pseudotime intervals, threshold, objective values, parameter trajectories,
#'     loss history, best epoch, best parameters, best prediction, stop cause, training time, and total epochs.
#' }
#'
#' @details For each cell type and branch, the function independently trains a model for each TF-TG pair:
#' \enumerate{
#'   \item Initializes genetic algorithm population with pre-set parameter values.
#'   \item Iterates through training epochs with a two-stage hybrid optimization:
#'     \itemize{
#'       \item Stage 1: Genetic algorithm for global parameter space exploration, initialized with the optimal solution from the previous epoch.
#'       \item Stage 2: Gradient ascent for local refinement based on the best GA solution.
#'     }
#'   \item Checks early stopping conditions in priority order after each epoch:
#'     maximum epochs reached, no loss improvement for consecutive epochs, predicted theta converges to observed theta, negligible loss improvement.
#'   \item Integrates results across all branches: binary regulatory relationships use union logic,
#'     continuous values are weighted by branch cell counts.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and temporarily switches to the \code{4.2 Build Prediction Model} directory,
#'     and restores the original working directory on exit.
#'   \item Each regulon is trained independently, enabling parallel computing at the regulon level via \code{ncores}.
#'   \item When \code{parallel_ga = TRUE}, an additional serial GA iteration is performed after parallel GA
#'     to mitigate result instability.
#'   \item Adjusting population size, maximum iterations and mutation probability can help the algorithm
#'     escape local optima but will increase runtime.
#' }
#'
#' @seealso \code{\link{f_train_cal}}, \code{\link{gradient_ascent_train}}, \code{\link{loss_cal}}, \code{\link{theta_p_cal}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_training_set = base::readRDS(
#'   "./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds"
#' )
#' interest_cell_type_branch_init_params = base::readRDS(
#'   "./4.2 Build Prediction Model/interest_cell_type_branch_init_params.rds"
#' )
#' interest_cell_type_group = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
#' interest_cell_type_genes_pseudotime_info = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#'
#' max_epochs = 50
#' early_stop = 5
#'
#' popSize_ga = 50
#' maxiter_ga = 150
#' pcrossover_ga = 0.9
#' pmutation_ga = 0.5
#' parallel_ga = FALSE
#' seed_ga = 123 # Reproducible
#' # seed_ga = stats::runif(1) # Not Reproducible
#'
#' iterations_grad = 30
#' ncores_grad = 1
#'
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#'
#' interest_cell_type_branch_model_train = model_train(
#'   interest_cell_type_branch_training_set = interest_cell_type_branch_training_set,
#'   interest_cell_type_branch_init_params = interest_cell_type_branch_init_params,
#'   interest_cell_type_group = interest_cell_type_group,
#'   interest_cell_type_tau = interest_cell_type_tau,
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   max_epochs = max_epochs,
#'   early_stop = early_stop,
#'   popSize_ga = popSize_ga,
#'   maxiter_ga = maxiter_ga,
#'   pcrossover_ga = pcrossover_ga,
#'   pmutation_ga = pmutation_ga,
#'   parallel_ga = parallel_ga,
#'   seed_ga = 123,
#'   iterations_grad = iterations_grad,
#'   ncores_grad = ncores_grad,
#'   ncores = ncores
#' )
#' base::saveRDS(
#'   interest_cell_type_branch_model_train,
#'   file = "./4.2 Build Prediction Model/interest_cell_type_branch_model_train.rds"
#' )
#' }
model_train = function (
    interest_cell_type_branch_training_set,
    interest_cell_type_branch_init_params,
    interest_cell_type_group,
    interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info,
    max_epochs = 1000, early_stop = 10,
    eps_theta = 1e-3, eps_loss = 1e-5,
    popSize_ga = 50, maxiter_ga = 150, pcrossover_ga = 0.9, pmutation_ga = 0.5, parallel_ga = FALSE, seed_ga = 123,
    iterations_grad = 30, ncores_grad = 1,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Training models for each branch. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "4.2 Build Prediction Model"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Process.
  interest_cell_type_branch_model_train = base::list()
  for (cell_type in base::names(interest_cell_type_branch_training_set)) {
    message("Starting to train prediction models for all branches in cell type: ", cell_type)

    message("Initializing model training iteration result storage.")
    interest_cell_type_branch_model_train[[cell_type]] = base::list()

    message("Initializing optimal prediction result storage.")
    n_branch = base::length(interest_cell_type_branch_training_set[[cell_type]]) - 1 # Subtract the "all" element to get the number of branches
    best_pred_result = interest_cell_type_branch_training_set[[cell_type]][["all"]]
    best_pred_result["theta_p"] = 0
    best_pred_result[base::paste0("theta_p_branch", 1:n_branch)] = 0
    best_pred_result$theta_i_bin = base::as.integer(best_pred_result$theta_i > 0)
    best_pred_result["theta_p_bin"] = 0
    best_pred_result[base::paste0("theta_p_bin_branch", 1:n_branch)] = 0
    best_pred_result["fit_loss"] = Inf
    best_pred_result[base::paste0("fit_loss_branch", 1:n_branch)] = Inf
    interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]] = best_pred_result

    message("Iterating through each branch and training each branch independently.")
    for (n in 1:n_branch) {
      message("Starting training for branch ", n, " of cell type ", cell_type)
      message("There are ", nrow(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]]), " regulons for training.")
      n_re = 0
      interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]] = base::list()

      # Train each regulon.
      message("Training each regulon.")
      model_train_process = function(re,
                                     n = n,
                                     n_re = base::match(re, base::names(interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]])),
                                     N_re = base::nrow(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]]),
                                     N_branch = n_branch,
                                     branch_params_list = interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]][[re]],
                                     train_data = base::as.matrix(interest_cell_type_branch_training_set[[cell_type]][[n]][re, , drop = FALSE]),
                                     Tslot_K = interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]],
                                     tau = interest_cell_type_tau[cell_type],
                                     popSize_ga = popSize_ga,
                                     max_epochs = max_epochs,
                                     pcrossover_ga = pcrossover_ga,
                                     pmutation_ga = pmutation_ga,
                                     maxiter_ga = maxiter_ga,
                                     parallel_ga = parallel_ga,
                                     seed_ga = seed_ga,
                                     best_pred_result_re = best_pred_result[re, , drop = FALSE]) {
        result = base::list()
        result[["re"]] = re

        signal = base::paste0(n_re, "/", N_re)
        message("Progress ", signal, "...")
        Signal = base::paste0("[Branch ", n, "/", N_branch, " - Train Regulon ", signal, "] ")
        message(Signal, "Extracting initial parameter values, upper bounds, lower bounds, and names.")
        params = base::unlist(base::lapply(branch_params_list, function(category) {
          base::lapply(category, function(param) param$initialization)
        }))
        params_upper = base::unlist(base::lapply(branch_params_list, function(category) {
          base::lapply(category, function(param) param$upper)
        }))
        params_lower = base::unlist(base::lapply(branch_params_list, function(category) {
          base::lapply(category, function(param) param$lower)
        }))
        params_init = base::unlist(base::lapply(branch_params_list, function(category) {
          base::lapply(category, function(param) param$initialization)
        }))
        params_names = base::names(params)

        message(Signal, "Extracting training set.")
        result[["train_data"]] = train_data

        message(Signal, "Extracting the pseudotime interval length for all cell groups.")
        result[["Tslot_K"]] = Tslot_K

        message(Signal, "Extracting the binary classification threshold for regulatory strength.")
        result[["tau"]] = tau

        message(Signal, "Initializing objective function with a large default value.")
        f_value_train = 1e+5 + 1

        message(Signal, "Initializing the genetic algorithm population.")
        prev_population = base::matrix(
          data = base::rep(base::as.numeric(params), popSize_ga),
          nrow = popSize_ga, ncol = base::length(params), byrow = TRUE
        )

        message(Signal, "Initializing the model parameters storage for each training round.")
        para_data_frame = base::as.data.frame(base::matrix(
          data = 0, nrow = max_epochs, ncol = base::length(params),
          dimnames = base::list(1:max_epochs, params_names)
        ))

        message(Signal, "Initializing the loss storage for each training round.")
        loss_train = base::rep(0, max_epochs)
        names(loss_train) = 1:max_epochs
        f_value_train_vec = base::rep(0, max_epochs)
        names(f_value_train_vec) = 1:max_epochs

        t1_re = base::Sys.time()
        for(epoch in 1:max_epochs) {
          # Each epoch workflow:
          # 1. Genetic algorithm for global exploration (initialized with previous epoch's optimal parameters)
          # 2. Gradient ascent for local refinement based on GA's best solution
          # 3. Use refined parameters as the initial population for next epoch's GA
          SIGNAL = base::paste0(base::gsub("] ", "", base::paste0(Signal, " - Epoch ", epoch), fixed = TRUE), "] ")
          if (epoch %% 10 == 1) {message(SIGNAL, "Epoch ", epoch, "/", max_epochs, " in progress")}

          # Genetic Algorithm Training
          result_ga = GA::ga(
            type = "real-valued", fitness = scDERI::f_train_cal, lower = params_lower, upper = params_upper,
            popSize = popSize_ga, pcrossover = pcrossover_ga, pmutation = pmutation_ga,
            elitism = base::ceiling(base::min(0.2 * maxiter_ga, popSize_ga / 5)),
            maxiter = maxiter_ga, suggestions = prev_population, parallel = parallel_ga, seed = seed_ga, monitor = FALSE,
            params_names = params_names, train_data = train_data, tau = tau, Tslot_K = Tslot_K,
            params_init = params_init, params_upper = params_upper, params_lower = params_lower
          )
          # Parallel GA may introduce result instability; rebuild population with the best solution,
          # run 1 serial GA iteration for stabilization, then free memory via garbage collection
          if (parallel_ga) {
            prev_population = base::matrix(
              data = base::rep(base::as.numeric(result_ga@solution[1, ]), popSize_ga),
              nrow = popSize_ga, ncol = base::length(params), byrow = TRUE
            )
            result_ga = GA::ga(
              type = "real-valued", fitness = scDERI::f_train_cal, lower = params_lower, upper = params_upper,
              popSize = popSize_ga, pcrossover = pcrossover_ga, pmutation = pmutation_ga,
              elitism = base::ceiling(base::min(0.2 * maxiter_ga, popSize_ga / 5)),
              maxiter = 1, suggestions = prev_population, parallel = FALSE, seed = seed_ga, monitor = FALSE,
              params_names = params_names, train_data = train_data, tau = tau, Tslot_K = Tslot_K,
              params_init = params_init, params_upper = params_upper, params_lower = params_lower
            ); base::gc()
          }
          f_value_train = -result_ga@fitnessValue

          # Gradient Ascent Refinement
          result_gradient = scDERI::gradient_ascent_train(
            iterations = iterations_grad, alpha_lower = 1e-5, alpha_upper = 1e-3, alpha_guess = 1e-4,
            init_theta = base::as.numeric(result_ga@solution[1, ]),
            h_max = 1e-9, h_min_percent = 0.01, ncores = ncores_grad,
            params_lower = params_lower, params_upper = params_upper, monitor = FALSE,
            params_names = params_names, train_data = train_data, tau = tau, Tslot_K = Tslot_K,
            params_init = params_init
          )
          f_value_train = -result_gradient$objective_value
          f_value_train_vec[epoch] = f_value_train
          result[["f_value_train"]] = f_value_train_vec
          prev_population = base::matrix(
            data = base::rep(result_gradient$parameters, popSize_ga),
            nrow = popSize_ga, ncol = base::length(params), byrow = TRUE
          )
          para_data_frame[epoch, ] = result_gradient$parameters

          result[["para_data_frame"]] = para_data_frame

          # Calculating training loss
          loss_train[epoch] = scDERI::loss_cal(
            data = train_data,
            params = result_gradient$parameters, params_names = params_names,
            tau = tau, Tslot_K = Tslot_K
          )

          # Early stopping criteria (checked in priority order):
          # 1. Reach maximum training epochs (max_epochs)
          # 2. No loss improvement for consecutive `early_stop` epochs
          # 3. Predicted regulatory strength (theta_p) converges to observed theta_i
          # 4. Loss improvement amplitude is below the threshold eps_loss
          if (epoch >= max_epochs) {
            message(SIGNAL, "Reaching the maximum number of training epochs.")
            loss_vec = loss_train[1:base::ifelse(
              base::length(base::which(loss_train != 0)) > 0, base::max(base::which(loss_train != 0)), 0
            )]
            result[["loss_train"]] = loss_vec

            best_epoch = base::as.numeric(base::names(base::which.min(loss_vec)))
            result[["best_epoch"]] = best_epoch

            best_params = base::as.numeric(para_data_frame[best_epoch, , drop = FALSE])
            base::names(best_params) = params_names
            result[["best_params"]] = best_params

            best_pred_result_re[re, base::paste0("fit_loss_branch", n)] = loss_vec[best_epoch]
            best_pred_result_re[re, base::paste0("theta_p_branch", n)] = scDERI::theta_p_cal(
              data = train_data,
              best_params = best_params, params_names = params_names,
              tau = tau, Tslot_K = Tslot_K
            )
            best_pred_result_re[re, base::paste0("theta_p_bin_branch", n)] = base::ifelse(best_pred_result_re[re, base::paste0("theta_p_branch", n)] >= tau, 1, 0)
            result[["best_pred_result_re"]] = best_pred_result_re

            stop_cause = "maximum training epochs reached"
            message(SIGNAL, "Meeting the early stop condition: ", stop_cause, ".")
            result[["stop_cause"]] = stop_cause
            break
          } else if (epoch >= early_stop) {
            loss_vec = loss_train[1:base::ifelse(
              base::length(base::which(loss_train != 0)) > 0, base::max(base::which(loss_train != 0)), 0
            )]
            result[["loss_train"]] = loss_vec

            best_epoch = base::as.numeric(base::names(base::which.min(loss_vec)))
            result[["best_epoch"]] = best_epoch

            best_params = base::as.numeric(para_data_frame[best_epoch, , drop = FALSE])
            base::names(best_params) = params_names
            result[["best_params"]] = best_params

            best_pred_result_re[re, base::paste0("fit_loss_branch", n)] = loss_vec[best_epoch]
            best_pred_result_re[re, base::paste0("theta_p_branch", n)] = scDERI::theta_p_cal(
              data = train_data,
              best_params = best_params, params_names = params_names,
              tau = tau, Tslot_K = Tslot_K
            )
            best_pred_result_re[re, base::paste0("theta_p_bin_branch", n)] = base::ifelse(best_pred_result_re[re, base::paste0("theta_p_branch", n)] >= tau, 1, 0)
            result[["best_pred_result_re"]] = best_pred_result_re

            if (best_epoch <= base::length(loss_vec) - early_stop + 1) {
              stop_cause = "no loss improvement for consecutive early_stop epochs"
              message(SIGNAL, "Meeting the early stop condition: ", stop_cause, ".")
              result[["stop_cause"]] = stop_cause
              break
            } else if (base::abs(best_pred_result_re[re, base::paste0("theta_p_branch", n)] - best_pred_result_re[re, "theta_i"]) < eps_theta) {
              stop_cause = "predicted theta_p converges to observed theta_i"
              message(SIGNAL, "Meeting the early stop condition: ", stop_cause, ".")
              result[["stop_cause"]] = stop_cause
              break
            } else if (base::max(loss_vec[(base::length(loss_vec) - early_stop + 1):base::length(loss_vec)]) - loss_vec[best_epoch] < eps_loss) {
              stop_cause = "loss improvement is below the threshold"
              message(SIGNAL, "Meeting the early stop condition: ", stop_cause, ".")
              result[["stop_cause"]] = stop_cause
              break
            }
          }
        }
        t2_re = base::Sys.time()
        message(Signal, "Training completed. Duration: ", base::format(t2_re - t1_re))

        result[["train_time"]] = base::sprintf("%.2f %s", base::as.numeric(t2_re - t1_re), base::units(t2_re - t1_re))
        result[["train_epoch"]] = epoch

        return(result)
      }

      interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]] = parallel::mclapply(
        X = base::rownames(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]]),
        FUN = function(re) {
          model_train_process(
            re,
            n = n,
            n_re = base::match(re, base::names(interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]])),
            N_re = base::nrow(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]]),
            N_branch = n_branch,
            branch_params_list = interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]][[re]],
            train_data = base::as.matrix(interest_cell_type_branch_training_set[[cell_type]][[n]][re, , drop = FALSE]),
            Tslot_K = interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]],
            tau = interest_cell_type_tau[cell_type],
            popSize_ga = popSize_ga,
            max_epochs = max_epochs,
            pcrossover_ga = pcrossover_ga,
            pmutation_ga = pmutation_ga,
            maxiter_ga = maxiter_ga,
            parallel_ga = parallel_ga,
            seed_ga = seed_ga,
            best_pred_result_re = best_pred_result[re, , drop = FALSE]
          )
        },
        mc.cores = ncores
      )
      base::names(interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]]) =
        base::rownames(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]])

      for (re in base::rownames(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]])) {
        best_pred_result[re, base::paste0("theta_p_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["best_pred_result_re"]][re, base::paste0("theta_p_branch", n)]
        best_pred_result[re, base::paste0("theta_p_bin_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["best_pred_result_re"]][re, base::paste0("theta_p_bin_branch", n)]
        best_pred_result[re, base::paste0("fit_loss_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["best_pred_result_re"]][re, base::paste0("fit_loss_branch", n)]
      }
    }
    message("All branch trainings have been completed.")

    # Multi-branch integration: binary regulatory relationships use union logic; continuous values are weighted by cell count per branch
    message("Integrating the results from all branches.")
    if (n_branch > 1) {
      best_pred_result$theta_p_bin = base::as.integer(base::rowSums(
        best_pred_result[, base::grep("^theta_p_bin_branch\\d+$", base::names(best_pred_result))]) > 0
      )

      weights = base::unlist(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]])
      best_pred_result$theta_p = base::apply(best_pred_result[, base::grep("^theta_p_branch([0-9]{1,2})$", base::names(best_pred_result), value = TRUE)], 1, function(x) {
        non_zero_idx = which(x != 0)
        if (base::length(non_zero_idx) == 0) 0 else {
          base::sum(x[non_zero_idx] * weights[non_zero_idx]) / base::sum(weights[non_zero_idx])
        }
      })
      best_pred_result$fit_loss = base::apply(best_pred_result[, base::grep("^fit_loss_branch([0-9]{1,2})$", base::names(best_pred_result), value = TRUE)], 1, function(x) {
        non_zero_idx = which(x != Inf)
        if (base::length(non_zero_idx) == 0) Inf else {
          base::sum(x[non_zero_idx] * weights[non_zero_idx]) / base::sum(weights[non_zero_idx])
        }
      })
    } else if (n_branch == 1) {
      best_pred_result$theta_p_bin = best_pred_result$theta_p_bin_branch1

      best_pred_result$theta_p = best_pred_result$theta_p_branch1
      best_pred_result$fit_loss = best_pred_result$fit_loss_branch1
    } else {
      stop("The number of branches is less than 1.")
    }

    best_pred_result = base::cbind(
      best_pred_result[, 1:4],
      base::abs(best_pred_result$theta_i - best_pred_result$theta_p),
      best_pred_result[, 5:base::ncol(best_pred_result)]
    )
    base::colnames(best_pred_result)[5] = "theta_difference"

    best_pred_result[base::paste0("theta_difference_branch", 1:n_branch)] = Inf
    best_pred_result[base::paste0("train_time_branch", 1:n_branch)] = 0
    best_pred_result[base::paste0("train_epoch_branch", 1:n_branch)] = 0
    best_pred_result[base::paste0("best_epoch_branch", 1:n_branch)] = 0
    best_pred_result[base::paste0("stop_cause_branch", 1:n_branch)] = NA
    for (n in 1:n_branch) {
      for (re in base::rownames(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]])) {
        best_pred_result[re, base::paste0("theta_difference_branch", n)] =
          base::abs(best_pred_result[re, "theta_i"] - best_pred_result[re, base::paste0("theta_p_branch", n)])
        best_pred_result[re, base::paste0("train_time_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["train_time"]]
        best_pred_result[re, base::paste0("train_epoch_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["train_epoch"]]
        best_pred_result[re, base::paste0("best_epoch_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["best_epoch"]]
        best_pred_result[re, base::paste0("stop_cause_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["stop_cause"]]
      }
    }

    interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]] = best_pred_result
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Training models for each branch. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_branch_model_train)
}

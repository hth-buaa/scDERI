#' @title Initialize Parameters for scDERI Gene Regulatory Network Model
#'
#' @description Generates initialized parameter structures for the scDERI gene regulatory network model for each TF-TG pair across cell types and branches.
#' Each parameter is assigned an initial value, upper bound and lower bound for subsequent optimization, with baseline expression parameters adaptively initialized from training set means.
#'
#' @param interest_cell_type_branch_training_set Named list.
#' Output from the \code{\link{select_training_set}} function,
#' containing branch-level training sets with expression and activity features.
#'
#' @return Nested list with the hierarchy \code{cell type → branch → TF-TG interaction}.
#' The innermost layer contains 7 categories of model parameters,
#' each with three elements: \code{initialization}, \code{upper} and \code{lower}.
#' The merged "all" dataset is not included in the output.
#'
#' @details The workflow is executed independently for each cell type:
#' \enumerate{
#'   \item Builds a hierarchical parameter structure for each branch, excluding the merged "all" dataset at the end of the input list.
#'   \item Initializes 7 categories of scDERI model parameters for each TF-TG pair:
#'     \itemize{
#'       \item \code{r}: effect efficiency function parameters (only r2 is trained)
#'       \item \code{alpha}: effect conversion rate parameters
#'       \item \code{beta}: linear correction term parameters
#'       \item \code{s}: power constant parameters
#'       \item \code{u}: Hill-type non-linear transformation parameters, including dissociation constants and Hill coefficients
#'       \item \code{E~}: environmental baseline expression level parameters
#'       \item \code{v}: rate parameters for maximum production, spontaneous degradation and usage-induced degradation
#'     }
#'   \item Adaptively initializes baseline parameters \code{E~_TG} and \code{E~_TF} using the mean TGE and TFE values across time points from the training set, and sets corresponding bounds accordingly.
#'   \item Presets initial values, upper bounds and lower bounds for all parameters to constrain the parameter space of subsequent optimization.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates a \code{4.2 Build Prediction Model} directory and temporarily switches the working directory into it; the original path is restored upon completion.
#'   \item Baseline expression parameters are data-adaptively initialized to improve optimization convergence efficiency.
#'   \item Parameter bounds are preset based on biological interpretability and numerical stability, and generally do not need manual adjustment.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_training_set = base::readRDS(
#'   "./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds"
#' )
#' interest_cell_type_branch_init_params = set_init_params(
#'   interest_cell_type_branch_training_set = interest_cell_type_branch_training_set
#' )
#' base::saveRDS(
#'   interest_cell_type_branch_init_params,
#'   file = "./4.2 Build Prediction Model/interest_cell_type_branch_init_params.rds"
#' )
#' }
set_init_params = function(interest_cell_type_branch_training_set) {
  ### Start.
  t_start = base::Sys.time()
  message("Run: Setting initial parameter values for all branches. Time: ", t_start, ".")
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

  ### Inner function: initialize parameters for a single TF-TG pair
  set_init_params_process = function(TFE_mean, TGE_mean) {
    message("Initializing with a list.")
    init_params = base::list()

    # message("Setting r.")
    init_params[["r"]] = base::list()
    message("Setting r2 in r.")
    init_params[["r"]][["r2"]] = base::list()
    init_params[["r"]][["r2"]][["initialization"]] = 1
    init_params[["r"]][["r2"]][["upper"]] = 1.113383
    init_params[["r"]][["r2"]][["lower"]] = 0.370506

    # message("Setting alpha.")
    init_params[["alpha"]] = base::list()
    message("Setting alpha1 in alpha.")
    init_params[["alpha"]][["alpha1"]] = base::list()
    init_params[["alpha"]][["alpha1"]][["initialization"]] = 1
    init_params[["alpha"]][["alpha1"]][["upper"]] = 50
    init_params[["alpha"]][["alpha1"]][["lower"]] = 0.1
    message("Setting alpha2 in alpha.")
    init_params[["alpha"]][["alpha2"]] = base::list()
    init_params[["alpha"]][["alpha2"]][["initialization"]] = 1
    init_params[["alpha"]][["alpha2"]][["upper"]] = 50
    init_params[["alpha"]][["alpha2"]][["lower"]] = 0.1

    # message("Setting beta.")
    init_params[["beta"]] = base::list()
    message("Setting beta1 in beta.")
    init_params[["beta"]][["beta1"]] = base::list()
    init_params[["beta"]][["beta1"]][["initialization"]] = 0
    init_params[["beta"]][["beta1"]][["upper"]] = 10
    init_params[["beta"]][["beta1"]][["lower"]] = -10
    message("Setting beta2 in beta.")
    init_params[["beta"]][["beta2"]] = base::list()
    init_params[["beta"]][["beta2"]][["initialization"]] = 0
    init_params[["beta"]][["beta2"]][["upper"]] = 10
    init_params[["beta"]][["beta2"]][["lower"]] = -10
    message("Setting beta3 in beta.")
    init_params[["beta"]][["beta3"]] = base::list()
    init_params[["beta"]][["beta3"]][["initialization"]] = 0
    init_params[["beta"]][["beta3"]][["upper"]] = 10
    init_params[["beta"]][["beta3"]][["lower"]] = -10

    # message("Setting s.")
    init_params[["s"]] = base::list()
    message("Setting s1 in s.")
    init_params[["s"]][["s1"]] = base::list()
    init_params[["s"]][["s1"]][["initialization"]] = 1
    init_params[["s"]][["s1"]][["upper"]] = 4.9
    init_params[["s"]][["s1"]][["lower"]] = 0.1
    message("Setting s2 in s.")
    init_params[["s"]][["s2"]] = base::list()
    init_params[["s"]][["s2"]][["initialization"]] = 1
    init_params[["s"]][["s2"]][["upper"]] = 4.9
    init_params[["s"]][["s2"]][["lower"]] = 0.1

    # message("Setting u.")
    init_params[["u"]] = base::list()
    message("Setting dissociation constant in u.")
    init_params[["u"]][["u11"]] = base::list()
    init_params[["u"]][["u11"]][["initialization"]] = 1
    init_params[["u"]][["u11"]][["upper"]] = 777600
    init_params[["u"]][["u11"]][["lower"]] = 0.1
    init_params[["u"]][["u21"]] = base::list()
    init_params[["u"]][["u21"]][["initialization"]] = 1
    init_params[["u"]][["u21"]][["upper"]] = 777600
    init_params[["u"]][["u21"]][["lower"]] = 0.1
    init_params[["u"]][["u311"]] = base::list()
    init_params[["u"]][["u311"]][["initialization"]] = 1
    init_params[["u"]][["u311"]][["upper"]] = 777600
    init_params[["u"]][["u311"]][["lower"]] = 0.1
    init_params[["u"]][["u321"]] = base::list()
    init_params[["u"]][["u321"]][["initialization"]] = 1
    init_params[["u"]][["u321"]][["upper"]] = 777600
    init_params[["u"]][["u321"]][["lower"]] = 0.1
    init_params[["u"]][["u331"]] = base::list()
    init_params[["u"]][["u331"]][["initialization"]] = 1
    init_params[["u"]][["u331"]][["upper"]] = 777600
    init_params[["u"]][["u331"]][["lower"]] = 0.1
    message("Setting Hill coefficient in u.")
    init_params[["u"]][["u12"]] = base::list()
    init_params[["u"]][["u12"]][["initialization"]] = 1
    init_params[["u"]][["u12"]][["upper"]] = 4.9
    init_params[["u"]][["u12"]][["lower"]] = 0.1
    init_params[["u"]][["u22"]] = base::list()
    init_params[["u"]][["u22"]][["initialization"]] = 1
    init_params[["u"]][["u22"]][["upper"]] = 4.9
    init_params[["u"]][["u22"]][["lower"]] = 0.1
    init_params[["u"]][["u312"]] = base::list()
    init_params[["u"]][["u312"]][["initialization"]] = 1
    init_params[["u"]][["u312"]][["upper"]] = 4.9
    init_params[["u"]][["u312"]][["lower"]] = 0.1
    init_params[["u"]][["u322"]] = base::list()
    init_params[["u"]][["u322"]][["initialization"]] = 1
    init_params[["u"]][["u322"]][["upper"]] = 4.9
    init_params[["u"]][["u322"]][["lower"]] = 0.1
    init_params[["u"]][["u332"]] = base::list()
    init_params[["u"]][["u332"]][["initialization"]] = 1
    init_params[["u"]][["u332"]][["upper"]] = 4.9
    init_params[["u"]][["u332"]][["lower"]] = 0.1

    # message("Setting E~.")
    init_params[["E~"]] = base::list()
    message("Setting E~_TG in E~.")
    init_params[["E~"]][["E~_TG"]] = base::list()
    init_params[["E~"]][["E~_TG"]][["initialization"]] = TGE_mean
    init_params[["E~"]][["E~_TG"]][["upper"]] = base::max(9, TGE_mean + 1)
    init_params[["E~"]][["E~_TG"]][["lower"]] = base::max(base::min(1, TGE_mean - 1), 0)
    message("Setting E~_TF in E~.")
    init_params[["E~"]][["E~_TF"]] = base::list()
    init_params[["E~"]][["E~_TF"]][["initialization"]] = TFE_mean
    init_params[["E~"]][["E~_TF"]][["upper"]] = base::max(9, TFE_mean + 1)
    init_params[["E~"]][["E~_TF"]][["lower"]] = base::max(base::min(1, TFE_mean - 1), 0)

    # message("Setting v.")
    init_params[["v"]] = base::list()
    message("Setting v1 in v.")
    init_params[["v"]][["v1"]] = base::list()
    init_params[["v"]][["v1"]][["initialization"]] = 1
    init_params[["v"]][["v1"]][["upper"]] = 200
    init_params[["v"]][["v1"]][["lower"]] = 0.1
    message("Setting v2 in v.")
    init_params[["v"]][["v2"]] = base::list()
    init_params[["v"]][["v2"]][["initialization"]] = 1
    init_params[["v"]][["v2"]][["upper"]] = 200
    init_params[["v"]][["v2"]][["lower"]] = 0.1
    message("Setting v3 in v.")
    init_params[["v"]][["v3"]] = base::list()
    init_params[["v"]][["v3"]][["initialization"]] = 1
    init_params[["v"]][["v3"]][["upper"]] = 200
    init_params[["v"]][["v3"]][["lower"]] = 0.1

    return(init_params)
  }

  ### Generate parameter lists for all cell types and branches
  interest_cell_type_branch_init_params = base::list()
  for (cell_type in base::names(interest_cell_type_branch_training_set)) {
    # Iterate over branch datasets only, exclude the last merged "all" dataset
    message("Setting initial parameter values for the data of cell type ", cell_type, ".")
    interest_cell_type_branch_init_params[[cell_type]] = base::list()
    for (n in 1:(base::length(interest_cell_type_branch_training_set[[cell_type]]) - 1)) {
      interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]] = base::list()
      for (re in rownames(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]])) {
        interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]][[re]] = set_init_params_process(
          TFE_mean = base::mean(base::as.numeric(dplyr::select(interest_cell_type_branch_training_set[[cell_type]][[n]], dplyr::matches("TFE"))[re, ]), na.rm = TRUE),
          TGE_mean = base::mean(base::as.numeric(dplyr::select(interest_cell_type_branch_training_set[[cell_type]][[n]], dplyr::matches("TGE"))[re, ]), na.rm = TRUE)
        )
      }
    }
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Setting initial parameter values for all branches. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_branch_init_params)
}

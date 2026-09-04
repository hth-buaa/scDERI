#' @title Optimize Regulatory Strength Theta for a Single TF-TG Pair
#'
#' @description Fixes all trained dynamic parameters and uses Brent's one-dimensional optimization
#' to find the optimal regulatory strength theta that minimizes prediction error of expression dynamics.
#' It is used to infer the regulatory strength of a given TF-TG pair.
#'
#' @param data Data frame. Input data for one TF-TG pair, containing a \code{theta_i} column
#'   and time-series columns prefixed with \code{TFE_}, \code{TGA_}, \code{TGE_}.
#' @param best_params Numeric vector. Trained and optimized model parameters.
#' @param params_names Character vector. Names corresponding to each parameter element.
#' @param tau Positive numeric. Regulatory binarization threshold.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths.
#'
#' @return Single numeric value. Optimized regulatory strength theta, bounded between 0.01 and 0.99.
#'
#' @details The function treats all dynamic model parameters as fixed and theta as the only optimization variable:
#' \enumerate{
#'   \item Uses the same weighted squared error calculation logic as \code{\link{loss_cal}} as the objective.
#'   \item Performs one-dimensional optimization via Brent's method within the interval \[0.01, 0.99\] via \code{\link[stats]{optim}}.
#'   \item Returns the theta value that minimizes prediction error.
#' }
#'
#' @note
#' \itemize{
#'   \item Typically used for regulatory strength inference after model training,
#'     or for performance evaluation and downstream analysis.
#'   \item The 0.01–0.99 boundary avoids extreme values that may cause numerical instability
#'     in subsequent calculations.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Infer optimal regulatory strength for a TF-TG pair
#' predicted_theta <- theta_p_cal(
#'   data = train_data,
#'   best_params = best_params,
#'   params_names = param_names,
#'   tau = tau,
#'   Tslot_K = Tslot_K
#' )
#' }
theta_p_cal = function (data = train_data,
                        best_params = best_params,
                        params_names = params_names,
                        tau = tau,
                        Tslot_K = Tslot_K)
{
  pa_best = base::as.numeric(best_params)
  base::names(pa_best) = params_names
  stats::optim(par = tau, fn = function(theta) {
    TFE = base::as.numeric(data[, base::grep("TFE_", base::colnames(data), value = TRUE)])
    TGA = base::as.numeric(data[, base::grep("TGA_", base::colnames(data), value = TRUE)])
    TGE = base::as.numeric(data[, base::grep("TGE_", base::colnames(data), value = TRUE)])

    Length = base::length(TFE)
    seqLength = base::seq_len(Length - 1)

    R = base::as.numeric(scDERI::R_cal(
      theta_TF_TG = theta,
      tau = tau,
      r1 = 200,
      r2 = pa_best["r.r2"],
      r3 = 1,
      r4 = 1000,
      r5 = 1
    ))

    alpha1 = pa_best["alpha.alpha1"]
    alpha2 = pa_best["alpha.alpha2"]

    beta1 = (pa_best[base::paste0("beta.beta1_TG_K-1.", seqLength)] + pa_best[base::paste0("beta.beta1_TF_K-1.", seqLength)])/2
    beta2 = (pa_best[base::paste0("beta.beta2_TG_K-1.", seqLength)] + pa_best[base::paste0("beta.beta2_TF_K-1.", seqLength)])/2
    beta3 = pa_best[base::paste0("beta.beta3_TG_K-1.", seqLength)]

    beta1 = pa_best["beta.beta1"]
    beta2 = pa_best["beta.beta2"]
    beta3 = pa_best["beta.beta3"]

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

    base::sum(base::sapply(2:Length, function(K) {
      prev_idx = K - 1
      2 * (TFE[prev_idx] + alpha1 * R * scDERI::S_cal(s = s1, U = U1[prev_idx], U_tilde = U1_tilde) + beta1 - TFE[K])^2 +
        2 * (TGA[prev_idx] + alpha2 * R * scDERI::S_cal(s = s2, U = U2[prev_idx], U_tilde = U2_tilde) + beta2 - TGA[K])^2 +
        (TGE[prev_idx] + R * (v1 * U31[prev_idx] - v2 * U32[prev_idx] - v3 * U33) * Tslot_K[K] + beta3 - TGE[K])^2
    }))
  }, lower = 0.01, upper = 0.99, method = "Brent")$par
}

#include "ctmc_moments_derivatives.h"

#include <utility>

#include "counting_tables.h"


//
// This source file defines the CTMC moment/derivative functions.
//


arma::cube ctmc_moments_Q_derivatives_aux(double t, const arma::mat& Q,
                                          const arma::mat& B, int max_order,
                                          Mode mode) {
  // Compute either CTMC moments or CTMC rate matrix derivatives associated with
  // the branch length `t`.
  arma::cube ctmc_mds(Q.n_rows, Q.n_cols, max_order + 1, arma::fill::zeros);

  // Initialize the counting tables.
  arma::vec factorial = factorial_table(max_order);
  arma::mat stirling_num = (mode == Mode::NSUBS_MOMENTS)
                               ? stirling_num_table(max_order)
                               : arma::mat();

  // Construct the auxiliary matrix `A`.
  // `A` is an upper block bidiagonal matrix of the form: (Q  B  0 ... 0)
  //                                                      (0  Q  B ... 0)
  //                                                      (0  0  Q ... 0)
  //                                                      (0  0  0 ... B)
  //                                                      (0  0  0 ... Q),
  // where the number of blocks is determined by `max_order`.
  arma::mat A(arma::size(Q) * (max_order + 1), arma::fill::zeros);
  A(0, 0, arma::size(Q)) = Q;

  for (int i = 0; i < max_order; ++i) {
    A(Q.n_rows * i, Q.n_cols * (i + 1), arma::size(Q)) = B;
    A(Q.n_rows * (i + 1), Q.n_cols * (i + 1), arma::size(Q)) = Q;
  }

  // Compute the matrix exponential integrals.
  // See (Van Loan, 1978) for an expression of the matrix exponential
  // `exp(A * t)`.
  arma::mat integrals_mat = arma::expmat(A * t).eval().head_rows(Q.n_rows);
  arma::cube integrals(integrals_mat.begin(), Q.n_rows, Q.n_cols,
                       max_order + 1);

  // Store the zeroth CTMC moment or rate matrix derivative (i.e. the transition
  // probability matrix).
  ctmc_mds.slice(0) = std::move(integrals.slice(0));

  for (int i = 1; i < max_order + 1; ++i) {
    // Multiply the matrix exponential integrals by the corresponding
    // factorials.
    integrals.slice(i) *= factorial(i);

    if (mode == Mode::NSUBS_MOMENTS) {
      // Mode 1: NSUBS_MOMENTS
      // Calculate the raw moments using the factorial moments and Stirling
      // numbers.
      for (int j = 1; j < i + 1; ++j) {
        ctmc_mds.slice(i) += integrals.slice(j) * stirling_num(i, j);
      }
    } else {
      // Mode 2: REWARD_MOMENTS
      // Mode 3: Q_DERIVATIVES
      ctmc_mds.slice(i) = std::move(integrals.slice(i));
    }
  }

  return ctmc_mds;
}


arma::cube ctmc_t_derivatives_aux(double t, const arma::mat& Q, int max_order) {
  // Compute CTMC branch length derivatives associated with the branch length
  // `t`.
  arma::cube ctmc_ds(Q.n_rows, Q.n_cols, max_order + 1, arma::fill::zeros);

  // Store the zeroth CTMC branch length derivative (i.e. the transition
  // probability matrix).
  ctmc_ds.slice(0) = arma::expmat(Q * t);

  for (int i = 1; i < max_order + 1; ++i) {
    ctmc_ds.slice(i) = Q * ctmc_ds.slice(i - 1);
  }

  return ctmc_ds;
}


arma::cube ctmc_moments_derivatives(double t, const arma::mat& Q,
                                    const arma::mat& B, int max_order,
                                    Mode mode) {
  if (mode != Mode::T_DERIVATIVES) {
    // Mode 1: NSUBS_MOMENTS
    // Mode 2: REWARD_MOMENTS
    // Mode 3: Q_DERIVATIVES
    return ctmc_moments_Q_derivatives_aux(t, Q, B, max_order, mode);
  } else {
    // Mode 4: T_DERIVATIVES
    return ctmc_t_derivatives_aux(t, Q, max_order);
  }
}


//' CTMC restricted moments of labeled substitution counts
//' 
//' Computes the CTMC restricted moment matrices of labeled substitution counts 
//' for all orders less than or equal to the user-specified maximum order of 
//' interest.
//' 
//' The zeroth-order restricted moment matrix is defined to be the transition 
//' probability matrix.
//' 
//' @param t A nonnegative numeric scalar representing the CTMC time interval 
//'   length.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param L An indicator matrix that labels the substitutions of interest.
//' @param max_order An integer specifying the maximum moment order of interest.
//'   
//' @return A three-dimensional array that stores the restricted moment 
//'   matrices, where the third dimension indexes the different moment orders.
//'   
//' @references Van Loan C (1978) \dQuote{Computing Integrals Involving the 
//'   Matrix Exponential}, \emph{IEEE Transactions on Automatic Control}, 
//'   23(3):395-404.
//'   
//'   Minin VN and Suchard MA (2008) \dQuote{Counting labeled transitions in 
//'   continuous-time Markov models of evolution}, \emph{Journal of Mathematical
//'   Biology}, 56(3):391-412.
//'   
//'   OUR PAPER!!
//'   
//' @seealso \code{\link{ctmc.reward.moments}}, 
//'   \code{\link{phylo.nsubs.moments}}
//'   
//' @export
// [[Rcpp::export(name = "ctmc.nsubs.moments")]]
arma::cube ctmc_nsubs_moments(double t, const Rcpp::List& subst_mod,
                              const arma::mat& L, int max_order) {
  if (t < 0.0) Rcpp::stop("'t' cannot be less than 0.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst_mod' must be an object of class 'substitution.model'.");
  if (arma::any(arma::abs(arma::vectorise(L) - 0) >= arma::datum::eps &&
                arma::abs(arma::vectorise(L) - 1) >= arma::datum::eps))
    Rcpp::stop("'L' must be an indicator matrix.");
  if (max_order < 0) Rcpp::stop("'max_order' cannot be less than 0.");

  const arma::mat& Q = subst_mod["Q"];

  if (arma::size(Q) != arma::size(L))
    Rcpp::stop("The rate matrix and 'L' must have the same dimensions.");

  return ctmc_moments_derivatives(t, Q, Q % L, max_order, Mode::NSUBS_MOMENTS);
}


//' CTMC restricted moments of total rewards
//' 
//' Computes the CTMC restricted moment matrices of total rewards for all orders
//' less than or equal to the user-specified maximum order of interest.
//' 
//' The zeroth-order restricted moment matrix is defined to be the transition 
//' probability matrix.
//' 
//' @param t A nonnegative numeric scalar representing the CTMC time interval 
//'   length.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param w A numeric vector that defines the reward per unit time (reward 
//'   rate) for each CTMC state.
//' @param max_order An integer specifying the maximum moment order of interest.
//'   
//' @return A three-dimensional array that stores the restricted moment 
//'   matrices, where the third dimension indexes the different moment orders.
//'   
//' @references Van Loan C (1978) \dQuote{Computing Integrals Involving the 
//'   Matrix Exponential}, \emph{IEEE Transactions on Automatic Control}, 
//'   23(3):395-404.
//'   
//'   Minin VN and Suchard MA (2008) \dQuote{Fast, Accurate and Simulation-Free 
//'   Stochastic Mapping}, \emph{Philosophical Transactions of the Royal Society
//'   B: Biological Sciences}, 363(1512):3985-3995.
//'   
//'   OUR PAPER!!
//'   
//' @seealso \code{\link{ctmc.nsubs.moments}}, 
//'   \code{\link{phylo.reward.moments}}
//'   
//' @export
// [[Rcpp::export(name = "ctmc.reward.moments")]]
arma::cube ctmc_reward_moments(double t, const Rcpp::List& subst_mod,
                               const arma::vec& w, int max_order) {
  if (t < 0.0) Rcpp::stop("'t' cannot be less than 0.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst_mod' must be an object of class 'substitution.model'.");
  if (max_order < 0) Rcpp::stop("'max_order' cannot be less than 0.");

  const arma::mat& Q = subst_mod["Q"];

  if (Q.n_rows != w.n_elem)
    Rcpp::stop("The rate matrix and 'w' must have compatible dimensions.");

  return ctmc_moments_derivatives(t, Q, arma::diagmat(w), max_order,
                                  Mode::REWARD_MOMENTS);
}


//' CTMC transition probability derivatives with respect to rate matrix 
//' parameters
//' 
//' Computes the CTMC transition probability derivative matrices with respect to
//' a given rate matrix parameter for all orders less than or equal to the 
//' user-specified maximum order of interest.
//' 
//' The zeroth-order transition probability derivative matrix is defined to be 
//' the transition probability matrix.
//' 
//' @param t A nonnegative numeric scalar representing the CTMC time interval 
//'   length.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param param_name The rate matrix parameter name of interest.
//' @param max_order An integer specifying the maximum derivative order of 
//'   interest.
//'   
//' @return A three-dimensional array that stores the transition probability 
//'   derivative matrices, where the third dimension indexes the different 
//'   derivative orders.
//'   
//' @references Van Loan C (1978) \dQuote{Computing Integrals Involving the 
//'   Matrix Exponential}, \emph{IEEE Transactions on Automatic Control}, 
//'   23(3):395-404.
//'   
//'   Kenney T and Gu H (2012) \dQuote{Hessian Calculation for Phylogenetic 
//'   Likelihood based on the Pruning Algorithm and its Applications}, 
//'   \emph{Statistical Applications in Genetics and Molecular Biology}, 11(4).
//'   
//'   OUR PAPER!!
//'   
//' @seealso \code{\link{ctmc.t.derivatives}}, \code{\link{phylo.Q.derivatives}}
//'   
//' @export
// [[Rcpp::export(name = "ctmc.Q.derivatives")]]
arma::cube ctmc_Q_derivatives(double t, const Rcpp::List& subst_mod,
                              const std::string& param_name, int max_order) {
  if (t < 0.0) Rcpp::stop("'t' cannot be less than 0.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst_mod' must be an object of class 'substitution.model'.");
  if (max_order < 0) Rcpp::stop("'max_order' cannot be less than 0.");
  
  const arma::mat& Q = subst_mod["Q"];
  std::string d_param_name = "d_" + param_name;
  
  if (!subst_mod.containsElementNamed(d_param_name.c_str()))
    Rcpp::stop("'param_name' is not a valid 'subst_mod' parameter name.");
  
  const arma::mat& dQ = subst_mod[d_param_name];
  
  return ctmc_moments_derivatives(t, Q, dQ, max_order, Mode::Q_DERIVATIVES);
}


//' CTMC transition probability derivatives with respect to time interval length
//' 
//' Computes the CTMC transition probability derivative matrices with respect to
//' time interval length for all orders less than or equal to the user-specified
//' maximum order of interest.
//' 
//' The zeroth-order transition probability derivative matrix is defined to be 
//' the transition probability matrix.
//' 
//' @param t A nonnegative numeric scalar representing the CTMC time interval 
//'   length.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param max_order An integer specifying the maximum derivative order of 
//'   interest.
//'   
//' @return A three-dimensional array that stores the transition probability 
//'   derivative matrices, where the third dimension indexes the different 
//'   derivative orders.
//'   
//' @references Kenney T and Gu H (2012) \dQuote{Hessian Calculation for 
//'   Phylogenetic Likelihood based on the Pruning Algorithm and its 
//'   Applications}, \emph{Statistical Applications in Genetics and Molecular 
//'   Biology}, 11(4).
//'   
//'   OUR PAPER!!
//'   
//' @seealso \code{\link{ctmc.Q.derivatives}}, \code{\link{phylo.t.derivatives}}
//'   
//' @export
// [[Rcpp::export(name = "ctmc.t.derivatives")]]
arma::cube ctmc_t_derivatives(double t, const Rcpp::List& subst_mod,
                              int max_order) {
  if (t < 0.0) Rcpp::stop("'t' cannot be less than 0.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst_mod' must be an object of class 'substitution.model'.");
  if (max_order < 0) Rcpp::stop("'max_order' cannot be less than 0.");

  const arma::mat& Q = subst_mod["Q"];

  return ctmc_moments_derivatives(t, Q, arma::mat(), max_order,
                                  Mode::T_DERIVATIVES);
}
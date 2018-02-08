#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


//   if (n < 0) Rcpp::stop("'n' cannot be less than 0.");
//   if (max_order < 0) Rcpp::stop("'max_order' cannot be less than 0.");
//   } else {
//     Rcpp::stop("'mode' must be set to either 'Mode::MOMENTS' or 'Mode::DERIVATIVES'.");
//   }


// [[Rcpp::export]]
arma::mat choose_table(int n) {
  // Initialize the combination table.
  arma::mat C(n + 1, n + 1, arma::fill::zeros);
  C.col(0).fill(1);

  // Fill the combination table.
  for (int i = 1; i < n + 1; ++i) {
    for (int j = 1; j < i + 1; ++j) {
      C(i, j) = C(i - 1, j - 1) + C(i - 1, j);
    }
  }

  return C;
}


// [[Rcpp::export]]
arma::mat stirling_num_table(int n) {
  // Initialize the Stirling number table.
  arma::mat S(n + 1, n + 1, arma::fill::zeros);
  S(0, 0) = 1;

  // Fill the Stirling number table.
  for (int i = 1; i < n + 1; ++i) {
    for (int j = 1; j < i + 1; ++j) {
      S(i, j) = S(i - 1, j - 1) + j * S(i - 1, j);
    }
  }

  return S;
}


// [[Rcpp::export]]
arma::vec factorial_table(int n) {
  // Initialize the factorial table.
  arma::vec F(n + 1, arma::fill::zeros);
  F(0) = 1;

  // Fill the factorial table.
  for (int i = 1; i < n + 1; ++i) {
    F(i) = i * F(i - 1);
  }

  return F;
}


enum class Mode { MOMENTS, Q_DERIVATIVES, T_DERIVATIVES };


arma::cube ctmc_moments_Q_derivatives_aux(double t, const arma::mat& Q,
                                          const arma::mat& B, int max_order,
                                          Mode mode) {
  // Compute either CTMC moments or CTMC rate matrix derivatives associated with
  // the branch length `t`.
  arma::cube ctmc_mds(Q.n_rows, Q.n_cols, max_order + 1, arma::fill::zeros);

  // Initialize the counting tables.
  arma::vec factorial = factorial_table(max_order);
  arma::mat stirling_num =
      (mode == Mode::MOMENTS) ? stirling_num_table(max_order) : arma::mat();

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

    if (mode == Mode::MOMENTS) {
      // Mode 1: MOMENTS
      // Calculate the raw moments using the factorial moments and Stirling
      // numbers.
      for (int j = 1; j < i + 1; ++j) {
        ctmc_mds.slice(i) += integrals.slice(j) * stirling_num(i, j);
      }
    } else {
      // Mode 2: Q_DERIVATIVES
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
  if (mode == Mode::MOMENTS || mode == Mode::Q_DERIVATIVES) {
    // Mode 1: MOMENTS
    // Mode 2: Q_DERIVATIVES
    return ctmc_moments_Q_derivatives_aux(t, Q, B, max_order, mode);
  } else {
    // Mode 3: T_DERIVATIVES
    return ctmc_t_derivatives_aux(t, Q, max_order);
  }
}







// [[Rcpp::export]]
arma::cube ctmc_moments_wrapper(double t, const arma::mat& rate_mat,
                                const arma::mat& meat_mat, int max_order) {
  return ctmc_moments_Q_derivatives_aux(t, rate_mat, meat_mat, max_order, Mode::MOMENTS);
}


// [[Rcpp::export]]
arma::cube ctmc_derivatives_wrapper(double t, const arma::mat& rate_mat,
                                    const arma::mat& meat_mat, int max_order) {
  return ctmc_moments_Q_derivatives_aux(t, rate_mat, meat_mat, max_order, Mode::Q_DERIVATIVES);
}

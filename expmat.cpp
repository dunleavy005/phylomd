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


Vector<arma::cube> ctmc_moments_Q_derivatives_aux(const arma::vec& t,
                                                  const arma::mat& Q,
                                                  const arma::mat& B,
                                                  const Vector<int>& max_orders,
                                                  Mode mode) {
  // Compute either CTMC moments or CTMC rate matrix derivatives for each entry
  // in `t`.
  Vector<arma::cube> outp;
  outp.reserve(t.n_elem);

  // Initialize the counting tables.
  int max_max_order = *std::max_element(max_orders.begin(), max_orders.end());
  arma::vec factorial = factorial_table(max_max_order);
  arma::mat stirling_num =
      (mode == Mode::MOMENTS) ? stirling_num_table(max_max_order) : arma::mat();

  for (arma::uword i = 0; i < t.n_elem; ++i) {
    // Construct the auxiliary matrix `A`.
    // `A` is an upper block bidiagonal matrix of the form: (Q  B  0 ... 0)
    //                                                      (0  Q  B ... 0)
    //                                                      (0  0  Q ... 0)
    //                                                      (0  0  0 ... B)
    //                                                      (0  0  0 ... Q),
    // where the number of blocks is determined by `max_orders[i]`.
    arma::mat A(arma::size(Q) * (max_orders[i] + 1), arma::fill::zeros);
    A(0, 0, arma::size(Q)) = Q;

    for (int j = 0; j < max_orders[i]; ++j) {
      A(Q.n_rows * j, Q.n_cols * (j + 1), arma::size(Q)) = B;
      A(Q.n_rows * (j + 1), Q.n_cols * (j + 1), arma::size(Q)) = Q;
    }

    // Compute the matrix exponential integrals.
    // See (Van Loan, 1978) for an expression of the matrix exponential
    // `exp(A * t(i))`.
    arma::mat integrals_mat = arma::expmat(A * t(i)).eval().head_rows(Q.n_rows);
    arma::cube integrals(integrals_mat.begin(), Q.n_rows, Q.n_cols,
                         max_orders[i] + 1);

    // Store the zeroth CTMC moment or rate matrix derivative (i.e. the
    // transition probability matrix).
    outp.emplace_back(arma::size(integrals), arma::fill::zeros);
    outp[i].slice(0) = std::move(integrals.slice(0));

    for (int j = 1; j < max_orders[i] + 1; ++j) {
      // Multiply the matrix exponential integrals by the corresponding
      // factorials.
      integrals.slice(j) *= factorial(j);

      if (mode == Mode::MOMENTS) {
        // Mode 1: MOMENTS
        // Calculate the raw moments using the factorial moments and Stirling
        // numbers.
        for (int k = 1; k < j + 1; ++k) {
          outp[i].slice(j) += integrals.slice(k) * stirling_num(j, k);
        }
      } else {
        // Mode 2: Q_DERIVATIVES
        outp[i].slice(j) = std::move(integrals.slice(j));
      }
    }
  }

  return outp;
}

Vector<arma::cube> ctmc_t_derivatives_aux(const arma::vec& t,
                                          const arma::mat& Q,
                                          const Vector<int>& max_orders) {
  // Compute CTMC branch length derivatives for each entry in `t`.
  Vector<arma::cube> outp;
  outp.reserve(t.n_elem);

  for (arma::uword i = 0; i < t.n_elem; ++i) {
    outp.emplace_back(Q.n_rows, Q.n_cols, max_orders[i] + 1, arma::fill::zeros);
    outp[i].slice(0) = arma::expmat(Q * t(i));

    for (int j = 1; j < max_orders[i] + 1; ++j) {
      outp[i].slice(j) = Q * outp[i].slice(j - 1);
    }
  }

  return outp;
}

Vector<arma::cube> ctmc_moments_derivatives(const arma::vec& t,
                                            const arma::mat& Q,
                                            const arma::mat& B,
                                            const Vector<int>& max_orders,
                                            Mode mode) {
  if (mode == Mode::MOMENTS || mode == Mode::Q_DERIVATIVES) {
    // Mode 1: MOMENTS
    // Mode 2: Q_DERIVATIVES
    return ctmc_moments_Q_derivatives_aux(t, Q, B, max_orders, mode);
  } else {
    // Mode 3: T_DERIVATIVES
    return ctmc_t_derivatives_aux(t, Q, max_orders);
  }
}







// [[Rcpp::export]]
Rcpp::List ctmc_moments_wrapper(const arma::vec& t, const arma::mat& rate_mat,
                                const arma::mat& meat_mat, const Vector<int>& max_orders) {
  std::vector<arma::cube> outp = ctmc_moments_Q_derivatives_aux(t, rate_mat, meat_mat, max_orders, Mode::MOMENTS);
  return Rcpp::List(outp.begin(), outp.end());
}


// [[Rcpp::export]]
Rcpp::List ctmc_derivatives_wrapper(const arma::vec& t, const arma::mat& rate_mat,
                                    const arma::mat& meat_mat, const Vector<int>& max_orders) {
  std::vector<arma::cube> outp = ctmc_moments_Q_derivatives_aux(t, rate_mat, meat_mat, max_orders, Mode::Q_DERIVATIVES);
  return Rcpp::List(outp.begin(), outp.end());
}

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]




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


enum class Mode { NSUBS_MOMENTS, REWARD_MOMENTS, Q_DERIVATIVES, T_DERIVATIVES };


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



// [[Rcpp::export]]
Rcpp::List JC69(double mu, bool scale = false) {
  if (mu < 0.0) Rcpp::stop("The rate parameter cannot be less than 0.");

  // Specify the CTMC states and create the JC69 rate matrix.
  Vector<std::string> states = {"A", "C", "G", "T"};
  arma::mat Q = {{-0.75 * mu, 0.25 * mu, 0.25 * mu, 0.25 * mu},
                 {0.25 * mu, -0.75 * mu, 0.25 * mu, 0.25 * mu},
                 {0.25 * mu, 0.25 * mu, -0.75 * mu, 0.25 * mu},
                 {0.25 * mu, 0.25 * mu, 0.25 * mu, -0.75 * mu}};

  // Differentiate the rate matrix with respect to the JC69 rate parameter.
  arma::mat d_mu = {{-0.75, 0.25, 0.25, 0.25},
                    {0.25, -0.75, 0.25, 0.25},
                    {0.25, 0.25, -0.75, 0.25},
                    {0.25, 0.25, 0.25, -0.75}};

  // Should time be specified in terms of the expected number of CTMC
  // substitutions per site?
  if (scale) {
    double scaler = arma::dot(arma::vec({0.25, 0.25, 0.25, 0.25}), -Q.diag());
    Q /= scaler;
    d_mu /= scaler;
  }

  // Construct the S3 substitution model object.
  Rcpp::List subst_mod = Rcpp::List::create(
      Rcpp::Named("states") = states, Rcpp::Named("Q") = Q,
      Rcpp::Named("pi") = Vector<double>({0.25, 0.25, 0.25, 0.25}),
      Rcpp::Named("d_mu") = d_mu);
  subst_mod.attr("class") = Vector<std::string>({"JC69", "substitution.model"});

  return subst_mod;
}


// [[Rcpp::export]]
Rcpp::List K80(double alpha, double beta, bool scale = false) {
  for (auto r : {alpha, beta}) {
    if (r < 0.0) Rcpp::stop("The rate parameters cannot be less than 0.");
  }

  // Specify the CTMC states and create the K80 rate matrix.
  Vector<std::string> states = {"A", "C", "G", "T"};
  arma::mat Q = {
      {-0.25 * alpha - 0.5 * beta, 0.25 * beta, 0.25 * alpha, 0.25 * beta},
      {0.25 * beta, -0.25 * alpha - 0.5 * beta, 0.25 * beta, 0.25 * alpha},
      {0.25 * alpha, 0.25 * beta, -0.25 * alpha - 0.5 * beta, 0.25 * beta},
      {0.25 * beta, 0.25 * alpha, 0.25 * beta, -0.25 * alpha - 0.5 * beta}};

  // Differentiate the rate matrix with respect to the K80 rate parameters.
  arma::mat d_alpha = {{-0.25, 0, 0.25, 0},
                       {0, -0.25, 0, 0.25},
                       {0.25, 0, -0.25, 0},
                       {0, 0.25, 0, -0.25}};
  arma::mat d_beta = {{-0.5, 0.25, 0, 0.25},
                      {0.25, -0.5, 0.25, 0},
                      {0, 0.25, -0.5, 0.25},
                      {0.25, 0, 0.25, -0.5}};

  // Should time be specified in terms of the expected number of CTMC
  // substitutions per site?
  if (scale) {
    double scaler = arma::dot(arma::vec({0.25, 0.25, 0.25, 0.25}), -Q.diag());
    Q /= scaler;
    d_alpha /= scaler;
    d_beta /= scaler;
  }

  // Construct the S3 substitution model object.
  Rcpp::List subst_mod = Rcpp::List::create(
      Rcpp::Named("states") = states, Rcpp::Named("Q") = Q,
      Rcpp::Named("pi") = Vector<double>({0.25, 0.25, 0.25, 0.25}),
      Rcpp::Named("d_alpha") = d_alpha, Rcpp::Named("d_beta") = d_beta);
  subst_mod.attr("class") = Vector<std::string>({"K80", "substitution.model"});

  return subst_mod;
}


// [[Rcpp::export]]
Rcpp::List F81(double mu, const arma::vec& pi, bool scale = false) {
  if (mu < 0.0) Rcpp::stop("The rate parameter cannot be less than 0.");
  if (pi.n_elem != 4 || arma::any(pi < 0.0) ||
      std::fabs(arma::sum(pi) - 1.0) >= arma::datum::eps)
    Rcpp::stop("'pi' must be an appropriate DNA frequency distribution.");

  // Specify the CTMC states and create the F81 rate matrix.
  Vector<std::string> states = {"A", "C", "G", "T"};
  arma::mat Q = {
      {-(pi(1) + pi(2) + pi(3)) * mu, pi(1) * mu, pi(2) * mu, pi(3) * mu},
      {pi(0) * mu, -(pi(0) + pi(2) + pi(3)) * mu, pi(2) * mu, pi(3) * mu},
      {pi(0) * mu, pi(1) * mu, -(pi(0) + pi(1) + pi(3)) * mu, pi(3) * mu},
      {pi(0) * mu, pi(1) * mu, pi(2) * mu, -(pi(0) + pi(1) + pi(2)) * mu}};

  // Differentiate the rate matrix with respect to the F81 rate parameter.
  arma::mat d_mu = {{-(pi(1) + pi(2) + pi(3)), pi(1), pi(2), pi(3)},
                    {pi(0), -(pi(0) + pi(2) + pi(3)), pi(2), pi(3)},
                    {pi(0), pi(1), -(pi(0) + pi(1) + pi(3)), pi(3)},
                    {pi(0), pi(1), pi(2), -(pi(0) + pi(1) + pi(2))}};

  // Should time be specified in terms of the expected number of CTMC
  // substitutions per site?
  if (scale) {
    double scaler = arma::dot(pi, -Q.diag());
    Q /= scaler;
    d_mu /= scaler;
  }

  // Construct the S3 substitution model object.
  Rcpp::List subst_mod = Rcpp::List::create(
      Rcpp::Named("states") = states, Rcpp::Named("Q") = Q,
      Rcpp::Named("pi") = arma::conv_to<Vector<double>>::from(pi),
      Rcpp::Named("d_mu") = d_mu);
  subst_mod.attr("class") = Vector<std::string>({"F81", "substitution.model"});

  return subst_mod;
}


// [[Rcpp::export]]
Rcpp::List HKY85(double alpha, double beta, const arma::vec& pi,
                 bool scale = false) {
  for (auto r : {alpha, beta}) {
    if (r < 0.0) Rcpp::stop("The rate parameters cannot be less than 0.");
  }
  if (pi.n_elem != 4 || arma::any(pi < 0.0) ||
      std::fabs(arma::sum(pi) - 1.0) >= arma::datum::eps)
    Rcpp::stop("'pi' must be an appropriate DNA frequency distribution.");

  // Specify the CTMC states and create the HKY85 rate matrix.
  Vector<std::string> states = {"A", "C", "G", "T"};
  arma::mat Q = {{-pi(2) * alpha - (pi(1) + pi(3)) * beta, pi(1) * beta,
                  pi(2) * alpha, pi(3) * beta},
                 {pi(0) * beta, -pi(3) * alpha - (pi(0) + pi(2)) * beta,
                  pi(2) * beta, pi(3) * alpha},
                 {pi(0) * alpha, pi(1) * beta,
                  -pi(0) * alpha - (pi(1) + pi(3)) * beta, pi(3) * beta},
                 {pi(0) * beta, pi(1) * alpha, pi(2) * beta,
                  -pi(1) * alpha - (pi(0) + pi(2)) * beta}};

  // Differentiate the rate matrix with respect to the HKY85 rate parameters.
  arma::mat d_alpha = {{-pi(2), 0, pi(2), 0},
                       {0, -pi(3), 0, pi(3)},
                       {pi(0), 0, -pi(0), 0},
                       {0, pi(1), 0, -pi(1)}};
  arma::mat d_beta = {{-(pi(1) + pi(3)), pi(1), 0, pi(3)},
                      {pi(0), -(pi(0) + pi(2)), pi(2), 0},
                      {0, pi(1), -(pi(1) + pi(3)), pi(3)},
                      {pi(0), 0, pi(2), -(pi(0) + pi(2))}};

  // Should time be specified in terms of the expected number of CTMC
  // substitutions per site?
  if (scale) {
    double scaler = arma::dot(pi, -Q.diag());
    Q /= scaler;
    d_alpha /= scaler;
    d_beta /= scaler;
  }

  // Construct the S3 substitution model object.
  Rcpp::List subst_mod = Rcpp::List::create(
      Rcpp::Named("states") = states, Rcpp::Named("Q") = Q,
      Rcpp::Named("pi") = arma::conv_to<Vector<double>>::from(pi),
      Rcpp::Named("d_alpha") = d_alpha, Rcpp::Named("d_beta") = d_beta);
  subst_mod.attr("class") =
      Vector<std::string>({"HKY85", "substitution.model"});

  return subst_mod;
}


// [[Rcpp::export]]
Rcpp::List GTR(double rAC, double rAG, double rAT, double rCG, double rCT,
               double rGT, const arma::vec& pi, bool scale = false) {
  for (auto r : {rAC, rAG, rAT, rCG, rCT, rGT}) {
    if (r < 0.0) Rcpp::stop("The rate parameters cannot be less than 0.");
  }
  if (pi.n_elem != 4 || arma::any(pi < 0.0) ||
      std::fabs(arma::sum(pi) - 1.0) >= arma::datum::eps)
    Rcpp::stop("'pi' must be an appropriate DNA frequency distribution.");

  // Specify the CTMC states and create the GTR rate matrix.
  Vector<std::string> states = {"A", "C", "G", "T"};
  arma::mat Q = {{-pi(1) * rAC - pi(2) * rAG - pi(3) * rAT, pi(1) * rAC,
                  pi(2) * rAG, pi(3) * rAT},
                 {pi(0) * rAC, -pi(0) * rAC - pi(2) * rCG - pi(3) * rCT,
                  pi(2) * rCG, pi(3) * rCT},
                 {pi(0) * rAG, pi(1) * rCG,
                  -pi(0) * rAG - pi(1) * rCG - pi(3) * rGT, pi(3) * rGT},
                 {pi(0) * rAT, pi(1) * rCT, pi(2) * rGT,
                  -pi(0) * rAT - pi(1) * rCT - pi(2) * rGT}};

  // Differentiate the rate matrix with respect to the GTR rate parameters.
  arma::mat d_rAC = {
      {-pi(1), pi(1), 0, 0}, {pi(0), -pi(0), 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  arma::mat d_rAG = {
      {-pi(2), 0, pi(2), 0}, {0, 0, 0, 0}, {pi(0), 0, -pi(0), 0}, {0, 0, 0, 0}};
  arma::mat d_rAT = {
      {-pi(3), 0, 0, pi(3)}, {0, 0, 0, 0}, {0, 0, 0, 0}, {pi(0), 0, 0, -pi(0)}};
  arma::mat d_rCG = {
      {0, 0, 0, 0}, {0, -pi(2), pi(2), 0}, {0, pi(1), -pi(1), 0}, {0, 0, 0, 0}};
  arma::mat d_rCT = {
      {0, 0, 0, 0}, {0, -pi(3), 0, pi(3)}, {0, 0, 0, 0}, {0, pi(1), 0, -pi(1)}};
  arma::mat d_rGT = {
      {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, -pi(3), pi(3)}, {0, 0, pi(2), -pi(2)}};

  // Should time be specified in terms of the expected number of CTMC
  // substitutions per site?
  if (scale) {
    double scaler = arma::dot(pi, -Q.diag());
    Q /= scaler;
    d_rAC /= scaler;
    d_rAG /= scaler;
    d_rAT /= scaler;
    d_rCG /= scaler;
    d_rCT /= scaler;
    d_rGT /= scaler;
  }

  // Construct the S3 substitution model object.
  Rcpp::List subst_mod = Rcpp::List::create(
      Rcpp::Named("states") = states, Rcpp::Named("Q") = Q,
      Rcpp::Named("pi") = arma::conv_to<Vector<double>>::from(pi),
      Rcpp::Named("d_rAC") = d_rAC, Rcpp::Named("d_rAG") = d_rAG,
      Rcpp::Named("d_rAT") = d_rAT, Rcpp::Named("d_rCG") = d_rCG,
      Rcpp::Named("d_rCT") = d_rCT, Rcpp::Named("d_rGT") = d_rGT);
  subst_mod.attr("class") = Vector<std::string>({"GTR", "substitution.model"});

  return subst_mod;
}






// [[Rcpp::export]]
arma::cube ctmc_nsubs_moments(double t, const Rcpp::List& subst_mod,
                              const arma::mat& L, int max_order) {
  if (t < 0.0) Rcpp::stop("'t' cannot be less than 0.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst.mod' must be an object of class 'substitution.model'.");
  if (arma::any(arma::abs(arma::vectorise(L) - 0) >= arma::datum::eps &&
                arma::abs(arma::vectorise(L) - 1) >= arma::datum::eps))
    Rcpp::stop("'L' must be an indicator matrix.");
  if (max_order < 0) Rcpp::stop("'max.order' cannot be less than 0.");

  const arma::mat& Q = subst_mod["Q"];

  if (arma::size(Q) != arma::size(L))
    Rcpp::stop("The rate matrix and 'L' must have the same dimensions.");

  return ctmc_moments_derivatives(t, Q, Q % L, max_order, Mode::NSUBS_MOMENTS);
}

// [[Rcpp::export]]
arma::cube ctmc_reward_moments(double t, const Rcpp::List& subst_mod,
                               const arma::vec& w, int max_order) {
  if (t < 0.0) Rcpp::stop("'t' cannot be less than 0.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst.mod' must be an object of class 'substitution.model'.");
  if (max_order < 0) Rcpp::stop("'max.order' cannot be less than 0.");

  const arma::mat& Q = subst_mod["Q"];

  if (Q.n_rows != w.n_elem)
    Rcpp::stop("The rate matrix and 'w' must have compatible dimensions.");

  return ctmc_moments_derivatives(t, Q, arma::diagmat(w), max_order,
                                  Mode::REWARD_MOMENTS);
}

// [[Rcpp::export]]
arma::cube ctmc_Q_derivatives(double t, const Rcpp::List& subst_mod,
                              const std::string& param_name, int max_order) {
  if (t < 0.0) Rcpp::stop("'t' cannot be less than 0.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst.mod' must be an object of class 'substitution.model'.");
  if (max_order < 0) Rcpp::stop("'max.order' cannot be less than 0.");

  std::string d_param_name = "d_" + param_name;

  if (!subst_mod.containsElementNamed(d_param_name.c_str()))
    Rcpp::stop("'param.name' is not a valid 'subst.mod' parameter name.");

  const arma::mat& Q = subst_mod["Q"];
  const arma::mat& dQ = subst_mod[d_param_name];

  return ctmc_moments_derivatives(t, Q, dQ, max_order, Mode::Q_DERIVATIVES);
}

// [[Rcpp::export]]
arma::cube ctmc_t_derivatives(double t, const Rcpp::List& subst_mod,
                              int max_order) {
  if (t < 0.0) Rcpp::stop("'t' cannot be less than 0.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst.mod' must be an object of class 'substitution.model'.");
  if (max_order < 0) Rcpp::stop("'max.order' cannot be less than 0.");

  const arma::mat& Q = subst_mod["Q"];

  return ctmc_moments_derivatives(t, Q, arma::mat(), max_order,
                                  Mode::T_DERIVATIVES);
}

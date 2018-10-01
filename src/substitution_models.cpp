#include "substitution_models.h"

#include <cmath>
#include <string>

#include "phylomd_types.h"


//
// This source file defines the CTMC substitution model functions.
//


//' The Jukes-Cantor (JC69) substitution model
//'
//' Creates a Jukes-Cantor (JC69) substitution model object.
//'
//' @param mu A nonnegative numeric scalar representing the overall substitution
//'   rate.
//' @param scale A boolean indicating whether to scale the time dimension.  If
//'   \code{TRUE}, then time is specified in terms of the expected number of
//'   CTMC substitutions per site.
//'
//' @return A list (S3 object of class \code{"substitution.model"}) with the
//'   following named entries:
//'   \describe{
//'     \item{\code{states}}{The vector of DNA states.}
//'     \item{\code{Q}}{The JC69 rate matrix.}
//'     \item{\code{pi}}{The JC69 stationary distribution.}
//'     \item{\code{d_mu}}{The JC69 rate matrix derivative with respect to
//'       \code{mu}.}
//'   }
//'
//' @references Jukes TH and Cantor CR (1969) \dQuote{Evolution of Protein
//'   Molecules}, \emph{Mammalian Protein Metabolism}, (3):21-132.
//'
//' @seealso \code{\link{K80}}, \code{\link{F81}}, \code{\link{HKY85}},
//'   \code{\link{GTR}}
//'
//' @export
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


//' The Kimura 2-parameter (K80) substitution model
//'
//' Creates a Kimura 2-parameter (K80) substitution model object.
//'
//' @param alpha A nonnegative numeric scalar representing the substitution rate
//'   of DNA transitions (\eqn{A <=> G} and \eqn{C <=> T}).
//' @param beta A nonnegative numeric scalar representing the substitution rate
//'   of DNA transversions (\eqn{A <=> C}, \eqn{A <=> T}, \eqn{C <=> G}, and
//'   \eqn{G <=> T}).
//' @param scale A boolean indicating whether to scale the time dimension.  If
//'   \code{TRUE}, then time is specified in terms of the expected number of
//'   CTMC substitutions per site.
//'
//' @return A list (S3 object of class \code{"substitution.model"}) with the
//'   following named entries:
//'   \describe{
//'     \item{\code{states}}{The vector of DNA states.}
//'     \item{\code{Q}}{The K80 rate matrix.}
//'     \item{\code{pi}}{The K80 stationary distribution.}
//'     \item{\code{d_alpha}}{The K80 rate matrix derivative with respect to
//'       \code{alpha}.}
//'     \item{\code{d_beta}}{The K80 rate matrix derivative with respect to
//'       \code{beta}.}
//'   }
//'
//' @references Kimura M (1980) \dQuote{A Simple Method for Estimating
//'   Evolutionary Rates of Base Substitutions Through Comparative Studies of
//'   Nucleotide Sequences}, \emph{Journal of Molecular Evolution},
//'   16(2):111-120.
//'
//' @seealso \code{\link{JC69}}, \code{\link{F81}}, \code{\link{HKY85}},
//'   \code{\link{GTR}}
//'
//' @export
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


//' The Felsenstein (F81) substitution model
//'
//' Creates a Felsenstein (F81) substitution model object.
//'
//' @param mu A nonnegative numeric scalar representing the overall substitution
//'   rate.
//' @param pi A numeric vector defining the F81 stationary distribution.
//' @param scale A boolean indicating whether to scale the time dimension.  If
//'   \code{TRUE}, then time is specified in terms of the expected number of
//'   CTMC substitutions per site.
//'
//' @return A list (S3 object of class \code{"substitution.model"}) with the
//'   following named entries:
//'   \describe{
//'     \item{\code{states}}{The vector of DNA states.}
//'     \item{\code{Q}}{The F81 rate matrix.}
//'     \item{\code{pi}}{The F81 stationary distribution.}
//'     \item{\code{d_mu}}{The F81 rate matrix derivative with respect to
//'       \code{mu}.}
//'   }
//'
//' @references Felsenstein J (1981) \dQuote{Evolutionary Trees from DNA
//'   Sequences: A Maximum Likelihood Approach}, \emph{Journal of Molecular
//'   Evolution}, 17(6):368-376.
//'
//' @seealso \code{\link{JC69}}, \code{\link{K80}}, \code{\link{HKY85}},
//'   \code{\link{GTR}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List F81(double mu, const arma::vec& pi, bool scale = false) {
  if (mu < 0.0) Rcpp::stop("The rate parameter cannot be less than 0.");
  if (pi.n_elem != 4 || arma::any(pi < 0.0) ||
      std::fabs(arma::sum(pi) - 1.0) >= EPS)
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


//' The Hasegawa-Kishino-Yano (HKY85) substitution model
//'
//' Creates a Hasegawa-Kishino-Yano (HKY85) substitution model object.
//'
//' @param alpha A nonnegative numeric scalar representing the substitution rate
//'   of DNA transitions (\eqn{A <=> G} and \eqn{C <=> T}).
//' @param beta A nonnegative numeric scalar representing the substitution rate
//'   of DNA transversions (\eqn{A <=> C}, \eqn{A <=> T}, \eqn{C <=> G}, and
//'   \eqn{G <=> T}).
//' @param pi A numeric vector defining the HKY85 stationary distribution.
//' @param scale A boolean indicating whether to scale the time dimension.  If
//'   \code{TRUE}, then time is specified in terms of the expected number of
//'   CTMC substitutions per site.
//'
//' @return A list (S3 object of class \code{"substitution.model"}) with the
//'   following named entries:
//'   \describe{
//'     \item{\code{states}}{The vector of DNA states.}
//'     \item{\code{Q}}{The HKY85 rate matrix.}
//'     \item{\code{pi}}{The HKY85 stationary distribution.}
//'     \item{\code{d_alpha}}{The HKY85 rate matrix derivative with respect to
//'       \code{alpha}.}
//'     \item{\code{d_beta}}{The HKY85 rate matrix derivative with respect to
//'       \code{beta}.}
//'   }
//'
//' @references Hasegawa M, Kishino H, and Yano TA (1985) \dQuote{Dating of the
//'   Human-Ape Splitting by a Molecular Clock of Mitochondrial DNA}, \emph{
//'   Journal of Molecular Evolution}, 22(2):160-174.
//'
//' @seealso \code{\link{JC69}}, \code{\link{K80}}, \code{\link{F81}},
//'   \code{\link{GTR}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List HKY85(double alpha, double beta, const arma::vec& pi,
                 bool scale = false) {
  for (auto r : {alpha, beta}) {
    if (r < 0.0) Rcpp::stop("The rate parameters cannot be less than 0.");
  }
  if (pi.n_elem != 4 || arma::any(pi < 0.0) ||
      std::fabs(arma::sum(pi) - 1.0) >= EPS)
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


//' The general time-reversible (GTR) substitution model
//'
//' Creates a general time-reversible (GTR) substitution model object.
//'
//' @param rAC,rAG,rAT,rCG,rCT,rGT Nonnegative numeric scalars representing the
//'   \eqn{A <=> C}, \eqn{A <=> G}, \eqn{A <=> T}, \eqn{C <=> G}, \eqn{C <=> T},
//'   and \eqn{G <=> T} substitution rates, respectively.
//' @param pi A numeric vector defining the GTR stationary distribution.
//' @param scale A boolean indicating whether to scale the time dimension.  If
//'   \code{TRUE}, then time is specified in terms of the expected number of
//'   CTMC substitutions per site.
//'
//' @return A list (S3 object of class \code{"substitution.model"}) with the
//'   following named entries:
//'   \describe{
//'     \item{\code{states}}{The vector of DNA states.}
//'     \item{\code{Q}}{The GTR rate matrix.}
//'     \item{\code{pi}}{The GTR stationary distribution.}
//'     \item{\code{d_rAC}}{The GTR rate matrix derivative with respect to
//'       \code{rAC}.}
//'     \item{\code{d_rAG}}{The GTR rate matrix derivative with respect to
//'       \code{rAG}.}
//'     \item{\code{d_rAT}}{The GTR rate matrix derivative with respect to
//'       \code{rAT}.}
//'     \item{\code{d_rCG}}{The GTR rate matrix derivative with respect to
//'       \code{rCG}.}
//'     \item{\code{d_rCT}}{The GTR rate matrix derivative with respect to
//'       \code{rCT}.}
//'     \item{\code{d_rGT}}{The GTR rate matrix derivative with respect to
//'       \code{rGT}.}
//'   }
//'
//' @references Tavare S (1986) \dQuote{Some Probabilistic and Statistical
//'   Problems in the Analysis of DNA Sequences}, \emph{Lectures on Mathematics
//'   in the Life Sciences}, 17(2):57-86.
//'
//' @seealso \code{\link{JC69}}, \code{\link{K80}}, \code{\link{F81}},
//'   \code{\link{HKY85}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List GTR(double rAC, double rAG, double rAT, double rCG, double rCT,
               double rGT, const arma::vec& pi, bool scale = false) {
  for (auto r : {rAC, rAG, rAT, rCG, rCT, rGT}) {
    if (r < 0.0) Rcpp::stop("The rate parameters cannot be less than 0.");
  }
  if (pi.n_elem != 4 || arma::any(pi < 0.0) ||
      std::fabs(arma::sum(pi) - 1.0) >= EPS)
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
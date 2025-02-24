% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bessel_Approximations.R
\name{Marginal_log_post_kappa_prime_prime}
\alias{Marginal_log_post_kappa_prime_prime}
\title{Marginal_log_post_kappa_prime_prime}
\usage{
Marginal_log_post_kappa_prime_prime(kappa, nu, norm_Y_bar, n)
}
\description{
Marginal_log_post_kappa_prime_prime
}
\examples{
Marginal_log_post_kappa_prime_prime(kappa=c(0, 0, 1, 2, 3), nu = 1, norm_Y_bar = .8, n = 10)
benchmark(Marginal_log_post_kappa_prime_prime(kappa=c(0, 0, 1, 2, 3), nu = 1, norm_Y_bar = .8, n = 10),
Marginal_log_post_kappa_prime_prime_old(kappa=c(0, 0, 1, 2, 3), nu = 1, norm_Y_bar = .8, n = 10),
replications = 10000 )
}

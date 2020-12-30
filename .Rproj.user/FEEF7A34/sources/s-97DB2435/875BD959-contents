# library(cmdstanr)
# library(posterior)
# library(bayesplot)
color_scheme_set("brightblue")
check_cmdstan_toolchain()
cmdstan_path()
set_cmdstan_path("OneDrive/Werk/Programming/cmdstan-2.25.0")
cmdstanr::cmdstan_version()

file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)
mod$print()
mod$exe_file()
# names correspond to the data block in the Stan program
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 2,
  refresh = 500
)
fit$summary()
# use a formula to summarize arbitrary functions, e.g. Pr(theta <= 0.5)
fit$summary("theta", pr_lt_half = ~ mean(. <= 0.235))
draws_array <- fit$draws()
str(draws_array)
draws_df <- as_draws_df(draws_array) # as_draws_matrix() for matrix
print(draws_df)
mcmc_hist(fit$draws("theta"))

# this is a draws_array object from the posterior package
str(fit$sampler_diagnostics())
tt <- fit$sampler_diagnostics()
data.frame(tt[1:10, 1, ])

diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
print(diagnostics_df)

fit$cmdstan_diagnose()

fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())
stanfit

fit_mle <- mod$optimize(data = data_list, seed = 123)

stanfit_mle <- rstan::read_stan_csv(fit_mle$output_files())



tt <- rstan::read_stan_csv("output.csv")
tt
ww <- posterior::as_draws(tt)
bayesplot::mcmc_dens(ww)
qq <- shinystan::as.shinystan(tt)

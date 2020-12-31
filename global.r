# call relevant packages
library(shiny)
library(shinystan)
library(bayesplot)
library(dplyr)
library(ggplot2)
library(shinyjs)
# sso <- get(".SHINYSTAN_OBJECT", envir = shinystan:::.sso_env) 

.inverse <- function(x) 1/x
.cloglog <- function(x) log(-log1p(-x))
.square <- function(x) x^2

# sso <- shinystan::eight_schools




### .max_t_sso

# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Trustees of Columbia University
# Copyright (C) 2018, 2019 Aki Vehtari, Paul Bürkner
# 
# Compute Geyer's initial positive sequence length for shinystan objects.
# From posterior package:
# https://github.com/jgabry/posterior/blob/4ac3054d029f5ddb6fadf1ef972344d195f530de/R/convergence.R
# Small adjustments made to get .max_t and .max_t_sso functions by Duco Veen.


# internal ----------------------------------------------------------------
#' To test is things are constant
#' 
#' @noRd
is_constant <- function(x, tol = .Machine$double.eps) {
  abs(max(x) - min(x)) < tol
}

#' Find the optimal next size for the FFT so that a minimum number of zeros
#' are padded.
fft_next_good_size <- function(N) {
  if (N <= 2)
    return(2)
  while (TRUE) {
    m <- N
    while ((m %% 2) == 0) m <- m / 2
    while ((m %% 3) == 0) m <- m / 3
    while ((m %% 5) == 0) m <- m / 5
    if (m <= 1)
      return(N)
    N <- N + 1
  }
}

#' Autocovariance estimates
#'
#' Compute autocovariance estimates for every lag for the specified
#' input sequence using a fast Fourier transform approach. The estimate
#' for lag t is scaled by N-t where N is the length of the sequence.
#' 
autocovariance <- function(x) {
  N <- length(x)
  M <- fft_next_good_size(N)
  Mt2 <- 2 * M
  yc <- x - mean(x)
  yc <- c(yc, rep.int(0, Mt2 - N))
  transform <- stats::fft(yc)
  ac <- stats::fft(Conj(transform) * transform, inverse = TRUE)
  # use "biased" estimate as recommended by Geyer (1992)
  ac <- Re(ac)[1:N] / (N^2 * 2)
  ac
}

#' Compute Geyer's initial positive sequence length.
#' From posterior package:
#' https://github.com/jgabry/posterior/blob/4ac3054d029f5ddb6fadf1ef972344d195f530de/R/convergence.R

.max_t <- function(x) {
  x <- as.matrix(x)
  nchains <- NCOL(x)
  niterations <- NROW(x)
  if (niterations < 3L || anyNA(x)) {
    return(NA)
  }
  if (any(!is.finite(x))) {
    return(NaN)
  }
  if (is_constant(x)) {
    return(NA)
  }
  acov_fun <- function(i) autocovariance(x[, i])
  acov <- lapply(seq_len(nchains), acov_fun)
  acov <- do.call(cbind, acov)
  chain_mean <- apply(x, 2, mean)
  mean_var <- mean(acov[1, ]) * niterations / (niterations - 1)
  var_plus <- mean_var * (niterations - 1) / niterations
  if (nchains > 1) {
    var_plus <- var_plus + var(chain_mean)
  }
  
  # Geyer's initial positive sequence
  rho_hat_t <- rep.int(0, niterations)
  t <- 0
  rho_hat_even <- 1
  rho_hat_t[t + 1] <- rho_hat_even
  rho_hat_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
  rho_hat_t[t + 2] <- rho_hat_odd
  while (t < NROW(acov) - 5 && !is.nan(rho_hat_even + rho_hat_odd) &&
         (rho_hat_even + rho_hat_odd > 0)) {
    t <- t + 2
    rho_hat_even = 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
    rho_hat_odd = 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat_t[t + 1] <- rho_hat_even
      rho_hat_t[t + 2] <- rho_hat_odd
    }
  }
  max_t <- t
  return(max_t)
}


#' Compute Geyer's initial positive sequence length for shinystan objects.
#' From posterior package:
#' https://github.com/jgabry/posterior/blob/4ac3054d029f5ddb6fadf1ef972344d195f530de/R/convergence.R

.max_t_sso <- function(sso, param, chains = 0) {
  
  max_t <- NULL
  for(i in 1:length(param)) {
    if(chains == 0){
      max_t[i] <- .max_t(sso@posterior_sample[, , param[i]])
    } else {
      max_t[i] <- .max_t(sso@posterior_sample[, chains, param[i]])
    }
    
  }
  return(max(max_t))
}



### mcmc_msce_hist
#' Add mcmc_mcse plot for use in shinystan. Based on Bayesplot package 
#' mcmc_neff and mcmc_rhat.
#'
# monte carlo standard error -------------------------------------------
mcmc_mcse_hist <- function(ratio, ..., binwidth = NULL, breaks = NULL) {
  check_ignored_arguments(...)
  data <- mcmc_mcse_data(ratio)
  
  ggplot(
    data,
    mapping = aes_(
      x = ~ value,
      color = ~ rating,
      fill = ~ rating)) +
    geom_histogram(
      size = .25,
      na.rm = TRUE,
      binwidth = binwidth,
      breaks = breaks) +
    scale_color_diagnostic("mcse") +
    scale_fill_diagnostic("mcse") +
    labs(x = expression(mcse/sd), y = NULL) +
    dont_expand_y_axis(c(0.005, 0)) +
    yaxis_title(FALSE) +
    yaxis_text(FALSE) +
    yaxis_ticks(FALSE) +
    bayesplot_theme_get()
}

mcmc_mcse_data <- function(ratio, ...) {
  check_ignored_arguments(...)
  ratio <- drop_NAs_and_warn(new_mcse_ratio(ratio))
  diagnostic_data_frame(ratio)
}


# internal ----------------------------------------------------------------


#' Convert numeric vector of diagnostic values to a factor
#'
diagnostic_factor <- function(x, breaks, ...) {
  UseMethod("diagnostic_factor")
}

diagnostic_factor.mcse_ratio <- function(x, breaks = c(0.05, 0.1)) {
  cut(x, breaks = c(-Inf, breaks, Inf),
      labels = c("low", "ok", "high"),
      ordered_result = FALSE)
}

diagnostic_data_frame <- function(x) {
  # x <- auto_name(sort(x))
  # quick fix getting right diagnostic factor because of class drop
  # at auto_name(sort(x))
  x <- structure(auto_name(sort(x)), class = c("mcse_ratio", "numeric"))
  stopifnot(!anyDuplicated(names(x)))
  diagnostic <- class(x)[1]
  
  d <- data.frame(
    diagnostic = diagnostic,
    parameter = factor(seq_along(x), labels = names(x)),
    value = as.numeric(x),
    rating = diagnostic_factor(x))
  
  labels <- diagnostic_color_labels[[diagnostic]]
  d$description <- as.character(labels[d$rating])
  d
}

auto_name <- function(xs) {
  if (is.null(names(xs))) {
    names(xs) <- zero_pad_int(seq_along(xs))
  }
  xs
}

# c(1, 2, 10, 20, 100) => c("001", "002", "010", "020", "100")
zero_pad_int <- function(xs) {
  formatter <- paste0("%0", max(nchar(xs)), "d")
  sprintf(formatter, xs)
}

diagnostic_points <- function(size = NULL) {
  args <- list(shape = 21, na.rm = TRUE)
  do.call("geom_point", c(args, size = size))
}


# Functions wrapping around scale_color_manual() and scale_fill_manual(), used to
# color the intervals by rhat value
scale_color_diagnostic <- function(diagnostic = c("rhat", "neff", "mcse")) {
  d <- match.arg(diagnostic)
  diagnostic_color_scale(d, aesthetic = "color")
}

scale_fill_diagnostic <- function(diagnostic = c("rhat", "neff", "mcse")) {
  d <- match.arg(diagnostic)
  diagnostic_color_scale(d, aesthetic = "fill")
}

diagnostic_color_scale <- function(diagnostic = c("rhat", "neff_ratio", "mcse_ratio"),
                                   aesthetic = c("color", "fill")) {
  diagnostic <- match.arg(diagnostic)
  aesthetic <- match.arg(aesthetic)
  dc <- diagnostic_colors(diagnostic, aesthetic)
  do.call(
    match.fun(paste0("scale_", aesthetic, "_manual")),
    list(
      name = NULL,
      drop = FALSE,
      values = dc$values,
      labels = dc$color_labels
    )
  )
}

diagnostic_colors <- function(diagnostic = c("rhat", "neff_ratio", "mcse_ratio"),
                              aesthetic = c("color", "fill")) {
  diagnostic <- match.arg(diagnostic)
  aesthetic <- match.arg(aesthetic)
  color_levels <- c("light", "mid", "dark")
  if (diagnostic == "neff_ratio") {
    color_levels <- rev(color_levels)
  }
  if (diagnostic == "mcse_ratio") {
    color_levels <- color_levels
  }
  if (aesthetic == "color") {
    color_levels <- paste0(color_levels, "_highlight")
  }
  
  color_labels <- diagnostic_color_labels[[diagnostic]]
  
  list(diagnostic = diagnostic,
       aesthetic = aesthetic,
       color_levels = color_levels,
       color_labels = color_labels,
       values = rlang::set_names(get_color(color_levels), c("low", "ok", "high")))
}

diagnostic_color_labels <- list(
  rhat = c(
    low  = expression(hat(R) <= 1.05),
    ok   = expression(hat(R) <= 1.10),
    high = expression(hat(R) > 1.10)
  ),
  neff_ratio = c(
    low  = expression(N[eff] / N <= 0.1),
    ok   = expression(N[eff] / N <= 0.5),
    high = expression(N[eff] / N > 0.5)
  ),
  mcse_ratio = c(
    low  = expression(mcse / sd <= 0.05),
    ok   = expression(mcse / sd <= 0.1),
    high = expression(mcse / sd > 0.1)
  )
)

# drop NAs from a vector and issue warning
drop_NAs_and_warn <- function(x) {
  is_NA <- is.na(x)
  if (anyNA(x)) {
    rlang::warn(paste0(
      "Dropped ", sum(is_NA), " NAs from '",
      deparse(substitute(x)), "'."
    ))
  }
  x[!is_NA]
}


#' Indexing method -- needed so that sort, etc. don't strip names.
new_mcse_ratio <- function(x) {
  # Convert a 1-d arrays to a vectors
  if (is.array(x) && length(dim(x)) == 1) {
    x <- as.vector(x)
  }
  as_mcse_ratio(validate_mcse_ratio(x))
}

validate_mcse_ratio <- function(x) {
  stopifnot(is.numeric(x), !is.list(x), !is.array(x))
  if (any(x < 0, na.rm = TRUE)) {
    rlang::abort("All mcse ratios must be positive.")
  }
  x
}

as_mcse_ratio <- function(x) {
  structure(x, class = c("mcse_ratio", "numeric"), names = names(x))
}

#' Indexing method -- needed so that sort, etc. don't strip names.
`[.mcse_ratio` <- function (x, i, j, drop = TRUE, ...) {
  as_mcse_ratio(NextMethod())
}



# Check for ignored arguments
check_ignored_arguments <- function(..., ok_args = character()) {
  dots <- list(...)
  if (length(dots)) {
    unrecognized <- if (!length(ok_args))
      names(dots) else setdiff(names(dots), ok_args)
    if (length(unrecognized)) {
      rlang::warn(paste(
        "The following arguments were unrecognized and ignored:",
        paste(unrecognized, collapse = ", ")
      ))
    }
  }
}


get_color <- function(levels) {
  levels <- full_level_name(levels)
  stopifnot(all(levels %in% scheme_level_names()))
  color_vals <- color_scheme_get()[levels]
  unlist(color_vals, use.names = FALSE)
}

full_level_name <- function(x) {
  map <- c(
    l = "light",
    lh = "light_highlight",
    m = "mid",
    mh = "mid_highlight",
    d = "dark",
    dh = "dark_highlight",
    light = "light",
    light_highlight = "light_highlight",
    mid = "mid",
    mid_highlight = "mid_highlight",
    dark = "dark",
    dark_highlight = "dark_highlight"
  )
  unname(map[x])
}

scheme_level_names <- function() {
  c("light",
    "light_highlight",
    "mid",
    "mid_highlight",
    "dark",
    "dark_highlight")
}

dont_expand_y_axis <- function(expand = c(0,0)) {
  scale_y_continuous(expand = expand)
}
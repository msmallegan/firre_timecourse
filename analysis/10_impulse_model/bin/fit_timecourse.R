#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(compiler)
library(tidyverse)
library(nls.multstart)
library(broom)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  batch <- args[1]
} else {
  stop("Too many args")
}


# Functions to fit
impulse_fun <- function(time, v_inter, v_final, t_rise, t_fall, rate) {
  (1/(1 + exp(-1*rate*(time - t_rise)))) * (v_final + (v_inter -
            v_final)*(1/(1 + exp(rate*(time - t_fall)))))
}

sigmoid_fun <- function(time, v_inter, t_rise, rate) {
  return(v_inter*(1/(1 + exp(-1*rate*(time - t_rise)))))
}

# Compile functions to speed up
cmp_impulse <- cmpfun(impulse_fun)
cmp_sigmoid <- cmpfun(sigmoid_fun)

# Wrapper functions to specify the ranges of the starts grid search
fit_sigmoid_wrapper <- function(y, x) {
  nls_multstart(y ~ cmp_sigmoid(x, v_inter, t_rise, rate),
                data = data.frame(x, y),
                iter = 500,
                start_lower = c(v_inter = 0, t_rise = 0, rate = 0.01),
                start_upper = c(v_inter = max(y), t_rise = 360, rate = 10),
                lower = c(v_inter = min(y), t_rise = 0, rate = 1e-3),
                supp_errors = "Y",
                convergence_count = 100)
}

fit_impulse_wrapper <- function(y, x) {
  nls_multstart(y ~ cmp_impulse(x, v_inter, v_final, t_rise, t_fall, rate), 
                data = data.frame(x, y),
                iter = 500,
                start_lower = c(v_inter = 0, v_final = 0, t_rise = 0, t_fall = 0, rate = 0.01),
                start_upper = c(v_inter = max(y), v_final = max(y), t_rise = 360, t_fall = 360, rate = 10),
                lower = c(v_inter = min(y), v_final = min(y), t_rise = 0, t_fall = 0, rate = 0),
                supp_errors = 'Y',
                convergence_count = 100)
}

# Compile these too -- because why not
cmp_fit_sigmoid_wrapper <- cmpfun(fit_sigmoid_wrapper)
cmp_fit_impulse_wrapper <- cmpfun(fit_impulse_wrapper)


# Import the data for this batch
rlog_counts <- read_rds(paste0("../results/rlog_counts/batch_", batch, ".rds"))

# Run the curve fits for lm, sigmoid, and imulse
cat("Beginning curve fits for", length(unique(rlog_counts$gene_id)), "genes\n")
curve_fits <- rlog_counts %>%
  dplyr::select(gene_id, gene_name, timepoint_minutes, 
                rlog_count, firre_ko, firre_induced) %>%
  unite(condition, firre_ko, firre_induced) %>%
  nest(timecourse = c(timepoint_minutes, rlog_count)) %>%
  mutate(lm_fit = map(timecourse, ~ lm(rlog_count ~ 0 + timepoint_minutes, data = .)),
         sigmoid_fit = map(timecourse, ~ cmp_fit_sigmoid_wrapper(.$rlog_count, .$timepoint_minutes)),
         impulse_fit = map(timecourse, ~ cmp_fit_impulse_wrapper(.$rlog_count, .$timepoint_minutes)))

# Save object to use in main Rmarkdown script
saveRDS(curve_fits, 
        file = paste0("../results/curve_fits/batch_", batch, ".rds"))



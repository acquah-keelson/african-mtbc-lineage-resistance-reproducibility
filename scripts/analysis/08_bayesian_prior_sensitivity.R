#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(brms); library(posterior); library(readr); library(dplyr); library(ggplot2)
})
options(mc.cores=min(4L,parallel::detectCores()))
fitdir <- file.path("outputs","bayesian","fits")
out <- file.path("outputs","bayesian","prior_sensitivity")
dir.create(file.path(out,'fits'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'tables'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'figures'),recursive=TRUE,showWarnings=FALSE)
models <- list(
  B1_primary=list(file='B1_primary_any_resistance.rds',sd=c(0.5,1.5),advanced=FALSE),
  B2_primary=list(file='B2_primary_MDR_or_higher.rds',sd=c(0.5,1.5),advanced=FALSE),
  B3_primary=list(file='B3_primary_advanced_resistance.rds',sd=c(0.5,1.5),advanced=TRUE),
  B3_geo=list(file='B3_geographic_sensitivity.rds',sd=c(0.5,1.5),advanced=TRUE)
)
make_priors <- function(sd,advanced) {
  p <- c(prior(normal(0,1),class='b'),prior(student_t(3,0,2.5),class='Intercept'))
  if (advanced) {
    for (lin in c('L1','L2','L3','L5','L6'))
      p <- c(p,set_prior(sprintf('normal(0,%s)',format(sd,scientific=FALSE)),class='b',coef=paste0('lineage',lin)))
  } else {
    p <- c(prior_string(sprintf('normal(0,%s)',format(sd,scientific=FALSE)),class='b'),
           prior(student_t(3,0,2.5),class='Intercept'),
           prior(student_t(3,0,2.5),class='sd'))
  }
  p
}
rows <- list(); diagnostics <- list()
for (nm in names(models)) {
  cfg <- models[[nm]]
  primary_path <- file.path(fitdir,cfg$file)
  if (!file.exists(primary_path)) stop('Missing original fit: ',primary_path)
  original <- readRDS(primary_path)
  prior_sd <- if(cfg$advanced) 0.75 else 1.0
  for (sd in c(prior_sd,cfg$sd)) {
    tag <- gsub('\\.','p',as.character(sd))
    dest <- file.path(out,'fits',paste0(nm,'_sd_',tag,'.rds'))
    if (identical(sd,prior_sd)) {
      fit <- original
    } else if (file.exists(dest)) {
      message('Loading saved sensitivity fit: ',dest)
      fit <- readRDS(dest)
    } else {
      message('Fitting ',nm,' prior SD=',sd)
      # update() retains the original response, covariates, data, family and model structure.
      fit <- update(original, prior=make_priors(sd,cfg$advanced),
                    chains=4, iter=if(cfg$advanced) 6000 else 4000,
                    warmup=if(cfg$advanced) 3000 else 2000,
                    seed=20261018+match(nm,names(models))*100+round(sd*10),
                    control=if(cfg$advanced) list(adapt_delta=0.995,max_treedepth=15)
                            else list(adapt_delta=0.99,max_treedepth=15),
                    backend='cmdstanr', refresh=500)
      saveRDS(fit,dest)
    }
    draws <- as_draws_df(fit)
    for (lin in c('L1','L2','L3','L5','L6')) {
      term <- paste0('b_lineage',lin)
      if (!term %in% names(draws)) stop('Missing coefficient ',term,' in ',nm)
      odds <- exp(draws[[term]])
      rows[[length(rows)+1]] <- data.frame(model=nm,prior_sd=sd,lineage=lin,
        OR_median=median(odds),CrI_low=unname(quantile(odds,.025)),
        CrI_high=unname(quantile(odds,.975)),Pr_OR_gt_1=mean(odds>1))
    }
    s <- posterior::summarise_draws(posterior::as_draws_array(fit))
    np <- brms::nuts_params(fit)
    diagnostics[[length(diagnostics)+1]] <- data.frame(model=nm,prior_sd=sd,
      max_Rhat=max(s$rhat,na.rm=TRUE),min_bulk_ESS=min(s$ess_bulk,na.rm=TRUE),
      min_tail_ESS=min(s$ess_tail,na.rm=TRUE),
      divergences=sum(np$Parameter=='divergent__' & np$Value>0))
    write_tsv(bind_rows(rows),file.path(out,'tables','prior_sensitivity_partial.tsv'))
    write_tsv(bind_rows(diagnostics),file.path(out,'tables','diagnostics_partial.tsv'))
  }
}
res <- bind_rows(rows); diag <- bind_rows(diagnostics)
write_tsv(res,file.path(out,'tables','Table_OCT2026_prior_sensitivity.tsv'))
write_tsv(diag,file.path(out,'tables','Table_OCT2026_prior_diagnostics.tsv'))

message("Prior-sensitivity modelling complete.")
message("Results written to: ", out)

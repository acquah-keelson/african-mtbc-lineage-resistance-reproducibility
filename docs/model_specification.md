# Statistical model specification

## Manuscript

**Study composition and geography shape lineage–drug resistance associations in African Mycobacterium tuberculosis complex**

This document defines the statistical models used for the lineage–drug resistance analyses. It is intended to accompany the executable scripts and frozen analysis inputs in this repository.

## Analysis populations

The principal analytical cohort contains 3,227 unmixed MTBC isolates assigned to Lineages 1–6 from 22 African countries.

Of these, 3,158 isolates have a resolved collection year and contribute to the primary collection-year-adjusted models. Sixty-nine isolates lack collection year; all 69 belong to Lineage 5.

The full cohort contains 112 contributing study accessions. The collection-year-resolved cohort contains 104 study accessions.

Lineage 4 is the reference lineage in all lineage contrast models.

## Outcomes

Three binary genomic resistance outcomes were analysed.

### Any genomic resistance

An isolate was coded as resistant if at least one resistance-associated genomic prediction was present.

### MDR or higher

This outcome included isolates classified as MDR-TB, pre-XDR-TB or XDR-TB.

### Advanced resistance

This outcome included pre-XDR-TB and XDR-TB.

There were 86 advanced-resistance isolates in the collection-year-resolved cohort.

## Collection year

Collection year was centred at the mean year of the complete-case analytical population and scaled in five-year units:

\[
Year5_i = \frac{Year_i - \bar{Year}}{5}.
\]

The mean collection year in the adjusted cohort was 2011.586.

## Frequentist models

### Crude lineage models

For isolate \(i\),

\[
\text{logit}\{P(Y_i=1)\}
=
\beta_0
+
\sum_{k \neq L4}\beta_k L_{ik}.
\]

Exponentiated lineage coefficients are reported as odds ratios relative to Lineage 4.

### Primary adjusted models: any resistance and MDR or higher

For isolate \(i\) from study \(j\),

\[
\text{logit}\{P(Y_{ij}=1)\}
=
\beta_0
+
\sum_{k \neq L4}\beta_k L_{ik}
+
\beta_y Year5_i
+
u_j,
\]

where

\[
u_j \sim N(0,\sigma^2_{\text{study}}).
\]

These models contain a random intercept for contributing study and were fitted to the 3,158 isolates with resolved collection year.

The reported lineage odds ratios are conditional odds ratios obtained by exponentiating the fixed lineage coefficients.

### Advanced-resistance model

Advanced resistance was sparse and could not reliably support the study-level random intercept used for the more common outcomes.

The primary model was therefore:

\[
\text{logit}\{P(Y_i=1)\}
=
\beta_0
+
\sum_{k \neq L4}\beta_k L_{ik}
+
\beta_y Year5_i.
\]

Bias-reduced logistic regression was used.

Lineage 5 contained no advanced-resistance events. Its frequentist contrast was therefore treated as not estimable rather than interpreted from the unstable fitted coefficient.

### Geographic sensitivity model for advanced resistance

African subregion was added as a fixed effect:

\[
\text{logit}\{P(Y_i=1)\}
=
\beta_0
+
\sum_{k \neq L4}\beta_k L_{ik}
+
\beta_y Year5_i
+
\sum_r \gamma_r Region_{ir}.
\]

This model assesses the extent to which the lineage association changes after accounting for broad geographic composition.

Study and subregion effects cannot be cleanly separated in this dataset because 101 of the 104 studies represented in the adjusted cohort contributed isolates from only one African subregion.

## Full-cohort sensitivity analysis

For any resistance and MDR-or-higher resistance, models were additionally fitted to all 3,227 isolates.

Because collection year is unavailable for 69 isolates, these models omit the year term:

\[
\text{logit}\{P(Y_{ij}=1)\}
=
\beta_0
+
\sum_{k \neq L4}\beta_k L_{ik}
+
u_j.
\]

Differences between complete-case and full-cohort estimates therefore reflect both restoration of isolates with missing year and removal of collection year from the model. They must not be interpreted as a pure missing-data effect.

## Bayesian models

Bayesian models used the same linear predictors as the corresponding frequentist models.

For the primary any-resistance and MDR-or-higher models, lineage and other population-level coefficients used weakly informative Normal(0, 1) priors.

The intercept used:

\[
t_3(0,2.5).
\]

The study random-effect standard deviation used a half-Student-t prior corresponding to:

\[
t_3(0,2.5)
\]

restricted to positive values.

For sparse advanced-resistance models, lineage coefficients used a more regularising Normal(0, 0.75) prior.

Models used four Markov chains.

The any-resistance and MDR-or-higher models used 4,000 iterations per chain with 2,000 warm-up iterations.

Advanced-resistance models used 6,000 iterations per chain with 3,000 warm-up iterations.

Sampling controls used `adapt_delta = 0.99` for the any-resistance and MDR-or-higher models and `adapt_delta = 0.995` for advanced-resistance models, with maximum tree depth 15.

Posterior lineage effects are reported as posterior median odds ratios, 95% credible intervals, and posterior probabilities that the odds ratio exceeds one.

## Bayesian prior sensitivity

Prior sensitivity was evaluated by varying prior standard deviations.

For the any-resistance and MDR-or-higher models, population-level coefficient priors were evaluated at Normal(0, 0.5), Normal(0, 1.0) and Normal(0, 1.5).

For the sparse advanced-resistance models, the lineage-specific coefficient prior was evaluated at Normal(0, 0.5), Normal(0, 0.75) and Normal(0, 1.5), while other fixed-effect priors retained their original specification.

This distinction reflects the implementation in the original validated analysis and is reproduced exactly in `08_bayesian_prior_sensitivity.R`.

The principal purpose of this analysis was to distinguish directionally stable associations from estimates strongly influenced by regularisation.

## Interpretation constraints

The assembled dataset was purposively derived from public genomic archives and is not a probability sample of tuberculosis cases or MTBC infections in Africa.

Accordingly:

- lineage frequencies must not be interpreted as population prevalence estimates;
- lineage odds ratios describe associations within the assembled genomic collection;
- adjustment for contributing study reduces but does not eliminate ascertainment bias;
- study and geography are strongly entangled;
- treatment history and several patient-level determinants of resistance were unavailable;
- an observed lineage association must not be interpreted as an intrinsic biological propensity to develop resistance.

The Lineage 5 MDR-or-higher estimate is explicitly considered prior-sensitive and sensitive to collection-year completeness.

The Lineage 5 advanced-resistance posterior is prior-driven because no Lineage 5 isolate had the outcome and must not be interpreted as evidence of protection.

The Lineage 2 advanced-resistance association remains positive across prior specifications but attenuates substantially after adjustment for African subregion. Its direction is therefore robust to prior specification, whereas its magnitude is sensitive to geographic composition.

## Software

Frequentist analyses were performed in R using `lme4` and `brglm2`.

Bayesian analyses were performed with `brms`, `rstan`/Stan infrastructure and `cmdstanr` for prior-sensitivity refitting.

Exact package and session information are provided under `environment/`.

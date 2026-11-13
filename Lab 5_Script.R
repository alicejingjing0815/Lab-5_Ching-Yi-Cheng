##Lab5_Count Regression: Workflow
##Poisson, negative binomial, and probing interactions with biochemist publications"

# Part 1: Configure the project and examine the data
### 1.Load packages and prepare the setup
install.packages("tidyverse")
install.packages("glmmTMB")
install.packages("DHARMa")
install.packages("ggeffects")
install.packages("emmeans")
install.packages("performance")
install.packages("broom.mixed")
install.packages("gt")
install.packages("pscl")
install.packages("AER")
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(emmeans)
library(performance)
library(broom.mixed)
library(gt)
library(pscl)
library(AER)

data("NMES1988", package = "AER")
dim(NMES1988)
nmes <- NMES1988 |>
  select(visits, chronic, health, gender, school, insurance) |>
  mutate(
    health = fct_relevel(health, "average"),
    gender = fct_relevel(gender, "male"),
    insurance = fct_relevel(insurance, "no")
  )

glimpse(nmes)
head(nmes)
colSums(is.na(nmes))
summary(nmes)

table(nmes$health)
table(nmes$gender)
table(nmes$insurance)

table(nmes$visits)
head(sort(nmes$visits, decreasing = TRUE), 10)

### 2.Explore the data
#| label: visits-summary
#| echo: false
visits_summary <- nmes |>
  summarise(
    n = sum(!is.na(visits)),
    mean_visits = mean(visits, na.rm = TRUE),
    variance_visits = var(visits, na.rm = TRUE),
    proportion_zero = mean(visits == 0, na.rm = TRUE),
    max_visits = max(visits, na.rm = TRUE)
  )

visits_summary |>
  gt::gt() |>
  gt::tab_header(title = "Descriptive Statistics for Physician Office Visits") |>
  gt::cols_label(
    n = "n",
    mean_visits = "Mean",
    variance_visits = "Variance",
    proportion_zero = "Proportion with zero visits",
    max_visits = "Maximum"
  ) |>
  gt::fmt_integer(columns = c(n, max_visits)) |>
  gt::fmt_number(columns = c(mean_visits, variance_visits), decimals = 2) |>
  gt::fmt_percent(columns = proportion_zero, decimals = 1)

###Bat chart
#| label: fig-visits-distribution
#| fig-cap: "Distribution of physician office visits."
#| echo: false

ggplot(nmes, aes(x = visits)) +
  geom_bar(fill = "#386B8C", width = 0.85, na.rm = TRUE) +
  labs(
    x = "Number of physician office visits",
    y = "Number of participants"
  ) +
  theme_minimal(base_size = 12)

# Part 2: Specify competing count models
## 2.Specify competing count models
## 2.1.Why not a linear model?
#| label: fig-gaussian-diagnostics
#| fig-cap: "DHARMa residual diagnostics for the Gaussian model."
#| echo: false

m_gaussian <- glmmTMB::glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  data = nmes,
  family = gaussian(link = "identity")
)

set.seed(9041)

res_gaussian <- DHARMa::simulateResiduals(
  fittedModel = m_gaussian,
  n = 500,
  seed = 9041,
  plot = FALSE
)

plot(res_gaussian)

## 2.2. Specify three nested Poisson models
#| label: poisson-models
#| echo: false

m0_poisson <- glmmTMB::glmmTMB(
  visits ~ 1,
  data = nmes,
  family = poisson(link = "log")
)

m1_poisson <- glmmTMB::glmmTMB(
  visits ~ chronic,
  data = nmes,
  family = poisson(link = "log")
)

m2_poisson <- glmmTMB::glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  data = nmes,
  family = poisson(link = "log")
)

#| label: poisson-model-comparison
#| echo: false

poisson_lrt <- anova(m0_poisson, m1_poisson, m2_poisson)

comparison_table <- tibble::tibble(
  model = c("M0 (Null)", "M1 (Chronic)", "M2 (Full)"),
  aic = poisson_lrt[["AIC"]],
  bic = poisson_lrt[["BIC"]],
  chi_sq = poisson_lrt[["Chisq"]],
  df = poisson_lrt[["Chi Df"]],
  p_value = poisson_lrt[["Pr(>Chisq)"]]
)

p_label <- function(p) {
  ifelse(p < .001, "< .001", paste0("= ", sprintf("%.3f", p)))
}

lrt_m1 <- comparison_table[2, ]
lrt_m2 <- comparison_table[3, ]

comparison_table |>
  gt::gt() |>
  gt::tab_header(title = "Comparison of Nested Poisson Models") |>
  gt::cols_label(
    model = "Model",
    aic = "AIC",
    bic = "BIC",
    chi_sq = "χ²",
    df = "df",
    p_value = "p"
  ) |>
  gt::fmt_number(columns = c(aic, bic, chi_sq), decimals = 2) |>
  gt::fmt(columns = p_value, fns = function(p) {
    ifelse(p < .001, "< .001", sprintf("%.3f", p))
  }) |>
  gt::sub_missing(
    columns = c(chi_sq, df, p_value),
    missing_text = "—"
  ) |>
  gt::tab_source_note(
    "Each likelihood-ratio test compares the model with the preceding row."
  )

# Part 3. Evaluate the Poisson model
## 3.1. Simulate residuals
#| label: fig-poisson-diagnostics
#| fig-cap: "DHARMa residual diagnostics for the full Poisson model (M2)."
#| fig-width: 10
#| fig-height: 5
#| echo: false

res_poisson <- DHARMa::simulateResiduals(
  fittedModel = m2_poisson,
  n = 500,
  seed = 9041,
  plot = FALSE
)

plot(res_poisson)

## 3.2. Test for overdispersion and excess zeros
#| label: poisson-diagnostic-tests
#| echo: false

disp_poisson <- DHARMa::testDispersion(res_poisson, plot = FALSE)
zeros_poisson <- DHARMa::testZeroInflation(res_poisson, plot = FALSE)
disp_poisson
zeros_poisson

## 3.3. Interpret the diagnostics

# Part 4: Fit and evaluate a negative binomial model
## 4.1. Estimate the NB2 model
#| label: fit-nb2
#| echo: false

m_nb2 <- glmmTMB::glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  data = nmes,
  family = glmmTMB::nbinom2(link = "log")
)

## 4.2. Compare the Poisson and NB2 models
#| label: compare-poisson-nb2
#| echo: false

nb2_lrt <- anova(m2_poisson, m_nb2)
nb2_lrt

## 4.3. Re-check the NB2 model with DHARMa
#| label: fig-nb2-diagnostics
#| fig-cap: "DHARMa residual diagnostics for the full negative binomial model."
#| fig-width: 10
#| fig-height: 5
#| echo: false

res_nb2 <- DHARMa::simulateResiduals(
  fittedModel = m_nb2,
  n = 500,
  seed = 9041,
  plot = FALSE
)

plot(res_nb2)

disp_nb2 <- DHARMa::testDispersion(res_nb2, plot = FALSE)
zeros_nb2 <- DHARMa::testZeroInflation(res_nb2, plot = FALSE)


## 4.4. Visualize predicted visit counts
#| label: fig-nb2-predicted-visits
#| fig-cap: "Expected physician office visits by number of chronic conditions from the NB2 model."
#| echo: false

pred_nb2 <- ggeffects::predict_response(
  m_nb2,
  terms = "chronic [all]",
  margin = "mean_reference"
)

plot(pred_nb2) +
  ggplot2::labs(
    x = "Number of chronic conditions",
    y = "Expected number of physician office visits"
  )

# Part 5. Decide whether the zeros need their own model
zeros_nb2

## 5.1. Make the decision
m_zinb <- update(
  m_nb2,
  ziformula = ~ chronic + health + gender + school + insurance
)

m_hurdle <- update(
  m_nb2,
  family = glmmTMB::truncated_nbinom2(link = "log"),
  ziformula = ~ chronic + health + gender + school + insurance
)

lrt_zinb <- anova(m_nb2, m_zinb)
lrt_zinb

AIC(m_nb2, m_zinb, m_hurdle)
BIC(m_nb2, m_zinb, m_hurdle)

## 5.2.Justify the choice
#| label: fig-hurdle-diagnostics
#| fig-cap: "DHARMa residual diagnostics for the hurdle negative binomial model."
#| fig-width: 10
#| fig-height: 5
#| echo: false

res_hurdle <- DHARMa::simulateResiduals(
  fittedModel = m_hurdle,
  n = 500,
  seed = 9041,
  plot = FALSE
)

plot(res_hurdle)

disp_hurdle <- DHARMa::testDispersion(res_hurdle, plot = FALSE)
zeros_hurdle <- DHARMa::testZeroInflation(res_hurdle, plot = FALSE)

# Part 6. Test and probe an interaction
## 6.1. State the prediction

## 6.2. Fit the interaction model
#| label: fit-hurdle-interaction
#| echo: false

m_hurdle_int <- update(
  m_hurdle,
  formula. = visits ~ chronic * insurance + health + gender + school,
  ziformula = ~ chronic + health + gender + school + insurance
)

## 6.3. Test the interaction
#| label: hurdle-interaction-lrt
#| echo: false

lrt_int <- anova(m_hurdle, m_hurdle_int)
lrt_int

## 6.4. Plot predicted counts
#| label: fig-hurdle-insurance-interaction
#| fig-cap: "Predictions from the hurdle count component by chronic conditions and private insurance. Shaded bands indicate 95% confidence intervals."
#| echo: false

pred_int <- ggeffects::predict_response(
  m_hurdle_int,
  terms = c("chronic [all]", "insurance"),
  type = "fixed",
  margin = "mean_reference"
)

pred_int_df <- as.data.frame(pred_int)

ggplot2::ggplot(
  pred_int_df,
  ggplot2::aes(x = x, y = predicted, colour = group, fill = group)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.16,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_colour_manual(values = c(no = "#376B8C", yes = "#A84B29")) +
  ggplot2::scale_fill_manual(values = c(no = "#376B8C", yes = "#A84B29")) +
  ggplot2::labs(
    x = "Number of chronic conditions",
    y = "Predicted count-component mean (before truncation)",
    colour = "Private insurance",
    fill = "Private insurance"
  ) +
  ggplot2::theme_minimal(base_size = 12)

## 6.5. Probe the interaction
#| label: hurdle-simple-slopes
#| echo: false

chronic_slopes <- emmeans::emtrends(
  m_hurdle_int,
  specs = ~ insurance,
  var = "chronic",
  component = "cond"
)

summary(chronic_slopes, infer = c(TRUE, TRUE))
summary(pairs(chronic_slopes), infer = c(TRUE, TRUE))

summary(m_hurdle_int)$coefficients$cond["chronic:insuranceyes", ]


#| label: hurdle-simple-effects
#| echo: false

chronic_mean <- mean(nmes$chronic, na.rm = TRUE)
chronic_sd <- sd(nmes$chronic, na.rm = TRUE)
chronic_values <- chronic_mean + c(-1, 0, 1) * chronic_sd

insurance_emm <- emmeans::emmeans(
  m_hurdle_int,
  specs = ~ insurance | chronic,
  at = list(chronic = chronic_values),
  component = "cond"
)

summary(insurance_emm, type = "response")

insurance_contrasts <- pairs(insurance_emm)
summary(insurance_contrasts, infer = c(TRUE, TRUE))


## 6.6. Report the simple effects
#| label: tbl-insurance-simple-effects
#| tbl-cap: "Simple effects of private insurance at selected values of chronic conditions."
#| echo: false

simple_effects <- as.data.frame(summary(insurance_contrasts))

probe_labels <- sprintf(
  "%s (%.3f)",
  c("−1 SD", "Mean", "+1 SD"),
  chronic_values
)

simple_effects |>
  dplyr::mutate(
    probe = probe_labels[match(round(chronic, 6), round(chronic_values, 6))],
    contrast = "No insurance − Private insurance"
  ) |>
  dplyr::arrange(chronic) |>
  dplyr::select(probe, contrast, estimate, SE, z.ratio, p.value) |>
  gt::gt() |>
  gt::cols_label(
    probe = "Chronic conditions",
    contrast = "Contrast",
    estimate = "Log contrast",
    SE = "SE",
    z.ratio = gt::md("*z*"),
    p.value = gt::md("*p*")
  ) |>
  gt::fmt_number(columns = c(estimate, SE, z.ratio), decimals = 3) |>
  gt::fmt(
    columns = p.value,
    fns = function(p) sub("^= ", "", p_label(p))
  ) |>
  gt::tab_source_note(
    source_note = "Tests were computed on the log scale of the hurdle model's count component."
  )


## 6.7. Interpret the interaction
# Part 7. Results and discussion
## 7.1. Results
#| label: results-values
#| echo: false

format_b <- function(b, lower, upper, p) {
  sprintf(
    "(*b* = %.3f, 95%% CI [%.3f, %.3f], *p* %s)",
    b, lower, upper, p_label(p)
  )
}

format_lrt <- function(result) {
  sprintf(
    "χ²(%d) = %.2f, *p* %s",
    result[["Chi Df"]][2],
    result[["Chisq"]][2],
    p_label(result[["Pr(>Chisq)"]][2])
  )
}

format_diagnostic <- function(result) {
  sprintf("%.2f (*p* %s)", result$statistic, p_label(result$p.value))
}

predictor_labels <- c(
  chronic = "each additional chronic condition",
  healthpoor = "poor versus average self-rated health",
  healthexcellent = "excellent versus average self-rated health",
  genderfemale = "female versus male gender",
  school = "each additional year of education",
  insuranceyes = "having versus lacking private insurance"
)

hurdle_coefficients <- summary(m_hurdle)$coefficients

describe_component <- function(coef_table, outcome) {
  terms <- names(predictor_labels)
  estimates <- coef_table[terms, "Estimate"]
  ses <- coef_table[terms, "Std. Error"]
  p_values <- coef_table[terms, "Pr(>|z|)"]
  critical_z <- qnorm(.975)
  
  sentences <- paste0(
    predictor_labels[terms], " was associated with ",
    ifelse(estimates >= 0, "a higher ", "a lower "), outcome, " ",
    format_b(
      estimates,
      estimates - critical_z * ses,
      estimates + critical_z * ses,
      p_values
    )
  )
  
  nonsignificant <- predictor_labels[terms[p_values >= .05]]
  significance_sentence <- if (length(nonsignificant) == 0) {
    "All six contrasts had p-values below .05."
  } else {
    paste0(
      "The following associations were not statistically significant: ",
      paste(nonsignificant, collapse = ", "), "."
    )
  }
  
  paste0(
    paste(sentences, collapse = "; "), ". ", significance_sentence
  )
}

count_results <- describe_component(
  hurdle_coefficients$cond,
  "count-component mean parameter"
)

zero_results <- describe_component(
  hurdle_coefficients$zi,
  "odds of having no physician office visits"
)

slope_results <- as.data.frame(
  summary(chronic_slopes, infer = c(TRUE, TRUE))
)
slope_results$report <- with(
  slope_results,
  format_b(chronic.trend, asymp.LCL, asymp.UCL, p.value)
)
slope_no <- slope_results$report[slope_results$insurance == "no"]
slope_yes <- slope_results$report[slope_results$insurance == "yes"]

effects_ordered <- simple_effects[order(simple_effects$chronic), ]
effect_results <- with(
  effects_ordered,
  sprintf("%.3f (*p* %s)", estimate, p_label(p.value))
)

## 7.2. Discussion















































































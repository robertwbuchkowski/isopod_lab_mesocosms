install.packages(c("tidyverse", "readxl", "glmmTMB", "DHARMa", "emmeans"))
library(tidyverse)
library(readxl)
library(glmmTMB)
library(DHARMa)
library(emmeans)

library(readxl)
feeding <- read_excel("Strip Assessments.xlsx")
View(feeding)

feeding <- feeding %>%
  mutate(
    Species = factor(Species),
    Env = factor(`Strip Env(M/D)`),
    Diet = factor(`Strip Diet (HN/LN)`),
    Arena = factor(`Arena #`),
    Trial = factor(`Trial #`),
    Area = `Area Consumed(mm2)`
  )

feeding <- feeding %>%
  mutate(
    Fed = ifelse(Area > 0, 1, 0)
  )

feeding <- feeding %>%
  mutate(
    Species = factor(Species),
    Species = relevel(Species, ref = "O. asellus")  # set reference if desired
  )


m_area_all <- glmmTMB(
  Area ~ Species*Env + Env * Diet + (1 | Arena),
  dispformula = ~Env + Species,
  data = feeding
)

summary(m_area_all)
plot(simulateResiduals(m_area_all))
testCategorical(simulateResiduals(m_area_all), catPred = feeding$Env)
testCategorical(simulateResiduals(m_area_all), catPred = feeding$Species)
testCategorical(simulateResiduals(m_area_all), catPred = feeding$Diet)
testCategorical(simulateResiduals(m_area_all), catPred = feeding$Trial)
testCategorical(simulateResiduals(m_area_all), catPred = feeding$Arena)














table(feeding$Fed)

m_fed <- glmmTMB(
  Fed ~ Species * Env * Diet + (1 | Trial) + (1 | Arena),
  family = binomial,
  data = feeding
)

summary(m_fed)

feeding_pos <- feeding %>% filter(Area > 0)

m_area <- glmmTMB(
  Area ~ Species * Env * Diet + (1 | Trial) + (1 | Arena),
  family = Gamma(link = "log"),
  data = feeding_pos
)


summary(m_area)

m_area_all <- glmmTMB(
  Area ~ Species * Env * Diet + (1 | Trial) + (1 | Arena),
  ziformula = ~ Env + Species,
  dispformula = ~Env + Species,
  data = feeding
)

summary(m_area_all)

emmeans(m_area_all,Area ~ Species, type = "response")

library(DHARMa)

plot(simulateResiduals(m_fed))
plot(simulateResiduals(m_area))
plot(simulateResiduals(m_area_all))
testCategorical(simulateResiduals(m_fed), catPred = feeding$Env)
testCategorical(simulateResiduals(m_fed), catPred = feeding$Diet)

emmeans(m_fed, ~ Species * Env * Diet, type = "response")
emmeans(m_area, ~ Species * Env * Diet, type = "response")

fed_emm <- emmeans(
  m_fed,
  ~ Species * Env * Diet,
  type = "response"
) %>% 
  as.data.frame()

ggplot(fed_emm, aes(x = Env, y = prob, color = Diet)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.2,
    position = position_dodge(width = 0.3)
  ) +
  facet_wrap(~ Species) +
  labs(
    y = "Probability of feeding",
    x = "Environment"
  ) +
  ylim(0, 1) +
  theme_classic()

head(fed_emm)
str(fed_emm)

summary(fed_emm$asymp.LCL)
summary(fed_emm$asymp.UCL)

table(fed_emm$Species, fed_emm$Diet)

area_emm <- emmeans(
  m_area,
  ~ Species * Env * Diet,
  type = "response"
) %>% 
  as.data.frame()

ggplot() +
  geom_jitter(
    data = feeding_pos,
    aes(x = Env, y = Area, color = Diet),
    width = 0.15,
    alpha = 0.6
  ) +
  geom_point(
    data = area_emm,
    aes(x = Env, y = response, color = Diet),
    size = 4,
    position = position_dodge(width = 0.3)
  ) +
  geom_errorbar(
    data = area_emm,
    aes(
      x = Env,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      color = Diet
    ),
    width = 0.2,
    position = position_dodge(width = 0.3)
  ) +
  facet_wrap(~ Species) +
  labs(
    y = expression("Area consumed (mm"^2*")"),
    x = "Environment"
  ) +
  theme_classic()


mdc_pairwise <- function(means, sds, ns,
                         alpha = 0.05,
                         power = 0.80,
                         var_method = c("welch", "pooled"),
                         mcomp = c("none", "bonferroni")) {
  # Arguments:
  # means: numeric vector of treatment means (length K)
  # sds:  numeric vector of treatment SDs  (length K)
  # ns:   numeric vector of treatment sample sizes (length K)
  # alpha: type I error (two-sided)
  # power: desired power (1 - beta)
  # var_method: "welch" (unequal variances; default) or "pooled" (equal variances)
  # mcomp: multiple-comparison adjustment for alpha across all pairwise tests:
  #       "none" (default) or "bonferroni"
  #
  # Returns: list with
  # - mdc: K x K symmetric matrix of minimum detectable absolute differences
  # - alpha_used: the per-comparison alpha actually used
  # - details: list with var_method and mcomp
  
  # Basic checks
  if (!(is.numeric(means) && is.numeric(sds) && is.numeric(ns))) {
    stop("means, sds, and ns must be numeric vectors.")
  }
  if (!(length(means) == length(sds) && length(sds) == length(ns))) {
    stop("means, sds, and ns must be the same length.")
  }
  if (any(ns <= 0) || any(sds < 0)) {
    stop("ns must be > 0 and sds must be >= 0.")
  }
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0,1).")
  if (power <= 0 || power >= 1) stop("power must be in (0,1).")
  
  var_method <- match.arg(var_method)
  mcomp <- match.arg(mcomp)
  
  K <- length(means)
  
  # Number of pairwise tests
  n_pairs <- K * (K - 1) / 2
  
  # Multiple-comparison-adjusted alpha (Bonferroni optional)
  if (mcomp == "bonferroni" && n_pairs > 0) {
    alpha_star <- alpha / n_pairs
  } else {
    alpha_star <- alpha
  }
  
  # Critical values (two-sided alpha*, and power)
  z_alpha <- qnorm(1 - alpha_star / 2)
  z_power <- qnorm(power)
  
  # Prepare output matrix
  mdc_mat <- matrix(0, nrow = K, ncol = K)
  rownames(mdc_mat) <- colnames(mdc_mat) <- paste0("Trt", seq_len(K))
  
  if (var_method == "pooled") {
    # Common (pooled) SD estimate:
    # s_p^2 = sum((n_i - 1) * s_i^2) / sum(n_i - 1)
    df_sum <- sum(ns - 1)
    if (df_sum <= 0) stop("Not enough degrees of freedom to compute pooled SD.")
    sp2 <- sum((ns - 1) * (sds^2)) / df_sum
    sp <- sqrt(sp2)
  }
  
  # Fill pairwise MDCs
  for (i in 1:K) {
    for (j in 1:K) {
      if (i == j) {
        mdc_mat[i, j] <- 0
      } else {
        if (var_method == "welch") {
          # Unequal variances: use SE_ij = sqrt(sd_i^2/n_i + sd_j^2/n_j)
          SE_ij <- sqrt((sds[i]^2) / ns[i] + (sds[j]^2) / ns[j])
        } else {
          # Equal variances: SE_ij = sp * sqrt(1/n_i + 1/n_j)
          SE_ij <- sp * sqrt(1 / ns[i] + 1 / ns[j])
        }
        mdc_mat[i, j] <- (z_alpha + z_power) * SE_ij
      }
    }
  }
  
  # Make symmetric (numerical stability)
  mdc_mat <- (mdc_mat + t(mdc_mat)) / 2
  
  return(list(
    mdc = mdc_mat,
    alpha_used = alpha_star,
    details = list(var_method = var_method,
                   mcomp = mcomp,
                   alpha_input = alpha,
                   power = power)
  ))
}


feeding_2 = feeding %>%
  group_by(Species, Diet, Env) %>%
  summarize(
    Mean = mean(Area),
    SD   = sd(Area),
    N    = n(),
    .groups = "drop"
  )

mdc_pairwise(
  means = feeding_2$Mean,
  sds   = feeding_2$SD,
  ns    = feeding_2$N
)

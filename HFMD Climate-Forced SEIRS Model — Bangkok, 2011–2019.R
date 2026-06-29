# ==============================================================================
# HFMD Climate-Forced SEIRS Model — Bangkok, 2011–2019
# ==============================================================================

library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(bbmle)


raw_data <- read.csv("Data-of-hfmd-cases-and-terraclimate-2011-to-2022raw.csv",
                     stringsAsFactors = FALSE)

bkk_data <- raw_data %>%
  filter(province == "Bangkok", year <= 2019) %>%
  mutate(
    date       = as.Date(paste(year, month, "01", sep = "-")),
    t_day      = as.numeric(date - min(date)),
    temp       = mean_tmax,
    RH         = mean_vap,
    rain       = mean_ppt,
    month_name = factor(month.abb[month], levels = month.abb)
  ) %>%
  arrange(date)

cat("Data loaded successfully: Bangkok", nrow(bkk_data), "months",
    min(bkk_data$year), "–", max(bkk_data$year), "\n")

# Global Parameters
sigma_val <- 1 / 3.5
gamma_val <- 1 / 7
omega_val <- 1 / (365 * 0.5)
N_pop     <- round(mean(bkk_data$population) * 0.07)
t_days    <- bkk_data$t_day

# Weather covariate lagging function
lag_fn <- function(x, k) {
  if (k == 0) return(x)
  c(rep(NA, k), head(x, -k))
}

# Apply lags (Temp = 3 months, RH = 2 months, Rain = 0 months)
temp_lag <- lag_fn(bkk_data$temp, 3)
rh_lag   <- lag_fn(bkk_data$RH,   2)
rain_lag <- lag_fn(bkk_data$rain, 0)

temp_lag[is.na(temp_lag)] <- mean(bkk_data$temp, na.rm = TRUE)
rh_lag[is.na(rh_lag)]     <- mean(bkk_data$RH,   na.rm = TRUE)
rain_lag[is.na(rain_lag)] <- mean(bkk_data$rain, na.rm = TRUE)

# Standardize covariates (z-scores)
temp_z <- as.numeric(scale(temp_lag))
rh_z   <- as.numeric(scale(rh_lag))
rain_z <- as.numeric(scale(rain_lag))


# ==============================================================================
# 1. SEASONAL PATTERNS ANALYSIS
# ==============================================================================

pal_seas <- setNames(
  colorRampPalette(c("#313695","#74add1","#fee090","#d73027"))(length(unique(bkk_data$year))),
  as.character(sort(unique(bkk_data$year)))
)

p_cases <- ggplot(bkk_data, aes(month_name, cases, group = factor(year), colour = factor(year))) +
  geom_line(linewidth = 0.8) + scale_colour_manual(values = pal_seas, name = "Year") +
  labs(title = "HFMD Cases", x = NULL, y = "Cases") + theme_bw(base_size = 11) +
  theme(legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 7))

p_temp <- ggplot(bkk_data, aes(month_name, temp, group = factor(year), colour = factor(year))) +
  geom_line(linewidth = 0.8) + scale_colour_manual(values = pal_seas, name = "Year") +
  labs(title = "Temperature (mean_tmax)", x = NULL, y = "°C") + theme_bw(base_size = 11) +
  theme(legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 7))

p_rh <- ggplot(bkk_data, aes(month_name, RH, group = factor(year), colour = factor(year))) +
  geom_line(linewidth = 0.8) + scale_colour_manual(values = pal_seas, name = "Year") +
  labs(title = "Vapour Pressure (mean_vap)", x = NULL, y = "kPa") + theme_bw(base_size = 11) +
  theme(legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 7))

p_rain <- ggplot(bkk_data, aes(month_name, rain, group = factor(year), colour = factor(year))) +
  geom_line(linewidth = 0.8) + scale_colour_manual(values = pal_seas, name = "Year") +
  labs(title = "Rainfall (mean_ppt)", x = NULL, y = "mm") + theme_bw(base_size = 11) +
  theme(legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 7))

fig_seas <- (p_cases | p_temp) / (p_rh | p_rain) +
  plot_annotation(
    title    = "Seasonal Patterns: Bangkok 2011–2019",
    subtitle = "Each line = one year (Jan–Dec)"
  )

print(fig_seas)
ggsave(filename = "Step1_Seasonal_Patterns_Bangkok.png", plot = fig_seas, 
       width = 10, height = 8, dpi = 300, bg = "white")


# ==============================================================================
beta0_base <- 1.5
rho_base   <- 0.003

S_star_base <- gamma_val / beta0_base
I_star_base <- omega_val * (1 - S_star_base) / (beta0_base * S_star_base + omega_val + gamma_val)
E_star_base <- gamma_val * I_star_base / sigma_val
R_star_base <- max(1 - S_star_base - E_star_base - I_star_base, 0)

state0_base <- c(S = round(S_star_base * N_pop),
                 E = round(E_star_base * N_pop),
                 I = round(I_star_base * N_pop),
                 R = round(R_star_base * N_pop),
                 C = 0)

seirs_ode_base <- function(t, state, parms) {
  S <- state[1]; E <- state[2]; I <- state[3]; R <- state[4]; C <- state[5]
  Neff   <- S + E + I + R
  beta_t <- parms[1]
  dS <- -beta_t * S * I / Neff + omega_val * R
  dE <-  beta_t * S * I / Neff - sigma_val * E
  dI <-  sigma_val * E         - gamma_val * I
  dR <-  gamma_val * I         - omega_val * R
  dC <-  sigma_val * E
  list(c(dS, dE, dI, dR, dC))
}

out_base <- as.data.frame(
  ode(y = state0_base, times = t_days, func = seirs_ode_base, 
      parms = c(beta0_base), method = "lsoda"))

cases_sim_base <- c(0, pmax(diff(out_base[, "C"]), 0)) * rho_base

fig_base <- ggplot() +
  geom_line(data = bkk_data, aes(date, cases, colour = "Actual Cases"), linewidth = 0.9) +
  geom_line(data = data.frame(date = bkk_data$date, y = cases_sim_base),
            aes(date, y, colour = "SEIRS β=1.5 (constant)"), linewidth = 0.8) +
  scale_colour_manual(values = c("Actual Cases" = "black", "SEIRS β=1.5 (constant)" = "red3"), name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Deterministic SEIRS: Constant β",
       subtitle = "σ=1/3.5 d⁻¹, γ=1/7 d⁻¹, ω=1/182 d⁻¹ | R₀=β/γ=10.5",
       x = NULL, y = "Cases") +
  theme_bw(base_size = 12) + theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

print(fig_base)
ggsave(filename = "Step2_SEIRS_Baseline_Constant_Beta.png", plot = fig_base, 
       width = 10, height = 6, dpi = 300, bg = "white")


# ==============================================================================

beta0_init <- 1.4
rho_init   <- 0.003

S_star_clim <- gamma_val / beta0_init
I_star_clim <- omega_val * (1 - S_star_clim) / (beta0_init * S_star_clim + omega_val + gamma_val)
E_star_clim <- gamma_val * I_star_clim / sigma_val
R_star_clim <- max(1 - S_star_clim - E_star_clim - I_star_clim, 0)

state0_clim <- c(S = round(S_star_clim * N_pop),
                 E = round(E_star_clim * N_pop),
                 I = round(I_star_clim * N_pop),
                 R = round(R_star_clim * N_pop),
                 C = 0)

seirs_ode_clim <- function(t, state, parms, Tz_vec, RHz_vec, Rz_vec, t_ref) {
  S <- state[1]; E <- state[2]; I <- state[3]; R <- state[4]; C <- state[5]
  Neff   <- S + E + I + R
  idx    <- which.min(abs(t_ref - t))
  beta_t <- max(parms[1] * exp(parms[2] * Tz_vec[idx] + parms[3] * RHz_vec[idx] + parms[4] * Rz_vec[idx]), 1e-9)
  dS <- -beta_t * S * I / Neff + omega_val * R
  dE <-  beta_t * S * I / Neff - sigma_val * E
  dI <-  sigma_val * E         - gamma_val * I
  dR <-  gamma_val * I         - omega_val * R
  dC <-  sigma_val * E
  list(c(dS, dE, dI, dR, dC))
}

run_seirs_clim <- function(b0, a1, a2, a3) {
  out <- tryCatch(
    as.data.frame(ode(y = state0_clim, times = t_days, func = seirs_ode_clim,
                      parms = c(b0, a1, a2, a3), method = "lsoda",
                      Tz_vec = temp_z, RHz_vec = rh_z, Rz_vec = rain_z, t_ref = t_days)),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA, length(t_days)))
  cases_raw <- c(0, pmax(diff(out[, "C"]), 0)) * rho_init
  sf        <- mean(bkk_data$cases) / max(mean(cases_raw), 1e-9)
  cases_raw * sf
}

pred_temp <- run_seirs_clim(beta0_init, 1.8,  0.0,  0.0)
pred_rh   <- run_seirs_clim(beta0_init, 0.0,  1.5,  0.0)
pred_rain <- run_seirs_clim(beta0_init, 0.0,  0.0, -1.4)
pred_full <- run_seirs_clim(beta0_init, 1.8,  1.5, -1.4)

df_clim_models <- data.frame(
  date      = rep(bkk_data$date, 4),
  cases_sim = c(pred_temp, pred_rh, pred_rain, pred_full),
  model     = rep(c("Temperature only", "Humidity only", "Rainfall only", "Temp+RH+Rain (Full)"), each = nrow(bkk_data))
)

df_clim_all <- bind_rows(data.frame(date = bkk_data$date, cases_sim = bkk_data$cases, model = "Actual Cases"), df_clim_models) %>%
  mutate(model = factor(model, levels = c("Actual Cases", "Temperature only", "Humidity only", "Rainfall only", "Temp+RH+Rain (Full)")))

colour_map_clim <- c("Actual Cases" = "black", "Temperature only" = "darkorange", 
                     "Humidity only" = "steelblue", "Rainfall only" = "seagreen", "Temp+RH+Rain (Full)" = "red3")

fig_clim <- ggplot(df_clim_all, aes(date, cases_sim, colour = model)) +
  geom_line(linewidth = 0.85) + facet_wrap(~ model, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_colour_manual(values = colour_map_clim, guide = "none") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "SEIRS with Climate Forcing β(t)", subtitle = "β(t) = β₀ × exp(α₁·Tz + α₂·RHz + α₃·Rainz) | lag: Temp=3mo, RH=2mo, Rain=0mo", x = NULL, y = "Cases") +
  theme_bw(base_size = 11) + theme(strip.text.y = element_text(size = 9, face = "bold"), strip.background = element_rect(fill = "grey90"), axis.text.x = element_text(angle = 45, hjust = 1))

print(fig_clim)
ggsave(filename = "Step3_SEIRS_Climate_Forcing_Beta.png", plot = fig_clim, 
       width = 10, height = 10, dpi = 300, bg = "white")


# ==============================================================================

run_ode_mle <- function(b0, a1, a2, a3) {
  out <- tryCatch(
    as.data.frame(ode(y = state0_clim, times = t_days, func = seirs_ode_clim,
                      parms = c(b0, a1, a2, a3), method = "lsoda",
                      Tz_vec = temp_z, RHz_vec = rh_z, Rz_vec = rain_z, t_ref = t_days)),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA, length(t_days)))
  c(0, pmax(diff(out[, "C"]), 0))
}

nll_mle <- function(log_b0, a1, a2, a3, log_rho) {
  b0      <- exp(log_b0)
  rho_fit <- plogis(log_rho)
  new_I   <- run_ode_mle(b0, a1, a2, a3)
  if (any(is.na(new_I))) return(1e9)
  mu <- pmax(new_I * rho_fit, 1e-9)
  -sum(dpois(bkk_data$cases, lambda = mu, log = TRUE))
}

cat("\nFitting SEIRS with mle2()...\n")
fit_mle <- tryCatch(
  mle2(nll_mle,
       start   = list(log_b0 = log(1.4), a1 = 0.30, a2 = 0.20, a3 = -0.10, log_rho = qlogis(0.03)),
       method  = "Nelder-Mead", control = list(maxit = 10000, reltol = 1e-6)),
  error = function(e) { cat("mle2 error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(fit_mle)) {
  cat("\n"); print(summary(fit_mle))
  cf_mle      <- coef(fit_mle) 
  b0_hat_mle  <- exp(cf_mle["log_b0"])
  a1_hat_mle  <- cf_mle["a1"]
  a2_hat_mle  <- cf_mle["a2"]
  a3_hat_mle  <- cf_mle["a3"]
  rho_hat_mle <- plogis(cf_mle["log_rho"])
} else {
  cat("mle2 failed → Using initial guesses\n")
  b0_hat_mle <- 1.4; a1_hat_mle <- 0.30; a2_hat_mle <- 0.20; a3_hat_mle <- -0.10; rho_hat_mle <- 0.03
}

cases_sim_mle <- run_ode_mle(b0_hat_mle, a1_hat_mle, a2_hat_mle, a3_hat_mle) * rho_hat_mle

df_mle_all <- bind_rows(
  data.frame(date = bkk_data$date, cases_sim = bkk_data$cases, model = factor("Actual Cases", levels = c("Actual Cases", "SEIRS MLE Fit"))),
  data.frame(date = bkk_data$date, cases_sim = cases_sim_mle, model = factor("SEIRS MLE Fit", levels = c("Actual Cases", "SEIRS MLE Fit")))
)

fig_mle <- ggplot(df_mle_all, aes(date, cases_sim, colour = model)) +
  geom_line(linewidth = 0.85) + facet_wrap(~ model, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_colour_manual(values = c("Actual Cases" = "black", "SEIRS MLE Fit" = "steelblue"), guide = "none") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "SEIRS Fitted Model (Poisson MLE)", subtitle = sprintf("β₀=%.3f, α₁=%.3f, α₂=%.3f, α₃=%.3f, ρ=%.5f", b0_hat_mle, a1_hat_mle, a2_hat_mle, a3_hat_mle, rho_hat_mle), x = NULL, y = "Cases") +
  theme_bw(base_size = 11) + theme(strip.text.y = element_text(size = 9, face = "bold"), strip.background = element_rect(fill = "grey90"), axis.text.x = element_text(angle = 45, hjust = 1))

print(fig_mle)
ggsave(filename = "Step4_SEIRS_Fitted_Poisson_MLE.png", plot = fig_mle, 
       width = 10, height = 8, dpi = 300, bg = "white")


# ==============================================================================

run_stoch_sim <- function(b0, a1, a2, a3, rho) {
  new_I <- run_ode_mle(b0, a1, a2, a3)
  if (any(is.na(new_I))) return(rep(NA, length(t_days)))
  rpois(length(new_I), lambda = pmax(new_I * rho, 1e-9))
}

set.seed(42)
stoch_matrix <- replicate(100, {
  b0_i  <- b0_hat_mle  * exp(rnorm(1, 0, 0.05))
  a1_i  <- a1_hat_mle  + rnorm(1, 0, 0.01)
  a2_i  <- a2_hat_mle  + rnorm(1, 0, 0.01)
  a3_i  <- a3_hat_mle  + rnorm(1, 0, 0.01)
  rho_i <- rho_hat_mle * exp(rnorm(1, 0, 0.05))
  run_stoch_sim(b0_i, a1_i, a2_i, a3_i, rho_i)
})

df_stoch <- data.frame(
  date = bkk_data$date,
  med  = apply(stoch_matrix, 1, median,   na.rm = TRUE),
  lo95 = apply(stoch_matrix, 1, quantile, 0.025, na.rm = TRUE),
  hi95 = apply(stoch_matrix, 1, quantile, 0.975, na.rm = TRUE)
)

fig_stoch <- ggplot() +
  geom_ribbon(data = df_stoch, aes(date, ymin = lo95, ymax = hi95, fill = "95% CI"), alpha = 0.25) +
  geom_line(data = df_stoch, aes(date, med, colour = "Sim median"), linewidth = 0.8) +
  geom_line(data = bkk_data, aes(date, cases, colour = "Actual Cases"), linewidth = 0.9) +
  scale_fill_manual(values = c("95% CI" = "steelblue"), name = NULL) +
  scale_colour_manual(values = c("Actual Cases" = "black", "Sim median" = "steelblue"), name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Stochastic SEIRS: 100 Simulations", subtitle = "Shaded = 95% envelope | Poisson observation noise | ±5% parameter perturbation", x = NULL, y = "Cases") +
  theme_bw(base_size = 12) + theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

print(fig_stoch)
ggsave(filename = "Step5_Stochastic_SEIRS_95CI.png", plot = fig_stoch, 
       width = 10, height = 6, dpi = 300, bg = "white")


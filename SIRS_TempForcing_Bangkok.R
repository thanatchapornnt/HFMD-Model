library(deSolve)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

raw <- read.csv("Data-of-hfmd-cases-and-terraclimate-2011-to-2022raw.csv",
                stringsAsFactors = FALSE)

bkk <- raw %>%
  filter(province == "Bangkok") %>%
  mutate(
    year  = as.integer(year),
    month = as.integer(month)
  ) %>%
  filter(year >= 2014 & year <= 2022) %>%
  arrange(year, month) %>%
  mutate(
    time_yr   = (year - 2014) + (month - 1) / 12,
    temp_mean = (mean_tmax + mean_tmin) / 2,
    semester  = case_when(
      month %in% c(5:9)        ~ "Semester 1",
      month %in% c(4, 10)      ~ "Break",
      month %in% c(1:3, 11:12) ~ "Semester 2"
    )
  )

cat("=== จำนวนแถวข้อมูล Bangkok 2014-2022 ===\n")
print(nrow(bkk))
cat("\n=== ตัวอย่างข้อมูล ===\n")
print(head(bkk[, c("year","month","cases","temp_mean","semester")]))

temp_monthly <- bkk %>%
  mutate(time_days = round(time_yr * 365)) %>%
  select(time_days, temp_mean)

temp_interpolator <- approxfun(x = temp_monthly$time_days,
                               y = temp_monthly$temp_mean,
                               rule = 2)

mean_T <- mean(temp_monthly$temp_mean, na.rm = TRUE)
sd_T   <- sd(temp_monthly$temp_mean,   na.rm = TRUE)

cat("\n=== อุณหภูมิ กทม. (ค่าเฉลี่ย/SD) ===\n")
cat("mean_T =", round(mean_T, 2), "°C  |  sd_T =", round(sd_T, 2), "°C\n")

MaxYears <- 9
shades_s1 <- data.frame()
shades_s2 <- data.frame()
shades_br <- data.frame()
for (y in 0:(MaxYears - 1)) {
  shades_s1 <- rbind(shades_s1, data.frame(xmin = y + 136/365, xmax = y + 273/365))
  shades_s2 <- rbind(shades_s2,
                     data.frame(xmin = y + 0/365,   xmax = y + 90/365),
                     data.frame(xmin = y + 305/365, xmax = y + 365/365))
  shades_br <- rbind(shades_br,
                     data.frame(xmin = y + 90/365,  xmax = y + 136/365),
                     data.frame(xmin = y + 273/365, xmax = y + 305/365))
}

SIRS_temp_forcing <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    temp_z <- (temp_interpolator(t) - mean_T) / sd_T
    beta_t <- max(Beta0 + B1 * temp_z, 0.01)
    
    dS <- mu - beta_t * S * I - mu * S + omega * R
    dI <- beta_t * S * I - gamma * I - mu * I
    dR <- gamma * I - mu * R - omega * R
    
    list(c(dS, dI, dR))
  })
}

beta0 <- 1.6
beta1 <- 0.4
gamma <- 1 / 7
mu    <- 1 / (70 * 365)
omega <- 1 / (365 * 1.5)

S0 <- 0.05; I0 <- 1e-4; R0_init <- 1 - S0 - I0
y0    <- c(S = S0, I = I0, R = R0_init)
times <- seq(0, MaxYears * 365, by = 1)
pars  <- c(Beta0 = beta0, B1 = beta1, gamma = gamma, mu = mu, omega = omega)

out <- as.data.frame(ode(y = y0, times = times, func = SIRS_temp_forcing, parms = pars))

max_actual       <- max(bkk$cases, na.rm = TRUE)
bkk$cases_scaled <- (bkk$cases / max_actual) * max(out$I) * 1.05

plot_sirs <- function(data, column, label, line_color) {
  p <- ggplot() +
    geom_rect(data = shades_s1, aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
              fill="steelblue", alpha=0.13) +
    geom_rect(data = shades_s2, aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
              fill="royalblue", alpha=0.07) +
    geom_rect(data = shades_br, aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
              fill="orange", alpha=0.10) +
    geom_line(data = data,
              aes(x = time / 365, y = .data[[column]]),
              linewidth = 0.8, color = line_color) +
    scale_x_continuous(limits = c(0, MaxYears), breaks = 0:MaxYears,
                       labels = 2014:(2014 + MaxYears), expand = c(0, 0)) +
    labs(title = label, x = "Year", y = "Proportion") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          axis.text.x  = element_text(angle = 45, hjust = 1))
  
  if (column == "I") {
    p <- p +
      geom_line(data = bkk,
                aes(x = time_yr, y = cases_scaled),
                color = "black", linewidth = 0.85)
  }
  p
}

p1 <- plot_sirs(out, "S", "Susceptible (S)", "green4") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p2 <- plot_sirs(out, "I", "Infectious (I) — Red = model, Black = actual cases", "red3") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p3 <- plot_sirs(out, "R", "Recovered (R)", "blue3")

final_plot <- (p1 / p2 / p3) +
  plot_annotation(
    title    = "HFMD SIRS Simulation — Bangkok (2014–2022)",
    subtitle = "Blue = Semester 1 | Light blue = Semester 2 | Orange = Break\nRed line = SIRS model | Black line = actual reported cases"
  )

ggsave("SIRS_TempForcing_Bangkok_plot.png", plot = final_plot, width = 10, height = 10, dpi = 300)
print(final_plot)

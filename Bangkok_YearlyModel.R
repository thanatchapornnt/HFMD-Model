library(deSolve)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)

raw <- read.csv("Data-of-hfmd-cases-and-terraclimate-2011-to-2022raw.csv",
                stringsAsFactors = FALSE)

bkk_raw <- raw %>%
  filter(province == "Bangkok") %>%
  mutate(year = as.integer(year), month = as.integer(month)) %>%
  arrange(year, month)

bkk_monthly_data <- split(bkk_raw$cases, bkk_raw$year)
bkk_monthly_data <- lapply(bkk_monthly_data, function(x) as.numeric(x))

SIR_vibrant <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    t_mod <- t %% 365
    is_term_time <- (t_mod >= 135 & t_mod <= 273) | (t_mod >= 304 | t_mod <= 60)
    force  <- if(is_term_time) 1.6 else -1.2
    beta_t <- Beta0 + B1 * force
    dS <- mu - beta_t * S * I - mu * S
    dI <- beta_t * S * I - gamma * I - mu * I
    dR <- gamma * I - mu * R
    list(c(dS, dI, dR))
  })
}

generate_yearly_plot <- function(year_str) {
  if (is.null(bkk_monthly_data[[year_str]])) return(NULL)
  actual_monthly <- bkk_monthly_data[[year_str]]
  if (length(actual_monthly) < 12) return(NULL)
  max_cases <- max(actual_monthly)

  pars <- c(Beta0 = 3.2, B1 = 1.5, gamma = 1/13, mu = 1/(3 * 365))
  y0   <- c(S = 0.35, I = 0.01, R = 0.64)

  out <- as.data.frame(ode(y = y0, times = seq(0, 365, by = 1),
                           func = SIR_vibrant, parms = pars))

  scale_factor <- max_cases / max(out$I)
  out$S_plot <- out$S * scale_factor * 1.8
  out$I_plot <- out$I * scale_factor
  out$R_plot <- out$R * scale_factor

  out_long <- pivot_longer(out, cols = c(S_plot, I_plot, R_plot),
                           names_to = "Group", values_to = "Count")

  ggplot() +
    geom_rect(aes(xmin = 5.5, xmax = 9.5, ymin = -Inf, ymax = Inf), fill = "skyblue", alpha = 0.2) +
    geom_rect(aes(xmin = 11,  xmax = 12.5, ymin = -Inf, ymax = Inf), fill = "skyblue", alpha = 0.2) +
    geom_rect(aes(xmin = 0.5, xmax = 2.5,  ymin = -Inf, ymax = Inf), fill = "skyblue", alpha = 0.2) +
    geom_line(data = out_long,
              aes(x = (time/365)*12 + 0.5, y = Count, color = Group), linewidth = 1) +
    geom_point(data = data.frame(m = 1:12, v = actual_monthly),
               aes(x = m, y = v), color = "black", size = 1.5) +
    scale_color_manual(values = c("S_plot" = "green4", "I_plot" = "red3", "R_plot" = "blue3"),
                       labels = c("S (Susceptible)", "I (Infectious)", "R (Recovered)")) +
    scale_x_continuous(breaks = 1:12, labels = month.abb, limits = c(0.5, 12.5)) +
    labs(title = paste("Year:", year_str), x = NULL, y = "Cases") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(size = 7))
}

plots <- lapply(names(bkk_monthly_data), generate_yearly_plot)
plots <- plots[!sapply(plots, is.null)]

if (length(plots) > 0) {
  final_layout <- wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title    = "Bangkok HFMD SIR Model Simulation (2011-2022)",
      subtitle = "Blue zones = School term | White zones = Holidays",
      caption  = "Red line: Model I | Black dots: Actual reported cases"
    )
  ggsave("Bangkok_YearlyModel_plot.png", plot = final_layout, width = 14, height = 16, dpi = 300)
  print(final_layout)
} else {
  print("No data to plot")
}

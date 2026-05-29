library(deSolve)
library(ggplot2)
library(patchwork)

term = function(t){
  t = t %% 365
  if((t > 90 & t < 136) | (t > 273 & t < 305)){
    return(-1) 
  } else {
    return(1)  
  }
}

SIR_term_forcing <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    beta_t <- Beta0 + B1 * term(t)
    
    dS <- mu - beta_t * S * I - mu * S
    dI <- beta_t * S * I - gamma * I - mu * I
    dR <- gamma * I - mu * R
    
    list(c(dS, dI, dR))
  })
}

beta0 <- 17 / 13 
beta1 <- 0.225 
gamma <- 1 / 13 
mu    <- 1 / (70 * 365)

S0 <- 1/17 
I0 <- 1e-4 
R0 <- 1 - S0 - I0
y0 <- c(S = S0, I = I0, R = R0)

MaxYears <- 10 
times <- seq(0, MaxYears * 365, by = 1)
pars  <- c(Beta0 = beta0, B1 = beta1, gamma = gamma, mu = mu)

out <- as.data.frame(ode(y = y0, times = times, func = SIR_term_forcing, parms = pars))

shades <- data.frame()
for(year in 0:(MaxYears-1)) {
  shades <- rbind(shades, data.frame(xmin = year + 0/365,   xmax = year + 90/365))
  shades <- rbind(shades, data.frame(xmin = year + 136/365, xmax = year + 273/365))
  shades <- rbind(shades, data.frame(xmin = year + 305/365, xmax = year + 365/365))
}

plot_sir_final <- function(data, column, label, line_color, show_years = 3) {
  ggplot() +
    geom_rect(data = shades, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "deepskyblue", alpha = 0.2) +
    geom_line(data = data, aes(x = time / 365, y = .data[[column]]),
              linewidth = 0.8, color = line_color) +
    coord_cartesian(xlim = c(0, show_years)) + 
    xlab("Time (years)") + ylab(label) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, MaxYears, 0.5)) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(face = "bold")
    )
}

p1 <- plot_sir_final(out, "S", "Susceptibles", "green4")
p2 <- plot_sir_final(out, "I", "Infectious", "red2")
p3 <- plot_sir_final(out, "R", "Recovereds", "blue3")

final_plot <- (p1 / p2 / p3) + 
  plot_annotation(
    title = "SIR Model with School Term Forcing",
    subtitle = "Light Blue Zones = School Term (High Transmission) | White Zones = Holidays (Low Transmission)",
    caption = "Showing first 3 years for clarity"
  )

ggsave("SIR_TermForcing_plot.png", plot = final_plot, width = 10, height = 8, dpi = 300)
print(final_plot)

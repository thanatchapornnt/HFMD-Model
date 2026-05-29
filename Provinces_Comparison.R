library(tidyr)
library(ggplot2)

all_provinces_data <- data.frame(
  Year = 2014:2024,
  Bangkok   = c(8227, 5295, 11406, 9705, 10193, 9309, 2285, 790, 5902, 7347, 14842),
  Chonburi  = c(979,  570,  1301,  834,  1633,  1651, 640,  115, 2160, 993,  6199),
  ChiangMai = c(2072, 1536, 3563,  2635, 3164,  4051, 1614, 1067,3460, 3404, 2624),
  Songkhla  = c(1019, 891,  1007,  1149, 1115,  1087, 1740, 571, 1177, 1567, 1081),
  Prachuap  = c(604,  709,  433,   548,  727,   521,  353,  170, 485,  536,  624),
  KhonKaen  = c(1489, 555,  1483,  1179, 1327,  1185, 939,  649, 3289, 1030, 1631)
)

all_data_long <- pivot_longer(all_provinces_data,
                              cols = -Year,
                              names_to = "Province",
                              values_to = "Cases")

ggplot(all_data_long, aes(x = Year, y = Cases, color = Province, group = Province)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2014:2024) +
  labs(
    title    = "Comparison of HFM Cases by Province (2014-2024)",
    subtitle = "Actual historical case totals from multi-province surveillance data",
    x        = "Year (A.D.)",
    y        = "Total Number of Cases"
  ) +
  theme_minimal() +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 14),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

p_final <- last_plot()
ggsave("Provinces_Comparison_plot.png", plot = p_final, width = 12, height = 7, dpi = 300)

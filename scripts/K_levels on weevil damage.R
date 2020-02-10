# Additional analysis to consider
## splitting the dataset at 5% internal damage and checking for differences in N levels

##  differences between k deficient and K sufficient fields 
library(data.table)
library(broom)
library(msme)
library(tidyverse)
library(MASS)
library(jtools)
library(here)
library(export)

# source custom functions
source(here::here("scripts", "Custom functions.R"))
# visual elements
axis_text <- element_text(color = "black", size = 12)

##  prepare data for analysis ####

# import& drop NAs
Dt <- data.table::fread(here::here("data", "raw", "Baseline_cleaned and combined2.csv")) %>%
  na.omit()

# change str
Dt <- Dt %>%
  .[, District := as.factor(District)] %>%
  .[, SubCounty := as.factor(SubCounty)] %>%
  .[, Village := as.factor(Village)] %>%
  .[, FL := as.integer(FL)] %>%
  .[, Yls := as.integer(Yls)]

# recode the sigatoka values to normal scale and sort  ### calculate insl instead
Dt2 <- Dt %>%
  mutate(., FL = ifelse(FL > Yls, FL, Yls)) %>%
  as.data.table() %>%
  # Recode negative FLs
  mutate(., sgtk := ifelse(Yls == 0 | Yls == FL, Yls / FL, FL - Yls)) %>%
  as.data.table() %>%
  .[, Sgtk := sgtk / FL]
Dt3 <- Dt2
Dt3 <- Dt3 %>%
  mutate(Region = recode(District,
                         Isingiro = "Southwestern", Mbarara = "Southwestern",
                         Luwero = "Central", Nakaseke = "Central", Kabarole = "Western")) %>%
  as.data.table() %>%
  mutate(., N_level = as.factor(ifelse(N > 0.2, "N_sufficient", "N_deficient"))) %>%
  as.data.table()
Dt3 <- Dt3 %>%
  mutate(K_level = case_when(
    K < 0.96 ~ "K_deficient",   K >= 0.96 ~ "K_sufficient",
    Region == "Western" & K < 0.93 ~ "K_deficient",
    Region == "Western" & K >= 0.93 ~ "K_sufficient",
    Region == "Southwestern" & K < 0.81 ~ "K_deficient",
    Region == "Southwestern" & K >= 0.81 ~ "K_sufficient",
    Region == "Central" & K < 1.5 ~ "K_deficient",
    Region == "Central" & K >= 1.5 ~ "K_sufficient")) %>%
  as.data.table(.) %>%
  .[, K_level := as.factor(K_level)]

# use negative binomial model on a data subset with non-zero damage levels
Dt4 <- Dt3[XT > 0]

# influence of K_levels on weevil damage
New_model7 <- glm.nb(XT ~ K_level, data = Dt4)
tidy_New_model7 <- tidy(New_model7)
tidy_New_model7
k <- tidy(P__disp(New_model7))
l <- bind_rows(tidy_New_model7, k)
fwrite(l, here::here("Results", "Tables", "positve XT between K levels.csv"))
res_effectplot5 <- jtools::effect_plot(New_model7, pred = K_level, interval = TRUE, plot.points = TRUE)

g8 <- res_effectplot5 +
  labs(y = "Weevil damage (%)") +
  theme_bw(base_size = 12)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = axis_text,
        axis.ticks.x = element_blank(),
        axis.text.y = axis_text,
        title = axis_text)
ggsave(path = here::here("Results", "Figs"), filename = "Total weevil damage against K status_nonzero.eps", width = 5, height = 7)
g8
graph2doc(
  file = paste0(here::here("Results", "Figs"), "/", "Grand4.docx"),
  width = 4,
  height = 5
)

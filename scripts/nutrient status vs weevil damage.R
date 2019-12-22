##  packages and functions
library(data.table)       # enhanced version of data.frame
library(tidyverse)        # tidy data wrangling
library(msme)             # methods for statistical model estimation (P__disp)
library(jtools)           # model support for regression analyses (robust SE)
library(broom.mixed)      # tidying mixed model outputs
library(glmmTMB)          # fit generalised linear models
library(rticles)          # Article templates
library(xtable)           # convert an R object to xtable object printable as LaTeX table
library(stargazer)        # well formatted LaTeX and html code
library(here)             # manage director paths
source(here::here("scripts", "Custom functions.R"))

# import data
dat <- data.table::fread(here::here("data", "raw", "GT_trialdat2.csv"))
dat <- na.omit(dat) # remove NAs

# calculate weevil damage
data.table::setnames(dat, c("UXI", "LXI", "UXO", "LXO"), c("UI_damage", "LI_damage", "UO_damage", "LO_damage"))
get_weevil_damage(dat)

# change structure
dat <- dat %>%
  dplyr::mutate_if(is.character, as.factor)

## weevil damage against K and N individually with non_zero XT 
KNx2 <- glmmTMB(XT ~ K + N + (1 | Spray), data = dat["XT" > 0], family = nbinom2())
broom.mixed::tidy(KNx2)
tidy(KNx2) %>% fwrite(., here::here("Results", "Tables", "XT_vs_KN_exp.csv"))

## nutrient status vs weevil damage

fit4 <- glmmTMB(XT ~ Nutrient.status + (1 | Spray), data = dat["XT" > 0], family = nbinom2())
tidy(fit4) %>% fwrite(., here::here("Results", "Tables", "XT_vs_Nutrient_status_nonzero.csv"))
res_effectplot13 <- jtools::effect_plot(fit4, pred = Nutrient.status, interval = TRUE, plot.points = TRUE,allow.new.levels=TRUE)

g18 <- res_effectplot13 +
  labs(y = "Total weevil damage (%)") +
  ggtitle("Total weevil damage against Nutrient status_nonzero")
ggsave(path = here::here("Results", "Figs"), filename = "Total weevil damage against Nutrient status_nonzero.png", width = 6, height = 7)


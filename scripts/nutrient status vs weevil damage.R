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
library(MASS)             # statistical modeling
library(patchwork)        # combining gg objects
library(export)           # Export graphs to word
library(here)             # manage director paths
source(here::here("scripts", "Custom functions.R"))

# visual elements
axis_text <- element_text(color = "black", size = 12)

# data preparation------------
# import data
dat <- data.table::fread(here::here("data", "raw", "GT_trialdat2.csv"))
dat <- na.omit(dat) # remove NAs

# calculate weevil damage
data.table::setnames(dat, c("UXI", "LXI", "UXO", "LXO"), c("UI_damage", "LI_damage", "UO_damage", "LO_damage"))
get_weevil_damage(dat)

# change structure
dat <- dat %>%
  dplyr::mutate_if(is.character, as.factor)


## full data weevil damage vs K and N individually with non_zero XT-------------- 
KNx2 <- glmmTMB(XT ~ K + N + (1 | Spray), data = dat["XT" > 0], family = nbinom2())
broom.mixed::tidy(KNx2)
tidy(KNx2) %>% fwrite(., here::here("Results", "Tables", "XT_vs_KN_exp.csv"))

# interraction
KNx2_interaction <- glmmTMB(XT ~ K*N + (1 | Spray), data = dat["XT" > 0], family = nbinom2())
tidy_KNx2_interaction <- broom.mixed::tidy(KNx2_interaction)
fwrite(tidy_KNx2_interaction, here::here("Results", "Tables", "XT_vs_KN_exp_interaction.csv"))

models_compared2 <- tidy(myaddterm(KNx2, . ~ . + K:N, test="Chisq"))
models_compared2
fwrite(models_compared2, here::here("Results", "Tables", "models_compared.csv"))

## nutrient status vs weevil damage
fit4 <- glmmTMB(XT ~ Nutrient.status + (1 | Spray), data = dat["XT" > 0], family = nbinom2())
tidy(fit4) %>% fwrite(., here::here("Results", "Tables", "XT_vs_Nutrient_status_nonzero.csv"))
res_effectplot13 <- jtools::effect_plot(fit4, pred = Nutrient.status, interval = TRUE, plot.points = TRUE,allow.new.levels=TRUE)

g18 <- res_effectplot13 +
  labs(y = "Total weevil damage (%)") +
  ggtitle("Total weevil damage against Nutrient status_nonzero")
ggsave(path = here::here("Results", "Figs"), filename = "Total weevil damage against Nutrient status_nonzero.png", width = 6, height = 7)


# repeat without sprayed plots
# ## weevil damage against K and N individually with non_zero XT and no spray
# dat2 <- as.data.table (dat) %>% .[Spray==0]
# names(dat2)
# KN_model <- glm.nb(XT ~ K + N, data = dat2["XT" > 0])
# summary(KN_model)
# alias(KN_model) # the data is not sufficient
# #because K and N are perfectly correlated in this sample, the coeficients for N are defined due to singularities use alias() to check
# broom::tidy(KN_model)
# tidy(KN_model) %>% fwrite(., here::here("Results", "Tables", "XT_vs_KN_exp_no_spray.csv"))

## with spray , K and N individually and non_zero XT-------
dat3 <- as.data.table (dat) %>% .[Spray==1]
KN_model2 <- glm.nb(XT ~ K + N, data = dat3["XT" > 0])
summary(KN_model2)
tidy_KN_model2 <- broom::tidy(KN_model2)
P__disp(KN_model2)
fwrite(tidy_KN_model2, here::here("Results", "Tables", "XT_vs_KN_exp_spray.csv"))

# with interaction term
KN_model2_interaction <- glm.nb(XT ~ K*N, data = dat3["XT" > 0])
summary(KN_model2_interaction)
tidy_interaction <- tidy(KN_model2_interaction)
fwrite(tidy_interaction, here::here("Results", "Tables", "XT_vs_KN_exp_spray.csv"))

# compare the two using AIC
extractAIC(KN_model2,KN_model2_interaction)  # do  not include the interaction, AIC increases a lot.
models_compared <- tidy(myaddterm(KN_model2, . ~ . + K:N, test="Chisq"))
models_compared
fwrite(models_compared, here::here("Results", "Tables", "models_compared.csv"))


## both sprayed and non-sprayed plots: nutrient status vs weevil damage with of non-zero XT -------------------
dat$Nutrient.status <- relevel(dat$Nutrient.status, ref = "None")

fit4 <- glmmTMB(XT ~ Nutrient.status + (1 | Spray), data = dat["XT" > 0], family = nbinom2())
tidy_fit4 <- tidy(fit4) 
fwrite(tidy_fit4, here::here("Results", "Tables", "XT_vs_Nutrient_status_nonzero.csv"))
dat3 <- as.data.table (dat) %>% .[Spray==1]

fit5 <- glm.nb(XT ~ Nutrient.status, data = dat3)
tidy(fit5) %>% fwrite(., here::here("Results", "Tables", "XT_vs_Nutrient_status_nonzero_spray_alone.csv"))

# vector graphics---------
res_effectplot20 <- jtools::effect_plot(fit5, pred = Nutrient.status, interval = TRUE, plot.points = TRUE,allow.new.levels=TRUE)
g20 <- res_effectplot20 +
  geom_point( alpha = I(1) )+
  labs(y = "Weevil damage (%)", x= "Nutrient status") +
  ggtitle("B")+
  theme_bw(base_size = 12)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.text.x = axis_text,
        axis.ticks.x = element_blank(),
        axis.text.y = axis_text,
        title = axis_text)
ggsave(path = here::here("Results", "Figs"), filename = "Weevil damage vs Nutrient status_nonzero_no spray.eps", width = 4, height = 7)



# vector graphic---------------------------------
res_effectplot13 <- jtools::effect_plot(fit4, pred = Nutrient.status, interval = TRUE, plot.points = TRUE,allow.new.levels=TRUE)
g18 <- res_effectplot13 +
  geom_point( alpha = I(1) )+
  labs(y = "Weevil damage (%)") +
  ggtitle("A")+
  theme_bw(base_size = 12)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.text.x = axis_text,
        axis.ticks.x = element_blank(),
        axis.text.y = axis_text,
        title = axis_text)
ggsave(path = here::here("Results", "Figs"), filename = "Total weevil damage against Nutrient status_nonzero.png", width = 6, height = 7)

## on sprayed plots: nutrient status vs weevil damage with of non-zero XT -------------------
dat3 <- as.data.table (dat) %>% .[Spray==1]
KNx3 <- glm.nb(XT ~ K+N, data = dat3)
tidy_KNx3 <- tidy(KNx3)
summary(KNx3)
dat3$Nutrient.status <- relevel(dat3$Nutrient.status, ref = "None")


# combine graphs
Grand7 <- g18+g20
Grand7
graph2doc(
  file = paste0(here::here("Results", "Figs"), "/", "Grand7.docx"),
  width = 6,
  height = 8
)

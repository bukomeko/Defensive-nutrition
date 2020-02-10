##  load packages and functions 
library(tidyverse)            # tidy data wrangling
library(data.table)           # enhanced version of data.frame
library(MASS)                 # glm models
library(broom)                # tidying model outputs
library(msme)                 # methods for statistical model estimation (P__disp)
tlibrary(factoextra)           # multivariate data analysis and visualization
library(cluster)              # methods for cluster analysis
library(clValid)              # choose clustering alogarithm and number of clusters
library(jtools)               # model support for regression analyses (robust SE)
library(patchwork)            # combining plots
library(export)               # export graphs to word and office

# source custom functions
source(here::here("scripts", "Custom functions.R"))

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

# setting Key column for sorting
keycol <- c("District", "SubCounty", "Village", "BWT")
# sorted
dT <- Dt2 %>%
  setorderv(., keycol) %>%
  dplyr::select(., c(N, K))

##  choose cluster alogarithm, number of clusters, then cluster and visualize  ####

# alogarithm and number of clusters, say yes!

# Enhanced hierarchical clustering, cut in 2 groups
my_data <- scale(dT)
res.hc <- eclust(my_data, "hclust", k = 2, graph = FALSE)

# visual elements
axis_text <- element_text(color = "black", size = 12)

# Visualize
g10 <- fviz_dend(res.hc, rect = TRUE, show_labels = FALSE) +
  ggtitle("A")+
  theme_bw(base_size = 12)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = axis_text,
        title = axis_text)

dT1 <- cbind(Dt2, res.hc$cluster)
setnames(dT1, "V2", "Cluster"); g10

##  Describe the clusters ####
# Describe the cluster
dT1_a <- dT1[, .(N = mean(N), K = mean(K)), by = Cluster]

# melt to compare quickly
dT1_b <- dT1[, c("N", "K", "Cluster")]
dT1_b %>%
  data.table::melt(., c("Cluster"),
                   measure = c("N", "K"),
                   variable.name = "Nutrients", value.name = "Concetration"
  ) %>%
  dplyr::group_by(Nutrients) %>%
  do(broom::tidy(kruskal.test(x = .$Concetration, g = .$Cluster))) %>%
  .[, c("Nutrients", "p.value")] -> comparisons_NK

# transpose the data and keep row names
dcast(melt(comparisons_NK, id.vars = "Nutrients"), variable ~ Nutrients) -> comparisons_NK
rbind(dT1_a, comparisons_NK, use.names = FALSE) -> dT1_a
dT1_a
data.table::fwrite(dT1_a, here::here("Results", "Tables", "NK_new_clusters.csv"))

##  testing differences between clusters #####
# first convert cluster column to factor and test weevil damage fo non-zero weevil damage
dT1[, Cluster := as.factor(Cluster)]

New_model1b <- glm.nb(XT ~ Cluster, data = dT1[XT > 0])
tidy_New_model1b <- tidy(New_model1b)
f <- tidy(P__disp(New_model1b))
g <- bind_rows(tidy_New_model1b, f)
fwrite(g, here::here("Results", "Tables", "XT between NK_new_clusters_nonzero.csv"))
res_effectplot7 <- jtools::effect_plot(New_model1b, pred = Cluster, interval = TRUE, plot.points = TRUE)
g12 <- res_effectplot7 +
  labs(y = "Weevil damage (%)") +
  ggtitle("B")+
  theme_bw(base_size = 12)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = axis_text,
        axis.text.y = axis_text,
        title = axis_text)

ggsave(path = here::here("Results", "Figs"), filename = "Total weevil damage against farm types_nonzero.eps", width = 6, height = 7)

Grand2 <- g10 | g12
Grand2
graph2doc(
  file = paste0(here::here("Results", "Figs"), "/", "Grand2.docx"),
  width = 6,
  height = 4
)

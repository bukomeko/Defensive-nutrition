# Notes--------
# fix the fonts

##  Packages and custom functions used-----
library(tidyverse)            # the world of tidy data wrangling
library(data.table)           # data structure
library(MASS)                 # modelling negative binomial
library(broom)                # summaring model info into tidy tibbles
library(msme)                 # dispersion parameter for poisson family
library(patchwork)            # combining plots
library(export)               # export graphs to word and office
library(extrafont)            # fonts
library(formatR)              # code readability

# fetch my custom functions
source(here::here("scripts", "Custom functions.R"))

##  Data preparation    ####

# yield data
DT <- fread(here::here("data", "raw", "Baseline_Yield data.csv"))

# Choose variables of interest
yield <- DT[, c(
  "ID",
  "Region",
  "Farm ID",
  "Plot ID",
  "Dist (m)",
  "Mat No.",
  "Plant No.",
  "Girth_0 (cm)",
  "Girth_100 (cm)",
  "# Hands",
  "# Fingers"
)]
yield <- na.omit(yield)
setnames(
  yield,
  c(
    "Dist (m)",
    "Girth_0 (cm)",
    "Girth_100 (cm)",
    "# Hands",
    "# Fingers"
  ),
  c("Dist", "Girth0", "Girth100", "Hands", "Fingers")
)
get_yield(yield)


##  Distance vs yield   ####
# # visualize first
# ggplot(data = yield, aes(Dist, BWT)) +
#   geom_point(position = "jitter") +
#   theme_bw()

# Get boundary points
yield %>%
  .[BWT < 55] %>%
  .[Dist < 100] %>% # drop outliers
  .[, Rel_BWT := (BWT / (max(.[, BWT])))] -> yield # relative Bunch Weight

# split the data into top 10% yielders
yield_10 <- yield %>% dplyr::filter(Rel_BWT > 0.9)
ymax <- mean(yield_10$Rel_BWT)

# # visualize again
# ggplot(data = yield, aes(Dist, Rel_BWT)) +
#   geom_point(position = "jitter") +
#   theme_bw()
#
yield %>%
  setorder(-Dist) %>% # sorting
  .[, Bpts := get_bpts(Rel_BWT)] %>% # picking bpts
  unique(., by = "Bpts") %>% # picking unique ones
  na.omit() -> Bpts_F # Droping NAs $ renaming

# fit nls model to Boundary points
m1 <-
  nls(Bpts ~ ymax / (1 + k * exp(-R * Dist)), start = list(k = 1, R = 8.35e-5), Bpts_F)
z <- tidy(m1)
m2 <- nls(
  Bpts ~ 1 / (1 + k * exp(-R * Dist)),
  start = list(k = z$estimate[1], R = z$estimate[2]),
  data = Bpts_F
)
formula(m2)
tidy_m1 <- tidy(m2)
fwrite(tidy_m1, here::here("Results", "Tables", "Dist_vs_yield.csv"))

# use the model to predict
Bpts_F[, pred := predict(m2)]

# visual elements
axis_text <- element_text(color = "black", size = 12)

# Visualize g1---------
g1 <- ggplot(Bpts_F) +
  geom_point(aes(Dist, Bpts)) +
  geom_smooth(aes(Dist, pred), se = F, colour = "red") +
  geom_point(
    data = yield[BWT > 0],
    aes(Dist, Rel_BWT),
    colour = "black",
    position = "jitter"
  ) +
  expand_limits(x = 0, y = 0) +
  labs(x = "Distance from homestead (m)", y = "Relative Bunch weight") +
  ggtitle("A") +
  theme_bw(base_size = 12, base_family = "TT Courier New") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = axis_text,
        axis.text.y = axis_text,
        title = axis_text)
ggsave(
  path = here::here("Results", "Figs"),
  filename = "BWT_vs_Dist.eps",
  width = 7,
  height = 6
)

##  Distance vs weevil  ####

# weevil data
weevil <-
  fread(here::here("data", "raw", "Baseline_weevil data.csv"))
setnames(
  weevil,
  c(
    "Dist (m)",
    "UXI (% Damage)",
    "UXO (% Damage)",
    "LLI (% Damage)",
    "LXO (% Damage)"
  ),
  c("Distw", "UI_damage", "UO_damage", "LI_damage", "LO_damage"),
  skip_absent = T
)
# Choose the variables of interest
weevil <-
  weevil[, c("ID",
             "Region",
             "Distw",
             "UI_damage",
             "UO_damage",
             "LI_damage",
             "LO_damage")]
weevil <- na.omit(weevil)
get_weevil_damage(weevil)

# # Visualize
# ggplot(data = weevil, aes(Distw, XI, colour = Region, label = ID)) +
#   geom_point(position = "jitter") +
#   geom_text(aes(label = ID), hjust = 0, vjust = 0)+
#   theme_bw()
#
# ggplot(data= weevil[XT<30], aes(Distw ,XT)) +
#     geom_point(position = "jitter")+
#     theme_bw()

# drop outliers
weevil <-
  weevil[XI < 30] %>% .[ID != 237] %>% .[ID != 1] %>% .[Distw < 150] %>%
  # weevil <- weevil[XI < 40] %>%
  # .[ID != 237] %>%
  # .[ID != 1] %>%
  .[Distw < 150]

# split the data into top 10%
p <- max(weevil$XI) * 0.9
Weevil_10 <- weevil %>% dplyr::filter(XI > p)
ymax <- mean(Weevil_10$XI)

# get bpts
weevil %>%
  setorder(.,-Distw) %>% # Order Ascending
  .[, Bpts := get_bpts(XI)] %>% # calculating Boundary points
  unique(., by = "Bpts") %>% # retain unique values
  na.omit() -> Bpts_F # drop NAs

# fit nls
m3 <- nls(Bpts ~ ymax / (1 + k * exp(-R * Distw)),
          start = list(k = -0.624574, R = -0.005247),
          data = Bpts_F)
z <- tidy(m3)
m4 <- nls(Bpts ~ 1 / (1 + k * exp(-R * Distw)),
          start = list(k = -0.9702669, R = 0.0012308),
          data = Bpts_F)
tidy_m4 <- tidy(m4)
fwrite(tidy_m4, here::here("Results", "Tables", "XT_vs_dist.csv"))
formula(m4)

Bpts_F[, Pred := predict(m4)]

# visualize g2--------
g2 <- ggplot(Bpts_F) +
  geom_point(aes(Distw, Bpts)) +
  geom_smooth(aes(Distw, Pred), se = F, colour = "red") +
  geom_point(data = weevil, aes(Distw, XI), position = "jitter") +
  expand_limits(x = 0, y = 0) +
  labs(x = "Distance from homestead (m)", y = "weevil damage(%)") +
  ggtitle("B") +
  theme_bw(base_size = 12, base_family = "TT Courier New") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = axis_text,
        axis.text.y = axis_text,
        title = axis_text)


ggsave(
  path = here::here("Results", "Figs"),
  filename = "damage_vs_Dist.eps",
  width = 7,
  height = 6
)


##  Weevil vs N ####
data.table::fread(here::here("data", "raw", "Baseline_cleaned and combined2.csv")) %>%
  na.omit() -> Nitrogen1

# drop some points
Nitrogen1 <-
  Nitrogen1[N < 0.3] %>% .[XT < 9] %>% .[FarmID != "Mwesigye Alkarito"]

# # Visualize
# ggplot(data = Nitrogen1, aes(N,XT)) +
#   geom_point(position = "jitter") +
#   labs(x="Soil Nitrogen (%)", y= "Weevil damage(%)",
#        title="Soil nitrogen against weevil damage in Bananas")+
#   theme_bw(base_size = 16)+
#   expand_limits(x = 0, y = 0) +
#   theme_bw(base_size = 12)
#

# split the data into top 10%
p <- max(Nitrogen1$XT) * 0.9
Nitrogen1_10 <- Nitrogen1 %>% dplyr::filter(XT > p)
ymax <- mean(Nitrogen1_10$XT)

# get bpts
Nitrogen1 %>%
  setorder(., N) %>%             # Order Ascending
  .[, Bpts := get_bpts(XT)] %>%   # calculating Boundary points
  unique(., by = "Bpts") %>%      # retain unique values
  na.omit(.) -> Bpts_F            # drop NAs

# fit nls to the boundary points
mmm <- nls(Bpts ~ ymax / (1 + k * exp(-R * N)),
           start = list(k = 1, R = -0.05693),
           data = Bpts_F)
tidy_mmmm <- tidy(mmm)

mmm2 <- nls(
  Bpts ~ ymax / (1 + k * exp(-R * N)),
  start = list(k = tidy_mmmm$estimate[1], R = tidy_mmmm$estimate[2]),
  data = Bpts_F
)
tidy_mmmm2 <- tidy(mmm2)

# predict usig m5
Bpts_F[, Pred := predict(mmm2)]

# visualize ggg---------
ggg <- ggplot(Bpts_F) +
  geom_point(aes(N, XT)) +
  geom_smooth(aes(N, Pred), se = F, colour = "red") +
  geom_point(data = Nitrogen1, aes(N, XT), position = "jitter") +
  expand_limits(x = 0.1, y = 0) +
  labs(x = "Soil Nitrogen (%)", y = "Weevil damage(%)") +
  ggtitle("F") +
  theme_bw(base_size = 12, base_family = "TT Courier New") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = axis_text,
        axis.text.y = axis_text,
        title = axis_text)
ggsave(
  path = here::here("Results", "Figs"),
  filename = "Nitrogen vs XT.eps",
  width = 10,
  height = 6
)


##  Distance vs N ####
data.table::fread(here::here("data", "raw", "Baseline_cleaned and combined2.csv")) %>%
  na.omit() -> Nitrogen

# # Visualize
# ggplot(data = Nitrogen, aes(Distance, N, label = FarmID)) +
#   geom_point(position = "jitter") +
#   # geom_text(aes(label = FarmID),hjurst = 0,vjust=0)
#   theme_bw()
# drop some points
Nitrogen <- Nitrogen[N < 0.25] %>% .[Distance > 5]

# split the data into top 10%
p <- max(Nitrogen$N) * 0.9
Nitrogen_10 <- Nitrogen %>% dplyr::filter(N > p)
ymax <- mean(Nitrogen_10$XI)

# get bpts
Nitrogen %>%
  setorder(.,-Distance) %>% # Order Ascending
  .[, Bpts := get_bpts(N)] %>% # calculating Boundary points
  unique(., by = "Bpts") %>% # retain unique values
  na.omit(.) -> Bpts_F # drop NAs

# fit nls to the boundary points does not fit well
m5 <- nls(Bpts ~ ymax / (1 + k * exp(-R * Distance)),
          start = list(k = 0.02024, R = -0.05693),
          data = Bpts_F)
summary(m5)

# predict usig m5
Bpts_F[, Pred := predict(m5)]

# fit a polynomial
p <- lm(Bpts ~ poly(Distance, 2, raw = T), data = Bpts_F)
N_vs_dist <- tidy(p)
fwrite(N_vs_dist, here::here("Results", "Tables", "N_vs_dist.csv"))
p_augment <- augment(p)
setnames(p_augment, "poly.Distance..2..raw...T.", "poly")

# Visualize g3-----
g3 <- ggplot(Bpts_F) +
  geom_point(aes(Distance, Bpts), position = "jitter") +
  theme(panel.background = element_rect(fill = "white"))+
  geom_smooth(data = p_augment, aes(poly[, 1], .fitted), se = FALSE, colour = "red") +
  geom_point(data = Nitrogen, aes(Distance, N), position = "jitter") +
  expand_limits(x = 7, y = 0.1) +
  labs(x = "Distance from homestead (m)", y = "Soil Nitrogen (%)") +
  ggtitle("C") +
  theme_bw(base_size = 12, base_family = "TT Courier New") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = axis_text,
        axis.text.y = axis_text,
        title = axis_text)
ggsave(
  path = here::here("Results", "Figs"),
  filename = "N_vs_Dist.eps",
  width = 7,
  height = 6
)


##  Distance vs K ####
data.table::fread(here::here("data", "raw", "Baseline_cleaned and combined2.csv")) %>%
  na.omit() -> Potassium

# Visualize
# ggplot(data = Potassium, aes(Distance, K, label = FarmID)) +
#   geom_point(position = "jitter") +
#   # geom_text(aes(label = FarmID),hjurst = 0,vjust=0)
#   theme_bw()
#
# drop some points
Potassium <- Potassium[K < 3] %>% .[Distance < 60]

# visualize g6---------
g6 <- ggplot(data = Potassium, aes(Distance, K, label = FarmID)) +
  geom_point(position = "jitter") +
  labs(x = "Distance from homestead (m)", y = "Soil Potassium \n (cmol(+)/kg soil)") +
  ggtitle("D") +
  theme_bw(base_size = 12, base_family = "TT Courier New") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = axis_text,
        axis.text.y = axis_text,
        title = axis_text)
ggsave(
  path = here::here("Results", "Figs"),
  filename = "K_vs_Dist.eps",
  width = 7,
  height = 6
)


##  K vs weevil ####
data.table::fread(here::here("data", "raw", "Baseline_cleaned and combined2.csv")) %>%
  na.omit() -> dat
dat1 <- dat[XT > 0] %>%
  .[XT < 10] %>%
  .[K < 3] %>%
  .[FarmID != "Barenwa Burandina"] %>%
  .[FarmID != "Mangarina Elvaida"] %>%
  .[FarmID != "Mehangye Henry"] %>%
  .[FarmID != "Manyirwehi Bosco"] %>%
  .[FarmID != "Kalyegira Olive"] %>%
  .[FarmID != "Nakalanzi cissy"]


# # Visualize
# ggplot(data = dat1, aes(K, XT, label = FarmID)) +
#   geom_point(position = "jitter") +
#   # geom_text(aes(label = FarmID),hjurst = 0,vjust=0)+
#   # ggrepel::geom_text_repel()
#   theme_bw()
# #
#     ggplot(data= dat, aes(K,XI,label = FarmID)) +
#         geom_point(position = "jitter")+
#         # geom_text(aes(label = FarmID),hjurst = 0,vjust=0)
#         theme_bw()

# get bpts
dat1 %>%
  setorder(.,-K) %>% # Order Ascending
  .[, Bpts := get_bpts(XT)] %>% # calculating Boundary points
  unique(., by = "Bpts") %>% # retain unique values
  na.omit() -> Bpts_F # drop NAs

# split the data into top 10%
p <- max(dat1$XT) * 0.9
dat1_10 <- dat1 %>% dplyr::filter(XT > p)
ymax <- mean(dat1_10$XT)


# fit nls
m10 <-
  nls(Bpts ~ ymax / (1 + k * exp(-R * K)),
      start = list(k = 1, R = 0.002),
      data = Bpts_F)

z <- tidy(m10)
z
m11 <- nls(
  Bpts ~ (max(Bpts_F$Bpts)) / (1 + k * exp(-R * K)),
  start = list(k = z$estimate[1], R = z$estimate[2]),
  data = Bpts_F
)
z <- tidy(m11)
z
fwrite(z, here::here("Results", "Tables", "XT_vs_K.csv"))
formula(m11)

Bpts_F[, Pred := predict(m11)]

# visualize g4---------
g4 <- ggplot(Bpts_F) +
  geom_point(aes(K, Bpts), position = "jitter", colour = "red") +
  geom_smooth(aes(K, Bpts), colour = "red", se = F) +
  geom_point(data = dat1, aes(K, XT), position = "jitter") +
  expand_limits(x = 0, y = 0) +
  labs(x = "Potassium (cmol (+)/kg soil)", y = "Weevil damage (%)") +
  ggtitle("E") +
  theme_bw(base_size = 12, base_family = "TT Courier New") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.text.x = axis_text,
        axis.text.y = axis_text,
        title = axis_text)
ggsave(
  path = here::here("Results", "Figs"),
  filename = "XT_vs_K.eps",
  width = 7,
  height = 6
)

# Combine the graphs---------
Grand <- (g1 | g2) / (g3 | g6) / (g4 | ggg)
Grand
graph2doc(
  file = paste0(here::here("Results", "Figs"), "/", "Grand.docx"),
  width = 6,
  height = 8
)

windows()
Grand
citation()
                                     

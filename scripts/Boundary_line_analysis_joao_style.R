# Additional analysis to consider
## weevil damage vs Nitrogen

##  Packages and custom functions used
library(tidyverse) # the world of tidy data wrangling
library(data.table) # data structure
library(MASS) # modelling negative binomial
library(broom) # summaring model info into tidy tibbles
library(msme) # dispersion parameter for poisson family
library(TRADER) # boundaryline analysis

# fetch my custom functions
source(here::here("scripts", "Custom functions.R"))

##  Data preparation    ####

# yield data
DT <- fread(here::here("data", "raw", "Baseline_Yield data.csv"))

# Choose variables of interest
yield <- DT[, c(
  "ID", "Region", "Farm ID", "Plot ID",
  "Dist (m)", "Mat No.", "Plant No.", "Girth_0 (cm)",
  "Girth_100 (cm)", "# Hands", "# Fingers"
)]
yield <- na.omit(yield)
setnames(
  yield, c("Dist (m)", "Girth_0 (cm)", "Girth_100 (cm)", "# Hands", "# Fingers"),
  c("Dist", "Girth0", "Girth100", "Hands", "Fingers")
)
get_yield(yield)


##  Distance vs yield   ####
# visualize first
ggplot(data = yield, aes(Dist, BWT, colour = Region)) +
  geom_point(position = "jitter") +
  theme_bw()

# Get boundary points
yield %>%
  .[BWT < 55] %>%
  .[Dist < 100] %>% # drop outliers
  .[, Rel_BWT := (BWT / (max(.[, BWT])))] -> yield # relative Bunch Weight

# split the data into top 10% yielders
yield_10 <- yield %>% dplyr::filter(Rel_BWT > 0.9)
ymax <- mean(yield_10$Rel_BWT)

# visualize again
ggplot(data = yield, aes(Dist, Rel_BWT, colour = Region)) +
  geom_point(position = "jitter") +
  theme_bw()

yield %>%
  setorder(-Dist) %>% # sorting
  .[, Bpts := get_bpts(Rel_BWT)] %>% # picking bpts
  unique(., by = "Bpts") %>% # picking unique ones
  na.omit() -> Bpts_F # Droping NAs $ renaming

# fit nls model to Boundary points
m1 <- nls(Bpts ~ ymax / (1 + k * exp(-R * Dist)), start = list(k = 1, R = 8.35e-5), Bpts_F)
z <- tidy(m1)
m2 <- nls(Bpts ~ 1 / (1 + k * exp(-R * Dist)),
          start = list(k = z$estimate[1], R = z$estimate[2]),
          data = Bpts_F)
formula(m2)
tidy_m1 <- tidy(m2)
fwrite(tidy_m1, here::here("Results", "Tables", "Dist_vs_yield.csv"))

# use the model to predict
Bpts_F[, pred := predict(m2)]

# Visualize
g1 <- ggplot(Bpts_F) +
  geom_point(aes(Dist, Bpts)) +
  geom_smooth(aes(Dist, pred), se = F, colour = "red") +
  geom_point(data = yield, aes(Dist, Rel_BWT), colour = "black") +
  expand_limits(x = 0, y = 0) +
  labs(x = "Distance from homestead (m)", y = "Relative Bunch weight") +
  ggtitle("A") +
  theme_bw()
ggsave(path = here::here("Results", "Figs"), filename = "BWT_vs_Dist.png", width = 7, height = 6)


##  Distance vs weevil  ####

# weevil data
weevil <- fread(here::here("data", "raw", "Baseline_weevil data.csv"))
setnames(weevil, c(
  "Dist (m)", "UXI (% Damage)", "UXO (% Damage)",
  "LLI (% Damage)", "LXO (% Damage)"
),
c("Distw", "UI_damage", "UO_damage", "LI_damage", "LO_damage"),
skip_absent = T
)
# Choose the variables of interest
weevil <- weevil[, c("ID", "Region", "Distw", "UI_damage", "UO_damage", "LI_damage", "LO_damage")]
weevil <- na.omit(weevil)
get_weevil_damage(weevil)

# Visualize
ggplot(data = weevil, aes(Distw, XI, colour = Region, label = ID)) +
  geom_point(position = "jitter") +
  geom_text(aes(label = ID), hjust = 0, vjust = 0)+
  theme_bw()

# ggplot(data= weevil[XT<30], aes(Distw ,XT)) +
#     geom_point(position = "jitter")+
#     theme_bw()
# drop outliers
weevil <- weevil[XI<30] %>% .[ID!=237] %>% .[ID!=1] %>% .[Distw<150] %>% 
# weevil <- weevil[XI < 40] %>%
  # .[ID != 237] %>%
  # .[ID != 1] %>%
  .[Distw < 150]

# split the data into top 10% 
p <- max(weevil$XI)*0.9
Weevil_10 <- weevil %>% dplyr::filter(XI > p)
ymax <- mean(Weevil_10$XI)

# get bpts
weevil %>%
  setorder(., -Distw) %>% # Order Ascending
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
          data = Bpts_F
)
tidy_m4 <- tidy(m4)
fwrite(tidy_m4, here::here("Results", "Tables", "XT_vs_dist.csv"))
formula(m4)

Bpts_F[, Pred := predict(m4)]

# plot and save
g2 <- ggplot(Bpts_F) +
  geom_point(aes(Distw, Bpts)) +
  geom_smooth(aes(Distw, Pred), se = F, colour = "red") +
  geom_point(data = weevil, aes(Distw, XI), position = "jitter") +
  expand_limits(x = 0, y = 0) +
  labs(x = "Distance from homestead (m)", y = "weevil damage(%)") +
  ggtitle("B") +
  theme_bw()
ggsave(path = here::here("Results", "Figs"), filename = "damage_vs_Dist.png", width = 7, height = 6)


##  Weevil vs N ####
data.table::fread(here::here("data", "raw", "Baseline_cleaned and combined2.csv")) %>%
  na.omit() -> Nitrogen1

# drop some points
Nitrogen1 <- Nitrogen1[N < 0.3] %>% .[XT < 9] %>% .[FarmID != "Mwesigye Alkarito"] 

# Visualize
ggplot(data = Nitrogen1, aes(N,XT)) +
  geom_point(position = "jitter") +
  labs(x="Soil Nitrogen (%)", y= "Weevil damage(%)",
       title="Soil nitrogen against weevil damage in Bananas")+
  theme_bw(base_size = 16)+
  expand_limits(x = 0, y = 0) +
  theme_bw(base_size = 12)


# split the data into top 10% 
p <- max(Nitrogen1$XT)*0.9
Nitrogen1_10 <- Nitrogen1 %>% dplyr::filter(XT > p)
ymax <- mean(Nitrogen1_10$XT)

# get bpts
Nitrogen1 %>%
  setorder(., N) %>%             # Order Ascending
  .[, Bpts := get_bpts(XT)] %>%   # calculating Boundary points
  unique(., by = "Bpts") %>%      # retain unique values
  na.omit(.) -> Bpts_F            # drop NAs

# fit nls to the boundary points
mmm <- nls(Bpts~ymax/(1+k*exp(-R*N)),
          start=list(k=1,R=-0.05693),
          data = Bpts_F)
tidy_mmmm <- tidy(mmm)

mmm2 <- nls(Bpts~ymax/(1+k*exp(-R*N)),
           start=list(k= tidy_mmmm$estimate[1], R= tidy_mmmm$estimate[2]),
           data = Bpts_F)
tidy_mmmm2 <- tidy(mmm2)

# predict usig m5
Bpts_F[,Pred := predict(mmm2)]

# plot and save
ggg <- ggplot(Bpts_F) +
  geom_point(aes(N,XT)) +
  geom_smooth(aes(N, Pred), se = F, colour = "red") +
  geom_point(data = Nitrogen1, aes(N,XT), position = "jitter") +
  expand_limits(x = 0, y = 0) +
  labs(x = "Soil Nitrogen (%)", y = "weevil damage(%)") +
  ggtitle("E") +
  theme_bw(base_size = 12)
ggsave(path = here::here("Results", "Figs"), filename = "Nitrogen vs XT.png", width = 7, height = 6)


##  Distance vs N ####
data.table::fread(here::here("data", "raw", "Baseline_cleaned and combined2.csv")) %>%
  na.omit() -> Nitrogen

# Visualize
ggplot(data = Nitrogen, aes(Distance, N, label = FarmID)) +
  geom_point(position = "jitter") +
  # geom_text(aes(label = FarmID),hjurst = 0,vjust=0)
  theme_bw()
# drop some points
Nitrogen <- Nitrogen[N < 0.25] %>% .[Distance > 5]

# split the data into top 10% 
p <- max(Nitrogen$N)*0.9
Nitrogen_10 <- Nitrogen %>% dplyr::filter(N > p)
ymax <- mean(Nitrogen_10$XI)

# get bpts
Nitrogen %>%
  setorder(., -Distance) %>% # Order Ascending
  .[, Bpts := get_bpts(N)] %>% # calculating Boundary points
  unique(., by = "Bpts") %>% # retain unique values
  na.omit(.) -> Bpts_F # drop NAs

# fit nls to the boundary points does not fit well
m5 <- nls(Bpts~ymax/(1+k*exp(-R*Distance)),
          start=list(k=0.02024,R=-0.05693),
          data = Bpts_F)
summary(m5)

# predict usig m5
Bpts_F[,Pred := predict(m5)]

# fit a polynomial
p <- lm(Bpts ~ poly(Distance, 2, raw = T), data = Bpts_F)
N_vs_dist <- tidy(p)
fwrite(N_vs_dist, here::here("Results", "Tables", "N_vs_dist.csv"))
p_augment <- augment(p)
setnames(p_augment, "poly.Distance..2..raw...T.", "poly")

# plot the polynomial fit onto the general scatter
g3 <- ggplot(Bpts_F) +
  geom_point(aes(Distance, Bpts), position = "jitter") +
  geom_smooth(data = p_augment, aes(poly[, 1], .fitted), colour = "red") +
  geom_point(data = Nitrogen, aes(Distance, N), position = "jitter") +
  expand_limits(x = 7, y = 0.1) +
  labs(x = "Distance from homestead (m)", y = "Soil Nitrogen") +
  ggtitle("C") +
  theme_bw()
ggsave(path = here::here("Results", "Figs"), filename = "N_vs_Dist.png", width = 7, height = 6)


##  Distance vs K ####
data.table::fread(here::here("data", "raw", "Baseline_cleaned and combined2.csv")) %>%
  na.omit() -> Potassium

# Visualize
ggplot(data = Potassium, aes(Distance, K, label = FarmID)) +
  geom_point(position = "jitter") +
  # geom_text(aes(label = FarmID),hjurst = 0,vjust=0)
  theme_bw()

# drop some points
Potassium <- Potassium[K < 3] %>% .[Distance < 60]

# plot
g6 <- ggplot(data = Potassium, aes(Distance, K, label = FarmID)) +
  geom_point(position = "jitter") +
  # geom_smooth(aes(Distance, K), method = "lm", colour = "red")
  # geom_text(aes(label = FarmID),hjurst = 0,vjust=0)
  labs(x = "Distance from homestead (m)", y = "Soil Potassium") +
  ggtitle("F") +
  theme_bw()
ggsave(path = here::here("Results", "Figs"), filename = "K_vs_Dist.png", width = 7, height = 6)


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


# Visualize
ggplot(data = dat1, aes(K, XT, label = FarmID)) +
  geom_point(position = "jitter") +
  # geom_text(aes(label = FarmID),hjurst = 0,vjust=0)+
  # ggrepel::geom_text_repel()
  theme_bw()
#
#     ggplot(data= dat, aes(K,XI,label = FarmID)) +
#         geom_point(position = "jitter")+
#         # geom_text(aes(label = FarmID),hjurst = 0,vjust=0)
#         theme_bw()

# get bpts
dat1 %>%
  setorder(., -K) %>% # Order Ascending
  .[, Bpts := get_bpts(XT)] %>% # calculating Boundary points
  unique(., by = "Bpts") %>% # retain unique values
  na.omit() -> Bpts_F # drop NAs

# split the data into top 10% 
p <- max(dat1$XT)*0.9
dat1_10 <- dat1 %>% dplyr::filter(XT > p)
ymax <- mean(dat1_10$XT)


# fit nls
m10 <- nls(Bpts ~ ymax/ (1 + k * exp(-R * K)), start = list(k = 1, R = 0.002), data = Bpts_F)

z <- tidy(m10)
z
m11 <- nls(Bpts ~ (max(Bpts_F$Bpts)) / (1 + k * exp(-R * K)),
           start = list(k = z$estimate[1], R = z$estimate[2]),
           data = Bpts_F)
z <- tidy(m11)
z
fwrite(z, here::here("Results", "Tables", "XT_vs_K.csv"))
formula(m11)

Bpts_F[, Pred := predict(m11)]

g4 <- ggplot(Bpts_F) +
  geom_point(aes(K, Bpts), position = "jitter", colour = "red") +
  geom_smooth(aes(K, Bpts), colour = "red", se = F) +
  geom_point(data = dat1, aes(K, XT), position = "jitter") +
  expand_limits(x = 0, y = 0) +
  labs(x = "Potassium", y = "Weevil damaage") +
  ggtitle("D") +
  theme_bw()
ggsave(path = here::here("Results", "Figs"), filename = "XT_vs_K.png", width = 7, height = 6)


# fit line
p2 <- lm(Bpts ~ poly(K, 2, raw = T), data = Bpts_F)
# p2 <- lm(Bpts~K, data = Bpts_F)
K_vs_XI <- tidy(p2)
K_vs_XI
fwrite(K_vs_XI, here::here("Results", "Tables", "XT_vs_K.csv"))
p2_augment <- augment(p2)
p2_augment
setnames(p2_augment, "poly.K..2..raw...T.", "poly")

# plot the polynomial fit onto the general scatter
ggplot(Bpts_F) +
  geom_point(aes(K, Bpts), position = "jitter", colour = "red") +
  geom_smooth(data = p2_augment, aes(poly[, 1], .fitted), colour = "red") +
  geom_point(data = dat1, aes(K, XT), position = "jitter") +
  expand_limits(x = 0, y = 0) +
  labs(x = "Soil Potassium", y = "Weevil damaage") +
  ggtitle("Potassium suppresses Weevil damage") +
  theme_bw()
ggsave(path = here::here("Results", "Figs"), filename = "XT_vs_K_poly.png", width = 7, height = 6)


# weevil with yield
dat2 <- dat[XT < 8] %>% .[BWT < 31] %>% .[, Rel_BWT := (BWT / (max(.[, BWT])))] # relative Bunch Weight
ggplot(dat2) +
  geom_point(aes(XT, Rel_BWT), position = "jitter", colour = "black") +
  theme_bw()

# get bpts
dat2 %>%
  setorder(., -XT) %>% # Order Ascending
  .[, Bpts := get_bpts(Rel_BWT)] %>% # calculating Boundary points
  unique(., by = "Bpts") %>% # retain unique values
  na.omit() -> Bpts_F # drop NAs

# split the data into top 10% 
dat2_10 <- dat2 %>% dplyr::filter(Rel_BWT > 0.9)
ymax <- mean(dat2_10$Rel_BWT)

r <- Bpts_F[Bpts > 28]
s <- Bpts_F[Bpts < 25]
DT <- rbind(r, s)

# fit nls
mm3 <- nls(Bpts ~ ymax / (1 + k * exp(-R * XT)),
           start = list(k = 0.7, R = 0.01),
           data = DT)
summary(mm3)
z <- tidy(mm3)
z
mm4 <- nls(Bpts ~ (max(Bpts_F$Bpts)) / (1 + k * exp(-R * XT)),
           start = list(k = z$estimate[1], R = z$estimate[2]),
           data = DT
)
z <- tidy(mm4)
z
fwrite(z, here::here("Results", "Tables", "XT_vs_K.csv"))
formula(mm4)

DT[, Pred := predict(mm4)]

g5 <- ggplot(DT) +
  geom_point(aes(XT, Bpts), position = "jitter", colour = "red") +
  geom_smooth(aes(XT, Pred), colour = "red") +
  geom_point(data = dat2, aes(XT, Rel_BWT), position = "jitter") +
  expand_limits(x = 0, y = 0) +
  labs(y = "Relative Bunch Weigth", X = "Total weevil damaage") +
  ggtitle("E") +
  theme_bw()
ggsave(path = here::here("Results", "Figs"), filename = "BWT_vs_XT.png", width = 7, height = 6)


# fit line
p3 <- lm(Bpts ~ poly(XT, 2, raw = T), data = DT)
XT_vs_BWT <- tidy(p3)
XT_vs_BWT
fwrite(XT_vs_BWT, here::here("Results", "Tables", "XT_vs_BWT.csv"))
p3_augment <- augment(p3)
p3_augment
setnames(p3_augment, "poly.XT..2..raw...T.", "poly")

# plot the polynomial fit onto the general scatter
ggplot(DT) +
  geom_point(aes(XT, Bpts), position = "jitter", colour = "red") +
  geom_smooth(data = p3_augment, aes(poly[, 1], .fitted), colour = "red") +
  geom_point(data = dat2, aes(XT, Rel_BWT), position = "jitter") +
  expand_limits(x = 0, y = 0) +
  labs(y = "Fresh Bunch Weigth", x = "Weevil damaage") +
  ggtitle("E") +
  theme_bw()
ggsave(path = here::here("Results", "Figs"), filename = "BWT_vs_XT_poly.png", width = 7, height = 6)


#Load Data
getwd()
setwd("C:/Users/watts/OneDrive - East Carolina University/Documents/Thesis/CSV files for R/")
HMALL <- read.csv("C:/Users/watts/OneDrive - East Carolina University/Thesis/CSV files for R/HeavyMetalsALLBDL.csv")
Points <- read.csv("C:/Users/watts/OneDrive - East Carolina University/Thesis/CSV files for R/SamplingPointsAll.csv")
#Load packages
install.packages("sf")
install.packages("sfdep")
install.packages("spdep")
library(sf)
library(sfdep)
library(spdep)
library(dplyr)
library(tidyr)
library(ggplot2)

HM_sf <- st_as_sf(
  HMALL,
  coords = c("Longitude", "Latitude"),
  crs = 4326   # WGS84 lat/long
)

HM_sf$As <- as.numeric(HM_sf$As)
HM_sf$B <- as.numeric(HM_sf$B)
HM_sf$Ca <- as.numeric(HM_sf$Ca)
HM_sf$Cd <- as.numeric(HM_sf$Cd)
HM_sf$Co <- as.numeric(HM_sf$Co)
HM_sf$K <- as.numeric(HM_sf$K)
HM_sf$Mg <- as.numeric(HM_sf$Mg)
HM_sf$Mn <- as.numeric(HM_sf$Mn)
HM_sf$Na <- as.numeric(HM_sf$Na)
HM_sf$Ni <- as.numeric(HM_sf$Ni)
HM_sf$Zn <- as.numeric(HM_sf$Zn)
HM_sf$P <- as.numeric(HM_sf$P)
HM_sf$S <- as.numeric(HM_sf$S)

#visualize data
hist(HMALL$Cr, main = "Distribution of Cr Concentrations", xlab = "Cr Concentrations mg/kg", ylab = "Frequency")
hist(HMALL$Cu, main = "Distribution of Cu Concentrations", xlab = "Cu Concentrations mg/kg", ylab = "Frequency")
hist(HMALL$Fe, main = "Distribution of Fe Concentrations", xlab = "Fe Concentrations mg/kg", ylab = "Frequency")
hist(HMALL$Pb, main = "Distribution of Pb Concentrations", xlab = "Pb Concentrations mg/kg", ylab = "Frequency")
#Visualize Cr Conc across samples
ggplot(HM_sf) +
  geom_sf(aes(color = Cr),size = 3) +
  facet_wrap(~Sample.Year)+
  scale_color_gradient(name = "Cr Conc",
                      low = "white",
                      high = "darkgreen") +
  ggtitle("Cr Conc samples") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 
ggplot(HM_sf) +
  geom_sf(aes(color = Cu),size = 3) +
  facet_wrap(~Sample.Year)+
  scale_color_gradient(name = "Cu Conc",
                       low = "white",
                       high = "darkgreen") +
  ggtitle("Cu Conc samples") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 
ggplot(HM_sf) +
  geom_sf(aes(color = Fe),size = 3) +
  facet_wrap(~Sample.Year)+
  scale_color_gradient(name = "Fe Conc",
                       low = "white",
                       high = "darkgreen") +
  ggtitle("Fe Conc samples") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 
ggplot(HM_sf) +
  geom_sf(aes(color = Pb),size = 3) +
  facet_wrap(~Sample.Year)+
  scale_color_gradient(name = "Pb Conc",
                       low = "white",
                       high = "darkgreen") +
  ggtitle("Pb Conc samples") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 
coords <- st_coordinates(HM_sf)
n <- nrow(HM_sf)

nb <- vector("list", n)

for (i in 1:n) {
  nb[[i]] <- as.integer(c(
    if (i > 1) i - 1 else integer(0),
    if (i < n) i + 1 else integer(0)
  ))
}

class(nb) <- "nb"
attr(nb, "region.id") <- as.character(1:n)
attr(nb, "call") <- match.call()

HM_w_binary <- nb2listw(nb, style = "B") #converts to binary spatial weights
HM_lag_Cr <- lag.listw(HM_w_binary, HM_sf$Cr) #computes spatial lag
HM_lag_Cu <- lag.listw(HM_w_binary, HM_sf$Cu)
HM_lag_Pb <- lag.listw(HM_w_binary, HM_sf$Pb)
HM_lag_Fe <- lag.listw(HM_w_binary, HM_sf$Fe)
HM_lag_Al <- lag.listw(HM_w_binary, HM_sf$Al)
HM_lag_As <- lag.listw(HM_w_binary, HM_sf$As)
HM_lag_B <- lag.listw(HM_w_binary, HM_sf$B)
HM_lag_Ca <- lag.listw(HM_w_binary, HM_sf$Ca)
HM_lag_Cd <- lag.listw(HM_w_binary, HM_sf$Cd)
HM_lag_Co <- lag.listw(HM_w_binary, HM_sf$Co)
HM_lag_K <- lag.listw(HM_w_binary, HM_sf$K)
HM_lag_Mg <- lag.listw(HM_w_binary, HM_sf$Mg)
HM_lag_Mn <- lag.listw(HM_w_binary, HM_sf$Mn)
HM_lag_Na <- lag.listw(HM_w_binary, HM_sf$Na)
HM_lag_Ni <- lag.listw(HM_w_binary, HM_sf$Ni)
HM_lag_P <- lag.listw(HM_w_binary, HM_sf$P)
HM_lag_S <- lag.listw(HM_w_binary, HM_sf$S)
HM_lag_Zn <- lag.listw(HM_w_binary, HM_sf$Zn)

#test for global G stat of Cr, Cu, Fe, Pb
globalG.test(HM_sf$Cr, HM_w_binary)
globalG.test(HM_sf$Cu, HM_w_binary)
globalG.test(HM_sf$Fe, HM_w_binary)
globalG.test(HM_sf$Pb, HM_w_binary)
globalG.test(HM_sf$Al, HM_w_binary)
globalG.test(HM_sf$As, HM_w_binary)
globalG.test(HM_sf$B, HM_w_binary)
globalG.test(HM_sf$Ca, HM_w_binary)
globalG.test(HM_sf$Cd, HM_w_binary)
globalG.test(HM_sf$Co, HM_w_binary)
globalG.test(HM_sf$K, HM_w_binary)
globalG.test(HM_sf$Mg, HM_w_binary)
globalG.test(HM_sf$Mn, HM_w_binary)
globalG.test(HM_sf$Na, HM_w_binary)
globalG.test(HM_sf$Ni, HM_w_binary)
globalG.test(HM_sf$P, HM_w_binary)
globalG.test(HM_sf$S, HM_w_binary)
globalG.test(HM_sf$Zn, HM_w_binary)

#HM_sf <- HM_sf |> select(-nb, -wt)
#test for hotspots
#HM_sf$nb <- nb
wt <- st_weights(nb)
#HM_sf$CrLag <- st_lag(HM_sf$Cr, nb, wt)
# 1. Compute Gi* outside mutate() 
Gi_results <- local_g_perm(HM_sf$Cr, nb, wt, nsim = 999)
HMCr_HS <- HM_sf |>
  mutate(Gi = Gi_results) |>
  unnest(Gi)
Gi_resultsCu <- local_g_perm(HM_sf$Cu, nb, wt, nsim=999)
HMCu_HS <- HM_sf |>
  mutate(GiCu = Gi_resultsCu) |>
  unnest(GiCu) |>
  rename (
    GiCu= gi, 
    pCu = p_folded_sim,
  )
Gi_resultsFe <- local_g_perm(HM_sf$Fe, nb, wt, nsim=999)
HMFe_HS <- HM_sf |>
  mutate(GiFe = Gi_resultsFe) |>
  unnest(GiFe) |>
  rename (
    GiFe= gi, 
    pFe = p_folded_sim,
  )
Gi_resultsPb <- local_g_perm(HM_sf$Pb, nb, wt, nsim=999)
HMPb_HS <- HM_sf |>
  mutate(GiPb = Gi_resultsPb) |>
  unnest(GiPb) |>
  rename (
    GiPb= gi, 
    pPb = p_folded_sim,
  )
Gi_resultsAl <- local_g_perm(HM_sf$Al, nb, wt, nsim=999)
HMAl_HS <- HM_sf |>
  mutate(GiAl = Gi_resultsAl) |>
  unnest(GiAl) |>
  rename (
    GiAl= gi, 
    pAl = p_folded_sim,
  )
Gi_resultsAs <- local_g_perm(HM_sf$As, nb, wt, nsim=999)
HMAs_HS <- HM_sf |>
  mutate(GiAs = Gi_resultsAs) |>
  unnest(GiAs) |>
  rename (
    GiAs= gi, 
    pAs = p_folded_sim,
  )
Gi_resultsB <- local_g_perm(HM_sf$B, nb, wt, nsim=999)
HMB_HS <- HM_sf |>
  mutate(GiB = Gi_resultsB) |>
  unnest(GiB) |>
  rename (
    GiB= gi, 
    pB = p_folded_sim,
  )
Gi_resultsCa <- local_g_perm(HM_sf$Ca, nb, wt, nsim=999)
HMCa_HS <- HM_sf |>
  mutate(GiCa = Gi_resultsCa) |>
  unnest(GiCa) |>
  rename (
    GiCa= gi, 
    pCa = p_folded_sim,
  )
Gi_resultsCd <- local_g_perm(HM_sf$Cd, nb, wt, nsim=999)
HMCd_HS <- HM_sf |>
  mutate(GiCd = Gi_resultsCd) |>
  unnest(GiCd) |>
  rename (
    GiCd= gi, 
    pCd = p_folded_sim,
  )
Gi_resultsCo <- local_g_perm(HM_sf$Co, nb, wt, nsim=999)
HMCo_HS <- HM_sf |>
  mutate(GiCo = Gi_resultsCo) |>
  unnest(GiCo) |>
  rename (
    GiCo= gi, 
    pCo = p_folded_sim,
  )
Gi_resultsK <- local_g_perm(HM_sf$K, nb, wt, nsim=999)
HMK_HS <- HM_sf |>
  mutate(GiK = Gi_resultsK) |>
  unnest(GiK) |>
  rename (
    GiK= gi, 
    pK = p_folded_sim,
  )
Gi_resultsMg <- local_g_perm(HM_sf$Mg, nb, wt, nsim=999)
HMMg_HS <- HM_sf |>
  mutate(GiMg = Gi_resultsMg) |>
  unnest(GiMg) |>
  rename (
    GiMg= gi, 
    pMg = p_folded_sim,
  )
Gi_resultsMn <- local_g_perm(HM_sf$Mn, nb, wt, nsim=999)
HMMn_HS <- HM_sf |>
  mutate(GiMn = Gi_resultsMn) |>
  unnest(GiMn) |>
  rename (
    GiMn= gi, 
    pMn = p_folded_sim,
  )
Gi_resultsNa <- local_g_perm(HM_sf$Na, nb, wt, nsim=999)
HMNa_HS <- HM_sf |>
  mutate(GiNa = Gi_resultsNa) |>
  unnest(GiNa) |>
  rename (
    GiNa= gi, 
    pNa = p_folded_sim,
  )
Gi_resultsNi <- local_g_perm(HM_sf$Ni, nb, wt, nsim=999)
HMNi_HS <- HM_sf |>
  mutate(GiNi = Gi_resultsNi) |>
  unnest(GiNi) |>
  rename (
    GiNi= gi, 
    pNi = p_folded_sim,
  )
Gi_resultsP <- local_g_perm(HM_sf$P, nb, wt, nsim=999)
HMP_HS <- HM_sf |>
  mutate(GiP = Gi_resultsP) |>
  unnest(GiP) |>
  rename (
    GiP= gi, 
    pP = p_folded_sim,
  )
Gi_resultsS <- local_g_perm(HM_sf$S, nb, wt, nsim=999)
HMS_HS <- HM_sf |>
  mutate(GiS = Gi_resultsS) |>
  unnest(GiS) |>
  rename (
    GiS= gi, 
    pS = p_folded_sim,
  )
Gi_resultsZn <- local_g_perm(HM_sf$Zn, nb, wt, nsim=999)
HMZn_HS <- HM_sf |>
  mutate(GiZn = Gi_resultsZn) |>
  unnest(GiZn) |>
  rename (
    GiZn= gi, 
    pZn = p_folded_sim,
  )
#visualization of Gi values
HMCr_HS |>
  ggplot((aes(color = gi))) +
  facet_wrap(~ Sample.Year) +
  geom_sf(size = 1) +
  scale_color_gradient2()+
  ggtitle("Gi* Hotspots for Chromium")
HMCu_HS |>
  ggplot((aes(color = GiCu))) +
  facet_wrap(~ Sample.Year) +
  geom_sf(size = 3) +
  scale_color_gradient2()+
  ggtitle("Gi* Hotspots for Copper")
HMFe_HS |>
  ggplot((aes(color = GiFe))) +
  facet_wrap(~ Sample.Year) +
  geom_sf(size = 3) +
  scale_color_gradient2()+
  ggtitle("Gi* Hotspots for Iron")
HMPb_HS |>
  ggplot((aes(color = GiPb))) +
  facet_wrap(~ Sample.Year) +
  geom_sf(size = 3) +
  scale_color_gradient2()+
  ggtitle("Gi* Hotspots for Lead")
#create new data frame
HMCr_HS |>
  select(gi, p_folded_sim, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
      gi > 0 & p_folded_sim <= 0.05 ~ "Hot",
      gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
      gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
      gi < 0 & p_folded_sim <= 0.05 ~ "Cold",
      gi < 0 & p_folded_sim <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    color = "Hot Spot Classification",
    title = "Cr Concentrations: Accomac and Aowa"
  )
HMCr_HS2 <- HMCr_HS |> filter(Sample.Year == "23-Aug") |> mutate( Classification = case_when( gi > 0 & p_folded_sim <= 0.01 ~ "Very hot", gi > 0 & p_folded_sim <= 0.05 ~ "Hot", gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot", gi < 0 & p_folded_sim <= 0.01 ~ "Very cold", gi < 0 & p_folded_sim <= 0.05 ~ "Cold", gi < 0 & p_folded_sim <= 0.1 ~ "Somewhat cold", TRUE ~ "Insignificant" ), Classification = factor( Classification, levels = c( "Very hot", "Hot", "Somewhat hot", "Insignificant", "Somewhat cold", "Cold", "Very cold" ) ) )
ggplot(HMCr_HS2,aes(color = Classification)) +
  geom_sf(size = 2) +
  scale_fill_brewer(type = "div", palette = 5) +
  coord_sf(xlim = c(-77.27089, -77.268992), ylim = c(38.4700, 38.472), expand = FALSE)+
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Cr Concentrations: Inter-ship Aowa Samples", 
    
  )
HMCu_HS |>
  select(GiCu, pCu, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiCu > 0 & pCu <= 0.01 ~ "Very hot",
      GiCu > 0 & pCu <= 0.05 ~ "Hot",
      GiCu > 0 & pCu <= 0.1 ~ "Somewhat hot",
      GiCu < 0 & pCu <= 0.01 ~ "Very cold",
      GiCu < 0 & pCu <= 0.05 ~ "Cold",
      GiCu < 0 & pCu <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Cu Concentrations: Accomac and Aowa"
  )
HMCu_HS2 <- HMCu_HS |> filter(Sample.Year == "23-Aug") |> mutate( Classification = case_when( GiCu > 0 & pCu <= 0.01 ~ "Very hot", GiCu > 0 & pCu <= 0.05 ~ "Hot", GiCu > 0 & pCu <= 0.1 ~ "Somewhat hot", GiCu < 0 & pCu <= 0.01 ~ "Very cold", GiCu < 0 & pCu <= 0.05 ~ "Cold", GiCu < 0 & pCu <= 0.1 ~ "Somewhat cold", TRUE ~ "Insignificant" ), Classification = factor( Classification, levels = c( "Very hot", "Hot", "Somewhat hot", "Insignificant", "Somewhat cold", "Cold", "Very cold" ) ) )
ggplot(HMCu_HS2,aes(color = Classification)) +
  geom_sf(size = 2) +   
  geom_sf_text(aes(label = Sample.1), size = 2, nudge_y = 0.00005)+
  scale_fill_brewer(type = "div", palette = 5) +
  coord_sf(xlim = c(-77.27089, -77.268992), ylim = c(38.4700, 38.472), expand = FALSE)+
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Cu Concentrations: Inter-ship Aowa Samples", 
    
  )
HMFe_HS |>
  select(GiFe, pFe, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiFe > 0 & pFe <= 0.01 ~ "Very hot",
      GiFe > 0 & pFe <= 0.05 ~ "Hot",
      GiFe > 0 & pFe <= 0.1 ~ "Somewhat hot",
      GiFe < 0 & pFe <= 0.01 ~ "Very cold",
      GiFe < 0 & pFe <= 0.05 ~ "Cold",
      GiFe < 0 & pFe <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Fe Concentrations: Accomac and Aowa"
  )
HMFe_HS2 <- HMFe_HS |> filter(Sample.Year == "23-Aug") |> mutate( Classification = case_when( GiFe > 0 & pFe <= 0.01 ~ "Very hot", GiFe > 0 & pFe <= 0.05 ~ "Hot", GiFe > 0 & pFe <= 0.1 ~ "Somewhat hot", GiFe < 0 & pFe <= 0.01 ~ "Very cold", GiFe < 0 & pFe <= 0.05 ~ "Cold", GiFe < 0 & pFe <= 0.1 ~ "Somewhat cold", TRUE ~ "Insignificant" ), Classification = factor( Classification, levels = c( "Very hot", "Hot", "Somewhat hot", "Insignificant", "Somewhat cold", "Cold", "Very cold" ) ) )
ggplot(HMFe_HS2,aes(color = Classification)) +
  geom_sf(size = 2) +   
  geom_sf_text(aes(label = Sample.1), size = 2, nudge_y = 0.00005)+
  scale_fill_brewer(type = "div", palette = 5) +
  coord_sf(xlim = c(-77.27089, -77.268992), ylim = c(38.4700, 38.472), expand = FALSE)+
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Fe Concentrations: Inter-ship Aowa Samples", 
    
  )
HMPb_HS |>
  select(GiPb, pPb, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiPb > 0 & pPb <= 0.01 ~ "Very hot",
      GiPb > 0 & pPb <= 0.05 ~ "Hot",
      GiPb > 0 & pPb <= 0.1 ~ "Somewhat hot",
      GiPb < 0 & pPb <= 0.01 ~ "Very cold",
      GiPb < 0 & pPb <= 0.05 ~ "Cold",
      GiPb < 0 & pPb <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Pb Concentrations: Accomac and Aowa"
  )
HMPb_HS2 <- HMPb_HS |> filter(Sample.Year == "23-Aug") |> mutate( Classification = case_when( GiPb > 0 & pPb <= 0.01 ~ "Very hot", GiPb > 0 & pPb <= 0.05 ~ "Hot", GiPb > 0 & pPb <= 0.1 ~ "Somewhat hot", GiPb < 0 & pPb <= 0.01 ~ "Very cold", GiPb < 0 & pPb <= 0.05 ~ "Cold", GiPb < 0 & pPb <= 0.1 ~ "Somewhat cold", TRUE ~ "Insignificant" ), Classification = factor( Classification, levels = c( "Very hot", "Hot", "Somewhat hot", "Insignificant", "Somewhat cold", "Cold", "Very cold" ) ) )
ggplot(HMPb_HS2,aes(color = Classification)) +
  geom_sf(size = 2) +   
  scale_fill_brewer(type = "div", palette = 5) +   
  geom_sf_text(aes(label = Sample.1), size = 2, nudge_y = 0.00005)+
  coord_sf(xlim = c(-77.27089, -77.268992), ylim = c(38.4700, 38.472), expand = FALSE)+
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Pb Concentrations: Inter-ship Aowa Samples", 
    
  )

HMAl_HS |>
  select(GiAl, pAl, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiAl > 0 & pAl <= 0.01 ~ "Very hot",
      GiAl > 0 & pAl <= 0.05 ~ "Hot",
      GiAl > 0 & pAl <= 0.1 ~ "Somewhat hot",
      GiAl < 0 & pAl <= 0.01 ~ "Very cold",
      GiAl < 0 & pAl <= 0.05 ~ "Cold",
      GiAl < 0 & pAl <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Al Concentrations: Accomac and Aowa"
  )
HMAs_HS |>
  select(GiAs, pAs, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiAs > 0 & pAs <= 0.01 ~ "Very hot",
      GiAs > 0 & pAs <= 0.05 ~ "Hot",
      GiAs > 0 & pAs <= 0.1 ~ "Somewhat hot",
      GiAs < 0 & pAs <= 0.01 ~ "Very cold",
      GiAs < 0 & pAs <= 0.05 ~ "Cold",
      GiAs < 0 & pAs <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "As Concentrations: Accomac and Aowa"
  )
HMB_HS |>
  select(GiB, pB, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiB > 0 & pB <= 0.01 ~ "Very hot",
      GiB > 0 & pB <= 0.05 ~ "Hot",
      GiB > 0 & pB <= 0.1 ~ "Somewhat hot",
      GiB < 0 & pB <= 0.01 ~ "Very cold",
      GiB < 0 & pB <= 0.05 ~ "Cold",
      GiB < 0 & pB <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "B Concentrations: Accomac and Aowa"
  )
HMCa_HS |>
  select(GiCa, pCa, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiCa > 0 & pCa <= 0.01 ~ "Very hot",
      GiCa > 0 & pCa <= 0.05 ~ "Hot",
      GiCa > 0 & pCa <= 0.1 ~ "Somewhat hot",
      GiCa < 0 & pCa <= 0.01 ~ "Very cold",
      GiCa < 0 & pCa <= 0.05 ~ "Cold",
      GiCa < 0 & pCa <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Ca Concentrations: Accomac and Aowa"
  )
HMCd_HS |>
  select(GiCd, pCd, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiCd > 0 & pCd <= 0.01 ~ "Very hot",
      GiCd > 0 & pCd <= 0.05 ~ "Hot",
      GiCd > 0 & pCd <= 0.1 ~ "Somewhat hot",
      GiCd < 0 & pCd <= 0.01 ~ "Very cold",
      GiCd < 0 & pCd <= 0.05 ~ "Cold",
      GiCd < 0 & pCd <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Cd Concentrations: Accomac and Aowa"
  )
HMCo_HS |>
  select(GiCo, pCo, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiCo > 0 & pCo <= 0.01 ~ "Very hot",
      GiCo > 0 & pCo <= 0.05 ~ "Hot",
      GiCo > 0 & pCo <= 0.1 ~ "Somewhat hot",
      GiCo < 0 & pCo <= 0.01 ~ "Very cold",
      GiCo < 0 & pCo <= 0.05 ~ "Cold",
      GiCo < 0 & pCo <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Co Concentrations: Accomac and Aowa"
  )
HMK_HS |>
  select(GiK, pK, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiK > 0 & pK <= 0.01 ~ "Very hot",
      GiK > 0 & pK <= 0.05 ~ "Hot",
      GiK > 0 & pK <= 0.1 ~ "Somewhat hot",
      GiK < 0 & pK <= 0.01 ~ "Very cold",
      GiK < 0 & pK <= 0.05 ~ "Cold",
      GiK < 0 & pK <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "K Concentrations: Accomac and Aowa"
  )
HMMg_HS |>
  select(GiMg, pMg, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiMg> 0 & pMg <= 0.01 ~ "Very hot",
      GiMg > 0 & pMg <= 0.05 ~ "Hot",
      GiMg > 0 & pMg <= 0.1 ~ "Somewhat hot",
      GiMg < 0 & pMg <= 0.01 ~ "Very cold",
      GiMg < 0 & pMg <= 0.05 ~ "Cold",
      GiMg < 0 & pMg <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Mg Concentrations: Accomac and Aowa"
  )
HMMn_HS |>
  select(GiMn, pMn, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiMn > 0 & pMn <= 0.01 ~ "Very hot",
      GiMn > 0 & pMn <= 0.05 ~ "Hot",
      GiMn > 0 & pMn <= 0.1 ~ "Somewhat hot",
      GiMn < 0 & pMn <= 0.01 ~ "Very cold",
      GiMn < 0 & pMn <= 0.05 ~ "Cold",
      GiMn < 0 & pMn <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Mn Concentrations: Accomac and Aowa"
  )
HMNa_HS |>
  select(GiNa, pNa, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiNa > 0 & pNa <= 0.01 ~ "Very hot",
      GiNa > 0 & pNa <= 0.05 ~ "Hot",
      GiNa > 0 & pNa <= 0.1 ~ "Somewhat hot",
      GiNa < 0 & pNa <= 0.01 ~ "Very cold",
      GiNa < 0 & pNa <= 0.05 ~ "Cold",
      GiNa < 0 & pNa <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Na Concentrations: Accomac and Aowa"
  )
HMNi_HS |>
  select(GiNi, pNi, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiNi > 0 & pNi <= 0.01 ~ "Very hot",
      GiNi > 0 & pNi <= 0.05 ~ "Hot",
      GiNi > 0 & pNi <= 0.1 ~ "Somewhat hot",
      GiNi < 0 & pNi <= 0.01 ~ "Very cold",
      GiNi < 0 & pNi <= 0.05 ~ "Cold",
      GiNi < 0 & pNi <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Ni Concentrations: Accomac and Aowa"
  )
HMP_HS |>
  select(GiP, pP, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiP > 0 & pP <= 0.01 ~ "Very hot",
      GiP > 0 & pP <= 0.05 ~ "Hot",
      GiP > 0 & pP <= 0.1 ~ "Somewhat hot",
      GiP < 0 & pP <= 0.01 ~ "Very cold",
      GiP < 0 & pP <= 0.05 ~ "Cold",
      GiP < 0 & pP <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "P Concentrations: Accomac and Aowa"
  )
HMS_HS |>
  select(GiS, pS, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiS > 0 & pS <= 0.01 ~ "Very hot",
      GiS > 0 & pS <= 0.05 ~ "Hot",
      GiS > 0 & pS <= 0.1 ~ "Somewhat hot",
      GiS < 0 & pS <= 0.01 ~ "Very cold",
      GiS < 0 & pS <= 0.05 ~ "Cold",
      GiS < 0 & pS <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "S Concentrations: Accomac and Aowa"
  )
HMZn_HS |>
  select(GiZn, pZn, geometry, Sample.Year) |>
  mutate(
    # Add a new column called "classification"
    Classification = case_when(
      # Classify based on the following criteria:
      GiZn > 0 & pZn <= 0.01 ~ "Very hot",
      GiZn > 0 & pZn <= 0.05 ~ "Hot",
      GiZn > 0 & pZn <= 0.1 ~ "Somewhat hot",
      GiZn < 0 & pZn <= 0.01 ~ "Very cold",
      GiZn < 0 & pZn <= 0.05 ~ "Cold",
      GiZn < 0 & pZn <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    Classification = factor(
      Classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(color = Classification)) +
  geom_sf(size = 2) +
  facet_wrap(~ Sample.Year) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Zn Concentrations: Accomac and Aowa"
  )

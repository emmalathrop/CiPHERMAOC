# Load libraries and directories ------------------------------------------
#data manipulation and visualization
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
#statistics
library(lme4)
library(glmmTMB)
library(sjPlot)
library(DHARMa)
library(emmeans)
library(rstatix)
library(nlme)
rm(list=ls()) # clear the workspace

# Plot theme set ----------------------------------------------------------
ts_theme <- theme_classic() +
  theme(plot.title = element_text(size = 26, hjust = 0.5),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.key = element_rect(fill = "transparent"),
        legend.title=element_text(size=20),
        legend.text=element_text(size=16),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks=element_line(color = "black", linewidth=1),
        strip.text.x = element_text(size = 12, color = "black"),
        strip.background = element_rect(color="black", fill="grey83", size=1.5, linetype="solid"))
theme_set(ts_theme)

treatPal <- c("#808080","#377EB8", "#E41A1C") 
offsetPal <- c("#628B48", "#F0C808", "#1C448E")
fracPal <- c("#FFFFFF", "#80808025")

# Load data and calculate density fraction metrics ---------------------------------------------------------------
soil <- read.csv("datasets/soil.csv") %>% select(-X)
stockDf <- read.csv("datasets/stockDf.csv") %>% select(-X)

# Percent Carbon in each fraction -----------------------------------------

#What is the %C in each fraction throughout the profile? 
percCPerFrac <- soil %>% 
  mutate(LFracC = LFracC/10, 
         HFracC = HFracC/10) %>% 
  dplyr::select(year, fence, treatment, midDepth, LFracC, HFracC) %>% 
  tidyr::pivot_longer(cols=c('LFracC', 'HFracC'),
                      names_to='Fraction',
                      values_to='percC') 

percCSummary <- percCPerFrac %>% 
  group_by(Fraction, treatment) %>% 
  summarise(mean = mean(percC, na.rm= TRUE),
            se = plotrix::std.error(percC, na.rm = TRUE))

mod <- lmer(HFracC ~ treatment + (1|fence), data = soil)
tab_model(mod)
emmeans::emmeans(mod, pairwise ~ treatment, adjust = "none")

# Analysis of fraction by depth (non-ash corrected) ------------------------------------------------------
soilLong <- soil %>%
  select(-c("depth0", "depth1", "block", "depth.cat", "soil.stock", "coreNum")) %>%
  pivot_longer(cols = -c(year, fence, plot, midDepth, treatment),
               names_to = "variable",
               values_to = "value")

soilSumTreat <- soilLong %>%
  group_by(variable, midDepth, treatment) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE)/sqrt(sum(!is.na(value))),
            .groups = "drop")

propSum <- stockDf %>%
  select(-c("depth0", "depth1", "block", "depth.cat", "soil.stock", "coreNum")) %>%
  pivot_longer(cols = -c(year, fence, plot, midDepth, treatment),
               names_to = "variable",
               values_to = "value") %>%
  filter(variable %in% c("propLC", "propHC")) %>%
  group_by(variable, midDepth, treatment) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE)/sqrt(sum(!is.na(value))),
            .groups = "drop")

## Depth plots for manuscript ----------------------------------------------
a <- ggplot(data = filter(soilSumTreat, variable %in% c("HFracC") & (midDepth > 20 & midDepth < 80)),
            aes(x = midDepth, y = mean/10, color = treatment)) +
  geom_point(alpha = 0.9) +
  geom_line(linewidth = 1, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean/10 - se/10, ymax = mean/10 + se/10), linewidth = 1, alpha = 0.9, width = 1.5) +
  scale_x_reverse() + coord_flip()+
  scale_color_manual(values = treatPal, 
                     name = "Treatment",
                     labels = c("2009", "2022: Ambient", "2022: Warming"))+
  facet_wrap(~variable, scales = "free", nrow = 1, labeller = as_labeller(c(`HFracC` = "%C in MAOC")))+
  labs(y = "%C of MAOC", x = "Depth (cm)")+
  theme(legend.position = "top",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15))
a
#ggsave(paste0(manDir, "figures/suppHFracCDepth", Sys.Date(), ".png"), width = 15, height = 20, units = "cm", dpi = 300)

# Fraction stocks with equivalent ash normalization -----------------------
#Create dataframe with soil samples by core, and calculate the cumulative sum of properties for stock calculations and equivalent ash method
core <- stockDf %>%            
  group_by(year, block, fence, plot, treatment) %>% 
  mutate(coreNum = cur_group_id(), 
         cu.soil.stock = cumsum(soil.stock),
         cu.ash.stock = cumsum(ash.stock),
         cu.C.stock = cumsum(C.stock),
         cu.N.stock = cumsum(N.stock), 
         cu.LC.stock = cumsum(stockLC), 
         cu.HC.stock = cumsum(stockHC), 
         cu.LN.stock = cumsum(stockLN),
         cu.HN.stock = cumsum(stockHN),
         bulkAsh = bulk.density*ash) %>% 
  ungroup()

#calculate the cumulative ash content that corresponds to the minimum ash content to reach desired depths
#from plaza et al. 209 Nature Geoscience and Lathrop et al. 2025 Global Change Bio
depths <- c(35,55,75)

ash.indices <- c()
placemark <- 1
cores2009 <- subset(core, year == 2009) %>% 
  filter(!is.na(ash.stock))
for(i in depths){
  ash.indices[placemark] <- min(cores2009[which(cores2009$depth1 == i), "cu.ash.stock"], na.rm = TRUE)
  placemark <- placemark + 1
}

#function to calculate sum of fraction stocks
sum.soil.prop <- function(depth0, ash.stock, cu.ash.stock, soil.prop.stock, cu.soil.prop.stock, ash.index){
  index <- first(which(cu.ash.stock >= ash.index)) #find which depth reaches the cumulative ash index val
  
  if(length(index) == 0 | length(which(!is.na(index))) == 0){ #conditional when the sample doesn't go deep enough to accumulate enough ash
    soil.c.stock <- NA
  } else if(depth0[index]==0){ #if the ash index is reached in the first sample in the core
    perc.layer <- (ash.stock[index] - (cu.ash.stock[index] - ash.index))/ash.stock[index]
    soil.c.stock <- cu.soil.prop.stock[index]*perc.layer
  } else {
    perc.layer <- (ash.stock[index] - (cu.ash.stock[index] - ash.index))/ash.stock[index]
    soil.c.stock <- cu.soil.prop.stock[index-1] + (soil.prop.stock[index]*perc.layer)
  }
  
  return(soil.c.stock)
  
}

bulkStockCalcCore <- core %>%
  group_by(year, block, fence, plot, treatment) %>%
  dplyr::summarize(stock.o = sum.soil.prop(depth0, ash.stock, cu.ash.stock, C.stock, cu.C.stock, ash.index = ash.indices[1]),
                   stock.m = sum.soil.prop(depth0, ash.stock, cu.ash.stock, C.stock, cu.C.stock, ash.index = ash.indices[2]),
                   stock.dm = sum.soil.prop(depth0, ash.stock, cu.ash.stock, C.stock, cu.C.stock, ash.index = ash.indices[3]),
                   Nstock.o = sum.soil.prop(depth0, ash.stock, cu.ash.stock, N.stock, cu.N.stock, ash.index = ash.indices[1]),
                   Nstock.m = sum.soil.prop(depth0, ash.stock, cu.ash.stock, N.stock, cu.N.stock, ash.index = ash.indices[2]),
                   Nstock.dm = sum.soil.prop(depth0, ash.stock, cu.ash.stock, N.stock, cu.N.stock, ash.index = ash.indices[3])) %>%
  mutate(fraction = "Bulk") 

#Calculate H and L fraction stocks
LFracStockCalcCore <- core %>% 
  group_by(year, block, fence, plot, treatment) %>% 
  summarise(stock.o = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockLC, cu.LC.stock, ash.index = ash.indices[1]),
            stock.m = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockLC, cu.LC.stock, ash.index = ash.indices[2]),
            stock.dm = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockLC, cu.LC.stock, ash.index = ash.indices[3]),
            Nstock.o = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockLN, cu.LN.stock, ash.index = ash.indices[1]),
            Nstock.m = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockLN, cu.LN.stock, ash.index = ash.indices[2]),
            Nstock.dm = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockLN, cu.LN.stock, ash.index = ash.indices[3])) %>% 
  mutate(fraction = "Light",
         #make NA for the single core that was not fractioned
         across(starts_with("stock"), ~ if_else(year == 2022 & (fence == 3 & plot == 4), NA_real_, .)),
         across(starts_with("Nstock"), ~ if_else(year == 2022 & (fence == 3 & plot == 4), NA_real_, .)))

HFracStockCalcCore <- core %>% 
  group_by(year, block, fence, plot, treatment) %>% 
  dplyr::summarize(stock.o = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockHC, cu.HC.stock, ash.index = ash.indices[1]),
                   stock.m = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockHC, cu.HC.stock, ash.index = ash.indices[2]),
                   stock.dm = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockHC, cu.HC.stock, ash.index = ash.indices[3]),
                   Nstock.o = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockHN, cu.HN.stock, ash.index = ash.indices[1]),
                   Nstock.m = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockHN, cu.HN.stock, ash.index = ash.indices[2]),
                   Nstock.dm = sum.soil.prop(depth0, ash.stock, cu.ash.stock, stockHN, cu.HN.stock, ash.index = ash.indices[3])) %>%
  mutate(fraction = "Heavy",
         #make NA for the single core that was not fractioned
         across(starts_with("stock"), ~ if_else(year == 2022 & (fence == 3 & plot == 4), NA_real_, .)),
         across(starts_with("Nstock"), ~ if_else(year == 2022 & (fence == 3 & plot == 4), NA_real_, .)))

FracStockCalcCoreNAs <- rbind(LFracStockCalcCore, HFracStockCalcCore, bulkStockCalcCore) %>%
  ungroup() %>% 
  rowwise() %>% 
  mutate(stock.dm = stock.dm-stock.m,
         stock.m = stock.m-stock.o,
         Nstock.dm = Nstock.dm-Nstock.m,
         Nstock.m = Nstock.m-Nstock.o) %>% 
  tidyr::pivot_longer(cols  = -c( "year", "block", "plot", "fence", "treatment", "fraction"), names_to = c(".value", "depth"),
                      names_sep = "\\.")

#gap fill missing samples
LFracCMod <- lm(stock ~ depth * treatment, data = subset(FracStockCalcCoreNAs, fraction == "Light"))
HFracCMod <- lm(stock ~ depth * treatment, data = subset(FracStockCalcCoreNAs, fraction == "Heavy"))
bulkCMod <- lm(stock ~ depth * treatment, data = subset(FracStockCalcCoreNAs, fraction == "Bulk"))

LFracNMod <- lm(Nstock ~ depth * treatment, subset(FracStockCalcCoreNAs, fraction == "Light"))
HFracNMod <- lm(Nstock ~ depth * treatment, subset(FracStockCalcCoreNAs, fraction == "Heavy"))
bulkNMod <- lm(Nstock ~ depth * treatment, data = subset(FracStockCalcCoreNAs, fraction == "Bulk"))

FracStockCalcCore <- FracStockCalcCoreNAs %>% 
  mutate(stock = if_else(fraction == "Light" & is.na(stock), predict(LFracCMod, newdata = FracStockCalcCoreNAs), stock),
         stock = if_else(fraction == "Heavy" & is.na(stock), predict(HFracCMod, newdata = FracStockCalcCoreNAs), stock),
         stock = if_else(fraction == "Bulk" & is.na(stock), predict(bulkCMod, newdata = FracStockCalcCoreNAs), stock),
         Nstock = if_else(fraction == "Light" & is.na(Nstock), predict(LFracNMod, newdata = FracStockCalcCoreNAs), Nstock),
         Nstock = if_else(fraction == "Heavy" & is.na(Nstock), predict(HFracNMod, newdata = FracStockCalcCoreNAs), Nstock),
         Nstock = if_else(fraction == "Bulk" & is.na(Nstock), predict(bulkNMod, newdata = FracStockCalcCoreNAs), Nstock)) %>% 
  tidyr::pivot_wider(names_from = c("depth"), 
                     values_from = -c("year", "block", "plot", "fence", "treatment", "depth", "fraction")) %>% 
  rowwise() %>% 
  mutate(stock_mineral = stock_m + stock_dm,
         stock_total = sum(stock_o, stock_mineral),
         Nstock_mineral = Nstock_m + Nstock_dm,
         Nstock_total = sum(Nstock_o, Nstock_mineral)) %>% 
  tidyr::pivot_longer(cols  = -c( "year", "block", "plot", "fence", "treatment", "fraction"), names_to = c(".value", "depth"),
                      names_sep = "\\_") %>% 
  tidyr::pivot_wider(names_from = fraction, values_from = c("stock", "Nstock")) %>% 
  dplyr::rename(stockLC = stock_Light, stockHC = stock_Heavy, stockBulk = stock_Bulk, 
                stockLN = Nstock_Light, stockHN = Nstock_Heavy, stockBulkN = Nstock_Bulk)

#wide dataframe with ash norm'd C and N stocks
ashNormStock <- FracStockCalcCore %>% 
  mutate(treatment = factor(treatment, 
                            levels = c("i", "c", "w"), 
                            labels = c("2009", "2022: Control", "2022: Warming"))) %>% 
  dplyr::rename(Cstock.L = stockLC,
                Cstock.H = stockHC,
                Nstock.L = stockLN,
                Nstock.H = stockHN,
                Cstock.bulk = stockBulk,
                Nstock.bulk = stockBulkN) %>% 
  mutate(CN.L = Cstock.L/Nstock.L,
         CN.H = Cstock.H/Nstock.H,
         propLight = ifelse(Cstock.L > 0, Cstock.L/Cstock.bulk, NA),
         propHeavy = ifelse(Cstock.H > 0, Cstock.H/Cstock.bulk, NA))

#long data frame for plotting
ashNormStockLong <- ashNormStock %>% 
  select(-Cstock.bulk, -Nstock.bulk) %>% 
  tidyr::pivot_longer(cols  = -c( "year", "block", "plot", "fence", "treatment", "depth","propLight", "propHeavy"), names_to = c(".value", "fraction"),
                      names_sep = "\\.")  %>% 
  mutate(treatment = factor(treatment, 
                            levels = c("2009", "2022: Control", "2022: Warming"), 
                            labels = c("i", "c", "w")))

ashNormStockLongWBulk <- ashNormStock %>% 
  tidyr::pivot_longer(cols  = -c( "year", "block", "plot", "fence", "treatment", "depth","propLight", "propHeavy"), names_to = c(".value", "fraction"),
                      names_sep = "\\.")  %>% 
  mutate(treatment = factor(treatment, 
                            levels = c("2009", "2022: Control", "2022: Warming"), 
                            labels = c("i", "c", "w")))

#write.csv(FracStockCalcCore, file = paste0(soilDir, "/stocks/fractionStocks", format(Sys.Date(), "%Y%m%d"), ".csv"))

# Summary table of fraction stocks by layer -------------------------------
fracStockSummary <- ashNormStockLongWBulk %>% 
  group_by(year, treatment, depth, fraction) %>% 
  summarise(meanC = mean(Cstock, na.rm = TRUE),
            sdC = sd(Cstock, na.rm = TRUE),
            seC = plotrix::std.error(Cstock, na.rm = TRUE),
            meanN = mean(Nstock, na.rm = TRUE),
            sdN = sd(Nstock, na.rm = TRUE),
            seN = plotrix::std.error(Nstock, na.rm = TRUE)) %>% 
  ungroup()

#what proportion of mineral layer (below 35cm) SOC is light fraction? 
total_mineral <- fracStockSummary %>% 
  filter(depth == "mineral") %>% 
  group_by(treatment) %>% 
  summarise(
    total_bulkC = sum(meanC[fraction == "bulk"], na.rm = TRUE),
    
    total_lightC = sum(meanC[fraction == "L"], na.rm = TRUE),
    total_lightC_se = sqrt(
      sum(seC[fraction == "L"]^2, na.rm = TRUE)
    ),    .groups = "drop") %>% 
  mutate(percent_light = (total_lightC / total_bulkC) * 100)

#What percent of the total C pool is in the light fraction?
fracStorageSummary <- ashNormStockLongWBulk %>% 
  filter(depth == "total") %>% 
  select(c(year:depth), fraction, Cstock) %>% 
  tidyr::pivot_wider(names_from = fraction,
                     values_from = c(Cstock),
                     names_sep = "_") %>% 
  rowwise() %>% 
  mutate(totalCPool = sum(L, H),
         propCinLFrac = L/totalCPool,
         propCinHFrac = H/totalCPool) %>% 
  group_by(year, treatment) %>% 
  summarise(meanCinLFrac = mean(propCinLFrac, na.rm = TRUE),
            sdCinLFrac = sd(propCinLFrac, na.rm = TRUE),
            seCinLFrac = plotrix::std.error(propCinLFrac, na.rm = TRUE)) %>% 
  ungroup()

#Mean for all years: 
fracStorageSummaryGlobalMean <- fracStorageSummary %>% 
  summarise(mean = mean(meanCinLFrac, na.rm = TRUE)*100,
            se = plotrix::std.error(meanCinLFrac, na.rm = TRUE)*100)

bulkPercChange <- fracStockSummary %>% 
  filter(fraction == "bulk") %>% 
  dplyr::select(depth, treatment, meanC, seC) %>% 
  pivot_wider(names_from = treatment,
              values_from = c(meanC, seC)) %>% 
  mutate(pctChange_c = ((meanC_c - meanC_i) / meanC_i) * 100,
         sePctChange_c = 100 * sqrt((seC_c / meanC_i)^2 +
                                      ((seC_i * (meanC_c - meanC_i)) / (meanC_i^2))^2),
         pctChange_w = ((meanC_w - meanC_i) / meanC_i) * 100,
         sePctChange_w = 100 * sqrt((seC_w / meanC_i)^2 +
                                      ((seC_i * (meanC_w - meanC_i)) / (meanC_i^2))^2))

#calculate change over time (as stock differences and  % differences)
lightStockChange <- ashNormStockLong %>% 
  filter(fraction == "L") %>% 
  group_by(year, treatment, depth) %>% 
  #calculate mean and se for each group
  summarise(meanStock = mean(Cstock, na.rm = TRUE),
            seStock = plotrix::std.error(Cstock, na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(depth, treatment, meanStock, seStock) %>%
  #make it wide so we can find differences by row
  tidyr::pivot_wider(names_from = treatment,
                     values_from = c(meanStock, seStock),
                     names_sep = "_") %>% 
  mutate(
    #calculate differences in control and warming from 2009 (i)
    diff_c = meanStock_c - meanStock_i,
    se_diff_c = sqrt(seStock_c^2 + seStock_i^2),
    diff_w = meanStock_w - meanStock_i,
    se_diff_w = sqrt(seStock_w^2 + seStock_i^2),
    diff_i = 0,
    se_diff_i = seStock_i
  ) %>% 
  dplyr::select(depth, diff_c, se_diff_c, diff_w, se_diff_w, diff_i, se_diff_i) %>%
  tidyr::pivot_longer(cols = c(diff_c, diff_w, diff_i, se_diff_c, se_diff_w, se_diff_i),
                      names_to = c(".value", "treatment"),
                      names_pattern = "(diff|se_diff)_(c|w|i)") %>%
  rename(diffFrom2009 = diff,
         seFrom2009 = se_diff)

lightPercChange <- fracStockSummary %>% 
  filter(fraction == "L") %>% 
  dplyr::select(depth, treatment, meanC, seC) %>% 
  pivot_wider(names_from = treatment,
              values_from = c(meanC, seC)) %>% 
  mutate(pctChange_c = ((meanC_c - meanC_i) / meanC_i) * 100,
         sePctChange_c = 100 * sqrt((seC_c / meanC_i)^2 +
                                      ((seC_i * (meanC_c - meanC_i)) / (meanC_i^2))^2),
         pctChange_w = ((meanC_w - meanC_i) / meanC_i) * 100,
         sePctChange_w = 100 * sqrt((seC_w / meanC_i)^2 +
                                      ((seC_i * (meanC_w - meanC_i)) / (meanC_i^2))^2))

#calculate change over time (as stock differences and  % differences)
heavyStockChange <- ashNormStockLong %>% 
  filter(fraction == "H") %>% 
  group_by(year, treatment, depth) %>% 
  #calculate mean and se for each group
  summarise(meanStock = mean(Cstock, na.rm = TRUE),
            seStock = plotrix::std.error(Cstock, na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(depth, treatment, meanStock, seStock) %>%
  #make it wide so we can find differences by row
  tidyr::pivot_wider(names_from = treatment,
                     values_from = c(meanStock, seStock),
                     names_sep = "_") %>% 
  mutate(
    #calculate differences in control and warming from 2009 (i)
    diff_c = meanStock_c - meanStock_i,
    se_diff_c = sqrt(seStock_c^2 + seStock_i^2),
    diff_w = meanStock_w - meanStock_i,
    se_diff_w = sqrt(seStock_w^2 + seStock_i^2),
    diff_i = 0,
    se_diff_i = seStock_i
  ) %>% 
  dplyr::select(depth, diff_c, se_diff_c, diff_w, se_diff_w, diff_i, se_diff_i) %>%
  tidyr::pivot_longer(cols = c(diff_c, diff_w, diff_i, se_diff_c, se_diff_w, se_diff_i),
                      names_to = c(".value", "treatment"),
                      names_pattern = "(diff|se_diff)_(c|w|i)") %>%
  rename(diffFrom2009 = diff,
         seFrom2009 = se_diff)

heavyPercChange <- fracStockSummary %>% 
  filter(fraction == "H") %>% 
  dplyr::select(depth, treatment, meanC, seC) %>% 
  pivot_wider(names_from = treatment,
              values_from = c(meanC, seC)) %>% 
  mutate(pctChange_c = ((meanC_c - meanC_i) / meanC_i) * 100,
         sePctChange_c = 100 * sqrt((seC_c / meanC_i)^2 +
                                      ((seC_i * (meanC_c - meanC_i)) / (meanC_i^2))^2),
         pctChange_w = ((meanC_w - meanC_i) / meanC_i) * 100,
         sePctChange_w = 100 * sqrt((seC_w / meanC_i)^2 +
                                      ((seC_i * (meanC_w - meanC_i)) / (meanC_i^2))^2))

#How much of stock losses were driven by light fraction losses
lossesFromLight <- ashNormStockLongWBulk %>% 
  group_by(treatment, depth, fraction) %>% 
  summarise(meanC = mean(Cstock, na.rm = TRUE),
            .groups = "drop") %>% 
  filter(fraction %in% c("L", "bulk")) %>% 
  tidyr::pivot_wider(names_from = fraction,
                     values_from = meanC) %>% 
  group_by(depth) %>% 
  mutate(bulk_i  = bulk[treatment == "i"],
         light_i = L[treatment == "i"],
         delta_bulk  = bulk - bulk_i,
         delta_light = L - light_i) %>% 
  ungroup() %>% 
  filter(treatment %in% c("c", "w")) %>% 
  mutate(prop_bulk_decline_from_light = delta_light / delta_bulk)

# Statistical differences: ash normalized stocks --------------------------

## Light fraction stock differences  -----------------------------------------------------------
FracStockCalcCoreMod <- ashNormStockLong %>% 
  dplyr::select(-block, -plot, -year) %>% 
  filter(fraction == "L") %>% 
  mutate(depth = factor(depth),
         fence = factor(fence),
         fractionDepth = interaction(fraction, depth)) 

hist(FracStockCalcCoreMod$Cstock, breaks = 50)

#Use GLMM with gamma log distribution for positive right skewed data
m0 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              FracStockCalcCoreMod, family = gaussian())

m1 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              FracStockCalcCoreMod, family = Gamma(link = "log"))

AIC(m0, m1) # GLMM improves model fit

performance::check_model(m1)
DHARMa::simulateResiduals(fittedModel = m1, plot = TRUE) #over dispersion detected

m2 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ depth,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

m3 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ depth*fence,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

m4 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ fence,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

m5 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ treatment*depth,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

m6 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ treatment,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

m7 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ treatment*depth*fence,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

AIC(m1, m2, m3, m4, m5, m6, m7)
mod.fin <- m5
mod.light <- mod.fin

performance::check_model(mod.fin) #much better performance
DHARMa::simulateResiduals(fittedModel = mod.fin, plot = TRUE) #no overdispersion

compare <- emmeans::emmeans(mod.fin, pairwise ~ treatment | depth, adjust = "none")
compare

cld_light <- multcomp::cld(compare,
                           adjust = "none",
                           alpha = 0.1,
                           Letters = letters,
                           sort = FALSE) %>% 
  mutate(fraction = "Light") %>% 
  rename(lower.CL = asymp.LCL,
         upper.CL = asymp.UCL)

summary(mod.fin)
sjPlot::tab_model(mod.fin)

## Heavy fraction stock differences ----------------------------------------

####Total, mineral, deep mineral layers ----------------------------------------
FracStockCalcCoreMod <- ashNormStockLong %>% 
  #Make cores without heavy fraction stock have a very very small amount of heavy fraction stock
  mutate(Cstock = Cstock + 0.0001) %>% 
  dplyr::select(-block, -plot, -year) %>% 
  filter(!(depth %in% c("o", "organic"))) %>% 
  filter(fraction == "H") %>% 
  mutate(depth = factor(depth),
         fence = factor(fence),
         fractionDepth = interaction(fraction, depth))

mod <- lmer(Cstock ~ treatment*depth + (1|fence), FracStockCalcCoreMod)
performance::check_model(mod) #not great with model fit

#right skewed, non-normal
hist(FracStockCalcCoreMod$Cstock, breaks = 50)

#GLMM to account for non-linear data
m0 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              FracStockCalcCoreMod, family = gaussian())

m1 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              FracStockCalcCoreMod, family = Gamma(link = "log"))

AIC(m0, m1) # GLMM improves model fit

performance::check_model(m1)
DHARMa::simulateResiduals(fittedModel = m1, plot = TRUE) #over dispersion detected

m2 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ depth,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

m3 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ depth*fence,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

m4 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ fence,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

m5 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ treatment*depth,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

m6 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ treatment,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

m7 <- glmmTMB(Cstock ~ treatment*depth + (1|fence),
              dispformula = ~ treatment*depth*fence,
              FracStockCalcCoreMod, family = Gamma(link = "log"))

AIC(m1, m2, m3, m4, m5, m6, m7)

mod.fin <- m7
mod.heavy.min <- mod.fin
performance::check_model(mod.fin) #looks good
DHARMa::simulateResiduals(fittedModel = mod.fin, plot = TRUE)

compare <- emmeans::emmeans(mod.fin, pairwise ~ treatment | depth, adjust = "none")
compare

cld_heavy_min <- multcomp::cld(compare,
                               adjust = "none",
                               alpha = 0.1,
                               Letters = letters,
                               sort = FALSE) %>% 
  mutate(fraction = "Heavy") %>% 
  rename(lower.CL = asymp.LCL,
         upper.CL = asymp.UCL)

summary(mod.fin)
sjPlot::tab_model(mod.fin)

### Organic layer ----------------------------------------
FracStockCalcCoreMod <- ashNormStockLong %>% 
  dplyr::select(-block, -plot, -year) %>% 
  filter(fraction == "H" & depth %in% c("o")) %>% 
  mutate(depth = factor(depth),
         fence = factor(fence),
         fractionDepth = interaction(fraction, depth)) 

hist(FracStockCalcCoreMod$Cstock, breaks = 100) #highly 0 inflated

#Non-normal distribution! heavily right skewed
#Tweedie distribution deals with 0s and also positive skewed data

m1 <- glmmTMB(Cstock ~ treatment + (1|fence),
              FracStockCalcCoreMod, family = tweedie(link = "log"))

performance::check_model(m1)
DHARMa::simulateResiduals(fittedModel = m1, plot = TRUE) #heteroskedacity

m2 <- glmmTMB(Cstock ~ treatment + (1|fence),
              dispformula = ~ fence,
              FracStockCalcCoreMod, family = tweedie(link = "log"))

m3 <- glmmTMB(Cstock ~ treatment + (1|fence),
              dispformula = ~ treatment,
              FracStockCalcCoreMod, family = tweedie(link = "log"))

m4 <- glmmTMB(Cstock ~ treatment + (1|fence),
              dispformula = ~ treatment*fence,
              FracStockCalcCoreMod, family = tweedie(link = "log"))

AIC(m1, m2, m3, m4)
mod.fin <- m3
mod.heavy.org <- mod.fin
performance::check_model(mod.fin)

compare <- emmeans::emmeans(mod.fin, pairwise ~ treatment, adjust = "none")
compare
cld_heavy_org <- multcomp::cld(compare,
                               adjust = "none",
                               alpha = 0.1,
                               Letters = letters,
                               sort = FALSE) %>% 
  mutate(fraction = "Heavy",
         depth = "o") %>% 
  rename(lower.CL = asymp.LCL,
         upper.CL = asymp.UCL)

summary(mod.fin)
tab_model(mod.fin)

# C stock plots -----------------------------------------------------------

## Final stock comparison --------------------------------------------------
#compact letter displays
cldAll <- rbind(cld_light, cld_heavy_org, cld_heavy_min) %>% 
  mutate(.group = trimws(.group)) %>% 
  select(treatment, depth, fraction, emmean, SE, lower.CL, upper.CL, .group) %>%
  rename(mean_model = emmean,
         se_model = SE)

#make dataframe with summary
stockSum <- ashNormStockLong %>%
  group_by(fraction, depth, treatment) %>%
  summarise(meanC = mean(Cstock, na.rm = TRUE),
            sdC   = sd(Cstock, na.rm = TRUE),
            seC   = sdC / sqrt(sum(!is.na(Cstock))),
            n     = sum(!is.na(Cstock)),
            .groups = "drop") %>% 
  mutate(fraction = ifelse(fraction == "H", "Heavy", "Light")) %>% 
  full_join(cldAll,
            by = c("fraction", "depth", "treatment")) %>% 
  mutate(depth = factor(depth, 
                        levels = c("o", "m", "dm", "mineral", "total"),
                        labels = c("Organic (~0-35cm)",
                                   "Mineral (~35-55cm)",
                                   "Deep mineral (~55-75cm)",
                                   "Total mineral (~35-75cm)",
                                   "Total stocks")),
         fraction = factor(fraction, 
                           levels = c("Light", "Heavy")),
         labelPos = ifelse(depth == "Organic (~0-35cm)",
                           meanC + seC + 0.02,
                           meanC + seC + 0.4))


#Create dataframe for stacked plots
stockSumStack <- stockSum %>%
  filter(depth == "Total stocks") %>%
  mutate(fraction = factor(fraction, levels = c("Heavy", "Light"))) %>%
  arrange(treatment, fraction) %>%
  group_by(treatment) %>%
  mutate(ymin = lag(cumsum(meanC), default = 0),
         ymax = ymin + meanC,
         err_ymin = ymax - seC,
         err_ymax = ymax + seC) %>%
  ungroup()

#signficance levels and testing from Lathrop et al. 2025 GCB
totalCld <- stockSumStack %>% 
  filter(!fraction == "Heavy") %>%
  mutate(.group = c("a", "b", "c"), #Comes from Lathrop et al. 2025 GCB
         label_y = err_ymax + 1) 

totalStocksStacked <- ggplot(stockSumStack) +
  geom_rect(aes(xmin = as.numeric(treatment) - 0.4,
                xmax = as.numeric(treatment) + 0.4,
                ymin = ymin,
                ymax = ymax,
                fill = treatment,
                alpha = fraction), color = "black") +
  geom_errorbar(aes(x = treatment, ymin = err_ymin, ymax = err_ymax), width = 0.35) +
  scale_x_discrete(labels = stringr::str_wrap(c("i" = "2009: Initial",
                                                "c" = "2022: Ambient",
                                                "w" = "2022: Warming"),
                                              width = 10))+
  geom_text(data = totalCld, aes(x = treatment, y = label_y, label = .group, group = treatment),
            position = position_dodge(0.9),
            size = 5,
            inherit.aes = TRUE)+
  scale_alpha_manual(values = c("Heavy" = 0.7, "Light" = 0.2),
                     breaks = c("Light", "Heavy"),
                     labels = c("POC", "MAOC"),
                     name = "Fraction") +
  scale_fill_manual(values = treatPal, name = element_blank(), guide = "none")+
  labs(x = element_blank(), y = expression(paste("Carbon stock (kg ", m^-2, ")"))) +
  ts_theme +
  ggtitle("Total C stocks (~0–75 cm)")+
  theme(legend.position = "left",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 20))

totalStocksStacked

depthStocksLight <- ggplot(subset(stockSum, depth %in%c("Organic (~0-35cm)", "Total mineral (~35-75cm)") & fraction == "Light"), 
                           aes(x = treatment, y = meanC, fill = treatment))+
  geom_bar(stat = "identity", position = position_dodge(0.9), color = "black", alpha = 0.2)+
  geom_errorbar(aes(ymax = meanC + seC, ymin = meanC - seC), position = position_dodge(0.9), width = 0.35)+
  scale_fill_manual(values = treatPal, name = "Treatment", labels = c("2009: Initial", "2022: Ambient", "2022: Warming"))+
  facet_wrap(~depth, ncol = 1, scales = "free_y")+
  scale_x_discrete(labels = stringr::str_wrap(c("i" = "2009: Initial",
                                                "c" = "2022: Ambient",
                                                "w" = "2022: Warming"),
                                              width = 10))+
  geom_text(aes(x = treatment, y = meanC+2*seC, label = .group, group = treatment),
            position = position_dodge(0.9),
            size = 5,
            inherit.aes = TRUE)+
  labs(x = element_blank(), y = expression(paste("Carbon stock (kg ", m^-2, ")")))+
  ts_theme+
  theme(strip.background = element_rect(color="black", fill=fracPal[], size=1.5, linetype="solid"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 15),
        legend.position="none")+
  ggtitle("POC stocks by depth")
depthStocksLight
#ggsave(paste0(manDir, "figures/suppLightStocksDepth", Sys.Date(), ".png"), width = 10, height = 15, units = "cm", dpi = 300)

depthStocksHeavy <- ggplot(subset(stockSum, depth %in% c("Organic (~0-35cm)", "Total mineral (~35-75cm)") & fraction == "Heavy"), 
                           aes(x = treatment, y = meanC, fill = treatment))+
  geom_bar(stat = "identity", position = position_dodge(0.9), color = "black", alpha = 0.7)+
  geom_errorbar(aes(ymax = meanC + seC, ymin = meanC - seC), position = position_dodge(0.9), width = 0.35)+
  scale_fill_manual(values = treatPal, name = "Treatment", labels = c("2009: Initial", "2022: Ambient", "2022: Warming"))+
  facet_wrap(~depth, ncol = 1, scales = "free_y")+
  scale_x_discrete(labels = stringr::str_wrap(c("i" = "2009: Initial",
                                                "c" = "2022: Ambient",
                                                "w" = "2022: Warming"),
                                              width = 10))+
  geom_text(aes(x = treatment, y = labelPos, label = .group, group = treatment),
            position = position_dodge(0.9),
            size = 5,
            inherit.aes = TRUE)+
  labs(x = element_blank(), y = expression(paste("MAOC stock (kg ", m^-2, ")")))+
  ts_theme+
  theme(strip.background = element_rect(color="black", size=1.5, linetype="solid"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 15),
        legend.position="none")+
  ggtitle("MAOC stocks by depth")
depthStocksHeavy

depthComb <- depthStocksHeavy
depthComb

ggarrange(totalStocksStacked, depthComb, widths = c(1.2,1), labels = c("A)", "B)"))
#ggsave(paste0(manDir, "figures/stocksDepth", Sys.Date(), ".png"), width = 27, height = 20, units = "cm", dpi = 300)

# Radiocarbon of fraction analysis --------------------------------

##Fraction modern for calculating age analysis -----------------------------------------------------
dfFracMod <- soil %>% 
  dplyr::select(year, fence, treatment, coreNum, midDepth, fracModern, LFracF14C, HFracF14C) %>%
  tidyr::pivot_longer(cols = -c(year, treatment, coreNum, midDepth, fence)) %>% 
  dplyr::rename(groups = "name") 

fracModSum <- dfFracMod  %>% 
  filter((midDepth > 20 & midDepth < 80)) %>% 
  #first take the average of all depths (uneven sample sizes by depth)
  group_by(treatment, groups, midDepth) %>% 
  summarise(se = plotrix::std.error(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE)) %>% 
  #then average across the depths
  group_by(treatment, groups) %>% 
  summarise(se = plotrix::std.error(mean, na.rm = TRUE),
            mean = mean(mean, na.rm = TRUE)) %>% 
  mutate(meanAge = -8033*log(mean))

#age of radiocarbon
dfAge <- soil %>% 
  dplyr::select(year, fence, treatment, coreNum, midDepth, bulkAge, LFracAge, HFracAge) %>%
  tidyr::pivot_longer(cols = -c(year, treatment, coreNum, midDepth, fence)) %>% 
  dplyr::rename(groups = "name") 

ageSum <- dfAge  %>% 
  filter((midDepth > 20 & midDepth < 80)) %>% 
  #first take the average of all depths (uneven sample sizes by depth)
  group_by(treatment, groups, midDepth) %>% 
  summarise(se = plotrix::std.error(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE)) %>% 
  #then average across the depths
  group_by(treatment, groups) %>% 
  summarise(se = plotrix::std.error(mean, na.rm = TRUE),
            mean = mean(mean, na.rm = TRUE)) 

ageSum2 <- dfAge  %>% 
  filter((midDepth > 20 & midDepth < 80)) %>% 
  #first take the average of all depths (uneven sample sizes by depth)
  group_by(groups) %>% 
  summarise(se = plotrix::std.error(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE))

#delta 14C
dfDel14 <- soil %>% 
  dplyr::select(year, fence, treatment, coreNum, midDepth, LFracd14C, HFracd14C, 
                delta14 = d14C) %>%
  tidyr::pivot_longer(cols = -c(year, treatment, coreNum, midDepth, fence)) %>% 
  dplyr::rename(groups = "name") 

#test significant differences:
#is heavy fraction older than light fraction or bulk? 
df <- dfDel14 %>% 
  filter((midDepth > 20 & midDepth < 80))

mod <- lmer(value ~ groups*midDepth + (1|fence), data = df)
summary(mod)
pairs(emmeans::emmeans(mod, ~groups, adjust = "none"))

d14Sum <- dfDel14  %>% 
  filter((midDepth > 20 & midDepth < 80)) %>% 
  #first take the average of all depths (uneven sample sizes by depth)
  group_by(treatment, groups, midDepth) %>% 
  summarise(se = plotrix::std.error(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE)) %>% 
  #then average across the depths
  group_by(treatment, groups) %>% 
  summarise(se = plotrix::std.error(mean, na.rm = TRUE),
            mean = mean(mean, na.rm = TRUE))  %>% 
  mutate(groups = factor(groups, 
                         levels = c("delta14", "LFracd14C", "HFracd14C")))

d14Sum2 <- dfDel14  %>% 
  filter((midDepth > 20 & midDepth < 80)) %>% 
  #first take the average of all depths (uneven sample sizes by depth)
  group_by(groups) %>% 
  summarise(se = plotrix::std.error(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE))

ggplot(data = subset(d14Sum, (groups %in% c("HFracd14C", "LFracd14C", "delta14"))),
       aes(x = treatment, y = mean, fill = treatment))+
  geom_col(width = 0.7, alpha = 0.7, color = "black") +
  #geom_bar(stat = "identity", color = "black", alpha = 0.7)+
  geom_errorbar(aes(ymax = mean + se, ymin =  mean-se), width = 0.3)+
  scale_fill_manual(values = treatPal[])+
  facet_wrap(~groups, scales = "free", 
             labeller = as_labeller(c(`delta14` = "Bulk SOC",
                                      `LFracd14C` = "POC",
                                      `HFracd14C` = "MAOC")))+
  scale_x_discrete(labels=stringr::str_wrap(c("i" = "2009: Initial", "c" = "2022: Ambient",
                                              "w" = "2022: Warming"), width = 10))+  
  coord_cartesian(ylim = c(-500, -300)) +
  ylab(expression(paste(Delta^14,"C ‰")))+ xlab(element_blank())+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16))

#ggsave(paste0(manDir, "figures/c14BarSupp", Sys.Date(), ".png"), width = 28, height = 15, units = "cm", dpi = 300)

## Age offset analysis -----------------------------------------------------
dfRad <- soil %>% 
  dplyr::select(year, fence, treatment, coreNum, midDepth, LFracd14C, HFracd14C, 
                delta14 = d14C, LFracOffset, HFracOffset) %>%
  tidyr::pivot_longer(cols = -c(year, treatment, coreNum, midDepth, fence)) %>% 
  dplyr::rename(groups = "name") %>% 
  mutate(groups = factor(groups, 
                         levels = c("delta14", "LFracd14C", "HFracd14C", "LFracOffset", "HFracOffset" ), 
                         labels = c("Bulk sample", "Light fraction", "Heavy fraction", "Light fraction offset", "Heavy fraction offset")),
         treatment = factor(treatment, 
                            levels = c("i", "c", "w")))

radOffSum <- dfRad %>% 
  filter((midDepth >20 &midDepth < 90)) %>% 
  group_by(treatment, groups) %>% 
  summarise(se = plotrix::std.error(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE))

#test significant differences:
#are light offsets different from heavy offsets? 
df <- filter(subset(dfRad, !is.na(value)))
mod <- lmer(value ~ groups + (1|fence), data = df)
summary(mod)
pairs(emmeans::emmeans(mod, ~groups, adjust = "none"))

dfRadSumDepths <- dfRad %>% 
  group_by( midDepth, treatment, groups) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            se = plotrix::std.error(value, na.rm = TRUE))

#supplemental material
ggplot(data = subset(dfRadSumDepths, 
                     (groups %in% c("Bulk sample", "Light fraction", "Heavy fraction")) &
                       (midDepth >20 & midDepth < 80)),
       aes(x = midDepth, y = mean, color = groups))+
  geom_point()+
  geom_line(linewidth = 1)+
  geom_errorbar(aes(ymax = mean + se, ymin =  mean-se), width = 0.4, linewidth = 1)+
  scale_color_manual(values = offsetPal[], name = "Fraction", labels = c("Bulk SOC", "POC", "MAOC"))+
  facet_wrap(~treatment, scales = "free",labeller = as_labeller(c(`i` = "2009: initial",
                                                                  `c` = "2022: ambient",
                                                                  `w` = "2022: warming")))+
  scale_x_reverse()+coord_flip()+
  ylab(expression(paste(Delta^14,"C ‰")))+ xlab("Depth (cm)")+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))

#ggsave(paste0(manDir, "figures/suppFrac14CByTreat", Sys.Date(), ".png"), width = 20, height = 15, units = "cm", dpi = 300)

ggplot(data = subset(dfRadSumDepths, 
                     !(groups %in% c("Bulk sample", "Light fraction", "Heavy fraction")) & 
                       (midDepth >20 &midDepth < 80) ),
       aes(x = midDepth, y = mean, color = treatment))+
  geom_point()+
  geom_line(linewidth = 1)+
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_errorbar(aes(ymax = mean + se, ymin =  mean-se), width = 0.3, linewidth = 1)+
  scale_color_manual(values = treatPal[], name = "Treatment", labels = c("2009", "2022: control", "2022: warming"))+
  facet_wrap(~groups, scales = "free", labeller = as_labeller(c(`Light fraction offset` = "POC age offset",
                                                                `Heavy fraction offset` = "MAOC age offset")))+
  scale_x_reverse()+coord_flip()+
  ylab("Age offset (fraction age - bulk age)")+
  xlab("Depth (cm)")+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))
#ggsave(paste0(manDir, "figures/suppAgeOffset", Sys.Date(), ".png"), width = 20, height = 15, units = "cm", dpi = 300)


## Plot age offset manuscript figure ---------------------------------------
#Treatment difference in light fraction offsets
df <- filter(subset(dfRad, groups == "Light fraction offset" & (midDepth >20 & midDepth < 80)), !is.na(value))

mod <- lmer(value ~ treatment + (1|fence), data = df)
summary(mod)
modLight <- mod
compare <- emmeans(modLight, pairwise ~ treatment, adjust = "none")
compare

cld_light <- multcomp::cld(compare, adjust = "none", 
                           alpha = 0.1,
                           Letters = letters, 
                           sort = FALSE) %>% 
  mutate(groups = "Light fraction offset")

cld_light

#Treatment differences in heavy fraction offsets
df <- filter(subset(dfRad, groups == "Heavy fraction offset" & (midDepth > 20 & midDepth < 80)), !is.na(value))

mod <- lmer(value ~ treatment + (1|fence), data = df)
compare <- emmeans(mod, pairwise ~ treatment, adjust = "none")
compare

cld_heavy <- multcomp::cld(compare, adjust = "none", 
                           alpha = 0.1,
                           Letters = letters, 
                           sort = FALSE) %>% 
  mutate(groups = "Heavy fraction offset")

cldAll <- rbind(cld_light, cld_heavy) %>% 
  mutate(.group = trimws(.group),
         cld_label = as.character(.group))

dfRadSum <- dfRad %>% 
  filter(!(groups %in% c("Bulk sample", "Light fraction", "Heavy fraction"))) %>% 
  filter(midDepth >20 & midDepth < 80) %>% 
  group_by(treatment, groups) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            se = plotrix::std.error(value, na.rm = TRUE)) %>% 
  full_join(cldAll,
            by = c("groups", "treatment")) %>% 
  mutate(yPos = ifelse(groups == "Light fraction offset",
                       (mean - se) - 100,
                       mean + se + 100),
         groups = factor(groups, 
                         levels = c("Light fraction offset",
                                    "Heavy fraction offset")),
         groups = factor(groups, 
                         levels = c("Light fraction offset", "Heavy fraction offset")))

barFracOffset <- ggplot(data = subset(dfRadSum, !(groups %in% c("Bulk sample", "Light fraction", "Heavy fraction"))),
                        aes(x = treatment, y = mean, fill = treatment))+
  geom_bar(stat = "identity", color = "black", alpha = 0.7)+
  geom_errorbar(aes(ymax = mean + se, ymin =  mean-se), width = 0.3)+
  scale_fill_manual(values = treatPal[])+
  geom_text(aes(x = treatment, y = yPos, label = cld_label, group = treatment),
            position = position_dodge(0.9),
            size = 5,
            inherit.aes = TRUE)+
  facet_wrap(~groups, labeller = as_labeller(c(`Light fraction offset` = "POC",
                                               `Heavy fraction offset` = "MAOC")))+
  
  scale_x_discrete(labels=stringr::str_wrap(c("i" = "2009: Initial", "c" = "2022: Ambient",
                                              "w" = "2022: Warming"), width = 10))+
  ylab(expression(paste("Age offset (", ""^14,"C yrs)")))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16))
barFracOffset

# PRS Probe analysis ------------------------------------------------------
prs <- read.csv("datasets/prs.csv") %>% select(-X)
meltPrs <- reshape::melt(as.data.frame(prs), id=c("sampleID", "treatment", "treatmentFull", "fence", "plot", "block")) %>% 
  mutate(value = as.numeric(value))

## Plot all species --------------------------------------------------------
ggplot(data = meltPrs, aes(x = variable, y = value, fill = treatment))+
  geom_boxplot()+
  facet_wrap(~variable, scales = "free")+
  scale_fill_manual(values = treatPal[-1])

#elements of interest
redoxVars <- c("NH4N", "NH4N", "Fe", "Mn", "S", "meanRedox")
allVars <- c("NO3-N", "NH4-N", "Fe", "Mn", "S", "P", "Al", "Mean of redox")

plotDf <- meltPrs %>%
  filter(!(variable == "S" & value > 100)) %>% #remove outliers in sulfur plot
  mutate(variable = recode(variable,
                           "NO3N" = "NO3-N",
                           "NH4N" = "NH4-N",
                           "meanRedox" = "Mean of redox")) %>%
  subset(variable %in% allVars) %>%
  mutate(variable = as.character(variable),
         variable = factor(variable, levels = allVars))

#Add signficance levels
df = plotDf
mod <- lme(value ~ variable * treatment,
           random = ~1 | fence,
           weights = varIdent(form = ~1 | variable),
           data = df,
           control = lmeControl(maxIter = 100, msMaxIter = 100),
           method = "REML")

#performance::check_model(mod1)
emm_pairs <- pairs(emmeans::emmeans(mod, ~treatment | variable, adjust = "none"))
emm_pairs

emm <- emmeans(mod, ~ treatment | variable)

yPos <- plotDf %>%
  group_by(variable) %>%
  summarise(y = max(value, na.rm = TRUE) * 1.08)

pvals <- pairs(emm, adjust = "none") %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame() %>% 
  select(variable, contrast, p.value) %>% 
  left_join(yPos, by = "variable") %>% 
  mutate(
    label = case_when(
      p.value < 0.001 ~ "p < 0.001",
      TRUE ~ paste0("p = ", signif(p.value, 2))))

facet_labeller <- function(x) {
  ifelse(x %in% c(redoxVars, "NH4-N"),
         paste0(x, "*"),
         x)
}

ggplot(data = plotDf,
       aes(x = treatment, y = value, fill = treatment)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~variable, ncol = 4, scales = "free_y", drop = FALSE,
             labeller = labeller(variable = facet_labeller)) +
  scale_fill_manual( values = treatPal[-1], name = "Treatment",
                     labels = c("2022: Ambient", "2022: Warming")) +
  geom_text(data = pvals,
            aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE,
            size = 6)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = "bottom")+
  ylab(expression(paste("Concentration (", mu, "g / 10 cm"^-2, ")")))
#ggsave(paste0(manDir, "figures/suppPRS", Sys.Date(), ".png"), width = 20, height = 15, units = "cm", dpi = 300)

# make final figure with 
plotDf <- meltPrs %>%
  mutate(variable = recode(variable,
                           "NO3N" = "NO3-N",
                           "NH4N" = "NH4-N",
                           "meanRedox" = "Mean of redox")) %>%
  filter(variable %in% c("Fe", "Mean of redox")) %>%
  mutate(variable = as.character(variable),
         variable = factor(variable, levels = c("Fe", "Mean of redox")))

#Add signficance levels
df = plotDf
mod <- lme(value ~ variable * treatment,
           random = ~1 | fence,
           weights = varIdent(form = ~1 | variable),
           data = df,
           control = lmeControl(maxIter = 100, msMaxIter = 100),
           method = "REML")

#performance::check_model(mod1)
emm_pairs <- pairs(emmeans::emmeans(mod, ~treatment | variable, adjust = "none"))
emm_pairs

emm <- emmeans(mod, ~ treatment | variable)

yPos <- plotDf %>%
  group_by(variable) %>%
  summarise(y = max(value, na.rm = TRUE) * 1.08)

pvals <- pairs(emm, adjust = "none") %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame() %>% 
  select(variable, contrast, p.value) %>% 
  left_join(yPos, by = "variable") %>% 
  mutate(
    label = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01  ~ paste0("p = ",signif(p.value, 2)),
      p.value < 0.05  ~ "**",
      p.value < 0.1 ~ "**",
      TRUE ~ paste0("p = ", signif(p.value, 2))))

a <- ggplot(data = subset(plotDf, variable %in% c("Fe", "Mean of redox")),
            aes(x = treatment, y = value, fill = treatment)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~variable, ncol = 4, scales = "free_y", drop = FALSE) +
  scale_fill_manual( values = treatPal[-1], name = "Treatment",
                     labels = c("2022: Ambient", "2022: Warming")) +
  geom_text(data = pvals,
            aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE,
            size = 6)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom")+
  ylab(expression(paste("Concentration (", mu, "g / 10 cm"^-2, ")")))
a

## Water table depth and PRS data ------------------------------------------
annualDataset <- read.csv("datasets/allEnvData.csv") %>% select(-X)
prsEnv <- annualDataset %>% 
  filter(year == 2022) %>% 
  mutate(treatment = ifelse(plot > 4, "w", "c")) %>% 
  left_join(select(prs, -block, -treatmentFull), by = c("fence", "plot", "treatment"))

b <- ggplot(data = subset(prsEnv, meanRedox < 3), aes(x = wtd.mean, y = meanRedox))+
  geom_point()+
  #ylim(0, 3)+
  geom_vline(xintercept = 0, linetype = "dashed")+
  ggpmisc::stat_poly_line(aes(x = wtd.mean, y = meanRedox), color = "black", alpha = 0.25) +
  ggpmisc::stat_poly_eq(aes(x = wtd.mean, y = meanRedox,
                            label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x,
                        label.x = "right",
                        label.y = "top") +  scale_x_reverse()+
  ylab("Mean of redox-sensitive ion loadings (-)")+
  xlab("Depth to reach water table in 2022 (cm)")+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
b

mod <- lmer(meanRedox ~ wtd.mean + (1|fence), data = prsEnv)
drop1(mod, test = "Chisq")
tab_model(mod)

ggarrange(a, b, widths = c(1, 0.8), labels = c("A)", "B)"))
#ggsave(paste0(manDir, "figures/prsEnv", Sys.Date(), ".png"), width = 25, height = 15, units = "cm", dpi = 300)
# XRF Fe analysis  -----------------------------------------------

## Fe and wtd plots --------------------------------------------------------

meanWtdTreatment <- annualDataset %>% 
  filter(year == 2009 | year == 2022) %>% 
  mutate(treatment = ifelse(plot>4, "w", "c"),
         treatment = ifelse(year == 2009, "initial", treatment)) %>% 
  group_by(year, treatment) %>% 
  summarise(meanWtd = mean(wtd.mean, na.rm = TRUE),
            sdWtd = sd(wtd.mean, na.rm = TRUE),
            seWtd = plotrix::std.error(wtd.mean, na.rm = TRUE),
            meanWtdSd = mean(wtd.sd, na.rm = TRUE)) %>% 
  ungroup() 

coreWtdTreatment <- core %>% 
  filter(year == 2009 | year == 2022) %>% 
  mutate(treatment = ifelse(plot>4, "w", "c"),
         treatment = ifelse(year == 2009, "initial", treatment),
         CtoN = C/N) %>% 
  group_by(treatment, midDepth) %>% 
  summarise(meanFe = mean(Fe, na.rm = TRUE),
            sdFe = sd(Fe, na.rm = TRUE),
            seFe = plotrix::std.error(Fe, na.rm = TRUE),
            meanHFrac = mean(HFracC, na.rm = TRUE),
            sdHFrac = sd(HFracC, na.rm = TRUE),
            meanLFrac = mean(LFracC, na.rm = TRUE),
            sdLFrac = sd(LFracC, na.rm = TRUE),
            meanHFrac14C = mean(HFracd14C, na.rm = TRUE),
            sdHFrac14C= sd(HFracd14C, na.rm = TRUE),
            meanLFrac14C = mean(LFracd14C, na.rm = TRUE),
            sdLFrac14C= sd(LFracd14C, na.rm = TRUE),
            meanBulk14C = mean(delta14, na.rm = TRUE),
            sdBulk14C = sd(delta14, na.rm = TRUE)) %>% 
  mutate(treatment = factor(treatment, 
                            levels = c("initial", "c", "w"))) %>% 
  ungroup()

treat_scale <- scale_color_manual(name   = "Treatment",
                                  values = c(initial = "#808080", c = "#377EB8", w = "#E41A1C"),
                                  # breaks = treat_levels,
                                  # limits = treat_levels,
                                  labels = c(initial = "2009: Initial",
                                             c       = "2022: Ambient",
                                             w       = "2022: Warming"),
                                  drop   = FALSE)
legend_fix <- guides(color = guide_legend(override.aes = list(linetype = "solid",
                                                              shape    = 16,
                                                              size     = 1.2,
                                                              alpha    = 1)))

iPlot <- ggplot(data = subset(coreWtdTreatment, treatment =="initial" & midDepth < 80), 
                aes(x = midDepth, y = meanFe, group = treatment))+
  geom_vline(data = subset(meanWtdTreatment, treatment == "initial"), 
             aes(xintercept = meanWtd, color = treatment), alpha = 0.7, size = 1, show.legend = FALSE)+
  geom_vline(data = subset(meanWtdTreatment, treatment == "initial"), 
             aes(xintercept = meanWtd + meanWtdSd, color = treatment), alpha = 0.7, size = 1,  lty = "dashed", show.legend = FALSE)+
  geom_vline(data = subset(meanWtdTreatment, treatment == "initial"), 
             aes(xintercept = meanWtd - meanWtdSd, color = treatment), alpha = 0.7, size = 1, lty = "dashed", show.legend = FALSE)+
  geom_point(aes(color = treatment))+
  geom_line(aes(color = treatment), alpha = 0.7, size = 1)+
  geom_errorbar(aes(ymin = meanFe-seFe, ymax = meanFe+seFe, color = treatment), size = 1, width = 1.5)+
  treat_scale + legend_fix+
  scale_x_reverse(limits = c(71, -1),
                  breaks = c(0, 20, 40, 60))+ coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.position = "bottom",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  xlab("Depth (cm)")+
  ylab(expression(Fe~"("~mg~kg^{-1}~soil~")"))
iPlot 

cPlot <- ggplot(data = subset(coreWtdTreatment, treatment != "w" & midDepth < 80), 
                aes(x = midDepth, y = meanFe, group = treatment))+
  geom_vline(data = subset(meanWtdTreatment, treatment != "w" & treatment != "initial"), 
             aes(xintercept = meanWtd, color = treatment), alpha = 0.7, size = 1)+
  geom_vline(data = subset(meanWtdTreatment, treatment != "w" & treatment != "initial"), 
             aes(xintercept = meanWtd + meanWtdSd, color = treatment), alpha = 0.7, size = 1,  lty = "dashed")+
  geom_vline(data = subset(meanWtdTreatment, treatment != "w" & treatment != "initial"), 
             aes(xintercept = meanWtd - meanWtdSd, color = treatment), alpha = 0.7, size = 1, lty = "dashed")+
  geom_point(aes(color = treatment))+
  geom_line(aes(color = treatment), alpha = 0.7, linewidth = 1)+
  geom_errorbar(aes(ymin = meanFe-seFe, ymax = meanFe+seFe, color = treatment), size = 1, width = 1.5)+
  treat_scale + legend_fix+
  scale_x_reverse(limits = c(71, -1),
                  breaks = c(0, 20, 40, 60))+ coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.position = "bottom",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  xlab("Depth (cm)")+
  ylab(expression(Fe~"("~mg~kg^{-1}~soil~")"))
cPlot 

wPlot <- ggplot(data = subset(coreWtdTreatment, treatment != "c" & midDepth < 80), 
                aes(x = midDepth, y = meanFe, group = treatment))+
  geom_vline(data = subset(meanWtdTreatment, treatment != "c" & treatment != "initial"), 
             aes(xintercept = meanWtd, color = treatment), alpha = 0.7, size = 1)+
  geom_vline(data = subset(meanWtdTreatment, treatment != "c" & treatment != "initial"), 
             aes(xintercept = meanWtd + meanWtdSd, color = treatment), alpha = 0.7, size = 1,  lty = "dashed")+
  geom_vline(data = subset(meanWtdTreatment, treatment != "c" & treatment != "initial"), 
             aes(xintercept = meanWtd - meanWtdSd, color = treatment), alpha = 0.7, size = 1, lty = "dashed")+
  geom_point(aes(color = treatment))+
  geom_line(aes(color = treatment), alpha = 0.7, size = 1)+
  geom_errorbar(aes(ymin = meanFe-seFe, ymax = meanFe+seFe, color = treatment), 
                alpha = 0.7, size = 1, width = 1.5)+
  treat_scale + legend_fix+
  scale_x_reverse(limits = c(71, -1),
                  breaks = c(0, 20, 40, 60))+ coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.position = "bottom",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+  
  xlab("Depth (cm)")+
  ylab(expression(Fe~"("~mg~kg^{-1}~soil~")"))
wPlot 

legend_df <- data.frame(treatment = factor(c("initial", "c", "w"),
                                           levels = c("initial", "c", "w")),
                        midDepth = Inf,
                        meanFe   = Inf)

legend_point <- geom_point(data = legend_df,
                           aes(x = midDepth, y = meanFe, color = treatment),
                           size = 3,
                           alpha = 0,
                           inherit.aes = FALSE,
                           show.legend = TRUE)

legend_line <- geom_line(data = legend_df,
                         aes(x = midDepth, y = meanFe, color = treatment, group = treatment),
                         size = 1,
                         alpha = 0,
                         inherit.aes = FALSE,
                         show.legend = TRUE)

cPlot <- cPlot + legend_line + legend_point
wPlot <- wPlot + legend_line + legend_point 
iPlot <- iPlot + legend_line + legend_point

y_scale <- scale_y_continuous(limits = range(coreWtdTreatment$meanFe, na.rm = TRUE),
                              expand = expansion(mult = c(0.02, 0.02)))

ggarrange(iPlot + y_scale, ((cPlot + y_scale) + rremove("ylab")), ((wPlot + y_scale) + rremove("ylab")), 
          widths = c(1, 0.9, 0.9), nrow = 1, common.legend = TRUE, legend = "bottom")

#ggsave(paste0(manDir, "figures/FeWtd", Sys.Date(), ".png"), width = 20, height = 15, units = "cm", dpi = 300)

# Water table variability analysis -----------------------------------------
## Water table change over time --------------------------------------------
wtdSum <- annualDataset %>% 
  filter(year == 2009 | year == 2022) %>% 
  group_by(year, treatment2) %>% 
  summarise(meanWtd = mean(wtd.mean, na.rm = TRUE),
            seMeanWtd = plotrix::std.error(wtd.mean, na.rm = TRUE),
            meanSdWtd = mean(wtd.sd, na.rm = TRUE),
            seSdWtd = plotrix::std.error(wtd.sd, na.rm = TRUE))

#wtd mean change over time
a <- ggplot(annualDataset, aes(x = as.factor(year), y = wtd.mean, fill= treatment2))+
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1, alpha = 0.8)+
  geom_hline(yintercept = 25, linetype = "dashed", color = "brown", linewidth = 1, alpha = 0.8)+
  geom_boxplot(alpha = 0.6)+
  scale_y_reverse()+
  ylab("Water table depth (cm below surface)")+
  xlab("Year")+
  scale_fill_manual(values = treatPal[], name = "Treatment", 
                    labels = c("i" = "Initial",
                               "c" = "Ambient",
                               "w" = "Warming"))+
  theme(axis.text.x = element_text(size = 17, angle = 45, vjust = 0.5, hjust=0.5),
        legend.position = "top",
        axis.title.y = element_text(size = 18))
a

#wtd sd summarized
b <- ggplot(subset(annualDataset, year == 2009 | year == 2022), 
            aes(x = treatment2, y = wtd.sd, fill = treatment2))+
  geom_boxplot(alpha = 0.6)+
  scale_x_discrete(labels = stringr::str_wrap(c("2009: Initial", "2022: Ambient", "2022: Warming"), width = 10))+
  scale_fill_manual(values = treatPal, name = element_blank(), 
                    labels = c("2009 initial", "2022 ambient", "2022 warming"))+
  ylab("Water table variability (sd of cm)")+
  scale_y_reverse()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 18))
b

#wtd sd summarized
c <- ggplot(subset(annualDataset, year == 2009 | year == 2022), 
            aes(x = treatment2, y = t40.filled.mean, fill = treatment2))+
  geom_boxplot(alpha = 0.6)+
  scale_x_discrete(labels = stringr::str_wrap(c("2009: Initial", "2022: Ambient", "2022: Warming"), width = 10))+
  scale_fill_manual(values = treatPal, name = element_blank(), 
                    labels = c("2009 initial", "2022 ambient", "2022 warming"))+
  ylab("40 cm soil temperature (°C)")+
  #scale_y_reverse()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 18))
c

ggarrange(c, a, b, ncol = 3, widths = c(0.7, 1, 0.7), heights = c(1, 0.75, 0.75), labels = c("A)", "B)", "C)"))
#ggsave(paste0(manDir, "figures/experimentWtdChange", Sys.Date(), ".png"), width = 35, height = 15, units = "cm", dpi = 300)


# Environmental analysis --------------------------------------------------

#combine stock data with annual dataset (with just 2009 and 2022)
stockEnv <- ashNormStock %>% 
  filter(fence != 7) %>% 
  select(-treatment, -block) %>% 
  left_join(filter(annualDataset, (year == 2009 | year == 2022)), by = c("year", "fence", "plot")) 

stockEnv$wtd.mean.sc <- scale(stockEnv$wtd.mean) 
stockEnv$wtd.sd.sc <- scale(stockEnv$wtd.sd)
stockEnv$t40.filled.mean.sc <- scale(stockEnv$t40.filled.mean)
stockEnv$winter.t40.mean.sc <- scale(stockEnv$winter.t40.mean)

## Controls on spatial patterns of heavy fraction stock --------------------------------
#hypothesized vars: 
#wtd.mean.sc, wtd.sd.sc, t40.filled.mean.sc, winter.t40.mean.sc

#2022 data paired with measurements from 2022

#Mineral layer
df <- subset(stockEnv, (depth == "m" & year == 2022))
hist(df$Cstock.H, breaks = 20) #fairly normally distributed

mod <- lmer(Cstock.H ~ wtd.mean.sc + wtd.sd.sc + t40.filled.mean.sc + winter.t40.mean.sc + (1|fence), data =  df)
mod.full.m <- mod 
drop1(mod, test = "Chisq")
#no significant predictors

#3. Third layer: deep mineral 
df <- subset(stockEnv, (depth == "dm" & year == 2022))

mod <- lmer(Cstock.H ~ wtd.mean.sc + wtd.sd.sc + t40.filled.mean.sc + winter.t40.mean.sc + (1|fence), data =  df)
mod.full.dm <- mod 
drop1(mod, test = "Chisq")

mod <- lmer(Cstock.H ~ wtd.mean.sc + (1|treatment), data = df)
tab_model(mod)
performance::check_model(datawizard::standardize(mod))
mod.dm <- mod

#create a table of the full and final models
#full models
m_tbl <- parameters::model_parameters(mod.full.m, standardize = "refit") %>%
  mutate(Model = "2022 Mineral stocks (~35–55 cm)")

dm_tbl <- parameters::model_parameters(mod.full.dm, standardize = "refit") %>%
  mutate(Model = "2022 Deep Mineral stocks (~55–75 cm)")

model_colors <- c("2022 Mineral stocks (~35–55 cm)" = alpha("#8B775F", alpha = 0.65),
                  "2022 Deep Mineral stocks (~55–75 cm)" = alpha("#6B472C", alpha = 0.65))

library(gt)
full_table_gt <- bind_rows(m_tbl, dm_tbl) %>%
  filter(!grepl("Intercept", Parameter)) %>%
  dplyr::select(Model, Parameter, Coefficient, CI_low, CI_high, SE, p) %>% 
  mutate(Model = factor(Model, levels = c("2022 Mineral stocks (~35–55 cm)", 
                                          "2022 Deep Mineral stocks (~55–75 cm)"))) %>% 
  filter(Parameter != "SD (Observations)") %>% 
  gt(groupname_col = "Model") %>%
  text_case_match("wtd.sd.sc" ~ "Water table variability", 
                  "wtd.mean.sc"~"Water table depth",
                  "winter.t40.mean.sc"~"Winter soil temp (40cm)",
                  "t40.filled.mean.sc"~"Summer soil temp (40cm)") %>% 
  # Round numeric columns to 2 decimal places
  fmt_number(columns = vars(Coefficient, SE, CI_low, CI_high, p),
             decimals = 2) %>%
  # Define the column labels
  cols_label(Parameter = "Predictor",
             Coefficient = "Std. Estimate",
             SE = "Std. Error",
             CI_low = "CI Low",
             CI_high = "CI High",
             p = "p-value") %>%
  tab_header(title = "Full spatial model by soil layer") %>%
  # Apply custom colors to each row group (model layer)
  tab_style(style = list(
    cell_text(align = "center", weight = "bold"),
    cell_fill(color = model_colors["2022 Mineral stocks (~35–55 cm)"])),
    locations = cells_row_groups(groups = "2022 Mineral stocks (~35–55 cm)")) %>%
  tab_style(style = list(
    cell_text(align = "center", weight = "bold"),
    cell_fill(color = model_colors["2022 Deep Mineral stocks (~55–75 cm)"])),
    locations = cells_row_groups(groups = "2022 Deep Mineral stocks (~55–75 cm)")) %>%
  # Center the row group labels (titles)
  tab_spanner(label = "", columns = everything())
full_table_gt 
#gtsave(full_table_gt , filename = paste0(manDir, "figures/spatialModelFull.png"))

#Final model table

#no mineral layer parameters improved model fit

dm_tbl <- parameters::model_parameters(mod.dm, standardize = "refit") %>%
  mutate(Model = "2022 Deep Mineral stocks (~55–75 cm)")

r2_fin_dm <- performance::r2_nakagawa(mod.dm)
r2_fin_text_dm <- paste0("R²m = ", round(r2_fin_dm$R2_marginal, 2),
                         ", R²c = ", round(r2_fin_dm$R2_conditional, 2))

final_table_gt <- bind_rows(dm_tbl) %>%
  filter(!grepl("Intercept", Parameter)) %>%
  dplyr::select(Model, Parameter, Coefficient, CI_low, CI_high, SE, p) %>% 
  mutate(Model = factor(Model, levels = c("2022 Deep Mineral stocks (~55–75 cm)"))) %>% 
  filter(Parameter != "SD (Observations)") %>% 
  gt(groupname_col = "Model") %>%
  text_case_match("wtd.sd.sc" ~ "Water table variability", 
                  "wtd.mean.sc"~"Water table depth",
                  "winter.t40.mean.sc"~"Winter soil temp (40cm)",
                  "t40.filled.mean.sc"~"Summer soil temp (40cm)") %>% 
  # Round numeric columns to 2 decimal places
  fmt_number(columns = vars(Coefficient, SE, CI_low, CI_high, p),
             decimals = 2) %>%
  # Define the column labels
  cols_label(Parameter = "Predictor",
             Coefficient = "Std. Estimate",
             SE = "Std. Error",
             CI_low = "CI Low",
             CI_high = "CI High",
             p = "p-value") %>%
  tab_header(title = "Final spatial model by soil layer") %>%
  # Apply custom colors to each row group (model layer)
  tab_style(style = list(
    cell_text(align = "center", weight = "bold"),
    cell_fill(color = model_colors["2022 Deep Mineral stocks (~55–75 cm)"])),
    locations = cells_row_groups(groups = "2022 Deep Mineral stocks (~55–75 cm)")) %>%
  # Center the row group labels (titles)
  tab_spanner(label = "", columns = everything())

final_table_gt %>% 
  tab_source_note(source_note = paste0(
    "Model fit: ", r2_fin_text_dm))
#gtsave(final_table_gt , filename = paste0(manDir, "figures/spatialModelFin.png"))

## Dataset with change in stocks and mean plot environmental variables---------------------------------------

#create dataframe that averages all values by plot for the 13 years of CiPHER
env <- annualDataset %>%
  dplyr::select(-year, -treatment, -treatment2) %>% 
  group_by(fence, plot) %>% 
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>% 
  mutate(treatment = ifelse(plot <= 4,'c', 'w')) %>% 
  ungroup() %>%
  mutate(treatment = ifelse(plot <= 4,'Control', 'Warming'),
         treatmentFull = ifelse(plot == 2 | plot == 4,
                                'Control',
                                ifelse(plot == 1 | plot == 3,
                                       'Air Warming',
                                       ifelse(plot == 6 | plot == 8,
                                              'Soil Warming',
                                              'Air + Soil Warming'))),
         block = ifelse(fence <= 2,
                        'a',
                        ifelse(fence <= 4,
                               'b',
                               'c'))) %>% 
  rowwise() %>% 
  mutate(treatment = as.factor(treatment),
         treatmentNum = ifelse(treatment == "Control", 1, 2)) %>% 
  dplyr::select(-treatmentFull, -treatment) %>% 
  ungroup() 

## scale variables
env$wtd.mean.sc <- scale(env$wtd.mean) 
env$wtd.sd.sc <- scale(env$wtd.sd)
env$t40.filled.mean.sc <- scale(env$t40.filled.mean)
env$winter.t40.mean.sc <- scale(env$winter.t40.mean)

#calculate the change in stocks by fraction and depth
stockDiff <- ashNormStock %>%
  filter(year %in% c(2009, 2022)) %>%
  group_by(fence, depth) %>%
  mutate(L_2009 = Cstock.L[year == 2009][1],
         H_2009 = Cstock.H[year == 2009][1]) %>%
  ungroup() %>%
  # Compute change ONLY for 2022 rows
  mutate(L_change = ifelse(year == 2022, Cstock.L - L_2009, NA),
         H_change = ifelse(year == 2022, Cstock.H - H_2009, NA)) %>% 
  filter(year != 2009) %>% 
  #add in environmental variables 
  left_join(env, by = c("fence", "plot"))

## Controls on change in heavy fraction stocks ---------------------------------
#hypothesized variables
#wtd.mean.sc, wtd.sd.sc, daysSat25cm.sc, daysDry25cm.sc, tp.sc, t40.filled.mean.sc, winter.t40.mean.sc

#Mineral layer 
df <- subset(stockDiff,depth %in% c("m"))
mod <- lmer(H_change ~wtd.mean.sc + wtd.sd.sc + t40.filled.mean.sc + winter.t40.mean.sc + (1|fence), data =  df)
mod.full.m <- mod
drop1(mod, test = "Chisq")
mod <- lmer(H_change ~ wtd.sd.sc + t40.filled.mean.sc + (1|fence), data =  df)
sjPlot::tab_model(mod)
fin.m <- mod

#Deep mineral layer
df <- subset(stockDiff,depth %in% c("dm"))
mod <- lmer(H_change ~wtd.mean.sc + wtd.sd.sc + t40.filled.mean.sc + winter.t40.mean.sc + (1|fence), data =  df)
mod.full.dm <- mod
drop1(mod, test = "Chisq")
#no significant predictors in final model

#create a table of the full and final models
#full models
m_tbl <- parameters::model_parameters(mod.full.m, standardize = "refit") %>%
  mutate(Model = "Change in MAOC stocks in mineral layer (~35–55 cm)")

dm_tbl <- parameters::model_parameters(mod.full.dm, standardize = "refit") %>%
  mutate(Model = "Change in MAOC stocks in deep mineral layer (~55–75 cm)")

model_colors <- c("Change in MAOC stocks in mineral layer (~35–55 cm)" = alpha("#8B775F", alpha = 0.65),
                  "Change in MAOC stocks in deep mineral layer (~55–75 cm)"= alpha("#6B472C", alpha = 0.65))

full_table_gt <- bind_rows(m_tbl, dm_tbl) %>%
  filter(!grepl("Intercept", Parameter)) %>%
  dplyr::select(Model, Parameter, Coefficient, CI_low, CI_high, SE, p) %>% 
  mutate(Model = factor(Model, levels = c("Change in MAOC stocks in mineral layer (~35–55 cm)", 
                                          "Change in MAOC stocks in deep mineral layer (~55–75 cm)"))) %>% 
  filter(Parameter != "SD (Observations)") %>% 
  gt(groupname_col = "Model") %>%
  text_case_match("wtd.sd.sc" ~ "Water table variability",
                  "wtd.mean.sc"~"Water table depth",
                  "winter.t40.mean.sc"~"Winter soil temp (40cm)",
                  "t40.filled.mean.sc"~"Summer soil temp (40cm)") %>% 
  # Round numeric columns to 2 decimal places
  fmt_number(columns = vars(Coefficient, SE, CI_low, CI_high, p),
             decimals = 2) %>%
  # Define the column labels
  cols_label(Parameter = "Predictor",
             Coefficient = "Std. Estimate",
             SE = "Std. Error",
             CI_low = "CI Low",
             CI_high = "CI High",
             p = "p-value") %>%
  tab_header(title = "Full difference model by soil layer") %>%
  # Apply custom colors to each row group (model layer)
  tab_style(style = list(cell_text(align = "center", weight = "bold"),
                         cell_fill(color = model_colors["Change in MAOC stocks in mineral layer (~35–55 cm)"])),
            locations = cells_row_groups(groups ="Change in MAOC stocks in mineral layer (~35–55 cm)")) %>%
  tab_style(style = list(cell_text(align = "center", weight = "bold"),
                         cell_fill(color = model_colors["Change in MAOC stocks in deep mineral layer (~55–75 cm)"])),
            locations = cells_row_groups(groups = "Change in MAOC stocks in deep mineral layer (~55–75 cm)")) %>%
  # Center the row group labels (titles)
  tab_spanner(label = "", columns = everything())

full_table_gt 
#gtsave(full_table_gt , filename = paste0(manDir, "figures/differenceModelFull.png"))

m_tbl <- parameters::model_parameters(fin.m, standardize = "refit") %>%
  mutate(Model = "Change in MAOC stocks in mineral layer (~35–55 cm)")
r2_fin_m <- performance::r2_nakagawa(fin.m)
r2_fin_text_m <- paste0("R²m = ", round(r2_fin_m$R2_marginal, 2),
                        ", R²c = ", round(r2_fin_m$R2_conditional, 2))

final_table_gt <- bind_rows(m_tbl) %>%
  filter(!grepl("Intercept", Parameter)) %>%
  dplyr::select(Model, Parameter, Coefficient, CI_low, CI_high, SE, p) %>% 
  mutate(Model = factor(Model, levels = c("Change in MAOC stocks in mineral layer (~35–55 cm)", 
                                          "Change in MAOC stocks in deep mineral layer (~55–75 cm)"))) %>% 
  filter(Parameter != "SD (Observations)") %>% 
  gt(groupname_col = "Model") %>%
  text_case_match("wtd.sd.sc" ~ "Water table variability",
                  "wtd.mean.sc"~"Water table depth",
                  "winter.t40.mean.sc"~"Winter soil temp (40cm)",
                  "t40.filled.mean.sc"~"Summer soil temp (40cm)",
                  "SD (Observations)"~"SD (Observations)") %>% 
  # Round numeric columns to 2 decimal places
  fmt_number(columns = vars(Coefficient, SE, CI_low, CI_high, p),
             decimals = 2) %>%
  # Define the column labels
  cols_label(Parameter = "Predictor",
             Coefficient = "Std. Estimate",
             SE = "Std. Error",
             CI_low = "CI Low",
             CI_high = "CI High",
             p = "p-value") %>%
  tab_header(title = "Final difference model by soil layer") %>%
  # Apply custom colors to each row group (model layer)
  tab_style(style = list(
    cell_text(align = "center", weight = "bold"),
    cell_fill(color = model_colors["Change in MAOC stocks in mineral layer (~35–55 cm)"])),
    locations = cells_row_groups(groups = "Change in MAOC stocks in mineral layer (~35–55 cm)")) %>%
  tab_spanner(label = "", columns = everything()) 

final_table_gt %>% 
  tab_source_note(source_note = paste0(
    "Model fit: ",
    "Mineral (35–55 cm): ", r2_fin_text_m))

#gtsave(final_table_gt , filename = paste0(manDir, "figures/differenceModelFin.png"))


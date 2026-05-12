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

#Create dataframe with soil samples by core, and calculate the cumulative sum of properties for stock calculations and equivalent ash method
core <- soilDf %>%            
  group_by(year, block, fence, plot, treatment) %>% 
  mutate(coreNum = cur_group_id(), 
         cu.soil.stock = cumsum(soil.stock),
         cu.ash.stock = cumsum(ash.stock),
         cu.C.stock = cumsum(C.stock),
         cu.N.stock = cumsum(N.stock), 
         cu.ice.stock = cumsum(moisture), 
         cu.Fe.stock = cumsum(Fe.stock),
         cu.FeDCB.stock = cumsum(FeDCB.stock),
         cu.LC.stock = cumsum(stockLC), 
         cu.HC.stock = cumsum(stockHC), 
         cu.LN.stock = cumsum(stockLN),
         cu.HN.stock = cumsum(stockHN),
         bulkAsh = bulk.density*ash) %>% 
  ungroup()

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
ggsave(paste0(manDir, "figures/suppHFracCDepth", Sys.Date(), ".png"), width = 15, height = 20, units = "cm", dpi = 300)

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
ggsave(paste0(manDir, "figures/suppLightStocksDepth", Sys.Date(), ".png"), width = 10, height = 15, units = "cm", dpi = 300)

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
ggsave(paste0(manDir, "figures/stocksDepth", Sys.Date(), ".png"), width = 27, height = 20, units = "cm", dpi = 300)

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

ggsave(paste0(manDir, "figures/c14BarSupp", Sys.Date(), ".png"), width = 28, height = 15, units = "cm", dpi = 300)

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

ggsave(paste0(manDir, "figures/suppFrac14CByTreat", Sys.Date(), ".png"), width = 20, height = 15, units = "cm", dpi = 300)

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
ggsave(paste0(manDir, "figures/suppAgeOffset", Sys.Date(), ".png"), width = 20, height = 15, units = "cm", dpi = 300)


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
ggsave(paste0(manDir, "figures/ageOffsetBar", Sys.Date(), ".png"), width = 20, height = 15, units = "cm", dpi = 300)


##                    ##
##      ALL SPECIES   ##
##                    ##
##        ACi fits    ##
##                    ##
##   Fits & Params    ##
##                    ##

# This script analyses all species together. It reads of a 
# master dataset created for data repository.
# I have found some differences with fits from individual scripts.
# Let's continue with this script, because we are removing some datapoints.
#
# In here: 
#   1. removing extra 400CO2 points in monotonic and split protocols.
#   2. analyse all data
#   3. test the effect of the later leg on split protocols
#
# I think this could become the only script with all the analyses.
#
# I'm removing a lot of unnecessary lines... go find them in the individual scripts.

rm(list= ls(all.names = T))

require(dplyr)
require(ggplot2)
require(emmeans)
require(lme4)

myTheme <- {theme(panel.background = element_rect(fill = NA, color = "black"),
                  panel.grid = element_blank(),
                  panel.spacing.x = unit(0, "cm"),
                  strip.background = element_rect(fill = NA, color = "black"),
                  strip.placement = "outside",
                  strip.text = element_text(size = 8, color = "black"),
                  legend.background = element_blank(),
                  legend.title = element_blank(),
                  legend.key = element_blank(),
                  axis.ticks.length = unit(.25, "cm"),
                  axis.ticks = element_line(lineend = "round"), 
                  text = element_text(size = 12), 
                  axis.text.x = element_text(size = 12, color = "black", 
                                             margin = margin(t = 7)),
                  axis.text.x.top = element_text(size = 12, color = "black", 
                                                 margin = margin(b = 7)),
                  axis.text.y = element_text(size = 12, color = "black", 
                                             margin = margin(r = 7))
)}
source("//es.msu.edu/prl/labs/walker/General Project Storage/Mauri TN/Useful R-Codes/Functions/axisLabels.R")
ACiCol <- c('#d7191c','#fdae61','grey35','#abd9e9','#2c7bb6')

## I. Load Data ----
# setwd("C:/Users/mauri/Desktop/ACi Protocols")
noDAT_CO2levels <- c(50, 100, 150, 200, 250, 300, 350, 400, 420, 500, 600, 700, 800, 1000, 1250, 1500, 1750)
errCO2 = .12

AllACi.df <- read.csv("//es.msu.edu/prl/labs/walker/General Project Storage/Mauri TN/Data&Analyses/Collaborations/ACi Methods/ACiMethod_4.0 RepositoryData.csv") %>% 
  mutate(RunID = ifelse(Species == "Tobacco", mapply(paste, Species, "Mch", RunID, sep = "_"), RunID),
         
     ## Make CO2_r discrete
         CO2level = ifelse(Protocol != "DAT", round(CO2_r),
                           cut(x = CO2_r,
                               breaks = matrix(c(noDAT_CO2levels - noDAT_CO2levels*errCO2, noDAT_CO2levels + noDAT_CO2levels*errCO2), ncol = 2) %>% t() %>% c,
                               labels = as.matrix(data.frame(noDAT_CO2levels, NA), ncol = 2) %>% t() %>% c %>% .[-34]) %>% as.character %>% as.numeric()
  )) %>% 
    ## Correct some CO2 levels (it's ugly but works)
  mutate(CO2level = ifelse(CO2level %in% 49:51, 50,
                           ifelse(CO2level %in% 99:101, 100, 
                                  ifelse(CO2level == 199, 200, 
                                         ifelse(CO2level == 303, 300,
                                                ifelse(CO2level %in% 398:420, 400, 
                                                       ifelse(CO2level %in% 499:501, 500, 
                                                              ifelse(CO2level == 601, 600, 
                                                                     ifelse(CO2level == 799, 800, 
                                                                            ifelse(CO2level %in% 1200:1249, 1250, 
                                                                                   ifelse(CO2level == 1001, 1000, 
                                                                                          ifelse(CO2level == 1800, 1750, 
                                                                                                 CO2level)))))))))))) 


## Check Protocol structure
AllACi.df %>% with(table(Protocol, Species))
  ## Looks good

## _i) Data Clean-up ----

## __ *Apple ----
## Remove extra DATs done in the field, and an extra SplitUpDown mistakenly done instead of a UpDown
AllACi.df <- AllACi.df %>% 
  filter(!(RunID == "Apple_sch1_2_DAT" & Elapsed > 10500) &
           !(RunID == "Apple_wlk1_2_DAT" & Elapsed > 18000) &
           !(RunID == "Apple_wlk4_2_DAT" & Elapsed > 17500) &
           !(RunID == "Apple_wlk4_1_SplitDownUp" & Elapsed < 2500))

## __ *Remove 400 CO2 ----
## The Split Protocols have a first and last 400, and 3 400 in the middle, only keep the center CO2 in the
## middle of the sequence. Also, remove the first 400 in the monotonic at the beginning of the sequence.
## For Splits, keep 3th 400 in all species, and 4th 400 in Tobacco (It has 2 400 at the beginning)
AllACi.df <- AllACi.df  %>% 
  group_by(RunID) %>% 
  mutate(keep400_01 = ifelse(CO2level == 400, 1, 0),
         keep400 = cumsum(CO2level == 400)*keep400_01,
         keep = ifelse(keep400 == 0 | Protocol == "DAT" | Protocol == "Random", "Y",
                       ifelse((Protocol %in% c("SplitDownUp", "SplitUpDown") & 
                                Species == "Tobacco" & keep400 == 4) |
                                (Protocol %in% c("SplitDownUp", "SplitUpDown") & 
                                   Species != "Tobacco" & keep400 == 3) |
                                (Protocol %in% c("Monotonic") & keep400 == 2), "Y", "N"))) %>% 
  filter(keep == "Y") # %>%  
  # filter(Species == "Apple") %>% 
  # ggplot(aes(x = Elapsed, y = A, col = RunID)) + 
  # facet_grid(Protocol ~ Species)+
  # geom_point(data = . %>% filter(CO2level == 400),
  #            show.legend = F) + geom_line(show.legend = F)
  # ^ This works like a charm, 
  # use the plot to compared the plot with or without the filter and confirmed what we need.

##
## II. Anet comparison by CO2level ----
##

avgAllACi.df <-
  AllACi.df %>% 
  group_by(Species, Protocol, CO2level) %>%     
  summarize(#avgElapsed = mean(elapsed, na.rm = T), sdElapsed = sd(elapsed, na.rm = T),
    avgCO2_r = mean(CO2_r, na.rm = T), sdCO2_r = sd(CO2_r, na.rm = T),
    avgCi = mean(Ci, na.rm = T), seCi = sd(Ci, na.rm = T)/sqrt(length(Ci)),
    avgA = mean(A, na.rm = T),  seA = sd(A, na.rm = T)/sqrt(length(A))) 

## Check structure
avgAllACi.df %>% with(table(CO2level, Species))

## The analysis on Ci is just to have a mean value to where to plot the data points.
emm.Ci <- lm(Ci ~ as.factor(CO2level)*Species, data = AllACi.df) %>%
  emmeans(~ CO2level|Species) %>% multcomp::cld(Letters = letters, reversed = T) %>% 
  as.data.frame() %>% select(Species, CO2level, emmean) %>% dplyr::rename(avgCi = emmean)

lm(A ~ as.factor(CO2level)*Protocol + Ci, data = AllACi.df) %>% anova
emm.A <- AllACi.df %>% 
  split(., .$Species) %>% 
  lapply(function(spp.df){
    lm(A ~ as.factor(CO2level)*Protocol + Ci, data = spp.df) %>% 
      emmeans(~ Protocol|CO2level) %>%   multcomp::cld(Letters = letters, reversed = T)  %>% 
      as.data.frame() %>% 
      mutate(Species = spp.df$Species[1])
  }) %>% do.call(rbind,.)

## Create dataframe with signficant differences and coordinate for "*"
star.df <- emm.A %>% as.data.frame() %>% 
  filter(CO2level != 420 & !is.na(emmean)) %>% 
  group_by(Species, CO2level) %>% 
  mutate(star = ifelse(all(.group == ' a'), "", "*")) %>% 
  filter(Protocol == Protocol[emmean == max(emmean, na.rm = T)]) %>% as.data.frame() %>% 
  
  select(Species, CO2level, emmean, SE, star) %>% 
  left_join(emm.Ci) %>% 
  # mutate(Species = factor(Species, levels = c("Tobacco", "Arabidopsis", "Soybean", "Potato", "Apple"))) %>% 
  filter(star == "*") %>% 
  mutate(emmean = ifelse(Species == "Tobacco" & avgCi > 500, 21, emmean))

## __i) Figure 1: Average ACi Plot ----
Fig1.AllACi.plt <-
  avgAllACi.df %>% 
  filter(Protocol != "DAT") %>%
  filter(avgA > 0) %>% 
  # mutate(Species = factor(Species, levels = c("Tobacco", "Arabidopsis", "Apple", "Soybean", "Potato"))) %>% 
  
  ggplot(aes(x = avgCi, y = avgA, col = Protocol)) +
  scale_color_manual(values = ACiCol[-1]) +
  facet_wrap(~Species, ncol = 3, scales = "free") +
    
  geom_errorbar(data = . %>% filter(Protocol != "DAT"),
                aes(ymin = avgA - seA, ymax = avgA + seA), width = 25, lwd = 1) +
  geom_errorbarh(aes(xmin = avgCi - seCi, xmax = avgCi + seCi), height = .75, lwd = 1) +

  geom_point(data = . %>% filter(Protocol != "DAT"), size = 2) +
  geom_line(data = . %>% filter(Protocol != "DAT" & !(Protocol == "SplitUpDown" & avgA > 17 & CO2level == 420)),
            lwd = 1) +
  geom_smooth(inherit.aes = F, data = AllACi.df %>% filter(Protocol == "DAT"), 
              aes(x = Ci, y = A), 
              col = ACiCol[1], lwd = 1, fill = ACiCol[1], alpha = .15) +
  
  geom_text(data = star.df,
            aes(x = avgCi, y = emmean + SE + .25, label = star),
            col = "black" , cex = 7) +
  
  scale_x_continuous(name = "Ci (umolCo2 mol)") +
  scale_y_continuous(name = axisLabel("Anet")) + expand_limits(y=0) + 
  myTheme + theme(legend.position = c(.8, .2),
                  axis.ticks.x.top = element_blank(),
                  axis.ticks.length.x.top = unit(0, "cm"))
Fig1.AllACi.plt

## __ii) Table 1: length of protocol ----
lngthACidf <- AllACi.df %>% group_by(Species, Protocol, RunID) %>% 
  summarize(# minTime = min(elapsed),
    # iniTime = range(elapsed)[1],
    # endTime = range(elapsed)[2],
    lngthTime = range(Elapsed) %>% diff) %>% as.data.frame() %>% 
  mutate(lngthTime = ifelse(Protocol  == "DAT", lngthTime + 5*60, lngthTime))
      ## ^ 5min Add pre-ramp wait time

lngth.emm <- lm(lngthTime/60 ~ Species*Protocol, data = lngthACidf) %>% 
  # anova
  emmeans(~Protocol|Species) %>%
  multcomp::cld(Letters = letters, reversed = F) %>% 
  as.data.frame() %>% 
  mutate(Param = "Length",
         .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)})

## 30.7 / 14.7 = 2x faster
## 30.7 - 14.7  = 16 min faster

allSpp.length.tbl <- lngth.emm %>% 
  mutate(Species, Protocol,
            emmean = round(emmean,1),
            SE = round(SE, 2),
            .group = .group) %>% 
  mutate(`Length (min)` = mapply(paste, emmean, SE, .group, collapse = "+-")) %>% 
  select(-(emmean:Param)) %>%
  tidyr::pivot_wider(id_cols = Protocol, values_from = `Length (min)`, names_from = Species)  #write.table("clipboard", row.names = F, sep = "\t")

allSpp.length.tbl
# Apl.length.tbl %>% write.table("clipboard", sep = "\t", row.names = F)

## __iii) Check DAT slope 
AllACi.df %>% with(table(Protocol, CO2level, Species))
AllACi.df %>% 
  group_by(Species, RunID) %>% 
  filter(Protocol == "DAT") %>% 
  group_by(Species, RunID) %>% 
  summarize(totTime = range(Elapsed) %>% diff,
            totTime_min = totTime/60,
            slope = diff(range(CO2_r))/(totTime/60))

##
## ---- III. ACi Fit (All Data) ----
## 

## __ i) Fit all ACi individually ----
require(msuRACiFit)

AllACi.df$RunID  <- AllACi.df$RunID %>% as.factor

AllACi.df %>% with(table(RunID, Species))
AllACiFits <- AllACi.df  %>% filter(!is.na(A)) %>% 
  mutate(Tleaf = ifelse(Tleaf > 100, 25, Tleaf)) %>% 
  split(., as.character(.$RunID)) %>% 
  lapply(function(x) {
    Tleaf = mean(x$Tleaf, na.rm = T)  
    fitACi(data = x, tleaf = Tleaf) %>% .[[1]] %>% as.data.frame() %>% mutate(Tleaf = Tleaf)
  }) %>% 
  do.call(rbind, .) %>%
  mutate(.before = VcMax, 
         RunID = rownames(.),
         Species = strsplit(RunID, split ="_") %>% do.call(rbind, .) %>% .[,1],
         Machine = strsplit(RunID, split ="_") %>% do.call(rbind, .) %>% .[,2],
         Plant = strsplit(RunID, split ="_") %>% do.call(rbind, .) %>% .[,3],
         Protocol = RunID %>% strsplit(split = "_") %>% do.call(rbind, .) %>% .[,4])  %>% 
  mutate(Protocol = ifelse(Species == "Tobacco", RunID %>% strsplit(split = "_") %>% do.call(rbind, .) %>% .[,3], Protocol),
         Plant = ifelse(Species == "Tobacco", RunID %>% strsplit(split = "_") %>% do.call(rbind, .) %>% .[,4], Plant),
         Species = ifelse(Species == "At", "Arabidopsis", Species))
AllACiFits %>% with(table(Species, Plant))

## This script is very smooth, I don't think there is a need
## to export dataset.
# write.csv(AllACiFits, "//es.msu.edu/prl/labs/walker/General Project Storage/Mauri TN/Data&Analyses/Collaborations/ACi Protocols/ACiParams/ACiParamsAllSpp.csv",
          # row.names = F)
# AllACiFits <- read.csv("//es.msu.edu/prl/labs/walker/General Project Storage/Mauri TN/Data&Analyses/Collaborations/ACi Protocols/ACiParams/ACiParamsAllSpp.csv")

## __ ii) Protocol effect on ACi params ----
AllACiFits.lng <- AllACiFits %>% 
  tidyr::pivot_longer(cols = VcMax:as,
                      values_to = "Value",
                      names_to = "Param")
AllACiFits.lng %>% head

## By Species and parameters independently
AllACiFits.emm <- 
  split(AllACiFits.lng, AllACiFits.lng$Species) %>% 
  lapply(function(spp.df){
    split(spp.df, spp.df$Param) %>% 
      lapply(function(Param) {
        tryCatch({
          lmer(Value ~ Protocol +  (1|Machine:Plant), data = Param)  %>% 
            emmeans(~Protocol) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
            as.data.frame() %>% 
            mutate(.before = Protocol,
                   Species = Param$Species[1],
                   Param = Param$Param[1],
                   .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)})
        }, error = function(e) data.frame(Species = Param$Species[1],
                                          Param = Param$Param[1],
                                          Protocol = c("DAT", "Monotonic", "SplitDownUp", "SplitUpDown"),
                                          emmean = mean(Param$Value), SE = NA, df = NA, 
                                          lower.CL = NA, upper.CL = NA, .group = NA) 
        )
      }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)

## Look at param comparisons
split(AllACiFits.lng, AllACiFits.lng$Param) %>% 
  lapply(function(Param) {
    tryCatch({
      
      lmer(Value ~ Protocol +  (1|Machine:Plant), # (1|Machine)
         data = Param)  %>% 
      emmeans(~Protocol) %>% pwpm()
    }, error = function(e) e)
  }) #%>% do.call(rbind, .)

## __iii) Figure S3: Param Comparison ----
ann_text <- data.frame(emmean = 3, Param = 7, 
                       Species = factor(c("Tobacco", "Arabidopsis", "Soybean", "Potato", "Apple")),
                       panel = factor("Small",levels = c("Large", "Small"))) %>% 
  mutate(Label = Species)


ACiCol2 <- c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "grey35")
FigS3.ACiparm.plt <-
  AllACiFits.emm %>% 
  mutate(panel = ifelse(Param %in% c("VcMax","J"), "Large", "Small") %>%  factor(levels = c("Small", "Large")),
         Param = factor(Param, levels = c("VcMax","J","TPU","gm","rL","ag","as") %>% rev),
         xPos = ifelse(Param %in% c("rL", "ag", "as"),  emmean + 1,  emmean + SE + emmean*.05)) %>% 
  
  ggplot(aes(x = emmean, y = Param, col = Protocol)) +
  scale_color_manual(values = ACiCol2) +
  facet_wrap(~Species*panel, scales = "free_x", ncol = 4) +
  
  geom_point(position = position_dodge(width = .75)) +
  geom_errorbarh(aes(xmin = emmean - SE, xmax = emmean + SE), 
                 position = position_dodge(width = .75)) +
  
  geom_text(aes(label = .group, x = xPos, group = Protocol),
            position = position_dodge(width = .75),
            show.legend = F, col = "black", cex = 2.5) +
  
  geom_text(data = ann_text, aes(label = Label),
            col = "black") +
  
  geom_hline(yintercept = c(Inf, -Inf)) +
  geom_vline(data = data.frame(panel = c("Small", "Large") %>%  factor(levels = c("Small", "Large")),
                               xInt = c(-Inf, Inf)),
             aes(xintercept = xInt)) +
  geom_vline(xintercept = Inf, lty = 2, col = "grey20")+
  
  scale_y_discrete(name = NULL, label = c("VcMax\n(umol m-2 s-1)", "J\n(umol m-2 s-1)", "TPU\n(umol m-2 s-1)", 
                                          "gm\n(umol m-2 s-1 pa-1)", "rL\n(umol m-2 s-1)", "ag", "as") %>% rev) +
  scale_x_continuous(name = "Parameter value")  +

  myTheme +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0.0, "cm"),
        legend.position = c(.9, .15)
  )
FigS3.ACiparm.plt

## __iv) Steady ~ DAT Correlation ----
AllACiFits %>% head
corrACiparam.wd <- left_join(
  AllACiFits %>% filter(Protocol != "DAT")  %>% `names<-`(c(names(.)[1:5], paste(names(.)[6:13], "Steady", sep = "_"))),
  AllACiFits %>% filter(Protocol == "DAT") %>% `names<-`(c(names(.)[1:5], paste(names(.)[6:13], "DAT", sep = "_"))) %>% 
    select(-Protocol, -RunID ) 
) %>% filter(Protocol != "Random")

corrACiparam.lng <- corrACiparam.wd %>%  
  tidyr::pivot_longer(
    cols = c(VcMax_Steady:Tleaf_DAT),
    names_to = "Param",
    values_to = "Value"
  ) %>%
  tidyr::separate(Param, into = c("Param", "Type"), sep = "_") %>% 
  tidyr::pivot_wider(
    names_from = Type,
    values_from = Value
  )
## ___* Deming Regression (doi:10.21105/joss.04148) ----
# library(SimplyAgree)

Deming.Param.df <- corrACiparam.lng %>% 
  as.data.frame() %>% 
  split(., .$Param) %>% 
  lapply(function(x) {
    SimplyAgree::dem_reg(y = "Steady", x = "DAT", id = "Species", conf.level = .90,
            data = x)  %>% .$model %>% 
      mutate(Int.slope = row.names(.),
             Signif = ifelse(Int.slope == "Intercept", 
                             ifelse(lower.ci*upper.ci > 0, "*", "NS"),
                             ifelse(lower.ci > 1 | upper.ci < 1 , "*", "NS"))) %>% 
      # select(Int.slope, coef, se, Signif) %>% 
      cbind(Param = x$Param[1], .)
  }) %>% do.call(rbind,.)

## Create intercept and slope dataframe for geom_abline()
reg.coef.df <- Deming.Param.df %>% select(Param, Int.slope, coef) %>% 
  mutate(Param = factor(Param, levels = c("VcMax", "J", "TPU", "gm", "rL", "ag", "as", "Tleaf"))) %>% 
  tidyr::pivot_wider(id_cols = Param, names_from = "Int.slope", values_from = "coef")

## Create the text of  equation for plotting
Reg.txt <-Deming.Param.df %>% 
  mutate(param.txt = mapply(paste, round(coef, 2), Signif, sep = ""),
         Int.slope = factor(Int.slope, levels = c("Slope", "Intercept"))) %>% 
  mutate(Param = factor(Param, levels = c("VcMax", "J", "TPU", "gm", "rL", "ag", "as", "Tleaf"))) %>% 
  arrange(Param, Int.slope) %>% 
    ## this is important so then they are pasted in the right order.
  group_by(Param) %>% 
  summarize(Func.txt = paste0(param.txt, collapse = "x + ") %>% paste("y =",.)) %>% 
  mutate(xpos = c(100, 150, 10 , 15, 2, 0.25, .3, 25),
         ypos = c(60, 120, 8, 7, 6, .05, .3, 23))

## create dataframe to set limits of facetted plots (regression plots are better as squares)
blank.df <- expand.grid(Param = c("VcMax", "J", "TPU"),
                        Limits = c("Low", "Up")) %>% 
  mutate(Param = factor(Param, levels = c("VcMax", "J", "TPU"))) %>% 
  arrange(Param) %>%
  mutate(X = c(25, 250, 75, 375, 5, 25),
         Y = X)
blank.df

## ___* Figure 2: Vcmax, J and TPU correlations ----
Fig2.Param.corr.plt <- corrACiparam.lng %>% 
  filter(Protocol != "DAT") %>% 
  mutate(Param = factor(Param, levels = c("VcMax", "J", "TPU", "gm", "rL", "ag", "as", "Tleaf")),
         Protocol = factor(Protocol, levels = c("Monotonic", "SplitDownUp", "SplitUpDown", "Overall")),
  ) %>% 
  filter(Param %in% c("VcMax", "J", "TPU")) %>%
  filter(!is.na(Protocol)) %>% 
  ## ^ THis are Random regime observations
  
  ggplot(aes( x = DAT, y = Steady, col = Protocol))+
  facet_wrap(~Param, scales = "free", ncol = 2, strip.position = "bottom")  +
  scale_color_manual(values = c(ACiCol[c(2, 4:5)], "black")) +
  geom_abline(intercept = 0, slope = 1) +
  
  geom_text(inherit.aes = T,
            data = Reg.txt %>% filter(Param %in% c("VcMax", "J", "TPU")),
            aes(x = xpos, y = ypos, label = Func.txt), col = "black", show.legend = F, hjust = 0,
            cex = 3) +

  geom_point(aes(pch = Species), size = 2, show.legend = T) +
  geom_abline(data = reg.coef.df %>% filter(Param %in% c("VcMax", "J", "TPU")), 
              aes(intercept = Intercept, slope = Slope),
              col = 'black', show.legend = F, size = 1.5) +
  # geom_smooth(method = "lm", se = F , col = 'black', show.legend = F, size = 1.5) +
      ## ^In case we need to show the differences between the two methods.
  
  
  geom_blank(data = blank.df, inherit.aes = F,
             aes(x = X, y = Y)) +
  ## ^this is a cool trick to set limits for each facet.
  
  theme_bw() +   
  theme(panel.grid = element_blank(), strip.placement = "outside",  strip.background = element_blank(),
        legend.position = c(.7, .25), legend.title = element_blank())
Fig2.Param.corr.plt

## ___* Figure S4: gm, rL, ag, as correlations ----
FigS4.Param.corr.plt <- corrACiparam.lng %>% 
  filter(Protocol != "DAT") %>% 
  mutate(Param = factor(Param, levels = c("VcMax", "J", "TPU", "gm", "rL", "ag", "as", "Tleaf")),
         Protocol = factor(Protocol, levels = c("Monotonic", "SplitDownUp", "SplitUpDown", "Overall")),
  ) %>% # 
  filter(!(Param %in% c("VcMax", "J", "TPU", "Tleaf"))) %>%
  filter(!is.na(Protocol)) %>% 
  ## ^ THis are Random regime observations
  
  ggplot(aes( x = DAT, y = Steady, col = Protocol))+
  facet_wrap(~Param, scales = "free", ncol = 2, strip.position = "bottom")  +
  scale_color_manual(values = c(ACiCol[c(2, 4:5)], "black")) +
  geom_abline(intercept = 0, slope = 1) +
  
  geom_text(inherit.aes = T,
            data = Reg.txt %>% filter(!(Param %in% c("VcMax", "J", "TPU", "Tleaf"))),
            aes(x = xpos, y = ypos, label = Func.txt), col = "black", show.legend = F, hjust = 0,
            cex = 3) +
  
  geom_abline(data = reg.coef.df %>% filter(!(Param %in% c("VcMax", "J", "TPU", "Tleaf"))), 
              aes(intercept = Intercept, slope = Slope),
              col = 'black', show.legend = F, size = 1.5) +
  # geom_smooth(method = "lm", se = F , col = 'black', show.legend = F, size = 1.5) +
  ## ^In case we need to show the differences between the two methods.
  
  geom_point(aes(pch = Species), size = 2, show.legend = T) +
  
  geom_blank(data = blank.df %>% filter(!(Param %in% c("VcMax", "J", "TPU", "Tleaf"))), inherit.aes = F,
             aes(x = X, y = Y)) +
  ## ^this is a cool trick to set limits for each facet.
  
  theme_bw() +   
  theme(panel.grid = element_blank(), strip.placement = "outside",  strip.background = element_blank(),
        legend.title = element_blank())
FigS4.Param.corr.plt

##
## ---- IV. ACi Fit (w.o Last leg) ----
## 

## Remove third leg of split Protocols.
AllACi.df %>% 
  filter(Protocol %in% c("SplitDownUp", "SplitUpDown")) %>% 
  ggplot(aes(x = Elapsed, y = CO2level, col = RunID)) +
  facet_grid(Species ~ Protocol) +
  geom_point(show.legend = F) + geom_line(show.legend = F)

AllACi.df %>% 
  group_by(Species, Protocol, RunID) %>% 
  mutate(keepSplit_01 = if_else((Protocol == "SplitDownUp" & Elapsed > Elapsed[CO2level == max(CO2level)]) |
                                 (Protocol == "SplitUpDown" & Elapsed > Elapsed[CO2level == min(CO2level)]), "N", "Y")) %>% 
  
  filter(Protocol %in% c("SplitDownUp", "SplitUpDown")) %>%
  
  ggplot(aes(x = Elapsed, y = CO2level, col = RunID)) +
  facet_grid(Species ~ Protocol) +
  # geom_point(show.legend = F) + 
  geom_line(aes(lty = keepSplit_01 %>% factor), show.legend = F)
    
## Time difference between two- and three-phase
AllACi.df %>% 
  group_by(Species, Protocol, RunID) %>% 
  mutate(keepSplit_01 = if_else((Protocol == "SplitDownUp" & Elapsed > Elapsed[CO2level == max(CO2level)]) |
                                  (Protocol == "SplitUpDown" & Elapsed > Elapsed[CO2level == min(CO2level)]), "N", "Y")) %>% 
  filter(Protocol %in% c("SplitDownUp", "SplitUpDown")) %>%
  group_by(Protocol, RunID) %>% summarise(length = diff(range(Elapsed))/60) %>% 
  group_by(Protocol) %>% summarise(length = mean(length))


AllACi.df %>% 
  group_by(Species, Protocol, RunID) %>% 
  mutate(keepSplit_01 = if_else((Protocol == "SplitDownUp" & Elapsed > Elapsed[CO2level == max(CO2level)]) |
                                  (Protocol == "SplitUpDown" & Elapsed > Elapsed[CO2level == min(CO2level)]), "N", "Y")) %>% 
  filter(Protocol %in% c("SplitDownUp", "SplitUpDown")) %>%
  filter(keepSplit_01 == "Y") %>% 
  group_by(Protocol, RunID) %>% summarise(length = diff(range(Elapsed))/60) %>% 
  group_by(Protocol) %>% summarise(length = mean(length))


AllACi.NoAmb.df <- AllACi.df %>% 
  group_by(Species, Protocol, RunID) %>% 
  mutate(keepSplit_01 = if_else((Protocol == "SplitDownUp" & Elapsed > Elapsed[CO2level == max(CO2level)]) |
                                  (Protocol == "SplitUpDown" & Elapsed > Elapsed[CO2level == min(CO2level)]), "N", "Y")) %>% 
  filter(keepSplit_01 == "Y")

AllACi.NoAmb.df %>%
  filter(Protocol %in% c("SplitDownUp", "SplitUpDown")) %>%
  ggplot(aes(x = Elapsed, y = CO2level, col = RunID)) +
  facet_grid(Species ~ Protocol) +
  # geom_point(show.legend = F) + 
  geom_line(aes(lty = keepSplit_01 %>% factor), show.legend = F)
  ## ^ Looks great!!
  ## Well done!

## __ i) Fit all ACi individually ----
AllACi.NoAmb.df$RunID  <- AllACi.NoAmb.df$RunID %>% as.factor
AllACi.NoAmb.df %>% with(table(RunID, Species))
AllACi.NoAmb.Fits <- AllACi.NoAmb.df  %>% filter(!is.na(A)) %>% 
  mutate(Tleaf = ifelse(Tleaf > 100, 25, Tleaf)) %>% 
  split(., as.character(.$RunID)) %>% 
  lapply(function(x) {
    Tleaf = mean(x$Tleaf, na.rm = T)  
    fitACi(data = x, tleaf = Tleaf) %>% .[[1]] %>% as.data.frame() %>% mutate(Tleaf = Tleaf)
  }) %>% 
  do.call(rbind, .) %>% 
  mutate(.before = VcMax, 
         RunID = rownames(.),
         Species = strsplit(RunID, split ="_") %>% do.call(rbind, .) %>% .[,1],
         Machine = strsplit(RunID, split ="_") %>% do.call(rbind, .) %>% .[,2],
         Plant = strsplit(RunID, split ="_") %>% do.call(rbind, .) %>% .[,3],
         Protocol = RunID %>% strsplit(split = "_") %>% do.call(rbind, .) %>% .[,4])  %>% 
  mutate(Protocol = ifelse(Species == "Tobacco", RunID %>% strsplit(split = "_") %>% do.call(rbind, .) %>% .[,3], Protocol),
         Plant = ifelse(Species == "Tobacco", RunID %>% strsplit(split = "_") %>% do.call(rbind, .) %>% .[,4], Plant),
         Species = ifelse(Species == "At", "Arabidopsis", Species))
AllACi.NoAmb.Fits %>% with(table(Species, Plant))

# write.csv(AllACi.NoAmb.Fits, "//es.msu.edu/prl/labs/walker/General Project Storage/Mauri TN/Data&Analyses/Collaborations/ACi Protocols/ACiParams/ACiParamsAllSpp.csv",
# row.names = F)
# AllACi.NoAmb.Fits <- read.csv("//es.msu.edu/prl/labs/walker/General Project Storage/Mauri TN/Data&Analyses/Collaborations/ACi Protocols/ACiParams/ACiParamsAllSpp.csv")

## __ ii) Protocol effect on ACi params ----
AllACi.NoAmb.Fits.lng <- AllACi.NoAmb.Fits %>% 
  tidyr::pivot_longer(cols = VcMax:as,
                      values_to = "Value",
                      names_to = "Param")
AllACi.NoAmb.Fits.lng %>% head

AllACi.NoAmb.Fits.emm <- 
  split(AllACi.NoAmb.Fits.lng, AllACi.NoAmb.Fits.lng$Species) %>% 
  lapply(function(spp.df){
    split(spp.df, spp.df$Param) %>% 
      lapply(function(Param) {
        tryCatch({
          lmer(Value ~ Protocol +  (1|Machine:Plant), data = Param)  %>% 
            emmeans(~Protocol) %>% multcomp::cld(Letters = letters, reversed = TRUE) %>% 
            as.data.frame() %>% 
            mutate(.before = Protocol,
                   Species = Param$Species[1],
                   Param = Param$Param[1],
                   .group = if(all(.group == " a")) {""} else {gsub(" ", "", .group)})
        }, error = function(e) data.frame(Species = Param$Species[1],
                                          Param = Param$Param[1],
                                          Protocol = c("DAT", "Monotonic", "SplitDownUp", "SplitUpDown"),
                                          emmean = mean(Param$Value), SE = NA, df = NA, 
                                          lower.CL = NA, upper.CL = NA, .group = NA) 
        )
      }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)

## __iii) No Show: Param Comparison ----
ann_text <- data.frame(emmean = 3, Param = 7, 
                       Species = factor(c("Tobacco", "Arabidopsis", "Soybean", "Potato", "Apple"), 
                                        levels = c("Tobacco", "Arabidopsis", "Soybean", "Potato", "Apple")),
                       panel = factor("Small",levels = c("Large", "Small"))) %>% 
  mutate(Label = Species)

## This is a cool graph, but I won't be published
emmAllACi.NoAmb.Fits.plt <-
  AllACi.NoAmb.Fits.emm %>% 
  mutate(panel = ifelse(Param %in% c("VcMax","J"), "Large", "Small") %>%  factor(levels = c("Small", "Large")),
         Param = factor(Param, levels = c("VcMax","J","TPU","gm","rL","ag","as") %>% rev),
         xPos = ifelse(Param %in% c("rL", "ag", "as"),  emmean + 1,  emmean + SE + emmean*.05)) %>% 
  
  ggplot(aes(x = emmean, y = Param, col = Protocol)) +
  scale_color_manual(values = ACiCol) +
  facet_wrap(~Species*panel, scales = "free_x", ncol = 6) +
  
  geom_point(position = position_dodge(width = .75)) +
  geom_errorbarh(aes(xmin = emmean - SE, xmax = emmean + SE), 
                 position = position_dodge(width = .75)) +
  
  geom_text(aes(label = .group, x = xPos, group = Protocol),
            position = position_dodge(width = .75),
            show.legend = F, col = "black", cex = 2.5) +
  
  geom_text(data = ann_text, aes(label = Label),
            col = "black") +
  
  geom_hline(yintercept = c(Inf, -Inf)) +
  geom_vline(data = data.frame(panel = c("Small", "Large") %>%  factor(levels = c("Small", "Large")),
                               xInt = c(-Inf, Inf)),
             aes(xintercept = xInt)) +
  geom_vline(xintercept = Inf, lty = 2, col = "grey20")+
  
  scale_y_discrete(name = NULL, label = c("VcMax\n(umol m-2 s-1)", "J\n(umol m-2 s-1)", "TPU\n(umol m-2 s-1)", 
                                          "gm\n(umol m-2 s-1 pa-1)", "rL\n(umol m-2 s-1)", "ag", "as") %>% rev) +
  scale_x_continuous(name = "Parameter value")  +
  
  myTheme +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0.0, "cm"),
        legend.position = c(.9, .15)
  )

emmAllACi.NoAmb.Fits.plt

## __iv) Amb ~ NoAmb Correlation ----
AllACi.NoAmb.Fits %>% head
names(AllACi.NoAmb.Fits)[6:13] <- paste("NoAmb", names(AllACi.NoAmb.Fits)[6:13], sep = "_")
names(AllACiFits)[6:13] <- paste("Amb", names(AllACiFits)[6:13], sep = "_")

Amb.NoAmb.df <-
  left_join(AllACi.NoAmb.Fits %>% select(-Machine, -Plant), 
            AllACiFits %>% select(-Machine, -Plant))  %>% 
  tidyr::pivot_longer(
    cols = c(NoAmb_VcMax:NoAmb_Tleaf, Amb_VcMax:Amb_Tleaf),
    names_to = c("Group", "Param"),  # Two new columns for group (Spp or All) and parameter name
    names_pattern = "^(NoAmb|Amb)_(.*)$",  # Regular expression to extract group and parameter
    values_to = "Value"
  ) %>% 
  tidyr::pivot_wider(id_cols = c(RunID:Protocol, Param), names_from = Group,values_from = Value)

## ___* Deming Regression (doi:10.21105/joss.04148) ----
# library(SimplyAgree)
NoAmb_Deming.Param.df <- 
  Amb.NoAmb.df %>% 
  filter(Protocol %in% c("SplitDownUp", "SplitUpDown")) %>% 
  split(., .$Protocol) %>% 
  lapply(function(Prot.df) {
    split(Prot.df, Prot.df$Param) %>% 
      lapply(function(x) {
        SimplyAgree::dem_reg(y = "NoAmb", x = "Amb", id = "Species",
                             data = x)  %>% .$model %>% 
          mutate(Int.slope = row.names(.),
                 Signif = ifelse(Int.slope == "Intercept", 
                                 ifelse(lower.ci*upper.ci > 0, "*", "NS"),
                                 ifelse(lower.ci > 1 | upper.ci < 1 , "*", "NS"))) %>% 
          select(Int.slope, coef, se, Signif) %>%
          cbind(Protocol = x$Protocol[1],
                Param = x$Param[1], .)
      }) %>% do.call(rbind,.)
  }) %>% do.call(rbind,.)
 
## Create intercept and slope dataframe for geom_abline()
reg.coef.df <- NoAmb_Deming.Param.df %>% select(Protocol, Param, Int.slope, coef) %>% 
  mutate(Param = factor(Param, levels = c("VcMax", "J", "TPU", "gm", "rL", "ag", "as", "Tleaf"))) %>% 
  tidyr::pivot_wider(id_cols = Protocol:Param, names_from = "Int.slope", values_from = "coef")

## Create the text of  equation for plotting
Reg.txt <- NoAmb_Deming.Param.df %>% 
  mutate(param.txt = mapply(paste, round(coef, 2), Signif, sep = ""),
         Int.slope = factor(Int.slope, levels = c("Slope", "Intercept"))) %>% 
  mutate(Param = factor(Param, levels = c("VcMax", "J", "TPU", "gm", "rL", "ag", "as", "Tleaf"))) %>% 
  arrange(Protocol, Param, Int.slope) %>% 
  ## this is important so then they are pasted in the right order.
  group_by(Protocol, Param) %>% 
  summarize(Func.txt = paste0(param.txt, collapse = "x + ") %>% paste("y =",.)) %>% 
  as.data.frame() %>% 
  mutate(xpos = c(100, 150, 10 , 15, 0, 0.15, .0, 25, 
                  100, 150, 10 , 15, 0, 0.15, .0, 25),
         ypos = c(70, 130, 8, 7, 6, .6, .7, 23, 
                  50, 100, 6, 5, 4, .5, .5, 21))


## Create blank dataset to fix the limits of each indivudal facet
# blank.df <- expand.grid(Parameter = c("VcMax", "J", "TPU"),
#                         Limits = c("Low", "Up")) %>% 
#   mutate(Parameter = factor(Parameter, levels = c("VcMax", "J", "TPU"))) %>% 
#   arrange(Parameter) %>%
#   mutate(X = c(25, 250, 75, 375, 5, 25),
#          Y = X)
# blank.df

## ___* Figure S5: Amb ~ NoAmb Correlation ----
FigS5.AmbNoAmb.corr.plt <- Amb.NoAmb.df %>% 
  mutate(Param = factor(Param, levels = c("VcMax", "J", "TPU", "gm", "rL", "ag", "as", "Tleaf"))) %>% 
  filter(Protocol %in% c("SplitDownUp", "SplitUpDown")) %>% 
  # filter(Param %in% c("VcMax", "J", "TPU")) %>%
  filter(Param != "Tleaf") %>%
  # mutate(Param = factor(Param, levels = c("VcMax", "J", "TPU"))) %>% 
  
  ggplot(aes(x = Amb, y = NoAmb, col = Protocol)) +
  facet_wrap(~ Param, scales = "free", ncol = 2) +
  scale_color_manual(values = ACiCol[4:5]) +
  
  geom_abline(intercept = 0, slope = 1) +

  geom_abline(data = reg.coef.df %>% filter(Param != "Tleaf"),# %>% filter(Param %in% c("VcMax", "J", "TPU")), 
              aes(intercept = Intercept, slope = Slope, col = Protocol),
              show.legend = F, size = 1.5) +
  # geom_smooth(method = "lm", lty =2, se = F) +
  geom_point() +
  
  geom_text(inherit.aes = T,
            data = Reg.txt %>% filter(Param != "Tleaf") ,# %>% filter(Param %in% c("VcMax", "J", "TPU")),
            aes(x = xpos, y = ypos, label = Func.txt, col = Protocol), show.legend = F, hjust = 0,
            cex = 3) +
  
  # geom_blank(data = blank.df, inherit.aes = F,
  #            aes(x = X, y = Y)) +
    ## ^this is a cool trick to set limits for each facet.
    
  scale_x_continuous(name = "Three-phase Split") +
  scale_y_continuous(name = "Two-phase Split") +
  
  theme_bw() + theme(panel.grid = element_blank(),
                     legend.title = element_blank(),
                     legend.background = element_blank(),
                     legend.position = c(.765, .1))
    ## This is cool, if the last leg of an splitDown up is removed Vcmax is overestimated.

## __v) NoAmb ~ DAT Correlation ----
AllACi.NoAmb.Fits %>% head
corrNoAmbACiparam.wd  <-
    left_join(
      ## Merge DAT and steady only for Splits without the third leg(back to ambient)
      AllACi.NoAmb.Fits %>% filter(Protocol != "DAT")  %>% `names<-`(c(names(.)[1:5], paste(names(.)[6:13], "Steady", sep = "_"))),
      AllACi.NoAmb.Fits %>% filter(Protocol == "DAT") %>% `names<-`(c(names(.)[1:5], paste(names(.)[6:13], "DAT", sep = "_"))) %>% 
        select(-Protocol, -RunID ) 
    ) %>% filter(Protocol %in% c("SplitDownUp", "SplitUpDown")) %>% mutate(Phase = "TwoPhase") 
corrNoAmbACiparam.wd %>% head

corrNoAmbACiparam.lng <- 
  corrNoAmbACiparam.wd %>%  
  tidyr::pivot_longer(
    cols = c(NoAmb_VcMax_Steady:NoAmb_Tleaf_DAT),
    names_to = "Param",
    values_to = "Value"
  ) %>%
  tidyr::separate(Param, into = c(NA, "Param", "Type"), sep = "_") %>% 
  tidyr::pivot_wider(
    names_from = Type,
    values_from = Value, #values_to = "Steady"
  )

corr.AmbNoAmb.ACiparam.df <- 
  rbind(corrACiparam.lng %>% 
          filter(Protocol %in% c("SplitDownUp", "SplitUpDown")) %>%
          mutate(Phase = "ThreePhase") ,
        corrNoAmbACiparam.lng) %>% 
  mutate(Param = factor(Param, levels = c("VcMax", "J", "TPU", "gm", "rL", "ag", "as", "Tleaf")))
  

corr.AmbNoAmb.ACiparam.df %>% 
  with(table(Protocol, Phase))

## ___* Figure S6: NoAmb ~ DAT Correlation ----
# I'm not going to do the same in depth analysis.
FigS6.NoAmbParam.corr.plt <- corr.AmbNoAmb.ACiparam.df %>% 
  # filter(Param %in% c("VcMax", "J", "TPU")) %>%
  filter(Param != "Tleaf") %>%
  ggplot(aes( x = DAT, y = Steady, col = Protocol, lty = Phase))+
  facet_wrap(~Param, scales = "free", ncol = 2)  +
  scale_color_manual(values = ACiCol[4:5]) +
  geom_abline(intercept = 0, slope = 1) +
  
  # geom_text(inherit.aes = T,
  #           data = Reg.txt %>% filter(Param %in% c("VcMax", "J", "TPU")),
  #           aes(x = xpos, y = ypos, label = Func.txt), col = "black", show.legend = F, hjust = 0,
  #           cex = 3) +
  
  geom_point(aes(pch = Species), size = 2, show.legend = T) +
  # geom_abline(data = reg.coef.df %>% filter(Param %in% c("VcMax", "J", "TPU")), 
  #             aes(intercept = Intercept, slope = Slope),
  #             col = 'black', show.legend = F, size = 1.5) +
  geom_smooth(method = "lm", se = F , show.legend = T) +
  ## ^In case we need to show the differences between the two methods.
  
  # geom_blank(data = blank.df, inherit.aes = F,
  #            aes(x = X, y = Y)) +
  ## ^this is a cool trick to set limits for each facet.
  
  theme_bw() +   
  theme(panel.grid = element_blank(), legend.position = c(.7, .05),
        legend.title = element_blank())
FigS6.NoAmbParam.corr.plt
##
## ---- XX Create powerpoint with all figures and tables ----
##

# Create a new powerpoint document
require(officer)
require(rvg)

# Write the document to a file

read_pptx() %>% 
  
  add_slide("Blank", "Office Theme") %>%
  ph_with(dml(ggobj =  Fig1.AllACi.plt),
          location = ph_location(width = 7.5, height = 5.2)) %>% 
  
  add_slide("Blank", "Office Theme") %>%
  ph_with(dml(ggobj =  FigS3.ACiparm.plt),
          location = ph_location(width = 6.5, height = 8.3)) %>% 

  add_slide("Blank", "Office Theme") %>%
  ph_with(dml(ggobj = Fig2.Param.corr.plt),
          location = ph_location(width = 5.8, height = 5.8)) %>% 
  
  add_slide("Blank", "Office Theme") %>%
  ph_with(dml(ggobj = FigS4.Param.corr.plt),
          location = ph_location(width = 7.6, height = 6)) %>% 
  
  add_slide("Blank", "Office Theme") %>%
  ph_with(dml(ggobj = FigS5.AmbNoAmb.corr.plt),
          location = ph_location(width = 5.5, height = 8.5)) %>% 
  
  add_slide("Blank", "Office Theme") %>%
  ph_with(dml(ggobj = FigS6.NoAmbParam.corr.plt),
          location = ph_location(width = 5.5, height = 8.5)) %>% 

  print(target = "C:/Users/tejerani/Desktop/R_firtsfigures/ACiProtocols_FirstFigure.pptx")

### ---- X. Extra! Extra! ----
##
## III. Anet correlations across Protocols ----
##

# Need to select the Anet values from DAT that are within range of the 
# CO2 levels measured on the other Protocols
#
# This was done here, but then I moved it to the Compilation script
# is the cut() function

corrACi.plt <- avgAllACi.df %>% 
  select(Protocol, CO2level, avgA) %>%  
  tidyr::pivot_wider(id_cols = CO2level, 
                     names_from = Protocol,
                     values_from = avgA) %>% 
  GGally::ggpairs(columns = 2:6, 
                  # upper = NULL, mapping=ggplot2::aes(colour = CO2level %>% as.factor)
  )


## IV Hysteresis ----

Split.df <- AllACi.df %>% filter(Protocol %in% c("SplitDownUp", "SplitUpDown"))


Split.df %>%  with(table(Run, Protocol))

Split.df %>% #filter(Run == 1) %>% 
  select(Protocol, Run, elapsed, CO2_r, Ci, A) %>% 
  group_by(Run) %>% 
  mutate(iniCi = Ci, iniA = A,
         endCi = c(Ci[-1], NA), endA = c(A[-1], NA)) %>% 
  
  ggplot(aes(x = Ci, y = A, col = elapsed)) +
  facet_wrap(~Run) +
  # geom_point() +
  geom_segment(aes(x = iniCi, y = iniA, xend = endCi, yend = endA),
               arrow = arrow(length = unit(0.25,"cm")))


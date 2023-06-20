# PACKAGES ====================================================================
library(here)
library(ggplot2)
library(janitor)
library(ggplotify)
library(cowplot)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(psych)
library(lme4)
library(lmerTest)
library(performance)
library(lmPerm)
library(emmeans)
library(multcompView)
library(scales)
library(eulerr)

# INFO ========================================================================
# This R script is part of the statistical analysis published in
# Eckert et al. (2023) in the scientific journal Metabolomics
# To load the corresponding raw data use the following R code:
# install.packages(zen4R)
# library(zen4R)
# zen4R::download_zenodo(doi = "10.5281/zenodo.7687634",
#                        path = here(".")) # add target directory to save files
# # unzip downloaded files
# utils::unzip("Eckert_et_al_2023_RawData.zip")

# DATA ========================================================================
## compound data ==============================================================
compounds_WideFormat <- readr::read_table(here(".", # replace with working directory
                                                 "Eckert_et_al_2023_RawData_Compounds.txt"),
                                          na="NA",
                                          col_types=cols(unique_ID="c",
                                                         donor_lab="f",
                                                         recipient_lab="f",
                                                         chemotype="f",
                                                         treatment="f",
                                                         plant_ID="i",
                                                         FW_leaf_g="d",
                                                         .default=col_double())); compounds_WideFormat
# split by chemotype
compounds_WideFormat <- split(compounds_WideFormat,
                              compounds_WideFormat$chemotype); compounds_WideFormat

## peak data ==================================================================
peak_counts <- bind_rows(compounds_WideFormat) %>% 
  filter(chemotype != 'blank') %>%
  dplyr::select(-c(FW_leaf_g, bromodecane_1_EIC)) %>%
  mutate(peak_counts = as.integer(rowSums(.[,-c(1:6)] > 0)),
         .after=plant_ID) %>%
  .[ , !sapply(., is.double)] %>% 
  mutate(donor_lab = factor(donor_lab,
                            levels=c("L1", "L2", "L3", "L4","L5")),
         recipient_lab = factor(recipient_lab,
                                levels=c("L1","L2","L4","L5"))); peak_counts

## bromodecane data ===========================================================
bromodecane <- readr::read_table(here(".", # replace with working directory
                               "Eckert_et_al_2023_RawData_Bromodecane.txt"),
                          na="NA",
                          col_types=cols(unique_ID="c",
                                         donor_lab="f",
                                         recipient_lab="f",
                                         chemotype="f",
                                         treatment="f",
                                         plant_ID="i",
                                         FW_leaf_g="d",
                                         bromodecane_applied="l",
                                         bromodecane_1_EIC="i")); bromodecane
# remove blanks and samples where bromodecane was not applied
# normalize to headspace-sampled leaf fresh weight
# adjust order of factors donor laboratory and recipient laboratory
bromodecane <- bromodecane %>%
  filter(chemotype != "blank",
         bromodecane_applied == 1) %>% 
  mutate(bromodecane_n_normalized = 
           ifelse(is.na(bromodecane_1_EIC) == T, NA, bromodecane_1_EIC / FW_leaf_g),
         donor_lab = factor(donor_lab,
                            levels = c("L1", "L2", "L3", "L4", "L5")),
         recipient_lab = factor(recipient_lab,
                               levels = c("L1", "L2", "L4", "L5"))); bromodecane # 370 x 10
# count samples where bromodecane was detected
length(which(bromodecane$bromodecane_1_EIC != 0)) # 167
length(bromodecane$bromodecane_1_EIC) # 370
# split by chemotype
bromodecane <- split(bromodecane, bromodecane$chemotype)

## ICC data ===================================================================
# Mono-chemotype
# count obs per recipient lab
rlab_mono_obs <- bromodecane$mono %>%
  count(recipient_lab); rlab_mono_obs$n # 40-49
# add Sample_ID (sequence within donor laboratory)
bromodecane_mono_LongFormat <- bromodecane$mono %>%
  dplyr::select(donor_lab,
         recipient_lab,
         plant_ID,
         treatment,
         bromodecane_n_normalized) %>%
  relocate(recipient_lab, .before = donor_lab) %>% 
  arrange(recipient_lab, donor_lab, treatment, plant_ID) %>% 
  mutate(bromodecane_sample_ID = unlist(lapply(rlab_mono_obs$n, sequence))) %>% 
  dplyr::select(-one_of("plant_ID", "treatment")); bromodecane_mono_LongFormat
# format from long-format to wide-format
bromodecane_mono_WideFormat <- pivot_wider(data=bromodecane_mono_LongFormat[,-2], # remove donor_lab column
                                 names_from = recipient_lab,
                                 values_from = bromodecane_n_normalized); bromodecane_mono_WideFormat
# Mixed-chemotype
# count obs per recipient lab
rlab_mixed_obs <- bromodecane$mixed %>%
  count(recipient_lab); rlab_mixed_obs$n # 40-48
# add Sample_ID (sequence within donor laboratory)
bromodecane_mixed_LongFormat <- bromodecane$mixed %>%
  dplyr::select(donor_lab,
         recipient_lab,
         plant_ID,
         treatment,
         bromodecane_n_normalized) %>%
  relocate(recipient_lab, .before = donor_lab) %>% 
  arrange(recipient_lab, donor_lab, treatment, plant_ID) %>% 
  mutate(bromodecane_sample_ID = unlist(lapply(rlab_mixed_obs$n, sequence))) %>% 
  dplyr::select(-one_of("plant_ID", "treatment")); bromodecane_mixed_LongFormat
# format from long-format to wide-format
bromodecane_mixed_WideFormat <- pivot_wider(data=bromodecane_mixed_LongFormat[,-2], # remove donor_lab column
                                           names_from = recipient_lab,
                                           values_from = bromodecane_n_normalized); bromodecane_mixed_WideFormat

# STATISTICS ==================================================================
## Mann-whitney U test ========================================================
# prepare data
bromodecane_MWUtest_df <- bind_rows(bromodecane$mono,
                                    bromodecane$mixed) %>% 
  filter(bromodecane_applied==T) %>%
  split(f = as.factor(.$chemotype)); bromodecane_MWUtest_df
# test
# Mono-chemotype
wilcox.test(bromodecane_MWUtest_df$mono[bromodecane_MWUtest_df$mono$treatment == "C",
                                         "bromodecane_n_normalized"][[1]],
            bromodecane_MWUtest_df$mono[bromodecane_MWUtest_df$mono$treatment == "JA",
                                         "bromodecane_n_normalized"][[1]],
            paired=F)
# Mixed-chemotype
wilcox.test(bromodecane_MWUtest_df$mixed[bromodecane_MWUtest_df$mixed$treatment == "C",
                                        "bromodecane_n_normalized"][[1]],
            bromodecane_MWUtest_df$mixed[bromodecane_MWUtest_df$mixed$treatment == "JA",
                                        "bromodecane_n_normalized"][[1]],
            paired=F)

## ICC ========================================================================
# Mono-chemotype
set.seed(123); bromodecane_mono_ICC <- irr::icc(bromodecane_mono_WideFormat[,-1],
                                                model = "twoway",
                                                type = "agreement",
                                                unit = "single"); bromodecane_mono_ICC

# Mixed-chemotype
set.seed(123); bromodecane_mixed_ICC <- irr::icc(bromodecane_mixed_WideFormat[,-1],
                                                 model = "twoway",
                                                 type = "agreement",
                                                 unit = "single"); bromodecane_mixed_ICC

## LMPERM =====================================================================
# prepare data
bromodecane_freq <- bind_rows(bromodecane$mono,
                              bromodecane$mixed) %>%
  filter(bromodecane_applied == TRUE) %>%
  filter(is.na(bromodecane_1_EIC) == F) %>% 
  group_by(donor_lab, recipient_lab, chemotype) %>%
  mutate(max_count = n()) %>% 
  mutate(across(bromodecane_n_normalized, ~ sum(.x != 0))) %>%
  dplyr::select(donor_lab, recipient_lab, chemotype, bromodecane_n_normalized, max_count) %>% 
  group_by(donor_lab, recipient_lab, chemotype) %>%
  mutate(freq = bromodecane_n_normalized / max_count) %>%
  ungroup() %>%
  distinct(.keep_all = T); bromodecane_freq
# percentage of samples with bromodecane detected across all samples
bromodecane_freq %>% 
  summarize(sum_bromodecane=sum(bromodecane_n_normalized),
            sum_total=sum(max_count)) %>% 
  mutate(sum_freq = sum_bromodecane / sum_total * 100)
# calculate means and sd of percentage values per donor lab and chemotype
bromodecane_freq_descr <- bromodecane_freq %>%
  group_by(donor_lab, chemotype) %>%
  dplyr::summarize(mean=mean(freq),
                   sd=sd(freq)) %>% 
  ungroup(); bromodecane_freq_descr
# model with frequency values
set.seed(123); bromodecane_mod <- summary(aovp(freq ~ donor_lab +
                                                 chemotype,
                                               perm="Exact",
                                               data=bromodecane_freq)); bromodecane_mod

## LMM ========================================================================
### Donor laboratories ========================================================
# create baseline model
dlab_chem_mod0 <- lmer(peak_counts ~ donor_lab * chemotype +
                               (1|recipient_lab/plant_ID) +
                               (1|treatment),
                             REML=F,
                             data=peak_counts,
                             control = lmerControl(
                               optimizer ='Nelder_Mead')); summary(dlab_chem_mod0)
# check assumptions visually
check_model(dlab_chem_mod0)
# update model without interaction
dlab_chem_mod1 <- update(dlab_chem_mod0, .~. - donor_lab:chemotype)
# check assumptions visually
check_model(dlab_chem_mod1)
# LRT
anova(dlab_chem_mod0, dlab_chem_mod1)
# update model without factor donor laboratory
dlab_chem_mod2 <- update(dlab_chem_mod1, .~. -donor_lab)
# check assumptions visually
check_model(dlab_chem_mod2)
# LRT
dlab_mod_LRT0 <- anova(dlab_chem_mod1, dlab_chem_mod2); dlab_mod_LRT0
# update model without factor chemotype
dlab_chem_mod3 <- update(dlab_chem_mod1, .~. -chemotype)
# check assumptions visually
check_model(dlab_chem_mod3)
# LRT
dlab_mod_LRT1 <- anova(dlab_chem_mod1, dlab_chem_mod3); dlab_mod_LRT1
# LMM model (type II using Satterthwaite's method of p-values)
dlab_mod <- anova(dlab_chem_mod1, type=2); dlab_mod
# post-hoc testing
dlab_HSD <- emmeans(dlab_chem_mod1,
                    list(pairwise ~ donor_lab),
                    adjust = "tukey"); dlab_HSD$`pairwise differences of donor_lab`
# transform to complete p-value matrix
dlab_HSD_pvals <- dlab_HSD$`pairwise differences of donor_lab` %>% 
  data.frame() %>% 
  dplyr::select(any_of(c("X1", "p.value"))) %>% 
  separate(X1,
           into=c("V1","V2"),
           sep=" - ",
           convert=T,
           extra="merge") %>% 
  spread(V1, p.value, fill=0) %>%
  column_to_rownames(var = "V2") %>%
  mutate(L5 = c(.[4,][2:4], 0)) %>% 
  rbind(list(0, .[1,1], .[2,1], .[3,1], .[4,1])) %>%
  magrittr::set_rownames(colnames(.)[c(2:5,1)]) %>% 
  arrange(row.names(.)) %>% 
  data.matrix(); dlab_HSD_pvals
# generate letters
dlab_HSD_letters <- multcompLetters(dlab_HSD_pvals,
                                    threshold=0.05,
                                    Letters=letters,
                                    reversed = FALSE); dlab_HSD_letters
# count observations
peak_counts %>%
  dplyr::group_by(donor_lab, chemotype) %>% 
  dplyr::count(counts=chemotype) # 30-40

### Recipient laboratories ====================================================
# create baseline model
rlab_chem_mod0 <- lmer(peak_counts ~ recipient_lab * chemotype +
                         (1|donor_lab/plant_ID) +
                         (1|treatment),
                       REML=F,
                       data=peak_counts,
                       control = lmerControl(
                         optimizer ='Nelder_Mead')); summary(rlab_chem_mod0)
# check assumptions visually
check_model(rlab_chem_mod0)
# update model without interaction
rlab_chem_mod1 <- update(rlab_chem_mod0, .~. - recipient_lab:chemotype)
# check assumptions visually
check_model(rlab_chem_mod1)
# LRT
anova(rlab_chem_mod0, rlab_chem_mod1)
# update model without factor donor laboratory
rlab_chem_mod2 <- update(rlab_chem_mod1, .~. -recipient_lab)
# check assumptions visually
check_model(rlab_chem_mod2)
# LRT
rlab_mod_LRT0 <- anova(rlab_chem_mod1, rlab_chem_mod2); rlab_mod_LRT0
# update model without factor chemotype
rlab_chem_mod3 <- update(rlab_chem_mod1, .~. -chemotype)
# check assumptions visually
check_model(rlab_chem_mod3)
# LRT
rlab_mod_LRT1 <- anova(rlab_chem_mod1, rlab_chem_mod3); rlab_mod_LRT1
# LMM model (type II using Satterthwaite's method of p-values)
rlab_mod <- anova(rlab_chem_mod1, type=2); rlab_mod
# post-hoc testing
rlab_HSD <- emmeans(rlab_chem_mod1,
                    list(pairwise ~ recipient_lab),
                    adjust = "tukey"); rlab_HSD$`pairwise differences of recipient_lab`
# transform to complete p-value matrix
rlab_HSD_pvals <- rlab_HSD$`pairwise differences of recipient_lab` %>% 
  data.frame() %>% 
  dplyr::select(any_of(c("X1", "p.value"))) %>% 
  separate(X1,
           into=c("V1","V2"),
           sep=" - ",
           convert=T,
           extra="merge") %>% 
  spread(V1, p.value, fill=0) %>%
  column_to_rownames(var = "V2") %>%
  mutate(L5 = c(.[3,][2:3], 0)) %>% 
  rbind(list(0, .[1,1], .[2,1], .[3,1])) %>%
  magrittr::set_rownames(colnames(.)[c(2:4,1)]) %>% 
  arrange(row.names(.)) %>% 
  data.matrix(); rlab_HSD_pvals
# generate letters
rlab_HSD_letters <- multcompLetters(rlab_HSD_pvals,
                                    threshold=0.05,
                                    Letters=letters,
                                    reversed = FALSE); rlab_HSD_letters
# count observations
peak_counts %>%
  dplyr::group_by(recipient_lab, chemotype) %>% 
  dplyr::count(counts=chemotype) # 40-50

## Venn =======================================================================
### Mono-chemotype ============================================================
# split by recipient lab
mono_lab_list <- compounds_WideFormat$mono %>% 
  split(., .$recipient_lab); mono_lab_list # 50 x 78 each
# while loop to find unique compounds per recipient lab
mono <- 1
mono_compound_list <- list()
while(mono < 5){
  mono_compound_list[[names(mono_lab_list)[mono]]] <- mono_lab_list[[names(mono_lab_list)[mono]]] %>%
    dplyr::select(contains("EIC")) %>%
    dplyr::select(-one_of("bromodecane_1_EIC")) %>% 
    janitor::remove_constant(na.rm=T) %>%
    names()
  mono <- mono+1
  }; mono_compound_list
# number of uniquely detected peaks across recipient labs
length(unique(c(mono_compound_list$L1,
                mono_compound_list$L2,
                mono_compound_list$L4,
                mono_compound_list$L))) # 64)
# generate venn object
mono_VennObj <- venn(mono_compound_list,
                     ggplot=F, borders=F); mono_VennObj

### Mixed-chemotype ===========================================================
# split by recipient lab
mixed_lab_list <- compounds_WideFormat$mixed %>% 
  split(., .$recipient_lab); mixed_lab_list # 49 x 78 each
# while loop to find unique compounds per recipient lab
mixed <- 1
mixed_compound_list <- list()
while(mixed < 5){
  mixed_compound_list[[names(mixed_lab_list)[mixed]]] <- mixed_lab_list[[names(mixed_lab_list)[mixed]]] %>%
    dplyr::select(contains("EIC")) %>%
    dplyr::select(-one_of("bromodecane_1_EIC")) %>% 
    janitor::remove_constant(na.rm=T) %>%
    names()
  mixed <- mixed+1
}; mixed_compound_list
# number of uniquely detected peaks across recipient labs
length(unique(c(mixed_compound_list$L1,
                mixed_compound_list$L2,
                mixed_compound_list$L4,
                mixed_compound_list$L))) # 64)
# generate venn object
mixed_VennObj <- venn(mixed_compound_list,
                     ggplot=F, borders=F); mixed_VennObj

# FIGURE 3 ====================================================================
## Subfigure 3A ===============================================================
# prepare asterisk notation for factor donor laboratory
asterisk_3A_dlab <- if_else(bromodecane_mod[[1]]$`Pr(Prob)`[1] < 0.001, "***",
                            if_else(bromodecane_mod[[1]]$`Pr(Prob)`[1] < 0.01, "**",
                                    if_else(bromodecane_mod[[1]]$`Pr(Prob)`[1] < 0.05, "*",
                                            "n.s."))); asterisk_3A_dlab
# prepare asterisk notation for factor chemotype
asterisk_3A_chem <- if_else(bromodecane_mod[[1]]$`Pr(Prob)`[2] < 0.001, "***",
                            if_else(bromodecane_mod[[1]]$`Pr(Prob)`[2] < 0.01, "**",
                                    if_else(bromodecane_mod[[1]]$`Pr(Prob)`[2] < 0.05, "*",
                                            "n.s."))); asterisk_3A_chem
# plot
set.seed(123); subfigure_3A <- ggplot(bromodecane_freq_descr,
                               aes(x=donor_lab,
                                   y=mean,
                                   color=chemotype,
                                   fill=chemotype)) +
  geom_point(data=bromodecane_freq,
             inherit.aes=F,
             aes(x=donor_lab,
                 y=freq,
                 color=chemotype),
             position=position_jitterdodge(),
             size=2, shape=21,
             show.legend=F) +
  geom_pointrange(aes(ymin=mean-sd,
                      ymax=mean+sd,
                      fill=chemotype,
                      color=chemotype),
                  size=1.1,
                  position=position_dodge(.8),
                  color=rep(c("black","black"),5),
                  show.legend=F) +
  geom_point(aes(y=mean,
                 x=donor_lab,
                 fill=chemotype,
                 color=chemotype),
             size=4,
             position=position_dodge(.8),
             color=rep(c("#4F94CD","#CD69C9"),5)) +
  scale_y_continuous(limits=c(0,1.05), expand=c(0,0),
                     breaks=c(0,.25,.5,.75,1),
                     labels=c(0,25,50,75,100)) +
  labs(y=expression(paste("Samples with 1-bromodecane detected [%]")),
       x="Donor laboratory",
       title=bquote(bold('LM'['Perm']*':')~
                      "Donor laboratory"~
                      "["*italic(I)[.(bromodecane_mod[[1]]$Df[1])*","*
                                      .(bromodecane_mod[[1]]$Df[3])]~
                      '='~.(bromodecane_mod[[1]]$Iter[1])*']'~.(asterisk_3A_dlab)),
       subtitle=bquote(~~~~~~~~~~~~~
                      " Chemotype"~
                      "["*italic(I)[.(bromodecane_mod[[1]]$Df[2])*","*
                                      .(bromodecane_mod[[1]]$Df[3])]~
                      '='~.(bromodecane_mod[[1]]$Iter[2])*']'~.(asterisk_3A_dlab))) +
 scale_color_manual(values=c("black","black"),
                   labels = c("Mono", "Mixed")) +
  theme_cowplot() +
  theme(legend.position = "none",#c(0.75, 0.92),
        legend.title = element_blank(),
        plot.title = element_text(size=12, face="plain"),
        plot.subtitle = element_text(size=12, face="plain"),
        panel.border=element_rect(fill="transparent",
                                  color="black",
                                  linewidth=0.8)) +
  guides(fill = guide_legend(override.aes = list(color = c("#4F94CD","#CD69C9")))) +
  annotate("text",
           x=0.5, y=.97,
           hjust=0, vjust=0,
           label=bquote(italic(n)~"= 3-4 (Receiver laboratories)")); subfigure_3A

## Subfigure 3B ===============================================================
# the following value will be removed for convenience: 
max(bromodecane$mono$bromodecane_n_normalized, na.rm=T) # 14116209 or 1.41e+07
# prepare asterisk notation for mono-chemotype
asterisk_3B_mono <- if_else(bromodecane_mono_ICC$p.value < 0.001, "***",
                            if_else(bromodecane_mono_ICC$p.value < 0.01, "**",
                                    if_else(bromodecane_mono_ICC$p.value < 0.05, "*",
                                            "n.s."))); asterisk_3B_mono
# prepare asterisk notation for mixed-chemotype
asterisk_3B_mixed <- if_else(bromodecane_mixed_ICC$p.value < 0.001, "***",
                             if_else(bromodecane_mixed_ICC$p.value < 0.01, "**",
                                     if_else(bromodecane_mixed_ICC$p.value < 0.05, "*",
                                             "n.s."))); asterisk_3B_mixed
# plot
set.seed(123); subfigure_3B <- bind_rows(bromodecane$mono,
                                                    bromodecane$mixed) %>% 
  ggplot(aes(x=recipient_lab,
             y=bromodecane_n_normalized,
             color=chemotype,
             fill=chemotype)) +
  geom_boxplot(fill=rep(c("#4F94CD","#CD69C9"),4),
               notch=F,
               outlier.shape = NA) +
  geom_point(fill="transparent",
             position=position_jitterdodge(),
             size=2, shape=21,
             show.legend=F) +
  scale_y_continuous(limits=c(0,1.6e+6), expand=c(0,0),
                     breaks=c(0, 1e+5, 5e+5, 1e+6, 1.5e+6),
                     labels=c(0, 1e+5, 5e+5, 1e+6, "1.5e+6")) +
  labs(y=expression(paste("1-bromodecane (normalised peak area)")),
       x="Recipient laboratory",
       title=bquote(bold('ICC'['Mono'])~"="~
                      .(round(bromodecane_mono_ICC$value,2))~'['*
                      italic(F)[.(bromodecane_mono_ICC$df1)*','*
                                  .(round(bromodecane_mono_ICC$df2,0))]~'='~
                      .(round(bromodecane_mono_ICC$Fvalue,2))~"]"~.(asterisk_3B_mono)),
       subtitle=bquote(bold('ICC'['Mixed'])~"="~
                         .(round(bromodecane_mixed_ICC$value,2))~'['*
                         italic(F)[.(bromodecane_mixed_ICC$df1)*','*
                                     .(round(bromodecane_mixed_ICC$df2,0))]~'='~
                         .(round(bromodecane_mixed_ICC$Fvalue,2))~"]"~.(asterisk_3B_mixed))) +
  scale_color_manual(values=c("black","black"),
                     labels=c("Mono", "Mixed")) +
  theme_cowplot() +
  theme(legend.position = c(0.75,0.92),
        legend.title = element_blank(),
        plot.title = element_text(size=12, face="plain"),
        plot.subtitle = element_text(size=12, face="plain"),
        panel.border=element_rect(fill="transparent",
                                  color="black",
                                  linewidth=0.8)) +
  guides(color = guide_legend(override.aes = list(fill = c("#4F94CD","#CD69C9")))) +
  annotate("point",
           x=3.45, y=c(1.495e+6, 1.405e+6),
           shape=21,
           size=4,
           color="black",
           fill=c("#4F94CD","#CD69C9")) +
  annotate("text",
           x=0.5, y=1.48e+6,
           hjust=0, vjust=0,
           label=bquote(italic(n)~"= 40-49")); subfigure_3B

## Subfigure 3C ===============================================================
# prepare asterisk notation for factor donor laboratory
asterisk_3C_dlab <- if_else(dlab_mod$`Pr(>F)`[1] < 0.001, "***",
                            if_else(dlab_mod$`Pr(>F)`[1] < 0.01, "**",
                                    if_else(dlab_mod$`Pr(>F)`[1] < 0.05, "*",
                                            "n.s."))); asterisk_3C_dlab
# prepare asterisk notation for factor chemotype
asterisk_3C_chem <- if_else(dlab_mod$`Pr(>F)`[2] < 0.001, "***",
                            if_else(dlab_mod$`Pr(>F)`[2] < 0.01, "**",
                                    if_else(dlab_mod$`Pr(>F)`[2] < 0.05, "*",
                                            "n.s."))); asterisk_3C_chem 
# plot
set.seed(123); subfigure_3C <- peak_counts %>% 
  ggplot(aes(x=donor_lab,
             y=peak_counts,
             color=chemotype,
             fill=chemotype)) +
  geom_boxplot(fill=rep(c("#4F94CD","#CD69C9"),5),
               notch=F,
               outlier.shape = NA) +
  geom_point(show.legend=F,
             fill="transparent",
             position=position_jitterdodge(),
             size=2, shape=21) +
  scale_y_continuous(limits=c(0,42), expand=c(0,0)) +
  labs(y="Number of peaks per sample",
       x="Donor laboratory",
       title=bquote(bold('LMM'['LRT']*":")~"Donor laboratory:"~
                      "["*italic(chi)^2~'='~
                      .(round(dlab_mod_LRT0$Chisq[2],2))*']'~.(asterisk_3C_dlab)),
       subtitle=bquote(~~~~~~~~~~~~~~~" Chemotype:"~
                                   "["*italic(chi)^2~'='~
                                   .(round(dlab_mod_LRT1$Chisq[2],2))*']'~.(asterisk_3C_chem))) +
  scale_color_manual(values=rep(c("black","black"),4),
                     labels = rep(c("Mono", "Mixed"), 4)) +
  theme_cowplot() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(size=12, face="plain"),
        plot.subtitle = element_text(size=12, face="plain"),
        panel.border=element_rect(fill="transparent",
                                  color="black",
                                  linewidth=0.8)) +
  annotate(geom="text",
           x=1:5,
           y=c((peak_counts %>% group_by(donor_lab) %>% summarise(Value = max(peak_counts))))[[2]] + 1,
           hjust=0, vjust=0,
           label=toupper(dlab_HSD_letters$Letters)) +
  annotate("text",
           x=0.5, y=39,
           hjust=0, vjust=0,
           label=bquote(italic(n)~"= 30-40")); subfigure_3C
# maximum number of peaks detected per donor laboratory
peak_counts %>% 
  group_by(donor_lab, chemotype) %>% 
  summarise(peaknax = max(peak_counts)) # 30-34

## Subfigure 3D ===============================================================
# prepare asterisk notation for factor donor laboratory
asterisk_3D_dlab <- if_else(rlab_mod$`Pr(>F)`[1] < 0.001, "***",
                            if_else(rlab_mod$`Pr(>F)`[1] < 0.01, "**",
                                    if_else(rlab_mod$`Pr(>F)`[1] < 0.05, "*",
                                            "n.s."))); asterisk_3D_dlab
# prepare asterisk notation for factor chemotype
asterisk_3D_chem <- if_else(rlab_mod$`Pr(>F)`[2] < 0.001, "***",
                            if_else(rlab_mod$`Pr(>F)`[2] < 0.01, "**",
                                    if_else(rlab_mod$`Pr(>F)`[2] < 0.05, "*",
                                            "n.s."))); asterisk_3D_chem 
# plot
set.seed(123); subfigure_3D <- peak_counts %>% 
  ggplot(aes(x=recipient_lab,
             y=peak_counts,
             color=chemotype,
             fill=chemotype)) +
  geom_boxplot(fill=rep(c("#4F94CD","#CD69C9"),4),
               notch=F,
               outlier.shape = NA) +
  geom_point(show.legend=F,
             fill="transparent",
             position=position_jitterdodge(),
             size=2, shape=21) +
  scale_y_continuous(limits=c(0,42), expand=c(0,0)) +
  labs(y="",
       x="Recipient laboratory",
       title=bquote(bold('LMM'['LRT']*":")~"Recipient laboratory:"~
                      "["*italic(chi)^2~'='~
                      .(round(rlab_mod_LRT0$Chisq[2],2))*'] ***'),
       subtitle=bquote(~~~~~~~~~~~~~~~" Chemotype:"~
                         "["*italic(chi)^2~'='~
                         .(round(rlab_mod_LRT1$Chisq[2],2))*'] n.s.')) +
  scale_color_manual(values=rep(c("black","black"),4),
                    labels = rep(c("Mono", "Mixed"), 4)) +
  theme_cowplot() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(size=12, face="plain"),
        plot.subtitle = element_text(size=12, face="plain"),
        panel.border=element_rect(fill="transparent",
                                  color="black",
                                  linewidth=0.8)) +
  annotate(geom="text",
           x=1:4,
           y=c((peak_counts %>% group_by(recipient_lab) %>% summarise(Value = max(peak_counts))))[[2]] + 1,
           hjust=0, vjust=0,
           label=toupper(rlab_HSD_letters$Letters)) +
  annotate("text",
           x=0.5, y=39,
           hjust=0, vjust=0,
           label=bquote(italic(n)~"= 40-50")); subfigure_3D
# maximum number of peaks detected per donor laboratory
peak_counts %>% 
  group_by(recipient_lab, chemotype) %>% 
  summarise(peaknax = max(peak_counts)) # 11-34

## Subfigure 3E ===============================================================
mono_VennPlot <- plot(mono_VennObj,
                      edges=T,
                      quantities=list(cex=1.3,
                                      col=c(rep("black",14), "white")),
                      
                      labels=list(labels=c("L1", "L2", "L4", "L5"), cex=1.2),
                      lty=2,
                      col=rep("gray40", 4),
                      fill=c(rep("transparent",14),
                             "gray40"),
                      lwd=2); mono_VennPlot
# ggplotify Venn diagram
subfigure_3E <- as.ggplot(mono_VennPlot) +
  labs(title=expression(paste(bold('Mono')))) +
  theme(plot.title=element_text(hjust=0.5)); subfigure_3E

## Subfigure 3E ===============================================================
mixed_VennPlot <- plot(mixed_VennObj,
                      edges=T,
                      quantities=list(cex=1.3,
                                      col=c(rep("black",14), "white")),
                      
                      labels=list(labels=c("L1", "L2", "L4", "L5"), cex=1.2),
                      lty=2,
                      col=rep("gray40", 4),
                      fill=c(rep("transparent",14),
                             "gray40"),
                      lwd=2); mixed_VennPlot
# ggplotify Venn diagram
subfigure_3F <- as.ggplot(mixed_VennPlot) +
  labs(title=expression(paste(bold('Mixed')))) +
  theme(plot.title=element_text(hjust=0.5)); subfigure_3F

## COMBINE PLOTS =============================================================
plot_grid(subfigure_3A,
          subfigure_3B,
          subfigure_3C,
          subfigure_3D,
          subfigure_3E,
          subfigure_3F,
          labels=LETTERS[1:6],
          ncol=2, nrow=3)
# and save
ggsave(last_plot(),
       file=here::here(".", # replace with desired file path
                       "Eckert_et_al_2023_Fig3.tiff"),
       width=10, height=15,
       dpi=600,
       compression="lzw")

# SESSION INFO ================================================================
sessionInfo()

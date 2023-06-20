# PACKAGES ====================================================================
library(here)
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
library(tibble)
library(performance)
library(fitdistrplus)
library(bestNormalize)
library(car)
library(agricolae)
library(mosaic)

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
pheno_data <- read_table(here(".", # replace with working directory
                              "Eckert_et_al_2023_RawData_PhenotypicData.txt"),
                          na="NA",
                          col_types=cols(unique_ID="c",
                                         donor_lab="f",
                                         chemotype="f",
                                         plant_ID="i",
                                         plant_height_cm="d",
                                         expanded_leaves_n="d",
                                         FW_shoot_g="d",
                                         FW_leaf_g="d",
                                         FW_biomass_g="d")); pheno_data # 50 x 9

# STATISTICS ==================================================================
## PLANT HEIGHT ===============================================================
# check distribution of response variable using Cullen-Frey graph
descdist(pheno_data$plant_height_cm, discrete=F)
# code interaction explicitly for later post-hoc testing
pheno_data$dlab_chemotype <- with(pheno_data, interaction(donor_lab, chemotype))
# run linear model
height_mod <- lm(plant_height_cm ~ donor_lab +
                   chemotype +
                   donor_lab:chemotype,
                 data=pheno_data,
                 contrasts=list(donor_lab=contr.sum,
                                chemotype=contr.sum)); height_mod
# check assumptions visually
check_model(height_mod)
# calculate type-III analysis-of-variance table
height_mod_sum <- Anova(height_mod, type=3); height_mod_sum
# post-hoc testing
set.seed(123); height_mod.HSD <- HSD.test(update(height_mod, .~.
                                                 -donor_lab:chemotype +
                                                   dlab_chemotype),
                                          "dlab_chemotype"); height_mod.HSD

## EXPANDED LEAVES ============================================================
# check distribution of response variable using Cullen-Frey graph
descdist(pheno_data$expanded_leaves_n, discrete=T)
# normalize response
set.seed(123); leaves_norm <- bestNormalize(pheno_data$expanded_leaves_n); leaves_norm
pheno_data$expanded_leaves_ORDERNORM <- leaves_norm$x.t
# run linear model
leaves_mod <- lm(expanded_leaves_ORDERNORM ~ donor_lab+chemotype,
                 data=pheno_data,
                 contrasts=list(donor_lab=contr.sum,
                                chemotype=contr.sum)); leaves_mod
# check assumptions visually
check_model(leaves_mod)
# calculate type-II analysis-of-variance table
leaves_mod_sum <- Anova(leaves_mod, type=2); leaves_mod_sum
# post-hoc testing
set.seed(123); leaves_mod.HSD <- HSD.test(leaves_mod, "donor_lab"); leaves_mod.HSD

## FW BIOMASS =================================================================
# check distribution of response variable using Cullen-Frey graph
descdist(pheno_data$FW_biomass_g, discrete=F)
# run linear model
FW_shoot_mod <- lm(FW_biomass_g ~ donor_lab+chemotype,
                   data=pheno_data,
                   contrasts=list(donor_lab=contr.sum,
                                  chemotype=contr.sum)); FW_shoot_mod
# check assumptions visually
check_model(FW_shoot_mod)
# calculate type-II analysis-of-variance table
FW_shoot_mod_sum <- Anova(FW_shoot_mod, type=2); FW_shoot_mod_sum
# post-hoc testing
set.seed(123); FW_shoot_mod1.HSD <- HSD.test(FW_shoot_mod, "donor_lab"); FW_shoot_mod1.HSD
set.seed(123); FW_shoot_mod2.HSD <- HSD.test(FW_shoot_mod, "chemotype"); FW_shoot_mod2.HSD

## FW LEAF ====================================================================
# check distribution of response variable using Cullen-Frey graph
descdist(pheno_data$FW_leaf_g, discrete=F)
# normalize response
set.seed(1); FW_leaf_norm <- bestNormalize(pheno_data$FW_leaf_g); FW_leaf_norm
pheno_data$FW_leaf_BC <- FW_leaf_norm$x.t
# run linear model
FW_leaf_mod <- lm(FW_leaf_BC ~ donor_lab*chemotype,
                  data=pheno_data,
                  contrasts=list(donor_lab=contr.sum,
                                 chemotype=contr.sum))
# check assumptions visually
check_model(FW_leaf_mod)
# calculate type-III analysis-of-variance table
FW_leaf_mod_sum <- Anova(FW_leaf_mod, type=3); FW_leaf_mod_sum
# post-hoc testing
set.seed(123); FW_leaf_mod.HSD <- HSD.test(update(FW_leaf_mod, .~.
                                                  -donor_lab:chemotype +
                                                    dlab_chemotype),
                                           "dlab_chemotype"); FW_leaf_mod.HSD

# FIGURE 2 ====================================================================
## Subfigure 2A ===============================================================
# PLANT HEIGHT
# prepare asterisk notation for interaction term
asterisk_2A <- if_else(height_mod_sum$`Pr(>F)`[4] < 0.001, "***",
                       if_else(height_mod_sum$`Pr(>F)`[4] < 0.01, "**",
                               if_else(height_mod_sum$`Pr(>F)`[4] < 0.05, "*",
                                       "n.s."))); asterisk_2A
# plot
set.seed(123); subfigure_2A <- pheno_data %>%
  ggplot(aes(x=donor_lab,
             y=plant_height_cm,
             color=chemotype,
             fill=chemotype)) +
  geom_boxplot(fill=rep(c("#4F94CD","#CD69C9"),5),
               notch=F,
               outlier.shape = NA) +
  geom_point(fill="transparent",
             position=position_jitterdodge(),
             size=2, shape=21) +
  scale_y_continuous(limits=c(0,28),
                     breaks=c(0,5,10,15,20,25),
                     expand=c(0,0)) +
  labs(y="Plant height [cm]",
       x=NULL,
       title="",
       subtitle=bquote("Donor laboratory"%*%"Chemotype"~"["*italic(F)[.(height_mod_sum$Df[4])*","*
                                                                        .(height_mod_sum$Df[5])]~
                         "="~.(round(height_mod_sum$`F value`[4],2))*"]"~.(asterisk_2A))) +
  scale_color_manual(values=rep(c("black","black"),4)) +
  theme_cowplot() +
  theme(legend.position = "none",
        plot.title = element_text(size=12, face="plain"),
        panel.border=element_rect(fill="transparent",
                                  color="black",
                                  linewidth=0.8)) +
  stat_summary(geom='text',
               label=height_mod.HSD$groups$groups[c(7,10,9,8,5,6,4:1)],
               fun = max,
               vjust = -1, hjust = 0.5,
               position = position_dodge(width = 0.9)) +
  annotate(geom="text",
           x=.5, y= 1,
           hjust=0, vjust=0,
           label=bquote(italic(n)~"= 5")); subfigure_2A

## Subfigure 2B ===============================================================
# EXPANDED LEAVES
# prepare asterisk notation for factor donor lab
asterisk_2B_dlab <- if_else(leaves_mod_sum$`Pr(>F)`[1] < 0.001, "***",
                       if_else(leaves_mod_sum$`Pr(>F)`[1] < 0.01, "**",
                               if_else(leaves_mod_sum$`Pr(>F)`[1] < 0.05, "*",
                                       "n.s."))); asterisk_2B_dlab
# prepare asterisk notation for factor chemotype
asterisk_2B_chem <- if_else(leaves_mod_sum$`Pr(>F)`[2] < 0.001, "***",
                            if_else(leaves_mod_sum$`Pr(>F)`[2] < 0.01, "**",
                                    if_else(leaves_mod_sum$`Pr(>F)`[2] < 0.05, "*",
                                            "n.s."))); asterisk_2B_chem
# plot
set.seed(123); subfigure_2B <- pheno_data %>%
  ggplot(aes(x=donor_lab,
             y=expanded_leaves_n,
             color=chemotype,
             fill=chemotype)) +
  geom_boxplot(fill=rep(c("#4F94CD","#CD69C9"),5),
               notch=F,
               outlier.shape = NA) +
  geom_point(show.legend=F,
             fill="transparent",
             position=position_jitterdodge(),
             size=2, shape=21) +
  scale_y_continuous(limits=c(0,28),
                     breaks=c(0,5,10,15,20,25),
                     expand=c(0,0),
                     sec.axis = sec_axis(~ predict(leaves_norm, .),
                                         name = "Ordered-quantile normalised data")) +
  labs(y="Number of leaves",
       x=NULL,
       title=bquote("Donor laboratory"~"["*italic(F)[.(leaves_mod_sum$Df[1])*","*
                                                .(leaves_mod_sum$Df[3])]~
                      "="~.(round(leaves_mod_sum$`F value`[1],2))*"]"~.(asterisk_2B_dlab)),
       subtitle=bquote("Chemotype"~"["*italic(F)[.(leaves_mod_sum$Df[2])*","*
                                                                      .(leaves_mod_sum$Df[3])]~
                         "="~.(round(leaves_mod_sum$`F value`[2],2))*"]"~.(asterisk_2B_chem))) +
  scale_color_manual(values=c("black","black"),
                    labels=c("Mono", "Mixed")) +
  theme_cowplot() +
  theme(legend.position = c(0.75, 0.9),
        legend.title = element_blank(),
        plot.title = element_text(size=12, face="plain"),
        panel.border=element_rect(fill="transparent",
                                  color="black",
                                  linewidth=0.8)) +
  guides(color = guide_legend(override.aes = list(fill = c("#4F94CD","#CD69C9")))) +
  annotate(geom="text",
           x=1:5,
           y=c((pheno_data %>%
                  group_by(donor_lab) %>%
                  summarise(Value = max(expanded_leaves_n))))[[2]] + 0.5,
           hjust=0, vjust=0,
           label=toupper(leaves_mod.HSD$groups$groups[c(1,2,4,5,3)])) +
  annotate(geom="text",
           x=.5, y= 1,
           hjust=0, vjust=0,
           label=bquote(italic(n)~"= 5")); subfigure_2B

## Subfigure 2C ===============================================================
# FW BIOMASS
# prepare asterisk notation for factor donor lab
asterisk_2C_dlab <- if_else(FW_shoot_mod_sum$`Pr(>F)`[1] < 0.001, "***",
                            if_else(FW_shoot_mod_sum$`Pr(>F)`[1] < 0.01, "**",
                                    if_else(FW_shoot_mod_sum$`Pr(>F)`[1] < 0.05, "*",
                                            "n.s."))); asterisk_2C_dlab
# prepare asterisk notation for factor chemotype
asterisk_2C_chem <- if_else(FW_shoot_mod_sum$`Pr(>F)`[2] < 0.001, "***",
                            if_else(FW_shoot_mod_sum$`Pr(>F)`[2] < 0.01, "**",
                                    if_else(FW_shoot_mod_sum$`Pr(>F)`[2] < 0.05, "*",
                                            "n.s."))); asterisk_2C_chem
# plot
set.seed(123); subfigure_2C <- pheno_data %>%
  ggplot(aes(x=donor_lab,
             y=FW_biomass_g,
             color=chemotype,
             fill=chemotype)) +
  geom_boxplot(fill=rep(c("#4F94CD","#CD69C9"),5),
               notch=F,
               outlier.shape = NA) +
  geom_point(fill="transparent",
             position=position_jitterdodge(),
             size=2, shape=21) +
 scale_y_continuous(limits=c(0,11),
                     expand=c(0,0),
                     breaks=c(0,2,4,6,8,10)) +
  labs(y="Fresh weight of aboveground biomass [g]",
       x="Donor laboratory",
       title=bquote("Donor laboratory"~"["*italic(F)[.(FW_shoot_mod_sum$Df[1])*","*
                                                       .(FW_shoot_mod_sum$Df[3])]~
                      "="~.(round(FW_shoot_mod_sum$`F value`[1],2))*"]"~.(asterisk_2C_dlab)),
       subtitle=bquote("Chemotype"~"["*italic(F)[.(FW_shoot_mod_sum$Df[2])*","*
                                                   .(FW_shoot_mod_sum$Df[3])]~
                         "="~.(round(FW_shoot_mod_sum$`F value`[2],2))*"]"~.(asterisk_2C_chem))) +
  scale_color_manual(values=rep(c("black","black"),5),
                     labels=c("Mono", "Mixed")) +
  theme_cowplot() +
  theme(legend.position = "none",
        plot.title = element_text(size=12, face="plain"),
        panel.border=element_rect(fill="transparent",
                                  color="black",
                                  linewidth=0.8)) +
  annotate(geom="text",
           x=1:5,
           y=c((pheno_data %>%
                  group_by(donor_lab) %>%
                  summarise(Value = max(FW_shoot_g))))[[2]] + 1.5,
           hjust=0, vjust=0,
           label=toupper(FW_shoot_mod1.HSD$groups$groups[c(4,5,2,1,3)])) +
  annotate(geom="text",
           x=.5, y= .4,
           hjust=0, vjust=0,
           label=bquote(italic(n)~"= 5")); subfigure_2C

## Subfigure 2D ===============================================================
# FW LEAF
# prepare asterisk notation for interaction term
asterisk_2D <- if_else(FW_leaf_mod_sum$`Pr(>F)`[4] < 0.001, "***",
                       if_else(FW_leaf_mod_sum$`Pr(>F)`[4] < 0.01, "**",
                               if_else(FW_leaf_mod_sum$`Pr(>F)`[4] < 0.05, "*",
                                       "n.s."))); asterisk_2D
# plot
set.seed(123); subfigure_2D <- pheno_data %>%
  ggplot(aes(x=donor_lab,
             y=FW_leaf_g,
             color=chemotype,
             fill=chemotype)) +
  geom_boxplot(fill=rep(c("#4F94CD","#CD69C9"),5),
               notch=F,
               outlier.shape = NA) +
  geom_point(fill="transparent",
             position=position_jitterdodge(),
             size=2, shape=21) +
  scale_y_continuous(limits=c(0,1.7),
                     expand=c(0,0),
                     sec.axis = sec_axis(~ predict(FW_leaf_norm, .), name = "Box-Cox transformed data")) +
  labs(y="Fresh weight of sampled leaf [g]",
       x="Donor laboratory",
       title="",
       subtitle=bquote("Donor laboratory"%*%"Chemotype"~"["*italic(F)[.(FW_leaf_mod_sum$Df[4])*","*
                                                                        .(FW_leaf_mod_sum$Df[5])]~
                         "="~.(round(FW_leaf_mod_sum$`F value`[4],2))*"]"~.(asterisk_2D))) +
  scale_color_manual(values=rep(c("black","black"),5),
                     labels=c("Mono", "Mixed")) +
  theme_cowplot() +
  theme(legend.position = "none",
        plot.title = element_text(size=12, face="plain"),
        panel.border=element_rect(fill="transparent",
                                  color="black",
                                  linewidth=0.8)) +
  stat_summary(geom='text',
               label=FW_leaf_mod.HSD$groups$groups[c(4,3,10,9,2,7,1,8,5:6)],
               fun = max,
               vjust = -1, hjust = 0.5,
               position = position_dodge(width = 0.9)) +
  annotate(geom="text",
           x=.5, y= .05,
           hjust=0, vjust=0,
           label=bquote(italic(n)~"= 5")); subfigure_2D

## COMBINE PLOTS ==============================================================
plot_grid(subfigure_2A,
          subfigure_2B,
          subfigure_2C,
          subfigure_2D,
          ncol=2,
          rel_heights = c(0.95,1),
          rel_widths = c(0.9,1),
          labels=LETTERS[1:4]) +
  theme(plot.margin = margin(1,1,1,1, "cm"))
# and save
ggsave(plot=last_plot(),
       file=here::here(".", # replace with desired file path
                       "Eckert_et_al_2023_Fig2.tiff"),
       width=10, height=9,
       dpi=600, bg="transparent",
       compression = "lzw")

# SESSION INFO ================================================================
sessionInfo()

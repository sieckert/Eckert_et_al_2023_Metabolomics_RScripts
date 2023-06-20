# LOAD PACKAGES ===============================================================
library(here)
library(tibble)
library(BiodiversityR) # also loads package vegan
library(ggplot2)
library(cowplot)
library(readr)
library(dplyr)
library(tidyr)
library(psych)
library(ggthemes)

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

# FUNCTION ====================================================================
`%notin%` <- function(x,y) !(x %in% y) 

# LOAD BLANKS =================================================================
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
# remove samples and blanks where bromodecane was not applied
# adjust order of factors donor laboratory and recipient laboratory
bromodecane_blanks <- bromodecane %>%
  filter(chemotype == "blank",
         treatment == "C",
         bromodecane_applied == TRUE) %>% 
  dplyr::select(-c(plant_ID, FW_leaf_g,chemotype,treatment,bromodecane_applied)) %>%
  group_by(recipient_lab) %>% 
  mutate(ID = seq(recipient_lab),
         donor_lab = factor(donor_lab,
                            levels = c("L1", "L2", "L3", "L4", "L5")),
         recipient_lab = factor(recipient_lab,
                                levels = c("L1", "L2", "L4", "L5"))); bromodecane_blanks # 36 x 5
# count observations where bromodecane was detected
bromodecane_blanks %>%
  dplyr::select(!c(ID, unique_ID, donor_lab)) %>%
  group_by(recipient_lab) %>%
  summarise(detected = sum(bromodecane_1_EIC > 0),
            not_detected = sum(bromodecane_1_EIC == 0),
            total = length(bromodecane_1_EIC),
            perc_detected = detected/total*100)
# max value
max(bromodecane_blanks$bromodecane_1_EIC) #5302261 (5.3e+06)
# prepare data for ICC
# format from long-format to wide-format
bromodecane_blanks_WideFormat <- bromodecane_blanks %>% 
  dplyr::select(-c(unique_ID,donor_lab)) %>% 
  pivot_wider(names_from = recipient_lab,
              values_from = bromodecane_1_EIC); bromodecane_blanks_WideFormat

# STATISTICS ==================================================================
## ICC ========================================================================
set.seed(123); bromodecane_blanks_ICC <- irr::icc(bromodecane_blanks_WideFormat[,-1], # remove ID column
                                 model = "twoway",
                                 type = "agreement",
                                 unit = "single"); bromodecane_blanks_ICC

# SUPPLEMENTARY FIGURE 2 ======================================================
# prepare asterisk notation for p-value
asterisk_S4 <- if_else(bromodecane_blanks_ICC$p.value < 0.001, "***",
                       if_else(bromodecane_blanks_ICC$p.value < 0.01, "**",
                               if_else(bromodecane_blanks_ICC$p.value < 0.05, "*",
                                       "n.s."))); asterisk_S4
# plot
set.seed(123); figure_S4 <- ggplot(bromodecane_blanks,
                                   aes(recipient_lab,
                                       y=bromodecane_1_EIC)) +
  geom_point(size=3, shape=21, alpha=0.5) +
  geom_boxplot(outlier.shape=NA,
               notch=F,
               fill=NA) +
  scale_y_continuous(limits=c(0,2e+6), expand=c(0,0),
                     breaks=c(0,1e+6,2e+6),
                     labels=c(0,1e+6,2e+6))+
  labs(y=expression(paste("1-bromodecane (raw peak area)")),
       x="Recipient laboratory",
       title=bquote(bold('ICC'['Blanks'])~"="~
                      .(round(bromodecane_blanks_ICC$value,2))~'['*
                      italic(F)[.(bromodecane_blanks_ICC$df1)*','*
                                  .(round(bromodecane_blanks_ICC$df2),2)]~'='~
                      .(round(bromodecane_blanks_ICC$Fvalue,2))~"]"~.(asterisk_S4))) +
  theme_cowplot() +
  theme(legend.position = c(0.75,0.9),
        legend.title = element_blank(),
        plot.title = element_text(size=12, face="plain"),
        plot.subtitle = element_text(size=12, face="plain"),
        panel.border=element_rect(fill="transparent",
                                  color="black",
                                  linewidth=0.8)) +
  annotate("text",
           x=0.5, y=5.6e+6,
           hjust=0, vjust=0,
           label="n = 12-14"); figure_S4
# and save
ggsave(last_plot(),
       file=here::here(".", # replace with desired file path
                       "Eckert_et_al_2023_SuppFig2.tiff"),
       width=5, height=5,
       dpi=600,
       compression="lzw")

# SESSION INFO ================================================================
sessionInfo()

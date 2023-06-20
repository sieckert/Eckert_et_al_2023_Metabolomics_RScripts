# PACKAGES ====================================================================
library(BiodiversityR) # also loads vegan
library(ggplot2)
library(readr)
library(dplyr)
library(ggthemes)
library(cowplot)
library(vegan)
library(tidyr)
library(colorblindr) # devtools::install_github("clauswilke/colorblindr")
library(colorspace)
library(lemon)
library(ggh4x)

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
compound_intersection <- readr::read_table(here(".", # replace with working directory
                                                "Eckert_et_al_2023_RawData_CompoundsIntersection.txt"),
                                           na="NA",
                                           col_types=cols(unique_ID = "c",
                                                          donor_lab = "f",
                                                          recipient_lab = "f",
                                                          chemotype = "f",
                                                          treatment = "f",
                                                          plant_ID = "i",
                                                          FW_leaf_g = "d",
                                                          .default = "i")) %>% 
  dplyr::select(-bromodecane_1_EIC); compound_intersection

## normalize data =============================================================
# change to matrix
compound_intersection_EIC <- dplyr::select(compound_intersection,
                                    contains(c("FW", "EIC"))); compound_intersection_EIC
# check names
names(compound_intersection_EIC)
# transform to matrix
matrix_compound_intersection_EIC <- as.matrix(compound_intersection_EIC); matrix_compound_intersection_EIC
# apply normalization (divide by bromodecane and leaf fresh weight)
matrix_compound_intersection_normalized <- apply(matrix_compound_intersection_EIC, 1,
                                                 function(x) (x)/x[1]); matrix_compound_intersection_normalized
# transform matrix to data frame
compound_intersection_normalized <- as.data.frame(t(matrix_compound_intersection_normalized)); compound_intersection_normalized
# combine with meta data
intersection_normalized <- tibble(cbind(compound_intersection[,c(1:6)],
                                        compound_intersection_normalized)); intersection_normalized

## CLEAN DATA =================================================================
# remove blanks and zero-only rows
compound_intersection_NMDS <- intersection_normalized %>%
  filter(chemotype != "blank" & treatment == "C") %>%
  dplyr::select(!c("unique_ID", "donor_lab", "recipient_lab", "chemotype", "treatment",
            "plant_ID", "FW_leaf_g")); compound_intersection_NMDS
# add row numbers as ID (to keep sample order in later percentage calculation)
compound_intersection_NMDS$id <- c(1:dim(compound_intersection_NMDS)[1])
# use percentage without standardizing the data (range between 0 and 1)
compound_intersection_NMDS_percentage <- compound_intersection_NMDS %>%
  gather(variable, value, -id) %>%
  group_by(id) %>%
  mutate(percentage = value/sum(value)*100) %>%
  dplyr::select(-value) %>%
  spread(variable, percentage) %>%
  replace(is.na(.), 1e-10); compound_intersection_NMDS_percentage
# create environmental data
compound_intersection_NMDS_env <- intersection_normalized %>%
  filter(chemotype != "blank" & treatment == "C") %>%
  dplyr::select(!c("plant_ID", "FW_leaf_g")) %>%
  filter_all(any_vars(. != 0)) %>%
  dplyr::select(c("unique_ID", "donor_lab", "recipient_lab", "chemotype", "treatment")); compound_intersection_NMDS_env

# PERCENTAGE ==================================================================
# combine
intersection_percentage_list <- cbind(compound_intersection_NMDS_env,
                           compound_intersection_NMDS_percentage) %>% 
  split(f = as.factor(.$chemotype), drop=T); intersection_percentage_list
# plyr::join_all(intersection_percentage_list)

# STATISTICS ==================================================================
## NMDS =======================================================================
set.seed(123); NMDS_mod <- metaMDS(wisconsin(sqrt(compound_intersection_NMDS_percentage[, -1])),
                                          distance="kulczynski",
                                          k=2,
                                          trymax=100,
                                          autotransform=F,
                                          plot=F); NMDS_mod
# extract NMDS data
NMDS_mod_summary <- ordiplot(NMDS_mod, choices=c(1,2))
NMDS_mod_ordsites <- sites.long(NMDS_mod_summary,
                                env.data=compound_intersection_NMDS_env)
# change levels of column Receiver_lab
NMDS_mod_ordsites$recipient_lab <- factor(NMDS_mod_ordsites$recipient_lab,
                                         levels=c("L1", "L2", "L4", "L5"))
# change levels of column Receiver_lab
NMDS_mod_ordsites$donor_lab <- factor(NMDS_mod_ordsites$donor_lab,
                                      levels=c("L1", "L2", "L3", "L4", "L5"))

## BETADISPER =================================================================
# create distance matrix from compound data
intersection_dist <- vegdist(wisconsin(sqrt(compound_intersection_NMDS_percentage[, -1])),
                         method="kulczynski"); intersection_dist
# calculate betadisper for factor donor_lab
dlab_betadisper <- betadisper(intersection_dist,
                              sqrt.dist = F,
                              group=compound_intersection_NMDS_env$donor_lab)
# calculate significance
set.seed(123); dlab_betadisper_permutest <- permutest(dlab_betadisper,
                                                      permutations=9999); dlab_betadisper_permutest
# calculate betadisper for factor chemotype
chemotype_betadisper <- betadisper(intersection_dist,
                                   sqrt.dist = F,
                                   group=compound_intersection_NMDS_env$chemotype)
# calculate significance
set.seed(123); chemotype_betadisper_permutest <- permutest(chemotype_betadisper,
                                                           permutations=9999); chemotype_betadisper_permutest
# calculate betadisper for factor receiver_lab
rlab_betadisper <- betadisper(intersection_dist,
                              sqrt.dist = F,
                              group=compound_intersection_NMDS_env$recipient_lab)
# calculate significance
set.seed(123); rlab_betadisper_permutest <- permutest(rlab_betadisper,
                                                      permutations=9999); rlab_betadisper_permutest

## npMANOVA ===================================================================
# default npMANOVA by terms for factor chemotype
set.seed(123); chemotype_npmanova <- adonis2(intersection_dist ~
                                               chemotype,
                                             by="terms",
                                             data = compound_intersection_NMDS_env,
                                             strata=compound_intersection_NMDS_env$donor_lab:compound_intersection_NMDS_env$recipient_lab,
                                         permutations=9999); chemotype_npmanova

# FIGURE 4 ====================================================================
## Subfigure 4A ===============================================================
### MONO ======================================================================
# prepare data
mono_averages <- intersection_percentage_list$mono %>%
  dplyr::select(-c(unique_ID, chemotype, treatment)) %>% 
  group_by(id) %>% 
  mutate(perc_sum = sum(borneol_EIC,
                        camphor_EIC,
                        cymene_p_EIC,
                        eucalyptol_EIC,
                        limonene_EIC,
                        pinene_a_EIC,
                        thujone_a_EIC,
                        thujone_b_EIC),
         .after=id) %>% 
  filter(perc_sum == 100) %>%
  ungroup() %>% 
  dplyr::select(-one_of("perc_sum", "id")) %>% 
  pivot_longer(!c("donor_lab",
                  "recipient_lab"),
               names_to="compounds",
               values_to="peak_perc") %>% 
  group_by(donor_lab, recipient_lab, compounds) %>% 
  summarise(perc_mean = mean(peak_perc)) %>% 
  # add dummy variable for x-axis in stacked barplot
  mutate(treatment = "C", .after=recipient_lab); mono_averages
# plot
subfigure_4A_1 <- ggplot(mono_averages,
                                    aes(fill=compounds,
                                        y=perc_mean,
                                        x=treatment)) +
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous(expand=c(0,0),
                     breaks=c(0,.25,.5,.75,1),
                     labels=c(0,25,50,75,100)) +
  facet_rep_grid(recipient_lab ~ donor_lab,
                 scales="free",
                 switch="x",
                 repeat.tick.labels = FALSE) +
  labs(title="Mono",
       x=NULL,
       y="Normalized relative peak area [%]",
       fill="Compound") +
  scale_fill_brewer(labels = c("borneol",
                               "camphor",
                               bquote(italic("p")*"-cymene"),
                               "eucalyptol",
                               "limonene",
                               bquote(alpha*"-pinene"),
                               bquote(alpha*"-thujone"),
                               bquote(beta*"-thujone")),
                    palette = "Purples") +
  theme_cowplot() +
  theme(legend.position = "top",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing.x=unit(-1.6, "lines"),
        panel.spacing.y=unit(0.8, "lines"),
        strip.background = element_rect(colour="black", fill="transparent")); subfigure_4A_1

### MIXED =====================================================================
# calculate mean of data
mixed_averages <- intersection_percentage_list$mixed %>%
  dplyr::select(-c(unique_ID, chemotype, treatment)) %>% 
  group_by(id) %>% 
  mutate(perc_sum = sum(borneol_EIC,
                        camphor_EIC,
                        cymene_p_EIC,
                        eucalyptol_EIC,
                        limonene_EIC,
                        pinene_a_EIC,
                        thujone_a_EIC,
                        thujone_b_EIC),
         .after=id) %>% 
  filter(perc_sum == 100) %>%
  ungroup() %>% 
  dplyr::select(-one_of("perc_sum", "id")) %>% 
  pivot_longer(!c("donor_lab",
                  "recipient_lab"),
               names_to="compounds",
               values_to="peak_perc") %>% 
  group_by(donor_lab, recipient_lab, compounds) %>% 
  summarise(perc_mean = mean(peak_perc)) %>% 
  # add dummy variable for x-axis in stacked barplot
  mutate(treatment = "C", .after=recipient_lab); mixed_averages
# and plot
subfigure_4A_2 <- ggplot(mixed_averages,
                         aes(fill=compounds,
                             y=perc_mean,
                             x=treatment)) +
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous(expand=c(0,0),
                     breaks=c(0,.25,.5,.75,1),
                     labels=c(0,25,50,75,100),
                     sec.axis = dup_axis(guide=guide_axis(title="Recipient laboratory"))) +
  facet_rep_grid(recipient_lab ~ donor_lab,
                 scales="free",
                 switch="x",
                 repeat.tick.labels = FALSE) +
  labs(title="Mixed",
       #  subtitle="n = 5",
       x=NULL,
       y=NULL,
       fill="Compound") +
  scale_fill_brewer(labels = c("borneol",
                               "camphor",
                               bquote(italic("p")*"-cymene"),
                               "eucalyptol",
                               "limonene",
                               bquote(alpha*"-pinene"),
                               bquote(alpha*"-thujone"),
                               bquote(beta*"-thujone")),
                    palette = "Purples") +
  theme_cowplot() +
  theme(legend.position="top",
        legend.title=element_text(color="transparent"),
        legend.text=element_text(color="transparent"),
        axis.line.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        panel.spacing.x=unit(-1.8, "lines"),
        panel.spacing.y=unit(0.8, "lines"),
        strip.background = element_rect(colour="black", fill="transparent")) +
  guides(fill=guide_legend(override.aes=list(fill=NA))); subfigure_4A_2

## Subfigure 4B ===============================================================
# find convex hulls
dlab_hulls <- NMDS_mod_ordsites %>% 
  group_by(donor_lab) %>% 
  slice(chull(axis1, axis2))
# count oberversions
NMDS_mod_ordsites %>%
  count(donor_lab) # 30-40
# NMDS plot
subfigure_4B <- ggplot(NMDS_mod_ordsites,
                            aes(x=axis1, y=axis2,
                                colour=donor_lab,
                                fill=donor_lab)) + 
  geom_point(data=NMDS_mod_ordsites,
             size=4, alpha = 0.7,
             aes(shape=chemotype)) +
  geom_polygon(data = dlab_hulls,
               aes(),
               alpha = 0,
               show.legend = FALSE) +
  labs(color="Donor laboratory",
       x="NMDS1",
       y="NMDS2",
       title=bquote(bold("Betadisper:")~"Donor laboratory"~
                      "["*italic('F')[.(dlab_betadisper_permutest$tab$Df[1])*","*
                                        .(dlab_betadisper_permutest$tab$Df[2])]~
                      "="~.(round(dlab_betadisper_permutest$statistic,2))*"] **"),
       subtitle=bquote(atop(~~~~~~~~~~~~~~~~~~~~" Chemotype"~
                              "["*italic('F')[.(chemotype_betadisper_permutest$tab$Df[1])*","*
                                                .(chemotype_betadisper_permutest$tab$Df[2])]~
                              "="~.(round(chemotype_betadisper_permutest$statistic,2))*"] n.s.",
                            bold("npMANOVA:")~"Chemotype"~
                              "["*italic('F')[.(chemotype_npmanova$Df[1])*","*
                                                .(chemotype_npmanova$Df[3])]~
                              "="~.(round(chemotype_npmanova$F[1],2))*"] ***"))) +
  scale_x_continuous(breaks=c(-1.5,-1,-0.5,0,0.5,1),
                     labels=c(-1.5,-1,-0.5,0,0.5,1)) +
  scale_y_continuous(breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5),
                     labels=c(-1.5,-1,-0.5,0,0.5,1,1.5)) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  scale_shape_manual(values=c(1,16)) +
  theme_cowplot() +
  theme(panel.border=element_rect(colour="black", fill=NA, linewidth=0.8),
        legend.position="top",
        legend.title = element_text(size=12, face="plain"),
        legend.text = element_text(size=11, face="plain"),
        plot.title = element_text(size=12, face="plain"),
        axis.title = element_text(size=12, face="plain")) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         fill="none",
         shape="none") +
  scale_color_manual(values=c("#000000",
                              "#E69F00",
                              "#CC79A7",
                              "#56B4E9",
                              "#009E73"),
                     labels=c("L1", "L2", "L3", "L4", "L5")) +
  scale_fill_manual(values=rep("transparent",5)) +
  annotate("point",
           x=c(-1.76,-1.76),
           y=c(-1.6,-1.85),
           size=4,
           shape=c(1,16),
           color="gray50") +
  annotate("text",
           x=c(-1.66,-1.66),
           y=c(-1.6,-1.85),
           hjust=0, vjust=0.5,
           label=c("Mono","Mixed")) +
  annotate("text",
           x=.7, y=1.5,
           hjust=0, vjust=0,
           size=3.5,
           label=bquote("stress ="~.(round(NMDS_mod$stress,2)))); subfigure_4B

## Subfigure 4C ===============================================================
# find convex hulls
rlab_hulls <- NMDS_mod_ordsites %>% 
  group_by(recipient_lab) %>% 
  slice(chull(axis1, axis2))
# count oberversions
NMDS_mod_ordsites %>%
  count(recipient_lab) # 40-50
# NMDS plot
subfigure_4C <- ggplot(NMDS_mod_ordsites,
                               aes(x=axis1, y=axis2,
                                   colour=recipient_lab,
                                   fill=recipient_lab)) + 
  geom_point(data=NMDS_mod_ordsites,
             size=4, alpha = 0.7,
             aes(shape=chemotype)) +
  geom_polygon(data = rlab_hulls,
               aes(),
               alpha = 0,
               show.legend = FALSE) +
  labs(color="Recipient laboratory",
       x="NMDS1",
       y="NMDS2",
       title=bquote(bold("Betadisper:")~"Recipient laboratory"~
                      "["*italic('F')[.(rlab_betadisper_permutest$tab$Df[1])*","*
                                        .(rlab_betadisper_permutest$tab$Df[2])]~
                      "="~.(round(rlab_betadisper_permutest$statistic,2))*"] ***"))+#,
  scale_x_continuous(breaks=c(-1.5,-1,-0.5,0,0.5,1),
                     labels=c(-1.5,-1,-0.5,0,0.5,1)) +
  scale_y_continuous(breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5),
                     labels=c(-1.5,-1,-0.5,0,0.5,1,1.5)) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  scale_shape_manual(values=c(1,16)) +
  theme_cowplot() +
  theme(panel.border=element_rect(colour="black", fill=NA, linewidth=0.8),
        legend.position="top",
        legend.title = element_text(size=12, face="plain"),
        legend.text = element_text(size=11, face="plain"),
        plot.title = element_text(size=12, face="plain"),
        axis.title = element_text(size=12, face="plain")) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         fill="none",
         shape="none") +
  scale_color_manual(values=c("#000000",
                              "#E69F00",
                              "#56B4E9",
                              "#009E73"),
                     labels=c("L1", "L2","L4", "L5")) +
  scale_fill_manual(values=rep("transparent",4)) +
  annotate("point",
           x=c(-1.76,-1.76),
           y=c(-1.6,-1.85),
           size=4,
           shape=c(1,16),
           color="gray50") +
  annotate("text",
           x=c(-1.66,-1.66),
           y=c(-1.6,-1.85),
           hjust=0, vjust=0.5,
           label=c("Mono","Mixed")) +
  annotate("text",
           x=.7, y=1.5,
           hjust=0, vjust=0,
           size=3.5,
           label=bquote("stress ="~.(round(NMDS_mod$stress,2)))); subfigure_4C

## COMBINE PLOTS ==============================================================
# subfigures 4A_1 and 4A_2 combined
subfigure_4A <- plot_grid(subfigure_4A_1,
                          subfigure_4A_2,
                          ncol=2) +
  theme(plot.margin = unit(c(0,0,0.8,0), "cm")) +
  annotate(geom = "text",
           x = 0.4, y = -0.015,
           hjust=0, vjust=0,
           label = "Donor laboratory",
           size=4.5)
# subfigures 4B and 4C combined
subfigures_4B_4C <- plot_grid(subfigure_4B,
                              subfigure_4C,
                              labels=LETTERS[2:3],
                              ncol=1, nrow=2,
                              rel_heights = c(1,.89))
# all figures combined
plot_grid(subfigure_4A,
          subfigures_4B_4C,
          labels=c("A",NA),
          rel_widths = c(1,.65),
          ncol=2, nrow=1)
# save
ggsave(last_plot(),
       file=here::here(".",
                       paste0(format(Sys.time(), "%Y%m%d_"),
                              "Eckert_et_al_2023_Fig4.tiff")),
       width=12, height=11, dpi=300,
       compression="lzw")

# SESSION INFO ================================================================
sessionInfo()

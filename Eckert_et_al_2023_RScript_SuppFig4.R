# PACKAGES ====================================================================
library(here)
library(ggplot2)
library(cowplot)
library(dplyr)
library(readr)
library(lme4)
library(ggthemes)
library(scales)
library(forcats)

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
standards <- readr::read_table(here(".", # replace with working directory
                                    "Eckert_et_al_2023_RawData_Standards.txt"),
                                 na="NA",
                                 col_types=cols(unique_ID = "c",
                                                compound = "f",
                                                concentration_ng="d",
                                                recipient_lab="f",
                                                slope="l",
                                                saturated="l",
                                                base_peak="i",
                                                EIC_RT="d",
                                                EIC_RI="i",
                                                EIC_area="i")) %>% 
  filter(slope == 1); standards
# add scaled EIC peak area (rescaled to [0,1])
standards$EIC_area_rescaled <- rescale(standards$EIC_area)

# STATISTICS ==================================================================
# calculate response factor for each compound per recipient lab
rlab_names <- unique(standards$recipient_lab); rlab_names
std_names <- unique(standards$compound); std_names
# extract slope and r-squared from simple linear regression in a loop
compound_fit <- data.frame()
rf_fits <- for (j in std_names) {
  for (i in rlab_names) {
    mod <- lm(EIC_area_rescaled ~ 0 + concentration_ng,
              data=standards %>% 
                filter(recipient_lab == i,
                       standards$compound == j))
    # model estimates per compound and for each receiver lab separately
    df <- tibble(recipient_lab = as.factor(i),
                 compound = as.factor(j),
                 r2 = round(summary(mod)$r.squared, 4),
                 slope = mod$coefficients[[1]])
    # append to dataframe
    compound_fit <- rbind(df, compound_fit)
  }
  # change to tibble
  compound_fit <- tibble(compound_fit)
}; compound_fit
# trim to the same unit (1e-04)
compound_fit$slope <- round(compound_fit$slope, 5); compound_fit
# sort facotr levels
compound_fit$recipient_lab <- factor(compound_fit$recipient_lab,
                                       levels=c("L1","L2","L4","L5"))
# split to list
compound_fit_list <- split(compound_fit,
                          f=compound_fit$compound); compound_fit_list
# sort all tables in list by recipient_lab column
sorted_compound_fit_list <- list()
for(i in names(compound_fit_list)){
  sorted_compound_fit_list[[i]] <- compound_fit_list[[i]] %>% 
    arrange(recipient_lab)
  }; sorted_compound_fit_list
compound_fit_list <- sorted_compound_fit_list; compound_fit_list

# SUPPLEMENTARY FIGURE 4 ======================================================
## subfigure S4A ==============================================================
subfigure_S4A <- standards %>% 
  filter(compound == "camphor") %>% 
  ggplot(aes(x=concentration_ng,
             y=EIC_area_rescaled,
             color=fct_inorder(recipient_lab),
             fill=fct_inorder(recipient_lab))) +
  geom_hline(yintercept=0, linetype=2, color="gray40") +
  geom_point(size=3) +
  geom_smooth(method = "lm", se=T, na.rm = T, alpha=0.3, formula=y~x-1) +
  scale_y_continuous(limits=c(-.05,.35),
                     expand=c(0,0)) +
  scale_x_continuous(limits=c(0,110),
                     breaks=c(0,20,40,60,80,100),
                     labels=c(0,20,40,60,80,100),
                     expand=c(0,0)) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  labs(title=expression(paste(bold("Compound:")~plain("camphor"))),
       color="Receiver laboratory",
       y="Normalised peak area (EIC)",
       x="") +
  theme_cowplot() +
  theme(legend.position="top",
        panel.border = element_rect(fill=NULL, color="black", linewidth=0.8),
        panel.grid.major = element_line(color = "grey90")) +
  guides(color=guide_legend(override.aes=list(fill=NA, size=6, shape=NA)),
         fill="none") +
  annotate("text",
           size=3,
           x=rep(2,4),
           y=c(.32,.29,.26,.23),
           hjust=0, vjust=0, parse=T,
           label=c(bquote(bold(.(as.character(compound_fit_list$camphor$recipient_lab [4]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$camphor$r2[4]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$camphor$slope[4],
                                                      format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$camphor$recipient_lab [3]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$camphor$r2[3]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$camphor$slope[3],
                                                      format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$camphor$recipient_lab [2]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$camphor$r2[2]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$camphor$slope[2],
                                                      format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$camphor$recipient_lab [1]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$camphor$r2[1]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$camphor$slope[1],
                                                      format="f", digits=5))))); subfigure_S4A

## subfigure S4C ==============================================================
subfigure_S4B <- standards %>% 
  filter(compound == "caryophyllene_b") %>% 
  ggplot(aes(x=concentration_ng,
             y=EIC_area_rescaled,
             color=fct_inorder(recipient_lab),
             fill=fct_inorder(recipient_lab))) +
  geom_hline(yintercept=0, linetype=2, color="gray40") +
  geom_point(size=3) +
  geom_smooth(method = "lm", se=T, na.rm = T, alpha=0.3, formula=y~x-1) +
  scale_y_continuous(limits=c(-.05,.4),
                     expand=c(0,0)) +
  scale_x_continuous(limits=c(0,110),
                     breaks=c(0,20,40,60,80,100),
                     labels=c(0,20,40,60,80,100),
                     expand=c(0,0)) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  labs(title=expression(paste(bold("Compound:")~beta*plain("-caryophyllene"))),
       color="Receiver laboratory",
       y="",
       x="") +
  guides(color=guide_legend(override.aes=list(fill=NA, color=NA, size=0, shape=NA)),
         fill="none") +
  theme_cowplot() +
  theme(legend.position="top",
        legend.text=element_text(color="transparent"),
        legend.title=element_text(color="transparent"),
        panel.border = element_rect(fill=NULL, color="black", linewidth =0.8),
        panel.grid.major = element_line(color = "grey90")) +
  annotate("text",
           size=3,
           x=rep(2,4),
           y=c(.36,.325,.29,.255),
           #y=c(3.7e+7, 4.2e+7, 4.7e+7, 5.2e+7),
           hjust=0, vjust=0, parse=T,
           label=c(bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[4]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[4]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[4],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[3]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[3]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[3],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[2]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[2]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[2],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[1]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[1]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[1],
                                                            format="f", digits=5))))); subfigure_S4B

## subfigure S4C ==============================================================
subfigure_S4C <- standards %>% 
  filter(compound == "cymene_p") %>% 
  ggplot(aes(x=concentration_ng,
             y=EIC_area_rescaled,
             color=fct_inorder(recipient_lab),
             fill=fct_inorder(recipient_lab))) +
  geom_hline(yintercept=0, linetype=2, color="gray40") +
  geom_point(size=3) +
  geom_smooth(method = "lm", se=T, na.rm = T, alpha=0.3, formula=y~x-1) +
  scale_y_continuous(limits=c(-.05,1.32),
                     expand=c(0,0)) +
  scale_x_continuous(limits=c(0,110),
                     breaks=c(0,20,40,60,80,100),
                     labels=c(0,20,40,60,80,100),
                     expand=c(0,0)) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  labs(title=expression(paste(bold("Compound:")~italic("p")*plain("-cymene"))),
       color="Receiver laboratory",
       y="Normalised peak area (EIC)",
       x="") +
  guides(color=guide_legend(override.aes=list(fill=NA, color=NA, size=0, shape=NA)),
         fill="none") +
  theme_cowplot() +
  theme(legend.position="none",
        legend.text=element_text(color="transparent"),
        legend.title=element_text(color="transparent"),
        panel.border = element_rect(fill=NULL, color="black", linewidth=0.8),
        panel.grid.major = element_line(color = "grey90")) +
  annotate("text",
           size=3,
           x=rep(2,4),
           y=c(1.2,1.1,1.0,.9),
           #y=c(1.95e+8, 2.2e+8, 2.45e+8, 2.7e+8),
           hjust=0, vjust=0, parse=T,
           label=c(bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[4]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[4]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[4],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[3]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[3]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[3],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[2]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[2]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[2],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[1]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[1]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[1],
                                                            format="f", digits=5))))); subfigure_S4C

## subfigure S4D ==============================================================
subfigure_S4D <- standards %>% 
  filter(compound == "hexenyl_acetate_3_Z") %>% 
  ggplot(aes(x=concentration_ng,
             y=EIC_area_rescaled,
             color=fct_inorder(recipient_lab),
             fill=fct_inorder(recipient_lab))) +
  geom_hline(yintercept=0, linetype=2, color="gray40") +
  geom_point(size=3) +
  geom_smooth(method = "lm", se=T, na.rm = T, alpha=0.3, formula=y~x-1) +
  scale_y_continuous(limits=c(-.05,.42),
                     expand=c(0,0)) +
  scale_x_continuous(limits=c(0,110),
                     breaks=c(0,20,40,60,80,100),
                     labels=c(0,20,40,60,80,100),
                     expand=c(0,0)) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  labs(title=expression(paste(bold("Compound:")~italic("(Z)")*plain("-hex-3-enyl-acetate"))),
       color="Receiver laboratory",
       y="",
       x=bquote("Concentration [ng µl"^-1*"]")) + #\u00B7
  guides(color=guide_legend(override.aes=list(fill=NA, color=NA, size=0, shape=NA)),
         fill="none") +
  theme_cowplot() +
  theme(legend.position="none",
        legend.text=element_text(color="transparent"),
        legend.title=element_text(color="transparent"),
        panel.border = element_rect(fill=NULL, color="black", linewidth=0.8),
        panel.grid.major = element_line(color = "grey90")) +
  annotate("text",
           size=3,
           x=rep(2,4),
           y=c(.38,.345,.31,.275),
           hjust=0, vjust=0, parse=T,
           label=c(bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[4]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[4]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[4],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[3]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[3]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[3],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[2]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[2]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[2],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[1]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[1]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[1],
                                                            format="f", digits=5))))); subfigure_S4D

## subfigure S4E ===================================================================
subfigure_S4E <- standards %>% 
  filter(compound == "pinene_a") %>% 
  ggplot(aes(x=concentration_ng,
             y=EIC_area_rescaled,
             color=fct_inorder(recipient_lab),
             fill=fct_inorder(recipient_lab))) +
  geom_hline(yintercept=0, linetype=2, color="gray40") +
  geom_point(size=3) +
  geom_smooth(method = "lm", se=T, na.rm = T, alpha=0.3, formula=y~x-1) +
  scale_y_continuous(limits=c(-.05,.81),
                     expand=c(0,0)) +
  scale_x_continuous(limits=c(0,110),
                     breaks=c(0,20,40,60,80,100),
                     labels=c(0,20,40,60,80,100),
                     expand=c(0,0)) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  labs(title=expression(paste(bold("Compound:")~alpha*plain("-pinene"))),
       color="Receiver laboratory",
       y="Normalised peak area (EIC)",
       x=bquote("Concentration [ng µl"^-1*"]")) + # \u00B7
  guides(color=guide_legend(override.aes=list(fill=NA, color=NA, size=0, shape=NA)),
         fill="none") +
  theme_cowplot() +
  theme(legend.position="none",
        legend.text=element_text(color="transparent"),
        legend.title=element_text(color="transparent"),
        panel.border = element_rect(fill=NULL, color="black", size=0.8),
        panel.grid.major = element_line(color = "grey90")) +
  annotate("text",
           size=3,
           x=rep(2,4),
           y=c(.73,.66,.59,.52),
           hjust=0, vjust=0, parse=T,
           label=c(bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[4]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[4]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[4],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[3]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[3]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[3],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[2]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[2]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[2],
                                                            format="f", digits=5))),
                   bquote(bold(.(as.character(compound_fit_list$caryophyllene$recipient_lab[1]))*":")~
                            italic("R")^2~"="~.(sprintf('%.2f', compound_fit_list$caryophyllene$r2[1]))*":"~
                            italic("b")['Rf']~"="~.(formatC(compound_fit_list$caryophyllene$slope[1],
                                                            format="f", digits=5))))); subfigure_S4E

## COMBINE PLOTS ==============================================================
plot_grid(subfigure_S4A,
          subfigure_S4B,
          subfigure_S4C,
          subfigure_S4D,
          subfigure_S4E,
          ncol=2,
          nrow=3,
          labels=LETTERS[1:5],
          rel_heights = c(1,0.9,0.9))
# and save
ggsave(last_plot(),
       file=here::here(".", # replace with desired file path
                       "Eckert_et_al_2023_SuppFig4.tiff"),
       width=10, height=11,
       dpi=600,
       compression="lzw")

# SESSION INFO ================================================================
sessionInfo()

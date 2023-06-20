# PACKAGES ====================================================================
library(here)
library(readr)
library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(janitor)
library(ggplot2)
library(lme4)
library(lmerTest)
library(performance)

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
# Function to calculate log2FC values from nested data
log2FC <- function(count_data,
                   meta_data,
                   numerator,
                   denominator,
                   samples = "Sample",
                   response = "Metabolite",
                   treatment = "Treatment",
                   log_transform = F){
  meta <- meta_data
  # tidy data
  counts <- count_data %>% 
    column_to_rownames(var=response) %>% 
    dplyr::select(-any_of(response)) %>%
    mutate(across(everything(), as.numeric)) %>% 
    sjmisc::rotate_df() %>%
    rownames_to_column(var=samples) %>% 
    left_join(meta, by=samples) %>% 
    mutate(ID = as.factor(ID)) %>% 
    column_to_rownames(samples) %>% 
    relocate(c("Treatment","ID"), .before=everything())
  # check (and apply) log transformation
  ifelse(log_transform == T,
         counts_mod <- counts %>%
           mutate_if(., is.numeric, log),
         counts_mod <- counts)
  # calculate log2-fold changes
  counts_mod_groups <- counts_mod %>% 
    pivot_longer(cols=where(is.numeric),
                 names_to="Compound",
                 values_to="Value") %>% 
    group_by(Compound) %>% 
    pivot_wider(names_from = treatment,
                values_from = "Value") %>%
    filter(complete.cases(.data[[numerator]], .data[[denominator]])) %>% 
    mutate(FC = .data[[numerator]]/.data[[denominator]]) %>% 
    dplyr::select(-any_of(c(numerator, denominator))) %>%
    summarize(log2FC_median = median(log2(FC)))
  # summarize(log2FC_median = log2(median(FC)))
  print("Yap, it worked!")
  return(counts_mod_groups)
}

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

## clean data =================================================================
# remove blanks and fresh weight column
compound_intersection_log2FC <- intersection_normalized %>%
  filter(chemotype != "blank") %>%
  dplyr::select(-FW_leaf_g); compound_intersection_log2FC

# LOG2FC ======================================================================
## add thresholds =============================================================
# min and max without zeros (this function takes unusually long, so maybe safe the output or use caching)
min_max <- compound_intersection_log2FC %>% 
  dplyr::select(dplyr::contains(c("recipient", "EIC"))) %>% 
  pivot_longer(cols=contains("EIC"),
               names_to = "compound",
               values_to = "value") %>%
  naniar::replace_with_na_all(data = .,
                              condition = ~.x == 0) %>%
  drop_na(value) %>% 
  group_by(recipient_lab) %>% 
  summarise(max_value = max(value),
            min_value = min((value))) %>%
  ungroup(); min_max
# add thresholds as reduced values by two orders of magnitude
min_max$threshold <- c(10,10,1,10); min_max
# add reduced threshold values to all values (not only zeros)
foldchange_thresholds <- compound_intersection_log2FC %>% 
  pivot_longer(cols=contains("EIC"),
               names_to = "compound",
               values_to = "value") %>% 
  group_by(recipient_lab) %>% 
  mutate(value_new = ifelse(recipient_lab == min_max$recipient_lab[1],
                            value+min_max$threshold[1],
                            ifelse(recipient_lab == min_max$recipient_lab[2],
                                   value+min_max$threshold[2],
                                   ifelse(recipient_lab == min_max$recipient_lab[3],
                                          value+min_max$threshold[3],
                                          value+min_max$threshold[3])))) %>%
  ungroup() %>% 
  dplyr::select(-"value") %>% 
  pivot_wider(names_from = "compound",
              values_from = "value_new"); foldchange_thresholds
# split by chemotype
log2FC_mono_data <- foldchange_thresholds[foldchange_thresholds$chemotype=="mono",]
log2FC_mixed_data <- foldchange_thresholds[foldchange_thresholds$chemotype=="mixed",]
# check if zeros present
min(log2FC_mono_data[,-c(1:6)]) # 1
min(log2FC_mixed_data[,-c(1:6)]) # 1

## MONO =======================================================================
# (1) Calculate log2FC by recipient_lab within donor_lab
# (2) Do this for every receiver_lab x donor_lab combination
# (3) Average calculated log2-fold changes using the median --> sample size n = 3-4
log2FC_mono_list <- list()
for (i in unique(log2FC_mono_data$donor_lab)){
  for (j in unique(log2FC_mono_data$recipient_lab)){
    if(any(log2FC_mono_data$donor_lab == i &
           log2FC_mono_data$recipient_lab == j)){
      # create counts data
      counts <- log2FC_mono_data %>%
        filter(donor_lab == i,
               recipient_lab == j) %>%
        dplyr::select(unique_ID, borneol_EIC, camphor_EIC, cymene_p_EIC, eucalyptol_EIC,
               limonene_EIC, pinene_a_EIC, thujone_a_EIC, thujone_b_EIC) %>%
        sjmisc::rotate_df(rn="Metabolite", cn=T)
      # create meta data
      meta <- log2FC_mono_data %>%
        filter(donor_lab == i,
               recipient_lab == j) %>%
        dplyr::select(unique_ID, treatment) %>%
        group_by(treatment) %>%
        mutate(ID = seq(treatment)) %>%
        #mutate(ID = rep(1:5,2)) %>%
        rename(Sample = unique_ID,
               Treatment = treatment)
      # calculate log2FC
      tmp <- log2FC(count_data = counts,
                    meta_data = meta,
                    numerator="JA",
                    denominator="C",
                    samples="Sample",
                    response="Metabolite",
                    treatment="Treatment",
                    log_transform=F)
      log2FC_mono_list[[paste0(i,"_",j)]] <- tmp$log2FC_median
    }
  }
}
# unlist the list
mono_log2FCs <- log2FC_mono_list %>% 
  unlist() %>% 
  data.frame() %>% 
  rownames_to_column() %>% 
  rename(dlab_rlab = "rowname",
         log2FC = '.') %>%
  tibble() %>% 
  mutate(compounds = rep(counts$Metabolite, 19),
         .before = log2FC) %>% 
  separate(col = dlab_rlab, into = c("donor_lab", "receiver_lab"), sep = "_"); mono_log2FCs
# average log2FC values within factor donor_lab
mono_log2FCs_median <- mono_log2FCs %>% 
  group_by(donor_lab, compounds) %>% 
  summarize(log2FC_median = median(log2FC)); mono_log2FCs_median # 40 x 3
# check
min(mono_log2FCs_median$log2FC_median); max(mono_log2FCs_median$log2FC_median)

## MIXED ======================================================================
log2FC_mixed_list <- list()
for (i in unique(log2FC_mixed_data$donor_lab)){
  for (j in unique(log2FC_mixed_data$recipient_lab)){
    if(any(log2FC_mixed_data$donor_lab == i &
           log2FC_mixed_data$recipient_lab == j)){
      # create counts data
      counts <- log2FC_mixed_data %>%
        filter(donor_lab == i,
               recipient_lab == j) %>%
        dplyr::select(unique_ID, borneol_EIC, camphor_EIC, cymene_p_EIC, eucalyptol_EIC,
               limonene_EIC, pinene_a_EIC, thujone_a_EIC, thujone_b_EIC) %>%
        sjmisc::rotate_df(rn="Metabolite", cn=T)
      # create meta data
      meta <- log2FC_mixed_data %>%
        filter(donor_lab == i,
               recipient_lab == j) %>%
        dplyr::select(unique_ID, treatment) %>%
        group_by(treatment) %>%
        mutate(ID = seq(treatment)) %>%
        #mutate(ID = rep(1:5,2)) %>%
        rename(Sample = unique_ID,
               Treatment = treatment)
      # calculate log2FC
      tmp <- log2FC(count_data = counts,
                    meta_data = meta,
                    numerator="JA",
                    denominator="C",
                    samples="Sample",
                    response="Metabolite",
                    treatment="Treatment",
                    log_transform=F)
      log2FC_mixed_list[[paste0(i,"_",j)]] <- tmp$log2FC_median
    }
  }
}
# unlist the list
mixed_log2FCs <- log2FC_mixed_list %>% 
  unlist() %>% 
  data.frame() %>% 
  rownames_to_column() %>% 
  rename(dlab_rlab = "rowname",
         log2FC = '.') %>%
  tibble() %>% 
  mutate(compounds = rep(counts$Metabolite, 19),
         .before = log2FC) %>% 
  separate(col = dlab_rlab, into = c("donor_lab", "receiver_lab"), sep = "_"); mixed_log2FCs
# average log2FC values within factor donor_lab
mixed_log2FCs_median <- mixed_log2FCs %>% 
  group_by(donor_lab, compounds) %>% 
  summarize(log2FC_median = median(log2FC)); mixed_log2FCs_median # 40 x 3
# check
min(mixed_log2FCs_median$log2FC_median); max(mixed_log2FCs_median$log2FC_median)

# STATISTICS ==================================================================
# calculate individual log2FC values per plant to apply mixed-effects model
## MONO ======================================================================
log2FC_mono_plantID <- log2FC_mono_data %>% 
  mutate(sample_ID = paste0(donor_lab, "_", recipient_lab, "_", plant_ID),
         rlab_plant_ID = paste0(recipient_lab, "_", plant_ID),
         .before=borneol_EIC) %>%
  dplyr::select(-unique_ID) %>% 
  pivot_longer(cols = contains("EIC"),
               names_to = "compound",
               values_to = "value") %>% 
  mutate(donor_lab = as.factor(donor_lab),
         compound = as.factor(compound),
         recipient_lab = as.factor(recipient_lab),
         sample_ID = as.factor(sample_ID)) %>% 
  pivot_wider(names_from = treatment,
              values_from = value) %>% 
  group_by(compound, sample_ID) %>%  # 95 x 6 = 570
  mutate(log2FC = log2(JA/C)) %>% 
  ungroup(); log2FC_mono_plantID
# normalize data
set.seed(123); fitdistrplus::descdist(log2FC_mono_plantID$log2FC, discrete=F, boot=999)
set.seed(123); log2FC_mono_norm <- bestNormalize::bestNormalize(log2FC_mono_plantID$log2FC); log2FC_mono_norm
log2FC_mono_plantID$log2FC_norm <- log2FC_mono_norm$x.t
set.seed(123); fitdistrplus::descdist(log2FC_mono_plantID$log2FC_norm, discrete=F, boot=999)
# apply mixed-effects model
mono_mod <- lmer(log2FC_norm ~ donor_lab * compound +
                   (1|rlab_plant_ID) +
                   (1|recipient_lab),
                 data = log2FC_mono_plantID,
                 REML=F); mono_mod@optinfo$conv$lme4$messages
# check model
check_model(mono_mod)
# update model without interaction
mono_mod_noInt <- update(mono_mod, .~. -donor_lab:compound); summary(mono_mod_noInt)
# check interaction with likelihood-ratio test
mono_mod_LRT <- anova(mono_mod, mono_mod_noInt, test="Chisq"); mono_mod_LRT
# generate p-values for LMM (type III because interaction is significant)
anova(mono_mod, type=3, ddf="Satterthwaite")

## MIXED =====================================================================
log2FC_mixed_plantID <- log2FC_mixed_data %>% 
  mutate(sample_ID = paste0(donor_lab, "_", recipient_lab, "_", plant_ID),
         rlab_plant_ID = paste0(recipient_lab, "_", plant_ID),
         .before=borneol_EIC) %>%
  dplyr::select(-unique_ID) %>% 
  pivot_longer(cols = contains("EIC"),
               names_to = "compound",
               values_to = "value") %>% 
  mutate(donor_lab = as.factor(donor_lab),
         compound = as.factor(compound),
         recipient_lab = as.factor(recipient_lab),
         sample_ID = as.factor(sample_ID)) %>% 
  pivot_wider(names_from = treatment,
              values_from = value) %>% 
  group_by(compound, sample_ID) %>%  # 95 x 6 = 570
  mutate(log2FC = log2(JA/C)) %>% 
  ungroup(); log2FC_mixed_plantID
# normalise
#set.seed(123); fitdistrplus::descdist(log2FC_mixed_plantID$log2FC, discrete=F, boot=999)
set.seed(123); log2FC_mixed_norm <- bestNormalize::bestNormalize(log2FC_mixed_plantID$log2FC); log2FC_mixed_norm
log2FC_mixed_plantID$log2FC_norm <- log2FC_mixed_norm$x.t
#set.seed(123); fitdistrplus::descdist(log2FC_mixed_plantID$log2FC_norm, discrete=F, boot=999)
# apply mixed-effects model
mixed_mod <- lmer(log2FC_norm ~ donor_lab * compound +
                   (1|rlab_plant_ID) +
                   (1|recipient_lab),
                 data = log2FC_mixed_plantID,
                 REML=F); mixed_mod@optinfo$conv$lme4$messages
# check model
check_model(mixed_mod)
# update model without interaction
mixed_mod_noInt <- update(mixed_mod, .~. -donor_lab:compound); summary(mixed_mod_noInt)
# check interaction with likelihood-ratio test
mixed_mod_LRT <- anova(mixed_mod, mixed_mod_noInt, test="Chisq"); mixed_mod_LRT
# generate p-values for LMM (type III because interaction is significant)
anova(mixed_mod, type=3, ddf="Satterthwaite")

# FIGURE 5 ====================================================================
## Subfigure 5A ===============================================================
subfigure_5A <- ggplot(mono_log2FCs_median,
                            aes(y=compounds, x=donor_lab)) +
  geom_tile(aes(fill = log2FC_median), colour = "white") +
  scale_fill_gradient2(low = "blue", mid="white", high = "red",
                       limits=c(-3.5, 8.5), breaks=seq(-3,8.5,by=2)) +
  labs(y="Volatile organic compound",
       x="Donor laboratory",
       title="Mono",
       subtitle=bquote(bold('LMM'['LRT']*":")~"Donor laboratory"%*%"Compounds"~
                         "["*italic(chi)^2~'='~
                         .(round(mono_mod_LRT$Chisq[2],2))*'] **'),
       fill=bquote(plain("log"[2])*"FC")) +
  scale_y_discrete(labels=rev(c("borneol",
                                "camphor",
                                bquote(italic(p)*"-cymene"),
                                "eucalyptol",
                                "limonene",
                                bquote(alpha*"-pinene"),
                                bquote(alpha*"-thujone"),
                                bquote(beta*"-thujone"))),
                   expand=c(0,0),
                   limits=rev) +
  scale_x_discrete(labels=c("L1", "L2", "L3", "L4", "L5"),
                   expand=c(0,0)) +
  theme_cowplot() +
  theme(panel.border=element_rect(fill=NA, color="black"),
        legend.position = "right",
        plot.subtitle = element_text(hjust = 0)); subfigure_5A

## Subfigure 5B ===============================================================
subfigure_5B <- ggplot(mixed_log2FCs_median,
                            aes(y=compounds, x=donor_lab)) +
  geom_tile(aes(fill = log2FC_median), colour = "white") +
  scale_fill_gradient2(low = "blue", mid="white", high = "red",
                       limits=c(-3.5, 8.5), breaks=seq(-3,8.5,by=2)) +
  labs(y=" ", #"Compound",
       x="Donor laboratory",
       title="Mixed",
       subtitle=bquote(bold('LMM'['LRT']*":")~"Donor laboratory"%*%"Compounds"~
                         "["*italic(chi)^2~'='~
                         .(round(mixed_mod_LRT$Chisq[2],2))*'] ***'),
       fill=bquote(plain("log"[2])*"FC")) +
  scale_y_discrete(labels=rev(c("borneol",
                                "camphor",
                                bquote(italic(p)*"-cymene"),
                                "eucalyptol",
                                "limonene",
                                bquote(alpha*"-pinene"),
                                bquote(alpha*"-thujone"),
                                bquote(beta*"-thujone"))),
                   expand=c(0,0),
                   limits=rev) +
  scale_x_discrete(labels=c("L1", "L2", "L3", "L4", "L5"),
                   expand=c(0,0)) +
  theme_cowplot() +
  theme(panel.border=element_rect(fill=NA, color="black"),
        legend.position = "none",
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=0, hjust=1)); subfigure_5B

## COMBINE PLOTS ==============================================================
hmaps_combined <- plot_grid(subfigure_5A + theme(legend.position="none"),
                            subfigure_5B + theme(legend.position="none"),
                            align = 'h',
                            labels = c("       A", "     B"),
                            label_size = 18,
                            rel_widths = c(1,.855),
                            hjust = .6,
                            nrow = 1)
# extract legend
legend_b <- get_legend(subfigure_5A + theme(legend.position="right"))
# add legend
plot_grid(hmaps_combined, legend_b,
          ncol = 2,
          rel_widths = c(1, .1))
# and save
ggsave(last_plot(),
       file=here::here(".", # replace with desired file path
                       "Eckert_et_al_2023_Fig5.tiff"),
       width=11, height=4,
       dpi=600,
       compression="lzw")

# SESSION INFO ================================================================
sessionInfo()

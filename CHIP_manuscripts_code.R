# CHIP and BioVU Project

library(tidyverse)
library(gtsummary)
library(lubridate)
library(survival)
library(survminer)

# Loading VUMC CHIP data for VUMC
calls <- readxl::read_xlsx("7492_CHIP_calls.xlsx")

# CHIP plot of frequency (includes patienst with mutiple CHIP mutations)
calls %>%
  filter(!is.na(Gene)) %>%
  ggplot(aes(x = fct_infreq(Gene), fill = Gene)) +
  geom_bar() +
  labs(
    x = "Gene",
    y = "Number of mutations",
    title = "CHIP in Heart Transplant Recipients (BioVU)"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Manuscript/VUMC_chip_mutations.png")

# Loading CHIP data from VUMC for general analyses
data <- read_csv("BioVU_SD_chip_il6r_combined_6.23.22.csv")

# Adding column for VAF >0.1
data <- data %>%
  mutate(
    vaf_greater_0.1 = case_when(
      VAF_1 >= 0.1 | VAF_2 >= 0.1 | VAF_3 >= 0.1 ~ 1,
      TRUE ~ 0
    ))

# Adding Yes/No column for homozygous IL6R protection
data <- data %>%
  mutate(
    il6r = case_when(
      Count_IL6R_D358A == 2 ~ 1,
      TRUE ~ 0
    )
  )

# Number of NAs in the dataset
View(data %>%
       select(everything()) %>%
       summarise_all(funs(sum(is.na(.)))))

# Calculate age at DNA sample collection (missing n = 257)
data <- data %>%
  mutate(
    dob = mdy(dob),
    dna_collection = mdy(dna_collection),
    date_dna = decimal_date(dna_collection) - decimal_date(dob)
  )

# Tables
table_1 <- data %>% select(has_CHIP, age_at_htx, gender, ischemic, htn, diabetes, statin, acr, amr, dsa, cav, deceased)

t1 <- tbl_summary(table_1,
                  by = has_CHIP,
                  missing = "no",
                  label = list(age_at_htx ~ "Age at transplant",
                               gender ~ "Gender (male)",
                               ischemic ~ "Ischemic etiology",
                               htn ~ "Hypertension",
                               diabetes ~ "Diabetes mellitus",
                               statin ~ "Statin use",
                               acr ~ "Acute cellular rejection",
                               amr ~ "Antibody-mediated rejection",
                               dsa ~ "Class II donor specific antibodies",
                               cav ~ "Cardiac allograft vasculopathy",
                               deceased ~ "Deceased")) %>%
  add_n() %>% add_p() %>%
  modify_header(label = "**Demographics by CHIP status**",
                stat_1 = "**No CHIP (N = 527)**",
                stat_2 = "**Has CHIP (N = 82)**") %>% bold_labels()
t1

gt::gtsave(as_gt(t1), file = "Manuscript/VUMC_chip_demographics.png")

# Logistic regression analyses
# Logistic regression - all CAV
mylogit <- glm(cav ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
               data = data, family = "binomial")

tbl_regression(mylogit,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(has_CHIP ~ "Presence of CHIP",
                            age_at_htx ~ "Age at transplant",
                            ischemic ~ "Ischemic etiology",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            amr ~ "History of AMR",
                            dsa ~ "Presence of class II DSAs")) %>%
  modify_header(label = "**CAV**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/VUMC_all_CAV_log_regression.png")

# Logistic regression - CAV 2/3
mylogit_cav_2_3 <- glm(cav_high_grade ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
                       data = data, family = "binomial")
tbl_regression(mylogit_cav_2_3,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(has_CHIP ~ "Presence of CHIP",
                            age_at_htx ~ "Age at transplant",
                            ischemic ~ "Ischemic etiology",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            amr ~ "History of AMR",
                            dsa ~ "Presence of class II DSAs")) %>%
  modify_header(label = "**CAV grade 2/3**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/VUMC_CAV_2_3_log_regression.png")

# Mortality logistic regression
mylogit_mortality <- glm(deceased ~ has_CHIP + age_at_htx + cav + htn + diabetes,
                         data = data, family = "binomial")
tbl_regression(mylogit_mortality,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(has_CHIP ~ "Presence of CHIP",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/VUMC_mortality_log_regression.png")

# DNMT3A vs others
dnmt3a <- data %>%
  mutate(
    dnmt3a_mut = case_when(
      Gene_1 == "DNMT3A" | Gene_2 == "DNMT3A" | Gene_3 == "DNMT3A" ~ 1,
      TRUE ~ 0
    ))

table_dnmt3a <- dnmt3a %>% select(dnmt3a_mut, age_at_htx, gender, ischemic, htn, diabetes, statin, acr, amr, dsa, cav, deceased)

t_dnmt3a <- tbl_summary(table_dnmt3a,
                  by = dnmt3a_mut,
                  missing = "no",
                  label = list(age_at_htx ~ "Age at transplant",
                               gender ~ "Gender (male)",
                               ischemic ~ "Ischemic etiology",
                               htn ~ "Hypertension",
                               diabetes ~ "Diabetes mellitus",
                               statin ~ "Statin use",
                               acr ~ "Acute cellular rejection",
                               amr ~ "Antibody-mediated rejection",
                               dsa ~ "Class II donor specific antibodies",
                               cav ~ "Cardiac allograft vasculopathy",
                               deceased ~ "Deceased")) %>%
  add_n() %>% add_p() %>%
  modify_header(label = "**Demographics by DNMT3A CHIP status**",
                stat_1 = "**No DNMT3A CHIP (N = 565)**",
                stat_2 = "**Has DNMT3A CHIP (N = 44)**") %>% bold_labels()

t_dnmt3a

gt::gtsave(as_gt(t_dnmt3a), file = "Manuscript/VUMC_DNMT3A_demographics.png")

dnmt3a_logit <- glm(cav ~ dnmt3a_mut + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
               data = dnmt3a, family = "binomial")

tbl_regression(dnmt3a_logit,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(dnmt3a_mut ~ "Presence of DNMT3A CHIP",
                            age_at_htx ~ "Age at transplant",
                            ischemic ~ "Ischemic etiology",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            amr ~ "History of AMR",
                            dsa ~ "Presence of class II DSAs")) %>%
  modify_header(label = "**CAV**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/VUMC_DNMT3A_log_regression.png")

dnmt3a_mortality <- glm(deceased ~ dnmt3a_mut + age_at_htx + cav + htn + diabetes,
                         data = dnmt3a, family = "binomial")
tbl_regression(dnmt3a_mortality,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(dnmt3a_mut ~ "Presence of DNMT3A CHIP",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/VUMC_DNMT3A_mortality_log_regression.png")

# VAF >0.1
vaf_mortality <- glm(deceased ~ vaf_greater_0.1 + age_at_htx + cav + htn + diabetes,
                        data = data, family = "binomial")

tbl_regression(vaf_mortality,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(vaf_greater_0.1 ~ "Presence of VAF >0.1",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/VUMC_VAF_0.1_mortality_log_regression.png")

vaf_cav <- glm(cav ~ vaf_greater_0.1 + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
               data = data, family = "binomial")

tbl_regression(vaf_cav,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(vaf_greater_0.1 ~ "Presence of VAF >0.1",
                            age_at_htx ~ "Age at transplant",
                            ischemic ~ "Ischemic etiology",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            amr ~ "History of AMR",
                            dsa ~ "Presence of class II DSAs")) %>%
  modify_header(label = "**CAV**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/VUMC_VAF_0.1_CAV_log_regression.png")

#################Columbia dataset

# Load data
cumc <- readxl::read_xlsx("CUMC Data Dictionary for Kaushik3.xlsx")

# CHIP mutations
col_chip <- readxl::read_xlsx("CUMC_chip_calls.xlsx")

col_chip %>%
  filter(!is.na(Gene)) %>%
  ggplot(aes(x = fct_infreq(Gene), fill = Gene)) +
  geom_bar() +
  labs(
    x = "Gene",
    y = "Number of mutations",
    title = "CHIP in Heart Transplant Recipients (Columbia)"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
cumc_muts

ggsave("Manuscript/Columbia_chip_mutations.png")

# Combine VUMC and Columbia mutations
combined_calls <- readxl::read_xlsx("combined_chip_calls.xlsx")

# CHIP plot of frequency (includes patienst with mutiple CHIP mutations)
combined_calls %>%
  filter(!is.na(Gene)) %>%
  ggplot(aes(x = fct_infreq(Gene), fill = Gene)) +
  geom_bar() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    x = "Gene",
    y = "Number of mutations",
    title = "CHIP in Heart Transplant Recipients (BioVU & Columbia)"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Manuscript/Combined_chip_mutations.png")


# Changing M = 1 and F = 0 for table

cumc$gender <- ifelse(cumc$gender=="M",1,0)

# Table 1 demographics
table_1 <- cumc %>% select(has_CHIP, age_at_htx, gender, ischemic, htn, diabetes, statin, TwoRplus, amr, dsa, cav, deceased)

t1 <- tbl_summary(table_1,
                  by = has_CHIP,
                  missing = "no",
                  label = list(age_at_htx ~ "Age at transplant",
                               gender ~ "Gender (male)",
                               ischemic ~ "Ischemic etiology",
                               htn ~ "Hypertension",
                               diabetes ~ "Diabetes mellitus",
                               statin ~ "Statin use",
                               TwoRplus ~ "Acute cellular rejection",
                               amr ~ "Antibody-mediated rejection",
                               dsa ~ "Donor specific antibodies",
                               cav ~ "Cardiac allograft vasculopathy",
                               deceased ~ "Deceased")) %>%
  add_n() %>% add_p() %>%
  modify_header(label = "**Demographics by CHIP status**",
                stat_1 = "**No CHIP (N = 143)**",
                stat_2 = "**Has CHIP (N = 35)**") %>% bold_labels()
t1
gt::gtsave(as_gt(t1), file = "Manuscript/CUMC_Table_1_demographics.png")

# Logistic regression - all CAV
mylogit_cumc <- glm(cav ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
               data = cumc, family = "binomial")

tbl_regression(mylogit_cumc,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(has_CHIP ~ "Presence of CHIP",
                            age_at_htx ~ "Age at transplant",
                            ischemic ~ "Ischemic etiology",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            amr ~ "History of AMR",
                            dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**CAV**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/CUMC_all_CAV_log_regression.png")

# Mortality logistic regression
mylogit_mortality_cumc <- glm(deceased ~ has_CHIP + age_at_htx + cav + htn + diabetes,
                         data = cumc, family = "binomial")
tbl_regression(mylogit_mortality_cumc,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(has_CHIP ~ "Presence of CHIP",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/CUMC_mortality_log_regression.png")
warnings()

# DNMT3A status
dnmt3a_logit_cumc <- glm(cav ~ DNMT3A + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
                    data = cumc, family = "binomial")

tbl_regression(dnmt3a_logit_cumc,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(DNMT3A ~ "Presence of DNMT3A CHIP",
                            age_at_htx ~ "Age at transplant",
                            ischemic ~ "Ischemic etiology",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            amr ~ "History of AMR",
                            dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**CAV**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/CUMC_DNMT3A_log_regression.png")

dnmt3a_mortality_cumc <- glm(deceased ~ DNMT3A + age_at_htx + cav + htn + diabetes,
                        data = cumc, family = "binomial")
tbl_regression(dnmt3a_mortality_cumc,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(DNMT3A ~ "Presence of DNMT3A CHIP",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/CUMC_DNMT3A_mortality_log_regression.png")
warnings()

# CUMC VAF >0.1

cumc_vaf_mortality <- glm(deceased ~ VAF + age_at_htx + cav + htn + diabetes,
                     data = cumc, family = "binomial")

tbl_regression(cumc_vaf_mortality,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(VAF ~ "Presence of VAF >0.1",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/Columbia_VAF_0.1_mortality_log_regression.png")

cumc_vaf_cav <- glm(cav ~ VAF + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
               data = cumc, family = "binomial")

tbl_regression(cumc_vaf_cav,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(VAF ~ "Presence of VAF >0.1",
                            age_at_htx ~ "Age at transplant",
                            ischemic ~ "Ischemic etiology",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            amr ~ "History of AMR",
                            dsa ~ "Presence of class II DSAs")) %>%
  modify_header(label = "**CAV**") %>%
  as_gt() %>%
  gt::gtsave("Manuscript/Columbia_VAF_0.1_CAV_log_regression.png")

####Fixed effect meta analysis
library(metafor)
library(meta)
library(dmetar)

power.analysis(OR = 1.36, k = 2, n1 = 609, n2 = 178, p = 0.05, heterogeneity = "fixed")

meta_cav <- read_csv("meta_cav.csv")
meta_mortality <- read_csv("meta_mortality.csv")

m.bin.cav <- metabin(Ee, Ne, Ec, Nc,
                     data = meta_cav,
                     studlab = paste(Group),
                     comb.fixed = T, comb.random = F,
                     method = "MH", sm = "OR")
m.bin.cav
forest(m.bin.cav, leftcols = c('studlab'))

png("Manuscript/Meta_CAV.png", width = 2400, height = 800, res = 300)
forest(m.bin.cav, leftcols = c('studlab'))
dev.off()


m.bin.mortality <- metabin(Ee, Ne, Ec, Nc,
                     data = meta_mortality,
                     studlab = paste(Group),
                     comb.fixed = T, comb.random = F,
                     method = "MH", sm = "OR")
m.bin.mortality
forest(m.bin.mortality, leftcols = c('studlab'))

png("Manuscript/Meta_mortality.png", width = 2400, height = 800, res = 300)
forest(m.bin.mortality, leftcols = c('studlab'))
dev.off()

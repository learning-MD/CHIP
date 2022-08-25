# IL6R and BioVU Project

library(tidyverse)
library(gtsummary)
library(lubridate)
library(survival)
library(survminer)

# Loading data from VUMC for general analyses
data <- read_csv("BioVU_SD_chip_il6r_combined_6.23.22.csv")

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

# Tables
table_1 <- data %>% select(il6r, age_at_htx, gender, ischemic, htn, diabetes, statin, acr, amr, dsa, cav, deceased)

t1 <- tbl_summary(table_1,
                  by = il6r,
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
  modify_header(label = "**Demographics by IL6R signaling status**",
                stat_1 = "**Preserved IL6R Signaling (N = 538)**",
                stat_2 = "**Has Impaired IL6R Signaling (N = 71)**") %>% bold_labels()
t1

gt::gtsave(as_gt(t1), file = "IL6R/VUMC_IL6R_demographics.png")

# Logistic regression analyses
# Logistic regression - all CAV
mylogit <- glm(cav ~ il6r + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
               data = data, family = "binomial")

tbl_regression(mylogit,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(il6r ~ "IL6R p.Asp358Ala",
                            age_at_htx ~ "Age at transplant",
                            ischemic ~ "Ischemic etiology",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            amr ~ "History of AMR",
                            dsa ~ "Presence of class II DSAs")) %>%
  modify_header(label = "**CAV**") %>%
  as_gt() %>%
  gt::gtsave("IL6R/VUMC_all_CAV_IL6R.png")

# Logistic regression - CAV 2/3
mylogit_cav_2_3 <- glm(cav_high_grade ~ il6r + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
                       data = data, family = "binomial")
tbl_regression(mylogit_cav_2_3,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(il6r ~ "IL6R p.Asp358Ala",
                            age_at_htx ~ "Age at transplant",
                            ischemic ~ "Ischemic etiology",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            amr ~ "History of AMR",
                            dsa ~ "Presence of class II DSAs")) %>%
  modify_header(label = "**CAV grade 2/3**") %>%
  as_gt() %>%
  gt::gtsave("IL6R/VUMC_CAV_2_3_IL6R.png")

# Mortality logistic regression
mylogit_mortality <- glm(deceased ~ il6r + age_at_htx + cav + htn + diabetes,
                         data = data, family = "binomial")
tbl_regression(mylogit_mortality,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(il6r ~ "IL6R p.Asp358Ala",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  as_gt() %>%
  gt::gtsave("IL6R/VUMC_mortality_IL6R.png")

# ACR logistic regression
mylogit_acr <- glm(acr ~ il6r + age_at_htx + ischemic + htn + diabetes + dsa,
               data = data, family = "binomial")
tbl_regression(mylogit_acr,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(il6r ~ "IL6R p.Asp358Ala",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            dsa ~ "Class II donor-specific antibodies")) %>%
  modify_header(label = "**Grade ??? 2R Acute Cellular Rejection**") %>%
  as_gt() %>%
  gt::gtsave("IL6R/VUMC_ACR_IL6R.png")

# AMR logistic regression
mylogit_amr <- glm(amr ~ il6r + age_at_htx + ischemic + htn + diabetes + acr + dsa,
                   data = data, family = "binomial")
tbl_regression(mylogit_amr,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(il6r ~ "IL6R p.Asp358Ala",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            ischemic ~ "Ischemic etiology",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            dsa ~ "Class II donor-specific antibodies")) %>%
  modify_header(label = "**Antibody-Mediated Rejection**") %>%
  as_gt() %>%
  gt::gtsave("IL6R/VUMC_AMR_IL6R.png")

#################Columbia dataset

# Load data
cumc <- readxl::read_xlsx("CUMC Data Dictionary for Kaushik3.xlsx")

# Changing M = 1 and F = 0 for table
cumc$gender <- ifelse(cumc$gender=="M",1,0)

# Create IL6R column
cumc <- cumc %>%
  mutate(
    il6r = case_when(
      Count_IL6R_D358A == 2 ~ 1,
      TRUE ~ 0
    )
  )

# Table 1 demographics
table_1 <- cumc %>% select(il6r, age_at_htx, gender, ischemic, htn, diabetes, statin, TwoRplus, amr, dsa, cav, deceased)

t1 <- tbl_summary(table_1,
                  by = il6r,
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
  modify_header(label = "**Demographics by IL6R signaling status**",
                stat_1 = "**Preserved IL6R Signaling (N = 152)**",
                stat_2 = "**Has Impaired IL6R Signaling (N = 26)**") %>% bold_labels()
t1
gt::gtsave(as_gt(t1), file = "IL6R/Columbia_Table_1_demographics_IL6R.png")

# Logistic regression - all CAV
mylogit_cumc <- glm(cav ~ il6r + age_at_htx + ischemic + htn + diabetes + TwoRplus + amr + dsa,
                    data = cumc, family = "binomial")

tbl_regression(mylogit_cumc,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(il6r ~ "IL6R p.Asp358Ala",
                            age_at_htx ~ "Age at transplant",
                            ischemic ~ "Ischemic etiology",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            TwoRplus ~ "History of ACR",
                            amr ~ "History of AMR",
                            dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**CAV**") %>%
  as_gt() %>%
  gt::gtsave("IL6R/Columbia_all_CAV_IL6R.png")

# Mortality logistic regression
mylogit_mortality_cumc <- glm(deceased ~ il6r + age_at_htx + cav + htn + diabetes,
                              data = cumc, family = "binomial")
tbl_regression(mylogit_mortality_cumc,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(il6r ~ "IL6R p.Asp358Ala",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  as_gt() %>%
  gt::gtsave("IL6R/Columbia_mortality_IL6R.png")

# ACR logistic regression
mylogit_acr_cumc <- glm(acr ~ il6r + age_at_htx + ischemic + htn + diabetes + dsa,
                   data = cumc, family = "binomial")
tbl_regression(mylogit_acr_cumc,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(il6r ~ "IL6R p.Asp358Ala",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**Grade ??? 2R Acute Cellular Rejection**") %>%
  as_gt() %>%
  gt::gtsave("IL6R/Columbia_ACR_IL6R.png")

# AMR logistic regression
mylogit_amr_cumc <- glm(amr ~ il6r + age_at_htx + ischemic + htn + diabetes + acr + dsa,
                   data = cumc, family = "binomial")
tbl_regression(mylogit_amr_cumc,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(il6r ~ "IL6R p.Asp358Ala",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            ischemic ~ "Ischemic etiology",
                            diabetes ~ "Diabetes mellitus",
                            acr ~ "History of ACR",
                            dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**Antibody-Mediated Rejection**") %>%
  as_gt() %>%
  gt::gtsave("IL6R/Columbia_AMR_IL6R.png")

####Fixed effect meta analysis
library(metafor)
library(meta)
library(dmetar)

power.analysis(OR = 1.36, k = 2, n1 = 609, n2 = 178, p = 0.05, heterogeneity = "fixed")

meta_cav_il6r <- read_csv("meta_il6r_cav.csv")
meta_mortality_il6r <- read_csv("meta_il6r_mortality.csv")
meta_acr_il6r <- read_csv("meta_il6r_acr.csv")
meta_amr_il6r <- read_csv("meta_il6r_amr.csv")

m.bin.cav.il6r <- metabin(Ee, Ne, Ec, Nc,
                     data = meta_cav_il6r,
                     studlab = paste(Group),
                     comb.fixed = T, comb.random = F,
                     method = "MH", sm = "OR")
m.bin.cav.il6r
forest(m.bin.cav.il6r, leftcols = c('studlab'))

png("IL6R/Meta_CAV_IL6R.png", width = 2400, height = 800, res = 300)
forest(m.bin.cav.il6r, leftcols = c('studlab'))
dev.off()


m.bin.mortality_il6r <- metabin(Ee, Ne, Ec, Nc,
                           data = meta_mortality_il6r,
                           studlab = paste(Group),
                           comb.fixed = T, comb.random = F,
                           method = "MH", sm = "OR")
m.bin.mortality_il6r
forest(m.bin.mortality_il6r, leftcols = c('studlab'))

png("IL6R/Meta_mortality_IL6R.png", width = 2400, height = 800, res = 300)
forest(m.bin.mortality_il6r, leftcols = c('studlab'))
dev.off()

m.bin.acr_il6r <- metabin(Ee, Ne, Ec, Nc,
                                data = meta_acr_il6r,
                                studlab = paste(Group),
                                comb.fixed = T, comb.random = F,
                                method = "MH", sm = "OR")
m.bin.acr_il6r
forest(m.bin.acr_il6r, leftcols = c('studlab'))

png("IL6R/Meta_ACR_IL6R.png", width = 2400, height = 800, res = 300)
forest(m.bin.acr_il6r, leftcols = c('studlab'))
dev.off()

m.bin.amr_il6r <- metabin(Ee, Ne, Ec, Nc,
                          data = meta_amr_il6r,
                          studlab = paste(Group),
                          comb.fixed = T, comb.random = F,
                          method = "MH", sm = "OR")
m.bin.amr_il6r
forest(m.bin.amr_il6r, leftcols = c('studlab'))

png("IL6R/Meta_AMR_IL6R.png", width = 2400, height = 800, res = 300)
forest(m.bin.amr_il6r, leftcols = c('studlab'))
dev.off()

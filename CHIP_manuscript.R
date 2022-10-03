library(mice)
library(survival)
library(tidyverse)
library(gtsummary)
library(sjmisc)
library(survminer)
library(metafor)
library(meta)
library(dmetar)

# CHIP plot of frequency (includes patienst with mutiple CHIP mutations)

combined_calls <- readxl::read_xlsx("_.xlsx")

combined_calls %>%
  filter(!is.na(Gene)) %>%
  ggplot(aes(x = fct_infreq(Gene), fill = Gene)) +
  geom_bar() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    x = "Gene",
    y = "Number of mutations",
    title = "CHIP in Heart Transplant Recipients (Vanderbilt & Columbia)"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_discrete(limits = c("DNMT3A", "PPM1D", "TET2",
                                 "ASXL1", "TP53", "JAK2",
                                 "GNAS", "BRCC3", "CBL",
                                 "SF3B1", "SRSF2"))
ggsave("Combined_chip_mutations.png")


# Loading data from VUMC for general analyses
data <- read_csv("_.csv")

# Adding Yes/No column for homozygous IL6R protection
data <- data %>%
  mutate(
    il6r = case_when(
      Count_IL6R_D358A == 2 ~ 1,
      TRUE ~ 0
    )
  )

# Adding VAF and DNMT3A info
data <- data %>%
  mutate(
    vaf_greater_0.1 = case_when(
      VAF_1 >= 0.1 | VAF_2 >= 0.1 | VAF_3 >= 0.1 ~ 1,
      TRUE ~ 0
    ))

data <- data %>%
  mutate(
    dnmt3a_mut = case_when(
      Gene_1 == "DNMT3A" | Gene_2 == "DNMT3A" | Gene_3 == "DNMT3A" ~ 1,
      TRUE ~ 0
    ))

# MICE for imputations
imp <- mice(list, m = 16)
saveRDS(imp, "imputations_mice.rds")
imp <- readRDS("imputations_mice.rds")

################## Tables and imputed/complete case analyses
# Table 1 - complete case analysis
table_1 <- data %>% select(has_CHIP, age_at_htx, gender, ischemic, htn, diabetes, statin, acr, amr, dsa, cav, deceased)

t1.imp <- tbl_summary(table_1,
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
  add_p() %>%
  modify_header(label = "**Demographics by CHIP status**",
                stat_1 = "**No CHIP (N = 527)**",
                stat_2 = "**Has CHIP (N = 82)**") %>% bold_labels()

# Imputed demographics - Table 1
datalist <- complete(imp, action ="all")
mat <- complete(imp)

table_2 <- mat %>% select(has_CHIP, age_at_htx, gender, ischemic, htn, diabetes, statin, acr, amr, dsa, cav, deceased)

t2 <- tbl_summary(table_2,
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
  add_p() %>%
  modify_header(label = "**Demographics by CHIP status**",
                stat_1 = "**No CHIP**",
                stat_2 = "**Has CHIP**") %>% bold_labels()
# All CAV - complete case analysis
mylogit <- glm(cav ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
               data = list, family = "binomial")

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
  modify_header(label = "**CAV**")

# All CAV - imputed model analysis
mi.model.logit.cav <- imp %>%
  with(glm(cav ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa, family = "binomial"))

tbl_regression(mi.model.logit.cav,
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
  modify_header(label = "**CAV**")

# All mortality - complete case analysis
mylogit_mortality <- glm(deceased ~ has_CHIP + age_at_htx + cav + htn + diabetes,
                         data = list, family = "binomial")
tbl_regression(mylogit_mortality,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(has_CHIP ~ "Presence of CHIP",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**")

# Mortality - imputed model analysis

mi.model.logit.mort <- imp %>%
  with(glm(deceased ~ has_CHIP + age_at_htx + cav + htn + diabetes, family = "binomial"))

tbl_regression(mi.model.logit.mort,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(has_CHIP ~ "Presence of CHIP",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**")

# DNMt3A - complete case analysis
dnmt3a_logit <- glm(cav ~ dnmt3a_mut + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
                    data = list, family = "binomial")

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
  modify_header(label = "**CAV**")

# DNMT3A - imputed model analysis
mi.model.logit.dnmt3a <- imp %>%
  with(glm(cav ~ dnmt3a_mut + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa, family = "binomial"))

tbl_regression(mi.model.logit.dnmt3a,
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
  modify_header(label = "**CAV**")

mi.model.logit.dnmt3a.mort <- imp %>%
  with(glm(deceased ~ dnmt3a_mut + age_at_htx + cav + htn + diabetes, family = "binomial"))

tbl_regression(mi.model.logit.dnmt3a.mort,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(dnmt3a_mut ~ "Presence of DNMT3A CHIP",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**")

# VAF >0.1 - complete case analysis
vaf_mortality <- glm(deceased ~ vaf_greater_0.1 + age_at_htx + cav + htn + diabetes,
                     data = list, family = "binomial")

tbl_regression(vaf_mortality,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(vaf_greater_0.1 ~ "Presence of VAF >0.1",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**")

vaf_cav <- glm(cav ~ vaf_greater_0.1 + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
               data = list, family = "binomial")

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
  modify_header(label = "**CAV**")

# VAC >0.1 - imputed model analysis
mi.model.logit.vaf.cav <- imp %>%
  with(glm(cav ~ vaf_greater_0.1 + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa, family = "binomial"))

tbl_regression(mi.model.logit.vaf.cav,
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
  modify_header(label = "**CAV**")

mi.model.logit.vaf.mort <- imp %>%
  with(glm(deceased ~ vaf_greater_0.1 + age_at_htx + cav + htn + diabetes, family = "binomial"))

tbl_regression(mi.model.logit.vaf.mort,
               exponentiate = TRUE,
               pvalue_fun = ~style_pvalue(.x, digits = 2),
               label = list(vaf_greater_0.1 ~ "Presence of VAF >0.1",
                            age_at_htx ~ "Age at transplant",
                            htn ~ "Hypertension",
                            diabetes ~ "Diabetes mellitus",
                            cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**")

# Grade 2/3 CAV - complete case analysis
mod.cav <- glm(cav_high_grade ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
               data = list, family = "binomial")

tbl_regression(mod.cav,
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
  modify_header(label = "**Moderate-to-severe CAV**")

# Grade 2/3 CAV - imputed model
mod.cav.imp <- imp %>%
  with(glm(cav_high_grade ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
           family = "binomial"))

tbl_regression(mod.cav.imp,
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
  modify_header(label = "**CAV**")
#################Columbia Dataset
# Load data
cumc <- readxl::read_xlsx("_.xlsx")

# Changing M = 1 and F = 0 for table
cumc$gender <- ifelse(cumc$gender=="M",1,0)

# Table 1 demographics
cumc.table_1 <- cumc %>% select(has_CHIP, age_at_htx, gender, ischemic, htn, diabetes, statin, TwoRplus, amr, dsa, cav, deceased)
cumc.table_1 <- cumc.table_1 %>% rename(acr = TwoRplus)

cumc.t1 <- tbl_summary(cumc.table_1,
                       by = has_CHIP,
                       missing = "no",
                       label = list(age_at_htx ~ "Age at transplant (years)",
                                    gender ~ "Gender (male)",
                                    ischemic ~ "Ischemic etiology",
                                    htn ~ "Hypertension",
                                    diabetes ~ "Diabetes mellitus",
                                    statin ~ "Statin use",
                                    acr ~ "Acute cellular rejection",
                                    amr ~ "Antibody-mediated rejection",
                                    dsa ~ "Donor specific antibodies",
                                    cav ~ "Cardiac allograft vasculopathy",
                                    deceased ~ "Deceased")) %>%
  add_p() %>%
  modify_header(label = "**Demographics by CHIP status**",
                stat_1 = "**No CHIP (N = 143)**",
                stat_2 = "**Has CHIP (N = 35)**") %>% bold_labels()

# CHIP and CAV
mylogit_cumc.cav <- glm(cav ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
                        data = cumc, family = "binomial")

exp(cbind("Odds ratio" = coef(mylogit_cumc.cav),
          confint.default(mylogit_cumc.cav, level = 0.95)))

cumc.cav <- tbl_regression(mylogit_cumc.cav,
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
  modify_table_styling(
    column = c(p.value),
    hide = TRUE
  )

# Mortality logistic regression
mylogit_mortality_cumc <- glm(deceased ~ has_CHIP + age_at_htx + cav + htn + diabetes,
                              data = cumc, family = "binomial")

cumc.mort <- tbl_regression(mylogit_mortality_cumc,
                            exponentiate = TRUE,
                            pvalue_fun = ~style_pvalue(.x, digits = 2),
                            label = list(has_CHIP ~ "Presence of CHIP",
                                         age_at_htx ~ "Age at transplant",
                                         htn ~ "Hypertension",
                                         diabetes ~ "Diabetes mellitus",
                                         cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  modify_table_styling(
    column = c(p.value),
    hide = TRUE
  )

# DNMT3A status and CAV
dnmt3a_logit_cumc <- glm(cav ~ DNMT3A + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
                         data = cumc, family = "binomial")

cumc.dnmt3a.cav <- tbl_regression(dnmt3a_logit_cumc,
                                  exponentiate = TRUE,
                                  pvalue_fun = ~style_pvalue(.x, digits = 2),
                                  label = list(DNMT3A ~ "Presence of DNMT3A mutation",
                                               age_at_htx ~ "Age at transplant",
                                               ischemic ~ "Ischemic etiology",
                                               htn ~ "Hypertension",
                                               diabetes ~ "Diabetes mellitus",
                                               acr ~ "History of ACR",
                                               amr ~ "History of AMR",
                                               dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**CAV**") %>%
  modify_table_styling(
    column = c(p.value),
    hide = TRUE
  )

# DNMT3A and mortality
dnmt3a_mortality_cumc <- glm(deceased ~ DNMT3A + age_at_htx + cav + htn + diabetes,
                             data = cumc, family = "binomial")

cumc.dnmt3a.mort <- tbl_regression(dnmt3a_mortality_cumc,
                                   exponentiate = TRUE,
                                   pvalue_fun = ~style_pvalue(.x, digits = 2),
                                   label = list(DNMT3A ~ "Presence of DNMT3A mutation",
                                                age_at_htx ~ "Age at transplant",
                                                htn ~ "Hypertension",
                                                diabetes ~ "Diabetes mellitus",
                                                cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  modify_table_styling(
    column = c(p.value),
    hide = TRUE
  )

# CUMC VAF >0.1
cumc_vaf_mortality <- glm(deceased ~ VAF + age_at_htx + cav + htn + diabetes,
                          data = cumc, family = "binomial")

cumc.vaf.mort <- tbl_regression(cumc_vaf_mortality,
                                exponentiate = TRUE,
                                pvalue_fun = ~style_pvalue(.x, digits = 2),
                                label = list(VAF ~ "Presence of VAF >0.1",
                                             age_at_htx ~ "Age at transplant",
                                             htn ~ "Hypertension",
                                             diabetes ~ "Diabetes mellitus",
                                             cav ~ "Cardiac allograft vasculopathy")) %>%
  modify_header(label = "**Mortality**") %>%
  modify_table_styling(
    column = c(p.value),
    hide = TRUE
  )

cumc_vaf_cav <- glm(cav ~ VAF + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
                    data = cumc, family = "binomial")

cumc.vaf.cav <- tbl_regression(cumc_vaf_cav,
                               exponentiate = TRUE,
                               pvalue_fun = ~style_pvalue(.x, digits = 2),
                               label = list(VAF ~ "Presence of VAF >0.1",
                                            age_at_htx ~ "Age at transplant",
                                            ischemic ~ "Ischemic etiology",
                                            htn ~ "Hypertension",
                                            diabetes ~ "Diabetes mellitus",
                                            acr ~ "History of ACR",
                                            amr ~ "History of AMR",
                                            dsa ~ "Presence of DSAs")) %>%
  modify_header(label = "**CAV**") %>%
  modify_table_styling(
    column = c(p.value),
    hide = TRUE
  )

# Grade 2/3 CAV
cumc.mod.cav <- glm(cav_high_grade ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
                    data = cumc, family = "binomial")

cumc.modcav <- tbl_regression(cumc.mod.cav,
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
  modify_header(label = "**Moderate-to-severe CAV**") %>%
  modify_table_styling(
    column = c(p.value),
    hide = TRUE
  )


# Stacking tables - All CAV

mi.model.logit.cav <- imp %>%
  with(glm(cav ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa, family = "binomial"))

mi.model.cav <- tbl_regression(mi.model.logit.cav,
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
  modify_header(label = "**CAV**")

mylogit_cumc <- glm(cav ~ has_CHIP + age_at_htx + ischemic + htn + diabetes + acr + amr + dsa,
                    data = cumc, family = "binomial")

cumc.cav <- tbl_regression(mylogit_cumc,
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
  modify_header(label = "**CAV**")

tbl_stack_cav <-
  tbl_stack(list(mi.model.cav, cumc.cav), group_header = c("VUMC/BioVU Cohort", "CUIMC Cohort"))

tbl_merge_ex1 <-
  tbl_merge(
    tbls = list(mi.model.cav, cumc.cav),
    tab_spanner = c("**VUMC/BioVU Cohort**", "**CUIMC Cohort**")
  )

gt::gtsave(as_gt(tbl_merge_ex1), file = "Imputation/Stacked_CAV.png")

####Fixed effect meta analysis
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

###### Survival analysis at VUMC
# Kaplan-Meier survival analysis
dates <- data %>% select(grid, age_at_htx, diabetes, htn, ischemic, acr, amr, dsa, has_CHIP, il6r, vaf_greater_0.1, dob, date_of_death, htx_date, cav, cav_grade, cav_mention, last_fup)

dates$surv_death_status <- ifelse(is.na(dates$date_of_death), 1, 2)
dates$surv_death_fup_dates <- coalesce(dates$date_of_death, dates$last_fup)

km_analysis <- dates %>%
  dplyr::mutate(
    dob = as.Date(dob, format = "%m/%d/%Y"),
    date_of_death = as.Date(date_of_death, format = "%m/%d/%Y"),
    htx_date = as.Date(htx_date, format = "%m/%d/%Y"),
    cav_mention = as.Date(cav_mention, format = "%m/%d/%Y"),
    last_fup = as.Date(last_fup, format = "%m/%d/%Y"),
    surv_death_fup_dates = as.Date(surv_death_fup_dates, format = "%m/%d/%Y")
  )

# Time to death post-HTx
km_analysis <- km_analysis %>%
  mutate(
    os_days =
      as.numeric(
        difftime(surv_death_fup_dates,
                 htx_date,
                 units = "days"))
  )

# Survival curves for MORTALITY
fit <- survfit(Surv(os_days, surv_death_status) ~ has_CHIP, data = km_analysis)
print(fit)

survp <- ggsurvplot(fit,
                    pval = TRUE, conf.int = FALSE,
                    risk.table = TRUE,
                    risk.table.col = "strata",
                    linetype = "strata",
                    surv.median.line = "hv",
                    ggtheme = theme_bw(),
                    palette = c("red", "blue"),
                    legend.labs = c("No CHIP", "Has CHIP"),
                    size = 2)

survp

surv_cox <- coxph(
  Surv(os_days, surv_death_status) ~ has_CHIP,
  data = km_analysis
)
surv_cox
summary(surv_cox)
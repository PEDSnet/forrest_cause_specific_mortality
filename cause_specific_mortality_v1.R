# =============================================================================
# Cause-Specific Mortality of Children in the US versus Peer Nations, 1954–2023 
# Analysis Script
# Authors: Forrest CB, Koenigsberg LJ, Zelinski AM, Noble MK, Wolf ER
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================

library(dplyr)
library(purrr)
library(tidyr)
library(tidyverse)
library(scales)
library(patchwork)
library(fixest)
library(broom)
library(ggplot2)
library(ggh4x)

# ── File paths ────────────────────────────────────────────────────────────────
# Edit these to match your local directory structure.
PATH_HMD <- "data/raw/HMD"
PATH_WHO_ICD7 <- "data/raw/WHO/Morticd7.csv"
PATH_WHO_ICD8 <- "data/raw/WHO/Morticd8.csv"
PATH_WHO_ICD9 <- "data/raw/WHO/Morticd9.csv"
PATH_WHO_ICD10 <- list(
  part1 = "data/raw/WHO/Morticd10_part1.csv",
  part2 = "data/raw/WHO/Morticd10_part2.csv",
  part3 = "data/raw/WHO/Morticd10_part3.csv",
  part4 = "data/raw/WHO/Morticd10_part4.csv",
  part5 = "data/raw/WHO/Morticd10_part5.csv",
  part6 = "data/raw/WHO/Morticd10_part6.csv"
)
PATH_CODESETS <- "data/raw/codesets/harmonized_codesets_forrest_v1.xlsx"
PATH_OUTDIR <- "output"

JAMA_COLORS <- c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF")

dir.create(file.path(PATH_OUTDIR, "figures"),
           showWarnings = FALSE,
           recursive = TRUE
)
dir.create(file.path(PATH_OUTDIR, "tables"),
           showWarnings = FALSE,
           recursive = TRUE
)

# ── Age group canonical encodings ─────────────────────────────────────────────
AG_LT1 <- "Age_lt1" # <1 year
AG_1TO9 <- "Age_1to9" # 1-9 years
AG_10T14 <- "Age_10to14" # 10-14 years
AG_15T19 <- "Age_15to19" # 15-19 years
AG_0T19 <- "Age_0to19" # overall 0-19

# Canonical display order, overall first, then ascending age
AGE_ORDER <- c(AG_0T19, AG_LT1, AG_1TO9, AG_10T14, AG_15T19)

# Domain map (abbreviated cause -> domain)
DOMAIN_MAP <- tribble(
  ~cause_abbrev, ~domain, ~cluster,
  "RespInf", "Infections", "Respiratory Infections",
  "TB", "Infections", "Tuberculosis",
  "VaxPrev", "Infections", "Vaccine-Preventable Infections",
  "DiarrInf", "Infections", "Diarrheal/gastrointestinal infections",
  "InfOther", "Infections", "Infections, Other",
  "UnintInj", "Injury", "Unintentional Injury",
  "Homicide", "Injury", "Homicide",
  "Suicide", "Injury", "Suicide",
  "Anemia", "Non-Communicable Diseases", "Anemia",
  "Cancer", "Non-Communicable Diseases", "Cancer",
  "CardioVasc", "Non-Communicable Diseases", "Cardiovascular disorders",
  "CongenAnom", "Non-Communicable Diseases", "Congenital Anomalies",
  "Diabetes", "Non-Communicable Diseases", "Diabetes Mellitus",
  "Epilepsy", "Non-Communicable Diseases", "Epilepsy",
  "Neuro", "Non-Communicable Diseases", "Neuro-psychiatric disorders",
  "GI-Hepatic", "Non-Communicable Diseases", "Gastrointestinal/Hepatic disorders, Non-Infectious",
  "GenitoUrin", "Non-Communicable Diseases", "Genitourinary disorders, Non-Infectious",
  "NCD-Other", "Non-Communicable Diseases", "Non-Communicable Diseases, Other",
  "Perinatal", "Non-Communicable Diseases", "Perinatal disorders",
  "BirthAsphyx", "Non-Communicable Diseases", "Birth asphyxia/birth trauma",
  "SUID", "Non-Communicable Diseases", "Sudden Unexpected Infant Deaths (SUID)",
  "IllDef", "Other", "Ill-defined and unassigned",
  "MisCode", "Other", "Miscoded"
)

DOMAIN_ORDER <- c(
  "Infections", "Injury", "Non-Communicable Diseases", "Other", "Total"
)

# ── Country lookup ────────────────────────────────────────────────────────────
COUNTRY_LOOKUP <- tribble(
  ~pop_name,  ~who_code, ~country_name,                 ~region,
  "AUS",      5020L,     "Australia",                   "OECD13",
  "BEL",      4020L,     "Belgium",                     "OECD13",
  "CAN",      2090L,     "Canada",                      "OECD13",
  "DNK",      4050L,     "Denmark",                     "OECD13",
  "ESP",      4280L,     "Spain",                       "OECD13",
  "FIN",      4070L,     "Finland",                     "OECD13",
  "FRATNP",   4080L,     "France",                      "OECD13",
  "ITA",      4180L,     "Italy",                       "OECD13",
  "NLD",      4210L,     "Netherlands",                 "OECD13",
  "NOR",      4220L,     "Norway",                      "OECD13",
  "SWE",      4290L,     "Sweden",                      "OECD13",
  "CHE",      4300L,     "Switzerland",                 "OECD13",
  "GRB_NP",   4308L,     "United Kingdom",              "OECD13",
  "USA",      2450L,     "United States",               "US"
)

ALL_WHO_CODES <- COUNTRY_LOOKUP$who_code

# ── Study constants ───────────────────────────────────────────────────────────
STUDY_YEARS <- 1954:2023

# Approximate ICD transition years for figure reference lines
ICD_TRANSITIONS <- c(1968L, 1979L, 1999L)


# APC interval definitions (DAP Section 5.2.2)
APC_INTERVALS <- list(
  "1954-1974" = 1954:1974,
  "1975-1990" = 1975:1990,
  "1991-2005" = 1991:2005,
  "2006-2023" = 2006:2023
)

# Vertical lines for ICD transitions
icd_vlines <- function() {
  geom_vline(
    xintercept = ICD_TRANSITIONS, linetype = "dashed",
    color = "#80796BFF", linewidth = 0.5
  )
}

# Shared theme for all figures
theme_mortality <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(size = base_size - 1, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text        = element_text(size = base_size - 1)
    )
}

# =============================================================================
# DATA INGESTION
# =============================================================================

# ── HMD: Exposures (population) and Deaths ────────────────────────────────────
read_hmd_file <- function(file_name,
                          hmd_path = PATH_HMD) {
  
  country_dirs <- list.dirs(
    hmd_path,
    recursive = FALSE,
    full.names = TRUE
  )
  
  map_dfr(country_dirs, function(country_dir) {
    
    read_table(
      file.path(country_dir, file_name),
      skip = 2,
      show_col_types = FALSE
    ) %>%
      mutate(Age = as.integer(Age),
             PopName = basename(country_dir)) %>% 
      filter(!is.na(Age),
             Age <= 19,
             Year %in% STUDY_YEARS)
  })
}

deaths_all <- read_hmd_file("Deaths_1x1.txt")
exposures_all <- read_hmd_file("Exposures_1x1.txt")

# ── WHO-MD: cause-specific deaths ─────────────────────────────────────────────
read_who_file <- function(path) {
  read_csv(path, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
    mutate(
      Country = as.integer(Country),
      Year    = as.integer(Year),
      Sex     = as.integer(Sex),
      across(starts_with("Deaths"), as.numeric)
    ) %>% 
    filter(Country %in% ALL_WHO_CODES, Year %in% STUDY_YEARS, Sex %in% c(1:2))
}

who_icd7 <- read_who_file(PATH_WHO_ICD7) %>% mutate(icd_version = "7")
who_icd8 <- read_who_file(PATH_WHO_ICD8) %>% mutate(icd_version = "8")
who_icd9 <- read_who_file(PATH_WHO_ICD9) %>% mutate(icd_version = "9")
who_icd10 <- map_dfr(PATH_WHO_ICD10, read_who_file) %>%
  mutate(icd_version = "10")

who_all <- bind_rows(who_icd7, who_icd8, who_icd9, who_icd10)

# ── Forrest codesets ──────────────────────────────────────────────────────────
fk_codesets_orig <- readxl::read_excel(PATH_CODESETS, sheet = "Codeset", col_types = 'text') %>%
  mutate(
    icd_version  = as.character(`ICD version`),
    icd_code     = as.character(`WHO / ICD Code`),
    condition_cluster = as.character(`Condition Cluster`),
    broad_category = as.character(`Broad Category`),
    injury_cause = `Injury Cause`,
    age_specific = `Logic for only specific codes`
  ) %>%
  dplyr::select(icd_version, icd_code, condition_cluster, broad_category, injury_cause, age_specific) %>% 
  left_join(dplyr::select(DOMAIN_MAP, cause_abbrev, cluster), by = c('condition_cluster' = 'cluster'))

codeset_check <- fk_codesets_orig %>%
  filter(icd_version == 9, str_length(icd_code) > 3) %>%
  mutate(icd_code = substr(icd_code, 1, 3)) %>%
  distinct()
  
fk_codesets <-
  fk_codesets_orig %>%
  anti_join(codeset_check, by = c('icd_code', 'icd_version'))

# =============================================================================
# DATA PROCESSING
# =============================================================================
# ── HMD: assign canonical age group labels, reshape long, aggregate ───────────
assign_age_group <- function(age) {
  dplyr::case_when(
    age == 0 ~ AG_LT1,
    age >= 1 & age <= 9 ~ AG_1TO9,
    age >= 10 & age <= 14 ~ AG_10T14,
    age >= 15 & age <= 19 ~ AG_15T19,
    TRUE ~ NA_character_
  )
}

process_hmd <- function(df, value_col) {
  df %>%
    pivot_longer(c(Female, Male), names_to = "Sex", values_to = value_col) %>%
    mutate(
      Age_Group = assign_age_group(Age),
      Sex       = factor(Sex, levels = c("Male", "Female"))
    ) %>%
    filter(!is.na(Age_Group)) %>%
    group_by(PopName, Year, Sex, Age_Group) %>%
    summarise(!!value_col := sum(.data[[value_col]], na.rm = TRUE),
              .groups = "drop"
    ) %>%
    inner_join(COUNTRY_LOOKUP %>% dplyr::select(pop_name),
               by = c("PopName" = "pop_name")
    )
}

hmd_pop <- process_hmd(exposures_all, "Population")
hmd_deaths <- process_hmd(deaths_all, "Deaths_AllCause")

hmd_intermed <- hmd_pop %>% 
  left_join(hmd_deaths, by = c("PopName", "Year", "Age_Group", "Sex")) %>% 
  rename(Country = PopName)

# ── WHO-MD: aggregate Deaths columns to canonical age groups ──────────────────
# WHO column mapping:
#   Deaths2           -> Age_lt1
#   Deaths3-Deaths7   -> Age_1to9  (ages 1,2,3,4, and 5-9)
#   Deaths8           -> Age_10to14
#   Deaths9           -> Age_15to19
who_age_agg <- who_all %>%
  mutate(
    !!AG_LT1 := Deaths2,
    !!AG_1TO9 := rowSums(
      across(c(Deaths3, Deaths4, Deaths5, Deaths6, Deaths7)),
      na.rm = TRUE
    ),
    !!AG_10T14 := Deaths8,
    !!AG_15T19 := Deaths9
  ) %>%
  dplyr::select(
    Country, Year, Sex, Cause, icd_version,
    all_of(c(AG_LT1, AG_1TO9, AG_10T14, AG_15T19))
  )

# Reshape to long
who_long <- who_age_agg %>%
  pivot_longer(
    cols      = all_of(c(AG_LT1, AG_1TO9, AG_10T14, AG_15T19)),
    names_to  = "Age_Group",
    values_to = "deaths"
  ) %>%
  mutate(
    deaths = replace_na(deaths, 0),
    Sex = factor(
      case_when(
        Sex == 1L ~ "Male",
        Sex == 2L ~ "Female",
        TRUE ~ NA_character_
      ),
      levels = c("Male", "Female")
    )
  ) %>%
  filter(!is.na(Sex))

# ── Map cause codes -> codeset cause categories ───────────────────────────────
who_long_causes <- who_long %>%
  # Remove 3 character parent ICD-9 codes first
  inner_join(fk_codesets,
             by = c("icd_version" = "icd_version", "Cause" = "icd_code")
  ) %>% 
  mutate(orig_cause = cause_abbrev) 

  # Get totals columns for future work
totals <- who_long %>% 
  filter((Cause == 'A000' & icd_version == '7') |
           (Cause == 'A000' & icd_version == '8') |
           (Cause == 'B00' & icd_version == '9') |
           (Cause == 'AAA' & icd_version == '10'))


# ── Codes with age restriction ────────────────────────────────────────────────
who_long_restrict <- who_long_causes %>%
  mutate(
    reassigned_cause = case_when(
      is.na(orig_cause) ~ "IllDef",
      orig_cause %in% c("SUID", "BirthAsphyx") & Age_Group != "Age_lt1" ~ "MisCode",
      orig_cause == "Perinatal" & is.na(age_specific) & Age_Group != "Age_lt1" ~ "MisCode",
      orig_cause == "Perinatal" & !is.na(age_specific) & Age_Group != "Age_lt1" ~ "NCD-Other",
      orig_cause == "Suicide" & Age_Group == "Age_lt1" ~ "MisCode",
      orig_cause == "Diabetes" & Age_Group == "Age_lt1" ~ "NCD-Other",
      Cause %in% c("R99",
                   "W75",
                   "W750",
                   "W751",
                   "W752",
                   "W754",
                   "W755",
                   "W758",
                   "W759") & icd_version == "10" & Age_Group == "Age_lt1" ~ "SUID",
      TRUE ~ orig_cause
    )
  ) %>% 
  mutate(
    reassigned_cause = case_when(
      Cause == "P351" & icd_version == "10" & Age_Group != "Age_lt1" ~ "InfOther",
      Cause %in% c("P219", "P916", "P209", "G931", "R090",
                   "P90", "P942", "P919") & icd_version == "10" & Age_Group != "Age_lt1" ~ "Neuro",
      TRUE ~ reassigned_cause
    )
  ) %>% 
  dplyr::select(-condition_cluster, -broad_category, -cause_abbrev)

# ── Add broad categories and sum deaths ───────────────────────────────────────
who_intermed <- who_long_restrict %>% 
  left_join(DOMAIN_MAP, by = c('reassigned_cause' = 'cause_abbrev')) %>% 
  rename(Condition_Cluster = cluster, Broad_Category = domain) %>% 
  group_by(Country, Year, icd_version, Age_Group, Sex, Condition_Cluster, Broad_Category) %>% 
  summarise(Deaths = sum(deaths, na.rm = TRUE), .groups = "drop")

# =============================================================================
# IMPUTATION
# =============================================================================

# ── HMD Imputation ───────────────────────────────────────────────────────────
# Only imputing tails for countries without HMD data up until 2023

tail_extrapolate <- function(data,
                             value_col,
                             target_year = 2023) {
  
  value_col <- enquo(value_col)
  
  forecasts <- data %>%
    arrange(Country, Sex, Age_Group, Year) %>%
    group_by(Country, Sex, Age_Group) %>%
    group_modify(~{
      
      last_year <- max(.x$Year)
      
      if (last_year >= target_year) {
        return(tibble())
      }
      
      values <- pull(.x, !!value_col)
      
      avg_change <- mean(diff(tail(values, 5)))
      
      future_years <- seq(last_year + 1, target_year)
      
      future_values <- values[length(values)] +
        avg_change * seq_along(future_years)
      
      tibble(
        Year = future_years,
        !!quo_name(value_col) := future_values
      )
    }) %>%
    mutate(imputed_hmd = 1L) %>% 
    ungroup()
  
  bind_rows(data, forecasts) %>%
    mutate(imputed_hmd = replace_na(imputed_hmd, 0)) %>% 
    arrange(Country, Sex, Age_Group, Year)
}

hmd_pop_imputation <- tail_extrapolate(hmd_intermed, "Population") %>% 
  dplyr::select(-Deaths_AllCause)
hmd_deaths_imputation <- tail_extrapolate(hmd_intermed, "Deaths_AllCause") %>% 
  dplyr::select(-Population)

hmd_imputation_final <- hmd_pop_imputation %>% 
  left_join(hmd_deaths_imputation)

hmd_imputation_final_ages <- hmd_imputation_final %>% 
  mutate(
    Age_Group = case_when(
      Age_Group == "Age_lt1" ~ "Age_lt1",
      Age_Group %in% c("Age_1to9", "Age_10to14", "Age_15to19") ~ "Age_1to19",
      TRUE ~ Age_Group
    )
  ) %>%
  group_by(Country, Year, Age_Group, Sex, imputed_hmd) %>%
  summarise(
    Population = sum(Population, na.rm = TRUE),
    Deaths_AllCause = sum(Deaths_AllCause, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# ANALYTIC DATASET
# =============================================================================

analytic <- who_intermed %>%
  rename(who_code = Country) %>% 
  left_join(COUNTRY_LOOKUP, by = c("who_code")) %>%
  left_join(
    hmd_imputation_final,
    by = c("pop_name" = "Country", "Year", "Sex", "Age_Group")
  ) %>% 
  mutate(Age_Group = factor(Age_Group, levels = AGE_ORDER))
    
analytic_ages <- who_intermed %>%
  mutate(
    Age_Group = case_when(
      Age_Group == "Age_lt1" ~ "Age_lt1",
      Age_Group %in% c("Age_1to9", "Age_10to14", "Age_15to19") ~ "Age_1to19",
      TRUE ~ Age_Group
    )
  ) %>%
  group_by(Country, Year, icd_version, Age_Group, Sex, Condition_Cluster, Broad_Category) %>%
  summarise(
    Deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(who_code = Country) %>% 
  left_join(COUNTRY_LOOKUP, by = c("who_code")) %>%
  left_join(
    hmd_imputation_final_ages,
    by = c("pop_name" = "Country", "Year", "Sex", "Age_Group")
  ) %>% 
  mutate(Age_Group = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19')))

  
analytic_population <- hmd_imputation_final_ages %>% 
  left_join(COUNTRY_LOOKUP, by = c('Country' = 'pop_name'))

pop <- analytic_population %>% 
  filter(!(Country == 'GRB_NP' & Year == 2000)) %>% 
  filter(!(Country == 'AUS' & Year == 2005),
         !(Country == 'BEL' & Year == 2023),
         !(Country == 'DNK' & Year == 2023),
         !(Country == 'NOR' & Year %in% c(2017, 2018, 2019, 2020, 2021, 2022, 2023))) %>% 
  distinct(region, Year, Age_Group, Sex, Population) %>%
  group_by(region, Year, Age_Group) %>%
  summarise(population = sum(Population), .groups = "drop")

overall_pop <- pop %>% 
  group_by(Year, Age_Group) %>%
  summarise(population = sum(population), .groups = "drop")

# =============================================================================
# TABLE 1
# =============================================================================
analytic_population_2023 <- analytic_population %>% filter(Year %in% c(1954, 2023))

table1 <- bind_rows(
  
  # Total population by age group
  analytic_population_2023 %>%
    group_by(Age_Group, Year, region) %>%
    summarise(Population = sum(Population, na.rm = TRUE), .groups = "drop") %>%
    transmute(
      Year = Year,
      Region = region,
      Section = "Age group",
      Category = as.character(Age_Group),
      Population
    ),
  
  # Total population by sex
  analytic_population_2023 %>%
    group_by(Sex, Year, region) %>%
    summarise(Population = sum(Population, na.rm = TRUE), .groups = "drop") %>%
    transmute(
      Year = Year,
      Region = region,
      Section = "Sex",
      Category = as.character(Sex),
      Population
    ),
  
  # Total population by age group and sex
  analytic_population_2023 %>%
    group_by(Age_Group, Sex, Year, region) %>%
    summarise(Population = sum(Population, na.rm = TRUE), .groups = "drop") %>%
    transmute(
      Year = Year,
      Region = region,
      Section = "Age group by sex",
      Category = paste(Age_Group, Sex, sep = ": "),
      Population
    ),
  
  # Total population by country
  analytic_population_2023 %>%
    group_by(country_name, Year, region) %>%
    summarise(Population = sum(Population, na.rm = TRUE), .groups = "drop") %>%
    transmute(
      Year = Year,
      Region = region,
      Section = "Country",
      Category = as.character(country_name),
      Population
    )
) %>%
  group_by(Section, Year, Region) %>%
  mutate(
    Percent = Population / sum(Population),
    `Population, n (%)` = paste0(
      comma(Population),
      " (",
      percent(Percent, accuracy = 0.1),
      ")"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(Year, Region, Section, Category, Population) %>% 
  pivot_wider(names_from = Year, values_from = Population) %>% 
  pivot_wider(names_from = Region, values_from = c(`1954`, `2023`))

write.csv(table1, file.path(PATH_OUTDIR, "tables", "table1.csv"), row.names = FALSE)


# =============================================================================
# DATA QUALITY CHECKS
# =============================================================================

# ── HMD vs. WHO All-Cause Consistency ─────────────────────────────────────────
all_who_deaths <- analytic_ages %>% 
  group_by(country_name, Year, Age_Group, icd_version) %>% 
  summarise(WHO_deaths = sum(Deaths),
            .groups = "drop")

all_hmd_deaths <- analytic_ages %>% 
  distinct(country_name, Year, Age_Group, Deaths_AllCause, icd_version)  %>% 
  group_by(country_name, Year, Age_Group, icd_version) %>% 
  summarise(HMD_deaths = sum(Deaths_AllCause),
            .groups = "drop")

hmd_who_const <- all_who_deaths %>% 
  left_join(all_hmd_deaths) %>% 
  mutate(percent_difference = ((WHO_deaths - HMD_deaths)/HMD_deaths) * 100) %>% 
  arrange(country_name, as.numeric(icd_version))

write.csv(hmd_who_const, file.path(PATH_OUTDIR, "tables", "hmd_vs_who.csv"), row.names = FALSE)

hmd_who_fig <-   hmd_who_const %>% 
  ggplot(aes(x = WHO_deaths, y = HMD_deaths, color = Age_Group)) +
  geom_point(size = 2) +
  facet_wrap(~country_name, ncol = 3, scales = "free") +
  labs(x = "WHO Deaths", 
       y = "HMD Deaths") +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 5
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.background = element_rect(
      fill = "grey90",
    ),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  ) +
  scale_color_manual(name = NULL, values = c(
    "Age_1to19" = "#9E89B9",
    "Age_lt1" = "#79AF97FF"
  )) +
  geom_abline(linetype = "dashed") +labs(x="WHO Total Deaths", y="HMD Total Deaths")+ scale_x_log10() + scale_y_log10()

ggsave(file.path(PATH_OUTDIR, "figures", "HMDvsWHO_deaths.pdf"),
       hmd_who_fig,
       width = 20, height = 20, limitsize = FALSE
)

# ── Proportion of Other deaths (miscoded/ill-defined) ─────────────────────────
death_sum_region <- analytic_ages %>% 
  group_by(Year, region, Age_Group) %>% 
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  )

death_sum_region_age <- analytic_ages %>% 
  group_by(region, Age_Group) %>% 
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  )

all_deaths <- death_sum_region %>% 
  group_by(region) %>% 
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  )

age_deaths <- death_sum_region %>% 
  group_by(Age_Group) %>% 
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  )

country_deaths <- analytic_ages %>% 
  group_by(country_name) %>% 
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  )

icd_deaths <- analytic_ages %>% 
  group_by(icd_version, region) %>% 
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  )

region_age_other_deaths <- analytic_ages %>% 
  filter(Broad_Category == 'Other') %>% 
  group_by(region, Age_Group, Broad_Category) %>% 
  summarise(deaths = sum(Deaths), .groups = "drop") %>% 
  left_join(death_sum_region_age, by = c('region', 'Age_Group')) %>% 
  mutate(deaths = deaths.x, total_deaths = deaths.y) %>% 
  mutate(prop_other_deaths = (deaths/total_deaths)) %>% 
  dplyr::select(-deaths.x, -deaths.y)

region_other_deaths <- analytic_ages %>% 
  filter(Broad_Category == 'Other') %>% 
  group_by(region, Broad_Category) %>% 
  summarise(deaths = sum(Deaths), .groups = "drop") %>% 
  left_join(all_deaths, by = c('region')) %>% 
  mutate(deaths = deaths.x, total_deaths = deaths.y) %>% 
  mutate(prop_other_deaths = (deaths/total_deaths)) %>% 
  dplyr::select(-deaths.x, -deaths.y)

age_other_deaths <- analytic_ages %>% 
  filter(Broad_Category == 'Other') %>% 
  group_by(Age_Group, Broad_Category) %>% 
  summarise(deaths = sum(Deaths), .groups = "drop") %>% 
  left_join(age_deaths, by = c('Age_Group')) %>% 
  mutate(deaths = deaths.x, total_deaths = deaths.y) %>% 
  mutate(prop_other_deaths = (deaths/total_deaths)) %>% 
  dplyr::select(-deaths.x, -deaths.y)

country_other_deaths <- analytic_ages %>% 
  filter(Broad_Category == 'Other') %>% 
  group_by(country_name, Broad_Category) %>% 
  summarise(deaths = sum(Deaths), .groups = "drop") %>% 
  left_join(country_deaths, by = c('country_name')) %>% 
  mutate(deaths = deaths.x, total_deaths = deaths.y) %>% 
  mutate(prop_other_deaths = (deaths/total_deaths)) %>% 
  dplyr::select(-deaths.x, -deaths.y)

icd_other_deaths <- analytic_ages %>% 
  filter(Broad_Category == 'Other') %>% 
  group_by(region, icd_version, Broad_Category) %>% 
  summarise(deaths = sum(Deaths), .groups = "drop") %>% 
  left_join(icd_deaths, by = c('icd_version', 'region')) %>% 
  mutate(deaths = deaths.x, total_deaths = deaths.y) %>% 
  mutate(prop_other_deaths = (deaths/total_deaths)) %>% 
  dplyr::select(-deaths.x, -deaths.y)

write.csv(region_other_deaths, file.path(PATH_OUTDIR, "tables", "region_other_deaths.csv"), row.names = FALSE)
write.csv(age_other_deaths, file.path(PATH_OUTDIR, "tables", "age_other_deaths.csv"), row.names = FALSE)
write.csv(country_other_deaths, file.path(PATH_OUTDIR, "tables", "country_other_deaths.csv"), row.names = FALSE)
write.csv(icd_other_deaths, file.path(PATH_OUTDIR, "tables", "icd_other_deaths.csv"), row.names = FALSE)
write.csv(region_age_other_deaths, file.path(PATH_OUTDIR, "tables", "region_age_other_deaths.csv"), row.names = FALSE)


# =============================================================================
# RATE COMPUTATION
# =============================================================================

# ── Top 10 codes after 1999 ───────────────────────────────────────────────────
fk_codesets_orig_desc <- readxl::read_excel(PATH_CODESETS, sheet = "Codeset", col_types = 'text') %>%
  mutate(
    icd_version  = as.character(`ICD version`),
    icd_code     = as.character(`WHO / ICD Code`),
    condition_cluster = as.character(`Condition Cluster`),
    broad_category = as.character(`Broad Category`),
    injury_cause = `Injury Cause`,
    age_specific = `Logic for only specific codes`
  ) %>%
  dplyr::select(icd_version, icd_code, condition_cluster, broad_category, injury_cause, age_specific, Description) %>% 
  left_join(dplyr::select(DOMAIN_MAP, cause_abbrev, cluster), by = c('condition_cluster' = 'cluster'))

codeset_check_desc <- fk_codesets_orig_desc %>%
  filter(icd_version == 9, str_length(icd_code) > 3) %>%
  mutate(icd_code = substr(icd_code, 1, 3)) %>%
  distinct()

fk_codesets_desc <-
  fk_codesets_orig_desc %>%
  anti_join(codeset_check_desc, by = c('icd_code', 'icd_version'))

icd10_codes <- who_long_causes %>%
  mutate(
    reassigned_cause = case_when(
      is.na(orig_cause) ~ "IllDef",
      orig_cause %in% c("SUID", "BirthAsphyx") & Age_Group != "Age_lt1" ~ "MisCode",
      orig_cause == "Perinatal" & is.na(age_specific) & Age_Group != "Age_lt1" ~ "MisCode",
      orig_cause == "Perinatal" & !is.na(age_specific) & Age_Group != "Age_lt1" ~ "NCD-Other",
      orig_cause == "Suicide" & Age_Group == "Age_lt1" ~ "MisCode",
      orig_cause == "Diabetes" & Age_Group == "Age_lt1" ~ "NCD-Other",
      Cause %in% c("R99",
                   "W75",
                   "W750",
                   "W751",
                   "W752",
                   "W754",
                   "W755",
                   "W758",
                   "W759") & icd_version == "10" & Age_Group == "Age_lt1" ~ "SUID",
      TRUE ~ orig_cause
    )
  ) %>% 
  mutate(
    reassigned_cause = case_when(
      Cause == "P351" & icd_version == "10" & Age_Group != "Age_lt1" ~ "InfOther",
      Cause %in% c("P219", "P916", "P209", "G931", "R090",
                   "P90", "P942", "P919") & icd_version == "10" & Age_Group != "Age_lt1" ~ "Neuro",
      TRUE ~ reassigned_cause
    )
  ) %>% 
  mutate(injury_cause = case_when(
    injury_cause == 'Firearms' & reassigned_cause == "UnintInj" & Age_Group == 'Age_lt1' ~ 'Other',
    injury_cause == 'Poisoning/overdose' & Age_Group == 'Age_lt1' ~ 'Other',
    TRUE ~ injury_cause)) %>% 
  dplyr::select(-condition_cluster, -broad_category, -cause_abbrev) %>% 
  left_join(DOMAIN_MAP, by = c('reassigned_cause' = 'cause_abbrev')) %>% 
  rename(Condition_Cluster = cluster, Broad_Category = domain) %>% 
  mutate(Condition_Cluster = case_when(
    Condition_Cluster != 'Miscoded' & injury_cause == 'Firearms' ~ 'Firearms',
    !is.na(injury_cause) & Condition_Cluster == 'Unintentional Injury' ~ injury_cause,
    TRUE ~ Condition_Cluster
  )) %>% filter(Year >= 1999) %>% 
  mutate(
    Age_Group = case_when(
      Age_Group == "Age_lt1" ~ "Age_lt1",
      Age_Group %in% c("Age_1to9", "Age_10to14", "Age_15to19") ~ "Age_1to19",
      TRUE ~ Age_Group
    )
  ) %>%
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == "Drowning" ~ "Unintentional Injury: Drowning",
    Condition_Cluster == "Motor vehicle crashes" ~ "Unintentional Injury: Motor Vehicle Accidents",
    Condition_Cluster == "Poisoning/overdose" ~ "Unintentional Injury: Poisoning/Overdose",
    Condition_Cluster == "Other" ~ "Unintentional Injury: Other",
    TRUE ~ Condition_Cluster)) %>% 
  group_by(Country, Year, Cause, icd_version, Age_Group, Sex, Condition_Cluster) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  filter(!Condition_Cluster %in% c('Ill-defined and unassigned', 'Miscoded'))

cluster_all_deaths <- icd10_codes %>% 
  group_by(Age_Group, Condition_Cluster) %>% 
  summarise(
    total_deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  )

top10_icd10 <- icd10_codes %>% 
  group_by(Cause, icd_version, Age_Group, Condition_Cluster) %>% 
  summarise(deaths = sum(deaths), .groups = "drop") %>% 
  group_by(Age_Group, Condition_Cluster) %>% 
  slice_max(deaths, n=10) %>% 
  ungroup() %>% 
  left_join(dplyr::select(fk_codesets_desc, icd_code, icd_version, Description), 
            by = c('Cause' = 'icd_code', 'icd_version')) %>% 
  left_join(cluster_all_deaths, by = c('Condition_Cluster', 'Age_Group')) %>% 
  mutate(rate = deaths/total_deaths) %>% 
  arrange(Condition_Cluster, Age_Group, desc(deaths), rate)

write.csv(top10_icd10, file.path(PATH_OUTDIR, "tables", "top10_icd10.csv"), row.names = FALSE)

# ── Overall mortality ─────────────────────────────────────────────────────────
overall_rates <- analytic_ages %>% 
  dplyr::select(-Population) %>% 
  group_by(Year, Sex, Age_Group, Condition_Cluster, Broad_Category) %>%
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  )

overall_domain_rates <- overall_rates %>%
  group_by(Year, Age_Group, Broad_Category) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(overall_pop, by = c("Year", "Age_Group")) %>% 
  filter(population > 0) %>%
  mutate(
    rate_per_100k = deaths / population * 1e5,
    Broad_Category = factor(Broad_Category, levels = DOMAIN_ORDER),
    Age_Group     = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19'))
  )


# ── Rates by sex ──────────────────────────────────────────────────────────────
rates_by_sex <- analytic_ages %>%
  dplyr::select(-Population) %>% 
  group_by(region, Year, Sex, Age_Group, Condition_Cluster, Broad_Category) %>%
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  )

# ── Rates by cluster ──────────────────────────────────────────────────────────
cluster_rates <- rates_by_sex %>%
  group_by(region, Year, Age_Group, Condition_Cluster, Broad_Category) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pop, by = c("region", "Year", "Age_Group")) %>% 
  filter(population > 0) %>%
  mutate(
    rate_per_100k = deaths / population * 100000,
    Age_Group     = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19'))
  )

# ── Rates by domain ───────────────────────────────────────────────────────────
domain_rates <- rates_by_sex %>%
  group_by(region, Year, Age_Group, Broad_Category) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pop, by = c("region", "Year", "Age_Group")) %>% 
  filter(population > 0) %>%
  mutate(
    rate_per_100k = deaths / population * 1e5,
    Broad_Category = factor(Broad_Category, levels = DOMAIN_ORDER),
    Age_Group     = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19'))
  )

# ── Domain Rate Plots ─────────────────────────────────────────────────────────
# Domain crossover years
fig1_crossover <- domain_rates %>% filter(Broad_Category != "Other") %>% 
  dplyr::select(-deaths, -population) %>% 
  pivot_wider(
    names_from = region,
    values_from = rate_per_100k
  ) %>%
  mutate(
    us_higher = !is.na(US) & !is.na(OECD13) & US > OECD13
  ) %>%
  arrange(Age_Group, Broad_Category, Year) %>%
  group_by(Age_Group, Broad_Category) %>%
  mutate(
    new_sequence = (
      us_higher != lag(us_higher, default = first(us_higher))
    ) | (
      Year != lag(Year, default = first(Year)) + 1
    ),
    sequence_id = cumsum(replace_na(new_sequence, TRUE))
  ) %>%
  group_by(Age_Group, Broad_Category, sequence_id) %>%
  mutate(
    consecutive_years = n()
  ) %>%
  ungroup() %>%
  filter(us_higher, consecutive_years >= 3) %>%
  arrange(Age_Group, Broad_Category, Year) %>% 
  filter(us_higher, new_sequence) %>%
  filter(Year != 1954) %>% 
  mutate(crossover_flag = 1L) %>% 
  dplyr::select(Year, Age_Group, Broad_Category, crossover_flag)

fig1_table <- domain_rates %>% filter(Broad_Category != "Other") %>% 
  pivot_wider(
    names_from = region,
    values_from = c(deaths, population, rate_per_100k)
  ) %>% 
  mutate(abs_rate_diff = abs(rate_per_100k_US - rate_per_100k_OECD13),
         rate_ratio = rate_per_100k_US/rate_per_100k_OECD13) %>% 
  left_join(fig1_crossover, by = c('Year', 'Age_Group', 'Broad_Category')) %>% 
  mutate(crossover_flag = if_else(is.na(crossover_flag), 0L, crossover_flag))
  
write.csv(fig1_table, file.path(PATH_OUTDIR, "tables", "fig1_table.csv"), row.names = FALSE)

# Domain plot 1954-2023
fig1 <- domain_rates %>% filter(Broad_Category != "Other") %>% 
  ggplot(
    aes(
      x = Year,
      y = rate_per_100k,
      color = region,
      group = region
    )
  ) +
  geom_line(linewidth = 0.9) +
  icd_vlines() +
  scale_color_manual(
    name = NULL,
    values = c(
      "US" = "#00A1D5FF",
      "OECD13" = "#DF8F44FF"
    ),
    breaks = c("US", "OECD13")
  ) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_grid2(
    rows = vars(Broad_Category),
    cols = vars(Age_Group),
    independent = "y",
    scales = "free_y",
    labeller = labeller(
      Broad_Category = label_wrap_gen(width = 20),
      Age_Group = c(
        Age_lt1 = "A. <1-Year-Old Infants",
        Age_1to19 = "B. 1- to 19-Year-Old Individuals"
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Deaths per 100,000"
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 10
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_text(
      size = 14,
      hjust = 0,
      margin = margin(t = 5, r = 5, b = 5, l = 5)
    ),
    strip.text.y = element_text(size = 14),
    
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "fig1_domain_rates.pdf"),
  plot = fig1,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

# Domain plot 1999-2023
fig1b <- domain_rates %>%
  filter(Broad_Category != "Other",
         Year >= 1999) %>%
  ggplot(
    aes(
      x = Year,
      y = rate_per_100k,
      color = region,
      group = region
    )
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    name = NULL,
    values = c(
      "US" = "#00A1D5FF",
      "OECD13" = "#DF8F44FF"
    ),
    breaks = c("US", "OECD13")
  ) +
  scale_x_continuous(breaks = seq(2000, 2030, 5)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_grid2(
    rows = vars(Broad_Category),
    cols = vars(Age_Group),
    independent = "y",
    scales = "free_y",
    labeller = labeller(
      Broad_Category = label_wrap_gen(width = 20),
      Age_Group = c(
        Age_lt1 = "A. <1-Year-Old Infants",
        Age_1to19 = "B. 1- to 19-Year-Old Individuals"
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Deaths per 100,000"
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 10
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_text(
      size = 14,
      hjust = 0,
      margin = margin(t = 5, r = 5, b = 5, l = 5)
    ),
    strip.text.y = element_text(size = 14),
    
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "suppfig1_domain_rates.pdf"),
  plot = fig1b,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ── Infection Rate Plots ──────────────────────────────────────────────────────
supp_infection_table <- cluster_rates %>%
  pivot_wider(
    names_from = region,
    values_from = c(deaths, population, rate_per_100k)
  ) %>% 
  mutate(abs_rate_diff = abs(rate_per_100k_US - rate_per_100k_OECD13),
         rate_ratio = rate_per_100k_US/rate_per_100k_OECD13) %>% 
  filter(Broad_Category == "Infections")

write.csv(supp_infection_table, file.path(PATH_OUTDIR, "tables", "supp_infection_table.csv"), row.names = FALSE)

# Infections plot 1954-2023
supplfig5 <- cluster_rates %>%
  filter(Broad_Category == "Infections") %>%
  dplyr::select(-Broad_Category) %>% 
  mutate(Condition_Cluster = if_else(Condition_Cluster == "Diarrheal/gastrointestinal infections", 
                                     "Diarrheal/Gastrointestinal Infections", Condition_Cluster)
  ) %>% 
  ggplot(
    aes(
      x = Year,
      y = rate_per_100k,
      color = region,
      group = region
    )
  ) +
  geom_line(linewidth = 0.9) +
  icd_vlines() +
  scale_color_manual(
    name = NULL,
    values = c(
      "US" = "#00A1D5FF",
      "OECD13" = "#DF8F44FF"
    ),
    breaks = c("US", "OECD13")
  ) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_grid2(
    rows = vars(Condition_Cluster),
    cols = vars(Age_Group),
    independent = "y",
    scales = "free",
    render_empty = FALSE,
    labeller = ggplot2::labeller(
      Condition_Cluster = label_wrap_gen(width = 22),
      Age_Group = c(
        Age_lt1 = "A. <1-Year-Old Infants",
        Age_1to19 = "B. 1- to 19-Year-Old Individuals"
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Deaths per 100,000"
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 10
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.text.x = element_text(
      size = 14,
      hjust = 0,
      margin = margin(t = 5, r = 5, b = 5, l = 5)
    ),
    strip.text.y = element_text(size = 14),
    
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "supp_infections_rates.pdf"),
  plot = supplfig5,
  width = 12,
  height = 15,
  units = "in",
  dpi = 600,
  bg = "white"
)

# Infections plot 1999-2023
supplfig5_icd10 <- cluster_rates %>%
  filter(Broad_Category == "Infections",
         Year >= 1999) %>%
  dplyr::select(-Broad_Category) %>% 
  mutate(Condition_Cluster = if_else(Condition_Cluster == "Diarrheal/gastrointestinal infections", 
                                     "Diarrheal/Gastrointestinal Infections", Condition_Cluster)
  ) %>% 
  ggplot(
    aes(
      x = Year,
      y = rate_per_100k,
      color = region,
      group = region
    )
  ) +
  geom_line(linewidth = 0.9) +
  icd_vlines() +
  scale_color_manual(
    name = NULL,
    values = c(
      "US" = "#00A1D5FF",
      "OECD13" = "#DF8F44FF"
    ),
    breaks = c("US", "OECD13")
  ) +
  scale_x_continuous(breaks = seq(2000, 2030, 5)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_grid2(
    rows = vars(Condition_Cluster),
    cols = vars(Age_Group),
    independent = "y",
    scales = "free",
    render_empty = FALSE,
    labeller = ggplot2::labeller(
      Condition_Cluster = label_wrap_gen(width = 22),
      Age_Group = c(
        Age_lt1 = "A. <1-Year-Old Infants",
        Age_1to19 = "B. 1- to 19-Year-Old Individuals"
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Deaths per 100,000"
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 10
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.text.x = element_text(
      size = 14,
      hjust = 0,
      margin = margin(t = 5, r = 5, b = 5, l = 5)
    ),
    strip.text.y = element_text(size = 14),
    
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "supp_infections_rates_icd10.pdf"),
  plot = supplfig5,
  width = 12,
  height = 15,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ── Injury Rate Plots ─────────────────────────────────────────────────────────
# Rerun rates to split out Unintentional Injury into subgroups
# Make cluster-specific age restrictions
who_intermed_inury <- who_long_restrict %>%
  mutate(injury_cause = case_when(
    injury_cause == 'Firearms' & reassigned_cause == "UnintInj" & Age_Group == 'Age_lt1' ~ 'Other',
    injury_cause == 'Poisoning/overdose' & Age_Group == 'Age_lt1' ~ 'Other',
    TRUE ~ injury_cause)) %>% 
  left_join(DOMAIN_MAP, by = c('reassigned_cause' = 'cause_abbrev')) %>%
  rename(Condition_Cluster = cluster, Broad_Category = domain) %>%
  group_by(Country, Year, icd_version, Age_Group, Sex, Condition_Cluster, Broad_Category, injury_cause) %>%
  summarise(Deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
  filter(Condition_Cluster == "Unintentional Injury", injury_cause != "See coding logic")

who_intermed_injury_age <- who_intermed_inury %>% 
  mutate(
    Age_Group = case_when(
      Age_Group == "Age_lt1" ~ "Age_lt1",
      Age_Group %in% c("Age_1to9", "Age_10to14", "Age_15to19") ~ "Age_1to19",
      TRUE ~ Age_Group)
  ) %>%
  left_join(COUNTRY_LOOKUP, by = c('Country' = 'who_code')) %>% 
  group_by(Country, Year, icd_version, Age_Group, Sex, Condition_Cluster, Broad_Category, injury_cause) %>% 
  summarise(Deaths = sum(Deaths), .groups = "drop")

analytic_injury_age <- who_intermed_injury_age %>%
  rename(who_code = Country) %>% 
  left_join(COUNTRY_LOOKUP, by = c("who_code")) %>%
  left_join(
    hmd_imputation_final_ages,
    by = c("pop_name" = "Country", "Year", "Sex", "Age_Group")
  ) %>% 
  mutate(Age_Group = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19')))

injury_clusters_ages <- analytic_injury_age %>% 
  group_by(region, Year, Condition_Cluster, injury_cause, Age_Group) %>%
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pop, by = c("region", "Year", "Age_Group")) %>% 
  mutate(rate_per_100k = deaths / population * 100000) %>% 
  dplyr::select(-Condition_Cluster) %>% 
  rename(Condition_Cluster = injury_cause) %>% 
  mutate(Age_Group = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19')))

fig2_table <- cluster_rates %>%
  filter(Broad_Category == "Injury") %>%
  dplyr::select(-Broad_Category) %>% 
  dplyr::union(injury_clusters_ages) %>% 
  filter(Condition_Cluster != "Unintentional Injury") %>% 
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == "Drowning" ~ "Unintentional Injury: Drowning",
    Condition_Cluster == "Firearms" ~ "Unintentional Injury: Firearms",
    Condition_Cluster == "Motor vehicle crashes" ~ "Unintentional Injury: Motor Vehicle Accidents",
    Condition_Cluster == "Poisoning/overdose" ~ "Unintentional Injury: Poisoning/Overdose",
    Condition_Cluster == "Other" ~ "Unintentional Injury: Other",
    TRUE ~ Condition_Cluster
  )) %>% 
  pivot_wider(
    names_from = region,
    values_from = c(deaths, population, rate_per_100k)
  ) %>% 
  mutate(abs_rate_diff = abs(rate_per_100k_US - rate_per_100k_OECD13),
         rate_ratio = rate_per_100k_US/rate_per_100k_OECD13)

write.csv(fig2_table, file.path(PATH_OUTDIR, "tables", "fig2_table.csv"), row.names = FALSE)

final_all_cluster_rates <- cluster_rates %>% 
  dplyr::select(-Broad_Category) %>% 
  dplyr::union(injury_clusters_ages) %>% 
  filter(Condition_Cluster != "Unintentional Injury") %>% 
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == "Drowning" ~ "Unintentional Injury: Drowning",
    Condition_Cluster == "Firearms" ~ "Unintentional Injury: Firearms",
    Condition_Cluster == "Motor vehicle crashes" ~ "Unintentional Injury: Motor Vehicle Accidents",
    Condition_Cluster == "Poisoning/overdose" ~ "Unintentional Injury: Poisoning/Overdose",
    Condition_Cluster == "Other" ~ "Unintentional Injury: Other",
    TRUE ~ Condition_Cluster
  ))
  

injury_order <- c("Homicide", "Suicide",
                  "Unintentional Injury: Drowning",
                  "Unintentional Injury: Firearms",
                  "Unintentional Injury: Motor Vehicle Accidents", "Unintentional Injury: Poisoning/Overdose",
                  "Unintentional Injury: Other")

# Injury plot 1954-2023
fig2 <- cluster_rates %>%
  filter(Broad_Category == "Injury") %>%
  dplyr::select(-Broad_Category) %>% 
  dplyr::union(injury_clusters_ages) %>% 
  filter(Condition_Cluster != "Unintentional Injury") %>% 
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == "Drowning" ~ "Unintentional Injury: Drowning",
    Condition_Cluster == "Firearms" ~ "Unintentional Injury: Firearms",
    Condition_Cluster == "Motor vehicle crashes" ~ "Unintentional Injury: Motor Vehicle Accidents",
    Condition_Cluster == "Poisoning/overdose" ~ "Unintentional Injury: Poisoning/Overdose",
    Condition_Cluster == "Other" ~ "Unintentional Injury: Other",
    TRUE ~ Condition_Cluster
  )) %>% 
  mutate(Condition_Cluster = factor(Condition_Cluster, levels = injury_order)) %>% 
  ggplot(
    aes(
      x = Year,
      y = rate_per_100k,
      color = region,
      group = region
    )
  ) +
  geom_line(linewidth = 0.9) +
  icd_vlines() +
  scale_color_manual(
    name = NULL,
    values = c(
      "US" = "#00A1D5FF",
      "OECD13" = "#DF8F44FF"
    ),
    breaks = c("US", "OECD13")
  ) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_grid2(
    rows = vars(Condition_Cluster),
    cols = vars(Age_Group),
    independent = "y",
    scales = "free",
    render_empty = FALSE,
    labeller = ggplot2::labeller(
      Condition_Cluster = label_wrap_gen(width = 22),
      Age_Group = c(
        Age_lt1 = "A. <1-Year-Old Infants",
        Age_1to19 = "B. 1- to 19-Year-Old Individuals"
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Deaths per 100,000"
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 10
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.text.x = element_text(
      size = 14,
      hjust = 0,
      margin = margin(t = 5, r = 5, b = 5, l = 5)
    ),
    strip.text.y = element_text(size = 14),
    
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "fig2_injury_rates.pdf"),
  plot = fig2,
  width = 12,
  height = 17,
  units = "in",
  dpi = 600,
  bg = "white"
)

# For 1999-2023, break out firearms from homicide and suicide
# Labelled decomp because will also use for decomposition analysis
who_long_restrict_decomp <- who_long_causes %>%
  mutate(
    reassigned_cause = case_when(
      is.na(orig_cause) ~ "IllDef",
      orig_cause %in% c("SUID", "BirthAsphyx") & Age_Group != "Age_lt1" ~ "MisCode",
      orig_cause == "Perinatal" & is.na(age_specific) & Age_Group != "Age_lt1" ~ "MisCode",
      orig_cause == "Perinatal" & !is.na(age_specific) & Age_Group != "Age_lt1" ~ "NCD-Other",
      orig_cause == "Suicide" & Age_Group == "Age_lt1" ~ "MisCode",
      orig_cause == "Diabetes" & Age_Group == "Age_lt1" ~ "NCD-Other",
      Cause %in% c("R99",
                   "W75",
                   "W750",
                   "W751",
                   "W752",
                   "W754",
                   "W755",
                   "W758",
                   "W759") & icd_version == "10" & Age_Group == "Age_lt1" ~ "SUID",
      TRUE ~ orig_cause
    )
  ) %>% 
  mutate(
    reassigned_cause = case_when(
      Cause == "P351" & icd_version == "10" & Age_Group != "Age_lt1" ~ "InfOther",
      Cause %in% c("P219", "P916", "P209", "G931", "R090",
                   "P90", "P942", "P919") & icd_version == "10" & Age_Group != "Age_lt1" ~ "Neuro",
      TRUE ~ reassigned_cause
    )
  ) %>% 
  mutate(injury_cause = case_when(
    injury_cause == 'Firearms' & reassigned_cause == "UnintInj" & Age_Group == 'Age_lt1' ~ 'Other',
    injury_cause == 'Poisoning/overdose' & Age_Group == 'Age_lt1' ~ 'Other',
    TRUE ~ injury_cause)) %>% 
  dplyr::select(-condition_cluster, -broad_category, -cause_abbrev)

who_intermed_decomp <- who_long_restrict_decomp %>% 
  left_join(DOMAIN_MAP, by = c('reassigned_cause' = 'cause_abbrev')) %>% 
  rename(Condition_Cluster = cluster, Broad_Category = domain) %>% 
  mutate(Condition_Cluster = case_when(
    Condition_Cluster != 'Miscoded' & injury_cause == 'Firearms' ~ 'Firearms',
    !is.na(injury_cause) & Condition_Cluster == 'Unintentional Injury' ~ injury_cause,
    TRUE ~ Condition_Cluster
  )) %>% 
  group_by(Country, Year, icd_version, Age_Group, Sex, Condition_Cluster, Broad_Category) %>% 
  summarise(Deaths = sum(deaths, na.rm = TRUE), .groups = "drop")

analytic_ages_decomp <- who_intermed_decomp %>%
  mutate(
    Age_Group = case_when(
      Age_Group == "Age_lt1" ~ "Age_lt1",
      Age_Group %in% c("Age_1to9", "Age_10to14", "Age_15to19") ~ "Age_1to19",
      TRUE ~ Age_Group
    )
  ) %>%
  group_by(Country, Year, icd_version, Age_Group, Sex, Condition_Cluster, Broad_Category) %>%
  summarise(
    Deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(who_code = Country) %>% 
  left_join(COUNTRY_LOOKUP, by = c("who_code")) %>%
  left_join(
    hmd_imputation_final_ages,
    by = c("pop_name" = "Country", "Year", "Sex", "Age_Group")
  ) %>% 
  mutate(Age_Group = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19')))

rates_by_sex_decomp <- analytic_ages_decomp %>%
  dplyr::select(-Population) %>% 
  group_by(region, Year, Sex, Age_Group, Condition_Cluster, Broad_Category) %>%
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  )

cluster_rates_decomp <- rates_by_sex_decomp %>%
  group_by(region, Year, Age_Group, Condition_Cluster, Broad_Category) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pop, by = c("region", "Year", "Age_Group")) %>% 
  filter(population > 0) %>%
  mutate(
    rate_per_100k = deaths / population * 100000,
    Age_Group     = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19'))
  )

injury_orderb <- c("Firearms", "Homicide", "Suicide",
                  "Unintentional Injury: Drowning",
                  "Unintentional Injury: Motor Vehicle Accidents", "Unintentional Injury: Poisoning/Overdose",
                  "Unintentional Injury: Other")

suppfig2_table <- cluster_rates_decomp %>%
  filter(Broad_Category == "Injury",
         Year >= 1999) %>%
  dplyr::select(-Broad_Category) %>% 
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == "Drowning" ~ "Unintentional Injury: Drowning",
    Condition_Cluster == "Motor vehicle crashes" ~ "Unintentional Injury: Motor Vehicle Accidents",
    Condition_Cluster == "Poisoning/overdose" ~ "Unintentional Injury: Poisoning/Overdose",
    Condition_Cluster == "Other" ~ "Unintentional Injury: Other",
    TRUE ~ Condition_Cluster
  )) %>% pivot_wider(
    names_from = region,
    values_from = c(deaths, population, rate_per_100k)
  ) %>% 
  mutate(abs_rate_diff = abs(rate_per_100k_US - rate_per_100k_OECD13),
         rate_ratio = rate_per_100k_US/rate_per_100k_OECD13)

write.csv(suppfig2_table, file.path(PATH_OUTDIR, "tables", "supp_injury_icd10.csv"), row.names = FALSE)

# Injury plot 1999-2023
fig2b <- cluster_rates_decomp %>%
  filter(Broad_Category == "Injury",
         Year >= 1999) %>%
  dplyr::select(-Broad_Category) %>% 
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == "Drowning" ~ "Unintentional Injury: Drowning",
    Condition_Cluster == "Motor vehicle crashes" ~ "Unintentional Injury: Motor Vehicle Accidents",
    Condition_Cluster == "Poisoning/overdose" ~ "Unintentional Injury: Poisoning/Overdose",
    Condition_Cluster == "Other" ~ "Unintentional Injury: Other",
    TRUE ~ Condition_Cluster
  )) %>% 
  mutate(Condition_Cluster = factor(Condition_Cluster, levels = injury_orderb)) %>% 
  ggplot(
    aes(
      x = Year,
      y = rate_per_100k,
      color = region,
      group = region
    )
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    name = NULL,
    values = c(
      "US" = "#00A1D5FF",
      "OECD13" = "#DF8F44FF"
    ),
    breaks = c("US", "OECD13")
  ) +
  scale_x_continuous(breaks = seq(1999, 2030, 5)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_grid2(
    rows = vars(Condition_Cluster),
    cols = vars(Age_Group),
    independent = "y",
    scales = "free",
    render_empty = FALSE,
    labeller = ggplot2::labeller(
      Condition_Cluster = label_wrap_gen(width = 22),
      Age_Group = c(
        Age_lt1 = "A. <1-Year-Old Infants",
        Age_1to19 = "B. 1- to 19-Year-Old Individuals"
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Deaths per 100,000"
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 10
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_text(
      size = 14,
      hjust = 0,
      margin = margin(t = 5, r = 5, b = 5, l = 5)
    ),
    strip.text.y = element_text(size = 14),
    
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "suppfig2_injury_rates_icd10.pdf"),
  plot = fig2b,
  width = 12,
  height = 17,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ── NCD Rate Plots ────────────────────────────────────────────────────────────
fig3_table <- cluster_rates %>%
  filter(Broad_Category == "Non-Communicable Diseases") %>% 
  dplyr::select(-Broad_Category) %>% 
  pivot_wider(
    names_from = region,
    values_from = c(deaths, population, rate_per_100k)
  ) %>% 
  mutate(abs_rate_diff = abs(rate_per_100k_US - rate_per_100k_OECD13),
         rate_ratio = rate_per_100k_US/rate_per_100k_OECD13)

write.csv(fig3_table, file.path(PATH_OUTDIR, "tables", "fig3_table.csv"), row.names = FALSE)

# NCD plot <1 1954-2023
fig3a <- cluster_rates %>%
  filter(Broad_Category == "Non-Communicable Diseases",
         Age_Group == "Age_lt1") %>%
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == 'Cardiovascular disorders' ~ 'Cardiovascular Disorders',
    Condition_Cluster == 'Neuro-psychiatric disorders' ~ 'Neuro-Psychiatric Disorders',
    Condition_Cluster == 'Gastrointestinal/Hepatic disorders, Non-Infectious' ~ 'Gastrointestinal/Hepatic Disorders',
    Condition_Cluster == 'Genitourinary disorders, Non-Infectious' ~ 'Genitourinary Disorders',
    Condition_Cluster == 'Prematurity/neonatal disorders' ~ 'Prematurity/Neonatal Disorders',
    Condition_Cluster == 'Birth asphyxia/birth trauma' ~ 'Birth Asphyxia/Birth Trauma',
    TRUE ~ Condition_Cluster
  )) %>% 
  ggplot(
    aes(
      x = Year,
      y = rate_per_100k,
      color = region,
      group = region
    )
  ) +
  geom_line(linewidth = 0.9) +
  icd_vlines() +
  scale_color_manual(
    name = NULL,
    values = c(
      "US" = "#00A1D5FF",
      "OECD13" = "#DF8F44FF"
    ),
    breaks = c("US", "OECD13")
  ) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_wrap(
    ~ Condition_Cluster,
    ncol = 3,
    scales = "free_y",
    labeller = labeller(
      Age_Group = c(
        Age_lt1 = "A. <1-Year-Old Infants",
        Age_1to19 = "B. 1- to 19-Year-Old Individuals"
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Deaths per 100,000"
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 10
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.background = element_rect(
      fill = "grey90",
    ),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "fig3_ncd_rates.pdf"),
  plot = fig3a,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

# NCD plot 1-19 1954-2023
fig3b <- cluster_rates %>%
  filter(Broad_Category == "Non-Communicable Diseases",
         Age_Group == "Age_1to19") %>%
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == 'Cardiovascular disorders' ~ 'Cardiovascular Disorders',
    Condition_Cluster == 'Neuro-psychiatric disorders' ~ 'Neuro-Psychiatric Disorders',
    Condition_Cluster == 'Gastrointestinal/Hepatic disorders, Non-Infectious' ~ 'Gastrointestinal/Hepatic Disorders',
    Condition_Cluster == 'Genitourinary disorders, Non-Infectious' ~ 'Genitourinary Disorders',
    Condition_Cluster == 'Prematurity/neonatal disorders' ~ 'Prematurity/Neonatal Disorders',
    Condition_Cluster == 'Birth asphyxia/birth trauma' ~ 'Birth Asphyxia/Birth Trauma',
    TRUE ~ Condition_Cluster
  )) %>% 
  ggplot(
    aes(
      x = Year,
      y = rate_per_100k,
      color = region,
      group = region
    )
  ) +
  geom_line(linewidth = 0.9) +
  icd_vlines() +
  scale_color_manual(
    name = NULL,
    values = c(
      "US" = "#00A1D5FF",
      "OECD13" = "#DF8F44FF"
    ),
    breaks = c("US", "OECD13")
  ) +
  scale_x_continuous(breaks = seq(1950, 2030, 10)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_wrap(
    ~ Condition_Cluster,
    ncol = 3,
    scales = "free_y",
    labeller = labeller(
      Age_Group = c(
        Age_lt1 = "A. <1-Year-Old Infants",
        Age_1to19 = "B. 1- to 19-Year-Old Individuals"
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Deaths per 100,000"
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 10
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.background = element_rect(
      fill = "grey90",
    ),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "suppfig3b_ncd_rates_1-19.pdf"),
  plot = fig3b,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

# NCD plot <1 1999-2023
fig3a_icd10 <- cluster_rates %>%
  filter(Broad_Category == "Non-Communicable Diseases",
         Age_Group == "Age_lt1",
         Year >= 1999) %>%
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == 'Cardiovascular disorders' ~ 'Cardiovascular Disorders',
    Condition_Cluster == 'Neuro-psychiatric disorders' ~ 'Neuro-Psychiatric Disorders',
    Condition_Cluster == 'Gastrointestinal/Hepatic disorders, Non-Infectious' ~ 'Gastrointestinal/Hepatic Disorders',
    Condition_Cluster == 'Genitourinary disorders, Non-Infectious' ~ 'Genitourinary Disorders',
    Condition_Cluster == 'Prematurity/neonatal disorders' ~ 'Prematurity/Neonatal Disorders',
    Condition_Cluster == 'Birth asphyxia/birth trauma' ~ 'Birth Asphyxia/Birth Trauma',
    Condition_Cluster == 'Respiratory disorders, Non-Infectious' ~ 'Respiratory Disorders',
    TRUE ~ Condition_Cluster
  )) %>% 
  ggplot(
    aes(
      x = Year,
      y = rate_per_100k,
      color = region,
      group = region
    )
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    name = NULL,
    values = c(
      "US" = "#00A1D5FF",
      "OECD13" = "#DF8F44FF"
    ),
    breaks = c("US", "OECD13")
  ) +
  scale_x_continuous(breaks = seq(2000, 2030, 5)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_wrap(
    ~ Condition_Cluster,
    ncol = 3,
    scales = "free_y",
    labeller = labeller(
      Age_Group = c(
        Age_lt1 = "A. <1-Year-Old Infants",
        Age_1to19 = "B. 1- to 19-Year-Old Individuals"
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Deaths per 100,000"
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 10
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.background = element_rect(
      fill = "grey90",
    ),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "suppfig4a_ncd_rates_icd10.pdf"),
  plot = fig3a_icd10,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

# NCD plot 1-19 1999-2023
fig3b_icd10 <- cluster_rates %>%
  filter(Broad_Category == "Non-Communicable Diseases",
         Age_Group == "Age_1to19",
         Year >= 1999) %>%
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == 'Cardiovascular disorders' ~ 'Cardiovascular Disorders',
    Condition_Cluster == 'Neuro-psychiatric disorders' ~ 'Neuro-Psychiatric Disorders',
    Condition_Cluster == 'Gastrointestinal/Hepatic disorders, Non-Infectious' ~ 'Gastrointestinal/Hepatic Disorders',
    Condition_Cluster == 'Genitourinary disorders, Non-Infectious' ~ 'Genitourinary Disorders',
    Condition_Cluster == 'Prematurity/neonatal disorders' ~ 'Prematurity/Neonatal Disorders',
    Condition_Cluster == 'Birth asphyxia/birth trauma' ~ 'Birth Asphyxia/Birth Trauma',
    Condition_Cluster == 'Respiratory disorders, Non-Infectious' ~ 'Respiratory Disorders',
    TRUE ~ Condition_Cluster
  )) %>% 
  ggplot(
    aes(
      x = Year,
      y = rate_per_100k,
      color = region,
      group = region
    )
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    name = NULL,
    values = c(
      "US" = "#00A1D5FF",
      "OECD13" = "#DF8F44FF"
    ),
    breaks = c("US", "OECD13")
  ) +
  scale_x_continuous(breaks = seq(2000, 2030, 5)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_wrap(
    ~ Condition_Cluster,
    ncol = 3,
    scales = "free_y",
    labeller = labeller(
      Age_Group = c(
        Age_lt1 = "A. <1-Year-Old Infants",
        Age_1to19 = "B. 1- to 19-Year-Old Individuals"
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Deaths per 100,000"
  ) +
  guides(
    color = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.7,
        size = 10
      )
    )
  ) +
  theme_mortality() +
  theme(
    strip.background = element_rect(
      fill = "grey90",
    ),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "suppfig4b_ncd_rates_icd10.pdf"),
  plot = fig3b_icd10,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ── Country Specific Mortality Rates ──────────────────────────────────────────
# Domain specific country rates
pop_country <- analytic_population %>% 
  filter(!(Country == 'GRB_NP' & Year == 2000)) %>% 
  filter(!(Country == 'AUS' & Year == 2005),
         !(Country == 'BEL' & Year == 2023),
         !(Country == 'DNK' & Year == 2023),
         !(Country == 'NOR' & Year %in% c(2017, 2018, 2019, 2020, 2021, 2022, 2023))) %>% 
  distinct(country_name, Year, Age_Group, Sex, Population) %>%
  group_by(country_name, Year, Age_Group) %>%
  summarise(population = sum(Population), .groups = "drop")

rates_by_sex_country <- analytic_ages %>%
  dplyr::select(-Population) %>% 
  group_by(country_name, Year, Sex, Age_Group, Condition_Cluster, Broad_Category) %>%
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  )

domain_rates_country <- rates_by_sex_country %>%
  group_by(country_name, Year, Age_Group, Broad_Category) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pop_country, by = c("country_name", "Year", "Age_Group")) %>% 
  filter(population > 0) %>%
  mutate(
    rate_per_100k = deaths / population * 1e5,
    Broad_Category = factor(Broad_Category, levels = DOMAIN_ORDER),
    Age_Group     = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19'))
  )

country_death_rates <- domain_rates_country %>% 
  filter(Broad_Category != 'Other')
write.csv(country_death_rates, file.path(PATH_OUTDIR, "tables", "country_death_rates.csv"), row.names = FALSE)

make_slug <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("&", "and") %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_remove_all("^_|_$")
}

age_filename <- function(x) {
  case_when(
    x == "Age_lt1"   ~ "lt1",
    x == "Age_1to19" ~ "1to19",
    TRUE             ~ make_slug(x)
  )
}

missing_periods <- tribble(
  ~country_name,    ~start_year, ~end_year,
  "Australia",       2005,        2005,
  "Belgium",         2023,        2023,
  "Denmark",         2023,        2023,
  "United Kingdom",  2000,        2000,
  "Norway",          2017,        2023
)

# One function for all domain country-specific rate plots
make_country_plot <- function(broad_category, age_group) {
  
  country_plot <- country_death_rates %>%
    filter(
      Broad_Category == broad_category,
      Age_Group == age_group
    ) %>%
    ggplot(
      aes(
        x = Year,
        y = rate_per_100k,
        color = Age_Group,
        group = Age_Group
      )
    ) +
    geom_rect(
      data = missing_periods,
      aes(
        xmin = start_year - 0.5,
        xmax = end_year + 0.5,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "#F4C95D",
      alpha = 0.25
    ) +
    geom_line(linewidth = 0.9) +
    icd_vlines() +
    scale_color_manual(
      name = NULL,
      labels = c(
        "Age_lt1" = "<1-Year-Old Infants",
        "Age_1to19" = "1- to 19-Year-Old Individuals"
      ),
      values = c(
        "Age_lt1" = "#79AF97FF",
        "Age_1to19" = "#9E89B9"
      ),
      breaks = c("Age_lt1", "Age_1to19")
    ) +
    scale_x_continuous(
      breaks = seq(1950, 2030, 10)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    facet_wrap(
      ~ country_name,
      ncol = 3,
      scales = "free_y"
    ) +
    labs(
      x = "Year",
      y = "Deaths per 100,000",
      caption = paste0(
        "Note. Shaded areas indicate years with unavailable WHO-MD data: ",
        "Australia (2005), Belgium (2023), Denmark (2023), ",
        "United Kingdom (2000), and Norway (2017-2023)."
      )
    ) +
    guides(
      color = guide_legend(
        title.position = "left",
        title.hjust = 0.5,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = 1.7
        )
      )
    ) +
    theme_mortality() +
    theme(
      strip.background = element_rect(
        fill = "grey90"
      ),
      strip.text = element_text(
        face = "bold",
        size = 14
      ),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.text = element_text(size = 12),
      panel.grid = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      
      plot.caption = element_text(
        hjust = 0,
        size = 12,
        color = "grey30",
        margin = margin(t = 10)
      )
    )
  
  plot_name <- paste0(
    "country_",
    make_slug(broad_category),
    "_rates_",
    age_filename(age_group)
  )
  
  ggsave(
    filename = file.path(
      PATH_OUTDIR,
      "figures",
      paste0(plot_name, ".pdf")
    ),
    plot = country_plot,
    width = 14,
    height = 10,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  country_plot
}

plot_combinations <- country_death_rates %>%
  distinct(Broad_Category, Age_Group) %>%
  arrange(Broad_Category, Age_Group)

country_plots <- plot_combinations %>%
  mutate(
    plot_name = paste0(
      "country_",
      map_chr(Broad_Category, make_slug),
      "_rates_",
      map_chr(Age_Group, age_filename)
    ),
    plot = map2(
      Broad_Category,
      Age_Group,
      make_country_plot
    )
  )

# Cluster specific country rates
cluster_rates_country <- rates_by_sex_country %>%
  group_by(country_name, Year, Age_Group, Condition_Cluster) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pop_country, by = c("country_name", "Year", "Age_Group")) %>% 
  filter(population > 0) %>%
  mutate(
    rate_per_100k = deaths / population * 1e5,
    Age_Group     = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19'))
  )

country_injury_clusters_ages <- analytic_injury_age %>% 
  group_by(country_name, Year, Condition_Cluster, injury_cause, Age_Group) %>%
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pop_country, by = c("country_name", "Year", "Age_Group")) %>% 
  mutate(rate_per_100k = deaths / population * 100000) %>% 
  dplyr::select(-Condition_Cluster) %>% 
  rename(Condition_Cluster = injury_cause) %>% 
  mutate(Age_Group = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19')))

all_country_cluster_rates <- cluster_rates_country %>%
  dplyr::union(country_injury_clusters_ages) %>% 
  filter(Condition_Cluster != "Unintentional Injury") %>% 
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == "Drowning" ~ "Unintentional Injury: Drowning",
    Condition_Cluster == "Firearms" ~ "Unintentional Injury: Firearms",
    Condition_Cluster == "Motor vehicle crashes" ~ "Unintentional Injury: Motor Vehicle Accidents",
    Condition_Cluster == "Poisoning/overdose" ~ "Unintentional Injury: Poisoning/Overdose",
    Condition_Cluster == "Other" ~ "Unintentional Injury: Other",
    Condition_Cluster == "Diarrheal/gastrointestinal infections" ~ "Diarrheal/Gastrointestinal Infections",
    Condition_Cluster == 'Cardiovascular disorders' ~ 'Cardiovascular Disorders',
    Condition_Cluster == 'Neuro-psychiatric disorders' ~ 'Neuro-Psychiatric Disorders',
    Condition_Cluster == 'Gastrointestinal/Hepatic disorders, Non-Infectious' ~ 'Gastrointestinal/Hepatic Disorders',
    Condition_Cluster == 'Genitourinary disorders, Non-Infectious' ~ 'Genitourinary Disorders',
    Condition_Cluster == 'Prematurity/neonatal disorders' ~ 'Prematurity/Neonatal Disorders',
    Condition_Cluster == 'Birth asphyxia/birth trauma' ~ 'Birth Asphyxia/Birth Trauma',
    TRUE ~ Condition_Cluster
  ))

country_cluster_rates <- all_country_cluster_rates %>% 
  filter(!Condition_Cluster %in% c('Ill-defined and unassigned', 'Miscoded'))
write.csv(country_cluster_rates, file.path(PATH_OUTDIR, "tables", "country_cluster_rates.csv"), row.names = FALSE)

# One function for all cluster country-specific rate plots
make_country_plot_cluster <- function(condition_cluster, age_group) {
  
  country_plot <- country_cluster_rates %>%
    filter(
      Condition_Cluster == condition_cluster,
      Age_Group == age_group
    ) %>%
    ggplot(
      aes(
        x = Year,
        y = rate_per_100k,
        color = Age_Group,
        group = Age_Group
      )
    ) +
    geom_rect(
      data = missing_periods,
      aes(
        xmin = start_year - 0.5,
        xmax = end_year + 0.5,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "#F4C95D",
      alpha = 0.25
    ) +
    geom_line(linewidth = 0.9) +
    icd_vlines() +
    scale_color_manual(
      name = NULL,
      labels = c(
        "Age_lt1" = "<1-Year-Old Infants",
        "Age_1to19" = "1- to 19-Year-Old Individuals"
      ),
      values = c(
        "Age_lt1" = "#79AF97FF",
        "Age_1to19" = "#9E89B9"
      ),
      breaks = c("Age_lt1", "Age_1to19")
    ) +
    scale_x_continuous(
      breaks = seq(1950, 2030, 10)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    facet_wrap(
      ~ country_name,
      ncol = 3,
      scales = "free_y"
    ) +
    labs(
      x = "Year",
      y = "Deaths per 100,000",
      caption = paste0(
        "Note. Shaded areas indicate years with unavailable WHO-MD data: ",
        "Australia (2005), Belgium (2023), Denmark (2023), ",
        "United Kingdom (2000), and Norway (2017-2023).")
    ) +
    guides(
      color = guide_legend(
        title.position = "left",
        title.hjust = 0.5,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = 1.7
        )
      )
    ) +
    theme_mortality() +
    theme(
      strip.background = element_rect(
        fill = "grey90"
      ),
      strip.text = element_text(
        face = "bold",
        size = 14
      ),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.text = element_text(size = 12),
      panel.grid = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      
      plot.caption = element_text(
        hjust = 0,
        size = 12,
        color = "grey30",
        margin = margin(t = 10)
      )
    )
  
  plot_name <- paste0(
    "country_",
    make_slug(condition_cluster),
    "_rates_",
    age_filename(age_group)
  )
  
  ggsave(
    filename = file.path(
      PATH_OUTDIR,
      "figures",
      paste0(plot_name, ".pdf")
    ),
    plot = country_plot,
    width = 14,
    height = 10,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  country_plot
}

plot_combinations_cluster <- country_cluster_rates %>%
  distinct(Condition_Cluster, Age_Group) %>%
  arrange(Condition_Cluster, Age_Group)

country_plots_cluster <- plot_combinations_cluster %>%
  mutate(
    plot_name = paste0(
      "country_",
      map_chr(Condition_Cluster, make_slug),
      "_rates_",
      map_chr(Age_Group, age_filename)
    ),
    plot = map2(
      Condition_Cluster,
      Age_Group,
      make_country_plot_cluster
    )
  )

# ── Rate Ratios ───────────────────────────────────────────────────────────────
# Time period 1999-2023
hmd_imputation_final_ages_4cat <- hmd_imputation_final %>% 
  mutate(
    Age_Group_Broad = case_when(
      Age_Group == "Age_lt1" ~ "Age_lt1",
      Age_Group %in% c("Age_1to9", "Age_10to14", "Age_15to19") ~ "Age_1to19",
      TRUE ~ Age_Group
    )
  ) %>%
  group_by(Country, Year, Age_Group, Age_Group_Broad, Sex, imputed_hmd) %>%
  summarise(
    Population = sum(Population, na.rm = TRUE),
    Deaths_AllCause = sum(Deaths_AllCause, na.rm = TRUE),
    .groups = "drop"
  )

analytic_ages_4cat <- who_intermed %>%
  mutate(
    Age_Group_Broad = case_when(
      Age_Group == "Age_lt1" ~ "Age_lt1",
      Age_Group %in% c("Age_1to9", "Age_10to14", "Age_15to19") ~ "Age_1to19",
      TRUE ~ Age_Group
    )
  ) %>%
  group_by(Country, Year, icd_version, Age_Group, Age_Group_Broad, Sex, Condition_Cluster, Broad_Category) %>%
  summarise(
    Deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(who_code = Country) %>%
  left_join(COUNTRY_LOOKUP, by = c("who_code")) %>%
  left_join(
    hmd_imputation_final_ages_4cat,
    by = c("pop_name" = "Country", "Year", "Sex", "Age_Group", "Age_Group_Broad")
  ) %>%
  mutate(Age_Group = factor(Age_Group, levels = c('Age_lt1', 'Age_1to9', 'Age_10to14', 'Age_15to19')))

analytic_ages_firearms_4cat <- who_intermed_decomp %>%
  mutate(
    Age_Group_Broad = case_when(
      Age_Group == "Age_lt1" ~ "Age_lt1",
      Age_Group %in% c("Age_1to9", "Age_10to14", "Age_15to19") ~ "Age_1to19",
      TRUE ~ Age_Group
    )
  ) %>%
  group_by(Country, Year, icd_version, Age_Group, Sex, Condition_Cluster, Broad_Category) %>%
  summarise(
    Deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(who_code = Country) %>% 
  left_join(COUNTRY_LOOKUP, by = c("who_code")) %>%
  left_join(
    hmd_imputation_final_ages_4cat,
    by = c("pop_name" = "Country", "Year", "Sex", "Age_Group")
  ) %>% 
  mutate(Age_Group = factor(Age_Group, levels = c('Age_lt1', 'Age_1to9', 'Age_10to14', 'Age_15to19')))


analytic_population_4cat <- hmd_imputation_final_ages_4cat %>% 
  left_join(COUNTRY_LOOKUP, by = c('Country' = 'pop_name'))

pop_4cat <- analytic_population_4cat %>% 
  filter(!(Country == 'GRB_NP' & Year == 2000)) %>% 
  filter(!(Country == 'AUS' & Year == 2005),
         !(Country == 'BEL' & Year == 2023),
         !(Country == 'DNK' & Year == 2023),
         !(Country == 'NOR' & Year %in% c(2017, 2018, 2019, 2020, 2021, 2022, 2023))) %>% 
  distinct(region, Year, Age_Group, Sex, Population) %>%
  group_by(region, Year, Age_Group, Sex) %>%
  summarise(population = sum(Population), .groups = "drop")

overall_pop_4cat <- pop_4cat %>% 
  group_by(Year, Age_Group) %>%
  summarise(population = sum(population), .groups = "drop")

cluster_rates_firearms_4cat <- analytic_ages_firearms_4cat %>%
  dplyr::select(-Population) %>% 
  group_by(region, Year, Sex, Age_Group, Age_Group_Broad, Condition_Cluster, Broad_Category) %>%
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pop_4cat, by = c("region", "Year", "Age_Group", "Sex")) %>% 
  filter(population > 0) %>%
  mutate(
    rate_per_100k = deaths / population * 100000,
    Age_Group     = factor(Age_Group, levels = c('Age_lt1', 'Age_1to9', 'Age_10to14', 'Age_15to19'))
  )

cluster_rr_table <- cluster_rates_firearms_4cat %>%
  filter(Broad_Category != "Other") %>%
  dplyr::select(-Broad_Category) %>% 
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == "Drowning" ~ "Unintentional Injury: Drowning",
    Condition_Cluster == "Motor vehicle crashes" ~ "Unintentional Injury: Motor Vehicle Accidents",
    Condition_Cluster == "Poisoning/overdose" ~ "Unintentional Injury: Poisoning/Overdose",
    Condition_Cluster == "Other" ~ "Unintentional Injury: Other",
    TRUE ~ Condition_Cluster
  ))


domain_rates_by_sex_4cat <- analytic_ages_4cat %>%
  dplyr::select(-Population) %>% 
  group_by(region, Year, Sex, Age_Group, Age_Group_Broad, Broad_Category) %>%
  summarise(
    deaths = sum(Deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  left_join(pop_4cat, by = c("region", "Year", "Age_Group", "Sex")) %>% 
  filter(population > 0) %>%
  mutate(
    Age_Group     = factor(Age_Group, levels = c('Age_lt1', 'Age_1to9', 'Age_10to14', 'Age_15to19'))
  )

domain_model_data <- domain_rates_by_sex_4cat %>%
  mutate(
    Year = as.integer(Year),
    
    region = factor(region, levels = c("OECD13", "US")),
    Sex = factor(Sex),
  ) %>%
  filter(
    Broad_Category != 'Other',
    Year >= 1999,
    !is.na(Age_Group),
    !is.na(deaths),
    !is.na(population),
    population > 0
  ) %>% 
  mutate(
    rate_per_100k = deaths / population * 100000
  )

cluster_model_data <- cluster_rr_table %>%
  mutate(
    Year = as.integer(Year),
    
    region = factor(region, levels = c("OECD13", "US")),
    Sex = factor(Sex),
  ) %>%
  filter(
    Year >= 1999,
    !is.na(Age_Group),
    !is.na(deaths),
    !is.na(population),
    population > 0
  ) %>% 
  mutate(
    rate_per_100k = deaths / population * 100000
  )

fit_rr_model <- function(data) {
  
  data <- droplevels(data)
  
    number_age_groups <- data %>%
    summarise(
      n = n_distinct(
        as.character(Age_Group),
        na.rm = TRUE
      )
    ) %>%
    pull(n)
  
  if (number_age_groups == 1) {
    
    model <- fepois(
      deaths ~ region + Sex + factor(Year),
      offset = ~ log(population),
      cluster = ~ Year,
      data = data
    )
    
  } else if (number_age_groups > 1) {
    
    model <- fepois(
      deaths ~ region + Sex + Age_Group + factor(Year),
      offset = ~ log(population),
      cluster = ~ Year,
      data = data
    )
    
  } else {
    
    stop("No valid detailed age groups were found in this subset.")
  }
  
  return(model)
}


model_results <- domain_model_data %>%
  group_by(Broad_Category, Age_Group_Broad) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    model = map(data, fit_rr_model)
  )

rr_results <- model_results %>%
  mutate(
    result = map(
      model,
      ~ broom::tidy(
        .x,
        conf.int = TRUE
      ) %>%
        filter(term == "regionUS") %>%
        mutate(
          RR = exp(estimate),
          CI_lower = exp(conf.low),
          CI_upper = exp(conf.high)
        )
    )
  ) %>%
  dplyr::select(
    Broad_Category,
    Age_Group_Broad,
    result
  ) %>%
  unnest(result) %>%
  transmute(
    Broad_Category,
    Age_Group_Broad,
    RR,
    CI_lower,
    CI_upper,
    p_value = round(p.value, digits = 3)
  ) %>%
  arrange(Age_Group_Broad, Broad_Category)

write.csv(rr_results, file.path(PATH_OUTDIR, "tables", "broad_rr.csv"), row.names = FALSE)

cluster_results <- cluster_model_data %>%
  group_by(Condition_Cluster, Age_Group_Broad) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    model = map(data, fit_rr_model)
  )

cluster_rr_results <- cluster_results %>%
  mutate(
    result = map(
      model,
      ~ broom::tidy(
        .x,
        conf.int = TRUE
      ) %>%
        filter(term == "regionUS") %>%
        mutate(
          RR = exp(estimate),
          CI_lower = exp(conf.low),
          CI_upper = exp(conf.high)
        )
    )
  ) %>%
  dplyr::select(
    Condition_Cluster,
    Age_Group_Broad,
    result
  ) %>%
  unnest(result) %>%
  transmute(
    Condition_Cluster,
    Age_Group_Broad,
    RR,
    CI_lower,
    CI_upper,
    p_value = round(p.value, digits = 3)
  ) %>%
  arrange(Age_Group_Broad, Condition_Cluster)

write.csv(cluster_rr_results, file.path(PATH_OUTDIR, "tables", "cluster_rr.csv"), row.names = FALSE)

# Time period 2014-2023
domain_model_data_2014 <- domain_rates_by_sex_4cat %>%
  mutate(
    Year = as.integer(Year),
    
    region = factor(region, levels = c("OECD13", "US")),
    Sex = factor(Sex),
  ) %>%
  filter(
    Broad_Category != 'Other',
    Year >= 2014,
    !is.na(Age_Group),
    !is.na(deaths),
    !is.na(population),
    population > 0
  ) %>% 
  mutate(
    rate_per_100k = deaths / population * 100000
  )

cluster_model_data_2014 <- cluster_rr_table %>%
  mutate(
    Year = as.integer(Year),
    
    region = factor(region, levels = c("OECD13", "US")),
    Sex = factor(Sex),
  ) %>%
  filter(
    Year >= 2014,
    !is.na(Age_Group),
    !is.na(deaths),
    !is.na(population),
    population > 0
  ) %>% 
  mutate(
    rate_per_100k = deaths / population * 100000
  )


model_results_2014 <- domain_model_data_2014 %>%
  group_by(Broad_Category, Age_Group_Broad) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    model = map(data, fit_rr_model)
  )

rr_results_2014 <- model_results_2014 %>%
  mutate(
    result = map(
      model,
      ~ broom::tidy(
        .x,
        conf.int = TRUE
      ) %>%
        filter(term == "regionUS") %>%
        mutate(
          RR = exp(estimate),
          CI_lower = exp(conf.low),
          CI_upper = exp(conf.high)
        )
    )
  ) %>%
  dplyr::select(
    Broad_Category,
    Age_Group_Broad,
    result
  ) %>%
  unnest(result) %>%
  transmute(
    Broad_Category,
    Age_Group_Broad,
    RR,
    CI_lower,
    CI_upper,
    p_value = round(p.value, digits = 3)
  ) %>%
  arrange(Age_Group_Broad, Broad_Category)

write.csv(rr_results_2014, file.path(PATH_OUTDIR, "tables", "broad_rr_2014.csv"), row.names = FALSE)

cluster_results_2014 <- cluster_model_data_2014 %>%
  group_by(Condition_Cluster, Age_Group_Broad) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    model = map(data, fit_rr_model)
  )

cluster_rr_results_2014 <- cluster_results_2014 %>%
  mutate(
    result = map(
      model,
      ~ broom::tidy(
        .x,
        conf.int = TRUE
      ) %>%
        filter(term == "regionUS") %>%
        mutate(
          RR = exp(estimate),
          CI_lower = exp(conf.low),
          CI_upper = exp(conf.high)
        )
    )
  ) %>%
  dplyr::select(
    Condition_Cluster,
    Age_Group_Broad,
    result
  ) %>%
  unnest(result) %>%
  transmute(
    Condition_Cluster,
    Age_Group_Broad,
    RR,
    CI_lower,
    CI_upper,
    p_value = round(p.value, digits = 3)
  ) %>%
  arrange(Age_Group_Broad, Condition_Cluster)

write.csv(cluster_rr_results_2014, file.path(PATH_OUTDIR, "tables", "cluster_rr_2014.csv"), row.names = FALSE)


# =============================================================================
# APC INTERVALS
# =============================================================================

# ── APC computation ───────────────────────────────────────────────────────────
compute_apc <- function(df, rate_col = "rate_per_100k") {
  df <- df %>%
    filter(.data[[rate_col]] > 0) %>%
    arrange(Year)
  n <- nrow(df)
  if (n < 3L) {
    return(tibble(
      n_obs = n, slope = NA_real_, APC = NA_real_,
      APC_lo = NA_real_, APC_hi = NA_real_,
      flag = "too few observations"
    ))
  }
  df <- df %>% mutate(log_rate = log(.data[[rate_col]]))
  mod <- lm(log_rate ~ Year, data = df)
  sl <- coef(mod)["Year"]
  se <- summary(mod)$coefficients["Year", "Std. Error"]
  ci <- confint(mod, "Year", level = 0.95)
  tibble(
    n_obs  = n,
    slope  = sl,
    APC    = (exp(sl) - 1) * 100,
    APC_lo = (exp(ci[1]) - 1) * 100,
    APC_hi = (exp(ci[2]) - 1) * 100,
    se     = (APC_hi - APC_lo) / (2*1.96),
    flag   = if_else(n < 5L, "fewer than 5 observations", "ok")
  )
}

cluster_apc_results <- map_dfr(names(APC_INTERVALS), function(int_name) {
  yrs <- APC_INTERVALS[[int_name]]
  final_all_cluster_rates %>%
    filter(Year %in% yrs) %>%
    group_by(region, Condition_Cluster, Age_Group) %>%
    group_modify(~ compute_apc(.x)) %>%
    ungroup() %>%
    mutate(interval = int_name)
})

domain_apc_results <- map_dfr(names(APC_INTERVALS), function(int_name) {
  yrs <- APC_INTERVALS[[int_name]]
  domain_rates %>%
    filter(Year %in% yrs,
           Broad_Category != 'Other') %>%
    group_by(region, Broad_Category, Age_Group) %>%
    group_modify(~ compute_apc(.x)) %>%
    ungroup() %>%
    mutate(interval = int_name)
})


# Differential APC (OECD13 - US)
cluster_apc_diff <- cluster_apc_results %>%
  filter(!Condition_Cluster %in% c('Ill-defined and unassigned', 'Miscoded')) %>% 
  dplyr::select(region, Condition_Cluster, Age_Group, interval, APC, APC_lo, APC_hi, se) %>%
  pivot_wider(
    names_from  = region,
    values_from = c(APC, APC_lo, APC_hi, se),
    values_fill = 0
  ) %>%
  mutate(se_diff = sqrt(se_US^2 + se_OECD13^2)) %>% 
  mutate(
    diff_APC    = APC_OECD13 - APC_US,
    diff_APC_lo = diff_APC - (1.96 * se_diff),
    diff_APC_hi = diff_APC + (1.96 * se_diff)
  )

write.csv(cluster_apc_diff, file.path(PATH_OUTDIR, "tables", "cluster_apc.csv"), row.names = FALSE)

domain_apc_diff <- domain_apc_results %>%
  dplyr::select(region, Broad_Category, Age_Group, interval, APC, APC_lo, APC_hi, se) %>%
  pivot_wider(
    names_from  = region,
    values_from = c(APC, APC_lo, APC_hi, se),
    values_fill = 0
  ) %>%
  mutate(se_diff = sqrt(se_US^2 + se_OECD13^2)) %>% 
  mutate(
    diff_APC    = APC_OECD13 - APC_US,
    diff_APC_lo = diff_APC - (1.96 * se_diff),
    diff_APC_hi = diff_APC + (1.96 * se_diff)
  )

write.csv(domain_apc_diff, file.path(PATH_OUTDIR, "tables", "broad_apc.csv"), row.names = FALSE)


# =============================================================================
# DECOMPOSITION
# =============================================================================

# Get overall difference US - OEDC13
overall_diff <- cluster_rates_decomp %>%
  dplyr::select(-Broad_Category) %>% 
  filter(!Condition_Cluster %in% c("Unintentional Injury", "Miscoded", "Ill-defined and unassigned"),
         Year >= 1999) %>% 
  group_by(region, Age_Group, Year) %>% 
  summarise(deaths = sum(deaths, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(pop, by = c("region", "Year", "Age_Group")) %>% 
  mutate(rate_per_100k = deaths / population * 100000) %>% 
  dplyr::select(-deaths, -population) %>% 
  pivot_wider(names_from = region, values_from = rate_per_100k) %>% 
  mutate(overall_rate_diff = US - OECD13)

cluster_diff <- cluster_rates_decomp %>%
  dplyr::select(-Broad_Category) %>% 
  filter(!Condition_Cluster %in% c("Unintentional Injury", "Miscoded", "Ill-defined and unassigned"),
         Year >= 1999) %>% 
  group_by(region, Age_Group, Year, Condition_Cluster) %>% 
  summarise(deaths = sum(deaths, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(pop, by = c("region", "Year", "Age_Group")) %>% 
  mutate(rate_per_100k = deaths / population * 100000) %>% 
  dplyr::select(-deaths, -population) %>% 
  pivot_wider(names_from = region, values_from = rate_per_100k) %>% 
  mutate(cluster_rate_diff = US - OECD13) %>% 
  left_join(dplyr::select(overall_diff, Age_Group, Year, overall_rate_diff), by = c('Age_Group', 'Year')) %>% 
  mutate(cluster_pct = cluster_rate_diff / overall_rate_diff * 100)

lt1_decomp <- cluster_diff %>% 
  filter(Age_Group == 'Age_lt1') %>% 
  mutate(Condition_Cluster = if_else(!Condition_Cluster %in% 
                                       c('Sudden Unexpected Infant Deaths (SUID)',
                                         'Congenital Anomalies',
                                         'Perinatal disorders',
                                         'Cardiovascular disorders'), 'Other', Condition_Cluster)) %>% 
  group_by(Year, Condition_Cluster, Age_Group, overall_rate_diff) %>%
  summarise(
    US = sum(US),
    OECD13 = sum(OECD13),
    cluster_rate_diff = sum(cluster_rate_diff),
    cluster_pct = sum(cluster_pct),
    .groups = "drop"
  )

age1to19_decomp <- cluster_diff %>% 
  filter(Age_Group == 'Age_1to19') %>% 
  mutate(Condition_Cluster = if_else(!Condition_Cluster %in% 
                                       c('Firearms',
                                         'Motor vehicle crashes',
                                         'Poisoning/overdose',
                                         'Drowning'), 'Other', Condition_Cluster)) %>% 
  group_by(Year, Condition_Cluster, Age_Group, overall_rate_diff) %>%
  summarise(
    US = sum(US),
    OECD13 = sum(OECD13),
    cluster_rate_diff = sum(cluster_rate_diff),
    cluster_pct = sum(cluster_pct),
    .groups = "drop"
  )

all_decomp <- lt1_decomp %>% 
  dplyr::union(age1to19_decomp) %>% 
  mutate(Age_Group = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19'))) %>% 
  mutate(cluster_pct = if_else(cluster_pct < 0, 0, cluster_pct))


write.csv(all_decomp, file.path(PATH_OUTDIR, "tables", "fig4_table.csv"), row.names = FALSE)

plot_dat <- all_decomp %>%
  mutate(
    Age_Group = recode(
      Age_Group,
      "Age_lt1"  = "A. <1-Year-Old Infants",
      "Age_1to19" = "B. 1- to 19-Year-Old Individuals"
    )
  ) %>% 
  mutate(Condition_Cluster = case_when(
    Condition_Cluster == 'Cardiovascular disorders' ~ 'Cardiovascular Disorders',
    Condition_Cluster == 'Perinatal disorders' ~ 'Perinatal Disorders',
    Condition_Cluster == 'Motor vehicle crashes' ~ 'Motor Vehicle Accidents',
    Condition_Cluster == 'Poisoning/overdose' ~ 'Poisoning/Overdose',
    TRUE ~ Condition_Cluster)
  )

# Panel A
p_infants <- plot_dat %>%
  mutate(Condition_Cluster = factor(Condition_Cluster, levels = c('Other',
                                                                  'Cardiovascular Disorders',
                                                                  'Congenital Anomalies',
                                                                  'Sudden Unexpected Infant Deaths (SUID)',
                                                                  'Perinatal Disorders'))) %>% 
  filter(Age_Group == "A. <1-Year-Old Infants") %>%
  ggplot(
    aes(
      x = Year,
      y = cluster_pct,
      fill = Condition_Cluster
    )
  ) +
  geom_area(color = "white", linewidth = 0.2) +
  facet_wrap(~Age_Group) +
  scale_fill_manual(
    values = c(
      "Perinatal Disorders" = "#2F5763",
      "Sudden Unexpected Infant Deaths (SUID)" = "#F7941D",
      "Congenital Anomalies" = "#00ADEF",
      "Cardiovascular Disorders" = "#8C7F70",
      "Other" = "#87AAB9"
    )
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 25),
    labels = scales::label_percent(scale = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    ylim = c(0, 100),
    expand = FALSE
  ) +
  labs(
    x = "Year",
    y = "Share of Total Excess Deaths (%)",
    fill = NULL
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  theme_minimal(base_size = 12) +
    theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    axis.line.x.bottom = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.line.y.left = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.x.bottom = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.y.left = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.length = grid::unit(0.1, "cm"),
    
    strip.text = element_text(
      hjust = 0,
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    
    legend.position = "bottom",
    legend.justification = c(0, 0.5),
    legend.box.just = "left",
    legend.margin = margin(t = 8, r = 0, b = 0, l = 0)
  )

# Panel B
p_1to19 <- plot_dat %>%
  mutate(Condition_Cluster = factor(Condition_Cluster, levels = c('Poisoning/Overdose',
                                                                  'Drowning',
                                                                  'Other',
                                                                  'Motor Vehicle Accidents',
                                                                  'Firearms'))) %>% 
  filter(Age_Group == "B. 1- to 19-Year-Old Individuals") %>%
  ggplot(
    aes(
      x = Year,
      y = cluster_pct,
      fill = Condition_Cluster
    )
  ) +
  geom_area(color = "white", linewidth = 0.2) +
  facet_wrap(~Age_Group) +
  scale_fill_manual(
    values = c(
      "Firearms" = "#E6417B",
      "Motor Vehicle Accidents" = "#9E89B9",
      "Poisoning/Overdose" = "#8C922A",
      "Drowning" = "#D23C3D",
      "Other" = "#00744E"
    )
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 25),
    labels = scales::label_percent(scale = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    ylim = c(0, 100),
    expand = FALSE
  ) +
  labs(
    x = "Year",
    y = NULL,
    fill = NULL
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    axis.line.x.bottom = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.line.y.left = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.x.bottom = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.y.left = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.length = grid::unit(0.1, "cm"),
    
    strip.text = element_text(
      hjust = 0,
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    
      legend.position = "bottom",
      legend.justification = c(0, 0.5),
      legend.box.just = "left",
      legend.margin = margin(t = 8, r = 0, b = 0, l = 0)
  )

fig4 <- (p_infants | p_1to19) +
  plot_layout(
    ncol = 2,
    guides = "keep"
  ) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(
        hjust = 0,
        face = "bold",
        size = 14,
        lineheight = 1.1,
        margin = margin(b = 2)
      ),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  )

ggsave(
  file.path(PATH_OUTDIR, "figures", "fig4_decomp.pdf"),
  plot = fig4,
  width = 16,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

domain_rates_decomp <- rates_by_sex_decomp %>%
  group_by(region, Year, Age_Group, Broad_Category) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pop, by = c("region", "Year", "Age_Group")) %>% 
  filter(population > 0) %>%
  mutate(
    rate_per_100k = deaths / population * 100000,
    Age_Group     = factor(Age_Group, levels = c('Age_lt1', 'Age_1to19'))
  )

overall_domain_diff <- domain_rates_decomp %>%
  filter(Year >= 1999,
         Broad_Category != 'Other') %>% 
  group_by(region, Age_Group, Year) %>% 
  summarise(deaths = sum(deaths, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(pop, by = c("region", "Year", "Age_Group")) %>% 
  mutate(rate_per_100k = deaths / population * 100000) %>% 
  dplyr::select(-deaths, -population) %>% 
  pivot_wider(names_from = region, values_from = rate_per_100k) %>% 
  mutate(overall_rate_diff = US - OECD13)


domain_diff <- domain_rates_decomp %>%
  filter(Year >= 1999,
         Broad_Category != 'Other') %>% 
  group_by(region, Age_Group, Year, Broad_Category) %>% 
  summarise(deaths = sum(deaths, na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(pop, by = c("region", "Year", "Age_Group")) %>% 
  mutate(rate_per_100k = deaths / population * 100000) %>% 
  dplyr::select(-deaths, -population) %>% 
  pivot_wider(names_from = region, values_from = rate_per_100k) %>% 
  mutate(domain_rate_diff = US - OECD13) %>% 
  left_join(dplyr::select(overall_domain_diff, Age_Group, Year, overall_rate_diff), by = c('Age_Group', 'Year')) %>% 
  mutate(domain_pct = domain_rate_diff / overall_rate_diff * 100)


write.csv(domain_diff, file.path(PATH_OUTDIR, "tables", "broad_excess.csv"), row.names = FALSE)

# =============================================================================
# ADDITIONAL SUPPLEMENTAL TABLES
# =============================================================================

# ── ICD code counts ───────────────────────────────────────────────────────────
icd_counts <- who_long_restrict %>% 
  mutate(
    Age_Group = case_when(
      Age_Group == "Age_lt1" ~ "Age_lt1",
      Age_Group %in% c("Age_1to9", "Age_10to14", "Age_15to19") ~ "Age_1to19",
      TRUE ~ Age_Group)
  ) %>%
  left_join(COUNTRY_LOOKUP, by = c('Country' = 'who_code')) %>% 
  group_by(country_name, Year, Cause, icd_version, Age_Group, reassigned_cause, injury_cause) %>% 
  summarise(deaths = sum(deaths), .groups = "drop") %>% 
  left_join(DOMAIN_MAP, by = c('reassigned_cause' = 'cause_abbrev')) %>% 
  rename(Condition_Cluster = cluster, Broad_Category = domain) %>% 
  dplyr::select(-reassigned_cause) %>% 
  pivot_wider(names_from = 'country_name', values_from = deaths)

icd_counts_oecd13 <- who_long_restrict %>% 
  mutate(
    Age_Group = case_when(
      Age_Group == "Age_lt1" ~ "Age_lt1",
      Age_Group %in% c("Age_1to9", "Age_10to14", "Age_15to19") ~ "Age_1to19",
      TRUE ~ Age_Group)
  ) %>%
  left_join(COUNTRY_LOOKUP, by = c('Country' = 'who_code')) %>% 
  filter(region == 'OECD13') %>% 
  group_by(region, Year, Cause, icd_version, Age_Group, reassigned_cause, injury_cause) %>% 
  summarise(deaths = sum(deaths), .groups = "drop") %>% 
  left_join(DOMAIN_MAP, by = c('reassigned_cause' = 'cause_abbrev')) %>% 
  rename(Condition_Cluster = cluster, Broad_Category = domain) %>% 
  dplyr::select(-reassigned_cause) %>% 
  pivot_wider(names_from = 'region', values_from = deaths)

all_icd_counts <- icd_counts %>% 
  left_join(icd_counts_oecd13) %>% 
  left_join(dplyr::select(fk_codesets_desc, icd_code, icd_version, Description), 
            by = c('Cause' = 'icd_code', 'icd_version'))


write.csv(all_icd_counts, file.path(PATH_OUTDIR, "tables", "icd_counts.csv"), row.names = FALSE)

# ── Population counts ─────────────────────────────────────────────────────────
pop_year_country <- analytic_population %>% 
  distinct(country_name, Year, Age_Group, Sex, Population, imputed_hmd) %>%
  group_by(country_name, Year, imputed_hmd) %>%
  summarise(population = sum(Population), .groups = "drop")

pop_year_oecd <- analytic_population %>% 
  filter(region == 'OECD13') %>% 
  distinct(country_name = region, Year, Age_Group, Sex, Population) %>%
  group_by(country_name, Year) %>%
  summarise(population = sum(Population), .groups = "drop") %>% 
  mutate(imputed_hmd = NA)

pop_year <- pop_year_country %>% 
  rbind(pop_year_oecd)

write.csv(pop_year, file.path(PATH_OUTDIR, "tables", "pop_year.csv"), row.names = FALSE)

pop_year_age_country <- analytic_population %>% 
  distinct(country_name, Year, Age_Group, Sex, Population, imputed_hmd) %>%
  group_by(country_name, Year, Age_Group, imputed_hmd) %>%
  summarise(population = sum(Population), .groups = "drop")

pop_year_age_oecd <- analytic_population %>% 
  filter(region == 'OECD13') %>% 
  distinct(country_name = region, Year, Age_Group, Sex, Population) %>%
  group_by(country_name, Year, Age_Group) %>%
  summarise(population = sum(Population), .groups = "drop") %>% 
  mutate(imputed_hmd = NA)

pop_year_age <- pop_year_age_country %>% 
  rbind(pop_year_age_oecd)

write.csv(pop_year_age, file.path(PATH_OUTDIR, "tables", "pop_year_age.csv"), row.names = FALSE)

pop_year_age_sex_country <- analytic_population %>% 
  distinct(country_name, Year, Age_Group, Sex, Population, imputed_hmd) %>%
  group_by(country_name, Year, Age_Group, Sex, imputed_hmd) %>%
  summarise(population = sum(Population), .groups = "drop")

pop_year_age_sex_oecd <- analytic_population %>% 
  filter(region == 'OECD13') %>% 
  distinct(country_name = region, Year, Age_Group, Sex, Population) %>%
  group_by(country_name, Year, Age_Group, Sex) %>%
  summarise(population = sum(Population), .groups = "drop") %>% 
  mutate(imputed_hmd = NA)

pop_year_age_sex <- pop_year_age_sex_country %>% 
  rbind(pop_year_age_sex_oecd)

write.csv(pop_year_age_sex, file.path(PATH_OUTDIR, "tables", "pop_year_age_sex.csv"), row.names = FALSE)

# ── ICD transition years ──────────────────────────────────────────────────────
country_icd_transition <- analytic_ages %>% 
  group_by(country_name, icd_version) %>% 
  summarise(first_year = min(Year),
            last_year = max(Year), .groups = "drop") %>% 
  arrange(country_name, as.numeric(icd_version))

write.csv(country_icd_transition, file.path(PATH_OUTDIR, "tables", "country_icd_transition.csv"), row.names = FALSE)


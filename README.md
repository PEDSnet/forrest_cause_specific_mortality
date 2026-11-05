# **Cause-Specific Mortality of Children in the US vs Peer Nations, 1954–2023**

Christopher B. Forrest, MD, PhD¹

Allison M. Zelinski, BS, MS¹

Elizabeth R. Wolf, MD, MPH²

Mari K. Noble, BS¹

Lauren J. Koenigsberg, BA¹

¹ Applied Clinical Research Center, Children’s Hospital of Philadelphia, Philadelphia, PA

² Department of Pediatrics, Virginia Commonwealth University, Richmond, VA

**Corresponding Author:**

  Christopher B. Forrest

  Children’s Hospital of Philadelphia

  2716 South St, Suite 11-473

  Philadelphia, PA 19146

  E: [forrestc\@chop.edu](mailto:forrestc@chop.edu)

  M: 267-909-4627

## Scripts and directories

``` sh
forrest_cause_specific_mortality/
├── .gitignore
├── README.md
├── cause_specific_mortality_v1.R
├── data/
│   └── raw/
│       ├── codesets/
│       │   └── harmonized_codesets_forrest_v1.xlsx
│       ├── HMD/
│       │   ├── AUS/
│       │   │   └── .placeholder
│       │   ├── BEL/
│       │   │   └── .placeholder
│       │   ├── CAN/
│       │   │   └── .placeholder
│       │   ├── CHE/
│       │   │   └── .placeholder
│       │   ├── DNK/
│       │   │   └── .placeholder
│       │   ├── ESP/
│       │   │   └── .placeholder
│       │   ├── FIN/
│       │   │   └── .placeholder
│       │   ├── FRATNP/
│       │   │   └── .placeholder
│       │   ├── GRB_NP/
│       │   │   └── .placeholder
│       │   ├── ITA/
│       │   │   └── .placeholder
│       │   ├── NLD/
│       │   │   └── .placeholder
│       │   ├── NOR/
│       │   │   └── .placeholder
│       │   ├── SWE/
│       │   │   └── .placeholder
│       │   └── USA/
│       │       └── .placeholder
│       └── WHO/
│           └── .placeholder
└── output/
    ├── figures/
    │   └── .placeholder
    ├── tables/
    │   └── .placeholder
    └── tables.html
  
```

cause_specific_mortality_v1.R: R code for all analyses

data/raw/codesets/harmonized_codesets_forrest_v1.xlsx: ICD code classification

data/raw/HMD: Expected location of country-specific Human Mortality Database (HMD) files

data/raw/WHO: Expected location of cause-specific WHO Morality Database (WHO-MD) files

output/figures: Output location for figures

output/tables: Output location for tables

output/tables.html: HTML of interactive data tables used in manuscript

## How to use

1.  Obtain population denominators from [HMD website](https://www.mortality.org/) and place files in corresponding country folders in the HMD folder

    - Files for each country: Deaths_1x1.txt, Exposures_1x1.txt

2.  Obtain cause-specific mortality data from the [WHO-MD website](https://www.who.int/data/data-collection-tools/who-mortality-database), save as csvs, and place files in the WHO folder

    - Files: MortIcd7.csv, Morticd8.csv, Morticd9.csv, Morticd10_part1.csv, Morticd10_part2.csv, Morticd10_part3.csv, Morticd10_part4.csv, Morticd10_part5.csv, Morticd10_part6.csv

3.  Ensure harmonized_codesets_forrest_v1.xlsx is in the codesets folder

4.  Install the following packages:

    - dplyr, purrr, tidyr, tidyverse, scales, patchwork, fixest, broom, ggplot2, ggh4x

5.  Run the script: cause_specific_mortality_v1.R

6.  Output figures and tables will populate in their respective folders

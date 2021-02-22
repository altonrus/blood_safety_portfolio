args = commandArgs(TRUE)
iter_start = as.integer(args[1])
iter_end = as.integer(args[2])



library(data.table)
#library(googlesheets4)

## Read in and process NMC simulation PSA
# source('nmc_r_functions.r')
source('portfolio_functions.R')
source('PSA_functions.R')

# Read in files
params_PSA <- fread("param_inputs_unif.csv")
ZIKV_psa_nmc <- fread('PSA_ZIKV_NMC.csv')
WNV_psa_nmc <- fread('PSA_WNV_NMC.csv')

ZIKV_lims <- ZIKV_psa_nmc[ , quantile(NMC_per_donation, probs = c(0.01, 0.99)) ]
WNV_lims <- WNV_psa_nmc[ , quantile(NMC_per_donation, probs = c(0.01, 0.99)) ]

ZIKV_unif_nmc <- data.table(
  NMC_per_donation = runif(10000, min = ZIKV_lims[1], ZIKV_lims[2])
)

WNV_unif_nmc <- data.table(
  NMC_per_donation = runif(10000, min = WNV_lims[1], WNV_lims[2])
)

donorGroups<- fread('donor_groups_zip3.csv')

# Sys.time()
Process_and_save_PSA(iter_start, iter_end, params_PSA, ZIKV_unif_nmc, WNV_unif_nmc, donorGroup,
                     fname_prefix = "psa_zip3_unif_")
# Sys.time()


args = commandArgs(TRUE)
iter_start = as.integer(args[1])
iter_end = as.integer(args[2])

# iter_start = 2210
# iter_end = 2211

library(data.table)
#library(googlesheets4)

## Read in and process NMC simulation PSA
# source('nmc_r_functions.r')
source('portfolio_functions.R')
source('PSA_functions.R')

# Read in files
params_PSA <- fread("param_inputs.csv")
ZIKV_psa_nmc <- fread('PSA_ZIKV_NMC.csv')
WNV_psa_nmc <- fread('PSA_WNV_NMC.csv')
donorGroups<- fread('donor_groups_zip3.csv')

# Sys.time()
Process_and_save_PSA(iter_start, iter_end, params_PSA, ZIKV_psa_nmc, WNV_psa_nmc, donorGroup)
# Sys.time()


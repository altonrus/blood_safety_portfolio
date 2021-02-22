library(data.table)
library(readxl)

## Read in and process NMC simulation PSA
source('NMC_sim/nmc_r_functions.r')
source('portfolio_functions.R')


###
# Generate distributions for PSA from the non-simulation parameters
###


# REGULAR TRIANGULAR DISTRIBUTIONS
# Read in params
param_sheet <- data.table(read_excel("data/params.xlsx", sheet = "params"))

n_PSA = 10000

params_PSA <- gen_PSA_inputs(param_sheet[ key != 'c_k'], n_PSA, two_index = TRUE)

fwrite(params_PSA, "./param_inputs.csv")

# UNIF SCENARIO
# Read in params
param_sheet <- data.table(read_excel("data/params.xlsx", sheet = "params"))
param_sheet[ , Distribution := ifelse(Distribution == "Tri", "Unif", Distribution)]
n_PSA = 10000

params_PSA <- gen_PSA_inputs(param_sheet[ key != 'c_k'], n_PSA, two_index = TRUE)

fwrite(params_PSA, "./param_inputs_unif.csv")



###
# Generate net health benefit per iteration of the microsim PSAs
###
# ZIKV
ZIKV_psa_output <- fread('NMC_sim/output/PSA_ZIKV_combined.csv')
transmissability_RBC = as.vector(t(params_PSA[key == "transmissability_ZIKV_RBC", 4:10003]))
transmissability_PLT = as.vector(t(params_PSA[key == "transmissability_ZIKV_PLT", 4:10003]))
transmissability_FFP = as.vector(t(params_PSA[key == "transmissability_ZIKV_FFP", 4:10003]))
RBC_per_don = as.vector(t(params_PSA[key == "RBC_per_don", 4:10003]))
PLT_per_don = as.vector(t(params_PSA[key == "PLT_per_don", 4:10003]))
FFP_per_don = as.vector(t(params_PSA[key == "FFP_per_don", 4:10003]))

ZIKV_psa_nmc <- nmc_from_sim_output(sim_output = ZIKV_psa_output, transmissability_RBC, transmissability_PLT, transmissability_FFP,
                                           RBC_per_don, PLT_per_don, FFP_per_don,
                                           WTP = 1E6)
fwrite(ZIKV_psa_nmc, 'NMC_sim/output/PSA_ZIKV_NMC.csv')
# WNV
WNV_psa_output <- fread('NMC_sim/output/PSA_WNV_combined.csv')
transmissability_RBC = as.vector(t(params_PSA[key == "transmissability_WNV_RBC", 4:10003]))
transmissability_PLT = as.vector(t(params_PSA[key == "transmissability_WNV_PLT", 4:10003]))
transmissability_FFP = as.vector(t(params_PSA[key == "transmissability_WNV_FFP", 4:10003]))
RBC_per_don = as.vector(t(params_PSA[key == "RBC_per_don", 4:10003]))
PLT_per_don = as.vector(t(params_PSA[key == "PLT_per_don", 4:10003]))
FFP_per_don = as.vector(t(params_PSA[key == "FFP_per_don", 4:10003]))

WNV_psa_nmc <- nmc_from_sim_output(sim_output = WNV_psa_output, transmissability_RBC, transmissability_PLT, transmissability_FFP,
                                    RBC_per_don, PLT_per_don, FFP_per_don,
                                    WTP = 1E6)
fwrite(WNV_psa_nmc, 'NMC_sim/output/PSA_WNV_NMC.csv')



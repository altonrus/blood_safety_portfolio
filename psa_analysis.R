library(data.table)
library(googlesheets4)

## Read in and process NMC simulation PSA
source('NMC_sim/nmc_r_functions.r')
source('portfolio_functions.R')
###
# Generate distributions for PSA from the non-simulation parameters
###

# Get parameters from Google shet
# param_sheet <- data.table(read_sheet("1Yjfq0SINstVPrszWYx9uJEYcqChxvGM1RIoR8uZe-sw", sheet = "params"))
# 
# n_PSA = 10000
# 
# params_PSA <- gen_PSA_inputs(param_sheet[ key != 'c_k'], n_PSA, two_index = TRUE)
# 
# fwrite(params_PSA, "Results/param_inputs.csv")
params_PSA <- fread("Results/param_inputs.csv")

###
# Generate net health benefit per iteration of the microsim PSAs
###
# # ZIKV
# ZIKV_psa_output <- fread('NMC_sim/output/PSA_ZIKV_combined.csv')
# transmissability_RBC = as.vector(t(params_PSA[key == "transmissability_ZIKV_RBC", 4:10003]))
# transmissability_PLT = as.vector(t(params_PSA[key == "transmissability_ZIKV_PLT", 4:10003]))
# transmissability_FFP = as.vector(t(params_PSA[key == "transmissability_ZIKV_FFP", 4:10003]))
# RBC_per_don = as.vector(t(params_PSA[key == "RBC_per_don", 4:10003]))
# PLT_per_don = as.vector(t(params_PSA[key == "PLT_per_don", 4:10003]))
# FFP_per_don = as.vector(t(params_PSA[key == "FFP_per_don", 4:10003]))
# 
# ZIKV_psa_nmc <- nmc_from_sim_output(sim_output = ZIKV_psa_output, transmissability_RBC, transmissability_PLT, transmissability_FFP,
#                                            RBC_per_don, PLT_per_don, FFP_per_don,
#                                            WTP = 1E6)
# fwrite(ZIKV_psa_nmc, 'NMC_sim/output/PSA_ZIKV_NMC.csv')
# # WNV
# WNV_psa_output <- fread('NMC_sim/output/PSA_WNV_combined.csv')
# transmissability_RBC = as.vector(t(params_PSA[key == "transmissability_WNV_RBC", 4:10003]))
# transmissability_PLT = as.vector(t(params_PSA[key == "transmissability_WNV_PLT", 4:10003]))
# transmissability_FFP = as.vector(t(params_PSA[key == "transmissability_WNV_FFP", 4:10003]))
# RBC_per_don = as.vector(t(params_PSA[key == "RBC_per_don", 4:10003]))
# PLT_per_don = as.vector(t(params_PSA[key == "PLT_per_don", 4:10003]))
# FFP_per_don = as.vector(t(params_PSA[key == "FFP_per_don", 4:10003]))
# 
# WNV_psa_nmc <- nmc_from_sim_output(sim_output = WNV_psa_output, transmissability_RBC, transmissability_PLT, transmissability_FFP,
#                                     RBC_per_don, PLT_per_don, FFP_per_don,
#                                     WTP = 1E6)
# fwrite(WNV_psa_nmc, 'NMC_sim/output/PSA_WNV_NMC.csv')
ZIKV_psa_nmc <- fread('NMC_sim/output/PSA_ZIKV_NMC.csv')
WNV_psa_nmc <- fread('NMC_sim/output/PSA_WNV_NMC.csv')
donorGroups<- data.table(read_sheet("1Yjfq0SINstVPrszWYx9uJEYcqChxvGM1RIoR8uZe-sw", sheet = "donor_groups"))
n_groups = nrow(donorGroups)
n_PSA_iter = nrow(ZIKV_psa_nmc)

opt_by_group_percent <- data.table(
  group = donorGroups$donorGroup,
  year = donorGroups$year,
  z = 0,
  a1 = 0,
  a2 = 0,
  a3 = 0,
  a4 = 0,
  m1 = 0)

opt_by_iter_metrics_2017 <- data.table(
  iter = 1:10000,
  obj_cost = 0, 
  yield = 0, 
  test_cost = 0, 
  mod_cost = 0, 
  resid_risk1 = 0,
  resid_risk2 = 0,
  downsteam_NMC = 0
  )
opt_by_iter_metrics_2018 <- data.table(
  iter = 1:10000,
  obj_cost = 0, 
  yield = 0, 
  test_cost = 0, 
  mod_cost = 0, 
  resid_risk1 = 0,
  resid_risk2 = 0,
  downsteam_NMC = 0
)
opt_by_iter_metrics_2019 <- data.table(
  iter = 1:10000,
  obj_cost = 0, 
  yield = 0, 
  test_cost = 0, 
  mod_cost = 0, 
  resid_risk1 = 0,
  resid_risk2 = 0,
  downsteam_NMC = 0
)

#Variable names
vars = unique(params_PSA$key)

Sys.time()
for (PSA_iter in 1:2){#10000){
  params = params_PSA[, c(1,2,3, 3+PSA_iter), with = FALSE]
  
  # Translate variable sheet into matrices for analysis
  
  for (var in vars){
    if(is.na(params[key == var]$row_idx[1])){
      eval(call("<-", as.name(var), params[key == var, 4][[1]]))
    } else {
      eval(call("<-", as.name(var), matrix(t(params[key == var, 4]), 
                                           ncol = params[key == var, max(col_idx)])))
    }
  }
  
  c_k = matrix(c(ZIKV_psa_nmc[PSA_iter, NMC_per_donation], WNV_psa_nmc[PSA_iter, NMC_per_donation]))
  
  # Count params
  K = length(c_k) #number of TTIDs
  I = 1 #number of donor segments, but performing analysis separately for each jurisdiction to take advantage of linear separability
  J = length(R_jk[,1]) #tests
  L = length(c_mod_l) #Mods
  
  
  # Enumerate policies, using linear separability
  n_pols = 2^(L+J) + 1
  pol_length = I*(1+J+L)
  
  #create policies matrix with every feasible policy for a single donor group
  policies <- matrix(1, n_pols, pol_length)
  policies[1, 1:pol_length] <- rep(0, pol_length)
  policies[2:n_pols, 2:pol_length] <- as.matrix(expand.grid(replicate(pol_length - 1, 0:1, simplify=FALSE)))
  
  #Create groupPolicies datatable which has exhaustive list of all policies for each donor group
  groupPolicies <- data.table(
    group = rep(donorGroups$donorGroup, times = 1, each = n_pols),
    year = rep(donorGroups$year, times = 1, each = n_pols),
    p1 = rep(donorGroups$pZIKV, times = 1, each = n_pols),
    p2 = rep(donorGroups$pWNV, times = 1, each = n_pols),
    n = rep(donorGroups$numDonor, times = 1, each = n_pols),
    z = rep(policies[ , 1], times = n_groups),
    a1 = rep(policies[ , 2], times = n_groups),
    a2 = rep(policies[ , 3], times = n_groups),
    a3 = rep(policies[ , 4], times = n_groups),
    a4 = rep(policies[ , 5], times = n_groups),
    m1 = rep(policies[ , 6], times = n_groups),
    obj_cost = 0, 
    yield = 0, 
    test_cost = 0, 
    mod_cost = 0, 
    resid_risk1 = 0,
    resid_risk2 = 0,
    downsteam_NMC = 0
  )
  
  # Evaluate each portfolio for each donor group
  for (row in 1:nrow(groupPolicies)){
    groupPolicies[ row, 12:18] =  as.list(cost_portfolio_unopt(z_i = matrix(unlist(groupPolicies[row, "z"])),
                                                               M_li = matrix(unlist(groupPolicies[row, "m1"], I, L)),
                                                               A_ji = matrix(unlist(groupPolicies[row, c("a1", "a2", "a3", "a4")], I, J)),
                                                               c_k = c_k, 
                                                               P_ik =  t(matrix(unlist(groupPolicies[row, c("p1", "p2")], I, K))),
                                                               d_i = d, 
                                                               w_i = w, 
                                                               n_i = matrix(unlist(groupPolicies[row, "n"])),
                                                               c_test_j = c_test_j, 
                                                               R_jk = R_jk, 
                                                               Q_jk = Q_jk, 
                                                               g = g, 
                                                               c_mod_l = c_mod_l, 
                                                               H_lk = H_lk, 
                                                               full_output = 1))
  }
  #Get optimal polict by donor group
  opt_by_group <- groupPolicies[ , .SD[which.min(obj_cost)], by = group]
  
  opt_by_group_percent[, c(3:8)] = opt_by_group_percent[, c(3:8)] + opt_by_group[, 6:11]
    opt_by_iter_metrics_2017[PSA_iter, c(2:8)] = as.list(colSums(opt_by_group[year=="2017", 12:18]))
  opt_by_iter_metrics_2018[PSA_iter, c(2:8)] = as.list(colSums(opt_by_group[year=="2018", 12:18]))
  opt_by_iter_metrics_2019[PSA_iter, c(2:8)] = as.list(colSums(opt_by_group[year=="2019", 12:18]))
  
}
Sys.time()
opt_by_group_percent[, c(3:8)] = opt_by_group_percent[, c(3:8)]/855

fwrite(opt_by_group_percent, 'Results/opt_by_group_percent_first_855.csv')

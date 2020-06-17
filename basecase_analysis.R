library(data.table)
library(googlesheets4)

# setwd("G:/My Drive/Blood Transfusion/Optimal portfolio/blood portfolio r project")

####
#
# Basecase
#
####

# Get parameters from Google shet
params <- data.table(read_sheet("1Yjfq0SINstVPrszWYx9uJEYcqChxvGM1RIoR8uZe-sw", sheet = "params"))[, c("key", "row_idx", "col_idx", "Basecase")]

#donorGroups<- data.table(read_sheet("1Yjfq0SINstVPrszWYx9uJEYcqChxvGM1RIoR8uZe-sw", sheet = "donor_groups"))
donorGroups <- fread('data/donor_groups_zip3.csv')
n_groups = nrow(donorGroups)

# Translate variable sheet into matrices for analysis
vars = unique(params$key)
for (var in vars){
  if(is.na(params[key == var]$row_idx[1])){
    eval(call("<-", as.name(var), params[key == var, Basecase]))
  } else {
    eval(call("<-", as.name(var), matrix(params[key == var, Basecase], ncol = params[key == var, max(col_idx)])))
  }
}


K = length(c_k) #number of TTIDs
I = n_groups #number of donor segments
J = length(R_jk[,1]) #tests
L = length(c_mod_l) #Mods


### No intervention and test all strategies
source('portfolio_functions.R')


basecase_pol_compare <- data.table(
  year =  numeric(),
  policy = character(),
  obj_cost = numeric(),
  yield = numeric(),
  test_cost = numeric(),
  mod_cost = numeric(),
  downsteam_NMC = numeric(),
  resid_risk1 = numeric(),
  resid_risk2 =numeric()
)



I = n_groups/3
for(yr in unique(donorGroups$year)){
  val = cost_portfolio_unopt(z_i = matrix(1, I),
                             M_li = matrix(0, L, I),
                             A_ji = matrix(0, J, I),
                             c_k = c_k, 
                             P_ik =  as.matrix(donorGroups[year == yr, c("pZIKV", "pWNV")]),
                             d_i = matrix(d, I), 
                             w_i = matrix(w, I), 
                             n_i = as.matrix(donorGroups[year == yr, "numDonor"]),
                             c_test_j = c_test_j, 
                             R_jk = R_jk, 
                             Q_jk = Q_jk, 
                             g = g, 
                             c_mod_l = c_mod_l, 
                             H_lk = H_lk, 
                             full_output = 1)
  basecase_pol_compare <- rbind(basecase_pol_compare,
                   t(data.table(c(yr, "No intervention", unlist(val)))), use.names = FALSE)
  
  val = cost_portfolio_unopt(z_i = matrix(1, I),
                             M_li = matrix(0, L, I),
                             A_ji = matrix(c(0,1), J, I),
                             c_k = c_k, 
                             P_ik =  as.matrix(donorGroups[year == yr, c("pZIKV", "pWNV")]),
                             d_i = matrix(d, I), 
                             w_i = matrix(w, I), 
                             n_i = as.matrix(donorGroups[year == yr, "numDonor"]),
                             c_test_j = c_test_j, 
                             R_jk = R_jk, 
                             Q_jk = Q_jk, 
                             g = g, 
                             c_mod_l = c_mod_l, 
                             H_lk = H_lk, 
                             full_output = 1)
  
  basecase_pol_compare <- rbind(basecase_pol_compare,
                   t(data.table(c(yr, "Universal MP-NAT", unlist(val)))), use.names = FALSE)
  
}

fwrite(basecase_pol_compare, "results/basecase_policy_comparison_zip3.csv")



# Enumerate policies -----
I = 1 #number of donor segments, but performing analysis separately for each jurisdiction to take advantage of linear separability
# First, enumerate all policies with between 2 and I - 1 donor segments
n_pols = 2^(L+J) + 1
pol_length = I*(1+J+L)

# expand.grid(replicate(pol_length - 1, 0:1, simplify=FALSE))


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
  downsteam_NMC = 0,
  resid_risk1 = 0,
  resid_risk2 = 0
)



#iterate over rows and calculate outcome of each policy
Sys.time()
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

Sys.time()

#Get optimal polict by donor group
opt_by_group <- groupPolicies[ , .SD[which.min(obj_cost)], by = group]
opt_by_group[ , season := substr(group, 13, 13)]
fwrite(opt_by_group, "Results/opt_by_zip3_basecase.csv")


opt_by_group[z == 0 | test_cost + mod_cost > 0]

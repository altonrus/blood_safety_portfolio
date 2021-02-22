library(data.table)

Process_PSA <- function(iter_start, iter_end, params_PSA, ZIKV_psa_nmc, WNV_psa_nmc, donorGroups){
  n_groups = nrow(donorGroups)
  
  opt_by_group_percent <- data.table(
    group = donorGroups$donorGroup,
    year = donorGroups$year,
    z = 0,
    a1 = 0,
    a2 = 0,
    a3 = 0,
    a4 = 0,
    m1 = 0)
  
  setorder(opt_by_group_percent, group)
  
  metrics_opt <- data.table(
    iter = rep(iter_start:iter_end, each = 3),
    policy = "Optimal",
    year = 0,
    obj_cost = 0,
    yield = 0,
    test_cost = 0,
    mod_cost = 0,
    downsteam_NMC = 0,
    resid_risk1 = 0,
    resid_risk2 = 0
  )
  
  metrics_no_int <- data.table(
    iter = rep(iter_start:iter_end, each = 3),
    policy = "No intervention",
    year = 0,
    obj_cost = 0,
    yield = 0,
    test_cost = 0,
    mod_cost = 0,
    downsteam_NMC = 0,
    resid_risk1 = 0,
    resid_risk2 = 0
  )
  
  metrics_univ_MP <- data.table(
    iter = rep(iter_start:iter_end, each = 3),
    policy = "Universal MP-NAT",
    year = 0,
    obj_cost = 0,
    yield = 0,
    test_cost = 0,
    mod_cost = 0,
    downsteam_NMC = 0,
    resid_risk1 = 0,
    resid_risk2 = 0
  )
  
  
  #Variable names
  vars = unique(params_PSA$key)
  
  for (PSA_iter in iter_start:iter_end){
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
    
    donorGroups_with_risk <- donorGroups[pZIKV > 0 | pWNV > 0]
    donorGroups_no_risk <- donorGroups[!(pZIKV > 0 | pWNV > 0)]
    
    #Create groupPolicies datatable which has exhaustive list of all policies for each donor group
    n_groups = nrow(donorGroups_with_risk)
    groupPolicies <- data.table(
      group = rep(donorGroups_with_risk$donorGroup, times = 1, each = n_pols),
      year = rep(donorGroups_with_risk$year, times = 1, each = n_pols),
      p1 = rep(donorGroups_with_risk$pZIKV, times = 1, each = n_pols),
      p2 = rep(donorGroups_with_risk$pWNV, times = 1, each = n_pols),
      n = rep(donorGroups_with_risk$numDonor, times = 1, each = n_pols),
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
    
    n_groups = nrow(donorGroups_no_risk)
    groupPolicies <- rbind(groupPolicies,
                           data.table(
                             group = rep(donorGroups_no_risk$donorGroup, times = 1, each = 2),
                             year = rep(donorGroups_no_risk$year, times = 1, each = 2),
                             p1 = rep(donorGroups_no_risk$pZIKV, times = 1, each = 2),
                             p2 = rep(donorGroups_no_risk$pWNV, times = 1, each = 2),
                             n = rep(donorGroups_no_risk$numDonor, times = 1, each = 2),
                             z = rep(c(1,1), times = n_groups),
                             a1 = rep(c(0, 0), times = n_groups),
                             a2 = rep(c(0, 1), times = n_groups),
                             a3 = rep(c(0, 0), times = n_groups),
                             a4 = rep(c(0, 1), times = n_groups),
                             m1 = rep(c(0, 0), times = n_groups),
                             obj_cost = 0, 
                             yield = 0, 
                             test_cost = 0, 
                             mod_cost = 0, 
                             downsteam_NMC = 0,
                             resid_risk1 = 0,
                             resid_risk2 = 0
    ))
    
    setorder(groupPolicies, group)
    
    # Delete rows for which p1 == 0 & p2 == 0 &
    
    
    # Evaluate each portfolio for each donor group
    for (row in 1:nrow(groupPolicies)){
      groupPolicies[row, 12:18] =  as.list(cost_portfolio_main(z_i = matrix(unlist(groupPolicies[row, "z"])),
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
    
    #For most outcomes a simple sum works, but for not residual risk
    #opt_by_group[ , sum(.SD), .SD = 12:16 , by = year]
    
    
    metrics_opt[iter == PSA_iter, c(3:10)] <- cbind(opt_by_group[ , lapply(.SD, sum), .SD = 12:16 , by = year], 
                                                    opt_by_group[ , sum(resid_risk1*yield)/sum(yield), by = year][, 2],
                                                    opt_by_group[ , sum(resid_risk2*yield)/sum(yield), by = year][, 2])
    
    univ_MP <- groupPolicies[z == 1 & a1 == 0 & a2 == 1 & a3 == 0 & a4 == 1 & m1 == 0]
    
    metrics_univ_MP[iter == PSA_iter, c(3:10)] <- cbind(univ_MP[ , lapply(.SD, sum), .SD = 12:16 , by = year], 
                                                        univ_MP[ , sum(resid_risk1*yield)/sum(yield), by = year][, 2],
                                                        univ_MP[ , sum(resid_risk2*yield)/sum(yield), by = year][, 2])
 
    no_int <- groupPolicies[z == 1 & a1 == 0 & a2 == 0 & a3 == 0 & a4 == 0 & m1 == 0]
    
    metrics_no_int[iter == PSA_iter, c(3:10)] <- cbind(no_int[ , lapply(.SD, sum), .SD = 12:16 , by = year], 
                                                        no_int[ , sum(resid_risk1*yield)/sum(yield), by = year][, 2],
                                                        no_int[ , sum(resid_risk2*yield)/sum(yield), by = year][, 2])
    
    PSA_metrics <- rbind(metrics_opt, metrics_no_int, metrics_univ_MP)
    
  }
  return(list(opt_by_group_percent, PSA_metrics))
}


Process_and_save_PSA <- function(iter_start, iter_end, params_PSA, ZIKV_psa_nmc, WNV_psa_nmc, donorGroup,
                                 fname_prefix = "PSA_zip3_"){
  output <- Process_PSA(iter_start, iter_end, params_PSA, ZIKV_psa_nmc, WNV_psa_nmc, donorGroups)
  
  fwrite(output[[1]], paste0(fname_prefix, 'opt_by_group_', iter_start, 'to', iter_end, '.csv'))
  fwrite(output[[2]], paste0(fname_prefix, 'metrics_', iter_start, 'to', iter_end, '.csv'))
  
}





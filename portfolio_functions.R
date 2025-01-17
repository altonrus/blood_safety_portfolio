cost_portfolio_main <- function(z_i, M_li, A_ji, c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_l, H_lk, 
                                full_output = 0 #If set to 1, returns 5 other metrics along with objective cost value
                                ) {
  K = length(c_k) #number of TTIDs
  I = length(z_i) #number of donor segments, but performing analysis separately for each jurisdiction to take advantage of linear separability
  J = length(R_jk[,1]) #tests
  L = length(c_mod_l) #Mods
  
  
  B1 <- matrix(1L, I, K) #Will vectorize
  for (i in 1:I){
    for (k in 1:K){
      for (j in 1:J){
        B1[i,k] = B1[i,k]*(1-R_jk[j,k]*A_ji[j,i])
      }
    }
  }
  B2 <- matrix(1L, I, K) #Will vectorize
  for (i in 1:I){
    for (k in 1:K){
      for (l in 1:L){
        B2[i,k] = B2[i,k]*(1-M_li[l,i]*(1-H_lk[l,k]))
      }
    }
  }
  
  v2 <- matrix(1L, I) #Will vectorize
  for (i in 1:I){
    all_neg = 1
    for (k in 1:K){
      for (j in 1:J){
        all_neg = all_neg*((1-P_ik[i,k]) * (1+A_ji[j,i]*(Q_jk[j,k] - 1)) + P_ik[i,k] * (1 - A_ji[j,i]*R_jk[j,k]))
      }
    }
    v2[i] = 1 - all_neg
  }
  if (full_output == 0) {
    return(drop(t((1-z_i)*n_i) %*% d_i + t(z_i * n_i) %*% (w_i + t(t(c_mod_l) %*% M_li) + t(t(c_test_j) %*% A_ji) + t(t(c_mod_l) %*% M_li) + 
                                                             (1-v2) * (B1*B2*P_ik) %*% c_k + ((g+d_i)*v2))))
  } else {
    colnames(z_i) <- NULL
    accepted_donors_by_seg = z_i * n_i # I X 1 vector
    RR_matrix = (B1*B2*P_ik) # I X K matrix

    
    yield = drop(t(z_i * (1 - v2) ) %*% n_i) # single variable
    false_negs = t(RR_matrix) %*% (z_i * (1 - v2) * n_i)  # K X 1 vector
    avg_resid_risk = false_negs / yield # K X 1 vector
    test_cost = drop(t(accepted_donors_by_seg) %*% (t(A_ji) %*% c_test_j)) # single variable
    mod_cost = drop(t(accepted_donors_by_seg) %*% (t(M_li) %*% c_mod_l)) # single variable
    defer_cost = drop(t((1-z_i)*n_i) %*% d_i)
    obj_cost = defer_cost + drop(t(accepted_donors_by_seg) %*% (w_i + t(t(c_mod_l) %*% M_li) + t(t(c_test_j) %*% A_ji) + (1-v2) * ((B1*B2*P_ik)%*%c_k) + (g+d_i)*v2))
    downstream_net_mon_cost =  drop(t(accepted_donors_by_seg) %*% (B1*B2*P_ik)%*%c_k)
    
    return(c("obj_cost" = obj_cost, 
             "yield" = yield, 
             "test_cost" = test_cost, 
             "mod_cost" = mod_cost, 
             "downstream_net_mon_cost" = downstream_net_mon_cost, 
             "avg_resid_risk" = avg_resid_risk))
  }
}


cost_portfolio_linearP <- function(z_i, M_li, A_ji, c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, c_mod_l, H_lk) {
  K = length(c_k) #number of TTIDs
  I = length(z_i) #number of donor segments, but performing analysis separately for each jurisdiction to take advantage of linear separability
  J = length(R_jk[,1]) #tests
  L = length(c_mod_l) #Mods
  
  
  B1 <- matrix(1L, I, K) #Will vectorize
  for (i in 1:I){
    for (k in 1:K){
      for (j in 1:J){
        B1[i,k] = B1[i,k]*(1-R_jk[j,k]*A_ji[j,i])
      }
    }
  }
  B2 <- matrix(1L, I, K) #Will vectorize
  for (i in 1:I){
    for (k in 1:K){
      for (l in 1:L){
        B2[i,k] = B2[i,k]*(1-M_li[l,i]*(1-H_lk[l,k]))
      }
    }
  }

  return(drop(t((1-z_i)*n_i) %*% d_i + t(z_i * n_i) %*% (w_i + t(t(c_mod_l) %*% M_li) + t(t(c_test_j) %*% A_ji) + t(t(c_mod_l) %*% M_li) + (B1*B2*P_ik)%*%c_k)))

}


cost_portfolio_nodeferral <- function(M_li, A_ji, c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_l, H_lk) {
  B1 <- matrix(1L, I, K) #Will vectorize
  for (i in 1:I){
    for (k in 1:K){
      for (j in 1:J){
        B1[i,k] = B1[i,k]*(1-R_jk[j,k]*A_ji[j,i])
      }
    }
  }
  B2 <- matrix(1L, I, K) #Will vectorize
  for (i in 1:I){
    for (k in 1:K){
      for (n in 1:N){
        B2[i,k] = B2[i,k]*(1-M_li[n,i]*(1-H_lk[n,k]))
      }
    }
  }
  
  v2 <- matrix(1L, I) #Will vectorize
  for (i in 1:I){
    all_neg = 1
    for (k in 1:K){
      for (j in 1:J){
        all_neg = all_neg*((1-P_ik[i,k]) * (1+A_ji[j,i]*(Q_jk[j,k] - 1)) + P_ik[i,k] * (1 - A_ji[j,i]*R_jk[j,k]))
      }
    }
    v2[i] = 1 - all_neg
  }
  return(drop(t( n_i) %*% (w_i + t(t(c_mod_l) %*% M_li) + t(t(c_test_j) %*% A_ji) + t(t(c_mod_l) %*% M_li) + (B1*B2*P_ik)%*%c_k + ((g+d_i)*v2))))
}






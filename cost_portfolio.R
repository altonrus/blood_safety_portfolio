library(Rcpp)
library(RcppArmadillo)
setwd("G:/My Drive/Blood Transfusion/Optimal portfolio/blood portfolio r project")
sourceCpp("cost_portfolio.cpp")

#OBJECTIVE FUNCTION Aman -----------------
# This function does not give the same results. Need to correct before implementing
cost_portfolio_aman <- function(z_i, M_ni, A_ji, c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_n, H_nk) {
  B1 <- mat_prod_prod(A_ji, R_jk)
  B2 <- mat_prod_prod(M_ni, H_nk - 1)
  vec_pos_test <- get_vec_pos_test(P_ik, A_ji, R_jk, Q_jk)

  c_mod_n_by_M_ni <- t(crossprod(c_mod_n, M_ni))
  right <-   
    w_i + 
    c_mod_n_by_M_ni + 
    t(crossprod(c_test_j, A_ji)) + 
    c_mod_n_by_M_ni + 
    (B1 * B2 * P_ik) %*% c_k + 
    (g * vec_pos_test)
  c(crossprod((1-z_i) * n_i, d_i) + crossprod(z_i * n_i, right))
}

#OBJECTIVE FUNCTION unoptimized-----------------
cost_portfolio_unopt <- function(z_i, M_ni, A_ji, c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_n, H_nk, full_output = 0) {
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
        B2[i,k] = B2[i,k]*(1-M_ni[n,i]*(1-H_nk[n,k]))
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
    return(drop(t((1-z_i)*n_i) %*% d_i + t(z_i * n_i) %*% (w_i + t(t(c_mod_n) %*% M_ni) + t(t(c_test_j) %*% A_ji) + t(t(c_mod_n) %*% M_ni) + (B1*B2*P_ik)%*%c_k + (g*v2))))
  } else {
    accepted_donors_by_seg = z_i * n_i # I X 1 vector
    RR_matrix = (B1*B2*P_ik) # I X K matrix
    yield = drop(t(z_i * all_neg) %*% n_i) # single variable
    false_negs = t(RR_matrix) %*% (z_i * all_neg * n_i)  # K X 1 vector
    avg_resid_risk = false_negs / yield # K X 1 vector
    test_cost = drop(t(accepted_donors_by_seg) %*% (t(A_ji) %*% c_test_j)) # single variable
    mod_cost = drop(t(accepted_donors_by_seg) %*% (t(M_ni) %*% c_mod_n)) # single variable
    defer_cost = drop(t((1-z_i)*n_i) %*% d_i)
    obj_cost = defer_cost + drop(t(accepted_donors_by_seg) %*% (w_i + t(t(c_mod_n) %*% M_ni) + t(t(c_test_j) %*% A_ji) + t(t(c_mod_n) %*% M_ni) + (B1*B2*P_ik)%*%c_k + (g*v2)))
    
    return(c(obj_cost, yield, test_cost, mod_cost, avg_resid_risk))
  }
}

# library(microbenchmark)
# microbenchmark(new = cost_portfolio(z_i, M_ni, A_ji, c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_n, H_nk), 
#                old = cost_portfolio_old(z_i, M_ni, A_ji, c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_n, H_nk),
#                times = 1000)

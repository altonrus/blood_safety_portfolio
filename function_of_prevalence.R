
# Optimal policy as function of p_1 (Zika prevalence) and p_2 (WNV prevalence)
library(data.table)
library(googlesheets4)
library(ggplot2)
params <- data.table(read_sheet("1Yjfq0SINstVPrszWYx9uJEYcqChxvGM1RIoR8uZe-sw", sheet = "params"))[, c("key", "row_idx", "col_idx", "Basecase")]
theme_set(theme_bw())

# Translate variable sheet into matrices for analysis
vars = unique(params[nchar(key) < 10 , key])
for (var in vars){
  if(is.na(params[key == var]$row_idx[1])){
    eval(call("<-", as.name(var), params[key == var, Basecase]))
  } else {
    eval(call("<-", as.name(var), matrix(params[key == var, Basecase], ncol = params[key == var, max(col_idx)])))
  }
}



#Create policies
n_pols = 2^(L+J) + 1
pol_length = (1+J+L)
policies <- matrix(1L, n_pols, pol_length)
policies[1, 1:pol_length] <- rep(0, pol_length)
policies[2:n_pols, 2:pol_length] <- as.matrix(expand.grid(replicate(pol_length - 1, 0:1, simplify=FALSE)))

constants_linear <- function(pi,
                       d, 
                       w, 
                       c_test_j,
                       R_jk, 
                       Q_jk, 
                       g, 
                       c_mod_l, 
                       H_lk){
  
  J = length(c_test_j)
  K = length(c_k)
  L = length(c_mod_l)
  
  z = pi[1]
  a_j = pi[2:(1+J)]
  m_l = pi[(2+J):length(pi)]
  
  alpha = (1-z)*d + z*(w+ a_j %*% c_test_j + m_l %*% c_mod_l )
  
  v_test = matrix(1L, K)
  for (j in 1:J){
    v_test = v_test*(1 - R_jk[j, ]*a_j[j])
  }
  v_mod = matrix(1L, K)
  for(l in 1:L){
    v_mod = v_mod*(1 + m_l[1]*(H_lk[l, ]-1))
  }
  
  beta <- z*v_test*v_mod*c_k
  
  return(list(alpha = alpha, beta = beta))
}



a_j = c(1, 1, 1, 1)
m_l = c(1)
z = 1



m = 1
b = 0




ggplot() + geom_abline(slope = m, intercept = b) +ylim(c(0,1)) + xlim(c(0,1))








get_linear_system <- function(policies,
                              c_k, d, w, c_test_j,
                              R_jk, Q_jk, g,
                              c_mod_l, H_lk){
  
  pairs <- combn(1:nrow(policies), 2, simplify = FALSE)
  coef_p1_vec <- c()
  coef_p2_vec <- c()
  const_vec <- c()
  for(pair in pairs){
    pair <- unlist(pair)
    
    param_1 <- constants_linear(pi = policies[pair[1], ],
                                d, 
                                w, 
                                c_test_j,
                                R_jk, 
                                Q_jk, 
                                g, 
                                c_mod_l, 
                                H_lk)
    
    param_2 <- constants_linear(pi = policies[pair[2], ],
                                d, 
                                w, 
                                c_test_j,
                                R_jk, 
                                Q_jk, 
                                g, 
                                c_mod_l, 
                                H_lk)
    
    # m = ((param_2$alpha - param_1$alpha)/(param_1$beta[1] - param_2$beta[1]))[1]
    # b = ((param_2$beta[2] - param_1$beta[2])/(param_1$beta[1] - param_2$beta[1]))[1]
    coef_p1 = (param_1$beta[1] - param_2$beta[1])[1]
    coef_p2 = (param_1$beta[2] - param_2$beta[2])[1]
    const = (param_2$alpha - param_1$alpha)[1]
    
    if ( (const/coef_p1 >= 0) | (const/coef_p2 >= 0)){
      coef_p1_vec <- append(coef_p1_vec, coef_p1)
      coef_p2_vec <- append(coef_p2_vec, coef_p2)
      const_vec <- append(const_vec, const)
    }
  }
  return(linear_system_params = data.table(coef_p1 = coef_p1_vec, coef_p2 = coef_p2_vec, const = const_vec))
  
}


linear_system <- get_linear_system(policies, 
                  c_k, d, w, c_test_j,
                  R_jk, Q_jk, g,
                  c_mod_l, H_lk)

linear_system[ , p1_1 := const/coef_p1]
linear_system[ , p2_1 := 0]
linear_system[ , p1_2 := 0]
linear_system[ , p2_2 := const/coef_p2]
linear_system[ , ID := 1:nrow(linear_system)]

#Melting for ploting with ggplot using 2 points
dt <- cbind(melt(linear_system, id.vars = c("ID"), measure.vars = c("p1_1", "p1_2"), value.name = "p1"),
      p2 = melt(linear_system, id.vars = c("ID"), measure.vars = c("p2_1", "p2_2"), value.name = "p2")$p2)

#This plot doesn't do vertical and horizontals properly
ggplot()+
  geom_line(data = dt, aes(x = p1, y = p2, colour = factor(ID)))+
  theme(legend.position = "none")+
  scale_x_log10(name = "Zika prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_y_log10(name = "WNV prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
                expand = c(0, 0),
                minor_breaks = NULL)+
  coord_cartesian(xlim = c(1e-6, 1), ylim = c(1e-6, 1))

# Non-scaled axes
ggplot()+
  geom_line(data = dt, aes(x = p1, y = p2, colour = factor(ID)))+
  theme(legend.position = "none")+
  scale_x_continuous(name = "Zika prevalence", 
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_y_continuous(name = "WNV prevalence", 
                expand = c(0, 0),
                minor_breaks = NULL)+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))


yac_assign(alpha_1, (1-z)*d + z*(w+ a_j %*% c_test_j + m_l %*% c_mod_l ))



yac_str("Sum(j, 1, 4, (1-p1)(1+a(q-1)))")

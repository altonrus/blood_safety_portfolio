library(data.table)





#setwd("G:/My Drive/Blood Transfusion/Optimal portfolio/blood portfolio r project")

library(googlesheets4)
params <- data.table(read_sheet("1Yjfq0SINstVPrszWYx9uJEYcqChxvGM1RIoR8uZe-sw", sheet = "params"))[, c("key", "row_idx", "col_idx", "Basecase")]

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
I = 1 #number of donor segments
J = length(R_jk[,1]) #tests
L = length(c_mod_l) #Mods

cost_zikv = 1574
cost_wnv = 6594.910244
QALY_zikv = 0.008165121
QALY_wnv = 0.031396707

source('portfolio_functions.R')


# Enumerate policies -----
# First, enumerate all policies with between 2 and I - 1 donor segments
n_pols = 2^(L+J) + 1
pol_length = I*(1+J+L)

# expand.grid(replicate(pol_length - 1, 0:1, simplify=FALSE))


#create policies matrix with every feasible policy for a single donor group
policies <- matrix(1, n_pols, pol_length)
policies[1, 1:pol_length] <- rep(0, pol_length)
policies[2:n_pols, 2:pol_length] <- as.matrix(expand.grid(replicate(pol_length - 1, 0:1, simplify=FALSE)))




# 
# ##########################################
# WTP = 1e6
# c_k <- matrix(data = c(cost_zikv + QALY_zikv*WTP, cost_wnv + QALY_wnv*WTP),
#               ncol = 1)
# dt.prev <- data.table(expand.grid(exp(seq(log(1e-6), log(0.3), length.out = 100)),
#                                   exp(seq(log(1e-6), log(0.3), length.out = 100))))
# 
# colnames(dt.prev) <- c("p_ZIKV", "p_WNV")
# 
# 
# #Create groupPolicies datatable which has exhaustive list of all policies for each donor group
# groupPolicies <- data.table(
#   p1 = rep(dt.prev$p_ZIKV, each = n_pols),
#   p2 = rep(dt.prev$p_WNV, each = n_pols),
#   n = 1,
#   z = rep(policies[ , 1], times = nrow(dt.prev)),
#   a1 = rep(policies[ , 2], times = nrow(dt.prev)),
#   a2 = rep(policies[ , 3], times = nrow(dt.prev)),
#   a3 = rep(policies[ , 4], times = nrow(dt.prev)),
#   a4 = rep(policies[ , 5], times = nrow(dt.prev)),
#   m1 = rep(policies[ , 6], times = nrow(dt.prev)),
#   obj_cost = 0, 
#   yield = 0, 
#   test_cost = 0, 
#   mod_cost = 0,
#   downsteam_NMC = 0,
#   resid_risk1 = 0,
#   resid_risk2 = 0
# )
# 
# 
# 
# for (row in 1:nrow(groupPolicies)){
#   groupPolicies[ row, 10:16] =  as.list(cost_portfolio_unopt(z_i = matrix(unlist(groupPolicies[row, "z"])),
#                                                              M_li = matrix(unlist(groupPolicies[row, "m1"], I, L)),
#                                                              A_ji = matrix(unlist(groupPolicies[row, c("a1", "a2", "a3", "a4")], I, J)),
#                                                              c_k = c_k, 
#                                                              P_ik =  t(matrix(unlist(groupPolicies[row, c("p1", "p2")], I, K))),
#                                                              d_i = d, 
#                                                              w_i = w, 
#                                                              n_i = matrix(unlist(groupPolicies[row, "n"])),
#                                                              c_test_j = c_test_j, 
#                                                              R_jk = R_jk, 
#                                                              Q_jk = Q_jk, 
#                                                              g = g, 
#                                                              c_mod_l = c_mod_l, 
#                                                              H_lk = H_lk, 
#                                                              full_output = 1))
# }
# 
# dt <- groupPolicies[ , .SD[which.min(obj_cost)], by = c("p1", "p2")]
# 
# 
# fwrite(dt, "Results//prev_2way_wtp1e6.csv")
# ###############################################
# 
# 
# WTP = 1e5
# c_k <- matrix(data = c(cost_zikv + QALY_zikv*WTP, cost_wnv + QALY_wnv*WTP),
#               ncol = 1)
# dt.prev <- data.table(expand.grid(exp(seq(log(1e-6), log(0.3), length.out = 100)),
#                                   exp(seq(log(1e-6), log(0.3), length.out = 100))))
# 
# colnames(dt.prev) <- c("p_ZIKV", "p_WNV")
# 
# 
# #Create groupPolicies datatable which has exhaustive list of all policies for each donor group
# groupPolicies <- data.table(
#   p1 = rep(dt.prev$p_ZIKV, each = n_pols),
#   p2 = rep(dt.prev$p_WNV, each = n_pols),
#   n = 1,
#   z = rep(policies[ , 1], times = nrow(dt.prev)),
#   a1 = rep(policies[ , 2], times = nrow(dt.prev)),
#   a2 = rep(policies[ , 3], times = nrow(dt.prev)),
#   a3 = rep(policies[ , 4], times = nrow(dt.prev)),
#   a4 = rep(policies[ , 5], times = nrow(dt.prev)),
#   m1 = rep(policies[ , 6], times = nrow(dt.prev)),
#   obj_cost = 0, 
#   yield = 0, 
#   test_cost = 0, 
#   mod_cost = 0,
#   downsteam_NMC = 0,
#   resid_risk1 = 0,
#   resid_risk2 = 0
# )
# 
# 
# 
# for (row in 1:nrow(groupPolicies)){
#   groupPolicies[ row, 10:16] =  as.list(cost_portfolio_unopt(z_i = matrix(unlist(groupPolicies[row, "z"])),
#                                                              M_li = matrix(unlist(groupPolicies[row, "m1"], I, L)),
#                                                              A_ji = matrix(unlist(groupPolicies[row, c("a1", "a2", "a3", "a4")], I, J)),
#                                                              c_k = c_k, 
#                                                              P_ik =  t(matrix(unlist(groupPolicies[row, c("p1", "p2")], I, K))),
#                                                              d_i = d, 
#                                                              w_i = w, 
#                                                              n_i = matrix(unlist(groupPolicies[row, "n"])),
#                                                              c_test_j = c_test_j, 
#                                                              R_jk = R_jk, 
#                                                              Q_jk = Q_jk, 
#                                                              g = g, 
#                                                              c_mod_l = c_mod_l, 
#                                                              H_lk = H_lk, 
#                                                              full_output = 1))
# }
# 
# dt <- groupPolicies[ , .SD[which.min(obj_cost)], by = c("p1", "p2")]
# 
# fwrite(dt, "Results//prev_2way_wtp1e5.csv")
# ###############################################
# 
# 
# WTP = 1e7
# c_k <- matrix(data = c(cost_zikv + QALY_zikv*WTP, cost_wnv + QALY_wnv*WTP),
#               ncol = 1)
# dt.prev <- data.table(expand.grid(exp(seq(log(1e-6), log(0.3), length.out = 100)),
#                                   exp(seq(log(1e-6), log(0.3), length.out = 100))))
# 
# colnames(dt.prev) <- c("p_ZIKV", "p_WNV")
# 
# 
# #Create groupPolicies datatable which has exhaustive list of all policies for each donor group
# groupPolicies <- data.table(
#   p1 = rep(dt.prev$p_ZIKV, each = n_pols),
#   p2 = rep(dt.prev$p_WNV, each = n_pols),
#   n = 1,
#   z = rep(policies[ , 1], times = nrow(dt.prev)),
#   a1 = rep(policies[ , 2], times = nrow(dt.prev)),
#   a2 = rep(policies[ , 3], times = nrow(dt.prev)),
#   a3 = rep(policies[ , 4], times = nrow(dt.prev)),
#   a4 = rep(policies[ , 5], times = nrow(dt.prev)),
#   m1 = rep(policies[ , 6], times = nrow(dt.prev)),
#   obj_cost = 0, 
#   yield = 0, 
#   test_cost = 0, 
#   mod_cost = 0,
#   downsteam_NMC = 0,
#   resid_risk1 = 0,
#   resid_risk2 = 0
# )
# 
# 
# 
# for (row in 1:nrow(groupPolicies)){
#   groupPolicies[ row, 10:16] =  as.list(cost_portfolio_unopt(z_i = matrix(unlist(groupPolicies[row, "z"])),
#                                                              M_li = matrix(unlist(groupPolicies[row, "m1"], I, L)),
#                                                              A_ji = matrix(unlist(groupPolicies[row, c("a1", "a2", "a3", "a4")], I, J)),
#                                                              c_k = c_k, 
#                                                              P_ik =  t(matrix(unlist(groupPolicies[row, c("p1", "p2")], I, K))),
#                                                              d_i = d, 
#                                                              w_i = w, 
#                                                              n_i = matrix(unlist(groupPolicies[row, "n"])),
#                                                              c_test_j = c_test_j, 
#                                                              R_jk = R_jk, 
#                                                              Q_jk = Q_jk, 
#                                                              g = g, 
#                                                              c_mod_l = c_mod_l, 
#                                                              H_lk = H_lk, 
#                                                              full_output = 1))
# }
# 
# dt <- groupPolicies[ , .SD[which.min(obj_cost)], by = c("p1", "p2")]
# 
# 
# fwrite(dt, "Results//prev_2way_wtp1e7.csv")
# ###############################################







############################################
## ANALYSIS

library(ggplot2)
library(gridExtra)


dt <- fread("Results//prev_2way_wtp1e6.csv")

dt.all <- rbind(cbind(fread("Results//prev_2way_wtp1e6.csv"), WTP = 1e6),
                cbind(fread("Results//prev_2way_wtp1e5.csv"), WTP = 1e5),
                cbind(fread("Results//prev_2way_wtp1e7.csv"), WTP = 1e7))







policies <- as.data.table(policies)
colnames(policies) <- c("z", "a1", "a2", "a3", "a4", "m1")

policies[ , name := c("Defer", #1
                   "No intervention", #2
                   "Zika-ID", #3
                   "Zika-MP", #4
                   "Zika-ID-MP", #5
                   "WNV-ID", #6
                   "Zika-ID, WNV-ID", #7
                   "Zika-MP, WNV-ID", #8
                   "Zika-ID-MP, WNV-ID", #9
                   "WNV-MP", #10
                   "Zika-ID, WNV-MP", #11
                   "Zika-MP, WNV-MP", #12
                   "Zika-ID-MP, WNV-MP", #13
                   "WNV-ID-MP", #14
                   "Zika-ID, WNV-ID-MP", #15
                   "Zika-MP, WNV-ID-MP", #16
                   "Zika-ID-MP, WNV-ID-MP", #17
                   rep("FFP-PI etc", 16) #18 - 33
)
]


policies[ , pol_code := paste0(z, a1, a2, a3, a4, m1)]
dt[ , pol_code := paste0(z, a1, a2, a3, a4, m1)]
dt.all[ , pol_code := paste0(z, a1, a2, a3, a4, m1)]

setkey(policies, pol_code)
setkey(dt, pol_code)
setkey(dt.all, pol_code)
dt <- dt[policies[, c("pol_code", "name")], nomatch = 0]
dt.all <- dt.all[policies[, c("pol_code", "name")], nomatch = 0]

dt[ , pol_name := factor(name,
                                          levels = c(
                                            "Defer", #1
                                            "No intervention", #2
                                            "Zika-MP", #4
                                            "Zika-ID", 
                                            "Zika-ID-MP", #5
                                            "WNV-MP",
                                            "WNV-ID", #6
                                            "WNV-ID-MP", #14
                                            "Zika-MP, WNV-MP", #12
                                            "Zika-MP, WNV-ID", #8
                                            "Zika-MP, WNV-ID-MP", #16
                                            "Zika-ID, WNV-MP", #11
                                            "Zika-ID, WNV-ID", #7
                                            "Zika-ID, WNV-ID-MP", #15
                                            "Zika-ID-MP, WNV-MP", #13
                                            "Zika-ID-MP, WNV-ID", #9
                                            "Zika-ID-MP, WNV-ID-MP"
                                          ))]

dt.all[ , pol_name := factor(name,
                         levels = c(
                           "Defer", #1
                           "No intervention", #2
                           "Zika-MP", #4
                           "Zika-ID", 
                           "Zika-ID-MP", #5
                           "WNV-MP",
                           "WNV-ID", #6
                           "WNV-ID-MP", #14
                           "Zika-MP, WNV-MP", #12
                           "Zika-MP, WNV-ID", #8
                           "Zika-MP, WNV-ID-MP", #16
                           "Zika-ID, WNV-MP", #11
                           "Zika-ID, WNV-ID", #7
                           "Zika-ID, WNV-ID-MP", #15
                           "Zika-ID-MP, WNV-MP", #13
                           "Zika-ID-MP, WNV-ID", #9
                           "Zika-ID-MP, WNV-ID-MP"
                         ))]


colors <- c("Defer" = "black",
            "No intervention" = "cornsilk",
            "Zika-ID" = "firebrick3",
            "Zika-MP" = "firebrick1", #4
            "Zika-ID-MP" = "firebrick4", #5
            "WNV-ID" = "dodgerblue3", #6
            "Zika-ID, WNV-ID" = "violetred3", #7
            "Zika-MP, WNV-ID" = "purple3", #8
            "Zika-ID-MP, WNV-ID" = "indianred3", #9
            "WNV-MP" = "dodgerblue1", #10
            "Zika-ID, WNV-MP" = "violetred1", #11
            "Zika-MP, WNV-MP" = "purple1", #12
            "Zika-ID-MP, WNV-MP" = "indianred1", #13
            "WNV-ID-MP" = "dodgerblue4", #14
            "Zika-ID, WNV-ID-MP" = "violetred4", #15
            "Zika-MP, WNV-ID-MP" = "purple4", #16
            "Zika-ID-MP, WNV-ID-MP" = "indianred4"
            
)


WTP.labs <- c("WTP = $1,000,000 per QALY (basecase)", "WTP = $10,000,000 per QALY", "WTP = $100,000 per QALY")
names(WTP.labs) <- c(1e6, 1e7, 1e5)

#p_pol_by_WTP <- 
ggsave("Results/fig_opt_by_prev_WTP_scenarios.png",
       ggplot(data=dt.all)+
  facet_wrap(vars(WTP), ncol = 1, labeller = labeller(WTP = WTP.labs))+
  geom_tile(aes(x = p1, y = p2, fill = pol_name)) +
  scale_x_log10(name = "Zika prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_y_log10(name = "WNV prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  theme_bw()+
  theme(#legend.position = "none",
    panel.grid = element_blank())+
  scale_fill_manual(values = colors, name = "Policy"),
       width = 5,
       height = 6,
       unit = "in")











p_policy <- ggplot()+
  geom_tile(data =dt, aes(x = p1, y = p2, fill = pol_name)) +
  scale_x_log10(name = "Zika prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_y_log10(name = "WNV prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  theme_bw()+
  theme(#legend.position = "none",
        panel.grid = element_blank())+
  scale_fill_manual(values = colors, name = "Policy")
  
p_policy



p_obj_cost <- ggplot()+
  geom_tile(data=dt, aes(x = p1, y = p2, fill = obj_cost)) +
  scale_x_log10(name = "Zika prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_y_log10(name = "WNV prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_fill_gradient2(name = "Cost",
                       low = "white",
                       mid = "paleturquoise",
                       high = "navyblue",
                       labels = scales::dollar)+
  theme_bw()+
  theme(#legend.position = "none",
    panel.grid = element_blank())+
  ggtitle("D. Total net monetary cost (objective function)")


#p_obj_cost

p_rr1 <- ggplot()+
  geom_tile(data=dt, aes(x = p1, y = p2, fill = resid_risk1)) +
  scale_x_log10(name = "Zika prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_y_log10(name = "WNV prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_fill_gradient2(name = "Risk",
                       low = "white",
                       mid = "paleturquoise",
                       high = "navyblue")+
  theme_bw()+
  theme(#legend.position = "none",
    panel.grid = element_blank())+
  ggtitle("E. Zika residual risk per donation")

p_rr2 <- ggplot()+
  geom_tile(data=dt, aes(x = p1, y = p2, fill = resid_risk2)) +
  scale_x_log10(name = "Zika prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_y_log10(name = "WNV prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_fill_gradient2(name = "Risk",
                       low = "white",
                       mid = "paleturquoise",
                       high = "navyblue")+
  theme_bw()+
  theme(#legend.position = "none",
    panel.grid = element_blank())+
  ggtitle("F. WNV residual risk per donation")

p_test_cost <- ggplot()+
  geom_tile(data=dt, aes(x = p1, y = p2, fill = test_cost)) +
  scale_x_log10(name = "Zika prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_y_log10(name = "WNV prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_fill_gradient2(name = "Cost",
                       low = "white",
                       mid = "paleturquoise",
                       high = "navyblue",
                       labels = scales::dollar)+
  theme_bw()+
  theme(#legend.position = "none",
    panel.grid = element_blank())+
  ggtitle("A. Test cost per donation")

p_NMC <- ggplot()+
  geom_tile(data=dt, aes(x = p1, y = p2, fill = downsteam_NMC)) +
  scale_x_log10(name = "Zika prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_y_log10(name = "WNV prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_fill_gradient2(name = "Cost",
                       low = "white",
                       mid = "paleturquoise",
                       high = "navyblue",
                       labels = scales::dollar)+
  theme_bw()+
  theme(#legend.position = "none",
    panel.grid = element_blank())+
  ggtitle("C. Downstream net monetary cost")


p_pos_tests <- ggplot()+
  geom_tile(data=dt, aes(x = p1, y = p2, fill = cost_pos_tests)) +
  scale_x_log10(name = "Zika prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_y_log10(name = "WNV prevalence", 
                breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-1),
                expand = c(0, 0),
                minor_breaks = NULL) + 
  scale_fill_gradient2(name = "Cost",
                       low = "white",
                       mid = "paleturquoise",
                       high = "navyblue",
                       labels = scales::dollar)+
  theme_bw()+
  theme(#legend.position = "none",
    panel.grid = element_blank())+
  ggtitle("B. Cost due to donations testing positive")

p_pos_tests

ggsave("Results/outcomes_by_prevalence.png",
       grid.arrange(#p_policy + ggtitle("A. Optimal policy") + theme(legend.text=element_text(size=7)),
         p_test_cost, p_pos_tests,
         p_NMC, p_obj_cost,
         p_rr1, p_rr2,
             ncol = 2),
       width = 10,
       height = 11,
       units = "in"
)


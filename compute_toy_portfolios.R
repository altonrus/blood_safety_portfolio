### Optimal blood safety portfolio
library(data.table)
source("portfolio_functions.R")

### Variable definitions ---------------

K = 4 #number of TTIDs
I = 5 #number of donor segments
J = 3 #number of tests
N = 2 #number of risk reducing modifications

#Net present cost. 
c_k <- matrix(0L, K) 
c_k[1] <- 200 + 0.002*10E6
c_k[2] <- 2000 + 0.025*1E6
c_k[3] <- 350 + 0.11*1E6
c_k[4] <- 51000 + 0.07*1E6


#Probability of infectiousness by donor segment
P_ik <- matrix(0L, I,K) # probability of infectiousness for TTID k in donor segment i
P_ik[1,1] <- .001
P_ik[1,2] <- .0001
P_ik[1,3] <- .0002
P_ik[1,4] <- .00001
P_ik[2,1] <- .001
P_ik[2,2] <- .0004
P_ik[2,3] <- .0002
P_ik[2,4] <- .00002
P_ik[3,1] <- .001
P_ik[3,2] <- .0001
P_ik[3,3] <- .0002
P_ik[3,4] <- .0005
P_ik[4,1] <- .001
P_ik[4,2] <- .00015
P_ik[4,3] <- .004
P_ik[4,4] <- .00001
P_ik[5,1] <- .009
P_ik[5,2] <- .0001
P_ik[5,3] <- .0002
P_ik[5,4] <- .00001

#Replacement cost by donor segment
d_i <- matrix(0L, I) # cost of replacing donation from segment i if deferred
d_i[1] <- 90
d_i[2] <- 90
d_i[3] <- 90
d_i[4] <- 90
d_i[5] <- 90

#Processing cost by donor segment
w_i <- matrix(0L, I) # cost of processing donation from segment i if accepted
w_i[1] <- 20
w_i[2] <- 20
w_i[3] <- 20
w_i[4] <- 20
w_i[5] <- 20

#Number of donors in each segment
n_i <- matrix(0L, I) # number of donors in segment i
n_i[1] <- 75000
n_i[2] <- 12000
n_i[3] <- 8000
n_i[4] <- 4000
n_i[5] <- 1000

#Cost of disease marker tests
c_test_j <- matrix(0L, J) #per-donation cost of test j
c_test_j[1] <- 13
c_test_j[2] <- 17
c_test_j[3] <- 14

#Sensitivity of tests
R_jk <- matrix(0L, J,K) #sensitivity of test j for TTID k
R_jk[1,]<- c(0,.97,0,0) #no test = 0
R_jk[2,]<- c(0,0,.92,.93)
R_jk[3,]<- c(0,0,.98,0)

#Specificity of tests
Q_jk <- matrix(0L, J,K) #specificity of test j for TTID k
Q_jk[1,] <- c(1,.999,1,1) # no test = 1
Q_jk[2,] <- c(1,1,.99,.99)
Q_jk[3,] <- c(1,1,.98,1)

g = 60+90 #disposal cost, should be larger than deferral cost

#cost of modification 
c_mod_n <- matrix(0L, N) #per-donation cost of modification n
c_mod_n[1] <- 6
c_mod_n[2] <- 44

#Risk-reducing multiplier for modification
H_nk <- matrix(0L, N,K) #risk-reduction multiplier for TTID k from using modification n
H_nk[1,] <- c(.8, 1, .7, .6) #No effect = 1
H_nk[2,] <- c(1, .02, .04, .4)

#Decision variables
z_i <- matrix(1L, I) # decision variable. 1 if accepting donations from segment i
M_ni <- matrix(1L, N,I) # decision variable. 1 if using modification n in segment i
A_ji <- matrix(1L, J,I) # decision variable. 1 if using test j in segment i


# Extra stuff ------------------------
universal_only = 1

I = 5
N = 1
J = 4

#Count number of policies
n_pols = 0
if (universal_only == 0){
  for (i in 0:I){
    n_pols = n_pols + choose(n=I, k=i)*2^(i*(N+J))
  } 
  
  pol_length = I + J + N
} else {
  for (i in 0:I){
    n_pols = n_pols + choose(n=I, k=i)*2^(N+J)
  } 
  
  pol_length = I + J*I + N*I
}

# #Universal only ----------
# policies <- matrix(0 , n_pols, pol_length + 4 + K)
# dim(policies)
# 
# 
# n_segment_combos = 2^I
# segment_combos <- as.matrix(expand.grid(replicate(I, 0:1, simplify=FALSE))) 
# n_testmod_combos = 2^(N+J)
# testmod_combos <- as.matrix(expand.grid(replicate(N+J, 0:1, simplify=FALSE))) 
# 
# for (row_chunk in 1:n_segment_combos){
#   start = 1 + (row_chunk-1)*n_segment_combos
#   policies[start:(start+n_segment_combos-1), 1:I] <- segment_combos
#   policies[start:(start+n_segment_combos-1), (I+1):(I+I*N+I*J)] <- testmod_combos[row_chunk, ]
# }
# 
# for (row_chunk in 1:n_segment_combos){
#   start = 1 + (row_chunk-1)*n_segment_combos
#   policies[start:(start+n_segment_combos-1), 1:pol_length]  <- cbind(segment_combos, t(matrix(rep(rep(testmod_combos[row_chunk, ], each=(N+J)), n_segment_combos), , n_segment_combos)))
# }
# 
# for (row in 1:dim(policies)[1]) {
#   for (col in (I+1):(I+I*N+I*J)){
#     mod = ifelse((col - I) %% I == 0, I, (col - I) %% I)
#     
#     if(policies[row, mod] == 0){policies[row, col] = 0}
#   }
# }
# 
# for (row in 1:dim(policies)[1]) {
#   policies[row, (pol_length+1):(pol_length+4+K)] = cost_portfolio_unopt(z_i = matrix(policies[row, 1:I]), 
#                                                M_ni = t(matrix(policies[row, (I+1):(I+I*N)], I, N)), 
#                                                A_ji = t(matrix(policies[row, (I+I*N+1):(I+I*N+I*J)], I, J)), 
#                                                c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_n, H_nk, full_output = 1)
# }
# colnames(policies) <- c(paste0('z', 1:I),
#                         paste0('m', rep(1:N, each = I), rep(1:I, times = N)),
#                         paste0('a', rep(1:J, each = I), rep(1:I, times = J)),
#                         'obj cost', 'yield', 'test cost', 'mod cost', paste0('RR disease ', 1:K)
#                         )
# 
# 
# 
# write.csv(policies, "toyproblem1_universal2.csv")
# 
# 

# Tailored policies -----------

I = 5 #number of donor segments
J = 3 #number of tests
N = 2 #number of risk reducing modifications

# POLICIES WITH 2to4 donor segments
# First, enumerate all policies wieth between 2 and I - 1 donor segments
# n_pols = 0
# 
# for (i in 1:(I-1)){
#   n_pols = n_pols + choose(n=I, k=i)*2^(i*(N+J))
# } 
# 
# pol_length = I*(1+J+N)
# 
# policies <- matrix(0 , n_pols, pol_length + 1)
# 
# n_segment_combos = 2^I
# segment_combos <- as.data.table(expand.grid(replicate(I, 0:1, simplify=FALSE)))
# segment_combos[ , num_segs := Var1+Var2+Var3+Var4+Var5]
# segment_combos <- segment_combos[ order(num_segs)]
# segment_combos <- segment_combos[num_segs %in% 1:4]
# n_seg = table(segment_combos$num_segs)
# segment_combos <- as.matrix(segment_combos[ , 1:5])



# # Add all policies except those with 0 or all I segments
# pol_row = 1
# seg_row = 1
# for (seg_size in 1:(I-1)){ # Number of donor segments
#   tm_times_seg_size <- (J+N)*seg_size #Max possible number of test/mod - segment decisions for this number of donor segments
#   tm_combos <- as.matrix(expand.grid(replicate((J+N)*seg_size, 0:1, simplify=FALSE))) #All binary combo for (J+M)*(# seg)
#   num_segs_of_size = n_seg[seg_size][[1]] # Number of donor segment policies of size seg_size
#   for (tm_row in 1:(2^tm_times_seg_size)){
#     for (seg in seg_row : (seg_row + num_segs_of_size - 1)){
#       policies[pol_row, 1:I] = segment_combos[seg , ]
#       a <- rep(segment_combos[ seg ,  ], J+N)
#       policies[ pol_row, (I+1) : pol_length ] = replace(a, split(seq_along(a),a)$'1', tm_combos[ tm_row, ])
#       pol_row = pol_row + 1
#     }
#   }
#   #print(seg_row)
#   seg_row = seg_row + num_segs_of_size
# }
# rm(tm_combos)
# saveRDS(policies, "pols_1to4_nocost.rds")
# 
# policies <- readRDS("pols_1to4_nocost.rds")
# 
# t_0 <- Sys.time()
# for (row in 1:dim(policies)[1]) {
#   policies[row, pol_length+1] = cost_portfolio_unopt(z_i = matrix(policies[row, 1:I]),
#                                                      M_ni = t(matrix(policies[row, (I+1):(I+I*N)], I, N)), 
#                                                      A_ji = t(matrix(policies[row, (I+I*N+1):(I+I*N+I*J)], I, J)), 
#                                                      c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_n, H_nk, full_output = 0)
# }
# print(Sys.time() - t_0)
# 
# policies <- rbind(policies, 0)
# row = dim(policies)[1]
# policies[row, pol_length + 1 ] <- cost_portfolio_unopt(z_i = matrix(policies[row, 1:I]),
#                                              M_ni = t(matrix(policies[row, (I+1):(I+I*N)], I, N)), 
#                                              A_ji = t(matrix(policies[row, (I+I*N+1):(I+I*N+I*J)], I, J)), 
#                                              c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_n, H_nk, full_output = 0)
# saveRDS(policies, "pols_0to4.rds")
# rm(policies)


# Policies with 5 donor segments
n_pols = 2^(I*(N+J))
pol_length = I*(J+N)

policies <- as.matrix(expand.grid(replicate( pol_length,  0:1, simplify=FALSE)))
policies <- cbind(policies, 0)

t_0 <- Sys.time()
for (row in 1: n_pols){
  policies[row, pol_length+1] = portfolio_nodeferral(M_ni = t(matrix(policies[row, 1:(I*N)], I, N)), 
                                                     A_ji = t(matrix(policies[row, (I*N+1):(I*N+I*J)], I, J)), 
                                                     c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_n, H_nk)
}
print(Sys.time() - t_0)

saveRDS(policies, "pols_5.rds")


policies5 <- readRDS("pols_5.rds")
policies0to4 <- readRDS("pols_0to4.rds")

policies5 <- cbind(matrix(1, nrow = dim(policies5)[1], ncol = 5) , policies5)

policies <- rbind(policies0to4, policies5)
rm(policies5); rm(policies0to4)


policies <- as.data.table(policies)

colnames(policies) <- c(paste0('z', 1:I),
                        paste0('m', rep(1:N, each = I), rep(1:I, times = N)),
                        paste0('a', rep(1:J, each = I), rep(1:I, times = J)),
                        'obj_cost')
policies <- policies[order(obj_cost)]

summary(policies$obj_cost)

policies_top500 <- policies[1:500,]
#rm(policies)

policies_top500 <- as.matrix(policies_top500)
policies_top500 <- cbind( policies_top500, matrix(0, nrow = 500, ncol = 7))



for(row in 1:500){
  policies_top500[row, 31:38] <- cost_portfolio_unopt(z_i = matrix(policies_top500[row, 1:I]),
                                                      M_ni = t(matrix(policies_top500[row, (I+1):(I+I*N)], I, N)),
                                                      A_ji = t(matrix(policies_top500[row, (I+I*N+1):(I+I*N+I*J)], I, J)),
                                                      c_k, P_ik, d_i, w_i, n_i, c_test_j, R_jk, Q_jk, g, c_mod_n, H_nk, full_output = 1)

  
write.csv(policies_top500, "toyproblem1_tailored500.csv")

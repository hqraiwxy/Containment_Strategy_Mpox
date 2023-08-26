library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(readxl)

setwd('C:/Users/wxyrc/Desktop/mpox_sim')
source("./mpox_initialize.R")
source("./mpox_helperfunction.R")
## 实验参数设置
expriment_scenarios <- read_excel ("参数设置.xlsx",sheet = "Sheet3")

stochastic_NH <- function(parms, Ns, delta_t, t){
  beta = parms[["beta"]] 
  beta.HH = parms[["beta.HH"]] 
  sigma1 = parms[["sigma1"]] 
  sigma2 = parms[["sigma2"]] 
  gamma = parms[["gamma1"]]
  gamma = parms[["gamma2"]]
  I.C_time=parms[["I.C_time"]] 
  I.C = ifelse(t>I.C_time,parms[["I.C"]],0) # make I.C dependent on I.C_time
  k.HH = parms[["k.HH"]] 
  k.HL = parms[["k.HL"]]
  k.LH = parms[["k.LH"]]
  k.LL = parms[["k.LL"]] 
  int.r = parms[["int.r"]] 
  ID= parms[["id.I"]] # id.I= 2 （duration of pre-symptomatic transmission）
  int.r = parms[["int.r"]] # 1 if go to new_recovered.rNC, 0 if go to new_recovered.rC
  int.hcw = parms[["int.hcw"]] # 1 if recovered work with NC and 0 if work with C
  int.rooms = parms[["int.rooms"]] # if 1 put susceptible and recovered together in NC
  mu.C = parms[["mu.C"]] # COVID mortality
  mu.NC = parms[["mu.NC"]] # regular mortality
  prop_rhcwR = parms[["prop_rhcwR"]] 
  VL.PCR_threshold = parms[["VL.PCR_threshold"]] # threshold for detectable VL (put arbitrary # in for now)
  VL.Antigen_threshold = parms[["VL.Antigen_threshold"]]
  int_time = parms[["int_time"]] #int_time=21  这个指的是intervention开始的时间吗？
  V_doses = parms[["V_doses"]]
  boost_vax_coverage_r = parms[["boost_vax_coverage_r"]] #疫苗覆盖率？还是疫苗保护总效果？
  full_vax_coverage_r = parms[["full_vax_coverage_r"]] #疫苗覆盖率？还是疫苗保护总效果？
  refusal_r = parms[["refusal_r"]]
  boost_vax_coverage_hcw = parms[["boost_vax_coverage_hcw"]]
  full_vax_coverage_hcw = parms[["full_vax_coverage_hcw"]]
  refusal_hcw= parms[["refusal_hcw"]] #拒绝接种疫苗者
  S.hNM = Ns[["S.hNM"]]
  E.hNM = Ns[["E.hNM"]]
  A.hNM = Ns[["A.hNM"]]
  I.hNM = Ns[["I.hNM"]]
  R.hNM = Ns[["R.hNM"]]
  I.hM  = Ns[["I.hM"]]
  R.hM  = Ns[["R.hM"]]
  
  S.lNM = Ns[["S.lNM"]]
  E.lNM = Ns[["E.lNM"]]
  A.lNM = Ns[["A.lNM"]]
  I.lNM = Ns[["I.lNM"]]
  R.lNM = Ns[["R.lNM"]]
  I.lM  = Ns[["I.lM"]]
  R.lM  = Ns[["R.lM"]]
  
  S.hNM3  = Ns[["S.hNM3"]]
  S.lNM3  = Ns[["S.lNM3"]]
  
  death.r=Ns[["death.r"]]
  death.rC=Ns[["death.rC"]]
  NewE_h_total=Ns[["NewE_h_total"]] # h累计新增E
  NewE_l_total=Ns[["NewE_l_total"]] # l累计新增E
  NewA_h_total=Ns[["NewA_h_total"]] # h累计新增A
  NewI_h_total=Ns[["NewI_h_total"]] # h累计新增I
  NewA_l_total=Ns[["NewA_l_total"]] # l累计新增A
  NewI_l_total=Ns[["NewI_l_total"]] # l累计新增I
  
  inc_vax1_h=Ns[["inc_vax1_h"]] 
  inc_vax2_h=Ns[["inc_vax2_h"]]
  inc_vax1_l=Ns[["inc_vax1_l"]] 
  inc_vax2_l=Ns[["inc_vax2_l"]]
  total=Ns[["total"]]
  total_l=Ns[["total_l"]]

  N.hNM <- nrow(S.hNM) + nrow(E.hNM) + nrow(A.hNM) + nrow(I.hNM) + nrow(R.hNM)
  N.hM <- nrow(I.hM) + nrow(R.hM) 

  N.lNM <-  nrow(S.lNM) + nrow(E.lNM) + nrow(A.lNM) + nrow(I.lNM) + nrow(R.lNM)
  N.lM <-  nrow(I.lM) + nrow(R.lM)
  
  ###########################   A/I to R  ###########################################
  # move highrisk from A/I to R and NC to C
  # NM队列中的A
  recover.A.hNM <- recover_or_test_h(A.hNM, parms, "A")
  R.hNM <- rbind(R.hNM,recover.A.hNM[[1]])    # move A.hNM to R.hNM if recovered
  A.hNM <- recover.A.hNM[[3]]                 # keep rest of A.hNM in A.hNM
  A.hNM  <- rbind(A.hNM,recover.A.hNM[[2]])     # 未recovered的A.hNM仍然在原队列中
  
  # NM队列中的I
  recover.I.hNM <- recover_or_test_h(I.hNM, parms, "I")
  R.hNM <- rbind(R.hNM,recover.I.hNM[[1]])
  I.hM  <- rbind(I.hM,recover.I.hNM[[2]]) 
  I.hNM <- recover.I.hNM[[3]] 
  # M队列中的I
  recover.I.hM <- recover_I.M(I.hM, parms)
  if (int.r==1 & t>int_time){
    R.hNM <- rbind(R.hNM,recover.I.hM[[1]])   # 采取干预在M队列中recover的I移回到NM队列
  } else {
    R.hM  <- rbind(R.hM, recover.I.hM[[1]])   # 不采取干预的话不移回
  }
    I.hM <- recover.I.hM[[2]]                 # keep rest in I.hM

  # move LOWRISK from A/I to R and NC or C
  # NM队列中的A
  recover.A.lNM <- recover_or_test_l(A.lNM, parms, "A")
  R.lNM <- rbind(R.lNM,recover.A.lNM[[1]]) # move A.lNM to R.lNM if recover
  A.lNM <- recover.A.lNM[[3]] # keep rest of A.lNM in A.lNM
  A.lNM  <- rbind(A.lNM,recover.A.lNM[[2]]) # 未recovered的A.lNM仍然在原队列中
  
  # NM队列中的I
  recover.I.lNM <- recover_or_test_l(I.lNM, parms, "I")
  R.lNM <- rbind(R.lNM,recover.I.lNM[[1]]) # move I.lNM to R.lNM if recover
  I.lM <- rbind(I.lM,recover.I.lNM[[2]]) # move I.lNM to I.lM if test positive
  I.lNM <- recover.I.lNM[[3]] # keep rest of I.lNM in I.lNM
 
  
  # M队列中的I
  recover.I.lM <- recover_I.M(I.lM, parms)
  if (int.r==1 & t>int_time){
    R.lNM <- rbind(R.lNM,recover.I.lM[[1]])   # 采取干预在M队列中recover的I移回到NM队列
  } else {
    R.lM  <- rbind(R.lM, recover.I.lM[[1]])   # 不采取干预的话不移回
  }
  I.lM <- recover.I.lM[[2]]                 # 
  
  ###########################   E -> A/I  ###########################################
  # 高风险NM队列中的E
  
  infect_E.hNM <- E_to_I(E.hNM,parms)
  A.hNM <- rbind(A.hNM,infect_E.hNM[[1]])
  I.hNM <- rbind(I.hNM,infect_E.hNM[[2]])
  E.hNM <- infect_E.hNM[[3]]
  
  
  NewA_h <- nrow(infect_E.hNM[[1]]) # 每日新增A
  NewI_h <- nrow(infect_E.hNM[[2]]) # 每日新增I
  
  NewA_h_total = (NewA_h_total + NewA_h)
  NewI_h_total = (NewI_h_total + NewI_h)
 
  #cat("infect",(nrow(S.hNM) + nrow(E.hNM) + nrow(A.hNM) + nrow(I.hNM) + nrow(R.hNM) + nrow(I.hM) + nrow(R.hM)),t,"\n")
  # 低风险NM队列中的E
  infect_E.lNM <- E_to_I(E.lNM,parms)
  A.lNM <- rbind(A.lNM,infect_E.lNM[[1]])
  I.lNM <- rbind(I.lNM,infect_E.lNM[[2]])
  E.lNM <- infect_E.lNM[[3]]
  
  NewA_l <- nrow(infect_E.lNM[[1]]) # 每日新增A
  NewI_l <- nrow(infect_E.lNM[[2]]) # 每日新增I
  
  NewA_l_total = (NewA_l_total + NewA_l)
  NewI_l_total = (NewI_l_total + NewI_l)
  

  ###########################   S -> E  ###########################################
  ## NM队列中的high risk S
  S.hNM2 <- S.hNM

  if(nrow(S.hNM)>=1){
    for (i in 1:nrow(S.hNM)){
      S.hNM.exp <- rbinom(1,1, I.C) # S来自社会面的暴露
     
      # adjust betas based on vaccination status
      beta.HH.v <- beta.HH*(1-S.hNM$VE_s[i])
      beta.v <- beta*(1-S.hNM$VE_s[i]) # reduce susceptibility by VEs (assuming leaky here)

      exp.prob <- beta.v*S.hNM.exp +
        beta.HH.v*(ifelse(N.hNM>0, k.HH*nrow(I.hNM)/N.hNM,0)) +  #【*(nrow(I.hNM))/N.hNM】(nrow(I.hNM)+nrow(I.hM))这里任意队列的I都在传播
        beta.v*(ifelse(N.lNM>0,k.HL*nrow(I.lNM)/N.lNM,0)) #*nrow(I.lNM)/N.lNM
      if (exp.prob >= 1) {
        exp.prob <- 1
      } else if (exp.prob < 0) {
        exp.prob <- 0
      }

      S.hNM.exposed <- rbinom(1,1,exp.prob)  #确定这个S是否进展到E

      if (S.hNM.exposed==1){  #条件：这个S进展到E；给定Inc.pd和Days=0初始状态
        VE_post_h = ifelse(S.hNM$V_doses[i]==0, parms[["VE_p"]]*rbinom(1,1,parms[["post_vax_coverage_h"]]), 0)
        
        S.hNM3 <- rbind(S.hNM3,cbind(ID=S.hNM$ID[i], 
                                     VE_s = S.hNM$VE_s[i],
                                     VE_p = VE_post_h,
                                     V_doses = ifelse( VE_post_h !=0, 2 , S.hNM$V_doses[i]),
                                     Inc.pd= t+round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), 
                                     Days=t))
        E.hNM <- rbind(E.hNM,cbind(ID=S.hNM$ID[i], 
                                    VE_s = S.hNM$VE_s[i],
                                    VE_p = VE_post_h,
                                    V_doses = ifelse( VE_post_h !=0, 2 , S.hNM$V_doses[i]),
                                    Inc.pd= t+round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), 
                                    Days=t))
        S.hNM2 %>%
          subset(ID != S.hNM$ID[i]) -> S.hNM2 # 从S中去除变成E的
      }
    }
  }

  expose.S.hNM <- nrow(S.hNM) - nrow(S.hNM2) #新E
  S.hNM <- S.hNM2 #剩下的S
  
  ## NM队列中的low risk S
  S.lNM2 <- S.lNM
  if(nrow(S.lNM)>=1){
    for (i in 1:nrow(S.lNM)){

      S.lNM.exp <- rbinom(1,1, I.C) # S来自社会面的暴露
      
      # adjust betas based on vaccination status
      beta.v <- beta*(1-S.lNM$VE_s[i]) # reduce susceptibility by VEs (assuming leaky here)
      
      exp.prob <- beta.v*S.lNM.exp +
        beta.v*(ifelse(N.hNM>0, k.LH*(nrow(I.hNM))/N.hNM,0)) + #*(nrow(I.hNM))/N.hNM 
        beta.v*(ifelse(N.lNM>0,k.LL*nrow(I.lNM)/N.lNM,0)) #*nrow(I.lNM)/N.lNM
      if (exp.prob >= 1){
        exp.prob <- 1
      } else if (exp.prob <0){
        exp.prob <- 0
      }
      
      S.lNM.exposed <- rbinom(1,1,exp.prob)  #确定这个S是否进展到E
      #cat("expose",S.hNM.exposed,nrow(S.hNM),nrow(E.hNM),(nrow(S.hNM) + nrow(E.hNM) + nrow(A.hNM) + nrow(I.hNM) + nrow(R.hNM) + nrow(I.hM) + nrow(R.hM)),t,"\n")
      #cat(i,S.hNM.exposed,"\n")
      if (S.lNM.exposed==1){  #条件：这个S进展到E；给定Inc.pd和Days=0初始状态
        VE_post_l = ifelse(S.lNM$V_doses[i]==0, parms[["VE_p"]]*rbinom(1,1,parms[["post_vax_coverage_l"]]), 0)
        
        S.lNM3 <- rbind(S.lNM3,cbind(ID=S.lNM$ID[i], 
                                     VE_s = S.lNM$VE_s[i],
                                     VE_p = VE_post_l,
                                     V_doses = ifelse( VE_post_l!=0, 2 , S.lNM$V_doses[i]),
                                     Inc.pd= t+round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), 
                                     Days=t))
        E.lNM <- rbind(E.lNM,cbind(ID=S.lNM$ID[i], 
                                    VE_s = S.lNM$VE_s[i],
                                    VE_p = VE_post_l,
                                    V_doses = ifelse( VE_post_l!=0, 2 , S.lNM$V_doses[i]),
                                    Inc.pd= t+round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), 
                                    Days=t))
        S.lNM2 %>%
          subset(ID != S.lNM$ID[i]) -> S.lNM2 # 从S中去除变成E的
      }
    }
  }
  
  expose.S.lNM <- nrow(S.lNM) - nrow(S.lNM2) #新E
  S.lNM <- S.lNM2 #剩下的S
  
  NewE_h <- expose.S.hNM # 每日新增E
  NewE_l <- expose.S.lNM # 每日新增E
  
  NewE_h_total = (NewE_h_total + NewE_h)
  NewE_l_total = (NewE_l_total + NewE_l)

  final <- list("S.hNM"=S.hNM,
                "E.hNM"=E.hNM,
                "A.hNM"=A.hNM,
                "I.hNM"=I.hNM,
                "R.hNM"=R.hNM,
                "I.hM"=I.hM,
                "R.hM"=R.hM,
                "S.lNM"=S.lNM,
                "E.lNM"=E.lNM,
                "A.lNM"=A.lNM,
                "I.lNM"=I.lNM,
                "R.lNM"=R.lNM,
                "I.lM"=I.lM,
                "R.lM"=R.lM,
                "S.hNM3"=S.hNM3,
                "S.lNM3"=S.lNM3,
                "NewE_h"=NewE_h,
                "NewE_l"=NewE_l,
                "NewE_h_total"=NewE_h_total,
                "NewE_l_total"=NewE_l_total,
                "total"=total,
                "total_l"=total_l,
                "NewA_h_total"=NewA_h_total,
                "NewI_h_total"=NewI_h_total,
                "NewA_l_total"=NewA_l_total,
                "NewI_l_total"=NewI_l_total)
  return(final)
}
# 初始化人群数量
inits=c(S.hNM.init = 495,
        E.hNM.init = 5,
        A.hNM.init = 0,
        I.hNM.init = 0,
        R.hNM.init = 0,
        I.hM.init  = 0,
        R.hM.init  = 0,
        S.lNM.init = 4999,
        E.lNM.init = 1,
        A.lNM.init = 0,
        I.lNM.init = 0,
        R.lNM.init = 0,
        I.lM.init  = 0,
        R.lM.init  = 0)

nsim <-50
t_step <- 1
dt <- seq(0,180,t_step)

### run simulations ----
parms <- list(
           sigma1=5,      ## incubation period shortest
           sigma2=10,     ## incubation period longest
           gamma1=14,     ## infectious period shortest
           gamma2=28,     ## infectious period longest
           I.C=0.00002,    ## probability of infection from community
           beta = 0.0028, ## 0.0032，0.0076beta for LL LH HL
           beta.HH = 0.0132, ## 0.0132，0.0229，beta for HH
           k.HH=9,        ## n contacts between highrisk
           k.HL=3,       ## n lowrisk contacted by each highrisk
           k.LH=1,        ## n highrisk contacted by each lowrisk
           k.LL=2,       ## n daily contacts between lowrisk
           pre_vax_coverage_h = 0.3,
           post_vax_coverage_h = 1,
           pre_vax_coverage_l = 0,
           post_vax_coverage_l = 0,
           VE_p = 0.4,
           VE_s = 0.8,## pre疫苗的保护效果
           move_M_ratio_h = 0, # highrisk症状出现后移动到M队列的比例
           move_M_ratio_l = 0, # lowrisk症状出现后移动到M队列的比例
           id.I_lower = 1, # duration of pre-symptomatic transmission 1-4
           id.I_upper = 4, 
           int_time=0,      #队列转移开始日期 intervention time
           int.r =1,
           I.C_time=0)

#####loop over expriment scenarios
for (es in seq_along(1:nrow(expriment_scenarios))) {
  random_seed = as.numeric(expriment_scenarios[es,1])
  set.seed(random_seed)
  
  parms[["pre_vax_coverage_h"]] <- as.numeric(expriment_scenarios[es,2])
  parms[["post_vax_coverage_h"]] <- as.numeric(expriment_scenarios[es,3])
  parms[["pre_vax_coverage_l"]] <- as.numeric(expriment_scenarios[es,4])
  parms[["post_vax_coverage_l"]] <- as.numeric(expriment_scenarios[es,5])
  parms[["move_M_ratio_h"]] <- as.numeric(expriment_scenarios[es,6])
  parms[["move_M_ratio_h"]] <- as.numeric(expriment_scenarios[es,6])

res_master <- NULL

for (sim in 1:nsim) {
  cat(sim,  "\n")
  Ns <- initialize(inits, parms, dt)
  res <-
    as.data.frame(matrix(nrow = length(dt), ncol = length(Ns) + 1))
  VLs <- NULL
  
  for (i in 1:length(dt)) {
    cat(i, "\n")
    debug <- Ns
    final <- stochastic_NH(parms, Ns, t_step, (i - 1) * t_step)
    Ns <- final
    res[i,] <- c(
      i,
      nrow(Ns[["S.hNM"]]),
      nrow(Ns[["E.hNM"]]),
      nrow(Ns[["A.hNM"]]),
      nrow(Ns[["I.hNM"]]),
      nrow(Ns[["R.hNM"]]),
      nrow(Ns[["I.hM"]]),
      nrow(Ns[["R.hM"]]),
      nrow(Ns[["S.lNM"]]),
      nrow(Ns[["E.lNM"]]),
      nrow(Ns[["A.lNM"]]),
      nrow(Ns[["I.lNM"]]),
      nrow(Ns[["R.lNM"]]),
      nrow(Ns[["I.lM"]]),
      nrow(Ns[["R.lM"]]),
      nrow(Ns[["S.hNM3"]]),
      nrow(Ns[["S.lNM3"]]),
      Ns[["NewE_h"]],
      Ns[["NewE_h_total"]],
      Ns[["NewA_h_total"]],
      Ns[["NewI_h_total"]],
      Ns[["NewE_l"]],
      Ns[["NewE_l_total"]],
      Ns[["NewA_l_total"]],
      Ns[["NewI_l_total"]],
      Ns[["total"]],
      Ns[["total_l"]]
    )
  }
  res %>%
    add_column("Sim" = sim) %>%
    add_column("seed" = random_seed) %>%
    add_column("pre_h" = parms[["pre_vax_coverage_h"]]) %>%
    add_column("post_h" = parms[["post_vax_coverage_h"]]) %>%
    add_column("pre_l" = parms[["pre_vax_coverage_l"]]) %>%
    add_column("post_l" = parms[["post_vax_coverage_l"]]) %>%
    add_column("move" = parms[["move_M_ratio_h"]]) %>%
    bind_rows(res_master) -> res_master
}

colnames(res_master) <- c("time", 
                          "S.hNM", "E.hNM", "A.hNM", "I.hNM", "R.hNM",
                          "I.hM", "R.hM",
                          "S.lNM", "E.lNM", "A.lNM", "I.lNM", "R.lNM",
                          "I.lM", "R.lM", 
                          "S.hNM3","S.lNM3", 
                          "NewE_h", "NewE_h_total", "NewA_h_total", "NewI_h_total",
                          "NewE_l", "NewE_l_total", "NewA_l_total", "NewI_l_total",
                          "total","total_l",
                          "Sim","seed","pre_h","post_h","pre_l","post_l","move")
write.csv(res_master, paste0("0504 seed-",random_seed," mpox.csv"))
}


initialize <- function(inits,params,dt){
  S.hNM.init = inits["S.hNM.init"]
  E.hNM.init = inits["E.hNM.init"]
  A.hNM.init = inits["A.hNM.init"]
  I.hNM.init = inits["I.hNM.init"]
  R.hNM.init = inits["R.hNM.init"]
  I.hM.init  = inits["I.hM.init"]
  R.hM.init = inits["R.hM.init"]
  
  S.lNM.init = inits["S.lNM.init"]
  E.lNM.init = inits["E.lNM.init"]
  A.lNM.init = inits["A.lNM.init"]
  I.lNM.init = inits["I.lNM.init"]
  R.lNM.init = inits["R.lNM.init"]
  I.lM.init  = inits["I.lM.init"]
  R.lM.init = inits["R.lM.init"]
  
  total <- 0
  ## 初始化S.hNM
  if (S.hNM.init > 0){
    S.hNM=as.data.frame(cbind("ID"=1:S.hNM.init))
    S.hNM %>%
      mutate(V_doses = sample(c(1,0),
                              size=S.hNM.init,
                              prob=c(parms[["pre_vax_coverage_h"]],1-parms[["pre_vax_coverage_h"]]),
                              replace=TRUE)
      ) -> S.hNM
    
    S.hNM %>%
      mutate(VE_s = V_doses * parms[["VE_s"]]
      ) -> S.hNM
    
    total <- total + S.hNM.init
  } else{
    S.hNM=as.data.frame(cbind("ID"=as.numeric(), 
                              "V_doses" = as.numeric(),
                              "VE_s" = as.numeric()) ) 
  }
  ## 初始化E.hNM
  if (E.hNM.init > 0){
    E.hNM=as.data.frame(cbind("ID"=(total+1):(total + E.hNM.init), 
                              "Inc.pd"=round(runif(E.hNM.init, parms[["sigma1"]], parms[["sigma2"]])),
                              "Days"=rep(0,E.hNM.init)))
    E.hNM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=E.hNM.init,
                              prob=c(parms[["post_vax_coverage_h"]],parms[["pre_vax_coverage_h"]],
                                     max(0,1-parms[["post_vax_coverage_h"]]-parms[["pre_vax_coverage_h"]])),
                              replace=TRUE)
      ) -> E.hNM
    E.hNM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> E.hNM

    total <- total + E.hNM.init
  } else{
    E.hNM=as.data.frame(cbind("ID"=as.numeric(),
                              "V_doses" = as.numeric(),"VE_s" = as.numeric(),"VE_p" = as.numeric(),
                              "Inc.pd"=as.numeric(),"Days"=as.numeric())) 
  }
  ## 初始化A.hNM 
  if (A.hNM.init > 0){
    A.hNM=as.data.frame(cbind("ID"=(total+1):(total+ A.hNM.init),
                              "removal.pd"=rep(NA, A.hNM.init),
                              "Inc.pd"=rep(0,A.hNM.init), #新增
                              "Inf.pd"=round(runif(A.hNM.init, parms[["gamma1"]], parms[["gamma2"]])),
                              "Inf.days"=rep(0,A.hNM.init), # can change if we want to assume they are partway into inf period
                              "ID.days"=rep(0,A.hNM.init)))
    A.hNM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=A.hNM.init,
                              prob=c(parms[["post_vax_coverage_h"]],parms[["pre_vax_coverage_h"]],
                                     max(0,1-parms[["post_vax_coverage_h"]]-parms[["pre_vax_coverage_h"]])),
                              replace=TRUE) 
    ) -> A.hNM
    
    A.hNM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> A.hNM
    
    total <- total + A.hNM.init
  } else{
    A.hNM=as.data.frame(cbind("ID"=as.numeric(),
                              "removal.pd"=as.numeric(),
                              "V_doses" = as.numeric(),
                              "VE_s" = as.numeric(),
                              "VE_p" = as.numeric(),
                              "Inc.pd"=as.numeric(),#新增
                              "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                              "ID.days"=as.numeric())) 
  }
  
  if (I.hNM.init > 0){
    I.hNM=as.data.frame(cbind("ID"=(total+1):(total+ I.hNM.init),
                              "removal.pd"=rep(NA, I.hNM.init),
                              "VE_s" = rep(0,I.hNM.init),
                              "Inc.pd"=rep(0,I.hNM.init),#新增
                              "Inf.pd"=round(runif(I.hNM.init, parms[["gamma1"]], parms[["gamma2"]])),
                              "Inf.days"=rep(0,I.hNM.init), # can change if we want to assume they are partway into inf period
                              "ID.days"=rep(0,I.hNM.init)))
    I.hNM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=I.hNM.init,
                              prob=c(parms[["post_vax_coverage_h"]],parms[["pre_vax_coverage_h"]],
                                     max(0,1-parms[["post_vax_coverage_h"]]-parms[["pre_vax_coverage_h"]])),
                              replace=TRUE)
      ) -> I.hNM
    
    I.hNM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> I.hNM
    total <- total + I.hNM.init
  } else{
    I.hNM=as.data.frame(cbind("ID"=as.numeric(),
                              "removal.pd"=as.numeric(),
                              "VE_s" = as.numeric(),
                              "V_doses" = as.numeric(),
                              "VE_s" = as.numeric(),
                              "VE_p" = as.numeric(),
                              "Inc.pd"=as.numeric(),#新增
                              "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                              "ID.days"=as.numeric()))
  }
  
  if (R.hNM.init > 0){
    R.hNM=as.data.frame(cbind("ID"=(total+1):(total + R.hNM.init),
                              "Inc.pd"=rep(0,R.hNM.init),#新增 
                              "Inf.pd"=rep(0,R.hNM.init),#新增
                              "Rec.days"=rep(0,R.hNM.init)))
    
    R.hNM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=R.hNM.init,
                              prob=c(parms[["post_vax_coverage_h"]],parms[["pre_vax_coverage_h"]],
                                     max(0,1-parms[["post_vax_coverage_h"]]-parms[["pre_vax_coverage_h"]])),
                              replace=TRUE)
      ) -> R.hNM
    
    R.hNM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> R.hNM
    
    total <- total + R.hNM.init
  } else{
    R.hNM=as.data.frame(cbind("ID"=as.numeric(),
                              "V_doses" = as.numeric(), 
                              "VE_s" = as.numeric(),
                              "VE_p" = as.numeric(),
                              "Inc.pd"=as.numeric(), "Inf.pd"=as.numeric(),#新增
                              "Rec.days"=as.numeric()))
  }
  
  if (I.hM.init > 0){
    I.hM=as.data.frame(cbind("ID"=(total+1):(total+ I.hM.init),
                             "removal.pd"=rep(NA, I.hM.init),
                             "VE_s" = rep(0,I.hM.init),
                             "Inc.pd"=rep(0,I.hM.init),#新增
                             "Inf.pd"=round(runif(I.hM.init, parms[["gamma1"]], parms[["gamma2"]])),
                             "Inf.days"=rep(0,I.hM.init),
                             "ID.days"=rep(0,I.hM.init))) # can change if we want to assume they are partway into inf period
    I.hM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=I.hM.init,
                              prob=c(parms[["post_vax_coverage_h"]],parms[["pre_vax_coverage_h"]],
                                     max(0,1-parms[["post_vax_coverage_h"]]-parms[["pre_vax_coverage_h"]])),
                              replace=TRUE)
      ) -> I.hM
    I.hM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> I.hM
    total <- total + I.hM.init
  } else{
    I.hM=as.data.frame(cbind("ID"=as.numeric(),
                             "removal.pd"=as.numeric(),
                             "VE_s" = as.numeric(),
                             "VE_p" = as.numeric(),
                             "V_doses" = as.numeric(),
                             "Inc.pd"=as.numeric(),#新增
                             "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                             "ID.days"=as.numeric())) 
    
  }
  
  if (R.hM.init > 0){
    R.hM=as.data.frame(cbind("ID"=(total+1):(total + R.hM.init),
                             "Inc.pd"=rep(0,R.hM.init), 
                             "Inf.pd"=rep(0,R.hM.init),#新增
                             "Rec.days"=rep(0,R.hM.init)))
    R.hM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=R.hM.init,
                              prob=c(parms[["post_vax_coverage_h"]],parms[["pre_vax_coverage_h"]],
                                     max(0,1-parms[["post_vax_coverage_h"]]-parms[["pre_vax_coverage_h"]])),
                              replace=TRUE)
      ) -> R.hM
    
    R.hM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> R.hM
    
    total <- total + R.hM.init
  } else{
    R.hM=as.data.frame(cbind("ID"=as.numeric(),
                             "V_doses" = as.numeric(),
                             "VE_s" = as.numeric(),
                             "VE_p" = as.numeric(),
                             "Inc.pd"=as.numeric(), "Inf.pd"=as.numeric(),#新增
                             "Rec.days"=as.numeric())) 
  }
  
  # hNM <- c(S.hNM$ID, E.hNM$ID, A.hNM$ID, I.hNM$ID,R.hNM$ID)
  # hM <- c(I.hM$ID,R.hM$ID)###################################################新增，这个要不要？
  # 
  total_l <- 9999 # l的编号从10000开始，为什么要在这里赋值？
  if (S.lNM.init > 0){
    S.lNM=as.data.frame(cbind("ID"=(total_l + 1):(total_l + S.lNM.init)))
    S.lNM %>%
      mutate(V_doses = sample(c(1,0),
                              size=S.lNM.init,
                              prob=c(parms[["pre_vax_coverage_l"]],1-parms[["pre_vax_coverage_l"]]),
                              replace=TRUE)
      ) -> S.lNM
    
    S.lNM %>%
      mutate(VE_s = V_doses * parms[["VE_s"]]
      ) -> S.lNM
    
    total_l <- total_l + S.lNM.init
  } else{
    S.lNM=as.data.frame(cbind("ID"=as.numeric(),
                              "VE_s"=as.numeric(),
                              "V_doses" = as.numeric()))
  }
  
  if (E.lNM.init > 0){
    E.lNM=as.data.frame(cbind("ID"=(total_l+1):(total_l + E.lNM.init),
                              "Inc.pd"=round(runif(E.lNM.init, parms[["sigma1"]], parms[["sigma2"]])),
                              "Days"=rep(0,E.lNM.init)))
    E.lNM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=E.lNM.init,
                              prob=c(parms[["post_vax_coverage_l"]],parms[["pre_vax_coverage_l"]],
                                     max(0,1-parms[["post_vax_coverage_l"]]-parms[["pre_vax_coverage_l"]])),
                              replace=TRUE)
      ) -> E.lNM
    E.lNM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> E.lNM
    total_l <- total_l + E.lNM.init
  } else{
    E.lNM=as.data.frame(cbind("ID"=as.numeric(),
                              "V_doses" = as.numeric(),
                              "VE_s" = as.numeric(),"VE_p" = as.numeric(),
                              "Inc.pd"=as.numeric(),"Days"=as.numeric())) 
  }
  
  if (A.lNM.init > 0){
    A.lNM=as.data.frame(cbind("ID"=(total_l+1):(total_l+ A.lNM.init),
                              "removal.pd"=rep(NA, A.lNM.init),
                              "Inc.pd"=rep(0,A.lNM.init), #新增
                              "Inf.pd"=round(runif(A.lNM.init, parms[["gamma1"]], parms[["gamma2"]])),
                              "Inf.days"=rep(0,A.lNM.init), # can change if we want to assume they are partway into inf period
                              "ID.days"=rep(0,A.lNM.init)))
    A.lNM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=A.lNM.init,
                              prob=c(parms[["post_vax_coverage_l"]],parms[["pre_vax_coverage_l"]],
                                     max(0,1-parms[["post_vax_coverage_l"]]-parms[["pre_vax_coverage_l"]])),
                              replace=TRUE)
      ) -> A.lNM
    A.lNM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> A.lNM
    total_l <- total_l + A.lNM.init
  } else{
    A.lNM=as.data.frame(cbind("ID"=as.numeric(),
                              "removal.pd"=as.numeric(),
                              "V_doses" = as.numeric(),
                              "VE_s" = as.numeric(),
                              "VE_p" = as.numeric(),
                              "Inc.pd"=as.numeric(),#新增
                              "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                              "ID.days"=as.numeric())) 
  }
  
  if (I.lNM.init > 0){
    I.lNM=as.data.frame(cbind("ID"=(total_l+1):(total_l+ I.lNM.init),
                              "removal.pd"=rep(NA, I.lNM.init),
                              "VE_s" = rep(0,I.lNM.init),
                              "Inc.pd"=rep(0,I.lNM.init),#新增
                              "Inf.pd"=round(runif(I.lNM.init, parms[["gamma1"]], parms[["gamma2"]])),
                              "Inf.days"=rep(0,I.lNM.init), # can change if we want to assume they are partway into inf period
                              "ID.days"=rep(0,I.lNM.init)))
    I.lNM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=I.lNM.init,
                              prob=c(parms[["post_vax_coverage_l"]],parms[["pre_vax_coverage_l"]],
                                     max(0,1-parms[["post_vax_coverage_l"]]-parms[["pre_vax_coverage_l"]])),
                              replace=TRUE)
      ) -> I.lNM
    I.lNM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> I.lNM
    total_l <- total_l + I.lNM.init
  } else{
    I.lNM=as.data.frame(cbind("ID"=as.numeric(),
                              "removal.pd"=as.numeric(),
                              "VE_s" = as.numeric(),
                              "VE_p" = as.numeric(),
                              "V_doses" = as.numeric(),
                              "Inc.pd"=as.numeric(),#新增
                              "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                              "ID.days"=as.numeric()))
  }
  
  if (R.lNM.init > 0){
    R.lNM=as.data.frame(cbind("ID"=(total_l+1):(total_l + R.lNM.init),
                              "Inc.pd"=rep(0,R.lNM.init),#新增 
                              "Inf.pd"=rep(0,R.lNM.init),#新增
                              "Rec.days"=rep(0,R.lNM.init)))
    
    R.lNM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=R.lNM.init,
                              prob=c(parms[["post_vax_coverage_l"]],parms[["pre_vax_coverage_l"]],
                                     max(0,1-parms[["post_vax_coverage_l"]]-parms[["pre_vax_coverage_l"]])),
                              replace=TRUE)
      ) -> R.lNM
    
    R.lNM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> R.lNM
    
    total_l <- total_l + R.lNM.init
  } else{
    R.lNM=as.data.frame(cbind("ID"=as.numeric(),
                              "V_doses" = as.numeric(),
                              "VE_s" = as.numeric(),
                              "VE_p" = as.numeric(),
                              "Inc.pd"=as.numeric(), "Inf.pd"=as.numeric(),#新增
                              "Rec.days"=as.numeric()))
  }
  
  if (I.lM.init > 0){
    I.lM=as.data.frame(cbind("ID"=(total_l+1):(total_l+ I.lM.init),
                             "removal.pd"=rep(NA, I.lM.init),
                             "VE_s" = rep(0,I.lM.init),
                             "Inc.pd"=rep(0,I.lM.init),#新增
                             "Inf.pd"=round(runif(I.lM.init, parms[["gamma1"]], parms[["gamma2"]])),
                             "Inf.days"=rep(0,I.lM.init),
                             "ID.days"=rep(0,I.lM.init))) # can change if we want to assume they are partway into inf period
    I.lM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=I.lM.init,
                              prob=c(parms[["post_vax_coverage_l"]],parms[["pre_vax_coverage_l"]],
                                     max(0,1-parms[["post_vax_coverage_l"]]-parms[["pre_vax_coverage_l"]])),
                              replace=TRUE)
      ) -> I.lM
    I.lM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> I.lM
    
    total_l <- total_l + I.lM.init
  } else{
    I.lM=as.data.frame(cbind("ID"=as.numeric(),
                             "removal.pd"=as.numeric(),
                             "VE_s" = as.numeric(),
                             "VE_p" = as.numeric(),
                             "V_doses" = as.numeric(),
                             "Inc.pd"=as.numeric(),#新增
                             "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                             "ID.days"=as.numeric())) 
    
  }
  
  if (R.lM.init > 0){
    R.lM=as.data.frame(cbind("ID"=(total_l+1):(total_l + R.lM.init),
                             "Inc.pd"=rep(0,R.lM.init), 
                             "Inf.pd"=rep(0,R.lM.init),#新增
                             "Rec.days"=rep(0,R.lM.init)))
    R.lM %>%
      mutate(V_doses = sample(c(2,1,0),
                              size=R.lM.init,
                              prob=c(parms[["post_vax_coverage_l"]],parms[["pre_vax_coverage_l"]],
                                     max(0,1-parms[["post_vax_coverage_l"]]-parms[["pre_vax_coverage_l"]])),
                              replace=TRUE)
      ) -> R.lM
    
    R.lM %>%
      mutate(VE_s = case_when(V_doses==2 ~ 0,
                              V_doses==1 ~ parms[["VE_s"]],
                              V_doses==0 ~ 0),
             VE_p = case_when(V_doses==2 ~ parms[["VE_p"]],
                              V_doses==1 ~ 0,
                              V_doses==0 ~ 0)) -> R.lM
    
    total_l <- total_l + R.lM.init
  } else{
    R.lM=as.data.frame(cbind("ID"=as.numeric(),
                             "V_doses" = as.numeric(),
                             "VE_s" = as.numeric(),
                             "VE_p" = as.numeric(),
                             "Inc.pd"=as.numeric(), "Inf.pd"=as.numeric(),#新增
                             "Rec.days"=as.numeric())) 
  }
  
  ###########NM3是不是看有多少接种疫苗的？定义需明确  当时是为了统计疫苗接种数量引入的，在
  S.hNM3 = as.data.frame(cbind("ID"=as.numeric(), 
                               "VE_s"==as.numeric(), 
                               "VE_p"==as.numeric(), 
                               "V_doses" = as.numeric(),                        
                               "Inc.pd"=as.numeric(), 
                               "Days"=as.numeric())) 
  S.lNM3 = as.data.frame(cbind("ID"=as.numeric(), 
                               "VE_s"==as.numeric(), 
                               "VE_p"==as.numeric(), 
                               "V_doses" = as.numeric(),                        
                               "Inc.pd"=as.numeric(), 
                               "Days"=as.numeric())) 
  ##############免疫情况
  vax_count1_h <- nrow(S.hNM %>% subset(V_doses==1))+nrow(E.hNM %>% subset(V_doses==1))+nrow(A.hNM %>% subset(V_doses==1))+nrow(I.hNM %>% subset(V_doses==1))+nrow(R.hNM %>% subset(V_doses==1))
  +nrow(I.hM %>% subset(V_doses==1))+nrow(R.hM %>% subset(V_doses==1)) #完成pre接种
  vax_count2_h <- nrow(S.hNM %>% subset(V_doses==2))+nrow(E.hNM %>% subset(V_doses==2))+nrow(A.hNM %>% subset(V_doses==2))+nrow(I.hNM %>% subset(V_doses==2))+nrow(R.hNM %>% subset(V_doses==2))
  +nrow(I.hM %>% subset(V_doses==2))+nrow(R.hM %>% subset(V_doses==2)) #完成post接种
  
  vax_count1_l <- nrow(S.lNM %>% subset(V_doses==1))+nrow(E.lNM %>% subset(V_doses==1))+nrow(A.lNM %>% subset(V_doses==1))+nrow(I.lNM %>% subset(V_doses==1))+nrow(R.lNM %>% subset(V_doses==1))
  +nrow(I.lM %>% subset(V_doses==1))+nrow(R.lM %>% subset(V_doses==1)) #完成pre接种
  vax_count2_l <- nrow(S.lNM %>% subset(V_doses==2))+nrow(E.lNM %>% subset(V_doses==2))+nrow(A.lNM %>% subset(V_doses==2))+nrow(I.lNM %>% subset(V_doses==2))+nrow(R.lNM %>% subset(V_doses==2))
  +nrow(I.lM %>% subset(V_doses==2))+nrow(R.lM %>% subset(V_doses==2)) #完成post接种
  
  total_l <- total_l-9999
  
  Ns <- list("S.hNM"=S.hNM,
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
             "NewE_h"=0,
             "NewE_l"=0,
             "NewE_h_total"=0,
             "NewE_l_total"=0,
             "total"=total,
             "total_l"=total_l,
             "NewA_h_total"=0,
             "NewI_h_total"=0,
             "NewA_l_total"=0,
             "NewI_l_total"=0)
  
  return(Ns)
}
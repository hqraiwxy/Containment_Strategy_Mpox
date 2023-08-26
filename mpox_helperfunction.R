E_to_I <- function(df,parms){
  df %>% 
    subset(Days==Inc.pd) %>%
    dplyr::select(ID,V_doses,VE_s,VE_p,Inc.pd,Days) %>% 
    mutate(asympt = ifelse( V_doses == 2 & rbinom(length(ID),1,parms[["VE_p"]]), 1,0), #需要检查
           Inf.pd = Inc.pd+round(runif(1, parms[["gamma1"]], parms[["gamma2"]])),
           Inf.days = Inc.pd,
           ID.days = Inc.pd,
           removal.pd = NA) -> infected
  infected %>%
    subset(asympt==1) -> infected_asympt
  #保证输出的A结构不变
  infected_asympt %>%
    dplyr::select(ID,V_doses,VE_s,VE_p,Inc.pd,Inf.pd,Inf.days,ID.days,removal.pd) -> infected_asympt
  
  infected %>%
    subset(asympt==0)  -> infected_sympt
  
  #保证输出的I结构不变
  infected_sympt %>%
    dplyr::select(ID,V_doses,VE_s,VE_p,Inc.pd,Inf.pd,Inf.days,ID.days,removal.pd) -> infected_sympt
  
  df %>%
    subset(Inc.pd != Days) %>%
    mutate(Days=Days + 1) -> df
  
  list("infected_asympt"=infected_asympt,
       "infected_sympt"=infected_sympt,
       "df"=df)  #E
}

recover_I.M <- function(df,parms){
  
  df %>%
    subset(Inf.days == Inf.pd) %>%
    mutate(Rec.days = Inf.pd) %>%
    dplyr::select(ID,Rec.days,VE_s,VE_p,V_doses) -> recovered
  
  df %>%
    subset(!(ID %in% recovered$ID)) %>%
    mutate(Inf.days = Inf.days + 1) -> df  
  
  list("recovered"=recovered,
       "df"=df)
}

recover_or_test_h <- function(df,parms,symptoms){
  
  df %>%
    subset(Inf.days ==Inf.pd) %>% 
    mutate(Rec.days = Inf.pd) %>% 
    dplyr::select(ID,Rec.days,V_doses,VE_s,VE_p) -> recovered # recovered的人
  
  df %>%
    subset(!(ID %in% recovered$ID)) -> df # 去除掉recovered剩下的人
  
  df %>%  
    subset(ID.days == removal.pd) -> sym_onset # 症状出现被转走
  
  df %>% 
    subset(!(ID %in% sym_onset$ID) & !(ID %in% recovered$ID)) -> df 
  
  if(symptoms=="I"){
    df %>%
      subset(is.na(removal.pd) ) %>% 
      mutate(removal.pd = Inc.pd + 
                          round(runif(length(ID),parms[["id.I_lower"]],parms[["id.I_upper"]]))
             ) -> df1   ## symptomatic identified through symptoms
    # sym_onset中仅有20%被转走
    sym_onset_to_M <- sym_onset[sample(nrow(sym_onset),parms[["move_M_ratio_h"]]*nrow(sym_onset)),]
    sym_onset%>% 
      subset(!(ID %in% sym_onset_to_M$ID)) -> sym_onset_to_NM 
  }
  else{   # 对于无症状的，不会转到M队列，仍然放在NM队列中
    df %>%
      subset(is.na(removal.pd)) %>% 
      mutate(removal.pd = 10000 ) -> df1 #asymptomatic 不被移走 removal取一个大值
    
    sym_onset_to_NM <- sym_onset
    sym_onset%>% 
      subset(!(ID %in% sym_onset_to_NM$ID)) -> sym_onset_to_M
  }

  df %>%
    subset((removal.pd > ID.days ) & !(ID %in% df1$ID)) %>% 
    mutate(Inf.days = Inf.days + 1,
           ID.days = ID.days + 1
           ) -> df3   # 尚未表现出症状
  
  sym_onset_to_NM %>%
    mutate(Inf.days = Inf.days + 1
           # ID.days = ID.days + 1
    ) -> sym_onset_to_NM
  
  df1 %>%
    subset(ID.days != removal.pd) %>%
    mutate(ID.days = ID.days + 1,
           Inf.days = Inf.days + 1
           ) %>%
    bind_rows(df3) %>%
    bind_rows(sym_onset_to_NM) -> df
  
  list("recovered"=recovered,   # recovered  
       "sym_onset"=sym_onset_to_M,   # 出现症状后被移到M队列
       "df"=df)                 # 剩余的保持不变，但是本身时间和感染时间均加一
}




recover_or_test_l <- function(df,parms,symptoms){
  # df<- I.lNM  A.lNM
  df %>%
    subset(Inf.days ==Inf.pd) %>% 
    mutate(Rec.days = Inf.pd) %>% 
    dplyr::select(ID,Rec.days,V_doses,VE_s,VE_p) -> recovered # recovered的人
  
  df %>%
    subset(!(ID %in% recovered$ID)) -> df # 去除掉recovered剩下的人
  
  df %>%  
    subset(ID.days == removal.pd) -> sym_onset # 症状出现被转走
  
  df %>% 
    subset(!(ID %in% sym_onset$ID) & !(ID %in% recovered$ID)) -> df 
  
  if(symptoms=="I"){
    df %>%
      subset(is.na(removal.pd) ) %>% 
      mutate(removal.pd = Inc.pd + 
               round(runif(length(ID),parms[["id.I_lower"]],parms[["id.I_upper"]]))
      ) -> df1   ## symptomatic identified through symptoms
    
    # sym_onset中仅有move_M_ratio_l比例被转走
    sym_onset_to_M <- sym_onset[sample(nrow(sym_onset),parms[["move_M_ratio_l"]]*nrow(sym_onset)),]
    sym_onset%>% 
      subset(!(ID %in% sym_onset_to_M$ID)) -> sym_onset_to_NM 
  }
  else{   # 对于无症状的，不会转到M队列，仍然放在NM队列中
    df %>%
      subset(is.na(removal.pd)) %>% 
      mutate(removal.pd = 10000 ) -> df1 #asymptomatic 不被移走 removal取一个大值
    
    sym_onset_to_NM <- sym_onset
    sym_onset%>% 
      subset(!(ID %in% sym_onset_to_NM$ID)) -> sym_onset_to_M
  }
  
  df %>%
    subset((removal.pd > ID.days ) & !(ID %in% df1$ID)) %>% 
    mutate(Inf.days = Inf.days + 1
           # ID.days = ID.days + 1
    ) -> df3   # 尚未表现出症状
  
  sym_onset_to_NM %>%
    mutate(Inf.days = Inf.days + 1
           # ID.days = ID.days + 1
    ) -> sym_onset_to_NM
  
  df1 %>%
    subset(ID.days != removal.pd) %>%
    mutate(ID.days = ID.days + 1,
           Inf.days = Inf.days + 1
    ) %>%
    bind_rows(df3) %>%
    bind_rows(sym_onset_to_NM) -> df
  
  list("recovered"=recovered,   # recovered  
       "sym_onset"=sym_onset_to_M,   # 出现症状后被移到M队列
       "df"=df)                 # 剩余的保持不变，但是本身时间和感染时间均加一
}
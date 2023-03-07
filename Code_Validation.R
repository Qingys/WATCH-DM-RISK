###import packages
library(tidyverse)
library(lubridate)
library(survival)
library(mice)
library(riskRegression)
library(survminer)


### Main Cox model in original ----
###cox model
fit_cox <- coxph(
  Surv(survtime, survindicator) ~ 
    strata(sex)+
    sbp+age+BMI+
    log(FPG):FPG+
    I(FPG^2):I(age^0.5)+
    log(FPG):I(FPG^2):I(BMI^0.5)+
    CVD+family_diabetes+hypertension+smoke+drink,
  data = dt, x=TRUE)
summary(fit_cox)
###means
fit_cox$means
###means*Beta
sum(fit_cox$means*fit_cox$coefficients)
###baseline survival
predictCox(fit_cox, times = seq(3, 6, 1))


#### validation ----
Val_original <- Score(
  list('Cox model'=fit_cox),
  formula = Surv(survtime, survindicator) ~ 1,
  data = dt,
  summary = c("risks","IPA"),
  plots = c('Calibration'),
  times = c(3, 4, 5, 6),
  split.method = 'bootcv',
  B = 500,
  progress.bar = 3
)


##### calibration ----
cal_original_Brier <- Val_original$Brier$score %>% 
  filter(model=='Cox model') %>% 
  select(-se) %>% 
  rename(Brier.lb=lower, Brier.ub=upper)

cal_original_EO <- Val_original$Calibration$plotframe %>% 
  mutate(diff=(risk-pseudovalue)) %>% 
  nest_by(times) %>% 
  mutate(TotalEOratio=mean(data$risk)/mean(data$pseudovalue)) %>% 
  mutate(fit=list(
    t.test(data$diff) %>% broom::tidy()
  )) %>% 
  unnest(fit) %>% 
  select(-data, -statistic, -parameter, -alternative) %>% 
  rename(TotalEOdiff=estimate, 
         TotalEOdiff.pvalue=p.value, 
         TotalEOdiff.lb=conf.low,
         TotalEOdiff.ub=conf.high,
         TotalEOdiff.method=method)

cal_original_slope <- Val_original$Calibration$plotframe %>% 
  nest_by(times) %>% 
  mutate(fit=
           list(
             lm(pseudovalue ~ risk, data = data) %>% 
               broom::tidy(conf.int=T)
           )) %>% 
  unnest(fit) %>% 
  select(-data, -std.error, -statistic)

cal_original <- cal_original_Brier %>% 
  left_join(cal_original_EO, by='times') %>% 
  left_join(cal_original_slope, by='times')


##### calibration plots ----
fun_calibration <- function(timei, methodi, data=Val_original){
  calplot <- plotCalibration(data, 
                             time = timei, 
                             cens.method="local",
                             round = F,
                             method = methodi,
                             plot = F)
  if(methodi=='quantile'){
    calplot$plotFrames$`Cox model` %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'riskQuantile') %>% 
      mutate(riskGroup=1:10) %>% 
      mutate(time=timei) %>% 
      mutate(method=methodi) %>% 
      as_tibble() -> calplot_final
  } else {
    calplot$plotFrames$`Cox model` %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'riskPercentile') %>% 
      mutate(time=timei) %>% 
      mutate(method=methodi) %>% 
      as_tibble() -> calplot_final
  }
  return(calplot_final)
}

maplist <- list(
  timei=seq(3, 6, 1),
  methodi=rep('quantile', 4)
)
calplot_original_quantile <- pmap_dfr(maplist, fun_calibration)


###quantile plot
p_cal_quantile_original <- calplot_original_quantile %>% 
  mutate(time=paste0(time, '-years')) %>% 
  ggplot(aes(Obs, Pred)) + 
  geom_abline(intercept = 0, slope = 1, color='lightgrey', size=0.8) +
  geom_point(
    aes(color=time, shape=time),   ###multiple cal plot
    size=2.3, stroke=1.1, alpha=0.8) +
  #geom_line(aes(color=times), size=0.5, alpha=0.8, linetype='dashed') +
  scale_shape_manual(values=c(6,5,0,1), name='Time horizon') +  ###multiple cal plot
  scale_color_viridis_d(name='Time horizon') +  ###multiple cal plot
  labs(y='Predicted risk', x='Observed Kaplan-Meier risk',  
       title='Calibration plot of incident diabetes') +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw()

###E/O ratio by quantile
calplot_original_quantile %>% 
  mutate(EO=Pred/Obs)


### Validation in Tianjin data ----
#### import data ----
dt_tianjin <- readxl::read_excel('Tianjin_Data.xlsx')

###check missing data
dt_tianjin %>% 
  summarise(
    across(everything(), ~sum(is.na(.x)))
  )
dt_tianjin %>% 
  summarise(
    across(everything(), ~sum(is.na(.x))/length(.x))
  )


#### validation ----
Val_tianjin <- Score(
  list('Cox model'=fit_cox),
  formula = Surv(survtime, survindicator) ~ 1,
  data = dt_tianjin_imp,
  summary = c("risks","IPA"),
  plots = c('Calibration'),
  times = c(3, 4, 4.5)
)

##### calibration ----
cal_tianjin_Brier <- Val_tianjin$Brier$score %>% 
  filter(model=='Cox model') %>% 
  select(-se) %>% 
  rename(Brier.lb=lower, Brier.ub=upper)

cal_tianjin_EO <- Val_tianjin$Calibration$plotframe %>% 
  mutate(diff=(risk-pseudovalue)) %>% 
  nest_by(times) %>% 
  mutate(TotalEOratio=mean(data$risk)/mean(data$pseudovalue)) %>% 
  mutate(fit=list(
    t.test(data$diff) %>% broom::tidy()
  )) %>% 
  unnest(fit) %>% 
  select(-data, -statistic, -parameter, -alternative) %>% 
  rename(TotalEOdiff=estimate, 
         TotalEOdiff.pvalue=p.value, 
         TotalEOdiff.lb=conf.low,
         TotalEOdiff.ub=conf.high,
         TotalEOdiff.method=method)

cal_tianjin_slope <- Val_tianjin$Calibration$plotframe %>% 
  nest_by(times) %>% 
  mutate(fit=
           list(
             lm(pseudovalue ~ risk, data = data) %>% 
               broom::tidy(conf.int=T)
           )) %>% 
  unnest(fit) %>% 
  select(-data, -std.error, -statistic)

cal_tianjin <- cal_tianjin_Brier %>% 
  left_join(cal_tianjin_EO, by='times') %>% 
  left_join(cal_tianjin_slope, by='times')


##### calibration plots ----
fun_calibration <- function(timei, methodi, data=Val_tianjin){
  calplot <- plotCalibration(data, 
                             time = timei, 
                             cens.method="local",
                             round = F,
                             method = methodi,
                             plot = F)
  if(methodi=='quantile'){
  calplot$plotFrames$`Cox model` %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'riskQuantile') %>% 
    mutate(riskGroup=1:10) %>% 
    mutate(time=timei) %>% 
    mutate(method=methodi) %>% 
    as_tibble() -> calplot_final
  } else {
    calplot$plotFrames$`Cox model` %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'riskPercentile') %>% 
      mutate(time=timei) %>% 
      mutate(method=methodi) %>% 
      as_tibble() -> calplot_final
  }
  return(calplot_final)
}

maplist <- list(
  timei=c(3, 4, 4.5),
  methodi=rep('quantile', 3)
)
calplot_tianjin_quantile <- pmap_dfr(maplist, fun_calibration)


###quantile plot
p_cal_quantile_tianjin <- calplot_tianjin_quantile %>% 
  mutate(time=paste0(time, '-years')) %>% 
  ggplot(aes(Obs, Pred)) + 
  geom_abline(intercept = 0, slope = 1, color='lightgrey', size=0.8) +
  geom_point(
    aes(color=time, shape=time),   ###multiple cal plot
    size=2.3, stroke=1.1, alpha=0.8) +
  #geom_line(aes(color=times), size=0.5, alpha=0.8, linetype='dashed') +
  scale_shape_manual(values=c(6,5,0), name='Time horizon') +  ###multiple cal plot
  scale_color_viridis_d(name='Time horizon') +  ###multiple cal plot
  labs(y='Predicted risk', x='Observed Kaplan-Meier risk',  
       title='Calibration plot of incident diabetes') +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw()

###E/O by quantile
calplot_tianjin_quantile %>% 
  mutate(EOdiff=Pred-Obs,
         EOratio=Pred/Obs)


##### re-calibration plot ----
Lp_tianjin <- tibble(
  Lp=coxLP(fit_cox, data = dt_tianjin_imp, center = T),
) %>% bind_cols(dt_tianjin_imp)

cal_tianjin_Lp <- coxph(
  Surv(survtime, survindicator) ~ offset(Lp) + strata(sex), 
  x=T,
  data = Lp_tianjin)

fun_recal_plot <- function(timei, methodi){
  bsurv <- predictCox(cal_tianjin_Lp, times = timei)
  recal_tianjin <- Lp_tianjin %>% 
    mutate(resurv=
             case_when(
               sex=='male' ~ bsurv$survival[2]^exp(Lp),
               sex=='female' ~ bsurv$survival[1]^exp(Lp)
             )) %>% 
    mutate(timei=timei)
  reVal_tianjin <- Score(
    list('Recalibration model'=(1-recal_tianjin$resurv)),
    formula = Surv(survtime, survindicator) ~ 1,
    data = dt_tianjin_imp,
    summary = c("risks","IPA"),
    plots = c('Calibration'),
    times = timei
  )
  calplot <- plotCalibration(reVal_tianjin, 
                             time = timei, 
                             cens.method="local",
                             round = F,
                             method = methodi,
                             plot = F)
  if(methodi=='quantile'){
    calplot$plotFrames$`Recalibration model` %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'riskQuantile') %>% 
      mutate(riskGroup=1:10) %>% 
      mutate(time=timei) %>% 
      mutate(method=methodi) %>% 
      as_tibble() -> calplot_final
  } else {
    calplot$plotFrames$`Recalibration model` %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'riskPercentile') %>% 
      mutate(time=timei) %>% 
      mutate(method=methodi) %>% 
      as_tibble() -> calplot_final
  }
  return(calplot_final)
}

maplist <- list(
  timei=c(3, 4, 4.5),
  methodi=rep('quantile', 3)
)
recalplot_tianjin_quantile <- pmap_dfr(maplist, fun_recal_plot)


###quantile plot
p_recal_quantile_tianjin <- recalplot_tianjin_quantile %>% 
  mutate(time=paste0(time, '-years')) %>% 
  ggplot(aes(Obs, Pred)) + 
  geom_abline(intercept = 0, slope = 1, color='lightgrey', size=0.8) +
  geom_point(
    aes(color=time, shape=time),   ###multiple cal plot
    size=2.3, stroke=1.1, alpha=0.8) +
  #geom_line(aes(color=times), size=0.5, alpha=0.8, linetype='dashed') +
  scale_shape_manual(values=c(6,5,0), name='Time horizon') +  ###multiple cal plot
  scale_color_viridis_d(name='Time horizon') +  ###multiple cal plot
  labs(y='Predicted risk', x='Observed Kaplan-Meier risk',  
       title='Recalibration plot of incident diabetes') +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw()


###E/O by quantile
recalplot_tianjin_quantile %>% 
  mutate(EOdiff=Pred-Obs,
         EOratio=Pred/Obs)



##### re-calibration Total EO and slope ----
fun_recal <- function(timei){
  bsurv <- predictCox(cal_tianjin_Lp, times = timei)
  recal_tianjin <- Lp_tianjin %>% 
    mutate(resurv=
             case_when(
               sex=='male' ~ bsurv$survival[2]^exp(Lp),
               sex=='female' ~ bsurv$survival[1]^exp(Lp)
             )) %>% 
    mutate(timei=timei)
  reVal_tianjin <- Score(
    list('Recalibration model'=(1-recal_tianjin$resurv)),
    formula = Surv(survtime, survindicator) ~ 1,
    data = dt_tianjin_imp,
    summary = c("risks","IPA"),
    plots = c('Calibration'),
    times = timei
  )
  
  recal_tianjin_Brier <- reVal_tianjin$Brier$score %>% 
    filter(model=='Recalibration model') %>% 
    select(-se) %>% 
    rename(Brier.lb=lower, Brier.ub=upper)
  
  recal_tianjin_EO <- reVal_tianjin$Calibration$plotframe %>% 
    mutate(diff=(risk-pseudovalue)) %>% 
    nest_by(times) %>% 
    mutate(TotalEOratio=mean(data$risk)/mean(data$pseudovalue)) %>% 
    mutate(fit=list(
      t.test(data$diff) %>% broom::tidy()
    )) %>% 
    unnest(fit) %>% 
    select(-data, -statistic, -parameter, -alternative) %>% 
    rename(TotalEOdiff=estimate, 
           TotalEOdiff.pvalue=p.value, 
           TotalEOdiff.lb=conf.low,
           TotalEOdiff.ub=conf.high,
           TotalEOdiff.method=method)
  
  recal_tianjin_slope <- reVal_tianjin$Calibration$plotframe %>% 
    nest_by(times) %>% 
    mutate(fit=
             list(
               lm(pseudovalue ~ risk, data = data) %>% 
                 broom::tidy(conf.int=T)
             )) %>% 
    unnest(fit) %>% 
    select(-data, -std.error, -statistic)
  
  recal_tianjin <- recal_tianjin_Brier %>% 
    left_join(recal_tianjin_EO, by='times') %>% 
    left_join(recal_tianjin_slope, by='times')

  return(recal_tianjin)
}

maplist <- list(
  timei=c(3, 4, 4.5)
)

recal_tianjin <- pmap_dfr(maplist, fun_recal)



### Validation in RCHealthCare ----
#### import data ----
dt_RCHealthCare <- read_csv('RCHealthCare_Data.csv')


###check missing data
dt_RCHealthCare %>% 
  summarise(
    across(everything(), ~sum(is.na(.x)))
  )
dt_RCHealthCare %>% 
  summarise(
    across(everything(), ~sum(is.na(.x))/length(.x))
  )


#### validation ----
Val_RCHealthCare <- Score(
  list('Cox model'=fit_cox),
  formula = Surv(survtime, survindicator) ~ 1,
  data = dt_RCHealthCare_imp,
  summary = c("risks","IPA"),
  plots = c('Calibration'),
  times = c(3, 4, 5)
)


##### calibration ----
cal_RCHealthCare_Brier <- Val_RCHealthCare$Brier$score %>% 
  filter(model=='Cox model') %>% 
  select(-se) %>% 
  rename(Brier.lb=lower, Brier.ub=upper)

cal_RCHealthCare_EO <- Val_RCHealthCare$Calibration$plotframe %>% 
  mutate(diff=(risk-pseudovalue)) %>% 
  nest_by(times) %>% 
  mutate(TotalEOratio=mean(data$risk)/mean(data$pseudovalue)) %>% 
  mutate(fit=list(
    t.test(data$diff) %>% broom::tidy()
  )) %>% 
  unnest(fit) %>% 
  select(-data, -statistic, -parameter, -alternative) %>% 
  rename(TotalEOdiff=estimate, 
         TotalEOdiff.pvalue=p.value, 
         TotalEOdiff.lb=conf.low,
         TotalEOdiff.ub=conf.high,
         TotalEOdiff.method=method)

cal_RCHealthCare_slope <- Val_RCHealthCare$Calibration$plotframe %>% 
  nest_by(times) %>% 
  mutate(fit=
           list(
             lm(pseudovalue ~ risk, data = data) %>% 
               broom::tidy(conf.int=T)
           )) %>% 
  unnest(fit) %>% 
  select(-data, -std.error, -statistic)

cal_RCHealthCare <- cal_RCHealthCare_Brier %>% 
  left_join(cal_RCHealthCare_EO, by='times') %>% 
  left_join(cal_RCHealthCare_slope, by='times')



##### calibration plots ----
fun_calibration <- function(timei, methodi, data=Val_RCHealthCare){
  calplot <- plotCalibration(data, 
                             time = timei, 
                             cens.method="local",
                             round = F,
                             method = methodi,
                             plot = F)
  if(methodi=='quantile'){
    calplot$plotFrames$`Cox model` %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'riskQuantile') %>% 
      mutate(riskGroup=1:10) %>% 
      mutate(time=timei) %>% 
      mutate(method=methodi) %>% 
      as_tibble() -> calplot_final
  } else {
    calplot$plotFrames$`Cox model` %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'riskPercentile') %>% 
      mutate(time=timei) %>% 
      mutate(method=methodi) %>% 
      as_tibble() -> calplot_final
  }
  return(calplot_final)
}

maplist <- list(
  timei=seq(3, 5, 1),
  methodi=rep('quantile', 3)
)
calplot_RCHealthCare_quantile <- pmap_dfr(maplist, fun_calibration)


###quantile plot
p_cal_quantile_RCHealthCare <- calplot_RCHealthCare_quantile %>% 
  mutate(time=paste0(time, '-years')) %>% 
  ggplot(aes(Obs, Pred)) + 
  geom_abline(intercept = 0, slope = 1, color='lightgrey', size=0.8) +
  geom_point(
    aes(color=time, shape=time),   ###multiple cal plot
    size=2.3, stroke=1.1, alpha=0.8) +
  #geom_line(aes(color=times), size=0.5, alpha=0.8, linetype='dashed') +
  scale_shape_manual(values=c(6,5,0), name='Time horizon') +  ###multiple cal plot
  scale_color_viridis_d(name='Time horizon') +  ###multiple cal plot
  labs(y='Predicted risk', x='Observed Kaplan-Meier risk',  
       title='Calibration plot of incident diabetes') +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw()


###E/O by quantile
calplot_RCHealthCare_quantile %>% 
  mutate(EOdiff=Pred-Obs,
         EOratio=Pred/Obs)


##### re-calibration plot ----
Lp_RCHealthCare <- tibble(
  Lp=coxLP(fit_cox, data = dt_RCHealthCare_imp, center = T),
) %>% bind_cols(dt_RCHealthCare_imp)

cal_RCHealthCare_Lp <- coxph(
  Surv(survtime, survindicator) ~ offset(Lp) + strata(sex), 
  x=T,
  data = Lp_RCHealthCare)

fun_recal <- function(timei, methodi){
  bsurv <- predictCox(cal_RCHealthCare_Lp, times = timei)
  recal_RCHealthCare <- Lp_RCHealthCare %>% 
    mutate(resurv=
             case_when(
               sex=='male' ~ bsurv$survival[2]^exp(Lp),
               sex=='female' ~ bsurv$survival[1]^exp(Lp)
             )) %>% 
    mutate(timei=timei)
  reVal_RCHealthCare <- Score(
    list('Recalibration model'=(1-recal_RCHealthCare$resurv)),
    formula = Surv(survtime, survindicator) ~ 1,
    data = dt_RCHealthCare_imp,
    summary = c("risks","IPA"),
    plots = c('Calibration'),
    times = timei
  )
  calplot <- plotCalibration(reVal_RCHealthCare, 
                             time = timei, 
                             cens.method="local",
                             round = F,
                             method = methodi,
                             plot = F)
  if(methodi=='quantile'){
    calplot$plotFrames$`Recalibration model` %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'riskQuantile') %>% 
      mutate(riskGroup=1:10) %>% 
      mutate(time=timei) %>% 
      mutate(method=methodi) %>% 
      as_tibble() -> calplot_final
  } else {
    calplot$plotFrames$`Recalibration model` %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'riskPercentile') %>% 
      mutate(time=timei) %>% 
      mutate(method=methodi) %>% 
      as_tibble() -> calplot_final
  }
  return(calplot_final)
}

maplist <- list(
  timei=seq(3, 5, 1),
  methodi=rep('quantile', 3)
)
recalplot_RCHealthCare_quantile <- pmap_dfr(maplist, fun_recal)


###quantile plot
p_recal_quantile_RCHealthCare <- recalplot_RCHealthCare_quantile %>% 
  mutate(time=paste0(time, '-years')) %>% 
  ggplot(aes(Obs, Pred)) + 
  geom_abline(intercept = 0, slope = 1, color='lightgrey', size=0.8) +
  geom_point(
    aes(color=time, shape=time),   ###multiple cal plot
    size=2.3, stroke=1.1, alpha=0.8) +
  #geom_line(aes(color=times), size=0.5, alpha=0.8, linetype='dashed') +
  scale_shape_manual(values=c(6,5,0), name='Time horizon') +  ###multiple cal plot
  scale_color_viridis_d(name='Time horizon') +  ###multiple cal plot
  labs(y='Predicted risk', x='Observed Kaplan-Meier risk',  
       title='Recalibration plot of incident diabetes') +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw()


###E/O ratio by quantile
recalplot_RCHealthCare_quantile %>% 
  mutate(EOdiff=Pred-Obs,
         EOratio=Pred/Obs)



##### re-calibration Total EO and slope ----
fun_recal <- function(timei){
  bsurv <- predictCox(cal_RCHealthCare_Lp, times = timei)
  recal_RCHealthCare <- Lp_RCHealthCare %>% 
    mutate(resurv=
             case_when(
               sex=='male' ~ bsurv$survival[2]^exp(Lp),
               sex=='female' ~ bsurv$survival[1]^exp(Lp)
             )) %>% 
    mutate(timei=timei)
  reVal_RCHealthCare <- Score(
    list('Recalibration model'=(1-recal_RCHealthCare$resurv)),
    formula = Surv(survtime, survindicator) ~ 1,
    data = dt_RCHealthCare_imp,
    summary = c("risks","IPA"),
    plots = c('Calibration'),
    times = timei
  )
  
  recal_RCHealthCare_Brier <- reVal_RCHealthCare$Brier$score %>% 
    filter(model=='Recalibration model') %>% 
    select(-se) %>% 
    rename(Brier.lb=lower, Brier.ub=upper)
  
  recal_RCHealthCare_EO <- reVal_RCHealthCare$Calibration$plotframe %>% 
    mutate(diff=(risk-pseudovalue)) %>% 
    nest_by(times) %>% 
    mutate(TotalEOratio=mean(data$risk)/mean(data$pseudovalue)) %>% 
    mutate(fit=list(
      t.test(data$diff) %>% broom::tidy()
    )) %>% 
    unnest(fit) %>% 
    select(-data, -statistic, -parameter, -alternative) %>% 
    rename(TotalEOdiff=estimate, 
           TotalEOdiff.pvalue=p.value, 
           TotalEOdiff.lb=conf.low,
           TotalEOdiff.ub=conf.high,
           TotalEOdiff.method=method)
  
  recal_RCHealthCare_slope <- reVal_RCHealthCare$Calibration$plotframe %>% 
    nest_by(times) %>% 
    mutate(fit=
             list(
               lm(pseudovalue ~ risk, data = data) %>% 
                 broom::tidy(conf.int=T)
             )) %>% 
    unnest(fit) %>% 
    select(-data, -std.error, -statistic)
  
  recal_RCHealthCare <- recal_RCHealthCare_Brier %>% 
    left_join(recal_RCHealthCare_EO, by='times') %>% 
    left_join(recal_RCHealthCare_slope, by='times')
  
  return(recal_RCHealthCare)
}

maplist <- list(
  timei=seq(3, 5, 1)
)

recal_RCHealthCare <- pmap_dfr(maplist, fun_recal)


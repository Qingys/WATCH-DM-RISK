###import packages
library(tidyverse)
library(tidytable)
library(data.table)
library(broom)
library(survival)
library(riskRegression)


######omitting tidy process
###...

### Imputation ----------------------
library(mice)
library(miceadds)

dt %>%
  group_by(id) %>% 
  mutate(group_id=cur_group_id()) %>% 
  ungroup() %>%
  select(3:4, 7:18, 20:41, 46) -> dt_miss

predM <- make.predictorMatrix(dt_miss)

predM[, 'group_id'] <- -2
predM['group_id', 'group_id'] <- 0

dt_imp <- mice(dt_miss, m=1, method = '2l.pmm', predictorMatrix = predM, maxit = 10)

dt_imp <- complete(dt_imp)

###bind data
dt <- dt %>% 
  select(-c(3:4, 7:18, 20:41)) %>% 
  bind_cols(dt_imp)

dt %>% is.na() %>% sum()

dt <- dt[index==1]


### Boruta --------------------------------
### Feature selection
library(Boruta)

dt <- dt_index %>%
  select(-c(1:7, 11:23, 46:50))

###for total
fit_Boruta <- Boruta(
  x=dt[, -c('survtime', 'survindicator')],
  y=Surv(dt$survtime, dt$survindicator),
  doTrace = 3,
  getImp = function(x, y) getImpRfZ(x, y, 
                                    num.threads=64, num.trees = 128, splitrule='maxstat')
)
###selected feature
Boruta_selectedVars_total <- attStats(fit_Boruta) %>% 
  rownames_to_column(var = 'Vars') %>% arrange(-meanImp) %>% 
  mutate(population='total')


### Random forest for interaction -----
library(ranger)
library(randomForestExplainer)

dt_index %>% 
  mutate(
    across(c('drink', 'smoke', 'hypertension', 'CVD', 'family_diabetes', 'sex'), 
           as_factor)
  ) -> dt_index


###build model
###sex+FPG+age+BMI+sbp
###drink+smoke+hypertension+CVD+family_diabetes
form <- as.formula(Surv(survtime, survindicator) ~ sex+FPG+age+BMI+sbp+
                     drink+smoke+hypertension+CVD+family_diabetes)

fit_rf <- ranger(formula = form, 
                 num.trees = 500,
                 data = dt_index,
                 splitrule = 'maxstat',
                 importance='impurity',
                 num.threads = 32
                 )

###interaction
interactions <- min_depth_interactions(fit_rf)
interactions %>%
  filter(occurrences>400) %>% 
  arrange(mean_min_depth)


### Boosted model for non-linear ------
library(mboost)

dt_index %>% 
  mutate(
    across(c('drink', 'smoke', 'hypertension', 'CVD', 'family_diabetes', 'sex'), 
           as_factor)
  ) -> dt_index


###build model
###form non-linear
form <- as.formula(Surv(survtime, survindicator) ~ FP(FPG,p=c(-1,-0.5,0.5,1,2))*(
  FP(age,p=c(-1,-0.5,0.5,1,2))+
    FP(sbp,p=c(-1,-0.5,0.5,1,2))+
    FP(BMI,p=c(-1,-0.5,0.5,1,2))
  ) +
    drink+smoke+hypertension+CVD+family_diabetes +strata(sex) -1)

###models
fit_mboost <- glmboost(form, data = dt_index, family = CoxPH(),
                       control = boost_control(mstop = 100, trace = T))

###stable selection
fit_ss <- stabsel(fit_mboost, cutoff = 0.8, PFER = 1, mc.cores=16)

result_ss_total <- fit_ss$phat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'Vars') %>% 
  pivot_longer(-Vars, names_to = 'iters') %>% 
  mutate(iters=
           str_extract(iters, '[:digit:]+'),
         iters=as.numeric(iters)) %>% 
  mutate(population='total')



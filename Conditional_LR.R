#CONDITIONAL logistic regression FID.
library("tidyverse")
library("epitools")
library("ggplot2")
library("survival")

#CLR models. 
#NSAID + covariates.
conditional_1 <- clogit(Backpn_3m_v1_v2 ~ NSAID_v0 + AGE + SEX + group + count_pain_sites_acute_v0 + strata(FID), data=backpain_chronicized, method = 'efron')
summary(conditional_1)

#PARACETAMOL + covariates.
conditional_2 <- clogit(Backpn_3m_v1_v2 ~ PARACETAMOL + AGE + SEX + group + count_pain_sites_acute_v0 + strata(FID), data=backpain_chronicized, method = 'efron')
summary(conditional_2)

#OPIODS + covariates
conditional_3 <- clogit(Backpn_3m_v1_v2 ~ OPIOIDS + AGE + SEX + group + count_pain_sites_acute_v0 + strata(FID), data=backpain_chronicized, method = 'efron')
summary(conditional_3)

#ANTIDEPRESSANTS + covariates
conditional_4 <- clogit(Backpn_3m_v1_v2 ~ ANTIDEPRESSANTS + AGE + SEX + group + count_pain_sites_acute_v0 + strata(FID), data=backpain_chronicized, method = 'efron')
summary(conditional_4)

#GABAPENTIN/PREGABALIN + covariates.
conditional_5 <- clogit(Backpn_3m_v1_v2 ~ gabapentin_pregabalin + AGE + SEX + group + count_pain_sites_acute_v0 + strata(FID), data=backpain_chronicized, method = 'efron')
summary(conditional_5)

#ALL MEDICATION IN 1 MODEL. 
conditional_6 <- clogit(Backpn_3m_v1_v2 ~ NSAID_v0 + AGE + SEX + group + count_pain_sites_acute_v0 + PARACETAMOL + ANTIDEPRESSANTS + OPIOIDS + gabapentin_pregabalin + strata(FID), data=backpain_chronicized, method = 'efron')
summary(conditional_6)





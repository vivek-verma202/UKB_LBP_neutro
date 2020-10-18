options(device = "RStudioGD", digits = 4)
df <- readRDS("./backpain_chronicized.rds")
sapply(df, function(x)(sum(is.na(x))))
library(tidyverse)
# select acute back pain IDs:
df1 <- df  %>% filter(Back_pn_v0 == "1") # lost 372410 IDs
df2 <- df1 %>% filter(Back_pn_3.m_v0 == "0") # lost 89553 IDs

# combine follow-up visits
sum(is.na(df2$Back_pn_3.m_v1)) # 40102 out of 40531
sum(is.na(df2$Backpn_3m_v1)) # 39215 out of 40531
# if v1 missing fill with v2
df2$Backpn_3m_v1[is.na(df2$Backpn_3m_v1)] <-
    df2$Backpn_3m_v2[is.na(df2$Backpn_3m_v1)]
sum(is.na(df2$Backpn_3m_v1)) # 37907 out of 40531

# select acute back pain IDs with follow-up data:
df3 <- df2[complete.cases(df2$Backpn_3m_v1),]
names(df3)[9] <- "chronic_back_pain"
summary(factor(df3$chronic_back_pain))
#    0    1
# 2183  441

# fix type from int to factor
df3$group[df3$group == "WHITE"] <- 0
df3$group[df3$group == "ASIAN"] <- 1
df3$group[df3$group == "BLACK"] <- 2
df3$group[df3$group == "MIXED"] <- 3
df3$group[df3$group == "OTHER"] <- 4

df3 <- within(df3, {
    chronic_back_pain <- as.factor(chronic_back_pain)
    ANTIDEPRESSANTS <- as.factor(ANTIDEPRESSANTS)
    NSAID_v0_clean <- as.factor(NSAID_v0_clean)
    OPIOIDS <- as.factor(OPIOIDS)
    PARACETAMOL <- as.factor(PARACETAMOL)
    SEX <- as.factor(SEX)
    group <- as.factor(group)
})

library(survival)
# could not find pregabalin / gabapentin in the dataset
mod6 <- clogit(chronic_back_pain ~ AGE + SEX +
                   group + count_pain_sites_acute_v0 +
                   NSAID_v0_clean + PARACETAMOL + OPIOIDS +
                   ANTIDEPRESSANTS + strata(FID),
               data = df3, method = "breslow")
# fails:
# Error in coxph(formula = Surv(rep(1, 2624L), chronic_back_pain) ~ AGE +  :
#                   an id statement is required for multi-state models
sum(is.na(df3$FID)) # 2524 (96.19% data is missing)

# simple logistic regression:
GLM.1 <- glm(chronic_back_pain ~ AGE + SEX + group + ANTIDEPRESSANTS
             + NSAID_v0_clean + OPIOIDS + PARACETAMOL +
                 count_pain_sites_acute_v0,
             family=binomial(logit), data=df3)
summary(GLM.1)
exp(coef(GLM.1))  # Exponentiated coefficients ("odds ratios")

# without count_pain_sites:
GLM.2 <- glm(chronic_back_pain ~ AGE + SEX + group + ANTIDEPRESSANTS
             + NSAID_v0_clean + OPIOIDS + PARACETAMOL,
             family=binomial(logit), data=df3)
summary(GLM.2)
exp(coef(GLM.2))

# adding WBC data
df <- readRDS("./data/df.RDS")
df <- df[,c(1,3,6:11)]
df3 <- merge(df3,df,by="ID",all.x = T)

# sanity check
plot(df$count_pain_sites_acute_v0,df$pain_count_1m)

# model with neutrophils
GLM.3 <- glm(chronic_back_pain ~ AGE + SEX + group + ANTIDEPRESSANTS
             + NSAID_v0_clean + OPIOIDS + PARACETAMOL + neutro_p_v0,
             family=binomial(logit), data=df3)
summary(GLM.3)
exp(coef(GLM.3))

GLM.4 <- glm(chronic_back_pain ~ AGE + SEX + group + ANTIDEPRESSANTS
            + NSAID_v0_clean + OPIOIDS + PARACETAMOL + neutro_p_v0 +
                count_pain_sites_acute_v0,
            family=binomial(logit), data=df3)
summary(GLM.4)
exp(coef(GLM.4))
df3$FID_fac <- as.factor(df3$FID)
write.csv(data.frame(table(df3$FID_fac)),"./data/FID.csv",
            row.names = F, quote = F)


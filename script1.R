options(device = "RStudioGD", digits = 2)
df <- readRDS("./data/backpain_chronicized-2.rds")
sapply(df, function(x)(sum(is.na(x))))

df <- within(df, {
    Backpn_3m_v1_v2       <- as.factor(Backpn_3m_v1_v2)
    ANTIDEPRESSANTS       <- as.factor(ANTIDEPRESSANTS)
    NSAID_v0              <- as.factor(NSAID_v0)
    OPIOIDS               <- as.factor(OPIOIDS)
    PARACETAMOL           <- as.factor(PARACETAMOL)
    gabapentin_pregabalin <- as.factor(gabapentin_pregabalin)
    SEX                   <- as.factor(SEX)
    group                 <- as.factor(group)
})

df$group <- relevel(df$group, ref = "WHITE")

library(questionr)
library(interactions)
library(jtools)
# model 6
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group,
               family=binomial(logit), data=df))

# model 5
odds.ratio(glm(Backpn_3m_v1_v2 ~ gabapentin_pregabalin + AGE + SEX + group,
               family=binomial(logit), data=df))

# model 4
odds.ratio(glm(Backpn_3m_v1_v2 ~ ANTIDEPRESSANTS +
                + AGE + SEX + group,
               family=binomial(logit), data=df))

# get wbc data
wbc <- readRDS("./data/UKB_WBC_all_raw.RDS")
df <- merge(df,wbc,by="ID",all.x = T)

library(Rcmdr)

odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group,
               family=binomial(logit), data=df))


odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0*lympho_p_v1 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group + WBC_v1,
               family=binomial(logit), data=df))


odds.ratio(glm(Backpn_3m_v1_v2 ~ AGE + SEX + group + WBC_v1
 #+ baso_p_v1
 #+ eosino_p_v1
 #+ neutro_p_v1
 #+ lympho_p_v1
 #+ mono_p_v1
               ,family=binomial(logit), data=df))


odds.ratio(glm(Backpn_3m_v1_v2 ~ AGE + SEX + group + NSAID_v0
                   + lympho_p_v1
               ,family=binomial(logit), data=df))


fit <- glm(Backpn_3m_v1_v2 ~ AGE + SEX + group +
               + lympho_p_v1 : NSAID_v0
           ,family=binomial(logit), data=df)
summ(fit)
odds.ratio(fit)
interact_plot(fit, pred = lympho_p_v1, modx = NSAID_v0,
              interval = T, line.thickness = 2)



with(wbc, cor.test(lympho_c_v1, neutro_c_v1, alternative="two.sided",
                   method="pearson"))

library(Rcmdr)






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
library(mediation)
# model 6
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group + neutro_p_v1,
               family=binomial(logit), data=df))

summ(glm(Backpn_3m_v1_v2 ~ NSAID_v0 + neutro_p_v1 + AGE + SEX + group,
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



odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group
               #+ baso_p_v1
               #+ eosino_p_v1
               #+ neutro_p_v1
               #+ lympho_p_v1
               #+ mono_p_v1
               ,family=binomial(logit), data=df))


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


fit <- glm(Backpn_3m_v1_v2 ~ AGE + SEX + group + NSAID_v0+lympho_p_v1
           ,family=binomial(logit), data=df)
summary(fit)
round(odds.ratio(fit),2)
interact_plot(fit, pred = neutro_p_v1, modx = NSAID_v0,
              interval = T, line.thickness = 2,
              x.label = "Neutrophil %", y.label = "Probability of chronic back pain development")



with(wbc, cor.test(lympho_c_v1, neutro_c_v1, alternative="two.sided",
                   method="pearson"))

library(Rcmdr)


odds.ratio(glm(neutro_p_v1 ~ AGE + SEX + group
               ,family=binomial(logit), data=df))


summary(lm(neutro_p_v1 ~ AGE + SEX + group, data=df))

# mediation
med.fit <- lm(lympho_p_v1 ~ NSAID_v0 + AGE + SEX + group, data = df)
out.fit <- glm(Backpn_3m_v1_v2 ~ lympho_p_v1 + NSAID_v0 + AGE + SEX + group,
               data = df, family = binomial("probit"))
med.out <- mediate(med.fit, out.fit, treat = "NSAID_v0",
                   mediator = "lympho_p_v1",
                   robustSE = TRUE, sims = 100)
summary(med.out)
test.TMint(med.out, conf.level = .95)
sens.out <- medsens(med.out, rho.by = 0.1, effect.type = "indirect",
                    sims = 100)
plot(sens.out, sens.par = "rho", main = "neutrophil %", ylim = c(-0.1, 0.1),
     xlim = c(-1,1))

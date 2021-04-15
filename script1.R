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
# remove relatives:
df <- df[-c(918,1951),]
df$time_interval <- ifelse(is.na(df$time_interval_years_v1),
                           df$time_interval_years_v2,
                           df$time_interval_years_v1)
library(questionr)
library(interactions)
library(jtools)

# model 1-6
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0
               + AGE + SEX + group + time_interval,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~ PARACETAMOL
               + AGE + SEX + group + time_interval,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~ OPIOIDS
               + AGE + SEX + group + time_interval,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~ ANTIDEPRESSANTS
               + AGE + SEX + group + time_interval,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~  gabapentin_pregabalin
               + AGE + SEX + group + time_interval,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group + time_interval,
               family=binomial(logit), data=df))

# model 7
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group + time_interval
               + neutro_p_v1,
               family=binomial(logit), data=df))

# supplementary 1-4

odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group + time_interval
               + baso_p_v1,
               family=binomial(logit), data=df))
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group + time_interval
               + eosino_p_v1,
               family=binomial(logit), data=df))
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group + time_interval
               + lympho_p_v1,
               family=binomial(logit), data=df))
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group + time_interval
               + mono_p_v1,
               family=binomial(logit), data=df))

# sup 5-8
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   AGE + SEX + group + time_interval
               + neutro_p_v1,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   AGE + SEX + group + time_interval
               + lympho_p_v1,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~
                   AGE + SEX + group + time_interval
               + neutro_p_v1:NSAID_v0,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~
                   AGE + SEX + group + time_interval
               + lympho_p_v1:NSAID_v0,
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
wbc <- readRDS("./data/WBC.RDS")
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


fit <- glm(Backpn_3m_v1_v2 ~ AGE.x + SEX + group + NSAID_v0*neutro_p_v1
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
saveRDS(df,"df.RDS")

# On 14 APR 2021
df1 <- read.delim("C:/Users/Vivek/Desktop/UKB_LBP_neutro/depression_anxiety.tsv")
df1 <- df1[,c(1,2,5)]
names(df1)[1] <- "ID"
df <- merge(df,df1,by="ID",all.x = T)
df$psych_distress1 <- df$depression_v0 + df$anxiety_panic_v0
df$psych_distress  <- ifelse(df$psych_distress1 > 0, 1, 0)
df$psych_distress1 <- NULL
df$psych_distress <- as.factor(df$psych_distress)
saveRDS(df,"df.RDS")


# models
round(odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                         PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                         gabapentin_pregabalin
               + AGE + SEX + group + time_interval + psych_distress + mono_p_v1,
               family=binomial(logit), data=df)),2)

# sup 5-8
round(odds.ratio(glm(Backpn_3m_v1_v2 ~
                   AGE + SEX + group + time_interval+ psych_distress
               + neutro_p_v1:NSAID_v0,
               family=binomial(logit), data=df)),2)

odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   AGE + SEX + group + time_interval
               + lympho_p_v1,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~
                   AGE + SEX + group + time_interval
               + neutro_p_v1:NSAID_v0,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~
                   AGE + SEX + group + time_interval
               + lympho_p_v1:NSAID_v0,
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



l <- c(6:10,17)
for (i in l) {
    print(names(df)[i])
    print(table(df[,64]))
}


# APR 15
df <- readRDS("C:/Users/vverma3/Desktop/Repos/UKB_LBP_neutro/data/df.RDS")
UKB_pain <- readRDS("./data/UKB_pain.RDS")

df <- merge(df,UKB_pain,by="ID",all.x = T)
rm(UKB_pain)
saveRDS(df,"./data/df.RDS")


library(questionr)
# models 1 -5
drugs <- c(6,8,10,9,17)
for (i in drugs) {
    print(odds.ratio(glm(as.formula(paste0("Backpn_3m_v1_v2"," ~ ",names(df)[i],
                                     "+AGE+SEX+group+time_interval+psych_distress+pain_count_1m")),
                   family=binomial(logit), data=df)))
}

# model 6
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group +
                   time_interval + psych_distress + pain_count_1m,
               family=binomial(logit), data=df))

# model 7
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   PARACETAMOL + OPIOIDS + ANTIDEPRESSANTS +
                   gabapentin_pregabalin + AGE + SEX + group +
                   time_interval + psych_distress + pain_count_1m +
                   neutro_p_v1,
               family=binomial(logit), data=df))

# sup models 1 -4
cells <- c(50,47,38,41)

for (i in cells) {
    print(odds.ratio(glm(as.formula(paste0("Backpn_3m_v1_v2"," ~ ",
    "NSAID_v0+PARACETAMOL+OPIOIDS+ANTIDEPRESSANTS+gabapentin_pregabalin",
    "+AGE+SEX+group+time_interval+psych_distress+pain_count_1m+",
                                           names(df)[i])),
                         family=binomial(logit), data=df)))
}


# sup 5-6
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                         AGE + SEX + group + time_interval+ psych_distress
                     + pain_count_1m + neutro_p_v1,
                     family=binomial(logit), data=df))
odds.ratio(glm(Backpn_3m_v1_v2 ~ NSAID_v0 +
                   AGE + SEX + group + time_interval+ psych_distress
               + pain_count_1m + lympho_p_v1,
               family=binomial(logit), data=df))


# sup 7-8
odds.ratio(glm(Backpn_3m_v1_v2 ~
                   AGE + SEX + group + time_interval+ psych_distress
               + pain_count_1m + neutro_p_v1:NSAID_v0,
               family=binomial(logit), data=df))

odds.ratio(glm(Backpn_3m_v1_v2 ~
                   AGE + SEX + group + time_interval+ psych_distress
               + pain_count_1m + lympho_p_v1:NSAID_v0,
               family=binomial(logit), data=df))















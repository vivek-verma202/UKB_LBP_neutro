options(device = "RStudioGD", digits = 4, verbose = T)
# prepare df:
## get pain pheno
UKB_pain <- readRDS("I:/Intralab Communications/Audrey-Vivek/UKB_pain.RDS")
UKB_pain <- UKB_pain[,c(1,5,37,43,48)]
UKB_cov <- readRDS("S:/UKB/UKB-pheno-R/UKB_cov.RDS")
UKB_cov <- UKB_cov[,c(1,2,4)]
names(UKB_cov) <- c("ID","age","sex")

## get covs


## add ethnicity
eth <- read.table("./data/Eth_v0_nona.csv", sep = ",")
eth <- eth[-1,c(1,4)]
names(eth) <- c("ID","ethnicity")

## make UKB_WBC data
library(readr)
UKB <- read_csv("S:/win7_desktop/UKB_WBC/UKB_followup/data/UKBBall.csv")
UKB <- UKB[,c(1:34)]
saveRDS(UKB,"./data/UKB_WBC_all_raw.RDS")
# ------ cleaning and imputation, not performed <start>
#Checking prop miss
propmiss <- function(dataframe) {
    m <- sapply(dataframe, function(x) {
        data.frame(
            nmiss=sum(is.na(x)),
            n=length(x),
            propmiss=sum(is.na(x))/length(x)
        )
    })
    d <- data.frame(t(m))
    d <- sapply(d, unlist)
    d <- as.data.frame(d)
    d$variable <- row.names(d)
    row.names(d) <- NULL
    d <- cbind(d[ncol(d)],d[-ncol(d)])
    return(d[order(d$propmiss), ])
}
propmiss(UKB)
hist(UKB$lympho_p_v1)
UKB$lympho_p_v1[is.na(UKB$lympho_p_v1)] <-
    UKB$lympho_c_v1[is.na(UKB$lympho_p_v1)]*100/
    UKB$WBC_v1[is.na(UKB$lympho_p_v1)]
UKB$mono_p_v1[is.na(UKB$mono_p_v1)] <-
    UKB$mono_c_v1[is.na(UKB$mono_p_v1)]*100/
    UKB$WBC_v1[is.na(UKB$mono_p_v1)]
UKB$neutro_p_v1[is.na(UKB$neutro_p_v1)] <-
    UKB$neutro_c_v1[is.na(UKB$neutro_p_v1)]*100/
    UKB$WBC_v1[is.na(UKB$neutro_p_v1)]
UKB$baso_p_v1[is.na(UKB$baso_p_v1)] <-
    UKB$baso_c_v1[is.na(UKB$baso_p_v1)]*100/
    UKB$WBC_v1[is.na(UKB$baso_p_v1)]
UKB$eosino_p_v1[is.na(UKB$eosino_p_v1)] <-
    UKB$eosino_c_v1[is.na(UKB$eosino_p_v1)]*100/
    UKB$WBC_v1[is.na(UKB$eosino_p_v1)]



UKB$WBC_v2[is.na(UKB$WBC_v2)] <- UKB$WBC_v3[is.na(UKB$WBC_v2)]
UKB$lympho_c_v2[is.na(UKB$lympho_c_v2)] <- UKB$lympho_c_v3[is.na(UKB$lympho_c_v2)]
UKB$mono_c_v2[is.na(UKB$mono_c_v2)] <- UKB$mono_c_v3[is.na(UKB$mono_c_v2)]
UKB$neutro_c_v2[is.na(UKB$neutro_c_v2)] <- UKB$neutro_c_v3[is.na(UKB$neutro_c_v2)]
UKB$eosino_c_v2[is.na(UKB$eosino_c_v2)] <- UKB$eosino_c_v3[is.na(UKB$eosino_c_v2)]
UKB$baso_c_v2[is.na(UKB$baso_c_v2)] <- UKB$baso_c_v3[is.na(UKB$baso_c_v2)]
UKB$lympho_p_v2[is.na(UKB$lympho_p_v2)] <- UKB$lympho_p_v3[is.na(UKB$lympho_p_v2)]
UKB$mono_p_v2[is.na(UKB$mono_p_v2)] <- UKB$mono_p_v3[is.na(UKB$mono_p_v2)]
UKB$neutro_p_v2[is.na(UKB$neutro_p_v2)] <- UKB$neutro_p_v3[is.na(UKB$neutro_p_v2)]
UKB$eosino_p_v2[is.na(UKB$eosino_p_v2)] <- UKB$eosino_p_v3[is.na(UKB$eosino_p_v2)]
UKB$baso_p_v2[is.na(UKB$baso_p_v2)] <- UKB$baso_p_v3[is.na(UKB$baso_p_v2)]

UKB <- UKB[,c(1:3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33)]
names(UKB) <- sub("*_v1", "_v0", names(UKB))
names(UKB) <- sub("*_v2", "_v1", names(UKB))
# ------ cleaning and imputation, not performed <end>
UKB$WBC_v2[is.na(UKB$WBC_v2)] <- UKB$WBC_v3[is.na(UKB$WBC_v2)]
UKB$neutro_c_v2[is.na(UKB$neutro_c_v2)] <-
    UKB$neutro_c_v3[is.na(UKB$neutro_c_v2)]
UKB$neutro_p_v2[is.na(UKB$neutro_p_v2)] <-
    UKB$neutro_p_v3[is.na(UKB$neutro_p_v2)]
UKB <- UKB[,c(1:3,11,12,26,27)]
names(UKB) <- c("ID","WBC_v0","WBC_v1","neutro_c_v0","neutro_c_v1",
                "neutro_p_v0","neutro_p_v1")

df1 <- merge(UKB_pain,UKB,by="ID",all = T)
df2 <- merge(UKB_cov,eth,by="ID",all = T)
df  <- merge(df1,df2,by="ID",all = T)
df  <- df[-c(1:14),]
row.names(df) <- NULL
df$sex <- as.factor(df$sex)
df$ethnicity <- as.factor(df$ethnicity)
saveRDS(df,"./data/df.RDS")

# summary table
library(arsenal)
df <- readRDS("./data/df.RDS")
table(df$back_pain_1m)
# FALSE   TRUE
# 197152 130122
df <- df[complete.cases(df$back_pain_1m),]
df <- df[df$back_pain_1m == TRUE,]
summary(df$chrn_back_pain)
df <- df[complete.cases(df$chrn_back_pain),c(1,3:14)]
df1 <- df[,-1]
table_one <- tableby(chrn_back_pain ~ ., data = df1)
summary(table_one)

library(ggplot2)
ggplot(df, aes(x = WBC_v0, fill = chrn_back_pain)) +
    geom_density(alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = c("darkgreen","red")) +
    theme_classic() +
    xlab("WBC (10^9 cells/L)") +
    theme(plot.title = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x = element_text(color = "black", size = 15),
          axis.title.x = element_text(color = "black", size = 15),
          legend.position = "none")
ggplot(df, aes(x = neutro_c_v0, fill = chrn_back_pain)) +
    geom_density(alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = c("darkgreen","red")) +
    theme_classic() +
    xlab("Neutrophils (10^9 cells/L)") +
    theme(plot.title = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x = element_text(color = "black", size = 15),
          axis.title = element_text(color = "black", size = 15),
          legend.position = "none")
ggplot(df, aes(x = neutro_p_v0, fill = chrn_back_pain)) +
    geom_density(alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = c("darkgreen","red")) +
    theme_classic() +
    xlab("Neutrophils %") +
    theme(plot.title = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x = element_text(color = "black", size = 15),
          axis.title = element_text(color = "black", size = 15),
          legend.position = "none")

# get drug data
df1 <- read_delim("I:/Intralab Communications/Jonathan-vivek/drugs.melt",
                 "\t", escape_double = FALSE, col_types = cols(df = col_skip(),
                                                               itemnum = col_skip()), trim_ws = TRUE)
df1 <- df1[df1$instnum == 0,]
codes <- read.csv("./data/NSAID_codes.csv")
codes <- codes$ï..Code
df1 <- df1[df1$value %in% codes, ]

df$nsaid <- as.logical(ifelse(df$ID %in% df1$IID, TRUE, FALSE))
saveRDS(df,"./data/df1.RDS")

ggplot(df, aes(x = neutro_c_v0, fill = nsaid)) +
    geom_density(alpha = 0.4, size = 1.5) +
    scale_fill_manual(values = c("#A674A5", "#4638D9")) +
    theme_classic() +
    xlab("WBC (10^9 cells/L)") +
    theme(plot.title = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x = element_text(color = "black", size = 15),
          axis.title = element_text(color = "black", size = 15),
          legend.position = "none")
df1 <- df[,-1]
table_one <- tableby(nsaid ~ ., data = df1)
summary(table_one)
# modelling
df <- readRDS("./data/df1.RDS")

with(df, plotMeans(neutro_p_v0, as.factor(pain_count_chrn), error.bars="se", connect=TRUE))


summary(glm(formula = chrn_back_pain ~ age + sex + ethnicity + neutro_c_v0
            , family = binomial(logit), data = df))

# updated analysis
backpain_chronicized <- readRDS("./backpain_chronicized.rds")
backpain_chronicized$acute_only_v0 <- NULL
backpain_chronicized$acute_only_v0[(backpain_chronicized$Back_pn_v0 == "1")] <- 1
backpain_chronicized$acute_only_v0[(backpain_chronicized$Back_pn_3.m_v0 != "0")] <- NA
df <- backpain_chronicized[!is.na(backpain_chronicized$acute_only_v0),]
df$missing_followup <- NULL
df$missing_followup <- as.factor(is.na(df$Back_pn_3.m_v1) + is.na(df$Back_pn_3.m_v2))


dim(df[(!is.na(df$Back_pn_3.m_v1) | !is.na(df$Back_pn_3.m_v2)), ]





#
#
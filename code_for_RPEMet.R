

#加载环境和数据 ####
rm(list = ls())
library(haven)
data <- read_dta("无OCT的86014人.dta")
rownames(data) <- data$eid_ageing

#定义结局
data$t2d_time <- data$t2d_date_2104 - data$startdate
data$t2d_time <- ifelse(is.na(data$t2d_time) & data$t2d_2104 == 0 & !is.na(data$deathdate), data$deathdate - data$startdate, data$t2d_time)
data$t2d_time <- ifelse(is.na(data$t2d_time) & data$t2d_2104 == 0 & !is.na(data$lostdate), data$lostdate - data$startdate, data$t2d_time)
data$t2d_time <- ifelse(is.na(data$t2d_time), data$dateend_2104 - data$startdate, data$t2d_time) #t2d_time若为缺失值，则以dateend-startdate填补
data$t2d_time <- ifelse(data$t2d_time < 0, NA, data$t2d_time) #t2d_time小于0的替换为缺失值
data$t2d_2104 <- ifelse(is.na(data$t2d_time), NA, data$t2d_2104) #t2d_time为缺失值的行中的t2d_2104也替换为缺失值
table(data$t2d_2104)
summary(data$t2d_time/365)

#拆分数据集
set.seed(666)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3)) ##7：3的比例生成随机数，命名为ind
train_set <- data[ind == 1, ] #训练数据
test_set <- data[ind == 2, ] #测试数据

#train set
traindata <- cbind(train_set[c("eid_ageing","t2d_time","t2d_2104")],train_set[,2:267],
                   train_set[c("baselineage","whr","hyperlipidemia","q8")])
traindata <- na.omit(traindata)
traindata <- na.omit(traindata)
table(train_set$t2d_2104)

#test data
testdata <- cbind(test_set[c("eid_ageing","t2d_time","t2d_2104")],test_set[,2:267],
                  test_set[c("baselineage","whr","hyperlipidemia","q8")])
testdata <- na.omit(testdata)
table(test_set$t2d_2104)


x <- as.matrix(traindata[,met_filtered$Met])
library(survival)
y <- data.matrix(Surv(traindata$t2d_time, traindata$t2d_2104))
library(glmnet)
fit <- glmnet(x, y, family = "cox", maxit = 200000, alpha=1)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 200000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoMet = row.names(coef)[index]
lassoMet = c("t2d_time","t2d_2104",lassoMet)
#展示筛选结果
lassoMet
lassoSigMet=traindata[,lassoMet]
paste(lassoMet,collapse="+")
#m4+m6+m8+m13+m17+m21+m29+m32+m39+m48+m55+m56+m59+m79+m80+m81+m83+m89+m91+m92+m93+m94+m98+m99+m100+m102+m103+m104+m119+m126+m133+m135+m137+m138+m139+m169+m170+m217+m222+m245

library(survival)
multiCox=coxph(Surv(t2d_time, t2d_2104) ~ 
                 m4+m6+m8+m9+m13+m14+m17+m18+m21+m25+m26+m29+m30+m32+m39+
                 m48+m55+m56+m59+m77+m78+m79+m80+m81+m82+m83+m84+m85+m86+
                 m87+m88+m89+m90+m91+m92+m93+m94+m96+m97+m98+m99+m100+m101+
                 m102+m103+m104+m119+m126+m127+m128+m129+m130+m133+m134+m135+
                 m136+m137+m138+m139+m169+m170+m217+m222+m245, data = train_set)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox_纳入所有.xls",sep="\t",row.names=F,quote=F)


riskScore=predict(multiCox,type="risk",newdata=test_set)
coxMet=rownames(multiCoxSum$coefficients)
coxMet=gsub("`","",coxMet)
outCol=c("t2d_time","t2d_2104",coxMet)

#分为四分位数
quantile_risk=as.vector(ifelse(riskScore>=quantile(riskScore,0.75),"Quantile 1",
                               ifelse(riskScore<quantile(riskScore,0.25),"Quantile 4",
                                      ifelse(riskScore>=quantile(riskScore,0.25) & riskScore<quantile(riskScore,0.5),"Quantile 3","Quantile 2"))))
table(quantile_risk)

write.table(cbind(id=rownames(cbind(test_set[,outCol],riskScore,risk)),cbind(test_set[,outCol],riskScore,quantile_risk)),
            file="risk_纳入所有.txt",
            sep="\t",
            quote=F,
            row.names=F)



------------------------------------------------------------------------------------------------------------
  #  Kaplan-Meier生存分析-----------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------
  
library(survival)
library(survminer)
library(ggthemes)
rt=read.table("risk.txt",header=T,sep="\t")
diff=survdiff(Surv(t2d_time, t2d_2104) ~ quantile_risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,3)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(t2d_time, t2d_2104) ~ quantile_risk, data = rt)
p<-ggsurvplot(fit, 
              data=rt,
              legend.labs=c("Quantile 1", "Quantile 2","Quantile 3","Quantile 4"),
              legend.title="T2D",
              xlab="Time(days)",
              pval=paste0("P=",pValue),
              pval.size=4,
              pval.coord=c(0,0.2),
              break.time.by = 730,
              censor.size = 2,
              conf.int = T,
              conf.int.style = "ribbon",
              fun = "cumhaz",
              surv.scale="percent",
              palette = "npg",
              ggtheme = theme_few())
p

pdf(file="Kaplan-Meier曲线_纳入所有.pdf", onefile = FALSE, width = 8, height =8)
p
dev.off()

#计算累计发病率比值RR
cal<-rt[c(colnames(rt)[3],colnames(rt)[2],"quantile_risk")]
rt<-na.omit(rt)
a1<-subset(rt,rt$quintile_risk=="Quantile 1")
a3<-subset(rt,rt$quintile_risk=="Quantile 4")
a1[,3]<-factor(a1[,3])
a3[,3]<-factor(a3[,3])
aa<-summary(a1[,3])
a<-aa[2] #top25%发病数
b<-sum(a1[,2])/365 #top25%人年
cc<-summary(a3[,3])
c<-cc[2] #bottom25%发病数
d<-sum(a3[,2])/365 #bottom25%人年
e1=a/b*1000
e2=c/d*1000 # c为千人年发病率
RR=e1/e2 #相对发病率
d=log(RR)
lower=exp(d - 1.96*sqrt((1/a)+(1/c)))####上限
upper=exp(d + 1.96*sqrt((1/a)+(1/c)))####下限
cat(cat(cat(RR, lower, sep=" ("), upper, sep="-"), ")", sep="")


------------------------------------------------------------------------------------------------------------
  #  ROC曲线------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------
  
library(survivalROC)
library(survival)
library(pROC)
library(ggplot2)

#传统
cox_m1 <- coxph(Surv(t2d_time,t2d_2104) ~ a1+a10+whr+a13+a6+a7+a2+a8+q8+hyperlipidemia+a12, data = traindata)
cox_m11<-step(cox_m1,direction = "both")
testdata$risk_score1<-predict(cox_m11,type="risk",newdata=testdata)
roc_t2d_1 <- roc(t2d_2104 ~ risk_score1, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T,
                 smooth = T)
#加入代谢物
cox_m2 <- coxph(Surv(t2d_time,t2d_2104) ~ a1+a10+whr+a13+a6+a7+a2+a8+q8+hyperlipidemia+a12+
                  m4+m6+m8+m9+m13+m14+m17+m18+m21+m25+m26+m29+m30+m32+m39+
                  m48+m55+m56+m59+m77+m78+m79+m80+m81+m82+m83+m84+m85+m86+
                  m87+m88+m89+m90+m91+m92+m93+m94+m96+m97+m98+m99+m100+m101+
                  m102+m103+m104+m119+m126+m127+m128+m129+m130+m133+m134+m135+
                  m136+m137+m138+m139+m169+m170+m217+m222+m245, data = traindata)
cox_m22<-step(cox_m2,direction = "both")
testdata$risk_score2<-predict(cox_m22,type="risk",newdata=testdata)
roc_t2d_2 <- roc(t2d_2104 ~ risk_score2, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T,
                 smooth = T)

#单独代谢物
cox_m3 <- coxph(Surv(t2d_time,t2d_2104) ~ 
                  m4+m6+m8+m9+m13+m14+m17+m18+m21+m25+m26+m29+m30+m32+m39+
                  m48+m55+m56+m59+m77+m78+m79+m80+m81+m82+m83+m84+m85+m86+
                  m87+m88+m89+m90+m91+m92+m93+m94+m96+m97+m98+m99+m100+m101+
                  m102+m103+m104+m119+m126+m127+m128+m129+m130+m133+m134+m135+
                  m136+m137+m138+m139+m169+m170+m217+m222+m245, data = traindata)
cox_m33<-step(cox_m3,direction = "both")
testdata$risk_score3<-predict(cox_m3,type="risk",newdata=testdata)
roc_t2d_3 <- roc(t2d_2104 ~ risk_score3, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T,
                 smooth = T)

#Delong's test
roc.test(roc_t2d_2,roc_t2d_1) #D = 7.1411, boot.n = 2000, boot.stratified = 1, p-value = 9.258e-13
roc.test(roc_t2d_3,roc_t2d_1)


#合并三曲线图
p_roc_3<-ggroc(list(roc_t2d_1, roc_t2d_3, roc_t2d_2),legacy.axes = TRUE,alpha=1,linetype=1,size=1) + ggtitle("Type 2 diabetes") +
  theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               colour = "darkgrey",
               linetype = "dashed") +
  xlab("1-Specifility") + ylab("Sensitivity") +
  theme(legend.title=element_blank()) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = c(1,0),legend.justification = c(1,0))+
  scale_colour_discrete(labels = c("Clinical indicators","RPET metabolic state","Clinical indicators + RPET metabolic state"))
p_roc_3

pdf(file="AUC_t2d三模型smooth_纳入所有.pdf", onefile = FALSE, width = 8, height =8)
p_roc_3
dev.off()

#合并双曲线图
p_roc_3<-ggroc(list(roc_t2d_1,  roc_t2d_2),legacy.axes = TRUE,alpha=1,linetype=1,size=1) + ggtitle("Type 2 diabetes") +
  theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               colour = "darkgrey",
               linetype = "dashed") +
  xlab("1-Specifility") + ylab("Sensitivity") +
  theme(legend.title=element_blank()) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = c(1,0),legend.justification = c(1,0))+
  scale_color_manual(values = c("#F8736A","#659EFF"),
                     labels = c("Clinical indicators","Clinical indicators + RPET metabolic state"))
p_roc_3

pdf(file="AUC_t2d三模型smooth_纳入所有.pdf", onefile = FALSE, width = 8, height =8)
p_roc_3
dev.off()

#建立其他单独危险因素曲线
roc_riskScore <- roc(t2d_2104 ~ risk_score3, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
reverse_t2d_2104 <- chartr('01', '10', testdata$t2d_2104)
roc_a1 <- roc(t2d_2104 ~ a1, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
roc_a10 <- roc(t2d_2104 ~ a10, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
roc_whr <- roc(t2d_2104 ~ whr, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
roc_a13 <- roc(t2d_2104 ~ a13, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
roc_a6 <- roc(t2d_2104 ~ a6, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
roc_a7 <- roc(t2d_2104 ~ a7, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
roc_a2 <- roc(t2d_2104 ~ a2, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
roc_a8 <- roc(t2d_2104 ~ a8, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
roc_q8 <- roc(t2d_2104 ~ q8, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
roc_hyperlipidemia <- roc(reverse_t2d_2104 ~ hyperlipidemia, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)
roc_a12 <- roc(t2d_2104 ~ a12, data = testdata, print.thres=TRUE, print.auc=TRUE,ci =T, plot=T)

p_roc_12<- ggroc(list(roc_a1, roc_a2, roc_a8, roc_a6, roc_a7, roc_a10, roc_whr, roc_q8, roc_hyperlipidemia, roc_a13, roc_a12, roc_riskScore),legacy.axes = TRUE,alpha=1,linetype=5,size=1) + 
  ggtitle("Type 2 diabetes") +
  theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               colour = "darkgrey",
               linetype = "dashed") + 
  xlab("1-Specifility") + ylab("Sensitivity") +
  theme(legend.title=element_blank()) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = c(1,0),legend.justification = c(1,0))+
  scale_colour_discrete(labels = c("Age","Sex","Ethnicity","Smoking","Drinking","BMI","WHR","Hypertension","Dyslipidemia","BP-lowering medication","Lipid-lowering medication","RPET metabolic state"))+
  scale_linetype_manual(labels = c("dashed", "solid"), values = c(5,5,5,5,5,5,5,5,5,5,5,1))
p_roc_12

pdf(file="AUC_t2d多模型_纳入所有.pdf",onefile = FALSE,
    width = 8,
    height =8)
p_roc_12
dev.off()

roc.test(roc_riskScore,roc_a10) #Z = 12.734, p-value < 2.2e-16
roc.test(roc_riskScore,roc_whr) #Z = 12.678, p-value < 2.2e-16
roc.test(roc_riskScore,roc_a13) #Z = 20.676, p-value < 2.2e-16
roc.test(roc_riskScore,roc_q8) #Z = 19.371, p-value < 2.2e-16
roc.test(roc_riskScore,roc_hyperlipidemia) #D = 20.679, boot.n = 2000, boot.stratified = 1, p-value < 2.2e-16
roc.test(roc_riskScore,roc_a12) #Z = 17.625, p-value < 2.2e-16

------------------------------------------------------------------------------------------------------------
  #  IDI和NRI-----------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------
  
#IDI & NRI
library(survival)
library(survIDINRI)
mat1 <- as.matrix(testdata[,c("a1","a10","whr","a13","a6","a7","a2","a8","q8","a12")])
mat2 <- as.matrix(testdata[,c("a1","a10","whr","a13","a6","a7","a2",
                              "a8","q8","a12","m4","m8","m9","m13","m14",
                              "m17","m18","m25","m29","m32","m39","m48","m55","m59","m78",
                              "m79","m80","m81","m83","m84","m85","m87","m88","m89","m90",
                              "m92","m93","m94","m96","m97","m100","m101","m102","m103",
                              "m104","m127","m129","m130","m134","m135","m138","m169")])
x<-IDI.INF(testdata[,2:3],mat1, mat2, t=365*12, npert=1000)
IDI.INF.OUT(x)

pdf(file="IDI_纳入所有.pdf",onefile = FALSE,
    width = 8,
    height =8)
IDI.INF.GRAPH(x)
dev.off()


------------------------------------------------------------------------------------------------------------
  #    Hosmer-Lemeshow检验----------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------

library(ResourceSelection)
yhat2 <- predict(cox_m2, newdata=testdata, type = "lp")
hoslem.test(x=testdata$t2d_2104, y=yhat2, g=10)

------------------------------------------------------------------------------------------------------------
  #    DCA曲线----------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------
  
gc()
library(survival)
library(ggDCA)
library(ggprism)
library(ggthemes)
library(ggplot2)
testdata<-na.omit(testdata) #drop NA

#训练生存模型
model1 <- coxph(Surv(t2d_time,t2d_2104) ~ a1 + a10 + whr + a13 + 
                  a6 + a7 + a2 + a8 + q8 + a12, data = traindata)
model2 <- coxph(Surv(t2d_time,t2d_2104) ~ a1 + a10 + whr + a13 + a6 + a7 + a2 + 
                  a8 + q8 + hyperlipidemia + a12 + m4 + m8 + m9 + m13 + m14 + 
                  m17 + m18 + m25 + m29 + m32 + m39 + m48 + m55 + m59 + m78 + 
                  m79 + m80 + m81 + m83 + m84 + m85 + m87 + m88 + m89 + m90 + 
                  m92 + m93 + m94 + m96 + m97 + m100 + m101 + m102 + m103 + 
                  m104 + m127 + m129 + m130 + m134 + m135 + m138 + m169, data = traindata)

dca_coxph <- dca(model1,model2,
                 model.names = c("Clinical indicators","Clinical indicators + RPET metabolic state"),
                 new.data = testdata)

p_dca<-ggplot(dca_combined, aes(thresholds, NB, color = model, linetype = model, lwd = model)) +
  scale_x_continuous(limits = c(0, 0.9),guide = "prism_minor") +
  scale_y_continuous(limits = c(-0.075/5, 0.075),guide = "prism_minor")+
  scale_color_manual(values = c("#F8736A","#659EFF", "#CAD2C5", "#CAD2C5"))+
  scale_linetype_manual(values = c(1,1,4,3))+
  scale_linewidth_manual(values=c(1,1,1,1))+
  geom_line() +
  theme_few() +
  theme(legend.position = c(1,1),legend.justification = c(1,1)) +
  theme(legend.title=element_blank()) +
  xlab("Risk Threshold") +
  ylab("Net Benefit") +
  ggtitle("Decision curve")
p_dca

pdf(file="DCA曲线_纳入所有.pdf", onefile = FALSE, width = 8, height =8)
p_dca
dev.off()

------------------------------------------------------------------------------------------------------------
  #    拼图-------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------
  
library(patchwork)

pdf(file="T2D汇总图_2×纳入所有.pdf",onefile = FALSE,
    width = 33.87/1.1*0.4,
    height =13.5/1.1*0.5)
p_roc_3 + p_dca +
  plot_layout(ncol=2,nrow=1)
dev.off()


#———————————————————————————————————————————————————————————————————————————————————————————————————————————####
#     基线对比----------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------
  
data$dataset <- ind
data$dataset[ind == 1] <- 1
data$dataset[ind == 2] <- 2
table(data$dataset)
#排除基线糖尿病
data$t2d_time <- data$t2d_date_2104 - data$startdate
data$t2d_time <- ifelse(is.na(data$t2d_time), data$dateend_2104 - data$startdate, data$t2d_time) #t2d_time若为缺失值，则以dateend-startdate填补
data<-subset(data,data$t2d_time>0) #排除基线糖尿病

data_7824 <- read_dta("/Volumes/杨少鹏/杨少鹏/010 文章5-RNFL& GCIPL的代谢组/数据/phase 1 GCIPL相关代谢物/OCT and Metabolomics_YSP.dta") 
data_7824$dataset <- 0
library(plyr)
data_merge <- rbind.fill(data,data_7824)
table(data_merge$dataset)

##合并基线年龄和随访时间
alldata <- read_dta("/Volumes/杨少鹏/杨少鹏/010 文章5-RNFL& GCIPL的代谢组/数据/UKB event_ww20220418.dta")
age_baseline <- alldata[c("eid_ageing","age_baseline")]
data_merge <- merge(data_merge, age_baseline, by = "eid_ageing", all.x = TRUE)
alldata$followup_time <- as.Date("2021-04-28") - alldata$startdate
followup_time <- alldata[c("eid_ageing","followup_time")]
data_merge <- merge(data_merge, followup_time, by = "eid_ageing", all.x = TRUE)

mean(data_merge$age_baseline) #基线年龄
sd(data_merge$age_baseline)
table(data_merge$a2) #性别
mean(data_merge$followup_time)/365 #随访时间
sd(data_merge$followup_time)/365
quantile(data_merge$followup_time, c(0.25,0.75))/365 #四分位数


#table 1
library(tableone)
library(haven)

#OCT人群、训练集、验证集对比
dput(names(data_merge))
allVars <-c("a1", "a2", "a8", "a5","a9","a4","a10", "a6", "a7","o1","o3","q8", "hyperlipidemia","a12","a13") #所有变量定义为allVars
fvars<-c("a1", "a2", "a8", "a5","a9","a4","a10", "a6", "a7", "q8", "hyperlipidemia","a12","a13") #分类变量定义为fvars
tab2 <- CreateTableOne(vars = allVars, data = data_merge, factorVars=fvars, strata = "dataset", addOverall = T )
tab3<- print(tab2,  quote = FALSE, noSpaces = T, printToggle = F, showAllLevels = T, addOverall = T)
print(tab3)
write.csv(tab3, file = "总人群基线对比.csv")

#训练集与验证集对比
dput(names(data))
allVars <-c("a1", "a2", "a8", "a5","a9","a4","a10", "a6", "a7","o1","o3","q8", "hyperlipidemia","a12","a13") #所有变量定义为allVars
fvars<-c("a1", "a2", "a8", "a5","a9","a4","a10", "a6", "a7", "q8", "hyperlipidemia","a12","a13") #分类变量定义为fvars
tab2 <- CreateTableOne(vars = allVars, data = data, factorVars=fvars, strata = "dataset", addOverall = T )
tab3<- print(tab2,  quote = FALSE, noSpaces = T, printToggle = F, showAllLevels = T, addOverall = T)
print(tab3)
write.csv(tab3, file = "训练验证基线对比.csv")



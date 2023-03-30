# Title: Collective Forest Land Rights Facilitate Cooperative Behavior. 
# Authors: Komal Preet Kaur, Kimberlee Chang, Krister Andersson
# Published in Conservation Letters, March 2023


#* load libraries ----
library(pacman)
p_load(haven, visreg, tidyverse, clubSandwich, sandwich, gdata, xtable,ggplot2, lmtest, foreign, gtools, reshape2, lfe,stargazer, stm, tm, SnowballC, tibble, mediation,readr, arm, interplot, multiwayvcov, miceaddss, corrplot, Hmisc, readstata13, xlsx, RItools, optmatch, estimatr, coin, lme4, brms,nbpMatching, rbounds, experiment, splines, optimx,sensitivitymv, sensitivitymw, sensitivitymult, optimx)


#* sessionInfo() ----
rm(list=ls())
set.seed(123)

#* Load data ----
dat <- read_dta("ifri.dta")

dat <- as.data.frame(dat)
dat$UHHNUM <- as.numeric(dat$UHHNUM)
dat$REFOWNLAND_DUM <- as.factor(dat$REFOWNLAND_DUM)

#* Make DV in percentages
dat$SOCAP_OUT_MEAN <- ((dat$SOCAP_OUT_MEAN-1)/3)*100
dat$SOCAP_IN_MEAN <- ((dat$SOCAP_IN_MEAN-1)/3)*100

# Table A-1
dat1 <- dat %>% 
  dplyr::select(CID,SOCAP_IN_MEAN, SOCAP_OUT_MEAN,REFOWNLAND_DUM,GRPYR,UHHNUM,UWEALTHDIF, ULITPROP, UPROPFEM) %>% 
  drop_na()

table(dat1$CID, dat1$REFOWNLAND_DUM) 

###############################################################
#* base models for forest-related interactions  (Models 1 and 2 in table 2) ----
main1 <- lm(SOCAP_IN_MEAN ~ REFOWNLAND_DUM + GRPYR + UHHNUM +UWEALTHDIF + ULITPROP + UPROPFEM, data = dat)   ##just land dejure rights
main2 <- lm(SOCAP_IN_MEAN ~ RDEJURE_DUM + GRPYR + UHHNUM + UWEALTHDIF + ULITPROP + UPROPFEM, data = dat)   ##just product dejure rights


mods <- list(main1,main2)
for(i in 1:2){
  assign(paste("cl.robust.se.",i,sep=""), 
         sqrt(diag(cluster.vcov(mods[[i]], dat$CID)))) # cluster-robust SEs for ols
}

stargazer(
  main1, main2,
  se=list(cl.robust.se.1, cl.robust.se.2),
  type = "text", star.char = c("+", "*", "**", "***"), star.cutoffs = c(.1, .05, .01, .001),
  notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), notes.append=FALSE)

#* base models  (Models 3 and 4 in table 2) ----
main1 <- lm(SOCAP_OUT_MEAN ~ REFOWNLAND_DUM + GRPYR + UHHNUM +UWEALTHDIF + ULITPROP + UPROPFEM, data = dat)   ##just land dejure rights
main2 <- lm(SOCAP_OUT_MEAN ~ RDEJURE_DUM + GRPYR + UHHNUM + UWEALTHDIF + ULITPROP + UPROPFEM, data = dat)   ##just product dejure rights

mods <- list(main1,main2)
for(i in 1:2){
  assign(paste("cl.robust.se.",i,sep=""), 
         sqrt(diag(cluster.vcov(mods[[i]], dat$CID)))) # cluster-robust SEs for ols
}

stargazer(
  main1, main2,
  se=list(cl.robust.se.1, cl.robust.se.2),
  digits = 2,
  type = "text", star.char = c("+", "*", "**", "***"), star.cutoffs = c(.1, .05, .01, .001),
  notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), notes.append=FALSE)



##################### Robustness checks #########################
#* 1. Accounting for country variation (Model 1 in Table A-6) ----
cvar <- dat %>% 
  filter(CID %in% c("BOL", "GUA", "MEX", "IND", "UGA")) %>% 
  dplyr::select(SOCAP_OUT_MEAN,SOCAP_IN_MEAN,REFOWNLAND_DUM, GRPYR,UHHNUM,UWEALTHDIF,ULITPROP,UPROPFEM, CID) %>% 
  tidyr::drop_na()

table(cvar$CID, cvar$REFOWNLAND_DUM)

m1_cvar <- lm(SOCAP_OUT_MEAN ~ REFOWNLAND_DUM + GRPYR + UHHNUM +UWEALTHDIF + ULITPROP + UPROPFEM, data = cvar)   

#* 2. Accounting for country variation (Model 2 in Table A-6) ----
m2_cvar <- lm(SOCAP_OUT_MEAN ~ REFOWNLAND_DUM + GRPYR + UHHNUM +UWEALTHDIF + ULITPROP + UPROPFEM + CID, data = cvar)   

stargazer(m1_cvar, m2_cvar, out = "cvar.doc", 
          digits = 2,
          align=TRUE, omit.stat=c("LL","ser","f"), no.space=TRUE,
          covariate.labels = c("De jure land rights", "User group age", "User group size", "User group wealth difference (dummy)", "User group literacy rate", "Proportion of women in user group"), 
          notes.align = "l", notes.append = F,
          star.char = c("+", "*", "**", "***"), star.cutoffs = c(.1, .05, .01, .001), type="text",
          notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), 
          se = starprep(m1_cvar, m2_cvar,  clusters = cvar$CID,se_type="stata"))

#* 3. Accounting for country history (Model 3 in Table A-6) ----
#** creating a variable for designating country's history of de jure land rights. bottom-up countries: India, Uganda coded as 0. Top-down: Bolivia, Guatemala, Mexico coded as 1  ----
dat$chist <- NA
dat$chist <- ifelse(dat$CID %in% c("BOL", "GUA", "MEX"), 1, 0)

chist_dat <- dat %>% 
  dplyr::select(SOCAP_OUT_MEAN,SOCAP_IN_MEAN,REFOWNLAND_DUM, GRPYR,UHHNUM,UWEALTHDIF,ULITPROP,UPROPFEM, CID) %>% 
  mutate(chist = ifelse(dat$CID %in% c("BOL", "GUA", "MEX"), 1, 0)) %>% 
  drop_na() 

m1_chist <- lm(SOCAP_OUT_MEAN ~ REFOWNLAND_DUM*chist + GRPYR + UHHNUM +UWEALTHDIF + ULITPROP + UPROPFEM, data = chist_dat)   

stargazer(m1_chist, type="text", out = "chist.doc", 
          digits = 2,
          align=TRUE, omit.stat=c("LL","ser","f"), no.space=TRUE,
          notes.align = "l", notes.append = F,star.char = c("+", "*", "**", "***"), 
          star.cutoffs = c(.1, .05, .01, .001),
          notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), 
          se = starprep(m1_chist, clusters = chist_dat$CID,se_type="stata"))

#* Mode 5 in Table A-6 with additional control variables ----
dat_add <- read_dta("ifri2.dta") 


dat_add <- as.data.frame(dat_add)
dat_add$UHHNUM <- as.numeric(dat_add$UHHNUM)
dat_add$REFOWNLAND_DUM <- as.factor(dat_add$REFOWNLAND_DUM)

dat_add$SOCAP_OUT_MEAN <- ((dat_add$SOCAP_OUT_MEAN-1)/3)*100
dat_add$SOCAP_IN_MEAN <- ((dat_add$SOCAP_IN_MEAN-1)/3)*100


main1 <- lm(SOCAP_OUT_MEAN ~ REFOWNLAND_DUM + GRPYR + UHHNUM + UWEALTHDIF + ULITPROP + UPROPFEM + SMARKET + GHHCOMM_perc, data = dat_add) 

mods <- list(main1)
for(i in 1:1){
  assign(paste("cl.robust.se.",i,sep=""), 
         sqrt(diag(cluster.vcov(mods[[i]], dat_add$CID)))) # cluster-robust SEs for ols
}

stargazer(main1, 
  se=list(cl.robust.se.1),
  type = "text", star.char = c("+", "*", "**", "***"), star.cutoffs = c(.1, .05, .01, .001),
  notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), notes.append=FALSE)



#* 3. Propensity score matching (Model 4 in Table A-6, Fig A-1) ----
dat$CID <- as.factor(dat$CID)
vars <- c("GRPYR","UHHNUM","UWEALTHDIF","ULITPROP","UPROPFEM","CID","REFOWNLAND_DUM", "SOCAP_OUT_MEAN")
subdat <- dat[which(rowSums(is.na(dat[,vars]))==0),]

subdat <- subset(subdat, CID == "BOL" | CID == "IND" | CID == "GUA" |CID == "MEX" | CID =="UGA")
subdat$GRPYRlog <- log(subdat$GRPYR)

#** Create Trt variable ----
subdat$Trt<- as.numeric(subdat$REFOWNLAND_DUM)-1
table(subdat$CID, subdat$Trt)

#** Define the glm for the propensity score ----
theglm <- glm (Trt ~ GRPYR + UHHNUM + UWEALTHDIF + ULITPROP + UPROPFEM + 
                 strata(CID), data=subdat, family=binomial(link="logit"))
summary(theglm)

subdat$REFOWNLAND_DUM <- as.factor(subdat$REFOWNLAND_DUM)
levels(subdat$REFOWNLAND_DUM) <- c("No land rights","Has land rights")

thepsscor01 <- predict(theglm, type="response")
subdat$pr_score <- thepsscor01
ggplot(aes(x = pr_score), data=subdat) +
  geom_histogram(color = "white") +
  facet_wrap(~REFOWNLAND_DUM) +
  xlab("Propensity score, by actual land right status") +
  theme_bw()
ggsave("psscore.png")


# Match
psmatrix <- match_on(theglm, data=subdat) #, caliper = 0.1, method = "euclidean"

# Do the matches
#  full match
fm1 <- fullmatch(psmatrix,data=subdat, min.controls=0, max.controls=5)
summary(fm1)
subdat$fm1 <- fm1


# Plot the balance of usergroup characteristics
library(gridExtra)
thecovs <- c("GRPYR" ,"UHHNUM","UWEALTHDIF","ULITPROP","UPROPFEM")

plot_list = list()
for(i in 1:length(thecovs)) {
  bppm1 <- ggplot(subdat,aes_string(x="fm1",y=thecovs[i])) +
    geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="red", fill="red")
  
  bporig <- ggplot(subdat,aes_string(y=thecovs[i]))+
    geom_boxplot()
  
  p <- grid.arrange(bppm1,bporig,ncol=2,layout_matrix=matrix(c(1,1,1,1,2),nrow=1))
  plot_list[[i]] <- p
  ggsave(p, file=paste("fig/balance_fm2_", thecovs[i], ".png", sep=""))
  
  bppm2 <- ggplot(subdat,aes_string(x="fm2",y=thecovs[i])) +
    geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="red", fill="red")
  
  bporig <- ggplot(subdat,aes_string(y=thecovs[i]))+
    geom_boxplot()
  
  p <- grid.arrange(bppm2,bporig,ncol=2,layout_matrix=matrix(c(1,1,1,1,2),nrow=1))
  plot_list[[i]] <- p
  ggsave(p, file=paste("fig/balance_fm2_", thecovs[i], ".png", sep=""))
  
  bppm1c <- ggplot(subdat,aes_string(x="fm1.c",y=thecovs[i])) +
    geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="red", fill="red")
  
  bporig <- ggplot(subdat,aes_string(y=thecovs[i]))+
    geom_boxplot()
  
  p <- grid.arrange(bppm1c,bporig,ncol=2,layout_matrix=matrix(c(1,1,1,1,2),nrow=1))
  plot_list[[i]] <- p
  ggsave(p, file=paste("fig/balance_fm1c_", thecovs[i], ".png", sep=""))
}


# Balance check with XBalance (Tables A-7, A-8)
xbres <- xBalance(Trt~GRPYR + UHHNUM + UWEALTHDIF + ULITPROP + UPROPFEM,
                  strata=list(nostrat=NULL,
                              fm1 = ~fm1
                              ),
                  data=subdat,report="all")
xbres$results
xbres$overall

png(filename="balance.png")
plot(xbres)
dev.off()

#evaluate treatment effect
lm1 <- lm(SOCAP_OUT_MEAN ~  Trt + strata(fm1), data=subdat)
stargazer(lm1, type="text", digits = 2,
          star.char = c("+", "*", "**", "***"), star.cutoffs = c(.1, .05, .01, .001),
          notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), notes.append=FALSE)


#* mediation analysis (Fig A-2 and Tablw A-10 ) ----
res <- NA
vars <- c("SOCAP_OUT_MEAN","USANCTIN","REFOWNLAND_DUM","GRPYR","UHHNUM","UWEALTHDIF","ULITPROP","UPROPFEM", "CID")
subdat <- na.omit(dat[,vars])

med1.socap_in <- lm(SOCAP_OUT_MEAN ~ REFOWNLAND_DUM + GRPYR + UHHNUM + UWEALTHDIF + ULITPROP + UPROPFEM, data = subdat)
med2.socap_in <- lm(USANCTIN ~ REFOWNLAND_DUM + GRPYR + UHHNUM + UWEALTHDIF + ULITPROP + UPROPFEM, data = subdat)
med3.socap_in <- lm(SOCAP_OUT_MEAN ~ USANCTIN  + REFOWNLAND_DUM + GRPYR + UHHNUM + UWEALTHDIF + ULITPROP + UPROPFEM, data = subdat)

mods <- list(med1.socap_in, med2.socap_in,med3.socap_in)
for(i in 1:3){
  assign(paste("cl.robust.se.",i,sep=""), 
         sqrt(diag(cluster.vcov(mods[[i]], subdat$CID)))) 
}

stargazer(
  med1.socap_in, med2.socap_in, med3.socap_in, 
  se=list(cl.robust.se.1, cl.robust.se.2, cl.robust.se.3),
  type = "text", star.char = c("+", "*", "**", "***"), star.cutoffs = c(.1, .05, .01, .001),
  notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), notes.append=FALSE)

med.out <- mediate(med2.socap_in, med3.socap_in, treat = "REFOWNLAND_DUM", mediator = "USANCTIN",robustSE = TRUE, sims = 1000)
summary(med.out)

# Sensitivity analysis
sens.out <- medsens(med.out, rho.by = 0.1, effect.type = "indirect", sims = 1000)
png(filename = "plot_sens.png")
plot(sens.out, sens.par = "rho") 
summary(sens.out)
dev.off()

png(filename = "plot_sens2.png")
plot(sens.out, sens.par = "R2", r.type = "total", sign.prod = "positive")
dev.off()

################################# End of the script ########################################

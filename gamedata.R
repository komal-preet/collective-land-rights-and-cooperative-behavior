library(pacman)
p_load(haven, tidyverse, stargazer, estimatr)

dat <- read_dta("datgroupsurv_021221.dta") 

dat <- dat %>% 
  mutate(REFOWNLAND_DUM = ifelse(Site_id ==1, 1, REFOWNLAND_DUM))

m1 <- lm(Harvest_sd ~ REFOWNLAND_DUM +  age_avg + UHHNUM+ UWEALTHDIF+ Treatment+ female_perc+ edu_avg, data = dat, Round <6)
summary(m1)

m2 <- lm(speak_numper ~ REFOWNLAND_DUM +  age_avg + UHHNUM+ UWEALTHDIF+ Treatment+ female_perc+ edu_avg, data = dat, Round <6)

m3 <- lm(speak_perc ~REFOWNLAND_DUM +  age_avg + UHHNUM+ UWEALTHDIF+ Treatment+ female_perc+ edu_avg, data = dat, Round <6)

m4 <- lm(commtrust_avg ~ REFOWNLAND_DUM +  age_avg + UHHNUM+ UWEALTHDIF+ Treatment+ female_perc+ edu_avg, data = dat, Round ==1)

m5 <- lm(gentrust_avg ~ REFOWNLAND_DUM +  age_avg + UHHNUM+ UWEALTHDIF+ Treatment+ female_perc+ edu_avg, data = dat, Round ==1)

stargazer(m1,m2,m3,m4,m5, out = "game.doc", 
          digits = 3,
          align=TRUE, omit.stat=c("f"), no.space=TRUE,
          star.char = c("+", "*", "**", "***"), 
          star.cutoffs = c(.1, .05, .01, .001),
          notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), 
          notes.append=FALSE, 
          se = starprep(m1, m2,m3,m4,m5, se_type="stata"))

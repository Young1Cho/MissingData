####################################################
# Author: Young Won Cho
# Date: 2022-09-18
# Purpose: Figure out what missing pattern model is 
#          suitable for the DAS empirical data.
#          Here, patients's positive affect is a DV.
####################################################


# Read Data ----
emp <- read.csv("C:/Users/YoungwonCho/OneDrive - The Pennsylvania State University/Simulation/Empirical/Cleaned_day8.csv")
names(emp)
emp1 <- emp[,c("dyad", "studyday", "pa.pain.CWC", "pa.supp.CWC",
        "pa.emp.CWC", "pa.neg.CWC", "pa.pos.CWC", "opyrsKOA.GMC",
        "opage.GMC", "osgender.GMC", "marsat.GMC")]
# Insert NA ----
n <- length(unique(emp$dyad))
base <- data.frame(dyad = rep(unique(emp$dyad), each = 22),
                   studyday = rep(1:22, n))
emp2 <- merge(base, emp1, by = c("dyad", "studyday"), all.x = TRUE)
head(emp2, 20)
emp2$pos.na <- ifelse(is.na(emp2$pa.pos.CWC), 1, 0)

# Logistic Model Fitting ----
model <- glm(pos.na ~ studyday + pa.pain.CWC + pa.supp.CWC + opyrsKOA.GMC + opage.GMC + osgender.GMC + marsat.GMC,
data = emp2, family = binomial(link = "logit"))
summary(model)

model1 <- glm(pos.na ~  studyday,
data = emp2, family = binomial(link = "logit")); summary(model1)
model2 <- glm(pos.na ~  studyday + I(studyday^2),
data = emp2, family = binomial(link = "logit")); summary(model2)
plot(emp2$studyday, model2$fitted.values)
plot(emp2$studyday, exp(model$fitted.values) / (1 + exp(model$fitted.values)))

table(emp2$pos.na)

# Result
# It looks like only 'studyday' affect?

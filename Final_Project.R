library(readxl)
library(faraway)
library(MASS)
library(leaps)

dev.new()

broilers <- read_excel("PasturedPoultryFarms.xlsx",sheet = "Broilers")

fecalsoil <- read_excel("PasturedPoultryFarms.xlsx",sheet = "FecalSoil")

Compositional <- read_excel("PasturedPoultryFarms.xlsx",sheet = "Compositional")


##Question 1

model_Camplyobact <- lm("Campylobact~Listeria+Salmonella",data = broilers)
model_Salmonella <- lm("Salmonella~Campylobact+Listeria",data = broilers)
model_Listeria <- lm("Listeria~Salmonella+Campylobact",data = broilers)

summary(model_Camplyobact)
summary(model_Salmonella)
summary(model_Listeria)

model_Camplyobact_binom <- glm("Campylobact~Listeria+Salmonella",data = broilers,
                         family = binomial)
model_Salmonella_binom <- glm("Salmonella~Campylobact+Listeria",data = broilers,
                              family = binomial)
model_Listeria_binom <- glm("Listeria~Salmonella+Campylobact",data = broilers,
                            family = binomial)

summary(model_Camplyobact_binom)
summary(model_Salmonella_binom)
summary(model_Listeria_binom)

boxcox(lm("I(Campylobact+0.1)~Listeria+Salmonella",data = broilers))


cor(broilers[,c("Campylobact","Listeria","Salmonella")],
    use = "complete.obs")
plot(model_Camplyobact$fit,model_Camplyobact$res,pch = 19)




##Question 2
campyl_model <- lm("Campylobact~factor(SampleType)+factor(Month)",
                   data = broilers)
summary(campyl_model)
plot(campyl_model$fit, campyl_model$res,pch = 19,xlab = "Fitted values",
     ylab = "Residual values",main = "Residual vs Fitted Plot: Linear Regression")

boxcox(lm("I(Campylobact+0.1)~factor(SampleType)+factor(Month)",
          data = broilers))

campyl_model_1 <- lm("log(Campylobact+0.1)~factor(SampleType)+factor(Month)",
                     data = broilers)
summary(campyl_model_1)

plot(campyl_model_1$fit, campyl_model_1$res,pch = 19)


campyl_model_binom_logit <- glm("Campylobact~factor(SampleType)+factor(Month)",
                         data = broilers, family = binomial(link = "logit"))
summary(campyl_model_binom_logit)


plot(residuals(campyl_model_binom_logit, type = "deviance"), pch = 19,
     ylab="Deviance",main ="Deviance plot of Logit Model")
abline(h = 0, lty = 2)

#Dispersion factor
campyl_model_binom_logit$deviance/campyl_model_binom_logit$df.residual

##Question 3
fecalsoil$EcoliLog10 <- as.numeric(fecalsoil$EcoliLog10)

eco_model <- lm("I(EcoliLog10 + 0.1) ~ factor(SampleType)+ factor(AnimalSource) + factor(SampleType)*factor(AnimalSource)", 
                data = fecalsoil)
summary(eco_model)


par(mfrow = c(1,2))
plot(eco_model$fit,eco_model$res,pch = 19)
abline(h=0,lty = 2)
identify(eco_model$fit,eco_model$res,atpen = T, tolerance = 0.5)

boxcox(eco_model, lambda = seq(1,2,by = 0.05))
title("Boxcox Plot of linear model")
par(mfrow = c(1,1))

step(eco_model)

Cpplot(leaps(model.matrix(eco_model)[,-1],na.omit(fecalsoil$EcoliLog10)+0.1))

maxadjr(leaps(model.matrix(eco_model)[,-1],na.omit(fecalsoil$EcoliLog10)+0.1,
              method="adjr2"), best = 8)

eco_model_gamma <- glm("I(EcoliLog10+0.1) ~ factor(SampleType)+ factor(AnimalSource) + factor(SampleType)*factor(AnimalSource)", 
                 data = fecalsoil, family = Gamma(link = "inverse"))

summary(eco_model_gamma)


par(mfrow = c(1,2))
plot(na.omit(fecalsoil$EcoliLog10)+0.1,eco_model$fit,pch = 19,
     xlab = "Original values",ylab="Fitted response",main = "Linear Model: Original vs Fitted")
abline(a=0,b=1,lty = 2)

plot(na.omit(fecalsoil$EcoliLog10)+0.1,eco_model_gamma$fit,pch = 19,
     xlab = "Original values",ylab="Fitted response",main = "Gamma Model: Original vs Fitted")
abline(a=0,b=1,lty = 2)
par(mfrow = c(1,1))

par(mfrow = c(1,2))
plot(hat(model.matrix(eco_model)),pch = 19,ylab = 'Leverages',
     main = 'Leverage Plot of Model')
abline(h = 2*ncol(model.matrix(eco_model))/nrow(model.matrix(eco_model)),lty = 2, col = 2)
points(137,hat(model.matrix(eco_model))[137], col = 2,pch = 19)
points(385,hat(model.matrix(eco_model))[385], col = 2,pch = 19)
points(1541,hat(model.matrix(eco_model))[1541], col = 2,pch = 19)


plot(cooks.distance(eco_model),type = 'h',lwd = 3,ylab = 'Cooks distance',
     main = 'Cooks distance Plot for Model')
identify(1:nrow(model.matrix(eco_model)),cooks.distance(eco_model),
         tolerance = 0.5, atpen = TRUE)
par(mfrow = c(1,1))

eco_model_385 <- lm("I(EcoliLog10 + 0.1) ~ factor(SampleType)+ factor(AnimalSource) + factor(SampleType)*factor(AnimalSource)", 
                data = na.omit(fecalsoil)[-385,])

eco_model_137 <- lm("I(EcoliLog10 + 0.1) ~ factor(SampleType)+ factor(AnimalSource) + factor(SampleType)*factor(AnimalSource)", 
                    data = na.omit(fecalsoil)[-137,])

eco_model_1541 <- lm("I(EcoliLog10 + 0.1) ~ factor(SampleType)+ factor(AnimalSource) + factor(SampleType)*factor(AnimalSource)", 
                    data = na.omit(fecalsoil)[-1541,])


par(mfrow = c(1,3))
plot(cooks.distance(eco_model_385),type = 'h',lwd = 3,ylab = 'Cooks distance',
     main = 'Cooks distance Plot for Model w/o 385')
identify(1:nrow(model.matrix(eco_model_385)),cooks.distance(eco_model_385),
         tolerance = 0.5, atpen = TRUE)

plot(cooks.distance(eco_model_137),type = 'h',lwd = 3,ylab = 'Cooks distance',
     main = 'Cooks distance Plot for Model w/o 137')
identify(1:nrow(model.matrix(eco_model_137)),cooks.distance(eco_model_137),
         tolerance = 0.5, atpen = TRUE)

plot(cooks.distance(eco_model_1541),type = 'h',lwd = 3,ylab = 'Cooks distance',
     main = 'Cooks distance Plot for Model w/o 1541')
identify(1:nrow(model.matrix(eco_model_1541)),cooks.distance(eco_model_1541),
         tolerance = 0.5, atpen = TRUE)

par(mfrow = c(1,1))

##Question 4
Question_2_model_w_farm_flock <- glm("Campylobact~factor(SampleType)+factor(Month)",
                                     data = broilers[broilers$Farm == "A",], family = binomial(link = "logit"))
summary(glm("Campylobact~factor(SampleType)+factor(Month)",
            data = broilers[broilers$Farm == "A",], family = binomial(link = "logit")))

library(car)
compareCoefs(campyl_model_binom_logit,Question_2_model_w_farm_flock)
##Question 5

attach(Compositional)
Sample_dummy <- factor(Compositional$Sampletype, labels = c(0,1))

summary(glm(Sample_dummy~A + B + C + D, family = binomial))
summary(glm(Sample_dummy~A + C + D + E, family = binomial))
summary(glm(Sample_dummy~A + B + D + E, family = binomial))
summary(glm(Sample_dummy~A + B + C + E, family = binomial))

compositional_model <- glm(Sample_dummy~A+ B + C + D, family = binomial)
summary(compositional_model_A)

compositional_model_E <- glm(Sample_dummy~ B+ C+D+E, family = binomial)
summary(compositional_model_E)

#Predicted values
Pred_sample <- compositional_model_A$fit > 0.5

#Confusion matrix
table(Pred_sample, Sampletype)

##Cross-Validation
test_num <- as.integer(nrow(Compositional)/3)
index <- sample(nrow(Compositional), test_num)
Sample_dummy_cross_val <- factor(Sampletype[-index],labels = c(0,1))
compos_cross_val_model <- glm("factor(Sampletype,labels = c(0,1)) ~ A+B+C+D",
                              data = Compositional[-index,], family = binomial)
Pred_sample_cross_val <- predict(compos_cross_val_model, Compositional[index,c("A","B","C","D","E")],
                                 type = "response")>0.5
table(Pred_sample_cross_val,Sampletype[index])


##Cross-Validation multiple loops
accuracy_score_list <- c()

for (i in seq(1,10)){
  test_num <- as.integer(nrow(Compositional)/3)
  index <- sample(nrow(Compositional), test_num)
  Sample_dummy_cross_val <- factor(Sampletype[-index],labels = c(0,1))
  compos_cross_val_model <- glm("factor(Sampletype,labels = c(0,1)) ~ A+B+C+D",
                                data = Compositional[-index,], family = binomial)
  Pred_sample_cross_val <- predict(compos_cross_val_model, Compositional[index,c("A","B","C","D","E")],
                                   type = "response")>0.5
  correct_obs <- sum(diag(table(Pred_sample_cross_val,Sampletype[index])))
  accuracy_score_list <- append(accuracy_score_list,(correct_obs/test_num))
  
  
}

summary(1-accuracy_score_list)

1-c(1,2,3)
predict(compos_cross_val_model, Compositional[index,c("A","B","C","D","E")],
        type = "response")
detach(Compositional)

as.integer(101.43)

?sample

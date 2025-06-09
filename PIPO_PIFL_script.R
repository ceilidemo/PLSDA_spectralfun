#################################################
### PLSDA: PIPO_PIFL pre-treatment needle spectra 
### Predominantly sourced Anna Schweiger's script in her PLSDA_example repository
### Ceili DeMarais, 28 October 2024 

library(devtools)
library(spectrolab)
library(lattice)
library(ggplot2)
library(caret)
library(reshape)
library(corrplot)
library (agricolae)

##View spectra
dati<- read.csv("data_in/PretreatmentNeedleSpectra.csv")
head(dati) ##Check .csv for correct input


####PLSDA time
 
### Quick PLSDA - Get general gist and idea of # of components
classi <- as.factor(dati$Species) ## define your classes

wvl <- paste0("X",seq(400,2400, by=10)) ## selects the wavelength range that you want to focus in on 
spec <- dati[,which(colnames(dati)%in%as.character(wvl))] ## creates spec matrix

set.seed(1840) ## random number of observatins in each species
rndid <- with(dati, ave(1:nrow(dati), Species, FUN=function(x) {sample.int(length(x))}))

### Partition data into training and testing data
# table(classi) ### check no of obs per class
inTrain <- rndid<=15 ### e.g. 15 spectra per class for training 
training <- spec[inTrain,] #### resampling optional
testing <- spec[!inTrain,]
trainclass <- classi[inTrain]
testclass <- classi[!inTrain]

#table(testclass) ### checks
#table(trainclass)
#table(dati$Species)

### K(number) folds CV, 10 repeats
ctrl <- trainControl(method = "repeatedcv", repeats = 10, number=10,
                     summaryFunction = multiClassSummary)

# quick model fit, might take a minute
plsFit <- train(training, trainclass, method = "pls", tuneLength = 20,
                trControl = ctrl,probMethod="softmax")

# summary stats to find ncomp
#pdf("./PLSDA_ncomp_summary.pdf",width = 5,height = 4)
plot(plsFit,metric = c("Sensitivity"))
plot(plsFit,metric = c("Specificity"))
plot(plsFit,metric = c("Kappa"))
plot(plsFit,metric = c("Accuracy"))
#dev.off()

plsProbs <- predict(plsFit, newdata = testing, type = "prob")
plsClasses <- predict(plsFit, newdata = testing)

# confusion matrix
(confu <- confusionMatrix(data = plsClasses, testclass))
confutab <- as.data.frame.matrix(confu$table)


# probabilities
probis <- as.data.frame(plsProbs)
probis <- cbind(testclass, probis)

## ok nice, good quick PLSDA, but now we want to iterate (I chose 50 iterations)

#########################
#### ITERATE PLSDA ###### 
#########################
classi <- as.factor(dati$Species) ### define classes

wvl <- paste0("X",seq(400,2400, by=10)) ### again, definind the wvl range and saying select every 10th wvl
spec <- dati[,which(colnames(dati)%in%as.character(wvl))] ### make spec matrix

### PLSDA iterative models for finding number of components
### this takes a while, decrease iterations (nsims) for test run

nsims <- 50 ### number of iterations
rndid <- list() ### specify no samples per species for training
set.seed(1840)
for (i in 1:nsims){ 
  rndid[[i]] <- with(dati, ave(1:nrow(dati), Species, FUN=function(x) {sample.int(length(x))}))
}

mods <- list()
for (nsim in seq(nsims)){
  inTrain <- rndid[[nsim]]<= 15 ### training on 15 samples per species
  print(nsim)
  flush.console()
  set.seed(nsim)
  traini <- spec[inTrain,seq(1,length(wvl), by=10)] 
  testi <- spec[!(inTrain),seq(1,length(wvl), by=10)]
  trainclass <- classi[inTrain];testclass <- classi[!(inTrain)]
  plsFit <- train(traini, trainclass, method = "simpls", tuneLength = 20,
                  trControl = trainControl(method="boot")) ### in this line we may have probMethod="Bayes" for larger datasets
  mods[[nsim]] <- plsFit
}

head(mods[[2]]$results) ### check that it worked

### Number of components 
ncomps <- vector(length = nsims)
for (i in 1:nsims){
  ncomps[i]<-mods[[i]]$finalModel$ncomp
}
table(ncomps)

### Kappa statistics
kappas <- data.frame(ncomps= 1:20,matrix(NA, nrow = 20, ncol = length(mods)))
for (i in 1:length(mods)){
  kappas[,i+1] <- mods[[i]]$results$Kappa
}

### Tukey test
kapp <- as.data.frame(as.numeric(t(kappas[,-1])))
kapp <- cbind(kapp, rep(1:20, each=length(mods)))
names(kapp) <- c("Kappa", "ncomps")

kapp$ncomps <- as.factor(kapp$ncomps)

modi <- lm (Kappa~ncomps, kapp) 
tuk <- HSD.test(modi,"ncomps")

tuk_dat <- as.data.frame(tuk$groups)
tuk_dat$var <- as.numeric(row.names(tuk_dat))
tuk_dat <- tuk_dat[order(tuk_dat$var,decreasing = F),]
letters <- as.character(tuk_dat$groups)

#### Ok, now let's visualize
## Kappa Tukey plot, we want high Kappa

#pdf("./PLSDA_kappas.pdf",width = 5,height = 4)
par(bty="l")
boxplot(kapp$Kappa~kapp$ncomps,ylim=c(0,1),
        xlab="Number of components",ylab="Kappa")
text(x=1:20, y=rep(1,20),letters) ##looks like it may level off at highest Kappa around 5/6 components
# dev.off()

#####################
#### Final model ###
compi <- 5 ### select number of components based on Kappa plot / Tukey test from above - I chose 5 components
finmods <- list()
nsims=10

for (nsim in 1:nsims){
  print(nsim)
  flush.console()
  set.seed(nsim)
  inTrain <- rndid[[nsim]]<= 15
  # inTrain <- createDataPartition(y =classi, p = .70, list = FALSE)
  training <- spec[inTrain,seq(1,length(wvl), by=10)]
  testing <- spec[!(inTrain),seq(1,length(wvl), by=10)]
  trainclass <- as.factor(classi[inTrain]); testclass <- as.factor(classi[!(inTrain)])
  
  finalModel <- plsda(training,trainclass, ncomp=compi, method = "simpls")
  finmods[[nsim]] <- finalModel
}

# saveRDS(finmods, "./PLSDA_finmods.rds")

col_ord <- c('#313695','#4575b4','#74add1','#abd9e9',"darkorchid4",'yellow2','#a50026','#d73027','#f46d43',
             '#fdae61',"grey50","black")

out <- list()
for (i in 1:nsims) {
  out[[i]] <- as.matrix(finmods[[i]]$scores[1:length(trainclass), 1:ncol(finmods[[i]]$scores)])
  out[[i]] <- as.data.frame(out[[i]])
  colnames(out[[i]]) <- gsub(pattern = " ", "_",colnames(out[[i]]))
  out[[i]]$Species <- as.character(trainclass)
  out[[i]] <- out[[i]][order(out[[i]]$Species,decreasing = T),]
  col_ordi <- rep(col_ord[1:length(unique(dati$Species))],as.numeric(table(out[[i]]$Species)))
  out[[i]]$col_ordi <- rev(col_ordi)
}
# write.csv(out, "./PLSDA_scores.csv", row.names=F)

### Probabilities and confusion matrix
probis <- list()
confus <- list()

for (nsim in seq(nsims)){
  print(nsim)
  flush.console()
  set.seed(nsim)
  inTrain <- rndid[[nsim]]<= 15
  # inTrain <- createDataPartition(y = classi, p = .80, list = FALSE)
  testing <- spec[!(inTrain),seq(1,length(wvl), by=10)]
  testclass <- as.factor(classi[!(inTrain)])
  
  plsProbs <- predict(finmods[[nsim]], newdata = testing, type = "prob")
  plsClasses <- predict(finmods[[nsim]], newdata = testing)
  confus[[nsim]] <- confusionMatrix(data = plsClasses, testclass)
  
  probs <- as.data.frame(plsProbs)
  names(probs) <- sapply(strsplit(names(probs),split = ".n"),"[",1)
  probs <- cbind(testclass, probs)
  probis[[nsim]] <- probs 
}


### Probability plot time!! 
arr <- array(unlist(probis), dim = c(dim(probis[[1]]),nsims))
prob_mean <- apply(arr, 1:2, mean)
prob_mean <- as.data.frame(prob_mean)
prob_mean$V1 <- probis[[1]]$testclass
colnames(prob_mean) <- colnames(probis[[1]])
# write.csv(prob_mean,"./PLSDA_probmean.csv")

# Calculate accuray from probmean
pp <- melt(prob_mean, id="testclass")
pp$position <- ifelse(as.character(pp$testclass) == as.character(pp$variable), 2, 1)
pp$testclass <- factor(pp$testclass, levels=rev(levels(pp$testclass)))

coli <- c("darkorchid4",'#313695','#4575b4','#74add1')
# pdf("./plots/probability_plot.pdf",width = 6,height = 4)
ggplot(pp, aes(x=testclass, y=value, fill=variable, group=position))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values= alpha(coli,1))+
  theme(legend.title=element_blank())+
  labs(x = "Species", y= "")+
  coord_flip()
# dev.off()

### Confusion table plot - this is usually the mose useful for visualization
tabs <- list()
for(i in 1:length(confus)){
  tabs[[i]] <- confus[[i]]$table
}

tabsi <- Reduce('+', tabs)
tab_mean <- as.data.frame.matrix(tabsi/length(confus))
# write.csv(tab_mean,"./PLSDA_confumean.csv")

sums <- colSums(tab_mean)
tabs_perc <- matrix(NA, length(sums),length(sums))
for (i in 1:length(sums)){
  tabs_perc[,i] <- tab_mean[,i]/sums[i]
}

colnames(tabs_perc) <- colnames(confus[[1]]$table)
rownames(tabs_perc) <- rownames(confus[[1]]$table)
# write.csv(tabs_perc,"./PLSDA_confuperc.csv")

col <- colorRampPalette(c("black","black","brown","gold","forestgreen")) 

# pdf("./plots/PLSDA_corrplot.pdf",width = 7,height = 6,pointsize = 11)
corrplot(tabs_perc, p.mat = tabs_perc, insig = "p-value", sig.level = -1, addCoef.col = 1,
         tl.srt = 70,col = col(20),cl.lim = c(0, 1),tl.col = 1, tl.offset =1.5, 
         cl.ratio = 0.2, cl.align.text = "l", cl.cex = 0.9, 
         mar=c(1,3,3,3))
mtext("Prediction",2, line=0, cex=1.2)
mtext("Reference",at = 2, line = 1.5, cex=1.2)
# dev.off()

#### Importance of bands! The loadings plot helps visualise which bands are most important for differentiation 
lls <- list()
for(i in 1:length(finmods)){
  lls[[i]] <- abs(loadings(finmods[[i]])[1:dim(loadings(finmods[[1]]))[1],1:compi])
  sumis <- lapply(lls,rowSums)
}

mm <- apply(simplify2array(sumis), 1, mean)
ss <- apply(simplify2array(sumis), 1, sd)

mm <- as.data.frame(mm)
mm <- cbind(mm,ss)

# pdf("./plots/PLSDA_loadings.pdf",width = 6,height = 4)
plot(1:ncol(testing), mm$mm, type="n", bty="l", ylab="abs (loadings)", 
     xaxt="n", xlab="Wavelength (nm)", 
     ylim=c(min(mm$mm)-mm[which(mm$mm==min(mm$mm)),]$ss, 
            max(mm$mm)+mm[which(mm$mm==max(mm$mm)),]$ss))
polygon(x=c(1:ncol(testing),ncol(testing):1), y=c(mm$mm-(mm$ss), rev(mm$mm+(mm$ss))),col = "grey", border = "grey")
lines(1:ncol(testing), mm$mm, type="l")
axis(1, at=seq(1,ncol(testing),4), labels=seq(400,2400,400))
# dev.off()


### Modelsstatistics ##

## Note: Accumean and accuracies in mods are based on probablities
## In the confuplot we are showing accuracies based on the 0/1 decisions
## both are relevant for reporting
## we make our class decision based on 0/1 result
## but when accumean is low model generalizability is not good (e.g. a mean accuracy of 60% is not very trustworthy)

### Accuracy and Kappa statistic
accu <- numeric(length=length(confus))
for (i in 1:length(confus)){
  accu[i] <- confus[[i]]$overall[1]
}

(accmean <- mean(accu)) #mean accuracy
(accsd <- sd(accu))  #SD of accuracy

kappas <- numeric(length=length(confus))
for (i in 1:length(confus)){
  kappas[i] <- confus[[i]]$overall[2]
}

(kappmean <- mean(kappas)) ## Mean Kappa
(kappsd <- sd(kappas))     ## SD of Kappa


### Sensitivity, specificity, precision
specifi <- list()
sensi <- list()
preci <- list()
mean_specifi <- list()
mean_sensi <- list()


for (i in 1:nsims){
  (confu <- confus[[i]]$table)
  (n <- sum(confu)) # number of instances
  (nc <- nrow(confu)) # number of classes
  (diag <-  diag(as.matrix(confu))) # number of correctly classified instances per class
  (rowsums <- apply(confu, 1, sum)) # number of instances per class
  (colsums <-  apply(confu, 2, sum)) # number of predictions per class
  (p <- rowsums / n) # distribution of instances over the actual classes
  (q <- colsums / n) # distribution of instances over the predicted classes
  preci[[i]] <- diag / colsums ### TP/TP+FP
   accuracy <- sum(diag) / n ### accurcay
   accu[[i]] <- accuracy
   expAccuracy <- sum(p*q) ### kappa
   kappas[[i]] <- (accuracy - expAccuracy) / (1 - expAccuracy)
  
  ### sensitivity, true positive rate: TP/TP+FN, TP/P
  sensi[[i]] <- diag / rowsums
  ### specificity, true negative rate: TN/TN+FP TN/N
  speci <- numeric(length = nrow(confu))
  for (j in 1:nrow(confu)){
    speci[j] <- sum(confu[-j,-j])/sum(confu[,-j])
  }
  names(speci) <- colnames(confu)
  specifi[[i]] <- speci
  ##### Specificity and Sensitivity
  mean_specifi[[i]] <- sum(speci)/nc
  mean_sensi[[i]] <- sum(sensi[[i]])/nc
}

mean(unlist(sensi)) ### mean sensitivity
sd(unlist(sensi)) ### standard deviation sensitivity

mean(unlist(specifi)) ### same for specificity
sd(unlist(specifi))

mean(unlist(preci)) ### same for precision
sd(unlist(preci))

mean(unlist(accu)) ### Get that accuracy mean
sd(unlist(accu))

mean(unlist(kappas)) ### Finally, kappa mean
sd(unlist(kappas))

### Just putting these all into a table 
modelstatistics_table <- data.frame(
  Metric = c("Sensitivity", "Specificity", "Precision", "Accuracy", "Kappa"),
  Mean = c(mean(unlist(sensi)), mean(unlist(specifi)), mean(unlist(preci)), 
           mean(unlist(accu)), mean(unlist(kappas))),
  SD = c(sd(unlist(sensi)), sd(unlist(specifi)), sd(unlist(preci)), 
         sd(unlist(accu)), sd(unlist(kappas)))
)

# write.csv(modelstatistics_table,"./plots/modelstatistics.csv")

#### END ######

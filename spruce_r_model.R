#read in dataset

library(lme4)
library(lmerTest)
library(bbmle)
data<-read.csv(file.choose(),check.names=FALSE,header=TRUE)
head(data)
fp <-subset(data,Type=="FP")
ig <-subset(data,Type=="IP")

fp$Depth <-as.factor(freshpeat$Depth)
ig$Depth <-as.factor(ingrowth$Depth)

#check for normality
hist(fp$Moisture)
boxplot(fp$Moisture)
qqnorm(fp$Moisture)
qqline(fp$Moisture)
qqnorm(log10(fp$Moisture))
qqline(log10(fp$Moisture))
hist(log10(fp$Moisture))
boxplot(log10(fp$Moisture))
fpmoistnorm <-subset(fp,Moisture>85)
boxplot(fpmoistnorm$Moisture)
qqnorm(log10(fpmoistnorm$Moisture))
qqline(log10(fpmoistnorm$Moisture))
qqnorm(fpmoistnorm$Moisture)
qqline(fpmoistnorm$Moisture)


boxplot(fp$DOC)
qqnorm(fp$DOC)
qqline(fp$DOC)
qqnorm(log10(fp$DOC))
qqline(log10(fp$DOC))
hist(log10(fp$DOC))
boxplot(log10(fp$DOC))

fpDOCn <-subset(fp,DOC<7)
boxplot(fpDOCn$DOC)
qqnorm(fpDOCn$DOC)
qqline(fpDOCn$DOC)
qqnorm(log10(fpDOCn$DOC))
qqline(log10(fpDOCn$DOC))
hist(log10(fpDOCn$DOC))
boxplot(log10(fpDOCn$DOC))


boxplot(fp$TN)
qqnorm(fp$TN)
qqline(fp$TN)
qqnorm(log10(fp$TN))
qqline(log10(fp$TN))
hist(log10(fp$TN))
boxplot(log10(fp$TN))


boxplot(fp$MBC)
qqnorm(fp$MBC)
qqline(fp$MBC)
qqnorm(log10(fp$MBC))
qqline(log10(fp$MBC))
hist(log10(fp$MBC))
boxplot(log10(fp$MBC))


boxplot(fp$MBN)
qqnorm(fp$MBN)
qqline(fp$MBN)
qqnorm(log10(fp$MBN))
qqline(log10(fp$MBN))
hist(log10(fp$MBN))
boxplot(log10(fp$MBN))

boxplot(fp$MBCMBN)
qqnorm(fp$MBCMBN)
qqline(fp$MBCMBN)
qqnorm(log10(fp$MBCMBN))
qqline(log10(fp$MBCMBN))
hist(log10(fp$MBCMBN))
boxplot(log10(fp$MBCMBN))

boxplot(fp$BG)
qqnorm(fp$BG)
qqline(fp$BG)
qqnorm(log10(fp$BG))
qqline(log10(fp$BG))
hist(log10(fp$BG))
boxplot(log10(fp$BG))

boxplot(fp$BX)
qqnorm(fp$BX)
qqline(fp$BX)
qqnorm(log10(fp$BX))
qqline(log10(fp$BX))
hist(log10(fp$BX))
boxplot(log10(fp$BX))

fpBXn <-subset(fp,BX<10000)
boxplot(fpBXn$BX)
qqnorm(fpBXn$BX)
qqline(fpBXn$BX)
qqnorm(log10(fpBXn$BX))
qqline(log10(fpBXn$BX))
hist(log10(fpBXn$BX))
boxplot(log10(fpBXn$BX))


boxplot(fp$CB)
qqnorm(fp$CB)
qqline(fp$CB)
qqnorm(log10(fp$CB))
qqline(log10(fp$CB))
hist(log10(fp$CB))
boxplot(log10(fp$CB))

fpCBn <-subset(fp,CB<10000)
boxplot(fpCBn$CB)
qqnorm(fpCBn$CB)
qqline(fpCBn$CB)
qqnorm(log10(fpCBn$CB))
qqline(log10(fpCBn$CB))
hist(log10(fpCBn$CB))
boxplot(log10(fpCBn$CB))


boxplot(fp$AP)
qqnorm(fp$AP)
qqline(fp$AP)
qqnorm(log10(fp$AP))
qqline(log10(fp$AP))
hist(log10(fp$AP))
boxplot(log10(fp$AP))


boxplot(fp$NAG)
hist(fp$NAG)
qqnorm(fp$NAG)
qqline(fp$NAG)
qqnorm(log10(fp$NAG))
qqline(log10(fp$NAG))
hist(log10(fp$NAG))
boxplot(log10(fp$NAG))


boxplot(fp$Leu)
qqnorm(fp$Leu)
qqline(fp$Leu)
qqnorm(log10(fp$Leu))
qqline(log10(fp$Leu))
hist(log10(fp$Leu))
boxplot(log10(fp$Leu))


fpLeun <-subset(fp,Leu<3000)
boxplot(fpLeun$Leu)
qqnorm(fpLeun$Leu)
qqline(fpLeun$Leu)
qqnorm(log10(fpLeun$Leu))
qqline(log10(fpLeun$Leu))
hist(log10(fpLeun$Leu))
boxplot(log10(fpLeun$Leu))

boxplot(fp$Ala)
qqnorm(fp$Ala)
qqline(fp$Ala)
qqnorm(log10(fp$Ala))
qqline(log10(fp$Ala))
hist(log10(fp$Ala))
boxplot(log10(fp$Ala))


boxplot(fp$AAP)
qqnorm(fp$AAP)
qqline(fp$AAP)
qqnorm(log10(fp$AAP))
qqline(log10(fp$AAP))
hist(log10(fp$AAP))
boxplot(log10(fp$AAP))

##############################################
#These are the normalized data:  log10(fpDOCn), log10(fp$TN), log10(fp$MBC), fp$MBN, log10(fp$BG), 
#  log10(fpCBn$CB), log10(fp$AP), log10(fp$NAG), log10(fp$Ala), log10(fp$AAP)

#questionable are Moisture, MBCMBN, BX, Leu

###############################################################################

##    Determine if hum hollow matters for equivilent depths
#exclue != +10 depth
fpbelow<-subset(fp,Depth != "10")

dim(fp)
dim(fpbelow)

model <-lmer(log10(fpbelow$NAG)~Topography*Depth*Date + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(NAG)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=TRUE)
model3 <-lmer(log10(fpbelow$NAG)~as.factor(Depth) +(1|Location),data=fpbelow,REML=FALSE)
model4 <-lmer(log10(fpbelow$NAG)~as.factor(Depth)*Date +(1|Location),data=fpbelow,REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)

#######################################  ??


model2 <-lmer(log10(NAG)~as.factor(Depth) +(1|Location),data=fp,REML=FALSE)

anova(model2)
 
model <-lmer(log10(NAG)~Depth + (1|Location/Topography/Depth),data=fp,REML=FALSE)
anova(model)

model <-lmer(log10(NAG)~Topography*Depth + (1|Location/Topography/Depth),data=fp,REML=TRUE)
anova(model)

model <-lmer(log10(NAG)~Topography*Depth + (1|Location/Topography/Depth),data=fp,REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(NAG)~Topography*Depth + (1|Location/Topography/Depth),data=fp,REML=TRUE)
model3 <-lmer(log10(NAG)~as.factor(Depth) +(1|Location),data=fp,REML=FALSE)
model4 <-lmer(log10(NAG)~as.factor(Depth)*Date +(1|Location),data=fp,REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)

###########################################           ************                ##################
model <-lmer(log10(NAG)~Depth*Date*Topography + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
anova(model)
difflsmeans(model, ddf="Satterthwaite",type=3,method.grad="simple")

model <-lmer(log10(NAG)~Depth*Date*Topography + (1|Location/Topography/Depth),data=fp,REML=FALSE)
anova(model)
model <-lmer(log10(Ala)~Depth*Date*Topography + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
anova(model)

library(gplots)
library(ggplot2)
library(plyr)
library(reshape)


NAG.n<-tapply (fp$NAG, fp$Depth, length)
NAG.u<-tapply(fp$NAG, fp$Depth, mean)
NAG.SD<-tapply(fp$NAG, fp$Depth, sd)
NAG.SE<-NAG.SD/sqrt(NAG.n-1)

Depth.colors<-c("black", "gray", "white")
x.axis<-c(" ", " ", " ")
pNAGlot<-barplot2(NAG.u, col=Depth.colors, names.arg=x.axis, plot.ci=TRUE,ci.l=NAG.u-NAG.SE, ci.u=NAG.u+NAG.SE, width=1, ylim=c(0,15000), ci.lwd=1.5, cex.axis=1, font=2) 
box(lty="solid", lwd=2, col="black")

#NAG Point plot, Depth effect
Depth.names<-c("+10","-10","-20NAG")
NAG.Depth<-data.frame(Depth.names, NAG.u, NAG.SE)
NAG.Depth$Depth.names<-factor(Depth.names, levels=unique(as.character(Depth.names)))
colors.Depth<-c("goldenrod3","blueviolet","darkgreen")
ylab<-expression(paste("NAGNAG.Depth (mg C Kg"^-1," dry soil)"))

ggplot(NAG.Depth, aes(Depth.names, NAG.u))+geom_pointrange(aes(colour=crop.names, ymin=MBCg.u-MBCg.SE, ymax=MBCg.u+MBCg.SE, size="1"))+ylab(ylab)+coord_cartesian(ylim=c(0,650))+scale_color_manual(values=colors.crop)+
theme(axis.line=element_line(size=2), axis.ticks=element_line(size=2, colour="black"), legend.position="none", panel.background=element_blank(), panel.grid=element_blank(), axis.text=element_text(size=20, face="bold", colour="black"), axis.title.x=element_blank(), axis.title.y=element_text(size=24, face="bold"))+
annotate("text", label="B", x=1, y=300, cex=7)+annotate("text", label="A", x=2, y=560, cex=7)+annotate("text", label="A", x=3, y=620, cex=7)+annotate("text", label="P<0.0001", x=1, y=600, font=2, cex=7)




#MBC_g Point plot, crop effect
crop.names<-c("Corn","Prairie","Fertilized Prairie")
MBCg.crop<-data.frame(crop.names, MBCg.u, MBCg.SE)
MBCg.crop$crop.names<-factor(crop.names, levels=unique(as.character(crop.names)))
colors.crop<-c("goldenrod3","blueviolet","darkgreen")
ylab<-expression(paste("Microbial Biomass C (mg C Kg"^-1," dry soil)"))

ggplot(MBCg.crop, aes(crop.names, MBCg.u))+geom_pointrange(aes(colour=crop.names, ymin=MBCg.u-MBCg.SE, ymax=MBCg.u+MBCg.SE, size="1"))+ylab(ylab)+coord_cartesian(ylim=c(0,650))+scale_color_manual(values=colors.crop)+
theme(axis.line=element_line(size=2), axis.ticks=element_line(size=2, colour="black"), legend.position="none", panel.background=element_blank(), panel.grid=element_blank(), axis.text=element_text(size=20, face="bold", colour="black"), axis.title.x=element_blank(), axis.title.y=element_text(size=24, face="bold"))+
annotate("text", label="B", x=1, y=300, cex=7)+annotate("text", label="A", x=2, y=560, cex=7)+annotate("text", label="A", x=3, y=620, cex=7)+annotate("text", label="P<0.0001", x=1, y=600, font=2, cex=7)






model <-lmer(log10(NAG/MBC)~Depth*Date*Topography + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
anova(model)

boxplot(fp$NAG/fp$MBC)
qqnorm(fp$NAG/fp$MBC)
qqline(fp$NAG/fp$MBC)
qqnorm(log10(fp$NAG/fp$MBC))
qqline(log10(fp$NAG/fp$MBC))
hist(log10(fp$NAG/fp$MBC))
boxplot(log10(fp$NAG/fp$MBC))

model <-lmer(log10(NAG/MBC)~Depth*Date*Topography + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
anova(model)
model <-lmer(log10(BG/MBC)~Depth*Date*Topography + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
anova(model)
model <-lmer(log10(AP/MBC)~Depth*Date*Topography + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
anova(model)
model <-lmer(log10(Ala/MBC)~Depth*Date*Topography + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
anova(model)
model <-lmer(log10(AAP/MBC)~Depth*Date*Topography + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
anova(model)


#These are the normalized data:  log10(fpDOCn), log10(fp$TN), log10(fp$MBC), fp$MBN, log10(fp$BG), 
#  log10(fpCBn$CB), log10(fp$AP), log10(fp$NAG), log10(fp$Ala), log10(fp$AAP)

#############################################
model <-lmer(log10(BG)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(NAG)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=TRUE)
model3 <-lmer(log10(fpbelow$NAG)~as.factor(Depth) +(1|Location),data=fpbelow,REML=FALSE)
model4 <-lmer(log10(fpbelow$NAG)~as.factor(Depth)*Date +(1|Location),data=fpbelow,REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)

model <-lmer(log10(Ala)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(Ala)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=TRUE)
model3 <-lmer(log10(Ala)~as.factor(Depth) +(1|Location),data=fpbelow,REML=FALSE)
model4 <-lmer(log10(Ala)~as.factor(Depth)*Date +(1|Location),data=fpbelow,REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)

model <-lmer(log10(AAP)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(AAP)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=TRUE)
model3 <-lmer(log10(AAP)~as.factor(Depth) +(1|Location),data=fpbelow,REML=FALSE)
model4 <-lmer(log10(AAP)~as.factor(Depth)*Date +(1|Location),data=fpbelow,REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)

model <-lmer(log10(AP)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(AP)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=TRUE)
model3 <-lmer(log10(AP)~as.factor(Depth) +(1|Location),data=fpbelow,REML=FALSE)
model4 <-lmer(log10(AP)~as.factor(Depth)*Date +(1|Location),data=fpbelow,REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)

model <-lmer(log10(BG)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(BG)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=TRUE)
model3 <-lmer(log10(BG)~as.factor(Depth) +(1|Location),data=fpbelow,REML=FALSE)
model4 <-lmer(log10(BG)~as.factor(Depth)*Date +(1|Location),data=fpbelow,REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)

model <-lmer(log10(fpCBn$CB)~Topography*Depth + (1|Location/Topography/Depth),data=subset(fp,CB<10000 & Depth !="10"),REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(fpCBn$CB)~Topography*Depth + (1|Location/Topography/Depth),data=subset(fp,CB<10000 & Depth !="10"),REML=TRUE)
model3 <-lmer(log10(fpCBn$CB)~as.factor(Depth) +(1|Location),data=subset(fp,CB<10000 & Depth !="10"),REML=FALSE)
model4 <-lmer(log10(fpCBn$CB)~as.factor(Depth)*Date +(1|Location),data=subset(fp,CB<10000 & Depth !="10"),REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)

model <-lmer(log10(MBC)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(MBC)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=TRUE)
model3 <-lmer(log10(MBC)~as.factor(Depth) +(1|Location),data=fpbelow,REML=FALSE)
model4 <-lmer(log10(MBC)~as.factor(Depth)*Date +(1|Location),data=fpbelow,REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)

model <-lmer(MBN~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
#try also REML=TRUE
model2 <-lmer(MBN~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=TRUE)
model3 <-lmer(MBN~as.factor(Depth) +(1|Location),data=fpbelow,REML=FALSE)
model4 <-lmer(MBN~as.factor(Depth)*Date +(1|Location),data=fpbelow,REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)


model <-lmer(log10(TN)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(TN)~Topography*Depth + (1|Location/Topography/Depth),data=fpbelow,REML=TRUE)
model3 <-lmer(log10(TN)~as.factor(Depth) +(1|Location),data=fpbelow,REML=FALSE)
model4 <-lmer(log10(TN)~as.factor(Depth)*Date +(1|Location),data=fpbelow,REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)


model <-lmer(log10(DOC)~Topography*Depth + (1|Location/Topography/Depth),data=subset(fp,DOC<7 & Depth !="10"),REML=FALSE)
#try also REML=TRUE
model2 <-lmer(log10(DOC)~Topography*Depth + (1|Location/Topography/Depth),data=subset(fp,CB<10000 & Depth !="10"),REML=TRUE)
model3 <-lmer(log10(DOC)~as.factor(Depth) +(1|Location),data=subset(fp,CB<10000 & Depth !="10"),REML=FALSE)
model4 <-lmer(log10(DOC)~as.factor(Depth)*Date +(1|Location),data=subset(fp,CB<10000 & Depth !="10"),REML=FALSE)
anova(model)
anova(model2)
anova(model3)
anova(model4)

#These are the normalized data:  log10(fpDOCn), log10(fp$TN), log10(fp$MBC), fp$MBN, log10(fp$BG), 
#  log10(fpCBn$CB), log10(fp$AP), log10(fp$NAG), log10(fp$Ala), log10(fp$AAP)
###########################################################################################################
############################################################################################################
#exclue != +10 depth
freshpeatbelow<-subset(freshpeat,Depth != "10")
freshpeatbelow
data
freshpeat
str(freshpeat)
dim(freshpeat)
dim(data)
dim(freshpeatbelow)

hist(freshpeatbelow$BG_Activity)
hist(freshpeat$BG_Activity)

boxplot(freshpeat$NAG_Activity)
freshpeat.NAG.nooutliers <-subset(freshpeat,NAG_Activity<25000 & Depth !="10")
hist(freshpeat.NAG.nooutliers$NAG_Activity)
hist(log10(freshpeat.NAG.nooutliers$NAG_Activity))

qqnorm(freshpeat.NAG.nooutliers$NAG_Activity)
qqline(freshpeat.NAG.nooutliers$NAG_Activity)

qqnorm(log10(freshpeat.NAG.nooutliers$NAG_Activity))
qqline(log10(freshpeat.NAG.nooutliers$NAG_Activity))

freshpeat.NAG.nooutliers.below<-subset(freshpeat.NAG.nooutliers,Depth != "10")
dim(freshpeat.NAG.nooutliers.below)

model <-lmer(log10(NAG_Activity)~Topography*Depth + (1|Location/Topography/Depth),data=freshpeat.NAG.nooutliers.below,REML=FALSE)
#try also REML=TRUE

model
anova(model)

model2 <-lmer(log10(NAG_Activity)~as.factor(Depth) +(1|Location),data=freshpeat.NAG.nooutliers.below,REML=FALSE)
model2 <-lmer(log10(NAG_Activity)~as.factor(Depth)*Date +(1|Location),data=freshpeat.NAG.nooutliers.below,REML=FALSE)

anova(model2)


model2 <-lmer(log10(NAG_Activity)~as.factor(Depth) +(1|Location),data=freshpeat.NAG.nooutliers.below,REML=FALSE)




####################################################   BG   ####################################
#Because hum ho w/i each rep, we'll use a split plot design, and depth is attached to each topography, so the random term is (1|Location/Topography/Depth) ;  coding random term as (1|term)  /is for nesting

model <-lmer(log10(BG_Activity)~Topography*Depth + (1|Location/Topography/Depth),data=freshpeatbelow,REML=FALSE)
#try also REML=TRUE

model
anova(model)
Run BG_Activity/MB_activity to test for normality with subsetted data set.
model <-lmer(log10(BG_Activity/MB_activity)~Depth + (1|Location/Depth),data=subset(freshpeat,NAG_Activity<25000 & Depth !="10" & MB < 7),REML=FALSE)
#try also REML=TRUE

model
anova(model)
# topography doesn't affect bg belowground
freshpeat.BG.nooutliers <-subset(freshpeat,BG_Activity<30000)

model2 <-lmer(log10(BG_Activity)~as.factor(Depth) +(1|Location),data=freshpeat.BG.nooutliers,REML=FALSE)

anova(model2)

#MODEL BELOW WORKS; MODEL ABOVE MUST HAVE +10 REMOVED FROM HUMMOCK IN ORDER TO WORK
#MODEL BELOW JUST ASKS IF DEPTH DIFFERS IN FRESH PEAT AND TOPOGRAPHY
#Or just look at hummock and hollow and depth with that
model2 <-lmer(log10(BG_Activity)~as.factor(Depth) +(1|Location),data=subset(freshpeat.BG.nooutliers,Topography=="Hi"),REML=FALSE)

anova(model2)


#To test different models
#Many more models can be tested, here I just included a model with main effects

#compare with AIC and anova for fitted models 
AICtab(model,model2) #,model3)
#AICtab in bbmle
anova(model,model2)#,model3)
#clearly the main effects model was the best for the data as it had the lowest AIC score by at least 16 (you have lots of confidence > 10 and less when below)

#pvalues are generated for pairwise comparisons based on the specified model.  These are not corrected p-values. Satterthwaite approx . is used as you would do in proc mixed in SAS
difflsmeans(test.model.main, ddf="Satterthwaite",type=3,method.grad="simple")

#check this against plain old anova, and they have similar results.  The effect of month is more clear with the mixed model, however, the cropping effect is still not therer




#sélection du Dossier de travail
setwd("D:/Documents/Cours_Master_1_Bioinformatique/Statistique sur les grands échantillons/Bootstrap") 

#importation des données de "basegenes"
basegene<-read.table("basegenes.txt", header= FALSE, sep=";", dec=".", na.strings="NA")


#importation des données de "baseclinique"
baseclinique<-read.table("baseclinique.csv", header= TRUE, sep=";", dec=".", na.strings="NA")


#extraction de baseclinique des individu servant a l'apprentissage

apprentissage<-baseclinique[baseclinique$Analysis=="Training",]
#extraction de baseclinique des individu servant à la validatione
validation<-baseclinique[baseclinique$Analysis=="Validation",]

###
# MTP
##

library(Biobase)
library(multtest)
library(genefilter)

dim(basegene)

X<-basegene[basegene$V1 %in% apprentissage$LYM.number,]

X<-X[-1]

X<-t(X)

mb<-apprentissage$Status
table(mb)

dim(X)
mb <- mb[!is.na(mb)]
length(mb)

mb.boot<-MTP(X = X , Y=mb, test = "t.twosamp.equalvar", robust=TRUE, alpha=c(0.05), B = 1000, get.cutoff = TRUE)

X <- X[mb.boot@adjp<=0.05,]
summary(logit.model <- glm(mb ~ t(X) , binomial))$coef

S <- logit.model$linear.predictor
exp(summary(logit.model)$coef[2,1])

se<-function(x) { return(sum(1*(S[mb==1]>x))/sum(mb==1)) }
sp<-function(x) { return(sum(1*(S[mb==0]<=x))/sum(mb==0)) }
temp<-sort(unique(S))
resultats<-data.frame( sp=sapply(temp, sp),sp1=1-sapply(temp, sp), se=sapply(temp, se), seuils=temp)

plot(c(1, resultats$sp1, 0), c(1, resultats$se, 0),xlab="False Positive Fraction (1-specificity)",ylab="True Positive Fraction (sensitivity)", type="n")


lines(c(1, resultats$sp1, 0), c(1, resultats$se, 0), type="s", lwd=2, col="red3")
abline(0,1, lty=2)


# Calcul de l'AUC
resultats[1:3,]
resultats<-resultats[order(resultats$sp1, resultats$se),]
resultats[dim(resultats)[1]+1,1]<-1
resultats[dim(resultats)[1],2]<-1
resultats[dim(resultats)[1],3]<-resultats[dim(resultats)[1]-1,3]
resultats[1:3,]

AUC <- sum( (resultats$sp1[2:length(resultats$sp1)] - resultats$sp1[1:(length(resultats$sp1)-1)] ) * (resultats$se[2:length(resultats$se)]) )
AUC

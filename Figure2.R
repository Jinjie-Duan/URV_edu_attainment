##############################
######### SWEDEN-WES #########
##############################

disrHBE <- read.table("SWE.pHIEbrain.disruptive.score", stringsAsFactor=F)
synoA <- read.table("SWE.synonymous.score", stringsAsFactor=F)
synoHBE <- read.table("SWE.pHIEbrain.synonymous.score", stringsAsFactor=F)
damaHBE <- read.table("SWE.pHIEbrainM.damaging.score", stringsAsFactor=F)


disrABE <- read.table("SWE.PALLEBRAIN.disruptive.score", stringsAsFactor=F)
synoABE <- read.table("SWE.PALLEBRAIN.synonymous.score", stringsAsFactor=F)
damaABE <- read.table("SWE.PALLEBRAIN.damaging.score", stringsAsFactor=F)


disrLBE <- read.table("SWE.pLOEbrain.disruptive.score", stringsAsFactor=F)
synoLBE <- read.table("SWE.pLOEbrain.synonymous.score", stringsAsFactor=F)
damaLBE <- read.table("SWE.pLOEbrainM.damaging.score", stringsAsFactor=F)





scoreSWE <- cbind(disrHBE,synoA[,2],synoHBE[,2],damaHBE[,2],disrABE[,2],synoABE[,2],damaABE[,2],disrLBE[,2],synoLBE[,2],damaLBE[,2])
colnames(scoreSWE) <- c("ID","disrHBE","synoA","synoHBE","damaHBE","disrABE","synoABE","damaABE","disrLBE","synoLBE","damaLBE")

sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)

## Merging datasets ##
eduswe <- merge(scoreSWE, sam, by="ID")



###############################
######### ESTONIA-WGS #########
###############################

disrHBE <- read.table("EST.pHIEbrain.disruptive.score", stringsAsFactor=F)
synoA <- read.table("EST.synonymous.score", stringsAsFactor=F)
synoHBE <- read.table("EST.pHIEbrain.synonymous.score", stringsAsFactor=F)
damaHBE <- read.table("EST.pHIEbrainM.damaging.score", stringsAsFactor=F)


disrABE <- read.table("EST.PALLEBRAIN.disruptive.score", stringsAsFactor=F)
synoABE <- read.table("EST.PALLEBRAIN.synonymous.score", stringsAsFactor=F)
damaABE <- read.table("EST.PALLEBRAIN.damaging.score", stringsAsFactor=F)


disrLBE <- read.table("EST.pLOEbrain.disruptive.score", stringsAsFactor=F)
synoLBE <- read.table("EST.pLOEbrain.synonymous.score", stringsAsFactor=F)
damaLBE <- read.table("EST.pLOEbrainM.damaging.score", stringsAsFactor=F)


scoreEST <- cbind(disrHBE,synoA[,2],synoHBE[,2],damaHBE[,2],disrABE[,2],synoABE[,2],damaABE[,2],disrLBE[,2],synoLBE[,2],damaLBE[,2])
colnames(scoreEST) <- c("ID","disrHBE","synoA","synoHBE","damaHBE","disrABE","synoABE","damaABE","disrLBE","synoLBE","damaLBE")

sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)

## Merging datasets ##
eduest <- merge(scoreEST, sam, by="ID")



###############################
######### FINLAND-WGS #########
###############################

disrHBE <- read.table("FINW.pHIEbrain.disruptive.score", stringsAsFactor=F)
synoA <- read.table("FINW.synonymous.score", stringsAsFactor=F)
synoHBE <- read.table("FINW.pHIEbrain.synonymous.score", stringsAsFactor=F)
damaHBE <- read.table("FINW.pHIEbrainM.damaging.score", stringsAsFactor=F)


disrABE <- read.table("FINW.PALLEBRAIN.disruptive.score", stringsAsFactor=F)
synoABE <- read.table("FINW.PALLEBRAIN.synonymous.score", stringsAsFactor=F)
damaABE <- read.table("FINW.PALLEBRAIN.damaging.score", stringsAsFactor=F)


disrLBE <- read.table("FINW.pLOEbrain.disruptive.score", stringsAsFactor=F)
synoLBE <- read.table("FINW.pLOEbrain.synonymous.score", stringsAsFactor=F)
damaLBE <- read.table("FINW.pLOEbrainM.damaging.score", stringsAsFactor=F)


scoreFINW <- cbind(disrHBE,synoA[,2],synoHBE[,2],damaHBE[,2],disrABE[,2],synoABE[,2],damaABE[,2],disrLBE[,2],synoLBE[,2],damaLBE[,2])
colnames(scoreFINW) <- c("ID","disrHBE","synoA","synoHBE","damaHBE","disrABE","synoABE","damaABE","disrLBE","synoLBE","damaLBE")


sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)

## Merging datasets ##
edufinw <- merge(scoreFINW, sam, by="ID")



###############################
######### FINLAND-WES #########
###############################

disrHBE <- read.table("FIN.pHIEbrain.disruptive.score", stringsAsFactor=F)
synoA <- read.table("FIN.synonymous.score", stringsAsFactor=F)
synoHBE <- read.table("FIN.pHIEbrain.synonymous.score", stringsAsFactor=F)
damaHBE <- read.table("FIN.pHIEbrainM.damaging.score", stringsAsFactor=F)


disrABE <- read.table("FIN.PALLEBRAIN.disruptive.score", stringsAsFactor=F)
synoABE <- read.table("FIN.PALLEBRAIN.synonymous.score", stringsAsFactor=F)
damaABE <- read.table("FIN.PALLEBRAIN.damaging.score", stringsAsFactor=F)


disrLBE <- read.table("FIN.pLOEbrain.disruptive.score", stringsAsFactor=F)
synoLBE <- read.table("FIN.pLOEbrain.synonymous.score", stringsAsFactor=F)
damaLBE <- read.table("FIN.pLOEbrainM.damaging.score", stringsAsFactor=F)


scoreFIN <- cbind(disrHBE,synoA[,2],synoHBE[,2],damaHBE[,2],disrABE[,2],synoABE[,2],damaABE[,2],disrLBE[,2],synoLBE[,2],damaLBE[,2])
colnames(scoreFIN) <- c("ID","disrHBE","synoA","synoHBE","damaHBE","disrABE","synoABE","damaABE","disrLBE","synoLBE","damaLBE")


sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)

## Merging datasets ##
edufin <- merge(scoreFIN, sam, by="ID")




######################################
######### ASSOCIATION ANALYSIS #######
######################################

## SWEDEN ##
mod0SWEHBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoHBE+yob2+yob3+scz, data=eduswe)
mod1SWEHBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrHBE+yob2+yob3+scz, data=eduswe)
mod2SWEHBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaHBE+yob2+yob3+scz, data=eduswe)

mod0SWEABE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoABE+yob2+yob3+scz, data=eduswe)
mod1SWEABE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrABE+yob2+yob3+scz, data=eduswe)
mod2SWEABE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaABE+yob2+yob3+scz, data=eduswe)

mod0SWELBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoLBE+yob2+yob3+scz, data=eduswe)
mod1SWELBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrLBE+yob2+yob3+scz, data=eduswe)
mod2SWELBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaLBE+yob2+yob3+scz, data=eduswe)



## ESTONIA ##
mod0ESTHBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoHBE+yob2+yob3, data=eduest)
mod1ESTHBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrHBE+yob2+yob3, data=eduest)
mod2ESTHBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaHBE+yob2+yob3, data=eduest)

mod0ESTABE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoABE+yob2+yob3, data=eduest)
mod1ESTABE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrABE+yob2+yob3, data=eduest)
mod2ESTABE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaABE+yob2+yob3, data=eduest)

mod0ESTLBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoLBE+yob2+yob3, data=eduest)
mod1ESTLBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrLBE+yob2+yob3, data=eduest)
mod2ESTLBE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaLBE+yob2+yob3, data=eduest)



## FINLAND WGS ##
mod0FINWHBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoHBE+yob2+yob3+scz, data=edufinw)
mod1FINWHBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrHBE+yob2+yob3+scz, data=edufinw)
mod2FINWHBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaHBE+yob2+yob3+scz, data=edufinw)

mod0FINWABE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoABE+yob2+yob3+scz, data=edufinw)
mod1FINWABE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrABE+yob2+yob3+scz, data=edufinw)
mod2FINWABE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaABE+yob2+yob3+scz, data=edufinw)

mod0FINWLBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoLBE+yob2+yob3+scz, data=edufinw)
mod1FINWLBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrLBE+yob2+yob3+scz, data=edufinw)
mod2FINWLBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaLBE+yob2+yob3+scz, data=edufinw)



## FINLAND WES ##
mod0FINHBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoHBE+yob2+yob3, data=edufin)
mod1FINHBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrHBE+yob2+yob3, data=edufin)
mod2FINHBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaHBE+yob2+yob3, data=edufin)

mod0FINABE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoABE+yob2+yob3, data=edufin)
mod1FINABE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrABE+yob2+yob3, data=edufin)
mod2FINABE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaABE+yob2+yob3, data=edufin)

mod0FINLBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoLBE+yob2+yob3, data=edufin)
mod1FINLBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrLBE+yob2+yob3, data=edufin)
mod2FINLBE <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaLBE+yob2+yob3, data=edufin)



### META  HIGHLY CONSTRAINED ###

library(meta)
beta <- c(as.numeric(coef(mod0SWEHBE)[15]),as.numeric(coef(mod0ESTHBE)[15]),as.numeric(coef(mod0FINWHBE)[15]),as.numeric(coef(mod0FINHBE)[15]))
se <- c(coef(summary(mod0SWEHBE))[15,2],coef(summary(mod0ESTHBE))[15,2],coef(summary(mod0FINWHBE))[15,2],coef(summary(mod0FINHBE))[15,2])
synometH <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))

beta <- c(as.numeric(coef(mod1SWEHBE)[15]),as.numeric(coef(mod1ESTHBE)[15]),as.numeric(coef(mod1FINWHBE)[15]),as.numeric(coef(mod1FINHBE)[15]))
se <- c(coef(summary(mod1SWEHBE))[15,2],coef(summary(mod1ESTHBE))[15,2],coef(summary(mod1FINWHBE))[15,2],coef(summary(mod1FINHBE))[15,2])
distmetH <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))

beta <- c(as.numeric(coef(mod2SWEHBE)[15]),as.numeric(coef(mod2ESTHBE)[15]),as.numeric(coef(mod2FINWHBE)[15]),as.numeric(coef(mod2FINHBE)[15]))
se <- c(coef(summary(mod2SWEHBE))[15,2],coef(summary(mod2ESTHBE))[15,2],coef(summary(mod2FINWHBE))[15,2],coef(summary(mod2FINHBE))[15,2])
damametH <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))



beta <- c(as.numeric(coef(mod0SWELBE)[15]),as.numeric(coef(mod0ESTLBE)[15]),as.numeric(coef(mod0FINWLBE)[15]),as.numeric(coef(mod0FINLBE)[15]))
se <- c(coef(summary(mod0SWELBE))[15,2],coef(summary(mod0ESTLBE))[15,2],coef(summary(mod0FINWLBE))[15,2],coef(summary(mod0FINLBE))[15,2])
synometL <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))

beta <- c(as.numeric(coef(mod1SWELBE)[15]),as.numeric(coef(mod1ESTLBE)[15]),as.numeric(coef(mod1FINWLBE)[15]),as.numeric(coef(mod1FINLBE)[15]))
se <- c(coef(summary(mod1SWELBE))[15,2],coef(summary(mod1ESTLBE))[15,2],coef(summary(mod1FINWLBE))[15,2],coef(summary(mod1FINLBE))[15,2])
distmetL <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))

beta <- c(as.numeric(coef(mod2SWELBE)[15]),as.numeric(coef(mod2ESTLBE)[15]),as.numeric(coef(mod2FINWLBE)[15]),as.numeric(coef(mod2FINLBE)[15]))
se <- c(coef(summary(mod2SWELBE))[15,2],coef(summary(mod2ESTLBE))[15,2],coef(summary(mod2FINWLBE))[15,2],coef(summary(mod2FINLBE))[15,2])
damametL <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))


beta <- c(as.numeric(coef(mod0SWEABE)[15]),as.numeric(coef(mod0ESTABE)[15]),as.numeric(coef(mod0FINWABE)[15]),as.numeric(coef(mod0FINABE)[15]))
se <- c(coef(summary(mod0SWEABE))[15,2],coef(summary(mod0ESTABE))[15,2],coef(summary(mod0FINWABE))[15,2],coef(summary(mod0FINABE))[15,2])
synometB <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))

beta <- c(as.numeric(coef(mod1SWEABE)[15]),as.numeric(coef(mod1ESTABE)[15]),as.numeric(coef(mod1FINWABE)[15]),as.numeric(coef(mod1FINABE)[15]))
se <- c(coef(summary(mod1SWEABE))[15,2],coef(summary(mod1ESTABE))[15,2],coef(summary(mod1FINWABE))[15,2],coef(summary(mod1FINABE))[15,2])
distmetB <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))

beta <- c(as.numeric(coef(mod2SWEABE)[15]),as.numeric(coef(mod2ESTABE)[15]),as.numeric(coef(mod2FINW)[15]),as.numeric(coef(mod2FINABE)[15]))
se <- c(coef(summary(mod2SWEABE))[15,2],coef(summary(mod2ESTABE))[15,2],coef(summary(mod2FINWABE))[15,2],coef(summary(mod2FINABE))[15,2])
damametB <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))





lhci <- function(x)
{
  co <- coef(x)[15]
  hi <- coef(x)[15] + 1.96 * coef(summary(x))[15,2]
  li <- coef(x)[15] - 1.96 * coef(summary(x))[15,2]
  p <- coef(summary(x))[15,4]
  return(c(co,hi,li,p))
}


df <- data.frame(
  y=c(
    synometH$TE.fixed,synometL$TE.fixed,synometB$TE.fixed,
    distmetH$TE.fixed,distmetL$TE.fixed,distmetB$TE.fixed,
    damametH$TE.fixed,damametL$TE.fixed,damametB$TE.fixed
    ),
  ymin=c(
  	synometH$lower.fixed,synometL$lower.fixed,synometB$lower.fixed,
    distmetH$lower.fixed,distmetL$lower.fixed,distmetB$lower.fixed,
    damametH$lower.fixed,damametL$lower.fixed,damametB$lower.fixed
    ),
  ymax=c(
  	synometH$upper.fixed,synometL$upper.fixed,synometB$upper.fixed,
    distmetH$upper.fixed,distmetL$upper.fixed,distmetB$upper.fixed,
    damametH$upper.fixed,damametL$upper.fixed,damametB$upper.fixed
    ),
  x=rev(c(1,2,3)),
  group=c(rep(c("Brain-expressed HC genes","Non brain-expressed HC genes","All brain-expressed genes"),3)),
  group1=c(rep("Synonymous",3),rep("Disruptive",3),rep("Damaging",3)),
  pval=c(
    synometH$pval.fixed,synometL$pval.fixed,synometB$pval.fixed,
    distmetH$pval.fixed,distmetL$pval.fixed,distmetB$pval.fixed,
    damametH$pval.fixed,damametL$pval.fixed,damametB$pval.fixed
    )
  )



df$group1 <- factor(df$group1, level=c("Disruptive","Damaging","Synonymous"))
library(grid)
library(ggplot2)
pdf("Figure2.pdf",width=8, height=6)
ggplot(df, aes(x=x, y=y, ymin=ymin, ymax=ymax)) +
  geom_hline(yintercept=0, color='red') + # se vuoi una linea verticale nel plot, se no togli
  geom_point(shape=23,size=6) + 
  geom_errorbar(width = .2) + # questo usa ymin e ymax
  theme_bw() +
  coord_flip() + # mette tutto orizzontale
  scale_x_continuous('',breaks=rev(c(1,2,3)), labels=c("Brain-expressed HC genes","Non brain-expressed HC genes","All brain-expressed genes")) + # toglie l'x-axis label
  scale_y_continuous('Change in years of education for 1 additional mutation') + # toglie il y-axis label
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(pval, digits=2))), hjust=-.1, vjust=-.8, size=3, parse=TRUE) + facet_wrap(~group1,ncol=1)  + theme(legend.position="none", plot.title = element_text(size=20, face="bold",  lineheight=0.6)) + expand_limits(x = c(1, 3.5))
dev.off()



##############################
######### SWEDEN-WES #########
##############################

disr <- read.table("SWE.pHI09.disruptive.score", stringsAsFactor=F)
synoA <- read.table("SWE.synonymous.score", stringsAsFactor=F)
syno <- read.table("SWE.pHI09.synonymous.score", stringsAsFactor=F)
dama <- read.table("SWE.pMISS309.damaging.score", stringsAsFactor=F)


disrL <- read.table("SWE.pLOW01.disruptive.score", stringsAsFactor=F)
synoL <- read.table("SWE.pLOW01.synonymous.score", stringsAsFactor=F)
damaL <- read.table("SWE.pLOW01.damaging.score", stringsAsFactor=F)


scoreSWE <- cbind(disr,synoA[,2],syno[,2],dama[,2],disrL[,2],synoL[,2],damaL[,2])
colnames(scoreSWE) <- c("ID","disr","synoA","syno","dama","disrL","synoL","damaL")

sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)

## Merging datasets ##
eduswe <- merge(scoreSWE, sam, by="ID")



###############################
######### ESTONIA-WGS #########
###############################

disr <- read.table("EST.pHI09.disruptive.score", stringsAsFactor=F)
synoA <- read.table("EST.synonymous.score", stringsAsFactor=F)
syno <- read.table("EST.pHI09.synonymous.score", stringsAsFactor=F)
dama <- read.table("EST.pMISS309.damaging.score", stringsAsFactor=F)


disrL <- read.table("EST.pLOW01.disruptive.score", stringsAsFactor=F)
synoL <- read.table("EST.pLOW01.synonymous.score", stringsAsFactor=F)
damaL <- read.table("EST.pLOW01.damaging.score", stringsAsFactor=F)


scoreEST <- cbind(disr,synoA[,2],syno[,2],dama[,2],disrL[,2],synoL[,2],damaL[,2])
colnames(scoreEST) <- c("ID","disr","synoA","syno","dama","disrL","synoL","damaL")

sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)

## Merging datasets ##
eduest <- merge(scoreEST, sam, by="ID")



###############################
######### FINLAND-WGS #########
###############################

disr <- read.table("FINW.pHI09.disruptive.score", stringsAsFactor=F)
synoA <- read.table("FINW.synonymous.score", stringsAsFactor=F)
syno <- read.table("FINW.pHI09.synonymous.score", stringsAsFactor=F)
dama <- read.table("FINW.pMISS309.damaging.score", stringsAsFactor=F)


disrL <- read.table("FINW.pLOW01.disruptive.score", stringsAsFactor=F)
synoL <- read.table("FINW.pLOW01.synonymous.score", stringsAsFactor=F)
damaL <- read.table("FINW.pLOW01.damaging.score", stringsAsFactor=F)


scoreFINW <- cbind(disr,synoA[,2],syno[,2],dama[,2],disrL[,2],synoL[,2],damaL[,2])
colnames(scoreFINW) <- c("ID","disr","synoA","syno","dama","disrL","synoL","damaL")

sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)

## Merging datasets ##
edufinw <- merge(scoreFINW, sam, by="ID")



###############################
######### FINLAND-WES #########
###############################

disr <- read.table("FIN.pHI09.disruptive.score", stringsAsFactor=F)
synoA <- read.table("FIN.synonymous.score", stringsAsFactor=F)
syno <- read.table("FIN.pHI09.synonymous.score", stringsAsFactor=F)
dama <- read.table("FIN.pMISS309.damaging.score", stringsAsFactor=F)


disrL <- read.table("FIN.pLOW01.disruptive.score", stringsAsFactor=F)
synoL <- read.table("FIN.pLOW01.synonymous.score", stringsAsFactor=F)
damaL <- read.table("FIN.pLOW01.damaging.score", stringsAsFactor=F)


scoreFIN <- cbind(disr,synoA[,2],syno[,2],dama[,2],disrL[,2],synoL[,2],damaL[,2])
colnames(scoreFIN) <- c("ID","disr","synoA","syno","dama","disrL","synoL","damaL")

sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)

## Merging datasets ##
edufin <- merge(scoreFIN, sam, by="ID")




######################################
######### ASSOCIATION ANALYSIS #######
######################################

## SWEDEN ##
mod0SWE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+syno+yob2+yob3+scz, data=eduswe)
mod1SWE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disr+yob2+yob3+scz, data=eduswe)
mod2SWE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+dama+yob2+yob3+scz, data=eduswe)

mod0SWEL <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoL+yob2+yob3+scz, data=eduswe)
mod1SWEL <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrL+yob2+yob3+scz, data=eduswe)
mod2SWEL <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaL+yob2+yob3+scz, data=eduswe)


## ESTONIA ##
mod0EST <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+syno+yob2+yob3, data=eduest)
mod1EST <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disr+yob2+yob3, data=eduest)
mod2EST <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+dama+yob2+yob3, data=eduest)

mod0ESTL <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoL+yob2+yob3, data=eduest)
mod1ESTL <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrL+yob2+yob3, data=eduest)
mod2ESTL <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaL+yob2+yob3, data=eduest)


## FINLAND WGS ##
mod0FINW <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+syno+yob2+yob3+scz, data=edufinw)
mod1FINW <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disr+yob2+yob3+scz, data=edufinw)
mod2FINW <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+dama+yob2+yob3+scz, data=edufinw)

mod0FINWL <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoL+yob2+yob3+scz, data=edufinw)
mod1FINWL <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrL+yob2+yob3+scz, data=edufinw)
mod2FINWL <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaL+yob2+yob3+scz, data=edufinw)



## FINLAND WES ##
mod0FIN <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+syno+yob2+yob3, data=edufin)
mod1FIN <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disr+yob2+yob3, data=edufin)
mod2FIN <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+dama+yob2+yob3, data=edufin)

mod0FINL <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+synoL+yob2+yob3, data=edufin)
mod1FINL <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+disrL+yob2+yob3, data=edufin)
mod2FINL <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+damaL+yob2+yob3, data=edufin)



### META  HIGHLY CONSTRAINED ###

library(meta)
beta <- c(as.numeric(coef(mod0SWE)[15]),as.numeric(coef(mod0EST)[15]),as.numeric(coef(mod0FINW)[15]),as.numeric(coef(mod0FIN)[15]))
se <- c(coef(summary(mod0SWE))[15,2],coef(summary(mod0EST))[15,2],coef(summary(mod0FINW))[15,2],coef(summary(mod0FIN))[15,2])
synomet <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))

beta <- c(as.numeric(coef(mod1SWE)[15]),as.numeric(coef(mod1EST)[15]),as.numeric(coef(mod1FINW)[15]),as.numeric(coef(mod1FIN)[15]))
se <- c(coef(summary(mod1SWE))[15,2],coef(summary(mod1EST))[15,2],coef(summary(mod1FINW))[15,2],coef(summary(mod1FIN))[15,2])
distmet <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))

beta <- c(as.numeric(coef(mod2SWE)[15]),as.numeric(coef(mod2EST)[15]),as.numeric(coef(mod2FINW)[15]),as.numeric(coef(mod2FIN)[15]))
se <- c(coef(summary(mod2SWE))[15,2],coef(summary(mod2EST))[15,2],coef(summary(mod2FINW))[15,2],coef(summary(mod2FIN))[15,2])
damamet <- metagen(beta,se,studlab=c("SwedenCA","Estonian","FinWGS","FinWES"))



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
    lhci(mod0SWE)[1],lhci(mod0EST)[1],lhci(mod0FINW)[1],lhci(mod0FIN)[1],synomet$TE.fixed,
    lhci(mod1SWE)[1],lhci(mod1EST)[1],lhci(mod1FINW)[1],lhci(mod1FIN)[1],distmet$TE.fixed,
    lhci(mod2SWE)[1],lhci(mod2EST)[1],lhci(mod2FINW)[1],lhci(mod2FIN)[1],damamet$TE.fixed
    ),
  ymin=c(
    lhci(mod0SWE)[3],lhci(mod0EST)[3],lhci(mod0FINW)[3],lhci(mod0FIN)[3],synomet$lower.fixed,
    lhci(mod1SWE)[3],lhci(mod1EST)[3],lhci(mod1FINW)[3],lhci(mod1FIN)[3],distmet$lower.fixed,
    lhci(mod2SWE)[3],lhci(mod2EST)[3],lhci(mod2FINW)[3],lhci(mod2FIN)[3],damamet$lower.fixed
    ),
  ymax=c(
    lhci(mod0SWE)[2],lhci(mod0EST)[2],lhci(mod0FINW)[2],lhci(mod0FIN)[2],synomet$upper.fixed,
    lhci(mod1SWE)[2],lhci(mod1EST)[2],lhci(mod1FINW)[2],lhci(mod1FIN)[2],distmet$upper.fixed,
    lhci(mod2SWE)[2],lhci(mod2EST)[2],lhci(mod2FINW)[2],lhci(mod2FIN)[2],damamet$upper.fixed
    ),
  x=rev(c(1,2,3,4,5)),
  group=c(rep(c("Swedish-WES","Estonian-WGS","Finnish-WGS","Finnish-WES","Meta-analysis"),3)),
  group1=c(rep("Synonymous",5),rep("Disruptive",5),rep("Damaging",5)),
  shape=rep(c(15,15,15,15,23),3),
  size=rep(c(length(mod0SWE$y),length(mod0EST$y),length(mod0FINW$y),length(mod0FIN$y),(length(mod0SWE$y)+length(mod0EST$y)+length(mod0FINW$y)+length(mod0FIN$y))/2),3),
  pval=c(
    lhci(mod0SWE)[4],lhci(mod0EST)[4],lhci(mod0FINW)[4],lhci(mod0FIN)[4],synomet$pval.fixed,
    lhci(mod1SWE)[4],lhci(mod1EST)[4],lhci(mod1FINW)[4],lhci(mod1FIN)[4],distmet$pval.fixed,
    lhci(mod2SWE)[4],lhci(mod2EST)[4],lhci(mod2FINW)[4],lhci(mod2FIN)[4],damamet$pval.fixed
    )
  )



df$group1 <- factor(df$group1, level=c("Disruptive","Damaging","Synonymous"))
library(grid)
library(ggplot2)
pdf("Figure1.pdf",width=7)
ggplot(df, aes(x=x, y=y, ymin=ymin, ymax=ymax)) +
  geom_hline(yintercept=0, color='red') + # se vuoi una linea verticale nel plot, se no togli
  geom_point(aes(size=size, shape=shape)) + scale_shape_identity() +
  geom_errorbar(width = .2) + # questo usa ymin e ymax
  theme_bw() +
  coord_flip() + # mette tutto orizzontale
  scale_x_continuous('',breaks=rev(c(1,2,3,4,5)), labels=c("Swedish-WES","Estonian-WGS","Finnish-WGS","Finnish-WES","Meta-analysis")) + # toglie l'x-axis label
  scale_y_continuous('Change in years of education for 1 additional mutation') + # toglie il y-axis label
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(pval, digits=2))), hjust=-.1, vjust=-1.2, size=3, parse=TRUE) + facet_wrap(~group1,ncol=1)  + theme(legend.position="none", plot.title = element_text(size=20, face="bold",  lineheight=0.6)) + expand_limits(x = c(1, 5.7))
dev.off()



##############################
######### SWEDEN-WES #########
##############################

disr <- read.table("SWE.pHI09.disruptive.score", stringsAsFactor=F)
synoA <- read.table("SWE.synonymous.score", stringsAsFactor=F)
syno <- read.table("SWE.pHI09.synonymous.score", stringsAsFactor=F)
dama <- read.table("SWE.pMISS309.damaging.score", stringsAsFactor=F)

scoreSWE <- cbind(disr,synoA[,2],syno[,2],dama[,2])
colnames(scoreSWE) <- c("ID","disr","synoA","syno","dama")

sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)
polySWE <- read.table("common_var_polygene.profile", stringsAsFactor=F, header=T)
cnvSWE <- read.table("cnv_per_individual_swe", stringsAsFactor=F, header=T)
rohSWE <- read.table("roh_per_individual_swe", stringsAsFactor=F, header=T)


## Merging datasets ##
eduswe00 <- merge(cnvSWE,rohSWE,by="ID",all.x=T,all.y=T)
eduswe0 <- merge(eduswe00,polySWE,by="ID",all.x=T,all.y=T)
eduswe <- merge(eduswe0, sam, by="ID")



###############################
######### ESTONIA-WGS #########
###############################

disr <- read.table("EST.pHI09.disruptive.score", stringsAsFactor=F)
synoA <- read.table("EST.synonymous.score", stringsAsFactor=F)
syno <- read.table("EST.pHI09.synonymous.score", stringsAsFactor=F)
dama <- read.table("EST.pMISS309.damaging.score", stringsAsFactor=F)


scoreEST <- cbind(disr,synoA[,2],syno[,2],dama[,2])
colnames(scoreEST) <- c("ID","disr","synoA","syno","dama")

sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)
cnvEST <- read.table("cnv_per_individual_est", stringsAsFactor=F, header=T)
rohEST <- read.table("roh_per_individual_est", stringsAsFactor=F, header=T)


## Merging datasets ##
eduest0 <- merge(cnvEST,rohEST,by="ID",all.x=T,all.y=T)
eduest <- merge(eduest0, sam, by="ID")



###############################
######### FINLAND-WGS #########
###############################

disr <- read.table("FINW.pHI09.disruptive.score", stringsAsFactor=F)
synoA <- read.table("FINW.synonymous.score", stringsAsFactor=F)
syno <- read.table("FINW.pHI09.synonymous.score", stringsAsFactor=F)
dama <- read.table("FINW.pMISS309.damaging.score", stringsAsFactor=F)


scoreFINW <- cbind(disr,synoA[,2],syno[,2],dama[,2])
colnames(scoreFINW) <- c("ID","disr","synoA","syno","dama")

sam <-read.table("file_with_phenotypes.tsv", sep="\t", stringsAsFactor=F, header=T)
cnvFINW <- read.table("cnv_per_individual_finw", stringsAsFactor=F, header=T)
rohFINW <- read.table("roh_per_individual_finw", stringsAsFactor=F, header=T)


## Merging datasets ##
edufinw0 <- merge(cnvFINW,rohFINW,by="ID",all.x=T,all.y=T)
edufinw <- merge(edufinw0, sam, by="ID")





### Sweden ###
eduswe$delS <- scale(eduswe$disr+eduswe$dama)
eduswe$EDUPGSS <- -scale(eduswe$CNT)
eduswe$CNVS <- scale(eduswe$CNV)
eduswe$ROHS <- scale(eduswe$ROH)


mod0polySWE <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+delS+CNVS+ROHS+EDUPGSS+yob2+yob3+ARRAY_TYPE+scz, data=eduswe)


### Estonia ###
eduest$delS <- scale(eduest$disr+eduest$dama)
eduest$CNVS <- scale(eduest$CNV)
eduest$ROHS <- scale(eduest$ROH)


mod0polyEST <- glm(yearsofedu ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+delS+CNVS+ROHS+yob2+yob3, data=eduest)


### Finland ###
edufinw$delS <- scale(edufinw$disr+edufinw$dama)
edufinw$CNVS <- scale(edufinw$CNV)
edufinw$ROHS <- scale(edufinw$ROH)

mod0polyFIN <- glm(yearsofschool ~ sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+yob+synoA+delS+CNVS+ROHS+yob2+yob3+scz, data=edufinw)



library(meta)
beta <- c(as.numeric(coef(mod0polySWE)[15]),as.numeric(coef(mod0polyEST)[15]),as.numeric(coef(mod0polyFIN)[15]))
se <- c(coef(summary(mod0polySWE))[15,2],coef(summary(mod0polyEST))[15,2],coef(summary(mod0polyFIN))[16,2])
delSmet <- metagen(beta,se,studlab=c("Sweden","Estonian","Finland"))

beta <- c(as.numeric(coef(mod0polySWE)[16]),as.numeric(coef(mod0polyEST)[16]),as.numeric(coef(mod0polyFIN)[16]))
se <- c(coef(summary(mod0polySWE))[16,2],coef(summary(mod0polyEST))[16,2],coef(summary(mod0polyFIN))[16,2])
CNVsmet <- metagen(beta,se,studlab=c("Sweden","Estonian","Finland"))

beta <- c(as.numeric(coef(mod0polySWE)[17]),as.numeric(coef(mod0polyEST)[17]),as.numeric(coef(mod0polyFIN)[17]))
se <- c(coef(summary(mod0polySWE))[17,2],coef(summary(mod0polyEST))[17,2],coef(summary(mod0polyFIN))[17,2])
ROHSmet <- metagen(beta,se,studlab=c("Sweden","Estonian","Finland"))



lhci <- function(x)
{
  co <- coef(x)[18]
  hi <- coef(x)[18] + 1.96 * coef(summary(x))[18,2]
  li <- coef(x)[18] - 1.96 * coef(summary(x))[18,2]
  p <- coef(summary(x))[18,4]
  return(c(co,hi,li,p))
}

df <- data.frame(
  y=c(
    lhci(mod0polySWE)[1],delSmet$TE.fixed,ROHSmet$TE.fixed,CNVsmet$TE.fixed
    ),
  ymin=c(
    lhci(mod0polySWE)[3],delSmet$lower.fixed,ROHSmet$lower.fixed,CNVsmet$lower.fixed
    ),
  ymax=c(
    lhci(mod0polySWE)[2],delSmet$upper.fixed,ROHSmet$upper.fixed,CNVsmet$upper.fixed
    ),
  x=rev(c(1,2,3,4)),
  group=c("Polygenic score","Disruptive + damaging ultra-rare","Run of homozygosity","Pathogenic CNVs"),
  shape=c(15,23,23,23),
  pval=c(
    lhci(mod0polySWE)[4],delSmet$pval.fixed,ROHSmet$pval.fixed,CNVsmet$pval.fixed
    )
  )



df$group <- factor(df$group, level=c("Polygenic score","Disruptive + damaging ultra-rare","Run of homozygosity","Pathogenic CNVs"))
#library(grid)
library(ggplot2)
pdf("Figure3.pdf",width=7, height=4)
p2 <- ggplot(df, aes(x=x, y=y, ymin=ymin, ymax=ymax)) +
  geom_hline(yintercept=0, color='dark gray') + # se vuoi una linea verticale nel plot, se no togli
  geom_point(size=5, aes(shape=shape)) + scale_shape_identity() +
  geom_errorbar(width = .2) + # questo usa ymin e ymax
  theme_bw() +
  coord_flip() + # mette tutto orizzontale
  scale_x_continuous('Normalized scores',breaks=rev(c(1,2,3,4)), labels=c("Polygenic score","Disruptive + damaging \n ultra-rare variants","Run of homozygosity","Pathogenic CNVs")) + # toglie l'x-axis label
  scale_y_continuous('Change in years of education for 1 SD change \n in the normalized scores') + # toglie il y-axis label
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(pval, digits=2))), hjust=-.1, vjust=-.8, size=3, parse=TRUE) + theme(legend.position="none", plot.title = element_text(size=20, face="bold",  lineheight=0.6)) + expand_limits(x = c(1, 4.5)) + ggtitle('') + theme(plot.title=element_text(hjust=0))
  p2
 dev.off()

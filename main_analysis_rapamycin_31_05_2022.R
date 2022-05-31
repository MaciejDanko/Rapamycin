rm(list=ls())

set.seed(123)
 
library(xlsx) #for write.xlsx
library(survival) #for Surv(), survFit()
library(magicaxis) #for magaxis()
library(mgcv) # for gam(), predict.gam(), s(), ...
library(coin) # for logrank_test() and independence_test(), ansari_trafo()
library(parallel) # for detectCores()
library(itsadug) # for acf_resid()
library(rms) ## for npsurv()

load(file='Rapamycin.rda')

mod<-''
skip <- -5
maxX <- 150
res <- 1200
BS <- 'ps'
K <- -1
X <- 0 : maxX
NCores <- detectCores()-1

################################################################################
# Non-parametric analysis - stratified logranks
# https://www.statisticshowto.com/log-rank-test/
################################################################################

calc.stat.coin<-function(object, name, tails){
  list(name=name,
       tails=tails,
       method=object@method,
       Z=object@statistic@teststatistic,
       P=(1-pnorm(abs(object@statistic@teststatistic)))*2)
}

RapamycinData$survdat$izoline<-as.factor(as.character(RapamycinData$survdat$izoline))
(STRATA_SEX_LR<-coin::logrank_test(Surv(midday,status) ~ sex | treatment + izoline, data = RapamycinData$survdat, type='logrank',
                                   ties.method='mid-ranks'))
(STRATA_SEX_GB<-coin::logrank_test(Surv(midday,status) ~ sex | treatment + izoline, data = RapamycinData$survdat, type='Gehan-Breslow',
                                   ties.method='Hothorn-Lausen'))
(STRATA_TR_LR<-coin::logrank_test(Surv(midday,status) ~ treatment | sex + izoline, data = RapamycinData$survdat, type='logrank',
                                  ties.method='mid-ranks'))
(STRATA_TR_GB<-coin::logrank_test(Surv(midday,status) ~ treatment | sex + izoline, data = RapamycinData$survdat, type='Gehan-Breslow',
                                  ties.method='Hothorn-Lausen'))

STRATA_SEX_LR<-calc.stat.coin(STRATA_SEX_LR,'Sex stratified by treatment and izoline','mid-ranks')
STRATA_SEX_GB<-calc.stat.coin(STRATA_SEX_GB,'Sex stratified by treatment and izoline','Hothorn-Lausen')
STRATA_TR_LR<-calc.stat.coin(STRATA_TR_LR,'Treatment stratified by sex and izoline','mid-ranks')
STRATA_TR_GB<-calc.stat.coin(STRATA_TR_GB,'Treatment stratified by sex and izoline','Hothorn-Lausen')

TAB1<-rbind(STRATA_SEX_LR,STRATA_SEX_GB,STRATA_TR_LR,STRATA_TR_GB)
xlsx::write.xlsx(TAB1,file=paste('./tables/TAB1_nonparametric_tests',mod,'.xlsx',sep=''),sheetName='Tests', col.names = TRUE, append = FALSE, row.names = TRUE)

################################################################################
# Asymptotic test against crossing-curve alternatives (Shen and Le, 2000)
################################################################################

if (FALSE) {
  plot(survfit(Surv(midday,status) ~ treatment, data = RapamycinData$survdat), col=1:2)
  shen_trafo <- function(x)
    ansari_trafo(logrank_trafo(x, type = "Prentice"))
  
  independence_test(Surv(midday,status) ~ treatment | sex + izoline, data = RapamycinData$survdat,
                    ytrafo = function(data)
                      trafo(data, surv_trafo = shen_trafo))
  independence_test(Surv(midday,status) ~ treatment + izoline, data = RapamycinData$survdat,
                    ytrafo = function(data)
                      trafo(data, surv_trafo = shen_trafo))
  
  plot(survfit(Surv(midday,status) ~ sex, data = RapamycinData$survdat), col=1:2)
  independence_test(Surv(midday,status) ~ sex , data = RapamycinData$survdat,
                    ytrafo = function(data)
                      trafo(data, surv_trafo = shen_trafo))
}

################################################################################
# Kaplan Meier plots
################################################################################
library(Cairo)
options(bitmapType="cairo")

#old colors
#c_red<-'#df536b'; c_blue<-'#2297e6'; c_yellow<-'#cd8500'; c_green<-'#008b00';
#new colors
c_red<-'#f36f83'; c_blue<-'#1550b0'; c_yellow<-'#d79a20'; c_green<-'#006b00';


#tiff(filename=paste('FIG1_Survivorship_KM',mod,'.tiff',sep=''),width=res*8,height=res*7,compression ='lzw',res=res,units='px')
pdf(paste('./figures/FIG1_Survivorship_KM',mod,'.pdf',sep=''),width=8,height=7)
#options(bitmapType="cairo")
par(mfrow=c(2,2))
par(mar=c(4,4,0,0),oma=c(1,1,0.15,0.15))
csf <- npsurv(Surv(midday,status) ~ treatment+izoline, data = RapamycinData$survdat[which(RapamycinData$survdat$sex=='F'),])
plot(csf,col=adjustcolor(c(rep(c_blue,15),rep(c_red,15)),alpha.f = 0.5),mark.time = !TRUE, conf.int = FALSE, xlab='', ylab='Survivorship',xlim=c(skip,maxX))
csf <- npsurv(Surv(midday,status) ~ treatment, data = RapamycinData$survdat[which(RapamycinData$survdat$sex=='F'),])
par(new=TRUE)
plot(csf,col=c(c_blue,c_red),mark.time = !TRUE, conf.int = FALSE, xlab='', ylab='',xlim=c(skip,maxX),lwd=3)
legend('topright',c('Control','Rapamycin'),col=c(c_blue,c_red),lty=1,bty='n',lwd=3)
legend('bottomleft','Females',bty='n',inset=c(-0.05,0))
legend('topleft',expression(bold(a)),bty='n',inset=c(-0.075,0),cex=1.3)

csf <- npsurv(Surv(midday,status) ~ treatment+izoline, data = RapamycinData$survdat[which(RapamycinData$survdat$sex=='M'),])
plot(csf,col=adjustcolor(c(rep(c_blue,15),rep(c_red,15)),alpha.f = 0.5),mark.time = !TRUE, conf.int = FALSE, xlab='', ylab='',xlim=c(skip,maxX))
csf <- npsurv(Surv(midday,status) ~ treatment, data = RapamycinData$survdat[which(RapamycinData$survdat$sex=='M'),])
par(new=TRUE)
plot(csf,col=c(c_blue,c_red),mark.time = !TRUE, conf.int = FALSE, xlab='', ylab='Survivorship',xlim=c(skip,maxX),lwd=3)
legend('bottomleft','Males',bty='n',inset=c(-0.05,0))
legend('topright',c('Control','Rapamycin'),col=c(c_blue,c_red),lty=1,bty='n',lwd=3)
legend('topleft',expression(bold(b)),bty='n',inset=c(-0.075,0),cex=1.3)

csf <- npsurv(Surv(midday,status) ~ sex+izoline, data = RapamycinData$survdat[which(RapamycinData$survdat$treatment=='C'),])
plot(csf,col=adjustcolor(c(rep(c_green,15),rep(c_yellow,15)),alpha.f = 0.5),mark.time = !TRUE, conf.int = FALSE, xlab='', ylab='',xlim=c(skip,maxX))
csf <- npsurv(Surv(midday,status) ~ sex, data = RapamycinData$survdat[which(RapamycinData$survdat$treatment=='C'),])
par(new=TRUE)
plot(csf,col=c(c_green,c_yellow),mark.time = !TRUE, conf.int = FALSE, xlab='Age', ylab='Survivorship',xlim=c(skip,maxX),lwd=3)
legend('bottomleft','Control',bty='n',inset=c(-0.05,0))
legend('topright',c('Females','Males'),col=c(c_green,c_yellow),lty=1,bty='n',lwd=3)
legend('topleft',expression(bold(c)),bty='n',inset=c(-0.075,0),cex=1.3)

csf <- npsurv(Surv(midday,status) ~ sex+izoline, data = RapamycinData$survdat[which(RapamycinData$survdat$treatment=='R'),])
plot(csf,col=adjustcolor(c(rep(c_green,15),rep(c_yellow,15)),alpha.f = 0.5),mark.time = !TRUE, conf.int = FALSE, xlab='', ylab='',xlim=c(skip,maxX))
csf <- npsurv(Surv(midday,status) ~ sex, data = RapamycinData$survdat[which(RapamycinData$survdat$treatment=='R'),])
par(new=TRUE)
plot(csf,col=c(c_green,c_yellow),mark.time = !TRUE, conf.int = FALSE, xlab='Age', ylab='',xlim=c(skip,maxX),lwd=3)
legend('bottomleft','Rapamycin',bty='n',inset=c(-0.05,0))
legend('topright',c('Females','Males'),col=c(c_green,c_yellow),lty=1,bty='n',lwd=3)
legend('topleft',expression(bold(d)),bty='n',inset=c(-0.075,0),cex=1.3)

dev.off()

################################################################################
# Mortality analysis using GAM - Life table construction 
################################################################################

LTM<-data.frame(data.table::rbindlist(lapply(unique(RapamycinData$data$sextrizoday),function(k){
  tmp<-RapamycinData$data[which(RapamycinData$data$sextrizoday==k),]
  data.frame(sex=unique(tmp$sex), treatment=unique(tmp$treatment), izoline=unique(tmp$izoline), 
             midday=unique(tmp$midday),dead=sum(tmp$dead, na.rm=TRUE), 
             nx=sum(tmp$NX.V),exposures=sum(tmp$exposuresV,na.rm = TRUE),
             stringsAsFactors = FALSE) 
})),stringsAsFactors = FALSE)

LTM$sextr<-paste(LTM$sex, LTM$treatment,sep='_')
LTM$sextrizo<-paste(LTM$sex, LTM$treatment, LTM$izoline, sep='_')
LTM$sextr<-as.factor(LTM$sextr)
LTM$sex<-as.factor(LTM$sex)
LTM$treatment<-as.factor(LTM$treatment)
LTM$izoline<-as.factor(LTM$izoline)
LTM$sextrizo<-as.factor(LTM$sextrizo)

ALLIZO<-unique(LTM$izoline)

NI1=sapply(ALLIZO, function(k) max(LTM[LTM$sex=='F'&LTM$treatment=='C'&LTM$izoline==k,'nx']))
NI2=sapply(ALLIZO, function(k) max(LTM[LTM$sex=='F'&LTM$treatment=='R'&LTM$izoline==k,'nx']))
NI3=sapply(ALLIZO, function(k) max(LTM[LTM$sex=='M'&LTM$treatment=='C'&LTM$izoline==k,'nx']))
NI4=sapply(ALLIZO, function(k) max(LTM[LTM$sex=='M'&LTM$treatment=='R'&LTM$izoline==k,'nx']))

BI1=sapply(ALLIZO, function(k) max(LTM[LTM$sex=='F'&LTM$treatment=='C'&LTM$izoline==k,'midday']))
BI2=sapply(ALLIZO, function(k) max(LTM[LTM$sex=='F'&LTM$treatment=='R'&LTM$izoline==k,'midday']))
BI3=sapply(ALLIZO, function(k) max(LTM[LTM$sex=='M'&LTM$treatment=='C'&LTM$izoline==k,'midday']))
BI4=sapply(ALLIZO, function(k) max(LTM[LTM$sex=='M'&LTM$treatment=='R'&LTM$izoline==k,'midday']))

################################################################################
# Fit GAM and GAMM models
################################################################################

if (!file.exists('Models.rda')){
  
  M_izoline_0_gamreml<-gam(dead~treatment*sex+s(midday,by=treatment:sex,bs=BS,k=K),#+
                           #select=TRUE,
                           offset=log(exposures),
                           data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  M_izoline_FS_0_gamreml<-gam(dead~treatment*sex+s(midday,by=treatment:sex,bs=BS,k=K)+
                                s(midday,izoline, bs='fs',k=K, m=1),
                              # select=TRUE,
                              offset=log(exposures),
                              data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  M_izoline_FS_sex_gamreml<-gam(dead~treatment*sex+s(midday,by=treatment:sex,bs=BS,k=K)+
                                  s(midday,izoline, by=sex, bs='fs',k=K, m=1),
                                #select=TRUE,
                                offset=log(exposures),
                                data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  
  M_izoline_FS_treatment_gamreml<-gam(dead~treatment*sex+s(midday,by=treatment:sex,bs=BS,k=K)+
                                        s(midday,izoline, by=treatment, bs='fs',k=K, m=1),
                                      #select=TRUE,
                                      offset=log(exposures),
                                      data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  
  M_izoline_FS_sex_treatment_gamreml<-gam(dead~treatment*sex+
                                            s(midday,by=treatment:sex,bs=BS,k=K)+
                                            s(midday,izoline, by=sex, bs='fs',k=K, m=1) +
                                            s(midday,izoline, by=treatment, bs='fs',k=K, m=1),
                                          #select=TRUE,
                                          offset=log(exposures),
                                          data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  
  M_izoline_RE_sex_treatment_inter_gamreml<-gam(dead~treatment*sex+
                                                  s(midday,by=treatment:sex,bs=BS,k=K)+
                                                  s(izoline, bs='re',k=K) +
                                                  s(izoline, midday, bs='re',k=K) +
                                                  s(izoline, sextr, bs='re',k=K),
                                                offset=log(exposures),
                                                data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  M_izoline_RE_sex_treatment_gamreml<-gam(dead~treatment*sex+
                                            s(midday,by=treatment:sex,bs=BS,k=K)+
                                            s(izoline, bs='re',k=K) +
                                            s(izoline, midday, bs='re',k=K) +
                                            s(izoline, sex, bs='re',k=K) +
                                            s(izoline, treatment, bs='re',k=K),
                                          offset=log(exposures),
                                          data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  M_izoline_RE_0_gamreml<-gam(dead~treatment*sex+
                                s(midday,by=treatment:sex,bs=BS,k=K)+
                                s(izoline, bs='re',k=K),# +
                              #select=TRUE,
                              offset=log(exposures),
                              data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  M_izoline_RE_age_gamreml<-gam(dead~treatment*sex+
                                  s(midday,by=treatment:sex,bs=BS,k=K)+
                                  s(izoline, bs='re',k=K) +
                                  s(izoline, midday, bs='re',k=K),
                                #select=TRUE,
                                offset=log(exposures),
                                data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  
  M_izoline_RE_sex_gamreml<-gam(dead~treatment*sex+
                                  s(midday,by=treatment:sex,bs=BS,k=K)+
                                  s(izoline, bs='re',k=K) +
                                  s(izoline, midday, bs='re',k=K) +
                                  s(izoline, sex, bs='re',k=K),
                                #select=TRUE,
                                offset=log(exposures),
                                data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  M_izoline_RE_treatment_gamreml<-gam(dead~treatment*sex+
                                        s(midday,by=treatment:sex,bs=BS,k=K)+
                                        s(izoline, bs='re',k=K) +
                                        s(izoline, midday, bs='re',k=K) +
                                        s(izoline, treatment, bs='re',k=K),
                                      #select=TRUE,
                                      offset=log(exposures),
                                      data=LTM,family = nb(), control = gam.control(nthreads=NCores), method='REML')
  
  
  save(list=c(
    'M_izoline_0_gamreml',
    'M_izoline_FS_0_gamreml',
    'M_izoline_FS_sex_gamreml',
    'M_izoline_FS_treatment_gamreml',
    'M_izoline_FS_sex_treatment_gamreml', # most parsimonious
    'M_izoline_RE_sex_treatment_gamreml',
    'M_izoline_RE_sex_treatment_inter_gamreml',
    'M_izoline_RE_0_gamreml',
    'M_izoline_RE_age_gamreml',
    'M_izoline_RE_sex_gamreml',
    'M_izoline_RE_treatment_gamreml'
  ),file='Models.rda')#,
  
} else load('Models.rda')

################################################################################
# Model selection
################################################################################

getFormula<-function(model, simplify=TRUE){
  res<-paste(deparse(formula(model), width.cutoff=200),collapse = '')
  if (simplify){
    res<-gsub(' ','',res)
    res<-gsub('dead','',res)
    res<-gsub('BS','"ps"',res,fixed = TRUE)
    res<-gsub(',k=K','',res)
    #res<-gsub(',m=1','',res)
    res<-gsub('sextr','TrSexInteraction',res)
    res<-gsub('sex','Sex',res)
    res<-gsub('treatment','Tr',res)
    res<-gsub('izoline','Iso',res)
    res<-gsub('midday','Age',res)
    res<-gsub('(Tr:Sex)','Tr:Sex',res,fixed = TRUE)
  }
  res
}

getAIC<-function(...){
  A<-round(AIC(...),2)
  rownames(A)<-sapply(list(...),function(k) getFormula(k, TRUE))
  A
}

A<-getAIC(
  M_izoline_0_gamreml,
  M_izoline_FS_0_gamreml,
  M_izoline_FS_sex_gamreml,
  M_izoline_FS_treatment_gamreml,
  M_izoline_FS_sex_treatment_gamreml, # most parsimonious
  M_izoline_RE_0_gamreml,
  M_izoline_RE_age_gamreml,
  M_izoline_RE_sex_gamreml,
  M_izoline_RE_treatment_gamreml,
  M_izoline_RE_sex_treatment_gamreml,
  M_izoline_RE_sex_treatment_inter_gamreml
)#,

Model_Selection<-A[order(A[,2]),]
m0<-M_izoline_FS_sex_treatment_gamreml

convert_terms<-function(tnames){
  tnames<-gsub('s(midday):','s(Age):',tnames,fixed = TRUE)
  tnames<-gsub('treatment:sexC:F',':Control:Females',tnames,fixed = TRUE)
  tnames<-gsub('treatment:sexC:M',':Control:Males',tnames,fixed = TRUE)
  tnames<-gsub('treatment:sexR:F',':Rapamycin:Females',tnames,fixed = TRUE)
  tnames<-gsub('treatment:sexR:M',':Rapamycin:Males',tnames,fixed = TRUE)
  tnames<-gsub('s(midday,izoline):','s(Age,Iso):',tnames,fixed = TRUE)
  tnames<-gsub(':treatmentC',':Control',tnames,fixed = TRUE)
  tnames<-gsub(':treatmentR',':Rapamycin',tnames,fixed = TRUE)
  tnames<-gsub(':sexF',':Females',tnames,fixed = TRUE)
  tnames<-gsub(':sexM',':Males',tnames,fixed = TRUE)
  tnames
}

convert_terms_2<-function(tnames){
  tnames<-gsub('R','(Rapamycin)',tnames,fixed = TRUE)
  tnames<-gsub('C','(Control)',tnames,fixed = TRUE)
  tnames<-gsub('M','(Males)',tnames,fixed = TRUE)
  tnames<-gsub('F','(Females)',tnames,fixed = TRUE)
  tnames<-gsub('treatment','Tr',tnames,fixed = TRUE)
  tnames<-gsub('sex','Sex',tnames,fixed = TRUE)
  tnames
}

convert_terms_3<-function(tnames){
  tnames<-gsub('R','(Rapamycin)',tnames,fixed = TRUE)
  tnames<-gsub('C','(Control)',tnames,fixed = TRUE)
  tnames<-gsub('M','(Males)',tnames,fixed = TRUE)
  tnames<-gsub('F','(Females)',tnames,fixed = TRUE)
  tnames<-gsub('treatment=','',tnames,fixed = TRUE)
  tnames<-gsub('sex=','',tnames,fixed = TRUE)
  tnames
}

################################################################################
# Exporting models summary to a file
################################################################################
make_sheet<-function(objectsum, aic){
  L_1<-as.matrix(cbind(t(as.matrix(getFormula(objectsum,TRUE))),'','',''))
  rownames(L_1)<-'Formula'
  L_2<-as.matrix(t(c('','','','')))
  rownames(L_2)<-''
  L_2b<-as.matrix(t(c(round(aic,2),'','','')))
  rownames(L_2b)<-'AIC'
  L_3<-t(as.matrix(c('','','','')))
  rownames(L_3)<-'Non-smooth terms:'
  L_5<-round(objectsum$p.table,4)
  rownames(L_5)<-convert_terms_2(rownames(L_5))
  L_4<-t(as.matrix(colnames(L_5)))
  L_6<-t(as.matrix(c('','','','')))
  rownames(L_6)<-''
  L_7<-t(as.matrix(c('','','','')))
  rownames(L_7)<-'Smooth terms:'
  L_9<-round(objectsum$s.table,4)
  rownames(L_9)<-convert_terms(rownames(L_9))
  L_8<-colnames(L_9)
  L<-rbind(L_1,' '=L_2,L_2b,' '=L_2,L_3,' '=L_4,' '=L_5,' '=L_6,L_7,' '=L_8,L_9)
  rownames(L)[nchar(rownames(L))==0]<-' '
  L<-cbind(rownames(L),L)
  colnames(L)<-NULL  
  rownames(L)<-NULL  
  data.frame(L,check.names = FALSE, 
             fix.empty.names = FALSE,
             stringsAsFactors = FALSE, 
             check.rows = FALSE)
}

L1<-make_sheet(summary(m0),AIC(m0))

xlsx::write.xlsx(Model_Selection,file=paste('./tables/TAB1Sab_GAM_summary',mod,'.xlsx',sep=''),sheetName='(a) Model selection', col.names = TRUE, append = FALSE, row.names = TRUE)
Sys.sleep(1)
xlsx::write.xlsx(L1,file=paste('./tables/TAB1Sab_GAM_summary',mod,'.xlsx',sep=''),sheetName='(b) Most parsimonious model', col.names = FALSE, append = TRUE, row.names = FALSE)
Sys.sleep(1)

if(FALSE){ # other models
  L2<-make_sheet(summary(M_izoline_0_gamreml),AIC(M_izoline_0_gamreml))
  L3<-make_sheet(summary(M_izoline_FS_0_gamreml),AIC(M_izoline_FS_0_gamreml))
  L4<-make_sheet(summary(M_izoline_FS_sex_gamreml),AIC(M_izoline_FS_sex_gamreml))
  L5<-make_sheet(summary(M_izoline_FS_treatment_gamreml),AIC(M_izoline_FS_treatment_gamreml))
  L6<-make_sheet(summary(M_izoline_RE_0_gamreml),AIC(M_izoline_RE_0_gamreml))
  L7<-make_sheet(summary(M_izoline_RE_age_gamreml),AIC(M_izoline_RE_age_gamreml))
  L8<-make_sheet(summary(M_izoline_RE_sex_gamreml),AIC(M_izoline_RE_sex_gamreml))
  L9<-make_sheet(summary(M_izoline_RE_treatment_gamreml),AIC(M_izoline_RE_treatment_gamreml))
  L10<-make_sheet(summary(M_izoline_RE_sex_treatment_gamreml),AIC(M_izoline_RE_sex_treatment_gamreml))
  L11<-make_sheet(summary(M_izoline_RE_sex_treatment_inter_gamreml),AIC(M_izoline_RE_sex_treatment_inter_gamreml))
  
  xlsx::write.xlsx(L1,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='(a) Most parsimonious model', col.names = FALSE, append = FALSE, row.names = FALSE)
  Sys.sleep(1)
  xlsx::write.xlsx(Model_Selection,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='(b) Model selection', col.names = TRUE, append = TRUE, row.names = TRUE)
  Sys.sleep(1)
  xlsx::write.xlsx(L2,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='Losing non-random effect model', col.names = FALSE, append = TRUE, row.names = FALSE)
  Sys.sleep(1)
  xlsx::write.xlsx(L3,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='Losing FS model 1', col.names = FALSE, append = TRUE, row.names = FALSE)
  Sys.sleep(1)
  xlsx::write.xlsx(L4,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='Losing FS model 2', col.names = FALSE, append = TRUE, row.names = FALSE)
  Sys.sleep(1)
  xlsx::write.xlsx(L5,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='Losing FS model 3', col.names = FALSE, append = TRUE, row.names = FALSE)
  Sys.sleep(1)
  xlsx::write.xlsx(L6,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='Losing RE model 1', col.names = FALSE, append = TRUE, row.names = FALSE)
  Sys.sleep(1)
  xlsx::write.xlsx(L7,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='Losing RE model 2', col.names = FALSE, append = TRUE, row.names = FALSE)
  Sys.sleep(1)
  xlsx::write.xlsx(L8,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='Losing RE model 3', col.names = FALSE, append = TRUE, row.names = FALSE)
  Sys.sleep(1)
  xlsx::write.xlsx(L9,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='Losing RE model 4', col.names = FALSE, append = TRUE, row.names = FALSE)
  Sys.sleep(1)
  xlsx::write.xlsx(L10,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='Losing RE model 5', col.names = FALSE, append = TRUE, row.names = FALSE)
  Sys.sleep(1)
  xlsx::write.xlsx(L11,file=paste('./tables/TAB2_GAM_summary_internal',mod,'.xlsx',sep=''),sheetName='Losing RE model 6', col.names = FALSE, append = TRUE, row.names = FALSE)
}

################################################################################
# Model predictions for each izoline; mortality and survivorship
################################################################################

PI1<-lapply(ALLIZO,function(k) predict(m0,data.frame(midday=X, sex='F', treatment='C', izoline=k, sextr='F_C'),type='response',se.fit = TRUE))
LI1<-lapply(ALLIZO,function(k) predict(m0,data.frame(midday=X, sex='F', treatment='C', izoline=k, sextr='F_C'),type='link',se.fit = TRUE))
SI1<-lapply(seq_along(PI1), function(k) exp(-cumsum(PI1[[k]]$fit)))
PI2<-lapply(ALLIZO,function(k) predict(m0,data.frame(midday=X, sex='F', treatment='R', izoline=k, sextr='F_R'),type='response',se.fit = TRUE))
LI2<-lapply(ALLIZO,function(k) predict(m0,data.frame(midday=X, sex='F', treatment='R', izoline=k, sextr='F_R'),type='link',se.fit = TRUE))
SI2<-lapply(seq_along(PI2), function(k) exp(-cumsum(PI2[[k]]$fit)))
PI3<-lapply(ALLIZO,function(k) predict(m0,data.frame(midday=X, sex='M', treatment='C', izoline=k, sextr='M_C'),type='response',se.fit = TRUE))
LI3<-lapply(ALLIZO,function(k) predict(m0,data.frame(midday=X, sex='M', treatment='C', izoline=k, sextr='M_C'),type='link',se.fit = TRUE))
SI3<-lapply(seq_along(PI3), function(k) exp(-cumsum(PI3[[k]]$fit)))
PI4<-lapply(ALLIZO,function(k) predict(m0,data.frame(midday=X, sex='M', treatment='R', izoline=k, sextr='M_R'),type='response',se.fit = TRUE))
LI4<-lapply(ALLIZO,function(k) predict(m0,data.frame(midday=X, sex='M', treatment='R', izoline=k, sextr='M_R'),type='link',se.fit = TRUE))
SI4<-lapply(seq_along(PI4), function(k) exp(-cumsum(PI4[[k]]$fit)))

################################################################################
# bootstrapped marginals and their differences
################################################################################

marginals.basic<-function(mu, N){
  dx<-0
  Sx<-0
  for(k in seq_along(mu)){ #k = izolinia
    sx<-c(1,exp(-cumsum(mu[[k]]$fit)))
    sx<-sx[1:(length(sx)-1)]
    dx<-dx+N[k]*mu[[k]]$fit*sx
    Sx<-Sx+N[k]*sx
  }  
  est<-dx/Sx
  est
}

boot_marginals<-function(logmu, N, B=100000) {
  w<-t(sapply(1:B,function(y){
    newmu<-lapply(logmu, function(k) list(fit=exp(rnorm(length(k$fit),k$fit,k$se))))
    marginals.basic(newmu,N)
  }))
  t(apply(w,2,quantile, probs=c(0.025,0.5,0.975)))
}

boot_marginals_diff<-function(logmu1, N1, logmu2, N2, B=1e5) {
  w<-t(sapply(1:B,function(y){
    newmu1<-lapply(logmu1, function(k) list(fit=exp(rnorm(length(k$fit),k$fit,k$se))))
    newmu2<-lapply(logmu2, function(k) list(fit=exp(rnorm(length(k$fit),k$fit,k$se))))
    marginals.basic(newmu1,N1)-marginals.basic(newmu2,N2)
  }))
  t(apply(w,2,quantile, probs=c(0.025,0.5,0.975)))
}

if (file.exists('Marginals_Boot.rda')) {
  load('Marginals_Boot.rda')
} else {
  mPI1B<-boot_marginals(LI1,NI1) # FC
  mPI2B<-boot_marginals(LI2,NI2) # FR
  mPI3B<-boot_marginals(LI3,NI3) # MC
  mPI4B<-boot_marginals(LI4,NI4) # MR
  mdP12B<-boot_marginals_diff(LI1,NI1,LI2,NI2)
  mdP34B<-boot_marginals_diff(LI3,NI3,LI4,NI4)
  mdP13B<-boot_marginals_diff(LI1,NI1,LI3,NI3)
  mdP24B<-boot_marginals_diff(LI2,NI2,LI4,NI4)
  save(list=c('mPI1B','mPI2B','mPI3B','mPI4B','mdP12B','mdP24B','mdP13B','mdP34B'),file='Marginals_Boot.rda')
}

sig_range_boot<-function(x, lof, hif) {
  rng<-x[(hif<0) | (lof>0)]
  rnglo<-x[(hif<0)]
  rnghi<-x[(lof>0)]
  testr<-c(1,which(c(diff(rng),1)!=1),length(rng))
  sapply(seq_along(testr)[-length(testr)], function(j) {
    g<-range(rng[c((testr[j]+(j!=1)*1):testr[j+1])])
    col<-c(c_red,c_blue)[1+all(g%in%rnglo)]
    lines(g,c(0,0),col=col,lwd=3)
    abline(v=g,col=col,lty=3)
    g
  })
}

#tiff(filename=paste('SFIG3_mortality_differences_boot',mod,'.tiff',sep=''),width=res*8,height=res*7,compression ='lzw',res=res,units='px')
pdf(paste('./figures/SFIG3_mortality_differences_boot',mod,'.pdf',sep=''),width=8,height=7)
par(mfrow=c(2,2))
par(mar=c(4,4,0,0),oma=c(1,1,0.15,0.15))

ind<-X<=min(max(BI1),max(BI2))
plot(X[ind],mdP12B[,3][ind],type='l',lty=2,xlab='',ylab = 'Marginal mortality differences',ylim=c(-0.05,0.05),xlim=c(0,maxX))
lines(X[ind],mdP12B[,1][ind],type='l',lty=2)
lines(X[ind],mdP12B[,2][ind],type='l',lty=1)
legend('bottomleft',bty='n',legend='Control - Rapamycin (Females)',cex=1,inset=c(-0.05,0))
legend('topleft',expression(bold(a)),bty='n',inset=c(-0.075,0),cex=1.3)
abline(h=0,col='gray')
sig_range_boot(X[ind],mdP12B[,1][ind],mdP12B[,3][ind]) #0-17

ind<-X<=min(max(BI3),max(BI4))
plot(X[ind],mdP34B[,3][ind],type='l',lty=2,xlab='',ylab = '',ylim=c(-0.05,0.05),xlim=c(0,maxX))
lines(X[ind],mdP34B[,1][ind],type='l',lty=2)
lines(X[ind],mdP34B[,2][ind],type='l',lty=1)
abline(h=0,col='gray')
sig_range_boot(X[ind],mdP34B[,1][ind],mdP34B[,3][ind]) #0-27, 69-83
legend('bottomleft',bty='n',legend='Control - Rapamycin (Males)',cex=1,inset=c(-0.05,0))
legend('topleft',expression(bold(b)),bty='n',inset=c(-0.075,0),cex=1.3)

ind<-X<=min(max(BI1),max(BI3))
plot(X[ind],mdP13B[,3][ind],type='l',lty=2,xlab='Age',ylab = 'Marginal mortality differences',ylim=c(-0.05,0.05),xlim=c(0,maxX))
lines(X[ind],mdP13B[,1][ind],type='l',lty=2)
lines(X[ind],mdP13B[,2][ind],type='l',lty=1)
abline(h=0,col='gray')
sig_range_boot(X[ind],mdP13B[,1][ind],mdP13B[,3][ind])
legend('bottomleft',bty='n',legend='Females - Males (Control)',cex=1,inset=c(-0.05,0))
legend('topleft',expression(bold(c)),bty='n',inset=c(-0.075,0),cex=1.3)

ind<-X<=min(max(BI2),max(BI4))
plot(X[ind],mdP24B[,3][ind],type='l',lty=2,xlab='Age',ylab = '',ylim=c(-0.05,0.05),xlim=c(0,maxX))
lines(X[ind],mdP24B[,1][ind],type='l',lty=2)
lines(X[ind],mdP24B[,2][ind],type='l',lty=1)
abline(h=0,col='gray')
sig_range_boot(X[ind],mdP24B[,1][ind],mdP24B[,3][ind])
legend('bottomleft',bty='n',legend='Females - Males (Rapamycin)',cex=1,inset=c(-0.05,0))
legend('topleft',expression(bold(d)),bty='n',inset=c(-0.075,0),cex=1.3)

dev.off()

################
# Mortality
################

#tiff(filename=paste('FIG2_Mortality_boot',mod,'.tiff',sep=''),width=res*8,height=res*7,compression ='lzw',res=res,units='px')
pdf(paste('./figures/FIG2_Mortality_boot',mod,'.pdf',sep=''),width=8,height=7)
par(mfrow=c(2,2))
par(mar=c(4,4,0,0),oma=c(1,1,0.15,0.15))

colCR=adjustcolor(c(c_blue,c_red),alpha.f = 0.5)
colFM=adjustcolor(c(c_green,c_yellow),alpha.f = 0.5)
Ylim=c(-4,0)

plot(X[X<=max(BI1)],log10(mPI1B[,2])[X<=max(BI1)],type='l',col=c_blue,lwd=3,ylim=Ylim,axes=F, ylab='Mortality rate', xlab='',xlim=c(0,maxX))
axis(1);magicaxis::magaxis(2,unlog = TRUE,las=2)
for(k in seq_along(PI1)) lines(X[X<=BI1[k]],log10(PI1[[k]]$fit[X<=BI1[k]]),col=colCR[1])
for(k in seq_along(PI2)) lines(X[X<=BI2[k]],log10(PI2[[k]]$fit[X<=BI2[k]]),col=colCR[2])
lines(X[X<=max(BI2)],log10(mPI2B[,2])[X<=max(BI2)],type='l',col=c_red,lwd=3)
lines(X[X<=max(BI1)],log10(mPI1B[,2])[X<=max(BI1)],type='l',col=c_blue,lwd=3)
legend('bottomleft','Females',bty='n')
legend('bottomright',c('Control','Rapamycin'),col=c(c_blue,c_red),lty=1,bty='n',lwd=3)
legend('topleft',expression(bold(a)),bty='n',inset=c(-0.075,0),cex=1.3)
box(); box()

plot(X[X<=max(BI3)],log10(mPI3B[,2])[X<=max(BI3)],type='l',col=c_blue,lwd=3,ylim=Ylim,axes=F, ylab='', xlab='',xlim=c(0,maxX))
axis(1);magicaxis::magaxis(2,unlog = TRUE,las=2)
for(k in seq_along(PI3)) lines(X[X<=BI3[k]],log10(PI3[[k]]$fit[X<=BI3[k]]),col=colCR[1])
for(k in seq_along(PI4)) lines(X[X<=BI4[k]],log10(PI4[[k]]$fit[X<=BI4[k]]),col=colCR[2])
lines(X[X<=max(BI4)],log10(mPI4B[,2])[X<=max(BI4)],type='l',col=c_red,lwd=3)
lines(X[X<=max(BI3)],log10(mPI3B[,2])[X<=max(BI3)],type='l',col=c_blue,lwd=3)
legend('bottomleft','Males',bty='n')
legend('bottomright',c('Control','Rapamycin'),col=c(c_blue,c_red),lty=1,bty='n',lwd=3)
legend('topleft',expression(bold(b)),bty='n',inset=c(-0.075,0),cex=1.3)
box(); box()

plot(X[X<=max(BI1)],log10(mPI1B[,2])[X<=max(BI1)],type='l',col=c_green,lwd=3,ylim=Ylim,axes=F, ylab='Mortality rate', xlab='Age',xlim=c(0,maxX))
axis(1);magicaxis::magaxis(2,unlog = TRUE,las=2)
for(k in seq_along(PI1)) lines(X[X<=BI1[k]],log10(PI1[[k]]$fit[X<=BI1[k]]),col=colFM[1])
for(k in seq_along(PI3)) lines(X[X<=BI3[k]],log10(PI3[[k]]$fit[X<=BI3[k]]),col=colFM[2])
lines(X[X<=max(BI3)],log10(mPI3B[,2])[X<=max(BI3)],type='l',col=c_yellow,lwd=3)
lines(X[X<=max(BI1)],log10(mPI1B[,2])[X<=max(BI1)],type='l',col=c_green,lwd=3)
legend('bottomleft','Control',bty='n')
legend('bottomright',c('Females','Males'),col=c(c_green,c_yellow),lty=1,bty='n')
legend('topleft',expression(bold(c)),bty='n',inset=c(-0.075,0),cex=1.3)
box(); box()

plot(X[X<=max(BI2)],log10(mPI2B[,2])[X<=max(BI2)],type='l',col=c_green,lwd=3,ylim=Ylim,axes=F, ylab='', xlab='Age',xlim=c(0,maxX))
axis(1);magicaxis::magaxis(2,unlog = TRUE,las=2)
for(k in seq_along(PI2)) lines(X[X<=BI2[k]],log10(PI2[[k]]$fit[X<=BI2[k]]),col=colFM[1])
for(k in seq_along(PI4)) lines(X[X<=BI4[k]],log10(PI4[[k]]$fit[X<=BI4[k]]),col=colFM[2])
lines(X[X<=max(BI4)],log10(mPI4B[,2])[X<=max(BI4)],type='l',col=c_yellow,lwd=3)
lines(X[X<=max(BI2)],log10(mPI2B[,2])[X<=max(BI2)],type='l',col=c_green,lwd=3)
legend('bottomleft','Rapamycin',bty='n')
legend('bottomright',c('Females','Males'),col=c(c_green,c_yellow),lty=1,bty='n')
legend('topleft',expression(bold(d)),bty='n',inset=c(-0.075,0),cex=1.3)
box(); box()

dev.off()


################################################################################
# Model predictions ignoring izoline; mortality and survivorship - New 2022
################################################################################

ELI1<-predict.gam(m0,data.frame(midday=X, sex='F', treatment='C', izoline=ALLIZO[3]),type='link',se.fit = TRUE,
            exclude=c('s(midday,izoline):sexF','s(midday,izoline):sexM','s(midday,izoline):treatmentC','s(midday,izoline):treatmentR'))
ELI2<-predict.gam(m0,data.frame(midday=X, sex='F', treatment='R', izoline=ALLIZO[3]),type='link',se.fit = TRUE,
                  exclude=c('s(midday,izoline):sexF','s(midday,izoline):sexM','s(midday,izoline):treatmentC','s(midday,izoline):treatmentR'))
ELI3<-predict.gam(m0,data.frame(midday=X, sex='M', treatment='C', izoline=ALLIZO[3]),type='link',se.fit = TRUE,
                  exclude=c('s(midday,izoline):sexF','s(midday,izoline):sexM','s(midday,izoline):treatmentC','s(midday,izoline):treatmentR'))
ELI4<-predict.gam(m0,data.frame(midday=X, sex='M', treatment='R', izoline=ALLIZO[3]),type='link',se.fit = TRUE,
                  exclude=c('s(midday,izoline):sexF','s(midday,izoline):sexM','s(midday,izoline):treatmentC','s(midday,izoline):treatmentR'))

boot_excluded_diff<-function(logmu1, logmu2, B=1e5) {
  w<-t(sapply(1:B,function(y){
    logmu1b<-rnorm(length(logmu1$fit),logmu1$fit,logmu1$se)
    logmu2b<-rnorm(length(logmu2$fit),logmu2$fit,logmu2$se)
    (logmu1b)-(logmu2b)
  }))
  t(apply(w,2,quantile, probs=c(0.025,0.5,0.975)))
}

if (file.exists('Excluded_Boot.rda')) {
  load('Excluded_Boot.rda')
} else {
  EmdP12B<-boot_excluded_diff(ELI1,ELI2)
  EmdP34B<-boot_excluded_diff(ELI3,ELI4)
  EmdP13B<-boot_excluded_diff(ELI1,ELI3)
  EmdP24B<-boot_excluded_diff(ELI2,ELI4)
  save(list=c('EmdP12B','EmdP24B','EmdP13B','EmdP34B'),file='Excluded_Boot.rda')
}

# # get difference estimates:
# diff12 <- itsadug::get_difference(m0, comp=list(treatment=c('C', 'R'),sex=c('F','F')), 
#                        cond=list(midday=X), sim.ci = TRUE)
# 
# diff34 <- itsadug::get_difference(m0, comp=list(treatment=c('C', 'R'),sex=c('M','M')), 
#                                   cond=list(midday=X),sim.ci = TRUE)
# 
# diff13 <- itsadug::get_difference(m0, comp=list(treatment=c('C', 'C'),sex=c('F','M')), 
#                                   cond=list(midday=X),sim.ci = TRUE)
# 
# diff24 <- itsadug::get_difference(m0, comp=list(treatment=c('R', 'R'),sex=c('F','M')), 
#                                   cond=list(midday=X),sim.ci = TRUE)
# diff12[,4]<-diff12[,5]
# diff34[,4]<-diff34[,5]
# diff13[,4]<-diff13[,5]
# diff24[,4]<-diff24[,5]

#tiff(filename=paste('FIG3_logmortality_differences_excluded_rand_boot',mod,'.tiff',sep=''),width=res*8,height=res*7,compression ='lzw',res=res,units='px')
pdf(paste('./figures/FIG3_logmortality_differences_excluded_rand_boot',mod,'.pdf',sep=''),width=8,height=7)
par(mfrow=c(2,2))
par(mar=c(4,4,0,0),oma=c(1,1,0.15,0.15))

ind<- X<=min(max(BI1),max(BI2))
plot(X[ind],EmdP12B[,3][ind],type='l',lty=2,xlab='',ylab = 'log mortality differences',ylim=c(-2,2),xlim=c(0,maxX))
lines(X[ind],EmdP12B[,1][ind],type='l',lty=2)
lines(X[ind],EmdP12B[,2][ind],type='l',lty=1)
#lines(X[ind],diff12[,3][ind],type='l',lty=1,col=2)
#lines(X[ind],diff12[,4][ind],type='l',lty=2,col=2)
#lines(X[ind],2*diff12[,3][ind]-diff12[,4][ind],type='l',lty=2,col=2)

legend('bottomleft',bty='n',legend='Control - Rapamycin (Females)',cex=1,inset=c(-0.05,0))
legend('topleft',expression(bold(a)),bty='n',inset=c(-0.075,0),cex=1.3)
abline(h=0,col='gray')
sig_range_boot(X[ind],EmdP12B[,1][ind],EmdP12B[,3][ind]) 
#sig_range_boot(X[ind],2*diff12[,3][ind]-diff12[,4][ind],diff12[,4][ind]) 

ind<-X<=min(max(BI3),max(BI4))
plot(X[ind],EmdP34B[,3][ind],type='l',lty=2,xlab='',ylab = '',ylim=c(-2,2),xlim=c(0,maxX))
lines(X[ind],EmdP34B[,1][ind],type='l',lty=2)
lines(X[ind],EmdP34B[,2][ind],type='l',lty=1)
#lines(X[ind],diff34[,3][ind],type='l',lty=1,col=2)
#lines(X[ind],diff34[,4][ind],type='l',lty=2,col=2)
#lines(X[ind],2*diff34[,3][ind]-diff34[,4][ind],type='l',lty=2,col=2)

abline(h=0,col='gray')
sig_range_boot(X[ind],EmdP34B[,1][ind],EmdP34B[,3][ind]) 
legend('bottomleft',bty='n',legend='Control - Rapamycin (Males)',cex=1,inset=c(-0.05,0))
legend('topleft',expression(bold(b)),bty='n',inset=c(-0.075,0),cex=1.3)

ind<-X<=min(max(BI1),max(BI3))
plot(X[ind],EmdP13B[,3][ind],type='l',lty=2,xlab='Age',ylab = 'log mortality differences',ylim=c(-2,2),xlim=c(0,maxX))
lines(X[ind],EmdP13B[,1][ind],type='l',lty=2)
lines(X[ind],EmdP13B[,2][ind],type='l',lty=1)
#lines(X[ind],diff13[,3][ind],type='l',lty=1,col=2)
#lines(X[ind],diff13[,4][ind],type='l',lty=2,col=2)
#lines(X[ind],2*diff13[,3][ind]-diff13[,4][ind],type='l',lty=2,col=2)

abline(h=0,col='gray')
sig_range_boot(X[ind],EmdP13B[,1][ind],EmdP13B[,3][ind])
legend('bottomleft',bty='n',legend='Females - Males (Control)',cex=1,inset=c(-0.05,0))
legend('topleft',expression(bold(c)),bty='n',inset=c(-0.075,0),cex=1.3)

ind<-X<=min(max(BI2),max(BI4))
plot(X[ind],EmdP24B[,3][ind],type='l',lty=2,xlab='Age',ylab = '',ylim=c(-2,2),xlim=c(0,maxX))
lines(X[ind],EmdP24B[,1][ind],type='l',lty=2)
lines(X[ind],EmdP24B[,2][ind],type='l',lty=1)
#lines(X[ind],diff24[,3][ind],type='l',lty=1,col=2)
#lines(X[ind],diff24[,4][ind],type='l',lty=2,col=2)
#lines(X[ind],2*diff24[,3][ind]-diff24[,4][ind],type='l',lty=2,col=2)

abline(h=0,col='gray')
sig_range_boot(X[ind],EmdP24B[,1][ind],EmdP24B[,3][ind])
legend('bottomleft',bty='n',legend='Females - Males (Rapamycin)',cex=1,inset=c(-0.05,0))
legend('topleft',expression(bold(d)),bty='n',inset=c(-0.075,0),cex=1.3)

dev.off()


################################################################################
# Goodness of fit
################################################################################

#tiff(filename=paste('FIG2S_Goodness_of_fit',mod,'.tiff',sep=''),width=res*8,height=res*9,compression ='lzw',res=res,units='px')
pdf(paste('./figures/SFIG1_Goodness_of_fit',mod,'.pdf',sep=''),width=8,height=9)
par(mar=c(4,4,1,1))
csf <- npsurv(Surv(midday,status) ~ treatment+sex, data = RapamycinData$survdat)
par(mfrow=c(2,1))
plot(csf,col=c(1,c_blue,c_red,c_yellow),mark.time = !TRUE, conf.int = FALSE, xlab='Age', ylab='Survivorship',xlim=c(skip,maxX),lwd=1)
lines(c(X,NA), c(1,exp(-cumsum(mPI1B[,2]))),col=1,lwd=2)
lines(c(X,NA), c(1,exp(-cumsum(mPI3B[,2]))),col=c_blue,lwd=2)
lines(c(X,NA), c(1,exp(-cumsum(mPI2B[,2]))),col=c_red,lwd=2)
lines(c(X,NA), c(1,exp(-cumsum(mPI4B[,2]))),col=c_yellow,lwd=2)
legend('topright',bty='n',legend='Kaplan-Meier vs. model-estimated marginal survivroships')
legend('topleft',expression(bold(a)),bty='n',inset=c(-0.04,0),cex=1.3)
legend('bottomleft',legend=gsub(')','',fixed=TRUE,gsub('(','',fixed=TRUE,convert_terms_3(names(csf$strata)))),col=c(1,c_blue,c_red,c_yellow),lwd=2,bty='n')
acf_resid(m0,main='')
legend('topright',bty='n',legend='Autocorrelation test')
legend('topleft',expression(bold(b)),bty='n',inset=c(-0.04,0),cex=1.3)
dev.off()


############## plotted "manually" please check original smooth ylab if consistent with plot ylab
Sm<-plot(m0,select = 0)

#tiff(filename=paste('SFIG1A_Model_Smooths_effects_on_log_mortality',mod,'.tiff',sep=''),width=res*8,height=res*9,compression ='lzw',res=res,units='px')
pdf(paste('./figures/SFIG2A_Model_Smooths_effects_on_log_mortality',mod,'.pdf',sep=''),width=8,height=9)
#pdf(file='FIG1S_Model_Smooths_effects_on_log_mortality.pdf',8,9,onefile = TRUE)
par(mfrow=c(2,2))
par(mar=c(4,4,1,1.5),oma=c(0,0,0,0))

print(Sm[[1]]$ylab)
plot(Sm[[1]]$x,Sm[[1]]$fit,type='l',xlab='Age',ylab="Effect of s(Age) : Control : Females",xlim=c(0,maxX),ylim=c(-3,6))
lines(Sm[[1]]$x,Sm[[1]]$fit+1.96*Sm[[1]]$se,type='l',lty=2)
lines(Sm[[1]]$x,Sm[[1]]$fit-1.96*Sm[[1]]$se,type='l',lty=2)
legend('topright','Fixed effects on log mortality rate:\nModel smooths by age for females in control',bty='n',cex=0.9)
legend('topleft',expression(bold(a)),bty='n',inset=c(-0.04,0),cex=1.3)

print(Sm[[2]]$ylab)
plot(Sm[[2]]$x,Sm[[2]]$fit,type='l',xlab='Age',ylab="Effect of s(Age) : Control : Males",xlim=c(0,maxX),ylim=c(-3,6))
lines(Sm[[2]]$x,Sm[[2]]$fit+1.96*Sm[[1]]$se,type='l',lty=2)
lines(Sm[[2]]$x,Sm[[2]]$fit-1.96*Sm[[1]]$se,type='l',lty=2)
legend('topright','Fixed effects on log mortality rate:\nModel smooths by age for males in control',bty='n',cex=0.9)
legend('topleft',expression(bold(b)),bty='n',inset=c(-0.04,0),cex=1.3)

print(Sm[[3]]$ylab)
plot(Sm[[3]]$x,Sm[[3]]$fit,type='l',xlab='Age',ylab="Effect of s(Age) : Rapamycin : Females",xlim=c(0,maxX),ylim=c(-3,6))
lines(Sm[[3]]$x,Sm[[3]]$fit+1.96*Sm[[1]]$se,type='l',lty=2)
lines(Sm[[3]]$x,Sm[[3]]$fit-1.96*Sm[[1]]$se,type='l',lty=2)
legend('topright','Fixed effects on log mortality rate:\nModel smooths by age for females in rapamycin',bty='n',cex=0.9)
legend('topleft',expression(bold(c)),bty='n',inset=c(-0.04,0),cex=1.3)

print(Sm[[4]]$ylab)
plot(Sm[[4]]$x,Sm[[4]]$fit,type='l',xlab='Age',ylab="Effect of s(Age) : Rapamycin : Males",xlim=c(0,maxX),ylim=c(-3,6))
lines(Sm[[4]]$x,Sm[[4]]$fit+1.96*Sm[[1]]$se,type='l',lty=2)
lines(Sm[[4]]$x,Sm[[4]]$fit-1.96*Sm[[1]]$se,type='l',lty=2)
legend('topright','Fixed effects on log mortality rate:\nModel smooths by age for males in rapamycin',bty='n',cex=0.9)
legend('topleft',expression(bold(d)),bty='n',inset=c(-0.04,0),cex=1.3)
dev.off()

#tiff(filename=paste('SFIG2B_Model_Smooths_effects_on_log_mortality',mod,'.tiff',sep=''),width=res*8,height=res*9,compression ='lzw',res=res,units='px')
pdf(paste('./figures/SFIG2B_Model_Smooths_effects_on_log_mortality',mod,'.pdf',sep=''),width=8,height=9)
par(mfrow=c(2,2))
par(mar=c(4,4,1,1.5),oma=c(0,0,0,0))

print(Sm[[5]]$ylab)
L<-length(Sm[[5]]$x)
D<-length(Sm[[5]]$fit)/L-1
COL<-terrain.colors(D+1)
plot(Sm[[5]]$x,Sm[[5]]$fit[1:L],type='l',xlab='Age',ylab="Effect of s(Age, Iso) : Females",xlim=c(0,maxX),ylim=c(-3,6),col=1)#COL[1])
for (k in 0:D) lines(Sm[[5]]$x,Sm[[5]]$fit[(1:L)+k*L],col=1)#COL[k+2])
legend('topright','Random effects on log mortality rate:\nModel smooths by isoline for females',bty='n',cex=0.9)
legend('topleft',expression(bold(e)),bty='n',inset=c(-0.04,0),cex=1.3)

print(Sm[[6]]$ylab)
L<-length(Sm[[6]]$x)
D<-length(Sm[[6]]$fit)/L-1
COL<-terrain.colors(D+1)
plot(Sm[[6]]$x,Sm[[6]]$fit[1:L],type='l',xlab='Age',ylab="Effect of s(Age, Iso) : Males",xlim=c(0,maxX),ylim=c(-3,6),col=1)#COL[1])
for (k in 0:D) lines(Sm[[6]]$x,Sm[[6]]$fit[(1:L)+k*L],col=1)#COL[k+2])
legend('topright','Random effects on log mortality rate:\nModel smooths by isoline for males',bty='n',cex=0.9)
legend('topleft',expression(bold(f)),bty='n',inset=c(-0.04,0),cex=1.3)

print(Sm[[7]]$ylab)
L<-length(Sm[[7]]$x)
D<-length(Sm[[7]]$fit)/L-1
COL<-terrain.colors(D+1)
plot(Sm[[7]]$x,Sm[[7]]$fit[1:L],type='l',xlab='Age',ylab="Effect of s(Age, Iso) : Control",xlim=c(0,maxX),ylim=c(-3,6),col=1)#COL[1])
for (k in 0:D) lines(Sm[[7]]$x,Sm[[7]]$fit[(1:L)+k*L],col=1)#COL[k+2])
legend('topright','Random effects on log mortality rate:\nModel smooths by isoline control',bty='n',cex=0.9)
legend('topleft',expression(bold(g)),bty='n',inset=c(-0.04,0),cex=1.3)

print(Sm[[8]]$ylab)
L<-length(Sm[[8]]$x)
D<-length(Sm[[8]]$fit)/L-1
COL<-terrain.colors(D+1)
plot(Sm[[8]]$x,Sm[[8]]$fit[1:L],type='l',xlab='Age',ylab="Effect of s(Age, Iso) : Rapamycin",xlim=c(0,maxX),ylim=c(-3,6),col=1)#COL[1])
for (k in 0:D) lines(Sm[[8]]$x,Sm[[8]]$fit[(1:L)+k*L],col=1)#COL[k+2])
legend('topright','Random effects on log mortality rate:\nModel smooths by isoline for rapamycin',bty='n',cex=0.9)
legend('topleft',expression(bold(h)),bty='n',inset=c(-0.04,0),cex=1.3)

dev.off()


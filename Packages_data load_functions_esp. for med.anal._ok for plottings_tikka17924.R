
# Simplifying the analysis by making the loading of libraries, datas and functions as straightforward as possible, i.e.
# select all and press enter :)

# Because some of these you may need for your analysis, Pauli Tikka: 1.3.24


#Select libraries and data loads all at once:
#Library loads. If problem to load: delete the package and reinstall, if that does not work copy your (or from internet) earlier version of the package.
library(scales);library("plotrix");library(Maaslin2);  library(ggtext);library(lavaan);library(psych);library("xlsx");library(lsr);library(quantreg);library("readxl")
library(semPlot);library(mediation);require(lme4);require(reshape2);library(mlma);library(binilib);library(plyr);library("viridis");library(RColorBrewer) #library(mma); # library(d3heatmap);
library(magrittr); library(ggplot2);library("censReg" );library(ggsankey);library(dplyr);library(tidyverse);library(dmetar);library(meta);library(ggforestplot); 
## note! do not add: library(forestplot!!!not working with ggforestplot..);
library(mdatools);library(circlize);library(igraph);library('bigsnpr');library(rcompanion);library(scRNAseq);library(tibble);library(stringr);library(MOFA2);library('qpgraph') #ok
library("grid"); library("ggplotify");library(ggpubr);library(rstatix);library(datarium);library(RColorBrewer); library(ggh4x); library(effsize)
library(chorddiag);library(corrplot);library(scater);library(mdatools);options(scipen = 999); library(car);library(FSA);library(pathviewr)#jotkin vanhoista paketeista (esim. SparkL) sattaaa hidastaa
library("lmtest");library(PerformanceAnalytics);library(psych);library("readxl");library(ggforce);library(ComplexHeatmap) #these are ok to drive in start
library('Hmisc');library(correlation);library(ggcorrplot);library(pheatmap);library(mgcv);library(lpSolve);library(glmnet); library(pathviewr)
library(extrafont)
font_import() #this is important
loadfonts(device = "win") #this is important too
# Data Loads... (if you have these you may skip these as now): 
setwd("C:/Users/patati/Desktop/TurkuOW/RWork")
date='tikka78524' #eli tässä ne on...
set.seed(30)

# Steroid data, simple figures:# NAFLD=read.csv("NAFLD_SteroidStudy_v3.csv", header = TRUE, sep=";")
NAFLD=read_excel("NAFLD_SteroidStudy.xlsx",sheet = "LFAT_steroidsDATA") #l ei tästä
oknames=colnames(NAFLD); NAFLD=data.frame(NAFLD)
# rownames(NAFLD)=NAFLD[,1]
groups=read.csv("groups_17823.csv", header = TRUE, sep=";")
#P4 was found from elswhere to be:
NAFLD[,'P4'] = as.numeric(NAFLD[,'P4'])
NAFLD[,'P4'][is.na(NAFLD[,'P4'])] = 22557.3330346846#median(NAFLD[,'P4'], na.rm=TRUE) 
NAFLD[,5:7][NAFLD[,5:7]==0.01]=0; colnames(NAFLD)=oknames
MASLD=read_excel("Combined.Matrix.For.Pauli.2023.10.17.Excel.Formatv2.xlsx") #tästä!
oknames=colnames(MASLD); MASLD=data.frame(MASLD)
colnames(MASLD)=oknames
rownames(MASLD)=MASLD[,1]
#P4 was found from elswhere to be:
MASLD[,'P4'] = as.numeric(MASLD[,'P4'])
MASLD[,'P4'][is.na(MASLD[,'P4'])] = 22557.3330346846 
eva=c('Grade(0-3)', 'Stage(0-4)','Necroinflammation')
MASLD[,eva][MASLD[,eva]==0.01]=0; #MASLD[,eva]
td=c('11-KDHT','AN','DHT','17a-OHP5','E2','P5','DOC')
val=c(103,252,51,200,26.5,253,10); vale=c(100,250,50,200,25,250,10)#MASLD[1:30,td] 
for (i in 1:7) {MASLD[,td][i][MASLD[,td][i]==val[i]]=vale[i]} #MASLD[1:30,td] # tu=c('E','11-KA4') # val=c(106000) # vale=c(100) # MASLD[,tu]

# These (E) are ok as per lab:
ME=read.csv('E_tikka231023.csv',header=TRUE, sep=";")
ME2=rownames(MASLD[MASLD[,'E']==106000,]) 
to=ME[which(ME[,1] %in% ME2),'patient.number']
te=ME[which(ME[,1] %in% ME2),'E']
MASLD[as.character(to),'E']=te
#These (11-KA4) will change in the lab (sometime after 24.10.23):
M11=read.csv('11KA4_tikka231023.csv',header=TRUE, sep=";")
M11[,1][c(1:5,9)];MASLD[as.character(M11[,1][c(1:5,9)]),'11-KA4'] #these wer denoted with 'big interference'
MASLD[as.character(M11[,1][c(1:5,9)]),'11-KA4'] = NA#median(MASLD[!rownames(MASLD) %in% as.character(M11[,1][c(1:5,9)]),'11-KA4'])
a=MASLD[order(MASLD[,'BMI']),'BMI']
b=NAFLD[order(NAFLD[,'BMI']),'BMI']
them=unique(b[! b %in% a])
NAFLD=NAFLD[order(NAFLD[,'BMI']),] 
NAFLD=NAFLD[NAFLD[,'BMI']!=them,]
MASLD=MASLD[order(MASLD[,'BMI']),]
#https://appsilon.com/imputation-in-r/ #https://www.datasciencemadesimple.com/get-minimum-value-of-a-column-in-r-2/?expand_article=1
#New data import withouth changing the conames: https://readxl.tidyverse.org/articles/column-names.html
Bali=data.frame(read_excel("Liver_bile_acids_PFAS.xlsx",sheet = "Liver_BA",.name_repair = "minimal")); row.names(Bali)=Bali[,1]
Pfase=data.frame(read_excel("Liver_bile_acids_PFAS.xlsx",sheet = "PFAS_serum",.name_repair = "minimal")); rownames(Pfase)=as.vector(unlist(Pfase[,1]))
Base=data.frame(read_excel("Liver_bile_acids_PFAS.xlsx",sheet = "Serum_BA",.name_repair = "minimal"));rownames(Base)=as.vector(unlist(Base[,1]))
C4=data.frame(read_excel("Liver_bile_acids_PFAS.xlsx",sheet = "C4",.name_repair = "minimal")); rownames(C4)=as.vector(unlist(C4[,1]))
Clini=data.frame(read_excel("Matching clinical data_all.xlsx",sheet = "Sheet1",.name_repair = "minimal")); rownames(Clini)=as.vector(unlist(Clini[,1]));
# https://www.analyticsvidhya.com/blog/2021/06/hypothesis-testing-parametric-and-non-parametric-tests-in-statistics/
head(NAFLD);head(MASLD)
#The below ordering needs to be changed...
Bali=Bali[as.character(MASLD$PatientNumber),];Bali[1:3,] #https://stackoverflow.com/questions/54264980/r-how-to-set-row-names-attribute-as-numeric-from-character I did otherway around
Base=Base[as.character(MASLD$PatientNumber),];Base[1:3,]
Clini=Clini[as.character(MASLD$PatientNumber),];Clini[1:3,]
C4=C4[as.character(MASLD$PatientNumber),];C4[1:3,]
Pfase=Pfase[as.character(MASLD$PatientNumber),];Pfase[1:3,]
#Menopause markers:
menopause=read_excel("Putative_metabolic_markers_menopause.xlsx",sheet='menopause markers',.name_repair = "minimal"); #rownames(Clini)=as.vector(unlist(Clini[,1]));
menopause=menopause[8:dim(menopause)[1],]; menopause=menopause[,-15]; menopause[2,2:14]=menopause[1,2:14]; menopause=data.frame(menopause); menopause[2,13:14]=c('v1','v2'); dim(menopause)
colnames(menopause)=c('row_names',menopause[2,2:dim(menopause)[2]]); menopause=menopause[3:dim(menopause)[1],];rownames(menopause)=as.vector(unlist(menopause[,1]));
menopause=menopause[as.character(MASLD$PatientNumber),]
colnames(Pfase)[colnames(Pfase)=='PFHxA.1']='PFHxA_Branched'
Pfase=Pfase[,colnames(Pfase)!='Benzylparaben.1']
Pfase[Pfase[,'Benzylparaben']>10,'Benzylparaben']=NA #m

Jeihou=data.frame(read_excel("Copy of BA_liverfat_RawData.xls",.name_repair = "minimal")); row.names(Jeihou)=Jeihou[,1];Jeihou=Jeihou[as.character(MASLD$PatientNumber),]
u=Jeihou[Jeihou[,'GHDGA']=='<LLOQ',1]; a=u[!is.na(u)];length(a); b=rownames(Bali[Bali[,'GHDGA']==1,]); length(b)
intersect(a,b); uu=Jeihou[Jeihou[,'GHDGA']=='No Result',1]; aa=uu[!is.na(uu)]; intersect(aa,b); c(aa,a)[!c(aa,a) %in% b] #24140250313 24112081112  #2/25
Bali[as.character(a),'GHDGA']=min(Bali[,'GHDGA'],na.rm=TRUE)/2
heps=Bali[Bali[,'GHDGA']==1,1] #2476250110  2487010610 24111141210
Bali[as.character(heps),'GHDGA']=NA
#https://www.datasciencemadesimple.com/get-minimum-value-of-a-column-in-r-2/?expand_article=1
mat=Bali[,c('TbMCA','ToMCA','TDCA','TDHCA','TLCA')]
mat[!mat>1]=10000
mat[mat==2]=10000 #colmins ei toiminuyt ja käytin:
hip=do.call(pmin, lapply(1:nrow(mat), function(i)mat[i,])) #https://stackoverflow.com/questions/13676878/fastest-way-to-get-min-from-every-column-in-a-matrix
hou=c('TbMCA','ToMCA','TDCA','TDHCA','TLCA')
for (i in 1:5) {Bali[Bali[,hou[i]]==1,hou[i]]=hip[i]}
for (i in 1:5) {Bali[Bali[,hou[i]]==2,hou[i]]=hip[i]}# hepsa=Bali[Bali[,c('TbMCA','ToMCA','TDCA','TDHCA','TLCA')]==1,1] #2476250110  2487010610 24111141210
#An imputation for missing values:
C4[is.na(C4[,2]),2]=median(C4[!is.na(C4[,2]),2]) #assuming that these were not below quantitation and replacing with median
#https://www.geeksforgeeks.org/performing-logarithmic-computations-in-r-programming-log-log10-log1p-and-log2-functions/# https://stackoverflow.com/questions/50476717/i-want-to-align-match-two-unequal-columns
#Matching two unequal columns..# #match the names of one original column (dat2) to ones that are missing (dat1 with to other) # #not sure if this should be this diffucult...
tv=cbind(MASLD[,1],NAFLD[,2:7],Clini[,'HOMA.IR'],MASLD[,colnames(NAFLD[,8:27])],Bali[,2:dim(Bali)[2]], C4[,2:dim(C4)[2]],Base[,2:dim(Base)[2]],Pfase[,(2:(dim(Pfase)[2]))], MASLD[,'PFAS']);colnames(tv)#,C4[,2:dim(C4)[2]]). Clini[,'HOMA-IR']
head(tv) #non nans 
which(is.na(tv))
# MASLD[1:30,1:12] # NAFLD[1:30,7:20]
colnames(tv)[colnames(tv)=='C4[, 2:dim(C4)[2]]']='C4';colnames(tv)[colnames(tv)=='Clini[, \"HOMA.IR\"]']='HOMA-IR'
colnames(tv)[colnames(tv)=='MASLD[, \"PFAS\"]']='PFAS';
colnames(tv)[colnames(tv)=="MASLD[, 1]" ]='PatientNumber';colnames(tv)#
rownames(tv)=unlist(Bali[,1]); tv[1:5,1:11];#tv[1:5,12:55]; dim(tv[1:3,9:28]);tv[1:5,1:80]
hep=colnames(tv)[!colnames(tv) %in% c( "Benzylparaben" ,"Methylparaben")] 
#not sure when it is the best time to take not needed variables away, perhaps at the very end?
tv=tv[,hep]
tv=cbind(tv,MASLD[,(dim(MASLD)[2]-13):dim(MASLD)[2]]) #here I add the lipids. In the future, I need to divide all the groups in their own components e.g. dataframe called 'lipids' so
#that adding them will be more straighfoward
head(tv) #non nans , ok colnames
which(is.na(tv))
# rna=rownames(tv)
# cna=colnames(tv)
# tv <- mapply(tv, FUN=as.numeric)
# tv=data.frame(tv)
# rownames(tv)=rna
# colnames(tv)=cna
# tv[,'MUFA']=abs(min(tv[,'MUFA']))+tv[,'MUFA']+0.0001 #this was for one test with different calculation for the mufa
tve=tv[,2:dim(tv)[2]]; tve[tve == 0] <- NA; #print(tve, max=nrow(tve)*ncol(tve)); note, here the covariates have not been normalized or scaled/elaborated in any way
tv_half <- tve %>% mutate(replace(., is.na(.), min(., na.rm = T)/2)) #https://mdatools.com/docs/preprocessing--autoscaling.html
tv_half_log2 <- log2(tv_half);
# print(tv_half_log2, max=nrow(tv_half_log2)*ncol(tv_half_log2))
tv_auto <- prep.autoscale(tv_half_log2, center = TRUE, scale = TRUE); 
head(tv_auto) #non nans 
which(is.na(tv_auto))

# tv_auto[1:5,1:11] 
#usually this should be the log2 value 'tv_half_log2' & #https://svkucheryavski.gitbooks.io/mdatools/content/preprocessing/text.html
# Necroinflammation  HOMA-IR Steatosis.Grade.0.To.3 Fibrosis.Stage.0.to.4
tv_all=cbind(tv[,1],tv_auto); #tv_all[1:5,1:11]; note, here the covariates have not been normalized or scaled/elaborated in any way
x1=colnames(tv_all[,c(1:8)]); v2=dim(NAFLD)[2]+1
x2=colnames(tv_all[,9:v2]);v3=(dim(Bali)[2]+v2);x3=colnames(tv_all[,(v2+1):(v3)]);v4=(dim(Base)[2])+v3
x4=colnames(tv_all[,(v3+1):(v4-1)]);x5=colnames(tv_all[,(v4):(dim(tv_all)[2])]); 
x3 <- paste(x3, "_L", sep="") #https://stackoverflow.com/questions/6984796/how-to-paste-a-string-on-each-element-of-a-vector-of-strings-using-apply-in-r
x4=gsub("(-[0-9]*)*.1", "", x4) #https://stackoverflow.com/questions/18997297/remove-ending-of-string-with-gsub
x4 <- paste(x4, "_S", sep="")# https://rdrr.io/bioc/qpgraph/man/qpNrr.html
x5a=x5[1:9]
x6=x5[10:length(x5)] #dividing to lipids
x5=x5a  #making sure that PFAS are separate
nm = c(x1,x2,x3,x4,x5,x6); nm=c('PatientNumber','Gender','AGE','BMI','Steatosis Grade','Fibrosis Stage','Necroinflammation','HOMA-IR',nm[9:length(nm)])
colnames(tv_all)=nm; #tv_all[1:5,1:30]; #NAFLD[1:2,1:28];
colnames(tv_all)[colnames(tv_all)=='MASLD[, \"PFAS\"]']='PFAS';
# head(tv_all) #non nans 
which(is.na(tv_all))
colnames(tv_all)

#jälkeenpäin lienee jeesh
x5=x5[x5!='PFAS'];x5=x5[x5!='Perfluorodecyl.ethanoic.acid']; x6=x6[x6!='Total_TG'] # x1;x2;x3;x4;x5;
tv_all=tv_all[,!colnames(tv_all) %in% c('Total_TG','PFAS',"Perfluorodecyl.ethanoic.acid")]

tv_half_log22=cbind(tv[,1],tv_half_log2);
x1=colnames(tv_half_log22[,c(1:8)]); v2=dim(NAFLD)[2]+1
x2=colnames(tv_half_log22[,9:v2]);v3=(dim(Bali)[2]+v2);x3=colnames(tv_half_log22[,(v2+1):(v3)]);v4=(dim(Base)[2])+v3
x4=colnames(tv_half_log22[,(v3+1):(v4-1)]);
x5=colnames(tv_half_log22[,(v4):(dim(tv_half_log22)[2])]);
x3 <- paste(x3, "_L", sep="") #https://stackoverflow.com/questions/6984796/how-to-paste-a-string-on-each-element-of-a-vector-of-strings-using-apply-in-r
x4=gsub("(-[0-9]*)*.1", "", x4) #https://stackoverflow.com/questions/18997297/remove-ending-of-string-with-gsub
x4 <- paste(x4, "_S", sep="")# https://rdrr.io/bioc/qpgraph/man/qpNrr.html
x5a=x5[1:9]
x6=x5[10:length(x5)] #dividing to lipids
x5=x5a  #making sure that PFAS are separate
nm = c(x1,x2,x3,x4,x5,x6); nm=c('PatientNumber','Gender','AGE','BMI','Steatosis Grade','Fibrosis Stage','Necroinflammation','HOMA-IR',nm[9:length(nm)])
colnames(tv_half_log22)=nm; #tv_half_log22[1:5,1:30]; #NAFLD[1:2,1:28];
colnames(tv_half_log22)[colnames(tv_half_log22)=='MASLD[, \"PFAS\"]']='PFAS';
colnames(tv_half_log22)

#jälkeenpäin lienee jeesh
x5=x5[x5!='PFAS'];x5=x5[x5!='Perfluorodecyl.ethanoic.acid']; x6=x6[x6!='Total_TG'] # x1;x2;x3;x4;x5;
tv_half_log22=tv_half_log22[,!colnames(tv_half_log22) %in% c('Total_TG','PFAS',"Perfluorodecyl.ethanoic.acid")]
ccova=tv[,c("Steatosis.Grade.0.To.3" , "Fibrosis.Stage.0.to.4" ,"Necroinflammation" ,  "HOMA-IR")]
sick_group=rowSums(ccova)>4
file_names=c("Steatosis" , "Fibrosis" ,"Necroinflammation" ,  "HOMAIR", 'Menopause')
lkm=30; simss=100 # 
ccovae=tv[,c("Steatosis.Grade.0.To.3")]; sick_group=ccovae>0 #toth# # hist(ccovae,breaks=100) # hist(ccova[,'HOMA-IR'],breaks=100)
fn=file_names[1]; 
t.val='no'; name='Generic'; sick='no'; Group='All';joo='joo';ip=1 #total, total (i.e. all samples (or all samples not sick) lähtökohtaisesti...

Treatment=colnames(tv_all)[71:77];
Mediator=colnames(tv_all)[9:28];
Outcome=colnames(tv_all)[c(29:51,78:90)]; ##https://sparkbyexamples.com/r-programming/r-remove-from-vector-with-examples/
Treatment=Treatment[!Treatment %in% c('Perfluorodecyl.ethanoic.acid')]
tv_all=tv_all[,!colnames(tv_all) %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]
tv_all=tv_all[,!colnames(tv_all) %in% x4]
tv_all=tv_all[,!colnames(tv_all) %in% c('TG_SFA','MUFA','TG_PUFA')]
Outcome=Outcome[!Outcome %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]
Outcome=Outcome[! Outcome %in% x4] #https://sparkbyexamples.com/r-programming/r-remove-from-vector-with-examples/
Outcome=Outcome[! Outcome %in% c('TG_SFA','MUFA','TG_PUFA')] #

colnames(tv_all)
Treatment
Mediator
Outcome



#The basic hypothesis. All are variables (y~x+m;m~x)
loop_med_simplified1a=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  if (Group=='female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else   {tv_red=tv_all[cond,]} 
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  # cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment
  Data=data.frame(Data)
  # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  # this time the b and c are in the loop, and b is the model 'M', i.e. M~X (e.g. Allocholic acid ~ PFOA_L)
  # the c is the model 'Y', i.e. Y~M+X, e.g. CAR ~PFOA + Allocholic acid (note, just one mediator at the time... )
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2])) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()  
  j=1;i=1;z=1
  # simss=2; length(y)*length(m)*length(x)
  for (i in 1:length(y)) {
    for (j in 1:length(m)) { #control.value=mina[i]
      for (z in 1:length(x)) {
        fmla1 <- as.formula(paste(paste(m[j], collapse= "+")," ~ ", paste(x[z], collapse= "+"))) 
        b = lm(fmla1, Data)
        xm=paste(paste(c(x[z],m[j]), collapse= "+"))
        fmla2 <- as.formula(paste(y[i]," ~ ", xm)) #https://www.statology.org/glm-vs-lm-in-r/
        c = lm(fmla2, Data) 
        if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax') 
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )} 
        
        med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
        tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
        res <- rbind(res,tmp);
        rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
        remove(tmp) }}} #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)}

# https://www.degruyter.com/document/doi/10.1515/ijb-2019-0088/html?lang=en
# https://www.degruyter.com/document/doi/10.1515/ijb-2019-0146/html

# #Basic with covariates. All are variables (y~x+m+covas;m~x+covas)
# loop_med_simplified2a = function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
#   if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
#   tv_red=c(); 
#   if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else   {tv_red=tv_all[cond,]} 
#   X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
#   M <- tv_red[,Mediator]  #
#   Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
#   cova <- tv_red[,c('AGE','BMI','Gender')] 
#   Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
#   Data=data.frame(Data)
#   colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment# https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
#   control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
#   treat.value=colMaxs(as.matrix(X))
#   #M~X
#   x=c();m=c(); y=c();ye=c()
#   for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
#   for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
#   #Y~X+M
#   for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
#   med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()  
#   j=1;i=1;z=1
#   # simss=10; length(y)*length(m)*length(x)
#   for (i in 1:length(y)) {
#     for (j in 1:length(m)) { #control.value=mina[i]
#       for (z in 1:length(x)) {
#         if (Group=='All') {fmla1 <- as.formula(paste(paste(m[j], collapse= "+")," ~ ", paste(c(paste(x[z], collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")))} else {
#           fmla1 <- as.formula(paste(paste(m[j], collapse= "+")," ~ ", paste(c(paste(x[z], collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")))} 
#         b = lm(fmla1, Data)  
#         if (Group=='All') {xm=paste(c(paste(c(x[z],m[j]), collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")} else
#         {xm=paste(c(paste(c(x[z],m[j]), collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")} 
#         fmla2 <- as.formula(paste(y[i]," ~ ", xm))
#         c = lm(fmla2, Data) 
#         if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
#         {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax')
#         {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )}
#         
#         med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
#         tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
#         res <- rbind(res,tmp);
#         rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
#         remove(tmp) }}} #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
#   
#   rownames(res)=rn #write.csv(rn,'iii.csv')
#   colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
#   res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
#   rownames(res) <- gsub("X11", "11", rownames(res))
#   rownames(res) <- gsub("X17", "17", rownames(res))
#   write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE)  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
#   return(res)}

#Basic with covariates and ...variable selections not needed. :) All are variables (y~x+m+covas;m~x+covas)
loop_med_simplified2b = function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c();
  if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else   {tv_red=tv_all[cond,]}
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"
  cova <- tv_red[,c('AGE','BMI','Gender')]
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment# https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()
  j=1;i=1;z=1
  # simss=10; length(y)*length(m)*length(x)
  
  try({ 
  for (i in 1:length(y)) {
    for (j in 1:length(m)) { #control.value=mina[i]
      for (z in 1:length(x)) {
        qqq=paste(c("Data[, ",as.numeric(dim(Data)[2])-3,"]"),collapse="")
        qq=paste(c("Data[, ",as.numeric(dim(Data)[2])-2,"]"),collapse="")
        q=paste(c("Data[, ",as.numeric(dim(Data)[2])-1,"]"),collapse="")

        
        if (Group=='All') {fmla1 <- as.formula(paste(paste(m[j], collapse= "+")," ~ ", paste(c(paste(x[z], collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else {
          fmla1 <- as.formula(paste(paste(m[j], collapse= "+")," ~ ", paste(c(paste(x[z], collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")))}
      
        b = lm(fmla1, Data)
        if (Group=='All') {xm=paste(c(paste(c(x[z],m[j]), collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")} else
        {xm=paste(c(paste(c(x[z],m[j]), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")}
        fmla2 <- as.formula(paste(y[i]," ~ ", xm))
        c = lm(fmla2, Data)
        if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax')
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )}

        med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
        tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2])
        res <- rbind(res,tmp);
        rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
        remove(tmp) }}}} ,  {print('mediation_notOK')})  #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html

  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE)  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  
  return(res)}

#x is stable! (all other are variables)
loop_med_simplified4=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  tv_red=tv_all[cond,] # tv_red=tv_all
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()  
  j=1;i=1;z=1
  # simss=10; length(y)*length(m)*length(x)
  for (i in 1:length(y)) {
    for (j in 1:length(m)) { #control.value=mina[i]
      for (z in 1:length(x)) {
        if (Group=='All') {fmla1 <- as.formula(paste(paste(m[j], collapse= "+")," ~ ", paste(c(paste(x, collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")))} else {
          fmla1 <- as.formula(paste(paste(m[j], collapse= "+")," ~ ", paste(c(paste(x, collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")))} 
        b = lm(fmla1, Data)  
        if (Group=='All') {xm=paste(c(paste(c(x,m[j]), collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")} else
        {xm=paste(c(paste(c(x,m[j]), collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")} 
        fmla2 <- as.formula(paste(y[i]," ~ ", xm))
        c = lm(fmla2, Data) 
        if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax')
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )}
        med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
        tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
        res <- rbind(res,tmp);
        rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
        remove(tmp) }}} #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)}

#x AND m are stable! (y is variable): (musta paras, hypo2)
loop_med_simplified5=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  
  # Group='All'
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  tv_red=tv_all[cond,] # tv_red=tv_all
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  if (Group=='All') {xm=paste(c(paste(c(x,m), collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")} else
  {xm=paste(c(paste(c(x,m), collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")} 
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()  
  # e=c()
  # if (Group=='All') {e=sum(c(names(c$coefficients)=='Data[, 82]', names(c$coefficients)=='Data[, 83]', names(c$coefficients)=='Data[, 84]'))} else 
  # {e=sum(c(names(c$coefficients)=='Data[, 82]', names(c$coefficients)=='Data[, 83]')) }
  # simss=2; length(y)*length(m)*length(x)
  j=1;i=1;z=1 #only for outcome model (c) put the AIC
  
  for (i in 1:length(y)) {
    fmla2 <- as.formula(paste(y[i]," ~ ", xm))
    c = lm(fmla2, Data)
    c=stepAIC(c, k=2) #https://www.scribbr.com/statistics/akaike-information-criterion/
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9]
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8
    if (length(be)<1) next 
    if (length(ce)<1) next 
    if (length(be)<1 & sum(c$coefficients>0.0001)<3) {next} else {
      he=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>28]
      ii=gfg[as.numeric(unlist(gfg_numbers))>28] #c("Data[, 82]","Data[, 83]","Data[, 84]")
      # gfg_numbers2 <- regmatches(ii, gregexpr("[[:digit:]]+", ii)) # iii=c("Data[, 82]","Data[, 83]")
      if (length(he)>0) { fmla1 <- as.formula(paste(paste(m[ce], collapse= "+")," ~ ", paste(c(paste(x[be], collapse= "+"),paste(ii,collapse= "+")),collapse= "+")))} else 
      {fmla1 <- as.formula(paste(paste(m[ce], collapse= "+")," ~ ", paste(x[be], collapse= "+")))} 
      #     else {fmla1 <- as.formula(paste(paste(m[ce], collapse= "+")," ~ ", paste(c(paste(x[be], collapse= "+"),paste(iii,collapse= "+")),collapse= "+")))} } #if (Group=='All') {
      b = lm(fmla1, Data) 
      # b <- MASS::stepAIC(b, k = 2)
      # l=length(names(b$coefficients))
      # gfg=names(b$coefficients)[2:l]
      # gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
      # a=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<10] 
      # t=intersect(a,be)# > i# [1] 48# > j# [1] 15# > z# [1] 5
      try(for (j in ce) { #control.value=mina[i]
        for (z in be) {
          # z=9
          # print(c$coefficients)
          # if (sum(c(names(c$coefficients)==m[j], names(c$coefficients)==x[z], names(b$coefficients)==x[z])>2)) {
          if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
          {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax')
          {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )} #} else {next}
          med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
          tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],
                med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],
                med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
          res <- rbind(res,tmp);
          rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
          remove(tmp) }}) }} #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)}

#x AND m are stable! (y is variable): (musta paras, hypo2)
loop_med_simplified5ö=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  # Group='Female'
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  tv_red=tv_all[cond,] # tv_red=tv_all
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  # if (Group=='All') {xm=paste(c(paste(c(x,m), collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")} else
  # {xm=paste(c(paste(c(x,m), collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")} 
  # xm=paste(c(x,m), collapse= "+")
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()  
  # simss=2; length(y)*length(m)*length(x)
  # j=1;i=1;z=1 #only for outcome model (c) put the AIC
  bvif=c()
  cvif=c()
  
  for (i in 1:length(y)) {
    # xm=paste(c(x,m), collapse= "+"); fmla2 <- as.formula(paste(y[i]," ~ ", xm))
    # c = lm(fmla2, Data)
    # l=length(names(c$coefficients))
    # if (l>2) try({
    yy=as.matrix(Data[,Outcome[i]])
    xx=as.matrix(Data[,c(colnames(data.frame(X)),colnames(data.frame(M)))])
    cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); best_lambda <- cv_model$lambda.min
    best_model <- glmnet(xx, yy, alpha = 1, lambda = best_lambda) #coef(best_model)[abs(coef(best_model))>0.00001]
    if (sum(abs(coef(best_model))>0.0001)>2) {
      cbm=coef(best_model)[2:length(coef(best_model)),]
      ok=abs(cbm)>0.0001 
      äm=names(cbm)[ok]
      xm=c()
      for (i in which(colnames(Data) %in% äm)) {xm=append(xm,paste("Data[, ",i , "]", sep=""))}
      xm=paste(xm, collapse= "+"); 
      fmla2 <- as.formula(paste(y[i]," ~ ", xm))
      c = lm(fmla2, Data)} else next
    # },next) 
    
    # if (l>2) try({c <- stepAIC(c, k=1.5)},next)
    #https://www.scribbr.com/statistics/akaike-information-criterion/
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9]
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8
    if (length(be)<1) next 
    if (length(ce)<1) next 
    if (sum(abs((c$coefficients)[2:l])>0.0005)<3) next 
    
    if (Group=='All') {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")} else 
    {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")} 
    fmla2 <- as.formula(paste(y[i]," ~ ", xm))
    c = lm(fmla2, Data)
    cvif = append(cvif,vif(c))  
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9]
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8
    
    fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+")))))
    b = lm(fmla1, Data)
    l=length(names(b$coefficients))
    if (l>2) try({
      yy=rowMedians(as.matrix(Data[,colnames(data.frame(M))]))
      xx=as.matrix(Data[,colnames(data.frame(X))])
      cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); best_lambda <- cv_model$lambda.min
      best_model <- glmnet(xx, yy, alpha = 1, lambda = best_lambda) #coef(best_model)[abs(coef(best_model))>0.00001]
      if (sum(abs(coef(best_model))>0.0001)>2) {
        cbm=coef(best_model)[2:length(coef(best_model)),]
        ok=abs(cbm)>0.0001
        äm=names(cbm)[ok]
        x=c()
        for (i in which(colnames(Data) %in% äm)) {x=append(x,paste("Data[, ",i , "]", sep=""))}
        fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+")))))
        b = lm(fmla1, Data)} else {next}})
    #
    
    l=length(names(b$coefficients))
    gfg=names(b$coefficients)[2:l]
    # if (l>2) try({
    if (Group=='All') { fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")))} else
    {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")))};
    b = lm(fmla1, Data) #}))
    bvif = append(bvif,vif(b))
    lu=length(names(b$coefficients))
    gfgu=names(b$coefficients)[2:lu]
    gfgu_numbers <- regmatches(gfgu, gregexpr("[[:digit:]]+", gfgu))
    beu=as.numeric(unlist(gfgu_numbers))[as.numeric(unlist(gfgu_numbers))<9]
    be=intersect(be,beu)
    
    try(
      for (j in ce) { #control.value=mina[i]
        for (z in be) {
          # z=9
          # print(c$coefficients)
          # if (sum(c(names(c$coefficients)==m[j], names(c$coefficients)==x[z], names(b$coefficients)==x[z])>2)) {
          if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
          {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax')
          {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )} #} else {next}
          med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
          tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],
                med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],
                med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
          res <- rbind(res,tmp);
          rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
          remove(tmp) }}) 
  } #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)} #list(res,bvif,cvif)

#x AND m are stable! (y is variable): (musta paras, hypo2)
loop_med_simplified5ö2=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  Group='Female'
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  tv_red=tv_all[cond,] # tv_red=tv_all
  X <- tv_red[,Treatment] #Standard values did not five errors # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  # if (Group=='All') {xm=paste(c(paste(c(x,m), collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")} else
  # {xm=paste(c(paste(c(x,m), collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")} 
  # xm=paste(c(x,m), collapse= "+")
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()  
  # simss=2; length(y)*length(m)*length(x)
  j=1;i=1;z=1 #only for outcome model (c) put the AIC
  for (i in 1:length(y)) {
    
    xm=paste(c(x,m), collapse= "+")
    ye=paste(y, collapse= "+")
    fmla2 <- as.formula(paste(ye," ~ ", xm))
    c = lm(fmla2, Data)
    l=length(names(c$coefficients))
    if (l>2) try({c <- stepAIC(c, k=1.5)},next)
    #https://www.scribbr.com/statistics/akaike-information-criterion/
    
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9]
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8
    if (length(be)<1) next 
    if (length(ce)<1) next 
    if (sum(abs((c$coefficients)[2:l])>0.0005)<3) next 
    
    if (Group=='All') {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")} else
    {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")} 
    # ye=paste(y, collapse= "+")
    fmla2 <- as.formula(paste(y[i]," ~ ", xm))
    c = lm(fmla2, Data)
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9]
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8
    
    fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+")))))
    b = lm(fmla1, Data) 
    l=length(names(b$coefficients))
    if (l>2) {b <- stepAIC(b, k=1)}
    l=length(names(b$coefficients))
    gfg=names(b$coefficients)[2:l]
    if (Group=='All') { fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")))} else
    {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")))}
    b = lm(fmla1, Data)
    lu=length(names(b$coefficients))
    gfgu=names(b$coefficients)[2:lu]
    gfgu_numbers <- regmatches(gfgu, gregexpr("[[:digit:]]+", gfgu))
    beu=as.numeric(unlist(gfgu_numbers))[as.numeric(unlist(gfgu_numbers))<10]
    be=intersect(be,beu)

    try(
      for (j in ce) { #control.value=mina[i]
        for (z in be) {
          # z=9
          # print(c$coefficients)
          # if (sum(c(names(c$coefficients)==m[j], names(c$coefficients)==x[z], names(b$coefficients)==x[z])>2)) {
          if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
          {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax')
          {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )} #} else {next}
          med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
          tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],
                med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],
                med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
          res <- rbind(res,tmp);
          rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...colnames(Y)[i]
          remove(tmp) }}) 
  } #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)}

#x AND m are stable! (y is variable): (musta paras, hypo2)
loop_med_simplified5ö3=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  # Group='All'
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  tv_red=tv_all[cond,] # tv_red=tv_all
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()  
  j=1;i=1;z=1
  # simss=10; length(y)*length(m)*length(x)
  for (i in 1:length(y)) {
    for (j in 1:length(m)) { #control.value=mina[i]
      for (z in 1:length(x)) {
        if (Group=='All') {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")))} else {
          fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")))} 
        b = lm(fmla1, Data)  
        if (Group=='All') {xm=paste(c(paste(c(x,m), collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")} else
        {xm=paste(c(paste(c(x,m), collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")} 
        fmla2 <- as.formula(paste(y[i]," ~ ", xm))
        c = lm(fmla2, Data) 
        if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax')
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )}
        med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
        tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
        res <- rbind(res,tmp);
        rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
        remove(tmp) }}} #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  #https://stackoverflow.com/questions/21937830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)}

#x AND m are stable! (y is variable): (musta paras, hypo2)
loop_med_simplified5ö4=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick, sick_group) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  # Group='All'
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if 
  (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else   {tv_red=tv_all[cond,]} # tv_red=tv_all
  # tv_red=tv_all #male 22 91, female 26 91: sum(as.vector(sick_group))
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()
  bvif=c()
  cvif=c()
  # colnames(Data);head(Data)
  
  for (i in 1:length(y)) {
    xm=paste(c(x,m), collapse= "+"); fmla2 <- as.formula(paste(y[i]," ~ ", xm))
    c = lm(fmla2, Data)
    l=length(names(c$coefficients))
    if (l>2) try({
      yy=as.matrix(Data[,Outcome[i]])
      xx=as.matrix(Data[,c(colnames(data.frame(X)),colnames(data.frame(M)))])
      cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); best_lambda <- cv_model$lambda.min
      best_model <- glmnet(xx, yy, alpha = 1, lambda = best_lambda) #coef(best_model)[abs(coef(best_model))>0.00001]
      if (sum(abs(coef(best_model))>0.0001)>2) {
        cbm=coef(best_model)[2:length(coef(best_model)),]
        ok=abs(cbm)>0.0001 
        äm=names(cbm)[ok]
        xm=c()
        for (i in which(colnames(Data) %in% äm)) {xm=append(xm,paste("Data[, ",i , "]", sep=""))}
        xm=paste(xm, collapse= "+"); 
        fmla2 <- as.formula(paste(y[i]," ~ ", xm))
        c = lm(fmla2, Data)} else next},next) 
    
    # if (l>2) try({c <- stepAIC(c, k=1.5)},next)
    #https://www.scribbr.com/statistics/akaike-information-criterion/
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9]
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8
    if (length(be)<1) next 
    if (length(ce)<1) next 
    if (sum(abs((c$coefficients)[2:l])>0.0005)<3) next 
    
    if (Group=='All') {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")} else 
    {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")} 
    fmla2 <- as.formula(paste(y[i]," ~ ", xm))
    c = lm(fmla2, Data)
    cvif = append(cvif,vif(c))  
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9] #the contaminants
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8 #the mediators
    
    fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+")))))
    b = lm(fmla1, Data)
    l=length(names(b$coefficients))
    if (l>2) try({
      yy=rowMedians(as.matrix(Data[,colnames(data.frame(M))]))
      xx=as.matrix(Data[,colnames(data.frame(X))])
      cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); best_lambda <- cv_model$lambda.min
      best_model <- glmnet(xx, yy, alpha = 1, lambda = best_lambda) #coef(best_model)[abs(coef(best_model))>0.00001]
      if (sum(abs(coef(best_model))>0.0001)>2) {
        cbm=coef(best_model)[2:length(coef(best_model)),]
        ok=abs(cbm)>0.0001
        äm=names(cbm)[ok]
        xe=c()
        for (i in which(colnames(Data) %in% äm)) {xe=append(xe,paste("Data[, ",i , "]", sep=""))}
        fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(xe, collapse= "+")))))
        b = lm(fmla1, Data)} else {next}})
    #
    
    l=length(names(b$coefficients))
    gfg=names(b$coefficients)[2:l]
    # if (l>2) try({
    if (Group=='All') { fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")))} else
    {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")))};
    b = lm(fmla1, Data) #}))
    bvif = append(bvif,vif(b))
    lu=length(names(b$coefficients))
    gfgu=names(b$coefficients)[2:lu]
    gfgu_numbers <- regmatches(gfgu, gregexpr("[[:digit:]]+", gfgu))
    beu=as.numeric(unlist(gfgu_numbers))[as.numeric(unlist(gfgu_numbers))<9]
    be=intersect(be,beu) #contaminants
    
    # m=m[1:length(ce)]
    # x=x[1:length(be)]
    
    try(
      for (j in 1:length(ce)) { #control.value=mina[i]
        for (z in 1:length(be)) {
          # z=9
          # print(c$coefficients)
          # if (sum(c(names(c$coefficients)==m[j], names(c$coefficients)==x[z], names(b$coefficients)==x[z])>2)) {
          if (t.val=='no'){med_oute=mediate(b, c, treat =  x[be[z]], mediator = m[ce[j]],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
          {med_oute=mediate(b, c, treat =  x[be[z]], mediator = m[ce[j]],sims = simss,control.value=control.value[be[z]],treat.value=X[test,be[z]] )} else if (t.val=='minmax')
          {med_oute=mediate(b, c, treat =  x[be[z]], mediator = m[ce[j]],sims = simss,control.value=control.value[be[z]],treat.value=treat.value[be[z]])} #} else {next}
          med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
          tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],
                med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],
                med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
          res <- rbind(res,tmp);
          rn=append(rn,paste(colnames(X)[be[z]],colnames(M)[ce[j]],colnames(Y)[i], sep=" ")) #attaching two rownames...
          remove(tmp) }}) 
  } #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)} #list(res,bvif,cvif)

#x AND m are stable! (y is variable): (musta paras, hypo2)
loop_med_simplified5a=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  tv_red=tv_all[cond,] # tv_red=tv_all
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  # cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment
  Data=data.frame(Data)
  # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2])) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c() 
  j=1;i=1;z=1
  # simss=2; length(y)*length(m)*length(x)
  
  # https://fmch.bmj.com/content/8/1/e000262
  
  for (i in 1:length(y)) {
    for (j in 1:length(m)) { #control.value=mina[i]
      for (z in 1:length(x)) {
        fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(x, collapse= "+"))) 
        b = lm(fmla1, Data)
        xm=paste(paste(c(x,m), collapse= "+"))
        fmla2 <- as.formula(paste(y[i]," ~ ", xm)) #https://www.statology.org/glm-vs-lm-in-r/
        c = lm(fmla2, Data) 
        if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax') 
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )} 
        
        med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
        tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
        res <- rbind(res,tmp);
        rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
        remove(tmp) }}} #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)}

#x and y vaihtelee (tai ovat stabiileja, pointsi on se, että m on aina sama)
loop_med_simplified6=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  tv_red=tv_all[cond,] # tv_red=tv_all
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment# https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()  
  j=1;i=1;z=1
  # simss=10; length(y)*length(m)*length(x)
  for (i in 1:length(y)) {
    for (j in 1:length(m)) { #control.value=mina[i]
      for (z in 1:length(x)) {
        if (Group=='All') {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x[z], collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")))} else {
          fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x[z], collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")))} 
        b = lm(fmla1, Data)  
        if (Group=='All') {xm=paste(c(paste(c(x[z],m), collapse= "+"),paste(c("Data[, 82]","Data[, 83]","Data[, 84]"),collapse= "+")),collapse= "+")} else
        {xm=paste(c(paste(c(x[z],m), collapse= "+"),paste(c("Data[, 82]","Data[, 83]"),collapse= "+")),collapse= "+")} 
        fmla2 <- as.formula(paste(y[i]," ~ ", xm))
        c = lm(fmla2, Data) 
        if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax')
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )}
        
        med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
        tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
        res <- rbind(res,tmp);
        rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
        remove(tmp) }}} #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE)  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)}

#x and y vaihtelee (tai ovat stabiileja, pointsi on se, että m on aina sama)
loop_med_simplified6a=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  tv_red=tv_all[cond,] # tv_red=tv_all
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  # cova <- tv_red[,c('AGE','BMI','Gender')] 
  # Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment# https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2])) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()  
  # j=1;i=1;z=1
  # simss=10; length(y)*length(m)*length(x)
  for (i in 1:length(y)) {
    for (j in 1:length(m)) { #control.value=mina[i]
      for (z in 1:length(x)) {
        if (Group=='All') {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x[z], collapse= "+")),collapse= "+")))} else {
          fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x[z], collapse= "+")),collapse= "+")))} 
        b = lm(fmla1, Data)  
        if (Group=='All') {xm=paste(c(paste(c(x[z],m), collapse= "+")))} else
        {xm=paste(c(paste(c(x[z],m), collapse= "+")))} 
        fmla2 <- as.formula(paste(y[i]," ~ ", xm))
        c = lm(fmla2, Data) 
        if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax')
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )}
        
        med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
        tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
        res <- rbind(res,tmp);
        rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
        remove(tmp) }}} #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE)  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)}

#x is stable! (all other are variables), ver ilman cov. y ja m vaihtelee, x ei, ei sama kuin yllä, tässä EI ole kovariaatteja
loop_med_simplified7=function(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test) { 
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  tv_red=tv_all[cond,] # tv_red=tv_all
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  # cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2])) {y=append(y,paste("Data[, ",z , "]", sep="")) } 
  #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()  
  j=1;i=1;z=1
  # simss=10; length(y)*length(m)*length(x)
  for (i in 1:length(y)) {
    for (j in 1:length(m)) { #control.value=mina[i]
      for (z in 1:length(x)) {
        fmla1 <- as.formula(paste(paste(m[j], collapse= "+")," ~ ", paste(x, collapse= "+"))) 
        b = lm(fmla1, Data)
        xm=paste(paste(c(x,m[j]), collapse= "+"))
        fmla2 <- as.formula(paste(y[i]," ~ ", xm)) #https://www.statology.org/glm-vs-lm-in-r/
        c = lm(fmla2, Data) 
        if (t.val=='no'){med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss)} else if (t.val=='yes') # control.value=control.value[z],treat.value=treat.value[z]  
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=X[test,z] )} else if (t.val=='minmax') 
        {med_oute=mediate(b, c, treat =  x[z], mediator = m[j],sims = simss,control.value=control.value[z],treat.value=treat.value[z] )} 
        
        med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
        tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
        res <- rbind(res,tmp);
        rn=append(rn,paste(colnames(X)[z],colnames(M)[j],colnames(Y)[i], sep=" ")) #attaching two rownames...
        remove(tmp) }}} #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  
  rownames(res)=rn #write.csv(rn,'iii.csv')
  colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  rownames(res) <- gsub("X11", "11", rownames(res))
  rownames(res) <- gsub("X17", "17", rownames(res))
  write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)}


#Hypoteesi 3
loop_med_simplified5ö5=function(Treatment, Mediator, Outcome,tv_all,Group,name,date,simss,t.val,test,sick, sick_group) { # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  
  #muuten hyvä, mutta _S:t katoaa...
  #nyt saattanee toimia.. hiotaan vielä funktiota, tikka9124
  # Group='All'
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if 
  (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else  {tv_red=tv_all} #{tv_red=tv_all[cond & !as.vector(sick_group),]} # tv_red=tv_all; oh, uusix...
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()
  # simss=3; length(y)*length(m)*length(x)
  ih=1;il=1;i=1; j=1;z=1 #only for outcome model (c) put the AIC
  # Y[43:53]
  # bvif=c()  # cvif=c()
  # colnames(Data);head(Data) # dim(Data)
  # i=1
  
  for (i in 1:length(y)) {
    # print(Outcome[i])
    # print(i)
    yy=as.matrix(Data[,Outcome[i]])
    xx=as.matrix(Data[,c(colnames(data.frame(X)),colnames(data.frame(M)))])
    cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); 
    df <- data.frame(x = seq(1:length(cv_model$lambda)),y = cv_model$lambda)#and also: hei=elbow(data = df)
    yn=find_curve_elbow(df)#or: which(hei$data[,4]==min(hei$data[,4]))[1]
    df2=df[yn:dim(df)[1],]#hei$data[yn:dim(hei$data)[1],4]# df2 <- data.frame(x = seq(1:length(df2)),y = df2)
    yn2=find_curve_elbow(df2) #plot_curve = FALSE
    
    yti=yn*0.95
    jeps=cv_model$lambda[yti]
    best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps) #coef(best_model)[abs(coef(best_model))>0.00001]
    cbm=coef(best_model)[2:length(coef(best_model)),]
    ok=abs(cbm)!=0
    äm=names(cbm)[ok]
    xm=c()
    for (ih in which(colnames(Data) %in% äm)) {xm=append(xm,paste("Data[, ",ih , "]", sep=""))}
    xm=paste(xm, collapse= "+");
    fmla2 <- as.formula(paste(y[i]," ~ ", xm))
    c = lm(fmla2, Data)
    l=length(names(c$coefficients))
    if (l<3) {xm=paste(c(x,m), collapse= "+"); fmla2 <- as.formula(paste(y[i]," ~ ", xm))
    c = lm(fmla2, Data);c <- stepAIC(c, k=1.2)} #https://www.scribbr.com/statistics/akaike-information-criterion/
    
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9]
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8
    
    qqq=paste(c("Data[, ",as.numeric(dim(Data)[2])-3,"]"),collapse="")
    qq=paste(c("Data[, ",as.numeric(dim(Data)[2])-2,"]"),collapse="")
    q=paste(c("Data[, ",as.numeric(dim(Data)[2])-1,"]"),collapse="")
    if (Group=='All') {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")} else 
    {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")}  
    fmla2 <- as.formula(paste(y[i]," ~ ", xm))
    c = lm(fmla2, Data)
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    # print(gfg)
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9] #the contaminants
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8 #the mediators
    
    if(sum(is.na(summary(c)[4]$coefficients[,4]))==length(summary(c)[4]$coefficients[,4])) {
      res <- rbind(res,rep(-1000, 16));
      # print(Outcome[i])
      rn=append(rn,paste('no_X','no_M',Outcome[i], sep=" ")); print('hei_c');next}
    
    yy=rowMedians(as.matrix(Data[,colnames(data.frame(M))]))
    xx=as.matrix(Data[,colnames(data.frame(X))])
    cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); 
    df <- data.frame(x = seq(1:length(cv_model$lambda)),y = cv_model$lambda)#and also: hei=elbow(data = df)
    yn=find_curve_elbow(df)#or: which(hei$data[,4]==min(hei$data[,4]))[1]
    df2=df[yn:dim(df)[1],]#hei$data[yn:dim(hei$data)[1],4]# df2 <- data.frame(x = seq(1:length(df2)),y = df2)
    yn2=find_curve_elbow(df2) #plot_curve = FALSE
    yti=yn+ceiling(yn2)
    jeps=cv_model$lambda[yti]
    best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps)
    
    if (sum(abs(coef(best_model))>0.0001)>1) {
      cbm=coef(best_model)[2:length(coef(best_model)),]
      ok=abs(cbm)!=0; äm=names(cbm)[ok]; xe=c()
      for (il in which(colnames(Data) %in% äm)) {xe=append(xe,paste("Data[, ",il , "]", sep=""))}
      fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(xe, collapse= "+")))))
      b = lm(fmla1, Data)} else {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+")))))
      b = lm(fmla1, Data)
      b <- stepAIC(b, k=1.2)} #https://bookdown.org/egarpor/PM-UC3M/lm-ii-modsel.html
    
    l=length(names(b$coefficients))
    gfg=names(b$coefficients)[2:l]
    # print(gfg)
    if (length(gfg)>0) {
      if (Group=='All') 
      { fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else
      {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")))};
      b = lm(fmla1, Data)} else {
        if (Group=='All') { fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+"),
                                                                                             paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else
                                                                                             {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+"),
                                                                                                                                                               paste(c(qqq,qq),collapse= "+")),collapse= "+")))};
        b = lm(fmla1, Data);b <- stepAIC(b, k=0.4)}# bvif = append(bvif,vif(b))
    lu=length(names(b$coefficients))
    gfgu=names(b$coefficients)[2:lu]
    # print(gfgu)
    gfgu_numbers <- regmatches(gfgu, gregexpr("[[:digit:]]+", gfgu))
    beu=as.numeric(unlist(gfgu_numbers))[as.numeric(unlist(gfgu_numbers))<9]
    be=intersect(be,beu) #contaminants
    print(be);print(i)
    
    # if(sum(is.na(summary(b)[4]$coefficients[,4]))==length(summary(b)[4]$coefficients[,4])) {print('hei_b');next}
    
    if (length(be)<1) { res <- rbind(res,rep(-1000, 16));
    # print(Outcome[i])
    rn=append(rn,paste('no_X','no_M',Outcome[i], sep=" ")); print('hei_b');next}
    
    # x=x[1:length(be)]
    
    try(
      for (j in 1:length(ce)) { #control.value=mina[i]
        for (z in 1:length(be)) {
          med_oute=mediate(b, c, treat =  x[be[z]], mediator = m[ce[j]],sims = simss,control.value=control.value[be[z]],treat.value=treat.value[be[z]]) #} else {next} control.value=control.value[be[z]],treat.value=treat.value[be[z]]
          med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
          tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],
                med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],
                med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
          res <- rbind(res,tmp);
          
          # print(Outcome[i])
          rn=append(rn,paste(colnames(X)[be[z]],colnames(M)[ce[j]],Outcome[i], sep=" ")) #attaching two rownames...
          remove(tmp) }},  {print('medc');next})} 
  #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html 
  # rownames(res)=rn #write.csv(rn,'iii.csv')
  # colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  # res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  # rownames(res) <- gsub("X11", "11", rownames(res))
  # rownames(res) <- gsub("X17", "17", rownames(res))
  # head(res)
  # # Group='Male'
  # write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  try({res},{res=t(data.frame(rep(0,16)))})
  if (!isTRUE(res) & dim(res)[1]<2 & dim(res)[2]<2) {res=t(data.frame(rep(0,16)))}
  try({rownames(res)=rn}, {print('rn')}) #write.csv(rn,'iii.csv')
  try({colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  }, {print('col')})
  try({res=res[order(res[,2]),]}) #res=res[rev(order(res[,1])),]
  try({rownames(res) <- gsub("X11", "11", rownames(res))}, {print('xx')})
  try({rownames(res) <- gsub("X17", "17", rownames(res)) }, {print('xxx')})
  # # head(res)
  # # Group='Male'
  # tryCatch({write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) }, error = function(msg){
  #   return(write.csv(res,paste(name,Group,dhist(rowSums(ccova),breaks=100)ate,'.csv')))})
  try({write.xlsx(res, file = paste(name,Group,date,'test1.xlsx'), append = FALSE, row.names = TRUE)},{print('xls')})
  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)} #list(res,bvif,cvif)

# Group='All'
Group='Male'
sick='yes';joo='no';ip=1
ccovae=tv[,c("Steatosis.Grade.0.To.3")]; sick_group=ccovae>0 #toth# # hist(ccovae,breaks=100) # hist(ccova[,'HOMA-IR'],breaks=100)
file_names=c("Steatosis" , "Fibrosis" ,"Necroinflammation" ,  "HOMAIR", 'Menopause')
fn=file_names[1]; 
hoi1=paste("C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/",fn,sep='')
setwd(hoi1) #
simss=10



# Group,name,date,simss,sick,sick_group,joo,ip... on devi but seemingly ok:
loop_med_simplified_hyp4a=function(Group,tv_all,name,date,simss,sick,sick_group,joo,ip) { 
    # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
    # muuten hyvä, mutta _S:t katoaa... #nyt saattanee toimia.. hiotaan vielä funktiota, tikka1124.. lienee ok
    # library(pathviewr) #Group='Male'
    if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if 
    (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])} 
    tv_red=c(); 
    # !as.vector(sick_group)
    if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else   {
      if (joo =='joo') {tv_red=tv_all} else {tv_red=tv_all[cond & !as.vector(sick_group),]}} #tv_red=tv_all[cond & !as.vector(sick_group),]} # tv_red=tv_all, I wonder this healthy should be 'cond' non-sick rather than all...
    X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
    M <- ip*tv_red[,Mediator]  #
    Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
    cova <- tv_red[,c('AGE','BMI','Gender')] 
    Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
    Data=data.frame(Data)
    colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
    control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
    treat.value=colMaxs(as.matrix(X))
    #M~X
    x=c();m=c(); y=c();ye=c()
    for (ii in 1:length(colnames(X))) {x=append(x,paste("Data[, ",ii , "]", sep=""))}
    for (jn in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",jn , "]", sep="")) }
    #Y~X+M
    for (zr in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",zr , "]", sep="")) } #this dimension was essential for the loop names
    med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()
    # simss=3; length(y)*length(m)*length(x)
    ih=1;il=1;i=1; j=1;z=1 #only for outcome model (c) put the AIC
    # control=as.matrix(X)>colMedians(as.matrix(X))
  
  # i=4
  # simss=20
  for (i in 1:length(y)) {# 
    # print(Outcome[i])
    Data=c()
    if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if
    (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
    tv_red=c();
    if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else  {tv_red=tv_all} #{tv_red=tv_all[cond,]} # tv_red=tv_all
    X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
    M <- ip*tv_red[,Mediator]  #
    
    colnames(M) <- gsub("11", "X11", colnames(M));
    colnames(M) <- gsub("17", "X17", colnames(M));
    colnames(M) <- gsub("-", ".", colnames(M));
    colnames(M) <- gsub("/", ".", colnames(M));
    joi=dim(M)[2]+dim(X)[2]
    å=length(y); list=1:å; list1=list[1:å!=i];list2=list1+joi;
    
    Y <- data.frame(tv_red[,Outcome[i]]); colnames(Y)=Outcome[i]  #-?
    Y2 <- -tv_red[,Outcome[1:å!=i]]; 
    Y3 <- tv_red[,Outcome[1:å!=i]]    #"Steatosis.Grade.0.To.3","Fibrosis.Stage.0.to.4","Necroinflammation","HOMA-IR"
    cova <- tv_red[,c('AGE','BMI','Gender')]
    Data <- cbind(X,M,Y,Y2,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
    Data=data.frame(Data);
    colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
    
    Data2 <- cbind(X,M,Y,Y3,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
    Data2=data.frame(Data2);
    colnames(Data2) <- gsub(" ", "_", colnames(Data2)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
    
    yy=as.matrix(Data[,Outcome[i]])
    xx=as.matrix(Data[,c(colnames(data.frame(X)),colnames(data.frame(M)),colnames(data.frame(Y2)))]) #this is ok!
    
    # cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); 
    try({cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1)} ,  {print('glmnet_na');next}) 
    df <- data.frame(x = seq(1:length(cv_model$lambda)),y = cv_model$lambda) #and also: hei=elbow(data = df)
    #to get the optimal amount of variables this is needed:
    yn=find_curve_elbow(df)#or: which(hei$data[,4]==min(hei$data[,4]))[1]
    df2=df[yn:dim(df)[1],]#hei$data[yn:dim(hei$data)[1],4]# df2 <- data.frame(x = seq(1:length(df2)),y = df2)
    yn2=find_curve_elbow(df2) #plot_curve = FALSE
    yti=yn*0.95
    jeps=cv_model$lambda[yti] #coef(best_model)[abs(coef(best_model))>0.00001]
    # best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps) 
    try({best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps) } ,  {print('glmnet_na');next}) 
    cbm=coef(best_model)[2:length(coef(best_model)),]
    ok=abs(cbm)!=0
    äm=names(cbm)[ok]
    xmo=c(); for (ih in which(colnames(Data) %in% äm)) {xmo=append(xmo,paste("Data[, ",ih , "]", sep=""))}
    xmo=colnames(Data)[colnames(Data) %in% äm]
    if (sum(colnames(X) %in% xmo)==0) {xmo=c(xmo,colnames(X))}
    # if (length(c$coefficients)>0) {xmo=xmo[!is.na(c$coefficients[2:l])]} else {xmo=xmo}
    xmoa=paste(xmo, collapse= "+");
    fmla2 <- as.formula(paste(Outcome[i]," ~ ", xmoa))
    c = lm(fmla2, Data)
    l=length(names(c$coefficients))
    
    if (l<3) {
      xmo=paste(c(x,m,y[list1]), collapse= "+"); fmla2 <- as.formula(paste(Outcome[i]," ~ ", xmo)); 
      c = lm(fmla2, Data);
      if (sum(is.na(c$coefficients))>0) {
        l=length(names(c$coefficients))
        xmo=c(x,m,y[list1])[!is.na(c$coefficients[2:l])]
        xmoa=paste(xmo, collapse= "+");
        fmla2 <- as.formula(paste(Outcome[i]," ~ ", xmoa))
        c = lm(fmla2, Data)}; if (sum(c$residuals)>0) {c <- stepAIC(c, k=1.2)} else {
          äm=names(cbm)[ok]
          xmo=c()
          for (ih in which(colnames(Data) %in% äm)) {xmo=append(xmo,paste("Data[, ",ih , "]", sep=""))}
          
          xmoa=paste(xmo, collapse= "+");
          fmla2 <- as.formula(paste(Outcome[i]," ~ ", xmoa))
          c = lm(fmla2, Data)
          
        }}; #https://www.scribbr.com/statistics/akaike-information-criterion/
    
    if (sum(is.na(c$coefficients))>0) {
      
      xmo=xmo[!is.na(c$coefficients[2:l])]
      xmoa=paste(xmo, collapse= "+");
      fmla2 <- as.formula(paste(Outcome[i]," ~ ", xmoa))
      c = lm(fmla2, Data)}
    
    caux=c

    #ok...
    
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    be = gfg[gfg %in% colnames(X)] #gfg #as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<(length(Treatment) +1)] #the contaminants
    ce=gfg[gfg %in% colnames(M)]
    qqq='AGE';qq='BMI';q='Gender'
    
    # if (Group=='All') {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")} else 
    # {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")}  
    # fmla2 <- as.formula(paste(Outcome[i]," ~ ", xm))
    # c = lm(fmla2, Data)
    
    #not ok...
    
    if (Group=='All') {xm=c(gfg,qqq,qq,q)} else {xm=c(gfg,qqq,qq)}
    xemm=paste(xm,collapse= "+"); fmla2 <- as.formula(paste(Outcome[i]," ~ ", xemm))

    c = lm(fmla2, Data)
    
    # Y_model_without_adjustement
    if (sum(is.na(summary(c)$coefficients))>1 | sum(is.na(summary(c)))>1) {c=caux;xm=xmo;xemm=xmoa;print('colinearity in adjustement in Y')}
    
    
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l] # 
    
    be = gfg[gfg %in% colnames(X)] #gfg #as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<(length(Treatment) +1)] #the contaminants
    # print(be)
    ce=gfg[gfg %in% colnames(M)] #colnames(M) Mediator


    joi2=dim(M)[2]; j=1
    
    for (j in 1:(joi2)) {
      # j=15
      print(i);print(j)
      print(Outcome[i]);colnames(M)[j]
      des=c();ci=c()
      
    
      
      # if the same m is not in indirect as it is in the direct (ce)
    if (sum( ce %in% colnames(M)[j])<1) 
      {des='putative_m'; xemm=paste(c(paste(c(colnames(M)[j]), collapse= "+"),paste(c(xm),collapse= "+")),collapse= "+")
      fmla2 <- as.formula(paste(Outcome[i]," ~ ", xemm)); ci = lm(fmla2, Data)} else {des='ok';ci=c}
      
    if (sum(is.na(summary(ci)$coefficients))>1 | sum(is.na(summary(ci)))>1) {ci=caux;xm=xmo;xemm=xmoa;print('colinearity in putative Y')}
      
    # ce %in% colnames(M)[j]
   
    å=joi2 #length(M); #list=1:å; list1=list[1:å!=j];list2=list1+joi2;
    MM <- data.frame(tv_red[,Mediator[j]]); 
    colnames(MM)=colnames(Data2)[8:27][j] #-? hööty 9324
    M2 <- -tv_red[,Mediator[1:å!=j]]; 
    colnames(M2)=colnames(data.frame(tv_red[,Mediator[1:å!=j]])) #Be awhare of the names...
    
    yy=as.matrix(MM)#Data2[,colnames(data.frame(M))[j]]) #rowMedians(as.matrix(Data[,colnames(data.frame(M))]))
    xx=as.matrix(cbind(X,M2))
    
    # cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); 
    try({cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1)},  {print('glmnet_na');next}); 
    
    df <- data.frame(x = seq(1:length(cv_model$lambda)),y = cv_model$lambda)#and also: hei=elbow(data = df)
    yn=find_curve_elbow(df)#or: which(hei$data[,4]==min(hei$data[,4]))[1]
    df2=df[yn:dim(df)[1],]#hei$data[yn:dim(hei$data)[1],4]# df2 <- data.frame(x = seq(1:length(df2)),y = df2)
    yn2=find_curve_elbow(df2) #plot_curve = FALSE
    yti=yn+ceiling(yn2)
    jeps=cv_model$lambda[yti]
    # best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps)
    try({best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps)},  {print('glmnet_na');next}); 
    
    
    if (sum(abs(coef(best_model))>0.0001)>1) {
      cbm=coef(best_model)[2:length(coef(best_model)),]
      ok=abs(cbm)!=0; äm=names(cbm)[ok]; xe=c()
      for (il in which(colnames(Data2) %in% äm)) {xe=append(xe,paste("Data2[, ",il , "]", sep=""))}
      xe=colnames(Data2)[colnames(Data2) %in% äm]
      fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(xe, collapse= "+")))))
      b = lm(fmla1, Data2)} else {
        x=colnames(X)
        fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(x, collapse= "+")))))
        b = lm(fmla1, Data2)
        b <- stepAIC(b, k=1.2)} #https://bookdown.org/egarpor/PM-UC3M/lm-ii-modsel.html
    
    l=length(names(b$coefficients));gfg=names(b$coefficients)[2:l]
    bex=b
    qqq='AGE';qq='BMI';q='Gender'

    if (length(gfg)>0) {
      if (Group=='All') 
      { cc=paste(c(paste(c(gfg), collapse= "+"),paste(c(qqq,qq,qq),collapse= "+")),collapse= "+")#}  
    fmla1 <- as.formula(paste(colnames(MM)," ~ ", cc))
    b = lm(fmla1, data=Data2); if (sum(is.na(b$coefficients))>0) {b=bex}} else
    {cc=paste(c(paste(c(gfg), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")#}  
    fmla1 <- as.formula(paste(colnames(MM)," ~ ", cc))
    b = lm(fmla1, data=Data2); if (sum(is.na(b$coefficients))>0) {b=bex}}} #ei datasta
    
    if (length(gfg)==0) {#b = glm(fmla1, Data2,family=gaussian)
        if (Group=='All') 
          { fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(colnames(X), 
                                  collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else
 {fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(colnames(X), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")))};
        b = lm(fmla1, Data2); 
        if (sum(is.na(b$coefficients))>0) {b=bex}; #b = glm(fmla1, Data2,family=gaussian)
        b <- stepAIC(b, k=0.4)}# bvif = append(bvif,vif(b))


    l=length(names(b$coefficients))
    gfg=names(b$coefficients)[2:l]
    gfg=gfg[!gfg %in% c('AGE','BMI','Gender','(Intercept)')]
    
    
    uip=rownames(summary(b)$coefficients)[!rownames(summary(b)$coefficients) %in% c('AGE','BMI','Gender','(Intercept)')]
    uip2=rownames(summary(ci)$coefficients)[!rownames(summary(c)$coefficients) %in% c('AGE','BMI','Gender','(Intercept)')]

    jii=intersect(uip,uip2)
    jx=colnames(X)[colnames(X) %in% jii]
    jm=colnames(M)[colnames(M) %in% jii]
    print(jm);print(jx)

    #For mediators (in M):
    ez=summary(b)$coefficients;ez= ez[rownames(ez) %in% colnames(M),];ez=ez[order(ez[,4]),]; ez2=rownames(ez)
    #For contaminants (in M):
    ezz=summary(b)$coefficients;ezz=ezz[order(ezz[,4]),];  ey2=rownames(ezz)[which(rownames(ezz) %in% colnames(X))][1]
    #For contaminants (in Y):
    ezz2=summary(ci)$coefficients;ezz2=ezz2[order(ezz2[,4]),]; eyh2=rownames(ezz2)[which(rownames(ezz2) %in% colnames(X))][1]
    if (is.na(eyh2)) {eyh2=names(colSums(X)[rev(order(colSums(X)))][1])}
    #For mediators (in Y):
    ezz3=summary(ci)$coefficients;ezz3=ezz3[order(ezz2[,4]),]; eyh3=rownames(ezz3)[which(rownames(ezz3) %in% colnames(M))][1]
    if (is.na(eyh3)) {eyh3=names(colSums(M)[rev(order(colSums(M)))][1])}
    
    # jx;jm
    
    # M_model_with_forced_int_x
    if (length(ey2)<1 | is.na(ey2)) {ey2=eyh2;if (Group=='All') 
    { fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", 
                                paste(c(paste(gfg, collapse= "+"),paste(c(ey2), collapse= "+"), paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else
    {fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c(ey2), collapse= "+"),
                                                                                 paste(c(qqq,qq),collapse= "+")),collapse= "+")))};
    b = lm(fmla1, Data2); des=paste(c('M_model_with_forced_int_x',des),collapse= "_");jx=c(ey2);
    print('Contaminants not in M model') } 
    
    #M_model_with_forced_int_m
    if (length(ez2)<1 ) {ey2=eyh3;if (Group=='All')
    { fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c(colnames(M)[!colnames(M) %in% colnames(MM)]), collapse= "+"), paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else
    {fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),
                              paste(c(colnames(M)[!colnames(M) %in% colnames(MM)]), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")))};
    b = lm(fmla1, Data2);des=paste(c('M_model_with_forced_int_m',des),collapse= "_"); jm=c(ey2);
    print('Mediators not in M model') }
    
    #For mediators (in M):
    ez=summary(b)$coefficients;ez= ez[rownames(ez) %in% colnames(M),];ez=ez[order(ez[,4]),]; ez2=rownames(ez)
    #For contaminants (in M):
    ezz=summary(b)$coefficients;ezz=ezz[order(ezz[,4]),];  ey2=rownames(ezz)[which(rownames(ezz) %in% colnames(X))][1]
    
    #Y_model_with_forced_int_m_and_x
    if (length(jm)<1 & length(jx)<1) {ci=c();
    xemm=paste(c(paste(c(ez2[1]), collapse= "+"),paste(c(ey2), collapse= "+"),paste(c(xm),collapse= "+")),collapse= "+");
    fmla2 <- as.formula(paste(Outcome[i]," ~ ", xemm)); 
    ci = lm(fmla2, Data); des=paste(c('Y_model_with_forced_int_m_and_x',des),collapse= "_");
    jm=c(ez2[1]);jx=c(ey2);
    print('Mediators and contaminants between models not intersecting')}
    
    #Y_model_with_forced_int_m
    if (length(jm)<1) {ci=c();xemm=paste(c(paste(c(ez2[1]), collapse= "+"),paste(c(xm),collapse= "+")),collapse= "+");#xemm=paste(c(ea[1],xm,collapse= "+"))
    fmla2 <- as.formula(paste(Outcome[i]," ~ ", xemm)); ci = lm(fmla2, Data); des=paste(c('Y_model_with_forced_int_m',des),collapse= "_");jm=c(ez2[1]);
    print('Mediators between models not intersecting')}
    
    #Y_model_with_forced_int_x
    if (length(jx)<1) {ci=c();xemm=paste(c(paste(c(ey2), collapse= "+"),paste(c(xm),collapse= "+")),collapse= "+");#xemm=paste(c(ea[1],xm,collapse= "+"))
    fmla2 <- as.formula(paste(Outcome[i]," ~ ", xemm)); ci = lm(fmla2, Data);
    des=paste(c('Y_model_with_forced_int_x',des),collapse= "_"); jx=c(ey2); print('Contaminants between models not intersect')}
    
    # print(jm);print(jx)


    jo=1;
    zz=1
    
    if (sum(is.na(summary(ci)$coefficients))>1 | sum(is.na(summary(ci)))>1) {next}
    if (sum(is.na(summary(b)$coefficients))>1 | sum(is.na(summary(b)))>1) {next}
    
    try({
      for (jo in 1:(length(jm))) { #control.value=mina[i]
        for (zz in 1:(length(jx))) {
          
          med_oute=mediate(b, ci, treat =  jx[zz], mediator = jm[jo], #boot=TRUE, boot.ci.type='bca',
                      sims = simss,control.value=control.value[jx[zz]],treat.value=treat.value[jx[zz]]) #[[1]]:set ei auta tässä.. eikä , control=control[,jx[zz]]
          #} else {next} control.value=control.value[be[z]],treat.value=treat.value[be[z]]
          med_out = summary(med_oute) #youm need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
          tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],
                med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],
                med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
          res <- rbind(res,tmp);
          # print(Outcome[i])
          rn=append(rn,paste(jx[zz],jm[jo],Outcome[i], des, sep=" ")) #attaching two rownames...
          remove(tmp) }}} ,  {print('mediation_notOK')}) #https://intro2r.com/loops.html 
    #https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html
    } 
    }
  #pilkku puuttu trysta: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/try
  
  # rasaa=res
  
  # try({res},{res=t(data.frame(rep(NA,16)))})
  # if (length(res)==0)  {res=t(data.frame(rep(NA,16)))}
  try({rownames(res)=rn}, {print('rn');next}) #write.csv(rn,'iii.csv')
  try({colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u',
                       'Proportion Mediated', 'n1.p','n1.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  }, {print('col');next})
  try({res=res[rev(order(res[,1])),]},{print('order');next}) #res=res[rev(order(res[,1])),]
  try({rownames(res) <- gsub("X11", "11", rownames(res))}, {print('xx');next})
  try({rownames(res) <- gsub("X17", "17", rownames(res)) }, {print('xxx');next})
    try({rownames(res) <- gsub("11b.OHA4", "11b-OHA4", rownames(res)) }, {print('xxx');next})
    try({rownames(res) <- gsub("17aOH.P4", "17a-OHP4", rownames(res)) }, {print('xxx');next})
    try({rownames(res) <- gsub("17a.OHP4", "17a-OHP4", rownames(res)) }, {print('xxx');next})
    try({rownames(res) <- gsub("17a.OHP5", "17a-OHP5", rownames(res)) }, {print('xxx');next})
    try({rownames(res) <- gsub("11.KT", "11-KT", rownames(res)) }, {print('xxx');next})
    try({rownames(res) <- gsub("11.KA4", "11-KA4", rownames(res)) }, {print('xxx');next})
    try({rownames(res) <- gsub("T.Epi.T", "T/Epi-T", rownames(res)) }, {print('xxx');next})
  # try({write.xlsx(res, file = paste(name,Group,date,'asdf.xlsx'), append = FALSE, row.names = TRUE)},{print('xls')})
  
    # write.csv(res,paste(name,Group,date,'hyp4a_oki.csv'))
    
  tryCatch({write.xlsx(res, file = paste(name,sick,Group,date,'hyp4a_ok.xlsx'), append = FALSE, row.names = TRUE) }, error = function(msg){
    return(write.csv(res,paste(name,sick,Group,date,'test.csv')))})
  
  return(res)}

# fem_ns=res
# mal_ns=res

# uh7f=fem_ns
# uh7ma=tot_ns


#Tässä tää unelma funktio... hypo4, le uusix äk:
loop_med_simplified_hyp4bb=function(Group,name,date,simss,sick,sick_group,joo,ip) { 
  # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  # muuten hyvä, mutta _S:t katoaa... #nyt saattanee toimia.. hiotaan vielä funktiota, tikka1124.. lienee ok
  # library(pathviewr) #Group='Male'
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if 
  (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])} 
  tv_red=c(); 
  # !as.vector(sick_group)
  if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else   {
    if (joo =='joo') {tv_red=tv_all} else {tv_red=tv_all[cond & !as.vector(sick_group),]}} #tv_red=tv_all[cond & !as.vector(sick_group),]} # tv_red=tv_all, I wonder this healthy should be 'cond' non-sick rather than all...
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- ip*tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()
  # simss=3; length(y)*length(m)*length(x)
  ih=1;il=1;i=1; j=1;z=1 #only for outcome model (c) put the AIC
  
  # i=14
  
  for (i in 1:length(y)) {# 
    # print(Outcome[i])
    Data=c()
    if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if
    (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
    tv_red=c();
    if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else  {tv_red=tv_all} #{tv_red=tv_all[cond,]} # tv_red=tv_all
    X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
    M <- ip*tv_red[,Mediator]  #
    
    colnames(M) <- gsub("11", "X11", colnames(M));
    colnames(M) <- gsub("17", "X17", colnames(M));
    colnames(M) <- gsub("-", ".", colnames(M));
    colnames(M) <- gsub("/", ".", colnames(M));
    joi=dim(M)[2]+dim(X)[2]
    å=length(y); list=1:å; list1=list[1:å!=i];list2=list1+joi;
    
    Y <- data.frame(tv_red[,Outcome[i]]); colnames(Y)=Outcome[i]  #-?
    Y2 <- -tv_red[,Outcome[1:å!=i]]; Y3 <- tv_red[,Outcome[1:å!=i]]    #"Steatosis.Grade.0.To.3","Fibrosis.Stage.0.to.4","Necroinflammation","HOMA-IR"
    cova <- tv_red[,c('AGE','BMI','Gender')]
    Data <- cbind(X,M,Y,Y2,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
    Data=data.frame(Data);
    colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
    
    Data2 <- cbind(X,M,Y,Y3,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
    Data2=data.frame(Data2);
    colnames(Data2) <- gsub(" ", "_", colnames(Data2)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
    
    yy=as.matrix(Data[,Outcome[i]])
    xx=as.matrix(Data[,c(colnames(data.frame(X)),colnames(data.frame(M)),colnames(data.frame(Y2)))]) #this is ok!
    
    # cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); 
    try({cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1)} ,  {print('glmnet_na');next}) 
    df <- data.frame(x = seq(1:length(cv_model$lambda)),y = cv_model$lambda) #and also: hei=elbow(data = df)
    yn=find_curve_elbow(df)#or: which(hei$data[,4]==min(hei$data[,4]))[1]
    df2=df[yn:dim(df)[1],]#hei$data[yn:dim(hei$data)[1],4]# df2 <- data.frame(x = seq(1:length(df2)),y = df2)
    yn2=find_curve_elbow(df2) #plot_curve = FALSE
    yti=yn*0.95
    jeps=cv_model$lambda[yti] #coef(best_model)[abs(coef(best_model))>0.00001]
    # best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps) 
    try({best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps) } ,  {print('glmnet_na');next}) 
    cbm=coef(best_model)[2:length(coef(best_model)),]
    ok=abs(cbm)!=0
    äm=names(cbm)[ok]
    xm=c(); for (ih in which(colnames(Data) %in% äm)) {xm=append(xm,paste("Data[, ",ih , "]", sep=""))}
    xm=colnames(Data)[colnames(Data) %in% äm]
    if (sum(colnames(X) %in% xm)==0) {xm=c(xm,colnames(X))}
    # if (length(c$coefficients)>0) {xm=xm[!is.na(c$coefficients[2:l])]} else {xm=xm}
    xma=paste(xm, collapse= "+");
    fmla2 <- as.formula(paste(Outcome[i]," ~ ", xma))
    c = lm(fmla2, Data)
    l=length(names(c$coefficients))
    if (l<3) {
      xm=paste(c(x,m,y[list1]), collapse= "+"); fmla2 <- as.formula(paste(Outcome[i]," ~ ", xm)); 
      c = lm(fmla2, Data);
      if (sum(is.na(c$coefficients))>0) {
        l=length(names(c$coefficients))
        xm=c(x,m,y[list1])[!is.na(c$coefficients[2:l])]
        xma=paste(xm, collapse= "+");
        fmla2 <- as.formula(paste(Outcome[i]," ~ ", xma))
        c = lm(fmla2, Data)}; if (sum(c$residuals)>0) {c <- stepAIC(c, k=1.2)} else {
          äm=names(cbm)[ok]
          xm=c()
          for (ih in which(colnames(Data) %in% äm)) {xm=append(xm,paste("Data[, ",ih , "]", sep=""))}
          
          xma=paste(xm, collapse= "+");
          fmla2 <- as.formula(paste(Outcome[i]," ~ ", xma))
          c = lm(fmla2, Data)
          
        }}; #https://www.scribbr.com/statistics/akaike-information-criterion/
    
    if (sum(is.na(c$coefficients))>0) {
      
      xm=xm[!is.na(c$coefficients[2:l])]
      xma=paste(xm, collapse= "+");
      fmla2 <- as.formula(paste(Outcome[i]," ~ ", xma))
      c = lm(fmla2, Data)}
    # l=length(names(c$coefficients))
    # gfg=names(c$coefficients)[2:l]
    # gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    # vex=as.numeric(unlist(gfg_numbers))[as.vector(is.na(c$coefficients[2:l]))]
    
    # c$coefficients=c$coefficients[!is.na(c$coefficients)]
    
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    # gfg_numbers <- regmatches(Treatment, gregexpr("[[:digit:]]+", gfg))
    be = gfg[gfg %in% colnames(X)] #gfg #as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<(length(Treatment) +1)] #the contaminants
    # ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>length(Treatment) ];
    # ce=ce1[ce1<(length(c(colnames(X),colnames(M)))+1)]; ce=ce-length(Treatment) #the mediators
    # print(gfg)
    ce=gfg[gfg %in% colnames(M)]
    # print(ce)
    
    # qqq=paste(c("Data[, ",as.numeric(dim(Data)[2])-3,"]"),collapse="")
    qqq='AGE';qq='BMI';q='Gender'
    # qq=paste(c("Data[, ",as.numeric(dim(Data)[2])-2,"]"),collapse="")
    # q=paste(c("Data[, ",as.numeric(dim(Data)[2])-1,"]"),collapse="")
    if (Group=='All') {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")} else 
    {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")}  
    fmla2 <- as.formula(paste(Outcome[i]," ~ ", xm))
    c = lm(fmla2, Data)
    
    if (Group=='All') {xm=c(gfg,qqq,qq,q)} else
    {xm=c(gfg,qqq,qq)}
    # fmla2 <- as.formula(paste(Outcome[i]," ~ ", xm))
    
    # oy=as.numeric(regmatches(Outcome[i], gregexpr("[[:digit:]]+", Outcome[i])));yh=paste(colnames(Data)[oy],collapse= "+")
    # xmm=as.numeric(regmatches(xm, gregexpr("[[:digit:]]+", xm)));
    # xemm=paste(colnames(Data)[xmm],collapse= "+")
    xemm=paste(xm,collapse= "+")
    fmla2 <- as.formula(paste(Outcome[i]," ~ ", xemm))
    
    # c = gam(fmla2, Data,family=gaussian) #gam??? :)
    c = lm(fmla2, Data)
    
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l] # 
    # print(gfg)
    # gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    # be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<(length(Treatment) +1)] #the contaminants
    # ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>length(Treatment) ];
    # ce=ce1[ce1<(length(c(colnames(X),colnames(M)))+1)]; ce=ce-length(Treatment) #the mediators
    be = gfg[gfg %in% colnames(X)] #gfg #as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<(length(Treatment) +1)] #the contaminants
    # print(be)
    ce=gfg[gfg %in% colnames(M)] #colnames(M) Mediator
    # print(ce)
    
    # }
    
    joi2=dim(M)[2]; j=1
    
    for (j in 1:(joi2)) {
      des=c()
      if (sum( ce %in% colnames(M)[j])<1) {des='putative'} else {des='ok'}
      
      # ce %in% colnames(M)[j]
      
      å=joi2 #length(M); #list=1:å; list1=list[1:å!=j];list2=list1+joi2;
      MM <- data.frame(tv_red[,Mediator[j]]); 
      colnames(MM)=colnames(Data2)[8:27][j] #-? hööty 9324
      M2 <- -tv_red[,Mediator[1:å!=j]]; 
      colnames(M2)=colnames(data.frame(tv_red[,Mediator[1:å!=j]])) #Be awhare of the names...
      
      yy=as.matrix(MM)#Data2[,colnames(data.frame(M))[j]]) #rowMedians(as.matrix(Data[,colnames(data.frame(M))]))
      xx=as.matrix(cbind(X,M2))
      
      # cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); 
      try({cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1)},  {print('glmnet_na');next}); 
      
      df <- data.frame(x = seq(1:length(cv_model$lambda)),y = cv_model$lambda)#and also: hei=elbow(data = df)
      yn=find_curve_elbow(df)#or: which(hei$data[,4]==min(hei$data[,4]))[1]
      df2=df[yn:dim(df)[1],]#hei$data[yn:dim(hei$data)[1],4]# df2 <- data.frame(x = seq(1:length(df2)),y = df2)
      yn2=find_curve_elbow(df2) #plot_curve = FALSE
      yti=yn+ceiling(yn2)
      jeps=cv_model$lambda[yti]
      # best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps)
      try({best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps)},  {print('glmnet_na');next}); 
      
      
      if (sum(abs(coef(best_model))>0.0001)>1) {
        cbm=coef(best_model)[2:length(coef(best_model)),]
        ok=abs(cbm)!=0; äm=names(cbm)[ok]; xe=c()
        for (il in which(colnames(Data2) %in% äm)) {xe=append(xe,paste("Data2[, ",il , "]", sep=""))}
        xe=colnames(Data2)[colnames(Data2) %in% äm]
        fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(xe, collapse= "+")))))
        b = lm(fmla1, Data2)} else {
          x=colnames(X)
          fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(x, collapse= "+")))))
          b = lm(fmla1, Data2)
          b <- stepAIC(b, k=1.2)} #https://bookdown.org/egarpor/PM-UC3M/lm-ii-modsel.html
      
      l=length(names(b$coefficients))
      gfg=names(b$coefficients)[2:l]
      
      # print(gfg)
      
      qqq='AGE';qq='BMI';q='Gender'
      
      # 
      # uip=rownames(summary(b)$coefficients)[!rownames(summary(b)$coefficients) %in% c('AGE','BMI','Gender','(Intercept)')]
      # try({uip2=names(summary(c)[[1]])[!names(summary(c)[[1]]) %in% c('AGE','BMI','Gender','(Intercept)')]},
      #     {uip2=names(c$coefficients)[!names(c$coefficients) %in% c('AGE','BMI','Gender','(Intercept)')]})
      # 
      # jii=intersect(uip,uip2)
      # jx=colnames(X)[colnames(X) %in% jii]
      # jm=colnames(M)[colnames(M) %in% jii]
      # 
      # #be #the contaminants
      # #ce #the mediators
      # 
      # if (length(jx)<1) {print('Contaminants_na')}#jx=names(c$coefficients)[names(c$coefficients) %in% colnames(X)]
      # if (length(jm)<1) {print('Mediators_na')} #jm=names(b$coefficients)[names(b$coefficients) %in% colnames(M)]
      
      # print(gfg)
      if (length(gfg)>0) {
        if (Group=='All') 
        { fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),
                                                                                      paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else
                                                                                      {fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", 
                                                                                                                 paste(c(paste(gfg, collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")))};
        b = lm(fmla1, Data)} else {#b = glm(fmla1, Data2,family=gaussian)
          if (Group=='All') { fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(x, collapse= "+"),
                                                                                                          paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else
                                                                                                          {fmla1 <- as.formula(paste(paste(colnames(MM), collapse= "+")," ~ ", paste(c(paste(x, collapse= "+"),
                                                                                                                                                                                       paste(c(qqq,qq),collapse= "+")),collapse= "+")))};
          b = lm(fmla1, Data) ; #b = glm(fmla1, Data2,family=gaussian)
          b <- stepAIC(b, k=0.4)}# bvif = append(bvif,vif(b))
      
      
      # lu=length(names(b$coefficients))
      # gfgu=names(b$coefficients)[2:lu]
      # #
      # beu = gfgu[gfgu %in% colnames(X)]  #the contaminants
      # print(be);print(beu)
      # ceu = gfgu[gfgu %in% colnames(M)] #the mediators
      # print(ce);print(ceu)
      # 
      # # gfgu_numbers <- regmatches(gfgu, gregexpr("[[:digit:]]+", gfgu))
      # # beu=as.numeric(unlist(gfgu_numbers))[as.numeric(unlist(gfgu_numbers))<(length(Treatment) +1)]
      # 
      # be=intersect(be,beu) #contaminants
      # ce=intersect(ce,ceu) #mediators
      # 
      # jx=be
      # jm=ce
      # 
      # if (length(jx)<1) {print('Contaminants_na_intersect')}
      # if (length(jm)<1) {print('Mediators_na_intersect')} 
      
      # 
      uip=rownames(summary(b)$coefficients)[!rownames(summary(b)$coefficients) %in% c('AGE','BMI','Gender','(Intercept)')]
      uip2=rownames(summary(c)$coefficients)[!rownames(summary(c)$coefficients) %in% c('AGE','BMI','Gender','(Intercept)')]
      # try({uip2=names(summary(c)[[1]])[!names(summary(c)[[1]]) %in% c('AGE','BMI','Gender','(Intercept)')]},
      #     {uip2=names(c$coefficients)[!names(c$coefficients) %in% c('AGE','BMI','Gender','(Intercept)')]})
      
      # print(uip)
      # print(uip[uip %in% colnames(M)])
      
      jii=intersect(uip,uip2)
      jx=colnames(X)[colnames(X) %in% jii]
      jm=colnames(M)[colnames(M) %in% jii]
      
      if (length(jx)<1) {print('Contaminants between models not intersect')}
      if (length(jm)<1) {print('Mediators between models not intersect')}
      print(jm)
      
      jo=1;
      zz=1
      
      try({
        for (jo in 1:length(jm)) { #control.value=mina[i]
          for (zz in 1:length(jx)) {
            
            med_oute=mediate(b, c, treat =  jx[zz], mediator = jm[jo], #boot=TRUE, boot.ci.type='bca',
                             sims = simss,control.value=control.value[jx[zz]],treat.value=treat.value[jx[zz]]) 
            #} else {next} control.value=control.value[be[z]],treat.value=treat.value[be[z]]
            med_out = summary(med_oute) #youm need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
            tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],
                  med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],
                  med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
            res <- rbind(res,tmp);
            # print(Outcome[i])
            rn=append(rn,paste(jx[zz],jm[jo],Outcome[i], des, sep=" ")) #attaching two rownames...
            remove(tmp) }}} ,  {print('mediate_na');next}) #https://intro2r.com/loops.html 
      #https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html
    } 
  }
  #pilkku puuttu trysta: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/try
  
  # rasaa=res
  
  # try({res},{res=t(data.frame(rep(NA,16)))})
  # if (length(res)==0)  {res=t(data.frame(rep(NA,16)))}
  try({rownames(res)=rn}, {print('rn');next}) #write.csv(rn,'iii.csv')
  try({colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u',
                       'Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  }, {print('col');next})
  try({res=res[rev(order(res[,1])),]},{print('order');next}) #res=res[rev(order(res[,1])),]
  try({rownames(res) <- gsub("X11", "11", rownames(res))}, {print('xx');next})
  try({rownames(res) <- gsub("X17", "17", rownames(res)) }, {print('xxx');next})
  # try({write.xlsx(res, file = paste(name,Group,date,'asdf.xlsx'), append = FALSE, row.names = TRUE)},{print('xls')})
  
  tryCatch({write.xlsx(res, file = paste(name,Group,date,'fun_hyp4bb_jeah.xlsx'), append = FALSE, row.names = TRUE) }, error = function(msg){
    return(write.csv(res,paste(name,Group,date,'testa2.csv')))})
  
  return(res)}


# Group,name,date,simss,sick,sick_group,joo,ip
loop_med_simplified_hyp4=function(Treatment,Mediator,Outcome,tv_all,Group,name,date,simss,t.val,test,sick,sick_group,joo) { 
  # if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') # {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  # muuten hyvä, mutta _S:t katoaa... #nyt saattanee toimia.. hiotaan vielä funktiota, tikka1124.. lienee ok
  # library(pathviewr) #Group='Male'
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if 
  (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c(); 
  # !as.vector(sick_group)
  if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else   {
    if (joo =='joo') {tv_red=tv_all} else {tv_red=tv_all[cond & !as.vector(sick_group),]}} #tv_red=tv_all[cond & !as.vector(sick_group),]} # tv_red=tv_all, I wonder this healthy should be 'cond' non-sick rather than all...
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  M <- tv_red[,Mediator]  #
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  cova <- tv_red[,c('AGE','BMI','Gender')] 
  Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data)
  colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  #M~X
  x=c();m=c(); y=c();ye=c()
  for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  #Y~X+M
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2]-3)) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
  med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c()
  # simss=3; length(y)*length(m)*length(x)
  ih=1;il=1;i=1; j=1;z=1 #only for outcome model (c) put the AIC
  # i=1
  
  for (i in 1:length(y)) {# 
    print(Outcome[i])
    Data=c()
    if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if 
    (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
    tv_red=c(); 
    if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else {tv_red=tv_all}  #{tv_red=tv_all[cond,]} # tv_red=tv_all
    X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
    M <- tv_red[,Mediator]  #
    Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3","Fibrosis.Stage.0.to.4","Necroinflammation","HOMA-IR"   
    cova <- tv_red[,c('AGE','BMI','Gender')] 
    Data <- cbind(X,M,Y,cova);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
    Data=data.frame(Data); colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
    å=length(y); list=1:å; list1=list[1:å!=i]; ya=c();list2=list1+28; 
    Data[,list2]=-Data[,list2]
    yy=as.matrix(Data[,Outcome[i]])
    D1=data.frame(Data[,list2])
    Data[,(dim(Data)[2]+1)]=rowSums(D1)
    xx=as.matrix(Data[,c(colnames(data.frame(X)),colnames(data.frame(M)),colnames(Data)[dim(Data)[2]])])
    # cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1); 
    try({cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1)}{next}); 
    df <- data.frame(x = seq(1:length(cv_model$lambda)),y = cv_model$lambda)#and also: hei=elbow(data = df)
    yn=find_curve_elbow(df)#or: which(hei$data[,4]==min(hei$data[,4]))[1]
    df2=df[yn:dim(df)[1],]#hei$data[yn:dim(hei$data)[1],4]# df2 <- data.frame(x = seq(1:length(df2)),y = df2)
    yn2=find_curve_elbow(df2) #plot_curve = FALSE
    yti=yn*0.95
    jeps=cv_model$lambda[yti] #coef(best_model)[abs(coef(best_model))>0.00001]
    # best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps) 
    try({best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps)}{next}); 
    cbm=coef(best_model)[2:length(coef(best_model)),]
    ok=abs(cbm)!=0
    äm=names(cbm)[ok]
    xm=c(); for (ih in which(colnames(Data) %in% äm)) {xm=append(xm,paste("Data[, ",ih , "]", sep=""))}
    # if (length(c$coefficients)>0) {xm=xm[!is.na(c$coefficients[2:l])]} else {xm=xm}
    xma=paste(xm, collapse= "+");
    fmla2 <- as.formula(paste(y[i]," ~ ", xma))
    c = lm(fmla2, Data)
    l=length(names(c$coefficients))
    summary(c)
    if (l<3) {
      xm=paste(c(x,m,y[list1]), collapse= "+"); fmla2 <- as.formula(paste(y[i]," ~ ", xm)); c = lm(fmla2, Data);
      if (sum(is.na(c$coefficients))>0) {
        l=length(names(c$coefficients))
        xm=c(x,m,y[list1])[!is.na(c$coefficients[2:l])]
        xma=paste(xm, collapse= "+");
        fmla2 <- as.formula(paste(y[i]," ~ ", xma))
        c = lm(fmla2, Data)}; if (sum(c$residuals)>0) {c <- stepAIC(c, k=1.2)} else {
          äm=names(cbm)[ok]
          xm=c()
          for (ih in which(colnames(Data) %in% äm)) {xm=append(xm,paste("Data[, ",ih , "]", sep=""))}
          
          xma=paste(xm, collapse= "+");
          fmla2 <- as.formula(paste(y[i]," ~ ", xma))
          c = lm(fmla2, Data)
          
        }}; #https://www.scribbr.com/statistics/akaike-information-criterion/
    
    if (sum(is.na(c$coefficients))>0) {
      
      xm=xm[!is.na(c$coefficients[2:l])]
      xma=paste(xm, collapse= "+");
      fmla2 <- as.formula(paste(y[i]," ~ ", xma))
      c = lm(fmla2, Data)}
    # l=length(names(c$coefficients))
    # gfg=names(c$coefficients)[2:l]
    # gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    # vex=as.numeric(unlist(gfg_numbers))[as.vector(is.na(c$coefficients[2:l]))]
    
    # c$coefficients=c$coefficients[!is.na(c$coefficients)]
    
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9]
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8
    
    # summary(c)
    
    qqq=paste(c("Data[, ",as.numeric(dim(Data)[2])-3,"]"),collapse="")
    qq=paste(c("Data[, ",as.numeric(dim(Data)[2])-2,"]"),collapse="")
    q=paste(c("Data[, ",as.numeric(dim(Data)[2])-1,"]"),collapse="")
    if (Group=='All') {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")} else 
    {xm=paste(c(paste(c(gfg), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")} 
    fmla2 <- as.formula(paste(y[i]," ~ ", xm))
    c = lm(fmla2, Data)
    l=length(names(c$coefficients))
    gfg=names(c$coefficients)[2:l]
    # print(gfg)
    gfg_numbers <- regmatches(gfg, gregexpr("[[:digit:]]+", gfg))
    be=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))<9] #the contaminants
    ce1=as.numeric(unlist(gfg_numbers))[as.numeric(unlist(gfg_numbers))>8];
    ce=ce1[ce1<29]; ce=ce-8 #the mediators
    
    summary(c)
    
    yy=as.matrix(rowMedians(as.matrix(Data[,colnames(data.frame(M))]))) #needs to be as.matrix...
    xx=as.matrix(Data[,colnames(data.frame(X))])
    try({cv_model <- cv.glmnet(x=xx, y=yy, alpha = 1)}{next}); 
    df <- data.frame(x = seq(1:length(cv_model$lambda)),y = cv_model$lambda)#and also: hei=elbow(data = df)
    yn=find_curve_elbow(df)#or: which(hei$data[,4]==min(hei$data[,4]))[1]
    df2=df[yn:dim(df)[1],]#hei$data[yn:dim(hei$data)[1],4]# df2 <- data.frame(x = seq(1:length(df2)),y = df2)
    yn2=find_curve_elbow(df2) #plot_curve = FALSE
    yti=yn+ceiling(yn2)
    jeps=cv_model$lambda[yti]
    best_model <- glmnet(xx, yy, alpha = 1, lambda = jeps)
    
    if (sum(abs(coef(best_model))>0.0001)>1) {
      cbm=coef(best_model)[2:length(coef(best_model)),]
      ok=abs(cbm)!=0; äm=names(cbm)[ok]; xe=c()
      for (il in which(colnames(Data) %in% äm)) {xe=append(xe,paste("Data[, ",il , "]", sep=""))}
      fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(xe, collapse= "+")))))
      b = lm(fmla1, Data)} else {
        fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+")))))
        b = lm(fmla1, Data)
        b <- stepAIC(b, k=1.2)} #https://bookdown.org/egarpor/PM-UC3M/lm-ii-modsel.html
    
    l=length(names(b$coefficients))
    gfg=names(b$coefficients)[2:l]
    # print(gfg)
    if (length(gfg)>0) {
      if (Group=='All') 
      { fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else
      {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(gfg, collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")))};
      b = lm(fmla1, Data)} else {
        if (Group=='All') { fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+"),
                                                                                             paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else
                                                                                             {fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(c(paste(x, collapse= "+"),
                                                                                                                                                               paste(c(qqq,qq),collapse= "+")),
                                                                                                                                                             collapse= "+")))};
        b = lm(fmla1, Data);b <- stepAIC(b, k=0.4)}# bvif = append(bvif,vif(b))
    lu=length(names(b$coefficients))
    gfgu=names(b$coefficients)[2:lu]
    # print(gfgu)
    gfgu_numbers <- regmatches(gfgu, gregexpr("[[:digit:]]+", gfgu))
    beu=as.numeric(unlist(gfgu_numbers))[as.numeric(unlist(gfgu_numbers))<9]
    be=intersect(be,beu) #contaminants
    print(be);print(i)
    if (length(be)<1) {next}
    if (length(ce)<1) {next} # x=x[1:length(be)]
    
    try({
      for (j in 1:length(ce)) { #control.value=mina[i]
        for (z in 1:length(be)) {
          
          med_oute=mediate(b, c, treat =  x[be[z]], mediator = m[ce[j]],sims = simss,control.value=control.value[be[z]],treat.value=treat.value[be[z]]) #} .. check thiselse {next} control.value=control.value[be[z]],treat.value=treat.value[be[z]]
          med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
          tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],
                med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],
                med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
          # print(Outcome[i])
          rn=append(rn,paste(colnames(X)[be[z]],colnames(M)[ce[j]],Outcome[i], sep=" ")) #attaching two rownames...
          res <- rbind(res,tmp);
          remove(tmp) }} },  {print('medc');next}) #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html
  } 
  #pilkku puuttu trysta: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/try
  # rownames(res)=rn #write.csv(rn,'iii.csv')
  # colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  # res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
  # rownames(res) <- gsub("X11", "11", rownames(res))
  # rownames(res) <- gsub("X17", "17", rownames(res))
  # dim(res)
  if (length(res)==0)  {res=t(data.frame(rep(NA,16)))}
  # try({res},{res=t(data.frame(rep(0,16)))})
  try({colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
  }, {print('col')}) 
  try({rownames(res)=rn}, {rownames(res)=1}) #write.csv(rn,'iii.csv')
  try({res=res[order(res[,2]),]},{print('order')}) #res=res[rev(order(res[,1])),]
  # try({rownames(res) <- gsub("X11", "11", rownames(res))}, {print('xx')})
  tryCatch({rownames(res) <- gsub("X11", "11", rownames(res))}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # try({rownames(res) <- gsub("X17", "17", rownames(res)) }, {print('xxx')})
  tryCatch({rownames(res) <- gsub("X17", "17", rownames(res))}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # # head(res)
  # # Group='Male'
  # tryCatch({write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) }, error = function(msg){
  # return(write.csv(res,paste(name,Group,date,'.csv')))})
  try({write.xlsx(res, file = paste(name,Group,date,'hypo4.xlsx'), append = FALSE, row.names = TRUE)},{write.csv(res,paste(name,Group,date,'.csv'));next})
  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
  return(res)} #list(res,bvif,cvif)

# tryCatch({asdf}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# uh7f=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group)

lp_fn=function(tv_all,Group,sick,sick_group,fn,name) {
  
  if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if 
  (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=c();x=c();m=c(); y=c();ii=1;i=1; koko=data.frame() 
  if (sick=='yes') {tv_red=tv_all[cond & as.vector(sick_group),]} else   {tv_red=tv_all}#tv_red=tv_all[cond & !as.vector(sick_group),]} # 
  X <- tv_red[,Treatment] #Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
  control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
  treat.value=colMaxs(as.matrix(X))
  M <- tv_red[,Mediator]  ## Thinking of remove multiple BA values by list, here the liver ones seem to be also good: #https://sparkbyexamples.com/r-programming/r-remove-from-vector-with-examples/
  Y <- tv_red[,Outcome]   #"Steatosis.Grade.0.To.3","Fibrosis.Stage.0.to.4","Necroinflammation","HOMA-IR"   # cova <- tv_red[,c('AGE','BMI','Gender')] #dim(X); dim(Y); dim(M)
  cova <- tv_red[,c('AGE','BMI','Gender')] #dim(X); dim(Y); dim(M)
  Data <- cbind(X,M,Y);   ## colnames(Data)[which(names(Data) == "X")] <- "PFOA_L"
  Data=data.frame(Data); colnames(Data) <- gsub(" ", "_", colnames(Data)) # colnames(Data[,1:2])[1]=Treatment # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
  for (ii in 1:length(colnames(X))) {x=append(x,paste("Data[, ",ii , "]", sep=""))}
  for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
  for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2])) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names  
  
  for (i in 1:length(y)) { #this loop is for the lipids
    # i=1
    print(c(Outcome[i], i)) 
    pl=length(colnames(X))
    mil=length(colnames(M))
    zl=pl+length(colnames(M))
    å=length(y); list=1:å; list1=list[1:å!=i]; ya=c();list2=list1+zl;#print(list2)
    Datae=c()
    Datae=Data
    Datae[,list2]=-Datae[,list2]
    y_true=as.matrix(Datae[,Outcome[i]]) #these are the true values... sum(xtot[1,])
    y_calc=Datae[,!colnames(Datae) %in% Outcome[i]] 
    # y_calc[,1:28]=-y_calc[,1:28] #maybe easier like this...
    ## Set the coefficients of the decision variables -> C
    for_one=c();mMa_tot=c();mXa_tot=c()
    j=1;z=1;
    
    for (j in 1:dim(y_calc)[1]) { #this is for the samples/cases
      ## Set the coefficients of the decision variables -> C
      C <- as.matrix(y_calc[j,]-y_true[j]) #c(30, 40, 80)
      C=as.vector(C)
      
      #Defining the constraints:
      Xl <- matrix(rep(0,dim(y_calc)[2]), nrow=1, byrow=TRUE) #this should be ok, the expression values are '1's so you basically count the number of variables here :)
      Xu <- matrix(rep(0,dim(y_calc)[2]), nrow=1, byrow=TRUE) #if needed check: #https://stackoverflow.com/questions/55753409/implementing-additional-constraint-variables-in-integer-programming-using-lpsolv
      Ml <- matrix(rep(0,dim(y_calc)[2]), nrow=1, byrow=TRUE) # :) yesh
      Mu <- matrix(rep(0,dim(y_calc)[2]), nrow=1, byrow=TRUE)
      Yl <- matrix(rep(0,dim(y_calc)[2]), nrow=1, byrow=TRUE)
      Yu <- matrix(rep(0,dim(y_calc)[2]), nrow=1, byrow=TRUE)
      
      Xl[1:pl]=1;Xu[1:pl]=1; Ml[(pl+1):(zl)]=1;Mu[(pl+1):(zl)]=1
      Yl[(zl+1):dim(y_calc)[2]]=1 #*(1:å!=i)
      Yu[(zl+1):dim(y_calc)[2]]=1 # sum(Xl);sum(Ml);sum(Yl);dim(X)[2];dim(M)[2];dim(Y)[2]
      
      #The constraint matrix:
      A=matrix(c(Xl,Xu,Ml,Mu,Yl,Yu), nrow=6, byrow=TRUE)
      # Right hand side for the constraints
      # xd=min(y_calc[1,1:8][order(unlist(y_calc[1,1:8]))]) #this is probably ok
      # xu=y_calc[1,1:8][order(unlist(y_calc[1,1:8]))][8] #but this is probably not, so I need the constraints for the number of 
      
      B <- c(1, pl,1,mil,2,54)
      # Direction of the constraints
      constranints_direction  <- c(">=", "<=",">=","<=", ">=","<=")
      # Find the optimal solution
      optimum <-  lp(direction="min", objective.in = C, const.mat = A,const.dir = constranints_direction, const.rhs = B, all.bin = T)
      best_sol <- optimum$solution; #print(best_sol);sum(y_calc[1,as.logical(best_sol)]);colnames(Datae[1,as.logical(best_sol)])
      # => 'Y' alias c...ok :)
      
      min_M=Mediator[Mediator %in% colnames(Datae[1,as.logical(best_sol)])] #lisää nää m:n
      min_X=Treatment[Treatment %in% colnames(Datae[1,as.logical(best_sol)])] #nääki on ehkä lisättävä..
      # min_M;min_X
      
      #for the M; apparently needed for each M
      ho=1 
      for (ho in 1:mil) {
        Dataa <- cbind(X,M);   ## colnames(Dataa)[which(names(Dataa) == "X")] <- "PFOA_L"
        m_true=as.matrix(Dataa[,Mediator[ho]]) #these are the true values... sum(xtot[1,]) :)
        m_calc=Dataa[,!colnames(Dataa) %in% Mediator[ho]]
        m_calc[,(pl+1):(zl-1)]=-m_calc[,(pl+1):(zl-1)] #maybe easier like this...
        
        ## Set the coefficients of the decision variables -> C
        C <- as.matrix(m_calc[j,]-m_true[j]) #c(30, 40, 80)
        C=as.vector(C)
        
        #Defining the constraints:
        Xl <- matrix(rep(0,dim(m_calc)[2]), nrow=1, byrow=TRUE) #this should be ok, the expression values are '1's so you basically count the number of variables here :)
        Xu <- matrix(rep(0,dim(m_calc)[2]), nrow=1, byrow=TRUE) #if needed check: #https://stackoverflow.com/questions/55753409/implementing-additional-constraint-variables-in-integer-programming-using-lpsolv
        Ml <- matrix(rep(0,dim(m_calc)[2]), nrow=1, byrow=TRUE) # :) yesh
        Mu <- matrix(rep(0,dim(m_calc)[2]), nrow=1, byrow=TRUE)
        Xl[1:pl]=1;Xu[1:pl]=1; Ml[(pl+1):(zl-1)]=1;Mu[(pl+1):(zl-1)]=1
        # Xl[1:8]=1;Xu[1:8]=1; Ml[9:27]=1;Mu[9:27]=1# Yl[29:dim(m_calc)[2]]=1;Yu[29:dim(m_calc)[2]]=1 # sum(Xl);sum(Ml);sum(Yl);dim(X)[2];dim(M)[2];dim(Y)[2]
        
        #The constraint matrix:
        A=matrix(c(Xl,Xu,Ml,Mu), nrow=4, byrow=TRUE)
        # Right hand side for the constraints
        # xd=min(m_calc[1,1:8][order(unlist(m_calc[1,1:8]))]) #this is probably ok
        # xu=m_calc[1,1:8][order(unlist(m_calc[1,1:8]))][8] #but this is probably not, so I need the constraints for the number of
        
        B <- c(1, pl,1,mil-1)
        # Direction of the constraints
        constranints_direction  <- c(">=", "<=",">=","<=")
        # Find the optimal solution
        optimum <-  lp(direction="min", objective.in = C, const.mat = A,const.dir = constranints_direction, const.rhs = B, all.bin = T)
        best_sol <- optimum$solution; #print(best_sol);sum(m_calc[1,as.logical(best_sol)]);names(Dataa[1,as.logical(best_sol)])
        # colnames(Dataaa[1,as.logical(best_sol)])
        # => 'M' alias b...ok :)
        min_Ma=Mediator[Mediator %in% names(Dataa[1,as.logical(best_sol)])] #lisää nää m:n
        min_Xa=Treatment[Treatment %in% names(Dataa[1,as.logical(best_sol)])] #nääki on ehkä lisättävä..
        # min_Ma;min_Xa
        mMa_tot=append(mMa_tot,min_Ma)
        mXa_tot=append(mXa_tot,min_Xa)}
      
      sm=c();sx=c();mt=c();xt=c();
      mt=table(unlist(mMa_tot));#hist(mt,breaks=30) #all are there.... if all are taken...
      sel1=ceiling((median(mt)-sd(mt))*1.05)
      sm=mt[mt>sel1];#sm# 11-KA4  11-KDHT    11-KT 11b-OHA4 17a-OHP5 17aOH-P4        A       A4       AN        B     DHEA     ....            55       70       57       63       30 
      xt=table(unlist(mXa_tot)) ;#hist(xt,breaks=30)
      sel2=ceiling((median(xt)-sd(xt))*1.05)
      sx=xt[xt>sel2];#sx
      
      min_Maa=names(sm)
      min_Xaa=names(sx)
      # Mt=list(c(min_M,min_Ma)) # Xt=list(c(min_X,min_Xa)) # fo=c(Mt,Xt)
      
      # Intersections
      M_is=list(intersect(min_M,min_Maa));#M_is;min_M #ok
      X_is=list(intersect(min_X,min_Xaa));#X_is;min_X #ok :)
      fo2=c(M_is,X_is)
      # XM=c(M_is,X_is)
      for_one=append(for_one,fo2)
    } # }...this not here
    
    # for_one # declaring a vector # vec <- 5:20 #print("Extracting every 4th element")
    # extracted_vec <- vec[seq(1, length(vec), 4)] #https://www.geeksforgeeks.org/extract-every-nth-element-of-a-vector-in-r/
    for_one_m <- for_one[seq(1, length(for_one)-1, 2)]
    for_one_x <- for_one[seq(2, length(for_one), 2)]
    # for_one_m_is=intersect(unlist(for_one_m))
    for_one_m_is=unique(unlist(for_one_m))
    sm=c();sx=c();mt=c();xt=c();
    mt=table(unlist(for_one_m));#hist(mt,breaks=30) #all are there.... if all are taken...
    # sel1=ceiling((median(mt)-sd(mt))*1.05) #sel1=ceiling((median(mt)))
    if (Group=='All') {sel1=ceiling((median(mt)-sd(mt))*1.05)} else {sel1=ceiling(median(mt))}
    sm=mt[mt>sel1];sm# 11-KA4  11-KDHT    11-KT 11b-OHA4 17a-O     P4  ....          57       68       69       65       55       70       57       63       30 
    if (sum(mt>sel1)==0) {sel1=ceiling(median(mt)-1)};sm=mt[mt>sel1];sm
    for_one_x_is=unique(unlist(for_one_x))
    xt=table(unlist(for_one_x)) ;#hist(xt,breaks=30)
    if (Group=='All') {sel2=ceiling((median(xt)-sd(xt))*1.05)} else {sel2=ceiling(median(xt)+1)} 
    sx=xt[xt>sel2];sx
    if (sum(xt>sel2)==0) {sel2=ceiling(median(xt)-2)}; sx=xt[xt>sel2];sx
    
    Data_v= cbind(Datae,cova);
    qqq=colnames(cova)[1]#paste(c("Data[, ",as.numeric(dim(Data)[2])-2,"]"),collapse="")
    qq=colnames(cova)[2]##paste(c("Data[, ",as.numeric(dim(Data)[2])-1,"]"),collapse="")
    q=colnames(cova)[3]##paste(c("Data[, ",as.numeric(dim(Data)[2]),"]"),collapse="")
    XM=c(names(sm),names(sx))# XM=c(names(mt),names(xt))
    if (Group=='All') {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")} else 
    {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")} 
    fmla2 <- as.formula(paste(Outcome[i]," ~ ", xm));fmla2
    
    c = lm(fmla2, Data_v);
    # print(summary(c))
    
    if (sum(is.na(summary(c)[[4]][,2]))>0) {sel2=ceiling(median(xt)+1);sx=xt[xt>sel2]; 
    XM=c(names(sm),names(xt[rev(order(xt))[1:2]]))# XM=c(names(mt),names(xt))
    if (Group=='All') {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")} else 
    {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")} 
    fmla2 <- as.formula(paste(Outcome[i]," ~ ", xm));fmla2;c = lm(fmla2, Data_v)}
    
    if (sum(is.na(summary(c)[[4]][,2]))>0) {
      sela=c$coefficients[!is.na(c$coefficients)]
      sela2=sela[names(sela) %in% colnames(M)]
      sv1=names(sela2[order(sela2)][1:2]);sv2=names(sela2[rev(order(sela2))][1:2])
      jap=c(sv1,sv1)
      sela3=sela[names(sela) %in% colnames(X)]
      sv1=names(sela3[order(sela3)][1]);sv2=names(sela3[rev(order(sela3))][1])
      jap1=c(sv1,sv2); XM=c(jap,jap1);XM=XM[!is.na(XM)]
      if (sum(XM %in% colnames(X))==0) {sv2=names(xt[rev(order(xt))[1]]);jap1=c(sv1,sv2); XM=c(jap,jap1);XM=XM[!is.na(XM)]}
      if (Group=='All') {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")} else 
      {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")} 
      fmla2 <- as.formula(paste(Outcome[i]," ~ ", xm));c = lm(fmla2, Data_v);
      if (sum(is.na(summary(c)[[4]][,2]))>0) {
        
        if (Group=='All') {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq,q),collapse= "+")),collapse= "+")} else 
        {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq),collapse= "+")),collapse= "+")} 
        fmla2 <- as.formula(paste(Outcome[i]," ~ ", xm));c = lm(fmla2, Data_v);}
      
      if (sum(is.na(summary(c)[[4]][,2]))>0) {
        sela=c$coefficients[!is.na(c$coefficients)]
        sela2=sela[names(sela) %in% colnames(M)]
        sv1=names(sela2[order(sela2)][1]);sv2=names(sela2[rev(order(sela2))][1])
        jap=c(sv1)
        sela3=sela[names(sela) %in% colnames(X)]
        sv1=names(sela3[order(sela3)][1]);sv2=names(sela3[rev(order(sela3))][1])
        jap1=c(sv1,sv2); XM=c(jap,jap1);XM=XM[!is.na(XM)]
        if (sum(XM %in% colnames(X))==0) {sv2=names(xt[rev(order(xt))[1]]);jap1=c(sv1,sv2); XM=c(jap,jap1);XM=XM[!is.na(XM)]}
        
        if (Group=='All') {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq,q),collapse= "+")),collapse= "+")} else 
        {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq),collapse= "+")),collapse= "+")} 
        fmla2 <- as.formula(paste(Outcome[i]," ~ ", xm));c = lm(fmla2, Data_v);}
    }
    if (sum(is.na(c$coefficients))>0) {sela=c$coefficients[!is.na(c$coefficients)]
    sela2=sela[names(sela) %in% colnames(M)]
    sv1=names(sela2[order(sela2)][1:2]);sv2=names(sela2[rev(order(sela2))][1:2])
    jap=c(sv1,sv1)
    sela3=sela[names(sela) %in% colnames(X)]
    sv1=names(sela3[order(sela3)][1]);sv2=names(sela3[rev(order(sela3))][1])
    jap1=c(sv1,sv2); XM=c(jap,jap1);XM=XM[!is.na(XM)]
    if (sum(XM %in% colnames(X))==0) {sv2=names(xt[rev(order(xt))[1]]);jap1=c(sv1,sv2); XM=c(jap,jap1);XM=XM[!is.na(XM)]}
    if (Group=='All') {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq,qq,q),collapse= "+")),collapse= "+")} else 
    {xm=paste(c(paste(c(XM), collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")} 
    fmla2 <- as.formula(paste(Outcome[i]," ~ ", xm));c = lm(fmla2, Data_v);}
    
    if (sum(is.na(c$coefficients))>0) {next}
    
    #https://bookdown.org/egarpor/PM-UC3M/lm-ii-modsel.html# rm(Mn,Xn)  
    sela=c$coefficients[!is.na(c$coefficients)];sela2=sela[names(sela) %in% colnames(M)];sela3=sela[names(sela) %in% colnames(X)]
    Mn=names(sela2); Xn=names(sela3)
    if (Group=='All') { fmla1 <- as.formula(paste(paste(Mn, collapse= "+")," ~ ", paste(c(paste(Xn, collapse= "+"),
                                                                                          paste(c(qqq,qq,q),collapse= "+")),collapse= "+")))} else 
                                                                                          {fmla1 <- as.formula(paste(paste(Mn, collapse= "+")," ~ ", paste(c(paste(Xn, collapse= "+"),paste(c(qqq,qq),collapse= "+")),collapse= "+")))};
    
    b = lm(fmla1, Data_v);
    # print(summary(b))
    if (sum(is.na(summary(b)[[4]][,2]))>0) {next}
    
    rn=c();med_out=c();res=c(); tmp=c();rn=c();med_oute=c();ji=1 
    # try({
    for (ji in 1:length(Mn)) { #control.value=mina[i]
      for (z in 1:length(Xn)) {
        med_oute=mediate(b, c, treat =  colnames(X)[colnames(X)==Xn[z]], mediator = colnames(M)[colnames(M)==Mn[ji]],
                         sims = simss,control.value=control.value[colnames(X)==Xn[z]],treat.value=treat.value[colnames(X)==Xn[z]]) #} else {next} control.value=control.value[be[z]],treat.value=treat.value[be[z]]
        med_out = summary(med_oute) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
        tmp=c(med_out$d0, med_out$d0.p, med_out$d0.ci[1],med_out$d0.ci[2],
              med_out$z0, med_out$z0.p, med_out$z0.ci[1],med_out$z0.ci[2],med_out$n1, med_out$n1.p,med_out$n1.ci[1],
              med_out$n1.ci[2],med_out$tau.coef,med_out$tau.p,med_out$tau.ci[1],med_out$tau.ci[2]) 
        res <- rbind(res,tmp);
        # print(Outcome[i])
        rn=append(rn,paste(colnames(X)[colnames(X)==Xn[z]],colnames(M)[colnames(M)==Mn[ji]],Outcome[i], sep=" ")) #so this seems to be for all lipids...
        remove(tmp) }} #,  {print('medc');next})  #https://intro2r.com/loops.html https://www.benjaminbell.co.uk/2022/12/loops-in-r-nested-loops.html
    # try({res},{res=t(data.frame(rep(NA,16)))})
    # if (isTRUE(res))  {res=t(data.frame(rep(0,16)))}
    
    
    try({rownames(res)=rn}, {print('rn')}) #write.csv(rn,'iii.csv')
    try({rownames(res) <- gsub("X11", "11", rownames(res))}, {print('xx')})
    try({rownames(res) <- gsub("X17", "17", rownames(res)) }, {print('xxx')}) 
    koko=rbind(koko,res)
  } #:) # }
  
  try({colnames(koko)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') }, {print('col')})
  try({koko=koko[order(koko[,2]),]},{print('order')}) #res=res[rev(order(res[,1])),]
  tryCatch({write.xlsx(koko, file = paste(fn,Group,'sick_',sick,name,date,'teste_4.xlsx'), append = FALSE, row.names = TRUE) }, error = function(msg){return(write.csv(koko,paste(name,Group,date,'test_lpcaa.csv')))}) # :)
  # }
  return(koko) }

#This function needs some changes...:
plottings=function(uh7ma,uh7m,uh7f,test,take,date) {
  rtot=uh7ma #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot) # name=paste(simss,'basic hypothesis',take)# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model# https://github.com/MarioniLab/miloR# https://www.nature.com/articles/s41467-023-40458-9/figures/4
  name=paste('Contaminants_Steroids_BAs_or_Lipids_100sims_basic_all',test,take,date) # rtot=rtot_2000_mrct # rtot=uh5
  pdf(paste("ACME_p_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'d0.p']),breaks=100, xlab="P value of ACME", ylab="Frequency",main="");dev.off()
  pdf(paste("Mediated Proportion_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'Proportion Mediated']),breaks=300,xlim=c(-1.2,1.2),xlab="Proportion of Mediated", ylab="Frequency",main="");dev.off()
  med.freq.cut=0.4; med.min=0; med.prop=0.1; Group='All'; med='Mediated' #or not mediated e.g. no
  if (med.min=='All') {rt2=rtot} else if (med.min!='All')  {if (med=='Mediated') {rt2=rtot[rtot[,'ACME']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']>med.prop,]} else {rt2=rtot[rtot[,'ADE']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']<med.prop,]}}
  if (med=='Mediated') {rt2=rt2[rt2[,'d0.p']<med.freq.cut,]; rt2=rt2[rt2[,'z0.p']>med.freq.cut,];rt2=rt2[rt2[,'ACME']>med.min,]} else {rt2=rt2[rt2[,'d0.p']>med.freq.cut,]; rt2=rt2[rt2[,'z0.p']<med.freq.cut,];rt2=rt2[rt2[,'ADE']>med.min,]}
  # rt2=rtot
  hoi=c(); hoi=scan(text=rownames(rt2), what="")
  hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
  colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender') ## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
  df2 <- hoi %>%make_long('Contaminants','Steroids','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
  meda='Sankey plot of'
  pdf(paste(meda,name,med.freq.cut,med.prop ,".pdf"), width = 20, height = 20,  pointsize = 10);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  jpeg(paste(meda,name,med.freq.cut,med.prop ,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 12, res=500);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
          geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
          scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  
  
  rtot=uh7m #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot) # name=paste(simss,'basic hypothesis',take)# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model# https://github.com/MarioniLab/miloR# https://www.nature.com/articles/s41467-023-40458-9/figures/4
  name=paste('Contaminants_Steroids_BAs_or_Lipids_100sims_basic_malea',test,take,date) # rtot=rtot_2000_mrct # rtot=uh5
  pdf(paste("ACME_p_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'d0.p']),breaks=100, xlab="P value of ACME", ylab="Frequency",main="");dev.off()
  pdf(paste("Mediated Proportion_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'Proportion Mediated']),breaks=300,xlim=c(-1.2,1.2),xlab="Proportion of Mediated", ylab="Frequency",main="");dev.off()
  med.freq.cut=0.3; med.min=0; med.prop=0.1; Group='All'; med='Mediated' #or not mediated e.g. no
  if (med.min=='All') {rt2=rtot} else if (med.min!='All')  {if (med=='Mediated') {rt2=rtot[rtot[,'ACME']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']>med.prop,]} else {rt2=rtot[rtot[,'ADE']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']<med.prop,]}}
  if (med=='Mediated') {rt2=rt2[rt2[,'d0.p']<med.freq.cut,]; rt2=rt2[rt2[,'z0.p']>med.freq.cut,];rt2=rt2[rt2[,'ACME']>med.min,]} else {rt2=rt2[rt2[,'d0.p']>med.freq.cut,]; rt2=rt2[rt2[,'z0.p']<med.freq.cut,];rt2=rt2[rt2[,'ADE']>med.min,]}
  hoi=c(); hoi=scan(text=rownames(rt2), what="")
  hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
  colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender') ## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
  df2 <- hoi %>%make_long('Contaminants','Steroids','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
  meda='Sankey plot of'
  pdf(paste(meda,name,med.freq.cut,med.prop ,"male.pdf"), width = 20, height = 20,  pointsize = 10);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  jpeg(paste(meda,name,med.freq.cut,med.prop ,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 12, res=500);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
          geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
          scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  
  rtot=uh7f #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot) # name=paste(simss,'basic hypothesis',take)# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model# https://github.com/MarioniLab/miloR# https://www.nature.com/articles/s41467-023-40458-9/figures/4
  name=paste('Contaminants_Steroids_BAs_or_Lipids_100sims_basic_female',test,take,date) # rtot=rtot_2000_mrct # rtot=uh5
  pdf(paste("ACME_p_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'d0.p']),breaks=100, xlab="P value of ACME", ylab="Frequency",main="");dev.off()
  pdf(paste("Mediated Proportion_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'Proportion Mediated']),breaks=300,xlim=c(-1.2,1.2),xlab="Proportion of Mediated", ylab="Frequency",main="");dev.off()
  med.freq.cut=0.3; med.min=0; med.prop=0.1; Group='All'; med='Mediated' #or not mediated e.g. no
  if (med.min=='All') {rt2=rtot} else if (med.min!='All')  {if (med=='Mediated') {rt2=rtot[rtot[,'ACME']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']>med.prop,]} else {rt2=rtot[rtot[,'ADE']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']<med.prop,]}}
  # if (med=='Mediated') {rt2=rt2[rt2[,'d0.p']<med.freq.cut,]; rt2=rt2[rt2[,'z0.p']>med.freq.cut,];rt2=rt2[rt2[,'ACME']>med.min,]} else {rt2=rt2[rt2[,'d0.p']>med.freq.cut,]; rt2=rt2[rt2[,'z0.p']<med.freq.cut,];rt2=rt2[rt2[,'ADE']>med.min,]}
  hoi=c(); hoi=scan(text=rownames(rt2), what="")
  hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
  colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender') ## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
  df2 <- hoi %>%make_long('Contaminants','Steroids','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
  meda='Sankey plot of'
  pdf(paste(meda,name,med.freq.cut,med.prop ,"female.pdf"), width = 20, height = 20,  pointsize = 10);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  jpeg(paste(meda,name,med.freq.cut,med.prop ,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 12, res=500);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
          geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
          scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  
  rtot=rbind(uh7m,uh7f)
  name=paste('Contaminants_Steroids_BAs_or_Lipids_100sims_basic_males and females',test,date) # rtot=rtot_2000_mrct # rtot=uh5
  pdf(paste("ACME_p_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'d0.p']),breaks=100, xlab="P value of ACME", ylab="Frequency",main="");dev.off()
  pdf(paste("Mediated Proportion_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'Proportion Mediated']),breaks=300,xlim=c(-1.2,1.2), xlab="Proportion of Mediated", ylab="Frequency",main="");dev.off()
  med.freq.cut=0.3; med.min=0; med.prop=0.1; Group='All'; med='Mediated' #or not mediated e.g. no
  if (med.min=='All') {rt2=rtot} else if (med.min!='All')  {if (med=='Mediated') {rt2=rtot[rtot[,'ACME']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']>med.prop,]} else {rt2=rtot[rtot[,'ADE']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']<med.prop,]}}
  # if (med=='Mediated') {rt2=rt2[rt2[,'d0.p']<med.freq.cut,]; rt2=rt2[rt2[,'z0.p']>med.freq.cut,];rt2=rt2[rt2[,'ACME']>med.min,]} else {rt2=rt2[rt2[,'d0.p']>med.freq.cut,]; rt2=rt2[rt2[,'z0.p']<med.freq.cut,];rt2=rt2[rt2[,'ADE']>med.min,]}
  hoi=c(); hoi=scan(text=rownames(rt2), what="")
  hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
  hoi[,4]=rt2[,'Gender'];hoi=hoi[,c(4,1,2,3)]
  colnames(hoi)=c('Gender','Contaminants','Steroids','Bile Acids or Lipids')#,'Gender') #
  df2 <- hoi %>%make_long('Gender','Contaminants','Steroids','Bile Acids or Lipids')
  meda='Sankey plot of'
  pdf(paste(meda,name,med.freq.cut,med.prop ,".pdf"), width = 20, height = 20,  pointsize = 10);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  jpeg(paste(meda,name,med.freq.cut,med.prop ,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 12, res=500);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
          geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
          scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  
  rtot=uha #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot)
  name=paste('Contaminants_Steroids_BAs_or_Lipids_sims_basic direct_allaa',test,take,date) # rtot=rtot_2000_mrct # rtot=uh5
  pdf(paste("ACME_p_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'d0.p']),breaks=100, xlab="P value of ACME", ylab="Frequency",main="");dev.off()
  pdf(paste("Mediated Proportion_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'Proportion Mediated']),breaks=300,
                                                                                                         xlim=c(-1.2,1.2),xlab="Proportion of Mediated", ylab="Frequency",main="");dev.off()
  # med.freq.cut=0.4; #rt2=rt2[rt2[,'d0.p']>med.freq.cut,];
  med.freq.cut=0.2;  med.min=1; med.prop=0.3;  med='no' #or not mediated e.g. no
  rt2=rtot
  if (med.min=='All') {rt2=rtot} else if (med.min!='All')  {if (med=='Mediated') {rt2=rtot[rtot[,'ACME']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']>med.prop,]} else {rt2=rtot[rtot[,'ADE']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']<med.prop,]}}
  if (med=='Mediated') {rt2=rt2[rt2[,'d0.p']<med.freq.cut,]; rt2=rt2[rt2[,'z0.p']>med.freq.cut,];rt2=rt2[rt2[,'ACME']>med.min,]} else { rt2=rt2[rt2[,'z0.p']<med.freq.cut,];rt2=rt2[rt2[,'ADE']>med.min,]}
  hoi=c(); hoi=scan(text=rownames(rt2), what="")
  hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
  colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender')
  hoi=hoi[,c('Contaminants','Bile Acids or Lipids')]## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
  df2 <- hoi %>%make_long('Contaminants','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
  meda='Sankey plot of'
  pdf(paste(meda,name,med.freq.cut,med.prop ,"all.pdf"), width = 20, height = 20,  pointsize = 10);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  jpeg(paste(meda,name,med.freq.cut,med.prop ,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 12, res=500);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
          geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
          scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  
  rtot=rbind(uh7m,uh7f)
  name=paste('Contaminants_Steroids_BAs_or_Lipids_sims_basic_direct_males and females',test,date) # rtot=rtot_2000_mrct # rtot=uh5
  pdf(paste("ACME_p_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'d0.p']),breaks=100, xlab="P value of ACME", ylab="Frequency",main="");dev.off()
  pdf(paste("Mediated Proportion_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'Proportion Mediated']),breaks=300,xlim=c(-1.2,1.2), xlab="Proportion of Mediated", ylab="Frequency",main="");dev.off()
  med.freq.cut=0.3; med.min=0; med.prop=0.1;  med='no' #or not mediated e.g. no
  # rt2=rtot
  if (med.min=='All') {rt2=rtot} else if (med.min!='All')  
  {if (med=='Mediated') {rt2=rtot[rtot[,'ACME']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']>med.prop,]} else {rt2=rtot[rtot[,'ADE']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']<med.prop,]}}
  if (med=='Mediated') {rt2=rt2[rt2[,'d0.p']<med.freq.cut,]; rt2=rt2[rt2[,'z0.p']>med.freq.cut,];rt2=rt2[rt2[,'ACME']>med.min,]} else 
  {rt2=rt2[rt2[,'d0.p']>med.freq.cut,]; rt2=rt2[rt2[,'z0.p']<med.freq.cut,];rt2=rt2[rt2[,'ADE']>med.min,]}
  hoi=c(); hoi=scan(text=rownames(rt2), what="")
  hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
  colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender')
  hoi=hoi[,c('Contaminants','Bile Acids or Lipids')]## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
  df2 <- hoi %>%make_long('Contaminants','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
  meda='Sankey plot of'
  pdf(paste(meda,name,med.freq.cut,med.prop ,".pdf"), width = 20, height = 20,  pointsize = 10);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  jpeg(paste(meda,name,med.freq.cut,med.prop ,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 12, res=500);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
          geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") + 
          scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+
          theme(axis.title.x = element_blank()));dev.off()}

#Mediaatio tulosten ilmoitus ja oleellisten 'poiminta', reduced funktiolla (3:n yksiselitteisempi..)...
reduced2=function(u3,Group,name,lkm) {
  u3=uh7mae
  ACMEMedian=c();ACMEpval=c();ACMEVar=c()
  ADEMedian=c();ADEpval=c();ADEVar=c()
  PMEDMedian=c();PMEDpval=c();PMEDVar=c(); c1=c()
  ADEMedian=median(u3[,'ADE'])+1*sd(u3[,'ADE']) 
  ADEpval=quantile(u3[,'z0.p'],0.1,na.rm=TRUE) #just want to have the negative ADEs... so this is not needed
  # DV=scale(abs(u3[,'z0.ci_l']-u3[,'z0.ci_u'])) #not sure if this is needed 12224
  DV=u3[,'d0.ci_l']/abs(u3[,'z0.ci_u'])
  ok=median(DV[DV>-4])
  # which(u3[,'z0.p']<0.05  & DV>-1.5)
  ce=u3[,'z0.p']<ADEpval  & DV>ok
  c1a= u3[ce,]
  # ADEVar=round(quantile(DV,0.975,na.rm=TRUE),1)
  c1b= u3[u3[,'ADE']<ADEMedian,]
  c1=rbind(c1b)
  # print(ADEMedian)
  # c1 = u3[u3[,'ADE']<ADEMedian,]
  
  PMEDMedian=median(c1[,'Proportion Mediated'][c1[,'Proportion Mediated']>0]) #ADE and ACME values twice in p.med?
  # PMEDpval=quantile(c1[,'n1.p'],0.95,na.rm=TRUE)
  # PV=scale(abs(c1[,'n.ci_l']-c1[,'n1.ci_u']))
  # PMEDVar=round(quantile(PV,0.975,na.rm=TRUE),1) #not sure if this is needed
  # c1=c1[c1[,'Proportion Mediated']>PMEDMedian & c1[,'n1.p']<PMEDpval & PV<PMEDVar,];
  c1=c1[c1[,'Proportion Mediated']>PMEDMedian,];
  
  ACMEMedian=median(c1[,'ACME'][c1[,'ACME']>0])
  ACMEpval=quantile(c1[,'d0.p'],0.5,na.rm=TRUE) #should be less than median
  # AV=scale(abs(c1[,'d0.ci_l']-c1[,'d0.ci_u']))
  # ACMEVar=round(quantile(AV,0.975,na.rm=TRUE),1)
  # c1=c1[c1[,'ACME']>ACMEMedian & c1[,'d0.p']<ACMEpval & AV<ACMEVar,] 
  c1=c1[c1[,'ACME']>ACMEMedian & c1[,'d0.p']<ACMEpval,] 
  io=c(); for (i in 1:length(rownames(c1))) {io=append(io,grepl( '_L', rownames(c1)[i], fixed = TRUE))}
  try({c1=c1[!io,]}, {c1=c1})
  c1=c1[rev(order(c1[,1])),];
  c1=tryCatch({c1[1:lkm,]}, error = function(msg){return(c1)}) #https://cnuge.github.io/post/trycatch/
  try({c1},{c1=t(data.frame(rep(0,16)))})
  # write.xlsx(c1, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  return(c1)}

reduced3=function(u3,Group,name,lkm) {
  ACMEMedian=c();ACMEpval=c();ACMEVar=c()
  ADEMedian=c();ADEpval=c();ADEVar=c()
  PMEDMedian=c();PMEDpval=c();PMEDVar=c()
  c1=c(); c1=u3
  ACMEMedian=median(c1[,'ACME'][c1[,'ACME']>0])
  ACMEpval=quantile(c1[,'d0.p'],0.5,na.rm=TRUE) #should be less than median
  # AV=scale(abs(c1[,'d0.ci_l']-c1[,'d0.ci_u'])) #I would not take this... 12224
  # ACMEVar=round(quantile(AV,0.95,na.rm=TRUE),1) #I would not take this...
  c1=c1[c1[,'ACME']>ACMEMedian & c1[,'d0.p']<ACMEpval,] 
  # c1=c1[c1[,'ACME']>ACMEMedian & c1[,'d0.p']<ACMEpval & AV<ACMEVar,] 
  io=c()
  for (i in 1:length(rownames(c1))) {io=append(io,grepl( '_L', rownames(c1)[i], fixed = TRUE))}
  try({c1=c1[!io,]}, {c1=c1})
  c1=c1[rev(order(c1[,1])),]
  c1=tryCatch({c1[1:lkm,]}, error = function(msg){return(c1)}) #try({res},{res=t(data.frame(rep(NA,16)))}) return(c1[1:dim(c1)[1],])}
  try({c1},{c1=t(data.frame(rep(0,16)))})
  write.xlsx(c1, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  return(c1)}

direct_sankey=function(rtot, med.freq.cut=0.3, med='no',med.min=0, med.prop=0.4,name,lkm) {
  rt2=rtot
  if (med.min=='All') {rt2=rtot} else if (med.min!='All')  {if (med=='Mediated') 
  {rt2=rtot[rtot[,'ACME']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']>med.prop,];rt2=rt2[rt2[,'d0.p']<med.freq.cut,];rt2=rt2[rt2[,'z0.p']>med.freq.cut,]} else 
  {rt2=rtot[rtot[,'ADE']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']<med.prop,];rt2=rt2[rt2[,'d0.p']>med.freq.cut,];rt2=rt2[rt2[,'z0.p']<med.freq.cut,]}}
  rt2=rt2[rev(order(rt2[,'ADE'])),]; 
  rt2=rt2[1:lkm,]
  hoi=c(); hoi=scan(text=rownames(rt2), what="")
  hoi=as.data.frame(matrix(hoi, ncol = 4,  byrow = TRUE), stringsAsFactors = FALSE)
  colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender')
  hoi=hoi[,c('Contaminants','Bile Acids or Lipids')]
  ## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct
  # https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation
  # https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
  df2 <- hoi %>%make_long('Contaminants','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
  meda='Sankey plot of'
  pdf(paste(meda,name,med.freq.cut,med.prop ,"all.pdf"), width = 20, height = 20,  pointsize = 10);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  jpeg(paste(meda,name,med.freq.cut,med.prop ,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 12, res=500);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
          geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
          scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  return(rt2)}

# testing testing...
# rownames(rt2[grep( 'putative',rownames(rt2)),])
# dy <- replace(rownames(rt2[grep( 'putative',rownames(rt2)),]), 2, 'blueberry')
# df$Names <- replace(df$Names, df$Names %in% conditions, replacement_values)
# hop = substr(rownames(rt2),grep( 'putative',rownames(rt2)),-nchar('putative'))
# nchar('putative')
# library('stringi')
# rownames(rt2[grep( 'putative',rownames(rt2)),]) = stri_sub(rownames(rt2[grep( 'putative',rownames(rt2)),]), 1, -(nchar('putative')+1)) 
# jai=c()
# for (i in 1:length(rownames(rt2))) {jai=append(jai,strsplit(rownames(rt2)[i], " +"))}
# for (i in 1:length(rownames(rt2))) if(length(jai[i][[1]])==3) jai[i][[1]][4]='ok'
# jeo=t(data.frame(jai))
# rownames(jeo)=1:length(rownames(jeo))
# jeo=jeo[,1:3]

plottings_sf=function(uh7ma,date,sick,Group) {
  rt2=alma #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot) # name=paste(simss,'basic hypothesis',take)# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model# https://github.com/MarioniLab/miloR# https://www.nature.com/articles/s41467-023-40458-9/figures/4
  name=paste('Contaminants_Steroids_BAs_or_Lipids_sims',date) # rtot=rtot_2000_mrct # rtot=uh5
  hoi=c(); 
  hoi=scan(text=rt2[,1], what=" ")#rownames(rt2)# names(
  hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
  
  # jai=c()
  # for (i in 1:length(rownames(rt2))) {jai=append(jai,strsplit(rownames(rt2)[i], " +"))}
  # for (i in 1:length(rownames(rt2))) if(length(jai[i][[1]])==3) jai[i][[1]][4]='ok'
  # jeo=t(data.frame(jai))
  # rownames(jeo)=1:length(rownames(jeo))
  # jeo=jeo[,1:3]
  # hoi=as.data.frame(jeo)
  colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender') ## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
  
  hoi[,'Steroids' ][hoi[,'Steroids' ]=='17aOH.P4']='17a-OHP4'
  # hoi[,'Steroids' ][hoi[,'Steroids' ]=='17a-OHP5']='17-aOHP5'
  hoi[,'Steroids' ]  <- gsub("\\.", "-",  hoi[,'Steroids' ] ) #:)
  hoi[,'Steroids' ][ hoi[,'Steroids' ]=='T-Epi-T']='T/Epi-T'
  
  df2 <- hoi %>%make_long('Contaminants','Steroids','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
  # meda='Sankey plot of'
  # pdf(paste(meda,name,sick,Group,".pdf"), width = 20, height = 20,  pointsize = 18);
  # print(
  ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() +
          theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank())#; dev.off()
  # jpeg(paste(meda,name ,sick,Group,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 18, res=500);
  # print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  #         scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));
  # dev.off()
  }
# https://stackoverflow.com/questions/8370548/how-can-i-interrupt-a-running-code-in-r-with-a-keyboard-command

reduced2=function(u3,Group,name,lkm) {
  c1=c()
  ACMEMedian=c();ACMEpval=c();ACMEVar=c()
  ADEMedian=c();ADEpval=c();ADEVar=c()
  # PMEDMedian=c();PMEDpval=c();PMEDVar=c(); 
  
  # ADEMedian=0#median(u3[,'ADE'])+1.5*sd(u3[,'ADE'])
  
  # ADEpval=quantile(u3[,'z0.p'],0.975) #just want to have the negative ADEs
  # DV=scale(abs(u3[,'z0.ci_l']-u3[,'z0.ci_u']))
  # ADEVar=round(quantile(DV,0.975),1)
  c1= u3#[u3[,'ADE'] < ADEMedian  & DV<ADEVar,] #& u3[,'z0.p']<ADEpval
  
  # PMEDMedian=median(c1[,'Proportion Mediated'][c1[,'Proportion Mediated']>0]) #ADE and ACME values twice in p.med?
  # PMEDpval=quantile(c1[,'n1.p'],0.95)
  # PV=scale(abs(c1[,'n.ci_l']-c1[,'n1.ci_u']))
  # PMEDVar=round(quantile(PV,0.975),1)
  # c1=c1[c1[,'Proportion Mediated']>PMEDMedian & c1[,'n1.p']<PMEDpval & PV<PMEDVar,];
  ACMEMedian=0#median(c1[,'ACME'][c1[,'ACME']>0])
  # ACMEpval=quantile(c1[,'d0.p'],0.5) #should be less than median
  # AV=scale(abs(c1[,'d0.ci_l']-c1[,'d0.ci_u']))
  # ACMEVar=round(quantile(AV,0.975),1)
  # asdf=c1
  # c1=asdf
  # c1[,'ACME']>ACMEMedian & 
  c1=c1[c1[,'ACME']>ACMEMedian & ((c1[,'ACME']-c1[,'ADE']) > 0), ] #&(-median(c1[,'ACME']))) 
          # (c1[,'z0.p']-c1[,'d0.p'])>-0.1,]# &
          # (c1[,'z0.ci_u']-c1[,'d0.ci_l'])<2.5,] #&
          # (c1[,'z0.ci_l']/abs(c1[,'z0.ci_u']-c1[,'z0.ci_l'])-
             # c1[,'d0.ci_l']/abs(c1[,'d0.ci_u']-c1[,'d0.ci_l']))>-0.05,]#& c1[,'d0.p']<ACMEpval & AV<ACMEVar
  c1=tryCatch({c1[1:lkm,]}, error = function(msg){return(c1)})
  c1=c1[rev(order(c1[,'ACME'])),];
  write.xlsx(c1, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE)
  
  return(c1)}

# reduced2=function(u3,Group,name,lkm) {
#   # u3=all_all1
#   c1=c()
#   
#   DV=-u3[,'z0.ci_l']/(abs(u3[,'z0.ci_u']-u3[,'z0.ci_l']))
#   AV=-u3[,'d0.ci_l']/(abs(u3[,'d0.ci_u']-u3[,'d0.ci_l']))
#   OK=DV-AV
#   
#   pmeda=(u3[,'ACME']+abs(min(u3[,'ACME'])))/(u3[,'ACME']+abs(min(u3[,'ACME']))+u3[,'ADE']+abs(min(u3[,'ADE'])))
#   c1= u3 
#   c1=c1[c1[,'ACME']>quantile(c1[,'ACME'],0.75) & c1[,'d0.p']<quantile(c1[,'d0.p'],0.25) & 
#           (c1[,'ACME']-c1[,'ADE']) > quantile(c1[,'ACME']-c1[,'ADE'],0.2) & 
#           OK>quantile(OK,0.20) & pmeda > quantile(pmeda,0.20), ] #& 
#   c1=c1[rev(order(c1[,'ACME'])),];
#   
#   #putatatiivien (tai melkein niin) poisto
#   # rownames(c1)
#   # the_ones=grepl('ok', rownames(c1))
#   # c1=c1[the_ones,]
#   # rownames(c1)  <- gsub(" ok", "",  rownames(c1) ) #:)
#   # the_ones=grepl('Y_model_with_forced_int_m_and_x_ok', rownames(c1)); c1=c1[!the_ones,];
#   # rownames(c1)  <- gsub(" Y_model_with_forced_int_x_ok", "",  rownames(c1) );rownames(c1)  <- gsub(" Y_model_with_forced_int_m_ok", "",  rownames(c1) ) #:)
#   
#   #ei poistoa, mutta piilotus:
#   # rownames(c1)  <- gsub(" ok", "",  rownames(c1) ); 
#   # rownames(c1)  <- gsub(" Y_model_with_forced_int_x_ok", "",  rownames(c1) );
#   # rownames(c1)  <- gsub(" Y_model_with_forced_int_m_ok", "",  rownames(c1) ) 
#   # rownames(c1)  <- gsub(" putative_m", "",  rownames(c1) ); 
#   rt2=c1 #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot) # name=paste(simss,'basic hypothesis',take)# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model# https://github.com/MarioniLab/miloR# https://www.nature.com/articles/s41467-023-40458-9/figures/4
#   # name=paste('Contaminants_Steroids_BAs_or_Lipids_sims',date) # rtot=rtot_2000_mrct # rtot=uh5
#   
#   # hoi=c(); 
#   # hoi=rownames(rt2)#scan(text=names(rt2[,1]), what=" ")#rownames(rt2)#
#   # hoi=as.data.frame(matrix(hoi, ncol = 1,  byrow = TRUE), stringsAsFactors = FALSE)
#   # hoi=as.data.frame(matrix(hoi, ncol = 4,  byrow = TRUE), stringsAsFactors = FALSE)
#   # ug=names(table(hoi$V4))
#   
#   # for (i in 1:length(ug)) {rownames(c1)  <- gsub(paste('',ug[i]), '' ,rownames(c1) )}
#   
#   
#   c1=tryCatch({c1[1:lkm,]}, error = function(msg){return(c1)}) 
#   
#   write.xlsx(c1, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE)
#   
#   return(c1)}


plottings_sf=function(uh7ma,date,sick,Group) {
  rt2=alma #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot) # name=paste(simss,'basic hypothesis',take)# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model# https://github.com/MarioniLab/miloR# https://www.nature.com/articles/s41467-023-40458-9/figures/4
  name=paste('Contaminants_Steroids_BAs_or_Lipids_sims',date) # rtot=rtot_2000_mrct # rtot=uh5
  
  hoi=c(); 
  # hoi=scan(text=rownames(rt2), what=" ")#rownames(rt2)# names(rt2[,1])
  # hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
  
  hoi=scan(text=rt2[,1] , what=" ")#rownames(rt2)# names(rt2[,1]) rownames(rt2)
  hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
  
  # rownames(c1)
  # the_ones=grepl('ok', rownames(c1))
  # c1=c1[the_ones,]
  # rownames(c1)  <- gsub(" ok", "",  rownames(c1) ) #:)
  # the_ones=grepl('Y_model_with_forced_int_m_and_x_ok', rownames(c1)); c1=c1[!the_ones,];
  # rownames(c1)  <- gsub(" Y_model_with_forced_int_x_ok", "",  rownames(c1) );rownames(c1)  <- gsub(" Y_model_with_forced_int_m_ok", "",  rownames(c1) ) #:)  
  
  colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender') ## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
  
  hoi[,'Steroids' ][hoi[,'Steroids' ]=='17aOH.P4']='17a-OHP4'
  hoi[,'Steroids' ][hoi[,'Steroids' ]=='17aOH-P4']='17a-OHP4'
  hoi[,'Steroids' ]  <- gsub("\\.", "-",  hoi[,'Steroids' ] ) #:)
  hoi[,'Steroids' ][ hoi[,'Steroids' ]=='T-Epi-T']='T/Epi-T'
  
  df2 <- hoi %>%make_long('Contaminants','Steroids','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
  meda='Sankey plot of_KAIKKI'
  pdf(paste(meda,name,sick,Group,".pdf"), width = 20, height = 20,  pointsize = 18);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + 
          theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  
  windowsFonts(A = windowsFont("Calibri (Body)")) 
  # text(font = 2)
  jpeg(paste(meda,name ,sick,Group,".jpg"), width = 6000, height = 7300, quality = 100,pointsize = 20, res=500);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,
                        fill = factor(node),label = node)) +
          
          geom_sankey(flow.alpha = 1, node.color = 1) + 
          geom_sankey_label(size = 8.5, color = 1, fill = "white") +
          # scale_fill_viridis_d(option = "D", alpha = 0.95) + 
          theme_sankey(base_size = 30) + 
          scale_fill_grey(start = 0.5, end = 0.5)+
          theme(axis.text.x = element_text(hjust = 0.5, vjust=7,colour = 'black'))+ #https://stackoverflow.com/questions/38862303/customize-ggplot2-axis-labels-with-different-colors
          theme(legend.position = "none") +
          theme(axis.title.x = element_blank()));
  dev.off()
}

ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
  geom_sankey(flow.alpha = 1, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +
  # scale_fill_viridis_d(option = "J", alpha = 1) + 
  scale_fill_grey(start = 0.5, end = 0.5)+
  theme_sankey(base_size = 18) + 
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())

#Testing hyp 3 and 4 with a combined function:
the_essentials=function(Treatment, Mediator, Outcome,tv, Group,name,simss,t.val,test,sick,sick_group,fn,lkm,date,joo,ip) {
  hoi1=paste("C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/",fn,sep='')
  setwd(hoi1) #https://topepo.github.io/caret/parallel-processing.html
  # https://ds-pl-r-book.netlify.app/optimization-in-r.html
  # t.val='no';simss=100; test=''; take=''; name='100hypo4_yes'; sick='yes'
  # Group='All'
  # t.val='no'; simss=10;
  # test=''; take=''; 
  name=paste(simss,'hypo4_yes_sick'); sick='yes'
  # loop_med_simplified_hyp4a(Group,tv,name,date,simss,sick,sick_group,joo,ip)
  # loop_med_simplified_hyp4(Treatment,Mediator,Outcome,tv,Group,name,date,simss,t.val,test,sick,sick_group,joo)
  # loop_med_simplified2b=function(Group,name,date,simss,sick,sick_group,joo,ip)
  Group='All'; uh7ma=loop_med_simplified1a(Group,tv,name,date,simss,sick,sick_group,joo,ip);try({uh7ma},{uh7ma=data.frame(0)});  save(uh7ma,file=paste(fn,'sick all.RData'));#try({uh7ma},{uh7ma=data.frame(0)})
  # # sum(tv[,2]==1 & sick_group)# Treatment, Mediator, Outcome,tv,Group,name,date,simss,t.val,test,sick,sick_group
  Group='female'; uh7f=loop_med_simplified1a(Group,tv,name,date,simss,sick,sick_group,joo,ip);  try({uh7f},{uh7f=data.frame(0)});save(uh7f, file=paste(fn,'sick female.RData'));
  # # uh7f=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv,Group,name,date,simss,t.val,test,sick,sick_group);try({uh7f},{uh7f=data.frame(0)})
  Group='male';   uh7m=loop_med_simplified1a(Group,tv,name,date,simss,sick,sick_group,joo,ip);try({uh7m},{uh7m=data.frame(0)});save(uh7m, file=paste(fn,'sick male.RData'));
  # # > Group='Female'; uh7f=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv,Group,name,date,simss,t.val,test,sick,sick_group);
  # Error in eigen(sigma, symmetric = TRUE) :  infinite or missing values in 'x'# In addition: Warning messages: 1: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
  # 2: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
  # > Group='Male';   uh7m=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv,Group,name,date,simss,t.val,test,sick,sick_group);# Error in mediate(b, c, treat = x[be[z]], mediator = m[ce[j]], sims = simss,  : NA in model coefficients; rerun models with nonsingular design matrix;n addition: Warning messages: 1: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold

  t.val='no'; #simss=10;
  test=''; take=''; name=paste(simss,'hypo4_no_not sick'); sick='no'
  Group='All';    uh7mae=loop_med_simplified1a(Group,tv,name,date,simss,sick,sick_group,joo,ip);save(uh7mae, file=paste(fn,'healthy alla.RData'));try({uh7mae},{uh7mae=data.frame(0)})
  Group='female'; uh7fe=loop_med_simplified1a(Group,tv,name,date,simss,sick,sick_group,joo,ip);save(uh7fe, file=paste(fn,'healthy female.RData'));try({uh7fe},{uh7fe=data.frame(0)})
  Group='male';   uh7me=loop_med_simplified1a(Group,tv,name,date,simss,sick,sick_group,joo,ip);save(uh7me, file=paste(fn,'healthy male.RData'));try({uh7me},{uh7me=data.frame(0)})
  # # setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/")
  # load("100hypo4_steatosis_all_22124.RData")
  # uh7mae[1,1]; uh7ma[1,1]
  # drive_through(uh7ma,uh7mae,lkm,fn,date)

  # setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/Tiedostot") #check this if needed...
  # files <- list.files(pattern="*.RData")
  # ldf <- lapply(files, load)
  # list_of_files <- list() #create empty list
  # for (i in files) {print(i); list_of_files[[i]] <- get(load(paste0("", i)))}  #add files to list position#https://www.reddit.com/r/Rlanguage/comments/nq773b/reading_multiple_rdata_files_into_a_list/
  # names(list_of_files) <- files #https://stackoverflow.com/questions/38643000/naming-list-elements-in-r
  # 
  
  # uh7f = list_of_files$`Necroinflammation uh7f.RData` #
  # uh7m = list_of_files$`Necroinflammation uh7m.RData`
  # uh7ma= list_of_files$`Necroinflammation uh7ma.RData`
  #   
  # uh7fe = list_of_files$`Necroinflammation uh7fe.RData` #e on healthy :)
  # uh7me = list_of_files$`Necroinflammation uh7me.RData`
  # uh7mae= list_of_files$`Necroinflammation uh7mae.RData`
  
  
  # date=paste(date,'_fem_yesa')
  # Group='Female';
  # allkoe=reduced2(uh7f,Group,name,lkm);
  # allkoe = na.omit(allkoe)
  # plottings_sf(allkoe,date,sick,Group)
  # 
  # date=paste(date,'_male_yes')
  # Group='Male';
  # if (sum(is.na(uh7m))==0) {allkom=reduced2(uh7m,Group,name,lkm);
  # allkom = na.omit(allkom);
  # plottings_sf(allkom,date,sick,Group)} else print('oho')
  # 
  # date=paste(date,'_all_yesa')
  # Group='All';
  # allk=reduced2(uh7ma,Group,name,lkm);
  # allk = na.omit(allk)
  # plottings_sf(allk,date,sick,Group)
  # 
  # date=paste(date,'_fem_no')
  # Group='Female';
  # allkoe=reduced2(uh7fe,Group,name,lkm);
  # allkoe = na.omit(allkoe)
  # plottings_sf(allkoe,date,sick,Group)
  # 
  # date=paste(date,'_male_no')
  # Group='Male';
  # allkom=reduced2(uh7me,Group,name,lkm);
  # allkom = na.omit(allkom)
  # plottings_sf(allkom,date,sick,Group)
  # 
  # date=paste(date,'_all_no')
  # Group='All';
  # allk=reduced2(uh7mae,Group,name,lkm);
  # allk = na.omit(allk)
  # plottings_sf(allk,date,sick,Group)
  # 
  # save.image(paste(simss,fn,date,"hypo4_noniin.RData"))
  #for the steroids:(hypoteesi3)
  # setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo3/") #https://topepo.github.io/caret/parallel-processing.html
  # hoi2=paste("C:/Users/patati/Desktop/TurkuOW/RWork/hypo3/",fn,sep='')
  # setwd(hoi2)
  # https://ds-pl-r-book.netlify.app/optimization-in-r.html
  # date='tikka19124' #change this... :)
  # t.val='no';#simss=10;
  # test=''; take=''; name=paste(simss,'hypo3_yes'); sick='yes'
  # Group='All';    uh7ma=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,date,simss,t.val,test,sick,sick_group);try({uh7ma},{uh7ma=data.frame(0)})
  # Group='Female'; uh7f=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,date,simss,t.val,test,sick,sick_group);try({uh7f},{uh7f=data.frame(0)})
  # Group='Male';   uh7m=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,date,simss,t.val,test,sick,sick_group);try({uh7m},{uh7m=data.frame(0)})
  # t.val='no';#simss=10;
  # test=''; take=''; name=paste(simss,'hypo3_no'); sick='no'
  # Group='All';    uh7mae=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,date,simss,t.val,test,sick,sick_group);try({uh7mae},{uh7mae=data.frame(0)})
  # Group='Female'; uh7fe=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,date,simss,t.val,test,sick,sick_group);try({uh7fe},{uh7fe=data.frame(0)})
  # Group='Male';   uh7me=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,date,simss,t.val,test,sick,sick_group);try({uh7me},{uh7me=data.frame(0)})
  # # setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests11/")
  # drive_through(uh7ma,uh7f,uh7m,uh7mae,uh7fe,uh7me,lkm,fn,date)
  # save.image(paste(simss,fn,date,"hypo3_all.RData"))
  }


#Testing hyp 3 and 4 with a combined function:
the_essentials=function(Treatment, Mediator, Outcome,tv, Group,name,simss,t.val,test,sick,sick_group,fn,lkm,date,joo,ip) {
  hoi1=paste("C:/Users/patati/Desktop/TurkuOW/RWork/basica/",fn,sep='')
  setwd(hoi1) #https://topepo.github.io/caret/parallel-processing.html # https://ds-pl-r-book.netlify.app/optimization-in-r.html

  name=paste(simss,'basic_yes_sick'); sick='yes'
  # Group='All'; uh7ma=loop_med_simplified1a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);try({uh7ma},{uh7ma=data.frame(0)});  save(uh7ma,file=paste(fn,'sick all.RData'));#try({uh7ma},{uh7ma=data.frame(0)})
  # # sum(tv[,2]==1 & sick_group)# Treatment, Mediator, Outcome,tv,Group,name,date,simss,t.val,test,sick,sick_group
  Group='female'; uh7f=loop_med_simplified1a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);  try({uh7f},{uh7f=data.frame(0)});save(uh7f, file=paste(fn,'sick female.RData'));
  # # uh7f=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv,Group,name,date,simss,t.val,test,sick,sick_group);try({uh7f},{uh7f=data.frame(0)})
  Group='male';   uh7m=loop_med_simplified1a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);try({res},{res=data.frame(0)});save(res, file=paste(fn,'sick male.RData'));

  # t.val='no'; #simss=10;
  # test=''; take=''; name=paste(simss,'basic_no_not sick'); sick='no'
  # Group='All';    uh7mae=loop_med_simplified1a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);save(uh7mae, file=paste(fn,'healthy alla.RData'));try({uh7mae},{uh7mae=data.frame(0)})
  # Group='female'; uh7fe=loop_med_simplified1a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);save(uh7fe, file=paste(fn,'healthy female.RData'));try({uh7fe},{uh7fe=data.frame(0)})
  # Group='male';   uh7me=loop_med_simplified1a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);save(uh7me, file=paste(fn,'healthy male.RData'));try({uh7me},{uh7me=data.frame(0)})
  # # setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/")
}


# load("100hypo3_basicah_all_18124.RData") # lkm=10; drive_through(uh7ma,uh7f,uh7m,uh7mae,uh7fe,uh7me,lkm)
drive_through=function(uh7ma,uh7mae,lkm,fn,date) {
  Group='All';name='oleelliset redusoitu mediaatiosta_sairaat_h'
  alloe=reduced2(uh7ma,Group,name,lkm);
  # Group='Female';
  # femalo=reduced3(uh7f,Group,name,lkm);
  # Group='Male';
  # malo=reduced3(uh7m,Group,name,lkm);
  # uh7mae=all_all
  
  Group='All';name='oleelliset redusoitu mediaatiosta_terveet_h'
  alloe=reduced3(uh7mae,Group,name,lkm); 
  # Group='Female'; 
  # femaloe=reduced3(uh7fe,Group,name,lkm);
  # Group='Male'; 
  # maloe=reduced3(uh7me,Group,name,lkm); 
  # rtt=rbind(allo,NA,femalo,NA,malo,NA,alloe,NA,femaloe,NA,maloe)
  rtt=rbind(allo,NA,alloe)
  
  sick='no';Group='All'; 
  plottings_sf(alloe,date,sick,Group)
  # rtt=rbind(alloe,NA,femaloe,NA,maloe)
  # rtt=rbind(alloe)
  # rownames(rtt)[rownames(rtt)=='']=c('Cond_female','Cond_male','Healthy_all','Healthy_female','Healthy_male')
  # rtt <- cbind(Row.Names = rownames(rtt), rtt)
  # rtt[rownames(rtt)=='']=c('Cond_female','Cond_male','Healthy_all','Healthy_female','Healthy_male')
  colnames(rtt)[1]='Condition_alla'
  rtt[is.na(rtt)]='empty'
  name=paste(fn,'essentials1')
  write.xlsx(rtt, file = paste(name,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  
  # Group='All';name='oleelliset redusoitu mediaatiosta_sairaat_bc_h'
  allon=reduced2(uh7ma,Group,lkm);
  # Group='Female';
  # femalo=reduced2(uh7f,Group,lkm);
  # Group='Male';
  # malo=reduced2(uh7m,Group,lkm);
  # Group='All';name='oleelliset redusoitu mediaatiosta_terveet_bc_h'
  # alloe=reduced2(uh7mae,Group,lkm); 
  # Group='Female'; 
  # femaloe=reduced2(uh7fe,Group,lkm);
  # Group='Male'; 
  # maloe=reduced2(uh7me,Group,lkm); 
  # # rtt=rbind(allo,NA,femalo,NA,malo,NA,alloe,NA,femaloe,NA,maloe)
  # rtt=rbind(alloe,NA,femaloe,NA,maloe)
  # # rownames(rtt)[rownames(rtt)=='']=c('Cond_female','Cond_male','Healthy_all','Healthy_female','Healthy_male')
  # # rtt <- cbind(Row.Names = rownames(rtt), rtt)
  # # rtt[rownames(rtt)=='']=c('Cond_female','Cond_male','Healthy_all','Healthy_female','Healthy_male')
  # colnames(rtt)[1]=paste(fn,'_all')
  # rtt[is.na(rtt)]='empty'
  # name=paste(fn,'essentials_red2'); write.xlsx(rtt, file = paste(name,date,'.xlsx'), append = FALSE, row.names = TRUE) 
  
  # sick='yes';Group='All';try({plottings_sf(alloe,date,sick,Group)})
}

reduced2=function(u3,Group,name,lkm) {
  
  # u3=jups
  u3=u3[rev(order(u3[,'ACME'])),]; 
  cut=round(quantile(abs(u3[1:200,'d0.ci_l']-u3[1:200,'z0.ci_u']),0.25),1)
  ACMEMedian=c();ACMEpval=c();ACMEVar=c()
  ADEMedian=c();ADEpval=c();ADEVar=c()
  PMEDMedian=c();PMEDpval=c();PMEDVar=c(); c1=c()
  ADEMedian=median(u3[,'ADE'])+1.5*sd(u3[,'ADE'])
  # ADEpval=quantile(u3[,'z0.p'],5) #just want to have the negative ADEs
  # DV=scale(abs(u3[,'z0.ci_l']-u3[,'z0.ci_u']))
  # ADEVar=round(quantile(DV,0.975),1)
  # min(u3[1:100,'ACME'])
  # quantile(u3[1:100,'z0.ci_l'],0.95)
  c1= u3[u3[,'ADE'] < ADEMedian  ,] #& u3[,'z0.p']<ADEpval
  
  PMEDMedian=median(c1[,'Proportion Mediated'][c1[,'Proportion Mediated']>0]) #ADE and ACME values twice in p.med?
  # PMEDpval=quantile(c1[,'n1.p'],0.95)
  # PV=scale(abs(c1[,'n.ci_l']-c1[,'n1.ci_u']))
  # PMEDVar=round(quantile(PV,0.975),1)
  c1=c1[c1[,'Proportion Mediated']>PMEDMedian  ,];
  ACMEMedian=median(c1[,'ACME'][c1[,'ACME']>0])
  ACMEpval=quantile(c1[,'d0.p'],0.5) #should be less than median
  # AV=scale(abs(c1[,'d0.ci_l']-c1[,'d0.ci_u']))
  # ACMEVar=round(quantile(AV,0.975),1)
  c1=c1[c1[,'ACME']>ACMEMedian & c1[,'d0.p']<ACMEpval,]
  # asdf=round(min(u3[1:50,'ACME'])/max(u3[1:50,'ADE']),1)# asdf=2
  cute=c()
  for (i in 1:dim(c1)[1]) {if (c1[i,'d0.p']>0.1 &  !((c1[i,'d0.ci_l']-c1[i,'z0.ci_u']) > -cut) & (c1[i,'ADE'] > -(0.1)) ) {cute=append(cute,i)} else {next}}
  c1=c1[!1:dim(c1)[1] %in% cute,] #round(quantile(u3[1:100,'z0.ci_l'],0.75),2)
  # cute2=c() # for (i in 1:dim(c1)[1]) {if (c1[i,'z0.p']<0.1 &  c1[i,'d0.p']>0.1 &  (c1[i,'ACME']-c1[i,'ADE']) < -1)  {cute2=append(cute2,i)} else {next}}
  c1=tryCatch({c1[1:lkm,]}, error = function(msg){return(c1)})
  write.xlsx(c1, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE)
  
  return(c1)}


reduced2a=function(u3,name,lkm) {
  
  # u3=list_of_files[[1]]
  u3=u3[rev(order(u3[,'ACME'])),]; 
  cut=round(quantile(abs(u3[1:200,'d0.ci_l']-u3[1:200,'z0.ci_u']),0.25),1)
  ACMEMedian=c();ACMEpval=c();ACMEVar=c()
  ADEMedian=c();ADEpval=c();ADEVar=c()
  PMEDMedian=c();PMEDpval=c();PMEDVar=c(); c1=c()
  ADEMedian=median(u3[,'ADE'])+1.5*sd(u3[,'ADE'])
  c1=u3
  # c1= u3[u3[,'ADE'] < ADEMedian  ,] #& u3[,'z0.p']<ADEpval
  PMEDMedian=median(c1[,'Proportion Mediated'][c1[,'Proportion Mediated']>0]) #ADE and ACME values twice in p.med?
  c1=c1[c1[,'Proportion Mediated']>PMEDMedian  ,];
  ACMEMedian=median(c1[,'ACME'][c1[,'ACME']>0])
  ACMEpval=quantile(c1[,'d0.p'],0.1) #should be less than median
  # kit=abs(c1[i,'d0.ci_l'])/(c1[i,'d0.ci_u']+abs(c1[i,'d0.ci_l']))
  
  c1=c1[c1[,'ACME']>0,]
  cute=c()
  for (i in 1:dim(c1)[1]) {if ((c1[i,'d0.ci_l']) <= -0.1 &  
          abs(c1[i,'d0.ci_l'])/(c1[i,'d0.ci_u']+abs(c1[i,'d0.ci_l'])) > 0.1 & c1[i,'d0.p'] >= ACMEpval &  c1[i,'ADE'] > -0.1 ) {cute=append(cute,i)} else {next}} #& c1[i,'z0.p'] < ACMEpval & c1[i,'d0.ci_l'] < c1[i,'ADE']
  c1=c1[!1:dim(c1)[1] %in% cute,] #round(quantile(u3[1:100,'z0.ci_l'],0.75),2)
  
  cute2=c()
  for (i in 1:dim(c1)[1]) {if ((c1[i,'d0.ci_l'] < c1[i,'z0.ci_u']) &  abs(c1[i,'z0.ci_u'])/(c1[i,'z0.ci_u']+abs(c1[i,'z0.ci_l'])) > 0.5 &
      abs(c1[i,'d0.ci_l'])/(c1[i,'d0.ci_u']+abs(c1[i,'d0.ci_l'])) > 0.1 & c1[i,'d0.p'] >= ACMEpval &  c1[i,'ADE'] > c1[i,'ACME'] ) {cute2=append(cute2,i)} else {next}} #& c1[i,'z0.p'] < ACMEpval & c1[i,'d0.ci_l'] < c1[i,'ADE']
  c1=c1[!1:dim(c1)[1] %in% cute2,] 
  
  
  c1=tryCatch({c1[1:lkm,]}, error = function(msg){return(c1)})

  return(c1)}


# 
# reduced3=function(u3,Group,name) {
#   ACMEMedian=c();ACMEpval=c();ACMEVar=c()
#   ADEMedian=c();ADEpval=c();ADEVar=c()
#   PMEDMedian=c();PMEDpval=c();PMEDVar=c()
#   c1=c(); c1=u3
#   ACMEMedian=median(c1[,'ACME'][c1[,'ACME']>0])
#   ACMEpval=quantile(c1[,'d0.p'],0.5) #should be less than median
#   AV=scale(abs(c1[,'d0.ci_l']-c1[,'d0.ci_u']))
#   ACMEVar=round(quantile(AV,0.95),1)
#   c1=c1[c1[,'ACME']>ACMEMedian & c1[,'d0.p']<ACMEpval & AV<ACMEVar,] 
#   io=c(); for (i in 1:length(rownames(c1))) {io=append(io,grepl( '_L', rownames(c1)[i], fixed = TRUE))}
#   c1=c1[!io,]
#   write.xlsx(c1, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
#   return(c1)}


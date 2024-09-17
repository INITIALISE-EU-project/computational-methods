# Codes starting from 8.8.23. Pauli Tikka. 
# Updated 27.05.2024 (if in conflict, see the save date of the file)

#Installations:
# remotes::install_github("davidsjoberg/ggsankey");BiocManager::install("qpgraph")
# devtools::install_github("NightingaleHealth/ggforestplot") ...
# this works quite often: install.packages("draw", type="source")
#DONE (ggsankey, monocle3 and others) :)
# length(6:80) # The libraries:
library(scales);library("plotrix");library(Maaslin2);  library(ggtext);library(lavaan);library(psych);library("xlsx");library(lsr);library(quantreg);library("readxl")
library(semPlot);library(mediation);require(lme4);require(reshape2);library(mlma);library(binilib);library(plyr);library("viridis");library(RColorBrewer) #library(mma);# library(d3heatmap);
library(magrittr); library(ggplot2);library( "censReg" );library(ggsankey);library(dplyr);library(tidyverse);library(dmetar);library(meta);library(ggforestplot); ## do not add: library(forestplot!!!not working with ggforestplot..);
library(mdatools);library(circlize);library(igraph);library('bigsnpr');library(rcompanion);library(scRNAseq);library(tibble);library(stringr);library(MOFA2);library('qpgraph') #ok
library("grid"); library("ggplotify");library(ggpubr);library(rstatix);library(datarium);library(RColorBrewer); library(ggh4x); library(effsize)
library(chorddiag);library(corrplot);library(scater);library(mdatools);options(scipen = 999); library(car);library(FSA);library(pathviewr);library(glmnet)
library("lmtest");library(PerformanceAnalytics);library(psych);library("readxl");library(ggforce);library(ComplexHeatmap) #these are ok to drive in start
library('Hmisc');library(correlation);library(ggcorrplot);library(pheatmap);library(mgcv);library('ppcor')#install.packages("ppcor")
library(extrafont)
font_import() #this is important
loadfonts(device = "win") #this is important too
library(rmdformats);library(prettydoc);library(hrbrthemes);library(tint);library(tufte)

# install.packages(c('rmdformats','prettydoc','tint','tufte'))

#Correlations can be found from line 660 (17.5.24)

# library(ComplexHeatmap);# library("heatmaply")
setwd("C:/Users/patati/Desktop/TurkuOW/RWork") 
#jotkin vanhoista paketeista (esim. SparkL) sattaaa hidastaa
date='tikka9924' #change this...
# saving all data to the path# save.image("saveworkspace_boot3_all.RData") #"saveworkspace_2000.RData" "saveworkspace.RData"
# save.image("huutumahe3924.RData")# save.image("basic2e.RData") #ilman x:iä
# loading the workspace
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork")
# load("huutumahe13824.RData") #xien kaa

# #Trying to especially reproduced foresplots as in # https://www.nature.com/articles/s41467-022-30875-7#citeas (fig. 5) for Steroid data:
# Steroid data, simple figures:# NAFLD=read.csv("NAFLD_SteroidStudy_v3.csv", header = TRUE, sep=";")
NAFLD=read_excel("NAFLD_SteroidStudy.xlsx",sheet = "LFAT_steroidsDATA") #l ei tästä
oknames=colnames(NAFLD); NAFLD=data.frame(NAFLD)
# rownames(NAFLD)=NAFLD[,1]
groups=read.csv("groups_17823.csv", header = TRUE, sep=";")
groups=groups[,c('Group','Abbreviation')]
groups=groups[groups[,'Abbreviation']!='F',]
groups=groups[order(groups[,'Group']),]

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
# MASLD[1:30,1:12]
# NAFLD[1:30,7:20]
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
# tve=tv[,9:dim(tv)[2]]; tve[tve == 0] <- NA; #print(tve, max=nrow(tve)*ncol(tve)); note, here the covariates have not been normalized or scaled/elaborated in any way; maybe I need to do so (28524...)
tve=tv[,2:dim(tv)[2]]; tve[tve == 0] <- NA; #print(tve, max=nrow(tve)*ncol(tve)); note, here the covariates have not been normalized or scaled/elaborated in any way; maybe I need to do so (28524...)
tv_half <- tve %>% mutate(replace(., is.na(.), min(., na.rm = T)/2)) #https://mdatools.com/docs/preprocessing--autoscaling.html
tv_half_log2 <- log2(tv_half);
# print(tv_half_log2, max=nrow(tv_half_log2)*ncol(tv_half_log2))
tv_auto <- prep.autoscale(tv_half_log2, center = TRUE, scale = TRUE); 
head(tv_auto) #non nans 
which(is.na(tv_auto))

# tv_auto[1:5,1:11] 
#usually this should be the log2 value 'tv_half_log2' & #https://svkucheryavski.gitbooks.io/mdatools/content/preprocessing/text.html
# Necroinflammation  HOMA-IR Steatosis.Grade.0.To.3 Fibrosis.Stage.0.to.4
tv_all=cbind(tv[,1],tv_auto); #tv_all[1:5,1:11]; note, here the covariates have not been normalized or scaled/elaborated in any way;  maybe I need to do so (1/324 or 28524...); check 27524 the notation: tv_all=cbind(tv[,1:8],tv_auto);
# tv_all=cbind(tv[,1:8],tv_auto); #tv_all[1:5,1:11]; note, here the covariates have not been normalized or scaled/elaborated in any way;  maybe I need to do so (1/324 or 28524...); check 27524 the notation: tv_all=cbind(tv[,1:8],tv_auto); 

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

#jälkeenpäin lienee jeesh
# x5=x5[x5!='PFAS'];
# # x5=x5[x5!='Perfluorodecyl.ethanoic.acid'];
# x6=x6[x6!='Total_TG'] # x1;x2;x3;x4;x5;
# tv_all=tv_all[,!colnames(tv_all) %in% c('Total_TG','PFAS')]

tv_half_log22=cbind(tv[,1],tv_half_log2);
# tv_half_log22=cbind(tv[,1:8],tv_half_log2);
x1=colnames(tv_half_log22[,c(1:8)]); v2=dim(NAFLD)[2]+1
x2=colnames(tv_half_log22[,9:v2]);v3=(dim(Bali)[2]+v2);
x3=colnames(tv_half_log22[,(v2+1):(v3)]);v4=(dim(Base)[2])+v3
x3=x3[c(length(x3),1:(length(x3)-1))]
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
# tv_half_log22=tv_half_log22[,!colnames(tv_half_log22) %in% c('Total_TG','PFAS')]

colnames(tv)[colnames(tv)=='17aOH-P4']='17a-OHP4'
colnames(tv_half_log22)[colnames(tv_half_log22)=='17aOH-P4']='17a-OHP4'
colnames(tv_all)[colnames(tv_all)=='17aOH-P4']='17a-OHP4'

Treatment=colnames(tv_all)[71:77];
Mediator=colnames(tv_all)[9:28];
Outcome=colnames(tv_all)[c(29:51,78:90)]; ##https://sparkbyexamples.com/r-programming/r-remove-from-vector-with-examples/
Treatment=Treatment[!Treatment %in% c('Perfluorodecyl.ethanoic.acid')]
tv_all=tv_all[,!colnames(tv_all) %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]
tv_all=tv_all[,!colnames(tv_all) %in% x4]
# tv_all=tv_all[,!colnames(tv_all) %in% c('TG_SFA','MUFA','TG_PUFA')]
Outcome=Outcome[!Outcome %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]
Outcome=Outcome[! Outcome %in% x4] #https://sparkbyexamples.com/r-programming/r-remove-from-vector-with-examples/
# Outcome=Outcome[! Outcome %in% c('TG_SFA','MUFA','TG_PUFA')] #
Mediator[Mediator=="17aOH-P4"]="17a-OHP4"

tv_covscl=tv_all
tv_covNS=cbind(tv[,1:8],tv_all[,9:dim(tv_all)[2]])
tv_LOG_covscl=tv_half_log22
tv_LOG_covNS=cbind(tv[,1:8],tv_half_log22[,9:dim(tv_half_log22)[2]])

colnames(tv_covNS)[1:8]=colnames(tv_all)[1:8]
colnames(tv_LOG_covNS)[1:8]=colnames(tv_all)[1:8]

# ...

# This is how you do it with one or two cases:# print(a)+print(b)
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/forest.html
# https://cran.r-project.org/web/packages/bigsnpr/index.html
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula
# https://stats.stackexchange.com/questions/190763/how-to-decide-which-glm-family-to-use
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm
# https://www.r-bloggers.com/2015/09/how-to-perform-a-logistic-regression-in-r/
# Others: https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/forest.html
# https://www.khstats.com/blog/forest-plots/
# https://www.geeksforgeeks.org/how-to-create-a-forest-plot-in-r/ 
# https://whitlockschluter3e.zoology.ubc.ca/RLabs/R_tutorial_Contingency_analysis.html
# Outcome: Necroinflammation  HOMA-IR Steatosis.Grade.0.To.3 Fibrosis.Stage.0.to.4
# Out:'Steatosis' 'HOMA-IR' 'Necroinflammation' 'Fibrosis'
# Fibrosis.Stage.0.to.4 Steatosis.Grade.0.To.3 Fibrosis Stage
# NAFLD=cbind(tv[,1:28]) #these are raw, not log nor centered values, such as tv_half_log2[,1:20]) # -0.19608540/-0.14242200 

#You need to do this to get rid of Xs and other not recognizable for the function, in the case of reduced number of values, use ie instead of tv for 'NAFLD':
NAFLD=cbind(tv[,1:28])#,tv_half_log2[,1:20])#tv[,1:28]); #ie# tv_half_log22
colnames(NAFLD) #autoscaled; # NAFLD=cbind(tv[,1:8],tv_half_log2[,1:20]) #for raw # NAFLD=cbind(tv[,1:28])
NAFLD[NAFLD[,c(5)]>0,5]=1;NAFLD[NAFLD[,c(6)]>0,6]=1;NAFLD[NAFLD[,c(7)]>0,7]=1;
NAFLD[NAFLD[,c(8)] <= 1.5,8]=0;NAFLD[NAFLD[,c(8)]>1.5,8]=1; #note the order of calculation...median
colnames(NAFLD) <- gsub("-", ".", colnames(NAFLD))
colnames(NAFLD) <- gsub("/", ".", colnames(NAFLD))
colnames(NAFLD) <- gsub("11", "X11", colnames(NAFLD))
colnames(NAFLD) <- gsub("17", "X17", colnames(NAFLD))
colnames(NAFLD) <- gsub("#", ".", colnames(NAFLD))
colnames(NAFLD)[colnames(NAFLD)=='X17aOH.P4']='X17.aOHP4' #oh...

#this is with first(!!), use it
#is normal (non-log scale) image (and autoscaled), loge=1 is  log10 scale image (and raw data)
# Necroinflammation  HOMA-IR Steatosis.Grade.0.To.3 Fibrosis.Stage.0.to.4
date='22224'; loge=0; #xlim=c(0.1,1.7)
Outcome='Steatosis.Grade.0.To.3';Out='Steatosis'; oute='Steatosis';
first=TRUE; e='P4';ordera=c();Group='All';
name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); 
hel=pre_errors_2(NAFLD,Outcome,Group,name,ordera,oute,first,e,xlim) 
#Afterwards:
first=FALSE;
Group='Female';name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); 
pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim)
Group='Male'; name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim)
#Afterwards:
first=FALSE;
Outcome='Fibrosis.Stage.0.to.4'; Out='Fibrosis';oute='Fibrosis'; 
Group='All'; name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date);pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim) #not the very first though...
Group='Female';
name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim)
Group='Male'; name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim)

Outcome='Necroinflammation'; Out='Necroinflammation';oute='Necroinflammation'; 
Group='All'; 
name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim) #not the very first though...
Group='Female';name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim)
Group='Male'; name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim)

Outcome='HOMA.IR';Out='HOMA-IR';oute='HOMAIR';
# xlim=c(0.1,5.3)
Group='All';name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date);pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim) #not the very first though...
# xlim=c(0.1,1.7)
Group='Female';name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim)
Group='Male'; name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); pre_errors_2(NAFLD,Outcome,Group,name,ordera=hel,oute,first,e,xlim)

# hist(df[df[,1]==e,'result.1'],breaks=50); hist(SG0[,'P4'],breaks=50) #colnames(NAFLDo)[9:28][ps<0.05]
# https://bookdown.org/content/b472c7b3-ede5-40f0-9677-75c3704c7e5c/more-than-one-mediator.html
#This works with the autoscaled (raw if loge=1 and remove 1 in the means) data NAFLD...
pre_errors_2=function(NAFLD,Outcome,Group,name,ordera,oute,first,e,xlim) { # Group='Female'
  
  if (Group=='Male') {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==2,]} else if (Group=='Female') 
  {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==1,]} else if (Group=='All') {NAFLDo=NAFLD} 
  # e=e;# #  # hist(NAFLDo[,e],breaks=50)# mean(NAFLDo[,e])
  sample_data=c();n0=c();n1=c()
  # shapie=c();for (j in 9:28) {shapie=append(shapie,shapiro.test(NAFLDo[,j])$p.value)}; sum(shapie<0.05) 
  # medians ok; sum(shapie<0.05) = 17 and 18 (for males and females)
  
  for (i in 1:2) {
    # i=2
    if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] == 0,];n0=dim(SG0)[1]} else if (i==2) {SG0=NAFLDo[NAFLDo[,Outcome] > 0,];n1=dim(SG0)[1]}#Steatosis.Grade.0.To.3 Fibrosis.Stage.0.to.4
    # hist(SG0[,e],breaks=50);print(shapiro.test(SG0[,e])) #ok in A4 but not in P4 in all cases and controls for both female and male + we have definately outliers ad described above..
    #https://stats.stackexchange.com/questions/237256/is-the-shapiro-wilk-test-only-applicable-to-smaller-sample-sizes # https://statisticsbyjim.com/hypothesis-testing/nonparametric-parametric-tests/
    means=c();for (j in 9:28) {means=append(means,median(SG0[,j], na.rm=TRUE))} #
    # shapi=c();for (j in 9:28) {shapi=append(shapi,shapiro.test(SG0[,j])$p.value)};sum(shapi<0.05)
    #medians ok; sum(shapi<0.05) = 14/20 not normal (for female it is 18!)
    # shapi2=c();for (j in 9:28) {shapi2=append(shapi2,shapiro.test(SG0[,j])$p.value)};sum(shapi2<0.05) 
    #sum(shapi2<0.05) = 7/20 not normal (for females it is 8!)
    # means[is.na(means)]=quantile(means[!is.na(means)],0.25) #low number of males in steatosis!
    #https://www.statology.org/mean-standard-deviation-grouped-data/ # https://amsi.org.au/ESA_Senior_Years/SeniorTopic4/4h/4h_2content_11.html # https://www.themathdoctors.org/mean-and-standard-deviation-of-grouped-data/
    sds=c();for (j in 9:28) {sds=append(sds,sd(SG0[,j],na.rm=TRUE))} #tähän jäätiin...
    # sds[is.na(sds)]=quantile(sds[!is.na(sds)],0.25) #low number of males in steatosis! and necroinfl.
    error_lower=means-sds# error_lower[error_lower <= 0] = 0
    error_upper=means+sds; error=sds
    sample_data <- append(sample_data,data.frame(study=colnames(NAFLD[,9:28]),index=colnames(NAFLD[,9:28]),result=means,error=error))} # cate<- (str_extract(colnames(SG0[,9:28]), "[aA-zZ]+")) #https://stackoverflow.com/questions/29825537/group-categories-in-r-according-to-first-letters-of-a-string
    # https://datatofish.com/create-dataframe-in-r/#if you know: ga=groups[,'Abbreviation']; md=metad[,'name']; unique(c(ga,md)); ga[!(ga %in% md)]; ga[ga=="17a-OHP4"]="17aOH-P4"
  # https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot
  # http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  # https://blogs.sas.com/content/iml/2011/04/27/log-transformations-how-to-handle-negative-data-values.html
  # https://stats.stackexchange.com/questions/155429/how-to-transform-negative-values-to-logarithms
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1534033/
  df=data.frame(sample_data) #
  ps=c();for (j in 9:28) {xnam <- colnames(NAFLDo)[j]; fmla <- as.formula(paste(xnam, "~",Outcome));
  ps=append(ps,wilcox.test(fmla, data = NAFLDo,exact = FALSE)$p.value)}#kruskal.test
  # psa=c();for (j in 9:28) {xnam <- colnames(NAFLDo)[j]; fmla <- as.formula(paste(xnam, "~",Outcome));psa=append(psa,kruskal.test(fmla, data = NAFLDo)$p.value)}#kruskal.test, http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  # > ps
  # [1] 0.609954658 0.715570972 0.537679346 0.549032904 0.181524372 0.138385547 0.548245444 0.176248250 0.253563009 0.621346211 0.086042142 0.192710546 0.081028520 0.076172801 0.004956634 0.001198884 0.007008253
  # [18] 0.070960960 0.286592996 0.307597411# > psa
  # [1] 0.607121328 0.712550297 0.535007089 0.546336641 0.180204741 0.137311216 0.544330463 0.174957178 0.251883482 0.618489657 0.085304517 0.190922561 0.080324413 0.075504719 0.004894636 0.001181976 0.006923521
  # [18] 0.070330542 0.284764613 0.305681757
  print(mean(df[df[,1]==e,'result.1']));#hist(df[df[,1]==e,'result.1'],breaks=50)
  print(mean(df[df[,1]==e,'result']));#hist(df[df[,1]==e,'result'],breaks=50)
  # https://www.statisticshowto.com/probability-and-statistics/statistics-definitions/parametric-and-non-parametric-data/
  a=df[df[,1]==e,'result.1']/df[df[,1]==e,'result']; 
  print(a); print(log(a))
  v2=data.frame(log(df$result.1/df$result)) #this should give the right order of the variables if not the absolute change
  v2[,'result']=v2[,1];v2[,'name']=df$study;v2=v2[,2:3]
  v2[,'name'] <- gsub("\\.", "-", v2[,'name']) #https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot
  v2[,'name'] <- gsub("X11", "11", v2[,'name'])
  v2[,'name'] <- gsub("X17", "17", v2[,'name'])
  v2[,'name'][v2[,'name']=="T-Epi-T"]="T/Epi-T"
  v2[,'pval']=ps# df2=df[rev(order(df[,'error'])),]
  # https://www.bmj.com/content/312/7038/1079.full# https://stats.stackexchange.com/questions/589920/how-can-i-back-transform-a-log-data-to-interpret-t-test-and-get-original-ci
  # https://www.biostars.org/p/16481/
  # https://whitlockschluter3e.zoology.ubc.ca/RLabs/R_tutorial_Contingency_analysis.html # https://sphweb.bumc.bu.edu/otlt/mph-modules/ep/ep713_randomerror/ep713_randomerror6.html
  # https://www.r-bloggers.com/2015/01/easy-error-propagation-in-r/# https://www.biostars.org/p/342756/ https://en.wikipedia.org/wiki/Tukey%27s_range_test 
  #https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test
  v2[,'result_pure']=(df$result.1/df$result) #case and control
  # dy=df$error.1/n1*1.64;dx=df$error/n0*1.64 #check the numbers
  # v2[,'error']=v2$result_pure*sqrt((dy/df$result.1)^2+(dx/df$result)^2) #v2$result*(sqrt((df$error.1/as.numeric(df1[ax,'result']))^2+(df$error/as.numeric(df0[ix,'result'])^2)))
  v2[,'error']=(abs((1/df$result)*df$error.1)+abs((df$result.1/df$result^2)*df$error))/dim(NAFLDo)[1]*1.64
  # v2[,'error']=dim(NAFLDo)[1]
  
  #This error above is ok... it is from the Jukka Vaari 1993 Fysiikan Laboratoriotyöt  
  # v2[,'error2']=v2[,'error']/sqrt(dim(NAFLDo)[1])*1.64#1.64
  v2[,'error'][v2[,'error']>(median(v2[,'error'])+sd(v2[,'error']))]=median(v2[,'error'])*1.25
  v2[,'errord1a']=v2[,'result_pure']-v2[,'error']#(df$result.1-dy)/(df$result+dx)# #
  v2[,'errord2a']=v2[,'result_pure']+v2[,'error']#(df$result.1+dy)/(df$result-dx)##
  # v2$errord1[v2$errord1 == '-Inf']=min(v2$errord1[v2$errord1 != 'NaN'])*1.05
  # v2$errord1[order(v2$errord1)]<0
  # v2[,'errord1a'][v2[,'errord1a']<0]=0#quantile(v2[,'errord1a'][v2[,'errord1a']0],0.1)
  # v2[,'errord2a'][v2[,'errord2a']<0]=0.1#quantile(v2[,'errord2a'][v2[,'errord2a']>0],0.1)
  # v2$errord1a[v2$errord1a<0]=v2$errord1a[order(v2$errord1a)][2]*0.95
  # abs(max(v2[,'error'])-median(v2[,'error']))
  # abs(v2$errord1a-v2$errord2a)
  v2[,'errord1']=log(v2[,'errord1a'])
  v2[,'errord2']=log(v2[,'errord2a'])
  v2[,'result']=log(v2[,'result_pure'])
  v2[,'Control']=df$result
  v2[,'Case']=df$result.1 #
  #https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot

  v2[,'pval0']=v2[,'pval']#fp$pval
  v2[,'pval1']=v2[,'pval']#fp$pval # write.csv(v2n, 'scaled_fibrosis results.csv')# write.csv(v2c, 'raw_fibrosis results.csv')# v2$pval0[order(v2$pval0)]
  v2[,'Significance0']= v2[,'pval0']<0.1#(v2[,case1]  < 1 & v2[,case2] < 1) | (v2[,case1]  >= 1 & v2[,case2] >= 1)
  v2[,'Significance0'][v2[,'Significance0']==TRUE]='Yes'
  v2[,'Significance0'][v2[,'Significance0']==FALSE]='No'
  v2[,'Color0']=v2[,'pval0'] < 0.1#(v2[,case1]  < 1 & v2[,case2] < 1) | (v2[,case1]  >= 1 & v2[,case2] >= 1)
  v2[,'Color0'][v2[,'Color0']==TRUE]='blue'
  v2[,'Color0'][v2[,'Color0']==FALSE]='grey'
  v2[,'Significance1']= v2[,'pval1']<0.1#(v2[,case1]  < 1 & v2[,case2] < 1) | (v2[,case1]  >= 1 & v2[,case2] >= 1)
  v2[,'Significance1'][v2[,'Significance1']==TRUE]='Yes'
  v2[,'Significance1'][v2[,'Significance1']==FALSE]='No'
  v2[,'Color1']=v2[,'pval1'] < 0.1#(v2[,case1]  < 1 & v2[,case2] < 1) | (v2[,case1]  >= 1 & v2[,case2] >= 1)
  v2[,'Color1'][v2[,'Color1']==TRUE]='blue'
  v2[,'Color1'][v2[,'Color1']==FALSE]='grey'

  gn=groups[,c('Group','Abbreviation')]
  gn=gn[gn[,'Abbreviation']!='F',]
  gn=gn[order(data.frame(gn[,'Abbreviation'])[,1]),]
  v2=v2[order(v2[,'name']),]
  v2=cbind(v2,gn[order(data.frame(gn[,'Abbreviation'])[,1]),])
  v2=v2[rev(order(v2[,'result'])),]
  # loge=1; # if (round(max(v2$errord2)*1.0,2) <= 0) {eo=1} else if (round(max(v2$errord2)*1.0,2)>0) {eo=round(max(v2$errord2)*1.05,2)+0.1}
  # if (eo<1) {eo=1} # xlim = c(round(min(v2$errord1),2), eo); 
  xlab = "Autoscaled Concentrations (SE)" #xlab = "Raw Concentrations in Log10 Scale (SE)"}
  # https://phas.ubc.ca/~oser/p509/Lec_10.pdf
  # xlim=xlim
  # xlim=c(min(v2$result)-max(abs(v2$errord1-v2$errord2)/2)*1.13,max(v2$result)+max(abs(v2$errord1-v2$errord2)/2)*0.9)
  xlim=c(min(v2$errord1),max(v2$errord2))
  # xlim=c(min(v2$result)*1.1,max(v2$result)*1.1) # if (xlim[2]>1) {xlim[2]=1};# if (xlim[1] < -0.75) {xlim[1]=-0.75};
  plote2=forestplot(df = v2, #drive tba_example_v3_oh_tikka17823.R if not working via 'x' error #coef, #
                    estimate = result, 
                    # lower=v2$errord1,
                    # upper=v2$errord2,
                    # xmin=v2$errord1,
                    # xmax=v2$errord2,
                    se=0,#abs(errord1-errord2)/4, #sterr,##this makes the significant value:
                    pvalue = pval1,psignif = 0.1,
                    xlim=xlim, xlab = 'Logged Ratio between Raw Concentrations of Case and Control with 90% CI',ylab='Steroid Groups',
                    title='',colour = Significance1 ) +#,colour = Significance 
    ggforce::facet_col(facets = ~Group,scales = "free_y",space = "free", strip.position='left')+
  geom_errorbarh(aes(xmin = errord1, xmax = errord2,height = .0,colour=Significance1));#plote2
  #+aes(x = result, xmin = errord1, xmax = errord1)#Space was free
      
  # ordera=hel
  
  if (sum(v2[,'Significance1']=='Yes')==20) {hp=c('blue','blue')} else {hp=c('#999999','blue')};#plote2
  if (Group=='All' & first==TRUE) {ordera=v2$name[order(v2$result)]; #
  plote2[["data"]][["name"]]=factor(plote2[["data"]][["name"]], levels = ordera)} else if
  (Group=='All' & first==FALSE) {plote2[["data"]][["name"]]=factor(plote2[["data"]][["name"]], levels = ordera)} else if
  (Group=='Female') {plote2[["data"]][["name"]]=factor(plote2[["data"]][["name"]], levels = ordera)} else if
  (Group=='Male') {plote2[["data"]][["name"]]=factor(plote2[["data"]][["name"]], levels = ordera)}
  #https://www.r-bloggers.com/2020/03/how-to-standardize-group-colors-in-data-visualizations-in-r/
  plote2$layers[[1]]$aes_params$odd <- "#00000000" #https://stackoverflow.com/questions/71745719/how-to-control-stripe-transparency-using-ggforestplot-geom-stripes
  v2$Group2=v2$Group
  v2 <- transform(v2,Group2 = as.numeric(as.factor(Group2)))
  v2$facet_fill_color <- c("red", "green", "blue", "yellow", "brown")[v2$Group2]
  jopon=plote2  +theme(axis.text.y=element_blank()) +theme_classic2();   #theme(axis.text.y = element_text(lineheight=.05));
  jopon2=jopon+geom_point(aes(colour = factor(Significance1)),colour = v2[,'Color1']) +  
    scale_color_manual(values=hp)+theme(legend.position = "none")+theme(strip.text.y = element_text(size=-Inf)) #ggtext::element_markdown(size = 12)
  # geom_vline(xintercept = 1);
  # https://stackoverflow.com/questions/10547487/remove-facet-wrap-labels-completely
  g <- ggplot_gtable(ggplot_build(jopon2))
  stripr <- which(grepl('strip-l', g$layout$name)); fills <- c("red","green","blue","yellow",'brown'); k <- 1; 
  for (i in stripr) {j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]; k <- k+1}
  grid::grid.draw(g)
  #https://stackoverflow.com/questions/24169675/multiple-colors-in-a-facet-strip-background-in-ggplot
  #http://www.sthda.com/english/wiki/ggplot2-axis-scales-and-transformations #https://www.statology.org/geom_point-fill/
  # https://stackoverflow.com/questions/62093084/set-geom-vline-line-types-and-sizes-with-aes-mapping-in-ggplot2
  # https://www.google.com/search?q=calculate+standard+error&sca_esv=564912168&rlz=1C1GCEA_enFI1072FI1072&ei=qjsBZYbOL-WMxc8PjumCkAg&ved=0ahUKEwiGy8qF4aaBAxVlRvEDHY60AIIQ4dUDCBA&uact=5&oq=calculate+standard+error&gs_lp=Egxnd3Mtd2l6LXNlcnAiGGNhbGN1bGF0ZSBzdGFuZGFyZCBlcnJvcjIKEAAYRxjWBBiwAzIKEAAYRxjWBBiwAzIKEAAYRxjWBBiwAzIKEAAYRxjWBBiwAzIKEAAYRxjWBBiwAzIKEAAYRxjWBBiwAzIKEAAYRxjWBBiwAzIKEAAYRxjWBBiwAzIKEAAYigUYsAMYQzIKEAAYigUYsAMYQ0jVRFDZEFjXIHACeAGQAQCYAUSgAZMEqgEBObgBA8gBAPgBAcICBhAAGAcYHsICCBAAGAcYHhgKwgIIEAAYBxgeGBPCAgoQABgHGB4YExgK4gMEGAAgQeIDBRIBMSBAiAYBkAYK&sclient=gws-wiz-serp
  jpeg(paste(name ,"divi.jpg"), width = 7500, height = 11000, quality = 100,pointsize = 16, res=1000); print(grid::grid.draw(g));dev.off();
  return(ordera) } 

# https://csass.ucsc.edu/certification/power_goldman.pdf
library(pwr) #The reason why I used non-parametric test:
hist(unlist(e),breaks=50) #e from the effect size calculations
most_min=-0.65; qm=quantile(unlist(e),0.25); mean=mean(unlist(e))
pwr.t2n.test(n1 =69 , n2= 35, d = most_min, sig.level =0.05) #0.873409
pwr.t2n.test(n1 =69 , n2= 35, d = qm, sig.level =0.05)       #0.348547
pwr.t2n.test(n1 =69 , n2= 35, d = mean, sig.level =0.05)     #0.1804521
#All steatosis
pwr.t2n.test(n1 =dim(SG0)[1] , n2= dim(SG1)[1], d = qm, sig.level =0.05) #0.2661434, dim(SG0)[1];dim(SG1)[1] :83 21
# Female (stetosis)
pwr.t2n.test(n1 =dim(SG0)[1] , n2= dim(SG1)[1], d = qm, sig.level =0.05) #0.1666289 # ddim(SG0)[1];dim(SG1)[1] 58  11
#Male (steatosis)
pwr.t2n.test(n1 =dim(SG0)[1] , n2= dim(SG1)[1], d = qm, sig.level =0.05) #0.1367615 # ddim(SG0)[1];dim(SG1)[1] 25  10
# https://www.statmethods.net/stats/power.html #https://www.scribbr.com/statistics/statistical-power/
# https://www.google.com/search?q=good+power+test+amount&rlz=1C1GCEA_enFI1072FI1072&oq=good+power+test+amount&gs_lcrp=EgZjaHJvbWUyBggAEEUYOTIHCAEQIRigAdIBCDY0MDBqMGo3qAIAsAIA&sourceid=chrome&ie=UTF-8
# https://www.moresteam.com/whitepapers/download/power-stat-test.pdf (80 is about the limit)->
#i.e. use non-parametric tests..
# https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS6.html
# https://www.sciencedirect.com/science/article/pii/S200103701930025X


#This may need the final scaled values (tikka 28524...)
#The circular/chord plots: #The values for the chord plot:
double_circular=function(tv_covscl,n_level,date,x1,x2,x3,x4,x5,x6) {
tv_c=tv_covscl#cbind(tv[,1:8], tv_half_log2) #check also not logged and then the auto one
# x1=colnames(tv_c[,c(1:8)]); v2=dim(NAFLD)[2]+1
# x2=colnames(tv_c[,9:v2]);v3=(dim(Bali)[2]+v2)
# x3=colnames(tv_c[,(v2+1):(v3)]);v4=(dim(Base)[2])+v3
# x4=colnames(tv_c[,(v3+1):(v4-1)])
# x5=colnames(tv_c[,(v4):(dim(tv_c)[2])]); 
# x5a=x5[1:9]
# x6=x5[10:length(x5)] #dividing to lipids
# x5=x5a  #making sure that PFAS are separate
# x1;x2;x3;x4;x5; # tv_auto[1:5,1:11]
# x3 <- paste(x3, "_L", sep="") #https://stackoverflow.com/questions/6984796/how-to-paste-a-string-on-each-element-of-a-vector-of-strings-using-apply-in-r
# x4=gsub("(-[0-9]*)*.1", "", x4) #https://stackoverflow.com/questions/18997297/remove-ending-of-string-with-gsub
# x4 <- paste(x4, "_S", sep="")# https://rdrr.io/bioc/qpgraph/man/qpNrr.html
# nm = c(x1,x2,x3,x4,x5,x6) # adding the lipids, x6
# nm=c('Subject#','Gender','AGE','BMI','Steatosis Grade','Fibrosis Stage','Necroinflammation','HOMA-IR',nm[9:91])
# colnames(tv_c)=nm # tv_c[1:5,1:30]; NAFLD[1:2,1:28];
# colnames(tv_c)[colnames(tv_c)=='MASLD[, \"PFAS\"]']='PFAS';
#removing PFAS and lipid generals:
# x5=x5[x5!='PFAS'] ; x6=x6[x6!='Total_TG'] # x1;x2;x3;x4;x5;
tv_c=data.frame(tv_c)

# tv_c=tv_c[,3:dim(tv_c)[2]]

tv_c=tv_c[,!colnames(tv_c) %in% c('Total_TG','PFAS',"Perfluorodecyl.ethanoic.acid")]
tvf=tv_c[tv_c[,'Gender']==min(tv_c[,'Gender']),1:dim(tv_c)[2]]
# tv['Steatosis.Grade.0.To.3'==0,9:27]] #tv[tv[,'Necroinflammation']==0,9:80]; #SG0i=as.numeric(SG0i); check also: tv[tv[,'HOMA-IR']==0,9:80]
tvm=tv_c[tv_c[,'Gender']==max(tv_c[,'Gender']),1:dim(tv_c)[2]]

# ;tvf=data.frame(tvf);tvm=data.frame(tvm);

tvtest=list(tv_c,tvf,tvm)
for (i in 1:3) {
  colnames(tvtest[[i]]) <- gsub("\\.", "-", colnames(tvtest[[i]]))
colnames(tvtest[[i]]) <- gsub("X11", "11", colnames(tvtest[[i]]))
colnames(tvtest[[i]]) <- gsub("X17", "17", colnames(tvtest[[i]]))
colnames(tvtest[[i]])[colnames(tvtest[[i]])=="T-Epi-T"]="T/Epi-T"
colnames(tvtest[[i]])[colnames(tvtest[[i]])=="Steatosis-Grade"]="Steatosis Grade"
colnames(tvtest[[i]])[colnames(tvtest[[i]])=="Fibrosis-Stage"]="Fibrosis Stage"
colnames(tvtest[[i]])[colnames(tvtest[[i]])=="17aOH-P4"]="17a-OHP4"
colnames(tvtest[[i]])[colnames(tvtest[[i]])=="HOMA IR"]="HOMA-IR"}
tv_c=tvtest[[1]]; tvf=tvtest[[2]]; tvm=tvtest[[3]];
x4[x4=="X7.oxo.DCA_S"]="X7-oxo-DCA_S"

dat = tv_c; dat = dat %>% select(-c('PatientNumber')) #this is quite nice way to delete columns, please remember...
resulta <- (rcorr(as.matrix(dat), type = c('spearman')))$r #compare pearson
intersect(colnames(resulta), rownames(resulta)) #https://stackoverflow.com/questions/45271448/r-finding-intersection-between-two-vectors

dat=tvf; dat= dat %>% select(-c('PatientNumber','Gender')) #this is quite nice way to delete columns, please remember...
resultaf <- (rcorr(as.matrix(dat), type = c('spearman')))$r #compare pearson

intersect(colnames(resultaf), rownames(resultaf)) #https://stackoverflow.com/questions/45271448/r-finding-intersection-between-two-vectors
dat=tvm; dat= dat %>% select(-c('PatientNumber','Gender')) #this is quite nice way to delete columns, please remember...
resultam <- (rcorr(as.matrix(dat), type = c('spearman')))$r #compare pearson
intersect(colnames(resultam), rownames(resultam)) #https://stackoverflow.com/questions/45271448/r-finding-intersection-between-two-vectors

#Check the columns away
at=colnames(resulta)[1:(length(x1)-1)] #clinicals
bt=colnames(resulta)[(length(at)+1):(length(at)+length(x2))] #Steroids
ct=colnames(resulta)[(length(at)+length(bt)+1):(length(at)+length(bt)+length(x3))] #BA_l
dt=colnames(resulta)[(length(at)+length(bt)+length(ct)+1):(length(at)+length(bt)+length(ct)+length(x4))] #BA_s
et=colnames(resulta)[(length(at)+length(bt)+length(ct)+length(dt)+1):(length(at)+length(bt)+length(ct)+length(dt)+length(x5))] #PFAS: change here
ft=colnames(resulta)[(length(at)+length(bt)+length(ct)+length(dt)+length(et)+1):(length(at)+length(bt)+length(ct)+length(dt)+length(et)+length(x6))] #
atl=length(at);btl=length(bt);ctl=length(ct);dtl=length(dt);etl=length(et);ftl=length(ft)

n_level=0.2; ## muuta tätä
# hist(as.numeric(Nrr)) #https://www.geeksforgeeks.org/elementwise-matrix-multiplication-in-r/

Nrr=qpNrr(resulta, verbose=FALSE);Nrr[is.na(Nrr)]=1; cond=data.frame(as.matrix(Nrr<n_level));RN=data.frame(resulta);tes_t=cond*RN;tes_t=as.matrix(tes_t);resulta=tes_t 
Nrr=qpNrr(resultaf, verbose=FALSE);Nrr[is.na(Nrr)]=1;cond=data.frame(as.matrix(Nrr<n_level));RN=data.frame(resultaf);tes_t=cond*RN;tes_t=as.matrix(tes_t);resultaf=tes_t
Nrr=qpNrr(resultam, verbose=FALSE);Nrr[is.na(Nrr)]=1;cond=data.frame(as.matrix(Nrr<n_level));RN=data.frame(resultam);tes_t=cond*RN;tes_t=as.matrix(tes_t);resultam=tes_t

# sum(resulta>0):sum(resultaf>0);sum(resultam>0)

#Plotting to chord plots at the same time...
#This function takes quite a bit of variables..
# circos.clear(); dev.off()
vars=list(resultaf,resultam); #
big='No';title='Genders Separated'; #or 'Yes' for the big plot alone
rem=x4; modi=4; colt='black';#rem=x3; modi=5; colt='green';
# rem=0; 
# rem='';modi=''; colt=''; #check this..
fig_name=paste('Chord Diagrams_v8 with',title,'_NRR_',n_level,date,'.png') # removing any of the clusters, e.g. liver bas: rem=x4; and colt='black'
two_chords(vars,n_level,fig_name,big,rem,modi,colt) #this drives the function
#for both variables
circos.clear(); dev.off()
# vars=resulta
vars=list(resulta)
big='Yes';title='All Variables' # 
rem=x4; modi=5; colt='black'
fig_name=paste('Chord Diagrams_v7a with',title,'_NRR_',n_level,date,'.png')
two_chords(vars,n_level,fig_name,big,rem,modi,colt) #this drives the function


return(hist(as.numeric(Nrr)))}

n_level=0.2; 
circos.clear(); dev.off()
nne=double_circular(tv_half_log22,n_level,date,x1,x2,x3,x4,x5,x6)

resok=resulta

resulta=resok

#This function takes elements from: https://jokergoo.github.io/circlize_book/book/advanced-layout.html#combine-circular-plots
two_chords=function(vars,n_level,fig_name, big,rem,modi,colt) {
  classes=5;
  tot=rownames(resulta)[2:dim(resulta)[1]];
  # resulta=resulta[2:dim(resulta)[1],2:dim(resulta)[2]]
  a=length(x1)-2;b=length(x2);c=length(x3);
  d=length(x4);e=length(x5);f=length(x6);#Check inside function
  range=1:(a+b+c+e+f)
  
  if (big=='Yes') {layout(matrix(1:1, 1, 1)); genders=c('Both Genders'); colors=c('black');title='All Variables'
  } else {layout(matrix(1:2, 1, 2)) ; genders= c('Female','Male');colors=c('white','black');title='Gender'}
  
  for (i in 1:length(vars)) {
    i=1
    # big='no'
    tes_t=vars[[i]];
    # tes_t=resulta[1:dim(resulta)[1],1:dim(resulta)[2]]
    
    if (big=='Yes') {colnames(tes_t)=rownames(resulta);rownames(tes_t)=rownames(resulta)} else
    {colnames(tes_t)=rownames(resultaf);rownames(tes_t)=rownames(resultaf)}
    
    a=length(x1)-2;b=length(x2);c=length(x3);d=length(x4);e=length(x5);f=length(x6);
    g1=c(rep('Clinical', a),rep('Steroids', b), rep('BA_liver', c),rep('Contaminants', e),rep('Lipids', f)) #rep('BA_serum', d)
    # removing self-correlation
    tes_t[1:a,1:a]=0
    tes_t[(a+1):(a+b),(a+1):(a+b)]=0
    tes_t[(a+b+1):(a+b+c),(a+b+1):(a+b+c)]=0
    # tes_t[(a+b+c+1):(a+b+c+d),(a+b+c+1):(a+b+c+d)]=0
    tes_t[(a+b+c+1):(a+b+c+e),(a+b+c+1):(a+b+c+e)]=0
    tes_t[(a+b+c+e+1):(a+b+c+e+f),(a+b+c+e+1):(a+b+c+e+f)]=0 #if you have more groups... make this automatic, now it is not (18.1.23)
    
    group = structure(g1, names = colnames(tes_t));#group
    grid.col = structure(c(rep('blue', a),rep('red', b), rep('green', c),  rep('orange', e), rep('#756BB1', f)), 
                         names = rownames(tes_t)); #rep('black', d),
    #brewer.pal(5, 'Purples')[4] #https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/ 
    
    tes_t=tes_t[range,range];grid.col = grid.col[range] #tes_t=resulta
    g <- graph.adjacency(tes_t, mode="upper", weighted=TRUE, diag=FALSE)
    e <- get.edgelist(g); df <- as.data.frame(cbind(e,E(g)$weight)); #
    df[,3]=as.numeric(df[, 3])
    
    rango <- function(x){((x-min(x))/(max(x)-min(x)))*2-1} #just a function for the -1 to 1 thing..
    col_fun = colorRamp2(c(min(df$V3), 0,max(df$V3)), c("blue",'white', "red"))
    df=df[!df$V1 %in% rem,];df=df[!df$V2 %in% rem,] #e.g.rem=x4
    # df$V3=rango(df$V3);
    
    for (i in 1:2) {
    df[,i]=  gsub("\\.", "-", df[,i])
    df[,i] <- gsub("X11", "11", df[,i])
    df[,i] <- gsub("X17", "17", df[,i])
    df[,i][df[,i]=="T-Epi-T"]="T/Epi-T"
    df[,i][df[,i]=="Steatosis.Grade"]="Steatosis Grade"
    df[,i][df[,i]=="Steatosis-Grade"]="Steatosis Grade"
    df[,i][df[,i]=="Fibrosis.Stage"]="Fibrosis Stage"
    df[,i][df[,i]=="Fibrosis-Stage"]="Fibrosis Stage"
    df[,i][df[,i]=="17aOH.P4"]="17a-OHP4"
    df[,i][df[,i]=="HOMA.IR"]="HOMA-IR"}
    
    classes=modi #modi=4
    # grid.col=grid.col[grid.col!=colt]# ='black'
    # namesh=unique(g1)[c(1:6)[1:6 != modi]];
    # cola=unique(grid.col)#[c(1:6)[1:6 != modi]]
    
    # grid.col=grid.col[grid.col!=colt]# ='black'
    namesh=unique(g1)    #[c(1:6)[1:6 != modi]];
    cola=unique(grid.col)#[c(1:6)[1:6 != modi]]
    
    # lgd_group = Legend(at = genders[i], type = "points", legend_gp = gpar(col = colors[i]), 
                       # title_position = "topleft", title = title)
    
    lgd_group = Legend(at = 'Both Genders', type = "points", legend_gp = gpar(col = 'black'), 
                       title_position = "topleft", title = title)
    
    lgd_points = Legend(at = namesh, type = "points", legend_gp = gpar(col = cola), title_position = "topleft", title = "Class")
    lgd_lines = Legend(at = c("Positive", "Negative"), type = "points", legend_gp = gpar(col = c('red','blue')), title_position = "topleft", title = "Correlation")#round(min(df$V3)), round(max(df$V3)) #round(min(df$V3)), round(max(df$V3)), -1,1
    lgd_edges= Legend(at = c(round(min(df$V3),1), round(max(df$V3),1)), col_fun = col_fun,  title_position = "topleft", title = "Edges")
    lgd_list_vertical = packLegend(lgd_group,lgd_points,  lgd_lines,lgd_edges) #lgd_lines,
    length(unique(colnames(resulta)));length(unique(df[,1]));length(unique(df[,2]));dim(df); #unique(df[,1]);unique(df[,2]);sum(colnames(resulta)=='CA')
    chordDiagram(df, annotationTrack = c("grid"),  grid.col=grid.col, directional = FALSE, 
                 order = rownames(tes_t), preAllocateTracks = 1, col = col_fun,transparency = 0.5)
    # unique(c(names(table(df[,1])),names(table(df[,2]))))
    #link.visible = df$V3 > 0, #big.gap = 60, small.gap = 1,
    
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim"); ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
      circos.axis(h = "top", labels.cex = 0.000001, major.tick.length = 0.2, sector.index = sector.name, track.index = 2)}, bg.border = NA) #https://stackoverflow.com/questions/31943102/rotate-labels-in-a-chorddiagram-r-circlize
    
    windowsFonts(A = windowsFont("Calibri (Body)")) 
    
    
    if (i == 1) {draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))} else
      {draw(lgd_list_vertical, x = unit(243, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))}
  # 
  if (big=='Yes') {dev.copy(jpeg,fig_name,width=9, height=12, units="in", res=1000);dev.off()} else
  {dev.copy(jpeg,fig_name,width=18, height=12, units="in", res=1000);dev.off() } 
    
    }

# draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))
# dev.copy(jpeg,fig_name,width=9.2, height=12.5, units="in", res=1000);dev.off()

  }


tv_half_log2 #...

n_level=0.2; 
circos.clear(); dev.off()
nne=double_circular(tv_covscl,n_level,date,x1,x2,x3,x4,x5,x6)

# 

# corrplot(tes_t)
# More info regarding the function:
# https://jokergoo.github.io/circlize_book/book/legends.html
# https://cran.r-project.org/web/packages/ggplotify/vignettes/ggplotify.html
# https://bioinfo4all.wordpress.com/2021/03/13/tutorial-7-how-to-do-chord-diagram-using-r/
# https://jokergoo.github.io/circlize_book/book/advanced-usage-of-chorddiagram.html
# https://jokergoo.github.io/circlize_book/book/a-complex-example-of-chord-diagram.html

#Correlation matrices:
#http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
#squares are good for individual associations, because the order is the same:
# ccovae=tv[,c("AGE")]; c1=ccovae<56; c2=ccovae>44; sick_group=c1 & c2 #tv[,c("AGE")][c1 & c2]
# as.integer(sick_group)
# tv_c=cbind(tv[,1:8],tv_half_log2[,8:dim(tv_half_log2)[2]]) #check also not logged and then the auto one
# 
# x1=colnames(tv_c[,c(1:8)]); v2=dim(NAFLD)[2]+1
# x2=colnames(tv_c[,9:v2]);v3=(dim(Bali)[2]+v2)
# x3=colnames(tv_c[,(v2+1):(v3-1)]);v4=(dim(Base)[2])+v3
# x4=colnames(tv_c[,(v3+1):(v4-1)])
# x5=colnames(tv_c[,(v4):(dim(tv_c)[2])]); 
# x5a=x5[2:8] #check this..
# x6=x5[10:length(x5)] #dividing to lipids
# x5=x5a  #making sure that PFAS are separate
# # x1;x2;x3;x4;x5; # tv_auto[1:5,1:11]
# x3 <- paste(x3, "_L", sep="") #https://stackoverflow.com/questions/6984796/how-to-paste-a-string-on-each-element-of-a-vector-of-strings-using-apply-in-r
# # x4=gsub("(-[0-9]*)*.1", "", x4) #https://stackoverflow.com/questions/18997297/remove-ending-of-string-with-gsub
# x4 <- paste(x4, "_S", sep="")# https://rdrr.io/bioc/qpgraph/man/qpNrr.html
# nm = c(x1,x2,x3,x4,x5,x6) # adding the lipids, x6
# nm=c('PatientNumber','Gender','AGE','BMI','Steatosis Grade','Fibrosis Stage','Necroinflammation','HOMA-IR',nm[9:84])
# colnames(tv_c)=nm # tv_c[1:5,1:30]; NAFLD[1:2,1:28];
# # tv_c=tv_c[,c(1:50,(length(x4)+51+1):dim(tv_c)[2])]
# colnames(tv_c)[colnames(tv_c)=='MASLD[, \"PFAS\"]']='PFAS';
#removing PFAS and lipid generals:
tv_c=tv_covscl#cbind(tv[,1:8], tv_half_log2) #check also not logged and then the auto one

# tv_c=cbind(tv[,1:8], tv_all[,9:dim(tv_all)[2]]) #check also not logged and then the auto one
# x1=colnames(tv_c[,c(1:8)]); v2=dim(NAFLD)[2]+1
# x2=colnames(tv_c[,9:v2]);v3=(dim(Bali)[2]+v2)
# x3=colnames(tv_c[,(v2+1):(v3)]);v4=(dim(Base)[2])+v3
# x4=colnames(tv_c[,(v3+1):(v4-1)])
# x5=colnames(tv_c[,(v4):(dim(tv_c)[2])]); 
# x5a=x5[1:9]
# x6=x5[10:length(x5)] #dividing to lipids
# x5=x5a  #making sure that PFAS are separate
# # x1;x2;x3;x4;x5; # tv_auto[1:5,1:11]
# x3 <- paste(x3, "_L", sep="") #https://stackoverflow.com/questions/6984796/how-to-paste-a-string-on-each-element-of-a-vector-of-strings-using-apply-in-r
# x4=gsub("(-[0-9]*)*.1", "", x4) #https://stackoverflow.com/questions/18997297/remove-ending-of-string-with-gsub
# x4 <- paste(x4, "_S", sep="")# https://rdrr.io/bioc/qpgraph/man/qpNrr.html
# nm = c(x1,x2,x3,x4,x5,x6) # adding the lipids, x6
# nm=c('Subject#','Gender','AGE','BMI','Steatosis Grade','Fibrosis Stage','Necroinflammation','HOMA-IR',nm[9:91])
# colnames(tv_c)=nm # tv_c[1:5,1:30]; NAFLD[1:2,1:28];
# colnames(tv_c)[colnames(tv_c)=='MASLD[, \"PFAS\"]']='PFAS';
# #removing PFAS and lipid generals:
# x5=x5[x5!='PFAS'] ; x6=x6[x6!='Total_TG'] # x1;x2;x3;x4;x5;
# dat = tv_c; ## dat = dat %>% select(-colnames(tv_all)[1:8]) #this is quite nice way to delete columns, please remember... tv_c[,9:dim(tv_c)[2]];
# dat=dat[,!colnames(dat) %in% c('Gender','PatientNumber')]
# # dat= dat %>% select(-x3)
# resulta <- (rcorr(as.matrix(dat), type = c('spearman')))$r #compare pearson # intersect(colnames(resulta), rownames(resulta)) #https://stackoverflow.com/questions/45271448/r-finding-intersection-between-two-vectors
# 
# dat=tvf; #
# dat=dat[,!colnames(dat) %in% c('Gender','PatientNumber')]
# # dat= dat %>% select(-c('Gender')) #this is quite nice way to delete columns, please remember...
# # dat= dat %>% select(-x3)
# resultaf <- (rcorr(as.matrix(dat), type = c('spearman')))$r #compare pearson# intersect(colnames(resultaf), rownames(resultaf)) #https://stackoverflow.com/questions/45271448/r-finding-intersection-between-two-vectors
# 
# dat=tvm;  #
# dat=dat[,!colnames(dat) %in% c('Gender','PatientNumber')] #dat= dat %>% select(-c('Gender')) #this is quite nice way to delete columns, please remember...
# # dat= dat %>% select(-x3)
# resultam <- (rcorr(as.matrix(dat), type = c('spearman')))$r #compare pearson # intersect(colnames(resultam), rownames(resultam)) #https://stackoverflow.com/questions/45271448/r-finding-intersection-between-two-vectors


#The correlations (ok):
# ccovae=tv[,c("AGE")]; c1=ccovae<56; c2=ccovae>44; sick_group=c1 & c2

tv_c=tv_covscl#tv_half_log22 #cbind(tv[,1:8], tv_half_log2) #check also not logged and then the auto one
# tv_c[,'Menopause']=as.numeric(sick_group)
tv_c=tv_c[,c(1:3,4:(dim(tv_c)[2]))]#dim(tv_c)[2],
tv_c=tv_c[,!colnames(tv_c) %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]
tv_c=tv_c[,!colnames(tv_c) %in% x4]
colnames(tv_c)[colnames(tv_c)=="17aOH-P4"]="17a-OHP4"
tvf=tv_c[tv_c[,'Gender']==min(tv_c[,'Gender']),1:dim(tv_c)[2]] #tv['Steatosis.Grade.0.To.3'==0,9:27]] #tv[tv[,'Necroinflammation']==0,9:80]; #SG0i=as.numeric(SG0i); check also: tv[tv[,'HOMA-IR']==0,9:80]
tvm=tv_c[tv_c[,'Gender']==max(tv_c[,'Gender']),1:dim(tv_c)[2]]
# tvf=tv_c[tv_c[,"SEX.1F.2M"]==min(tv_c[,"SEX.1F.2M"]),1:dim(tv_c)[2]] #tv['Steatosis.Grade.0.To.3'==0,9:27]] #tv[tv[,'Necroinflammation']==0,9:80]; #SG0i=as.numeric(SG0i); check also: tv[tv[,'HOMA-IR']==0,9:80]
# tvm=tv_c[tv_c[,"SEX.1F.2M"]==max(tv_c[,"SEX.1F.2M"]),1:dim(tv_c)[2]]
# "SEX.1F.2M"

# colnames(tv_c)
# rango <- function(x){((x-min(x))/(max(x)-min(x)))*2-1} 
rango = function(x,mi,ma) {(ma-mi)/(max(x)-min(x))*(x-min(x))+mi}
dat = tv_c; 
dat=dat[,!colnames(dat) %in% c('Gender','PatientNumber')] #SEX.1F.2M
resulta <- (rcorr(as.matrix(dat), type = c('spearman')))$r #compare pearson # intersect(colnames(resulta), rownames(resulta)) #https://stackoverflow.com/questions/45271448/r-finding-intersection-between-two-vectors
p.mat.a=rcorr(as.matrix(dat), type = c('spearman'))$P; 
p.mat.a[is.na(p.mat.a)]=1; 
p.mat.aa=matrix(p.adjust(p.mat.a,method="BH"),nrow=dim(p.mat.a)[1],ncol=dim(p.mat.a)[2]); 
rownames(p.mat.aa)=rownames(p.mat.a);colnames(p.mat.aa)=colnames(p.mat.a)
# write.csv(resulta,'MASLD_steroid_study correlations with spearman_log_tikka12424.csv')
# resulta=dat
dat=tvf; # 
dat=dat[,!colnames(dat) %in% c('Gender','PatientNumber')] #SEX.1F.2M
resultaf <- (rcorr(as.matrix(dat), type = c('spearman')))$r #compare pearson# intersect(colnames(resultaf), rownames(resultaf)) #https://stackoverflow.com/questions/45271448/r-finding-intersection-between-two-vectors
p.mat.f=rcorr(as.matrix(dat), type = c('spearman'))$P
p.mat.f[is.na(p.mat.f)]=1; 
p.mat.ff=matrix(p.adjust(p.mat.f,method="BH"),nrow=dim(p.mat.f)[1],ncol=dim(p.mat.f)[2]); 
rownames(p.mat.ff)=rownames(p.mat.f);colnames(p.mat.ff)=colnames(p.mat.f)
# write.csv(resultaf,'MASLD_steroid_study female correlations with spearman_log_tikka12424.csv')
dat=tvm;  # 
dat=dat[,!colnames(dat) %in% c('Gender','PatientNumber')] #dat= dat %>% select(-c('Gender')) #this is quite nice way to delete columns, please remember...
resultam <- (rcorr(as.matrix(dat), type = c('spearman')))$r #compare pearson # intersect(colnames(resultam), rownames(resultam)) #https://stackoverflow.com/questions/45271448/r-finding-intersection-between-two-vectors
p.mat.m=rcorr(as.matrix(dat), type = c('spearman'))$P
p.mat.m[is.na(p.mat.m)]=1; 
p.mat.mm=matrix(p.adjust(p.mat.m,method="BH"),nrow=dim(p.mat.m)[1],ncol=dim(p.mat.m)[2]); 
rownames(p.mat.mm)=rownames(p.mat.m);colnames(p.mat.mm)=colnames(p.mat.m)
# write.csv(resultam,'MASLD_steroid_study male correlations with spearman_log_tikka12424.csv')
resulta[resulta==1]=0
resultam[resultam==1]=0
resultaf[resultaf==1]=0
min(resultaf)
max(resultaf)
n_level=0.2


# dat = tv_c; 
# dat=dat[,!colnames(dat) %in% c('Gender','PatientNumber')] #SEX.1F.2M
# resulta <- (rcorr(as.matrix(dat), type = c('spearman')))$r
# n_level=0.2
Nrr=qpNrr(resulta, verbose=FALSE);Nrr[is.na(Nrr)]=1;print(hist(as.numeric(Nrr),breaks=50)); cond=data.frame(as.matrix(Nrr<n_level))
RN=data.frame(resulta);tes_t=cond*RN;tes_t=as.matrix(tes_t);resulta=tes_t;colnames(resulta)=rownames(resulta) #https://www.geeksforgeeks.org/elementwise-matrix-multiplication-in-r/
# tv_ah <- prep.autoscale(resulta, center = TRUE, scale = TRUE);
# tv_ah=rango(resulta,-0.4,0.4);
# resulta=tv_ah;
# colnames(resulta)=rownames(resulta)
# hist(resulta, breaks=50)

# Nrr=qpNrr(resultaf, verbose=FALSE);Nrr[is.na(Nrr)]=1;print(hist(as.numeric(Nrr),breaks=50));cond=data.frame(as.matrix(Nrr<n_level))
# RN=data.frame(resultaf);tes_t=cond*RN;tes_t=as.matrix(tes_t);resultaf=tes_t;colnames(resultaf)=rownames(resultaf)
# tv_ah <- prep.autoscale(resultaf, center = TRUE, scale = TRUE);
# tv_ah=rango(resultaf,min(resulta),max(resulta));
# resultaf=tv_ah;
# colnames(resultaf)=rownames(resultaf)
# hist(resultaf, breaks=50)
# Nrr=qpNrr(resultam, verbose=FALSE);Nrr[is.na(Nrr)]=1;print(hist(as.numeric(Nrr),breaks=50)); cond=data.frame(as.matrix(Nrr<n_level))#n_level=0.4;
# RN=data.frame(resultam);tes_t=cond*RN;tes_t=as.matrix(tes_t);resultam=tes_t;colnames(resultam)=rownames(resultam)
# tv_ah <- prep.autoscale(resultam, center = TRUE, scale = TRUE);
# tv_ah=rango(resultam,min(resultaf),max(resultaf));
# resultam=tv_ah;
# colnames(resultam)=rownames(resultam)
# hip1='h_with_L only_oke'
# jpeg(paste("Circle Correlation Plot of All Variables",hip1,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 16, res=300); 
# corrplot(resulta, method="square",order = "hclust");
# dev.off()
# 
# jpeg(paste("Circle Correlation Plot of Females",hip1,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 16, res=300);  
# corrplot(resultaf, method="square");
# dev.off()

# jpeg(paste("Circle Correlation Plot of Males",hip1,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 16, res=300); 
# corrplot(resultam, method="square") 
# dev.off()

# # jpeg(paste("Circle Correlation Plot of All Variables NRR checked.jpg"), width = 20000, height = 20000, quality = 100,pointsize = 16, res=1000); 
# # corrplot(tes_t, method="circle");dev.off() 
# #Triangular matrix reveal associations of a compound to many other faster
# # rev(COL2('RdBu'))
# https://www.rdocumentation.org/packages/corrplot/versions/0.92/topics/corrplot
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
order="original" #alphabet, hclust, original #https://stackoverflow.com/questions/51115495/how-to-keep-order-of-the-correlation-plot-labels-as-same-in-the-datafile
range='orig';corre='re_renorma'; method='color' #color square
jpeg(paste("Correlations with Full Plot of All_vok_nes",n_level,order,range,corre,method,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 23, res=300);
corrplot(resulta, type = "lower", order = order,method=method, tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu')),is.corr = FALSE) #,is.corr = FALSE
dev.off()
jpeg(paste("Correlations with Full Plot of Female_vok",n_level,order,range,corre,method,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 23, res=300);
corrplot(resultaf, type = "lower", order = order,method=method,tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu')),is.corr = FALSE)
dev.off()
jpeg(paste("Correlations with Full Plot of Male_voka",n_level,order,range,corre,method,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 23, res=300); 
corrplot(resultam, type = "lower", order = order, method=method,tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu')),is.corr = FALSE) #order = "alphabet",  order = "hclust",
dev.off()

# install.packages("RcmdrMisc")
# library('RcmdrMisc')
# oi=rcorr.adjust(dat, type = c("spearman"))
# oi2=as.matrix(oi$R$r)
# jpeg(paste("Correlations with Full Plot of All_v9_altb.jpg"), width = 8000, height = 8000, quality = 100,pointsize = 21, res=300);
# corrplot(oi2, type = "full", method='square', tl.col = "black", tl.srt = 45, diag = FALSE,col = COL2('RdBu'),is.corr = FALSE) #is.corr = FALSE is.corr = FALSE,col.lim = c(-1, 1)
# dev.off()
# n_level=0.4

Nrr=qpNrr(resulta, verbose=FALSE);Nrr[is.na(Nrr)]=1;print(hist(as.numeric(Nrr),breaks=50)); cond=data.frame(as.matrix(Nrr<n_level))
RN=data.frame(resulta);tes_t=cond*RN;tes_t=as.matrix(tes_t);resulta=tes_t;colnames(resulta)=rownames(resulta) #https://www.geeksforgeeks.org/elementwise-matrix-multiplication-in-r/
tv_ah <- prep.autoscale(resulta, center = TRUE, scale = TRUE); tv_ah=rango(tv_ah,-1,1); resulta=tv_ah;colnames(resulta)=rownames(resulta)
#PFAS vs. steroids; PFAS vs. lipids/bas; steroids vs. lipds/bas -kuvaajat
x5=x5[!x5=='Perfluorodecyl.ethanoic.acid']
colnames(resulta)[colnames(resulta)=="17aOH-P4"]="17a-OHP4"
x2[x2=="17aOH-P4"]="17a-OHP4"

resulta1=resultaf[x5,x2];x2
resulta2=resultaf[x5,c(x3,x6)];c(x3,x6)
resulta3=resultaf[x2,c(x3,x6)]
order="original" #alphabet, hclust, original #https://stackoverflow.com/questions/51115495/how-to-keep-order-of-the-correlation-plot-labels-as-same-in-the-datafile
range='orig';corre='no_renorm'; method='color' #color square
type='full'
hip1='h_r1_with_L_ok_fem_rnga4'
jpeg(paste("Circle Correlation Plot of PFAS vs. steroids NRR and -1to1 checked_female",hip1,".jpg"), width = 8000, height = 3000, quality = 100,pointsize = 30, res=300); 
corrplot(resulta1, type = type, order = order,method=method, tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu'))) #
dev.off()
jpeg(paste("Circle Correlation Plot of PFAS vs. lipids_bas NRR and -1to1 checked_female",hip1,".jpg"), width = 8000, height = 4000, quality = 100,pointsize = 30, res=300);  
corrplot(resulta2, type = type, order = order,method=method, tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu'))) #,is.corr = FALSE
dev.off()
jpeg(paste("Circle Correlation Plot of steroids vs. lipds_bas and -1to1 checked_female",hip1,".jpg"), width = 8000, height = 6000, quality = 100,pointsize = 30, res=300); 
corrplot(resulta3, type = type, order = order,method=method, tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu'))) #
dev.off()

resulta1=resultam[x5,x2];x2
resulta2=resultam[x5,c(x3,x6)];c(x3,x6)
resulta3=resultam[x2,c(x3,x6)]
hip1='h_r1_with_L_ok_male_rnga4'
jpeg(paste("Square Correlation Plot of PFAS vs. steroids NRR and -1to1 checked_male",hip1,".jpg"), width = 8000, height = 3000, quality = 100,pointsize = 30, res=300); 
corrplot(resulta1, type = type, order = order,method=method, tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu'))) #
dev.off()
jpeg(paste("Square Correlation Plot of PFAS vs. lipids_bas NRR and -1to1 checked_male",hip1,".jpg"), width = 8000, height = 4000, quality = 100,pointsize = 30, res=300);  
corrplot(resulta2, type = type, order = order,method=method, tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu'))) #
dev.off()
jpeg(paste("Square Correlation Plot of steroids vs. lipds_bas and -1to1 checked_male",hip1,".jpg"), width = 8000, height = 6000, quality = 100,pointsize = 30, res=300); 
corrplot(resulta3, type = type, order = order,method=method, tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu'))) #,is.corr = FALSE
dev.off()

resulta1=resulta[x5,x2];x2
resulta2=resulta[x5,c(x3,x6)];c(x3,x6)
resulta3=resulta[x2,c(x3,x6)]
hip1='h_r1_with_L_ok_all_rnga4'
tv_ah=rango(resulta1,-0.2,0.2); resulta1=tv_ah;
tv_ah=rango(resulta2,-0.2,0.2); resulta2=tv_ah;
tv_ah=rango(resulta2,-0.2,0.2); resulta2=tv_ah;
jpeg(paste("Square Correlation Plot of PFAS vs. steroids NRR and -1to1 checked",hip1,".jpg"), width = 8000, height = 3000, quality = 100,pointsize = 30, res=300); 
corrplot(resulta1, type = type, order = order,method=method, tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu')),is.corr = FALSE) #
dev.off()
jpeg(paste("Square Correlation Plot of PFAS vs. lipids_bas NRR and -1to1 checked",hip1,".jpg"), width = 8000, height = 4000, quality = 100,pointsize = 30, res=300);  
corrplot(resulta2, type = type, order = order,method=method, tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu')),is.corr = FALSE) #
dev.off()
jpeg(paste("Square Correlation Plot of steroids vs. lipds_bas and -1to1 checked",hip1,".jpg"), width = 8000, height = 6000, quality = 100,pointsize = 30, res=300); 
corrplot(resulta3, type = type, order = order,method=method, tl.col = "black", tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu')),is.corr = FALSE) #,is.corr = FALSE
dev.off()

#The ok ones:
x1=colnames(resulta)[c(1:6)]
x1=c(x1[3:6],x1[1],x1[2])#x1[2],
x5=x5[!x5=='Perfluorodecyl.ethanoic.acid']
colnames(resulta)[colnames(resulta)=="17aOH-P4"]="17a-OHP4"
colnames(p.mat.a)[colnames(p.mat.a)=="17aOH-P4"]="17a-OHP4"
x2[x2=="17aOH-P4"]="17a-OHP4"
x2=x2[order(match(x2,groups[,2]))] #https://stackoverflow.com/questions/1568511/how-do-i-sort-one-vector-based-on-values-of-another

x5=x5[!x5=='Perfluorodecyl.ethanoic.acid']

resulta1=resulta[c(x1,x2),x5];p.mat.a1=p.mat.aa[c(x1,x2),x5]
resulta2=resultaf[c(x1,x2),x5];p.mat.f1=p.mat.ff[c(x1,x2),x5]
resulta3=resultam[c(x1,x2),x5];p.mat.m1=p.mat.mm[c(x1,x2),x5]
# tv_ah=rango(resulta3,(min(resulta2)),max(resulta2)); resulta3=tv_ah;#
hip1='transposesa_kaikki scale';width = 2400;height=6000;pch.cex=1.2;
ho='PFAS vs. clinical factors and steroids'
resulta1=t(resulta1);resulta2=t(resulta2);resulta3=t(resulta3)
p.mat.a1=t(p.mat.a1);p.mat.f1=t(p.mat.f1);p.mat.m1=t(p.mat.m1)
width = 6000;height=2800;

resulta1=resulta[c(x3,x6),x5];p.mat.a1=p.mat.aa[c(x3,x6),x5]
resulta2=resultaf[c(x3,x6),x5];p.mat.f1=p.mat.ff[c(x3,x6),x5]
resulta3=resultam[c(x3,x6),x5];p.mat.m1=p.mat.mm[c(x3,x6),x5]
tv_ah=rango(resulta3,(min(resulta2)),max(resulta2)); resulta3=tv_ah;#
hip1='transpose';width = 2400;height=6000;pch.cex=1.2;ho='PFAS vs. BAs and lipids'
resulta1=t(resulta1);resulta2=t(resulta2);resulta3=t(resulta3)
p.mat.a1=t(p.mat.a1);p.mat.f1=t(p.mat.f1);p.mat.m1=t(p.mat.m1)
width = 9000;height=2800;

resulta1=resulta[c(x2),x5];p.mat.a1=p.mat.aa[c(x2),x5]
resulta2=resultaf[c(x2),x5];p.mat.f1=p.mat.ff[c(x2),x5]
resulta3=resultam[c(x2),x5];p.mat.m1=p.mat.mm[c(x2),x5]
tv_ah=rango(resulta3,(min(resulta2)),max(resulta2)); resulta3=tv_ah;#
hip1='transposes';width = 2800;height=6000;pch.cex=1.2;ho='PFAS vs. steroids_v2'
resulta1=t(resulta1);resulta2=t(resulta2);resulta3=t(resulta3)
p.mat.a1=t(p.mat.a1);p.mat.f1=t(p.mat.f1);p.mat.m1=t(p.mat.m1)

resulta1

resulta1=resulta[c(x1,x3,x5,x6),x2];  p.mat.a1=p.mat.aa[c(x1,x3,x5,x6),x2]
resulta2=resultaf[c(x1,x3,x5,x6),x2]; p.mat.f1=p.mat.ff[c(x1,x3,x5,x6),x2]
resulta3=resultam[c(x1,x3,x5,x6),x2]; p.mat.m1=p.mat.mm[c(x1,x3,x5,x6),x2]
tv_ah=rango(resulta3,(min(resulta2)),max(resulta2)); resulta3=tv_ah;#
hip1='transpose'; width = 4000;pch.cex=0.7;ho='steroids vs. all others'

resulta1=resulta[c(x1,x3,x6),x2];  p.mat.a1=p.mat.aa[c(x1,x3,x6),x2]
resulta2=resultaf[c(x1,x3,x6),x2]; p.mat.f1=p.mat.ff[c(x1,x3,x6),x2]
resulta3=resultam[c(x1,x3,x6),x2]; p.mat.m1=p.mat.mm[c(x1,x3,x6),x2]
# tv_ah=rango(resulta3,(min(resulta2)),max(resulta2)); resulta3=tv_ah;#
hip1='transpose'; width = 3700;height=6300;ho='steroids vs. all others except PFAS';ps=28 #pch=10;
min(c(resulta2)); max(c(resulta2)) #These are around -0.4 and 0.4

tv_ah=rango(resulta1,-0.5,0.5); resulta1=tv_ah;#
tv_ah=rango(resulta2,-0.5,0.5); resulta2=tv_ah;#
tv_ah=rango(resulta3,-0.5,0.5); resulta3=tv_ah;# # tv_ah=rango(resulta3,(min(resulta2)),max(resulta2)); resulta3=tv_ah;#

#https://www.rdocumentation.org/packages/corrplot/versions/0.92/topics/corrplot
#https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
#https://statisticsglobe.com/change-font-size-corrplot-r
#order can be: alphabet, hclust, original #https://stackoverflow.com/questions/51115495/how-to-keep-order-of-the-correlation-plot-labels-as-same-in-the-datafile

order="original"; range='orig';corre='no_renorm'; type='full'; method='color';ga='All';gf='Female';gm='Male' #color square
cl.offset=1.0;cl.length=5;cl.cex = 1.3;pch.cex=1.3;pch=20;cl.pos = 'r';#cl.pos = 'b' ;#pch.cex=0.95,1.3; height=6300; pos 'b' cl.pos = 'b' 
jpeg(paste("Square Correlation Plot of",ho,ga,hip1,"2.jpg"), width = width, height = height, quality = 100,pointsize = 30, res=300);# par( ps=ps)# par(cex.lab=90)
corrplot(resulta1, type = type, order = order,method=method, p.mat=p.mat.a1, tl.col = "black", #sum(COL2('RdBu')=="#FF7417")
         cl.cex = cl.cex, pch.cex=pch.cex, pch.col='black',pch=pch,#pitikö vain pch lisätä pch väriin väriin... mystistä...'#FEE12B'
sig.level = c(.001,.05, .2),cl.pos = cl.pos, insig = "label_sig", cl.offset=cl.offset,cl.length=cl.length,
tl.srt = 90, diag = TRUE,col = rev(COL2('RdBu')[25:(length(COL2('RdBu'))-25)]),is.corr = FALSE) #only in age...0.001,
dev.off()
pch.cex=1.3;
jpeg(paste("Square Correlation Plot of",ho,gf,hip1,"2.jpg"), width = width, height = height, quality = 100,pointsize = 30, res=300);  
corrplot(resulta2, type = type, order = order,method=method, p.mat=p.mat.f1,tl.col = "black",
         cl.cex = cl.cex,  pch.cex=pch.cex,pch.col='black',pch=pch,
sig.level = c(.001, .05, .2), cl.pos = cl.pos, insig = "label_sig",cl.offset=cl.offset,cl.length=cl.length,
tl.srt = 90, diag = TRUE,col = rev(COL2('RdBu')[25:(length(COL2('RdBu'))-25)]),is.corr = FALSE) #
dev.off()
# pch.cex=2.9;
jpeg(paste("Square Correlation Plot of",ho,gm,hip1,"2.jpg"), width = width, height = height, quality = 100,pointsize = 30, res=300); 
corrplot(resulta3, type = type, order = order,method=method, p.mat=p.mat.m1, tl.col = "black", cl.cex = cl.cex,pch.cex=pch.cex,
         pch.col='black',pch=pch,
sig.level = c(.001, .05, .2),cl.pos = cl.pos, insig = "label_sig",cl.offset=cl.offset,cl.length=cl.length,
tl.srt = 90, diag = TRUE,col = rev(COL2('RdBu')[25:(length(COL2('RdBu'))-25)]),is.corr = FALSE) #,is.corr = FALSE
dev.off()




# https://stats.stackexchange.com/questions/484935/how-should-i-interpret-different-p-values-when-all-of-them-are-far-from-the-sign

# aputaulukkoa v2
# resulta1=resulta[c(x1,x2),x5];p.mat.a1=p.mat.aa[c(x1,x2),x5]
rt1=matrix(as.vector(resulta1)*as.numeric(p.mat.a1<0.2),nrow=dim(resulta1)[1],ncol=dim(resulta1)[2])
rownames(rt1)=rownames(resulta1);colnames(rt1)=colnames(resulta1)
# resulta2=resultaf[c(x1,x2),x5];p.mat.f1=p.mat.ff[c(x1,x2),x5]
rt2=matrix(as.vector(resulta2)*as.numeric(p.mat.f1<0.2),nrow=dim(resulta2)[1],ncol=dim(resulta2)[2])
rownames(rt2)=rownames(resulta2);colnames(rt2)=colnames(resulta2)
# resulta3=resultam[c(x1,x2),x5];p.mat.m1=p.mat.mm[c(x1,x2),x5];tv_ah=rango(resulta3,(min(resulta2)),max(resulta2)); resulta3=tv_ah;#
rt3=matrix(as.vector(resulta3)*as.numeric(p.mat.f1<0.2),nrow=dim(resulta3)[1],ncol=dim(resulta3)[2])
rownames(rt3)=rownames(resulta3);colnames(rt3)=colnames(resulta3)
rt1=t(rt1);rt2=t(rt2);rt3=t(rt3)

rt1=rt1[,!colnames(rt1) %in% Outcome]; rt1=rt1[,c(4:7,1:3,8:(dim(rt1)[2]-3))]
rt2=rt2[,!colnames(rt2) %in% Outcome];rt2=rt2[,c(4:7,1:3,8:(dim(rt2)[2]-3))]
rt3=rt3[,!colnames(rt3) %in% Outcome];rt3=rt3[,c(4:7,1:3,8:(dim(rt3)[2]-3))]

groupse=groups[groups[,'Abbreviation']!='F',]
groupse=groupse[order(groupse[,'Group']),]
rt1=rt1[groupse[,'Abbreviation'],];rt2=rt2[groupse[,'Abbreviation'],];rt3=rt3[groupse[,'Abbreviation'],]; #rt1=data.frame(cbind(groupse[,'Group'],rt1))

totaalis=c();rtt=c()
for (i in 1:dim(rt1)[2])
     {rtt=cbind(rt1[,i],rt2[,i],rt3[,i]); 
     totaalis=append(totaalis,rtt)}

tto=matrix(totaalis,nrow=dim(rt1)[1],ncol=dim(rt1)[2]*3)
rownames(tto)=rownames(rt1);
colnames(tto)=rep(colnames(rt1),times=rep(3,dim(rt1)[2]))
tto=data.frame(cbind(groupse[,'Group'],tto))
colnames(tto)[1]='Steroid Group'# c(3:7,1,2,8:dim(rt1)[2])
tta=matrix(as.numeric(unlist(tto)),nrow=dim(tto)[1],ncol=dim(tto)[2])
rownames(tta)=rownames(tto);
colnames(tta)=colnames(tto)
# tta=tta[,2:dim(tta)[2]]

colnames(tta)=sub('.1' , ' Female', colnames(tta))
colnames(tta)=sub('.2' , ' Male', colnames(tta))
colnames(tta)=sub('\\.' , ' ', colnames(tta))
tte=tta

tte[tte > 0] <-  1 
tte[tte < 0] <- -1 



#The testosterone-to-estradiol ratio (T/E2 ratio) 
median(tv[,'T/Epi-T'])/median(tv[,'E2'])
wo=tv[,'SEX.1F.2M']==1;me=tv[,'SEX.1F.2M']==2
wom=tv[wo,];men=tv[me,]
median(wom[,'T/Epi-T'])/median(wom[,'E2']*10)
median(men[,'T/Epi-T'])/median(men[,'E2']*10)
hist(tv[,'E2'],breaks=50)

tvaux=tv

mean(tv[,'T/Epi-T'])/mean(tv[,'E2']*10)

wo=tv[,'SEX.1F.2M']==1;me=tv[,'SEX.1F.2M']==2
wom=tv[wo,];men=tv[me,]

median(wom[,'T/Epi-T'])/median(wom[,'E2']*10)
median(men[,'T/Epi-T'])/median(men[,'E2']*10)


# the_cases_div=function(tv) {

tvaux=tv
asdf=tvaux
# case = 'Necroinflammation'

lope=c()


for (i in 1:15) {
  
  cases=c('Total',colnames(tv)[5:8],'Menopause','AGE', 'BMI',x5)  
  # i=1
  case=cases[i]
  
  if (case == 'Steatosis.Grade.0.To.3') {ccovae=tv[,c("Steatosis.Grade.0.To.3")]; sick_group=ccovae>0 } else if 
  (case == 'Fibrosis.Stage.0.to.4') {ccovae=tv[,c("Fibrosis.Stage.0.to.4")]; sick_group=ccovae>0} else if 
  (case == 'Necroinflammation') {ccovae=tv[,c("Necroinflammation")]; sick_group=ccovae>0} else if 
  (case == 'HOMA-IR'){ccovae=tv[,c("HOMA-IR")]; sick_group=ccovae>1.5} else if
  (case == 'Menopause') {ccovae=tv[,c("AGE")]; c1=ccovae<56; c2=ccovae>44; sick_group=c1 & c2} else if 
  (case == 'AGE') {ccovae=tv[,c("AGE")];sick_group=ccovae>median(ccovae)} else if
  (case == 'BMI') {ccovae=tv[,c("BMI")];sick_group=ccovae>median(ccovae)} else if
  (case == 'Total') {
    ccovae=as.logical(rep(1,104));sick_group=ccovae
    tvs=tvaux[sick_group,]
    wo=tvs[,'SEX.1F.2M']==1
    me=tvs[,'SEX.1F.2M']==2
    wom=tvs[wo,] 
    men=tvs[me,]
    all=tvs
    
    a=mean(all[,'T/Epi-T'])/mean(all[,'E2']*10)
    b=mean(wom[,'T/Epi-T'])/mean(wom[,'E2']*10)
    c=mean(men[,'T/Epi-T'])/mean(men[,'E2']*10)
    tot_sick=data.frame(a,b,c);tot_ok=data.frame(0,0,0);pvals=data.frame(0,0,0)
    colnames(tot_ok)=c('All','Females','Males');
    rownames(tot_ok)='Control (Average)'
    colnames(tot_sick)=c('All','Females','Males');rownames(tot_sick)='Case (Average)'
    colnames(pvals)=c('All','Females','Males');rownames(pvals)='p value'
    tot_tot=rbind(tot_ok,tot_sick,pvals);
    rownames(tot_tot)=c(paste(case,rownames(tot_tot)[1]),paste(case,rownames(tot_tot)[2]),paste(case,rownames(tot_tot)[3]))
    colnames(tot_tot)=c('Ratio of All Subjects','Ratio of Female Subjects','Ratio of Male Subjects')
    lope=rbind(lope,tot_tot)
    next
    } else if# (case == 'PFAS') {ccovae=tv[,c("PFAS")];sick_group=ccovae>median(ccovae)} else if
  (case == 'PFHpA') {ccovae=tv[,c("PFHpA")];sick_group=ccovae>median(ccovae)} else if
  (case == 'PFHxA') {ccovae=tv[,c("PFHxA")];sick_group=ccovae>median(ccovae)} else if
  (case == 'PFHxA_Branched') {ccovae=tv[,c("PFHxA_Branched")];sick_group=ccovae>median(ccovae)} else if
  (case == 'PFHxS') {ccovae=tv[,c("PFHxS")];sick_group=ccovae>median(ccovae)} else if
  (case == 'PFNA') {ccovae=tv[,c("PFNA")];sick_group=ccovae>median(ccovae)} else if
  (case == 'PFOA') {ccovae=tv[,c("PFOA")];sick_group=ccovae>median(ccovae)} else if
  (case == 'PFOS') {ccovae=tv[,c("PFOS")];sick_group=ccovae>median(ccovae)} 
  
  
  tvs=tvaux
  wo=tvs[,'SEX.1F.2M']==1
  me=tvs[,'SEX.1F.2M']==2
  wom=tvs[wo,] 
  men=tvs[me,]
  all=tvs
  
  # tot_sick=data.frame(a,b,c)
  a=all[,'T/Epi-T']/(all[,'E2']*10) #sick
  b=wom[,'T/Epi-T']/(wom[,'E2']*10)
  c=men[,'T/Epi-T']/(men[,'E2']*10)
  
  
  
  tvs=tvaux[sick_group,]
  wo=tvs[,'SEX.1F.2M']==1
  me=tvs[,'SEX.1F.2M']==2
  wom=tvs[wo,] 
  men=tvs[me,]
  all=tvs
  
  # a=mean(all[,'T/Epi-T'])/(mean(all[,'E2'])*10)
  # b=mean(wom[,'T/Epi-T'])/(mean(wom[,'E2'])*10)
  # c=mean(men[,'T/Epi-T'])/(mean(men[,'E2'])*10)
  # tot_sick=data.frame(a,b,c)
  a=all[,'T/Epi-T']/(all[,'E2']*10) #sick
  b=wom[,'T/Epi-T']/(wom[,'E2']*10)
  c=men[,'T/Epi-T']/(men[,'E2']*10)

  tv_ok=tvaux[!sick_group,]
  wo=tv_ok[,'SEX.1F.2M']==1
  me=tv_ok[,'SEX.1F.2M']==2
  wom=tv_ok[wo,] 
  men=tv_ok[me,]
  all=tv_ok
  
  aa=all[,'T/Epi-T']/(all[,'E2']*10) #non sick
  bb=wom[,'T/Epi-T']/(wom[,'E2']*10)
  cc=men[,'T/Epi-T']/(men[,'E2']*10)
  
  # print(a);print(aa);print(b);print(bb);print(c);print(cc)
  
  jaa=t.test(a,aa)
  baa=t.test(b,bb)
  caa=t.test(c,cc)
  
  # a=mean(all[,'T/Epi-T'])/(mean(all[,'E2'])*10)
  # b=mean(wom[,'T/Epi-T'])/(mean(wom[,'E2'])*10)
  # c=mean(men[,'T/Epi-T'])/(mean(men[,'E2'])*10)
  
  tot_ok=data.frame(jaa$estimate[2],baa$estimate[2],caa$estimate[2]);
  colnames(tot_ok)=c('All','Females','Males');
  rownames(tot_ok)='Control (Average)'
  tot_sick=data.frame(jaa$estimate[1],baa$estimate[1],caa$estimate[1]);
  colnames(tot_sick)=c('All','Females','Males');rownames(tot_sick)='Case (Average)'
  pvals=data.frame(jaa$p.value,baa$p.value,caa$p.value);
  colnames(pvals)=c('All','Females','Males');rownames(pvals)='p value'
  
  tot_tot=rbind(tot_ok,tot_sick,pvals)
  rownames(tot_tot)=c(paste(case,rownames(tot_tot)[1]),paste(case,rownames(tot_tot)[2]),paste(case,rownames(tot_tot)[3]))
  
  colnames(tot_tot)=c('Ratio of All Subjects','Ratio of Female Subjects','Ratio of Male Subjects')
  
  lope=rbind(lope,tot_tot)
  
}



lope[1,]=colMeans(lope[seq(4, nrow(lope), n),])
lope[2,]=colMeans(lope[seq(5, nrow(lope), n),])

lope[3,]=c(t.test(lope[seq(4, nrow(lope), n),1],lope[seq(5, nrow(lope), n),1])$p.val,
t.test(lope[seq(4, nrow(lope), n),2],lope[seq(5, nrow(lope), n),2])$p.val,t.test(lope[seq(4, nrow(lope), n),3],lope[seq(5, nrow(lope), n),3])$p.val)

lope=t(lope)
totaalis2=c();rtt=c()
for (i in 1:dim(lope)[2])
{rtt=cbind(lope[1,i],lope[2,i],lope[3,i]); 
totaalis2=append(totaalis2,rtt)}


ttm=matrix(totaalis2,nrow=1,ncol=dim(lope)[2]*3)
rownames(ttm)='all';
colnames(ttm)=rep(colnames(lope),times=rep(3,dim(lope)[2]))

# jou=data.frame(lope)

write.csv(t(lope), 'ratiosaan3.csv',col.names = TRUE)  

ttm2=t(data.frame(ttm))
ttm2=matrix(t(data.frame(ttm)),nrow=9) #ncol=length(ttm)/3

totaalis3=c();rtt=c()
for (i in 1:dim(ttm2)[2])
{rtt=rbind(ttm2[1:3,i],ttm2[4:6,i],ttm2[7:9,i]); 
totaalis3=append(totaalis3,rtt)}
ttm3=matrix(totaalis3,nrow=3)
rownames(ttm3)=c('Case avg','Control avg', 'p val');
colnames(ttm3)=c('tot','fem','mal',colnames(tte))#rep(colnames(lope),times=rep(3,dim(lope)[2]))

write.csv(ttm3, 'ratiosaan8.csv',col.names = TRUE)  


# return(jou)}

colnames(tto)=sub('.1' , ' Female', colnames(tto))
colnames(tto)=sub('.2' , ' Male', colnames(tto))
colnames(tto)=sub('\\.' , ' ', colnames(tto))

write.csv(tto, 'aputaulukko_revised2.csv')  

# data("ToothGrowth")
# df <- ToothGrowth
# df %>% t_test(len ~ 1, mu = 0)
# df %>% t_test(len ~ supp)
# df %>% t_test (len ~ supp, paired = TRUE)



#testing the variance explained...
#this is it! https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html
SG0i;resultam
hep=tv
tv2=t(hep[,9:80])
sce=SingleCellExperiment(tv2)
logcounts(sce)=tv2
sce@colData=DataFrame(hep[,1:8]) #it is DataFrame with big Ds and Fs
vars <- getVarianceExplained(sce,variables=names(colData(sce))[1:8])
plotExplanatoryVariables(vars) #this is it! https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html

#or
# https://stats.stackexchange.com/questions/79399/calculate-variance-explained-by-each-predictor-in-multiple-regression-using-r
# https://rdrr.io/github/MRCIEU/TwoSampleMR/man/get_r_from_pn.html
# https://onlinestatbook.com/2/effect_size/variance_explained.html
#https://stackoverflow.com/questions/10441437/why-am-i-getting-x-in-my-column-names-when-reading-a-data-frame
#https://stackoverflow.com/questions/27044727/removing-characters-from-string-in-r

#Boxplots:https://r-graph-gallery.com/265-grouped-boxplot-with-ggplot2.html
# tv_half_log22[,9:dim(tv_half_log22)] == min(tv_half_log22[,9:dim(tv_half_log22)])
tv_half_log225=tv_half_log22
# tv_half_log22=tv_half_log225 #you may need this
tv_half_log22[,'11-KA4'][tv_half_log22[,'11-KA4']==min(tv_half_log22[,'11-KA4'])]=median(tv_half_log22[,'11-KA4'])
# tv_half_log22 %>%  mutate_all(~ replace(.x, which.min(.x), 0))
other='23824e';ie=tv_half_log22## 'Steatosis Grade','Fibrosis Stage','Necroinflammation','HOMA-IR'
# hist(as.numeric(unlist(tv_half_log22[,9:dim(tv_half_log22)[2]])),breaks=50)

Outcome='Steatosis Grade';Out='Steatosis'; oute='Steatosis';num=0;Group='All';boxplots(ie,Group,Outcome,Out,oute,other);Group='Female';boxplots(ie,Group,Outcome,Out,oute,other);Group='Male';boxplots(ie,Group,Outcome,Out,oute,other)
Outcome='Fibrosis Stage';Out='Fibrosis'; oute='Fibrosis Stage';num=0;Group='All';boxplots(ie,Group,Outcome,Out,oute,other);Group='Female';boxplots(ie,Group,Outcome,Out,oute,other);Group='Male';boxplots(ie,Group,Outcome,Out,oute,other) 
Outcome='Necroinflammation';Out='Necroinflammation'; oute='Level';num=0;Group='All';boxplots(ie,Group,Outcome,Out,oute,other);Group='Female';boxplots(ie,Group,Outcome,Out,oute,other);Group='Male';boxplots(ie,Group,Outcome,Out,oute,other)
Outcome='HOMA-IR';Out='HOMA-IR'; oute='Level';num=1.5 ;Group='All';boxplots(ie,Group,Outcome,Out,oute,other);Group='Female';boxplots(ie,Group,Outcome,Out,oute,other);Group='Male';boxplots(ie,Group,Outcome,Out,oute,other)

#check the imputed 'ie' dataframe from the beginning... as well as num #you may want to use logs...
Group='All';boxplots(ie,Group,Outcome,Out,oute,other);Group='Female';boxplots(ie,Group,Outcome,Out,oute,other);Group='Male';boxplots(ie,Group,Outcome,Out,oute,other)

library(extrafont)
font_import() #this is important
loadfonts(device = "win") #this is important too

# Eow check what you have:
fonts()
windowsFonts()
# install.packages('ggpval')
# library(superb)
# library(ggpval)
# tvt=ie

boxplots=function(tvt,Group,Outcome,Out,oute,other) {
if (Group=='Male') {tvt=tvt[tvt[,'Gender']==2,]} else if (Group=='Female') 
{tvt=tvt[tvt[,'Gender']==1,]} else if (Group=='All') {tvt=tvt}
Steroid=rep(colnames(tvt[,9:28]), each=dim(tvt)[1])
data2=rep('Control',dim(tvt)[1])
data2[tvt[,Outcome]>num]='Case' #'Steatosis.Grade.0.To.3' #
Treatment=data2
note=unlist(tvt[,9:28])
Concentration=as.vector(note)
data=data.frame(Steroid, Treatment ,  Concentration)
data[,'Group'] = 0
# data$Steroid #check the Xs etc out..
data$Steroid [data$Steroid  == '17aOH-P4']='17a-OHP4'
# groups$Abbreviation[groups$Abbreviation == '17a-OHP4']='17aOH-P4'
for (i in 1:21) {data[data$Steroid %in% groups$Abbreviation[i],'Group']=groups$Group[i]}
title = paste(Out,"'s Effect to Concentrations of Steroids", ' in ',Group,sep="")
if (Group=='Male') {lep=theme(legend.position = "none")} else if (Group=='Female') 
{lep=theme(legend.position = "none")} else if (Group=='All') {lep=theme_classic2()+theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2,size = 25))}
# lep=theme(legend.position = "none")
e1=paste('Case (>',num,')',sep="")
e2=paste('Control (=',num,')',sep="") #≤ 
data=data[!is.na(data$Concentration),]
# grouped boxplot: https://stackoverflow.com/questions/32539222/group-boxplot-data-while-keeping-their-individual-x-axis-labels-in-ggplot2-in-r
plot = ggplot(data, aes(x=Steroid, y=Concentration, fill=Treatment))+
  geom_boxplot(notch=F, notchwidth=0.5,outlier.shape=1,outlier.size=2, coef=1.5)+
  # coord_cartesian(ylim = c(0, 60000))+
  theme(axis.text=element_text(color="black"))+
  theme_classic2()+#https://stackoverflow.com/questions/34522732/changing-fonts-in-ggplot2
  # scale_x_discrete(guide = guide_axis(angle = 90))+ https://stackoverflow.com/questions/37488075/align-axis-label-on-the-right-with-ggplot2
  theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2))+#annotate(geom="text",family="Broadway",size=20)+#annotate(geom="text",family="Calibri",size = 14)+
  theme(panel.grid.minor=element_blank())+ #http://www.sthda.com/english/wiki/ggplot2-rotate-a-graph-reverse-and-flip-the-plot
  labs(size= "Type",x = "Steroids",y = "Concentration (Log2 pM)", title=title,size = 25)+ #log2 Autoscaled
  scale_fill_manual(values=c("orange","blue"),name=oute,labels=c(e1,e2))+ #abels=c("Case (>5)", "Control (=<5)"))
  #showSignificance( c(1,1), 13.5, 0, "**")+
  facet_grid(~Group, scales = "free_x", space = "free")+lep+
  theme(text=element_text(size=12.5,family="Calibri"), #change font size of all text
          axis.text=element_text(size=25), #change font size of axis text
          axis.title=element_text(size=25), #change font size of axis titles
          plot.title=element_text(size=25), #change font size of plot title
          legend.text=element_text(size=25), #change font size of legend text
          legend.title=element_text(size=25))+theme(axis.text=element_text(color="black"))#+annotate(geom="text",family="Broadway",size=20) # add_pval(plot, pairs = list(c(1, 2)), test='wilcox.test');plot
  # theme(text=element_text(size=20,  family="Calibri"))#change font size of legend title
  #+coord_flip()+scale_y_log10() https://datavizpyr.com/horizontal-boxplots-with-ggplot2-in-r/
jpeg(paste(title ,other,".jpg"), width = 15000, height = 6500, quality = 100,pointsize = 30, res=1000); 
print(plot);dev.off()}

#Menopause correlation with age:
plot(tv[,'AGE'],tv[,'E2'],xlab="Age", ylab="E2")
abline(lm(tv[,'E2'] ~ tv[,'AGE']), col = "red", lwd = 3) #
text(paste("Correlation:", round(cor(tv[,'AGE'], tv[,'E2']), 2)), x = 32, y = 500)

#Matching a criteria and other:
m1=rownames(tv[tv[,'AGE']<55 & tv[,'E2']>600,])
m2=rownames(tv[rev(order(tv[,''])),])[1:10]
m3=rownames(tv[tv[,'E']>100000,])
m3[order(m3)]
m1[order(m1)]# unique(c(m1,m2))
setequal(m1,m3)# which(m1 %in% m2)
Reduce(intersect, list(m1,m3)) #https://intellipaat.com/community/7100/how-to-find-common-elements-from-multiple-vectors

#Heatmaps:
#https://www.intechopen.com/chapters/52527
mat <- matrix(as.vector(unlist(tv_all[,9:28])),ncol=20)
# heatmap with the defaults parameters
pheatmap(mat)
# with "rainbow" colors
pheatmap(mat, color=rainbow(50))
# with blue to red (middle color white)
pheatmap(mat, color=colorRampPalette(c("blue", "white", "red"))(50))
# https://biocorecrg.github.io/CRG_RIntroduction/pheatmap-function-from-the-pheatmap-package.html
# add column names to mat
colnames(mat) <- paste0("Sample", 1:6)
colnames(mat) =colnames(tv_all[,9:28])

# create data frame for annotation (in the case of samples, information about the experiment, for example)
annot_cols = data.frame(
  Group = c(rep("WT", 3), rep("KO", 3)), 
  TimePoint = rep(c(0, 5, 10), each=2),
  row.names = colnames(mat))

Outcome='Steatosis.Grade.0.To.3';Out='Steatosis'; oute='Steatosis';#jap1=sa; jap2=sm; jap3=sf #
Outcome='Fibrosis.Stage.0.to.4';Out='Fibrosis';oute='Fibrosis'; jap1=fa; jap2=fm; jap3=ff#Tikka26923_v1
Outcome='Necroinflammation';Out='Necroinflammation';oute='Necroinflammation'; jap1=na; jap2=nm; jap3=nf  #
Outcome='HOMA-IR';Out='HOMA-IR';oute='HOMAIR';jap1=ha; jap2=hm; jap3=hf #Tikka26923_v1

# date='tikka31123'; 
Group='All'; 
Outcome=Out
heatmaps(tv_all,Group,Outcome,Out,date)
Group='Female'; heatmaps(tv_all,Group,Outcome,Out,date)
Group='Male'; heatmaps(tv_all,Group,Outcome,Out,date)

colnames(tf)[1:28]

heatmaps=function(tv_all,Group,Outcome,Out,date) {
tf=tv_all[,1:28]
tf[tf[,'11-KA4']>5,]=mean(tf[,'11-KA4']<5)
if (Group=='Male') {NAFLDo=tf[tf[,'SEX.1F.2M']==2,]} else if (Group=='Female') 
{NAFLDo=tf[tf[,'SEX.1F.2M']==1,]} else if (Group=='All') {NAFLDo=tf}
if (Out!='HOMA-IR') {for (i in 1:2) {if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] == 0,]} else {SG1=NAFLDo[NAFLDo[,Outcome] > 0,]}}} else
    {for (i in 1:2) {if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] < 5,]} else {SG1=NAFLDo[NAFLDo[,Outcome] > 5,]}}}
sgt=rbind(SG0,SG1)
mat <- matrix(as.vector(unlist(sgt[,c(9:28)])),ncol=20)
colnames(mat) =colnames(tf[,9:28])#
rownames(mat)=rownames(sgt)#
# define the annotation
annotation_row = data.frame(Clinical = 
  factor(rep(c("Control", "Case"),c(dim(SG0)[1], dim(SG1)[1]))))
colnames(annotation_row)[colnames(annotation_row) == "Clinical"] =Out
rownames(annotation_row) = rownames(sgt)#paste("Gene", 1:105, sep = "")
#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
# https://jokergoo.github.io/ComplexHeatmap-reference/book/integrate-with-other-packages.html
color_palette=brewer.pal(n = 11, name = 'RdBu')
brks_heatmap <- function(mat, color_palette){
  rng <- range(mat, na.rm = TRUE)
  lpal <- length(color_palette)
  c(seq(rng[1], 0, length.out=ceiling(lpal/2) + 1),
    seq(rng[2]/dim(mat)[1], rng[2], length.out=floor(lpal/2)))}
# https://shanguangyu.com/articles/customized-color-key-for-pheatmap/
hope=pheatmap(mat, annotation_row = annotation_row, #annotation_col = annot_cols, 
         cluster_rows=FALSE, cluster_cols=TRUE,column_names_side = c("top"), angle_col = c("90"),
         breaks = brks_heatmap(mat, color_palette),color=color_palette) #toimi! :)
other=Out;title=Group
jpeg(paste(title ,other,date,".jpg"), width = 6500 , height =15000 , 
     quality = 100,pointsize = 16, res=1000); print(hope);dev.off()}


Outcome='Steatosis Grade';Out='Steatosis'; oute='Steatosis';jap1=sa; jap2=sm; jap3=sf #
Outcome='Fibrosis.Stage.0.to.4';Out='Fibrosis';oute='Fibrosis'; jap1=fa; jap2=fm; jap3=ff#Tikka26923_v1
Outcome='Necroinflammation';Out='Necroinflammation';oute='Necroinflammation'; jap1=na; jap2=nm; jap3=nf  #
Outcome='HOMA-IR';Out='HOMA-IR';oute='HOMAIR';jap1=ha; jap2=hm; jap3=hf #Tikka26923_v1

heatmaps=function(tv_all,Group,Outcome,Out,date) {
  tf=tv_all[,1:28]
  # tf[tf[,'11-KA4']>5,]=mean(tf[,'11-KA4']<5)
  if (Group=='Male') {NAFLDo=tf[tf[,'Gender']==2,]} else if (Group=='Female') 
  {NAFLDo=tf[tf[,'Gender']==1,]} else if (Group=='All') {NAFLDo=tf}
  if (Out!='HOMA-IR') {for (i in 1:2) {if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] == 0,]} else {SG1=NAFLDo[NAFLDo[,Outcome] > 0,]}}} else
  {for (i in 1:2) {if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] < 5,]} else {SG1=NAFLDo[NAFLDo[,Outcome] > 5,]}}}
  sgt=rbind(SG0,SG1)
  mat <- matrix(as.vector(unlist(sgt[,c(9:28)])),ncol=20) # mat=t(mat)
  colnames(mat) =colnames(tf[,9:28])#
  rownames(mat)=rownames(sgt)#
  mat=t(mat)
  # define the annotation
  annotation_row = data.frame(Clinical = factor(rep(c("Control", "Case"),c(dim(SG0)[1], dim(SG1)[1]))))
  colnames(annotation_row)[colnames(annotation_row) == "Clinical"] =Out
  rownames(annotation_row) = rownames(sgt) #paste("Gene", 1:105, sep = "")
  #https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/integrate-with-other-packages.html
  color_palette=brewer.pal(n = 11, name = 'RdBu')
  # https://shanguangyu.com/articles/customized-color-key-for-pheatmap/
  hope=pheatmap(mat, annotation_row = annotation_row, #annotation_col = annot_cols, 
               cluster_rows=FALSE, cluster_cols=FALSE,column_names_side = c("top"), angle_col = c("90"),
                breaks = brks_heatmap(mat, color_palette),color=color_palette) #toimi! :)
  
  
  # pheatmap(mat, annotation_row = annotation_row,  breaks,
  #           cluster_rows=FALSE, cluster_cols=FALSE,column_names_side = c("top"), angle_col = c("90"),color=color_palette)
  # 
  # 
  other=Out;title=Group
  jpeg(paste(title ,other,date,".jpg"), width = 6500 , height =15000 , 
       quality = 100,pointsize = 16, res=1000); print(hope);dev.off()}

d3heatmap(mat, dendrogram = 'none', key = TRUE, col = 'RdYlGn',
          scale = 'column', key.title = "Legend", print.values = T,
          notecol = 'white') %>% 
  hmAxis("x", title = "test", location = 'top') %>% 
  hmAxis("y", title = "test", location = 'left') %>% 
  hmCells(font.size = 8, color = 'blue') %>% 
  hmLegend(show = T, title = "Title", location = "tl")
# https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
# https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
heatmaply::heatmaply(mat, Rowv = FALSE,
                     Colv = FALSE,
                     xlab = "Steroids",
                     ylab = "Patient Id", 
                     scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                       low = "blue", 
                       high = "red", 
                       midpoint = (max(mat)-min(mat))/2+min(mat),
                       limits = c(min(mat), max(mat))
                     )) #%>% layout(xaxis = list(side = "top"))
#http://www.sthda.com/english/wiki/ggplot2-rotate-a-graph-reverse-and-flip-the-plot
#https://stackoverflow.com/questions/10547487/remove-facet-wrap-labels-completely

# https://r-graph-gallery.com/79-levelplot-with-ggplot2.html

#The modelling:
# SG0=tv[tv[,'Steatosis.Grade.0.To.3']==0,c(5,9:28)]
Group='male'
name="Forest plot of All Steroid Ratios in HOMA-IR"
# pre_errors=function(NAFLD,Outcome,Group,name,ordera) {
if (Group=='male') {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==2,]} else if (Group=='female') 
{NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==1,]} else if (Group=='All') {NAFLDo=NAFLD}

SG0=NAFLDo[,c(8:28)]#tv[,'Steatosis.Grade.0.To.3']==0# colnames(SG0[,2:21]) <- gsub(" ", ".", colnames(SG0[,2:21])) #https://stackoverflow.com/questions/10688137/how-to-fix-spaces-in-column-names-of-a-data-frame-remove-spaces-inject-dots
oknames=colnames(SG0)
SG0=data.frame(SG0)
colnames(SG0)
colnames(SG0[,2:21]) <- gsub("-", ".", colnames(SG0[,2:21]))
colnames(SG0[,2:21]) <- gsub("/", ".", colnames(SG0[,2:21]))
# SG0=destroyX(data.frame(SG0))

xnam <- colnames(SG0[,2:21])#paste("x", 1:25, sep="")
fmla <- as.formula(paste("HOMA.IR ~ ", paste(xnam, collapse= "+")))
#https://stats.stackexchange.com/questions/190763/how-to-decide-which-glm-family-to-use
poissone=glm( fmla, data=SG0,family = gaussian) 
anova(poissone);
summary(poissone);
poissone$coefficients;
hist(poissone$coefficients)
valuees=summary(poissone)
ms=poissone$coefficients[2:21]
error=valuees$coefficients[2:21,2]
pval=anova(poissone)[1:20,5]
error_lower=ms-error# error_lower[error_lower <= 0] = 0
error_upper=ms+error
hip=colnames(tv[,9:28])
sample_data <- data.frame(study=colnames(tv[,9:28]),index=colnames(tv[,9:28]),result=ms,error_lower=error_lower,error_upper=error_upper)#,pval=pval)
sample_data %>%
  mutate(study = fct_reorder(study, result)) %>%
  ggplot(aes(y=study, x=ms,xmin=error_lower,xmax=error_upper)) +
  geom_point(color= "black", pch= 1, size=3) +
  # coord_flip() +
  geom_errorbarh(height=.5, color= "black", lwd=0.5) +
  labs(title='Steatosis Grade Effect Forest Plot', x='LM Estimates (SD)', y = 'Steroids')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = 
          element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))#+

metad=data.frame(colnames(SG0)[2:21],'cate') #https://datatofish.com/create-dataframe-in-r/
colnames(metad) =c('name','group'); #if you know: ga=groups[,'Abbreviation']; md=metad[,'name']; unique(c(ga,md)); ga[!(ga %in% md)]; ga[ga=="17a-OHP4"]="17aOH-P4"
groups[,'Abbreviation'][groups[,'Abbreviation']=="17a-OHP4"]="17aOH-P4"
groups=groups[groups[,'Abbreviation']!='F',];rownames(groups)=1:20
groups=groups[order(groups[,'Abbreviation']),]; metad=metad[order(metad[,'name']),]
metad[,'group']=groups[,'Group']# groups[,'Abbreviation'] == metad[,'name']
cate=metad[,'group']
hup=names(table(cate)[rev(order(table(cate)))])#[c(1:15)]) #2:3, 1, 4, 5:8, 9-19
metad2=metad[metad[,'group'] %in% hup,]
df_grouping=tibble(metad2) #should be ok
sample_data2=sample_data[rev(order(sample_data[,'result'])),]
sample_data2=tibble(sample_data2)
colnames(sample_data2) =c('name','name2','result','nerror','error','ok')# Join the association data frame df_compare_traits with group data # use right_join, with df_grouping on the right, to preserve the order of biomarkers it specifies.
df_compare_traits_groups <- sample_data2 %>% dplyr::right_join(., df_grouping, by = "name") %>% dplyr::mutate(group = factor(.data$group, levels = unique(.data$group)))
df_compare_traits_groups=data.frame(df_compare_traits_groups)
df_compare_traits_groups=df_compare_traits_groups[rev(order(df_compare_traits_groups[,'result'])),]
df_compare_traits_groups=tibble(df_compare_traits_groups)

gn=df_compare_traits_groups[,c('group','name')]
# gn[order(data.frame(gn[,'name'])[,1]),]
gn[order(data.frame(gn[,'name'])[,1]),]
sample_data=sample_data[order(sample_data[,'study']),]
sample_data=cbind(sample_data,gn[order(data.frame(gn[,'name'])[,1]),])
# sample_data[,'Group']='Group'
# groups$Abbreviation[groups$Abbreviation == '17a-OHP4']='17aOH-P4'
# for (i in 1:21) {sample_data[sample_data$study %in% groups$Abbreviation[i],'Group']=groups$Group[i]}
# groups=groups[groups[,'Abbreviation']!='F',]
# # for (i in 1:20) {groups[groups$Abbreviation %in% sample_data$study[i],'Group']=sample_data$Group[i]}
sample_data=sample_data[order(sample_data[,'study']),]
sample_data=cbind(sample_data,groups[order(data.frame(groups[,'Abbreviation'])[,1]),])

cm=sign(sample_data[,'error_lower']) < 1 & sign(sample_data[,'error_upper']) < 1
cp=sign(sample_data[,'error_lower']) > 0 & sign(sample_data[,'error_upper']) > 0
ct= cm | cp
sample_data[,'Significance']=ct
sample_data[,'Significance'][sample_data[,'Significance']==TRUE]='Yes'
sample_data[,'Significance'][sample_data[,'Significance']==FALSE]='No'
sample_data[,'Color']=ct#!(sample_data[,'error_lower'] >= 0 | sample_data[,'error_upper'] >= 0)
sample_data[,'Color'][sample_data[,'Color']==TRUE]='blue'
sample_data[,'Color'][sample_data[,'Color']==FALSE]='grey'
sample_data[,'pval']=ct#!(sample_data[,'error_lower'] >= 0 | sample_data[,'error_upper'] >= 0)
sample_data[,'pval'][sample_data[,'pval']==TRUE]=0.01
sample_data[,'pval'][sample_data[,'pval']==FALSE]=0.06
sample_data2=sample_data
# sample_data2=sample_data2[,c(7,6,3:5,9:11)]
sample_data2[,unique(colnames(sample_data2))]
plot=forestplot(df = sample_data, #drive tba_example_v3_oh_tikka17823.R if not working via 'x' error
                estimate = result,
                se= abs(error_lower-error_upper)/4,
                pvalue = pval, #this makes the significant value..
                psignif = 0.1,
                xlab = "LM  Estimates (SD)",
                ylab='Steroid Groups',
                title='Forest Plot of Linear Model Estimates for Male HOMAIRs',
                colour = Significance) +
  ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "fixed",
    strip.position='left') 
if (Group=='All') {ordera=sample_data$study[order(sample_data$result)] #this should be ok
plot[["data"]][["study"]]=factor(plot[["data"]][["study"]], levels = sample_data$study[order(sample_data$result)] )} else if 
(Group!='All') {plot[["data"]][["study"]]=factor(plot[["data"]][["study"]], levels = sample_data$study[order(sample_data$result)] )}
plot$layers[[1]]$aes_params$odd <- "#00000000" 
#https://stackoverflow.com/questions/71745719/how-to-control-stripe-transparency-using-ggforestplot-geom-stripes
jop=plot #+theme(axis.text.y=element_blank());
jop2=jop+geom_point(aes(colour = factor(Significance)),colour = sample_data[,'Color']) +
  scale_color_manual(values=c('#999999','blue'));jop2 

# plot$layers[[1]]$aes_params$odd <- "#00000000" #https://stackoverflow.com/questions/71745719/how-to-control-stripe-transparency-using-ggforestplot-geom-stripes
#https://rdrr.io/rforge/CALIBERdatamanage/man/multiforest.html

#https://www.nature.com/articles/s41467-022-30875-7#citeas #with https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html
#other # https://academic.oup.com/bioinformatics/article/26/18/2336/208507?login=false#393579283 http://locuszoom.sph.umich.edu//
#https://www.cell.com/ajhg/fulltext/S0002-9297(07)61352-4#secd28548230e6399
#https://www.mdpi.com/2218-273X/13/6/917

#The tba:
includedx=read.csv("Lipid Data for Pauli.csv", header = FALSE,sep=","); #note the dot and not ;
included=read.csv("Lipid_Data_with_Meta_Data_For_Pauli (1).csv", header = FALSE,sep=";");
included[8:12,1:24]
inc2=includedx[8:2199,]
lipids=inc2[2,12:369]
lipids=unname(lipids)[1:358]
row.names(lipids) <- NULL
Start=inc2[1,1:11]; Start=unname(Start)
row.names(Start) <- NULL
colnames(inc2)=c(Start,lipids)
inc22=inc2[5:2192,1:368]
# rownames(inc22)=inc22[,'GUPi'] #duplicates occur..
jip=inc22[,7]=='Mild'
inc33=inc22[inc22[,7]=='Mild',]
inc33[is.na(inc33)]='NaN'
inc34a <- inc33[rowSums(inc33=='NaN') < 2, ] # hist(rowSums(inc33=='NaN'))
inc34a[40:50,c(8,12:13)]#Do I need the fragment values?

ic2=included[8:2199,]
lipidsi=ic2[2,24:380]
lipidsi=unname(lipidsi)[1:357]
row.names(lipidsi) <- NULL
start=ic2[1,1:23]; Start=unname(start)
row.names(start) <- NULL
colnames(ic2)=c(start,lipidsi)
#https://stackoverflow.com/questions/10688137/how-to-fix-spaces-in-column-names-of-a-data-frame-remove-spaces-inject-dots
names(ic2) <- gsub(" ", ".", names(ic2))
names(ic2) <- gsub("_", ".", names(ic2))
names(ic2) <- gsub("-", ".", names(ic2))
ic2[20:30,1:25]
ic22=ic2[5:2192,1:368]
ic33=ic22[jip,]
ic33[70:80,c(1:24)]#Do I need the fragment values?
ic33[is.na(ic33)]='NaN'
ic34 <- ic33[rowSums(ic33=='NaN') < 5, ] #Interesting is same as with the other
# hist(rowSums(ic33=='NaN'))
ic34[,24:368] = lapply(ic34[,24:368],as.numeric)
sum(ic34[,24:368]=='NaN');sum(ic34[,24:368]==0)
ic34[70:80,c(1:24)]#Do I need the fragment values?
ic34[ic34[,'Duration.of.Amnesia']=='NaN',]=0
# sum(ic34[,'WorstGCS']=='NaN')
ic34=ic34[!rownames(ic34) %in% rownames(ic34[ic34[,'WorstGCS']=='NaN',]),]
ic34=ic34[!rownames(ic34) %in% rownames(ic34[ic34[,'Education']=='NaN',]),]
ic34=ic34[!rownames(ic34) %in% rownames(ic34[ic34[,'Gose.7.8.Good.outcome']=='NaN',]),]

cate<- (str_extract(colnames(ic34[,24:368]), "[aA-zZ]+")) #https://stackoverflow.com/questions/29825537/group-categories-in-r-according-to-first-letters-of-a-string
metad=data.frame(colnames(ic34[,24:368]),cate,'trait') #https://datatofish.com/create-dataframe-in-r/

ic34[40:50,c(1:24)]
SG0=ic34[ic34[,'Gose.7.8.Good.outcome']=='Poor Outcome',]#Steatosis.Grade.0.To.3
SG0[40:50,c(24:28)]
sum(SG0[,24:368]=='NaN');sum(SG0[,24:368]==0)
means=c();for (i in 24:368) {means=append(means,mean(SG0[,i]))}
sds=c();for (i in 24:368) {sds=append(sds,sd(SG0[,i]))}
error_lower=means-sds
error_upper=means+sds
sample_data <- data.frame(study=colnames(SG0[,24:368]),index=colnames(SG0[,24:368]),
                          result=means,error_lower=error_lower,error_upper=error_upper)
sample_data %>%
  mutate(study = fct_reorder(study, result)) %>%
  ggplot(aes(y=study, x=result,xmin=error_lower,xmax=error_upper)) +
  geom_point(color= "black", pch= 1, size=3) +
  # coord_flip() +
  geom_errorbarh(height=.5, color= "black", lwd=0.5) +
  labs(title='Forest Plot', x='Log Concentration Averages (SD)', y = 'Lipids')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))#+

colnames(metad) =c('name','group','trait')
hup=names(table(cate)[rev(order(table(cate)))])[c(1:19)] #2:3, 1, 4, 5:8, 9-19
metad2=metad[metad[,'group'] %in% hup,]
df_grouping=tibble(metad2) #should be ok
sample_data2=sample_data[rev(order(sample_data[,'result'])),]
# sample_data2=sample_data2[1:300,]
sample_data2=tibble(sample_data2)
colnames(sample_data2) =c('name','name2','result','nerror','error')
# table(cate)[rev(order(table(cate)))]
# Join the association data frame df_compare_traits with group data # use right_join, with df_grouping on the right, to preserve the order of biomarkers it specifies.
df_compare_traits_groups <-
  sample_data2 %>%
  dplyr::right_join(., df_grouping, by = "name") %>% dplyr::mutate(group = factor(.data$group, levels = unique(.data$group)))
df_compare_traits_groups=data.frame(df_compare_traits_groups)
df_compare_traits_groups=df_compare_traits_groups[rev(order(df_compare_traits_groups[,'result'])),]
df_compare_traits_groups=tibble(df_compare_traits_groups)
print(df_compare_traits_groups, n=80)
# df_compare_traits_groups=df_compare_traits_groups[1:135,] #135 is too many..
df_compare_traits_groups=data.frame(df_compare_traits_groups)
df_compare_traits_groups=df_compare_traits_groups %>% group_by(group) %>% slice_max(order_by = result, n = 10)
df_compare_traits_groups=tibble(df_compare_traits_groups)
# Draw a forestplot of cross-sectional, linear associations.
plot=forestplot(
  df = df_compare_traits_groups,
  estimate = result,
  se=error,
  pvalue = NULL,
  psignif = 0.1,
  xlab = "Log Concentration (SD)",
  ylab='Lipids',
  colour = group) +
  ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "fixed")
jpeg(paste("Lipid groups's log concentrations_all_poor outcome_v" ,".jpg"), 
     width = 12000, height = 44000, quality = 100,pointsize = 16, res=1200); plot;dev.off() #height = 44000 in limit..

#https://nightingalehealth.github.io/ggforestplot/articles/ggforestplot.html
SG0m <- SG0[,c(5,6,13,18,21,24:368)]
SG0m= SG0m[rowSums(SG0m=='NaN') ==0, ] 

#https://stats.stackexchange.com/questions/190763/how-to-decide-which-glm-family-to-use


for (i in 1:3) {SG0m[,i]=as.numeric(SG0m[,i])}
xnam <- colnames(SG0m[,c(1,3:350)])#paste("x", 1:25, sep="")
colnames(SG0m)[1:5]
fmla <- as.formula(paste("WorstGCS ~ ", paste(xnam, collapse= "+")));fmla
#https://stats.stackexchange.com/questions/190763/how-to-decide-which-glm-family-to-use

poissone=lm( fmla, data=SG0m)#, family = binomial) 
ap=anova(poissone);
sp=summary(poissone);
ap[rev(order(ap$'F value')),]
table=sp[[4]][rev(order(sp[[4]][,1])),]
error=table[,'Std. Error']#valuees$coefficients[2:21,2]
ms=table[,'Estimate'] #poissone$coefficients[2:21]
error_lower=ms-error
error_upper=ms+error
# hip=colnames(NAFLD[,8:27])
sample_data <- data.frame(study=names(ms),index=names(ms),
                          result=ms,
                          error_lower=error_lower,
                          error_upper=error_upper)
plote=sample_data %>%
  mutate(study = fct_reorder(study, result)) %>%
  ggplot(aes(y=study, x=ms,xmin=error_lower,xmax=error_upper)) +
  geom_point(color= "black", pch= 1, size=3) +
  geom_errorbarh(height=.5, color= "black", lwd=0.5) +
  labs(title='WorstGCS Forest Plot', x='LM Estimates (SD)', y = 'Lipids')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))#+
# scale_fill_discrete(breaks = sample_data$study) # scale_fill_discrete(breaks = sample_data$study) 
jpeg(paste("WorstGCS LM Estimates",".jpg"), 
     width = 12000, height = 40000, quality = 100,pointsize = 16, res=1200); plote;dev.off()

#Now a simple forest plot as in the example:
library(ggplot2)
includede=read.csv("Meta_Data_GOSe2.csv", sep=";", header = TRUE); 
hist(rowSums(includede=='NaN'))
includedee <- includede[rowSums(includede=='NaN') < 2, ] #note the sign..
includedee[1:3,1:11]
rownames(includedee)=includedee[,1]
includedee=includedee[,3:11]

#Sankey plots:
#https://r-charts.com/flow/sankey-diagram-ggplot2/
# df <- mtcars %>% make_long(cyl, vs, am, gear, carb)
# ggplot(df, aes(x = x,next_x = next_x, 
hop=includedee[ic34[,'GUPi'],]# hop=includedee
hop[1:5,1:5]
hop2=hop[,5:7]
hop2=hop2[rowSums(hop2=='NaN') == 0, ] 
hop2=hop2[!is.na(hop2[,1]), ] 

colnames(hop2)=c('3','6','12')
hop2[hop2==2]='2/3'
# hop2 = na.omit(hop2)
df2 <- hop2 %>% make_long(colnames(hop2))
top_group <-c(1,'2/3',4,5,6,7,8)#c(1:8)# df %>% top_n(3, sample.size) %>% pull(group)
#https://stackoverflow.com/questions/70197763/only-displaying-some-legend-contents
ggplot(df2, aes(x = x, next_x = next_x, 
                node = node, next_node = next_node,
                fill = factor(node),
                label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white")+
  scale_fill_viridis_d(option = "A", alpha = 0.95) +
  theme_sankey(base_size = 16)+
  labs(y= "GOSe", x = "Time (months)") +
  guides(fill = guide_legend(title = "GOSe"))+
  scale_fill_discrete(breaks = top_group)
#https://rcompanion.org/rcompanion/e_07.html
#https://cran.r-project.org/web/packages/bigsnpr/index.html

# Sitten omalla:
# Data.omit = na.omit(Data)
inc35=inc34[,12:368]
hop=includedee[inc34[,'GUPi'],]
# hop=includedee
hop[1:5,1:5]
hop2=hop[,5:7]
hop2=hop2[rowSums(hop2=='NaN') == 0, ] 
colnames(hop2)=c('3','6','12')
hop2[hop2==2]='2/3'
inc35=cbind(status=hop2[,1], inc35)
inc35=inc35[!is.na(inc35$status),]
inc35[inc35$status < 5,1]=0
inc35[inc35$status > 4,1]=1

hop=includedee[ic34[,'GUPi'],]# hop=includedee
hop[1:5,1:5]
hop2=hop[,5:7]
hop2=hop2[rowSums(hop2=='NaN') == 0, ] 
hop2=hop2[!is.na(hop2[,1]), ] 

colnames(hop2)=c('3','6','12')
hop2[hop2==2]='2/3'
# hop2 = na.omit(hop2)
df2 <- hop2 %>% make_long(colnames(hop2))
top_group <-c(1,'2/3',4,5,6,7,8)#c(1:8)# df %>% top_n(3, sample.size) %>% pull(group)
#https://stackoverflow.com/questions/70197763/only-displaying-some-legend-contents
ggplot(df2, aes(x = x, next_x = next_x, 
                node = node, next_node = next_node,
                fill = factor(node),
                label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white")+
  scale_fill_viridis_d(option = "A", alpha = 0.95) +
  theme_sankey(base_size = 16)+
  labs(y= "GOSe", x = "Time (months)") +
  guides(fill = guide_legend(title = "GOSe"))+
  scale_fill_discrete(breaks = top_group)



NAFLD=cbind(tv[,1:28])#,tv_half_log2[,1:20]); #ie# tv_all=cbind(tv[,1:8],tv_auto); let's take tv 
colnames(NAFLD) #autoscaled; # NAFLD=cbind(tv[,1:8],tv_half_log2[,1:20]) #for raw # NAFLD=cbind(tv[,1:28])
NAFLD[NAFLD[,c(5)]>0,5]=1;NAFLD[NAFLD[,c(6)]>0,6]=1;NAFLD[NAFLD[,c(7)]>0,7]=1;
NAFLD[NAFLD[,c(8)] <= 5,8]=0;NAFLD[NAFLD[,c(8)]>5,8]=1; #note the order of calculation...
colnames(NAFLD) <- gsub("-", ".", colnames(NAFLD))
colnames(NAFLD) <- gsub("/", ".", colnames(NAFLD))
colnames(NAFLD) <- gsub("11", "X11", colnames(NAFLD))
colnames(NAFLD) <- gsub("17", "X17", colnames(NAFLD))
colnames(NAFLD) <- gsub("#", ".", colnames(NAFLD))

# Necroinflammation  HOMA-IR Steatosis.Grade.0.To.3 Fibrosis.Stage.0.to.4
# date='tikka31123'; 
loge=0 
Outcome='Steatosis.Grade.0.To.3';Out='Steatosis'; oute='Steatosis';jap1=sa; jap2=sm; jap3=sf #
Outcome='Steatosis Grade';
Outcome='Fibrosis.Stage.0.to.4';Out='Fibrosis';oute='Fibrosis'; jap1=fa; jap2=fm; jap3=ff#Tikka26923_v1
Outcome='Necroinflammation';Out='Necroinflammation';oute='Necroinflammation'; jap1=na; jap2=nm; jap3=nf  #
Outcome='HOMA.IR';Out='HOMA-IR';oute='HOMAIR';jap1=ha; jap2=hm; jap3=hf #Tikka26923_v1

#this is with first(!!), use it
#is normal (non-log scale) image (and autoscaled), loge=1 is  log10 scale image (and raw data)
first=TRUE; ordera=c();
Group='All';name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); 
hel=pre_errors(NAFLD,Outcome,Group,name,ordera,loge,jap1,oute,first) 

#Afterwards:
first=FALSE
Group='All';name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date);
Group='Female'; name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); 
Group='Male'; name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); 

#This works with the normal raw data
# https://www.statology.org/cohens-d-in-r/
cohd=function(NAFLD, tv,Group,Outcome) {
if (Group=='Male') {NAFLDo=NAFLD[NAFLD[,'Gender']==max(NAFLD[,'Gender']),];tva=tv[tv[,'SEX.1F.2M']==max(tv[,'SEX.1F.2M']),]} else if (Group=='Female') 
{NAFLDo=NAFLD[NAFLD[,'Gender']==min(NAFLD[,'Gender']),];tva=tv[tv[,'SEX.1F.2M']==min(tv[,'SEX.1F.2M']),]} else if (Group=='All') {NAFLDo=NAFLD;tva=tv}
if (Outcome != 'HOMA-IR') {SG0=NAFLDo[NAFLDo[,Outcome] == min(NAFLDo[,Outcome]),]; SG1=NAFLDo[NAFLDo[,Outcome] > min(NAFLDo[,Outcome]),]} else if (Outcome == 'HOMA-IR') 
  {SG0=NAFLDo[tva[,'HOMA-IR'] <= 1.5,]; SG1=NAFLDo[tva[,'HOMA-IR'] > 1.5,]}
 
  cd=c();for (i in 1:20) {cd=append(cd,cohen.d(SG1[,i+8],SG0[,i+8])$estimate)} 
 
 return(cd)}

colnames(x) <- gsub(".1", "ho", colnames(x))

d=c()
NAFLD=tv_all
Outcome='Steatosis Grade';Group='All'; a=cohd(NAFLD, tv,Group,Outcome);Group='Female'; b=cohd(NAFLD, tv,Group,Outcome);Group='Male'; c=cohd(NAFLD, tv,Group,Outcome);d=cbind(d,cbind(a,b,c))
Outcome='Fibrosis Stage';Group='All'; a=cohd(NAFLD, tv,Group,Outcome);Group='Female'; b=cohd(NAFLD, tv,Group,Outcome);Group='Male'; c=cohd(NAFLD, tv,Group,Outcome);d=cbind(d,cbind(a,b,c))
Outcome='Necroinflammation';Group='All'; a=cohd(NAFLD, tv,Group,Outcome);Group='Female'; b=cohd(NAFLD, tv,Group,Outcome);Group='Male'; c=cohd(NAFLD, tv,Group,Outcome);d=cbind(d,cbind(a,b,c))
Outcome='HOMA-IR';Group='All'; a=cohd(NAFLD, tv,Group,Outcome);Group='Female'; b=cohd(NAFLD, tv,Group,Outcome);Group='Male'; c=cohd(NAFLD, tv,Group,Outcome);d=cbind(d,cbind(a,b,c))
rownames(d)=colnames(tv_all[,9:28])
colnames(d)=rep(c('All','Female','Male'),4)
write.csv(d,'cohens_da_tikka_v11924a.csv')
n=d



# heattiin=function(NAFLD,Group,Outcome){
# if (Group=='Male') {NAFLDo=NAFLD[NAFLD[,'Gender']==max(NAFLD[,'Gender']),];tva=tv[tv[,'SEX.1F.2M']==max(tv[,'SEX.1F.2M']),]} else if (Group=='Female') 
#   {NAFLDo=NAFLD[NAFLD[,'Gender']==min(NAFLD[,'Gender']),];tva=tv[tv[,'SEX.1F.2M']==min(tv[,'SEX.1F.2M']),]} else if (Group=='All') {NAFLDo=NAFLD;tva=tv}
#   
# sample_data=c()
# for (i in 1:2) {
#   if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] == 0,]} else if (i==2) {SG0=NAFLDo[NAFLDo[,Outcome] > 0,]}
#   means=c();for (j in 9:28) {means=append(means,mean(SG0[,j], na.rm=TRUE))}
#   means[is.na(means)]=quantile(means[!is.na(means)],0.5) #low number of males in steatosis!
#   sds=c();for (j in 9:28) {sds=append(sds,sd(SG0[,j],na.rm=TRUE))} #tähän jäätiin...
#   sds[is.na(sds)]=quantile(sds[!is.na(sds)],0.25) #low number of males in steatosis! and necroinfl.
#   error_lower=means-sds # error_lower[error_lower <= 0] = 0
#   error_upper=means+sds; error=sds
#   sample_data <- append(sample_data,data.frame(study=colnames(NAFLD[,9:28]),index=colnames(NAFLD[,9:28]),result=means,error=error))} # cate<- (str_extract(colnames(SG0[,9:28]), "[aA-zZ]+")) #https://stackoverflow.com/questions/29825537/group-categories-in-r-according-to-first-letters-of-a-string
# df=data.frame(sample_data)
# v2=data.frame(log(df$result.1/df$result))
# colnames(v2)=Group
# rownames(v2)=colnames(tv_all[,9:28])
# return(v2)}
# 
e=c()
# Outcome='Steatosis Grade';Group='All'; a=heattiin(NAFLD, Group,Outcome);Group='Female'; b=heattiin(NAFLD, Group,Outcome);Group='Male'; c=heattiin(NAFLD, Group,Outcome);e=append(e,cbind(a,b,c))
# Outcome='Fibrosis Stage';Group='All'; a=heattiin(NAFLD, Group,Outcome);Group='Female'; b=heattiin(NAFLD, Group,Outcome);Group='Male'; c=heattiin(NAFLD, Group,Outcome);e=append(e,cbind(a,b,c))
# Outcome='Necroinflammation';Group='All'; a=heattiin(NAFLD, Group,Outcome);Group='Female'; b=heattiin(NAFLD, Group,Outcome);Group='Male'; c=heattiin(NAFLD, Group,Outcome);e=append(e,cbind(a,b,c))
# Outcome='HOMA-IR';Group='All'; a=heattiin(NAFLD, Group,Outcome);Group='Female'; b=heattiin(NAFLD, Group,Outcome);Group='Male'; c=heattiin(NAFLD, Group,Outcome);e=append(e,cbind(a,b,c))
x=data.frame(n)
row.names(x)=colnames(tv_all[,9:28])

colnames(x)[1:3]=c("All_St.","Female_St.","Male_St.")
colnames(x) <- gsub(".1", "_Fib.", colnames(x))
colnames(x) <- gsub(".2", "_Nec.", colnames(x))
colnames(x) <- gsub(".3", "_HI", colnames(x))

groups=groups
groups[groups=="17a-OHP4"]="17aOH-P4"
op=groups[order(groups$Group),'Abbreviation']
op=op[op %in% row.names(x)]
x=x[op,]
brks_heatmap <- function(mat, color_palette){
     
       rng <- range(mat, na.rm = TRUE)
       lpal <- length(color_palette)
       
         c(seq(rng[1], 0, length.out=ceiling(lpal/2) + 1),
               seq(rng[2]/dim(mat)[1], rng[2], length.out=floor(lpal/2)))
 }
color_palette=brewer.pal(n = 11, name = 'RdBu')
color_palette=rev(color_palette)
library(grid)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(x,cluster_cols=FALSE,cluster_rows=FALSE,breaks = brks_heatmap(x, color_palette),color=color_palette,column_names_side = c("bottom"), angle_col = c("90"))
setHook("grid.newpage", NULL, "replace")
grid.text("Steatosis, Fibrosis, Necroinflammation, HOMA-IR", y=-0.07, x=0.4,gp=gpar(fontsize=16))
grid.text("Androgens, Estrogens, Gluc., Mineraloc., Progestogens ", x=-0.07, rot=90, gp=gpar(fontsize=16))

#this is with first(!!), use it
#is normal (non-log scale) image (and autoscaled), loge=1 is  log10 scale image (and raw data)
first=TRUE; ordera=c();
Group='All';name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); 
hel=pre_errors(NAFLD,Outcome,Group,name,ordera,loge,jap1,oute,first) 

#Afterwards:
first=FALSE
Group='All';name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date);
Group='Female'; name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); 
Group='Male'; name=paste("Forest plot of",Group, "Steroid Ratios in",Out,date); 







#This works with the normal raw data
# https://www.statology.org/cohens-d-in-r/
cohd=function(NAFLD, Group,Outcome) {
  if (Group=='Male') {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==2,]} else if (Group=='Female') 
  {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==1,]} else if (Group=='All') {NAFLDo=NAFLD}
  if (Outcome!='HOMA.IR') {for (i in 1:2) {if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] == 0,]} else if (i==2) {SG1=NAFLDo[NAFLDo[,Outcome] > 0,]}}} else
  for (i in 1:2) {if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] <= 1.5,]} else if (i==2) {SG1=NAFLDo[NAFLDo[,Outcome] > 1.5,]}}
  cd=c();for (i in 1:20) {cd=append(cd,cohen.d(SG1[,i+8],SG0[,i+8])$estimate)} 
  return(cd)}

d=c()
NAFLD=tv
Outcome='Steatosis.Grade.0.To.3';Group='All';a=cohd(NAFLD, Group,Outcome);Group='Female'; b=cohd(NAFLD, Group,Outcome);Group='Male'; c=cohd(NAFLD, Group,Outcome);d=cbind(d,cbind(a,b,c))
Outcome='Fibrosis.Stage.0.to.4';Group='All';a=cohd(NAFLD, Group,Outcome);Group='Female'; b=cohd(NAFLD, Group,Outcome);Group='Male'; c=cohd(NAFLD, Group,Outcome);d=cbind(d,cbind(a,b,c))
Outcome='Necroinflammation';Group='All';a=cohd(NAFLD, Group,Outcome);Group='Female'; b=cohd(NAFLD, Group,Outcome);Group='Male'; c=cohd(NAFLD, Group,Outcome);d=cbind(d,cbind(a,b,c))
Outcome='HOMA-IR';Group='All'; a=cohd(NAFLD, Group,Outcome);Group='Female'; b=cohd(NAFLD, Group,Outcome);Group='Male'; c=cohd(NAFLD, Group,Outcome);d=cbind(d,cbind(a,b,c))
d
rownames(d)=colnames(tv_all[,9:28])
colnames(d)=rep(c('All','Female','Male'),4)
write.csv(d,'cohens_d_tikka_v2asdfasfd.csv')
n=d

Group='All';
Outcome='Steatosis.Grade.0.To.3';
Outcome='Fibrosis.Stage.0.to.4';
Outcome='Necroinflammation';
Outcome='HOMA.IR'

heattiin=function(NAFLD,Group,Outcome){
  if (Group=='Male') {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==2,]} else if (Group=='Female') 
  {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==1,]} else if (Group=='All') {NAFLDo=NAFLD}
  sample_data=c()
  for (i in 1:2) {
    if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] == 0,]} else if (i==2) {SG0=NAFLDo[NAFLDo[,Outcome] > 0,]}
    means=c();for (j in 9:28) {means=append(means,mean(SG0[,j], na.rm=TRUE))}
    means[is.na(means)]=quantile(means[!is.na(means)],0.5) #low number of males in steatosis!
    sds=c();for (j in 9:28) {sds=append(sds,sd(SG0[,j],na.rm=TRUE))} #tähän jäätiin...
    sds[is.na(sds)]=quantile(sds[!is.na(sds)],0.25) #low number of males in steatosis! and necroinfl.
    error_lower=means-sds # error_lower[error_lower <= 0] = 0
    error_upper=means+sds; error=sds
    sample_data <- append(sample_data,data.frame(study=colnames(NAFLD[,9:28]),index=colnames(NAFLD[,9:28]),result=means,error=error))} # cate<- (str_extract(colnames(SG0[,9:28]), "[aA-zZ]+")) #https://stackoverflow.com/questions/29825537/group-categories-in-r-according-to-first-letters-of-a-string
  df=data.frame(sample_data)
  v2=data.frame(log(df$result.1/df$result))
  colnames(v2)=Group
  rownames(v2)=colnames(tv_all[,9:28])
  return(v2)}



e=c()
Outcome='Steatosis.Grade.0.To.3';Group='All'; a=heattiin(NAFLD, Group,Outcome);Group='Female'; b=heattiin(NAFLD, Group,Outcome);Group='Male'; c=heattiin(NAFLD, Group,Outcome);e=append(e,cbind(a,b,c))
Outcome='Fibrosis.Stage.0.to.4';Group='All'; a=heattiin(NAFLD, Group,Outcome);Group='Female'; b=heattiin(NAFLD, Group,Outcome);Group='Male'; c=heattiin(NAFLD, Group,Outcome);e=append(e,cbind(a,b,c))
Outcome='Necroinflammation';Group='All'; a=heattiin(NAFLD, Group,Outcome);Group='Female'; b=heattiin(NAFLD, Group,Outcome);Group='Male'; c=heattiin(NAFLD, Group,Outcome);e=append(e,cbind(a,b,c))
Outcome='HOMA.IR';Group='All'; a=heattiin(NAFLD, Group,Outcome);Group='Female'; b=heattiin(NAFLD, Group,Outcome);Group='Male'; c=heattiin(NAFLD, Group,Outcome);e=append(e,cbind(a,b,c))
x=data.frame(e)
row.names(x)=colnames(tv_all[,9:28])

groups=groups
groups[groups=="17a-OHP4"]="17aOH-P4"
op=groups[order(groups$Group),'Abbreviation']
op=op[op %in% row.names(x)]
x=x[op,]
brks_heatmap <- function(mat, color_palette){
  
  rng <- range(mat, na.rm = TRUE)
  lpal <- length(color_palette)
  
  c(seq(rng[1], 0, length.out=ceiling(lpal/2) + 1),
    seq(rng[2]/dim(mat)[1], rng[2], length.out=floor(lpal/2)))
}
color_palette=brewer.pal(n = 11, name = 'RdBu')
color_palette=rev(color_palette)
pheatmap(x,cluster_cols=FALSE,cluster_rows=FALSE,breaks = brks_heatmap(x, color_palette),
         color=color_palette,column_names_side = c("bottom"), angle_col = c("90"))



#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
# https://jokergoo.github.io/ComplexHeatmap-reference/book/integrate-with-other-packages.html
# color_palette=brewer.pal(n = 11, name = 'RdBu')
# # https://shanguangyu.com/articles/customized-color-key-for-pheatmap/

#https://stats.oarc.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
# https://lavaan.ugent.be/tutorial/mediation.html
#https://cran.r-project.org/web/packages/mma/mma.pdf
# install.packages('mma')

#Mediation Analysis
# date='tikka31123'; 
Outcome='Steatosis Grade';Group='All';
res1a=loop_med(Group,Outcome,date,tv_all);Group='Female'; res1f=loop_med(Group,Outcome,date,tv_all);Group='Male'; res1m=loop_med(Group,Outcome,date,tv_all)
Outcome='Fibrosis Stage';Group='All'; res2a=loop_med(Group,Outcome,date,tv_all);Group='Female'; res2f=loop_med(Group,Outcome,date,tv_all);Group='Male'; res2m=loop_med(Group,Outcome,date,tv_all)
Outcome='Necroinflammation';Group='All'; res3a=loop_med(Group,Outcome,date,tv_all);Group='Female'; res3f=loop_med(Group,Outcome,date,tv_all);Group='Male'; res3m=loop_med(Group,Outcome,date,tv_all)
Outcome='HOMA-IR';Group='All'; res4a=loop_med(Group,Outcome,date,tv_all);Group='Female'; res4f=loop_med(Group,Outcome,date,tv_all);Group='Male'; res4m=loop_med(Group,Outcome,date,tv_all)
# #https://stackoverflow.com/questions/19535996/avoid-rbind-cbind-conversion-from-numeric-to-factor
res1f=cbind(res1f,'Steatosis','Female');colnames(res1f)[9]='Case';colnames(res1f)[10]='Gender';
res1m=cbind(res1m,'Steatosis','Male');colnames(res1m)[9]='Case';colnames(res1m)[10]='Gender';
res2f=cbind(res2f,'Fibrosis','Female');colnames(res2f)[9]='Case';colnames(res2f)[10]='Gender';
res2m=cbind(res2m,'Fibrosis','Male');colnames(res2m)[9]='Case';colnames(res2m)[10]='Gender';
res3f=cbind(res3f,'Necroinflammation','Female');colnames(res3f)[9]='Case';colnames(res3f)[10]='Gender';
res3m=cbind(res3m,'Necroinflammation','Male');colnames(res3m)[9]='Case';colnames(res3m)[10]='Gender';
res4f=cbind(res4f,'HOMA-IR','Female');colnames(res4f)[9]='Case';colnames(res4f)[10]='Gender';
res4m=cbind(res4m,'HOMA-IR','Male');colnames(res4m)[9]='Case';colnames(res4m)[10]='Gender';# sum(is.na(tv))

# https://www.rdocumentation.org/packages/mediation/versions/4.5.0/topics/mediate

loop_med=function(Group,Outcome,date,tv_all) {
  if (Group=='Female') {cond=tv_all[,'Gender']==1} else if (Group=='Male') 
  {cond=tv_all[,'Gender']==2} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
  tv_red=tv_all[cond,]
  X <- tv_red[,71:(length(colnames(tv_all))-1)]#Standard values did not five erros
  M <- tv_red[,9:28]#
  Y <- tv_red[,Outcome] #"Steatosis.Grade.0.To.3"       "Fibrosis.Stage.0.to.4"       "Necroinflammation"            "HOMA-IR"   
  Data <- cbind(X, Y , M ); colnames(Data) <- gsub(" ", ".", colnames(Data))
  #multiple predictors (categorical and continuous predictors) https://rdrr.io/cran/mlma/man/data.org.html
  # https://cran.r-project.org/web/packages/mma/mma.pdf
  hiio=c();hoi=c(); 
  for (i in 1:dim(X)[2]) {hiio=append(hiio,paste("Data[, ",i , "]", sep=""));}
  for (j in (dim(X)[2]+2):dim(Data)[2]) {hoi=append(hoi,paste("Data[, ",j , "]", sep="")); }
  fmla1 <- as.formula(paste(paste(hoi, collapse= "+")," ~ ", paste(hiio, collapse= "+"))); e=paste('Data[,',which(colnames(Data)=='Y'),']');en=paste(c(hoi,hiio), collapse= "+")
  fmla2 <- as.formula(paste(e," ~ ", en))
  b = lm(fmla1, Data)  
  c = lm(fmla2, Data)  
  # Another linear model for the association (pheno ~ mediator + cause)
  med_out=c();res=c(); tmp=c();rn=c()
  # mediate(b, c, treat =  "Data[, 2]", mediator = "Data[, 12]", sims = 30) #ok...
  for (i in 1:dim(X)[2]) {
    for (j in 1:dim(M)[2]) {
      med_out = summary(mediate(b, c, treat = hiio[i], mediator = hoi[j],sims = 100)) #you need sims=100 min for the paper, maybe more like 1000... 10 was too little, but can get you results fast..
      tmp=c(med_out$d0, med_out$d0.p, med_out$z0, med_out$z0.p, med_out$tau.coef, med_out$tau.p, med_out$n0, med_out$n0.p,med_out$d1, med_out$d1.p, med_out$z1, med_out$z1.p, med_out$n1, med_out$n1.p) #med_out$d1,med_out$d1.p,
      #med_out$z0.p is for the ADE and  med_out$d0.p for the ACME
      res <- rbind(res,tmp);rn=append(rn,paste(colnames(Data)[i],colnames(Data)[j+11]))
      remove(tmp)}}
  rownames(res)=rn
  colnames(res)=c('d0', 'd0.p', 'z0', 'z0.p', 'tau.coef', 'tau.p', 'n0', 'n0.p','d1', 'd1.p', 'z1', 'z1.p', 'n1', 'n1.p') #d0 and d1 are the same as.. 'd1', 'd1.p', 
  res=res[order(res[,2]),]
  # write.xlsx(res, file = paste(Outcome,Group,'loop_mediation_',date,'.xlsx'), append = FALSE)
  # write.csv(res,paste(Outcome,Group,'loop_mediation_',date,'.csv'))
  return(res)}


# https://rdrr.io/cran/mediation/src/R/mediate.R
# https://m-clark.github.io/mixed-models-with-R/random_intercepts.html
# https://rdrr.io/cran/mediation/src/R/mediate.R
# https://library.virginia.edu/data/articles/introduction-to-mediation-analysis
library(ggsankey)
# rtote=rtot
rtot=rbind(res1f,res1m,res2f,res2m,res3f,res3m,res4f,res4m)
colnames(rtot)=c(colnames(rtot)[1:8],'d1','d1.p',colnames(rtot)[11:14],colnames(rtot)[9:10])
# rtot_pb=rtot
# rtot_pb1=rtot
# rtot_nok=rtot
# rtot_pbr=rtot
# rtot_smooth_10=rtot
rtot_smooth_100=rtot
# write.csv(data.frame(rtot_smooth_10),'mediationbackup_rtot_smooth_10_all_tikka271023.csv')
write.csv(data.frame(rtot_smooth_100),'mediationbackup_rtot_smooth_100_all_tikka291023.csv')
# rtot=rtot_smooth_100

sankeys=function(med,rtot,cut,Group,size,div) {
  if (size=='All') {rt2=rtot} else if (size!='All')  {rt2=rtot[rtot[,'d0']>size,]}
  if (med=='Mediated') {rt2=rt2[rt2[,'d0.p']<pcut,]; rt2=rt2[rt2[,'z0.p']>pcut,]} else
  {rt2=rt2[rt2[,'d0.p']>pcut,]; rt2=rt2[rt2[,'z0.p']<pcut,]} 
  hoi=c();for (i in 1:dim(rt2)[1]) {hoi=append(hoi,scan(text=rownames(rt2)[i], what=""))}
  hoi=as.data.frame(matrix(hoi, ncol = 2,  byrow = TRUE), stringsAsFactors = FALSE)
  hoi[,3]=rt2[,15]
  hoi[,4]=rt2[,16]
  hoi[,5]=rt2[,1]
  colnames(hoi)=c('Contaminants','Steroids','Case','Gender','Mediation')# median(as.numeric(hoi[,'Mediation']),breaks=50)# hist(as.numeric(hoi[,'Mediation']),breaks=50)
  # https://stats.stackexchange.com/questions/271155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct
  # https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation
  hoi=hoi[,1:4]
  if (Group=='All') {hoi=hoi} else if  
  (Group=='Female') {hoi=hoi[hoi[,'Gender']=='Female',]} else if 
  (Group=='Male') {hoi=hoi[hoi[,'Gender']=='Male',]};hoi
  df2 <- hoi %>%make_long('Gender',"Contaminants",'Steroids',"Case" )
  jpeg(paste(med,Group,div ,"e.jpg"), width = 7500, height = 11000, quality = 100,pointsize = 16, res=1000);
  print(
    ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
      geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
      scale_fill_viridis_d() +
      theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank())
  );dev.off()}

# https://stackoverflow.com/questions/64829727/mediation-r-package-p-values-workaround-to-get-more-significant-digits
# https://bookdown.org/mike/data_analysis/mediation.html

hist(as.numeric(rtot[,'d0.p']),breaks=50)
pcut=0.4
med='Mediated'
size=0
Group='All'
div='Yeshaanet4'
sankeys(med,rtot,cut,Group,size,div) 

#Load MaAsLin 2 package into the R environment
# input_data = system.file("extdata", "HMP2_taxonomy - Copy.tsv", package="Maaslin2");input_data # The abundance table file
# input_metadata = system.file("extdata", "HMP2_metadata - Copy.tsv", package="Maaslin2") 
# The metadata table file, ("https://raw.githubusercontent.com/biobakery/biobakery_workflows/master/examples/tutorial/stats_vis/input/pathabundance_relab.tsv", "./pathabundance_relab.tsv")

#Load MaAsLin 2 package into the R environment
# C:\Users\patati\AppData\Local\Programs\R\R-4.3.1\library\Maaslin2\extdata\
input_data = system.file("extdata", "HMP2_taxonomy.tsv", package="Maaslin2"); input_data # The abundance table file... new one  :) 15.1.15
input_metadata = system.file("extdata", "HMP2_metadata.tsv", package="Maaslin2") #


# pop=c(colnames(tv_all)[2:8],x5)#c('Gender','AGE','BMI','Steatosis Grade','Fibrosis Grade','Necroinflammation','HOMA-IR','PFAS')
# df_input_metadatae=tv_all[,pop]
# # write.csv(df_input_metadatae,'HMP2_metadata_v2e.csv')
# write.table(df_input_metadatae, file='HMP2_metadata_v4.tsv', quote=FALSE, sep='\t')

# C:\Users\patati\AppData\Local\Programs\R\R-4.3.1\library\Maaslin2\extdata\ https://stackoverflow.com/questions/17108191/how-to-export-proper-tsv
# The metadata table file, ("https://raw.githubusercontent.com/biobakery/biobakery_workflows/master/examples/tutorial/stats_vis/input/pathabundance_relab.tsv", "./pathabundance_relab.tsv")
df_input_data = read.table(file = input_data,header  = TRUE,sep = "\t", row.names  = 1, stringsAsFactors = FALSE)
df_input_metadata = read.table(file  = input_metadata, header           = TRUE, sep              = "\t", row.names        = 1,stringsAsFactors = FALSE) ##Check the the names, no # or - are needed..
#Scope
ui=1:20; df_input_data=df_input_data[,ui]
# As I gather now, this is not ok: df_input_data=cbind(df_input_data,df_input_metadata[,c('AGE' ,'BMI','Gender')]) # ui=1:23 #steroids +others
#Genders:
# df_input_data=cbind(df_input_data,df_input_metadata[,c('AGE' ,'BMI')])
df_input_dataf=df_input_data[df_input_metadata[,'Gender']==1,ui]
df_input_metadataf=df_input_metadata[df_input_metadata[,'Gender']==1, ]
df_input_datam=df_input_data[df_input_metadata[,'Gender']==2, ui]
df_input_metadatam=df_input_metadata[df_input_metadata[,'Gender']==2, ]
df_input_data[1:2, 1:20]; 
df_input_metadata[1:2, ];
# tv_all[1:2,9:28]; #tv_all[1:5,1:8];min(df_input_metadata[,'HOMAIR'])
# hist(as.numeric(df_input_metadata[,'HOMAIR']),xlim=c(0,18), breaks=40)

#Steatosis.Grade ei hitsi näiden nimien kaa sarakkeissa... :)

#These models work equally well without the zero reference. I am assuming that the reference for the age, bmi and gender are from the min values.

#Steatosis models:
fit_data_sa = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform ="NONE",
                       output         = "Steatosis_all",
                       fixed_effects  = c('Steatosis','AGE' ,'BMI','Gender'), #'Gender', 
                       reference      = c("Steatosis,0")) #fit_data_sa$results
#,random_effects=c() ok, not I think to be evaluated, since at least subjects? https://github.com/biobakery/biobakery/wiki/maaslin2
# https://lost-stats.github.io/
# https://julianfaraway.github.io/brinlabook/chaglmm.html
fit_data_sf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform =    "NONE",
                       output         = "Steatosis_female",
                       fixed_effects  = c('Steatosis','AGE' ,'BMI'), #'Gender',
                       reference      = c("Steatosis,0")) #,random_effects=c() ok
fit_data_sm = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform =    "NONE",
                       output         = "Steatosis_male_adj",fixed_effects  = c('Steatosis','AGE' ,'BMI'), #'Gender',
                       reference      = c("Steatosis,0")) #,random_effects=c()

#Fibrosis Models
fit_data_fa = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform =    "NONE",
                       output         = "Fibrosis_all_adj",fixed_effects  = c('Fibrosis','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("Fibrosis,0")) #,random_effects=c()
fit_data_ff = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform =    "NONE",
                       output         = "Fibrosis_female_adj",fixed_effects  = c('Fibrosis','AGE' ,'BMI'), #'Gender', #need ti cgecj what is the effect...
                       reference      = c("Fibrosis,0")) #,random_effects=c()
fit_data_fm = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform =    "NONE",
                       output         = "Fibrosis_male_adj",fixed_effects  = c('Fibrosis','AGE' ,'BMI'), #'Gender',
                       reference      = c("Fibrosis,0")) #,random_effects=c()

#HOMAIR Models
fit_data_ha = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform = "NONE",
                       output         = "HOMAIR_all_adj",fixed_effects  = c('HOMAIR','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("HOMAIR")) #,random_effects=c()
fit_data_hf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform =    "NONE",
                       output         = "HOMAIR_female_adj",fixed_effects  = c('HOMAIR','AGE' ,'BMI'), #'Gender',
                       reference      = c("HOMAIR")) #,random_effects=c()
#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_hm = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam,
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
                       output         = "HOMAIR_male_adj",fixed_effects  = c('HOMAIR','AGE' ,'BMI'), #'Gender',
                       reference      = c("HOMAIR")) #,random_effects=c()

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_na = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
                       output         = "Necroinflammation_all_adj",fixed_effects  = c('Necroinflammation','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("Necroinflammation,0")) #,random_effects=c()
fit_data_naf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                        normalization  = "NONE",transform = "NONE", 
                        min_prevalence = 0.05,min_abundance=0.0,
                        output         = "Necroinflammation_female_adj",fixed_effects  = c('Necroinflammation','AGE' ,'BMI'), #'Gender',
                        reference      = c("Necroinflammation,0")) #,random_effects=c()
fit_data_nam = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                        min_prevalence = 0.05,min_abundance=0.0,
                        normalization  = "NONE",transform = "NONE",
                        output         = "Necroinflammation_male_adj",fixed_effects  = c('Necroinflammation','AGE' ,'BMI'), #'Gender',
                        reference      = c("Necroinflammation,0")) #,random_effects=c()

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_pa = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
                       output         = "PFAS_all_adj",fixed_effects  = c('PFAS','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("PFAS,0")) #,random_effects=c()
fit_data_paf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                        normalization  = "NONE",transform = "NONE", 
                        min_prevalence = 0.05,min_abundance=0.0,
                        output         = "PFAS_female_adj",fixed_effects  = c('PFAS','AGE' ,'BMI'), #'Gender',
                        reference      = c("PFAS,0")) #,random_effects=c()
fit_data_pam = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                        min_prevalence = 0.05,min_abundance=0.0,
                        normalization  = "NONE",transform = "NONE",
                        output         = "PFAS_male_adj",fixed_effects  = c('PFAS','AGE' ,'BMI'), #'Gender',
                        reference      = c("PFAS,0")) #,random_effects=c()



# PFHpA     PFHxA PFHxA_Branched      PFHxS      PFNA       PFOA      PFOS

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_PFHpA_a = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                       min_prevalence = 0.05,min_abundance=0.0,
                       normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
                       output         = "PFHpA_all_adj",fixed_effects  = c('PFHpA','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("PFHpA,0")) #,random_effects=c()
fit_data_PFHpA_f = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                        normalization  = "NONE",transform = "NONE", 
                        min_prevalence = 0.05,min_abundance=0.0,
                        output         = "PFHpA_female_adj",fixed_effects  = c('PFHpA','AGE' ,'BMI'), #'Gender',
                        reference      = c("PFHpA,0")) #,random_effects=c()
fit_data_PFHpA_m = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                        min_prevalence = 0.05,min_abundance=0.0,
                        normalization  = "NONE",transform = "NONE",
                        output         = "PFHpA_male_adj",fixed_effects  = c('PFHpA','AGE' ,'BMI'), #'Gender',
                        reference      = c("PFHpA,0")) #,random_effects=c()

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_PFHxA_a = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFHxA_all_adj",fixed_effects  = c('PFHxA','AGE' ,'BMI','Gender'), #'Gender',
                            reference      = c("PFHxA,0")) #,random_effects=c()
fit_data_PFHxA_f = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                            normalization  = "NONE",transform = "NONE", 
                            min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFHxA_female_adj",fixed_effects  = c('PFHxA','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFHxA,0")) #,random_effects=c()
fit_data_PFHxA_m = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE",
                            output         = "PFHxA_male_adj",fixed_effects  = c('PFHxA','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFHxA,0")) #,random_effects=c()

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_PFHxA_Branched_a = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFHxA_Branched_all_adj",fixed_effects  = c('PFHxA_Branched','AGE' ,'BMI','Gender'), #'Gender',
                            reference      = c("PFHxA_Branched,0")) #,random_effects=c()
fit_data_PFHxA_Branched_f = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                            normalization  = "NONE",transform = "NONE", 
                            min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFHxA_Branched_female_adj",fixed_effects  = c('PFHxA_Branched','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFHxA_Branched,0")) #,random_effects=c()
fit_data_PFHxA_Branched_m = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE",
                            output         = "PFHxA_Branched_male_adj",fixed_effects  = c('PFHxA_Branched','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFHxA_Branched,0")) #,random_effects=c()

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_PFHxS_a = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFHxS_all_adj",fixed_effects  = c('PFHxS','AGE' ,'BMI','Gender'), #'Gender',
                            reference      = c("PFHxS,0")) #,random_effects=c()
fit_data_PFHxS_f = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                            normalization  = "NONE",transform = "NONE", 
                            min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFHxS_female_adj",fixed_effects  = c('PFHxS','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFHxS,0")) #,random_effects=c()
fit_data_PFHxS_m = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE",
                            output         = "PFHxS_male_adj",fixed_effects  = c('PFHxS','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFHxS,0")) #,random_effects=c()

# PFNA       PFOA      PFOS

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_PFNA_a = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFNA_all_adj",fixed_effects  = c('PFNA','AGE' ,'BMI','Gender'), #'Gender',
                            reference      = c("PFNA,0")) #,random_effects=c()
fit_data_PFNA_f = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                            normalization  = "NONE",transform = "NONE", 
                            min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFNA_female_adj",fixed_effects  = c('PFNA','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFNA,0")) #,random_effects=c()
fit_data_PFNA_m = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE",
                            output         = "PFNA_male_adj",fixed_effects  = c('PFNA','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFNA,0")) #,random_effects=c()

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_PFNA_a = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFNA_all_adj",fixed_effects  = c('PFNA','AGE' ,'BMI','Gender'), #'Gender',
                            reference      = c("PFNA,0")) #,random_effects=c()
fit_data_PFNA_f = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                            normalization  = "NONE",transform = "NONE", 
                            min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFNA_female_adj",fixed_effects  = c('PFNA','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFNA,0")) #,random_effects=c()
fit_data_PFNA_m = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE",
                            output         = "PFNA_male_adj",fixed_effects  = c('PFNA','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFNA,0")) #,random_effects=c()

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_PFOS_a = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFOS_all_adj",fixed_effects  = c('PFOS','AGE' ,'BMI','Gender'), #'Gender',
                            reference      = c("PFOS,0")) #,random_effects=c()
fit_data_PFOS_f = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
                            normalization  = "NONE",transform = "NONE", 
                            min_prevalence = 0.05,min_abundance=0.0,
                            output         = "PFOS_female_adj",fixed_effects  = c('PFOS','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFOS,0")) #,random_effects=c()
fit_data_PFOS_m = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
                            min_prevalence = 0.05,min_abundance=0.0,
                            normalization  = "NONE",transform = "NONE",
                            output         = "PFOS_male_adj",fixed_effects  = c('PFOS','AGE' ,'BMI'), #'Gender',
                            reference      = c("PFOS,0")) #,random_effects=c()



fit_data_sa; sa=data.frame(fit_data_sa[[1]])
fit_data_sf ;sf=data.frame(fit_data_sf[[1]])
fit_data_sm ; sm=data.frame(fit_data_sm[[1]])
fit_data_fa ;fa=data.frame(fit_data_fa[[1]])
fit_data_ff;ff=data.frame(fit_data_ff[[1]])
fit_data_fm;fm=data.frame(fit_data_fm[[1]])
fit_data_ha ;ha=data.frame(fit_data_ha[[1]])
fit_data_hf ;hf=data.frame(fit_data_hf[[1]])
fit_data_hm ;hm=data.frame(fit_data_ha[[1]])
# hm=read.csv("C:/Users/patati/Desktop/TurkuOW/RWork/HOMAIR_male_adj/all_results.tsv", header = TRUE,sep="\t");
# htm=hm[hm$metadata=='HOMAIR',]# hta=ha[ha$metadata=='HOMAIR',]
fit_data_na; na=data.frame(fit_data_na[[1]]) 
fit_data_naf  ;nf=data.frame(fit_data_naf[[1]]) 
fit_data_nam; nm=data.frame(fit_data_nam[[1]]) #malelta puuttu joku steroidi... ja piti lisätä:min_prevalence = -0.000001,min_abundance=-0.000001,
# sa;sf;sm;fa;ff;fm;ha;hf;hm;na;nf;nm

dim(sa[sa[,2]=='Steatosis',])
dim(sf[sf[,2]=='Steatosis',])#? without minimun prevalnce, I do not have all the steroids... particularly T/Epi-T here...
dim(sm[sm[,2]=='Steatosis',])#?

dim(fa[fa[,2]=='Fibrosis',])
dim(ff[ff[,2]=='Fibrosis',])#? 
dim(fm[fm[,2]=='Fibrosis',])#

dim(na[na[,2]=='Necroinflammation',])
dim(nf[nf[,2]=='Necroinflammation',])#? 
dim(nm[nm[,2]=='Necroinflammation',])#

dim(ha[ha[,2]=='HOMAIR',])
dim(hf[hf[,2]=='HOMAIR',])#? 
dim(hm[hm[,2]=='HOMAIR',])#

# > dim(sa[sa[,2]=='Steatosis',])[1] 20 10
# > dim(sf[sf[,2]=='Steatosis',])[1] 19 10 #checkt the min prevalence closer to zero
# > dim(sm[sm[,2]=='Steatosis',])#?[1] 20 10
# > dim(fa[fa[,2]=='Fibrosis',])[1] 20 10
# > dim(ff[ff[,2]=='Fibrosis',])#? [1] 19 10
# > dim(fm[fm[,2]=='Fibrosis',])#[1] 20 10
# > dim(na[na[,2]=='Necroinflammation',])[1] 20 10
# > dim(nf[nf[,2]=='Necroinflammation',])#? [1] 19 10
# > dim(nm[nm[,2]=='Necroinflammation',])#[1] 20 10
# > dim(ha[ha[,2]=='HOMAIR',])[1] 20 10
# > dim(hf[hf[,2]=='HOMAIR',])#? [1] 19 10
# > dim(hm[hm[,2]=='HOMAIR',])#[1] 20 10

totis=rbind('Steatosis_all',sa,'Steatosis_female',sf,'Steatosis_male',sm,
        'Fibrosis_all',fa,'Fibrosis_female',ff,'Fibrosis_male',fm,
        'Necroinflammation_all',na,'Necroinflammation_female',nf,'Necroinflammation_male',nm,
        'HOMAIR_all',ha,'HOMAIR_female',hf,'HOMAIR_male',hm)

# totis=data.frame(totis)

# totis2=totis

write.xlsx(totis, file ='lms_tikka16324_nonadj.xlsx', append = FALSE, row.names = TRUE) 

#Steatosis models:
# fit_data_sa = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
#                        min_prevalence = 0.05,min_abundance=0.0,
#                        normalization  = "NONE",transform ="NONE",
#                        output         = "Steatosis_all",
#                        fixed_effects  = c('Steatosis'), #'Gender', 
#                        reference      = c("Steatosis,0")) #fit_data_sa$results
# fit_data_sf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
#                        min_prevalence = 0.05,min_abundance=0.0,
#                        normalization  = "NONE",transform =    "NONE",
#                        output         = "Steatosis_female",
#                        fixed_effects  = c('Steatosis'), #'Gender',
#                        reference      = c("Steatosis,0")) #,random_effects=c() ok
# fit_data_sm = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
#                        min_prevalence = 0.05,min_abundance=0.0,
#                        normalization  = "NONE",transform =    "NONE",
#                        output         = "Steatosis_male",fixed_effects  = c('Steatosis'), #'Gender',
#                        reference      = c("Steatosis,0")) #,random_effects=c()
# #Fibrosis Models
# fit_data_fa = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
#                        min_prevalence = 0.05,min_abundance=0.0,
#                        normalization  = "NONE",transform =    "NONE",
#                        output         = "Fibrosis_all",fixed_effects  = c('Fibrosis'), #'Gender',
#                        reference      = c("Fibrosis,0")) #,random_effects=c()
# fit_data_ff = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
#                        min_prevalence = 0.05,min_abundance=0.0,
#                        normalization  = "NONE",transform =    "NONE",
#                        output         = "Fibrosis_female",fixed_effects  = c('Fibrosis'), #'Gender', #need ti cgecj what is the effect...
#                        reference      = c("Fibrosis,0")) #,random_effects=c()
# fit_data_fm = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
#                        min_prevalence = 0.05,min_abundance=0.0,
#                        normalization  = "NONE",transform =    "NONE",
#                        output         = "Fibrosis_male",fixed_effects  = c('Fibrosis'), #'Gender',
#                        reference      = c("Fibrosis,0")) #,random_effects=c()
# #HOMAIR Models
# fit_data_ha = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
#                        min_prevalence = 0.05,min_abundance=0.0,
#                        normalization  = "NONE",transform = "NONE",
#                        output         = "HOMAIR_all",fixed_effects  = c('HOMAIR'), #'Gender',
#                        reference      = c("HOMAIR")) #,random_effects=c()
# fit_data_hf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
#                        min_prevalence = 0.05,min_abundance=0.0,
#                        normalization  = "NONE",transform =    "NONE",
#                        output         = "HOMAIR_female",fixed_effects  = c('HOMAIR'), #'Gender',
#                        reference      = c("HOMAIR")) #,random_effects=c()
# #note the prevalences and abudances ahve to be little bit below zero to get all the variables:
# fit_data_hm = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam,
#                        min_prevalence = 0.05,min_abundance=0.0,
#                        normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
#                        output         = "HOMAIR_male",fixed_effects  = c('HOMAIR'), #'Gender',
#                        reference      = c("HOMAIR")) #,random_effects=c()
# #note the prevalences and abudances have to be little bit below zero to get all the variables:
# fit_data_na = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata, 
#                        min_prevalence = 0.05,min_abundance=0.0,
#                        normalization  = "NONE",transform = "NONE", #min_prevalence = 0.05,min_abundance=0.0,
#                        output         = "Necroinflammation_all",fixed_effects  = c('Necroinflammation'), #'Gender',
#                        reference      = c("Necroinflammation,0")) #,random_effects=c()
# fit_data_naf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf, 
#                         normalization  = "NONE",transform = "NONE", 
#                         min_prevalence = 0.05,min_abundance=0.0,
#                         output         = "Necroinflammation_female",fixed_effects  = c('Necroinflammation'), #'Gender',
#                         reference      = c("Necroinflammation,0")) #,random_effects=c()
# fit_data_nam = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam, 
#                         min_prevalence = 0.05,min_abundance=0.0,
#                         normalization  = "NONE",transform = "NONE",
#                         output         = "Necroinflammation_male",fixed_effects  = c('Necroinflammation'), #'Gender',
#                         reference      = c("Necroinflammation,0")) #,random_effects=c()

sa[sa[,2]=='Steatosis',][1:12,];
sm[sm[,2]=='Steatosis',][1:12,]
sf[sf[,2]=='Steatosis',]

fa[fa[,2]=='Fibrosis',][1:12,];fm[fm[,2]=='Fibrosis',][1:12,]
hm[hm[,2]=='HOMAIR',][1:12,]

# https://lost-stats.github.io/Model_Estimation/OLS/fixed_effects_in_linear_regression.html
# One variable was 'missing' in males without min_prevalence = -0.000001,min_abundance=-0.000001,
# https://huttenhower.sph.harvard.edu/maaslin/
# https://datatofish.com/create-dataframe-in-r/




ok_forestplot=function(Outcome, Group, name, data_frame,grouping,ordera,first) { 
  valuees=data_frame;sample_data=valuees
  sample_data=sample_data[!sample_data$metadata %in% c('AGE', 'BMI', 'Gender'),] # sample_data=sample_data[sample_data$metadata %in% c('AGE'),]
  gn=groups[,c('Group','Abbreviation')]
  gn=gn[gn[,'Abbreviation']!='F',]
  gn[order(data.frame(gn[,'Abbreviation'])[,1]),]
  sample_data=sample_data[,!colnames(sample_data) %in% c('metadata','value','name','N.not.zero','N','N.not.0')]
  names(sample_data)[names(sample_data) == "feature"] <- "name" #"feature"] <- "namee"
  sample_data=sample_data[!sample_data[,'name' ] %in% c('Gender','AGE','BMI'),]
  rownames(sample_data)=1:dim(sample_data)[1]
  # sample_data=sample_data[order(sample_data[,'name']),] #stud feature
  # ch=c("T.Epi.T","X11.KA4",   "X11.KDHT",  "X11.KT",    "X11b.OHA4", "X17a.OHP5", "X17aOH.P4")
  # sample_data[,'name' ] [sample_data[,'name' ] %in% ch] =c('T/Epi-T','11-KA4','11-KDHT','11-KT','11b-OHA4','17a-OHP5','17aOH-P4')
  # # sample_data=sample_data[order(sample_data[,'name']),] #stud feature
  ch=c("T.Epi.T","X11.KA4",   "X11.KDHT",  "X11.KT",    "X11b.OHA4", "X17a.OHP5", "X17.aOHP4")
  sample_data[,'name' ][sample_data[,'name' ]=='X17aOH.P4']='X17.aOHP4'
  # # for (i in 1:length(sample_data[,'name' ])) {
  # #   for (j in 1:length(ch)) {#   if (sample_data[i,'name' ]==ch[j]) }}
  # # if (all((sample_data[,'name' ] %in% ch)[14:20])==FALSE) {break}
  # # sample_data[,'name' ] [sample_data[,'name' ] %in% ch] =c('T/Epi-T','11-KA4','11-KDHT','11-KT','11b-OHA4','17a-OHP5','17aOH-P4')
  sample_data[,'name' ] [sample_data[,'name' ] %in% ch] <- gsub("\\.", "-", sample_data[,'name' ] [sample_data[,'name' ] %in% ch])
  sample_data[,'name' ]  <- gsub("X", "", sample_data[,'name' ] )
  sample_data[,'name' ][sample_data[,'name' ]=='T-Epi-T']='T/Epi-T'
  sample_data=sample_data[order(sample_data[,'name']),]
  
  sample_data=cbind(sample_data,gn[order(data.frame(gn[,'Abbreviation'])[,1]),])
  sample_data[,'Significance']=sample_data[,'pval'] < 0.1
  sample_data[,'Significance'][sample_data[,'Significance']==TRUE]='Yes'
  sample_data[,'Significance'][sample_data[,'Significance']==FALSE]='No'
  sample_data[,'Color']=sample_data[,'pval'] < 0.1
  sample_data[,'Color'][sample_data[,'Color']==TRUE]='blue'
  sample_data[,'Color'][sample_data[,'Color']==FALSE]='grey'# Rename column where names is "..."
  
  # name    index      result  error_lower    error_upper      error        pval              Group Abbreviation Significance Color
  plot=forestplot(df = sample_data, #drive tba_example_v3_oh_tikka17823.R if not working via 'x' error
                  estimate = coef, #coef, # result
                  se=stderr/2,#abs(error_lower -error_upper)/2, #abs(error-nerror)/4, #stderr,#abs(error_lower -error_upper)/4 
                  pvalue = pval, #this makes the significant value..
                  psignif = 0.1,
                  xlim = c(-0.6, 0.6),
                  xlab = "Linear Model Coefficients (SE)",
                  ylab='Steroid Groups',
                  title='',#name,
                  colour = Significance ) +#,colour = Significance 
    ggforce::facet_col(
      facets = ~Group,
      scales = "free_y",
      space = "free",
      strip.position='left');plot #standard error s/sqrt(n)
  
  # #It is not the same as previously factor(...) as such
  if (sum(sample_data[,'Significance']=='Yes')==20) {hp=c('blue','blue')} else {hp=c('#999999','blue')}
  if (Group=='All' & first==TRUE) {ordera=sample_data$name[order(sample_data$coef)] #this should be ok
  plot[["data"]][["name"]]=factor(plot[["data"]][["name"]], levels = ordera)} else if
  (Group=='All' & first==FALSE) {plot[["data"]][["name"]]=factor(plot[["data"]][["name"]], levels = ordera)} else if
  (Group=='female') {plot[["data"]][["name"]]=factor(plot[["data"]][["name"]], levels = ordera)} else if
  (Group=='male') {plot[["data"]][["name"]]=factor(plot[["data"]][["name"]], levels = ordera)}#this does it.. yesh :)
  plot$layers[[1]]$aes_params$odd <- "#00000000" #https://stackoverflow.com/questions/71745719/how-to-control-stripe-transparency-using-ggforestplot-geom-stripes
  
  sample_data$Group2=sample_data$Group
  sample_data <- transform(sample_data,Group2 = as.numeric(as.factor(Group2)))
  sample_data$facet_fill_color <- c("red", "green", "blue", "yellow", "brown")[sample_data$Group2]
  jopon=plot  +theme(axis.text.y=element_blank()) +theme_classic2();   #theme(axis.text.y = element_text(lineheight=.05));
  jopon2=jopon+geom_point(aes(colour = factor(Significance)),colour = sample_data[,'Color']) +  #facet_wrap2(~ Group, strip = strip)+ #axes='x',strip.position='left'
    #element_rect(colour = c("black",'yellow','red','green','brown'), fill = c("black",'yellow','red','green','brown')))+
    # geom_col(aes(x = result, y = name, fill = Group))+scale_fill_manual(values = c('red','blue','black','yellow','orange'))+ #geom_vline(xintercept = 1,size = 1)+
    scale_color_manual(values=hp)+theme(legend.position = "none")+theme(strip.text.y = element_text(size=-Inf)); jopon2 #ggtext::element_markdown(size = 12)
  # https://stackoverflow.com/questions/10547487/remove-facet-wrap-labels-completely
  
  g <- ggplot_gtable(ggplot_build(jopon2))
  stripr <- which(grepl('strip-l', g$layout$name))
  fills <- c("red","green","blue","yellow",'brown')
  k <- 1
  for (i in stripr) {j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1}; grid::grid.draw(g)
  jpeg(paste("Forest plot of Grouped Steroids Showing Ratios of ",Outcome,'with',Group,date,".jpg"), 
       width = 6500, height = 11000, quality = 100,pointsize = 16, res=1000); print(grid::grid.draw(g));dev.off();return(ordera)} #remember print :)

median(tv[1:20,c('11-KT')])


#Maaslin2, https://github.com/biobakery/biobakery/wiki/maaslin2, "https://raw.githubusercontent.com/biobakery/biobakery_workflows/master/examples/tutorial/stats_vis/input/pathabundance_relab.tsv", "./pathabundance_relab.tsv")
Outcome='HOMA-IR'
Group='All';data_frame=ha
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples')
grouping=groups #elsewhere defined

# Outcome='Steatosis'
Group='male';data_frame=hm
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,ordera=hels)


data_frame=sf

# ok_forestplot=function(Outcome, Group, name, data_frame,grouping,ordera,first)
# See below the ok one...


#Tryings:
Outcome='Steatosis';Group='All';first=TRUE
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
data_frame=sa
grouping=groups #elsewhere defined
hels=ok_forestplot(Outcome, Group, name, data_frame,grouping,ordera,first);
first=FALSE

Outcome='Steatosis'
Group='male';data_frame=sm
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first)

Outcome='Steatosis'
Group='female';data_frame=sf
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first)

Outcome='Fibrosis'
Group='All';data_frame=fa
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first)

Outcome='Fibrosis'
Group='female';data_frame=ff
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first)

Outcome='Fibrosis'
Group='male';data_frame=fm
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first)

Outcome='Necroinflammation'
Group='All';data_frame=na
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first) #helh=

Outcome='Necroinflammation'
Group='female';data_frame=nf
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first)

Outcome='Necroinflammation'
Group='male';data_frame=nm
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first)

Outcome='HOMA-IR'
Group='All';data_frame=ha
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first)#helh=

Outcome='HOMA-IR'
Group='female';data_frame=hf
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first)

Outcome='HOMA-IR'
Group='male';data_frame=hm
name=paste("Forest plot of Steroid Ratios in", Outcome, 'with',Group,'samples', date)
ok_forestplot(Outcome, Group, name, data_frame,grouping,hels,first)

# ha[order(hm[,'metadata']),]

# 
# #Steatosis models:
fit_data_sa = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata,
                       normalization  = "NONE",transform ="NONE", #
                       min_prevalence = 0.00001,
                       output         = "Steatosis_all",
                       fixed_effects  = c('Steatosis','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("Steatosis,0"))
# #,random_effects=c() ok, not I think to be evaluated, since at least subjects? https://github.com/biobakery/biobakery/wiki/maaslin2
# # https://lost-stats.github.io/
# # https://julianfaraway.github.io/brinlabook/chaglmm.html
fit_data_sf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf,
                       normalization  = "NONE",transform =    "NONE", #min_prevalence = 0,
                       min_prevalence = 0.00001,
                       output         = "Steatosis_female",
                       fixed_effects  = c('Steatosis','AGE' ,'BMI'), #'Gender',
                       reference      = c("Steatosis,0")) #,random_effects=c() ok
fit_data_sm = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam,
                       normalization  = "NONE",transform =    "NONE", #min_prevalence = -0.000001,min_abundance=-0.000001,
                       output         = "Steatosis_male_adj",fixed_effects  = c('Steatosis','AGE' ,'BMI'), #'Gender',
                       reference      = c("Steatosis,0")) #,random_effects=c()

#Fibrosis Models
fit_data_fa = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata,
                       normalization  = "NONE",transform =    "NONE",#min_prevalence = 0,
                       min_prevalence = 0.00001,
                       output         = "Fibrosis_all_adj",fixed_effects  = c('Fibrosis','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("Fibrosis,0")) #,random_effects=c()
fit_data_ff = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf,
                       normalization  = "NONE",transform =    "NONE",#min_prevalence = 0,
                       min_prevalence = 0.00001,
                       output         = "Fibrosis_female_adj",fixed_effects  = c('Fibrosis','AGE' ,'BMI'), #'Gender', #need ti cgecj what is the effect...
                       reference      = c("Fibrosis,0")) #,random_effects=c()
fit_data_fm = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam,
                       min_prevalence = 0.00001,
                       normalization  = "NONE",transform =    "NONE", #min_prevalence = -0.000001,min_abundance=-0.000001,
                       output         = "Fibrosis_male_adj",fixed_effects  = c('Fibrosis','AGE' ,'BMI'), #'Gender',
                       reference      = c("Fibrosis,0")) #,random_effects=c()

#HOMAIR Models
fit_data_ha = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata,
                       normalization  = "NONE",transform = "NONE", #
                       # min_prevalence = 0, #
                       min_prevalence = 0.00001,
                       output         = "HOMAIR_all_adj",fixed_effects  = c('HOMAIR','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("HOMAIR")) #,random_effects=c()
fit_data_hf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf,
                       normalization  = "NONE",transform =    "NONE", #min_prevalence = 0,
                       min_prevalence = 0.00001,
                       output         = "HOMAIR_female_adj",fixed_effects  = c('HOMAIR','AGE' ,'BMI'), #'Gender',
                       reference      = c("HOMAIR")) #,random_effects=c()
#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_hm = Maaslin2(input_data     = df_input_datam,
                       min_prevalence = 0.00001,
                       input_metadata = df_input_metadatam,normalization  = "NONE",transform = "NONE", #min_prevalence = -0.000001,min_abundance=-0.000001,
                       output         = "HOMAIR_male_adj",fixed_effects  = c('HOMAIR','AGE' ,'BMI'), #'Gender',
                       reference      = c("HOMAIR")) #,random_effects=c()

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_na = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata,
                       min_prevalence = 0.00001,
                       normalization  = "NONE",transform = "NONE", #min_prevalence = -0.000001,min_abundance=-0.000001,
                       output         = "Necroinflammation_all_adj",fixed_effects  = c('Necroinflammation','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("Necroinflammation,0")) #,random_effects=c()
fit_data_naf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf,
                        normalization  = "NONE",transform = "NONE",
                        min_prevalence = 0.00001,
                        # min_prevalence = -0.000001,min_abundance=-0.000001,
                        output         = "Necroinflammation_female_adj",fixed_effects  = c('Necroinflammation','AGE' ,'BMI'), #'Gender',
                        reference      = c("Necroinflammation,0")) #,random_effects=c()
fit_data_nam = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam,
                        # min_prevalence = -0.000001,min_abundance=-0.000001,
                        min_prevalence = 0.00001,
                        normalization  = "NONE",transform = "NONE",
                        output         = "Necroinflammation_male_adj",fixed_effects  = c('Necroinflammation','AGE' ,'BMI'), #'Gender',
                        reference      = c("Necroinflammation,0")) #,random_effects=c()

# #Steatosis models:
fit_data_sa = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata,
                       normalization  = "NONE",transform ="NONE", #
                       min_prevalence = 0.05,
                       output         = "Steatosis_all",
                       fixed_effects  = c('Steatosis','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("Steatosis"))

# #,random_effects=c() ok, not I think to be evaluated, since at least subjects? https://github.com/biobakery/biobakery/wiki/maaslin2
# # https://lost-stats.github.io/
# # https://julianfaraway.github.io/brinlabook/chaglmm.html

fit_data_sf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf,
                       normalization  = "NONE",transform =    "NONE", #min_prevalence = 0,
                       min_prevalence = 0.01,
                       output         = "Steatosis_female",
                       fixed_effects  = c('Steatosis','AGE' ,'BMI'), #'Gender',
                       reference      = c("Steatosis")) #,random_effects=c() ok
fit_data_sm = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam,
                       normalization  = "NONE",transform =    "NONE", #min_prevalence = -0.000001,min_abundance=-0.000001,
                       output         = "Steatosis_male_adj",fixed_effects  = c('Steatosis','AGE' ,'BMI'), #'Gender',
                       reference      = c("Steatosis")) #,random_effects=c()

#Fibrosis Models
fit_data_fa = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata,
                       normalization  = "NONE",transform =    "NONE",#min_prevalence = 0,
                       min_prevalence = 0.01,
                       output         = "Fibrosis_all_adj",fixed_effects  = c('Fibrosis','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("Fibrosis")) #,random_effects=c()
fit_data_ff = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf,
                       normalization  = "NONE",transform =    "NONE",#min_prevalence = 0,
                       min_prevalence = 0.01,
                       output         = "Fibrosis_female_adj",fixed_effects  = c('Fibrosis','AGE' ,'BMI'), #'Gender', #need ti cgecj what is the effect...
                       reference      = c("Fibrosis")) #,random_effects=c()
fit_data_fm = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam,
                       min_prevalence = 0.01,
                       normalization  = "NONE",transform =    "NONE", #min_prevalence = -0.000001,min_abundance=-0.000001,
                       output         = "Fibrosis_male_adj",fixed_effects  = c('Fibrosis','AGE' ,'BMI'), #'Gender',
                       reference      = c("Fibrosis")) #,random_effects=c()

#HOMAIR Models
fit_data_ha = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata,
                       normalization  = "NONE",transform = "NONE", #
                       # min_prevalence = 0, #
                       min_prevalence = 0.01,
                       output         = "HOMAIR_all_adj",fixed_effects  = c('HOMAIR','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("HOMAIR")) #,random_effects=c()
fit_data_hf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf,
                       normalization  = "NONE",transform =    "NONE", #min_prevalence = 0,
                       min_prevalence = 0.01,
                       output         = "HOMAIR_female_adj",fixed_effects  = c('HOMAIR','AGE' ,'BMI'), #'Gender',
                       reference      = c("HOMAIR")) #,random_effects=c()
#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_hm = Maaslin2(input_data     = df_input_datam,
                       min_prevalence = 0.01,
                       input_metadata = df_input_metadatam,normalization  = "NONE",transform = "NONE", #min_prevalence = -0.000001,min_abundance=-0.000001,
                       output         = "HOMAIR_male_adj",fixed_effects  = c('HOMAIR','AGE' ,'BMI'), #'Gender',
                       reference      = c("HOMAIR")) #,random_effects=c()

#note the prevalences and abudances ahve to be little bit below zero to get all the variables:
fit_data_na = Maaslin2(input_data     = df_input_data, input_metadata = df_input_metadata,
                       min_prevalence = 0.01,
                       normalization  = "NONE",transform = "NONE", #min_prevalence = -0.000001,min_abundance=-0.000001,
                       output         = "Necroinflammation_all_adj",fixed_effects  = c('Necroinflammation','AGE' ,'BMI','Gender'), #'Gender',
                       reference      = c("Necroinflammation")) #,random_effects=c()
fit_data_naf = Maaslin2(input_data     = df_input_dataf, input_metadata = df_input_metadataf,
                        normalization  = "NONE",transform = "NONE",
                        min_prevalence = 0.01,
                        # min_prevalence = -0.000001,min_abundance=-0.000001,
                        output         = "Necroinflammation_female_adj",fixed_effects  = c('Necroinflammation','AGE' ,'BMI'), #'Gender',
                        reference      = c("Necroinflammation")) #,random_effects=c()
fit_data_nam = Maaslin2(input_data     = df_input_datam, input_metadata = df_input_metadatam,
                        # min_prevalence = -0.000001,min_abundance=-0.000001,
                        min_prevalence = 0.01,
                        normalization  = "NONE",transform = "NONE",
                        output         = "Necroinflammation_male_adj",fixed_effects  = c('Necroinflammation','AGE' ,'BMI'), #'Gender',
                        reference      = c("Necroinflammation")) #,random_effects=c()

#Test comparison to similar model with lm: # for (i in 1:3) {SG0m[,i]=as.numeric(SG0m[,i])}
colnames(tv_all)[5] ='Steatosis'
colnames(tv_all)[8] ='HOMAIR'
xnam <- colnames(tv_all[,c(2:5)])#paste("x", 1:25, sep="") # colnames(SG0m)[1:5]
fmla <- as.formula(paste("P4 ~ ", paste(xnam, collapse= "+")));fmla
#https://stats.stackexchange.com/questions/190763/how-to-decide-which-glm-family-to-use

poissone=lm( fmla, data=tv_all)#, family = binomial) 
ap=anova(poissone);
sp=summary(poissone);

write.csv(tv_all, 'tv_all_upd.csv')


#... sample sizes in cases... use previous functions...

#this is with first(!!), use it
#is normal (non-log scale) image (and autoscaled), loge=1 is  log10 scale image (and raw data)
# Necroinflammation  HOMA-IR Steatosis.Grade.0.To.3 Fibrosis.Stage.0.to.4

# hist(df[df[,1]==e,'result.1'],breaks=50); hist(SG0[,'P4'],breaks=50) #colnames(NAFLDo)[9:28][ps<0.05]
# https://bookdown.org/content/b472c7b3-ede5-40f0-9677-75c3704c7e5c/more-than-one-mediator.html
#This works with the autoscaled (raw if loge=1 and remove 1 in the means) data NAFLD...


#If you have case/control as >0/0 (homa.ir was converted as per 1.5), or then see the other function if you have the division
sample_size=function(NAFLD,Outcome,Group) { 
  if (Group=='Male') {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==2,]} else if (Group=='Female') {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==1,]} else if (Group=='All') {NAFLDo=NAFLD} 
  sample_data=c();n0=c();n1=c()
  for (i in 1:2) {
    if (i==1) 
    {SG0=NAFLDo[NAFLDo[,Outcome] == 0,];n0=dim(SG0)[1]; sample_data=append(sample_data,n0)} else if (i==2) {SG0=NAFLDo[NAFLDo[,Outcome] > 0,];n1=dim(SG0)[1];sample_data=append(sample_data,n1)}
  }
  return(sample_data)} 

sample_size2=function(NAFLD,Outcome,Group,sick_groupe) { 
  sample_data=c();n0=c();n1=c();NAFLDo=data.frame()
  for (i in 1:2) {
    if (i==1) {if (Group=='Male') {NAFLDo=NAFLD[!sick_groupe & NAFLD[,'SEX.1F.2M']==2,]} else if (Group=='Female') 
        {NAFLDo=NAFLD[!sick_groupe & NAFLD[,'SEX.1F.2M']==1,]} else if (Group=='All') {NAFLDo=NAFLD[!sick_groupe,]}; n0=dim(NAFLDo)[1];sample_data=append(sample_data,n0)} 
    else if (i==2) {if (Group=='Male') {NAFLDo=NAFLD[sick_groupe & NAFLD[,'SEX.1F.2M']==2,]} else if (Group=='Female') 
          {NAFLDo=NAFLD[sick_groupe & NAFLD[,'SEX.1F.2M']==1,]} else if (Group=='All') {NAFLDo=NAFLD[sick_groupe,]}; n1=dim(NAFLDo)[1];sample_data=append(sample_data,n1)}}
  return(sample_data)} 


  #   #%%
  #   means=c();for (j in 9:28) {means=append(means,median(SG0[,j], na.rm=TRUE))} #
  #   sds=c();for (j in 9:28) {sds=append(sds,sd(SG0[,j],na.rm=TRUE))} #tähän jäätiin...
  #   error_lower=means-sds 
  #   # error_lower[error_lower <= 0] = 0
  #   error_upper=means+sds; error=sds
  #   sample_data <- append(sample_data, data.frame(study=colnames(NAFLD[,9:28]),index=colnames(NAFLD[,9:28]),result=means,error=error))} 
  # df=data.frame(sample_data) #
  # ps=c();for (j in 9:28) {xnam <- colnames(NAFLDo)[j]; fmla <- as.formula(paste(xnam, "~",Outcome));
  # ps=append(ps,wilcox.test(fmla, data = NAFLDo,exact = FALSE)$p.value)}#kruskal.test
  # print(mean(df[df[,1]==e,'result.1']));#hist(df[df[,1]==e,'result.1'],breaks=50)
  # print(mean(df[df[,1]==e,'result']));#hist(df[df[,1]==e,'result'],breaks=50)
  # a=df[df[,1]==e,'result.1']/df[df[,1]==e,'result']; 
  # print(a); print(log(a)) }

 # do this to get rid of Xs and other not recognizable for the function, in the case of reduced number of values, use ie instead of tv for 'NAFLD':
NAFLD=cbind(tv[,1:28])#,tv_half_log2[,1:20])#tv[,1:28]); #ie# tv_half_log22
colnames(NAFLD) #autoscaled; # NAFLD=cbind(tv[,1:8],tv_half_log2[,1:20]) #for raw # NAFLD=cbind(tv[,1:28])
NAFLD[NAFLD[,c(5)]>0,5]=1;NAFLD[NAFLD[,c(6)]>0,6]=1;NAFLD[NAFLD[,c(7)]>0,7]=1;
NAFLD[NAFLD[,c(8)] <= 1.5,8]=0;NAFLD[NAFLD[,c(8)]>1.5,8]=1; #note the order of calculation...median
colnames(NAFLD) <- gsub("-", ".", colnames(NAFLD))
colnames(NAFLD) <- gsub("/", ".", colnames(NAFLD))
colnames(NAFLD) <- gsub("11", "X11", colnames(NAFLD))
colnames(NAFLD) <- gsub("17", "X17", colnames(NAFLD))
colnames(NAFLD) <- gsub("#", ".", colnames(NAFLD))
colnames(NAFLD)[colnames(NAFLD)=='X17aOH.P4']='X17.aOHP4' #oh...

jappend=c()
date='26224';  #xlim=c(0.1,1.7)
Outcome='Steatosis.Grade.0.To.3';
Group='All'; jappend=c(jappend,sample_size(NAFLD,Outcome,Group)); 
Group='Female';jappend=c(jappend,sample_size(NAFLD,Outcome,Group))
Group='Male'; jappend=c(jappend,sample_size(NAFLD,Outcome,Group))

Outcome='Fibrosis.Stage.0.to.4'; 
Group='All'; jappend=c(jappend,sample_size(NAFLD,Outcome,Group)) #
Group='Female';jappend=c(jappend,sample_size(NAFLD,Outcome,Group))
Group='Male'; jappend=c(jappend,sample_size(NAFLD,Outcome,Group))

Outcome='Necroinflammation'; Group='All'; 
jappend=c(jappend,sample_size(NAFLD,Outcome,Group)) #not the very first though...
Group='Female';jappend=c(jappend,sample_size(NAFLD,Outcome,Group))
Group='Male'; jappend=c(jappend,sample_size(NAFLD,Outcome,Group))

Outcome='HOMA.IR';# xlim=c(0.1,5.3)
Group='All';jappend=c(jappend,sample_size(NAFLD,Outcome,Group)) #not the very first though...# xlim=c(0.1,1.7)
Group='Female';jappend=c(jappend,sample_size(NAFLD,Outcome,Group))
Group='Male'; jappend=c(jappend,sample_size(NAFLD,Outcome,Group))

matrix(unlist(jappend),ncol=8)

# For menopause
# ccovae=tv[,c("AGE")]; c1=ccovae<56; c2=ccovae>44; sick_groupe=c1 & c2 
# #tv[,c("AGE")][c1 & c2] # fn=file_names[5]; # the_essentials(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_groupe,fn,lkm,date,joo)
# Group='All';jappend=c(jappend,sample_size2(NAFLD,Outcome,Group,sick_groupe)) #not the very first though...# xlim=c(0.1,1.7)
# Group='Female';jappend=c(jappend,sample_size2(NAFLD,Outcome,Group,sick_groupe))
# Group='Male'; jappend=c(jappend,sample_size2(NAFLD,Outcome,Group,sick_groupe))
# 
# # #For masld
# ccova=tv[,c("Steatosis.Grade.0.To.3" , "Fibrosis.Stage.0.to.4" ,"Necroinflammation" ,  "HOMA-IR")] 
# fn='MASLD';sick_groupe=rowSums(ccova)>4 #toth#
# Group='All';jappend=c(jappend,sample_size2(NAFLD,Outcome,Group,sick_groupe)) #not the very first though...# xlim=c(0.1,1.7)
# Group='Female';jappend=c(jappend,sample_size2(NAFLD,Outcome,Group,sick_groupe))
# Group='Male'; jappend=c(jappend,sample_size2(NAFLD,Outcome,Group,sick_groupe))


#...
cM=c()
csd=c()






Group='Male'
if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
tv_red=c(); 
tv_red=tv[cond,]
#Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
M <- tv_red[,Mediator]  #
cM=append(cM,colMeans(M,na.rm = TRUE))
csd=append(csd,apply(M, 2, sd,na.rm = TRUE))

# cM2=cM
# csd2=csd
tot=cbind(cM,csd)
tot2=cbind(tot[1:20,],tot[21:40,],tot[41:60,])

tot3=round(tot2,0)

write.xlsx(tot3, file ='okestd3.xlsx', append = FALSE, row.names = TRUE)  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r




#...
cM=c()
csd=c()
tot=c()
tot2=c()
tot3=c()

Group='Male'
if (Group=='Female') {cond=tv_all[,'Gender']==min(tv_all[,'Gender'])} else if (Group=='Male') {cond=tv_all[,'Gender']==max(tv_all[,'Gender'])} else if (Group=='All') {cond=rep(TRUE,dim(tv_all)[1])}
tv_red=c(); 
tv_red=tv[cond,]
#Standard values did not five erros # hep=colnames(X)[!colnames(X) %in% c( "Benzylparaben" ,"Methylparaben")] # X=X[,hep]
M <- tv_red[,Mediator]  #
cM=round(apply(M, 2, median,na.rm = TRUE),0)
quants <- c(0.25,0.75)
csd=round(apply( M , 2 , quantile , probs = quants , na.rm = TRUE ),0)
tot=rbind(cM,csd)

tot2=cbind(tot2,tot)

# t4=matrix(tot2,nrow = 100, ncol = 3)

# round(apply( M , 2 , quantile , probs = quants , na.rm = TRUE ))


# cM2=cM
# csd2=csd
# tot=cbind(cM,csd)
tot3=rbind(tot2[,1:20],tot2[,21:40],tot2[,41:60])

tot4=t(round(tot3,0))

write.xlsx(tot4, file ='okestd5.xlsx', append = FALSE, row.names = TRUE)  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r



# csd
# cM


rownames(res)=rn #write.csv(rn,'iii.csv')
colnames(res)=c('ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
res=res[order(res[,2]),] #res=res[rev(order(res[,1])),]
rownames(res) <- gsub("X11", "11", rownames(res))
rownames(res) <- gsub("X17", "17", rownames(res))
write.xlsx(res, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE)  #https://stackoverflow.com/questions/21847830/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r




hyy=matrix(data = jappend, nrow = 18, ncol = 2, byrow = TRUE, dimnames = NULL)

#Other

CSVfile <- file.path('\\\\shared', 'data', 'abc.csv')
read.csv(CSVfile, header=T)`

list.files(path="('\\\\utu.fi\\taltio\\Bioscience-Sysmed\\000 Project archive\\Orion\\")


hupse=list.files(path="\\\\utu.fi\\taltio\\Bioscience-Sysmed\\000 Project archive\\Orion\\", pattern=".", all.files=TRUE, 
           full.names=FALSE)

write.csv(hupse, 'the_names_orion_tikka27224.csv')


# CSVfile <- file.path('\\\\utu.fi\\taltio\\Bioscience-Sysmed\\000 Project archive\\Orion', '\\ORION-10102-B.abf')
# asdf=read.csv(CSVfile, header=T)`

#And now the 'figure 5', forest plot info: 
#https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/forest.html
#https://rpubs.com/mbounthavong/forest_plots_r
#https://www.khstats.com/blog/forest-plots/
#https://www.geeksforgeeks.org/how-to-create-a-forest-plot-in-r/
#multiforest: https://rdrr.io/rforge/CALIBERdatamanage/man/multiforest.html
# https://www.tutorialspoint.com/how-to-change-the-linetype-for-geom-vline-in-r
#https://samuel-marsh.github.io/scCustomize/articles/Gene_Expression_Plotting.html # Alternative way to add the scale could be  to copy paste it in to image

#

# Houdaa
# Pitää tehdä lineaariset mallit semiuudestaan...
# The modelling, retesting:
# SG0=tv[tv[,'Steatosis.Grade.0.To.3']==0,c(5,9:28)]
# name="Forest plot of All Steroid Ratios in HOMA-IR"
# tv_alls=tv_all


lin_model_adjusted=function(NAFLD,Outcome,Group,name,ordera) {
Group='female';
if (Group=='male') {NAFLDo=tv_all[tv_all[,'Gender']==max(tv_all[,'Gender']),]} else if (Group=='female') 
{NAFLDo=tv_all[tv_all[,'Gender']==min(tv_all[,'Gender']),]} else if (Group=='All') {NAFLDo=tv_all}
SG0=NAFLDo[,c(2:dim(tv_all)[2])]
#tv[,'Steatosis.Grade.0.To.3']==0# colnames(SG0[,2:21]) <- gsub(" ", ".", colnames(SG0[,2:21])) 
#https://stackoverflow.com/questions/10688137/how-to-fix-spaces-in-column-names-of-a-data-frame-remove-spaces-inject-dots
oknames=colnames(SG0)
SG0=data.frame(SG0)
colnames(SG0)
colnames(SG0[,8:27]) <- gsub("-", ".", colnames(SG0[,8:27]))
colnames(SG0[,8:27]) <- gsub("/", ".", colnames(SG0[,8:27]))
# SG0=destroyX(data.frame(SG0))
hesh=c()
xnam <- colnames(SG0)[c(4:7)]
xnam2 <- colnames(SG0)[c(70:76)]
j=1;i=1
for (i in 1:length(xnam)) {
  for (j in 1:length(xnam2)) {
    if (Group!='All')  {fmla <- as.formula(paste(paste(c(xnam2[j]," ~ "), collapse= ""), paste(c(xnam[i],'BMI','AGE'), collapse= "+")))} else if (Group=='All') 
    {fmla <- as.formula(paste(paste(c(xnam2[j]," ~ "), collapse= ""), paste(c(xnam[i],'BMI','AGE','Gender'), collapse= "+")))} #https://stats.stackexchange.com/questions/190763/how-to-decide-which-glm-family-to-use
    poissone=lm( fmla, data=SG0) # anova(poissone);# poissone
    ps=summary(poissone);
    pss=ps[[4]] # fmla <- as.formula(paste(paste(c(colnames(SG0[,8:27])[j]," ~ "), collapse= ""), paste(c(xnam[i],'AGE','BMI'), collapse= "+")))
    hesh=rbind(hesh,c(xnam2[j],xnam[i],Group,pss[2,1],pss[2,4]))}}
j=1;i=1
xnam <- colnames(SG0)[c(2)]
xnam2 <- colnames(SG0)[c(70:76)]
for (i in 1:length(xnam)) {
  for (j in 1:length(xnam2)) {
    if (Group!='All')  {fmla <- as.formula(paste(paste(c(xnam2[j]," ~ "), collapse= ""), paste(c(xnam[i],'BMI'), collapse= "+")))} else if (Group=='All') 
    {fmla <- as.formula(paste(paste(c(xnam2[j]," ~ "), collapse= ""), paste(c(xnam[i],'BMI','Gender'), collapse= "+")))} #https://stats.stackexchange.com/questions/190763/how-to-decide-which-glm-family-to-use
    poissone=lm( fmla, data=SG0) # anova(poissone);# poissone
    ps=summary(poissone);
    pss=ps[[4]] 
    hesh=rbind(hesh,c(xnam2[j],xnam[i],Group,pss[2,1],pss[2,4]))}}
j=1;i=1
xnam <- colnames(SG0)[c(3)]
for (i in 1:length(xnam)) {
  for (j in 1:length(xnam2)) {
    if (Group!='All')  {fmla <- as.formula(paste(paste(c(xnam2[j]," ~ "), collapse= ""), paste(c(xnam[i],'AGE'), collapse= "+")))} else if (Group=='All') 
    {fmla <- as.formula(paste(paste(c(xnam2[j]," ~ "), collapse= ""), paste(c(xnam[i],'AGE','Gender'), collapse= "+")))} #https://stats.stackexchange.com/questions/190763/how-to-decide-which-glm-family-to-use
    poissone=lm( fmla, data=SG0) # anova(poissone);# poissone
    ps=summary(poissone);
    pss=ps[[4]] # fmla <- as.formula(paste(paste(c(colnames(SG0[,8:27])[j]," ~ "), collapse= ""), paste(c(xnam[i],'AGE','BMI'), collapse= "+")))
    hesh=rbind(hesh,c(xnam2[j],xnam[i],Group,pss[2,1],pss[2,4]))}}

# hoi=c();for (i in 1:dim(rt2)[1]) {hoi=append(hoi,scan(text=rownames(rt2)[i], what=""))}
j=1;i=1
xnam <- colnames(SG0)[c(70:76)]#paste("x", 1:25, sep="") 28:length(colnames(SG0)
for (i in 1:length(xnam)) {
  for (j in 1:length(colnames(SG0[,8:27]))) {
    if (Group!='All')  {fmla <- as.formula(paste(paste(c(colnames(SG0[,8:27])[j]," ~ "), collapse= ""), paste(c(xnam[i],'AGE','BMI'), collapse= "+")))} else if (Group=='All') 
    {fmla <- as.formula(paste(paste(c(colnames(SG0[,8:27])[j]," ~ "), collapse= ""), paste(c(xnam[i],'AGE','BMI','Gender'), collapse= "+")))} #https://stats.stackexchange.com/questions/190763/how-to-decide-which-glm-family-to-use
    poissone=lm( fmla, data=SG0) # anova(poissone);# poissone
    ps=summary(poissone);
    pss=ps[[4]] # fmla <- as.formula(paste(paste(c(colnames(SG0[,8:27])[j]," ~ "), collapse= ""), paste(c(xnam[i],'AGE','BMI'), collapse= "+")))
    hesh=rbind(hesh,c(colnames(SG0[,8:27])[j],xnam[i],Group,pss[2,1],pss[2,4]))}}

hesa=hesh
hoi=as.data.frame(hesh)
rsa=c()
for (i in 4:5) {
rs=hoi[,c(i)] # rs=data.frame(rs)
rs=matrix(as.numeric(rs),nrow=7)
# hoi=as.data.frame(matrix(as.numeric(rs), ncol = 20,  byrow = TRUE), stringsAsFactors = TRUE)
cn=c(colnames(SG0)[c(4:7)],'AGE','BMI',colnames(SG0)[c(8:27)])
rn=colnames(SG0)[c(70:76)]; rs=as.data.frame(rs)
rownames(rs)=rn; colnames(rs)=cn
colnames(rs) <- gsub("\\.", "-", colnames(rs))
colnames(rs) <- gsub("X11", "11", colnames(rs))
colnames(rs) <- gsub("X17", "17", colnames(rs))
colnames(rs)["Steatosis-Grade"]
colnames(rs)[colnames(rs)=="T-Epi-T"]="T/Epi-T"
colnames(rs)[colnames(rs)=="Steatosis-Grade"]="Steatosis Grade"
colnames(rs)[colnames(rs)=="Fibrosis-Stage"]="Fibrosis Stage"
colnames(rs)[colnames(rs)=="17aOH-P4"]="17a-OHP4"
rsa=rbind(rsa,rs) }
rs1a=rsa[1:7,];
rs2a=rsa[8:14,]
rs1=rs1a;rs2=rs2a
rownames(rs2)=str_sub(rownames(rs2), end = -2)
# https://scales.arabpsychology.com/stats/how-to-remove-the-last-character-from-a-string-in-r-2-examples/
# rs<0.05
rs1=as.data.frame(rs1[,colnames(resulta1)]); 
rs2=as.data.frame(rs2[,colnames(resulta1)])
# colnames(resulta3)[!colnames(resulta3) %in% colnames(rs1)]
tv_ah=rango(rs1,-0.5,0.5); rs1=tv_ah;#
order="original"; range='orig';corre='no_renorm'; type='full'; method='color';ga='All';gf='Female';gm='Male' #color square
cl.offset=1.0;cl.length=5;cl.cex = 1.3;pch.cex=1.3;pch=20;cl.pos = 'r';#cl.pos = 'b' ;#pch.cex=0.95,1.3; height=6300; pos 'b' cl.pos = 'b' 
ho=Group;hip1='hurraa'
rs1=as.matrix(rs1)
rs2=as.matrix(rs2)
tv_ah=rango(rs1,-0.5,0.5); rs1=tv_ah;#
corrplot(rs1, type = type, order = order,method=method, p.mat=rs2)
width=6000; height=2800

jpeg(paste("Square Correlation Plot of_e",ho,gm,hip1,".jpg"), width = width, height = height, quality = 100,pointsize = 30, res=300); 
corrplot(rs1, type = type, order = order,method=method, p.mat=rs2, tl.col = "black", cl.cex = cl.cex,pch.cex=pch.cex,pch.col='black',pch=pch,
         sig.level = c(.001, .05, .2),cl.pos = cl.pos, insig = "label_sig",cl.offset=cl.offset,cl.length=cl.length,
         tl.srt = 90, diag = TRUE,col = rev(COL2('RdBu')[25:(length(COL2('RdBu'))-25)]),is.corr = FALSE) #,is.corr = FALSE
dev.off()
}

# jpeg(paste("Square Correlation Plot ofx3",ho,gm,hip1,".jpg"), width = width, height = height, quality = 100,pointsize = 30, res=300);
# corrplot(rs1, type = type, order = order,method=method, p.mat=rs2, tl.col = "black", cl.cex = cl.cex,pch.cex=pch.cex,
#          pch.col='black',pch=pch,
#          sig.level = c(.001, .05, .2),cl.pos = cl.pos, insig = "label_sig",cl.offset=cl.offset,cl.length=cl.length,
#          tl.srt = 90, diag = TRUE,col = rev(COL2('RdBu')[25:(length(COL2('RdBu'))-25)]),is.corr = FALSE) #,is.corr = FALSE
# dev.off()
# rs1=as.data.frame(rs1[,colnames(resulta3)])
# rs2=as.data.frame(rs2[,colnames(resulta3)])
# colnames(resulta3)[!colnames(resulta3) %in% colnames(rs1)]
# 
# tv_ah=rango(rs1,-0.5,0.5); rs1=tv_ah;#
# 
# rs1 <- cor(rs1)
# # corrplot(M, method="number")
# 
# all(sapply(rs1, is.finite))



# devtools::install_github('cran/ggplot2') 


jepsm=tvm[,c('DHEA','DHT','T/Epi-T','PFHpA')];jepsm=data.frame(jepsm)
a='T.Epi.T'; b='PFHpA'  
r <-  round(cor(jeps[,c('PFHpA')], jeps[,a]), 3);       p  <- cor.test(jeps[,c('PFHpA')], jeps[,a])$p.value
rf <- round(cor(jepsf[,c('PFHpA')], jepsf[,c(a)]), 3); pf <- cor.test(jepsf[,c('PFHpA')], jepsf[,c(a)])$p.value
rm <- round(cor(jepsm[,c('PFHpA')], jepsm[,c(a)]), 3); pm <- cor.test(jepsm[,c('PFHpA')], jepsm[,c(a)])$p.value

jpeg(paste("Correlations plotted_all",a,b,".jpg"), width = 1000, height = 1000, quality = 100,pointsize = 20, res=300); 
ggplot(jeps, aes(y=T.Epi.T, x=PFHpA)) + 
  geom_point() +
  geom_smooth(method="lm", col="black") +
  annotate("text", x=min(jeps[,c('PFHpA')]), y=max(jeps[,c(a)]), label=paste0("r = ", r), hjust=0) +
  annotate("text", x=min(jeps[,c('PFHpA')])+0.07, y=(max(jeps[,c(a)])-1.3), label=paste0("p = ", round(p, 5)), hjust=0) +
  theme_classic()
dev.off()



order="original"; range='orig';corre='no_renorm'; type='full'; method='color' #color square
jpeg(paste("Square Correlation Plot of PFAS vs. clinical factors and steroids_alle_0.2ee.jpg"), width = 10000, height = 7000, quality = 100,pointsize = 30, res=300); 
corrplot(tte, type = type, order = order,method=method, tl.col = "black",  cl.pos = 'r', ,cl.offset=1.5,cl.length=5,tl.srt = 90, diag = FALSE,col = c('blue','white','red'),is.corr = FALSE) #,is.corr = FALSE
# corrplot(resulta3, type = type, order = order,method=method,  tl.col = "black",  pch.cex=2,cl.pos = 'b',cl.offset=1.5,cl.length=5,tl.srt = 90, diag = FALSE,col = rev(COL2('RdBu')),is.corr = FALSE) #,is.corr = FALSE
dev.off()



# install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
# my_data <- mtcars[, c(1,3,4,5,6,7)]
# chart.Correlation(my_data, histogram=TRUE, pch=19)
# chart.Correlation(jeps, histogram=FALSE, pch=16)
# chart.Correlation(jepsf, histogram=FALSE, pch=16)
# chart.Correlation(jepsm, histogram=FALSE, pch=16)

tv_c=tv #cbind(tv[,1:8], tv[,9:dim(tv)[2]]) #check also not logged and then the auto one
# colnames(tv_c) <- gsub("-", ".", colnames(tv_c))
# colnames(tv_c) <- gsub("/", ".", colnames(tv_c))
# colnames(tv_c) <- gsub("11", "X11", colnames(tv_c))
# colnames(tv_c) <- gsub("17", "X17", colnames(tv_c))
# colnames(tv_c) <- gsub("#", ".", colnames(tv_c))
x1=colnames(tv_c[,c(1:8)]); v2=dim(NAFLD)[2]+1
x2=colnames(tv_c[,9:v2]);v3=(dim(Bali)[2]+v2)
x3=colnames(tv_c[,(v2+1):(v3)]);v4=(dim(Base)[2])+v3
x4=colnames(tv_c[,(v3+1):(v4-1)])
x5=colnames(tv_c[,(v4):(dim(tv_c)[2])]); 
x5a=x5[1:9]
x6=x5[10:length(x5)] #dividing to lipids
x5=x5a  #making sure that PFAS are separate
# x1;x2;x3;x4;x5; # tv_auto[1:5,1:11]
x3 <- paste(x3, "_L", sep="") #https://stackoverflow.com/questions/6984796/how-to-paste-a-string-on-each-element-of-a-vector-of-strings-using-apply-in-r
x4=gsub("(-[0-9]*)*.1", "", x4) #https://stackoverflow.com/questions/18997297/remove-ending-of-string-with-gsub
x4 <- paste(x4, "_S", sep="")# https://rdrr.io/bioc/qpgraph/man/qpNrr.html
nm = c(x1,x2,x3,x4,x5,x6) # adding the lipids, x6
nm=c('Subject#','Gender','AGE','BMI','Steatosis Grade','Fibrosis Stage','Necroinflammation','HOMA-IR',nm[9:93])
colnames(tv_c)=nm # tv_c[1:5,1:30]; NAFLD[1:2,1:28];
colnames(tv_c)[colnames(tv_c)=='MASLD[, \"PFAS\"]']='PFAS';
#removing PFAS and lipid generals:
x5=x5[x5!='PFAS']; x6=x6[x6!='Total_TG'] 
# x1;x2;x3;x4;x5;
# The correlations (ok):
# ccovae=tv[,c("AGE")]; c1=ccovae<56; c2=ccovae>44; sick_group=c1 & c2
# tv_c=tv#tv_half_log22 #cbind(tv[,1:8], tv_half_log2) #check also not logged and then the auto one
# tv_c[,'Menopause']=as.numeric(sick_group)
# tv_c=tv_c[,c(1:3,dim(tv_c)[2],4:(dim(tv_c)[2]-1))]
colnames(tv_c) <- gsub("-", ".", colnames(tv_c))
colnames(tv_c) <- gsub("/", ".", colnames(tv_c))
colnames(tv_c) <- gsub("11", "X11", colnames(tv_c))
colnames(tv_c) <- gsub("17", "X17", colnames(tv_c))
colnames(tv_c) <- gsub("#", ".", colnames(tv_c))
colnames(tv_c)[colnames(tv_c)=="X17aOH.P4"]="X17a.OHP4"

tv_c=tv_c[,!colnames(tv_c) %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]
tv_c=tv_c[,!colnames(tv_c) %in% x4]
colnames(tv_c)[colnames(tv_c)=="17aOH-P4"]="17a-OHP4"
tvf=tv_c[tv_c[,'Gender']==min(tv_c[,'Gender']),1:dim(tv_c)[2]] #tv['Steatosis.Grade.0.To.3'==0,9:27]] #tv[tv[,'Necroinflammation']==0,9:80]; #SG0i=as.numeric(SG0i); check also: tv[tv[,'HOMA-IR']==0,9:80]
tvm=tv_c[tv_c[,'Gender']==max(tv_c[,'Gender']),1:dim(tv_c)[2]]
# tvf=tv_c[tv_c[,"SEX.1F.2M"]==min(tv_c[,"SEX.1F.2M"]),1:dim(tv_c)[2]] #tv['Steatosis.Grade.0.To.3'==0,9:27]] #tv[tv[,'Necroinflammation']==0,9:80]; #SG0i=as.numeric(SG0i); check also: tv[tv[,'HOMA-IR']==0,9:80]
# tvm=tv_c[tv_c[,"SEX.1F.2M"]==max(tv_c[,"SEX.1F.2M"]),1:dim(tv_c)[2]]
# jeps=tv_c[,c('DHEA','DHT','T/Epi-T','PFHpA')];jeps=data.frame(jeps)
# jepsf=tvf[,c('DHEA','DHT','T/Epi-T','PFHpA')];jepsf=data.frame(jepsf)
# jepsm=tvm[,c('DHEA','DHT','T/Epi-T','PFHpA')];jepsm=data.frame(jepsm)
# colnames(tvf)[colnames(tvf)=="17aOH-P4"]="17a-OHP4"
# colnames(tvm)[colnames(tvm)=="17aOH-P4"]="17a-OHP4"
jeps  = tv_c; #jeps=data.frame(jeps)
jepsf = tvf;  #jepsf=data.frame(jepsf)
jepsm = tvm;  #jepsm=data.frame(jepsm)

# colnames(tv_c) <- gsub("-", ".", colnames(tv_c))
# colnames(tv_c) <- gsub("/", ".", colnames(tv_c))
# colnames(tv_c) <- gsub("11", "X11", colnames(tv_c))
# colnames(tv_c) <- gsub("17", "X17", colnames(tv_c))
# colnames(tv_c) <- gsub("#", ".", colnames(tv_c))
# colnames(tv_c)[colnames(tv_c)=="17aOH.P4"]="17a.OHP4"
# colnames(tv_c)[colnames(tv_c)=="17aOH-P4"]="17a-OHP4"



#...
Treatment
Mediator
Outcome

Mediator[!Mediator %in% colnames(jeps)]

# Treatment=colnames(tv_all)[71:77];
# Mediator=colnames(tv_all)[9:28];
# Outcome=colnames(tv_all)[c(29:51,78:90)]; ##https://sparkbyexamples.com/r-programming/r-remove-from-vector-with-examples/
# Treatment=Treatment[!Treatment %in% c('Perfluorodecyl.ethanoic.acid')]
# tv_all=tv_all[,!colnames(tv_all) %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]
# tv_all=tv_all[,!colnames(tv_all) %in% x4]
# tv_all=tv_all[,!colnames(tv_all) %in% c('TG_SFA','MUFA','TG_PUFA')]
# Outcome=Outcome[!Outcome %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]
# Outcome=Outcome[! Outcome %in% x4] #https://sparkbyexamples.com/r-programming/r-remove-from-vector-with-examples/
# Outcome=Outcome[! Outcome %in% c('TG_SFA','MUFA','TG_PUFA')] #

# Nyt löyty looppi menetelmä
# a='T.Epi.T'; b='PFHpA'  
# a='T.Epi.T'
ioppi=Mediator2
pp0=Mediator
Treatmenta=Treatment


#(y=x); steroid = pfas, covar =pfas, #ok pfas=covar
#(y=x=; covar= steroid, BA/lipid ~ steroid #steroid=covar  

#uudet, pitäis olla ok: 
#pfas=covar; steroid=covar; BA/lipid=covar 
#....: BA/lipid ~ steroid, steroid ~ PFAS

# 1) pfas=covar; steroid=covar; BA/lipid=covar; Mediator2=colnames(tv_all)[9:28];
# Mediator2=c(colnames(jeps)[3:8]) #colnames(tv_all)[9:28]
# Mediator=c(colnames(jeps)[3:8])
#x eli treatment
Treatment = colnames(jeps)[3:8]

# Mediator2=c(colnames(jeps)[3:8],Outcome,c("TG_SFA", "MUFA", "TG_PUFA"))
#responses
Mediator2=c(x2,x5,x3,x6)
Mediator = Mediator2

#2) BA/lipid ~ steroid, steroid ~ PFAS
#a)
Treatment = c(x2)#colnames(jeps)[3:8]
Mediator= c(x3,x6)#Mediator2[1:20]
Mediator2=c(x3,x6)#Mediator2[1:20]
#b)
Treatment = c(x5)#colnames(jeps)[3:8]
Mediator= c(x2)#Mediator2[1:20]
Mediator2=c(x2)#Mediator2[1:20]

Mediator=c(colnames(jeps)[3:8],Outcome,c("TG_SFA", "MUFA", "TG_PUFA"))





Mediator <- gsub("-", ".", Mediator)
Mediator <- gsub("/", ".", Mediator)
Mediator <- gsub("11", "X11", Mediator)
Mediator <- gsub("17", "X17", Mediator)
Mediator <- gsub("#", ".", Mediator)
Mediator[Mediator=="X17aOH.P4"]="X17a.OHP4"
# Mediator[Mediator=="17aOH-P4"]="17a.OHP4"

Mediator2 <- gsub("-", ".", Mediator)
Mediator2 <- gsub("/", ".", Mediator)
Mediator2 <- gsub("11", "X11", Mediator)
Mediator2 <- gsub("17", "X17", Mediator)
Mediator2 <- gsub("#", ".", Mediator)
Mediator2[Mediator2=="X17aOH.P4"]="X17a.OHP4"

# Treatment

Treatment <- gsub("-", ".", Treatment)
Treatment <- gsub("/", ".", Treatment)
Treatment <- gsub("11", "X11", Treatment)
Treatment <- gsub("17", "X17", Treatment)
Treatment <- gsub("#", ".", Treatment)
Treatment[Treatment=="X17aOH.P4"]="X17a.OHP4"



Treatment = c(x2)#colnames(jeps)[3:8]

eh=colnames(tv_all)[2:8]
Mediator=c(eh,x3,x5,x6)
Mediator = Mediator2




# pcor(x, method = c("pearson"))

tryCatch(pm  <- round(pcor(cbind(jepsm[,Treatment[j]], jepsm[,Mediator[i]],jepsm[,"AGE"],jepsm[,"BMI"]),method="spearman")$p.value[2,1],5),
                 error = function(cond) {pm <- cor.test(jepsm[,Treatment[j]], jepsm[,Mediator[i]],method="spearman")$p.value})

pm = tryCatch({round(pcor(cbind(jepsm[,Treatment[j]], jepsm[,Mediator[i]],jepsm[,"AGE"],jepsm[,"BMI"]),method="spearman")$p.value[2,1],5)}, warning = function(w) {round(cor.test(jepsm[,Treatment[j]], jepsm[,Mediator[i]],method="spearman")$p.value,5)})

# fs::dir_create( path = sub_dir)

# unlink("C:/Users/patati/Desktop/TurkuOW/RWork/males/",recursive = TRUE, force = TRUE)
ld=list.dirs(path = "C:/Users/patati/Desktop/TurkuOW/RWork/males/", full.names = TRUE, recursive = TRUE)
for (i in 1:length(ld)) {unlink(paste0(ld[i],"/*",sep=""))} #https://www.geeksforgeeks.org/obtain-list-of-directories-in-r/
# do.call(file.remove, list(list.files("C:/Users/patati/Desktop/TurkuOW/RWork/males/*", full.names = TRUE)))
# lda=list.dirs(path = "C:/Users/patati/Desktop/TurkuOW/RWork/alla/", full.names = TRUE, recursive = TRUE)
# ldf=list.dirs(path = "C:/Users/patati/Desktop/TurkuOW/RWork/females/", full.names = TRUE, recursive = TRUE)
# ldm=list.dirs(path = "C:/Users/patati/Desktop/TurkuOW/RWork/males/", full.names = TRUE, recursive = TRUE)
# sub_dir<-Treatment
# file.exists(paste("C:/Users/patati/Desktop/TurkuOW/RWork/alla/",sub_dir[i],sep=""))
# setting up the main directory
main_dir <- "C:\\Users\\patati\\Desktop\\TurkuOW\\RWork\\males\\"
# setting up the sub directory

for (i in 1:length(Treatment)) {sub_dir <- Treatment[i]
# check if sub directory exists 
if (file.exists(sub_dir)){# specifying the working directory
  setwd(file.path(main_dir, sub_dir))} else {
  # create a new sub directory inside# the main path
  dir.create(file.path(main_dir, sub_dir))
  # specifying the working directory
  setwd(file.path(main_dir, sub_dir))}}


# jeps=...
                
i=1;j=1


for (i in 1:length(Mediator)) {
  for (j in 1:length(Treatment)) {
    r = tryCatch({c(round(pcor(cbind(jeps[,Treatment[j]], jeps[,Mediator[i]],jeps[,"AGE"],jeps[,"BMI"],jeps[,"Gender"]),method="spearman")$estimate[2,1],5),f='')}, 
                 error = function(w) {c(round(cor.test(jeps[,Treatment[j]], jeps[,Mediator[i]],method="spearman")$estimate,5),f='no adj.')}, 
                 warning = function(w) {c(round(cor.test(jeps[,Treatment[j]], jeps[,Mediator[i]],method="spearman")$estimate,5),f='no adj.')}); 
    r=as.numeric(r[1])
    
    
    p = tryCatch({c(round(pcor(cbind(jeps[,Treatment[j]], jeps[,Mediator[i]],jeps[,"AGE"],jeps[,"BMI"],jeps[,"Gender"]),method="spearman")$p.value[2,1],5),f='')}, 
  error = function(w) {c(round(cor.test(jeps[,Treatment[j]], jeps[,Mediator[i]],method="spearman")$p.value,5),f='no adj.')}, 
  warning = function(w) {c(round(cor.test(jeps[,Treatment[j]], jeps[,Mediator[i]],method="spearman")$p.value,5),f='no adj.')}); u=p[2];p=as.numeric(p[1]);#f='no adj.'
    
    
    rf <- tryCatch({c(round(pcor(cbind(jepsf[,Treatment[j]], jepsf[,Mediator[i]],jepsf[,"AGE"],jepsf[,"BMI"]),method="spearman")$estimate[2,1],5),f='')}, 
                   error = function(w) {c(round(cor.test(jepsf[,Treatment[j]], jepsf[,Mediator[i]],method="spearman")$estimate,5),f='no adj.')}, 
                   warning = function(w) {c(round(cor.test(jepsf[,Treatment[j]], jepsf[,Mediator[i]],method="spearman")$estimate,5),f='no adj.')}); 
    rf=as.numeric(rf[1]) 
    
    pf = tryCatch({c(round(pcor(cbind(jepsf[,Treatment[j]], jepsf[,Mediator[i]],jepsf[,"AGE"],jepsf[,"BMI"]),method="spearman")$p.value[2,1],5),f='')}, 
                  error = function(w) {c(round(cor.test(jepsf[,Treatment[j]], jepsf[,Mediator[i]],method="spearman")$p.value,5),f='no adj.')},
                  warning = function(w) {c(round(cor.test(jepsf[,Treatment[j]], jepsf[,Mediator[i]],method="spearman")$p.value,5),f='no adj.')}); 
    f=pf[2];pf=as.numeric(pf[1]);#f='no adj.'
    
    rm <- tryCatch({c(round(pcor(cbind(jepsm[,Treatment[j]], jepsm[,Mediator[i]],jepsm[,"AGE"],jepsm[,"BMI"]),method="spearman")$estimate[2,1],5),f='')}, 
                   error = function(w) {c(round(cor.test(jepsm[,Treatment[j]], jepsm[,Mediator[i]],method="spearman")$estimate,5),f='no adj.')}, 
                   warning = function(w) {c(round(cor.test(jepsm[,Treatment[j]], jepsm[,Mediator[i]],method="spearman")$estimate,5),f='no adj.')}); 
    rm=as.numeric(rm[1]) 
    
    pm = tryCatch({c(round(pcor(cbind(jepsm[,Treatment[j]], 
                                      jepsm[,Mediator[i]],jepsm[,"AGE"],
                                      jepsm[,"BMI"]),method="spearman")$p.value[2,1],5),f='')}, 
                  error = function(w) {c(round(cor.test(jepsm[,Treatment[j]], jepsm[,Mediator[i]],method="spearman")$p.value,5),f='no adj.')},
                  warning = function(w) {c(round(cor.test(jepsm[,Treatment[j]], jepsm[,Mediator[i]],method="spearman")$p.value,5),f='no adj.')}); 
    m=pm[2];pm=as.numeric(pm[1]); #f='no adj.'

    setwd(paste("C:/Users/patati/Desktop/TurkuOW/RWork/alla/",Treatment[j],sep=""))
    jpeg(paste("Correlations plotted_all",Treatment[j],Mediator2[i],".jpg"), width = 1000, height = 1000, quality = 100,pointsize = 20, res=300); 
    a=ggplot(jeps, aes(y=jeps[,Mediator[i]], x=jeps[,Treatment[j]])) + 
      geom_point() +
      xlab(Treatment[j]) +
      ylab(Mediator2[i]) +
      geom_smooth(method="lm", col="black") +
      annotate("text", x=min(jeps[,Treatment[j]]), y=jeps[rev(order(jeps[, Mediator[i]])),Mediator[i]][2], label=paste0("r = ", r), hjust=0) +
      annotate("text", x=min(jeps[,Treatment[j]])+0.02, y=jeps[rev(order(jeps[, Mediator[i]])),Mediator[i]][2]-sd(jeps[, Mediator[i]])/2.3, label=paste0("p = ", round(p, 5)," ",u), hjust=0) +
      theme_classic()
    print(a)
    dev.off()
    setwd(paste("C:/Users/patati/Desktop/TurkuOW/RWork/females/",Treatment[j],sep=""))
    jpeg(paste("Correlations plotted_female",Treatment[j],Mediator2[i],".jpg"), width = 1000, height = 1000, quality = 100,pointsize = 20, res=300); 
    b=ggplot(jepsf, aes(y=jepsf[,Mediator[i]], x=jepsf[,Treatment[j]])) + 
      geom_point() +
      xlab(Treatment[j]) +
      ylab(Mediator2[i]) +
      geom_smooth(method="lm", col="black") +
      annotate("text", x=min(jepsf[,Treatment[j]]), y=jepsf[rev(order(jepsf[, Mediator[i]])),Mediator[i]][2], label=paste0("r = ", rf), hjust=0) +
      annotate("text", x=min(jepsf[,Treatment[j]])+0.02, y=jepsf[rev(order(jepsf[, Mediator[i]])),Mediator[i]][2]-sd(jepsf[, Mediator[i]])/2.3, label=paste0("p = ", round(pf, 5)," ",f), hjust=0) +
      theme_classic()
      print(b)
    dev.off()
    
    setwd(paste("C:/Users/patati/Desktop/TurkuOW/RWork/males/",Treatment[j],sep=""))
    jpeg(paste("Correlations plotted_male",Treatment[j],Mediator2[i],".jpg"), width = 1000, height = 1000, quality = 100,pointsize = 20, res=300); 
    c=ggplot(jepsm, aes(y=jepsm[,Mediator[i]], x=jepsm[,Treatment[j]])) + 
      geom_point() +
      xlab(Treatment[j]) +
      ylab(Mediator2[i]) +
      geom_smooth(method="lm", col="black") +
      annotate("text", x=min(jepsm[,Treatment[j]]), y=jepsm[rev(order(jepsm[, Mediator[i]])),Mediator[i]][2], label=paste0("r = ", rm), hjust=0) +
      annotate("text", x=min(jepsm[,Treatment[j]])+0.02, y=jepsm[rev(order(jepsm[, Mediator[i]])),Mediator[i]][2]-sd(jepsm[, Mediator[i]])/2.3, label=paste0("p = ", round(pm, 5)," ",m), hjust=0) +
      theme_classic()
      print(c)
    dev.off()
    
  }}





i=1;j=1

for (i in 1:length(Mediator)) {
  for (j in 1:length(Treatment)) {
    r = c(round(cor.test(jeps[,Treatment[j]], jeps[,Mediator[i]],method="spearman")$estimate,5),f='no adj.')
    r=as.numeric(r[1])
    p = c(round(cor.test(jeps[,Treatment[j]], jeps[,Mediator[i]],method="spearman")$p.value,5),f='no adj.')
    u=p[2];p=as.numeric(p[1]);#f='no adj.'
    
    
    rf <- c(round(cor.test(jepsf[,Treatment[j]], jepsf[,Mediator[i]],method="spearman")$estimate,5),f='no adj.')
    rf=as.numeric(rf[1]) 
    pf = c(round(cor.test(jepsf[,Treatment[j]], jepsf[,Mediator[i]],method="spearman")$p.value,5),f='no adj.')
    f=pf[2];pf=as.numeric(pf[1]); #f='no adj.'
    
    rm <- c(round(cor.test(jepsm[,Treatment[j]], jepsm[,Mediator[i]],method="spearman")$estimate,5),f='no adj.')
    rm=as.numeric(rm[1]) 
    pm = c(round(cor.test(jepsm[,Treatment[j]], jepsm[,Mediator[i]],method="spearman")$p.value,5),f='no adj.') 
    m=pm[2];pm=as.numeric(pm[1]); #f='no adj.'
    
    setwd(paste("C:/Users/patati/Desktop/TurkuOW/RWork/alla/",Treatment[j],sep=""))
    jpeg(paste("Correlations plotted_all",Treatment[j],Mediator2[i],".jpg"), width = 1000, height = 1000, quality = 100,pointsize = 20, res=300); 
    a=ggplot(jeps, aes(y=jeps[,Mediator[i]], x=jeps[,Treatment[j]])) + 
      geom_point() +
      xlab(Treatment[j]) +
      ylab(Mediator2[i]) +
      geom_smooth(method="lm", col="black") +
      annotate("text", x=-Inf,y=-Inf,hjust=0,vjust=0, label=paste0("r = ", r)) +
      annotate("text", x=-Inf,y=-Inf,hjust=0,vjust=0, label=paste0("p = ", round(p, 5)," ",u)) +
      theme_classic()
    print(a)
    dev.off()
    setwd(paste("C:/Users/patati/Desktop/TurkuOW/RWork/females/",Treatment[j],sep=""))
    jpeg(paste("Correlations plotted_female",Treatment[j],Mediator2[i],".jpg"), width = 1000, height = 1000, quality = 100,pointsize = 20, res=300); 
    b=ggplot(jepsf, aes(y=jepsf[,Mediator[i]], x=jepsf[,Treatment[j]])) + 
      geom_point() +
      xlab(Treatment[j]) +
      ylab(Mediator2[i]) +
      geom_smooth(method="lm", col="black") +
      annotate("text", x=-Inf,y=-Inf,hjust=0,vjust=0, label=paste0("r = ", rf)) +
      annotate("text", x=-Inf,y=-Inf,hjust=0,vjust=0, label=paste0("p = ", round(pf, 5)," ",f)) +
      theme_classic()
    print(b)
    dev.off()
    
    setwd(paste("C:/Users/patati/Desktop/TurkuOW/RWork/males/",Treatment[j],sep=""))
    jpeg(paste("Correlations plotted_male",Treatment[j],Mediator2[i],".jpg"), width = 1000, height = 1000, quality = 100,pointsize = 20, res=300); 
    c=ggplot(jepsm, aes(y=jepsm[,Mediator[i]], x=jepsm[,Treatment[j]])) + 
      geom_point() +
      xlab(Treatment[j]) +
      ylab(Mediator2[i]) +
      geom_smooth(method="lm", col="black") +
      annotate("text", x=-Inf,y=-Inf,hjust=0,vjust=0, label=paste0("r = ", rm)) +
      annotate("text", x=-Inf,y=-Inf,hjust=0,vjust=0, label=paste0("p = ", round(pm, 5)," ",m)) +
      theme_classic()
    print(c)
    dev.off()
    
  }}

# ggplot(jeps, aes(y=T.Epi.T, x=PFHpA))

# jpeg(paste("Correlations plotted_female",a,b,".jpg"), width = 1000, height = 1000, quality = 100,pointsize = 20, res=300); 
# ggplot(jepsf, aes(y=T.Epi.T, x=PFHpA)) + 
#   geom_point() +
#   geom_smooth(method="lm", col="black") +
#   annotate("text", x=min(jepsf[,c('PFHpA')]), y=max(jepsf[,c(a)]), label=paste0("r = ", rf), hjust=0) +
#   annotate("text", x=min(jepsf[,c('PFHpA')])+0.07, y=(max(jepsf[,c(a)])-1.3), label=paste0("p = ", round(pf, 5)), hjust=0) +
#   theme_classic()
# dev.off()
# 
# jpeg(paste("Correlations plotted_male",a,b,".jpg"), width = 1000, height = 1000, quality = 100,pointsize = 20, res=300);
# ggplot(jepsm, aes(y=T.Epi.T, x=PFHpA)) +
#   geom_point() +
#   geom_smooth(method="lm", col="black") +
#   annotate("text", x=min(jepsm[,c('PFHpA')]), y=max(jepsm[,c(a)]), label=paste0("r = ", rm), hjust=0) +
#   annotate("text", x=min(jepsm[,c('PFHpA')])+0.07, y=(max(jepsm[,c(a)])-1.3), label=paste0("p = ", round(pm, 5)), hjust=0) +
#   theme_classic()
# dev.off()

# DHEA

# jeps1=jeps[,1:12]
# jeps2=jeps[,13:24]
# jeps3=jeps[,25:37]
# jeps4=jeps[,38:50]
# jeps5=jeps[,51:63]
jpeg(paste("Correlations plotted1.jpg"), width = 8000, height = 8000, quality = 100,pointsize = 23, res=300); 
chart.Correlation(jeps, histogram=FALSE, pch=16)
# chart.Correlation(jeps1, histogram=TRUE, pch=12)
dev.off()
jpeg(paste("Correlations plotted2.jpg"), width = 8000, height = 8000, quality = 100,pointsize = 23, res=300); 
chart.Correlation(jepsf, histogram=FALSE, pch=16)
# chart.Correlation(jeps2, histogram=TRUE, pch=12)
dev.off()
jpeg(paste("Correlations plotted3.jpg"), width = 8000, height = 8000, quality = 100,pointsize = 23, res=300); 
chart.Correlation(jepsm, histogram=FALSE, pch=16)
# chart.Correlation(jeps3, histogram=TRUE, pch=12)
dev.off()
jpeg(paste("Correlations plotted4.jpg"), width = 8000, height = 8000, quality = 100,pointsize = 23, res=300); 
chart.Correlation(jeps4, histogram=TRUE, pch=12)
dev.off()
jpeg(paste("Correlations plotted5.jpg"), width = 8000, height = 8000, quality = 100,pointsize = 23, res=300); 
chart.Correlation(jeps5, histogram=TRUE, pch=12)
dev.off()



#The network plots: :)

#With correlations
tv_c=tv_covscl#[,9:dim(tv_covscl)[2]] #tv_half_log22 #cbind(tv[,1:8], tv_half_log2) #check also not logged and then the auto one
tv_c=tv_c[,!colnames(tv_c) %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]; tv_c=tv_c[,!colnames(tv_c) %in% x4]
colnames(tv_c)[colnames(tv_c)=="17aOH-P4"]="17a-OHP4"
dat = tv_c; 
dat=dat[,!colnames(dat) %in% c('Gender','PatientNumber')] #SEX.1F.2M
resulta <- (rcorr(as.matrix(dat), type = c('spearman')))$r
n_level=0.9
Nrr=qpNrr(resulta, verbose=FALSE);Nrr[is.na(Nrr)]=1;print(hist(as.numeric(Nrr),breaks=50)); cond=data.frame(as.matrix(Nrr<n_level))
RN=data.frame(resulta);tes_t=cond*RN;tes_t=as.matrix(tes_t);resulta=tes_t;colnames(resulta)=rownames(resulta) #https://www.geeksforgeeks.org/elementwise-matrix-multiplication-in-r/
# tv_ah <- prep.autoscale(resulta, center = TRUE, scale = TRUE);
# tv_c=tv_covscl#tv_half_log22 #cbind(tv[,1:8], tv_half_log2) #check also not logged and then the auto one
# # tv_c[,'Menopause']=as.numeric(sick_group)
# tv_c=tv_c[,c(1:3,4:(dim(tv_c)[2]))]#dim(tv_c)[2],
# tv_c=tv_c[,!colnames(tv_c) %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]
# tv_c=tv_c[,!colnames(tv_c) %in% x4]
# colnames(tv_c)[colnames(tv_c)=="17aOH-P4"]="17a-OHP4"
# dat = tv_c; 
# dat=dat[,!colnames(dat) %in% c('Gender','PatientNumber')] #SEX.1F.2M
# resulta <- (rcorr(as.matrix(dat), type = c('spearman')))$r
tes_t=resulta

a=length(x1)-2;b=length(x2);c=length(x3);d=length(x4);e=length(x5);f=length(x6);
# # removing self-correlation
# tes_t[1:b,1:b]=0
# tes_t[(b+1):(c+b),(c+1):(c+b)]=0
# tes_t[(c+b+1):(d+b+c),(c+b+1):(d+b+c)]=0
# tes_t[(a+b+c+1):(a+b+c+d),(a+b+c+1):(a+b+c+d)]=0
# tes_t[(a+b+c+d+1):(a+b+c+d+e),(a+b+c+d+1):(a+b+c+d+e)]=0
# tes_t[(a+b+c+d+e+1):(a+b+c+d+e+f),(a+b+c+d+e+1):(a+b+c+d+e+f)]=0 
# removing self-correlation
tes_t[1:a,1:a]=0
tes_t[(a+1):(a+b),(a+1):(a+b)]=0
tes_t[(a+b+1):(a+b+c),(a+b+1):(a+b+c)]=0
# tes_t[(a+b+c+1):(a+b+c+d),(a+b+c+1):(a+b+c+d)]=0
tes_t[(a+b+c+1):(a+b+c+e),(a+b+c+1):(a+b+c+e)]=0
tes_t[(a+b+c+e+1):(a+b+c+e+f),(a+b+c+e+1):(a+b+c+e+f)]=0
# resulta=resulta[colnames(resulta)[7:69],colnames(resulta)[7:69]]
# tes_t=tes_t[colnames(tes_t)[7:69],colnames(tes_t)[7:69]]

tes_t=tes_t[colnames(tes_t)[7:66],colnames(tes_t)[7:66]]
colnames(tes_t)
# tes_t=resulta #ks. correlations
g <- graph_from_adjacency_matrix(tes_t, mode="upper", weighted=TRUE, diag=FALSE)
e <- as_edgelist(g); df <- as.data.frame(cbind(e,E(g)$weight)); #
df[,3]=as.numeric(df[, 3])
hoi=df
hoi[rev(order(hoi[,3])),]
# colnames(hoi)=c('V1','V2','V4')
hoi=hoi[!duplicated(hoi[,c(1,2)]),]


# g.b <- betweenness(g, directed = FALSE)



#With ACMEs
# rtot=all_all#rbind(rtot1,rtot2,rtot3,rtot4) #ks. mediation all_all list_of_files[[1]]
# head(rtot)
# # rt2=rtot[rtot[,'d0.p']<0.42,]
# # rt2=rt2[rt2[,'Proportion Mediated']>0.3,]
# rt2=rtot[rtot[,1]>0,]
# rt2=rt2[rt2[,1]>quantile(rt2[,1],0.75),]
# rt2=rt2[rt2[,10]>quantile(rt2[,10],0.5),]
# # strsplit(rownames(rt2)[116], ".", fixed = TRUE)[[1]]
# hoi=c();for (i in 1:dim(rt2)[1]) {hoi=append(hoi,scan(text=rownames(rt2)[i], what=""))}
# # hoi=as.data.frame(matrix(hoi, ncol = 4,  byrow = TRUE), stringsAsFactors = FALSE) #check this..
# hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE) #check this..
# 
# hoi[,2] <- gsub("\\.", "-", hoi[,2])
# # hoi[,3]=rt2[,9]
# hoi[,4]=rt2[,1]/2
# hoi[,5]=rt2[,1]/2
# # hoi=hoi[hoi[,3]=='Steatosis',c(1,2,4)] #hoi[,3]=='Steatosis' tibble
# hoi1=hoi[,c(1,3,4)]
# hoi2=hoi[,c(1,2,4)]
# colnames(hoi1)=c('V1','V2','V4')
# colnames(hoi2)=c('V1','V2','V4')
# hoi=rbind(hoi1,hoi2)
# hoi=hoi[!duplicated(hoi[,c(1,2)]),]


# library, https://r-graph-gallery.com/249-igraph-network-map-a-color.html
library(igraph)
# # create data:
links <- data.frame(
  source=c("A","A", "A", "A", "A","J", "B", "B", "C", "C", "D","I"),
  target=c("B","B", "C", "D", "J","A","E", "F", "G", "H", "I","I"),
  importance=(sample(1:4, 12, replace=T)))
colnames(hoi)=colnames(links)
links=hoi
# nodes <- data.frame(
#   name=LETTERS[1:10],
#   carac=c( rep("young",3),rep("adult",2), rep("old",5)))
# sources=hoi %>% distinct('V1') %>% rename(source='label')
# destinations=hoi %>% distinct(V2) %>% rename(target ='label')
sources=hoi %>% distinct(source) %>% rename(source='label')
destinations=hoi %>% distinct(target) %>% rename(target ='label')

# destinations$label=c(xb,xl,xs)

nodess <- full_join(sources, destinations, by = "label")

# as.vector(nodess)$label

xc=x5[x5 %in% nodess$label]
xb=x3[x3 %in% nodess$label]
xl=x6[x6 %in% nodess$label]
xs=x2[x2 %in% nodess$label]

# x2[x2 =='17aOH-P4']='17a-OHP4' #Next time check these early on... :)

nodess$label=c(xc,xb,xl,xs)

# str_detect(a, regex(u, ignore_case = TRUE)) #not like this... does not work the way you want
# a <- "_L"
# u=as.vector(nodess)$label
# ie=grep(a,u,fixed=TRUE)
# check=c(8:27)
# hips=check[! check %in% ie]
# u2=u[ie];u3=u[hips]
# u4=c(u2,u3)
# nodess[check,1]=u4
# ls=length(u)-max(check)

nodes <- data.frame(
  name=nodess[,1], #  carac=c( rep("Contaminants",length(7)),rep("Bile Acids",length(u2)),rep("Lipids",length(u3)),rep("Steroids",ls))) #range on kaikki +1
  carac=( c(rep("Contaminants",length(xc)),rep("Bile Acids",length(xb)),rep("Lipids",length(xl)),rep("Steroids",length(xs))))) #range on kaikki +1
# nodess[,1]=c(x2,x3,x5,x6[1:(length(x6)-3)])
# nodes <- data.frame(name=nodess[,1], carac=c(  rep("Steroids",20),rep("Bile Acids",23),rep('Contaminants',7),rep("Lipids",10))) #range on kaikki +1, -rep("Covariates",6), some lipids
# nodes[,1]=c(x2,x3,x5,x6[1:(length(x6)-3)])

# nodes[42:48,2]='Contaminants'
# nodes[49:60,2]='Lipids'

# Turn it into igraph object
network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
# Make a palette of 3 colors
library(RColorBrewer)
coul  <- c('#B2BEB5','Green','Red','Orange') # brewer.pal(2, "Set1") ,'#ADD8E6','#faf0e6'
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network)$carac))]
# Make the plot
# plot(network, vertex.color=my_color)
# Add a legend
# legend("bottomleft", legend=levels(as.factor(V(network)$carac))  , col = coul , bty = "n", pch=20 , pt.cex = 3, cex = 0.5, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))
# Plot
# plot(network, mode = "circle",vertex.color=my_color, vertex.size = 07, 
#      edge.arrow.size = 0.8,vertex.label.cex = 0.65,edge.width=as.numeric(E(network)$importance)*30 )
# legend("topright", legend=levels(as.factor(V(network)$carac))  , 
#        col = coul , bty = "n", pch=20 , pt.cex = 1.3, cex = 1.3, text.col=coul , 
#        horiz = FALSE, inset = c(0.01, 0.01))



plot(network, mode = "circle",vertex.color=my_color, vertex.size = 10, 
     edge.arrow.size = 0.8,vertex.label.cex = 0.45,edge.width=as.numeric(E(network)$importance)*6.00 )

legend("topright", legend=levels(as.factor(V(network)$carac))  , 
       col = coul , bty = "n", pch=20 , pt.cex = 1.3, cex = 1.3, text.col=coul , 
       horiz = FALSE, inset = c(0.01, 0.01))


#make an edgelesti:
# pheatmap(jeg,  cluster_rows=FALSE, cluster_cols=FALSE)
# nodes <- full_join('V1', 'V2', by = "label") #https://www.jessesadler.com/post/network-analysis-with-r/

# sources      =      hoi %>% distinct(V1) %>% rename(V1='label')
# destinations =      hoi %>% distinct(V2) %>% rename(V2='label')
# nodess <- full_join(sources, destinations, by = "label")
# nodes <- rowid_to_column(nodess, "id")
# nodes
# per_route <- hoi %>%  group_by('V1','V2') %>% summarise(hoi,n= n(), .groups = "drop")
# per_route
# # %>%
# #   summarise(weight = n(), .groups = "drop")
# hoi
# hoi[order(hoi[,1]),]
# 
# edges <- select(hoi, V1, V2, V4)
# edges <- edges %>% 
#   left_join(nodes, by = c("V1" = "label")) %>% 
#   rename(from = id)
# edges <- edges %>% 
#   left_join(nodes, by = c("V2" = "label")) %>% 
#   rename(to = id)
# library(network)
# routes_network <- network(tibble(edges),
#                           vertex.attr = nodes,
#                           matrix.type = "edgelist",
#                           ignore.eval = TRUE, multiple=TRUE)
# 
# library(ggsankey)
# df <- mtcars %>% make_long(cyl, vs, am, gear, carb)
# df=hoi
# 
# rango <- function(x){((x-min(x))/(max(x)-min(x)))*2-1} #just a function for the -1 to 1 thing..
# col_fun = colorRamp2(as.numeric(c(min(df$V4), 0,max(df$V4))), c("blue",'white', "red"))
# df=df[!df$V1 %in% rem,];df=df[!df$V2 %in% rem,] #e.g.rem=x4
# df$V4=rango(as.numeric(df$V4));
# classes=modi #modi=4
# grid.col=grid.col[grid.col!=colt]# ='black'
# namesh=unique(g1)[1:classes];cola=unique(grid.col)[1:classes]
# lgd_group = Legend(at = genders[i], type = "points", legend_gp = gpar(col = colors[i]), title_position = "topleft", title = title)
# lgd_points = Legend(at = namesh, type = "points", legend_gp = gpar(col = cola), title_position = "topleft", title = "Class")
# lgd_lines = Legend(at = c("Positive", "Negative"), type = "points", legend_gp = gpar(col = c('red','blue')), title_position = "topleft", title = "Correlation")
# lgd_edges= Legend(at = c(min(df$V3), max(df$V3)), col_fun = col_fun,  title_position = "topleft", title = "Edges")
# lgd_list_vertical = packLegend(lgd_group,lgd_points,  lgd_lines,lgd_edges) #lgd_lines,
# chordDiagram(df,  transparency = 0.5,col = col_fun,)
# circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
#   xlim = get.cell.meta.data("xlim"); ylim = get.cell.meta.data("ylim")
#   sector.name = get.cell.meta.data("sector.index")
#   circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
#   circos.axis(h = "top", labels.cex = 0.000001, major.tick.length = 0.2, sector.index = sector.name, track.index = 2)}, bg.border = NA) 



# poissone$coefficients;
# hist(poissone$coefficients)
# poissone$coefficients

# valuees=summary(poissone)
# ms=poissone$coefficients[2:21]
# error=valuees$coefficients[2:21,2]
# pval=anova(poissone)[1:20,5]
# error_lower=ms-error# error_lower[error_lower <= 0] = 0
# error_upper=ms+error
# hip=colnames(tv[,9:28])
# sample_data <- data.frame(study=colnames(tv[,9:28]),index=colnames(tv[,9:28]),result=ms,error_lower=error_lower,error_upper=error_upper)#,pval=pval)
# sample_data %>%
#   mutate(study = fct_reorder(study, result)) %>%
#   ggplot(aes(y=study, x=ms,xmin=error_lower,xmax=error_upper)) +
#   geom_point(color= "black", pch= 1, size=3) +
#   # coord_flip() +
#   geom_errorbarh(height=.5, color= "black", lwd=0.5) +
#   labs(title='Steatosis Grade Effect Forest Plot', x='LM Estimates (SD)', y = 'Steroids')+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = 
#           element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))#+
# 

# %%%%....


# metad=data.frame(colnames(SG0)[2:21],'cate') #https://datatofish.com/create-dataframe-in-r/
# colnames(metad) =c('name','group'); #if you know: ga=groups[,'Abbreviation']; md=metad[,'name']; unique(c(ga,md)); ga[!(ga %in% md)]; ga[ga=="17a-OHP4"]="17aOH-P4"
# groups[,'Abbreviation'][groups[,'Abbreviation']=="17a-OHP4"]="17aOH-P4"
# groups=groups[groups[,'Abbreviation']!='F',];rownames(groups)=1:20
# groups=groups[order(groups[,'Abbreviation']),]; metad=metad[order(metad[,'name']),]
# metad[,'group']=groups[,'Group']# groups[,'Abbreviation'] == metad[,'name']
# cate=metad[,'group']
# hup=names(table(cate)[rev(order(table(cate)))])#[c(1:15)]) #2:3, 1, 4, 5:8, 9-19
# metad2=metad[metad[,'group'] %in% hup,]
# df_grouping=tibble(metad2) #should be ok
# sample_data2=sample_data[rev(order(sample_data[,'result'])),]
# sample_data2=tibble(sample_data2)
# colnames(sample_data2) =c('name','name2','result','nerror','error','ok')# Join the association data frame df_compare_traits with group data # use right_join, with df_grouping on the right, to preserve the order of biomarkers it specifies.
# df_compare_traits_groups <- sample_data2 %>% dplyr::right_join(., df_grouping, by = "name") %>% dplyr::mutate(group = factor(.data$group, levels = unique(.data$group)))
# df_compare_traits_groups=data.frame(df_compare_traits_groups)
# df_compare_traits_groups=df_compare_traits_groups[rev(order(df_compare_traits_groups[,'result'])),]
# df_compare_traits_groups=tibble(df_compare_traits_groups)
# 
# gn=df_compare_traits_groups[,c('group','name')]
# # gn[order(data.frame(gn[,'name'])[,1]),]
# gn[order(data.frame(gn[,'name'])[,1]),]
# sample_data=sample_data[order(sample_data[,'study']),]
# sample_data=cbind(sample_data,gn[order(data.frame(gn[,'name'])[,1]),])
# # sample_data[,'Group']='Group'
# # groups$Abbreviation[groups$Abbreviation == '17a-OHP4']='17aOH-P4'
# # for (i in 1:21) {sample_data[sample_data$study %in% groups$Abbreviation[i],'Group']=groups$Group[i]}
# # groups=groups[groups[,'Abbreviation']!='F',]
# # # for (i in 1:20) {groups[groups$Abbreviation %in% sample_data$study[i],'Group']=sample_data$Group[i]}
# sample_data=sample_data[order(sample_data[,'study']),]
# sample_data=cbind(sample_data,groups[order(data.frame(groups[,'Abbreviation'])[,1]),])
# 
# cm=sign(sample_data[,'error_lower']) < 1 & sign(sample_data[,'error_upper']) < 1
# cp=sign(sample_data[,'error_lower']) > 0 & sign(sample_data[,'error_upper']) > 0
# ct= cm | cp
# sample_data[,'Significance']=ct
# sample_data[,'Significance'][sample_data[,'Significance']==TRUE]='Yes'
# sample_data[,'Significance'][sample_data[,'Significance']==FALSE]='No'
# sample_data[,'Color']=ct#!(sample_data[,'error_lower'] >= 0 | sample_data[,'error_upper'] >= 0)
# sample_data[,'Color'][sample_data[,'Color']==TRUE]='blue'
# sample_data[,'Color'][sample_data[,'Color']==FALSE]='grey'
# sample_data[,'pval']=ct#!(sample_data[,'error_lower'] >= 0 | sample_data[,'error_upper'] >= 0)
# sample_data[,'pval'][sample_data[,'pval']==TRUE]=0.01
# sample_data[,'pval'][sample_data[,'pval']==FALSE]=0.06
# sample_data2=sample_data
# # sample_data2=sample_data2[,c(7,6,3:5,9:11)]
# sample_data2[,unique(colnames(sample_data2))]
# plot=forestplot(df = sample_data, #drive tba_example_v3_oh_tikka17823.R if not working via 'x' error
#                 estimate = result,
#                 se= abs(error_lower-error_upper)/4,
#                 pvalue = pval, #this makes the significant value..
#                 psignif = 0.1,
#                 xlab = "LM  Estimates (SD)",
#                 ylab='Steroid Groups',
#                 title='Forest Plot of Linear Model Estimates for Male HOMAIRs',
#                 colour = Significance) +
#   ggforce::facet_col(
#     facets = ~group,
#     scales = "free_y",
#     space = "fixed",
#     strip.position='left') 
# if (Group=='All') {ordera=sample_data$study[order(sample_data$result)] #this should be ok
# plot[["data"]][["study"]]=factor(plot[["data"]][["study"]], levels = sample_data$study[order(sample_data$result)] )} else if 
# (Group!='All') {plot[["data"]][["study"]]=factor(plot[["data"]][["study"]], levels = sample_data$study[order(sample_data$result)] )}
# plot$layers[[1]]$aes_params$odd <- "#00000000" 
# #https://stackoverflow.com/questions/71745719/how-to-control-stripe-transparency-using-ggforestplot-geom-stripes
# jop=plot #+theme(axis.text.y=element_blank());
# jop2=jop+geom_point(aes(colour = factor(Significance)),colour = sample_data[,'Color']) +
#   scale_color_manual(values=c('#999999','blue'));jop2 

# plot$layers[[1]]$aes_params$odd <- "#00000000" #https://stackoverflow.com/questions/71745719/how-to-control-stripe-transparency-using-ggforestplot-geom-stripes
#https://rdrr.io/rforge/CALIBERdatamanage/man/multiforest.html

#https://www.nature.com/articles/s41467-022-30875-7#citeas #with https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html
#other # https://academic.oup.com/bioinformatics/article/26/18/2336/208507?login=false#393579283 http://locuszoom.sph.umich.edu//
#https://www.cell.com/ajhg/fulltext/S0002-9297(07)61352-4#secd28548230e6399
#https://www.mdpi.com/2218-273X/13/6/917


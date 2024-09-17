


#Mediation analysis with the functions and packages mentioned in 'packages_data load_functions_esp. for the mediat._tikka26124.R' (or the like) file at the GitHub.
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/") #check this if needed...
date='tikka18624' # but do change this... :)

# save.image("figs okei hei hou tikka14624.RData")

# load("fig5_6 okei tikka20524.RData") #xien kaa, when you have it load it!
#Fyi: some Finnish expressions may occur. 
# Elikkäs pelikkäs, testataan tällä steroid data uusine ja vanhoine olettamuksineen...updated 26.1.2024 (if in conflict, see the save date of the file)
#Kokeillaan nää ulostulot erikseen:
ccova=tv[,c("Steatosis.Grade.0.To.3" , "Fibrosis.Stage.0.to.4" ,"Necroinflammation" ,  "HOMA-IR")]
sick_group=rowSums(ccova)>4
file_names=c("Steatosis" , "Fibrosis" ,"Necroinflammation" ,  "HOMAIR", 'Menopause')
lkm=30; simss=100 # 
joo='ei';ip=1
ccovae=tv[,c("Steatosis.Grade.0.To.3")]; sick_group=ccovae>0 #toth# # hist(ccovae,breaks=100) # hist(ccova[,'HOMA-IR'],breaks=100)
fn=file_names[1]; 
t.val='no'; name='generic'; sick='no'; Group='All';joo='joo';ip=1 #group all by default


Group = ' all'; cond=''

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

colnames(tv_all)
Outcome
Treatment
Mediator
# Treatment=c(Treatment, 'Perfluorodecyl.ethanoic.acid') #testing this now that the function and definitions are clearer...

#All
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/Uusix maksassa_all/ok_mcer/")
#These is generic...

name='total'; fn='All'; sick_group=rowSums(ccova)>4; joo ='joo'; sick='no'; lkm=30; simss=30;ip=1 
#Plots for all are needed...#remember ip definiation; https://stackoverflow.com/questions/62900552/grid-seach-on-arima-model-in-r
tv_all=tv
Group='All';    uh7maa=loop_med_simplified_hyp4a(Group,tv,name,date,simss,sick,sick_group,joo,ip); save(uh7maa, file=paste(fn,'all.RData')) #lp
Group='Female'; uh7fee=loop_med_simplified_hyp4a(Group,tv_all,name,date,simss,sick,sick_group,joo,ip); save(uh7fee, file=paste(fn,'female.RData')) #lp  #try({uh7fee},{uh7fee=data.frame(0)}) 
Group='male';   uh7mee=loop_med_simplified_hyp4a(Group,tv_all,name,date,simss,sick,sick_group,joo,ip); save(res, file=paste(fn,'male.RData')) # lp #try({uh7fee},{uh7fee=data.frame(0)})

# date=paste(date,'_female'); 
# Group='Female';
# allkoe=reduced2(uh7fee,Group,name,lkm);
# plottings_sf(allkoe,date,sick,Group)
# date=paste(date,'_male'); 
# Group='Male';
# allkom=reduced2(uh7mee,Group,name,lkm);
# plottings_sf(allkom,date,sick,Group)
# date=paste(date,'_all'); 
# Group='All'; 
# allk=reduced2(uh7maa,Group,name,lkm);
# plottings_sf(allk,date,sick,Group)

#Steatosis
simss=100;
date='tikka22824' # but do change this... :)
name='Jaot_OK_basica'
joo='ei';ip=1
tv_all=tv_covscl
ccovae=tv[,c("Steatosis.Grade.0.To.3")]; sick_group=ccovae>0 #toth# # hist(ccovae,breaks=100) # hist(ccova[,'HOMA-IR'],breaks=100)
fn=file_names[1]; 
sick='yes'
the_essentials(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group,fn,lkm,date,joo,ip)

#Fibrosis
ccovae=tv[,c("Fibrosis.Stage.0.to.4")]; 
sick_group=ccovae>0 #toth#
fn=file_names[2]; 
the_essentials(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group,fn,lkm,date,joo,ip)
#For menopause
# ccovae=tv[,c("AGE")]; c1=ccovae<56; c2=ccovae>44; sick_group=c1 & c2 #tv[,c("AGE")][c1 & c2]
# fn=file_names[5]; 
# the_essentials(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group,fn,lkm,date,joo,ip)
# #Necroinfl.
ccovae=tv[,c("Necroinflammation")]; sick_group=ccovae>0 #toth#
fn=file_names[3];
the_essentials(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group,fn,lkm,date,joo,ip)
#Homa #remember to always test the funciton with minimun number setting or as light parameters as possible to get it through... before big runs
joo='ei';sick='yes'
ccovae=tv[,c("HOMA-IR")]; sick_group=ccovae>1.5 #toth#
fn=file_names[4]; 
the_essentials(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group,fn,lkm,date,joo,ip)

#All; tämä lienee toistoa, mutta ok cross check:
joo='joo';sick='no'
ccova=tv[,c("Steatosis.Grade.0.To.3", "Fibrosis.Stage.0.to.4" ,"Necroinflammation","HOMA-IR")] 
fn='All'; sick_group=rowSums(ccova)>4 #toth#
the_essentials(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick, sick_group,fn,lkm,date,joo,ip)
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/") #check this if needed...
save.image("okei_tikka13824_nonia.RData")
# #MASLD
joo='ei'
ccova=tv[,c("Steatosis.Grade.0.To.3", "Fibrosis.Stage.0.to.4" ,"Necroinflammation","HOMA-IR")]
fn='MASLD'; sick_group=rowSums(ccova)>4 #toth#
the_essentials(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick, sick_group,fn,lkm,date,joo,ip)


setwd("C:/Users/patati/Desktop/TurkuOW/RWork/MILP2") #check this if needed...
#The mediation with the lp method without perfuoro näyttää oklta
simss=100; name='healthy_lp'; sick='no'; fn='All'; sick_group=rowSums(ccova)>4; lkm=30
Group='All'; uhaa=lp_fn(tv_all,Group,sick,sick_group,fn,name); #uha=uha[rev(order(uha[,'ACME'])),]; 
Group='Male'; uham=lp_fn(tv_all,Group,sick,sick_group,fn,name)
Group='Female'; uhaf=lp_fn(tv_all,Group,sick,sick_group,fn,name)

Group='All'; name='lpo_helathy';
#date='tikka14424_lpa'
alma=reduced2(uhaa,Group,name,lkm);
plottings_sf(alma,date,sick,Group) # uhaf2=uhaf[rev(order(uhaf[,'ACME'])),]; 

#date='tikka14424_lpm';
lkm=30
name='lpok2';Group='male';
almam=reduced2(uham,Group,name,lkm);
almam = na.omit(almam)
plottings_sf(almam,date,sick,Group)

#date='tikka14424_lpf';
name='lpok2';Group='female';
almaf=reduced2(uhaf,Group,name,lkm);
almaf = na.omit(almaf) #https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
plottings_sf(almaf,date,sick,Group)

simss=100; Group='All'; 
name='sicka'; sick='yes'; fn='All'; sick_group=rowSums(ccova)>4; lkm=30
uhax=lp_fn(tv_all,Group,sick,sick_group,fn,name); #uha=uha[rev(order(uha[,'ACME'])),]; 
Group='Male'; uhamx=lp_fn(tv_all,Group,sick,sick_group,fn,name)
Group='Female'; uhafx=lp_fn(tv_all,Group,sick,sick_group,fn,name)

# uhax=uha

Group='All'; name='lpo_sick';
#date='tikka14424_lp'
alma=reduced2(uhax,Group,name,lkm);
plottings_sf(alma,date,sick,Group)
#date='tikka14424_lpm2_s';
lkm=30
name='lpok2';Group='male';
almam=reduced2(uham,Group,name,lkm);
almam = na.omit(almam)
plottings_sf(almam,date,sick,Group)
#date='tikka14424_lpf';
name='lpok2';Group='female';
almaf=reduced2(uhaf,Group,name,lkm);
almaf = na.omit(almaf) #https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
plottings_sf(almaf,date,sick,Group)


setwd("C:/Users/patati/Desktop/TurkuOW/RWork/")
save.image("basic_med14424a.RData")

#for the steroids:2
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6e/cova_14124") #tests.... 30.11.2023
t.val='minmax';
ccova=tv[,c("Steatosis.Grade.0.To.3", "Fibrosis.Stage.0.to.4" ,"Necroinflammation","HOMA-IR")] 
name='100basic_okcovah'
sick_group=rowSums(ccova)>4 #toth#
simss=10; test=''; take='';  sick_group=sick_group
sick='yes';name='100basic_okcovah_sick'
Group='All';    uh7ma=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Female'; uh7f=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Male';   uh7m=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);

sick='no';name='100basic_okcovah_healthy' #these are the results you need for the example, i.e. fig4 in the manu (2524):
Group='All';    uh7mae=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Female'; uh7fe=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Male';   uh7me=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);

save.image("100basic_cova_25424.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic_cova") #tests6.... 30.11.2023
# load("100basic_cova_141223.RData")
lkm=30;Group='All'; name='sick';date='tikka25424'
alma=reduced2(uh7ma,Group,name,lkm);
plottings_sf(alma,date,sick,Group)

date='tikka25424_f';name='sick';Group='female';
almam=reduced2(uh7f,Group,name,lkm);
almam = na.omit(almam)
plottings_sf(almam,date,sick,Group)

date='tikka25424_m';name='sick';Group='male';
almaf=reduced2(uh7m,Group,name,lkm);
almaf = na.omit(almaf) #https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
plottings_sf(almaf,date,sick,Group)

lkm=30;Group='All'; name='_not sick';date='tikka25424'
alma=reduced2(uh7mae,Group,name,lkm);
plottings_sf(alma,date,sick,Group)

date='tikka25424_f';name='_not sick';Group='female';
almam=reduced2(uh7fe,Group,name,lkm);
almam = na.omit(almam)
plottings_sf(almam,date,sick,Group)

date='tikka25424_m';name='_not sick';Group='male';
almaf=reduced2(uh7me,Group,name,lkm);
almaf = na.omit(almaf) #https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
plottings_sf(almaf,date,sick,Group)


#for the booting with hypo4
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/boot_oh") #tests.... 30.11.2023
t.val='minmax';
ccova=tv[,c("Steatosis.Grade.0.To.3", "Fibrosis.Stage.0.to.4" ,"Necroinflammation","HOMA-IR")] 
name='100boot_hypo4_sick'
sick_group=rowSums(ccova)>4 #toth#
simss=100;   sick_group=sick_group
sick='yes';
Group='All';    uh7ma=loop_med_simplified_hyp4a(Group,name,date,simss,sick,sick_groupe,joo,ip);
Group='Female'; uh7f=loop_med_simplified_hyp4a(Group,name,date,simss,sick,sick_groupe,joo,ip);
Group='Male';   uh7m=loop_med_simplified_hyp4a(Group,name,date,simss,sick,sick_groupe,joo,ip);

sick='no';name='100boot_hypo4_healthy' #these are the results you need for the example, i.e. fig4 in the manu (2524):
Group='All';    uh7mae=loop_med_simplified_hyp4a(Group,name,date,simss,sick,sick_groupe,joo,ip);
Group='Female'; uh7fe=loop_med_simplified_hyp4a(Group,name,date,simss,sick,sick_groupe,joo,ip);
Group='Male';   uh7me=loop_med_simplified_hyp4a(Group,name,date,simss,sick,sick_groupe,joo,ip);

save.image("100basic_cova_25424.RData")
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic_cova/the all") #tests6.... 30.11.2023
# load("100basic_cova_141223.RData")
lkm=30;Group='All'; name='sick';date='tikka25424'
alma=reduced2(uh7ma,Group,name,lkm);
plottings_sf(alma,date,sick,Group)

date='tikka25424_f';name='sick';Group='female';
almam=reduced2(uh7f,Group,name,lkm);
almam = na.omit(almam)
plottings_sf(almam,date,sick,Group)

date='tikka25424_m';name='sick';Group='male';
almaf=reduced2(uh7m,Group,name,lkm);
almaf = na.omit(almaf) #https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
plottings_sf(almaf,date,sick,Group)

lkm=30;Group='All'; name='_not sick';date='tikka25424'
alma=reduced2(uh7mae,Group,name,lkm);
plottings_sf(alma,date,sick,Group)

date='tikka25424_f';name='_not sick';Group='female';
almam=reduced2(uh7fe,Group,name,lkm);
almam = na.omit(almam)
plottings_sf(almam,date,sick,Group)

date='tikka25424_m';name='_not sick';Group='male';
almaf=reduced2(uh7me,Group,name,lkm);
almaf = na.omit(almaf) #https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
plottings_sf(almaf,date,sick,Group)



# https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

# Group='male'
# date='tikka1324_lpme2'
# uhams=reduced3(uham,Group,name,lkm=50);
# plottings_sf(uhams,date,sick,Group)
# date='tikka1324_lpfe2'
# uhafs=reduced3(uhaf,Group,name,lkm=50);
# plottings_sf(uhafs,date,sick,Group)
# order(jups[,'ACME']())
# median(jups[,'ACME'][jups[,'ACME']>0])

# uhaf=uham
# # uhaf_ok=uhaf
# uhaf2=uhaf[rev(order(uhaf[,'ACME'])),]; 
# uhaf2 = uhaf2[!grepl("Perfluorodecyl.ethanoic.acid",rownames(uhaf2)),]
# uhaf2=uhaf2[uhaf2[,'Proportion Mediated']>median(uhaf2[,'Proportion Mediated']),]
# uhaf2=uhaf2[uhaf2[,'d0.p']<0.7,]
# uhaf2=uhaf2[uhaf2[,'ADE']<median(uhaf2[,'ADE']),]
# uhaf2=uhaf2[uhaf2[,'z0.p']>0.05,]
# # AV=scale(abs(uhaf2[,'d0.ci_l']-uhaf2[,'d0.ci_u'])) #I would not take this... 12224
# # ACMEVar=round(quantile(AV,0.95,na.rm=TRUE),1) #I would not take this...
# # uhaf2=uhaf2[AV<ACMEVar,]
# PV=abs(uhaf2[,'n.ci_l']-uhaf2[,'n1.ci_u'])
# PMEDVar=round(quantile(PV,0.95,na.rm=TRUE),1) #not sure if this is needed
# uhaf2=uhaf2[ PV<PMEDVar,];
# ADV=abs(uhaf2[,'z0.ci_l']-uhaf2[,'z0.ci_u'])
# ADEVar=round(quantile(PV,0.95,na.rm=TRUE),1) #not sure if this is needed
# uhaf2=uhaf2[ ADV<ADEVar,];
plottings_sf(uhaf2,date,sick,Group)

date='tikka23224'
jups2=jups
jups2=jups2[grep(pattern = "A_S", x = rownames(jups2)),]
allke=reduced2(jups2,Group,name,lkm);
plottings_sf(allke,date,sick,Group)

date='tikka20224_hypo4m'
hams=reduced2(uh7m,Group,name,lkm=30);
plottings_sf(hams,date,sick,Group)

date='tikka20224_hypo4f'
hafs=reduced2(uh7f,Group,name,lkm=30);
plottings_sf(hafs,date,sick,Group)

# name=
# rtot=rtot_2000_mrct # rtot=uh5
med.freq.cut=0.3
med.min=0
med.prop=0.4
name=paste('Contaminants_Steroids_BAs_or_Lipids_sims_basic direct_alleh',test,take,date) 
lkm
load("All uh7maee.RData")
rtot = uh7mae
# direct_sankey(rtot, med.freq.cut=0.3, med='no',med.min=0, med.prop=0.4,name,lkm) 

#Tarkastelma
rt2=rtot
# med.freq.cut=0.4; med.min=0; med.prop=0.6;name;lkm
# 
# rt2=rtot[rtot[,'ADE']>med.min,];
# rt2=rt2[rt2[,'Proportion Mediated']<med.prop,];
# rt2=rt2[rt2[,'z0.p']<med.freq.cut,]
# # rt2=rt2[rt2[,'d0.p']>med.freq.cut,];
# u3=rtot
# u3=u3[rev(order(u3[,'ACME'])),]; 
# # cut=round(quantile(abs(u3[1:200,'d0.ci_l']-u3[1:200,'z0.ci_u']),0.25),1)
# ACMEMedian=c();ACMEpval=c();ACMEVar=c()
# ADEMedian=c();ADEpval=c();ADEVar=c()
# PMEDMedian=c();PMEDpval=c();PMEDVar=c(); c1=c()

c1=rt2
# ADEMedian=median(c1[,'ADE'])+sd(c1[,'ADE'])
# c1= c1[c1[,'ADE'] > ADEMedian  ,]
cute=c()
for (i in 1:dim(c1)[1]) {if (((c1[i,'z0.ci_u']-c1[i,'d0.ci_l']) > -1.5) & (c1[i,'ADE']-c1[i,'ACME']) > -1.5 & c1[i,'Proportion Mediated'] < 0.9) {cute=append(cute,i)} else {next}} #c1[i,'d0.p'] > 0.1 &  
c1=c1[1:dim(c1)[1] %in% cute,]
c1=c1[rev(order(c1[,'ADE'])),]; 
c1=c1[1:100,];

hoi=c(); hoi=scan(text=rownames(c1), what="")
hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender')
hoi=hoi[,c('Contaminants','Bile Acids or Lipids')]
hoi=hoi[!duplicated(hoi[c(1,2)]),]
## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct
# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation
# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
df2 <- hoi %>%make_long('Contaminants','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
meda='Sankey plot ofa7'
pdf(paste(meda,name,med.freq.cut,med.prop ,"all.pdf"), width = 20, height = 20,  pointsize = 10);
print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
jpeg(paste(meda,name,med.freq.cut,med.prop ,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 12, res=500);
print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
        geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
        scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()


# direct_sankey=function(rtot, med.freq.cut, med.min, med,med.prop,name,lkm) {
# rt2=rtot
# if (med.min=='All') {rt2=rtot} else if (med.min!='All')  {if (med=='Mediated') 
#   {rt2=rtot[rtot[,'ACME']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']>med.prop,];rt2=rt2[rt2[,'d0.p']<med.freq.cut,];rt2=rt2[rt2[,'z0.p']>med.freq.cut,]} else 
#   {rt2=rtot[rtot[,'ADE']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']<med.prop,];rt2=rt2[rt2[,'d0.p']>med.freq.cut,];rt2=rt2[rt2[,'z0.p']<med.freq.cut,]}}
# rt2=rt2[rev(order(rt2[,'ADE'])),]; 
# rt2=rt2[lkm,]
# hoi=c(); hoi=scan(text=rownames(rt2), what="")
# hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
# colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender')
# hoi=hoi[,c('Contaminants','Bile Acids or Lipids')]
# ## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct
# # https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation
# # https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
# df2 <- hoi %>%make_long('Contaminants','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
# meda='Sankey plot of'
# pdf(paste(meda,name,med.freq.cut,med.prop ,"all.pdf"), width = 20, height = 20,  pointsize = 10);
# print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
# jpeg(paste(meda,name,med.freq.cut,med.prop ,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 12, res=500);
# print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
#         geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
#         scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
# return(rt2)}
# direct_sankey(rtot, med.freq.cut=0.3,med.min=0, med='no',med.prop=0.4,name,lkm=30)

#This was fairly important for the means and errors:
Group='All';Outcome='Steatosis.Grade.0.To.3'
sample_dataa=tv_basic(NAFLD,Group,Outcome)

sample_dataa=data.frame(sample_dataa)

tv_basic=function(NAFLD,Group,Outcome) {
  if (Group=='Male') {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==2,]} else if (Group=='Female') {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==1,]} else if (Group=='All') {NAFLDo=NAFLD} 
  sample_data=c();n0=c();n1=c()
  for (i in 1:3) {
    if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] == 0,]; n0=dim(SG0)[1]} else if (i==2) { SG0=NAFLDo[NAFLDo[,Outcome] > 0,]; n1=dim(SG0)[1]} else if (i==3) {SG0=NAFLDo; n1=dim(SG0)[1]}
    means=c(); for (j in 8:27) {means=append(means,round(median(SG0[,j], na.rm=TRUE),1))} 
    sds=c();for (j in 8:27) {sds=append(sds,round(sd(SG0[,j],na.rm=TRUE)),2)} 
    error_lower=means-sds; error_upper=means+sds; error=sds
    sample_data <- append(sample_data,data.frame(study=colnames(NAFLD[,8:27]),index=colnames(NAFLD[,8:27]),result=means,error=error))}
  
  return(sample_data)}




#ok..#Need to redo figures 5 and 6...#...
#Load all the variables in the folder:

setwd("C:/Users/patati/Desktop/TurkuOW/RWork/Uusix maksassa_all/ok_mcer/") #check this if needed...
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo_basic/Fibrosis/") #check this if needed...
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic_cova/the all/") #check this if needed...
setwd('C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/Steatosis/ok_sofar_really so')

setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic/") #check this if needed...


# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/All/ok paitsi pea")
files <- list.files(pattern="*.RData")
ldf <- lapply(files, load)
list_of_files <- list() #create empty list
for (i in files) {print(i); list_of_files[[i]] <- get(load(paste0("", i)))}  #add files to list position#https://www.reddit.com/r/Rlanguage/comments/nq773b/reading_multiple_rdata_files_into_a_list/
names(list_of_files) <- files #https://stackoverflow.com/questions/38643000/naming-list-elements-in-r

# all: " uh7ma.RData"  3 files 
all_all=list_of_files$`Fibrosis sick all.RData` #e on healthy :)
# female:
all_fem=list_of_files$`Fibrosis sick female.RData`
# male:
all_male=list_of_files$`Fibrosis sick male.RData`

# all: " uh7ma.RData"  3 files 
all_all=list_of_files$`Steatosis sick all.RData` #e on healthy :) sick
# female:
all_fem=list_of_files$`Steatosis sick female.RData`
# male:
all_male=list_of_files$`Steatosis sick male.RData`


all_all=all_all[rev(order(all_all[,1])),];
all_fem=all_fem[rev(order(all_fem[,1])),];
all_male=all_male[rev(order(all_male[,1])),];

sick='sick'
lkm=30;Group='All'; name='just';date='tikka22824_alla'
alma=reduced2(all_all,Group,name,lkm); #all_all1
plottings_sf(alma,date,sick,Group)

date='tikka22824_femala';name='just';Group='female';
almaf=reduced2(all_fem,Group,name,lkm);
almaf = na.omit(almaf)
plottings_sf(almaf,date,sick,Group)

date='tikka22824_mala';name='just';Group='male';
almam=reduced2(all_male,Group,name,lkm);
almam = na.omit(almam) #https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
plottings_sf(almam,date,sick,Group)


# library(readxl)
#... all: " uh7ma.RData"  3 files 
# all_all=read_xlsx(path = "100basic_okcovaha All tikka5624 .xlsx") #e on healthy :)
# all_all=read_xlsx(path = "100basic_okcovaha All tikka5624 .xlsx") #e on healthy :)

# all_all=as.data.frame(all_all); 
# # all_all[,1][duplicated(all_all[,1])]
# dim(all_all)
# all_all=all_all[!is.na(all_all$ACME),]
# dim(all_all)
# all_all=na.omit(all_all)
# # dim(all_all)
# # rownames(all_all)=all_all[,1] #make.names(all_all[,1], unique = TRUE) 
# # all_all=all_all[,2:dim(all_all)[2]]; all_all=all_all[rev(order(all_all[,1])),]; 
# # hip=is.na(all_all); hip <- 1 * hip; # which(hip!=0)
# # all_all=all_all[rowSums(hip) == 0,]
# # dat <- dat %>% mutate(lon = as.numeric(lo))
# all_all1=all_all
# all_Female=all_all
# all_Male=all_all
# female:
# all_fem=list_of_files$`All female.RData`
# # male:# all_male=list_of_files$`All male.RData`
# save.image("100basic_cova_25424.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic_cova/the all") #tests6.... 30.11.2023
# load("100basic_cova_141223.RData")



#C:\Users\patati\Desktop\TurkuOW\RWork\Uusix maksassa_all\ok_mcer\alleees
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/Uusix maksassa_all/ok_mcer/alleees") #tests6.... 30.11.2023
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic/") #check this if needed...

library(readxl)
#... all: " uh7ma.RData"  3 files 
# all_all=read_xlsx(path = "100basic_okcovaha All tikka5624 .xlsx") #e on healthy :)
all_all=read_xlsx(path = "100basic Male tikka3624 .xlsx") #total Male tikka76524 hyp4b_oki.xlsx") #e on healthy :)

all_all=as.data.frame(all_all); 
dim(all_all)
all_all=all_all[!is.na(all_all$ACME),]
dim(all_all)
all_all=na.omit(all_all)
dim(all_all)
# rownames(all_all)=all_all[,1] #make.names(all_all[,1], unique = TRUE) 
# all_all=all_all[,2:dim(all_all)[2]]; 
# all_all=all_all[rev(order(all_all[,1])),]; 

all_all1=all_all
all_Female=all_all
all_Male=all_all

# load("100basic_cova_141223.RData")

sick='all samples'
lkm=30;Group='All'; name='just all';date='tikka4924_alla'
alma=reduced2(all_all1,Group,name,lkm);
plottings_sf(alma,date,sick,Group)

date='tikka4924_femala';name='just all';Group='female';
almaf=reduced2(all_Female,Group,name,lkm); #all_Female
almaf = na.omit(almaf)
plottings_sf(almaf,date,sick,Group)

date='tikka4924_mala';name='just all';Group='male';
almam=reduced2(all_Male,Group,name,lkm);
almam = na.omit(almam) #https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
plottings_sf(almam,date,sick,Group)



reduced2=function(u3,Group,name,lkm) {
  # u3=all_all1
  c1=c()
  
  DV=-u3[,'z0.ci_l']/(abs(u3[,'z0.ci_u']-u3[,'z0.ci_l']))
  AV=-u3[,'d0.ci_l']/(abs(u3[,'d0.ci_u']-u3[,'d0.ci_l']))
  OK=DV-AV
  
  pmeda=(u3[,'ACME']+abs(min(u3[,'ACME'])))/(u3[,'ACME']+abs(min(u3[,'ACME']))+u3[,'ADE']+abs(min(u3[,'ADE'])))
  c1= u3 
  c1=c1[c1[,'ACME']>quantile(c1[,'ACME'],0.75) & c1[,'d0.p']<quantile(c1[,'d0.p'],0.25) & 
          (c1[,'ACME']-c1[,'ADE']) > quantile(c1[,'ACME']-c1[,'ADE'],0.2) & 
          OK>quantile(OK,0.20) & pmeda > quantile(pmeda,0.20), ] #& 
  c1=c1[rev(order(c1[,'ACME'])),];
  
  #putatatiivien (tai melkein niin) poisto
  # rownames(c1)
  # the_ones=grepl('ok', rownames(c1))
  # c1=c1[the_ones,]
  # rownames(c1)  <- gsub(" ok", "",  rownames(c1) ) #:)
  # the_ones=grepl('Y_model_with_forced_int_m_and_x_ok', rownames(c1)); c1=c1[!the_ones,];
  # rownames(c1)  <- gsub(" Y_model_with_forced_int_x_ok", "",  rownames(c1) );rownames(c1)  <- gsub(" Y_model_with_forced_int_m_ok", "",  rownames(c1) ) #:)
  
  #ei poistoa, mutta piilotus:
  # rownames(c1)  <- gsub(" ok", "",  rownames(c1) ); 
  # rownames(c1)  <- gsub(" Y_model_with_forced_int_x_ok", "",  rownames(c1) );
  # rownames(c1)  <- gsub(" Y_model_with_forced_int_m_ok", "",  rownames(c1) ) 
  # rownames(c1)  <- gsub(" putative_m", "",  rownames(c1) ); 
  rt2=c1 #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot) # name=paste(simss,'basic hypothesis',take)# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model# https://github.com/MarioniLab/miloR# https://www.nature.com/articles/s41467-023-40458-9/figures/4
  name=paste('Contaminants_Steroids_BAs_or_Lipids_sims',date) # rtot=rtot_2000_mrct # rtot=uh5
  
  hoi=c(); 
  # hoi=scan(text=names(rt2[,1]), what=" ")#rownames(rt2)#check
  hoi=scan(text=rownames(rt2), what=" ")#rownames(rt2)#
  hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE) #check ncol
  ug=names(table(hoi$V3)) #check
  
  # for (i in 1:length(ug)) {rownames(c1)  <- gsub(paste('',ug[i]), '' ,rownames(c1) )}
  
  
  c1=tryCatch({c1[1:lkm,]}, error = function(msg){return(c1)}) 
  write.xlsx(c1, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE)
  
  return(c1)}


plottings_sf=function(uh7ma,date,sick,Group) {
  rt2=alma #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot) # name=paste(simss,'basic hypothesis',take)# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model# https://github.com/MarioniLab/miloR# https://www.nature.com/articles/s41467-023-40458-9/figures/4
  name=paste('Contaminants_Steroids_BAs_or_Lipids_sims',date) # rtot=rtot_2000_mrct # rtot=uh5
  
  hoi=c(); 
  hoi=scan(text=rt2[,1], what=" ") #rownames(rt2)#
  # hoi=scan(text=rownames(rt2), what=" ") #rownames(rt2)#
  hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE) #check this...
  
  
  # rownames(c1)
  # the_ones=grepl('ok', rownames(c1))
  # c1=c1[the_ones,]
  # rownames(c1)  <- gsub(" ok", "",  rownames(c1) ) #:)
  # the_ones=grepl('Y_model_with_forced_int_m_and_x_ok', rownames(c1)); c1=c1[!the_ones,];
  # rownames(c1)  <- gsub(" Y_model_with_forced_int_x_ok", "",  rownames(c1) );rownames(c1)  <- gsub(" Y_model_with_forced_int_m_ok", "",  rownames(c1) ) #:)  
  
  colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender') ## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
  
  hoi[,'Steroids' ][hoi[,'Steroids' ]=='17aOH.P4']='17a-OHP4'
  hoi[,'Steroids' ]  <- gsub("\\.", "-",  hoi[,'Steroids' ] ) #:)
  hoi[,'Steroids' ][ hoi[,'Steroids' ]=='T-Epi-T']='T/Epi-T'
  
  df2 <- hoi %>%make_long('Contaminants','Steroids','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
  meda='Sankey plot of_KAIKKI'
  pdf(paste(meda,name,sick,Group,".pdf"), width = 20, height = 20,  pointsize = 18);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + 
          theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
  jpeg(paste(meda,name ,sick,Group,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 18, res=500);
  print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
          scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));
  dev.off()
}
# https://stackoverflow.com/questions/8370548/how-can-i-interrupt-a-running-code-in-r-with-a-keyboard-command


# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/All/ok paitsi pea")
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/Uusix maksassa_all/ok_mcer/") #check this if needed...

name='all1111'; Group='All'; lkm=30; sick='just allaaa'
all_sankey=reduced2(all_all,Group,name,lkm)
plottings_sf(all_sankey,date,sick,Group)

name='female11'; Group='Female'; lkm=30; sick='just all'
all_fem_sankey=reduced2(all_fem,Group,name,lkm)
plottings_sf(all_fem_sankey,date,sick,Group)

name='male11'; Group='Male'; lkm=30; sick='just all'
all_male_sankey=reduced2(all_male,Group,name,lkm)
plottings_sf(all_male_sankey,date,sick,Group)




# name='female'; Group='Female'; lkm=30; sick='just all'
# all_fem_sankey=reduced2(all_fem,Group,name,lkm)
# Read 1532 items
# plottings_sf(all_fem_sankey,date,sick,Group)
# Read 90 items
# null device 
# name='male'; Group='Male'; lkm=30; sick='just all'
# all_male_sankey=reduced2(all_male,Group,name,lkm)
# Read 1548 items
# plottings_sf(all_male_sankey,date,sick,Group)




#Directit
u3=all_Male; Group='male'
c1=c()
ACMEMedian=c();ACMEpval=c();ACMEVar=c()
ADEMedian=c();ADEpval=c();ADEVar=c()
PMEDMedian=c();PMEDpval=c();PMEDVar=c();
DV=-u3[,'z0.ci_l']/(abs(u3[,'z0.ci_u']-u3[,'z0.ci_l']))
AV=-u3[,'d0.ci_l']/(abs(u3[,'d0.ci_u']-u3[,'d0.ci_l']))
OK=AV-DV
pmeda=(u3[,'ACME']+abs(min(u3[,'ACME'])))/(u3[,'ACME']+abs(min(u3[,'ACME']))+u3[,'ADE']+abs(min(u3[,'ADE'])))
c1= u3 #[u3[,'ADE'] < ADEMedian  & DV<ADEVar,] #& u3[,'z0.p']<ADEpval
c1=c1[c1[,'ADE']>quantile(c1[,'ADE'],0.75) & c1[,'z0.p']<quantile(c1[,'z0.p'],0.25) & (c1[,'ADE']-c1[,'ACME']) > quantile(c1[,'ADE']-c1[,'ACME'],0.2) & 
        OK>quantile(OK,0.20) & pmeda < quantile(pmeda,0.5), ] #& 

c1=c1[rev(order(c1[,'ADE'])),];

#putatiivin kanssa parempi... indirectinä
# rt2=c1 #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot) # name=paste(simss,'basic hypothesis',take)# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model# https://github.com/MarioniLab/miloR# https://www.nature.com/articles/s41467-023-40458-9/figures/4
# name=paste('Contaminants_Steroids_BAs_or_Lipids_sims',date) # rtot=rtot_2000_mrct # rtot=uh5
# hoi=c();
# hoi=scan(text=names(rt2[,1]), what=" ")#rownames(rt2)#
# hoi=as.data.frame(matrix(hoi, ncol = 4,  byrow = TRUE), stringsAsFactors = FALSE)
# ug=names(table(hoi$V4))
# for (i in 1:length(ug)) {rownames(c1)  <- gsub(paste('',ug[i]), '' ,rownames(c1) )}

#muttei direktinä..., eli sitten tämä
# rownames(c1)
# the_ones=grepl('ok', rownames(c1))
# c1=c1[the_ones,]
# 
# rownames(c1)  <- gsub(" ok", "",  rownames(c1) ) #:)
# the_ones=grepl('Y_model_with_forced_int_m_and_x_ok', rownames(c1)); c1=c1[!the_ones,];
# rownames(c1)  <- gsub(" Y_model_with_forced_int_x_ok", "",  rownames(c1) );rownames(c1)  <- gsub(" Y_model_with_forced_int_m_ok", "",  rownames(c1) ) #:)
# rownames(c1)  <- gsub(" putative_m", "",  rownames(c1) );

c1=tryCatch({c1[1:lkm,]}, error = function(msg){return(c1)}) 

# direct_sankey=function(rtot, 
med.freq.cut=0.3; med='no';med.min=0; med.prop=0.4;#name,lkm) {
rt2=c1

# if (med.min=='All') {rt2=rtot} else if (med.min!='All')  {if (med=='Mediated') 
# {rt2=rtot[rtot[,'ACME']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']>med.prop,];rt2=rt2[rt2[,'d0.p']<med.freq.cut,];rt2=rt2[rt2[,'z0.p']>med.freq.cut,]} else 
# {rt2=rtot[rtot[,'ADE']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']<med.prop,];rt2=rt2[rt2[,'d0.p']>med.freq.cut,];rt2=rt2[rt2[,'z0.p']<med.freq.cut,]}}
# rt2=rt2[rev(order(rt2[,'ADE'])),]; 
# rt2=rt2[1:lkm,]

# hoi=c(); hoi=scan(text=rownames(rt2), what="")
# hoi=as.data.frame(matrix(hoi, ncol = 4,  byrow = TRUE), stringsAsFactors = FALSE)

hoi=scan(text=rownames(rt2), what=" ") #rownames(rt2)#
hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE) #check this...

# colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids','Desig.')#,'Gender')
colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender')
hoi=hoi[,c('Contaminants','Bile Acids or Lipids')]## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct
# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation
# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
df2 <- hoi %>%make_long('Contaminants','Bile Acids or Lipids') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
meda='Sankey plot of_directi ilman putatiiveinaa_male'
pdf(paste(meda,name,med.freq.cut,med.prop ,Group,".pdf"), width = 20, height = 20,  pointsize = 10);
print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) + geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 5.5, color = 1, fill = "white") +scale_fill_viridis_d() + theme_sankey(base_size = 30) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
jpeg(paste(meda,name,med.freq.cut,med.prop,Group,".jpg"), width = 8000, height = 8000, quality = 100,pointsize = 12, res=500);
print(ggplot(df2, aes(x = x,  next_x = next_x, node = node,  next_node = next_node,fill = factor(node),label = node)) +
        geom_sankey(flow.alpha = 0.5, node.color = 1) + geom_sankey_label(size = 3.5, color = 1, fill = "white") +
        scale_fill_viridis_d() + theme_sankey(base_size = 16) + theme(legend.position = "none")+theme(axis.title.x = element_blank()));dev.off()
# return(rt2)}


#Tehdään tässä välissä kuvan '4' esimerkkejä: ok, eli tää on tää 6.6.24 ok jei! :)
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/Uusix maksassa_all/ok_mcer/") #check this if needed...
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/All/ok paitsi pea")
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic_cova/the all/") #check this if needed...
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/Steatosis")

#tiedostoilla
files <- list.files(pattern="*.RData")
ldf <- lapply(files, load)
list_of_files <- list() #create empty list
for (i in files) {print(i); list_of_files[[i]] <- get(load(paste0("", i)))}  #add files to list position#https://www.reddit.com/r/Rlanguage/comments/nq773b/reading_multiple_rdata_files_into_a_list/
names(list_of_files) <- files #https://stackoverflow.com/questions/38643000/naming-list-elements-in-r
# all: " uh7ma.RData"  3 files 
all_all=list_of_files$`All all.RData` #e on healthy :)
# female:
all_all=list_of_files$`All female.RData`
# male:
all_all=list_of_files$`All male.RData`


# list_of_files <- list() #create empty list
# for (i in files) {print(i); list_of_files[[i]] <- get(load(paste0("", i)))}
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic_cova/the all/") #check this if needed...
library(readxl)
#... all: " uh7ma.RData"  3 files 
all_all=read_xlsx(path = "100basic_okcovaha All tikka5624 .xlsx") #e on healthy :) #_2 
all_all=as.data.frame(all_all); all_all=all_all[!is.na(all_all[,1]),];rownames(all_all)=all_all[,1]; all_all=all_all[,2:dim(all_all)[2]]; all_all=all_all[rev(order(all_all[,1])),]
all_all=all_all; #all_Female=all_all;all_Male=all_all, this is important: all_all=all_all[!is.na(all_all[,1]),]
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/Covat/cova_14124/") #check this if needed...
# library(readxl)
# all: " uh7ma.RData"  3 files 
# all_all=read_xlsx(path = "100basic_okcovah_healthy All tikka25424 .xlsx") #e on healthy :)
# all_all=as.data.frame(all_all)
# rownames(all_all)=all_all[,1]
# all_all=all_all[,2:dim(all_all)[2]]
# all_all=all_all[rev(order(all_all[,1])),]
all_all=all_all[all_all[,1]>0,]

# > names(list_of_files) 
# [1] "Steatosis healthy all.RData"    "Steatosis healthy female.RData" "Steatosis healthy male.RData"   "Steatosis sick all.RData"      
# [5] "Steatosis sick female.RData"    "Steatosis sick male.RData" 

#This is for many files/folders with also non-basic...
# for (hii in 1:6) {
# hii=3
# rt2=list_of_files[[hii]]
# hoi=c(); hoi=scan(text=rownames(rt2), what="")
# hoi=matrix(hoi, ncol = 4,  byrow = TRUE) #or 3
# colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids','Desig.')#,'Gender') ,'Desig.')
# # hoi=hoi[,c('Contaminants','Steroids','Bile Acids or Lipids','Desig.')]## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct
# # https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation
# # https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
# # df2 <- hoi %>%make_long('Contaminants','Steroids','Bile Acids or Lipids','Desig.') #('Gender','Contaminants','Steroids','Bile Acids or Lipids')
# # hoi[1:100,c(1,2)]
# hoi[,c(2)]  <- gsub("\\.", "-",  hoi[,c(2)]  ) #:)
# # length(table(hoi[1:585,c(2)])) #min mukaan; Mediator[!Mediator %in% names(table(hoi[1:3774,c(2)]))]; "17aOH-P4" was mising
# Mediator_ok=Mediator[Mediator %in% names(table(hoi[1:dim(hoi)[1],c(2)]))]
# 

# all_all

#check if it is basic or not...
# rt2=data.frame(all_all,row.names = 1) #klassinen aikakilleri: row.names = 1



rt2=all_all
rt2=rt2[rt2[,1]>0,]
hoi=c(); hoi=scan(text=rownames(rt2), what="")
hoi=matrix(hoi, ncol = 3,  byrow = TRUE) #or 3
colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids')#,'Gender') ,'Desig.')
hoi[,c(2)]  <- gsub("\\.", "-",  hoi[,c(2)]  ) #:)
hoi[,'Steroids' ][hoi[,'Steroids' ]=='17aOH-P4']='17a-OHP4'

# #oh...
# Mediator=colnames(tv_all)[9:28];
# Treatment=colnames(tv_all)[71:77];
# Outcome=colnames(tv_all)[c(29:51,78:90)]; 
# Treatment=Treatment[!Treatment %in% c('Perfluorodecyl.ethanoic.acid')]
# Outcome=Outcome[!Outcome %in% c('Total_TG','PFAS','Perfluorodecyl.ethanoic.acid')]
# 
# # Mediator[Mediator=='17aOH-P4']='17a-OHP4'
# # Mediator_ok
# Mediator_ok=Mediator[Mediator %in% names(table(hoi[1:dim(hoi)[1],c(2)]))]
# # all_all
# 
# # print(Mediator[!Mediator %in% names(table(hoi[1:dim(hoi)[1],c(2)]))])
# Mediator_ok=Outcome[Outcome %in% names(table(hoi[1:dim(hoi)[1],c(3)]))]
# # length(table(hoi[1:100,c(1)]))
# # tf=hoi[1:585,c(1:2)]
# 




houdees=function(hoi, Outcome, Mediator,switch) {
  # rt2=all_all
  indir=c(); dir=c(); ip=c();rn=c();rn2=c()
  #direct...
  
  if (switch==1) {
    Mediator_ok=Outcome[Outcome %in% names(table(hoi[1:dim(hoi)[1],c(3)]))]
    for (i in 1:7) {for (j in 1:length(Mediator_ok)) {
      indir=append(indir,rt2[which(hoi[,1]==Treatment[i] & hoi[,3]==Mediator_ok[j])[1],c(5)]) #or c(1) hoi 1
      ip=append(ip,rt2[which(hoi[,1]==Treatment[i] & hoi[,3]==Mediator_ok[j])[1],c(6)])
      rn=append(rn,hoi[,3][which(hoi[,1]==Treatment[i] & hoi[,3]==Mediator_ok[j])[1]]) #change this...
      rn2=append(rn2,hoi[,1][which(hoi[,1]==Treatment[i] & hoi[,3]==Mediator_ok[j])[1]])
    }}} else if (switch==0) {
      # indir:
      Mediator_ok=Mediator[Mediator %in% names(table(hoi[1:dim(hoi)[1],c(2)]))]
      for (i in 1:7) {for (j in 1:length(Mediator_ok)) {
        indir=append(indir,rt2[which(hoi[,1]==Treatment[i] & hoi[,2]==Mediator_ok[j])[1],c(1)]) #or c(1) hoi 1
        ip=append(ip,rt2[which(hoi[,1]==Treatment[i] & hoi[,2]==Mediator_ok[j])[1],c(2)])
        rn=append(rn,hoi[,2][which(hoi[,1]==Treatment[i] & hoi[,2]==Mediator_ok[j])[1]]) #change this...
        rn2=append(rn2,hoi[,1][which(hoi[,1]==Treatment[i] & hoi[,2]==Mediator_ok[j])[1]])
      }}
    } 
  
  # rn2=rep(rn,times=2,each=2,len=length(indir))
  tot=cbind(rn2,rn,indir) #or indir or dir
  tot=tot[!is.na(tot[,1]),]
  tot=as.data.frame(tot)# tot %>%make_short('rn2','rn','indir')
  library(reshape2)
  jops=dcast(tot, rn2~rn, value.var='indir')
  jops[is.na(jops)]=0
  rownames(jops)=jops[,1]
  jops=jops[,2:dim(jops)[2]]
  jops=as.data.frame(jops)
  jopsr=matrix(as.numeric(unlist(jops)),nrow=dim(jops)[1],ncol=dim(jops)[2])
  colnames(jopsr)=colnames(jops);rownames(jopsr)=rownames(jops)
  # jopsr=jopsr[,Outcome] #groups[,'Abbreviation'][groups[,'Abbreviation'] %in% colnames(jopsr)]
  # jopsr=cbind(jopsr,matrix(rep(jopsr[,1]*0,9),ncol=9))
  # groups=groups[groups[,'Abbreviation']!='F',]
  # colnames(jopsr)[12:20]=groups[,'Abbreviation'][!groups[,'Abbreviation'] %in% colnames(jopsr)]
  
  #colnames(jopsr); muutoksia tähän
  
  # groups[,'Abbreviation']
  
  if (switch==1) {
    #for direct:
    jopsr=jopsr[,Outcome[Outcome %in% colnames(jopsr)]] #colnames(jopsr)][colnames(jopsr) %in% Outcome]
  } else if (switch==0) {
    #for indirect:
    ums=groups[order(groups[,'Group']),'Abbreviation']
    # ums=ums[ums!="17a-OHP4"]
    # jopsr[,'Abbreviation']
    # jopsr=jopsr[,groups[,'Abbreviation'][groups[,'Abbreviation'] %in% colnames(jopsr)]] #groups[,'Abbreviation'][groups[,'Abbreviation'] %in% colnames(jopsr)]
    jopsr=jopsr[,ums[ums %in% colnames(jopsr)]] 
  }
  # print(jopsr[rev(order(as.numeric(jopsr)))]) 
  
  # jopsr[jopsr>2]=2
  # jopsr[jopsr<-0.5]=-0.5
  
  # if (hii==3) {jopsr[rev(order(as.numeric(jopsr)))][1:8]=c(5.7,5.4,5,4.6,4.2,3.8,3.2,3)}
  # if (hii==3) {jopsr[rev(order(as.numeric(jopsr)))][1:6]=c(5,4.6,4.2,3.8,3.2,3)}
  # if (hii==3) {jopsr[rev(order(as.numeric(jopsr)))][1:3]=c(2,1.5,1)}
  
  # jopsr=rango(jopsr,0,0.5);
  
  tot=cbind(rn2,rn,ip)
  tot=tot[!is.na(tot[,1]),]
  tot=as.data.frame(tot)# tot %>%make_short('rn2','rn','indir')
  library(reshape2)
  jopsa=dcast(tot, rn2~rn, value.var='ip')
  jopsa[is.na(jopsa)]=0
  rownames(jopsa)=jopsa[,1]
  jopsa=jopsa[,2:dim(jopsa)[2]]
  jopsra=matrix(as.numeric(unlist(jopsa)),nrow=dim(jopsa)[1],ncol=dim(jopsa)[2])
  colnames(jopsra)=colnames(jopsa);rownames(jopsra)=rownames(jopsa)
  # jopsra=matrix(p.adjust(jopsra,method="BH"),nrow=dim(jopsra)[1],ncol=dim(jopsra)[2]); #bonferroni
  colnames(jopsra)=colnames(jopsa);rownames(jopsra)=rownames(jopsa)
  # jopsra=jopsra[,Outcome]#groups[,'Abbreviation'][groups[,'Abbreviation'] %in% colnames(jopsra)]
  
  
  if (switch==1) {
    #for direct:
    # jopsra=jopsra[,Outcome] 
    jopsra=jopsra[,Outcome[Outcome %in% colnames(jopsra)]]
  } else if (switch==0) {
    #for indirect
    jopsra=jopsra[,groups[,'Abbreviation'][groups[,'Abbreviation'] %in% colnames(jopsra)]]#groups[,'Abbreviation'][groups[,'Abbreviation'] %in% colnames(jopsra)]
    jopsra=jopsra[,ums[ums %in% colnames(jopsr)]] 
  }
  
  
  setwd("C:/Users/patati/Desktop/TurkuOW/RWork/") #check this if needed...
  hip1='transpose';width = 8000;height=2800;pch.cex=1.2;
  ho=paste('PFAS vs. steroids_ for the hypo4_colors_stea', switch)
  resulta1=jopsr
  p.mat.a1=jopsra
  # tv_ah=rango(resulta1,-0.5,0.5); resulta1=tv_ah;#
  #https://www.rdocumentation.org/packages/corrplot/versions/0.92/topics/corrplot
  #https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
  #https://statisticsglobe.com/change-font-size-corrplot-r
  #order can be: alphabet, hclust, original #https://stackoverflow.com/questions/51115495/how-to-keep-order-of-the-correlation-plot-labels-as-same-in-the-datafile
  order="original"; range='orig';corre='no_renorm'; type='full'; method='color';ga='All';gf='Female';gm='Male' #color square
  cl.offset=1.0;cl.length=5;cl.cex = 1.3;pch.cex=3.2;pch=20;cl.pos = 'r';
  
  
  if (switch==1) {rbo=rev(COL2('RdBu')[25:(length(COL2('RdBu'))-25)]) } else if (switch==0) {rbo=rev(COL2('RdBu')[25:100]) } 
  #rainbow(30)[1:24];rev(heat.colors(100))#cl.pos = 'b' ;#pch.cex=0.95,1.3; height=6300; pos 'b' cl.pos = 'b' 
  
  jpeg(paste("Square Correlation Plot of_7624_3_korkeatACME",ho,ga,hip1,".jpg"), width = width, height = height, quality = 100,pointsize = 30, res=300);# par( ps=ps)# par(cex.lab=90)
  col = brewer.pal(n = 9, name = "YlOrRd")
  corrplot(resulta1, type = type, order = order,method=method, p.mat=p.mat.a1, tl.col = "black", #sum(COL2('RdBu')=="#FF7417")
           cl.cex = cl.cex, pch.cex=pch.cex, pch.col='black',pch=pch,#pitikö vain pch lisätä pch väriin väriin... mystistä...'#FEE12B'
           sig.level = c(.05),cl.pos = cl.pos, insig = "label_sig", cl.offset=cl.offset,cl.length=cl.length,
           tl.srt = 90, diag = TRUE,col=rbo,is.corr = FALSE) #only in age...0.001, #rev(COL2('RdBu')[25:(length(COL2('RdBu'))-25)])
  dev.off()
  return(resulta1)
}
# }

uliulie=houdees(hoi, Outcome, Mediator,switch=0)


# install.packages("wesanderson")
# library(wesanderson) #http://www.sthda.com/english/wiki/colors-in-r

#Load all the variables in the folder:

setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/Tiedostot/") #check this if needed...

files <- list.files(pattern="*.RData")
ldf <- lapply(files, load)
list_of_files <- list() #create empty list
#loop through the files
for (i in files) {print(i); list_of_files[[i]] <- get(load(paste0("", i)))}  #add files to list position
#https://www.reddit.com/r/Rlanguage/comments/nq773b/reading_multiple_rdata_files_into_a_list/
names(list_of_files) <- files #https://stackoverflow.com/questions/38643000/naming-list-elements-in-r

#Steatosis all: "Steatosis uh7ma.RData"  28 files 
# name=paste(simss,'hypo4_yes_sick'); sick='yes'
# Group='All'; uh7ma=loop_med_boot_myj(Group,name,date,simss,sick,sick_group,joo,ip);try({uh7ma},{uh7ma=data.frame(0)})
# test=''; take=''; name=paste(simss,'hypo4_no_not sick'); sick='no'
# Group='All';    uh7mae=loop_med_boot_myj(Group,name,date,simss,sick,sick_group,joo,ip);try({uh7mae},{uh7mae=data.frame(0)})
# list_of_files[[28]]
s_all=list_of_files$`MASLD sick all.RData` #e on healthy :)
#Steatosis female:
s_fem=list_of_files$`MASLD sick female.RData`
#Steatosis male:
s_male=list_of_files$`MASLD sick male.RData`
#Fibrosis all:
f_all=list_of_files$`MASLD healthy all.RData`
#Fibrosis female:
f_fem=list_of_files$`MASLD healthy female.RData`
#Fibrosis male:
f_male=list_of_files$`MASLD healthy male.RData`

head(s_all)

s_all=s_all[rev(order(s_all[,1])),]; s_all20=s_all[1:10,]
s_fem=s_fem[rev(order(s_fem[,1])),]; s_fem20=s_fem[1:10,]
s_male=s_male[rev(order(s_male[,1])),]; s_male20=s_male[1:10,]

f_all=f_all[rev(order(f_all[,1])),]; f_all20=f_all[1:10,]
f_fem=f_fem[rev(order(f_fem[,1])),]; f_fem20=f_fem[1:10,]
f_male=f_male[rev(order(f_male[,1])),]; f_male20=f_male[1:10,]

totees=rbind(s_all20,s_fem20,s_male20,f_all20,f_fem20,f_male20) #https://stackoverflow.com/questions/38643000/naming-list-elements-in-r
# totees2=totees
# tot3=totees
names <- rownames(totees); rownames(totees) <- NULL;totees <- cbind(names,totees);totees=data.frame(totees)
try({colnames(totees)=c('Mediation','ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
}, {print('col');next})
# try({totees=totees[rev(order(totees[,2])),]},{print('order');next}) #totees=totees[rev(order(totees[,1])),]
totees[, c(2:17)] <- sapply(totees[, c(2:17)], as.numeric)
fn='noni'
write.xlsx(totees, file = paste(fn,Group,'sick_alle_s',sick,name,date,'.xlsx'), append = FALSE, row.names = TRUE)
write.csv(f_all,'.csv')

#loop through the files
for (i in files) {print(i); list_of_files[[i]] <- get(load(paste0("", i)))}  #add files to list position
#https://www.reddit.com/r/Rlanguage/comments/nq773b/reading_multiple_rdata_files_into_a_list/

names(list_of_files) <- files #https://stackoverflow.com/questions/38643000/naming-list-elements-in-r

#let's take it with sick...
#Steatosis all: "Steatosis sick malea.RData"  28 files 
s_all=list_of_files$`Steatosis healthy all.RData`
#Steatosis female:
s_fem=list_of_files$`Steatosis healthy female.RData`
#Steatosis male:
s_male=list_of_files$`Steatosis healthy male.RData`
#Fibrosis all:
f_all=list_of_files$`Fibrosis healthy all.RData`
#Fibrosis female:
f_fem=list_of_files$`Fibrosis healthy female.RData`
#Fibrosis male:
f_male=list_of_files$`Fibrosis healthy male.RData`

#HOMAIR all: "HOMAIR healthy all.RData"   files 
h_all=list_of_files$`HOMAIR healthy all.RData`
#HOMAIR female:
h_fem=list_of_files$`HOMAIR healthy female.RData`
#HOMAIR male:
h_male=list_of_files$`HOMAIR healthy male.RData`
# #Menopause all:
# m_all=list_of_files$`Menopause healthy all.RData`
# #Menopause female:
# m_fem=list_of_files$`Menopause healthy female.RData`
# #Menopause male:
# m_male=list_of_files$`Menopause healthy male.RData`

#Necroinflammation  all: "Necroinflammation  healthy all.RData"  .. files 
n_all=list_of_files$`Necroinflammation healthy all.RData`
#Steatosis female:
n_fem=list_of_files$`Necroinflammation healthy female.RData`
#Steatosis male:
n_male=list_of_files$`Necroinflammation healthy male.RData`

head(s_all)

s_all=s_all[rev(order(s_all[,1])),]; s_all_20=s_all[1:10,]
s_fem=s_fem[rev(order(s_fem[,1])),]; s_fem_20=s_fem[1:10,]
s_male=s_male[rev(order(s_male[,1])),]; s_male_20=s_male[1:10,]
f_all=f_all[rev(order(f_all[,1])),]; f_all_20=f_all[1:10,]
f_fem=f_fem[rev(order(f_fem[,1])),]; f_fem_20=f_fem[1:10,]
f_male=f_male[rev(order(f_male[,1])),]; f_male_20=f_male[1:10,]

h_all=h_all[rev(order(h_all[,1])),]; h_all_20=h_all[1:10,]
h_fem=h_fem[rev(order(h_fem[,1])),]; h_fem_20=h_fem[1:10,]
h_male=h_male[rev(order(h_male[,1])),]; h_male_20=h_male[1:10,]
m_all=m_all[rev(order(m_all[,1])),]; m_all_20=m_all[1:10,]
m_fem=m_fem[rev(order(m_fem[,1])),]; m_fem_20=m_fem[1:10,]
m_male=m_male[rev(order(m_male[,1])),]; m_male_20=m_male[1:10,]

n_all=n_all[rev(order(n_all[,1])),]; n_all_20=n_all[1:10,]
n_fem=n_fem[rev(order(n_fem[,1])),]; n_fem_20=n_fem[1:10,]
n_male=n_male[rev(order(n_male[,1])),]; n_male_20=n_male[1:10,]




totees=rbind(rep('s_all_10', 16),s_all_20,rep('s_fem_10', 16),s_fem_20,rep('s_male_10', 16),s_male_20,
             rep('f_all_10', 16),f_all_20,rep('f_fem_10', 16),f_fem_20,rep('f_male_10', 16),f_male_20,
             rep('n_all_10', 16),n_all_20,rep('n_fem_10', 16),n_fem_20,rep('n_male_10', 16),n_male_20,
             rep('h_all_10', 16),h_all_20,rep('h_fem_10', 16),h_fem_20,rep('h_male_10', 16),h_male_20,
             rep('m_all_10', 16),m_all_20,rep('m_fem_10', 16),m_fem_20,rep('m_male_10', 16),m_male_20) #https://stackoverflow.com/questions/38643000/naming-list-elements-in-r
# totees2=totees
# tot3=totees
names <- rownames(totees); rownames(totees) <- NULL;totees <- cbind(names,totees);totees=data.frame(totees)
try({colnames(totees)=c('Mediation','ACME', 'd0.p', 'd0.ci_l','d0.ci_u','ADE', 'z0.p', 'z0.ci_l','z0.ci_u','Proportion Mediated', 'n1.p','n.ci_l','n1.ci_u','Total Effect','tau.p','tau.ci_l','tau.ci_u') #d0 and d1 are the same as.. 'd1', 'd1.p',
}, {print('col');next})
# try({totees=totees[rev(order(totees[,2])),]},{print('order');next}) #totees=totees[rev(order(totees[,1])),]
totees[, c(2:17)] <- sapply(totees[, c(2:17)], as.numeric)


fn='table s2 redo_v5eea'
write.xlsx(totees, file = paste(fn,Group,'healthy_afmeee',sick,name,date,'testaas_4.xlsx'), append = FALSE, row.names = TRUE)


totees_sick=totees
totees_healthy=totees

#f totees_sick[13:22,1]
#m totees_sick[24:33,1]

hoi=c(); hoi=scan(text=totees_sick[13:22,1], what="")
hoi=as.data.frame(matrix(hoi, ncol = 4,  byrow = TRUE), stringsAsFactors = FALSE)
colnames(hoi)=c('Contaminants','Steroids','Bile Acids or Lipids','Desig.')#,'Gender')
hoi_sick=hoi[,c('Contaminants','Steroids','Bile Acids or Lipids')]## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct

table(hoi_sick[,1])[rev(order(table(hoi_sick[,1])))];table(hoi_healthy[,1])[rev(order(table(hoi_healthy[,1])))]
table(hoi_sick[,2])[rev(order(table(hoi_sick[,2])))];
table(hoi_healthy[,2])[rev(order(table(hoi_healthy[,2])))]
table(hoi_sick[,3])[rev(order(table(hoi_sick[,3])))];table(hoi_healthy[,3])[rev(order(table(hoi_healthy[,3])))]


table(hoi_sick[,1])[rev(order(table(hoi_sick[,1])))];

#Changing to excels: this was because the rownames...
a=list.dirs(path = "C:/Users/patati/Desktop/TurkuOW/RWork/hypo4/", full.names = TRUE, recursive = FALSE)
for (j in 1:(length(a)))
{
  filenames <- list.files(a[j], pattern="*.csv", full.names=TRUE)
  for(i in filenames) {
    b <- read.csv(i)
    new_name <- sub('.csv', '.xlsx', i, fixed = TRUE)
    write.xlsx(b, new_name)
  }
}


# I want to check the steroids with highest ACMEs between healthy (all, all) and sick (all, all).
# I need to have the cutoffs

setwd("C:/Users/patati/Desktop/TurkuOW/RWork/hypo_basic/Tiedostot/") #check this if needed... hypo4

files <- list.files(pattern="*.RData")
ldf <- lapply(files, load)
list_of_files <- list() #create empty list
#loop through the files
for (i in files) {print(i); list_of_files[[i]] <- get(load(paste0("", i)))}  #add files to list position
#https://www.reddit.com/r/Rlanguage/comments/nq773b/reading_multiple_rdata_files_into_a_list/
names(list_of_files) <- files #https://stackoverflow.com/questions/38643000/naming-list-elements-in-r
names(list_of_files)
library(stringr)

the_combos=function(list_of_files,Group,cond) { #cond='' vastaa kaikkia
#General categories of females or males (without 'All, 'MASLD, and 'Menopause')
u <- names(list_of_files)
a <- Group #note the writing, yes with " " in between. " female"  or " male" or " all"
ie=grep(a,u,fixed=TRUE)
u2=u[ie]
del=c(grep("All",u2,fixed=TRUE),grep("MASLD",u2,fixed=TRUE),grep("Menopause",u2,fixed=TRUE))
yl=1:length(u2); lop=yl[!yl %in% del]
u=u2[lop] #general male

ie=grep(cond,u,fixed=TRUE); u=u[ie]

a <- "sick"; ie=grep(a,u,fixed=TRUE); u3=u[ie]# sick ones, male or females
u_sick=list_of_files[u3]#https://www.tutorialspoint.com/how-to-extract-strings-that-contains-a-particular-substring-in-an-r-vector

a='healthy';ie=grep(a,u,fixed=TRUE);u4=u[ie]
u_healthy=list_of_files[u4] # u_healthy=list_of_files[grep(a,u,fixed=TRUE)] 
# plot(as.numeric(u_sick[[1]][,1]))# plot(1:500,as.numeric(u_sick[[1]][1:500,1]))

tcross=function(u_sick) { 
  all_names = c(); i=1
  for (i in (1:length(u_sick))) { #length(u_sick) is 18
    us=u_sick[[i]]
    aux=c();
    if (dim(us)[1]>200) {aux=200} else {aux=dim(us)[1]}
    plot(1:aux,as.numeric(us[1:aux,1]))
    values=(1:(length(as.numeric(us[,1]))-1))
    coo=c(); z=0
    for (z in values) {coo=append(coo,abs(us[z,1]-us[(z+1),1]))}  
    pss=which(coo>quantile(coo,0.95));ro=round(length(coo)/3)# pss[pss<ro]#round(length(which(coo>quantile(coo,0.95)))/2)]
    dpp=diff(pss[pss<ro]) #https://stackoverflow.com/questions/13602170/how-do-i-find-the-difference-between-two-values-without-knowing-which-is-larger
    dpp_sort=sort(dpp,decreasing = TRUE)
    if (length(dpp_sort)<5) { for_comp=length(dpp)+1} else {
      if (sum(dpp_sort[1:5]>5)==5) {cf=dpp_sort[5];cff=which(dpp>cf)[1]} else {cf=dpp_sort[2];cff=which(dpp>cf)[1]}
      cff=which(dpp>(cf-1))[1]
      for_comp=pss[cff]; }
    
    if (aux<30) {for_comp  =max(which(as.vector(us[,1]>0)))}
    if (aux>30) {if (max(pss[pss<ro])<30) (for_comp=30)} #due the small amount of good ones that can be like 4...
    
    rt2=us[1:for_comp,]; j=0
    hoi=c();for (j in 1:dim(rt2)[1]) {hoi=append(hoi,scan(text=rownames(rt2)[j], what=""))}
    hoi=as.data.frame(matrix(hoi, ncol = 4,  byrow = TRUE), stringsAsFactors = FALSE)
    hoi[,2] <- gsub("\\.", "-",  hoi[,2] ) #:)
    xz=round(quantile(table(hoi[,2]),0.25)); if (xz<2) {xz=0} else {xz=xz} #let's start with 25% and if not ok go to like 5%
    names=c();names=names(table(hoi[,2])[table(hoi[,2])>xz]) 
    # print(all_names)
    all_names=append(names,all_names)
  }
  return(sort(table(all_names),decreasing = TRUE))  
}

# sort(table(all_names))# diff(sort(table(all_names)))
tc_sick=tcross(u_sick);tc_healthy=tcross(u_healthy) # table(the_cross)# hist(table(the_cross),breaks=10)

cae1=as.numeric(names(sort(table(tc_sick),decreasing = TRUE)))[1];cae2=as.numeric(names(sort(table(tc_sick),decreasing = TRUE)))[2]
if (max(tc_sick)==cae1 | max(tc_sick)==cae2)
cfn=cae2-1 #sometimes as with females no differnces, i.. tc_sick;rev(as.numeric(names(table(tc_sick))))[2]-1#
if (is.na(cfn)) {steroid_sick=names(tc_sick)} else {steroid_sick=names(tc_sick[tc_sick>(cfn)])}

cae1=as.numeric(names(sort(table(tc_healthy),decreasing = TRUE)))[1];cae2=as.numeric(names(sort(table(tc_healthy),decreasing = TRUE)))[2]
if (max(tc_healthy)==cae1 | max(tc_healthy)==cae2)
cfn=cae2-1 
if (is.na(cfn)) {steroid_healthy=names(tc_healthy)} else {steroid_healthy=names(tc_healthy[tc_healthy>(cfn)])}
# steroid_healthy=names(tc_healthy[tc_healthy>cfn])
# https://stackoverflow.com/questions/45271448/r-finding-intersection-between-two-vectors
tbe=intersect(steroid_healthy, steroid_sick)

totaali_sh_all=steroid_sick[!steroid_sick %in% tbe] #https://www.geeksforgeeks.org/difference-between-two-vectors-in-r/ 

return(list(steroid_sick,tbe,totaali_sh_all)) } #"17a-OHP4" "DHT"      "DOC"      'P4'


names(list_of_files)


Group = ' all'; cond=''
all_all=the_combos(list_of_files,Group,cond) #ekat allit ('All...') oli deletoitu, yes, ja käytetty vain spesifisiä alleja...

Group = ' female'; cond=''
female_all=the_combos(list_of_files,Group,cond) 

Group = ' male'; cond=''
male_all=the_combos(list_of_files,Group,cond) 


Group = ' all'; cond='Steatosis'
all_steatosis=the_combos(list_of_files,Group,cond) 

Group = ' female'; cond='Steatosis'
female_steatosis=the_combos(list_of_files,Group,cond) 

Group = ' male'; cond='Steatosis'
male_steatosis=the_combos(list_of_files,Group,cond) 


Group = ' all'; cond='Fibrosis'
all_Fibrosis=the_combos(list_of_files,Group,cond) 

Group = ' female'; cond='Fibrosis'
female_Fibrosis=the_combos(list_of_files,Group,cond) 

Group = ' male'; cond='Fibrosis'
male_Fibrosis=the_combos(list_of_files,Group,cond) 


Group = ' all'; cond='Necroinflammation'
all_Necroinflammation=the_combos(list_of_files,Group,cond) 

Group = ' female'; cond='Necroinflammation'
female_Necroinflammation=the_combos(list_of_files,Group,cond) 

Group = ' male'; cond='Necroinflammation'
male_Necroinflammation=the_combos(list_of_files,Group,cond) 


Group = ' all'; cond='HOMAIR'
all_HOMAIR=the_combos(list_of_files,Group,cond) 

Group = ' female'; cond='HOMAIR'
female_HOMAIR=the_combos(list_of_files,Group,cond) 

Group = ' male'; cond='HOMAIR'
male_HOMAIR=the_combos(list_of_files,Group,cond) 



pottees=c(  all_all[3],female_all[3],male_all[3],
            all_steatosis[3],female_steatosis[3],male_steatosis[3],
            all_Fibrosis[3],female_Fibrosis[3],male_Fibrosis[3],
            all_Necroinflammation[3],female_Necroinflammation[3],male_Necroinflammation[3],
            all_HOMAIR[3],female_HOMAIR[3],male_HOMAIR[3])

pt1=pottees
pt2=pottees

mylist <- pt2

file <- "myfile_ok13824.txt"
conn <- file(description=file, open="w")

newlist <- lapply(seq_len(length(mylist)), function(i){
  
  lapply(seq_len(length(mylist[[i]])), function(j) {
    temp <- c(i, j, mylist[[i]][[j]])
    writeLines(text=paste(temp, collapse=","), con=conn, sep="\r\n")  
  })
  
})
close(conn)

pottees=c(  all_all[3],female_all[3],male_all[3],
            all_steatosis[3],female_steatosis[3],male_steatosis[3],
            all_Fibrosis[3],female_Fibrosis[3],male_Fibrosis[3],
            all_Necroinflammation[3],female_Necroinflammation[3],male_Necroinflammation[3],
            all_HOMAIR[3],female_HOMAIR[3],male_HOMAIR[3])

pottees_fem=c( 
           female_Fibrosis[3],
            female_Necroinflammation[3],
            female_HOMAIR[3])

pottees_male=c( 
  male_Fibrosis[3],
  male_Necroinflammation[3],
  male_HOMAIR[3])

pottees=c( 
            female_steatosis[3],male_steatosis[3],
            female_Fibrosis[3],male_Fibrosis[3],
            female_Necroinflammation[3],male_HOMAIR[3],
            female_HOMAIR[3],male_HOMAIR[3])


# unlist(pottees)[rev(order(unlist(pottees))]

# table(unlist(pottees))[rev(order(table(unlist(pottees))))]

table(unlist(pottees))[rev(order(table(unlist(pottees))))]

# earlier 
# all: DOC
# female"17a-OHP4" "DHT"      "DOC"      "E"        "E2" or
# male: #ei eroa..., joten ehkä vaan nämä: 
# "E2"      "T/Epi-T" "11-KDHT" "E"   ja verrattuna naisiin E ja E2 

# totaali_sh=steroid_sick[!steroid_sick %in% tbe] #https://www.geeksforgeeks.org/difference-between-two-vectors-in-r/ #DOC
#DOC is an important mediator in MASLD
#this was all, and females and males..., separately maybe needed...
totaali_sh_all=steroid_sick[!steroid_sick %in% tbe] #https://www.geeksforgeeks.org/difference-between-two-vectors-in-r/ 
totaali_sh_all                                         #"17a-OHP4" "DHT"      "DOC"      'P4'
fem_sh=steroid_sick[!steroid_sick %in% tbe] #naisilla:  "17a-OHP4" "DHT"      "DOC"      "E1"       "11-KT"    "11b-OHA4"   "P4"         
fem_sh
mal_sh=steroid_sick[!steroid_sick %in% tbe] #           "17a-OHP4" "DOC"      "DHT"                 "11-KT"    "11b-OHA4"   "P4"       "S"     
mal_sh

cfn=as.numeric(names(sort(table(tc_healthy),decreasing = TRUE)[2]))-1# rev(as.numeric(names(table(tc_healthy))))[2]-1
if (is.na(cfn)) {steroid_healthy=names(tc_healthy)} else {steroid_healthy=names(tc_healthy[tc_healthy>(cfn)])}


#Testit...
# values=(1:(length(as.numeric(u_sick[[1]][,1]))-1))
# coo=c()
# for (i in values) {coo=append(coo,abs(u_sick[[1]][i,1]-u_sick[[1]][(i+1),1]))}  
# pss=which(coo>quantile(coo,0.95));ro=round(length(coo)/3)# pss[pss<ro]#round(length(which(coo>quantile(coo,0.95)))/2)]
# dpp=diff(pss[pss<ro]) #https://stackoverflow.com/questions/13602170/how-do-i-find-the-difference-between-two-values-without-knowing-which-is-larger
# dpp_sort=sort(dpp,decreasing = TRUE)
# if (sum(dpp_sort[1:5]>5)==5) {cf=dpp_sort[5];cff=which(dpp>cf)[1]} else {cf=dpp_sort[2];cff=which(dpp>cf)[1]}
# 
# cff=which(dpp>cf)[1]
# for_comp=pss[cff]
# u_sick[[1]][1:for_comp,]



#https://stat.ethz.ch/R-manual/R-devel/library/base/html/sort.html

# dfddfdf=reduced2(list_of_files[[1]],name,lkm)

#For the top mediated:
BiocManager::install("colMap")
library(topGO)
library(ALL)
library(hgu95av2.db)
data(ALL)
data(geneList)


hoee=read.csv('topmed2_t14324.csv',header = FALSE, sep=';')
hoee=unique(hoee)
rownames(hoee)=1:17

# write.csv(hoee,'topmed2_t14324.csv')

hoee2=hoee[,2]
names(hoee2)=hoee[,1]
hoee2=1/(2.5*hoee2)

library(OmnipathR)
#fyi: https://pkgdown.r-lib.org/

d <- data.frame(protein_name = names(hoee2))
d <- translate_ids(d, protein_name, genesymbol)

sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = d, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)

sampleGOdata
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultFisher
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
pValue.classic <- score(resultKS)
pValue.elim    <- score(resultKS.elim)[names(pValue.classic)]
gstat          <- termStat(sampleGOdata, names(pValue.classic))
gSize          <- gstat$Annotated / max(gstat$Annotated) * 4

colMap <- function(x) {.col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
return(.col[match(1:length(x), order(x))]) } #https://rpubs.com/aemoore62/TopGo_colMap_Func_Troubleshoot

gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",pch = 19, cex = gSize, col = gCol)

sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]

cbind(termStat(sampleGOdata, sel.go),
      elim = pValue.elim[sel.go],
      classic = pValue.classic[sel.go])

showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')

#for the steroids:1
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic") #tests6.... 30.11.2023
t.val='minmax';simss=100; test=''; take=''; sick='no'; sick_group=sick_group

name='100basic'; tv_all=tv_covscl
Group='All';uh7ma=loop_med_simplified1a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Female'; uh7f=loop_med_simplified1a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Male';   uh7m=loop_med_simplified1a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';

# save.image("100basic_3624.RData")
# # setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic") #tests6.... 30.11.2023
# load("100basic_141223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

#...



#for the steroids:2
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic_cova") #tests.... 30.11.2023
t.val='minmax';
ccova=tv[,c("Steatosis.Grade.0.To.3", "Fibrosis.Stage.0.to.4" ,"Necroinflammation","HOMA-IR")] 
name='100basic_okcovaha'
sick_group=rowSums(ccova)>4 #toth#
simss=100; test=''; take='';  sick_group=sick_group
t.val='minmax';simss=100; test=''; take=''; sick='no'; sick_group=sick_group
# name='100basic'; 
tv_all=tv_covscl

sick='yes'
Group='All';    uh7ma=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Female'; uh7f=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Male';   uh7m=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);

sick='no'
Group='All';    uh7mae=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Female'; uh7fe=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Male';   uh7me=loop_med_simplified2b(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);

save.image("100basic_cova_14424.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic_cova") #tests6.... 30.11.2023
# load("100basic_cova_141223.RData")
lkm=30;Group='All'; name='ok3';date='tikka14424'
alma=reduced2(uh7ma,Group,name,lkm);
plottings_sf(alma,date,sick,Group)

date='tikka14424_m2';name='ok2';Group='male';
almam=reduced2(uh7f,Group,name,lkm);
almam = na.omit(almam)
plottings_sf(almam,date,sick,Group)

date='tikka14424_f3_s';name='ok2';Group='female';
almaf=reduced2(uh7m,Group,name,lkm);
almaf = na.omit(almaf) #https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
plottings_sf(almaf,date,sick,Group)

lkm=30;Group='All'; name='_not sick';date='tikka14424'
alma=reduced2(uh7mae,Group,name,lkm);
plottings_sf(alma,date,sick,Group)

date='tikka14424_m';name='_not sick';Group='male';
almam=reduced2(uh7fe,Group,name,lkm);
almam = na.omit(almam)
plottings_sf(almam,date,sick,Group)

date='tikka14424_f';name='ok2_not sick';Group='female';
almaf=reduced2(uh7me,Group,name,lkm);
almaf = na.omit(almaf) #https://www.tutorialspoint.com/how-to-remove-rows-from-data-frame-in-r-that-contains-nan
plottings_sf(almaf,date,sick,Group)

plottings(uh7ma,uh7m,uh7f,test,take,date)
plottings(uh7mae,uh7me,uh7fe,test,take,date)
#for the steroids:3... not ok/excessive
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo1") #tests6.... 30.11.2023
# Group='All';uh7ma=loop_med_simplified3(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo1") #tests6.... 30.11.2023
# load("100basic_cova_191223.RData")
# plottings(uh7ma,uh7m,uh7f,t.val,test,take,date)

#Mediation explanation:
data(jobs)
####################################################
# Example 1: Linear Outcome and Mediator Models
####################################################
b <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
c <- lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data=jobs)
# Estimation via quasi-Bayesian approximation
contcont <- mediate(b, c, sims=50, treat="treat", mediator="job_seek")

summary(contcont)

plot(contcont)

mean(as.vector(unlist(summary(contcont)[7])))

heps=as.vector(unlist(summary(contcont)[7]))

wilcox.test(heps)
t.test(heps, mu = 100)


contcont.boot <- mediate(b, c, boot=TRUE, sims=50, treat="treat", mediator="job_seek")
summary(contcont.boot)


t.test(rep(1,50),heps)
summary(contcont)



# https://rdrr.io/cran/mediation/src/R/mediate.R

# summary(med_oute)
# 
# Causal Mediation Analysis 
# 
# Quasi-Bayesian Confidence Intervals
# 
# Estimate 95% CI Lower 95% CI Upper p-value
# ACME            -0.3667      -1.0206         0.25    0.26
# ADE              0.2777      -1.2545         1.78    0.70
# Total Effect    -0.0891      -1.3952         1.27    0.92
# Prop. Mediated   0.0103      -9.0591         6.51    0.98
# 
# Sample Size Used: 35 


# Simulations: 100 

# pval(med_oute$d0.sims,med_oute$d0)
# [1] 0.26 #ok :)



#no ni...
# https://stackoverflow.com/questions/64829727/mediation-r-package-p-values-workaround-to-get-more-significant-digits


# contcont.boot
# 
# require(parallel)
# require(MASS)
# data(jobs)
# a <- lm(comply ~ treat + sex + age + marital + nonwhite + educ + income,
#         data = jobs)
# b <- glm(job_dich ~ comply + treat + sex + age + marital + nonwhite + educ + income,
#          data = jobs, family = binomial)
# # 10 jobs
# c <- lm(depress2 ~ job_dich * (comply + treat) + sex + age + marital + nonwhite + educ + income,
#         data = jobs)
# out <- ivmediate(a, b, c, sims = 50, boot = FALSE, enc = "treat", treat = "comply", mediator = "job_dich")
# 
# out
# 
# wilcox.test(out$dc0.sims,out$dc1.sims)
# 
# 
# x = cumsum(c(0, runif(100, -1, +1)))
# y = cumsum(c(0, runif(100, -1, +1)))
# fit = lm(y ~ x)
# summary(fit)
# 
# lmp <- function (modelobject) {
#   if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
#   f <- summary(modelobject)$fstatistic
#   p <- pf(f[1],f[2],f[3],lower.tail=F)
#   attributes(p) <- NULL
#   return(p)
# }
# 
# lmp(fit)
# 
# summary(out)
# plot(out)
# lmp(out)





#for the steroids:'4' (eka uusiks)
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo1_ok") #tests6.... 30.11.2023
t.val='no';simss=100; test=''; name='hypo1_ok'
Group='All';uh7ma=loop_med_simplified4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Female'; uh7f=loop_med_simplified4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Male';   uh7m=loop_med_simplified4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100hypo1_cova_ok_141223.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo1_ok") #tests6.... 30.11.2023
# load("100hypo1_cova_ok_141223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

#for the steroids:'5' (kakkonen)
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/t7_lasso/") #tests6.... 30.11.2023
t.val='no';simss=10; test=''; take=''; 
# Outcome=colnames(tv_all)[c(29:70)]; Treatment=colnames(tv_all)[82:78]; Mediator=colnames(tv_all)[9:28]; X <- tv_all[,Treatment]; 
name='10hypo3_lassoana2'
Group='All';uh7ma=loop_med_simplified5ö(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
# write.xlsx(uh7ma, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
# name='100hypo2_aiced2a'
Group='Female'; uh7f=loop_med_simplified5ö(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
# write.xlsx(uh7f, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
# name='100hypo2_aiced2b'
Group='Male';   uh7m=loop_med_simplified5ö(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
# write.xlsx(uh7m, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100hypo3_lassoa_191223.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo2") #tests6.... 30.11.2023
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/t7_lasso/") #tests6.... 30.11.2023
# load("100hypo2_cova_ok_aiced2_131223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

#3rd plots... for elucidating the important factors for visualizations and excel files (esp. related to mediation)
# rtot=uh7ma
# p.med<- rtot[,'Proportion Mediated']
# ACME<- rtot[,'ACME']
# ACME_pval<- rtot[,'d0.p']
# scatter3d(x = p.med, y = ACME, z = ACME_pval)

#for the steroids:'5' (kakkonen)
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo2_aiced") #tests6.... 30.11.2023
t.val='no';simss=100; test=''; 
name='100hypo2_aiced'
Group='All';uh7ma=loop_med_simplified5(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
# name='100hypo2_aiced2a'
Group='Female'; uh7f=loop_med_simplified5(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
# name='100hypo2_aiced2b'
Group='Male';   uh7m=loop_med_simplified5(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100hypo2_cova_ok_aiced_131223.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo2") #tests6.... 30.11.2023
# load("100hypo2_cova_ok_141223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

#for the steroids:'5' (kakkonen)
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo4_aiced") #tests6.... 30.11.2023
name='100hypo4_aiced2e'
Group='All';
uh7ma=loop_med_simplified5ö2(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
# write.xlsx(uh7ma, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
# name='100hypo2_aiced2a'
Group='Female'; 
uh7f=loop_med_simplified5ö2(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
# write.xlsx(uh7f, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
# name='100hypo2_aiced2b'
Group='Male';   
uh7m=loop_med_simplified5ö2(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
# write.xlsx(uh7m, file = paste(name,Group,date,'.xlsx'), append = FALSE, row.names = TRUE) 
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100hypo4_cova_ok_aiced_131223.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo2") #tests6.... 30.11.2023
# load("100hypo2_cova_ok_141223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

#for the steroids:'5' (kakkonen)
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests8/") #https://topepo.github.io/caret/parallel-processing.html
# https://ds-pl-r-book.netlify.app/optimization-in-r.html
t.val='no';simss=100; test=''; take=''; 
name='100hypo3_healthyee'; sick='no'
Group='All';uh7man=loop_med_simplified5ö4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,toth);
Group='Female'; uh7fn=loop_med_simplified5ö4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,toth);
Group='Male';   uh7mn=loop_med_simplified5ö4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,toth);
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100ehhypo3ee_sick_8124.RData")

#for the steroids:'5' (kakkonen) ILMAN adjustementtia
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo2_ilman") #tests6.... 30.11.2023
t.val='no';simss=100; test=''; 
name='100hypo2_ilman'
Group='All';uh7ma=loop_med_simplified5a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Female'; uh7f=loop_med_simplified5a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Male';   uh7m=loop_med_simplified5a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100hypo2_ilmancova_ok_141223.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo2") #tests6.... 30.11.2023
# load("100hypo2_ilmancova_ok_141223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

#for the steroids:'6' (kolmonen) #ämmä vapaana eli siinä kaikki variaabelit
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo3") #tests6.... 30.11.2023
t.val='no';simss=100; test=''; take=''; 
name='100hypo3'
Group='All';uh7ma=loop_med_simplified6(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Female'; uh7f=loop_med_simplified6(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Male';   uh7m=loop_med_simplified6(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100hypo3_cova_ok_141223.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo3") #tests6.... 30.11.2023
# load("100hypo3_cova_ok_141223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

#for the steroids:'6' (kolmonen) #ämmä vapaana eli siinä kaikki variaabelit ja ilman adj
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo3_ilman") #tests6.... 30.11.2023
t.val='no';simss=100; test=''; take=''; 
name='100hypo3_ilman'
Group='All';uh7ma=loop_med_simplified6a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Female'; uh7f=loop_med_simplified6a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Male';   uh7m=loop_med_simplified6a(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100hypo3_ilmancova_ok_141223.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo3") #tests6.... 30.11.2023
# load("100hypo3_cova_ok_141223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

#for the steroids:'7' (ykis) #ämmä vapaana eli siinä kaikki variaabelit
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo1_ilman") #tests6.... 30.11.2023
t.val='no';simss=100; test=''; take=''; 
name='hypo1_ilman'
Group='All';uh7ma=loop_med_simplified7(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Female'; uh7f=loop_med_simplified7(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Male';   uh7m=loop_med_simplified7(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100hypo1_ilmancova_ok_141223.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo1_ilman") #tests6.... 30.11.2023
# load("100hypo1_ilmancova_ok_141223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

#for the steroids:'8' (ykis) #ämmä vapaana eli siinä kaikki variaabelit
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/basic_minmax") #tests6.... 30.11.2023
t.val='minmax';simss=100; test=''; take=''; 
name='100basic_minmax'
Group='All';uh7ma=loop_med_simplified7(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Female'; uh7f=loop_med_simplified7(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Male';   uh7m=loop_med_simplified7(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100basic_minmax_141223.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/basic_minmax") #tests6.... 30.11.2023
# load("100basic_minmax_141223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

#for the steroids:'...' (kakkonen) JA MINMAX
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo2") #tests6.... 30.11.2023
t.val='minmax';simss=100; test=''; 
name='100hypo2'
Group='All'; uh7ma=loop_med_simplified5ö3(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Female'; uh7f=loop_med_simplified5ö3(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
Group='Male';   uh7m=loop_med_simplified5ö3(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test);
uh7m=cbind(uh7m,'Male');colnames(uh7m)[dim(uh7m)[2]]='Gender';
uh7f=cbind(uh7f,'Female');colnames(uh7f)[dim(uh7f)[2]]='Gender';
save.image("100hypo2_cova_ok_141223.RData")
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/hypo2") #tests6.... 30.11.2023
# load("100hypo2_cova_ok_141223.RData")
plottings(uh7ma,uh7m,uh7f,test,take,date)

# save.image("100hypo2asdfasdf_cova_ok_141223.RData")

#Hypoteesi4:
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests11/") #https://topepo.github.io/caret/parallel-processing.html
# https://ds-pl-r-book.netlify.app/optimization-in-r.html
# t.val='no';simss=100; test=''; take=''; name='100hypo3_yes'; sick='yes'
# Group='All'
t.val='no'; simss=100; test=''; take=''; name='100hypo4_yes'; sick='yes'
Group='All'; uh7ma=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Female'; uh7f=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Male';   uh7m=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
t.val='no'; simss=100; test=''; take=''; name='100hypo4_no'; sick='no'
Group='All'; uh7mae=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Female'; uh7fe=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Male';   uh7me=loop_med_simplified_hyp4(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
save.image("100hypo4_basicaaaha_all_18124.RData")
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests11/")
# load("100hypo4_basicaaaha_all_10124.RData") 
uh7mae[1,1]; uh7ma[1,1]

#for the steroids:(hypoteesi3)
setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests11/") #https://topepo.github.io/caret/parallel-processing.html
# https://ds-pl-r-book.netlify.app/optimization-in-r.html

# date='tikka19124' #change this... :)
t.val='no';simss=100; test=''; take=''; name='100hypo3_yes'; sick='yes'
Group='All';    uh7ma=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Female'; uh7f=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Male';   uh7m=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
t.val='no';simss=100; test=''; take=''; name='100hypo3_no'; sick='no'
Group='All';    uh7mae=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Female'; uh7fe=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
Group='Male';   uh7me=loop_med_simplified5ö5(Treatment, Mediator, Outcome,tv_all,Group,name,simss,t.val,test,sick,sick_group);
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests11/")
# save.image("100hypo3_basice_all_19124.RData")
# load("100hypo3_basicah_all_18124.RData")

#Steroidien lukumäärän pienennys mediaatiota varten, testailua:  
head(res[rev(order(res[,'ACME'])),])
install.packages('pathviewr')
devtools::install_github("ahasverus/elbow")
library(pathviewr); library(nlstools);library(elbow)
df <- data.frame(x = seq(1:length(cv_model$lambda)),y = cv_model$lambda)
hei=elbow(data = df)
plot(hei$data[,4])
yn=which(hei$data[,4]==min(hei$data[,4]))[1]
df2=hei$data[yn:dim(hei$data)[1],4]
df2 <- data.frame(x = seq(1:length(df2)),y = df2)
yn2=find_curve_elbow(df2)
yti=yn+yn2
jeps=cv_model$lambda[yti]
plot(cv_model)
abline(v=log(jeps), col="blue") #http://www.sthda.com/english/wiki/abline-r-function-an-easy-way-to-add-straight-lines-to-a-plot-using-r-software
df <- data.frame(x = seq(1:length(cv_model$lambda)),y = cv_model$lambda)
hei=elbow(data = df)
yn=which(hei$data[,4]==min(hei$data[,4]))[1]
df2=hei$data[yn:dim(hei$data)[1],4]
df2 <- data.frame(x = seq(1:length(df2)),y = df2)
yn2=find_curve_elbow(df2)
yti=yn+yn2
jeps=cv_model$lambda[yti]
plot(cv_model)
abline(v=log(jeps), col="blue")

# Mediaation testailua
# å=length(y); list=1:å; list1=list[1:å!=i];ya=c();
# list2=list1+28; DataA=Data; DataA[,list1]=-DataA[,list1]
# # ya=y; ya[list1[1]]=paste('-',ya[list1[1]],sep='')
# # yt1=paste(ya[list1],collapse= "-")
# # dim(Data)[2]
# # Data <- cbind(X,M,Y,cova); 
# # list=28
# xm=paste(c(x,m), collapse= "+")
# all_right=paste(xm,yt1,sep='')
# fmla2 <- as.formula(paste(y[i]," ~ ", all_right))
# c = lm(fmla2, Data)

#Here is some tests for adjusting the variables with dummy codes....
# Install and load the caret package
# install.packages("caret")
# library(caret)
# # Create a sample dataset with categorical variables
# data <- data.frame(
#   Gender = c("Male", "Female", "Male", "Female"),
#   Region = c("North", "South", "East", "West"),
#   Age = c(25, 30, 22, 35),
#   Income = c(50000, 60000, 45000, 70000),
#   Outcome = c(1, 0, 1, 0))
# # Display the original dataset
# print("Original Dataset:")
# print(data)
# # Specify the indices of categorical variables
# cat_cols <- c("Gender", "Region")
# # Create a dummy variable transformer
# dummy_transformer <- dummyVars("~.", data = data, levelsOnly = TRUE, cat_cols = cat_cols)
# # Apply the transformer to the original dataset
# dummy_data <- predict(dummy_transformer, newdata = data)
# # Display the dataset with dummy variables
# print("Dataset with Dummy Variables:"); print(dummy_data)

# https://stats.stackexchange.com/questions/218401/testing-mediation-and-moderation-can-one-variable-function-as-both-mediator-and
# # Load the mediation package 
# library(mediation) 
# # Summarize the mediation analysis results # summary(mediation_model) 
# https://stackoverflow.com/questions/75802827/how-to-colour-skankey-node-with-all-colours-of-first-node-that-it-is-related-to
# p <- ggplot(df2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
# p +geom_rect(data = dat, aes(xmin = x, xmax = x + width,

# Import data
#Sankey with graduate change of color:
# library(Director)
# data(ovca)# Level1 <- createList(poorprog$Level1) # miRNA-gene data Level2 <- createList(poorprog$Level2) # gene-pathway data. List <- append2List(Level1, Level2) # append lists.
# Figure 2A with default settings
#initSankey() # Generates necessary javascript and CSS templates.
#template1 <- makeSankey(List, averagePath = TRUE)
#Sankey1 <- drawSankey(template1, caption="Top SPIA pathways with pGFdr < 0.1.<br><font size=-1>Nodes are red
#if upregulated in good diagnostic samples.</font>",width=1100, height=500,legendfont="lato, helvetica, sans-serif",legendsize=20, nodeValue="Fold change")
# remotes::install_github("ssp3nc3r/ggSankeyGrad")
# library(ggSankeyGrad)
# c1     <- c("A", "A", "B", "B")
# c2     <- c("C", "D", "C", "D")
# values <- c(2L, 5L, 8L, 3L)
# col1   <- c("red", "red", "green", "green")
# col2   <- c("blue", "orange", "blue", "orange")
# ggSankeyGrad(c1, c2, col1, col2, values, label = TRUE)
# # https://www.informationisbeautifulawards.com/showcase/204-nobels-no-degrees
# d5 <- read.csv(text = '
# "Category","University","n","col1","col2"
# "Chemistry","Berkeley",5,"#cc5b47","#003262"
# "Chemistry","Caltech",4,"#cc5b47","#FF6C0C"
# "Chemistry","Cambridge",3,"#cc5b47","#A3C1AD"
# "Chemistry","Columbia",3,"#cc5b47","#B9D9EB" #... color sankeys..
# rtot=uh7ma #[,1:17]# rtot=rtot[,1:17]# rtot=data.frame(rtot) # name=paste(simss,'basic hypothesis',t.val,take)# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model# https://github.com/MarioniLab/miloR# https://www.nature.com/articles/s41467-023-40458-9/figures/4
# name=paste('Contaminants_Steroids_BAs_or_Lipids_100sims_basic_all',t.val,test,take,date) # rtot=rtot_2000_mrct # rtot=uh5
# pdf(paste("ACME_p_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'d0.p']),breaks=100, xlab="P value of ACME", ylab="Frequency",main="");dev.off()
# pdf(paste("Mediated Proportion_histogram",name,".pdf"), width = 80, height = 80, pointsize = 70); hist(as.numeric(rtot[,'Proportion Mediated']),breaks=300,xlim=c(-1.2,1.2),xlab="Proportion of Mediated", ylab="Frequency",main="");dev.off()
# med.freq.cut=0.3; med.min=0; med.prop=0.1; Group='All'; med='Mediated' #or not mediated e.g. no
# if (med.min=='All') {rt2=rtot} else if (med.min!='All')  {if (med=='Mediated') {rt2=rtot[rtot[,'ACME']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']>med.prop,]} else {rt2=rtot[rtot[,'ADE']>med.min,];rt2=rt2[rt2[,'Proportion Mediated']<med.prop,]}}
# if (med=='Mediated') {rt2=rt2[rt2[,'d0.p']<med.freq.cut,]; rt2=rt2[rt2[,'z0.p']>med.freq.cut,];rt2=rt2[rt2[,'ACME']>med.min,]} else {rt2=rt2[rt2[,'d0.p']>med.freq.cut,]; rt2=rt2[rt2[,'z0.p']<med.freq.cut,];rt2=rt2[rt2[,'ADE']>med.min,]}
# hoi=c(); hoi=scan(text=rownames(rt2), what="")
# hoi=as.data.frame(matrix(hoi, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
# hoi[,4]=rt2[,1]; colnames(hoi)=c('Contaminants','Steroids','BAs_Lipids','Strength')#,'Gender') ## https://stats.stackexchange.com/questions/282155/causal-mediation-analysis-negative-indirect-and-total-effect-positive-direct# https://www.researchgate.net/post/How_can_I_interpret_a_negative_indirect_effect_for_significant_mediation# https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot replacing dot
# hoi[,'col1']=string.to.colors(hoi[,'Steroids'])
# hoi[,'col2']=string.to.colors(hoi[,'BAs_Lipids'])
# hoi=hoi[1:20,]; hoi[,4]=hoi[,4]*20; hoi[1:5,4]=hoi[1:5,4]*20
# with(hoi, ggSankeyGrad(c1 = Steroids,c2 = BAs_Lipids,col1 = col1[1:17], col2 = col2[1:17],values = Strength, label = TRUE))

#For vif/bif/other characteristics of the model, e.g. AUC
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/") #tests6.... 30.11.2023
#...
# # https://rdrr.io/cran/mlma/man/data.org.html # https://cran.r-project.org/web/packages/mma/mma.pdf
# control.value=colMins(as.matrix(X)) #test also with colMedians colMins -abs(colMins(as.matrix(X))*2
# treat.value=colMaxs(as.matrix(X))
# #M~X
# x=c();m=c(); y=c();ye=c()
# for (i in 1:length(colnames(X))) {x=append(x,paste("Data[, ",i , "]", sep=""))}
# for (j in (dim(X)[2]+1):(length(colnames(M))+dim(X)[2])) {m=append(m,paste("Data[, ",j , "]", sep="")) }
# #Y~X+M
# for (z in (dim(M)[2]+dim(X)[2]+1):(dim(Data)[2])) {y=append(y,paste("Data[, ",z , "]", sep="")) } #this dimension was essential for the loop names
# med_out=c();res=c(); tmp=c();rn=c();med_oute=c();med_sense=c();resa=c() 
# j=1;i=1;z=1
# simss=2; length(y)*length(m)*length(x)
# # https://fmch.bmj.com/content/8/1/e000262
# bvif=c()# cvif=c()
# for (i in 1:length(y)) {
#   # for (j in 1:length(m)) { #control.value=mina[i]
#     # for (z in 1:length(x)) {
#       fmla1 <- as.formula(paste(paste(m, collapse= "+")," ~ ", paste(x, collapse= "+"))) 
#       b = lm(fmla1, Data)
#       bvif = append(bvif,vif(b))  
#       xm=paste(paste(c(x,m), collapse= "+"))
#       fmla2 <- as.formula(paste(y[i]," ~ ", xm)) #https://www.statology.org/glm-vs-lm-in-r/
#       c = lm(fmla2, Data) 
#       cvif = append(cvif,vif(c)) }
# par(mfrow=c(2,2)); plot(c)
# bvif = append(bvif,vif(b))
# cvif = append(cvif,vif(c))
# hist(bvif); hist(cvif)
# library(GGally)
# ggpairs(data.frame(tv_all[,18:28]))
# #define multiple linear regression model
# model <- lm(rating ~ points + assists + rebounds, data=df)
# #calculate the VIF for each predictor variable in the model
# vif(model)      

# Lasso best model option 1:
# library(MASS); library(glmnet)
# Outcome=colnames(tv_all)[c(29:70,80:90)]; Treatment=colnames(tv_all)[71:78]; Mediator=colnames(tv_all)[9:28]; X <- tv_all[,Treatment]; X <- tv_all[,Treatment]; 
# lasso_fit <- glmnet(x = as.matrix(tv_all[,c(71:78)]), y = tv_all[,Outcome[1]], alpha = 1)
# #lambda = 0.5
# coef(lasso_fit,s=0.1)
# plot(lasso_fit, xvar = "lambda", label = TRUE)
# library(plotmo) # for plot_glmnet; install.packages('plotmo') https://www.rdocumentation.org/packages/plotmo/versions/3.6.2
# plot_glmnet(lasso_fit, label=TRUE) 
# plot_glmnet(lasso_fit, label=8, xvar ="norm")   
# #use 5-fold cross validation to pick lambda
# cv_lasso_fit <- cv.glmnet(x = as.matrix(tv_all[,c(71:78)]), y = tv_all[,Outcome[1]], alpha = 1, nfolds = 5)
# plot(cv_lasso_fit)
# cv_lasso_fit$lambda.min
# Boston_IS_pred <- predict(lasso_fit, as.matrix(tv_all[,c(71:78)]), s = cv_lasso_fit$lambda.min) #https://xiaorui.site/Data-Mining-R/lecture/3.C_LASSO.html
# df.comp <- data.frame(beta    = beta,Linear  = li.eq$coefficients,Lasso   = la.eq$beta[,1],Ridge   = ri.eq$beta[,1])
#perform k-fold cross-validation to find optimal lambda value
# cv_model <- cv.glmnet(X, y, alpha = 1)
# #find optimal lambda value that minimizes test MSE
# best_lambda <- cv_model$lambda.min
# best_lambda
# # [1] 5.616345
# #produce plot of test MSE by lambda value
# # plot(cv_model) 
# #find coefficients of best model
# best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)
# # coef(best_model)
# coef(best_model)[abs(coef(best_model))>0.0001]
# abs(coef(best_model))>0.0001 
# cv_model <- cv.glmnet(X, y, alpha = 1); best_lambda <- cv_model$lambda.min
# best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda) #coef(best_model)[abs(coef(best_model))>0.00001]
# abs(coef(best_model))>0.0001 

# https://mathstat.slu.edu/~speegle/Spring2020/4870/_book/variable-selection.html
# https://helda.helsinki.fi/server/api/core/bitstreams/e4842c4f-fd4f-4c2e-8226-aac8b03f5ac5/content
# https://bioconductor.org/packages/devel/bioc/vignettes/Director/inst/doc/vignette.pdf


# These tests are almost directly from previous works, see the links below: 
# #WQS test
# # install.packages('wqs')
# library(wqs)
# install.packages('gWQS')
# # https://cran.r-project.org/web/packages/gWQS/vignettes/gwqs-vignette.html
# library(gWQS)
# 
# library(ggplot2)
# library(knitr)
# library(kableExtra)
# library(reshape2)
# 
# # we save the names of the mixture variables in the variable "PCBs"
# PCBs <- names(wqs_data)[1:34]
# # we run the model and save the results in the variable "results2i"
# results2i <- gwqs(yLBX ~ pwqs + nwqs, mix_name = PCBs, data = wqs_data, 
#                   q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
#                   family = "gaussian", seed = 2016)
# 
# summary(results2i)
# 
# # bar plot
# gwqs_barplot(results2i)
# # scatter plot y vs wqs
# gwqs_scatterplot(results2i)
# # scatter plot residuals vs fitted values
# gwqs_fitted_vs_resid(results2i)
# # boxplot of the weights estimated at each repeated holdout step
# gwqs_boxplot(results2i)
# 
# head(results2i$final_weights)
# 
# gwqs_summary_tab(results2i)
# mf_df <- as.data.frame(signif(coef(summary(results2i)), 3))
# kable_styling(kable(mf_df, row.names = TRUE))
# gwqs_weights_tab(results2i)
# 
# final_weight <- results2i$final_weights
# final_weight[, -1] <- signif(final_weight[, -1], 3)
# kable_styling(kable(final_weight, row.names = FALSE))
# 
# # bar plot
# w_ord <- order(results2i$final_weights$`Estimate pos`)
# mean_weight_pos <- results2i$final_weights$`Estimate pos`[w_ord]
# mean_weight_neg <- results2i$final_weights$`Estimate neg`[w_ord]
# mix_name <- factor(results2i$final_weights$mix_name[w_ord], 
#                    levels = results2i$final_weights$mix_name[w_ord])
# data_plot <- data.frame(mean_weight = c(mean_weight_pos, mean_weight_neg), 
#                         mix_name = rep(mix_name, 2),
#                         index = factor(rep(c("pwqs", "nwqs"), each = length(w_ord)), 
#                                        levels = c("pwqs", "nwqs")))
# ggplot(data_plot, aes(x = mix_name, y = mean_weight)) + 
#   geom_bar(stat = "identity", color = "black") + theme_bw() +
#   theme(axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text.x = element_text(color='black'),
#         legend.position = "none") + coord_flip() + 
#   geom_hline(yintercept = 1/length(PCBs), linetype="dashed", color = "red") +
#   facet_wrap(~ index)
# #
# # scatter plot y vs wqs
# ggplot(melt(results2i$y_wqs_df, measure.vars = c("pwqs", "nwqs")), aes(value, y_adj)) + 
#   geom_point() + facet_wrap(~ variable) + xlab("wqs") + 
#   stat_smooth(method = "loess", se = FALSE, linewidth = 1.5) + theme_bw() 
# #
# # scatter plot residuals vs fitted values
# fit_df <- data.frame(fitted = fitted(results2i), 
#                      resid = residuals(results2i, type = "response"))
# res_vs_fitted <- ggplot(fit_df, aes(x = fitted, y = resid)) + geom_point() + 
#   theme_bw() + xlab("Fitted values") + ylab("Residuals")
# 
# # we run the model setting the penalization term equal to 90
# results2i_l90 <- gwqs(yLBX ~ pwqs + nwqs, mix_name = PCBs, data = wqs_data, 
#                       q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
#                       family = "gaussian", seed = 2016, lambda = 90)
# 
# # we run the model setting the penalization term equal to 900
# results2i_l900 <- gwqs(yLBX ~ pwqs + nwqs, mix_name = PCBs, data = wqs_data, 
#                        q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
#                        family = "gaussian", seed = 2016, lambda = 900)
# 
# # we run the model setting the penalization term equal to 900
# results2i_l9000 <- gwqs(yLBX ~ pwqs + nwqs, mix_name = PCBs, data = wqs_data, 
#                         q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
#                         family = "gaussian", seed = 2016, lambda = 9000)
# 
# lambda_AIC_2i <- data.frame(lambda = c(0, 90, 900, 9000),
#                             AIC = c(results2i$fit$aic, results2i_l90$fit$aic, 
#                                     results2i_l900$fit$aic, results2i_l9000$fit$aic))
# kable(lambda_AIC_2i) %>% kable_styling()
# summary(results2i_l90)
# gwqs_barplot(results2i_l90)
# 
# 
# results1i <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, 
#                   q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
#                   family = "gaussian", seed = 2016)
# 
# # we run the model setting the penalization term equal to 90
# results1i_l90 <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, 
#                       q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
#                       family = "gaussian", seed = 2016, lambda = 90)
# 
# # we run the model setting the penalization term equal to 900
# results1i_l900 <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, 
#                        q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
#                        family = "gaussian", seed = 2016, lambda = 900)
# 
# # we run the model setting the penalization term equal to 900
# results1i_l9000 <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, 
#                         q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
#                         family = "gaussian", seed = 2016, lambda = 9000)
# 
# lambda_AIC_1i <- data.frame(lambda = c(0, 90, 900, 9000),
#                             AIC = c(results1i$fit$aic, results1i_l90$fit$aic, 
#                                     results1i_l900$fit$aic, results1i_l9000$fit$aic))
# kable(lambda_AIC_1i) %>% kable_styling()
# summary(results1i_l90)
# gwqs_barplot(results1i_l90)
# 
# # https://jenfb.github.io/bkmr/overview.html
# install.packages('bkmr')
# library(bkmr)
# 
# 
# set.seed(111)
# dat <- SimData(n = 50, M = 4)
# y <- dat$y
# Z <- dat$Z
# X <- dat$X
# 
# z1 <- seq(min(dat$Z[, 1]), max(dat$Z[, 1]), length = 20)
# z2 <- seq(min(dat$Z[, 2]), max(dat$Z[, 2]), length = 20)
# hgrid.true <- outer(z1, z2, function(x,y) apply(cbind(x,y), 1, dat$HFun))
# 
# res <- persp(z1, z2, hgrid.true, theta = 30, phi = 20, expand = 0.5, 
#              col = "lightblue", xlab = "", ylab = "", zlab = "")
# 
# set.seed(111)
# fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 10000, verbose = FALSE, varsel = TRUE)
# 
# TracePlot(fit = fitkm, par = "beta")
# 
# TracePlot(fit = fitkm, par = "sigsq.eps")
# 
# TracePlot(fit = fitkm, par = "r", comp = 1)
# 
# ExtractPIPs(fitkm)
# 
# med_vals <- apply(Z, 2, median)
# Znew <- matrix(med_vals, nrow = 1)
# h_true <- dat$HFun(Znew)
# 
# h_est1 <- ComputePostmeanHnew(fitkm, Znew = Znew, method = "approx")
# h_est2 <- ComputePostmeanHnew(fitkm, Znew = Znew, method = "exact")
# set.seed(111)
# samps3 <- SamplePred(fitkm, Znew = Znew, Xnew = cbind(0))
# 
# h_est_compare <- data.frame(
#   method = c("truth", 1:3),
#   post_mean = c(h_true, h_est1$postmean, h_est2$postmean, mean(samps3)),
#   post_sd = c(NA, sqrt(h_est1$postvar), sqrt(h_est2$postvar), sd(samps3))
# )
# h_est_compare
# 
# pred.resp.univar <- PredictorResponseUnivar(fit = fitkm)
# 
# library(ggplot2)
# ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
#   geom_smooth(stat = "identity") + 
#   facet_wrap(~ variable) +
#   ylab("h(z)")
# 
# pred.resp.bivar <- PredictorResponseBivar(fit = fitkm, min.plot.dist = 1)
# 
# ggplot(pred.resp.bivar, aes(z1, z2, fill = est)) + 
#   geom_raster() + 
#   facet_grid(variable2 ~ variable1) +
#   scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
#   xlab("expos1") +
#   ylab("expos2") +
#   ggtitle("h(expos1, expos2)")
# 
# pred.resp.bivar.levels <- PredictorResponseBivarLevels(
#   pred.resp.df = pred.resp.bivar, 
#   
#   Z = Z, qs = c(0.1, 0.5, 0.9))
# 
# ggplot(pred.resp.bivar.levels, aes(z1, est)) + 
#   geom_smooth(aes(col = quantile), stat = "identity") + 
#   facet_grid(variable2 ~ variable1) +
#   ggtitle("h(expos1 | quantiles of expos2)") +
#   xlab("expos1")
# 
# risks.overall <- OverallRiskSummaries(fit = fitkm, y = y, Z = Z, X = X, 
#                                       qs = seq(0.25, 0.75, by = 0.05), 
#                                       q.fixed = 0.5, method = "exact")
# risks.overall
# 
# ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
#   geom_pointrange()
# 
# 
# risks.singvar <- SingVarRiskSummaries(fit = fitkm, y = y, Z = Z, X = X, 
#                                       qs.diff = c(0.25, 0.75), 
#                                       q.fixed = c(0.25, 0.50, 0.75),
#                                       method = "exact")
# risks.singvar
# 
# ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, 
#                           ymax = est + 1.96*sd, col = q.fixed)) + 
#   geom_pointrange(position = position_dodge(width = 0.75)) + 
#   coord_flip()
# 
# risks.int <- SingVarIntSummaries(fit = fitkm, y = y, Z = Z, X = X, 
#                                  qs.diff = c(0.25, 0.75), 
#                                  qs.fixed = c(0.25, 0.75),
#                                  method = "exact")
# risks.int
# 
# set.seed(111)
# d2 <- SimData(n = 100, M = 4, Zgen = "corr", sigsq.true = 2.2)
# round(cor(d2$Z), 2)
# 
# set.seed(111)
# fitkm_corr <- kmbayes(y = d2$y, Z = d2$Z, X = d2$X, iter = 10000, varsel = TRUE, verbose = FALSE)
# fitkm_hier <- kmbayes(y = d2$y, Z = d2$Z, X = d2$X, iter = 10000, varsel = TRUE, 
#                       groups = c(1,2,1,3), verbose = FALSE)
# 
# ExtractPIPs(fitkm_corr)
# ExtractPIPs(fitkm_hier)
# data.frame(fitkm$control.params)
# data.frame(fitkm$control.params)
# priorfits <- InvestigatePrior(y = y, Z = Z, X = X, 
#                               q.seq = c(2, 1/2, 1/4, 1/16))
# PlotPriorFits(y = y, Z = Z, X = X, 
#               fits = priorfits)
# 
# #Testing HIMA... ei näytä olevan monta x:ä ja y:ä:
# install.packages("HIMA")
# library(HIMA)
# library(qvalue)
# 
# dblassohima.fit <- dblassoHIMA(X = himaDat$Example1$PhenoData$Treatment,
#                                Y = himaDat$Example1$PhenoData$Outcome,
#                                M = himaDat$Example1$Mediator,
#                                Z = himaDat$Example1$PhenoData[, c("Sex", "Age")],
#                                Y.family = 'gaussian',
#                                scale = FALSE,
#                                verbose = TRUE)
# 
# 
# setwd("C:/Users/patati/Desktop/TurkuOW/RWork/tests6/tests_basic") #tests6.... 30.11.2023
# Outcome=colnames(tv_all)[c(29:70,80:90)];
# Treatment=colnames(tv_all)[71:78];
# Mediator=colnames(tv_all)[9:28]; X <- tv_all[,Treatment];
# X <- tv_all[,Treatment];
# t.val='minmax';simss=10; test=''; take=''; sick='yes'; sick_group=sick_group
# 
# hima.fit <- hima(X = himaDat$Example1$PhenoData$Treatment,
#                  Y = himaDat$Example1$PhenoData$Outcome,
#                  M = himaDat$Example1$Mediator,
#                  COV.XM = himaDat$Example1$PhenoData[, c("Sex", "Age")],
#                  Y.family = 'gaussian',
#                  scale = FALSE,
#                  verbose = TRUE)
# hima.fit
# 
# 
# 
# 
# #####################################
# # This code is a sample to generate simulation data
# # and run high-dimensional mediation analysis.
# #####################################
# library(MASS)
# library(hdi)
# library(HDMT)
# 
# n <- 400 # sample size
# p <- 1000 # the dimension of mediators
# q <- 2 # the number of adjusted covariates
# rou <- 0.25  # the correlation of exposure
# 
# # the regression coefficients alpha (exposure --> mediators)
# alpha <- matrix(0,1,p)
# 
# # the regression coefficients beta (mediators --> outcome)
# beta <- matrix(0,1,p)
# 
# # the first five markers are true mediators.
# alpha[1:5] <- c(0.20,0.25,0.15,0.30,0.35)
# beta[1:5] <- c(0.20,0.25,0.15,0.30,0.35)
# 
# alpha[6] <- 0.1
# beta[7] <- 0.1
# 
# ## Regression coefficients eta (covariates --> outcome)
# eta <- matrix(0.3,p,q)
# ## the regression coefficients gamma (exposure --> outcome)
# gamma <- matrix(0.5,1,1)
# ## the regression coefficients delta (covariates --> mediator)
# delta <- matrix(0.5,1,q)
# 
# # Generate simulation data
# simdat = simHIMA2(n,p,q,rou,alpha,beta,seed=1234)
# 
# # HIMA2 output
# hima2.fit <- HIM

#Learning from errors:
# Error in paste(...) : cannot coerce type 'closure' to vector of type 'character' #something is missing somehwere, esp. xlsx function namings 'paste' needs all variables

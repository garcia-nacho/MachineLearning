#MyHerritage analyzer
#nacho garcia 2019
#garcia.nacho@gmail.com

library(ggplot2)
library(SNPediaR)

#data loading skipping the first 12 rows because they are comments in the file
df <- read.csv("/home/nacho/MyHerr/MyHeritage_raw_dna_dataMV.csv", sep = ",", skip = 12)

#Generation of files for www.snp-nexus.org containing 100K SNPS
n.files <- round(nrow(df)/100000) + 1

for (i in 1:n.files){
  start<-((i-1)*100000)+1
  end<- min(i* 100000,nrow(df))
  dummy<-as.data.frame(rep("dbsnp", end-start+1))
  dummy$rsid<-df[start:end,1]

  write.table(dummy, file = paste("/home/nacho/MyHerr/SNPS_DB",i,".txt",sep = ""),
              col.names=FALSE,
              quote = FALSE,
              row.names = FALSE,
              sep = "\t")
}

#Populations loading 
path <- "/home/nacho/MyHerr/data/populations/"
files<-list.files(path)
setwd(path)

pop<-list()

pb <- txtProgressBar(min = 1, max = length(files), initial = 1) 
for (i in 1:length(files)) {
  setTxtProgressBar(pb,i)
  pop[[i]]<-read.csv(paste(path,files[i],sep = ""), sep = "\t")
}

#Subset populations 
files <- gsub("-.*","",files)

#1000genomes
code1000genomes<-c("afr", "eur", "sas", "amr", "eas")
pop1000genomes<-pop[which(files %in% code1000genomes)]
pop<-pop[-which(files %in% code1000genomes)]
code1000genomes<-files[which(files %in% code1000genomes)]

#EAC
codeEAC<-c("AFR", "AMR", "OTH", "EAS", "FIN", "NFE", "SAS")
popEAC<-pop[which(files %in% codeEAC)]
pop<-pop[-which(files %in% codeEAC)]
codeEAC<-files[which(files %in% codeEAC)]

##### 1000genomes analysis 
#Merge data with sample 
priors1000genomes<- c(0.17, #afr
                      0.15, #amr
                      0.27, #eas
                      0.1, #eur
                      0.31 #sas
                      )
#Check that all of them are equal
rsid<-pop1000genomes[[1]]$SNP
for (i in 2:length(code1000genomes)) {
  rsid<-Reduce(intersect, list(pop1000genomes[[i]]$SNP, rsid ))  
}
rsid<- Reduce(intersect, list(df$RSID, rsid ))

df1000g<-subset(df, df$RSID %in% rsid )
df1000g<- df1000g[,c(1,4)]    
#df1000g$Allele1<- gsub("^.","",df1000g$RESULT)
#df1000g$Allele2<- gsub(".$","",df1000g$RESULT)
#df1000g$RESULT<-NULL
#colnames(df1000g)<-c("RSID","Allele1","Allele2")


  for (i in 1:length(code1000genomes)) {
    pop1000genomes[[i]]<-subset(pop1000genomes[[i]], pop1000genomes[[i]]$SNP %in% rsid)
    #changes zeros and ones
    pop1000genomes[[i]]$Frequency<-as.numeric(as.character(pop1000genomes[[i]]$Frequency))
    pop1000genomes[[i]]$Frequency[pop1000genomes[[i]]$Frequency==0]<-0.00005
    pop1000genomes[[i]]$Frequency[pop1000genomes[[i]]$Frequency==1]<-0.99995
    pop1000genomes[[i]]<-pop1000genomes[[i]][order(pop1000genomes[[i]]$SNP),]  
    }

#Repeated elemnts

for (i in 1:length(code1000genomes)) {
  pop1000genomes[[i]][,c(2,3,7)]<-NULL
  
  colnames(pop1000genomes[[i]])<-paste(colnames(pop1000genomes[[i]]), code1000genomes[i], sep = "-")
  colnames(pop1000genomes[[i]])<- gsub("SNP-.*", "RSID", colnames(pop1000genomes[[i]]))
  
  #merging #Population size 5000
  if (i==1) df1000g <- merge( pop1000genomes[[i]],df1000g, by="RSID")
  if (i>1) df1000g <- cbind(df1000g,pop1000genomes[[i]] )

  list<- c(1:length(code1000genomes))
  list <- list[-i]
  
  Rest <- 0
  RestAlt1N<-0
  RestRef<-0
  
  for (h in 1:length(list)) {
#Modify rests
    Rest <- ((as.numeric(as.character(pop1000genomes[[list[h]]]$Frequency))^2) * priors1000genomes[[list[h]]]) + Rest
    
    RestAlt1N <- ((2*as.numeric(as.character(pop1000genomes[[list[h]]]$Frequency))-
                   2*as.numeric(as.character(pop1000genomes[[list[h]]]$Frequency))^2) * 
                   priors1000genomes[[list[h]]]) + RestAlt1N
    
    RestRef <- (1-2*as.numeric(as.character(pop1000genomes[[list[h]]]$Frequency))+
                  as.numeric(as.character(pop1000genomes[[list[h]]]$Frequency))^2) *
                  priors1000genomes[[list[h]]] + RestRef
  }
  
  df1000g$last<-Rest
  colnames(df1000g)[ncol(df1000g)] <- paste("Rest.2N.Alt",code1000genomes[i], sep = "-")
  df1000g$last<-RestAlt1N
  colnames(df1000g)[ncol(df1000g)] <- paste("Rest.1N.Alt",code1000genomes[i], sep = "-")
  df1000g$last<-RestRef
  colnames(df1000g)[ncol(df1000g)] <- paste("Rest.2N.Ref",code1000genomes[i], sep = "-")

  }

df1000g <- df1000g[complete.cases(df1000g),]

#df1000g <- df1000g[,-c(7,8,12,13,17,18,22,23)]
df1000g <- subset(df1000g,df1000g$RESULT!="--")
df1000g<-df1000g[which(nchar(as.character(df1000g$`Ref_Allele-afr`))==1),]
df1000g<-df1000g[which(nchar(as.character(df1000g$`Alt_Allele-afr`))==1),]
df1000g<-df1000g[which(nchar(as.character(df1000g$`Alt_Allele-eur`))==1),]
df1000g<-df1000g[which(nchar(as.character(df1000g$`Alt_Allele-sas`))==1),]
df1000g<-df1000g[which(nchar(as.character(df1000g$`Alt_Allele-amr`))==1),]
df1000g<-df1000g[which(nchar(as.character(df1000g$`Alt_Allele-eas`))==1),]

df1000g$Ref.2N <- paste(df1000g$`Ref_Allele-afr`, df1000g$`Ref_Allele-afr`, sep = "")

#AFR
df1000g$Alt.1NA.afr <- paste(df1000g$`Ref_Allele-afr`, df1000g$`Alt_Allele-afr`, sep = "")
df1000g$Alt.1NB.afr <- paste( df1000g$`Alt_Allele-afr`, df1000g$`Ref_Allele-afr`, sep = "")
df1000g$Alt.2N.afr<-  paste(df1000g$`Alt_Allele-afr`, df1000g$`Alt_Allele-afr`, sep = "")
df1000g$Ref.2N.afr.Frq <- (1+ df1000g$`Frequency-afr`^ 2 - 2*df1000g$`Frequency-afr`)
df1000g$Alt.2N.afr.Frq <- df1000g$`Frequency-afr`^ 2

#EUR
df1000g$Alt.1NA.eur <- paste(df1000g$`Ref_Allele-eur`, df1000g$`Alt_Allele-eur`, sep = "")
df1000g$Alt.1NB.eur <- paste( df1000g$`Alt_Allele-eur`, df1000g$`Ref_Allele-eur`, sep = "")
df1000g$Alt.2N.eur<-  paste(df1000g$`Alt_Allele-eur`, df1000g$`Alt_Allele-eur`, sep = "")
df1000g$Ref.2N.eur.Frq <- (1+ df1000g$`Frequency-eur`^ 2 - 2*df1000g$`Frequency-eur`)
df1000g$Alt.2N.eur.Frq <- df1000g$`Frequency-eur`^ 2

#SAS
df1000g$Alt.1NA.sas <- paste(df1000g$`Ref_Allele-sas`, df1000g$`Alt_Allele-sas`, sep = "")
df1000g$Alt.1NB.sas <- paste( df1000g$`Alt_Allele-sas`, df1000g$`Ref_Allele-sas`, sep = "")
df1000g$Alt.2N.sas<-  paste(df1000g$`Alt_Allele-sas`, df1000g$`Alt_Allele-sas`, sep = "")
df1000g$Ref.2N.sas.Frq <- (1+ df1000g$`Frequency-sas`^ 2 - 2*df1000g$`Frequency-sas`)
df1000g$Alt.2N.sas.Frq <- df1000g$`Frequency-sas`^ 2

#AMR
df1000g$Alt.1NA.amr <- paste(df1000g$`Ref_Allele-amr`, df1000g$`Alt_Allele-amr`, sep = "")
df1000g$Alt.1NB.amr <- paste( df1000g$`Alt_Allele-amr`, df1000g$`Ref_Allele-amr`, sep = "")
df1000g$Alt.2N.amr<-  paste(df1000g$`Alt_Allele-amr`, df1000g$`Alt_Allele-amr`, sep = "")
df1000g$Ref.2N.amr.Frq <- (1+ df1000g$`Frequency-amr`^ 2 - 2*df1000g$`Frequency-amr`)
df1000g$Alt.2N.amr.Frq <- df1000g$`Frequency-amr`^ 2

#EAS
df1000g$Alt.1NA.eas <- paste(df1000g$`Ref_Allele-eas`, df1000g$`Alt_Allele-eas`, sep = "")
df1000g$Alt.1NB.eas <- paste( df1000g$`Alt_Allele-eas`, df1000g$`Ref_Allele-afr`, sep = "")
df1000g$Alt.2N.eas<-  paste(df1000g$`Alt_Allele-eas`, df1000g$`Alt_Allele-eas`, sep = "")
df1000g$Ref.2N.eas.Frq <- (1+ df1000g$`Frequency-eas`^ 2 - 2*df1000g$`Frequency-eas`)
df1000g$Alt.2N.eas.Frq <- df1000g$`Frequency-eas`^ 2

###Probs

#AFR subsetting and calculation
df.afr.Ref.2N<- subset(df1000g[,c(1,4,5:8,37:42)], df1000g$RESULT==df1000g$Ref.2N)
df.afr.Alt.1N<- subset(df1000g[,c(1,4,5:8,37:42)], df1000g$RESULT==df1000g$Alt.1NA.afr |  df1000g$RESULT==df1000g$Alt.1NB.afr)  
df.afr.Alt.2N <- subset(df1000g[,c(1,4,5:8,37:42)], df1000g$RESULT==df1000g$Alt.2N.afr)

#EUR subsetting and calculation
df.eur.Ref.2N<- subset(df1000g[,c(5,12,13,14,15,42,43,44,46,47)], df1000g$RESULT==df1000g$Ref.2N)
df.eur.Alt.1N<- subset(df1000g[,c(5,12,13,14,15,42,43,44,46,47)], df1000g$RESULT==df1000g$Alt.1NA.eur |  df1000g$RESULT==df1000g$Alt.1NB.eur)  
df.eur.Alt.2N<- subset(df1000g[,c(5,12,13,14,15,42,43,44,46,47)], df1000g$RESULT==df1000g$Alt.2N.eur)

#SAS subsetting and calculation
df.sas.Ref.2N<- subset(df1000g[,c(5,19,20,21,22,48,49,50,51,52)], df1000g$RESULT==df1000g$Ref.2N)
df.sas.Alt.1N<- subset(df1000g[,c(5,19,20,21,22,48,49,50,51,52)], df1000g$RESULT==df1000g$Alt.1NA.sas |  df1000g$RESULT==df1000g$Alt.1NB.sas)  
df.sas.Alt.2N<- subset(df1000g[,c(5,19,20,21,22,48,49,50,51,52)], df1000g$RESULT==df1000g$Alt.2N.sas)

#EAS subsetting and calculation
df.eas.Ref.2N<- subset(df1000g[,c(5,33,34,35,36,61,62)], df1000g$RESULT==df1000g$Ref.2N)
df.eas.Alt.1N<- subset(df1000g[,c(5,33,34,35,36,61,62)], df1000g$RESULT==df1000g$Alt.1NA.eas |  df1000g$RESULT==df1000g$Alt.1NB.eas)  
df.eas.Alt.2N<- subset(df1000g[,c(5,33,34,35,36,61,62)], df1000g$RESULT==df1000g$Alt.2N.sas)

#AMR subsetting and calculation
df.amr.Ref.2N<- subset(df1000g[,c(5,26,27,28,29,56,57)], df1000g$RESULT==df1000g$Ref.2N)
df.amr.Alt.1N<- subset(df1000g[,c(5,26,27,28,29,56,57)], df1000g$RESULT==df1000g$Alt.1NA.sas |  df1000g$RESULT==df1000g$Alt.1NB.sas)  
df.amr.Alt.2N<- subset(df1000g[,c(5,26,27,28,29,56,57)], df1000g$RESULT==df1000g$Alt.2N.sas)

#Total probability for 

#Probability calculator
afr.priors <- 0.1
eur.priors <- 0.1
sas.priors<-0.3
amr.priors<-0.2
eas.priors<-0.3

N.elements<-150
#Probability calculator


stats<- as.data.frame(array(data = NA, dim = c(N.elements,6)))
colnames(stats)<-c("N",code1000genomes)
rnd<-round(runif(1, min = N.elements, max = 50000))
  
for(n in 1:N.elements){
  stats$N[n]<-n*3

  #EUR
  num.eur<-eur.priors*
  prod(df.eur.Ref.2N$Ref.2N.eur.Frq[1:n+rnd])*
  prod(2*df.eur.Alt.1N$`Frequency-eur`[1:n+rnd] - 2*df.eur.Alt.1N$`Frequency-eur`[1:n+rnd]^2)*
  prod(df.eur.Alt.2N$Alt.2N.eur.Frq[1:n+rnd])

den.eur<-num.eur+
  ((1-eur.priors)*
     prod(df.eur.Ref.2N$`Rest.2N.Ref-eur`[1:n+rnd])*
     prod(df.eur.Alt.1N$`Rest.1N.Alt-eur`[1:n+rnd])*
     prod(df.eur.Alt.2N$`Rest.2N.Alt-eur`[1:n+rnd]))

prob.eur<- num.eur/den.eur
  
  #AFR
  num.afr<-afr.priors*
  prod(df.afr.Ref.2N$Ref.2N.afr.Frq[1:n+rnd])*
  prod(2*df.afr.Alt.1N$`Frequency-afr`[1:n+rnd] - 2*df.afr.Alt.1N$`Frequency-afr`[1:n+rnd]^2)*
  prod(df.afr.Alt.2N$Alt.2N.afr.Frq[1:n+rnd])

den.afr<-num.afr+
  ((1-afr.priors)*
     prod(df.afr.Ref.2N$`Rest.2N.Ref-afr`[1:n+rnd])*
     prod(df.afr.Alt.1N$`Rest.1N.Alt-afr`[1:n+rnd])*
     prod(df.afr.Alt.2N$`Rest.2N.Alt-afr`[1:n+rnd]))

prob.afr<- num.afr/den.afr
 
 #SAS
  num.sas<-sas.priors*
  prod(df.sas.Ref.2N$Ref.2N.sas.Frq[1:n+rnd])*
  prod(2*df.sas.Alt.1N$`Frequency-sas`[1:n+rnd] - 2*df.sas.Alt.1N$`Frequency-sas`[1:n+rnd]^2)*
  prod(df.sas.Alt.2N$Alt.2N.sas.Frq[1:n+rnd])

den.sas<-num.sas+
  ((1-sas.priors)*
     prod(df.sas.Ref.2N$`Rest.2N.Ref-sas`[1:n+rnd])*
     prod(df.sas.Alt.1N$`Rest.1N.Alt-sas`[1:n+rnd])*
     prod(df.sas.Alt.2N$`Rest.2N.Alt-sas`[1:n+rnd]))
prob.sas<- num.sas/den.sas

  #EAS
  num.eas<-eas.priors*
  prod(df.eas.Ref.2N$Ref.2N.eas.Frq[1:n+rnd])*
  prod(2*df.eas.Alt.1N$`Frequency-eas`[1:n+rnd] - 2*df.eas.Alt.1N$`Frequency-eas`[1:n+rnd]^2)*
  prod(df.eas.Alt.2N$Alt.2N.eas.Frq[1:n+rnd])

den.eas<-num.eas+
  ((1-eas.priors)*
     prod(df.eas.Ref.2N$`Rest.2N.Ref-eas`[1:n+rnd])*
     prod(df.sas.Alt.1N$`Rest.1N.Alt-eas`[1:n+rnd])*
     prod(df.eas.Alt.2N$`Rest.2N.Alt-eas`[1:n+rnd]))

prob.eas<- num.eas/den.eas

  #AMR
  num.amr<-amr.priors*
  prod(df.amr.Ref.2N$Ref.2N.amr.Frq[1:n+rnd])*
  prod(2*df.amr.Alt.1N$`Frequency-amr`[1:n+rnd] - 2*df.amr.Alt.1N$`Frequency-amr`[1:n+rnd]^2)*
  prod(df.amr.Alt.2N$Alt.2N.amr.Frq[1:n+rnd])

den.amr<-num.amr+
  ((1-amr.priors)*
     prod(df.amr.Ref.2N$`Rest.2N.Ref-amr`[1:n+rnd])*
     prod(df.amr.Alt.1N$`Rest.1N.Alt-amr`[1:n+rnd])*
     prod(df.amr.Alt.2N$`Rest.2N.Alt-amr`[1:n+rnd]))

prob.amr<- num.amr/den.amr

#Probabilities

stats$afr[n]<-prob.afr
stats$eur[n]<-prob.eur
stats$sas[n]<-prob.sas
stats$eas[n]<-prob.eas
stats$amr[n]<-prob.amr

}

ggplot(stats)+
  geom_line(aes(x=N, y= afr), colour="red")+
  geom_line(aes(x=N, y= eur), colour="blue")+
  geom_line(aes(x=N, y= sas), colour="green")+
  geom_line(aes(x=N, y= amr), colour="black")+
  geom_line(aes(x=N, y= eas), colour="yellow")+
  theme_minimal()+
  xlab("Number of SNPs")+
  ylab("Probability")

for (i in 1:nrow(stats)) {
  stats$AFR[i]<-stats$afr[i]/sum(stats[i,c(2:6)])
  stats$EUR[i]<-stats$eur[i]/sum(stats[i,c(2:6)])
  stats$AMR[i]<-stats$amr[i]/sum(stats[i,c(2:6)])
  stats$EAS[i]<-stats$eas[i]/sum(stats[i,c(2:6)])
  stats$SAS[i]<-stats$sas[i]/sum(stats[i,c(2:6)])
  }

ggplot(stats)+
  geom_line(aes(x=N, y= afr), colour="red")+
  geom_line(aes(x=N, y= eur), colour="blue")+
  geom_line(aes(x=N, y= sas), colour="green")+
  geom_line(aes(x=N, y= amr), colour="black")+
  geom_line(aes(x=N, y= eas), colour="yellow")+
  theme_minimal()+
  xlab("Number of SNPs")+
  ylab("Probability")


#Random bootstrap and see

#Clasical bayes

for(i in 1:ncol()){
  
  
  
}



#MySNP

listsnp<-read.csv("/home/nacho/MyHerr/rsidlist.txt")
listsnp<-gsub("Rs","rs", listsnp$I1000001)
df.list<-df$RSID

listsnp<-listsnp[listsnp %in% df.list]

listsnp[1]

pg<-getPages (titles = listsnp[2])

extractSnpTags (pg[[1]])

#Atopic eczema

atopic<-read.csv("/home/nacho/MyHerr/atopic eczema.csv", sep = " ", header = FALSE)
atopic$RSID<-gsub("\t.*", "", atopic$V1)
atopic$Effect<-gsub("/.*", "", atopic$V1)
atopic$Effect<-gsub(".*\t","",atopic$Effect)
atopic$Other<-gsub(".*/", "", atopic$V1)
atopic$Other<-gsub("\t.*","",atopic$Other)
atopic<-atopic[,c(5:7)]
df.atopic<-df[which(df$RSID %in% atopic$RSID),c(1,4)]
df.atopic<-merge(df.atopic, atopic[which(atopic$RSID %in% df.atopic$RSID),])


#AS

#Atopic eczema

as<-read.csv("/home/nacho/MyHerr/arthritis", sep = " ", header = FALSE)
as<-as$V1
as<-gsub("\t.*","",as)
df.as<-df[which(df$RSID %in% as),c(1,4)]

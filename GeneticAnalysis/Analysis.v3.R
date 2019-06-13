#MyHerritage analyzer
#nacho garcia 2019
#garcia.nacho@gmail.com

library(ggplot2)
library(SNPediaR)
library(chromoMap)

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

##### 1000genomes analysis 
#Merge data with sample 
priors1000genomes<- c(0.17, #afr
                      0.1, #eur
                      0.31, #sas
                      0.15, #amr
                      0.27 #eas
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
df1000g$Alt.1NA.eur <- paste(df1000g$`Ref_Allele-afr`, df1000g$`Alt_Allele-eur`, sep = "")
df1000g$Alt.1NB.eur <- paste( df1000g$`Alt_Allele-eur`, df1000g$`Ref_Allele-afr`, sep = "")
df1000g$Alt.2N.eur<-  paste(df1000g$`Alt_Allele-eur`, df1000g$`Alt_Allele-eur`, sep = "")
df1000g$Ref.2N.eur.Frq <- (1+ df1000g$`Frequency-eur`^ 2 - 2*df1000g$`Frequency-eur`)
df1000g$Alt.2N.eur.Frq <- df1000g$`Frequency-eur`^ 2

#SAS
df1000g$Alt.1NA.sas <- paste(df1000g$`Ref_Allele-afr`, df1000g$`Alt_Allele-sas`, sep = "")
df1000g$Alt.1NB.sas <- paste( df1000g$`Alt_Allele-sas`, df1000g$`Ref_Allele-afr`, sep = "")
df1000g$Alt.2N.sas<-  paste(df1000g$`Alt_Allele-sas`, df1000g$`Alt_Allele-sas`, sep = "")
df1000g$Ref.2N.sas.Frq <- (1+ df1000g$`Frequency-sas`^ 2 - 2*df1000g$`Frequency-sas`)
df1000g$Alt.2N.sas.Frq <- df1000g$`Frequency-sas`^ 2

#AMR
df1000g$Alt.1NA.amr <- paste(df1000g$`Ref_Allele-afr`, df1000g$`Alt_Allele-amr`, sep = "")
df1000g$Alt.1NB.amr <- paste( df1000g$`Alt_Allele-amr`, df1000g$`Ref_Allele-afr`, sep = "")
df1000g$Alt.2N.amr<-  paste(df1000g$`Alt_Allele-amr`, df1000g$`Alt_Allele-amr`, sep = "")
df1000g$Ref.2N.amr.Frq <- (1+ df1000g$`Frequency-amr`^ 2 - 2*df1000g$`Frequency-amr`)
df1000g$Alt.2N.amr.Frq <- df1000g$`Frequency-amr`^ 2

#EAS
df1000g$Alt.1NA.eas <- paste(df1000g$`Ref_Allele-afr`, df1000g$`Alt_Allele-eas`, sep = "")
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
colnames(stats)<-c("N","AFR","EUR","SAS","EAS","AMR")

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

stats$AFR[n]<-prob.afr
stats$EUR[n]<-prob.eur
stats$SAS[n]<-prob.sas
stats$EAS[n]<-prob.eas
stats$AMR[n]<-prob.amr

}

ggplot(stats)+
  geom_line(aes(x=N, y= AFR), colour="red")+
  geom_line(aes(x=N, y= EUR), colour="blue")+
  geom_line(aes(x=N, y= SAS), colour="green")+
  geom_line(aes(x=N, y= AMR), colour="black")+
  geom_line(aes(x=N, y= EAS), colour="yellow")+
  theme_minimal()+
  xlab("Number of SNPs")+
  ylab("Probability")

for (i in 1:nrow(stats)) {
  stats$AFR[i]<-stats$AFR[i]/sum(stats[i,c(2:6)])
  stats$EUR[i]<-stats$EUR[i]/sum(stats[i,c(2:6)])
  stats$AMR[i]<-stats$AMR[i]/sum(stats[i,c(2:6)])
  stats$EAS[i]<-stats$EAS[i]/sum(stats[i,c(2:6)])
  stats$SAS[i]<-stats$SAS[i]/sum(stats[i,c(2:6)])
  }

ggplot(stats)+
  geom_line(aes(x=N, y= AFR), colour="red")+
  geom_line(aes(x=N, y= EUR), colour="blue")+
  geom_line(aes(x=N, y= SAS), colour="green")+
  geom_line(aes(x=N, y= AMR), colour="black")+
  geom_line(aes(x=N, y= EAS), colour="yellow")+
  theme_minimal()+
  xlab("Number of SNPs")+
  ylab("Probability")


#Random bootstrap---------------
sample.N <- 20



N.2N<-nrow(df.eur.Alt.2N)
N.1N<-nrow(df.eur.Alt.1N)
N.Ref<-nrow(df.eas.Ref.2N)

boot <- as.data.frame(array(data = NA, dim = c(sample.N, 5))) 

pb <- txtProgressBar(min = 1, max = sample.N, initial = 1) 

for (i in 1:sample.N) {
  setTxtProgressBar(pb,i)
  samp.1N <- sample(nrow(df1000g), N.1N, replace=FALSE)
  samp.2N <- sample(nrow(df1000g[-samp.1N,]), N.2N, replace = FALSE)
  samp.Ref <- sample(nrow(df1000g[-c(samp.1N,samp.2N)]), N.Ref, replace = FALSE)
  
df.boot.1N <- df1000g[samp.1N,]
df.boot.2N <- df1000g[samp.2N,]
df.boot.Ref <- df1000g[samp.Ref,]
  
#Probability calculator
afr.priors <- 0.1
eur.priors <- 0.1
sas.priors<-0.3
amr.priors<-0.2
eas.priors<-0.3

n<-60
#Probability calculator


stats<- as.data.frame(array(data = NA, dim = c(N.elements,6)))
colnames(stats)<-c("N","AFR","EUR","SAS","EAS","AMR")

rnd<-round(runif(1, min = N.elements, max = 50000))

  #EUR
  num.eur<-eur.priors*
    prod(df.boot.Ref$Ref.2N.eur.Frq[1:n+rnd])*
    prod(2*df.boot.1N$`Frequency-eur`[1:n+rnd] - 2*df.boot.1N$`Frequency-eur`[1:n+rnd]^2)*
    prod(df.boot.2N$Alt.2N.eur.Frq[1:n+rnd])
  
  den.eur<-num.eur+
    ((1-eur.priors)*
       prod(df.boot.Ref$`Rest.2N.Ref-eur`[1:n+rnd])*
       prod(df.boot.1N$`Rest.1N.Alt-eur`[1:n+rnd])*
       prod(df.boot.2N$`Rest.2N.Alt-eur`[1:n+rnd]))
  
  prob.eur<- num.eur/den.eur
  
  #AFR
  num.afr<-afr.priors*
    prod(df.boot.Ref$Ref.2N.afr.Frq[1:n+rnd])*
    prod(2*df.boot.1N$`Frequency-afr`[1:n+rnd] - 2*df.boot.1N$`Frequency-afr`[1:n+rnd]^2)*
    prod(df.boot.2N$Alt.2N.afr.Frq[1:n+rnd])
  
  den.afr<-num.afr+
    ((1-afr.priors)*
       prod(df.boot.Ref$`Rest.2N.Ref-afr`[1:n+rnd])*
       prod(df.boot.1N$`Rest.1N.Alt-afr`[1:n+rnd])*
       prod(df.boot.2N$`Rest.2N.Alt-afr`[1:n+rnd]))
  
  prob.afr<- num.afr/den.afr
  
  #SAS
  num.sas<-sas.priors*
    prod(df.boot.Ref$Ref.2N.sas.Frq[1:n+rnd])*
    prod(2*df.boot.1N$`Frequency-sas`[1:n+rnd] - 2*df.boot.1N$`Frequency-sas`[1:n+rnd]^2)*
    prod(df.boot.2N$Alt.2N.sas.Frq[1:n+rnd])
  
  den.sas<-num.sas+
    ((1-sas.priors)*
       prod(df.boot.Ref$`Rest.2N.Ref-sas`[1:n+rnd])*
       prod(df.boot.1N$`Rest.1N.Alt-sas`[1:n+rnd])*
       prod(df.boot.2N$`Rest.2N.Alt-sas`[1:n+rnd]))
  prob.sas<- num.sas/den.sas
  
  #EAS
  num.eas<-eas.priors*
    prod(df.boot.Ref$Ref.2N.eas.Frq[1:n+rnd])*
    prod(2*df.boot.1N$`Frequency-eas`[1:n+rnd] - 2*df.boot.1N$`Frequency-eas`[1:n+rnd]^2)*
    prod(df.boot.2N$Alt.2N.eas.Frq[1:n+rnd])
  
  den.eas<-num.eas+
    ((1-eas.priors)*
       prod(df.boot.Ref$`Rest.2N.Ref-eas`[1:n+rnd])*
       prod(df.boot.1N$`Rest.1N.Alt-eas`[1:n+rnd])*
       prod(df.boot.2N$`Rest.2N.Alt-eas`[1:n+rnd]))
  
  prob.eas<- num.eas/den.eas
  
  #AMR
  num.amr<-amr.priors*
    prod(df.boot.Ref$Ref.2N.amr.Frq[1:n+rnd])*
    prod(2*df.boot.1N$`Frequency-amr`[1:n+rnd] - 2*df.boot.1N$`Frequency-amr`[1:n+rnd]^2)*
    prod(df.boot.2N$Alt.2N.amr.Frq[1:n+rnd])
  
  den.amr<-num.amr+
    ((1-amr.priors)*
       prod(df.boot.Ref$`Rest.2N.Ref-amr`[1:n+rnd])*
       prod(df.boot.1N$`Rest.1N.Alt-amr`[1:n+rnd])*
       prod(df.boot.2N$`Rest.2N.Alt-amr`[1:n+rnd]))
  
  prob.amr<- num.amr/den.amr
  
  #Probabilities

  boot[i, 1]<-prob.afr
  boot[i, 2]<-prob.eur
  boot[i, 3]<-prob.sas
  boot[i, 4]<-prob.eas
  boot[i, 5]<-prob.amr

  
}

colnames(boot)<-c("AFR","EUR","SAS","EAS","AMR")
boxplot(log(boot[,1:5]), ylab="log(P)")

for (i in 1:nrow(boot)){
  boot[i,]<-boot[i,]/sum(boot[i,])
  
}

colnames(boot)<-c("AFR","EUR","SAS","EAS","AMR")
boxplot(log(boot[,1:5]), ylab="log(P)")



#Clasical bayes----------------------#######
prior.AFR <- 0.1
prior.AMR <- 0.2
prior.EUR <- 0.1
prior.SAS <- 0.3
prior.EAS <- 0.3

df1000g$PSNP <- prior.AFR*df1000g$`Frequency-afr`+
                prior.AMR*df1000g$`Frequency-amr`+
                prior.EUR*df1000g$`Frequency-eur`+
                prior.SAS*df1000g$`Frequency-sas`+
                prior.EAS*df1000g$`Frequency-eas`

df1000g$PSNP.Alt.1N <- 2*df1000g$PSNP - 2*df1000g$PSNP^2
df1000g$PSNP.Alt.2N <- df1000g$PSNP^2
df1000g$PSNP.Ref.2N <- 1- 2*df1000g$PSNP + df1000g$PSNP^2
df1000g$Bayes.AFR<-NA
df1000g$Bayes.EUR<-NA
df1000g$Bayes.AMR<-NA
df1000g$Bayes.SAS<-NA
df1000g$Bayes.EAS<-NA


size<-nrow(df1000g)
size<-1000

pb <- txtProgressBar(min = 1, max = size, initial = 1) 


for(i in 1:size){
  setTxtProgressBar(pb,i)  
  
  
  prior.AFR <- 0.1
  prior.AMR <- 0.2
  prior.EUR <- 0.1
  prior.SAS <- 0.3
  prior.EAS <- 0.3
  
  
  rnd<-round(runif(1,min = 1,max = 3))
    
  #Homozigous Reference
    if(df1000g$RESULT[i]==df1000g$Ref.2N[i]){
    #if(rnd==1){
        
    prior.AFR <- prior.AFR*(1- (2*df1000g$`Frequency-afr`[i])+ (df1000g$`Frequency-afr`[i]^ 2))/
                 ((prior.AFR*(1- (2*df1000g$`Frequency-afr`[i])+ (df1000g$`Frequency-afr`[i]^ 2)))+
                    ((1-prior.AFR)*df1000g$`Rest.2N.Ref-afr`[i]))
    
    prior.EUR <- prior.EUR*(1- (2*df1000g$`Frequency-eur`[i])+ (df1000g$`Frequency-eur`[i]^ 2))/
      ((prior.EUR*(1- (2*df1000g$`Frequency-eur`[i])+ (df1000g$`Frequency-eur`[i]^ 2)))+
         ((1-prior.EUR)*df1000g$`Rest.2N.Ref-eur`[i]))
    
    prior.AMR <- prior.AMR*(1- (2*df1000g$`Frequency-amr`[i])+ (df1000g$`Frequency-amr`[i]^ 2))/
      ((prior.AMR*(1- (2*df1000g$`Frequency-amr`[i])+ (df1000g$`Frequency-amr`[i]^ 2)))+
         ((1-prior.AMR)*df1000g$`Rest.2N.Ref-amr`[i]))
    
    prior.SAS <- prior.SAS*(1- (2*df1000g$`Frequency-sas`[i])+ (df1000g$`Frequency-sas`[i]^ 2))/
      ((prior.SAS*(1- (2*df1000g$`Frequency-sas`[i])+ (df1000g$`Frequency-sas`[i]^ 2)))+
         ((1-prior.SAS)*df1000g$`Rest.2N.Ref-sas`[i]))
    
    prior.EAS <- prior.EAS*(1- (2*df1000g$`Frequency-eas`[i])+ (df1000g$`Frequency-eas`[i]^ 2))/
      ((prior.EAS*(1- (2*df1000g$`Frequency-eas`[i])+ (df1000g$`Frequency-eas`[i]^ 2)))+
         ((1-prior.EAS)*df1000g$`Rest.2N.Ref-eas`[i]))
    
  }
  
  #Homozigous Alt
  if(df1000g$RESULT[i]==df1000g$Alt.2N.afr[i]){
  #if(rnd==2){
      
    prior.AFR <- prior.AFR*(df1000g$`Frequency-afr`[i]^ 2)/
      (prior.AFR*(df1000g$`Frequency-afr`[i]^ 2)+((1-prior.AFR)*df1000g$`Rest.2N.Alt-afr`[i]))
      
    prior.EUR <- prior.EUR*(df1000g$`Frequency-eur`[i]^ 2)/
      (prior.EUR*(df1000g$`Frequency-eur`[i]^ 2)+((1-prior.EUR)*df1000g$`Rest.2N.Alt-eur`[i]))
    
    prior.AMR <- prior.AMR*(df1000g$`Frequency-amr`[i]^ 2)/
      (prior.AMR*(df1000g$`Frequency-amr`[i]^ 2)+((1-prior.AMR)*df1000g$`Rest.2N.Alt-amr`[i]))
    
    prior.SAS <- prior.SAS*(df1000g$`Frequency-sas`[i]^ 2)/
      (prior.SAS*(df1000g$`Frequency-sas`[i]^ 2)+((1-prior.SAS)*df1000g$`Rest.2N.Alt-sas`[i]))
    
    prior.EAS <- prior.EAS*(df1000g$`Frequency-eas`[i]^ 2)/
      (prior.EAS*(df1000g$`Frequency-eas`[i]^ 2)+((1-prior.EAS)*df1000g$`Rest.2N.Alt-eas`[i]))
    
  }
  
  #Heterozigous Alt
  if(df1000g$RESULT[i]!=df1000g$Ref.2N[i] & df1000g$RESULT[i]!=df1000g$Ref.2N[i] ){
  #if(rnd==3){
      
    prior.AFR <- prior.AFR*((2*df1000g$`Frequency-afr`[i])-(2*df1000g$`Frequency-afr`[i]^ 2))/
    ((prior.AFR*((2*df1000g$`Frequency-afr`[i])-(2*df1000g$`Frequency-afr`[i]^ 2)))+(1-prior.AFR)*df1000g$`Rest.1N.Alt-afr`[i])
    
    prior.EUR <- prior.EUR*((2*df1000g$`Frequency-eur`[i])-(2*df1000g$`Frequency-eur`[i]^ 2))/
      ((prior.EUR*((2*df1000g$`Frequency-eur`[i])-(2*df1000g$`Frequency-eur`[i]^ 2)))+(1-prior.EUR)*df1000g$`Rest.1N.Alt-eur`[i])
    
    prior.AMR <- prior.AMR*((2*df1000g$`Frequency-amr`[i])-(2*df1000g$`Frequency-amr`[i]^ 2))/
      ((prior.AMR*((2*df1000g$`Frequency-amr`[i])-(2*df1000g$`Frequency-amr`[i]^ 2)))+(1-prior.AMR)*df1000g$`Rest.1N.Alt-amr`[i])
    
    prior.SAS <- prior.SAS*((2*df1000g$`Frequency-sas`[i])-(2*df1000g$`Frequency-sas`[i]^ 2))/
      ((prior.SAS*((2*df1000g$`Frequency-sas`[i])-(2*df1000g$`Frequency-sas`[i]^ 2)))+(1-prior.SAS)*df1000g$`Rest.1N.Alt-sas`[i])
    
    
  }
  
  
  df1000g$Bayes.AFR[i]<-prior.AFR
  df1000g$Bayes.EUR[i]<-prior.EUR
  df1000g$Bayes.AMR[i]<-prior.AMR
  df1000g$Bayes.SAS[i]<-prior.SAS
  df1000g$Bayes.EAS[i]<-prior.EAS
  
}
check<-df1000g[,c(1:5,12,19,26,33,63:71)]
boxplot(check[,14:18])

colnames(check)<- paste(colnames(check),"Samp")
check$ID<-"Sample"
randomized$ID<-"Random"

a<-check[,14:19]
b<-randomized[,14:19]

colnames(a)<-c("AFR","EUR","AMR","SAS","EAS","ID")
colnames(b)<-c("AFR","EUR","AMR","SAS","EAS","ID")

toplot <- rbind(a,b)

ggplot(toplot)+
  geom_density(aes(x=SAS, fill=ID), alpha=0.5)+
  theme_minimal()

#AFR-SNP very frequent
prior.AFR <- 0.1
prior.AMR <- 0.2
prior.EUR <- 0.1
prior.SAS <- 0.3
prior.EAS <- 0.3

#Ratios

df1000g$Overall.Freq<- df1000g$`Frequency-afr`*prior.AFR+
                       df1000g$`Frequency-amr`*prior.AMR+
                       df1000g$`Frequency-eur`*prior.EUR+
                       df1000g$`Frequency-sas`*prior.SAS+
                       df1000g$`Frequency-eas`*prior.EAS
  
  
df1000g$AFR.ratio<-df1000g$`Frequency-afr`/df1000g$Overall.Freq
df1000g$EUR.ratio<-df1000g$`Frequency-eur`/df1000g$Overall.Freq
df1000g$AMR.ratio<-df1000g$`Frequency-amr`/df1000g$Overall.Freq
df1000g$SAS.ratio<-df1000g$`Frequency-sas`/df1000g$Overall.Freq
df1000g$EAS.ratio<-df1000g$`Frequency-eas`/df1000g$Overall.Freq

#Check Ratios
SNPs.N <- 1000

df1000g<-df1000g[order(-df1000g$AFR.ratio),]
AFR.SNP<-as.character(df1000g$RSID[1:SNPs.N])

df1000g<-df1000g[order(-df1000g$EUR.ratio),]
EUR.SNP<-as.character(df1000g$RSID[1:SNPs.N])

df1000g<-df1000g[order(-df1000g$SAS.ratio),]
SAS.SNP<-as.character(df1000g$RSID[1:SNPs.N])

df1000g<-df1000g[order(-df1000g$EAS.ratio),]
EAS.SNP<-as.character(df1000g$RSID[1:SNPs.N])

df1000g<-df1000g[order(-df1000g$AMR.ratio),]
AMR.SNP<-as.character(df1000g$RSID[1:SNPs.N])

AFR.SNP.list <- as.character(df1000g$RSID [df1000g$RSID%in% AFR.SNP & df1000g$RESULT != df1000g$Ref.2N])
AFR.SNP.list1N <- as.character(df1000g$RSID [df1000g$RSID%in% AFR.SNP.list & df1000g$RESULT != df1000g$Alt.2N.afr]) 

EUR.SNP.list <- as.character(df1000g$RSID [df1000g$RSID%in% EUR.SNP & df1000g$RESULT != df1000g$Ref.2N])
EUR.SNP.list1N <- as.character(df1000g$RSID [df1000g$RSID%in% EUR.SNP.list & df1000g$RESULT != df1000g$Alt.2N.eur]) 

SAS.SNP.list <- as.character(df1000g$RSID [df1000g$RSID%in% SAS.SNP & df1000g$RESULT != df1000g$Ref.2N])
SAS.SNP.list1N <- as.character(df1000g$RSID [df1000g$RSID%in% SAS.SNP.list & df1000g$RESULT != df1000g$Alt.2N.sas]) 

EAS.SNP.list <- as.character(df1000g$RSID [df1000g$RSID%in% EAS.SNP & df1000g$RESULT != df1000g$Ref.2N])
EAS.SNP.list1N <- as.character(df1000g$RSID [df1000g$RSID%in% EAS.SNP.list & df1000g$RESULT != df1000g$Alt.2N.eas]) 

AMR.SNP.list <- as.character(df1000g$RSID [df1000g$RSID%in% AMR.SNP & df1000g$RESULT != df1000g$Ref.2N])
AMR.SNP.list1N <- as.character(df1000g$RSID [df1000g$RSID%in% AMR.SNP.list & df1000g$RESULT != df1000g$Alt.2N.amr]) 

AFR.corrector<-mean(df1000g$`Frequency-afr`[df1000g$RSID %in% AFR.SNP.list])-mean(df1000g$Overall.Freq[df1000g$RSID %in% AFR.SNP.list])
EUR.corrector<-mean(df1000g$`Frequency-eur`[df1000g$RSID %in% EUR.SNP.list])-mean(df1000g$Overall.Freq[df1000g$RSID %in% EUR.SNP.list])
SAS.corrector<-mean(df1000g$`Frequency-sas`[df1000g$RSID %in% SAS.SNP.list])-mean(df1000g$Overall.Freq[df1000g$RSID %in% SAS.SNP.list])
EAS.corrector<-mean(df1000g$`Frequency-eas`[df1000g$RSID %in% EAS.SNP.list])-mean(df1000g$Overall.Freq[df1000g$RSID %in% EAS.SNP.list])
AMR.corrector<-mean(df1000g$`Frequency-amr`[df1000g$RSID %in% AMR.SNP.list])-mean(df1000g$Overall.Freq[df1000g$RSID %in% AMR.SNP.list])

AFR.perc <- AFR.corrector*(length(AFR.SNP.list1N)+(length(AFR.SNP.list)-length(AFR.SNP.list1N))*2)/(length(AFR.SNP)*2)
EUR.perc <- EUR.corrector*(length(EUR.SNP.list1N)+(length(EUR.SNP.list)-length(EUR.SNP.list1N))*2)/(length(EUR.SNP)*2)
AMR.perc <- SAS.corrector*(length(AMR.SNP.list1N)+(length(AMR.SNP.list)-length(AMR.SNP.list1N))*2)/(length(AMR.SNP)*2)
SAS.perc <- EAS.corrector*(length(SAS.SNP.list1N)+(length(SAS.SNP.list)-length(SAS.SNP.list1N))*2)/(length(SAS.SNP)*2)
EAS.perc <- AMR.corrector*(length(EAS.SNP.list1N)+(length(EAS.SNP.list)-length(EAS.SNP.list1N))*2)/(length(EAS.SNP)*2)

WW<-sum(AFR.perc,EUR.perc,AMR.perc,SAS.perc,EAS.perc)

AFR.perc<-AFR.perc/WW
EUR.perc<-EUR.perc/WW
SAS.perc<-SAS.perc/WW
EAS.perc<-EAS.perc/WW
AMR.perc<-AMR.perc/WW

barplot(c(AFR.perc,EUR.perc,SAS.perc,EAS.perc,AMR.perc), names.arg =c("AFR","EUR","SAS","EAS","AMR"), xlab = "Population",
        ylab = "Percentage Normalized")



##CLustering method 

library("cluster")
library("factoextra")
library("magrittr")

SampleSize<-20
RelevantSNP <- c(AFR.SNP,
                 EUR.SNP,
                 AMR.SNP,
                 SAS.SNP,
                 EAS.SNP)

#Divide genome in 200 regions to plot it afterwards

RelevantSNP<-sample(as.character(df1000g$RSID),size = 500000, replace = FALSE)

population<-as.data.frame(matrix(data = NA, nrow = (SampleSize*5)+1, ncol = length(RelevantSNP) ))

dfpop<-subset(df1000g, df1000g$RSID %in% RelevantSNP)

dfpop$ploidy<-1
dfpop$ploidy[which(dfpop$RESULT==dfpop$Ref.2N)]<-0
dfpop$ploidy[which(dfpop$RESULT==dfpop$Alt.2N.afr)]<-2

population[1,]<-dfpop$ploidy
rownames(population)[1]<-"Sample"

pb <- txtProgressBar(min = 1, max = SampleSize*nrow(dfpop), initial = 1) 
count<-0

 #AFR generation
for (i in 2:(SampleSize+1)) {
  dfpop$dummy<-0
  for (h in 1:nrow(dfpop)) {
    count<-count+1
    setTxtProgressBar(pb,count)  
    rnd<-runif(1,min = 0, max = 1)
    if(rnd<dfpop$`Frequency-afr`[h]^2) dfpop$dummy[h]<-2
    if(rnd>(dfpop$`Frequency-afr`[h]^2) & rnd<(2*dfpop$`Frequency-afr`[h]-2*dfpop$`Frequency-afr`[h]^2)) dfpop$dummy[h]<-1     
  }
  population[i,]<-dfpop$dummy
  rownames(population)[i]<-paste("AFR",i,sep = "")

}

pb <- txtProgressBar(min = 1, max = SampleSize*nrow(dfpop), initial = 1) 
count<-0
#EUR generation
for (i in (SampleSize+2):(2*SampleSize+1)) {
  dfpop$dummy<-0
  for (h in 1:nrow(dfpop)) {
    count<-count+1
    setTxtProgressBar(pb,count)  
    rnd<-runif(1,min = 0, max = 1)
    if(rnd<dfpop$`Frequency-eur`[h]^2) dfpop$dummy[h]<-2
    if(rnd>(dfpop$`Frequency-eur`[h]^2) & rnd<(2*dfpop$`Frequency-eur`[h]-2*dfpop$`Frequency-eur`[h]^2)) dfpop$dummy[h]<-1     
  }
  population[i,]<-dfpop$dummy
  rownames(population)[i]<-paste("EUR",i,sep = "")
  
}

pb <- txtProgressBar(min = 1, max = SampleSize*nrow(dfpop), initial = 1) 
count<-0
#AMR generation
for (i in (2*SampleSize+2):(3*SampleSize+1)) {
  dfpop$dummy<-0
  for (h in 1:nrow(dfpop)) {
    count<-count+1
    setTxtProgressBar(pb,count)  
    rnd<-runif(1,min = 0, max = 1)
    if(rnd<dfpop$`Frequency-amr`[h]^2) dfpop$dummy[h]<-2
    if(rnd>(dfpop$`Frequency-amr`[h]^2) & rnd<(2*dfpop$`Frequency-amr`[h]-2*dfpop$`Frequency-amr`[h]^2)) dfpop$dummy[h]<-1     
  }
  population[i,]<-dfpop$dummy
  rownames(population)[i]<-paste("AMR",i,sep = "")
  
}

pb <- txtProgressBar(min = 1, max = SampleSize*nrow(dfpop), initial = 1) 
count<-0
#SAS generation
for (i in (3*SampleSize+2):(4*SampleSize+1)) {
  dfpop$dummy<-0
  for (h in 1:nrow(dfpop)) {
    count<-count+1
    setTxtProgressBar(pb,count)  
    rnd<-runif(1,min = 0, max = 1)
    if(rnd<dfpop$`Frequency-sas`[h]^2) dfpop$dummy[h]<-2
    if(rnd>(dfpop$`Frequency-sas`[h]^2) & rnd<(2*dfpop$`Frequency-sas`[h]-2*dfpop$`Frequency-sas`[h]^2)) dfpop$dummy[h]<-1     
  }
  population[i,]<-dfpop$dummy
  rownames(population)[i]<-paste("SAS",i,sep = "")
  
}

pb <- txtProgressBar(min = 1, max = SampleSize*nrow(dfpop), initial = 1) 
count<-0
#EAS generation
for (i in (4*SampleSize+2):(5*SampleSize+1)) {
  dfpop$dummy<-0
  for (h in 1:nrow(dfpop)) {
    count<-count+1
    setTxtProgressBar(pb,count)  
    rnd<-runif(1,min = 0, max = 1)
    if(rnd<dfpop$`Frequency-eas`[h]^2) dfpop$dummy[h]<-2
    if(rnd>(dfpop$`Frequency-eas`[h]^2) & rnd<(2*dfpop$`Frequency-eas`[h]-2*dfpop$`Frequency-eas`[h]^2)) dfpop$dummy[h]<-1     
  }
  population[i,]<-dfpop$dummy
  rownames(population)[i]<-paste("EAS",i,sep = "")
  
}

k.clust<- kmeans(population, 5)
k.result<-k.clust$cluster
k.result<-as.data.frame(k.result)




df.clust<-  dist(population,method = "binary")
df.clust<-hclust(df.clust,method = "ward.D2")  
fviz_dend(df.clust, k = 10,
          cex = 0.5, # label size
          
          color_labels_by_k = TRUE, 
          rect = TRUE 
)

#Cluster in chromomap

#Load positions

positions<-read.csv("/home/nacho/MyHerr/data/populations/afr-set.txt",sep = "\t")
positions<-positions[,1:3]
colnames(positions)<-c("RSID","Chr","nt")
positions<-positions[!duplicated(positions$RSID),]
df1000g<-subset(df1000g,df1000g$RSID %in% unique(positions$RSID))
df1000g<-df1000g[!duplicated(df1000g$RSID),]
positions<-subset(positions, positions$RSID %in% unique(df1000g$RSID))

df1000g<-merge(df1000g,positions, by="RSID")

df1000g$nt<-as.numeric(as.character(df1000g$nt))

ggplot(subset(df1000g,df1000g$Chr!="X" & df1000g$Chr!="Y"))+
  geom_jitter(aes(x=as.character(Chr), y=nt))+
  theme_minimal()


genomic.size<-150000
chr<-read.csv("/home/nacho/MyHerr/genome.txt.csv", sep = "\t",header = FALSE)
chr$chunks<-round(chr$V3/genomic.size)+1
chr$V1<-gsub("X","23",chr$V1)
chr$V1<-gsub("Y","24",chr$V1)


for (i in 1:nrow(chr)) {
  chr.tolook<-paste("chr",i,sep = "")
  chunks<-chr$chunks[chr$V1==chr.tolook]  
  
  pb <- txtProgressBar(min = 1, max = chunks, initial = 1) 
  for (h in 1:chunks) {
  setTxtProgressBar(pb,h)   
  RelevantSNP<-as.character(df1000g$RSID[as.numeric(df1000g$nt)>(h-1)*genomic.size & as.numeric(df1000g$nt) < h*genomic.size ])
  
  population<-as.data.frame(matrix(data = NA, nrow = (SampleSize*5)+1, ncol = length(RelevantSNP) ))
  dfpop<-subset(df1000g, df1000g$RSID %in% RelevantSNP)
  
  dfpop$ploidy<-1
  dfpop$ploidy1[which(dfpop$RESULT==dfpop$Ref.2N)]<-0
  dfpop$ploidy2[which(dfpop$RESULT==dfpop$Ref.2N)]<-0
  dfpop$ploidy1[which(dfpop$RESULT==dfpop$Alt.2N.afr)]<-1
  
  #Test etniticy in SNP.chunk     
  }

  
  
  }

#ChromoMAP
library(chromoMap)
chromoMap(c("/home/nacho/MyHerr/genome.txt.csv","/home/nacho/MyHerr/genome.txt.csv"),
          c("/home/nacho/MyHerr/anotation.csv","/home/nacho/MyHerr/anotation.csv"), ploidy = 2)



#write.csv(check,"/home/nacho/MyHerr/test.csv")

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

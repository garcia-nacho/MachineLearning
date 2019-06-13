  #Chromosome plotter
  #Nacho garcia 2019
  #garcia.nacho@gmail.com
  #GPL v3
  
  library("cluster")
  library("factoextra")
  library("magrittr")
  library(ggplot2)
  library(tcltk)
  
  
  #data loading skipping the first 12 rows because they are comments in the file
  df <- read.csv("/home/nacho/MyHerr/MyHeritage_raw_dna_dataMV.csv", sep = ",", skip = 12)
  
  
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
    pop1000genomes[[i]]$Frequency[pop1000genomes[[i]]$Frequency==0]<-0.0005
    pop1000genomes[[i]]$Frequency[pop1000genomes[[i]]$Frequency==1]<-0.9995
    pop1000genomes[[i]]<-pop1000genomes[[i]][order(pop1000genomes[[i]]$SNP),]  
  }
  
  
  
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
  
  df1000g<-df1000g[,c(1:5,12,19,26,33)]
  
  df1000g$Ref2N<-paste(df1000g$`Ref_Allele-afr`,df1000g$`Ref_Allele-afr`,sep = "")
  df1000g$Alt2N<-paste(df1000g$`Alt_Allele-afr`,df1000g$`Alt_Allele-afr`,sep = "")
  
  #Load positions
  positions<-read.csv("/home/nacho/MyHerr/data/populations/afr-set.txt",sep = "\t")
  positions<-positions[,1:3]
  colnames(positions)<-c("RSID","Chr","nt")
  positions<-positions[!duplicated(positions$RSID),]
  df1000g<-subset(df1000g,df1000g$RSID %in% unique(positions$RSID))
  df1000g<-df1000g[!duplicated(df1000g$RSID),]
  positions<-subset(positions, positions$RSID %in% unique(df1000g$RSID))
  
  df1000g<-merge(df1000g,positions, by="RSID")
  
  genomic.size<-15000000
  chr<-read.csv("/home/nacho/MyHerr/genome.txt.csv", sep = "\t",header = FALSE)
  chr$chunks<-round(chr$V3/genomic.size)+1
  chr$V1<-gsub("X","23",chr$V1)
  chr$V1<-gsub("Y","24",chr$V1)
  
  df1000g$nt<-as.numeric(as.character(df1000g$nt))
  df1000g$Eth<-"ND"
  df1000g$withinn<-NA
  df1000g$cluster.score<-NA
  df1000g$AFR.score<-NA
  df1000g$EUR.score<-NA
  df1000g$SAS.score<-NA
  df1000g$AMR.score<-NA
  df1000g$EAS.score<-NA
  
  
  df1000g$Chr<-gsub("X","23", df1000g$Chr)
  df1000g$Chr<-gsub("Y","24", df1000g$Chr)
  df1000g$Chr<-as.numeric(df1000g$Chr)
  
  #Priors
  priors1000genomes<- c(0.17, #afr
                        0.1, #eur
                        0.31, #sas
                        0.15, #amr
                        0.27 #eas
  )
  
  df1000g$Total <- df1000g$`Frequency-afr`*priors1000genomes[1]+
                   df1000g$`Frequency-eur`*priors1000genomes[2]+
                   df1000g$`Frequency-sas`*priors1000genomes[3]+
                   df1000g$`Frequency-amr`*priors1000genomes[4]+
                   df1000g$`Frequency-eas`*priors1000genomes[5]
  
  df1000g$Anti.AFR <-   (df1000g$`Frequency-eur`*priors1000genomes[2]+
    df1000g$`Frequency-sas`*priors1000genomes[3]+
    df1000g$`Frequency-amr`*priors1000genomes[4]+
    df1000g$`Frequency-eas`*priors1000genomes[5])/sum(priors1000genomes[2:5])
  
  df1000g$Anti.EUR <- (df1000g$`Frequency-afr`*priors1000genomes[1]+
    df1000g$`Frequency-sas`*priors1000genomes[3]+
    df1000g$`Frequency-amr`*priors1000genomes[4]+
    df1000g$`Frequency-eas`*priors1000genomes[5])/sum(priors1000genomes[c(1,3:5)])
  
  df1000g$Anti.AMR <- (df1000g$`Frequency-afr`*priors1000genomes[1]+
    df1000g$`Frequency-eur`*priors1000genomes[2]+
    df1000g$`Frequency-sas`*priors1000genomes[3]+
    df1000g$`Frequency-eas`*priors1000genomes[5])/sum(priors1000genomes[c(1,2,3,5)])
  
  df1000g$Anti.SAS <- (df1000g$`Frequency-afr`*priors1000genomes[1]+
    df1000g$`Frequency-eur`*priors1000genomes[2]+
    df1000g$`Frequency-amr`*priors1000genomes[4]+
    df1000g$`Frequency-eas`*priors1000genomes[5])/sum(priors1000genomes[c(1,2,4,5)])
  
  df1000g$Anti.EAS <- (df1000g$`Frequency-afr`*priors1000genomes[1]+
    df1000g$`Frequency-eur`*priors1000genomes[2]+
    df1000g$`Frequency-sas`*priors1000genomes[3]+
    df1000g$`Frequency-amr`*priors1000genomes[4])/sum(priors1000genomes[c(1:4)])
  
  
  #Strategy: check if the sample is equal to result and provide an score 
  
  # multipy it by a regularization factor ((P(local)-P(Global-L)))
  #Genetic score function
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  
  GenScore <- function(ploidy, Freq, Rest.Freq){
    Score<-rep(NA,length(Freq))
    Dip.Var <- Freq^2
    Dip.Ref <- (1-Freq)^2
    Hap<- 2*Freq*(1-Freq)
    
    Rest.Dip.Var <- Rest.Freq^2
    Rest.Dip.Ref <- (1-Rest.Freq)^2
    Rest.Hap<- 2*Rest.Freq*(1-Rest.Freq)
    
    Score[which(ploidy==2)]<- Dip.Var[which(ploidy==2)] - Rest.Dip.Var[which(ploidy==2)]
    Score[which(ploidy==0)]<- Dip.Ref[which(ploidy==0)] - Rest.Dip.Ref[which(ploidy==0)]
    Score[which(ploidy==1)]<- Hap[which(ploidy==1)] - Rest.Hap[which(ploidy==1)]
    
    return(Score)
    
  }
  
  #Cluster parallel
  SampleSize<-25 #number of individuals for the clustering
  
  library(doParallel)
  library(parallel)
  library(foreach)
  cores<-detectCores()-1
  cluster.cores<-makeCluster(cores)
  registerDoParallel(cluster.cores)
  
  
  output<-foreach(i2=1:nrow(chr), .verbose = TRUE, .packages = "cluster") %dopar%{
  
    chr.tolook<-paste("chr",i2,sep = "")
    chunks<-chr$chunks[chr$V1==chr.tolook]  
    #print(paste("Analyzing Chr",i2,sep=""))
    #pb <- txtProgressBar(min = 1, max = chunks, initial = 1) 
    for (h2 in 1:chunks) {
      #setTxtProgressBar(pb,h2)   
      RelevantSNP<-as.character(df1000g$RSID[df1000g$nt > (h2-1)*genomic.size & df1000g$nt < h2*genomic.size & df1000g$Chr==i2])
      if(length(RelevantSNP)>1){
      dfpop<-subset(df1000g, df1000g$RSID %in% RelevantSNP)
      dfpop$ploidy<-1
      dfpop$ploidy[which(dfpop$RESULT==dfpop$Ref2N)]<-0
      dfpop$ploidy[which(dfpop$RESULT==dfpop$Alt2N)]<-2
      
      population<-as.data.frame(matrix(data = NA, nrow = (SampleSize*5)+1, ncol = length(RelevantSNP) ))
      
      population[1,]<-dfpop$ploidy
      rownames(population)[1]<-"Sample"
      
      EthScore<-rep(NA, length(code1000genomes))
    
      #AFR generation
      for (i in 2:(SampleSize+1)) {
        dfpop$dummy<-0
        for (h in 1:nrow(dfpop)) {
          rnd<-runif(1,min = 0, max = 1)
          if(rnd<dfpop$`Frequency-afr`[h]^2) dfpop$dummy[h]<-2
          if(rnd>(dfpop$`Frequency-afr`[h]^2) & rnd<(2*dfpop$`Frequency-afr`[h]-2*dfpop$`Frequency-afr`[h]^2)) dfpop$dummy[h]<-1     
        }
        population[i,]<-dfpop$dummy
        rownames(population)[i]<-paste("AFR",i,sep = "")
        
      }
  
      #EUR generation
      for (i in (SampleSize+2):(2*SampleSize+1)) {
        dfpop$dummy<-0
        for (h in 1:nrow(dfpop)) {
          rnd<-runif(1,min = 0, max = 1)
          if(rnd<dfpop$`Frequency-eur`[h]^2) dfpop$dummy[h]<-2
          if(rnd>(dfpop$`Frequency-eur`[h]^2) & rnd<(2*dfpop$`Frequency-eur`[h]-2*dfpop$`Frequency-eur`[h]^2)) dfpop$dummy[h]<-1     
        }
        population[i,]<-dfpop$dummy
        rownames(population)[i]<-paste("EUR",i,sep = "")
        
      }
      
      #AMR generation
      for (i in (2*SampleSize+2):(3*SampleSize+1)) {
        dfpop$dummy<-0
        for (h in 1:nrow(dfpop)) {
          rnd<-runif(1,min = 0, max = 1)
          if(rnd<dfpop$`Frequency-amr`[h]^2) dfpop$dummy[h]<-2
          if(rnd>(dfpop$`Frequency-amr`[h]^2) & rnd<(2*dfpop$`Frequency-amr`[h]-2*dfpop$`Frequency-amr`[h]^2)) dfpop$dummy[h]<-1     
        }
        population[i,]<-dfpop$dummy
        rownames(population)[i]<-paste("AMR",i,sep = "")
        
      }
      
      #SAS generation
      for (i in (3*SampleSize+2):(4*SampleSize+1)) {
        dfpop$dummy<-0
        for (h in 1:nrow(dfpop)) {
          rnd<-runif(1,min = 0, max = 1)
          if(rnd<dfpop$`Frequency-sas`[h]^2) dfpop$dummy[h]<-2
          if(rnd>(dfpop$`Frequency-sas`[h]^2) & rnd<(2*dfpop$`Frequency-sas`[h]-2*dfpop$`Frequency-sas`[h]^2)) dfpop$dummy[h]<-1     
        }
        population[i,]<-dfpop$dummy
        rownames(population)[i]<-paste("SAS",i,sep = "")
        
      }
      
      #EAS generation
      for (i in (4*SampleSize+2):(5*SampleSize+1)) {
        dfpop$dummy<-0
        for (h in 1:nrow(dfpop)) {
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
      
      cluster<-rownames(k.result)[k.result$k.result== k.result[1,1]]
      
      cluster<-cluster[-1]
      cluster<-gsub("\\d","",cluster)
  
      df1000g$Eth[df1000g$RSID %in% RelevantSNP]<-getmode(cluster)
      df1000g$withinn[df1000g$RSID %in% RelevantSNP] <- k.clust$withinss[as.numeric(k.clust$cluster[1])]
      df1000g$cluster.score[df1000g$RSID %in% RelevantSNP] <- k.clust$size[k.clust$cluster[1]]-SampleSize+1
  
      rm(cluster)
      
      }#If there are SNP in the region
      
    }#Chunks
    
    return(subset(df1000g[c("RSID","Eth","cluster.score","withinn")], df1000g$Chr==i2))
    
  }#chr
  
  stopCluster(cluster.cores)

  out<-rbind(output[[1]],output[[2]])
  for (i in 3:length(output)){
    out<-rbind(out,output[[i]])
  }
  
  df1000g<-merge(df1000g,out, by="RSID")  



#Plotting Chr

#Check point
#write.csv(df1000g, "/home/nacho/MyHerr/df1000g.csv")

  #Scrambling
  df1000g$Eth.y<-sample(df1000g$Eth.y)
  
ggplot(subset(df1000g))+
  geom_jitter(aes(x=as.numeric(as.character(Chr)), y=nt, colour=Eth.y),size=0.05)+
  #geom_rect(aes(xmin=0.5, xmax=1.5, ymin=0, ymax=10000000), color="black", alpha=0) +
  ggtitle("Ethniticity by Chromosome")+
  xlab("Chromosome")+
  ylab(" ")+
  scale_x_discrete(limits=c(1:24), labels=c(1:22,"X","Y"))+
  scale_y_continuous(breaks=NULL)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(subset(df1000g))+
  geom_jitter(aes(x=as.numeric(as.character(Chr)), y=nt, colour=cluster.score.y),size=0.05)+
  scale_colour_gradient2(midpoint = 0, mid = "white", high = "red", low = "blue")+
  ggtitle("Clustering Score by Chromosome")+
  xlab("Chromosome")+
  ylab(" ")+
  scale_x_discrete(limits=c(1:24), labels=c(1:22,"X","Y"))+
  scale_y_continuous(breaks=NULL)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(subset(df1000g))+
  geom_jitter(aes(x=as.numeric(as.character(Chr)), y=nt, colour=withinn.y),size=0.05)+
  scale_colour_gradient2(midpoint = mean(df1000g$withinn.y, na.rm = TRUE), mid = "white", high = "red", low = "blue")+
  ggtitle("Within Score by Chromosome")+
  xlab("Chromosome")+
  ylab(" ")+
  scale_x_discrete(limits=c(1:24), labels=c(1:22,"X","Y"))+
  scale_y_continuous(breaks=NULL)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(subset(df1000g))+
  geom_jitter(aes(x=as.numeric(as.character(Chr)), y=nt, colour=SAS.score),size=0.05)+
  scale_colour_gradient2(midpoint = 0, mid = "white", high = "red", low = "blue")+
  ggtitle("South-Asian Score by Chromosome")+
  xlab("Chromosome")+
  ylab(" ")+
  scale_x_discrete(limits=c(1:24), labels=c(1:22,"X","Y"))+
  scale_y_continuous(breaks=NULL)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(subset(df1000g))+
  geom_jitter(aes(x=as.numeric(as.character(Chr)), y=nt, colour=AFR.score),size=0.05)+
  scale_colour_gradient2(midpoint = 0, mid = "white", high = "red", low = "blue")+
  ggtitle("African Score by Chromosome")+
  xlab("Chromosome")+
  ylab(" ")+
  scale_x_discrete(limits=c(1:24), labels=c(1:22,"X","Y"))+
  scale_y_continuous(breaks=NULL)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))


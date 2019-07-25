#Protein autoencoder
#Nacho Garcia 2019

library(ggplot)
library(keras)
library(rcdk)

SC<-read.csv("D:/Proteomes/cerevisiae.csv")
CA<-read.csv("D:/Proteomes/auris.csv")
CG<-read.csv("D:/Proteomes/glabrata.csv")
HS<-read.csv("D:/Proteomes/Human.csv")

SC$GN<-paste("SC_",SC$GN)
CA$GN<-paste("CA_",CA$GN)
CG$GN<-paste("CG_",CG$GN)
HS$GN<-paste("HS_",HS$GN)

prot<-rbind(SC,CA)
prot<-rbind(prot,CG)
prot<-rbind(prot,HS)

pad<-2000

size<-vector()
for(i in 1:nrow(prot)){

  size[i]<-length(strsplit(as.character(prot$sequence[i]), "")[[1]])  
}

plot(density(size))

genes<-prot$GN

df<-matrix(data = "X", nrow = nrow(prot), ncol = pad )

for(i in 1:nrow(prot)){
  uplimit<-length(strsplit(as.character(prot$sequence[i]), "")[[1]])
  if(uplimit>2000) uplimit<-2000
                 
  df[i,1:uplimit] <- strsplit(as.character(prot$sequence[1]), "")[[1]][1:uplimit]
  
}

df[is.na(df)] <- "X"


#Create similarity table for amino acids

AA<-read.csv("C:/Users/AutophagyCrusher/Documents/Proteomes/aminoacids.csv", sep = "\t", header = FALSE)
AA<-AA[,c(3,5,8,10,11,12)]
colnames(AA)<-c("Code","SMILES","Hydropaty","pI","pK1","pK2" )

#Getting descriptors of amino acids from mordred

desc<-read.csv("C:/Users/AutophagyCrusher/Documents/Proteomes/AAdescriptors.csv")

desc<-desc[vapply(desc, function(x) length(unique(x)) > 1, logical(1L))]

clean.desc<-desc[,!colSums(is.na(desc)) > 0]

pca<-prcomp(as.matrix(clean.desc[2:ncol(clean.desc)]))
plot(pca,type="l")
summary(pca)

components.n<-6

#PC6 captures 99.053 of the variance
aminoacids<-cbind(as.character(AA$Code), pca$x[,1:components.n])

#Scale between 0-1 (sigmoid function)

for (i in 2:ncol(aminoacids)) {
  aminoacids[,i]<- (as.numeric(aminoacids[,i])-as.numeric(min(aminoacids[,i])))/(as.numeric(max(aminoacids[,i]))-as.numeric(min(aminoacids[,i])))
  aminoacids[,i]<-as.numeric(as.character(aminoacids[,i]))
}

#Scale between -1 and 1 (tanh)

for (i in 2:ncol(aminoacids)) {
  aminoacids[,i]<- 2*as.numeric(aminoacids[,i]) - 1
  aminoacids[,i]<-as.numeric(as.character(aminoacids[,i]))
}


#Padding Amino acid

#aminoacids <- rbind(aminoacids, c("X", 0, 0, 0, 0, 0, 0))


#Apply PCA.aminoacids to the matrices 

#New matrix creation
aminoacids<-as.data.frame(aminoacids)

df3D <- array(data = 0, dim = c(dim(df),ncol(aminoacids)-1) ) 

pb <- txtProgressBar(min = 1, max = ncol(df)*nrow(df), initial = 1) 

k<-0
for (i in 1:nrow(df)) {
  for (h in 1:ncol(df)) {
    k<-k+1
    
    if(df[i,h]!="X"){
    
    setTxtProgressBar(pb,k)
    df3D[i,h,1:6]<-as.numeric(as.character(aminoacids[aminoacids$V1==df[i,h],2:ncol(aminoacids)]))
    }
  }
  
}

#Keras model 

#Input

Input

#Parameters
n.filters <- 32
dim.l1 <- 500
dim.l2 <- 100
dim.l3 <- 20
dim.ls <- 2

model <-keras_model_sequential()

model %>% 
  layer_conv_2d(filters = n.filters, kernel_size = c(6,6),  strides = c(1, 1), padding = "same", activation = "tanh") %>% 
  layer_conv_2d(filters = n.filters*2, kernel_size = c(6,6),  strides = c(1, 1), padding = "same", activation = "tanh") %>% 
  layer_conv_2d(filters = n.filters*4, kernel_size = c(6,6),  strides = c(1, 1), padding = "same", activation = "tanh") %>% 
  layer_flatten() %>%
  layer_dense(units = dim.l1, activation = "tanh") %>% 
  layer_dense(units = dim.l2, activation = "tanh") %>%
  layer_dense(units = dim.l3, activation = "tanh") %>%
  layer_dense(units = dim.ls, activation = "tanh") %>%
  layer_dense(units = dim.l3, activation = "tanh") %>%
  layer_dense(units = dim.l2, activation = "tanh") %>%
  layer_dense(units = dim.l1, activation = "tanh") %>%
  layer_reshape(target_shape = c(150,3,n.filters*8)) %>%
  layer_conv_2d_transpose(filters = n.filters, kernel_size = c(6,6),  strides = c(1, 1), padding = "same", activation = "tanh") %>% 
  layer_conv_2d_transpose(filters = n.filters*2, kernel_size = c(6,6),  strides = c(1, 1), padding = "same", activation = "tanh") %>% 
  layer_conv_2d_transpose(filters = n.filters*4, kernel_size = c(6,6),  strides = c(1, 1), padding = "same", activation = "tanh") 
  
  

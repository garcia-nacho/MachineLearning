#### Computer vision in R

library(fitdistrplus)
library(jpeg)
library(countcolors)
library(keras)

df<-read.csv("C:/Users/AutophagyCrusher/Desktop/CellCV/Results.csv")
df$image<-gsub("C.-","",df$Label)
df$image<-gsub(":.*",".jpg",df$image)
df$time<-gsub(".*t:","",df$Label)
df$time<-gsub("C.*","12/12",df$time)
df$time<-gsub("/.*","",df$time)
df$image<-paste("T",df$time,df$image,sep = "")

path<-"C:/User/Autophagy/Cruser"
w<-1584 #width
h<-1584 #height

img<-readJPEG("C:/Users/AutophagyCrusher/Desktop/CellCV/PICS/T10TimeCourseNoise.lsm - TimeCourseNoise #1.jpg")

df.image<-subset(df,df$image=="T10TimeCourseNoise.lsm - TimeCourseNoise #1.jpg")

df.image<-df.image[,c(7:12)]
df.image$class<-1
df.noise<-df.image

df.noise$class<-0 
df.noise$X<-runif(nrow(df.image), min=31, max = w-31)
df.noise$Y<-runif(nrow(df.image), min=31, max = h-31)


df.total<-rbind(df.image,df.noise)




#Extract cells for training
plot(density(c(df.total$Width,df.total$Height))) #60x60
sizeX<-60
sizeY<-60
#Training array

training<-array(data = 0, dim = c(nrow(df.total),sizeY+1,sizeX+1,3))

pb<-txtProgressBar(min = 1, max = nrow(df.total), initial = 1)
cellID<-rep("KO",nrow(df.total))

for (i in 1:nrow(df.total)) {
  setTxtProgressBar(pb,i)
  
  if(df.total$Width[i]<sizeX+1 &
     df.total$Height[i]<sizeY+1 &
     df.total$X[i]>30 &
     df.total$Y[i]>30 &
     df.total$Y[i]<h-30 &
     df.total$X[i]<w-30){
    
    
  cellID[i]<-"OK"

  #Bounding boxes
  center<-vector()
  center[1]<-df.total$X[i]
  center[2]<-df.total$Y[i]


  dummyimage<-img[c(center[2]-30):(center[2]+30) , c(center[1]-30):(center[1]+30) , ]
  dummyimage<-array(dummyimage, dim=dim(dummyimage))}
  
  #plotArrayAsImage(dummyimage)
  

  
  training[i,,, ]<-dummyimage[,,]
  
  }




#Model 
model <-keras_model_sequential()

model %>%
  #
  layer_conv_2d(filter=32,kernel_size=c(3,3),padding="same",input_shape=c(sizeY+1,sizeX+1,3) ) %>%  
  layer_activation("relu") %>%  
  layer_conv_2d(filter=32 ,kernel_size=c(3,3))  %>%  
  layer_activation("relu") %>%
  layer_max_pooling_2d(pool_size = c(2,2)) %>% 
  
  layer_flatten() %>% 
  layer_dense(1024) %>%
  layer_activation("relu") %>% 
  layer_dense(128) %>% 
  layer_activation("relu") %>% 
  layer_dropout(0.3) %>% 
  layer_dense(40) %>% 
  layer_dropout(0.2) %>%
  layer_dense(2) %>% 
  layer_activation("softmax")

summary(model)

opt<-optimizer_adam( lr= 0.0001 , decay = 1e-6 )

compile(model,optimizer = opt, loss = 'categorical_crossentropy' )

#Ytrain
training<-training[which(cellID=="OK"),,,,drop=FALSE]
y_train<-df.total$class[which(cellID=="OK")]
y_train<-to_categorical(y_train)

date<-as.character(date())
logs<-gsub(" ","_",date)
logs<-gsub(":",".",logs)
logs<-paste("logs/",logs,sep = "")



history<-model %>% fit(training, y_train,
                       batch_size=10,
                       epoch=30,
                       validation_spit=0.1,
                       callbacks = callback_tensorboard(logs),
                       
                       view_metrics=TRUE,
                       shuffle=TRUE)

tensorboard(logs)


a<-predict(model,training)

###NO padding




# Face recognition 
# Nacho garcia 2019
# garcia.nacho@gmail.com

library(imager)
library(keras)

Path<-"/home/nacho/Downloads/faces/"

files<-list.files(Path)
setwd(Path)

df.temp <- load.image(files[1])
df.temp<-as.matrix(df.temp)

df<-array(data = 0, dim = c(length(files),nrow(df.temp),ncol(df.temp),1))

for (i in 1:length(files)) {
  df.temp <- load.image(files[i])
  df.temp<-as.matrix(df.temp)
  df[i,,,1]<-df.temp
   }

Y <- gsub(" .*","",files)
Y<- gsub("s","",Y)

  image(df[3,,,],
        useRaster = TRUE,
        axes=FALSE,
        col = gray.colors(256, start = 0, end = 1, gamma = 2.2, alpha = NULL))
  
  #dim(df)

#Model 
model <-keras_model_sequential()

model %>%
  #
  layer_conv_2d(filter=32,kernel_size=c(3,3),padding="same",input_shape=c(92,112,1) ) %>%  
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
  layer_activation("softmax")

summary(model)

opt<-optimizer_adam( lr= 0.0001 , decay = 1e-6 )

#Training and test
TestSeq<-(c(1:40)*10)+round(runif(1, min = 1, max = 10))-10

x_test <- df[TestSeq,,,1]
y_test <- Y[TestSeq]  
x_train <- df[-TestSeq,,,1]
y_train <- Y[-TestSeq]

x_test <- array(x_test, dim = c(dim(x_test),1))
x_train <- array(x_train, dim = c(dim(x_train),1))

y_train <- to_categorical(y_train)
y_test <- to_categorical(y_test)

y_train<-y_train[,2:41]
y_test<-y_test[,2:41]

#image(x_test[11,,,], useRaster = TRUE, axes=FALSE)

#Model fit

compile(model,optimizer = opt, loss = 'categorical_crossentropy' )

history<-model %>% fit(x_train, y_train,
              batch_size=10,
              epoch=10,
              validation_data = list(x_test, y_test),
              callbacks = callback_tensorboard("logs/run_a"),
              view_metrics=TRUE,
              shuffle=TRUE)

model %>% predict_classes(x_test)
y_test_pred <- Y[TestSeq]  

tensorboard("logs/run_a")

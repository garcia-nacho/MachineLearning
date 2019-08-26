#Sheep drawer
#Ignacio Garcia 2019

# library(reticulate)
# np <- import("numpy")
# 
# 
# df<-np$load("/home/nacho/VAE_Faces/full_numpy_bitmap_sheep.npy")
# df<-array(df, dim = c(dim(df)[1], 28,28))
# df<-df/255
# write.csv(df, "/home/nacho/VAE_Faces/datasheep.csv")
# 
# detach("package:reticulate", unload=TRUE)

library(keras)
library(imager)
library(ggplot2)
library(R6)

df<-read.csv("/home/nacho/VAE_Faces/datasheep.csv")
df<-as.matrix(df)
df<-array(as.numeric(df[,2:785]), dim = c(nrow(df),28,28,1))

image(df[10,,,1],
      useRaster = TRUE,
      axes=FALSE,
      col = gray.colors(256, start = 0, end = 1, gamma = 2.2, alpha = NULL))


#Model
reset_states(vae)
filters <- 10
intermediate_dim<-100
latent_dim<-2
epsilon_std <- 1
batch_size <- 8
epoch<-15
activation<-"relu"

dimensions<-dim(df)
dimensions<-dimensions[-1]

Input <- layer_input(shape = dimensions)

faces<- Input %>%
  layer_conv_2d(filters=filters, kernel_size=c(4,4), activation=activation, padding='same',strides=c(1,1),data_format='channels_last')%>% 
  layer_conv_2d(filters=filters*2, kernel_size=c(4,4), activation=activation, padding='same',strides=c(1,1),data_format='channels_last')%>%
  layer_conv_2d(filters=filters*4, kernel_size=c(4,4), activation=activation, padding='same',strides=c(1,1),data_format='channels_last')%>%
  layer_conv_2d(filters=filters*8, kernel_size=c(4,4), activation=activation, padding='same',strides=c(1,1),data_format='channels_last')%>%
  layer_flatten()

hidden <- faces %>% layer_dense( units = intermediate_dim, activation = activation) %>% 
  layer_dropout(0.1) %>% 
  layer_batch_normalization() %>% 
  layer_dense( units = round(intermediate_dim/2), activation = activation) %>% 
  layer_dropout(0.1) %>% 
  layer_batch_normalization() %>% 
  layer_dense( units = round(intermediate_dim/4), activation = activation)

z_mean <- hidden %>% layer_dense( units = latent_dim)
z_log_var <- hidden %>% layer_dense( units = latent_dim)


sampling <- function(args) {
  z_mean <- args[, 1:(latent_dim)]
  z_log_var <- args[, (latent_dim + 1):(2 * latent_dim)]
  
  epsilon <- k_random_normal(
    shape = c(k_shape(z_mean)[[1]]),
    mean = 0.,
    stddev = epsilon_std
  )
  z_mean + k_exp(z_log_var) * epsilon
}

z <- layer_concatenate(list(z_mean, z_log_var)) %>% layer_lambda(sampling, name="LatentSpace")

#Initialization of layers
Output1<- layer_dense( units = round(intermediate_dim/4), activation = activation)
Output2<- layer_dropout(rate=0.1)
Output3<- layer_batch_normalization()
Output4<- layer_dense( units = round(intermediate_dim/2), activation = activation)
Output5<- layer_dense(units = intermediate_dim, activation = activation)
Output6<- layer_dense(units = prod(28,28,filters*8), activation = activation)
Output7<- layer_reshape(target_shape = c(28,28,filters*8))
Output8<- layer_conv_2d_transpose(filters=filters*8, kernel_size=c(4,4), activation=activation, padding='same',strides=c(1,1),data_format='channels_last')
Output9<- layer_conv_2d_transpose(filters=filters*4, kernel_size=c(4,4), activation=activation, padding='same',strides=c(1,1),data_format='channels_last')
Output10<-layer_conv_2d_transpose(filters=filters*2, kernel_size=c(4,4), activation=activation, padding='same',strides=c(1,1),data_format='channels_last')
Output11<-layer_conv_2d_transpose(filters=filters, kernel_size=c(4,4), activation=activation, padding='same',strides=c(2,2),data_format='channels_last') 
Output12<-layer_conv_2d(filters=1, kernel_size=c(4,4), activation="sigmoid", padding='same',strides=c(2,2),data_format='channels_last') 

#Concatenation of layers for the VAE
O1<-Output1(z)
O2<-Output2(O1)
O3<-Output3(O2)
O4<-Output4(O3)
O5<-Output5(O4)
O6<-Output6(O5)
O7<-Output7(O6)
O8<-Output8(O7)
O9<-Output9(O8)
O10<-Output10(O9)
O11<-Output11(O10)
O12<-Output12(O11)

#Concatenation of layers for the Decoder
decoder_input <- layer_input(shape = latent_dim)
O1D<-Output1(decoder_input)
O2D<-Output2(O1D)
O3D<-Output3(O2D)
O4D<-Output4(O3D)
O5D<-Output5(O4D)
O6D<-Output6(O5D)
O7D<-Output7(O6D)
O8D<-Output8(O7D)
O9D<-Output9(O8D)
O10D<-Output10(O9D)
O11D<-Output11(O10D)
O12D<-Output12(O11D)

## variational autoencoder
vae <- keras_model(Input, O12)
summary(vae)

## Decoder
decoder<- keras_model(decoder_input, O12D)
summary(decoder)

#KL-Annealing

weight <- k_variable(0)
epochs <-30
kl.start <-10
kl.steep <- 10


# custom loss function
loss_vae<- function(weight){
  l.f<-function(x, x_decoded_mean){
    
    x <- k_flatten(x)
    x_decoded_mean <- k_flatten(x_decoded_mean)
    xent_loss <- 2.0 * dimensions[1]* dimensions[2]*loss_mean_squared_error(x, x_decoded_mean)
    kl_loss <- -0.5*k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1L)
    xent_loss + weight*kl_loss}
  return(l.f)
}


#Custom callback KL-Annealing

KL.Ann <- R6::R6Class("KL.Ann",
                      inherit = KerasCallback,
                      
                      public = list(
                        
                        losses = NULL,
                        params = NULL,
                        model = NULL,
                        weight = NULL,
                        
                        set_context = function(params = NULL, model = NULL) {
                          self$params <- params
                          self$model <- model
                          self$weight<-weight
                        },                        
                        
                        on_epoch_end = function(epoch, logs = NULL) {
                          
                          if(epoch>kl.start){
                            new_weight<- min((epoch-kl.start)/kl.steep,1)
                            k_set_value(self$weight, new_weight)
                            print(paste("     ANNEALING KLD:", k_get_value(self$weight), sep = " "))
                            
                          }
                        }
                        
                      ))




vae %>% compile(optimizer = "rmsprop", loss = loss_vae(weight))


date<-as.character(date())
logs<-gsub(" ","_",date)
logs<-gsub(":",".",logs)
logs<-paste("logs/",logs,sep = "")

#Callbcks
tb<-callback_tensorboard(logs)
annealing <- KL.Ann$new()

history<-vae %>% fit(x= df[1:10000,,,,drop=FALSE],
                y= df[1:10000,,,,drop=FALSE],
                batch_size=batch_size,
                epoch=epochs,
                callbacks = list(annealing,tb),
                #initial_epoch=30,
                view_metrics=FALSE,
                shuffle=TRUE)

tensorboard(logs)

load_model_weights_hdf5(vae, "/home/nacho/VAE_Faces/Baeh.h5")
#save_model_weights_hdf5(vae, "/home/nacho/VAE_Faces/Baeh.h5")

Samp1.gen <- predict(vae, df[24,,,,drop=FALSE], batch_size = 10)
Samp1.gen<-as.matrix(Samp1.gen)
Samp1.gen<-matrix(Samp1.gen, ncol = 28)
sharp.kern<- matrix(-1/9, ncol = 3, nrow = 3)
sharp.kern[2,2]<--1

library(OpenImageR)

Samp1.gen <- convolution(Samp1.gen, sharp.kern, mode="same")


image(Samp1.gen,
      useRaster = TRUE,
      axes=FALSE, 
      col = gray.colors(256, start = 0, end = 1, gamma = 2.2, alpha = NULL))

layer_name <- 'LatentSpace'
intermediate_layer_model <- keras_model(inputs = vae$input,
                                        outputs = get_layer(vae, layer_name)$output)
intermediate_output <- as.data.frame(predict(intermediate_layer_model, df))
plot(intermediate_output$V1,intermediate_output$V2,xlim = c(-5,5), ylim=c(-5,5))

#write.csv(intermediate_output,"/home/nacho/VAE_Faces/BaehhhLS.csv")

ggplot(intermediate_output)+
  geom_point(aes(V1,V2), colour="blue", alpha=0.05)+
  xlim(-2.5,2.5)+
  ylim(-2.5,2.5)+
  theme_minimal()

to.remove<-which(is.infinite(intermediate_output[,1]))

clusters<-kmeans(intermediate_output[-to.remove,],centers = 196)

intermediate_output$Cluster<-NA
intermediate_output$Cluster[-to.remove]<-clusters$cluster

withnn<-vector()
for (i in 1:250) {
  dummy<-kmeans(intermediate_output[-to.remove,],centers = 5)
  withnn[i]<-mean(dummy$withinss)
  
}

plot(withnn)

#DBSCAN

db<-as.data.frame(matrix(data = NA, ncol = 3, nrow = 100))

pb<-txtProgressBar(min = 1, max = 100, initial = 1)
for (i in 1:100) {
  setTxtProgressBar(pb,i)
  eps<-runif(1, min = 0.5, max =0.5 )
  mp<-round(runif(1,min = 50, max = 1000))
  dummy<-dbscan(intermediate_output[-to.remove,1:2], eps = eps, MinPts = mp)
  db[i,1]<-eps
  db[i,2]<-mp
  db[i,3]<-max(dummy$cluster)
}

library(dbscan)
db.s<-dbscan(intermediate_output[-to.remove,1:2], eps = 0.04270631, MinPts = 377)
intermediate_output[-to.remove,]$Cluster<-db.s$cluster

ggplot(intermediate_output)+
  geom_point(aes(V1,V2,colour=as.character(Cluster)),  alpha=0.5)+
  xlim(-2.5,2.5)+
  ylim(-2.5,2.5)+
  theme_minimal()

cluster6<-which(intermediate_output$Cluster==6)
cluster8<-which(intermediate_output$Cluster==8)

sample6<-df[sample(cluster6,1),,,]
sample8<-df[sample(cluster8,1),,,]

cbind(abs(sample6-1),abs(sample8-1)) %>% as.raster() %>% plot()

#Annomaly detection 
MSE<-vector()
pb<-txtProgressBar(min = 1, max = nrow(df), initial = 1)
for (i in 1:nrow(df)) {
  setTxtProgressBar(pb,i)
  dummy<-predict(vae, df[i,,,,drop=FALSE], batch_size = 10)
  MSE[i] <- sum((df[i,,,] - dummy[1,,,])^2) 
}

plot(MSE)
abline(h =300, col="red", lwd=3, lty=2)

anomalies<-which(MSE>300)
an.plot<-df[anomalies[1],,,]
for (i in 2:length(anomalies)) {
  an.plot<-cbind(an.plot, df[anomalies[i],,,])
}

abs(an.plot-1) %>% as.raster() %>% plot()

# generator, from latent space to reconstructed inputs

Samp<-array(data = c(rnorm(1, mean = 1, sd=0.5),rnorm(1, mean = 1, sd=0.5)), dim = c(1,2))
Samp<- predict(decoder, generator.arrray, batch_size = 10)
Samp <- Samp[1,,,1]

for (i in 1:2) {
  generator.arrray<-array(data = c(rnorm(1, mean = 1, sd=0.5),rnorm(1, mean = 1, sd=0.5)), dim = c(1,2))
  Gen1 <- predict(decoder, generator.arrray, batch_size = 10)
  Gen1 <- Gen1[1,,,1]
  Samp<-cbind(Samp,Gen1)
  
}
abs(Samp-1) %>% as.raster() %>% plot()

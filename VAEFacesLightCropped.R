#VAE faces
#Nacho Garcia 2019
#garcia.nacho@gmail.com

#Library loading
library(keras)
library(imager)
library(matlab)

#File loading


Path<-"/home/nacho/VAE_Faces/ImagesPGM/"

files<-list.files(Path)
setwd(Path)

#Testing purposes 10%
files<-sample(files, round(0.2*length(files)), replace = FALSE)

df.temp <- load.image(files[1])
df.temp<-as.matrix(df.temp)

df<-array(data = 0, dim = c(length(files),nrow(df.temp),ncol(df.temp),1))

for (i in 1:length(files)) {
  df.temp <- load.image(files[i])
  df.temp<-as.matrix(df.temp)
  #Padding smaller pics
  if(ncol(df.temp)<dim(df)[3]) df.temp<-padarray(df.temp,c(0,1+round((dim(df)[3]-ncol(df.temp))/2)))
  if(nrow(df.temp)<dim(df)[2]) df.temp<-padarray(df.temp,c(1,1+round((dim(df)[2]-nrow(df.temp))/2)))
  
  #Cropping bigger pics
  if(ncol(df.temp)>dim(df)[3]) df.temp<-df.temp[,1:dim(df)[2]]
  if(nrow(df.temp)>dim(df)[2]) df.temp<-df.temp[1:dim(df)[2],]
  
  df[i,,,1]<-df.temp
}

df2<- array(data = 0, dim = c(dim(df)[1],300,300,1))
df2[,1:299,1:299,]<-df
df<-df2
rm(df2)

#Image cropper based on data obtained
df<-as.array(df[,51:250,51:250,])
dim(df)<-c(dim(df)[1],200,200,1)

#Image visualization
image(df[1,,,],
      useRaster = TRUE,
      axes=FALSE,
      col = gray.colors(256, start = 0, end = 1, gamma = 2.2, alpha = NULL))

#Cleaning file names
files<-gsub("_.*","",files)
files<-gsub(".pgm","",files)
files<-gsub("face.*","",files)
files<-gsub("\\d","",files)


#Model
reset_states(vae)
filters <- 4
intermediate_dim<-64
latent_dim<-2
epsilon_std <- 1.0
batch_size <- 8
epoch<-30

dimensions<-dim(df)
dimensions<-dimensions[-1]

Input <- layer_input(shape = dimensions)

faces<- Input %>%
  layer_conv_2d(filters=1, kernel_size=c(5,5), activation='sigmoid', padding='same',strides=c(1,1),data_format='channels_last')%>% 
  layer_conv_2d(filters=filters*2, kernel_size=c(5,5), activation='sigmoid', padding='same',strides=c(2,2),data_format='channels_last')%>%
  layer_conv_2d(filters=filters*4, kernel_size=c(5,5), activation='sigmoid', padding='same',strides=c(1,1),data_format='channels_last')%>%
  layer_conv_2d(filters=filters*8, kernel_size=c(5,5), activation='sigmoid', padding='same',strides=c(2,2),data_format='channels_last')%>%
  layer_flatten()

hidden <- faces %>% layer_dense( units = intermediate_dim, activation = "sigmoid") %>% 
  layer_dense( units = round(intermediate_dim/2), activation = "sigmoid") %>% 
  layer_dense( units = round(intermediate_dim/4), activation = "sigmoid")

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

z <- layer_concatenate(list(z_mean, z_log_var)) %>% layer_lambda(sampling)

Output<- z %>%
  layer_dense( units = round(intermediate_dim/4), activation = "sigmoid") %>% 
  layer_dense( units = round(intermediate_dim/2), activation = "sigmoid") %>%
  layer_dense(units = intermediate_dim, activation = "sigmoid") %>%
  layer_dense(units = prod(50,50,filters*8), activation = "sigmoid") %>%
  layer_reshape(target_shape = c(50,50,filters*8)) %>%
  layer_conv_2d_transpose(filters=filters*8, kernel_size=c(5,5), activation='sigmoid', padding='same',strides=c(2,2),data_format='channels_last')%>%
  layer_conv_2d_transpose(filters=filters*4, kernel_size=c(5,5), activation='sigmoid', padding='same',strides=c(1,1),data_format='channels_last')%>%
  layer_conv_2d_transpose(filters=filters*2, kernel_size=c(5,5), activation='sigmoid', padding='same',strides=c(2,2),data_format='channels_last')%>%
  layer_conv_2d_transpose(filters=1, kernel_size=c(5,5), activation='sigmoid', padding='same',strides=c(1,1),data_format='channels_last')

# custom loss function
vae_loss <- function(x, x_decoded_mean_squash) {
  
  x <- k_flatten(x)
  x_decoded_mean_squash <- k_flatten(x_decoded_mean_squash)
  
  xent_loss <- 5.0 * dimensions[1]* dimensions[1] *
    loss_binary_crossentropy(x, x_decoded_mean_squash)
  
  kl_loss <- -0.5 * k_mean(1 + z_log_var - k_square(z_mean) -
                             k_exp(z_log_var), axis = -1L)
  
  k_mean(xent_loss + kl_loss)
}

## variational autoencoder
vae <- keras_model(Input, Output)
vae %>% compile(optimizer = "rmsprop", loss = vae_loss)
summary(vae)

## Split for Training 
reset_states(vae)


date<-as.character(date())
logs<-gsub(" ","_",date)
logs<-gsub(":",".",logs)
logs<-paste("logs/",logs,sep = "")

history<-vae %>% fit(x= df,
                     y=df,
                     batch_size=batch_size,
                     epoch=epoch,
                     callbacks = callback_tensorboard(logs),
                     view_metrics=FALSE,
                     shuffle=TRUE)

tensorboard(logs)


Samp1.gen <- predict(vae, array(df[15,,,], dim = c(1,200,200,1)), batch_size = 10)
Samp1.gen<-as.matrix(Samp1.gen)
Samp1.gen<-matrix(Samp1.gen, ncol = 200)

image(Samp1.gen,
      useRaster = TRUE,
      axes=FALSE, 
      col = gray.colors(256, start = 0, end = 1, gamma = 2.2, alpha = NULL))

#Latent space

layer_name <- 'lambda_3'
intermediate_layer_model <- keras_model(inputs = vae$input,
                                        outputs = get_layer(vae, layer_name)$output)
intermediate_output <- as.data.frame(predict(intermediate_layer_model, df))
plot(intermediate_output$V1,intermediate_output$V2)

#Looking for outliers to be removed, next pull 
for (i in 1:dim(df)[1]) {
  
}

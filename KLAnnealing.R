weight <- k_variable(0)
epochs <-30
kl.start <-2
kl.steep <- 4

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
                               print(paste("                      ANNEALING KLD:", k_get_value(self$weight), sep = " "))
                               
                               }
                               }
                             
                           ))


loss<- function(weight){
          l.f<-function(x, x_decoded_mean){
           xent_loss <- (original_dim/1.0)*loss_binary_crossentropy(x, x_decoded_mean)
           kl_loss <- -0.5*k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1L)
            xent_loss + weight*kl_loss}
          return(l.f)
}

vae %>% compile(optimizer = "rmsprop", loss = loss(weight) )

# Model training ----------------------------------------------------------

history <- KL.Ann$new()

vae %>% fit(
  x_train, x_train, 
  shuffle = TRUE, 
  epochs = 10, 
  batch_size = batch_size,
  callbacks = list(history),
  validation_data = list(x_test, x_test)
)

#loading a pre-trained framework called Xception
pretrained_framework <- application_xception(weights = 'imagenet', 
                                             include_top = FALSE, 
                                             input_shape = c(fish_width, fish_height, 3))
#freeze all weights from Xception
freeze_pretrained_framework <- freeze_weights(pretrained_framework)

#writing a layer on top of the pre-trained framework including setting parameters to tune
fish_model_function <- function(learning_rate = 0.001, dropout_rate = 0.2,
                                n_dense = 1024) {

fish_model <- keras_model_sequential() %>%
    pretrained_framework %>% 
    layer_global_average_pooling_2d() %>% 
    layer_dense(units = n_dense) %>%
    layer_activation("relu") %>%
    layer_dropout(dropout_rate) %>%
    layer_dense(units = fish_label_length, activation = "softmax")
  
fish_model %>% compile(
    loss = "categorical_crossentropy",
    optimizer = optimizer_adam(lr = learning_rate),
    metrics = "accuracy")

return(fish_model)
}

fish_model <- fish_model_function()

#training the model 
batch_size <- 32 #standard amount 
epochs <- 12

#create path to save checkpoints during and after training
checkpoint_path <- "C:/Users/rikac/OneDrive/Documents/R/IBC 2"
#create checkpoint callback
checkpoint_callback <- callback_model_checkpoint(
  filepath = checkpoint_path,
  save_weights_only = TRUE,
  verbose = 0
)

#val_acc: 0.8125 - latest trained
training_history1 <- fish_model %>% fit(
                          fish_train_images,
                          steps_per_epoch = fish_train_images$n %/% batch_size, 
                          epochs = 12, 
                          validation_data = fish_val_images,
                          validation_steps = fish_val_images$n %/% batch_size,
                          verbose = 2)

#evaluate on test set (accuracy: 0.75 - trained on training_history1)
test_eval <- fish_model %>% 
  evaluate(fish_test_images, steps = fish_test_images$n %/% batch_size)





















tune_grid <- data.frame("learning_rate" = c(0.01, 0.001, 0.0001, 0.00001),
                        "dropout_rate" = c(0.5, 0.6, 0.7, 0.8),
                        "n_dense" = c(512, 256, 128, 64))
tuning_results2 <- NULL

for (i in 1:length(tune_grid$learning_rate)){
  for (j in 1:length(tune_grid$dropout_rate)){
    for (k in 1:length(tune_grid$n_dense)){
    
      fish_model <- fish_model_function(  
        learning_rate = tune_grid$learning_rate[i],
        dropout_rate = tune_grid$dropout_rate[j],
        n_dense = tune_grid$n_dense[k])
      
      training_history <- fish_model %>% fit(
        fish_train_images,
        steps_per_epoch = fish_train_images$n %/% batch_size, 
        epochs = epochs, 
        validation_data = fish_val_images,
        validation_steps = fish_val_images$n %/% batch_size,
        verbose = 2
      )
      
      #Save model configurations
      tuning_results2 <- rbind(tuning_results2,
        c("learning_rate" = tune_grid$learning_rate[i],
          "dropout_rate" = tune_grid$dropout_rate[j],
          "n_dense" = tune_grid$n_dense[k],
          "val_accuracy" = training_history$metrics$val_accuracy))
    }
  }
}

#lr: 0.001, dr: 0.6, n_dense: 1024 (high val acc: 0.78125)
best_results2 <- tuning_results2[which( 
  tuning_results2[,ncol(tuning_results2)] == 
    max(tuning_results2[,ncol(tuning_results2)])),]

#lr: 0.001, dr: 0.7, n_dense: 1024 (high val acc: 0.78125)
best_results3 <- tuning_results3[which( 
  tuning_results3[,ncol(tuning_results3)] == 
    max(tuning_results3[,ncol(tuning_results3)])),]






















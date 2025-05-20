fish_model <- fish_model_function(
               learning_rate = best_results2["learning_rate"],
               dropout_rate = best_results2["dropout_rate"],
               n_dense = best_results2["n_dense"])

training_history1 <- fish_model %>% fit(
  fish_train_images,
  steps_per_epoch = fish_train_images$n %/% batch_size, 
  epochs = 12, 
  validation_data = fish_val_images,
  validation_steps = fish_val_images$n %/% batch_size,
  verbose = 2)

fish_model %>% save_model_tf("fish_model")



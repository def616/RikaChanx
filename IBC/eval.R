#checking to see which species the model predicts best
library(keras)
library(readr)
library(tensorflow)
library(stringr)
library(purrr)
library(base)
library(ggplot2)
library(imager)
library(tidyverse)

fishID <- load_model_tf("www/fish_model")

predictions <- fishID %>% 
  predict(fish_test_images,
    steps = fish_test_images$n
  ) %>% as.data.frame

names(predictions) <- paste0("Class", 0:13)

predictions$predicted_class <- 
  paste0("Class", apply(predictions, 1, which.max) - 1)
predictions$true_class <- paste0("Class",fish_test_images$classes)

predictions %>% group_by(true_class) %>% 
  summarise(percentage_true = 100*sum(predicted_class == 
                                        true_class)/n()) %>% 
  left_join(data.frame(fish = names(fish_test_images$class_indices), 
                       true_class=paste0("Class", 0:13)), by="true_class") %>%
  select(fish, percentage_true) %>% 
  mutate(fish = fct_reorder(fish, percentage_true)) %>%
  ggplot(aes(x=fish,y=percentage_true,fill=percentage_true, 
             label=percentage_true)) +
  geom_col() + theme_minimal() + coord_flip() +
  geom_text(nudge_y = 3) + 
  ggtitle("Percentage correct classifications by fish species")

#checking tuned model accuracy on test images
#train accuracy: 0.9797980
fishID %>% evaluate(fish_train_images, steps = fish_train_images$n)
#validation accuracy: 0.8181818
fishID %>% evaluate(fish_val_images, steps = fish_val_images$n)
#validation accuracy: 0.6507937 
fishID %>% evaluate(fish_test_images, steps = fish_test_images$n)




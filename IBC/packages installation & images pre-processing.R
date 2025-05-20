library(keras)
library(readr)
library(tensorflow)
library(stringr)
library(purrr)
library(base)
library(ggplot2)
library(imager)
library(tidyverse)

#install tensorflow
install_tensorflow()
install_tensorflow(version = "gpu")
install_keras(tensorflow = "gpu")
install_tensorflow(extra_packages="pillow") #python package

#photos directory 
working_dir <- setwd("C:/Users/rikac/OneDrive/Documents/R/IBC 2/For project - FINAL (Train & Test)")
fish_label <- dir("C:/Users/rikac/OneDrive/Documents/R/IBC 2/For project - FINAL (Train & Test)/TRAIN")
fish_label_length <- length(fish_label)
save(fish_label, file="fish_label.RData")

#setting images scaling (in pixels)
fish_width <- 256
fish_height <- 256
target_fish_size <- c(fish_width, fish_height)
fish_color <- 3 #RGB

#create train and validation sets
fish_dir <- "C:/Users/rikac/OneDrive/Documents/R/IBC 2/For project - FINAL (Train & Test)/TRAIN"
#generate fish images for train and validation set
fish_img_gen <- image_data_generator(rescale = 1/255, validation_split = 0.2)

#train set: 198 images
fish_train_images <- flow_images_from_directory(fish_dir,
                                                fish_img_gen,
                                                subset = 'training',
                                                target_size = target_fish_size,
                                                class_mode = "categorical",
                                                shuffle = F,
                                                classes = fish_label,
                                                seed = 2021)
#validation set: 44 images
fish_val_images <- flow_images_from_directory(fish_dir,
                                              fish_img_gen,
                                              subset = 'validation',
                                              target_size = target_fish_size,
                                              class_mode = "categorical",
                                              shuffle = F,
                                              classes = fish_label,
                                              seed = 2021)

#test set
fish_dir_test <- "C:/Users/rikac/OneDrive/Documents/R/IBC 2/For project - FINAL (Train & Test)/TEST"
fish_img_gen_test <- image_data_generator(rescale = 1/255)

#test set: 63 images
fish_test_images <- flow_images_from_directory(fish_dir_test,
                                               fish_img_gen_test,
                                               target_size = target_fish_size,
                                               class_mode = "categorical",
                                               shuffle = F,
                                               classes = fish_label,
                                               seed = 2021)

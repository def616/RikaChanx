# Fish Image-Based Classification

**Tools Used:** R, tensorflow, keras, Xception (pre-trained framework)
**Goal:** To classify 12 families of fish species based on their underwater photos

## Methods
- Preprocessed images
- Initated and fine-tuned a pre-trained model
- Created a simple Shiny app for UI
- Showed evaluation metrics

## Folder structure
- `IBC/`
  - `README.md`: description of the project
  - `images pre-processing`: pre-processing of obtained fish images
  - `model training & evaluating.R`: model training on train dataset and evaluation on train, validation, and test sets
  - `model tuning`: tuning the model for better accuracy
  - `model training with tuned parameters:` retraining the model with tuned parameters
  - `eval.R`: final evaluation of the tuned model on test dataset
  - `result.jpeg`: accuracy result across fish species
  - `app.R`: script for creating the app
  - `shinyapp.R`: script for launching the app

## Description
This project aims to create a simple UI interface and use a pre-trained framework to create a classification model of fish images. There were a total of 14 fish families, known by their common names. The fish photos were taken underwater under poor visibility; hence, the challenge of the project was to tune th model so that it would recognize poor quality picture.






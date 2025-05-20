library(keras)
library(tensorflow)
library(shiny)
library(rsconnect)
library(tidyverse)
library(shinydashboard)

fishID <- load_model_tf("www/fish_model")
load("www/fish_label.RData")
target_sizeID <- c(256, 256, 3)
options(scipen = 999)


library(shiny)

  ui <- dashboardPage(
    skin="blue",
    
    #(1) Header
    
    dashboardHeader(title = tags$h1("Fish Identification", 
                                    style = "font-size: 120%; font-weight: bold; color: white"),
                    titleWidth = 350,
                    tags$li(class = "dropdown"),
                    dropdownMenu(type = "notifications", icon = icon("question-circle", "fa-1x"), 
                                 badgeStatus = NULL,
                                 headerText="",
                                 tags$li(a(href = "https://forloopsandpiepkicks.wordpress.com",
                                           target = "_blank",
                                           tagAppendAttributes(icon("icon-circle"), class = "info"),
                                           "Created by"))
                    )),
    
    #(2) Sidebar
    
    dashboardSidebar(
      width = 350,
      fileInput("input_image", "File" , accept = c('.jpg','.jpeg')), 
      tags$br(),
      tags$p("Upload image")
    ),
    
    #(3) Body
    
    dashboardBody(
      h4("Instruction:"),
      tags$br(),tags$p("1. Take a picture of a fish species."),
      tags$p("2. Crop image so that the fish fills out most of the image."),
      tags$p("3. Upload image with menu on the left."),
      tags$br(),
      
      fluidRow(
        column(h4("Image:"),imageOutput("output_image"), width=6),
        column(h4("Result:"),tags$br(),textOutput("warntext",), tags$br(),
               tags$p("This species is probably a:"),tableOutput("text"),width=6)
      ),tags$br()
      
    ))

server <- function(input, output) {
 
    
    image <- reactive({image_load(input$input_image$datapath, 
                                  target_size = target_sizeID[1:2])})
    
    
    prediction <- reactive({
      if(is.null(input$input_image)){return(NULL)}
      x <- image_to_array(image())
      x <- array_reshape(x, c(1, dim(x)))
      x <- x/255
      pred <- fishID %>% predict(x)
      pred <- data.frame("Fish" = fish_label, "Prediction" = t(pred))
      pred <- pred[order(pred$Prediction, decreasing = TRUE),][1:5,]
      pred$Prediction <- paste(format(100 * pred$Prediction, 2), "%")
      pred
    })
    
    output$text <- renderTable({
      prediction()
    })
    
    output$warntext <- renderText({
      req(input$input_image)
      
      if(as.numeric(substr(prediction()[1,2],1,4)) >= 30){return(NULL)}
      warntext <- "Warning: I am not sure about this species!"
      warntext
    })
    
    
    output$output_image <- renderImage({
      req(input$input_image)
      
      outfile <- input$input_image$datapath
      contentType <- input$input_image$type
      list(src = outfile,
           contentType=contentType,
           width = 400)
    }, deleteFile = TRUE)
  
}

shinyApp(ui, server)
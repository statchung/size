library(shiny)
library(shinythemes)
library(fpow) 

fluidPage( tagList(
    
    navbarPage(
      # theme = "cerulean",  # <--- To use a theme, uncomment this
      "Design options",
      tabPanel("Factorial Design",
               sidebarPanel(
                 numericInput('nf','Number of Factor',2,min=2,max=10),
                 textInput("fl", 'Factor Levels', value = "2,2"),
                 checkboxGroupInput("checkGroup", "order", 
                                    choices = list("Main" = 1, "Two-way Interactions" = 2),
                                    selected = 1),
                  
                  radioButtons(
                      "delta_type", "Type of effect size", inline = TRUE, 
                      c("SD" = 1, "Range of effect" =2)
                    ),
                    conditionalPanel("input.delta_type==1", 
                                     
                                       numericInput('de1','SD(main)',1,min=0,max=10),
                                       numericInput('de2','SD(Interaction)',1,min=0,max=10),
                                       numericInput('de3','SD(Noise)',1,min=0,max=10)
                                      
                                     )  ,
                 conditionalPanel("input.delta_type==2", 
                                  
                                  numericInput('de11','Range of effect(main)',1,min=0,max=10),
                                  numericInput('de12','Range of effect(Interaction)',1,min=0,max=10),
                                  numericInput('de13','Range of effect(Noise)',1,min=0,max=10)
                                  
                 )  ,
                    
                   
                 numericInput('a','Type I error',0.05,min=0,max=1),
                 numericInput('b','Power at dectectable error',0.8,min=0,max=1)
                 ,
                 actionButton("do", "Click Me" , icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
               ),
                
               mainPanel(
                 tabsetPanel(
                   tabPanel("Result",
                            h4('Sample Size'),
                            verbatimTextOutput("Size1"),
                            h4('Dectable effect size'),
                            verbatimTextOutput("Size2")
                   ),
                   tabPanel("Delta vs Power Plot", 
                            sidebarLayout(      
                              
                              # Define the sidebar with one input
                              sidebarPanel(
                                selectInput("plot_order", "List of factors:", 
                                            choices=c(toupper(letters)[1:2])) 

                              ),
                              
                            mainPanel(plotOutput("Delta_graph")
                                      ,
                                      tableOutput("values")
                                      
                                      )) ) ,
                   tabPanel("Delta vs Size Plot",  plotOutput("Size_graph")) ,
                   tabPanel("Power vs Size Plot",
                            sidebarLayout(      
                              
                              # Define the sidebar with one input
                              sidebarPanel(
                                selectInput("plot_delta", "SD", 
                                            choices=c(1,1.5,2.0)) 
                                
                              ),
                              
                              mainPanel(plotOutput("power_graph")
                                        
                              ))) 
                 )
               )
      ),
      tabPanel("Fractional FD", "This panel is intentionally left blank"),
      tabPanel("Randomized Block", "This panel is intentionally left blank"),
      tabPanel("Split Design", "This panel is intentionally left blank")
    )
  ))

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
                 checkboxGroupInput("checkGroup", "Order", 
                                    choices = list("Main" = 1, "Two-way Interactions" = 2),
                                    selected = 1),
                  
                  radioButtons(
                      "delta_type", "Type of effect size", inline = TRUE, 
                      c("SD" = 1, "Range of effect" =2)
                    ),
                    conditionalPanel("input.delta_type==1", 
                                     
                                       numericInput('de1','SD(Main)',1,min=0,max=10),
                                       numericInput('de2','SD(Interaction)',1,min=0,max=10),
                                       numericInput('de3','SD(Noise)',1,min=0,max=10)
                                      
                                     )  ,
                 conditionalPanel("input.delta_type==2", 
                                  
                                  numericInput('de11','Range of effect(Main)',1,min=0,max=10),
                                  numericInput('de12','Range of effect(Interaction)',1,min=0,max=10),
                                  numericInput('de13','SD(Noise)',1,min=0,max=10)
                                  
                 )  ,
                    
                   
                 numericInput('a','Type I error',0.05,min=0,max=1),
                 numericInput('b','Power',0.8,min=0,max=1)
                 ,
                 actionButton("do", "Click Me" , icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
               ),
                
               mainPanel(
                 tabsetPanel(
                   tabPanel("Result",
                            h4('Model'),
                            verbatimTextOutput("list1" ),
                            h4('Sample Size'),
                            verbatimTextOutput("Size1"),
                            h4('Detectable effect size'),
                            verbatimTextOutput("Size2")
                            
                   ),
                   tabPanel("Delta vs Power Plot", 
                            sidebarLayout(      
                              
                              # Define the sidebar with one input
                              sidebarPanel(
                                selectInput("plot_order", "List of factors:", 
                                            choices=c(toupper(letters)[1:2],"ALL")) 

                              ),
                              
                            mainPanel(plotOutput("Delta_graph")
                                      ,
                                      tableOutput("values")
                                      
                                      )) ) ,
                   tabPanel("Size vs Delta Plot",  plotOutput("Size_graph")) ,
                   tabPanel("Size vs Power Plot",
                            sidebarLayout(      
                              
                              # Define the sidebar with one input
                              sidebarPanel(
                                selectInput("plot_delta", "Delta", 
                                            choices=c(1,1.5,2.0)) 
                                
                              ),
                              
                              mainPanel(plotOutput("power_graph")
                                        
                              ))) 
                 )
               )
      ),
      tabPanel("Fractional FD", sidebarPanel(
        numericInput('nf2','Number of 2 level Factor',2,min=2,max=10),
        numericInput('fr2','Fraction p (eg:2^(k-p))',1,min=1,max=10),
     
        radioButtons(
          "check_res", "Resolution", inline = TRUE, 
          c("III" = 3, "Higher resolution" =4)
        ),
        
        radioButtons(
          "delta_type2", "Type of effect size", inline = TRUE, 
          c("SD" = 1, "Range of effect" =2)
        ),
        conditionalPanel("input.delta_type2==1", 
                         
                         numericInput('de1_2','SD(Main)',1,min=0,max=10),
                        
                         numericInput('de3_2','SD(Noise)',1,min=0,max=10)
                         
        )  ,
        conditionalPanel("input.delta_type2==2", 
                         
                         numericInput('de11_2','Range of effect(Main)',1,min=0,max=10),
                         
                         numericInput('de13_2','SD(Noise)',1,min=0,max=10)
                         
        )  ,
        
        
        numericInput('a2','Type I error',0.05,min=0,max=1),
        numericInput('b2','Power',0.8,min=0,max=1)
        ,
        actionButton("do2", "Click Me" , icon("paper-plane"), 
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel("Result",
                   h4('Model'),
                   verbatimTextOutput("list1_2" ),
                   h4('Sample Size'),
                   verbatimTextOutput("Size12"),
                   h4('Detectable effect size'),
                   verbatimTextOutput("Size22")
          ),
          tabPanel("Delta vs Power Plot", 
                   
                       plotOutput("Delta_graph2")
                             ,
                               tableOutput("values2"))
                               
                      ,
          tabPanel("Delta vs Size Plot",  plotOutput("Size_graph2")) ,
          tabPanel("Power vs Size Plot",
                   sidebarLayout(      
                     
                     # Define the sidebar with one input
                     sidebarPanel(
                       selectInput("plot_delta2", "Delta", 
                                   choices=c(1,1.5,2.0)) 
                       
                     ),
                     
                     mainPanel(plotOutput("power_graph2")
                               
                     ))) 
        )
      )),
      tabPanel("Randomized Block",
               sidebarPanel(
                 numericInput('nf3','Number of Factor',2,min=2,max=10),
                 textInput("fl3", 'Factor Levels', value = "2,2"),
                 checkboxGroupInput("checkGroup3", "Order", 
                                    choices = list("Main" = 1, "Two-way Interactions" = 2),
                                    selected = 1),
                 
                 radioButtons(
                   "delta_type3", "Type of effect size", inline = TRUE, 
                   c("SD" = 1, "Range of effect" =2)
                 ),
                 conditionalPanel("input.delta_type3==1", 
                                  
                                  numericInput('de1_3','SD(Main)',1,min=0,max=10),
                                  numericInput('de2_3','SD(Interaction)',1,min=0,max=10),
                                  numericInput('de3_3','SD(Noise)',1,min=0,max=10)
                                  
                 )  ,
                 conditionalPanel("input.delta_type3==2", 
                                  
                                  numericInput('de11_3','Range of effect(Main)',1,min=0,max=10),
                                  numericInput('de12_3','Range of effect(Interaction)',1,min=0,max=10),
                                  numericInput('de13_3','SD(Noise)',1,min=0,max=10)
                                  
                 )  ,
                 
                 
                 numericInput('a3','Type I error',0.05,min=0,max=1),
                 numericInput('b3','Power',0.8,min=0,max=1)
                 ,
                 actionButton("do3", "Click Me" , icon("paper-plane"), 
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("Result",
                            
                            h4('Model'),
                            verbatimTextOutput("list1_3" ),
                            h4('Sample Size'),
                            verbatimTextOutput("Size1_3"),
                            h4('Detectable effect size'),
                            verbatimTextOutput("Size2_3")
                   ),
                   tabPanel("Delta vs Power Plot", 
                            sidebarLayout(      
                              
                              # Define the sidebar with one input
                              sidebarPanel(
                                selectInput("plot_order3", "List of factors:", 
                                            choices=c(toupper(letters)[1:2],"ALL")) 
                                
                              ),
                              
                              mainPanel(plotOutput("Delta_graph3")
                                        ,
                                        tableOutput("values3")
                                        
                              )) ) ,
                   tabPanel("Size vs Delta Plot",  plotOutput("Size_graph3")) ,
                   tabPanel("Size vs Power Plot",
                            sidebarLayout(      
                              
                              # Define the sidebar with one input
                              sidebarPanel(
                                selectInput("plot_delta3", "Delta", 
                                            choices=c(1,1.5,2.0)) 
                                
                              ),
                              
                              mainPanel(plotOutput("power_graph3")
                                        
                              ))) 
                 )
               )
      ),
      tabPanel("Split-Plot Design",
               sidebarPanel(
                 numericInput('wf','Number of Whole Factor',2,min=2,max=10),
                 textInput("wfl", 'Whole Factor Levels', value = "2,2"),
                 numericInput('sf','Number of Split Factor',2,min=2,max=10),
                 textInput("sfl", 'Split Factor Levels', value = "2,2"),
                 checkboxGroupInput("checkGroup4", "Order", 
                                    choices = list("Main" = 1, "Two-way Interactions" = 2),
                                    selected = 1),
                 
                 radioButtons(
                   "delta_type4", "Type of effect size", inline = TRUE, 
                   c("SD" = 1, "Range of effect" =2)
                 ),
                 conditionalPanel("input.delta_type4==1", 
                                  
                                  numericInput('de1_4','SD(Main)',1,min=0,max=10),
                                  numericInput('de2_4','SD(Interaction)',1,min=0,max=10),
                                  numericInput('de3_4','SD(Whole Noise)',1,min=0,max=10),
                                  numericInput('de4_4','SD(Noise)',1,min=0,max=10)
                                  
                 )  ,
                 conditionalPanel("input.delta_type4==2", 
                                  
                                  numericInput('de11_4','Range of effect(Main)',1,min=0,max=10),
                                  numericInput('de12_4','Range of effect(Interaction)',1,min=0,max=10),
                                  numericInput('de13_4','SD(Whole Noise)',1,min=0,max=10),
                                  numericInput('de14_4','SD(Noise)',1,min=0,max=10)
                                  
                 )  ,
                 
                 
                 numericInput('a4','Type I error',0.05,min=0,max=1),
                 numericInput('b4','Power',0.8,min=0,max=1)
                 ,
                 actionButton("do4", "Click Me" , icon("paper-plane"), 
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("Result",
                            
                            h4('Model'),
                            verbatimTextOutput("list1_4" ),
                            h4('Sample Size'),
                            verbatimTextOutput("Size1_4"),
                            h4('Detectable effect size'),
                            verbatimTextOutput("Size2_4")
                   ),
                   tabPanel("Delta vs Power Plot", 
                            sidebarLayout(      
                              
                              # Define the sidebar with one input
                              sidebarPanel(
                                selectInput("plot_order4", "List of factors:", 
                                            choices=c(toupper(letters)[1:2],"ALL")) 
                                
                              ),
                              
                              mainPanel(plotOutput("Delta_graph4")
                                        ,
                                        tableOutput("values4")
                                        
                              )) ) ,
                   tabPanel("Size vs Delta Plot",  plotOutput("Size_graph4")) ,
                   tabPanel("Size vs Power Plot",
                            sidebarLayout(      
                              
                              # Define the sidebar with one input
                              sidebarPanel(
                                selectInput("plot_delta4", "Delta", 
                                            choices=c(1,1.5,2.0)) 
                                
                              ),
                              
                              mainPanel(plotOutput("power_graph4")
                                        
                              ))) 
                 )
               )
      )
    )
  ))

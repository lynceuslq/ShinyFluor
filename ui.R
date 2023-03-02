library(shiny)
library(shinyBS)
library(shinycssloaders)
library(shinydashboard)
library(dashboardthemes)
library(shinyjs)
library(shinyWidgets)
library(ggpubr)
library(survminer)
library(tidyverse)
library(rlang)
library(plyr)
library(tippy)
library(plotly)
library(kableExtra)


js <- '.nav-tabs-custom .nav-tabs li.active {
    border-top-color: #fac002;
}"'

ui <-    dashboardPage(title = "ShinyFluor",
                       skin = "yellow",
                       dashboardHeader(title = "ShinyFlor", titleWidth = 200),
                       dashboardSidebar(
                           collapsed = T,
                           width = 200,
                           sidebarMenu(
                               menuItem("Florensent analyses", icon = icon("bar-chart-o"), tabName = "menu_top"),
                               menuItem("Characterised gates", icon = icon("table"), tabName = "menu_tab"),
                               menuItem("Github", icon = icon("github"), href = "https://github.com/lynceuslq"),
                               menuItem("About", icon = icon("question-circle-o"), tabName = "menu_about")
                           )
                       ),
                       dashboardBody(
                           tags$style(js),
                           tabItems(
                               tabItem(tabName = "menu_top",
                                       fluidRow(
                                           box(title = "Parameter", width = 4, solidHeader = T, status = "primary",height = "400px",
                                               fluidRow(
                                                   column(3,
                                                          p(HTML("<b>Gate model</b>"),span(shiny::icon("info-circle"), id = "info_uu"),
                                                            selectInput(
                                                                "gatemodel",
                                                                NULL,
                                                                choices = c(
                                                                    "Sensor" =1,
                                                                    "NOT gate" =2,
                                                                    "AND gate" =3
                                                                    # "NAND gate"=4
                                                                ),
                                                                selected=NULL

                                                            ),
                                                            tippy::tippy_this(elementId = "info_uu",tooltip = "Select a gate model for characterisation",placement = "right")
                                                          )
                                                   ),

                                                column(3,
                                                       p(HTML("<b>Steady state equation</b>"),span(shiny::icon("info-circle"), id = "info_se"),
                                                       tippy::tippy_this(elementId = "info_se",tooltip = "Output signals estimated from input signals at the steady state",placement = "right")
                                                       ),
                                                       uiOutput("gatefunc")

                                                )
                                               ),
                                               fluidRow(
                                                   column(6, fileInput("file1", "Upload florenscent results", accept = ".csv"), ),
                                                   column(3,
                                                          p(HTML("<b>Test Period</b>"),span(shiny::icon("info-circle"), id = "info_test_period"),
                                                            uiOutput("nameControl"),
                                                            tippy::tippy_this(elementId = "info_test_period",tooltip = "Please name your gate element to store it to the record'",placement = "right")
                                                          )
                                                   )
                                                 #  column(3, numericInput('power', "power = 1 - β", 0.80, min = 0.80, max = 0.99, step = 0.01))
                                               ),
                                               uiOutput("andControl"),
                                               fluidRow(
                                                   column(6, actionBttn("add", "Add record",color ="primary", style = "jelly")),
                                                   column(6, actionBttn("addto", "Add to resource",color ="primary", style = "jelly"))
                                               )
                                           ),
                                           box(title = "Summary of florenscent results",
                                               width = 4,
                                               solidHeader = T,
                                               status = "primary",
                                               DT::dataTableOutput("flo_tab"),
                                               height = "400px"
                                           ),
                                           uiOutput("moreControls"),

                                           box(title = "Table", width = 12, solidHeader = T, status = "success",
                                               DT::dataTableOutput("recordtab"),
                                               fluidRow(
                                                   column(1, actionBttn("btn_remove", "Remove Record",color = "success",style = "jelly")),
                                                      column(1, uiOutput("ui_dlbtn"))
                                                   #   column(10, uiOutput("ui_unvisible_columns"))
                                               )
                                           ),
                                           uiOutput("predControl1"),

                                           uiOutput("predControl2")


                                       )
                               ),
                               tabItem(tabName = "menu_tab",
                                       fluidRow(
                                           box(title = "Retrive characterised gate resources", width = 6,solidHeader = T, status = "primary",
                                               selectInput(
                                                   "gatetype",
                                                   "Select gate type from",
                                                   choices = c(
                                                       "Sensor" = 1,
                                                       "NOT gate" =2,
                                                       "AND gate" =3
                                                       # "NAND gate"=4
                                                   ), selected=NULL

                                               ),
                                               box(title = "Characterised gate resource", width = 12, solidHeader = T, status = "success",
                                                   DT::dataTableOutput("gate_resources")
                                               )
                                               ),
                                           tabBox(title = "Compare selected gates", width = 6,
                                                  id = "tabset_menu_tab",
                                                  tabPanel("Radar plot",
                                                           plotlyOutput(outputId="gate_compare_radar",
                                                                        width = "100%",
                                                                        height = "600px")
                                                  ),
                                                  tabPanel("Output-Input relation",
                                                           plotlyOutput(outputId="gate_compare_io",
                                                                        width = "100%",
                                                                        height = "600px")
                                                           #  plotlyOutput("pmf_plot") %>% withSpinner(type = 5)
                                                  )
                                           )
                                       )
                               ),
                               tabItem(tabName = "menu_about",
                                           includeMarkdown("www/helppage.Rmd")
                               )
                           )
                       )
)

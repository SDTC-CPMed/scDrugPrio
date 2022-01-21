library(shiny)
library(shinyjs)
ui <- fluidPage(
  actionButton("reset", "Highlight all"),
  plotlyOutput("network"),
)
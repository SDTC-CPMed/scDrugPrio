library(shiny)
library(shinyjs)
library(igraph)
library(plotly)
source("ui.R")

#Load all data created in 'network.R'
load("shiny.RData")

#Factorize the nodes
point_data$col_DEGs = as.factor(point_data$col_DEGs)
levels(point_data$col_DEGs) = c('DEGs', 'Unknown drugs', 'Known drugs')

server <- function(input, output) {
  
  idx = 0
  values <- reactiveValues(col_DEGs_ord = col_DEGs[order(col_DEGs)],
                           size = node_size_DEGs,
                           idx = 0,
                           connections = 0,
                           cam = cam)
  
  #If the highlight button is clicked
  observeEvent(input$reset,{
    values$idx = 0
    values$cam = event_data('plotly_relayout')
    })
    
    
  #Set up the graph 
  observeEvent(values$idx,{ #If idx is changes, create a new output
    output$network <- renderPlotly({
      
      #To turn up the whole network
      if(values$idx == 0){
        fig = plot_ly()
        fig = fig %>% add_trace(
          x=line_data$Xe,
          y=line_data$Ye,
          z=line_data$Ze,
          mode='lines',
          line=list(color='grey',
                    width=0.3),
          hoverinfo = 'none',
          opacity = 0.4,
          name = 'PP connections'
          
        )
        
        
        axis=list(
          showticklabels=FALSE,
          showspikes = FALSE,
          title=''
        )
        

        
        fig = fig %>% layout(
          title="Interaction between drug candidates and AIA plasma cell DEGs",
          width=1000,
          height=1000,
          showlegend=TRUE,
          scene=list(
            zaxis = list(range = c(-10, 10),
                         showticklabels=FALSE,
                         showspikes = FALSE,
                         title=''),
            xaxis = list(showticklabels=FALSE,
                         showspikes = FALSE,
                         title=''),
            yaxis = list(showticklabels=FALSE,
                         showspikes = FALSE,
                         title=''),
            camera = values$cam$scene.camera
            

          ),
          margin=list(
            t=100
          ),
          hovermode='closest'
        )

        
        
        
        
        
        fig = fig %>% add_trace( data = point_data,
                                 x=~Xn,
                                 y=~Yn,
                                 z=~Zn,
                                 type = 'scatter3d',
                                 mode='markers',
                                 marker=list(symbol='circle',
                                             size=values$size[order(col_DEGs)],
                                            
                                             color=values$col_DEGs_ord,
                                             colors= c('blue', 'red', 'yellow'),

                                             line=list(color='black', width=0.001)
                                             
                                 ),
                                 color= ~col_DEGs,
                                 text=hover_name_DEGs,
                                 hoverinfo='text',
                                 opacity = 1
        )
      }
      

    
      # In case one of the nodes was clicked
      else{
        fig = plot_ly()
        
        
        fig = fig %>% add_trace(
          x=line_data$Xe[values$connections],
          y=line_data$Ye[values$connections],
          z=line_data$Ze[values$connections],
          name = 'PP connections',
          mode='lines',
          line=list(color='grey',
                    width=1),
          hoverinfo = 'none',
          opacity = 1,
          showlegend = TRUE
          
        )
        
        fig = fig %>% add_trace(
          x=line_data$Xe[-values$connections],
          y=line_data$Ye[-values$connections],
          z=line_data$Ze[-values$connections],
          name = 'PP connections',
          mode='lines',
          line=list(color='grey',
                    width=0.3),
          hoverinfo = 'none',
          opacity = 0.1,
          showlegend = FALSE
          
        )

        fig = fig %>% add_trace(data = point_data,
                                x=~Xn[-values$idx],
                                y=~Yn[-values$idx],
                                z=~Zn[-values$idx],
                                type = 'scatter3d',
                                mode='markers',
                                marker=list(symbol='circle',
                                            size=values$size[-values$idx],
                                            color=col_DEGs[-values$idx][order(c(col_DEGs[-values$idx]))],
                                            colors= c('blue', 'red', 'yellow'),
                                            
                                            line=list(color='black', width=0.001)
                                            
                                ),
                                color= ~col_DEGs[-values$idx],
                                hoverinfo = 'none',
                                opacity = 0.1,
                                showlegend = FALSE
        )
        
      }
      

      fig
    })
  })
  
  #Light up the neighboors after the click
  observeEvent(event_data("plotly_click"),{
    d = event_data("plotly_click")
    print(d)
    if ((d$curveNumber > 0 & values$idx == 0) | (d$curveNumber > 3 & values$idx > 0)){
      values$cam = event_data('plotly_relayout')

      values$idx = c(unlist(Neighboors[[which(point_data$Xn == d$x)]]))
      values$size = node_size_DEGs[order(col_DEGs)]
      con = c(unlist(Connections[[which(point_data$Xn == d$x)]]))
      values$connections = rep(0, 3*length(con))
      for(i in 1:length(con)){
        values$connections[3*i-2] = con[i]*3-2
        values$connections[3*i-1] = con[i]*3-1
        values$connections[3*i] = con[i]*3
      }
      
    }

    
  })
  
  
  
  
}

shinyApp(ui, server)

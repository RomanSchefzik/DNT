library(shiny)


ui1 <- basicPage(
  title = "DNT",

  navbarPage(
    "Differential network test",

    #1. Panel: Comparison -----------------------
    tabPanel(
      "Comparison",
      fluidRow(
        column(3,
               sliderInput("vSize", "Vertice Size", min=8, max=50, value=16),
               sliderInput("tSize", "Text Size", min=0.4, max=2.5, value=0.8),
        ),
        column(3,
               textInput("negcol", "Color of negative edges", value = "blue"),
               textInput("poscol", "Color of positive edges", value = "red"),
               numericInput("multiplier", "Multiplier with edge weight", value=4)
        ),
        column(3,
               selectInput("layout", "Layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve", "Edge curved", choices=list(TRUE, FALSE))
        ),
        column(3,
               selectInput("method", "Method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method", "Method 2", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic", "Number for EBICglasso",min = 0,max = 1, step = 0.01, value = 0)#,Distcorr = list("DistCorr"))),
        )
      ),
      mainPanel(width = 12,splitLayout(plotOutput("plot1"),
                            plotOutput("plot2")))
    ),
    tabPanel(
      "Survivors Admission",
      sliderInput("t", "T", min = 0, max = 1, value = 0),
      mainPanel(width = 12, plotOutput("plot3"))
    ),
    tabPanel(
      "Survivors Event",
      sliderInput("t", "T", min = 0, max = 1, value = 0),
      mainPanel(width = 12, plotOutput("plot4"))
    ),
    tabPanel(
      "Without Clusters",
      fluidRow(
        column(3,
               sliderInput("vSize5", "Vertice Size", min=8, max=50, value=16),
               sliderInput("tSize5", "Text Size", min=0.4, max=2.5, value=0.8),
        ),
        column(3,
               textInput("negcol5", "Color of negative edges", value = "blue"),
               textInput("poscol5", "Color of positive edges", value = "red"),
               numericInput("multiplier5", "Multiplier with edge weight", value=4)
        ),
        column(3,
               selectInput("layout5", "Layout", choices=c(d = list(layout.auto), e = list(layout.fruchterman.reingold))),
               selectInput("curve5", "Edge curved", choices=list(TRUE, FALSE))
        ),
        column(3,
               selectInput("method5", "Method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method5", "Method 2", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic5", "Number for EBICglasso",min = 0, step = 0.3, value = 0)#,Distcorr = list("DistCorr"))),
        )
      ),
      mainPanel(width = 12,splitLayout(plotOutput("plot5"),
                                       plotOutput("plot6")))
    )
  )
)

server1 <- function(input,output, session){
  observe({
    x <- input$method

    if(x == "EBICglasso"){
      updateSelectInput(session, "adj.method",
                        label = paste("Select input label"),
                        choices = list("pearson", "kendall", "spearman"),
                        selected = "pearson")
    }else{
      updateSelectInput(session, "adj.method",
                        label = paste("Select input label"),
                        choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                        selected = "holm")
    }
  })
  observe({
    x <- input$method5

    if(x == "EBICglasso"){
      updateSelectInput(session, "adj.method5",
                        label = paste("Select input label"),
                        choices = list("pearson", "kendall", "spearman"),
                        selected = "pearson")
    }else{
      updateSelectInput(session, "adj.method5",
                        label = paste("Select input label"),
                        choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                        selected = "holm")
    }
  })


  output$plot1 <- renderPlot({
    g1 <- create.Igraphclustering(x1, methodlist = list(input$method, input$adj.method, input$tun.ebic))

    E(g1[[2]])[which(E(g1[[2]])$weight<0)]$color <- input$negcol
    E(g1[[2]])[which(E(g1[[2]])$weight>0)]$color <- input$poscol

    E(g1[[2]])$weight <- abs(E(g1[[2]])$weight)

    if(input$layout == "layout.auto"){
      layout1 = layout.auto
    }
    else if(input$layout == "layout.circle"){
      layout1 = layout.circle
    }
    else if(input$layout == "layout.davidson.harel"){
      layout1 = layout.davidson.harel
    }
    else if(input$layout == "layout.drl"){
      layout1 = layout.drl
    }
    else if(input$layout == "layout.fruchterman.reingold"){
      layout1 = layout.fruchterman.reingold
    }
    else if(input$layout == "layout.fruchterman.reingold.grid"){
      layout1 = layout.fruchterman.reingold.grid
    }
    else if(input$layout == "layout.gem"){
      layout1 = layout.gem
    }
    else if(input$layout == "layout.graphopt"){
      layout1 = layout.graphopt
    }
    else if(input$layout == "layout.grid"){
      layout1 = layout.grid
    }
    else if(input$layout == "layout.kamada.kawai"){
      layout1 = layout.kamada.kawai
    }
    else if(input$layout == "layout.lgl"){
      layout1 = layout.lgl
    }
    else if(input$layout == "layout.mds"){
      layout1 = layout.mds
    }
    else if(input$layout == "layout.reingold.tilford"){
      layout1 = layout.reingold.tilford
    }
    else if(input$layout == "layout.spring"){
      layout1 = layout.spring
    }
    else if(input$layout == "layout.star"){
      layout1 = layout.star
    }
    else{
      layout1 = layout.svd
    }

    plot(g1[[1]], g1[[2]], main="Survivors Admission",edge.curved=input$curve,edge.width=E(g1[[2]])$weight*input$multiplier, edge.color = E(g1[[2]])$color, layout = layout1, vertex.label.cex=input$tSize, vertex.size=input$vSize)
  })

  output$plot2 <- renderPlot({
    g2 <- create.Igraphclustering(x2, methodlist = list(input$method, input$adj.method, input$tun.ebic))

    E(g2[[2]])[which(E(g2[[2]])$weight<0)]$color <- input$negcol
    E(g2[[2]])[which(E(g2[[2]])$weight>0)]$color <- input$poscol

    E(g2[[2]])$weight <- abs(E(g2[[2]])$weight)

    if(input$layout == "layout.auto"){
      layout1 = layout.auto
    }
    else if(input$layout == "layout.circle"){
      layout1 = layout.circle
    }
    else if(input$layout == "layout.davidson.harel"){
      layout1 = layout.davidson.harel
    }
    else if(input$layout == "layout.drl"){
      layout1 = layout.drl
    }
    else if(input$layout == "layout.fruchterman.reingold"){
      layout1 = layout.fruchterman.reingold
    }
    else if(input$layout == "layout.fruchterman.reingold.grid"){
      layout1 = layout.fruchterman.reingold.grid
    }
    else if(input$layout == "layout.gem"){
      layout1 = layout.gem
    }
    else if(input$layout == "layout.graphopt"){
      layout1 = layout.graphopt
    }
    else if(input$layout == "layout.grid"){
      layout1 = layout.grid
    }
    else if(input$layout == "layout.kamada.kawai"){
      layout1 = layout.kamada.kawai
    }
    else if(input$layout == "layout.lgl"){
      layout1 = layout.lgl
    }
    else if(input$layout == "layout.mds"){
      layout1 = layout.mds
    }
    else if(input$layout == "layout.reingold.tilford"){
      layout1 = layout.reingold.tilford
    }
    else if(input$layout == "layout.spring"){
      layout1 = layout.spring
    }
    else if(input$layout == "layout.star"){
      layout1 = layout.star
    }
    else{
      layout1 = layout.svd
    }

    plot(g2[[1]], g2[[2]], main="Survivors Event",edge.curved=input$curve,edge.width=E(g2[[2]])$weight*input$multiplier, edge.color = E(g2[[2]])$color, layout = layout1, vertex.label.cex=input$tSize, vertex.size=input$vSize)
  })


  output$plot3 <- renderPlot({
    df <- create.Igraphclustering(x1, list("Spearman"))
  })

  output$plot5 <- renderPlot({
    #df <- create.Igraphclustering(x1, list("Spearman"))
    g1 <- create.Igraphclustering(x1, methodlist = list(input$method5, input$adj.method5, input$tun.ebic5))

    E(g1[[2]])[which(E(g1[[2]])$weight<0)]$color <- input$negcol5
    E(g1[[2]])[which(E(g1[[2]])$weight>0)]$color <- input$poscol5

    E(g1[[2]])$weight <- abs(E(g1[[2]])$weight)

    if(input$layout5 == "layout.auto"){
      layout1 = layout.auto
    }
    else{
      layout1 = layout.fruchterman.reingold
    }

    plot(g1[[2]], main="Survivors Admission",edge.curved=input$curve5,edge.width=E(g1[[2]])$weight*input$multiplier5, edge.color = E(g1[[2]])$color, layout = layout1, vertex.label.cex=input$tSize5, vertex.size=input$vSize5)
  })

  output$plot6 <- renderPlot({
    g2 <- create.Igraphclustering(x2, methodlist = list(input$method5, input$adj.method5, input$tun.ebic5))

    E(g2[[2]])[which(E(g2[[2]])$weight<0)]$color <- input$negcol5
    E(g2[[2]])[which(E(g2[[2]])$weight>0)]$color <- input$poscol5

    E(g2[[2]])$weight <- abs(E(g2[[2]])$weight)

    if(input$layout5 == "layout.auto"){
      layout1 = layout.auto
    }
    else{
      layout1 = layout.fruchterman.reingold
    }

    plot(g2[[2]], main="Survivors Event",edge.curved=input$curve5,edge.width=E(g2[[2]])$weight*input$multiplier5, edge.color = E(g2[[2]])$color, layout = layout1, vertex.label.cex=input$tSize5, vertex.size=input$vSize5)
  })
}
#shinyApp(ui1, server1)


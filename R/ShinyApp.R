ui <- basicPage(
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
      mainPanel(width = 12,plotOutput("plot1"))
    ),
    #2. Survivors Admission -----------------------
    tabPanel(
      "Survivors Admission",
      fluidRow(
        column(3,
               sliderInput("vSize2", "Vertice Size", min=8, max=50, value=16),
               sliderInput("tSize2", "Text Size", min=0.4, max=2.5, value=0.8),
        ),
        column(3,
               textInput("negcol2", "Color of negative edges", value = "blue"),
               textInput("poscol2", "Color of positive edges", value = "red"),
               numericInput("multiplier2", "Multiplier with edge weight", value=4)
        ),
        column(3,
               selectInput("layout2", "Layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve2", "Edge curved", choices=list(TRUE, FALSE))
        ),
        column(3,
               selectInput("method2", "Method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method2", "Method 2", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic2", "Number for EBICglasso",min = 0,max = 1, step = 0.01, value = 0)#,Distcorr = list("DistCorr"))),
        )
      ),
      mainPanel(width = 12, plotOutput("plot2"))
    ),
    #3. Panel: Survivors Event -----------------------
    tabPanel(
      "Survivors Event",
      fluidRow(
        column(3,
               sliderInput("vSize3", "Vertice Size", min=8, max=50, value=16),
               sliderInput("tSize3", "Text Size", min=0.4, max=2.5, value=0.8),
        ),
        column(3,
               textInput("negcol3", "Color of negative edges", value = "blue"),
               textInput("poscol3", "Color of positive edges", value = "red"),
               numericInput("multiplier3", "Multiplier with edge weight", value=4)
        ),
        column(3,
               selectInput("layout3", "Layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve3", "Edge curved", choices=list(TRUE, FALSE))
        ),
        column(3,
               selectInput("method3", "Method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method3", "Method 2", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic3", "Number for EBICglasso",min = 0,max = 1, step = 0.01, value = 0)#,Distcorr = list("DistCorr"))),
        )
      ),
      mainPanel(width = 12, plotOutput("plot3"))
    ),
    #4. Panel: Without Clusters -----------------------
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
               selectInput("layout5", "Layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve5", "Edge curved", choices=list(TRUE, FALSE))
        ),
        column(3,
               selectInput("method5", "Method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method5", "Method 2", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic5", "Number for EBICglasso",min = 0, step = 0.3, value = 0)#,Distcorr = list("DistCorr"))),
        )
      ),
      mainPanel(width = 12,plotOutput("plot5"))
    )
  )
)

#' Server
#' @description
#' Server for the ShinyApp.
#' @export
#'
server <- function(input,output, session){
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
    x <- input$method2

    if(x == "EBICglasso"){
      updateSelectInput(session, "adj.method2",
                        label = paste("Select input label"),
                        choices = list("pearson", "kendall", "spearman"),
                        selected = "pearson")
    }else{
      updateSelectInput(session, "adj.method2",
                        label = paste("Select input label"),
                        choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                        selected = "holm")
    }
  })
  observe({
    x <- input$method3

    if(x == "EBICglasso"){
      updateSelectInput(session, "adj.method3",
                        label = paste("Select input label"),
                        choices = list("pearson", "kendall", "spearman"),
                        selected = "pearson")
    }else{
      updateSelectInput(session, "adj.method3",
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
    graph.plot(A = x1, B = x2, methodlist = list(input$method, input$adj.method, input$tun.ebic), thresh = 0.05, cluster = TRUE, negcol = input$negcol, poscol = input$poscol, multiplier = input$multiplier, curved = input$curve, layout = layout.func(input$layout), vSize = input$vSize, tSize = input$tSize)
  })

  output$plot2 <- renderPlot({
    graph.plot(A = x1, B = x2, methodlist = list(input$method2, input$adj.method2, input$tun.ebic2), admission = TRUE, event = FALSE, thresh = 0.05, cluster = TRUE, negcol = input$negcol2, poscol = input$poscol2, multiplier = input$multiplier2, curved = input$curve2, layout = layout.func(input$layout2), vSize = input$vSize2, tSize = input$tSize2)
  })

  output$plot3 <- renderPlot({
    graph.plot(A = x1, B = x2, methodlist = list(input$method3, input$adj.method3, input$tun.ebic3), admission = FALSE, event = TRUE, thresh = 0.05, cluster = TRUE, negcol = input$negcol3, poscol = input$poscol3, multiplier = input$multiplier3, curved = input$curve3, layout = layout.func(input$layout3), vSize = input$vSize3, tSize = input$tSize3)
  })

  output$plot5 <- renderPlot({
    graph.plot(A = x1, B = x2, methodlist = list(input$method5, input$adj.method5, input$tun.ebic5), thresh = 0.05, cluster = FALSE, negcol = input$negcol5, poscol = input$poscol5, multiplier = input$multiplier5,  curved = input$curve5, layout = layout.func(input$layout5), vSize = input$vSize5, tSize = input$tSize5)
  })
}

#' Layout function
#' @description
#' Layout function turns the layout-String into the suitable function.
#' @param layout The layout as a String.
#' @details
#' If the String is not a existing layout the layout function is automatically layout.auto.
#' @return The suitable layout function.
#' @examples
#' layout.func("layout.fruchterman.reingold")
#' @export
layout.func <- function(layout){
  if(layout == "layout.circle"){
    layout1 = layout.circle
  }
  else if(layout == "layout.davidson.harel"){
    layout1 = layout.davidson.harel
  }
  else if(layout == "layout.drl"){
    layout1 = layout.drl
  }
  else if(layout == "layout.fruchterman.reingold"){
    layout1 = layout.fruchterman.reingold
  }
  else if(layout == "layout.fruchterman.reingold.grid"){
    layout1 = layout.fruchterman.reingold.grid
  }
  else if(layout == "layout.gem"){
    layout1 = layout.gem
  }
  else if(layout == "layout.graphopt"){
    layout1 = layout.graphopt
  }
  else if(layout == "layout.grid"){
    layout1 = layout.grid
  }
  else if(layout == "layout.kamada.kawai"){
    layout1 = layout.kamada.kawai
  }
  else if(layout == "layout.lgl"){
    layout1 = layout.lgl
  }
  else if(layout == "layout.mds"){
    layout1 = layout.mds
  }
  else if(layout == "layout.reingold.tilford"){
    layout1 = layout.reingold.tilford
  }
  else if(layout == "layout.spring"){
    layout1 = layout.spring
  }
  else if(layout == "layout.star"){
    layout1 = layout.star
  }
  else if(layout == "layout.svd"){
    layout1 = layout.svd
  }
  else{
    layout1 = layout.auto
  }
  return(layout1)
}

#' Create App
#' @description
#' Creates the Shiny App of the Differential Network Tests
#' @export
create.app <- function(){
  shinyApp(ui, server)
}

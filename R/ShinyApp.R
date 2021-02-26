ui <- basicPage(
  title = "DNT",

 navbarPage(
    "Differential network tests",

    #1. Panel: Comparison -----------------------
    tabPanel(
      "Comparison",
      shinybusy::add_busy_spinner(spin = "fading-circle"),
      fluidRow(
        column(3,
               fileInput("fileA", "select .RData, .txt, .csv or .xls file for network A", accept = c(".RData", ".txt", ".csv", ".xls")),
               fileInput("fileB", "select .RData, .txt, .csv or .xls File for network B", accept = c(".RData", ".txt", ".csv", ".xls"))
        ),
        column(3,
               textInput("negcol", "color of edges with negative association", value = "blue"),
               textInput("poscol", "color of edges with positive association", value = "red"),
               numericInput("multiplier", "multiplier for edge width", value=4)
        ),
        column(3,
               selectInput("layout", "layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve", "curved edges", choices=list(TRUE, FALSE)),
               sliderInput("vSize", "node size", min=8, max=50, value=16),
               sliderInput("tSize", "text size", min=0.4, max=2.5, value=0.8)
        ),
        column(3,
               selectInput("method", "network estimation method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method", "p-value adjustment method in network estimation", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic", "tuning parameter for EBICglasso",min = 0,max = 1, step = 0.01, value = 0)#,Distcorr = list("DistCorr"))),
        )
      ),
      mainPanel(width = 12,plotOutput("plot1"))
    ),
    #2. Network A -----------------------
    tabPanel(
      "Network A",
      shinybusy::add_busy_spinner(spin = "fading-circle"),
      fluidRow(
        column(3,
               p("Takes file for network A from 'Comparison' tab if no file is selected here."),
               fileInput("fileA2", "select .RData, .txt, .csv or .xls file", accept = c(".RData", ".txt", ".csv", ".xls"))
        ),
        column(3,
               textInput("negcol2", "color of edges with negative association", value = "blue"),
               textInput("poscol2", "color of edges with positive association", value = "red"),
               numericInput("multiplier2", "multiplier for edge width", value=4)
        ),
        column(3,
               selectInput("layout2", "layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve2", "curved edges", choices=list(TRUE, FALSE)),
               sliderInput("vSize2", "node size", min=8, max=50, value=16),
               sliderInput("tSize2", "text size", min=0.4, max=2.5, value=0.8)
        ),
        column(3,
               selectInput("method2", "network estimation method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method2", "p-value adjustment method in network estimation", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic2", "tuning parameter for EBICglasso",min = 0,max = 1, step = 0.01, value = 0)#,Distcorr = list("DistCorr"))),
        )
      ),
      mainPanel(width = 12, plotOutput("plot2"))
    ),
    #3. Network B -----------------------
    tabPanel(
      "Network B",
      shinybusy::add_busy_spinner(spin = "fading-circle"),
      fluidRow(
        column(3,
               p("Takes file for network B from 'Comparison' tab if no file is selected here."),
               fileInput("fileB3", "select .RData, .txt, .csv or .xls file", accept = c(".RData", ".txt", ".csv", ".xls"))
        ),
        column(3,
               textInput("negcol3", "color of edges with negative association", value = "blue"),
               textInput("poscol3", "color of edges with positive association", value = "red"),
               numericInput("multiplier3", "multiplier for edge width", value=4)
        ),
        column(3,
               selectInput("layout3", "layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve3", "curved edges", choices=list(TRUE, FALSE)),
               sliderInput("vSize3", "node size", min=8, max=50, value=16),
               sliderInput("tSize3", "text size", min=0.4, max=2.5, value=0.8)
        ),
        column(3,
               selectInput("method3", "Method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method3", "p-value adjustment method in network estimation", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic3", "tuning parameter for EBICglasso",min = 0,max = 1, step = 0.01, value = 0)#,Distcorr = list("DistCorr"))),
        )
      ),
      mainPanel(width = 12, plotOutput("plot3"))
    ),
    #4. Panel: Without Clusters -----------------------
    tabPanel(
      "Without Clusters",
      shinybusy::add_busy_spinner(spin = "fading-circle"),
      fluidRow(
        column(3,
               p("Takes files for networks A and B from 'Comparison' tab if no files are selected here."),
               fileInput("fileA5", "select .RData, .txt, .csv or .xls file for network A", accept = c(".RData", ".txt", ".csv", ".xls")),
               fileInput("fileB5", "select .RData, .txt, .csv or .xls file for network B", accept = c(".RData", ".txt", ".csv", ".xls"))
        ),
        column(3,
               textInput("negcol5", "color of edges with negative association", value = "blue"),
               textInput("poscol5", "color of edges with positive association", value = "red"),
               numericInput("multiplier5", "multiplier for edge width", value=4)
        ),
        column(3,
               selectInput("layout5", "layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve5", "curved edges", choices=list(TRUE, FALSE)),
               sliderInput("vSize5", "node size", min=8, max=50, value=16),
               sliderInput("tSize5", "text size", min=0.4, max=2.5, value=0.8)
        ),
        column(3,
               selectInput("method5", "network estimation method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method5", "p-value adjustment method in network estimation", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic5", "tuning parameter for EBICglasso",min = 0, step = 0.3, value = 0)#,Distcorr = list("DistCorr"))),
        )
      ),
      mainPanel(width = 12,plotOutput("plot5"))
    )
  )
)

#' Server function for Shiny app
#' @description
#' Server function for Shiny app
#' @param input list-like object containing all the input data sent from the browser, named according to the input ID
#' @param output list-like object, named according to the output ID
#' @param session environment that can be employed to access information and functionality with respect to the session
#' @details As one never calls the server function by oneself, one will never create the three argument objects by oneself. Instead, they are created by Shiny when the session starts, connecting back to a specific session.
#' @return Server function for Shiny app
#' @export
#'
server <- function(input,output, session){
  observe({
    x <- input$method

    if(x == "EBICglasso"){
      updateSelectInput(session, "adj.method",
                        choices = list("pearson", "kendall", "spearman"),
                        selected = "pearson")
    }else{
      updateSelectInput(session, "adj.method",
                        choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                        selected = "holm")
    }
  })
  observe({
    x <- input$method2

    if(x == "EBICglasso"){
      updateSelectInput(session, "adj.method2",
                        choices = list("pearson", "kendall", "spearman"),
                        selected = "pearson")
    }else{
      updateSelectInput(session, "adj.method2",
                        choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                        selected = "holm")
    }
  })
  observe({
    x <- input$method3

    if(x == "EBICglasso"){
      updateSelectInput(session, "adj.method3",
                        choices = list("pearson", "kendall", "spearman"),
                        selected = "pearson")
    }else{
      updateSelectInput(session, "adj.method3",
                        choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                        selected = "holm")
    }
  })
  observe({
    x <- input$method5

    if(x == "EBICglasso"){
      updateSelectInput(session, "adj.method5",
                        choices = list("pearson", "kendall", "spearman"),
                        selected = "pearson")
    }else{
      updateSelectInput(session, "adj.method5",
                        choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                        selected = "holm")
    }
  })



  output$plot1 <- renderPlot({

    fileA <- input$fileA
    fileB <- input$fileB

    req(fileA)
    req(fileB)

    extA <- tools::file_ext(fileA$datapath)
    validate(need(extA %in% c("RData", "csv", "xls", "txt"), "please upload a .RData, .xls, .csv or .txt file for network A"))

    extB <- tools::file_ext(fileB$datapath)
    validate(need(extB %in% c("RData", "csv", "xls", "txt"), "please upload a .RData, .xls, .csv or .txt file for network B"))

    shinybusy::show_spinner()

    dataA <- import.table(fileA$datapath)
    dataB <- import.table(fileB$datapath)

    comp.plot(A = dataA, B = dataB, methodlist = list(input$method, input$adj.method, input$tun.ebic), thresh = 0.05, cluster = TRUE, negcol = input$negcol, poscol = input$poscol, multiplier = input$multiplier, curved = input$curve, layout = layout.func(input$layout), vSize = input$vSize, tSize = input$tSize)

    shinybusy::hide_spinner()
  })

  output$plot2 <- renderPlot({

    file <- input$fileA2

    if(is.null(file)){
      file <- input$fileA
    }

    req(file)

    ext <- tools::file_ext(file$datapath)
    validate(need(ext %in% c("RData", "csv", "xls", "txt"), "please upload a .RData, .xls, .csv or .txt file"))

    shinybusy::show_spinner()

    dataA <- import.table(file$datapath)
    dataB <- matrix(0, nrow = length(dataA), ncol = length(dataA))

    comp.plot(A = dataA, B = dataB, methodlist = list(input$method2, input$adj.method2, input$tun.ebic2), networkA = TRUE, networkB = FALSE, thresh = 0.05, cluster = TRUE, negcol = input$negcol2, poscol = input$poscol2, multiplier = input$multiplier2, curved = input$curve2, layout = layout.func(input$layout2), vSize = input$vSize2, tSize = input$tSize2)

    shinybusy::hide_spinner()
  })

  output$plot3 <- renderPlot({

    file <- input$fileB3

    if(is.null(file)){
      file <- input$fileB
    }

    req(file)

    ext <- tools::file_ext(file$datapath)
    validate(need(ext %in% c("RData", "csv", "xls", "txt"), "please upload a .RData, .xls, .csv or .txt file"))

    shinybusy::show_spinner()

    dataB <- import.table(file$datapath)
    dataA <- matrix(0, nrow = length(dataB), ncol = length(dataB))

    comp.plot(A = dataA, B = dataB, methodlist = list(input$method3, input$adj.method3, input$tun.ebic3), networkA = FALSE, networkB = TRUE, thresh = 0.05, cluster = TRUE, negcol = input$negcol3, poscol = input$poscol3, multiplier = input$multiplier3, curved = input$curve3, layout = layout.func(input$layout3), vSize = input$vSize3, tSize = input$tSize3)

    shinybusy::hide_spinner()
  })

  output$plot5 <- renderPlot({

    fileA <- input$fileA5
    fileB <- input$fileB5

    if(is.null(fileA)){
      fileA <- input$fileA
    }
    if(is.null(fileB)){
      fileB <- input$fileB
    }

    req(fileA)
    req(fileB)

    extA <- tools::file_ext(fileA$datapath)
    validate(need(extA %in% c("RData", "csv", "xls", "txt"), "please upload a .RData, .xls, .csv or .txt file for network A"))

    extB <- tools::file_ext(fileB$datapath)
    validate(need(extB %in% c("RData", "csv", "xls", "txt"), "please upload a .RData, .xls, .csv or .txt file for network B"))

    shinybusy::show_spinner()

    dataA <- import.table(fileA$datapath)
    dataB <- import.table(fileB$datapath)

    comp.plot(A = dataA, B = dataB, methodlist = list(input$method5, input$adj.method5, input$tun.ebic5), thresh = 0.05, cluster = FALSE, negcol = input$negcol5, poscol = input$poscol5, multiplier = input$multiplier5,  curved = input$curve5, layout = layout.func(input$layout5), vSize = input$vSize5, tSize = input$tSize5)

    shinybusy::hide_spinner()
  })
}

#' Creates layout function
#' @description
#' This function transforms a layout string from the igraph package into a corresponding layout function.
#' @param layout the layout as a string
#' @details
#' This function transforms a layout string from the igraph package into a corresponding layout function. If the string does not represent an existing layout from the igraph package, the layout function is automatically set to layout.auto.
#' @return the corresponding layout function
#' @export
layout.func <- function(layout){
  if(layout == "layout.circle"){
    layout1 = igraph::layout.circle
  }
  else if(layout == "layout.davidson.harel"){
    layout1 = igraph::layout.davidson.harel
  }
  else if(layout == "layout.drl"){
    layout1 = igraph::layout.drl
  }
  else if(layout == "layout.fruchterman.reingold"){
    layout1 = igraph::layout.fruchterman.reingold
  }
  else if(layout == "layout.fruchterman.reingold.grid"){
    layout1 = igraph::layout.fruchterman.reingold.grid
  }
  else if(layout == "layout.gem"){
    layout1 = igraph::layout.gem
  }
  else if(layout == "layout.graphopt"){
    layout1 = igraph::layout.graphopt
  }
  else if(layout == "layout.grid"){
    layout1 = igraph::layout.grid
  }
  else if(layout == "layout.kamada.kawai"){
    layout1 = igraph::layout.kamada.kawai
  }
  else if(layout == "layout.lgl"){
    layout1 = igraph::layout.lgl
  }
  else if(layout == "layout.mds"){
    layout1 = igraph::layout.mds
  }
  else if(layout == "layout.reingold.tilford"){
    layout1 = igraph::layout.reingold.tilford
  }
  else if(layout == "layout.spring"){
    layout1 = igraph::layout.spring
  }
  else if(layout == "layout.star"){
    layout1 = igraph::layout.star
  }
  else if(layout == "layout.svd"){
    layout1 = igraph::layout.svd
  }
  else{
    layout1 = igraph::layout.auto
  }
  return(layout1)
}


#' Imports an external table into R
#' @description
#' This function imports an external table into R.
#' @param directory the directory to the external table that should be imported into R
#' @details
#' This function imports an external table into R. It causes that only columns with integers or doubles are kept and works with .RData, .txt, .xls and .csv files as input.
#' @return
#' the imported data table
#' @export
#'
import.table <- function(directory){
  ending <- tools::file_ext(directory)
  stopifnot("this data format is not supported" = ending %in% c("xls", "RData", "csv", "txt"))
  table <- 0
  if(ending == "RData"){
    table <- get(load(directory))
  }
  else if(ending == "xls"){
    table <- readxl::read_excel(directory)
  }
  else if(ending == "csv"){
    table <- utils::read.csv2(directory)
  }
  else if(ending == "txt"){
    table <- utils::read.table(directory)
  }
  table <- dplyr::select(table, where(is.double), where(is.integer))
  return(table)
}
utils::globalVariables("where")

#' Creates Shiny app
#' @description
#' Creates the Shiny app corresponding to the DNT package
#' @export
create.app <- function(){
  shinyApp(ui, server)
}


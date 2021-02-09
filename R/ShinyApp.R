ui <- basicPage(
  title = "DNT",

 navbarPage(
    "Differential network test",

    #1. Panel: Comparison -----------------------
    tabPanel(
      "Comparison",
      shinybusy::add_busy_spinner(spin = "fading-circle"),
      fluidRow(
        column(3,
               fileInput("fileA", "Choose RData, TXT, CSV or XLS File for Network A", accept = c(".RData", ".txt", ".csv", ".xls")),
               fileInput("fileB", "Choose RData, TXT, CSV or XLS File for Network B", accept = c(".RData", ".txt", ".csv", ".xls"))
        ),
        column(3,
               textInput("negcol", "Color of negative edges", value = "blue"),
               textInput("poscol", "Color of positive edges", value = "red"),
               numericInput("multiplier", "Multiplier with edge weight", value=4)
        ),
        column(3,
               selectInput("layout", "Layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve", "Edge curved", choices=list(TRUE, FALSE)),
               sliderInput("vSize", "Vertice Size", min=8, max=50, value=16),
               sliderInput("tSize", "Text Size", min=0.4, max=2.5, value=0.8)
        ),
        column(3,
               selectInput("method", "Method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method", "Adjust method", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic", "Number for EBICglasso",min = 0,max = 1, step = 0.01, value = 0)#,Distcorr = list("DistCorr"))),
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
               p("Takes file for Network A from 'Comparison'-Tab if no file is selected here."),
               fileInput("fileA2", "Choose RData, TXT, CSV or XLS File", accept = c(".RData", ".txt", ".csv", ".xls"))
        ),
        column(3,
               textInput("negcol2", "Color of negative edges", value = "blue"),
               textInput("poscol2", "Color of positive edges", value = "red"),
               numericInput("multiplier2", "Multiplier with edge weight", value=4)
        ),
        column(3,
               selectInput("layout2", "Layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve2", "Edge curved", choices=list(TRUE, FALSE)),
               sliderInput("vSize2", "Vertice Size", min=8, max=50, value=16),
               sliderInput("tSize2", "Text Size", min=0.4, max=2.5, value=0.8)
        ),
        column(3,
               selectInput("method2", "Method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method2", "Adjust method", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic2", "Number for EBICglasso",min = 0,max = 1, step = 0.01, value = 0)#,Distcorr = list("DistCorr"))),
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
               p("Takes file for Network B from 'Comparison'-Tab if no file is selected here."),
               fileInput("fileB3", "Choose RData, TXT, CSV or XLS File", accept = c(".RData", ".txt", ".csv", ".xls"))
        ),
        column(3,
               textInput("negcol3", "Color of negative edges", value = "blue"),
               textInput("poscol3", "Color of positive edges", value = "red"),
               numericInput("multiplier3", "Multiplier with edge weight", value=4)
        ),
        column(3,
               selectInput("layout3", "Layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve3", "Edge curved", choices=list(TRUE, FALSE)),
               sliderInput("vSize3", "Vertice Size", min=8, max=50, value=16),
               sliderInput("tSize3", "Text Size", min=0.4, max=2.5, value=0.8)
        ),
        column(3,
               selectInput("method3", "Method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method3", "Adjust method", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
               numericInput("tun.ebic3", "Number for EBICglasso",min = 0,max = 1, step = 0.01, value = 0)#,Distcorr = list("DistCorr"))),
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
               p("Takes file for Network A and B from 'Comparison'-Tab if no file is selected here."),
               fileInput("fileA5", "Choose RData, TXT, CSV or XLS File for Network A", accept = c(".RData", ".txt", ".csv", ".xls")),
               fileInput("fileB5", "Choose RData, TXT, CSV or XLS File for Network B", accept = c(".RData", ".txt", ".csv", ".xls"))
        ),
        column(3,
               textInput("negcol5", "Color of negative edges", value = "blue"),
               textInput("poscol5", "Color of positive edges", value = "red"),
               numericInput("multiplier5", "Multiplier with edge weight", value=4)
        ),
        column(3,
               selectInput("layout5", "Layout", choices=list("layout.auto", "layout.circle", "layout.davidson.harel", "layout.drl", "layout.fruchterman.reingold", "layout.fruchterman.reingold.grid", "layout.gem", "layout.graphopt", "layout.grid", "layout.kamada.kawai", "layout.lgl", "layout.mds", "layout.reingold.tilford", "layout.spring", "layout.star", "layout.svd")),
               selectInput("curve5", "Edge curved", choices=list(TRUE, FALSE)),
               sliderInput("vSize5", "Vertice Size", min=8, max=50, value=16),
               sliderInput("tSize5", "Text Size", min=0.4, max=2.5, value=0.8)
        ),
        column(3,
               selectInput("method5", "Method to create adjacency matrix", choices = list("Spearman","PCSpearman","Spearman.adj","PCSpearman.adj", "DistCorr", "DistCorr.adj", "EBICglasso")),
               selectInput("adj.method5", "Adjust method", choices = list("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
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
    validate(need(extA %in% c("RData", "csv", "xls", "txt"), "Please upload a RData, xls, csv or txt file for Network A"))

    extB <- tools::file_ext(fileB$datapath)
    validate(need(extB %in% c("RData", "csv", "xls", "txt"), "Please upload a RData, xls, csv or txt file for Network B"))

    shinybusy::show_spinner()

    dataA <- import.table(fileA$datapath)
    dataB <- import.table(fileB$datapath)

    graph.plot(A = dataA, B = dataB, methodlist = list(input$method, input$adj.method, input$tun.ebic), thresh = 0.05, cluster = TRUE, negcol = input$negcol, poscol = input$poscol, multiplier = input$multiplier, curved = input$curve, layout = layout.func(input$layout), vSize = input$vSize, tSize = input$tSize)

    shinybusy::hide_spinner()
  })

  output$plot2 <- renderPlot({

    file <- input$fileA2

    if(is.null(file)){
      file <- input$fileA
    }

    req(file)

    ext <- tools::file_ext(file$datapath)
    validate(need(ext %in% c("RData", "csv", "xls", "txt"), "Please upload a RData, xls, csv or txt file"))

    shinybusy::show_spinner()

    dataA <- import.table(file$datapath)
    dataB <- matrix(0, nrow = length(dataA), ncol = length(dataA))

    graph.plot(A = dataA, B = dataB, methodlist = list(input$method2, input$adj.method2, input$tun.ebic2), networkA = TRUE, networkB = FALSE, thresh = 0.05, cluster = TRUE, negcol = input$negcol2, poscol = input$poscol2, multiplier = input$multiplier2, curved = input$curve2, layout = layout.func(input$layout2), vSize = input$vSize2, tSize = input$tSize2)

    shinybusy::hide_spinner()
  })

  output$plot3 <- renderPlot({

    file <- input$fileB3

    if(is.null(file)){
      file <- input$fileB
    }

    req(file)

    ext <- tools::file_ext(file$datapath)
    validate(need(ext %in% c("RData", "csv", "xls", "txt"), "Please upload a RData, xls, csv or txt file"))

    shinybusy::show_spinner()

    dataB <- import.table(file$datapath)
    dataA <- matrix(0, nrow = length(dataB), ncol = length(dataB))

    graph.plot(A = dataA, B = dataB, methodlist = list(input$method3, input$adj.method3, input$tun.ebic3), networkA = FALSE, networkB = TRUE, thresh = 0.05, cluster = TRUE, negcol = input$negcol3, poscol = input$poscol3, multiplier = input$multiplier3, curved = input$curve3, layout = layout.func(input$layout3), vSize = input$vSize3, tSize = input$tSize3)

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
    validate(need(extA %in% c("RData", "csv", "xls", "txt"), "Please upload a RData, xls, csv or txt file for Network A"))

    extB <- tools::file_ext(fileB$datapath)
    validate(need(extB %in% c("RData", "csv", "xls", "txt"), "Please upload a RData, xls, csv or txt file for Network B"))

    shinybusy::show_spinner()

    dataA <- import.table(fileA$datapath)
    dataB <- import.table(fileB$datapath)

    graph.plot(A = dataA, B = dataB, methodlist = list(input$method5, input$adj.method5, input$tun.ebic5), thresh = 0.05, cluster = FALSE, negcol = input$negcol5, poscol = input$poscol5, multiplier = input$multiplier5,  curved = input$curve5, layout = layout.func(input$layout5), vSize = input$vSize5, tSize = input$tSize5)

    shinybusy::hide_spinner()
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


#' Import table
#' @description
#' The function loads a table out of an extern file into R.
#' @param directory The directory to the file that should be loaded into R as a tibble.
#' @details
#' The function will only keep the columns with integers or doubles and the function can work with RData, text, excel, csv files.
#' @return
#' The data loaded out of the data file.
#' @export
#'
import.table <- function(directory){
  ending <- tools::file_ext(directory)
  stopifnot("This data format is not supported" = ending %in% c("xls", "RData", "csv", "txt"))
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

#' Create App
#' @description
#' Creates the Shiny App of the Differential Network Tests
#' @export
create.app <- function(){
  shinyApp(ui, server)
}


# server.R

# the following packages are required before commerce the program
# install.packages(c("shinydashboard", "DT","shiny", "ggplot2", "gplots", "tidyverse", "heatmaply",
#               "RColorBrewer", "plotly", "networkD3", "igraph", "reshape2"))
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("limma")


## Load required packages
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(gplots)
library(limma)
library(heatmaply)
library(RColorBrewer)
library(tidyverse)
library(plotly)
library(networkD3)
library(igraph)
library(reshape2)


## Start shiny server
shinyServer(function(input, output, session) {
  
  nTopGenes <- reactive({
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:input$topGenes]
    highly_variable_lcpm <- logcounts_heatmap[select_var,]
    highly_variable_lcpm <- as.data.frame(highly_variable_lcpm)
    colnames(highly_variable_lcpm) <- targets$sample_name
    return(highly_variable_lcpm)
  })
  
  
  # plot MDS
  output$mdsplot <- renderPlotly({
    
    mp <- plotMDS(nTopGenes(), top=dim(nTopGenes())[1], ndim=3, plot=FALSE)
    dat <- as.data.frame(mp$cmdscale.out)
    rownames(dat) <- targets$sample_ID
    colnames(dat) <- c("dim_1", "dim_2", "dim_3")
    dat1 <- cbind(dat, targets)
    colors <- c('#4AC6B7', '#1972A4')
    mdsplot <- plot_ly(dat1, x = ~dim_1, y = ~dim_2, z = ~dim_3, color=~group, colors=colors, 
            text = ~paste('Name:', sample_name, '<br>Group:', mds_cluster)) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'dim 1'),
                          yaxis = list(title = 'dim 2'),
                          zaxis = list(title = 'dim 3'),
                          autosize=T
                          # width=800,
                          # height=800,
                          ),
             title=paste("Top", dim(nTopGenes())[1], "most variables genes", sep=" "))
    mdsplot
  })
  
  # color for clusters
  mypalette <- brewer.pal(9,"Spectral")
  
  # output dendrogram
  output$dendro <- renderDendroNetwork({
    hc <- hclust(dist(t(nTopGenes())), "complete")
    dendroNetwork(hc, 
                  treeOrientation="vertical", 
                  textColour=mypalette[cutree(hc, input$ncluster)], 
                  fontSize=13)
  })
  
  #output sample information
  output$sampleTable <- DT::renderDataTable({
    targets %>%
      transmute(sampleID = sample_ID,
                sampleName = sample_name,
                disease = group,
                MDS_cluster = group2)
  }, options = list(searching = TRUE, pageLength = 5))
  
  
  ## top gene for heatmap tab
  nTopGeneH <- reactive({
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:input$topGeneH]
    highly_variable_lcpm <- logcounts_heatmap[select_var,]
    highly_variable_lcpm <- as.data.frame(highly_variable_lcpm)
    colnames(highly_variable_lcpm) <- targets$sample_name
    return(highly_variable_lcpm)
  })
  
  # ## output heatmap
  # output$heatmap <- renderD3heatmap({
  #   d3heatmap(x = nTopGeneH(),
  #             scale="row",
  #             cexCol = 0.2)
  # })
  
  # plot heatmap using plotly
  output$heatmap <- renderPlotly({
    heatmaply(nTopGeneH(), 
              scale = "row",# col=rev(morecols(50)), 
              #ColSideColors=targets$group,
              showticklabels = c(TRUE, TRUE),
              fontsize_row = input$heatmapFontsize,
              fontsize_col = input$heatmapSampleFontsize,
              #fontsize_col = 7,
              #labRow = highly_variable_lcpm$Symbol,
              plot_method = 'plotly')
  })
  
  
  # dataset chosen
  GEP <- reactive({
    if(is.null(input$dataset))
      return()
    
    dat <- get(input$dataset)
    dat <- as.data.frame(dat)
    
    selectedFdr <- as.numeric(input$fdr)
    selectedFC <- as.numeric(input$fc)
    
    dat %>%
      mutate(
        type=factor(ifelse((logFC > selectedFC & adj.P.Val < selectedFdr), "sigUp",
                           ifelse((logFC < (selectedFC * -1) & adj.P.Val < selectedFdr), "sigDown", "NotSig"))),
        minusLog10Pvalue = -log10(P.Value),
        tooltip = paste(Symbol, "logFC:", signif(logFC, digits = 3)
                        , "FDR :", signif(adj.P.Val,digits=2), sep = "\n"))
    
  })
  
  ## scatter plot
  output$volcanoPlot <- renderPlotly({

    plot <- GEP() %>%
      ggplot(aes(x = logFC,
                 y = minusLog10Pvalue,
                 colour = type,
                 text = tooltip,
                 key = row.names(GEP()))) +
      geom_point() +
      xlab("log fold change") +
      ylab("-log10(P-value)") +
      scale_fill_manual("grey", "green", "red")

    plot %>%
      ggplotly(tooltip = "tooltip") %>%
      layout(dragmode = "select")


  })
  
  #plot table
  output$selectedGeneTable <- DT::renderDataTable({

    eventData <- event_data("plotly_click")

    selectedData <- GEP() %>% slice(0)
    if (!is.null(eventData)) selectedData <- GEP()[eventData$key,]

    selectedData %>%
      transmute(
        Gene = Symbol,
        `log fold change` = signif(logFC, digits = 3),
        `p-value` = signif(P.Value, digits = 2),
        `FDR p-value` = signif(adj.P.Val, digits = 2)
      )
  },
  options = list(dom = "tip", pageLength = 10, searching = FALSE))
  
  
  ## plot stripchart
  output$genechart <- renderPlotly({
    eventData <- event_data("plotly_click")
    
    selectedData1 <- GEP() %>% slice(0)
    if (!is.null(eventData)) {
      selectedData1 <- GEP()[eventData$key,]
      
      targeted_symbol <- as.character(selectedData1$Symbol)
      plot <- ccat %>%
        filter(Symbol==targeted_symbol) %>%
        ggplot(aes(x = group2,
                   y = value,
                   colour = group2,
                   text = tooltip)) +
        geom_boxplot()+
        geom_jitter(position=position_jitter(0.1), size=2) +
        xlab("group") +
        ylab("gene expression (log2)") +
        ggtitle(paste(targeted_symbol, "gene expression", sep=" ")) +
        scale_color_brewer(palette="Dark2")
      
      plot %>%
        ggplotly(tooltip = "tooltip")
      
    }
    else {
      plotly_empty(type="scatter", mode = 'markers')
    }
  })
  
  
  # dataset selected for Genesearch
  GEP1 <- reactive({
    if(is.null(input$dataset_gs))
      return()
    
    dat1 <- get(input$dataset_gs)
    dat1 <- as.data.frame(dat1)
    
    selectedFdr_gs <- as.numeric(input$fdr_gs)
    selectedFC_gs <- as.numeric(input$fc_gs)
    
    dat1 %>%
      transmute(`Gene ID` = GeneID, 
                chromosome = Chr, 
                `Gene Symbol` = Symbol, 
                `log fold change` = signif(logFC, digits = 3),
                `Average expression` = signif(AveExpr, digits = 3),
                `p-value` = signif(P.Value, digits = 2),
                `FDR p-value` = signif(adj.P.Val, digits = 2)) %>%
      arrange(`FDR p-value`, `p-value`) %>%
      filter(`FDR p-value` < selectedFdr_gs,
             abs(`log fold change`) > selectedFC_gs)
    
  })
  
  # output significant gene table
  output$sigTable <- DT::renderDataTable({
    GEP1()
  }, selection="single")
  
  
  
  ## plot stripchart for Genesearch
  output$genechart_gs <- renderPlotly({
    s = input$sigTable_rows_selected
    if (length(s)) {
      selectedData <- GEP1()
      selectedData1 <- selectedData[as.numeric(s),]
      
      targeted_symbol <- as.character(selectedData1$`Gene Symbol`)
      plot <- ccat %>%
        filter(Symbol==targeted_symbol) %>%
        ggplot(aes(x = group2,
                   y = value,
                   colour = group2,
                   text = tooltip)) +
        geom_boxplot()+
        geom_jitter(position=position_jitter(0.1), size=2) +
        xlab("group") +
        ylab("gene expression (log2)") +
        ggtitle(paste(targeted_symbol, "gene expression", sep=" ")) +
        scale_color_brewer(palette="Dark2")
      
      plot %>%
        ggplotly(tooltip = "tooltip")
      
    }
    else {
      plotly_empty(type="scatter", mode = 'markers')
    }
    
  })
  options(warn = -1)
  network_gep <- reactive({
    if(is.null(input$dataset_network))
      return()
    
    dat_net <- get(input$dataset_network)
    dat_net <- as.data.frame(dat_net)
    
    selectedFdr_net <- as.numeric(input$fdr_network)
    selectedFC_net <- as.numeric(input$fc_network)
    
    res <- dat_net %>%
      select(Symbol, logFC, adj.P.Val) %>%
      filter(abs(logFC) > selectedFC_net,
             adj.P.Val < selectedFdr_net)
    
    res <- res[1:input$topnetworkgenes,]
    res
  })
  
  
  #plot network analysis
  output$force <- renderForceNetwork({
    res <- network_gep()
    res_sym <- data.frame(Symbol=toupper(res$Symbol),stringsAsFactors = F)
    datA <- inner_join(hi2,res_sym,by=c("Symbol.A"="Symbol"))
    datB <- inner_join(hi2,res_sym,by=c("Symbol.B"="Symbol"))
    datx <- rbind(datA,datB)
    
    # remove duplicates
    datx2 <- unique(datx[c("Symbol.A","Symbol.B")])
    
    vertices <- data.frame("name"=unique(unlist(datx2)))
    net <- graph_from_data_frame(datx2,directed=F,vertices=vertices)
    #
    ## http://www.r-graph-gallery.com/87-interactive-network-with-networkd3-package/
    vertices$group <- edge.betweenness.community(net)$membership  # betweeness centrality for each node for grouping
    
    
    # remove cluster with few nodes
    groupSummary <- table(vertices$group)
    KeepIdx <- which(groupSummary>= input$minCluster)
    KeepVertices <- vertices[vertices$group %in% KeepIdx,1]
    datx3 <- datx2[datx2$Symbol.A %in% KeepVertices & datx2$Symbol.B %in% KeepVertices,]
    vertices <- data.frame("name"=unique(unlist(datx3)))
    net <- graph_from_data_frame(datx3,directed=F,vertices=vertices)
    vertices$group <- edge.betweenness.community(net)$membership  # betweeness centrality for each node for grouping
    
    # choose the node with the largest betweeness value from each group
    verticesBetweeness <- data.frame(Vertices=V(net)$name,bw=betweenness(net,directed=F))
    dfGroup <- inner_join(vertices,verticesBetweeness,by=c("name"="Vertices")) %>%
      group_by(group) %>%
      filter(rank(-bw,ties.method="first")==1) %>% # minus since rank sorts from min to max
      arrange(group) %>%
      mutate(groupName=paste0(group,' - ',name,' (',round(bw),')')) %>%
      select(group,groupName)
    
    vertices <- left_join(vertices,dfGroup,by=c('group'))
    
    relations <- datx3
    relations$source.index = match(relations$Symbol.A, vertices$name)-1
    relations$target.index = match(relations$Symbol.B, vertices$name)-1
    
    forceNetwork(Links = relations, Nodes = vertices,
                 Source = "source.index", Target = "target.index",
                 NodeID = "name",
                 Group = "groupName", # color nodes by betweeness calculated earlier
                 charge = -70, # node repulsion
                 linkDistance = 25,
                 fontSize = 16,
                 legend=T,
                 zoom = T, 
                 opacity = input$opacity)
  })
  
})

## End shiny server

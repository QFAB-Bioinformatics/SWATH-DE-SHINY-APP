require(shiny)
require(ggplot2)
require(preprocessCore)
require(ggrepel)
require(svglite)
## function performing the differential expression analysis
DEprot <- function(data,nbCond=NULL, rep,condName=NULL, normalization="MEAN"){
  
  ##applying different normalization depending on the normalization parameter 
  normalization=toupper(normalization)
  if(identical(normalization, "NULL")){
  }
  else if(identical(normalization, "MEAN")){
    data<-scale(log2(data))
  }
  else if(identical(normalization, "QUANTILE")){
    require(preprocessCore)
    data <- log2(data)
    data <- normalize.quantiles(as.matrix(data))
  } 
  else if(identical(normalization, "MEDIAN")){
    data<-log2(data)
    for (x in 1:length(data)){
      data[,x]<-(data[,x]-median(data[,x],na.rm = T))/(mad(data[,x],na.rm = T))
    }
  }
  else{
    data<-scale(log2(data))
  }
  
  col<-1
  ##cond is a list containing the columns number of each different conditions.
  cond<-list()
  ## the p value is stored in a 3d array, the two first dimmention beeing the comparated conditions and the third dimension the values for each protein.
  p.value<-array(dim = c(nbCond,nbCond,nrow(data))) 
  ## the adjusted p value is stored in a 3d array, the two first dimmention beeing the comparated conditions and the third dimension the values for each protein.
  adjust.p.value<-array(dim = c(nbCond,nbCond,nrow(data)))
  ## the fold change is stored in a 3d array, the two first dimmention beeing the comparated conditions and the third dimension the values for each protein.
  fc<-array(dim = c(nbCond,nbCond,nrow(data)))
  
  ## resultsTable is a 2D matrix of data frames containing the differents values calculated for each comparison. 
  resultsTable<-matrix(data.frame(), nrow = nbCond, ncol = nbCond)
  
  ## if the parameter condName is null, or there is an error in the number of name in the vector, the name are set numbers.
  if(is.null(condName) | length(condName) != nbCond){
    condName<-c(1:nbCond)
  }
  
  ## Setting the cond list depending on the number of conditions and replicates in parameters
  for(i in 1:nbCond){
    cond[[i]]<-c(as.numeric(rep[[i]]))
  }
  data<-as.data.frame(data[complete.cases(data),])
  
  ## Calcul of the p value and fold change for each condition, and setting the resultsTable for each comparisons
  for(i in 1:nbCond){
    for(j in 1:nbCond){
      ##student's test
      p.value[i,j,]<-apply(data,1,function(x){t.test(as.numeric(x[c(cond[[i]])]),as.numeric(x[c(cond[[j]])]), alternative = "t") $p.value})
      ## Benjamini hochberg correction of the p value
      adjust.p.value[i,j,]<-p.adjust(p.value[i,j,], method = "BH")
      
      ##fold change
      fc[i,j,]<-rowMeans(data[,cond[[i]]])-rowMeans(data[,cond[[j]]])
      
      resultsTable[[i,j]]<-data.frame(protein=c(row.names(data)))
      resultsTable[[i,j]][paste("p.value.",condName[i],".vs.",condName[j],sep="")]=c(p.value[i,j,])
      resultsTable[[i,j]][paste("adjust.p.value.",condName[i],".vs.",condName[j],sep="")]=c(adjust.p.value[i,j,])
      resultsTable[[i,j]][paste("fc.",condName[i],".vs.",condName[j],sep="")]=c(fc[i,j,])
    }
  }
  ##returns a list containing a data frame with the normalized values and a matrix containing the results of the differnt comparisons
  results <-list( "normdata" = data,"dataprocessed"=resultsTable)
  return(results)
}

norm <- function(data, normalization="MEAN"){
  
  if(!is.data.frame(data) & !is.matrix(data)){print("the data must be a data.frame or a matrix")}
  
  normalization=toupper(normalization)
  if(identical(normalization, "MEAN")){
    data<-scale(log2(data))
  }
  else if(identical(normalization, "QUANTILE")){
    data <- log2(data)
    nomCol <- colnames(data)
    data <- as.data.frame(normalize.quantiles(as.matrix(data)), row.names = row.names(data))
    colnames(data) <- nomCol
  }
  else if(identical(normalization, "MEDIAN")){
    data<-log2(data)
    for (x in 1:length(data)){
      data[,x]<-(data[,x]-median(data[,x],na.rm = T))/(mad(data[,x],na.rm = T))
    }
  }
  else{
    stop("normalization must be : 'MEAN', 'QUANTILE' or 'MEDIAN'")
  }
  return(data)
}

test.stat <- function(data, stat = "t.test" , design , contrast){
  if(!is.data.frame(data) & !is.matrix(data)){print("the data must be a data.frame or a matrix")}
  if(!is.data.frame(design) & !is.matrix(design)){print("the experience design must be a matrix generated from model.matrix")}
  if(!is.data.frame(contrast) & !is.matrix(contrast)){print("the contrast matrix must be  a matrix generated from makeContrasts()")}
  
  
  nbComp<-ncol(contrast)
  listResults<-list()
  
  for (i in 1:nbComp){
    
    op<-row.names(contrast)[contrast[,i] == 1 ]
    ref<-row.names(contrast)[contrast[,i] == -1 ]
    
    samplesop <- row.names(design)[row(as.matrix(design[,op]))[design[,op]==1]]
    samplesref <- row.names(design)[row(as.matrix(design[,ref]))[design[,ref]==1]]
    
    colop <- which(colnames(data) %in% samplesop)
    colref <- which(colnames(data) %in% samplesref)
    
    p.value<-apply(data,1,function(x){t.test(as.numeric(x[colref]),as.numeric(x[colop]), alternative = "t") $p.value})
    adjust.p.value<-p.adjust(p.value, method = "BH")
    
    fc<-rowMeans(data[,colop])-rowMeans(data[,colref])
    
    listResults[[i]]<-data.frame(protein=c(row.names(data)))
    listResults[[i]][paste("p.value.",colnames(contrast)[i],sep="")]=c(p.value)
    listResults[[i]][paste("adjust.p.value.",colnames(contrast)[i],sep="")]=c(adjust.p.value)
    listResults[[i]][paste("fc.",colnames(contrast)[i],sep="")]=c(fc)
    
  }
  return(listResults)
}


shinyServer(function(input, output, session) {
  ##rendering the text input to enter the conditions name depending on the number of conditions 
  output$text<-renderUI({
    lapply(1:input$nbCond, function(i) {
      textInput(paste0('c', i), h5(paste0('name of condition ',i," : ")))
    })
  })
  
  output$select<-renderUI({
    lapply(1:input$nbCond, function(i) {
      selectizeInput(paste0('r', i), h5(paste0('replicates columns of confitions ',i," : ")), choices = list('please choose a dataset'="NULL"),multiple =T )
    })
  }) 
  
  lcol <- reactiveValues()
  observeEvent(input$dataFile,{
    lcol$name = colnames(read.table(input$dataFile$datapath, header=T, sep=",",row.names = 1))
    for(i in 1:input$nbCond){
      updateSelectizeInput(session, paste0('r', i), choices = setNames(c(1:length(lcol$name)), lcol$name) )
    }
  })
  observeEvent(input$nbCond,{
    for(i in 1:input$nbCond){
      if (!is.null(input$dataFile)){
        updateSelectizeInput(session, paste0('r', i), choices = setNames(c(1:length(lcol$name)), lcol$name) )
      }
    }
  })
  
  
  ## when the user clicks on the button submit :
  observeEvent(input$submit, {
    
    ##Hides the data selection panel and show the results 
    hide("dselect")
    show("cselect")
    show("plots")
    
    
    dataFile <- input$dataFile
    ##loading the file choosed by the user via the file input
    data <- read.table(dataFile$datapath, header=T, sep=",",row.names = 1)
    
    
    if(input$manual == FALSE){
      
      designFile <- input$designFile
      design <-as.matrix(read.csv(designFile$datapath, row.names = 1, check.names = F))
      contrastFile <- input$contrastFile
      contrast <- as.matrix(read.csv(contrastFile$datapath, row.names = 1, check.names = F))
      nbComp<-ncol(contrast)
      updateNumericInput(session, "comp",max = nbComp)
      
      norm.data <- norm(data, normalization = input$norm)
      results <- test.stat(norm.data, design = design, contrast = contrast)
      
      
      #rendering the data frame containing all the results for the selected comparison
      output$allProteins <- renderDataTable({
        as.data.frame(results[[input$comp]])
      })
      #rendering the button allowing to download all the results for the selected comparison
      output$dld <- downloadHandler(
        filename = function(){paste0("DEproteins.",substr(colnames(results[[input$comp]])[2], 9 ,nchar(colnames(results[[input$comp]])[2])), ".csv")},
        content = function(file) {
          write.csv(results[[input$comp]], file)
        }
      )
      
      ##rendering the Volcano plot for the selected comparison
      ggplotInput <- reactive({
        plotTitle <- substr(colnames(results[[input$comp]])[2], 9 ,nchar(colnames(results[[input$comp]])[2]))
        values <- as.data.frame(results[[input$comp]])
        forplot <- data.frame(x=as.numeric(values[,4]), y=-log10(values[,3]), id=as.character(values[,1]))
        tmp <- forplot[as.numeric(forplot$y)>=-log10(input$p) & abs(forplot$x)>input$fc,]
        p <- ggplot(forplot) + geom_point(aes(x, y , label = id, color = ifelse(y>=-log10(input$p) & abs(x)>=input$fc, "DE proteins", "not significant")),show.legend = F) +
          scale_color_manual(values = c("red", "blue")) +
          geom_text_repel(data = subset(forplot, abs(forplot$x)>=input$fc & forplot$y>=-log10(input$p)),
                          aes(x,y,label = id),
                          size = 4) +
          geom_vline(xintercept = input$fc ) +
          geom_vline(xintercept = -input$fc) + 
          geom_hline(yintercept = -log10(input$p)) + 
          labs(title = plotTitle,x="log2(Fold-change)", y="-log10(P.Value)") + theme_bw()  
      })
      
      output$volcanoPlot<- renderPlotly({
        ggplotly(ggplotInput(), tooltip = c("x","y","id")) %>% layout()
      })
      
      
      ##rendering the button allowing to download the Volcano plot for the selected comparison
      output$dlvp <- downloadHandler(
        filename = function() { paste0("VolcanoPlot-",substr(colnames(results[[input$comp]])[2], 9 ,nchar(colnames(results[[input$comp]])[2])), '.svg')},
        content = function(file) {
          plotTitle <- substr(colnames(results[[input$comp]])[2], 9 ,nchar(colnames(results[[input$comp]])[2]))
          values <- as.data.frame(results[[input$comp]])
          forplot <- data.frame(x=as.numeric(values[,4]), y=-log10(values[,3]), id=as.character(values[,1]))
          tmp <- forplot[as.numeric(forplot$y)>=-log10(input$p) & abs(forplot$x)>input$fc,]
          p <- ggplot(forplot) + geom_point(aes(x, y , color = ifelse(y>=-log10(input$p) & abs(x)>=input$fc, "not signi", "FC")),show.legend = F) +
            scale_color_manual(values = c("blue", "red")) +
            geom_text_repel(data = subset(forplot, abs(forplot$x)>=input$fc & forplot$y>=-log10(input$p)),
                            aes(x,y,label = id),
                            size = 2) +
            geom_vline(xintercept = input$fc ) +
            geom_vline(xintercept = -input$fc) + 
            geom_hline(yintercept = -log10(input$p)) + 
            labs(title = plotTitle,x="log2(Fold-change)", y="-log10(P.Value)") + theme_bw() 
          ggsave(file,plot = p)
        }
      )
      
      ##rendering the data frame containing the results for the proteins above the threshold (significant proteins) for the selected comparison
      output$significantProteins <- renderDataTable({
        fc = as.data.frame(results[[input$comp]])[,4]
        p = as.data.frame(results[[input$comp]])[,3]
        dt<-as.data.frame(results[[input$comp]])
        if(length(which(p<=input$p & abs(fc)>=input$fc))!=0){
          return(dt[which(p<=input$p & abs(fc)>=input$fc),])
        }
      })
      ##rendering the button allowing to dowload the results for the proteins above the threshold (significant proteins) for the selected comparison
      output$dlsd <- downloadHandler(
        filename = function(){ paste0("significant.DEproteins.",substr(colnames(results[[input$comp]])[2], 9 ,nchar(colnames(results[[input$comp]])[2])), ".csv")},
        content = function(file) {
          fc = as.data.frame(results[[input$comp]])[,4]
          p = as.data.frame(results[[input$comp]])[,3]
          dt<-as.data.frame(results[[input$comp]])
          write.csv(dt[which(p<=input$p & abs(fc)>=input$fc),], file)
        }
      )
      
      output$normalizedPlot <- renderPlot({
        boxplot(norm.data ,main = "boxplot of the normalized data",ylab="normalized Intensity)", xlab="samples",names=c(1:ncol(data)))
      })
      
    }
    
    
    
    
    
    if(input$manual == TRUE) {
      
      nbCond <- input$nbCond
      
      updateNumericInput(session, "cond1",max = nbCond)
      updateNumericInput(session, "cond2",max = nbCond)
      
      rep<-list()
      for(i in 1:nbCond){
        rep[[i]]<-eval(parse(text = paste0("input$r",i)))
      }
      
      
      
      ##Setting the names with the values of the text inputs
      names<-vector()
      for(i in 1:nbCond){
        if(eval(parse(text = paste0("input$c",i)))==""){
          names[i]<- as.character(i)
        }else {names[i] <- eval(parse(text = paste0("input$c",i)))}
      }
      
      if (is.null(dataFile))
        return(NULL)
      
      
      
      
      ##applying the DEprot function to the data
      tryCatch(results <- DEprot(data,nbCond, rep = rep,condName = names,normalization = as.character(input$norm)),
               warning = function(w){shinyjs::info(w)},
               error = function(e){
                 if (grepl("essentially constant", as.character(e))){
                   shinyjs::info("error : the chosen normalization technique makes the data too constant for a student's test. Try to run another normalization.")
                   show("dselect")
                   hide("cselect")
                   hide("plots")
                 }else if (grepl("not enough", as.character(e))){
                   shinyjs::info("error : please choose the right number of conditions and replicate. \nThere should be at least 2 replicates and two samples")
                   show("dselect")
                   hide("cselect")
                   hide("plots")
                 }else if (grepl("column", as.character(e))){
                   shinyjs::info(paste0("number of column in the file : ", ncol(data),"\nnumber of samples selected :", (nbCond), "\nplease enter a valid number of samples"))
                   show("dselect")
                   hide("cselect")
                   hide("plots")
                 }else{shinyjs::info(e)}
               })
      
      ##rendering the data frame containing all the results for the selected comparison
      output$allProteins <- renderDataTable({
        as.data.frame(results$dataprocessed[input$cond1,input$cond2])
      })
      ##rendering the button allowing to download all the results for the selected comparison
      output$dld <- downloadHandler(
        filename = function(){paste0("DEproteins.",names[input$cond1],".vs.",names[input$cond2], ".csv")},
        content = function(file) {
          write.csv(results$dataprocessed[input$cond1,input$cond2], file)
        }
      )
      
      ##rendering the Volcano plot for the selected comparison
      ggplotInput <- reactive({
        values <- as.data.frame(results$dataprocessed[input$cond1,input$cond2])
        forplot <- data.frame(x=as.numeric(values[,4]), y=-log10(values[,3]), id=as.character(values[,1]))
        tmp <- forplot[as.numeric(forplot$y)>=-log10(input$p) & abs(forplot$x)>input$fc,]
        p <- ggplot(forplot) + geom_point(aes(x, y, color = ifelse(y>=-log10(input$p) & abs(x)>=input$fc, "DE proteins", "not significant")),show.legend = F) +
          scale_color_manual(values = c("blue", "red")) +
          geom_text_repel(data = subset(forplot, abs(forplot$x)>=input$fc & forplot$y>=-log10(input$p)),
                          aes(x,y,label = id)) +
          geom_vline(xintercept = input$fc ) +
          geom_vline(xintercept = -input$fc) + 
          geom_hline(yintercept = -log10(input$p)) + 
          labs(title=paste0("VolcanoPlot.",names[input$cond1],".vs.",names[input$cond2]),x="log2(Fold-change)", y="-log10(P.Value)") + theme_bw() 
        
      })
      
      output$volcanoPlot<- renderPlotly({
        ggplotly(ggplotInput(), tooltip = c("x","y","id")) %>% layout()
      })
      
      
      ##rendering the button allowing to download the Volcano plot for the selected comparison
      output$dlvp <- downloadHandler(
        filename = function() { paste0("VolcanoPlot.",names[input$cond1],".vs.",names[input$cond2], '.svg')},
        content = function(file) {
          values <- as.data.frame(results$dataprocessed[input$cond1,input$cond2])
          forplot <- data.frame(x=as.numeric(values[,4]), y=-log10(values[,3]), id=as.character(values[,1]))
          tmp <- forplot[as.numeric(forplot$y)>=-log10(input$p) & abs(forplot$x)>input$fc,]
          p <- ggplot(forplot) + geom_point(aes(x, y, color = ifelse(y>=-log10(input$p) & abs(x)>=input$fc, "not signi", "FC")),show.legend = F) +
            scale_color_manual(values = c("blue", "red")) +
            geom_text_repel(data = subset(forplot, abs(forplot$x)>=input$fc & forplot$y>=-log10(input$p)),
                            aes(x,y,label = id),
                            size = 2) +
            geom_vline(xintercept = input$fc ) +
            geom_vline(xintercept = -input$fc) + 
            geom_hline(yintercept = -log10(input$p)) + 
            labs(title=paste0("VolcanoPlot.",names[input$cond1],".vs.",names[input$cond2]), x="log2(Fold-change)", y="-log10(P.Value)") +theme_bw() 
          ggsave(file,plot = p)
        }
      )
      ##rendering the data frame containing the results for the proteins above the threshold (significant proteins) for the selected comparison
      output$significantProteins <- renderDataTable({
        fc = as.data.frame(results$dataprocessed[input$cond1,input$cond2])[,4]
        p = as.data.frame(results$dataprocessed[input$cond1,input$cond2])[,3]
        dt<-as.data.frame(results$dataprocessed[input$cond1,input$cond2])
        if(length(which(p<=input$p & abs(fc)>=input$fc))!=0){
          return(dt[which(p<=input$p & abs(fc)>=input$fc),])
        }
      })
      ##rendering the button allowing to dowload the results for the proteins above the threshold (significant proteins) for the selected comparison
      output$dlsd <- downloadHandler(
        filename = function(){ paste0("significant.DEproteins.",names[input$cond1],".vs.",names[input$cond2], ".csv")},
        content = function(file) {
          fc = as.data.frame(results$dataprocessed[input$cond1,input$cond2])[,4]
          p = as.data.frame(results$dataprocessed[input$cond1,input$cond2])[,3]
          dt<-as.data.frame(results$dataprocessed[input$cond1,input$cond2])
          write.csv(dt[which(p<=input$p & abs(fc)>=input$fc),], file)
        }
      )
      
      output$normalizedPlot <- renderPlot({
        boxplot(results$normdata ,main = "boxplot of the normalized data",ylab="normalized Intensity)", xlab="samples",names=c(1:ncol(data)))
      })
      
    }
    
    
    
    ##rendering the boxplots to visualize the normalization of the data
    output$unnormalizedPlot <- renderPlot({
      boxplot(log2(data) ,main = "boxplot of the log2 data",ylab="log2(Intensity)", xlab="samples",names=c(1:ncol(data)))
    })
    
    
  })
  
  ##Hides the results and show the data selection panel when clicking on the "back" button
  observeEvent(input$back,{
    show("dselect")
    hide("cselect")
    hide("plots")
  })
  
})
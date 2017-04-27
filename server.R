library(shiny)

DEprot <- function(data,nbCond=NULL,nbRep=NULL,condName=NULL, normalization="SCALING"){
  
  
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
    cond[[i]]<-c(col:(col+nbRep-1))
    col<-col+nbRep
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

shinyServer(function(input, output, session) {

  ##rendering the text input to enter the conditions name depending on the number of conditions 
  output$text<-renderUI({
    lapply(1:input$nbCond, function(i) {
      textInput(paste0('c', i), paste0('name of condition ',i," : "))
    })
  })
  

  ## when the user clicks on the button submit :
  observeEvent(input$submit, {
    
    dataFile <- input$dataFile
    nbCond <- input$nbCond
    nbRep <- input$nbRep
    
    updateNumericInput(session, "cond1",max = nbCond)
    updateNumericInput(session, "cond2",max = nbCond)
    
    
    ##Hides the data selection panel and show the results 
    hide("dselect")
    show("cselect")
    show("plots")
    
    ##Setting the names with the values of the text inputs
    names<-vector()
    for(i in 1:nbCond){
      if(eval(parse(text = paste0("input$c",i)))==""){
        names[i]<- as.character(i)
      }else {names[i] <- eval(parse(text = paste0("input$c",i)))}
    }
    
    if (is.null(dataFile))
      return(NULL)
    
    ##loading the file choosed by the user via the file input
    data <- read.table(dataFile$datapath, header=T, sep=",",row.names = 1)
    
    
    ##applying the DEprot function to the data
    tryCatch(results <- DEprot(data,nbCond,nbRep,condName = names,normalization = as.character(input$norm)),
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
                 shinyjs::info(paste0("number of column in the file : ", ncol(data),"\nnumber of samnples selected :", (nbRep*nbCond), "\nplease enter a valid number of samples"))
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
      output$volcanoPlot <- renderPlot({
        fc = as.data.frame(results$dataprocessed[input$cond1,input$cond2])[,4]
        p = as.data.frame(results$dataprocessed[input$cond1,input$cond2])[,3]
        plot(fc,
             -log10(p),
             main = 'Volcano plot',
             ylab ="-log10(adjusted p value)",
             xlab = paste0("fold change ",names[input$cond1]," vs ",names[input$cond2]),
             abline(h = -log10(input$p), v = c(-input$fc,input$fc)),
             pch = 20,
             col = ifelse(p<=input$p & abs(fc)>=input$fc, "red", "blue"))
         N<-which(p<=input$p & abs(fc)>=input$fc)##threshold
        if(length(N)!=0){
          text(fc[N], -log10(p)[N], labels = as.data.frame(results$dataprocessed[input$cond1,input$cond2])[N,1], cex= 0.5, pos = 2)##adding text on plot
        }
      })
      ##rendering the button allowing to download the Volcano plot for the selected comparison
      output$dlvp <- downloadHandler(
        filename = function() { paste0("VolcanoPlot.",names[input$cond1],".vs.",names[input$cond2], '.svg')},
        content = function(file) {
          fc = as.data.frame(results$dataprocessed[input$cond1,input$cond2])[,4]
          p = as.data.frame(results$dataprocessed[input$cond1,input$cond2])[,3]
          svg(file)
          plot(fc,
               -log10(p),
               main = 'Volcano plot',
               ylab ="-log10(adjusted p value)",
               xlab = paste0("fold change ",names[input$cond1]," vs ",names[input$cond2]),
               abline(h = -log10(input$p), v = c(-input$fc,input$fc)),
               pch = 20,
               col = ifelse(p<=input$p & abs(fc)>=input$fc, "red", "blue"))
          N<-which(p<=input$p & abs(fc)>=input$fc)
          if(length(N)!=0){
            text(fc[N], -log10(p)[N], labels = as.data.frame(results$dataprocessed[input$cond1,input$cond2])[N,1], cex= 0.5, pos = 2)##row(as.matrix(fc))[N]
          }
          dev.off()
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
    
      
    ##rendering the boxplots to visualize the normalization of the data
    output$unnormalizedPlot <- renderPlot({
      boxplot(log2(data) ,main = "boxplot of the log2 data",ylab="log2(Intensity)", xlab="samples",names=c(1:(nbRep*nbCond)))
    })
    output$normalizedPlot <- renderPlot({
      boxplot(results$normdata ,main = "boxplot of the normalized data",ylab="normalized Intensity)", xlab="samples")#,names=c(1:(nbRep*nbCond))
    })
    
  })
  
  ##Hides the results and show the data selection panel when clicking on the "back" button
  observeEvent(input$back,{
    show("dselect")
    hide("cselect")
    hide("plots")
  })
  
})
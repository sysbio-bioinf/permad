my.plot <- function(result, name = "", ylim = c(0,105), cex.axis = 1.8)
{
    train <- result$train.performance
    test  <- result$test.performance
    colnames(train) <- colnames(test) <- c("Accuracy\n(overall)","Specificity\n(not in 100 d.)","Sensitivity\n(within 100 d.)","SS2\n(overall)")
    
    index <- c(1,3,2,4)
    
    train <- train[,index]
    test  <- test[,index]
    
    forg = c("#a6cee3","#b2df8a","#fb9a99","#fdbf6f")
    backg = c("#1f78b4","#33a02c","#e31a1c","#ff7f00")

    pdf(paste(name,".pdf",sep  =""), height = 9, width = 14)
    par(mar = c(1,1,1,1)*7)
    plot( 0,0, 
          type = "n",
          xlim = c(0,5),
          ylim = ylim,
          xaxs = 'i', 
          yaxs = 'i', 
          xlab  ="", 
          yaxt = "n", ylab = "",
          xaxt = "n")
    
    axis(2, seq(0,100,10), paste(seq(0,100,10),"%"), las = 2,  cex.axis = cex.axis)
    
    abline(h = seq(0,100,10), col = "gray")
    
    #boxplot(train*100, col = forg,  xaxt  ="n", yaxt = "n",add = T, at = c(1,4,7,10), box = "n")
    #boxplot(test*100,  col = backg, xaxt  ="n", yaxt = "n",add = T, at = c(2,5,8,11), box = "n")
    boxplot(test*100,  col = backg, xaxt  ="n", yaxt = "n",add = T, at = 1:4, box = "n")
    
    #axis(1, 1:4 , rep("TE",4), las = 0, cex.axis = cex.axis)
    #axis(1, c(1,2,4,5,7,8,10,11), rep(c("TR","TE"),4), las = 0)
    axis(3, 1:4, colnames(train), las = 0, cex.axis = cex.axis)
    dev.off()
    
    system(paste("pdfcrop ",name,".pdf", sep  =""))
    system(paste("mv ",name,"-crop.pdf ",name,".pdf", sep  =""))
    
    
     nms <- c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
    
    print(train)
    
    #sink(paste0(name,"_summary.txt"))
        #print("Train")
            tmp <- apply(train*100,2,function(x){
                round(as.vector(summary(x)), digits = 2)
            })
            rownames(tmp) <- nms
         #   print(tmp)
        #print("Test")
            tmp1 <- apply(test*100,2,function(x){
                round(as.vector(summary(x)), digits = 2)
            })
            rownames(tmp1) <- nms
            #print(tmp)
            
            tmp <- cbind(tmp,tmp1)
            
            write.csv(tmp, file = paste0(name,"_summary.csv") )
    #sink()
    
}



load("RESULT_TT_KNN.RData")
my.plot(result.knn, name = "TT_KNN", ylim = c(20,105) )


load("RESULT_TT_RF.RData")
my.plot(result.rf, name = "TT_RF",  ylim = c(40,105))


load("RESULT_TT_SVM.RData")
my.plot(result.svm, name = "TT_SVM",  ylim = c(0,105))



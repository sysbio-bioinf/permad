#####################################################################################
###
### inner.knn

inner.knn <- function(data, labs, ntimes=5,nfold=5, ks=c(1,3,5,7))
{
    require(doParallel)
    require(class)
    require(TunePareto)
    
    labs = as.character(labs)
    
    #######################################
    #generate foldList
    foldList <- generateCVRuns(labs, ntimes = ntimes, nfold = nfold, stratified = TRUE)
    
    #######################################
    #generate combs
    
    pred <- do.call("rbind",foreach::foreach(k=ks) %dopar% {
                mat <- sapply(foldList, function(run){
                    prediction <- unlist(lapply(run, function(fold){
                        knn(train = data[-fold,,drop=FALSE],
                            test  = data[fold,,drop=FALSE],
                            cl    = labs[-fold],
                            k     = k)
                    }))
    
                    rearranged.labs <- labs[unlist(run)]
                   
                    index.one  <- rearranged.labs == "1"
                    index.two  <- rearranged.labs == "2"
                
                    correct <- prediction == rearranged.labs
                
                    return( c(  mean(correct), #ACC
                                mean(correct[index.one]), #Sens 1
                                mean(correct[index.two]), #Sens 2
                                (mean(correct[index.one])+mean(correct[index.two]))/2 ) # SS2
                            )
                })
                
                if(!is.matrix(mat))
                    mat <- matrix(mat, nrow = 4, ncol = ntimes)
                
                rownames(mat) <- c("Acc","Sens1","Sens2","SS2")
                colnames(mat) <- paste("Run", 1:ntimes,sep = "")
                
                return(c(mean(mat["Acc",]),
                         mean(mat["Sens1",]),
                         mean(mat["Sens2",]),
                         mean(mat["SS2",]),
                         mat["Acc",],
                         mat["Sens1",],
                         mat["Sens2",],
                         mat["SS2",]))
               
    })

    colnames(pred) <- c(    "Acc","Sens1","Sens2","SS2",
                            paste("Acc_Run", 1:ntimes, sep = ""),
                            paste("Sens1_Run", 1:ntimes, sep = ""),
                            paste("Sens2_Run", 1:ntimes, sep = ""),
                            paste("SS1_Run", 1:ntimes, sep = ""))
                            
    rownames(pred) <- as.character(ks)
    
    pred <- pred[order(pred[,"SS2"], decreasing = TRUE)[1],,drop = FALSE]
    
    return(pred)
}

#####################################################################################
###
### outer.knn

outer.knn <- function(dataset, foldList, itimes = 5, ifold=5, ks=c(1,3,5,7) ){
    require(class)
    require(doParallel)
    require(TunePareto)
    
    numRun <- length(foldList)
    numFol <- length(foldList[[1]])

    result <- lapply(foldList, function(run){
        lapply(run, function(fold){
    
            #######################################
            #generate test samples
    
            test.samples <- fold
    
            train.data <- t(dataset$data[,-test.samples, drop = FALSE])
            test.data  <- t(dataset$data[,test.samples, drop = FALSE])

            train.labs <- as.factor(dataset$labs)[-test.samples, drop = FALSE]
            test.labs  <- as.factor(dataset$labs)[test.samples, drop = FALSE]

            index.one  <- test.labs == "1"
            index.two  <- test.labs == "2"
    
            #######################################
            #inner CV
            in.cv <- inner.knn(train.data, train.labs, ntimes = itimes, nfold=ifold, ks=ks)
            k <-as.numeric(rownames(in.cv))
            names(k) <- "k"
    
            prediction = knn(   train = train.data,
                                test  = test.data,
                                cl    = train.labs,
                                k     = k)
                        
            names(prediction) = test.samples

            correct <- prediction == test.labs
                
            pred <- matrix( c(  mean(correct), #ACC
                                mean(correct[index.one]), #Sens 1
                                mean(correct[index.two]), #Sens 2
                                (mean(correct[index.one])+mean(correct[index.two]))/2 ) # SS2
                                , nrow=1)

            colnames(pred) <- c("Acc","Sens1","Sens2","SS2")
            rownames(pred) <- as.character(k)
    
            list(   test.performance = pred,
                    train.performance= in.cv[,1:4,drop=FALSE],
                    params           = k,
                    test.prediction  = prediction)
        })
    })
    
    tasks <- expand.grid(1:numFol,1:numRun)[,2:1]
    
    nms <- apply(tasks, 1, function(x){paste("R",x[1],"F",x[2],sep="")})
    
    test.performance <- do.call("rbind",lapply(result, function(run){
        mat <- t(sapply(run,function(fold){
            fold$test.performance
        }))
        colnames(mat) <- c("Acc","Sens1","Sens2","SS2")
        
        return(mat)
    }))
    
    train.performance <- do.call("rbind",lapply(result, function(run){
        mat <- t(sapply(run,function(fold){
            fold$train.performance
        }))
        colnames(mat) <- c("Acc","Sens1","Sens2","SS2")
        
        return(mat)
    }))
    
    params <- lapply(1:nrow(tasks), function(i){
        result[[tasks[i,1]]][[tasks[i,2]]]$params
    })
    
    test.prediction <- lapply(1:nrow(tasks), function(i){
        result[[tasks[i,1]]][[tasks[i,2]]]$test.prediction
    })
    
    rownames(train.performance) <- nms
    rownames(test.performance)  <- nms
    names(params)               <- nms
    names(test.prediction)      <- nms
    
    return(list(train.performance=train.performance,
                test.performance = test.performance,
                test.prediction  = test.prediction,
                params = params,
                experimentName = "cvKNN",
                dataName = dataset$name))
}


#####################################################################################
###
### inner.svm

inner.svm <- function(data, labs, ntimes=5,nfold=5, costs=c(0.01,0.1,1,10,100))
{
    require(doParallel)
    require(e1071)
    require(TunePareto)
    
    labs = as.character(labs)
    
    #######################################
    #generate foldList
    foldList <- generateCVRuns(labs, ntimes = ntimes, nfold = nfold, stratified = TRUE)
    
    #######################################
    #generate combs
    
    pred <- do.call("rbind",foreach::foreach(cost=costs) %dopar% {
                mat <- sapply(foldList, function(run){
                    prediction <- unlist(lapply(run, function(fold){
                        model <- svm(x = data[ -fold,,drop=FALSE],
                                     y = labs[ -fold],
                                     type = "C-classification",
                                     kernel = "linear",
                                     scale = FALSE,
                                     cost   = cost)
                            
                        predict(model, data[fold,,drop=FALSE])
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
                            
    rownames(pred) <- as.character(costs)
    
    pred <- pred[order(pred[,"SS2"], decreasing = TRUE)[1],,drop = FALSE]
    
    return(pred)
}

#####################################################################################
###
### outer.svm

outer.svm <- function(dataset, foldList, itimes = 5, ifold=5, costs=c(0.01,0.1,1,10,100) ){
    require(e1071)
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
            in.cv <- inner.svm(train.data, train.labs, ntimes = itimes, nfold=ifold, costs=costs)
                   
            cost <-as.numeric(rownames(in.cv))
            names(cost) <- "cost"
                                
            model <- svm(x = train.data,
                         y = train.labs,
                         type = "C-classification",
                         kernel = "linear",
                         scale = FALSE,
                         cost   = cost)
                                    
            prediction <- predict(model, test.data)
            
            names(prediction) = test.samples

            correct <- prediction == test.labs
                
            pred <- matrix( c(  mean(correct), #ACC
                                mean(correct[index.one]), #Sens 1
                                mean(correct[index.two]), #Sens 2
                                (mean(correct[index.one])+mean(correct[index.two]))/2 ) # SS2
                                , nrow=1)
            
    

            colnames(pred) <- c("Acc","Sens1","Sens2","SS2")
            rownames(pred) <- as.character(cost)
    
            list(   test.performance = pred,
                    train.performance= in.cv[,1:4,drop=FALSE],
                    params           = cost,
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
                experimentName = "cvSVM",
                dataName = dataset$name))
}


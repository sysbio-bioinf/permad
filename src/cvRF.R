#####################################################################################
###
### inner.svm

inner.rf <- function(data, labs, ntimes=5,nfold=5, ntree = seq(25,2000,by=25), nodesize = 1:5, classwt = seq(0.05, 0.95, by = 0.05))
{
    require(doParallel)
    require(randomForest)
    require(TunePareto)
    
    set.seed(123)
    
    labs = as.factor(labs)
    
    #######################################
    #generate foldList
    foldList <- generateCVRuns(labs, ntimes = ntimes, nfold = nfold, stratified = TRUE)
    
    #######################################
    #generate combs
    combs <- expand.grid(ntree = ntree, nodesize = nodesize, classwt = classwt)
    
    pred <- do.call("rbind",foreach::foreach(i =1:nrow(combs)) %dopar% {
                mat <- sapply(foldList, function(run){
                    prediction <- unlist(lapply(run, function(fold){
                        
                        set.seed(123)
                        
                        model <- randomForest(  x=data[-fold,,drop=FALSE],
                                                y=labs[-fold],
                                                ntree=combs[i,"ntree"],
                                                classwt=c(combs[i, "classwt"], 1-combs[i, "classwt"]),
                                                nodesize = combs[i, "nodesize"])
                        
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
                            
    rownames(pred) <- apply(combs,1,function(x){paste(x, collapse = "_")})
    
    pred <- pred[order(pred[,"SS2"], decreasing = TRUE)[1],,drop = FALSE]
    
    return(pred)
}

#####################################################################################
###
### outer.svm

outer.rf <- function(dataset, foldList, itimes = 5, ifold=5, ntree = seq(25,2000,by=25), nodesize = 1:5, classwt = seq(0.05, 0.95, by = 0.05) ){
    require(randomForest)
    require(doParallel)
    require(TunePareto)
    
    numRun <- length(foldList)
    numFol <- length(foldList[[1]])

    result <- lapply(foldList, function(run){
        lapply(run, function(fold){
            print("Tut")
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
            in.cv <- inner.rf(train.data, train.labs, ntimes = itimes, nfold=ifold, ntree = ntree, nodesize = nodesize, classwt = classwt)
                   
            cv <- rownames(in.cv)
            cv <-as.numeric(strsplit(cv, "_")[[1]])
                   
            names(cv) <- c("ntree", "nodesize", "classwt")
                                
            model <- randomForest(  x=train.data,
                                    y=train.labs,
                                    ntree=cv["ntree"],
                                    classwt=c(unname(cv["classwt"]), unname(1-cv["classwt"])),
                                    nodesize = cv["nodesize"],
                                    importance = TRUE)
                                    
            prediction <- predict(model, test.data)
            
            names(prediction) = test.samples

            correct <- prediction == test.labs
                
            pred <- matrix( c(  mean(correct), #ACC
                                mean(correct[index.one]), #Sens 1
                                mean(correct[index.two]), #Sens 2
                                (mean(correct[index.one])+mean(correct[index.two]))/2 ) # SS2
                                , nrow=1)
            
    

            colnames(pred) <- c("Acc","Sens1","Sens2","SS2")
            rownames(pred) <- paste(cv, collapse = "_")
    
            list(   test.performance = pred,
                    train.performance= in.cv[,1:4,drop=FALSE],
                    params           = cv,
                    test.prediction  = prediction,
                    importance       = importance(model,type = 2))
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
    
    importance <- do.call("cbind",lapply(1:nrow(tasks), function(i){
        result[[tasks[i,1]]][[tasks[i,2]]]$importance
    }))

    rownames(train.performance) <- nms
    rownames(test.performance)  <- nms
    names(params)               <- nms
    names(test.prediction)      <- nms
    colnames(importance)        <- nms
    
    return(list(train.performance=train.performance,
                test.performance = test.performance,
                test.prediction  = test.prediction,
                params = params,
                importance = importance,
                experimentName = "cvRF",
                dataName = dataset$name))
}


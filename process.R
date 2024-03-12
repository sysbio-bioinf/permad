## Simualtions for PERMAD data
## KNN runs on a MacBook Pro M1 Max, 64GB RAM in a couple of seconds
## SVM runs on a 16 core server in a couple of hours
## RFF takes a week on a 96 core AMD Epyc server

## load data
load("data/PERMAD.RData")
data <- t(PERMAD$data)
labs <- PERMAD$labs

## create foldList
library(TunePareto)
set.seed(54321)
ntimes = 5
nfold  = 10
foldList <- generateCVRuns(PERMAD$labs, ntimes = ntimes, nfold = nfold, stratified = TRUE)

## KNN
## initiate parallel cores
library(doParallel)
numCpu = 100
cl <- makeCluster(numCpu)
registerDoParallel(cores = cl)

## KNN
source("src/cvKNN.R")
result.knn <- outer.knn(PERMAD, foldList, itimes = ntimes,  ifold=nfold)
save(result.knn, file = "results/RESULT_KNN.RData")

## SVM
source("src/cvSVM.R")
result.svm <- outer.svm(PERMAD, foldList, itimes = ntimes, ifold=nfold)
save(result.svm, file = "results/RESULT_SVM.RData")

## RF
source("src/cvRF.R")
result.rf <- outer.rf(PERMAD, foldList, itimes = ntimes, ifold=nfold)
save(result.rf, file = "results/RESULT_RF.RData")

## RF importance ranking
importance.rank <- apply(result.rf$importance,2,rank)

## RF again on top10 features
top10 <- rownames(importance.rank)[1:10]
PERMAD.top10 <- PERMAD
PERMAD.top10$data <- PERMAD.top10$data[match(top10, rownames(PERMAD$data)),]
PERMAD.top10$features <- PERMAD$features[match(top10, PERMAD$features)]
PERMAD.top10$name <- "PERMAD_top10"
save(PERMAD.top10, file="data/PERMAD_top10.RData")
result.rf.top10 <- outer.rf(PERMAD.top10, foldList, itimes = ntimes, ifold=nfold)
save(result.rf.top10, file = "results/RESULT_RF_top10.RData")

## stop parallel cores
stopCluster(cl)

## RF top10 all samples reclass
model <- randomForest(x=t(PERMAD.top10$data), y=as.factor(PERMAD$labs), nodesize=1, ntree=41)
table(predict(model, t(PERMAD.top10$data)))

load("RESULT_RF.RData")
#result.rf

importance.rank <- apply(result.rf$importance,2,rank)




importance.rank <- (importance.rank[order(print(apply(importance.rank,1,median)), decreasing = TRUE),,drop=FALSE]-1)/nrow(importance.rank)
#importance.rank <- (importance.rank[order(print(apply(importance.rank,1,mean)), decreasing = TRUE),,drop=FALSE]-1)/nrow(importance.rank)
#print(importance.rank)
  
  
pdf("Importance-Ranking.pdf", width = 16, height = 12)
par(mar = c(5,5,25,5))
plot(   x=0,y=0,
        xlim = c(0, nrow(importance.rank))+c(-1,+1),
        ylim = c(1, 0),
        type = "n", las = 2, ylab = "Importance ranking", xaxt = "n",
        #main = "Importance ranking (RF-Experiments, 10x10 CV)", 
        yaxs = 'i',
xaxs = 'i', xlab  ="", yaxt = "n")
axis(2, seq(0,1,0.1),seq(0,1,0.1),las = 0)

nr <- 1:nrow(importance.rank)
axis(3, nr, paste(nr,". ",rownames(importance.rank), sep = ""),las = 2, cex.axis = 0.8)

abline(h=seq(0,1,0.1), col = "gray")
abline(v=seq(10.5,200,10), col = "gray")
boxplot(t(importance.rank), las  = 2, col = "orange", cex.axis = 0.5, add = TRUE, yaxt = "n", xaxt = "n", at=1:nrow(importance.rank) )
dev.off()

system(paste("pdfcrop Importance-Ranking.pdf", sep  =""))
system(paste("mv Importance-Ranking-crop.pdf Importance-Ranking.pdf", sep  =""))

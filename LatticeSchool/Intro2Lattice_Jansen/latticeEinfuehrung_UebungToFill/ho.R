library(boot)

######### analysis ###################

hodata <- read.table("Routput", skip=1, col.names=c("npath","x","x2","acceptance"),colClasses =c("integer","double","double","double"))

xmean <- mean(hodata$x)
xse <- sqrt(var(hodata$x)/length((hodata$x-1)))
print(paste0("<x>= ", xmean))
#sprintf("<x>= %14.7e\n", xmean)
print(paste0("d<x>= ", xse))

xmean2 <- mean(hodata$x2)
xse2 <- sqrt(var(hodata$x2)/length((hodata$x2-1)))
print(paste0("<x^2>= ", xmean2))
print(paste0("d<x^2>= ", xse2))

accmean <- mean(hodata$acceptance)
accse <- sqrt(var(hodata$acceptance)/length((hodata$acceptance-1)))
print(paste0("Acceptance= ", accmean))
print(paste0("dAcceptance= ", accse))

#
set.seed(1717)
#
bmeanx <- function(data,indices){
         d <- data[indices,]
         average <- mean(d$x)
         return(average)
}
xboot <- boot(data=hodata, statistic=bmeanx, R=1000)
#print(xboot)   
#print(xboot$t0)
#print(xboot$t)
cat("<x> = ",xboot$t0," Standard deviation = ", sd(xboot$t[,1]), "\n")
ci <- boot.ci(xboot, type="basic")
cat("95% CI from ", ci$basic[1,4], " - ", ci$basic[1,5],"\n")

bmeanx2 <- function(data,indices){
           d <- data[indices,]
           average <- mean(d$x2)
           return(average)
}
xboot <- boot(data=hodata, statistic=bmeanx2, R=1000)
cat("<x^2> = ",xboot$t0," Standard deviation = ", sd(xboot$t[,1]), "\n")
ci = boot.ci(xboot, type="basic")
cat("95% CI from ", ci$basic[1,4], " - ", ci$basic[1,5], "\n")




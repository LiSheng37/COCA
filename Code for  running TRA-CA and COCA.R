########## load function and data, set path====================================
library(rsatscan)
library(dplyr)
library(ggplot2)
library(cowplot)
library(sp)
setwd("E:/COCA/cov")
load(file = "Simulationdataset.RData")
###TRA-CA+COCA####
js <- c(1,4,7,10,13,16,19,22,25)
for(j in js){
  #1 ycluster，j is 1，4，7，10，13，16，19，22，25
  Rcluster2 <- c("36047","36081","34017","36061","36085","36005","34013","34039")
  #2 ycluster
  #Rcluster2 <- c("25025","25021","25009","25017","25023","25005","33015","44007","36047", 
  #               "36081","34017","36061","36085","36005","34013","34039")
  #3 ycluster
  #Rcluster2 <- c("25025","25021","25009","25017","25023","25005","33015","44007","36047","36081","34017", "36061",
  #               "36085","36005","34013","34039","36065","36053","36043","36049","36067","36075","36077","36017")
simcase <- simlist.beta1[[j]]
xsim <- xlist.beta1[[j]]
{
simudata <- read.table("simdata.cas",stringsAsFactors = F,header = F,
                       colClasses = c("character","numeric","numeric","numeric","character",
                                      "numeric","numeric","numeric","numeric","character",
                                      "numeric","numeric"))
colnames(simudata) <- c("code","x1","x","var","name","year","pop","expect","simupatient","name","geo.x","geo.y")
simudata <- simudata[,c(1,3,5,6,7,8,9,11,12)]
simudata2 <- simudata
casepop1 <- read.table("Cases.cas",stringsAsFactors = F,header = F,
                       colClasses = c("character","numeric"))
expectpop1 <- read.table("Population.pop",stringsAsFactors = F,header = F,
                         colClasses = c("character","character","numeric"))
expectadjpop1 <- read.table("Population.pop",stringsAsFactors = F,header = F,
                            colClasses = c("character","character","numeric"))
}
res <- data.frame(matrix(nrow=0,ncol=13))
simutime <- 2000
colnames(res) <- c("sensi1","sensi2","speci1","speci2","ppv1","ppv2","ydi1","ydi2","misclass1","misclass2","power1","power2","beteadj")
for(i in 1:simutime){
  case <- simcase1[,i]
  casepop1$V2 <- case
  write.table(casepop1,file = "Cases.cas",row.names = F,col.names = F)
  simudata$x  <- xsim[,i]
  model <- glm(case ~ offset(log(simudata$pop)) + simudata$x,family = poisson())
  expectpop1$V3 <- model$fitted.values 
  write.table(expectpop1,file = "expect1.pop",row.names = F,col.names = F)
  td <- setwd("E:/COCA/cov") #the location of data
  sslocation <- "C:\\Program Files\\SaTScan" ## the installation location of the SaTScan software
  invisible(ss.options(reset=TRUE)) 
  ss.options(list(CaseFile="Cases.cas",StartDate="2000/1/1",EndDate="2000/12/31", 
                  CoordinatesFile="Coordinates.geo", PopulationFile="expect1.pop",
                  PrecisionCaseTimes=0,CoordinatesType=0, AnalysisType=1, ModelType=0))
  #adjust the file of population
  ss.options(c("NonCompactnessPenalty=0", "ReportGiniClusters=n", "LogRunToHistoryFile=n"))
  write.ss.prm(td,"xsimu0706") 
  cluster <- satscan(td, "xsimu0706", sslocation="C:\\Program Files\\SaTScan")
  clua1 <- cluster[['gis']]
  clu1 <- subset(clua1,P_VALUE<=0.05)
  ##calculate a
  cluboth <- intersect(Rcluster2,clu1$LOC_ID)
  cluboth<- as.data.frame(cluboth)
  colnames(cluboth) <- c("code")
  clubothpop <- merge(cluboth,simudata) 
  apop <- sum(clubothpop$pop)
  #calculate b
  clutrue <- intersect(Rcluster2,setdiff(Rcluster2,clu1$LOC_ID))
  clutrue <- as.data.frame(clutrue)
  colnames(clutrue) <- c("code")
  clutruepop <- merge(clutrue,simudata)
  bpop <- sum(clutruepop$pop)
  #calculate c
  cludetec <- intersect(clu1$LOC_ID,setdiff(clu1$LOC_ID,Rcluster2))
  cludetec <- as.data.frame(cludetec)
  colnames(cludetec) <- c("code")
  cludetecpop <- merge(cludetec,simudata) 
  cpop <- sum(cludetecpop$pop)
  #calculate d
  cluno <- setdiff(simudata$code,union(Rcluster2,clu1$LOC_ID))
  cluno <- as.data.frame(cluno)
  colnames(cluno) <- c("code")
  clunopop <- merge(simudata,cluno)
  dpop <- sum(clunopop$pop)
  #performance
  sensi1 <- apop/(apop+bpop)
  speci1 <- dpop/(cpop+dpop)
  ppv1 <- apop/(apop+cpop)
  ydi1 <- sensi1+speci1-1
  misclass1 <- (bpop+cpop)/(apop+bpop+cpop+dpop)
  res[i,1] <- sensi1 
  res[i,3] <- speci1 
  res[i,5] <- ppv1 
  res[i,7] <- ydi1
  res[i,9] <- misclass1
  if(nrow(clu1)>0){
    res[i,11] <- 1
  } else{
    res[i,11] <- 0
  }
  ########remove clustery in TRA-CA
  ycluster <- subset(clua1,P_VALUE<=0.05)
  yclustercode <- ycluster$LOC_ID
  Fclustercod <- !simudata$cod %in% yclustercode
  fcase <- case[Fclustercod];fpop <- simudata$pop[Fclustercod];fx <- simudata$x[Fclustercod]
  fmodel <- glm(fcase ~ offset(log(fpop)) + fx,family = poisson())
  beta2 <- fmodel$coefficients["fx"]
  beta1 <- model$coefficients["simudata$x"]
  expectadjpop1$V3 <- predict(fmodel,newdata = data.frame(fpop = simudata$pop,fx = simudata$x),type = "response")
  write.table(expectadjpop1,file = "expectadj.pop",row.names = F,col.names = F)
  td <- setwd("E:/COCA/cov")
  sslocation <- "C:\\Program Files\\SaTScan" ## the installation location of the SaTScan software  
  invisible(ss.options(reset=TRUE)) 
  ss.options(list(CaseFile="Cases.cas",StartDate="2000/1/1",EndDate="2000/12/31", 
                  CoordinatesFile="Coordinates.geo", PopulationFile="expectadj.pop",
                  PrecisionCaseTimes=0,CoordinatesType=0, AnalysisType=1, ModelType=0))
  ss.options(c("NonCompactnessPenalty=0", "ReportGiniClusters=n", "LogRunToHistoryFile=n"))
  write.ss.prm(td,"xsimuadj") 
  cluster2 <- satscan(td, "xsimuadj", sslocation="C:\\Program Files\\SaTScan")
  clu2 <- cluster2[['gis']]
  clu2 <- subset(clu2,P_VALUE<=0.05)
  ##calculate a
  cluboth2 <- intersect(Rcluster2,clu2$LOC_ID)
  cluboth2 <- as.data.frame(cluboth2)
  colnames(cluboth2) <- c("code")
  clubothpop2 <- merge(cluboth2,simudata2) 
  apop2 <- sum(clubothpop2$pop)
  #calculate b
  clutrue2 <- intersect(Rcluster2,setdiff(Rcluster2,clu2$LOC_ID))
  clutrue2 <- as.data.frame(clutrue2)
  colnames(clutrue2) <- c("code")
  clutruepop2 <- merge(clutrue2,simudata2)
  bpop2 <- sum(clutruepop2$pop)
  #calculate c
  cludetec2 <- intersect(clu2$LOC_ID,setdiff(clu2$LOC_ID,Rcluster2))
  cludetec2 <- as.data.frame(cludetec2)
  colnames(cludetec2) <- c("code")
  cludetecpop2 <- merge(cludetec2,simudata2) 
  cpop2 <- sum(cludetecpop2$pop)
  #calculate d
  cluno2 <- setdiff(simudata2$code,union(Rcluster2,clu2$LOC_ID))
  cluno2 <- as.data.frame(cluno2)
  colnames(cluno2) <- c("code")
  clunopop2 <- merge(simudata2,cluno2)
  dpop2 <- sum(clunopop2$pop)
  #指标
  sensi2 <- apop2/(apop2+bpop2)
  speci2 <- dpop2/(cpop2+dpop2)
  ppv2 <- apop2/(apop2+cpop2)
  ydi2 <- sensi2+speci2-1
  misclass2 <- (bpop2+cpop2)/(apop2+bpop2+cpop2+dpop2)
  res[i,2] <- sensi2 
  res[i,4] <- speci2 
  res[i,6] <- ppv2 
  res[i,8] <- ydi2
  res[i,10] <- misclass2
  res[i,13] <- beta2
  if(nrow(clu2)>0){
    res[i,12] <- 1
  } else{
    res[i,12] <- 0
  }
}
res
ress <- colMeans(res,na.rm = TRUE)
f1 <- paste("TRA-CA and COCA ","res",substring(names(xlist[j]),2,8),"1beta",".csv",sep = "")
write.csv(res,f1)
f2 <- paste("TRA-CA and COCA ","meanres",substring(names(xlist[j]),2,8),"1beta",".csv",sep = "")
write.csv(ress,f2)
}


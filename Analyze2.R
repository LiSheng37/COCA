### Load the R packages to enable multithreading
library('parallel') 

### the function to complete the simulation study for COCA.
test_funciton <- function(j){
  ### Load the required R packages
  library("DClusterm")
  library("sp")
  library("dplyr")
  
  setwd("E:/COCA/cov") ### set the working directory to the location of data. All data is contained in the cov folder. 
  
  ### When there is 1 true cluster, j is 1，4，7，10，13，16，19，22，25
  ### The location codes of spatial units contained in the true cluster.
  Rcluster2 <- c("36047","36081","34017","36061","36085","36005","34013","34039")
  
  ### When there are 2 true clusters, j is 2，5，8，11，14，17，20，23，26
  ### The location codes of spatial units contained in the true cluster.
  # Rcluster2 <- c("25025","25021","25009","25017","25023","25005","33015","44007","36047", 
  #               "36081","34017","36061","36085","36005","34013","34039")
  
  ### When there are 3 true clusters, j is 3，6，9，12，15，18，21，24，27
  ### The location codes of spatial units contained in the true cluster.
  # Rcluster2 <- c("25025","25021","25009","25017","25023","25005","33015","44007","36047","36081","34017", "36061",
  #               "36085","36005","34013","34039","36065","36053","36043","36049","36067","36075","36077","36017")
  
  ### Read the real dataset to ensure that the format of the simulated data matches that of the real dataset.
  ap <- st_read("COUNTIES.SHP")
  map <- as(ap, "Spatial")
  map$x <- coordinates(map)[, 1]
  map$y <- coordinates(map)[, 2]
  map2 <- data.frame(map)
  map2 <- map2[,c(1,5)]
  colnames(map2) <- c("name","code")
  load(file = "Simulationdataset.RData")
  {
    simudata <- read.table("simdata.cas",stringsAsFactors = F,header = F,
                           colClasses = c("character","numeric","numeric","numeric","character",
                                          "numeric","numeric","numeric","numeric","character",
                                          "numeric","numeric"))
    colnames(simudata) <- c("code","x1","x","var","name","year","pop","expect","simupatient","name","geo.x","geo.y")
    simudata <- simudata[,c(1,3,5,6,7,8,9,11,12)]
    simudata2 <- simudata
    simu <- simudata[,c(1,3,5)]
  }
  
  ### read the simulated dataset corresponding to a simulation scenario
  simcase1 <- simlist.beta1[[j]]  ### Simulated case data
  xsim <- xlist.beta1[[j]] ### Simulated covariate data
  
  ### Define a data frame to save the results.
  res <- data.frame(matrix(nrow=0,ncol=6))
  colnames(res) <- c("sensi","speci","ppv","ydi","misclass","power")
  
  ### Given the simulation scenario, repeat the simulation 2000 times.
  simutime <- 2001
  for(i in 2:simutime){
    ### Read i-th simulated dataset and initially preprocess the data.
    simu$simupatient <-  simcase1[,i]
    simu$observed <- round(simu$simupatient)
    simu$Expected <- simu$pop * sum(simu$observed) / sum(simu$pop)
    simu$x1  <- xsim[,i]
    map3 <- left_join(map2,simu)
    map$pop <- map3$pop
    map$x1 <- map3$x1
    map$Observed <- map3$observed
    map$Expected <- map3$Expected
    
    ### Run DGLM
    ny.m0 <- glm(Observed ~ offset(log(Expected)) + x1, family = "poisson",data = map)
    ny.cl0 <- DetectClustersModel(map, thegrid = as.data.frame(map)[c("x", "y")], fractpop = 0.50, alpha = 0.05, radius = Inf,
                                  step = NULL, typeCluster = "S", R = 999, model0 = ny.m0,
                                  ClusterSizeContribution = "pop")
    if(ny.cl0 == "No clusters found"){
      m2 <- data.frame(code= "00000")
    }else{
      hh <- slimknclusters(map, ny.cl0, 1)
      map$clu <- get.allknclusters(map, hh)
      m <- map[which(map$clu == "CLUSTER"),]
      m1 <- data.frame(m) 
      m2 <- data.frame(m1[,5])
      colnames(m2) <- "code"
    }
    
    ### process the results from SaTScan and calculate evaluation indices for DGLM
    ### calculate a
    cluboth3 <- intersect(Rcluster2,m2$code)
    cluboth3<- as.data.frame(cluboth3)
    colnames(cluboth3) <- c("code")
    clubothpop3 <- merge(cluboth3,simudata) 
    apop3 <- sum(clubothpop3$pop)
    
    ### calculate b
    clutrue3 <- intersect(Rcluster2,setdiff(Rcluster2,m2$code))
    clutrue3 <- as.data.frame(clutrue3)
    colnames(clutrue3) <- c("code")
    clutruepop3 <- merge(clutrue3,simudata)
    bpop3 <- sum(clutruepop3$pop)
    
    ### calculate c
    cludetec3 <- intersect(m2$code,setdiff(m2$code,Rcluster2))
    cludetec3 <- as.data.frame(cludetec3)
    colnames(cludetec3) <- c("code")
    cludetecpop3 <- merge(cludetec3,simudata) 
    cpop3 <- sum(cludetecpop3$pop)
    
    ### calculate d
    cluno3 <- setdiff(simudata$code,union(Rcluster2,m2$code))
    cluno3 <- as.data.frame(cluno3)
    colnames(cluno3) <- c("code")
    clunopop3 <- merge(simudata,cluno3)
    dpop3 <- sum(clunopop3$pop)
    
    ### calculate evaluation index
    sensi3 <- apop3/(apop3+bpop3)
    speci3 <- dpop3/(cpop3+dpop3)
    ppv3 <- apop3/(apop3+cpop3)
    ydi3 <- sensi3+speci3-1
    misclass3 <- (bpop3+cpop3)/(apop3+bpop3+cpop3+dpop3)
    
    ### store evaluation index values 
    res[i,1] <- sensi3 
    res[i,2] <- speci3 
    res[i,3] <- ppv3 
    res[i,4] <- ydi3
    res[i,5] <- misclass3
    if(ny.cl0 == "No clusters found"){
      res[i,6] <- 0
    } else{
      res[i,6] <- 1
    }
    print(i)
  }
  res
  ress <- colMeans(res,na.rm = TRUE)
  f1 <- paste("DGLM","Dres",substring(names(xlist[j]),2,8),"1beta",".csv",sep = "")
  write.csv(res,f1)
  f2 <- paste("DGLM","meanDres",substring(names(xlist[j]),2,8),"1beta",".csv",sep = "")
  write.csv(ress,f2)
}

### Complete the simulation study of DLGM using multithreading
detectCores() ### Check the number of cores on the computer for running the code
cl <- makeCluster(9) ### Set the number of cores for running the code to 9
system.time({ ### Calculate the runtime
  res <- parLapply(cl, c(1,4,7,10,13,16,19,22,25),test_funciton)
})
stopCluster(cl) ### Release the memory


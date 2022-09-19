####################################################
# Author: Young Won Cho
# Date: 2022-08-30
# Purpose: Data simulation for univariate OSC model and
#          Missing data handling
####################################################

rm(list=ls(all=TRUE)) #To remove all stored objects
username<-Sys.info()[7]

# Packages ----
library(deSolve)
library(Matrix)   # Containing data generating function for a state-space model
library(ggplot2)
require(plyr)     # for 'ddply' function
library(dplyr)
library(MASS)
library(reshape2) # for long to wide format
library(foreach)  # Use foreach to estimate models in parallel
library(doParallel)
library(mice)

# Functions ----
# 1. Time delay embedding function
Embed = function(x, E, tau) {  
  len = length(x)
  out = x[1:(len-(E*tau)+tau)]
  for(i in 2:E) { out = cbind(out,x[(1+((i-1)*tau)):(len-(E*tau)+(i*tau))]) }
  return(as.matrix(out))
}
# 2. Use column bind (cbind) the output of embedding each ID's data in turn
embedMultiPerson = function(x,E,tau,variableNames){
  for (j in 1:length(variableNames)){
    x_temp = x[,c("ID",variableNames[j])]
    colnames(x_temp)[2] = 'y'
    if (j==1){
      out = as.matrix(ddply(x_temp, .(ID), here(summarise), Embed(y,E,tau)))}
    else{
      out = cbind(out,as.matrix(ddply(x_temp, .(ID), here(summarise), Embed(y,E,tau))[,2]))}
  }
  return(out)
}

# # Edited fuction using "dplyr"
# head(ddply(long,.(ID),summarize, Embed(x1,E=6,tau=1)))
# long %>% group_by(ID) %>% summarise(Embed(x1,E=6,tau=1)) # %>% arrange(desc(total))
# # 1. Time delay embedding function
# Embed = function(x, E, tau) {  
#   len = length(x)
#   out = x[1:(len-(E*tau)+tau)]
#   for(i in 2:E) { out = cbind(out,x[(1+((i-1)*tau)):(len-(E*tau)+(i*tau))]) }
#   return(as.matrix(out))
# }
# # 2. Use column bind (cbind) the output of embedding each ID's data in turn
# embedMultiPerson = function(x,E,tau,variableNames){
#   for (j in 1:length(variableNames)){
#     x_temp = x[,c("ID",variableNames[j])]
#     colnames(x_temp)[2] = 'y'
#     if (j==1){
#       out = as.matrix(
#         ddply(x_temp, .(ID), here(summarise), Embed(y,E,tau))
#         )}
#     else{
#       out = cbind(out,as.matrix(ddply(x_temp, .(ID), here(summarise), Embed(y,E,tau))[,2]))}
#   }
#   return(out)
# }



# 3. Common model
DLOmodel <- function(t, prevState, parms) {
  x   <- prevState[1] # x[t]
  dx  <- prevState[2]
  
  with(as.list(parms), {
    dx  <- dx
    d2x <- parms[1]*x + parms[2]*dx
    res <-c(dx,d2x)
    list(res)
  })
}


# Data generation ----
set.seed(20220830)
Model <- "Osc" # Multilevel damped oscillator model
## Conditions ----
ny = 1        #Number of manifest indicators
neWithin = 1  #Number of within-person factors
neBetween = 1 #Number of between-person factors

## Edit Sample Configurations HERE ----
rep.start = 1
rep = 1#500
Samplesize <- 1#c(150, 300)   # 300 people only with 10 time-points
SeriesLength <- 99#c(9, 19, 99) # T = 10, 20, 100

Type='Raw'

# Cluster
n.cores <- parallel::detectCores() - 6
w.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(w.cluster)
# Add global variables to the clusters
parallel::clusterExport(w.cluster,
                        varlist = c("Model", "ny",'neWithin','neBetween','rep','Samplesize','SeriesLength','Type'))

setwd("C:/Users/YoungwonCho/OneDrive - The Pennsylvania State University/Simulation/MI/MIwithLoading")

setwd("C:/Users/YoungwonCho/OneDrive - The Pennsylvania State University/Simulation/MI/MNAR")

maxT=99
totalPersons=150
r=1

# ___________________________________________________________________________
runtime<-system.time(
  for(totalPersons in Samplesize){
    for(maxT in SeriesLength){
      if(totalPersons==300 & maxT!=9){next} # run 300-sample only when maxT==9
      #for(r in 1:10){
      foreach::foreach(r = rep.start:rep, #maxT = SeriesLength,
                       .packages=c("dplyr",'deSolve','Matrix','plyr','MASS','reshape2','mice')) %dopar% {
                         
                         ## begin ----                    
                         TimeSeq <- seq(0, maxT, by=1)  # the measurement occasions
                         # Note: because the time starts from 0, the # of time points we have is 'maxT+1' (not maxT)
                         deltaT <- unique(diff(TimeSeq))
                         
                         
                         ## Within-person dynamics ----
                         
                         # time-invariant covariate (e.g., MarSat, gender)
                         # Only 'cov1' will be used in the simulation
                         cov1 <- runif(totalPersons,-3,3)      #Continuous variable 
                         cov2 <- rbinom(totalPersons,1,.5)-0.5 #Binary variable (e.g., Gender): 0.5,-0.5
                         intcov<-runif(totalPersons,-2,2)
                         
                         # Empty Matrix
                         nn=1 # x
                         withinFactor <- matrix(NA, nrow=totalPersons, ncol=nn*(maxT+1)+4) 
                         # Because I'll assign [ID, Cov1, Cov2, intCov] into the first 4 columns, the # of columns is "4+2*(maxT+1)"
                         # one (maxT+1) for cov3, one (maxT+1) for y1
                         for (id in 1:totalPersons) {
                           # Assigning ID, Cov1, Cov2 into the first three columns
                           withinFactor[id,nn*(maxT+1)+1] <- id
                           withinFactor[id,nn*(maxT+1)+2] <- cov1[id]
                           withinFactor[id,nn*(maxT+1)+3] <- cov2[id]
                           withinFactor[id,nn*(maxT+1)+4] <- intcov[id]
                           
                           # Generating one PA items and assigning them into the 'withinFactor' matrix
                           # using lsoda(initial, times, function, parms) # parms: list of parameters used in function
                           if (Model=="Osc"){  
                             eta1 = -0.5 
                             zeta1= -0.04
                             
                             eta2 = -0.5
                             zeta2= -0.04
                             
                             parms <- c(eta1, zeta1, eta2, zeta2)
                             
                             #withinFactor[id,1:(maxT+1)+0*(maxT+1)] <- xcov3[1:(maxT+1)+(id-1)*(maxT+1)]
                             #withinFactor[id,1:(maxT+1)+1*(maxT+1)] <- ycov3[1:(maxT+1)+(id-1)*(maxT+1)]
                             
                             xstart <- c(x=rnorm(1,0,1), dx=rnorm(1,0,1)) # initial values
                             out1 <- as.data.frame(lsoda(xstart, TimeSeq, DLOmodel, parms))
                             withinFactor[id,1:(maxT+1)+0*(maxT+1)] <- out1$x #x: repeated measure/ID (wide)
                             #withinFactor[id,1:(maxT+1)+1*(maxT+1)] <- out1$y #y: repeated measure/ID (wide)
                           }}#End of loop through id
                         colnames(withinFactor) = c(paste0("x",0:maxT),
                                                    "ID","cov1","cov2",'intcov')
                         
                         ## Add between part ----
                         x2i = matrix(NA,totalPersons,ny)  #Empty cell for n-Intercept(Between-person)
                         x1it = matrix(NA,totalPersons,ny) #Empty cell for n-WithinFactors
                         #y2i = matrix(NA,totalPersons,ny)  #Empty cell for n-Intercept(Between-person)
                         #y1it = matrix(NA,totalPersons,ny) #Empty cell for n-WithinFactors
                         
                         xdataMatrix = matrix(NA,totalPersons*(maxT+1),ny) #Empty Matrix for the final result: y2i+y1it
                         #ydataMatrix = matrix(NA,totalPersons*(maxT+1),ny) 
                         
                         nu = matrix(c(0),ncol=1) # item-specific averages: arbitrary values. Here, it is just zero.
                         betweenFactorLoadings = matrix(c(1),ncol=neBetween)
                         
                         withinFactorLoadings = matrix(c(1),ncol=neWithin)
                         R = diag(c(1))
                         
                         # element for PA1
                         InterceptFactors = rep(0,totalPersons)    #no variance for intercept
                         
                         for (id in 1:totalPersons){
                           #Between-person effects  
                           #there can be a between person factor, summarizing systematic covariations across PA1, PA2, & PA3 across individuals
                           x2i[id,] = nu + betweenFactorLoadings %*% InterceptFactors[id]
                           #y2i[id,] = nu + betweenFactorLoadings %*% InterceptFactors[id]
                           # by multiplying by between_Factor_Loadings, we're creating the other two values(PA2, PA3).
                           # intercept_i = nu0 + nu1(=0.3)*covariate2[ID] + random_i(=rnorm(totalPersons,0,1))
                           for (t in 1:(maxT+1)){
                             x1it[id,] = withinFactorLoadings %*% withinFactor[id,t] + chol(R)%*%rnorm(ny,0,1)
                             #y1it[id,] = withinFactorLoadings %*% withinFactor[id,3*(maxT+1)+t] + chol(R)%*%rnorm(ny,0,1)
                             
                             # by multiplying by within_Factor_Loadings, we're creating the other two values...
                             xdataMatrix[t+(id-1)*(maxT+1),] = x2i[id,] + x1it[id,]
                             #ydataMatrix[t+(id-1)*(maxT+1),] = y2i[id,] + y1it[id,]
                           }
                           }
                         
                         
                         ## Long data ----
                         long = matrix(NA,totalPersons*(maxT+1),2+nn) #nn=4;c(x,xTVcov,y,yTVcov)
                         for (id in 1:totalPersons){ 
                           long[(1:(maxT+1))+(id-1)*(maxT+1),1] = id
                           long[(1:(maxT+1))+(id-1)*(maxT+1),2] = TimeSeq
                           
                           long[(1:(maxT+1))+(id-1)*(maxT+1),3] = xdataMatrix[(1:(maxT+1))+(id-1)*(maxT+1),] #x1
                     
                           #s<-summary(long[(1:(maxT+1))+(id-1)*(maxT+1),3:6])
                           #print(s[4,])
                         }
                         
                         # For the simulation, we are going to use the raw data before detrending, so we can model individual differences in intercepts.
                         colnames(long) = c("ID","Time", paste0("x",1:ny))
                         # long = data.frame(long, cov1=rep(withinFactor[,"cov1"],each=(maxT+1)),
                         #                   cov2=rep(withinFactor[,"cov2"],each=(maxT+1)),
                         #                   intcov=rep(withinFactor[,"intcov"],each=(maxT+1)))
                         head(long)
                        # write.table(round(long,5), file=paste0(Model,"Raw_long.txt"), append=F,
                        #             row.names=FALSE,col.names=TRUE,sep=",")
                         
                         ## Insert Missing (AMPUTE) ----
                         # ampute(
                         #   longdata,
                         #   prop = 0.5,
                         #   patterns = NULL,
                         #   freq = NULL,
                         #   mech = "MAR",
                         #   weights = NULL,
                         #   std = TRUE,
                         #   cont = TRUE,
                         #   type = NULL,
                         #   odds = NULL,
                         #   bycases = TRUE,
                         #   run = TRUE
                         # )
                         
                         # missing
                         mprob=c(.40, .20, .10, .05, 0) #missing ratio
                         for (p in mprob){
                         #amputed<- ampute(data.frame(long[,c('Time','x1')]), mech = "MCAR", prop=p)$amp[,2]
                         #newlong <- data.frame(long[,c('ID','Time')], x1=amputed)
                         
mypattern <- matrix(c(1,1,0), nrow = 1)
myweight <- matrix(c(0,0,1), nrow = 1)
newlong <- ampute(
   data = long,
   prop = p,
   patterns = mypattern,
   freq = NULL,
   mech = "MNAR",
   weights = myweight)$amp
                         #str(newlong)
                         #md.pattern(newlong)
                         
                         print(paste("proportion of missing is", round(1-nrow(na.omit(newlong))/nrow(newlong),2))) # proportion of missing
                         
                         newlong.comp <- na.omit(newlong)
                         #str(newlong.comp)
                         
                         # Add delta T ----
                         dat<-newlong.comp
                         #dat<-long
                        
                         dat$deltaT <- c(1, diff(dat$Time)) # End point was set to '1' here (replacing NA)
                         dat$deltaT <-ifelse(dat$deltaT<0, 1, dat$deltaT)
                         print("Delta t distribution is")
                         print(table(dat$deltaT))
                         
                
                         ## Time Delay Embedding ----
                         tau = 1     # The lag between subsequent columns in the embedded matrix
                         deltaT = 1  # The amount of time elapsed between subsequent observations
                         embedD = 6  # The number of columns in the time-delay embedded matrix
                         
                         ## DATA WITH NAs ##
                         OscData = data.frame(newlong)
                         lnames = c("x1","Time")

                         embedAll = embedMultiPerson(OscData[,c("ID",lnames)], embedD, tau, lnames)
                         colnames(embedAll) = c("ID", paste0("x1lag",(embedD-1):0),paste0("t",(embedD-1):0))
                         
                         #dim(embedAll); dim(OscData)
                         ifelse(dim(OscData)[1]-dim(embedAll)[1]==(embedD-1)*totalPersons, print('good'), print('stop'))
                         
                         embedAll = data.frame(embedAll)
                         # new loadings
                         first <- embedAll[,paste0("t",(embedD-1):0)]-rowMeans(embedAll[,paste0("t",(embedD-1):0)])
                         second <- (first)^2/factorial(2)
                         third  <- (first)^3/factorial(3)
                         fourth <- (first)^4/factorial(4)
                         
                         final<-data.frame(embedAll, first, second, third, fourth)
                         colnames(final)<- c("ID", paste0("x1lag",(embedD-1):0),
                         paste0("t",(embedD-1):0),
                         paste0("a",(embedD-1):0),paste0("b",(embedD-1):0))
                         
                         filename=paste0("UniOSCRaw_","w_NA","_TDE",embedD,"_n",totalPersons,'_t',maxT+1,'_','M',p*100,'_L',formatC(r, width=4, flag="0"),".txt")
                         write.table(round(final,5),
                                     file=paste0(getwd(),"/",filename),
                                     append=F, row.names=FALSE, col.names=FALSE, sep=",",
                                     na = "-999")
                         
                         
                         ## DATA WITHOUT NAs ##
                         OscData = data.frame(dat) #; Type='Raw'
                         lnames = c("x1","Time")
                         embedAll = embedMultiPerson(OscData[,c("ID",lnames)], embedD, tau, lnames)
                         colnames(embedAll) = c("ID", paste0("x1lag",(embedD-1):0),paste0("t",(embedD-1):0))
                         
                         embedAll = data.frame(embedAll)
                         # new loadings
                         first <- embedAll[,paste0("t",(embedD-1):0)]-rowMeans(embedAll[,paste0("t",(embedD-1):0)])
                         second<- (first)^2/factorial(2)
                         third <- (first)^3/factorial(3)
                         fourth<- (first)^4/factorial(4)
                         
                         final<-data.frame(embedAll, first, second, third, fourth)
                         colnames(final)<- c("ID", paste0("x1lag",(embedD-1):0),
                         paste0("t",(embedD-1):0),
                         paste0("a",(embedD-1):0),paste0("b",(embedD-1):0))
                         
                         # ##
                         # for(t in 1:nrow(dat)){
                         #  slope<- dat$x1[t+1]-dat$x1[t]/dat$deltaT[t]
                         #  print(paste(t,"th", round(slope,3)))
                         # }
                         
                         
                         
                         # ##
                         # deltaT<-data.frame(embedAll[,8:13])
                         # colnames(deltaT)<-c('dt1','dt2','dt3','dt4','dt5','dt6')
                         # 
                         # #dx loadings    
                         # deltaT$a1 = -2.5*deltaT$dt1;
                         # deltaT$a2 = -1.5*deltaT$dt2;
                         # deltaT$a3 = -0.5*deltaT$dt3;
                         # deltaT$a4 = 0.5* deltaT$dt4;
                         # deltaT$a5 = 1.5* deltaT$dt5;
                         # deltaT$a6 = 2.5* deltaT$dt6;
                         # 
                         # #d2x loadings
                         # deltaT$b1 = (-2.5*deltaT$dt1)**2/2;
                         # deltaT$b2 = (-1.5*deltaT$dt2)**2/2;
                         # deltaT$b3 = (-0.5*deltaT$dt3)**2/2;
                         # deltaT$b4 = (0.5*deltaT$dt4)**2/2;
                         # deltaT$b5 = (1.5*deltaT$dt5)**2/2;
                         # deltaT$b6 = (2.5*deltaT$dt6)**2/2;
                         # 
                         # final<-data.frame(embedAll[,1:7], deltaT)
                         # 
                         # psych::describe(final)

                         
                         #dataDir=paste0("C:/Users/",username,"/OneDrive - The Pennsylvania State University/Simulation/Data");setwd(dataDir)
                         #r=1
                         filename=paste0("UniOSCRaw_","wo_NA","_TDE",embedD,"_n",totalPersons,'_t',maxT+1,'_','M',p*100,'_L',formatC(r, width=4, flag="0"),".txt")
                         
                         write.table(round(final,5),
                                     file=paste0(getwd(),"/",filename),
                                     append=F, row.names=FALSE, col.names=FALSE, sep=",",
                                     na = "-999")
                       }
                         ## end ----
                       }}}
  #}
)

# Check runtime
lubridate::as.duration(runtime[3])

# Stop Cluster
parallel::stopCluster(w.cluster)

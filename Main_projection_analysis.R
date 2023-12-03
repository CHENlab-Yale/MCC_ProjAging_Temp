###################################################################################################################
## Sample R code for the MCC_ProjAging project
## "Impact of population aging on future temperature-related mortality at different global warming levels" 
## Ana Vicedo and Kai Chen
## August, 2023
###################################################################################################################

library(dlnm); library(lubridate); library(MASS); library(tsModel)

# DEFINE DIRECTORY
dir <- "D:/MCC/MCC_ProjAging/Analysis"

# LOAD the exposure_response functions estimated from a standard 2-stage time-series analysis with DLNM.
# load(paste(dir,"Output/MCCProjAging_2ndtstage19952014.RData",sep="/"))

# Load the age factors, computed elsewehre 
# load(paste(dir,"Output/Agefactor.RData",sep="/"))
# agefactor

####################################################################################
# LABELS DEFINING THE PROJECTION PERIODS
## Warming levels: 1.5 °C, 2.0 °C, and 3.0 °C global-mean annual-mean warming with respect to 1850 to 1900.
## Based on IPCC AR6 WG1 report: Cross-Section Box TS.1, Table 1: Assessment results for 20-year averaged change in global surface temperature
## based on multiple lines of evidence. The change is displayed in °C relative to the 13 1850–1900 reference period for selected time periods, 
## and as the first 20-year period during which the average global surface temperature change exceeds the specified level relative to
## the period 1850–1900. 
## Here we use the central estimate as the time period

period.ssp370 <- rbind(c(1995, 2014), c(2021, 2040), c(2037, 2056), c(2066, 2085))
period.ssp585 <- rbind(c(1995, 2014), c(2018, 2037), c(2032, 2051),c(2055, 2074))

period.list <- list(period.ssp370,period.ssp585)
####################################################################################

# CMIP6 CLIMATE SCENARIOS
rcp <- c("ssp370","ssp585")

# WARMING TARGETS
wtg <- c("historical", "1.5C", "2C", "3C")

# GENERAL CIRCULATION MODELS 
gcm <- c("awi_cm_1_1_mr", "miroc6", "miroc_es2l", "mri_esm2_0", "noresm2_mm", 
         "cnrm_cm6_1", "access_cm2", "bcc_csm2_mr", "cnrm_cm6_1_hr", "ipsl_cm6a_lr", 
         "cesm2","cnrm_esm2_1","gfdl_esm4", 
         "iitm_esm", "inm_cm5_0", "inm_cm4_8", "mpi_esm1_2_lr", "ukesm1_0_ll")

# NUMBER OF ITERATION IN THE MONTE-CARLO SIMULATION OF ATTRIBUTABLE RISK
nsim <- 1000 

# AGE CATEGORIES
agecat <- c("age064","age6574","age75plus")

# AGE NO AGE
age_noage <- c("age","noage")

# TEMPERATURE RAGE
cold_heat <- c("cold","heat")

# CREATE ARRAYS TO STORE INFO ON PROJECTED TEMPERATURE
tmeanperiod <- array(NA,c(nrow(cities),3,length(rcp),length(wtg),length(gcm)),
                     dimnames=list(cities$city,c("median","P25","P75"),names(rcp),wtg,gcm))

# CREATE A MATRIX TO STORE THE PROJECTED MORTALITY BY PERIOD
deathperiod <- array(0,dim=c(nrow(cities),length(rcp),length(wtg),length(agecat),2),
                     dimnames=list(cities$city,rcp,wtg,agecat,c("age","noage")))

# CREATE TEMPORARY OBJECT TO STORE SIMULATED RESULTS USED FOR UNCERTAINTY
# STRATIFIED ALSO BY GCM AND ITERATION, BUT NOT RCP
ancitysim <- array(0,dim=c(nrow(cities),length(rcp),length(wtg),length(gcm),2,length(agecat),2,
                           nsim+1),dimnames=list(cities$city,rcp,wtg,gcm,
                                                 c("cold","heat"),agecat,
                                                 c("age","noage"),
                                                 c("est",paste0("sim",seq(nsim)))))



################################################################################
# RUN THE LOOPS

################################
#### LOOP 1 - ACROSS CITIES ####
################################

for (i in c(1:nrow(cities))){
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dlist[[i]]
  
  # DEFINE THE OUTCOME
  data$y <- if(indnonext[i]) as.integer(data$nonext) else data[[out]]
  
  # DEFINE ARGVAR (AS IN ESTIMATION), CENTERING AND COEF-VCOV
  argvar <- list(fun=varfun,knots=quantile(data$tmean,varper/100,na.rm=T),
                 Bound=range(data$tmean,na.rm=T))
  if(!is.null(vardegree)) argvar$degree <- vardegree
  cen <- mintempcity[i]
  
  ################################
  ####   LOOP 2 - ACROSS RCP  ####
  ################################
  
  for (k in seq(length(rcp))){
    
    # LOAD PROJECTED TEMPERATURE SERIES FOR i CITY AND k RCP
    tmeanproj <- read.csv(paste0(dir,"Data/Tmean/",rcp[k],"/",names(dlist)[i],
                                 "_daily_Tmean_calproj.csv"))
    tmeanproj$X <- NULL
    tmeanproj$Date <- as.Date(tmeanproj$Date)
    ind1 <- substr(tmeanproj$Date,6,10)!="02-29" & tmeanproj$Date < "2086-01-01"
    tmeanproj <- tmeanproj[ind1,]
    
    # LOAD PROJECTED MORTALITY
    mortproj <- readRDS(paste0(dir,"Data/PopMort/",rcp[k],"/",names(dlist)[i],".RDS"))
    
    # DEFINE WEIGHTS FOR AGE GROUPS IN THE OBSERVED PERIOD
    deathobs <- if(indnonext[i]) as.integer(data$nonext) else 
      data[[out]]
    inddate <- with(data,date[date%in%mortproj$date & !is.na(deathobs)])
    wage <- colSums(mortproj[mortproj$date%in%inddate,c("u65","o65_74","o75")]) /
      sum(mortproj[mortproj$date%in%inddate,c("u65","o65_74","o75")])
    
    ####################################
    ####   LOOP 3 - ACROSS AGE/NOAGE  ##
    ####################################
    
    for (j in c(1,2)){  ## J==1  AGING // J=2 NO AGING 
      
      ################################
      ####   LOOP 4 - ACROSS WTG  ####
      ################################
      
      for (w in seq(length(wtg))){
        
        # DEFINE PERIOD BY SCENARIO
        selperiod <- period.list[[k]]
        selperiod <- seq(selperiod[w,1], selperiod[w,2])
        
        # SELECT TEMPERATURE SERIES
        indperiod <- substr(tmeanproj$Date,1,4)%in%selperiod
        tmeansel <- tmeanproj[indperiod,]
        
        if(j==1){
          # SELECT MORTALITY PROJECTIONS
          indperiod <- substr(mortproj$date,1,4)%in%selperiod
          mortprojsel <- mortproj[indperiod,c("u65","o65_74","o75")]
        } else {
          selperiod <- period.list[[k]]
          selperiod <- seq(selperiod[1,1], selperiod[1,2])
          indperiod <- substr(mortproj$date,1,4)%in%selperiod
          mortprojsel <- mortproj[indperiod,c("u65","o65_74","o75")]
        }
        
        deathperiod[i,k,w,,j] <- round(colSums(mortprojsel)/(nrow(mortprojsel)/365),0)
        
        ################################
        ####   LOOP 5 - ACROSS gcm  ####
        ################################
        
        for (m in seq(length(gcm))){
          
          # SELECT MODEL
          tmeanmodelsel <- tmeansel[,which(names(tmeansel)==gcm[m])]
          while(any(isna <- is.na(tmeanmodelsel)))
            tmeanmodelsel[isna] <- rowMeans(Lag(tmeanmodelsel,c(-20,20)),na.rm=T)[isna]
          # SUMMARY
          tmeanperiod[i,,k,w,m] <- quantile(tmeanmodelsel,c(0.5,0.25,0.75), na.rm=T)
          
          # DERIVE THE CENTERED BASIS AND EXTRACT PARAMETERS
          bvar <- do.call(onebasis,c(list(x=tmeanmodelsel), argvar))
          cenvec <- do.call(onebasis,c(list(x=cen),argvar))
          bvarcen <- scale(bvar,center=cenvec,scale=F)
          coef <- blupall[[i]]$blup
          vcov <- blupall[[i]]$vcov
          
          # INDICATOR FOR COLD/HEAT DAYS
          indheat <- tmeanmodelsel>cen
          
          # COMPUTE AGE-SPECIFIC LOG EXPOSURE-RESPONSE - USE WAGE (AGE-SPECIFIC WEIGHT BASED ON OBSERVED PERIOD)
          rr <- exp(drop(bvarcen%*%coef))
          rr064 <- (rr-1) / (wage[1] + wage[2]*agefactor[,1][indheat+1] + 
                               wage[3]*agefactor[,2][indheat+1]) + 1
          rr6574 <- agefactor[,1][indheat+1] * (rr064-1) + 1
          rr75plus <- agefactor[,2][indheat+1] * (rr064-1) + 1
          
          # ATTRIBUTABLE DEATHS 
          rrall <- cbind(rr064,rr6574,rr75plus)
          an <- (rrall-1)/rrall*as.matrix(mortprojsel)
          ancitysim[i,k,w,m,"cold",,j,1] <- colSums(an[!indheat,], na.rm=T)
          ancitysim[i,k,w,m,"heat",,j,1] <- colSums(an[indheat,], na.rm=T)
          
          # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
          set.seed(13041975)
          coefsim <- mvrnorm(nsim,coef,vcov)
          
          # LOOP ACROSS ITERATIONS
          for(s in seq(nsim)) {
            
            # COMPUTE AGE-SPECIFIC LOG EXPOSURE-RESPONSE
            rr <- exp(drop(bvarcen%*%coefsim[s,]))
            rr064 <- (rr-1) / (wage[1] + wage[2]*agefactor[,1][indheat+1] + 
                                 wage[3]*agefactor[,2][indheat+1]) + 1
            rr6574 <- agefactor[,1][indheat+1] * (rr064-1) + 1
            rr75plus <- agefactor[,2][indheat+1] * (rr064-1) + 1
            
            # ATTRIBUTABLE DEATHS
            rrall <- cbind(rr064,rr6574,rr75plus)
            an <- (rrall-1)/rrall*as.matrix(mortprojsel)
            ancitysim[i,k,w,m,"cold",,j,1+s] <- colSums(an[!indheat,], na.rm=T)
            ancitysim[i,k,w,m,"heat",,j,1+s] <- colSums(an[indheat,], na.rm=T)
            
          }
        }
      }  
    }        
  }
}            

##Save the results 
# save.image(paste(dir,"Output/03projMCCAging.RData",sep="/"))

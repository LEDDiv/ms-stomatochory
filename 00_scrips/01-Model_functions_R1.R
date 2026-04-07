#==============================================================================
## Function of Seed Dispersal by Stomatochory
#==============================================================================

rm(list=ls())
getwd()

## Loading packages ----
#install.packages("renv")
#renv::init()
library(devtools)
library(adehabitatHR)
library(rio)
#install_github('MarieAugerMethe/CCRWvsLW') # run with R 4.0 version
#You will probably need to install an older version of R to
#run this package. It works with version 4.0.5
library(CCRWvsLW)
library(VGAMdata)

# Parameters for competition and stomatochory probability functions:
# P_base <- baseline prob. of stomatochory affected by factors other than competition
# P_max <- upper limit for the probability of stomatochory.

## Function ----
Macaw_Kernel <- function(ind= c(100, 3000), rep= 100, n= (0.6*3000), b= 3, P_base= 0.1, P_max= 0.7,
                        mean.fru.car=11, p.fru.car.pred=0.74, p.landing= 0.6, a=3, mu=2.3,
                        steps= 600, k=20,walk= 10, movement="levy"){
  #-----------------------------------------------------------
  #1. Stomatochory Events or Indiv doing Stomatochory
  #-----------------------------------------------------------
  S <- length(ind) # Scenarios
  
  Esto.ev <- list() # n events
  Competition_n <- list() # intensity of competition in each scenario
  Probability_n <- list() # Prob of stomatochory affected by competition in each scenario
  for(s in 1:S){
    C_n <- (ind[s]^b) / ((n^b) + (ind[s]^b))
    p.esto <- P_base + (P_max - P_base) * C_n
    n_esto <- rbinom(n= rep,size= ind[s],prob= p.esto)
    Esto.ev[[s]]<- n_esto
    Competition_n[s] <- C_n
    Probability_n[s] <- p.esto
  }
  
  #------------------------------------------------------------
  #2. Fruits Carried per Stomato Events 
  #------------------------------------------------------------
  
  All.fru.car.ev <- list()
  AUX<- list()
  for(s in 1:S){
    for(r in 1:length(Esto.ev[[s]])){
      fru_car_ev <- VGAMdata::rpospois(n= Esto.ev[[s]][r], lambda= mean.fru.car)
      AUX[[r]] <- fru_car_ev
    }
    All.fru.car.ev[[s]]<-AUX
  }
  
  #each element of the list is a replicate,
  #each element of the vector is the number of seeds per estomat event
  
  #-------------------------------------------------------------
  #3. Seeds dispersed per Stomato Events 
  #-------------------------------------------------------------
  
  All.fru.disp.ev <- list()
  AUX <- list()
  AUX2 <- list()
  deposit.time.all <- list(NA)
  for(s in 1:S){
    for (r in 1:length(All.fru.car.ev[[s]])){
      fru_disp_ev <- rbinom(n = length(All.fru.car.ev[[s]][[r]]),
                            size= All.fru.car.ev[[s]][[r]],
                            prob= 1-p.fru.car.pred) #prob of seed predation
      AUX[[r]] <- fru_disp_ev
      All.fru.disp.ev[[s]] <- AUX
      
      
      # replace 0 by NA
      replace <- function(vector) {
        ifelse(vector == 0, NA, vector)
      }
      
      # apllying fc
      All.fru.disp.ev.na <- lapply(All.fru.disp.ev, function(inside_list) {
        lapply(inside_list, replace)
      })
      
      # removing NA
      remove_nna <- function(vector) {
        vector[!is.na(vector)]
      }
      
      # applying to each vector of the sublist of lists
      All.fru.disp.ev.sna <- lapply(All.fru.disp.ev.na, function(sublist) {
        lapply(sublist, remove_nna)
      })
      
      #---------------------------------------------------------
      #4. The moment the fruits are dispersed (number of steps until stop)
      #---------------------------------------------------------
      # considering that the weight of the propagules has no influence on the
      # flight distance
      
      dt <- matrix(NA,1,length(All.fru.disp.ev.sna[[s]][[r]]))
      for (ev in 1:length(All.fru.disp.ev.sna[[s]][[r]])){
        dt[0:1,ev] <- rgeom(1,p.landing)
      }
      AUX2[[r]] <- round(dt)
    }
    deposit.time.all[[s]] <- AUX2
  }
  
  # Add one value to each vector element
  add_one <- function(vector) {
    vector + 1
  }
  
  # Applying
  deposit.time.all.sum <- lapply(deposit.time.all, function(sub_list) {
    lapply(sub_list, add_one)
  })
  
  
  #each element of the list is a replicate
  #in each matrix, the column represents a stomatochory event,
  #the row represents each fruit
  #the matrix element represents the time step a given seed has been deposited
  
  #------------------------------------------------------
  #5. Simulating Movement
  #------------------------------------------------------
  dist.all <- list()  #distances for individuals of all spp
  for(s in 1:S){
    ev_disp <- unlist(lapply(deposit.time.all.sum[[s]], ncol))
    for (r in 1:length(ev_disp)){ #for each replicate
      
      st.dist <- matrix(NA,steps+2,ev_disp[r]) #if others matrix has n+2 rows
      
      
      step.size <- 1/walk #step size for species s
      
      # one movement simulation per stomatochory event
      
      for (i in 1:ev_disp[r]){ #Simulating individual walks
        if (movement=="brown"){
          w <- simmBW(n=steps, l=step.size, a=a)
        }else if(movement=="levy"){
          w<- simmLW(n=steps, mu=mu, a=a)
        }else if(movement=="corr"){
          w <- simmCRW(n=steps,l=step.size,k=k,a=a)
        }
        mov<-w[[1]]
        st.dist[,i] <- sqrt(mov[,1]^2+mov[,2]^2)
      }
      AUX[[r]]<- st.dist
    }
    dist.all[[s]]<-AUX #Distance from source at each time step
  }
  
  #---------------------------------------------------
  #6. Combining deposit time and movement
  #---------------------------------------------------
  # replaces the time each seed was deposited by the distance traveled by the
  # individual at that time step
  
  dist.time.all <- list(NA)
  for(s in 1:S){
    ev_disp <- unlist(lapply(deposit.time.all.sum[[s]], ncol))
    
    for(r in 1:length(deposit.time.all.sum[[s]])){
      
      deposit.time <- deposit.time.all.sum[[s]][[r]]#matrix representing  each replica
      # Calls each of the arrays in this object
      dist.time <- matrix(NA,nrow(deposit.time),ev_disp[r])
      dist <- dist.all[[s]][[r]] #distance matrix at each step of each replica
      deposit.time[deposit.time>steps] <- steps
      for(i in 1:ev_disp[r]){
        dist.time[,i] <- dist[deposit.time[,i]+1,i]
      }
      dist.time <- rep(dist.time, All.fru.disp.ev.sna[[s]][[r]]) #repeats the
      #deposition distance for each fruit in the bunch
      AUX2[[r]] <- dist.time
    }
    dist.time.all[[s]] <- AUX2 #Distance each seed was deposited for each event
  }
  return(dist.time.all)
}




#==============================================================================
## Simulations
#==============================================================================

## movement models ----
## Simulating Population Scenarios and parameters ---

ind <- c(100, 500, 1000, 2000, 3000)
S <-length(ind) # different scenarios
replicas <-10
steps <- 600 #,3600)
walk <- 10

#set.seed(42)
LW_sim <- Macaw_Kernel(ind, rep= replicas, steps= steps,walk= walk,
                       movement="levy")


BW_sim <- Macaw_Kernel(ind, rep= replicas,steps= steps,walk= walk,
                       movement="brown")


CRW_sim <- Macaw_Kernel(ind, rep= replicas,steps= steps,walk= walk,
                        movement="corr")


#Saving simulations in .rds (read by R)
saveRDS(LW_sim, "../output/LW_sim.rds", compress = FALSE)
saveRDS(BW_sim, "../output/BW_sim.rds", compress = FALSE)
saveRDS(CRW_sim, "../output/CRW_sim.rds", compress = FALSE)

## Sensitivity analyses ----

## Landing-flight probability
# variation in landing probability as a proxy for perch density around feeding areas
# simulations for the levy walks movement pattern

# p.landing=0.2
LW_sim_p2 <- Macaw_Kernel(ind, rep= replicas, p.landing= 0.2,
                          steps= steps, walk= walk,movement="levy")


# p.landing=0.4
LW_sim_p4 <- Macaw_Kernel(ind, rep= replicas, p.landing= 0.4,
                          steps= steps, walk= walk,movement="levy")

# p.landing=0.8
LW_sim_p8 <- Macaw_Kernel(ind, rep= replicas,p.landing= 0.8,
                          steps= steps, walk= walk,movement="levy")



#Saving
saveRDS(LW_sim_p2, "../output/LW_sim_p2.rds", compress = FALSE)
saveRDS(LW_sim_p4, "../output/LW_sim_p4.rds", compress = FALSE)
saveRDS(LW_sim_p8, "../output/LW_sim_p8.rds", compress = FALSE)


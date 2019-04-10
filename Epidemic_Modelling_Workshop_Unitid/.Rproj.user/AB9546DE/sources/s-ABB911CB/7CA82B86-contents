# Pulliam_Silkey2016.R
# Day 2: Models for study design
# SACEMA-UNITID Epidemiological Modelling Workshop
# University of Narobi, Nairobi, Kenya, 8-11 April 2019
#
# Juliet R.C. Pulliam, 2019
#
# MODIFIED FROM: 
#
## SWCRT_ARMA
## R code for simulating an SWCRT for an intervention against an infectious disease
## Program TransmissionSimulator.R (copyright Mariabeth Silkey and Thomas Smith, 2013)
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU GPL version 3 as published by the Free Software Foundation.
## <http://www.gnu.org/licenses/
##
## Original code was provided as Additional File 2 as part of this publication:
## Silkey et al. 2016 Trials 17:278
## https://trialsjournal.biomedcentral.com/articles/10.1186/s13063-016-1378-1
##
## OVERVIEW
##
## first for loop populates the area with individuals and randomly assigning initial infection status
## second for loop introduces a pure ARMA (autoregressive moving average) process under which steady state of infection is reached;
## once the intervention is introduced a more complex ARMA process takes over, including direct and indirect effects of the intervention
## cluster sequences for each design type are generated in separate code and stored in external databases "cluster_order.csv"
## datasets distmatrix and hh created using function spDists from package sp to calculate the pairwise distances in km across all households
## first line loads dataset distmatrix and data set hh which contains demographic information on each household
rm(list = ls())
load('rusingahh.RData')

# Define missing parameters -----------------------------------------------
initialIncidence <- 0.5
clseq <- 296
radius <- 1#00
suffix <- 'ACG'
initialPrevalence <- 0.22
efficacy <- 0.3
#===============================================================
# set constants 
endInitial  = 10; endOfBurnIn   = 40  # simple ARMA process endpoint
maxTimeSteps= 160; nIterations   = 4
nhh         = nrow(distmatrix)      # nhh is the number of households
underlying  = -log(1-initialIncidence)/initialIncidence
#  underlying rate constrained so that the ARMA process will stabilize and continuously 
#  match values  from connecting sections (run-in, rollout, post-intervention)

#===============================================================
# get cluster sequence    
cluster_o = read.csv('cluster_order.csv', header=F, nrow=1, col.names='clseq', skip=clseq-1)[,1]
# calculate the number of and identify neighboring structures within radius 
neighborHS = apply(distmatrix,1, function(X) which(X<radius)) # neighboring houses
# initialize data structure for results
write.table('clseq, radius, efficacy, initialPrevalence, iteration, t, treatCC, compare1CC
  , compare2CC, compare3CC, intervenedno, comparator1no, comparator2no, comparator3no'
            , file=paste("out", suffix, sep='_'), quote=F, row.names=FALSE, col.names=FALSE)

#===============================================================
# first for loop populates the area with individuals and randomly assigning initial infection status
for (iteration in 1:nIterations){   
  # clear containers at the beginning of each iteration
  CCMatrix  = matrix(data=0, nrow = maxTimeSteps, ncol=nhh)   # clinical case status   
  coverageMatrix = matrix(data=0, nrow = maxTimeSteps, ncol=nhh) # average intervention status  
  aveCCnhood = matrix(data=0, nrow = maxTimeSteps, ncol=nhh)   # average of clinical cases  
  reservoir  = matrix(data=0, nrow = maxTimeSteps, ncol=nhh)   # weighted mean of coverage  
  boundaryNeighbors  = matrix(data=0, nrow = maxTimeSteps, ncol=nhh) # households on the border 
  intervened = rep(0, nhh)                                   # households with intervention
  output = matrix(data=0, nrow = maxTimeSteps, ncol = 10 )   # to hold clinical case counts
  
  # second for loop introduces a pure ARMA (autoregressive moving average) process under which steady state of infection is reached;
  for(  t in 1:endInitial){ 
    # assignment of original disease state. CC is clinical cases per household
    CCMatrix[t,] = rbinom(nhh, hh$N, initialPrevalence)
    aveCCnhood[t,] = sapply(neighborHS, function(X) {sum(CCMatrix[t,unlist(X)])/sum(hh$N[unlist(X)])}) # generates error
    output[t,] =  c(iteration, t, 0, sum(CCMatrix[t,]), sum(CCMatrix[t,]), 0, 0, sum(hh$N), sum(hh$N), 0 )      
  }   
  
  # once the intervention is introduced a more complex ARMA process takes over, including direct and indirect effects of the intervention
  for (t in (endInitial+1):maxTimeSteps) {
    reservoir[t,] =  
      0.1* aveCCnhood[t-10, ] + 0.1* aveCCnhood[t-8 , ] +
      0.2* aveCCnhood[t-7 , ] + 0.2* aveCCnhood[t-9 , ] + 0.4* aveCCnhood[t-8 , ] 
    
    # length of intervened is the number of households
    intervened        =   cluster_o[hh$cluster] <= (t - endOfBurnIn)
    coverageMatrix[t,]= sapply(neighborHS, function(X) weighted.mean(intervened[X], w=hh$N[X],na.rm=T))
    
    # infection status update - lovely strange syntax of rbinom 
    CCMatrix [t,]  = rbinom(rep(1,nhh), hh$N, (1-exp(-(1-coverageMatrix[t,]*efficacy)*reservoir[t,]*underlying))) 
    aveCCnhood[t,] = sapply(neighborHS, function(X) {sum(CCMatrix[t,unlist(X)])/sum(hh$N[unlist(X)])})     
    
    # comparator1 is all households who have not received the intervention
    comparator1    = 1 - intervened   # households of mixed comparator group
    
    # neighboring households of the intervened         
    neighborhh    = unique(unlist(neighborHS[intervened]))
    boundary      = setdiff(neighborhh, (1:nhh)[intervened])
    boundaryNeighbors[t,] = rep(0,nhh) # initialize vector
    boundaryNeighbors[t,boundary] = 1
    
    # comparator2 is the subset of comparator 1x that are far from the intervention  
    comparator2   = comparator1 - boundaryNeighbors[t,] 
    
    # output components for effectiveness measures		
    treatCC       = sum(CCMatrix [t, ]*intervened) 
    compare1CC    = sum(CCMatrix [t, ]*comparator1 )
    compare2CC    = sum(CCMatrix [t, ]*comparator2 ) 
    compare3CC    = sum(CCMatrix [t, ]*boundaryNeighbors[t,])
    intervenedno  = sum(intervened*hh$N)
    comparator1no = sum(comparator1*hh$N)
    comparator2no = sum(comparator2*hh$N)
    
    # comparator3 is the complement of comparator 2 (i.e., the set of non-intervened neighbors)
    comparator3no = sum(boundaryNeighbors[t,]*hh$N) 
    output[t,] =  c(iteration, t, treatCC, compare1CC, compare2CC,compare3CC, intervenedno, comparator1no, comparator2no, comparator3no)    
  }  
  output = cbind(clseq, radius, efficacy, initialIncidence, output)
  write.table(output, file=paste("out", suffix, sep='_'), sep = ",", append=TRUE,row.names=FALSE, col.names=FALSE)
  gc() # garbage collection to release as much memory as possible before returning to top of loop
}
#=============================================================== 
#   End of program
#===============================================================
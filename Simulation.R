rm(list = ls())

### set current directory


library(deSolve)
library(DEoptim)
library(pracma)
library(Hmisc)
library(scales)
library(RColorBrewer)
#library(ComplexHeatmap)
library(lattice)      
library(latticeExtra) 
library(grid)
library(circlize)
#library(scater)
library(ggpubr)
library(rapport)
library(dplyr)
library(gmodels)
library(ggplot2)
library(raster)
library(dplyr)
library(EnvStats)
library(deSolve)

################################
# Model Function
################################

SARS_COV_model <- function(t,x,params){
  with(as.list(x),{   
    
    ddt_S = -beta*V*S
    ddt_I = beta*V*S - delta*I^k*I - m*E^r*I/(E^r+E50^r) 
    ddt_V = p*I-c*V
    
    ddt_M1 = w*I*M1 - q*M1
    ddt_M2 = q*(M1-M2)
    ddt_E = q*M2- de*E
    
    der <- c(ddt_S,ddt_I,ddt_V,ddt_E,ddt_M1,ddt_M2)
    
    list(der)
  })       
}


################################   
# AUC calculation function
################################  
#AUC calculation function
AUC <- function(t,y) {
  y1=na.omit(y)
  mn=length(y1)
  val_AUC=(trapz(t,y))
  return(val_AUC)
}
################################
# Main
################################




################################
# 1. Read Population parameters
################################


Parameters = read.csv("SimulatedParameters.txt",header=TRUE)

Three_params = data.frame(read.csv(file="alpha_beta_contact_combos_p4_p2.csv",as.is=TRUE))

three_params_options_length=length(Three_params$alpha)

########################################################
########################################################

max_replicates=1000 ### how many individuals simulated

dT=1/6 # this has to be less than 1 --- every 4 hours
no_parameter_set <- c(1)
dispersion_options <- c(40/dT)

############################################ NEW
perc_people_wearing_mask= 0.0 # 0 mean nobody wearing mask , 1 means all wearing mask
mask_efficacy=0.0  ## 0 mean mask is useless, 1 means mask is great
perc_time_wearing_mask=0.0  ## 0 means no mask at all, 1 means all the time
############################################ END NEW

########## Specify inputs for the following

perc_infector_wearing_mask_options=c(seq(0,1,by=0.2))
perc_infectee_wearing_mask_options=c(seq(0,1,by=0.2))
perc_time_infector_wearing_mask_options=c(seq(0,1,by=0.2))
perc_time_infectee_wearing_mask_options=c(seq(0,1,by=0.2))
mask_infector_efficacy_options=c(seq(0,1,by=0.2))
mask_infectee_efficacy_options=c(seq(0,1,by=0.2))


#############################

PAR_proj_R0_SI_mask = matrix(0,nrow=length(perc_infector_wearing_mask_options)*length(perc_infectee_wearing_mask_options)*length(perc_time_infector_wearing_mask_options)*length(perc_time_infectee_wearing_mask_options)*length(mask_infector_efficacy_options)*length(mask_infectee_efficacy_options),ncol=17)
colnames(PAR_proj_R0_SI_mask)=c("perc_infector_wearing_mask","perc_infectee_wearing_mask", "perc_time_infector_wearing_mask", 
                                "perc_time_infectee_wearing_mask ",  "mask_infector_efficacy_options","mask_infectee_efficacy_options" ,
                                "R0","SI","GT","perc_transmitters","perc_SSE_5","perc_SSE_10","perc_SSE_20",
                                "trans_exposed_mask","trans_mask","exposed_mask","none_mask")

### inputs: perc_infector_wearing_mask, perc_infectee_wearing_mask , perc_time_infector_wearing_mask , perc_time_infectee_wearing_mask, mask_infector_efficacy_options,  mask_infectee_efficacy_options
### outputs: R0 (mean R0) ; SI ( Serial Interval) ; GT (Generation Time) ; perc_SSE_5 (percentage superspreaders of all infected individuals, where superspreaders are those infecting more than 5)
### outputs: perc_SSE_10 (percentage superspreaders of all infected individuals, where superspreaders are those infecting more than 10) ;  perc_SSE_20 (percentage superspreaders of all infected individuals, where superspreaders are those infecting more than 20)
### outputs: trans_exposed_mask (transmission from masked to masked individuals) ; trans_mask (transmission from masked to unmasked individuals) 
### outputs: exposed_mask (transmission from unmasked to masked individuals) ; none_mask (transmission from unmasked to unmasked individuals) 


par_ind_mask=1

for (perc_infector_wearing_mask in perc_infector_wearing_mask_options){
  for (perc_infectee_wearing_mask in perc_infectee_wearing_mask_options){
    for (perc_time_infector_wearing_mask in perc_time_infector_wearing_mask_options){
      for (perc_time_infectee_wearing_mask in perc_time_infectee_wearing_mask_options){
        for (mask_infector_efficacy in mask_infector_efficacy_options){ 
          for (mask_infectee_efficacy in mask_infectee_efficacy_options){
            

            print(length(perc_infector_wearing_mask_options)*length(perc_infectee_wearing_mask_options)*length(perc_time_infector_wearing_mask_options)*length(perc_time_infectee_wearing_mask_options)*length(mask_infector_efficacy_options)*length(mask_infectee_efficacy_options))
            print(par_ind_mask)
            
            
            mean_R0_iter <- vector()
            mean_SI_iter<-vector()
            mean_GT_iter<-vector()
            mean_Trans_iter <- vector()
            mean_SSE5_of_all_iter<-vector()
            mean_SSE10_of_all_iter<-vector()
            mean_SSE20_of_all_iter<-vector()
            trans_exposed_mask<-vector()
            trans_mask<- vector()
            exposed_mask<-vector()
            none_mask<-vector()
            
            no_iterations <- c(seq(1,3,by=1))  ## 3 replicates to get mean R0 (to overcome stochasticity)
            for (iter in no_iterations){
  
              ijk=1
              
              R0_all_infections = matrix(0,nrow=length(dispersion_options)*max_replicates*150,ncol=4) # 15 is for a maximum 15 time points
              colnames(R0_all_infections)=c("ID","time","dispersionx","new_infections")
            
              par_ind=1
              
              par_ID=1
              
              
              R0_more_than_10<-vector()
              
              lambda_options_superspreaders<-vector()
              lambda_options_non_superspreaders<-vector()
              
              contact_options_superspreaders<-vector()
              contact_options_non_superspreaders<-vector()
              
              alpha_options_superspreaders<-vector()
              alpha_options_non_superspreaders<-vector()
              
              AUC_risk_superspreaders<-vector()
              AUC_risk_non_superspreaders<-vector()
              
              PeakVL_superspreaders<-vector()
              PeakVL_non_superspreaders<-vector()
     
              for (dispersion in seq(1,length(dispersion_options),by=1)) { ### we only choose one value but the user can specify more and store this as well
             
                first_loop=seq(1,max_replicates,by=1) # we simulate 1000 transmitters but the user could simulate 10000 transmitters 
                ijx=1
                
                R0_temp_all_infections = matrix(0,nrow=max_replicates*31/dT,ncol=9) # 15 is for a maximum 15 time points
                colnames(R0_temp_all_infections)=c("ID","time","VL","contact","lambda","alpha","R0","R0_cumul","Pre")
                
                itx=1
                
                R0_infections<- vector()
                SI_infections<- vector()
                GT_infections<- vector()
                
                contact_options<-vector()
                lambda_options<-vector()
                alpha_options<-vector()
                AUC_risk <- vector()
                
                R0_asym_infections<- vector()
                
                asymptomatic_incidence<-vector()
                
                presymptomatic_incidence<- vector()
                symptomatic_incidence<- vector()
                
                contact_options<-vector()
                lambda_options<-vector()
                alpha_options<-vector()
                AUC_risk <- vector()
                peakVL_options<-vector()
                
                
                for (yu in first_loop){
  
                  tmin=min(Parameters$tzero[c(seq(1,max_replicates,by=1))])
                  
                  No_to_select=1 ### We select lambda, alpha, theta from the file
                   
                  three_mtp_options<- c(no_parameter_set) #5th parameter set
                  tzero=Three_params$tzero[three_mtp_options[1]]
                  AUC_factor=Three_params$betat[three_mtp_options[1]] 
                  alpha=Three_params$alpha[three_mtp_options[1]] 
                  no_contact_per_day=Three_params$no_contact[three_mtp_options[1]] 
                  
                  mtp=raster::sampleInt(10^4, 1, replace = FALSE) ### We select viral kinetic parameters from the file
                  beta=Parameters$beta[mtp]
                  delta=Parameters$delta[mtp]
                  k=Parameters$k[mtp]
                  p=Parameters$p[mtp]
                  m=Parameters$m[mtp]
                  w=Parameters$w[mtp]
                  E50=Parameters$E50[mtp]
                  r=Parameters$r[mtp]
                  q=Parameters$q[mtp]
                  de=Parameters$de[mtp]
                  c=Parameters$c[mtp]
                  S_0 = 1e7
                  I_0 = 1
                  V_0 = p*I_0/c
                  E_0 = 0
                  M1_0 = 1
                  M2_0 = 0
                  
                  
                  
                  ## the incubation period (time from exosure to symtptoms) are gamma distributed -- an incubation period with a mean of 6.4 and SD of 2.3 days [4], and a mean of 4.8 and a SD 2.6 days [7]  --- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/
                  mean_incubation_period=5.2 # mean from https://www.acpjournals.org/doi/10.7326/M20-0504 and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/
                  std_incubation_period=2.8 # days --  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/
                  incubation_period_infector=pmax(0,rgamma(n = 1, shape = 3.45,
                                                           scale=(1/0.66))) # rate=0.66, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/
                  
                  
                  Non_replication_phase_duration_infector=tzero
                  
                  #### For transmitter
                  tzero2=tzero
                  
                  incubation_period_infectee=pmax(0,rgamma(n = 1, shape = 3.45,
                                                           scale=(1/0.66))) # rate=0.66, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/
                  
               
                  Non_replication_phase_duration_infectee=tzero2
                  
                  ########################################################################
                  ################## Checking if transmission occur in asym phase or not with viral loads of ~1 copies/mL
                  ########################################################################
                  
                  Non_replication_phase_duration_infector_round=round(Non_replication_phase_duration_infector)
                  
                  Non_replication_phase_duration_infector_round_discrete=c(seq(0,Non_replication_phase_duration_infector_round,by=dT))
                  V_asym=rep((V_0),length(Non_replication_phase_duration_infector_round_discrete))
                  
                  
                  ## check for successful transmissions (since transmission risk is not successful transmission)
                  no_contact_per_day_options<-c()
                  no_contact_per_day_options=rgamma(n = length(Non_replication_phase_duration_infector_round_discrete), shape = (no_contact_per_day/dispersion_options[dispersion]),
                                                    scale=dispersion_options[dispersion]) # Random variable X is distributed X=gamma(A,B) with mean =AB 
                  
                  
                  ############################ New Infector and infectee wearing mask
                  is_infector_wearing_mask=rbinom(n=1,size=1,prob=(perc_infector_wearing_mask)) # 1 wearing mask --- we decide this once
                  
                  V_asym1 = V_asym ## when only trnamitter wearing mask
                  V_asym2=V_asym ## when both wearing mask
                  V_asym3=V_asym ## when only exposed contact wearing mask
                  V_asym4=V_asym ## when none is wearing mask
                  
                  perc_infectee_wearing_mask_at_all_times=perc_infector_wearing_mask*perc_time_infectee_wearing_mask
                  if(perc_infectee_wearing_mask_at_all_times>0){
                    frac_infectee_wearing_mask_average = (1-(perc_infectee_wearing_mask_at_all_times))
                    frac_infectee_wearing_mask_at_sampled_times = rbinom(n=length(no_contact_per_day_options),size=1,prob=(frac_infectee_wearing_mask_average))
                    number_infectee_wearing_mask_at_sampled_times =  no_contact_per_day_options * frac_infectee_wearing_mask_at_sampled_times 
                    number_infectee_not_wearing_mask_at_sampled_times =  no_contact_per_day_options * (1-frac_infectee_wearing_mask_at_sampled_times)
                  } else{
                    number_infectee_wearing_mask_at_sampled_times =  no_contact_per_day_options * 0 
                    number_infectee_not_wearing_mask_at_sampled_times =  no_contact_per_day_options * 1
                    
                  }
                  
                  random_numbers_asym1= runif(n=length(V_asym1),min=0,max=1)  # uniform distribution on the interval from min to max.
                  random_numbers_asym2= runif(n=length(V_asym2),min=0,max=1)  # uniform distribution on the interval from min to max.
                  random_numbers_asym3= runif(n=length(V_asym3),min=0,max=1)  # uniform distribution on the interval from min to max.
                  random_numbers_asym4= runif(n=length(V_asym4),min=0,max=1)  # uniform distribution on the interval from min to max.
                  
                  
                  if(is_infector_wearing_mask==1){
                    
                    infector_wearing_mask_when= rbinom(n=length(Non_replication_phase_duration_infector_round_discrete),size=1,prob=(perc_time_infector_wearing_mask)) # 1 means wearing mask and will not transmit at that time
                    
                    ## times when only transmitter wearing mask
                    V_asym1[which(infector_wearing_mask_when==1)] = (1-mask_infector_efficacy)*
                      V_asym1[which(infector_wearing_mask_when==1)]  ## Viral load for infectiousness/contagiousness at the time of wearing mask should be reduced with mask efficacy
                    V_asym1[which(infector_wearing_mask_when==0)] = 0
                    V_asym1[is.na(V_asym1)]=0
                    
                    Prob_V_asym1=(pmax(0,(V_asym1)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_asym1)))^alpha)
                    
                    ## times when both wearing mask
                    V_asym2[which(infector_wearing_mask_when==1)] = (1-mask_infectee_efficacy)*
                      (1-mask_infector_efficacy)*
                      V_asym2[which(infector_wearing_mask_when==1)]
                    V_asym2[which(infector_wearing_mask_when==0)]=0
                    V_asym2[is.na(V_asym2)]=0
                    
                    Prob_V_asym2=(pmax(0,(V_asym2)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_asym2)))^alpha)
                    
                    ## when only exposed contact wearing mask and transmisio is from tramitter at time when not wearing mask
                    V_asym3[which(infector_wearing_mask_when==0)] = (1-mask_infectee_efficacy)*
                      V_asym3[which(infector_wearing_mask_when==0)]
                    V_asym3[which(infector_wearing_mask_when==1)]=0
                    V_asym3[is.na(V_asym3)]=0
                    
                    Prob_V_asym3=(pmax(0,(V_asym3)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_asym3)))^alpha)
                    
                    ## when none is wearing mask and transmisio is from tramitter at time when not wearing mask
                    V_asym4[which(infector_wearing_mask_when==0)] = V_asym4[which(infector_wearing_mask_when==0)]
                    V_asym4[which(infector_wearing_mask_when==1)]=0
                    V_asym4[is.na(V_asym4)]=0
                    
                    Prob_V_asym4=(pmax(0,(V_asym4)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_asym4)))^alpha)
                    
                    # From masked transmitter at times when wearing mask to unmasked contacts
                    number_of_new_infections_actual_asym1=rep(0,length=length(Non_replication_phase_duration_infector_round_discrete))
                    number_of_new_infections_actual_asym1[Prob_V_asym1>random_numbers_asym1] = dT*number_infectee_not_wearing_mask_at_sampled_times[Prob_V_asym1>random_numbers_asym1]*Prob_V_asym1[Prob_V_asym1>random_numbers_asym1]
                    # From masked transmitter at times when wearing mask  to masked contacts
                    number_of_new_infections_actual_asym2=rep(0,length=length(Non_replication_phase_duration_infector_round_discrete))
                    number_of_new_infections_actual_asym2[Prob_V_asym2>random_numbers_asym2] = dT*number_infectee_wearing_mask_at_sampled_times[Prob_V_asym2>random_numbers_asym2]*Prob_V_asym2[Prob_V_asym2>random_numbers_asym2]
                    # From masked transmitter at times when not wearing mask to masked contacts
                    number_of_new_infections_actual_asym3=rep(0,length=length(Non_replication_phase_duration_infector_round_discrete))
                    number_of_new_infections_actual_asym3[Prob_V_asym3>random_numbers_asym3] = dT*number_infectee_wearing_mask_at_sampled_times[Prob_V_asym3>random_numbers_asym3]*Prob_V_asym3[Prob_V_asym3>random_numbers_asym3]
                    # From masked transmitter at times when not wearing mask  to unmasked contacts
                    number_of_new_infections_actual_asym4=rep(0,length=length(Non_replication_phase_duration_infector_round_discrete))
                    number_of_new_infections_actual_asym4[Prob_V_asym4>random_numbers_asym4] = dT*number_infectee_not_wearing_mask_at_sampled_times[Prob_V_asym4>random_numbers_asym4]*Prob_V_asym4[Prob_V_asym4>random_numbers_asym4]
                    
                    
                    number_of_new_infections_actual_asym=number_of_new_infections_actual_asym1+
                      number_of_new_infections_actual_asym2+
                      number_of_new_infections_actual_asym3+
                      number_of_new_infections_actual_asym4
                    
                    type_mask_to_unmask_asym=sum(number_of_new_infections_actual_asym1)
                    type_mask_to_mask_asym=sum(number_of_new_infections_actual_asym2)
                    type_unmask_to_mask_asym = sum(number_of_new_infections_actual_asym3)
                    type_unmask_to_unmask_asym = sum(number_of_new_infections_actual_asym4)
                    
                  } else {
                    ## when only exposed contact wearing mask
                    V_asym3 = (1-mask_infectee_efficacy)*V_asym3
                    V_asym3[is.na(V_asym3)]=0
                    Prob_V_asym3=(pmax(0,(V_asym3)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_asym3)))^alpha)
                    ## when none is wearing mask
                    V_asym4 = V_asym4
                    V_asym4[is.na(V_asym4)]=0
                    Prob_V_asym4=(pmax(0,(V_asym4)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_asym4)))^alpha)
                    
                    # From unmasked transmitter to masked contacts
                    number_of_new_infections_actual_asym3=rep(0,length=length(Non_replication_phase_duration_infector_round_discrete))
                    number_of_new_infections_actual_asym3[Prob_V_asym3>random_numbers_asym3] = dT*number_infectee_wearing_mask_at_sampled_times[Prob_V_asym3>random_numbers_asym3]*Prob_V_asym3[Prob_V_asym3>random_numbers_asym3]
                    # From unmasked transmitter to unmasked contacts
                    number_of_new_infections_actual_asym4=rep(0,length=length(Non_replication_phase_duration_infector_round_discrete))
                    number_of_new_infections_actual_asym4[Prob_V_asym4>random_numbers_asym4] = dT*number_infectee_not_wearing_mask_at_sampled_times[Prob_V_asym4>random_numbers_asym4]*Prob_V_asym4[Prob_V_asym4>random_numbers_asym4]
                    
                    number_of_new_infections_actual_asym= number_of_new_infections_actual_asym3+
                      number_of_new_infections_actual_asym4
                    
                    type_mask_to_unmask_asym=0
                    type_mask_to_mask_asym=0
                    type_unmask_to_mask_asym = sum(number_of_new_infections_actual_asym3)
                    type_unmask_to_unmask_asym = sum(number_of_new_infections_actual_asym4)
                    
                  }
                  
                  
                  
                  time_new_infections_actual_asym = Non_replication_phase_duration_infector_round_discrete[number_of_new_infections_actual_asym>0]
                  

                  R0_asym_infections[ijx]=round(sum(number_of_new_infections_actual_asym))
                  asymptomatic_incidence[ijx]=round(sum(number_of_new_infections_actual_asym))
                  
                  if(R0_asym_infections[ijx]>0){  ### if infector started spreading infection in asym phase
                    
                    #print("if infector started spreading infection in asym phase")
                    
                    ##### Now for period after viral load kicks off    ---- if sum(number_of_new_infections_actual_asym)==0 else we get serial interval and generation time directly from asymtpmatic phase
                    init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
                    t.out <- seq(tzero,tzero+30,by=dT)
                    #length(t.out)
                    params=c()
                    out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
                    
                    t.out2 <- seq(tzero,tzero+30,by=0.01)
                    out2 <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
                    #print(out$V)
                    
                    out$V[out$time>tzero+20]=0
                    V_sym=pmax(0,(out$V))
                    
                    no_contact_per_day_options_mid<-c()
                    no_contact_per_day_options_mid=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[dispersion]),
                                                          scale=dispersion_options[dispersion]) # Random variable X is distributed X=gamma(A,B) with mean =AB 
                    
                    no_contact_per_day_options<-c()
                    no_contact_per_day_options=no_contact_per_day_options_mid[seq(1,length(is.finite(V_sym)),by=1)]
      
                    
                    V_sym1 = V_sym ## when only transmitter wearing mask
                    V_sym2=V_sym ## when both wearing mask
                    V_sym3=V_sym ## when only exposed contact wearing mask
                    V_sym4=V_sym ## when none is wearing mask
                    
                    perc_infectee_wearing_mask_at_all_times=perc_infector_wearing_mask*perc_time_infectee_wearing_mask

                    if(perc_infectee_wearing_mask_at_all_times>0){
                      frac_infectee_wearing_mask_average = (1-(perc_infectee_wearing_mask_at_all_times))
                      frac_infectee_wearing_mask_at_sampled_times = rbinom(n=length(no_contact_per_day_options),size=1,prob=(frac_infectee_wearing_mask_average))
                      number_infectee_wearing_mask_at_sampled_times =  no_contact_per_day_options * frac_infectee_wearing_mask_at_sampled_times 
                      number_infectee_not_wearing_mask_at_sampled_times =  no_contact_per_day_options * (1-frac_infectee_wearing_mask_at_sampled_times)

                    } else{
                      number_infectee_wearing_mask_at_sampled_times =  no_contact_per_day_options * 0 
                      number_infectee_not_wearing_mask_at_sampled_times =  no_contact_per_day_options * 1
                      
                    }

                    
                    random_numbers_sym1= runif(n=length(V_sym1),min=0,max=1)  # uniform distribution on the interval from min to max.
                    random_numbers_sym2= runif(n=length(V_sym2),min=0,max=1)  # uniform distribution on the interval from min to max.
                    random_numbers_sym3= runif(n=length(V_sym3),min=0,max=1)  # uniform distribution on the interval from min to max.
                    random_numbers_sym4= runif(n=length(V_sym4),min=0,max=1)  # uniform distribution on the interval from min to max.
                    
                    
                    if(is_infector_wearing_mask==1){
                      
                      infector_wearing_mask_when= rbinom(n=length(no_contact_per_day_options),size=1,prob=(perc_time_infector_wearing_mask)) # 1 means wearing mask and will not transmit at that time
                      
                      ## times when only transmitter wearing mask
                      V_sym1[which(infector_wearing_mask_when==1)] = (1-mask_infector_efficacy)*
                        V_sym1[which(infector_wearing_mask_when==1)]  ## Viral load for infectiousness/contagiousness at the time of wearing mask should be reduced with mask efficacy
                      V_sym1[which(infector_wearing_mask_when==0)] = 0
                      V_sym1[is.na(V_sym1)]=0
                      
                      Prob_V_sym1=(pmax(0,(V_sym1)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym1)))^alpha)
                      
                      
                      ## times when both wearing mask
                      V_sym2[which(infector_wearing_mask_when==1)] = (1-mask_infectee_efficacy)*
                        (1-mask_infector_efficacy)*
                        V_sym2[which(infector_wearing_mask_when==1)]
                      V_sym2[which(infector_wearing_mask_when==0)]=0
                      V_sym2[is.na(V_sym2)]=0
                      
                      Prob_V_sym2=(pmax(0,(V_sym2)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym2)))^alpha)
                      
                      ## when only exposed contact wearing mask and transmisio is from tramitter at time when not wearing mask
                      V_sym3[which(infector_wearing_mask_when==0)] = (1-mask_infectee_efficacy)*
                        V_sym3[which(infector_wearing_mask_when==0)]
                      V_sym3[which(infector_wearing_mask_when==1)]=0
                      V_sym3[is.na(V_sym3)]=0
                      
                      Prob_V_sym3=(pmax(0,(V_sym3)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym3)))^alpha)
                      
                      ## when none is wearing mask and transmisio is from tramitter at time when not wearing mask
                      V_sym4[which(infector_wearing_mask_when==0)] = V_sym4[which(infector_wearing_mask_when==0)]
                      V_sym4[which(infector_wearing_mask_when==1)]=0
                      V_sym4[is.na(V_sym4)]=0
                      Prob_V_sym4=(pmax(0,(V_sym4)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym4)))^alpha)
                      
                      # From masked transmitter at times when wearing mask to unmasked contacts
                      number_of_new_infections_actual_sym1=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym1[Prob_V_sym1>random_numbers_sym1] = dT*number_infectee_not_wearing_mask_at_sampled_times[Prob_V_sym1>random_numbers_sym1]*Prob_V_sym1[Prob_V_sym1>random_numbers_sym1]
                      # From masked transmitter at times when wearing mask  to masked contacts
                      number_of_new_infections_actual_sym2=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym2[Prob_V_sym2>random_numbers_sym2] = dT*number_infectee_wearing_mask_at_sampled_times[Prob_V_sym2>random_numbers_sym2]*Prob_V_sym2[Prob_V_sym2>random_numbers_sym2]
                      # From masked transmitter at times when not wearing mask to masked contacts
                      number_of_new_infections_actual_sym3=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym3[Prob_V_sym3>random_numbers_sym3] = dT*number_infectee_wearing_mask_at_sampled_times[Prob_V_sym3>random_numbers_sym3]*Prob_V_sym3[Prob_V_sym3>random_numbers_sym3]
                      # From masked transmitter at times when not wearing mask  to unmasked contacts
                      number_of_new_infections_actual_sym4=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym4[Prob_V_sym4>random_numbers_sym4] = dT*number_infectee_not_wearing_mask_at_sampled_times[Prob_V_sym4>random_numbers_sym4]*Prob_V_sym4[Prob_V_sym4>random_numbers_sym4]
                      
                      
                      number_of_new_infections_actual=number_of_new_infections_actual_sym1+
                        number_of_new_infections_actual_sym2+
                        number_of_new_infections_actual_sym3+
                        number_of_new_infections_actual_sym4
                      
                      type_mask_to_unmask_sym=sum(number_of_new_infections_actual_sym1)
                      type_mask_to_mask_sym=sum(number_of_new_infections_actual_sym2)
                      type_unmask_to_mask_sym = sum(number_of_new_infections_actual_sym3)
                      type_unmask_to_unmask_sym = sum(number_of_new_infections_actual_sym4)
                      
                    } else {
                      ## when only exposed contact wearing mask
                      V_sym3 = (1-mask_infectee_efficacy)*V_sym3
                      V_sym3[is.na(V_sym3)]=0
                      Prob_V_sym3=(pmax(0,(V_sym3)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym3)))^alpha)
                      
                      ## when none is wearing mask
                      V_sym4 = V_sym4
                      V_sym4[is.na(V_sym4)]=0
                      Prob_V_sym4=(pmax(0,(V_sym4)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym4)))^alpha)
                      
                      # From unmasked transmitter to masked contacts
                      number_of_new_infections_actual_sym3=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym3[Prob_V_sym3>random_numbers_sym3] = dT*number_infectee_wearing_mask_at_sampled_times[Prob_V_sym3>random_numbers_sym3]*Prob_V_sym3[Prob_V_sym3>random_numbers_sym3]
                      # From unmasked transmitter to unmasked contacts
                      number_of_new_infections_actual_sym4=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym4[Prob_V_sym4>random_numbers_sym4] = dT*number_infectee_not_wearing_mask_at_sampled_times[Prob_V_sym4>random_numbers_sym4]*Prob_V_sym4[Prob_V_sym4>random_numbers_sym4]
                      
                      number_of_new_infections_actual= number_of_new_infections_actual_sym3+
                        number_of_new_infections_actual_sym4
                      
                      type_mask_to_unmask_sym=0
                      type_mask_to_mask_sym=0
                      type_unmask_to_mask_sym = sum(number_of_new_infections_actual_sym3)
                      type_unmask_to_unmask_sym = sum(number_of_new_infections_actual_sym4)
                      
                    }
                    
                    
                    time_new_infections_actual = t.out[number_of_new_infections_actual>0]
                    

                    if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){
                      
                      total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
                      ks_options=seq(1,length(number_of_new_infections_actual),by=1)
                      for (ks in ks_options) {
                        if(ks==1){ 
                          total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
                        } else {
                          total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
                        }
                        
                      }
                      
                      total_number_of_new_infections=round(total_number_of_new_infections)
                      
                    }
                    
                    R0_infections[ijx]=max(total_number_of_new_infections) + R0_asym_infections[ijx]
                    GT_infections[ijx]=time_new_infections_actual_asym[1]
                    
                    time_from_infection_to_symptoms_infector=incubation_period_infector
                    time_from_infection_to_symptoms_infectee=incubation_period_infectee
                    
                    SI_infections[ijx]=round(time_new_infections_actual_asym[1]
                                             +time_from_infection_to_symptoms_infectee
                                             -time_from_infection_to_symptoms_infector)
                    
                    ########################################################################
                  }  else { ###  infector does not spread infection in asym phase
                    
                    
                    #print("infector does not spread infection in asym phase")
                    
                    init.x <- c(S=S_0,I=I_0,V=V_0,E=E_0,M1=M1_0,M2=M2_0)
                    t.out <- seq(tzero,tzero+30,by=dT)
                    #print(t.out)
                    params=c()
                    out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
                    #print(out$V)
                    #print(out)
                    out$V[out$time>tzero+20]=0
                    V_sym=pmax(0,(out$V))
                    
                    no_contact_per_day_options_mid<-c()
                    no_contact_per_day_options_mid=rgamma(n = length(t.out), shape = (no_contact_per_day/dispersion_options[dispersion]),
                                                          scale=dispersion_options[dispersion]) # Random variable X is distributed X=gamma(A,B) with mean =AB 
                    
                    no_contact_per_day_options<-c()
                    no_contact_per_day_options=no_contact_per_day_options_mid[seq(1,length(is.finite(V_sym)),by=1)]
                    
                    
                    V_sym1 = V_sym ## when only trnamitter wearing mask
                    V_sym2=V_sym ## when both wearing mask
                    V_sym3=V_sym ## when only exposed contact wearing mask
                    V_sym4=V_sym ## when none is wearing mask
                    
                    perc_infectee_wearing_mask_at_all_times=perc_infector_wearing_mask*perc_time_infectee_wearing_mask

                    if(perc_infectee_wearing_mask_at_all_times>0){
                      frac_infectee_wearing_mask_average = (1-(perc_infectee_wearing_mask_at_all_times))
                      frac_infectee_wearing_mask_at_sampled_times = rbinom(n=length(no_contact_per_day_options),size=1,prob=(frac_infectee_wearing_mask_average))
                      number_infectee_wearing_mask_at_sampled_times =  no_contact_per_day_options * frac_infectee_wearing_mask_at_sampled_times 
                      number_infectee_not_wearing_mask_at_sampled_times =  no_contact_per_day_options * (1-frac_infectee_wearing_mask_at_sampled_times)

                    } else{
                      number_infectee_wearing_mask_at_sampled_times =  no_contact_per_day_options * 0 
                      number_infectee_not_wearing_mask_at_sampled_times =  no_contact_per_day_options * 1
                      
                    }
                  
                    random_numbers_sym1= runif(n=length(V_sym1),min=0,max=1)  # uniform distribution on the interval from min to max.
                    random_numbers_sym2= runif(n=length(V_sym2),min=0,max=1)  # uniform distribution on the interval from min to max.
                    random_numbers_sym3= runif(n=length(V_sym3),min=0,max=1)  # uniform distribution on the interval from min to max.
                    random_numbers_sym4= runif(n=length(V_sym4),min=0,max=1)  # uniform distribution on the interval from min to max.
                    
                    
                    if(is_infector_wearing_mask==1){
                      
                      infector_wearing_mask_when= rbinom(n=length(no_contact_per_day_options),size=1,prob=(perc_time_infector_wearing_mask)) # 1 means wearing mask and will not transmit at that time
                      
                      ## times when only transmitter wearing mask
                      V_sym1[which(infector_wearing_mask_when==1)] = (1-mask_infector_efficacy)*
                        V_sym1[which(infector_wearing_mask_when==1)]  ## Viral load for infectiousness/contagiousness at the time of wearing mask should be reduced with mask efficacy
                      V_sym1[which(infector_wearing_mask_when==0)] = 0
                      V_sym1[is.na(V_sym1)]=0
                      
                      Prob_V_sym1=(pmax(0,(V_sym1)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym1)))^alpha)
                      
                      
                      ## times when both wearing mask
                      V_sym2[which(infector_wearing_mask_when==1)] = (1-mask_infectee_efficacy)*
                        (1-mask_infector_efficacy)*
                        V_sym2[which(infector_wearing_mask_when==1)]
                      V_sym2[which(infector_wearing_mask_when==0)]=0
                      V_sym2[is.na(V_sym2)]=0
                      
                      Prob_V_sym2=(pmax(0,(V_sym2)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym2)))^alpha)
                      
                      ## when only exposed contact wearing mask and transmisio is from tramitter at time when not wearing mask
                      V_sym3[which(infector_wearing_mask_when==0)] = (1-mask_infectee_efficacy)*
                        V_sym3[which(infector_wearing_mask_when==0)]
                      V_sym3[which(infector_wearing_mask_when==1)]=0
                      V_sym3[is.na(V_sym3)]=0
                      
                      Prob_V_sym3=(pmax(0,(V_sym3)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym3)))^alpha)
                      
                      ## when none is wearing mask and transmisio is from tramitter at time when not wearing mask
                      V_sym4[which(infector_wearing_mask_when==0)] = V_sym4[which(infector_wearing_mask_when==0)]
                      V_sym4[which(infector_wearing_mask_when==1)]=0
                      V_sym4[is.na(V_sym4)]=0
                      Prob_V_sym4=(pmax(0,(V_sym4)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym4)))^alpha)
                      
                      # From masked transmitter at times when wearing mask to unmasked contacts
                      number_of_new_infections_actual_sym1=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym1[Prob_V_sym1>random_numbers_sym1] = dT*number_infectee_not_wearing_mask_at_sampled_times[Prob_V_sym1>random_numbers_sym1]*Prob_V_sym1[Prob_V_sym1>random_numbers_sym1]
                      # From masked transmitter at times when wearing mask  to masked contacts
                      number_of_new_infections_actual_sym2=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym2[Prob_V_sym2>random_numbers_sym2] = dT*number_infectee_wearing_mask_at_sampled_times[Prob_V_sym2>random_numbers_sym2]*Prob_V_sym2[Prob_V_sym2>random_numbers_sym2]
                      # From masked transmitter at times when not wearing mask to masked contacts
                      number_of_new_infections_actual_sym3=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym3[Prob_V_sym3>random_numbers_sym3] = dT*number_infectee_wearing_mask_at_sampled_times[Prob_V_sym3>random_numbers_sym3]*Prob_V_sym3[Prob_V_sym3>random_numbers_sym3]
                      # From masked transmitter at times when not wearing mask  to unmasked contacts
                      number_of_new_infections_actual_sym4=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym4[Prob_V_sym4>random_numbers_sym4] = dT*number_infectee_not_wearing_mask_at_sampled_times[Prob_V_sym4>random_numbers_sym4]*Prob_V_sym4[Prob_V_sym4>random_numbers_sym4]
                      
                      
                      number_of_new_infections_actual=number_of_new_infections_actual_sym1+
                        number_of_new_infections_actual_sym2+
                        number_of_new_infections_actual_sym3+
                        number_of_new_infections_actual_sym4
                      
                      type_mask_to_unmask_sym=sum(number_of_new_infections_actual_sym1)
                      type_mask_to_mask_sym=sum(number_of_new_infections_actual_sym2)
                      type_unmask_to_mask_sym = sum(number_of_new_infections_actual_sym3)
                      type_unmask_to_unmask_sym = sum(number_of_new_infections_actual_sym4)
                      
                    } else {
                      ## when only exposed contact wearing mask
                      V_sym3 = (1-mask_infectee_efficacy)*V_sym3
                      V_sym3[is.na(V_sym3)]=0
                      Prob_V_sym3=(pmax(0,(V_sym3)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym3)))^alpha)
                      
                      ## when none is wearing mask
                      V_sym4 = V_sym4
                      V_sym4[is.na(V_sym4)]=0
                      Prob_V_sym4=(pmax(0,(V_sym4)))^alpha/((10^AUC_factor)^alpha   +   (pmax(0,(V_sym4)))^alpha)
                      
                      # From unmasked transmitter to masked contacts
                      number_of_new_infections_actual_sym3=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym3[Prob_V_sym3>random_numbers_sym3] = dT*number_infectee_wearing_mask_at_sampled_times[Prob_V_sym3>random_numbers_sym3]*Prob_V_sym3[Prob_V_sym3>random_numbers_sym3]
                      # From unmasked transmitter to unmasked contacts
                      number_of_new_infections_actual_sym4=rep(0,length=length(t.out))
                      number_of_new_infections_actual_sym4[Prob_V_sym4>random_numbers_sym4] = dT*number_infectee_not_wearing_mask_at_sampled_times[Prob_V_sym4>random_numbers_sym4]*Prob_V_sym4[Prob_V_sym4>random_numbers_sym4]
                      
                      number_of_new_infections_actual= number_of_new_infections_actual_sym3+
                        number_of_new_infections_actual_sym4
                      
                      type_mask_to_unmask_sym=0
                      type_mask_to_mask_sym=0
                      type_unmask_to_mask_sym = sum(number_of_new_infections_actual_sym3)
                      type_unmask_to_unmask_sym = sum(number_of_new_infections_actual_sym4)
                      
                    }
                    
                    
                    time_new_infections_actual = t.out[number_of_new_infections_actual>0]
         
                    if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5){
                      
                      total_number_of_new_infections<- vector() # the addition is done to avoid classify 10^-6 as new infection
                      ks_options=seq(1,length(number_of_new_infections_actual),by=1)
                      for (ks in ks_options) {
                        if(ks==1){ 
                          total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]
                        } else {
                          total_number_of_new_infections[ks]=number_of_new_infections_actual[ks]+total_number_of_new_infections[ks-1] 
                        }
                        
                      }
                      
             
                      total_number_of_new_infections=round(total_number_of_new_infections)
                      
                      #### as zero do not play any role, we remove them
                      
                      time_new_infections_actual_after_randomization=t.out[total_number_of_new_infections>0]
                      total_number_of_new_infections=total_number_of_new_infections[total_number_of_new_infections>0]
                      
                      index_where_new_infections_happen=match(unique(total_number_of_new_infections),total_number_of_new_infections)
                      
                      time_when_new_infection_is_added=time_new_infections_actual_after_randomization[index_where_new_infections_happen]  ##
                      
                      time_when_new_infection_is_added_relative_to_tzero=time_when_new_infection_is_added  
                      
                      cumulative_incidence_temp=total_number_of_new_infections[index_where_new_infections_happen]
                      
                      ### finding incidence every day
                      incidence<-vector()
                      y_length=seq(1,length(cumulative_incidence_temp),by=1)
                      for (yid in y_length){
                        if(yid==1){
                          incidence[yid]=cumulative_incidence_temp[yid]
                        } else{
                          incidence[yid]=cumulative_incidence_temp[yid]-cumulative_incidence_temp[yid-1]  ################ this is imp
                        }
                      }
                      
                      #### presymptomatic vs symptomatic incidence
                      
                      
                      ind_presympt = which(time_when_new_infection_is_added+tzero-incubation_period_infector<=0)

                      presymptomatic_incidence[ijx]=sum(incidence[ind_presympt])   ################ this is imp
                      
                      ind_sympt = which(time_when_new_infection_is_added+tzero-incubation_period_infector>0)    ################ this is imp
                      symptomatic_incidence[ijx]=sum(incidence[ind_sympt])
                   
                      #### Determine reproduction number
                      R0_infections[ijx] = sum(incidence)
                      
                      #### time when the first secondary infection happens
                      if(isempty(incidence)){
                        print("AA")
                        SI_infections[ijx]= NA # maximum time period if no new infections happen
                        R0_infections[ijx]=0
                        GT_infections[ijx]=NA
                        
                        
                        ID_temp <-matrix(ijx,nrow=length(t.out),ncol=1)
                        incidence_temp<- matrix(0,nrow=length(t.out),ncol=1)
                        R0_temp<- matrix(0,nrow=length(t.out),ncol=1)
                        R0_cumulative_temp<- matrix(0,nrow=length(t.out),ncol=1)
                        alpha_temp<- matrix(alpha_options[ijx],nrow=length(t.out),ncol=1)
                        theta_temp<- matrix(0,nrow=length(t.out),ncol=1)
                        lambda_temp<- matrix(lambda_options[ijx],nrow=length(t.out),ncol=1)
                        
                        pre_temp<-  matrix(0,nrow=length(t.out),ncol=1)
                        pre_temp= t.out<round(incubation_period_infector,digits=1)
  
                        theta_temp[seq(1,length(no_contact_per_day_options))]=(round(no_contact_per_day_options))

                        V_temp<- matrix(0,nrow=length(t.out),ncol=1)
                        V_temp[seq(1,length(V_sym))]=pmax(0,log10(V_sym))
   
                        R0_temp_all_infections[c(seq(1+(ijx-1)*length(t.out),length(t.out)+(ijx-1)*length(t.out),by=1)),]=c(ID_temp,t.out,V_temp,theta_temp,lambda_temp,alpha_temp,R0_temp,R0_cumulative_temp,pre_temp)

                        
                      } else{
                        
                        
                        R0_infections[ijx] = sum(incidence)
                        
                        time_from_infection_to_VL_kick_off_infector=Non_replication_phase_duration_infector
                        time_from_VL_kick_off_infector_to_secondary_infection= time_when_new_infection_is_added[1] 
                        
                        GT_infections[ijx]=round(time_from_infection_to_VL_kick_off_infector
                                                 +time_from_VL_kick_off_infector_to_secondary_infection)
                        
                        time_from_infection_to_symptoms_infector=incubation_period_infector
                        time_from_infection_to_symptoms_infectee=incubation_period_infectee
                        
                        SI_infections[ijx]=round(time_from_infection_to_VL_kick_off_infector
                                                 +time_from_VL_kick_off_infector_to_secondary_infection
                                                 +time_from_infection_to_symptoms_infectee
                                                 -time_from_infection_to_symptoms_infector)
                        

                        IND_a=is.element(t.out,time_when_new_infection_is_added_relative_to_tzero)
                        
                        ID_temp <-matrix(ijx,nrow=length(t.out),ncol=1)
                        incidence_temp<- matrix(0,nrow=length(t.out),ncol=1)
                        R0_temp<- matrix(0,nrow=length(t.out),ncol=1)
                        R0_cumulative_temp<- matrix(sum(incidence),nrow=length(t.out),ncol=1)
                        alpha_temp<- matrix(alpha_options[ijx],nrow=length(t.out),ncol=1)
                        theta_temp<- matrix(0,nrow=length(t.out),ncol=1)
                        lambda_temp<- matrix(lambda_options[ijx],nrow=length(t.out),ncol=1)
                        
                        
                        pre_temp<-  matrix(0,nrow=length(t.out),ncol=1)
                        pre_temp= t.out<round(incubation_period_infector,digits=1)
                        

                        theta_temp[seq(1,length(no_contact_per_day_options))]=(round(no_contact_per_day_options))
                        R0_temp[IND_a]=incidence
                        
                        V_temp<- matrix(0,nrow=length(t.out),ncol=1)
                        V_temp[seq(1,length(V_sym))]=pmax(0,log10(V_sym))
                        
                        
                        R0_temp_all_infections[c(seq(1+(ijx-1)*length(t.out),length(t.out)+(ijx-1)*length(t.out),by=1)),]=c(ID_temp,t.out,V_temp,theta_temp,lambda_temp,alpha_temp,R0_temp,R0_cumulative_temp,pre_temp)
                        
                        
                      }
                      
                    }  else { ### this checks whether infections happened in symptomatic phase or not
                      SI_infections[ijx]= NA # maximum time period if no new infections happen
                      R0_infections[ijx]=0
                      GT_infections[ijx]=NA
                      
                      
                      ID_temp <-matrix(ijx,nrow=length(t.out),ncol=1)
                      incidence_temp<- matrix(0,nrow=length(t.out),ncol=1)
                      R0_temp<- matrix(0,nrow=length(t.out),ncol=1)
                      R0_cumulative_temp<- matrix(0,nrow=length(t.out),ncol=1)
                      alpha_temp<- matrix(alpha_options[ijx],nrow=length(t.out),ncol=1)
                      theta_temp<- matrix(0,nrow=length(t.out),ncol=1)
                      lambda_temp<- matrix(lambda_options[ijx],nrow=length(t.out),ncol=1)
                      
                      
                      pre_temp<-  matrix(0,nrow=length(t.out),ncol=1)
                      pre_temp= t.out<round(incubation_period_infector,digits=1)
                      

                      theta_temp[seq(1,length(no_contact_per_day_options))]=(round(no_contact_per_day_options))

                      V_temp<- matrix(1,nrow=length(t.out),ncol=1)
                      V_temp[seq(1,length(V_sym))]=pmax(0,log10(V_sym))

                      
                      R0_temp_all_infections[c(seq(1+(ijx-1)*length(t.out),length(t.out)+(ijx-1)*length(t.out),by=1)),]=c(ID_temp,t.out,V_temp,theta_temp,lambda_temp,alpha_temp,R0_temp,R0_cumulative_temp,pre_temp)

                      
                      
                    } # end of  else of if(length(number_of_new_infections_actual)>0  && sum(number_of_new_infections_actual)>0.5)
                    
                    
                    
                    
                    
                    
                    ijx=ijx+1
                    
                    
                  } # here is the end of else of if(R0_asym_infections[ijx]>0)
                  
                  
                }  #### end of yu --- 1000 hypothetical individuals
                
                
                mean_R0=round(mean(R0_infections),digits = 1)
                mean_R0_iter[iter] = mean_R0
                
                
                SI_infections_new=na.omit(SI_infections)
                #print(SI_infections)
                if(length(SI_infections_new)==0){
                  mean_SI=0
                } else { 
                  mean_SI=mean(SI_infections_new,na.rm=TRUE)
                }
                mean_SI_iter[iter] = mean_SI
                
                
                GT_infections_new=na.omit(GT_infections)
                #print(SI_infections)
                if(length(GT_infections_new)==0){
                  mean_GT=0
                } else { 
                  mean_GT=mean(GT_infections_new,na.rm=TRUE)
                }
                mean_GT_iter[iter] = mean_GT
                
                
                mean_Trans_iter[iter] = 100*length(R0_infections[which(R0_infections>0)])/max_replicates
                mean_SSE5_of_all_iter[iter] = 100*length(R0_infections[which(R0_infections>=5)])/max_replicates
                mean_SSE10_of_all_iter[iter] = 100*length(R0_infections[which(R0_infections>=10)])/max_replicates
                mean_SSE20_of_all_iter[iter] = 100*length(R0_infections[which(R0_infections>=20)])/max_replicates
                
                type_mask_to_unmask_sym=0
                type_mask_to_mask_sym=0
                type_unmask_to_mask_sym = sum(number_of_new_infections_actual_sym3)
                type_unmask_to_unmask_sym = sum(number_of_new_infections_actual_sym4)
                
                
                trans_exposed_mask[iter]= 100*(type_mask_to_mask_sym + type_mask_to_mask_asym) / (type_mask_to_mask_sym + type_mask_to_mask_asym+
                                                                                                type_mask_to_unmask_sym + type_mask_to_unmask_asym+
                                                                                                type_unmask_to_mask_sym + type_unmask_to_mask_asym+
                                                                                                type_unmask_to_unmask_sym + type_unmask_to_unmask_asym)
                trans_mask[iter]=100*(type_mask_to_unmask_sym + type_mask_to_unmask_asym) / (type_mask_to_mask_sym + type_mask_to_mask_asym+
                                                                                       type_mask_to_unmask_sym + type_mask_to_unmask_asym+
                                                                                       type_unmask_to_mask_sym + type_unmask_to_mask_asym+
                                                                                       type_unmask_to_unmask_sym + type_unmask_to_unmask_asym)
                exposed_mask[iter]=100*(type_unmask_to_mask_sym + type_unmask_to_mask_asym) / (type_mask_to_mask_sym + type_mask_to_mask_asym+
                                                                                         type_mask_to_unmask_sym + type_mask_to_unmask_asym+
                                                                                         type_unmask_to_mask_sym + type_unmask_to_mask_asym+
                                                                                         type_unmask_to_unmask_sym + type_unmask_to_unmask_asym)
                none_mask[iter]=100*(type_unmask_to_unmask_sym + type_unmask_to_unmask_asym) / (type_mask_to_mask_sym + type_mask_to_mask_asym+
                                                                                      type_mask_to_unmask_sym + type_mask_to_unmask_asym+
                                                                                      type_unmask_to_mask_sym + type_unmask_to_mask_asym+
                                                                                      type_unmask_to_unmask_sym + type_unmask_to_unmask_asym)
                
                
                ijk=ijk+1
                
                par_ID=par_ID+1
                
                
                
              }
              
              
            }
            
            
            Final_R0=mean(mean_R0_iter)
            Final_SI=mean(mean_SI_iter)
            Final_GT=mean(mean_GT_iter)
            Final_transmitters=mean(mean_Trans_iter)
            Final_SSE5_of_all=mean(mean_SSE5_of_all_iter)
            Final_SSE10_of_all=mean(mean_SSE10_of_all_iter)
            Final_SSE20_of_all=mean(mean_SSE20_of_all_iter)
            Final_trans_exposed_mask=mean(trans_exposed_mask)
            Final_trans_mask=mean(trans_mask)
            Final_exposed_mask=mean(exposed_mask)
            Final_none_mask=mean(none_mask)
            
            #"trans_exposed_mask","trans_mask","exposed_mask","none_mask")
            
            PAR_proj_R0_SI_mask[par_ind_mask,]=c(perc_infector_wearing_mask,perc_infectee_wearing_mask,
                                                 perc_time_infector_wearing_mask,perc_time_infectee_wearing_mask,
                                                 mask_infector_efficacy,mask_infectee_efficacy,Final_R0,Final_SI,Final_GT,
                                                 Final_transmitters,Final_SSE5_of_all,Final_SSE10_of_all,Final_SSE20_of_all,
                                                 Final_trans_exposed_mask,Final_trans_mask,Final_exposed_mask,Final_none_mask)
            
            write.csv(PAR_proj_R0_SI_mask,file="Proj_COVID_R0_SI_mask_use.csv",row.names = FALSE) 
            
            par_ind_mask=par_ind_mask+1
            
            
          }
          
        }
        
      }
    }
  }
}

write.csv(PAR_proj_R0_SI_mask,file="Proj_COVID_R0_SI_mask_use.csv",row.names = FALSE) 

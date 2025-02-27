# LIP model (latent variable model, see Krucien et al. 2017)
# Eye-tracking comparison experiment
# Universidad de Chile - National University of Singapore
# Code by Bastian Henríquez-Jara
rm(list = ls())
#install.packages('apollo')
library(apollo)
library(dplyr)
library(plyr)
detach(package:plyr)    
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cases=c('monitor','laptop','tablet','mobile')
exp='exp1'

normalize<-function(x){
  x=as.numeric(x)
  x[is.infinite(x)]=NA
  x=(x-mean(x,na.rm = TRUE))/sd(x,na.rm = TRUE) #run again without sd!
  return(x)}
max_min<-function(x){
  x=as.numeric(x)
  x[is.infinite(x)]=NA
  x=(x-min(x,na.rm = TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm = TRUE))
  return(x)}
mat_power=function(A,b){
  for(i in 1:ncol(A)){
    A[,i]=A[,i]^b[i]  
  }
  return(A)
}


########
ETs=c('WG','WL')
VST_summary=c()
names=c()
IPs=c()
alpha_beta=c()
LL_list=c()
Cross_validation<-c()


for (ET in ETs){
i=0
for (case in cases){

  path=paste0('../Compiled_data/data_merged_',case,'_',exp,'_',ET,'.csv')
  data=read.csv(path) 
  data=data[!(data$pid %in% c(20,21)),]
  
  data=data[!((data$pid %in% c(17,36) & case=="mobile")|
                (data$pid==28 & case=="tablet")|
                (data$pid %in% c(5,8) & case=="monitor")|
                (data$pid %in% c(2,22) & case=="laptop"))
                ,]
  data$Choice_id=1*(data$Choice=="Bus_label")+2*(data$Choice=="Metro_label")+3*(data$Choice=="RH_label")
  
  
  data=data%>%group_by(pid)%>%mutate(max_freq=max(sum(Choice_id==1),
                                                  sum(Choice_id==2),
                                                  sum(Choice_id==3)))%>%filter(max_freq<12)

  data=data%>%group_by(pid)%>%mutate(inertia_bus=c(0,1*((Choice_id[-length(Choice_id)])==1)),
                                     inertia_metro=c(0,1*((Choice_id[-length(Choice_id)])==2)),
                                     inertia_RH=c(0,1*((Choice_id[-length(Choice_id)])==3)))

  #   data=data%>%mutate(RT=rowSums(cbind(Travel_FT,Cost_FT,Comfort_FT,fixation_counts.Other_WL),na.rm=TRUE))
#   data=data%>%filter(RT>0)
 
  data=data%>%
    mutate(
      Comfort_FT_sum_norm=((Comfort_FT_count)/Tot_counts),#/sum(FT_sum)
      Cost_FT_sum_norm=((Cost_FT_count)/Tot_counts),
      Travel_FT_sum_norm=((Travel_FT_count)/Tot_counts))

  temp=data %>%
    filter(Trial == 1) %>% select(Comfort_FT_count,Cost_FT_count,Travel_FT_count,pid,
                                  Comfort_FT_time,Cost_FT_time,Travel_FT_time,Tot_counts)
  temp$Comfort_FT_norm = normalize(temp$Comfort_FT_count)#/temp$Tot_counts
  temp$Cost_FT_norm = normalize(temp$Cost_FT_count)#/temp$Tot_counts
  temp$Travel_FT_norm = normalize(temp$Travel_FT_count)#/temp$Tot_counts
 
  data=merge(data,temp,by=c('pid'))  

  
  data<-data%>%group_by(pid)%>%mutate(Exp_duration=sum(TrialDuration))
  
  set.seed(i)
  
  Ne=round(length(unique(data$pid))*0.8)
  Pid_estimate=sample(data$pid,Ne)
  
  ## split data
  estimation_data=data[data$pid %in% Pid_estimate,]
  test_data=data[!(data$pid %in% Pid_estimate),]
  
  apollo_initialise()
  #for(i in 101:150){
  ### Initialise code
  
  
  ### Set core controls
  apollo_control = list(
    modelName       = paste0(exp,'_',case,'_',ET,'_LIP_counts'),
    modelDescr      = "LIP for transport choices - ET comparison",
    indivID         = "pid", 
    outputDirectory = "output",
    #workInLogs=TRUE,
    #mixing=TRUE,
    nCores=5,
    panelData=TRUE
  )
  
    
  database = data
  
  apollo_beta = c(b_bus=0,
                  b_metro=0,
                  b_RH=0,
                  b_time=0,
                  b_cost=0,
                  b_comfort=0,
                  b_inertia=0,
                  gamma_IP=0,
                  
                  #gamma_RT=0,
                  gamma_time=1,
                  gamma_cost=1,
                  gamma_comfort=1,
                  
                  gamma0_comfort=0,
                  gamma0_cost=0,
                  gamma0_time=0,
                  
                  alpha_time=0,
                  alpha_cost=0,
                  alpha_comfort=0,
                  
                  sigma_IP=1,
                  sigma_cost=1,
                  sigma_time=1,
                  sigma_comfort=1
  )
  apollo_fixed=c('sigma_IP','b_bus','gamma0_cost')
  
  apollo_draws = list(
    interDrawsType="halton", 
    interNDraws=2000,          
    interUnifDraws=c(),      
    interNormDraws=c("eta")
    
    # intraDrawsType="",
    # intraNDraws=0,          
    # intraUnifDraws=c(),     
    # intraNormDraws=c()      
  )
  
  ### Create random parameters
  apollo_randCoeff=function(apollo_beta, apollo_inputs){
    randcoeff = list()
    
    randcoeff[["IP"]] = gamma_IP+sigma_IP*eta

    
    randcoeff[['b_cost_LV']] = b_cost+alpha_cost*randcoeff[["IP"]]
    randcoeff[['b_time_LV']]=b_time+alpha_time*randcoeff[["IP"]]
    randcoeff[['b_comfort_LV']]=b_comfort+alpha_comfort*randcoeff[["IP"]]
    
    return(randcoeff)
  }
  
  
  apollo_inputs = apollo_validateInputs()
  
  apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
    
    ### Attach inputs and detach after function exit
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    
    ### Create list of probabilities P
    P = list()
    
    V=list()
    V[['bus']]=b_bus+(b_cost_LV)*Bus_cost+(b_time_LV)*Bus_travel_time+(b_comfort_LV)*Bus_Comfort+b_inertia*inertia_bus
    
    V[['metro']]=b_metro+(b_cost_LV)*metro_cost+(b_time_LV)*metro_travel_time+(b_comfort_LV)*metro_Comfort+b_inertia*inertia_metro
    
    V[['RH']]=b_RH+(b_cost_LV)*RH_cost+(b_time_LV)*RH_travel_time+(b_comfort_LV)*RH_Comfort+b_inertia*inertia_RH
  
    mnl_settings = list(
      utilities=V,
      alternatives  = c(bus=1,metro=2,RH=3),
      avail=1,
      choiceVar = Choice_id
    )
  
    normalDensity_settings1 = list(outcomeNormal = (Cost_FT_norm), #valor real (observado)
                                   xNormal       = gamma0_cost+IP*gamma_cost, #valor estiamdo
                                   mu            = 0,
                                   sigma         = abs(sigma_cost), #error estandar de la ecuacion de medicion
                                   rows          = (Trial==1),
                                   componentName = "cost")
    normalDensity_settings2 = list(outcomeNormal = (Travel_FT_norm),
                                   xNormal       = IP*gamma_time,
                                   mu            = gamma0_time,
                                   sigma         = abs(sigma_time),
                                   rows          = (Trial==1),
                                   componentName = "time")
    normalDensity_settings3 = list(outcomeNormal = (Comfort_FT_norm),
                                   xNormal       = IP*gamma_comfort,
                                   mu            = gamma0_comfort,
                                   sigma         = abs(sigma_comfort),
                                   rows          = (Trial==1),
                                   componentName = "comfort")
    
    P[["cost"]]     = apollo_normalDensity(normalDensity_settings1, functionality)
    P[["time"]] = apollo_normalDensity(normalDensity_settings2, functionality)
    P[["comfort"]]      = apollo_normalDensity(normalDensity_settings3, functionality)
    
      ### Compute within-class choice probabilities using MNL model
    P[['choice']] = apollo_mnl(mnl_settings, functionality)
      
    
    P = apollo_combineModels(P, apollo_inputs, functionality)
    ### Take product across observation for same individual
    #P = apollo_avgIntraDraws(P, apollo_inputs, functionality)
    P = apollo_panelProd(P, apollo_inputs, functionality)
    ### Average across inter-individual draws
    P = apollo_avgInterDraws(P, apollo_inputs, functionality)
    ### Prepare and return outputs of function
    P = apollo_prepareProb(P, apollo_inputs, functionality)
    
    return(P)
  }
  
  
  
  # ################################################################# #
  #### MODEL ESTIMATION                                            ####
  # ################################################################# #
  
  model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                          estimate_settings=list(maxIterations=400))
  
  #
  apollo_modelOutput(model)
  
  # ----------------------------------------------------------------- #
  #---- FORMATTED OUTPUT (TO FILE, using model name)               ----
  # ----------------------------------------------------------------- #
  # betas<-rbind(betas,model$estimate)
  # robse<-rbind(robse,model$robse)
  # }
  
  apollo_saveOutput(model)
  
#}

## save LL of choices
  

LL_0<-apollo_llCalc(
  apollo_beta = apollo_beta, 
  apollo_probabilities = model$apollo_probabilities, 
  apollo_inputs = apollo_inputs
)

LL_final<-apollo_llCalc(
  apollo_beta = model$estimate, 
  apollo_probabilities = model$apollo_probabilities, 
  apollo_inputs = apollo_inputs
)

LL_list<-rbind(LL_list,cbind(exp=exp,device=case,ET=ET,
                 LL0=LL_0$choice,LLfinal=LL_final$choice))

## save posteriors of random parameters
cond=apollo_conditionals(model,apollo_probabilities,apollo_inputs)
b_cost_post=cond$b_cost_LV
b_time_post=cond$b_time_LV
IP_temp=cond$IP
temp=cbind(exp=exp,device=case,ET=ET,IP_temp,
           IP_alpha_time=model$estimate['alpha_time']*IP_temp$post.mean,
           IP_alpha_cost=model$estimate['alpha_cost']*IP_temp$post.mean,
           IP_alpha_comfort=model$estimate['alpha_comfort']*IP_temp$post.mean)
IPs=rbind(IPs,temp)
## calculate VST
VST=(b_time_post$post.mean/b_cost_post$post.mean)

VST=data.frame(VST=VST,exp=exp,device=case,ET=ET)
VST_summary=rbind(VST_summary,VST)
#alpha/beta
delta_settings=list(expression=c('alpha_time/b_time',
                                 'alpha_cost/b_cost',
                                 'alpha_comfort/b_comfort'))

temp=apollo_deltaMethod(model,deltaMethod_settings = delta_settings)
temp=cbind(temp,alpha_mean=c(model$estimate['alpha_time'],
                        model$estimate['alpha_cost'],
                        model$estimate['alpha_comfort']),
           alpha_se=c(model$robse['alpha_time'],
                      model$robse['alpha_cost'],
                      model$robse['alpha_comfort']),
           var=c('time','cost','comfort'),
           exp=exp,device=case,ET=ET)
alpha_beta=rbind(alpha_beta,temp)


### Replace `database` with the test set
# database = test_data
# apollo_inputs = apollo_validateInputs()

### Cross-validation
# predictions = apollo_prediction(model, apollo_probabilities, apollo_inputs)
# LL<-mean(log(predictions$choice$chosen))
# 
# FPR<- sum(predictions$choice[,6]==apply(predictions$choice[,c(3:5)],1,max))/nrow(predictions$choice)
# 
# Cross_validation<-rbind(Cross_validation,
#                         cbind(exp=exp,device=case,ET=ET,LL=LL,FPR=FPR))
# outOfSample_settings=list(nRep=5,
#                           validationSize=0.1)
# apollo_outOfSample(apollo_beta,
#                    apollo_fixed,
#                    apollo_probabilities,
#                    apollo_inputs,
#                    estimate_settings=list(maxIterations=400),outOfSample_settings)

i=i+1
}
}


###### out of sample ########


# 
# ############### Baseline ###############
# # 
# ET='Baseline'
# i=0
# for (case in cases){
#   path=paste0('../Compiled_data/data_merged_',case,'_',exp,'_','WL','.csv')
#   data=read.csv(path)
#   data=data[!(data$pid %in% c(20,21)),]
#   data=data[!((data$pid %in% c(17,36) & case=="mobile")|
#                 (data$pid==28 & case=="tablet")|
#                 (data$pid %in% c(5,8) & case=="monitor")|
#                 (data$pid %in% c(2,20,21,22) & case=="laptop"))
#             ,]
#   data$Choice_id=1*(data$Choice=="Bus_label")+2*(data$Choice=="Metro_label")+3*(data$Choice=="RH_label")
#   #robse=c()
#   data=data%>%group_by(pid)%>%mutate(max_freq=max(sum(Choice_id==1),
#                                                   sum(Choice_id==2),
#                                                   sum(Choice_id==3)))%>%filter(max_freq<12)
#   
#   data=data%>%group_by(pid)%>%mutate(inertia_bus=c(0,1*((Choice_id[-length(Choice_id)])==1)),
#                                      inertia_metro=c(0,1*((Choice_id[-length(Choice_id)])==2)),
#                                      inertia_RH=c(0,1*((Choice_id[-length(Choice_id)])==3)))
#   
#   x=data%>%group_by(pid)%>%reframe(x=var(Choice_id))
#   print(x,n=40)
#   set.seed(27012025)
#   Ne=round(length(unique(data$pid))*0.8)
#   Pid_estimate=sample(data$pid,Ne)
#   
#   ## split data
#   estimation_data=data[data$pid %in% Pid_estimate,]
#   test_data=data[!(data$pid %in% Pid_estimate),]
#   
#   #data=merge(data,temp,by=c('pid'))
# 
#   apollo_initialise()
#   #for(i in 101:150){
#   ### Initialise code
# 
# 
#   ### Set core controls
#   apollo_control = list(
#     modelName       = paste0(exp,'_',case,'_',ET,'_MixedLogit'),
#     modelDescr      = "LIP for transport choices - ET comparison",
#     indivID         = "pid",
#     outputDirectory = "output",
#     #workInLogs=TRUE,
#     mixing=TRUE,
#     nCores=5,
#     panelData=TRUE
#   )
# 
# 
#   database = data
# 
#   apollo_beta = c(b_bus=0,
#                   b_metro=0,
#                   b_RH=0,
#                   b_time=0,
#                   b_cost=0,
#                   b_comfort=0,
#                   b_inertia=0,
# 
#                   alpha_time=0,
#                   alpha_cost=0,
#                   alpha_comfort=0
# 
#   )
#   apollo_fixed=c('b_bus')
# 
#   apollo_draws = list(
#     interDrawsType="halton",
#     interNDraws=1600,
#     interUnifDraws=c(),
#     interNormDraws=c("eta1")
# 
#     # intraDrawsType="",
#     # intraNDraws=0,
#     # intraUnifDraws=c(),
#     # intraNormDraws=c()
#   )
# 
#   ### Create random parameters
#   apollo_randCoeff=function(apollo_beta, apollo_inputs){
#     randcoeff = list()
# 
#     randcoeff[['b_cost_LV']] = b_cost+(alpha_cost)*eta1
#     randcoeff[['b_time_LV']]=b_time+(alpha_time)*eta1
#     randcoeff[['b_comfort_LV']]=b_comfort+(alpha_comfort)*eta1
#     # randcoeff[['VA_cost']]=randcoeff[["IP"]]*gamma_cost+sigma_cost*eta1
#     # randcoeff[['VA_time']]=randcoeff[["IP"]]*gamma_time+sigma_time*eta2
#     # randcoeff[['VA_comfort']]=randcoeff[["IP"]]*gamma_comfort+sigma_comfort*eta3
# 
#     return(randcoeff)
#   }
# 
# 
#   apollo_inputs = apollo_validateInputs()
# 
#   apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
# 
#     ### Attach inputs and detach after function exit
#     apollo_attach(apollo_beta, apollo_inputs)
#     on.exit(apollo_detach(apollo_beta, apollo_inputs))
# 
#     ### Create list of probabilities P
#     P = list()
# 
#     V=list()
#     V[['bus']]=b_bus+(b_cost_LV)*Bus_cost+(b_time_LV)*Bus_travel_time+(b_comfort_LV)*Bus_Comfort+b_inertia*inertia_bus
#     
#     V[['metro']]=b_metro+(b_cost_LV)*metro_cost+(b_time_LV)*metro_travel_time+(b_comfort_LV)*metro_Comfort+b_inertia*inertia_metro
#     
#     V[['RH']]=b_RH+(b_cost_LV)*RH_cost+(b_time_LV)*RH_travel_time+(b_comfort_LV)*RH_Comfort+b_inertia*inertia_RH
#     
# 
#     mnl_settings = list(
#       utilities=V,
#       alternatives  = c(bus=1,metro=2,RH=3),
#       avail=1,
#       choiceVar = Choice_id
#     )
# 
# 
# 
# 
# 
#     ### Compute within-class choice probabilities using MNL model
#     P[['model']] = apollo_mnl(mnl_settings, functionality)
# 
# 
#     #P = apollo_combineModels(P, apollo_inputs, functionality)
#     ### Take product across observation for same individual
#     #P = apollo_avgIntraDraws(P, apollo_inputs, functionality)
#     P = apollo_panelProd(P, apollo_inputs, functionality)
#     ### Average across inter-individual draws
#     P = apollo_avgInterDraws(P, apollo_inputs, functionality)
#     ### Prepare and return outputs of function
#     P = apollo_prepareProb(P, apollo_inputs, functionality)
# 
#     return(P)
#   }
# 
# 
# 
#   # ################################################################# #
#   #### MODEL ESTIMATION                                            ####
#   # ################################################################# #
# 
#   model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
# 
#   #,estimate_settings=list(estimationRoutine="bfgs")
#   apollo_modelOutput(model)
# 
#   # ----------------------------------------------------------------- #
#   #---- FORMATTED OUTPUT (TO FILE, using model name)               ----
#   # ----------------------------------------------------------------- #
#   # betas<-rbind(betas,model$estimate)
#   # robse<-rbind(robse,model$robse)
#   # }
# 
#   apollo_saveOutput(model)
# 
#   #}
# 
# 
#   LL_0<-apollo_llCalc(
#     apollo_beta = apollo_beta,
#     apollo_probabilities = model$apollo_probabilities,
#     apollo_inputs = apollo_inputs
#   )
# 
#   LL_final<-apollo_llCalc(
#     apollo_beta = model$estimate,
#     apollo_probabilities = model$apollo_probabilities,
#     apollo_inputs = apollo_inputs
#   )
# 
#   LL_list<-rbind(LL_list,cbind(exp=exp,device=case,ET=ET,
#                                LL0=LL_0$model,LLfinal=LL_final$model))
# 
# 
#   cond=apollo_conditionals(model,apollo_probabilities,apollo_inputs)
#   b_cost_post=cond$b_cost_LV
#   b_time_post=cond$b_time_LV
#   VST=(b_time_post$post.mean/b_cost_post$post.mean)
# 
#   VST=data.frame(VST=VST,exp=exp,device=case,ET=ET)
#   VST_summary=rbind(VST_summary,VST)
#   delta_settings=list(expression=c('alpha_time/b_time',
#                                    'alpha_cost/b_cost',
#                                    'alpha_comfort/b_comfort'))
# 
#   temp=apollo_deltaMethod(model,deltaMethod_settings = delta_settings)
#   temp=cbind(temp,exp=exp,device=case,ET=ET)
#   alpha_beta=rbind(alpha_beta,temp)
#   
#   
#   ### Replace `database` with the test set
#   # database = test_data
#   # apollo_inputs = apollo_validateInputs()
#   # 
#   ### Cross-validation
#   # predictions = apollo_prediction(model, apollo_probabilities, apollo_inputs)
#   # 
#   # LL<-mean(log(predictions$chosen))
#   # 
#   # FPR<-sum(predictions[,6]==apply(predictions[,c(3:5)],1,max))/nrow(predictions)
#   # 
#   # Cross_validation<-rbind(Cross_validation,
#   #                         cbind(exp=exp,device=case,ET=ET,LL=LL,FPR=FPR))
#   # 
#   # 
#   # outOfSample_settings=list(nRep=5,
#   #                           validationSize=0.1)
#   # apollo_outOfSample(apollo_beta,
#   #                    apollo_fixed,
#   #                    apollo_probabilities,
#   #                    apollo_inputs,
#   #                    estimate_settings=list(maxIterations=400),outOfSample_settings)
#   # i=i+1
# }


################# now compare models ############

rm(ET)

VST_summary$ET[VST_summary$ET=='WG']='Web-Gazer'
VST_summary$ET[VST_summary$ET=='WL']='Eye-Link'
VST_summary<-VST_summary%>%filter(ET=='Eye-Link'|ET=='Web-Gazer')

VST_summary=VST_summary%>%group_by(device,ET)%>%mutate(mean_VST=mean(VST),
                                                       sd=sd(VST))

## VST distribution
VST_kstest <- VST_summary %>%
  group_by(device) %>%
  do({
    data = .
    combn(unique(data$ET), 2, function(x) {
      test <- ks.test(data$VST[data$ET == x[1]], data$VST[data$ET == x[2]])
      t_test<- t.test(data$VST[data$ET == x[1]], data$VST[data$ET == x[2]])
      mean_qst1<-mean(data$VST[data$ET == x[1]])
      mean_qst2<-mean(data$VST[data$ET == x[2]])
      se_qst1<-sd(data$VST[data$ET == x[1]])/sqrt(sum(data$ET == x[1]))
      se_qst2<-sd(data$VST[data$ET == x[2]])/sqrt(sum(data$ET == x[2]))
      ci1=paste0(round(mean_qst1,3)," (",round(mean_qst1,3)-round(1.96*se_qst1,3),", ",round(mean_qst1,3)+round(1.96*se_qst1,3),")")
      ci2=paste0(round(mean_qst2,3)," (",round(mean_qst2,3)-round(1.96*se_qst2,3),", ",round(mean_qst2,3)+round(1.96*se_qst2,3),")")
      
      overlap<-overlap(list(data$VST[data$ET == x[1]], data$VST[data$ET == x[2]]))
      data.frame(Group1 = x[1], Group2 = x[2], D.statistic = test$statistic, P.value = test$p.value,
                 overlap=overlap$OV,mean1=ci1,mean2=ci2,t_test=t_test$statistic,p_value_t=t_test$p.value
      )
    }, simplify = FALSE) %>% bind_rows()
  }) %>%
  bind_rows() %>%
  mutate(Comparison=paste0(Group1,"-",Group2),
         Label = paste("KS:",round(D.statistic, 3), "p:", round(P.value, 3)),
         Label_ov = paste0("OV ",Comparison,": ",round(overlap*100, 2),"%"))

write.csv(VST_kstest,'output/VST_ttest.csv')

desired_order2 <- c("monitor", "laptop", "tablet","mobile")  # Adjust as needed

# Reorder the 'attribute' within the ggplot call or before it
VST_kstest$Comparison=as.factor(VST_kstest$Comparison)
VST_kstest$Comparison=as.numeric(VST_kstest$Comparison)
VST_kstest$device<-factor(VST_kstest$device,levels=desired_order2)
VST_summary$device<-factor(VST_summary$device,levels=desired_order2)

# Define a color palette
palette <- brewer.pal(n = 6, name = "Dark2")


p1<-ggplot(VST_summary, aes(x=VST,fill=ET,color=ET)) +
  geom_density(alpha=0.3,position='identity') +
  geom_vline(aes(xintercept = mean_VST,color=ET), show.legend = FALSE) +
  facet_wrap(~device, ncol=2,scales="free") +
  scale_fill_manual(values = palette) +  # Apply the custom color palette
  scale_color_manual(values = palette) +# Define a color palette
  theme_minimal() +
  xlab('QWTP [SGD/qual]') + ylab('density') +
  guides(fill=guide_legend(title="ET"), color="none") +
  # geom_text(data=VST_kstest, aes(label=Label_ov, y=Inf, x=Inf,color=NULL,fill=NULL),
  #           vjust=VST_kstest$Comparison+(VST_kstest$Comparison-1)*1.5, hjust=1, size=2.5, color="black",
  #           show.legend = FALSE,
  #           show_guide  = FALSE)  +
  theme(legend.position = "bottom") 


VST_kstest <- VST_kstest %>%
  arrange(overlap,device) %>%
  mutate(Group = factor(paste0(Group1, "-", Group2), levels = unique(paste0(Group1, "-", Group2))))


p2 <- ggplot(VST_kstest,aes(x=device,y=overlap,fill=device))+
  geom_bar(alpha=0.8,stat="identity",position="dodge")+
  theme_minimal()+
  scale_fill_manual(values = palette) + 
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  xlab('Device')+ylab('% Overlapping')+
  theme(legend.position = "bottom") 


# Combine plots horizontally
combined_plot_horizontal <- p1 + p2 + plot_layout(ncol = 2)

# Or, combine plots vertically
combined_plot_vertical <- p1 / p2 + plot_layout(nrow = 2)

# Print one of the combined layouts
png(filename="QWTP_v1_transport.png",width=20,height=15,units="cm",res=300)
print(combined_plot_horizontal)
dev.off()

## VST - across devices
VST_kstest2 <- VST_summary %>%
  group_by(ET) %>%
  do({
    data = .
    combn(unique(data$device), 2, function(x) {
      test <- ks.test(data$VST[data$device == x[1]], data$VST[data$device == x[2]])
     
      overlap<-overlap(list(data$VST[data$device == x[1]], data$VST[data$device == x[2]]))
      data.frame(Group1 = x[1], Group2 = x[2], D.statistic = test$statistic, P.value = test$p.value,
                 overlap=overlap$OV
      )
    }, simplify = FALSE) %>% bind_rows()
  }) %>%
  bind_rows() %>%
  mutate(Comparison=paste0(Group1,"-",Group2),
         Label = paste("KS:",round(D.statistic, 3), "p:", round(P.value, 3)),
         Label_ov = paste0("OV ",Comparison,": ",round(overlap*100, 2),"%"))


VST_kstest2$Comparison=as.factor(VST_kstest2$Comparison)
VST_kstest2$Comparison=as.numeric(VST_kstest2$Comparison)

p1<-ggplot(VST_summary, aes(x=VST,fill=device,color=device)) +
  geom_density(alpha=0.4,position='identity') +
  #geom_vline(aes(xintercept = mean_VST,color=device), show.legend = FALSE) +
  facet_wrap(~ET, ncol=2,scales="free") +
  scale_fill_manual(values = palette) +  # Apply the custom color palette
  scale_color_manual(values = palette) +# Define a color palette
  theme_minimal() +
  xlab('QWTP [SGD/qual]') + ylab('density') +
  guides(fill=guide_legend(title="Device"), color="none") +
  # geom_text(data=VST_kstest2, aes(label=Label_ov, y=Inf, x=Inf,color=NULL,fill=NULL),
  #           vjust=VST_kstest2$Comparison+(VST_kstest2$Comparison-1)*1.5, hjust=1, size=2.5, color="black",
  #           show.legend = FALSE,
  #           show_guide  = FALSE)  +
  theme(legend.position = "bottom") 


VST_kstest2 <- VST_kstest2 %>%
  arrange(ET, overlap) %>%
  mutate(Group = factor(paste0(Group1, "-", Group2), levels = unique(paste0(Group1, "-", Group2))))


p2 <- ggplot(VST_kstest2,aes(x=Group,y=overlap,fill=ET))+
  geom_bar(alpha=0.8,stat="identity",position="dodge")+
  facet_wrap(~ET)+theme_minimal()+
  scale_fill_manual(values = palette) + 
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  xlab('Comparison')+ylab('% Overlapping')+
  theme(legend.position = "bottom") 


# Combine plots horizontally
combined_plot_horizontal <- p1 + p2 + plot_layout(ncol = 2)

# Or, combine plots vertically
combined_plot_vertical <- p1 / p2 + plot_layout(nrow = 2)

# Print one of the combined layouts
png(filename="QWTP_v2_transport.png",width=20,height=15,units="cm",res=300)
print(combined_plot_horizontal)
dev.off()



#### alpha_beta#### alpha_betscale_y_binned()a
# alpha_beta$ET[alpha_beta$ET=='WG']='Web-Gazer'
# alpha_beta$ET[alpha_beta$ET=='WL']='Eye-Link'
# 
# ggplot(alpha_beta,aes(x=var,y=alpha_mean,col=device,fill=device))+geom_point( position=position_dodge(width=0.25))+geom_errorbar(aes(ymin=alpha_mean-qt(0.95,1000)*alpha_se,ymax=alpha_mean+qt(0.95,1000)*alpha_se),position=position_dodge(width=0.25))+
#   facet_wrap(~ET)+theme_minimal()+geom_hline(yintercept = 0,linetype="dashed")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## IP distribution

summary_data1 <- IPs %>%
  group_by(device, ET) %>%
  reframe(
    Variable="time",
    mean_prob =(IP_alpha_time)
  )


summary_data2 <- IPs %>%
  group_by(device, ET) %>%
  reframe(
    Variable= "cost",
    mean_prob = (IP_alpha_cost)
    
  )
summary_data3 <- IPs %>%
  group_by(device, ET) %>%
  reframe(
    Variable= "comfort",
    mean_prob = (IP_alpha_comfort)
    
  )

summary_data<-rbind(summary_data1,summary_data2,summary_data3)
summary_data$ET[summary_data$ET=='WG']='Web-Gazer'
summary_data$ET[summary_data$ET=='WL']='Eye-Link'


# IPxalpha

# device x variable

Summary_kstest <- summary_data %>%
  group_by(Variable,device) %>%
  do({
    data = .
    combn(unique(data$ET), 2, function(x) {
      test <- ks.test(data$mean_prob[data$ET == x[1]], data$mean_prob[data$ET == x[2]])
      t_test<- t.test(data$mean_prob[data$ET == x[1]], data$mean_prob[data$ET == x[2]])
      mean_qst1<-mean(data$mean_prob[data$ET == x[1]])
      mean_qst2<-mean(data$mean_prob[data$ET == x[2]])
      se_qst1<-sd(data$mean_prob[data$ET == x[1]])/sqrt(sum(data$ET == x[1]))
      se_qst2<-sd(data$mean_prob[data$ET == x[2]])/sqrt(sum(data$ET == x[2]))
      ci1=paste0(round(mean_qst1,3)," (",round(mean_qst1,3)-round(1.96*se_qst1,3),", ",round(mean_qst1,3)+round(1.96*se_qst1,3),")")
      ci2=paste0(round(mean_qst2,3)," (",round(mean_qst2,3)-round(1.96*se_qst2,3),", ",round(mean_qst2,3)+round(1.96*se_qst2,3),")")
      
      overlap<-overlap(list(data$mean_prob[data$ET == x[1]], data$mean_prob[data$ET == x[2]]))
      data.frame(Group1 = x[1], Group2 = x[2], D.statistic = test$statistic, P.value = test$p.value,
                 overlap=overlap$OV,mean1=ci1,mean2=ci2,t_test=t_test$statistic,p_value_t=t_test$p.value)
    }, simplify = FALSE) %>% bind_rows()
  }) %>%
  bind_rows() %>%
  mutate(Comparison=paste0(Group1,"-",Group2),
         Label = paste("KS:", round(D.statistic, 3), "p:", round(P.value, 3)),
         Label_ov = paste0("OV ",Comparison,": ", round(overlap*100, 2),"%"))

write.csv(Summary_kstest,'output/IP_test_transport.csv')
Summary_kstest$Comparison=as.factor(Summary_kstest$Comparison)
Summary_kstest$Comparison=as.numeric(Summary_kstest$Comparison)

desired_order2 <- c("monitor", "laptop", "tablet","mobile")  # Adjust as needed

# Reorder the 'attribute' within the ggplot call or before it
Summary_kstest$device<-factor(Summary_kstest$device,levels=desired_order2)
summary_data$device<-factor(summary_data$device,levels=desired_order2)

# Create the plot
p1<-ggplot(summary_data, aes(x=mean_prob, fill=ET, color=ET)) +
  geom_density(alpha=0.6, position='identity') +  # Ensure binwidth is set
  scale_fill_manual(values = palette) +  
  scale_color_manual(values = palette) +
  facet_wrap(~device+Variable, ncol=3,scales="free") +
  theme_minimal() +
  xlab(TeX("Posterior mean $IP_n\\cdot \\alpha_k$")) +
  ylab('Density') +theme(legend.position = "bottom")+
  theme(axis.text.x=element_text(size=7,
                                 angle=45,hjust=1),
        )
# geom_text(data=Summary_kstest, aes(label=Label_ov, y=Inf, x=Inf, color=NULL, fill=NULL),
#           vjust=Summary_kstest$Comparison+(Summary_kstest$Comparison-1)*1.5, hjust=1, size=2.5, color="black", show.legend = FALSE) +
# theme(legend.position = "right",  # Adjust legend positioning
#       legend.title = element_text(size = 10),  # Adjust legend title size
#       legend.text = element_text(size = 8)) 

Summary_kstest <- Summary_kstest %>%
  arrange(Variable, overlap) %>%
  mutate(Group = factor(device, levels = unique(device)))

p2 <- ggplot(Summary_kstest,aes(x=device,y=overlap,fill=Variable))+
  geom_bar(alpha=0.7,stat="identity",position="dodge")+
  scale_fill_manual(values = palette) + 
  ylab('% Overlapping')+
  #facet_wrap(~device,ncol=2)+theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  xlab('Device')+
  theme_minimal()+theme(legend.position = "bottom")

# Combine plots horizontally
combined_plot_horizontal <- p1 + p2 + plot_layout(ncol = 2)

# Or, combine plots vertically
combined_plot_vertical <- p1 / p2 + plot_layout(nrow = 2)

# Print one of the combined layouts
png(filename="IP_alpha_v1_transport.png",width=20,height=15,units="cm",res=300)
print(combined_plot_horizontal)
dev.off()

### variable + ET grouping

Summary_kstest2 <- summary_data %>%
  group_by(Variable,ET) %>%
  do({
    data = .
    combn(unique(data$device), 2, function(x) {
      test <- ks.test(data$mean_prob[data$device == x[1]], data$mean_prob[data$device == x[2]])
      overlap<-overlap(list(data$mean_prob[data$device == x[1]], data$mean_prob[data$device == x[2]]))
      data.frame(Group1 = x[1], Group2 = x[2], D.statistic = test$statistic, P.value = test$p.value,
                 overlap=overlap$OV)
    }, simplify = FALSE) %>% bind_rows()
  }) %>%
  bind_rows() %>%
  mutate(Comparison=paste0(Group1,"-",Group2),
         Label = paste("KS:", round(D.statistic, 3), "p:", round(P.value, 3)),
         Label_ov = paste0("OV ",Comparison,": ", round(overlap*100, 2),"%"))

Summary_kstest2$Comparison=as.factor(Summary_kstest2$Comparison)
Summary_kstest2$Comparison=as.numeric(Summary_kstest2$Comparison)

p<-ggplot(summary_data, aes(x=mean_prob, fill=device, color=device)) +
  geom_density(alpha=0.6, position='identity') +  # Ensure binwidth is set
  scale_fill_manual(values = palette) +  
  scale_color_manual(values = palette) +
  facet_wrap(~ET+Variable, ncol=3,scales="free") +
  theme_minimal() +
  xlab(TeX("Posterior mean $IP_n\\cdot \\alpha_k$")) +
  ylab('Density') +theme(legend.position = "bottom")+
  theme(axis.text.x=element_text(size=7,
                                 angle=45,hjust=1))
# geom_text(data=Summary_kstest2, aes(label=Label_ov, y=Inf, x=Inf, color=NULL, fill=NULL),
#           vjust=Summary_kstest2$Comparison+(Summary_kstest2$Comparison-1)*1.5, hjust=1, size=2.5, color="black", show.legend = FALSE) +
# theme(legend.position = "right",  # Adjust legend positioning
#       legend.title = element_text(size = 10),  # Adjust legend title size
#       legend.text = element_text(size = 8)) 

Summary_kstest2$facet=as.numeric(as.factor(paste(Summary_kstest2$ET,Summary_kstest2$Variable,sep="-")))

Summary_kstest2 <- Summary_kstest2 %>%
  arrange(ET, Variable, overlap) %>%
  mutate(Group = factor(paste0(Group1, "-", Group2), levels = unique(paste0(Group1, "-", Group2))))

p2 <- ggplot(Summary_kstest2,aes(x=Group,y=overlap,fill=Variable))+
  geom_bar(alpha=0.7,stat="identity",position="dodge")+
  facet_wrap(~ET)+theme_minimal()+
  scale_fill_manual(values = palette) + 
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  xlab('Comparison')+ylab('% Overlapping')+theme(legend.position = "bottom")

# Combine plots horizontally
combined_plot_horizontal <- p + p2 + plot_layout(ncol = 2)

# Or, combine plots vertically
combined_plot_vertical <- p / p2 + plot_layout(nrow = 2)

# Print one of the combined layouts
png(filename="IP_alpha_v2_transport.png",width=20,height=15,units="cm",res=300)
print(combined_plot_horizontal)
dev.off()
# 
# #### log-likelihood plot ######
# LL_list<-data.frame(LL_list)
# LL_list$LLfinal<-as.numeric(LL_list$LLfinal)
# write.csv(LL_list,'LL_counts.csv')
# 
# LL_list_baseline<-LL_list%>%filter(ET=="Baseline")
# LL_list_model<-LL_list%>%filter(ET!="Baseline")
# LL_list_model<-LL_list_model%>%group_by(ET,device)%>%mutate(
#   LRT=-2*(LL_list_baseline$LLfinal[LL_list_baseline$device==device]-
#             LLfinal)
# )
# 
# ggplot(LL_list_model, aes(x = ET, y = LRT)) +
#   geom_bar(stat = 'identity') +
#   facet_wrap(~ device,scale="free") +theme_minimal()
# 
# ##### Cross Validation ####
# Cross_validation2=as.data.frame(Cross_validation)
# 
# Cross_validation2=Cross_validation2%>%mutate(LL=as.numeric(LL),
#                                            FPR=as.numeric(FPR))
# 
# ggplot(Cross_validation2, aes(x = ET, y = exp(LL))) +
#   geom_bar(stat = 'identity') +
#   facet_wrap(~ device,scale="free") +theme_minimal()
# 
# ggplot(Cross_validation2, aes(x = ET, y = FPR)) +
#   geom_bar(stat = 'identity') +
#   facet_wrap(~ device,scale="free") +theme_minimal()

# 
# cases=c('monitor','laptop','tablet','mobile')
# exp='exp1'
# ETs=c('WG','WL')
# 
# modelNames=
#   c(apply(expand.grid(exp,cases,'WL','MixedLogit'),1,function(x){paste0(x,collapse = "_")}),   
#     apply(expand.grid(exp,cases,ETs,'LIP_counts'),1,function(x){paste0(x,collapse = "_")}))
# 
# combineResults_settings=list(estimateDigits=3,
#                              modelNames=modelNames
# )
# 
# 
# apollo_combineResults(combineResults_settings)


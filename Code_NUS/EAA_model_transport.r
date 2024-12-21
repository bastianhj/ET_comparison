# EAA model (latent classes without ET data)
# Eye-tracking comparison experiment
# Universidad de Chile - National University of Singapore
# Bastian Henríquez-Jara
library(apollo)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cases=c('monitor','laptop','tablet','mobile')
exp='EXPERIMENT1'



########3

for (case in cases){
  
  path=paste0('../Compiled_data/',exp,'_',case,'_','choices','.csv')
  data=read.csv(path) 
  data$Choice_id=1*(data$Choice=="Bus_label")+2*(data$Choice=="Metro_label")+3*(data$Choice=="RH_label")
  #robse=c()
  apollo_initialise()
  #for(i in 101:150){
  ### Initialise code
  
  
  ### Set core controls
  apollo_control = list(
    modelName       = paste0(exp,'_',case,'_EAA'),
    modelDescr      = "MNL for transport choices - ET comparison",
    indivID         = "pid", 
    outputDirectory = "output"
  )
  
  
  database = data
  
  apollo_beta = c(b_bus=0,
                  b_metro=0,
                  b_RH=0,
                  b_time=0,
                  b_cost=0,
                  b_comfort=0,
                  b_time_att=0,
                  b_cost_att=0,
                  b_comfort_att=0
  )
  apollo_fixed=c('b_bus')
  
  apollo_lcPars=function(apollo_beta, apollo_inputs){
    lcpars = list()
    # lcpars[["b_brand1"]] = list(b_brand1)
    # lcpars[["b_brand2"]] = list(b_brand2)
    # lcpars[["b_brand3"]] = list(b_brand3)
    # lcpars[["b_material1"]] = list(b_material1)
    # lcpars[["b_material2"]] = list(b_material2)
    # lcpars[["b_system1"]] = list(b_system1)
    # lcpars[["b_design1"]] = list(b_design1)
    # lcpars[["b_design2"]] = list(b_design2)
    # lcpars[["b_design3"]] = list(b_design3)
    # lcpars[["b_pricecup"]] = list(b_pricecup)
    # lcpars[["b_price"]] = list(b_price)
    ### Utilities of class allocation model
    
    V_cost=b_cost_att#+b_va*va1
    V_time=b_time_att#+b_va*va2
    V_comfort=b_comfort_att#+b_va*va3
    
    V_attributes=cbind(V_cost,V_time,V_comfort)
    class_matrix=as.matrix(expand.grid(c(1,0), c(1,0), c(1,0)))
    
    pi_attribute=c(1/(1+exp(-V_attributes)))
    list_pi=list()
    for(i in 1:8){
      A=prod((pi_attribute^class_matrix[i,])*(1-pi_attribute)^(1-class_matrix[i,]))
      list_pi[[i]]=A#apply(A,1,prod)#apply(matrix_power(A,(class_matrix[i,])),1,prod)
    }
    
    
    lcpars[["pi_values"]] = list_pi
    
    return(lcpars)
  }
  
  apollo_inputs = apollo_validateInputs()
  
  apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
    
    ### Attach inputs and detach after function exit
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    
    ### Create list of probabilities P
    P = list()
    
    mnl_settings = list(
      alternatives  = c(bus=1,metro=2,RH=3),
      avail=1,
      choiceVar = Choice_id
    )

    class_matrix=(expand.grid(c(1,0), c(1,0), c(1,0)))
    
    for(s in 1:8){
      
      ### Compute class-specific utilities
      V=list()
      V[['bus']]=b_bus+b_cost*Bus_cost*class_matrix[s,1]+b_time*Bus_travel_time*class_matrix[s,2]+b_comfort*Bus_Comfort*class_matrix[s,3]
      
      V[['metro']]=b_metro+b_cost*metro_cost*class_matrix[s,1]+b_time*metro_travel_time*class_matrix[s,2]+b_comfort*metro_Comfort*class_matrix[s,3]
      
      V[['RH']]=b_RH+b_cost*RH_cost*class_matrix[s,1]+b_time*RH_travel_time*class_matrix[s,2]+b_comfort*RH_Comfort*class_matrix[s,3]
      
      mnl_settings$utilities     = V
      mnl_settings$componentName = paste0("Class_",s)
      
      ### Compute within-class choice probabilities using MNL model
      P[[paste0("Class_",s)]] = apollo_mnl(mnl_settings, functionality)
      
      ### Take product across observation for same individual
      P[[paste0("Class_",s)]] = apollo_panelProd(P[[paste0("Class_",s)]], apollo_inputs ,functionality)
      
    }
    
    lc_settings  = list(inClassProb = P, classProb=pi_values)
    
    P[["model"]] = apollo_lc(lc_settings, apollo_inputs, functionality)
    
    #P = apollo_combineModels(P, apollo_inputs, functionality)
    ### Take product across observation for same individual
    #P = apollo_avgIntraDraws(P, apollo_inputs, functionality)
    #P = apollo_panelProd(P, apollo_inputs, functionality)
    ### Average across inter-individual draws
    # P = apollo_avgInterDraws(P, apollo_inputs, functionality)
    ### Prepare and return outputs of function
    P = apollo_prepareProb(P, apollo_inputs, functionality)
    
    return(P)
  }
  
  
  
  # ################################################################# #
  #### MODEL ESTIMATION                                            ####
  # ################################################################# #
  
  model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
  
  
  apollo_modelOutput(model)
  
  # ----------------------------------------------------------------- #
  #---- FORMATTED OUTPUT (TO FILE, using model name)               ----
  # ----------------------------------------------------------------- #
  # betas<-rbind(betas,model$estimate)
  # robse<-rbind(robse,model$robse)
  # }
  
  apollo_saveOutput(model)
  
}




### now compare models

combineResults_settings=list(estimateDigits=3,
                             modelNames=paste0(exp,'_',cases,'_EAA'))

apollo_combineResults(combineResults_settings)

?apollo_combineResults

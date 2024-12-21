# pizzas baseline model
# Eye-tracking comparison experiment
# Universidad de Chile - National University of Singapore
# Bastian Henríquez-Jara
library(apollo)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cases=c('monitor','laptop','tablet','mobile')
exp='EXPERIMENT2'

for (case in cases){
  
  path=paste0('../Compiled_data/',exp,'_',case,'_','choices','.csv')
  data=read.csv(path) 
  data$Choice_id=1*(data$Choice=="Option1_label")+2*(data$Choice=="Option2_label")
  #robse=c()
  apollo_initialise()
  #for(i in 101:150){
  ### Initialise code
  
  
  ### Set core controls
  apollo_control = list(
    modelName       = paste0(exp,'_',case,'_baseline_MNL'),
    modelDescr      = "MNL for pizzass choices - ET comparison",
    indivID         = "pid", 
    outputDirectory = "output"
  )
  
  database = data
  
  apollo_beta = c(b_opt1=0,
                  b_opt2=0,
                  b_qual=0,
                  b_cost=0)
  
  apollo_fixed=c('b_opt2')
  
  
  apollo_inputs = apollo_validateInputs()
  
  apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
    
    ### Attach inputs and detach after function exit
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    
    ### Create list of probabilities P
    P = list()
    V = list()
    
    
    V[['opt1']]=b_opt1+b_cost*cost1+b_qual*qual1
    
    V[['opt2']]=b_opt2+b_cost*cost2+b_qual*qual2
  
    mnl_settings = list(
      alternatives  = c(opt1=1,opt2=2),
      avail=1,
      choiceVar = Choice_id,
      utilities = V
    )
    
    P[["model"]] = apollo_mnl(mnl_settings, functionality)
    
    #P = apollo_combineModels(P, apollo_inputs, functionality)
    ### Take product across observation for same individual
    #P = apollo_avgIntraDraws(P, apollo_inputs, functionality)
    P = apollo_panelProd(P, apollo_inputs, functionality)
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
                             modelNames=paste0(exp,'_',cases,'_baseline_MNL'))

apollo_combineResults(combineResults_settings)

?apollo_combineResults

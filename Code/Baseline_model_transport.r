# transport base model
# Eye-tracking comparison experiment
# Universidad de Chile - National University of Singapore
# Bastian Henríquez-Jara
library(apollo)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cases=c('monitor','laptop','tablet','mobile')
exp='EXPERIMENT1'
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
  modelName       = paste0(exp,'_',case,'_baseline_MNL'),
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
                b_comfort=0)

apollo_fixed=c('b_bus')


apollo_inputs = apollo_validateInputs()

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  V = list()

  
  V[['bus']]=b_bus+b_cost*Bus_cost+b_time*Bus_travel_time+b_comfort*Bus_Comfort
  
  V[['metro']]=b_metro+b_cost*metro_cost+b_time*metro_travel_time+b_comfort*metro_Comfort
  
  V[['RH']]=b_RH+b_cost*RH_cost+b_time*RH_travel_time+b_comfort*RH_Comfort

  
  mnl_settings = list(
    alternatives  = c(bus=1,metro=2,RH=3),
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

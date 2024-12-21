# compile data for each case
# Eye-tracking comparison experiment
# Universidad de Chile - National University of Singapore
# Bastian Henríquez-Jara
rm(list = ls())
library(dplyr)
library(plyr) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## change this if needed
cases=c('monitor','laptop','tablet','mobile')
N_participants=40

exp='EXPERIMENT2'

filenames <- list.files(paste0('../Data_NUS/',exp,'/'),pattern="*.csv", full.names=TRUE)


### function to detect if file is empty
is_csv_empty <- function(file_path) {
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop("File does not exist.")
  }
  
  # Check if the file size is zero
  file_size <- file.info(file_path)$size
  if (file_size == 0) {
    return(TRUE) # File is empty
  }
  
  # Open a connection to the file and check its content
  con <- file(file_path, "r")
  first_line <- readLines(con, n = 1, warn = FALSE)
  close(con)
  
  # Check if the first line exists (non-empty header)
  if (length(first_line) == 0 || all(trimws(first_line) == "")) {
    return(TRUE) # File is effectively empty
  }
  
  # Attempt to read a few rows to check for actual data
  tryCatch({
    data <- read.csv(file_path, nrows = 5)
    if (nrow(data) == 0) {
      return(TRUE) # No data rows found
    }
    return(FALSE) # File has data
  }, error = function(e) {
    warning("File might have structural issues: ", e$message)
    return(TRUE) # Assume empty or problematic file
  })
}
########## compile data ###############
errors_id=c()
for (case in cases){
  

devices<-sub(".*_", "",filenames)
ids<-sub(".*/", "",filenames)
ids<-sub("_.*", "",ids)[devices==paste0(case,'.csv')]

temp<-sub(paste0("_",case,'.*'), "", filenames)
temp=sub(".*_", "", temp)

block=temp[devices==paste0(case,'.csv')]


compiled=c()
#compile data for all participants
for(id in 1:N_participants){
 # print(length(block[ids==paste0('pct',id)]))
  if (length(block[ids==paste0('pct',id)])==1){
    path=paste0('../Data_NUS/',exp,'/pct',id,'_',block[ids==paste0('pct',id)],'_',case,'.csv')
    if(!is_csv_empty(path)){
      data=read.csv(path)  
      data=data[!is.na(data$Trial),]
      data$pid=id
      compiled=rbind.fill(compiled,data)
    }else{
      errors_id=c(errors_id,id)
    }
  }else{
    errors_id=c(errors_id,id)
  }
}

if(exp=='EXPERIMENT1'){
  compiled_choices<-compiled[,c('pid',"Trial","TrialDuration","Block","Device","RH_cost" ,"RH_travel_time","RH_Comfort","metro_cost","metro_travel_time","metro_Comfort","Bus_cost","Bus_travel_time","Bus_Comfort","Choice")]
}else{
  compiled_choices<-compiled[,c('pid',"Trial","TrialDuration","Block","Device","cost1","qual1","cost2","qual2","Choice")]
}

write.csv(compiled_choices,paste0('../Compiled_data/',exp,'_',case,'_choices.csv'))

}

errors_id



##################################################################################################
# MultiLabel Similarities Measures
# Copyright (C) 2021                                                                             #
#                                                                                                #
# This code is free software: you can redistribute it and/or modify it under the terms of the    #
# GNU General Public License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version. This code is distributed in the hope       #
# that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    #
# more details.                                                                                  #
#                                                                                                #
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin                     #
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus Sao Carlos           #
# Computer Department (DC: https://site.dc.ufscar.br/)                                           #
# Program of Post Graduation in Computer Science (PPG-CC: http://ppgcc.dc.ufscar.br/)            #
# Bioinformatics and Machine Learning Group (BIOMAL: http://www.biomal.ufscar.br/)               #
#                                                                                                #
##################################################################################################

rm(list=ls())

##################################################################################################
# Configures the workspace according to the operating system                                     #
##################################################################################################
sistema = c(Sys.info())
FolderRoot = ""
if (sistema[1] == "Linux") {
  FolderRoot = paste("/home/", sistema[7], "/MultiLabelSimilaritiesMeasures", sep = "")
  setwd(FolderRoot)
} else {
  FolderRoot = paste("C:/Users/", sistema[7], "/MultiLabelSimilaritiesMeasures", sep = "")
  setwd(FolderRoot)
}
setwd(FolderRoot)
FolderScripts =  paste(FolderRoot, "/R", sep = "")


##################################################################################################
# Options Configuration                                                                          #
##################################################################################################
cat("\nopções de configurações")
options(java.parameters = "-Xmx32g")
options(show.error.messages = TRUE)
options(scipen = 10)


##################################################################################################
# Read the dataset file with the information for each dataset                                    #
##################################################################################################
setwd(FolderRoot)
datasets <- data.frame(read.csv("datasets-2022.csv"))


##################################################################################################
# ARGS COMMAND LINE                                                                              #
##################################################################################################
args <- commandArgs(TRUE)


##################################################################################################
# Get dataset information                                                                        #
##################################################################################################
number_dataset = as.numeric(args[1])
cat("\nMLSM: DS \t ", number_dataset)


##################################################################################################
# Get dataset information                                                                        #
##################################################################################################
ds <- datasets[number_dataset,]


##################################################################################################
# Get the number of cores                                                                        #
##################################################################################################
number_cores <- as.numeric(args[2])
cat("\nMLSM: cores \t ", number_cores)


##################################################################################################
# Get the number of folds                                                                        #
##################################################################################################
number_folds <- as.numeric(args[3])
cat("\nMLSM: folds \t ", number_folds)


##################################################################################################
# Get the number of folds                                                                        #
##################################################################################################
FolderResults  <- toString(args[4])
cat("\nMLSM: folder \t ", FolderResults)


##################################################################################################
# Get dataset name                                                                               #
##################################################################################################
dataset_name  <- toString(ds$Name) 
cat("\nMLSM: nome \t ", dataset_name)


##################################################################################################
# DON'T RUN -- it's only for test the code
 ds <- datasets[23,]
 dataset_name = ds$Name
 number_dataset = ds$Id
 number_cores = 10             
 number_folds = 10
 FolderResults = "/dev/shm/teste"
##################################################################################################


##################################################################################################
# CONFIG THE FOLDER RESULTS                                                                      #
##################################################################################################
if(dir.exists(FolderResults)==FALSE){
  dir.create(FolderResults)
  cat("\n")
}


##################################################################################################
# LOAD SOURCES                                                                                     #
##################################################################################################
setwd(FolderScripts)
source("runCV.R") 

setwd(FolderScripts)
source("runNCV.R") 

setwd(FolderScripts)
source("libraries.R") 

setwd(FolderScripts)
source("utils.R") 


##################################################################################################
# GET THE DIRECTORIES                                                                            #
##################################################################################################
setwd(FolderRoot)
cat("\nGet directories\n")
folder = diretorios(dataset_name, FolderResults)


##################################################################################################
cat("\nCopy FROM google drive \n")
destino = folder$FolderDS
origem = paste("cloud:Datasets/CrossValidation_WithValidation/", dataset_name, sep="")
comando = paste("rclone -v copy ", origem, " ", destino, sep="")
cat("\n", comando, "\n") 
a = print(system(comando))
a = as.numeric(a)
if(a != 0) {
  stop("Erro RCLONE")
  quit("yes")
}
 
##################################################################################################
cat("\nCopy FROM google drive \n")
destino = folder$FolderDS
origem = paste("cloud:Datasets/Originais/", dataset_name, ".arff", sep="")
comando = paste("rclone -v copy ", origem, " ", destino, sep="")
cat("\n", comando, "\n") 
a = print(system(comando))
a = as.numeric(a)
if(a != 0) {
 stop("Erro RCLONE")
 quit("yes")
}
 
##################################################################################################
cat("\nCopy FROM google drive \n")
destino = folder$FolderDS
origem = paste("cloud:Datasets/Originais/", dataset_name, ".xml", sep="")
comando = paste("rclone -v copy ", origem, " ", destino, sep="")
cat("\n", comando, "\n") 
a = print(system(comando))
a = as.numeric(a)
if(a != 0) {
 stop("Erro RCLONE")
 quit("yes")
}


##################################################################################################
# execute the code and get the total execution time                                              #
# n_dataset, number_cores, number_folds                                                          #
##################################################################################################

if(number_folds==1){
  cat("\nExecute MLSM without Cross Validation \n")
  timeFinal <- system.time(results <- executeMLSM_NCV(number_dataset, FolderResults))
  print(timeFinal)
  
} else {
  cat("\nExecute MLSM with Cross Validation (number folds > 1) \n")
  timeFinal <- system.time(results <- executeMLSM_CV(ds, number_dataset, number_cores, number_folds, FolderResults))
  print(timeFinal)
}


 
###################################################################################################
cat("\nCopy to google drive \n")
origem = folder$FolderRD
destino = paste("cloud:[2022]ResultadosExperimentos/Similaridades/", dataset_name, sep="")
comando = paste("rclone -v copy ", origem, " ", destino, sep="")
cat("\n", comando, "\n") 
a = print(system(comando))
a = as.numeric(a)
if(a != 0) {
  stop("Erro RCLONE")
  quit("yes")
}


##################################################################################################
# del                                                                                      #
##################################################################################################
cat("\nDelete folder \n")
str5 = paste("rm -r ", folder$FolderResults, sep="")
print(system(str5))

cat("\nDelete folder \n")
str5 = paste("rm -r ", FolderRoot, "/datasets/", dataset_name, sep="")
print(system(str5))

cat("\nDelete folder \n")
str5 = paste("rm -r ", folder$FolderRD, sep="")
print(system(str5))

gc()

rm(list=ls())

##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
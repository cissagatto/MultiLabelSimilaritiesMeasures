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

cat("\n\n###################################################################")
cat("\n# ====> MLSM:  SET WORK SPACE                                       #")
cat("\n#####################################################################\n\n")
FolderRoot = "~/MultiLabelSimilaritiesMeasures"
FolderScripts = paste(FolderRoot, "/R/", sep="")


cat("\n\n###################################################################")
cat("\n# ====> MLSM: OPTIONS CONFIGURATIONS                                #")
cat("\n#####################################################################\n\n")
options(java.parameters = "-Xmx32g")
options(show.error.messages = TRUE)
options(scipen=20)


cat("\n\n###################################################################")
cat("\n# ====> MLSM: READ DATASETS-2022                                    #")
cat("\n#####################################################################\n\n")
setwd(FolderRoot)
datasets <- data.frame(read.csv("datasets-2022.csv"))


cat("\n\n###################################################################")
cat("\n# ====> MLSM: GET THE ARGUMENTS COMMAND LINE                        #")
cat("\n#####################################################################\n\n")
args <- commandArgs(TRUE)


cat("\n\n###################################################################")
cat("\n# ====> MLSM: DATASET INFORMATION                                   #")
cat("\n#####################################################################\n\n")
ds <- datasets[args[1],]
print(ds)


number_dataset <- as.numeric(args[1])
cat("\n\n###################################################################")
cat("\n# ====> MLSM: NUMBER DATASET: ", number_dataset, "                  #")
cat("\n#####################################################################\n\n")


number_cores <- as.numeric(args[2])
cat("\n\n###################################################################")
cat("\n# ====> MLSM: GET THE NUMBER CORES: ", number_cores, "              #")
cat("\n#####################################################################\n\n")


number_folds <- as.numeric(args[3])
cat("\n\n###################################################################")
cat("\n# ====> MLSM: GET THE NUMBER FOLDS: ", number_folds, "              #")
cat("\n#####################################################################\n\n")


FolderResults <- toString(args[4])
cat("\n\n###################################################################")
cat("\n# ====> MLSM: GET THE STRING FOLDER RESULTS: ", FolderResults, "    #")
cat("\n#####################################################################\n\n")


dataset_name <- toString(ds$Name)
cat("\n\n###################################################################")
cat("\n# ====> MLSM: GET THE DATASET NAME: ", dataset_name, "              #")
cat("\n#####################################################################\n\n")


#ds <- datasets[22,]
#number_dataset = 22
#number_cores = 10
#number_folds = 10
#FolderResults = "/dev/shm/res"
#dataset_name = ds$Name



cat("\n\n###################################################################")
cat("\n# ====> MLSM: CREATING FOLDER RESULTS TEMP                          #")
cat("\n#####################################################################\n\n")
if(dir.exists(FolderResults)==FALSE){ dir.create(FolderResults)}
cat("\n")


cat("\n\n###################################################################")
cat("\n# ====> MLSM: LOAD SOURCES R                                        #")
cat("\n#####################################################################\n\n")

setwd(FolderScripts)
source("runCV.R") 

#setwd(FolderScripts)
#source("runNCV.R") 

setwd(FolderScripts)
source("libraries.R") 

setwd(FolderScripts)
source("utils.R") 


cat("\n\n###################################################################")
cat("\n# ====> MLSM: CREATING DIRECTORIES                                  #")
cat("\n#####################################################################\n\n")
setwd(FolderRoot)
cat("\nGet directories\n")
folder = diretorios(dataset_name, FolderResults)



#cat("\n\n###################################################################")
#cat("\n# ====> MLSM: COPY DATA SETS FROM GOOGLE DRIVE                      #")
#cat("\n#####################################################################\n\n")
#destino = folder$FolderDS
#origem = paste("cloud:Datasets/CrossValidation_WithValidation/", dataset_name, sep="")
#comando = paste("rclone -P copy ", origem, " ", destino, sep="")
#cat("\n\n\n", comando, "\n\n\n") 
#a = print(system(comando))
#a = as.numeric(a)
#if(a != 0) {
#  stop("Erro RCLONE")
#  quit("yes")
#}

 
#cat("\n\n###################################################################")
#cat("\n# ====> MLSM: COPY DATA SETS FROM GOOGLE DRIVE                      #")
#cat("\n#####################################################################\n\n")
#destino = folder$FolderDS
#origem = paste("cloud:Datasets/Originais/", dataset_name, ".arff", sep="")
#comando = paste("rclone -P copy ", origem, " ", destino, sep="")
#cat("\n\n\n", comando, "\n\n\n") 
#a = print(system(comando))
#a = as.numeric(a)
#if(a != 0) {
# stop("Erro RCLONE")
# quit("yes")
#}



cat("\n\n###################################################################")
cat("\n# ====> TCP-KNN-H: ORGANIZANDO AS COISAS                            #")
cat("\n#####################################################################\n\n")

cat("\nCOPIANDO DATASETS")
str26 = paste("cp ~/MultiLabelSimilaritiesMeasures/Datasets/", ds$Name, ".tar.gz ",
              folder$FolderDatasets, sep="")
res=system(str26)
if(res!=0){break}else{cat("\ncopiou")}


cat("\nDESCOMPACTANDO DATASETS")
str27 = paste("tar xzf ", folder$FolderResults,
              "/Datasets/", ds$Name, ".tar.gz -C ",
              folder$FolderResults, "/Datasets", sep="")
res=system(str27)
if(res!=0){break}else{cat("\ndescompactou")}


cat("\n APAGANDO TAR")
str28 = paste("rm ", folder$FolderResults, "/Datasets/",
              ds$Name, ".tar.gz", sep="")
res=system(str28)
if(res!=0){break}else{cat("\napagou")}



cat("\n\n###################################################################")
cat("\n# ====> MLSM: START                                                 #")
cat("\n#####################################################################\n\n")

if(number_folds==1){
  #cat("\nExecute MLSM without Cross Validation \n")
  #timeFinal <- system.time(results <- executeMLSM_NCV(number_dataset, FolderResults))
  #print(timeFinal)
  cat("\nnot available for now")
  
} else {
  cat("\nExecute MLSM with Cross Validation (number folds > 1) \n")
  timeFinal <- system.time(results <- executeMLSM_CV(ds, number_dataset, number_cores, number_folds, FolderResults))
  print(timeFinal)
  gc()
  
  cat("\n\n###################################################################")
  cat("\n# ====> MLSM: SAVE RUNTIME                                          #")
  cat("\n#####################################################################\n\n")
  result_set <- t(data.matrix(timeFinal))
  setwd(folder$FolderRD)
  write.csv(result_set, "Runtime.csv")
  print(timeFinal)
  cat("\n")
  
}

cat("\n\n###################################################################")
cat("\n# ====> MLSM: DELETE DATASETS FOLDER                                #")
cat("\n#####################################################################\n\n")
print(system(paste("rm -r ", folder$FolderDS, sep="")))
print(system(paste("rm -r ", folder$FolderDatasets, sep="")))

 
#cat("\n\n###################################################################")
#cat("\n# ====> MLSM: COPY TO GOOGLE DRIVE                                  #")
#cat("\n#####################################################################\n\n")
#origem = folder$FolderRD
#destino = paste("cloud:[2022]ResultadosExperimentos/Similaridades/", dataset_name, sep="")
#comando = paste("rclone -P copy ", origem, " ", destino, sep="")
#cat("\n\n\n", comando, "\n\n\n") 
#a = print(system(comando))
#a = as.numeric(a)
#if(a != 0) {
#  stop("Erro RCLONE")
#  quit("yes")
#}


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-NH: COPY TO HOME                                     #")
cat("\n#####################################################################\n\n")

str0 = "~/MultiLabelSimilaritiesMeasures/Reports"
if(dir.exists(str0)==FALSE){dir.create(str0)}

str2 = paste("cp -r ", folder$FolderResults, "/", ds$Name, "/ ", str0, sep="")
print(system(str2))


cat("\n\n###################################################################")
cat("\n# ====> MLSM: CLEAN                                                 #")
cat("\n#####################################################################\n\n")
str5 = paste("rm -r ", folder$FolderResults, sep="")
print(system(str5))
gc()
rm(list=ls())

##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
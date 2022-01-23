##################################################################################################

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

##################################################################################################
# Configures the workspace according to the operating system                                     #
##################################################################################################
sistema = c(Sys.info())
FolderRoot = ""
if (sistema[1] == "Linux"){
  FolderRoot = paste("/home/", sistema[7], "/MultiLabelSimilaritiesMeasures", sep="")
  setwd(FolderRoot)
} else {
  FolderRoot = paste("C:/Users/", sistema[7], "/MultiLabelSimilaritiesMeasures", sep="")
  setwd(FolderRoot)
}
setwd(FolderRoot)
FolderScripts = paste(FolderRoot, "/Scripts", sep="")

library(stringr)

setwd(FolderRoot)
datasets = data.frame(read.csv("datasets-2022.csv"))
n = nrow(datasets)

FolderJob = paste(FolderRoot, "/Jobs-Euler-2", sep="")
if(dir.exists(FolderJob)==FALSE){dir.create(FolderJob)}

i = 1
while(i<=n){
  
  dataset = datasets[i,]
  cat("\n\tDataset:", dataset$Name)
  
  nome = paste(FolderJob, "/l2-", dataset$Name, ".sh", sep="")
  output.file <- file(nome, "wb")
  
  write("#!/bin/bash", file = output.file, append = TRUE)
  
  write("eval \"$(conda shell.bash hook)\" ", file = output.file, append = TRUE)
  
  write("conda activate AmbienteTeste", file = output.file, append = TRUE)
  
  str1 = paste("rclone -v copy nuvem:Datasets/CrossValidation_WithValidation/", dataset$Name, " /mnt/nfs/home/elaine/MultiLabelSimilaritiesMeasures/Datasets/", dataset$Name, sep="")
  write(str1, file = output.file, append = TRUE)
  
  str0 = paste("rclone -v copy nuvem:Datasets/Originais/", dataset$Name, ".arff /mnt/nfs/home/elaine/MultiLabelSimilaritiesMeasures/Datasets/", dataset$Name, sep="")
  write(str0, file = output.file, append = TRUE)
  
  str10 = paste("rclone -v copy nuvem:Datasets/Originais/", dataset$Name, ".xml /mnt/nfs/home/elaine/MultiLabelSimilaritiesMeasures/Datasets/", dataset$Name, sep="")
  write(str10, file = output.file, append = TRUE)
  
  str2 = paste("JOBID=$(qsub s1-", dataset$Name, ".sh)", sep="")
  write(str2, file = output.file, append = TRUE)
  
  write("echo \" O job é $JOBID\" ", file = output.file, append = TRUE)
  
  str3 = paste("while [[ -n $(qstat -an1u elaine | grep \"$JOBID\") ]] ; do", sep="")
  write(str3, file = output.file, append = TRUE)
  write("sleep 5", file = output.file, append = TRUE)
  write("done", file = output.file, append = TRUE)
  
  str4 = paste("rclone -v copy /mnt/nfs/home/elaine/MultiLabelSimilaritiesMeasures/Reports/", dataset$Name ," nuvem:[2022]ResultadosExperimentos/Similaridades/", dataset$Name, sep="")
  write(str4, file = output.file, append = TRUE)
  
  str5 = paste("rm -r /lustre/elaine/s1-", dataset$Name, "/*", sep="")
  write(str5, file = output.file, append = TRUE)
  
  write("conda deactivate", file = output.file, append = TRUE)
  
  close(output.file)
  
  i = i + 1
  gc()
  
}
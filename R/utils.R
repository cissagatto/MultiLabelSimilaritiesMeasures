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

FolderRoot = "~/MultiLabelSimilaritiesMeasures"
FolderScripts = paste(FolderRoot, "/R/", sep="")

##################################################################################################

diretorios <- function(dataset_name, FolderResults){
  retorno = list()
  
  if(dir.exists(FolderResults)==FALSE){dir.create(FolderResults)}
  
  FolderScripts =  paste(FolderRoot, "/R", sep="")
  if(dir.exists(FolderScripts)==FALSE){dir.create(FolderScripts)}
  
  FolderRD =  paste(FolderResults, "/", dataset_name, sep="")
  if(dir.exists(FolderRD)==FALSE){dir.create(FolderRD)}
  
  FolderDatasets = paste(FolderResults, "/Datasets", sep="")
  if(dir.exists(FolderDatasets)==FALSE){dir.create(FolderDatasets)}
  
  FolderDS = paste(FolderDatasets, "/", dataset_name, sep="")
  if(dir.exists(FolderDS)==FALSE){dir.create(FolderDS)}
  
  FolderCV = paste(FolderDS, "/CrossValidation", sep="")
  if(dir.exists(FolderCV)==FALSE){dir.create(FolderCV)}
  
  FolderNamesLabels = paste(FolderDS, "/NamesLabels", sep="")
  if(dir.exists(FolderNamesLabels)==FALSE){dir.create(FolderNamesLabels)}
  
  FolderLabelSpace = paste(FolderDS, "/LabelSpace", sep="")
  if(dir.exists(FolderLabelSpace)==FALSE){dir.create(FolderLabelSpace)}
  
  FolderTR = paste(FolderCV, "/Tr", sep="")
  if(dir.exists(FolderTR)==FALSE){dir.create(FolderTR)}
  
  FolderTS = paste(FolderCV, "/Ts", sep="")
  if(dir.exists(FolderTS)==FALSE){dir.create(FolderTS)}
  
  FolderVL = paste(FolderCV, "/Vl", sep="")
  if(dir.exists(FolderVL)==FALSE){dir.create(FolderVL)}

  
  retorno$FolderDS = FolderDS
  retorno$FolderResults = FolderResults
  retorno$FolderScripts = FolderScripts
  retorno$FolderRD = FolderRD
  retorno$FolderDatasets = FolderDatasets
  retorno$FolderCV = FolderCV
  retorno$FolderNamesLabels = FolderNamesLabels
  retorno$FolderLabelSpace = FolderLabelSpace
  retorno$FolderTR = FolderTR
  retorno$FolderTS = FolderTS
  retorno$FolderVL = FolderVL
  
  return(retorno)
  gc()
  
}



##################################################################################################
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
FolderScripts = paste("~/MultiLabelSimilaritiesMeasures", "/R", sep="")

##################################################################################################
# Load sources
setwd(FolderScripts)
source("functions_contingency_table_multilabel.R")

setwd(FolderScripts)
source("functions_measures_binary_data.R")

setwd(FolderScripts)
source("functions_multilabel_binary_measures.R")

##################################################################################################
# load packages
library(progress)
library(dplyr)

##################################################################################################

# open dataset
setwd("/home/elaine/MultiLabelSimilaritiesMeasures/GpositiveGO/CrossValidation/Tr")
dados = foreign::read.arff("GpositiveGO-Split-Tr-1.arff")

# open information about dataset
setwd(FolderRoot)
datasets = read.csv("datasets.csv")
ds = datasets[29,]
start.label = ds$LabelStart
end.label = ds$LabelEnd
num.labels = ds$Labels
num.instances = nrow(dados)
num.attributes = ncol(dados)

# get labels
labels = data.frame(dados[,start.label:end.label])

# get instances
instances = data.frame(dados[, ds$AttStart:ds$AttEnd])

# get instances columns names 
names_instances = colnames(instances)

# get labels columns names
names_labels = colnames(labels)

# compute a, b, c and d
res1 = compute.cont.table(labels, num.labels)

# results from res1
a = res1$ma
b = res1$mb
c = res1$mc
d = res1$md

# compute ab, ac, ad, bc, bd and cd
res2 = compute.marg.probs(labels, num.labels, a, b, c, d)

# results from res2
ab = res2$mab
ac = res2$mac
ad = res2$mad
bc = res2$mbc
bd = res2$mbd
cd = res2$mcd
n = res2$mn

# covariance
res3 = compute.covar(labels, num.labels, ad, bc)
res3

# get the names of similarities functions
funs = getNamesListFunctions()
funs

res4 = compute.measure(labels, num.labels, a, b, c, d, n, ample.e)
res4








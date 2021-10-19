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
# Options Configuration                                                                          #
##################################################################################################
options(java.parameters = "-Xmx32g")
options(show.error.messages = TRUE)
options(scipen=10)

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

compute.measure(labels, num.labels, a, b, c, d, n, funs[1], ample.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[2], anderberg.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[3], baroni.urbani.buser.1.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[4], baroni.urbani.buser.2.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[5], braun.banquet.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[6], bray.curtis.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[7], canberra.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[8], chord.e )
compute.measure(labels, num.labels, a, b, c, d, n, funs[9], cityblock.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[10], cole.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[11], cosine.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[12], czekanowski.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[13], dennis.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[14], dice.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[15], disperson.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[16], driver.kroeber.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[17], euclidean.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[18], eyraud.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[19], fager.mcgowan.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[20], faith.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[21], forbes.2.e )
compute.measure(labels, num.labels, a, b, c, d, n, funs[22], forbesi.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[23], fossum.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[24], gilbert.well.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[25], goodman.kruskal.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[26], gower.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[27], gower.legendre.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[28], hamann.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[29], hamming.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[30], hellinger.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[31], inner.product.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[32], intersection.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[33], jaccard.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[34], johnson.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[35], kulczynski.1.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[36], kulczynski.2.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[37], lance.williams.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[38], manhattan.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[39], mcconnaughey.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[40], mean.manhattan.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[41], michael.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[42], minowski.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[43], mountford.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[44], nei.li.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[45], ochiai.2.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[46], otsuka.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[47], pattern.difference.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[48], pearson.1.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[49], pearson.2.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[50], pearson.3.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[51], pearson.heron.1.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[52], pearson.heron.2.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[53], peirce.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[54], roger.tanimoto.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[55], russel.rao.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[56], shape.differnece.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[57], simpson.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[58], size.difference.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[59], sokal.michener.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[60], sokal.sneath.1.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[61], sokal.sneath.2.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[62], sokal.sneath.3.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[63], sokal.sneath.4.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[64], sokal.sneath.5.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[65], sorgenfrei.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[66], square.euclidean.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[67], stiles.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[68], tanimoto.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[69], tarantula.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[70], tarwid.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[71], three.w.jaccard.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[72], vari.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[73], yule.w.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[74], yuleq.1.e)
compute.measure(labels, num.labels, a, b, c, d, n, funs[75], yuleq.2.e)




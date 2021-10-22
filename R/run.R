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


##################################################################################################
# Folders
##################################################################################################
FolderScripts =  paste(FolderRoot, "/R", sep="")
FolderDocs =  paste(FolderRoot, "/Docs", sep="")
FolderReports =  paste(FolderRoot, "/Reports", sep="")
FolderDatasets = paste(FolderRoot, "/Datasets", sep="")


##################################################################################################
# Options Configuration                                                                          #
##################################################################################################
options(java.parameters = "-Xmx32g")
options(show.error.messages = TRUE)
options(scipen=10)


###############################################################################
# Load sources
###############################################################################
setwd(FolderScripts)
source("functions_contingency_table_multilabel.R")

setwd(FolderScripts)
source("functions_measures_binary_data.R")

setwd(FolderScripts)
source("functions_multilabel_binary_measures.R")


###############################################################################
# load packages
###############################################################################
library(progress)
library(dplyr)


###############################################################################
# 
###############################################################################

# open dataset
Folder = paste(FolderDatasets, "/GpositiveGO/CrossValidation/Tr", sep="")
setwd(Folder)
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
head(labels)

# get instances
instances = data.frame(dados[, ds$AttStart:ds$AttEnd])

# get instances columns names 
names_instances = colnames(instances)

# get labels columns names
names_labels = colnames(labels)

# compute a, b, c and d
res1 = compute.cont.table(labels, num.labels)
res1

# results from res1
m.a = res1$ma
m.b = res1$mb
m.c = res1$mc
m.d = res1$md

# compute ab, ac, ad, bc, bd and cd
res2 = compute.marg.probs(labels, num.labels, m.a, m.b, m.c, m.d)

# results from res2
m.ab = res2$mab
m.ac = res2$mac
m.ad = res2$mad
m.bc = res2$mbc
m.bd = res2$mbd
m.cd = res2$mcd
m.n = res2$mn

# covariance
res3 = compute.covar(labels, num.labels, m.ad, m.bc)
res3

# get the names of similarities functions
funs = getNamesListFunctions()
funs

setwd(FolderReports)

################################################################################
# 1. AMPLE
# retornou números positivos e infinitos
result1 = compute.measure(Fun = ample.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
write.csv(result1, "ample-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[1], ample.e.2)
write.csv(result2, "ample-2.csv")
rm(result2)

################################################################################
# 2. ANDERBERG
# deu diferente os resultados das funções
# retornou números negativos e positivos
result1 = compute.measure(Fun = anderberg.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
write.csv(result1, "anderberg-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[2], anderberg.e.2)
write.csv(result2, "anderberg-2.csv")
rm(result2)

################################################################################
# 3. BARONI URBANI BUSER 1
# retornou valores positivos entre 0 e 1
result1 = compute.measure(Fun = baroni.urbani.buser.1.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
write.csv(result1, "baroni-urbani-user-1-1.csv")
rm(result1)

results2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[3], baroni.urbani.buser.1.e.2)
write.csv(result2, "baroni-urbani-user-1-2.csv")
rm(result2)

################################################################################
# 4. BARONI URBANI BUSER 2
# retornou valores positivos e negativos. 
# A faixa deve estar entre -1 e +1
result1 = compute.measure(Fun = baroni.urbani.buser.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
write.csv(result1, "baroni-urbani-user-2-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[4], baroni.urbani.buser.2.e.2)
write.csv(result2, "baroni-urbani-user-2-2.csv")
rm(result2)

################################################################################
# 5. BRAUN BANQUET
# deram resultados diferentes
# retornou valores positivos entre 0 e 1
result1 = compute.measure(Fun = braun.banquet.e.1, a = m.a, m.b = m.b, c = m.c) 
write.csv(result1, "baroni-branquet-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[5], braun.banquet.e.2)
write.csv(result2, "baroni-branquet-2.csv")
rm(result2)

################################################################################
# 6. BRAY CURTIS
# retornou valores positivos entre 0 e 1
result1 = compute.measure(Fun = bray.curtis.e.1, a = m.a, b = m.b, c = m.c)
write.csv(result1, "bray-curtis-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[6], bray.curtis.e.2)
write.csv(result2, "bray-curtis-2.csv")
rm(result2)

################################################################################
# 7. CANBERRA
# retornou valores positivos dentro dos valores da tabela de contingência
result1 = compute.measure(Fun = canberra.e.1, b = m.b, c = m.c)
write.csv(result1, "canberra-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[7], canberra.e.2)
write.csv(result2, "canberra-2.csv")
rm(result2)

################################################################################
# 8. CHORD
# retornou valores positivos 
result1 = compute.measure(Fun = chord.e.1, a = m.a, b = m.b, c = m.c)
write.csv(result1, "chord-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[8], chord.e.2)
write.csv(result2, "chord-2.csv")
rm(result2)

################################################################################
# 9. CITYBLOCK
# retornou valores positivos dentro dos valores da tabela de contingência
result1 = compute.measure(Fun = cityblock.e.1, b = m.b, c = m.c)
write.csv(result1, "cityblock-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[9], cityblock.e.2)
write.csv(result2, "cityblock-2.csv")
rm(result2)

################################################################################
# 10. COLE
# retornou apenas INFs and NANs
result1 = compute.measure(Fun = cole.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
write.csv(result1, "cole-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[10], cole.e.2)
write.csv(result2, "cole-2.csv")
rm(result2)

################################################################################
# 11. Cosine
# retornou valores positivos
result1 = compute.measure(Fun = cosine.e.1, a = m.a, b = m.b, c = m.c)
write.csv(result1, "cosine-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[11], cosine.e.2)
write.csv(result2, "cosine-2.csv")
rm(result2)

################################################################################
# CZEKANOWSKI
# retornou valores positivos entre 0 e 1
result1 = compute.measure(Fun = czekanowski.e.1, a = m.a, b = m.b, c = m.c)
write.csv(result1, "czekanowski-1.csv")
rm(result1)

compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[12], czekanowski.e.2)
write.csv(result2, "czekanowski-2.csv")
rm(result2)

################################################################################
# Dennis
# retornou valores positivos e negativos
result1 = compute.measure(Fun = dennis.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
write.csv(result1, "dennis-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[13], dennis.e.2)
write.csv(result2, "dennis-2.csv")
rm(result2)

################################################################################
# Dice
# retornou valores positivos entre 0 e 1
result1 = compute.measure(Fun = dice.e.1, a = m.a, b = m.b, c = m.c)
write.csv(result1, "dice-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[14], dice.e.2)
write.csv(result2, "dice-2.csv")
rm(result2)

################################################################################
# Dispersion
# retornou valores positivos e negativos
result1 = compute.measure(Fun = disperson.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
write.csv(result1, "disperson-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[15], disperson.e.2)
write.csv(result2, "disperson-2.csv")
rm(result2)

################################################################################
# Driver Kroeber
# retornou valores positivos entre 0 e 1
result1 = compute.measure(Fun =  driver.kroeber.e.1, a = m.a, b = m.b, c = m.c)
write.csv(result1, "driver-kroeber-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[16], driver.kroeber.e.2)
write.csv(result2, "driver-kroeber-2.csv")
rm(result2)

################################################################################
# Euclidean
# retornou valores positivos
result1 = compute.measure(Fun = euclidean.e.1, b = m.b, c = m.c)
write.csv(result1, "euclidean-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[17], euclidean.e.2)
write.csv(result2, "euclidean-2.csv")
rm(result2)

################################################################################
# Eyraud
# retornou valores positivos e negativos
result1 = compute.measure(Fun = eyraud.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
write.csv(result1, "eyraud-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[18], eyraud.e.2)
write.csv(result2, "eyraud-2.csv")
rm(result2)

################################################################################
# Fager Mcgowman
# deu diferente
# retornou valores positivos e negativos
result1 = compute.measure(Fun = fager.mcgowan.e.1, a = m.a, b = m.b, c = m.c)
write.csv(result1, "fager-mcgowan-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[19], fager.mcgowan.e.2)
write.csv(result2, "fager-mcgowan-2.csv")
rm(result2)

################################################################################
# Faith
# retornou valores positivos
result1 = compute.measure(Fun = faith.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
write.csv(result1, "faith-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[20], faith.e.2)
write.csv(result2, "faith-2.csv")
rm(result2)

################################################################################
# Forbes
# retornou valores positivos
result1 = compute.measure(Fun = forbes.2.e.1, a = m.a, b = m.b, c = m.c, n = m.n)
write.csv(result1, "forbes-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[21], forbes.2.e.2)
write.csv(result2, "forbes-2.csv")
rm(result2)

################################################################################
# Forbesi
# deu diferente
# retornou valores positivos
result1 = compute.measure(Fun = forbesi.e.1, a = m.a, b = m.b, c = m.c)
write.csv(result1, "forbesi-1.csv")
rm(result1)

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[22], forbesi.e.2)
write.csv(result2, "forbesi-2.csv")
rm(result2)

################################################################################
# deu difernete
# retornou valores positivos
result1 = compute.measure(Fun = fossum.e.1, a = m.a, b = m.b, c = m.c)
write.csv(result1, "fossum.e-1.csv")

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,
                            funs[23], fossum.e.2)
write.csv(result2, "fossum.e-2.csv")

################################################################################
# retornou valores positivos, negativos e INFs
result1 = compute.measure(Fun = gilbert.well.e.1, a = m.a, b = m.b, 
                          c = m.c, n = m.n)
write.csv(result1, "gilbert-well-1.csv")

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,
                            funs[24], gilbert.well.e.2)
write.csv(result2, "gilbert-well-2.csv")

################################################################################
# deu diferente
# retornou valores positivos e negativos
result1 = compute.measure(Fun = goodman.kruskal.e.1, a = m.a, b = m.b, 
                          c = m.c, d = m.d, n = m.n)
write.csv(result1, "goodman-kruskal-1.csv")

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,
                            funs[25], goodman.kruskal.e.2)
write.csv(result2, "goodman-kruskal-2.csv")

################################################################################
# retornou valores positivos

result1 = compute.measure(Fun = gower.e.2, a = m.a, b = m.b, 
                          c = m.c, d = m.d, n = m.n)
write.csv(result1, "goodman-kruskal-1.csv")

result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,
                            funs[26], gower.e.2)
write.csv(result2, "goodman-kruskal-2.csv")

compute.measure(labels, num.labels, a, b, c, d, n, )

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[27], gower.legendre.e)

################################################################################
# retornou valores positivos e negativos
compute.measure(labels, num.labels, a, b, c, d, n, funs[28], hamann.e)

################################################################################
# retornou valores positivos dentro dos valores da tabela de contingência
compute.measure(labels, num.labels, a, b, c, d, n, funs[29], hamming.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[30], hellinger.e)

################################################################################
# retornou valores positivos dentro dos valores da tabela de contingência
compute.measure(labels, num.labels, a, b, c, d, n, funs[31], inner.product.e)

################################################################################
# retornou valores positivos dentro dos valores da tabela de contingência
compute.measure(labels, num.labels, a, b, c, d, n, funs[32], intersection.e)

################################################################################
# retornou valores positivos entre o e 1
compute.measure(Fun = jaccard.e, a = m.a, b = m.b, c = m.c)
compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[33], jaccard.e.2)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[34], johnson.e)

################################################################################
# retornou valores positivos e INFs
compute.measure(labels, num.labels, a, b, c, d, n, funs[35], kulczynski.1.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[36], kulczynski.2.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[37], lance.williams.e)

################################################################################
# retornou valores positivos dentro dos valores da tabela de contingência
compute.measure(labels, num.labels, a, b, c, d, n, funs[38], manhattan.e)

################################################################################
# retornou valores positivos e negativos
compute.measure(labels, num.labels, a, b, c, d, n, funs[39], mcconnaughey.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[40], mean.manhattan.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[41], michael.e)

################################################################################
# retornou valores positivos dentro dos valores da tabela de contingência
compute.measure(labels, num.labels, a, b, c, d, n, funs[42], minowski.e)

################################################################################
# retornou valores positivos e INFs
compute.measure(labels, num.labels, a, b, c, d, n, funs[43], mountford.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[44], nei.li.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[45], ochiai.2.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[46], otsuka.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[47], pattern.difference.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[48], pearson.1.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[49], pearson.2.e)

################################################################################
# retornou valores positivos e NANs
compute.measure(labels, num.labels, a, b, c, d, n, funs[50], pearson.3.e)

################################################################################
# retornou valores positivos e negativos
compute.measure(labels, num.labels, a, b, c, d, n, funs[51], pearson.heron.1.e)

################################################################################
# retornou valores positivos, negativos e NANs
compute.measure(labels, num.labels, a, b, c, d, n, funs[52], pearson.heron.2.e)

################################################################################
# retornou valores positivos e NANs
compute.measure(labels, num.labels, a, b, c, d, n, funs[53], peirce.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[54], roger.tanimoto.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[55], russel.rao.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[56], shape.differnece.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[57], simpson.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[58], size.difference.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[59], sokal.michener.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[60], sokal.sneath.1.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[61], sokal.sneath.2.e)

################################################################################
# retornou valores positivos e INFs
compute.measure(labels, num.labels, a, b, c, d, n, funs[62], sokal.sneath.3.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[63], sokal.sneath.4.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[64], sokal.sneath.5.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[65], sorgenfrei.e)

################################################################################
# retornou valores positivos dentro dos valores da tabela de contingência
compute.measure(labels, num.labels, a, b, c, d, n, funs[66], square.euclidean.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[67], stiles.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[68], tanimoto.e)

################################################################################
# retornou valores positivos e INFs
compute.measure(labels, num.labels, a, b, c, d, n, funs[69], tarantula.e)

################################################################################
# retornou valores positivos, negativos e INFs
compute.measure(labels, num.labels, a, b, c, d, n, funs[70], tarwid.e)

################################################################################
# retornou valores positivos entre 0 e 1
compute.measure(labels, num.labels, a, b, c, d, n, funs[71], three.w.jaccard.e)

################################################################################
# retornou valores positivos
compute.measure(labels, num.labels, a, b, c, d, n, funs[72], vari.e)

################################################################################
# retornou valores positivos, negativos e INFs
compute.measure(labels, num.labels, a, b, c, d, n, funs[73], yule.w.e)

################################################################################
# retornou valores positivos e negativos entre -1 e +1
compute.measure(labels, num.labels, a, b, c, d, n, funs[74], yuleq.1.e)

################################################################################
# retornou valores positivos e INFs
compute.measure(labels, num.labels, a, b, c, d, n, funs[75], yuleq.2.e)




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

###############################################################################
# Load sources
###############################################################################
cat("\nCarregando scripts")

setwd(FolderScripts)
source("libraries.R")

setwd(FolderScripts)
source("utils.R")

setwd(FolderScripts)
source("functions_contingency_table_multilabel.R")

setwd(FolderScripts)
source("functions_measures_binary_data.R")

setwd(FolderScripts)
source("functions_multilabel_binary_measures.R")


##################################################################################################
# GET THE DIRECTORIES                                                                            #
##################################################################################################
cat("\nGet directories\n")
folder <- diretorios(dataset_name, FolderResults)


###############################################################################
# 
###############################################################################

executeMLSM_NCV <- function(number_dataset, FolderResults){
  
  # open dataset
  setwd(folder$FolderDS)
  dados = foreign::read.arff(paste(dataset_name, ".arff", sep=""))
  
  # open information about dataset
  setwd(FolderRoot)
  datasets = read.csv("datasets.csv")
  ds = datasets[number_dataset,]
  start.label = ds$LabelStart
  end.label = ds$LabelEnd
  start.att  = as.numeric(ds$AttStart)
  end.att  = as.numeric(ds$AttEnd)
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
  
  cat("\nCompute a, b, c and d\n\n")
  res1 = compute.cont.table(labels, num.labels)
  setwd(folder$FolderRD)
  write.csv(res1, paste(dataset_name, "-contingecy-table.csv"))
  
  # results from res1
  m.a = res1$ma
  m.b = res1$mb
  m.c = res1$mc
  m.d = res1$md
  
  cat("\nCompute compute ab, ac, ad, bc, bd and cd\n\n")
  res2 = compute.marg.probs(labels, num.labels, m.a, m.b, m.c, m.d)
  setwd(folder$FolderRD)
  write.csv(res2, paste(dataset_name, "-marginal-probabilities.csv"))
  
  # results from res2
  m.ab = res2$mab
  m.ac = res2$mac
  m.ad = res2$mad
  m.bc = res2$mbc
  m.bd = res2$mbd
  m.cd = res2$mcd
  m.n = res2$mn
  
  cat("\nCompute covariance\n\n")
  res3 = compute.covar(labels, num.labels, m.ad, m.bc)
  setwd(folder$FolderRD)
  write.csv(res3, paste(dataset_name, "-covariance.csv", sep=""))
  
  cat("\nCompute total classes and conditional probabilities\n\n")
  res4 <- computeInitialClassProbabilitiesTotals(labels, num.labels)
  setwd(folder$FolderRD)
  write.csv(res4$class.probabilities, paste(dataset_name, "-conditional-probabilities.csv", sep=""))
  write.csv(res4$class.totals, paste(dataset_name, "-total-per-labels.csv", sep=""))
  
  # get the names of similarities functions
  funs = getNamesListFunctions()
  
  setwd(folder$FolderRD)
  
  ################################################################################
  # 1. AMPLE
  # retornou números positivos e infinitos
  result1 = compute.measure(Fun = ample.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "ample-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[1], ample.e.2)
  write.csv(result2, "ample-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 2. ANDERBERG
  # ------------->>>>>>>>>>>> DEU DIFERENTE!!!
  # retornou números negativos e positivos
  result1 = compute.measure(Fun = anderberg.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
  write.csv(result1, "anderberg-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[2], anderberg.e.2)
  write.csv(result2, "anderberg-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 3. BARONI URBANI BUSER 1
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = baroni.urbani.buser.1.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "baroni-urbani-user-1-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[3], baroni.urbani.buser.1.e.2)
  write.csv(result2, "baroni-urbani-user-1-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 4. BARONI URBANI BUSER 2
  # retornou valores positivos e negativos. 
  # A faixa deve estar entre -1 e +1
  result1 = compute.measure(Fun = baroni.urbani.buser.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "baroni-urbani-user-2-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[4], baroni.urbani.buser.2.e.2)
  write.csv(result2, "baroni-urbani-user-2-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 5. BRAUN BANQUET
  # -------------------->>>>>>>>>>>>> deram resultados diferentes
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = braun.banquet.e.1, a = m.a, m.b = m.b, c = m.c) 
  write.csv(result1, "baroni-branquet-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[5], braun.banquet.e.2)
  write.csv(result2, "baroni-branquet-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 6. BRAY CURTIS
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = bray.curtis.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "bray-curtis-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[6], bray.curtis.e.2)
  write.csv(result2, "bray-curtis-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 7. CANBERRA
  # retornou valores positivos dentro dos valores da tabela de contingência
  result1 = compute.measure(Fun = canberra.e.1, b = m.b, c = m.c)
  write.csv(result1, "canberra-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[7], canberra.e.2)
  write.csv(result2, "canberra-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 8. CHORD
  # retornou valores positivos 
  result1 = compute.measure(Fun = chord.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "chord-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[8], chord.e.2)
  write.csv(result2, "chord-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 9. CITYBLOCK
  # retornou valores positivos dentro dos valores da tabela de contingência
  result1 = compute.measure(Fun = cityblock.e.1, b = m.b, c = m.c)
  write.csv(result1, "cityblock-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[9], cityblock.e.2)
  write.csv(result2, "cityblock-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 10. COLE
  # retornou apenas INFs and NANs
  result1 = compute.measure(Fun = cole.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "cole-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[10], cole.e.2)
  write.csv(result2, "cole-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 11. Cosine
  # retornou valores positivos
  result1 = compute.measure(Fun = cosine.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "cosine-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[11], cosine.e.2)
  write.csv(result2, "cosine-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 12. CZEKANOWSKI
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = czekanowski.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "czekanowski-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[12], czekanowski.e.2)
  write.csv(result2, "czekanowski-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 13. Dennis
  # retornou valores positivos e negativos
  result1 = compute.measure(Fun = dennis.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
  write.csv(result1, "dennis-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[13], dennis.e.2)
  write.csv(result2, "dennis-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 14. Dice
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = dice.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "dice-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[14], dice.e.2)
  write.csv(result2, "dice-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 15. Dispersion
  # retornou valores positivos e negativos
  result1 = compute.measure(Fun = disperson.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "disperson-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[15], disperson.e.2)
  write.csv(result2, "disperson-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 16. Driver Kroeber
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun =  driver.kroeber.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "driver-kroeber-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[16], driver.kroeber.e.2)
  write.csv(result2, "driver-kroeber-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 17. Euclidean
  # retornou valores positivos
  result1 = compute.measure(Fun = euclidean.e.1, b = m.b, c = m.c)
  write.csv(result1, "euclidean-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[17], euclidean.e.2)
  write.csv(result2, "euclidean-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 18. Eyraud
  # retornou valores positivos e negativos
  result1 = compute.measure(Fun = eyraud.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
  write.csv(result1, "eyraud-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[18], eyraud.e.2)
  write.csv(result2, "eyraud-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 19. Fager Mcgowman
  # --------------->>>>>>>>>>>>> deu diferente
  # retornou valores positivos e negativos
  result1 = compute.measure(Fun = fager.mcgowan.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "fager-mcgowan-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[19], fager.mcgowan.e.2)
  write.csv(result2, "fager-mcgowan-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 20. Faith
  # retornou valores positivos
  result1 = compute.measure(Fun = faith.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "faith-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[20], faith.e.2)
  write.csv(result2, "faith-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 21. Forbes
  # -------------------->>>>>>>>>>>>>>>>>> DEU DIFERENTE
  # retornou valores positivos
  result1 = compute.measure(Fun = forbes.2.e.1, a = m.a, b = m.b, c = m.c, n = m.n)
  write.csv(result1, "forbes-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[21], forbes.2.e.2)
  write.csv(result2, "forbes-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 22. Forbesi
  # --------------------->>>>>>>>>>>>>> deu diferente
  # retornou valores positivos
  result1 = compute.measure(Fun = forbesi.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "forbesi-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[22], forbesi.e.2)
  write.csv(result2, "forbesi-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 23. Fossum
  # retornou valores positivos
  result1 = compute.measure(Fun = fossum.e.1, a = m.a, b = m.b, c = m.c, n = m.n)
  write.csv(result1, "fossum.e-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[23], fossum.e.2)
  write.csv(result2, "fossum.e-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 24. Gilber Well
  # retornou valores positivos, negativos e INFs
  result1 = compute.measure(Fun = gilbert.well.e.1, a = m.a, b = m.b, c = m.c, n = m.n)
  write.csv(result1, "gilbert-well-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[24], gilbert.well.e.2)
  write.csv(result2, "gilbert-well-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 25. Goodman Kruskal 
  # ------->>>>>>>>>>>>>>>>>>>> deu diferente
  # retornou valores positivos e negativos
  result1 = compute.measure(Fun = goodman.kruskal.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
  write.csv(result1, "goodman-kruskal-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[25], goodman.kruskal.e.2)
  write.csv(result2, "goodman-kruskal-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 26. Gower 
  # retornou valores positivos
  result1 = compute.measure(Fun = gower.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "gower-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[26], gower.e.2)
  write.csv(result2, "gower-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 27. Gower Legender 
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = gower.legendre.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "gower-legendre-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[27], gower.legendre.e.2)
  write.csv(result2, "gower-legendre-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 28. Hamann
  # retornou valores positivos e negativos
  result1 = compute.measure(Fun = hamann.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "hamann-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[28], hamann.e.2)
  write.csv(result2, "hamann-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 29. Hamming 
  # retornou valores positivos dentro dos valores da tabela de contingência
  result1 = compute.measure(Fun = hamming.e.1, b = m.b, c = m.c)
  write.csv(result1, "hamming-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[29], hamming.e.2)
  write.csv(result2, "hamming-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 30. Hellinger
  # retornou valores positivos
  result1 = compute.measure(Fun = hellinger.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "hellinger-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[30], hellinger.e.2)
  write.csv(result2, "hellinger-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 31. Inner Product 
  # retornou valores positivos dentro dos valores da tabela de contingência
  result1 = compute.measure(Fun = inner.product.e.1, a = m.a, d = m.d)
  write.csv(result1, "inner-product-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[31], inner.product.e.2)
  write.csv(result2, "inner-product-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 32. Intersection 
  # retornou valores positivos dentro dos valores da tabela de contingência
  result1 = compute.measure(Fun = intersection.e.1, a = m.a)
  write.csv(result1, "intersection-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[32], intersection.e.2)
  write.csv(result2, "intersection-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 33. Jaccard
  # retornou valores positivos entre o e 1
  result1 = compute.measure(Fun = jaccard.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "jaccard-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[33], jaccard.e.2)
  write.csv(result2, "jaccard-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 34. Johnson
  # retornou valores positivos
  result1 = compute.measure(Fun = jonhson.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "jonhson-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[34], jonhson.e.2)
  write.csv(result2, "jonhson-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 35. kulczynski 1
  # retornou valores positivos e INFs
  result1 = compute.measure(Fun = kulczynski.1.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "kulczynski-1-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[35], kulczynski.1.e.2)
  write.csv(result2, "kulczynski-1-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 36. kulczynski 2
  # retornou valores positivos
  result1 = compute.measure(Fun = kulczynski.2.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "kulczynski-2-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[36], kulczynski.2.e.2)
  write.csv(result2, "kulczynski-2-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 37. Lance William
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = lance.williams.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "lance-williams-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[37], lance.williams.e.2)
  write.csv(result2, "lance-williams-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 38. Manhatan
  # retornou valores positivos dentro dos valores da tabela de contingência
  result1 = compute.measure(Fun = manhattan.e.1, b = m.b, c = m.c)
  write.csv(result1, "manhattan-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[38], manhattan.e.2)
  write.csv(result2, "manhattan-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 39. Mcconnaughey
  # retornou valores positivos e negativos
  result1 = compute.measure(Fun = mcconnaughey.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "mcconnaughey-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[39], mcconnaughey.e.2)
  write.csv(result2, "mcconnaughey-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 40. Mean Manhatan
  # retornou valores positivos
  result1 = compute.measure(Fun = mean.manhattan.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "mean-manhattan-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[40], mean.manhattan.e.2)
  write.csv(result2, "mean-manhattan-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 41. Michael
  # retornou valores positivos
  result1 = compute.measure(Fun = michael.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "michael-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[41], michael.e.2)
  write.csv(result2, "michael-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 42. Minowski
  # retornou valores positivos dentro dos valores da tabela de contingência
  result1 = compute.measure(Fun = minowski.e.1, b = m.b, c = m.c)
  write.csv(result1, "minowski-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[42], minowski.e.2)
  write.csv(result2, "minowski-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 43. Mountford
  # retornou valores positivos e INFs
  result1 = compute.measure(Fun = mountford.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "mountford-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[43], mountford.e.2)
  write.csv(result2, "mountford-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 44. Nei Li 
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = nei.li.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "nei-li-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[44], nei.li.e.2)
  write.csv(result2, "nei-li-1.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 45. Ochiai 2
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = ochiai.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "ochiai-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[45], ochiai.2.e.2)
  write.csv(result2, "ochiai-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 46. Otsuka
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = otsuka.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "otsuka-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[46], otsuka.e.2)
  write.csv(result2, "otsuka-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 47. Pattern Difference
  # retornou valores positivos
  result1 = compute.measure(Fun = pattern.difference.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "pattern-difference-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[47], pattern.difference.e.2)
  write.csv(result2, "pattern-difference-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 48. Pearson 1
  # retornou valores positivos
  result1 = compute.measure(Fun = pearson.1.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
  write.csv(result1, "pearson-1-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[48], pearson.1.e.2)
  write.csv(result2, "pearson-1-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 49. Pearson 2
  # retornou valores positivos
  result1 = compute.measure(Fun = pearson.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
  write.csv(result1, "pearson-2-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[49], pearson.2.e.2)
  write.csv(result2, "pearson-2-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 50. Pearson 3
  # retornou valores positivos e NANs
  result1 = compute.measure(Fun = pearson.3.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
  write.csv(result1, "pearson-3-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[50], pearson.3.e.2)
  write.csv(result2, "pearson-3-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 51. Pearson Heron 1
  # retornou valores positivos e negativos
  result1 = compute.measure(Fun = pearson.heron.1.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "pearson-heron-1-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[51], pearson.heron.1.e.2)
  write.csv(result2, "pearson-heron-1-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 52. Pearson Heron 2
  # retornou valores positivos, negativos e NANs
  result1 = compute.measure(Fun = pearson.heron.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "pearson-heron-2-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[52], pearson.heron.2.e.2)
  write.csv(result2, "pearson-heron-2-.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 53. Peirce
  # retornou valores positivos e NANs
  result1 = compute.measure(Fun = peirce.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "peirce-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[53], peirce.e.2)
  write.csv(result2, "peirce-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 54. Roger
  # retornou valores positivos
  result1 = compute.measure(Fun = roger.tanimoto.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "roger-tanimoto-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[54], roger.tanimoto.e.2)
  write.csv(result2, "roger-tanimoto-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 55. Russel Rao
  # retornou valores positivos
  result1 = compute.measure(Fun = russel.rao.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "russel-rao-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[55], russel.rao.e.2)
  write.csv(result2, "russel-rao-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 56. Shape Difference
  # retornou valores positivos
  result1 = compute.measure(Fun = shape.differnece.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
  write.csv(result1, "shape-difference-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[56], shape.differnece.e.2)
  write.csv(result2, "shape-difference-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 57. Simpson
  # -------------->>>>>>>>>>> deu diferente
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = simpson.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "simpson-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[57], simpson.e.2)
  write.csv(result2, "simpson-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 58. Size Difference
  # retornou valores positivos
  result1 = compute.measure(Fun = size.difference.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "size-difference-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[58], size.difference.e.2)
  write.csv(result2, "size-difference-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 59. Sokal Michener
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = sokal.michener.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "sokal-michener-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[59], sokal.michener.e.2)
  write.csv(result2, "sokal-michener-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 60. Sokal Sneath 1
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = sokal.sneath.1.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "sokal-sneath-1-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[60], sokal.sneath.1.e.2)
  write.csv(result2, "sokal-sneath-1-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 61. Sokal Sneath 2
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = sokal.sneath.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "sokal-sneath-2-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[61], sokal.sneath.2.e.2)
  write.csv(result2, "sokal-sneath-2-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 62. Sokal Sneath 3
  # retornou valores positivos e INFs
  result1 = compute.measure(Fun = sokal.sneath.3.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "sokal-sneath-3-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[62], sokal.sneath.3.e.2)
  write.csv(result2, "sokal-sneath-3-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 63. Sokal Sneath 4
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = sokal.sneath.4.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "sokal-sneath-4-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[63], sokal.sneath.4.e.2)
  write.csv(result2, "sokal-sneath-4-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 64. Sokal Sneath 5 
  # retornou valores positivos
  result1 = compute.measure(Fun = sokal.sneath.5.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "sokal-sneath-5-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[64], sokal.sneath.5.e.2)
  write.csv(result2, "sokal-sneath-5-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 65. Sorgenfrei
  # retornou valores positivos
  result1 = compute.measure(Fun = sorgenfrei.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "songerfrei-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[65], sorgenfrei.e.2)
  write.csv(result2, "songerfrei-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 66. Square Euclidean
  # retornou valores positivos dentro dos valores da tabela de contingência
  result1 = compute.measure(Fun = square.euclidean.e.1, b = m.b, c = m.c)
  write.csv(result1, "square-euclidean-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[66], square.euclidean.e.2)
  write.csv(result2, "square-euclidean-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 67. Stiles
  # retornou valores positivos
  result1 = compute.measure(Fun = stiles.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
  write.csv(result1, "stiles-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[67], stiles.e.2)
  write.csv(result2, "stiles-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 68. Tanimoto
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = tanimoto.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "tanimoto-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[68], tanimoto.e.2)
  write.csv(result2, "tanimoto-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 69. Tarantula
  # retornou valores positivos e INFs
  result1 = compute.measure(Fun = tarantula.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "tarantula-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[69], tarantula.e.2)
  write.csv(result2, "tarantula-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 70. Tarwid
  # retornou valores positivos, negativos e INFs
  result1 = compute.measure(Fun = tarwid.e.1, a = m.a, b = m.b, c = m.c, n = m.n)
  write.csv(result1, "tarwid-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[70], tarwid.e.2)
  write.csv(result2, "tarwid-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 71. 3W Jaccard
  # retornou valores positivos entre 0 e 1
  result1 = compute.measure(Fun = three.w.jaccard.e.1, a = m.a, b = m.b, c = m.c)
  write.csv(result1, "three-w-jaccard-2.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[71], three.w.jaccard.e.2)
  write.csv(result2, "three-w-jaccard-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 72. Vari
  # retornou valores positivos
  result1 = compute.measure(Fun = vari.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "vari-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[72], vari.e.2)
  write.csv(result2, "vari-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 73. Yule W
  # retornou valores positivos, negativos e INFs
  result1 = compute.measure(Fun = yule.w.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "yule-w-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[73], yule.w.e.2)
  write.csv(result2, "yule-w-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 74. Yuleq 1
  # retornou valores positivos e negativos entre -1 e +1
  result1 = compute.measure(Fun = yuleq.1.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "yuleq-1-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[74], yuleq.1.e.2)
  write.csv(result2, "yuleq-1-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  ################################################################################
  # 75. Yuleq 2
  # retornou valores positivos e INFs
  result1 = compute.measure(Fun = yuleq.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
  write.csv(result1, "yuleq-2-1.csv")
  
  result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[75], yuleq.2.e.2)
  write.csv(result2, "yuleq-2-2.csv")
  
  result1 == result2
  rm(result1, result2)
  
  
  ################################################################################
  # measure from cerri and mauri
  setwd(folder$FolderDS)
  dados.2 = read.csv(paste(dataset_name, ".csv", sep=""), stringsAsFactors = FALSE)
  result0 = compute.cerri.ferrandin(dados.2, labels, ds, k=6)
  setwd(folder$FolderRD)
  write.csv(result0, paste(dataset_name, "cerri-ferrandin.csv"), row.names = FALSE)  

  gc()  
  
}


################################################################################
# any errors, please, contact me: elainececiliagatto@gmail.com                 #
################################################################################
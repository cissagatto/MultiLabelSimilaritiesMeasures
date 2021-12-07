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
if (sistema[1] == "Linux") {
  FolderRoot = paste("/home/", sistema[7], "/MultiLabelSimilaritiesMeasures", sep = "")
  setwd(FolderRoot)
} else {
  FolderRoot = paste("C:/Users/", sistema[7], "/MultiLabelSimilaritiesMeasures", sep = "")
  setwd(FolderRoot)
}
setwd(FolderRoot)
FolderScripts =  paste(FolderRoot, "/R", sep = "")

##################

executeMLSM_CV <- function(ds, number_dataset, number_cores, number_folds, FolderResults){
  
  dataset_name = ds$Name
  
  ##################################################################################################
  if(number_cores == 0){
    cat("\nZero is a disallowed value for number_cores. Please choose a value greater than or equal to 1.")
  } else {
    cl <- parallel::makeCluster(number_cores)
    doParallel::registerDoParallel(cl)
    print(cl)
    
    if(number_cores==1){
      cat("\n\n################################################################################################")
      cat("\n# Running Sequentially!                                                                          #")
      cat("\n##################################################################################################\n\n") 
    } else {
      cat("\n\n################################################################################################")
      cat("\n# Running in parallel with ", number_cores, " cores!                                             #")
      cat("\n##################################################################################################\n\n") 
    }
  }
  cl = cl
  
  f = 1
  mlsmParalel <- foreach(f = 1:number_folds) %dopar% {
    
    cat("\nFold ", f)
    
    ##################################################################################################
    # Configures the workspace according to the operating system                                     #
    ##################################################################################################
    sistema = c(Sys.info())
    FolderRoot = ""
    if (sistema[1] == "Linux") {
      FolderRoot = paste("/home/", sistema[7], "/MultiLabelSimilaritiesMeasures", sep = "")
      setwd(FolderRoot)
      print(FolderRoot)
    } else {
      FolderRoot = paste("C:/Users/", sistema[7], "/MultiLabelSimilaritiesMeasures", sep = "")
      setwd(FolderRoot)
      print(FolderRoot)
    }
    setwd(FolderRoot)
    FolderScripts =  paste(FolderRoot, "/R", sep = "")
    
    
    ###############################################################################
    # Load sources
    ###############################################################################
    cat("\nCarregando scripts")
    
    cat("\nLibraries")
    setwd(FolderScripts)
    source("libraries.R")
    
    cat("\nUtils")
    setwd(FolderScripts)
    source("utils.R")
    
    cat("\nContingency Table")
    setwd(FolderScripts)
    source("functions_contingency_table_multilabel.R")
    
    cat("\nBinary Data")
    setwd(FolderScripts)
    source("functions_measures_binary_data.R")
    
    cat("\nBinary Measures")
    setwd(FolderScripts)
    source("functions_multilabel_binary_measures.R")
    
    ###############################################################################
    cat("\nGet folders")
    folder <- diretorios(ds$Name, FolderResults)
    
    ##############################################################
    cat("\nSplit Folder")
    FolderSplit = paste(folder$FolderResults, "/Split-", f, sep="")
    if(dir.exists(FolderSplit)==FALSE){dir.create(FolderSplit)}
    setwd(FolderSplit)
    
    ###############################################################################
    cat("\nOpen information about dataset")
    setwd(FolderRoot)
    datasets = read.csv("datasets.csv")
    start.label = ds$LabelStart
    end.label = ds$LabelEnd
    start.att  = as.numeric(ds$AttStart)
    end.att  = as.numeric(ds$AttEnd)
    num.labels = ds$Labels
    
    ###############################################################################
    cat("\nOpen dataset")
    setwd(folder$FolderTR)
    dados = foreign::read.arff(paste(dataset_name, "-Split-Tr-", f, ".arff", sep=""))
    
    # instances number
    num.instances = nrow(dados)
    
    # attributes number
    num.attributes = ncol(dados)
    
    # get labels
    labels = data.frame(dados[,start.label:end.label])
    
    # get instances
    instances = data.frame(dados[,ds$AttStart:ds$AttEnd])
    
    # get instances columns names 
    names_instances = colnames(instances)
    
    # get labels columns names
    names_labels = colnames(labels)
    
    cat("\ncompute a, b, c and d")
    res1 = compute.cont.table(labels, num.labels)
    setwd(FolderSplit)
    write.csv(res1, paste(dataset_name, "-contingecy-table.csv", sep=""))
    
    # results from res1
    m.a = res1$ma
    m.b = res1$mb
    m.c = res1$mc
    m.d = res1$md
    
    cat("\ncompute ab, ac, ad, bc, bd and cd")
    res2 = compute.marg.probs(labels, num.labels, m.a, m.b, m.c, m.d)
    setwd(FolderSplit)
    write.csv(res2, paste(dataset_name, "-marginal-probabilities.csv", sep=""))
    
    # results from res2
    m.ab = res2$mab
    m.ac = res2$mac
    m.ad = res2$mad
    m.bc = res2$mbc
    m.bd = res2$mbd
    m.cd = res2$mcd
    m.n = res2$mn
    
    cat("\ncompute covariance")
    res3 = compute.covar(labels, num.labels, m.ad, m.bc)
    setwd(FolderSplit)
    write.csv(res3, paste(dataset_name, "-covariance.csv", sep=""))
    
    cat("\ncompute total classes and conditional probabilities")
    res4 <- computeInitialClassProbabilitiesTotals(labels, num.labels, names_labels)
    setwd(FolderSplit)
    write.csv(res4$class.probabilities, paste(dataset_name, "-conditional-probabilities.csv", sep=""))
    write.csv(res4$class.totals, paste(dataset_name, "-total-per-labels.csv", sep=""))
    
    # get the names of similarities functions
    funs = getNamesListFunctions()
    
    setwd(FolderSplit)
    
    
    ################################################################################
    # 1. AMPLE
    # retornou números positivos e infinitos
    #result1 = compute.measure(Fun = ample.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    #result1[which(!is.finite(result1))] <- 0
    #write.csv(result1, "ample-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[1], ample.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "ample.csv")
    rm(result2)
    
    ################################################################################
    # 2. ANDERBERG
    # ------------->>>>>>>>>>>> DEU DIFERENTE!!!
    # retornou números negativos e positivos
    #result1 = compute.measure(Fun = anderberg.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
    #write.csv(result1, "anderberg-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[2], anderberg.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "anderberg.csv")
    rm(result2)
    
    ################################################################################
    # 3. BARONI URBANI BUSER 1
    # retornou valores positivos entre 0 e 1
    #result1 = compute.measure(Fun = baroni.urbani.buser.1.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    #write.csv(result1, "baroni-urbani-user-1-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[3], baroni.urbani.buser.1.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "baroni-urbani-user-1.csv")
    rm(result2)
    
    ################################################################################
    # 4. BARONI URBANI BUSER 2
    # retornou valores positivos e negativos. 
    # A faixa deve estar entre -1 e +1
    # result1 = compute.measure(Fun = baroni.urbani.buser.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "baroni-urbani-user-2-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[4], baroni.urbani.buser.2.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "baroni-urbani-user-2.csv")
    rm(result2)
    
    ################################################################################
    # 5. BRAUN BANQUET
    # -------------------->>>>>>>>>>>>> deram resultados diferentes
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = braun.banquet.e.1, a = m.a, m.b = m.b, c = m.c) 
    # write.csv(result1, "baroni-branquet-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[5], braun.banquet.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "baroni-branquet.csv")
    rm(result2)
    
    ################################################################################
    # 6. BRAY CURTIS
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = bray.curtis.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "bray-curtis-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[6], bray.curtis.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "bray-curtis.csv")
    rm(result2)
    
    ################################################################################
    # 7. CANBERRA
    # retornou valores positivos dentro dos valores da tabela de contingência
    # result1 = compute.measure(Fun = canberra.e.1, b = m.b, c = m.c)
    # write.csv(result1, "canberra-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[7], canberra.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "canberra.csv")
    rm(result2)
    
    ################################################################################
    # 8. CHORD
    # retornou valores positivos 
    # result1 = compute.measure(Fun = chord.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "chord-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[8], chord.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "chord.csv")
    rm(result2)
    
    ################################################################################
    # 9. CITYBLOCK
    # retornou valores positivos dentro dos valores da tabela de contingência
    # result1 = compute.measure(Fun = cityblock.e.1, b = m.b, c = m.c)
    # write.csv(result1, "cityblock-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[9], cityblock.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "cityblock.csv")
    rm(result2)
    
    ################################################################################
    # 10. COLE
    # retornou apenas INFs and NANs
    # result1 = compute.measure(Fun = cole.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "cole-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[10], cole.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "cole.csv")
    rm(result2)
    
    ################################################################################
    # 11. Cosine
    # retornou valores positivos
    # result1 = compute.measure(Fun = cosine.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "cosine-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[11], cosine.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "cosine.csv")
    rm(result2)
    
    ################################################################################
    # 12. CZEKANOWSKI
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = czekanowski.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "czekanowski-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[12], czekanowski.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "czekanowski.csv")
    rm(result2)
    
    ################################################################################
    # 13. Dennis
    # retornou valores positivos e negativos
    # result1 = compute.measure(Fun = dennis.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
    # write.csv(result1, "dennis-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[13], dennis.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "dennis.csv")
    rm(result2)
    
    ################################################################################
    # 14. Dice
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = dice.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "dice-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[14], dice.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "dice.csv")
    rm(result2)
    
    ################################################################################
    # 15. Dispersion
    # retornou valores positivos e negativos
    # result1 = compute.measure(Fun = disperson.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "disperson-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[15], disperson.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "disperson.csv")
    rm(result2)
    
    ################################################################################
    # 16. Driver Kroeber
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun =  driver.kroeber.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "driver-kroeber-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[16], driver.kroeber.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "driver-kroeber.csv")
    rm(result2)
    
    ################################################################################
    # 17. Euclidean
    # retornou valores positivos
    # result1 = compute.measure(Fun = euclidean.e.1, b = m.b, c = m.c)
    # write.csv(result1, "euclidean-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[17], euclidean.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "euclidean.csv")
    rm(result2)
    
    ################################################################################
    # 18. Eyraud
    # retornou valores positivos e negativos
    # result1 = compute.measure(Fun = eyraud.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
    # write.csv(result1, "eyraud-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[18], eyraud.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "eyraud.csv")
    rm(result2)
    
    ################################################################################
    # 19. Fager Mcgowman
    # --------------->>>>>>>>>>>>> deu diferente
    # retornou valores positivos e negativos
    # result1 = compute.measure(Fun = fager.mcgowan.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "fager-mcgowan-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[19], fager.mcgowan.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "fager-mcgowan.csv")
    rm(result2)
    
    ################################################################################
    # 20. Faith
    # retornou valores positivos
    # result1 = compute.measure(Fun = faith.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "faith-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[20], faith.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "faith.csv")
    rm(result2)
    
    ################################################################################
    # 21. Forbes
    # -------------------->>>>>>>>>>>>>>>>>> DEU DIFERENTE
    # retornou valores positivos
    # result1 = compute.measure(Fun = forbes.2.e.1, a = m.a, b = m.b, c = m.c, n = m.n)
    # write.csv(result1, "forbes-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n,funs[21], forbes.2.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "forbes.csv")
    rm(result2)
    
    ################################################################################
    # 22. Forbesi
    # --------------------->>>>>>>>>>>>>> deu diferente
    # retornou valores positivos
    # result1 = compute.measure(Fun = forbesi.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "forbesi-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[22], forbesi.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "forbesi.csv")
    rm(result2)
    
    ################################################################################
    # 23. Fossum
    # retornou valores positivos
    # result1 = compute.measure(Fun = fossum.e.1, a = m.a, b = m.b, c = m.c, n = m.n)
    # write.csv(result1, "fossum-e-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[23], fossum.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "fossum-e.csv")
    rm(result2)
    
    ################################################################################
    # 24. Gilber Well
    # retornou valores positivos, negativos e INFs
    # result1 = compute.measure(Fun = gilbert.well.e.1, a = m.a, b = m.b, c = m.c, n = m.n)
    # write.csv(result1, "gilbert-well-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[24], gilbert.well.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "gilbert-well.csv")
    rm(result2)
    
    ################################################################################
    # 25. Goodman Kruskal 
    # ------->>>>>>>>>>>>>>>>>>>> deu diferente
    # retornou valores positivos e negativos
    # result1 = compute.measure(Fun = goodman.kruskal.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
    # write.csv(result1, "goodman-kruskal-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[25], goodman.kruskal.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "goodman-kruskal.csv")
    rm(result2)
    
    ################################################################################
    # 26. Gower 
    # retornou valores positivos
    # result1 = compute.measure(Fun = gower.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "gower-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[26], gower.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "gower.csv")
    rm(result2)
    
    ################################################################################
    # 27. Gower Legender 
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = gower.legendre.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "gower-legendre-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[27], gower.legendre.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "gower-legendre.csv")
    rm(result2)
    
    ################################################################################
    # 28. Hamann
    # retornou valores positivos e negativos
    # result1 = compute.measure(Fun = hamann.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "hamann-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[28], hamann.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "hamann.csv")
    rm(result2)
    
    ################################################################################
    # 29. Hamming 
    # retornou valores positivos dentro dos valores da tabela de contingência
    # result1 = compute.measure(Fun = hamming.e.1, b = m.b, c = m.c)
    # write.csv(result1, "hamming-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[29], hamming.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "hamming.csv")
    rm(result2)
    
    ################################################################################
    # 30. Hellinger
    # retornou valores positivos
    # result1 = compute.measure(Fun = hellinger.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "hellinger-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[30], hellinger.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "hellinger.csv")
    rm(result2)
    
    ################################################################################
    # 31. Inner Product 
    # retornou valores positivos dentro dos valores da tabela de contingência
    # result1 = compute.measure(Fun = inner.product.e.1, a = m.a, d = m.d)
    # write.csv(result1, "inner-product-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[31], inner.product.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "inner-product.csv")
    rm(result2)
    
    ################################################################################
    # 32. Intersection 
    # retornou valores positivos dentro dos valores da tabela de contingência
    # result1 = compute.measure(Fun = intersection.e.1, a = m.a)
    # write.csv(result1, "intersection-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[32], intersection.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "intersection.csv")
    rm(result2)
    
    ################################################################################
    # 33. Jaccard
    # retornou valores positivos entre o e 1
    # result1 = compute.measure(Fun = jaccard.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "jaccard-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[33], jaccard.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "jaccard.csv")
    rm(result2)
    
    ################################################################################
    # 34. Johnson
    # retornou valores positivos
    # result1 = compute.measure(Fun = jonhson.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "jonhson-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[34], jonhson.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "jonhson.csv")
    rm(result2)
    
    ################################################################################
    # 35. kulczynski 1
    # retornou valores positivos e INFs
    # result1 = compute.measure(Fun = kulczynski.1.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "kulczynski-1-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[35], kulczynski.1.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "kulczynski-1.csv")
    rm(result2)
    
    ################################################################################
    # 36. kulczynski 2
    # retornou valores positivos
    # result1 = compute.measure(Fun = kulczynski.2.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "kulczynski-2-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[36], kulczynski.2.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "kulczynski-2.csv")
    rm(result2)
    
    ################################################################################
    # 37. Lance William
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = lance.williams.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "lance-williams-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[37], lance.williams.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "lance-williams.csv")
    rm(result2)
    
    ################################################################################
    # 38. Manhatan
    # retornou valores positivos dentro dos valores da tabela de contingência
    # result1 = compute.measure(Fun = manhattan.e.1, b = m.b, c = m.c)
    # write.csv(result1, "manhattan-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[38], manhattan.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "manhattan.csv")
    rm(result2)
    
    ################################################################################
    # 39. Mcconnaughey
    # retornou valores positivos e negativos
    # result1 = compute.measure(Fun = mcconnaughey.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "mcconnaughey-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[39], mcconnaughey.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "mcconnaughey.csv")
    rm(result2)
    
    ################################################################################
    # 40. Mean Manhatan
    # retornou valores positivos
    # result1 = compute.measure(Fun = mean.manhattan.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "mean-manhattan-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[40], mean.manhattan.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "mean-manhattan.csv")
    rm(result2)
    
    ################################################################################
    # 41. Michael
    # retornou valores positivos
    # result1 = compute.measure(Fun = michael.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "michael-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[41], michael.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "michael.csv")
    rm(result2)
    
    ################################################################################
    # 42. Minowski
    # retornou valores positivos dentro dos valores da tabela de contingência
    # result1 = compute.measure(Fun = minowski.e.1, b = m.b, c = m.c)
    #write.csv(result1, "minowski-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[42], minowski.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "minowski.csv")
    rm(result2)
    
    ################################################################################
    # 43. Mountford
    # retornou valores positivos e INFs
    # result1 = compute.measure(Fun = mountford.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "mountford-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[43], mountford.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "mountford.csv")
    rm(result2)
    
    ################################################################################
    # 44. Nei Li 
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = nei.li.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "nei-li-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[44], nei.li.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "nei-li.csv")
    rm(result2)
    
    ################################################################################
    # 45. Ochiai 1
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = ochiai.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "ochiai-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[45], ochiai.2.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "ochiai-1.csv")
    rm(result2)
    
    ################################################################################
    # 45. Ochiai 2
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = ochiai.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "ochiai-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[46], ochiai.2.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "ochiai-2.csv")
    rm(result2)
    
    ################################################################################
    # 46. Otsuka
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = otsuka.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "otsuka-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[47], otsuka.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "otsuka.csv")
    rm(result2)
    
    ################################################################################
    # 47. Pattern Difference
    # retornou valores positivos
    # result1 = compute.measure(Fun = pattern.difference.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "pattern-difference-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[48], pattern.difference.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pattern-difference.csv")
    rm(result2)
    
    ################################################################################
    # 48. Pearson 1
    # retornou valores positivos
    # result1 = compute.measure(Fun = pearson.1.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
    # write.csv(result1, "pearson-1-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[49], pearson.1.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pearson-1.csv")
    rm(result2)
    
    ################################################################################
    # 49. Pearson 2
    # retornou valores positivos
    # result1 = compute.measure(Fun = pearson.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
    # write.csv(result1, "pearson-2-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[50], pearson.2.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pearson-2.csv")
    rm(result2)
    
    ################################################################################
    # 50. Pearson 3
    # retornou valores positivos e NANs
    # result1 = compute.measure(Fun = pearson.3.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
    # write.csv(result1, "pearson-3-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[51], pearson.3.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pearson-3.csv")
    rm(result2)
    
    ################################################################################
    # 51. Pearson Heron 1
    # retornou valores positivos e negativos
    # result1 = compute.measure(Fun = pearson.heron.1.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "pearson-heron-1-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[52], pearson.heron.1.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pearson-heron-1.csv")
    rm(result2)
    
    ################################################################################
    # 52. Pearson Heron 2
    # retornou valores positivos, negativos e NANs
    # result1 = compute.measure(Fun = pearson.heron.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "pearson-heron-2-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[53], pearson.heron.2.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pearson-heron-2.csv")
    rm(result2)
    
    ################################################################################
    # 53. Peirce
    # retornou valores positivos e NANs
    # result1 = compute.measure(Fun = peirce.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "peirce-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[54], peirce.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "peirce.csv")
    rm(result2)
    
    ################################################################################
    # 54. Roger
    # retornou valores positivos
    # result1 = compute.measure(Fun = roger.tanimoto.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "roger-tanimoto-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[55], roger.tanimoto.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "roger-tanimoto.csv")
    rm(result2)
    
    ################################################################################
    # 55. Russel Rao
    # retornou valores positivos
    # result1 = compute.measure(Fun = russel.rao.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "russel-rao.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[56], russel.rao.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "russel-rao-2.csv")
    rm(result2)
    
    ################################################################################
    # 56. Shape Difference
    # retornou valores positivos
    # result1 = compute.measure(Fun = shape.differnece.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
    # write.csv(result1, "shape-difference-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[57], shape.differnece.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "shape-difference.csv")
    rm(result2)
    
    ################################################################################
    # 57. Simpson
    # -------------->>>>>>>>>>> deu diferente
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = simpson.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "simpson-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[58], simpson.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "simpson.csv")
    rm(result2)
    
    ################################################################################
    # 58. Size Difference
    # retornou valores positivos
    # result1 = compute.measure(Fun = size.difference.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "size-difference-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[59], size.difference.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "size-difference.csv")
    rm(result2)
    
    ################################################################################
    # 59. Sokal Michener
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = sokal.michener.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "sokal-michener-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[60], sokal.michener.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-michener.csv")
    rm(result2)
    
    ################################################################################
    # 60. Sokal Sneath 1
    # retornou valores positivos entre 0 e 1
    result1 = compute.measure(Fun = sokal.sneath.1.e.1, a = m.a, b = m.b, c = m.c)
    write.csv(result1, "sokal-sneath-1-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[61], sokal.sneath.1.e.2)
    write.csv(result2, "sokal-sneath-1.csv")
    
    rm(result2)
    
    ################################################################################
    # 61. Sokal Sneath 2
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = sokal.sneath.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "sokal-sneath-2-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[62], sokal.sneath.2.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-2.csv")
    rm(result2)
    
    ################################################################################
    # 62. Sokal Sneath 3
    # retornou valores positivos e INFs
    # result1 = compute.measure(Fun = sokal.sneath.3.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "sokal-sneath-3-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[63], sokal.sneath.3.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-3.csv")
    rm(result2)
    
    ################################################################################
    # 63. Sokal Sneath 4
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = sokal.sneath.4.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "sokal-sneath-4-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[64], sokal.sneath.4.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-4.csv")
    rm(result2)
    
    ################################################################################
    # 64. Sokal Sneath 5 
    # retornou valores positivos
    # result1 = compute.measure(Fun = sokal.sneath.5.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "sokal-sneath-5-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[65], sokal.sneath.5.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-5.csv")
    rm(result2)
    
    ################################################################################
    # 65. Sorgenfrei
    # retornou valores positivos
    # result1 = compute.measure(Fun = sorgenfrei.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "songerfrei-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[66], sorgenfrei.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "songerfrei.csv")
    rm(result2)
    
    ################################################################################
    # 66. Square Euclidean
    # retornou valores positivos dentro dos valores da tabela de contingência
    # result1 = compute.measure(Fun = square.euclidean.e.1, b = m.b, c = m.c)
    # write.csv(result1, "square-euclidean-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[67], square.euclidean.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "square-euclidean.csv")
    rm(result2)
    
    ################################################################################
    # 67. Stiles
    # retornou valores positivos
    # result1 = compute.measure(Fun = stiles.e.1, a = m.a, b = m.b, c = m.c, d = m.d, n = m.n)
    # write.csv(result1, "stiles-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[68], stiles.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "stiles.csv")
    rm(result2)
    
    ################################################################################
    # 68. Tanimoto
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = tanimoto.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "tanimoto-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[69], tanimoto.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "tanimoto.csv")
    rm(result2)
    
    ################################################################################
    # 69. Tarantula
    # retornou valores positivos e INFs entre 0 e 1
    # result1 = compute.measure(Fun = tarantula.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "tarantula-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[70], tarantula.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "tarantula.csv")
    rm(result2)
    
    ################################################################################
    # 70. Tarwid
    # retornou valores positivos, negativos e INFs
    # entre 0 e  .....
    # result1 = compute.measure(Fun = tarwid.e.1, a = m.a, b = m.b, c = m.c, n = m.n)
    # write.csv(result1, "tarwid-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[71], tarwid.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "tarwid.csv")
    rm(result2)
    
    ################################################################################
    # 71. 3W Jaccard
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = three.w.jaccard.e.1, a = m.a, b = m.b, c = m.c)
    # write.csv(result1, "three-w-jaccard-2.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[72], three.w.jaccard.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "three-w-jaccard.csv")
    rm(result2)
    
    ################################################################################
    # 72. Vari
    # retornou valores positivos entre 0 e 1
    # result1 = compute.measure(Fun = vari.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "vari-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[73], vari.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "vari.csv")
    rm(result2)
    
    ################################################################################
    # 73. Yule W
    # retornou valores positivos, negativos e INFs
    # entre 0 e ......
    # result1 = compute.measure(Fun = yule.w.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "yule-w-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[74], yule.w.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "yule-w.csv")
    rm(result2)
    
    ################################################################################
    # 74. Yuleq 1
    # retornou valores positivos e negativos entre -1 e +1
    # result1 = compute.measure(Fun = yuleq.1.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "yuleq-1-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[75], yuleq.1.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "yuleq-1.csv")
    rm(result2)
    
    ################################################################################
    # 75. Yuleq 2
    # retornou valores positivos e INFs entre 0 e .....
    # result1 = compute.measure(Fun = yuleq.2.e.1, a = m.a, b = m.b, c = m.c, d = m.d)
    # write.csv(result1, "yuleq-2-1.csv")
    
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, funs[76], yuleq.2.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "yuleq-2.csv")
    rm(result2)
    
    ################################################################################
    # measure from cerri and mauri
    Folder2 = paste(folder$FolderDatasets, "/", dataset_name ,"/CrossValidation/Tr", sep="")
    setwd(Folder2)
    dados.2 = read.csv(paste(dataset_name, "-Split-Tr-", f, ".csv", sep=""), stringsAsFactors = FALSE)
    result0 = compute.cerri.ferrandin(dados.2, labels, ds, k=6)
    setwd(FolderSplit)
    write.csv(result0, paste(dataset_name, "split-", f ,"-cerri-ferrandin.csv"))  
    
    ##################################################################################################
    # copy file                                                                                      #
    ##################################################################################################
    cat("\nCopy files \n")
    str4 = paste("cp -r ", FolderSplit, "/ ",  folder$FolderRD, sep="")
    print(system(str4))
    
    gc()
    
  } # fim do for each
  
  gc()
  
} # fim da função



################################################################################
# any errors, please, contact me: elainececiliagatto@gmail.com                 #
################################################################################
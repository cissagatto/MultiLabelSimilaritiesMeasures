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
    datasets = read.csv("datasets-2022.csv")
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
    
    cat("\ncompute independence")
    res6 = compute.indep(labels, num.labels, m.ad, m.bc)
    setwd(FolderSplit)
    write.csv(res6, paste(dataset_name, "-independence.csv", sep=""))
    
    cat("\ncompute total classes and conditional probabilities")
    res4 <- computeInitialClassProbabilitiesTotals(labels, num.labels, names_labels)
    setwd(FolderSplit)
    write.csv(res4$class.probabilities, paste(dataset_name, "-conditional-probabilities.csv", sep=""))
    write.csv(res4$class.totals, paste(dataset_name, "-total-per-labels.csv", sep=""))
    
    # get the names of similarities functions
    funs = getNamesListFunctions()
    
    setwd(FolderSplit)
    
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs$nomes_funcoes[1], ample.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "ample.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[2], anderberg.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "anderberg.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[3], baroni.urbani.buser.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "baroni-urbani-user-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[4], baroni.urbani.buser.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "baroni-urbani-user-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[5], braun.blanquet.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "braun-blanquet.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[6], bray.curtis.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "bray-curtis.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[7], canberra.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "canberra.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[8], chord.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "chord.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[9], cityblock.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "cityblock.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[10], clement.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "clement.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n,funs[11], cohen.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "cohen.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n,funs[12], cole.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "cole-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n,funs[13], cole.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "cole-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n,funs[14], cole.e.3)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "cole-3.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n,funs[15], cosine.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "cosine.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n,funs[16], czekanowski.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "czekanowski.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n,funs[17], dennis.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "dennis.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n,funs[18], dibby.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "dibby.csv")
    rm(result2)
    
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[19], dice.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "dice-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[20], dice.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "dice-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[21], dice.e.3)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "dice-3.csv")
    rm(result2)
    
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[22], disperson.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "disperson.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[23], doolittle.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "doolittle.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[24], driver.kroeber.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "driver-kroeber-1.csv")
    rm(result2)
    
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[25], driver.kroeber.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "driver-kroeber-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[26], euclidean.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "euclidean.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[27], eyraud.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "eyraud.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n,funs[28], fager.mcgowan.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "fager-mcgowan-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n,funs[28], fager.mcgowan.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "fager-mcgowan-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[30], faith.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "faith.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[31], fleiss.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "fleiss.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[32], forbes.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "forbes-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[33], forbes.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "forbes-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[34], forbesi.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "forbesi.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[35], fossum.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "fossum.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[36], gilbert.well.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "gilbert-well.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[37], goodman.kruskal.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "goodman-kruskal-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[38], goodman.kruskal.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "goodman-kruskal-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[39], gower.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "gower.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[40], gower.legendre.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "gower-legendre.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[41], hamann.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "hamann-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[42], hamann.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "hamann-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[43], hamming.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "hamming.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[44], harris.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "harris.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[45], hawkins.dotson.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "hawkins.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[46], hellinger.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "hellinger.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[47], inner.product.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "inner-product.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[48], intersection.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "intersection.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[49], jaccard.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "jaccard.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[50], jonhson.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "jonhson.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[51], kent.foster.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "kent-foster-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[52], kent.foster.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "kent-foster-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[53], kuder.richardson.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "kuder.richardson..csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[54], kulczynski.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "kulczynski-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[55], kulczynski.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "kulczynski-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[56], kulczynski.e.3)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "kulczynski-3.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[57], kulczynski.e.4)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "kulczynski-4.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[58], lance.williams.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "lance-williams.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[59], loevinger.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "loevinger.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, 
                                m.d, m.n, funs[60], manhattan.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "manhattan.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, 
                                m.d, m.n, funs[61], maxwell.pilliner.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "maxwell-pilliner.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[62], mcconnaughey.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "mcconnaughey.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[63], mean.manhattan.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "mean-manhattan.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[64], michael.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "michael.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[65], minowski.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "minowski.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[66], mountford.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "mountford-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[67], mountford.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "mountford-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[68], nei.li.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "nei-li.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[69], ochiai.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "ochiai-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[70], ochiai.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "ochiai-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[71], ochiai.e.3)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "ochiai-3.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[72], otsuka.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "otsuka.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[73], pattern.difference.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pattern-difference.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[74], pearson.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pearson-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[75], pearson.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pearson-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[76], pearson.e.3)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pearson-3.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, 
                                funs[77], pearson.heron.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pearson-heron-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, m.n, 
                                funs[78], pearson.heron.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "pearson-heron-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[79], peirce.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "peirce-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[80], peirce.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "peirce-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[81], peirce.e.3)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "peirce-3.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[82], roger.tanimoto.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "roger-tanimoto.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[83], rogot.goldberd.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, " rogot-goldberd.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[84], russel.rao.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "russel-rao.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[85], scott.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "scott.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[86], shape.difference.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "shape-difference.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d,
                                m.n, funs[87], simpson.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "simpson.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[88], size.difference.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "size-difference.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[89], sokal.michener.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-michener.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[90], sokal.sneath.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-1.csv")
    rm(result2)
    
    ################################################################################
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[91], sokal.sneath.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[92], sokal.sneath.e.3)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-3.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[93], sokal.sneath.e.4a)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-4a.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[94], sokal.sneath.e.4b)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-4b.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[95], sokal.sneath.e.5a)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-5a.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[96], sokal.sneath.e.5b)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "sokal-sneath-5b.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[97], sorgenfrei.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "songerfrei.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[98], square.euclidean.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "square-euclidean.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[99], stiles.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "stiles.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[100], tanimoto.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "tanimoto.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[101], tarantula.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "tarantula.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[102], tarwid.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "tarwid.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[103], three.w.jaccard.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "three-w-jaccard.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[104], vari.e)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "vari.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[105], yule.e.1)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "yule-1.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[106], yule.e.2)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "yule-2.csv")
    rm(result2)
    
    ################################################################################
    result2 = compute.measure.2(labels, num.labels, m.a, m.b, m.c, m.d, 
                                m.n, funs[107], yule.e.3)
    result2[which(!is.finite(result2))] <- 0
    write.csv(result2, "yule-3.csv")
    rm(result2)
    
    ##################################################################################################
    # copy file                                                                                      #
    ##################################################################################################
    cat("\nCopy files \n")
    str4 = paste("cp -r ", FolderSplit, "/ ",  folder$FolderRD, sep="")
    print(system(str4))
    
    gc()
    
  } # fim do for each
  
  
  cat("\n# Run: Stop Parallel                                                                              #")
  parallel::stopCluster(cl) 
  gc()
  
} # fim da função



################################################################################
# any errors, please, contact me: elainececiliagatto@gmail.com                 #
################################################################################
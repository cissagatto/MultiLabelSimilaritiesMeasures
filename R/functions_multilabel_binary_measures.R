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
FolderScripts = paste(FolderRoot, "/R", sep="")

##################################################################################################

# BUILD CORRELATION MATRIX ----------------------------------------------------------------
#' Measures for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param labels Label Space from dataset
#' @param num.labels Number labels from Label Space
#' @return correlation matrix num.labels X num.labels
#' @references
#'  
#' @export
#'
#' @examples
#' 
build.matrix.corr <- function(num.labels, labels){
  matrix.corr <- matrix(nrow=num.labels, ncol=num.labels, data=0)
  colnames(matrix.corr) <- colnames(labels)
  rownames(matrix.corr) <- colnames(labels)
  return(matrix.corr)
  gc()
}

# CONTINGENCY TABLE ----------------------------------------------------------------
#' Measures for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param labels Label Space from dataset
#' @param num.labels Number labels from Label Space
#' @return contingence table for all labels
#' @references
#'  
#' @export
#'
#' @examples
#' 
compute.cont.table <- function(labels, num.labels){
  
  retorno = list()
  
  ma <- build.matrix.corr(num.labels, labels)
  mb <- build.matrix.corr(num.labels, labels)
  mc <- build.matrix.corr(num.labels, labels)
  md <- build.matrix.corr(num.labels, labels)
  
  u = (num.labels*num.labels)
  pb <- progress_bar$new(total = u)
  
  #i = 1
  #j = 1
  
  for (i in 1:num.labels){
    
    for (j in 1:num.labels){
      
      x = labels[,i]
      y = labels[,j]
      
      x = as.numeric(levels(x))[x]
      y = as.numeric(levels(y))[y]
      
      ma[i,j] = compute.a(x,y)
      mb[i,j] = compute.b(x,y)
      mc[i,j] = compute.c(x,y)
      md[i,j] = compute.d(x,y)
      
      pb$tick()
      Sys.sleep(1/u)
      gc()
    } # end intern for
    
    gc()
  } # enf extern for    
  
  retorno$ma = ma
  retorno$mb = mb
  retorno$mc = mc 
  retorno$md = md
  return(retorno)
  
  gc()
}

# COMPUTE MARGINAL PROBABILITIES ----------------------------------------------------------------
#' Measures for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param labels Label Space from dataset
#' @param num.labels Number labels from Label Space
#' @param res A object from function compute.cont.table(labels, num.labels)
#' @return marginal probabilities for all labels
#' @references
#'  
#' @export
#'
#' @examples
#' 
compute.marg.probs <- function(labels, num.labels, a, b, c, d){
  
  retorno = list()
  
  mab <- build.matrix.corr(num.labels, labels)
  mac <- build.matrix.corr(num.labels, labels)
  mad <- build.matrix.corr(num.labels, labels)
  mbc <- build.matrix.corr(num.labels, labels)
  mbd <- build.matrix.corr(num.labels, labels)
  mcd <- build.matrix.corr(num.labels, labels)
  mn <- build.matrix.corr(num.labels, labels)
  
  u = (num.labels*num.labels)
  pb <- progress_bar$new(total = u)
  
  #i = 1
  #j = 1
  for (i in 1:num.labels){
    for (j in 1:num.labels){
      
      x = a[i,j]
      y = b[i,j]
      mab[i,j] = compute.ab(x,y)
      #cat("\nAB ", mab[i,j])
      
      w = a[i,j]
      v = c[i,j]
      mac[i,j] = compute.ac(w,v)
      #cat("\nAC ", mac[i,j])
      
      e = a[i,j]
      f = d[i,j]
      mad[i,j] = compute.ad(e,f)
      #cat("\nAD ", mad[i,j])
      
      g = b[i,j]
      h = c[i,j]
      mbc[i,j] = compute.bc(g,h)
      #cat("\nBC ", mbc[i,j])
      
      k = b[i,j]
      l = d[i,j]
      mbd[i,j] = compute.bd(k,l)
      #cat("\nBD ", mbd[i,j])
      
      m = c[i,j]
      n = d[i,j]
      mcd[i,j] = compute.cd(m,n)
      #cat("\nCD ", mcd[i,j])
      
      o = a[i,j]
      p = b[i,j]
      q = c[i,j]
      r = d[i,j]
      mn[i,j] = compute.n(o,p,q,r)
      #cat("\nN ", mn[i,j])
      
      pb$tick()
      Sys.sleep(1/u)
      
      #i = i + 1
      gc()
    } # end intern for
    #j = j + 1
    gc()
  } # enf extern for    
  
  retorno$mab = mab
  retorno$mac = mac
  retorno$mad = mad
  retorno$mbc = mbc
  retorno$mbd = mbd
  retorno$mcd = mcd
  retorno$mn = mn
  return(retorno)
  
  gc()
}

# COVARIANCE ----------------------------------------------------------------
#' Measures for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param labels Label Space from dataset
#' @param num.labels Number labels from Label Space
#' @param res A object from function compute.marg.probs(labels, num.labels, res)
#' @return covariance for all labels
#' @references
#'  
#' @export
#'
#' @examples
#' 
compute.covar <- function(labels, num.labels, ad, bc){
  
  mco <- build.matrix.corr(num.labels, labels)
  u = (num.labels*num.labels)
  pb <- progress_bar$new(total = u)
  for (i in 1:num.labels){
    for (j in 1:num.labels){
      x = ad[i,j]
      y = bc[i,j]
      mco[i,j] = covariance(x,y)
      pb$tick()
      Sys.sleep(1/u)
      gc()
    } # end intern for
    gc()
  } # enf extern for  
  return(mco)
  gc()
}


# ALL ----------------------------------------------------------------
#' Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param labels Label Space from dataset
#' @param num.labels Number labels from Label Space
#' @param  a
#' @param  b
#' @param  c
#' @param  d
#' @param  n
#' @return values for all measures
#' @references
#'  
#' @export
#'
#' @examples
#' 
compute.measure.2 <- function(labels, num.labels, a, b, c, d, n, name, FUN){
  
  retorno = list()
  
  m <- build.matrix.corr(num.labels, labels)
  u = (num.labels*num.labels) # tamanho da matriz 
  pb <- progress_bar$new(total = u) # barra de progresso
  
  for (i in 1:num.labels){
    for (j in 1:num.labels){
      x = as.numeric(a[i,j])
      y = as.numeric(b[i,j])
      w = as.numeric(c[i,j])
      z = as.numeric(d[i,j])
      k = as.numeric(n[i,j])
      m[i,j] = FUN(x,y,w,z,k)
      pb$tick()
      Sys.sleep(1/u)
      gc()
    } # end intern for
    gc()
  } # enf extern for  
  
  return(m)
  gc()
}

compute.measure <- function(Fun, ...) {
  Args <- list(...)
  Args <- lapply(Args, function(M) apply(X = M, MARGIN = c(1,2), FUN = as.numeric))
  Fun(Args)
}

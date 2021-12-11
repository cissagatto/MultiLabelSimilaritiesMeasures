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
# proporção de 1s que ambas as variáveis compartilham nas mesmas posições
# correspondências positivas entre x e y: x e y == 1 
compute.a <- function(x, y){
  return(sum(x == 1 & y == 1))
}

# proporção de 0s na primeira variável e 1s na segunda variável nas mesmas posições
# x ausente: x == 0 e y == 1 
compute.b <- function(x, y){
  return(sum(x == 0 & y == 1))
}

# proporção de 1s na primeira variável e 0s na segunda variável nas mesmas posições
# y ausente: x == 1 e y == 0 
compute.c <- function(x, y){
  return(sum(x == 1 & y == 0))
}

# proporção de zeros que ambas as variáveis compartilham
# correspondências positivas entre x e y: x and y == 0
compute.d <- function(x,y,m){
  return(sum(x == 0 & y == 0))
}

# marginal probabilities

# p1 = a + b --> proporção de uns na primeira variável
compute.ab <- function(a, b){
  return(a+b)
}

# p2 = a + c --> proporção de uns na segunda variável
compute.ac <- function(a, c){
  return(a+c)
}

# p3 = a + d --> diagonal (11) (00)
compute.ad <- function(a, d){
  return(a+d)
}

# p4 = b + c --> diagonal (10) (01) 
compute.bc <- function(b, c){
  return(b+c)
}

# p5 = b + d --> proporção de zeros na segunda variável
compute.bd <- function(b, d){
  return(b+d)
}

# p6 = c + d --> proporção de zeros na primeira variável
compute.cd <- function(c, d){
  return(c+d)
}

compute.n <- function(a,b,c,d){
  return(a+b+c+d)
}


# Na teoria da probabilidade, duas variáveis binárias são chamadas 
# de não correlacionadas se elas  compartilha covariância zero, 
# ou seja, ad - bc = 0 
# returns 0 if (ad - bc) != 0
# returns 1 if (ad - bc) == 0
covariance <- function(x,y){
  return((x-y)==0)
}

# Duas variáveis binárias satisfazem a independência estatística se
# ad / bc = 1
# returns 0 if (ad / bc) != 1
# returns 1 if (ad / bc) == 1
independence <- function(x,y){
  return((x/y)==1)
}


################################################################################
# any errors, please, contact me: elainececiliagatto@gmail.com                 #
################################################################################
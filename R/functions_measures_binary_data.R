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

###############################################################################
# Ample ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param d
ample.e.1 <- function(l){
  return(abs((l$a*(l$c+l$d))/(l$c*(l$a+l$b))))
} 

ample.e.2 <- function(a,b,c,d,n){
  return(abs((a*(c+d))/(c*(a+b))))
} 

###############################################################################
# Anderberg ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param d
#' @param n
anderberg.e.1 <- function(l){
  w = max(l$a,l$b) + max(l$c,l$d) + max(l$a,l$c) + max(l$b,l$d) 
  z = max((l$a+l$c),(l$b+l$d)) + max((l$a+l$b),(l$c+l$d))
  x = ((w-z)/(2*l$n))
  return(x)
}

anderberg.e.2 <- function(a,b,c,d,n){
  w = max(a,b) + max(c,d) + max(a,c) + max(b,d) 
  z = max( (a+c), (b+d) ) + max( (a+b) , (c+d) )
  x = ((w-z)/(2*n))
  return(x)
}

###############################################################################
# Baroni Urbani Buser 1 ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param d
baroni.urbani.buser.1.e.1 <- function(l){
  return((sqrt((l$a*l$d)) + l$a)/((sqrt((l$a*l$d))) + l$a + l$b + l$c))
}

baroni.urbani.buser.1.e.2 <- function(a,b,c,d,n){
  return((sqrt((a*d)) + a)/((sqrt((a*d))) + a + b + c))
}

###############################################################################
# Baroni Urbani Buser 2 ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param d
baroni.urbani.buser.2.e.1 <- function(l){
  return((sqrt((l$a*l$d))+l$a-(l$b+l$c))/((sqrt((l$a*l$d))) + l$a + l$b + l$c))
}

baroni.urbani.buser.2.e.2 <- function(a,b,c,d,n){
  return((sqrt((a*d))+a-(b+c))/((sqrt((a*d))) + a + b + c))
}

###############################################################################
# Braun Banquet ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
braun.banquet.e.1 <- function(l){
  return(l$a/max((l$a+l$b),(l$a+l$c)))
}

braun.banquet.e.2 <- function(a,b,c,d,n){
  return(a/max((a+b),(a+c)))
}

###############################################################################
# Bray Curtis ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
bray.curtis.e.1 <- function(l){
  return((l$b+l$c)/((2*l$a)+l$b+l$c))
}

bray.curtis.e.2 <- function(a,b,c,d,n){
  return((b+c)/((2*a)+b+c))
}

###############################################################################
# Canberra ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param b
#' @param c
canberra.e.1 <- function(l){
  return((l$b+l$c)^(2/2))
}

canberra.e.2 <- function(a,b,c,d,n){
  return((b+c)^(2/2))
}

##############################################################################
# Chord ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
chord.e.1 <- function(l){
  return(sqrt(2*(1-((l$a)/(sqrt((l$a+l$b)*(l$a+l$c)))))))
}

chord.e.2 <- function(a,b,c,d,n){
  return(sqrt(2*(1-((a)/(sqrt((a+b)*(a+c)))))))
}

##############################################################################
# Cityblock ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param b
#' @param c
cityblock.e.1 <- function(l){
  return(l$b+l$c)
}

cityblock.e.2 <- function(a,b,c,d,n){
  return(b+c)
}

##############################################################################
# Cole ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param d
cole.e.1 <- function(l){
  d1 = sqrt(2) * ((l$a*l$d)-(l$b*l$c))
  d2 = (((l$a*l$d)-(l$b*l$c))^2) - ((l$a+l$b)*(l$a+l$c)*(l$b+l$d)*(l$c+l$d))
  d3 = sqrt(d2)
  d4 = d1/d3
  return(d4)
}

cole.e.2 <- function(a,b,c,d,n){
  d1 = sqrt(2) * ((a*d)-(b*c))
  d2 = (((a*d)-(b*c))^2) - ((a+b)*(a+c)*(b+d)*(c+d))
  d3 = sqrt(d2)
  d4 = d1/d3
  return(d4)
}

##############################################################################
# Cosine ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
cosine.e.1 <- function(l){
  return(l$a / sqrt((((l$a+l$b) * (l$a+l$c))))^2)
}

cosine.e.2 <- function(a,b,c,d,n){
  return(a / sqrt((((a+b) * (a+c))))^2)
}

##############################################################################
# Czekanowski ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
czekanowski.e.1 <- function(l){
  return((2*l$a)/((2*l$a)+l$b+l$c))
}

czekanowski.e.2 <- function(a,b,c,d,n){
  return((2*a)/((2*a)+b+c))
}

##############################################################################
# Dennis ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param d
#' @param n
dennis.e.1 <- function(l){
  return(((l$a*l$d)-(l$b*l$c))/sqrt(l$n*(l$a+l$b)*(l$a+l$c)))
}

dennis.e.2 <- function(a,b,c,d,n){
  return(((a*d)-(b*c))/sqrt(n*(a+b)*(a+c)))
}

##############################################################################
# Dice ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
dice.e.1 <- function(l){
  return((2*l$a)/((2*l$a)+l$b+l$c))
}

dice.e.2 <- function(a,b,c,d,n){
  return((2*a)/((2*a)+b+c))
}

##############################################################################
# Dispersion ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param d
disperson.e.1 <- function(l){
  return(((l$a+l$d)-(l$b+l$c))/((l$a+l$b+l$c+l$d)^2))
}

disperson.e.2 <- function(a,b,c,d,n){
  return(((a+d)-(b+c))/((a+b+c+d)^2))
}

##############################################################################
#  Driver Kroeber ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param d
driver.kroeber.e.1 <- function(l){
  return((l$a/2) * ((1/(l$a+l$b)) + (1/(l$a+l$c))))
}

driver.kroeber.e.2 <- function(a,b,c,d,n){
  return((a/2) * ((1/(a+b)) + (1/(a+c))))
}

##############################################################################
# Euclidean ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param b
#' @param c
euclidean.e.1 <- function(l){
  return(sqrt(l$b+l$c))
}

euclidean.e.2 <- function(a,b,c,d,n){
  return(sqrt(b+c))
}

##############################################################################
# Eyraud ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param d
#' @param n
eyraud.e.1 <- function(l){
  return((l$n^2)*((l$n*l$a)-(l$a+l$b)*(l$a+l$c))/(l$a+l$b)*(l$a+l$c)*(l$b+l$d)*(l$c+l$d))
}

eyraud.e.2 <- function(a,b,c,d,n){
  return((n^2)*((n*a)-(a+b)*(a+c))/(a+b)*(a+c)*(b+d)*(c+d))
}

##############################################################################
# Fager Mcgowan ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
fager.mcgowan.e.1 <- function(l){
  return((l$a/sqrt((l$a+l$b)+(l$a+l$c))) - (max((l$a+l$b),(l$a+l$c))/2))
}

fager.mcgowan.e.2 <- function(a,b,c,d,n){
  return((a/sqrt((a+b)+(a+c))) - (max((a+b),(a+c))/2))
}

##############################################################################
# Faith ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param d
faith.e.1 <- function(l){
  return((l$a+(0.5*l$d))/(l$a+l$b+l$c+l$d))
}

faith.e.2 <- function(a,b,c,d,n){
  return((a+(0.5*d))/(a+b+c+d))
}

##############################################################################
# Forbes ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
#' @param n
forbes.2.e.1 <- function(l){
  return((l$n*l$a)-((l$a+l$b)*(l$a+l$c))/l$n*(min((l$a+l$b),(l$a+l$c))-((l$a+l$b)*(l$a+l$c))))
}

forbes.2.e.2 <- function(a,b,c,d,n){
  return((n*a)-((a+b)*(a+c))/n*(min((a+b),(a+c))-((a+b)*(a+c))))
}

##############################################################################
# Forbesi ----------------------------------------------------------------
#' Similaritie Measure for categorial data
#'
#' @family Multi_Label_Binary_Measures
#' @param a
#' @param b
#' @param c
forbesi.e.1 <- function(l){
  return((l$n*l$a)/((l$a+l$b)*(l$a+l$c)))
}

forbesi.e.2 <- function(a,b,c,d,n){
  return((n*a)/((a+b)*(a+c)))
}

###############################################################################
# Fossum
fossum.e.1 <- function(l){
  return((l$n*(l$a-0.5)^2)/(l$a+l$b)*(l$a+l$c))
}

fossum.e.2 <- function(a,b,c,d,n){
  return((n*(a-0.5)^2)/(a+b)*(a+c))
}

###############################################################################
# Gilbert Well
gilbert.well.e.1 <- function(l){
  return(log(l$a)-log(l$n)-log((l$a+l$b)/l$n)-log((l$a+l$c)/l$n))
}

gilbert.well.e.2 <- function(a,b,c,d,n){
  return(log(a) - log(n) - log((a+b)/n) - log((a+c)/n))
}

###############################################################################
# Goodman Kruskal
goodman.kruskal.e.1 <- function(l){
  z = max(l$a,l$b)+max(l$c,l$d)+max(l$a,l$c)+max(l$b,l$c)
  w = max((l$a+l$c),(l$b+l$d))+max((l$a+l$b),(l$c+l$d))
  p = (z-w)/((2*l$n)-w)
  return(p)
}

goodman.kruskal.e.2 <- function(a,b,c,d,n){
  z = max(a,b)+max(c,d)+max(a,c)+max(b,c)
  w = max((a+c),(b+d))+max((a+b),(c+d))
  p = (z-w)/((2*n)-w)
  return(p)
}

###############################################################################
# Gower
gower.e.1 <- function(l){
  return((l$a+l$d)/sqrt((l$a+l$b)*(l$a+l$c)*(l$b+l$d)*(l$c+l$d)))
}

gower.e.2 <- function(a,b,c,d,n){
  return((a+d)/sqrt((a+b)*(a+c)*(b+d)*(c+d)))
}

###############################################################################
# Gower Legendre
gower.legendre.e.1 <- function(l){
  return((l$a+l$d)/(l$a+(0.5*(l$b+l$c))+l$d))
}

gower.legendre.e.2 <- function(a,b,c,d,n){
  return((a+d)/(a+(0.5*(b+c))+d))
}

###############################################################################
# Hamann
hamann.e.1 <- function(l){
  return(((l$a+l$d)-(l$b+l$c))/(l$a+l$b+l$c+l$d))
}

hamann.e.2 <- function(a,b,c,d,n){
  return(((a+d)-(b+c))/(a+b+c+d))
}

###############################################################################
# Hamming 
hamming.e.1 <- function(l){
  return(l$b+l$c)
}

hamming.e.2 <- function(a,b,c,d,n){
  return(b+c)
}

###############################################################################
# Helllinger
hellinger.e.1 <- function(l){
  return(2 * sqrt(1 - ((l$a)/(sqrt((l$a+l$b)*(l$a+l$c))))))
}

hellinger.e.2 <- function(a,b,c,d,n){
  return(2 * sqrt(1 - ((a)/(sqrt((a+b)*(a+c))))))
}

###############################################################################
# Inner Product
inner.product.e.1 <- function(l){
  return(l$a+l$d)
}

inner.product.e.2 <- function(a,b,c,d,n){
  return(a+d)
}

###############################################################################
# Intersection
intersection.e.1 <- function(l){
  return(l$a)
}

intersection.e.2 <- function(a,b,c,d,n){
  return(a)
}

###############################################################################
# Jaccard
jaccard.e.1 <- function(l){
  return(l$a/(l$a+l$b+l$c))
}

jaccard.e.2 <- function(a,b,c,d,n){
  return(a/(a+b+c))
}

###############################################################################
# Jonhson
jonhson.e.1 <- function(l){
  return((l$a/(l$a+l$b)) + (l$a/(l$a+l$c)))
}

jonhson.e.2 <- function(a,b,c,d,n){
  return((a/(a+b)) + (a/(a+c)))
}

###############################################################################
# kulczynski 1
kulczynski.1.e.1 <- function(l){
  return(l$a/(l$b+l$c))
}

kulczynski.1.e.1 <- function(a,b,c,d,n){
  return(a/(b+c))
}

###############################################################################
# kulczynski 2
kulczynski.2.e.1 <- function(l){
  return((l$a/2) * ((2*l$a)+l$b+l$c))/(l$a+l$b)*(l$a+l$c)
}

kulczynski.2.e.2 <- function(a,b,c,d,n){
  return((a/2) * ((2*a)+b+c))/(a+b)*(a+c)
}

###############################################################################
# Lance Williams
lance.williams.e.1 <- function(l){
  return((l$b+l$c)/((2*l$a)+l$b+l$c))
}

lance.williams.e.2 <- function(a,b,c,d,n){
  return((b+c)/((2*a)+b+c))
}

###############################################################################
# Manhattan 
manhattan.e.1 <- function(l){
  return(l$b+l$c)
}

manhattan.e.2 <- function(a,b,c,d,n){
  return(b+c)
}

###############################################################################
# Mcconnaughey
mcconnaughey.e.1 <- function(l){
  return(((l$a^2)-(l$b-l$c))/(l$a+l$b)*(l$a+l$c))
}

mcconnaughey.e.2 <- function(a,b,c,d,n){
  return(((a^2) - (b-c))/(a+b)*(a+c))
}

###############################################################################
# Mean Manhattan
mean.manhattan.e.1 <- function(l){
  return((l$b+l$c)/(l$a+l$b+l$c+l$d))
}

mean.manhattan.e.2 <- function(a,b,c,d,n){
  return((b+c)/(a+b+c+d))
}

###############################################################################
# Michael
michael.e.1 <- function(l){
  return(4*((l$a*l$d)-(l$b*l$c))/((l$a+l$d)^2)+((l$b+l$c)^2))
}

michael.e.2 <- function(a,b,c,d,n){
  return(4*((a*d)-(b*c))/((a+d)^2)  + ((b+c)^2))
}

###############################################################################
# Minowski
minowski.e.1 <- function(l){
  return((l$b+l$c)^(1/1))
}

minowski.e.2 <- function(a,b,c,d,n){
  return((b+c)^(1/1))
}

###############################################################################
# Mountford
mountford.e.1 <- function(l){
  return(l$a/(0.5*(((l$a*l$b)+(l$a*l$c))+(l$b*l$c))))
}

mountford.e.2 <- function(a,b,c,d,n){
  return(a/(0.5*(((a*b)+(a*c))+(b*c))))
}

###############################################################################
# Mean Manhattan
nei.li.e.1 <- function(l){
  return((2*l$a)/((l$a+l$b)+(l$a+l$c)))
}

nei.li.e.2 <- function(a,b,c,d,n){
  return((2*a)/((a+b)+(a+c)))
}

###############################################################################
# Ochiai 2
ochiai.2.e.1 <- function(l){
  return((l$a*l$d)/sqrt((l$a+l$b)*(l$a+l$c)*(l$b+l$d)*(l$c+l$d)))
}

ochiai.2.e.1 <- function(a,b,c,d,n){
  return((a*d)/sqrt((a+b)*(a+c)*(b+d)*(c+d)))
}

###############################################################################
# Otsuka
otsuka.e.1 <- function(l){
  return(l$a/((l$a+l$b)*(l$a+l$c))^0.5)
}

otsuka.e.2 <- function(a,b,c,d,n){
  return(a/((a+b)*(a+c))^0.5)
}

###############################################################################
# Pattern Difference
pattern.difference.e.1 <- function(l){
  return((4*l$b*l$c)/((l$a+l$b+l$c+l$d)^2))
}

pattern.difference.e.2 <- function(a,b,c,d,n){
  return((4*b*c)/((a+b+c+d)^2))
}

###############################################################################
# Pearson 
pearson.1.e.1 <- function(l){
  return(l$n*(((l$a*l$d)-(l$b*l$c))^2)/(l$a+l$b)*(l$a+l$c)*(l$c+l$d)*(l$b+l$d))
}

pearson.1.e.2 <- function(a,b,c,d,n){
  return(n*(((a*d)-(b*c))^2)/(a+b)*(a+c)*(c+d)*(b+d))
}

pearson.2.e.1 <- function(l){
  z = pearson.1.e(l$a,l$b,l$c,l$d,l$n)
  w = (z/(l$n*z))^(1/2)
  return(w)
}

pearson.2.e.2 <- function(a,b,c,d,n){
  z = pearson.1.e(a,b,c,d,n)
  w = (z/(n*z))^(1/2)
  return(w)
}

pearson.3.e.1 <- function(l){
  z = ((l$a*l$d)-(l$b*l$c))/sqrt((l$a+l$b)*(l$a+l$c)*(l$b+l$d)*(l$c+l$d))
  w = (z/(l$n+z))^(1/2)
  return(w)
}

pearson.3.e.2 <- function(a,b,c,d,n){
  z = ((a*d)-(b*c))/sqrt((a+b)*(a+c)*(b+d)*(c+d))
  w = (z/(n+z))^(1/2)
  return(w)
}

pearson.heron.1.e.1 <- function(l){
  return((l$a*l$d)-(l$b*l$c)/(l$a+l$b)*(l$a+l$c)*(l$b+l$d)*(l$c+d))
}

pearson.heron.1.e.2 <- function(a,b,c,d,n){
  return((a*d)-(b*c)/(a+b)*(a+c)*(b+d)*(c+d))
}

pearson.heron.2.e.1 <- function(l){
  return(cos(pi*(sqrt(l$b*l$c))/sqrt(l$a*l$d)+sqrt(l$b*l$c)))
}

pearson.heron.2.e.2 <- function(a,b,c,d,n){
  return(cos(pi*(sqrt(b*c))/sqrt(a*d)+sqrt(b*c)))
}

###############################################################################
# Pearson 
peirce.e <- function(a,b,c,d,n){
  return((a*b)+(b*c)/((a*b)+(2*b*c)+(c*d)))
}

###############################################################################
# Roger Tanimoto 
roger.tanimoto.e.1 <- function(l){
  return((l$a+l$d)/(2*(l$b+l$c)+l$d))
}

roger.tanimoto.e.2 <- function(a,b,c,d,n){
  return((a+d)/(2*(b+c)+d))
}

###############################################################################
# Russel Rao 
russel.rao.e.1 <- function(l){
  return(l$a/(l$a+l$b+l$c+l$d))
}

russel.rao.e.2 <- function(a,b,c,d,n){
  return(a/(a+b+c+d))
}

###############################################################################
# Shape Difference 
shape.differnece.e.1 <- function(l){
  return((l$n*(l$b+l$c)-((l$b-l$c)^2))/((l$a+l$b+l$c+l$d)^2))
}

shape.differnece.e.2 <- function(a,b,c,d,n){
  return((n*(b+c)-((b-c)^2))/((a+b+c+d)^2))
}

###############################################################################
# Simpson
simpson.e.1 <- function(l){
  return(l$a/min((l$a+l$b),(l$a+l$c)))
}

simpson.e.2 <- function(a,b,c,d,n){
  return(a/min((a+b),(a+c)))
}

###############################################################################
# Size Difference 
size.difference.e.1 <- function(l){
  return(((l$b+l$c)^2)/((l$a+l$b+l$c+l$d)^2))
}

size.difference.e.2 <- function(a,b,c,d,n){
  return(((b+c)^2)/((a+b+c+d)^2))
}

###############################################################################
# Sokal 
sokal.michener.e.1 <- function(l){
  return((l$a+l$d)/(l$a+l$b+l$c+l$d))
}

sokal.michener.e.2 <- function(a,b,c,d,n){
  return((a+d)/(a+b+c+d))
}

sokal.sneath.1.e.1 <- function(l){
  return(l$a/(l$a+(2*l$b)+(2*l$c)))
}

sokal.sneath.1.e.2 <- function(a,b,c,d,n){
  return(a/(a+(2*b)+(2*c)))
}

sokal.sneath.2.e.1 <- function(l){
  return(2*(l$a+l$d)/((2*l$a)+l$b+l$c+(2*l$d)))
}

sokal.sneath.2.e.2 <- function(a,b,c,d,n){
  return(2*(a+d)/((2*a)+b+c+(2*d)))
}

sokal.sneath.3.e.1 <- function(l){
  return((l$a+l$d)/(l$b+l$c))
}

sokal.sneath.3.e.2 <- function(a,b,c,d,n){
  return((a+d)/(b+c))
}

sokal.sneath.4.e.1 <- function(l){
  return(((l$a/(l$a+l$b))+(l$a/(l$a+l$c))+(l$d/(l$b+l$d))+(l$d/(l$b+l$d)))/4)
}

sokal.sneath.4.e.2 <- function(a,b,c,d,n){
  return(((a/(a+b))+(a/(a+c))+(d/(b+d))+(d/(b+d)))/4)
}

sokal.sneath.5.e.1 <- function(l){
  return((l$a*l$d)/(l$a+l$b)*(l$a+l$c)*(l$b+l$d)*((l$c*l$d)^0.5))
}

sokal.sneath.5.e.2 <- function(a,b,c,d,n){
  return((a*d)/(a+b)*(a+c)*(b+d)*((c*d)^0.5))
}

###############################################################################
# Sorgenfrei
sorgenfrei.e.1 <- function(l){
  return((l$a^2)/(l$a+l$b)*(l$a+l$c))
}

sorgenfrei.e.2 <- function(a,b,c,d,n){
  return((a^2)/(a+b)*(a+c))
}

###############################################################################
# Square Euclidean
square.euclidean.e.1 <- function(l){
  return(sqrt((l$b+l$c)^2))
}

square.euclidean.e.2 <- function(a,b,c,d,n){
  return(sqrt((b+c)^2))
}

###############################################################################
# Stiles
stiles.e.1 <- function(l){
  return(log10((l$n*(abs((l$a*l$d)-(l$b*l$c))-(l$n/2))^2) /((l$a+l$b)*(l$a+l$c)*(l$b+l$d)*(l$c+l$d))))
}

stiles.e.2 <- function(a,b,c,d,n){
  return(log10((n*(abs((a*d)-(b*c))-(n/2))^2) /((a+b)*(a+c)*(b+d)*(c+d))))
}

###############################################################################
# Tanimoto
tanimoto.e.1 <- function(l){
  return(l$a/((l$a+l$b)+(l$a+l$c)-l$a))
}

tanimoto.e.2 <- function(a,b,c,d,n){
  return(a/((a+b)+(a+c)-a))
}

###############################################################################
# Tarantula
tarantula.e.1 <- function(l){
  return((l$a*(l$c+l$d))/(l$c*(l$a+l$b)))
}

tarantula.e.2 <- function(a,b,c,d,n){
  return((a*(c+d))/(c*(a+b)))
}

###############################################################################
# Tarwid
tarwid.e.1 <- function(l){
  return((l$n*l$a)-(l$a+l$b)*(l$a+l$c)/(l$n*l$a)+(l$a+l$b)*(l$a+l$c))
}

tarwid.e.2 <- function(a,b,c,d,n){
  return((n*a) - (a+b) * (a+c) / (n*a) + (a+b) * (a+c))
}

###############################################################################
# 3W Jaccard
three.w.jaccard.e.1 <- function(l){
  return((3*l$a)/((3*l$a)+l$b+l$c))
}

three.w.jaccard.e.2 <- function(a,b,c,d,n){
  return((3*a)/((3*a)+b+c))
}

###############################################################################
# Vari
vari.e.1 <- function(l){
  return((l$b+l$c)/(4*(l$a+l$b+l$c+l$d)))
}

vari.e.2 <- function(a,b,c,d,n){
  return((b+c)/(4*(a+b+c+d)))
}

###############################################################################
# Yuleq
yuleq.1.e.1 <- function(l){
  return(((l$a*l$d)-(l$b*l$c))/((l$a*l$d)+(l$b*l$c)))
}

yuleq.1.e.2 <- function(a,b,c,d,n){
  return(((a*d)-(b*c))/((a*d)+(b*c)))
}

yuleq.2.e.1 <- function(l){
  return((2*l$b*l$c)/(l$a*l$d)+(l$b*l$c))
}

yuleq.2.e.2 <- function(a,b,c,d,n){
  return((2*b*c)/(a*d) + (b*c))
}

yule.w.e.1 <- function(l){
  return(sqrt((l$a*l$d))-sqrt((l$b*l$c))/sqrt((l$a*l$d))+sqrt((l$b*l$c)))
}

yule.w.e.2 <- function(a,b,c,d,n){
  return(sqrt((a*d))-sqrt((b*c))/sqrt((a*d))+sqrt((b*c)))
}


###############################################################################
# get name list
getNamesListFunctions <- function(){
  retorno = list()
  names_function = c("ample.e",
                   "anderberg.e ",
                   "baroni.urbani.buser.1.e",
                   "baroni.urbani.buser.2.e",
                   "braun.banquet.e",
                   "bray.curtis.e",
                   "canberra.e",
                   "chord.e ",
                   "cityblock.e",
                   "cole.e",
                   "cosine.e",
                   "czekanowski.e ",
                   "dennis.e",
                   "dice.e",
                   "disperson.e",
                   "driver.kroeber.e",
                   "euclidean.e ",
                   "eyraud.e",
                   "fager.mcgowan.e",
                   "faith.e ",
                   "forbes.2.e", 
                   "forbesi.e",
                   "fossum.e ",
                   "gilbert.well.e",
                   "goodman.kruskal.e",
                   "gower.e ",
                   "gower.legendre.e ",
                   "hamann.e ",
                   "hamming.e",
                   "hellinger.e",
                   "inner.product.e",
                   "intersection.e",
                   "jaccard.e ",
                   "johnson.e",
                   "kulczynski.1.e",
                   "kulczynski.2.e",
                   "lance.williams.e",
                   "manhattan.e",
                   "mcconnaughey.e",
                   "mean.manhattan.e",
                   "michael.e",
                   "minowski.e ",
                   "mountford.e ",
                   "nei.li.e",
                   "ochiai.2.e",
                   "otsuka.e",
                   "pattern.difference.e",
                   "pearson.1.e",
                   "pearson.2.e",
                   "pearson.3.e",
                   "pearson.heron.1.e ",
                   "pearson.heron.2.e ",
                   "peirce.e",
                   "roger.tanimoto.e",
                   "russel.rao.e ",
                   "shape.differnece.e",
                   "simpson.e",
                   "size.difference.e ",
                   "sokal.michener.e",
                   "sokal.sneath.1.e ",
                   "sokal.sneath.2.e ",
                   "sokal.sneath.3.e",
                   "sokal.sneath.4.e ",
                   "sokal.sneath.5.e ",
                   "sorgenfrei.e",
                   "square.euclidean.e",
                   "stiles.e",
                   "tanimoto.e",
                   "tarantula.e",
                   "tarwid.e ",
                   "tw.jaccard.e ",
                   "vari.e",
                   "yule.w.e ",
                   "yuleq.1.e",
                   "yuleq.2.e")
  retorno$nomes_funcoes = names_function
  retorno$tamanho = length(names_function)
  return(retorno)
  gc()
}

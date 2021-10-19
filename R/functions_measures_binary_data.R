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
FolderScripts = paste(FolderRoot, "/scripts", sep="")

##################################################################################################
ample.e <- function(a,b,c,d,n){
  return(abs((a*(c+d))/(c*(a+b))))
}

anderberg.e <- function(a,b,c,d,n){
  z = max(a,b)+max(c,d)+max(a,c)+max(b,c)
  w = max((a+c),(b+d))+max((a+b),(c+d))
  p = (z-w)/(2*n)
  return(p)
}

baroni.urbani.buser.1.e <- function(a,b,c,d,n){
  return((sqrt((a*d)) + a)/((sqrt((a*d))) + a + b + c))
}

baroni.urbani.buser.2.e <- function(a,b,c,d,n){
  return((sqrt((a*d))+a-(b+c))/((sqrt((a*d))) + a + b + c))
}

braun.banquet.e <- function(a,b,c,d,n){
  return(a/max((a+b),(a+c)))
}

bray.curtis.e <- function(a,b,c,d,n){
  return((b+c)/((2*a)+b+c))
}

canberra.e <- function(a,b,c,d,n){
  return((b+c)^(2/2))
}

chord.e <- function(a,b,c,d,n){
  return(sqrt(2*(1-((a)/(sqrt((a+b)*(a+c)))))))
}

cityblock.e <- function(a,b,c,d,n){
  return(b+c)
}

cole.e <- function(a,b,c,d,n){
  d1 = sqrt(2) * ((a*d)-(b*c))
  d2 = (((a*d)-(b*c))^2) - ((a+b)*(a+c)*(b+d)*(c+d))
  d3 = sqrt(d2)
  d4 = d1/d3
  return(d4)
}

cosine.e <- function(a,b,c,d,n){
  return(a / sqrt((((a+b) * (a+c))))^2)
}

czekanowski.e <- function(a,b,c,d,n){
  return((2*a)/((2*a)+b+c))
}

dennis.e <- function(a,b,c,d,n){
  return(((a*d)-(b*c))/sqrt(n*(a+b)*(a+c)))
}

dice.e <- function(a,b,c,d,n){
  return((2*a)/((2*a)+b+c))
}

disperson.e <- function(a,b,c,d,n){
  return(((a+d)-(b+c))/((a+b+c+d)^2))
}

driver.kroeber.e <- function(a,b,c,d,n){
  return((a/2) * ((1/(a+b)) + (1/(a+c))))
}

euclidean.e <- function(a,b,c,d,n){
  return(sqrt(b+c))
}

eyraud.e <- function(a,b,c,d,n){
  return((n^2)*((n*a)-(a+b)*(a+c))/(a+b)*(a+c)*(b+d)*(c+d))
}

fager.mcgowan.e <- function(a,b,c,d,n){
  return((a/sqrt((a+b)+(a+c))) - (max((a+b),(a+c))/2))
}

faith.e <- function(a,b,c,d,n){
  return((a+(0.5*d))/(a+b+c+d))
}

forbes.2.e <- function(a,b,c,d,n){
  return((n*a)-((a+b)*(a+c))/n*(min((a+b),(a+c))-((a+b)*(a+c))))
}

forbesi.e <- function(a,b,c,d,n){
  return((n*a)/((a+b)*(a+c)))
}

fossum.e <- function(a,b,c,d,n){
  return((n*(a-0.5)^2)/(a+b)*(a+c))
}

gilbert.well.e <- function(a,b,c,d,n){
  return(log(a) - log(n) - log((a+b)/n) - log((a+c)/n))
}

goodman.kruskal.e <- function(a,b,c,d,n){
  z = max(a,b)+max(c,d)+max(a,c)+max(b,c)
  w = max((a+c),(b+d))+max((a+b),(c+d))
  p = (z-w)/((2*n)-w)
  return(p)
}

gower.e <- function(a,b,c,d,n){
  return((a+d)/sqrt((a+b)*(a+c)*(b+d)*(c+d)))
}

gower.legendre.e <- function(a,b,c,d,n){
  return((a+d)/(a+(0.5*(b+c))+d))
}

hamann.e <- function(a,b,c,d,n){
  return(((a+d)-(b+c))/(a+b+c+d))
}

hamming.e <- function(a,b,c,d,n){
  return(b+c)
}

hellinger.e <- function(a,b,c,d,n){
  return(2 * sqrt(1 - ((a)/(sqrt((a+b)*(a+c))))))
}

inner.product.e <- function(a,b,c,d,n){
  return(a+d)
}

intersection.e <- function(a,b,c,d,n){
  return(a)
}

jaccard.e <- function(a,b,c,d,n){
  return(a/(a+b+c))
}

johnson.e <- function(a,b,c,d,n){
  return((a/(a+b)) + (a/(a+c)))
}

kulczynski.1.e <- function(a,b,c,d,n){
  return(a/(b+c))
}

kulczynski.2.e <- function(a,b,c,d,n){
  return((a/2) * ((2*a)+b+c))/(a+b)*(a+c)
}

lance.williams.e <- function(a,b,c,d,n){
  return((b+c)/((2*a)+b+c))
}

manhattan.e <- function(a,b,c,d,n){
  return(b+c)
}

mcconnaughey.e <- function(a,b,c,d,n){
  return(((a^2) - (b-c))/(a+b)*(a+c))
}

mean.manhattan.e <- function(a,b,c,d,n){
  return((b+c)/(a+b+c+d))
}

michael.e <- function(a,b,c,d,n){
  return(4*((a*d)-(b*c))/((a+d)^2)  + ((b+c)^2))
}

minowski.e <- function(a,b,c,d,n){
  return((b+c)^(1/1))
}

mountford.e <- function(a,b,c,d,n){
  return(a/(0.5*(((a*b)+(a*c))+(b*c))))
}

nei.li.e <- function(a,b,c,d,n){
  return((2*a)/((a+b)+(a+c)))
}

ochiai.2.e <- function(a,b,c,d,n){
  return((a*d)/sqrt((a+b)*(a+c)*(b+d)*(c+d)))
}

otsuka.e <- function(a,b,c,d,n){
  return(a/((a+b)*(a+c))^0.5)
}

pattern.difference.e <- function(a,b,c,d,n){
  return((4*b*c)/((a+b+c+d)^2))
}

pearson.1.e <- function(a,b,c,d,n){
  return(n*(((a*d)-(b*c))^2)/(a+b)*(a+c)*(c+d)*(b+d))
}

pearson.2.e <- function(a,b,c,d,n){
  z = pearson.1.e(a,b,c,d,n)
  w = (z/(n*z))^(1/2)
  return(w)
}

pearson.3.e <- function(a,b,c,d,n){
  z = ((a*d)-(b*c))/sqrt((a+b)*(a+c)*(b+d)*(c+d))
  w = (z/(n+p))^(1/2)
  return(w)
}

pearson.heron.1.e <- function(a,b,c,d,n){
  return((a*d)-(b*c)/(a+b)*(a+c)*(b+d)*(c+d))
}

pearson.heron.2.e <- function(a,b,c,d,n){
  return(cos(pi*(sqrt(b*c))/sqrt(a*d)+sqrt(b*c)))
}

roger.tanimoto.e <- function(a,b,c,d,n){
  return((a+d)/(2*(b+c)+d))
}

russel.rao.e <- function(a,b,c,d,n){
  return(a/(a+b+c+d))
}

shape.differnece.e <- function(a,b,c,d,n){
  return((n*(b+c)-((b-c)^2))/((a+b+c+d)^2))
}

simpson.e <- function(a,b,c,d,n){
  return(a/min((a+b),(a+c)))
}

size.difference.e <- function(a,b,c,d,n){
  return(((b+c)^2)/((a+b+c+d)^2))
}

sokal.michener.e <- function(a,b,c,d,n){
  return((a+d)/(a+b+c+d))
}

sokal.sneath.1.e <- function(a,b,c,d,n){
  return(a/(a+(2*b)+(2*c)))
}

sokal.sneath.2.e <- function(a,b,c,d){
  return(2*(a+d)/((2*a)+b+c+(2*d)))
}

sokal.sneath.3.e <- function(a,b,c,d,n){
  return((a+d)/(b+c))
}

sokal.sneath.4.e <- function(a,b,c,d,n){
  return(((a/(a+b))+(a/(a+c))+(d/(b+d))+(d/(b+d)))/4)
}

sokal.sneath.5.e <- function(a,b,c,d,n){
  return((a*d)/(a+b)*(a+c)*(b+d)*((c*d)^0.5))
}

sorgenfrei.e <- function(a,b,c,d,n){
  return((a^2)/(a+b)*(a+c))
}

square.euclidean.e <- function(a,b,c,d,n){
  return(sqrt((b+c)^2))
}

stiles.e <- function(a,b,c,d,n){
  return(log10((n*(abs((a*d)-(b*c))-(n/2))^2) /((a+b)*(a+c)*(b+d)*(c+d))))
}

tanimoto.e <- function(a,b,c,d,n){
  return(a/((a+b)+(a+c)-a))
}

tarantula.e <- function(a,b,c,d,n){
  return((a*(c+d))/(c*(a+b)))
}

tarwid.e <- function(a,b,c,d,n){
  return((n*a) - (a+b) * (a+c) / (n*a) + (a+b) * (a+c))
}

three.w.jaccard.e <- function(a,b,c,d,n){
  return((3*a)/((3*a)+b+c))
}

vari.e <- function(a,b,c,d,n){
  return((b+c)/(4*(a+b+c+d)))
}

yuleq.1.e <- function(a,b,c,d,n){
  return(((a*d)-(b*c))/((a*d)+(b*c)))
}

yuleq.2.e <- function(a,b,c,d,n){
  return((2*b*c)/(a*d) + (b*c))
}

yule.w.e <- function(a,b,c,d,n){
  return(sqrt((a*d))-sqrt((b*c))/sqrt((a*d))+sqrt((b*c)))
}

# get name list
getNamesListFunctions <- function(){
  names_function = list("ample.e ",
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
                   "ochiai.e ",
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
  return(names_function)
}

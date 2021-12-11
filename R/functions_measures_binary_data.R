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
    d1 = l$a*(l$c+l$d)
    d2 = l$c*(l$a+l$b)
    d3 = d1/d2
    d4 = abs(d3)
    return(d4)
  } 
  
  ample.e <- function(a,b,c,d,n){
    d1 = a*(c+d)
    d2 = c*(a+b)
    d3 = d1/d2
    d4 = abs(d3)
    return(d4)
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
    delta_1 = max(l$a,l$b) + max(l$c,l$d) + max(l$a,l$c) + max(l$b,l$d) 
    delta_2 = max((l$a+l$c),(l$b+l$d)) + max((l$a+l$b),(l$c+l$d))
    delta = ((delta_1-delta_2)/(2*l$n))
    return(delta)
  }
  
  anderberg.e <- function(a,b,c,d,n){
    delta_1 = max(a,b) + max(c,d) + max(a,c) + max(b,d) 
    delta_2 = max((a+c),(b+d)) + max((a+b),(c+d))
    delta = ((delta_1-delta_2)/(2*n))
    return(delta)
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
  baroni.urbani.buser.e.11 <- function(l){
    d1 = sqrt(abs(l$a*l$d)) + l$a
    d2 = sqrt(abs(l$a*l$d)) + a + l$b + l$c
    d3 = d1/d2
    return(d3)
  }
  
  baroni.urbani.buser.e.1 <- function(a,b,c,d,n){
    d1 = sqrt(abs(a*d)) + a
    d2 = sqrt(abs(a*d)) + a + b + c
    d3 = d1/d2
    return(d3)
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
  baroni.urbani.buser.e.22 <- function(l){
    d1 = l$a - l$b - l$c + sqrt(abs(l$a*l$d)) 
    d2 = l$a + l$b + l$c + sqrt(abs(l$a*l$d)) 
    d3 = d1/d2
    return(d3)
  }
  
  baroni.urbani.buser.e.2 <- function(a,b,c,d,n){
    d1 = a - b - c + sqrt(abs(a*d)) 
    d2 = a + b + c + sqrt(abs(a*d)) 
    d3 = d1/d2
    return(d3)
  }
  
  ###############################################################################
  # Braun Banquet ----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  braun.blanquet.e.1 <- function(l){
    return(l$a/max((l$a+l$b),(l$a+l$c)))
  }
  
  braun.blanquet.e <- function(a,b,c,d,n){
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
  
  bray.curtis.e <- function(a,b,c,d,n){
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
  
  canberra.e <- function(a,b,c,d,n){
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
  cohen.e.1 <- function(l){
    d1 = 2 * ((l$a*l$d) - (l$b*l$c))
    d2 = ((l$a+l$b)*(l$b+l$d)) + ((l$a+l$c)*(l$c+l$d))
    d3 = d1/d2  
    return(d3)
  }
  
  cohen.e <- function(a,b,c,d,n){
    d1 = 2 * ((a*d) - (b*c))
    d2 = ((a+b)*(b+d)) + ((a+c)*(c+d))
    d3 = d1/d2  
    return(d3)
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
    d1 = l$a/sqrt(abs((l$a+l$b)*(l$a+l$c)) ) 
    d2 = 1 - d1
    d3 = sqrt(abs(2 * d2))  
    return(d3)
  }
  
  chord.e <- function(a,b,c,d,n){
    d1 = a/sqrt(abs((a+b)*(a+c)) ) 
    d2 = 1 - d1
    d3 = sqrt(abs(2 * d2))  
    return(d3)
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
  
  cityblock.e <- function(a,b,c,d,n){
    return(b+c)
  }
  
  
  ##############################################################################
  # Cityblock ----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param b
  #' @param c
  clement.e.1 <- function(l){
    d1 = l$a*(l$c+l$d) / (l$a+l$b)
    d2 = l$d*(l$a+l$b) / (l$c+l$d)
    d3 = d1 + d2
    return(d3)
  }
  
  clement.e <- function(a,b,c,d,n){
    d1 = a*(c+d) / (a+b)
    d2 = d*(a+b) / (c+d)
    d3 = d1 + d2
    return(d3)
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
  cole.e.11 <- function(l){
    d1 = sqrt(2) * ((l$a*l$d)-(l$b*l$c))
    d2 = (((l$a*l$d)-(l$b*l$c))^2)
    d3 = ((l$a+l$b)*(l$a+l$c)*(l$b+l$d)*(l$c+l$d))
    d4 = d2 - d3
    d5 = sqrt(abs(d4))
    return(d5)
  }
  
  cole.e.1 <- function(a,b,c,d,n){
    d1 = sqrt(2) * ((a*d)-(b*c))
    d2 = (((a*d)-(b*c))^2)
    d3 = ((a+b)*(a+c)*(b+d)*(c+d))
    d4 = d2 - d3
    d5 = sqrt(abs(d4))
    return(d5)
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
  cole.e.22 <- function(l){
    d1 = (l$a*l$d) - (l$b*l$c)
    d2 = (l$a+l$b) * (l$b+l$d)
    d3 = d1/d2
    return(d3)
  }
  
  cole.e.2 <- function(a,b,c,d,n){
    d1 = (a*d) - (b*c)
    d2 = (a+b) * (b+d)
    d3 = d1/d2
    return(d3)
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
  cole.e.33 <- function(l){
    d1 = (l$a*l$d) - (l$b*l$c)
    d2 = (l$a+l$c) * (l$c+l$d) 
    d3 = d1/d2
    return(d3)
  }
  
  cole.e.3 <- function(a,b,c,d,n){
    d1 = (a*d) - (b*c)
    d2 = (a+c) * (c+d) 
    d3 = d1/d2
    return(d3)
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
    return(l$a / sqrt(abs((((l$a+l$b) * (l$a+l$c))))^2))
  }
  
  cosine.e <- function(a,b,c,d,n){
    return(a/sqrt(abs((((a+b) * (a+c))))^2))
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
  
  czekanowski.e <- function(a,b,c,d,n){
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
    return(((l$a*l$d)-(l$b*l$c))/sqrt(abs(l$n*(l$a+l$b)*(l$a+l$c))))
  }
  
  dennis.e <- function(a,b,c,d,n){
    return(((a*d)-(b*c))/sqrt(abs(n*(a+b)*(a+c))))
  }
  
  ##############################################################################
  # Digby ----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  #' @param n
  digby.e.1 <- function(l){
    d1 = (l$a*l$d)^(3/4) - (l$b*l$c)^(3/4)
    d2 = (l$a*l$d)^(3/4) - (l$b*l$c)^(3/4)
    d3 = d1/d2
    return(d3)
  }
  
  dibby.e <- function(a,b,c,d,n){
    d1 = (a*d)^(3/4) - (b*c)^(3/4)
    d2 = (a*d)^(3/4) - (b*c)^(3/4)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Dice ----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  dice.e.11 <- function(l){
    return((2*l$a)/((2*l$a)+l$b+l$c))
  }
  
  dice.e.1 <- function(a,b,c,d,n){
    return((2*a)/((2*a)+b+c))
  }
  
  
  ##############################################################################
  # Dice ----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  dice.e.22 <- function(l){
    return(l$a/l$a+l$b)
  }
  
  dice.e.2 <- function(a,b,c,d,n){
    return(a/a+b)
  }
  
  ##############################################################################
  # Dice ----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  dice.e.33 <- function(l){
    return(l$a/l$a+l$c)
  }
  
  dice.e.3 <- function(a,b,c,d,n){
    return(a/a+c)
  }
  
  
  ##############################################################################
  # Dice ----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  doolittle.e.1 <- function(l){
    d1 = ((l$a*l$d) - (l$b*l$c))^2
    d2 = (l$a+l$b)*(l$a+l$c)*(l$c+l$d)*(l$b+l$d)
    d3 = d1/d2
    return(d3)
  }
  
  doolittle.e <- function(a,b,c,d,n){
    d1 = ((a*d) - (b*c))^2
    d2 = (a+b)*(a+c)*(c+d)*(b+d)
    d3 = d1/d2
    return(d3)
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
  
  disperson.e <- function(a,b,c,d,n){
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
  driver.kroeber.e.11 <- function(l){
    return((l$a/2) * ((1/(l$a+l$b)) + (1/(l$a+l$c))))
  }
  
  driver.kroeber.e.1 <- function(a,b,c,d,n){
    return((a/2) * ((1/(a+b)) + (1/(a+c))))
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
  driver.kroeber.e.22 <- function(l){
    return(l$a / sqrt(abs((l$a+l$b)*(l$a+l$c))))
  }
  
  driver.kroeber.e.2 <- function(a,b,c,d,n){
    return(a / sqrt(abs((a+b)*(a+c))))
  }
  
  
  
  ##############################################################################
  # Euclidean ----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param b
  #' @param c
  euclidean.e.1 <- function(l){
    return(l$b+l$c)
  }
  
  euclidean.e <- function(a,b,c,d,n){
    return(b+c)
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
    d1 = (l$n^2) * ((l$n*l$a) - (l$a+l$b) * (l$a+l$c))
    d2 = (l$a+l$b) * (l$a+l$c) * (l$b+l$d) * (l$c+l$d)
    d3 = d1/d2
    return(d3)
  }
  
  eyraud.e <- function(a,b,c,d,n){
    d1 = (n^2) * ((n*a) - (a+b) * (a+c))
    d2 = (a+b) * (a+c) * (b+d) * (c+d)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Fager Mcgowan ----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  fager.mcgowan.e.11 <- function(l){
    d1 = l$a/sqrt(abs((l$a+l$b)*(l$a+l$c)))
    d2 = max((l$a+l$b),(l$a+l$c))/2
    d3 = d1-d2
    return(d3)
  }
  
  fager.mcgowan.e.1 <- function(a,b,c,d,n){
    d1 = a/sqrt(abs((a+b)*(a+c)))
    d2 = max((a+b),(a+c))/2
    d3 = d1-d2
    return(d3)
  }
  
  ##############################################################################
  # Fager Mcgowan ----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  fager.mcgowan.e.22 <- function(l){
    d1 = l$a/sqrt(abs((l$a+l$b)*(l$a+l$c)))
    d2 = 1/ (2 * sqrt(abs(max((l$a+l$b),(l$a+l$c)))))
    d3 = d1-d2
    return(d3)
  }
  
  fager.mcgowan.e.2 <- function(a,b,c,d,n){
    d1 = a/sqrt(abs((a+b)*(a+c)))
    d2 = 1/ (2 * sqrt(abs(max((a+b),(a+c)))))
    d3 = d1-d2
    return(d3)
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
  
  faith.e <- function(a,b,c,d,n){
    return((a+(0.5*d))/(a+b+c+d))
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
  fleiss.e.1 <- function(l){
    d1 = ((l$a*l$d) - (l$b*l$c)) * ((l$a+l$b)*(l$b+l$d)) + ((l$a+l$c)*(l$c+l$d))
    d2 = 2 * (l$a+l$b)*(l$a+l$c)*(l$c+l$d)*(l$b+l$d)
    d3 = d1/d2
    return(d3)
  }
  
  fleiss.e <- function(a,b,c,d,n){
    d1 = ((a*d) - (b*c)) * ((a+b)*(b+d)) + ((a+c)*(c+d))
    d2 = 2 * (a+b)*(a+c)*(c+d)*(b+d)
    d3 = d1/d2
    return(d3)
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
  forbes.e.22 <- function(l){
    d1 = (l$n*l$a) - (l$a+l$b) * (l$a+l$c)
    d2 = l$n * ( min((l$a+l$b),(l$a+l$c)) - ((l$a+l$b) * (l$a+l$c)) )
    d3 = d1/d2
    return(d3)
  }
  
  forbes.e.2 <- function(a,b,c,d,n){
    d1 = (n*a) - (a+b) * (a+c)
    d2 = n * ( min((a+b),(a+c)) - ((a+b) * (a+c)) )
    d3 = d1/d2
    return(d3)
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
  forbes.e.11 <- function(l){
    return((l$n*l$a)/ (l$a+l$b) * (l$a+l$c))
  }
  
  forbes.e.1 <- function(a,b,c,d,n){
    return((n*a)/ (a+b) * (a+c))
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
  
  forbesi.e <- function(a,b,c,d,n){
    return((n*a)/((a+b)*(a+c)))
  }
  
  ##############################################################################
  # Fossum ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param n
  fossum.e.1 <- function(l){
    return((l$n*(l$a-0.5)^2)/(l$a+l$b)*(l$a+l$c))
  }
  
  fossum.e <- function(a,b,c,d,n){
    return((n*(a-0.5)^2)/(a+b)*(a+c))
  }
  
  ##############################################################################
  # Gilbert Well ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param n
  gilbert.well.e.1 <- function(l){
    return(log(l$a)-log(l$n)-log((l$a+l$b)/l$n)-log((l$a+l$c)/l$n))
  }
  
  gilbert.well.e <- function(a,b,c,d,n){
    return(log(a) - log(n) - log((a+b)/n) - log((a+c)/n))
  }
  
  ##############################################################################
  # Goodman Kruskal ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  #' @param n
  goodman.kruskal.e.11 <- function(l){
    z = max(l$a,l$b)+max(l$c,l$d)+max(l$a,l$c)+max(l$b,l$c)
    w = max((l$a+l$c),(l$b+l$d))+max((l$a+l$b),(l$c+l$d))
    p = (z-w)/((2*l$n)-w)
    return(p)
  }
  
  goodman.kruskal.e.1 <- function(a,b,c,d,n){
    z = max(a,b)+max(c,d)+max(a,c)+max(b,c)
    w = max((a+c),(b+d))+max((a+b),(c+d))
    p = (z-w)/((2*n)-w)
    return(p)
  }
  
  ##############################################################################
  # Goodman Kruskal ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  #' @param n
  goodman.kruskal.e.22 <- function(l){
    d1 = 2* min(l$a,l$d) - l$b - l$c
    d2= 2* min(l$a,l$d) + l$b + l$c
    d3 = d1/d2
    return(d3)
  }
  
  goodman.kruskal.e.2 <- function(a,b,c,d,n){
    d1 = 2* min(a,d) - b - c
    d2= 2* min(a,d) + b + c
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Gower ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  gower.e.1 <- function(l){
    d1 = l$a+l$d
    d2 = sqrt(abs((l$a+l$b) * (l$a+l$c) * (l$b+l$d) * (l$c+l$d))) 
    d3 = d1/d2
    return(d3)
  }
  
  gower.e <- function(a,b,c,d,n){
    d1 = a+d
    d2 = sqrt(abs((a+b) * (a+c) * (b+d) * (c+d))) 
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Gower Legendre ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  gower.legendre.e.1 <- function(l){
    return((l$a+l$d)/(l$a+(0.5*(l$b+l$c))+l$d))
  }
  
  gower.legendre.e <- function(a,b,c,d,n){
    return((a+d)/(a+(0.5*(b+c))+d))
  }
  
  ##############################################################################
  # Hamann ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  hamann.e.11 <- function(l){
    return(((l$a+l$d)-(l$b+l$c))/(l$a+l$b+l$c+l$d))
  }
  
  hamann.e.1 <- function(a,b,c,d,n){
    return(((a+d)-(b+c))/(a+b+c+d))
  }
  
  ##############################################################################
  # Harris and Lahey ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  harris.e.1 <- function(l){
    d1 = (l$a * (l$c+l$d) + (l$b+l$d) ) / (2 * (l$a+l$b+l$c))
    d2 = (l$d * (l$a+l$b) + (l$a+l$c)) / (2 * (l$b+l$c+l$d))
    d3 = d1+d2
    return(d3)
  }
  
  harris.e <- function(a,b,c,d,n){
    d1 = (a * (c+d) + (b+d) ) / (2 * (a+b+c))
    d2 = (d * (a+b) + (a+c)) / (2 * (b+c+d))
    d3 = d1+d2
    return(d3)
  }
  
  ##############################################################################
  # Hamann ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  hamann.e.22 <- function(l){
    return((l$a-l$b-l$c+l$d)/(l$a+l$b+l$c+l$d))
  }
  
  hamann.e.2 <- function(a,b,c,d,n){
    return((a-b-c+d)/(a+b+c+d))
  }
  
  ##############################################################################
  # Hamming  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param b
  #' @param c
  hamming.e.1 <- function(l){
    return(l$b+l$c)
  }
  
  hamming.e <- function(a,b,c,d,n){
    return(b+c)
  }
  
  ##############################################################################
  # Hamming  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param b
  #' @param c
  hawkins.dotson.e.1 <- function(l){
    return((1/2) * ((l$a/l$a+l$b+l$c) + (l$d/l$b+l$c+l$d)))
  }
  
  hawkins.dotson.e <- function(a,b,c,d,n){
    return((1/2) * ((a/a+b+c) + (d/b+c+d)))
  }
  
  ##############################################################################
  # Helllinger  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  hellinger.e.1 <- function(l){
    d1 = l$a / sqrt(abs((l$a+l$b) * (l$a+l$c)))
    d2 = 1 - d1
    d3 = sqrt(abs(d2))
    d4 = 2*d3
    return(d4)
  }
  
  hellinger.e <- function(a,b,c,d,n){
    d1 = a / sqrt(abs((a+b) * (a+c)))
    d2 = 1 - d1
    d3 = sqrt(abs(d2))
    d4 = 2*d3
    return(d4)
  }
  
  ##############################################################################
  # Inner Product  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param d
  inner.product.e.1 <- function(l){
    return(l$a+l$d)
  }
  
  inner.product.e <- function(a,b,c,d,n){
    return(a+d)
  }
  
  ##############################################################################
  # Intersection  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  intersection.e.1 <- function(l){
    return(l$a)
  }
  
  intersection.e <- function(a,b,c,d,n){
    return(a)
  }
  
  ##############################################################################
  # Jaccard ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  jaccard.e.1 <- function(l){
    return(l$a/(l$a+l$b+l$c))
  }
  
  jaccard.e <- function(a,b,c,d,n){
    return(a/(a+b+c))
  }
  
  ##############################################################################
  # Jonhson ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  jonhson.e.1 <- function(l){
    return((l$a/(l$a+l$b)) + (l$a/(l$a+l$c)))
  }
  
  jonhson.e <- function(a,b,c,d,n){
    return((a/(a+b)) + (a/(a+c)))
  }
  
  ##############################################################################
  # kulczynski 1 ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  kent.foster.e.11 <- function(l){
    d1 = -(l$b*l$c)
    d2 = l$b*(l$a+l$b) + l$c*(l$a+l$c) + (l$b*l$c)
    d3 = d1/d2
    return(d3)
  }
  
  kent.foster.e.1 <- function(a,b,c,d,n){
    d1 = -(b*c)
    d2 = b*(a+b) + c*(a+c) + (b*c)
    d3 = d1/d2
    return(d3) 
  }
  
    ##############################################################################
  # kulczynski 1 ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  kent.foster.e.22 <- function(l){
    d1 = -(l$b*l$c)
    d2 = l$b*(l$c+l$d) + l$c*(l$b+l$d) + (l$b*l$c)
    d3 = d1/d2
    return(d3) 
  }
  
  kent.foster.e.2 <- function(a,b,c,d,n){
    d1 = -(b*c)
    d2 = b*(c+d) + c*(b+d) + (b*c)
    d3 = d1/d2
    return(d3) 
  }
  
  ##############################################################################
  # kulczynski 1 ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  kulczynski.e.11 <- function(l){
    return(l$a/(l$b+l$c))
  }
  
  kulczynski.e.1 <- function(a,b,c,d,n){
    return(a/(b+c)) 
  }
  
  ##############################################################################
  # kulczynski 2 ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  kulczynski.e.22 <- function(l){
    d1 = (l$a/2) * ((2*l$a)+l$b+l$c)
    d2 = (l$a+l$b) * (l$a+l$c)
    d3 = d1/d2
    return(d3)
  }
  
  kulczynski.e.2 <- function(a,b,c,d,n){
    d1 = (a/2) * ((2*a)+b+c)
    d2 = (a+b) * (a+c)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # kulczynski 2 ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  kulczynski.e.33 <- function(l){
    return((1/2) * ((l$a/l$a+l$b) + (l$a/l$a+l$c)))
  }
  
  kulczynski.e.3 <- function(a,b,c,d,n){
    return((1/2) * ((a/a+b) + (a/a+c)))
  }
  
  ##############################################################################
  # kulczynski 2 ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  kulczynski.e.44 <- function(l){
    return(l$a/l$b+l$c)
  }
  
  kulczynski.e.4 <- function(a,b,c,d,n){
    return(a/b+c)
  }
  
  ##############################################################################
  # kulczynski 2 ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  kuder.richardson.e.1 <- function(l){
    d1 = 4 * ((l$a*l$d) - (l$b*l$c))
    d2 = ((l$a+l$b)*(l$c+l$d)) + ((l$a+l$c)*(l$b+l$d)) + (2 *((l$a*l$d)-(l$b*l$c)))
    d3 = d1/d2
    return(d3)
  }
  
  kuder.richardson.e <- function(a,b,c,d,n){
    d1 = 4 * ((a*d) - (b*c))
    d2 = ((a+b)*(c+d)) + ((a+c)*(b+d)) + (2 *((a*d)-(b*c)))
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Lance Williams ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  lance.williams.e.1 <- function(l){
    return((l$b+l$c)/((2*l$a)+l$b+l$c))
  }
  
  lance.williams.e <- function(a,b,c,d,n){
    return((b+c)/((2*a)+b+c))
  }
  
  ##############################################################################
  # Lance Williams ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  loevinger.e.1 <- function(l){
    d1 = (l$a*l$d) - (l$b*l$c)
    d2 = min(((l$a+l$b)*(l$b+l$d)),((l$a+l$c)*(l$c+l$d)))
    d3 = d1/d2
    return(d3)
  }
  
  loevinger.e <- function(a,b,c,d,n){
    d1 = (a*d) - (b*c)
    d2 = min(((a+b)*(b+d)),((a+c)*(c+d)))
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Manhattan ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param b
  #' @param c
  maxwell.pilliner.e.1 <- function(l){
    d1 = 2 * ((l$a*l$d) - (l$b*l$c))
    d2 = (l$a+l$b)*(l$c+l$d) + (l$c+l$d)*(l$b+l$d)
    d3 = d1/d2
    return(d3)
  }
  
  maxwell.pilliner.e <- function(a,b,c,d,n){
    d1 = 2 * ((a*d) - (b*c))
    d2 = (a+b)*(c+d) + (c+d)*(b+d)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Manhattan ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param b
  #' @param c
  manhattan.e.1 <- function(l){
    return(l$b+l$c)
  }
  
  manhattan.e <- function(a,b,c,d,n){
    return(b+c)
  }
  
  ##############################################################################
  # Mcconnaughey ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  mcconnaughey.e.1 <- function(l){
    return(((l$a^2)-(l$b-l$c))/(l$a+l$b)*(l$a+l$c))
  }
  
  mcconnaughey.e <- function(a,b,c,d,n){
    return(((a^2) - (b-c))/(a+b)*(a+c))
  }
  
  ##############################################################################
  # Mean Manhattan ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  mean.manhattan.e.1 <- function(l){
    return((l$b+l$c)/(l$a+l$b+l$c+l$d))
  }
  
  mean.manhattan.e <- function(a,b,c,d,n){
    return((b+c)/(a+b+c+d))
  }
  
  ##############################################################################
  # Michael ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  michael.e.1 <- function(l){
    return(4*((l$a*l$d)-(l$b*l$c))/((l$a+l$d)^2)+((l$b+l$c)^2))
  }
  
  michael.e <- function(a,b,c,d,n){
    return(4*((a*d)-(b*c))/((a+d)^2) + ((b+c)^2))
  }
  
  ##############################################################################
  # Minowski ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param b
  #' @param c
  minowski.e.1 <- function(l){
    return((l$b+l$c)^(1/1))
  }
  
  minowski.e <- function(a,b,c,d,n){
    return((b+c)^(1/1))
  }
  
  ##############################################################################
  # Mountford ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  mountford.e.11 <- function(l){
    return(l$a/(0.5*(((l$a*l$b)+(l$a*l$c))+(l$b*l$c))))
  }
  
  mountford.e.1 <- function(a,b,c,d,n){
    return(a/(0.5*(((a*b)+(a*c))+(b*c))))
  }
  
  ##############################################################################
  # Mountford ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  mountford.e.22 <- function(l){
    return((2*l$a) / (a*(l$a+l$c) + 2*l$b*l$c))
  }
  
  mountford.e.2 <- function(a,b,c,d,n){
    return((2*a) / (a*(a+c) + 2*b*c))
  }
  
  ##############################################################################
  # Nei Li  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  nei.li.e.1 <- function(l){
    return((2*l$a)/((l$a+l$b)+(l$a+l$c)))
  }
  
  nei.li.e <- function(a,b,c,d,n){
    return((2*a)/((a+b)+(a+c)))
  }
  
  ##############################################################################
  # Ochiai 1  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  ochiai.e.11 <- function(l){
    return(l$a/sqrt(abs((l$a+l$b)*(l$a+l$c))))
  }
  
  ochiai.e.1 <- function(a,b,c,d,n){
    return(a/sqrt(abs((a+b)*(a+c))))
  }
  
  ##############################################################################
  # Ochiai 2  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  ochiai.e.22 <- function(l){
    return((l$a*l$d)/sqrt(abs((l$a+l$d)*(l$a+l$c)*(l$d+l$b)*(l$d+l$c))))
  }
  
  ochiai.e.2 <- function(a,b,c,d,n){
    return((a*d)/sqrt(abs((a+d)*(a+c)*(d+b)*(d+c))))
  }
  
  ##############################################################################
  # Ochiai 2  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  ochiai.e.33 <- function(l){
    return((l$a*l$d)/sqrt(abs((l$a+l$b)*(l$a+l$c)*(l$b+l$d)*(l$c+l$d))))
  }
  
  ochiai.e.3 <- function(a,b,c,d,n){
    return((a*d)/sqrt(abs((a+b)*(a+c)*(b+d)*(c+d))))
  }
  
  ##############################################################################
  # Otsuka  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  otsuka.e.1 <- function(l){
    return(l$a/((l$a+l$b)*(l$a+l$c))^0.5)
  }
  
  otsuka.e <- function(a,b,c,d,n){
    return(a/((a+b)*(a+c))^0.5)
  }
  
  ##############################################################################
  # Pattern Difference  ---------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  pattern.difference.e.1 <- function(l){
    return((4*l$b*l$c)/((l$a+l$b+l$c+l$d)^2))
  }
  
  pattern.difference.e <- function(a,b,c,d,n){
    return((4*b*c)/((a+b+c+d)^2))
  }
  
  ##############################################################################
  # Pearson 1 ------------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  pearson.e.11 <- function(l){
    d1 = l$n * (((l$a*l$d) - (l$b*l$c))^2)
    d2 = (l$a+l$b) * (l$a+l$c) * (l$c+l$d) * (l$b+l$d)
    d3 = d1/d2
    return(d3)
  }
  
  pearson.e.1 <- function(a,b,c,d,n){
    d1 = n * (((a*d) - (b*c))^2)
    d2 = (a+b) * (a+c) * (c+d) * (b+d)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Pearson 2 ------------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  #' @param n
  pearson.e.22 <- function(l){
    d1 = l$n * (((l$a*l$d) - (l$b*l$c))^2)
    d2 = (l$a+l$b) * (l$a+l$c) * (l$c+l$d) * (l$b+l$d)
    d3 = d1/d2
    d4 = (d3/n+d3)^(1/2)
    return(d4)
  }
  
  pearson.e.2 <- function(a,b,c,d,n){
    d1 = n * (((a*d) - (b*c))^2)
    d2 = (a+b) * (a+c) * (c+d) * (b+d)
    d3 = d1/d2
    d4 = (d3/n+d3)^(1/2)
    return(d4)
  }
  
  ##############################################################################
  # Pearson 3 ------------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  pearson.e.33 <- function(l){
    d1 = (l$a*l$d) - (l$b*l$c)
    d2 = sqrt( abs( (l$a+l$b)*(l$a+l$c)*(l$b+l$c)*(l$c+l$d) ) )
    d3 = d1/d2
    d4 = (d3/b+d3)^(1/2)
    return(d4)
  }
  
  pearson.e.3 <- function(a,b,c,d,n){
    d1 = (a*d) - (b*c)
    d2 = sqrt( abs( (a+b)*(a+c)*(b+c)*(c+d) ) )
    d3 = d1/d2
    d4 = (d3/b+d3)^(1/2)
    return(d4)
  }
  
  ##############################################################################
  # Pearson Heron 1------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  pearson.heron.e.11 <- function(l){
    d1 = (l$a*l$d) - (l$b*l$c)
    d2 = (l$a+l$b) * (l$a+l$c) * (l$b+l$d) * (l$c+l$d)
    d3 = sqrt(abs(d2))
    d4 = d1/d2
    return(d4)
  }
  
  pearson.heron.e.1 <- function(a,b,c,d,n){
    d1 = (a*d) - (b*c)
    d2 = (a+b) * (a+c) * (b+d) * (c+d)
    d3 = sqrt(abs(d2))
    d4 = d1/d2
    return(d4)
  }
  
  
  ##############################################################################
  # Pearson Heron 2------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  pearson.heron.e.22 <- function(l){
    d1 = pi * sqrt(abs(l$b*l$c))
    d2  = sqrt(abs(l$a*l$d)) + sqrt(abs(l$b*l$c))
    d3 = d1/d2
    d4 = cos(d3)
    return(d4)
  }
  
  pearson.heron.e.2 <- function(a,b,c,d,n){
    d1 = pi * sqrt(abs(b*c))
    d2  = sqrt(abs(a*d)) + sqrt(abs(b*c))
    d3 = d1/d2
    d4 = cos(d3)
    return(d4)
  }
  
  ##############################################################################
  # Pierce --------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  peirce.e.11 <- function(l){
    d1 = (l$a*l$b) + (l$b*l$c)
    d2 = (l$a*l$b) + (2*l$b*l$c) + (l$c*l$d)
    d3 = d1/d2
    return(d3)
  }
  
  peirce.e.1 <- function(a,b,c,d,n){
    d1 = (a*b) + (b*c)
    d2 = (a*b) + (2*b*c) + (c*d)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Pierce --------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  peirce.e.22 <- function(l){
    d1 = (l$a*l$d) - (l$b*l$c)
    d2 = (l$a+l$b) * (l$c+l$d)
    d3 = d1/d2
    return(d3)
  }
  
  peirce.e.2 <- function(a,b,c,d,n){
    d1 = (a*d) - (b*c)
    d2 = (a+b) * (c+d)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Pierce --------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  peirce.e.33 <- function(l){
    d1 = (l$a*l$d) - (l$b*l$c)
    d2 = (l$a+l$c) * (l$b+l$d)
    d3 = d1/d2
    return(d3)
  }
  
  peirce.e.3 <- function(a,b,c,d,n){
    d1 = (a*d) - (b*c)
    d2 = (a+c) * (b+d)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Roger Tanimoto -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  roger.tanimoto.e.1 <- function(l){
    d1 = l$a + l$d
    d2 = l$a + (2*(l$b+l$c)) + l$d
    d3 = d1/d2
    return(d3)
  }
  
  roger.tanimoto.e <- function(a,b,c,d,n){
    d1 = a + d
    d2 = a + (2*(b+c)) + d
    d3 = d1/d2
    return(d3)
  }
  
  
  ##############################################################################
  # Roger Tanimoto -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  rogot.goldberd.e.1 <- function(l){
    d1 = l$a/(l$a+l$b)+(l$a+l$c)
    d2 = l$d/(l$c+l$d)+(l$b+l$d)
    d3 = d1+d2
    return(d3)
  }
  
  rogot.goldberd.e <- function(a,b,c,d,n){
    d1 = a/(a+b)+(a+c)
    d2 = d/(c+d)+(b+d)
    d3 = d1+d2
    return(d3)
  }
  
  ##############################################################################
  # Russel Rao -----------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  russel.rao.e.1 <- function(l){
    return(l$a/(l$a+l$b+l$c+l$d))
  }
  
  russel.rao.e <- function(a,b,c,d,n){
    return(a/(a+b+c+d))
  }
  
  ##############################################################################
  # Shape Difference -----------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  #' @param n
  shape.difference.e.1 <- function(l){
    d1 = l$n * (l$b+l$c) - ((l$b-l$c)^2)
    d2 = (l$a + l$b + l$c + l$d)^2
    d3 = d1/d2
    return(d3)
  }
  
  shape.difference.e <- function(a,b,c,d,n){
    d1 = n * (b+c) - ((b-c)^2)
    d2 = (a + b + c + d)^2
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Simpson --------------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  simpson.e.1 <- function(l){
    return(l$a/min((l$a+l$b),(l$a+l$c)))
  }
  
  simpson.e <- function(a,b,c,d,n){
    return(a/min((a+b),(a+c)))
  }
  
  ##############################################################################
  # Size Difference ------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  size.difference.e.1 <- function(l){
    d1 = (l$b+l$c)^2
    d2 = (l$a+l$b+l$c+l$d)^2
    d3 = d1/d2
    return(d3)
  }
  
  size.difference.e <- function(a,b,c,d,n){
    d1 = (b+c)^2
    d2 = (a+b+c+d)^2
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Sokal Michener -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  scott.e.1 <- function(l){
    d1 = 4*l$a*l$d - ((l$b+l$c)^2)
    d2 = ((l$a+l$b) + (l$a+l$c)) * ((l$c+l$d)+(l$b+l$d))
    d3 = d1/d2
    return(d3)
  }
  
  scott.e <- function(a,b,c,d,n){
    d1 = 4*a*d - ((b+c)^2)
    d2 = ((a+b) + (a+c)) * ((c+d)+(b+d))
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Sokal Michener -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  sokal.michener.e.1 <- function(l){
    return((l$a+l$d)/(l$a+l$b+l$c+l$d))
  }
  
  sokal.michener.e <- function(a,b,c,d,n){
    return((a+d)/(a+b+c+d))
  }
  
  ##############################################################################
  # Sokal Michener -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  sokal.sneath.e.11 <- function(l){
    return(l$a/(l$a+(2*l$b)+(2*l$c)))
  }
  
  sokal.sneath.e.1 <- function(a,b,c,d,n){
    return(a/(a+(2*b)+(2*c)))
  }
  
  ##############################################################################
  # Sokal Michener -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  sokal.sneath.e.21 <- function(l){
    return(2*(l$a+l$d)/((2*l$a)+l$b+l$c+(2*l$d)))
  }
  
  sokal.sneath.e.2 <- function(a,b,c,d,n){
    return(2*(a+d)/((2*a)+b+c+(2*d)))
  }
  
  ##############################################################################
  # Sokal Michener -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  sokal.sneath.e.31 <- function(l){
    return((l$a+l$d)/(l$b+l$c))
  }
  
  sokal.sneath.e.3 <- function(a,b,c,d,n){
    return((a+d)/(b+c))
  }
  
  ##############################################################################
  # Sokal Michener -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  sokal.sneath.e.4b1 <- function(l){
    return((l$a/l$a+l$b) * (l$a/l$a+l$c) * (l$d/l$b+l$c) * (l$d/l$b+l$d))
  }
  
  sokal.sneath.e.4b <- function(a,b,c,d,n){
    return((a/a+b) * (a/a+c) * (d/b+c) * (d/b+d))
  }
  
  ##############################################################################
  # Sokal Michener -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  sokal.sneath.e.4a1 <- function(l){
    d1 = 1/4
    d2 = (l$a/l$a+l$b) * (l$a/l$a+l$c) * (l$a/l$c+l$d) * (l$a/l$b+l$d)
    d3 = d1/d2
    return(d3)
  }
  
  sokal.sneath.e.4a <- function(a,b,c,d,n){
    d1 = 1/4
    d2 = (a/a+b) * (a/a+c) * (a/c+d) * (a/b+d)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Sokal Michener -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  sokal.sneath.e.5b1 <- function(l){
    d1 = l$a*l$d
    d2 = sqrt(abs((l$a+l$d) * (l$a+l$c) * (l$b+l$d) * (l$c+l$d)))
    d3 = d1/d2
    return(d3)
  }
  
  sokal.sneath.e.5b <- function(a,b,c,d,n){
    d1 = a*d
    d2 = sqrt(abs((a+d) * (a+c) * (b+d) * (c+d)))
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Sokal Michener -------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  sokal.sneath.e.5a1 <- function(l){
    d1 = l$a*l$d
    d2 = ((l$a+l$d) * (l$a+l$c) * (l$b+l$d) * (l$c+l$d))^0.5
    d3 = d1/d2
    return(d3)
  }
  
  sokal.sneath.e.5a <- function(a,b,c,d,n){
    d1 = a*d
    d2 = ((a+d) * (a+c) * (b+d) * (c+d))^0.5
    d3 = d1/d2
    return(d3)
  }
  
  ###############################################################################
  # Sorgenfrei
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  sorgenfrei.e.1 <- function(l){
    d1 = l$a^2
    d2 = (l$a+l$b) * (l$a+l$b)
    d3 = d1/d2
    return(d3)
  }
  
  sorgenfrei.e <- function(a,b,c,d,n){
    d1 = a^2
    d2 = (a+b) * (a+b)
    d3 = d1/d2
    return(d3)
  }
  
  ###############################################################################
  # Square Euclidean
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  square.euclidean.e.1 <- function(l){
    return(sqrt(abs((l$b+l$c)^2)))
  }
  
  square.euclidean.e <- function(a,b,c,d,n){
    return(sqrt(abs((b+c)^2)))
  }
  
  ###############################################################################
  # Stiles
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  stiles.e.1 <- function(l){
    d1 = l$n * ( abs((l$a*l$d) - (l$n/2)) )^2
    d2 = (l$a+l$b) * (l$a+l$c) * (l$b+l$d) * (l$c+l$d)
    d3 = d1/d2
    d4 = log10(d3)
    return(d4)
  }
  
  stiles.e <- function(a,b,c,d,n){
    d1 = n * ( abs((a*d) - (n/2)) )^2
    d2 = (a+b) * (a+c) * (b+d) * (c+d)
    d3 = d1/d2
    d4 = log10(d3)
    return(d4)
  }
  
  ###############################################################################
  # Tanimoto
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  tanimoto.e.1 <- function(l){
    d1 = l$a
    d2 = (l$a+l$b) + (l$a+l$c) - l$a
    d3 = d1/d2
    return(d3)
  }
  
  tanimoto.e <- function(a,b,c,d,n){
    d1 = a
    d2 = (a+b) + (a+c) - a
    d3 = d1/d2
    return(d3)
  }
  
  ###############################################################################
  # Tarantula
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  tarantula.e.1 <- function(l){
    d1 = l$l$a* (l$c+l$d)
    d2 = l$c * (l$a+l$d)
    d3 = d1/d2
    return(d3)
  }
  
  tarantula.e <- function(a,b,c,d,n){
    d1 = a* (c+d)
    d2 = c* (a+d)
    d3 = d1/d2
    return(d3)
  }
  
  ###############################################################################
  # Tarwid
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  tarwid.e.1 <- function(l){
    d1 = l$n*l$a - (l$a+l$b) * (l$a+l$c)
    d2 = l$n*l$a + (l$a+l$b) * (l$a+l$c)
    d3 = d1/d2
    return(d3)
  }
  
  tarwid.e <- function(a,b,c,d,n){
    d1 = n*a - (a+b) * (a+c)
    d2 = n*a + (a+b) * (a+c)
    d3 = d1/d2
    return(d3)
  }
  
  ###############################################################################
  # 3W Jaccard
  #' Similaritie Measure for categorial data
  #'
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  three.w.jaccard.e.1 <- function(l){
    d1 = 3 * l$a
    d2 = 3 * l$a + b + c
    d3 = d1/d2
    return(d3)
  }
  
  three.w.jaccard.e <- function(a,b,c,d,n){
    d1 = 3 * l$a
    d2 = 3 * l$a + b + c
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Vari -----------------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  vari.e.1 <- function(l){
    d1 = l$b+l$c
    d2 = 4 * (l$a + l$b + l$c + l$d)
    d3 = d1/d2
    return(d3)
  }
  
  vari.e <- function(a,b,c,d,n){
    d1 = b+c
    d2 = 4 * (a + b + c + d)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Yuleq 1 --------------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  yule.e.11 <- function(l){
    d1 = (l$a*l$d) - (l$b*l$c)
    d2 = (l$a*l$d) - (l$b*l$c)
    d3 = d1/d2
    return(d3)
  }
  
  yule.e.1 <- function(a,b,c,d,n){
    d1 = (a*d) - (b*c)
    d2 = (a*d) - (b*c)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Yuleq 2 --------------------------------------------------------------------
  #' Similaritie Measure for categorial data
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  yule.e.21 <- function(l){
    d1 = 2*l$b*l$c
    d2 = (l$a*l$d) + (l$b*l$c)
    d3 = d1/d2
    return(d3)
  }
  
  yule.e.2 <- function(a,b,c,d,n){
    d1 = 2*b*c
    d2 = (a*d) + (b*c)
    d3 = d1/d2
    return(d3)
  }
  
  ##############################################################################
  # Yule ---------------------------------------------------------------------
  #' Similaritie Measure for categorial data 
  #' @family Multi_Label_Binary_Measures
  #' @param a
  #' @param b
  #' @param c
  #' @param d
  yule.e.31 <- function(l){
    d1 = sqrt(abs(l$a*l$d)) - sqrt(abs(l$b*l$c))
    d2 = sqrt(abs(l$a*l$d)) - sqrt(abs(l$b*l$c))
    d3 = d1/d2
    return(d3)
  }
  
  yule.e.3 <- function(a,b,c,d,n){
    d1 = sqrt(abs(a*d)) - sqrt(abs(b*c))
    d2 = sqrt(abs(a*d)) - sqrt(abs(b*c))
    d3 = d1/d2
    return(d3)
  }
  
  
  ###############################################################################
  # get name list
  getNamesListFunctions <- function(){
    retorno = list()
    names_function = c("ample.e",
                       "anderberg.e",
                       "baroni.urbani.buser.e.1",
                       "baroni.urbani.buser.e.2",
                       "braun.blanquet.e",
                       "bray.curtis.e",
                       "canberra.e",
                       "chord.e",
                       "cityblock.e",
                       "clement.e",
                       "cohen.e",
                       "cole.e.1",
                       "cole.e.2",
                       "cole.e.3",
                       "cosine.e", 
                       "czekanowski.e",
                       "dennis.e", 
                       "dibby.e",
                       "dice.e.1",
                       "dice.e.2",
                       "dice.e.3",
                       "disperson.e",
                       "doolittle.e",
                       "driver.kroeber.e.1",
                       "driver.kroeber.e.2",
                       "euclidean.e",
                       "eyraud.e",
                       "fager.mcgowan.e.1",
                       "fager.mcgowan.e.2",
                       "faith.e",
                       "fleiss.e",
                       "forbes.e.1",
                       "forbes.e.2",
                       "forbesi.e",
                       "fossum.e",
                       "gilbert.well.e",
                       "goodman.kruskal.e.1",
                       "goodman.kruskal.e.2",
                       "gower.e",
                       "gower.legendre.e", 
                       "hamann.e.1",
                       "hamann.e.2",
                       "hamming.e",
                       "harris.e",
                       "hawkins.dotson.e",
                       "hellinger.e", 
                       "inner.product.e", 
                       "intersection.e",
                       "jaccard.e",
                       "jonhson.e",
                       "kent.foster.e.1",
                       "kent.foster.e.2",
                       "kuder.richardson.e",
                       "kulczynski.e.1",
                       "kulczynski.e.2",
                       "kulczynski.e.3",
                       "kulczynski.e.4",
                       "lance.williams.e",
                       "loevinger.e",
                       "manhattan.e",
                       "maxwell.pilliner.e",
                       "mcconnaughey.e",
                       "mean.manhattan.e",
                       "michael.e",
                       "minowski.e",
                       "mountford.e.1",
                       "mountford.e.2",
                       "nei.li.e",
                       "ochiai.e.1",
                       "ochiai.e.2",
                       "ochiai.e.3",
                       "otsuka.e",
                       "pattern.difference.e",
                       "pearson.e.1",
                       "pearson.e.2",
                       "pearson.e.3",
                       "pearson.heron.e.1",
                       "pearson.heron.e.2",
                       "peirce.e.1",
                       "peirce.e.2",
                       "peirce.e.3",
                       "roger.tanimoto.e",
                       "rogot.goldberd.e",
                       "russel.rao.e",
                       "scott.e",
                       "shape.difference.e",
                       "simpson.e",
                       "size.difference.e",
                       "sokal.michener.e",
                       "sokal.sneath.e.1",
                       "sokal.sneath.e.2",
                       "sokal.sneath.e.3",
                       "sokal.sneath.e.4a",
                       "sokal.sneath.e.4b",
                       "sokal.sneath.e.5a",
                       "sokal.sneath.e.5b",
                       "sorgenfrei.e",
                       "square.euclidean.e",
                       "stiles.e",
                       "tanimoto.e",
                       "tarantula.e",
                       "tarwid.e",
                       "three.w.jaccard.e",
                       "vari.e",
                       "yule.e.1",
                       "yule.e.2",
                       "yule.e.3")
    retorno$nomes_funcoes = names_function
    retorno$tamanho = length(names_function)
    return(retorno)
    gc()
  }
  
  
  ################################################################################
  # any errors, please, contact me: elainececiliagatto@gmail.com                 #
  ################################################################################
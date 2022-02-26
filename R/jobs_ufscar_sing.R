##################################################################################################

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
FolderRoot = "~/Documents/DOUTORADO/[2022]CODIGOS/MultiLabelSimilaritiesMeasures"
FolderScripts = paste(FolderRoot, "/R/", sep="")

library(stringr)

setwd(FolderRoot)
datasets = data.frame(read.csv("datasets-2022.csv"))
n = nrow(datasets)

FolderJob = paste(FolderRoot, "/Jobs-Singularity", sep="")
if(dir.exists(FolderJob)==FALSE){dir.create(FolderJob)}

  i = 1
  while(i<=n){

    dataset = datasets[i,]
    cat("\n\tDataset:", dataset$Name)

    nome_job = paste("mlsm-", dataset$Name, sep="")

    nome = paste(FolderJob, "/", nome_job, ".sh", sep="")
    output.file <- file(nome, "wb")

    write("#!/bin/bash", file = output.file)

    str1 = paste("#SBATCH -J ", nome_job, sep="")
    write(str1, file = output.file, append = TRUE)

    write("#SBATCH -o %j.out", file = output.file, append = TRUE)

    write("#SBATCH -n 1", file = output.file, append = TRUE)

    write("#SBATCH -c 10", file = output.file, append = TRUE)

    write("#SBATCH --partition slow", file = output.file, append = TRUE)

    write("#SBATCH -t 720:00:00", file = output.file, append = TRUE)

    #write("#SBATCH -t 128:00:00", file = output.file, append = TRUE)

    write("#SBATCH --mail-user=elainegatto@estudante.ufscar.br",
          file = output.file, append = TRUE)

    write("#SBATCH --mail-type=ALL", file = output.file, append = TRUE)
    write("", file = output.file, append = TRUE)

    str2 = paste("local_job=",  "\"/scratch/", nome_job, "\"", sep="")
    write(str2, file = output.file, append = TRUE)

    write("function clean_job(){", file = output.file, append = TRUE)

    str3 = paste(" echo", "\"Limpando ambiente...\"", sep=" ")
    write(str3, file = output.file, append = TRUE)

    str4 = paste(" rm -rf ", "\"${local_job}\"", sep="")
    write(str4, file = output.file, append = TRUE)
    write("}", file = output.file, append = TRUE)

    write("trap clean_job EXIT HUP INT TERM ERR",
          file = output.file, append = TRUE)
    write("", file = output.file, append = TRUE)

    write("set -eE", file = output.file, append = TRUE)
    write("umask 077", file = output.file, append = TRUE)
    str5 = paste("echo RUN ", nome_job, sep="")
    write(str5, file = output.file, append = TRUE)

    str6 = paste("sbcast /home/u704616/experimentos-R.simg /scratch/",
                 nome_job, sep="")
    write(str6, file = output.file, append = TRUE)

    str7 = paste("singularity instance start --bind ~/.config/rclone/:/root/.config/rclone experimentos-R.simg ", nome_job ,sep="")
    write(str7, file = output.file, append = TRUE)

    str8 = paste("singularity run --app Rscript instance://",
                 nome_job, " /MultiLabelSimilaritiesMeasures/R/mlsm.R ",
                 dataset$Id , " 10 10 ",
                 "\"/scratch/", nome_job, "\"", sep="")
    write(str8, file = output.file, append = TRUE)

    str9 = paste("singularity instance stop ", nome_job, sep="")
    write(str9, file = output.file, append = TRUE)
    write("", file = output.file, append = TRUE)

    close(output.file)

    i = i + 1
    gc()
  }
# MultiLabel Similarities Measures

Compute similarities measures (categorical data) for all labels in label space for a multilabel dataset

## Tutorial

https://rpubs.com/cissagatto/MultiLabelSimilaritiesMeasures

## How to cite 
@misc{Gatto2021, author = {Gatto, E. C.}, title = {Compute Similarities Measures for MultiLabel Classification}, year = {2021}, publisher = {GitHub}, journal = {GitHub repository}, howpublished = {\url{https://github.com/cissagatto/MultiLabelSimilaritiesMeasures}}}

# Scripts
This code has the following script in the R folder

1. functions_contingency_table_multilabel.R
2. functions_measures_binary_data.R
3. functions_multilabel_binary_measures.R
4. libraries.R
5. run.R


## Folder Path
Place a copy of this code in _"C:/Users/[username]/BellPartitionsMultiLabel"_ or _"/home/username/BellPartitionsMultiLabel"_. You can change the path in the code if you want.

## File "datasets.csv"
A file called *datasets.csv* must be in the root project folder. This file is used to read information about the datasets and they are used in the code. All 74 datasets available in cometa (https://cometa.ujaen.es/) are in this file. If you want to use another dataset, please, add the following information about the dataset in the file:

*Id, Name, Domain, Labels, Instances, Attributes, Inputs, Labelsets, Single, Max freq, Card, Dens, MeanIR, Scumble, TCS, AttStart, AttEnd, LabelStart, LabelEnd, xn, yn, gridn*

The _"Id"_ of the dataset is a mandatory parameter in the command line to run all code. The fields are used in a lot of internal functions. Please, make sure that this information is available before running the code. *xn* and *yn* correspond to a dimension of the quadrangular map for kohonen, and *gridn* is *xn* * *yn*. Example: xn = 4, yn = 4, gridn = 16.

## Datasets
You will need the originals datasets to run this code. Please, make sure that the folder *datasets* contain the arrf and xml files from cometa.

## RUN
The script RUN is a example of how to test each measure function.




## Acknowledgment
This study is financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES) - Finance Code 001

## Links

[Post-Graduate Program in Computer Science](http://ppgcc.dc.ufscar.br/pt-br)

[Biomal](http://www.biomal.ufscar.br/)

[Computer Department](https://site.dc.ufscar.br/)

[CAPES](https://www.gov.br/capes/pt-br)

[Embarcados](https://www.embarcados.com.br/author/cissa/)

[Linkedin](https://www.linkedin.com/in/elainececiliagatto/)

[Linkedin](https://www.linkedin.com/company/27241216)

[Instagram](https://www.instagram.com/professoracissa/)

[Facebook](https://www.facebook.com/ProfessoraCissa/)

[Twitter](https://twitter.com/professoracissa)

# Thanks



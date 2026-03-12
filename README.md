<b><h3>uem915</h3></b>


<b>Install uem915 R package:</b>
```r
source("https://raw.githubusercontent.com/fdumbioinfo/rtools/main/uem915-demo/0-uem915-install-r-universe.r")
```

<b>Cours introduction à R</b>

```r


#
# UEM915 analyse omique
#
# introduction a R
#
# installation du package UEM915
#
if(!require("moalannotgene",quietly=TRUE)){install.packages("moalannotgene",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalannotensg",quietly=TRUE)){install.packages("moalannotensg",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalannotenst",quietly=TRUE)){install.packages("moalannotenst",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalannotensp",quietly=TRUE)){install.packages("moalannotensp",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbhs",quietly=TRUE)){install.packages("moalstringdbhs",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbmm",quietly=TRUE)){install.packages("moalstringdbmm",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbrn",quietly=TRUE)){install.packages("moalstringdbrn",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbdr",quietly=TRUE)){install.packages("moalstringdbdr",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbss",quietly=TRUE)){install.packages("moalstringdbss",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("uem915",quietly=TRUE)){install.packages("uem915",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
# chargement du package
library(uem915)
uem915::env()
#
# working directory: dossier dans lequel les fichiers sont cr??s
setwd("~/Desktop/IPSIT/enseignements/UEM915/2026/cours/rwd")
#
# les objets de donnees R
#
# vecteurs
# 
10+15
c(10)
letters
c(1,4,6)
5:10
# operateur d'affectation <- ->
a <- 10
10 -> a
# affichage
a
print(a)
# operateurs pipe %>% 
if( !require( "magrittr", quietly=TRUE ) ){ install.packages("magrittr") }
library(magrittr)
# matrice : tableau a 2 dimension
rnorm(1000) %>% matrix(ncol=10) -> m0
m0 %>% head
m0 %>% dim
m0 %>% str
m0 %>% graphics::hist(breaks=100)
rep(c("WT","KO"),c(5,5)) %>% as.factor -> GROUP
m0 %>% uem915:::boxplot(GROUP)
# data.frame
rnorm(1000) %>% matrix(ncol=10) %>% data.frame -> m0
# list : peut contenir des object de type diff?rents
# permet faire des listes de listes pour ranger des donn?es
list(
  list( nom = "gzerg", prenom = "fzef" , age = 10 ),
  list( nom = "dazd", prenom = "erge" , age = 40 ) )
# Fonction de description des objets
a %>% mode # donne le type de l'objet (character, integer, numeric, boolean)
a %>% class 
a %>% length
GROUP %>% class
m0 %>% dim
m0 %>% head
m0 %>% str
# -----
# anova CHAUSSURE
# -----
paste("s",1:20,sep="") -> SampleID
c( 170.3,175.3, 180.3,182.1, 158.8,172.1, 185.3,186.1, 155.3,158.1,
   186.3,188.3, 156.7,160.1, 178.8,183.1, 175.3,176.5, 181.3,182.8) -> TAILLE
rep(c("F","F","H","H"),5) %>% as.factor -> GENRE
rep(c("SANS","AVEC"),10) %>% ordered(c("SANS","AVEC")) -> CHAUSSURE
CHAUSSURE
paste("IND",rep(1:10,rep(2,10)),sep="") -> CASE
CASE %>% ordered(CASE %>% unique) -> CASE
CASE
paste(SampleID,CHAUSSURE,GENRE,CASE,sep="_") -> SampleName
data.frame(SampleID,CHAUSSURE,GENRE,CASE,SampleName,TAILLE) -> m0
m0 %>% head
m0 %>% dim
m0$CHAUSSURE %>% "=="("SANS") %>% which %>% m0$TAILLE[.] %>% mean -> meanS
m0$CHAUSSURE %>% "=="("AVEC") %>% which %>% m0$TAILLE[.] %>% mean -> meanA
meanA/meanS
m0$TAILLE %>% hist(breaks=10)
# CHAUSSURE
data.frame(SANS=m0$CHAUSSURE %>% "=="("SANS") %>% which %>% m0$TAILLE[.],
           AVEC=m0$CHAUSSURE %>% "=="("AVEC") %>% which %>% m0$TAILLE[.]) -> m1
m1 %>% uem915::boxplot(c("SANS","AVEC") %>% as.factor)
m1 %>% uem915::boxplot(m0$CHAUSSURE %>% levels %>% as.factor)
# CASE
m1 <- foreach(i=1:length(m0$CASE %>% unique),.combine="cbind") %do%
  {
    m0$CASE %>% unique %>% "["(i) %>% "=="(m0$CASE) %>% which %>% m0$TAILLE[.] %>% data.frame
  }
m1
m1 %>% data.frame %>% uem915::boxplot(m0$CASE %>% unique %>% ordered(m0$CASE %>% unique))
# GENRE
m1 <- foreach(i=1:length(m0$GENRE %>% table),.combine="cbind") %do%
  {
    m0$GENRE %>% levels %>% "["(i) %>% grep(m0$GENRE) %>% m0$TAILLE[.]
  }
m1 %>% uem915::boxplot(m0$GENRE %>% levels %>% as.factor)
# anova
lm(formula = "TAILLE~CHAUSSURE", data = m0) %>% aov %>% summary
lm(formula = "TAILLE~CHAUSSURE+CASE", data = m0) %>% aov %>% summary
lm(formula = "TAILLE~CHAUSSURE*GENRE+CASE", data = m0) %>% aov %>% summary
#
```


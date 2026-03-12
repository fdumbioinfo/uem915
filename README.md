<b><h3>uem915</h3></b>


<b>Install:</b>

```r
# -----
# UMS IPSIT BIOINFO - Licence GPL-3
# https://github.com/fdumbioinfo/uem915
# title: installation uem915
# date: 12-03-2026
# -----
#
options(pkgType = "binary")
# annotation packages
if(!require("moalannotgene",quietly=TRUE)){install.packages("moalannotgene",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalannotensg",quietly=TRUE)){install.packages("moalannotensg",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalannotenst",quietly=TRUE)){install.packages("moalannotenst",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalannotensp",quietly=TRUE)){install.packages("moalannotensp",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbhs",quietly=TRUE)){install.packages("moalstringdbhs",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbmm",quietly=TRUE)){install.packages("moalstringdbmm",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbrn",quietly=TRUE)){install.packages("moalstringdbrn",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbdr",quietly=TRUE)){install.packages("moalstringdbdr",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbss",quietly=TRUE)){install.packages("moalstringdbss",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
# depend packages
if(!require("BiocManager",quietly=TRUE)){install.packages("BiocManager")}
if(!require("broom",quietly=TRUE)){BiocManager::install("broom",update=F)}
if(!require("dendextend",quietly=TRUE)){BiocManager::install("dendextend",update=F)}
if(!require("doParallel",quietly=TRUE)){BiocManager::install("doParallel",update=F)}
if(!require("dplyr",quietly=TRUE)){BiocManager::install("dplyr",update=F)}
if(!require("forcats",quietly=TRUE)){BiocManager::install("forcats",update=F)}
if(!require("foreach",quietly=TRUE)){BiocManager::install("foreach",update=F)}
if(!require("ggforce",quietly=TRUE)){BiocManager::install("ggforce",update=F)}
if(!require("ggpubr",quietly=TRUE)){BiocManager::install("ggpubr",update=F)}
if(!require("gplots",quietly=TRUE)){BiocManager::install("gplots",update=F)}
if(!require("ggrepel",quietly=TRUE)){BiocManager::install("ggrepel",update=F)}
if(!require("grDevices",quietly=TRUE)){BiocManager::install("grDevices",update=F)}
if(!require("gridExtra",quietly=TRUE)){BiocManager::install("gridExtra",update=F)}
if(!require("igraph",quietly=TRUE)){BiocManager::install("igraph",update=F)}
if(!require("parallel",quietly=TRUE)){BiocManager::install("parallel",update=F)}
if(!require("plyr",quietly=TRUE)){BiocManager::install("plyr",update=F)}
if(!require("Rgraphviz",quietly=TRUE)){BiocManager::install("Rgraphviz",update=F)}
if(!require("graphics",quietly=TRUE)){BiocManager::install("graphics",update=F)}
if(!require("stringr",quietly=TRUE)){BiocManager::install("stringr",update=F)}
if(!require("tidyselect",quietly=TRUE)){BiocManager::install("tidyselect",update=F)}
if(!require("utils",quietly=TRUE)){BiocManager::install("utils",update=F)}
if(!require("colourvalues",quietly=TRUE)){BiocManager::install("colourvalues",update=F)}
if(!require("fgsea",quietly=TRUE)){BiocManager::install("fgsea",update=F)}
if(!require("limma",quietly=TRUE)){BiocManager::install("limma",update=F)}
# moal package
if(!require("moal",quietly=TRUE)){install.packages("moal",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
#
```

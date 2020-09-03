#!/bin/bash

# script for installing BiocManager "by hand", as changes in Travis' code makes automatic
# installation impossible
# inspired by the code of this package: https://github.com/bimberlabinternal/OOSAP

mkdir $HOME/.R

export R_BIOC_VERSION='3.11'

# Verify config, add Bioconductor repos to default repos, so install.packages will pick them up.
if [ -e ~/.Rprofile ];then
    echo 'Rprofile exists:'
    cat ~/.Rprofile
fi

if [ -e ~/.Rprofile.site ];then
    echo 'Rprofile.site exists:'
    cat ~/.Rprofile.site
fi

Rscript -e "install.packages(c('devtools', 'BiocManager', 'remotes'), dependencies=TRUE, ask = FALSE)"

# See: https://stackoverflow.com/questions/26042751/cannot-install-package-xml-to-r
echo "options(repos = c(BiocManager::repositories(version = '${R_BIOC_VERSION}'), c('xml' = 'http://www.omegahat.net/R')));" >> ~/.Rprofile.site

# Log repos to ensure Bioconductor used:
echo 'R repos:'
Rscript -e "getOption('repos')"

# Bioconductor
Rscript -e "BiocManager::install(version='${R_BIOC_VERSION}')"

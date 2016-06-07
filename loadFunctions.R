# The following program loads all R utility functions stored in:
# ./utilities/Rcode/

utilitiesDir = './utilities/Rcode/';
rfiles = dir(utilitiesDir, pattern = '*.R$');
for (i in 1:length(rfiles)){
    loadPath = paste(utilitiesDir,rfiles[i],sep='');
    source(loadPath);
}


# installing/loading pandoc
pacs<-c("lubridate", "data.table", "roxygen2", "tidyverse", "data.table", "docstring", "here", 
        "installr", "MASS", "broom", "broom.mixed", "papaja", "msme", "kableExtra", "gridExtra", 
        "officer", "flextable", "knitr", "mclust", "factoextra", "cluster", "clValid", "jtools",
        "knitr", "stargazer", "gridExtra", "glmmTMB", "rticles", "xtable", "usethis","TRADER")
for(i in 1:length(pacs)){
  pac<-pacs[i]
  if(length((find.package(pac,quiet=TRUE)))==0){
    install.packages(pac)
  }else{
    print(paste0("Package ",pac, " was already installed"))
  }
}
#load pacs
invisible(lapply(pacs, function(x)require(x,character.only = T,quietly = T)))




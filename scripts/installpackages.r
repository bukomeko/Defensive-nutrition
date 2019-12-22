# installing/loading pandoc
pacs<-c("lubridate", "data.table", "roxygen2", "tidyverse", "data.table", "docstring", "here", 
        "installr", "MASS", "broom", "broom.mixed", "papaja", "msme", "kableExtra", "gridExtra", 
        "officer", "flextable", "knitr", "mclust", "factoextra", "cluster", "clValid", "jtools",
        "knitr", "stargazer", "gridExtra", "glmmTMB", "rticles", "xtable", "usethis")
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



# token: 6a41b6a8621490907762dddb203e62f0ce7f1dbd

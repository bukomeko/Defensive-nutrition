  # custom functions
  get_bpts            <- function(x) {
    q <- list()
    for (i in 1:length(x)) {
      p <- ifelse(x[i + 1] < x[i], x[i + 1] <- x[i], x[i + 1] <- x[i + 1])
      q[[i]] <- p
    }
    return(unlist(q))
  }
  
  get_dispersion      <- function(model) {
    rdf <- df.residual(model) # residual degrees of freedom
    rp <- residuals(model, type = "pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq / rdf
    pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
    c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
  }
  
  get_clustering      <- function(data, n) {
    intern <- scale(data) %>%
      clValid::clValid(data,
                       nClust = 2:n,
                       clMethods = c("hierarchical", "kmeans", "pam"),
                       validation = "internal"
      )
    return(summary(intern))
  }
  
  get_yield           <- function(data) {
    k <- -6.832
    a <- 0.755
    b <- 1.059
    c <- 0.569
    d <- 0.478
    data[, BWT := exp(k + (a * log(Girth0)) + (b * log(Girth100)) + (c * log(Hands)) + (d * (log(Fingers))))]
    write.csv(data, here::here("data", "temp", "Baseline_yield.csv"))
    return(data)
  }
  
  get_weevil_damage   <- function(data) {
    data %>%
      .[, XI := ((UI_damage + LI_damage) / 2)] %>%
      .[, XO := ((UO_damage + LO_damage) / 2)] %>%
      .[, XT := ((UI_damage + LI_damage + UO_damage + LO_damage) / 4)] -> data
    return(data)
  }
  
  get_insl            <- function(data, FL, YLS) {
    sl <- FL
    yls <- YLS
    insl <- ((yls-1)/sl)*100
    data[, I_nsl := insl]
    return(data)
  }
  
  get_necrosis        <- function(data,r1,r2,r3,r4,r5){
    nec <- r1+r2+r3+r4+r5
    data[,necrosis := nec]
    return(data)
  }
  
  get_density         <- function(data, d1, d2, d3, d4){
    data[, sqmeandist := (rowMeans(c(d1,d2,d3,d4),na.rm = TRUE)^2)]
    data[, mat_density:= 10000/sqmeandist]
  }

  get_leaf_area       <- function(data, Girth0, Height, FL) {
    a <- 0.404
    b <- 0.381
    c <- 0.411
    data <- data[, leaf_area := (a + (b*Height) + (c*Girth0)) * FL]
  }
  
  get_biomass <- function(Girth0, c = 0.085, a = 0.066){
    data[, AGB := exp(c + (a * log(Girth0)))] 
  }
  
  
  # Sorting variables of interest if Positively correlated
  org  <- function (x,y) {
    
    z  <- x [ order (y),]
    return (z)
  }
  
  # Function for calculation of Boundary Points
  bla  <- function (x) {# function for selecting the boundary points!!!
    b <- list()
    for (i in 1:length(x)) {
      a <-    ifelse (x [i+1] < x [i], x [i+1]  <-  x [i], 
                      x [i+1]  <-  x [i+1] )
      b [[i]] <- a    
    }
    return (b) 
  } 
  
  Boundline <- function(data,Fig.name){
    #'@title Boundline
    #'@param data The data containing the variables of interest
    #'@param Ymax is local maximum value of y variable
    # Selecting the Attainable Yield
    data$Ymax <-  max(data$RBWT,na.rm=T)
    
    # Calculating boundary points for positive correlations
    #   Use a loop for calculation boundary points 
    # find ways to generalise x and y variable names and provide them as arguments
    pf = which(!is.na(match(names(data), c("POT")))) # positive factors
    pf_names  <- c("POT")
    y_names = c("RBWTpts")
    
    m <- NULL
    n <- NULL
    L <- NULL 
    for (i in 1:length (pf)) {
      pf_i = data [,pf [i]] 
      data <- org(data,pf_i)
      
      m <-data$RBWT
      L <-m [1]
      m <- bla (data$RBWT )
      m[1] <- L
      n <- unlist(m)
      data$Y  <- n 
      names(data)[names(data) == "Y"] <- y_names[i]
    } 
    
    # Plots with the Boundary Points, we need to generalise the x and y variables
    y <- data$RBWTpts %>% na.omit()
    x <- data$POT
    summary(y)
    summary(x)
    
    Model <- broom::tidy(nls(y~1/(1+(K*exp(-R*x))), start=list(K=0.41504, R=0.1504168)))
    
    Model_summary <- summary(Model)
    
    data$Ypred_POT <- 1/(1+(Model$estimate[[1]]*exp(-Model$estimate[[2]]*x)))
    
    data$RBWTpts2=data$RBWTpts
    
    crit <- approx(y=data$POT, x=data$Ypred_POT, xout=0.9) # find the corresponding x value for which yield in 0.9
    
    bla_plot <-  ggplot (data = data, aes (x=POT, y=RBWTpts2)) +
      geom_point (aes(y=RBWT , colour="Y obs"),colour="black") + 
      geom_line(aes(y=Ypred_POT,colour="BP"),colour="red",lwd=1,lty=2)+
      theme_bw()+
      scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1))+
      scale_x_continuous(breaks=c(0,1,2,3,4))+ 
      geom_hline(yintercept = 0.9,lty=2,colour="grey")+
      geom_vline(xintercept = crit$y,lwd=0.75,lty=7,colour="grey")+
      labs(x = expression("K concentration"(cmol+kg^-1)),
           y = "Relative fresh bunch weight", title = dplyr::ensym(Fig.name))
    
return(bla_plot)
    
  }
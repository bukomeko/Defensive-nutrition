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
  
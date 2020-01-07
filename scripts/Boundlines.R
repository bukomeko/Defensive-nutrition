# Written by Dennis Ochola and Hannington Bukomeko
# functions used
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

Model <- nls(y~1/(1+(K*exp(-R*x))), start=list(K=0.41504, R=0.1504168))

Model_summary <<- summary(Model)

data$Ypred_POT <- 1/(1+(0.5216*exp(-1.035 *x)))

data$RBWTpts2=data$RBWTpts

approx(y=data$POT, x=data$Ypred_POT, xout=0.9) # find the corresponding x value for which yield in 0.9

bla_plot <-  ggplot (data = data, aes (x=POT, y=RBWTpts2)) +
  geom_point (aes(y=RBWT , colour="Y obs"),colour="black") + 
  geom_line(aes(y=Ypred_POT,colour="BP"),colour="red",lwd=1,lty=2)+
  theme_minimal()+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1))+
  scale_x_continuous(breaks=c(0,1,2,3,4))+
  labs(x = expression("K concentration"(cmol+kg^-1)), 
       y = "Relative fresh bunch weight", title = dplyr::ensym(Fig.name))+ 
  geom_hline(yintercept = 0.9,lty=2,colour="grey")

}

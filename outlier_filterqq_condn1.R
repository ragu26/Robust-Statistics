rm(list = ls())
setwd("E:/Mstat Courses/Robust Statistics/Assignment")
cat("\f")  #clear Console
library(ggplot2)
set.seed(12345)

nSample <- 2000
epsilon <- c(0.1, 0.2, 0.3, 0.4)


# function to generate normal dist
ge.norm <- function(gamma, nSample, epsilon) {
  
  # nSample <- 2000 epsilon <- c(0.1,0.2,0.3,0.4)
  ne <- length(epsilon)
  
  # Generating the Experiments with STD.NORMAL.DIST ####
  exp1 <- rep(NA, ne)
  exp2 <- rep(NA, ne)
  data.norm1 <- matrix(NA, nSample, ne)
  data.norm2 <- matrix(NA, nSample, ne)
  
  for (i in 1:4) {
    # 1st Experiment
    set.seed(12345)
    
    # sorting the data so as to contaminate
    data.norm <- sort(rnorm(nSample))
    exp1[i] <- nSample * epsilon[i]
    N1 <- nSample - exp1[i]
    # the data was contaminated by adding a constant gamma to the largest epsilon_n
    # data values.
    data.norm[N1:nSample] <- data.norm[N1:nSample] + gamma
    data.norm1[, i] <- data.norm
    # 2nd Experiment
    set.seed(12345)
    data.norm <- rnorm(nSample)
    exp2[i] <- (nSample * epsilon[i])/2
    N2 <- nSample - exp2[i]
    
    # we add the gamma value to the largest epsilon_n/2 data values while
    # simultaneously subtracting gamma from the smallest epsilon_n/2 data values
    data.norm[1:exp2[i]] <- data.norm[1:exp2[i]] - gamma
    data.norm[N2:nSample] <- data.norm[N2:nSample] + gamma
    data.norm2[, i] <- data.norm
  }
  return(list(data.norm1 = data.norm1, data.norm2 = data.norm2))
}
# function ends here

ge.cauchy <- function(gamma, nSample, epsilon) {
  
  # Generating the Experiments with CAUCHY.DIST ####
  
  ne <- length(epsilon)
  
  # Generating the Experiments with CAUCHY.DIST ####
  exp1 <- rep(NA, ne)
  exp2 <- rep(NA, ne)
  epsilon <- c(0.1, 0.2, 0.3, 0.4)
  data.cauchy1 <- matrix(NA, nSample, ne)
  data.cauchy2 <- matrix(NA, nSample, ne)
  for (i in 1:4) {
    # 1st Experiment
    set.seed(12345)
    data.cauchy <- sort(rcauchy(nSample))
    exp1[i] <- nSample * epsilon[i]
    N1 <- nSample - exp1[i]
    # the data was contaminated by adding a constant gamma to the largest epsilon_n
    # data values.
    data.cauchy[N1:nSample] <- data.cauchy[N1:nSample] + gamma
    data.cauchy1[, i] <- data.cauchy
    # 2nd Experiment
    set.seed(12345)
    data.cauchy <- rcauchy(nSample)
    exp2[i] <- (nSample * epsilon[i])/2
    N2 <- nSample - exp2[i]
    
    # we add the gamma value to the largest epsilon_n/2 data values while
    # simultaneously subtracting gamma from the smallest epsilon_n/2 data values
    data.cauchy[1:exp2[i]] <- data.cauchy[1:exp2[i]] - gamma
    data.cauchy[N2:nSample] <- data.cauchy[N2:nSample] + gamma
    data.cauchy2[, i] <- data.cauchy
  }
  return(list(data.cauchy1 = data.cauchy1, data.cauchy2 = data.cauchy2))
  
}
# function ends here

# The method flags all the data points for which |zi| > a and zi /qi > a, with a
# being a tuning constant.  In order to know what value for a is needed to have
# 1% of outliers in the uncontaminated data,

# tuning parameter for std normal dist
flag.tuning.norm <- function(temp.norm, nSample) {
  
  # tuning parameter values
  a1.norm <- matrix(0, 4, 1)
  
  # outlier info
  outlier.norm1 <- matrix(0, 4, 1)
  
  # tuning parameter values
  a2.norm <- matrix(0, 4, 1)
  
  # outlier info
  outlier.norm2 <- matrix(0, 4, 1)
  
  # dataset
  data.norm1 <- temp.norm$data.norm1
  data.norm2 <- temp.norm$data.norm2
  
  # calculate tuning parameter such that 1% of the outliers are detected in
  # uncontaminated data
  for (col in 1:4) {
    a = 1
    
    uncont = 1000
    
    # measure of location and scale
    m1 <- median(data.norm1[, col])
    s1 <- mad(data.norm1[, col])
    m2 <- median(data.norm2[, col])
    s2 <- mad(data.norm2[, col])
    
    i = 1:2000  #NSample is 2000
    
    # calculating robust median
    p = (i - 1/3)/(nSample + 1/3)
    q = qnorm(p)
    z1 = (data.norm1[, col] - m1)/s1
    z2 = (data.norm2[, col] - m2)/s2
    if (col == 1) {
      # for epsilon 0.1 and experiment 1
      uSample = nSample * 0.9
      repeat {
        uncont = sum((abs(z1[1:uSample]) > a) & ((z1[1:uSample]/q[1:uSample]) > 
                                                   a))
        if (uncont <= 18) 
          break
        a = a + 0.001
      }
      a1.norm[col] = a
      outlier.norm1[col] = sum((abs(z1[(uSample + 1):nSample]) > a) & ((z1[(uSample + 
                                                                              1):nSample]/q[(uSample + 1):nSample]) > a))
      
      # for epsilon 0.1 and experiment 2
      a = 1
      uSample = nSample * 0.95
      lSample = nSample * 0.05
      seq <- c(1:100, 1901:2000)
      repeat {
        uncont = sum((abs(z2[lSample:uSample]) > a) & ((z2[lSample:uSample]/q[lSample:uSample]) > 
                                                         a))
        if (uncont <= 18) 
          break
        a = a + 0.001
      }
      a2.norm[col] = a
      outlier.norm2[col] = sum((abs(z2[seq]) > a) & ((z2[seq]/q[seq]) > a))
    }
    if (col == 2) {
      # for epsilon 0.2 and experiment 1
      uSample = nSample * 0.8
      repeat {
        uncont = sum((abs(z1[1:uSample]) > a) & ((z1[1:uSample]/q[1:uSample]) > 
                                                   a))
        if (uncont <= 16) 
          break
        a = a + 0.001
      }
      a1.norm[col] = a
      outlier.norm1[col] = sum((abs(z1[(uSample + 1):nSample]) > a) & ((z1[(uSample + 
                                                                              1):nSample]/q[(uSample + 1):nSample]) > a))
      # for epsilon 0.2 and experiment 2
      a = 1
      uSample = nSample * 0.9
      lSample = nSample * 0.1
      seq <- c(1:200, 1801:2000)
      repeat {
        uncont = sum((abs(z2[lSample:uSample]) > a) & ((z2[lSample:uSample]/q[lSample:uSample]) > 
                                                         a))
        if (uncont <= 16) 
          break
        a = a + 0.001
      }
      a2.norm[col] = a
      outlier.norm2[col] = sum((abs(z2[seq]) > a) & ((z2[seq]/q[seq]) > a))
    }
    if (col == 3) {
      # for epsilon 0.3 and experiment 1
      uSample = nSample * 0.7
      repeat {
        uncont = sum((abs(z1[1:uSample]) > a) & ((z1[1:uSample]/q[1:uSample]) > 
                                                   a))
        if (uncont <= 14) 
          break
        a = a + 0.001
      }
      a1.norm[col] = a
      outlier.norm1[col] = sum((abs(z1[(uSample + 1):nSample]) > a) & ((z1[(uSample + 
                                                                              1):nSample]/q[(uSample + 1):nSample]) > a))
      # for epsilon 0.3 and experiment 2
      a = 1
      uSample = nSample * 0.85
      lSample = nSample * 0.15
      seq <- c(1:300, 1701:2000)
      repeat {
        uncont = sum((abs(z2[lSample:uSample]) > a) & ((z2[lSample:uSample]/q[lSample:uSample]) > 
                                                         a))
        if (uncont <= 14) 
          break
        a = a + 0.001
      }
      a2.norm[col] = a
      outlier.norm2[col] = sum((abs(z2[seq]) > a) & ((z2[seq]/q[seq]) > a))
    }
    if (col == 4) {
      # for epsilon 0.4 and experiment 1
      uSample = nSample * 0.6
      repeat {
        uncont = sum((abs(z1[1:uSample]) > a) & ((z1[1:uSample]/q[1:uSample]) > 
                                                   a))
        if (uncont <= 12) 
          break
        a = a + 0.001
      }
      a1.norm[col] = a
      outlier.norm1[col] = sum((abs(z1[(uSample + 1):nSample]) > a) & ((z1[(uSample + 
                                                                              1):nSample]/q[(uSample + 1):nSample]) > a))
      # for epsilon 0.4 and experiment 2
      a = 1
      uSample = nSample * 0.85
      lSample = nSample * 0.2
      seq <- c(1:400, 1601:2000)
      repeat {
        uncont = sum((abs(z2[lSample:uSample]) > a) & ((z2[lSample:uSample]/q[lSample:uSample]) > 
                                                         a))
        if (uncont <= 12) 
          break
        a = a + 0.001
      }
      a2.norm[col] = a
      outlier.norm2[col] = sum((abs(z2[seq]) > a) & ((z2[seq]/q[seq]) > a))
    }
  }
  return(list(a1.norm = a1.norm, a2.norm = a2.norm, outlier.norm1 = outlier.norm1, 
              outlier.norm2 = outlier.norm2))
}

# tuning parameter for cauchy dist
flag.tuning.cauchy <- function(temp.cauchy, nSample) {
  
  
  uncont = 1000
  a1.cauchy <- matrix(0, 4, 1)
  a2.cauchy <- matrix(0, 4, 1)
  data.cauchy1 <- temp.cauchy$data.cauchy1
  data.cauchy2 <- temp.cauchy$data.cauchy2
  outlier.cauchy1 <- matrix(0, 4, 1)
  outlier.cauchy2 <- matrix(0, 4, 1)
  for (col in 1:4) {
    a = 10
    
    # robust estimates of scale and location
    m1 <- median(data.cauchy1[, col])
    s1 <- mad(data.cauchy1[, col])
    m2 <- median(data.cauchy2[, col])
    s2 <- mad(data.cauchy2[, col])
    i = 1:2000
    
    # robust median, to be used as theoretical quantile
    p = (i - 1/3)/(nSample + 1/3)
    q = qnorm(p)
    
    # standardizing the data
    z1 = (data.cauchy1[, col] - m1)/s1
    z2 = (data.cauchy2[, col] - m2)/s2
    if (col == 1) {
      a = 10
      # for epsilon 0.1 and experiment 1
      uSample = nSample * 0.9
      repeat {
        uncont = sum((abs(z1[1:uSample]) > a) & ((z1[1:uSample]/q[1:uSample]) > 
                                                   a))
        if (uncont <= 18) 
          break
        a = a + 0.01
      }
      a1.cauchy[col] = a
      outlier.cauchy1[col] = sum((abs(z1[(uSample + 1):nSample]) > a) & ((z1[(uSample + 
                                                                                1):nSample]/q[(uSample + 1):nSample]) > a))
      # for epsilon 0.1 and experiment 2
      a = 10
      uSample = nSample * 0.95
      lSample = nSample * 0.05
      seq <- c(1:100, 1901:2000)
      repeat {
        uncont = sum((abs(z2[lSample:uSample]) > a) & ((z2[lSample:uSample]/q[lSample:uSample]) > 
                                                         a))
        if (uncont <= 18) 
          break
        a = a + 0.01
      }
      a2.cauchy[col] = a
      
      outlier.cauchy2[col] = sum((abs(z2[seq]) > a) & ((z2[seq]/q[seq]) > a))
      
    }
    if (col == 2) {
      a = 10
      # for epsilon 0.2 and experiment 1
      uSample = nSample * 0.8
      repeat {
        uncont = sum((abs(z1[1:uSample]) > a) & ((z1[1:uSample]/q[1:uSample]) > 
                                                   a))
        if (uncont <= 16) 
          break
        a = a + 0.01
      }
      a1.cauchy[col] = a
      outlier.cauchy1[col] = sum((abs(z1[(uSample + 1):nSample]) > a) & ((z1[(uSample + 
                                                                                1):nSample]/q[(uSample + 1):nSample]) > a))
      # for epsilon 0.2 and experiment 2
      a = 10
      uSample = nSample * 0.9
      lSample = nSample * 0.1
      seq <- c(1:200, 1801:2000)
      repeat {
        uncont = sum((abs(z2[lSample:uSample]) > a) & ((z2[lSample:uSample]/q[lSample:uSample]) > 
                                                         a))
        if (uncont <= 16) 
          break
        a = a + 0.01
      }
      a2.cauchy[col] = a
      
      outlier.cauchy2[col] = sum((abs(z2[seq]) > a) & ((z2[seq]/q[seq]) > a))
    }
    if (col == 3) {
      a = 10
      # for epsilon 0.3 and experiment 1
      uSample = nSample * 0.7
      repeat {
        uncont = sum((abs(z1[1:uSample]) > a) & ((z1[1:uSample]/q[1:uSample]) > 
                                                   a))
        if (uncont <= 14) 
          break
        a = a + 0.01
      }
      a1.cauchy[col] = a
      outlier.cauchy1[col] = sum((abs(z1[(uSample + 1):nSample]) > a) & ((z1[(uSample + 
                                                                                1):nSample]/q[(uSample + 1):nSample]) > a))
      # for epsilon 0.3 and experiment 2
      a = 10
      uSample = nSample * 0.85
      lSample = nSample * 0.15
      seq <- c(1:300, 1701:2000)
      repeat {
        uncont = sum((abs(z2[lSample:uSample]) > a) & ((z2[lSample:uSample]/q[lSample:uSample]) > 
                                                         a))
        if (uncont <= 14) 
          break
        a = a + 0.01
      }
      a2.cauchy[col] = a
      
      outlier.cauchy2[col] = sum((abs(z2[seq]) > a) & ((z2[seq]/q[seq]) > a))
    }
    if (col == 4) {
      a = 10
      # for epsilon 0.4 and experiment 1
      uSample = nSample * 0.6
      repeat {
        uncont = sum((abs(z1[1:uSample]) > a) & ((z1[1:uSample]/q[1:uSample]) > 
                                                   a))
        if (uncont <= 12) 
          break
        a = a + 0.01
      }
      a1.cauchy[col] = a
      outlier.cauchy1[col] = sum((abs(z1[(uSample + 1):nSample]) > a) & ((z1[(uSample + 
                                                                                1):nSample]/q[(uSample + 1):nSample]) > a))
      # for epsilon 0.4 and experiment 2
      a = 10
      uSample = nSample * 0.8
      lSample = nSample * 0.2
      seq <- c(1:400, 1601:2000)
      repeat {
        uncont = sum((abs(z2[lSample:uSample]) > a) & ((z2[lSample:uSample]/q[lSample:uSample]) > 
                                                         a))
        if (uncont <= 12) 
          break
        a = a + 0.01
      }
      a2.cauchy[col] = a
      
      outlier.cauchy2[col] = sum((abs(z2[seq]) > a) & ((z2[seq]/q[seq]) > a))
    }
  }
  return(list(a1.cauchy = a1.cauchy, a2.cauchy = a2.cauchy, outlier.cauchy1 = outlier.cauchy1, 
              outlier.cauchy2 = outlier.cauchy2))
}


# function to convert wide form of data to long form. inbuilt reshape function is
# not customized for the project needs
reshape.plot <- function(datax, gamma_val) {
  plot.data <- data.frame((cbind(gamma_val, datax)))
  melted <- matrix(0, 800, 3)
  
  # gamma sequence
  melted[1:200, 1] <- gamma_val[1:200]
  melted[201:400, 1] <- gamma_val[1:200]
  melted[401:600, 1] <- gamma_val[1:200]
  melted[601:800, 1] <- gamma_val[1:200]
  
  # contamination sequence
  melted[1:200, 2] <- rep(0.1, 200, 1)
  melted[201:400, 2] <- rep(0.2, 200, 1)
  melted[401:600, 2] <- rep(0.3, 200, 1)
  melted[601:800, 2] <- rep(0.4, 200, 1)
  
  # proportion of outliers
  melted[1:200, 3] <- datax[, 1]/nSample
  melted[201:400, 3] <- datax[, 2]/nSample
  melted[401:600, 3] <- datax[, 3]/nSample
  melted[601:800, 3] <- datax[, 4]/nSample
  melted[, 1] <- as.numeric(melted[, 1])
  melted[, 2] <- as.numeric(melted[, 2])
  melted[, 3] <- as.numeric(melted[, 3])
  melted <- data.frame(melted)
  
  # return set
  names(melted) <- c("gamma", "contamination", "outliers")
  return(list(melted = melted))
}

# finding pattern of tuning constant coressponding to gamma pattern
a.norm1 <- matrix(0, 200, 4)
a.norm2 <- matrix(0, 200, 4)
outlier.cont.norm1 <- matrix(0, 200, 4)
outlier.cont.norm2 <- matrix(0, 200, 4)
gamma_val <- matrix(0, 200, 1)
for (k in 1:200) {
  # gamma values from 0 to 20 incremnted by 0.1
  gamma_val[k] <- k * 0.1
  # experiments
  temp.norm <- ge.norm(gamma_val[k], nSample, epsilon)
  temp.a <- flag.tuning.norm(temp.norm, nSample)
  
  # store tuning constants
  a.norm1[k, 1] <- temp.a$a1.norm[1]
  a.norm1[k, 2] <- temp.a$a1.norm[2]
  a.norm1[k, 3] <- temp.a$a1.norm[3]
  a.norm1[k, 4] <- temp.a$a1.norm[4]
  a.norm2[k, 1] <- temp.a$a2.norm[1]
  a.norm2[k, 2] <- temp.a$a2.norm[2]
  a.norm2[k, 3] <- temp.a$a2.norm[3]
  a.norm2[k, 4] <- temp.a$a2.norm[4]
  
  # store outlier values in contaminated data i.e., true positives
  outlier.cont.norm1[k, 1] <- temp.a$outlier.norm1[1]
  outlier.cont.norm1[k, 2] <- temp.a$outlier.norm1[2]
  outlier.cont.norm1[k, 3] <- temp.a$outlier.norm1[3]
  outlier.cont.norm1[k, 4] <- temp.a$outlier.norm1[4]
  
  outlier.cont.norm2[k, 1] <- temp.a$outlier.norm2[1]
  outlier.cont.norm2[k, 2] <- temp.a$outlier.norm2[2]
  outlier.cont.norm2[k, 3] <- temp.a$outlier.norm2[3]
  outlier.cont.norm2[k, 4] <- temp.a$outlier.norm2[4]
}

# Plot Tuning Constant vs gamma for Standard Normal
longdata <- reshape.plot(outlier.cont.norm1, gamma_val)$melted
ggplot(longdata, aes(x = gamma, y = tuning_constant, group = contamination, colour = contamination), 
       main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + 
  xlab("Gamma") + ylab("Outlier Count for contaminated data") + ggtitle("Experiment 1 with Standard Normal data")

longdata <- reshape.plot(outlier.cont.norm2, gamma_val)$melted
ggplot(longdata, aes(x = gamma, y = tuning_constant, group = contamination, colour = contamination), 
       main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + 
  xlab("Gamma") + ylab("Outlier Count for contaminated data") + ggtitle("Experiment 2 with Standard Normal data")

# finding pattern of tuning constant coressponding to gamma pattern
outlier.cont.cauchy1 <- matrix(0, 200, 4)
outlier.cont.cauchy2 <- matrix(0, 200, 4)
a.cauchy1 <- matrix(0, 200, 4)
a.cauchy2 <- matrix(0, 200, 4)
gamma_val <- matrix(0, 200, 1)
for (k in 1:200) {
  
  # gamma values from 0 to 20 incremnted by 0.1
  gamma_val[k] <- k * 0.1
  # experiments
  temp.cauchy <- ge.cauchy(gamma_val[k], nSample, epsilon)
  temp.a <- flag.tuning.cauchy(temp.cauchy, nSample)
  a.cauchy1[k, 1] <- temp.a$a1.cauchy[1]
  a.cauchy1[k, 2] <- temp.a$a1.cauchy[2]
  a.cauchy1[k, 3] <- temp.a$a1.cauchy[3]
  a.cauchy1[k, 4] <- temp.a$a1.cauchy[4]
  a.cauchy2[k, 1] <- temp.a$a2.cauchy[1]
  a.cauchy2[k, 2] <- temp.a$a2.cauchy[2]
  a.cauchy2[k, 3] <- temp.a$a2.cauchy[3]
  a.cauchy2[k, 4] <- temp.a$a2.cauchy[4]
  
  # outlier.cont.cauchy1[k,1]<-temp.a$outlier.cauchy1[1]
  # outlier.cont.cauchy1[k,2]<-temp.a$outlier.cauchy1[2]
  # outlier.cont.cauchy1[k,3]<-temp.a$outlier.cauchy1[3]
  # outlier.cont.cauchy1[k,4]<-temp.a$outlier.cauchy1[4]
  # outlier.cont.cauchy2[k,1]<-temp.a$outlier.cauchy2[1]
  # outlier.cont.cauchy2[k,2]<-temp.a$outlier.cauchy2[2]
  # outlier.cont.cauchy2[k,3]<-temp.a$outlier.cauchy2[3]
  # outlier.cont.cauchy2[k,4]<-temp.a$outlier.cauchy2[4]
}

# #Gamma vs contaminated outliers #Plot Tuning Constant vs gamma for Cauchy
# longdata<-reshape.plot(outlier.cont.cauchy1,gamma_val)$melted ggplot(longdata,
# aes(x = gamma, y = tuning_constant,group=contamination, colour =
# contamination), main='Normal Dist') + geom_line()+
# scale_color_gradientn(colours = rainbow(4))+ xlab('Gamma')+ylab('Outlier count
# for contaminated data')+ggtitle('Experiment 1 with Cauchy data') #Plot Tuning
# Constant vs gamma for Cauchy
# longdata<-reshape.plot(outlier.cont.cauchy2,gamma_val)$melted ggplot(longdata,
# aes(x = gamma, y = tuning_constant,group=contamination, colour =
# contamination), main='Normal Dist') + geom_line()+
# scale_color_gradientn(colours = rainbow(4))+ xlab('Gamma')+ylab('Outlier count
# for contaminated data')+ggtitle('Experiment 2 with Cauchy data')


# Plotting tuning constants

# Plot Tuning Constant vs gamma for Standard Normal
longdata <- reshape.plot(a.norm1, gamma_val)$melted
ggplot(longdata, aes(x = gamma, y = tuning_constant, group = contamination, colour = contamination), 
       main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + 
  xlab("Gamma") + ylab("Tuning Constant") + ggtitle("Experiment 1 with Standard Normal data")

# Plot Tuning Constant vs gamma for Standard Normal
longdata <- reshape.plot(a.norm2, gamma_val)$melted
ggplot(longdata, aes(x = gamma, y = tuning_constant, group = contamination, colour = contamination), 
       main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + 
  xlab("Gamma") + ylab("Tuning Constant") + ggtitle("Experiment 2 with Standard Normal data")

# Plot Tuning Constant vs gamma for Cauchy
longdata <- reshape.plot(a.cauchy1, gamma_val)$melted
ggplot(longdata, aes(x = gamma, y = tuning_constant, group = contamination, colour = contamination), 
       main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + 
  xlab("Gamma") + ylab("Tuning Constant") + ggtitle("Experiment 1 with Cauchy data")

# Plot Tuning Constant vs gamma for Cauchy
longdata <- reshape.plot(a.cauchy2, gamma_val)$melted
ggplot(longdata, aes(x = gamma, y = tuning_constant, group = contamination, colour = contamination), 
       main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + 
  xlab("Gamma") + ylab("Tuning Constant") + ggtitle("Experiment 2 with Cauchy data")

# To exemplify the result given above a simulation was performed. We can see that
# when the value of the tuning constant is set equal to 2.3 and with a relatively
# small value of gamma, with which we can actually have outliers in the
# uncontaminated data, 1% of the uncontaminated data are flagged.  It is worth
# mentioning that the simulated data consists of 2000 data points and the
# simulation was done using the contamination procedure of the first experiment.
# Hence for the 10%,20%,30% and 40% of contaminated data we know how many
# outliers we should have, i.e. 18, 16, 14 and 12 respectively which seems to be
# satisfied. Of course when � becomes larger then, as mentioned before, it is
# easier to identify the outliers, therefore in the uncontaminated data we cannot
# find more outliers. 
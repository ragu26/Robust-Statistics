rm(list = ls())
cat("\014")  #clear Console
library(ggplot2)
set.seed(12345)


# Two experiments were conducted, with the difference between the experiments being the way
# the data was contaminated. In the first experiment, the data was contaminated by adding a
# constant gamma to the largest epsilon_n data values. In the second experiment, we add the
# gamma value to the largest epsilon_n/2 data values while simultaneously subtracting gamma
# from the smallest epsilon_n/2 data values, with epsilon being the fraction (or level) of
# contamination. Each experiment contains four levels of contamination of the data: 10, 20,
# 30 and 40 percent.
gamma <- 2

# total number of samples to be generated
nSample <- 2000

# fraction of contamination
epsilon <- c(0.1, 0.2, 0.3, 0.4)
ne <- length(epsilon)

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
    # the data was contaminated by adding a constant gamma to the largest epsilon_n data values.
    data.norm[N1:nSample] <- data.norm[N1:nSample] + gamma
    data.norm1[, i] <- data.norm
    # 2nd Experiment
    set.seed(12345)
    data.norm <- rnorm(nSample)
    exp2[i] <- (nSample * epsilon[i])/2
    N2 <- nSample - exp2[i]
    
    # we add the gamma value to the largest epsilon_n/2 data values while simultaneously
    # subtracting gamma from the smallest epsilon_n/2 data values
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
    # the data was contaminated by adding a constant gamma to the largest epsilon_n data values.
    data.cauchy[N1:nSample] <- data.cauchy[N1:nSample] + gamma
    data.cauchy1[, i] <- data.cauchy
    # 2nd Experiment
    set.seed(12345)
    data.cauchy <- rcauchy(nSample)
    exp2[i] <- (nSample * epsilon[i])/2
    N2 <- nSample - exp2[i]
    
    # we add the gamma value to the largest epsilon_n/2 data values while simultaneously
    # subtracting gamma from the smallest epsilon_n/2 data values
    data.cauchy[1:exp2[i]] <- data.cauchy[1:exp2[i]] - gamma
    data.cauchy[N2:nSample] <- data.cauchy[N2:nSample] + gamma
    data.cauchy2[, i] <- data.cauchy
  }
  return(list(data.cauchy1 = data.cauchy1, data.cauchy2 = data.cauchy2))
  
}
# function ends here

# function to detect outliers
flag_outlier <- function(datax) {
  
  # initialize indicator error flag
  ind <- matrix(0, nrow(datax), ncol(datax))
  # chisq with 99% probability and 1 df
  c <- sqrt(qchisq(0.9, 1))
  
  # loop to flag outliers for all the columns
  for (count in 1:ncol(datax)) {
    # location and scale for a particular column
    m <- median(datax[, count])
    s <- mad(datax[, count])
    # indicator storing outliers
    ind[, count] <- ifelse(abs((datax[, count] - m)/s) > c, TRUE, FALSE)
  }
  return(ind = ind)
  
}

# function to convert wide form of data to long form. inbuilt reshape function is not
# customized for the project needs
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

# Question 1 For Normal Dist;clea
gamma_val <- rep(0, 200)

# to store the proportion
outlier.uncont.norm1 <- matrix(0, 200, ne)
outlier.cont.norm1 <- matrix(0, 200, ne)
outlier.uncont.norm2 <- matrix(0, 200, ne)
outlier.cont.norm2 <- matrix(0, 200, ne)

# sequence of contamination in experiment 2
seq1 <- c(1:100, 1901:2000)
seq2 <- c(1:200, 1801:2000)
seq3 <- c(1:300, 1701:2000)
seq4 <- c(1:400, 1601:2000)

# collect the outlier data
for (k in 1:200) {
  # indicator to store if the particular cell is an outlier
  indq1.norm1 <- matrix(NA, nSample, ne)
  indq1.norm2 <- matrix(NA, nSample, ne)
  
  # series of gamma values from 0 to 20 incremented by 0.1
  gamma_val[k] <- k * 0.1
  # experiments
  temp.norm <- ge.norm(gamma_val[k], nSample, epsilon)
  data.norm1 <- temp.norm$data.norm1
  data.norm2 <- temp.norm$data.norm2
  
  # call function to flag outlier
  indq1.norm1 <- flag_outlier(data.norm1)
  indq1.norm2 <- flag_outlier(data.norm2)
  # outlier details outlier.norm1[k,] <- colSums(indq1.norm1) outlier.norm2[k,] <-
  # colSums(indq1.norm2)
  
  
  # summing the indicators for false positives
  outlier.uncont.norm1[k, 1] <- sum(indq1.norm1[1:1800, 1])
  outlier.uncont.norm1[k, 2] <- sum(indq1.norm1[1:1600, 2])
  outlier.uncont.norm1[k, 3] <- sum(indq1.norm1[1:1400, 3])
  outlier.uncont.norm1[k, 4] <- sum(indq1.norm1[1:1200, 4])
  outlier.cont.norm1[k, 1] <- sum(indq1.norm1[1801:2000, 1])
  outlier.cont.norm1[k, 2] <- sum(indq1.norm1[1601:2000, 2])
  outlier.cont.norm1[k, 3] <- sum(indq1.norm1[1401:2000, 3])
  outlier.cont.norm1[k, 4] <- sum(indq1.norm1[1201:2000, 4])
  
  # summing the indicators for true positives
  outlier.uncont.norm2[k, 1] <- sum(indq1.norm2[101:1900, 1])
  outlier.uncont.norm2[k, 2] <- sum(indq1.norm2[201:1800, 2])
  outlier.uncont.norm2[k, 3] <- sum(indq1.norm2[301:1700, 3])
  outlier.uncont.norm2[k, 4] <- sum(indq1.norm2[401:1600, 4])
  outlier.cont.norm2[k, 1] <- sum(indq1.norm2[seq1, 1])
  outlier.cont.norm2[k, 2] <- sum(indq1.norm2[seq2, 2])
  outlier.cont.norm2[k, 3] <- sum(indq1.norm2[seq3, 3])
  outlier.cont.norm2[k, 4] <- sum(indq1.norm2[seq4, 4])
}

library(gridExtra)
# Plot outlier vs gamma for Standard Normal convert wide to long
longdata.uncont <- reshape.plot(outlier.uncont.norm1, gamma_val)$melted
longdata.cont <- reshape.plot(outlier.cont.norm1, gamma_val)$melted
# store the plot
p1 <- ggplot(longdata.uncont, aes(x = gamma, y = outliers, group = contamination, colour = contamination), 
             main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + xlab("Gamma") + 
  ylab("Outlier Proportion") + ggtitle("Experiment 1 with Standard Normal Uncontaminated data")
p2 <- ggplot(longdata.cont, aes(x = gamma, y = outliers, group = contamination, colour = contamination), 
             main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + xlab("Gamma") + 
  ylab("Outlier Proportion") + ggtitle("Experiment 1 with Standard Normal Contaminated data")
grid.arrange(p1, p2)

# convert wide to long
longdata.uncont <- reshape.plot(outlier.uncont.norm2, gamma_val)$melted
longdata.cont <- reshape.plot(outlier.cont.norm2, gamma_val)$melted
p1 <- ggplot(longdata.uncont, aes(x = gamma, y = outliers, group = contamination, colour = contamination), 
             main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + xlab("Gamma") + 
  ylab("Outlier Proportion") + ggtitle("Experiment 2 with Standard Normal Uncontaminated data")
p2 <- ggplot(longdata.cont, aes(x = gamma, y = outliers, group = contamination, colour = contamination), 
             main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + xlab("Gamma") + 
  ylab("Outlier Proportion") + ggtitle("Experiment 2 with Standard Normal Contaminated data")
grid.arrange(p1, p2)

# Question 1 For Cauchy Dist;
outlier.uncont.cauchy1 <- matrix(0, 200, ne)
outlier.cont.cauchy1 <- matrix(0, 200, ne)
outlier.uncont.cauchy2 <- matrix(0, 200, ne)
outlier.cont.cauchy2 <- matrix(0, 200, ne)

# collect the outlier data
for (k in 1:200) {
  # indicator to store if the particular cell is an outlier
  indq1.cauchy1 <- matrix(NA, nSample, ne)
  indq1.cauchy2 <- matrix(NA, nSample, ne)
  # series of gamma values from 0 to 20 incremented by 0.1
  gamma_val[k] <- k * 0.1
  # experiments
  temp.cauchy <- ge.cauchy(gamma_val[k], nSample, epsilon)
  data.cauchy1 <- temp.cauchy$data.cauchy1
  data.cauchy2 <- temp.cauchy$data.cauchy2
  
  # call to flag outliers
  indq1.cauchy1 <- flag_outlier(data.cauchy1)
  indq1.cauchy2 <- flag_outlier(data.cauchy2)
  # outlier details
  
  # to count the number of false positives
  outlier.uncont.cauchy1[k, 1] <- sum(indq1.cauchy1[1:1800, 1])
  outlier.uncont.cauchy1[k, 2] <- sum(indq1.cauchy1[1:1600, 2])
  outlier.uncont.cauchy1[k, 3] <- sum(indq1.cauchy1[1:1400, 3])
  outlier.uncont.cauchy1[k, 4] <- sum(indq1.cauchy1[1:1200, 4])
  outlier.cont.cauchy1[k, 1] <- sum(indq1.cauchy1[1801:2000, 1])
  outlier.cont.cauchy1[k, 2] <- sum(indq1.cauchy1[1601:2000, 2])
  outlier.cont.cauchy1[k, 3] <- sum(indq1.cauchy1[1401:2000, 3])
  outlier.cont.cauchy1[k, 4] <- sum(indq1.cauchy1[1201:2000, 4])
  
  # to count the number of true positives
  outlier.uncont.cauchy2[k, 1] <- sum(indq1.cauchy2[101:1900, 1])
  outlier.uncont.cauchy2[k, 2] <- sum(indq1.cauchy2[201:1800, 2])
  outlier.uncont.cauchy2[k, 3] <- sum(indq1.cauchy2[301:1700, 3])
  outlier.uncont.cauchy2[k, 4] <- sum(indq1.cauchy2[401:1600, 4])
  outlier.cont.cauchy2[k, 1] <- sum(indq1.cauchy2[seq1, 1])
  outlier.cont.cauchy2[k, 2] <- sum(indq1.cauchy2[seq2, 2])
  outlier.cont.cauchy2[k, 3] <- sum(indq1.cauchy2[seq3, 3])
  outlier.cont.cauchy2[k, 4] <- sum(indq1.cauchy2[seq4, 4])
}

# Plot outlier vs gamma for Cauchy convert long to wide
longdata.uncont <- reshape.plot(outlier.uncont.cauchy1, gamma_val)$melted
longdata.cont <- reshape.plot(outlier.cont.cauchy1, gamma_val)$melted

# store plot
p1 <- ggplot(longdata.uncont, aes(x = gamma, y = outliers, group = contamination, colour = contamination), 
             main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + xlab("Gamma") + 
  ylab("Outlier Proportion") + ggtitle("Experiment 1 with Cauchy Uncontaminated data")
# store plot
p2 <- ggplot(longdata.cont, aes(x = gamma, y = outliers, group = contamination, colour = contamination), 
             main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + xlab("Gamma") + 
  ylab("Outlier Proportion") + ggtitle("Experiment 1 with Cauchy Contaminated data")
grid.arrange(p1, p2)

# convert long to wide
longdata.uncont <- reshape.plot(outlier.uncont.cauchy2, gamma_val)$melted
longdata.cont <- reshape.plot(outlier.cont.cauchy2, gamma_val)$melted
# store plot
p1 <- ggplot(longdata.uncont, aes(x = gamma, y = outliers, group = contamination, colour = contamination), 
             main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + xlab("Gamma") + 
  ylab("Outlier Proportion") + ggtitle("Experiment 2 with Cauchy Uncontaminated data")
# store plot
p2 <- ggplot(longdata.cont, aes(x = gamma, y = outliers, group = contamination, colour = contamination), 
             main = "Normal Dist") + geom_line() + scale_color_gradientn(colours = rainbow(4)) + xlab("Gamma") + 
  ylab("Outlier Proportion") + ggtitle("Experiment 2 with Cauchy Contaminated data")
grid.arrange(p1, p2) 
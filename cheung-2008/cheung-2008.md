The following set of analyses explore the article by Cheung (2008)

http://www.statmodel.com/download/MCheung.pdf

Cheung, M.W.L. (2008). A model for integrating fixed-, random-, and mixed-effects meta-analyses into structural equation modeling. Psychological Methods, 13, 182.

# Fixed effects meta-analysis
Let's create some data from the model for a fixed-effects meta-analysis.


```r
library(metafor)
```






```r
k <- 10
beta_fixed <- 0.2
meta_data <- data.frame(k = seq(k), sigma = runif(k, 0.05, 0.2))
meta_data$sigma_sq <- meta_data$sigma^2
meta_data$y <- rnorm(k, mean = beta_fixed, sd = meta_data$sigma)
```







```r
fixed_effects_meta_analysis <- function(y, sigma_sq) {
    # checks
    if (length(y) != length(sigma_sq)) {
        stop("length of y and sigma_sq should be equal")
    }
    if (any(sigma_sq <= 0)) {
        stop("non-positive sigma_sq: All sigma_sq must be positive")
    }
    
    
    k <- length(y)
    w <- 1/sigma_sq
    beta_hat_fixed <- sum(w * y)/sum(w)
    s_squared_hat_fixed <- 1/sum(w)
    s_hat_fixed <- sqrt(s_squared_hat_fixed)
    
    Z1 <- beta_hat_fixed/s_hat_fixed
    Z1_p_two_tail <- 1 - pnorm(abs(Z1))
    
    q <- w * (y - beta_hat_fixed)^2
    Q <- sum(q)
    Q_df <- k - 1
    Q_p <- 1 - pchisq(Q, df = Q_df)
    
    meta_table <- data.frame(y, sigma_sq, sigma = sqrt(sigma_sq), w, q)
    results <- list(meta_table = meta_table, beta_hat_fixed = beta_hat_fixed, 
        s_squared_hat_fixed = s_squared_hat_fixed, Z1 = c(Z1 = Z1, p_two_tail = Z1_p_two_tail), 
        Q = c(Q = Q, df = Q_df, p = Q_p))
    class(results) <- "fixed_effects_meta_analysis"
    results
}

print.fixed_effects_meta_analysis <- function(x, digits = 3, ...) {
    print(lapply(x, function(X) round(X, digits = digits)))
}

plot.fixed_effects_meta_analysis <- function(x, ...) {
    require(metafor)
    forest.default(x = x$meta_table$y, se = x$meta_table$sigma, refline = x$beta_hat_fixed, 
        ...)
}

meta1 <- fixed_effects_meta_analysis(y = meta_data$y, sigma_sq = meta_data$sigma_sq)
```




* The above code implements equations 2 to 5 in Cheung (2008).
* I created an S3 class. A `print` method is used to make display of the results a little easier to read. This is called by default when the return object printed to the console.
* A `plot` method produces a forest plot using the `forest.default` function in the `metafor` pakcage.



```r
plot(meta1)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

```r
meta1
```

```
## $meta_table
##        y sigma_sq sigma      w     q
## 1  0.119    0.005 0.069 209.07 1.378
## 2  0.152    0.022 0.147  45.99 0.109
## 3  0.328    0.005 0.072 194.87 3.148
## 4  0.192    0.019 0.139  51.88 0.004
## 5  0.144    0.005 0.067 222.01 0.717
## 6  0.226    0.006 0.079 160.13 0.103
## 7  0.231    0.006 0.077 167.80 0.153
## 8  0.137    0.008 0.087 131.76 0.531
## 9  0.380    0.009 0.096 108.67 3.503
## 10 0.156    0.003 0.059 286.32 0.567
## 
## $beta_hat_fixed
## [1] 0.201
## 
## $s_squared_hat_fixed
## [1] 0.001
## 
## $Z1
##         Z1 p_two_tail 
##      7.968      0.000 
## 
## $Q
##      Q     df      p 
## 10.213  9.000  0.334 
## 
```




# Models with covariates


```r
k <- 10
beta <- c(0.2, 0.3)
meta_data_covariates <- data.frame(k = seq(k), sigma = runif(k, 0.05, 
    0.2), x1 = c(rep(1, k/2), rep(0, k/2)), x2 = c(rep(0, k/2), rep(1, k/2)))
meta_data_covariates$sigma_sq <- meta_data$sigma^2
meta_data_covariates$beta <- as.matrix(meta_data_covariates[, c("x1", 
    "x2")]) %*% beta
meta_data_covariates$y <- rnorm(k, mean = meta_data_covariates$beta, 
    sd = meta_data_covariates$sigma)
```






```r
fixed_effects_covariates_meta_analysis <- function(y, X, sigma_sq) {
    # y: vector of effect sizes of length k, X: design matrix (k rows by p
    # columns) Ve: vector of effect size variances of length k
    p <- ncol(X)
    k <- length(y)
    
    Ve <- diag(sigma_sq)
    beta_hat <- solve(t(X) %*% solve(Ve) %*% X) %*% t(X) %*% solve(Ve) %*% y
    
    V_beta_hat <- solve(t(X) %*% solve(Ve) %*% X)
    
    Q_tilde <- t(beta_hat) %*% solve(V_beta_hat) %*% beta_hat
    Q_df <- p - 1
    Q_p <- 1 - pchisq(Q_tilde, df = Q_df)
    
    # significance test for H0: Beta_i = 0
    beta_hat <- as.vector(beta_hat)
    Z2 <- sapply(seq(beta_hat), function(X) beta_hat[X]/sqrt(V_beta_hat[X, X]))
    Z2_p_two_tail <- 1 - pnorm(abs(Z2))
    
    
    results <- list(meta_table = data.frame(y, X), V_beta_hat = V_beta_hat, 
        beta_hat = data.frame(beta_hat, z = Z2, p_two_tailed = Z2_p_two_tail), 
        Q = c(Q = Q_tilde, df = Q_df, p = Q_p))
    class(results) <- c("fixed_effects_covariates_meta_analysis", "fixed_effects_meta_analysis")
    results
    
}


fixed_effects_covariates_meta_analysis(y = meta_data_covariates$y, 
    X = as.matrix(meta_data_covariates[, c("x1", "x2")]), sigma_sq = meta_data_covariates$sigma_sq)
```

```
## $meta_table
##        y x1 x2
## 1  0.181  1  0
## 2  0.085  1  0
## 3  0.257  1  0
## 4  0.039  1  0
## 5  0.177  1  0
## 6  0.170  0  1
## 7  0.202  0  1
## 8  0.273  0  1
## 9  0.418  0  1
## 10 0.240  0  1
## 
## $V_beta_hat
##       x1    x2
## x1 0.001 0.000
## x2 0.000 0.001
## 
## $beta_hat
##   beta_hat     z p_two_tailed
## 1    0.184 4.950            0
## 2    0.247 7.227            0
## 
## $Q
##     Q    df     p 
## 76.74  1.00  0.00 
## 
```






# A modified version of findbeta ----------------------------------------------
findbeta2 <- function (themean = NULL, themedian = NULL, themode = NULL, percentile = 0.95, 
          lower.v = F, percentile.value) 
{
  stopifnot((is.null(themean) + is.null(themedian) + is.null(themode)) == 
              2)
  if (is.null(themode) && is.null(themedian)) {
    stopifnot((lower.v == T && themean <= percentile.value) | 
                (lower.v == F && themean >= percentile.value))
  }
  else if (is.null(themean) && is.null(themode)) {
    stopifnot((lower.v == T && themedian <= percentile.value) | 
                (lower.v == F && themedian >= percentile.value))
  }
  else {
    stopifnot((lower.v == T && themode <= percentile.value) | 
                (lower.v == F && themode >= percentile.value))
  }
  a = runif(1, 1, 10)
  if (lower.v == T) {
    pr_n = percentile
  }
  else {
    pr_n = 1 - percentile
  }
  if (is.null(themode) && is.null(themedian)) {
    to.minimize <- function(a) {
      abs(qbeta(pr_n, shape1 = a, shape2 = a * (1 - themean)/themean) - 
            percentile.value)
    }
  }
  else if (is.null(themean) && is.null(themode)) {
    to.minimize <- function(a) {
      abs(qbeta(pr_n, shape1 = a, shape2 = (3 * a * (1 - 
                                                       themedian) + 2 * themedian - 1)/(3 * themedian)) - 
            percentile.value)
    }
  }
  else {
    to.minimize <- function(a) {
      abs(qbeta(pr_n, shape1 = a, shape2 = (a * (1 - themode) + 
                                              2 * themode - 1)/themode) - percentile.value)
    }
  }
  estimate <- optim(runif(1, 1, 10), to.minimize, lower = 0.1, 
                    upper = 10^4, method = "Brent")
  finalshape1 = estimate$par
  if (is.null(themode) && is.null(themedian)) {
    finalshape2 = finalshape1 * (1 - themean)/themean
  }
  else if (is.null(themean) && is.null(themode)) {
    finalshape2 = (3 * finalshape1 * (1 - themedian) + 2 * 
                     themedian - 1)/(3 * themedian)
  }
  else {
    finalshape2 = (finalshape1 * (1 - themode) + 2 * themode - 
                     1)/themode
  }
  shapes <- c(round(finalshape1, 2), round(finalshape2, 2))
  
  return(shapes)
}

# A modified version of findbetaqq --------------------------------------------

findbetaqq2 <- function (percentile.value1, percentile1, percentile.value2, percentile2) 
{
  findcentiles <- function(x) {
    c(F1 = qbeta(percentile1, x[1], x[2]) - percentile.value1, 
      F2 = qbeta(percentile2, x[1], x[2]) - percentile.value2)
  }
  
  ss <- multiroot(f = findcentiles, start = c(1, 1))
  finalshape1 = ss$root[1]
  finalshape2 = ss$root[2]
  sample_beta = rbeta(10000, finalshape1, finalshape2)
  shapes <- c(round(finalshape1, 2), round(finalshape2, 2))
  
  return(shapes)
}

# Test
findbetaqq2(percentile.value1=0.30, percentile1=0.20, percentile.value2=0.60, percentile2=0.90)

# logit / inv_logit-----

logit <- function(x) log(x/(1-x))
inv_logit <- function(x) exp(x)/(1+exp(x))

# normal prior on logit Se/Sp ----

normal_prior <- function(lo95 = NULL, hi95 = NULL) {
  a <-
    matrix(c(1, 1, -1.96, 1.96),
           nrow = 2)
  b <-
    matrix(logit(c(lo95, hi95)),
           nrow = 2)
  d <- as.vector(solve(a, b))
  names(d) <- c("logit_mu", "logit_se")
  
  return(d)
}

# Test
parms <- normal_prior(lo95=0.95, hi95=0.99)
parms

bm_nondif <- " model {
for (i in 1:m) {
                # likelihood
                y[i,1] ~ dbin(p[i], N[i])
                p[i] <- tp[i]*Se + (1-tp[i])*(1-Sp)
                }

A <- N[1] * tp[1]
C <- N[1] * (1-tp[1])
B <- N[2] * tp[2]
D <- N[2] * (1-tp[2])

# priors for prevalence parameters
                tp[1] ~ dbeta(1, 1)
                tp[2] ~ dbeta(1, 1)
                
# priors for sens and spec
logit_Se ~ dnorm(logit_Se_mu, logit_Se_tau0)
Se <- exp(logit_Se)/(1+exp(logit_Se))

logit_Sp ~ dnorm(logit_Sp_mu, logit_Sp_tau0)
Sp <- exp(logit_Sp)/(1+exp(logit_Sp))

RR <- tp[1]/tp[2]
RD <- tp[1]-tp[2]
OR <- (tp[1]/(1-tp[1])) / (tp[2]/(1-tp[2]))
logRR <- log(RR)
logOR <- log(OR)

#data# m, N, y, logit_Se_mu, logit_Sp_mu, logit_Se_tau0, logit_Sp_tau0
#inits#
#monitor# tp, RR, RD, OR, logRR, logOR
}
"

bm_perf <- " model {
for (i in 1:m) {
                # likelihood
                y[i,1] ~ dbin(p[i], N[i])
                p[i] <- tp[i]
                }

A <- N[1] * tp[1]
C <- N[1] * (1-tp[1])
B <- N[2] * tp[2]
D <- N[2] * (1-tp[2])

# priors for prevalence parameters
                tp[1] ~ dbeta(1, 1)
                tp[2] ~ dbeta(1, 1)

RR <- tp[1]/tp[2]
RD <- tp[1]-tp[2]
OR <- (tp[1]/(1-tp[1])) / (tp[2]/(1-tp[2]))
logRR <- log(RR)
logOR <- log(OR)

#data# m, N, y
#inits#
#monitor# Se, Sp, tp, RR, RD, OR, logRR, logOR, A, B, C, D
}
"


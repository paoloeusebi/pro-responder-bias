library(tidyverse)
library(runjags)
library(rjags)
library(rootSolve)
# library(kableExtra)
testjags()
source("functions.R")

# prior calibration (manually)
logit_Se_mu <- logit(0.97)
logit_Se_sd <-  0.27
logit_Se_tau0 <- 1/(logit_Se_sd^2)
logit_Se <- rnorm(100000, mean=logit_Se_mu, sd=logit_Se_sd)
Se <- exp(logit_Se)/(1+exp(logit_Se))
quantile(Se, c(0.025, 0.975))
#hist(Se)

logit_Sp_mu <- logit(0.88)
logit_Sp_sd <- 0.32
logit_Sp_tau0 <- 1/(logit_Sp_sd^2)
logit_Sp <- rnorm(100000, mean=logit_Sp_mu, sd=logit_Sp_sd)
Sp <- exp(logit_Sp)/(1+exp(logit_Sp))
quantile(Sp, c(0.025, 0.975))
#hist(Sp)

y <- c(54, 58, 
       16, 99)
y <- matrix(y, nrow = 2, byrow = T); y

N <- apply(y, 1, sum); N

m <- 2


## initial values
inits1 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 100022)
inits2 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 300022)
inits3 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 500022)

res_bm <- run.jags(
  bm_nondif,
  burnin = 50000, 
  sample = 150000,
  n.chains = 3,
  inits = list(inits1, inits2, inits3)
)
res_bm

plot(
  res_bm,
  plot.type = c("h","t", "au"),
  vars = c("tp","logRR", "RR"),
  layout = c(4, 3)
)


res_bm_p <- run.jags(
  bm_perf,
  burnin = 50000, 
  sample = 150000,
  n.chains = 3,
  inits = list(inits1, inits2, inits3)
)
res_bm_p

plot(
  res_bm_p,
  plot.type = c("h","t", "au"),
  vars = c("tp","logRR", "RR"),
  layout = c(4, 3)
)


library(ggpubr)
library(scales)
theme_set(theme_minimal())
theme_update(axis.title.y = element_text(size=12),
             plot.background = element_rect(fill="#262626", colour = NA),
             panel.background = element_rect(fill="#262626", colour = NA),
             panel.grid = element_blank(),
             axis.text = element_text(colour = "white"),
             text = element_text(colour = "white", family="Lato")
)

pt <- ggplot(data = tibble(x = 0:1), aes(x)) +
  stat_function(fun = dbeta, n = 101, args = list(1, 1), col = "white") +
  labs(title = "Response in Sildenafil", 
       y = "density", x = "") +
  scale_y_continuous(labels = label_number(accuracy = 0.01),
                     limits = c(0, 1))
pt
ggsave(filename = "output/priors_pt.png", plot = pt, width = 3, height = 3)

pc <- ggplot(data = tibble(x = 0:1), aes(x)) +
  stat_function(fun = dbeta, n = 101, args = list(1, 1), col = "white") +
  labs(title = "Response in Placebo", 
       y = "density", x = "") +
  scale_y_continuous(labels = label_number(accuracy = 0.01),
                     limits = c(0, 1))
pc
ggsave(filename = "output/priors_pc.png", plot = pc, width = 3, height = 3)

se <- ggplot(data = tibble(x = c(0.5, 1)), aes(x)) +
  stat_function(fun = function(x) (dnorm(logit(x), mean = logit_Se_mu, sd=logit_Se_sd)), n = 1e06, col = "white") +
  labs(title = "Threshold sensitivity", 
       y = "density", x = "") +
  scale_y_continuous(labels = label_number(accuracy = 0.01))
se
ggsave(filename = "output/priors_se.png", plot = se, width = 3, height = 3)

sp <- ggplot(data = tibble(x = c(0.5,1)), aes(x)) +
  stat_function(fun = function(x) inv_logit(dnorm(logit(x), mean = logit_Sp_mu, sd=logit_Sp_sd)), n = 1e06, col = "white") +
  labs(title = "Threshold specificity", 
       y = "density", x = "") +
  scale_y_continuous(labels = label_number(accuracy = 0.01))
sp
ggsave(filename = "output/priors_sp.png", plot = sp, width = 3, height = 3)

priors <- ggarrange(pt, pc, se, sp,
                    # labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
ggsave(filename = "output/priors.png", plot = priors, width = 6, height = 6)

d_p <- res_bm_p$mcmc[[3]] %>% 
  as.data.frame() %>%
  mutate(logRR=log(RR))
d_b <- res_bm$mcmc[[3]] %>% 
  as.data.frame() %>%
  mutate(logRR=log(RR))

names(d_b)[3] <- "tp_1"
names(d_b)[4] <- "tp_2"
d_b <- d_b %>% 
  mutate(prior_logRR=log(tp_1/tp_2))

posterior <- ggplot() +
  geom_histogram(data = d_p, aes(x=logRR), col="grey90", fill="grey95", alpha=0.5) +
  geom_histogram(data = d_b, aes(x=logRR), col="#5B9BD5", fill="#5B9BD5", alpha=0.5, ) +
  labs(title = "",
       xlab = "log RR") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  scale_x_continuous(limits = c(0, 7.5))
posterior
ggsave(filename = "output/posterior.png", plot = sp, width = 4, height = 4)

post <- ggplot() +
  geom_density(data = d_b, aes(x=logRR), col="#5B9BD5", fill="#5B9BD5", alpha=0.75) +
  geom_density(data = d_p, aes(x=logRR), col="white", fill="grey90", alpha=0.25) +
  labs(title = "",
       x = "log RR") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  scale_x_continuous(limits = c(0, 10))
posterior
ggsave(filename = "output/posterior.png", plot = post, width = 4, height = 4)

d_p <- res_bm_p$mcmc[[3]] %>% 
  as.data.frame() %>%
  mutate(logRR=log(RR))
d_b <- res_bm$mcmc[[3]] %>% 
  as.data.frame() %>%
  mutate(logRR=log(RR))

d_pp <- data.frame("x" = c(0 , 1),
                   "y" = c(res_bm$summaries[c("RR"),c("Median")],
                           res_bm_p$summaries[c("RR"),c("Median")]),
                   "ymin" = c(res_bm$summaries[c("RR"),c("Lower95")],
                              res_bm_p$summaries[c("RR"),c("Lower95")]),
                   "ymax" = c(res_bm$summaries[c("RR"),c("Upper95")],
                              res_bm_p$summaries[c("RR"),c("Upper95")]))
post2 <- ggplot() +
  geom_pointrange(data = d_pp, 
                  aes(x=x, y=y, ymin=ymin, ymax=ymax))
post2
ggsave(filename = "output/posterior.png", plot = post, width = 4, height = 4)


ggplot() +
  geom_histogram(data = d_b, aes(x=RR), col="#5B9BD5", fill="#5B9BD5", alpha=0.5) +
  geom_histogram(data = d_p, aes(x=logRR), col="white", fill="grey95", alpha=0.5) +
  labs(title = "") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  scale_x_continuous(limits = c(0, 5))




posterior <- ggplot() +
  geom_density(data = d_b, aes(x=logRR), col="blue") +
  #geom_density(data = d_p, aes(x=logRR), col="black") +
  # geom_density(data = d_f, aes(x=logRR), col="blue") +
  labs(title = "") +
  scale_y_continuous(labels = label_number(accuracy = 0.01))
fig2






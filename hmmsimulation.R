library(tileHMM); library(lattice); library(ggplot2); library(caret)
pinull <- 0.05
pisignal <- 0.15
tau <- 1
sigma <- 1
state.names <- c("null","signal")
transition <- c(pinull,pisignal)
location <- c(0,tau)
scale <- c(sigma,sigma)
df <- c(Inf,Inf)

hmm <- getHMM(list(a = transition, mu = location, sigma = scale, nu = df), state.names)
set.seed(1994)
obs <- sampleSeq(hmm,1000, return.states = T)
states <- obs$states
obs <- obs$observation

n2n <- 0; n2s <- 0; s2n <- 0; s2s <- 0
for(i in 1:(length(states)-1)) {
    if(states[i] == 'null' & states[i+1] == 'null') {
        n2n <- n2n + 1
    } else if(states[i] == 'null' & states[i+1] == 'signal') {
        n2s <- n2s + 1
    } else if(states[i] == 'signal' & states[i+1] == 'null') {
        s2n <- s2n + 1
    } else {
        s2s <- s2s + 1
    }
}

reject <- function(pt,t){
  return(ifelse(pt > t, "signal","null"))
}

simulatemethod <- function(n,t = 0.95){
  vec <- character(n)
  for(i in 1:1000){
      seq <- obs[1:i]
      posteriorvec <- posterior(seq,hmm, log = F)
      pt <- posteriorvec[2,i]
      vec[i] <- reject(pt,t)
    }
  return(vec)
}

# calculates P(s_t = 0 | z_1, ..., z_(t-1))
# note P(s_t = 1 | z_1, ..., z_(t-1)) = 
# 1 - P(s_t = 0 | z_1, ..., z_(t-1))
p_next_null = function(posteriorm1){
  return(posteriorm1[1] * pinull + posteriorm1[2] * (1 - pisignal))
}

classifications <- simulatemethod(1000)
#Doesn't control overall error rate
errorrate <- mean(classifications != states)

plot(0, xlim=c(0, length(errorLocations)), type='n')
errorLocations <- classifications != states
for(i in 1:length(errorLocations)) {
    if(errorLocations[i] == TRUE)
        abline(v=i)
}
confusionMatrix(classifications, states)

table(classifications)

#Rina's Code
alpha_inv = function(P,alpha){
  # alpha-investing by Foster & Stine (2007)
  # W(0) = alpha
  # at time i reject if P[i]<=alpha[i]
  # choose alpha[i]<=W(i-1)/(1+W(i-1))
  # if reject, W(i) = W(i-1) + alpha
  # if not, W(i) = W(i-1) - alpha[i]/(1-alpha[i])
  
  # scheme: choose alpha[i] = 0.5*W(i-1)/(1+W(i-1))
  
  n = length(P)
  W = rep(0,n+1)
  alphas = rep(0,n)
  W[1] = alpha
  discoveries = rep(0,n)
  
  for(i in 1:n){
    alphas[i] = 0.5 * W[i] / (1+W[i])
    if(P[i]<=alphas[i]){
      discoveries[i] = 1
      W[i+1] = W[i]+alpha
    }else{
      W[i+1] = W[i] - alphas[i]/(1-alphas[i])
    }		
  }
  output = list()
  output$discoveries = discoveries
  output$alphas = alphas
  return(output)
}

LBOND = function(P,alpha){
  # "level based on number of discoveries", Javanmard & Montanari 2015
  # at time i reject if P[i]<=alpha[i]
  # set alpha[i] = beta[i]*max{1,D(i-1)}
  # where D(i-1) is the number of discoveries up to time i-1
  # and beta is a sequence with sum_{i=1}^{infty} beta[i] = alpha
  
  # scheme: beta[i] proportional to i^{-1.5}
  
  n = length(P)
  alphas = rep(0,n)
  beta = (1:n)^(-1.5); beta = beta * alpha / sum(beta)
  discoveries = rep(0,n)
  ndisc = 0
  
  for(i in 1:n){
    alphas[i] = beta[i] * max(1,ndisc)
    if(P[i]<=alphas[i]){
      discoveries[i] = 1
      ndisc = ndisc + 1
    }
  }
  
  output = list()
  output$discoveries = discoveries
  output$alphas = alphas
  return(output)	
}
LBORD = function(P,alpha){
  # "level based on recent discoveries", Javanmard & Montanari 2015
  # at time i reject if P[i]<=alpha[i]
  # set alpha[i] = beta[i-tau(i)]
  # where tau(i) is the time of the most recent discovery
  # and beta is a sequence with sum_{i=1}^{infty} beta[i] = alpha
  
  # scheme: beta[i] proportional to i^{-1.5}
  
  n = length(P)
  alphas = rep(0,n)
  beta = (1:n)^(-1.5); beta = beta * alpha / sum(beta)
  discoveries = rep(0,n)
  tau = 0
  
  for(i in 1:n){
    alphas[i] = beta[i - tau]
    if(P[i]<=alphas[i]){
      discoveries[i] = 1
      tau = i
    }
  }
  
  output = list()
  output$discoveries = discoveries
  output$alphas = alphas
  return(output)	
}

#Alpha investing only doing one sided testing
#This might not be the right thing to do
pvals <- rep(0,1000)
pvals[states == "signal"] <- 1 - pnorm(obs[states == "signal"])
pvals[states == "null"] <- 1 - pnorm(obs[states == "null"])

alphainv <- alpha_inv(pvals,0.10)
#Only makes one discovery, which is obviously not that good
alphainv$discoveries

lbond <- LBOND(pvals,0.10)
lbond$discoveries
#However, this simulation doesn't really model the use case we have in mind. 

lbord <- LBORD(pvals,0.10)
lbord$discoveries


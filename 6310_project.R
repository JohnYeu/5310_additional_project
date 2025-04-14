
# number of sample points
n1 <- 26 # complete data
n2 <- 15 # missing B
n3 <- 17 # missing A
n = n1 + n2 + n3 # total number

# index: p=complete, q=missing B, r=missing A
p <- 1:n1  # [1,2...26]
q <- (n1+1):(n1+n2)  # [27...41]
r <- (n1+n2+1):(n1+n2+n3)  # [42...58]


# initial guess for parameters
dataA = interexp[c(p,q),1]
dataB = interexp[c(p,r),2]
lenA = length(dataA)
lenB = length(dataB)
muA.old = mean(dataA)
muB.old = mean(dataB)
sigma11.old = var(dataA) * (lenA - 1) / lenA
sigma22.old = var(dataB) * (lenB - 1) / lenB
sigma12.old = cov(interexp[p,1], interexp[p,2])
sigma21.old = cov(interexp[p,1], interexp[p,2])


# EM algorithm
epsilon = 1e-6
iter <- 0
repeat {
  #complete data
  mean1.1 <- mean(interexp[p,1])
  mean1.2 <- mean(interexp[p,2])
  sig1 <- cov(interexp[p,]) * (n1 - 1)
  
  # data missing B
  mean2.1 <- mean(interexp[q,1])
  z.2 <- (muA.old + sigma12.old*(interexp[q,1] - mean2.1)/sigma11.old)
  z2.2 <- (sigma22.old - sigma12.old^2/sigma11.old) + z.2^2
  
  # data missing A
  mean3.2 <- mean(interexp[r,2])
  z.1 <- (muB.old + sigma12.old*(interexp[r,2] - mean3.2)/sigma22.old)
  z2.1 <- (sigma11.old - sigma12.old^2/sigma22.old) + z.1^2
  
  # M-step: update muA, muB
  muA.new <- mean(c(interexp[p,1], interexp[q,1], z.1))
  muB.new <- mean(c(interexp[p,2], z.2, interexp[r,2]))
  #print(length(c(interexp[p,1], z.1, interexp[r,1])))
  
  # M-step: update sigma
  sig2 <- matrix(c(sum((interexp[q,1] - muA.new)*(interexp[q,1] - muA.new)),
                   rep(sum((z.2 - muB.new)*(interexp[q,1] - muA.new)),2),
                   sum(z2.2 - 2*muB.new*z.2 + muB.new^2)), nrow=2, ncol=2)
  sig3 <- matrix(c(sum(z2.1 - 2*muA.new*z.1 + muA.new^2),
                   rep(sum((z.1 - muA.new)*(interexp[r,2] - muB.new)),2),
                   sum((interexp[r,2] - muB.new)*(interexp[r,2] - muB.new))), nrow=2, ncol=2)
  sigma11.new <-((sig1 + sig2 + sig3)/n)[1,1]
  sigma12.new <-((sig1 + sig2 + sig3)/n)[1,2]
  sigma21.new <-((sig1 + sig2 + sig3)/n)[2,1]
  sigma22.new <-((sig1 + sig2 + sig3)/n)[2,2]
  iter <- iter + 1
  
  d.muA <- max(abs(muA.new - muA.old))
  d.muB <- max(abs(muB.new - muB.old))
  d.sigma11 <- max(abs(sigma11.new - sigma11.old))
  d.sigma12 <- max(abs(sigma12.new - sigma12.old))
  d.sigma21 <- max(abs(sigma21.new - sigma21.old))
  d.sigma22 <- max(abs(sigma22.new - sigma22.old))
  
  if(max(d.muA, d.muB, d.sigma11, d.sigma12, d.sigma21, d.sigma22) < epsilon || iter > 500) break else {
    muA.old <- muA.new
    muB.old <- muB.new
    sigma11.old <- sigma11.new
    sigma12.old <- sigma12.new
    sigma21.old <- sigma21.new
    sigma22.old <- sigma22.new
  }
}
cat("iter = ", iter, "\n")
cat("muA = ", muA.old, "\n")
cat("muB = ", muB.old, "\n")
cat("sigma11 = ", sigma11.old, "\n")
cat("sigma12 = ", sigma12.old, "\n")
cat("sigma21 = ", sigma21.old, "\n")
cat("sigma22 = ", sigma22.old, "\n")

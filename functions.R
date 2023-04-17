SumAllComb <- function(x, m, n){
  # this function is to get sum of products of m elements from x
  # x is vector 
  # n is length(x)
  if(m==0){
    return(1)
  }else{
    tmp <- rep(1, (n-m+1))
    for(k in 1:m){
      tmp <- cumsum(tmp*x[k:(n-m+k)])
    }
    return(tmp[(n-m+1)])
  }
}

icpwfunc <- function(exb, a, ni){
  # calculate cipw in a culster
  # exb is vector of (exp(X_i1%*%beta),...,exp(X_ini%*%beta))
  # ni is sample size of cluster i
  # a is vector of (A_i1,..,A_ini)
  t <- sum(a)
  if(t==ni | ni==1 | t==0){
    return(rep(1,ni))
  }else{
    denom <- SumAllComb(x=exb, m=t, n=ni)
    numer <- sapply(c(1:ni), FUN=function(j){SumAllComb(x=exb[-j],m=(t-a[j]),n=(ni-1))})
    numer[a==1] <- numer[a==1]*exb[a==1]
    return(denom/numer)
  }
}
  
  p1_IPWfunc <- function(A, Y, prob){
    mean(Y*A/prob)
  }
  
  p0_IPWfunc <- function(A, Y, prob){
    mean(Y*(1-A)/(1-prob))
  }
  
  p1_SWfunc <- function(A, Y, prob.num, prob.den){
    sum(Y*A*prob.num/prob.den)/sum(prob.num/prob.den*A)
  }
  
  p0_SWfunc <- function(A, Y, prob.num, prob.den){
    sum(Y*(1-A)*(1-prob.num)/(1-prob.den))/sum((1-prob.num)/(1-prob.den)*(1-A))
  }
  
  OR_func <- function(p0, p1){p1/(1-p1)/(p0/(1-p0))}
  
  
  p1_DR_funct <- function(A, Y, prob, m1) {
    mean(A*Y/prob-(A-prob)/prob*m1)
  }
  
  p0_DR_funct <- function(A, Y, prob, m0) {
    mean((1-A)*Y/(1-prob)+(A-prob)/(1-prob)*m0)
  }
  

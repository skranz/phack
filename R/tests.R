#####################################################################
# Tests for Applications 1 and 2
# Paper: Detecting p-hacking
# Original Authors: G. Elliott, N. Kudrin, K. Wuthrich
# Some cosmetic adaptions: Sebastian Kranz
#####################################################################

warn_p_rounding = function(p, min_bunch=3) {
  if (min_bunch <= 0) return()
  p = p[p > 0 & p < 1]
  if (length(p) < min_bunch) return()

  count = sort(table(p), decreasing = TRUE)
  if (count[1] < min_bunch) return()
  num_bunch = sum(count > min_bunch)

  msg = paste0("
Your p-values seem to face an unsolved rounding problem.\nFor example, you have ", count[1], " p-values that all have the value ", names(count)[1],".\nA likely reason for such bunching are rounding errors.\n(The rounding errors might be present in the \ncorresponding t-values or reported coefficients and standard errors\n from which p-values were computed.).\n This test for p-hacking is not valid if p-values face rounding errors.\nPlease first apply an appropriate de-rounding procedure for your p-values.")

  warning(msg)
  return(invisible(msg))

}

#' Binomial test for p-hacking
#'
#' Elliot et al. (2022) show that absent rounding errors, p-hacking and
#' publication bias, the density of p-values across many different tests
#' should be decreasing in the p-value.
#'
#' This means if we split any interval [p_min, p_max] in the center, a
#' significantly higher proportion in the right half than the left half
#' suggests a violation of the no p-hacking and no publication bias assumption
#' (if there are no rounding errors or p-values are approbriately de-rounded).
#'
#' This function tests this via a Binomial test.
#'
#' @param p vector of p-values (should be derounded if there are rounding errors)
#' @param p_min lower bound of the interval used for the test
#' @param p_max upper bound of the interval used for the test
#' @param open_interval if TRUE take p-values from open interval (p_min, p_max).
#' @param min_bunch minimal number of elements of p that have exactly the same p-value in order to show a warning that there seems to be a rounding problem.
phack_test_binomial <-function(p,p_min=0.04, p_max=0.05, open_interval=FALSE, min_bunch=3){
  P = p
  # closed interval
  if (!open_interval){
    P = P[P<=p_max & P>=p_min]
  } else {
    P = P[P<p_max & P>p_min]
  }
  warn_p_rounding(P, min_bunch)

  nn = length(P)
  kk = sum(P>(p_max+p_min)/2)
  return(1 - pbinom(kk-1, nn, 0.5))
}


#' Fisher's test
#'
#' Similar to phack_binomial_test but using Fisher's test instead
#' of the binomial test.
#'
#' For this test the half-open interval [p_min, p_max) is used.
#'
#' @param p vector of p-values (should be derounded if there are rounding errors)
#' @param p_min lower bound of the interval used for the test
#' @param p_max upper bound of the interval used for the test
#' @param min_bunch minimal number of elements of p that have exactly the same p-value in order to show a warning that there seems to be a rounding problem.
phack_test_fisher <- function(p, p_min, p_max, min_bunch=3){
  P = p
  P = P[P<p_max & P>=p_min]
  warn_p_rounding(P, min_bunch)

  nn = length(P)
  statFM = -2*sum(log(1 - (P-p_min)/(p_max-p_min)))
  return(1 - pchisq(statFM, df = 2*nn))
}


#' Simulate Brownian Bridge (BB) and ||LCM(BB)-BB||
#'
#' @param M -- number of repetitions
SimBB = function(M){
  N = 10000
  x = 1:N
  x = x/N
  BBsup = matrix(0, nrow=M, ncol = 1)
  for (m in 1:M){
    eps <- rnorm(n = N, sd = 1, mean = 0)
    eps = eps/sqrt(N)
    W = cumsum(eps)
    B = W - x*W[N]
    C = c(0,x)
    B = c(0,B)
    lcmaj = gcmlcm(C,B, type="lcm")
    f_lcm = approxfun(lcmaj$x.knots, lcmaj$y.knots)
    y = f_lcm(C)
    BBsup[m,1] = max(abs(y - B))
  }
  return(BBsup)
}

get.phack.F_LCMsup = function() {
  res = getOption("phack.F_LCMsup")
  if (is.null(res)) {
    res = ecdf(SimBB(10000))
    options(phack.F_LCMsup = res)
  }
  res
}

#' LCM test on [p_min, p_max]
#' @param p -- vector of p-values
#' @param p_min lower bound of the interval used for the test
#' @param p_max upper bound of the interval used for the test
#' @param F_LCMsup cdf for LCM test
#' @param min_bunch minimal number of elements of p that have exactly the same p-value in order to show a warning that there seems to be a rounding problem.
phack_test_lcm <-function(p,p_min, p_max, F_LCMsup =get.phack.F_LCMsup(), min_bunch=3){
  P = p
  P = P[P<=p_max & P>=p_min]
  warn_p_rounding(P, min_bunch)

  nn = length(P)
  f<- ecdf(P)
  x = seq(0,1,length=1000)
  y = f(x*(p_max-p_min)+p_min)
  lcmaj = gcmlcm(x,y, type="lcm")
  ff_lcm = approxfun(lcmaj$x.knots, lcmaj$y.knots)
  z = as.numeric(ff_lcm(x))
  Test_stat = sqrt(nn)*max(abs(y-z))
  return(1 - F_LCMsup(Test_stat))
}

#' Discontinuity test
#'
#' @param p vector of p-values
#' @param c potential discontinuity point
#' @param min_bunch minimal number of elements of p that have exactly the same p-value in order to show a warning that there seems to be a rounding problem.

phack_test_discontinuity <- function(p, c, min_bunch=3){
  warn_p_rounding(p, min_bunch)
  res = rddensity(p, c = c)
  return(res$test$p_jk)
}

##### Cox&Shi test #####

# Critical value function for 2-sided t-test
cv2 <- function(p){
  cv = qnorm(1 - p/2)
  return(cv)
}

lambda2 <- function(x1, x2, h){
  lambda = pnorm(cv2(x1) - h) - pnorm(cv2(x2) - h)+pnorm(cv2(x1) + h) - pnorm(cv2(x2) + h)
  return(lambda)
}

#Bounds on the proportions (two-sided t-test);
#[p_min, p_max] interval with J bins
Bound0 <- function(p_min, p_max, J){
  h = seq(0,100,by = 0.001)
  X = linspace(p_min, p_max, J+1)
  B = matrix(0, J, 1)
  for (j in 1:(J)){
    Obj1 = lambda2(X[j], X[j+1], h)
    B[j] = max(Obj1)
  }
  if (p_min==0){
  B[1]=1
  }
  return(B)
}

#Bounds on the first differences of proportions (two-sided t-test);
#[p_min, p_max] interval with J bins
Bound1 <- function(p_min, p_max, J){
  h = seq(0,100, by = 0.001)
  X = linspace(p_min, p_max, J+1)
  B = matrix(0, J-1, 1)
  for (j in 1:(J-1)){
    Obj1 = lambda2(X[j], X[j+1], h)
    Obj2 = lambda2(X[j+1], X[j+2], h)
    A = Obj2-Obj1
    B[j] = max(abs(A))
  }
  if (p_min==0){
  B[1]=1
  }
  return(B)
}

#Bounds on the second differences of proportions (two-sided t-test);
#[p_min, p_max] interval with J bins
Bound2 <- function(p_min, p_max, J){
  h = seq(0,100, by = 0.001)
  X = linspace(p_min, p_max, J+1)
  B = matrix(0, J-2, 1)
  for (j in 1:(J-2)){
   Obj1 = lambda2(X[j], X[j+1], h)
   Obj2 = lambda2(X[j+1], X[j+2], h)
   Obj3 = lambda2(X[j+2], X[j+3], h)
   A = Obj3 - 2*Obj2 + Obj1
   B[j] = max(abs(A))
  }
  if (p_min==0){
  B[1]=1
  }
  return(B)
}



#' Cox-Shi histogram test and more general test for K-monotonicity and bounds on [p_min, p_max] interval
#'
#' For the defaults K=1 and use_bounds=FALSE we have a basic histogram test.
#'
#' @param p vector of p-values (make sure rounding problems are dealt with)
#' @param article vector of unique article ids for approbriate clustuering
#' @param J number of subintervals
#' @param K degree of K-monotonicity (see Section 4.3 in Elliot et al. 2022)
#' @param use_bounds use bounds or test without bounds (see Appendix A in Elliot et al. 2022)
#' @param min_bunch minimal number of elements of p that have exactly the same p-value in order to show a warning that there seems to be a rounding problem.
phack_test_cox_shi <- function(p, article=NA, p_min=0, p_max=0.15, J=30, K=1, use_bounds=FALSE, min_bunch=3){


  Q = p

  if (length(article)>1){
    article = article[Q<=p_max & Q>=p_min]
    uni_article = unique(article)
    ind = match(article, uni_article)
    indu = unique(ind)
  } else {
    ind = indu = NA
  }


  B = use_bounds*1L

  B0 = Bound0(p_min, p_max, J)
  B1 = Bound1(p_min, p_max, J)
  B2 = Bound2(p_min, p_max, J)

  P = Q[Q<=p_max & Q>=p_min]
  warn_p_rounding(P, min_bunch)


  Bnd_adj = length(P)/length(Q)
  N = length(P)
  bin = seq(p_min,p_max,length=J+1)
  Phat = matrix(0, nrow=J-1, ncol = 1)

  for (s in 1:(J-1)){
    Phat[s] = sum((P>bin[s])*(P<=bin[s+1]))/N
  }
  Phat[1] = Phat[1]+sum(P==bin[1])/N
  if (B==0){
    B0 = -matrix(1, nrow = J, ncol = 1)
  }
  if (B==1){
    B0=-B0/Bnd_adj
    B1=-B1/Bnd_adj
    B2=-B2/Bnd_adj
    if (p_min==0){
      B0[1]=-1
      B1[1]=-1
      B2[1]=-1
    }
  }

  if (length(ind)>1){
    Omega = matrix(0, J-1, J-1)
    for (i in c(indu)){
      X = P[ind==i]

      mq = repmat(matrix(X, 1, length(X)), J-1,1)
      a1 = (mq<= bin[2:J])
      a2 = (mq>bin[1:(J-1)])
      mq0 = 1*(mq==0)
      mq0[2:(J-1),] = 0
      mq = 1*(a1*a2) + mq0
      mq = mq - repmat((Phat), 1,length(X))

      Omega = Omega + (mq)%*%matrix(1, length(X), length(X))%*%t(mq)
    }
    Omega = Omega/length(P)
  }
  if (length(ind)==1){
    if (min(Phat)==0){
      Qhat = Phat*N/(N+1) + 1/(J*(N+1))
      Omega = diag(c(Qhat)) - Qhat%*%t(Qhat)
    }
    if (min(Phat)>0){
      Omega = diag(c(Phat)) - Phat%*%t(Phat)
    }
  }
  D = matrix(0, J-1, J)
  for (i in 1:(J-1)){
    for (j in 1:J){
      if (i==j){
        D[i,j] = -1
      }
      if (i+1==j){
        D[i,j] = 1
      }

    }
  }
  Dk = -D
  if (K>1){
    d = D
    for (k in 2:K){
      d = D[1:(J-k), 1:(J-k+1)]%*%d
      Dk = rbind(Dk, (-1)^k*d)
    }
  }
  if (B==0){
    Dk = rbind(-diag(J), diag(J), Dk)
  }
  if (B==1){
    Dk = rbind(-diag(J),-Dk, diag(J), Dk)
  }
  eJ = matrix(0, J, 1)
  eJ[J] = 1

  F1 = rbind(-diag(J-1), matrix(1, 1,J-1))
  c = matrix(0, (K+1)*(J-K/2), 1)

  if (B==0){
    c = rbind(matrix(-1, J,1), c)
  }
  if (B==1){

    if (K==0){
      c=rbind(B0, c)
    }
    if (K==1){
      c=rbind(B0, B1, c)
    }
    if (K==2){
      c=rbind(B0, B1, B2, c)
    }
  }
  A = Dk%*%F1
  b = Dk%*%eJ - c

  if (abs(det(Omega))>0){
  myQ = solve(Omega)
  fn <- function(t){
    Obj = N*(t(Phat - t))%*%myQ%*%((Phat - t))
    return(Obj)
  }
  t0 = matrix(1, (J-1), 1)/(J)
  res = fmincon(t0, fn, A = A, b = b)
  t_opt = t(t(res$par))
  #print(t_opt)
  T = fn(t_opt)
  Ba = A[which(res$info$lambda$ineqlin>0),]
  JX = qr(Ba)$rank
  if (res$convergence==0){
    return(as.numeric(
      1 - pchisq(T, df = JX)*(JX>0)
    ))
  }
  else {
    return(999)
  }
  }
  else {
    return (888)
  }
}


CoxShi <- function(Q, ind, p_min, p_max, J, K, B){
  restorepoint::restore.point("coxshi")
  B0 = Bound0(p_min, p_max, J)
  B1 = Bound1(p_min, p_max, J)
  B2 = Bound2(p_min, p_max, J)

  P = Q[Q<=p_max & Q>=p_min]
  if (length(ind)>1){
    ind = ind[Q<=p_max & Q>=p_min]
    indu = unique(ind)
  }
  Bnd_adj = length(P)/length(Q)
  N = length(P)
  bin = seq(p_min,p_max,length=J+1)
  Phat = matrix(0, nrow=J-1, ncol = 1)

  for (s in 1:(J-1)){
    Phat[s] = sum((P>bin[s])*(P<=bin[s+1]))/N
  }
  Phat[1] = Phat[1]+sum(P==bin[1])/N
  if (B==0){
    B0 = -matrix(1, nrow = J, ncol = 1)
  }
  if (B==1){
    B0=-B0/Bnd_adj
    B1=-B1/Bnd_adj
    B2=-B2/Bnd_adj
    if (p_min==0){
      B0[1]=-1
      B1[1]=-1
      B2[1]=-1
    }
  }

  if (length(ind)>1){
    Omega = matrix(0, J-1, J-1)
    for (i in c(indu)){
      X = P[ind==i]

      mq = repmat(matrix(X, 1, length(X)), J-1,1)
      a1 = (mq<= bin[2:J])
      a2 = (mq>bin[1:(J-1)])
      mq0 = 1*(mq==0)
      mq0[2:(J-1),] = 0
      mq = 1*(a1*a2) + mq0
      mq = mq - repmat((Phat), 1,length(X))

      Omega = Omega + (mq)%*%matrix(1, length(X), length(X))%*%t(mq)
    }
    Omega = Omega/length(P)
  }
  if (length(ind)==1){
    if (min(Phat)==0){
      Qhat = Phat*N/(N+1) + 1/(J*(N+1))
      Omega = diag(c(Qhat)) - Qhat%*%t(Qhat)
    }
    if (min(Phat)>0){
      Omega = diag(c(Phat)) - Phat%*%t(Phat)
    }
  }
  D = matrix(0, J-1, J)
  for (i in 1:(J-1)){
    for (j in 1:J){
      if (i==j){
        D[i,j] = -1
      }
      if (i+1==j){
        D[i,j] = 1
      }

    }
  }
  Dk = -D
  if (K>1){
    d = D
    for (k in 2:K){
      d = D[1:(J-k), 1:(J-k+1)]%*%d
      Dk = rbind(Dk, (-1)^k*d)
    }
  }
  if (B==0){
    Dk = rbind(-diag(J), diag(J), Dk)
  }
  if (B==1){
    Dk = rbind(-diag(J),-Dk, diag(J), Dk)
  }
  eJ = matrix(0, J, 1)
  eJ[J] = 1

  F1 = rbind(-diag(J-1), matrix(1, 1,J-1))
  c = matrix(0, (K+1)*(J-K/2), 1)

  if (B==0){
    c = rbind(matrix(-1, J,1), c)
  }
  if (B==1){

    if (K==0){
      c=rbind(B0, c)
    }
    if (K==1){
      c=rbind(B0, B1, c)
    }
    if (K==2){
      c=rbind(B0, B1, B2, c)
    }
  }
  A = Dk%*%F1
  b = Dk%*%eJ - c

  if (abs(det(Omega))>0){
  myQ = solve(Omega)
  fn <- function(t){
    Obj = N*(t(Phat - t))%*%myQ%*%((Phat - t))
    return(Obj)
  }
  t0 = matrix(1, (J-1), 1)/(J)
  res = fmincon(t0, fn, A = A, b = b)
  t_opt = t(t(res$par))
  print(t_opt)
  T = fn(t_opt)
  Ba = A[which(res$info$lambda$ineqlin>0),]
  JX = qr(Ba)$rank
  if (res$convergence==0){
    return(1 - pchisq(T, df = JX)*(JX>0))
  }
  else {
    return(999)
  }
  }
  else {
    return (888)
  }
}

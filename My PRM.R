

my_prm<- function (X, y, a, fairct = 4, opt = "l1m", usesvd = FALSE) 
{
  
  n = nrow(X)
  p = ncol(X)
  
  # Compute l1median and median
  mx <- l1median(X) ###
  my <- median(y) ###
  
  # Center and scale the data
  Xmc <- as.matrix(scale(X, center = mx, scale = FALSE))
  ymc <- as.vector(scale(y, center = my, scale = FALSE))
  
  #W_i_X;
  wx <- sqrt(apply(Xmc^2, 1, sum))
  wx <- wx/median(wx)
  wx <- 1/((1 + abs(wx/fairct))^2)
  
  #W_i_r;
  wy<- abs(ymc)
  wy <-  wy/median( wy)
  wy <- 1/((1 + abs( wy/fairct))^2)
  
  # Combine weights
  #W_i;  
  w <- wx * wy
  Xw <- Xmc * sqrt(w)
  yw <- ymc * sqrt(w)
  
  pls1_nipals<-function (Xw, yw, a, it = 50, tol = 1e-08) 
  {
    Xh <- Xw
    yh <- yw
    TT <- NULL
    P <- NULL
    C <- NULL
    W <- NULL
    for (h in 1:a) {
      wh <- t(Xh) %*% yh
      wh <- wh/as.vector(sqrt(t(wh) %*% wh))
      th <- Xh %*% wh
      ch <- as.numeric(crossprod(yh, th)/crossprod(th))
      ph <- t(Xh) %*% th/as.vector(t(th) %*% th)
      Xh <- Xh - th %*% t(ph)
      yh <- yh - th * ch
      TT <- cbind(TT, th)
      P <- cbind(P, ph)
      C <- c(C, ch)
      W <- cbind(W, wh)
    }
    b <- W %*% solve(t(P) %*% W) %*% C
    
    list(P = P, TT = TT, W = W, C = C, b = b)
  }
  
  loops <- 1
  ngamma <- 10^5
  difference <- 1
  while ((difference > 0.01) && loops < 30) {
    ngammaold <- ngamma
    pls <- pls1_nipals(Xw, yw, a)
    b <- pls$b
    
    gamma <- t(t(yw) %*% pls$TT)
    T <- pls$TT/sqrt(w)
    r <- ymc - T %*% gamma
    
    # W_i_r Update;
    rc <- r - median(r)
    r <- rc/median(abs(rc))            ### r <- ymc - (Xmc)%*% pls$b 
    wy <- 1/((1 + abs(r/fairct))^2)
    
    mt <- l1median(T) ###
    
    dt <- T - mt
    wt <- sqrt(sum(dt^2))
    wt <- wt/median(wt)
    wt <- 1/((1 + abs(wt/fairct))^2)
    
    ngamma <- sqrt(sum(pls$b^2))
    difference <- abs(ngamma - ngammaold)/ngamma
    w <- drop(wy * wt)
    w0 <- which(w == 0)
    if (length(w0) != 0) {
      w <- replace(w, list = w0, values = 10^(-6))
    }
    
    Xw <- Xmc * sqrt(w)
    yw <- ymc * sqrt(w)
    
    loops <- loops + 1
  }
  
  list(coef = b, wy = wy, wt = wt, w = w, scores = T, 
       mx = mx, my = my, yscores = pls$C, xweights=pls$W ) ###
}


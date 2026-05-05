library(numDeriv)

# data generation function

MCCR.data.new2 = function(n,c,alp1,lam1,alp2,lam2,b){
  
  y = rep(NA,n)
  d = rep(NA,n)
  
  U1 = runif(n,0,1)
  R = runif(n,0,b)
  
  T_cure = rweibull(n,shape=alp2, scale = lam2^(-1/alp2))
  
  T1_unc = rweibull(n,shape=alp1, scale = lam1^(-1/alp1))
  T2_unc = rweibull(n,shape=alp2, scale = lam2^(-1/alp2))
  T_unc = pmin(T1_unc, T2_unc)
  d = apply(cbind(T1_unc,T2_unc),1,which.min) 
  
  
  T = ifelse(U1<=c,T_cure, T_unc)
  d[U1<=c]  = 2   # cure case
  
  y = pmin(T,R) # censoring
  d[T>R] = 0    # Cesored case
  return(data.frame(cbind(y,d)))
}

# calculation of observed log-likelihod function

log.lik.fn = function(param=c(cure,alp1,lam1,alp2,lam2),y0,y1,y2){
  
  h1.0 = param[2]*param[3]*(y0^(param[2]-1))
  h1.1 = param[2]*param[3]*(y1^(param[2]-1)) 
  h1.2 = param[2]*param[3]*(y2^(param[2]-1))
  
  h2.0 = param[4]*param[5]*(y0^(param[4]-1))
  h2.1 = param[4]*param[5]*(y1^(param[4]-1)) 
  h2.2 = param[4]*param[5]*(y2^(param[4]-1))
  
  H1.0 = param[3]*(y0^param[2])
  H1.1 = param[3]*(y1^param[2])
  H1.2 = param[3]*(y2^param[2])
  
  H2.0 = param[5]*(y0^param[4])
  H2.1 = param[5]*(y1^param[4])
  H2.2 = param[5]*(y2^param[4])
  
  l.obs = sum(log(1-param[1])+log(h1.1)-H1.1-H2.1)+sum(log((param[1]*exp(-H2.0))+((1-param[1])*exp(-H1.0-H2.0))))+
          sum(log(h2.2))+sum(log((param[1]*exp(-H2.2))+((1-param[1])*exp(-H1.2-H2.2))))
  return(l.obs)
}


# EM algorithm with 2 competing risks

EM.MCCR2=function(data,tol,maxit,c,alpha1,lambda1,alpha2,lambda2){
  
  data0 = data[data$d==0,]
  data1 = data[data$d==1,]
  data2 = data[data$d==2,]
  
  y0 = data0$y
  y1 = data1$y
  y2 = data2$y
  
  p.new=matrix(0,ncol=1,nrow=5)
  p.old=matrix(0,ncol=1,nrow=5)
  
  p.old[1,1]=c
  p.old[2,1]=alpha1
  p.old[3,1]=lambda1
  p.old[4,1]=alpha2
  p.old[5,1]=lambda2
  
  continue = TRUE
  iter=1
  
  while(continue){
    
    #E step:
    
    h1.0 = p.old[2,1]*p.old[3,1]*(y0^(p.old[2,1]-1))
    h1.1 = p.old[2,1]*p.old[3,1]*(y1^(p.old[2,1]-1)) 
    h1.2 = p.old[2,1]*p.old[3,1]*(y2^(p.old[2,1]-1))
    
    h2.0 = p.old[4,1]*p.old[5,1]*(y0^(p.old[4,1]-1))
    h2.1 = p.old[4,1]*p.old[5,1]*(y1^(p.old[4,1]-1)) 
    h2.2 = p.old[4,1]*p.old[5,1]*(y2^(p.old[4,1]-1))
    
    H1.0 = p.old[3,1]*(y0^p.old[2,1])
    H1.1 = p.old[3,1]*(y1^p.old[2,1])
    H1.2 = p.old[3,1]*(y2^p.old[2,1])
    
    H2.0 = p.old[5,1]*(y0^p.old[4,1])
    H2.1 = p.old[5,1]*(y1^p.old[4,1])
    H2.2 = p.old[5,1]*(y2^p.old[4,1])
    
    w0 = ((1-p.old[1,1])*exp(-H1.0-H2.0))/((p.old[1,1]*exp(-H2.0)) + ((1-p.old[1,1])*exp(-H1.0-H2.0))) 
    w2 = ((1-p.old[1,1])*exp(-H1.2-H2.2))/((p.old[1,1]*exp(-H2.2)) + ((1-p.old[1,1])*exp(-H1.2-H2.2))) 
    
    # defining the function involving the cure rate 
    
    Q = function(par=c(c1)){
      
      #c1 = exp(par[1])/(1+exp(par[1]))
      #res =  (length(y1)*log(1-c1)) + sum(((1-w0)*log(c1))+(w0*log(1-c1))) + sum(((1-w2)*log(c1))+(w2*log(1-c1))) 
      res =  (length(y1)*log(1-par[1])) + sum(((1-w0)*log(par[1]))+(w0*log(1-par[1]))) + sum(((1-w2)*log(par[1]))+(w2*log(1-par[1]))) 
      return(res)
      
    }# end of Q
    
    
    # defining the function involving parameters of risk 1
    
    Q1 = function(par1=c(a1,l1)){
      
      haz1.1 = par1[1]*par1[2]*(y1^(par1[1]-1))
      chaz1.0 = par1[2]*(y0^par1[1])
      chaz1.1 = par1[2]*(y1^par1[1])
      chaz1.2 = par1[2]*(y2^par1[1])
      
      res1 = sum(log(haz1.1) - chaz1.1) - sum(w0*chaz1.0) - sum(w2*chaz1.2)
      #res1 = sum(log(par1[1])+log(par1[2])+((par1[1]-1)*log(y1)) - chaz1.1) - sum(w0*chaz1.0) - sum(w2*chaz1.2)
      return(-res1)
      
    }# end of Q1
    
    # defining the function involving parameters of risk 2 
    
    Q2 = function(par2=c(a2,l2)){
      
      haz2.2 = par2[1]*par2[2]*(y2^(par2[1]-1))
      chaz2.0 = par2[2]*(y0^par2[1])
      chaz2.1 = par2[2]*(y1^par2[1])
      chaz2.2 = par2[2]*(y2^par2[1])
      
      res2 = sum(log(haz2.2)) - sum(chaz2.1) - sum(chaz2.0) - sum(chaz2.2)
      #res2 = sum(log(par2[1])+log(par2[2])+((par2[1]-1)*log(y2))) - sum(chaz2.1) - sum(chaz2.0) - sum(chaz2.2)
      return(-res2)
      
    }# end of Q2
    
    c.new = optimize(Q,c(p.old[1,1]-0.01,p.old[1,1]+0.01),maximum=T)$maximum
    risk1.new = optim(par=c(p.old[2,1],p.old[3,1]),fn=Q1, method="Nelder-Mead")$par
    risk2.new = optim(par=c(p.old[4,1],p.old[5,1]),fn=Q2, method="Nelder-Mead")$par
    
    p.new = matrix(c(c.new,risk1.new,risk2.new)) # updated estimates
    
    #print(p.new)
    
    iter = iter + 1
    
    continue=(abs((p.new[1,1]-p.old[1,1])/p.old[1,1])>tol | abs((p.new[2,1]-p.old[2,1])/p.old[2,1])>tol |
                abs((p.new[3,1]-p.old[3,1])/p.old[3,1])>tol | abs((p.new[4,1]-p.old[4,1])/p.old[4,1])>tol | 
                abs((p.new[5,1]-p.old[5,1])/p.old[5,1])>tol) & (iter<maxit)
    
    p.old[1,1]=p.new[1,1]
    p.old[2,1]=p.new[2,1]
    p.old[3,1]=p.new[3,1]
    p.old[4,1]=p.new[4,1]
    p.old[5,1]=p.new[5,1]
    
    
  }#end of while
  
  #p.est = ifelse(iter<maxit,p.new,rep(0,5))
  
  #print("iteration")
  #print(iter)
  
  #return(p.est)
  
  hess.mat = hessian(log.lik.fn,p.new,y0=y0,y1=y1,y2=y2,method="Richardson")
  FI = solve(-1*hess.mat)
  se.c = sqrt(FI[1,1])
  se.alp1 = sqrt(FI[2,2])
  se.lam1 = sqrt(FI[3,3])
  se.alp2 = sqrt(FI[4,4])
  se.lam2 = sqrt(FI[5,5])
  
  se.est = c(se.c,se.alp1,se.lam1,se.alp2,se.lam2)
  
  EM.out = rep(NA,10)
  
  if(iter<maxit & sum(is.na(se.est))==0){
    EM.out = c(p.new,se.est)
  }else{
    EM.out = rep(0,10)
  }
  
  return(EM.out)
  
}#end of the EM function


# calling the functions

  
  N = 500
  n = 600
  c.true = 0.6
  alp1.true = 1.5
  lam1.true = 2
  alp2.true = 0.5
  lam2.true = 1
  b = 2
  incr = 0.2
  
  c.hat = rep(NA,N)
  alp1.hat = rep(NA,N)
  lam1.hat = rep(NA,N)
  alp2.hat = rep(NA,N)
  lam2.hat = rep(NA,N)
  
  std.c = rep(NA,N)
  std.alp1 = rep(NA,N)
  std.lam1 = rep(NA,N)
  std.alp2 = rep(NA,N)
  std.lam2 = rep(NA,N)
  
  c.temp = rep(NA,N)
  alp1.temp = rep(NA,N)
  lam1.temp = rep(NA,N)
  alp2.temp = rep(NA,N)
  lam2.temp = rep(NA,N)
  
  lcl_c_95=rep(NA,N)
  lcl_alp1_95=rep(NA,N)
  lcl_lam1_95=rep(NA,N)
  lcl_alp2_95=rep(NA,N)
  lcl_lam2_95=rep(NA,N)
  ucl_c_95=rep(NA,N)
  ucl_alp1_95=rep(NA,N)
  ucl_lam1_95=rep(NA,N)
  ucl_alp2_95=rep(NA,N)
  ucl_lam2_95=rep(NA,N)
  
  t_c_95=0
  t_alp1_95=0
  t_lam1_95=0
  t_alp2_95=0
  t_lam2_95=0
  
  j = 1
  count = 0
  
  while(j <= N){
    print(j)
    data = MCCR.data.new2(n=n,c=c.true,alp1=alp1.true,lam1=lam1.true,alp2=alp2.true,lam2=lam2.true,b=b)
    
    c.init = sample(seq((c.true-(incr*abs(c.true))),(c.true+(incr*abs(c.true))),by=0.01),1)
    alp1.init = sample(seq((alp1.true-(incr*abs(alp1.true))),(alp1.true+(incr*abs(alp1.true))),by=0.01),1)
    lam1.init = sample(seq((lam1.true-(incr*abs(lam1.true))),(lam1.true+(incr*abs(lam1.true))),by=0.01),1)
    alp2.init = sample(seq((alp2.true-(incr*abs(alp2.true))),(alp2.true+(incr*abs(alp2.true))),by=0.01),1)
    lam2.init = sample(seq((lam2.true-(incr*abs(lam2.true))),(lam2.true+(incr*abs(lam2.true))),by=0.01),1)
    
    EM = EM.MCCR2(data=data,tol=0.001,maxit=500,c=c.init,alpha1=alp1.init,lambda1=lam1.init,alpha2=alp2.init,lambda2=lam2.init)
    
    if(all(EM==0)){
      j = j + 0
      count=count+1
    }else{
    
    c.hat[j] = EM[1]
    alp1.hat[j] = EM[2]
    lam1.hat[j] = EM[3]
    alp2.hat[j] = EM[4]
    lam2.hat[j] = EM[5]
    std.c[j] = EM[6]
    std.alp1[j] = EM[7]
    std.lam1[j] = EM[8]
    std.alp2[j] = EM[9]
    std.lam2[j] = EM[10]
    
    c.temp[j] = c.hat[j] - c.true
    alp1.temp[j] = alp1.hat[j] - alp1.true
    lam1.temp[j] = lam1.hat[j] - lam1.true
    alp2.temp[j] = alp2.hat[j] - alp2.true
    lam2.temp[j] = lam2.hat[j] - lam2.true
    
    lcl_c_95[j] = c.hat[j] - (1.96*std.c[j])
    lcl_alp1_95[j] = alp1.hat[j] - (1.96*std.alp1[j])
    lcl_lam1_95[j] = lam1.hat[j] - (1.96*std.lam1[j])
    lcl_alp2_95[j] = alp2.hat[j] - (1.96*std.alp2[j])
    lcl_lam2_95[j] = lam2.hat[j] - (1.96*std.lam2[j])
    
    ucl_c_95[j] = c.hat[j] + (1.96*std.c[j])
    ucl_alp1_95[j] = alp1.hat[j] + (1.96*std.alp1[j])
    ucl_lam1_95[j] = lam1.hat[j] + (1.96*std.lam1[j])
    ucl_alp2_95[j] = alp2.hat[j] + (1.96*std.alp2[j])
    ucl_lam2_95[j] = lam2.hat[j] + (1.96*std.lam2[j])
    
    if(c.true > lcl_c_95[j] & c.true < ucl_c_95[j]){
      t_c_95 = t_c_95 + 1
    }
    if(alp1.true > lcl_alp1_95[j] & alp1.true < ucl_alp1_95[j]){
      t_alp1_95 = t_alp1_95 + 1
    }
    if(lam1.true > lcl_lam1_95[j] & lam1.true < ucl_lam1_95[j]){
      t_lam1_95 = t_lam1_95 + 1
    }
    if(alp2.true > lcl_alp2_95[j] & alp2.true < ucl_alp2_95[j]){
      t_alp2_95 = t_alp2_95 + 1
    }
    if(lam2.true > lcl_lam2_95[j] & lam2.true < ucl_lam2_95[j]){
      t_lam2_95 = t_lam2_95 + 1
    }
    
  j = j+1
    
   }#end of else
  }#end of while
  
  avg.c = sum(c.hat)/N
  avg.alp1 = sum(alp1.hat)/N
  avg.lam1 = sum(lam1.hat)/N
  avg.alp2 = sum(alp2.hat)/N
  avg.lam2 = sum(lam2.hat)/N
  
  avg.std.c = sum(std.c)/N
  avg.std.alp1 = sum(std.alp1)/N
  avg.std.lam1 = sum(std.lam1)/N
  avg.std.alp2 = sum(std.alp2)/N
  avg.std.lam2 = sum(std.lam2)/N
  
  bias.c = sum(c.temp)/N
  bias.alp1 = sum(alp1.temp)/N
  bias.lam1 = sum(lam1.temp)/N
  bias.alp2 = sum(alp2.temp)/N
  bias.lam2 = sum(lam2.temp)/N
  
  rmse.c = sqrt(sum(c.temp^2)/(N-1))
  rmse.alp1 = sqrt(sum(alp1.temp^2)/(N-1))
  rmse.lam1 = sqrt(sum(lam1.temp^2)/(N-1))
  rmse.alp2 = sqrt(sum(alp2.temp^2)/(N-1))
  rmse.lam2 = sqrt(sum(lam2.temp^2)/(N-1))
  
  cp_c_95 = t_c_95/N
  cp_alp1_95 = t_alp1_95/N
  cp_lam1_95 = t_lam1_95/N
  cp_alp2_95 = t_alp2_95/N
  cp_lam2_95 = t_lam2_95/N
  
  
  avg.c
  avg.alp1
  avg.lam1
  avg.alp2
  avg.lam2
  
  bias.c
  bias.alp1
  bias.lam1
  bias.alp2
  bias.lam2
  
  rmse.c
  rmse.alp1
  rmse.lam1
  rmse.alp2
  rmse.lam2
  
  avg.std.c
  avg.std.alp1
  avg.std.lam1
  avg.std.alp2
  avg.std.lam2
  
  cp_c_95 
  cp_alp1_95 
  cp_lam1_95 
  cp_alp2_95 
  cp_lam2_95 
  
  par.est = round(c(avg.c,avg.alp1,avg.lam1,avg.alp2,avg.lam2),4)
  par.bias = round(c(bias.c,bias.alp1,bias.lam1,bias.alp2,bias.lam2),4)
  par.rmse = round(c(rmse.c,rmse.alp1,rmse.lam1,rmse.alp2,rmse.lam2),4)
  par.se = round(c(avg.std.c,avg.std.alp1,avg.std.lam1,avg.std.alp2,avg.std.lam2),4)
  par.cp.95 = round(c(cp_c_95,cp_alp1_95,cp_lam1_95,cp_alp2_95,cp_lam2_95),4)
  
  data.frame(par.est,par.bias,par.rmse,par.se,par.cp.95)
  


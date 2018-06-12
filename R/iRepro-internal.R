.beta.i <-
function(z, params, sd.t, sd.tb, lower.limit, upper.limit){
  return(params[1]*exp(-(z/sd.tb)^2)*(pnorm((upper.limit*sd.tb - params[1]^2*z/sd.tb)/(params[1]*sd.t)) - pnorm((lower.limit*sd.tb - params[1]^2*z/sd.tb)/(params[1]*sd.t)))/sd.tb  )
}
.cfun1 <-
function(theta,ratings.1,ratings.2,i.limits){  
  
  n_classes = nrow(i.limits);
  n = length(ratings.1);
  f = n*log(theta[1]);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  
  f.table <- mat.or.vec(n_classes, n_classes)
  for (lb1 in 1:n_classes){
   for (lb2 in lb1:n_classes){
      i1.lb <- i.limits[lb1,1]
      i1.ub <- i.limits[lb1,2]
      i2.lb <- i.limits[lb2,1]
      i2.ub <- i.limits[lb2,2]      
                
      intg <- .f.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .f.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .f.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .f.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
      
      if(intg < 0){
        stop('f.table < 0')}else{
                f.table[lb1, lb2] <- - log(intg)
        }
      f.table[lb2, lb1] <- f.table[lb1, lb2]
    }
  } 
    
  nonzero <- table(ratings.1,ratings.2)>0
  f = f + sum(table(ratings.1,ratings.2)[nonzero]*f.table[nonzero]);   
  return(f);
}
.cfun2 <-
function(theta,ratings.1,ratings.2,i.limits){  
  
  n_classes = nrow(i.limits);
  n = length(ratings.1);
  f = n*log(theta[1]);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  
  f.table <- mat.or.vec(n_classes, n_classes)
  for (lb1 in 1:n_classes){
   for (lb2 in lb1:n_classes){
      i1.lb <- i.limits[lb1,1]
      i1.ub <- i.limits[lb1,2]
      i2.lb <- i.limits[lb2,1]
      i2.ub <- i.limits[lb2,2]      
                
      if((i1.ub - i2.ub)<(i1.lb - i2.lb)){
          intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i1.lb,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
        }else if((i1.ub - i2.ub)>(i1.lb - i2.lb)){
          intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i2.ub,1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
        }else{
          intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
        }
      
      if(intg < 0){
        stop('f.table < 0')}else{
                f.table[lb1, lb2] <- - log(intg)
        }
      f.table[lb2, lb1] <- f.table[lb1, lb2]
    }
  } 
    
  nonzero <- table(ratings.1,ratings.2)>0
  f = f + sum(table(ratings.1,ratings.2)[nonzero]*f.table[nonzero]);   
  return(f);
}
.cfun3 <-
function(theta,ratings.1,ratings.2){  
    
  n = nrow(ratings.1);
  f = n*log(theta[1]);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  
  for (i.n in 1:n){
    
    i1.lb <- ratings.1[i.n,1]
    i1.ub <- ratings.1[i.n,2]
    i2.lb <- ratings.2[i.n,1]
    i2.ub <- ratings.2[i.n,2]     
                
    intg <- .f.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .f.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .f.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .f.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
    
    if(intg < 0){
      stop('f.table < 0')}else{
        f <- f - log(intg)
      }
  }    
  return(f);
}
.cfun4 <-
function(theta,ratings.1,ratings.2){  
    
  n = nrow(ratings.1);
  f = n*log(theta[1]);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  
  for (i.n in 1:n){
    
    i1.lb <- ratings.1[i.n,1]
    i1.ub <- ratings.1[i.n,2]
    i2.lb <- ratings.2[i.n,1]
    i2.ub <- ratings.2[i.n,2]     
                
    if((i1.ub - i2.ub)<(i1.lb - i2.lb)){
      intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i1.lb,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
      }else if((i1.ub - i2.ub)>(i1.lb - i2.lb)){
        intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i2.ub,1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
        }else{
          intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
          }
    
    if(intg < 0){
      stop('f.table < 0')}else{
        f <- f - log(intg)
      }
  }    
  return(f);
}
.chi2.int.test <-
function(d.lower,d.upper,bins = 10,mu,stdev){
  
  # test if interval-censored data follows normal distribution with mean mu and standard deviation stdev
  # using chi-squared test
  #
  # INPUT:
  # d.lower - lower limits of data censoring intervals
  # d.upper - upper limits of data censoring intervals
  # bins - number of categories in chi-squared test
  # mu - mean of baseline distribution
  # stdev - standard deviation of baseline distribution
  
  breaks <- qnorm(seq(0,1,1/bins),mean=mu,sd=stdev)
  nresp <- length(d.lower) 
  lims <- sort(union(d.lower,d.upper))
  dis_lim <- sort(union(lims, breaks))
  
  n_counts <- length(dis_lim)-1
  dis_prob <- rep(0, n_counts)
  mids <- rep(0, n_counts)
  for(i in 1:n_counts){
    dis_prob[i] <- pnorm(dis_lim[i+1], mean=mu, sd=stdev) - pnorm(dis_lim[i], mean=mu, sd=stdev)  
    mids[i] <- 0.5*(dis_lim[i] + dis_lim[i+1])
  }

  counts <- rep(0,n_counts)
  vals <- rep(0,nresp)

  for(i in 1:nresp){
    ind1 <- which(dis_lim==d.lower[i])
    ind2 <- which(dis_lim==d.upper[i])
    normf <- sum(dis_prob[ind1:(ind2-1)])
    for(j in ind1:(ind2-1)){
      counts[j] <- counts[j] + dis_prob[j]/normf
    }  
  }
  
  n_counts.eq <- length(breaks)-1
  counts.eq <- rep(0, n_counts.eq)
  dis_prob.eq <- rep(0, n_counts.eq)
  mids.eq <- rep(0, n_counts.eq)
  
  for(i in 1:n_counts.eq){
    ind1 <- which(dis_lim == breaks[i])
    ind2 <- which(dis_lim == breaks[i+1])
    mids.eq[i] <- 0.5*(breaks[i] + breaks[i+1])
    counts.eq[i] <- sum(counts[ind1:(ind2-1)])
    dis_prob.eq[i] <- sum(dis_prob[ind1:(ind2-1)])
  }  
  
  t <- chisq.test(c(counts.eq),p = c(dis_prob.eq))  
  return(t);
}
.f.int <-
function(z, sgn, params, sd.t, lower.limit, upper.limit){  
  I <- function(u){
    return(exp(-(u/params[1])^2)*pnorm((sgn*u + z - params[3])/sd.t));
  }
  return(integrate(I,lower.limit,upper.limit,stop.on.error=1)$value);
}
.fmod.int <-
function(z1, sgn1, z2, sgn2, params, sd.t, lower.limit, upper.limit){   
  I <- function(u){
    return(exp(-(u/params[1])^2)*(pnorm((sgn1*u + z1 - params[3])/sd.t) - pnorm((sgn2*u + z2 - params[3])/sd.t)));
  } 
  return(integrate(I,lower.limit,upper.limit,stop.on.error=1)$value);
}
.gradcfun1 <-
function(theta,ratings.1,ratings.2,i.limits){
  
  n_classes = nrow(i.limits);
  n = length(ratings.1);
  Jf = c(n/theta[1], 0, 0);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  sdev.thb = sqrt(2*(theta[1]^2 + theta[2]^2));  
 
  Jf1.table <- mat.or.vec(n_classes, n_classes)
  Jf2.table <- mat.or.vec(n_classes, n_classes)
  Jf3.table <- mat.or.vec(n_classes, n_classes)
  
  for (lb1 in 1:n_classes){
   for (lb2 in lb1:n_classes){
      i1.lb <- i.limits[lb1,1]
      i1.ub <- i.limits[lb1,2]
      i2.lb <- i.limits[lb2,1]
      i2.ub <- i.limits[lb2,2]      
      
      corr.n <- (.f.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .f.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .f.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .f.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)));

      Jf1.table[lb1, lb2] <- - (.sigma.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigma.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigma.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigma.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf1.table[lb2, lb1] <- Jf1.table[lb1, lb2]
      Jf2.table[lb1, lb2] <- - (.sigmab.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigmab.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigmab.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigmab.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf2.table[lb2, lb1] <- Jf2.table[lb1, lb2]
      Jf3.table[lb1, lb2] <- (.beta.i(theta[3]-i2.ub,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .beta.i(i1.ub-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .beta.i(i1.lb-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .beta.i(theta[3]-i2.lb,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf3.table[lb2, lb1] <- Jf3.table[lb1, lb2]
    }
  }
  nonzero <- table(ratings.1,ratings.2)>0
  Jf[1] = Jf[1] + sum(table(ratings.1,ratings.2)[nonzero]*Jf1.table[nonzero]);
  Jf[2] = Jf[2] + sum(table(ratings.1,ratings.2)[nonzero]*Jf2.table[nonzero]);
  Jf[3] = Jf[3] + sum(table(ratings.1,ratings.2)[nonzero]*Jf3.table[nonzero]);
  
  return(Jf);
}
.gradcfun2 <-
function(theta,ratings.1,ratings.2,i.limits){
  
  n_classes = nrow(i.limits);
  n = length(ratings.1);
  Jf = c(n/theta[1], 0, 0);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  sdev.thb = sqrt(2*(theta[1]^2 + theta[2]^2));  
 
  Jf1.table <- mat.or.vec(n_classes, n_classes)
  Jf2.table <- mat.or.vec(n_classes, n_classes)
  Jf3.table <- mat.or.vec(n_classes, n_classes)
  
  for (lb1 in 1:n_classes){
   for (lb2 in lb1:n_classes){
      i1.lb <- i.limits[lb1,1]
      i1.ub <- i.limits[lb1,2]
      i2.lb <- i.limits[lb2,1]
      i2.ub <- i.limits[lb2,2]      
      
      if((i1.ub - i2.ub)<(i1.lb - i2.lb)){
        corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i1.lb,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
        }else if((i1.ub - i2.ub)>(i1.lb - i2.lb)){
          corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i2.ub,1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
          }else{
            corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
          }
          

      Jf1.table[lb1, lb2] <- - (.sigma.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigma.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigma.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigma.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf1.table[lb2, lb1] <- Jf1.table[lb1, lb2]
      Jf2.table[lb1, lb2] <- - (.sigmab.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigmab.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigmab.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigmab.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf2.table[lb2, lb1] <- Jf2.table[lb1, lb2]
      Jf3.table[lb1, lb2] <- (.beta.i(theta[3]-i2.ub,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .beta.i(i1.ub-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .beta.i(i1.lb-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .beta.i(theta[3]-i2.lb,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf3.table[lb2, lb1] <- Jf3.table[lb1, lb2]
    }
  }
  nonzero <- table(ratings.1,ratings.2)>0
  Jf[1] = Jf[1] + sum(table(ratings.1,ratings.2)[nonzero]*Jf1.table[nonzero]);
  Jf[2] = Jf[2] + sum(table(ratings.1,ratings.2)[nonzero]*Jf2.table[nonzero]);
  Jf[3] = Jf[3] + sum(table(ratings.1,ratings.2)[nonzero]*Jf3.table[nonzero]);  
  
  return(Jf);
}
.gradcfun3 <-
function(theta,ratings.1,ratings.2){
  
  n = nrow(ratings.1);
  Jf = c(n/theta[1], 0, 0);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  sdev.thb = sqrt(2*(theta[1]^2 + theta[2]^2));  
    
  for (i.n in 1:n){
      i1.lb <- ratings.1[i.n,1]
      i1.ub <- ratings.1[i.n,2]
      i2.lb <- ratings.2[i.n,1]
      i2.ub <- ratings.2[i.n,2]
      
      corr.n <- (.f.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .f.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .f.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .f.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)));

      Jf[1] <- Jf[1] - (.sigma.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigma.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigma.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigma.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf[2] <- Jf[2] - (.sigmab.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigmab.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigmab.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigmab.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf[3] <- Jf[3] + (.beta.i(theta[3]-i2.ub,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .beta.i(i1.ub-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .beta.i(i1.lb-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .beta.i(theta[3]-i2.lb,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
       
  }  
  return(Jf);
}
.gradcfun4 <-
function(theta,ratings.1,ratings.2){
  
  n = nrow(ratings.1);
  Jf = c(n/theta[1], 0, 0);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  sdev.thb = sqrt(2*(theta[1]^2 + theta[2]^2));  
    
  for (i.n in 1:n){
      i1.lb <- ratings.1[i.n,1]
      i1.ub <- ratings.1[i.n,2]
      i2.lb <- ratings.2[i.n,1]
      i2.ub <- ratings.2[i.n,2]
      
      if((i1.ub - i2.ub)<(i1.lb - i2.lb)){
        corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i1.lb,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
        }else if((i1.ub - i2.ub)>(i1.lb - i2.lb)){
          corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i2.ub,1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
          }else{
            corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
          }
          
      Jf[1] <- Jf[1] - (.sigma.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigma.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigma.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigma.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf[2] <- Jf[2] - (.sigmab.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigmab.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigmab.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigmab.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf[3] <- Jf[3] + (.beta.i(theta[3]-i2.ub,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .beta.i(i1.ub-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .beta.i(i1.lb-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .beta.i(theta[3]-i2.lb,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;       
  }  
  return(Jf);
}
.intervalICC.est1 <-
function(ratings, classes, c.limits, theta0){

  costFunc <- function(theta){
    return(.cfun1(theta,ratings$t1,ratings$t2,c.limits))
  }
  grad.costFunc <- function(theta){
    return(.gradcfun1(theta,ratings$t1,ratings$t2,c.limits))
  }
  
  th_optim <- tryCatch({constrOptim(theta=theta0, f=costFunc, grad=grad.costFunc, ui=rbind(c(1,0,0),c(0,1,0)), ci=c(.Machine$double.eps,.Machine$double.eps), method = "BFGS")}, error = function(e) {print("Numerical optimization failed; try optim.method=2"); return(0)})
    
  sigma2.b <- th_optim$par[2]^2
  sigma2.w <- th_optim$par[1]^2
  mu <- th_optim$par[3]
  icc <- sigma2.b/(sigma2.b + sigma2.w)
  loglikelihood <- -th_optim$value - nrow(ratings)*log(pi)/2
  
  list(icc = icc,
       sigma2.b = sigma2.b,
       sigma2.w = sigma2.w,
       mu = mu,
       loglikelihood = loglikelihood)
}
.intervalICC.est2 <-
function(ratings, classes, c.limits, theta0){

  costFunc <- function(theta){
    return(.cfun2(theta,ratings$t1,ratings$t2,c.limits))
  }
  grad.costFunc <- function(theta){
    return(.gradcfun2(theta,ratings$t1,ratings$t2,c.limits))
  }
  
  th_optim <- tryCatch({constrOptim(theta=theta0, f=costFunc, grad=grad.costFunc, ui=rbind(c(1,0,0),c(0,1,0)), ci=c(.Machine$double.eps,.Machine$double.eps), method = "BFGS")}, error = function(e) {print("Numerical optimization failed; try optim.method=1"); return(0)})
    
  sigma2.b <- th_optim$par[2]^2
  sigma2.w <- th_optim$par[1]^2
  mu <- th_optim$par[3]
  icc <- sigma2.b/(sigma2.b + sigma2.w)
  loglikelihood <- -th_optim$value - nrow(ratings)*log(pi)/2
  
  list(icc = icc,
       sigma2.b = sigma2.b,
       sigma2.w = sigma2.w,
       mu = mu,
       loglikelihood = loglikelihood)
}
.intervalICC.est3 <-
function(ratings1, ratings2, theta0){

  costFunc <- function(theta){
    return(.cfun3(theta,ratings1,ratings2))
  }
  grad.costFunc <- function(theta){
    return(.gradcfun3(theta,ratings1,ratings2))
  }
  
  th_optim <- tryCatch({constrOptim(theta=theta0, f=costFunc, grad=grad.costFunc, ui=rbind(c(1,0,0),c(0,1,0)), ci=c(.Machine$double.eps,.Machine$double.eps), method = "BFGS")}, error = function(e) {print("Numerical optimization failed; try optim.method=1"); return(0)})
    
  sigma2.b <- th_optim$par[2]^2
  sigma2.w <- th_optim$par[1]^2
  mu <- th_optim$par[3]
  icc <- sigma2.b/(sigma2.b + sigma2.w)
  loglikelihood <- -th_optim$value - nrow(ratings1)*log(pi)/2
  
  list(icc = icc,
       sigma2.b = sigma2.b,
       sigma2.w = sigma2.w,
       mu = mu,
       loglikelihood = loglikelihood)
}
.intervalICC.est4 <-
function(ratings1, ratings2, theta0){

  costFunc <- function(theta){
    return(.cfun4(theta,ratings1,ratings2))
  }
  grad.costFunc <- function(theta){
    return(.gradcfun4(theta,ratings1,ratings2))
  }
  
  th_optim <- tryCatch({constrOptim(theta=theta0, f=costFunc, grad=grad.costFunc, ui=rbind(c(1,0,0),c(0,1,0)), ci=c(.Machine$double.eps,.Machine$double.eps), method = "BFGS")}, error = function(e) {print("Numerical optimization failed; try optim.method=1"); return(0)})
    
  sigma2.b <- th_optim$par[2]^2
  sigma2.w <- th_optim$par[1]^2
  mu <- th_optim$par[3]
  icc <- sigma2.b/(sigma2.b + sigma2.w)
  loglikelihood <- -th_optim$value - nrow(ratings1)*log(pi)/2
  
  list(icc = icc,
       sigma2.b = sigma2.b,
       sigma2.w = sigma2.w,
       mu = mu,
       loglikelihood = loglikelihood)
}
.Random.seed <-
c(403L, 10L, 1952910500L, -1912918390L, 829351157L, -1253074097L, 
-747428778L, 1555670784L, -585063869L, -1794967071L, 1510579536L, 
-196422738L, -445361831L, -1573023893L, 430866138L, 817205148L, 
887403423L, -40019467L, 1656052876L, 286418274L, -410846915L, 
-2096229961L, -1836018082L, -477937768L, -12906757L, -2145524071L, 
312020776L, -1225415210L, 176695889L, 75497731L, 1391664274L, 
-1083452892L, -867767065L, -146931395L, 740357556L, -1630557670L, 
-1298673307L, 1311064063L, -1379985786L, 1305570832L, -1539689293L, 
1193693009L, -917609728L, 1715567070L, -810540503L, 389042651L, 
1456226026L, 494636172L, 864711631L, -1627124859L, -2069993412L, 
417425682L, 1716108813L, 1891642695L, 1832920878L, 1845156488L, 
874466187L, -201290903L, 1738065208L, -1982674970L, 1119035265L, 
261191443L, 224239810L, -1382073484L, 2019542903L, -1343545235L, 
1086750468L, -8110550L, -911182955L, -1945053137L, -1492815562L, 
-231894688L, -1902891101L, -745449279L, 1494878384L, 1261606862L, 
453577529L, 170754891L, 1669986554L, 233808380L, 748438079L, 
1236272917L, 2050372652L, -1065251582L, 2091460061L, -1543121705L, 
-1068490818L, -1634764296L, 421772763L, 92263481L, -760588920L, 
-1027423946L, -1417233551L, -337757789L, 2062451954L, -1498269628L, 
725064455L, 1817672413L, -980261932L, -1121082630L, -1569130939L, 
-1946431841L, 23229478L, 2072158384L, -1344349869L, -1674984335L, 
1644067232L, 2076620094L, -102476727L, 373022971L, 1966278858L, 
-670773716L, 1487794735L, -1309185307L, -737841252L, -1789916686L, 
200008813L, -710441561L, -166203442L, 1578593576L, -429781205L, 
-474090167L, 956021336L, 1960833798L, -358768927L, -1397540237L, 
-1649163550L, -1045979564L, 1488299863L, -641459379L, 489107172L, 
440233674L, -757518539L, -86145649L, 836854422L, -1309541824L, 
-1825027325L, -985785567L, 1244781584L, 242404334L, 104717209L, 
-1724873685L, -1541299942L, -901278756L, -1239706145L, -1469177419L, 
1636036428L, -117832286L, 2116977917L, -2036446729L, 270731038L, 
1759902168L, -1840147653L, -691448615L, -1008789400L, -1090359914L, 
493040913L, -456755517L, -204384558L, -377843868L, 1768125095L, 
-2060149763L, 907513332L, 748178138L, -1096627163L, -517250753L, 
51533510L, -456410288L, -1195777165L, -825617007L, -1902725184L, 
1576960414L, 1951525737L, -1628020837L, 271541930L, -1234038580L, 
-615214705L, 430466501L, 215689596L, -323003950L, 1656931661L, 
2067126279L, 1829173230L, 875036360L, -1687457333L, -1807842263L, 
-1363316104L, -1577287770L, 614250945L, -1599156653L, -499989246L, 
-1575738316L, 1961784759L, -2118067795L, 1300888260L, -1143693590L, 
457503061L, -1859393809L, 362697974L, 654992416L, -1500244509L, 
-1602845055L, 1280651504L, -2080723058L, -870932231L, 342837899L, 
1349274298L, 1404151228L, -1521840385L, -480574891L, 243676524L, 
167984578L, 1654308381L, 1772635799L, 590204926L, -1952121928L, 
-1811400805L, 47021561L, -412533432L, -1678182282L, 449627569L, 
-198358045L, 63150258L, 1109865732L, 153196103L, 354763037L, 
1382010772L, -2145176774L, 1101692805L, -857496592L, -1750993220L, 
-1484147436L, -1484156734L, 683801696L, 859502972L, 307423184L, 
-516609758L, -1274134872L, -629666932L, 397290684L, -1152819406L, 
-1579249088L, 920872116L, -983608248L, -141690438L, -58814960L, 
346355948L, 1212164436L, -996071662L, 1210189216L, -243569972L, 
1152788064L, -2048454718L, -1162820824L, -17662548L, -12578500L, 
1819963250L, -1116962512L, -1363697884L, 1595733304L, -1460086054L, 
-764962928L, -2058251524L, 1960992660L, 1602525250L, 1325616384L, 
1727851580L, 1477180752L, 1951197346L, 1856243528L, 1397102348L, 
1356790492L, 1816544082L, -1597881728L, 1358374996L, -2007745560L, 
-1765378310L, 1930045776L, 1699103532L, -1182861452L, 172946322L, 
-1329146272L, 220905356L, -334372000L, -850393118L, -138481752L, 
-1683628852L, -1420800516L, 93866930L, -1437619728L, -1118156348L, 
-2033667752L, -1323948678L, -147263056L, -448606276L, 764596308L, 
-1126661182L, 654599200L, -574362436L, 732326480L, -1119399966L, 
-1210911128L, 2110704204L, -1269186308L, 874049522L, -1903739648L, 
1068210036L, -663199736L, 1434449338L, -1212564016L, -2016246804L, 
-2121097324L, -84259886L, 898774880L, -1585826484L, -1188258528L, 
1492849538L, 579963624L, -2007075604L, 1478343228L, -1576920206L, 
-700757392L, -2037102556L, -790366408L, 312720538L, 1119990480L, 
-1473470276L, 1514284692L, -1883380094L, -1136763712L, 1076949308L, 
1538064528L, -982525406L, -1812613624L, 1702564108L, 1133109916L, 
1761951250L, -643240576L, -269652460L, -1026202072L, -428632774L, 
458483664L, -4832660L, -2059614412L, -1716850030L, 571767136L, 
2000970060L, 1693237216L, -1143050334L, -1871120472L, 1929075596L, 
65459900L, -2075387022L, 759171312L, -1613965244L, -717685416L, 
-12712902L, -1921735824L, 121963580L, 704887572L, -1661302334L, 
-1109913376L, -1282078852L, 1120284368L, -1900501342L, -2093939672L, 
579176716L, -1937939524L, -80810702L, 2078310336L, -575704524L, 
902728648L, -599554886L, -965300336L, -25762324L, -1582849964L, 
-1738920302L, -893363424L, 590372300L, 2082753504L, -479178046L, 
1782835752L, -809984340L, -215797060L, -919059982L, 1581861168L, 
-1232765788L, -1552076104L, -512914598L, -1630605936L, 834756988L, 
1424956692L, -1755551422L, -1423457536L, -918996036L, 1025010512L, 
-894836446L, 1538592072L, 632010764L, -991635364L, -1902750510L, 
-2115803136L, -1210061996L, 736956904L, 1253935098L, -762198320L, 
-489843412L, -929072908L, 1090893714L, -1979274144L, 944696588L, 
462030560L, -210793118L, 354255528L, -1376261556L, -1403394180L, 
-976508110L, 109276784L, -1924899516L, -1628623912L, -1415757702L, 
-1826990544L, -1474569540L, -1518053548L, -1972766142L, -2138809440L, 
2035071164L, 185741776L, -1998633630L, -1483021336L, -1367791668L, 
2107598588L, 769691634L, -1547322240L, 152389876L, 1596232840L, 
-1714244294L, -1047985328L, 1187458156L, -2079854444L, 787859794L, 
-807706912L, -832724788L, 13080736L, 1178261890L, 1279424744L, 
-1271240468L, -209983428L, -1779572622L, -1587812240L, 250602020L, 
-1765162440L, 1198979610L, -473669424L, -1775925700L, -147551084L, 
-864556243L, 166188954L, -825776808L, 606324929L, 736483187L, 
-1358218092L, 494717218L, 1467119671L, -1489524479L, 705412102L, 
697801284L, 1530655605L, -1223969057L, -272453064L, -453772778L, 
1141190179L, -1993782811L, -1177823230L, 20075120L, 1102641177L, 
73060875L, 1594931452L, 222151978L, 1043757919L, -2127770903L, 
1439257086L, -778756628L, -1436349123L, -1494485753L, 388691184L, 
272074894L, 1763410619L, -903996771L, 1322630730L, 1846953064L, 
-798090159L, 1268541251L, -852925884L, 772460786L, -1767859001L, 
1984792209L, -1708959722L, 1829151348L, -299272955L, -628724913L, 
-2122554040L, -1962108698L, 513707539L, 169947509L, -1280326318L, 
-457276000L, -1736248183L, -1855556421L, -1394571508L, -251498854L, 
8442895L, -432172967L, 823188718L, -1431218308L, -695581715L, 
-1056215913L, 1294841696L, -2036256962L, 941330443L, -1379592051L, 
388145594L, 414593016L, 1251494817L, 618833235L, 1573231732L, 
990834178L, 229892247L, -67538335L, -1672459610L, -1800072540L, 
-1682262187L, -1817838593L, -112170344L, -443838538L, -482941885L, 
-94749307L, 870131234L, -1101144816L, 36572345L, 1000943979L, 
1391452252L, -1121383030L, -796599937L, 714245321L, 823933278L, 
1883552268L, 796736989L, -1695750937L, 1958183184L, -804331474L, 
-726231333L, 1156378557L, 744180010L, 1972280776L, 1902437105L, 
-2114409373L, -636372892L, 2021471122L, 1865848935L, -411162575L, 
-1250150026L, 1735050900L, -959591259L, -1005317713L, 330730088L, 
-243273146L, 585848819L, 1876879189L, 1489754418L, -1481233152L, 
569667049L, 536901723L, -2130498900L, 182539450L, 598070895L, 
-1713486855L, -838571378L, 2035182428L, 367052749L, 1556883255L, 
2108831168L, -388425058L, 138875755L, -1821947155L, -1675509926L, 
315946008L, 1200345473L, -1462524109L, 1531220L, 1222776674L, 
940714999L, -1666919359L, -1921318714L, -1319490044L, 1728432565L, 
1002354335L, 556259064L, 694557398L, 1446219619L, -230344923L, 
599338178L, -814167632L, 176079321L, -1741172917L, -317807172L, 
1560073706L, 1375728159L, 750360745L, 1450119230L, -1273011924L, 
2133707645L, 1860633031L, 1907766448L, -653712690L, 1012406267L, 
-1301061539L, -1459411062L, 1625270184L, 829017361L, -106507261L, 
1723604740L, 1229517618L, 1191727495L, 2070091729L, 530071134L
)
.sigma.int <-
function(z, sgn, params, sd.t, lower.limit, upper.limit){  
  I <- function(u){
    arg.in = sgn*u + z - params[3];    
    return(exp(-(u/params[1])^2)*(2*(u^2/params[1]^3)*pnorm(arg.in/sd.t) - params[1]*arg.in/(2*sd.t^3)*dnorm(arg.in/sd.t)));
  }
  return(integrate(I,lower.limit,upper.limit,stop.on.error=1)$value);
}
.sigmab.int <-
function(z, sgn, params, sd.t, lower.limit, upper.limit){ 
  I <- function(u){
    arg.in = sgn*u + z - params[3];    
    return(exp(-(u/params[1])^2)*(-params[2]*arg.in*dnorm(arg.in/sd.t)/sd.t^3));
  }
  return(integrate(I,lower.limit,upper.limit,stop.on.error=1)$value);
}

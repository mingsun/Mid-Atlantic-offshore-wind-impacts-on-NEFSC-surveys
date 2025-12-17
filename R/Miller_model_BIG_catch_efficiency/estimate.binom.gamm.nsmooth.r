estimate.binom.gamm.nsmooth <- function(
  dat,
  mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ s(length, bs = 'cr', k = 10), 
  mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~ 1,
  lambda.station.factor,
  beta.station.var.factor,
  do_pred = FALSE, 
  use_beta_re = 0,
  use_beta_re_cor = 0,
  use_mean_smooth_re = 0,
  use_station_smooth_re = 0,
  new_dat,
  fit = TRUE,
  input
)
{
  require(mgcv)
  cl = match.call()

  if(missing(input))
  {
    temp.dat = gam(mu.mean.form, family = binomial, data = dat, fit = FALSE, control = list(scalePenalty=FALSE))  
    temp.fit = gam(mu.mean.form, family = binomial, data = dat, control = list(scalePenalty=FALSE))  
    pred.dat = predict(temp.fit, newdata = new_dat, type = "lpmatrix")
    
    mu.dat <- list(X = temp.dat$X)
    mu.dat$Dplus_diag = numeric()
    mu.dat$fp_mean_smooth_re = integer()
    mu.dat$lp_mean_smooth_re = integer()
    mu.dat$Xpred = pred.dat
    not_mean_smooth_re_ind = 1:NCOL(mu.dat$X)
    mean_smooth_re_ind = integer()
    station_smooth_re_ind = integer()
    
    if(length(temp.dat$smooth))#only 1 smooth allowed
    {
      for(i in 1:length(temp.dat$smooth))
      {
        S <- temp.dat$smooth[[i]]$S[[1]]
        UDU <- eigen(S)
        mu.dat$Dplus_diag <- c(mu.dat$Dplus_diag, UDU$val[which(UDU$val > 1e-10)])
        UF <- cbind(UDU$vec[,which(UDU$val <= 1e-10)])
        UR <- cbind(UDU$vec[,which(UDU$val > 1e-10)])
        ind = temp.dat$smooth[[i]]$first.para:temp.dat$smooth[[i]]$last.para
        mu.dat$X[,ind] = cbind(cbind(mu.dat$X[,ind]) %*% UF, cbind(mu.dat$X[,ind]) %*% UR)
        mu.dat$Xpred[,ind] = cbind(cbind(mu.dat$Xpred[,ind]) %*% UF, cbind(mu.dat$Xpred[,ind]) %*% UR)
        ind = ind[-(1:NCOL(UF))]
        mu.dat$fp_mean_smooth_re = c(mu.dat$fp_mean_smooth_re, min(ind))
        mu.dat$lp_mean_smooth_re = c(mu.dat$lp_mean_smooth_re, max(ind))
        mean_smooth_re_ind = c(mean_smooth_re_ind, ind)
      }
      not_mean_smooth_re_ind = not_mean_smooth_re_ind[-mean_smooth_re_ind]
    }
    
    temp.dat = gam(mu.station.form, family = binomial, data = dat, fit = FALSE, control = list(scalePenalty=FALSE))
    temp.fit = gam(mu.station.form, family = binomial, data = dat, control = list(scalePenalty=FALSE))  
    pred.dat = predict(temp.fit, newdata = new_dat, type = "lpmatrix")
    
    mu.dat$Z = temp.dat$X
    mu.dat$fp_station_smooth_re = integer()
    mu.dat$lp_station_smooth_re = integer()
    mu.dat$Zpred = pred.dat
    mu.dat$Dplus_diag_re = numeric()  
    not_station_smooth_re_ind = 1:NCOL(mu.dat$Z)
    
    
    if(length(temp.dat$smooth)) #only 1 smooth allowed
    {
      for(i in 1:length(temp.dat$smooth))
      {
        S <- temp.dat$smooth[[1]]$S[[1]]
        UDU <- eigen(S)
        mu.dat$Dplus_diag_re <- c(mu.dat$Dplus_diag_re, UDU$val[which(UDU$val > 1e-10)])
        UF <- cbind(UDU$vec[,which(UDU$val <= 1e-10)])
        UR <- cbind(UDU$vec[,which(UDU$val > 1e-10)])
        ind = temp.dat$smooth[[i]]$first.para:temp.dat$smooth[[i]]$last.para
        mu.dat$Z[,ind] = cbind(cbind(mu.dat$Z[,ind]) %*% UF, cbind(mu.dat$Z[,ind]) %*% UR)
        mu.dat$Zpred[,ind] = cbind(cbind(mu.dat$Zpred[,ind]) %*% UF, cbind(mu.dat$Zpred[,ind]) %*% UR)
        ind = ind[-(1:NCOL(UF))]
        mu.dat$fp_station_smooth_re = c(mu.dat$fp_station_smooth_re, min(ind))
        mu.dat$lp_station_smooth_re = c(mu.dat$lp_station_smooth_re, max(ind))
        station_smooth_re_ind = c(station_smooth_re_ind, ind)
      }
      not_station_smooth_re_ind = not_station_smooth_re_ind[-station_smooth_re_ind]
    }
    nstation = length(unique(dat$id))
    dat = list(
      id = as.integer(dat$id)-1, #because of C++ 
      n = apply(temp.dat$y,1,sum),
      y = temp.dat$y[,1],
      offst = dat$offst,
      X = cbind(mu.dat$X),
      beta_ns_ind = not_mean_smooth_re_ind - 1, #because of C++
      fp_mean_smooth_re = mu.dat$fp_mean_smooth_re-1, #because of C++
      lp_mean_smooth_re = mu.dat$lp_mean_smooth_re-1, #because of C++
      Z = cbind(mu.dat$Z),
      beta_re_ns_ind = not_station_smooth_re_ind - 1, #because of C++
      fp_station_smooth_re = mu.dat$fp_station_smooth_re-1, #because of C++
      lp_station_smooth_re = mu.dat$lp_station_smooth_re-1, #because of C++
      Dplus_diag = mu.dat$Dplus_diag,
      Dplus_diag_station = mu.dat$Dplus_diag_re,
      id_Dplus_diag_station = matrix(nrow = 0, ncol = 0),
      use_beta_re = use_beta_re,
      use_beta_re_cor = use_beta_re_cor,
      use_mean_smooth_re = use_mean_smooth_re,
      use_station_smooth_re = use_station_smooth_re,
      do_pred = do_pred,
      Xpred = mu.dat$Xpred,
      Zpred = mu.dat$Zpred 
    )
    if(length(mu.dat$Dplus_diag_re)) 
      dat$id_Dplus_diag_station = matrix(0:(length(mu.dat$Dplus_diag_re)-1), nstation, length(mu.dat$Dplus_diag_re), byrow = TRUE) #because of C++
    
    par = list(
      beta = rep(0, length(not_mean_smooth_re_ind)),
      mean_smooth_re = rep(0, length(mean_smooth_re_ind)),
      station_smooth_re = matrix(0, length(unique(dat$id)), length(station_smooth_re_ind)),
      lambda_par = rep(0,length(dat$fp_mean_smooth_re)),
      lambda_par_station = matrix(0,length(unique(dat$id)),length(dat$fp_station_smooth_re)),
      beta_re = matrix(0, length(unique(dat$id)), length(not_station_smooth_re_ind)),
      beta_re_var_pars = matrix(0, length(unique(dat$id)), length(not_station_smooth_re_ind)),
      beta_re_cor_pars = matrix(0, length(unique(dat$id)), length(not_station_smooth_re_ind)*(length(not_station_smooth_re_ind)-1)/2)
    )

    map = list()
    if(use_mean_smooth_re == 0) 
    {
      map$lambda_par = factor(rep(NA,length(par$lambda_par)))
      map$mean_smooth_re = factor(rep(NA, length(par$mean_smooth_re)))
    }
    
    if(use_station_smooth_re == 0) 
    {
      map$lambda_par_station = factor(rep(NA, length(par$lambda_par_station)))
      map$station_smooth_re = factor(rep(NA, length(par$station_smooth_re)))
    }
    else
    {
      if(!missing(lambda.station.factor))
      {
        un = unique(lambda.station.factor)
        lambda.station.map = cbind(match(lambda.station.factor, un))
        if(length(dat$fp_station_smooth_re)>1) 
        {
          for(i in 2:length(dat$fp_station_smooth_re)) 
            lambda.station.map = cbind(lambda.station.map, max(lambda.station.map) + match(lambda.station.factor, un))
        }
        map$lambda_par_station = factor(lambda.station.map)
      }
      else map$lambda_par_station = factor(rep(1:NCOL(par$lambda_par_station), each = NROW(par$lambda_par_station)))
    }
    
    map$beta_re_cor_pars =factor(rep(NA, length(par$beta_re_cor_pars)))
    if(use_beta_re == 0) 
    {
      map$beta_re_var_pars =factor(rep(NA, length(par$beta_re_var_pars)))
      map$beta_re = factor(rep(NA, length(par$beta_re)))
    }
    else 
    {
      if(!missing(beta.station.var.factor))
      {
        un = unique(beta.station.var.factor)      
        beta.station.var.map = cbind(match(beta.station.var.factor, un))
        if(NCOL(par$beta_re_var_pars)>1) 
        {
          for(i in 2:NCOL(par$beta_re_var_pars))
            beta.station.var.map = cbind(beta.station.var.map, max(beta.station.var.map) + match(beta.station.var.factor, un))
        }
        map$beta_re_var_pars = factor(beta.station.var.map)
        if(use_beta_re_cor == 1)
        {
          beta.station.cor.map = cbind(match(beta.station.var.factor, un))
          if(NCOL(par$beta_re_cor_pars)>1) 
          {
            for(i in 2:NCOL(par$beta_re_cor_pars))
              beta.station.cor.map = cbind(beta.station.cor.map, max(beta.station.cor.map) + match(beta.station.cor.factor, un))
          }
          map$beta_re_cor_pars = factor(beta.station.cor.map)
        }
      }
      else 
      {
        map$beta_re_var_pars = factor(rep(1:NCOL(par$beta_re_var_pars), each = NROW(par$beta_re_var_pars)))
        if(use_beta_re_cor == 1 & NCOL(par$beta_re_cor_pars)) map$beta_re_cor_pars = factor(rep(1:NCOL(par$beta_re_cor_pars), each = NROW(par$beta_re_cor_pars)))
      }
    }

    random = character()
    if(dat$use_beta_re == 1) random = c(random, "beta_re")
    if(dat$use_mean_smooth_re == 1) random = c(random, "mean_smooth_re")
    if(dat$use_station_smooth_re == 1) random = c(random, "station_smooth_re")

    temp = list(dat = dat, par = par, map = map, random = random)
    
  }
  else temp = input
  if(fit)
  {
    mod = MakeADFun(temp$dat, temp$par, DLL="binom_gamm_nsmooth", map = temp$map, random= temp$random, silent = TRUE)
    #return(list(mod,temp))
    mod = fit_tmb(mod, do.sdrep = FALSE)
    return(list(model.res=mod, input = temp, call = cl))
  }
  else return(list(model.res = NULL, input = temp, call = cl))
}

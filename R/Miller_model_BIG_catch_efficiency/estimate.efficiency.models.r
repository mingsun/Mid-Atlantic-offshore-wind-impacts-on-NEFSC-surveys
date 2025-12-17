estimate.efficiency.models <- function(dat, k=10, fit.models = c(paste0("bi",0:4),paste0("bb", 0:7)), 
  do_pred = TRUE, bootstrap = FALSE, predict.length)
{
  #bi 0-4 and bb 0-7 are the models estimated in Miller 2013 CJFAS.
  #bi5 is the same as bi4, but with day/night effects also on mean relative efficiency 
  require(mgcv)
  if(bootstrap)
  {
    ids = unique(dat$id)
    ids.boot = sample(ids, replace = TRUE)
    n.per.sta.boot = sapply(ids.boot, function(x) sum(dat$id == x))
    id.boot = rep(1:length(ids.boot), n.per.sta.boot)
    index.boot = unlist(sapply(ids.boot, function(x) which(dat$id == x)))
    dat = dat[index.boot,]
    dat$id = id.boot
  }
  #n <- c(NROW(cl.less), length(predict.length))
  if(missing(predict.length)) predict.length = seq(min(dat$length),max(dat$length),0.01)
  out = list(pred.length = predict.length)

  new.dat = data.frame(dn = factor(rep(c("day","night"), each = length(predict.length))), length = rep(predict.length,2))

  #if(do_pred) n_pred = n[2]
  #else n_pred = 0

  model.fits = list()
  if("bi0" %in% fit.models)
  {
    ############################
    #binomial model, no smoother
    model.fits$bi0 = estimate.binom.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ 1, 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~ 1, 
      do_pred = do_pred, 
      use_beta_re = 0,
      use_mean_smooth_re = 0,
      use_station_smooth_re = 0,
      new_dat = new.dat
    )
    print("bi0 done")
  }
  if("bi1" %in% fit.models)
  {
    ############################
    #binomial model, no smoother, with station-specific random effects
    model.fits$bi1 = estimate.binom.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ 1, 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~ 1, 
      do_pred = do_pred, 
      use_beta_re = 1,
      use_mean_smooth_re = 0,
      use_station_smooth_re = 0,
      new_dat = new.dat
    )
    print("bi1 done")
  }
  if("bi2" %in% fit.models)
  {
    ############################
    #normal gam with constant random effects for smoother
    model.fits$bi2 = estimate.binom.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ s(length, bs = 'cr', k = 10), 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~ 1, 
      do_pred = do_pred, 
      use_beta_re = 0,
      use_mean_smooth_re = 1,
      use_station_smooth_re = 0,
      new_dat = new.dat
    )
    print("bi2 done")
  }
  if("bi3" %in% fit.models)
  {
    ############################
    #random effects on intercept of mu for each station
    model.fits$bi3 = estimate.binom.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ s(length, bs = 'cr', k = 10), 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~ 1, 
      do_pred = do_pred, 
      use_beta_re = 1,
      use_mean_smooth_re = 1,
      use_station_smooth_re = 0,
      new_dat = new.dat
    )
    print("bi3 done")
  }  
  if("bi4" %in% fit.models)
  {
    ############################
    #random effects on smoothers for each station, but lambda is constant
    model.fits$bi4 = estimate.binom.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ s(length, bs = 'cr', k = 10), 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~  s(length, bs = 'cr', k = 10), 
      do_pred = do_pred, 
      use_beta_re = 1,
      use_beta_re_cor = 1,
      use_mean_smooth_re = 1,
      use_station_smooth_re = 1,
      new_dat = new.dat
    )
    print("bi4 done")
  }
  if("bi5" %in% fit.models)
  {
    ############################
    #day/night effects on relative efficiency and lambda, random effects on smoothers for each station
    model.fits$bi5 = estimate.binom.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ dn + s(length, by = dn, bs = 'cr', k = 10), 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~  s(length, bs = 'cr', k = 10), 
      do_pred = do_pred, 
      use_beta_re = 1,
      use_beta_re_cor = 1,
      use_mean_smooth_re = 1,
      use_station_smooth_re = 1,
      new_dat = new.dat
    )
    print("bi5 done")
  }

  if("bb0" %in% fit.models)
  {
    ############################
    #no smoother in mu or phi parameter
    model.fits$bb0 = estimate.betabin.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ 1, 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~  1, 
      phi.mean.form = cbind(recnumlen.ch,recnumlen.rh)  ~ 1, 
      do_pred = do_pred, 
      use_beta_re = 0,
      use_beta_re_cor = 0,
      use_mean_smooth_re = 0,
      use_station_smooth_re = 0,
      use_phi_mean_smooth_re = 0,
      new_dat = new.dat
    )
    print("bb0 done")
  }
  if("bb1" %in% fit.models)
  {
    ############################
    #no smoother in mu or phi parameter, station-specific random intercepts in mu
    model.fits$bb1 = estimate.betabin.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ 1, 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~  1, 
      phi.mean.form = cbind(recnumlen.ch,recnumlen.rh)  ~ 1, 
      do_pred = do_pred, 
      use_beta_re = 1,
      use_beta_re_cor = 0,
      use_mean_smooth_re = 0,
      use_station_smooth_re = 0,
      use_phi_mean_smooth_re = 0,
      new_dat = new.dat
    )
    print("bb1 done")
  }
  if("bb2" %in% fit.models)
  {
    ############################
    #mean smoother for mu, no smoother in phi parameter
    model.fits$bb2 = estimate.betabin.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ s(length, bs = 'cr', k = 10), 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~  1, 
      phi.mean.form = cbind(recnumlen.ch,recnumlen.rh)  ~ 1, 
      do_pred = do_pred, 
      use_beta_re = 0,
      use_beta_re_cor = 0,
      use_mean_smooth_re = 1,
      use_station_smooth_re = 0,
      use_phi_mean_smooth_re = 0,
      new_dat = new.dat
    )
    print("bb2 done")
  }
  if("bb3" %in% fit.models)
  {
    ############################
    #with smoother in phi parameter
    model.fits$bb3 = estimate.betabin.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ s(length, bs = 'cr', k = 10), 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~  1, 
      phi.mean.form = cbind(recnumlen.ch,recnumlen.rh)  ~ s(length, bs = 'cr', k = 10), 
      do_pred = do_pred, 
      use_beta_re = 0,
      use_beta_re_cor = 0,
      use_mean_smooth_re = 1,
      use_station_smooth_re = 0,
      use_phi_mean_smooth_re = 1,
      new_dat = new.dat
    )
    print("bb3 done")
  }
  if("bb4" %in% fit.models)
  {
    ############################
    #with station-specific random effects on intercept of mu, no smoother on phi
    model.fits$bb4 = estimate.betabin.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ s(length, bs = 'cr', k = 10), 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~  1, 
      phi.mean.form = cbind(recnumlen.ch,recnumlen.rh)  ~ 1, 
      do_pred = do_pred, 
      use_beta_re = 1,
      use_beta_re_cor = 0,
      use_mean_smooth_re = 1,
      use_station_smooth_re = 0,
      use_phi_mean_smooth_re = 0,
      new_dat = new.dat
    )
    print("bb4 done")
  }
  if("bb5" %in% fit.models)
  {
    ############################
    #with station-specific random effects on intercept of mu, smoother on phi
    model.fits$bb5 = estimate.betabin.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ s(length, bs = 'cr', k = 10), 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~  1, 
      phi.mean.form = cbind(recnumlen.ch,recnumlen.rh)  ~ s(length, bs = 'cr', k = 10), 
      do_pred = do_pred, 
      use_beta_re = 1,
      use_beta_re_cor = 0,
      use_mean_smooth_re = 1,
      use_station_smooth_re = 0,
      use_phi_mean_smooth_re = 1,
      new_dat = new.dat
    )
    print("bb5 done")
  }
  if("bb6" %in% fit.models)
  {
    ############################
    #with station-specific random effects on fixed and random effects portion of mu smoother, no smoother on phi
    model.fits$bb6 = estimate.betabin.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ s(length, bs = 'cr', k = 10), 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~  s(length, bs = 'cr', k = 10), 
      phi.mean.form = cbind(recnumlen.ch,recnumlen.rh)  ~ 1, 
      do_pred = do_pred, 
      use_beta_re = 1,
      use_beta_re_cor = 1,
      use_mean_smooth_re = 1,
      use_station_smooth_re = 1,
      use_phi_mean_smooth_re = 0,
      new_dat = new.dat
    )
    print("bb6 done")
  }
  if("bb7" %in% fit.models)
  {
    ############################
    #with station-specific random effects on fixed and random effects portion of mu smoother, smoother on phi
    model.fits$bb7 = estimate.betabin.gamm.nsmooth(
      dat = dat,
      mu.mean.form = cbind(recnumlen.ch,recnumlen.rh) ~ s(length, bs = 'cr', k = 10), 
      mu.station.form = cbind(recnumlen.ch,recnumlen.rh)  ~  s(length, bs = 'cr', k = 10), 
      phi.mean.form = cbind(recnumlen.ch,recnumlen.rh)  ~ s(length, bs = 'cr', k = 10), 
      do_pred = do_pred, 
      use_beta_re = 1,
      use_beta_re_cor = 1,
      use_mean_smooth_re = 1,
      use_station_smooth_re = 1,
      use_phi_mean_smooth_re = 1,
      new_dat = new.dat
    )
    print("bb7 done")
  }
  
  out$dat = dat
  out$new.dat = new.dat
  out$model.fits = model.fits
  return(out) 
}

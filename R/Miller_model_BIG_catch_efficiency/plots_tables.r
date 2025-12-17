library(TMB)
library(Hmisc)
load("data/survey/stock.strata.RData")
source("code/get.best.r")
source("code/make.all.AIC.table.fn.r")
parentdir = getwd()

#datasets by species are already compiled by code in .gitignored code/gather_data

species.order = c(1:6,10,7:9)
sp.info$sp.pretty.names[species.order]

for(i in 1:NROW(sp.info))
{
  sp = sp.info$sp.names[i]
  combined.data = readRDS(paste0(parentdir, "/data/", sp, "_data.RDS"))
  dat = combined.data$catch.length  
  x = aggregate(cbind(expnumlen.rh,expnumlen.ch) ~ length * day.night, data = combined.data$catch.length, FUN = sum) 
  colnames(x)[3:4] = c("expnumlen.rockhopper", "expnumlen.chainsweep")
  write.csv(x, file = paste0(parentdir,"/results/", sp, "_twin_trawl_length_frequencies.csv"))
}


source("code/make.nobs.table.fn.r")
nobs.table = make.nobs.table.fn()
nobs.table = nobs.table[species.order,]
nobs.table = format.df(nobs.table, big.mark = ",", numeric.dollar=F)
x = latex(nobs.table, file = paste0(parentdir,"/paper/data_table.tex"), 
  cgroup = c("Paired Tows", "Captured", "\\thead{Both Gears \\\\ Measured}", "\\thead{Chainsweep \\\\ Measured}", "\\thead{Rockhopper \\\\ Measured}"), first.hline.double = FALSE,
  n.cgroup = c(3,1,3,3,3), 
  col.just = rep("r",NCOL(nobs.table)),
  cgroupTexCmd="normalfont",#numeric.dollar = FALSE, #cellTexCmds = cell.format,
  rowlabel = 'Species', table.env = FALSE)#, rowlabel.just = "c")

#make summary table of AIC for each model/species
species.get.best = lapply(sp.info$sp.names, get.best, pdir = parentdir)
all.models = unique(unlist(sapply(species.get.best, function(y) names(y$np))))

n.rep = c(5,8)
first.names = paste0(rep(c("bi", "bb"), n.rep), c(0:4,0:7))

out = cbind.data.frame(
  Model = paste0(rep(c("BI", "BB"), n.rep), "$_{", c(0:(n.rep[1]-1),0:(n.rep[2]-1)), "}$"),
  "$\\log\\left(\\rho\\right)$" = paste0("$\\sim$ ", c(
  1, 
  "1 + 1$|$pair", 
  "s(length)", 
  "s(length) + 1$|$pair",
  "s(length) + s(length)$|$pair",
  1,
  "1 + 1$|$pair", 
  "s(length)", 
  "s(length)", 
  "s(length) + 1$|$pair",
  "s(length) + 1$|$pair",
  "s(length) + s(length)$|$pair",
  "s(length) + s(length)$|$pair")),
  "$\\log\\left(\\phi\\right)$" = c(
  "--", 
  "--", 
  "--", 
  "--", 
  "--", 
  paste0("$\\sim$ ", c(
  1, 
  1, 
  1, 
  "s(length)", 
  1, 
  "s(length)", 
  1, 
  "s(length)"))),  
#  "log-likelihood" = round(nll.rep,2),
  "$n_p$" = species.get.best[[1]]$np[first.names]#,
  #aic = round(aic.rep,2)
#  "log-likelihood" = round(-sapply(x$model.fits, function(x) x$model.res$opt$obj),2),
#  "$P$" = sapply(x$model.fits, function(x) length(x$model.res$opt$par)),
#  aic = round(sapply(x$model.fits, function(x) 2*(x$model.res$opt$obj + length(x$model.res$opt$par))),2)
)
#out[["$\\Delta$AIC"]] = round(out$aic - min(out$aic, na.rm = TRUE), 2)
out$Description = c(
  "population-level mean for all observations",
  "population- and random station-level $\\rho$",
  "population-level smooth size effect on $\\rho$",
  "population-level smooth size effect and random station-level intercept for $\\rho$",
  "population-level and random station-level smooth size effects for $\\rho$",
  "population-level $\\rho$ and $\\phi$",
  "population-level and random station-level intercept for $\\rho$ and population-level $\\phi$",
  "population-level smooth size effect on $\\rho$ and population-level $\\phi$",
  "population-level smooth size effect on $\\rho$ and $\\phi$",
  "population-level smooth size effect and random station-level intercept for $\\rho$ and population-level $\\phi$",
  "population-level smooth size effect on $\\rho$ and $\\phi$ and random station-level intercepts for $\\rho$",
  "population-level and random station-level smooth size effects on $\\rho$ and population-level $\\phi$",
  "population-level and random station-level smooth size effects on $\\rho$ and population-level smooth size effects on $\\phi$")#,
#  "population-level smooth size and day/night effects and random station-level smooth size effects on $\\rho$ and population-level smooth size effects on $\\phi$",
#  "BB$_8$ with no correlation of station-specific random effects",
#  "BB$_9$ with interaction of population-level day/night and smooth size effects")

x = latex(out, file = paste0(parentdir, "/paper/base_model_description.tex"), 
  table.env = FALSE, col.just = c("l","l","l","r","p{0.5\\textwidth}"), rowname = NULL)

aic.table = make.all.AIC.table.fn(sp.info$sp.names, sp.info$sp.pretty.names, pdir = parentdir)
aic.table = t(aic.table[,-12])
aic.table = aic.table[c(2:7,11,8:10),]
aic.table[11,] = round(aic.table[11,], 0)
cell.format = matrix("", NROW(aic.table), NCOL(aic.table))
cell.format[cbind(1:10, apply(aic.table[1:10,], 1, function(x) which(x == 0)))] = "bfseries"
x = aic.table
x[] = as.character(x)
x = latex(x, file = paste0(parentdir,"/paper/base_model_compare.tex"), 
  #rgroup = c("", ""), first.hline.double = FALSE,
  #n.rgroup = c(length(sp.info$sp.names)-1,1), 
  col.just = rep("r",NCOL(x)),
  numeric.dollar = FALSE, cellTexCmds = cell.format,
  rowlabel = '', table.env = FALSE)#, rowlabel.just = "c")


source("code/make_best_model_table.r") #makes a data.frame called finalout
ind = c(1:6,10,7:9)
temp = finalout[,-1]
cell.format = matrix("", NROW(temp), NCOL(temp))
cell.format[temp[[5]] == "0",5] = "bfseries"
x = latex(temp, file = paste0(parentdir,"/paper/best_model_compare.tex"), 
  rgroup = sp.info$sp.pretty.names[ind],
  n.rgroup = rep(3,length(ind)), 
  col.just = c("l", "l", "l", rep("r",2)),
  numeric.dollar = FALSE, cellTexCmds = cell.format,
  rowlabel = '', rowname = rep("",NROW(temp)), table.env = FALSE)#, rowlabel.just = "c")



#make plots of results for each species including bootstrap-based confidence intervals.
#source("get.best.r")
n.good.boot = c()
for(i in species.order) 
{
  boot.pred.eta = readRDS(paste0("results/", sp.info$sp.names[i], "_boot_pred_eta_0.RDS"))
  for(j in 1:9) boot.pred.eta = rbind(boot.pred.eta, readRDS(paste0("results/", sp.info$sp.names[i], "_boot_pred_eta_",j,".RDS")))
  n.good.boot[i] = sum(!is.na(boot.pred.eta[,1]))
}
names(n.good.boot) = sp.info$sp.names[species.order]

boot.pred.eta = readRDS(paste0("results/", sp.name, "_boot_pred_eta_0.RDS"))
#temp = boot.pred.eta
for(j in 1:9) boot.pred.eta = rbind(boot.pred.eta, readRDS(paste0("results/", sp.name, "_boot_pred_eta_",j,".RDS")))
temp = readRDS(paste0("/home/tmiller2/work/paired_tow_studies/R/2020/",sp.name,"/boot.pred.eta.RDS"))
apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.025, na.rm = TRUE)
apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.5, na.rm = TRUE)
apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.975, na.rm = TRUE)

ymax.sp = c(6,6,15,6,10,10,6,50,10,10)
#png(filename = paste0("paper/sp_rho_plot.png"), width = 8*144, height = 12*144, res = 144, pointsize = 12, family = "Times")
cairo_pdf('paper/sp_rho_plot.pdf', family = "Times", height = 12, width = 8)
source("code/plot.results.r")
par(oma = c(3,3,0,0), mar = c(2,2,3,1), mfcol = c(5,2))
for(i in species.order) 
{
  print(i)
  dn.legend.loc = NULL
  if(i == species.order[length(species.order)]) dn.legend.loc = "topright"
  plot.results(sp.info$sp.names[i], pdir = parentdir, ymax = ymax.sp[i], dn.legend.loc = dn.legend.loc) 
  mtext(side = 3, sp.info$sp.pretty.names[i], line = 0, cex = 1.5)
}
mtext(side = 1, 'Length (cm)', line = 1, cex = 1.5, outer = TRUE)
mtext(side = 2, "Relative Catch Efficiency (Chain:Rockhopper)", line = 1, cex = 1.5, outer = TRUE)
dev.off()

##########BIOMASS INDICES

setwd(parentdir)
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)
for( i in 1:length(stocks)) {
  load(paste0("results/", stocks[i], "_N.W.RData"))
  s.ind = match(as.character(survey.years),colnames(all.spring.N.W))
  f.ind = match(as.character(survey.years),colnames(all.fall.N.W))
  y = cbind(spring = all.spring.N.W[2,s.ind], fall = all.fall.N.W[2,f.ind])/1000
  write.csv(y, file = paste0("results/", stocks[i], "_bigelow_biomass.csv"))
  y = cbind(spring = all.spring.N.W[3,s.ind], fall = all.fall.N.W[3,f.ind])
  write.csv(y, file = paste0("results/", stocks[i], "_bigelow_n_per_tow.csv"))
  y = cbind(spring = all.spring.N.W[4,s.ind], fall = all.fall.N.W[4,f.ind])
  write.csv(y, file = paste0("results/", stocks[i], "_bigelow_kg_per_tow.csv"))
}

source("code/survey/get_biomass_boots.r")
survey.years = 2009:2019
all.boots = lapply(1:length(x), function(y) get_biomass_boots(sp.i = x[y], stock.i = y,years = survey.years))

sp.name = "fluke"
boot.pred.eta = readRDS(paste0("results/", sp.name, "_boot_pred_eta_0.RDS"))
for(j in 1:9) boot.pred.eta = rbind(boot.pred.eta, readRDS(paste0("results/", sp.name, "_boot_pred_eta_",j,".RDS")))
temp  = all.boots[[1]][[1]][1,]
which(temp == max(temp)) #710
apply(boot.pred.eta[,1:20], 2, function(x) which(x == max(x,na.rm=T))) #710

#remove the extreme value for fluke relative efficiency
all.boots[[1]] = lapply(all.boots[[1]], function(y) y[,-710])
#fall 2017 survey is very poorly sampled for fluke strata 
all.boots[[1]]$fall.calib.boots[which(survey.years==2017),] = NA
all.boots[[1]]$fall.uncalib.boots[which(survey.years==2017),] = NA


cv.fn = function(x) {
  sd(x,na.rm=T)/mean(x,na.rm=T)
  }

biomass.cvs = lapply(all.boots, function(y){
  fall.calib.cv = apply(y$fall.calib.boots,1,cv.fn)
  fall.uncalib.cv = apply(y$fall.uncalib.boots,1,cv.fn)
  spring.calib.cv = apply(y$spring.calib.boots,1,cv.fn)
  spring.uncalib.cv = apply(y$spring.uncalib.boots,1,cv.fn)
  return(cbind(fall.calib.cv, fall.uncalib.cv,spring.calib.cv,spring.uncalib.cv))
  })

lapply(1:length(all.boots), function(y){

  x = all.boots[[y]]
  out = cbind(apply(x$fall.calib.boots,1, cv.fn)) #cv
  out = cbind(out, apply(x$fall.calib.boots,1, quantile, probs = 0.025, na.rm = TRUE)/1000)
  out = cbind(out, apply(x$fall.calib.boots,1, quantile, probs = 0.975, na.rm = TRUE)/1000)
  out = cbind(out, apply(x$spring.calib.boots,1, cv.fn)) #cv
  out = cbind(out, apply(x$spring.calib.boots,1, quantile, probs = 0.025, na.rm = TRUE)/1000)
  out = cbind(out, apply(x$spring.calib.boots,1, quantile, probs = 0.975, na.rm = TRUE)/1000)
  print(dim(out))
  colnames(out) = c("cv.fall", "lo.fall", "hi.fall", "cv.spring", "lo.spring", "hi.spring")
  print(out)
  write.csv(round(out,2), file = paste0("results/", stocks[y], ".cv.biomass.csv"))
  temp = x$fall.uncalib.boots/x$fall.calib.boots
  out = cbind(apply(temp,1, cv.fn)) #cv
  out = cbind(out, apply(temp,1, quantile, probs = 0.025, na.rm = TRUE))
  out = cbind(out, apply(temp,1, quantile, probs = 0.975, na.rm = TRUE))
  temp = x$spring.uncalib.boots/x$spring.calib.boots  
  out = cbind(out, apply(temp,1, cv.fn)) #cv
  out = cbind(out, apply(temp,1, quantile, probs = 0.025, na.rm = TRUE))
  out = cbind(out, apply(temp,1, quantile, probs = 0.975, na.rm = TRUE))
  colnames(out) = c("cv.fall", "lo.fall", "hi.fall", "cv.spring", "lo.spring", "hi.spring")
  write.csv(round(out,2), file = paste0("results/", stocks[y], ".cv.biomass.ratio.csv"))
  })




use.stock.names.2 = c(
  "Summer flounder", 
  "American plaice", 
  "GB-GOM windowpane", 
  "SNE-MAB windowpane", 
  "GB winter flounder", 
  "GOM winter flounder", 
  "SNE winter flounder", 
  "GB yellowtail flounder", 
  "SNE-MA yellowtail flounder", 
  "CC-GOM yellowtail flounder", 
  "Witch flounder",
  "Northern goosefish",
  "Southern goosefish",
  "Barndoor skate",
  "Thorny skate",
  "Northern red hake",
  "Southern red hake",
  "GOM cod",
  "GB cod"
  )
stock.definition.table = cbind(Stock = c(
  "Summer flounder",
  "American Plaice",
  "Georges Bank-Gulf of Maine (GB-GOM) windowpane",
  "Southern New England-Mid-Atlantic Bight (SNE-MAB) windowpane",
  "Georges Bank (GB) winter flounder",
  "Gulf of Maine (GOM) winter flounder",
  "Southern New England (SNE) winter flounder",
  "GB yellowtail flounder",
  "Southern New England-Mid-Atlantic (SNE-MA) yellowtail flounder",
  "Cape Cod-Gulf of Maine (CC-GOM) yellowtail flounder",
  "Witch flounder",
  "Northern red hake",
  "Southern red hake",
  "Northern goosefish",
  "Southern goosefish",
  "Barndoor skate",
  "Thorny skate"))
x = latex(stock.definition.table, file = paste0("paper/stock_definition_table.tex"), 
  first.hline.double = FALSE,
  col.just = rep("l",2),
  table.env = FALSE)#, rowlabel.just = "c")


stock.order = c(1:11,16,17,12:15)
cv.ratios = t(sapply(biomass.cvs, function(x) c(mean(x[,1]/x[,2],na.rm=T), mean(x[,3]/x[,4],na.rm=T))))
rownames(cv.ratios) = use.stock.names.2
colnames(cv.ratios) = c("Fall","Spring")
x = latex(round(cv.ratios[stock.order,2:1],2), file = paste0(parentdir,"/paper/stock_cv_ratios.tex"), 
  cgroup = c("\\thead{Average CV Ratio \\\\ Calibrated:Uncalibrated}"), first.hline.double = FALSE,
  n.cgroup = 2,
  collabel.just= rep("r",2), 
  col.just = rep("r",2),
  cgroupTexCmd="normalfont",#numeric.dollar = FALSE, #cellTexCmds = cell.format,
  rowlabel = 'Stock', table.env = FALSE)#, rowlabel.just = "c")

mean.cor = t(sapply(all.boots, function(y) {
  sapply(y, function(z) {
    res = cor(t(z))
    res = mean(res[lower.tri(res)], na.rm = TRUE)
    return(res)
    })
}))
rownames(mean.cor) = use.stock.names.2

x = latex(round(mean.cor[stock.order,c(3,1)],2), file = paste0(parentdir,"/paper/stock_mean_cor.tex"), 
  first.hline.double = FALSE,
  n.cgroup = 2,
  colheads = c("Spring","Fall"),
#  collabel.just= rep("r",2), 
#  col.just = rep("r",2),
  rowlabel = 'Stock', table.env = FALSE)#, rowlabel.just = "c")


stock.order = c(1:11,15:17,14,12:13) #change to get similar magnitude in the same rows
stock.order = c(1,3:10,2,11,15:17,14,12:13) #change to get similar magnitude in the same rows
cairo_pdf('paper/stock_biomass_plot.pdf', family = "Times", height = 12, width = 8)
source("code/survey/plot_biomass.r")
plot_biomass(stocks[stock.order], use.stock.names[stock.order])
dev.off()

source("code/survey/plot_biomass.r")
plot_biomass(stocks[stock.order], use.stock.names[stock.order])

cairo_pdf('paper/stock_biomass_efficiency_plot.pdf', family = "Times", height = 12, width = 8)
source("code/survey/plot_biomass_efficiency.r")
plot_biomass_efficiency(stocks[stock.order], use.stock.names[stock.order])
dev.off()

make.nobs.table.fn = function(sp.dat=sp.info)
{

  library(Hmisc)
  
  out = matrix(NA, NROW(sp.dat), 13)
  print(out)
  print(dim(out))
  rownames(out) = sp.dat$sp.pretty.names
  colnames(out) = c("Total", "Day", "Night", "Total", rep(c("Total", "Day", "Night"),3))
  for(i in 1:NROW(out)) 
  {
    tt.dat = readRDS(paste0("data/", sp.info$sp.names[i], "_data.RDS"))$catch.length

    #number of stations
    out[i,1] = aggregate(id ~ 1, data = tt.dat, FUN = function(x) length(unique(x)))[1,1]
    #number of stations day, night
    out[i,2:3] = aggregate(id ~ day.night, data = tt.dat, FUN = function(x) length(unique(x)))[,2]
    #number of fish caught across all stations
    out[i,4] = aggregate(expnumlen.ch + expnumlen.rh ~ 1, data = tt.dat, FUN = sum)[1,1]
    #number of fish measured for length across all stations
    out[i,5] = aggregate(recnumlen.ch + recnumlen.rh ~ 1, data = tt.dat, FUN = sum)[1,1]
    #number of fish measured for length across all stations day, night
    out[i,6:7] = aggregate(recnumlen.ch + recnumlen.rh ~ day.night, data = tt.dat, FUN = sum)[,2]
    #number of fish measured for length across all stations in chainsweep gear
    out[i,8] = aggregate(recnumlen.ch ~ 1, data = tt.dat, FUN = sum)[1,1]
    #number of fish measured for length across all stations in chainsweep gear day, night
    out[i,9:10] = aggregate(recnumlen.ch ~ day.night, data = tt.dat, FUN = sum)[,2]
    #number of fish measured for length across all stations in rockhopper gear
    out[i,11] = aggregate(recnumlen.rh ~ 1, data = tt.dat, FUN = sum)[1,1]
    #number of fish measured for length across all stations in rockhopper gear day, night
    out[i,12:13] = aggregate(recnumlen.rh ~ day.night, data = tt.dat, FUN = sum)[,2]
  }
  
  #x = latex(out, file = '', 
  #  rowlabel = 'Species', table.env = FALSE, rowlabel.just = "c")
  #write.csv(out, file= "nobs_table.csv")
  return(out)
}

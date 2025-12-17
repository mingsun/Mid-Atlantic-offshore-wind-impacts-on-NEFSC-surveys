

get_stMean <- function(towData, stMeanCol){
  
  # Calculate the unweighted means by stratum
  meanByStratum <- towData %>%
    
    # Group the data according to station information and sum up at that level
    # of aggregation. If length data are input, this will put the data on
    # a tow level; if tow level data are input it will simply change the
    # column of interest to X
    group_by(comname, year, Season, stratum, StratumAreaSqNMiles, station) %>%
    summarize(X = sum({{stMeanCol}}),
              .groups = 'drop') %>%
    
    # remove the tow-level grouping and group by stratum
    group_by(comname, year, Season, stratum, StratumAreaSqNMiles) %>%
    
    # Calculate the unweighted mean and variance by stratum
    summarize(meanX = mean(X),
              n = n(),
              varX = sum((X - meanX)^2) / (n - 1),
              .groups = 'drop')

  
  # Calculate the stratified mean
  stMean <- meanByStratum %>%
    
    # Remove the stratum-level grouping
    group_by(comname, year, Season) %>%
    
    # Calculate the stratum weighting as a proportion
    mutate(W = StratumAreaSqNMiles / sum(StratumAreaSqNMiles)) %>%
    
    # Calculate the weighted mean and variance according to Cochran (1977)
    # methods
    summarize(stMeanX = weighted.mean(meanX, w = W),
              stVarX = sum(W^2 * varX / n, na.rm = TRUE),
              stSeX = sqrt(stVarX),
              .groups = 'drop') %>%
    
    # Calculate the confidence bounds and the CV
    mutate(up95X = qnorm(0.975, mean = stMeanX, sd = stSeX),
           lo95X = qnorm(0.025, mean = stMeanX, sd = stSeX),
           lo95X = ifelse(lo95X < 0, yes = 0, no = lo95X),
           CVX = stSeX / stMeanX) %>%
    
    # Rename columns for easy reading
    rename(Year = year, `Stratified Mean` = stMeanX) %>%
    
    # Make sure any missing years (e.g., 2020) are included in the data
    # set and list all values other than species and season as NA
    complete(Year = min(Year):max(Year), 
             fill = list(comname = unique(.$comname),
                         Season = unique(.$Season)))

  return(stMean)
  
}





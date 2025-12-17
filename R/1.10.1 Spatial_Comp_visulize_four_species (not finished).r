
### some metrics from  Li et al. 2017, consider Taylor diagram

# correlation coefficient

# Average relative error
consider 0 value here

# average absolute error

# root mean square error

# reliability index

# modeling efficiency



# Taylor diagram

library(openair)

TaylorDiagram(subset(full_df, OWF == "INSIDE"), obs = "obs", mod = "mod", group = "model")
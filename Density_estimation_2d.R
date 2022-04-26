#This program runs the SMC algorithm under the settings of Section 5.1.1
#The package "SMCMP" needs to be installed before running the following codes

library(SMCMP)

#The model settings need to be imported from "models.R"
source("LOCATION/HMPT_experiments/models.R")

#Parameters
#We use one combination of parameter values for illustration
model_name = "2D_Uniforms" #model name
d = 2 #dimension
n = 1000 #sample size 
max_K = 10 #maximum resolution
eta = 10^(-1) #eta in the location prior
NL = 32 #NL-1 = # grid points for the location prior
M = 1000 # # particles
ran_seed = 1 # random seed

#In the paper, the following parameter values are used\
# model_name = "2D_Uniforms", "2D_Clusters", and "2D_Smooth"
# n = 500, 600, ..., 1000
# eta = 0.1, 0.01
#ran_seed = 1, 2, ..., 100

#Generate data
set.seed(ran_seed)
data = Generate_data(model_name, d, n)

#Make the grid points to visualize the estimated densities and compute the KL divergence
#In the manuscript, a part of the boundaries are removed for better visualization
num_grid_density_2D = 100
by = 1/num_grid_density_2D
start = by/2
end = 1 - by/2

x.grid = seq(start,end,by=by)
y.grid = seq(start,end,by=by)
grid_points = as.matrix(expand.grid(x.grid,y.grid))

#Run the SMC
out_APT = smc.apt(data, grid_points, max.resol = max_K, n.particle = M, n.grid.L = NL-1, eta.R = eta)
densities = out_APT$pred.density

#Plot the estimated density
library(ggplot2)
library(viridis)
densities.df = data.frame(x1=grid_points[,1],x2=grid_points[,2],density = densities)
p = ggplot() + geom_tile(data = densities.df, aes(x=x1, y=x2, fill=density)) + scale_fill_viridis(discrete=FALSE)
print(p)

#Plot the MAP tree
tree_left = out_APT$MAPtree.left
tree_right = out_APT$MAPtree.right

n_all_node = ncol(tree_left)
N = (n_all_node-1)/2

x_all = numeric(N)
y_all = numeric(N)
x_end_all = numeric(N)
y_end_all = numeric(N)

for(j in 1:N){
  start = tree_left[,2*j+1]
  end   = tree_right[,2*j]
  
  x_all[j] = start[1]
  y_all[j] = start[2]
  
  x_end_all[j] = end[1]
  y_end_all[j] = end[2]
}

x_st = c(0,0,1,1)
x_end = c(0,1,1,0)
y_st = c(0,1,1,0)
y_end = c(1,1,0,0)

print(p + geom_segment(aes(x = x_all, y = y_all, xend = x_end_all, yend = y_end_all)) +
        geom_segment(aes(x = x_st, y = y_st, xend = x_end, yend = y_end)) + xlab("x1") + ylab("x2"))

#Compue the KL divergence from the true distribution
true_densities = True_density(model_name, d, grid_points)
SMALL_NUMBER = 1e-100
KL = mean(true_densities * (log(true_densities + SMALL_NUMBER) - log(out_APT$pred.density + SMALL_NUMBER)))

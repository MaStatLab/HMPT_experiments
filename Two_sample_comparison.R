#This program runs the SMC algorithm under the settings of Section 5.2
#The package "SMCMP" needs to be installed before running the following codes

library(SMCMP)

#The model settings need to be imported from "models.R"
source("LOCATION/HMPT_experiments/models.R")

#Parameters
#We use one combination of parameter values for illustration
model_name = "Location_shift" #model name
# model names: "Location_shift", "Dispersion_difference", "Correlation"
d = 50 #dimension
n = 4000 #total sample size (= n1 + n2)
max_K = 15 #maximum resolution
eta = 10^(-1) #eta in the location prior (0.1 or 0.01 in the experiment)
NL = 32 #NL-1 = # grid points for the location prior
M = 1000 # # particles
ran_seed = 1 # random seed (set to 1,...,200 in the experiment)
n_groups = 2 # # groups
same_or_not = 0 #0: two groups are generated from the same distribution. 1: different distributions

#Also we can fix the value of L(A) to 0.5 by setting eta to an extremely large value such as 10^100

#Generate data
set.seed(ran_seed)
groups = rep(c(1,n_groups), each = n/2) #group labels
data = Generate_data_two_sample(model_name, d, groups, n, n_groups, same_or_not)

#Run the SMC
out_MRS = smc.mrs(data, groups, max.resol = max_K, n.particle = M, n.grid.L = NL, eta.R = eta)

print("Posterior prob of the two groups being the same")
print(out_MRS$post.null.global)

#Draw the nodes in which the largest difference is found

#preparation
Indices_included = function(X, left, right){
  n = nrow(X)
  out = c()
  
  for(i in 1:n){
    X_i = X[i,]
    
    if(all(left < X_i) & all(X_i < right)){
      out = c(out, i)
    }
  }
  
  return(out)
}

X1 = data[1:(n/n_groups),]
X2 = data[(n/n_groups + 1):n,]

index_top = which.min(out_MRS$post.null.nodewise)

tree_left = out_MRS$MAPtree.left
tree_right = out_MRS$MAPtree.right

children_IDs = out_MRS$children.IDs

left  = tree_left[,index_top]
right = tree_right[,index_top]

left_l = tree_left[,children_IDs[1,index_top]]
left_r = tree_right[,children_IDs[1,index_top]]

right_l = tree_left[,children_IDs[2,index_top]]
right_r = tree_right[,children_IDs[2,index_top]]

dim = which(left_l != right_l)
is_odd = dim %% 2 == 1

if(is_odd){
  another_dim = dim + 1
}else{
  dim = dim - 1
  another_dim = dim + 1
}

X1_included = X1[Indices_included(X1, left, right),]
X2_included = X2[Indices_included(X2, left, right),]

dim_name = paste("x", dim, sep = "")
another_dim_name = paste("x", another_dim, sep = "")

range_x = range(data[,dim], left[dim], right[dim])
range_y = range(data[,another_dim], left[another_dim], right[another_dim])

# plot the observed points
plot(data[,dim], data[,another_dim], col="azure3", pch=20, xlab = dim_name, ylab = another_dim_name, main = model_name,
     xlim = range_x, ylim = range_y)

points(X1_included[,dim], X1_included[,another_dim], pch=21, col="blue")
points(X2_included[,dim], X2_included[,another_dim], pch=24, col="red")


#Draw the window
x_st = c(left[dim],left[dim],right[dim],right[dim])
x_end = c(left[dim],right[dim],right[dim],left[dim])
y_st = c(left[another_dim],right[another_dim],right[another_dim],left[another_dim])
y_end = c(right[another_dim],right[another_dim],left[another_dim],left[another_dim])

segments(x_st, y_st, x_end, y_end, lwd = 2)

#Draw the boundary
if(is_odd){
  x_st = x_end = left_r[dim]
  y_st = left[another_dim]
  y_end = right[another_dim]
}else{
  y_st = y_end = left_r[another_dim]
  x_st = left[dim]
  x_end = right[dim]
}
segments(x_st, y_st, x_end, y_end, lwd = 2)

legend("topleft", 
       legend = c("Group 1","Group 2"),
       pch = c(21,24),
       col = c("blue", 'red'))

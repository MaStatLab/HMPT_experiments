Generate_data = function(model_name, d, n){

  data = matrix(NA, nrow = n, ncol = d)

 if(model_name == "2D_Uniforms"){

    w = c(1/3, 1/3, 1/3)
    b1_x = c(0.1, 0.45)
    b1_y = c(0.35, 0.9)

    b2_x = c(0.2, 0.8)
    b2_y = c(0.45, 0.5)

    b3_x = c(0.7, 0.9)
    b3_y = c(0.05, 0.6)

    for(i in 1:n){
      u = runif(1)
      if(u < w[1]){
        x_temp = b1_x[1] + runif(1) * (b1_x[2] - b1_x[1])
        y_temp = b1_y[1] + runif(1) * (b1_y[2] - b1_y[1])
      }else if(u < w[1]+w[2]){
        x_temp = b2_x[1] + runif(1) * (b2_x[2] - b2_x[1])
        y_temp = b2_y[1] + runif(1) * (b2_y[2] - b2_y[1])
      }else{
        x_temp = b3_x[1] + runif(1) * (b3_x[2] - b3_x[1])
        y_temp = b3_y[1] + runif(1) * (b3_y[2] - b3_y[1])
      }

      data[i,1] = x_temp
      data[i,2] = y_temp
    }

  }else if(model_name == "2D_Clusters"){ 

    w = c(1/10, 3/10, 3/10, 3/10)
    cum_w = cumsum(w)

    parameters_x = matrix(NA, ncol = 4, nrow = 4)
    parameters_x[1,] = c(1, 1)
    parameters_x[2,] = c(15,45)
    parameters_x[3,] = c(45,15)
    parameters_x[4,] = c(37.5,22.5)

    parameters_y = matrix(NA, ncol = 4, nrow = 4)
    parameters_y[1,] = c(1, 1)
    parameters_y[2,] = c(15,45)
    parameters_y[3,] = c(22.5,37.5)
    parameters_y[4,] = c(45,15)

    for(i in 1:n){
      u = runif(1)
      sampled = 0

      for(j in 1:4){
        if((u < cum_w[j]) & (sampled == 0)){
          data[i,1] = rbeta(1, parameters_x[j,1], parameters_x[j,2])
          data[i,2] = rbeta(1, parameters_y[j,1], parameters_y[j,2])

          sampled = 1
        }
      }
    }

    #High-dimensional models
  }else if(model_name == "2D_Smooth"){

    parameter_x = 10
    parameter_y = 20


    for(i in 1:n){
      data[i,1] = rbeta(1, parameter_x, parameter_y)
      data[i,2] = rbeta(1, parameter_x, parameter_y)
    }


  }else if(model_name == "Multi_beta"){
    
    for(j in 1:(d/2)){
      
      dims = c(2*(j-1)+1, 2*j)
      p = 0.25 + 0.7/j
      
      for(i in 1:n){
        if(runif(1) < p){
          data[i, dims] = rbeta(2, 0.25, 1)
        }else{
          data[i, dims] = rbeta(2, 50/j, 50/j)
        }
      }
    }
    
  }

  return(data)
}

True_density = function(model_name, d, grid_points){

  if(d == 1){
    num_grid_density = length(grid_points)
    true_densities = numeric(num_grid_density)

  }else if(d==2){
    num_grid_density = nrow(grid_points)
    true_densities = numeric(num_grid_density)
  }

  if(model_name == "2D_Uniforms"){

    w = c(1/3, 1/3, 1/3)
    b1_x = c(0.1, 0.45)
    b1_y = c(0.35, 0.9)

    b2_x = c(0.2, 0.8)
    b2_y = c(0.45, 0.5)

    b3_x = c(0.7, 0.9)
    b3_y = c(0.05, 0.6)

    for(i in 1:nrow(grid_points)){
      if((grid_points[i,1] > b1_x[1]) & (grid_points[i,1] < b1_x[2])){
        if((grid_points[i,2] > b1_y[1]) & (grid_points[i,2] < b1_y[2])){
          true_densities[i] = true_densities[i] + w[1] * (b1_x[2] - b1_x[1])^(-1) * (b1_y[2] - b1_y[1])^(-1)
        }
      }

      if((grid_points[i,1] > b2_x[1]) & (grid_points[i,1] < b2_x[2])){
        if((grid_points[i,2] > b2_y[1]) & (grid_points[i,2] < b2_y[2])){
          true_densities[i] = true_densities[i] + w[2] * (b2_x[2] - b2_x[1])^(-1) * (b2_y[2] - b2_y[1])^(-1)
        }
      }

      if((grid_points[i,1] > b3_x[1]) & (grid_points[i,1] < b3_x[2])){
        if((grid_points[i,2] > b3_y[1]) & (grid_points[i,2] < b3_y[2])){
          true_densities[i] = true_densities[i] +w[3] * (b3_x[2] - b3_x[1])^(-1) * (b3_y[2] - b3_y[1])^(-1)
        }
      }
    }

  }else if(model_name == "2D_Clusters"){

    w = c(1/10, 3/10, 3/10, 3/10)
    cum_w = cumsum(w)

    parameters_x = matrix(NA, ncol = 4, nrow = 4)
    parameters_x[1,] = c(1, 1)
    parameters_x[2,] = c(15,45)
    parameters_x[3,] = c(45,15)
    parameters_x[4,] = c(37.5,22.5)

    parameters_y = matrix(NA, ncol = 4, nrow = 4)
    parameters_y[1,] = c(1, 1)
    parameters_y[2,] = c(15,45)
    parameters_y[3,] = c(22.5,37.5)
    parameters_y[4,] = c(45,15)


    for(i in 1:nrow(grid_points)){
      for(j in 1:4){
        true_densities[i] = true_densities[i] + w[j] * dbeta(grid_points[i,1], parameters_x[j,1], parameters_x[j,2]) *
          dbeta(grid_points[i,2], parameters_y[j,1], parameters_y[j,2])
      }
    }

  }else if(model_name == "2D_Smooth"){

    parameter_x = 10
    parameter_y = 20

    for(i in 1:nrow(grid_points)){
      true_densities[i] = dbeta(grid_points[i,1], parameter_x, parameter_y) *
                           dbeta(grid_points[i,2], parameter_x, parameter_y)
    }

  }

  return(true_densities)
}


Generate_data_two_sample = function(model_name, d, groups, n, G, null_or_not){

  if(model_name == "Location_shift"){

    data = matrix(NA, nrow=n, ncol=d)

    w = c(1/3, 1/3, 1/3)
    cum_w = cumsum(w)

    shift = -0.5

    means1 = c(-2.5,1,2)
    means2 = c(1,-2,2.5)

    sd1 = 0.5
    sd2 = 0.7

    for(j in 1:(d/2)){
      for(i in 1:n){

        group_i = groups[i]
        u = runif(1)

        if(u < cum_w[1]){

          data[i,(j-1)*2+1] = rnorm(1, means1[1], sd1)
          data[i,j*2]       = rnorm(1, means2[1], sd2)

        }else if(u < cum_w[2]){

          data[i,(j-1)*2+1] = rnorm(1, means1[2], sd1)
          data[i,j*2]       = rnorm(1, means2[2], sd2)

        }else{

          if((group_i == 2) & (null_or_not == 0) & (j <= (d/10))){

            data[i,(j-1)*2+1] = rnorm(1, means1[3] + shift, sd1)
            data[i,j*2]       = rnorm(1, means2[3] + shift, sd2)

          }else{

            data[i,(j-1)*2+1] = rnorm(1, means1[3], sd1)
            data[i,j*2]       = rnorm(1, means2[3], sd2)

          }

        }
      }
    }
  }else if(model_name == "Dispersion_difference"){

    data = matrix(NA, nrow=n, ncol=d)

    w = c(1/3, 1/3, 1/3)
    cum_w = cumsum(w)

    means1 = c(-2.5,1,2)
    means2 = c(1,-2,2.5)

    sd1 = 0.5
    sd2 = 0.7

    sd_sub = 0.40

    for(j in 1:(d/2)){
      for(i in 1:n){

        group_i = groups[i]
        u = runif(1)

        if(u < cum_w[1]){

          data[i,(j-1)*2+1] = rnorm(1, means1[1], sd1)
          data[i,j*2]       = rnorm(1, means2[1], sd2)

        }else if(u < cum_w[2]){

          data[i,(j-1)*2+1] = rnorm(1, means1[2], sd1)
          data[i,j*2]       = rnorm(1, means2[2], sd2)

        }else{

          if((group_i == 2) & (null_or_not == 0) & (j == 1) & (j <= (d/10))){

            data[i,(j-1)*2+1] = rnorm(1, means1[3], sd1-sd_sub)
            data[i,j*2]       = rnorm(1, means2[3], sd2-sd_sub)

          }else{

            data[i,(j-1)*2+1] = rnorm(1, means1[3], sd1)
            data[i,j*2]       = rnorm(1, means2[3], sd2)

          }

        }
      }
    }
  }else if(model_name == "Correlation"){

    data = matrix(NA, nrow=n, ncol=d)

    means = c(0.0, 0.0)
    sd = c(1.0, 1.0)
    cor = 0.75

    Cov = matrix(NA, nrow=2, ncol=2)
    Cov[1,1] = sd[1]
    Cov[2,1] = sqrt(sd[1] * sd[2]) * cor
    Cov[1,2] = sqrt(sd[1] * sd[2]) * cor
    Cov[2,2] = sd[2]

    R = chol(Cov)

    t(R) %*% R

    for(j in 1:(d/2)){
      for(i in 1:n){

        group_i = groups[i]

        obs = numeric(2)

        obs[1] = rnorm(1, means[1], sd[1])
        obs[2] = rnorm(1, means[2], sd[2])

        if((group_i == 2) & (null_or_not == 0)  & (j <= (d/10))){
          obs = t(R) %*% matrix(obs, nrow=2,ncol=1)
        }

        data[i,(j-1)*2+1] = obs[1]
        data[i,j*2]       = obs[2]
      }
    }
  }

  return(data)
}

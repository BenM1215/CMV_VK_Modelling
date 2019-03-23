functions{                                                // ODEs defined here
  real[] ode(real time,
                  real[] y,
                  real[] theta,
                  real[] x_r,
                  int[] x_i) {
        
        //define variables;
        real u = y[1];
        real i = y[2];
        real v = y[3];
        
        //define thetas
        real beta = theta[1];
        real lambda = theta[2];
        real delta = theta[3];
        real rho = theta[4];
        real c = theta[5];

        //define equations
        real du_dt = -beta * u * v + lambda;            // uninfected
        real di_dt = beta * u * v - delta * i;          // infected
        real dv_dt = rho * i - c * v;                   // virus
        
        //return numerical values for each equation;
        return { du_dt, di_dt, dv_dt };
  }
}


// Read-in our data from R.
data {                          // Read-in data here
        int K;                  // Number of patients
        int N;                  // Number of observations
        int numODEs;            // Number of ODEs
        int s[K];               // Group sizes

        real Y0s[K];            // Y0 for each ODE
        real t0;                // t0
        
        real time[N];           // Vector - all timepoints
        real obs[N];            // Vector - all VL observations
      }


transformed data {              // Transform data here.
      real x_r[0];              // Empty - covariates (real)
      int x_i[0];               // Empty - covariates (integer)
}


parameters {                                    // Parameters to estimate go here

      // Define uniform parameter bounds
      real<lower=0,upper=1> beta;               // Portion of cells infected
      real<lower=0,upper=1000> lambda;          // Influx of new healthy cells
      real<lower=0,upper=1> delta;              // Loss of infected cells
      real<lower=0,upper=10> rho;               // Secretion of virus from infected cells
      real<lower=0,upper=10> c;                 // Loss of virus
      real<lower=0,upper=1000000> startingU;    // Starting number of uninfected cells

      real<lower=0,upper=10> sigma;             // Variance^2
}


transformed parameters {                        // Evaluate the ODEs and transform parameters here
      real theta[5];                            // Vector - thetas
      real Y0[3];                               // Vector - Y0 for each ODE
      real z[N, 3];                             // Stores solutions for model
      
      theta[1] = beta;
      theta[2] = lambda;
      theta[3] = delta;
      theta[4] = rho;
      theta[5] = c;

      Y0[1] = startingU;
      Y0[3] = Y0s[1];                           // Set starting VL for each Ps
      Y0[2] = Y0s[1]/rho;                       // Set starting Infected Cells for each Ps

      {                                         // Anything in here is seperate to transformed paramters
                                                
      int pos;                                  // Loops over flat vector of measurements
      pos = 1;
      for (k in 1:K) {
        real yhat[s[k], 3];                     // Multidimensional vector - predictions
        yhat = 
          integrate_ode_bdf(ode, Y0, t0, 
          //integrate_ode_rk45(ode, Y0, t0, 
          time[pos:pos+s[k]-1], theta, 
          x_r, x_i);
        
        z[pos:pos+s[k]-1] = yhat;
        pos = pos + s[k];
      }
      
      }
}


model {                                       // Define error model here
                                              
                                              //STRONG PRIOR SEEMS ESSENTIAL FOR DECENT RUN TIMES
      theta[1:5] ~ normal(0,2);               // Prior on thetas
      obs ~ lognormal(z[,3], sigma);
      
}


generated quantities{                       // For posterior predictive checks
      
      //gen_time
      
}




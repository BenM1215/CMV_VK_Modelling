functions{                                              // ODEs defined here
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
        real beta = theta[1];                           // Portion of cells infected
        real lambda = theta[2];                         // Influx of new healthy cells
        real delta = theta[3];                          // Loss of infected cells
        real rho = theta[4];                            // Secretion of virus from infected cells
        real c = theta[5];                              // Loss of virus

        //define equations
        real du_dt = -beta * u * v + lambda;            // uninfected
        real di_dt = beta * u * v - delta * i;          // infected
        real dv_dt = rho * i - c * v;                   // virus
        
        //return numerical values for each equation;
        return { du_dt, di_dt, dv_dt };
  }
}


// Read-in our data from R.
data {                            // Read-in data here
        int nSubjects;            // Number of patients
        int nObs;                 // Number of observations
        int nODEs;                // Number of ODEs
        int s[nSubjects];         // Group sizes
        int gen_time_N;           // Length of time vector for generating posterior predictive checks

        real Y0s[nSubjects];      // Y0 for each ODE
        real t0;                  // t0
        
        real time[nObs];           // Vector - all timepoints
        real obs[nObs];            // Vector - all VL observations
        
        real gen_time[gen_time_N]; //Vector - time points to integrate over for posterior predictives.
      }


transformed data {                 // Transform data here.
      real x_r[0];                 // Empty - covariates (real)
      int x_i[0];                  // Empty - covariates (integer)
}


parameters {                                          // Parameters to estimate go here

      // Define uniform parameter bounds
      //real<lower=0> theta[5];                       // Vector - thetas {beta, lambda, delta, rho, c}
      
      real<lower = 0> beta;                           // Portion of cells infected
      real<lower = 0> lambda;                         // Influx of new healthy cells
      real<lower = 0> delta;                          // Loss of infected cells
      real<lower = 0> rho;                            // Secretion of virus from infected cells
      real<lower = 0> c;                              // Loss of virus
      
      real<lower=0> U0[nSubjects];                    // Vector - Starting number of uninfected cells
      //real<lower=0> U0;                             // Vector - Starting number of uninfected cells
      real<lower=0> I0[nSubjects];                    // Vector - Starting number of infected cells

      //real<lower=0,upper=100000> U0;
      //real<lower=0,upper=100000> I0;
      
      real<lower=0> sigma;                            // Variance^2
      
}


transformed parameters {                              // Evaluate the ODEs and transform parameters here
      real z[nObs, 3];                                // Vector - solutions for ODEs
      
      real thetaHat[5];
      real theta[5];
      
      thetaHat[1] = beta;
      thetaHat[2] = lambda;
      thetaHat[3] = delta;
      thetaHat[4] = rho;
      thetaHat[5] = c;
      
      theta = thetaHat;                          //To be replaced with matrix multiplication for IIV
      
      {                                         // Anything in here is seperate to transformed paramters
                                                
      int pos;                                  // Loops over flat vector of measurements
      real Y0[3];                               // Vector - Y0 for each ODE
      pos = 1;
      for (subject in 1:nSubjects) {
        real yhat[s[subject], 3];                     // Multidimensional vector - predictions
        
        Y0[1] = U0[subject];                      // Set starting Infected Cells for each patient
        //Y0[1] = U0-I0[subject];                      // Set starting Infected Cells for each patient
        Y0[2] = I0[subject];                      // Set starting Infected Cells for each patient
        //Y0[1] = U0;
        //Y0[2] = I0;
        Y0[3] = Y0s[subject];                         // Set starting VL for each patient
        
        yhat = 
          integrate_ode_bdf(ode, Y0, t0, 
          time[pos:pos+s[subject]-1], theta, 
          x_r, x_i);
        
        z[pos:pos+s[subject]-1] = yhat;
        pos = pos + s[subject];
      }
      
      }
}


model {                                       // Define error model here
                                              //STRONG PRIOR SEEMS ESSENTIAL FOR DECENT RUN TIMES
      //theta[1:5] ~ normal(0,2);               // Prior on thetas
      
      beta ~ normal(0, 2);
      lambda ~ normal(0, 2);
      delta ~ normal(0, 2);
      rho ~ normal(0, 2);
      c ~ normal(0, 2);
      
      U0[1:nSubjects] ~ normal(0,1e6);
      //U0 ~ normal(0,1e6);
      I0[1:nSubjects] ~ normal(0,1e6);
      obs ~ lognormal(z[,3], sigma);
      
}


generated quantities{                       // For posterior predictive checks
      real gen_pred[gen_time_N, 3];
        
      real Y0[3];
      //for (subject in 1:nSubjects) {
      //  
      //  Y0[1] = U0[subject];
      //  Y0[2] = I0[subject];
      //  Y0[3] = Y0s[subject];
      //  
      //  gen_pred = integrate_ode_bdf(ode, Y0, t0, gen_time, theta, x_r, x_i);
        ////rng
        ////for each patient
        ////etas
        
        
      //}
      Y0[1] = U0[1];
      //Y0[1] = U0;
      //Y0[1] = U0-I0[1];
      Y0[2] = I0[1];
      Y0[3] = Y0s[1];
      gen_pred = integrate_ode_bdf(ode, Y0, t0, gen_time, theta, x_r, x_i);
}




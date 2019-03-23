
functions{                                                    // ODEs defined here
  real[] ode(real time,
                  real[] y,
                  real[] theta,
                  real[] x_r,
                  int[] x_i) {
        real dydt[3];

        dydt[1] = -theta[1] * y[1] * y[3] + theta[2];         // uninfected
        dydt[2] = theta[1] * y[1] * y[3] - theta[3] * y[2];   // infected
        dydt[3] = theta[4] * y[2] - theta[5] * y[3];          // virus
        return dydt;
  }
}


// Read-in our data from R.

data {                          // Read-in data here
        int numPs;              // Number of Ps
        int totNumObs;          // Total number of observations
        int numODEs;            // Number of ODEs

        real Y0s[numPs];        // Y0 for each ODE
        real t0;                // t0
        
        real time[totNumObs];   // Vector - All timepoints
        real obs[totNumObs];    // Vector - All VL observations
      }


transformed data {              // Transform read-in data here.
      real x_r[0];              // Empty - covariates (real)
      int x_i[0];               // Empty - covariates (integer)

}


parameters {                                    // Parameters to estimate go here
      real<lower=0,upper=1> beta;               // Portion of cells infected
      real<lower=0,upper=1000> lambda;          // Influx of new healthy cells
      real<lower=0,upper=1> delta;              // Loss of infected cells
      real<lower=0,upper=10> rho;               // Secretion of virus from infected cells
      real<lower=0,upper=10> c;                 // Loss of virus
      real<lower=0,upper=1000000> startingU;    // Starting number of uninfected cells

      real<lower=0,upper=10> sigma;             // Variance squared
}


transformed parameters {                        // Evaluate the ODEs and transform parameters here
      real theta[5];                            // Vector - thetas
      real Y0[3];                               // Vector - Y0 for each ODE


      theta[1] = beta;
      theta[2] = lambda;
      theta[3] = delta;
      theta[4] = rho;
      theta[5] = c;

      Y0[1] = startingU;
      Y0[3] = Y0s[1];                          // Set starting VL for each Ps
      Y0[2] = Y0s[1]/rho;                      // Set starting Infected Cells for each Ps

      
}


model {                                         // Define data error model here
      real yhat[totNumObs, 3];                  // Multidimensional vector - predictions
      
      yhat = integrate_ode_bdf(ode, Y0, t0, time, theta, x_r, x_i);
      obs ~ normal(yhat[,3], sigma);
      // print(yhat)
}


generated quantities{
      real ypred[totNumObs, 3];
      ypred = integrate_ode_bdf(ode, Y0, t0, time, theta, x_r, x_i);;
}


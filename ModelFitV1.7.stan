functions{                                              // ODE function defined here
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

        //define ODEs - basic Nowak model
        real du_dt = -(beta * u * v) + lambda;            // uninfected
        real di_dt = (beta * u * v) - (delta * i);          // infected
        real dv_dt = (rho * i) - (c * v);                   // virus
        
        //return numerical values for each equation;
        return { du_dt, di_dt, dv_dt };
  }
}

data {                              // Read-in data here
        int nSubjects;              // Number of subjects
        int nObs;                   // Number of observations
        int nODEs;                  // Number of ODEs
        int s[nSubjects];           // Group sizes for ragged array iteration
        int gen_time_N;             // Length of time vector for generating posterior predictive checks

        real Y0s[nSubjects];        // Y0 for each ODE
        real t0;                    // t0
        
        real time[nObs];            // Vector - all observed timepoints
        real obs[nObs];             // Vector - all observed viral loads
        
        real gen_time[gen_time_N];  //Vector - time points to integrate over for posterior predictive checks
      }


transformed data {                // Define transformed data here
      
      int nRandom = 5;            // Number of random variables/thetas
      real x_r[0];                // Empty - covariates for ODEs (real)
      int x_i[0];                 // Empty - covariates for ODEs (integer)
}


parameters {                                          // Parameters and bounds are defined here

      real<lower = 0> betaHat;                        // Population - Portion of cells infected
      real<lower = 0> lambdaHat;                      // Population - Influx of new healthy cells
      real<lower = 0> deltaHat;                       // Population - Loss of infected cells
      real<lower = 0> rhoHat;                         // Population Secretion of virus from infected cells
      real<lower = 0> cHat;                           // Population - Loss of virus
      
      real<lower=0> U0[nSubjects];                    // Vector - Starting number of uninfected cells
      real<lower=0> I0[nSubjects];                    // Vector - Starting number of infected cells
      
      real<lower=0> sigma;                            // Variance^2
      
      cholesky_factor_corr[nRandom] L;                // Cholesky Factor for correlation matrix
      vector[nRandom] omega;                          // Vector of length of random variables
      matrix[nRandom, nSubjects] eta;                 // Matrix of R.V. for each subject and theta
      
}


transformed parameters {                              // Evaluate the ODEs and transform parameters here
      real z[nObs, 3];                                // Vector - solutions for all ODEs for each subject
      
      vector[nRandom] thetaHat;                       // Vector - Population theta estimates
      real theta[nRandom];                            // ODE solver requires variable 'theta'
      matrix[nSubjects, nRandom] theta_matrix;        // Matrix of theta estimates for each subject
      
      real<lower = 0> beta[nSubjects];                // Vector - Individual beta estimates
      real<lower = 0> lambda[nSubjects];              // Vector - Individual lambda estimates
      real<lower = 0> delta[nSubjects];               // Vector - Individual delta estimates
      real<lower = 0> rho[nSubjects];                 // Vector - Individual rho estimates
      real<lower = 0> c[nSubjects];                   // Vector - Individual c estimates
      
      thetaHat[1] = betaHat;                          // Fill vector of population estimates (beta)
      thetaHat[2] = lambdaHat;                        // Fill vector of population estimates (lambda)
      thetaHat[3] = deltaHat;                         // Fill vector of population estimates (delta)
      thetaHat[4] = rhoHat;                           // Fill vector of population estimates (rho)
      thetaHat[5] = cHat;                             // Fill vector of population estimates (c)
      
      //print("lp = ", target());
      
      //print("thetaHat = ", thetaHat);
      
      //print("eta = ", eta);
      
      // Define a matrix of individual theta estimates - matrix of population estimates repeated by the number of subjects undergoing elementwise multiplication of the exponent of a diagonal matrix of cholesky factors * their respective etas
      theta_matrix = (rep_matrix(thetaHat, nSubjects) .*  
        exp(diag_pre_multiply(omega, L * eta)))';
        
        //print("theta_matrix = ", theta_matrix);

      { // Anything in here is seperate from transformed paramters
                                                
      int pos;                                        // Loops over flat vector of measurements
      real Y0[3];                                     // Vector - Y0 for each ODE
      pos = 1;
      for (subject in 1:nSubjects) {                  // Loop over each subject in ragged array
        real yhat[s[subject], 3];                     // Vector for ODE output
        
        beta[subject] = theta_matrix[subject, 1];     // Individual level beta
        lambda[subject] = theta_matrix[subject, 2];   // Individual level lambda
        delta[subject] = theta_matrix[subject, 3];    // Individual level delta
        rho[subject] = theta_matrix[subject, 4];      // Individual level rho
        c[subject] = theta_matrix[subject, 5];        // Individual level c
        
        theta[1] = theta_matrix[subject, 1];          // Pass theta 1 to ODE solver
        theta[2] = theta_matrix[subject, 2];          // Pass theta 2 to ODE solver
        theta[3] = theta_matrix[subject, 3];          // Pass theta 3 to ODE solver
        theta[4] = theta_matrix[subject, 4];          // Pass theta 4 to ODE solver
        theta[5] = theta_matrix[subject, 5];          // Pass theta 5 to ODE solver
        
        //print("theta_ode = ", theta);
        
        Y0[1] = U0[subject];                          // Starting Infected Cells for individual
        Y0[2] = I0[subject];                          // Starting Infected Cells for individual
        Y0[3] = Y0s[subject];                         // Starting VL for individual
        
        //print("Y0 = ", Y0);
        
        yhat =                                        // Compute ODE given parameters and time vector
          integrate_ode_bdf(ode, Y0, t0, 
          time[pos:pos+s[subject]-1], theta, 
          x_r, x_i);
        
        //print("yhat = ", yhat);
        z[pos:pos+s[subject]-1] = yhat;               // Store predictions for each ind. in vector
        pos = pos + s[subject];                       // Update position for loop
      }
      }
}

model {                                       // Define error model and priors here
      
      betaHat ~ lognormal(0, 2);                 // Prior on beta population estimate
      //print("betaHat = ", betaHat);
      lambdaHat ~ lognormal(0, 2);               // Prior on lambda population estimate
      //print("lambdaHat = ", lambdaHat);
      deltaHat ~ cauchy(0, 1);                // Prior on delta population estimate
      //print("deltaHat = ", deltaHat);
      rhoHat ~ lognormal(0, 2);                  // Prior on rho population estimate
      //print("rhoHat = ", rhoHat);
      cHat ~ lognormal(0, 2);                    // Prior on c population estimate
      //print("cHat = ", cHat);
      
      U0[1:nSubjects] ~ lognormal(0,1e6);        // Prior on each ind. starting uninfected cells
      //print("U0 = ", U0);
      I0[1:nSubjects] ~ lognormal(0,1e6);        // Prior on each ind. starting infected cells
      //print("I0 = ", I0);
      //U0[1:nSubjects] ~ cauchy(100000,1);        // Prior on each ind. starting uninfected cells
      //I0[1:nSubjects] ~ cauchy(100000,1);        // Prior on each ind. starting infected cells
      to_vector(eta) ~ normal(0, 0.1);          // Prior on etas for IIV
      //print("eta_model = ", eta);
      obs ~ lognormal(z[,3], sigma);          // Association between obs. data and predictions
      //print("obs = ", obs);
}

generated quantities{                       // Posterior predictive checks - runs every sample
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
      
      // Trial output for a single patient (patient 1)
      Y0[1] = U0[1];
      Y0[2] = I0[1];
      Y0[3] = Y0s[1];
      gen_pred = integrate_ode_bdf(ode, Y0, t0, gen_time, theta, x_r, x_i);
      //print("gen_pred", gen_pred);
}




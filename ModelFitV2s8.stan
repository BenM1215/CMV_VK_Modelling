functions{
  real[] ode(     real time,    // time
                  real[] y,     // system state {Uninfected, Infected, Virus}
                  real[] theta, // parameters
                  real[] x_r,   // data
                  int[] x_i) {
        
      real U = y[1];            // uninfected
      real I = y[2];            // infected
      real V = y[3];            // free virus
      
      real dydt[3];
      
      
      real lambda = theta[1];   // uninfected cell birth
      real d = theta[2];        // uninfected cell death
      real k = theta[3];        // viral infection rate
      real delta = theta[4];    // infected cell loss
      real rho = theta[5];      // viral production
      real c = theta[6];        // viral clearance rate
      
      dydt[1] = lambda - d * U - k * V * U;
      dydt[2] = k * V * U - delta * I;
      // (1-drug) * rho * I - cl * V
      dydt[3] = rho * I - c * V;
      
      return dydt;
  }
}

data {
  int<lower = 1> nSubjects;     // number of patients
  int<lower = 1> nObs;          // number of measurements
  int<lower = 1> nODEs;         // number of ODEs
  int<lower = 1> s[nSubjects];  // size of viral load vectors
  int<lower = 1> gen_time_N;    // size of time vectors for gen_quantities

  real Y0s[nSubjects];          // initial measured viral loads {y0}
  real t0;                      // initial time {t0}
  real time[nObs];              // measurement times {t1 -> tn}
  real obs[nObs];               // measured viral loads {y1 -> yn}
  
  real gen_time[gen_time_N];    // time for gen_quantities{t1 -> tn}
}


transformed data {
  real x_r[0];
  int x_i[0];
}


parameters {
  real <lower = 0> lambda_mu;
  real <lower = 0> d_mu;
  real <lower = 0> k_mu;
  real <lower = 0> delta_mu;
  real <lower = 0> rho_mu;
  real <lower = 0> c_mu;
  
  real <lower = 0> lambda_phi[nSubjects];
  real <lower = 0> d_phi[nSubjects];
  real <lower = 0> k_phi[nSubjects];
  real <lower = 0> delta_phi[nSubjects];
  real <lower = 0> rho_phi[nSubjects];
  real <lower = 0> c_phi[nSubjects];
  
  real <lower = 0> sigma;             // viral measurement error
  real <lower = 0> Y0_hat[nSubjects, nODEs];
}


transformed parameters {
    real z[nObs, nODEs];            // stores integrator output
    real pop_theta[6];              // fixed effects
    real ind_theta[6, nSubjects];   // fixed effects + mixed effects
    
    // *** Fixed parameters initiation ***
    real lambda_fix = 10000;
    real d_fix = 0.03;
    real k_fix = 2.6e-8;
    real delta_fix = 0.3;
    real rho_fix = 400;
    real c_fix = 5;
    
    // *** Define mu and phi ***
    real mu[6];
    real phi[6, nSubjects];
    
    // *** Population parameter fixing *** // uncomment parameters out to fix them
    //mu[1] = lambda_fix;
    //mu[2] = d_fix;
    //mu[3] = k_fix;
    //mu[4] = delta_fix;
    //mu[5] = rho_fix;
    //mu[6] = c_fix;
    // *** Individual parameter fixing *** // uncomment parameters out to fix them
    //phi[1, ] = rep_array(lambda_fix, nSubjects);
    //phi[2, ] = rep_array(d_fix, nSubjects);
    //phi[3, ] = rep_array(k_fix, nSubjects);
    //phi[4, ] = rep_array(delta_fix, nSubjects);
    //phi[5, ] = rep_array(rho_fix, nSubjects);
    //phi[6, ] = rep_array(c_fix, nSubjects);
    
    // *** Estimated population parameters *** // uncomment parameters to estimate them
    mu[1] = lambda_mu;
    mu[2] = d_mu;
    mu[3] = k_mu;
    mu[4] = delta_mu;
    mu[5] = rho_mu;
    mu[6] = c_mu;
    // *** Estimated individual parameters *** // uncomment parameters to estimate them
    phi[1, ] = lambda_phi;
    phi[2, ] = d_phi;
    phi[3, ] = k_phi;
    phi[4, ] = delta_phi;
    phi[5, ] = rho_phi;
    phi[6, ] = c_phi;
    
    
    { // Loop over patients
    int position; 
    position = 1;
    
    pop_theta = mu;                   // population level effects
    
    for (subject in 1:nSubjects) {    // loop over patients
      real yhat[s[subject], nODEs];   // store all ODE output for the patient
      //real Y0[nODEs];                 // store initial states in vector
      
      ind_theta[, subject] = phi[, subject];     // individual level effects

      //Y0[1] = 1000000;                          // fixed uninfected y0
      //Y0[2] = 0;                                // fixed infected y0
      //Y0[3] = Y0s[subject];                     // measured viral y0 - estimate?
        
      yhat = 
        integrate_ode_bdf(ode, Y0_hat[subject,], t0, 
        time[position:position+s[subject]-1], ind_theta[ ,subject], 
        x_r, x_i);
      
      z[position:position+s[subject]-1] = yhat; // update z vector
      
      position = position + s[subject];         // update position in loop

      }
    }
}


model { // comment parameters out if fixed
  // *** Literature priors ***
  lambda_mu ~ normal(10000, 2);     // lambda
  //d_mu ~ normal(0.03, 0.01);   // d
  k_mu ~ normal(2.6e-8, 0.01); // k
  //delta_mu ~ normal(0.3, 0.2);     // delta
  //rho_mu ~ normal(400, 0.2);     // rho
  //c_mu ~ normal(5, 0.2);       // c
  
  Y0_hat[,1] ~ normal(100000, 100);
  Y0_hat[,2] ~ normal(0, 2);
  Y0_hat[,3] ~ normal(1000, 10);
  
  // *** Weakly informative priors 1 ***
  //lambda_mu ~ normal(0, 10000);
  //d_mu ~ normal(0, 2);
  //k_mu ~ normal(0, 2);
  //delta_mu ~ normal(0, 2);
  //rho_mu ~ normal(0, 1000);
  //c_mu ~ normal(0, 10);

  // *** Weakly informative priors 2 ***
  //lambda_mu ~ normal(0, 100);       // lambda
  d_mu ~ normal(0, 100);       // d
  //k_mu ~ normal(0, 100);       // k
  delta_mu ~ normal(0, 100);       // delta
  rho_mu ~ normal(0, 100);       // rho
  c_mu ~ normal(0, 100);       // c
  
  // *** Non-informative priors ***
  //lambda_mu ~ normal(0, 1e6);
  //d_mu ~ normal(0, 1e6);
  //k_mu ~ normal(0, 1e6);
  //delta_mu ~ normal(0, 1e6);
  //rho_mu ~ normal(0, 1e6);
  //c_mu ~ normal(0, 1e6);
  
  // *** Center individual parameters on population means ***
  lambda_phi ~ normal(lambda_mu, 0.5);
  d_phi ~ normal(d_mu, 0.5);
  k_phi ~ normal(k_mu, 0.5);
  delta_phi ~ normal(delta_mu, 0.5);
  rho_phi ~ normal(rho_mu, 0.5);
  c_phi ~ normal(c_mu, 0.5);

  
  // *** Measurement error ***
  //sigma ~ cauchy(0, 1);
  sigma ~ lognormal(-1, 1);
  
  { // loop over patients
    int position;
    position = 1;
    
    for (subject in 1:nSubjects) {
      
      Y0s[subject] ~ lognormal(log(Y0_hat[subject, 3]), sigma);
      
      obs[position:position+s[subject]-1] ~ lognormal(log(z[position:position+s[subject]-1, nODEs]), sigma);
      
      position = position + s[subject];
    }
  }
}


generated quantities{
  real pop_pred[nObs, nODEs];
  real ind_pred[nObs, nODEs];
  
  {
  int position;
  position = 1;
  
  for (subject in 1:nSubjects) {
    real yhat_pop[s[subject], nODEs];
    real yhat_ind[s[subject], nODEs];
    //real Y0[nODEs];
    
    //Y0[1] = 1000000;
    //Y0[2] = 0;
    //Y0[3] = Y0s[subject];
    
    yhat_pop =
      integrate_ode_bdf(ode, Y0_hat[subject, ], t0, 
      time[position:position+s[subject]-1], pop_theta, 
      x_r, x_i);
    
    yhat_ind = 
      integrate_ode_bdf(ode, Y0_hat[subject, ], t0, 
      time[position:position+s[subject]-1], ind_theta[,subject], 
      x_r, x_i);
      
    pop_pred[position:position+s[subject]-1] = yhat_pop;
    ind_pred[position:position+s[subject]-1] = yhat_ind;
    
    position = position + s[subject];
  
    }
  }
}




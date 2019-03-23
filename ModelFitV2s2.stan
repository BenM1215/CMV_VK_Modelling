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
  // *** Thetas ***
  //real <lower = 0> lambda;      // uninfected cell birth
  //real <lower = 0> d;           // uninfected cell death
  real <lower = 0> k;             // viral infection rate
  //real <lower = 0> delta;       // infected cell loss
  //real <lower = 0> rho;         // viral production
  //real <lower = 0> c;           // viral clearance rate
  
  // *** Mixed effects ***
  //vector[nSubjects] eta_lambda;
  //vector[nSubjects] eta_d;
  vector[nSubjects] eta_k;
  //vector[nSubjects] eta_delta;
  //vector[nSubjects] eta_rho;
  //vector[nSubjects] eta_c;
  
  real <lower = 0> sigma[nSubjects];         // measurement error
  // initial population
  // means {mu}
}


transformed parameters {
    real z[nObs, nODEs];            // stores integrator output
    real pop_theta[6];              // fixed effects
    real ind_theta[nSubjects, 6];   // fixed effects + mixed effects
    
    // *** Fixed parameters initiation ***
    real lambda;
    real d;
    //real k;
    real delta;
    real rho;
    real c;
    
    { // Loop over patients
    int position; 
    position = 1;
    
    // *** Fixed parameters initiation ***
    lambda = 10000;
    d = 0.03;
    //k = 2.6e-8;
    delta = 0.3;
    rho = 400;
    c = 5;
    
    pop_theta[1] = lambda;
    pop_theta[2] = d;
    pop_theta[3] = k;
    pop_theta[4] = delta;
    pop_theta[5] = rho;
    pop_theta[6] = c;

    for (subject in 1:nSubjects) {    // loop over patients
      real yhat[s[subject], nODEs];   // store all ODE output for the patient
      real Y0[nODEs];                 // store initial states in vector
      
      // *** Fixed - no random effects ***
      ind_theta[subject, 1] = lambda;
      ind_theta[subject, 2] = d;
      ind_theta[subject, 3] = k;
      ind_theta[subject, 4] = delta;
      ind_theta[subject, 5] = rho;
      ind_theta[subject, 6] = c;
      
      // *** Estimated - w/ random effects ***
      //ind_theta[subject, 1] = lambda * exp(eta_lambda[subject]);
      //ind_theta[subject, 2] = d * exp(eta_d[subject]);
      //ind_theta[subject, 3] = k * exp(eta_k[subject]);
      //ind_theta[subject, 4] = delta * exp(eta_delta[subject]);
      //ind_theta[subject, 5] = rho * exp(eta_rho[subject]);
      //ind_theta[subject, 6] = c * exp(eta_c[subject]);
    
      Y0[1] = 1000000;                          // fixed uninfected y0
      Y0[2] = 0;                                // fixed infected y0
      Y0[3] = Y0s[subject];                     // measured viral y0 - estimate?
        
      yhat = 
        integrate_ode_bdf(ode, Y0, t0, 
        time[position:position+s[subject]-1], ind_theta[subject], 
        x_r, x_i);
      
      z[position:position+s[subject]-1] = yhat; // update z vector
      
      position = position + s[subject];         // update position in loop

      }
    }
}


model { 
  // *** Literature priors ***
  //lambda ~ normal(10000, 2);
  //d ~ normal(0.03, 0.01);
  //k ~ normal(2.6e-8, 0.01);
  //delta ~ normal(0.3, 0.2);
  //rho ~ normal(400, 0.2);
  //c ~ normal(5, 0.2);
  
  // *** Weakly informative priors 1 ***
  //lambda ~ normal(0, 10000);
  //d ~ normal(0, 2);
  k ~ normal(0, 2);
  //delta ~ normal(0, 2);
  //rho ~ normal(0, 1000);
  //c ~ normal(0, 10);
    
  // *** Weakly informative priors 2 ***
  //lambda ~ normal(0, 100);
  //d ~ normal(0, 100);
  //k ~ normal(0, 100);
  //delta ~ normal(0, 100);
  //rho ~ normal(0, 100);
  //c ~ normal(0, 100);
  
  // *** Non-informative priors ***
  //lambda ~ normal(0, 1e6);
  //d ~ normal(0, 1e6);
  //k ~ normal(0, 1e6);
  //delta ~ normal(0, 1e6);
  //rho ~ normal(0, 1e6);
  //c ~ normal(0, 1e6);
  
  // *** Measurement error ***
  sigma ~ cauchy(0, 1);
  
  // *** Eta priors ***
  //eta_lambda ~ normal(0, 0.1);
  //eta_d ~ normal(0, 0.1);
  //eta_k ~ normal(0, 0.1);
  //eta_delta ~ normal(0, 0.1);
  //eta_rho ~ normal(0, 0.1);
  //eta_c ~ normal(0, 0.1);
  
  { // loop over patients
    int position;
    position = 1;
    
    for (subject in 1:nSubjects) {
    
      obs[position:position+s[subject]-1] ~ lognormal(log(z[position:position+s[subject]-1, nODEs]), sigma[subject]);
      
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
    real Y0[nODEs];
    
    Y0[1] = 1000000;
    Y0[2] = 0;
    Y0[3] = Y0s[subject];
    
    yhat_pop =
      integrate_ode_bdf(ode, Y0, t0, 
      time[position:position+s[subject]-1], pop_theta, 
      x_r, x_i);
    
    yhat_ind = 
      integrate_ode_bdf(ode, Y0, t0, 
      time[position:position+s[subject]-1], ind_theta[subject], 
      x_r, x_i);
      
    pop_pred[position:position+s[subject]-1] = yhat_pop;
    ind_pred[position:position+s[subject]-1] = yhat_ind;
    
    position = position + s[subject];
  
    }
  }
}




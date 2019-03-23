functions{
  real[] ode(     real time,    // time
                  real[] y,     // system state {Uninfected, Infected, Virus}
                  real[] theta, // parameters
                  real[] x_r,   // data
                  int[] x_i) {
        
      real U = y[1];            // uninfected
      real I = y[2];            // infected
      real V = y[3];            // free virus
      int trigger;
      
      real dydt[3];
      
      matrix[2,(size(x_r)/2)] ode_data;
      
      real drugs[size(x_r)/2];          // data - drug
      real t[size(x_r)/2];          // data - time
      
      real lambda = theta[1];   // uninfected cell birth
      real d = theta[2];        // uninfected cell death
      real k = theta[3];        // viral infection rate
      real delta = theta[4];    // infected cell loss
      real rho = theta[5];      // viral production
      real c = theta[6];        // viral clearance rate
      real drugE = theta[7];        // drug effect
      
      ode_data = to_matrix(x_r, 2, (size(x_r)/2), 0);
      
      drugs = to_array_1d(ode_data[1,]);
      t = to_array_1d(ode_data[2,]);
      
      dydt[1] = lambda - d * U - k * V * U;
      dydt[2] = k * V * U - delta * I;
      
      trigger = 0;
      
      //With x_i
      //for (i in 1:min(x_i)){    //loop through all time points
      //  if (time == t[i]){  //time - matching index
      //    if (drugs[i] == 1){   //drug - matching index
      //      dydt[3] = (1 - drugE) * rho * I - c * V;
      //      trigger = 1;
      //    }
      //  }
      //}
      
      //Without x_i
      for (i in 1:size(t)){     //loop through all time points
        if (time == t[i]){        //time - matching index
          if (drugs[i] == 1){     //drug - matching index
            print(trigger);
            dydt[3] = (1 - drugE) * rho * I - c * V;
            trigger = 1;
          }
        }
      }
      
      if (trigger == 0){
        dydt[3] = rho * I - c * V;
      }
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
  real drugs[nObs];             // measured drug administrations
  
  real x_r[nObs*2];             // measured drug administrations
  
  real gen_time[gen_time_N];    // time for gen_quantities{t1 -> tn}
  
}


transformed data {
  //real x_r[2,nObs];
  //real x_r[nObs*2];
  int x_i[0];
  //int p;
  //int x_i[nSubjects];
    
    // forms array [drug_1, t_1, drug_2, t_2, ..., drug_n, t_n]  
    //p = 1;
    //for (subject in 1:nSubjects) {    // loop over patients
    
      //x_r[p:p+s[subject]-1] = drugs;

      //p = p + s[subject];         // update position in loop
      
      //x_r[p:p+s[subject]-1] = time;
      
      //p = p + s[subject];         // update position in loop

    //}

  
  
  //x_i = s;
}


parameters {
  real <lower = 0> lambda_mu;
  real <lower = 0> d_mu;
  real <lower = 0> k_mu;
  real <lower = 0> delta_mu;
  real <lower = 0> rho_mu;
  real <lower = 0> c_mu;
  real <lower = 0, upper = 1> drugE_mu;
  
  real <lower = 0> lambda_phi[nSubjects];
  real <lower = 0> d_phi[nSubjects];
  real <lower = 0> k_phi[nSubjects];
  real <lower = 0> delta_phi[nSubjects];
  real <lower = 0> rho_phi[nSubjects];
  real <lower = 0> c_phi[nSubjects];
  real <lower = 0, upper = 1> drugE_phi[nSubjects];
  
  real <lower = 0> U0_hat[nSubjects];
  real <lower = 0> I0_hat[nSubjects];
  real <lower = 0> V0_hat[nSubjects];
  
  real <lower = 0> sigma;             // viral measurement error
}


transformed parameters {
    real z[nObs, nODEs];            // stores integrator output
    real pop_theta[7];              // fixed effects
    real ind_theta[7, nSubjects];   // fixed effects + mixed effects
    real Y0_hat[3, nSubjects];      // Y0s
    
    // *** Fixed parameters initiation ***
    //real lambda_fix = 10000;
    //real d_fix = 0.03;
    //real k_fix = 2.6e-8;
    //real delta_fix = 0.3;
    //real rho_fix = 400;
    //real drugE_fix = 1;
    //real c_fix = 5;
    //real U0_fix = 100000;
    //real I0_fix = 0;
    //real V0_fix = 1000;
    
    // *** Define mu and phi ***
    real mu[7];
    real phi[7, nSubjects];
    
    // *** Population parameter fixing *** // uncomment parameters out to fix them
    //mu[1] = lambda_fix;
    //mu[2] = d_fix;
    //mu[3] = k_fix;
    //mu[4] = delta_fix;
    //mu[5] = rho_fix;
    //mu[6] = c_fix;
    //mu[7] = drugE_fix;
    
    // *** Individual parameter fixing *** // uncomment parameters out to fix them
    //phi[1, ] = rep_array(lambda_fix, nSubjects);
    //phi[2, ] = rep_array(d_fix, nSubjects);
    //phi[3, ] = rep_array(k_fix, nSubjects);
    //phi[4, ] = rep_array(delta_fix, nSubjects);
    //phi[5, ] = rep_array(rho_fix, nSubjects);
    //phi[6, ] = rep_array(c_fix, nSubjects);
    //phi[6, ] = rep_array(drugE_fix, nSubjects);
    
    // *** Y0s fixed *** // uncomment parameters out to fix them
    //Y0_hat[1, ] = rep_array(U0_fix, nSubjects);
    //Y0_hat[2, ] = rep_array(I0_fix, nSubjects);
    //Y0_hat[3, ] = rep_array(V0_fixed, nSubjects);
    
    // *** Estimated population parameters *** // uncomment parameters to estimate them
    mu[1] = lambda_mu;
    mu[2] = d_mu;
    mu[3] = k_mu;
    mu[4] = delta_mu;
    mu[5] = rho_mu;
    mu[6] = c_mu;
    mu[7] = drugE_mu;
    
    // *** Estimated individual parameters *** // uncomment parameters to estimate them
    phi[1, ] = lambda_phi;
    phi[2, ] = d_phi;
    phi[3, ] = k_phi;
    phi[4, ] = delta_phi;
    phi[5, ] = rho_phi;
    phi[6, ] = c_phi;
    phi[7, ] = drugE_phi;
    
    // *** Estimated Y0s *** //
    Y0_hat[1, ] = U0_hat;
    Y0_hat[2, ] = I0_hat;
    Y0_hat[3, ] = V0_hat;
    
    
    { // Loop over patients
    int position;
    int p2;
    position = 1;
    p2 = 1;
    
    pop_theta = mu;                   // population level effects
    
    for (subject in 1:nSubjects) {    // loop over patients
      real yhat[s[subject], nODEs];   // store all ODE output for the patient

      ind_theta[, subject] = phi[, subject];     // individual level effects
      
      yhat = 
        integrate_ode_bdf(ode, Y0_hat[,subject], t0, 
        time[position:position+s[subject]-1], ind_theta[ ,subject], 
        x_r[p2:(p2+(s[subject]*2)-1)], x_i);
        //x_r, x_i);
      
      z[position:position+s[subject]-1] = yhat; // update z vector
      
      position = position + s[subject];         // update position in loop
      p2 = p2 + s[subject] + s[subject];         // update p2 in loop

      }
    }
}


model { // comment parameters out if fixed
  // *** Literature priors ***
  //lambda_mu ~ normal(10000, 2);     // lambda
  //d_mu ~ normal(0.03, 0.01);   // d
  //k_mu ~ normal(2.6e-8, 0.01); // k
  //delta_mu ~ normal(0.3, 0.2);     // delta
  //rho_mu ~ normal(400, 0.2);     // rho
  //c_mu ~ normal(5, 0.2);       // c
  
  // *** Y0 priors *** //
  U0_hat ~ normal(100000, 100);
  I0_hat ~ normal(0, 2);
  V0_hat ~ normal(Y0s, 1);
  
  // *** Weakly informative priors 1 ***
  lambda_mu ~ normal(10000, 1000);
  d_mu ~ normal(0.03, 0.003);
  k_mu ~ normal(2.6e-8, 1e-8);
  delta_mu ~ normal(0.3, 0.03);
  rho_mu ~ normal(400, 40);
  c_mu ~ normal(5, 0.5);
  drugE_mu ~ normal(0.5, 0.2);

  // *** Weakly informative priors 2 ***
  //lambda_mu ~ normal(0, 100);       // lambda
  //d_mu ~ normal(0, 100);       // d
  //k_mu ~ normal(0, 100);       // k
  //delta_mu ~ normal(0, 100);       // delta
  //rho_mu ~ normal(0, 100);       // rho
  //c_mu ~ normal(0, 100);       // c
  
  // *** Non-informative priors ***
  //lambda_mu ~ normal(0, 1e6);
  //d_mu ~ normal(0, 1e6);
  //k_mu ~ normal(0, 1e6);
  //delta_mu ~ normal(0, 1e6);
  //rho_mu ~ normal(0, 1e6);
  //c_mu ~ normal(0, 1e6);
  
  // *** Center individual parameters on population means ***
  lambda_phi ~ normal(lambda_mu, 100);
  d_phi ~ normal(d_mu, 0.003);
  k_phi ~ normal(k_mu, 1e-8);
  delta_phi ~ normal(delta_mu, 0.03);
  rho_phi ~ normal(rho_mu, 40);
  c_phi ~ normal(c_mu, 0.5);
  drugE_phi ~ normal(drugE_mu, 2);

  
  // *** Measurement error ***
  //sigma ~ cauchy(0, 1);
  sigma ~ lognormal(-1, 1);
  
  { // loop over patients
    int position;
    position = 1;
    
    for (subject in 1:nSubjects) {
      
      Y0s[subject] ~ lognormal(log(Y0_hat[3, subject]), sigma);
      
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
  int p2;
  position = 1;
  p2 = 1;
  
  for (subject in 1:nSubjects) {
    real yhat_pop[s[subject], nODEs];
    real yhat_ind[s[subject], nODEs];
    //real x_r[(s[subject]*2)];        // stores x_r data
    
    //x_r[1:s[subject]] = drugs[position:position+s[subject]-1];
    //x_r[(s[subject]+1):(s[subject]*2)] = time[position:position+s[subject]-1];

    yhat_pop =
      integrate_ode_bdf(ode, Y0_hat[,subject], t0, 
      time[position:position+s[subject]-1], pop_theta, 
      //x_r[,position:position+s[subject]-1], x_i[subject]);
      x_r[p2:(p2+(s[subject]*2)-1)], x_i);    
      
    yhat_ind = 
      integrate_ode_bdf(ode, Y0_hat[,subject], t0, 
      time[position:position+s[subject]-1], ind_theta[,subject], 
      //x_r[,position:position+s[subject]-1], x_i[subject]);
      x_r[p2:(p2+(s[subject]*2)-1)], x_i);   
      
    pop_pred[position:position+s[subject]-1] = yhat_pop;
    ind_pred[position:position+s[subject]-1] = yhat_ind;
    
    position = position + s[subject];
    p2 = p2 + s[subject] + s[subject];         // update p2 in loop

    }
  }
}




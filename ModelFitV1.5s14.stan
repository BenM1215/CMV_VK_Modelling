//runs, 1 patient, all params. = 1
//MWE 
//  + prior per theta & theta initialisation
//  + looped observation fitting
//  + multiple IDs (same values and time points)
//  + Y0s
//  + all etas
//  + generated quantities

functions{
  real[] ode(     real time,
                  real[] y,
                  real[] theta,
                  real[] x_r,
                  int[] x_i) {
        
        real dydt[3];
        
        dydt[1] = -theta[1] * y[1] * y[3] + theta[2];         //uninfected
        dydt[2] = theta[1] * y[1] * y[3] - theta[3] * y[2];   //infected
        dydt[3] = theta[4] * y[2] - theta[5] * y[3];          //virus
        
        return dydt;
  }
}

data {
        int<lower = 1> nSubjects;
        int<lower = 1> nObs;
        int<lower = 1> nODEs;
        int<lower = 1> s[nSubjects];
        int<lower = 1> gen_time_N;

        real Y0s[nSubjects];
        real t0;
        real time[nObs];
        real obs[nObs];
        
        real gen_time[gen_time_N]; 
      }


transformed data {
      real x_r[0];
      int x_i[0];
}


parameters {

      real <lower = 0> beta;
      real <lower = 0> lambda;
      real <lower = 0> delta;
      real <lower = 0> rho;
      //real <lower = 0> c;
      
      real <lower = 0> mu_sigma;
      real <lower = 0> sigma[nSubjects];
      real <lower = 0> xi;
      
      vector[nSubjects] eta_beta;
      vector[nSubjects] eta_lambda;
      vector[nSubjects] eta_delta;
      vector[nSubjects] eta_rho;
      //vector[nSubjects] eta_c;
}


transformed parameters {
  
      real z[nObs, nODEs];
      real ind_theta[nSubjects, 5];
      real pop_theta[5];
      
      // Fixed parameters
      real c;
      
      {
      int position;
      position = 1;
      
      c = 5;
      
      pop_theta[1] = beta;
      pop_theta[2] = lambda;
      pop_theta[3] = delta;
      pop_theta[4] = rho;
      pop_theta[5] = c;
      
      for (subject in 1:nSubjects) {
        real yhat[s[subject], nODEs];
        real Y0[nODEs];
        
        ind_theta[subject, 1] = beta * exp(eta_beta[subject]);
        ind_theta[subject, 2] = lambda * exp(eta_lambda[subject]);
        ind_theta[subject, 3] = delta * exp(eta_delta[subject]);
        ind_theta[subject, 4] = rho * exp(eta_rho[subject]);
        //ind_theta[subject, 5] = c * exp(eta_c[subject]); //Forces model instability
        ind_theta[subject, 5] = c;
      
        Y0[1] = 100;
        Y0[2] = 0;
        Y0[3] = Y0s[subject];
          
        yhat = 
          integrate_ode_bdf(ode, Y0, t0, 
          time[position:position+s[subject]-1], ind_theta[subject], 
          x_r, x_i);
        
        z[position:position+s[subject]-1] = yhat;
        position = position + s[subject];

        }
      }
}


model { 
      beta ~ normal(1, 0.2);
      lambda ~ normal(1, 0.2);
      delta ~ normal(1, 0.2);
      rho ~ normal(1, 0.2);
      //c ~ normal(1, 0.2);
      mu_sigma ~ normal(0, 0.01);
      sigma ~ normal(mu_sigma, xi);
      xi ~ cauchy(0, 2.5);
      
      eta_beta ~ normal(0, 0.1);
      eta_lambda ~ normal(0, 0.1);
      eta_delta ~ normal(0, 0.1);
      eta_rho ~ normal(0, 0.1);
      //eta_c ~ normal(0, 0.1);
      
      {
      int position;
      position = 1;
      
      for (subject in 1:nSubjects) {
      
        obs[position:position+s[subject]-1] ~ normal(z[position:position+s[subject]-1, nODEs], sigma[subject]);
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
        
        Y0[1] = 100;
        Y0[2] = 0;
        Y0[3] = Y0s[subject];;
        
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




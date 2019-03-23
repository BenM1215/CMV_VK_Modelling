//runs, 1 patient, all params. = 1
//MWE 
//  + prior per theta & theta initialisation
//  + looped observation fitting
//  + multiple IDs (same values and time points)
//  + Y0s
//  + all etas
//  + generated quantities
//  + real data

functions{
  real[] ode(     real time,
                  real[] y,
                  real[] theta,
                  real[] x_r,
                  int[] x_i) {
        
        real dydt[3];
        
        //dydt[1] = -theta[1] * y[1] * y[3] + theta[2];         //uninfected
        //dydt[2] = theta[1] * y[1] * y[3] - theta[3] * y[2];   //infected
        //dydt[3] = theta[4] * y[2] - theta[5] * y[3];          //virus
        
        // lambda - d * U - k * V * U
        dydt[1] = theta[1] - theta[2] * y[1] - theta[3] * y[3] * y[1];   //uninfected
        //dydt[1] = 10000 - theta[2] * y[1] - 2.6e-8 * y[3] * y[1];   //uninfected
        // k * V * U - delta * I
        dydt[2] = theta[3] * y[1] * y[3] - theta[4] * y[2];   //infected
        //dydt[2] = 2.6e-8 * y[1] * y[3] - 0.3 * y[2];   //infected
        // (1-drug) * rho * I - cl * V
        // rho * I - cl * V
        dydt[3] = theta[5] * y[2] - theta[6] * y[3];          //virus
        //dydt[3] = 400 * y[2] - 5 * y[3];          //virus
        
        
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


      //real <lower = 0> lambda;
      //real <lower = 0> d;
      real <lower = 0> k;
      //real <lower = 0> delta;
      //real <lower = 0> rho;
      //real <lower = 0> c;
      
      //real <lower = 0> startingI;

      real <lower = 0> sigma;
      
      //vector[nSubjects] eta_lambda;
      //vector[nSubjects] eta_d;
      vector[nSubjects] eta_k;
      //vector[nSubjects] eta_delta;
      //vector[nSubjects] eta_rho;
      //vector[nSubjects] eta_c;
}


transformed parameters {
  
      real z[nObs, nODEs];
      real ind_theta[nSubjects, 6];
      real pop_theta[6];
      
      // FIXED
      real lambda;
      real d;
      //real k;
      real delta;
      real rho;
      real c;
      
      {
      int position;
      position = 1;
      
      //FIXED
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

      for (subject in 1:nSubjects) {
        real yhat[s[subject], nODEs];
        real Y0[nODEs];
        
        //Fixed - no random effects
        ind_theta[subject, 1] = lambda;
        ind_theta[subject, 2] = d;
        //ind_theta[subject, 3] = k;
        ind_theta[subject, 4] = delta;
        ind_theta[subject, 5] = rho;
        ind_theta[subject, 6] = c;
        
        //Estimated - w/ random effects
        //ind_theta[subject, 1] = lambda * exp(eta_lambda[subject]);
        //ind_theta[subject, 2] = d * exp(eta_d[subject]);
        ind_theta[subject, 3] = k * exp(eta_k[subject]);
        //ind_theta[subject, 4] = delta * exp(eta_delta[subject]);
        //ind_theta[subject, 5] = rho * exp(eta_rho[subject]);
        //ind_theta[subject, 6] = c * exp(eta_c[subject]);
        
      
        Y0[1] = 1000000;
        //Y0[2] = startingI;
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
      
      //proper
      //lambda ~ normal(10000, 2);
      //d ~ normal(0.03, 0.01);
      k ~ normal(2.6e-8, 0.01);
      //delta ~ normal(0.3, 0.2);
      //rho ~ normal(400, 0.2);
      //c ~ normal(5, 0.2);
      
      //middle ground - nope
      //lambda ~ normal(10000, 100);
      //d ~ normal(0.03, 10);
      //k ~ normal(2.6e-8, 10);
      //delta ~ normal(0.3, 10);
      //rho ~ normal(400, 10);
      //c ~ normal(5, 20);
      
      //non-informative - less
      //lambda ~ normal(0, 10000);
      //d ~ normal(0, 2);
      //k ~ normal(0, 2);
      //delta ~ normal(0, 2);
      //rho ~ normal(0, 1000);
      //c ~ normal(0, 10);
        
      //startingI ~ normal(0, 1e6);
      
      //works
      //lambda ~ normal(0, 100);
      //d ~ normal(0, 100);
      //k ~ normal(0, 100);
      //delta ~ normal(0, 100);
      //rho ~ normal(0, 100);
      //c ~ normal(0, 100);

      
      //non-informative
      //lambda ~ normal(0, 1e6);
      //d ~ normal(0, 1e6);
      //k ~ normal(0, 1e6);
      //delta ~ normal(0, 1e6);
      //rho ~ normal(0, 1e6);
      //c ~ normal(0, 1e6);
      
      //mu_sigma ~ normal(0, 0.01);
      //sigma ~ normal(mu_sigma, xi);
      sigma ~ cauchy(0, 1);
      //xi ~ cauchy(0, 2.5);
      
      //eta_lambda ~ normal(0, 0.1);
      //eta_d ~ normal(0, 0.1);
      //eta_k ~ normal(0, 0.1);
      //eta_delta ~ normal(0, 0.1);
      //eta_rho ~ normal(0, 0.1);
      //eta_c ~ normal(0, 0.1);
      
      {
      int position;
      position = 1;
      
      for (subject in 1:nSubjects) {
      
        obs[position:position+s[subject]-1] ~ normal(z[position:position+s[subject]-1, nODEs], sigma);
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
        //Y0[2] = startingI;
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




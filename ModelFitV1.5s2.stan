//runs, 1 patient, all params. = 1
//MWE 
//  + prior per theta & theta initialisation

functions{
  real[] ode(     real time,
                  real[] y,
                  real[] theta,
                  real[] x_r,
                  int[] x_i) {
        
        real dydt[3];
        
        dydt[1] = -theta[1] * y[1] * y[3] + theta[2];
        dydt[2] = theta[1] * y[1] * y[3] - theta[3] * y[2];
        dydt[3] = theta[4] * y[2] - theta[5] * y[3];
        
        return dydt;
  }
}

data {
        int<lower = 1> nSubjects;
        int<lower = 1> nObs;
        int<lower = 1> nODEs;
        int<lower = 1> s;
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
      real <lower = 0> c;
      real <lower = 0> sigma;
}


transformed parameters {
  
      real yhat[s, 3];
      real Y0[nODEs];
      real theta[5];
      
      theta[1] = beta;
      theta[2] = lambda;
      theta[3] = delta;
      theta[4] = rho;
      theta[5] = c;

      Y0[1] = 100;
      Y0[2] = 0;
      Y0[3] = 100;
        
      yhat = 
        integrate_ode_rk45(ode, Y0, t0, 
        time, theta, 
        x_r, x_i);
}


model { 
                                              
      beta ~ normal(0.5, 0.2);
      lambda ~ normal(0.5, 0.2);
      delta ~ normal(0.5, 0.2);
      rho ~ normal(0.5, 0.2);
      c ~ normal(0.5, 0.2);
      sigma ~ cauchy(0, 1);
      
      obs ~ normal(yhat[,3], sigma);
}


generated quantities{

}




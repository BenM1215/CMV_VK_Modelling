//runs, 1 patient, all params. = 1
//MWE

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

      real <lower = 0> theta[5];
      real <lower = 0> sigma;
}


transformed parameters {
  
      real yhat[s, 3];
      real Y0[nODEs];

      Y0[1] = 100;
      Y0[2] = 0;
      Y0[3] = 100;
        
      yhat = 
        integrate_ode_rk45(ode, Y0, t0, 
        time, theta, 
        x_r, x_i);
}


model { 
                                              
      theta ~ normal(0.5, 0.2);
      sigma ~ cauchy(0, 1);
      
      obs ~ normal(yhat[,3], sigma);
}


generated quantities{

}




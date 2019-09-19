data{
//Sample sizes
 int<lower=1> N; // total number of obs (long format)
 int<lower=1> nstud; // # students
 int<lower=1> nsec; // # sections (problems nested in sections)

// indices
 int<lower=1,upper=nstud> studentM[N]; //for each obs, which student?
 int<lower=1,upper=nsec> section[N]; // for each obs, which section?

// data data
 int<lower=0,upper=1> hint[N]; // for each obs, was a hint requested?

}
parameters{

 vector[nstud] studEff; // ability parameter
 real<lower=0,upper=1> p1; // mixing proportion
 //ordered[2] mu;             // locations of mixture components
 real<upper=0> mu0; //location of mixture component 0
 vector<lower=0>[2] sigma;  // scales of mixture components

 real secEff[nsec]; // section difficulty parameter
 //real<lower=0> sigSec; // scale of difficulty parameters

}
transformed parameters{
 real<lower=0> mu1; // location of mixture component 1
 mu1= -(1-p1)/p1*mu0;
}
model{
 // logs of mixing proportions to save time
 real log_p1=log(p1);
 real log_p0=log(1-p1);

// linear predictor for logit model
 real linPred[N];

// priors
 sigma ~ lognormal(0, 2);
 mu0 ~ normal(0, 2);
// sigSec ~ lognormal(0,2);
 p1~ beta(2,2); // I'd rather the groups have similar sizes

// model for section difficulty parameters
 secEff~normal(0,3);//sigSec);

// Rasch model for hint requests
 for(i in 1:N)
  linPred[i]= studEff[studentM[i]]-secEff[section[i]];

 hint~bernoulli_logit(linPred);


// mixture model for ability parameters
 for(n in 1:nstud)
  target += log_sum_exp(log_p0+normal_lpdf(studEff[n]|mu0,sigma[1]),
			log_p1+normal_lpdf(studEff[n]|mu1,sigma[2]));

}
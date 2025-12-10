data {
	int<lower=0> N; // sample size
	vector[N] X; // intervals
	vector[N] Y; // rates
	int<lower=0> pop[N]; // population indexes
	int<lower=0> J; // number of populations
}

parameters {
	vector[J] alpha; // pop intercepts
	vector[J] beta; // pop slope
	real<lower=0> sigma; // residual error
	real mu[2]; // mean for alphas and betas
	real<lower=0> sd[2]; // sd for alphas and betas
}

model{
	// linear model and likelihood
	for(i in 1:N){
		Y[i] ~ normal(alpha[pop[i]] + beta[pop[i]] * X[i], sigma);
	}
	// hierarchical priors
	alpha ~ normal(mu[1], sd[1]);
	beta ~ normal(mu[2], sd[2]);
	// hyperpriors
	mu ~ normal(0,10);
	sd ~ normal(0,10);

	// other priors
	sigma ~ normal(0,10);
}


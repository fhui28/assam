data {
    int<lower=1> num_sites;
    int<lower=1> num_species;
    int<lower=1> num_X;
    int<lower=1> num_archetypes;
    
    array[num_sites, num_species] int<lower=0, upper=1> y; 
    matrix[num_sites, num_X] X;
    vector<lower=0>[num_archetypes] alpha_dirichlet;
}

parameters {
    //positive_ordered[num_archetypes] mixprop_raw; ; // Forced to be: mixprop_raw[1] < mixprop_raw[2] < ...
    simplex[num_archetypes] mixprop;             // Mixing proportions (sums to 1)
    matrix[num_archetypes, num_X] beta;     // Archetype-level coefficients
    vector[num_species] alpha;              // Species-level intercepts
}

model {
    // Priors
    mixprop ~ dirichlet(alpha_dirichlet);
    alpha ~ normal(0, 5);                   // Weakly informative
    for (g in 1:num_archetypes) {
        beta[g] ~ normal(0, 5);
    }
    
    // Likelihood: Marginalizing out the discrete latent labels z[j]
    for (j in 1:num_species) {
        vector[num_archetypes] lps; // log-probabilities for each archetype
        
        for (g in 1:num_archetypes) {
            // Start with the log-prior probability of being in archetype g
            lps[g] = log(mixprop[g]);
            
            // Add the log-likelihood of the data for species j given archetype g
            // logit(p) = alpha + X * beta
            lps[g] += bernoulli_logit_lpmf(y[, j] | alpha[j] + X * beta[g]');
    }
    
    // Log-Sum-Exp trick to marginalize: log(sum(exp(lps)))
    target += log_sum_exp(lps);
  }
}

generated quantities {
  // Posterior membership probabilities (Soft clustering)
  matrix[num_species, num_archetypes] z_prob;
  
  for (j in 1:num_species) {
    vector[num_archetypes] log_probs;
    for (g in 1:num_archetypes) {
      log_probs[g] = log(mixprop[g]) + bernoulli_logit_lpmf(y[, j] | alpha[j] + X * beta[g]');
        }
    z_prob[j, ] = to_row_vector(softmax(log_probs));
    }
}


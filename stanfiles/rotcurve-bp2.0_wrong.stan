/*
*
*
*
*
*
*
*
*
*
*
*/

functions{
  array[] vector predicted_proper_motions(vector plx, array[] vector r, vector sunpos, 
      real h_param, real p_param, array[] vector p, array[] vector q, real Av) {
    
    array[size(plx)] vector[2] predicted_pm;
    vector[3] starpos;
    real Rstar;
    real vphistar;
    real phi;
    vector[3] vdiff;

    for (i in 1:size(plx)) {
      starpos = (1000.0/plx[i])*r[i] + sunpos;
      Rstar = sqrt(starpos[1]^2+starpos[2]^2);
      vphistar = Rstar/h_param * (1 + (Rstar/h_param)^2) ^ ((p_param-2)/4);
      phi = atan2(starpos[2], starpos[1]);
      vdiff = [-vphistar*sin(phi), vphistar*cos(phi), 0.0]'; // In andere stond ook nog -vsun, is dat hier nog nodig? (Dit gaat met name om peculiar motion en dit was in het vorige model ook een parameter.
      predicted_pm[i][1] = dot_product(p[i], vdiff) * plx[i] / Av;
      predicted_pm[i][2] = dot_product(q[i], vdiff) * plx[i] / Av;
    }

    return predicted_pm;
  }
}

data {
  int<lower=0> N;
  vector[N] galon;
  vector[N] galat;
  vector[N] pml_obs;
  vector[N] pml_obs_unc;
  vector[N] pmb_obs;
  vector[N] pmb_obs_unc;
  vector[N] pml_pmb_corr;
  vector[N] plx_obs;
  real Rsun;
  real Zsun;
}

transformed data {
  real Ysun = 0.0;                // Sun galactocentric Cartesian y-coordinate (0 by definition)
  // parameters for priors
  real h_param_prior_alpha = 1.1; // Breder maken oud: 2.0
  real h_param_prior_beta = 10; // Breder maken oud: 1.0
  real p_param_mean = -0.5;
  real p_param_sigma = 1.5;

  real auInMeter = 149597870700.0;
  real julianYearSeconds = 365.25 * 86400.0;
  real auKmYearPerSec = auInMeter/(julianYearSeconds*1000.0);

  array[N] cov_matrix[2] cov_pm;             // Covariance matrix of the proper motions (auxiliary variable only)
  array[N] vector[2] pm_obs;                 // Observed proper motions
  array[N] vector[3] pvec;
  array[N] vector[3] qvec;
  array[N] vector[3] rvec;

  for (n in 1:N) {
    cov_pm[n][1,1] = pml_obs_unc[n]^2;
    cov_pm[n][2,2] = pmb_obs_unc[n]^2;
    cov_pm[n][1,2] = pml_obs_unc[n]*pmb_obs_unc[n]*pml_pmb_corr[n];
    cov_pm[n][2,1] = cov_pm[n][1,2];

    pm_obs[n][1] = pml_obs[n];
    pm_obs[n][2] = pmb_obs[n];

    pvec[n][1] = -sin(galon[n]);
    pvec[n][2] = cos(galon[n]);
    pvec[n][3] = 0.0;
    
    qvec[n][1] = -sin(galat[n])*cos(galon[n]);
    qvec[n][2] = -sin(galat[n])*sin(galon[n]);
    qvec[n][3] = cos(galat[n]);
    
    rvec[n][1] = cos(galat[n])*cos(galon[n]);
    rvec[n][2] = cos(galat[n])*sin(galon[n]);
    rvec[n][3] = sin(galat[n]);
  }
}

parameters {
  real <lower=0> h_param;
  real p_param;
}

transformed parameters {
  array[N] vector[2] model_pm;
  model_pm = predicted_proper_motions(plx_obs, rvec, [-Rsun, Ysun, Zsun]', h_param, p_param, pvec, qvec, auKmYearPerSec);
}

model {
  h_param ~ gamma(h_param_prior_alpha, h_param_prior_beta);
  p_param ~ normal(p_param_mean, p_param_sigma);
  for (i in 1:N) {
    pm_obs[i] ~ multi_normal(model_pm[i], cov_pm[i]); // Wat doet multi normal precies en is cov_pm nu wel normal?
// Klopt dit model voor pm_obs nu?
  }
}

generated quantities {
  vector[N] pred_pml;
  vector[N] pred_pmb;
  vector[2] pred_pm;

  for (i in 1:N) {
    pred_pm = multi_normal_rng(model_pm[i], cov_pm[i]); // Klopt dit nog? Kijken naar cov_pm (dcov in ..._old.stan) het is in old een combi met de vdisp parameter en hier is het alleen een observable (en miss niet eens een normal function!)
    pred_pml[i] = pred_pm[1];
    pred_pmb[i] = pred_pm[2];
  }
}


// Zelf wiskundige achtergrond van het model even afleiden voor het goede inzicht. Nogmaals paper BP erbij houden.
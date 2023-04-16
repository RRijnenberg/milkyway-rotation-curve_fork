/*
* Model with a function for the velocity dispersion of giants.
* This model will only be accurate for the giant subsample. When another subsample is fitted to this stan model the wrong velocity dispersion behaviour
* will be modelled.
*
* The v_0 parameter is removed from this model, as from previous fits we can conclude that there is no decent constraint for this parameter.
*
*
*
*
*
*/

functions{
  array[] vector predicted_proper_motions(vector plx, array[] vector r, vector sunpos, 
      real h_param, real p_param, array[] vector p, array[] vector q, real Av, vector Vsun_pec, real v0) {
    
    array[size(plx)] vector[2] predicted_pm;
    vector[3] starpos;
    real Rstar;
    real vphistar;
    real phi;
    vector[3] vdiff;

    // Toegevoegd:
    real Rsun = sqrt(sunpos[1]^2+sunpos[2]^2)/1000.0; // In kpc
    real Vcirc_sun = v0 * Rsun/h_param * (1 + (Rsun/h_param)^2) ^ ((p_param-2)/4);
    vector[3] Vsun = [0.0, Vcirc_sun, 0.0]' + Vsun_pec;

    for (i in 1:size(plx)) {
      starpos = (1.0/plx[i])*r[i] + sunpos/1000.0;  // In kpc
      Rstar = sqrt(starpos[1]^2+starpos[2]^2);
      vphistar = v0*Rstar/h_param * (1 + (Rstar/h_param)^2) ^ ((p_param-2)/4);
      phi = atan2(starpos[2], starpos[1]);
      vdiff = [-vphistar*sin(phi), vphistar*cos(phi), 0.0]' - Vsun;
      predicted_pm[i][1] = dot_product(p[i], vdiff) * plx[i] / Av;
      predicted_pm[i][2] = dot_product(q[i], vdiff) * plx[i] / Av;
    }

    return predicted_pm;
  }


  // real xy_disp_giants(real plx, vector r, vector sunpos, real amplitude, real R_scale) {
	// /* sigma(R)
  // * Part one: R <= constant --> return constant
	// * Part two: R > constant --> return value from function
  // */

  // real disp_mean;
  // vector[3] starpos;
  // real Rstar;

  // starpos = (1.0/plx)*r + sunpos/1000.0;  // In kpc
  // Rstar = sqrt(starpos[1]^2+starpos[2]^2);

  // if (Rstar < 2.5) {
  //   disp_mean = amplitude;
  // }
  // else {
  //   disp_mean = amplitude * exp(-(Rstar-2.5)/R_scale);
  // }
  // return disp_mean;
  // }


  // real z_disp_giants(real plx, vector r, vector sunpos, real amplitude_z, real R_scale) {
  
  // real disp_mean;
  // vector[3] starpos;
  // real Rstar;

  // starpos = (1.0/plx)*r + sunpos/1000.0;  // In kpc
  // Rstar = sqrt(starpos[1]^2+starpos[2]^2);
  // disp_mean = amplitude_z * exp(-Rstar/R_scale);
  // return disp_mean;
  // }

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
  real Ysun = 0.0;                     // Sun galactocentric Cartesian y-coordinate (0 by definition)
  // parameters for priors
  real h_param_prior_alpha = 1.1;      // In kpc
  real h_param_prior_beta = 10;        // In kpc
  real p_param_mean = -0.5; 
  real p_param_sigma = 1.5;
  real Vsun_pec_x_prior_mean = 11.0;
  real Vsun_pec_y_prior_mean = 12.0;
  real Vsun_pec_z_prior_mean = 7.0;
  real Vsun_pec_x_prior_sigma = 20.0;
  real Vsun_pec_y_prior_sigma = 20.0;
  real Vsun_pec_z_prior_sigma = 20.0;
  real v0_prior_mean = 234.0;      // In km/s
  real v0_prior_sigma = 20;
  real vdisp_prior_alpha = 2.0; // Used for both xy and z dispersion
  real vdisp_prior_beta = 0.1;  //

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
  real <lower=0> h_param;      // Scalelength of the rotational velocity function introduced by B&P
  real p_param;                // Model parameter allowing study of different versions of the introduced model type by B&P
  real Vsun_pec_x;             // Peculiar velocity of Sun in Galactocentric Cartesian X
  real Vsun_pec_y;             // Peculiar velocity of Sun in Galactocentric Cartesian Y
  real Vsun_pec_z;             // Peculiar velocity of Sun in Galactocentric Cartesian Z
  real vdispxy;
  real vdispz;
  real v0;

}

transformed parameters {
  array[N] vector[2] model_pm;
  cov_matrix[3] scov;          // Model covariance matrix for velocity dispersions
  array[N] cov_matrix[2] dcov;       // Total covariance matrix (observational uncertainties plus velocity dispersion)

  scov[1,1] = vdispxy^2;
  scov[2,2] = vdispxy^2;
  scov[3,3] = vdispz^2;
  scov[1,2] = 0.0;
  scov[1,3] = 0.0;
  scov[2,3] = 0.0;
  scov[2,1] = scov[1,2];
  scov[3,1] = scov[1,3];
  scov[3,2] = scov[2,3];

  for (n in 1:N) {
    dcov[n][1,1] = dot_product(pvec[n], pvec[n]); 
    dcov[n][2,1] = dot_product(qvec[n], pvec[n]); 
    dcov[n][1,2] = dot_product(pvec[n], qvec[n]); 
    dcov[n][2,2] = dot_product(qvec[n], qvec[n]);
    dcov[n] = cov_pm[n] + (plx_obs[n]/auKmYearPerSec)^2 * dcov[n];
  }

  model_pm = predicted_proper_motions(plx_obs, rvec, [-Rsun, Ysun, Zsun]', h_param, p_param, pvec, qvec, auKmYearPerSec, [Vsun_pec_x, Vsun_pec_y, Vsun_pec_z]', v0);
}

model {
  h_param ~ gamma(h_param_prior_alpha, h_param_prior_beta);
  p_param ~ normal(p_param_mean, p_param_sigma);
  Vsun_pec_x ~ normal(Vsun_pec_x_prior_mean, Vsun_pec_x_prior_sigma);
  Vsun_pec_y ~ normal(Vsun_pec_y_prior_mean, Vsun_pec_y_prior_sigma);
  Vsun_pec_z ~ normal(Vsun_pec_z_prior_mean, Vsun_pec_z_prior_sigma);
  vdispxy ~ gamma(vdisp_prior_alpha, vdisp_prior_beta);
  vdispz ~ gamma(vdisp_prior_alpha, vdisp_prior_beta);
  v0 ~ normal(v0_prior_mean, v0_prior_sigma);

  for (i in 1:N) {
    pm_obs[i] ~ multi_normal(model_pm[i], dcov[i]);
  }
}

generated quantities {
  vector[N] pred_pml;
  vector[N] pred_pmb;
  vector[2] pred_pm;

  for (i in 1:N) {
    pred_pm = multi_normal_rng(model_pm[i], dcov[i]);
    pred_pml[i] = pred_pm[1];
    pred_pmb[i] = pred_pm[2];
  }
}
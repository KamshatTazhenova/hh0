functions{
    real h2n(real h, vector a) {
        return a[1] + a[2] * h + a[3] * h ^ 2 + a[4] * h ^ 3 + a[5] * h ^ 4;
    }
    real hq2n(real h, real q, vector a) {
        return a[1] + a[2] * h + a[3] * q + a[4] * h ^ 2 + 
               a[5] * q ^ 2 + a[6] * h * q + a[7] * h ^ 3 + 
               a[8] * q ^ 3 + a[9] * q * h ^ 2 + a[10] * h * q ^ 2 + 
               a[11] * h ^ 4 + a[12] * q ^ 4 + a[13] * q * h ^ 3 + 
               a[14] * h * q ^ 3 + a[15] * (h * q) ^ 2;
    }
    //straight line (triangular) distribution function
    real new_straight_line_lpdf(real y, real a, real b) {
	return log(a*y + b);
    }
    //function following the histogram shape
    real new_hist_lpdf(real y, int n_bins, vector z_bins, vector histogram) {
	int index=1;	
	for(i in 1:n_bins-1) {
	    if (y >= z_bins[i] && y < z_bins[i+1]) {
		index = i;
	    }
	}  	
	return log(histogram[index]);  
    }
}
data {
    int<lower=0, upper=1> ntlo;        // use next-to-leading-order expansion
    int<lower=0, upper=1> vary_m_c;    // variable chirp mass
    int<lower=0, upper=1> z_dep_rate;  // redshift merger rate
    int<lower=0, upper=1> fixed_n_bns; // assume sample size known
    int<lower=0> n_bns;                // total number of mergers
    vector[n_bns] obs_amp_plus;        // measured plus amplitude
    vector[n_bns] obs_amp_cross;       // measured cross amplitude
    
    real a_coeff;                      // a coefficient of straight line dist   
    real b_coeff;                      // b coefficient of straight line dist
    real z_min;                        // minimum prior redshift
    real z_max;                        // maximum prior redshift

    int n_bins;                        //number of bins for histogram
    vector[n_bins-1] histogram;        //values of histogram
    vector[n_bins] z_bins;             //redshift bins of histogram

    
    vector[n_bns] obs_m_c_z;           // measured redshifted chirp masses
    real amp_s;                        // intrinsic GW amplitude
    real amp_n;                        // GW noise level
    real sig_v_pec;                    // std of true peculiar velocities
    
    real sig_z;                        // noise on observed redshifts
    real mu_m_c;                       // chirp mass prior mean
    real sig_m_c;                      // chirp mass prior std
    real sig_obs_m_c_z;                // noise on observed redshifted chirp masses
    int<lower=0> n_coeffs;             // number of coefficients of polynomial fit to \bar{N}(H_0,q_0)
    vector[n_coeffs] n_bar_det_coeffs; // coefficients of polynomial fit to \bar{N}(H_0,q_0)
}
transformed data {
    real c;                            // c in km/s
    real g;                            // g in Mpc / M_sol / (km/s)^2
    int n_m_c;                         // number of chirp masses to sample
    c = 2.99792458e5;
    g = 4.301e-9;
    n_m_c = vary_m_c * n_bns;
}
parameters {
    real<lower=20.0, upper=140.0> h_0;
    vector<lower=0.07296, upper=0.3>[n_bns] z_cos;  //limits on z_cos come from the initial z_cos distribution used for GW simulation
    vector<lower=-1.0, upper=1.0>[n_bns] cos_i;
    vector[n_bns] v_pec;    
}

transformed parameters {

    vector<lower=z_min, upper=z_max>[n_bns] true_z;
    vector<lower=0.0>[n_bns] true_d;
    vector[n_bns] true_amp_plus;
    vector[n_bns] true_amp_cross;
    real<lower=0.0> n_bar_det;

    for(i in 1:n_bns) {
        
        // pick order-appropriate distance-redshift relation
        if (ntlo) {
            true_z[i] = z_cos[i] + (1.0 + z_cos[i]) * 
                        v_pec[i] / c;
            
        } else {
            true_z[i] = z_cos[i] + v_pec[i] / c;
            true_d[i] = c * z_cos[i] / h_0;
        }
        if (true_d[i] < 0) {
            if (ntlo) {
                print("D BAD! ", true_d[i], " ", z_cos[i], " ", h_0, " ");
            } else {
                print("D BAD! ", n_bar_det, " ", z_cos[i], " ", h_0);
            }
        }

        // different amplitudes if sampling chirp masses
        if (n_m_c > 0) {
            true_amp_plus[i] = amp_s; 
        } else {
            true_amp_plus[i] = amp_s * (1.0 + cos_i[i] ^ 2) / 
                               2.0 / true_d[i];
            true_amp_cross[i] = -amp_s * cos_i[i] / true_d[i];
        }

    }
    if (ntlo) {
        n_bar_det = hq2n(h_0, h_0, n_bar_det_coeffs);
    } else {
        n_bar_det = h2n(h_0, n_bar_det_coeffs);
    }
    if (n_bar_det < 0) {
        if (ntlo) {
            print("N BAD! ", n_bar_det, " ", h_0, " ");
        } else {
            print("N BAD! ", n_bar_det, " ", h_0);
        }
    }
    
    
}

model {

    // priors on true parameters. priors on cos(i) are uniform
    
    v_pec ~ normal(0.0, sig_v_pec);
    //v_pec ~ normal(0.0, sig_v_pec + c*sig_z); // to take into account redshift measurement uncertainty
    for(i in 1:n_bns) {    
	//target += new_straight_line_lpdf(true_z[i] | a_coeff, b_coeff); //for straight line distribution
	target += new_hist_lpdf(true_z[i] | n_bins, z_bins, histogram);   //for histogram 
    }

    // pick order-appropriate volume element. NB: as addition acts 
    // per vector element, the statement below (correctly) applies a 
    // factor of 1/H_0^3 per object
    target += -3.0 * log(h_0) + 2.0 * log(z_cos);
    

    // Poisson exponent
    //if (fixed_n_bns) {
    //    target += -n_bns * log(n_bar_det);
    //} else {
    //    target += -n_bar_det;
    //}

    // GW likelihoods
    obs_amp_plus ~ normal(true_amp_plus, amp_n);
    obs_amp_cross ~ normal(true_amp_cross, amp_n);
   

}


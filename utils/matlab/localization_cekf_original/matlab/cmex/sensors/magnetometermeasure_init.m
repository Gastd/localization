function magnetometermeasure = magnetometermeasure_init(flagnoise)

magnetometermeasure.mx = 0;
magnetometermeasure.my = 0;
magnetometermeasure.mz = 0;

variance_mx = 1*(0.0624)^2;
variance_my = 1*(0.0274)^2;
variance_mz = 1*(0.0546)^2;

magnetometermeasure.mxvariance = flagnoise*variance_mx;
magnetometermeasure.myvariance = flagnoise*variance_my;
magnetometermeasure.mzvariance = flagnoise*variance_mz;

magnetometermeasure.flagvalidmeasure = 1;

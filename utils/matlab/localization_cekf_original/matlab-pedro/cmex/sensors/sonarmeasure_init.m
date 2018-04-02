% Função que simula o sonar
function sonarmeasure = sonarmeasure_init(flagnoise)

variance_range = (0.01)^2;
sonarmeasure.rangevariance = flagnoise*variance_range;
sonarmeasure.flagvalidmeasure = 1;
sonarmeasure.range = 0;
sonarmeasure.R_s2b = eye(3);
sonarmeasure.t_s2b = zeros(3,1);
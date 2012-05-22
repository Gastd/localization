% Fun��o que simula o GPS
function gpsmeasure = gpsmeasure_init(flagnoise)

P_p = diag([0.3 0.3 0.3].^2); 
P_v = diag([0.05 0.05 05].^2); 

gpsmeasure.P_p = flagnoise*P_p;
gpsmeasure.P_v = flagnoise*P_v;
gpsmeasure.flagvalidpmeasure = 1;
gpsmeasure.flagvalidvmeasure = 1;
gpsmeasure.p = zeros(3,1);
gpsmeasure.v = zeros(3,1);

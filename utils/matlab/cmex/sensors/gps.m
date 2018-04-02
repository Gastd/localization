% Função que simula o GPS
function gpsmeasure = gps(vehiclestate, flagnoise)

gpsmeasure = gpsmeasure_init(flagnoise);
if isstruct(vehiclestate)==0
    return;
end

%calcula os supostas valores de medida do GPS
gpsmeasure.p = [vehiclestate.x; vehiclestate.y; vehiclestate.z];
gpsmeasure.v = [vehiclestate.dx_dt; vehiclestate.dy_dt; vehiclestate.dz_dt];

%adiciona o ruido caso habilitado
if flagnoise
    A = chol(gpsmeasure.P_p)';
    gpsmeasure.p = gpsmeasure.p + flagnoise*A*randn(3,1);
    A = chol(gpsmeasure.P_v)';
    gpsmeasure.v = gpsmeasure.v + flagnoise*A*randn(3,1);
end
gpsmeasure.flagvalidpmeasure = 1;
gpsmeasure.flagvalidvmeasure = 1;

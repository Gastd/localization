%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_out = localization_filter_model_f(procedure_name,x,u,G,T,flagestimateaccelerometerbias, flagestimategyrometerbias)
% x : [q; p; v]
% u : [imu.ax; imu.ay; imu.az; imu.wx; imu.wy; imu.wz]
% w : process noise

if nargin == 6
    flagestimategyrometerbias = 0;
end

% vari�veis
imumeasure.ax = u(1);
imumeasure.ay = u(2);
imumeasure.az = u(3);
imumeasure.wx = u(4);
imumeasure.wy = u(5);
imumeasure.wz = u(6);

if flagestimategyrometerbias
    sx = (imumeasure.wx-x(14))*T;
    sy = (imumeasure.wy-x(15))*T;
    sz = (imumeasure.wz-x(16))*T;
else
    sx = imumeasure.wx*T;
    sy = imumeasure.wy*T;
    sz = imumeasure.wz*T;
end
v = sqrt(sx^2 + sy^2 + sz^2);

switch procedure_name
    case 'evaluate'
        imumeasure.ax = u(1);
        imumeasure.ay = u(2);
        imumeasure.az = u(3);
        imumeasure.wx = u(4);
        imumeasure.wy = u(5);
        imumeasure.wz = u(6);

        %Atualiza a atitude
        if(v > 0)
            q(1) = cos(v/2)*x(1) - sin(v/2)/v*( sx*x(2) + sy*x(3) + sz*x(4));
            q(2) = cos(v/2)*x(2) - sin(v/2)/v*(-sx*x(1) - sz*x(3) + sy*x(4));
            q(3) = cos(v/2)*x(3) - sin(v/2)/v*(-sy*x(1) + sz*x(2) - sx*x(4));
            q(4) = cos(v/2)*x(4) - sin(v/2)/v*(-sz*x(1) - sy*x(2) + sx*x(3));
        else
            q(1) = x(1);
            q(2) = x(2);
            q(3) = x(3);
            q(4) = x(4);
        end
        %preserva a atitude, posi��o e velocidade pr�via
        if flagestimateaccelerometerbias
            a = quaternions2dcm(x(1:4))*[imumeasure.ax-x(11); imumeasure.ay-x(12); imumeasure.az-x(13)] + G; % gravity compensation
        else
            a = quaternions2dcm(x(1:4))*[imumeasure.ax; imumeasure.ay; imumeasure.az] + G; % gravity compensation
        end
        p = x(5:7) + x(8:10)*T + a*(T^2)/2;
        v = x(8:10) + T*a;
        % v = (rotation_from_quaternions(x(1:4),[imumeasure.ax; imumeasure.ay; imumeasure.az]) + G)*T + x(8:10);
        % p = x(8:10)*T + x(5:7);

        %Armazena a posi��o
        var_out(1:4,1)  = q;
        var_out(5:7,1)  = p;
        var_out(8:10,1) = v;
        
        % par�metros:
        if flagestimateaccelerometerbias
            var_out(11:13,1) = x(11:13,1);
        end
        
        if flagestimategyrometerbias
            var_out(14:16,1) = x(14:16,1);
        end
        
    case 'df_dx'
        df_dx = zeros(length(x),length(x));

        if(v > 0)
            df_dx(1,1) = cos(v/2);
            df_dx(1,2) = -(sin(v/2)/v)*sx;
            df_dx(1,3) = -(sin(v/2)/v)*sy;
            df_dx(1,4) = -(sin(v/2)/v)*sz;

            df_dx(2,1) = (sin(v/2)/v)*sx;
            df_dx(2,2) = cos(v/2);
            df_dx(2,3) = (sin(v/2)/v)*sz;
            df_dx(2,4) = -(sin(v/2)/v)*sy;

            df_dx(3,1) = -(sin(v/2)/v)*sy;
            df_dx(3,2) = (sin(v/2)/v)*sz;
            df_dx(3,3) = cos(v/2);
            df_dx(3,4) = -(sin(v/2)/v)*sx;

            df_dx(4,1) = -(sin(v/2)/v)*sz;
            df_dx(4,2) = -(sin(v/2)/v)*sy;
            df_dx(4,3) = (sin(v/2)/v)*sx;
            df_dx(4,4) = cos(v/2);

            R = quaternions2dcm(x(1:4));
            dR_dx1 = 2*[(2*x(1))/2    (-x(4))       (x(3));
                        ( x(4))       (2*x(1))/2    (-x(2));
                        (-x(3))       (x(2))        (2*x(1))/2];
            dR_dx2 = 2*[(2*x(2))/2    (x(3))        (x(4));
                        (x(3))        (-2*x(2))/2   (-x(1));
                        (x(4))        (x(1))        (-2*x(2))/2];
            dR_dx3 = 2*[(-2*x(3))/2   (x(2))        (x(1));
                        (x(2))        (2*x(3))/2    (x(4));
                        (-x(1))       (x(4))        (-2*x(3))/2];
            dR_dx4 = 2*[(-2*x(4))/2   -x(1)         (x(2));
                        (x(1))        (-2*x(4))/2    (x(3));
                        (x(2))        (x(3))         (2*x(4))/2];

            da_dx = zeros(3,4); % somente quaternions
            da_dx(:,1) = dR_dx1*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
            da_dx(:,2) = dR_dx2*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
            da_dx(:,3) = dR_dx3*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
            da_dx(:,4) = dR_dx4*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions

            df_dx(5:7,1:4)  = da_dx*(T^2)/2;
            df_dx(5:7,5:7)  = eye(3);
            df_dx(5:7,8:10) = eye(3)*T;
            if flagestimateaccelerometerbias
                df_dx(5:7,11:13) = -R*(T^2)/2;
            end

            df_dx(8:10,1:4)  = da_dx*(T);
            df_dx(8:10,8:10) = eye(3);
            if flagestimateaccelerometerbias
                df_dx(8:10,11:13) = -R*(T);
            end

            if flagestimateaccelerometerbias
                df_dx(11:13,11:13) = eye(3);
            end

            if flagestimategyrometerbias
                df_dx(14:16,14:16) = eye(3);
            end
        end
        
        var_out = df_dx;
        
    case 'df_du'
        df_du = zeros(length(x),length(u));
    
        if(v > 0)
            dv_dwx = (1/2)*(1/v)*(2*sx*T);
            dv_dwy = (1/2)*(1/v)*(2*sy*T);
            dv_dwz = (1/2)*(1/v)*(2*sz*T);

            df_du(1,3) = -sin(v/2)*(1/2)*dv_dwx*x(1) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwx * ( sx*x(2) + sy*x(3) + sz*x(4)) - sin(v/2)/v*( T*x(2) ); % df_dwx
            df_du(1,4) = -sin(v/2)*(1/2)*dv_dwy*x(1) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwy * ( sx*x(2) + sy*x(3) + sz*x(4)) - sin(v/2)/v*( T*x(3) ); % df_dwy
            df_du(1,5) = -sin(v/2)*(1/2)*dv_dwz*x(1) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwz * ( sx*x(2) + sy*x(3) + sz*x(4)) - sin(v/2)/v*( T*x(4) ); % df_dwz

            df_du(2,3) = -sin(v/2)*(1/2)*dv_dwx*x(2) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwx * (-sx*x(1) - sz*x(3) + sy*x(4)) - sin(v/2)/v*(-T*x(1) ); % df_dwx
            df_du(2,4) = -sin(v/2)*(1/2)*dv_dwy*x(2) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwy * (-sx*x(1) - sz*x(3) + sy*x(4)) - sin(v/2)/v*( T*x(4) ); % df_dwy
            df_du(2,5) = -sin(v/2)*(1/2)*dv_dwz*x(2) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwz * (-sx*x(1) - sz*x(3) + sy*x(4)) - sin(v/2)/v*(-T*x(3) ); % df_dwz

            df_du(3,3) = -sin(v/2)*(1/2)*dv_dwx*x(3) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwx * (-sy*x(1) + sz*x(2) - sx*x(4)) - sin(v/2)/v*(-T*x(4) ); % df_dwx
            df_du(3,4) = -sin(v/2)*(1/2)*dv_dwy*x(3) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwy * (-sy*x(1) + sz*x(2) - sx*x(4)) - sin(v/2)/v*(-T*x(1) ); % df_dwy
            df_du(3,5) = -sin(v/2)*(1/2)*dv_dwz*x(3) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwz * (-sy*x(1) + sz*x(2) - sx*x(4)) - sin(v/2)/v*( T*x(2) ); % df_dwz

            df_du(4,3) = -sin(v/2)*(1/2)*dv_dwx*x(4) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwx * (-sz*x(1) - sy*x(2) + sx*x(3)) - sin(v/2)/v*( T*x(3) ); % df_dwx
            df_du(4,4) = -sin(v/2)*(1/2)*dv_dwy*x(4) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwy * (-sz*x(1) - sy*x(2) + sx*x(3)) - sin(v/2)/v*(-T*x(2) ); % df_dwy
            df_du(4,5) = -sin(v/2)*(1/2)*dv_dwz*x(4) + (-cos(v/2)*(1/2)*v + sin(v/2))/(v*v)* dv_dwz * (-sz*x(1) - sy*x(2) + sx*x(3)) - sin(v/2)/v*(-T*x(1) ); % df_dwz

            R = quaternions2dcm(x(1:4));
            df_du(5:7 ,1:3)  = R*(T^2)/2;
            df_du(8:10,1:3)  = R*(T);
        end
        
        var_out = df_du;
        
end
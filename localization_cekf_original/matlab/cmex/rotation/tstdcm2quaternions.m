clear all;
close all;

N = 10000;
for n=1:N
    roll = pi*(2*rand(1)-1);
    pitch = pi/2*(2*rand(1)-1);
    yaw = pi*(2*rand(1)-1);
    
    a = roll;
    Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
    a = pitch;
    Ry = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)];
    a = yaw;
    Rz = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
    
    R = Rz*Ry*Rx;
    [q1] = matdcm2quaternions(R);
    [q2] = mexdcm2quaternions(R);
    
    d = q1-q2;
    if sqrt(d*d') > 1e-7
        disp('Erro para o caso seguinte:');
        roll
        pitch
        yaw
        disp('Com a seguinte diferença:');
        d
        return;
    end
end
disp('Teste concluido com sucesso');

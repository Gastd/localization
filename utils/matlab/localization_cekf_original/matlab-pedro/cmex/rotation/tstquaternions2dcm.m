clear all;
close all;

N = 10000;
for n=1:N
    roll = pi*(2*rand(1)-1);
    pitch = pi/2*(2*rand(1)-1);
    yaw = pi*(2*rand(1)-1);

    [q] = mexeuler2quaternions(roll,pitch,yaw);
    
    [R1] = mexquaternions2dcm(q);
    [R2] = matquaternions2dcm(q);
    
    d = R1-R2;
    if sqrt(d*d') > 1e-7
        disp('Erro para o caso seguinte:');
        roll
        pitch
        yaw
        q1
        q2
        disp('Com a seguinte diferenša:');
        d
        return;
    end
end
disp('Teste concluido com sucesso');

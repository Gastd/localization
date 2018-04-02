clear all;
close all;

N = 10000;
for n=1:N
    roll = pi*(2*rand(1)-1);
    pitch = pi/2*(2*rand(1)-1);
    yaw = pi*(2*rand(1)-1);
    
    [R1] = mexeuler2dcm([roll,pitch,yaw]');
    [R2] = mateuler2dcm([roll,pitch,yaw]');
    
    d = R1-R2;
    if sqrt(d*d') > 1e-7
        disp('Erro para o caso seguinte:');
        roll
        pitch
        yaw
        R1
        R2
        disp('Com a seguinte diferença:');
        d
        return;
    end
end
disp('Teste concluido com sucesso');

%funcao que organiza os dados para que sejam desenhado pela funcao chamada

function [hvehicle] = vehicle_draw(vehiclestate,hvehicle)

x = vehiclestate.x;
y = vehiclestate.y;
z = vehiclestate.z;
pitch = vehiclestate.pitch;
yaw = vehiclestate.yaw;
roll = vehiclestate.roll;

%passa os angulos para o sistema de referencia da terra
R = euler2dcm([roll,pitch,yaw]');
P = [x;y;z];
Length = 3; % length of each arrow.

if nargin == 1 %verifica se o numero de argumento eh 1
   %%% first creation
   hvehicle = coordinate_system_draw(P,R,Length,'X_b','Y_b','Z_b','k');
else
   %%% update
   hvehicle = coordinate_system_draw(P,R,Length,'X_b','Y_b','Z_b','k',hvehicle);
end
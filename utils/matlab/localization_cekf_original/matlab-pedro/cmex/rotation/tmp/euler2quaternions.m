%Observar que na linha 1 coluna 2 há um erro na Matriz do Padilha
function q = euler2quaternions(varargin)

if nargin==1
    v = varargin{1};
    roll = v(1);
    pitch = v(2);
    yaw = v(3);
end
if nargin==3
    roll = varargin{1};
    pitch = varargin{2};
    yaw = varargin{3};
end

R = euler2dcm(roll,pitch,yaw);
q = dcm2quaternions(R);
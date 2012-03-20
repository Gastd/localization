%Observar que na linha 1 coluna 2 há um erro na Matriz do Padilha
function R = euler2dcm(rpy)

% if nargin==1
%     v = varargin{1};
%     roll = v(1);
%     pitch = v(2);
%     yaw = v(3);
% end
% if nargin==3
%     roll = varargin{1};
%     pitch = varargin{2};
%     yaw = varargin{3};
% end

[R] = mexeuler2dcm(rpy);
%[R] = mateuler2dcm(rpy);
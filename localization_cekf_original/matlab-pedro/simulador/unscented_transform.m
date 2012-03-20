%Transformada Unscented

%   Parametros de entrada
%       Xm -> vetor em torno do qual são calculados os pontos sigma
%       Px -> matriz de covariancia associada a Xm
%       functionname -> funcao nao-linear
%
%
%
%   Parametros de saida
%       Ym -> media transformada
%       Py -> matriz de covariancia transformada

function [Ym,Py] = unscented_transform(Xm,Px,functionname,Xanglemask,Yanglemask)

%%% Propagação de incertezas por Unscented Transform.
DimX = length(Xm);
DimY = length(eval(sprintf('%s(Xm);',functionname)));

if nargin < 4,
    Xanglemask = zeros(length(Xm),1);
end
if nargin < 5,
    Yanglemask = zeros(length(Ym),1);
end

kappa = 1;
UtNSamples = 2*DimX+1;
XSamples = zeros(DimX,UtNSamples);
YSamples = zeros(DimY,UtNSamples);
WSamples = zeros(1,UtNSamples);
Q = sqrtm((DimX+kappa)*Px)';
for i=1:DimX,
   WSamples(1,i)   = 1/(2*(DimX+kappa));
   WSamples(1,i+DimX) = 1/(2*(DimX+kappa));
   XSamples(:,i)   = Xm + Q(:,i);
   XSamples(:,i+DimX) = Xm - Q(:,i);
   for j=1:DimX
       if Xanglemask(j)
           XSamples(j,i) = atan2(sin(XSamples(j,i)),cos(XSamples(j,i)));
           XSamples(j,i+DimX) = atan2(sin(XSamples(j,i+DimX)),cos(XSamples(j,i+DimX)));
       end
   end
   YSamples(:,i)   = eval(sprintf('%s(XSamples(:,i));',functionname));
   YSamples(:,i+DimX) = eval(sprintf('%s(XSamples(:,i+DimX));',functionname));
end,
WSamples(1,UtNSamples) = kappa/(DimX+kappa);
XSamples(:,UtNSamples) = Xm;
YSamples(:,UtNSamples) = eval(sprintf('%s(XSamples(:,UtNSamples));',functionname)); 

Ym 	= zeros(DimY,1);
Py  = zeros(DimY,DimY);
for ns=1:UtNSamples
   Ym = Ym + WSamples(1,ns)*YSamples(:,ns);
end,
if(sum(isreal(Ym)==0)>0)
    eig(Px)
    disp('Ym tem componente complexo');
    pause
end
for j=1:DimY
    if Yanglemask(j)
        Ym(j) = atan2(sin(Ym(j)),cos(Ym(j)));
    end
end
for ns=1:UtNSamples
    dy = YSamples(:,ns)-Ym;
    for j=1:DimY
        if Yanglemask(j)
            dy(j) = atan2(sin(dy(j)),cos(dy(j)));
        end
    end
   Py = Py + WSamples(1,ns)*(dy)*(dy)';
end,

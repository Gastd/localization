%funcao que cria e atualiza as setas orientadas que caracterizam a posição
%do veiculo
function [hcoordinatesystem] = coordinate_system_draw(P,R,Length,Xstring,Ystring,Zstring,Color,hcoordinatesystem)

if nargin == 7
   %%% first creation
   hcoordinatesystem = zeros(6,1);
   Start = P; %posição do inicio das setas
   Stop = P+Length*R(:,1); %posição final da seta
   %desenha uma seta orientada
   hcoordinatesystem(1) = arrow('Start',Start,'Stop',Stop,'LineStyle','-','BaseAngle',45,'Length',Length*2); hold on;
   hcoordinatesystem(2) = text(Stop(1),Stop(2),Stop(3),Xstring);
   Stop = P+Length*R(:,2);
   hcoordinatesystem(3) = arrow('Start',Start,'Stop',Stop,'LineStyle','-','BaseAngle',45,'Length',Length*2);
   hcoordinatesystem(4) = text(Stop(1),Stop(2),Stop(3),Ystring);
   Stop = P+Length*R(:,3);
   hcoordinatesystem(5) = arrow('Start',Start,'Stop',Stop,'LineStyle','-','BaseAngle',45,'Length',Length*2);
   hcoordinatesystem(6) = text(Stop(1),Stop(2),Stop(3),Zstring);
else
   %%% update
   Start = P;
   Stop = P+Length*R(:,1);
   arrow(hcoordinatesystem(1),'Start',Start,'Stop',Stop,'LineStyle','-','BaseAngle',45);
   set(hcoordinatesystem(2),'Position',Stop);
   Stop = P+Length*R(:,2);
   arrow(hcoordinatesystem(3),'Start',Start,'Stop',Stop,'LineStyle','-','BaseAngle',45);
   set(hcoordinatesystem(4),'Position',Stop);
   Stop = P+Length*R(:,3);
   arrow(hcoordinatesystem(5),'Start',Start,'Stop',Stop,'LineStyle','-','BaseAngle',45);
   set(hcoordinatesystem(6),'Position',Stop);
end
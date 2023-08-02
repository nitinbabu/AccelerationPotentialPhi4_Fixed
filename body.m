
function [bodyX,bodyY,centx,centy] = body(side,R)
side = 0.5;
R = 0.064;
nop = 11;
nor = 10; 
deficit = side-(R);
xg = 0;
yg = 0;

topX =  linspace(R,deficit,nop) ; topY =  zeros(1,length(topX));
topX(1) = []; topX(end) = [];topY(1) = []; topY(end) = [];

rightY =  linspace(-R,-deficit,nop) ; rightX =  side+zeros(1,length(rightY));
rightY(1) = [];rightY(end) = [];    rightX(1) = [];rightX(end) = [];

bottomX =  linspace(deficit,R,nop) ;
bottomY =  -side+zeros(1,length(bottomX));

bottomY(1) = [];bottomY(end) = [];    bottomX(1) = [];bottomX(end) = [];



leftY =  linspace(-deficit,-R,nop) ; leftX =  zeros(1,length(leftY));
leftY(1) = [];  leftY(end) = [];    leftX(1) = []; leftX(end) = [];


bilgeTRx = (side-R)+R*cos(linspace(pi/2,0,nor));
bilgeTRy = -R+R*sin(linspace(pi/2,0,nor));


bilgeLRx = (side-R)+R*cos(linspace(0,-pi/2,nor));
bilgeLRy = (-side+R)+R*sin(linspace(0,-pi/2,nor));


bilgeLLx = (R)+R*cos(linspace(-pi/2,-pi,nor));
bilgeLLy = (-side+R)+R*sin(linspace(-pi/2,-pi,nor));


bilgeTLx = R+R*cos(linspace(pi,pi/2,nor));
bilgeTLy = -R+R*sin(linspace(pi,pi/2,nor));

% figure;
% plot(topX,topY,'.-')
% hold on
% plot(rightX,rightY,'.-')
% plot(bottomX,bottomY,'.-')
% plot(leftX,leftY,'.-')
% axis equal
% grid minor
% 
% plot(side-R,-R,'.r')
% plot(bilgeTRx,bilgeTRy,'o-b');
% hold on
% plot(R,-side+R,'.r')
% plot(bilgeLRx,bilgeLRy,'o-b');
% 
% 
% plot(side-R,-side+R,'.r')
% plot(bilgeLLx,bilgeLLy,'o-b');
% 
% 
% plot(R,-R,'.r')
% plot(bilgeTLx,bilgeTLy,'o-b');


bodyX = [topX bilgeTRx rightX bilgeLRx bottomX bilgeLLx leftX bilgeTLx];
bodyY = [topY bilgeTRy rightY bilgeLRy bottomY bilgeLLy leftY bilgeTLy];

centx = side/2; centy = -side/2;
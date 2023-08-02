% manbody
%% Manipulates body poition
clc
clear
close all
[bodyX,bodyY,centx,centy] = body(0.5,0.064); % Body at false origin

xg = -centx ; yg = -centy; 
xbody = bodyX+xg; % Shift to designated origin at xoz
ybody = bodyY+yg;
figure;
hold on
grid minor
axis equal
plot(xbody,ybody);
theta = 05;
rxbody = xbody.*cosd(theta)-ybody.*sind(theta);
rybody = xbody.*sind(theta)+ybody.*cosd(theta);

hold on
plot(rxbody,rybody,'.-k');

testx = 1*cos(0:pi/120:2*pi)
curver = curvspace([rxbody; rybody ]',1000);
curx = curver(:,1); cury = curver(:,2);
% cury = curvspace(rybody,100);

plot(curx,cury,'.-m')

%%
[newx,newy,~] = CSregrid(curx,cury,0.*curx)


plot(newx,newy,'.-.');
legend('old','shift','rotate','spline')



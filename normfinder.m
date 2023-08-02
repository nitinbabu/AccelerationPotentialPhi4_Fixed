

clc
clear
close all

%% This bit of code keeps the node on the wetted parts. Answers how the water knows the body is present there
[bodyX,bodyY,centx,centy] = body(0.5,0.064); % Body at false origin

xg = -centx ; yg = -centy; 
xbody = bodyX+xg; % Shift to designated origin at xoz
ybody = bodyY+yg;
figure;
hold on
grid minor
axis equal
plot(xbody,ybody);
[cx,cy]=CenterMass(xbody,ybody)

% This was in body fixed coordinates
%%

cx = 8;
theta = 0;
% rxbody = xbody.*cosd(theta)-ybody.*sind(theta);
% rybody = xbody.*sind(theta)+ybody.*cosd(theta);
sway = 0; heave = 0;
rbody= [cx;cy]+[sway;heave]+[cosd(theta)  sind(theta); -sind(theta) cosd(theta)]*[xbody; ybody];

rxbody = rbody(1,:);
rybody = rbody(2,:);


hold on
plot(rxbody,rybody,'.-k');
%%
% testx = 1*cos(0:pi/120:2*pi)
curver = curvspace([rxbody; rybody ]',1000);
curx = curver(:,1); cury = curver(:,2);
% cury = curvspace(rybody,100);

plot(curx,cury,'.-m')



% x = curx;
% y = cury;
% 
% 
% 
% xbody = x;
% ybody = y;
% plot(x,y,'-x')

%% Left

Px = -.25+0.1; Py = 0.05;

hold on
plot(Px,Py,'bo')

%% Right
Kx = .25; Ky = -0.05;
hold on
plot(Kx,Ky,'bo')
%%
newpoints = curvspace([curx,cury],1000)
newrx =newpoints(:,1)';
newry =newpoints(:,2)';
plot(newrx,newry,'.-r')
%% Left Wetted Redist
rdist = sqrt((curx-Px).^2+(cury-Py).^2);
[val,locl]=min(rdist);
Qx = newrx(locl);
Qy = newry(locl);
plot(Qx,Qy,'sr')

%% Right Wetted Redist 



rdist = sqrt((curx-Kx).^2+(cury-Ky).^2);
[val,locr]=min(rdist);
Qx = newrx(locr);
Qy = newry(locr);
plot(Qx,Qy,'sr')

%%

wetx = newrx(locr:locl);
wety = newry(locr:locl);

% plot(wetx,wety,'.k')


wetx= flip(wetx);
wety= flip(wety);


figure;
plot(wetx,wety,'.-k')


hold on

newnodes = curvspace([wetx',wety'],30); %% This must be output
x = newnodes(:,1);
y = newnodes(:,2);

plot(x,y,'o-')
axis equal

[cx,cy]=CenterMass(xbody,ybody)
%% Manipulation


% rxbody = xbody.*cosd(theta)-ybody.*sind(theta);
% rybody = xbody.*sind(theta)+ybody.*cosd(theta);


% Calculate the area
theta = 5;
sway = -0.05; heave = -0.05;



rbody= 1.*[cx;cy]+[sway;heave]+[cosd(theta)  sind(theta); -sind(theta) cosd(theta)]*[xbody; ybody];

rxbody = rbody(1,:);
rybody = rbody(2,:);
plot(rxbody,rybody,'d-r')
hold on
axis([-2 2 -1 1]);
axis equal
% grid minor

[cx,cy]=CenterMass(rxbody,rybody);
%% Redist Projection
newpoints = curvspace([rxbody;rybody]',1000);
newrx =newpoints(:,1)';
newry =newpoints(:,2)';
plot(newrx,newry,'.-r')

%% Left Redist
rdist = sqrt((newrx-Px).^2+(newry-Py).^2);
[val,locl]=min(rdist);
Qx = newrx(locl);
Qy = newry(locl);
plot(Qx,Qy,'sr')

%% Right Redist



rdist = sqrt((newrx-Kx).^2+(newry-Ky).^2);
[val,locr]=min(rdist);
Qx = newrx(locr);
Qy = newry(locr);
plot(Qx,Qy,'sr')

%%

wetx = newrx(locr:locl);
wety = newry(locr:locl);

plot(wetx,wety,'.k')


wetx= flip(wetx);0
wety= flip(wety);


figure;
plot(wetx,wety,'.-k')

hold on

newnodes = curvspace([wetx',wety'],30);
x = newnodes(:,1);
y = newnodes(:,2);

plot(x,y,'o')
axis equal
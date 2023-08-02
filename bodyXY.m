function [xbody,ybody,wetx,wety,cx,cy] = bodyXY(L,wpoints,theta,sway,heave)


    side = 0.5;
    [bodyX,bodyY,centx,centy] = body(side,0.064); % Body at false origin
    xg = -centx ; yg = -centy; 
    xbody = bodyX+xg; % Shift to designated origin at xoz
    ybody = bodyY+yg;



% figure;
% hold on
% grid minor
% axis equal
% plot(xbody,ybody);
[cx,cy]=CenterMass(xbody,ybody);

% This was in body fixed coordinates
%% This Body is placed....relative to the origin.

cx = L/2;  % Location relative to the Wave tank origin
cy; 
theta = 0;  sway = 0; heave = 0;

rbody= [cx;cy]+[sway;heave]+[cosd(theta)  sind(theta); -sind(theta) cosd(theta)]*[xbody; ybody];
rxbody = rbody(1,:);    rybody = rbody(2,:);

% hold on
% plot(rxbody,rybody,'.-k');
%% 
curver = curvspace([rxbody; rybody ]',1000);
curx = curver(:,1); cury = curver(:,2);
% plot(curx,cury,'.-m')


%% Left
Px = cx-side/2; Py = 0;
% hold on
% plot(Px,Py,'bo')

%% Right
Kx = cx+side/2; Ky = 0;
hold on
% plot(Kx,Ky,'bo')
%%
newpoints = curvspace([curx,cury],1000)
newrx =newpoints(:,1)';
newry =newpoints(:,2)';
% plot(newrx,newry,'.-r')
%% Left Wetted Redist
rdist = sqrt((curx-Px).^2+(cury-Py).^2);
[val,locl]=min(rdist);
Qx = newrx(locl);
Qy = newry(locl);
% plot(Qx,Qy,'sr')

%% Right Wetted Redist 



rdist = sqrt((curx-Kx).^2+(cury-Ky).^2);
[val,locr]=min(rdist);
Qx = newrx(locr);
Qy = newry(locr);
% plot(Qx,Qy,'sr')

%% Plots wetted body surface.. 

wetx = newrx(locr:locl);
wety = newry(locr:locl);


wetx= flip(wetx);
wety= flip(wety);

% figure;
% plot(wetx,wety,'.-k')

% hold on
newnodes = curvspace([wetx',wety'],wpoints); %% This must be output
wetx = newnodes(:,1);
wety = newnodes(:,2);

% plot(x,y,'o-')
% axis equal
[cx,cy]=CenterMass(xbody,ybody);
xbody = rxbody;
ybody = rybody;
nodePortX = wetx(1);    nodePortY = wety(1);
nodeStarX = wetx(end) ; nodeStarY = wety(end);








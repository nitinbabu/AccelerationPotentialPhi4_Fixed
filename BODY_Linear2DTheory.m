
clc
close all
clear
% tic

% Works for these combinations|| not sure why                                      
% especially when fx is 20/40/60/80/100
% fx = 40; L = 100; h = 20;
% fx = 60; L = 300; h = 30;
% fx = 40; L = 300; h = 60;


% fx =130;  h = 1; Amp = 0.01; lambda = 12.57; L = 6*lambda; T = 4.17; w = 2*pi/T; k = 2*pi/lambda; g = 9.81;

fx =51;  h = 1; Amp = 0.01; lambda = 1; L = 6*lambda; T = 0.8005; w = 2*pi/T; k = 2*pi/lambda; g = 9.81;
Vn = 0;
t0=0;
time = t0:0.1:10;
%%  eta = 0.5.*sin(0.09.*xf-1*t); Lambda = 70; T = 6;

%% Code breaks if there is some mismatch between designated elements and
%% actual elements due to positions
time = t0; % for any instant
t = t0
for m = 1
    t = time(m);
    
% [x,y,xf,xd,xb,xp,yf,yd,yb,yp,uin,uout,pot] = boundary(fx,h,L,t,w,k,Amp); %dx changed
            BL = L/2-0.2*L/2; BR = L/2+0.2*L/2;
            [x,y,xfl, xbody, xfr,xd,xb,xp, yfl,ybody,yfr,yd,yb,yp] = nwtbody(fx,h,L,w,k,Amp); %dx changed; %dx changed
%             [x,y,xf,xd,xb,xp,yf,yd,yb,yp,uin,uout,pot] = boundary(fx,h,L,t,w,k,Amp); %dx changed

            potL = 1.*(g*Amp/w)*(cosh(k*(1.*yfl+h))/cosh(k*h)).*sin(k.*xfl-w*t);
            potR = 1.*(g*Amp/w)*(cosh(k*(1.*yfr+h))/cosh(k*h)).*sin(k.*xfr-w*t);
            
            uin = -(g*k*Amp/w) * (cosh(k*(yp+h))/cosh(k*h)).*cos(k.*xp-w*t) ; % Need to check the multiple fluxes value at the very last element 0;
            uL = -(g*k*Amp/w) * (cosh(k*(yfl(1)+h))/cosh(k*h)).*cos(k.*xfl(1)-w*t);
            
            uout = 1.*(g*k*Amp/w) * (cosh(k*(yd+h))/cosh(k*h)).*cos(k.*xd-w*t) ;

[C1, C2, C3, C4, C5, C6] = bodycorners(xfl,xbody,xfr,xd,xb,xp);
nnode = length(x);  nelem  = ((length(x))/2);
% corner4 = [C1 C2 C3 C4];
corner6 = [C1 C2 C3 C4 C5 C6];


% Meshing
qnode =zeros(nelem,3);
for n = 1:nelem
     if n ~= nelem
            for i = 1:3
                qnode(n,1) = 2*(n-1)+1;
                qnode(n,2) = 2*(n-1)+2;
                qnode(n,3) = 2*(n-1)+3;
            end
     else
            for i = 1:3
                qnode(n,1) = 2*(n-1)+1;
                qnode(n,2) = 2*(n-1)+2;
                qnode(n,3) = 1;%2*(n-1)+3; % Closing final node for geometry
            end
     end
end

xint = x; yint = y;


 yint(C1:C2) = 0;  xint(C1:C2) = xfl; %
 yint(C3:C4) = 0;  xint(C3:C4) = xfr; %
figure(1)
plot(x,y,'.-'); grid on; title('2D Domain'); xlabel('x'); ylabel('y');
hold on; grid on; title('2D NWT');  axis([0 L -h 2*h])
[ktvec,knvec]=bcsB( xfl,xbody,xfr, xd, xb, xp,potL,Vn,potR,uout,uin)

[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEMBODY(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner6,uL);

toc


if m <2 
figure(2);
sgtitle('HOBEM Results')
subplot(1,2,1)
Handle1 = plot(1:length(x),p,'.-'); hold on

title('Potential ')
xlabel('nodes')
% pa = .2.*y+5;
% pa = [5 5 5 4 3 3 3 4];
hold on
% plot(pa,'o r')

% legend('mf-HOBEM','Analytical')
grid on
subplot(1,2,2)
Handle2 = plot(1:length(x),dpdn,'.-'); hold on
grid on
title('Gradient ')
hold on

% plot(qa,'o r')
xlabel('nodes')
% legend('mf-HOBEM','Analytical')

else
    set(Handle1,'XData',1:length(x),'YData', p);
                set(Handle2,'XData',1:length(x),'YData', dpdn);
                drawnow;
end
end
%%
% [U,V,NX,NY] = sderivatives(xf,yf,x,y,p,dpdn,qnode);
[UL,VL,NXL,NYL,detadxL] = sderivatives2LEFT(xfl,yfl,xint,yint,p,dpdn,qnode,corner6);
[UR,VR,NXR,NYR,detadxR] = sderivatives2RIGHT(xfr,yfr,xint,yint,p,dpdn,qnode,corner6);
[UB,VB,NXB,NYB,detadxB] = sderivatives2BODY(xbody,ybody,xint,yint,p,dpdn,qnode,corner6);
figure(11)
plot(C2+1:C3-1,NXB,'^-r')
hold on
figure(11)
plot(C2+1:C3-1,NYB,'v-k')

t = t0;
% clc
figure;
sgtitle('Seabed potential')

% subplot(1,2,1)
plot(xb,p(C5:C6),'-x')
grid on; hold on
xlabel('bottom nodes')
% subplot(1,2,2)
grid on
potb = (g*Amp/w)*(cosh(k*(yb+h))/cosh(k*h)).*sin(k.*xb-w*t);
plot(xb,potb,'o-k')
xlabel('bottom nodes')
grid on
legend('BEM','Airy')


uaL = (g*k*Amp/w) * (cosh(k.*(yfl+h))/cosh(k*h)).*cos(k.*xfl-w*t);
vaL = (g*k*Amp/w) * (sinh(k.*(yfl+h))/cosh(k*h)).*sin(k.*xfl-w*t);


uaR = (g*k*Amp/w) * (cosh(k.*(yfr+h))/cosh(k*h)).*cos(k.*xfr-w*t);
vaR = (g*k*Amp/w) * (sinh(k.*(yfr+h))/cosh(k*h)).*sin(k.*xfr-w*t);



figure;
sgtitle('Free Surface Velocities')
% subplot(1,2,1)
title('BEM')
plot(xfl,UL,'o-'); hold on; grid 
plot(xfl,uaL,'-x');
legend('U','ua')
xlabel('fs nodes')



figure;
sgtitle('Free Surface Velocities')
% subplot(1,2,1)
title('BEM')
plot(xfr,UR,'o-'); hold on; grid 
plot(xfr,uaR,'-x');
legend('U','ua')
xlabel('fs nodes')



% subplot(1,2,2)
 figure;
title('Airy')
plot(xfl,VL,'o-b');hold on; grid on
plot(xfl,vaL,'x-r')

legend('V','va')
xlabel('fs nodes')
% legend('U','ua','V','va')

 figure;
title('Airy')
plot(xfr,VR,'o-b');hold on; grid on
plot(xfr,vaR,'x-r')

legend('V','va')
xlabel('fs nodes')

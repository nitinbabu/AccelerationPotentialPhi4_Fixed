
clc
close all
clear
% tic

% Works for these combinations|| not sure why                                      
% especially when fx is 20/40/60/80/100
% fx = 40; L = 100; h = 20;
% fx = 60; L = 300; h = 30;
% fx = 40; L = 300; h = 60;


fx =60*1; L = 6; h = 1; Amp = 0.03; lambda = 1; T = 0.8005; w = 2*pi/T; k = 2*pi/lambda; g = 9.81;
t0=0;
time = t0:0.1:10;
%%  eta = 0.5.*sin(0.09.*xf-1*t); Lambda = 70; T = 6;

%% Code breaks if there is some mismatch between designated elements and
%% actual elements due to positions
time = t0; % for any instant
t = t0
for m = 1
    t = time(m);
    
[x,y,xf,xd,xb,xp,yf,yd,yb,yp,uin,uout,pot] = boundary(fx,h,L,t,w,k,Amp); %dx changed
            pot = 1.*(g*Amp/w)*(cosh(k*(1.*yf+h))/cosh(k*h)).*sin(k.*xf-w*t);
            uin = -(g*k*Amp/w) * (cosh(k*(yp+h))/cosh(k*h)).*cos(k.*xp-w*t) ; % Need to check the multiple fluxes value at the very last element
%             uL = -(g*k*Amp/w)*cos(k*x(1)-w*t); %uout = 0;
            uL = -(g*k*Amp/w) * (cosh(k*(yf(1)+h))/cosh(k*h)).*cos(k.*x(1)-w*t);
            uout = 0%(g*k*Amp/w) * (cosh(k*(yd+h))/cosh(k*h)).*cos(k.*xd-w*t) ;
% pot = 5; uin = 0; 
% uout = 0;
% uL = 0
[ktvec,knvec]=bcs(xf, xd, xb, xp,pot,uout,uin);
[C1, C2, C3, C4] = corners(xf,xd,xb,xp);
nnode = length(x);  nelem  = ((length(x))/2);
corner4 = [C1 C2 C3 C4];

figure(1)
plot(x,y,'.-'); grid on; title('2D Domain'); xlabel('x'); ylabel('y');

%%
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

[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(x,y,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
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
[U,V,NX,NY] = sderivatives(xf,yf,x,y,p,dpdn,qnode);
t = t0;
% clc
figure;
sgtitle('Seabed potential')

% subplot(1,2,1)
plot(xb,p(C3:C4),'-x')
grid on; hold on
xlabel('bottom nodes')
% subplot(1,2,2)
grid on
potb = (g*Amp/w)*(cosh(k*(yb+h))/cosh(k*h)).*sin(k.*xb-w*t);
plot(xb,potb,'o-k')
xlabel('bottom nodes')
grid on
legend('BEM','Airy')


ua = (g*k*Amp/w) * (cosh(k.*(yf+h))/cosh(k*h)).*cos(k.*xf-w*t);
va = (g*k*Amp/w) * (sinh(k.*(yf+h))/cosh(k*h)).*sin(k.*xf-w*t);
figure;
sgtitle('Free Surface Velocities')
% subplot(1,2,1)
title('BEM')
plot(xf,U,'o-'); hold on; grid 
plot(xf,ua,'-x');
legend('U','ua')
xlabel('fs nodes')
% subplot(1,2,2)
 figure;
title('Airy')
plot(xf,V,'o-b');hold on; grid on
plot(xf,va,'x-r')

legend('V','va')
xlabel('fs nodes')
% legend('U','ua','V','va')

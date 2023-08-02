%% New Time Marching code
clc
close all
clear

tic
g = 9.81;
fx =130*1;          % FS nodes
h = 1;   
Amp = 0.05;         % Amplitude in meters
 
T = 4.17; 
Amp = 0.035;
T  = 1.1318;
w = 2*pi/T;  % Wave period int Depth

lambda = 12.57;
% [lambda]=dispersioncalc(w,h)
lambda = 2;
h = 2
eph =0 ;
% Amp = [0.0111    0.0124    0.0136    0.0146    0.0157    0.0192]; 
% T = [4.5152    4.0386    3.6867    3.4132    3.1928    2.6069];
w = 2*pi./T;  
% eph = [0 pi/8 pi/4 pi/2 pi/3 pi/6];% Wave period Deepwater CASE 1
% lambda = [13.6756   12.1272   10.9751   10.0725    9.3392    7.3557];
           % Tank depth
waveComp = length(Amp);




k = 2*pi./lambda

requiredAmp = 0.142*tanh (k*h)*lambda/2

H_L = 2*Amp./lambda           % Steepness
ka = k.*Amp                   % KA
kh = k.*h
disp('d/L > 0.5?'); h/max(lambda)
L = 8*max(lambda);                       % Tank Length
t0 = 0; tend =12*max(T);           % Simulation timing
dt = min(T)/50                   % Time Step
% m = 1;


% Amp = 0.01274;         % Amplitude in meters
% T = 0.9363; w = 2*pi/T;  % Wave period
% h = 0.1; 

timeaxis = t0:dt:tend;
steps = 1:length(timeaxis);
airy = zeros(1,length(steps));
stepi = 1;
mely = airy;
%% Preprocessing    
[x,y,xf,xd,xb,xp,yf,yd,yb,yp] = tanknwt(fx,h,L,w,k,Amp); %dx changed
[C1, C2, C3, C4] = corners(xf,xd,xb,xp);
nnode = length(x);  nelem  = ((length(x))/2);
corner4 = [C1 C2 C3 C4];

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
Handle1 = fill(x,y,'.-'); hold on; grid on; title('2D NWT');  axis([0 L -h 2*h])
% dtmax = sqrt((8/pi)*( -(xf(1)-xf(2)) /g))

%% Damping zone definition

n=2*max(lambda);     % Damping zone length in meters as a function of wavelength
[val,idx]=min(abs(xf-n));
minVal=xf(idx); % Finds close enough values then looks for the index
xpos=L-xf(idx); % x0
alphad = 1; % alpha;
dzone = (fx+1)-idx;
nu = zeros(1,length(xf));
nu(dzone:end) = alphad*max(w)*((xf(dzone:end)-xpos)/max(lambda)).^2;

counter = 0;
% m = 1;
pot = 0.*xf; yf = 0.*xf; eta = yf;
%%
yp0 = yp;
while t0<tend

[r_t] = rampfunction(t0,min(T));
%% Step 1
[yp] = semiLateral(yp0,t0,g,Amp,w,h,k,r_t);
    yint = y; % intermediate Y
    xint = x;
%     uin = zeros(1,length(xp)); uL = 0;
%     for m = 1:2
%     uin = -r_t*( uin+(g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yp+h))/cosh(k(m)*h)).*cos(k(m).*xp-w(m)*t0+eph(m))  ); % Need to check the multiple fluxes value at the very last element
%      uL = -r_t*( uL+(g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yf(1)+h))/cosh(k(m)*h)).*cos(k(m).*0-w(m)*t0+eph(m))  );% -r_t*(g*k*Amp/w) * cos(k*x(1)-w*t0);
%     end
    uout = 0;
   
%     uin = -r_t*( (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yp+h))/cosh(k(m)*h)).*cos(k(m).*xp-w(m)*t0+eph(m)) + (g*k(m+1)*Amp(m+1)/w(m+1)) * (cosh(k(m+1)*(yp+h))/cosh(k(m+1)*h)).*cos(k(m+1).*xp-w(m+1)*t0+eph(m+1)) ); % Need to check the multiple fluxes value at the very last element
%      uL = -r_t*( (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yf(1)+h))/cosh(k(m)*h)).*cos(k(m).*0-w(m)*t0+eph(m)) + (g*k(m+1)*Amp(m+1)/w(m+1)) * (cosh(k(m+1)*(yf(1)+h))/cosh(k(m+1)*h)).*cos(k(m+1).*0-w(m+1)*t0+eph(m+1))  );% -r_t*(g*k*Amp/w) * cos(k*x(1)-w*t0);

%%    Inlet Velocity % MNA eta(x= 0) and x = 0 assumed see if that is a problem or not
  [uin, uL] =     StokesIrregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,yf,t0,h); % without interacton
%     [uin,uL] = nonlinearStokesVel(r_t,waveComp,Amp,k,w,eph,xp,yp,yf,t0,h); % with frequency interaction
%%

    mely(stepi) = eta((fx/2)+1);
    stepi = stepi+1;

disp('Time = ');disp(t0);

pot; eta;
[ktvec,knvec]=bcs(xf, xd, xb, xp,pot,uout,uin);
[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(x,yint,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
[U,V,NX,NY,detadx] = sderivatives2(xf,yf,x,yint,p,dpdn,qnode);

%% Semi-Lagrange
% K1p = -g.*eta - 0.5*(U.^2+V.^2) + V.*(V-U.*detadx) -nu.*p(1:fx+1);
% K1e = (V-U.*detadx-nu.*eta);
%% Fully-Lagrange
K1p = -g.*eta + 0.5*(U.^2+V.^2) -nu.*p(1:fx+1);
K1e = (V-nu.*eta);
K1x = U;

%% Step 2
poti = pot+0.5*K1p*dt;
etai = eta+0.5*K1e*dt;
xfi = xf+0.5*K1x*dt;

ti = t0+dt/2;
[yp] = semiLateral(yp0,ti,g,Amp,w,h,k,r_t);
%      uin = zeros(1,length(xp)); uL = 0;
%     for m = 1:2
%     uin = -r_t*( uin+(g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yp+h))/cosh(k(m)*h)).*cos(k(m).*xp-w(m)*ti+eph(m))    ); % Need to check the multiple fluxes value at the very last element
%      uL = -r_t*( uL+(g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(etai(1)+h))/cosh(k(m)*h)).*cos(k(m).*0-w(m)*ti+eph(m))  );% -r_t*(g*k*Amp/w) * cos(k*x(1)-w*t0);
%     end
%      uin = zeros(1,length(xp)); uL = 0;
%     for m = 1:2
%     uin = -r_t*( (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yp+h))/cosh(k(m)*h)).*cos(k(m).*xp-w(m)*ti+eph(m)) + (g*k(m+1)*Amp(m+1)/w(m+1)) * (cosh(k(m+1)*(yp+h))/cosh(k(m+1)*h)).*cos(k(m+1).*xp-w(m+1)*ti+eph(m+1))    ); % Need to check the multiple fluxes value at the very last element
%      uL = -r_t*( (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(etai(1)+h))/cosh(k(m)*h)).*cos(k(m).*0-w(m)*ti+eph(m)) +  (g*k(m+1)*Amp(m+1)/w(m+1)) * (cosh(k(m+1)*(etai(1)+h))/cosh(k(m+1)*h)).*cos(k(m+1).*0-w(m+1)*ti+eph(m+1))  );% -r_t*(g*k*Amp/w) * cos(k*x(1)-w*t0);
%     end
%      uout = (Amp*g*k/w) * (cosh(k.*yd+h)/sinh(k*h)).*cos( k*xd-w*ti);
%%  Inlet Velocity
[uin, uL] = StokesIrregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etai,ti,h); % without interacton
% [uin,uL] = nonlinearStokesVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etai,ti,h); % with frequency interaction
%%
[ktvec,knvec]=bcs(xf, xd, xb, xp,poti,uout,uin);
yint(1:fx+1) = etai; %
xint(1:fx+1) = xfi; %


[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
[U,V,NX,NY,detadx] = sderivatives2(xfi,etai,xint,yint,p,dpdn,qnode);

% K2p = -g.*etai - 0.5*(U.^2+V.^2 ) + V.*(V-U.*detadx)-nu.*p(1:fx+1);
% K2e = (V-U.*detadx-nu.*etai);

K2p = -g.*etai + 0.5*(U.^2+V.^2 ) -nu.*p(1:fx+1);
K2e = (V-nu.*etai);
K2x = U;

%% Step 3
potii = pot+0.5*K2p*dt;
etaii = eta+0.5*K2e*dt;
xfii = xf+0.5*K2x*dt;


tii = t0+dt/2;
[yp] = semiLateral(yp0,tii,g,Amp,w,h,k,r_t);
%  uin = zeros(1,length(xp)); uL = 0;
%     for m = 1:2
%     uin = -r_t*( uin+(g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yp+h))/cosh(k(m)*h)).*cos(k(m).*xp-w(m)*tii+eph(m))    ); % Need to check the multiple fluxes value at the very last element
%      uL = -r_t*( uL+(g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(etaii(1)+h))/cosh(k(m)*h)).*cos(k(m).*0-w(m)*tii+eph(m))  );% -r_t*(g*k*Amp/w) * cos(k*x(1)-w*t0);
%     end
    
%  uin = zeros(1,length(xp)); uL = 0;
%     for m = 1:2
%     uin = -r_t*( (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yp+h))/cosh(k(m)*h)).*cos(k(m).*xp-w(m)*tii+eph(m)) +(g*k(m+1)*Amp(m+1)/w(m+1)) * (cosh(k(m+1)*(yp+h))/cosh(k(m+1)*h)).*cos(k(m+1).*xp-w(m+1)*tii+eph(m+1))    ); % Need to check the multiple fluxes value at the very last element
%      uL = -r_t*( (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(etaii(1)+h))/cosh(k(m)*h)).*cos(k(m).*0-w(m)*tii+eph(m)) + (g*k(m+1)*Amp(m+1)/w(m+1)) * (cosh(k(m+1)*(etaii(1)+h))/cosh(k(m+1)*h)).*cos(k(m+1).*0-w(m+1)*tii+eph(m+1))   );% -r_t*(g*k*Amp/w) * cos(k*x(1)-w*t0);
%     end
%      uout = (Amp*g*k/w) * (cosh(k.*yd+h)/sinh(k*h)).*cos( k*xd-w*tii);
%% Inlet Velocity
[uin, uL] = StokesIrregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etaii,tii,h);% without interacton
% [uin,uL] = nonlinearStokesVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etaii,tii,h); % with frequency interaction
%%
[ktvec,knvec]=bcs(xf, xd, xb, xp,potii,uout,uin);
yint(1:fx+1) = etaii; %
xint(1:fx+1) = xfii; %

[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
[U,V,NX,NY,detadx] = sderivatives2(xfii,etaii,xint,yint,p,dpdn,qnode);

% K3p = -g.*etaii - 0.5*(U.^2+V.^2) + V.*(V-U.*detadx)-nu.*p(1:fx+1);
% K3e = (V-U.*detadx-nu.*etaii);

K3p = -g.*etaii + 0.5*(U.^2+V.^2) -nu.*p(1:fx+1);
K3e = (V-nu.*etaii);
K3x = U;
%% Step 4
potiii = pot+K3p*dt;
etaiii = eta+K3e*dt;
xfiii = xf+K3x*dt;

tiii = t0+dt;

[yp] = semiLateral(yp0,tiii,g,Amp,w,h,k,r_t);
%  uin = zeros(1,length(xp)); uL = 0;
%     for m = 1:2
%     uin = -r_t*( uin+(g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yp+h))/cosh(k(m)*h)).*cos(k(m).*xp-w(m)*tiii+eph(m))   ); % Need to check the multiple fluxes value at the very last element
%      uL = -r_t*( uL+(g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(etaiii(1)+h))/cosh(k(m)*h)).*cos(k(m).*0-w(m)*tiii+eph(m)) );% -r_t*(g*k*Amp/w) * cos(k*x(1)-w*t0);
%     end

%  uin = zeros(1,length(xp)); uL = 0;
%     for m = 1:2
%     uin = -r_t*( (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yp+h))/cosh(k(m)*h)).*cos(k(m).*xp-w(m)*tiii+eph(m))  +(g*k(m+1)*Amp(m+1)/w(m+1)) * (cosh(k(m+1)*(yp+h))/cosh(k(m+1)*h)).*cos(k(m+1).*xp-w(m+1)*tiii+eph(m+1))   ); % Need to check the multiple fluxes value at the very last element
%      uL = -r_t*( (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(etaiii(1)+h))/cosh(k(m)*h)).*cos(k(m).*0-w(m)*tiii+eph(m))  +(g*k(m+1)*Amp(m+1)/w(m+1)) * (cosh(k(m+1)*(etaiii(1)+h))/cosh(k(m+1)*h)).*cos(k(m+1).*0-w(m+1)*tiii+eph(m+1)) );% -r_t*(g*k*Amp/w) * cos(k*x(1)-w*t0);
%     end
%      uout = (Amp*g*k/w) * (cosh(k.*yd+h)/sinh(k*h)).*cos( k*xd-w*tiii);

[uin, uL] = StokesIrregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etaiii,tiii,h); % without interacton
% [uin,uL] = nonlinearStokesVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etaiii,tiii,h); % with frequency interaction

[ktvec,knvec]=bcs(xf, xd, xb, xp,potiii,uout,uin);
yint(1:fx+1) = etaiii; %
xint(1:fx+1) = xfiii; %

[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
[U,V,NX,NY,detadx] = sderivatives2(xfiii,etaiii,xint,yint,p,dpdn,qnode);

K4p = -g.*etaiii + 0.5*(U.^2+V.^2) -nu.*p(1:fx+1);
K4e = (V-nu.*etaiii);
K4x = U;

potN = pot +(dt/6)*(K1p+2*K2p+2*K3p+K4p);
etaN = eta +(dt/6)*(K1e+2*K2e+2*K3e+K4e);
xN = xf +(dt/6)*(K1x+2*K2x+2*K3x+K4x);

% Prepare for next iteration
xN(1) = 0;
% etaN(1) = r_t*( Amp*cos(-w.*tiii)+ (3*k*Amp*Amp/(4*(k*h)^3))*cos(2*((-w.*tiii))) ) ;
% Not working
% potN(1) = r_t*( -Amp*g*w^-1 *cosh(k*( h+etaN(1) ))*cosh(k*h)^-1 * sin(-w*tiii)    -(3/8)*Amp^2*w*cosh(2*k*(h+etaN(1) ) )/sinh(k*h)^4 *sin(2*(-w*tiii)) );
[xN,etaN,potN] = CSregrid(xN,etaN,potN);
if counter == 10
    [etaN] = FreeSurfSmooth5pt(xf,etaN);
    [potN] = FreeSurfSmooth5pt(xf,potN);
    counter = 0;
else
    counter = counter+1;
end

pot = potN;
eta = etaN; y(1:fx+1) = eta;
t0 = tiii;
yf = eta;
xf = xN;
y = [eta yd yb yp];
x = [xf xd xb xp];
set(Handle1,'XData',x,'YData', y); grid on ;  drawnow; 
end

toc
%% Post Process
figure(20);

% % mely4DT = mely
% eta4DT = eta
plot(timeaxis,mely,'.r','LineWidth',1)
grid on; hold on
set(gca,'FontSize',18)
% airyy  = Amp*cos(k.*xf(((fx/2)+1))-w.*timeaxis);
% plot(timeaxis,airyy,'-r')
% legend('MEL','Steep-Airy')
stokey = 0;
for m = 1:waveComp
stokey  = stokey+Amp(m)*cos(k(m).*xf(((fx/2)+1))-w(m).*timeaxis+eph(m) ) + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xf(((fx/2)+1))-w(m).*timeaxis+eph(m) )))  ;
end
plot(timeaxis,stokey,'-k','LineWidth',1)
xprobe = xf(((fx/2)+1));
% [TSetaS] = NLTSelevationStokes(waveComp,Amp,k,w,eph,xprobe,timeaxis,h)
% plot(timeaxis,TSetaS,'-k','LineWidth',0.1)

legend('MEL','Theoretical','Theoretical Interaction')
% legend('MEL','Theoretical Interaction')
%%

title('Comparing \eta at mid tank time-history')
ylabel('\eta (m) ')
xlabel('t (s)')
%%
figure(30);

plot(xf,eta,'.-r','LineWidth',1)
grid on; hold on
set(gca,'FontSize',18)
% etaairyy  = Amp*cos(k.*xf-w.*tend);
% plot(xf,etaairyy,'--r')
% legend('MEL','Steep-Airy')
hold on;
%     plot(xf(dzone)./L,0,'rd')
etastokey = 0;
for m = 1:waveComp
etastokey  = etastokey+Amp(m)*cos(k(m).*xf-w(m).*tiii+eph(m)) + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xf-w(m).*tiii+eph(m)))) ;
end
plot(xf,etastokey,'--k','LineWidth',0.5)

[etaS] = elevationStokes(waveComp,Amp,k,w,eph,xf,tiii,h)
% plot(xf,etaS,'-k','LineWidth',0.8)

legend('MEL','Theoretical','Theoretical Interaction')
% legend('MEL','Theoretical Interaction')
% axis square
%%
str = sprintf('Comparing elevation at %d seconds', tend)
title(str)

ylabel('\eta (m)')
xlabel('x (m)')
%%
figure(100)
plot(xf,pot)
hold on;
grid on
xlabel('Velocity Potential')
ylabel('Velocity Potential')
xlabel('x')
% 
% %  umean = mean(U(1:dzone))
%  
% % Ueul= -Amp*Amp*k*pi*coth(k*h)/(T*k*h)
% 
% 
% 

%% Post Process

figure(25);
% % mely4DT = mely
% eta4DT = eta
plot(timeaxis,mely./sum(Amp),'.')
grid on; hold on
set(gca,'FontSize',18)
% airyy  = Amp*cos(k.*xf(((fx/2)+1))-w.*timeaxis);
% plot(timeaxis,airyy,'-r')
% legend('MEL','Steep-Airy')
stokey = 0;
for m = 1:waveComp
stokey  = stokey +Amp(m)*cos(k(m).*xf(((fx/2)+1))-w(m).*timeaxis+eph(m))  + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xf(((fx/2)+1))-w(m).*timeaxis+eph(m))))  ;
end
plot(timeaxis,stokey./sum(Amp),'-r')
legend('MEL','Theoretical')
%%

title('Comparing \eta at mid tank time-history')
ylabel('\eta (m)')
xlabel('t (s)')
%%
figure(35);

plot(xf,eta./sum(Amp))
grid on; hold on
set(gca,'FontSize',18)
% etaairyy  = Amp*cos(k.*xf-w.*tend);
% plot(xf,etaairyy,'--r')
% legend('MEL','Steep-Airy')
hold on;
%     plot(xf(dzone)./L,0,'rd')
etastokey = 0;
for m = 1:waveComp
etastokey  = etastokey+Amp(m)*cos(k(m).*xf-w(m).*tiii+eph(m))+ (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xf-w(m).*tiii+eph(m))))  ;
end
plot(xf,etastokey./sum(Amp),'--r')
legend('MEL','Theoretical')
% axis square
%%
str = sprintf('Comparing elevation at %d seconds', tend)
title(str)

ylabel('\eta (m)')
xlabel('x (m)')

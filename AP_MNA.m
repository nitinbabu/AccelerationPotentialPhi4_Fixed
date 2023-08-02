%% New Time Marching code
clc
close all
clear

% tic
% g = 9.81;
% fx =91*1;          % FS nodes
% wpoints = 41;
% Amp = 0.05;         % Amplitude in meters
% T = 4.17; w = 2*pi/T;  % Wave period int Depth
% lambda = 12.57;


caseno = 1:12;
Ampk = (0.07/2)+zeros(1,length(caseno));
Tcase = [0.800300000000000,0.876700000000000,0.946900000000000,1.01230000000000,1.07370000000000,1.13180000000000,1.18700000000000,1.23980000000000,1.33920000000000,1.49720000000000,1.56010000000000,1.64010000000000];
lambdacase = [1,1.20000000000000,1.40000000000000,1.60000000000000,1.80000000000000,2,2.20000000000000,2.40000000000000,2.80000000000000,3.50000000000000,3.80000000000000,4.20000000000000];
zetacase = [1.57080000000000,1.30900000000000,1.12200000000000,0.981700000000000,0.872700000000000,0.785400000000000,0.714000000000000,0.654500000000000,0.561000000000000,0.448800000000000,0.413400000000000,0.374000000000000];
    g = 9.81;
    fx =91*1;          % FS nodes
    wpoints = 41;
% Amp = 0.01;         % Amplitude in meters
% T = 0.8005; w = 2*pi/T;  % Wave period int Depth
% lambda = 1.0005;
%% No.01 MDPI paper
% Amp = 0.035;         % Amplitude in meters
% T = 0.8003; w = 2*pi/T;  % Wave period int Depth
% lambda = 1.0;
% %% No.06 MDPI paper
% Amp = 0.035;         % Amplitude in meters
% T = 1.1318; w = 2*pi/T;  % Wave period int Depth
% lambda = 2.0;
% 
% %% No.12 MDPI paper
% Amp = 0.035;         % Amplitude in meters
% T = 1.6401; w = 2*pi/T;  % Wave period int Depth
% lambda = 4.2;

for casenoi = 1:12
    tic
   
    
    
    Amp = Ampk(casenoi);
    T = Tcase(casenoi);w = 2*pi/T;  % Wave period int Depth
    lambda = lambdacase(casenoi);
    

zetafreq =  w^2*0.5/2/g;
h = lambda;
eph =0 ;
waveComp = length(Amp);
  
k = 2*pi./lambda;
H_L = 2*Amp./lambda ;          % Steepness
ka = k.*Amp  ;                 % KA
kh = k.*h;
disp('d/L > 0.5?'); h/max(lambda);
L = 8*max(lambda);                       % Tank Length
% t0 = 0; tend =15*T;           % Simulation timing
t0 = 0; 
if casenoi>=6
    
tend =20*T;           % Simulation timing
elseif casenoi < 6 && casenoi > 3
        tend = 15*T;
else
        tend = 10*T;
end

dt = min(T)/60                  % Time Step


Vn = 0; % Body Boundary Condition

timeaxis = t0:dt:tend;
steps = 1:length(timeaxis);
% Fx = steps; Fz = Fx; My = Fz;
airy = zeros(1,length(steps));
stepi = 1;
mely = airy; melyC1 = mely; melyC3 = mely; melyC2 = mely;
%% Preprocessing    
% [x,y,xf,xd,xb,xp,yf,yd,yb,yp] = t anknwt(fx,h,L,w,k,Amp); %dx changed
% BL = L/2-0.2*L/2; BR = L/2+0.2*L/2;


[x,y,xfl, nodexbody,xfr,xd,xb,xp,yfl,nodeybody,yfr,yd,yb,yp,cx,cy,xbody,ybody,wpoints] = INITnwtbody(fx,h,L,w,k,Amp,wpoints);
[C1, C2, C3, C4, C5, C6] = bodycorners(xfl,nodexbody,xfr,xd,xb,xp);
nnode = length(x);  nelem  = ((length(x))/2);
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
% plot(x,y,'.-k')
Handle1 = fill(x,y,'.-'); hold on; grid on; title('2D NWT');  axis([0 L -h 2*h]);
% axis([6 10 -0.5 0.5]);

%% Damping zone definition


% Inlet Damping Zone-
nl=2*max(lambda);     % Damping zone length in meters as a function of wavelength
[val,idxl]=min(abs(xfl-nl));
minVal=xfl(idxl); % Finds close enough values then looks for the index
xposl=xfl(idxl); % x0
alphadl = 1; % alpha;
dzonel = length(xfl)-idxl;
nul = zeros(1,length(xfl));
nul(1:dzonel) = alphadl*max(w)*((xfl(1:dzonel)-xposl)/max(lambda)).^2;
plot(xposl,yfl(dzonel),'rd')



%Damping Zone Right-Downstream
n=L-2*max(lambda);     % Damping zone length in meters as a function of wavelength
[val,idx]=min(abs(xfr-n));
minVal=xfr(idx); % Finds close enough values then looks for the index
xpos=xfr(idx); % x0
alphad = 1; % alpha;
dzone = length(xfr)-idx;
nu = zeros(1,length(xfr));
nu(dzone:end) = alphad*max(w)*((xfr(dzone:end)-xpos)/max(lambda)).^2;
plot(xpos,yfr(dzone),'rd')


counter = 0;
potL = 0.*xfl; potR = 0.*xfr; lyfl = 0.*xfl;  yfr = 0.*xfr; etaL = yfl; etaR = yfr;

%%
yp0 = yp;
while t0<tend

[r_t] = rampfunction(t0,min(T));
%% Step 1
[yp] = semiLateral(yp0,t0,g,Amp,w,h,k,r_t);
    yint = y; % intermediate Y
    xint = x;

    uout = 0;
   
%%    Inlet Velocity % MNA eta(x= 0) and x = 0 assumed see if that is a problem or not
  [uin, uL] =     StokesIrregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,yfl,t0,h); % without interacton
%     [uin,uL] = nonlinearStokesVel(r_t,waveComp,Amp,k,w,eph,xp,yp,yf,t0,h); % with frequency interaction


%    mely(stepi) = etaL(21);
%     melyC1(stepi) = etaL(21);   xfl(21);
%     melyC2(stepi) = etaL(93-6);  xfl(93-6);
%     melyC3(stepi) = etaR(7);   xfr(7);
   

str = sprintf('Time =  %.3g T seconds', t0/T);
disp(str)
% pot; eta;
% [ktvec,knvec]=bcs(xf, xd, xb, xp,pot,uout,uin);
[ktvec,knvec]=bcsB( xfl,nodexbody,xfr, xd, xb, xp,potL,Vn,potR,uout,uin);

% [p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(x,yint,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEMBODY(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner6,uL);

% [U,V,NX,NY,detadx] = sderivatives2(xf,yf,x,yint,p,dpdn,qnode);
 [UL,VL,NXL,NYL,detadxL] = sderivatives2LEFT(xfl,yfl,xint,yint,p,dpdn,qnode,corner6);
 [Ub,Vb,NXb,NYb,detadxb] = sderivatives2BODY(nodexbody,nodeybody,xint,yint,p,dpdn,qnode,corner6);
 [UR,VR,NXR,NYR,detadxR] = sderivatives2RIGHT(xfr,yfr,xint,yint,p,dpdn,qnode,corner6);

 
%% Semi-Lagrange
% K1p = -g.*eta - 0.5*(U.^2+V.^2) + V.*(V-U.*detadx) -nu.*p(1:fx+1);
% K1e = (V-U.*detadx-nu.*eta);
%% Fully-Lagrange
pa = 0;
for m = 1:waveComp
pa  = pa +r_t*(g*Amp(m))/w(m).*(cosh(k(m).*(yfl+h))/cosh(k(m)*h)).*sin(k(m).*xfl-w(m).*t0+eph(m)) + (3/8*Amp(m)*Amp(m)*w(m))*( cosh(2*k(m).*(yfl+h) )/sinh(k(m)*h)^4 ).*sin(2*((k(m).*xfl-w(m).*t0+eph(m)))) ;
end

etaanalytical = 0;
for m = 1:waveComp
etaanalytical  = etaanalytical+r_t*Amp(m)*cos(k(m).*xfl-w(m).*t0+eph(m)) + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xfl-w(m).*t0+eph(m)))) ;
end

K1pL = -g.*etaL + 0.5*(UL.^2+VL.^2) -nul.*(p(C1:C2)-pa);
K1eL = (VL-nul.*(yfl-etaanalytical));
K1xL = UL;


K1pR = -g.*etaR + 0.5*(UR.^2+VR.^2) -nu.*p(C3:C4);
K1eR = (VR-nu.*etaR);
K1xR = UR;

%% Step 2
potiL = potL+0.5*K1pL*dt;
etaiL = etaL+0.5*K1eL*dt;
xfil = xfl+0.5*K1xL*dt;


potiR = potR+0.5*K1pR*dt;
etaiR = etaR+0.5*K1eR*dt;
xfir = xfr+0.5*K1xR*dt;

ti = t0+dt/2;
[yp] = semiLateral(yp0,ti,g,Amp,w,h,k,r_t);

%%  Inlet Velocity
[uin, uL] = StokesIrregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etaiL,ti,h); % without interacton
% [uin,uL] = nonlinearStokesVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etai,ti,h); % with frequency interaction
%%
% [ktvec,knvec]=bcs(xf, xd, xb, xp,poti,uout,uin);
% yint(1:fx+1) = etai; %
% xint(1:fx+1) = xfi; %
% 
% 
% [p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
% [U,V,NX,NY,detadx] = sderivatives2(xfi,etai,xint,yint,p,dpdn,qnode);


[ktvec,knvec]=bcsB( xfil,nodexbody,xfir, xd, xb, xp,potiL,Vn,potiR,uout,uin);
 yint(C1:C2) = etaiL;  xint(C1:C2) = xfil; %
 yint(C3:C4) = etaiR;  xint(C3:C4) = xfir; %
 
[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEMBODY(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner6,uL);
[UL,VL,NXL,NYL,detadxL] = sderivatives2LEFT(xfil,etaiL,xint,yint,p,dpdn,qnode,corner6);
 [Ub,Vb,NXb,NYb,detadxb] = sderivatives2BODY(nodexbody,nodeybody,xint,yint,p,dpdn,qnode,corner6);
[UR,VR,NXR,NYR,detadxR] = sderivatives2RIGHT(xfir,etaiR,xint,yint,p,dpdn,qnode,corner6);

pa = 0;
for m = 1:waveComp
pa  = pa +r_t*(g*Amp(m))/w(m).*(cosh(k(m).*(etaiL+h))/cosh(k(m)*h)).*sin(k(m).*xfil-w(m).*ti+eph(m)) + (3/8*Amp(m)*Amp(m)*w(m))*( cosh(2*k(m).*(etaiL+h) )/sinh(k(m)*h)^4 ).*sin(2*((k(m).*xfil-w(m).*ti+eph(m)))) ;
end

etaanalytical = 0;
for m = 1:waveComp
etaanalytical  = etaanalytical+r_t*Amp(m)*cos(k(m).*xfil-w(m).*ti+eph(m)) + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xfil-w(m).*ti+eph(m)))) ;
end




K2pL = -g.*etaiL + 0.5*(UL.^2+VL.^2 ) -nul.*(p(C1:C2)-pa);
K2eL = (VL-nul.*(etaiL-etaanalytical));
K2xL = UL;

K2pR = -g.*etaiR + 0.5*(UR.^2+VR.^2 ) -nu.*p(C3:C4);
K2eR = (VR-nu.*etaiR);
K2xR = UR;

%% Step 3

potiiL = potL+0.5*K2pL*dt;
etaiiL = etaL+0.5*K2eL*dt;
xfiil = xfl+0.5*K2xL*dt;

potiiR = potR+0.5*K2pR*dt;
etaiiR = etaR+0.5*K2eR*dt;
xfiir = xfr+0.5*K2xR*dt;


tii = t0+dt/2;
[yp] = semiLateral(yp0,tii,g,Amp,w,h,k,r_t);
%% Inlet Velocity
[uin, uL] = StokesIrregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etaiiL,tii,h);% without interacton
% [uin,uL] = nonlinearStokesVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etaii,tii,h); % with frequency interaction
%% OLD
% [ktvec,knvec]=bcs(xf, xd, xb, xp,potii,uout,uin);
% yint(1:fx+1) = etaii; %
% xint(1:fx+1) = xfii; %
% 
% [p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
% [U,V,NX,NY,detadx] = sderivatives2(xfii,etaii,xint,yint,p,dpdn,qnode);
% 
% % K3p = -g.*etaii - 0.5*(U.^2+V.^2) + V.*(V-U.*detadx)-nu.*p(1:fx+1);
% % K3e = (V-U.*detadx-nu.*etaii);
% 
% K3p = -g.*etaii + 0.5*(U.^2+V.^2) -nu.*p(1:fx+1);
% K3e = (V-nu.*etaii);
% K3x = U;
% NEW

[ktvec,knvec]=bcsB( xfiil,nodexbody,xfiir, xd, xb, xp,potiiL,Vn,potiiR,uout,uin);
 yint(C1:C2) = etaiiL;  xint(C1:C2) = xfiil; %
 yint(C3:C4) = etaiiR;  xint(C3:C4) = xfiir; %
 
[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEMBODY(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner6,uL);
[UL,VL,NXL,NYL,detadxL] = sderivatives2LEFT(xfiil,etaiiL,xint,yint,p,dpdn,qnode,corner6);
 [Ub,Vb,NXb,NYb,detadxb] = sderivatives2BODY(nodexbody,nodeybody,xint,yint,p,dpdn,qnode,corner6);
[UR,VR,NXR,NYR,detadxR] = sderivatives2RIGHT(xfiir,etaiiR,xint,yint,p,dpdn,qnode,corner6);


pa = 0;
for m = 1:waveComp
pa  = pa +r_t*(g*Amp(m))/w(m).*(cosh(k(m).*(etaiiL+h))/cosh(k(m)*h)).*sin(k(m).*xfiil-w(m).*tii+eph(m)) + (3/8*Amp(m)*Amp(m)*w(m))*( cosh(2*k(m).*(etaiiL+h) )/sinh(k(m)*h)^4 ).*sin(2*((k(m).*xfiil-w(m).*tii+eph(m)))) ;
end

etaanalytical = 0;
for m = 1:waveComp
etaanalytical  = etaanalytical+r_t*Amp(m)*cos(k(m).*xfiil-w(m).*tii+eph(m)) + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xfiil-w(m).*tii+eph(m)))) ;
end

K3pL = -g.*etaiiL + 0.5*(UL.^2+VL.^2 ) -nul.*(p(C1:C2)-pa);
K3eL = (VL-nul.*(etaiiL-etaanalytical));
K3xL = UL;

K3pR = -g.*etaiiR + 0.5*(UR.^2+VR.^2 ) -nu.*p(C3:C4);
K3eR = (VR-nu.*etaiiR);
K3xR = UR;

%% Step 4

potiiiL = potL+K3pL*dt;
etaiiiL = etaL+K3eL*dt;
xfiiil = xfl+K3xL*dt;

potiiiR = potR+K3pR*dt;
etaiiiR = etaR+K3eR*dt;
xfiiir = xfr+K3xR*dt;

tiii = t0+dt;

[yp] = semiLateral(yp0,tiii,g,Amp,w,h,k,r_t);
%%
[uin, uL] = StokesIrregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etaiiiL,tiii,h); % without interacton
% [uin,uL] = nonlinearStokesVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etaiii,tiii,h); % with frequency interaction
%%
% [ktvec,knvec]=bcs(xf, xd, xb, xp,potiii,uout,uin);
% yint(1:fx+1) = etaiii; %
% xint(1:fx+1) = xfiii; %
% 
% [p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
% [U,V,NX,NY,detadx] = sderivatives2(xfiii,etaiii,xint,yint,p,dpdn,qnode);
% 
% K4p = -g.*etaiii + 0.5*(U.^2+V.^2) -nu.*p(1:fx+1);
% K4e = (V-nu.*etaiii);
% K4x = U;

[ktvec,knvec]=bcsB( xfiiil,nodexbody,xfiiir, xd, xb, xp,potiiiL,Vn,potiiiR,uout,uin);
 yint(C1:C2) = etaiiiL;  xint(C1:C2) = xfiiil; %
 yint(C3:C4) = etaiiiR;  xint(C3:C4) = xfiiir; %
 
[p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEMBODY(xint,yint,nnode,nelem,qnode,ktvec,knvec,corner6,uL);
[UL,VL,NXL,NYL,detadxL] = sderivatives2LEFT(xfiiil,etaiiiL,xint,yint,p,dpdn,qnode,corner6);
 [Ub,Vb,NXb,NYb,detadxb] = sderivatives2BODY(nodexbody,nodeybody,xint,yint,p,dpdn,qnode,corner6);
[UR,VR,NXR,NYR,detadxR] = sderivatives2RIGHT(xfiiir,etaiiiR,xint,yint,p,dpdn,qnode,corner6);

pa = 0;
for m = 1:waveComp
pa  = pa +r_t*(g*Amp(m))/w(m).*(cosh(k(m).*(etaiiiL+h))/cosh(k(m)*h)).*sin(k(m).*xfiiil-w(m).*tiii+eph(m)) + (3/8*Amp(m)*Amp(m)*w(m))*( cosh(2*k(m).*(etaiiiL+h) )/sinh(k(m)*h)^4 ).*sin(2*((k(m).*xfiiil-w(m).*tiii+eph(m)))) ;
end

etaanalytical = 0;
for m = 1:waveComp
etaanalytical  = etaanalytical+r_t*Amp(m)*cos(k(m).*xfiiil-w(m).*tiii+eph(m)) + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xfiiil-w(m).*tiii+eph(m)))) ;
end


K4pL = -g.*etaiiL + 0.5*(UL.^2+VL.^2 ) -nul.*(p(C1:C2)-pa);
K4eL = (VL-nul.*(etaiiiL-etaanalytical));
K4xL = UL;

K4pR = -g.*etaiiR + 0.5*(UR.^2+VR.^2 ) -nu.*p(C3:C4);
K4eR = (VR-nu.*etaiiR);
K4xR = UR;


%%


potNL = potL +(dt/6)*(K1pL+2*K2pL+2*K3pL+K4pL);
etaNL = etaL +(dt/6)*(K1eL+2*K2eL+2*K3eL+K4eL);
xNl = xfl +(dt/6)*(K1xL+2*K2xL+2*K3xL+K4xL);


potNR = potR +(dt/6)*(K1pR+2*K2pR+2*K3pR+K4pR);
etaNR = etaR +(dt/6)*(K1eR+2*K2eR+2*K3eR+K4eR);
xNr = xfr +(dt/6)*(K1xR+2*K2xR+2*K3xR+K4xR);
%%
% Prepare for next iteration
xNl(1) = 0;
xNl(end)= x(C2);

xNr(1) = x(C3);
xNr(end) = x(C4);
% figure;
% plot(xNl,etaNL);
% hold on
% plot(xNr,etaNR)


% etaNL(1) = r_t*( Amp*cos(-w.*tiii)+ (3*k*Amp*Amp/(4*(k*h)^3))*cos(2*((-w.*tiii))) ) ;
% Not working
% potN(1) = r_t*( -Amp*g*w^-1 *cosh(k*( h+etaN(1) ))*cosh(k*h)^-1 * sin(-w*tiii)    -(3/8)*Amp^2*w*cosh(2*k*(h+etaN(1) ) )/sinh(k*h)^4 *sin(2*(-w*tiii)) );

[xNl,etaNL,potNL] = CSregrid(xNl,etaNL,potNL);
[xNr,etaNR,potNR] = CSregrid(xNr,etaNR,potNR);

% dddfg
% figure;
% plot(xNl,etaNL);
% hold on
% plot(xNr,etaNR)

if counter == 10
    [etaNL] = FreeSurfSmooth5pt(xNl,etaNL);
    [potNL] = FreeSurfSmooth5pt(xNl,potNL);
    
    [etaNR] = FreeSurfSmooth5pt(xNr,etaNR);
    [potNR] = FreeSurfSmooth5pt(xNr,potNR);
    
    counter = 0;
else
    counter = counter+1;
end
% sss
potL = potNL;
etaL = etaNL; y(C1:C2) = etaL;

potR = potNR;
etaR = etaNR; y(C3:C4) = etaR;



t0 = tiii;
yfl = etaL;
xfl = xNl;

yfr = etaR;
xfr = xNr;

roll = 0; sway = 0; heave = 0; %% Manipulates the body
[x,y,xfl, nodexbody,xfr,xd,xb,xp,yfl,nodeybody,yfr,yd,yb,yp,xbody,ybody] = updatebodymesh(cx,cy,xfl,yfl,xbody,ybody,xfr,yfr,roll,sway,heave,xd,xb,xp,yd,yb,yp,wpoints);
[~,~,p(C2+1:C3-1)] = CSregrid(nodexbody,nodeybody,p(C2+1:C3-1)); %%

y = [etaL nodeybody etaR yd yb yp];
x = [xfl nodexbody xfr xd xb xp];
set(Handle1,'XData',x,'YData', y); grid on ;  drawnow; 


%% Acceleration field
%% Diffraction
[uin, uL] = StokesIrregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etaL,t0,h); % without interacton
uout = 0;
[ktvecdiff,knvecdiff]=bcsB( xfl,nodexbody,xfr, xd, xb, xp,potL,Vn,potR,uout,uin);
[pdiff,dpdndiff,~,~,~,~]=HOBEMBODY(x,y,nnode,nelem,qnode,ktvecdiff,knvecdiff,corner6,uL);

[ULdiff,VLdiff,~,~,~] = sderivatives2LEFT(xfl,etaL,x,y,pdiff,dpdndiff,qnode,corner6);

bodynodexbody = [xfl(end) nodexbody xfr(1)];
bodynodeybody = [yfl(end) nodeybody yfr(1)];

[Ubdiff,Vbdiff,NXbdiff,NYbdiff,detadx,d2pxi2,d2pxin,dpxi] = sderivatives2BODYQ( bodynodexbody, bodynodeybody ,x,y,pdiff,dpdndiff,qnode,corner6);

[URdiff,VRdiff,~,~,~] = sderivatives2RIGHT(xfr,etaR,x,y,pdiff,dpdndiff,qnode,corner6);


[ain, aL] = StokesIrregular_inletACC(r_t,waveComp,Amp,k,w,eph,xp,yp,yfl,t0,h); aout = 0;


potLDiffr = -g*etaL-0.5*(ULdiff.^2+VLdiff.^2);
% potRDiffr = -g*etaL-0.5*(URdiff.^2+VRdiff.^2); %% error
potRDiffr = -g*etaR-0.5*(URdiff.^2+VRdiff.^2); %% fix


qBod = dpdn(C2:C3).*d2pxi2-dpxi.*d2pxin;
[ktvecAccDiff,knvecAccDiff]=bcsB( xfl,nodexbody,xfr, xd, xb, xp,potLDiffr,qBod(2:end-1),potRDiffr,aout,ain);
[ptAccDiff,dpdntAccDiff,~,~,~,~]=HOBEMBODY(x,y,nnode,nelem,qnode,ktvecAccDiff,knvecAccDiff,corner6,aL);

phi_t4 = ptAccDiff(C2:C3);
%%  ptAccDiff % Phi_4 to get phi_t
for it = 2:length(bodynodexbody)
ds(it) = sqrt( (bodynodexbody(it-1)-bodynodexbody(it))^2 + (bodynodeybody(it-1)-bodynodeybody(it))^2  );

end

for it = 1:length(bodynodexbody)
rx_c(it) = (bodynodexbody(it)-L/2);
rz_c(it) = (bodynodeybody(it)-(-0.115)); %Changed gravity G ref. Fig 1 Huang et. al. 2008
end

sumpot4x = 0; 
sumpot4z = 0;
sumpot4m = 0;

% rx-rz is the position vector from C.G. to any location
for kcount =1:length(bodynodexbody)
   sumpot4x = sumpot4x+ -1000*(phi_t4(kcount)+g*bodynodeybody(kcount)+0.5*(Ubdiff(kcount)^2+Vbdiff(kcount)^2) )*NXbdiff(kcount)*ds(kcount);
   sumpot4z = sumpot4z+ -1000*(phi_t4(kcount)+g*bodynodeybody(kcount)+0.5*(Ubdiff(kcount)^2+Vbdiff(kcount)^2) )*NYbdiff(kcount)*ds(kcount);
   sumpot4m = sumpot4m+ -1000*(phi_t4(kcount)+g*bodynodeybody(kcount)+0.5*(Ubdiff(kcount)^2+Vbdiff(kcount)^2) )*(rx_c(kcount)*NYbdiff(kcount)  - rz_c(kcount)*NXbdiff(kcount))*ds(kcount);
   
end

Fx(stepi) = sumpot4x;
Fz(stepi) = sumpot4z;
My(stepi) = sumpot4m;
 stepi = stepi+1;
% Fz
% Fx
end

%%
% figure;
% plot(x(C2:C3),ptAccDiff(C2:C3))
% hold on
% plot(x(C2+1:C3-1),NXbdiff)
%%

toc

% %%
% 
% 
% 
% %     disp(Fx)
%     ampFx = abs(fft(Fx(220:end),2^10));
% 
% %     disp(Fz)
%     ampFz = abs(fft(Fz(220:end),2^10));
% 
% 
% % figure(66)
% % plot(ampFx)
% normxforce =    1000*9.81*0.25*Amp;
% % figure(67)
% % plot(ampFz)
% normzforce =    1000*9.81*0.5*Amp;
% 
% 
% figure(65);
% sgtitle('Forces on Surface Piercing body B=0.5m d = 0.25m, \lambda = 2m')
% subplot(3,1,1)
% title('Horizontal Force')
% plot(timeaxis(1:end-1), Fx,'r-^');
% hold on
% legend('Fx')
% grid minor
% ylabel('Force-x')
% xlabel('t(s)')
% subplot(3,1,2)
% 
% title('Vertical Force')
% plot(timeaxis(1:end-1), Fz,'b-^');
% legend('Fz')
% hold on
% grid minor
% ylabel('Force-z')
% xlabel('t(s)')
% %
% 
% subplot(3,1,3)
% 
% title('Moment')
% plot(timeaxis(1:end-1), My,'b-^');
% legend('My')
% hold on
% grid minor
% ylabel('Moment-y')
% xlabel('t(s)')
% 
% 
% 
% Fs = 1/dt;
% figure(66);
% subplot(1,2,2)
% % Step 2: Calculate the FFT
% fft_dataz = fft(Fz(400:end)./normzforce);
% % Step 3: Compute the amplitudes and frequencies
% N = length(fft_dataz);
% amplitudes = 2/N * abs(fft_dataz(1:N/2+1)); % Compute the one-sided amplitudes
% frequencies = (0:N/2) * Fs / N; % Compute the corresponding frequencies
% 
% % Plot the amplitude spectrum
% % title('')
% plot(frequencies(2:end), amplitudes(2:end));
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% title('Amplitude Spectrum');
% 
% legend('Fz')
% % Note: The amplitudes are scaled by 2/N to account for the one-sided amplitude spectrum.
% % The amplitudes are doubled to represent the energy in both the positive and negative frequency components.
% 
% % You can also find the dominant frequency (peak frequency) using the following:
% [max_amplitude, index] = max(amplitudes);
% dominant_frequency = frequencies(index);
% disp(['Dominant frequency: ', num2str(dominant_frequency), ' Hz']);
% 
% subplot(1,2,1)  %% X force
% 
% % Step 2: Calculate the FFT
% fft_dataX = fft(Fx(400:end)./normxforce);
% % Step 3: Compute the amplitudes and frequencies
% N = length(fft_dataX);
% amplitudes = 2/N * abs(fft_dataX(1:N/2+1)); % Compute the one-sided amplitudes
% frequencies = (0:N/2) * Fs / N; % Compute the corresponding frequencies
% 
% % Plot the amplitude spectrum
% 
% plot(frequencies(2:end), amplitudes(2:end));
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% title('Amplitude Spectrum');
% legend('Fx')
% % Note: The amplitudes are scaled by 2/N to account for the one-sided amplitude spectrum.
% % The amplitudes are doubled to represent the energy in both the positive and negative frequency components.
% 
% % You can also find the dominant frequency (peak frequency) using the following:
% [max_amplitude, index] = max(amplitudes);
% dominant_frequency = frequencies(index);
% disp(['Dominant frequency: ', num2str(dominant_frequency), ' Hz']);
% 
% 
% %%
% 
% 
% %%
% figure;
% quiver(bodynodexbody, bodynodeybody, NXbdiff, NYbdiff)
% hold on
% plot(bodynodexbody, bodynodeybody,'o-r')
% axis equal
% grid minor

% scatter(bodynodexbody, bodynodeybody,50,ptAccDiff(C2:C3),'o-')














% %% Post Process
% figure(20);
% plot(timeaxis,mely,'.r','LineWidth',1)
% grid on; hold on
% set(gca,'FontSize',18)
% stokey = 0;
% for m = 1:waveComp
% stokey  = stokey+Amp(m)*cos(k(m).*xfl(21)-w(m).*timeaxis+eph(m) ) + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xfl(21)-w(m).*timeaxis+eph(m) )))  ;
% end
% plot(timeaxis,stokey,'-k','LineWidth',1)
% xprobe = xfl(idx+5);
% % [TSetaS] = NLTSelevationStokes(waveComp,Amp,k,w,eph,xprobe,timeaxis,h)
% % plot(timeaxis,TSetaS,'-k','LineWidth',0.1)
% 
% legend('MEL','Theoretical','Theoretical Interaction')
% % legend('MEL','Theoretical Interaction')
% %%
% 
% title('Comparing \eta at mid tank time-history')
% ylabel('\eta (m) ')
% xlabel('t (s)')
% %%
% figure(30);
% 
% plot(xfl,etaL,'.-r','LineWidth',1)
% grid on; hold on
% set(gca,'FontSize',18)
% hold on;
% etastokey = 0;
% for m = 1:waveComp
% etastokey  = etastokey+Amp(m)*cos(k(m).*xfl-w(m).*tiii+eph(m)) + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xfl-w(m).*tiii+eph(m)))) ;
% end
% plot(xfl,etastokey,'--k','LineWidth',0.5)
% 
% [etaS] = elevationStokes(waveComp,Amp,k,w,eph,xfl,tiii,h)
% % plot(xf,etaS,'-k','LineWidth',0.8)
% 
% legend('MEL','Theoretical','Theoretical Interaction')
% % legend('MEL','Theoretical Interaction')
% % axis square
% %%
% str = sprintf('Comparing elevation at %d seconds', tend)
% title(str)
% 
% ylabel('\eta (m)')
% xlabel('x (m)')
% 
% %%
% figure(50);
% 
% plot(xfr,etaR,'.-r','LineWidth',1)
% grid on; hold on
% set(gca,'FontSize',18)
% hold on;
% etastokey = 0;
% for m = 1:waveComp
% etastokey  = etastokey+Amp(m)*cos(k(m).*xfr-w(m).*tiii+eph(m)) + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xfr-w(m).*tiii+eph(m)))) ;
% end
% plot(xfr,etastokey,'--k','LineWidth',0.5)
% 
% [etaS] = elevationStokes(waveComp,Amp,k,w,eph,xfr,tiii,h)
% % plot(xf,etaS,'-k','LineWidth',0.8)
% 
% legend('MEL','Theoretical','Theoretical Interaction')
% % legend('MEL','Theoretical Interaction')
% % axis square
% %%
% str = sprintf('Comparing elevation at %d seconds', tend)
% title(str)
% 
% ylabel('\eta (m)')
% xlabel('x (m)')
% 
% %%
% figure(70);
% hold on
% %
% subplot(3,1,1)
% 
% plot(timeaxis,melyC1,'.-r','LineWidth',1)
% grid minor; hold on
% hold on
% ylabel('\eta (m)');
% xlabel('t(s)')
% % set(gca,'FontSize',5)
% stokey = 0;
% for m = 1:waveComp
% stokey  = stokey+Amp(m)*cos(k(m).*xfl(21)-w(m).*timeaxis+eph(m) ) + (3*k(m)*Amp(m)*Amp(m)/(4*(k(m)*h)^3))*cos(2*((k(m).*xfl(21)-w(m).*timeaxis+eph(m) )))  ;
% end
% plot(timeaxis,stokey,'-k','LineWidth',1)
% hold on;
% set(gca,'FontSize',15)
% set(gca, 'fontname','times' )
% ylabel('\eta (m)');
% xlabel('t(s)')
% axis([1 15 -0.1 0.1])
% str = sprintf(' Elevation at G1 x = %1.3g  m', xfl(21));
% title(str);
% legend('MEL','Analytical')
% %
% subplot(3,1,2)
% 
% plot(timeaxis,melyC2,'--r','LineWidth',1)
% grid minor; 
% hold on;
% axis([0 15 -0.1 0.1])
% legend('MEL')
% str = sprintf(' Elevation at G2 x = %1.3g  m', xfl(93-6));
% title(str);
% ylabel('\eta (m)');
% xlabel('t(s)')
% set(gca,'FontSize',15)
% set(gca, 'fontname','times' )
% axis([5 15 -0.1 0.1])
% %
% subplot(3,1,3)
% figure;
% plot(timeaxis,melyC3,'--r','LineWidth',1)
% grid minor; 
% hold on;
% axis([5 15 -0.03 0.03])
% legend('MEL')
% ylabel('\eta (m)');
% xlabel('t(s)')
% str = sprintf(' Elevation at G3 x = %1.3g  m ', xfr(7));
% title(str);
% load('G3etamdpi.mat')
% load('G3txmdpi.mat')
% 
% hold on
% plot(G3txmdpi,G3etamdpi,'--b')
% set(gca,'FontSize',15)
% set(gca, 'fontname','times' )
% %%
% % str = sprintf('Comparing elevation at %d seconds', tend)
% % title(str)
% 
% ylabel('\eta (m)')
% xlabel('x (m)')


%% Save Data
% sprintf('%1.3g',1/3)
% timeaxiscomp(casenoi,:)=timeaxis
% ForceXwComp(casenoi,:) = Fx
% ForceZwComp(casenoi,:) = Fz
% MomenYwComp(casenoi,:) = My
str = sprintf(' Simulation Data Full mod case no.   %1.3g %1.3g  wavelength  %s ', casenoi, lambda)
mainfolder = pwd;
mkdir (str)
cd (str)
save('WorkspaceVar')
cd(mainfolder)


close all
clc


end
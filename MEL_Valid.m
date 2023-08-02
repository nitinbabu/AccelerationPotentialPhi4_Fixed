clc
clear all
close all

tic

%% Wave Tank Parameters
dt = 0.008; t0 = 0; tend = 10;
fx = 100;
Lm = 6; h = 1;
Amp = 0.03;
linearSW = 1;  % 1: NL ; 0: Linearized
T= 0.8005;  w = 2*pi/T;
% Constants
g = 9.81; 
% [k] = fdispersion(w,g,h);
lambda = 2*pi/k;

steepness = 2*Amp/lambda
%% Initializations for Linear wave
[x,y,xf,xd,xb,xp,yf,yd,yb,yp,uin,uout,pot] = tank(fx,h,Lm,t0,w,k,Amp);
yp0 = yp;
[C1, C2, C3, C4] = corners(xf,xd,xb,xp);
nnode = length(x);  nelem  = ((length(x))/2);
corner4 = [C1 C2 C3 C4];
[qnode] = meshing(nelem);
tm = t0;

Handle1 = plot(x,y,'.-'); 
grid on; title('2D Domain, Velocity Inlet. Zero Potential on FS '); xlabel('x'); ylabel('y'); axis([-Lm/20 Lm -h h]) ;drawnow
hold on
%%
xo = xf;
eta = Amp*cos(k*xo-w*t0);
Handle2 = plot(xo,eta,'.-k');
legend('MEL2D','Benchmark')

%% Damping zone definition

n=1*lambda;
[val,idx]=min(abs(xf-n));
minVal=xf(idx); % Finds close enough values then looks for the index
xpos=Lm-xf(idx); % x0
alphad = 1; % alpha;
dzone = (fx+1)-idx;
nu = zeros(1,length(xf));
nu(dzone:end) = alphad*w*((xf(dzone:end)-xpos)/lambda).^2;
%%


counter = 0;
timeaxis = t0:dt:tend;
steps = 1:length(timeaxis);
airy = zeros(1,length(steps));
mely = airy;
stepi = 2; % initial index for final figure;
%%
while tm <tend+dt
    if tm == t0
        
        [r_t] = rampfunction(tm,T);  % Single Ramp coefficient applied to Inlet BC
        airy(1) = eta((fx/2)+1);
        %%
        % for uL yp != 0 yp is elevation itself BIG Finding 
        
        [yp] = semiLateral(yp0,t0,g,Amp,w,h,k);
        
        pot = 0.*(g*Amp/w)*(cosh(k*(1.*yf+h))/cosh(k*h)).*sin(k.*xf-w*t0);
        uin = -r_t*(g*k*Amp/w) * (cosh(k*(yp+h))/cosh(k*h)).*cos(k.*xp-w*t0) ; % Need to check the multiple fluxes value at the very last element
        uL = -r_t*(g*k*Amp/w) * (cosh(k*(yf(1)+h))/cosh(k*h)).*cos(k.*x(1)-w*t0);% -r_t*(g*k*Amp/w) * cos(k*x(1)-w*t0);
        uout = 0;

        [ktvec,knvec]=bcs(xf, xd, xb, xp,pot,uout,uin);

        
    else
        
        if counter == 10
            [yu] = FreeSurfSmooth5pt(xu,yu);
            counter = 0;
        else
            counter = counter+1;
        end
        
        xf = xu; yf = yu; pot = phiu;
        [r_t] = rampfunction(tm,T);
         
        [yp] = semiLateral(yp0,tm,g,Amp,w,h,k);
        
        uin = -r_t*(g*k*Amp/w) * (cosh(k*(yp+h))/cosh(k*h)).*cos(k.*xp-w*tm) ; % Need to check the multiple fluxes value at the very last element
        uL =  -r_t*(g*k*Amp/w) * (cosh(k*(yf(1)+h))/cosh(k*h)).*cos(k.*x(1)-w*tm);%-r_t*(g*k*Amp/w) * cos(k*x(1)-w*tm);
        uout = 0;
       
        [ktvec,knvec]=bcs(xf, xd, xb, xp,pot,uout,uin);
        
        set(Handle1,'XData',x,'YData', y); grid on ;  drawnow; %axis([-Amp Lm -h h]) ;

    eta = Amp*cos(k*xo-w*tm);
    set(Handle2,'XData',xo,'YData', eta); drawnow;
    mely(stepi) = yf((fx/2)+1);
    airy(stepi) = eta((fx/2)+1);
    
    stepi = stepi+1;
    end
    disp('time: ')
    disp(tm)
 
    %% STEP 1
    x = [xf xd xb xp];           y = [yf yd yb yp];
       
    %   stage 1
    [p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(x,y,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
    [U,V,NX,NY,detadx] = sderivatives2(xf,yf,x,y,p,dpdn,qnode); %changed
    
    %   stage 2
%     k0x = dt*U; 
    k0y = dt*(V-linearSW.*U.*detadx-nu.*yf); Vel2 = (U+V).^2; %mod
    k0p = dt*(-linearSW.*0.5*Vel2 - g*yf + linearSW.*V.*(V-U.*detadx)-nu.*p(1:fx+1)); %mode
    
    %   stage 3
%     x_i = xf + 0.5*k0x;
    y_i = yf + 0.5*k0y;
    phi_i = pot + 0.5*k0p;
    
    %% STEP 2
    %   BC need to be updated for intermediate substeps
    ti = tm+ 0.5*dt; 
    [yp] = semiLateral(yp0,ti,g,Amp,w,h,k);
        
    x = [xf xd xb xp];           y = [y_i yd yb yp];
        uin = -r_t*(g*k*Amp/w) * (cosh(k*(yp+h))/cosh(k*h)).*cos(k.*xp-w*ti) ; % Need to check the multiple fluxes value at the very last element
        uL =  -r_t*(g*k*Amp/w) * (cosh(k*(yf(1)+h))/cosh(k*h)).*cos(k.*x(1)-w*ti) ;%-r_t*(g*k*Amp/w) * cos(k*x(1)-w*ti);
        uout = 0;

        [ktvec,knvec]=bcs(xf, xd, xb, xp,phi_i,uout,uin);
    
    %   stage 1
    [p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(x,y,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
    [U,V,NX,NY,detadx] = sderivatives2(xf,y_i,x,y,p,dpdn,qnode); % changed
    
    %   stage 2

    k1y = dt*(V-linearSW.*U.*detadx-nu.*y_i); Vel2 = (U+V).^2; %mod
%     k1p = dt*(linearSW.*0.5*Vel2 - g*y_i + linearSW.*V.*(V-U.*detadx) -nu.*p(1:fx+1)); %mod
    k1p = dt*(-linearSW.*0.5*Vel2 - g*y_i + linearSW.*V.*(V-U.*detadx) -nu.*p(1:fx+1)); %mod
%             ^ mistake    
    %   stage 3
%     x_ii = xf + 0.5*k1x;
    y_ii = yf + 0.5*k1y;
    phi_ii = pot + 0.5*k1p;
    
    %% STEP 3
    %   BC need to be updated for intermediate substeps
    tii = tm+ 0.5*dt;
        [yp] = semiLateral(yp0,tii,g,Amp,w,h,k);
        
    x = [xf xd xb xp];           y = [y_ii yd yb yp];
        uin = -r_t*(g*k*Amp/w) * (cosh(k*(yp+h))/cosh(k*h)).*cos(k.*xp-w*tii) ; % Need to check the multiple fluxes value at the very last element
        uL =  -r_t*(g*k*Amp/w) * (cosh(k*(yf(1)+h))/cosh(k*h)).*cos(k.*x(1)-w*tii) ;%-r_t*(g*k*Amp/w) * cos(k*x(1)-w*tii);
        uout = 0;
        [ktvec,knvec]=bcs(xf, xd, xb, xp,phi_ii,uout,uin);
    
    %   stage 1
    [p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(x,y,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
    [U,V,NX,NY,detadx] = sderivatives2(xf,y_ii,x,y,p,dpdn,qnode); %changed
    
    %   stage 2
%     k2x = dt*U;  
    k2y = dt*(V-linearSW.*U.*detadx-nu.*y_ii); Vel2 = (U+V).^2; %mod
%     k2p = dt*(linearSW.*0.5*Vel2 - g*y_ii + linearSW.*V.*(V-U.*detadx) -nu.*p(1:fx+1));%mod
    k2p = dt*(-linearSW.*0.5*Vel2 - g*y_ii + linearSW.*V.*(V-U.*detadx) -nu.*p(1:fx+1));%mod
   %          ^ mistake      
    %   stage 3   % FOUND ONE MISTAKE did not solve the issue though
    y_iii = yf + k2y;
    phi_iii = pot + k2p;
    %% STEP 4
    %   BC need to be updated for intermediate substeps
    tiii = tm + dt; 
        [yp] = semiLateral(yp0,tiii,g,Amp,w,h,k);
        
    x = [xf xd xb xp];           y = [y_iii yd yb yp];
        uin = -r_t*(g*k*Amp/w) * (cosh(k*(yp+h))/cosh(k*h)).*cos(k.*xp-w*tiii) ; % Need to check the multiple fluxes value at the very last element
        uL =  -r_t*(g*k*Amp/w) * (cosh(k*(yf(1)+h))/cosh(k*h)).*cos(k.*x(1)-w*tiii) ;% -r_t*(g*k*Amp/w) * cos(k*x(1)-w*tiii);
          uout = 0;
        [ktvec,knvec]=bcs(xf, xd, xb, xp,phi_iii,uout,uin);
    
    %   stage 1
    [p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(x,y,nnode,nelem,qnode,ktvec,knvec,corner4,uL);
    [U,V,NX,NY,detadx] = sderivatives2(xf,y_iii,x,y,p,dpdn,qnode); %changed
    
    %   stage 2
%     k3x = dt*U;  
    k3y = dt*(V-linearSW.*U.*detadx-nu.*y_iii); Vel2 = (U+V).^2; %mod
% k3p = dt*(linearSW.*0.5*Vel2 - g*y_iii + linearSW.*V.*(V-U.*detadx) -nu.*p(1:fx+1)); %mod
k3p = dt*(-linearSW.*0.5*Vel2 - g*y_iii + linearSW.*V.*(V-U.*detadx) -nu.*p(1:fx+1)); %mod
%         ^ Mistake
    %   stage 3
    xu = xf ;%+ (k0x + 2*k1x + 2*k2x + k3x)/6;
    yu = yf + (k0y + 2*k1y + 2*k2y + k3y)/6;
    phiu = pot + (k0p + 2*k1p + 2*k2p + k3p)/6; 
    
    tm = tiii; 
    
  
end

%%
figure(2);
timeaxis = t0:dt:tend%+dt;
plot(timeaxis,airy,'--r')
hold on;plot(timeaxis,mely,'-b')
legend('Airy','MEL')
axis([0 tend+1 -(Amp+0.5*Amp) Amp+0.5*Amp])

grid on
title('MEL vs Airy at x = 3m  dt = 0.01 dzone 1x \lambda alpha = 1 smoothing = 10*dt')
xlabel('time (s)')  
ylabel('\eta (m)')


%% CFL
figure(3)
plot(xf,eta);hold on
plot(xf,yf);

hold on
grid on
xlabel('Tank Length')
ylabel('\eta (m)')

axis([0 Lm -(Amp+0.5*Amp) Amp+0.5*Amp])
legend('Airy','MEL')

Dt1 = 0.25/max((U)/2.5)

%% vonNeumann
Dt2 = sqrt( (8*[xf(2)-xf(1)])/(pi*9.81) )
toc
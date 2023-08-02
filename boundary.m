function [x,y,xf,xd,xb,xp,yf,yd,yb,yp,uin,uout,pot] = boundary(fx,h,L,t,w,k,Amp)
        xpp = 0;
        Lm = L-xpp; 
        dx = Lm/fx; 
%         lambda = 200/6;
%         dy = 1*dx;% h/4/10;
        
        dy = 1*dx;% h/4/10;
%         dy = dx;


%         xf = [0:dx:Lm];  yf = zeros(1,length(xf)); %yf(2) = 5;   yf(3) = 10;  yf(4) = 5; 
%         
% %         yf = -0.*Amp.*cos(k.*xf-w*t);
%         yd = -dy:-dy:-h+dy; xd = Lm+zeros(1,length(yd)); 
% 
%         xb = [ Lm:-dx*1:0 ];    yb = -h+zeros(1,length(xb));         %Double nodes
%         
%         yp = -h+dy:dy:-dy; xp = zeros(1,length(yp));    %Double elements
%     

        xf = [0:dx:Lm];  yf = zeros(1,length(xf)); %yf(2) = 5;   yf(3) = 10;  yf(4) = 5; 
%         yf =Amp.*cos(k.*xf-w*t);
%         yf = -0.*Amp.*cos(k.*xf-w*0);
        yd = -dy:-dy*1:-h+dy; xd = Lm+zeros(1,length(yd)); 
        xb = [ Lm:-dx*1:0 ];    yb = -h+zeros(1,length(xb));              
        yp = [-h+dy:dy:-dy     ];% -(Amp+Amp) -Amp] ;
%           yp = [ -25 -20 -15 -10 -9 -8 -7 -6 -5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 -0.25];% -0.1 0];
%         yp = [-h -20 -15 -10 -9:0.5:-2]; 

%         yp = [-0.8 -0.6 -0.4 -0.3   -0.2 -0.1 -0.09 -0.08 -0.06 -0.05 -0.045 -0.04 ];%-0.035 -0.03 -0.025 -0.02 -0.015];% -0.04 -0.03 -0.035 -0.025 -0.02 -0.017]; % Amp = 0.01
%         yp = [-0.8 -0.6 -0.4 -0.3   -0.2 -0.1 -0.08 -0.06  -0.055 -0.05 -0.045  -0.037 -0.035 ];  % Amp = 0.03
        xp = zeros(1,length(yp));    
        
x = [xf xd xb xp];           y = [yf yd yb yp];
uin = 0;uout=0;pot = 0;
 
end
function [lambda]=dispersioncalc(w,d)
%% Code from V Sundar
T=2*pi/w;
d;
g=9.81;
L0= 1.56*T^2;
L0_old=L0;
for i = 1:10000
L_new =((g*T^2)/(2*pi))*(tanh((2*pi)*d/L0_old));
err = abs(L_new-L0_old);
if err <1e-6
L=L_new;
break;
else
L0_old = L_new;
end
end
lambda = L;

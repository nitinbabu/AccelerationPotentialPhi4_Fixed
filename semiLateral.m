function  [yp] = semiLateral(yp0,tx,g,Amp,w,h,k,r_t)
%% xp =0;
%     dyp = yp0;
%     yp = yp0;
    dyp = r_t.*Amp* (sinh(k*(h+yp0))/sinh(k*h)) * (cos(k*0-w*tx));
    yp = yp0+dyp;
%     dyp(1:end-2)   = Amp* (sinh(k*(h+yp0(1:end-2)))/sinh(k*h)) * (cos(k*0-w*tx)); %% Displacement
%     yp(1:end-2) = yp0(1:end-2)+dyp(1:end-2);
%     yp(end-1) = yp0(end-1)-2*Amp;
%     yp(end) = yp0(end)-1*Amp;
end
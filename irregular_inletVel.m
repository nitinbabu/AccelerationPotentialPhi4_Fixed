function [uin, uL] = irregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etax,tx,h)
% etax is eta at some step tx
     uin = 0; uL = 0; g = 9.81;
     
     for m = 1:waveComp
        uin = uin+(  (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(yp+h))/cosh(k(m)*h)).*cos(k(m).*xp-w(m)*tx+eph(m)) );
        uL =  uL+( (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(etax(1)+h))/cosh(k(m)*h)).*cos(k(m).*0-w(m)*tx+eph(m)) );
     end
     
    uin = -r_t*uin; uL  = -r_t*uL;

end
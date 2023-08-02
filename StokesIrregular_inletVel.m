function [uin, uL] = StokesIrregular_inletVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etax,tx,h)
% etax is eta at some step tx
     uin = 0; uL = 0; g = 9.81;
     
     for m = 1:waveComp
        uin = uin+(  (g*k(m)*Amp(m)/w(m)) * (cosh(k(m).*(yp+h))/cosh(k(m)*h)).*cos(k(m).*xp-w(m)*tx+eph(m)) + (3/8)*Amp(m)*Amp(m)*w(m)*(2*k(m)).*cosh(2*k(m)*(yp+h)).*cos(2*(k(m).*xp-w(m)*tx +eph(m)))./(sinh(k(m)*h))^4  );
        uL =  uL+( (g*k(m)*Amp(m)/w(m)) * (cosh(k(m)*(etax(1)+h))/cosh(k(m)*h)).*cos(k(m).*0-w(m)*tx+eph(m)) + (3/8)*Amp(m)*Amp(m)*w(m)*(2*k(m)).*cosh(2*k(m)*(etax(1)+h)).*cos(2*(k(m).*0-w(m)*tx +eph(m)))./(sinh(k(m)*h))^4  );
     end
     
    uin = -r_t*uin; uL  = -r_t*uL;

end
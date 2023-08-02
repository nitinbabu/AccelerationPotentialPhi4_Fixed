function [ain, aL] = StokesIrregular_inletACC(r_t,waveComp,Amp,k,w,eph,xp,yp,etax,tx,h)
% etax is eta at some step tx
     ain = 0; aL = 0; g = 9.81;
     
     for m = 1:waveComp
        ain = ain+(  (g*k(m)*Amp(m)) * (cosh(k(m).*(yp+h))/cosh(k(m)*h)).*sin(k(m).*xp-w(m)*tx+eph(m)) + 2*w(m)*(3/8)*Amp(m)*Amp(m)*w(m)*(2*k(m)).*cosh(2*k(m)*(yp+h)).*sin(2*(k(m).*xp-w(m)*tx +eph(m)))./(sinh(k(m)*h))^4  );
        aL =  aL+( (g*k(m)*Amp(m)) * (cosh(k(m)*(etax(1)+h))/cosh(k(m)*h)).*sin(k(m).*0-w(m)*tx+eph(m)) + 2*w(m)*(3/8)*Amp(m)*Amp(m)*w(m)*(2*k(m)).*cosh(2*k(m)*(etax(1)+h)).*sin(2*(k(m).*0-w(m)*tx +eph(m)))./(sinh(k(m)*h))^4  );
       % WAS Ul INSTEAD OF aL
     end
     
    ain = -r_t*ain; aL  = -r_t*aL;

end
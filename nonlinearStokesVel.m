%analytical ETA
function [uin,uL] = nonlinearStokesVel(r_t,waveComp,Amp,k,w,eph,xp,yp,etax,t,d)
 u1 = 0; u2a = 0; u2b = 0; u2c = 0; g = 9.81;
 u1L = 0; u2aL = 0; u2bL = 0; u2cL = 0;
 
 x0 = 0; t0 = 0;

 for i = 1:waveComp
    u1 = u1 + ( g*Amp(i)*k(i)/w(i) ) * cosh(k(i).*(yp+d))/cosh(k(i)*d)   .* cos(k(i).*(xp-x0) - w(i)*(t-t0) + eph(i));
    u1L = u1L + ( g*Amp(i)*k(i)/w(i) ) * cosh(k(i)*(etax(1)+d))/cosh(k(i)*d)   * cos(k(i)*(0-x0) - w(i)*(t-t0) + eph(i));
     for j =1:waveComp % SUM FREQUENCIES
          if j>i
%                   HP(w(i),w(j),k(i),k(j),d)
%                   GP(w(i),w(j),k(i),k(j),d)
              u2a = u2a +   Amp(i)*Amp(j)*(k(i)+k(j))*(  ( GP(w(i),w(j),k(i),k(j),d)/DP(w(i),w(j),k(i),k(j),d) ) .* ( cosh((k(i)+k(j)).*(yp+d)) /  cosh((k(i)+k(j))*d)  )    ).*cos( (k(i)+k(j)).*(xp-x0)-(w(i)+w(j))*(t-t0) + (eph(i)+eph(j))) ;
              u2aL = u2aL +   Amp(i)*Amp(j)*(k(i)+k(j))*(  ( GP(w(i),w(j),k(i),k(j),d)/DP(w(i),w(j),k(i),k(j),d) ) * ( cosh((k(i)+k(j))*(etax(1)+d)) /  cosh((k(i)+k(j))*d)  )    )*cos( (k(i)+k(j))*(0-x0)-(w(i)+w(j))*(t-t0) + (eph(i)+eph(j))) ;
        
          end           
     end

      for j =1:waveComp % DIFFERENCE FREQUENCIES
          if j>i
                u2b = u2b +   Amp(i)*Amp(j)*(k(i)-k(j))*(  ( GN(w(i),w(j),k(i),k(j),d)/DN(w(i),w(j),k(i),k(j),d) ) * ( cosh((k(i)-k(j)).*(yp+d)) /  cosh((k(i)-k(j))*(d))  )    ).*cos( (k(i)-k(j)).*(xp-x0)-(w(i)-w(j))*(t-t0) + (eph(i)-eph(j))) ;
                u2bL = u2bL +   Amp(i)*Amp(j)*(k(i)-k(j))*(  ( GN(w(i),w(j),k(i),k(j),d)/DN(w(i),w(j),k(i),k(j),d) ) * ( cosh((k(i)-k(j))*(etax(1)+d)) /  cosh((k(i)-k(j))*(d))  )    )*cos( (k(i)-k(j))*(0-x0)-(w(i)-w(j))*(t-t0) + (eph(i)-eph(j))) ;
          
          end           
      end        


    u2c = u2c +   Amp(i)^2*k(i)* ( GP(w(i),w(j),k(i),k(j),d)/DP(w(i),w(j),k(i),k(j),d) ) *( cosh((2*k(i)*(yp+d)))/cosh(2*k(i)*d)  )  .*sin( 2*(k(i).*(xp-x0)-w(i)*(t-t0) + eph(i)) ) ;
    u2cL = u2cL +   Amp(i)^2*k(i)* ( GP(w(i),w(j),k(i),k(j),d)/DP(w(i),w(j),k(i),k(j),d) ) *( cosh((2*k(i)*(etax(1)+d)))/cosh(2*k(i)*d)  )  *sin( 2*(k(i)*(0-x0)-w(i)*(t-t0) + eph(i)) ) ;

 end

     
uS = u1+u2a+u2b+u2c;
uSL = u1L+u2aL+u2bL+u2cL;
% This was tested against a simple stokes inlet velocity model

uin = -r_t.*uS;
uL = -r_t*uSL;
     
     
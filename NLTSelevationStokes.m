%analytical ETA
function [TSetaS] = NLTSelevationStokes(waveComp,Amp,k,w,eph,xprobe,timeseries,d)
 eta1 = 0; eta2a = 0; eta2b = 0; eta2c = 0; g = 9.81;
 x0 = 0; t = timeseries;t0 = 0;
 
 for i = 1:waveComp
    eta1 = eta1 + Amp(i)*cos(k(i).*(xprobe-x0) - w(i)*(t-t0) + eph(i));
     for j =1:waveComp
          if j>i
%                   HP(w(i),w(j),k(i),k(j),d)
%                   GP(w(i),w(j),k(i),k(j),d)
              eta2a = eta2a +   Amp(i)*Amp(j)*HP(w(i),w(j),k(i),k(j),d)*cos( (k(i)+k(j)).*(xprobe-x0)-(w(i)+w(j)).*(t-t0) + (eph(i)+eph(j))) ;
          end           
     end

      for j =1:waveComp
          if j>i
              eta2b = eta2b +   Amp(i)*Amp(j)*HN(w(i),w(j),k(i),k(j),d)*cos( (k(i)-k(j)).*(xprobe-x0)-(w(i)-w(j)).*(t-t0) + (eph(i)-eph(j))) ;
          end           
      end        


    eta2c = eta2c +   Amp(i)^2*HP(w(i),w(j),k(i),k(j),d)*sin( 2*(k(i).*(xprobe-x0)-w(i).*(t-t0) + eph(i)) ) ;

 end

     
     TSetaS = eta1+eta2a+eta2b+eta2c;
     
     
end
     
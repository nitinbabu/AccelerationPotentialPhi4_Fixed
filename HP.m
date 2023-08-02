function [Hp]=HP(wi,wj,ki,kj,d)
 
   Hp = ((wi+wj)/9.81) * (GP(wi,wj,ki,kj,d)/DP(wi,wj,ki,kj,d)) + 0.5*(ki*tanh(ki*d) + kj*tanh(kj*d)) -...
...         +              +                 +
           (9.81/2)*(ki*kj/wi/wj) * (cosh((ki-kj)*d)/(cosh(ki*d)*cosh(kj*d)));
%                                            -
    
end
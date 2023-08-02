function [Hn]=HN(wi,wj,ki,kj,d)
 
   Hn = ((wi-wj)/9.81) * (GN(wi,wj,ki,kj,d)/DN(wi,wj,ki,kj,d)) + 0.5*(ki*tanh(ki*d) + kj*tanh(kj*d)) -...
...         -              -                 -
           (9.81/2)*(ki*kj/wi/wj) * (cosh((ki+kj)*d)/(cosh(ki*d)*cosh(kj*d)));
%                                            +
    
end
function [Dn]=DN(wi,wj,ki,kj,d)
 
   Dn = -(wi-wj)^2 + 9.81*(ki-kj)*tanh((ki-kj)*d);
%           -                -            -
    
end
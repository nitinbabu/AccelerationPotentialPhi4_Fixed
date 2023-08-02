function [Dp]=DP(wi,wj,ki,kj,d)
 
   Dp = -(wi+wj)^2 + 9.81*(ki+kj)*tanh((ki+kj)*d);
%           +                +            +
    
end
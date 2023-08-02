function [Gn]=GN(wi,wj,ki,kj,d)
 
    Gn = -9.81^2*( (ki^2/(2*wi*cosh(ki*d)^2))  - (kj^2/(2*wj*cosh(kj*d)^2))   +...
                ...                            -
                    ( (ki*kj)/(wi*wj) * (wi-wj)*(1+tanh(ki*d)*tanh(kj*d))  )  ) ;
%                                          -      +
    
end
function [N,Nd] = shapefunction(zeta)

N = [(0.5*zeta*(zeta-1)), 1-zeta^2 , (0.5*zeta*(zeta+1))];
    
Nd = [ zeta-0.5, -2*zeta , zeta+0.5 ];

% end
function [J ,nx ,ny ] = jacobian(xc,yc,Nd) %dependent on zeta starting shapefunction.m
    
    dxdeta = Nd(1).*xc(1)  +Nd(2).*xc(2)  +Nd(3).*xc(3) ;
    dydeta = Nd(1).*yc(1) + Nd(2).*yc(2) + Nd(3).*yc(3);
    
    J = sqrt( (dydeta)^2 + (dxdeta)^2 );
    nx=-dydeta/J;
    ny=dxdeta/J;
    
    
% end
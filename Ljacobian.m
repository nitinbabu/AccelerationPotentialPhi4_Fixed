function [J ,nx ,ny ] = Ljacobian(xc,yc,NLd) %dependent on eta starting shapefunction.m
    %% Do I need normal vector or unit normal vector? CHECK
%     ny = -NLd(1).*xc(1) + -NLd(2).*xc(2) + -NLd(3).*xc(3) ;
%     nx = NLd(1).*yc(1) + NLd(2).*yc(2) + NLd(3).*yc(3);
%     
%     J = sqrt( (nx)^2 + (ny)^2 );
    
    dxdeta = NLd(1).*xc(1)  +NLd(2).*xc(2)  +NLd(3).*xc(3) ; %03.07.21
    dydeta = NLd(1).*yc(1) + NLd(2).*yc(2) + NLd(3).*yc(3);
    
    J = sqrt( (dydeta)^2 + (dxdeta)^2 );
    nx=-dydeta/J;
    ny=dxdeta/J;
    
end
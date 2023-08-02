function [intA, intB] = lineintK(P,~,x,y,c,xc,yc)
[wg, xig,~,~] = GuassianQuad8points();
xp = x(P); yp = y(P);
% xq = x(Q); yq = y(Q);


tempA = 0; tempB = 0;

for ig = 1:length(xig)
    zeta = xig(ig);
    [N,Nd] = shapefunction(zeta);

    [J ,nx ,ny ] = jacobian(xc,yc,Nd);
%     nx = -nx; ny = -ny; %Grilli 1989 Eq 31
    
    xq = N(1)*xc(1) + N(2)*xc(2) + N(3)*xc(3);
    yq = N(1)*yc(1) + N(2)*yc(2) + N(3)*yc(3);
    
    
    rx = xq-xp ;
    ry = yq-yp;
 
    
    r = sqrt((xp-xq)^2 + (yp-yq)^2) ;
    K1 = -(rx*nx+ry*ny)/r^2/(2*pi) ;
    
    
    K2 = log(1/r)*1/(2*pi);
    
    tempA = tempA + K1*N(c)*J*wg(ig);
    tempB = tempB + K2*N(c)*J*wg(ig);
    

end

intA = tempA;
intB = tempB;

% end
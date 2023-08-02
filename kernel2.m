function intB = kernel2(~,~,~,~,c,xc,yc,in)
%Term with Normal Gaussan Quadrature

[wg, xig,~,~] = GuassianQuad8points();
% xp = x(P); yp = y(P);
% xq = x(Q); yq = y(Q);
xxq = 0; yyq = 0;
tempB = 0;

for ig = 1:length(xig)
    zeta = xig(ig); %!
    [N,Nd] = shapefunction(zeta);

    [J ,~ ,~ ] = jacobian(xc,yc,Nd);

    if in ==1
        xxq = (zeta-2)*xc(1)+2*(1-zeta)*xc(2)+zeta*xc(3);
        yyq = (zeta-2)*yc(1)+2*(1-zeta)*yc(2)+zeta*yc(3);
    elseif in == 2
        xxq = -0.5*(zeta-1)*xc(1)+zeta*xc(2)-0.5*(zeta+1)*xc(3);
        yyq = -0.5*(zeta-1)*yc(1)+zeta*yc(2)-0.5*(zeta+1)*yc(3);
    elseif in == 3
        xxq = -zeta*xc(1)+2*(1+zeta)*xc(2)-(zeta+2)*xc(3);
        yyq = -zeta*yc(1)+2*(1+zeta)*yc(2)-(zeta+2)*yc(3);
    end
    
    r = sqrt(xxq^2+yyq^2);

    K2 = log(1/r)*1/(2*pi);
    

    tempB = tempB + K2*N(c)*J*wg(ig);
    
    
    

end


intB = tempB; %1/2*pi in K2

% end
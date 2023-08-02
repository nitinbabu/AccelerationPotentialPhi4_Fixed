function intB = singular(P,Q,x,y,c,xc,yc,in)
%Term with Logarithmic Gaussan Quadrature

[wg, xig,wgl8,etag8] = GuassianQuad8points();
xp = x(P); yp = y(P);
% xq = x(Q); yq = y(Q);

tempB = 0;

for ig = 1:length(etag8)
    eta = etag8(ig); %!

    if in ==1
        ic =1; %to Check how ic fits in J and xc yc
        [NL,NLd] = Lshapefunction(eta,ic);
        [J ,nx ,ny ] = Ljacobian(xc,yc,NLd);%needs ic   
        dzde = 2;
        tempB = tempB + wgl8(ig)*NL(c)*J*dzde;
    elseif in == 2
        ic =2;
        [NL,NLd] = Lshapefunction(eta,ic);
        [J ,nx ,ny ] = Ljacobian(xc,yc,NLd);%needs ic   
        dzde = 1;
        tempB = tempB + wgl8(ig)*NL(c)*J*dzde;
        ic =3;
        [NL,NLd] = Lshapefunction(eta,ic);
        [J ,nx ,ny ] = Ljacobian(xc,yc,NLd);%needs ic 
        dzde =1;
        tempB = tempB + wgl8(ig)*NL(c)*J*dzde;
    elseif in == 3
        ic = 4;
        [NL,NLd] = Lshapefunction(eta,ic);
        [J ,nx ,ny ] = Ljacobian(xc,yc,NLd);%needs ic   
        dzde = 2;
        tempB = tempB + wgl8(ig)*NL(c)*J*dzde;

    end
    
end

intB = tempB/(2*pi);

% end
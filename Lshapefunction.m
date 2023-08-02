function[NL,NLd] = Lshapefunction(eta,ic)

if ic == 1
    zeta = 2*eta-1;  
elseif ic == 2
    zeta = eta;
elseif ic == 3
    zeta = -eta;
elseif ic == 4
    zeta = 1-2*eta;
end

NL = [(0.5*zeta*(zeta-1)), 1-zeta^2 , (0.5*zeta*(zeta+1))];
    
NLd = [ zeta-0.5, -2*zeta , zeta+0.5 ];

end
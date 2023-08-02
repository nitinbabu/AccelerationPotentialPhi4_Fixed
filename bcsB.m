function [ktvec,knvec]=bcsB( xfl,xbody,xfr, xd, xb, xp,phi_jL,Vn,phi_jR,uout,uin)

%%       Free Surface Upstream          BODY                        Free Surface Downstream     Right Boundary          Bottom Boundary         Inlet/Paddle
ktvec = [1+zeros(1,length(xfl))         2+zeros(1,length(xbody))    1+zeros(1,length(xfr))      2+zeros(1,length(xd)) 2+zeros(1,length(xb))  2+zeros(1,length(xp))];
knvec = [phi_jL+zeros(1,length(xfl))    Vn+zeros(1,length(xbody))   phi_jR+zeros(1,length(xfl)) uout+zeros(1,length(xd)) 0+zeros(1,length(xb))  uin+zeros(1,length(xp))];
% 

 
end
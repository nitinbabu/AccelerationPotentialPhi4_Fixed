function [ktvec,knvec]=bcs( xf, xd, xb, xp,phi_j,uout,uin)

%%       Free Surface           Right Boundary          Bottom Boundary         Inlet/Paddle
ktvec = [1+zeros(1,length(xf))  2+zeros(1,length(xd)) 2+zeros(1,length(xb))  2+zeros(1,length(xp))];
knvec = [phi_j+zeros(1,length(xf))  uout+zeros(1,length(xd)) 0+zeros(1,length(xb))  uin+zeros(1,length(xp))];
% 

 
end
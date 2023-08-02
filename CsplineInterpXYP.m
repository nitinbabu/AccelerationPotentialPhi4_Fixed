function [Bx,By, Bp] = CsplineInterpXYZ(nodes,t,x,y,pot)
Px = x;
Py = y;
Ppot = pot;

%% Equation 5-15 % Rvector  # RHS vector [M][P'] = {R}
Rvecx = zeros(1,nodes); Rvecy = Rvecx;
Rvecx(1)  = 3/(2*t(2))*(x(2)-x(1));
Rvecy(1)  = 3/(2*t(2))*(y(2)-y(1));
Rvecp(1)  = 3/(2*t(2))*(pot(2)-pot(1));

for k = 2:nodes-1 
   Rvecx(k)= 3*( t(k)^2*(x(k+1)-x(k)) + t(k+1)^2*(x(k)-x(k-1)) )/(t(k)*t(k+1));
   Rvecy(k)= 3*( t(k)^2*(y(k+1)-y(k)) + t(k+1)^2*(y(k)-y(k-1)) )/(t(k)*t(k+1));
   Rvecp(k)= 3*( t(k)^2*(pot(k+1)-pot(k)) + t(k+1)^2*(pot(k)-pot(k-1)) )/(t(k)*t(k+1));
   
end
Rvecx(nodes)  = (6/t(nodes)) * (x(nodes)-x(nodes-1)) ;
Rvecy(nodes)  = (6/t(nodes)) * (y(nodes)-y(nodes-1)) ;
Rvecp(nodes)  = (6/t(nodes)) * (pot(nodes)-pot(nodes-1)) ;
%% Equation 5-15 % MatM

MatM = zeros(nodes,nodes);
for kd = 1:nodes
    for ka = 1:nodes

        if kd == 1 || kd == nodes  % If first or last node
%             MatM(kd,kd) = 1; %% I am not sure about this
            
        elseif kd>1 && kd<nodes
            MatM(kd,kd-1)    =   t(kd+1);
            MatM(kd,kd)    =   2*(t(kd)+t(kd+1));
            MatM(kd,kd+1)    =   t(kd);
        end
%         MatM(kd,ka) = ;

    end
end
MatM(1,1) = 1;      MatM(1,2) = 0.5; 
MatM(nodes,nodes) = 4;  MatM(nodes,nodes-1) = 2;


%% [P'] =inv[M]{R} for x and y each
pprimex = (MatM)\Rvecx';
pprimey = (MatM)\Rvecy';
pprimep = (MatM)\Rvecp';


PPx  = pprimex; PPy = pprimey; PPpot = pprimep;


for k = 1:nodes-1 % Each segment
    B1x(k)= Px(k);
    B2x(k) = PPx(k);
    B3x(k) = 3*(Px(k+1)-Px(k))/t(k+1)^2   - 2*PPx(k)/t(k+1) - PPx(k+1)/t(k+1);
    B4x(k) = 2*(Px(k)-Px(k+1))/t(k+1)^3   + PPx(k)/t(k+1)^2 + PPx(k+1)/t(k+1)^2;
    Bx(1:4,k) = [B1x(k) B2x(k) B3x(k) B4x(k)];

    B1y(k)= Py(k);
    B2y(k) = PPy(k);
    B3y(k) = 3*(Py(k+1)-Py(k))/t(k+1)^2   - 2*PPy(k)/t(k+1) - PPy(k+1)/t(k+1);
    B4y(k) = 2*(Py(k)-Py(k+1))/t(k+1)^3   + PPy(k)/t(k+1)^2 + PPy(k+1)/t(k+1)^2;
    By(1:4,k) = [B1y(k) B2y(k) B3y(k) B4y(k)];
    
    B1p(k)= Ppot(k);
    B2p(k) = PPpot(k);
    B3p(k) = 3*(Ppot(k+1)-Ppot(k))/t(k+1)^2   - 2*PPpot(k)/t(k+1) - PPpot(k+1)/t(k+1);
    B4p(k) = 2*(Ppot(k)-Ppot(k+1))/t(k+1)^3   + PPpot(k)/t(k+1)^2 + PPpot(k+1)/t(k+1)^2;
    Bp(1:4,k) = [B1p(k) B2p(k) B3p(k) B4p(k)];


end


end
function [p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEMBODY(x,y,nnode,nelem,qnode,ktvec,knvec,corner6,uL)
p = zeros(1,nnode); dpdn = zeros(1,nnode);
% [SysA,b,MatA,MatB]=SysEquations(x,y,nnode,nelem,qnode,knvec,ktvec); %31.01.22

Jmid = 0; Ivec = 0;
% [SysA,b,MatA,MatB,Ivec,Jmid]=SysEquationsMFBODY(x,y,nnode,nelem,qnode,knvec,ktvec,corner6,uL); %08.09.22

[SysA,b,MatA,MatB,Ivec,Jmid]=SysEquationsMFBODY_mex(x,y,nnode,nelem,qnode,knvec,ktvec,corner6,uL); %27.05.2023
zvec = SysA\b; %BeckerFormul.
MatA; MatB;
% *****Recollecting all BCs from z ****
for Q = 1:nnode     
    if ktvec(Q) == 1
       p(Q) = knvec(Q);
       dpdn(Q) = zvec(Q);
    else
        dpdn(Q)= knvec(Q);
        p(Q) = zvec(Q);
    end
end

end
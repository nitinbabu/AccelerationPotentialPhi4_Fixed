function [p,dpdn,MatA,MatB,Ivec,Jmid]=HOBEM(x,y,nnode,nelem,qnode,ktvec,knvec,corner4,uL)
p = zeros(1,nnode); dpdn = zeros(1,nnode);
% [SysA,b,MatA,MatB]=SysEquations(x,y,nnode,nelem,qnode,knvec,ktvec); %31.01.22

Jmid = 0; Ivec = 0;
[SysA,b,MatA,MatB,Ivec,Jmid]=SysEquationsMF_mex(x,y,nnode,nelem,qnode,knvec,ktvec,corner4,uL); %08.09.22
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
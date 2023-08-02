function [A,b,testA,testB]=SysEquations(x,y,nnode,nelem,qnode,knvec,ktvec)
A = zeros(nnode,nnode);
b = zeros(nnode,1);
M = nelem;
MatA = zeros(nnode,nnode); MatB = zeros(nnode,nnode);
for P = 1:nnode
    for m = 1:M
    
        for c = 1:3
            Q = qnode(m,c);
            if Q == P
%                 MatA(P,Q) = 0;
                xc = x(qnode(m,:));
                yc = y(qnode(m,:));
                in = c;
                if in == 1
                MatB(P,Q) = MatB(P,Q)+singular(P,Q,x,y,c,xc,yc,in)+kernel2(P,Q,x,y,c,xc,yc,in);                                
                elseif in == 2
                MatB(P,Q) = MatB(P,Q)+singular(P,Q,x,y,c,xc,yc,in)+kernel2(P,Q,x,y,c,xc,yc,in);           
                elseif in == 3
                MatB(P,Q) = MatB(P,Q)+singular(P,Q,x,y,c,xc,yc,in)+kernel2(P,Q,x,y,c,xc,yc,in);
                end
                
            else
                xc = x(qnode(m,:));
                yc = y(qnode(m,:));
                [intA, intB] = lineintK(P,Q,x,y,c,xc,yc);
%          
                MatA(P,Q) = MatA(P,Q)+intA;
                MatB(P,Q) = MatB(P,Q)+intB;
                
            end
   
        end
   
    end
end

diagA = 0; %!! Checked and validated
for P = 1:nnode %Source Point
    for Q = 1:nnode %Field Point
        if P ~= Q
        diagA = diagA-MatA(P,Q);
        else
        end
    end
    MatA(P,P) = diagA;
    diagA = 0;
end

MatA;
MatB; %Fine (maybe need further scrutiny)

% Applying the BCs following Becker formulation, a.k.a basic arithmetic

A = zeros(size(MatA));

% * * * * * * * * * * * * * * * * * * * * *

for P = 1:nnode
    sumB = 0; stor =1;
    for Q = 1:nnode
        
        if ktvec(Q) == 1 %if potential type B.C.
            sumB = sumB -MatA(P,Q)*knvec(Q);%MatA(P,Q)*knvector(m,c):
%             disp(sumB)
            A(P,Q) = -MatB(P,Q);
        elseif ktvec(Q) == 3
           
            
            sumB = sumB + MatB(P,Q)*beta1;%(beta); %(which is beta1) ;
%             disp(sumB)
            A(P,Q) = MatA(P,Q) -alpha1*MatB(P,Q);%- alpha*MatB(P,Q);
        else
            sumB = sumB + MatB(P,Q)*knvec(Q);%MatB(P,Q)*knvector(m,c);
%             disp(sumB)
            A(P,Q) = MatA(P,Q);
        end

    end
     b(P) = sumB;
end

testA = MatA;
testB = MatB;
end


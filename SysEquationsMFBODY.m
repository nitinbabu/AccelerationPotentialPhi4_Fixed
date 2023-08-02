function [SysA,b,MatA,MatB,Ivec,Jmid]=SysEquationsMFBODY(x,y,nnode,nelem,qnode,knvec,ktvec,corner6,uL)
%%
% corner4 = [C1 C2 C3 C4];
%% Made faster by filling MatA in the same loop. Simple fix 08/01/23
%% Almost as fast as default HOBEM
C1 = corner6(1); C2 = corner6(2); C3 = corner6(3); C4 = corner6(4); C5 = corner6(5); C6 = corner6(6);

SysA = zeros(nnode,nnode);
b = zeros(nnode,1);
M = nelem;
MatA = zeros(nnode,nnode); MatB = zeros(nnode,3*nelem); testB = MatB;
% for P = 1:nnode
%     for m = 1:M
%         for c = 1:3
%             Q = qnode(m,c);
%             if Q ~= P
%                 xc = x(qnode(m,:));                yc = y(qnode(m,:));
%                 [intA, ~] = lineintK(P,Q,x,y,c,xc,yc);
%                 MatA(P,Q) = MatA(P,Q)+intA;                
%             end
%         end   
%     end
% end

% diagA = 0; %!! Checked and validated
% for P = 1:nnode %Source Point
%     for Q = 1:nnode %Field Point
%         if P ~= Q
%         diagA = diagA-MatA(P,Q);
%         else
%         end
%     end
%     MatA(P,P) = diagA;
%     diagA = 0;
% end

%% Took 3 Hrs on 07/09/22. How long will it take for 3D?

Ivec = zeros(1,3*nelem);

Mvec = 1:nelem; Qmid = zeros(size(Mvec)); Jmid = Qmid;%Jmid = zeros(length(nelem));
count = 1;
for j = 1:nnode
    if rem(j,2)==0  ;Qmid(count) = j   ; count = count+1;
    end
end
for j = 1:nelem
    Jmid(j) = Qmid(j) + Mvec(j)-1;
end


for P = 1:nnode
   e = 1; c= 0;
    for QM = 1:3*nelem


        if QM == 3*nelem

            Q = 1;  c = 3;  xc = x(qnode(e,:));     yc = y(qnode(e,:));  
            if Q == P             
                in = c;
                if in == 1
                testB(P,QM) =singular(P,Q,x,y,c,xc,yc,in)+kernel2(P,Q,x,y,c,xc,yc,in);                                
                elseif in == 2
                testB(P,QM) = singular(P,Q,x,y,c,xc,yc,in)+kernel2(P,Q,x,y,c,xc,yc,in);           
                elseif in == 3
                testB(P,QM) = singular(P,Q,x,y,c,xc,yc,in)+kernel2(P,Q,x,y,c,xc,yc,in);
                end
                
            else
            
            [intA, intB] = lineintK(P,Q,x,y,c,xc,yc);
            testB(P,QM) = intB;
            MatA(P,Q) = MatA(P,Q)+intA;  
%             
            end

%             testB(P,QM) = 10*QM;    testC(P,QM) = c;
            Ivec(QM) = 1;


        else
            Q = QM-(e-1);   % Logical relation between MFlux nodes and local element nodes, they are an arithematic progression
            c= c+1;
            xc = x(qnode(e,:));                yc = y(qnode(e,:)); 
            if Q == P

                in = c;
                if in == 1
                testB(P,QM) = singular(P,Q,x,y,c,xc,yc,in)+kernel2(P,Q,x,y,c,xc,yc,in);                                
                elseif in == 2
                testB(P,QM) = singular(P,Q,x,y,c,xc,yc,in)+kernel2(P,Q,x,y,c,xc,yc,in);           
                elseif in == 3
                testB(P,QM) = singular(P,Q,x,y,c,xc,yc,in)+kernel2(P,Q,x,y,c,xc,yc,in);
                end  
            else
            [intA, intB] = lineintK(P,Q,x,y,c,xc,yc);
            testB(P,QM) = intB;
            MatA(P,Q) = MatA(P,Q)+intA;  
%             testB(P,QM) =  10*QM; testC(P,QM) = c;
            Ivec(QM) = Q;
            end
        end
        if rem(QM,3) == 0 % Reached element last node? Repeat node
            e =  e+1;
            c=0;
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


MatB = testB;
%%
% Applying the BCs following Becker formulation, a.k.a basic arithmetic
Jvec = 1:3*nelem;

sktvec = zeros(1,3*nelem); sknvec = zeros(1,3*nelem); 
j = 1;
while j<  3*nelem
    icount = Ivec(j); % introduce a condition check for corner reps
    if icount == C2 || icount == C4 || icount == C6
    
                sktvec(j) = ktvec(icount);
                sknvec(j) = knvec(icount);

                sktvec(j+1) = ktvec(icount+1);
                sknvec(j+1) = knvec(icount+1);
                j = j+2;

    elseif icount == C3 || icount == C5 
                sktvec(j) = ktvec(icount-1);
                sknvec(j) = knvec(icount-1);    
                
                j = j+1;

                sktvec(j) = ktvec(icount);
                sknvec(j) = knvec(icount);

                j = j +1;
    
    else
                sktvec(j) = ktvec(icount);
                sknvec(j) = knvec(icount);

                j = j+1;
    end
end

sktvec(3*nelem) = 2;
sknvec(3*nelem) = uL;

swapA = zeros(size(MatA));          swapB = zeros(size(MatB)); 
swapB = MatB;                       swapA = MatA;



[swapA,swapB] = swapperAB1(MatA,MatB,sktvec,Ivec,nnode,nelem,Jmid);



SysA = swapA;
b = swapB*sknvec';

end











            


       

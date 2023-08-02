function [swapA,swapB] = swapperAB1(MatA,MatB,sktvec,Ivec,nnode,nelem,Jmid)
%% HAGUE+SWAN SWAPPING ALGO

%% Seems to be working on paper now at 1551Hrs 14.09.2022 not in the laplace code
%% There is some error in the swapping algorithm in the Neumann Neumann case and specifically Dirichlet to neumann case

%% 1050 15.09.20222 There seems to be a real bug in N-N Case

%% 1332 15.09.2022  There was an error in the Main Ktvec Knvec aka becker BC vectors
swapB = MatB;
swapA = MatA;

%MatB(1,1) = 10;   % Uncomment only for testing test_swap.m
for P = 1:nnode
    for J = 2:3*nelem
        I = Ivec(J);

        [index] = binarySearch(Jmid, nelem, J); [index2] = binarySearch(Jmid, nelem, J+1);

        if index ~=0 && J == Jmid(index) % Is J one of the mid nodes
           
            if sktvec(J) == 2       % Known gradient (Keep on Right)
                    swapB(P,J) =  MatB(P,J);
                    swapA(P,I) =  MatA(P,I);
            elseif sktvec(J) == 1   % Known potential (Go to right with sign change)    
                    swapB(P,J) =  -MatA(P,I);
                    swapA(P,I) =  -MatB(P,J);
            end



        else     % element End nodes
                if J ~= 3*nelem %All but last
                    if sktvec(J) == 2 && sktvec(J+1) == 2
                       % No change Neumann Neighbours                   
                    
                    elseif sktvec(J) == 2 && sktvec(J+1) ==1
                        % Neumann Dirichlet Corner
                        swapA(P,I) = -MatB(P,J+1); swapB(P,J+1) = -MatA(P,I);
                    elseif sktvec(J) == 1 && sktvec(J+1) ==1 &&  index2 ==0  
                        % Dirichlet Neighbours
                        swapA(P,I) = -MatB(P,J)-MatB(P,J+1);   swapB(P,J) = -MatA(P,I); swapB(P,J+1) = 0;
                    elseif sktvec(J) == 1 && sktvec(J+1) == 2
                        % Dirichlet Neumann Corner
                        swapA(P,I) = -MatB(P,J); swapB(P,J) = -MatA(P,I);                  
                    end

                else    % LAST NODE
                    
                    if sktvec(J) == 2 && sktvec(1) ==   1
                        % Neumann Dirichlet Corner
                        swapA(P,I) = -MatB(P,1); swapB(P,1) = -MatA(P,I);
                    elseif sktvec(J) == 1 && sktvec(1) == 2
                        % Dirichlet Neumann Corner  (For 2D NWT Might be unnesessary for our convention)
                        swapA(P,I) = -MatB(P,J); swapB(P,J) = -MatA(P,I);
                        
                    end

                end
         
                



        end
    end
end




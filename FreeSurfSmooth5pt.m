%% Input a dirty free surface, clean it out
function [yf] = FreeSurfSmooth5pt(x,y)
N   = length(y);
f = y;
F = zeros(1,length(y));


for j  = 3:N-2

    dxjm2 = x(j-2) - x(j);
    dxjm1 = x(j-1) - x(j);
    dxjp1 = x(j+1) - x(j);
    dxjp2 = x(j+2) - x(j);
    
    a0 = [  (-dxjm1*dxjp1)/(2*dxjm2*(dxjp2-dxjm2))   ]*f(j-2) + ...
...
            [   (-dxjp1)/(2*(dxjm1-dxjp1))  ]*f(j-1) + ...
...
                [ (dxjm1*dxjp1)/(2*dxjm2*dxjp2)   + 0.5 ]*f(j) + ...
...
                    [  (dxjm1)/(2*(dxjm1-dxjp1))   ]*f(j+1) + ...
...
                        [  (dxjm1*dxjp1)/(2*dxjp2*(dxjp2-dxjm2))   ]*f(j+2) ;
            





    F(j) = a0;

    %% Not sure how to use these four points  \\ 17.11.22 Just Plug it using n-2 and 3 number nodes it was simple
    if j ==3
        
            a1 = [  (dxjp2 + dxjp1 + dxjm1)/(2*dxjm2*(dxjp2-dxjm2))   ]*f(j-2) + ...
...
            [   (-1)/(2*(dxjp1-dxjm1))  ]*f(j-1) + ...
...
                [ -(dxjm2+dxjm1+dxjp1+dxjp2)/(2*dxjm2*dxjp2) ]*f(j) + ...
...
                    [  (1)/(2*(dxjp1-dxjm1))   ]*f(j+1) + ...
...
                        [  -(dxjp1+dxjm1+dxjm2)/(2*dxjp2*(dxjp2-dxjm2))   ]*f(j+2) ;



    a2 = [-1/(dxjm2*(dxjp2-dxjm2))]*f(j-2) + [1/(dxjm2*dxjp2)]*f(j) + [1/(dxjp2*(dxjp2-dxjm2))]*f(j+2);

        F1 = a0     +    a1*(x(1)-x(3)) +    a2*(x(1)-x(3))^2;
        F2 = a0     +    a1*(x(2)-x(3)) +    a2*(x(2)-x(3))^2;
    end

% F2 = y(2)     +    a1*(x(2)-x(3)) +    a2*(x(2)-x(3))^2;
% F1 = y(1)     +    a1*(x(1)-x(3)) +    a2*(x(1)-x(3))^2;
% 
if j == N-2
    
        a1 = [  (dxjp2 + dxjp1 + dxjm1)/(2*dxjm2*(dxjp2-dxjm2))   ]*f(j-2) + ...
...
            [   (-1)/(2*(dxjp1-dxjm1))  ]*f(j-1) + ...
...
                [ -(dxjm2+dxjm1+dxjp1+dxjp2)/(2*dxjm2*dxjp2) ]*f(j) + ...
...
                    [  (1)/(2*(dxjp1-dxjm1))   ]*f(j+1) + ...
...
                        [  -(dxjp1+dxjm1+dxjm2)/(2*dxjp2*(dxjp2-dxjm2))   ]*f(j+2) ;



        a2 = [-1/(dxjm2*(dxjp2-dxjm2))]*f(j-2) + [1/(dxjm2*dxjp2)]*f(j) + [1/(dxjp2*(dxjp2-dxjm2))]*f(j+2);

        FN_1 = a0   +   a1*(x(N-1)-x(N-2)) +    a2*(x(N-1)-x(N-2))^2;
        FN   = a0   +   a1*(x(N)-x(N-2)) +      a2*(x(N)-x(N-2))^2;

end
end


F(1) = F1; F(2) = F2;
F(end-1) = FN_1; F(end) = FN;

yf = F; 
end


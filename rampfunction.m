function [r_t] = rampfunction(tx,T)
    if tx<=2*T
       r_t = (1-cos(pi*tx/(2*T)))./2; 
    else
        r_t = 1;
    end
    
%     r_t = 0;
end
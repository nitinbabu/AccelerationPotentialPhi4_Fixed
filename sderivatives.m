function [U,V,NX,NY] = sderivatives(xf,yf,x,y,p,dpdn,qnode)
pf = p(1:length(xf)); qf = dpdn(1:length(xf)); % get fs with double nodes


%% Prefer to keep odd number of nodes on the free surface  
M = length(xf)/2-0.5;
fnode = qnode(1:M,:); % Free surface nodal address
%% xz and yz starting and ending there seems to be a problem. not due to corner problem
zeta = [-1 0 1];

for M = 1:M
    for n = 1:3
            xi  = zeta(n);
            xnode(M,n) = -(xi/2)*(1-xi)*xf(fnode(M,n)) +  (1+xi)*(1-xi)*xf(fnode(M,n)) + (xi/2)*(1+xi)*xf(fnode(M,n)) ;
            ynode(M,n) = -(xi/2)*(1-xi)*yf(fnode(M,n)) +  (1+xi)*(1-xi)*yf(fnode(M,n)) + (xi/2)*(1+xi)*yf(fnode(M,n)) ;
            pnode(M,n) = -(xi/2)*(1-xi)*pf(fnode(M,n)) +  (1+xi)*(1-xi)*pf(fnode(M,n)) + (xi/2)*(1+xi)*pf(fnode(M,n)) ;
            dpnode(M,n) = -(xi/2)*(1-xi)*qf(fnode(M,n)) +  (1+xi)*(1-xi)*qf(fnode(M,n)) + (xi/2)*(1+xi)*qf(fnode(M,n)) ;
    end
end

for M = 1:M
    for n = 1:3
            xi  = zeta(n);
            dxe(M,n) = xi*(xnode(M,1) -2*xnode(M,2)+xnode(M,3))    +   (xnode(M,3) - xnode(M,1))*0.5;
            dye(M,n) = xi*(ynode(M,1) -2*ynode(M,2)+ynode(M,3))    +   (ynode(M,3) - ynode(M,1))*0.5;
            dpe(M,n) = xi*(pnode(M,1) -2*pnode(M,2)+pnode(M,3))    +   (pnode(M,3) - pnode(M,1))*0.5;
            jacobian(M,n) = sqrt((dxe(M,n))^2+dye(M,n)^2);
            mx(M,n) = dxe(M,n)/jacobian(M,n); my(M,n) = dye(M,n)/jacobian(M,n);
%             nx(M,n) = dye(M,n)/jacobian(M,n); ny(M,n) = -dxe(M,n)/jacobian(M,n);
                    % Mistake ny-nx sign change fro becker to outward not
                    % inwards
            nx(M,n) = -dye(M,n)/jacobian(M,n); ny(M,n) = dxe(M,n)/jacobian(M,n);

    end
end



%% Reaffirm values
for M = 1:M
    for n = 1:3
        Q = fnode(M,n); xz(Q) = xnode(M,n); yz(Q) = ynode(M,n); pz(Q) = pnode(M,n);
        
        dxxi(Q) = dxe(M,n); dyxi(Q) = dye(M,n);    dpxi(Q) = dpe(M,n);
        nxc(Q) = nx(M,n);   nyc(Q) = ny(M,n); % (Changed before)Grilli 1989 Eq 31
%         dpds(Q) = dpe(M,n); % Not correct-placeholder
    end
end


%***********************************************************
% %%  Using the local derivatives %        Ning & Teng
for m = 1:length(xf)
A = [dxxi(m) dyxi(m);nxc(m) nyc(m)];
b = [dpxi(m); qf(m)];
Vvec = A\b;
U(m) = Vvec(1,:);   V(m) = Vvec(2,:);
end
NX = nxc;   NY = nyc;

%*********************************************************



end





































%% Last FS element if FS has even number of nodes. Avoid it!! If using use direcly before %% Reaffirm values section
% for M = length(xf)/2
%     for n = 1:2
%             xi  = zeta(n);
%             xnode(M,n) = -(xi/2)*(1-xi)*xf(fnode(M,n)) +  (1+xi)*(1-xi)*xf(fnode(M,n)) + (xi/2)*(1+xi)*xf(fnode(M,n)) ;
%             ynode(M,n) = -(xi/2)*(1-xi)*yf(fnode(M,n)) +  (1+xi)*(1-xi)*yf(fnode(M,n)) + (xi/2)*(1+xi)*yf(fnode(M,n)) ;
%     end
% end
% 
% for M = length(xf)/2
%     for n = 1:2
%             xi  = zeta(n);
%             dxe(M,n) = xi*(xnode(M-1,2) -2*xnode(M,1)+xnode(M,2))    +   (xnode(M,2) - xnode(M-1,2))*0.5;
%             dye(M,n) = xi*(ynode(M-1,2) -2*ynode(M,1)+ynode(M,2))    +   (ynode(M,2) - ynode(M-1,2))*0.5;
%             jacobian(M,n) = sqrt((dxe(M,n))^2+dye(M,n)^2);
%             mx(M,n) = dxe(M,n)/jacobian(M,n); my(M,n) = dye(M,n)/jacobian(M,n);
%             nx(M,n) = dye(M,n)/jacobian(M,n);   ny(M,n) = -dxe(M,n)/jacobian(M,n);
%     end
% end


% clear xz yz nxc nyc

%% Contrived recovery of values of FS has even number of nodes.
%% Even number of nodes means even nodes are to be arranged in M 3 noded elements hence last element will only have 2 nodal values to fill in
%% Not a problem if FS is odd number of nodes with/without double nodes
% 
% xz = [xz xnode(length(xf)/2,2)];
% yz = [yz ynode(length(xf)/2,2)];
% nxc = [nxc nx(length(xf)/2,2)];
% nyc = [nyc ny(length(xf)/2,2)];

% xz = xz(1:end-1);
% yz = yz(1:end-1);
% nxc = nxc(1:end-1);
% nyc = nyc(1:end-1);






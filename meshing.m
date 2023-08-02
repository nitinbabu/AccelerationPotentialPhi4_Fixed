function [qnode] = meshing(nelem)
qnode =zeros(nelem,3);
for n = 1:nelem
     if n ~= nelem
            for i = 1:3
                qnode(n,1) = 2*(n-1)+1;
                qnode(n,2) = 2*(n-1)+2;
                qnode(n,3) = 2*(n-1)+3;
            end
     else
            for i = 1:3
                qnode(n,1) = 2*(n-1)+1;
                qnode(n,2) = 2*(n-1)+2;
                qnode(n,3) = 1;%2*(n-1)+3; % Closing final node for geometry
            end
     end
end
function [C1, C2, C3, C4] = corners(xf,xd,xb,xp)
G1 = 1:length(xf);              % Free surface | Dirichlet Segment
G2 = G1(end)+(1:length(xd));    % Downstream   | Neumann Segment
G3 = G2(end)+ (1:length(xb));   % Bottom       | Dirichlet Segment (Laplace)
% G3 = G2(end)+ (1:length(xb)); % Bottom       | Neumann Segment
G4 = G3(end)+ (1:length(xp));   % Upstream     | Dirichlet Segment



C1 = G1(1); C2 = G1(end);   C3 = G3(1) ; C4 = G3(end);
end
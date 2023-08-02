function [C1, C2, C3, C4, C5, C6] = bodycorners(xfl,xbody,xfr,xd,xb,xp)
G1 = 1:length(xfl);              % Free surface L | Dirichlet Segment
G2 = G1(end)+(1:length(xbody));    % BODY   | Neumann Segment
G3 = G2(end)+(1:length(xfr));    % Free surface R | Dirichlet Segment

G4 = G3(end)+(1:length(xd));    % Downstream   | Neumann Segment
G5 = G4(end)+ (1:length(xb));   % Bottom       | Dirichlet Segment (Laplace)
% G3 = G2(end)+ (1:length(xb)); % Bottom       | Neumann Segment
G6 = G5(end)+ (1:length(xp));   % Upstream     | Dirichlet Segment



C1 = G1(1); C2 = G1(end);   C3 = G3(1) ; C4 = G3(end); C5 = G5(1); C6 = G5(end);
end
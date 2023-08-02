
function [cx,cy]=CenterMass(x,y)
% Define the x and y coordinates of the closed curve points

% Append the first point at the end to close the curve
x(end+1) = x(1);
y(end+1) = y(1);

% Calculate the centroid using the formula for the centroid of a polygon
n = length(x);
area = 0;
cx = 0;
cy = 0;

for i = 1:(n-1)
    term = (x(i) * y(i+1)) - (x(i+1) * y(i));
    area = area + term;
    cx = cx + (x(i) + x(i+1)) * term;
    cy = cy + (y(i) + y(i+1)) * term;
end

area = 0.5 * area;
cx = cx / (6 * area);
cy = cy / (6 * area);



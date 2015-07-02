function contained = contains(e1,e2)
%% Determing whether one ellipse, e1, is fully contained within another, e2
%% e1, e2: lists of ellipse descriptors: [xpos, ypos, major, minor rot]
%% Calculate Coordinate Axes
% for e1
i1 = [cosd(e1(5)), sind(e1(5))];
j1 = [-i1(2),i1(1)];
% for e2
i2 = [cosd(e2(5)), sind(e2(5))];
j2 = [-i2(2),i2(1)];
%% Calculate Extrema of Perimeter Points for e1
p = [    i1*e1(3);
        -i1*e1(3);
         j1*e1(4);
        -j1*e1(4)];
% Covert to Global System Origin
p(:,1) = p(:,1) + e1(1);
p(:,2) = p(:,2) + e1(2);
%% Transform to Axes of e2
p(:,1) = p(:,1) - e2(1);
p(:,2) = p(:,2) - e2(2);
% Calculate Transformation Matrix
T = [i2', j2'];
% Transform
p = p * T;
%% Calculate Relative Distance of Points to Perimeter of e2
Q = p(:,1).^2 ./ e2(3)^2  +  p(:,2).^2 ./ e2(4)^2;

% The below criteria is neccessary and, almost, sufficient
if sum(Q <= 1) == 4
    contained = true;
else
    contained = false;
end
end
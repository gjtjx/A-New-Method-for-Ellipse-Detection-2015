function B = opterode(A,se)
%% Modified Version of Matlabs Morphop Code
% Specific for erosion, assumes one flat SE & greyscale or binary image
%% Padding
[m,n] = size(A);
[sm,sn] = size(se);
B = zeros(m+(2*sm),n+2*sn);
B(sm+(1:m),sn+(1:n)) = A;

%% Apply the sequence of dilations/erosions.
height = zeros([sm,sn]);
if islogical(A)
    B = logical(B);
    B = images.internal.morphmex('erode_binary_twod', B, se, height, -1);
else
    B = images.internal.morphmex('erode_gray_flat', B, se, height, -1);
end

%% Extract the "middle" of the result; it should be the same size as the input image.
B = B(sm+(1:m),sn+(1:n));
end
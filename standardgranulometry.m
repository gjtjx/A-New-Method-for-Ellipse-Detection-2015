radii=200;
sek = cell(radii+1,1);
%% Strict Definition
% se = strel('disk',1);
% sek{1} = 1;
% sek{2} = se.getnhood;
% for r=2:radii
%     sek{r+1} = imdilate(sek{r},se,'full');
% end
%% MATLAB Definition
%delete the fourth one!
for r=0:radii
    sek{r+1} = strel('disk',r,0).getnhood;
end
%% Opening
intensity_area = zeros(radii,1);
for r=0:radii-1
    remain = imopen(im,sek{r+1});
    intensity_area(r+1) = sum(remain(:));
end
intensity_area_prime = diff(intensity_area);
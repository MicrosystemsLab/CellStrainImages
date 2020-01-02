function watershedAvgIm(FNS,cellNperLength)
% WATERSHEDAVGIM(FNS,CELLNPERLENGTH) finds cell boundaries in a series of
%  microscope images of a cell monolayer.
%  

fprintf(1,'\nwatershedAvgIm\n');

% load results from previous steps
load('tform.mat','tformCell')
%watershedSelected = {};
%waterImTSeries = {};

% load the fixed ("zero-strain") image
fixedAd = imgprocess2(FNS{1},1);


%% create "average image"
% Create an image which is the average of all the images in the set:
%  The "fixed image", and all the stretch images with the inverse transform
%  applied. This image will be used to find cell boundaries.

% start average registered images
avgImg = fixedAd/size(FNS,2);

fprintf(1,'watershedAvgIm: read all images and average\n');
% loop over all images
for i = 2:size(FNS,2)
	% read the next image
    fprintf(1,' Image: %d / %d: %s\n',i,size(FNS,2),FNS{i});
	movingAd = imgprocess2(FNS{i},1);
	
	% apply transform to this moving image
	tform = tformCell{i};
	movingRegistered = imwarp(movingAd,tform,'OutputView',imref2d(size(fixedAd)));

	% add result to average registered image
	avgImg = avgImg + movingRegistered/size(FNS,2);
end
save('avgImg.mat','avgImg')


%% apply watershed algorithm to find cell boundaries
rawIm = avgImg;
I2 = im2uint16(rawIm);
adjIm = imadjust(I2);
adaIm = adapthisteq(adjIm);
objWidth = 100;
regMinThreshold = 0.09;
% Convert to double-precision and re-scale.
rawIm = im2double(adaIm);
%rawIm = im2double(avgImg);
% Calculate local signal-to-noise ratios.
%  (In a way this corrects for uneven illumination)
snrIm = rawIm ./ imfilter(rawIm, fspecial('average', objWidth),'replicate');
% Enhance contrast by adaptive histogram equalization.
eqIm = adapthisteq(mat2gray(snrIm));
% Filter speckle noise using a 3x3 median filter.
medIm = medfilt2(eqIm,[3 3]);
% Suppress regional minima to mark the basins for the watershed algorithm.
% The cell-cell contacts are supposed to be the ridges.
imposedIm = imhmin(medIm, regMinThreshold);
figure, imshow(imposedIm,'InitialMagnification','fit')
title('Watershed Input')
% Apply watershed algorithm.
waterIm = watershed(imposedIm);

% Create binary mask of cell-cell contacts.
bwIm = waterIm == false;
bwIm = imdilate(bwIm, strel('sq', 3));
% Display results.
[m, n] = size(bwIm);
greenIm = cat(3, zeros(m, n), bwIm, zeros(m, n));
fusedIm = imfuse(rawIm, greenIm, 'blend');
% rgb = label2rgb(waterIm,'jet',[.5 .5 .5]);


%% filter cells based on size
% WARNING: Uses magic numbers based on pixels/cell 
AreaIm = regionprops(waterIm,'Area','BoundingBox','Centroid');
Areas = cat(1,AreaIm.Area);
cellSize = n*m/cellNperLength^2;
fprintf(1,'watershedAvgIm: using typical cell area of %.0f pixels^2\n',cellSize);
c_lower = cellSize/2.84; % lower limit of cell size
c_upper = cellSize*4.3;  % upper limit of cell size
fprintf(1,'watershedAvgIm: ignoring cells < %.0f and > %0.f pixels^2\n',c_lower,c_upper);
ind = find([AreaIm.Area] >= c_lower & [AreaIm.Area] <= c_upper);
% show cell boundaries
figure, imshow(fusedIm,'InitialMagnification','fit')
title(sprintf('Cell Boundaries (%d cells)',length(AreaIm)))
fprintf(1,'watershedAvgIm: identified %d cells\n',length(ind));
% show size distribution
figure
histogram(Areas,500);
hold on
grid on
plot(gca,[c_lower c_lower],ylim(gca),'-r')
plot(gca,[c_upper c_upper],ylim(gca),'-r')
title('Cell Size Histogram')

% show filtered cells
Iout = ismember(waterIm,ind);
rgb = label2rgb(Iout,'jet',[.5 .5 .5]);
figure, imshow(rgb,'InitialMagnification','fit')
title(sprintf('Cells Identified (%d cells)',length(ind)))
drawnow;


%% save results
watershedSelected = ind;
waterImTSeries = waterIm;


save('waterImTSeries.mat','waterImTSeries','watershedSelected')


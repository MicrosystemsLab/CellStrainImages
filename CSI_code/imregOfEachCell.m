function imregOfEachCell(FNS)
% IMREGOFEACHCELL(FNS) estimates strain in each cell by comparing it to the
%  the same cell in the zero-strain image.

fprintf(1,'\nimregOfEachCell\n');
warning('off','images:regmex:registrationFailedException');

% load results from previous steps
load('waterImTSeries.mat','watershedSelected','waterImTSeries')
load('tform.mat','tformCell')
ind = watershedSelected;
waterIm = waterImTSeries;

% load the fixed ("zero-strain") image
fixedAd = imgprocess2(FNS{1},2);

%% estimate strain for each cell in each image
%

TFORM = {};
fprintf(1,'imregOfEachCell: register cells in each image file\n');
% loop over all images
for i = 2:size(FNS,2)

	tstart = now;
	
	% read the next image
	fprintf(1,' Image: %d / %d start at %s: %s\n',i,size(FNS,2),datestr(tstart),FNS{i});	
	movingAd = imgprocess2(FNS{i},2);
    
	movingRegistered = imwarp(movingAd,tformCell{i},'OutputView',imref2d(size(fixedAd)));
% 	figure;
% 	imshowpair(fixedAd, movingRegistered,'Scaling','joint');
% 	title(FNS{i},'Interpreter','none')

	%% find deformation of each cell
	fprintf(1,' local imregistration for %d cells\n',length(ind));
	wb = waitbar(0,sprintf('%s (%d)',FNS{i},length(ind)),'Name',sprintf('Register Cells (%d / %d)',i,size(FNS,2)));
	wb.Children.Title.Interpreter = 'none';
	% loop over cells
	for j = 1:length(ind)
		
		% load indexed image ROI
		selectedReg = im2uint8(waterIm == ind(j));

		% make binary mask matrix
		level = graythresh(selectedReg);
		sbw = im2bw(selectedReg,level); %#ok<IM2BW>

		% dilate mask matrix
		se = strel('diamond',18);
		diIm = imdilate(sbw,se);

		% Make 1s-element matrix of doubles
		diImd = im2double(diIm);
		diImd = diImd/max(max(diImd));

		% find rectangle bounding box
		s = regionprops(diImd,'BoundingBox');
		cFixedIm = imcrop(fixedAd,s.BoundingBox);
		cMovingIm = imcrop(movingRegistered,s.BoundingBox);

		% imregister
		[optimizer, metric] = imregconfig('multimodal');
		optimizer.InitialRadius = 0.0005;
		optimizer.Epsilon = 1.5e-8;
		optimizer.GrowthFactor = 1.005;
		optimizer.MaximumIterations = 500;

		tform = imregtform(cMovingIm,cFixedIm,'affine',optimizer,metric);
		TFORM{i,j} = tform.T; %#ok<AGROW>
		% rmovingIm = imwarp(cMovingIm,tform,'OutputView',imref2d(size(cMovingIm)));
		% figure, imshowpair(rmovingIm,cFixedIm)

		waitbar(j/length(ind),wb)
		%drawnow;

	end
	close(wb)
	tend = now;
	fprintf(1,'imregOfEachCell: Loop elapsed time: %.0f seconds. End at %s\n',(tend-tstart)*86400,datestr(tend));
	drawnow;
	
end
fprintf(1,'imregOfEachCell: analyzed %d images.\n',length(tformCell)-1);
save('tformLocals.mat','TFORM')


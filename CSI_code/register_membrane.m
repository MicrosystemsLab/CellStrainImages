function register_membrane(FNS)
% REGISTER_MEMBRANE(FNS) registers the deformation of a series of
%  microscope images relative to the first image in the series.
%

% values related to registration quality
REGIST_QUALITY_COEFF = 150;
Rm = [];
regqual0 = 0;


fprintf(1,'\nregister_membrane\n');

tformCell = {};

if exist('tformAdj.mat', 'file') == 2
	fprintf(1,'register_membrane: Using manual adjustment values\n');
    load('tformAdj.mat','tformCell')
    [m,n] = size(tformCell) %#ok<NOPRT>
    % n is the size of the tformCell, upto which images are processed
else
    n = 1;
end


% load the fixed ("zero-strain") image
fprintf(1,'register_membrane: zero-strain image is "%s"\n',FNS{1});
fixedAd = imgprocess2(FNS{1},2);

figure, imshow(fixedAd,'InitialMagnification','fit');
title(['zero strain: ' FNS{1}],'Interpreter','none')
drawnow;

%% process all images
% Determine the deformation of each image relative to the first
%  ("zero-strain") image. This provides a good estimate of average
%  strain in the monolayer, but does not estimate strain in each cell.

% save strains estimated from image
imexx = zeros(1,size(FNS,2));
imeyy = zeros(1,size(FNS,2));

% loop over all images
fprintf(1,'register_membrane: register each image\n');
for i = n+1:size(FNS,2)
	tloop = now;
	% load an image
	fprintf(1,'Image: %d / %d: %s ',i,size(FNS,2),FNS{i});
	movingAd = imgprocess2(FNS{i},2);
	
	% use imregtform to estimate the distortion of the image relative to
	%  the zero-strain image
	if n == 1
		[optimizer, metric] = imregconfig('multimodal');
		optimizer.InitialRadius = 0.0004;
		optimizer.Epsilon = 1.5e-8;
		optimizer.GrowthFactor = 1.005;
		optimizer.MaximumIterations = 500;

		if isempty(Rm)
			tformi = imregtform(movingAd,imref2d(size(movingAd)),fixedAd,imref2d(size(fixedAd)),'affine',optimizer,metric);
		else
			tformi = imregtform(movingAd,Rm,fixedAd,imref2d(size(fixedAd)),'affine',optimizer,metric);
		end
		
	else
		tformInit = tformCell{2};
		movingRegInit = imwarp(movingAd,tformInit,'OutputView',imref2d(size(fixedAd)));

		[optimizer, metric] = imregconfig('multimodal');
		optimizer.InitialRadius = 0.0004;
		optimizer.Epsilon = 1.5e-8;
		optimizer.GrowthFactor = 1.005;
		optimizer.MaximumIterations = 500;
		
		tformInter = imregtform(movingRegInit,fixedAd,'affine',optimizer,metric);
		tformi = tformInit;
		tformi.T = tformInter.T*tformInit.T;
	end

	tformCell{i} = tformi; %#ok<AGROW>
	imexx(i) = tformi.T(1,1);
	imeyy(i) = tformi.T(2,2);

	% Test if the registration is failed
	if i ~= n+1
		regqual = tformCell{i-1}.T/(tformCell{i}.T);
		idm = [1 0 0;0 1 0;0 0 1];
		regqual0 = [ones(2,3)*1e3;ones(1,3)].*abs(idm - regqual);
	end
	if sum(sum(regqual0)) <= REGIST_QUALITY_COEFF  % if NOT failed
		if isempty(Rm)
			[movingRegistered,Rm] = imwarp(movingAd,imref2d(size(fixedAd)),tformi,'OutputView',imref2d(size(fixedAd)));
		else
			[movingRegistered,Rm] = imwarp(movingAd,Rm,tformi,'OutputView',imref2d(size(fixedAd)));
		end
	else % try again
		fprintf(1,'\n Re-trying image registration... ');
		tformInit = tformCell{i-1};
		movingRegInit = imwarp(movingAd,tformInit,'OutputView',imref2d(size(fixedAd)));

		[optimizer, metric] = imregconfig('multimodal');
		optimizer.InitialRadius = 0.0005;
		optimizer.Epsilon = 1.5e-8;
		optimizer.GrowthFactor = 1.006;
		optimizer.MaximumIterations = 600;

		tformInter = imregtform(movingRegInit,fixedAd,'affine',optimizer,metric);
		tformi = tformInit;
		tformi.T = tformInter.T*tformInit.T;
		tformCell{i} = tformi; %#ok<AGROW>
	end
	
%  	% optional: show comparison of the fixed and moving image
%  	imshowpair(fixedAd, movingRegistered,'Scaling','joint');
% 	title(sprintf('Image Pair %d',i),'Interpreter','none')
% 	drawnow;
    
	fprintf(1,'(%.0f sec)\n',(now-tloop)*86400);
end

drawnow;
fprintf(1,'register_membrane: registered %d images.\n',length(tformCell)-1);

% estimate strains from image registration
imexx = 1./imexx-1; imexx(1) = 0;
imeyy = 1./imeyy-1; imeyy(1) = 0;
pri = 0:1:size(FNS,2)-1;
figure
plot(pri,imexx,'o-',pri,imeyy,'s-');
title('Strains from Image Registration')
grid on

% save results
save('tform.mat','tformCell')
save('image_reg_strains.mat','imexx','imeyy','pri')


function [If,testOut] = segmentCellBackground(img,background_seg,pStruct,frames)
testOut=struct();
%         channel = alterChanName(channel);
nucDiameter = pStruct.(background_seg).nucDiameter;
threshFactor = pStruct.(background_seg).threshFactor;
sigmaScaledToParticle = pStruct.(background_seg).sigmaScaledToParticle;
fnames = fieldnames(pStruct.(background_seg));
if sum(strcmpi(fnames,'denoise'))>0
    wienerP=pStruct.(background_seg).denoise;
else
    wienerP=5;
end
kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize

percentSmoothed = pStruct.(background_seg).percentSmoothed;
psm = 100-percentSmoothed;

ogimg = img;
            imgRaw = img;
            imgRaw(imgRaw>prctile(imgRaw(:),psm)) = prctile(imgRaw(:),psm);
            img = regionfill(img,img>prctile(imgRaw(:),psm));
efiltsize = round(nucDiameter./10);
if mod(efiltsize,2)==0
    efiltsize = efiltsize+1;
end
% output = entropyfilt(uint16(img),true(efiltsize));
output = stdfilt(img,true(efiltsize));
output(output>prctile(output(:),90)) = prctile(output(:),90);
[~,~,~,threshLocation] = method4(output(:));
outputlog = false(size(output));
outputlog(:) = false;
outputlog(output<(threshLocation*threshFactor))=true;
Ih = ~outputlog;
areaOfSegmentation = sum(Ih(:));
percentageOfImageSegmented = round(100*(areaOfSegmentation./(size(img,1)*size(img,2))));
if percentageOfImageSegmented < 95
    %     Ih = imdilate(Ih,strel('disk',floor(nucDiameter./2)));
end

%USE IF IMAGE IS NOT FLAT
% imgLowPass = gaussianBlurz(img,sigma,kernelgsize);
% intimg = regionfill(imgLowPass,Ih);
% flatten = intimg./prctile(intimg(:),95);
% flatten = wiener2(flatten,[wienerP wienerP]);
% flatten(flatten>1.2) = 1.2;
% img = img./flatten;

% flatten = intimg;
% flatten = wiener2(flatten,[wienerP wienerP]);
% % flatten(flatten>1.2) = 1.2;
% img = img-flatten;


%%
%identify background pixels by finding and grouping low fluorescence pixels in the non
%segmented area
fillIh = imfill(Ih,'holes');
newlog = img<prctile(img(~fillIh),20); %mask -->nonsegmentated areas, 20th prctile fluorescence
newlogclose = imclose(newlog,strel('disk',20)); %mask -->cluster the regions together
newlogclose(fillIh) = false; %mask -->make sure clustered regions don't overlap with original segmentation

%this should identify areas above background
%finds pixels in img whose fluorescence is higher than background pixels
%(as identified by nelogclose mask previously)
closertest = false(size(newlogclose));
closertest(img>prctile(img(newlogclose),50)) = true; %mask -->areas of img greater than 50% prctile of prevmask areas

%group non-cell areas (imclose of inverted cell areas), 
%then open/erode cell areas (inverted image)
%then dilate the cell areas
closertestclose = imdilate(imopen(~imclose(~closertest,strel('disk',2)),strel('disk',2)),strel('disk',5));

%find pixels greater 20 prctile of segmented area and set them false (these
%are incorrect background areas)
%group the true regions of newlogclose (which is the background area pixels)
%then make sure the areas identified as cells above are excluded
% newlogclose(img>prctile(img(Ih),20)) = false;
newlogclose(img>prctile(img(closertestclose),20)) = false;
newlogcloseclose = imclose(newlogclose,strel('disk',5)); %group the true regions of newlogclose (which is the background area pixels)
newlogcloseclose(closertestclose) = false;
If = ~newlogcloseclose;


% f = figure(9999);
% f.Position =[495 204 1200 663];
% subplot(2,4,1);imagesc(ogimg);title('ogimg')
% subplot(2,4,2);imagesc(img);title('img')
% subplot(2,4,3);imagesc(Ih);title('Ih')
% subplot(2,4,4);imagesc(flatten);title('flatten')
% subplot(2,4,5);imagesc(newlog);title('newlog')
% subplot(2,4,6);imagesc(closertestclose);title('closertestclose')
% subplot(2,4,7);imagesc(newlogclose);title('newlogclose')
% subplot(2,4,8);imagesc(If);title('If')

%%

% If = ~Ihf;

%
% figure;
% subplot(2,3,1);imagesc(ogimg);
% subplot(2,3,2);imagesc(img);
% subplot(2,3,3);imagesc(flatten);
% subplot(2,3,4);imagesc(Ih);
% subplot(2,3,4);imagesc(newlog);
% subplot(2,3,5);imagesc(newlogclose);
% subplot(2,3,5);imagesc(closertest);


Im= zeros(size(img));

if frames==1
    testOut.img = img;
    testOut.I = -1.*Im;
    testOut.imgRawDenoised = Im;
    testOut.imgLowPass = Im;
    testOut.rawMinusLP = Im;
    testOut.rawMinusLPScaled = Im;
    testOut.Ih = Ih;
    testOut.Inew = zeros([512 512]);
    testOut.Im = Im;
    %         testOut.Ihcd = Ihcd;
    testOut.L = zeros([512 512]);
    %         testOut.gradmag = gradmag;
    testOut.gradmag = zeros(size(img));
    %         testOut.gradmag2 =  gradmag2;
    testOut.gradmag2 = zeros(size(img));
    %         testOut.Ie = Ie;
    testOut.Ie = zeros(size(img));
    %         testOut.fgm4 = fgm4;
    testOut.fgm4 = zeros(size(img));
    %         testOut.Ieg = Ieg;
    testOut.Ieg = zeros(size(img));
    testOut.Shapes = zeros(size(img));
    %         testOut.waterBoundary = waterBoundary;
    
end

end


function bw = gaussianBlurz(im,sigma,kernelgsize,varargin)

filtersize = [kernelgsize kernelgsize];
kernelg = fspecial('gaussian',filtersize,sigma);

gFrame = imfilter(im,kernelg,'repl');

if ~isempty(varargin)
    bw=gFrame.*uint16(varargin{1}>0);
else
    bw=gFrame;
end
end

function [sdfdone,fraction,bincenters,threshLocation]= method3(vec)
lowperc = prctile(vec,0.01);
highperc = prctile(vec,100);
[numbers,bincenters] = hist(vec,lowperc:(highperc-lowperc)/500:highperc);
%     numbersone = movmean(numbers,10,2,'Endpoints','shrink');
%     numberstwo = movmean(numbersone,100,2,'Endpoints','shrink');
%     numbersone = movmean(numbers,10,2);
%     numbersone = smooth(numbers,'rlowess');
numbersone = movmean(numbers,5,2,'Endpoints','fill')';
movmeanmagnitude= 50;
numberstwo = movmean(numbersone,movmeanmagnitude,1,'Endpoints','fill');
fraction = numberstwo./max(numberstwo);
diffFraction = diff(fraction,1,1);
diffFraction = diffFraction./max(diffFraction);
%     sdf = smooth(diffFraction,'lowess');
sdf = movmean(diffFraction,5,1,'EndPoints','fill');
%     sdf = movmean(diffFraction,5,1);
%     sdfdone = movmean(sdf,50,1);
sdfdone = movmean(sdf,5,1,'Endpoints','fill');
%     sdfdone(1:movmeanmagnitude)=NaN; sdfdone(end-movmeanmagnitude:end) =NaN;

leftedge = find(sdfdone == max(sdfdone),1,'first');
insideslopedown = find(sdfdone(leftedge:end) == min(sdfdone(leftedge:end)),1,'first');
slopeup = find(sdfdone(leftedge+insideslopedown-1:end) == max(sdfdone(leftedge+insideslopedown-1:end)),1,'first');

threshLocation=[];
%conditionals for determining threshLocation
if isempty(leftedge)
    stp=1;
elseif leftedge==1
    whatup=1;
elseif isempty(insideslopedown) && (sum(logvec==0)>100)
    threshLocation = bincenters(leftedge);
elseif isempty(slopeup) || slopeup==1
    threshLocation = bincenters(leftedge);
    threshFactor = 0.5;
elseif sdfdone(leftedge)<0.3
    threshLocation = bincenters(leftedge);
    threshFactor = 0.5;
elseif sdfdone(leftedge+insideslopedown-1)>0
    threshLocation = bincenters(leftedge);
else
    threshLocation = bincenters(leftedge+insideslopedown-1);
end

end


function [sdfdone,fraction,bincenters,threshLocation]= method4(vec)
lowperc = prctile(vec,0.01);
highperc = prctile(vec,100);
[numbers,bincenters] = hist(vec,lowperc:(highperc-lowperc)/500:highperc);


numbersone = movmean(numbers,5,2,'Endpoints','shrink')';
movmeanmagnitude= 50;
numberstwo = movmean(numbersone,movmeanmagnitude,1,'Endpoints','shrink');
fraction = numberstwo./max(numberstwo);
diffFraction = diff(fraction,1,1);
diffFraction = diffFraction./max(diffFraction);

sdf = movmean(diffFraction,5,1,'EndPoints','shrink');
sdfdone = movmean(sdf,5,1,'Endpoints','shrink');

firstrightedge = find(sdfdone == min(sdfdone),1,'first');
firsttrough = find(sdfdone(firstrightedge:end) > -0.15,1,'first')+firstrightedge-1;
firstleftedge = find(sdfdone(firsttrough:end) == max(sdfdone(firsttrough:end)),1,'first')+firsttrough-1;

%     figure(22)
%     bvec = bincenters(1:end-1);
%     plot(bvec,sdfdone);hold on
%     scatter(bvec(firstrightedge),sdfdone(firstrightedge))
%     scatter(bvec(firsttrough),sdfdone(firsttrough));hold off
% 
%     plot(bincenters,fraction);hold on
%     scatter(bincenters(firstrightedge),fraction(firstrightedge));
%     scatter(bincenters(firsttrough),fraction(firsttrough));
%     scatter(bincenters(firstleftedge),fraction(firstleftedge))


threshLocation = bincenters(firsttrough);
if isempty(threshLocation)
    threshLocation = bincenters(firstrightedge);
end
end



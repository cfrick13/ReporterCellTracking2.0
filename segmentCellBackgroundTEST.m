
function [If,testOut] = segmentCellBackground(img,background_seg,pStruct,frames)
testOut=struct();    
%         channel = alterChanName(channel);
nucDiameter = pStruct.(background_seg).nucDiameter;
threshFactor = pStruct.(background_seg).threshFactor;
sigmaScaledToParticle = pStruct.(background_seg).sigmaScaledToParticle;
kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize
percentSmoothed = pStruct.(background_seg).percentSmoothed;


%initial segmentation to determine how much of image is covered by cells
% img = FinalImage(:,:,1); 
imgW = wiener2(img,[1 20]);
imgWW = wiener2(imgW,[20 1]);
imgWWW = wiener2(imgWW,[5 5]);
imgRawDenoised = imgWWW;
% denoiseVec = single(reshape(imgRawDenoised,size(imgRawDenoised,1)^2,1));
% highpoints = prctile(imgRawDenoised(:),80);
% imgRawDenoised(imgRawDenoised>highpoints) = highpoints;
%
imgLowPass = gaussianBlurz(single(imgRawDenoised),sigma,kernelgsize);
rawMinusLP = single(imgRawDenoised) -single(imgLowPass);%%%%%%% key step!
rawMinusLPvec = reshape(rawMinusLP,size(rawMinusLP,1)^2,1);
globalMinimaValues = prctile(rawMinusLPvec,0.01);
globalMinimaIndices = find(rawMinusLP < globalMinimaValues);
LPscalingFactor = imgRawDenoised(globalMinimaIndices)./imgLowPass(globalMinimaIndices);
imgLPScaled = imgLowPass.*nanmedian(LPscalingFactor);
rawMinusLPScaled = single(imgRawDenoised) - single(imgLPScaled);


%rescale the image
    lcontrast = 0;
    tcontrast = 100;
    lprcntl = prctile(rawMinusLPScaled(:),lcontrast);
    prcntl = prctile(rawMinusLPScaled(:),tcontrast);
    scaleFactor = 1./(prcntl - lprcntl);
    rawMinusLPScaledContrasted = rawMinusLPScaled.*scaleFactor;
    rawMinusLPScaledContrasted = rawMinusLPScaledContrasted-(lprcntl.*scaleFactor);
                
vecOG = single(reshape(rawMinusLPScaledContrasted,size(rawMinusLPScaledContrasted,1)^2,1));
logvecpre = vecOG; logvecpre(logvecpre==0)=[];
logvec = log10(logvecpre);
vec = logvec;

[~,~,~,threshLocation] = method3(vec);
subtractionThreshold = threshLocation;

if size(subtractionThreshold,1)==size(subtractionThreshold,2)
    else
     subtractionThreshold = mean(threshLocation);
end
subtractionThresholdScaled = (10.^subtractionThreshold).*threshFactor;
subtracted = single(rawMinusLPScaledContrasted)-subtractionThresholdScaled;
subzero = (subtracted<0);
Ih = ~subzero;
    
    imgW = wiener2(img,[1 20]);
    imgWW = wiener2(imgW,[20 1]);
    imgWWW = wiener2(imgWW,[5 5]);
    imgRawDenoised = imgWWW;
    highpoints = prctile(imgWWW(Ih),percentSmoothed);
    imgRawDenoised(imgRawDenoised>highpoints) = highpoints;

    If = imgRawDenoised;
    mmIf = max(If(:)) ;
    If(If<mmIf)=0;
    If(If == mmIf)=1;
    If = logical(If);
    
    areaOfSegmentation = sum(If(:));
    percentageOfImageSegmented = round(100*(areaOfSegmentation./(size(img,1)*size(img,2))));
    if percentageOfImageSegmented > 99
        percentageOfImageSegmented = 99;
    elseif percentageOfImageSegmented == 0
        percentageOfImageSegmented = 1;
    end
    
    
    areamin = max([(100-percentageOfImageSegmented)./10 1]);
    areamax = max([(100-percentageOfImageSegmented)./2 1]);
    imgarea = (size(If,1).*size(If,2));

   width = 30;
   se = strel('disk',width);
   Ig = imdilate(If,se);
   a = ((imgarea-sum(Ig(:)))./imgarea).*100;
   while  a>areamax
       width = ceil(width*1.5);
       se = strel('disk',width);
       Ig = imdilate(If,se);
       a = ((imgarea-sum(Ig(:)))./imgarea).*100;
   end
   
   while  a<areamin
       width = ceil(width./2);
       se = strel('disk',width);
       Ig = imdilate(If,se);
       a = ((imgarea-sum(Ig(:)))./imgarea).*100;
   end

   if sum(~Ig(:))>0
       If = Ig;
   else
       If=If;
   end
   
   
  Im= zeros(size(imgW));
  
       if frames==1
        testOut.img = img;
        testOut.I = -1.*rawMinusLPScaled;
        testOut.imgRawDenoised = imgRawDenoised;
        testOut.imgLowPass = imgLowPass;
        testOut.rawMinusLP = rawMinusLP;
        testOut.rawMinusLPScaled = rawMinusLPScaled;
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

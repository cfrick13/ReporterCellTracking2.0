
function [If,testOut] = segmentCellBackground(img,background_seg,pStruct,frames)
testOut=struct();    
%         channel = alterChanName(channel);
nucDiameter = pStruct.(background_seg).nucDiameter;
threshFactor = pStruct.(background_seg).threshFactor;
sigmaScaledToParticle = pStruct.(background_seg).sigmaScaledToParticle;
kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize



%initial segmentation to determine how much of image is covered by cells
% img = FinalImage(:,:,1); 
imgW = wiener2(img,[1 20]);
imgWW = wiener2(imgW,[20 1]);
imgWWW = wiener2(imgWW,[5 5]);
imgRawDenoised = imgWWW;
denoiseVec = single(reshape(imgRawDenoised,size(imgRawDenoised,1)^2,1));
highpoints = prctile(denoiseVec,95);
imgRawDenoised(imgRawDenoised>highpoints) = highpoints;
%
imgLowPass = gaussianBlurz(single(imgRawDenoised),sigma,kernelgsize);
rawMinusLP = single(imgRawDenoised) -single(imgLowPass);%%%%%%% key step!
rawMinusLPvec = reshape(rawMinusLP,size(rawMinusLP,1)^2,1);
globalMinimaValues = prctile(rawMinusLPvec,0.01);
globalMinimaIndices = find(rawMinusLP < globalMinimaValues);
LPscalingFactor = imgRawDenoised(globalMinimaIndices)./imgLowPass(globalMinimaIndices);
imgLPScaled = imgLowPass.*nanmedian(LPscalingFactor);
rawMinusLPScaled = single(imgRawDenoised) - single(imgLPScaled);

rawMinusLPScaledvec = reshape(rawMinusLPScaled,size(rawMinusLPScaled,1)^2,1);
high_in = prctile(rawMinusLPScaledvec,99);
rawMinusLPScaledContrasted = imadjust(rawMinusLPScaled./high_in,[0.1; 0.99],[0; 1]);

vecOG = single(reshape(rawMinusLPScaledContrasted,size(rawMinusLPScaledContrasted,1)^2,1));
logvecpre = vecOG; logvecpre(logvecpre==0)=[];
logvec = log10(logvecpre);
vec = logvec;
[numbers,bincenters] = hist(vec,prctile(vec,1):(prctile(vec,99)-prctile(vec,1))/1000:max(vec));
numbersone = medfilt1(numbers, 10); %smooths curve
numberstwo = medfilt1(numbersone, 100); %smooths curve
fraction = numberstwo./sum(numberstwo);
mf = max(fraction);
    %%%%%%%%%%%%%%%%%%%% Important parameters for finding minima of
    %%%%%%%%%%%%%%%%%%%% histogram
    left=0.5*mf;
    slopedown=0.4*mf;
    %%%%%%%%%%%%%%%%%%%%%
leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
threshLocation = bincenters(leftedge+insideslopedown-1);
subtractionThreshold = threshLocation;

if size(subtractionThreshold,1)==size(subtractionThreshold,2)
    else
     subtractionThreshold = mean(threshLocation);
end
subtractionThresholdScaled = (10.^subtractionThreshold).*threshFactor;
subtracted = single(rawMinusLPScaledContrasted)-subtractionThresholdScaled;
subzero = (subtracted<0);
Ih = ~subzero;
Ih = imclose(Ih,strel('disk',20));
areaOfSegmentation = sum(sum(Ih));
%
percentageOfImageSegmented = round(100*(areaOfSegmentation./(size(img,1)*size(img,2))));
if percentageOfImageSegmented > 99
    percentageOfImageSegmented = 99;
elseif percentageOfImageSegmented == 0
    percentageOfImageSegmented = 1;
end

    
    imgW = wiener2(img,[1 20]);
    imgWW = wiener2(imgW,[20 1]);
    imgWWW = wiener2(imgWW,[5 5]);
    imgRawDenoised = imgWWW;
    denoiseVec = single(reshape(imgRawDenoised,size(imgRawDenoised,1)^2,1));
    highpoints = prctile(denoiseVec,percentageOfImageSegmented);
    imgRawDenoised(imgRawDenoised>highpoints) = highpoints;

    
    %Based on algorithm of Fast and accurate automated cell boundary determination for fluorescence microscopy by Arce et al (2013)   
    %LOW PASS FILTER THE IMAGE (scale the gaussian filter to diameter of
    %nuclei -- diameter of nuclei is about 50 to 60))
    
    imgLowPass = gaussianBlurz(single(imgRawDenoised),sigma,kernelgsize);
    rawMinusLP = single(imgRawDenoised) -single(imgLowPass);%%%%%%% key step!
    rawMinusLPvec = reshape(rawMinusLP,size(rawMinusLP,1)^2,1);
    globalMinimaValues = prctile(rawMinusLPvec,0.01);
    globalMinimaIndices = find(rawMinusLP < globalMinimaValues);
    LPscalingFactor = imgRawDenoised(globalMinimaIndices)./imgLowPass(globalMinimaIndices);
    imgLPScaled = imgLowPass.*nanmedian(LPscalingFactor);
    rawMinusLPScaled = single(imgRawDenoised) - single(imgLPScaled);


    %determine the threshold by looking for minima in log-scaled histogram
    %of pixels from rawMinusLPScaled
    rawMinusLPScaledvec = reshape(rawMinusLPScaled,size(rawMinusLPScaled,1)^2,1);
    high_in = prctile(rawMinusLPScaledvec,99);
    rawMinusLPScaledContrasted = imadjust(rawMinusLPScaled./high_in,[0.1; 0.99],[0; 1]);
    
    vecOG = single(reshape(rawMinusLPScaledContrasted,size(rawMinusLPScaledContrasted,1)^2,1));
    logvecpre = vecOG; logvecpre(logvecpre==0)=[];
    logvec = log10(logvecpre);
    vec = logvec;
    [numbers,bincenters] = hist(vec,prctile(vec,1):(prctile(vec,99)-prctile(vec,1))/1000:max(vec));
    numbersone = medfilt1(numbers, 10); %smooths curve
    numberstwo = medfilt1(numbersone, 100); %smooths curve
    fraction = numberstwo./sum(numberstwo);
    mf = max(fraction);
        %%%%%%%%%%%%%%%%%%%% Important parameters for finding minima of
        %%%%%%%%%%%%%%%%%%%% histogram
        left=0.5*mf;
        slopedown=0.4*mf;
        %%%%%%%%%%%%%%%%%%%%%
    leftedge = find(fraction > left,1,'first');
    insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
    threshLocation = bincenters(leftedge+insideslopedown-1);
    subtractionThreshold = threshLocation;

    if size(subtractionThreshold,1)==size(subtractionThreshold,2)
        else
         subtractionThreshold = mean(threshLocation);
    end


    subtractionThresholdScaled = (10.^subtractionThreshold).*threshFactor;
    if ~(length(subtractionThresholdScaled)==1)
        stophere=1;
    end
    subtracted = single(rawMinusLPScaledContrasted)-subtractionThresholdScaled;
    subzero = (subtracted<0);
    Ih = ~subzero;

    width = 10;
    Ihc = imclose(Ih,strel('disk',width));
    Im=Ihc;


    If = imgRawDenoised;
    mmIf = max(max(If)) ;
    If(If<mmIf)=0;
    If(If == mmIf)=1;
    If = logical(If);
    arealimit = (100-percentageOfImageSegmented)./8;
    imgarea = (size(If,1).*size(If,2));

   a = length(If==0);
   width = 10;
   Ig= If;
   while  1
       se = strel('disk',width);
       a = ((imgarea-sum(sum(Ig)))./imgarea).*100;
       if a<arealimit
           break
       else
            Ig = imdilate(Ig,se);
       end

   end

    If=Ig;
    
       if frames==1
        testOut.img = img;
        testOut.I = -1.*rawMinusLPScaled;
        testOut.imgRawDenoised = imgRawDenoised;
        testOut.imgLowPass = imgLowPass;
        testOut.rawMinusLP = rawMinusLP;
        testOut.rawMinusLPScaled = rawMinusLPScaled;
        testOut.Ih = Ih;
%         testOut.Ihc = Ihc;
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

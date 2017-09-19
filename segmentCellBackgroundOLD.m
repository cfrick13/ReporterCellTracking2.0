
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
Ih = imclose(Ih,strel('disk',20));
areaOfSegmentation = sum(Ih(:));
%
percentageOfImageSegmented = round(100*(areaOfSegmentation./(size(img,1)*size(img,2))));
while percentageOfImageSegmented > 99
    threshFactor = threshFactor*1.1;
    subtractionThresholdScaled = (10.^subtractionThreshold).*threshFactor;
    subtracted = single(rawMinusLPScaledContrasted)-subtractionThresholdScaled;
    subzero = (subtracted<0);
    Ih = ~subzero;
    Ih = imclose(Ih,strel('disk',20));
    areaOfSegmentation = sum(Ih(:));
    %
    percentageOfImageSegmented = round(100*(areaOfSegmentation./(size(img,1)*size(img,2))));
end

newlog = img<prctile(img(~Ih),0.5);
newlogclose = imclose(newlog,strel('disk',20));
newlogclose(Ih) = false;
newlogclose(img>prctile(img(Ih),0.5)) = false;
newlogcloseclose = imclose(newlogclose,strel('disk',5));
If = ~newlogcloseclose;



output = entropyfilt(uint16(img),true(21));
[~,~,~,threshLocation] = method4(output(:));
outputlog = false(size(output));
outputlog(:) = false;
outputlog(output<threshLocation)=true;
Ih = ~outputlog;


newlog = img<prctile(img(~Ih),0.5);
newlogclose = imclose(newlog,strel('disk',20));
newlogclose(Ih) = false;
newlogclose(img>prctile(img(Ih),0.5)) = false;
newlogcloseclose = imclose(newlogclose,strel('disk',5));
Ihf = ~newlogcloseclose;


% figure(9999)
% subplot(1,2,1);imagesc(If);
% subplot(1,2,2);imagesc(Ihf)

If = Ihf;

stophere=1;

%% identify regions where standard deviation of image is small
% output = zeros(size(img));
% csize = 10;
% for i = 1:size(img,1)
%     for j = 1:size(img,2)
%         x1 = max([1 i-csize]);
%         x2 = min([size(img,1) i+csize]);
%         y1 = max([1 j-csize]);
%         y2 = min([size(img,2) j+csize]);
%         
%         input = img(x1:x2,y1:y2);
%         compval = nanstd(input(:))./nanmean(input(:));
%         output(i,j) = compval;
%         
%     end
% end


% 
% 
% %try someting new
% %first inmeanvals
% inmean0 = prctile(img(~Ih),95);
% instd0 = nanstd(img(~Ih));
% 
% %next iteratively erode and dilate
% origlog = ~Ih;
% inmeanvals = [];
% instdvals =[];
% senumvec =[];
% sevec = 5:5:100;
% for senum = 1:length(sevec)
%     se = strel('disk',sevec(senum));
%     Ihe = imdilate(origlog,se);
%     sumihe = sum(~Ihe);
%     if sumihe<100
%         break
%     end
%     inmean = prctile(img(Ihe),95);
%     instd = nanstd(img(Ihe));
%     inmeanvals = [inmeanvals(:)' inmean(:)'];
%     instdvals = [instdvals(:)' instd(:)'];
%     senumvec = [senumvec(:)' -sevec(senum)'];
% end
% cycle1 = senum-1;
% inmean1 = inmeanvals;
% instd1 = instdvals;
% senum1 = senumvec;
% 
% origlog = ~Ih;
% inmeanvals = [];
% instdvals =[];
% senumvec =[];
% sevec = 5:5:100;
% for senum = 1:length(sevec)
%     se = strel('disk',sevec(senum));
%     Ihe = imerode(origlog,se);
%     sumihe = sum(Ihe);
%     if sumihe<100
%         break
%     end
%     inmean = prctile(img(Ihe),95);
%     instd = nanstd(img(Ihe));
%     inmeanvals = [inmeanvals(:)' inmean(:)'];
%     instdvals = [instdvals(:)' instd(:)'];    
%     senumvec = [senumvec(:)' sevec(senum)'];
% end
% cycle2=senum-1;
% inmean2 = inmeanvals;
% instd2 = instdvals;
% senum2 = senumvec;
% 
% 
% 
% meanvals = [inmean1(:)' inmean0 inmean2(:)'];
% stdvals = [instd1(:)' instd0 instd2(:)'];
% senumvals = [senum1(:)' 0 senum2(:)'];
% 
% 
% midx = meanvals < min(meanvals)*1.05;
% sidx = stdvals < min(stdvals)*1.05;
% nidx = midx & sidx;
% seidx = senumvals(nidx);
% 
% senum = seidx(1);
% if senum<0
%     se = strel('disk',abs(senum));
%     Ihnew = imdilate(origlog,se);
%     %imdilate
% else
%     se = strel('disk',senum);
%     Ihnew = imerode(origlog,se);
% end
%     
% If = Ihnew;
    
%     imgW = wiener2(img,[1 20]);
%     imgWW = wiener2(imgW,[20 1]);
%     imgWWW = wiener2(imgWW,[5 5]);
%     imgRawDenoised = imgWWW;
%     denoiseVec = single(reshape(imgRawDenoised,size(imgRawDenoised,1)^2,1));
%     highpoints = prctile(imgWWW(Ih),percentSmoothed);
% %     highpoints = prctile(denoiseVec,percentageOfImageSegmented);
% 
% 
% %what is this for?????
%     a = sum(imgRawDenoised(:)>highpoints)>(size(imgW,1)*size(imgW,2).*0.85);
%     stepup=1;
%     newperc = 1;
%     areacriteria = size(imgW,1)*size(imgW,2).*0.85;
%     while a==1
%         newperc = newperc+stepup;
%         highpoints = prctile(imgWWW(Ih),newperc);
%         a = sum(imgRawDenoised(:)>highpoints)>areacriteria;
%     end
%     imgRawDenoised(imgRawDenoised>highpoints) = highpoints;
% 
%  
%     
% 
% 
%     If = imgRawDenoised;
% %     If = Im;
%     mmIf = max(If(:)) ;
%     If(If<mmIf)=0;
%     If(If == mmIf)=1;
%     If = logical(If);
%     arealimit = min([(100-percentageOfImageSegmented) 5]);
%     imgarea = (size(If,1).*size(If,2));
% 
%    a = length(If==0);
%    width = 10;
%    Ig= If;
%    while  1
%        se = strel('disk',width);
%        a = ((imgarea-sum(Ig(:)))./imgarea).*100;
%        if a<arealimit
%            break
%        else
%             Ig = imdilate(Ig,se);
%        end
% 
%    end
% 
%    if sum(~Ig(:))>0
%        If = Ig;
%    else
%        %If=If;
%    end
%    
   
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


function [sdfdone,fraction,bincenters,threshLocation]= method4(vec)
    lowperc = prctile(vec,0.01);
    highperc = prctile(vec,100);
    [numbers,bincenters] = hist(vec,lowperc:(highperc-lowperc)/500:highperc);

    
    numbersone = movmean(numbers,5,2,'Endpoints','fill')';
    movmeanmagnitude= 50;
    numberstwo = movmean(numbersone,movmeanmagnitude,1,'Endpoints','fill');
    fraction = numberstwo./max(numberstwo);
    diffFraction = diff(fraction,1,1);
    diffFraction = diffFraction./max(diffFraction);

    sdf = movmean(diffFraction,5,1,'EndPoints','fill');
    sdfdone = movmean(sdf,5,1,'Endpoints','fill');
    
    firstrightedge = find(sdfdone == min(sdfdone),1,'first');
    firsttrough = find(sdfdone(firstrightedge:end) > 0,1,'first')+firstrightedge-1;
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
    
    threshLocation = bincenters(firstrightedge);
    
end

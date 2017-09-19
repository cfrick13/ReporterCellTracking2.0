
function [If,testOut] = segmentCellBackgroundNEW(img,background_seg,pStruct,frames)
testOut=struct();    
%         channel = alterChanName(channel);
nucDiameter = pStruct.(background_seg).nucDiameter;
threshFactor = pStruct.(background_seg).threshFactor;
sigmaScaledToParticle = pStruct.(background_seg).sigmaScaledToParticle;
kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize
percentSmoothed = pStruct.(background_seg).percentSmoothed;

efiltsize = nucDiameter;
if mod(efiltsize,2)==0
    efiltsize = efiltsize+1;
end
% output = entropyfilt(uint16(img),true(efiltsize));
output = stdfilt(img,true(efiltsize));
[~,~,~,threshLocation] = method4(output(:));
outputlog = false(size(output));
outputlog(:) = false;
outputlog(output<(threshLocation*threshFactor))=true;
Ih = ~outputlog;
areaOfSegmentation = sum(Ih(:));
percentageOfImageSegmented = round(100*(areaOfSegmentation./(size(img,1)*size(img,2))));
if percentageOfImageSegmented < 95
    Ih = imdilate(Ih,strel('disk',floor(nucDiameter./2)));
end



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


%    
   
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

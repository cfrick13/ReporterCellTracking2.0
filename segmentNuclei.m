
function [If,testOut] = segmentNuclei(img,nucleus_seg,pStruct,frames)
    nucDiameter = pStruct.(nucleus_seg).nucDiameter;
    threshFactor = pStruct.(nucleus_seg).threshFactor;
    sigmaScaledToParticle = pStruct.(nucleus_seg).sigmaScaledToParticle;
    fnames = fieldnames(pStruct.(nucleus_seg));
    if sum(strcmp(fnames,'metthresh'))>0
        metthresh = pStruct.(nucleus_seg).metthresh;
    else
        metthresh=0.1;
    end
    fnames = fieldnames(pStruct.(nucleus_seg));
    if sum(strcmpi(fnames,'denoise'))>0
        wienerP=pStruct.(nucleus_seg).denoise;
    else
        wienerP=5;
    end
    
    percentSmoothed = pStruct.(nucleus_seg).percentSmoothed;
    testOut = struct();
           
            imgRaw = img;
            imgRawDenoised = wiener2(imgRaw,[wienerP wienerP]);


            %Based on algorithm of Fast and accurate automated cell boundary determination for fluorescence microscopy by Arce et al (2013)   
            %LOW PASS FILTER THE IMAGE (scale the gaussian filter to diameter of
            %nuclei -- diameter of nuclei is about 50 to 60))
            kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
            sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize
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

            

            %rescale image
                lcontrast = 0;
                tcontrast = 100;
                lprcntl = prctile(rawMinusLPScaled(:),lcontrast);
                prcntl = prctile(rawMinusLPScaled(:),tcontrast);
                scaleFactor = 1./(prcntl - lprcntl);
                rawMinusLPScaledContrasted = rawMinusLPScaled.*scaleFactor;
                rawMinusLPScaledContrasted = rawMinusLPScaledContrasted-(lprcntl.*scaleFactor);

            
            
            vecOG = rawMinusLPScaledContrasted(:);
            logvecpre = vecOG; logvecpre(~(logvecpre>0))=[];
            logvec = log10(logvecpre);
            vec = logvec;
           

%             [sdf1,frac1,bc1,threshLocation1] = method1(vec);
%             [sdf2,frac2,bc2,threshLocation2] = method2(vec);
%             [sdf3,frac3,bc3,threshLocation3] = method3(vec);
            [~,~,~,threshLocation] = method3(vec);
            
%             frameNum=2;
%             if frames == frameNum
%             figure(99)
%             subplot(1,3,1);plot(bc1,frac1);hold on;plot(bc1(2:end),sdf1);stem(threshLocation1,1);
%             title(num2str(frameNum));
%             subplot(1,3,2);plot(bc2,frac2);hold on;plot(bc2(2:end),sdf2);stem(threshLocation2,1);
%             subplot(1,3,3);plot(bc3,frac3);hold on;plot(bc3(2:end),sdf3);stem(threshLocation3,1);
%             end

            
            


            subtractionThreshold = threshLocation;

            if size(subtractionThreshold,1)==size(subtractionThreshold,2)
                else
                 subtractionThreshold = mean(threshLocation);
            end
            

            subtractionThresholdScaled = (10.^subtractionThreshold).*threshFactor;
            subtracted = single(rawMinusLPScaledContrasted)-subtractionThresholdScaled;
            subzero = (subtracted<0);
            Ih = ~subzero;
            Ihe = imerode(Ih,strel('disk',2));
            Ihed = imdilate(Ihe,strel('disk',2));
            Ihc = imclose(Ihed,strel('disk',2));
            Ihcf = imfill(Ihc,'holes');
            Im=Ihcf;
%             Im=Ih;


            %%%% this is the ultimate addition for watershed segmentation!!!
            see = strel('disk',1);
            Isum = Im;
            Ier = Isum;
            for i=1:round((nucDiameter/2))
                Ier = imerode(Ier,see);
                Isum = Isum+Ier;
            end
            Isum(Isum>nucDiameter) = nucDiameter;
            waterBoundary = imerode(Im,strel('disk',1));

            %BEGIN THE WATERSHET ALGORITHM
%             I = imgRawDenoised;
%             I = gaussianBlurz(rawMinusLPScaled,sigma./4,kernelgsize);
%             I = rawMinusLPScaledContrasted;

            Inew = rawMinusLPScaledContrasted;
            pxvals = Inew(waterBoundary); 
            px10 = prctile(pxvals,percentSmoothed); 
            Inew(Inew>px10)=px10;

            

            %gradmag
            hy = fspecial('sobel');
            hx = hy';
            Iy = imfilter(single(Inew), hy, 'replicate');
            Ix = imfilter(single(Inew), hx, 'replicate');
            gradmag = sqrt(Ix.^2 + Iy.^2);

            %Smoothing and identification of regional maxima (seeding watershed)
            %fgm4
            I = single(Isum);
            width = round(nucDiameter./10);
            se = strel('disk', width);
            Io = imopen(I, se);
            Ie = imerode(Io, se);
            Ieg = gaussianBlurz(Ie,round(sigma./2),round(kernelgsize./2));
            fgm = imregionalmax(Ieg);
            width = round(nucDiameter./10);
            fgm4 = imdilate(fgm,strel('disk',width));

            %bgm
            bw = Im;
            D = bwdist(bw);
            DL = watershed(D,4);
            bgm = DL == 0;
%             gradmag2 = uint16(imimposemin(gradmag, bgm | fgm4));
            gradmag2 = imimposemin(gradmag, bgm | fgm4);

            %L
            L = watershed(gradmag2,8);
            L(waterBoundary<1) = 0;
            If = L>0;

            %remove incorrect nuclei
            CellObjects = bwconncomp(If,8);

            nm = CellObjects.NumObjects;

            %determine nuclei that meet roundness criteria
            stats = regionprops(CellObjects,'Area','Perimeter');
            areavec = horzcat(stats.Area);
            perimetervec = horzcat(stats.Perimeter);
            metric = 4.*pi.*areavec./(perimetervec.^2);
            %metric = 4*pi*area/perimeter^2.
%             metthresh = 0.1;
            metriclog = (metric>metthresh);



            %determine nuclei that meet size criteria
            PX = CellObjects.PixelIdxList;
            pxl = cellfun(@length,PX,'UniformOutput',1);
            obRad = (nucDiameter./2);
            objectArea = pi.*(obRad.^2);

            pxlogS = pxl>(objectArea./8); %only keep if area is bigger than small limit
%             pxlogL = pxl<(objectArea.*8); %only keep if area is smaller than large limit
            pxlogL = pxl<(objectArea.*16); %only keep if area is smaller than large limit

            %apply logicals that remove small, large, and non-round segmented objects
            PXX = PX(~(pxlogS & pxlogL & metriclog));
            If(vertcat(PXX{:})) = 0;
            
            

            %%%%%%%%%remove segmentation if it overlaps with the border of theimage%%%%%%%%%%%%%%%%%%%
%             testImage = false(size(If));
%             dim = size(If);
%             bW = 1; %bW = borderWidth
%             testImage(1:(1+bW),1:dim(2)) =1; testImage((dim(1)-bW):dim(1),1:dim(2)) =1; testImage(1:dim(1),(1:1+bW)) =1; testImage(1:dim(1),(dim(2)-bW):dim(2)) =1;
%             borderpixels = find(testImage == 1);
%             for pidx = 1:length(PX)
%                 px = PX{pidx};
%                 testlog = ismember(px,borderpixels);
%                 if sum(testlog)>(nucDiameter/2)
%                     If(px) = 0;
%                 end
% 
%             end

            
            if frames ==7
                eo=1;
            end
            if frames==1
                testOut.img = img;
                testOut.imgRawDenoised = imgRawDenoised;
                testOut.imgLowPass = imgLowPass;
                testOut.rawMinusLP = rawMinusLP;
                testOut.rawMinusLPScaled = rawMinusLPScaled;
                testOut.Ih = Ih;
                testOut.Ihcf = Ihcf;
                testOut.Im = Im;
                testOut.Ieg = Ieg;
                testOut.fgm4 = fgm4;
                testOut.Ie = Ie;
                testOut.Inew = Inew;
                testOut.L = L;
                testOut.gradmag = gradmag;
                testOut.gradmag2 = gradmag2;
                %testOut.waterBoundary = waterBoundary;
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


function [sdfdone,fraction,bincenters,threshLocation] = method1(vec)
    lowperc = prctile(vec,1);
    highperc = prctile(vec,100);
    [numbers,bincenters] = hist(vec,lowperc:(highperc-lowperc)/500:highperc);
%     numbersone = movmean(numbers,10,2,'Endpoints','shrink');
%     numberstwo = movmean(numbersone,100,2,'Endpoints','shrink');
%     numbersone = movmean(numbers,10,2);
%     numbersone = smooth(numbers,'rlowess');
    numbersone = movmean(numbers,5,2,'Endpoints','fill')';
    movmeanmagnitude= 50;
    numberstwo = movmean(numbersone,movmeanmagnitude,1,'Endpoints','shrink');
    fraction = numberstwo./max(numberstwo);
    diffFraction = diff(fraction,1,1);
    diffFraction = diffFraction./max(diffFraction);
%     sdf = smooth(diffFraction,'lowess');
    sdf = movmean(diffFraction,5,1,'EndPoints','fill');
%     sdf = movmean(diffFraction,5,1);
%     sdfdone = movmean(sdf,50,1);
    sdfdone = movmean(sdf,5,1,'Endpoints','fill');
    sdfdone(1:movmeanmagnitude)=NaN; sdfdone(end-movmeanmagnitude:end) =NaN;
    
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

function [sdfdone,fraction,bincenters,threshLocation] = method2(vec)
    lowperc = prctile(vec,1);
    highperc = prctile(vec,100);
    [numbers,bincenters] = hist(vec,lowperc:(highperc-lowperc)/250:highperc);
%     numbersone = movmean(numbers,10,2,'Endpoints','shrink');
%     numberstwo = movmean(numbersone,100,2,'Endpoints','shrink');
%     numbersone = movmean(numbers,10,2);
%     numbersone = smooth(numbers,'rlowess');
    numbersone = movmean(numbers,5,2,'Endpoints','fill')';
    movmeanmagnitude= 25;
    numberstwo = movmean(numbersone,movmeanmagnitude,1,'Endpoints','shrink');
    fraction = numberstwo./max(numberstwo);
    diffFraction = diff(fraction,1,1);
    diffFraction = diffFraction./max(diffFraction);
%     sdf = smooth(diffFraction,'lowess');
    sdf = movmean(diffFraction,5,1,'EndPoints','fill');
%     sdf = movmean(diffFraction,5,1);
%     sdfdone = movmean(sdf,50,1);
    sdfdone = movmean(sdf,5,1,'Endpoints','fill');
    sdfdone(1:movmeanmagnitude)=NaN; sdfdone(end-movmeanmagnitude:end) =NaN;
    
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
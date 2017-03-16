
function [If,testOut] = segmentNuclei(img,nucleus_seg,pStruct,frames)
    nucDiameter = pStruct.(nucleus_seg).nucDiameter;
    threshFactor = pStruct.(nucleus_seg).threshFactor;
    sigmaScaledToParticle = pStruct.(nucleus_seg).sigmaScaledToParticle;
    wienerP=5;
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
            rawMinusLPScaledvec = reshape(rawMinusLPScaled,size(rawMinusLPScaled,1)^2,1);
            high_in = prctile(rawMinusLPScaledvec,99);
            rawMinusLPScaledContrasted = imadjust(rawMinusLPScaled./high_in,[0.1; 0.99],[0; 1]);
            vecOG = single(reshape(rawMinusLPScaledContrasted,size(rawMinusLPScaledContrasted,1)^2,1));
            logvecpre = vecOG; logvecpre(logvecpre==0)=[];
            logvec = log10(logvecpre);
            vec = logvec;
            lowperc = prctile(vec,1);
            highperc = prctile(vec,100);
            [numbers,bincenters] = hist(vec,lowperc:(highperc-lowperc)/1000:highperc);
            numbersone = medfilt1(numbers, 10); %smooths curve
            numberstwo = medfilt1(numbersone, 100); %smooths curve
            fraction = numberstwo./sum(numberstwo);
            mf = max(fraction);
                %%%%%%%%%%%%%%%%%%% Important parameters for finding minima of
                %%%%%%%%%%%%%%%%%%% histogram
                left=0.5*mf;
                slopedown=0.4*mf;
                %%%%%%%%%%%%%%%%%%%%
            leftedge = find(fraction > left,1,'first');
            insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
            threshLocation = bincenters(leftedge+insideslopedown);
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
            I = imgRawDenoised;
            I = gaussianBlurz(rawMinusLPScaled,sigma./4,kernelgsize);
            I = rawMinusLPScaledContrasted;

            %gradmag
            hy = fspecial('sobel');
            hx = hy';
            Iy = imfilter(single(I), hy, 'replicate');
            Ix = imfilter(single(I), hx, 'replicate');
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
            gradmag2 = uint16(imimposemin(gradmag, bgm | fgm4));
            gradmag2 = imimposemin(gradmag, bgm | fgm4);

            %L
            L = watershed(gradmag2,8);
            L(waterBoundary<1) = 0;
            If = L>1;

            %remove incorrect nuclei
            CellObjects = bwconncomp(If,8);

            %determine nuclei that meet roundness criteria
            stats = regionprops(CellObjects,'Area','Perimeter');
            areavec = horzcat(stats.Area);
            perimetervec = horzcat(stats.Perimeter);
            metric = 4.*pi.*areavec./(perimetervec.^2);
            %metric = 4*pi*area/perimeter^2.
            metthresh = 0.6;
            metriclog = (metric>metthresh);



            %determine nuclei that meet size criteria
            PX = CellObjects.PixelIdxList;
            pxl = cellfun(@length,PX,'UniformOutput',1);
            obRad = (nucDiameter./2);
            objectArea = pi.*(obRad.^2);

            pxlogS = pxl>(objectArea./8); %only keep if area is bigger than small limit
            pxlogL = pxl<(objectArea.*8); %only keep if area is smaller than large limit

            %apply logicals that remove small, large, and non-round segmented objects
            PXX = PX(~(pxlogS & pxlogL & metriclog));
            If(vertcat(PXX{:})) = 0;
            
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

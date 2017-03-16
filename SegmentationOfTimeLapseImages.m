function  SegmentationOfTimeLapseImages(~,datename,BACKGROUND,segInstruct)

nucleus_seg = segInstruct.nucleus;
cell_seg = segInstruct.cell;
background_seg = segInstruct.background;

%%%%
%determine path to .m file being executed
    mdir = mfilename('fullpath');
        [~,b] = regexp(mdir,'Tracking\w*/');
            if isempty(b)
                [~,b] = regexp(mdir,'Tracking\w*\');
            end
    parentdir = mdir(1:b);
    exportdir = strcat(parentdir,'Export/');

%determine path to gparent folder
    [~,b ] = regexp(parentdir,'/');
        if isempty(b)
            [~,b] = regexp(parentdir,'\');
        end
        gparentdir = parentdir(1:b(end-1));

    %assign path to experiment folder
    cd(gparentdir)
    experimentdir = strcat(gparentdir,datename);   
    cd(experimentdir)


    %assign path for flatfield corrected images and segmentation folder
    mstackName = 'flat mstack';
    mstackPath = strcat(experimentdir,'/',mstackName);
    segmentName = 'segment mstack';
    segmentPath = strcat(experimentdir,'/',segmentName);
    
%make segmentation dir of not present
    cd(experimentdir)
    dirlist = dir((segmentName));
    if isempty(dirlist)
        mkdir(strcat(segmentName));
    end




%define parameter structure and default parameter values
    pStruct = struct();
    parameterDefaults.EGFP = [106 1 15 1];
    parameterDefaults.CFP = [20 1 2 1];
    parameterDefaults.mKate = [40 1 2 1];
    parameterDefaults.Hoechst = [20 1 2 1];
    parameterDefaults.DIC = [20 1 2 1];
    parameterStrings = {'nucDiameter','threshFactor','sigmaScaledToParticle','noparametercurrently'};
    channelList = fieldnames(parameterDefaults);
    for p = 1:length(parameterStrings)
        pString = char(parameterStrings{p});
        for c = 1:length(channelList)
            cstr = char(channelList{c});
            cString = alterChanName(cstr);
            pd = parameterDefaults.(cString);
            pStruct.(cString).(pString) = pd(p); 
        end
    end
    pStruct = loadSegmentParameters(pStruct,datename,exportdir); %loads saved value of pStruct


bkarray = bkarraymaker(BACKGROUND);
bkinputs =channelregexpmaker(bkarray);

cd(mstackPath)
nucleusFileList = dir(mstackPath);
[~,~,~,d] = regexp({nucleusFileList.name},'s[0-9]+');
dlog = ~cellfun(@isempty,d,'UniformOutput',1); 
dcell = d(dlog);
SceneList = unique(cellfun(@(x) x{1},dcell,'UniformOutput',0));
    
[~,~,~,d] = regexp(SceneList,bkinputs);
bkgscenelog = cellfun(@isempty,d,'UniformOutput',1);
SceneList = SceneList(bkgscenelog);


    
        for i=1:length(SceneList)
            tic
            sceneStr = SceneList{i};
            sceneName = char(sceneStr);

            cd(mstackPath)
            nucleusFileList = dir(strcat('*',sceneName,'*',nucleus_seg,'*.mat'));
            if isempty(nucleusFileList)
                error(strcat('the channel set for "nucleus_seg" ("',nucleus_seg,'") does not exist'))
            end
          
            
            nucleusFileName = char(nucleusFileList.name);
            fileObject = matfile(nucleusFileName);
            FinalImage = fileObject.flatstack;
            disp(sceneName)
            [~,~] = segmentationNucleus(FinalImage,segmentPath,nucleus_seg,nucleusFileName,pStruct);
                

            % %load image of segmented nuclei
            % cd(segmentPath)
            % nucleusFileList = dir(strcat('*',sceneName,'*',nucleus_seg,'*.mat'));
            % nucleusFileName = char(nucleusFileList.name);
            % fileObject = matfile(nucleusFileName);
            % FinalImage = fileObject.flatstack;

            %segmentationREPORTERBKG(FinalImage,'mKate',scenename,fname,'BKGbinary_flat',pStruct);
            %segmentationRFP(FinalImage,subdirname,scenename,fname,'mKatebinary_flat');
            %segmentationMNG(FinalImage,subdirname,scenename,fname,'EGFPbinary_flat',pStruct);
            
            
            cd(mstackPath)
            backgroundFileList = dir(strcat('*',sceneName,'*',background_seg,'*.mat'));
            if isempty(backgroundFileList)
                error(strcat('the channel set for "background_seg" ("',background_seg,'") does not exist'))
            end
            backgroundFileName = char(backgroundFileList.name);
            fileObject = matfile(backgroundFileName);
            FinalImage = fileObject.flatstack;
            disp(sceneName)
            [~,~] = segmentationImageBackground(FinalImage,segmentPath,background_seg,backgroundFileName,pStruct);
            toc
            
            cd ..
        end
end


function pStruct = loadSegmentParameters(pStruct,datename,exportdir)


cd(exportdir)
filename = strcat('*',datename,'*segmentParameters*');

filelist = dir(filename);
if ~isempty(filelist)
loadname = char((filelist.name));
A = load(loadname); %load pstruct values
pStruct = A.pStruct;
else
    disp('RUN uiSegmentTimeLapseImages to set segmentation parameters')
end

    

end

function FinalImage=loadStack(FileTif)
% [a,b] = uigetfile;
% FileTif = a;
% cd (b)
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   FinalImage(:,:,i)=TifLink.read();
end
TifLink.close();
end






function [IfFinal,testOut] = segmentationNucleus(FinalImage,segmentPath,nucleus_seg,nucleusFileName,pStruct)
                           
    testOut = struct();
    frames = 1;
    img = FinalImage(:,:,frames); 
    [~,testOut] = segmentNuclei(img,nucleus_seg,pStruct,frames);
    
    IfFinal = false(size(FinalImage));
    
    
    %determine number of parallel cores to employ based on memory requirements
        dim = size(img);
        memoryrequired = size(dim(1)*dim(2)*2);%background stack + imagestack
            disp(strcat('memory required for one imgae =',num2str(round(memoryrequired./(1e6),1,'decimals')),'MegaBytes'))
        detWorkers = 4.5./(memoryrequired./1e9); %determine the number of workers you can use while staying under 6 GB
        nWorkers = floor(detWorkers);
        possibleWorkers = feature('numcores');
        if nWorkers>possibleWorkers
            nWorkers = possibleWorkers;
        elseif nWorkers <1
            nWorkers = 1;
        end
        
     % Enter parallel loop
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            poolobj = parpool(nWorkers);
        else
        nw = poolobj.NumWorkers;
            if nw ==nWorkers
            else
                delete(poolobj)
                poolobj = parpool(nWorkers);
            end
        end
    
    % start segmentation
    parfor frames = 1:size(FinalImage,3)
        img = FinalImage(:,:,frames); 
        [If,~] = segmentNuclei(img,nucleus_seg,pStruct,frames);
        IfFinal(:,:,frames)=If;
    end
                
    savethatimagestack(IfFinal,nucleusFileName,segmentPath)
    stophere=1;
end

function If = segmentationCell(FinalImage,subdirname,scenename,filename,channel,pStruct)
global nucleus_seg foldernameglobal
fig=1;

dirlist = dir(strcat(channel));
if isempty(dirlist)
mkdir(strcat(channel));
end

% parameters
left = 0.004;
slopedown = 0.003;

dimdiff = 2048./size(FinalImage(:,:,1),1);

zerostrel = 2;
firststrel = round(30./dimdiff);
sigmafirst = firststrel.*3;
kernelgsizefirst = firststrel.*6;
fracsmoothing = 0.5.*dimdiff;

dirlist = dir(nucleus_seg);
if isempty(dirlist)
    dirlist = dir(foldernameglobal);
    foldername = foldernameglobal;
else
    foldername = nucleus_seg;
end


tsn = determineTimeFrame(foldername);


% start
for frames = 1:size(FinalImage,3)
img = FinalImage(:,:,frames); 
imgorig = img;

img = wiener2(img,[5 5]);
se =strel('disk',zerostrel);
Ie = imerode(imgorig,se);
Iobr = imreconstruct(Ie,img);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
IobrcbrF = imcomplement(Iobrcbr);
gaus = double(IobrcbrF);

se =strel('disk',firststrel);
Ie = imerode(gaus,se);
Iobr = imreconstruct(Ie,gaus);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
gaus = double(Iobrcbr);


sigma = sigmafirst;
kernelgsize = kernelgsizefirst;
gaustwo = gaussianBlurz(double(gaus),sigma,kernelgsize);

sub = double(gaus) -double(gaustwo);%%%%%%% key step!
b = find(sub == min(min(sub)),1,'first');
rattio = gaustwo(b)./gaus(b);
gaustwocorr = gaustwo./rattio;
sub_scale_corr = double(gaus) - double(gaustwocorr);




subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
[numbers,bincenters] = hist(double(vec),0:fracsmoothing:10000);



numbers = medfilt1(numbers, 10); %smooths curve
fraction = numbers./sum(numbers);

mf = max(fraction);

%%%%%%%%%%%%%%%%%%%%
left=0.3*mf;
slopedown=0.2*mf;
%%%%%%%%%%%%%%%%%%%%%

leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
insideslopeup = find(fraction(leftedge+insideslopedown:end) >0.0012,1,'first');
trough = min(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup));
troughindex = find(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup) == trough);
troughindexrounded = round(median(troughindex));



threshlocation = bincenters(leftedge+insideslopedown);
% 
% figure(22)
%     bar(bincenters,fraction);hold on
%     xlim([-500 1000])
%     ylim([0 0.1])
% stem ([threshlocation threshlocation],[0 1]);hold off
% drawnow

subtractionthreshold = threshlocation;

if size(subtractionthreshold,1)==size(subtractionthreshold,2)
else
     subtractionthreshold = mean(threshlocation);
end

% subtractionthreshold = graythresh(subtractionref);

subtracted = sub_scale_corr-subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
a = find(Ie>0);
submax = zeros(size(Ie));
Ie(a)=50;





Ih = Ie>0;
Igclose = imclose(Ih,strel('disk',round(30./dimdiff)));
Igclosemax = imclose(Ih,strel('disk',round(80./dimdiff)));
Igcopenmax = imopen(Igclosemax,strel('disk',round(10./dimdiff)));
Igcopen = imopen(Igclose,strel('disk',2));
Igcofill = imfill(Igcopen,'holes');
Igcfopen = bwareaopen(Igcofill,round(5000./(dimdiff.^2)));
Igcfopendil = imerode(Igcfopen,strel('disk',round(5)));

% finalsigma=20;
% finalkernelgsize=40;
gaus = gaussianBlurz(IobrcbrF,round(sigmafirst./dimdiff),(kernelgsizefirst./dimdiff));

% sigma=40;
% kernelgsize=80;
% gaus = gaussianBlurz(gaus,sigmafirst.*2,kernelgsizefirst.*2);

imgt = -double(gaus);
% imgt(~(Igcfopen>0)) = -Inf;
imgt(~(Igcopenmax>0)) = -Inf;

L=watershed(imgt);

L(Igcfopendil<1) = 0;
% imagesc(L)
% colormap parula
If = L>1;

stophere=1;
% time = settimecharacter(frames);
time = tsn{frames};
tim = time(2:end);
savethatimagestack(scenename,time,If.*255,frames,filename,channel)
end


end

function [IfFinal,testOut] = segmentationImageBackground(FinalImage,segmentPath,background_seg,backgroundFileName,pStruct)
% parameters
testOut = struct();
nucDiameter = pStruct.(background_seg).nucDiameter;
threshFactor = pStruct.(background_seg).threshFactor;
sigmaScaledToParticle = pStruct.(background_seg).sigmaScaledToParticle;
kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize



%initial segmentation to determine how much of image is covered by cells
img = FinalImage(:,:,1); 
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
end
% disp(percentageOfImageSegmented);

    
IfFinal = false(size(FinalImage));
parfor frames = 1:size(FinalImage,3)

    img = FinalImage(:,:,frames); 
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
    IfFinal(:,:,frames)=If;
end

savethatimagestack(IfFinal,backgroundFileName,segmentPath)
end


function If = segmentationREPORTERBKG(FinalImage,channel,scenename,filename,segchannel,pStruct)
global   nucleus_seg foldernameglobal
mkdir(strcat(segchannel));

foldername = foldernameglobal;
% foldername = foldernameglobal;
tsn = determineTimeFrame(foldername);

% parameters
nucDiameter = pStruct.(channel).nucDiameter;
threshFactor = pStruct.(channel).threshFactor;
sigmaScaledToParticle = pStruct.(channel).sigmaScaledToParticle;
kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize


for frames=1:size(FinalImage,3)
img = FinalImage(:,:,frames);
If = logical(img);
se = strel('disk',20);
Ifd = imdilate(If,se);
If = ~Ifd;


    
    time = tsn{frames};
    tim = time(2:end);
%     IfFinal(:,:,frames)=If;
    savethatimagestack(scenename,time,If.*255,frames,filename,segchannel)
end


stophere=1;
end

function savethatimagestack(IfFinal,filename,segmentPath)
olddir  = pwd;
cd(segmentPath)
save(filename,'IfFinal','-v7.3');
cd (olddir)


end

function bw = gaussianBlurz(im,sigma,kernelgsize,varargin)

filtersize = [kernelgsize kernelgsize];
kernelg = fspecial('gaussian',filtersize,sigma);

%% image filtering
gFrame = imfilter(im,kernelg,'repl');

if ~isempty(varargin)
    bw=gFrame.*uint16(varargin{1}>0);
else
    bw=gFrame;
end
end


function bw = logMasked(im,ksize,varargin)
% Discrete Laplacian
kernel = chooseKernel(ksize);
%% image filtering
lapFrame = imfilter(im,kernel,'repl');
if ~isempty(varargin)
    bw=lapFrame.*uint16(varargin{1}>0);
else
    bw=lapFrame;
end
end

function kernel = chooseKernel(ksize)
if ksize ==5
kernel = [-4 -1  0 -1 -4;...
     -1  2  3  2 -1;...
      0  3  4  3  0;...
     -1  2  3  2 -1;...
     -4 -1  0 -1 -4];


% % % -4 -1  0 -1 -4
% % % -1  2  3  2 -1
% % % 0  3  4  3  0
% % % -1  2  3  2 -1
% % % -4 -1  0 -1 -4

elseif ksize == 7
kernel =[-10 -5 -2 -1 -2 -5 -10;... 
    -5  0  3  4  3  0  -5;... 
    -2  3  6  7  6  3  -2;... 
    -1  4  7  8  7  4  -1;... 
    -2  3  6  7  6  3  -2;... 
    -5  0  3  4  3  0  -5;... 
    -10 -5 -2 -1 -2 -5 -10];... 
    
% % % -10 -5 -2 -1 -2 -5 -10 
% % % -5  0  3  4  3  0  -5 
% % % -2  3  6  7  6  3  -2 
% % % -1  4  7  8  7  4  -1 
% % % -2  3  6  7  6  3  -2 
% % % -5  0  3  4  3  0  -5
% % % -10 -5 -2 -1 -2 -5 -10
end
end


function  LoGstack = LaplacianOfGaussianStack(imgstack,dims,ksize)
        LoGstack = zeros(dims(1),dims(2),dims(3));
        for i = 1:size(imgstack,3)
        LoGstack(:,:,i) = logMasked(imgstack(:,:,i),ksize);
        end
        
end

function tsn = determineTimeFrame(foldername)
cd(foldername)
fflist = dir(strcat('*.tif'));
ffnames = {fflist.name};
[a,b,c,d] = regexp(ffnames,'_t[0-9]+');
[a,b,c,d] = regexp(ffnames,'_t[0-9]++');
tnames = cellfun(@(x) x{1},d,'UniformOutput',0);
[a,b,c,d] = regexp(tnames,'t[0-9]++');%added extra step because sometimes the time frame parsin was mistaken
tnames = cellfun(@(x) x{1},d,'UniformOutput',0);
tsn = sort(tnames);
cd ..
end

function channelinputs =channelregexpmaker(channelstoinput)
    channelinputs = '(';
    for i=1:length(channelstoinput) % creates a string of from '(c1|c2|c3|c4)' for regexp functions
        if i ==1
        channelinputs = strcat(channelinputs,channelstoinput{i});
        elseif i < length(channelstoinput)
            channelinputs = strcat(channelinputs,'|',channelstoinput{i});
        else
            channelinputs = strcat(channelinputs,'|',channelstoinput{i},')');
        end
    end
end
function bkarray = bkarraymaker(BACKGROUND)
    for i = 1:length(BACKGROUND)
        bkstr = num2str(BACKGROUND(i)); 
        if length(bkstr)>1
            bkarray{i} = strcat('s',bkstr);
        else
            bkarray{i} = strcat('s0',bkstr); 
        end
    end
end

function chanstruct = alterChanName(chan)
[a,~] = regexp(chan,'(_|\W|\s|[0-9])'); %remove underscore, dashes, or whitespace
chanstruct = chan;
chanstruct(a) = [];

% if isempty(chanstruct)
%     chanstruct = 'CFP';
% end
end

function Seg_Flat_and_Sub_TimeLapseImages(~,datename,channelstoinput,BACKGROUND)
% datename = '2015_12_15 smad3g smFISH';
E={datename}; %beginning of file name for saving

%determine path to .m file being executed
mdir = mfilename('fullpath');
[~,b] = regexp(mdir,'Tracking\w*/');
if isempty(b)
    [~,b] = regexp(mdir,'Tracking\w*\');
end
parentdir = mdir(1:b);

%determine path to gparent folder
[~,b ] = regexp(parentdir,'/');
if isempty(b)
    [~,b] = regexp(parentdir,'\');
end
gparentdir = parentdir(1:b(end-1));

%assign path to experiment folder
cd(gparentdir)
expPath = strcat(gparentdir,datename);
cd(expPath)

%determine date and experiment number
[a,b] = regexp(datename,'exp[0-9]');
expnum = datename(a:b);
dateofexp = datename(1:10);

%assign path to .mat image stack exports
mstackdir = 'mstack images';
flatfielddirname = 'flat mstack';
mstackPath = strcat(expPath,'/',mstackdir);

%assign path to folder for metadata export files
metadir = strcat(parentdir,'Export/');
cd(metadir)



%find associated extracted metadata
% filelist  = dir(strcat('*',dateofexp,'*',expnum,'*metaData.mat'));
filelist  = dir(strcat('*',datename,'*metaData.mat'));
metadatafile = char(filelist.name);
A = load(metadatafile);
dim = A.dimensions;
timeCount = A.timeCount;
channelinputs =channelregexpmaker(channelstoinput);

% load excel-extracted data
datequery = strcat(datename,'*DoseAndScene*');
cd(metadir)
filelist = dir(datequery);
if isempty(filelist)
    error(strcat('need to run ExtractMetadata for-',FileName));
    %dosestruct = makeDoseStruct; %run function to make doseStruct
else
    dosestructstruct = load(char(filelist.name));
    dosestruct = dosestructstruct.dosestruct;
    segInstruct = dosestructstruct.segInstruct;
end
conditionsListArray = {dosestruct.conditions};
conditionslog = ~cellfun(@isempty,conditionsListArray,'UniformOutput',1);
conditionsListArrayCells = conditionsListArray(conditionslog);
conditionsTruncListArray = cellfun(@(x) x{1},conditionsListArrayCells,'UniformOutput',0);
conditionsListArray(~conditionslog)={'empty'};
conditionsListArray(conditionslog) = conditionsTruncListArray;

for i = 1:length(conditionsListArray)
    cond = conditionsListArray{i};
    numsinthere = regexp(cond,'[0-9]');
    if ~isempty(numsinthere)
        wordstr = num2word(str2num(cond(numsinthere)));
        cond(numsinthere)=[];
        conditionsListArray{i} = [cond wordstr];
    else
    end
end
conditionsList = unique(conditionsListArray(conditionslog));

cd(mstackPath)
dirlist = dir('*.mat');
[~,~,~,channelsListed] = regexp([dirlist.name],channelinputs);
channelList = unique(channelsListed);

cd(expPath)
dirlist = dir(strcat('*',flatfielddirname,'*'));
flatPath = strcat(expPath,'/',flatfielddirname);
if isempty(dirlist)
    mkdir(flatPath);
end




sceneVector = A.sceneVector;
msv = max(sceneVector);
backarray = cell(1,length(BACKGROUND));
cycle=1;
for g = BACKGROUND
    if msv>99
        backstr = 's000';
    else
        backstr = 's00';
    end
    backnumstr = num2str(g);
    backstr(end-(length(backnumstr)-1):end) = backnumstr;
    backarray{cycle} = backstr;
    cycle=cycle+1;
end
backinputs = channelregexpmaker(backarray);
cd(mstackPath)
dirlist = dir('*.mat');
dirlistarray = {dirlist.name};
[~,~,~,SceneFileNames] = regexp(dirlistarray,backinputs);
sfnlog = ~(cellfun(@isempty,SceneFileNames,'UniformOutput',1));
%logical for only background


[~,~,~,SceneFileNames] = regexp(dirlistarray,'s[0-9]++');
sfnarray= cellfun(@(x) x{1},SceneFileNames,'UniformOutput',0);
sceneListArray = unique(sfnarray);





%create the median stack of media only images for all timepoints for
%each channel -- these will be used for flatfield correction in next
%step
% channelList = channelList(2);
for condition = conditionsList
    for channel = channelList %cycle through one channel at a time
        cond = char(condition);
        chan = char(channel);
        [~,~,~,channelsListed] = regexp(dirlistarray,chan);
        cllog = ~(cellfun(@isempty,channelsListed,'UniformOutput',1));
        %         bkgForChannel = dirlistarray(sfnlog & cllog);
        bkgForChannel = dirlistarray(cllog);
        
        [~,~,~,condList] = regexp(conditionsListArray,cond);
        condlog = ~(cellfun(@isempty,condList,'UniformOutput',1));
        %         bkgForChannel = dirlistarray(sfnlog & cllog);
        bkgForhannelForCond = bkgForChannel(condlog);
        
        clear A dirlist
        
        
        bki = zeros(dim(1),dim(2),timeCount,length(bkgForhannelForCond),'uint16');
        for j=1:length(bkgForhannelForCond)
            bkname = char(bkgForhannelForCond{j});
            imgs = load(bkname);
            imgsz = imgs.img_stack_test;
            disp(j)
            for ii = 1:size(imgsz,3)
%             for ii = 1
                if strcmpi(chan,'DIC')
%                     imgsz(:,:,ii) = gaussianBlurz(imgsz(:,:,ii),20,50);
                else
%                     imgsz(:,:,ii) = segimg_reconstruct(imgsz(:,:,ii));
%                     imgsz(:,:,ii) = segimg_reconstruct(imgsz(:,:,ii));
                end
            end

            bki(:,:,:,j) = imgsz;
        end
        clear imgs
        
        bki = single(bki);
        bki(bki==0) = NaN;
%         bkflat = nanmedian(bki,4);
        bkflat = prctile(bki,25,4);       
        clear bki
        
        bkall = prctile(bkflat,25,3);
%         mask = false(size(bkall));
%         mask(isnan(bkall)) = true;
%         mask = imdilate(mask,strel('disk',10));
%         bkfill = regionfill(bkall,mask);
%         bkfill(bkfill>(prctile(bkall(:),50))*1.1) = prctile(bkall(:),50)*1.1;
%         
        for ii = 1:size(bkflat,3)
            bkflat(:,:,ii) = bkall;
        end
        bkflat = single(bkflat);
        
        rsbflat = reshape(bkflat,dim(1).*dim(2),timeCount);
        pp = prctile(rsbflat,99.99,1);
        ppp = prctile(rsbflat,99.999,1);
        clear rsbflat
        p = mean([pp' ppp'],2);
        normBkgimgMat = zeros(dim(1),dim(2),timeCount,'single');
        for k = 1:timeCount
            normBkgimgMat(:,:,k) = bkflat(:,:,k)./p(k);
        end
        
        chanstruct = alterChanName(chan);
        
        bkgimgstruct.(cond).(chanstruct)=normBkgimgMat;
    end
end
clear normBkgimgMat






memoryrequired = dim(1)*dim(2)*timeCount*(length(bkgForChannel)+1)*2;%background stack + imagestack

disp(strcat('memory required for one set =',num2str(round(memoryrequired./(1e6),1,'decimals')),'MegaBytes'))

detWorkers = 7./(memoryrequired./1e9); %determine the number of workers you can use while staying under 6 GB
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

cd(mstackPath)

scenesCorr = sceneListArray([1:15 41:45]);
for channel = channelList %cycle through one channel at a time
    filelist = dir(strcat('*',char(channel),'*.mat'));
    cfile = {filelist.name};
    
    tic
    for j = 1:length(cfile) %load image_stack for each scene // and determine scene and channel
        %     parfor j = 1:length(cfile) %load image_stack for each scene // and determine scene and channel
        filename = char(cfile(j));
        imgstruct = load(filename);
        imgstack = single(imgstruct.img_stack_test);
        
        [a,b] = regexp(filename,'s[0-9]++'); %determine scene
        scene = filename(a:b);
        sceneid = find(strcmpi(scene,sceneListArray));
        
        cond = conditionsListArray{sceneid};
        if strcmpi(cond,'empty')
            cond = conditionsList{1};
        end
        
        
        [e,f] = regexp(filename,channelinputs); %determine channel
        chan = filename(e:f);
        chanstruct = alterChanName(chan);
        normBkgimgMat = bkgimgstruct.(cond).(chanstruct);
        
        if sum(strcmpi(scenesCorr,scene))>0
            flatstack = imgstack./normBkgimgMat; %flatten the experimental image with the background image
        else
            flatstack = imgstack;
        end
        %         flatstack = imgstack./normBkgimgMat; %flatten the experimental image with the background image
        imgstack = [];
        normBkgimgMat = [];
        savename = strcat(E{1},'_flat_',scene,'_','_',chan,'.mat');
        SAVdir = strcat(expPath,'/',flatfielddirname,'/');
        savethatimage(savename,SAVdir,mstackPath,flatstack,j);
        %                 stophere=1;
    end
    toc
    stopehre=1;
end
stophere=1;
end



function savethatimage(savename,SAVdir,mstackPath,flatstack,j)
disp(strcat(savename,'...',num2str(j)));
cd (SAVdir);
%     imwrite(flatstack,char(savename),'tiff');
%     save(flatstack,char(savename),'tiff');
%     save(savename,'flatstack','-v7.3');
save(savename,'flatstack','-v6');
cd (mstackPath);
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

function chanstruct = alterChanName(chan)
[a,~] = regexp(chan,'(_|\W|\s|[0-9])'); %remove underscore, dashes, or whitespace
chanstruct = chan;
chanstruct(a) = [];

if isempty(chanstruct)
    chanstruct = 'CFP';
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

function wordstr = num2word(num)
nDigits = dec2base(num,10) - '0';
% wordArray = [' ';length(nDigits);' '];
% wordArray = cellstr(wordArray);
wordArray = cell(1,length(nDigits));
%wordArray = cell(1,20); tried this but no..
%these variables are switches, one the of statement containing them is ran
%they change values meaning that the switch has been turned off
a = 0;
b = 0;
c = 0;
d = 0;
e = 0;
f = 0;
for i=1:(length(nDigits))
    %millions
    if length(nDigits) >= 7 && i == 1
        switch nDigits(i)
            case 0
                wordArray{i} = ' ';
            case 1
                wordArray{i} = 'One Million';
            case 2
                wordArray{i} = 'Two Million';
            case 3
                wordArray{i} = 'Three Million';
            case 4
                wordArray{i} = 'Four Million';
            case 5
                wordArray{i} = 'Five Million';
            case 6
                wordArray{i} = 'Six Million';
            case 7
                wordArray{i} = 'Seven Million';
            case 8
                wordArray{i} = 'Eight Million';
            case 9
                wordArray{i} = 'Nine Million';
        end
        %hundred thousands
    elseif length(nDigits) >= 6 && a == 0
        a = 1;
        switch nDigits(i)
            case 0
                wordArray{i} = ' ';
            case 1
                wordArray{i} = 'One Hundred and';
            case 2
                wordArray{i} = 'Two Hundred and';
            case 3
                wordArray{i} = 'Three Hundred and';
            case 4
                wordArray{i} = 'Four Hundred and';
            case 5
                wordArray{i} = 'Five Hundred and';
            case 6
                wordArray{i} = 'Six Hundred and';
            case 7
                wordArray{i} = 'Seven Hundred and';
            case 8
                wordArray{i} = 'Eight Hundred and';
            case 9
                wordArray{i} = 'Nine Hundred and';
        end
        % Ten thousands
    elseif length(nDigits) >= 5 && b ==0
        b = 1;
        switch nDigits(i)
            case 0
                wordArray{i} = ' ';
            case 1
                % if wordArray{i+1} > 0 (teens); the c variable is so that
                % we can leave the next cell blank ''
                switch nDigits(i+1)
                    case 0
                        wordArray{i} = 'Ten Thousand';
                        c = 2;
                    case 1
                        wordArray{i} = 'Eleven Thousand';
                        c = 2;
                    case 2
                        wordArray{i} = 'Twelve Thousand';
                        c = 2;
                    case 3
                        wordArray{i} = 'Thirteen Thousand';
                        c = 2;
                    case 4
                        wordArray{i} = 'Fourteen Thousand';
                        c = 2;
                    case 5
                        wordArray{i} = 'Fifteen Thousand';
                        c = 2;
                    case 6
                        wordArray{i} = 'Sixteen Thousand';
                        c = 2;
                    case 7
                        wordArray{i} = 'Seventeen Thousand';
                        c = 2;
                    case 8
                        wordArray{i} = 'Eighteen Thousand';
                        c = 2;
                    case 9
                        wordArray{i} = 'Nineteen Thousand';
                        c = 2;
                end
            case 2
                wordArray{i} = 'Twenty ';
            case 3
                wordArray{i} = 'Thrirty ';
            case 4
                wordArray{i} = 'Fourty ';
            case 5
                wordArray{i} = 'Fivty ';
            case 6
                wordArray{i} = 'Sixty ';
            case 7
                wordArray{i} = 'Seventy ';
            case 8
                wordArray{i} = 'Eighty ';
            case 9
                wordArray{i} = 'Ninety ';
        end
        % Ten thousands
    elseif length(nDigits) >= 4 && c ==0
        c = 1;
        switch nDigits(i)
            case 0
                wordArray{i} = ' ';
            case 1
                wordArray{i} = 'One';
            case 2
                wordArray{i} = 'Two Thousand';
            case 3
                wordArray{i} = 'Three Thousand';
            case 4
                wordArray{i} = 'Four Thousand';
            case 5
                wordArray{i} = 'Five Thousand';
            case 6
                wordArray{i} = 'Six Thousand';
            case 7
                wordArray{i} = 'Seven Thousand';
            case 8
                wordArray{i} = 'Eight Thousand';
            case 9
                wordArray{i} = 'Nine Thousand';
        end
        % Ten thousands for teens
    elseif length(nDigits) >= 4 && c ==2
        c = 1;
        wordArray{i} = '';
        % Hundreds
    elseif length(nDigits) >= 3 && d ==0
        d = 1;
        switch nDigits(i)
            case 0
                wordArray{i} = ' ';
            case 1
                wordArray{i} = 'One Hundred';
            case 2
                wordArray{i} = 'Two Hundred';
            case 3
                wordArray{i} = 'Three Hundred';
            case 4
                wordArray{i} = 'Four Hundred';
            case 5
                wordArray{i} = 'Five Hundred';
            case 6
                wordArray{i} = 'Six Hundred';
            case 7
                wordArray{i} = 'Seven Hundred';
            case 8
                wordArray{i} = 'Eight Hundred';
            case 9
                wordArray{i} = 'Nine Hundred';
        end
    elseif length(nDigits) >= 2 && e ==0
        e = 1;
        nDigits(i+1)
        switch nDigits(i)
            case 1
                % if wordArray{i+1} > 0 (teens); the f variable is so that
                % we can leave the next cell blank ''
                
                switch nDigits(i+1)
                    case 0
                        wordArray{i} = 'Ten ';
                        f = 2;
                    case 1
                        wordArray{i} = 'Eleven ';
                        f = 2;
                    case 2
                        wordArray{i} = 'Twelve ';
                        f = 2;
                    case 3
                        wordArray{i} = 'Thirteen ';
                        f = 2;
                    case 4
                        wordArray{i} = 'Fourteen ';
                        f = 2;
                    case 5
                        wordArray{i} = 'Fifteen ';
                        f = 2;
                    case 6
                        wordArray{i} = 'Sixteen ';
                        f = 2;
                    case 7
                        wordArray{i} = 'Seventeen ';
                        f = 2;
                    case 8
                        wordArray{i} = 'Eighteen ';
                        f = 2;
                    case 9
                        wordArray{i} = 'Nineteen ';
                        f = 2;
                end
            case 2
                wordArray{i} = 'Twenty ';
            case 3
                wordArray{i} = 'Thrirty ';
            case 4
                wordArray{i} = 'Fourty ';
            case 5
                wordArray{i} = 'Fivty ';
            case 6
                wordArray{i} = 'Sixty ';
            case 7
                wordArray{i} = 'Seventy ';
            case 8
                wordArray{i} = 'Eighty ';
            case 9
                wordArray{i} = 'Ninety ';
        end
        
    elseif length(nDigits) >= 1 && f ==0
        f = 1;
        switch nDigits(i)
            case 0
                wordArray{i} = 'zero';
            case 1
                wordArray{i} = 'one';
            case 2
                wordArray{i} = 'two';
            case 3
                wordArray{i} = 'three';
            case 4
                wordArray{i} = 'four';
            case 5
                wordArray{i} = 'five';
            case 6
                wordArray{i} = 'six';
            case 7
                wordArray{i} = 'seven';
            case 8
                wordArray{i} = 'eight';
            case 9
                wordArray{i} = 'nine';
        end
        % double digits for teens
    elseif length(nDigits) >= 4 && f ==2
        f = 1;
        wordArray{i} = '';
    end
end
ilog = cellfun(@isempty,wordArray,'UniformOutput',1);
wordArrayz = horzcat(wordArray{~ilog});
blanklog = regexp(wordArrayz,'\s');
wordstr = wordArrayz;
wordstr(blanklog)=[];
end

function [flatten] = segimg_reconstruct(img)
%         channel = alterChanName(channel);
%%
% profile on
nucDiameter =40;
threshFactor = 1;
sigmaScaledToParticle = 1;

wienerP=15;
kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize

imgRaw = img;
imgRaw = wiener2(imgRaw,[wienerP wienerP]);
imgRaw(imgRaw>prctile(imgRaw(:),95)) = prctile(imgRaw(:),95);
imgRaw(imgRaw<prctile(imgRaw(:),20)) = prctile(imgRaw(:),20);

efiltsize = round(nucDiameter./10);
if mod(efiltsize,2)==0
    efiltsize = efiltsize+1;
end
% output = entropyfilt(uint16(img),true(efiltsize));
output = stdfilt(imgRaw,true(efiltsize));
output(output>prctile(output(:),90)) = prctile(output(:),90);
[~,~,~,threshLocation] = method4(output(:));
outputlog = false(size(output));
outputlog(:) = false;
outputlog(output<(threshLocation*threshFactor))=true;
Ih = imclose(imdilate(~outputlog,strel('disk',10)),strel('disk',5));

flatten = img;
flatten(Ih) = NaN;
% flatten = regionfill(img,Ih);
% flatten(flatten>(prctile(img(~Ih),50)).*1.2) = prctile(img(~Ih),50).*1.2;
% profile off; profile report
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

threshLocation = bincenters(firstrightedge);
threshLocation = bincenters(firsttrough);

end

function ExtractMetadataAndImages(B)
% reader = bfGetReader('/Volumes/Seagate Backup Plus Drive/FrickData/2017_01_27 plate/exp1 timelapse of all snail clones in NGSmad3 cells.czi');

%set directory to location of code being used (generally external harddrive
mdir = mfilename('fullpath');
    [~,b] = regexp(mdir,'Tracking\w*/');
        if isempty(b)
            [~,b] = regexp(mdir,'Tracking\w*\');
        end
parentdir = mdir(1:b);
exportdir = strcat(parentdir,'Export/');
readerpath = strcat(parentdir,'bfmatlab/');
    addpath(readerpath);
cd(parentdir)

[~,b ] = regexp(parentdir,'/');
    if isempty(b)
        [~,b] = regexp(parentdir,'\');
    end
    gparentdir = parentdir(1:b(end-1));

cd(gparentdir)
expdir = strcat(gparentdir,'/',B);
cd(expdir)
PathName = pwd;

filelist = dir('*plate exp*.czi');
FileName = char(filelist.name);

cd(PathName)
disp(PathName)
export_stack_dir = 'mstack images';
dirlist = dir(strcat('*',export_stack_dir,'*'));
if isempty(dirlist)
    mkdir(export_stack_dir) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct a Bio-Formats reader decorated with the Memoizer wrapper
maxjava = java.lang.Runtime.getRuntime.maxMemory;
freejava = java.lang.Runtime.getRuntime.freeMemory;
percentJavaMemFree = round(freejava./maxjava,1,'decimals');
disp(strcat('% java mem free=',num2str(percentJavaMemFree)))


reader = loci.formats.Memoizer(bfGetReader(), 0);
% Initialize the reader with an input file to cache the reader
reader.setId(FileName);
omeMeta = reader.getMetadataStore();



maxjava = java.lang.Runtime.getRuntime.maxMemory;
freejava = java.lang.Runtime.getRuntime.freeMemory;
percentJavaMemFree = round(freejava./maxjava,1,'decimals');
disp(strcat('% java mem free=',num2str(percentJavaMemFree)))
% data = bfopen(FileName);
% seriesCount = size(data, 1); %number of scenes
[a,~] = regexp(FileName,'.czi');
expName = FileName(1:a-1);
% reader = bfGetReader('/Volumes/Seagate Backup Plus Drive/FrickData/2017_01_25 plate/exp1 timelapse of caga and pctgf in cmvNGSmad3ex1 cells.czi');
% omeMeta = reader.getMetadataStore();

% omeMeta = data{1,4};

%determine dimensions of images
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
dimensions = [stackSizeX stackSizeY stackSizeZ]; %[x y z];

timeCount = omeMeta.getPixelsSizeT(0).getValue();
planeCount = omeMeta.getPlaneCount(0); % a plane is number of images for a certain scene (each channel at each timepoint)
channelCount = omeMeta.getChannelCount(0); %number of channels
imageCount = omeMeta.getImageCount(); %number of scenes in experiment
% sceneCount = size(data,1);



% change all channel names of .czi to EGFP or CFP or mKate, etc
channelNames = cell(1,channelCount);
for cycle = 1:channelCount
channelNames{cycle} = omeMeta.getChannelName(0,cycle-1).char;
end

channelNameSwapArray = cell(1,channelCount);
for cycle = 1:channelCount
    
    cname = channelNames{cycle};
    aname = regexp(cname,'DIC');
    if isempty(aname)
        emwave = double(omeMeta.getChannelEmissionWavelength(0,cycle-1).value);
        exwave = double(omeMeta.getChannelExcitationWavelength(0,cycle-1).value);
    else
        emwave = 1;
    end
    if ~isempty(aname)
        cswap = 'DIC';
    elseif emwave>505 && emwave<560
        cswap = 'EGFP';
    elseif emwave>600 && emwave<650
        cswap = 'mKate';
    elseif (emwave<480 && emwave>420) && exwave>400
        cswap = 'CFP';
    elseif exwave<370
        cswap = 'Hoechst';
    else
        a = regexp(cname,'\s');
        cswap = cname;
        cswap(a) = [];
    end
channelNameSwapArray{cycle} = cswap;
end
disp(strcat('Channel Names changed from--`',channelNames','` to `',channelNameSwapArray','`'))



% Close reader
reader.close()
clear reader omeMeta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save all images as mStack files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine amount of memory required to load one stack
memoryrequired = stackSizeX*stackSizeY*planeCount;
totalmemorysize = memoryrequired*imageCount;
disp(strcat('memory required for one set =',num2str(round(memoryrequired./(1e6),1,'decimals')),'MegaBytes'))
disp(strcat('memory required for full set =',num2str(round(totalmemorysize./(1e6),1,'decimals')),'MegaBytes'))
detWorkers = 7./(memoryrequired./1e9); %determine the number of workers you can use while staying under 6 GB
nWorkers = floor(detWorkers);
possibleWorkers = feature('numcores');
if nWorkers>possibleWorkers
    nWorkers = possibleWorkers;
elseif nWorkers <1
    nWorkers = 1;
end

% Initialize parallel loop
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



%split scenes across workers as evenly as possible
sceneVector = 1:imageCount;
fractions = length(sceneVector)./nWorkers;
firstlengths = ceil(fractions);
    if firstlengths*nWorkers > imageCount
        firstlengths = firstlengths-1;
    end
sceneArray = cell(1,nWorkers);
    for nw = 1:nWorkers
        if nw == nWorkers
            sceneArray{nw} = sceneVector(((nw-1)*firstlengths)+1:end);
        else
            sceneArray{nw} = sceneVector(((nw-1)*firstlengths)+1:firstlengths*nw);
        end
    end
cd(PathName)
tic

%run parallel image extraction
disp('...running extraction...')
disp(strcat('total images=',num2str(planeCount),'--',num2str(stackSizeX),'x',num2str(stackSizeY),'pixels'));
% for iNW = 1 : nWorkers
parfor iNW = 1 : nWorkers
    
    % Initialize logging at INFO level
    bfInitLogging('INFO');
    % Initialize a new reader per worker as Bio-Formats is not thread safe
    r2 = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
    % Initialization should use the memo file cached before entering the
    % parallel loop
    r2.setId(FileName);
    %do work
    %imageCount is number of scenes/series
    sceneArrayVec = sceneArray{iNW};
    for iSeries = sceneArrayVec
        r2.setSeries(iSeries - 1);
        disp(iSeries)
        for iC = 1:channelCount

            if max(sceneVector)>99
                sceneNumber = 's000';
            else
                sceneNumber = 's00';
            end
            scenestr = num2str(iSeries);
            sceneNumber(end-(length(scenestr)-1):end) = scenestr;
            channelName = channelNameSwapArray{iC};

            img_stack_test = zeros(stackSizeX,stackSizeY,timeCount,'uint16');
            for iT = 1:timeCount
                iPlane = r2.getIndex(0, iC -1, iT - 1) + 1;
                img_stack_test(:,:,iT) = bfGetPlane(r2, iPlane);
            end

            savename = strcat(expName,'_',channelName,'_',sceneNumber,'.mat');

%             tic
%             savefilefunc(savename,img_stack_test,'img_stack_test',export_stack_dir,PathName)
%             toc
%             tic
%             savename = strcat(expName,'_',channelName,'_',sceneNumber,'fast.mat');
            savefastfilefunc(savename,img_stack_test,'img_stack_test',export_stack_dir,PathName)
%             toc
        end
    end
    
    % Close the reader
    r2.close()
end
% cd(PathName)
disp('...extractionFinished')
toc


reader = bfGetReader(FileName);
omeMeta = reader.getMetadataStore();

t=[];
thet=[];
thec=[];
timeVec=[];

planeCountMat = 0:planeCount-1;
for i = 0:1:imageCount-1
    for pp = 1:timeCount
        p = planeCountMat(pp*channelCount);
    
        deltaT_Hash = omeMeta.getPlaneDeltaT(i,p); %Time getPlaneDeltaT(int imageIndex,int planeIndex)
        deltaT = double(deltaT_Hash.value);
        theT_Hash = omeMeta.getPlaneTheT(i,p); %getPlaneTheT(int imageIndex, int planeIndex)
        theT = double(theT_Hash.getValue);
        theC_Hash = omeMeta.getPlaneTheC(i,p); %getPlaneTheT(int imageIndex, int planeIndex)
        theC = double(theC_Hash.getValue);
        
        deltaTarray(pp) = deltaT;
        thet(pp)=theT+1;
        thec(pp) = theC+1;
        if theT==0
            t(theT+1) = deltaT;
        else
%             t(theT+1) = t(theT)+deltaT;
            t(theT+1) = deltaT;
        end
        %extract the times

    end
    timeVec(i+1,:) = t./60; %save in units of minutes
    clear t
end




datastruct  = struct();

datastruct.imageDimensions = dimensions;
datastruct.timeCount= timeCount;
datastruct.channelCount = channelCount;
datastruct.sceneCount = imageCount; 
datastruct.channelNames = channelNameSwapArray;


[~,b] = regexp(FileName,'exp[0-9]');
savename = FileName(1:b);
cd(exportdir)
savename = strcat(savename,'-metaData.mat');
clear omeMeta reader deltaT_Hash theC_Hash theT_Hash FileName
save(savename)





% series1 = data{1, 1};
% series_label_array = cell(1,size(series1,1));
% for p = 1 : size(series1,1)
% series_label_array{p} = series1{p,2};
% end
% 
% export_stack_dir = 'mstack images';
% dirlist = dir(strcat('*',export_stack_dir,'*'));
% if isempty(dirlist)
%     mkdir(export_stack_dir) 
% end
% 
% cd(PathName)
% cd(export_stack_dir)
% for scene = 1 : imageCount
%     sceneNumber = 's00';
%     scenestr = num2str(scene);
%     sceneNumber(end-(length(scenestr)-1):end) = scenestr;
%     
%     img_series = data{scene,1};
%     for channel = 1 : channelCount    
%         channelName = char(omeMeta.getChannelName(0,channel-1));
%         timepoints = channel:channelCount:planeCount;
%         timepoint = timepoints;
%         img_stack = cat(3,img_series{timepoints});
%         savename = strcat(expName,'_',channelName,'_',sceneNumber,'.mat');
%         save(savename,'img_stack')
% %         savefilefunc(savename,img_stack)
%     end
% end
% cd(PathName)
 
end

function savefilefunc(savename,img_stack_test,img_stack_str,export_stack_dir,PathName)
    cd(export_stack_dir)
    save(savename,img_stack_str,'-v7.3');
    cd(PathName)
end

function savefastfilefunc(savename,img_stack_test,img_stack_str,export_stack_dir,PathName)
    cd(export_stack_dir)
    save(savename,img_stack_str,'-v6');
    cd(PathName)
end

function chanstruct = alterChanName(chan)
[a,~] = regexp(chan,'(_|\W|\s|[0-9])'); %remove underscore, dashes, or whitespace
chanstruct = chan;
chanstruct(a) = [];

if isempty(chanstruct)
    chanstruct = 'CFP';
end
end
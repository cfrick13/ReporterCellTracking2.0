
%% fullbkg
function BleachCorrectionOfTimeLapseImages(~,datename,channelstoinput,BACKGROUND)
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
mstackdir = 'flat mstack';
segmentdir = 'segment mstack';
flatfielddirname = 'bleach mstack';
mstackPath = strcat(expPath,'/',mstackdir);
segmentPath = [expPath '/' segmentdir];

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


unstimulatedwellstrings = {'0.00','mock','ctrl','none'};
unstimulatedwellstringsinputs = channelregexpmaker(unstimulatedwellstrings);


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
%%
conditions = dosestructstruct.conditions;
doses = dosestructstruct.doses;
dosestructstruct.doseToScene

condarray = dosestructstruct.conditionsToScenez;
dosearray = dosestructstruct.doseToScenez;
scenearray = dosestructstruct.sceneList;

newstruct = struct();
for i = 1:length(conditions) %conditions
    ii = str2num(condarray{i});
[newstruct(ii).conditions] = deal(conditions(i)) ;
end
for i = 1:length(doses) %doses
    ii = str2num(dosearray{i});
[newstruct(ii).doses] = deal(doses(i)) ;
end
for i = 1:length(scenearray) %scenes
    ii = i;
[newstruct(ii).scenes] = deal(scenearray(i)) ;
dv = newstruct(ii).doses;
cv =newstruct(ii).conditions;
if ~isempty(dv)
newstruct(ii).doseAndCondition = [dv{1} cv{1}];
else
    newstruct(ii).doseAndCondition = 'empty';
    newstruct(ii).conditions ={''};
    newstruct(ii).doses = {''};
end
end

dosearray = [newstruct.doses];
[~,~,~,d] = regexp(dosearray,unstimulatedwellstringsinputs);
ddose = ~cellfun(@isempty,d);
dnamescells = d(ddose);
dnames = cellfun(@(x) x{1},dnamescells,'UniformOutput',0);
unstimulatedidx = find(ddose);

conditionsArray = [newstruct.conditions];
unstimcond = conditions;

%%

%bleach struct
bleachstruct = struct();
%bleachnames

% bleachstruct.(fname) = bleachstackref;

for i = 1:length(unstimcond)
    condstr = unstimcond(i);
    condidx = strcmpi(conditionsArray,condstr);
    cdidx = ddose&condidx;
    sceneArraySub = sceneListArray(cdidx);
    sceneregexp = channelregexpmaker(sceneArraySub);
    [~,~,~,sceneListed] = regexp(dirlistarray,sceneregexp);
    sllog = ~(cellfun(@isempty,sceneListed,'UniformOutput',1));
    
    
    
    
    %%
    
    
    cd(mstackPath)
    %create the median stack of media only images for all timepoints for
    %each channel -- these will be used for flatfield correction in next
    %step
    for channel = channelList %cycle through one channel at a time
        chan = char(channel);
        disp([char(channel) ' - ' condstr])
        [~,~,~,channelsListed] = regexp(dirlistarray,chan);
        cllog = ~(cellfun(@isempty,channelsListed,'UniformOutput',1));
        unstimulatedForGivenChannelArray = dirlistarray(sllog & cllog);
        
        
        
        clear A dirlist
        
        
        cd(segmentPath)
        if strcmpi(chan,'EGFP')
            segdir = dir('*background.mat');
            disp('b')
        else
            segdir = dir('*nucleus*.mat');
            disp('nuc')
        end
        segdirnames = {segdir.name};
        [~,~,~,sceneListed2] = regexp(segdirnames,sceneregexp);
        sidx = ~cellfun(@isempty,sceneListed2);
        backForGivenChannelArray = segdirnames(sidx);
        bkibkg = zeros(dim(1),dim(2),timeCount,length(sceneArraySub),'logical');
        %     bkibkg = nan(dim(1),dim(2),timeCount,length(sceneArraySub),'single');
        for j=1:length(backForGivenChannelArray)
            bkname = char(backForGivenChannelArray{j});
            imgs = load(bkname);
            bkibkg(:,:,:,j) = imgs.IfFinal;
        end
        disp(sceneregexp)
        
        cd(mstackPath)
        
        
        bki = nan(timeCount,length(unstimulatedForGivenChannelArray),'single');
        for j=1:length(unstimulatedForGivenChannelArray)
            bkname = char(unstimulatedForGivenChannelArray{j});
            imgs = load(bkname);
            bkbk = imgs.flatstack;
            bkgb = bkibkg(:,:,:,j);
            for jki = 1:size(bkbk,3)
                bksun = bkbk(:,:,jki);
                bkg = bkgb(:,:,jki);
                if sum(bkg(:))>0
                    bki(jki,j) = nanmedian(bksun(bkg)) - nanmedian(bksun(~bkg));
                    %                 bki(jki,j) = nanmedian(bksun(:));
                elseif sum(bkg(:))==0 && jki==1
                    bki(jki,j) = nanmedian(bksun(:));
                    disp('no background')
                    disp(jki)
                    disp(j)
                else
                    bki(jki,j) = bki(jki-1,j);
                end
            end
        end
        clear imgs
        bkflat = nanmedian(bki,2);
        
        clear bki
        bkflat = single(bkflat);
        
        
        tframe = dosestructstruct.tgfFrames;
        normBkgimgMat = zeros(dim(1),dim(2),timeCount,1,'single');
        for k = 1:timeCount
            bknew = bkflat(k)./bkflat(tframe);
%             bknew = bkflat(k);
            normBkgimgMat(:,:,k) = ones(dim(1),dim(2),'single').*bknew;
        end
        
        chanstruct = alterChanName(chan);
        
        bleachstruct(i).(chanstruct)=normBkgimgMat;
    end
end
clear normBkgimgMat

%%




memoryrequired = dim(1)*dim(2)*timeCount*(length(unstimulatedForGivenChannelArray)+1)*2;%background stack + imagestack

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

% scenesCorr = sceneListArray([1:15 41:45]);
for channel = channelList %cycle through one channel at a time
    filelist = dir(strcat('*',char(channel),'*.mat'));
    cfile1 = {filelist.name};
    [~,~,~,SceneFileNames] = regexp(cfile1,backinputs);
    sfnlog = ~(cellfun(@isempty,SceneFileNames,'UniformOutput',1));
    cfile = cfile1(~sfnlog);
    
    tic
%     for j = 1:length(cfile) %load image_stack for each scene // and determine scene and channel
    parfor j = 1:length(cfile) %load image_stack for each scene // and determine scene and channel
        filename = char(cfile(j));
        imgstruct = load(filename);
        imgstack = single(imgstruct.flatstack);
        
        [a,b] = regexp(filename,'s[0-9]++'); %determine scene
        scene = filename(a:b);
        
        [e,f] = regexp(filename,channelinputs); %determine channel
        chan = filename(e:f);
        chanstruct = alterChanName(chan);
%         normBkgimgMat = bkgimgstruct.(chanstruct);
        
        sidx = strcmpi(sceneListArray,scene);
%         condstr = conditionsArray(sidx);
        condcell = newstruct(sidx).conditions;
        condstr = char(condcell);
        condidx = find(strcmpi(conditions,condstr));
        
        
        normBkgimgMat = bleachstruct(condidx).(chanstruct);
        
        
%         if sum(strcmpi(scenesCorr,scene))>0 && ~strcmpi(chan,'DIC')
%             flatstack = imgstack;
%         else
%             flatstack = imgstack./normBkgimgMat; %flatten the experimental image with the background image
%             
%             flatstack = imgstack./normBkgimgMat; %flatten the experimental image with the background image
%             imgstack = [];
%             normBkgimgMat = [];
%             savename = strcat(E{1},'_flat_',scene,'_','_',chan,'.mat');
%             SAVdir = strcat(expPath,'/',flatfielddirname,'/');
%             savethatimage(savename,SAVdir,mstackPath,flatstack,j);
%         end
        flatstack = imgstack./normBkgimgMat; %flatten the experimental image with the background image
        imgstack = [];
        normBkgimgMat = [];
        savename = strcat(E{1},'_flat_',scene,'_','_',chan,'.mat');
        SAVdir = strcat(expPath,'/',flatfielddirname,'/');
        savethatimage(savename,SAVdir,mstackPath,flatstack,j);
        
    end
    toc
    stopehre=1;
end
stophere=1;

        olddir = pwd;
        SAVdir = strcat(expPath,'/bleachcorrcode/');
        if isdir(SAVdir)
        else
            mkdir(SAVdir)
        end
        cd (SAVdir);
        
        
        savedirfull = pwd;
        copyfile([mdir '.m'],[savedirfull '/bleachcorrCODE.m'])
        cd(olddir)
% lineageyo(A,B);
% cd('D:/Users/zeiss/Documents/MATLAB')
% % bleachcorrection(A,B)
% bleachcorrectionNew(A,B)
% cd('D:/Users/zeiss/Documents/MATLAB')
% toughsegmentationforstacks(A,B)
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



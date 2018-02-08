
%% fullbkg
function FlatfieldCorrectionOfTimeLapseImagesTIFFs(~,datename,channelstoinput,BACKGROUND)
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
mstackdir = 'tifftack images';
flatfielddirname = 'tiflat tifftack';
mstackPath = strcat(expPath,'/',mstackdir);

%assign path to folder for metadata export files
metadir = strcat(parentdir,'Export/');
cd(metadir)



%find associated extracted metadata
filelist  = dir(strcat('*',dateofexp,'*',expnum,'*metaData.mat'));
metadatafile = char(filelist.name);
A = load(metadatafile);
dim = A.dimensions;
timeCount = A.timeCount;
channelinputs =channelregexpmaker(channelstoinput);

cd(mstackPath)
    dirlist = dir('*.tif');
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
    dirlist = dir('*.tif');
    dirlistarray = {dirlist.name};
    [~,~,~,SceneFileNames] = regexp(dirlistarray,backinputs);
    sfnlog = ~(cellfun(@isempty,SceneFileNames,'UniformOutput',1));
    %logical for only background
    

    
    
    %create the median stack of media only images for all timepoints for
    %each channel -- these will be used for flatfield correction in next
    %step
        chan = char(channelList{1});
        [~,~,~,channelsListed] = regexp(dirlistarray,chan);
        cllog = ~(cellfun(@isempty,channelsListed,'UniformOutput',1));
        backgroundForGivenChannelArray = dirlistarray(sfnlog & cllog);

        
        
        clear A dirlist

    
    memoryrequired = dim(1)*dim(2)*timeCount*(length(backgroundForGivenChannelArray)+1)*2;%background stack + imagestack

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

    for channel = channelList %cycle through one channel at a time
%     for channel = {'EGFP'} %cycle through one channel at a time
    filelist = dir(strcat('*',char(channel),'*.tif'));
    cfile = {filelist.name};

        tic
        parfor j = 1:length(cfile) %load image_stack for each scene // and determine scene and channel
            filename = char(cfile(j));
            
            
             info = imfinfo(filename);
            sz1 = info.Height;
            sz2 = info.Width;
            num_images = numel(info);
            imgs = zeros(sz1,sz2,num_images);
            for k = 1:num_images
                imgs(:,:,k) = imread(filename, k, 'Info', info);
            end
            imgs = double(imgs);
%             gimgs = imgs;
            gmimgs = imgs;
            for k = 1:num_images
                gimg = imgaussfilt(imgs(:,:,3),50);
                gmax = prctile(gimg(:),90)*1.5;
                gmimgs(:,:,k) = gimg./gmax;
%                 gimgs(:,:,k) = gimg;
            end
            bimgs = imgs./gmimgs;
%             simgs = imgs-gimg;


            [a,b] = regexp(filename,'s[0-9]++'); %determine scene
            scene = filename(a:b);

            [e,f] = regexp(filename,channelinputs); %determine channel
            chan = filename(e:f);

%                 normBkgimgMat = [];
            flatstack = bimgs;
            savename = strcat(E{1},'_flat_',scene,'_','_',chan,'.tif');
            SAVdir = strcat(expPath,'/',flatfielddirname,'/');
            savethatimage(savename,SAVdir,mstackPath,uint16(flatstack),j);
            
%             flatstack = simgs;
%             savename = strcat(E{1},'_sflat_',scene,'_','_',chan,'.tif');
%             SAVdir = strcat(expPath,'/',flatfielddirname,'/');
%             savethatimage(savename,SAVdir,mstackPath,uint16(flatstack),j);
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
% save(savename,'flatstack','-v6');



imwrite( flatstack(:,:,1), savename,'WriteMode','overwrite');
for ind = 2:size(flatstack,3)
    imwrite( flatstack(:,:,ind), savename,'WriteMode','append');
end




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



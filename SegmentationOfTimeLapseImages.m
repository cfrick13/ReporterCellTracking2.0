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
    
    
    fnames = fieldnames(segInstruct);
    segInstructList = cell(1,length(fnames));
    for i = 1:length(segInstructList)
        str = fnames{i};
        segInstructList{i} = str;
    end




%define parameter structure and default parameter values
pStruct = defaultpStructFunc(segInstructList);

pStruct = loadSegmentParameters(pStruct,datename,exportdir); %loads saved value of pStruct


bkarray = bkarraymaker(BACKGROUND);
bkinputs =channelregexpmaker(bkarray);

cd(mstackPath)
nucleusFileList = dir(mstackPath);
[~,~,~,d] = regexp({nucleusFileList.name},'s[0-9]++');
dlog = ~cellfun(@isempty,d,'UniformOutput',1); 
dcell = d(dlog);
SceneList = unique(cellfun(@(x) x{1},dcell,'UniformOutput',0));
nucleusFileListArrayPre = {nucleusFileList.name};
nucleusFileListArray = nucleusFileListArrayPre(dlog);
    
[~,~,~,d] = regexp(SceneList,bkinputs);
bkgscenelog = cellfun(@isempty,d,'UniformOutput',1);
SceneList = SceneList(bkgscenelog);

fileobject = matfile(nucleusFileListArray{1});
imgstack = fileobject.flatstack;
    %determine number of parallel cores to employ based on memory requirements
        dim = size(imgstack);
        if max(size(dim))>2
        else
            newdim = [dim(:)' 1];
            dim = newdim;
        end
        memoryrequired = dim(1)*dim(2)*dim(3)*2*2;%background stack + imagestack
        memoryrequired = dim(1)*dim(2)*dim(3)*2*2*4;%background stack + imagestack + extra processes
%             disp(strcat('memory required for one image =',num2str(round(memoryrequired./(1e6),2,'decimals')),'MegaBytes'))
        detWorkers = 6./(memoryrequired./1e9); %determine the number of workers you can use while staying under 6 GB
        nWorkers = floor(detWorkers);
        possibleWorkers = feature('numcores');
        if nWorkers>possibleWorkers
            nWorkers = possibleWorkers;
        elseif nWorkers <1
            nWorkers = 1;
        end
        

        disp(['memory required for imagestack = ' num2str(round(memoryrequired./(1e6),0,'decimals')) ' MB'])

        
%      Enter parallel loop
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
    sList = SceneList;
%     exclude={'s01','s02','s03','s08','s09','s10','s16','s17','s22','s23','s24','s29','s30','s31'};
%     smem = ismember(sList,exclude);
%     sList(~smem)=[];
    sceneVector = 1:length(sList);
    fractions = length(sceneVector)./nWorkers;
    firstlengths = ceil(fractions);
    if firstlengths*nWorkers > length(sList)
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
    
    parfor iNW=1:nWorkers
%     for iNW=1:nWorkers
        sceneArrayVec = sceneArray{iNW};
        sub_sList = sList(sceneArrayVec);
        for i = 1:length(sceneArrayVec)
            scenestr = char(sub_sList{i});

        %locate nuclei data set
            nuctic = tic;
            cd(mstackPath)
            nucleusFileList = dir(strcat('*',scenestr,'*',nucleus_seg,'*.mat'));
            if isempty(nucleusFileList)
                error(strcat('the channel set for "nucleus_seg" ("',nucleus_seg,'") does not exist'))
            end
            %load nuclei data set
            nucleusFileName = char(nucleusFileList.name);
            fileObject = matfile(nucleusFileName);
            FinalImage = fileObject.flatstack;
            %segment nuclei
            [~,~] = segmentationNucleus(FinalImage,segmentPath,'nucleus',nucleusFileName,pStruct);
            nuctoc = num2str(round(toc(nuctic),0,'decimals'));


        %locate cell fluorescence data set for background segmentation
            backtic = tic;
            cd(mstackPath)
            backgroundFileList = dir(strcat('*',scenestr,'*',background_seg,'*.mat'));
            if isempty(backgroundFileList)
                error(strcat('the channel set for "background_seg" ("',background_seg,'") does not exist'))
            end
            %load cell fluorescence data set
            backgroundFileName = char(backgroundFileList.name);
            fileObject = matfile(backgroundFileName);
            FinalImage = fileObject.flatstack;
            %segment background
            [~,~] = segmentationImageBackground(FinalImage,segmentPath,'background',backgroundFileName,pStruct);  
            cd ..
            backtoc = num2str(round(toc(backtic),0,'decimals'));
            
            disp([scenestr ' nucSeg time= ' nuctoc ' s , backSeg time= ' backtoc ' s'])
            ssss=1;
            %%
        end
    end
end


function pStruct = defaultpStructFunc(segInstructList)
    pStruct = struct();
    parameterDefaults.background = [30 1 2 0.5 10 10];
    parameterDefaults.nucleus = [30 1 2 0.5 10 10];
    parameterDefaults.cell = [40 1 2 0.5 10 10];
    parameterStrings = {'nucDiameter','threshFactor','sigmaScaledToParticle','metthresh','percentSmoothed','denoise'};
    for p = 1:length(parameterStrings)
        pString = char(parameterStrings{p});
        for c = 1:length(segInstructList)
            cstr = char(segInstructList{c});
            cString = alterChanName(cstr);
            pd = parameterDefaults.(cString);
            pStruct.(cString).(pString) = pd(p); 
        end
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

function [IfFinal,testOut] = segmentationNucleus(FinalImage,segmentPath,nucleus_seg,nucleusFileName,pStruct)
    testOut = struct();                           
    img = FinalImage(:,:,1);
    
    
%     %determine number of parallel cores to employ based on memory requirements
%         dim = size(img);
%         memoryrequired = dim(1)*dim(2)*2;%background stack + imagestack
% %             disp(strcat('memory required for one image =',num2str(round(memoryrequired./(1e6),2,'decimals')),'MegaBytes'))
%         detWorkers = 6./(memoryrequired./1e9); %determine the number of workers you can use while staying under 6 GB
%         nWorkers = floor(detWorkers);
%         possibleWorkers = feature('numcores');
%         if nWorkers>possibleWorkers
%             nWorkers = possibleWorkers;
%         elseif nWorkers <1
%             nWorkers = 1;
%         end
%         
%      % Enter parallel loop
%         poolobj = gcp('nocreate');
%         if isempty(poolobj)
%             poolobj = parpool(nWorkers);
%         else
%         nw = poolobj.NumWorkers;
%             if nw ==nWorkers
%             else
%                 delete(poolobj)
%                 poolobj = parpool(nWorkers);
%             end
%         end
    
%     start segmentation
%     parfor frames = 1:size(FinalImage,3)
    IfFinal = false(size(FinalImage));
    for frames = 1:size(FinalImage,3)
        img = FinalImage(:,:,frames); 
        [If,~] = segmentNuclei(img,nucleus_seg,pStruct,frames);
        IfFinal(:,:,frames)=If;
    end
                
    nfilename = [nucleusFileName(1:end-4) '_' nucleus_seg '.mat'];
    savethatimagestack(IfFinal,nfilename,segmentPath)
end


function [IfFinal,testOut] = segmentationImageBackground(FinalImage,segmentPath,background_seg,backgroundFileName,pStruct)
    testOut = struct();                           
    img = FinalImage(:,:,1);
    
%     %determine number of parallel cores to employ based on memory requirements
%         dim = size(img);
%         memoryrequired = dim(1)*dim(2)*2;%background stack + imagestack
% %             disp(strcat('memory required for one image =',num2str(round(memoryrequired./(1e6),1,'decimals')),'MegaBytes'))
%         detWorkers = 6./(memoryrequired./1e9); %determine the number of workers you can use while staying under 6 GB
%         nWorkers = floor(detWorkers);
%         possibleWorkers = feature('numcores');
%         if nWorkers>possibleWorkers
%             nWorkers = possibleWorkers;
%         elseif nWorkers <1
%             nWorkers = 1;
%         end
%         
%      % Enter parallel loop
%         poolobj = gcp('nocreate');
%         if isempty(poolobj)
%             poolobj = parpool(nWorkers);
%         else
%         nw = poolobj.NumWorkers;
%             if nw ==nWorkers
%             else
%                 delete(poolobj)
%                 poolobj = parpool(nWorkers);
%             end
%         end

    IfFinal = false(size(FinalImage));
%     parfor frames = 1:size(FinalImage,3)
    for frames = 1:size(FinalImage,3)
%         disp(frames)
        img = FinalImage(:,:,frames); 
        [If,~] = segmentCellBackground(img,background_seg,pStruct,frames);
        IfFinal(:,:,frames)=If;
    end
bfilename = [backgroundFileName(1:end-4) '_' background_seg '.mat'];                
savethatimagestack(IfFinal,bfilename,segmentPath)


% parameters

end


function savethatimagestack(IfFinal,filename,segmentPath)
olddir  = pwd;
cd(segmentPath)
save(filename,'IfFinal','-v7.3');
cd (olddir)


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

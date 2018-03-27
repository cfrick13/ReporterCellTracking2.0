
%% fullbkg
function MovieMakerBaby(~,datename,channelstoinput,stimulationFrame)
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
moviedirname = 'moviestack';

mstackPath = strcat(expPath,'/',mstackdir);
flatstackPath = [expPath '/' flatfielddirname];
moviestackPath = [expPath '/' moviedirname];

%assign path to folder for metadata export files
metadir = strcat(parentdir,'Export/');
cd(metadir)



%find associated extracted metadata
filelist  = dir(strcat('*',dateofexp,'*',expnum,'*metaData.mat'));
metadatafile = char(filelist.name);
A = load(metadatafile);
dim = A.dimensions;
timeCount = A.timeCount;
timeMatrix = A.timeVec;
channelinputs =channelregexpmaker(channelstoinput);

cd(flatstackPath)
dirlist = dir('*.tif');
[~,~,~,channelsListed] = regexp([dirlist.name],channelinputs);
[~,~,~,sceneListArray] = regexp([dirlist.name],'s[0-9]+');
channelList = unique(channelsListed);
sceneList = unique(sceneListArray);


memoryrequired = dim(1)*dim(2)*timeCount*(length(timeCount)+1)*2;%background stack + imagestack

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

cd(flatstackPath)

sceneVec = {'s02','s14','s20','s32','s40','s46'};
conditionsArray = {'no transfection','mCherry-Smad3 mRNA','mCherry-Smad3 mRNA','mock transfection','mCherry-Smad5 mRNA','mCherry-Smad5 mRNA'};
LigandArray = {'Tgf','Tgf','no ligand','Tgf','Tgf','BMP-4'};


parfor j = 1:length(sceneVec)
%     for c = 1:length(channelList) %cycle through one channel at a time
    tic
    for c = 1:3 %cycle through one channel at a time
        channel = channelList{c};
     
        scenestr = sceneVec{j};
        filelist = dir(strcat('*',scenestr,'*',channel,'*.tif'));
        filename = char(filelist.name);
        
        [a,b] = regexp(filename,'s[0-9]++'); %determine scene
        scene = filename(a:b);
        sidx = strcmpi(sceneList,scene);
        
        [e,f] = regexp(filename,channelinputs); %determine channel
        chan = filename(e:f);
        
        
        
        info = imfinfo(filename);
        sz1 = info.Height;
        sz2 = info.Width;
        num_images = numel(info);
        imgs = zeros(sz1,sz2,num_images);
%         for k = 1:num_images
%             imgs(:,:,k) = imread(filename, k, 'Info', info);
%         end
%         imgs = double(imgs);
        
        

        
        timeVec = round(timeMatrix(sidx,:));
        timeVec = timeVec - timeVec(stimulationFrame);
        
        expdate = '2018_03_07 plate james exp1';
        fluor = 'NG-Smad3 (endogenous)';
        condition = conditionsArray{j};
        CHANNEL = char(channel);
        if strcmp(channel,'EGFP')
            CHANNEL = 'endogenous NG-Smad3';
        elseif strcmp(channel,'mKate') && ~isempty(regexpi(condition,'Smad5'))
            CHANNEL = 'mCherry-Smad5 (mRNA)';
            timeVec = round(timeMatrix(sidx,:));
            timeVec = timeVec - timeVec(15);
        elseif strcmp(channel,'mKate') && ~isempty(regexpi(condition,'Smad3'))
            CHANNEL = 'mCherry-Smad3(mRNA)';
        elseif strcmp(channel,'DIC')
            CHANNEL = 'DIC';
        end
        
        ligand = LigandArray{j};
        str = scenestr;
        detailsArray = {str,['(' condition ')'],ligand,CHANNEL};
%         detailsArray = {str,['(' condition ')'],ligand,CHANNEL,fluor};
        detY = linspace(340,480,length(detailsArray))-80;
        fontSize = 34;
        
        %% save DetailsImage
%         f89 = figure(89);
%         f89.Units = 'pixels';
%         f89.Position = [500 500 dim(1) dim(2)];
%         colormap(f89,'gray');
%         f89.Color = 'k';
%         f89.InvertHardcopy = 'off';
        imgs = zeros(sz1,sz2,num_images);
        for k = 1:num_images
            %%
            tnum=k;
            img = zeros(dim(1),dim(2),3,'uint16');
            
%             II = imagesc(img);
% 
%             ax = II.Parent;
%             ax.Units = 'pixels';
%             ax.YTick =[];
%             ax.XTick = [];
%             ax.Position = [0 0 dim(1) dim(2)];
%             II.delete;
%             ax.Color = 'k';

            
            for jj = 1:length(detailsArray)
                detailstr = detailsArray{jj};
%                 t = text(ax,0,0,detailstr);
%                 t.Units = 'pixels';
%                 t.Position = [5 detY(jj)];
%                 t.Color = 'w';
%                 t.FontSize = fontSize;
%                 t.FontWeight ='bold';
                img = insertText(img,[5 dim(2)-detY(jj)],detailstr,'TextColor','white','FontSize',fontSize,'BoxColor','black','AnchorPoint','LeftBottom','Font','Consolas','BoxOpacity',0);
            end
            
            if timeVec(tnum)<0
                hr = ceil(timeVec(tnum)./60);
                mint = abs(round(timeVec(tnum) - hr*60));
                isneg = true;
            else
                hr = floor(timeVec(tnum)./60);
                mint = round(timeVec(tnum) - hr*60);
                isneg = false;
            end
            
            minstr = '00';
            %     minsub = num2str(round(timeVec(tnum)./60,1,'decimals'));
            minsub = num2str(mint);
            minidx = false(size(minstr));
            minidx(end-(length(minsub)-1):end) = true;
            if timeVec(tnum)<0
                minstr(end-(length(minsub)-1):end) = minsub;
                %         minstr(~minidx) = ' ';
            else
                minstr(end-(length(minsub)-1):end) = minsub;
                
                %         minstr(~minidx)=' ';
            end
            
            hrstr = '00';
            %     minsub = num2str(round(timeVec(tnum)./60,1,'decimals'));
            minsub = num2str(hr);
            minidx = false(size(hrstr));
            minidx(end-(length(minsub)-1):end) = true;
            if timeVec(tnum)<0
                hrstr(end-(length(minsub)-1):end) = minsub;
                hrstr(~minidx) = ' ';
                hrstr = num2str(abs(hr),'%02d');
                hrstr(1) = '-';
            else
                hrstr(end-(length(minsub)-1):end) = minsub;
                hrstr = num2str((hr),'%02d');
                if ~(hr>=10)
                    hrstr(1) = [' '];
                end
                %         hrstr(~minidx)='';
            end
            
            timestr = [hrstr ':' minstr ' hr:min'];
%             frametext= text(ax,0,0,[hrstr ':' minstr ' hr:min']);
            img = insertText(img,[5 dim(2)-10],timestr,'TextColor','white','FontSize',fontSize,'BoxColor','black','AnchorPoint','LeftBottom','Font','Consolas');
%             frametext.Units = 'pixels';
%             frametext.Position = [5 30];
%             frametext.FontSize = fontSize;
%             frametext.Color = 'w';
%             frametext.FontWeight = 'bold';
%             frametext.HorizontalAlignment = 'left';
%             frametext.BackgroundColor = 'k';
            
            if ~(timeVec(tnum)<0)
                ligandstr = [ligand];
%                 frametext= text(ax,0,0,ligandstr);
%                 frametext.Units = 'pixels';
%                 frametext.Position = [5 dim(2)-80];
%                 frametext.FontSize = fontSize;
%                 frametext.Color ='w';
%                 frametext.FontWeight = 'bold';
%                 frametext.HorizontalAlignment = 'left';
%                 frametext.BackgroundColor = 'k';
                img = insertText(img,[5 dim(2)-60],ligandstr,'TextColor','white','FontSize',fontSize,'BoxColor','black','AnchorPoint','LeftBottom','Font','Consolas');
            end
            
            imgflat = uint16(mean(img,3));
            imgs(:,:,k) = imgflat;
%             texxx = findobj('Type','Text');
%             [texxx.FontName] = deal('consolas');
            tstr = num2str(k,'t%04d');
%            SAVdir = strcat(expPath,'/',moviedirname,'/');
%            savename = strcat(E{1},'_movie_',scene,'_',channel,'_',tstr,'.tif');
% %            savethatimage(savename,SAVdir,moviestackPath,uint16(img),j);
%            olddir = pwd;
%            if isdir(SAVdir)
%            else
%                mkdir(SAVdir)
%            end
%            cd (SAVdir);
           
%            disp([savename tstr])
%            saveas(f89,[savename '.tif'],'tif')
           
%            cd(olddir)
        end
        
        olddir = pwd;
        SAVdir = strcat(expPath,'/',moviedirname,'/');
        if isdir(SAVdir)
        else
            mkdir(SAVdir)
        end
        cd (SAVdir);

        

        savename = strcat(E{1},'_movie_',scene,'_',channel,'.tif');
        savethatimage(savename,SAVdir,moviestackPath,uint16(imgs),j);
        cd(olddir)
        

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



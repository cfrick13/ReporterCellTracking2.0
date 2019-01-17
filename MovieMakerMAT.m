
%% fullbkg
function MovieMakerMAT(~,datename,channelstoinput,stimulationFrame)
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
dirlist = dir('*.mat');
[~,~,~,channelsListed] = regexp([dirlist.name],channelinputs);
[~,~,~,sceneListArray] = regexp([dirlist.name],'s[0-9]+');
channelList = unique(channelsListed);
sceneList = unique(sceneListArray);



cd(flatstackPath)

% sceneVec = {'s29','s36'};
% conditionsArray = {'CMV-NG-Smad3 (NMuMG)','NG-Smad3(NMuMG)'};
% LigandArray = {'Tgf','Tgf'};

sceneVec = {'s07','s12','s26','s33','s52'};
conditionsArray = {'C2C12','C2C12','C2C12','C2C12','C2C12'};
LigandArray = {'BMP-10','BMP-10','BMP-10','BMP-10','BMP-10'};
smadArray = {'Citrine-Smad5','Citrine-Smad1','Citrine-Smad1','Citrine-Smad5','Citrine-Smad3'};
scarlettArray = {'Scarlett-Smad1','Scarlett-Smad5','Scarlett-Smad1','Scarlett-Smad5','Scarlet-Smad1'};


for j = 1:length(sceneVec)
%     for c = 1:length(channelList) %cycle through one channel at a time
    tic
    for c = 1:length(channelList) %cycle through one channel at a time
        channel = channelList{c};
     
        scenestr = sceneVec{j};
        filelist = dir(strcat('*',scenestr,'*',channel,'*.mat'));
        filename = char(filelist.name);
        
        [a,b] = regexp(filename,'s[0-9]++'); %determine scene
        scene = filename(a:b);
        sidx = strcmpi(sceneList,scene);
        
        [e,f] = regexp(filename,channelinputs); %determine channel
        chan = filename(e:f);
        
        
        mfile = matfile(filename);
        imgs1 = mfile.flatstack;
%         imgs1 = imgs1(:,:,1:num_images);
%         imgs1 = double(imgs1);
        
        
        
        

        
        timeVec = round(timeMatrix(sidx,:));
        timeVec = timeVec - timeVec(stimulationFrame);
        
        datename = '2018_04_25 plate citrine scarlett james mrna exp1';
%         fluor = 'NG-Smad3 (endogenous)';
        condition = conditionsArray{j};
        CHANNEL = char(channel);
        if strcmp(channel,'EGFP')
            CHANNEL = [smadArray{j} ' (mRNA)'];
        elseif strcmp(channel,'mKate')
            CHANNEL = [scarlettArray{j} ' (mRNA)'];
        end
        

        ligand = LigandArray{j};
        str = scenestr;
        detailsArray = {str,['(' condition ')'],ligand,CHANNEL};
%         detailsArray = {str,['(' condition ')'],ligand,CHANNEL,fluor};
        detY = linspace(340,480,length(detailsArray))-80;
        fontSize = 34;
        
        dim = size(imgs1);
        imgs = zeros(dim);
        for k = 1:dim(3)
            %%
            tnum=k;
            img = zeros(dim(1),dim(2),3,'uint16');
            
            for jj = 1:length(detailsArray)
                detailstr = detailsArray{jj};
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
             minsub = num2str(mint);
            minidx = false(size(minstr));
            minidx(end-(length(minsub)-1):end) = true;
            if timeVec(tnum)<0
                minstr(end-(length(minsub)-1):end) = minsub;
 
            else
                minstr(end-(length(minsub)-1):end) = minsub;

            end
            
            hrstr = '00';
 
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
            img = insertText(img,[5 dim(2)-10],timestr,'TextColor','white','FontSize',fontSize,'BoxColor','black','AnchorPoint','LeftBottom','Font','Consolas');

            if ~(timeVec(tnum)<0)
                ligandstr = [ligand];
                img = insertText(img,[5 dim(2)-60],ligandstr,'TextColor','white','FontSize',fontSize,'BoxColor','black','AnchorPoint','LeftBottom','Font','Consolas');
            end
            
            imgflat = uint16(mean(img,3));
            imgs(:,:,k) = imgflat;
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
        savedirfull = pwd;
        
        %%
        imgsnew = horzcat(imgs1,imgs);
        imgsnew(imgsnew>prctile(imgs1(:),98)) = prctile(imgs1(:),98);
        imgsnew(imgsnew<prctile(imgs1(:),1)) = prctile(imgs1(:),1);
        
        specialdir = [expPath '/AVIs'];
        if ~isdir(specialdir)
            mkdir(specialdir)
        end
        savename = [specialdir '/' scene '_' channel '.avi'];
        profile = 'Motion JPEG AVI';
        v = VideoWriter(savename,profile);
        % t = frames/(frames/s);
        frameRate = 5;
        % v.Duration = length(nannum(1):nannum(end))./frameRate;
        v.FrameRate = frameRate;
        v.Quality = 95;
        
        open(v);
        
        f997 = figure(997);
        h = axes();
        colormap(f997,'gray');
        h.Position = [0 0 1 1];

%         h.NextPlot = 'replacechildren';
        h.YDir='reverse';
        f997.Units = 'pixels';
        f997.Position = [200 200 size(imgsnew,2) size(imgsnew,1)];
        for i = 1:size(imgsnew,3)
            imagesc(h,imgsnew(:,:,i));
            h.XTick=[];
            h.YTick=[];
            h.CLim = [prctile(imgsnew(:),2) prctile(imgsnew(:),99)];
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        
        close(v);
        %%



        

        savename1 = strcat(E{1},'_movie_',scene,'_',channel,'.tif');
        savethatimage(savename1,SAVdir,moviestackPath,uint16(imgs),j);
        copyfile([mdir '.m'],[savedirfull 'moviemakerCode.m'])
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



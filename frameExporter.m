
%% fullbkg
function frameExporter(~,datename,channelstoinput,stimulationFrame)
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
[a,b] = regexp(datename,'exp[0-9]+');
expnum = datename(a:b);
dateofexp = datename(1:10);

%assign path to .mat image stack exports
mstackdir = 'mstack images';
flatfielddirname = 'flat mstack';
moviedirname = 'moviestack';
segdirname = 'segment mstack';

mstackPath = strcat(expPath,'/',mstackdir);
flatstackPath = [expPath '/' flatfielddirname];
segstackpath = [expPath '/' segdirname];
moviestackPath = [expPath '/' moviedirname];

%assign path to folder for metadata export files
metadir = strcat(parentdir,'Export/');
cd(metadir)


%%
FileName =datename;

[a,~] = regexp(datename,'_nuclei');
datequery = strcat(datename,'*metaData.mat');
filedatestr = datename(1:a-1);
cd(metadir)
load([datename '_nuclei_export.mat'])
filelist = dir(datequery);
if length({filelist.name}) ==1
    metaData = load(char(filelist.name));
else
    filename = uigetfile();
    metaData = load(filename);
end
%determine the timeVector from metaData [dim 1 is scene#, dim 2 is each time frame]
timeVec = metaData.timeVec;

%load information regarding doses and scenes and tgfbeta addition
% [a,~] = regexp(FileName,'_nuclei');
datequery = strcat(FileName,'*DoseAndScene*');
cd(metadir)
filelist = dir(datequery);
if isempty(filelist)
    dosestruct = makeDoseStruct; %run function to make doseStruct
else
    dosestructstruct = load(char(filelist.name));
    dosestruct = dosestructstruct.dosestruct;
end



%finished collecting metaData
%now combine metaDeta together



%determine the scenes present in the experiment in order to combine metadata
scenestr = 'scene';
sceneListArray = vertcat({exportNucleiStruct.(scenestr)});
sceneList = unique(sceneListArray);
sceneListArrayTwo = vertcat({dosestruct.(scenestr)});

%combine the exportStruct information with dosestruct information
for i=1:length(sceneList)
    sceneChoice=sceneList{i};
    indices = strcmp(sceneListArray,sceneChoice);
    indicestwo = strcmp(sceneListArrayTwo,sceneChoice);
    
    fnamesz = fieldnames(dosestruct);
    for p = 1:length(fnamesz)
        fnstr = fnamesz{p};
        var = dosestruct(indicestwo).(fnstr);
        if isnumeric(var)
            [exportNucleiStruct(indices).(fnstr)] = deal(var);
        else
            [exportNucleiStruct(indices).(fnstr)] = deal(char(var));
        end
    end
end

%combine the exportStruct information with dosestruct information
for i=1:length(sceneList)
    sceneChoice=sceneList{i};
    indices = strcmp(sceneListArray,sceneChoice);
    indicestwo = strcmp(sceneListArrayTwo,sceneChoice);
    
    
    expdate = FileName(1:a-1);
    fnamesz = fieldnames(dosestruct);
    for p = 1:length(fnamesz)
        fnstr1 = 'dosestr';
        fnstr2 = 'conditions';
        fnstr3 = 'expdate';
        var1 = dosestruct(indicestwo).(fnstr1);
        var2 = dosestruct(indicestwo).(fnstr2);
        [exportNucleiStruct(indices).doseAndCondition] = deal([char(var1) char(var2)]);
        [exportNucleiStruct(indices).doseconddate] = deal([expdate ' ' char(var1) ' ' char(var2)]);
        [exportNucleiStruct(indices).(fnstr3)] = deal(expdate);
    end
end
doseListArray = vertcat({exportNucleiStruct.dosestr});



%%

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

sceneVec = {'s005','s020','s035','s050','s065','s080','s095','s110'};
sceneVec = {'s05','s10','s15','s20','s25','s30'};
% conditionsArray = {'C2C12','C2C12','C2C12','C2C12','C2C12'};
% LigandArray = {'BMP-10','BMP-10','BMP-10','BMP-10','BMP-10'};
% smadArray = {'Citrine-Smad5','Citrine-Smad1','Citrine-Smad1','Citrine-Smad5','Citrine-Smad3'};
% scarlettArray = {'Scarlett-Smad1','Scarlett-Smad5','Scarlett-Smad1','Scarlett-Smad5','Scarlet-Smad1'};
%%
close all

%scaled individually
f997 = figure(997);
colormap(f997,'gray');
f997.Units = 'pixels';
f997.Position = [200 200 2045 512];
% 
cycle=0;
for j = 1:length(sceneVec)
    for c = 1:length(channelList) %cycle through one channel at a time
        cycle=cycle+1;
        h(cycle) = axes();
        pos = zeros(1,4);
        
        pos(1) = (j-1)./length(sceneVec);
        pos(2) = (c-1)./length(channelList);
        pos(3) = 1./length(sceneVec);
        pos(4) = 1./length(channelList);
        
        h(cycle).Position =pos;
        
    end
end
[h.XTick]=deal([]);
[h.YTick]=deal([]);
[h.YDir] = deal('reverse');
[h.NextPlot] = deal('replacechildren');


%scaled identically
f998 = figure(998);
colormap(f998,'gray');
f998.Units = 'pixels';
f998.Position = [200 200 2045 512];
cycle=0;
for j = 1:length(sceneVec)
    for c = 1:length(channelList) %cycle through one channel at a time
        cycle=cycle+1;
        hh(cycle) = axes();
        pos = zeros(1,4);
        
        pos(1) = (j-1)./length(sceneVec);
        pos(2) = (c-1)./length(channelList);
        pos(3) = 1./length(sceneVec);
        pos(4) = 1./length(channelList);
        
        hh(cycle).Position =pos;
        
    end
end
[hh.XTick]=deal([]);
[hh.YTick]=deal([]);
[hh.YDir] = deal('reverse');
[hh.NextPlot] = deal('replacechildren');


%
cd(flatstackPath)
cycle =0;
for j = 1:length(sceneVec)
    tic
    for c = 1:length(channelList) %cycle through one channel at a time
        cycle =cycle+1;
        channel = channelList{c};
        
        scenestr = sceneVec{j};
        scenenum = str2double(scenestr(2:end));
        dandc=exportNucleiStruct(scenenum).doseAndCondition;
        
        
        disp(dandc);
        filelist = dir(strcat('*',scenestr,'*',channel,'*.mat'));
        filename = char(filelist.name);
        
        
        mfile = matfile(filename);
        imgs1 = mfile.flatstack;
        
        cd(segstackpath)
         filelist = dir(strcat('*',scenestr,'*',channel,'*background*.mat'));
        filename = char(filelist.name);
        
        
        mfile = matfile(filename);
        imgs1back = mfile.IfFinal;
        
        for kk = 1:size(imgs1,3)
            imm = imgs1(:,:,kk);
            immb = imgs1back(:,:,kk);
            imgs1(:,:,kk) = imm-median(median(imm(~immb)));
        end
        cd(flatstackPath)
        
        olddir = pwd;
        SAVdir = strcat(expPath,'/',moviedirname,'/');
        if isdir(SAVdir)
        else
            mkdir(SAVdir)
        end
        cd (SAVdir);
        
        %
        imgsnew = imgs1(100:400,100:400);
        imgsnew(imgsnew>prctile(imgs1(:),98)) = prctile(imgs1(:),98);
        imgsnew(imgsnew<prctile(imgs1(:),1)) = prctile(imgs1(:),1);
        
        
        
        specialdir = [expPath '/FocusFrames'];
        if ~isdir(specialdir)
            mkdir(specialdir)
        end
        
        imagesc(h(cycle),imgsnew(:,:,end));
        h(cycle).CLim = [prctile(imgsnew(:),2) prctile(imgsnew(:),99)];
        h(cycle).Title.String = scenestr;
        h(cycle).XLim = [0 size(imgsnew,1)];
        h(cycle).YLim = [0 size(imgsnew,2)];
        colorv = [0 1];
        posv = [0.05 0];
        for jo = 1:2
            t=text(h(cycle),0,0,'r');
            t.Position=[0,posv(jo)];
            t.String = [scenestr '-' dandc];
            t.Color = [colorv(jo) 0 0];
            t.FontSize = 15;
            t.VerticalAlignment = 'top';
            t.FontWeight = 'bold';
        end
        cd(olddir)
        
        
        
        imgs2 = imgs1(100:400,100:400,:);
        %background subtraction
        imagesc(hh(cycle),imgs2(:,:,1));
        
%         h(cycle).CLim = [prctile(imgsnew(:),2) prctile(imgsnew(:),99)];
        hh(cycle).CLim = [0 10000];%set CLIM
        hh(cycle).Title.String = scenestr;
        hh(cycle).XLim = [0 size(imgs2,1)];
        hh(cycle).YLim = [0 size(imgs2,2)];
        colorv = [0 1];
        posv = [0.05 0];
        for jo = 1:2
            t=text(hh(cycle),0,0,'r');
            t.Position=[0,posv(jo)];
            t.String = [scenestr '-' dandc];
            t.Color = [colorv(jo) 0 0];
            t.FontSize = 15;
            t.VerticalAlignment = 'top';
            t.FontWeight = 'bold';
        end
        cd(olddir)
        
        
    end
    toc
end
% cd(specialdir)
cd('I:\Frick')
saveas(f997,[E{1} '-focustest.tif'],'tiff');
saveas(f998,[E{1} '-scaledidentically.tif'],'tiff');

end



function savethatimage(savename,SAVdir,mstackPath,flatstack,j)
disp(strcat(savename,'...',num2str(j)));
cd (SAVdir);
%     imwrite(flatstack,char(savename),'tiff');
%     save(flatstack,char(savename),'tiff');
%     save(savename,'flatstack','-v7.3');
% save(savename,'flatstack','-v6');



imwrite( flatstack(:,:,1), savename,'WriteMode','overwrite');





cd (mstackPath);
end



function channelinputs =channelregexpmaker(channelstoinput)
channelinputs = '(';
for i=1:length(channelstoinput) % creates a string of from '(c1|c2|c3|c4)' for regexp functions
    if i ==1 & length(channelstoinput)==1
        channelinputs = channelstoinput{i};
    elseif i ==1
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



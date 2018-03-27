function uiTrackCellzNEW(FileDate,AutoTrackStr)
global  stackstruct channelstoinput medmeantot mainfig stimulationFrame Tracked chanbkgmedianmat nucbkgmedianmat pStruct timeVec expDetailsStruct dirStruct timeSteps DivisionStruct segStruct foStruct xAxisLimits DICimgstack pathStruct segmentimgstack channelimgstack ExportNameKey ExportName plottingTotalOrMedian channelinputs tcontrast lcontrast ThirdPlotAxes SecondPlotAxes togStruct PlotAxes cmap  timeFrames frameToLoad ImageDetails MainAxes SceneList

close all

pathStruct = struct();
DivisionStruct = struct();
segStruct = struct();
foStruct = struct();
dirStruct = struct();
expDetailsStruct = struct();
ImageDetails = struct();
togStruct = struct();

%determine matfile directory
mdir = mfilename('fullpath');
[~,b ] = regexp(mdir,'/');
if isempty(b)
    [~,b] = regexp(mdir,'\');
end
dirStruct.parentdir = mdir(1:b(end)); %folder in which matfile exists
dirStruct.exportdir = strcat(dirStruct.parentdir,'Export/');

[~,b ] = regexp(dirStruct.parentdir,'/');
if isempty(b)
    [~,b] = regexp(dirStruct.parentdir,'\');
end
dirStruct.gparentdir = dirStruct.parentdir(1:b(end-1)); %folder in which matfile dirStruct.parentdir exists


%specify important directory names
mstackName = 'flat mstack';
trackName = 'tracking files';
segmentName = 'segment mstack';


%set contrast and plotting details
togStruct.updateContrast=false;
togStruct.runIterate =false;
togStruct.trackUpdated = true;
togStruct.plotSettingsToggle =false;
togStruct.displayTrackingToggle = false;
togStruct.SmallMovie = false;
togStruct.SaveMovie = false;

plottingTotalOrMedian = 'median';
tcontrast = 99;
lcontrast = 1;
medmeantot = 'total';


%determine export details
ExportNameKey = 'tsi_';
% ExportNameKey = 'final';
if ~strcmp(ExportNameKey,'final')
    disp(strcat('Export name key is "',ExportNameKey,'" not FINAL'))
end
ExportName = 'fricktrack_num_mmt';


%initialize global variables
channelimgstack =[];
segmentimgstack =[];
DICimgstack=[];




% user selects experiment directory to be analyzed
cd(dirStruct.gparentdir)
if isempty(FileDate)
    pathStruct.expDirPath = uigetdir;
else
    cd(dirStruct.gparentdir)
    expdir = strcat(dirStruct.gparentdir,'/',FileDate);
    cd(expdir)
    pathStruct.expDirPath = pwd;
end

cd(pathStruct.expDirPath)
experimentdir = pathStruct.expDirPath;
pathStruct.mstackPath = strcat(experimentdir,'/',mstackName);
pathStruct.segmentPath = strcat(experimentdir,'/',segmentName);
pathStruct.trackingPath = strcat(experimentdir,'/',trackName);
pathStruct.mdirPath = mdir;


%make tracking file folder if it does not exist
cd(experimentdir)
dirlist = dir(trackName);
if isempty(dirlist)
    mkdir(trackName);
end


%determine date of experiment
[a,b] = regexp(pathStruct.expDirPath,'201[0-9]');
[~,d] = regexp(pathStruct.expDirPath,'exp[0-9]');
expDetailsStruct.expDate = pathStruct.expDirPath(a:b+6); expDetailsStruct.expDateStr = pathStruct.expDirPath(a:d);
[a,~] = regexp(expDetailsStruct.expDate,'_');expDetailsStruct.expDate(a) = '-';


%subdirectories should include
%> [ flat mstack ]
%> [ mstack images ]
%> [ segment mstack ]
%> [ tracking files ]


%load helpful metadata
cd(dirStruct.exportdir)
FileName = expDetailsStruct.expDateStr;
datequery = strcat(FileName,'*DoseAndScene*');
cd(dirStruct.exportdir)
filelist = dir(datequery);
if isempty(filelist)
    error(strcat('need to run ExtractMetadata for-',FileName));
    %        dosestruct = makeDoseStruct; %run function to make doseStruct
else
    dosestructstruct = load(char(filelist.name));
    dosestruct = dosestructstruct.dosestruct;
end
segInstruct = dosestructstruct.segInstruct;

segStruct.nucleus_seg = segInstruct.nucleus;
% segStruct.cell_seg = segInstruct.cell;
segStruct.background_seg = segInstruct.background;
channelstoinput = dosestructstruct.channelNameSwapArray;
channelinputs =channelregexpmaker(channelstoinput);
bkg = dosestructstruct.BACKGROUND;
ImageDetails.ImgSize = dosestructstruct.dimensions;

BACKGROUND = bkg{1};
bkarray = cell(1,length(BACKGROUND));
for i = 1:length(BACKGROUND)
    bkstr = num2str(BACKGROUND(i));
    if length(bkstr)>1
        bkarray{i} = strcat('s',bkstr);
    else
        bkarray{i} = strcat('s0',bkstr);
    end
end
bkinputs =channelregexpmaker(bkarray);


%determine how many scenes are present
dirlist = dir(pathStruct.mstackPath);
[~,~,~,d] = regexp({dirlist.name},'s[0-9]++');
dlog = ~cellfun(@isempty,d,'UniformOutput',1);
dcell = d(dlog);
SceneList = unique(cellfun(@(x) x{1},dcell,'UniformOutput',0));


%remove scenes that are background images
[~,~,~,d] = regexp(SceneList,bkinputs);
bkgscenelog = cellfun(@isempty,d,'UniformOutput',1);
SceneList = SceneList(bkgscenelog);


%determine the number of frames per scene
cd(pathStruct.mstackPath)
dirlist = dir('*.mat');
filearray = {dirlist.name};
filename = filearray{1};
fileObject = matfile(filename);
dim = size(fileObject,'flatstack');
    if max(size(dim))>2
        timeFrames = dim(3);
    else
        timeFrames = 1;
    end
chanbkgmedianmat = zeros(1,timeFrames);
nucbkgmedianmat = zeros(1,timeFrames);
frameToLoad = 1;

xAxisLimits = [0 timeFrames];
%set colormap
cd(dirStruct.parentdir)
addpath('Colormaps')
cmaplength = min([max(timeFrames) 100]);
cmap = vertcat(gray(256),jet(cmaplength));

%% set up user interface
fnum=22;
mainfig = figure(fnum);
initialsize = [1 1 2560 1080];
mainfig.Units = 'pixels';
figdim = [initialsize(3) initialsize(4)];
fpos = [1 1 figdim(1) figdim(2)];
mainfig.Position = fpos;
fW = fpos(3);
fH = fpos(4);


bW = initialsize(3)./11; %button width
bH = initialsize(4)/40; %button height
yP = sort(110:bH:(fH-(fH/9)),'descend');
xP = ones(1,length(yP)).*(fW-(bW*2));
fontSize = 10;




%instruct user how they can change to view different channels
mmm=1;
uicontrol('Style','text','String','To choose channel push 1, 2, or 3',...
    'Position',[xP(mmm),yP(mmm)+bH./1.5,bW,bH*2]);


%set axis limits button
uicontrol('Style','pushbutton',...
    'String','plotAxisLimits [x]',...
    'Position',[fW-(bW*5),yP(mmm),bW,bH],...
    'Callback',@plotAxis_callback);

%move forward frame or previuos frame
mmm=2;
uicontrol('Style','pushbutton',...
    'String','NextFrame [f]',...
    'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
    'Callback',@nextbutton_callback);
uicontrol('Style','pushbutton',...
    'String','PreviousFrame [a]',...
    'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
    'Callback',@prevbutton_callback);
uicontrol('Style','pushbutton',...
    'String','JumpToFrame [j]',...
    'Position',[xP(mmm)-(bW./2)-bW,yP(mmm),bW,bH],...
    'Callback',@jumpToFrame);

%next row is final Frame and first Frame commands
mmm=3;
uicontrol('Style','pushbutton',...
    'String','FinalFrame [h]',...
    'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
    'Callback',@finalbutton_callback);
uicontrol('Style','pushbutton',...
    'String','FirstFrame [z]',...
    'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
    'Callback',@firstbutton_callback);
uicontrol('Style','pushbutton',...
    'String','StimulationFrame [g]',...
    'Position',[xP(mmm)-bW - bW./2,yP(mmm),bW,bH],...
    'Callback',@firstbutton_callback);

%choose scene text and popup menu
mmm=5;
uicontrol('Style','text','String','Choose Scene',...
    'Position',[xP(mmm),yP(mmm),bW,bH]);

uicontrol('Style','popupmenu',...
    'String',SceneList',...
    'Position',[xP(mmm),yP(mmm)-bH./2,bW,bH./1.5],...
    'Callback',@popup_menu_Callback);

%add area and remove area
mmm=7;
uicontrol('Style','pushbutton','String','AddArea [v]',...
    'Position',[xP(mmm)-bW,yP(mmm),bW,bH],...
    'BackgroundColor',[0.7 0.7 1],...
    'Callback',@addareabutton_Callback);
uicontrol('Style','pushbutton','String','ReplaceArea [y]',...
    'Position',[xP(mmm)+bW,yP(mmm),bW,bH],...
    'BackgroundColor',[0.7 0.7 1],...
    'Callback',@replaceareabutton_Callback);
uicontrol('Style', 'pushbutton', 'String', 'Remove area [r]',...
    'Position',[xP(mmm),yP(mmm),bW,bH],...
    'BackgroundColor',[0.7 0.7 1],...
    'Callback',@removeArea_Callback);

%delete commands
mmm=8;
uicontrol('Style','pushbutton','String','DeleteArea',...
    'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
    'BackgroundColor',[0.7 0.7 1],...
    'Callback',@deletebutton_Callback);
uicontrol('Style','pushbutton','String','DeleteAllOnFrame',...
    'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
    'BackgroundColor',[0.7 0.7 1],...
    'Callback',@deleteAllonFrame_Callback);

%erode or dilate
mmm=9;
uicontrol('Style','pushbutton','String','Erode Selected Nuclei',...
    'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
    'BackgroundColor',[0.7 0.7 1],...
    'Callback',@erodeOnes_Callback);
uicontrol('Style','pushbutton','String','Dilate Selected Nuclei',...
    'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
    'BackgroundColor',[0.7 0.7 1],...
    'Callback',@dilateOnes_Callback);



%destroy commands
mmm=11;
uicontrol('Style','pushbutton','String','Destroy [q]',...
    'Position',[xP(mmm),yP(mmm),bW,bH],...
    'Callback',@destroybutton_Callback);
hDestroy = uicontrol('Style','pushbutton','String','DestroyPrevious [w]',...
    'Position',[xP(mmm)-bW,yP(mmm),bW,bH],...
    'Callback',@destroybuttonPrevious_Callback);
hDestroy.FontSize=fontSize-2;
hDestroy = uicontrol('Style','pushbutton','String','DestroySubsequent [e]',...
    'Position',[xP(mmm)+bW,yP(mmm),bW,bH],...
    'Callback',@destroybuttonSubsequent_Callback);
hDestroy.FontSize=fontSize-2;

%destroy commands
mmm=12;

hDestroy = uicontrol('Style','pushbutton','String','DestroyAllFramePrevious',...
    'Position',[xP(mmm)-bW,yP(mmm),bW+bW./2,bH],...
    'Callback',@destroyAllFramePrevious_Callback);
hDestroy.FontSize=fontSize-2;
hDestroy = uicontrol('Style','pushbutton','String','DestroyAllFrameSubsequent',...
    'Position',[xP(mmm)+bW./2,yP(mmm),bW+bW./2,bH],...
    'Callback',@destroyAllFrameSubsequent_Callback);
hDestroy.FontSize=fontSize-2;


%chosen ones commands
mmm=13;
uicontrol('Style','pushbutton','String','Chosen Ones',...
    'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
    'Callback',@chosenOnes_Callback);
uicontrol('Style','pushbutton','String','Chosen OnesAllOnFrame',...
    'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
    'Callback',@chosenOnesAllOnFrame_Callback);



%tracking commands
mmm=15;
uicontrol('Style','pushbutton','String','Break link [d]',...
    'Position',[xP(mmm)-bW/2,yP(mmm),bW,bH],...
    'BackgroundColor',[0.7 1 0.7],...
    'Callback',@breaklink_callback);
uicontrol('Style','pushbutton','String','autoBreakEdge [9]',...
    'Position',[xP(mmm)+bW/2,yP(mmm),bW,bH],...
    'BackgroundColor',[0.7 1 0.1],...
    'Callback',@breakEdge_callback);

mmm=16;
uicontrol('Style','pushbutton','String','LinkCells [c]',...
    'Position',[xP(mmm),yP(mmm),bW,bH],...
    'BackgroundColor',[0.7 1 0.7],...
    'Callback',@linkCells_Callback);

mmm=17;
uicontrol('Style','pushbutton',...
    'String','Run Tracking [t]',...
    'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
    'Callback',@trackbutton_Callback);
uicontrol('Style','pushbutton',...
    'String','divisionTrack [y]',...
    'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
    'BackgroundColor',[1 0.7 0.7],...
    'Callback',@divisionTrack_Callback);

mmm=18;
uicontrol('Style','pushbutton',...
    'String','LoadTracking',...
    'Position',[xP(mmm)-bW/2,yP(mmm),bW,bH],...
    'Callback',@loadTrackingFile_callback);
uicontrol('Style','pushbutton',...
    'String','updateTracking',...
    'Position',[xP(mmm)+bW/2,yP(mmm),bW,bH],...
    'Callback',@updatetracking_callback);




%display commands
mmm=20;
uicontrol('Style','pushbutton',...
    'String','Set Contrast [k]',...
    'Position',[xP(mmm),yP(mmm),bW,bH],...
    'Callback',@contrast_Callback);
uicontrol('Style','pushbutton',...
    'String','togStruct.displayTrackingToggle [m]',...
    'Position',[xP(mmm),yP(mmm)-bH./2,bW,bH./1.5],...
    'Callback',@displayTrackingButton_Callback);


%save commands
mmm=22;
uicontrol('Style','pushbutton',...
    'String','SaveTrackingAs',...
    'Position',[xP(mmm)-bW./2,yP(mmm),bW+bW,bH],...
    'Callback',@saveTrackingFileAs_callback);

%autotrack commands
mmm=23;
uicontrol('Style','pushbutton',...
    'String','trackSaveIterate',...
    'Position',[xP(mmm)+bW./2,yP(mmm),bW+bW./2,bH./1.5],...
    'Callback',@trackSaveIterate_callback);
uicontrol('Style','pushbutton',...
    'String','TSIchosen',...
    'Position',[xP(mmm)-bW,yP(mmm),bW+bW./2,bH./1.5],...
    'Callback',@trackSaveIterateChosen_callback);


%plot and plot settings commands
mmm=25;
uicontrol('Style','pushbutton',...
    'String','PLOT!',...
    'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
    'Callback',@Plot_callback);
uicontrol('Style','pushbutton',...
    'String','Plot Specific Cell!',...
    'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
    'Callback',@PlotSpecificCell_callback);
mmm=26;
uicontrol('Style','pushbutton',...
    'String','Plot Settings!',...
    'Position',[xP(mmm)-bW./2,yP(mmm),bW+bW./2,bH./1.5],...
    'Callback',@PlotSettings_callback);

%export commands
mmm=28;
uicontrol('Style','pushbutton',...
    'String','ExportTrackedCells',...
    'Position',[xP(mmm),yP(mmm),bW,bH*2],...
    'Callback',@exportTrackedCells);
uicontrol('Style','pushbutton',...
    'String','ExportAllCells',...
    'Position',[xP(mmm)-bW,yP(mmm),bW,bH./1.5],...
    'Callback',@exportAllCells);
uicontrol('Style','pushbutton',...
    'String',{'ExportSegmentedCells'},...
    'Position',[xP(mmm)+bW,yP(mmm),bW,bH*2],...
    'Callback',@exportSegmentedCells);

%label and comment commands
mmm=30;
uicontrol('Style','pushbutton',...
    'String','Label Cells',...
    'Position',[xP(mmm)-bW,yP(mmm),bW,bH./1.5],...
    'Callback',@labelCells);
uicontrol('Style','pushbutton','String','Comments',...
    'Position',[xP(mmm),yP(mmm),bW,bH./1.5],...
    'Callback',@comment_Callback);
uicontrol('Style','pushbutton',...
    'String','ExportLabels',...
    'Position',[xP(mmm)+bW,yP(mmm),bW,bH./1.5],...
    'Callback',@exportLabels);

mmm=31;
uicontrol('Style','pushbutton',...
    'String','MakeMovie',...
    'Position',[xP(mmm)-bW/2,yP(mmm),bW,bH./1.5],...
    'Callback',@makeMovieCropped);

mmm=31;
uicontrol('Style','pushbutton',...
    'String','exportFrames',...
    'Position',[xP(mmm)+bW/2 + bW,yP(mmm),bW,bH./1.5],...
    'Callback',@saveFramesCropped);

uicontrol('Style','pushbutton',...
    'String','SaveMovie',...
    'Position',[xP(mmm)+bW/2,yP(mmm),bW,bH./1.5],...
    'Callback',@saveMovieCropped);

mmm=32;
uicontrol('Style','pushbutton',...
    'String','CellMovie [h]',...
    'Position',[xP(mmm)-bW/2,yP(mmm),bW,bH./1.5],...
    'Callback',@displaySmallMovieFunct);
uicontrol('Style','pushbutton',...
    'String','SaveCellMovie',...
    'Position',[xP(mmm)+bW/2,yP(mmm),bW,bH./1.5],...
    'Callback',@saveSmallMovieFunct);

%refine figure details
mainfig.Visible = 'on'   ;
mainfig.Units = 'normalized';
for i = 1:length(mainfig.Children)
    hhh = mainfig.Children(i);
    hhh.Units = 'normalized';
end

MainAxes = axes;
MainAxes.Units = 'pixels';
MainAxes.XTick=[];
MainAxes.YTick = [];
imgdim = (fW./2).*0.95;
Position = [10 10 imgdim imgdim];
% Position = [0.1 0.3 0.65 0.65];
MainAxes.Position = Position;
MainAxes.Units = 'normalized';
colormap(MainAxes,cmap);
c = colorbar(MainAxes);
c.Limits = [256 size(cmap,1)];
c.Location = 'eastoutside';
c.TickLabels = c.Ticks -255;
c.Label.String = 'track length, # of frames';


PlotAxes = axes;
PlotAxes.Units = 'pixels';
pos = [imgdim+10+200   0.6605    0.1500    0.1500];
PlotAxes.Position = pos;
PlotAxes.Units = 'normalized';
pos = PlotAxes.Position;
pos(2:4) = [0.6605    0.1500    0.1500];
PlotAxes.Position = pos;

SecondPlotAxes = axes;
SecondPlotAxes.Units = 'pixels';
pos = [imgdim+10+200   0.6605    0.1500    0.1500];
SecondPlotAxes.Position = pos;
SecondPlotAxes.Units = 'normalized';
pos = SecondPlotAxes.Position;
pos(2:4) = [0.4605    0.1500    0.1500];
SecondPlotAxes.Position = pos;

ThirdPlotAxes = axes;
ThirdPlotAxes.Units = 'pixels';
pos = [imgdim+10+200   0.6605    0.1500    0.1500];
ThirdPlotAxes.Position = pos;
ThirdPlotAxes.Units = 'normalized';
pos = ThirdPlotAxes.Position;
pos(2:4) = [0.2605    0.1500    0.1500];
ThirdPlotAxes.Position = pos;

mainfig.Units = 'pixels';
scrnsize = get(0,'screensize');
fpos = scrnsize;
fpos(1) = scrnsize(3)/2;
fpos(3) = scrnsize(3)/2;
fpos(4) = scrnsize(4).*0.75;
mainfig.Position = fpos;

%constrain image axis to be square initially
MainAxes.Units = 'pixels';
pos = MainAxes.Position;
pos(3:4) = [min([pos(3) pos(4)]) min([pos(3) pos(4)])];
MainAxes.Position = pos;
MainAxes.Units = 'normalized';
pos = MainAxes.Position;
pos(2)= 0.5 - pos(4)/2;
MainAxes.Position = pos;
imgsize =ImageDetails.ImgSize;
MainAxes.YAxis.Limits = [0 imgsize(1)];
MainAxes.XAxis.Limits = [0 imgsize(2)];

mainfig.Color = 'w';
set(mainfig,'KeyPressFcn',@keypress);


pStruct = loadSegmentParameters([],FileName,dirStruct.exportdir); %loads saved value of pStruct
timeMat = loadMetadata(FileName,dirStruct.exportdir);
rtimeMat = round(timeMat);
stimulationFrame = dosestruct(1).tgfFrame + 1;
timeVec = rtimeMat(1,:)-rtimeMat(1,stimulationFrame);
timeSteps = diff(round(timeVec),[],2);


stackstruct = struct(); %place up top and make global
foStruct = struct();
for ic = 1:length(channelstoinput)
    cstr = channelstoinput{ic};
    %choose the channel image
    imgstack = zeros([imgsize([1 2]) timeFrames]);
    
    stackstruct.(cstr) = imgstack;
    foStruct.(cstr) = [];
end

segStruct.nucleus_seg = segInstruct.nucleus;
% segStruct.cell_seg = segInstruct.cell;
segStruct.background_seg = segInstruct.background;
segfnames = fieldnames(segStruct);
for ic = 1:length(segfnames)
    cstr = segfnames{ic};
    %choose the channel image
    imgstack = zeros([imgsize([1 2]) timeFrames]);
    
    stackstruct.(cstr) = logical(imgstack);
    foStruct.(cstr) = [];
end


if strcmpi(AutoTrackStr,'AutoTrackCells')
    mainfig.Visible = 'off';
    trackSaveIterateChosen_callback([],[])
    % exportNuclei([],[])
    close(mainfig)
elseif strcmpi(AutoTrackStr,'AutoExportNuclei')
    exportSegmentedCells([],[])
    close(mainfig)
elseif strcmpi(AutoTrackStr,'AutoExportTracks')
    mainfig.Visible = 'off';
    exportTrackedCells([],[])
    close(mainfig)
else
    Tracked = makeTrackingFile(timeFrames);
    ImageDetails.Scene = SceneList{1};
    ImageDetails.Channel = segStruct.nucleus_seg;
    ImageDetails.Frame = 1;
    ImageDetails.ImgSize = ImageDetails.ImgSize;
    setSceneAndTime
end


end

function timeVec = loadMetadata(FileName,loaddir)

%load metadata associated with the experiment (requires manual input if there is ambiguity)
datequery = strcat(FileName,'*metaData.mat');
cd(loaddir)
filelist = dir(datequery);
if length({filelist.name}) ==1
    metaData = load(char(filelist.name));
else
    filename = uigetfile();
    metaData = load(filename);
end

%determine the timeVector from metaData [dim 1 is scene#, dim 2 is each time frame]
timeVec = metaData.timeVec;
timeVec = round(timeVec);
end
function pStruct = loadSegmentParameters(pStruct,datename,loaddir)

cd(loaddir)
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


%% uifunctions
function keypress(fig_obj,~)
global  ImageDetails segStruct channelstoinput
key = get(fig_obj,'CurrentKey');

switch key
    case '~'
        
    case '1'
        ImageDetails.Channel = segStruct.background_seg;
        setSceneAndTime
    case '2'
        ImageDetails.Channel = segStruct.nucleus_seg;
        setSceneAndTime
    case '3'
        ImageDetails.Channel = 'DIC';
        setSceneAndTime
    case '4'
        ImageDetails.Channel = 'BKGbinary';
        setSceneAndTime
    case '5'
        ImageDetails.Channel = 'overlay';
        setSceneAndTime
    case '6'
        ImageDetails.Channel = 'FluorOnlyOverlay';
        setSceneAndTime
    case 'hyphen'
        ci = strcmp(ImageDetails.Channel,channelstoinput);
        if sum(ci)>0
            cidx = find(ci==1)-1;
        else
            cidx = 1;
        end
        if cidx<1
            cidx = 1;
        end
        cstr = channelstoinput{cidx};
        ImageDetails.Channel = cstr;
        setSceneAndTime
    case 'equal'
        ci = strcmp(ImageDetails.Channel,channelstoinput);
        if sum(ci)>0
            cidx = find(ci==1)+1;
        else
            cidx = 1;
        end
        if cidx>length(channelstoinput)
            cidx = 1;
        end
        cstr = channelstoinput{cidx};
        ImageDetails.Channel = cstr;
        setSceneAndTime
    case 'q'
        destroybutton_Callback([],[])
    case 'w'
        destroybuttonPrevious_Callback([],[])
    case 'e' 
        destroybuttonSubsequent_Callback([],[])
    case 'a'
        prevbutton_callback([],[])
    case 'f'
        nextbutton_callback([],[])
    case 'leftbracket'
        prevScene_callback([],[])
    case 'rightbracket'
        nextScene_callback([],[])
        
    case 'd'
        breaklink_callback([],[]);
    case '9'
        breakEdge_callback([],[]);
    case 'h'
        displaySmallMovieFunct;
    case '8'
        togglemmstr
    case '7'
        offplot
    case 't'
        trackbutton_Callback([],[]);
    case 'v'
        addareabutton_Callback([],[]);
    case 'r'
        removeArea_Callback([],[]);
    case 'y'
        replaceareabutton_Callback([],[]);
    case 'c'
        linkCells_Callback([],[]);
    case 'm'
        displayTrackingButton_Callback([],[])
    case 'j'
        jumpToFrame([],[])
    case 'g'
        stimulationFramebutton_callback([],[])
    case 'z'
        firstbutton_callback([],[])
    case 's'
        saveTrackingFileAs_callback([],[])
    case 'l'
        loadTrackingFile_callback([],[])
    case 'p'
        Plot_callback([],[]);
    case 'comma'
        makeFCmovie([],[])
    case 'o'
        PlotSpecificCell_callback([],[])
    case 'i'
        PlotSpecificCellIteratively_callback([],[])
    case '0'
        makeMovieAllChannels([],[])
        
%     case 'o'
%         labelCells;
%     case 'u'
%         %         if displaycomments==1
%         %             displaycomments=0;
%         %         else
%         displaycomments=1;
%         [~] = getxy([],[]);
%         %         end
    case 'k'
        contrast_Callback([],[])
%     case 'k'
%         comment_Callback([],[])
%     case 'j'
%         comment_CallbackJ([],[])
%     case 'n'
%         PlotCFPnorm_callback([],[])
%     case 'b'
%         PlotCFPnotnorm_callback([],[])
    case 'x'
        plotAxis_callback([],[])
%     case '0'
%         displaycomments=1;
%         xy = getxy([],[]);
%         [~,comments,commentpos,~]=updatecomments(xy);
%         setcommentsTracking(comments,commentpos)
%         dispxy(xy)
end

end

%% division functions
function togglemmstr
global medmeantot

if strcmp(medmeantot,'median')
    medmeantot = 'mean';
elseif strcmp(medmeantot,'mean')
    medmeantot = 'total';
elseif strcmp(medmeantot,'total')
    medmeantot = 'nuccyto';
elseif strcmp(medmeantot,'nuccyto')
    medmeantot = 'median';
end
Plot_callback([],[])
end
function offplot
global mainfig
h = findobj(mainfig,'Type','Axes');
o = findobj(h,'Type','Line');
o.delete;
end

function divisionTrack_Callback(~,~)
global ImageDetails Tracked MainAxes DivisionStruct

t = ImageDetails.Frame;
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;


% Display mesh plot of the currently selected data.

MainAxes.Title.String = 'click mother cell...';
[cellxx,cellyy] = ginput(1);

cellx = round(cellxx);
celly = round(cellyy);

mothercellind = sub2ind(ImageDetails.ImgSize,celly,cellx);

for j = 1:length(cellxx)
    if j==1
        idxs = cellfun(@(x) isempty(find(x==mothercellind(j),1)),PX,'UniformOutput',1);
    else
        idxs = idxs & cellfun(@(x) isempty(find(x==mothercellind(j),1)),PX,'UniformOutput',1);
    end
end



nextbutton_callback([],[]);


MainAxes.Title.String = 'click daughter cells...';
[cellxx,cellyy] = ginput(2);

cellx = round(cellxx);
celly = round(cellyy);

daughtercellinds = sub2ind(ImageDetails.ImgSize,celly,cellx);

for j = 1:length(cellxx)
    if j==1
        idxs = cellfun(@(x) isempty(find(x==daughtercellinds(j),1)),PX,'UniformOutput',1);
    else
        idxs = idxs & cellfun(@(x) isempty(find(x==daughtercellinds(j),1)),PX,'UniformOutput',1);
    end
end



divStruct.motherLocation = mothercellind;
divStruct.frame = t;
divStruct.daughterOneLocation = daughtercellinds(1);
divStruct.daughterTwoLocation = daughtercellinds(2);

ds = checkDivision(divStruct);
DivisionStruct = ds;


setSceneAndTime;
end
function ds = checkDivision(divStruct)
global Tracked DivisionStruct

%open the pixels from the scene of the recently selected dividing cell
t = divStruct.frame;
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;

%determine the cellID of the recently selected dividing cell
mothercellind = divStruct.motherLocation;
idxs = cellfun(@(x) ~isempty(find(x==mothercellind,1)),PX,'UniformOutput',1);
%idxs is the recently divided cell

olddivStruct = DivisionStruct;
if ~isempty(olddivStruct)&&isstruct(olddivStruct)
    %find repeat cells
    fnames = fieldnames(olddivStruct)';
    cycle=0;
    newdivStruct = struct();
    for i = 1:length(olddivStruct)
        cycle=cycle+1;
        tdiv = olddivStruct(i).frame;
        mothercellind = olddivStruct(i).motherLocation;
        CC = Tracked{tdiv}.Cellz;
        PX = CC.PixelIdxList;
        idxold = cellfun(@(x) ~isempty(find(x==mothercellind,1)),PX,'UniformOutput',1);
        %if the cells do match, then replace the old cell (that is, do not keep the old cell)
        iold = find(idxold==1);
        inew = find(idxs==1);
        if ~(iold == inew) % if the cells don't match then keep the old cell
            for fname = fnames
                newdivStruct(cycle).(char(fname)) = olddivStruct(i).(char(fname));
            end
        end
    end
    
    %add the old cell
    lnds = length(newdivStruct)+1;
    if lnds ==0
        lnds=1;
    end
    for fname = fnames
        newdivStruct(lnds).(char(fname)) = divStruct.(char(fname));
    end
    
    divisionStructure = newdivStruct;
else
    divisionStructure=struct();
    divisionStructure(1).motherLocation = divStruct.motherLocation;
    divisionStructure(1).frame = t;
    divisionStructure(1).daughterOneLocation = divStruct.daughterOneLocation;
    divisionStructure(1).daughterTwoLocation = divStruct.daughterTwoLocation;
end

ds = divisionStructure;
end

%choose scenes
function nextScene_callback(~,~)
global ImageDetails Tracked timeFrames SceneList

if isempty(ImageDetails.Scene)
    ImageDetails.Scene = SceneList{1};
    
else
    sceneval = ImageDetails.Scene;
    [~,~,~,d] = regexp(SceneList,sceneval);
    didx = ~cellfun(@isempty,d,'UniformOutput',true);
    dval = find(didx==1)+1;
    dval2 = min([dval length(SceneList)]);
    ImageDetails.Scene = SceneList{dval2};
    

end

Trackedz = makeTrackingFile(timeFrames);
Tracked=Trackedz;
setSceneAndTime

end
function prevScene_callback(~,~)
global ImageDetails Tracked timeFrames SceneList

if isempty(ImageDetails.Scene)
    ImageDetails.Scene = SceneList{1};
    
else
    sceneval = ImageDetails.Scene;
    [~,~,~,d] = regexp(SceneList,sceneval);
    didx = ~cellfun(@isempty,d,'UniformOutput',true);
    dval = find(didx==1)-1;
    dval2 = max([dval 1]);
    ImageDetails.Scene = SceneList{dval2};
    

end

Trackedz = makeTrackingFile(timeFrames);
Tracked=Trackedz;
setSceneAndTime

end

%choose frames
function nextbutton_callback(~,~) %next buttong
global frameToLoad ImageDetails timeFrames

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = frameToLoad;
end

frameToLoad = ImageDetails.Frame + 1;
if frameToLoad>timeFrames
    frameToLoad = timeFrames;
end

ImageDetails.Frame = frameToLoad;
setSceneAndTime
end
function prevbutton_callback(~,~) %previous button
global frameToLoad ImageDetails

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = frameToLoad;
end

frameToLoad = ImageDetails.Frame - 1;
if frameToLoad<1
    frameToLoad = 1;
end
ImageDetails.Frame = frameToLoad;
setSceneAndTime
end
function finalbutton_callback(~,~) %final button
global frameToLoad ImageDetails timeFrames

frameToLoad = timeFrames;
ImageDetails.Frame = frameToLoad;
setSceneAndTime
end
function stimulationFramebutton_callback(~,~)
global frameToLoad ImageDetails timeFrames stimulationFrame

frameToLoad = timeFrames;
ImageDetails.Frame = stimulationFrame;
setSceneAndTime
end
function firstbutton_callback(~,~) %first button
global frameToLoad ImageDetails

frameToLoad = 1;
ImageDetails.Frame = frameToLoad;
setSceneAndTime
end
function jumpToFrame(~,~)
global frameToLoad ImageDetails

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = frameToLoad;
end

prompt = {'Go to which frame'};
dlg_title = 'Go to frame...';
idx = str2double(cell2mat(inputdlg(prompt,dlg_title)));

frameToLoad = idx;
ImageDetails.Frame = frameToLoad;
setSceneAndTime
end


function popup_menu_Callback(source,~)
global ImageDetails Tracked timeFrames


Trackedz = makeTrackingFile(timeFrames);
Tracked=Trackedz; 

%Determine the selected data set.
str = source.String;
val = source.Value;
pvalue = char(str{val});

ImageDetails.Scene = pvalue;
setSceneAndTime

end

%% these functions alter segmented cell areas/pixels
%add cells and link cells
function addareabutton_Callback(~,~)
global  ImageDetails Tracked 

ArrayStruct = Tracked.arrayStruct;
trackmatrix = Tracked.trackmatrix;
segment_Pixels_array = ArrayStruct.pixels;
PXarray = [];
tvec = [];
dropArray = [];
newArray = [];
keepArray=[];
allArray=[];

button=1;
while button == 1
    [polyx,polyy,button] = ginput();
    button = round(mean(button));
    
    if button ==1 && length(polyx)>2
        M = zeros(1,length(polyx)*2);
        M(1:2:end) = polyx;
        M(2:2:end) = polyy;
        zeroImage = zeros(ImageDetails.ImgSize);
        zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
        zerogray = rgb2gray(zeroImage);
        
        if ~isempty(segment_Pixels_array)
            t = ImageDetails.Frame;
            %determine set of pixels
            px = find((zerogray>0)==1);
            PX = segment_Pixels_array{t};
            %find pixel overlaps
            idxs = cellfun(@(x) sum(ismember(x,px)),PX,'UniformOutput',1);
            index = find(idxs>1);
            
            if ~isempty(index)
                %find location in trackmatrix [trackmatrix(frame,tracknum)]
                alltrackframe = trackmatrix(t,:);
                trackidx = zeros(size(index));
                tracklength = zeros(size(index));
                for j = 1:length(index)
                    indexnum = index(j);
                    ti = find(alltrackframe == indexnum);
                    if ~isempty(ti)
                        trackidx(j) = ti;
                        tracklength(j) = sum(~isnan(trackmatrix(:,j)));
                    end
                end
                
                [~,maxtrackidxnum] = max(tracklength);
                maxtrackidx = false(size(index));
                maxtrackidx(maxtrackidxnum) = true;
                
                pxGroup = vertcat(PX{index}); %join pixel values from all overlapping cells
                pxnew = {unique(vertcat(pxGroup,px))}; %combine all pixel values into set
                
                %define keepidx,dropidx,newidx
                keepidx = index(maxtrackidx);
                dropidx = index(~maxtrackidx);
                newidx = [];
                allidx = keepidx;
            else
                pxnew = {px};
                
                %define keepidx,dropidx,newidx
                keepidx = [];
                dropidx = [];
                newidx = length(PX)+1;
                allidx = newidx;
            end
            
            trackmatrix = updateTrackMatrix(trackmatrix,t,dropidx,keepidx,newidx);
            dropArray = horzcat(dropArray,{dropidx});
            keepArray = horzcat(keepArray,{keepidx});
            newArray = horzcat(newArray,{newidx});
            allArray = horzcat(allArray,{allidx});
            PXarray = horzcat(PXarray,{pxnew});
            tvec = horzcat(tvec,t);
        end
        nextbutton_callback([],[]);
    else
        if isempty(tvec)
            return
        else
            break
        end
    end
end
    %array struct does not work here because it runs twice!
ArrayStruct = updateArrayStruct(ArrayStruct,PXarray,tvec,dropArray,newArray,keepArray,allArray,trackmatrix,ImageDetails.ImgSize); %this is not ideal for updating trackmatrix or removing
Tracked.arrayStruct = ArrayStruct;
Tracked.trackmatrix = trackmatrix;
ImageDetails.Frame = tvec(1);
ucell = 1:size(trackmatrix,2);
trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
setSceneAndTime;
end
%add cells and link cells
function replaceareabutton_Callback(~,~)
global  ImageDetails Tracked 

ArrayStruct = Tracked.arrayStruct;
trackmatrix = Tracked.trackmatrix;
segment_Pixels_array = ArrayStruct.pixels;
PXarray = [];
tvec = [];
dropArray = [];
newArray = [];
keepArray=[];
allArray=[];

button=1;
while button == 1
    [polyx,polyy,button] = ginput();
    button = round(mean(button));
    
    if button ==1 && length(polyx)>2
        M = zeros(1,length(polyx)*2);
        M(1:2:end) = polyx;
        M(2:2:end) = polyy;
        zeroImage = zeros(ImageDetails.ImgSize);
        zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
        zerogray = rgb2gray(zeroImage);
        
        if ~isempty(segment_Pixels_array)
            t = ImageDetails.Frame;
            %determine set of pixels
            px = find((zerogray>0)==1);
            PX = segment_Pixels_array{t};
            %find pixel overlaps
            idxs = cellfun(@(x) sum(ismember(x,px)),PX,'UniformOutput',1);
            index = find(idxs>1);
            
            if ~isempty(index)
                %find location in trackmatrix [trackmatrix(frame,tracknum)]
                alltrackframe = trackmatrix(t,:);
                trackidx = zeros(size(index));
                tracklength = zeros(size(index));
                for j = 1:length(index)
                    indexnum = index(j);
                    ti = find(alltrackframe == indexnum);
                    if ~isempty(ti)
                        trackidx(j) = ti;
                        tracklength(j) = sum(~isnan(trackmatrix(:,j)));
                    end
                end
                
                [~,maxtrackidxnum] = max(tracklength);
                maxtrackidx = false(size(index));
                maxtrackidx(maxtrackidxnum) = true;
                
%                 pxGroup = vertcat(PX{index}); %join pixel values from all overlapping cells
%                 pxnew = {unique(vertcat(pxGroup,px))}; %combine all pixel values into set
                pxnew = {px}; %combine all pixel values into set
                
                %define keepidx,dropidx,newidx
                keepidx = index(maxtrackidx);
                dropidx = index(~maxtrackidx);
                newidx = [];
                allidx = keepidx;
            else
                pxnew = {px};
                
                %define keepidx,dropidx,newidx
                keepidx = [];
                dropidx = [];
                newidx = length(PX)+1;
                allidx = newidx;
            end
            
            trackmatrix = updateTrackMatrix(trackmatrix,t,dropidx,keepidx,newidx);
            dropArray = horzcat(dropArray,{dropidx});
            keepArray = horzcat(keepArray,{keepidx});
            newArray = horzcat(newArray,{newidx});
            allArray = horzcat(allArray,{allidx});
            PXarray = horzcat(PXarray,{pxnew});
            tvec = horzcat(tvec,t);
        end
        nextbutton_callback([],[]);
    else
        if isempty(tvec)
            return
        else
            break
        end
    end
end
    %array struct does not work here because it runs twice!
ArrayStruct = updateArrayStruct(ArrayStruct,PXarray,tvec,dropArray,newArray,keepArray,allArray,trackmatrix,ImageDetails.ImgSize); %this is not ideal for updating trackmatrix or removing
Tracked.arrayStruct = ArrayStruct;
Tracked.trackmatrix = trackmatrix;
ImageDetails.Frame = tvec(1);
ucell = 1:size(trackmatrix,2);
trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
setSceneAndTime;
end
%removeArea
function removeArea_Callback(~,~)
global  ImageDetails Tracked 

ArrayStruct = Tracked.arrayStruct;
trackmatrix = Tracked.trackmatrix;
segment_Pixels_array = ArrayStruct.pixels;
PXarray = [];
tvec = [];
dropArray = [];
newArray = [];
keepArray=[];
allArray =[];

button=1;
while button==1
    [polyx,polyy,button] = ginput();
    button = round(mean(button));
    
    if button ==1 && length(polyx)>2
        M = zeros(1,length(polyx)*2);
        M(1:2:end) = polyx;
        M(2:2:end) = polyy;
        zeroImage = zeros(ImageDetails.ImgSize);
        zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
        zerogray = rgb2gray(zeroImage);
        
        if ~isempty(segment_Pixels_array)
            t = ImageDetails.Frame;
            %identify pixels
            px = find((zerogray>0)==1);
            PX = segment_Pixels_array{t};
            %find overlap
            idxs = cellfun(@(x) sum(ismember(x,px)),PX,'UniformOutput',1);
            index = find(idxs>1);
            if ~isempty(index)
                %the possibilites are: empty, all cells shaved only, all cells bisected only, some bisected and some shaved (which can lead to same number but different cells)
                
                %define oldPX
                oldPX = PX(index);
                %define newPX
                oldMass = vertcat(oldPX{:});
                overlap = ismember(oldMass,px);
                oldMass(overlap)=[];
                imagio = zeros(ImageDetails.ImgSize);
                imagio(oldMass)=1;
                cc = bwconncomp(imagio);
                newPX = cc.PixelIdxList;
                %                 S = regionprops(cc,'Centroid');
                %                 centroidMat = vertcat(S.Centroid);
%                 subindeces = nan(1,length(newPX));
                subindeces = [];
                for k = 1:length(newPX)
                    px = newPX{k};
                    %find overlap
                    idxs = cellfun(@(x) sum(ismember(x,px)),oldPX,'UniformOutput',1);
                    idxsnum = find(idxs>1);
                    if ~isempty(idxsnum)
                        subindeces = [subindeces(:)' idxsnum(:)'];
                    end
                end
                [nn, histidx] = histc(subindeces(:)',1:length(index));
                samecells = nn==1;
                newcells = nn>1;
                gonecells = nn<1;
                
                %define dropidx,keepidx,newidx
                keepidx = index(samecells);
                searchidx = index(subindeces);
                newnum = find(newcells==1);
                newspotsidx = [];
                searchupdateidx = ismember(histidx,newnum);
                for k = 1:length(newnum)
                    newn = searchidx(find(newnum(k) == histidx));
                    newspotsidx = horzcat(newspotsidx,newn);
                end
                newidx = length(PX) + (1:length(newspotsidx));
                dropidx = horzcat(index(gonecells),unique(newspotsidx)); %good
%                 allidx = [newidx(:)' keepidx(:)'];
                searchidx(searchupdateidx) = newidx;
                allidx = searchidx;
                pxnew = newPX; %good
            else
                dropidx = [];
                newidx = [];
                keepidx = [];
                allidx = [];
            end
            
            trackmatrix = updateTrackMatrix(trackmatrix,t,dropidx,keepidx,newidx);
            dropArray = horzcat(dropArray,{dropidx});
            keepArray = horzcat(keepArray,{keepidx});
            newArray = horzcat(newArray,{newidx});
            allArray = horzcat(allArray,{allidx});
            PXarray = horzcat(PXarray,{pxnew});
            tvec = horzcat(tvec,t);
        end
        nextbutton_callback([],[]);
    else
        if isempty(tvec)
            return
        else
            break
        end
    end
end
ArrayStruct = updateArrayStruct(ArrayStruct,PXarray,tvec,dropArray,newArray,keepArray,allArray,trackmatrix,ImageDetails.ImgSize); %this is not ideal for updating trackmatrix or removing
Tracked.arrayStruct = ArrayStruct;
ucell = 1:size(trackmatrix,2);
trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
Tracked.trackmatrix = trackmatrix;
setSceneAndTime
end
%delete cells
function deletebutton_Callback(~,~)
global ImageDetails  Tracked 

ArrayStruct = Tracked.arrayStruct;
trackmatrix = Tracked.trackmatrix;
t = ImageDetails.Frame;
segment_Pixels_array = ArrayStruct.pixels;


% Display mesh plot of the currently selected data.
[cellxx,cellyy] = ginput();
cellx = round(cellxx);
celly = round(cellyy);
cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);


PX = segment_Pixels_array{t};
for j = 1:length(cellxx)
    if j==1
        idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
    else
        idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
    end
end
index = find((~idxs) == 1);
dropidx = index;
dropArray = {dropidx};
newidx = [];
newArray = {newidx};
keepidx = [];
keepArray = {keepidx};
allArray = {index};
PXarray = {[]};
tvec = t;
trackmatrix = updateTrackMatrix(trackmatrix,t,dropidx,keepidx,newidx);
ArrayStruct = updateArrayStruct(ArrayStruct,PXarray,tvec,dropArray,newArray,keepArray,allArray,trackmatrix,ImageDetails.ImgSize); %this is not ideal for updating trackmatrix or removing
ucell = 1:size(trackmatrix,2);
trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
Tracked.arrayStruct = ArrayStruct;
Tracked.trackmatrix = trackmatrix;
setSceneAndTime;
end
function deleteAllonFrame_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked
ArrayStruct = Tracked.arrayStruct;
trackmatrix = Tracked.trackmatrix;
t = ImageDetails.Frame;
segment_Pixels_array = ArrayStruct.pixels;

PX = segment_Pixels_array{t};
idxs = true(size(PX));
index = find(idxs == 1);
dropidx = index;
dropArray = {dropidx};
newidx = [];
newArray = {newidx};
keepidx = [];
keepArray = {keepidx};
allArray = {index};
PXarray = {[]};
tvec = t;
trackmatrix = updateTrackMatrix(trackmatrix,t,dropidx,keepidx,newidx);
ArrayStruct = updateArrayStruct(ArrayStruct,PXarray,tvec,dropArray,newArray,keepArray,allArray,trackmatrix,ImageDetails.ImgSize); %this is not ideal for updating trackmatrix or removing
ucell = 1:size(trackmatrix,2);
trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
Tracked.arrayStruct = ArrayStruct;
Tracked.trackmatrix = trackmatrix;
setSceneAndTime;
end


function [cellxx,cellyy,timingkeeper,button] = identifyCellbyClick(t,ginputnum)
cellxx  =    [];
cellyy  =    [];
timingkeeper = [];
if ~isempty(ginputnum)
    [xx,yy,button] = ginputmod(ginputnum); %record each click
else
    [xx,yy,button] = ginput(); %record each click
end
if ~isempty(button)
    cellxx  =    xx;
    cellyy  =    yy;
    timingkeeper = t;
else
    return
end
end
function [idx,idxlog,notACell] = findCellByCoordinate(Tracked,ImgSize,cxx,cyy,tkeeper)
ArrayStruct = Tracked.arrayStruct;
segment_Pixels_array = ArrayStruct.pixels;
%convert coordinates to indices
cellx = round(cxx);
celly = round(cyy);
cellind = sub2ind(ImgSize,celly,cellx);

%find the cell
PX = segment_Pixels_array{tkeeper};
idx = nan(1,length(cellx));
for i =1:length(cellind)
    cellindi = cellind(i);
    idxlog = ~cellfun(@(x) sum(x==cellindi)>0,PX,'UniformOutput',1);
    if sum(~idxlog)>0
        idx(i) = find(idxlog==0);
    end
end
idx(isnan(idx))=[];
notACell = isempty(idx);
% if ~notACell && ~isempty(idx)
%     plotTheChosenCells(idx);
% end
end
function linkCells_Callback(~,~)
%fast link!
global ImageDetails Tracked  timeFrames


trackmatrix = Tracked.trackmatrix;

%set accumulating variables equal to empty
cellxx = [];
cellyy = [];
cellnn = [];
idxnn = [];
timingkeeper=[];


%choose the initial cell
t = ImageDetails.Frame;
notACell = true;
while notACell %keep running until a cell is clicked
    ginputnum=1;
    [cxx,cyy,tkeeper,button] = identifyCellbyClick(t,ginputnum);
    if isempty(button)
        return
    end
    %find cell
    [idx,~,notACell] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,tkeeper);
    if ~notACell && ~isempty(idx)
        plotTheChosenCells(idx);
    end
end
cellxx =[cellxx(:)' cxx];
cellyy =[cellyy(:)' cyy];
timingkeeper =[timingkeeper(:)' tkeeper];
%update cellnum
alltrackframe = trackmatrix(t,:);
cellnum = find(alltrackframe == idx);
if isempty(cellnum)
    newtrack = nan(size(trackmatrix,1),1);
    trackmatrix = horzcat(trackmatrix,newtrack);
    cellnum = size(trackmatrix,2);
    trackmatrix(t,cellnum) = idx;
end
cellnn = [cellnn(:)' cellnum];
idxnn = [idxnn(:)' idx];

a=1;
while a ==1
    broken = false;
    if button == 1
        for t = tkeeper:timeFrames
            ImageDetails.Frame = t;
            setSceneAndTime
%             pause(0.005)
            nanidx = isnan(trackmatrix(t,cellnum));
            if nanidx
                broken = true;
                break
            end
        end
    elseif button ==3
        t=t+1;
        broken = true;
        ImageDetails.Frame = t;
        setSceneAndTime
    end
    
    if broken %==true
        notACell = true;
        while notACell %keep running until a cell is clicked
            ginputnum=1;
            [cxx,cyy,tkeeper,button] = identifyCellbyClick(t,ginputnum);
            if isempty(button)
                break
            end
            %find cell
            [idx,~,notACell] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,tkeeper);
            if ~notACell && ~isempty(idx)
                plotTheChosenCells(idx);
            end
        end
        if isempty(button)
            break
        end
        cellxx =[cellxx(:)' cxx];
        cellyy =[cellyy(:)' cyy];
        timingkeeper =[timingkeeper(:)' tkeeper];
        %update cellnum
        alltrackframe = trackmatrix(t,:);
        cellnum = find(alltrackframe == idx);
        if isempty(cellnum)
            newtrack = nan(size(trackmatrix,1),1);
            trackmatrix = horzcat(trackmatrix,newtrack);
            cellnum = size(trackmatrix,2);
            trackmatrix(t,cellnum) = idx;
        end
        cellnn = [cellnn(:)' cellnum];
        idxnn = [idxnn(:)' idx];
    end
    
    if ~broken
        break
    end
end

timestart = horzcat(1,timingkeeper(2:end));
timeend = horzcat(timingkeeper(2:end)-1,timeFrames);
% disp(horzcat(timestart',timeend'))
cellnnz = horzcat(cellnn,cellnn(end));
cellLoc = cellnnz(1);
newtrackmatrix = trackmatrix;
for i = 1:length(timestart)
    trange = timestart(i):timeend(i);
    newtrackmatrix(trange,cellLoc) = trackmatrix(trange,cellnnz(i));
    if ~(cellnnz(i)==cellLoc)
        newtrackmatrix(trange,cellnnz(i)) = NaN;
    end
end

% ucell = unique(cellnnz);
ucell = 1:size(newtrackmatrix,2);
trackmatrix = refineTrackingMatrix(newtrackmatrix,ucell);
ImageDetails.Frame = max([1 ImageDetails.Frame-1]);
Tracked.trackmatrix = trackmatrix;
setSceneAndTime
alltrackframe = trackmatrix(timingkeeper(1),:);
idx = idxnn(1);
cellnum = find(alltrackframe == idx);
alltrackframe = trackmatrix(ImageDetails.Frame,:);
idx = alltrackframe(cellnum);
    if  ~isempty(idx)
        plotTheChosenCells(idx);
    end
end

function trackmatrix = refineTrackingMatrix(newtrackmatrix,ucell)
global togStruct
tracklengthcuttoff = 3;
% disp('correcting tracking...')
%now move old tracks to end
%(run refine tracking)
newtracker = [];
%remove empty tracks
logmat = newtrackmatrix>0;
dlogmat = diff(logmat);
for j = 1:length(ucell)
    cellnum = ucell(j);
    idxo = logmat(:,cellnum);
    didxo = dlogmat(:,cellnum);
    endoftrack = find(didxo(:)' == -1);
    beginoftrack = find(didxo(:)' == 1)+1;
    if idxo(1)==1
        beginoftrack = [1 beginoftrack(:)'];
    end
    
    %add last frame if track exists at last frame
    if idxo(length(idxo)) == 1
        endoftrack = [endoftrack(:)' length(idxo)];
    end
    
    tracks = zeros(length(beginoftrack),2);
    for i =1:length(beginoftrack)
        tracks(i,:)  = [beginoftrack(i) endoftrack(i)];
    end
    
    
    tracklengths = tracks(:,2) - tracks(:,1) + 1;
    for trackidx = 1:size(tracks,1) %start from the second track
        
        beginoftrack = tracks(trackidx,1);
        endoftrack = tracks(trackidx,2);
        
        newtrack = nan(size(idxo));
        %if tracklength is too short, then remove it
        if tracklengths(trackidx)<tracklengthcuttoff
            newtrackmatrix(beginoftrack:endoftrack,cellnum) = NaN; %and remove
        elseif trackidx>1 %elseif track length is long enough and there exists a subsequent track
            %then move the subsequent track to the end of the structure
            %make the same changes to trt
            newtrack(beginoftrack:endoftrack) = newtrackmatrix(beginoftrack:endoftrack,cellnum); %move to end
            newtrackmatrix(beginoftrack:endoftrack,cellnum) = NaN; %and remove
            newtracker = horzcat(newtracker,newtrack);
        else %if track length is long enough and the only track,
            %then do nothing
        end
    end
end
newtrackmatrix = horzcat(newtrackmatrix,newtracker);


%remove empty tracks
logmat = newtrackmatrix>0;
startlogmat = diff(logmat);
startlogmat(1,logmat(1,:)==1)=1;
startvec = sum(startlogmat==1,1);
endlogmat = diff(logmat);
endlogmat(end,logmat(end,:)==1)=-1;
endvec = sum(endlogmat==-1,1);
emptylog = (startvec==0)|(endvec==0);

if sum(emptylog)>0
    newtrackmatrix(:,emptylog)=[];
end
trackmatrix = newtrackmatrix;
togStruct.trackUpdated=true;
end

%destroy and chosen ones
function destroyAllFramePrevious_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked

%   determine the frame to load
t = ImageDetails.Frame;
trackmatrix = Tracked.trackmatrix;
trackmatrix(1:t-1,:) = NaN;
ucell = 1:size(trackmatrix,2);
trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
Tracked.trackmatrix = trackmatrix;
setSceneAndTime;

end
function destroyAllFrameSubsequent_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked

%   determine the frame to load
t = ImageDetails.Frame;
trackmatrix = Tracked.trackmatrix;
trackmatrix(t+1:end,:) = NaN;
ucell = 1:size(trackmatrix,2);
trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
Tracked.trackmatrix = trackmatrix;
setSceneAndTime;
end
function destroybuttonPrevious_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked
trackmatrix = Tracked.trackmatrix;
t = ImageDetails.Frame;
ginputnum=[];
[cxx,cyy,~,button] = identifyCellbyClick(t,ginputnum);
if isempty(button)
    return
end
%find cell
[idx,~,~] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,t);
alltrackframe = trackmatrix(t,:);
cellnum = find(ismember(alltrackframe,idx));
if ~isempty(cellnum)
    trackmatrix(1:t-1,cellnum) = NaN;
    ucell = 1:size(trackmatrix,2);
    trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
    Tracked.trackmatrix = trackmatrix;
    if ~isempty(idx)
        plotTheChosenCells(idx);
    end
end

setSceneAndTime;
end
function destroybuttonSubsequent_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked  
trackmatrix = Tracked.trackmatrix;
t = ImageDetails.Frame;
ginputnum=[];
[cxx,cyy,~,button] = identifyCellbyClick(t,ginputnum);
if isempty(button)
    return
end
%find cell
[idx,~,~] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,t);
alltrackframe = trackmatrix(t,:);
cellnum = find(ismember(alltrackframe,idx));
if ~isempty(cellnum)
    trackmatrix(t+1:end,cellnum) = NaN;
    ucell = 1:size(trackmatrix,2);
    trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
    Tracked.trackmatrix = trackmatrix;
    if ~isempty(idx)
        plotTheChosenCells(idx);
    end
end
setSceneAndTime;
plotTheChosenCells(idx);
end
function destroybutton_Callback(~,~)
%remove entire track from all frames
global ImageDetails Tracked
trackmatrix = Tracked.trackmatrix;
t = ImageDetails.Frame;
ginputnum=[];
[cxx,cyy,~,button] = identifyCellbyClick(t,ginputnum);
if isempty(button)
    return
end
%find cell
[idx,~,~] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,t);
alltrackframe = trackmatrix(t,:);
cellnum = find(ismember(alltrackframe,idx));
if ~isempty(cellnum)
    trackmatrix(:,cellnum) = NaN;
    ucell = 1:size(trackmatrix,2);
    trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
    Tracked.trackmatrix = trackmatrix;
end
setSceneAndTime;
end
function breaklink_callback(~,~)
%remove entire track from all frames
global ImageDetails Tracked
trackmatrix = Tracked.trackmatrix;
t = ImageDetails.Frame;
ginputnum=[];
[cxx,cyy,~,button] = identifyCellbyClick(t,ginputnum);
if isempty(button)
    return
end
%find cell
[idx,~,~] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,t);
alltrackframe = trackmatrix(t,:);
cellnum = find(ismember(alltrackframe,idx));
if ~isempty(cellnum)
    trackmatrix(t,cellnum) = NaN;
    ucell = 1:size(trackmatrix,2);
    trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
    Tracked.trackmatrix = trackmatrix;
end
setSceneAndTime;
end
function breakEdge_callback(~,~)
%remove entire track from all frames
global ImageDetails Tracked timeFrames
trackmatrix = Tracked.trackmatrix;
ArrayStruct = Tracked.arrayStruct;
segment_Pixels_array = ArrayStruct.pixels;

t = ImageDetails.Frame;
imgsize =ImageDetails.ImgSize;
pxborder = true(imgsize(1),imgsize(2));
pxborder(2:imgsize(1)-1,2:imgsize(2)-1) = false;
edgepixels = find(pxborder);
for t = 1:timeFrames
    %find the cells
    PX = segment_Pixels_array{t};
    idxlog = cellfun(@(x) sum(ismember(edgepixels,x))>0,PX,'UniformOutput',1);
    idx = find(idxlog);
    alltrackframe = trackmatrix(t,:);
    cellnum = find(ismember(alltrackframe,idx));
    if ~isempty(cellnum)
        trackmatrix(t,cellnum) = NaN;
    end
end
ucell = 1:size(trackmatrix,2);
trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
Tracked.trackmatrix = trackmatrix;
setSceneAndTime
end

%choose the cells you want
function chosenOnes_Callback(~,~)
%keep the entire track of selected cells
global ImageDetails Tracked
trackmatrix = Tracked.trackmatrix;
t = ImageDetails.Frame;
ginputnum=[];
[cxx,cyy,~,button] = identifyCellbyClick(t,ginputnum);
if isempty(button)
    return
end
%find cell
[idx,~,~] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,t);
alltrackframe = trackmatrix(t,:);
cellnum = find(ismember(alltrackframe,idx));
if ~isempty(cellnum)
    cellmatlog = true(1,size(trackmatrix,2));
    cellmatlog(cellnum) = false;
    trackmatrix(:,cellmatlog) = NaN;
    ucell = 1:size(trackmatrix,2);
    trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
    Tracked.trackmatrix = trackmatrix;
end
setSceneAndTime;
end
function chosenOnesAllOnFrame_Callback(~,~)
%keep the entire tracks from only the cells on the frame
global ImageDetails Tracked
trackmatrix = Tracked.trackmatrix;
t = ImageDetails.Frame;
alltrackframe = trackmatrix(t,:);
idx = find(~isnan(alltrackframe));
if ~isempty(idx)
    cellnum = idx;
    cellmatlog = true(1,size(trackmatrix,2));
    cellmatlog(cellnum) = false;
    trackmatrix(:,cellmatlog) = NaN;
    ucell = 1:size(trackmatrix,2);
    trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
    Tracked.trackmatrix = trackmatrix;
end
setSceneAndTime;
end
function erodeOnes_Callback(~,~)
global  ImageDetails Tracked


ArrayStruct = Tracked.arrayStruct;
trackmatrix = Tracked.trackmatrix;
segment_Pixels_array = ArrayStruct.pixels;
PXarray = [];
tvec = [];
dropArray = [];
newArray = [];
keepArray=[];
allArray=[];

t = ImageDetails.Frame;
ginputnum=[];
[cxx,cyy,~,button] = identifyCellbyClick(t,ginputnum);
if isempty(button)
    return
end
%find cell
[idxall,~,notACell] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,t);
if notACell
    return
end


idxu = unique(idxall);
erodemag = zeros(1,length(idxu));
for j = 1:length(idxu)
    erodemag(j) = sum(ismember(idxall,idxu(j)));
end

for idxi = 1:length(idxu)
    idx = idxu(idxi);
    emag = erodemag(idxi);
    %loop everythihng in here
    
    alltrackframe = trackmatrix(t,:);
    cellnum = find(ismember(alltrackframe,idx));
    celltrack = trackmatrix(:,cellnum);
    timesyo = find(~isnan(celltrack));
    
    for ti = timesyo'
        index = trackmatrix(ti,cellnum);
        
        if ~isempty(index)
            PX = segment_Pixels_array{ti};
            pxGroup = vertcat(PX{index}); %join pixel values from all overlapping cells
            img = false(ImageDetails.ImgSize);
            img(pxGroup) = true;
            imge =imerode(img,strel('disk',emag));
            pxnew = {find(imge==true)};
            
            %define keepidx,dropidx,newidx
            keepidx = index;
            dropidx = [];
            newidx = [];
            allidx = keepidx;
        else
            %         pxnew = {px};
            %
            %         %define keepidx,dropidx,newidx
            %         keepidx = [];
            %         dropidx = [];
            %         newidx = length(PX)+1;
            %         allidx = newidx;
        end
        
        %         trackmatrix = updateTrackMatrix(trackmatrix,ti,dropidx,keepidx,newidx);
        dropArray = horzcat(dropArray,{dropidx});
        keepArray = horzcat(keepArray,{keepidx});
        newArray = horzcat(newArray,{newidx});
        allArray = horzcat(allArray,{allidx});
        PXarray = horzcat(PXarray,{pxnew});
        tvec = horzcat(tvec,ti);
        
        if isempty(tvec)
            return
        else
        end
        disp(ti)
    end
end    
    %array struct does not work here because it runs twice!
    tic
    ArrayStruct = updateArrayStruct(ArrayStruct,PXarray,tvec,dropArray,newArray,keepArray,allArray,trackmatrix,ImageDetails.ImgSize); %this is not ideal for updating trackmatrix or removing
    Tracked.arrayStruct = ArrayStruct;
    Tracked.trackmatrix = trackmatrix;
    toc
    ucell = 1:size(trackmatrix,2);
    trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
    
    ImageDetails.Frame = t;
setSceneAndTime;

end
function dilateOnes_Callback(~,~)
global  ImageDetails Tracked


ArrayStruct = Tracked.arrayStruct;
trackmatrix = Tracked.trackmatrix;
segment_Pixels_array = ArrayStruct.pixels;
PXarray = [];
tvec = [];
dropArray = [];
newArray = [];
keepArray=[];
allArray=[];

t = ImageDetails.Frame;
ginputnum=[];
[cxx,cyy,~,button] = identifyCellbyClick(t,ginputnum);
if isempty(button)
    return
end
%find cell
[idxall,~,notACell] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,t);
if notACell
    return
end


idxu = unique(idxall);
erodemag = zeros(1,length(idxu));
for j = 1:length(idxu)
    erodemag(j) = sum(ismember(idxall,idxu(j)));
end

for idxi = 1:length(idxu)
    idx = idxu(idxi);
    emag = erodemag(idxi);
    %loop everythihng in here
    
    alltrackframe = trackmatrix(t,:);
    cellnum = find(ismember(alltrackframe,idx));
    celltrack = trackmatrix(:,cellnum);
    timesyo = find(~isnan(celltrack));
    
    for ti = timesyo'
        index = trackmatrix(ti,cellnum);
        
        if ~isempty(index)
            PX = segment_Pixels_array{ti};
            pxGroup = vertcat(PX{index}); %join pixel values from all overlapping cells
            img = false(ImageDetails.ImgSize);
            img(pxGroup) = true;
            imge = imdilate(img,strel('disk',emag));
            pxnew = {find(imge==true)};
            
            %define keepidx,dropidx,newidx
            keepidx = index;
            dropidx = [];
            newidx = [];
            allidx = keepidx;
        else
            %         pxnew = {px};
            %
            %         %define keepidx,dropidx,newidx
            %         keepidx = [];
            %         dropidx = [];
            %         newidx = length(PX)+1;
            %         allidx = newidx;
        end
        
        %         trackmatrix = updateTrackMatrix(trackmatrix,ti,dropidx,keepidx,newidx);
        dropArray = horzcat(dropArray,{dropidx});
        keepArray = horzcat(keepArray,{keepidx});
        newArray = horzcat(newArray,{newidx});
        allArray = horzcat(allArray,{allidx});
        PXarray = horzcat(PXarray,{pxnew});
        tvec = horzcat(tvec,ti);
        
        if isempty(tvec)
            return
        else
        end
        disp(ti)
    end
end    
    %array struct does not work here because it runs twice!
    tic
    ArrayStruct = updateArrayStruct(ArrayStruct,PXarray,tvec,dropArray,newArray,keepArray,allArray,trackmatrix,ImageDetails.ImgSize); %this is not ideal for updating trackmatrix or removing
    Tracked.arrayStruct = ArrayStruct;
    Tracked.trackmatrix = trackmatrix;
    toc
    ucell = 1:size(trackmatrix,2);
    trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
    
    ImageDetails.Frame = t;
setSceneAndTime;

end


%plot your cells!
function psettings = PlotSettings_callback(~,~)
global dirStruct expDetailsStruct frameToLoad
cd(dirStruct.exportdir)

queryName = strcat(expDetailsStruct.expDateStr,'*DoseAndScene*.mat');
filelist = dir(queryName);
if isempty(queryName)
    prompt = {'tgfbeta frame','last frame'};
    dlg_title = 'frames where cells must be tracked...';
    inputdlgOutput = inputdlg(prompt,dlg_title);
    framesThatMustBeTracked = cellfun(@num2str,inputdlgOutput,'UniformOutput',1);
else
    A = load(char(filelist.name));
    dosestruct = A.dosestruct;
    if length(frameToLoad)<10
        framesThatMustBeTracked = [dosestruct(1).tgfFrame 1+dosestruct(1).tgfFrame];
    else
        framesThatMustBeTracked = [dosestruct(1).tgfFrame 10+dosestruct(1).tgfFrame];
    end
end
psettings.framesThatMustBeTracked = framesThatMustBeTracked;
end
function [fc,abs,pidx] = Plot_callback(~,~)
global medmeantot segStruct xAxisLimits cmaplz SecondPlotAxes Tracked ImageDetails pathStruct timeFrames PlotAxes togStruct psettings cmap

trackmatrix = Tracked.trackmatrix;
ArrayStruct = Tracked.arrayStruct;
segment_Pixels_array = ArrayStruct.pixels;
if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
end
framesThatMustBeTracked = psettings.framesThatMustBeTracked;


for jy = 1
    PX = segment_Pixels_array{jy};
    makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
end

makeIMGidx = find(makeIMG==1);

smooththat=0;
[plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,pathStruct,ImageDetails,timeFrames,smooththat);


    
    %second axes, plot fold-change
    mmsstr = medmeantot;
    fcplotstr = 'FC';
    ylabelstr = 'fold change';
    if strcmp(ImageDetails.Channel,segStruct.background_seg)
        plotstr = ['Smad' fcplotstr];
    elseif strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
        plotstr = ['mkate' fcplotstr];
    elseif strcmp()
        
    else
        plotstr = ['Smad' fcplotstr];
    end
    plotMat = plotStructUI.(mmsstr).(plotstr);
    ylimit = [prctile(plotMat(:),0) prctile(plotMat(:),100)];
    
    idx = true(size(plotMat,1),1);
    fc = plotMat;
    pidx = idx;
    cmapl = cmaplz;
    idxa = find(idx==1);
    h = plot(SecondPlotAxes,plotMat(idx,:)','LineWidth',2);
    %label with nvalue
    nvalue = max(sum(~isnan(plotMat),1));
    nstr = num2str(nvalue);
    tstr = ['n = ' nstr ' tracks'];
    t = text(SecondPlotAxes,0,0,tstr);
    t.Units = 'normalized';
    t.Position = [0.01 0.95];
    t.FontSize = 8;
    %update color traces to match image
    if togStruct.displayTrackingToggle ==1
        for i=1:length(h)
            h(i).Color = cmapl(idxa(i),:);
        end
%         colormap(cmap);
    end
    
    SecondPlotAxes.XLim = ([xAxisLimits(1) xAxisLimits(2)]);
    SecondPlotAxes.YLim = (ylimit);
    SecondPlotAxes.YLabel.String = [ylabelstr ' ' mmsstr];
    SecondPlotAxes.XLabel.String = 'frames';
    SecondPlotAxes.XGrid = 'on';
    SecondPlotAxes.YGrid = 'on';
    SecondPlotAxes.Color = [0.95 0.95 0.95];
    SecondPlotAxes.FontSize = 7;
    
    
    %main axes, plot abundance
    fcplotstr = '';
    ylabelstr = 'abundance';
    if strcmp(ImageDetails.Channel,segStruct.background_seg)
        plotstr = ['Smad' fcplotstr];
    elseif strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
        plotstr = ['mkate' fcplotstr];
    else
        plotstr = ['Smad' fcplotstr];
    end
    plotMat = plotStructUI.(mmsstr).(plotstr);
    ylimit = [prctile(plotMat(:),0) prctile(plotMat(:),100)];
    abs = plotMat;
    h = plot(PlotAxes,plotMat(idx,:)','LineWidth',2);
    nvalue = sum(idx);
    nstr = num2str(nvalue);
    tstr = ['n = ' nstr ' tracks'];
    t = text(PlotAxes,0,0,tstr);
    t.Units = 'normalized';
    t.Position = [0.01 0.95];
    t.FontSize = 8;
    if togStruct.displayTrackingToggle ==1
        for i=1:length(h)
            h(i).Color = cmapl(i,:);
        end
%         colormap(cmap);
    end
    
    PlotAxes.XLim = ([xAxisLimits(1) xAxisLimits(2)]);
    PlotAxes.YLim = ([prctile(plotMat(:),0.5)./1.2 prctile(plotMat(:),99.5).*1.2]);
    PlotAxes.YLabel.String = [ylabelstr ' ' mmsstr];
    PlotAxes.XLabel.String = 'frames';
    PlotAxes.XGrid = 'on';
    PlotAxes.YGrid = 'on';
    PlotAxes.Color = [0.95 0.95 0.95];
    PlotAxes.Title.String = [ImageDetails.Channel ' in tracked cells'];
    PlotAxes.FontSize = 7;

end
function plotAxis_callback(~,~)
global xAxisLimits
prompt = {'xmin','xmax'};
dlg_title = 'set x axis limits...';
inputdlgOutput = inputdlg(prompt,dlg_title);
xAxisLimits = cellfun(@str2num,inputdlgOutput,'UniformOutput',1);
Plot_callback([],[])
end
function PlotSpecificCell_callback(~,~)
global Tracked ImageDetails
    ginputnum=[];
    t = ImageDetails.Frame;
    [cxx,cyy,tkeeper,button] = identifyCellbyClick(t,ginputnum);
    if isempty(button)
        return
    end
    %find cell
    [idx,~,notACell] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,tkeeper);
    if ~notACell && ~isempty(idx)
        plotTheChosenCells(idx);
    end
end
function PlotSpecificCellIteratively_callback(~,~)
global Tracked ImageDetails
button = 1;
while button == 1
    ginputnum=1;
    t = ImageDetails.Frame;
    [cxx,cyy,tkeeper,button] = identifyCellbyClick(t,ginputnum);
    if isempty(button)
        return
    end
    %find cell
    [idx,~,notACell] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,tkeeper);
    if ~notACell && ~isempty(idx)
        plotTheChosenCells(idx);
    end
end
end
function plotTheChosenCells(idxs) % Display mesh plot of the currently selected data.
global togStruct cmaplz xAxisLimits ThirdPlotAxes Tracked ImageDetails pathStruct timeFrames psettings segStruct medmeantot

if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
end


trackmatrix = Tracked.trackmatrix;
if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
end
framesThatMustBeTracked = psettings.framesThatMustBeTracked;


t = ImageDetails.Frame;
smooththat=0;
[plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,pathStruct,ImageDetails,timeFrames,smooththat);

alltrackframe = trackmatrix(t,:);
makeIMGidx = find(ismember(alltrackframe,idxs));
basalnanidx = ~isnan(trackmatrix(framesThatMustBeTracked(1),makeIMGidx));

lrstrArray = {'left', 'right'};
mmsstrArray = {medmeantot,'total'};
for lrnum = [1 2]
    lrstr = lrstrArray{lrnum}; %determine left or right axis
    mmsstr = mmsstrArray{lrnum};
    
    if isempty(find(basalnanidx==0))
        fcplotstr = 'FC';
        ylabelstr = 'fold change';
    else
        fcplotstr = '';
        ylabelstr = 'abundance';
    end
    
    if strcmp(ImageDetails.Channel,segStruct.background_seg)
        plotstr = ['Smad' fcplotstr];
    elseif strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
        plotstr = ['mkate' fcplotstr];
    else
        plotstr = ['Smad' fcplotstr];
    end
    plotMat = plotStructUI.(mmsstr).(plotstr);
    if lrnum == 1
        ylimit = [prctile(plotMat(:),0) max([4 prctile(plotMat(:),100)])];
    end
    
    if togStruct.displayTrackingToggle == 1
        %set colormap
        cmapl = cmaplz;
        if lrnum >1
            darkenfactor = 1.5;
        else
            darkenfactor = 1;
        end
        cmapnew = cmapl(makeIMGidx,:)/darkenfactor;
    else
        cmapnew = feval('lines',length(makeIMGidx));
    end
        cbut = num2cell(cmapnew,2);
    
    
    yyaxis(ThirdPlotAxes,lrstr)
    toplot = plotMat(makeIMGidx,:);
    if lrnum == 1
        h1 = plot(ThirdPlotAxes,toplot','LineWidth',2);
        [h1.LineStyle] = deal('-');
        [h1.DisplayName] = deal(mmsstr);
        [h1.Marker] = deal('none');
        if togStruct.displayTrackingToggle ==1
            [h1.Color] = cbut{:};
        end
    else
        h2 = plot(ThirdPlotAxes,toplot','LineWidth',2);
        [h2.LineStyle] = deal(':');
        [h2.DisplayName] = deal(mmsstr);
        [h2.Marker] = deal('none');
        if togStruct.displayTrackingToggle ==1
            [h2.Color] = cbut{:};
        end
    end

    ThirdPlotAxes.XLim = ([xAxisLimits(1) xAxisLimits(2)]);
    ThirdPlotAxes.YLim = (ylimit);
    ThirdPlotAxes.YLabel.String = [ylabelstr ''  mmsstr];
    ThirdPlotAxes.XLabel.String = 'frames';
    ThirdPlotAxes.XGrid = 'on';
    ThirdPlotAxes.YGrid = 'on';
    ThirdPlotAxes.Color = [0.95 0.95 0.95];
    ThirdPlotAxes.FontSize = 7;
    
    

    if lrnum ==1
        children2 = findobj(ThirdPlotAxes,'Type','Legend');
        delete(children2)
    end
end
if ~isempty(h1)
    leg = legend(ThirdPlotAxes,[h1(1) h2(1)],mmsstrArray{1},mmsstrArray{2});
    % leg.String = mmsstrArray;
    pos = ThirdPlotAxes.Position;
    leg.Location = 'southoutside';
    ThirdPlotAxes.Position =pos;
    leg.FontSize = 8;
    leg.Orientation = 'horizontal';
end
% l.Box = 'off';
end
function [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,pathStruct,ImageDetails,timeFrames,smooththat)


trackmatrix = Tracked.trackmatrix;
ArrayStruct = Tracked.arrayStruct;

%extract pixel intensities
segment_Pixels_array = ArrayStruct.pixels;
segment_cellFluor_array = ArrayStruct.cellFluor;
segment_nucFluor_array = ArrayStruct.nucFluor;

segment_cellcyto_array = ArrayStruct.cellcyto;
segment_nuccyto_array = ArrayStruct.nuccyto;

cellFluorInit = ArrayStruct.cellFluor{1};
mmsstrArray = {'median','mean','total','nuccyto'};
% for k = 1:size(cellFluorInit,2)
for k = 1:length(mmsstrArray)
    mmsstr = mmsstrArray{k};

    
    %now update all the fields
    plotTracesCell = cell(size(trackmatrix));
    plotTracesCell(:) = {NaN};
    nucFluorMatrix = nan(size(trackmatrix,2),timeFrames);
    cellFluorMatrix = nan(size(trackmatrix,2),timeFrames);
    
    for i = 1:timeFrames
        a=i;
        trackvals=trackmatrix(i,:);
        trackidx = ~isnan(trackvals);
        tracknums = trackvals(trackidx);
        
        px = segment_Pixels_array{a};
        nucFluor = segment_nucFluor_array{a};
        cellFluor = segment_cellFluor_array{a};
        if ~isempty(px)
            %first pixels
            pxpx = plotTracesCell(i,:);
            pxpx(trackidx) = px(tracknums);
            plotTracesCell(i,:) = pxpx;
            
            if strcmp(mmsstr,'nuccyto')
                k=1;
                nuccyto = segment_nuccyto_array{a}';
                cellcyto = segment_cellcyto_array{a}';
                nucFluorMatrix(trackidx,i) = nucFluor(tracknums,k)./nuccyto(tracknums); %choose median, mean or total (k=1,2,3 respectively)
                cellFluorMatrix(trackidx,i) = cellFluor(tracknums,k)./cellcyto(tracknums);
            else
                nucFluorMatrix(trackidx,i) = nucFluor(tracknums,k); %choose median, mean or total (k=1,2,3 respectively)
                cellFluorMatrix(trackidx,i) = cellFluor(tracknums,k);
            end
        end
    end
    
    Smad = cellFluorMatrix;
    mkate = nucFluorMatrix;
    
    basalSUB = max([1 framesThatMustBeTracked(1)-3]);
    
    basalsmad = nanmean(Smad(:,basalSUB:framesThatMustBeTracked(1)),2);
    SmadFC = zeros(size(Smad),'single');
    for i = 1:size(Smad,2)
        SmadFC(:,i) = Smad(:,i)./basalsmad;
    end
    
    basalmkate = nanmean(mkate(:,basalSUB:framesThatMustBeTracked(1)),2);
    mkateFC = zeros(size(mkate),'single');
    for i = 1:size(mkate,2)
        mkateFC(:,i) = mkate(:,i)./basalmkate;
    end
    
    nanvec= Smad(:,framesThatMustBeTracked(1));
    nanidx = isnan(nanvec);
    
    SmadFC(nanidx,:) = NaN;
    mkateFC(nanidx,:) = NaN;
    
    %     Smad,Cfp,mkate,CfpFC,SmadFC,mkateFC,Smadbkg,Cfpbkg,mkatebkg
    plotStructUI.(mmsstr).Smad = Smad;
    plotStructUI.(mmsstr).mkate = mkate;
    plotStructUI.(mmsstr).SmadFC = SmadFC;
    plotStructUI.(mmsstr).mkateFC = mkateFC;
end

end




function ncpx = nuccytopxspitter(px,imdim)

newim = zeros(imdim(1),imdim(2));
for j = 1:length(px)
    subpx = px{j};
    newim(subpx) =j;
end

se = strel('disk',3);
newimd2 = imdilate(newim,se);
se = strel('disk',1);
newimd1 = imdilate(newim,se);
newimd2log = newimd2>0;
newimd1log = newimd1>0;
newimlog = ~newimd1log & newimd2log;
newimnew = newimd2;
newimnew(~newimlog) = 0;

for j = 1:length(px)
    newpx = find(newimnew==j);
    px{j} = newpx;
end
ncpx = px;

end


function plotStruct = plotthemfunctionToStructure(Tracked,idScene,pathStruct,timeFrames,segStruct,ctinput)
cd(pathStruct.mstackPath)
%no bleach correction option yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   open the image files   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  open smad img  %
cd(pathStruct.mstackPath)
ff = dir(strcat('*',idScene,'*',segStruct.background_seg,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
filename = char(ff.name);
channelfileObject = matfile(filename);
cellQ_imgstack = channelfileObject.flatstack;
imdim = size(cellQ_imgstack);
clear cellQ_imgstack

plotStruct = struct();
trackmatrix = Tracked.trackmatrix;
ArrayStruct = Tracked.arrayStruct;

%extract pixel intensities
segment_Pixels_array = ArrayStruct.pixels;
segment_cellFluor_array = ArrayStruct.cellFluor;
segment_nucFluor_array = ArrayStruct.nucFluor;
segment_Centroid_array = ArrayStruct.centroid;
centroids = segment_Centroid_array{1};
%now update all the fields
plotTracesCell = cell(size(trackmatrix));
plotTracesCell(:) = {NaN};
nucFluorMatrix = nan(size(trackmatrix,2),timeFrames);
cellFluorMatrix = nan(size(trackmatrix,2),timeFrames);
if isempty(centroids)
    centnew = nan(size(trackmatrix,2),2,timeFrames);
else
    centnew = nan(size(trackmatrix,2),size(centroids,2),timeFrames);
end
k=1; %choose median
for i = 1:timeFrames
    a=i;
    trackvals=trackmatrix(i,:);
    trackidx = ~isnan(trackvals);
    tracknums = trackvals(trackidx);
    
    px = segment_Pixels_array{a};
    nucFluor = segment_nucFluor_array{a};
    cellFluor = segment_cellFluor_array{a};
    centroids = segment_Centroid_array{a};
    if ~isempty(px)
        %first pixels
        pxpx = plotTracesCell(i,:);
        pxpx(trackidx) = px(tracknums);
        plotTracesCell(i,:) = pxpx;
        nucFluorMatrix(trackidx,i) = nucFluor(tracknums,k);
        cellFluorMatrix(trackidx,i) = cellFluor(tracknums,k);
        centnew(trackidx,:,i) = centroids(tracknums,:);
    end
end


ncRatio = cell(size(trackmatrix));
ncRatio(:) = {NaN};
% nucFluorMatrix = nan(size(trackmatrix,2),timeFrames);
% cellFluorMatrix = nan(size(trackmatrix,2),timeFrames);
% centnew = nan(size(trackmatrix,2),size(centroids,2),timeFrames);
k=1; %choose median

for i = 1:timeFrames
    a=i;
    trackvals=trackmatrix(i,:);
    trackidx = ~isnan(trackvals);
    tracknums = trackvals(trackidx);
    px = segment_Pixels_array{a};
    
    if ~isempty(px) 
        ncpx = nuccytopxspitter(px,imdim);
        px = ncpx;

            
        %first pixels
        pxpx = ncRatio(i,:);
        pxpx(trackidx) = px(tracknums);
        ncRatio(i,:) = pxpx;
    end
    
end

ncRatio = ncRatio';
plotTracesCell = plotTracesCell';



cd(pathStruct.mstackPath)
%no bleach correction option yet
%% 
nctinput = ctinput;
for i = 1:length(ctinput)
    ctstr = ctinput{i};
    a = regexp(ctstr,'\s');
    astr = regexp(ctstr,'DIC');
    if ~isempty(astr)
        ctstr = 'DIC';
    else
        ctstr(a) = [];
    end
    nctinput{i} = ctstr;
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   open the image files   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chanstruct = struct();
for i = 1:length(nctinput)
    nctstr = nctinput{i};
%  open smad img  %
cd(pathStruct.mstackPath)
ff = dir(['*' idScene '*' nctstr '*']);
%         ff = dir(strcat(ImageDetails.Channel,'*'));
filename = char(ff.name);
channelfileObject = matfile(filename);
% cellQ_imgstack = channelfileObject.flatstack;
imgstack = channelfileObject.flatstack;
chanstruct.(nctstr) = imgstack;
end


%% open segmentation images
% open background Logical img  %
cd(pathStruct.segmentPath)
ff = dir(strcat('*',idScene,'*',segStruct.background_seg,'*'));
if length(ff)>1
    ff = dir(strcat('*',idScene,'*',segStruct.background_seg,'*background*'));
end
%         ff = dir(strcat(ImageDetails.Channel,'*'));
filename = char(ff.name);
channelfileObject = matfile(filename);
bkglogimgstack = channelfileObject.IfFinal;

% open nuclear Logical img  %
cd(pathStruct.segmentPath)
ff = dir(strcat('*',idScene,'*',segStruct.nucleus_seg,'*'));
if length(ff)>1
    ff = dir(strcat('*',idScene,'*',segStruct.nucleus_seg,'*nucleus*'));
end
%         ff = dir(strcat(ImageDetails.Channel,'*'));
filename = char(ff.name);
channelfileObject = matfile(filename);
nuclogimgstack = channelfileObject.IfFinal;
bkglogimgstack(nuclogimgstack) = true;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subimgstruct = chanstruct;
bkgstruct = struct();
nctnames = fieldnames(chanstruct);
for i = 1:length(nctnames)
    nctstr = nctnames{i};
    imgstack = chanstruct.(nctstr);
    subimgstack = imgstack;
    
    %perform bkg subtraction
    imgBKG = zeros(1,timeFrames,'single');
    for k=1:timeFrames
        bkglog = ~bkglogimgstack(:,:,k);
        img = single(imgstack(:,:,k));
        
        %background subtraction is just subtraction with a value
        imgBKG(k) = nanmedian(img(bkglog));
        subimgstack(:,:,k) = img-imgBKG(k);
    end
    bkgstruct.(nctstr) = imgBKG;
    subimgstruct.(nctstr) = subimgstack;
end




%extract pixel intensities
for n = 1:length(nctnames)
    nstr = nctnames{n};
    
    img_pxls = cell(size(plotTracesCell,1),size(plotTracesCell,2));
    img_nc_pxls = img_pxls;
    imgstack = subimgstruct.(nstr);
    cellBKG =  bkgstruct.(nctstr);
    
    for i = 1:size(plotTracesCell,2)
        img_img = single(squeeze(imgstack(:,:,i)));
        for j=1:size(plotTracesCell,1)
            pxidx = plotTracesCell{j,i};
            ncidx = ncRatio{j,i};
            
            if ~isnan(pxidx)
                img_pxls(j,i) = {img_img(pxidx)};
                img_nc_pxls(j,i) = {img_img(ncidx)};
            else
                img_pxls(j,i) = {single(13579)};
                img_nc_pxls(j,i) = {single(13579)};
            end
            
        end
    end
    
    %determine median pxl intensities
    img = cellfun(@nanmedian,img_pxls,'UniformOutput',1);
    img(img==single(13579)) = NaN;
    for i = 1:size(img_pxls,1)
        plotStruct(i).(['medianNuc' nstr]) = img(i,:);
    end
    
    %determine median pxl intensities
    img_nc = cellfun(@nanmedian,img_nc_pxls,'UniformOutput',1);
    img_nc(img_nc==single(13579)) = NaN;
    for i = 1:size(img_nc_pxls,1)
        plotStruct(i).(['nuccytoNuc' nstr]) = img(i,:)./img_nc(i,:);
    end
    
    
    %determine total pxl intensities
    img = cellfun(@nansum,img_pxls,'UniformOutput',1);
    img(img==single(13579)) = NaN;
    for i = 1:size(img_pxls,1)
        plotStruct(i).(['totalNuc' nstr]) = img(i,:);
    end
    
    %determine mean pxl intensities
    img = cellfun(@nanmedian,img_pxls,'UniformOutput',1);
    img(img==single(13579)) = NaN;
    for i = 1:size(img_pxls,1)
        plotStruct(i).(['medianNuc' nstr]) = img(i,:);
        plotStruct(i).(['medianBKG' nstr])= cellBKG;
    end
    
    
    %determine mean pxl intensities
    img = cellfun(@nanmean,img_pxls,'UniformOutput',1);
    img(img==single(13579)) = NaN;
    for i = 1:size(img_pxls,1)
        plotStruct(i).(['meanNuc' nstr]) = img(i,:);
        plotStruct(i).(['meanBKG' nstr])= cellBKG;
    end
    
    
    %determine total pxl intensities
    img = cellfun(@length,img_pxls,'UniformOutput',1);
    img(img==single(13579)) = NaN;
    for i = 1:size(img_pxls,1)
        plotStruct(i).(['NucArea' nstr]) = img(i,:);
    end
    
    
    
    for i = 1:size(img_pxls,1)
        plotStruct(i).Centroid = centnew(i,:,:);
    end
    
    for i = 1:size(img_pxls,1)
        plotStruct(i).CellNum = trackmatrix(:,i);
    end
    
    stophere=1;
end
end


%%%%comments!!!
function xy = getxy(~,~)
global Tracked togStruct psettings ImageDetails
if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     t=length(Tracked);
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


framesThatMustBeTracked = psettings.framesThatMustBeTracked;
PX = Tracked{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
makeIMG = false(length(framesThatMustBeTracked),length(PX));
for jy = 1:length(framesThatMustBeTracked)
    PX = Tracked{framesThatMustBeTracked(jy)}.Cellz.PixelIdxList;
    %                 makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)==1,PX,'UniformOutput',1)); %choose only the cells without NAN
    makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
end
makeIMG = makeIMG(1,:)&makeIMG(2,:);
makeIMGidx = find(makeIMG==1);



PXX = Tracked{t}.Cellz.PixelIdxList;
PX = PXX(makeIMG);
mx = nan(1,length(PX));
my = nan(1,length(PX));
for j = 1:length(PX)
    px = PX{j};
    y = rem(px-1,ImageDetails.ImgSize(1))+1; %these two lines replace ind2sub
    x = (px-y)/ImageDetails.ImgSize(2) + 1;  %these two lines replace ind2sub
    %%%% sort and select middle index is faster than median %%%
    sx = sort(x);
    sy = sort(y);
    pseudomean = round(length(sx)./2);
    if pseudomean == 0
        mx(j) = NaN;
        my(j) = NaN;
    else
        mx(j) = sx(pseudomean);  %use the make Centroids index to keep the centroids the same color when plotting
        my(j) = sy(pseudomean);
    end
end
xy = horzcat(mx',my');


end
function dispxy(xy)
global frameToLoad displaycomments Tracked



for i = 1:size(xy,1)
    text(xy(i,1)+1,xy(i,2)+1,num2str(i),'Color',[0 0 0],'FontSize',12);hold on
    text(xy(i,1),xy(i,2),num2str(i),'Color',[1 1 0],'FontSize',12);hold on
end

t = length(Tracked);
xxyy=xy;
xxyy(:,1) = xxyy(:,1)+10;
xxyy(:,2) = xxyy(:,2)+10;
if t==length(frameToLoad)
    if displaycomments==1
        comments = Tracked{t}.comments;
        commentpos = Tracked{t}.commentpos;
        %         for i = 1:size(xxyy,1)
        %         for i = 1:length(comments)
        lll = length(comments);
        lllx = size(xxyy,1);
        if lll <lllx
            lcc = lll;
        else
            lcc =lllx;
        end
        for i = 1:lcc
            text(xxyy(i,1)+1,xxyy(i,2)+1,comments{i},'Color',[0 0 0],'FontSize',8);hold on
            text(xxyy(i,1),xxyy(i,2),comments{i},'Color',[1 1 0],'FontSize',8);hold on
        end
    end
end
end



function [cellnums,comments,commentpos,cellidx]=updatecomments(xy)
%choose the cells you want to comment on
global ImageDetails Tracked
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;



ix = cellfun(@length,PX,'UniformOutput',1);
idxs = ix>1;

idx = find(idxs==1);


fnames = fieldnames(Tracked{t});
if sum(strcmp(fnames,'comments'))
    oldcomments = Tracked{t}.comments;
    oldcommentpos = Tracked{t}.commentpos;
    comments=[];
    commentpos =[];
else
    oldcomments = [];
    oldcommentpos = [];
    comments = [];
    commentpos =[];
end

comments = cell(1,size(xy,1));
commentpos = zeros(1,size(xy,1));
cellidx =zeros(1,size(xy,1));

indies = sub2ind([ImageDetails.ImgSize],xy(:,2),xy(:,1));

if ~isempty(oldcommentpos)
    for jim = idx
        cycle=0;
        for i = 1:length(oldcommentpos)
            if i==5
                stophere=1;
            end
            cycle=cycle+1;
            px = PX{jim};
            alreadycommented = ismember(oldcommentpos(i),px);
            indiidx = find(ismember(indies,px)==1);
            if alreadycommented ==1
                comments{indiidx} = oldcomments{i};
                commentpos(indiidx) = oldcommentpos(i);
                cellidx(indiidx) = jim;
            else
                
            end
        end
    end
    
    
    
else
    
    
    
    for i=idx
        px = PX{i};
        indiidx = find(ismember(indies,px)==1);
        commentpos(indiidx) = indies(indiidx);
    end
    
end

%try this out 2016_02_17
aa = find(commentpos==0);
commentpos(aa) = indies(aa);


for i=1:length(comments)
    cellnums{i} = strcat('cell#',num2str(i));
    if isempty(comments{i})
        comments{i}='';
    end
end

end
function comments = nowaddcomments(xy,cellnums,comments,commentpos)
global Tracked
t=length(Tracked);



prompt = cellnums;
dlg_title = 'Comment for selected cells:';
num_lines = 1;
defaultans = comments;
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

if ~isempty(answer)
    comments = answer;
    
    Tracked{t}.comments= comments;
    Tracked{t}.commentpos = commentpos;
end
% displaycomments =1;

end
function comment_Callback(~,~)
global displaycomments Tracked


fnames = fieldnames(Tracked{length(Tracked)});
if sum(strcmp(fnames,'comments'))
    displaycomments =1;
else
    displaycomments =0;
end
xy = getxy([],[]);
dispxy(xy)
[cellnums,comments,commentpos,cellidx]=updatecomments(xy);
comments = nowaddcomments(xy,cellnums,comments,commentpos);
displaycomments =1;
dispxy(xy)
%    displaycomments =0;


end



function exportTrackedCells(~,~)
global segStruct pathStruct ExportNameKey ExportName dirStruct expDetailsStruct SceneList timeFrames togStruct psettings channelstoinput
exportStruct = struct();



if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
end

%remove all global variables before parfor loop
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
eNameKey = ExportNameKey;
eName = ExportName;
paStruct = pathStruct;
sStruct = segStruct;
tPath = pathStruct.trackingPath;
sList = SceneList;
tFrames = timeFrames;
ctinput = channelstoinput;
cd(tPath)
cd ..

plotStructArray = cell(1,length(sList));
% for scenenumber = 1:length(sList)
parfor scenenumber = 1:length(sList)
    cd(tPath)
    sceneN = sList{scenenumber};
    disp(sceneN)
    idScene = sceneN;
    
    trackfile = dir(strcat(eNameKey,'*',sceneN,'*',eName,'.mat'));
    trackfilename = char({trackfile.name});
    
    if ~isempty(trackfilename)
        trackedArray = loadTrackedArray(trackfilename); %trackedArray = Tracked
        plotStruct = plotthemfunctionToStructure(trackedArray,idScene,paStruct,tFrames,sStruct,ctinput);
        plotStructArray{scenenumber} = plotStruct;
    end
end

idx = ~cellfun(@isempty,plotStructArray,'UniformOutput',1);
idxa = find(idx==1);

for scenenumber = idxa
    sceneN = sList{scenenumber};
    plotStruct = plotStructArray{scenenumber};
    
    fnames = fieldnames(plotStruct);
    if isempty(fieldnames(exportStruct)) %if exportStruct is empty
        
        idx = 0;
        for i = 1:length(plotStruct)
            for j = 1:length(fnames)
                exportStruct(idx+i).(fnames{j}) = plotStruct(i).(fnames{j});

            end
            exportStruct(idx+i).scene = sceneN;
            exportStruct(idx+i).cellID = i;
        end
        
    else    %if fields are defined, append the cell data to the next available index
        
        idx = length(exportStruct);
        for i = 1:length(plotStruct)
            for j = 1:length(fnames)
                exportStruct(idx+i).(fnames{j}) = plotStruct(i).(fnames{j});
            end
            exportStruct(idx+i).scene = sceneN;
            exportStruct(idx+i).cellID = i;
        end
        
    end
    
end

%save the exportStruct
cd(dirStruct.exportdir)
filename = strcat(expDetailsStruct.expDateStr,'_tracking_export.mat');
save(filename,'exportStruct','-v6');
end
function trackedArray = loadTrackedArray(trackfilename)
load(trackfilename)
trackedArray = Tracked;
end
function exportSegmentedCells(~,~)
global segStruct pathStruct dirStruct expDetailsStruct SceneList timeFrames


exportNucleiStruct=struct();
%remove all global variables before parfor loop
paStruct = pathStruct;
sStruct = segStruct;
tPath = pathStruct.trackingPath;
sList = SceneList;
tFrames = timeFrames;
cd(tPath)
cd ..

nucleiStructArray = cell(1,length(sList));
% for scenenumber = 1:length(sList)
parfor scenenumber = 1:length(sList)
    cd(tPath)
    sceneN = sList{scenenumber};
    disp(sceneN)
    idScene = sceneN;
    
    exportNucleiz = nucleiToStructure(idScene,paStruct,tFrames,sStruct);
    nucleiStructArray{scenenumber} = exportNucleiz;
end


for scenenumber = 1:length(sList)
    sceneN = sList{scenenumber};
    exportNuclei = nucleiStructArray{scenenumber};
    %             exportNucleiStruct.(sceneN) = exportNuclei;
    
    fnames = fieldnames(exportNuclei);
    if isempty(fieldnames(exportNucleiStruct)) %if exportStruct is empty
        idx = 0;
        for i = 1
            for j = 1:length(fnames)
                exportNucleiStruct(idx+i).(fnames{j}) = {exportNuclei.(fnames{j})};
                exportNucleiStruct(idx+i).scene = sceneN;
                exportNucleiStruct(idx+i).cellID = i;
            end
        end
    else    %if fields are defined, append the cell data to the next available index
        idx = length(exportNucleiStruct);
        for i = 1
            for j = 1:length(fnames)
                exportNucleiStruct(idx+i).(fnames{j}) = {exportNuclei.(fnames{j})};
                exportNucleiStruct(idx+i).scene = sceneN;
                exportNucleiStruct(idx+i).cellID = i;
            end
        end
    end
end

%save the exportStruct
cd(dirStruct.exportdir)
filename = strcat(expDetailsStruct.expDateStr,'_nuclei_export.mat');
save(filename,'exportNucleiStruct');


end
function exportNucleiStruct = nucleiToStructure(idScene,pathStruct,timeFrames,segStruct)
cd(pathStruct.mstackPath)
%no bleach correction option yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   open the image files   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  open smad img  %
cd(pathStruct.mstackPath)
ff = dir(strcat('*',idScene,'*',segStruct.background_seg,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
filename = char(ff.name);
channelfileObject = matfile(filename);
cellQ_imgstack = channelfileObject.flatstack;

%    open cfp img  %
cd(pathStruct.mstackPath)
ff = dir(strcat('*',idScene,'*',segStruct.nucleus_seg,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
filename = char(ff.name);
channelfileObject = matfile(filename);
nuc_imgstack = channelfileObject.flatstack;

% open background Logical img  %
cd(pathStruct.segmentPath)
ff = dir(strcat('*',idScene,'*',segStruct.background_seg,'*'));
if length(ff)>1
    ff = dir(strcat('*',idScene,'*',segStruct.background_seg,'*background*'));
end
%         ff = dir(strcat(ImageDetails.Channel,'*'));
filename = char(ff.name);
channelfileObject = matfile(filename);
bkglogimgstack = channelfileObject.IfFinal;

% open nuclear Logical img  %
cd(pathStruct.segmentPath)
ff = dir(strcat('*',idScene,'*',segStruct.nucleus_seg,'*'));
if length(ff)>1
    ff = dir(strcat('*',idScene,'*',segStruct.nucleus_seg,'*nucleus*'));
end
%         ff = dir(strcat(ImageDetails.Channel,'*'));
filename = char(ff.name);
channelfileObject = matfile(filename);
nuclogimgstack = channelfileObject.IfFinal;
bkglogimgstack(nuclogimgstack) = true;  %make sure that nuclei are not counted in background image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%perform bkg subtraction
cellBKG = zeros(1,timeFrames,'single');
nucBKG = zeros(1,timeFrames,'single');
for k=1:timeFrames
    bkglog = ~bkglogimgstack(:,:,k);
    cellQ_img = single(cellQ_imgstack(:,:,k));
    nuc_img = single(nuc_imgstack(:,:,k));
    
    %background subtraction is just subtraction with a value
    cellBKG(k) = nanmedian(cellQ_img(bkglog));
    nucBKG(k) = nanmedian(nuc_img(bkglog));
    cellQ_imgstack(:,:,k) = cellQ_img-cellBKG(k);
    nuc_imgstack(:,:,k) = nuc_img-nucBKG(k);
    
    %background subtraction is subtraction with an interpolated image
    %     smadbkgimg = regionfill(cellQ_img,~bkglog);
    %     cfpbkgimg = regionfill(nuc_img,~bkglog); %fill in the regions where bkglog is 0
    %     cellQ_imgstack(:,:,k) = cellQ_imgstack(:,:,k)-smadbkgimg;
    %     nuc_imgstack(:,:,k) = nuc_imgstack(:,:,k)-cfpbkgimg;
end


%extract pixel intensities
exportNucleiStruct = struct();
exportNucleiStruct(size(cellQ_imgstack,3)).medianReporter = [];
for i = 1:size(cellQ_imgstack,3)
    cellQ_img = single(squeeze(cellQ_imgstack(:,:,i)));
    nuc_img = single(squeeze(nuc_imgstack(:,:,i)));
    nuclog = squeeze(nuclogimgstack(:,:,i));
    CC = bwconncomp(nuclog);
    PX = CC.PixelIdxList;
    
    imdim = CC.ImageSize;
    ncpx = nuccytopxspitter(PX,imdim);
    nucIntensities = cell(1,length(PX));
    cellIntensities = nucIntensities;
    cytonucIntensities = nucIntensities;
    cytocellIntensities = nucIntensities;
    for j=1:length(PX)
        pxidx = PX{j};
        nxpxidx = ncpx{j};
        nucIntensities{j} = nuc_img(pxidx);
        cellIntensities{j} = cellQ_img(pxidx);
        
        cytonucIntensities{j} = nuc_img(nxpxidx);
        cytocellIntensities{j} = cellQ_img(nxpxidx);
    end
    
    %store centroid values
    stats = regionprops(CC,'Centroid');
    centroidarray = {stats.Centroid};
    Centroid = zeros(2,length(centroidarray));
    for ijk = 1:length(centroidarray)
        Centroid(:,ijk) = centroidarray{ijk};
    end
    
    exportNucleiStruct(i).medianReporter = cellfun(@nanmedian, nucIntensities,'UniformOutput',1);
    exportNucleiStruct(i).medianSmad = cellfun(@nanmedian, cellIntensities,'UniformOutput',1);
    exportNucleiStruct(i).meanReporter = cellfun(@nanmean, nucIntensities,'UniformOutput',1);
    exportNucleiStruct(i).meanSmad = cellfun(@nanmean, cellIntensities,'UniformOutput',1);
    exportNucleiStruct(i).totalReporter = cellfun(@nansum, nucIntensities,'UniformOutput',1);
    exportNucleiStruct(i).totalSmad = cellfun(@nansum, cellIntensities,'UniformOutput',1);
    exportNucleiStruct(i).areaReporter = cellfun(@length, nucIntensities,'UniformOutput',1);
    exportNucleiStruct(i).areaSmad = cellfun(@length, cellIntensities,'UniformOutput',1);
    exportNucleiStruct(i).Centroid = Centroid;
    
    nuccytoRep = zeros(1,length(cytonucIntensities));
    nuccytoSmad = nuccytoRep;
    for j = 1:length(cytonucIntensities)
        nuccytoRep(j) = nanmedian(nucIntensities{j})./nanmedian(cytonucIntensities{j});
        nuccytoSmad(j) = nanmedian(cellIntensities{j})./nanmedian(cytocellIntensities{j});
    end
    
    exportNucleiStruct(i).nuccytoReporter = nuccytoRep;
    exportNucleiStruct(i).nuccytoSmad = nuccytoSmad;
    
end
end
function exportAllCells(~,~)
global ExportNameKey ExportName dirStruct segStruct  pathStruct expDetailsStruct SceneList timeFrames togStruct psettings channelstoinput
exportStruct = struct();

if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
end

%remove all global variables before parfor loop
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
eNameKey = ExportNameKey;
eName = ExportName;
tPath = pathStruct.trackingPath;
mPath = pathStruct.mstackPath;
sList = SceneList;
tFrames = timeFrames;
cSeg = segStruct.background_seg;
nSeg = segStruct.nucleus_seg;
bSeg = segStruct.background_seg;
sPath = pathStruct.segmentPath;
cd(tPath)
cd ..

for scenenumber = 1:length(sList)
    cd(tPath)
    cd ..
    sceneN = char(sList{scenenumber});
    disp(sceneN)
    scenedir = dir(strcat('*',sceneN,'*'));
    scenedirname = char({scenedir.name});
    cd(scenedirname)
    SceneDirPath = char({pwd});
    
    trackfile = dir(strcat(eNameKey,eName,'.mat'));
    trackfilename = char({trackfile.name});
    if ~isempty(trackfilename)
        load(trackfilename)
        PX = trackedArray{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
        makeIMG = zeros(length(framesThatMustBeTracked),length(PX));
        
        makeIMG = makeIMG(1,:)& makeIMG(2,:);
        makeIMG = true(size(makeIMG));
        makeIMGidx = find(makeIMG==1);
        smooththat=0;
        %                 [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,pathStruct.expDirPath,ImageDetails,SceneDirPath,timeFrames,frameToLoad,PlotAxes,ImageDetails.ImgSize,togStruct.plotSettingsToggle,psettings,makeIMG,makeIMGidx,smooththat);
        %                 plotStruct = plotthemfunctionToStructure(trackedArray,idScene,mPath,tFrames,makeIMG,makeIMGidx);
        plotStruct = plotthemfunctionToStructure(Tracked,idScene,pathStruct,timeFrames,makeIMG,makeIMGidx,segStruct,channelstoinput);
        
        fnames = fieldnames(plotStruct);
        if isempty(fieldnames(exportStruct)) %if exportStruct is empty
            idx = 0;
            for i = 1:length(plotStruct)
                for j = 1:length(fnames)
                    exportStruct(idx+i).(fnames{j}) = plotStruct(i).(fnames{j});
                end
                exportStruct(idx+i).scene = sceneN;
                exportStruct(idx+i).cellID = i;
            end
        else    %if fields are defined, append the cell data to the next available index
            idx = length(exportStruct);
            for i = 1:length(plotStruct)
                for j = 1:length(fnames)
                    exportStruct(idx+i).(fnames{j}) = plotStruct(i).(fnames{j});
                end
                exportStruct(idx+i).scene = sceneN;
                exportStruct(idx+i).cellID = i;
            end
        end
    end
end
cd(dirStruct.exportdir)
filename = strcat(expDetailsStruct.expDateStr,'_tracking_allcells_export.mat');
save(filename,'exportStruct');
end


function exportLabels(~,~)
global ExportNameKey ExportName ImageDetails togStruct SceneList pathStruct psettings

if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
end

framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(pathStruct.trackingPath)
cd ..
CENTROID = struct();
for scenenumber = 1:length(SceneList)
    togStruct.updateContrast=1;
    cd(pathStruct.trackingPath)
    cd ..
    sceneN = char(SceneList{scenenumber});
    disp(sceneN)
    ImageDetails.Scene = sceneN;
    scenedir = dir(strcat('*',sceneN,'*'));
    scenedirname = char({scenedir.name});
    cd(scenedirname)
    
    trackfile = dir(strcat(ExportNameKey,ExportName,'.mat'));
    %         trackfile = dir('finalfricktrack.mat');
    trackfilename = char({trackfile.name});
    
    %%% save figure image and centroids %%%
    if ~isempty(trackfilename)
        load(trackfilename)
        setSceneAndTime;
        togStruct.displayTrackingToggle =1;
        finalbutton_callback([],[]);
        xy = labelCells([],[]);
        set(gcf,'Units','Pixels');
        set(gca,'Units','Pixels');
        P = get(gca,'pos');
        F = getframe(gcf,P);
        [X,Map] = frame2im(F);
        set(gcf,'Units','normalized');
        set(gca,'Units','normalized');
        cd(pathStruct.trackingPath)
        cd ..
        cd ..
        cd('ANNOTATIONS')
        imwrite(X,strcat(sceneN,'.jpg'),'JPEG'); %save image
        cd ..
        cd('CENTROIDS')
        
        CENTROIDS.ImageDetails.ImgSize = ImageDetails.ImgSize;
        CENTROIDS.centroids = xy;
        CENTROIDS.scene = sceneN;
        save(strcat('CENTROIDS-',sceneN,'.mat'),'CENTROIDS'); %save CENTROIDS
    end
    %%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
end
togStruct.updateContrast=0;
end


%% tracking functions


function Tracked = makeTrackingFile(timeFrames)

Tracked = struct();
Tracked.trackmatrix = [];
Tracked.arrayStruct = [];

end
function trackbutton_Callback(~,~)
global  pStruct timeSteps Tracked ImageDetails togStruct pathStruct segStruct
pvalue = ImageDetails.Scene;

nucleiDist = pStruct.nucleus.nucDiameter;
tsteps = timeSteps;

tic
trackfilelist = {'yes','no'};
[S,~] = listdlg('PromptString','Are you sure you want to run tracking?',...
    'SelectionMode','single',...
    'ListSize',[200 300],...
    'ListString',trackfilelist);

if S==1
    h=waitbar(0,'running tracking algorithm...');
    Tracked = FrickTrackCellsYeah(pathStruct,pvalue,segStruct,h,nucleiDist,tsteps);
else
end

close(h);
setSceneAndTime;
toc
end

function Tracked = loadTrackedStructure
global pathStruct timeFrames togStruct ImageDetails MainAxes expDetailsStruct trackfilestr ExportName

scenestr = ImageDetails.Scene;
if togStruct.runIterate ==0
    cd(pathStruct.trackingPath)
    trackfile = dir(strcat('*',scenestr,'*num_mmt.mat'));
    if ~isempty(trackfile)
        trackfilelist = {trackfile.name};
        Selection=[];
        [Selection,~] = listdlg('PromptString','Select a tracking file:',...
            'SelectionMode','single',...
            'ListSize',[500 300],...
            'ListString',trackfilelist);
        if ~isempty(Selection)
            load(trackfilelist{Selection}); %load Tracked
            trackfilename = trackfilelist{Selection};
            [a,~] = regexp(trackfilename,ExportName);
            trackfilesub = trackfilename(1:a-1);
            [a,~] = regexp(trackfilesub,'s[0-9]++');
            trackfilestr = trackfilesub(1:a(1)-2);
            
        else
            Tracked = makeTrackingFile(timeFrames);
        end

        
    else
        Tracked = makeTrackingFile(timeFrames);
    end
else
    Tracked = makeTrackingFile(timeFrames);
end

expDetailsStruct.expDateTitleStr = expDetailsStruct.expDateStr;
[a,~] = regexp(expDetailsStruct.expDateTitleStr,'_');expDetailsStruct.expDateTitleStr(a) = '-';

ttl = get(MainAxes,'Title');
set(ttl,'String',[expDetailsStruct.expDateTitleStr ' ' scenestr]);
set(ttl,'FontSize',12);
end
function loadTrackingFile_callback(~,~)
global  Tracked togStruct

Tracked = loadTrackedStructure;
togStruct.trackUpdated = true;
end

function updatetracking_callback(~,~)

        Selection=[];
        trackfilelist = {'yes','no'};
        [Selection,~] = listdlg('PromptString','Are you sure you want to update tracking? this could be slow:',...
            'SelectionMode','single',...
            'ListSize',[500 300],...
            'ListString',trackfilelist);
        if ~isempty(Selection)
            if strcmpi(trackfilelist{Selection},'yes')
                updateFrickTrackCellsYeah
            end
        end
            


end
%make trajectories for overlay of tracking

function trackSaveIterate_callback(~,~)
global togStruct pStruct timeSteps SceneList ImageDetails Tracked pathStruct ExportName timeFrames segStruct

nucleiDist = pStruct.nucleus.nucDiameter;
tsteps = timeSteps;

togStruct.runIterate =1;
togStruct.parallel = 1;
if isempty(h) % Enter parallel loop
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
end

for i=1:length(SceneList)
    %iteratively set scenes
    Trackedz = makeTrackingFile(timeFrames);
    Tracked=Trackedz;
    % Determine the selected data set.
    str = SceneList;
    val = i;
    pvalue = char(str{val});
    ImageDetails.Scene = pvalue;
    setSceneAndTime
    disp(pvalue)
    
    
    %run tracking
    pvalue = ImageDetails.Scene;
    
    Tracked = FrickTrackCellsYeah(pathStruct,pvalue,segStruct,[],nucleiDist,tsteps);
    %             Tracked = FrickTrackCellsYeah(pathStruct.expDirPath,frameToLoad,pvalue,[]);
    setSceneAndTime;
    
    
    %save
    % saveTrackingFileAs_callback([],[])
    cd(pathStruct.trackingPath)
    filename = 'initial';
    save(strcat(filename,'_',ImageDetails.Scene,'_',ExportName,'.mat'),'Tracked')
    %             save(strcat(filename,ExportName,'.mat'),'Tracked')
    
end
togStruct.runIterate =0;
togStruct.parallel = 0;

end
function trackSaveIterateChosen_callback(~,~)
global SceneList pathStruct ExportName segStruct pStruct timeSteps timeFrames

%initialize this structure outside of the parfor loop
nucleiDist = pStruct.nucleus.nucDiameter;
tsteps = timeSteps;

%initialize this structure outside of the parfor loop
psettings = PlotSettings_callback([],[]);
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
framesThatMustBeTracked(2) =  min([timeFrames framesThatMustBeTracked(2)+10]);

%remove all global variables before parfor loop
paStruct = pathStruct;
sStruct = segStruct;
tPath = pathStruct.trackingPath;
sList = SceneList;
eName = ExportName;



possibleWorkers = feature('numcores');
nWorkers = possibleWorkers-2;
%create parallel pool
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool(nWorkers);
else
    nw = poolobj.NumWorkers;
    if ~(nw==nWorkers)
        delete(poolobj)
        poolobj = parpool(nWorkers);
    end
end


%split scenes across workers as evenly as possible
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


%     for iNW=1:nWorkers
parfor iNW=1:nWorkers
    sceneArrayVec = sceneArray{iNW};
    sub_sList = sList(sceneArrayVec);
    for i = 1:length(sceneArrayVec)
        %iteratively set scenes
        scenestr = char(sub_sList{i});
        
        %run tracking
        tf1 = tic;
        Tracked = FrickTrackCellsYeah(paStruct,scenestr,sStruct,[],nucleiDist,tsteps);
        tp1 = toc(tf1);
        
        %display time it takes to track and refine
        
        
        cd(tPath)
        filename = 'tsi';
        saveTrackingIteratively(filename,scenestr,eName,Tracked)
        
        
        %run chosen at two specific frames
        tf2 = tic;
        for ftmbt = framesThatMustBeTracked
            frame = ftmbt;
            Tracked = specialchosenOnesAllOnFrame_Callback(frame,Tracked);
        end
        tp2 = toc(tf2);
        
        disp([scenestr ' - track=' num2str(round(tp1,0,'decimals')) 's - refine=' num2str(round(tp2,1,'decimals')) 's'])
        
        cd(tPath)
        filename = 'tsichosen';
        saveTrackingIteratively(filename,scenestr,eName,Tracked)
    end
end
end
function Tracked = specialchosenOnesAllOnFrame_Callback(frame,Tracked)
trackmatrix = Tracked.trackmatrix;
t = frame;
alltrackframe = trackmatrix(t,:);
idx = find(~isnan(alltrackframe));
if ~isempty(idx)
    cellnum = idx;
    cellmatlog = true(1,size(trackmatrix,2));
    cellmatlog(cellnum) = false;
    trackmatrix(:,cellmatlog) = NaN;
    ucell = 1:size(trackmatrix,2);
    trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
    Tracked.trackmatrix = trackmatrix;
end
end
function saveTrackingIteratively(filename,scenestr,eName,Tracked)
save(strcat(filename,'_',scenestr,'_',eName,'.mat'),'Tracked');
end



%% these functions are for updating segmentation areas/pixels and trackmatrix
function trackmatrix = updateTrackMatrix(trackmatrix,t,dropidx,keepidx,newidx)
% first need to update the cellNumber values so they correspond to the indices of arraystruct;

trackmatrixorig = trackmatrix;
% need to remove tracks from trackmatrix corresponding to dropidx
alltrackframe = trackmatrix(t,:); %pull track idx nums from trackmatrix for given frame
if ~isempty(dropidx)
    
    dropidxorig = dropidx;
    
    if ~isempty(dropidx)
        dropnum = zeros(1,length(dropidx));
        for k = 1:length(dropidx)
            dr = find(dropidx(k) == alltrackframe);
            if ~isempty(dr)
                dropnum(k) = dr;
            else
                %                 disp('fail')
            end
        end
        dropnum(dropnum==0)=[];
    end
    
    %     dropidx(isnan(dropnum)) =[];
    
    
    alltrackframe = trackmatrix(t,:); %pull track idx nums from trackmatrix for given frame
    updatetrackframe = nan(size(alltrackframe)); %make all values = nan
    postdroptrackidx = sort(alltrackframe(~isnan(alltrackframe))); %find all numbers
    
    %determine updated indices
    dropnumz = [];
    for k = dropidx
        dr = find(k == postdroptrackidx);
        if ~isempty(dr)
            dropnumz = horzcat(dropnumz,dr);
        else
            %             disp('fail')
        end
    end
    dropidxz = 1:length(dropnumz);
    dropidx = dropnumz(dropidxz);
    
    postdroptrackidx(dropidx) = []; %remove dropped cells
    if ~isempty(postdroptrackidx)
        newalltrackidx = postdroptrackidx;
        dropidxsub = dropidxorig;
        for k = 1:length(dropidxorig)
            knum = dropidxsub(k);
            newalltrackidx(newalltrackidx>knum) = newalltrackidx(newalltrackidx>knum)-1;
            dropidxsub = dropidxsub-1;
        end
        %         newalltrackidxz = (1:max(postdroptrackidx)-length(dropidx)); %update numbers of non-dropped cells to reflect new numbers in ArrayStruct;
        %need to renumber postdroptrack to be continuous integer set
        for k = 1:length(newalltrackidx)
            allidxnum = find(postdroptrackidx(k) == alltrackframe); %find index where previous numbers matched
            if ~isempty(allidxnum)
                updatetrackframe(allidxnum) = newalltrackidx(k); % if a match is found, update number
            end
        end
    end
    trackmatrix(t,:) = updatetrackframe;
    
    %add the remaining portions of the removed track (beginning first, then ending portion next) to the end of trackmatrix
    newtrackmat = nan(size(trackmatrix,1),length(dropidx)*2);
    newtrackmat(1:(t-1),1:length(dropidx)) = trackmatrix(1:(t-1),dropnum); %first add the beginning portion
    newtrackmat((t+1):size(trackmatrix,1),(length(dropidx)+1):(length(dropidx)*2)) = trackmatrix((t+1):size(trackmatrix,1),dropnum); %next add the subsequent portion
    newtrackmatrix = horzcat(trackmatrix,newtrackmat);
    trackmatrix = newtrackmatrix;
    trackmatrix(:,dropnum) = [];
end



% and need to add tracks to trackmatrix corresponding to newidx
if ~isempty(newidx)
    %     newtrackmatrix = nan(size(trackmatrix,1),size(trackmatrix,2)+length(newidx));
    %     oldtracklog = false(size(newtrackmatrix));
    %     oldtracklog(1:size(trackmatrix,1),1:size(trackmatrix,2)) = true;
    %     newtrackmatrix(oldtracklog) = trackmatrix;
    %     better method
    
    newtrackmat = nan(size(trackmatrix,1),length(newidx));
    updatenewidx = newidx;
    newtrackmat(t,:) = updatenewidx;
    newtrackmatrix = horzcat(trackmatrix,newtrackmat);
    trackmatrix = newtrackmatrix;
end

% ucell = 1:size(trackmatrix,2);
% trackmatrix = refineTrackingMatrix(trackmatrix,ucell);

end
function ArrayStruct = updateArrayStruct(ArrayStruct,PXarray,tvec,dropArray,newArray,keepArray,allArray,trackmatrix,imgsize)
global nucleusimgstack channelimgstack stackstruct segStruct


%update image stack
%background subtract nucleusimgstack
nucleusimgstackz = stackstruct.(segStruct.nucleus_seg);
cellimgstackz = stackstruct.(segStruct.background_seg);



segment_Area_array = ArrayStruct.area;
segment_nucFluor_array = ArrayStruct.nucFluor;
segment_cellFluor_array = ArrayStruct.cellFluor;
segment_Stdev_array = ArrayStruct.stdev;
segment_Centroid_array = ArrayStruct.centroid;
segment_Ellipt_array = ArrayStruct.ellipt;
segment_Pixels_array = ArrayStruct.pixels;
segment_Cellnum_array = ArrayStruct.cellnum;

segment_PixelsCyto_array = ArrayStruct.pixelscyto;
segment_cellCyto_array = ArrayStruct.cellcyto;
segment_nucCyto_array = ArrayStruct.nuccyto;

cnnec = 8;
CC = struct();
CC.ImageSize = [imgsize(1) imgsize(2)];
CC.Connectivity = cnnec;
CC.PixelIdxList =	[];
CC.NumObjects = [];


for i = 1:length(tvec)
    
    dropidx = dropArray{i};
    newidx = newArray{i};
    keepidx = keepArray{i};
    allidx = allArray{i};
    PXnew = PXarray{i};
    tnum = tvec(i);
    nucI = nucleusimgstackz(:,:,tnum);
    cellI = cellimgstackz(:,:,tnum);
    


    %cytoplasmic annulus
    imdim = CC.ImageSize;
    if ~isempty(PXnew)
        ncpx = nuccytopxspitter(PXnew,imdim);
    else
        ncpx = PXnew;
    end
    
    cellnumvec = trackmatrix(tnum,:);
    cellnum = cellnumvec(~isnan(cellnumvec))';
    if ~isempty(PXnew)
        CC.PixelIdxList = PXnew;
        CC.NumObjects = numel(PXnew);
        S = regionprops(CC,'Centroid','Area','Perimeter');
        areavec = vertcat(S.Area);
        perimetervec = vertcat(S.Perimeter);
        elliptvec = 4.*pi.*areavec./(perimetervec.^2);
        centroidMat = vertcat(S.Centroid);
        
        
%         nucFluor = cellfun(@(x) nanmedian(nucI(x)),PXnew,'UniformOutput',1)';
        nucFluor = horzcat(cellfun(@(x) nanmedian(nucI(x)),PXnew,'UniformOutput',1)',cellfun(@(x) nanmean(nucI(x)),PXnew,'UniformOutput',1)',cellfun(@(x) nansum(nucI(x)),PXnew,'UniformOutput',1)');
%         cellFluor = cellfun(@(x) nanmedian(cellI(x)),PXnew,'UniformOutput',1)';
        cellFluor = horzcat(cellfun(@(x) nanmedian(cellI(x)),PXnew,'UniformOutput',1)',cellfun(@(x) nanmean(cellI(x)),PXnew,'UniformOutput',1)',cellfun(@(x) nansum(cellI(x)),PXnew,'UniformOutput',1)');
        stdeval = cellfun(@(x) nanstd(nucI(x)),PXnew,'UniformOutput',1)';
        
        
        cellcyto = cellfun(@(x) nanmedian(cellI(x)),ncpx,'UniformOutput',1);
        nuccyto = cellfun(@(x) nanmedian(nucI(x)),ncpx,'UniformOutput',1);

        
        
        if ~isempty(newidx)
            newnum = ismember(allidx,newidx);
           
            segment_Area_array{tnum} = vertcat(segment_Area_array{tnum},areavec(newnum));
            segment_Centroid_array{tnum} = vertcat(segment_Centroid_array{tnum},centroidMat(newnum,:));
            segment_Pixels_array{tnum} = horzcat(segment_Pixels_array{tnum},PXnew(newnum));
            segment_Cellnum_array{tnum} = cellnum;
            segment_Ellipt_array{tnum} = vertcat(segment_Ellipt_array{tnum},elliptvec(newnum));
            segment_nucFluor_array{tnum} = vertcat(segment_nucFluor_array{tnum},nucFluor(newnum,:));
            segment_cellFluor_array{tnum} = vertcat(segment_cellFluor_array{tnum},cellFluor(newnum,:));
            segment_Stdev_array{tnum} = vertcat(segment_Stdev_array{tnum},stdeval(newnum));
            
            segment_PixelsCyto_array{tnum} = horzcat(segment_PixelsCyto_array{tnum},ncpx(newnum));
            segment_cellCyto_array{tnum} = horzcat(segment_cellCyto_array{tnum},cellcyto(newnum));
            segment_nucCyto_array{tnum} = horzcat(segment_nucCyto_array{tnum},nuccyto(newnum));
            


        end
        
        if ~isempty(keepidx)
            newnum = ismember(allidx,keepidx);
            
            segment_Area_array{tnum}(keepidx) = areavec(newnum);
            segment_Centroid_array{tnum}(keepidx,true(1,size(centroidMat,2))) = centroidMat(newnum,:);
            segment_Pixels_array{tnum}(keepidx) = PXnew(newnum);
            segment_Cellnum_array{tnum}= cellnum;
            segment_Ellipt_array{tnum}(keepidx) = elliptvec(newnum);
            segment_nucFluor_array{tnum}(keepidx,true(1,size(nucFluor,2))) = nucFluor(newnum,:);
            segment_cellFluor_array{tnum}(keepidx,true(1,size(cellFluor,2)))= cellFluor(newnum,:);
            segment_Stdev_array{tnum}(keepidx) = stdeval(newnum);
            
            segment_PixelsCyto_array{tnum}(keepidx) = ncpx(newnum);
            segment_cellCyto_array{tnum}(keepidx) = cellcyto(newnum);
            segment_nucCyto_array{tnum}(keepidx) = nuccyto(newnum); 
        end
    end
    
    if ~isempty(dropidx)
        centold = segment_Centroid_array{tnum};
        centidx = 1:size(centold,1);
        centidx(dropidx)=[];
        centnew = centold(centidx,:);
        
        nucFluorold = segment_nucFluor_array{tnum};
        nucFluornew = nucFluorold(centidx,:);
        
        cellFluorold = segment_cellFluor_array{tnum};
        cellFluornew = cellFluorold(centidx,:);
        
        segment_Area_array{tnum}(dropidx) = [];
        segment_Centroid_array{tnum} = centnew;
        segment_Pixels_array{tnum}(dropidx) = [];
        segment_Cellnum_array{tnum}= cellnum;
        segment_Ellipt_array{tnum}(dropidx) = [];
        segment_nucFluor_array{tnum} = nucFluornew;
        segment_cellFluor_array{tnum} = cellFluornew;
        segment_Stdev_array{tnum}(dropidx) = [];
        
        
        segment_PixelsCyto_array{tnum}(dropidx) = [];
        segment_cellCyto_array{tnum}(dropidx) = [];
        segment_nucCyto_array{tnum}(dropidx) = [];
    end
end
ArrayStruct.area    =   segment_Area_array;
ArrayStruct.ellipt     =   segment_Ellipt_array;
ArrayStruct.nucFluor     =   segment_nucFluor_array;
ArrayStruct.cellFluor     =   segment_cellFluor_array;
ArrayStruct.stdev     =   segment_Stdev_array;
ArrayStruct.centroid     =   segment_Centroid_array;
ArrayStruct.pixels     =   segment_Pixels_array;
ArrayStruct.cellnum     =   segment_Cellnum_array;

ArrayStruct.pixelscyto = segment_PixelsCyto_array;
ArrayStruct.cellcyto = segment_cellCyto_array;
ArrayStruct.nuccyto = segment_nucCyto_array;
end

%% functions for tracking cells
function updateFrickTrackCellsYeah
global pathStruct segStruct SceneList ImageDetails


for snum = 1:length(SceneList)
% for snum = 1
disp(snum)
    pvalue = SceneList{snum};
%     pvalue = ImageDetails.Scene;
    cd(pathStruct.trackingPath)
    trackfile = dir(strcat('*',pvalue,'*_num_mmt.mat'));
    trackfile = dir(strcat('*final*',pvalue,'*_num_mmt.mat'));
    tfilelist = {trackfile.name};
    if ~isempty(trackfile)
        
        
        
        %% load the images
        
        %set details for opening images necessary for tracking
        imgstruct = struct();
        imgstruct(1).path = pathStruct.segmentPath;
        imgstruct(2).path = pathStruct.segmentPath;
        imgstruct(3).path = pathStruct.mstackPath;
        imgstruct(4).path = pathStruct.mstackPath;
        
        imgstruct(1).seginstruct = segStruct.nucleus_seg;
        imgstruct(2).seginstruct = segStruct.background_seg;
        imgstruct(3).seginstruct = segStruct.nucleus_seg;
        imgstruct(4).seginstruct = segStruct.background_seg;
        
        imgstruct(1).segstr = 'nucleus';
        imgstruct(2).segstr = 'background';
        imgstruct(3).segstr = 'nucleus';
        imgstruct(4).segstr = 'cell';
        
        imgstruct(1).openstr = 'IfFinal';
        imgstruct(2).openstr = 'IfFinal';
        imgstruct(3).openstr = 'flatstack';
        imgstruct(4).openstr = 'flatstack';
        
        imgstruct(1).image =[];
        imgstruct(1).filename =[];
        
        
        
        %load nucleus segmentation binary as segmentimgstack
        for i = 1:length(imgstruct) %parfor is slower with only 4 images
            pathstr = imgstruct(i).path;
            seginstruct = imgstruct(i).seginstruct;
            segstr = imgstruct(i).segstr;
            openstr = imgstruct(i).openstr;
            cd(pathstr)
            ff = dir(strcat('*',pvalue,'*',seginstruct,'*'));
            if length(ff)>1
                ff = dir(strcat('*',pvalue,'*',seginstruct,'*',segstr,'*'));
            end
            fimgname = char(ff.name);
            nucleusFileObject = matfile(fimgname);
            indimage = nucleusFileObject.(openstr);
            imgstruct(i).image = indimage;
            imgstruct(i).filename = fimgname;
        end
        
        segmentimgstack = imgstruct(1).image;
        backgroundimgstack = imgstruct(2).image;
        nosub_nucleusimgstack = imgstruct(3).image;
        nosub_cellimgstack = imgstruct(4).image;
        filename = imgstruct(1).filename;
        
        backgroundimgstack(segmentimgstack) = true; %update backgruond segmentation to exclude nuclei segmentation
        
        %background subtract nucleusimgstack
        nucleusimgstack = zeros(size(nosub_nucleusimgstack));
        cellimgstack = nucleusimgstack;
        for i=1:size(nosub_nucleusimgstack,3)
            nucI = nosub_nucleusimgstack(:,:,i);
            cellI = nosub_cellimgstack(:,:,i);
            bkgI = backgroundimgstack(:,:,i);
            nucleusimgstack(:,:,i) = nucI-nanmedian(nucI(~bkgI));
            cellimgstack(:,:,i) = cellI - nanmedian(cellI(~bkgI));
        end
        
        
        
        %% load previous tracking
        for tnum = 1:length(tfilelist)
            tfilestr = tfilelist{tnum};
            disp(pvalue)
            cd(pathStruct.trackingPath)
            load(tfilestr); %load Tracked
            
            % ArrayStruct run through based on track'
            
            %determine parameters of nuclei (%area,centroid,fluorescence,elllipticity,velocities?)
            ArrayStruct = Tracked.arrayStruct;
            segment_Area_array = ArrayStruct.area;
            segment_nucFluor_array = ArrayStruct.nucFluor;
            segment_cellFluor_array = ArrayStruct.cellFluor;
            segment_Stdev_array = ArrayStruct.stdev;
            segment_Centroid_array = ArrayStruct.centroid;
            segment_Ellipt_array = ArrayStruct.ellipt;
            segment_Pixels_array = ArrayStruct.pixels;
            segment_Cellnum_array = ArrayStruct.cellnum;
            
            
            fnames = fieldnames(ArrayStruct);
            if sum(strcmp(fnames,'pixelscyto'))>0
                segment_PixelsCyto_array = ArrayStruct.pixelscyto;
                segment_cellCyto_array = ArrayStruct.cellcyto;
                segment_nucCyto_array = ArrayStruct.nuccyto;
            end

            segmentsequence = 1:size(segmentimgstack,3); %track from first frame to last frame
            for i = segmentsequence
                segI = segmentimgstack(:,:,i);
                nucI = nucleusimgstack(:,:,i);
                cellI = cellimgstack(:,:,i);
                CC = bwconncomp(segI);
                %
                PX = segment_Pixels_array{i}; %use the previous pixels
                CC.PixelIdxList = PX;
                CC.NumObjects = length(PX);
                imdim = CC.ImageSize;
                %
                if ~sum(strcmp(fnames,'pixelscyto'))>0
                    ncpx = nuccytopxspitter(PX,imdim);
                else
                    ncpx = segment_PixelsCyto_array{i};
                end
                
                
                S = regionprops(CC,'Centroid','Area','Perimeter');
                areavec = vertcat(S.Area);
                perimetervec = vertcat(S.Perimeter);
                elliptvec = 4.*pi.*areavec./(perimetervec.^2);
                centroidMat = vertcat(S.Centroid);
                
                segment_Area_array{i} = areavec;
                segment_Ellipt_array{i} = elliptvec;
                %     segment_nucFluor_array{i} = cellfun(@(x) nanmedian(nucI(x)),PX,'UniformOutput',1)';
                segment_nucFluor_array{i} = horzcat(cellfun(@(x) nanmedian(nucI(x)),PX,'UniformOutput',1)',cellfun(@(x) nanmean(nucI(x)),PX,'UniformOutput',1)',cellfun(@(x) nansum(nucI(x)),PX,'UniformOutput',1)');
                segment_cellFluor_array{i} = horzcat(cellfun(@(x) nanmedian(cellI(x)),PX,'UniformOutput',1)',cellfun(@(x) nanmean(cellI(x)),PX,'UniformOutput',1)',cellfun(@(x) nansum(cellI(x)),PX,'UniformOutput',1)');
                %     segment_cellFluor_array{i} = cellfun(@(x) nanmedian(cellI(x)),PX,'UniformOutput',1)';
                segment_Stdev_array{i} = cellfun(@(x) nanstd(nucI(x)),PX,'UniformOutput',1)';
                segment_Centroid_array{i} = centroidMat;
                segment_Pixels_array{i} = PX;
                segment_Cellnum_array{i} = (1:length(PX))';
                
                segment_PixelsCyto_array{i} = ncpx;
                segment_cellCyto_array{i} = cellfun(@(x) nanmedian(cellI(x)),ncpx,'UniformOutput',1);
                segment_nucCyto_array{i} = cellfun(@(x) nanmedian(nucI(x)),ncpx,'UniformOutput',1);
                
                
            end
            
            ArrayStruct.area     =   segment_Area_array;
            ArrayStruct.ellipt     =   segment_Ellipt_array;
            ArrayStruct.nucFluor     =   segment_nucFluor_array;
            ArrayStruct.cellFluor     =   segment_cellFluor_array;
            ArrayStruct.stdev     =   segment_Stdev_array;
            ArrayStruct.centroid     =   segment_Centroid_array;
            ArrayStruct.pixels     =   segment_Pixels_array;
            ArrayStruct.cellnum     =   segment_Cellnum_array;
            
            ArrayStruct.pixelscyto = segment_PixelsCyto_array;
            ArrayStruct.cellcyto = segment_cellCyto_array;
            ArrayStruct.nuccyto = segment_nucCyto_array;
            Tracked.arrayStruct = ArrayStruct;
            
            
            
            %             [~,b] = regexp(tfilestr,'fricktrack');
            %             newtfilestr = [tfilestr(1:b) '_num_mmt.mat'];
            newtfilestr = tfilestr;
            cd(pathStruct.trackingPath)
            save(newtfilestr,'Tracked')
        end
    end
end
end

function [ Tracked ] = FrickTrackCellsYeah(pathStruct,pvalue,segStruct,h,nucleiDist,tsteps)
%% intialize
%set details for opening images necessary for tracking
imgstruct = struct();
imgstruct(1).path = pathStruct.segmentPath;
imgstruct(2).path = pathStruct.segmentPath;
imgstruct(3).path = pathStruct.mstackPath;
imgstruct(4).path = pathStruct.mstackPath;

imgstruct(1).seginstruct = segStruct.nucleus_seg;
imgstruct(2).seginstruct = segStruct.background_seg;
imgstruct(3).seginstruct = segStruct.nucleus_seg;
imgstruct(4).seginstruct = segStruct.background_seg;

imgstruct(1).segstr = 'nucleus';
imgstruct(2).segstr = 'background';
imgstruct(3).segstr = 'nucleus';
imgstruct(4).segstr = 'cell';

imgstruct(1).openstr = 'IfFinal';
imgstruct(2).openstr = 'IfFinal';
imgstruct(3).openstr = 'flatstack';
imgstruct(4).openstr = 'flatstack';

imgstruct(1).image =[];
imgstruct(1).filename =[];



%load nucleus segmentation binary as segmentimgstack
for i = 1:length(imgstruct) %parfor is slower with only 4 images
    pathstr = imgstruct(i).path;
    seginstruct = imgstruct(i).seginstruct;
    segstr = imgstruct(i).segstr;
    openstr = imgstruct(i).openstr;
    cd(pathstr)
    ff = dir(strcat('*',pvalue,'*',seginstruct,'*'));
    if length(ff)>1
        ff = dir(strcat('*',pvalue,'*',seginstruct,'*',segstr,'*'));
    end
    fimgname = char(ff.name);
    nucleusFileObject = matfile(fimgname);
    indimage = nucleusFileObject.(openstr);
    imgstruct(i).image = indimage;
    imgstruct(i).filename = fimgname;
end

segmentimgstack = imgstruct(1).image;
backgroundimgstack = imgstruct(2).image;
nosub_nucleusimgstack = imgstruct(3).image;
nosub_cellimgstack = imgstruct(4).image;
filename = imgstruct(1).filename;

backgroundimgstack(segmentimgstack) = true; %update backgruond segmentation to exclude nuclei segmentation

%background subtract nucleusimgstack
nucleusimgstack = zeros(size(nosub_nucleusimgstack));
cellimgstack = nucleusimgstack;
for i=1:size(nosub_nucleusimgstack,3)
    nucI = nosub_nucleusimgstack(:,:,i);
    cellI = nosub_cellimgstack(:,:,i);
    bkgI = backgroundimgstack(:,:,i);
    nucleusimgstack(:,:,i) = nucI-nanmedian(nucI(~bkgI));
    cellimgstack(:,:,i) = cellI - nanmedian(cellI(~bkgI));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section involves calculating the nearest neighbor to the centroid
%and organizing PixelLists and Centroid lists to match nearest neighbor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine parameters of nuclei (%area,centroid,fluorescence,elllipticity,velocities?)
segment_Area_array = cell(1,size(segmentimgstack,3));
segment_nucFluor_array = cell(1,size(segmentimgstack,3));
segment_cellFluor_array = cell(1,size(segmentimgstack,3));
segment_Stdev_array = cell(1,size(segmentimgstack,3));
segment_Centroid_array = cell(1,size(segmentimgstack,3));
segment_Ellipt_array = cell(1,size(segmentimgstack,3));
segment_Pixels_array = cell(1,size(segmentimgstack,3));
segment_Cellnum_array = cell(1,size(segmentimgstack,3));

segment_PixelsCyto_array{i} = cell(1,size(segmentimgstack,3));
segment_cellCyto_array{i} = cell(1,size(segmentimgstack,3));
segment_nucCyto_array{i} = cell(1,size(segmentimgstack,3));
segmentsequence = 1:size(segmentimgstack,3); %track from first frame to last frame
for i = segmentsequence
    segI = segmentimgstack(:,:,i);
    nucI = nucleusimgstack(:,:,i);
    cellI = cellimgstack(:,:,i);
    CC = bwconncomp(segI);
    PX = CC.PixelIdxList;
    


    
    S = regionprops(CC,'Centroid','Area','Perimeter');
    areavec = vertcat(S.Area);
    perimetervec = vertcat(S.Perimeter);
    elliptvec = 4.*pi.*areavec./(perimetervec.^2);
    centroidMat = vertcat(S.Centroid);
    
    segment_Area_array{i} = areavec;
    segment_Ellipt_array{i} = elliptvec;
%     segment_nucFluor_array{i} = cellfun(@(x) nanmedian(nucI(x)),PX,'UniformOutput',1)';
    segment_nucFluor_array{i} = horzcat(cellfun(@(x) nanmedian(nucI(x)),PX,'UniformOutput',1)',cellfun(@(x) nanmean(nucI(x)),PX,'UniformOutput',1)',cellfun(@(x) nansum(nucI(x)),PX,'UniformOutput',1)');
    segment_cellFluor_array{i} = horzcat(cellfun(@(x) nanmedian(cellI(x)),PX,'UniformOutput',1)',cellfun(@(x) nanmean(cellI(x)),PX,'UniformOutput',1)',cellfun(@(x) nansum(cellI(x)),PX,'UniformOutput',1)');
%     segment_cellFluor_array{i} = cellfun(@(x) nanmedian(cellI(x)),PX,'UniformOutput',1)';
    segment_Stdev_array{i} = cellfun(@(x) nanstd(nucI(x)),PX,'UniformOutput',1)';
    segment_Centroid_array{i} = centroidMat;
    segment_Pixels_array{i} = PX;
    segment_Cellnum_array{i} = (1:length(PX))';
    
    
    %cytoplasmic annulus
    imdim = CC.ImageSize;
    if ~isempty(PX)
        ncpx = nuccytopxspitter(PX,imdim);
    else
        ncpx = PX;
    end
    segment_PixelsCyto_array{i} = ncpx;
    segment_cellCyto_array{i} = cellfun(@(x) nanmedian(cellI(x)),ncpx,'UniformOutput',1);
    segment_nucCyto_array{i} = cellfun(@(x) nanmedian(nucI(x)),ncpx,'UniformOutput',1);
    
end

ArrayStruct = struct();
ArrayStruct.area     =   segment_Area_array;
ArrayStruct.ellipt     =   segment_Ellipt_array;
ArrayStruct.nucFluor     =   segment_nucFluor_array;
ArrayStruct.cellFluor     =   segment_cellFluor_array;
ArrayStruct.stdev     =   segment_Stdev_array;
ArrayStruct.centroid     =   segment_Centroid_array;
ArrayStruct.pixels     =   segment_Pixels_array;
ArrayStruct.cellnum     =   segment_Cellnum_array;

ArrayStruct.pixelscyto = segment_PixelsCyto_array;
ArrayStruct.cellcyto = segment_cellCyto_array;
ArrayStruct.nuccyto = segment_nucCyto_array;
%         ArrayStruct.If  = IfPerimFunction();


%
%% tracking

%track based on assigning probabilities determined by minimizing changes to measured parameters
segmentsequence = 2:size(segmentimgstack,3); %track from first frame to last frame
birtharray = cell(1,size(segmentimgstack,3));
nidxarray = cell(1,size(segmentimgstack,3));
birtharray(1) = segment_Cellnum_array(1);


possibleWorkers = feature('numcores');
nWorkers = possibleWorkers;

%split scenes across workers as evenly as possible
sceneVector = segmentsequence;
fractions = length(sceneVector)./nWorkers;
firstlengths = ceil(fractions);
if firstlengths*nWorkers > length(segmentsequence)
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

birtharrayz = cell(1,nWorkers);
nidxarrayz = cell(1,nWorkers);

for iNW = 1:nWorkers
    sceneArrayVec = sceneArray{iNW};
    birtharraypre = cell(size(sceneArrayVec));
    nidxarraypre = cell(size(sceneArrayVec));
    for ii = 1:length(sceneArrayVec)
        i = sceneArrayVec(ii);
        if ~isempty(h)
            waitbar(i./max(segmentsequence),h)
        end
        
        

        a=i-1;
        b=i;
        
        xcorrim0 = min(1,segmentimgstack(:,:,a));
        xcorrim1 = min(1,segmentimgstack(:,:,b));
         nres=4;
        if max(max(xcorrim0))==0 || max(max(xcorrim1))==0
            %no cells were segmented
            regshift=[0,0];
        else
            xc=xcorr2(xcorrim0(1:nres:end,1:nres:end),xcorrim1(1:nres:end,1:nres:end)); %compute low resolution cross correlation
            [I,J]=find(xc==max(xc(:)),1);
            regshift=nres*[I,J]-[imdim(2),imdim(1)]; %determine how far the maxima in the image moves, i think by computing cross correlation      
        end
        % if regshift is small just ignore it
        if max(abs(regshift))<1
            regshift=[0,0];
        end
        
        %%

        centroids = segment_Centroid_array{b};
        centroidsPrev = segment_Centroid_array{a};
        area = segment_Area_array{b};
        areaPrev = segment_Area_array{a};
        nucfluor = segment_nucFluor_array{b};
        nucfluorPrev = segment_nucFluor_array{a};
        cellfluor = segment_cellFluor_array{b};
        cellfluorPrev = segment_cellFluor_array{a};
        
        cytocellfluor = segment_cellCyto_array{b};
        cytocellfluorPrev = segment_cellCyto_array{a};
        cytonucfluor = segment_nucCyto_array{b};
        cytonucfluorPrev = segment_nucCyto_array{a};
%         nucfluor = segment_nucFluor_array{b}(:,1); %choose median
%         nucfluorPrev = segment_nucFluor_array{a}(:,1); %choose median
%         cellfluor = segment_cellFluor_array{b}(:,1); %choose median
%         cellfluorPrev = segment_cellFluor_array{a}(:,1); %choose median

        pixels = segment_Pixels_array{b};
        pixelsPrev = segment_Pixels_array{a};
        
        displacementCutoff = (nucleiDist)*(tsteps(i-1)./10);
        
        if ~(isempty(pixels) || isempty(pixelsPrev))
            %nearest neighbor distances to determine probability that
            %two cells are the same
            knnnum = 5;
            
            areachoice = 1.2; %1.2
            displacementCutoff = (nucleiDist)*(tsteps(i-1)./4); %(nucleiDist)*(tsteps(i-1)./10);
            nucFluorchoice = 3; %1.2
            cellFluorchoice = 3; %1.2
            
            [distProb,~,distcut] = probSpitter(centroids-regshift,centroidsPrev,knnnum,displacementCutoff,[]);

            %added multiFactor to improve track lengths and accuracy
            multiFactorPrev = [nucfluorPrev(:,1) cellfluorPrev(:,1) areaPrev cytocellfluorPrev' cytonucfluorPrev'];
            multiFactor = [nucfluor(:,1) cellfluor(:,1) area cytocellfluor' cytonucfluor'];
            [distProb2,~,distcut2] = probSpitter(multiFactor,multiFactorPrev,knnnum,displacementCutoff,[]);
            
            [areaProbette,~,areacut] = probSpitter(area,areaPrev,size(area,1),[],areachoice);
            [nucFluorProbette,~,nuccut] = probSpitter(nucfluor(:,1),nucfluorPrev(:,1),size(nucfluor,1),[],nucFluorchoice); %choose to track based on median
            [cellFluorProbette,~,cellcut] = probSpitter(cellfluor(:,1),cellfluorPrev(:,1),size(cellfluor,1),[],cellFluorchoice); %choose to track based on median
%             newProb = distProb.*areaProbette.*nucFluorProbette.*cellFluorProbette;
            newProb = distProb.*distProb2;
            %distProb dimesionas are length(centroidsPrev) x length(centroids)
            
            
            %centroids [37x2] centroidPrev[34x2] idx[34x3] distProb[34x37];
            trackProb = distProb;
            trackProb(distcut|nuccut|areacut)=0;
%             trackProb(distcut)=0;
%             trackProb(distcut|nuccut|areacut|cellcut)=0;
            [maxvals,idx] = max(trackProb,[],2); %idx is index of input that matches to inputPrev such that input(idx) = inputPrev;
            %idx has length(idx) = length(inputPrev);
            %idx has max(idx) = length(input);
            
            
            %now some cells are likely to be assigned twice. Correct this based on highest probabilities
            num_cells_currentFrame = 1:size(trackProb,2); %size(trackProb,2) = length(input); %num_cells_set should be = the size of the current NOT the prev
            %num_cells_prevFrame = 1:size(trackProb,1); %size(trackProb,1) = length(inputPrev)
            
            arbitrateProb = newProb;
            arbitrateProb(distcut|nuccut|areacut|cellcut)=0;
            [arbmaxvals,~] = max(arbitrateProb,[],2); %idx is index of input that matches to inputPrev such that input(idx) = inputPrev;
            arbmaxvals(arbmaxvals==0)=NaN;
            idx(maxvals==0) = NaN;
            idx(isnan(maxvals)) = NaN;
            
            [n, bin] = histc(idx, num_cells_currentFrame);
            multiple = find(n>1);%the same cell is called closest to two previous cells
            missers = find(n<1);
            loserz = [];
            if ~isempty(multiple)
                for loop = multiple' %loop is the cellID in current frame
                    testidx    = find(ismember(bin, loop)); %cell IDs from prev frame
                    winneridx = find(arbmaxvals(testidx) == max(arbmaxvals(testidx))); %find the cell with the highest probability match
                    winneridx(isnan(arbmaxvals(winneridx)))=[];%if the match is an artifact (if probability=NaN, then remove)
                    if ~isempty(winneridx)
                        testidx(winneridx(1))=[]; %must be winneridx(1) because only one cell can win
                    end
                    lidx=testidx;
                    loseridx = find(true(size(testidx))==1);
                    
                    for ll = 1:length(lidx)
                        lval = lidx(ll);
                        arbval = zeros(1,length(missers));
                        for u = 1:length(missers)
                            arbval(u) = arbitrateProb(lval,missers(u));
                        end
                        [maxarbval,maxidx] = max(arbval);
                        if maxarbval>0
                            idx(lval) = missers(maxidx);
                            missers(maxidx)=[];
                            loseridx(ll)=NaN;
                        end
                    end
                    
                    loseridx(isnan(loseridx))=[];
                    lidx=testidx(loseridx);
                    losern = lidx;
                    loserz = [loserz(:)' losern(:)'];
                end
            end
            %the index of loserz represents the cell number of
            %previous frame and the new indexed numbers of current
            %frame
            %i.e. loserz =3 means cell 3 from previous frame
            %matches with index 3 of newidx
            nidx = idx;
            nidx(loserz)=NaN; %remove duplicates (loserz)// that is to say, tracks that merge onto one cell
            
            [nn, ~] = histc([nidx(:)' missers(:)'], num_cells_currentFrame);
            updatemissers = find(nn<1);
            newidx = [nidx(:)' missers(:)' updatemissers(:)']';
            
        elseif isempty(pixelsPrev) %no cells in the previous frame
            idx = nan(1);
            nidx = idx;
            missers = [];
            num_cells_currentFrame = length(area);
            [nn, ~] = histc([nidx(:)' missers(:)'], num_cells_currentFrame);
            updatemissers = find(nn<1);
            newidx = [nidx(:)' missers(:)' updatemissers(:)']';
        else %if no cells are segmented on the current frame
            idx = nan(size(areaPrev));
            nidx = idx;
            newidx = idx;
            missers = [];
            updatemissers = [];
        end
        
        [nn, ~] = histc(newidx, 1:nanmax(newidx));
        if max(nn)>1
            error('cells called twice!!!')
        end
        
        birtharraypre(ii) = {[missers(:)' updatemissers(:)']'};
        nidxarraypre(ii) = {nidx};
    end
    birtharrayz{iNW} = birtharraypre;
    nidxarrayz{iNW} = nidxarraypre;
end
birtharray(2:end) = horzcat(birtharrayz{:});
nidxarray(2:end) = horzcat(nidxarrayz{:}) ;
%THIS ENDS THE TRACKING PORTION

%% assemby portion

%determine the number of births (new tracks) to determine the dimensions of the track vector
numbirths = zeros(1,length(birtharray));
birthdistmat = zeros(1,length(birtharray));
for i = 1:length(nidxarray)
    births = birtharray{i};
    numbirths(i) = length(births);
    if i>1
        birthdistmat(i) = sum(numbirths(1:i-1));
    else
        birthdistmat(i) = 0;
    end
end

%now move tracks so that they each involve their own row in a matrix
number_of_tracks = sum(numbirths);
trackmatrix = nan(size(segmentimgstack,3),number_of_tracks);
for i = 1:size(segmentimgstack,3)
    births = birtharray{i};
    nidx = nidxarray{i}';
    nidvec = nidx(~isnan(nidx));
    bdist = birthdistmat(i);
    nb = numbirths(i);
    
    if i>1
        prevvecinit = trackmatrix(i-1,1:bdist);
        prevvec = prevvecinit(~isnan(prevvecinit));
        %they are always the same size
        %        disp(size(prevvec))
        %        disp(size(nidx))
        %        disp(max(prevvec))
    else
        prevvec = 1:length(nidvec);
    end
    
    if ~isempty(nidx(~isnan(nidx)))
        jvec = 1:length(prevvec);
        jvectwo = zeros(size(jvec));
        for k = 1:length(prevvec)
            nval = find(prevvecinit == prevvec(k));
            jvectwo(k) = nval;
        end
        
        for j = jvec
            trackmatrix(i,jvectwo(j)) = nidx(prevvec(j));
        end
    end
    
    if ~isempty(births)
        trackmatrix(i,bdist+1:bdist+nb) = births;
    end
end

ucell = 1:size(trackmatrix,2);
trackmatrix = refineTrackingMatrix(trackmatrix,ucell);
Tracked = struct();
Tracked.trackmatrix = trackmatrix;
Tracked.arrayStruct = ArrayStruct;
end



function [distProb,idx,cutoffidx] = probSpitter(input,inputPrev,knnnum,displacementCutoff,relativeCutoff)
%make a matrix that tells you the probability a value in one vector is the same as anothe





if ~isempty(displacementCutoff)
    [idx,eps] = knnsearch(input,inputPrev,'K',knnnum);
    eps_scaled = false(size(eps));
elseif ~isempty(relativeCutoff)
    [idx,eps_scaled] = knnsearch(input,inputPrev,'K',knnnum,'Distance',@distfun);
    eps = eps_scaled;
% %     [idx,eps_scaled] = knnsearch(input,inputPrev,'K',knnnum,'Distance','seuclidean');
% %     valPrev = input(idx);
% %     epsmin = min(inputPrev,valPrev);
% %     epsmax = max(inputPrev,valPrev);
% %     eps_scaled = epsmax./epsmin;
%     
%     valCurr = input(idx);
%     valPrev = nan(size(inputPrev,1),knnnum);
%     for ki = 1:knnnum
%         valPrev(:,ki) = inputPrev(:,1);
%     end
%     epsmin = min(valPrev,valCurr);
%     epsmax = max(valPrev,valCurr);
%     %     epsmax2 = max(input,valPrev')';
%     eps_scaled = epsmax./epsmin;
end
distvec = eps;
distsub = 1./(distvec./max(distvec(:)));
%     distsub = distvec./max(distvec(:));
%     distProbValues = 1-distsub;
if (sum(distsub(:))==0) || isnan(nansum(distsub(:)))
    distsuber = 1./(distvec);
    distsuber(distsuber==Inf)=1;
    distsub = distsuber./(max(distsuber(:)));
    disp('one cell')
end
distsub(distsub==Inf) = max(distsub(~(distsub==Inf)));
distProbValues = distsub./max(distsub(:));

distProb = zeros(size(inputPrev,1),size(input,1));
epsMat = distProb;
for di = 1:size(idx,1)
    didx = idx(di,:);
    distProb(di,didx) = distProbValues(di,:);
    if ~isempty(displacementCutoff)
        epsMat(di,didx) = eps(di,:);
    elseif ~isempty(relativeCutoff)
        epsMat(di,didx) = eps_scaled(di,:);
    end
end

if ~isempty(displacementCutoff)
    cutoffidx = epsMat>displacementCutoff;
elseif ~isempty(relativeCutoff)
    cutoffidx = epsMat>relativeCutoff;
end

end
function D2 = distfun(ZI,ZJ)
% xs = ZI;
% yt = ZJ;
% 
%ZJ can be a vector

% disp(ZI)
% disp(ZJ)
x1 = min(ZI(:)',ZJ(:)');
y1 = max(ZI(:)',ZJ(:)');

% D2 = max([ZI ZJ])./min([ZI ZJ]);
D2 = (y1./x1)';
% disp(D2)
% D2 = 1;
end
function saveTrackingFileAs_callback(~,~)
global  pathStruct Tracked ExportName ImageDetails trackfilestr
cd(pathStruct.trackingPath)
prompt = 'filename of tracking structure to be saved?';
dlg_title = 'save tracking structure as...specific filename';
filename = char(inputdlg(prompt,dlg_title));
trackfilename = strcat(filename,'_',ImageDetails.Scene,'_',ExportName,'.mat');
trackfilestr = filename;
save(trackfilename,'Tracked')
end

%% Image Display functions
function setSceneAndTime
global  channelstoinput stackstruct bkgmedianmat Ifstack togStruct foStruct  segStruct  pathStruct frameToLoad ImageDetails Tracked SceneList imgfile
trackmatrix = Tracked.trackmatrix;

%check for empty variables
% cd(pathStruct.mstackPath)
if isempty(ImageDetails.Scene)
    ImageDetails.Scene = SceneList{1};
end
if isempty(ImageDetails.Channel)
    ImageDetails.Channel = segStruct.nucleus_seg;
end

%determine the frame to load
if isempty(ImageDetails.Frame)
    ImageDetails.Frame = frameToLoad;
    t = ImageDetails.Frame;
else
    t = ImageDetails.Frame;
end




%% iterate through all channels
% stackstruct = struct(); %place up top and make global
for ic = 1:length(channelstoinput)
    cstr = channelstoinput{ic};
    imgstack = stackstruct.(cstr);
    %choose the channel image
    if ~isempty(foStruct.(cstr)) %if channelfileObject has been made, check to see if the scene has changed.
        [a,~] = regexp(foStruct.(cstr),ImageDetails.Scene);
        if isempty(a) %if the scene has changed load the new channelimgstack
            cd(pathStruct.mstackPath)
            ff = dir(strcat('*',ImageDetails.Scene,'*',cstr,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            imgstack = channelfileObject.flatstack;
            foStruct.(cstr) = char(channelfileObject.Properties.Source);%update foStruct.(cstr)
        elseif ~isempty(a) && isempty(imgstack)  %if the scene is same but unloaded
            cd(pathStruct.mstackPath)
            ff = dir(strcat('*',ImageDetails.Scene,'*',cstr,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            imgstack = channelfileObject.flatstack;
            foStruct.(cstr) = char(channelfileObject.Properties.Source);%update foStruct.(cstr)
        else
        end
    else %if no foStruct.(cstr), then
        cd(pathStruct.mstackPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',cstr,'*'));
        filename = char(ff.name);
        channelfileObject = matfile(filename);
        imgstack = channelfileObject.flatstack;
        foStruct.(cstr) = char(channelfileObject.Properties.Source);%update foStruct.(cstr)
    end
    
    stackstruct.(cstr) = imgstack;
end


%% iterate through all segmentation channels
updatebkgmedianmat = true;

segfnames = fieldnames(segStruct);
for ic = 1:length(segfnames)
    cstr = segfnames{ic};
    imgstack = stackstruct.(cstr);
    %choose the channel image
    
    aa = regexp(cstr,'_seg');
    segstr = cstr(1:aa-1);

    if ~isempty(foStruct.(cstr)) %if channelfileObject has been made, check to see if the scene has changed.
        [a,~] = regexp(foStruct.(cstr),ImageDetails.Scene);
        if isempty(a) %if the scene has changed load the new channelimgstack
            cd(pathStruct.segmentPath)
            ff = dir(strcat('*',ImageDetails.Scene,'*',[segStruct.(cstr) '_' segstr],'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            imgstack = channelfileObject.IfFinal;
            foStruct.(cstr) = char(channelfileObject.Properties.Source);%update foStruct.(cstr)
        elseif ~isempty(a) && isempty(imgstack)  %if the scene is same but unloaded
            cd(pathStruct.segmentPath)
            ff = dir(strcat('*',ImageDetails.Scene,'*',[segStruct.(cstr) '_' segstr],'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            imgstack = channelfileObject.IfFinal;
            foStruct.(cstr) = char(channelfileObject.Properties.Source);%update foStruct.(cstr)
        else
            updatebkgmedianmat = false;
        end
    else %if no foStruct.(cstr), then
        cd(pathStruct.segmentPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',[segStruct.(cstr) '_' segstr],'*'));
        filename = char(ff.name);
        channelfileObject = matfile(filename);
        imgstack = channelfileObject.IfFinal;
        foStruct.(cstr) = char(channelfileObject.Properties.Source);%update foStruct.(cstr)
    end
    
    stackstruct.(cstr) = imgstack;
end




if updatebkgmedianmat %udpdate background subtraction if new files loaded
    bkgmedianmat = struct();
    for ic = 1:length(channelstoinput)
        cstr = channelstoinput{ic};
        medianmat = zeros(size(imgstack,3),1);
        backgroundimgstack = stackstruct.background_seg;
        segmentimgstack = stackstruct.nucleus_seg;
        chanimgstack = stackstruct.(cstr);
        if ~strcmpi('dic',cstr)
            for k = 1:size(imgstack,3)
                
                backgroundimg = backgroundimgstack(:,:,k);
                segimg = segmentimgstack(:,:,k);
                chanimg = chanimgstack(:,:,k);
                
                backgroundimg(segimg) = true;
                bkgpixels = chanimg(~backgroundimg);
                chanbkg = nanmedian(bkgpixels);
                
                
                medianmat(k) = chanbkg;
                chanimgstack(:,:,k) = chanimg - single(chanbkg);
            end
        end
        stackstruct.(cstr) = chanimgstack;
        bkgmedianmat.(cstr) = medianmat;
    end
    
      clear segimg chanimg backgroundimg chanimgstack segmentimgstack backgroundimgstack
end









%choose the segmentation image
trackexists = true;
if isempty(trackmatrix) %Tracked is empty try to load tracking or make new tracking structure
    cd(pathStruct.segmentPath)
    imgfile = dir(strcat('*',ImageDetails.Scene,'*',segStruct.nucleus_seg,'*'));
    
    Tracked = loadTrackedStructure;
    togStruct.trackUpdated = true;
    if isempty(Tracked.trackmatrix)
%         segmentimgstack = stackstruct.nucleus_seg(:,:,t);
        If = stackstruct.nucleus_seg(:,:,t);
        trackexists = false;
    else
    end
else
    %if there exists segmenttracking already...then use that.
end

if togStruct.trackUpdated && trackexists
    Ifstack = IfPerimFunction(Tracked);
    If = Ifstack(:,:,t);
elseif ~togStruct.trackUpdated && trackexists
    If = Ifstack(:,:,t);
else
%     If = bwperim(segmentimgstack(:,:,t));
    If = bwperim(stackstruct.nucleus_seg(:,:,t));   
    
end


fnames = fieldnames(stackstruct);
cstr = ImageDetails.Channel;

if sum(strcmpi(fnames,cstr))>0
    imgstack = stackstruct.(cstr);
    channelimg = imgstack(:,:,t);
    bkgmedian = 0;
elseif strcmp(ImageDetails.Channel,'BKGbinary')
    ncstr = segStruct.background_seg;
    imgstack = stackstruct.(ncstr);
    channelimg = imgstack(:,:,t);

    segmentimgstack = stackstruct.nucleus_seg;
    backgroundimgstack = stackstruct.background_seg;
    segmentimg = segmentimgstack(:,:,t);
    backgroundimg = backgroundimgstack(:,:,t);
    backgroundimg(segmentimg) = true;
    prim = imdilate(bwperim(~logical(backgroundimg)),strel('square',1));
    channelimg(prim) = max(max(channelimg));
    channelimg = channelimg+1000;
    bkgmedian = 0;
    
elseif strcmp(ImageDetails.Channel,'overlay')
    channelimg = zeros(size(imgstack,1),size(imgstack,2),3);
    
    channelimg(:,:,1) = stackstruct.(segStruct.nucleus_seg)(:,:,t);
    channelimg(:,:,2) = stackstruct.(segStruct.background_seg)(:,:,t);
    channelimg(:,:,3) = stackstruct.DIC(:,:,t);
    bkgmedian = 0;
elseif strcmpi(ImageDetails.Channel,'FluorOnlyOverlay')
    channelimg = zeros(size(imgstack,1),size(imgstack,2),3);
    channelimg(:,:,1) = stackstruct.(segStruct.nucleus_seg)(:,:,t);
    channelimg(:,:,2) = stackstruct.(segStruct.background_seg)(:,:,t);
    bkgmedian = 0;
else
    channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
    channelimg(:,:,1) = stackstruct.(segStruct.nucleus_seg)(:,:,t);
    channelimg(:,:,2) = stackstruct.(segStruct.background_seg)(:,:,t);
    bkgmedian = 0;
end


if togStruct.runIterate ==0
    displayImageFunct(If,channelimg,bkgmedian);
elseif togStruct.SmallMovie == true
    displaySmallMovieFunct
else
    
end


end
function contrast_Callback(~,~)
global tcontrast lcontrast togStruct
prompt = {'High Contrast Percent Limit','Low Contrast Percent Limit'};
dlg_title = 'Contrast limits from 0 to 100';
num_lines = 1;
defaultans = {num2str(tcontrast),num2str(lcontrast)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
tcontrast = str2double(answer{1});
lcontrast = str2double(answer{2});

togStruct.updateContrast =1;
setSceneAndTime
togStruct.updateContrast =0;
end
function displayTrackingButton_Callback(~,~)
global togStruct Tracked


if isempty(Tracked.trackmatrix)
    disp('need to run tracking!!!')
else
    if togStruct.displayTrackingToggle == 0
        togStruct.displayTrackingToggle = 1;
    else
        togStruct.displayTrackingToggle =0;
    end
    setSceneAndTime
end


end
function Ifstack = IfPerimFunction(Tracked)
global ImageDetails
imgsize = ImageDetails.ImgSize;
trackmatrix = Tracked.trackmatrix;
ArrayStruct = Tracked.arrayStruct;

segment_Pixels_array = ArrayStruct.pixels;

trackbrackets = 1:10:size(trackmatrix,1);
Ifstack = zeros(imgsize(1),imgsize(2),size(trackmatrix,1));
for t = 1:size(trackmatrix,1)
    imagio = zeros(imgsize(1),imgsize(2));
    PX = segment_Pixels_array{t};
%     imdim = [imgsize(1) imgsize(2)];
%     PX = nuccytopxspitter(PX,imdim);
    alltrackframe = trackmatrix(t,:);
    cellindices = find(~isnan(alltrackframe));
    
    if ~isempty(PX)
        whiteiter = 1:length(PX);
        if ~isempty(cellindices) %if cells are tracked on the frame
            cellnums = alltrackframe(cellindices);
            celltracks = trackmatrix(:,cellindices);
            tracklengths = sum(~isnan(celltracks),1);
            cellnumlog = false(size(PX));
            cellnumlog(cellnums) = true;
            coloriter = cellnums;
            whiteiter = find(~cellnumlog);
            for j = 1:length(coloriter)
                cellnum = coloriter(j);
                px = PX{cellnum};
                imagio(px) = tracklengths(j);
            end
        end
        
        for j = 1:length(whiteiter)
            cellnum = whiteiter(j);
            px = PX{cellnum};
            imagio(px) = 1;
        end
    end
    
    
    %erode, compute logical, then dilate, then remove logical
    se = strel('disk',1);
    imagiolog = imagio>0;
    %erode the removal logical
    imagiolog = imerode(imagiolog,se);
    %dilate imagio
    imagio = imdilate(imagio,se);
    %remove logical.
    imagio(imagiolog) = 0;
    Ifstack(:,:,t) = imagio;
    
end

%Display the annulus for cytoplasm
% Ifstacknuc = Ifstack;
% 
% segment_PixelsCyto_array = ArrayStruct.pixelscyto;
% trackbrackets = 1:10:size(trackmatrix,1);
% Ifstack = zeros(imgsize(1),imgsize(2),size(trackmatrix,1));
% for t = 1:size(trackmatrix,1)
%     imagio = zeros(imgsize(1),imgsize(2));
%     PX = segment_PixelsCyto_array{t};
%     alltrackframe = trackmatrix(t,:);
%     cellindices = find(~isnan(alltrackframe));
%     
%     if ~isempty(PX)
%         whiteiter = 1:length(PX);
%         if ~isempty(cellindices) %if cells are tracked on the frame
%             cellnums = alltrackframe(cellindices);
%             celltracks = trackmatrix(:,cellindices);
%             tracklengths = sum(~isnan(celltracks),1);
%             cellnumlog = false(size(PX));
%             cellnumlog(cellnums) = true;
%             coloriter = cellnums;
%             whiteiter = find(~cellnumlog);
%             for j = 1:length(coloriter)
%                 cellnum = coloriter(j);
%                 px = PX{cellnum};
%                 imagio(px) = tracklengths(j);
%             end
%         end
%         
%         for j = 1:length(whiteiter)
%             cellnum = whiteiter(j);
%             px = PX{cellnum};
%             imagio(px) = 1;
%         end
%     end
%     
%     
% %     %erode, compute logical, then dilate, then remove logical
% %     se = strel('disk',1);
% %     imagiolog = imagio>0;
% %     %erode the removal logical
% %     imagiolog = imerode(imagiolog,se);
% %     %dilate imagio
% %     imagio = imdilate(imagio,se);
% %     %remove logical.
% %     imagio(imagiolog) = 0;
% %     Ifstack(:,:,t) = imagio;
%     Ifstack(:,:,t) = imagio;
%     
% nucim = Ifstacknuc(:,:,t);
% cytoim = Ifstack(:,:,t);
% nucim(cytoim>2) = 200;
% 
% Ifstack(:,:,t) = nucim;
%     
% end



end
function traject = trackingTrajectories(timeFrames)
global Tracked

trackmatrix = Tracked.trackmatrix;
ArrayStruct = Tracked.arrayStruct;
segment_Pixels_array = ArrayStruct.pixels;
segment_Centroid_array = ArrayStruct.centroid;

%now update all the fields
centroids = segment_Centroid_array{1};
if ~isempty(centroids)
    centnew = nan(size(trackmatrix,2),size(centroids,2),timeFrames);
    for i = 1:timeFrames
        a=i;
        trackvals=trackmatrix(i,:);
        trackidx = ~isnan(trackvals);
        tracknums = trackvals(trackidx);
        px = segment_Pixels_array{a};
        centroids = segment_Centroid_array{a};
        if ~isempty(px)
            centnew(trackidx,:,i) = centroids(tracknums,:);
        end
    end
    traject = centnew;
end
end
function displayImageFunct(If,channelimg,bkgmedian)
global mainfig oldscenename togStruct timeFrames tcontrast lcontrast MainAxes ImageDetails prcntl lprcntl cmap cmaplz expDetailsStruct

scenestr = ImageDetails.Scene;
tnum = ImageDetails.Frame;
figure(mainfig)

%delete old images to keep memory usage low
set(MainAxes,'NextPlot','replacechildren')

% children1 = findobj(MainAxes,'Type','image');
% children2 = findobj(MainAxes,'Type','Line');
% children3 = findobj(MainAxes,'Type','Text');
% delete(children1);
% delete(children2);
% delete(children3);
% set(MainAxes,'NextPlot','add')

%constrain image axis to be square initially
set(MainAxes,'Units','pixels');
pos = get(MainAxes,'Position');
pos(3:4) = [min([pos(3) pos(4)]) min([pos(3) pos(4)])];
set(MainAxes,'Position',pos);
set(MainAxes,'Units','normalized');
pos = get(MainAxes,'Position');
pos(2)= 0.5 - pos(4)/2;
set(MainAxes,'Position',pos);



% important for updating the contrast
%update contrast if time =1
togStruct.changeSceneOrChannel=0;
if tnum==1
    oldscenename='new';
    togStruct.changeSceneOrChannel=1;
end

%update contrast if channel has changed
if ~strcmp(ImageDetails.Channel,oldscenename)
    togStruct.changeSceneOrChannel=1;
    oldscenename = ImageDetails.Channel;
end


%scripts for displaying contrasted image
alog = regexpi(ImageDetails.Channel,'overlay'); %when overlay display is desired
if ~isempty(alog)
    disprgb = zeros(size(channelimg));
    channelimgrgb = channelimg;

    Ifrgb = ind2rgb(uint8(If),vertcat([0 0 0],cmap(257:end,:)));
    
    for i = 1:size(channelimg,3)
        if i==3
            channelimg = channelimgrgb(:,:,i);
            lprcntl = prctile(channelimg(:),lcontrast);
            prcntl = prctile(channelimg(:),tcontrast);
            scaleFactor = 1./(prcntl - lprcntl);
            dispimg = channelimg.*scaleFactor;
            dispimg = dispimg-(lprcntl.*scaleFactor);
            dispimg(dispimg> 1) =1;
            dispimg(dispimg<0) = 0;
        else
            channelimg = channelimgrgb(:,:,i);
            lprcntl = 0;
            prcntl = prctile(channelimg(:),tcontrast);
            scaleFactor = 1./(prcntl - lprcntl);
            dispimg = channelimg.*scaleFactor;
            dispimg = dispimg-(lprcntl.*scaleFactor);
            dispimg(dispimg> 1) =1;
            dispimg(dispimg<0) = 0;
        end
        
        
        if i==1
            disprgb(:,:,i) = dispimg;
        elseif i==2
            disprgb(:,:,i) = dispimg;
        elseif i==3
            if sum(dispimg(:))>100
                for j = 1:3
                    cimg =  disprgb(:,:,j);
                    cimg = dispimg.*cimg;
                    disprgb(:,:,j) = cimg;
                end
                
                for j = 1:3
                    cimg =  disprgb(:,:,j);
                    cimg = cimg + dispimg./3;
                    cimg(cimg>1) = 1;
                    disprgb(:,:,j) = cimg;
                end
            end
        end
    end
    
    newimg = zeros(size(disprgb));
    for i = 1:size(disprgb,3)
        Ifsub = Ifrgb(:,:,i);
        Ifnew = zeros(size(Ifsub));
        Ifnew(If>0) = Ifsub(If>0);
        dispimg = disprgb(:,:,i);
        dispimg(If>0) = Ifsub(If>0);
        disprgb(:,:,i) = dispimg;
        newimg(:,:,i) = Ifnew;
    end
    
    
    dispimg = disprgb;
    image(MainAxes,dispimg);
    
else  %under normal circumstances
    if togStruct.changeSceneOrChannel ==1
        prcntl = prctile(channelimg(:),tcontrast);
        togStruct.changeSceneOrChannel =0;
    end
    lprcntl = bkgmedian;
    scaleFactor = 255./(prcntl - lprcntl);
    dispimg = channelimg.*scaleFactor;
    dispimg = dispimg-(lprcntl.*scaleFactor);
    dispimg(dispimg > 250) =253;
    dispimg(dispimg < 0) = 0;
    %     dispimg(If>0)=255;
    dispimg(If>0)=If(If>0)+257;
    imagesc(MainAxes,int16(dispimg));
%     colormap(MainAxes,cmap);
%     MainAxes.ColorOrder
    set(MainAxes,'CLim',[0 size(cmap,1)]);
end

%title the displayed image
ttl = get(MainAxes,'Title');
set(ttl,'String',{ImageDetails.Channel,[expDetailsStruct.expDateTitleStr] ,scenestr});
set(ttl,'FontSize',12); 

frametext= text(MainAxes,0,0,[' frame ' num2str(tnum) '/' num2str(timeFrames)]);
% frametext.Units = 'normalized';
imgsize = ImageDetails.ImgSize;
textspace = round(imgsize(1)/100);
texty = imgsize(2)-textspace*2;
frametext.Position = [textspace texty];
frametext.FontSize = 12;
frametext.Color =[1 0 0];

%plot cell tracking "tails"
if ~(tnum==1)
    if togStruct.displayTrackingToggle==1
        trajectForPlot = trackingTrajectories(timeFrames);
        mainplotX = squeeze(trajectForPlot(:,1,:)); %28x22 means 28 cells on frame 22;
        mainplotY = squeeze(trajectForPlot(:,2,:));
        
        if size(mainplotY,2) == 1
            mainplotY=mainplotY';
            mainplotX=mainplotX';
        end
        
        %only plot if the cell is currently tracked/segmented in this frame
        idx = ~isnan(mainplotY(:,tnum));
        idxa = find(idx==1);
        if ~isempty(idxa)
            
            plotFCcolor=0;
            if plotFCcolor==0
                if togStruct.trackUpdated
                    cmaplz = colorcubemodified(length(idx),'jet');
                elseif ~(size(cmaplz,1)==length(idx))
                    cmaplz = colorcubemodified(length(idx),'jet');
                end
                cmapl = cmaplz;
                plotcmap = zeros(length(idxa),3);
                for i=1:length(idxa)
                    plotcmap(i,:) =  cmapl(idxa(i),:);
                end
                
                
                
                x = mainplotX(idx,:);
                y = mainplotY(idx,:);
                t0 =max([1 tnum-20]);
                h = line(MainAxes,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',2);
                set(h, {'color'}, num2cell(plotcmap,2));
                % only plot if the cell is currently tracked/segmented in this frame
                
            else
                if togStruct.trackUpdated
                    cmaplz = colorcubemodified(length(idx),'colorcube');
                elseif ~(size(cmaplz,1)==length(idx))
                    cmaplz = colorcubemodified(length(idx),'colorcube');
                end
                %                 cmaplz = feval('jet',length(idx));
                cmapl = cmaplz;
                
                [fc,abs,pidx] = Plot_callback([],[]);
                [ss,sidx] = sort(fc(:,tnum+2));
                
                
                

                plotcmap = zeros(length(idxa),3);
                for i=1:length(idxa)
                    plotcmap(i,:) =  cmapl(idxa(i),:);
                end
                
                
                x = mainplotX(idx,:);
                y = mainplotY(idx,:);
                fcy = fc(idx,tnum+20);
                fcy = fc(idx,tnum+20)-fc(idx,tnum);
                [~,sidx] = sort(fcy);
                plotcmap = flipud(feval('viridis',length(sidx)));
                %             x = mainplotX(sidx,:);
                %             y = mainplotY(sidx,:);
                t0 =max([1 tnum-20]);
                h = line(MainAxes,x(sidx,t0:tnum)',y(sidx,t0:tnum)','LineWidth',2);
                set(h, {'color'}, num2cell(plotcmap,2));
                togStruct.trackUpdated=1;
            end


        end
    end
end

togStruct.trackUpdated = false;
drawnow
end

function saveMovieCropped(~,~)
global togStruct
togStruct.SaveMovie = true;
makeMovieCroppedYEAH([],[])
togStruct.SaveMovie = false;
end

function makeMovieCropped(~,~)
global togStruct
togStruct.SaveMovie = false;
makeMovieCroppedYEAH([],[])
togStruct.SaveMovie = false;
end

function saveFramesCropped(~,~)
global togStruct
togStruct.SaveMovie = true;
exportFramesmakeMovieCroppedYEAH([],[])
togStruct.SaveMovie = true;
end



function makeMovieAllChannels(~,~)
global trackfilestr expDetailsStruct dirStruct timeSteps timeVec pathStruct psettings medmeantot Tracked segStruct Ifstack channelimgstack nucleusimgstack DICimgstack togStruct timeFrames tcontrast lcontrast  ImageDetails prcntl lprcntl cmap cmaplz



contrastFrame = ImageDetails.Frame;
% screen_capture(movie_name,recording_time)
bkgmedian=0;


% click the desired cell
ArrayStruct = Tracked.arrayStruct;
trackmatrix = Tracked.trackmatrix;
%choose the initial cell




%load traces to be plotted
if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
end
framesThatMustBeTracked = psettings.framesThatMustBeTracked;


smooththat=0;
[plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,pathStruct,ImageDetails,timeFrames,smooththat);


stimulationFrame = find(~(timeVec<0),1,'first');


cropbuffer = 10;
imgdim = ImageDetails.ImgSize;
celltrack = ~isnan(trackmatrix);
nanvec = sum((celltrack),2);
% nannum1 = max([2 find(nanvec == max(nanvec),1,'first')]);
nannum1 = max([2 find(nanvec > 0,1,'first')]);
% nannum1 = 1;
nannum2 = find(nanvec == max(nanvec),1,'last');
segment_Centroids_array = ArrayStruct.centroid;

%determine crop region
xvals = nan(nannum2,2);
yvals = nan(nannum2,2);
cycle = 0;
for tnum = nannum1:nannum2
    cycle = cycle+1;
    cellindices = trackmatrix(tnum,:);
    centroids = segment_Centroids_array{tnum};
    cellcentroids = centroids(cellindices(~isnan(cellindices)),:);
    yvals(cycle,:) = round([min(cellcentroids(:,1)) max(cellcentroids(:,1))]);
    xvals(cycle,:) = round([min(cellcentroids(:,2)) max(cellcentroids(:,2))]);
end
cropxmin = max([min(xvals(:,1))-cropbuffer 1]);
cropxmax = min([max(xvals(:,2))+cropbuffer imgdim(1)]);
cropymin = max([min(yvals(:,1))-cropbuffer 1]);
cropymax = min([max(yvals(:,2))+cropbuffer imgdim(2)]);
xdist = cropxmax-cropxmin;
ydist = cropymax-cropymin;

ddiff = abs(ydist-xdist);
dd2x = floor(ddiff/2);
dd2y = ceil(ddiff/2);
if ddiff>0
    if ydist>xdist
        x1 = cropxmin-1;
        x2 = 512-cropxmax;
        if x1 < dd2x
            dd2x = x1;
            dd2y = ddiff-x1;
        elseif x2 < dd2y
            dd2y = x2;
            dd2x = ddiff-x2;
        end
        cropxmin = cropxmin-dd2x;
        cropxmax = cropxmax+dd2y;
    else
        x1 = cropymin-1;
        x2 = 512-cropymax;
        if x1 < dd2x
            dd2x = x1;
            dd2y = ddiff-x1;
        elseif x2 < dd2y
            dd2y = x2;
            dd2x = ddiff-x2;
        end
        cropymin = cropymin-dd2x;
        cropymax = cropymax+dd2y;
    end
end

xdist = cropxmax-cropxmin;
ydist = cropymax-cropymin;



% set up figure
f33 = figure(33);
close 33

%%%% make movie by iterating through frames
f33 = figure(33);
f33.Color = 'w';
f33.Units = 'pixels';
f33.Position = [100 100 (xdist*6) ydist];

ax1 = axes();
ax1.Units = 'pixels';
ax1.Position = [0 0 xdist ydist];

ax2 = axes();
ax2.Units = 'pixels';
ax2.Position = [xdist 0 xdist ydist];

ax3 = axes();
ax3.Units = 'pixels';
ax3.Position = [(xdist*2) 0 xdist ydist];

ax4 = axes();
ax4.Units = 'pixels';
ax4.Position = [(xdist*3) 0 xdist ydist];

ax5 = axes();
ax5.Units = 'pixels';
ax5.Position = [(xdist*4) 0 xdist ydist];

ax6 = axes();
ax6.Units = 'pixels';
ax6.Position = [(xdist*5) 0 xdist ydist];





%set contrast
t = contrastFrame;
contrastFrameSmad = find(timeVec>30,1,'first');
contrastFrameSnail = find(timeVec>300,1,'first');
If = Ifstack(:,:,t);
nucleusImg = nucleusimgstack(:,:,contrastFrameSnail);
cellImg = channelimgstack(:,:,contrastFrameSmad);
DICImg = DICimgstack(:,:,t);

ImageDetails.Channel = 'overlay';
channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
channelimg(:,:,1) = nucleusImg;
channelimg(:,:,2) = cellImg;
channelimg(:,:,3) = DICImg;


lc = nan(1,3);
tc = nan(1,3);
%scripts for displaying contrasted image
channelimgrgb = channelimg;
for i = 1:size(channelimgrgb,3)
    if i==3
        channelimg = channelimgrgb(:,:,i);
        lc(i) = prctile(channelimg(:),lcontrast);
        tc(i) = prctile(channelimg(:),tcontrast);
    else
        channelimg = channelimgrgb(:,:,i);
        lc(i) = 0;
        tc(i) = prctile(channelimg(:),tcontrast);
    end
end


%display images
for tnum = nannum1:nannum2
    t = tnum;
    
    If = Ifstack(:,:,t);
    nucleusImg = nucleusimgstack(:,:,t);
    cellImg = channelimgstack(:,:,t);
    DICImg = DICimgstack(:,:,t);
    
    ImageDetails.Channel = 'overlay';
        channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
        channelimg(:,:,1) = nucleusImg;
        channelimg(:,:,2) = cellImg;
        channelimg(:,:,3) = DICImg;
    
    if length(size(channelimg))>2
        cropim = channelimg(cropxmin:cropxmax,cropymin:cropymax,:);
    else
        cropim = channelimg(cropxmin:cropxmax,cropymin:cropymax);
    end
    If = If(cropxmin:cropxmax,cropymin:cropymax);
    channelimg = cropim;
    
    If(If<2)=0;
    If(If>2) = 40;
%     If(If>1) = 0;
     
    
    %scripts for displaying contrasted image
    disprgb = zeros(size(cropim));
    channelimgrgb = cropim;
    Ifrgb = ind2rgb(uint8(If),vertcat([0 0 0],cmap(257:end,:)));
    for i = 1:size(cropim,3)
        if i==3
            channelimg = channelimgrgb(:,:,i);
            lc(i) = prctile(channelimg(:),lcontrast);
            tc(i) = prctile(channelimg(:),tcontrast);
            scaleFactor = 1./(tc(i) - lc(i));
            dispimg = channelimg.*scaleFactor;
            dispimg = dispimg-(lc(i).*scaleFactor);
            dispimg(dispimg> 1) =1;
            dispimg(dispimg<0) = 0;
        else
            channelimg = channelimgrgb(:,:,i);
            scaleFactor = 1./(tc(i) - lc(i));
            dispimg = channelimg.*scaleFactor;
            dispimg = dispimg-(lc(i).*scaleFactor);
            dispimg(dispimg> 1) =1;
            dispimg(dispimg<0) = 0;
        end
        
        
        if i==1
            disprgb(:,:,i) = dispimg;
        elseif i==2
            disprgb(:,:,i) = dispimg;
        elseif i==3
            OGrgb = disprgb;

            DICrgb = channelimgrgb;
            for j = 1:3
                cimg =  disprgb(:,:,j);
                cimgpre = cimg;
%                 cimg(dispimg>(cimg./a)) = dispimg(dispimg>(cimg./a));
                cimg = cimg.*dispimg;
                disprgb(:,:,j) = cimg;
                DICrgb(:,:,j) = dispimg;
            end
            cimgmax =  max(disprgb,[],3);
            a = 0.1;
            for j = 1:3
                cimg =  disprgb(:,:,j);
                cimg = cimg + dispimg./10;
                cimg(cimg>1) = 1;
                disprgb(:,:,j) = cimg;
            end
            
        end
    end
    
    
    dispimg = disprgb;
    im1 = OGrgb;
    im2 = OGrgb;
    im3 = OGrgb;
    im4 = OGrgb;
    
   
    im1(:,:,2:3) = 0;
    im2(:,:,[1 3]) = 0;
    im3 = DICrgb;
    im4 = dispimg;
    
    disprgb = im2;
    %add nuclei segmentation overlay???
    newimg = zeros(size(disprgb));
    for i = 1:size(disprgb,3)
        Ifsub = Ifrgb(:,:,i);
        Ifnew = zeros(size(Ifsub));
        Ifnew(If>0) = Ifsub(If>0);
        dispimg = disprgb(:,:,i);
        dispimg(If>0) = Ifsub(If>0);
        disprgb(:,:,i) = dispimg;
        newimg(:,:,i) = Ifnew;
    end
    im2 = disprgb;

    disprgb = im1;
    %add nuclei segmentation overlay???
    newimg = zeros(size(disprgb));
    for i = 1:size(disprgb,3)
        Ifsub = Ifrgb(:,:,i);
        Ifnew = zeros(size(Ifsub));
        Ifnew(If>0) = Ifsub(If>0);
        dispimg = disprgb(:,:,i);
%         dispimg(If>0) = Ifsub(If>0);
        disprgb(:,:,i) = dispimg;
        newimg(:,:,i) = Ifnew;
    end
    im1 = disprgb;
    
    
    
    image(ax4,im2);
    image(ax3,im1);
    image(ax1,im3);
    image(ax2,im4);


    minstr = '0000';
    minsub = num2str(round(timeVec(tnum)./60,1,'decimals'));
    minidx = false(size(minstr));
    minidx(end-(length(minsub)-1):end) = true;
    if timeVec(tnum)<0
        minstr(end-(length(minsub)-1):end) = minsub;
        minstr(~minidx) = ' ';
    else
        minstr(end-(length(minsub)-1):end) = minsub;

        minstr(~minidx)=' ';
    end
    
    frametext= text(ax4,0,0,[minstr ' hr']);
    frametext.Units = 'pixels';
    frametext.Position = [xdist-5 15];
    frametext.FontSize = 15;
    frametext.Color =[1 1 0];
    frametext.FontWeight = 'bold';
    frametext.HorizontalAlignment = 'right';
    
    if ~(timeVec(tnum)<0)
        frametext= text(ax4,0,0,'+TGFbeta');
        frametext.Units = 'pixels';
        frametext.Position = [xdist-5 35];
        frametext.FontSize = 15;
        frametext.Color =[1 1 0];
        frametext.FontWeight = 'bold';
        frametext.HorizontalAlignment = 'right';
    end
    
    
    
    f=gcf;
    for h = f.Children'
        if strcmp(h.Type,'axes')
            h.XTick = [];
            h.YTick = [];
            set(h,'NextPlot','Add')
        end
    end
    
    

    %plot cell tracking "tails"
    if tnum>0
        
        trajectForPlot = trackingTrajectories(timeFrames);
        mainplotX = squeeze(trajectForPlot(:,1,:)); %28x22 means 28 cells on frame 22;
        mainplotY = squeeze(trajectForPlot(:,2,:));
        
        if size(mainplotY,2) == 1
            mainplotY=mainplotY';
            mainplotX=mainplotX';
        end
        %only plot if the cell is currently tracked/segmented in this frame
        idx = ~isnan(mainplotY(:,tnum));
        idxa = find(idx);
        
        cmapl = cool(size(mainplotY,1));
        cmapl = feval('lines',(size(mainplotY,1)));
        plotcmap = zeros(length(idxa),3);
        for i=1:length(idxa)
            plotcmap(i,:) =  cmapl(idxa(i),:);
        end
        %
        %             x = mainplotX(idxa,:) - cellcentroid(1) + cropsize;
        %             y = mainplotY(idxa,:) - cellcentroid(2) + cropsize;
        x = mainplotX(idxa,:) - cropymin;
        y = mainplotY(idxa,:) - cropxmin;
        t0 =max([max([1 nannum1-1]) tnum-20]);
        %             h = plot(MainAxes,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',2);
        if togStruct.displayTrackingToggle==1 && tnum>1
            h = line(ax4,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',3);
            set(h, {'color'}, num2cell(plotcmap,2));
            h = line(ax3,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',3);
            set(h, {'color'}, num2cell(plotcmap,2));
        end
    end
    
    
    
    
    
    
    
        
    %plot traces
    mmsstr = 'median';
    fcplotstr = 'FC';
    ylabelstr = 'fold change';
    
%     mmsstr = 'median';
%     fcplotstr = '';
%     ylabelstr = 'abundance';
    
%     mmsstr = 'nuccyto';
%     fcplotstr = '';
%     ylabelstr = 'nuccyto';
    %first axes, plot abundance
    

    plotstr = ['Smad' fcplotstr];
    plotMat = plotStructUI.(mmsstr).(plotstr);
    plotMatStat = plotMat(:,nannum1:nannum2);
    ylimit = [prctile(plotMatStat(:),0)*0.95 prctile(plotMatStat(:),100)*1.2];
    stem(ax5,0,10,'LineStyle',':','Color',[0.7 0.7 0.7],'LineWidth',3); hold on
    
    h = plot(ax5,timeVec(nannum1:tnum),plotMat(idx,nannum1:tnum)','LineWidth',2.5);hold off
    if togStruct.displayTrackingToggle ==1
        for i=1:length(h)
            h(i).Color = cmapl(i,:);
        end
        colormap(cmap);
    end
    
    ax5.XLim = timeVec([nannum1 nannum2]);
    ax5.YLim = (ylimit);
    ax5.YLabel.String = ['Smad3 ' ylabelstr ' ' mmsstr];
    ax5.XLabel.String = 'minutes';
    ax5.XGrid = 'on';
    ax5.YGrid = 'on';
    ax5.Color = [0.95 0.95 0.95];
    ax5.GridColor = [1 1 1];
    ax5.GridAlpha = 1;
    ax5.Title.String = ['mNG-Smad3'];
    ax5.FontSize = 10;
    
    xtick = [0:60:timeVec(nannum2)+240]-240;
    ax5.XTick = xtick;
    xticklabel = mat2cell(xtick,1,ones(1,length(xtick)));
    xtidx = false(size(xtick));
    xtidx(1:floor(length(xtick)/5):end) = true;
%     xtidx(xtick==0) = true;
    xticklabel(~xtidx) = {''};
    
    ax5.XTickLabel = xticklabel;
    ytick = [0:0.5:20]-5;
    ax5.YTick = ytick;
    yticklabel = mat2cell(ytick,1,ones(1,length(ytick)));
    yticklabel(2:2:end) = {''};
    ax5.YTickLabel = yticklabel;
    ax5.Box = 'off';
    
    
  
    xshift = 40;
    yshift = 50;
    ax5.Position = [(xdist*4)+xshift yshift xdist-xshift-15 ydist-yshift-25];
    
    
    
    
        %plot traces
    mmsstr = 'mean';
    %first axes, plot abundance
    fcplotstr = '';
    ylabelstr = 'abundance';
    plotstr = ['mkate' fcplotstr];
    plotMat = plotStructUI.(mmsstr).(plotstr);
    pm = plotMat(:,stimulationFrame);
    plotMat = plotMat./nanmedian(pm(:));
    plotMatStat = plotMat(:,nannum1:nannum2);
    ylimit = [prctile(plotMatStat(:),0)*0.95 prctile(plotMatStat(:),100)*1.2];
    stem(ax6,0,10,'LineStyle',':','Color',[0.7 0.7 0.7],'LineWidth',3); hold on
    h = plot(ax6,timeVec(nannum1:tnum),plotMat(idx,nannum1:tnum)','LineWidth',2.5);hold off
    if togStruct.displayTrackingToggle ==1
        for i=1:length(h)
            h(i).Color = cmapl(i,:);
        end
        colormap(cmap);
    end
    
    ax6.XLim = timeVec([nannum1 nannum2]);
    ax6.YLim = (ylimit);
    ax6.YLabel.String = ['mCherry ' ylabelstr ' ' mmsstr];
    ax6.XLabel.String = 'minutes';
    ax6.XGrid = 'on';
    ax6.YGrid = 'on';
    ax6.Color = [0.95 0.95 0.95];
    ax6.GridColor = [1 1 1];
    ax6.GridAlpha = 1;
    ax6.Title.String = ['Snail-mCherry'];
    ax6.FontSize = 10;
    
    xtick = [0:60:timeVec(nannum2)+240]-240;
    ax6.XTick = xtick;
    xticklabel = mat2cell(xtick,1,ones(1,length(xtick)));
    xtidx = false(size(xtick));
    xtidx(1:floor(length(xtick)/5):end) = true;
%     xtidx(xtick==0) = true;
    xticklabel(~xtidx) = {''};
    ax6.XTickLabel = xticklabel;
    
    ytick = [0:0.25:20]-5;
    ax6.YTick = ytick;
    yticklabel = mat2cell(ytick,1,ones(1,length(ytick)));
    yticklabel(2:2:end) = {''};
    ax6.YTickLabel = yticklabel;
    ax6.Box = 'off';
    
    
    
  

    ax6.Position = [(xdist*5)+xshift yshift xdist-xshift-15 ydist-yshift-25];
    
        f=gcf;
    for h = f.Children'
        if strcmp(h.Type,'axes')
            set(h,'NextPlot','Replace')
        end
    end
    
    
    
    
    
    
    if tnum==1
        tstepdelay = 0; 
    else
        tstepdelay = timeSteps(tnum-1);
    end
%     pause(tstepdelay*0.015)
    pause(tstepdelay*0.005)
    drawnow
    
    
%     channelimg = horzcat(im3,im2,im1,im4);
    
%     togStruct.SaveMovie = false;
    if togStruct.SaveMovie
        specialdir = [dirStruct.parentdir 'LookingAtData/AnimatedTrackingPlottingMovie'];
        olddir = pwd;
        saveimg = channelimg;
        if ~isdir(specialdir)
            mkdir(specialdir)
        end
        cd(specialdir)
        savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' tif frames'];
        if ~isdir(savedir)
            mkdir(savedir)
        end
        cd(savedir)
        savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.tif'];
        imwrite(uint16(saveimg*30000),savename,'tif');
        
        cd(specialdir)
        savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' figure tiff frames'];
        if ~isdir(savedir)
            mkdir(savedir)
        end
        cd(savedir)
        savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.tif'];
        saveas(f33,savename,'tif')
        disp(num2str(ImageDetails.Frame));
        cd(olddir)
    end
end

end
function exportFramesmakeMovieCroppedYEAH(~,~)
%FOR IMMUNO!
global stackstruct channelstoinput ExportName ExportNameKey trackfilestr expDetailsStruct dirStruct timeVec pathStruct psettings medmeantot Tracked segStruct Ifstack channelimgstack nucleusimgstack DICimgstack togStruct timeFrames tcontrast lcontrast  ImageDetails prcntl lprcntl cmap cmaplz

bkgmedian=0;
%LONG TRACK



%choose the initial cell
contrastVal = 4000;
% scenesTo = {'s005','s015','s026','s035','s043','s054','s064','s072','s085','s098','s102','s114','s126','s138'};
scenesTo = {'s002','s004','s011','s019','s022','s026','s027','s182','s192'};
scenesTo = {'s03','s64'}
scenesTo = {'s02','s07','s11','s15','s19'};
% scenesTo = {'s03','s08','s25','s26'};
% scenesTo = {'s01','s05'};

trange = [3];




%load traces to be plotted
if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
end
framesThatMustBeTracked = psettings.framesThatMustBeTracked;



for io = 1:length(scenesTo)
    scenestr = scenesTo{io};
    for jimmy = 1:length(channelstoinput)
        ImageDetails.Channel = channelstoinput{jimmy};
        
        ImageDetails.Scene = scenestr;
        setSceneAndTime
        
        
        
        
        if togStruct.plotSettingsToggle == 0
            psettings = PlotSettings_callback([],[]);
            togStruct.plotSettingsToggle=1;
        end
        
        % %remove all global variables before parfor loop
        eNameKey = ExportNameKey;
        eName = ExportName;
        tPath = pathStruct.trackingPath;
        
        
        cd(tPath)
        trackfile = dir(strcat(eNameKey,'*',scenestr,'*',eName,'.mat'));
        trackfilename = char({trackfile.name});
        
        if ~isempty(trackfilename)
            trackedArray = loadTrackedArray(trackfilename); %trackedArray = Tracked
        end
        Tracked = trackedArray;
        % click the desired cell
        
        ArrayStruct = trackedArray.arrayStruct;
        trackmatrix = trackedArray.trackmatrix;
        smooththat=0;
        [plotStructUI] = plotthemfunction(framesThatMustBeTracked,trackedArray,pathStruct,ImageDetails,timeFrames,smooththat);
        
        
        % set up figure
        f33 = figure(33);
        close 33
        
        %%%% make movie by iterating through frames
        f33 = figure(33);
        f33.Color = 'w';
        f33.Units = 'pixels';
        f33.Position = [593 156 1460 760];
        
        ax1 = axes();
        ax1.Units = 'pixels';
        ax1.Position = [25 25 720 720];
        
        ax2 = axes();
        ax2.Units = 'pixels';
        ax2.Position = [1000 425 335 305];
        
        ax3 = axes();
        ax3.Units = 'pixels';
        ax3.Position = [1000 50 335 305];
        
        f33.Units = 'normalized';
        h = findobj(f33,'Type','Axes');
        [h.Units] = deal('normalized');
        f33.Position = [0.05 0.05 0.4 0.5];
        
        
        cropbuffer = 20;
        imgdim = ImageDetails.ImgSize;
        celltrack = ~isnan(trackmatrix);
        nanvec = sum((celltrack),2);
        % nannum1 = max([2 find(nanvec == max(nanvec),1,'first')]);
        nannum1 = max([2 find(nanvec > 0,1,'first')]);
        % nannum1 = 1;
        
        % nannum2 = find(nanvec == max(nanvec),1,'last');
        nannum2 = find(nanvec > 0,1,'last');
        segment_Centroids_array = ArrayStruct.centroid;
        
        %determine crop region
        xvals = nan(nannum2,2);
        yvals = nan(nannum2,2);
        cycle = 0;
        for tnum = nannum1:nannum2
            cycle = cycle+1;
            cellindices = trackmatrix(tnum,:);
            centroids = segment_Centroids_array{tnum};
            cellcentroids = centroids(cellindices(~isnan(cellindices)),:);
            yvals(cycle,:) = round([min(cellcentroids(:,1)) max(cellcentroids(:,1))]);
            xvals(cycle,:) = round([min(cellcentroids(:,2)) max(cellcentroids(:,2))]);
        end
        cropxmin = max([min(xvals(:,1))-cropbuffer 1]);
        cropxmax = min([max(xvals(:,2))+cropbuffer imgdim(1)]);
        cropymin = max([min(yvals(:,1))-cropbuffer 1]);
        cropymax = min([max(yvals(:,2))+cropbuffer imgdim(2)]);
        xdist = cropxmax-cropxmin;
        ydist = cropymax-cropymin;
        
        ddiff = abs(ydist-xdist);
        dd2x = floor(ddiff/2);
        dd2y = ceil(ddiff/2);
        if ddiff>0
            if ydist>xdist
                x1 = cropxmin-1;
                x2 = 512-cropxmax;
                if x1 < dd2x
                    dd2x = x1;
                    dd2y = ddiff-x1;
                elseif x2 < dd2y
                    dd2y = x2;
                    dd2x = ddiff-x2;
                end
                cropxmin = cropxmin-dd2x;
                cropxmax = cropxmax+dd2y;
            else
                x1 = cropymin-1;
                x2 = 512-cropymax;
                if x1 < dd2x
                    dd2x = x1;
                    dd2y = ddiff-x1;
                elseif x2 < dd2y
                    dd2y = x2;
                    dd2x = ddiff-x2;
                end
                cropymin = cropymin-dd2x;
                cropymax = cropymax+dd2y;
            end
        end
        
        xdist = cropxmax-cropxmin;
        ydist = cropymax-cropymin;
        
        
        %     trange = 1:nannum2;
        for tnum = trange
            t = tnum;
            
            If = Ifstack(:,:,t);
            nucleusImg = stackstruct.(segStruct.nucleus_seg)(:,:,t);
            cellImg = stackstruct.(segStruct.background_seg)(:,:,t);
            
            
            if strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
                channelimg = nucleusImg;
            elseif strcmp(ImageDetails.Channel,segStruct.background_seg)
                channelimg = cellImg;
            elseif strcmp(ImageDetails.Channel,'overlay')
                DICImg = DICimgstack(:,:,t);
                channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
                channelimg(:,:,1) = nucleusImg;
                %         channelimg(:,:,1) = zeros(size(nucleusImg));
                channelimg(:,:,2) = cellImg;
                channelimg(:,:,3) = DICImg;
            elseif strcmpi(ImageDetails.Channel,'FluorOnlyOverlay')
                DICImg = stackstruct.('DIC')(:,:,t);
                channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
                channelimg(:,:,1) = nucleusImg;
                channelimg(:,:,2) = cellImg;
            elseif ~isempty(regexp(ImageDetails.Channel,'647'))
                channelimg = stackstruct.('AlexaFluor647')(:,:,t);
            elseif strcmpi(ImageDetails.Channel,'DIC')
                channelimg = stackstruct.('DIC')(:,:,t);
            end
            
            if length(size(channelimg))>2
                cropim = channelimg(cropxmin:cropxmax,cropymin:cropymax,:);
            else
                cropim = channelimg(cropxmin:cropxmax,cropymin:cropymax);
            end
            If = If(cropxmin:cropxmax,cropymin:cropymax);
            channelimg = cropim;
            
            If(If<2)=0;
            If(If>2) = 100;
            %     If(If>1) = 0;
            
            
            
            
            
            
            %scripts for displaying contrasted image
            alog = regexpi(ImageDetails.Channel,'overlay'); %when overlay display is desired
            if ~isempty(alog)
                %scripts for displaying contrasted image
                disprgb = zeros(size(cropim));
                channelimgrgb = cropim;
                Ifrgb = ind2rgb(uint8(If),vertcat([0 0 0],cmap(257:end,:)));
                for i = 1:size(cropim,3)
                    if i==3
                        channelimg = channelimgrgb(:,:,i);
                        lprcntl = prctile(channelimg(:),lcontrast);
                        prcntl = prctile(channelimg(:),tcontrast);
                        scaleFactor = 1./(prcntl - lprcntl);
                        dispimg = channelimg.*scaleFactor;
                        dispimg = dispimg-(lprcntl.*scaleFactor);
                        dispimg(dispimg> 1) =1;
                        dispimg(dispimg<0) = 0;
                    else
                        channelimg = channelimgrgb(:,:,i);
                        lprcntl = 0;
                        prcntl = prctile(channelimg(:),tcontrast);
                        scaleFactor = 1./(prcntl - lprcntl);
                        dispimg = channelimg.*scaleFactor;
                        dispimg = dispimg-(lprcntl.*scaleFactor);
                        dispimg(dispimg> 1) =1;
                        dispimg(dispimg<0) = 0;
                    end
                    
                    
                    if i==1
                        disprgb(:,:,i) = dispimg;
                    elseif i==2
                        disprgb(:,:,i) = dispimg;
                    elseif i==3
                        OGrgb = disprgb;
                        
                        DICrgb = channelimgrgb;
                        for j = 1:3
                            cimg =  disprgb(:,:,j);
                            cimgpre = cimg;
                            %                 cimg(dispimg>(cimg./a)) = dispimg(dispimg>(cimg./a));
                            cimg = cimg.*dispimg;
                            disprgb(:,:,j) = cimg;
                            DICrgb(:,:,j) = dispimg;
                        end
                        cimgmax =  max(disprgb,[],3);
                        a = 0.1;
                        for j = 1:3
                            cimg =  disprgb(:,:,j);
                            cimg = cimg + dispimg./5;
                            cimg(cimg>1) = 1;
                            disprgb(:,:,j) = cimg;
                        end
                        
                    end
                end
                
                newimg = zeros(size(disprgb));
                for i = 1:size(disprgb,3)
                    Ifsub = Ifrgb(:,:,i);
                    Ifnew = zeros(size(Ifsub));
                    Ifnew(If>0) = Ifsub(If>0);
                    %             Ifnew(If>0) = 0;
                    dispimg = disprgb(:,:,i);
                    dispimg(If>0) = Ifsub(If>0);
                    %             dispimg(If>0) = 0;
                    disprgb(:,:,i) = dispimg;
                    newimg(:,:,i) = Ifnew;
                end
                
                
                dispimg = disprgb;
                %     dispimg = uint8(dispimg);
                image(ax1,dispimg);
                %     colormap(MainAxes,cmap);
                %                 set(MainAxes,'CLim',[0 size(cmap,1)]);
                %             set(MainAxes,'YDir',[0 size(cmap,1)]);
                
            else  %under normal circumstances
                if togStruct.changeSceneOrChannel ==1
                    prcntl = prctile(channelimg(:),tcontrast);
                    prcntl = 1800;
                    prcntl = contrastVal;
                    togStruct.changeSceneOrChannel =0;
                end
                %     lprcntl = bkgmedian.*0.90;
                lprcntl = bkgmedian;
                prcntl = contrastVal;
                scaleFactor = 255./(prcntl - lprcntl);
                dispimg = channelimg.*scaleFactor;
                dispimg = dispimg-(lprcntl.*scaleFactor);
                dispimg(dispimg > 250) =253;
                dispimg(dispimg < 0) = 0;
                %     dispimg(If>0)=255;
                %         dispimg(If>0)=If(If>0)+257;
                imagesc(ax1,int16(dispimg));
                colormap(ax1,cmap);
                set(ax1,'CLim',[0 size(cmap,1)]);
                set(ax1,'YDir','normal');
            end
            
            %title the displayed image
            
            
            
            %     frametext= text(ax1,0,0,['frame ' num2str(tnum) '/' num2str(timeFrames)]);
            %     frametext.Units = 'normalized';
            %     frametext.Position = [0.01 0.95];
            %     frametext.FontSize = 12;
            %     frametext.Color =[1 1 0];
            minstr = '0000';
            minsub = num2str(round(timeVec(tnum)./60,1,'decimals'));
            minidx = false(size(minstr));
            minidx(end-(length(minsub)-1):end) = true;
            if timeVec(tnum)<0
                minstr(end-(length(minsub)-1):end) = minsub;
                minstr(~minidx) = ' ';
            else
                minstr(end-(length(minsub)-1):end) = minsub;
                
                minstr(~minidx)=' ';
            end
            
            frametext= text(ax1,0,0,[minstr ' hr']);
            frametext.Units = 'normalized';
            frametext.Position = [0.99 0.02];
            frametext.FontSize = 15;
            frametext.Color =[1 1 0];
            frametext.FontWeight = 'bold';
            frametext.HorizontalAlignment = 'right';
            
            if ~(timeVec(tnum)<0)
                frametext= text(ax1,0,0,'+TGFbeta');
                frametext.Units = 'normalized';
                frametext.Position = [0.99 0.05];
                frametext.FontSize = 15;
                frametext.Color =[1 1 0];
                frametext.FontWeight = 'bold';
                frametext.HorizontalAlignment = 'right';
            end
            
            
            stimulationFrame = find(~(timeVec<0),1,'first');
            
            set(ax1,'NextPlot','Add')
            ax1.XTick =[];
            ax1.YTick = [];
            %plot cell tracking "tails"
            if tnum>0
                
                trajectForPlot = trackingTrajectories(timeFrames);
                mainplotX = squeeze(trajectForPlot(:,1,:)); %28x22 means 28 cells on frame 22;
                mainplotY = squeeze(trajectForPlot(:,2,:));
                
                if size(mainplotY,2) == 1
                    mainplotY=mainplotY';
                    mainplotX=mainplotX';
                end
                %only plot if the cell is currently tracked/segmented in this frame
                idx = ~isnan(mainplotY(:,stimulationFrame));
                idxa = find(idx);
                
                cmapl = lines(size(mainplotY,1));
                plotcmap = zeros(length(idxa),3);
                for i=1:length(idxa)
                    plotcmap(i,:) =  cmapl(idxa(i),:);
                end
                %
                %             x = mainplotX(idxa,:) - cellcentroid(1) + cropsize;
                %             y = mainplotY(idxa,:) - cellcentroid(2) + cropsize;
                x = mainplotX(idxa,:) - cropymin;
                y = mainplotY(idxa,:) - cropxmin;
                t0 =max([nannum1-1 tnum-20]);
                %             h = plot(MainAxes,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',2);
                if togStruct.displayTrackingToggle==1
                    h = line(ax1,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',4);
                    set(h, {'color'}, num2cell(plotcmap,2));
                    nanan = isnan(y(:,tnum));
                    [h(nanan).Visible] = deal('off');
                end
            end
            set(ax1,'NextPlot','Replace')
            
            
            
            
            %plot traces
            mmsstr = medmeantot;
            %first axes, plot abundance
            fcplotstr = '';
            ylabelstr = 'abundance';
            if strcmp(ImageDetails.Channel,segStruct.background_seg)
                plotstr = ['Smad' fcplotstr];
            elseif strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
                plotstr = ['mkate' fcplotstr];
            else
                plotstr = ['Smad' fcplotstr];
            end
            plotMatOG = plotStructUI.(mmsstr).(plotstr);
            plotMatSub = plotMatOG(idx,:);
            % % % %     plotMat = nan(size(plotMatSub));
            % % % %     for j = 1:size(plotMatSub,1)
            % % % %         plotMat(j,:) = smooth(plotMatSub(j,:),0.05,'loess');
            % % % %     end
            % % % %     plotMat(isnan(plotMatSub)) = NaN;
            plotMat = medfilt1(plotMatSub,3,[],2);
            plotMat = plotMatSub./(nanmedian(plotMatSub(:))*0.01);
            plotMatStat = plotMat(:,nannum1:nannum2);
            ylimit = [prctile(plotMatStat(:),0)*0.95 prctile(plotMatStat(:),100)*1.2];
            
            h = plot(ax2,timeVec(nannum1:tnum),plotMat(:,nannum1:tnum)','LineWidth',2);
            if togStruct.displayTrackingToggle ==1
                for i=1:length(h)
                    h(i).Color = cmapl(i,:);
                end
                colormap(cmap);
            end
            
            if sum(diff(timeVec([nannum1 nannum2])))>0
                xlim = timeVec([nannum1 nannum2]);
            else
                xlim = [0 1];
            end
            
            ax2.XLim = xlim;
            ax2.YLim = (ylimit);
            ax2.YLabel.String = ['Smad3 ' ylabelstr ' ' mmsstr];
            ax2.YLabel.String = ['Smad3 ' ylabelstr];
            ax2.XLabel.String = 'minutes';
            ax2.XGrid = 'on';
            ax2.YGrid = 'on';
            ax2.Color = [0.95 0.95 0.95];
            ax2.Title.String = ['Nuclear mNG-Smad3 in tracked cells'];
            ax2.FontSize = 12;
            
            
            %second axes, plot fold-change
            
            fcplotstr = 'FC';
            ylabelstr = 'fold change';
            if strcmp(ImageDetails.Channel,segStruct.background_seg)
                plotstr = ['Smad' fcplotstr];
            elseif strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
                plotstr = ['mkate' fcplotstr];
            else
                plotstr = ['Smad' fcplotstr];
            end
            plotMatOG = plotStructUI.(mmsstr).(plotstr);
            plotMatSub = plotMatOG(idx,:);
            % % %     plotMat = nan(size(plotMatSub));
            % % %     for j = 1:size(plotMatSub,1)
            % % %         plotMat(j,:) = smooth(plotMatSub(j,:),0.05,'loess');
            % % %     end
            % % %     plotMat(isnan(plotMatSub)) = NaN;
            plotMat = medfilt1(plotMatSub,3,[],2);
            plotMat = plotMatSub;
            plotMatStat = plotMat(:,nannum1:nannum2);
            ylimit = [prctile(plotMatStat(:),0)*0.95 prctile(plotMatStat(:),100)*1.2];
            
            idx = true(size(plotMat,1),1);
            idxa = find(idx==1);
            h = plot(ax3,timeVec(nannum1:tnum),plotMat(idx,nannum1:tnum)','LineWidth',2);
            %update color traces to match image
            if togStruct.displayTrackingToggle ==1
                for i=1:length(h)
                    h(i).Color = cmapl(idxa(i),:);
                end
                colormap(cmap);
            end
            
            ax3.XLim = xlim;
            ax3.YLim = (ylimit);
            %     ax3.YLabel.String = ['Smad3 ' ylabelstr ' ' mmsstr];
            ax3.YLabel.String = ['Smad3 ' ylabelstr];
            ax3.XLabel.String = 'minutes';
            ax3.XGrid = 'on';
            ax3.YGrid = 'on';
            ax3.Color = [0.95 0.95 0.95];
            ax3.FontSize = 12;
            
            pause(0.15)
            drawnow
            
            
            if togStruct.SaveMovie
                specialdir = [dirStruct.parentdir 'LookingAtData/ExportImageFrames'];
                olddir = pwd;
                saveimg = channelimg;
                if ~isdir(specialdir)
                    mkdir(specialdir)
                end
                cd(specialdir)
                savedir = [specialdir '/' expDetailsStruct.expDateStr '/'];
                if ~isdir(savedir)
                    mkdir(savedir)
                end
                cd(savedir)
                savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.tif'];
                imwrite(uint16(saveimg),savename,'tif');
                
                %SAVE FIGURE
                %             cd(specialdir)
                %             %         savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' figure png frames'];
                %             savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' figure tiff frames'];
                %             if ~isdir(savedir)
                %                 mkdir(savedir)
                %             end
                %             cd(savedir)
                %             %         savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.png'];
                %             %         saveas(f33,savename,'png')
                %             savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.tif'];
                %             saveas(f33,savename,'tif')
                %             disp(num2str(ImageDetails.Frame));
                %             cd(olddir)
                
                cd(specialdir)
                savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' movieCode'];
                if ~isdir(savedir)
                    mkdir(savedir)
                end
                cd(savedir)
                savenamestr = ['uitrackcellz_' ImageDetails.Channel '_' ImageDetails.Scene '_' num2str(tnum) '.m'];
                mdir = pathStruct.mdirPath;
                copyfile([mdir '.m'],[savedir '/' savenamestr '.m'])
                cd(olddir)
            end
        end
    end
end
end
% 
% 
% function makeFCmovie(~,~)
% global stackstruct trackfilestr expDetailsStruct dirStruct timeVec pathStruct psettings medmeantot Tracked segStruct Ifstack channelimgstack nucleusimgstack DICimgstack togStruct timeFrames tcontrast lcontrast  ImageDetails prcntl lprcntl cmap cmaplz
% 
% bkgmedian=0;
% %LONG TRACK
% 
% 
% % click the desired cell
% ArrayStruct = Tracked.arrayStruct;
% trackmatrix = Tracked.trackmatrix;
% %choose the initial cell
% 
% 
% 
% 
% scenesTo = {'s02','s07','s11','s15','s19'}
% 
% %load traces to be plotted
% if togStruct.plotSettingsToggle == 0
%     psettings = PlotSettings_callback([],[]);
%     togStruct.plotSettingsToggle=1;
% end
% framesThatMustBeTracked = psettings.framesThatMustBeTracked;
% 
% 
% smooththat=0;
% [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,pathStruct,ImageDetails,timeFrames,smooththat);
% 
% 
% % set up figure
% f33 = figure(33);
% close 33
% 
% %%%% make movie by iterating through frames
% f33 = figure(33);
% f33.Color = 'w';
% f33.Units = 'pixels';
% f33.Position = [593 156 1210 760];
% 
% ax1 = axes();
% ax1.Units = 'pixels';
% ax1.Position = [25 25 720 720];
% 
% ax2 = axes();
% ax2.Units = 'pixels';
% ax2.Position = [830 445 335 285];
% 
% ax3 = axes();
% ax3.Units = 'pixels';
% ax3.Position = [830 80 335 285];
% 
% f33.Units = 'normalized';
% h = findobj(f33,'Type','Axes');
% [h.Units] = deal('normalized');
% pos = f33.Position;
% pos(1:2) = [0.05 0.05];
% pos(3:4) = pos(3:4)./1.2;
% f33.Position = pos;
% 
% 
% cropbuffer = 20;
% imgdim = ImageDetails.ImgSize;
% celltrack = ~isnan(trackmatrix);
% nanvec = sum((celltrack),2);
% % nannum1 = max([2 find(nanvec == max(nanvec),1,'first')]);
% nannum1 = max([2 find(nanvec > 0,1,'first')]);
% % nannum1 = 1;
% 
% % nannum2 = find(nanvec == max(nanvec),1,'last');
% nannum2 = find(nanvec > 0,1,'last');
% segment_Centroids_array = ArrayStruct.centroid;
% 
% %determine crop region
% xvals = nan(nannum2,2);
% yvals = nan(nannum2,2);
% cycle = 0;
% for tnum = nannum1:nannum2
%     cycle = cycle+1;
%     cellindices = trackmatrix(tnum,:);
%     centroids = segment_Centroids_array{tnum};
%     cellcentroids = centroids(cellindices(~isnan(cellindices)),:);
%     yvals(cycle,:) = round([min(cellcentroids(:,1)) max(cellcentroids(:,1))]);
%     xvals(cycle,:) = round([min(cellcentroids(:,2)) max(cellcentroids(:,2))]);
% end
% cropxmin = max([min(xvals(:,1))-cropbuffer 1]);
% cropxmax = min([max(xvals(:,2))+cropbuffer imgdim(1)]);
% cropymin = max([min(yvals(:,1))-cropbuffer 1]);
% cropymax = min([max(yvals(:,2))+cropbuffer imgdim(2)]);
% xdist = cropxmax-cropxmin;
% ydist = cropymax-cropymin;
% 
% ddiff = abs(ydist-xdist);
% dd2x = floor(ddiff/2);
% dd2y = ceil(ddiff/2);
% if ddiff>0
%     if ydist>xdist
%         x1 = cropxmin-1;
%         x2 = 512-cropxmax;
%         if x1 < dd2x
%             dd2x = x1;
%             dd2y = ddiff-x1;
%         elseif x2 < dd2y
%             dd2y = x2;
%             dd2x = ddiff-x2;
%         end
%         cropxmin = cropxmin-dd2x;
%         cropxmax = cropxmax+dd2y;
%     else
%         x1 = cropymin-1;
%         x2 = 512-cropymax;
%         if x1 < dd2x
%             dd2x = x1;
%             dd2y = ddiff-x1;
%         elseif x2 < dd2y
%             dd2y = x2;
%             dd2x = ddiff-x2;
%         end
%         cropymin = cropymin-dd2x;
%         cropymax = cropymax+dd2y;
%     end
% end
% 
% xdist = cropxmax-cropxmin;
% ydist = cropymax-cropymin;
% 
% 
% for tnum = nannum1:nannum2
%     t = tnum;
%     
%     If = Ifstack(:,:,t);
%     CC = bwconncomp(If);
%     px = CC.PixelIdxList;
%     x = mainplotX(:,:);
%     y = mainplotY(:,:);
%     idx = ~isnan(mainplotY(:,tnum));
%     idxa = find(idx==1);
%     xyind = sub2ind(size(If),round(x(:,tnum)),round(y(:,tnum)));
%     ii=1;
%     for ii = 1:length(px)
%         midx = sum(ismember(px{ii},xyind));
%         if midx>0
%             fcval = fc(idx,:);
%             Ifn(px{ii}) = fcval;
%         end
%     end
%         
%     end
%     
%     [fc,abb,pidx] = Plot_callback([],[]);
%     for ui = 1:size(fc,1)
%         ui;
%         
%     end
%     plotcmap = flipud(feval('viridis',length(sidx)));
%     
%     trajectForPlot = trackingTrajectories(timeFrames);
%     mainplotX = squeeze(trajectForPlot(:,1,:)); %28x22 means 28 cells on frame 22;
%     mainplotY = squeeze(trajectForPlot(:,2,:));
%     if size(mainplotY,2) == 1
%         mainplotY=mainplotY';
%         mainplotX=mainplotX';
%     end
%     
%     %only plot if the cell is currently tracked/segmented in this frame
% 
%     if togStruct.trackUpdated
%         cmaplz = colorcubemodified(length(idx),'jet');
%     elseif ~(size(cmaplz,1)==length(idx))
%         cmaplz = colorcubemodified(length(idx),'jet');
%     end
%     cmapl = cmaplz;
%     plotcmap = zeros(length(idxa),3);
%     for i=1:length(idxa)
%         plotcmap(i,:) =  cmapl(idxa(i),:);
%     end
%     
%     
%     
%     x = mainplotX(idx,:);
%     y = mainplotY(idx,:);
%     t0 =max([1 tnum-20]);
%     h = line(MainAxes,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',2);
%     set(h, {'color'}, num2cell(plotcmap,2));
%     % only plot if the cell is currently tracked/segmented in this frame
% 
%     
%     nucleusImg = stackstruct.(segStruct.nucleus_seg)(:,:,t);
%     cellImg = stackstruct.(segStruct.background_seg)(:,:,t);
%             
%     
%     
%     if strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
%         channelimg = nucleusImg;
%     elseif strcmp(ImageDetails.Channel,segStruct.background_seg)
%         channelimg = cellImg;
%     elseif strcmp(ImageDetails.Channel,'overlay')
%         DICImg = DICimgstack(:,:,t);
%         channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
%         channelimg(:,:,1) = nucleusImg;
% %         channelimg(:,:,1) = zeros(size(nucleusImg));
%         channelimg(:,:,2) = cellImg;
%         channelimg(:,:,3) = DICImg;
%     elseif strcmpi(ImageDetails.Channel,'FluorOnlyOverlay')
%         DICImg = DICimgstack(:,:,t);
%         channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
%         channelimg(:,:,1) = nucleusImg;
%         channelimg(:,:,2) = cellImg;
%     elseif strcmpi(ImageDetails.Channel,'647')
%         channelimg = DICimgstack(:,:,t);
%     elseif strcmpi(ImageDetails.Channel,'DIC')
%         channelimg = stackstruct.('DIC')(:,:,t);
%     end
%     
%     if length(size(channelimg))>2
%         cropim = channelimg(cropxmin:cropxmax,cropymin:cropymax,:);
%     else
%         cropim = channelimg(cropxmin:cropxmax,cropymin:cropymax);
%     end
%     If = If(cropxmin:cropxmax,cropymin:cropymax);
%     channelimg = cropim;
%     
%     If(If<2)=0;
%     If(If>2) = 100;
% %     If(If>1) = 0;
%     
%     
%     
%     
%     
%     
%     %scripts for displaying contrasted image
%     alog = regexpi(ImageDetails.Channel,'overlay'); %when overlay display is desired
%     if ~isempty(alog)
%             %scripts for displaying contrasted image
%             disprgb = zeros(size(cropim));
%             channelimgrgb = cropim;
%             Ifrgb = ind2rgb(uint8(If),vertcat([0 0 0],cmap(257:end,:)));
%             for i = 1:size(cropim,3)
%                 if i==3
%                     channelimg = channelimgrgb(:,:,i);
%                     lprcntl = prctile(channelimg(:),lcontrast);
%                     prcntl = prctile(channelimg(:),tcontrast);
%                     scaleFactor = 1./(prcntl - lprcntl);
%                     dispimg = channelimg.*scaleFactor;
%                     dispimg = dispimg-(lprcntl.*scaleFactor);
%                     dispimg(dispimg> 1) =1;
%                     dispimg(dispimg<0) = 0;
%                 else
%                     channelimg = channelimgrgb(:,:,i);
%                     lprcntl = 0;
%                     prcntl = prctile(channelimg(:),tcontrast);
%                     scaleFactor = 1./(prcntl - lprcntl);
%                     dispimg = channelimg.*scaleFactor;
%                     dispimg = dispimg-(lprcntl.*scaleFactor);
%                     dispimg(dispimg> 1) =1;
%                     dispimg(dispimg<0) = 0;
%                 end
%                 
%                 
%                 if i==1
%                     disprgb(:,:,i) = dispimg;
%                 elseif i==2
%                     disprgb(:,:,i) = dispimg;
%                 elseif i==3
%                     OGrgb = disprgb;
%                     
%                     DICrgb = channelimgrgb;
%                     for j = 1:3
%                         cimg =  disprgb(:,:,j);
%                         cimgpre = cimg;
%                         %                 cimg(dispimg>(cimg./a)) = dispimg(dispimg>(cimg./a));
%                         cimg = cimg.*dispimg;
%                         disprgb(:,:,j) = cimg;
%                         DICrgb(:,:,j) = dispimg;
%                     end
%                     cimgmax =  max(disprgb,[],3);
%                     a = 0.1;
%                     for j = 1:3
%                         cimg =  disprgb(:,:,j);
%                         cimg = cimg + dispimg./5;
%                         cimg(cimg>1) = 1;
%                         disprgb(:,:,j) = cimg;
%                     end
%                     
%                 end
%             end
%         
%         newimg = zeros(size(disprgb));
%         for i = 1:size(disprgb,3)
%             Ifsub = Ifrgb(:,:,i);
%             Ifnew = zeros(size(Ifsub));
%             Ifnew(If>0) = Ifsub(If>0);
% %             Ifnew(If>0) = 0;
%             dispimg = disprgb(:,:,i);
%             dispimg(If>0) = Ifsub(If>0);
% %             dispimg(If>0) = 0;
%             disprgb(:,:,i) = dispimg;
%             newimg(:,:,i) = Ifnew;
%         end
%         
%         
%         dispimg = disprgb;
%         %     dispimg = uint8(dispimg);
%         image(ax1,dispimg);
%         %     colormap(MainAxes,cmap);
%         %     set(MainAxes,'CLim',[0 size(cmap,1)]);
%         
%     else  %under normal circumstances
%         if togStruct.changeSceneOrChannel ==1
%             prcntl = prctile(channelimg(:),tcontrast);
%             togStruct.changeSceneOrChannel =0;
%         end
%         %     lprcntl = bkgmedian.*0.90;
%         lprcntl = bkgmedian;
%         scaleFactor = 255./(prcntl - lprcntl);
%         dispimg = channelimg.*scaleFactor;
%         dispimg = dispimg-(lprcntl.*scaleFactor);
%         dispimg(dispimg > 250) =253;
%         dispimg(dispimg < 0) = 0;
% %         dispimg(If>0)=255;
%         dispimg(If>0)=If(If>0)+257;
%         imagesc(ax1,int16(dispimg));
%         colormap(ax1,cmap);
%         set(ax1,'CLim',[0 size(cmap,1)]);
%     end
%     
%     %title the displayed image
%     
%     
%     
% %     frametext= text(ax1,0,0,['frame ' num2str(tnum) '/' num2str(timeFrames)]);
% %     frametext.Units = 'normalized';
% %     frametext.Position = [0.01 0.95];
% %     frametext.FontSize = 12;
% %     frametext.Color =[1 1 0];
%     minstr = '0000';
%     minsub = num2str(timeVec(tnum));
%     minidx = false(size(minstr));
%     minidx(end-(length(minsub)-1):end) = true;
%     if timeVec(tnum)<0
%         minstr(end-(length(minsub)-1):end) = minsub;
%         minstr(~minidx) = ' ';
%     else
%         minstr(end-(length(minsub)-1):end) = minsub;
% 
%         minstr(~minidx)=' ';
%     end
%     
%     frametext= text(ax1,0,0,[minstr ' min']);
%     frametext.Units = 'normalized';
%     frametext.Position = [0.99 0.02];
%     frametext.FontSize = 15;
%     frametext.Color =[1 1 0];
%     frametext.FontWeight = 'bold';
%     frametext.HorizontalAlignment = 'right';
%     
%     if ~(timeVec(tnum)<0)
%         frametext= text(ax1,0,0,'+TGFbeta');
%         frametext.Units = 'normalized';
%         frametext.Position = [0.99 0.05];
%         frametext.FontSize = 15;
%         frametext.Color =[1 1 0];
%         frametext.FontWeight = 'bold';
%         frametext.HorizontalAlignment = 'right';
%     end
%     
%     
%     stimulationFrame = find(~(timeVec<0),1,'first');
%     
%     set(ax1,'NextPlot','Add')
%     ax1.XTick =[];
%     ax1.YTick = [];
%     %plot cell tracking "tails"
%     if tnum>0
%         
%         trajectForPlot = trackingTrajectories(timeFrames);
%         mainplotX = squeeze(trajectForPlot(:,1,:)); %28x22 means 28 cells on frame 22;
%         mainplotY = squeeze(trajectForPlot(:,2,:));
%         
%         if size(mainplotY,2) == 1
%             mainplotY=mainplotY';
%             mainplotX=mainplotX';
%         end
%         %only plot if the cell is currently tracked/segmented in this frame
%         idx = ~isnan(mainplotY(:,stimulationFrame));
%         idxa = find(idx);
%         
%         cmapl = lines(size(mainplotY,1));
%         plotcmap = zeros(length(idxa),3);
%         for i=1:length(idxa)
%             plotcmap(i,:) =  cmapl(idxa(i),:);
%         end
%         %
%         %             x = mainplotX(idxa,:) - cellcentroid(1) + cropsize;
%         %             y = mainplotY(idxa,:) - cellcentroid(2) + cropsize;
%         x = mainplotX(idxa,:) - cropymin;
%         y = mainplotY(idxa,:) - cropxmin;
%         t0 =max([nannum1-1 tnum-20]);
%         %             h = plot(MainAxes,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',2);
%         if togStruct.displayTrackingToggle==1
%             h = line(ax1,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',4);
%             set(h, {'color'}, num2cell(plotcmap,2));
%             nanan = isnan(y(:,tnum));
%             [h(nanan).Visible] = deal('off');
%         end
%     end
%     set(ax1,'NextPlot','Replace')
%     
%     
%     
%     
%     %plot traces
%     mmsstr = medmeantot;
%     %first axes, plot abundance
%     fcplotstr = '';
%     ylabelstr = 'abundance';
%     if strcmp(ImageDetails.Channel,segStruct.background_seg)
%         plotstr = ['Smad' fcplotstr];
%     elseif strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
%         plotstr = ['mkate' fcplotstr];
%     else
%         plotstr = ['Smad' fcplotstr];
%     end
%     plotMatOG = plotStructUI.(mmsstr).(plotstr);
%     plotMatSub = plotMatOG(idx,:);
% % % % %     plotMat = nan(size(plotMatSub));
% % % % %     for j = 1:size(plotMatSub,1)
% % % % %         plotMat(j,:) = smooth(plotMatSub(j,:),0.05,'loess');
% % % % %     end
% % % % %     plotMat(isnan(plotMatSub)) = NaN;
%     plotMat = medfilt1(plotMatSub,3,[],2);
%     plotMat = plotMatSub./(nanmedian(plotMatSub(:))*0.01);
%     plotMatStat = plotMat(:,nannum1:nannum2);
%     ylimit = [prctile(plotMatStat(:),0)*0.95 prctile(plotMatStat(:),100)*1.2];
%     ylimit = [0 prctile(plotMatStat(:),100)*1.2];
%     
%     h = plot(ax2,timeVec(nannum1:tnum),plotMat(:,nannum1:tnum)','LineWidth',2);
%     if togStruct.displayTrackingToggle ==1
%         for i=1:length(h)
%             h(i).Color = cmapl(i,:);
%         end
%         colormap(cmap);
%     end
%     
%     
%     if sum(diff(timeVec([nannum1 nannum2])))==0
%         xlim = [0 1];
%     else
%         xlim = timeVec([nannum1 nannum2]);
%     end
%         
%         
%     ax2.XLim = xlim;
%     ax2.YLim = (ylimit);
%     ax2.YLabel.String = ['Smad3 ' ylabelstr ' ' mmsstr];
%     ax2.YLabel.String = ['Smad3 ' ylabelstr];
%     ax2.XLabel.String = 'minutes';
%     ax2.XGrid = 'on';
%     ax2.YGrid = 'on';
%     ax2.Color = [0.95 0.95 0.95];
%     ax2.Title.String = ['Nuclear mNG-Smad3 in tracked cells'];
%     ax2.FontSize = 12;
%     
%     
%     %second axes, plot fold-change
%     
%     fcplotstr = 'FC';
%     ylabelstr = 'fold change';
%     if strcmp(ImageDetails.Channel,segStruct.background_seg)
%         plotstr = ['Smad' fcplotstr];
%     elseif strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
%         plotstr = ['mkate' fcplotstr];
%     else
%         plotstr = ['Smad' fcplotstr];
%     end
%     plotMatOG = plotStructUI.(mmsstr).(plotstr);
%     plotMatSub = plotMatOG(idx,:);
% % % %     plotMat = nan(size(plotMatSub));
% % % %     for j = 1:size(plotMatSub,1)
% % % %         plotMat(j,:) = smooth(plotMatSub(j,:),0.05,'loess');
% % % %     end
% % % %     plotMat(isnan(plotMatSub)) = NaN;
%     plotMat = medfilt1(plotMatSub,3,[],2);
%     plotMat = plotMatSub;
%     plotMatStat = plotMat(:,nannum1:nannum2);
%     ylimit = [prctile(plotMatStat(:),0)*0.95 prctile(plotMatStat(:),100)*1.2];
%     ylimit = [0 prctile(plotMatStat(:),100)*1.2];
%     
%     idx = true(size(plotMat,1),1);
%     idxa = find(idx==1);
%     h = plot(ax3,timeVec(nannum1:tnum),plotMat(idx,nannum1:tnum)','LineWidth',2);
%     %update color traces to match image
%     if togStruct.displayTrackingToggle ==1
%         for i=1:length(h)
%             h(i).Color = cmapl(idxa(i),:);
%         end
%         colormap(cmap);
%     end
%     
%     ax3.XLim = xlim;
%     ax3.YLim = (ylimit);
% %     ax3.YLabel.String = ['Smad3 ' ylabelstr ' ' mmsstr];
%     ax3.YLabel.String = ['Smad3 ' ylabelstr];
%     ax3.XLabel.String = 'minutes';
%     ax3.XGrid = 'on';
%     ax3.YGrid = 'on';
%     ax3.Color = [0.95 0.95 0.95];
%     ax3.FontSize = 12;
%     
%     
%     
%     if togStruct.SaveMovie
%     else
%     pause(0.15)
%     end
%     drawnow
%     
%     
%     if togStruct.SaveMovie
%         specialdir = [dirStruct.parentdir 'LookingAtData/AnimatedTrackingPlottingMovie'];
%         olddir = pwd;
%         saveimg = channelimg;
%         if ~isdir(specialdir)
%             mkdir(specialdir)
%         end
%         cd(specialdir)
%         savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' tif frames'];
%         if ~isdir(savedir)
%             mkdir(savedir)
%         end
%         cd(savedir)
%         savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.tif'];
%         imwrite(uint16(saveimg),savename,'tif');
%         
%         cd(specialdir)
% %         savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' figure png frames'];
%         savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' figure tiff frames'];
%         if ~isdir(savedir)
%             mkdir(savedir)
%         end
%         cd(savedir)
% %         savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.png'];
% %         saveas(f33,savename,'png')
%         savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.tif'];
%         saveas(f33,savename,'tif')
%         disp(num2str(ImageDetails.Frame));
%         cd(olddir)
%         
%         cd(specialdir)
%         savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' movieCode'];
%         if ~isdir(savedir)
%             mkdir(savedir)
%         end
%         cd(savedir)
%         savenamestr = ['uitrackcellz_' ImageDetails.Channel '_' ImageDetails.Scene '_' num2str(tnum) '.m'];
%         mdir = pathStruct.mdirPath;
%         copyfile([mdir '.m'],[savedir '/' savenamestr '.m'])
%         cd(olddir)
%     end
% end
% end
% 


%excellent
function makeMovieCroppedYEAH(~,~)
global stackstruct trackfilestr expDetailsStruct dirStruct timeVec pathStruct psettings medmeantot Tracked segStruct Ifstack channelimgstack nucleusimgstack DICimgstack togStruct timeFrames tcontrast lcontrast  ImageDetails prcntl lprcntl cmap cmaplz

bkgmedian=0;
%LONG TRACK


% click the desired cell
ArrayStruct = Tracked.arrayStruct;
trackmatrix = Tracked.trackmatrix;
%choose the initial cell




scenesTo = {'s02','s07','s11','s15','s19'}

%load traces to be plotted
if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
end
framesThatMustBeTracked = psettings.framesThatMustBeTracked;


smooththat=0;
[plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,pathStruct,ImageDetails,timeFrames,smooththat);


% set up figure
f33 = figure(33);
close 33

%%%% make movie by iterating through frames
f33 = figure(33);
f33.Color = 'w';
f33.Units = 'pixels';
f33.Position = [593 156 1210 760];

ax1 = axes();
ax1.Units = 'pixels';
ax1.Position = [25 25 500 500];
% ax1.Position = [25 25 720 720];

ax2 = axes();
ax2.Units = 'pixels';
% ax2.Position = [830 445 335 285];
ax2.Position = [830 1000 1000 285];

ax3 = axes();
ax3.Units = 'pixels';
% ax3.Position = [830 80 335 285];
ax3.Position = [650 100 400 350];

f33.Units = 'normalized';
h = findobj(f33,'Type','Axes');
[h.Units] = deal('normalized');
pos = f33.Position;
pos(1:2) = [0.05 0.05];
pos(3:4) = pos(3:4)./1.2;
f33.Position = pos;


cropbuffer = 20;
imgdim = ImageDetails.ImgSize;
celltrack = ~isnan(trackmatrix);
nanvec = sum((celltrack),2);
% nannum1 = max([2 find(nanvec == max(nanvec),1,'first')]);
nannum1 = max([2 find(nanvec > 0,1,'first')]);
% nannum1 = 1;

% nannum2 = find(nanvec == max(nanvec),1,'last');
nannum2 = find(nanvec > 0,1,'last');
segment_Centroids_array = ArrayStruct.centroid;

%determine crop region
xvals = nan(nannum2,2);
yvals = nan(nannum2,2);
cycle = 0;
for tnum = nannum1:nannum2
    cycle = cycle+1;
    cellindices = trackmatrix(tnum,:);
    centroids = segment_Centroids_array{tnum};
    cellcentroids = centroids(cellindices(~isnan(cellindices)),:);
    yvals(cycle,:) = round([min(cellcentroids(:,1)) max(cellcentroids(:,1))]);
    xvals(cycle,:) = round([min(cellcentroids(:,2)) max(cellcentroids(:,2))]);
end
cropxmin = max([min(xvals(:,1))-cropbuffer 1]);
cropxmax = min([max(xvals(:,2))+cropbuffer imgdim(1)]);
cropymin = max([min(yvals(:,1))-cropbuffer 1]);
cropymax = min([max(yvals(:,2))+cropbuffer imgdim(2)]);
xdist = cropxmax-cropxmin;
ydist = cropymax-cropymin;

ddiff = abs(ydist-xdist);
dd2x = floor(ddiff/2);
dd2y = ceil(ddiff/2);
if ddiff>0
    if ydist>xdist
        x1 = cropxmin-1;
        x2 = 512-cropxmax;
        if x1 < dd2x
            dd2x = x1;
            dd2y = ddiff-x1;
        elseif x2 < dd2y
            dd2y = x2;
            dd2x = ddiff-x2;
        end
        cropxmin = cropxmin-dd2x;
        cropxmax = cropxmax+dd2y;
    else
        x1 = cropymin-1;
        x2 = 512-cropymax;
        if x1 < dd2x
            dd2x = x1;
            dd2y = ddiff-x1;
        elseif x2 < dd2y
            dd2y = x2;
            dd2x = ddiff-x2;
        end
        cropymin = cropymin-dd2x;
        cropymax = cropymax+dd2y;
    end
end

xdist = cropxmax-cropxmin;
ydist = cropymax-cropymin;

for tnum = nannum1:nannum2
    t = tnum;
    
    If = Ifstack(:,:,t);
    nucleusImg = stackstruct.(segStruct.nucleus_seg)(:,:,t);
    cellImg = stackstruct.(segStruct.background_seg)(:,:,t);
            
    
    
    if strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
        channelimg = nucleusImg;
    elseif strcmp(ImageDetails.Channel,segStruct.background_seg)
        channelimg = cellImg;
    elseif strcmp(ImageDetails.Channel,'overlay')
        DICImg = stackstruct.('DIC')(:,:,t);
        channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
        channelimg(:,:,1) = nucleusImg;
%         channelimg(:,:,1) = zeros(size(nucleusImg));
        channelimg(:,:,2) = cellImg;
        channelimg(:,:,3) = DICImg;
    elseif strcmpi(ImageDetails.Channel,'FluorOnlyOverlay')
        DICImg = stackstruct.('DIC')(:,:,t);
        channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
        channelimg(:,:,1) = nucleusImg;
        channelimg(:,:,2) = cellImg;
    elseif strcmpi(ImageDetails.Channel,'647')
        channelimg = stackstruct.('DIC')(:,:,t);
    elseif strcmpi(ImageDetails.Channel,'DIC')
        channelimg = stackstruct.('DIC')(:,:,t);
    end
    
    if length(size(channelimg))>2
        cropim = channelimg(cropxmin:cropxmax,cropymin:cropymax,:);
    else
        cropim = channelimg(cropxmin:cropxmax,cropymin:cropymax);
    end
    If = If(cropxmin:cropxmax,cropymin:cropymax);
    channelimg = cropim;
    
    If(If<2)=0;
    If(If>2) = 100;
%     If(If>1) = 0;
    
    
    
    
    
    
    %scripts for displaying contrasted image
    alog = regexpi(ImageDetails.Channel,'overlay'); %when overlay display is desired
    if ~isempty(alog)
            %scripts for displaying contrasted image
            disprgb = zeros(size(cropim));
            channelimgrgb = cropim;
            Ifrgb = ind2rgb(uint8(If),vertcat([0 0 0],cmap(257:end,:)));
            for i = 1:size(cropim,3)
                if i==3
                    channelimg = channelimgrgb(:,:,i);
                    lprcntl = prctile(channelimg(:),lcontrast);
                    prcntl = prctile(channelimg(:),tcontrast);
                    scaleFactor = 1./(prcntl - lprcntl);
                    dispimg = channelimg.*scaleFactor;
                    dispimg = dispimg-(lprcntl.*scaleFactor);
                    dispimg(dispimg> 1) =1;
                    dispimg(dispimg<0) = 0;
                else
                    channelimg = channelimgrgb(:,:,i);
                    lprcntl = 0;
                    prcntl = prctile(channelimg(:),tcontrast);
                    scaleFactor = 1./(prcntl - lprcntl);
                    dispimg = channelimg.*scaleFactor;
                    dispimg = dispimg-(lprcntl.*scaleFactor);
                    dispimg(dispimg> 1) =1;
                    dispimg(dispimg<0) = 0;
                end
                
                
                if i==1
                    disprgb(:,:,i) = dispimg;
                elseif i==2
                    disprgb(:,:,i) = dispimg;
                elseif i==3
                    OGrgb = disprgb;
                    
                    DICrgb = channelimgrgb;
                    for j = 1:3
                        cimg =  disprgb(:,:,j);
                        cimgpre = cimg;
                        %                 cimg(dispimg>(cimg./a)) = dispimg(dispimg>(cimg./a));
                        cimg = cimg.*dispimg;
                        disprgb(:,:,j) = cimg;
                        DICrgb(:,:,j) = dispimg;
                    end
                    cimgmax =  max(disprgb,[],3);
                    a = 0.1;
                    for j = 1:3
                        cimg =  disprgb(:,:,j);
                        cimg = cimg + dispimg./5;
                        cimg(cimg>1) = 1;
                        disprgb(:,:,j) = cimg;
                    end
                    
                end
            end
        
        newimg = zeros(size(disprgb));
        for i = 1:size(disprgb,3)
            Ifsub = Ifrgb(:,:,i);
            Ifnew = zeros(size(Ifsub));
            Ifnew(If>0) = Ifsub(If>0);
%             Ifnew(If>0) = 0;
            dispimg = disprgb(:,:,i);
            dispimg(If>0) = Ifsub(If>0);
%             dispimg(If>0) = 0;
            disprgb(:,:,i) = dispimg;
            newimg(:,:,i) = Ifnew;
        end
        
        
        dispimg = disprgb;
        %     dispimg = uint8(dispimg);
        image(ax1,dispimg);
        %     colormap(MainAxes,cmap);
        %     set(MainAxes,'CLim',[0 size(cmap,1)]);
        
    else  %under normal circumstances
        if togStruct.changeSceneOrChannel ==1
            prcntl = prctile(channelimg(:),tcontrast);
            togStruct.changeSceneOrChannel =0;
        end
        %     lprcntl = bkgmedian.*0.90;
        lprcntl = bkgmedian;
        scaleFactor = 255./(prcntl - lprcntl);
        dispimg = channelimg.*scaleFactor;
        dispimg = dispimg-(lprcntl.*scaleFactor);
        dispimg(dispimg > 250) =253;
        dispimg(dispimg < 0) = 0;
%         dispimg(If>0)=255;
%         dispimg(If>0)=If(If>0)+257;
        imagesc(ax1,int16(dispimg));
        colormap(ax1,cmap);
        set(ax1,'CLim',[0 size(cmap,1)]);
    end
    
    %title the displayed image
    
    
    
%     frametext= text(ax1,0,0,['frame ' num2str(tnum) '/' num2str(timeFrames)]);
%     frametext.Units = 'normalized';
%     frametext.Position = [0.01 0.95];
%     frametext.FontSize = 12;
%     frametext.Color =[1 1 0];
if timeVec(tnum)<0
    hr = ceil(timeVec(tnum)./60);
    mint = abs(round(timeVec(tnum) - hr*60));
else
    hr = floor(timeVec(tnum)./60);
    mint = round(timeVec(tnum) - hr*60);
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
    else
        hrstr(end-(length(minsub)-1):end) = minsub;

%         hrstr(~minidx)='';
    end
    
    
    frametext= text(ax1,0,0,[hrstr ':' minstr ' hr:min']);
    frametext.Units = 'normalized';
    frametext.Position = [0.99 0.05];
    frametext.FontSize = 26;
    frametext.Color =[1 1 0];
    frametext.FontWeight = 'bold';
    frametext.HorizontalAlignment = 'right';
    frametext.BackgroundColor = 'k';
    
    if ~(timeVec(tnum)<0)
        frametext= text(ax1,0,0,'+TGFbeta');
        frametext.Units = 'normalized';
        frametext.Position = [0.99 0.12];
        frametext.FontSize = 26;
        frametext.Color =[1 1 0];
        frametext.FontWeight = 'bold';
        frametext.HorizontalAlignment = 'right';
        frametext.BackgroundColor = 'k';
    end
    
    
    stimulationFrame = find(~(timeVec<0),1,'first');
    
    set(ax1,'NextPlot','Add')
    ax1.XTick =[];
    ax1.YTick = [];
    %plot cell tracking "tails"
    if tnum>0
        
        trajectForPlot = trackingTrajectories(timeFrames);
        mainplotX = squeeze(trajectForPlot(:,1,:)); %28x22 means 28 cells on frame 22;
        mainplotY = squeeze(trajectForPlot(:,2,:));
        
        if size(mainplotY,2) == 1
            mainplotY=mainplotY';
            mainplotX=mainplotX';
        end
        %only plot if the cell is currently tracked/segmented in this frame
        idx = ~isnan(mainplotY(:,stimulationFrame));
        idxa = find(idx);
        
        cmapl = lines(size(mainplotY,1));
        plotcmap = zeros(length(idxa),3);
        for i=1:length(idxa)
            plotcmap(i,:) =  cmapl(idxa(i),:);
        end
        %
        %             x = mainplotX(idxa,:) - cellcentroid(1) + cropsize;
        %             y = mainplotY(idxa,:) - cellcentroid(2) + cropsize;
        x = mainplotX(idxa,:) - cropymin;
        y = mainplotY(idxa,:) - cropxmin;
        t0 =max([nannum1-1 tnum-20]);
        %             h = plot(MainAxes,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',2);
        if togStruct.displayTrackingToggle==1
            h = line(ax1,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',4);
            set(h, {'color'}, num2cell(plotcmap,2));
            nanan = isnan(y(:,tnum));
            [h(nanan).Visible] = deal('off');
        end
    end
    set(ax1,'NextPlot','Replace')
    
    
    
    
    %plot traces
    mmsstr = medmeantot;
    %first axes, plot abundance
    fcplotstr = '';
    ylabelstr = 'abundance';
    if strcmp(ImageDetails.Channel,segStruct.background_seg)
        plotstr = ['Smad' fcplotstr];
    elseif strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
        plotstr = ['mkate' fcplotstr];
    else
        plotstr = ['Smad' fcplotstr];
    end
    plotMatOG = plotStructUI.(mmsstr).(plotstr);
    plotMatSub = plotMatOG(idx,:);
% % % %     plotMat = nan(size(plotMatSub));
% % % %     for j = 1:size(plotMatSub,1)
% % % %         plotMat(j,:) = smooth(plotMatSub(j,:),0.05,'loess');
% % % %     end
% % % %     plotMat(isnan(plotMatSub)) = NaN;
    plotMat = medfilt1(plotMatSub,3,[],2);
    plotMat = plotMatSub./(nanmedian(plotMatSub(:))*0.01);
    plotMatStat = plotMat(:,nannum1:nannum2);
    ylimit = [prctile(plotMatStat(:),0)*0.95 prctile(plotMatStat(:),100)*1.2];
    ylimit = [0 prctile(plotMatStat(:),100)*1.2];
    
    h = plot(ax2,timeVec(nannum1:tnum),plotMat(:,nannum1:tnum)','LineWidth',4);
    if togStruct.displayTrackingToggle ==1
        for i=1:length(h)
            h(i).Color = cmapl(i,:);
        end
        colormap(cmap);
    end
    
    
    if sum(diff(timeVec([nannum1 nannum2])))==0
        xlim = [0 1];
    else
        xlim = timeVec([nannum1 nannum2]);
    end
        
        fSize = 24;
    ax2.XLim = xlim;
    ax2.YLim = (ylimit);
    ax2.YLabel.String = ['Smad3 ' ylabelstr ' ' mmsstr];
    ax2.YLabel.String = ['Smad3 ' ylabelstr];
    ax2.XLabel.String = 'minutes';
    ax2.XGrid = 'on';
    ax2.YGrid = 'on';
    ax2.Color = [0.95 0.95 0.95];
    ax2.Title.String = ['NG-Smad3 response'];
    ax2.FontSize = fSize;
    
    
    %second axes, plot fold-change
    
    fcplotstr = 'FC';
    ylabelstr = 'fold change';
    if strcmp(ImageDetails.Channel,segStruct.background_seg)
        plotstr = ['Smad' fcplotstr];
    elseif strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
        plotstr = ['mkate' fcplotstr];
    else
        plotstr = ['Smad' fcplotstr];
    end
    plotMatOG = plotStructUI.(mmsstr).(plotstr);
    plotMatSub = plotMatOG(idx,:);
% % %     plotMat = nan(size(plotMatSub));
% % %     for j = 1:size(plotMatSub,1)
% % %         plotMat(j,:) = smooth(plotMatSub(j,:),0.05,'loess');
% % %     end
% % %     plotMat(isnan(plotMatSub)) = NaN;
    plotMat = medfilt1(plotMatSub,3,[],2);
    plotMat = plotMatSub;
    plotMatStat = plotMat(:,nannum1:nannum2);
    ylimit = [prctile(plotMatStat(:),0)*0.95 prctile(plotMatStat(:),100)*1.2];
    ylimit = [0 prctile(plotMatStat(:),100)*1.2];
    
    idx = true(size(plotMat,1),1);
    idxa = find(idx==1);
    h = plot(ax3,timeVec(nannum1:tnum),plotMat(idx,nannum1:tnum)','LineWidth',4);
    %update color traces to match image
    if togStruct.displayTrackingToggle ==1
        for i=1:length(h)
            h(i).Color = cmapl(idxa(i),:);
        end
        colormap(cmap);
    end
    
    ax3.XLim = xlim;
    ax3.YLim = (ylimit);
%     ax3.YLabel.String = ['Smad3 ' ylabelstr ' ' mmsstr];
    ax3.YLabel.String = ['Smad3 ' ylabelstr];
    ax3.XLabel.String = 'minutes';
    ax3.XGrid = 'on';
    ax3.YGrid = 'on';
    ax3.Color = [0.95 0.95 0.95];
    ax3.Title.String = ['NG-Smad3 response'];
    ax3.FontSize = fSize;
    
    
    
    if togStruct.SaveMovie
    else
    pause(0.15)
    end
    drawnow
    
    
    f33.InvertHardcopy = 'off';
    if togStruct.SaveMovie
        specialdir = [dirStruct.parentdir 'LookingAtData/AnimatedTrackingPlottingMovie'];
        olddir = pwd;
        saveimg = channelimg;
        if ~isdir(specialdir)
            mkdir(specialdir)
        end
        cd(specialdir)
        savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' tif frames'];
        if ~isdir(savedir)
            mkdir(savedir)
        end
        cd(savedir)
        savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.tif'];
        imwrite(uint16(saveimg),savename,'tif');
        
        cd(specialdir)
%         savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' figure png frames'];
        savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' figure tiff frames'];
        if ~isdir(savedir)
            mkdir(savedir)
        end
        cd(savedir)
%         savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.png'];
%         saveas(f33,savename,'png')
        savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.tif'];
        saveas(f33,savename,'tif')
        disp(num2str(ImageDetails.Frame));
        cd(olddir)
        
        cd(specialdir)
        savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel ' - ' trackfilestr ' movieCode'];
        if ~isdir(savedir)
            mkdir(savedir)
        end
        cd(savedir)
        savenamestr = ['uitrackcellz_' ImageDetails.Channel '_' ImageDetails.Scene '_' num2str(tnum) '.m'];
        mdir = pathStruct.mdirPath;
        copyfile([mdir '.m'],[savedir '/' savenamestr '.m'])
        cd(olddir)
    end
end
end

%BEST

function saveSmallMovieFunct(~,~)
global togStruct
togStruct.SaveMovie = true;
displaySmallMovieFunct([],[])
togStruct.SaveMovie = false;
end
function displaySmallMovieFunct(~,~)
global timeVec stackstruct trackfilestr dirStruct expDetailsStruct Tracked segStruct Ifstack channelimgstack nucleusimgstack DICimgstack oldscenename togStruct timeFrames tcontrast lcontrast  ImageDetails prcntl lprcntl cmap cmaplz

bkgmedian=0;

scenestr = ImageDetails.Scene;

% click the desired cell
ArrayStruct = Tracked.arrayStruct;
trackmatrix = Tracked.trackmatrix;
%choose the initial cell
t = ImageDetails.Frame;
notACell = true;
while notACell %keep running until a cell is clicked
    ginputnum=1;
    [cxx,cyy,tkeeper,button] = identifyCellbyClick(t,ginputnum);
    if isempty(button)
        return
    end
    %find cell
    [idx,~,notACell] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,tkeeper);
    if ~notACell && ~isempty(idx)
        plotTheChosenCells(idx);
    end
end
alltrackframe = trackmatrix(t,:);
cellidx = alltrackframe == idx;


%%%% make movie by iterating through frames
f33 = figure(33);
axes();
MainAxes = gca;


cropsize = 75;
imgdim = ImageDetails.ImgSize;
celltrack = trackmatrix(:,cellidx);
nannum = find(~isnan(celltrack));
segment_Centroids_array = ArrayStruct.centroid;
tval=t;
for tnum = nannum(1):nannum(end)
    t = tnum;
    
    If = Ifstack(:,:,t);
    nucleusImg = stackstruct.(segStruct.nucleus_seg)(:,:,t);
    cellImg = stackstruct.(segStruct.background_seg)(:,:,t);
    DICImg = stackstruct.('DIC')(:,:,t);
    
    if strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
        channelimg = nucleusImg;
    elseif strcmp(ImageDetails.Channel,segStruct.background_seg)
        channelimg = cellImg;
    elseif strcmp(ImageDetails.Channel,'DIC')
        channelimg = DICImg;
    elseif strcmp(ImageDetails.Channel,'overlay')
        channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
        channelimg(:,:,1) = nucleusImg;
        channelimg(:,:,2) = cellImg;
        channelimg(:,:,3) = DICImg;
    elseif strcmpi(ImageDetails.Channel,'FluorOnlyOverlay')
        channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
        channelimg(:,:,1) = nucleusImg;
        channelimg(:,:,2) = cellImg;
    end
    
    
    
    cellnum = celltrack(tnum);
    centroids = round(segment_Centroids_array{tnum});
    cellcentroid = centroids(cellnum,:);
    
    
    cropymin = max([1 cellcentroid(1)-cropsize]);
    cropymax = min([cellcentroid(1)+cropsize imgdim(1)]);
    cropxmin = max([1 cellcentroid(2)-cropsize]);
    cropxmax = min([cellcentroid(2)+cropsize imgdim(2)]);
    
    xdist = (cropxmax+1 - cropxmin);
    ydist = (cropymax+1 - cropymin);
    xl = cropsize*2 + 1;
    xw = cropsize*2 + 1;
    
    xoff = xl - xdist;
    x1 = ceil(xoff/2);
    x2 = xoff-x1;
    yoff = xw - ydist;
    y1 = ceil(yoff/2);
    y2 = yoff-y1;
    
    if length(size(channelimg))>2
        cropim = zeros(cropsize*2 + 1,cropsize*2 + 1,size(channelimg,3));

        cropim(x1+1 :(xl-x2),y1+1 :(xw-y2),:) = channelimg(cropxmin:cropxmax,cropymin:cropymax,:);
    else
        cropim = zeros(xl,xw);
        cropim(x1+1 :(xl-x2),y1+1 :(xw-y2)) = channelimg(cropxmin:cropxmax,cropymin:cropymax);
    end
    newIf = zeros(xl,xw);
    newIf(x1+1 :(xl-x2),y1+1 :(xw-y2)) = If(cropxmin:cropxmax,cropymin:cropymax);
    If = newIf;
    channelimg = cropim;
    
    
    


    
    
    
    %scripts for displaying contrasted image
    alog = regexpi(ImageDetails.Channel,'overlay'); %when overlay display is desired
    if ~isempty(alog)
        disprgb = zeros(size(cropim));
        channelimgrgb = cropim;
        Ifrgb = ind2rgb(uint8(If),vertcat([0 0 0],cmap(257:end,:)));
        
        for i = 1:size(cropim,3)
            if i==3
                channelimg = channelimgrgb(:,:,i);
                lprcntl = prctile(channelimg(:),lcontrast);
                prcntl = prctile(channelimg(:),tcontrast);
                scaleFactor = 0.6./(prcntl - lprcntl);
                dispimg = channelimg.*scaleFactor;
                dispimg = dispimg-(lprcntl.*scaleFactor);
                dispimg(dispimg> 1) =1;
                dispimg(dispimg<0) = 0;
            else
                channelimg = channelimgrgb(:,:,i);
                lprcntl = -100;
                prcntl = prctile(channelimg(:),tcontrast);
                scaleFactor = 1./(prcntl - lprcntl);
                dispimg = channelimg.*scaleFactor;
                dispimg = dispimg-(lprcntl.*scaleFactor);
                dispimg(dispimg> 1) =1;
                dispimg(dispimg<0) = 0;
            end
            
            
            if i==1
                disprgb(:,:,i) = dispimg;
            elseif i==2
                disprgb(:,:,i) = dispimg;
            elseif i==3
                for j = 1:3
                    cimg =  disprgb(:,:,j);
                    cimg(dispimg>cimg) = dispimg(dispimg>cimg);
                    disprgb(:,:,j) = cimg;
                end
            end
        end
        
        newimg = zeros(size(disprgb));
        for i = 1:size(disprgb,3)
            Ifsub = Ifrgb(:,:,i);
            Ifnew = zeros(size(Ifsub));
            Ifnew(If>0) = Ifsub(If>0);
            dispimg = disprgb(:,:,i);
            dispimg(If>0) = Ifsub(If>0);
            disprgb(:,:,i) = dispimg;
            newimg(:,:,i) = Ifnew;
        end
        
        
        dispimg = disprgb;
        %     dispimg = uint8(dispimg);
        image(MainAxes,dispimg);
        %     colormap(MainAxes,cmap);
        %     set(MainAxes,'CLim',[0 size(cmap,1)]);
        
    else  %under normal circumstances
        if togStruct.changeSceneOrChannel ==1
            prcntl = prctile(channelimg(:),tcontrast);
            togStruct.changeSceneOrChannel =0;
        end
        %     lprcntl = bkgmedian.*0.90;
        lprcntl = bkgmedian;
        scaleFactor = 255./(prcntl - lprcntl);
        dispimg = channelimg.*scaleFactor;
        dispimg = dispimg-(lprcntl.*scaleFactor);
        dispimg(dispimg > 250) =253;
        dispimg(dispimg < 0) = 0;
        %     dispimg(If>0)=255;
        
        
%         dispimg(If>0)=If(If>0)+257;


        imagesc(MainAxes,int16(dispimg));
        colormap(MainAxes,cmap);
        set(MainAxes,'CLim',[0 size(cmap,1)]);
    end
    
    %title the displayed image
    
    
    
%     frametext= text(MainAxes,0,0,[' frame ' num2str(tnum) '/' num2str(timeFrames)]);
%     % frametext.Units = 'normalized';
%     imgsize = size(dispimg);
%     textspace = round(imgsize(1)/100);
%     texty = imgsize(2)-textspace*2;
%     frametext.Position = [textspace texty];
%     frametext.FontSize = 12;
%     frametext.Color =[1 0 0];
%     
    %     frametext= text(ax1,0,0,['frame ' num2str(tnum) '/' num2str(timeFrames)]);
%     frametext.Units = 'normalized';
%     frametext.Position = [0.01 0.95];
%     frametext.FontSize = 12;
%     frametext.Color =[1 1 0];
    minstr = '0000';
    minsub = num2str(round(timeVec(tnum)./60,1,'decimals'));
    minidx = false(size(minstr));
    minidx(end-(length(minsub)-1):end) = true;
    if timeVec(tnum)<0
        minstr(end-(length(minsub)-1):end) = minsub;
        minstr(~minidx) = ' ';
    else
        minstr(end-(length(minsub)-1):end) = minsub;

        minstr(~minidx)=' ';
    end
    
    frametext= text(MainAxes,0,0,[minstr ' hr']);
    frametext.Units = 'normalized';
    frametext.Position = [0.99 0.035];
    frametext.FontSize = 15;
    frametext.Color =[1 1 0];
    frametext.FontWeight = 'bold';
    frametext.HorizontalAlignment = 'right';
    
    if ~(timeVec(tnum)<0)
        frametext= text(MainAxes,0,0,'+TGFbeta');
        frametext.Units = 'normalized';
        frametext.Position = [0.99 0.10];
        frametext.FontSize = 15;
        frametext.Color =[1 1 0];
        frametext.FontWeight = 'bold';
        frametext.HorizontalAlignment = 'right';
    end




set(MainAxes,'NextPlot','Add')
%plot cell tracking "tails"
if ~(tnum==1)
    if togStruct.displayTrackingToggle==1
        trajectForPlot = trackingTrajectories(timeFrames);
        mainplotX = squeeze(trajectForPlot(:,1,:)); %28x22 means 28 cells on frame 22;
        mainplotY = squeeze(trajectForPlot(:,2,:));
        
        if size(mainplotY,2) == 1
            mainplotY=mainplotY';
            mainplotX=mainplotX';
        end
        %only plot if the cell is currently tracked/segmented in this frame
        idx = ~isnan(mainplotY(cellidx,tnum));
        idxa = find(cellidx);

            cmapl = cmaplz;
            plotcmap = zeros(length(idxa),3);
            for i=1:length(idxa)
                plotcmap(i,:) =  cmapl(idxa(i),:);
            end
            
            x = mainplotX(cellidx,:) - cellcentroid(1) + cropsize;
            y = mainplotY(cellidx,:) - cellcentroid(2) + cropsize;
            t0 =max([1 tnum-20]);
            %             h = plot(MainAxes,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',2);
            h = line(MainAxes,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',4);
            set(h, {'color'}, num2cell(plotcmap,2));

    end
end
drawnow
set(MainAxes,'NextPlot','Replace')
if ~togStruct.SaveMovie
% pause(1/15)
end

cellidxstr = num2str(find(cellidx));

    if togStruct.SaveMovie
        specialdir = [dirStruct.parentdir 'LookingAtData/SingleCellMovies'];
        olddir = pwd;
        saveimg = channelimg;
        if ~isdir(specialdir)
            mkdir(specialdir)
        end
        cd(specialdir)
        savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel 'cell' cellidxstr ' - ' trackfilestr ' tif'];
        if ~isdir(savedir)
            mkdir(savedir)
        end
        cd(savedir)
        savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.tif'];
        imwrite(uint16(saveimg),savename,'tif');
        
        cd(specialdir)
        savedir = [specialdir '/' expDetailsStruct.expDateStr '/' ImageDetails.Scene '-' ImageDetails.Channel 'cell' cellidxstr ' - ' trackfilestr ' png'];
        if ~isdir(savedir)
            mkdir(savedir)
        end
        cd(savedir)
        savename = [ImageDetails.Channel ' ' ImageDetails.Scene ' ' num2str(tnum) '.tif'];
%         saveas(f33,savename,'tif')
        disp(num2str(ImageDetails.Frame));
        cd(olddir)
    end





end
end

function cmaplz = colorcubemodified(cmaplength,cmapstr)
%default darkmin = 0.2
%default brightmax = 2
darkmin = 0.3;
colorationmin = 0.4;
brightmax = 2.4;

%generate colormap based on number of cells tracked
cnew=[];
%                 ccc = vertcat(colormap('summer'),colormap('autumn'),colormap('winter'),colormap('spring'));
%                 ccc = vertcat(colormap('hsv'),colormap('hot'));
%     cmapccc = colorcube();

numunique = max([floor(cmaplength^(1/3)) 3]);
a = linspace(0,1,numunique);
cycle=0;
la = length(a);
ccc = zeros(la*la*la,3);
for i=1:la
    for j=1:la
        for k=1:la
            cycle=cycle+1;
            ccc(cycle,:) = [a(i) a(j) a(k)];
        end
    end
end
cmapccc = ccc;
%     cmapccc = zeros(1000,3);
%     for j = 1:size(cmapccc,2)
%         x = linspace(0,1,size(ccc,1));
%         v = ccc(:,j);
%         xq = linspace(0,1,1000);
%         cmapccc(1:length(xq),j) = interp1(x,v,xq);
%     end

%                 ccc = colormap(jet(1000));
cccyc = 0;
for k = 1:size(cmapccc,1)
    cvec = cmapccc(k,:);
    colortest = abs([cvec(1)-cvec(2) cvec(2)-cvec(3) cvec(3)-cvec(1)]);
    if sum(cvec)>darkmin && max(colortest)>colorationmin && sum(cvec)<brightmax
        cccyc = cccyc+1;
        cnew(cccyc,:) = cvec;
    end
end

%     [~,cnewidx] = sort(cnew(:,1));
%     cnew = cnew(cnewidx,:);
ccnew = zeros(cmaplength,3);
for j = 1:size(ccnew,2)
    x = linspace(0,1,size(cnew,1));
    v = cnew(:,j);
    xq = linspace(0,1,cmaplength);
    ccnew(1:length(xq),j) = interp1(x,v,xq);
end

% randidx = randi(size(ccnew,1),1,size(ccnew,1));
randidx = randperm(size(ccnew,1));
ccnew = ccnew(randidx,:);
cmaplz = ccnew;
end

%% functions for determining variables

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

%% ginput modified function
function [out1,out2,out3] = ginputmod(arg1)
%GINPUT Graphical input from mouse.
%   [X,Y] = GINPUT(N) gets N points from the current axes and returns
%   the X- and Y-coordinates in length N vectors X and Y.  The cursor
%   can be positioned using a mouse.  Data points are entered by pressing
%   a mouse button or any key on the keyboard except carriage return,
%   which terminates the input before N points are entered.
%
%   [X,Y] = GINPUT gathers an unlimited number of points until the
%   return key is pressed.
%
%   [X,Y,BUTTON] = GINPUT(N) returns a third result, BUTTON, that
%   contains a vector of integers specifying which mouse button was
%   used (1,2,3 from left) or ASCII numbers if a key on the keyboard
%   was used.
%
%   Examples:
%       [x,y] = ginput;
%
%       [x,y] = ginput(5);
%
%       [x, y, button] = ginput(1);
%
%   See also GTEXT, WAITFORBUTTONPRESS.

%   Copyright 1984-2015 The MathWorks, Inc.

out1 = []; out2 = []; out3 = []; y = [];

if ~matlab.ui.internal.isFigureShowEnabled
    error(message('MATLAB:hg:NoDisplayNoFigureSupport', 'ginput'))
end

% Check Inputs
if nargin == 0
    how_many = -1;
    b = [];
else
    how_many = arg1;
    b = [];
    if  ~isPositiveScalarIntegerNumber(how_many)
        error(message('MATLAB:ginput:NeedPositiveInt'))
    end
    if how_many == 0
        % If input argument is equal to zero points,
        % give a warning and return empty for the outputs.
        warning (message('MATLAB:ginput:InputArgumentZero'));
    end
end

% Get figure
fig = gcf;
figure(gcf);

% Make sure the figure has an axes
gca(fig);

% Setup the figure to disable interactive modes and activate pointers.
initialState = setupFcn(fig);

% onCleanup object to restore everything to original state in event of
% completion, closing of figure errors or ctrl+c.
c = onCleanup(@() restoreFcn(initialState));

drawnow
char = 0;

while how_many ~= 0
    waserr = 0;
    try
        keydown = wfbp;
    catch %#ok<CTCH>
        waserr = 1;
    end
    if(waserr == 1)
        if(ishghandle(fig))
            cleanup(c);
            error(message('MATLAB:ginput:Interrupted'));
        else
            cleanup(c);
            error(message('MATLAB:ginput:FigureDeletionPause'));
        end
    end
    % g467403 - ginput failed to discern clicks/keypresses on the figure it was
    % registered to operate on and any other open figures whose handle
    % visibility were set to off
    figchildren = allchild(0);
    if ~isempty(figchildren)
        ptr_fig = figchildren(1);
    else
        error(message('MATLAB:ginput:FigureUnavailable'));
    end
    %         old code -> ptr_fig = get(0,'CurrentFigure'); Fails when the
    %         clicked figure has handlevisibility set to callback
    if(ptr_fig == fig)
        if keydown
            char = get(fig, 'CurrentCharacter');
            button = abs(get(fig, 'CurrentCharacter'));
        else
            button = get(fig, 'SelectionType');
            if strcmp(button,'open')
                button = 1;
            elseif strcmp(button,'normal')
                button = 1;
            elseif strcmp(button,'extend')
                button = 2;
            elseif strcmp(button,'alt')
                button = 3;
            else
                error(message('MATLAB:ginput:InvalidSelection'))
            end
        end
        axes_handle = gca;
        drawnow;
        pt = get(axes_handle, 'CurrentPoint');
        
        how_many = how_many - 1;
        
        if(char == 13) % & how_many ~= 0)
            % if the return key was pressed, char will == 13,
            % and that's our signal to break out of here whether
            % or not we have collected all the requested data
            % points.
            % If this was an early breakout, don't include
            % the <Return> key info in the return arrays.
            % We will no longer count it if it's the last input.
            break;
        end
        
        out1 = [out1;pt(1,1)]; %#ok<AGROW>
        y = [y;pt(1,2)]; %#ok<AGROW>
        b = [b;button]; %#ok<AGROW>
    end
end

% Cleanup and Restore
cleanup(c);

if nargout > 1
    out2 = y;
    if nargout > 2
        out3 = b;
    end
else
    out1 = [out1 y];
end

end

function valid = isPositiveScalarIntegerNumber(how_many)
valid = ~ischar(how_many) && ...            % is numeric
    isscalar(how_many) && ...           % is scalar
    (fix(how_many) == how_many) && ...  % is integer in value
    how_many >= 0;                      % is positive
end

function key = wfbp
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

fig = gcf;
current_char = []; %#ok<NASGU>

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
try
    h=findall(fig,'Type','uimenu','Accelerator','C');   % Disabling ^C for edit menu so the only ^C is for
    set(h,'Accelerator','');                            % interrupting the function.
    keydown = waitforbuttonpress;
    current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
    if~isempty(current_char) && (keydown == 1)          % If the character was generated by the
        if(current_char == 3)                           % current keypress AND is ^C, set 'waserr'to 1
            waserr = 1;                                 % so that it errors out.
        end
    end
    
    set(h,'Accelerator','C');                           % Set back the accelerator for edit menu.
catch %#ok<CTCH>
    waserr = 1;
end
drawnow;
if(waserr == 1)
    set(h,'Accelerator','C');                          % Set back the accelerator if it errored out.
    error(message('MATLAB:ginput:Interrupted'));
end

if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function initialState = setupFcn(fig)

% Store Figure Handle.
initialState.figureHandle = fig;

% Suspend figure functions
initialState.uisuspendState = uisuspend(fig);

% Disable Plottools Buttons
initialState.toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
if ~isempty(initialState.toolbar)
    initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
        uigettool(initialState.toolbar,'Plottools.PlottoolsOn')];
    initialState.ptState = get (initialState.ptButtons,'Enable');
    set (initialState.ptButtons,'Enable','off');
end

%Setup empty pointer
cdata = NaN(16,16);
% cdata = ones(16,16)*2;
hotspot = [8,8];
set(gcf,'Pointer','custom','PointerShapeCData',cdata,'PointerShapeHotSpot',hotspot)



% Adding this to enable automatic updating of currentpoint on the figure
% This function is also used to update the display of the fullcrosshair
% pointer and make them track the currentpoint.
set(fig,'WindowButtonMotionFcn',@(o,e) dummy()); % Add dummy so that the CurrentPoint is constantly updated
% Create uicontrols to simulate fullcrosshair pointer.
initialState.CrossHair = createCrossHair(fig);
% updateCrossHair(fig,initialState.CrossHair)
initialState.MouseListener = addlistener(fig,'WindowMouseMotion', @(o,e) updateCrossHair(o,initialState.CrossHair));

% Get the initial Figure Units
initialState.fig_units = get(fig,'Units');
end

function restoreFcn(initialState)
if ishghandle(initialState.figureHandle)
    delete(initialState.CrossHair);
    
    % Figure Units
    set(initialState.figureHandle,'Units',initialState.fig_units);
    
    set(initialState.figureHandle,'WindowButtonMotionFcn','');
    delete(initialState.MouseListener);
    
    % Plottools Icons
    if ~isempty(initialState.toolbar) && ~isempty(initialState.ptButtons)
        set (initialState.ptButtons(1),'Enable',initialState.ptState{1});
        set (initialState.ptButtons(2),'Enable',initialState.ptState{2});
    end
    
    % UISUSPEND
    uirestore(initialState.uisuspendState);
end
end

function updateCrossHair(fig, crossHair)
% update cross hair for figure.
gap = 5; % 3 pixel view port between the crosshairs
cp = hgconvertunits(fig, [fig.CurrentPoint 0 0], fig.Units, 'pixels', fig);
cp = cp(1:2);
figPos = hgconvertunits(fig, fig.Position, fig.Units, 'pixels', fig.Parent);
figWidth = figPos(3);
figHeight = figPos(4);

% Early return if point is outside the figure
if cp(1) < gap || cp(2) < gap || cp(1)>figWidth-gap || cp(2)>figHeight-gap
    return
end

set(crossHair, 'Visible', 'on');
thickness = 2; % 1 Pixel thin lines.
set(crossHair(1), 'Position', [0 cp(2) cp(1)-gap thickness]);
set(crossHair(2), 'Position', [cp(1)+gap cp(2) figWidth-cp(1)-gap thickness]);
set(crossHair(3), 'Position', [cp(1) 0 thickness cp(2)-gap]);
set(crossHair(4), 'Position', [cp(1) cp(2)+gap thickness figHeight-cp(2)-gap]);

offset = 50;
boxpoint = 4;
set(crossHair(1), 'Position', [cp(1)-offset cp(2) offset-gap thickness]);
set(crossHair(2), 'Position', [cp(1)+gap cp(2) offset-gap thickness]);
set(crossHair(3), 'Position', [cp(1) cp(2)-offset thickness offset-gap]);
set(crossHair(4), 'Position', [cp(1) cp(2)+gap thickness offset-gap]);
set(crossHair(5), 'Position', [cp(1)-1 cp(2)-1 boxpoint boxpoint]);
end

function crossHair = createCrossHair(fig)
% Create thin uicontrols with black backgrounds to simulate fullcrosshair pointer.
% 1: horizontal left, 2: horizontal right, 3: vertical bottom, 4: vertical top
cmapz=jet(4);
for k = 1:5
    crossHair(k) = uicontrol(fig, 'Style', 'text', 'Visible', 'off', 'Units', 'pixels', 'BackgroundColor', [0.2 0.9 0], 'HandleVisibility', 'off', 'HitTest', 'off'); %#ok<AGROW>
end
end

function cleanup(c)
if isvalid(c)
    delete(c);
end
end

function dummy(~,~)
end


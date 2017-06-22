function uiTrackCellz
global xAxisLimits DICimgstack dfoName cfoName trackingPath background_seg bfoName nfoName sfoName cell_seg nucleus_seg segmentimgstack channelimgstack segmentPath mstackPath runIterate ExportNameKey ExportName exportdir plottingTotalOrMedian channelinputs adjuster cmapper tcontrast lcontrast ThirdPlotAxes SecondPlotAxes OGExpDate plottingON PlotAxes cmap TC expDirPath  timeFrames frameToLoad ImageDetails MainAxes SceneList displaytracking imgsize ExpDate


%determine matfile directory
    mdir = mfilename('fullpath');
        [~,b ] = regexp(mdir,'/');
            if isempty(b)
                [~,b] = regexp(mdir,'\');
            end
    parentdir = mdir(1:b(end)); %folder in which matfile exists
    exportdir = strcat(parentdir,'Export/');

    [~,b ] = regexp(parentdir,'/');
        if isempty(b)
            [~,b] = regexp(parentdir,'\');
        end
    gparentdir = parentdir(1:b(end-1)); %folder in which matfile parentdir exists
    
    
%specify important directory names
    mstackName = 'flat mstack';
    trackName = 'tracking files';
    segmentName = 'segment mstack';
    

%set contrast and plotting details
    adjuster=0;
    plottingTotalOrMedian = 'median';
    tcontrast = 99;
    lcontrast = 1;
    runIterate =0;
    TC = 1;
    ImageDetails = InitializeImageDetails;
    displaytracking = 0;
    plottingON =0;
    xAxisLimits = [0 50];
    
    
%determine export details
    ExportNameKey = 'final';
    if strcmp(ExportNameKey,'final')
    else
    disp(strcat('Export name key is "',ExportNameKey,'" not FINAL'))
    end
    ExportName = 'fricktrack';
   
    
%initialize global variables    
    cfoName = [];
    channelimgstack =[];
    segmentimgstack =[];
    sfoName =[];
    bfoName = [];
    nfoName = [];
    DICimgstack=[];
    dfoName=[];



%set colormap
    cd(parentdir)
    addpath('Colormaps')
        cmap = colormap(gray(255));
%         cmap = colormap(viridis(255));
        % cmap = colormap(magma(255));
        % cmap = colormap(inferno(255));
        % cmap = colormap(plasma(255));
        cmap(255,:)=[1 0 0];
        cmapper = cmap;
        close all

        
% user selects experiment directory to be analyzed
    cd(gparentdir)
    expDirPath = uigetdir;

    cd(expDirPath)
    experimentdir = expDirPath;
    mstackPath = strcat(experimentdir,'/',mstackName);
    segmentPath = strcat(experimentdir,'/',segmentName);
    trackingPath = strcat(experimentdir,'/',trackName);

%make tracking file folder if it does not exist
    cd(experimentdir)
    dirlist = dir(trackName);
    if isempty(dirlist)
        mkdir(trackName);
    end

%determine date of experiment
    [a,b] = regexp(expDirPath,'201[0-9]');
    [~,d] = regexp(expDirPath,'exp[0-9]');
    ExpDate = expDirPath(a:b+6);OGExpDate = expDirPath(a:d); [a,~] = regexp(ExpDate,'_');ExpDate(a) = '-';


%subdirectories should include
%> [ flat mstack ]
%> [ mstack images ]
%> [ segment mstack ]
%> [ tracking files ]


%load helpful metadata
    cd(exportdir)
    FileName = OGExpDate;
    datequery = strcat(FileName,'*DoseAndScene*');
    cd(exportdir)
    filelist = dir(datequery);
        if isempty(filelist)
            error(strcat('need to run ExtractMetadata for-',FileName));
    %        dosestruct = makeDoseStruct; %run function to make doseStruct 
        else
            dosestructstruct = load(char(filelist.name));
            dosestruct = dosestructstruct.dosestruct;
        end
        segInstruct = dosestructstruct.segInstruct;

    nucleus_seg = segInstruct.nucleus;
    cell_seg = segInstruct.cell;
    background_seg = segInstruct.background;
    channelstoinput = dosestructstruct.channelNames;
    channelinputs =channelregexpmaker(channelstoinput);
    bkg = dosestructstruct.BACKGROUND;
    imgsize = dosestructstruct.dimensions;

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

timeVec = dosestructstruct.timeVec;
% timeMatrix(i,:) = timeVec(idxtwo,1:finalFrame)-timeVec(1,stimulationFrame+1); %subtract by time closest to Tgfbeta addition 


%determine how many scenes are present
    dirlist = dir(mstackPath);
    [~,~,~,d] = regexp({dirlist.name},'s[0-9]+');
    dlog = ~cellfun(@isempty,d,'UniformOutput',1); 
    dcell = d(dlog);
    SceneList = unique(cellfun(@(x) x{1},dcell,'UniformOutput',0));
    
%remove scenes that are background images
    [~,~,~,d] = regexp(SceneList,bkinputs);
    bkgscenelog = cellfun(@isempty,d,'UniformOutput',1);
    SceneList = SceneList(bkgscenelog);


%determine the number of frames per scene
    cd(mstackPath)
    dirlist = dir('*.mat');
    filearray = {dirlist.name};
    filename = filearray{1};
    fileObject = matfile(filename);
    dim = size(fileObject,'flatstack');
    timeFrames = dim(3);
    frameToLoad = 1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Set up  user interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;
f.Units = 'pixels';
fpos = [1 1 1800 900];
f.Position =[1 1 1800 900];
fW = fpos(3);
fH = fpos(4);


bW = 90;
bH = 25;
yP = sort(100:bH:(fH-fH./9),'descend');
xP = ones(1,length(yP)).*(fW-fW./9);
fontSize = 10;
%instruct user how they can change to view different channels
    mmm=1;
    chooseChannelText = uicontrol('Style','text','String','To choose channel push 1, 2, or 3',...
              'Position',[xP(mmm),yP(mmm)+bH./1.5,bW,bH*2]);


%set axis limits button

    hplotAxisLimits = uicontrol('Style','pushbutton',...
        'String','plotAxisLimits [x]',...
        'Position',[fW-(bW*6),yP(mmm),bW,bH],...
        'Callback',@plotAxis_callback);

%move forward frame or previuos frame
    mmm=2;
    hNextFrame = uicontrol('Style','pushbutton',...
        'String','NextFrame [f]',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@nextbutton_callback);
    hPreviousFrame = uicontrol('Style','pushbutton',...
        'String','PreviousFrame [a]',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@prevbutton_callback);
    hGoToFrame = uicontrol('Style','pushbutton',...
        'String','GoToFrame',...
        'Position',[xP(mmm)-(bW./2)-bW,yP(mmm),bW,bH],...
        'Callback',@gotobutton_callback);

%next row is final Frame and first Frame commands   
    mmm=3;
    hFinalFrame = uicontrol('Style','pushbutton',...
        'String','FinalFrame [g]',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@finalbutton_callback);
    hFirstFrame = uicontrol('Style','pushbutton',...
        'String','FirstFrame [z]',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@firstbutton_callback);

%choose scene text and popup menu    
    mmm=5;
    chooseSceneText = uicontrol('Style','text','String','Choose Scene',...
            'Position',[xP(mmm),yP(mmm),bW,bH]);

    hpopup = uicontrol('Style','popupmenu',...
        'String',SceneList',...
        'Position',[xP(mmm),yP(mmm)-bH./2,bW,bH./1.5],...
        'Callback',@popup_menu_Callback);

%add area and remove area
    mmm=7;
    hAddArea = uicontrol('Style','pushbutton','String','AddArea [v]',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@addareabutton_Callback);
    hRemoveArea = uicontrol('Style', 'pushbutton', 'String', 'Remove area',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@removeArea_Callback);
       
%delete commands
    mmm=8; 
    hDelete = uicontrol('Style','pushbutton','String','Delete',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@deletebutton_Callback);
    hDeleteAllOnFrame = uicontrol('Style','pushbutton','String','DeleteAllOnFrame',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@deleteAllonFrame_Callback);


%destroy commands    
    mmm=9;
    hDestroy = uicontrol('Style','pushbutton','String','Destroy [d]',...
        'Position',[xP(mmm),yP(mmm),bW,bH],...
        'Callback',@destroybutton_Callback);
    hDestroy = uicontrol('Style','pushbutton','String','DestroyPrevious',...
        'Position',[xP(mmm)-bW,yP(mmm),bW,bH],...
        'Callback',@destroybuttonAllPrevious_Callback);
        hDestroy.FontSize=fontSize-2;
    hDestroy = uicontrol('Style','pushbutton','String','DestroySubsequent',...
        'Position',[xP(mmm)+bW,yP(mmm),bW,bH],...
        'Callback',@destroybuttonAllSubsequent_Callback);
        hDestroy.FontSize=fontSize-2;
        
%destroy commands    
    mmm=10;

    hDestroy = uicontrol('Style','pushbutton','String','DestroyAllFramePrevious',...
        'Position',[xP(mmm)-bW,yP(mmm),bW+bW./2,bH],...
        'Callback',@destroyAllFramePrevious_Callback);
        hDestroy.FontSize=fontSize-2;
    hDestroy = uicontrol('Style','pushbutton','String','DestroyAllFrameSubsequent',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW+bW./2,bH],...
        'Callback',@destroyAllFrameSubsequent_Callback);
        hDestroy.FontSize=fontSize-2;


%chosen ones commands    
    mmm=11;
    hchosenOnes = uicontrol('Style','pushbutton','String','Chosen Ones',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@chosenOnes_Callback);
    hchosenOnes = uicontrol('Style','pushbutton','String','Chosen OnesAllOnFrame',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@chosenOnesAllOnFrame_Callback);



  



%erode or dilate nuclei       
    mmm=12;
    hErode = uicontrol('Style','pushbutton','String','Erode Selected Nuclei',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@erodeOnes_Callback);
    hDilate = uicontrol('Style','pushbutton','String','Dilate Selected Nuclei',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@dilateOnes_Callback);


%tracking commands
mmm=14;
hLinkCells = uicontrol('Style','pushbutton','String','LinkCells [r]',...
    'Position',[xP(mmm),yP(mmm),bW,bH],...
    'Callback',@linkCells_Callback);
mmm=15;
hTrack = uicontrol('Style','pushbutton',...
    'String','Run Tracking [t]',...
    'Position',[xP(mmm)-bW./2,yP(mmm),bW*2,bH],...
    'Callback',@trackbutton_Callback);    
mmm=16;
hLoadTracking = uicontrol('Style','pushbutton',...
    'String','LoadTracking',...
    'Position',[xP(mmm),yP(mmm),bW,bH],...
    'Callback',@loadTrackingFile_callback);


        

%display commands
    mmm=18;
    hContrast = uicontrol('Style','pushbutton',...
        'String','contrast user',...
        'Position',[xP(mmm),yP(mmm),bW,bH],...
        'Callback',@contrast_Callback);
    hDisplayTracking = uicontrol('Style','pushbutton',...
        'String','DisplayTracking [m]',...
        'Position',[xP(mmm),yP(mmm)-bH./2,bW,bH./1.5],...
        'Callback',@displayTrackingButton_Callback);
        

%save commands
    mmm=20;
    hSaveTrackingAs = uicontrol('Style','pushbutton',...
        'String','SaveTrackingAs',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW+bW,bH],...
        'Callback',@saveTrackingFileAs_callback);
    

%autotrack commands
    mmm=21;
    hTrackSaveIterate = uicontrol('Style','pushbutton',...
            'String','trackSaveIterate',...
            'Position',[xP(mmm)+bW./2,yP(mmm),bW+bW./2,bH./1.5],...
            'Callback',@trackSaveIterate_callback);
    hTrackSaveIterate = uicontrol('Style','pushbutton',...
        'String','TSIchosen',...
        'Position',[xP(mmm)-bW,yP(mmm),bW+bW./2,bH./1.5],...
        'Callback',@trackSaveIterateChosen_callback);


%plot and plot settings commands
    mmm=23;  
    hPlot = uicontrol('Style','pushbutton',...
        'String','PLOT!',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@Plot_callback);
    hPlotSpecificCell = uicontrol('Style','pushbutton',...
        'String','Plot Specific Cell!',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@Plot_SpecificCell_callback);
    mmm=24;
    hPlotSettings = uicontrol('Style','pushbutton',...
        'String','Plot Settings!',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW+bW./2,bH./1.5],...
        'Callback',@PlotSettings_callback);

%export commands  
    mmm=26;
    hExportCells = uicontrol('Style','pushbutton',...
        'String','ExportTrackedCells',...
        'Position',[xP(mmm),yP(mmm),bW,bH*2],...
        'Callback',@exportTrackedCells);
    hExportCells = uicontrol('Style','pushbutton',...
        'String','ExportAllCells',...
        'Position',[xP(mmm)-bW,yP(mmm),bW,bH./1.5],...
        'Callback',@exportAllCells);
    hExportCells = uicontrol('Style','pushbutton',...
        'String',{'ExportSegmentedCells'},...
        'Position',[xP(mmm)+bW,yP(mmm),bW,bH*2],...
        'Callback',@exportSegmentedCells);

%label and comment commands
    mmm=28;
    hLabelCells = uicontrol('Style','pushbutton',...
        'String','Label Cells',...
        'Position',[xP(mmm)-bW,yP(mmm),bW,bH./1.5],...
        'Callback',@labelCells); 
    hcomment = uicontrol('Style','pushbutton','String','Comments',...
        'Position',[xP(mmm),yP(mmm),bW,bH./1.5],...
        'Callback',@comment_Callback);
    hExportLabelsCells = uicontrol('Style','pushbutton',...
        'String','ExportLabels',...
        'Position',[xP(mmm)+bW,yP(mmm),bW,bH./1.5],...
        'Callback',@ exportLabels);

        

f.Visible = 'on'   ;
f.Units = 'normalized';
for i = 1:length(f.Children)
   hhh = f.Children(i);
   hhh.Units = 'normalized';
end

MainAxes = axes;
MainAxes.Units = 'pixels';
MainAxes.XTick=[];
MainAxes.YTick = [];
imgdim = 512.*1.5;
Position = [25 25 imgdim imgdim];
% Position = [0.1 0.3 0.65 0.65];
MainAxes.Position = Position;
MainAxes.Units = 'normalized';


PlotAxes = axes;
% Position = [0.1 0.05 0.15 0.15];
Position = [0.6440    0.6605    0.1500    0.1500];
% Position = [822.7440 465.9920 191.4000 105.6000]
PlotAxes.Position = Position;

SecondPlotAxes = axes;
% Position = [0.3 0.05 0.15 0.15];
Position = [0.6440    0.4605    0.1500    0.1500];
SecondPlotAxes.Position = Position;

ThirdPlotAxes = axes;
% Position = [0.5 0.05 0.15 0.15];
Position = [0.6440    0.2605    0.1500    0.1500];
ThirdPlotAxes.Position = Position;

% f.Position =[0.1,0.1,0.7,0.8];
% f.Position = [0.1461 0.1370 0.4457 0.7315];
f.Units = 'pixels';
f.Position = fpos;
f.Color = 'w';
set(f,'KeyPressFcn',@keypress);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%% uifunctions

function PMthreshslider(source,callbackdata)
    global threshinput Imagez nucleus_seg
    str = source.String;
     threshinput.(str) =source.Value;
     
    if strcmp(source.String,'channel1')
    elseif strcmp(source.String,'channel2')
        Imagez.NucSeg = segmentNucleus(Imagez.NucImage); %semgnet the image
    elseif strcmp(source.String,'channel3')
        Imagez.CellSeg = segmentCell(Imagez.CellImage);
    end
updatePMseg
updateImage
end



function keypress(fig_obj,~)
global  ImageDetails displaycomments nucleus_seg cell_seg
key = get(fig_obj,'CurrentKey');

switch key
    case '1'
        ImageDetails.Channel = cell_seg;
        setSceneAndTime
    case '2'
        ImageDetails.Channel = nucleus_seg;
        setSceneAndTime    
    case '3'
        ImageDetails.Channel = 'DIC';
        setSceneAndTime
    case '4'
        ImageDetails.Channel = 'EGFP';
        setSceneAndTime
    case '5'
        ImageDetails.Channel = 'BKGbinary';
        setSceneAndTime
    case '6'
        ImageDetails.Channel = 'overlay';
        setSceneAndTime
    case 'q'
        prevscenebutton_Callback([],[])
    case 'w'
        nextscenebutton_Callback([],[])
    case 'a'
        prevbutton_callback([],[])
    case 'f'
        nextbutton_callback([],[])
    case 'd'
        destroybutton_Callback([],[]);
    case 't'
        trackbutton_Callback([],[]);
    case 'e'
        eliminatebutton_Callback([],[]);
    case 'v'
        addareabutton_Callback([],[]);
    case 'r'
        linkCells_Callback([],[]);  
    case 'm'
        displayTrackingButton_Callback([],[])
    case 'g'
        finalbutton_callback([],[])
    case 'z'
        firstbutton_callback([],[])
    case 's'
        saveTrackingFileAs_callback([],[])
    case 'l'
        loadTrackingFile_callback([],[])
    case 'p'
        Plot_callback([],[])
    case 'o'
        labelCells;
    case 'u'
%         if displaycomments==1
%             displaycomments=0;
%         else
            displaycomments=1;
            xy = getxy([],[]);
%         end
    case 'c'
        contrast_Callback([],[])
    case 'k'
        comment_Callback([],[])
    case 'j'
        comment_CallbackJ([],[])
    case 'n'
        PlotCFPnorm_callback([],[])
    case 'b'
        PlotCFPnotnorm_callback([],[])
    case 'x'
        plotAxis_callback([],[])
    case '0'
        displaycomments=1;
    xy = getxy([],[]);
    [~,comments,commentpos,cellidx]=updatecomments(xy);
    setcommentsTracking(comments,commentpos)
    dispxy(xy)
end

end

function plotAxis_callback(~,~)
global xAxisLimits 


prompt = {'xmin','xmax'};
dlg_title = 'set x axis limits...';
inputdlgOutput = inputdlg(prompt,dlg_title);
xAxisLimits = cellfun(@str2num,inputdlgOutput,'UniformOutput',1);
Plot_callback([],[])
end



%choose frames
function nextbutton_callback(~,~)
global frameToLoad ImageDetails timeFrames


if isempty(ImageDetails.Frame)
    ImageDetails.Frame = frameToLoad;
end

frameToLoad = ImageDetails.Frame + 1;

if frameToLoad>timeFrames
    frameToLoad = timeFrames; 
end

ImageDetails.Frame = frameToLoad;
% disp(frameToLoad)
setSceneAndTime
end
function prevbutton_callback(~,~) 
global frameToLoad ImageDetails timeFrames

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = frameToLoad;
end

frameToLoad = ImageDetails.Frame - 1;

if frameToLoad<1
    frameToLoad = 1;
end
% disp(frameToLoad)
ImageDetails.Frame = frameToLoad;
setSceneAndTime
end

function finalbutton_callback(~,~)
global frameToLoad ImageDetails timeFrames

frameToLoad = timeFrames;
ImageDetails.Frame = frameToLoad;
setSceneAndTime
end
function firstbutton_callback(~,~)
global frameToLoad ImageDetails 

frameToLoad = 1;
ImageDetails.Frame = frameToLoad;
setSceneAndTime
end
function gotobutton_callback(~,~)
global frameToLoad ImageDetails 

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = frameToLoad;
end

prompt = {'Go to which frame'};
dlg_title = 'Go to frame...';
idx = str2num(cell2mat(inputdlg(prompt,dlg_title)));

frameToLoad = idx;
ImageDetails.Frame = frameToLoad;
setSceneAndTime
end


%choose scenes
function nextscenebutton_Callback(~,~) 
global   ImageDetails SceneList expDirPath trackPath




if isempty(ImageDetails.Scene)
    ImageDetails.Scene = SceneList{1};
end

Idx = strcmp(ImageDetails.Scene,SceneList);
idx = find(Idx == 1);
if idx == length(SceneList)
else
idx = idx + 1;
end
ImageDetails.Scene = SceneList{idx};


loadTrackingFile_callback([],[])
setSceneAndTime
end
function prevscenebutton_Callback(~,~) 
global   ImageDetails SceneList expDirPath trackPath

if isempty(ImageDetails.Scene)
    ImageDetails.Scene = SceneList{1};
end

Idx = strcmp(ImageDetails.Scene,SceneList);
idx = find(Idx == 1);
if idx ==1
else
idx = idx - 1;
end
ImageDetails.Scene = SceneList{idx};


loadTrackingFile_callback([],[])
setSceneAndTime
end
function popup_menu_Callback(source,~) 
global ImageDetails Tracked timeFrames

Trackedz = makeTrackingFile(timeFrames);
Tracked=Trackedz;

% Determine the selected data set.
 str = source.String;
 val = source.Value;
 pvalue = char(str{val});

ImageDetails.Scene = pvalue;
setSceneAndTime

end

%removeArea
function removeArea_Callback(~,~)
global  ImageDetails frameToLoad Tracked imgsize
 % choose cell
%       [cellx,celly] = ginput(1);
       % construct a polygon to add

       If = 1;
       button=1;
       while button==1
      [polyx,polyy,button] = ginput();
      button = round(mean(button));
      
      if button ==1
          M = zeros(1,length(polyx)*2);
          M(1:2:end) = polyx;
          M(2:2:end) = polyy;
          zeroImage = zeros(imgsize);
          zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
          zerogray = rgb2gray(zeroImage);

          if isempty(Tracked{1}.Cellz)


          else  %if there exists segmenttracking already...then load that. 

          imagio = zeros(imgsize);
          imagio(zerogray>0)=1;
          cc = bwconncomp(imagio);
          px = cc.PixelIdxList;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   determine the frame to load
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t = ImageDetails.Frame;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    CC = Tracked{t}.Cellz;
    PX = CC.PixelIdxList;        

    idxs = cellfun(@(x) sum(ismember(x,px{1})),PX,'UniformOutput',1);
    index = find(idxs>1);
        if ~isempty(index)
            for abc = index
            oldMass = PX{abc};
            overlap = ismember(oldMass,px{1});
            oldMass(overlap)=[];
            imagio = zeros(imgsize);
            imagio(oldMass)=1;
            cc = bwconncomp(imagio);
            numcells = cc.NumObjects;
                if numcells>1
                    splitcells = cc.PixelIdxList;
                    PX(abc) = {NaN};
                        for nums = 1:numcells
                            PX{end+1} = splitcells{nums};
                        end
                else
                    PX{abc} = oldMass;
                end
            end
%             PX{min(index)} = unique(vertcat(oldMass,px{1}));
%                 if length(index)>1
%                 index(find(index == min(index)))=[];
%                 PX(index) = {NaN};
%                 end
            CC.PixelIdxList = PX;
%             end
        end
    CC.NumObjects = length(CC.PixelIdxList);
        S = regionprops(CC,'Centroid');
        Smat = vertcat(S.Centroid);
        CC.Centroid = Smat;
    Tracked{t}.Cellz = CC;
          end



        nextbutton_callback([],[]);
      end 
      end

end

%delete cells
function deletebutton_Callback(~,~) 
    global imgsize ImageDetails frameToLoad Tracked
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
     [cellxx,cellyy] = ginput();
       
% for j = 1:length(cellxx)
%     cellx = cellxx(j);
%     celly = cellyy(j);
%     
    cellx = round(cellxx);
    celly = round(cellyy);
    
    cellind = sub2ind(imgsize,celly,cellx);
    
    for j = 1:length(cellxx)
    if j==1
        idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
    else
        idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
    end
    end
    
PX(~idxs) = {NaN};
CC.PixelIdxList = PX;
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Tracked{t}.Cellz = CC;
    
setSceneAndTime;
   
end
function eliminatebutton_Callback(~,~)

global ImageDetails frameToLoad Tracked imgsize

button=1;
while button == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
     [cellxx,cellyy,button] = ginput(1);

     

       
% for j = 1:length(cellxx)
%     cellx = cellxx(j);
%     celly = cellyy(j);
%     
    cellx = round(cellxx);
    celly = round(cellyy);
    
    cellind = sub2ind(imgsize,celly,cellx);
    
    for j = 1:length(cellxx)
    if j==1
        idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
    else
        idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
    end
    end
    
PX(~idxs) = {NaN};
CC.PixelIdxList = PX;
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Tracked{t}.Cellz = CC;

if button==1
    nextbutton_callback([],[])
elseif button == 3
    prevbutton_callback([],[])
else 
   setSceneAndTime
end

end
end

function destroyAllFramePrevious_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked imgsize

%   determine the frame to load
    t = ImageDetails.Frame;
    CC = Tracked{t}.Cellz;
    PX = CC.PixelIdxList;   

    idxs = false(size(PX));   %choose all cells on frame
    Trackedz = crushThem(Tracked,idxs,[],t); %delete from t=t backwar
    Tracked = Trackedz;
       setSceneAndTime

end
function destroyAllFrameSubsequent_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked imgsize

%   determine the frame to load
    t = ImageDetails.Frame;
    CC = Tracked{t}.Cellz;
    PX = CC.PixelIdxList;   

    idxs = false(size(PX));   %choose all cells on frame
    Trackedz = crushThem(Tracked,idxs,t,[]); %delete from t=t onward
    Tracked = Trackedz;
       setSceneAndTime
end
function destroybuttonAllPrevious_Callback(~,~)
%delete a cell from all frames
global ImageDetails frameToLoad Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
Trackedz = crushThem(Tracked,idxs,[],t);      
Tracked = Trackedz;


   setSceneAndTime



end
function destroybuttonAllSubsequent_Callback(~,~)
%delete a cell from all frames
global ImageDetails frameToLoad Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
Trackedz = crushThem(Tracked,idxs,t,[]);      
Tracked = Trackedz;


   setSceneAndTime



end

function deleteAllonFrame_Callback(~,~)
%delete all cells from one frame
  global  ImageDetails  Tracked
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.

       
% for j = 1:length(cellxx)
%     cellx = cellxx(j);
%     celly = cellyy(j);
    
    idxs = false(size(PX));
    
PX(~idxs) = {NaN};
CC.PixelIdxList = PX;
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Tracked{t}.Cellz = CC;
    
setSceneAndTime;

end



function destroybutton_Callback(~,~)
%delete a cell from all frames
global ImageDetails frameToLoad Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
Trackedz = crushThem(Tracked,idxs,[],[]);      
Tracked = Trackedz;


   setSceneAndTime


end
%choose the cells you want
function chosenOnes_Callback(~,~)
%choose the cells you want
global ImageDetails frameToLoad Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      

      
% Trackedz = crushThem(Tracked,~idxs,length(Tracked),length(Tracked));
Trackedz = crushThem(Tracked,~idxs,1,length(Tracked)); 
Tracked = Trackedz;


   setSceneAndTime


end
function chosenOnesAllOnFrame_Callback(~,~)
%choose the cells you want
global ImageDetails frameToLoad Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
  idxs = ~cellfun(@(x) length(x)>1,PX,'UniformOutput',1); %choose all cells on frame
  
% Trackedz = crushThem(Tracked,~idxs,length(Tracked),length(Tracked));
Trackedz = crushThem(Tracked,~idxs,1,length(Tracked)); 
Tracked = Trackedz;


   setSceneAndTime


end
function chosenOnesEnd_Callback(~,~)
%choose the cells you want
global ImageDetails frameToLoad Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
Trackedz = crushThem(Tracked,~idxs,length(Tracked),length(Tracked));
% Trackedz = crushThem(Tracked,~idxs,1,length(Tracked)); 
Tracked = Trackedz;


   setSceneAndTime

end

function erodeOnes_Callback(~,~)
%choose the cells you want
global ImageDetails frameToLoad Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
% Trackedz = crushThem(Tracked,~idxs,length(Tracked),length(Tracked));
% Trackedz = crushThem(Tracked,~idxs,1,length(Tracked)); 


initialframe = 1;
finalframe = length(Tracked);
Stacked = Tracked;



if isempty(initialframe)
    initialframe=1;
end

if isempty(finalframe)
    finalframe=length(Stacked);
end

se = strel('disk',1);
idxf = find((~idxs) ==1);
for i=initialframe:finalframe 
CC = Stacked{i}.Cellz;
PX = CC.PixelIdxList;   
    for jim = idxf
    pixxies = PX{jim};
        if ~isnan(pixxies)
        imdub = zeros(imgsize);
        imdub(pixxies) = 1;
        imdub = imerode(imdub,se);
        pxtwo = find(imdub==1);
            if ~isempty(pxtwo)
            PX{jim} = pxtwo;
            else
            PX{jim} = NaN;
            end
        end
    end
CC.PixelIdxList = PX;
CC.NumObjects = length(PX);
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Stacked{i}.Cellz = CC;

end
Trackedz=Stacked;








Tracked = Trackedz;


   setSceneAndTime


end



function dilateOnes_Callback(~,~)
%choose the cells you want
global ImageDetails frameToLoad Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
% Trackedz = crushThem(Tracked,~idxs,length(Tracked),length(Tracked));
% Trackedz = crushThem(Tracked,~idxs,1,length(Tracked)); 


initialframe = 1;
finalframe = length(Tracked);
Stacked = Tracked;



if isempty(initialframe)
    initialframe=1;
end

if isempty(finalframe)
    finalframe=length(Stacked);
end

se = strel('disk',1);
idxf = find((~idxs) ==1);
for i=initialframe:finalframe 
CC = Stacked{i}.Cellz;
PX = CC.PixelIdxList;   
    for jim = idxf
    pixxies = PX{jim};
%     disp(pixxies)
        if ~isnan(pixxies)
        imdub = zeros(imgsize);
        imdub(pixxies) = 1;
        imdub = imdilate(imdub,se);
        pxtwo = find(imdub==1);
            if ~isempty(pxtwo)
            PX{jim} = pxtwo;
            else
            PX{jim} = NaN;
            end
        end
    end
CC.PixelIdxList = PX;
CC.NumObjects = length(PX);
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Stacked{i}.Cellz = CC;

end
Trackedz=Stacked;








Tracked = Trackedz;


   setSceneAndTime


end

%internal function for chosenOnes and destroy
function Trackedz = crushThem(Stacked,idxs,initialframe,finalframe)  

if isempty(initialframe)
    initialframe=1;
end

if isempty(finalframe)
    finalframe=length(Stacked);
end

for i=initialframe:finalframe 
CC = Stacked{i}.Cellz;
PX = CC.PixelIdxList;   
PX(~idxs) = {NaN};
CC.PixelIdxList = PX;
CC.NumObjects = length(PX);
Stacked{i}.Cellz = CC;
end
Trackedz=Stacked;
end




%add cells and link cells
function addareabutton_Callback(~,~) 
 global  ImageDetails frameToLoad Tracked imgsize
 % choose cell
%       [cellx,celly] = ginput(1);
       % construct a polygon to add

       If = 1;
       button=1;
       while button==1
      [polyx,polyy,button] = ginput();
      button = round(mean(button));
      
      if button ==1
      M = zeros(1,length(polyx)*2);
      M(1:2:end) = polyx;
      M(2:2:end) = polyy;
      zeroImage = zeros(imgsize);
      zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
      zerogray = rgb2gray(zeroImage);

      if isempty(Tracked{1}.Cellz)

          
      else  %if there exists segmenttracking already...then load that. 
        
      imagio = zeros(imgsize);
      imagio(zerogray>0)=1;
      cc = bwconncomp(imagio);
      px = cc.PixelIdxList;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   determine the frame to load
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;        
        
idxs = cellfun(@(x) sum(ismember(x,px{1})),PX,'UniformOutput',1);
index = find(idxs>1);
    if ~isempty(index)
    newMass = vertcat(PX{index});
    PX{min(index)} = unique(vertcat(newMass,px{1}));
        if length(index)>1
        index(find(index == min(index)))=[];
        PX(index) = {NaN};
        end
    CC.PixelIdxList = PX;
    else
    CC.PixelIdxList = horzcat(PX,px);    
    end
CC.NumObjects = length(CC.PixelIdxList);
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Tracked{t}.Cellz = CC;
% 
% makeIMG = cellfun(@(x) length(x)==1,PX,'UniformOutput',1); %choose only the cells without NAN
% newpxlist = horzcat(PX(~makeIMG),px);
% CC.PixelIdxList = newpxlist;
% CC.NumObjects = length(newpxlist);
% 
% segmentimgL = labelmatrix(CC);
% axes(MainAxes)
% imagesc(segmentimgL);
% segmentimg = zeros(size(segmentimgL));
% segmentimg(segmentimgL>0)=1;
% If = segmentimg;
    
      end
      
      
      
    nextbutton_callback([],[]);
      end 
      end
end
function linkCells_Callback(~,~)
global ImageDetails frameToLoad Tracked TC imgsize timeFrames

%this is a quick workaround to get linking working with tracking
%trajectories on
% TC =1;
% setSceneAndTime

button=1;
i=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while button==1
TC=0;
    [xx,yy,button] = ginput(1); %record each click
    
    if button ==1
    t = ImageDetails.Frame;
    cellxx(i)  =    xx; 
    cellyy(i)  =    yy;
    timingkeeper(i) = t;
        if t==timeFrames
            i;
        else
            i=i+1;
        end
    nextbutton_callback([],[]);
    end
    
end
TC=1;

cellx = round(cellxx);
celly = round(cellyy);
cellind = sub2ind(imgsize,celly,cellx);

for i = 1:length(timingkeeper)
t = timingkeeper(i);
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;  
idxlog = cellfun(@(x) isempty(find(x==cellind(i),1)),PX,'UniformOutput',1);
idx(i) = find(idxlog==0);
end
% 
% cellindx = 1:length(frameToLoad);
% timeindx = 1:length(frameToLoad);
cellindx = 1:timeFrames;
timeindx = 1:timeFrames;
cellindx(1:timingkeeper(1))=idx(1);
cellindx(timingkeeper(end):end) = idx(end);

for i=timingkeeper
cellindx(i) = idx(i-(min(timingkeeper)-1));  
end


for i=1:timeFrames

CC = Tracked{i}.Cellz;
PX = CC.PixelIdxList;
PX{cellindx(1)} = PX{cellindx(i)};
if cellindx(1)~=cellindx(i)
PX{cellindx(i)} = NaN;
end
CC.PixelIdxList = PX;
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Tracked{i}.Cellz = CC;
end

setSceneAndTime

end

%plot your cells!
function psettings = PlotSettings_callback(~,~)
global exportdir OGExpDate frameToLoad
cd(exportdir)

queryName = strcat(OGExpDate,'*DoseAndScene*.mat');
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


function PlotCFPnotnorm_callback(~,~)
global toggleCFPnorm
% tcontrast =99;
% lcontrast =1;
toggleCFPnorm = 0;
Plot_callback([],[]);
end
function PlotCFPnorm_callback(~,~)
global toggleCFPnorm
% tcontrast =99;
% lcontrast =1;
toggleCFPnorm = 1;
Plot_callback([],[]);
end

function Plot_callback(~,~)
global cell_seg nucleus_seg xAxisLimits trunccmaplz toggleCFPnorm SecondPlotAxes Tracked ImageDetails expDirPath trackingPath timeFrames mstackPath frameToLoad PlotAxes imgsize plottingON psettings cmaplz displaytracking cmap

    if plottingON == 0
        psettings = PlotSettings_callback([],[]);
        plottingON=1;
    end
    framesThatMustBeTracked = psettings.framesThatMustBeTracked;


    for jy = 1
        PX = Tracked{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
        makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
    end
    [comments,commentpos,cellidx,plotidx]=commentsforplot(Tracked);
    makeIMGidx = find(makeIMG==1);
    
    if ~isempty(plotidx)
        iidd = find(~ismember(makeIMGidx,plotidx));
    else
        iidd=[];
    end


    smooththat=0;
    [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,expDirPath,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);

    % plotStructUI.EGFP = Smad;
    % plotStructUI.mKate = mkate;
    % plotStructUI.Cfp = Cfp;
    % plotStructUI.CfpFC = CfpFC;
    % plotStructUI.SmadFC = SmadFC;
    % plotStructUI.mkateFC = mkateFC;
    % plotStructUI.Smadbkg = Smadbkg;
    % plotStructUI.Cfpbkg = Cfpbkg;
    % plotStructUI.mkatebkg = mkatebkg;

    if strcmp(ImageDetails.Channel,cell_seg)
        plotMat = plotStructUI.Smad;
        plotMatFC = plotStructUI.SmadFC;
        ylimit =[0 6];
    elseif strcmp(ImageDetails.Channel,nucleus_seg)
        plotMat = plotStructUI.mkatetotal;
        plotMatFC = plotStructUI.mkateFCtotal;
        ylimit =[0 3];
    else
        plotMat = plotStructUI.Smad;
        plotMatFC = plotStructUI.SmadFC;
        ylimit = [0 10];
    end
    
    xmin = 0;
    toplot = plotMatFC;
    idx = true(size(toplot,1),1);
    cmapl = trunccmaplz;
    idxa = find(idx==1);
    h = plot(SecondPlotAxes,toplot(idx,:)','LineWidth',2);
            
    if displaytracking ==1
        for i=1:length(h)
            h(i).Color = cmapl(idxa(i),:);
        end
        colormap(cmap);    
    end
    
    for po = iidd  %fade out the commented cells
        h(po).LineStyle = ':';
    end

    SecondPlotAxes.XLim = ([xAxisLimits(1) xAxisLimits(2)]);
    SecondPlotAxes.YLim = (ylimit);
    SecondPlotAxes.YLabel.String = 'fold change';
    SecondPlotAxes.XLabel.String = 'frames';
    SecondPlotAxes.XGrid = 'on';
    SecondPlotAxes.YGrid = 'on';
    SecondPlotAxes.Color = [0.95 0.95 0.95];


    toplot = plotMat;
    h = plot(PlotAxes,toplot(idx,:)','LineWidth',2);
    
    if displaytracking ==1
        for i=1:length(h)
            h(i).Color = cmapl(i,:);
        end
        colormap(cmap);    
    end
    
    for po = iidd  %fade out the commented cells
        h(po).LineStyle = ':';
    end
        

    PlotAxes.XLim = ([xAxisLimits(1) xAxisLimits(2)]);
    PlotAxes.YLim = ([prctile(toplot(:),0.5)./1.2 prctile(toplot(:),99.5).*1.2]);
    PlotAxes.YLabel.String = 'abundance';
    PlotAxes.XLabel.String = 'frames';
    PlotAxes.XGrid = 'on';
    PlotAxes.YGrid = 'on';
    PlotAxes.Color = [0.95 0.95 0.95];

end


function Plot_SpecificCell_callback(~,~)
global displaytracking trunccmaplz xAxisLimits plottingTotalOrMedian ThirdPlotAxes toggleCFPnorm Tracked ImageDetails expDirPath mstackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings

if plottingON == 0
psettings = PlotSettings_callback([],[]);
plottingON=1;
end




framesThatMustBeTracked = psettings.framesThatMustBeTracked;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
      
% Display mesh plot of the currently selected data.
    [cellxx,cellyy] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    cellind = sub2ind(imgsize,celly,cellx);
        for j = 1:length(cellxx)
            if j==1
                idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
            else
                idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
            end
        end
makeIMG = ~idxs;        
makeIMGidx = find(makeIMG==1);


smooththat=0;
[plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,expDirPath,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
smooththat=toggleCFPnorm;
SmadFC = plotStructUI.SmadFC;
    
    if strcmp(ImageDetails.Channel,'EGFP')
        plotMat = plotStructUI.Smad;
        plotMatFC = plotStructUI.SmadFC;
        ylimit =[0 6];
    elseif strcmp(ImageDetails.Channel,'mKate')
        plotMat = plotStructUI.mkatetotal;
        plotMatFC = plotStructUI.mkateFCtotal;
        ylimit =[0 3];
    else
        plotMat = plotStructUI.Smad;
        plotMatFC = plotStructUI.SmadFC;
    end
    
toplot = plotMatFC(makeIMGidx,:);
h = plot(ThirdPlotAxes,toplot','LineWidth',3);
ThirdPlotAxes.XLim = ([xAxisLimits(1) xAxisLimits(2)]);
ThirdPlotAxes.YLim = (ylimit);
ThirdPlotAxes.YLabel.String = 'fold change';
ThirdPlotAxes.XLabel.String = 'frames';
ThirdPlotAxes.XGrid = 'on';
ThirdPlotAxes.YGrid = 'on';
ThirdPlotAxes.Color = [0.95 0.95 0.95];


    cmapl = trunccmaplz;
            
    if displaytracking ==1
        for i=1:length(h)
            h(i).Color = cmapl(makeIMGidx(i),:);
        end 
    end

end

function [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,expDirPath,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat)
global plottingTotalOrMedian cell_seg nucleus_seg background_seg segmentPath

PXX = Tracked{1}.Cellz.PixelIdxList;
plotTracesCell = cell(length(PXX),length(Tracked));
    for i = 1:length(Tracked)
        PXX = Tracked{i}.Cellz.PixelIdxList;
        plotTracesCell(:,i) = PXX;
    end
    

cd(mstackPath)
%no bleach correction option yet
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%   open the image files   %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %  open smad img  %
            cd(mstackPath)         
            ff = dir(strcat('*',ImageDetails.Scene,'*',cell_seg,'*'));
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            cellQ_imgstack = channelfileObject.flatstack;
            
                     %  open mKate img  %
            cd(mstackPath)         
            ff = dir(strcat('*',ImageDetails.Scene,'*',nucleus_seg,'*'));
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            mkateimgstack = channelfileObject.flatstack;


                   %    open cfp img  %
            cd(mstackPath)         
            ff = dir(strcat('*',ImageDetails.Scene,'*',nucleus_seg,'*'));
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            nuc_imgstack = channelfileObject.flatstack;
            
                    % open background Logical img  %
            cd(segmentPath)         
            ff = dir(strcat('*',ImageDetails.Scene,'*',background_seg,'*'));
            if length(ff)>1
                ff = dir(strcat('*',ImageDetails.Scene,'*',background_seg,'*background*'));
            end
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            bkglogimgstack = channelfileObject.IfFinal;
            
                    % open nuclear Logical img  %
            cd(segmentPath)         
            ff = dir(strcat('*',ImageDetails.Scene,'*',nucleus_seg,'*'));
            if length(ff)>1
                ff = dir(strcat('*',ImageDetails.Scene,'*',nucleus_seg,'*nucleus*'));
            end
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            nuclogimgstack = channelfileObject.IfFinal;
            bkglogimgstack(nuclogimgstack) = true; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           

%perform bkg subtraction
Smadbkg = zeros(1,timeFrames,'single');
Cfpbkg = zeros(1,timeFrames,'single');
mkatebkg = zeros(1,timeFrames,'single');
bkgstd = zeros(1,timeFrames,'single');
for k=1:timeFrames
    bkglog = ~bkglogimgstack(:,:,k);
    cellQ_img = single(cellQ_imgstack(:,:,k));
    nuc_img = single(nuc_imgstack(:,:,k));
    mkateimg = single(mkateimgstack(:,:,k));
    %background subtraction is just subtraction with a value
    Smadbkg(k) = nanmedian(cellQ_img(bkglog));
    Cfpbkg(k) = nanmedian(nuc_img(bkglog));
    mkatebkg(k) = nanmedian(mkateimg(bkglog));
    bkgstd(k) = nanstd(mkateimg(bkglog));
    cellQ_imgstack(:,:,k) = cellQ_img-Smadbkg(k);
    nuc_imgstack(:,:,k) = nuc_img-Cfpbkg(k);
    mkateimgstack(:,:,k) = mkateimg - mkatebkg(k);
    %background subtraction is subtraction with an interpolated image
%     smadbkgimg = regionfill(cellQ_img,~bkglog);
%     cfpbkgimg = regionfill(nuc_img,~bkglog); %fill in the regions where bkglog is 0
%     cellQ_imgstack(:,:,k) = cellQ_imgstack(:,:,k)-smadbkgimg;
%     nuc_imgstack(:,:,k) = nuc_imgstack(:,:,k)-cfpbkgimg;
end


%extract pixel intensities
cellQ_pxls = cell(size(plotTracesCell,1),size(plotTracesCell,2));
nuc_pxls = cell(size(plotTracesCell,1),size(plotTracesCell,2));

for i = 1:size(plotTracesCell,2)
    cellQ_img = single(squeeze(cellQ_imgstack(:,:,i)));
    nuc_img = single(squeeze(nuc_imgstack(:,:,i)));
    mkateimg = single(squeeze(mkateimgstack(:,:,i)));
    for j=1:size(plotTracesCell,1)
    pxidx = plotTracesCell{j,i};
        if ~isnan(pxidx)
        cellQ_pxls(j,i) = {cellQ_img(pxidx)};
        nuc_pxls(j,i) = {nuc_img(pxidx)};
        mkatepxls(j,i) = {mkateimg(pxidx)};
        else
        cellQ_pxls(j,i) = {single(13579)};
        nuc_pxls(j,i) = {single(13579)};
        mkatepxls(j,i) = {single(13579)};
        end
    end
end


    Smadtotal = cellfun(@nansum,cellQ_pxls,'UniformOutput',1);
    % Smad = cellfun(@nanmean,cellQ_pxls,'UniformOutput',1);
    Smadtotal(Smadtotal==single(13579)) = NaN;
    Cfptotal = cellfun(@nansum,nuc_pxls,'UniformOutput',1);
    Cfptotal(Cfptotal==single(13579)) = NaN;
    mkatetotal = cellfun(@nansum,mkatepxls,'UniformOutput',1);
    mkatetotal(mkatetotal==single(13579)) = NaN;

    Smad = cellfun(@nanmedian,cellQ_pxls,'UniformOutput',1);
    % Smad = cellfun(@nanmean,cellQ_pxls,'UniformOutput',1);
    Smad(Smad==single(13579)) = NaN;
    Cfp = cellfun(@nanmedian,nuc_pxls,'UniformOutput',1);
    Cfp(Cfp==single(13579)) = NaN;
    mkate = cellfun(@nanmedian,mkatepxls,'UniformOutput',1);
    mkate(mkate==single(13579)) = NaN;


basalSUB = framesThatMustBeTracked(1)-7;
if basalSUB<1
    basalSUB = framesThatMustBeTracked-1;
end
basalcfp = nanmean(Cfp(:,framesThatMustBeTracked(1)-basalSUB:framesThatMustBeTracked(1)),2);

CfpFC = zeros(size(Cfp),'single');
    for i = 1:size(Cfp,2)
       CfpFC(:,i) = Cfp(:,i)./basalcfp; 
    end
    
    if smooththat==1;
    %%%%%%%%%%%%%%%%%%%
    Smad = Smad./CfpFC;
    %%%%%%%%%%%%%%%%%%%
    end
basalsmad = nanmean(Smad(:,framesThatMustBeTracked(1)-basalSUB:framesThatMustBeTracked(1)),2);
SmadFC = zeros(size(Smad),'single');
    for i = 1:size(Smad,2)
       SmadFC(:,i) = Smad(:,i)./basalsmad; 
    end
    
basalmkate = nanmean(mkate(:,framesThatMustBeTracked(1)-basalSUB:framesThatMustBeTracked(1)),2);
mkateFC = zeros(size(mkate),'single');
    for i = 1:size(mkate,2)
       mkateFC(:,i) = mkate(:,i)./basalmkate; 
    end
    
basalmkatetotal = nanmean(mkatetotal(:,framesThatMustBeTracked(1)-basalSUB:framesThatMustBeTracked(1)),2);
mkateFCtotal = zeros(size(mkatetotal),'single');
    for i = 1:size(mkatetotal,2)
       mkateFCtotal(:,i) = mkatetotal(:,i)./basalmkatetotal; 
    end
    
%     Smad,Cfp,mkate,CfpFC,SmadFC,mkateFC,Smadbkg,Cfpbkg,mkatebkg
plotStructUI.Smad = Smad;
plotStructUI.mkate = mkate;
plotStructUI.mkatetotal = mkatetotal;
plotStructUI.Cfp = Cfp;
plotStructUI.CfpFC = CfpFC;
plotStructUI.SmadFC = SmadFC;
plotStructUI.mkateFC = mkateFC;
plotStructUI.mkateFCtotal = mkateFCtotal;
plotStructUI.Smadbkg = Smadbkg;
plotStructUI.Cfpbkg = Cfpbkg;
plotStructUI.mkatebkg = mkatebkg;
    
end
function plotStruct = plotthemfunctionToStructure(Tracked,idScene,mstackPath,timeFrames,makeIMG,makeIMGidx,cell_seg,nucleus_seg,background_seg,segmentPath)

plotStruct = struct();



plotTracesCell = cell(length(makeIMGidx),length(Tracked));
centroidarray = cell(length(makeIMGidx),length(Tracked));
    for i = 1:length(Tracked)
        CC = Tracked{i}.Cellz;
        stats = regionprops(CC,'Centroid');
        centroids = {stats.Centroid};
        centroidarray(:,i) = centroids(makeIMG);
        PXX = Tracked{i}.Cellz.PixelIdxList;
        plotTracesCell(:,i) = PXX(makeIMG);
    end
    

cd(mstackPath)
%no bleach correction option yet
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%   open the image files   %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %  open smad img  %
            cd(mstackPath)         
            ff = dir(strcat('*',idScene,'*',cell_seg,'*'));
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            cellQ_imgstack = channelfileObject.flatstack;

                   %    open cfp img  %
            cd(mstackPath)         
            ff = dir(strcat('*',idScene,'*',nucleus_seg,'*'));
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            nuc_imgstack = channelfileObject.flatstack;
            
                    % open background Logical img  %
            cd(segmentPath)         
            ff = dir(strcat('*',idScene,'*',background_seg,'*'));
            if length(ff)>1
                ff = dir(strcat('*',idScene,'*',background_seg,'*background*'));
            end
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            bkglogimgstack = channelfileObject.IfFinal;
            
                    % open nuclear Logical img  %
            cd(segmentPath)         
            ff = dir(strcat('*',idScene,'*',nucleus_seg,'*'));
            if length(ff)>1
                ff = dir(strcat('*',idScene,'*',nucleus_seg,'*nucleus*'));
            end
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            nuclogimgstack = channelfileObject.IfFinal;
            bkglogimgstack(nuclogimgstack) = true; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
cellQ_pxls = cell(size(plotTracesCell,1),size(plotTracesCell,2));
nuc_pxls = cell(size(plotTracesCell,1),size(plotTracesCell,2));

for i = 1:size(plotTracesCell,2)
    cellQ_img = single(squeeze(cellQ_imgstack(:,:,i)));
    nuc_img = single(squeeze(nuc_imgstack(:,:,i)));
    for j=1:size(plotTracesCell,1)
    pxidx = plotTracesCell{j,i};
        if ~isnan(pxidx)
        cellQ_pxls(j,i) = {cellQ_img(pxidx)};
        nuc_pxls(j,i) = {nuc_img(pxidx)};
        else
        cellQ_pxls(j,i) = {single(13579)};
        nuc_pxls(j,i) = {single(13579)};
        end
    end
end

%determine median pxl intensities
    cellQ = cellfun(@nanmedian,cellQ_pxls,'UniformOutput',1);
        cellQ(cellQ==single(13579)) = NaN;
    nucQ = cellfun(@nanmedian,nuc_pxls,'UniformOutput',1);
        nucQ(nucQ==single(13579)) = NaN;
for i = 1:size(cellQ_pxls,1)
    plotStruct(i).medianNucEGFP = cellQ(i,:);
    plotStruct(i).medianNucRFP = nucQ(i,:);
end

%determine total pxl intensities
    cellQ = cellfun(@nansum,cellQ_pxls,'UniformOutput',1);
        cellQ(cellQ==single(13579)) = NaN;
    nucQ = cellfun(@nansum,nuc_pxls,'UniformOutput',1);
        nucQ(nucQ==single(13579)) = NaN;
for i = 1:size(cellQ_pxls,1)
    plotStruct(i).totalNucEGFP = cellQ(i,:);
    plotStruct(i).totalNucRFP = nucQ(i,:);
end

%determine mean pxl intensities
    cellQ = cellfun(@nanmean,cellQ_pxls,'UniformOutput',1);
    cellQ(cellQ==single(13579)) = NaN;
    nucQ = cellfun(@nanmean,nuc_pxls,'UniformOutput',1);
    nucQ(nucQ==single(13579)) = NaN;
for i = 1:size(cellQ_pxls,1)
    plotStruct(i).meanNucEGFP = cellQ(i,:);
    plotStruct(i).meanNucRFP = nucQ(i,:);
    plotStruct(i).medianCfpbkg = nucBKG;
    plotStruct(i).medianSmadbkg = cellBKG;
end


%determine mean pxl intensities
    cellQ = cellfun(@nanmean,cellQ_pxls,'UniformOutput',1);
    cellQ(cellQ==single(13579)) = NaN;
    nucQ = cellfun(@nanmean,nuc_pxls,'UniformOutput',1);
    nucQ(nucQ==single(13579)) = NaN;
for i = 1:size(cellQ_pxls,1)
    plotStruct(i).meanNucEGFP = cellQ(i,:);
    plotStruct(i).meanNucRFP = nucQ(i,:);
    plotStruct(i).medianCfpbkg = nucBKG;
    plotStruct(i).medianSmadbkg = cellBKG;
end

for i = 1:size(centroidarray,1)
    plotStruct(i).Centroid = centroidarray(i,:);
end





end



%%%%comments!!!
function xy = getxy(~,~)
global Tracked plottingON psettings imgsize ImageDetails frameToLoad
    
    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
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
        y = rem(px-1,imgsize(1))+1; %these two lines replace ind2sub
        x = (px-y)/imgsize(2) + 1;  %these two lines replace ind2sub
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     t=length(Tracked);
t = ImageDetails.Frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
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


function [comments,commentpos,cellidx,plotidx]=commentsforplot(Tracked)
%choose the cells you want to comment on
global ImageDetails frameToLoad imgsize 

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
if sum(strcmp(fnames,'comments'));
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

comments = cell(1,length(oldcommentpos));
commentpos = zeros(1,length(oldcommentpos));
cellidx =zeros(1,length(oldcommentpos));

indies = oldcommentpos;

if ~isempty(oldcommentpos)
    for jim = idx
        cycle=0;
        for i = 1:length(oldcommentpos)
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


for i=1:length(comments)
    cellnums{i} = strcat('cell#',num2str(i));
    if isempty(comments{i})
        comments{i}='';
    end
end

stringsarray = {'overdriver','no cfp','nocfp','dimmer',...
            'overlap','remove','doublenuc','saturated','big'};
        for ubi = 1:length(stringsarray)
            if ubi==1
            stringer = strcat('(',stringsarray{ubi},'|');
            elseif ubi==length(stringsarray)
                stringer = horzcat(stringer,strcat(stringsarray{ubi},')'));
            else
                stringer = horzcat(stringer,strcat(stringsarray{ubi},'|'));
            end
        end
        
[o,oo,~,~] = regexp(comments,stringer);
overlay = cellfun(@length,o,'UniformOutput',1);
plotidx = cellidx(~overlay);
stophere=1;

end
function setcommentsTracking(comments,commentpos)
global Tracked
t=length(Tracked);
Tracked{t}.comments= comments;
Tracked{t}.commentpos = commentpos;
end
function [cellnums,comments,commentpos,cellidx]=updatecomments(xy)
%choose the cells you want to comment on
global ImageDetails frameToLoad Tracked imgsize 

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
if sum(strcmp(fnames,'comments'));
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

indies = sub2ind([imgsize],xy(:,2),xy(:,1));

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

function comments = nowaddcommentsJ(xy,cellnums,comments,commentpos)
global Tracked
t=length(Tracked);



prompt = cellnums;
dlg_title = 'Comment for selected cells:';
num_lines = 1;
defaultans = comments;
% answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
answer = defaultans;

if ~isempty(answer)
comments = answer;

Tracked{t}.comments= comments;
Tracked{t}.commentpos = commentpos;
end
% displaycomments =1;

end

function comment_CallbackJ(~,~)
global ExportNameKey ExportName OGExpDate displaycomments SceneList Tracked ImageDetails expDirPath trackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(trackPath)
cd ..
    for scenenumber = 1:length(SceneList)
        cd(trackPath)
        cd ..
        sceneN = char(SceneList{scenenumber});
        disp(sceneN)
        scenedir = dir(strcat('*',sceneN,'*'));
        scenedirname = char({scenedir.name});
        cd(scenedirname)
        SceneDirPath = char({pwd});
        
        
        trackfile = dir(strcat(ExportNameKey,ExportName,'.mat'));
%         trackfile = dir('finalfricktrack.mat');
        trackfilename = char({trackfile.name});
        
            if ~isempty(trackfilename)
                load(trackfilename)
                fnames = fieldnames(Tracked{length(Tracked)});
                    if sum(strcmp(fnames,'comments'))
                        displaycomments =1;
                    else
                        displaycomments =0;
                    end
                xy = getxy([],[]);
                dispxy(xy)
                [cellnums,comments,commentpos,cellidx]=updatecomments(xy);
                comments = nowaddcommentsJ(xy,cellnums,comments,commentpos);
                displaycomments =1;
                dispxy(xy)   
                saveTrackingFileAs_callbackJ([],[],SceneDirPath)
            else
                stopmehere=1;
            end
    end
      
end

function exportTrackedCells(~,~)
global cell_seg nucleus_seg background_seg segmentPath ExportNameKey ExportName exportdir mstackPath OGExpDate SceneList ImageDetails expDirPath trackingPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings
exportStruct = struct();

    if plottingON == 0
        psettings = PlotSettings_callback([],[]);
        plottingON=1;
    end
    
%remove all global variables before parfor loop
    framesThatMustBeTracked = psettings.framesThatMustBeTracked;
    eNameKey = ExportNameKey;
    eName = ExportName;
    tPath = trackingPath;
    mPath = mstackPath;
    sList = SceneList;
    tFrames = timeFrames;
    cSeg = cell_seg;
    nSeg = nucleus_seg;
    bSeg = background_seg;
    sPath = segmentPath;
    cd(tPath)
    cd ..
    
        plotStructArray = cell(1,length(sList));
%         for scenenumber = 1:length(sList)
        parfor scenenumber = 1:length(sList)
            cd(tPath)
            sceneN = sList{scenenumber};
            disp(sceneN)
            idScene = sceneN;

            trackfile = dir(strcat(eNameKey,'*',sceneN,'*',eName,'.mat'));
            trackfilename = char({trackfile.name});

                if ~isempty(trackfilename)
                    trackedArray = loadTrackedArray(trackfilename); %trackedArray = Tracked
                    PX = trackedArray{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
                    makeIMG = false(length(framesThatMustBeTracked),length(PX));

                        for jy = 1:length(framesThatMustBeTracked)
                        PX = trackedArray{framesThatMustBeTracked(jy)}.Cellz.PixelIdxList;
                        makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
                        end

                    makeIMG = makeIMG(1,:)&makeIMG(2,:);
                    makeIMGidx = find(makeIMG==1);

                    plotStruct = plotthemfunctionToStructure(trackedArray,idScene,mPath,tFrames,makeIMG,makeIMGidx,cSeg,nSeg,bSeg,sPath);
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
                                exportStruct(idx+i).scene = sceneN;
                                exportStruct(idx+i).cellID = i;
                            end
                        end

                    else    %if fields are defined, append the cell data to the next available index

                        idx = length(exportStruct);
                        for i = 1:length(plotStruct)
                            for j = 1:length(fnames)
                                exportStruct(idx+i).(fnames{j}) = plotStruct(i).(fnames{j});
                                exportStruct(idx+i).scene = sceneN;
                                exportStruct(idx+i).cellID = i;
                            end
                        end

                    end

        end
        
        %save the exportStruct
            cd(exportdir)
            filename = strcat(OGExpDate,'_tracking_export.mat'); 
            save(filename,'exportStruct');
end


function trackedArray = loadTrackedArray(trackfilename)
% tic
load(trackfilename)
trackedArray = Tracked;
% toc
% tic
% tfileobject = matfile(trackfilename);
% trackedArray = tfileobject.Tracked;
% toc
end

function exportSegmentedCells(~,~)
global cell_seg nucleus_seg background_seg segmentPath   exportdir mstackPath OGExpDate SceneList  trackingPath timeFrames


    exportNucleiStruct=struct();
    %remove all global variables before parfor loop
        tPath = trackingPath;
        mPath = mstackPath;
        sList = SceneList;
        tFrames = timeFrames;
        cSeg = cell_seg;
        nSeg = nucleus_seg;
        bSeg = background_seg;
        sPath = segmentPath;
        cd(tPath)
        cd ..
    
        nucleiStructArray = cell(1,length(sList));
        parfor scenenumber = 1:length(sList)
            cd(tPath)
            sceneN = sList{scenenumber};
            disp(sceneN)
            idScene = sceneN;

                    exportNucleiz = nucleiToStructure(idScene,mPath,tFrames,cSeg,nSeg,bSeg,sPath);
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
            cd(exportdir)
            filename = strcat(OGExpDate,'_nuclei_export.mat'); 
            save(filename,'exportNucleiStruct');
            

end

function exportNucleiStruct = nucleiToStructure(idScene,mstackPath,timeFrames,cell_seg,nucleus_seg,background_seg,segmentPath)
cd(mstackPath)
%no bleach correction option yet
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%   open the image files   %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %  open smad img  %
        cd(mstackPath)         
        ff = dir(strcat('*',idScene,'*',cell_seg,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
        filename = char(ff.name);
        channelfileObject = matfile(filename);
        cellQ_imgstack = channelfileObject.flatstack;

               %    open cfp img  %
        cd(mstackPath)         
        ff = dir(strcat('*',idScene,'*',nucleus_seg,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
        filename = char(ff.name);
        channelfileObject = matfile(filename);
        nuc_imgstack = channelfileObject.flatstack;

                % open background Logical img  %
        cd(segmentPath)         
        ff = dir(strcat('*',idScene,'*',background_seg,'*'));
        if length(ff)>1
            ff = dir(strcat('*',idScene,'*',background_seg,'*background*'));
        end
%         ff = dir(strcat(ImageDetails.Channel,'*'));
        filename = char(ff.name);
        channelfileObject = matfile(filename);
        bkglogimgstack = channelfileObject.IfFinal;

                % open nuclear Logical img  %
        cd(segmentPath)         
        ff = dir(strcat('*',idScene,'*',nucleus_seg,'*'));
        if length(ff)>1
            ff = dir(strcat('*',idScene,'*',nucleus_seg,'*nucleus*'));
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

        nucIntensities = cell(1,length(PX));
        cellIntensities = nucIntensities;
        for j=1:length(PX)
            pxidx = PX{j};
            nucIntensities{j} = nuc_img(pxidx);
            cellIntensities{j} = cellQ_img(pxidx);
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
                    
    end
end

function exportAllCells(~,~)
global ExportNameKey ExportName exportdir  cell_seg nucleus_seg background_seg segmentPath OGExpDate SceneList ImageDetails expDirPath mstackPath timeFrames trackingPath frameToLoad PlotAxes imgsize plottingON psettings
exportStruct = struct();

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
        
%remove all global variables before parfor loop
    framesThatMustBeTracked = psettings.framesThatMustBeTracked;
    eNameKey = ExportNameKey;
    eName = ExportName;
    tPath = trackingPath;
    mPath = mstackPath;
    sList = SceneList;
    tFrames = timeFrames;
    cSeg = cell_seg;
    nSeg = nucleus_seg;
    bSeg = background_seg;
    sPath = segmentPath;
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
%                 [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,expDirPath,ImageDetails,SceneDirPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
                plotStruct = plotthemfunctionToStructure(trackedArray,idScene,mPath,tFrames,makeIMG,makeIMGidx);
                
                fnames = fieldnames(plotStruct);
                if isempty(fieldnames(exportStruct)) %if exportStruct is empty 
                    idx = 0;
                    for i = 1:length(plotStruct)
                        for j = 1:length(fnames)
                            exportStruct(idx+i).(fnames{j}) = plotStruct(i).(fnames{j});
                            exportStruct(idx+i).scene = sceneN;
                            exportStruct(idx+i).cellID = i;
                        end
                    end
                else    %if fields are defined, append the cell data to the next available index
                    idx = length(exportStruct);
                    for i = 1:length(plotStruct)
                        for j = 1:length(fnames)
                            exportStruct(idx+i).(fnames{j}) = plotStruct(i).(fnames{j});
                            exportStruct(idx+i).scene = sceneN;
                            exportStruct(idx+i).cellID = i;
                        end
                    end
                end
            end
    end
        cd(exportdir)
        filename = strcat(OGExpDate,'_tracking_allcells_export.mat'); 
        save(filename,'exportStruct');
end

function xy = labelCells(~,~)


global ExportNameKey ExportName displaycomments  SceneList Tracked ImageDetails expDirPath trackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
t = ImageDetails.Frame;
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(trackPath)
cd ..

sceneN = char(ImageDetails.Scene);
disp(sceneN)
scenedir = dir(strcat('*',sceneN,'*'));
scenedirname = char({scenedir.name});
cd(scenedirname)
SceneDirPath = char({pwd});

trackfile = dir(strcat(ExportNameKey,ExportName,'.mat'));
% trackfile = dir('finalfricktrack.mat');
trackfilename = char({trackfile.name});

    if ~isempty(trackfilename)
        load(trackfilename)
        
        PX = Tracked{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
        makeIMG = false(length(framesThatMustBeTracked),length(PX));
            for jy = 1:length(framesThatMustBeTracked)
                PX = Tracked{framesThatMustBeTracked(jy)}.Cellz.PixelIdxList;
%                 makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)==1,PX,'UniformOutput',1)); %choose only the cells without NAN
                makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
            end
        makeIMG = makeIMG(1,:)&makeIMG(2,:);
        makeIMGidx = find(makeIMG==1);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PXX = Tracked{t}.Cellz.PixelIdxList;
    PX = PXX(makeIMG);
    mx = nan(1,length(PX));
    my = nan(1,length(PX));
        for j = 1:length(PX)
        px = PX{j};
        y = rem(px-1,imgsize(1))+1; %these two lines replace ind2sub
        x = (px-y)/imgsize(2) + 1;  %these two lines replace ind2sub
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(xy,1)
text(xy(i,1)+1,xy(i,2)+1,num2str(i),'Color',[0 0 0],'FontSize',12);hold on
text(xy(i,1),xy(i,2),num2str(i),'Color',[1 1 1],'FontSize',12);hold on
end

   xxyy=xy;
   xxyy(:,1) = xxyy(:,1)+10;
   xxyy(:,2) = xxyy(:,2)+10;
fnames = fieldnames(Tracked{t});
if sum(strcmp(fnames,'comments'))
    if t==length(frameToLoad)
        if displaycomments==1
            comments = Tracked{t}.comments;
            commentpos = Tracked{t}.commentpos;
%             for i = 1:size(xxyy,1)
%             for i = 1:length(comments)
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

end

function exportFrames(~,~)
global ExportNameKey ExportName  imgsize displaytracking SceneList Tracked ImageDetails expDirPath trackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
    
%   determine the frame to load
t = ImageDetails.Frame;


framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(trackPath)
cd ..
CENTROID = struct();
    
        cd(trackPath)
        cd ..
%         sceneN = char(SceneList{scenenumber});
%         disp(sceneN)
%         ImageDetails.Scene = sceneN;
        sceneN = ImageDetails.Scene;
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
            displaytracking =1;
%             finalbutton_callback([],[]);
            for i=1:length(frameToLoad)
                if i==1
                    firstbutton_callback([],[]);
                else
                    nextbutton_callback([],[]);
                end
                
%        xy = labelCells([],[]);
            set(gcf,'Units','Pixels');
            set(gca,'Units','Pixels');
            P = get(gca,'pos');
            F = getframe(gcf,P);
            [X,Map] = frame2im(F);
            set(gcf,'Units','normalized');
            set(gca,'Units','normalized');
            cd(trackPath)
            cd ..
            cd ..
            cd('ANNOTATIONS')
            imwrite(X,strcat('s',sceneN,'-t',num2str(i),'.jpg'),'JPEG'); %save image
            cd .. 
            cd('CENTROIDS')
            cd ..
            
%             CENTROIDS.imgsize = imgsize;
%             CENTROIDS.centroids = xy;
%             CENTROIDS.scene = sceneN;
%             save(strcat('CENTROIDS-',sceneN,'.mat'),'CENTROIDS'); %save CENTROIDS
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        

        
        
        
 
end
function exportLabels(~,~)
global ExportNameKey ExportName  imgsize displaytracking SceneList Tracked ImageDetails expDirPath trackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings adjuster

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(trackPath)
cd ..
CENTROID = struct();
    for scenenumber = 1:length(SceneList)
        adjuster=1;
        cd(trackPath)
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
            displaytracking =1;
            finalbutton_callback([],[]);
       xy = labelCells([],[]);
            set(gcf,'Units','Pixels');
            set(gca,'Units','Pixels');
            P = get(gca,'pos');
            F = getframe(gcf,P);
            [X,Map] = frame2im(F);
            set(gcf,'Units','normalized');
            set(gca,'Units','normalized');
            cd(trackPath)
            cd ..
            cd ..
            cd('ANNOTATIONS')
            imwrite(X,strcat(sceneN,'.jpg'),'JPEG'); %save image
            cd .. 
            cd('CENTROIDS')
            
            CENTROIDS.imgsize = imgsize;
            CENTROIDS.centroids = xy;
            CENTROIDS.scene = sceneN;
            save(strcat('CENTROIDS-',sceneN,'.mat'),'CENTROIDS'); %save CENTROIDS
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        

        
        
        
    end
    adjuster=0;
end


%% tracking functions
function displayTrackingButton_Callback(~,~)
global    displaytracking 

if displaytracking == 0
displaytracking = 1;
else
displaytracking =0;
end

setSceneAndTime
end

function Tracked = makeTrackingFile(timeFrames)

Tracked = cell(1,length(timeFrames));
Framez = struct();
Framez.Cellz =[];
Framez.imgsize=[];
Framez.filename=[];


for i = 1:length(Tracked)
    Tracked{i} = Framez;
end

end
function trackbutton_Callback(~,~)
global  Tracked ImageDetails TC segmentPath nucleus_seg
pvalue = ImageDetails.Scene;

    trackfilelist = {'yes','no'};
    [S,~] = listdlg('PromptString','Are you sure you want to run tracking?',...
                'SelectionMode','single',...
                'ListSize',[200 300],...
                'ListString',trackfilelist);
            
            if S==1
Tracked = FrickTrackCellsYeah(segmentPath,pvalue,nucleus_seg);
            else
            end

            TC =1;
setSceneAndTime;
end

function Tracked = loadTrackedStructure
global trackingPath timeFrames TC ExportName runIterate ImageDetails
if runIterate ==0
    cd(trackingPath)
    trackfile = dir(strcat('*',ImageDetails.Scene,'*.mat'));
    if ~isempty(trackfile)
        trackfilelist = {trackfile.name};
        Selection=[];
        [Selection,~] = listdlg('PromptString','Select a tracking file:',...
                    'SelectionMode','single',...
                    'ListSize',[500 300],...
                    'ListString',trackfilelist);
        if ~isempty(Selection)
        load(trackfilelist{Selection}); %load Tracked
        else
        Tracked = makeTrackingFile(timeFrames);
        end

        if isempty(Tracked{1}.Cellz)
        TC=0;
        else
        TC =1;
        end

    else
        Tracked = makeTrackingFile(timeFrames);
    end
else
    Tracked = makeTrackingFile(timeFrames);
end


end
function loadTrackingFile_callback(~,~)
global  Tracked TC

Tracked = loadTrackedStructure;
TC =1;
end

%make trajectories for overlay of tracking
function traject = trackingTrajectories(frameToLoad,ImageDetails,timeFrames)
global Tracked imgsize psettings


%   determine the frame to load
    t = timeFrames;
    framesThatMustBeTracked = psettings.framesThatMustBeTracked;

    xy = cell(1,t);
    lxy = zeros(1,t);
    CC = Tracked{framesThatMustBeTracked(1)}.Cellz;
    PX = CC.PixelIdxList;
    makeIMG = false(1,length(PX));

    for i = 1:t %determine centroids from 1:t for plotting a tracking tail 
        CC = Tracked{i}.Cellz;

        PXX =  CC.PixelIdxList;
        makeCentroids = find((~makeIMG)==1); %you don't want to do this calculation through all the NaN, so index for non NAN
        PX = PXX(~makeIMG);
        mx = nan(1,length(PXX));
        my = nan(1,length(PXX));
            for j = 1:length(PX)
            px = PX{j};
            y = rem(px-1,imgsize(1))+1; %these two lines replace ind2sub
            x = (px-y)/imgsize(2) + 1;  %these two lines replace ind2sub

            sx = sort(x);
            sy = sort(y);

            pseudomean = round(length(sx)./2);
                if pseudomean == 0
                mx(makeCentroids(j)) = NaN; 
                my(makeCentroids(j)) = NaN;
                else    
                mx(makeCentroids(j)) = sx(pseudomean);  %use the make Centroids index to keep the centroids the same color when plotting
                my(makeCentroids(j)) = sy(pseudomean);  
                end
            end

        xy{i} = horzcat(mx',my');
        lxy(i) = size(xy{i},1);
    end
    
    traject = nan(max(lxy),2,t);

    for i = 1:t
        traject(1:lxy(i),1:2,i) = xy{i};
    end
end
function trt = calculateTrackingLogical(Stacked)

    %make a matrix where 1 means yes there is a segmented cell and 0 means
    %there is no cell there (just NaN value)
    MAR1 = Stacked{1}.Cellz.PixelIdxList;
    trt = false(length(MAR1),length(Stacked));
    for i=1:length(Stacked)
        MAR = Stacked{i}.Cellz.PixelIdxList; %MAR contains the pixel indices for each nuclei
        logx = false(1,length(MAR));
        for j = 1:length(MAR)
           x = MAR{j}; 
           logx(j)  = ~isnan(x(1));
        end
        trt(:,i) = logx;
    end
end
function trt = calculateTrackingLogicalOLD(Stacked)

%make a matrix where 0 means yes there is a segmented cell and 1 means
%there is no cell there (just NaN value)

mar = cell(1,length(Stacked));
for i=1:length(Stacked)
MAR = Stacked{i}.Cellz.PixelIdxList;
    % mar{i} = cellfun(@(x) isnan(x(1)),MAR,'UniformOutput',1);
    % CELLFUN IS SLOW!!!!!
    logx = zeros(1,length(MAR));
    for j = 1:length(MAR)
       x = MAR{j}; 
       if isempty(x)
       x=NaN;    
       end
       logx(j)  = isnan(x(1));
    end
mar{i}=logx;
end

ml = cellfun(@length,mar);
trt = zeros(max(ml),length(Stacked));
for i=1:length(Stacked)
trt(1:ml(i),i) = mar{i};
end
end
function visualizeTrackedStructure(Tracked)

mar = cell(1,length(Tracked));
for i=1:length(Tracked)
MAR = Tracked{i}.Cellz.PixelIdxList;
mar{i} = cellfun(@(x) nansum(isnan(x)),MAR,'UniformOutput',1);
end
ml = cellfun(@length,mar);
trt = zeros(max(ml),length(Tracked));
for i=1:length(Tracked)
trt(1:ml(i),i) = mar{i};
end

figure(11)
strt = sum(trt,2);
[~,itrt] = sort(strt);
sortedtrt = trt(itrt,:);
startwithCELLidx = (sortedtrt(:,1) ==0);
endwithCELLidx = (sortedtrt(:,end) ==0);
subplot(2,1,1);plot(sortedtrt(startwithCELLidx,:)');
title('should only see rising lines')
ylim([-0.5 1.5])
subplot(2,1,2);plot(sortedtrt(endwithCELLidx,:)');
title('should only see falling lines')
ylim([-0.5 1.5])


end
function Trackedz = trackingCosmetics(Stacked)

%this identifies the maximum length to identify the total number of cells
%MARlength contains the number of cells in each frame
    MARlength = zeros(1,length(Stacked));
    for i=1:length(Stacked)
    MAR = Stacked{i}.Cellz.PixelIdxList;
    MARlength(i) = length(MAR);
    end

    %this makes a cell for each cell in all of the frames so each frame has
    %the same total number of cells. (When track is ended NaN is added)
    for i=1:length(Stacked)
    MAR = Stacked{i}.Cellz.PixelIdxList;
    PX = cell(1,max(MARlength));
    PX(1:MARlength(i)) =  MAR;
    PX(MARlength(i)+1:max(MARlength)) =  {NaN};
    Stacked{i}.Cellz.PixelIdxList = PX;
    end

    
trt = calculateTrackingLogical(Stacked);
idxo = trt;
%add the number of NaN remaining
%the goal is that each track will only have one continous segment tracked
didxo = diff(idxo,[],2);
    for j=1:size(didxo,1) %iterate through each cell
        endoftrack = find(didxo(j,:) == -1);
        beginoftrack = find(didxo(j,:) == 1)+1;
        
        %add frame 1 if track exists at frame 1
        if idxo(j,1)==1
            beginoftrack = [1 beginoftrack];
        end
        
        %add last frame if track exists at last frame
        if idxo(j,size(idxo,2)) == 1
            endoftrack = [endoftrack size(idxo,2)];
        end
        
        if length(beginoftrack)==length(endoftrack)
        else
            error(('mis-identification of track length?'))
        end
        
        tracks = zeros(length(beginoftrack),2);
        for i =1:length(beginoftrack)
            tracks(i,:)  = [beginoftrack(i) endoftrack(i)];
        end

            %move the post-gap cells to the end
            for trackidx = 2:size(tracks,1) %start from the second track
                beginoftrack = tracks(trackidx,1);
                endoftrack = tracks(trackidx,2);
                for frame = 1:length(Stacked)
                    PX = Stacked{frame}.Cellz.PixelIdxList;
                    px = cell(1,length(PX)+1);
                    px(1:length(PX)) = PX;
                    if frame>beginoftrack-1 && frame<endoftrack+1 % if frame is where cell is tracked
                        px(length(PX)+1) = PX(j); %add the track to the end
                        px(j) = {NaN}; %remove the track data from where it was previously located
                    else %if the frame is not where the cell is located fill the track with NaN
                        px(length(PX)+1) = {NaN};
                    end
                    Stacked{frame}.Cellz.PixelIdxList = px;
                end
                
                %make the same changes to trt
                Trt = false(size(idxo,1)+1,size(idxo,2));
                Trt(1:size(idxo,1),:) = idxo; %set up new matrix
                Trt(j,beginoftrack:endoftrack) = false; %make the track at the current position false
                Trt(size(idxo,1)+1,beginoftrack:endoftrack) = true; %make the track true at the new ending position
                idxo = Trt;
            end
    end



%remove tracks that have no cells
    index = find(idxo(:,1)==0); %find cells in the first frame that are not tracked
    strt = sum(~idxo(index,:),2); %determine how 
    istrt = (strt == length(1:length(Stacked)));
    pxidxremove = index(istrt);
    pxidx = 1:length(Stacked{1}.Cellz.PixelIdxList);
    pxidx(pxidxremove) = [];
    %update the fields of Stacked.Cellz to be accurate
        for j=1:length(Stacked)
        PX = Stacked{j}.Cellz.PixelIdxList;
        px = PX(pxidx);
        Stacked{j}.Cellz.PixelIdxList = px;
        Stacked{j}.Cellz.NumObjects = length(px);
        CC.Connectivity = Stacked{1}.Cellz.Connectivity;
        CC.ImageSize = Stacked{1}.Cellz.ImageSize;
        CC.NumObjects = Stacked{1}.Cellz.NumObjects;
        CC.PixelIdxList = px;
        stats = regionprops(CC,'Centroid');
        centroidarray = {stats.Centroid};
        Centroid = zeros(2,length(centroidarray));
        for ijk = 1:length(centroidarray)
            Centroid(:,ijk) = centroidarray{ijk};
        end
        CC.Centroid = Centroid;
        Stacked{j}.Cellz = CC;
        end
    Trackedz=Stacked;

    stophere=1;
end
function Trackedz = trackingCosmeticsOld(Stacked)


%this identifies the maximum length to identify the total number of cells
    MARlength = zeros(1,length(Stacked));
    for i=1:length(Stacked)
    MAR = Stacked{i}.Cellz.PixelIdxList;
    MARlength(i) = length(MAR);
    end

    %this makes a cell for each cell in all of the frames
    for i=1:length(Stacked)
    MAR = Stacked{i}.Cellz.PixelIdxList;
    PX = cell(1,max(MARlength));
    PX(1:MARlength(i)) =  MAR;
    PX(MARlength(i)+1:max(MARlength)) =  {NaN};
    Stacked{i}.Cellz.PixelIdxList = PX;
    end

    %this makes a cell for each cell in all of the frames
    for i=1:length(Stacked)
    MAR = Stacked{i}.Cellz.PixelIdxList;
%     llist = Stacked{i}.Cellz.PixelIdxList.lengthlist;
    MARi = cellfun(@(x) length(x),MAR,'UniformOutput',1);
    MARidx = ~logical(MARi<2);
    MARidx = ~logical(cellfun(@(x) length(x)<2,MAR,'UniformOutput',1));
    PX = cell(1,max(MARlength));
    PX(MARidx) =  MAR(MARidx);
    PX(~MARidx) =  {NaN};
    Stacked{i}.Cellz.PixelIdxList = PX;
    end


trt = calculateTrackingLogical(Stacked);
%add the number of NaN remaining
idxo = ~logical(trt);  %set value of 1 where trt = 0;
didxo = diff(idxo,[],2);
    for j=1:size(didxo,1) 
    %if you ever find a negative 1, you can never have a positive 1 afterward
    beginoftrack = find(idxo(j,:)==1,1,'first');%first frame of first tracked cell
    endoftrack = find(didxo(j,:)==-1,1,'last');%last frame of last tracked cell
    fni = find(didxo(j,:)== -1); %frame at which tracking ends
    fi = find(didxo(j,:)==1); %frame where tracking begins again

    fi = fi(fi>beginoftrack);
%     fni = fni(fni==endoftrack | fni<endoftrack); %original code
    fni = fni(fni<endoftrack);

    if j==120
        stophere=1;
    end
    
        if ~(sum([isempty(fi) isempty(fni)])>0) %if a track has one cell then none then another cell (track with gap)
        idxtest = fi>fni; 
        fiToRemove = fi(idxtest)+1; %the extra cell begins at (fi+1)
            %move the post-gap cells to the end
            for jy = 1:length(fiToRemove)
                for ky = 1:length(Stacked)
                    PX = Stacked{ky}.Cellz.PixelIdxList;
                    px = cell(1,length(PX)+1);
                    px(1:length(PX)) = PX;
                    if ky<fiToRemove(jy)  %if it is before the extra cell, then set the end values to NaN and do nothing to the tracking
                    %PX modifications
                    px(length(PX)+1) = {NaN};
                    Stacked{ky}.Cellz.PixelIdxList = px;                    
                    else %once the frame corresponds to the beginning of the extra cell, move it to the end and change the orignial tracking to NaN;
                    px(length(PX)+1) = PX(j);
                    px(j) = {NaN};
                    Stacked{ky}.Cellz.PixelIdxList = px;
                    end
                end
            end
            %make the same changes to trt
                for jy = 1:length(fiToRemove)
                Trt = zeros(size(trt,1)+1,size(trt,2));
                Trt(1:size(trt,1),:) = trt; %set up new matrix
                Trt(size(trt,1)+1,1:fiToRemove(jy)-1) = 1; %for all times before gap cell, make = NaN
                Trt(size(trt,1)+1,fiToRemove(jy):end) = trt(j,fiToRemove(jy):end); %for all times following, make equal to gap cell
                Trt(j,fiToRemove(jy):end) = 1;
                trt = Trt;
                end
            
        else    
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% trtz = calculateTrackingLogical(Stacked); %1 is no cell, 0 is cell

%if strt = 1 all the way through = length(1:length(Tracked))
index = find(trt(:,1)==1);
strt = sum(trt(index,:),2);
    istrt = (strt == length(1:length(Stacked)));
    pxidxremove = index(istrt);
    pxidx = 1:length(Stacked{1}.Cellz.PixelIdxList);
    pxidx(pxidxremove) = [];
%     pxidx(pxidxremove) = {NaN};
        for j=1:length(Stacked)
        PX = Stacked{j}.Cellz.PixelIdxList;
        px = PX(pxidx);
        Stacked{j}.Cellz.PixelIdxList = px;
        Stacked{j}.Cellz.NumObjects = length(px);
        end
    Trackedz=Stacked;
stophere=1;
end

%function for tracking cells
function [ Tracked ] = FrickTrackCellsYeah(segmentPath,pvalue,nucleus_seg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% cd('D:\Users\zeiss\Documents\MATLAB')
minimum_nucleus_area=20;
% minimum_nucleus_area=35;
Tracked = [];
Frame = struct();


       %load segmented image
        cd(segmentPath)
        ff = dir(strcat('*',pvalue,'*',nucleus_seg,'*'));
        if length(ff)>1
            ff = dir(strcat('*',pvalue,'*',nucleus_seg,'*nucleus*')); 
        end
%         ff = dir(strcat(ImageDetails.Channel,'*'));
        filename = char(ff.name);
        segmentfileObject = matfile(filename);
        segmentimgstack = segmentfileObject.IfFinal;
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section involves calculating the nearest neighbor to the centroid
%and organizing PixelLists and Centroid lists to match nearest neighbor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
for i = 1:size(segmentimgstack,3)
    img = segmentimgstack(:,:,i);
    
%     %%%incorporate mKate signal into the segmentation
%     cd .. 
%     cd('mKatebinary_flat')
%         imgmkate = imread(char(cfilemkate{i}));
%     If = uint8(logical(img) | logical(imgmkate));
%     img = If;
    
CC = bwconncomp(img);
PX = CC.PixelIdxList;
%choose segmented nuclei above a certain size only
segment_Area = cellfun(@length,PX);
true_cells = segment_Area > minimum_nucleus_area;
PX = PX(true_cells);
CC.PixelIdxList = PX;
CC.NumObjects = numel(PX);

PXarray{i} = PX;

if i==1
S = regionprops(CC,'Centroid');
Smat = vertcat(S.Centroid);

B{i}=Smat;
Frame.filename = filename;
Frame.Cellz = CC;
Tracked{i} = Frame;
else
S = regionprops(CC,'Centroid');
Smat = vertcat(S.Centroid); 

if ~isempty(Smat)
[Idx,Eps] = knnsearch(Smat,B{i-1},'K',1); %B{i-1} = Smat(Idx)
SameCellPX = PX(Idx);
else
SameCellPX = Tracked{i-1}.Cellz.PixelIdxList;
end


num_cells_set = 1:size(Smat,1);
[n, bin] = histc(Idx, num_cells_set);
multiple = find(n>1); %the same cell is called closest to two previous cells
missers = find(n<1); %these are likely new cells
    if ~isempty(multiple)
        for loop = multiple'
            index    = find(ismember(bin, loop));
            winnerofrepeatidx =  find(Eps == min(Eps(index)));
            loserz = setdiff(index,winnerofrepeatidx);
%             if length(loserz)>100
%                 stophere=1;
%             end
%             for lsrz = loserz'
%             SameCellPX(end+1) = SameCellPX(lsrz);
%             end
            
            SameCellPX(loserz) = {NaN}; %remove multiple links to same cell from previous frame so that cell is only linked to one previous
            
        end
    end


AllCellsPX = horzcat(SameCellPX,PX(missers));
CC.PixelIdxList = AllCellsPX;
CC.NumObjects = numel(AllCellsPX);
S = regionprops(CC,'Centroid');
Smat = vertcat(S.Centroid);
CC.Centroid = Smat;
B{i}=Smat;

Frame.filename = filename;
Frame.Cellz = CC;
Tracked{i} = Frame;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section involves calculating the nearest neighbor to the centroid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stophere=1;




% for i = 1:length(B)
% xy = B{i};
% c = colormap(lines(length(xy)));
% scatter(xy(:,1),xy(:,2),[],c);hold on
% end

stophere=1;
end

%% Image Display functions
function setSceneAndTime
global TC DICimgstack dfoName  nucleus_seg backgroundimgstack bfoName nfoName background_seg cell_seg nucleusimgstack sfoName segmentimgstack  channelimgstack cfoName segmentPath frameToLoad ImageDetails  Tracked SceneList  trackPath imgfile mstackPath

%check for empty variables
    cd(mstackPath)
    if isempty(ImageDetails.Scene)
        ImageDetails.Scene = SceneList{1};
    end
    if isempty(ImageDetails.Channel)
        ImageDetails.Channel = nucleus_seg;
    end
    
%determine the frame to load
    if isempty(ImageDetails.Frame)
       ImageDetails.Frame = frameToLoad;
       t = ImageDetails.Frame;
    else
       t = ImageDetails.Frame;
    end


%choose the channel image
        cd(mstackPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',cell_seg,'*'));
        filename = char(ff.name);
        if ~isempty(cfoName) %if channelfileObject has been made, check to see if the scene has changed. 
            [a,~] = regexp(cfoName,ImageDetails.Scene);
            if isempty(a) %if the scene has changed load the new channelimgstack
                 channelfileObject = matfile(filename);
                 channelimgstack = channelfileObject.flatstack;
                 cfoName = char(channelfileObject.Properties.Source);%update cfoName
            elseif ~isempty(a) && isempty(channelimgstack)  %if the scene is same but unloaded
                 channelfileObject = matfile(filename);
                 channelimgstack = channelfileObject.flatstack;
                 cfoName = char(channelfileObject.Properties.Source);%update cfoName
            else
            end
        else %if no cfoName, then 
                 channelfileObject = matfile(filename);
                 channelimgstack = channelfileObject.flatstack;
                 cfoName = char(channelfileObject.Properties.Source);%update cfoName
        end
        
        
        cellImg = channelimgstack(:,:,t);
        

        cd(mstackPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',nucleus_seg,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
        filename = char(ff.name);
        if ~isempty(nfoName) %if channelfileObject has been made, check to see if the scene has changed. 
            [a,~] = regexp(nfoName,ImageDetails.Scene);
            if isempty(a) %if the scene has changed load the new channelimgstack
                 nucleusfileObject = matfile(filename);
                 nucleusimgstack = nucleusfileObject.flatstack;
                 nfoName = char(nucleusfileObject.Properties.Source);%update cfoName
%                  disp('if -> if')
            elseif ~isempty(a) && isempty(nucleusimgstack)  %if the scene is same but unloaded
                 nucleusfileObject = matfile(filename);
                 nucleusimgstack = nucleusfileObject.flatstack;
                 nfoName = char(nucleusfileObject.Properties.Source);%update cfoName
%                  disp('if -> elseif')
            else
%                 disp('if -> else')
                %dont do anything
            end
        else %if no cfoName, then 
                 nucleusfileObject = matfile(filename);
                 nucleusimgstack = nucleusfileObject.flatstack;
                 nfoName = char(nucleusfileObject.Properties.Source);%update cfoName
%                  disp('else')
        end
        nucleusImg = nucleusimgstack(:,:,t);
        
        
        
         cd(mstackPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*','DIC','*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
        filename = char(ff.name);
        if ~isempty(dfoName) %if channelfileObject has been made, check to see if the scene has changed. 
            [a,~] = regexp(dfoName,ImageDetails.Scene);
            if isempty(a) %if the scene has changed load the new channelimgstack
                 DICfileObject = matfile(filename);
                 DICimgstack = DICfileObject.flatstack;
                 dfoName = char(DICfileObject.Properties.Source);%update cfoName
%                  disp('if -> if')
            elseif ~isempty(a) && isempty(DICimgstack)  %if the scene is same but unloaded
                 DICfileObject = matfile(filename);
                 DICimgstack = DICfileObject.flatstack;
                 dfoName = char(DICfileObject.Properties.Source);%update cfoName
%                  disp('if -> elseif')
            else
%                 disp('if -> else')
                %dont do anything
            end
        else %if no cfoName, then 
            flist = dir(filename);
            if ~isempty(flist)
                 DICfileObject = matfile(filename);
                 DICimgstack = DICfileObject.flatstack;
                 dfoName = char(DICfileObject.Properties.Source);%update cfoName
%                  disp('else')
            else
                DICimgstack = false(size(nucleusimgstack));
                dfoName = filename;
            end
        end
        DICImg = DICimgstack(:,:,t);
        
         
       %load nucleus segmented image
        cd(segmentPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',nucleus_seg,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
        if length(ff)>1
           ff = dir(strcat('*',ImageDetails.Scene,'*',nucleus_seg,'*nucleus*'));
        end
        filename = char(ff.name);
        
        if ~isempty(sfoName) %if channelfileObject has been made, check to see if the scene has changed. 
            [a,~] = regexp(sfoName,ImageDetails.Scene);
            if isempty(a) %if the scene has changed load the new channelimgstack
                 segmentfileObject = matfile(filename);
                 segmentimgstack = segmentfileObject.IfFinal;
                 sfoName = char(segmentfileObject.Properties.Source);%update cfoName
%                  disp('if -> if')
            elseif ~isempty(a) && isempty(segmentimgstack)  %if the scene is same but unloaded
                 segmentfileObject = matfile(filename);
                 segmentimgstack = segmentfileObject.IfFinal;
                 sfoName = char(segmentfileObject.Properties.Source);%update cfoName
%                  disp('if -> elseif')
            else
%                 disp('if -> else')
                %dont do anything
            end
        else %if no cfoName, then 
                 segmentfileObject = matfile(filename);
                 segmentimgstack = segmentfileObject.IfFinal;
                 sfoName = char(segmentfileObject.Properties.Source);%update cfoName
%                  disp('else')
        end
        segmentimg = segmentimgstack(:,:,t);
        
        
        
        %load background segmented image
        cd(segmentPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',background_seg,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
        if length(ff)>1
            ff = dir(strcat('*',ImageDetails.Scene,'*',background_seg,'*background*'));
        end
        filename = char(ff.name);
        if ~isempty(bfoName) %if channelfileObject has been made, check to see if the scene has changed. 
            [a,~] = regexp(bfoName,ImageDetails.Scene);
            if isempty(a) %if the scene has changed load the new channelimgstack
                 backgroundfileObject = matfile(filename);
                 backgroundimgstack = backgroundfileObject.IfFinal;
                 bfoName = char(backgroundfileObject.Properties.Source);%update cfoName
%                  disp('if -> if')
            elseif ~isempty(a) && isempty(backgroundimgstack)  %if the scene is same but unloaded
                 backgroundfileObject = matfile(filename);
                 backgroundimgstack = backgroundfileObject.IfFinal;
                 bfoName = char(backgroundfileObject.Properties.Source);%update cfoName
%                  disp('if -> elseif')
            else
%                 disp('if -> else')
                %dont do anything
            end
        else %if no cfoName, then 
                 backgroundfileObject = matfile(filename);
                 backgroundimgstack = backgroundfileObject.IfFinal;
                 bfoName = char(backgroundfileObject.Properties.Source);%update cfoName
%                  disp('else')
        end
        backgroundimg = backgroundimgstack(:,:,t);
         
        
        
        
        
        
        
        
        
        
%         channelimg = double(loadUpTiffStackFrame(char(ff.name),t));
    %     channelimg = double(imread(char(imgfile.name)));    %load normal image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   choose the segmentation image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(Tracked{1}.Cellz) %Tracked is empty try to load tracking or make new tracking structure
        TC = 0;
        cd(segmentPath)
        imgfile = dir(strcat('*',ImageDetails.Scene,'*',nucleus_seg,'*'));
        If = segmentimg;
        Tracked = loadTrackedStructure;
    else  %if there exists segmenttracking already...then load that. 
        CC = Tracked{t}.Cellz;
        PX = CC.PixelIdxList;
    %     makeIMG = cellfun(@(x) length(x)==1,PX,'UniformOutput',1); %choose only the cells without NAN
        makeIMG = cellfun(@(x) length(x)<2,PX,'UniformOutput',1); %choose only the cells without NAN
        CC.PixelIdxList = PX(~makeIMG);
        CC.NumObjects = length(PX(~makeIMG));

        segmentimgL = labelmatrix(CC);
        segmentimgz = false(size(segmentimgL));
        segmentimgz(segmentimgL>0)=1;
        If = segmentimgz;
    end
    




% disp(TC)
if TC == 1
Tracked = trackingCosmetics(Tracked);
end
TC=0;


if strcmp(ImageDetails.Channel,nucleus_seg)
    channelimg = nucleusImg;
elseif strcmp(ImageDetails.Channel,cell_seg)
    channelimg = cellImg;
elseif strcmp(ImageDetails.Channel,'EGFP')
    channelimg = cellImg;
elseif strcmp(ImageDetails.Channel,'DIC')
    channelimg = DICImg;
elseif strcmp(ImageDetails.Channel,'BKGbinary')
     channelimg = cellImg;
     backgroundimg(segmentimg) = true; 
     prim = imdilate(bwperim(~logical(backgroundimg)),strel('square',1));
     channelimg(prim) = max(max(channelimg));
elseif strcmp(ImageDetails.Channel,'reporter_quantify')
elseif strcmp(ImageDetails.Channel,'overlay')
    channelimg = zeros(size(cellImg,1),size(cellImg,2),3);
    channelimg(:,:,1) = nucleusImg;
    channelimg(:,:,2) = cellImg;
    channelimg(:,:,3) = DICImg;
end

backgroundimg(segmentimg) = true; 
bkgpixels = channelimg(backgroundimg);
bkgmedian = nanmedian(bkgpixels);
displayImageFunct(If,channelimg,bkgmedian);

end


function contrast_Callback(~,~)
global tcontrast lcontrast adjuster
prompt = {'High Contrast Percent Limit','Low Contrast Percent Limit'};
dlg_title = 'Contrast limits from 0 to 100';
num_lines = 1;
defaultans = {num2str(tcontrast),num2str(lcontrast)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
tcontrast = str2double(answer{1});
lcontrast = str2double(answer{2});

adjuster =1;
setSceneAndTime
adjuster =0;
end

function contrast_Callbacknew(~,~)
global tcontrast lcontrast adjuster
prompt = {'High Contrast Limit','Low Contrast Limit'};
dlg_title = 'Contrast limits from 0 to 100';
num_lines = 1;
defaultans = {num2str(tcontrast),num2str(lcontrast)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
tcontrast = str2double(answer{1});
lcontrast = str2double(answer{2});

adjuster =1;
setSceneAndTime
adjuster =0;
end

function displayImageFunct(If,channelimg,bkgmedian)
global displaycomments psettings plottingON trunccmaplz timeFrames lprcntlt prcntlt tcontrast lcontrast MainAxes displaytracking ImageDetails frameToLoad prcntlz lprcntlz prcntlk lprcntlk prcntl lprcntl D ExpDate cmap cmaplz adjuster


%determine current time Frame
    t = frameToLoad;

%delete old images to keep memory usage low
    axes(MainAxes);
    children = findobj(MainAxes,'Type','image');
    delete(children);

% important for updating the contrast
    %update contrast if time =1
    ifCHANGEofCHANNELorSCENE=0;
    if t==1
    D='new';
    ifCHANGEofCHANNELorSCENE=1;
    end

%update contrast if channel has changed
    if ~strcmp(ImageDetails.Channel,D)
    ifCHANGEofCHANNELorSCENE=1;
    D = ImageDetails.Channel;
    end

%update contrast if contrast values are updated
    if adjuster ==1
        ifCHANGEofCHANNELorSCENE = 1;
        D = ImageDetails.Channel;
    end

%scripts for displaying contrasted image
    if strcmp(ImageDetails.Channel,'overlay') %when overlay display is desired
        disprgb = zeros(size(channelimg));
        channelimgrgb = channelimg;
        for i = 1:size(channelimg,3)
            if i==3
                channelimg = channelimgrgb(:,:,i);
                lprcntl = prctile(channelimg(:),lcontrast);
                prcntl = prctile(channelimg(:),tcontrast);
                scaleFactor = 150./(prcntl - lprcntl);
                dispimg = channelimg.*scaleFactor;
                dispimg = dispimg-(lprcntl.*scaleFactor);
                dispimg(dispimg> 255) =254;
                dispimg(dispimg<0) = 0;
            else
                channelimg = channelimgrgb(:,:,i);
                lprcntl = prctile(channelimg(:),lcontrast);
                prcntl = prctile(channelimg(:),tcontrast);
                scaleFactor = 255./(prcntl - lprcntl);
                dispimg = channelimg.*scaleFactor;
                dispimg = dispimg-(lprcntl.*scaleFactor);
                dispimg(dispimg> 255) =254;
                dispimg(dispimg<0) = 0;
            end
            
            if i==1
                If = bwperim(If) | bwperim(imdilate(If,strel('disk',1)));
                dispimg(If>0)=255;  
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
            colormap(cmap);
        end
        
        dispimg = disprgb;
        
    else  %under normal circumstances
        if ifCHANGEofCHANNELorSCENE==1
            lprcntl = prctile(channelimg(:),lcontrast);
            prcntl = prctile(channelimg(:),tcontrast);
            scaleFactor = 255./(prcntl - lprcntl);
            ifCHANGEofCHANNELorSCENE=0;
        end
        lprcntl = bkgmedian.*0.90;
        scaleFactor = 255./(prcntl - lprcntl);
        dispimg = channelimg.*scaleFactor;
        dispimg = dispimg-(lprcntl.*scaleFactor);
        dispimg(dispimg> 255) =254;
        colormap(cmap);
        If = bwperim(If) | bwperim(imdilate(If,strel('disk',1)));
        dispimg(If>0)=255;
    end

%title the displayed image
    himg = imagesc(uint8(dispimg));
    himgax = get(himg,'Parent');
    himgax.CLim = [0 256];
    ttl = get(himgax,'Title');
    t = ImageDetails.Frame;
    set(ttl,'String',[ExpDate ' ' ImageDetails.Scene ' frame ' num2str(t) ' out of ' num2str(timeFrames)]);
    set(ttl,'FontSize',12);
    
    if plottingON == 0
        psettings = PlotSettings_callback([],[]);
        plottingON=1;
    end
    framesThatMustBeTracked = psettings.framesThatMustBeTracked;

    if ~(t==1)
        if displaytracking==1
            trajectForPlot = trackingTrajectories(t,ImageDetails,timeFrames);
            himgax.NextPlot = 'add';
            mainplotX = squeeze(trajectForPlot(:,1,:)); %28x22 means 28 cells on frame 22;
            mainplotY = squeeze(trajectForPlot(:,2,:));

            if size(mainplotY,2) == 1
                mainplotY=mainplotY';
                mainplotX=mainplotX';
            end
            %only plot if the cell is currently tracked/segmented in this frame
                idx = ~isnan(mainplotY(:,t));
%                 h = plot(mainplotX(idx,1:t)',mainplotY(idx,1:t)','LineWidth',1,'Marker','s','MarkerSize',8);
                h = plot(mainplotX(idx,1:t)',mainplotY(idx,1:t)','LineWidth',2);

            %generate colormap based on number of cells tracked
                cnew=[];
%                 ccc = vertcat(colormap('summer'),colormap('autumn'),colormap('winter'),colormap('spring'));
%                 ccc = vertcat(colormap('hsv'),colormap('hot'));
                ccc = colormap('colorcube');
                cccyc = 0;
                for k = 1:size(ccc,1)
                    cvec = ccc(k,:);
                    if sum(cvec)>0.5 && sum(abs(diff(cvec)))>0.2 && sum(cvec)<2
                        cccyc = cccyc+1;
                        cnew(cccyc,:) = cvec;
                    end
                end
                ccnew = zeros(size(mainplotX,1),size(cnew,2));
                for j = 1:size(cnew,2)
                    x = linspace(0,1,size(cnew,1));
                    v = cnew(:,j);
                    xq = linspace(0,1,size(mainplotX,1));
                    ccnew(1:length(xq),j) = interp1(x,v,xq);
                end
                
                cmaplz = ccnew;
                cmapl = cmaplz;
                idxa = find(idx==1);
                trunccmaplz = cmaplz;
                
                plotcmap = zeros(length(idxa),3);
                for i=1:length(h)
%                     h(i).Color = cmapl(idxa(i),:);
                    plotcmap(i,:) =  cmapl(idxa(i),:);
%                     h(i).MarkerFaceColor = cmapl(idxa(i),:);
%                     h(i).MarkerEdgeColor = cmapl(idxa(i),:)./1.2;
                end
                set(h, {'color'}, num2cell(plotcmap,2));
                colormap(cmap);%return colormap so images display properly
                hax = h.Parent;
                hax.Color = 'none';
                himgax.CLim = [0 256];
                himgax.NextPlot = 'replace';
        end
    end
    himgax.YTick = [];
    himgax.XTick = [];
end

function displayImageFunctold(If,channelimg)
global displaycomments timeFrames lprcntlt prcntlt tcontrast lcontrast MainAxes displaytracking ImageDetails frameToLoad prcntlz lprcntlz prcntlk lprcntlk prcntl lprcntl D ExpDate cmap cmaplz adjuster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(frameToLoad,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   display the images overlayed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(MainAxes);
children = findobj(MainAxes,'Type','image');
delete(children);

ifCHANGEofCHANNELorSCENE=0;
if t==1
D='new';
ifCHANGEofCHANNELorSCENE=1;
end

if ~strcmp(ImageDetails.Channel,D)
ifCHANGEofCHANNELorSCENE=1;
D = ImageDetails.Channel;
end

if adjuster ==1
    ifCHANGEofCHANNELorSCENE = 1;
    D = ImageDetails.Channel;
end

if strcmp(ImageDetails.Channel,'overlay') %when overlay display is desired
    imgone = channelimg(:,:,1);
    imgtwo = channelimg(:,:,2);
    imgthree = channelimg(:,:,3);
    imgfour = channelimg(:,:,4);
    channelimg = zeros(size(channelimg,1),size(channelimg,2),3);
        if ifCHANGEofCHANNELorSCENE==1
        
        cimgline = reshape(imgone,[1 size(imgone,1).*size(imgone,2)]);
        lprcntlz = prctile(cimgline,lcontrast);
%         prcntlz = prctile(cimgline,tcontrast) - lprcntlz;
        prcntlz = prctile(cimgline-lprcntlz,tcontrast);
        
        cimgline = reshape(imgtwo,[1 size(imgtwo,1).*size(imgtwo,2)]);
        lprcntl = prctile(cimgline,lcontrast);
%         prcntl = prctile(cimgline,tcontrast)-lprcntl;
        prcntl = prctile(cimgline-lprcntl,tcontrast);
        
        cimgline = reshape(imgthree,[1 size(imgtwo,1).*size(imgtwo,2)]);
        lprcntlk = prctile(cimgline,lcontrast);
%         prcntlk = prctile(cimgline,tcontrast)-lprcntlk;
        prcntlk = prctile(cimgline-lprcntlk,tcontrast);
        
        cimgline = reshape(imgfour,[1 size(imgtwo,1).*size(imgtwo,2)]);
        lprcntlt = prctile(cimgline,lcontrast);
%         prcntlk = prctile(cimgline,tcontrast)-lprcntlk;
        prcntlt = prctile(cimgline-lprcntlt,tcontrast);
        ifCHANGEofCHANNELorSCENE=0;
        end
    imgone = uint8(((imgone-lprcntlz)./prcntlz).*255);
    imgtwo = uint8(((imgtwo-lprcntl)./prcntl).*255);
    imgthree = uint8(((imgthree-lprcntlk)./prcntlk).*255);
    imgfour = uint8(((imgfour-lprcntlt)./prcntlt).*255);
    
        imgone(imgone<imgfour) = imgfour(imgone<imgfour);
        imgtwo(imgtwo<imgfour) = imgfour(imgtwo<imgfour);
        imgthree(imgthree<imgfour) = imgfour(imgthree<imgfour);
    
    channelimg = uint8(channelimg);
    channelimg(:,:,2) = imgone;
    If = bwperim(If);
    imgtwo(If) = 255;
    channelimg(:,:,3) = imgtwo;
    imgthree(If) = 255;
    channelimg(:,:,1) = imgthree;
    
%     If = bwperim(If);
%     channelimg(:,:,3)=uint8(If.*255);

else  %under normal circumstances
        if ifCHANGEofCHANNELorSCENE==1
        cimgline = reshape(channelimg,[1 size(channelimg,1).*size(channelimg,2)]);
        lprcntl = prctile(cimgline,lcontrast);
%         prcntl = prctile(cimgline,tcontrast)-lprcntl;
        prcntl = prctile(cimgline-lprcntl,tcontrast);
        ifCHANGEofCHANNELorSCENE=0;
        end

    channelimg = uint8(((channelimg-lprcntl)./prcntl).*255);
    channelimg(channelimg == 255) =254;
    colormap(cmap);
    If = bwperim(If) | bwperim(imdilate(If,strel('disk',1)));
%     If = imdilate(If,strel('disk',1));
    channelimg(If>0)=255;
end


himg = imagesc(channelimg);
himgax = get(himg,'Parent');
himgax.CLim = [0 256];
% SLOW 
%himgax.Title.String = strcat(ExpDate,'...',ImageDetails.Scene,'...frame ',num2str(t),' out of', num2str(length(frameToLoad)));
%himgax.Title.FontSize = 12;
ttl = get(himgax,'Title');
t = ImageDetails.Frame;
set(ttl,'String',strcat(ExpDate,'...',ImageDetails.Scene,'...frame ',num2str(t),' out of', num2str(timeFrames)));
set(ttl,'FontSize',12);


    if ~(t==1)
        if displaytracking==1
            traject = trackingTrajectories(frameToLoad,ImageDetails);
            
            himgax.NextPlot = 'add';
            % rgbhax.NextPlot = 'replace';
            mainX = squeeze(traject(:,1,:));
            mainY = squeeze(traject(:,2,:));
            
            %only plot if the cell is currently tracked/segmented in this frame
            idx = ~isnan(mainY(:,t));
            h = plot(mainX(idx,:)',mainY(idx,:)','LineWidth',3);
            
            cmaplz = colormap(colorcube(size(mainX,1).*1.5));
            cmapl = cmaplz;
            idxa = find(idx==1);
                for i=1:length(h)
                    h(i).Color = cmapl(idxa(i),:);
                end
            colormap(cmap);
            hax = h.Parent;
            hax.Color = 'none';
            himgax.CLim = [0 256];
            himgax.NextPlot = 'replace';
        end
    end
himgax.YTick = [];
himgax.XTick = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saveChannelFiveImages
end

%% saving functions
function saveTrackingFileAs_callback(~,~)
global  trackingPath Tracked ExportName ImageDetails
cd(trackingPath)
prompt = 'filename of tracking structure to be saved?';
dlg_title = 'save tracking structure as...specific filename';
filename = char(inputdlg(prompt,dlg_title));
save(strcat(filename,'_',ImageDetails.Scene,'_',ExportName,'.mat'),'Tracked')
end

function trackSaveIterate_callback(~,~)
global runIterate SceneList ImageDetails TC expDirPath frameToLoad Tracked trackingPath ExportName timeFrames segmentPath nucleus_seg

runIterate =1;
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
            
            Tracked = FrickTrackCellsYeah(segmentPath,pvalue,nucleus_seg);
%             Tracked = FrickTrackCellsYeah(expDirPath,frameToLoad,pvalue,[]);
            TC =1;
            setSceneAndTime;


        %save
            % saveTrackingFileAs_callback([],[])
            cd(trackingPath)
            filename = 'initial';
            save(strcat(filename,'_',ImageDetails.Scene,'_',ExportName,'.mat'),'Tracked')
%             save(strcat(filename,ExportName,'.mat'),'Tracked')
            
    end
runIterate =0;

end


function trackSaveIterateChosen_callback(~,~)
global plottingON psettings runIterate SceneList ImageDetails TC expDirPath frameToLoad Tracked trackingPath ExportName timeFrames segmentPath nucleus_seg

runIterate =1;
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
            
            Tracked = FrickTrackCellsYeah(segmentPath,pvalue,nucleus_seg);
%             Tracked = FrickTrackCellsYeah(expDirPath,frameToLoad,pvalue,[]);
            TC =1;
            setSceneAndTime;
            
    if plottingON == 0
        psettings = PlotSettings_callback([],[]);
        plottingON=1;
    end
    
        %run chosen at two specific frames
        framesThatMustBeTracked = psettings.framesThatMustBeTracked;
        framesThatMustBeTracked(2) =  framesThatMustBeTracked(2)+10;
        for ftmbt = framesThatMustBeTracked
            ImageDetails.Frame =ftmbt;
            chosenOnesAllOnFrame_Callback
        end


        %save
            % saveTrackingFileAs_callback([],[])
            cd(trackingPath)
            filename = 'tsichosen';
            save(strcat(filename,'_',ImageDetails.Scene,'_',ExportName,'.mat'),'Tracked')
%             save(strcat(filename,ExportName,'.mat'),'Tracked')
            
    end
runIterate =0;

end

function saveTrackingFileAs_callbackJ(~,~,SceneDirPath)
global   Tracked
cd(SceneDirPath)
% prompt = 'filename of tracking structure to be saved?';
% dlg_title = 'save tracking structure as...specific filename';
% filename = char(inputdlg(prompt,dlg_title));
disp(iixixixi)
%why am i saving a specific file name without asking if I want to
%overwrite?
save(strcat('finalfricktrack.mat'),'Tracked')
end
function saveTrackingFile
global  trackPath
cd(trackPath)
save('generalTrackingSavedfricktrack.mat','Tracked')
end
function saveChannelFiveImages
global If imgfile trackPath channelinputs
cd(trackPath)
fname = imgfile.name;
[a,b] = regexp(fname,channelinputs);
fname(a:b) = 'c5';

cd('c5_flat')
imwrite(If,fname,'Tiff');
cd ..

end

%% functions for determining variables
function img_stack = loadImageStack(fname)
info = imfinfo(fname);
num_images = numel(info);
img_stack = zeros([info(1).Width info(1).Height num_images]);
for k = 1:num_images
   img_stack(:,:,k) = imread('NucleusBinary_flat.tif', k, 'Info', info);
end
end

function [timeFrames,frameToLoad] = determineTimeFrames(spec_directory)
dirlist = dir(spec_directory);
if isempty(dirlist)
    foldername = '_mKate_flat';
else
    foldername = char(dirlist.name);
end
cd (foldername)

files = dir('*.tif');
[~,~,~,chlist] = regexp([files.name],'_t[0-9]+');
numbrsCell = cellfun(@(x) str2double(x(3:end)),chlist,'UniformOutput',0);
numbrsMat = cellfun(@(x) str2double(x(3:end)),chlist,'UniformOutput',1);
numJump = min(numbrsMat)-1;
maxnumJump = max(numbrsMat);
timeFrames = cellfun(@(x) doTime(x,numJump,maxnumJump,1,chlist),numbrsCell,'UniformOutput',0);
frameToLoad = cellfun(@(x) doTime(x,numJump,maxnumJump,0,chlist),numbrsCell,'UniformOutput',0);
timeJump = numJump;
end
function time = doTime(x,numJump,maxnumJump,opt,chlist)
if opt ==1
digitz = length(chlist{1})-1;%determineNumber of digits in max time number;
time = num2str(zeros(digitz,1))';
time(1) = 't';
if length(num2str(x-numJump))>1
% time(2:end) = num2str(x-numJump);
time(length(time)-length(num2str(x-numJump))+1:end) = num2str(x-numJump);
else
   time(end) = num2str(x-numJump);
end

else
digitz = length(chlist{1})-1;%determineNumber of digits in max time number;
time = num2str(zeros(digitz,1))';
time(1) = 't';
if length(num2str(x))>1
% time(2:end) = num2str(x);
time(length(time)-length(num2str(x))+1:end) = num2str(x);
else
   time(end) = num2str(x);
end

end

end
function ImageDetails = InitializeImageDetails

ImageDetails.Scene=[];
ImageDetails.Channel=[];
ImageDetails.Frame=[];

end


%% load up images
function FinalImage=loadUpFinalImageOfStack(filenames)
cfile = filenames;
FileTif = char(cfile);
InfoImage=imfinfo(FileTif);

mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
bb=1;
i = 1;

i = NumberImages;
   TifLink.setDirectory(i);
   FinalImage=TifLink.read();
   
 
end


function FinalImage=loadUpTiffStack(filenames)
cfile = filenames;
FileTif = char(cfile);
InfoImage=imfinfo(FileTif);

mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
    for i = 1:NumberImages;
       TifLink.setDirectory(i);
       FinalImage(:,:,i)=TifLink.read();
    end
end

function FinalImage=loadUpTiffStackFrame(filenames,frame)
cfile = filenames;
FileTif = char(cfile);
InfoImage=imfinfo(FileTif);

mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,1,'uint16');
 
TifLink = Tiff(FileTif, 'r');
    i = frame;
       TifLink.setDirectory(i);
       FinalImage=TifLink.read();
    
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




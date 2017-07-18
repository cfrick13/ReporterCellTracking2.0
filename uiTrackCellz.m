function uiTrackCellz(FileDate,AutoTrackStr)
global Tracked pStruct timeVec timeSteps DivisionStruct xAxisLimits DICimgstack dfoName cfoName trackingPath background_seg bfoName nfoName sfoName cell_seg nucleus_seg segmentimgstack channelimgstack segmentPath mstackPath runIterateToggle ExportNameKey ExportName exportdir plottingTotalOrMedian channelinputs updateContrastToggle cmapper tcontrast lcontrast ThirdPlotAxes SecondPlotAxes expDateStr plotSettingsToggle PlotAxes cmap refineTrackingToggle expDirPath  timeFrames frameToLoad ImageDetails MainAxes SceneList displayTrackingToggle imgsize ExpDate

    DivisionStruct = struct();
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
    updateContrastToggle=0;
    plottingTotalOrMedian = 'median';
    tcontrast = 99;
    lcontrast = 1;
    runIterateToggle =0;
    refineTrackingToggle = 1;
    ImageDetails = InitializeImageDetails;
    displayTrackingToggle = 0;
    plotSettingsToggle =0;
    xAxisLimits = [0 50];
    
    
%determine export details
    ExportNameKey = 'tsi_';
    if ~strcmp(ExportNameKey,'final')
        disp(strcat('Export name key is "',ExportNameKey,'" not FINAL'))
    end
    ExportName = 'fricktrack';
   
    
%initialize global variables  
    channelimgstack =[];
    segmentimgstack =[];
    cfoName = [];
    sfoName =[];
    bfoName = [];
    nfoName = [];
    dfoName=[];
    DICimgstack=[];




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
    if isempty(FileDate)
    expDirPath = uigetdir;
    else
    cd(gparentdir)
    expdir = strcat(gparentdir,'/',FileDate);
    cd(expdir)
    expDirPath = pwd;
    end

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
    ExpDate = expDirPath(a:b+6);expDateStr = expDirPath(a:d); 
    [a,~] = regexp(ExpDate,'_');ExpDate(a) = '-';


%subdirectories should include
    %> [ flat mstack ]
    %> [ mstack images ]
    %> [ segment mstack ]
    %> [ tracking files ]


%load helpful metadata
    cd(exportdir)
    FileName = expDateStr;
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
    f = figure(1);
    initialsize = [1 1 2560 1080];
    f.Units = 'pixels';
    figdim = [initialsize(3) initialsize(4)];
    fpos = [1 1 figdim(1) figdim(2)];
    f.Position = fpos;
    fW = fpos(3);
    fH = fpos(4);


    bW = initialsize(3)./11; %button width
    bH = initialsize(4)/36; %button height
    yP = sort(100:bH:(fH-(fH/9)),'descend');
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
        'String','GoToFrame',...
        'Position',[xP(mmm)-(bW./2)-bW,yP(mmm),bW,bH],...
        'Callback',@gotobutton_callback);

%next row is final Frame and first Frame commands   
    mmm=3;
    uicontrol('Style','pushbutton',...
        'String','FinalFrame [g]',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@finalbutton_callback);
    uicontrol('Style','pushbutton',...
        'String','FirstFrame [z]',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
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
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@addareabutton_Callback);
    uicontrol('Style', 'pushbutton', 'String', 'Remove area',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@removeArea_Callback);
       
%delete commands
    mmm=8; 
    uicontrol('Style','pushbutton','String','Delete',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@deletebutton_Callback);
    uicontrol('Style','pushbutton','String','DeleteAllOnFrame',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@deleteAllonFrame_Callback);


%destroy commands    
    mmm=9;
    uicontrol('Style','pushbutton','String','Destroy [d]',...
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
    uicontrol('Style','pushbutton','String','Chosen Ones',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@chosenOnes_Callback);
    uicontrol('Style','pushbutton','String','Chosen OnesAllOnFrame',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@chosenOnesAllOnFrame_Callback);


%erode or dilate nuclei       
    mmm=12;
    uicontrol('Style','pushbutton','String','Erode Selected Nuclei',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@erodeOnes_Callback);
    uicontrol('Style','pushbutton','String','Dilate Selected Nuclei',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@dilateOnes_Callback);


%tracking commands
    mmm=14;
    uicontrol('Style','pushbutton','String','LinkCells [r]',...
        'Position',[xP(mmm),yP(mmm),bW,bH],...
        'Callback',@linkCells_Callback);
    mmm=15;
h=    uicontrol('Style','pushbutton',...
        'String','Run Tracking [t]',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@trackbutton_Callback);    
h=    uicontrol('Style','pushbutton',...
        'String','divisionTrack [d]',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'BackgroundColor',[1 0.7 0.7],...
        'Callback',@divisionTrack_Callback);        
    
    mmm=16;
    uicontrol('Style','pushbutton',...
        'String','LoadTracking',...
        'Position',[xP(mmm),yP(mmm),bW,bH],...
        'Callback',@loadTrackingFile_callback);


%display commands
    mmm=18;
    uicontrol('Style','pushbutton',...
        'String','contrast user',...
        'Position',[xP(mmm),yP(mmm),bW,bH],...
        'Callback',@contrast_Callback);
    uicontrol('Style','pushbutton',...
        'String','displayTrackingToggle [m]',...
        'Position',[xP(mmm),yP(mmm)-bH./2,bW,bH./1.5],...
        'Callback',@displayTrackingButton_Callback);
        

%save commands
    mmm=20;
    uicontrol('Style','pushbutton',...
        'String','SaveTrackingAs',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW+bW,bH],...
        'Callback',@saveTrackingFileAs_callback);
    
%autotrack commands
    mmm=21;
    uicontrol('Style','pushbutton',...
            'String','trackSaveIterate',...
            'Position',[xP(mmm)+bW./2,yP(mmm),bW+bW./2,bH./1.5],...
            'Callback',@trackSaveIterate_callback);
    uicontrol('Style','pushbutton',...
        'String','TSIchosen',...
        'Position',[xP(mmm)-bW,yP(mmm),bW+bW./2,bH./1.5],...
        'Callback',@trackSaveIterateChosen_callback);


%plot and plot settings commands
    mmm=23;  
    uicontrol('Style','pushbutton',...
        'String','PLOT!',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW,bH],...
        'Callback',@Plot_callback);
    uicontrol('Style','pushbutton',...
        'String','Plot Specific Cell!',...
        'Position',[xP(mmm)+bW./2,yP(mmm),bW,bH],...
        'Callback',@Plot_SpecificCell_callback);
    mmm=24;
    uicontrol('Style','pushbutton',...
        'String','Plot Settings!',...
        'Position',[xP(mmm)-bW./2,yP(mmm),bW+bW./2,bH./1.5],...
        'Callback',@PlotSettings_callback);

%export commands  
    mmm=26;
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
    mmm=28;
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

        
%refine figure details
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
    imgdim = (fW./2).*0.95;
    Position = [10 10 imgdim imgdim];
    % Position = [0.1 0.3 0.65 0.65];
    MainAxes.Position = Position;
    MainAxes.Units = 'normalized';


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

    f.Units = 'pixels';
    scrnsize = get(0,'screensize');
    fpos = scrnsize;
    fpos(1) = scrnsize(3)/2;
    fpos(3) = scrnsize(3)/2;
    fpos(4) = scrnsize(4).*0.75;
    f.Position = fpos;
    
%constrain image axis to be square initially
    MainAxes.Units = 'pixels';
    pos = MainAxes.Position;
    pos(3:4) = [min([pos(3) pos(4)]) min([pos(3) pos(4)])];
    MainAxes.Position = pos;
    MainAxes.Units = 'normalized'; 
    pos = MainAxes.Position;
    pos(2)= 0.5 - pos(4)/2;
    MainAxes.Position = pos;
    
    
    f.Color = 'w';
    set(f,'KeyPressFcn',@keypress);
    
    
    pStruct = loadSegmentParameters([],FileName,exportdir); %loads saved value of pStruct
    timeMat = loadMetadata(FileName,exportdir);
    rtimeMat = round(timeMat);
    stimulationFrame = dosestruct(1).tgfFrame + 1;
    timeVec = rtimeMat(1,:)-rtimeMat(1,stimulationFrame);
    timeSteps = diff(round(timeVec),[],2);
    
    if strcmpi(AutoTrackStr,'AutoTrackCells')
        f.Visible = 'off';
        trackSaveIterateChosen_callback([],[])
        exportTrackedCells([],[])
        % exportNuclei([],[])
        close(f)
    elseif strcmpi(AutoTrackStr,'AutoExportNuclei')
        exportSegmentedCells([],[])
        close(f)
    else
        Tracked = makeTrackingFile(timeFrames);
        ImageDetails.Scene = SceneList{1};
        ImageDetails.Channel = nucleus_seg;
        ImageDetails.Frame = 1;
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


%% uifunctions

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
%         destroybutton_Callback([],[]);
        divisionTrack_Callback([],[]);
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
            [~] = getxy([],[]);
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
    [~,comments,commentpos,~]=updatecomments(xy);
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


function divisionTrack_Callback(~,~)
global imgsize ImageDetails Tracked MainAxes DivisionStruct

t = ImageDetails.Frame;
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.

    MainAxes.Title.String = 'click mother cell...';
     [cellxx,cellyy] = ginput(1);
       
    cellx = round(cellxx);
    celly = round(cellyy);
    
    mothercellind = sub2ind(imgsize,celly,cellx);
    
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
    
    daughtercellinds = sub2ind(imgsize,celly,cellx);
    
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
    global  ImageDetails Tracked DivisionStruct 

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
global frameToLoad ImageDetails 

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
idx = str2double(cell2mat(inputdlg(prompt,dlg_title)));

frameToLoad = idx;
ImageDetails.Frame = frameToLoad;
setSceneAndTime
end


%choose scenes
function nextscenebutton_Callback(~,~) 
global   ImageDetails SceneList  




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
global   ImageDetails SceneList  

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
global  ImageDetails  Tracked imgsize  refineTrackingToggle
 % choose cell
%       [cellx,celly] = ginput(1);
       % construct a polygon to add

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


     refineTrackingToggle = 1;
        nextbutton_callback([],[]);
      end 
      end

end

%delete cells
function deletebutton_Callback(~,~) 
    global imgsize ImageDetails  Tracked refineTrackingToggle
    

t = ImageDetails.Frame;
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
    
PX(~idxs) = {NaN};
CC.PixelIdxList = PX;
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Tracked{t}.Cellz = CC;
    
refineTrackingToggle = 1; 
setSceneAndTime;
   
end
function eliminatebutton_Callback(~,~)

global ImageDetails Tracked imgsize refineTrackingToggle

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
   refineTrackingToggle = 1; 
    setSceneAndTime;
end

end
end

function destroyAllFramePrevious_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked refineTrackingToggle

%   determine the frame to load
    t = ImageDetails.Frame;
    CC = Tracked{t}.Cellz;
    PX = CC.PixelIdxList;   

    idxs = false(size(PX));   %choose all cells on frame
    Trackedz = crushThem(Tracked,idxs,[],t); %delete from t=t backwar
    Tracked = Trackedz;
    refineTrackingToggle = 1; 
    setSceneAndTime;

end
function destroyAllFrameSubsequent_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked  refineTrackingToggle

%   determine the frame to load
    t = ImageDetails.Frame;
    CC = Tracked{t}.Cellz;
    PX = CC.PixelIdxList;   

    idxs = false(size(PX));   %choose all cells on frame
    Trackedz = crushThem(Tracked,idxs,t,[]); %delete from t=t onward
    Tracked = Trackedz;
    refineTrackingToggle = 1; 
    setSceneAndTime;
end
function destroybuttonAllPrevious_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked imgsize refineTrackingToggle

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


    refineTrackingToggle = 1; 
    setSceneAndTime;



end
function destroybuttonAllSubsequent_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked imgsize refineTrackingToggle

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


    refineTrackingToggle = 1; 
    setSceneAndTime;



end

function deleteAllonFrame_Callback(~,~)
%delete all cells from one frame
  global  ImageDetails  Tracked refineTrackingToggle
    
    
    
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
    
    refineTrackingToggle = 1; 
    setSceneAndTime;

end



function destroybutton_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked imgsize refineTrackingToggle

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


    refineTrackingToggle = 1; 
    setSceneAndTime;


end
%choose the cells you want
function chosenOnes_Callback(~,~)
%choose the cells you want
global ImageDetails Tracked imgsize refineTrackingToggle

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


    refineTrackingToggle = 1; 
    setSceneAndTime;


end
function chosenOnesAllOnFrame_Callback(~,~)
%choose the cells you want
global ImageDetails Tracked refineTrackingToggle

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


    refineTrackingToggle = 1; 
    setSceneAndTime;


end

function erodeOnes_Callback(~,~)
%choose the cells you want
global ImageDetails Tracked imgsize refineTrackingToggle

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

    refineTrackingToggle = 1; 
    setSceneAndTime;


end



function dilateOnes_Callback(~,~)
%choose the cells you want
global ImageDetails Tracked imgsize refineTrackingToggle

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


    refineTrackingToggle = 1; 
    setSceneAndTime;

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
 global  ImageDetails Tracked imgsize refineTrackingToggle
 % choose cell
%       [cellx,celly] = ginput(1);
       % construct a polygon to add

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
      end
      
      
      refineTrackingToggle = 1;
    nextbutton_callback([],[]);
      end 
      end
end
function linkCells_Callback(~,~)
global ImageDetails Tracked refineTrackingToggle imgsize timeFrames

%this is a quick workaround to get linking working with tracking
%trajectories on
% refineTrackingToggle =1;
% setSceneAndTime

button=1;
i=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while button==1
refineTrackingToggle=0;
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
refineTrackingToggle=1;

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
global exportdir expDateStr frameToLoad
cd(exportdir)

queryName = strcat(expDateStr,'*DoseAndScene*.mat');
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
global cell_seg nucleus_seg xAxisLimits trunccmaplz SecondPlotAxes Tracked ImageDetails expDirPath timeFrames mstackPath frameToLoad PlotAxes imgsize plotSettingsToggle psettings cmaplz displayTrackingToggle cmap

    if plotSettingsToggle == 0
        psettings = PlotSettings_callback([],[]);
        plotSettingsToggle=1;
    end
    framesThatMustBeTracked = psettings.framesThatMustBeTracked;


    for jy = 1
        PX = Tracked{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
        makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
    end
    [~,~,~,plotidx]=commentsforplot(Tracked);
    makeIMGidx = find(makeIMG==1);
    
    if ~isempty(plotidx)
        iidd = find(~ismember(makeIMGidx,plotidx));
    else
        iidd=[];
    end


    smooththat=0;
    [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,expDirPath,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plotSettingsToggle,psettings,makeIMG,makeIMGidx,smooththat);

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
%         plotMat = plotStructUI.mkate;
%         plotMatFC = plotStructUI.mkateFC;
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
            
    if displayTrackingToggle ==1
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
    
    if displayTrackingToggle ==1
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
    PlotAxes.Title.String = [ImageDetails.Channel ' in tracked cells'];

end


function Plot_SpecificCell_callback(~,~)
global displayTrackingToggle trunccmaplz xAxisLimits ThirdPlotAxes toggleCFPnorm Tracked ImageDetails expDirPath mstackPath timeFrames frameToLoad PlotAxes imgsize plotSettingsToggle psettings

if plotSettingsToggle == 0
psettings = PlotSettings_callback([],[]);
plotSettingsToggle=1;
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
[plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,expDirPath,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plotSettingsToggle,psettings,makeIMG,makeIMGidx,smooththat);
    
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
            
    if displayTrackingToggle ==1
        for i=1:length(h)
            h(i).Color = cmapl(makeIMGidx(i),:);
        end 
    end

end

function [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,expDirPath,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plotSettingsToggle,psettings,makeIMG,makeIMGidx,smooththat)
global cell_seg nucleus_seg background_seg segmentPath

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
Nucbkg = zeros(1,timeFrames,'single');
for k=1:timeFrames
    bkglog = ~bkglogimgstack(:,:,k);
    cellQ_img = single(cellQ_imgstack(:,:,k));
    nuc_img = single(nuc_imgstack(:,:,k));
    %background subtraction is just subtraction with a value
    Smadbkg(k) = nanmedian(cellQ_img(bkglog));
    Nucbkg(k) = nanmedian(nuc_img(bkglog));
    cellQ_imgstack(:,:,k) = cellQ_img-Smadbkg(k);
    nuc_imgstack(:,:,k) = nuc_img-Nucbkg(k);
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

    
    Smadtotal = zeros(size(cellQ_pxls),'single');
    mkatetotal = Smadtotal;
    Smad = Smadtotal;
    mkate = Smadtotal;
    for i = 1:size(cellQ_pxls,1)
        for j = 1:size(cellQ_pxls,2)
            Smadtotal(i,j) = sum(cellQ_pxls{i,j});
            mkatetotal(i,j) = sum(nuc_pxls{i,j});
            Smad(i,j) = median(cellQ_pxls{i,j});
            mkate(i,j) = median(nuc_pxls{i,j});
        end
    end
    Smadtotal(Smadtotal==single(13579)) = NaN;
    mkate(mkate==single(13579)) = NaN;
    mkatetotal(mkatetotal==single(13579)) = NaN;
    Smad(Smad==single(13579)) = NaN;


basalSUB = framesThatMustBeTracked(1)-7;
if basalSUB<1
    basalSUB = framesThatMustBeTracked-1;
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
    
    if smooththat==1
    %%%%%%%%%%%%%%%%%%%
    Smad = Smad./mkateFC;
    %%%%%%%%%%%%%%%%%%%
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
plotStructUI.SmadFC = SmadFC;
plotStructUI.mkateFC = mkateFC;
plotStructUI.mkateFCtotal = mkateFCtotal;
plotStructUI.Smadbkg = Smadbkg;
plotStructUI.mkatebkg = Nucbkg;
    
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
global Tracked plotSettingsToggle psettings imgsize ImageDetails     
    if plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    plotSettingsToggle=1;
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
global ImageDetails  

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


end
function setcommentsTracking(comments,commentpos)
global Tracked
t=length(Tracked);
Tracked{t}.comments= comments;
Tracked{t}.commentpos = commentpos;
end
function [cellnums,comments,commentpos,cellidx]=updatecomments(xy)
%choose the cells you want to comment on
global ImageDetails  Tracked imgsize 

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
global ExportNameKey ExportName  displaycomments SceneList Tracked   trackPath     plotSettingsToggle psettings

    if plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    plotSettingsToggle=1;
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
global cell_seg nucleus_seg background_seg segmentPath ExportNameKey ExportName exportdir mstackPath expDateStr SceneList  trackingPath timeFrames plotSettingsToggle psettings
exportStruct = struct();

    if plotSettingsToggle == 0
        psettings = PlotSettings_callback([],[]);
        plotSettingsToggle=1;
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
            filename = strcat(expDateStr,'_tracking_export.mat'); 
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
global cell_seg nucleus_seg background_seg segmentPath   exportdir mstackPath expDateStr SceneList  trackingPath timeFrames


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
            filename = strcat(expDateStr,'_nuclei_export.mat'); 
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
global ExportNameKey ExportName exportdir  cell_seg nucleus_seg background_seg segmentPath expDateStr SceneList   mstackPath timeFrames trackingPath    plotSettingsToggle psettings
exportStruct = struct();

    if plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    plotSettingsToggle=1;
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
%                 [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,expDirPath,ImageDetails,SceneDirPath,timeFrames,frameToLoad,PlotAxes,imgsize,plotSettingsToggle,psettings,makeIMG,makeIMGidx,smooththat);
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
        filename = strcat(expDateStr,'_tracking_allcells_export.mat'); 
        save(filename,'exportStruct');
end

function xy = labelCells(~,~)


global ExportNameKey ExportName displaycomments   Tracked ImageDetails  trackPath  frameToLoad  imgsize plotSettingsToggle psettings

    if plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    plotSettingsToggle=1;
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
global ExportNameKey ExportName  displayTrackingToggle   ImageDetails  trackPath  frameToLoad   plotSettingsToggle psettings

    if plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    plotSettingsToggle=1;
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
            displayTrackingToggle =1;
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
global ExportNameKey ExportName  imgsize displayTrackingToggle SceneList  ImageDetails  trackPath  plotSettingsToggle psettings updateContrastToggle

    if plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    plotSettingsToggle=1;
    end
    
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(trackPath)
cd ..
CENTROID = struct();
    for scenenumber = 1:length(SceneList)
        updateContrastToggle=1;
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
            displayTrackingToggle =1;
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
    updateContrastToggle=0;
end


%% tracking functions
function displayTrackingButton_Callback(~,~)
global    displayTrackingToggle Tracked


if isempty(Tracked{1}.Cellz)
    disp('need to run tracking!!!')
else
    if displayTrackingToggle == 0
    displayTrackingToggle = 1;
    else
    displayTrackingToggle =0;
    end
    setSceneAndTime
end


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
global  Tracked ImageDetails refineTrackingToggle segmentPath nucleus_seg mstackPath cell_seg background_seg
pvalue = ImageDetails.Scene;
    tic
    trackfilelist = {'yes','no'};
    [S,~] = listdlg('PromptString','Are you sure you want to run tracking?',...
                'SelectionMode','single',...
                'ListSize',[200 300],...
                'ListString',trackfilelist);
            
            if S==1
                h=waitbar(0,'running tracking algorithm...');
            Tracked = FrickTrackCellsYeah(segmentPath,mstackPath,pvalue,nucleus_seg,cell_seg,background_seg,h);
            else
            end

            refineTrackingToggle =1;
            close(h);
setSceneAndTime;
    toc
end

function Tracked = loadTrackedStructure
global trackingPath timeFrames refineTrackingToggle runIterateToggle ImageDetails
if runIterateToggle ==0
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
        refineTrackingToggle=0;
        else
        refineTrackingToggle =1;
        end

    else
        Tracked = makeTrackingFile(timeFrames);
    end
else
    Tracked = makeTrackingFile(timeFrames);
end


end
function loadTrackingFile_callback(~,~)
global  Tracked refineTrackingToggle

Tracked = loadTrackedStructure;
refineTrackingToggle =1;
end

%make trajectories for overlay of tracking
function traject = trackingTrajectories(timeFrames)
global Tracked


%   determine the frame to load
    t = timeFrames;

    xy = cell(1,t);
    lxy = zeros(1,t);

    for i = 1:t %determine centroids from 1:t for plotting a tracking tail 
        CC = Tracked{i}.Cellz;
        Centroids  = CC.Centroid;
        xy{i} = Centroids';
        lxy(i) = length(xy{i});
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
function Trackedz = trackingCosmetics(Stacked)
global displayTrackingToggle
%this identifies the maximum length to identify the total number of cells
%MARlength contains the number of cells in each frame
    MARlength = zeros(1,length(Stacked));
    for i=1:length(Stacked)
    MAR = Stacked{i}.Cellz.PixelIdxList;
    MARlength(i) = length(MAR);
    end

    %this makes a cell for each cell in all of the frames so each frame has
    %the same total number of cells. (When track is ended NaN is added)
    mmarl = max(MARlength);
    PXBOX = cell(mmarl,length(Stacked));
    for i=1:length(Stacked)
    MAR = Stacked{i}.Cellz.PixelIdxList;
    PX = cell(1,mmarl);
    PX(1:MARlength(i)) =  MAR;
    PX(MARlength(i)+1:max(MARlength)) =  {NaN};
    Stacked{i}.Cellz.PixelIdxList = PX;
    PXBOX(:,i) = PX;
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
                disp('correcting tracking...')
                beginoftrack = tracks(trackidx,1);
                endoftrack = tracks(trackidx,2);
                for frame = 1:length(Stacked)
                    S = Stacked{frame};
                    C = S.Cellz;
                    PX = C.PixelIdxList;
%                     PX = Stacked{frame}.Cellz.PixelIdxList;
                    px = cell(1,length(PX)+1);
                    px(1:length(PX)) = PX;
                    if frame>beginoftrack-1 && frame<endoftrack+1 % if frame is where cell is tracked
                        px(length(PX)+1) = PX(j); %add the track to the end
                        px(j) = {NaN}; %remove the track data from where it was previously located
                    else %if the frame is not where the cell is located fill the track with NaN
                        px(length(PX)+1) = {NaN};
                    end
                    C.PixelIdxList = px;
                    S.Cellz = C;
                    Stacked{frame} = S;
%                     Stacked{frame}.Cellz.PixelIdxList = px;
                end
                
                
                %try tracking nucleus number (as opposed to all pixel data)
                %and then use the changes below in TRT to track nuclear
                %number instead of just ones and zeros and then you can
                %easily update the PX with the new TRT. 
                
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
    connectivity = Stacked{1}.Cellz.Connectivity;
    imgsize = Stacked{1}.Cellz.ImageSize;
        for j=1:length(Stacked)
        S = Stacked{j};
        CC = S.Cellz;
        PX = CC.PixelIdxList;
%         PX = Stacked{j}.Cellz.PixelIdxList;
        px = PX(pxidx);
%         CC = Stacked{j}.Cellz;
        CC.Connectivity = connectivity;
        CC.ImageSize = imgsize;
        CC.NumObjects = length(px);
        CC.PixelIdxList = px;
        
        if displayTrackingToggle==1
            stats = regionprops(CC,'Centroid');
            centroidarray = {stats.Centroid};
            Centroid = zeros(2,length(centroidarray));
            for ijk = 1:length(centroidarray)
                Centroid(:,ijk) = centroidarray{ijk};
            end
            CC.Centroid = Centroid;
        end
        
%         CC = checkDivision()
        S.Cellz = CC;
        Stacked{j} = S;
        end
    Trackedz=Stacked;


end
function trackSaveIterate_callback(~,~)
global parellelToggle runIterateToggle SceneList ImageDetails refineTrackingToggle   Tracked trackingPath ExportName timeFrames segmentPath nucleus_seg mstackPath cell_seg background_seg

runIterateToggle =1;
parellelToggle = 1;
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
            
            Tracked = FrickTrackCellsYeah(segmentPath,mstackPath,pvalue,nucleus_seg,cell_seg,background_seg,[]);
%             Tracked = FrickTrackCellsYeah(expDirPath,frameToLoad,pvalue,[]);
            refineTrackingToggle =1;
            setSceneAndTime;


        %save
            % saveTrackingFileAs_callback([],[])
            cd(trackingPath)
            filename = 'initial';
            save(strcat(filename,'_',ImageDetails.Scene,'_',ExportName,'.mat'),'Tracked')
%             save(strcat(filename,ExportName,'.mat'),'Tracked')
            
    end
runIterateToggle =0;
parellelToggle = 0;

end
function trackSaveIterateChosen_callback(~,~)
global parellelToggle plotSettingsToggle psettings runIterateToggle SceneList ImageDetails refineTrackingToggle   Tracked trackingPath ExportName timeFrames segmentPath nucleus_seg mstackPath cell_seg background_seg

runIterateToggle =1;
parellelToggle = 1;
    for i=1:length(SceneList)
%     for i=22
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
            
            Tracked = FrickTrackCellsYeah(segmentPath,mstackPath,pvalue,nucleus_seg,cell_seg,background_seg,[]);
%             Tracked = FrickTrackCellsYeah(expDirPath,frameToLoad,pvalue,[]);
            refineTrackingToggle =1;
            setSceneAndTime;
            
            cd(trackingPath)
            filename = 'tsi';
            save(strcat(filename,'_',ImageDetails.Scene,'_',ExportName,'.mat'),'Tracked')
            
            
            
    if plotSettingsToggle == 0
        psettings = PlotSettings_callback([],[]);
        plotSettingsToggle=1;
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
runIterateToggle =0;
parellelToggle = 0;

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
function valProbette = detProbette(val,valPrev,fidx,ijk)
   valprobsub = abs(val(fidx) - valPrev(ijk));
   valprosubvec =abs(valPrev - val(fidx));
   valprobette = valprobsub/max(valprosubvec(:));
   valProb = 1-valprobette;
   valProbette = valProb';
end
function valnew = updatevaluefunc(val,oldidx,newidx)
% validx = true(1,size(val,1));
% valnew = val(validx,idx);
% valnew(validx,loserz) = NaN;
% valnew = horzcat(valnew,val(validx,missers));
validx = true(1,size(val,2));
if iscell(val)
valnew = cell(length(oldidx),size(val,2));
else
    valnew = nan(length(oldidx),size(val,2));
end

idx1 = ~isnan(oldidx)&~(isnan(newidx));
valnew(idx1,validx) = val(newidx(idx1),validx);

idx2 = isnan(oldidx)&~isnan(newidx);
valnew(idx2,validx) = val(newidx(idx2),validx);

if iscell(val)
idx3 = cellfun(@isempty,valnew,'UniformOutput',1);
valnew(idx3,validx) = {NaN};
end

end


function [distProb,idx,cutoffidx] = probSpitter(input,inputPrev,knnnum,displacementCutoff)
%make a matrix that tells you the probability a value in one vector is the same as anothe

    
    [idx,eps] = knnsearch(input,inputPrev,'K',knnnum); 
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
    
%     distProbValues = nan(size(eps));
%     for i = 1:size(distinputvec,1)
% %        valprobsub = abs(val(fidx) - valPrev(ijk));
% %        valprobette = valprobsub/max(valprobsub);
% %        valProbette = 1-valprobette;
%         distvec_norm = (distvec(i,:) - min(distvec(i,:)))./(max(distvec(i,:)) - min(distvec(i,:)));
%         distProbValues(i,:) = 1 - distvec_norm;
%     end


    
    distProb = zeros(size(inputPrev,1),size(input,1));
    epsMat = distProb;
    for di = 1:size(idx,1)
       didx = idx(di,:);
       distProb(di,didx) = distProbValues(di,:);
       epsMat(di,didx) = eps(di,:);
    end
    
    cutoffidx = epsMat>displacementCutoff;
%     distProb(cutoffidx)=0;
    
    
end
%function for tracking cells
function [ Tracked ] = FrickTrackCellsYeah(segmentPath,mstackPath,pvalue,nucleus_seg,cell_seg,background_seg,h)
global pStruct timeVec timeSteps
%function for tracking cells
%   Detailed explanation goes here


    %determine the size of the nuclei in each experiment
    nucleiDist = pStruct.nucleus.nucDiameter;
    


    Frame = struct();
    %load nucleus segmentation binary as segmentimgstack
    cd(segmentPath)
    ff = dir(strcat('*',pvalue,'*',nucleus_seg,'*'));
    if length(ff)>1
        ff = dir(strcat('*',pvalue,'*',nucleus_seg,'*nucleus*')); 
    end
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
    filename = char(ff.name);
    nucleusFileObject = matfile(filename);
    segmentimgstack = nucleusFileObject.IfFinal;


    %load background segmentation binary as segmentimgstack
    cd(segmentPath)
    ff = dir(strcat('*',pvalue,'*',background_seg,'*'));
    if length(ff)>1
        ff = dir(strcat('*',pvalue,'*',background_seg,'*background*')); 
    end
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
    filename = char(ff.name);
    backgroundFileObject = matfile(filename);
    backgroundimgstack = backgroundFileObject.IfFinal;

    %load nucleus images as nosub_nucleusimgstack
    cd(mstackPath)
    ff = dir(strcat('*',pvalue,'*',nucleus_seg,'*'));
    if length(ff)>1
        ff = dir(strcat('*',pvalue,'*',nucleus_seg,'*nucleus*')); 
    end
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
    filename = char(ff.name);
    nucleusFileObject = matfile(filename);
    nosub_nucleusimgstack = nucleusFileObject.flatstack;   

    %load cell images as nosub_cellimgstack
    cd(mstackPath)
    ff = dir(strcat('*',pvalue,'*',cell_seg,'*'));
    if length(ff)>1
        ff = dir(strcat('*',pvalue,'*',cell_seg,'*cell*')); 
    end
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
    filename = char(ff.name);
    cellFileObject = matfile(filename);
    nosub_cellimgstack = cellFileObject.flatstack;      
    

    %background subtract nucleusimgstack
    nucleusimgstack = zeros(size(nosub_nucleusimgstack));
    cellimgstack = nucleusimgstack;
    for i=1:size(nosub_nucleusimgstack,3)
        nucI = nosub_nucleusimgstack(:,:,i);
        cellI = nosub_cellimgstack(:,:,i);
        bkgI = backgroundimgstack(:,:,i);
        nucleusimgstack(:,:,i) = nucI-nanmedian(nucI(bkgI));
        cellimgstack(:,:,i) = cellI - nanmedian(cellI(bkgI));
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
%     segmentsequence = fliplr(1:size(segmentimgstack,3)); %track from last frame to first frame
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
            segment_nucFluor_array{i} = cellfun(@(x) nanmedian(nucI(x)),PX,'UniformOutput',1)';
            segment_cellFluor_array{i} = cellfun(@(x) nanmedian(cellI(x)),PX,'UniformOutput',1)';
            segment_Stdev_array{i} = cellfun(@(x) nanstd(nucI(x)),PX,'UniformOutput',1)';
            segment_Centroid_array{i} = centroidMat;
            segment_Pixels_array{i} = PX;
            segment_Cellnum_array{i} = [1:length(PX)]';
        end
        
    %track based on assigning probabilities determined by minimizing changes to measured parameters
    
        
%     segmentsequence = fliplr(1:size(segmentimgstack,3)-1); %track from last frame to first frame
    segmentsequence = 2:size(segmentimgstack,3); %track from first frame to last frame
    birtharray = cell(1,size(segmentimgstack,3));
    nidxarray = cell(1,size(segmentimgstack,3));
    birtharray(1) = segment_Cellnum_array(1);

               
    possibleWorkers = feature('numcores');
    nWorkers = possibleWorkers;
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

        parfor i = segmentsequence
            a=i-1;
            b=i;
            centroids = segment_Centroid_array{b};
            centroidsPrev = segment_Centroid_array{a};
               if isempty(centroids) %that is, no cells were segmented
                   b=i-1;
                   centroids = segment_Centroid_array{b};
               end
               if isempty(centroidsPrev) %that is, no cells were tracked on previous frame
                   a=i;
                   centroidsPrev = segment_Centroid_array{a};
               end
            area = segment_Area_array{b};
            areaPrev = segment_Area_array{a};
            nucfluor = segment_nucFluor_array{b};
            nucfluorPrev = segment_nucFluor_array{a};
            cellfluor = segment_cellFluor_array{b};
            cellfluorPrev = segment_cellFluor_array{a};            
            pixels = segment_Pixels_array{b};
            pixelsPrev = segment_Pixels_array{a};
            ellipt = segment_Ellipt_array{b};
            elliptPrev = segment_Ellipt_array{a};
            
            displacementCutoff = (nucleiDist)*(timeSteps(i-1)./10);

            
               %nearest neighbor distances to determine probability that
               %two cells are the same
                knnnum = 5;
               [distProb,~,distcut] = probSpitter(centroids,centroidsPrev,knnnum,displacementCutoff);
               [areaProbette,~,areacut] = probSpitter(area,areaPrev,size(area,1),nanmedian(area));
               [nucFluorProbette,~,nuccut] = probSpitter(nucfluor,nucfluorPrev,size(nucfluor,1),nanmedian(nucfluor)./2);
               [cellFluorProbette,~,cellcut] = probSpitter(cellfluor,cellfluorPrev,size(cellfluor,1),Inf);
               newProb = distProb.*areaProbette.*nucFluorProbette.*cellFluorProbette;     
               %distProb dimesionas are length(centroidsPrev) x length(centroids)
               
              %you should assign weights to the probilities as well

               %centroids [37x2] centroidPrev[34x2] idx[34x3] distProb[34x37];
               %idx(1,:) = [1 5 9] means that centroids(1,:)  centroids(5,:)and centroids(9.:) are closest to centroidPrev(1,:)
               %now you need to build a probability matrix that has rows
               % for all cells in centroidsPrev (size(CentroidsPrev,1) and columns for all cells in centroids (size(centroids,1))
                trackProb = distProb;
                trackProb(distcut|nuccut)=0;
               [maxvals,idx] = max(trackProb,[],2); %idx is index of input that matches to inputPrev such that input(idx) = inputPrev;
               
                %now some cells are assigned twice. Correct this based on highest probabilities
                %num_cells_set should be = the size of the current NOT the prev
                num_cells_currentFrame = 1:size(trackProb,2); %size(cellProb,2) = length(input)
                num_cells_prevFrame = 1:size(trackProb,1); %size(cellProb,2) = length(input)
                
                arbitrateProb = newProb;
                arbitrateProb(distcut|nuccut|areacut)=0;
                [arbmaxvals,~] = max(arbitrateProb,[],2); %idx is index of input that matches to inputPrev such that input(idx) = inputPrev;
                arbmaxvals(arbmaxvals==0)=NaN;
                idx(maxvals==0) = NaN;
                idx(isnan(maxvals)) = NaN;
                [n, bin] = histc(idx, num_cells_currentFrame);
                multiple = find(n>1); %the same cell is called closest to two previous cells
                missers = find(n<1); %these are likely new cells
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
                            loserz = [loserz losern'];
                        end
                    end 
                    %the index of loserz represents the cell number of
                    %previous frame and the new indexed numbers of current
                    %frame
                    %i.e. loserz =3 means cell 3 from previous frame
                    %matches with index 3 of newidx
                    

                nidx = idx; 
                nidx(loserz)=NaN; %remove duplicates (loserz)// that is to say, tracks that merge onto one cell
                [nn, ~] = histc(vertcat(nidx,missers), num_cells_currentFrame);
                updatemissers = find(nn<1);
                
                
                newidx = vertcat(nidx,missers,updatemissers);
                oldidx = [1:length(newidx)]'; oldidx(oldidx>length(pixelsPrev))=NaN;
                oldidx((isnan(areaPrev)==1))=NaN;
            
            %display a table showing which cells become which
                dispidx = zeros(length(newidx),length([num2str(length(newidx)) ' -> ' num2str(max(newidx))]));
                dispidx = char(dispidx);
                for iter=1:size(dispidx,1)
                    origcellstr = num2str(oldidx(iter));
                    vec  = [origcellstr ' -> ' num2str(newidx(iter))];
                    dispidx(iter,1:length(vec)) = vec;
                end
%                 disp({'',['i = ' num2str(i)]})
%                 disp(dispidx)
            
            [nn, ~] = histc(newidx, 1:nanmax(newidx));
            if max(nn)>1
                error('cells called twice!!!')
            end

            
            birtharray(i) = {vertcat(missers,updatemissers)};
            nidxarray(i) = {nidx};
        end
  

%determine the number of births (new tracks) to determine the dimensions of
%the track vector
numbirths = zeros(1,length(birtharray));
for i = 1:length(nidxarray)
    births = birtharray{i};
    numbirths(i) = length(births);
end

number_of_tracks = sum(numbirths);
trackmatrix = nan(size(segmentimgstack,3),number_of_tracks);
birthdistold = 0;
for i = 1:size(segmentimgstack,3)
    births = birtharray{i};
    nidx = nidxarray{i}';
    nidvec = nidx(~isnan(nidx));

    if i>1
       prevvecinit = trackmatrix(i-1,1:birthdistold);
       prevvec = prevvecinit(~isnan(prevvecinit));
       %they are always the same size
        %        disp(size(prevvec))
        %        disp(size(nidx))
        %        disp(max(prevvec))
    else
       prevvec = 1:length(nidvec);
    end

    if ~isempty(nidx)
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
        trackmatrix(i,birthdistold+1:birthdistold+length(births)) = births;
        birthdistnew = max(birthdistold+1:birthdistold+length(births));
        birthdistold=birthdistnew;
    end
end


%now update all the fields
newpx = cell(size(trackmatrix));
newpx(:) = {NaN};
segI = segmentimgstack(:,:,1); %choose the first frame of segmentimgstack to start with
CC = bwconncomp(segI);
Tracked = cell(1,size(segmentimgstack,3));
    for i = 1:size(segmentimgstack,3)
        px = segment_Pixels_array{i};
        trackvals=trackmatrix(i,:);
        trackidx = ~isnan(trackvals);
        tracknums = trackvals(trackidx);
        newpx(i,trackidx) = px(tracknums);
        
        centroids = segment_Centroid_array{i};
        centnew = nan(length(trackidx),size(centroids,2));
        centnew(trackidx,:) = centroids(tracknums,:);
        

        AllCellsPX = newpx(i,:);
        CC.PixelIdxList = AllCellsPX;
        CC.NumObjects = numel(AllCellsPX);
%             stats = regionprops(CC,'Centroid');
%             Smat = vertcat(stats.Centroid);
        CC.Centroid = centnew;
        CC.Division = [];
        Frame.filename = filename;
        Frame.Cellz = CC;
        Tracked{i} = Frame;
    end
end
function saveTrackingFileAs_callback(~,~)
global  trackingPath Tracked ExportName ImageDetails
cd(trackingPath)
prompt = 'filename of tracking structure to be saved?';
dlg_title = 'save tracking structure as...specific filename';
filename = char(inputdlg(prompt,dlg_title));
save(strcat(filename,'_',ImageDetails.Scene,'_',ExportName,'.mat'),'Tracked')
end

%% Image Display functions
function setSceneAndTime
global runIterateToggle displayTrackingToggle refineTrackingToggle DICimgstack dfoName  nucleus_seg backgroundimgstack bfoName nfoName background_seg cell_seg nucleusimgstack sfoName segmentimgstack  channelimgstack cfoName segmentPath frameToLoad ImageDetails  Tracked SceneList  trackPath imgfile mstackPath

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
        refineTrackingToggle = 0;
        cd(segmentPath)
        imgfile = dir(strcat('*',ImageDetails.Scene,'*',nucleus_seg,'*'));
        If = segmentimg;
        Tracked = loadTrackedStructure;
        centroidDisagreement=false;
    else  %if there exists segmenttracking already...then load that. 
        CC = Tracked{t}.Cellz;
        PX = CC.PixelIdxList;
        Centroids  = CC.Centroid;
        centroidDisagreement = ~(size(Centroids,2)==CC.NumObjects);
    %     makeIMG = cellfun(@(x) length(x)==1,PX,'UniformOutput',1); %choose only the cells without NAN
        makeIMG = cellfun(@(x) length(x)<2,PX,'UniformOutput',1); %choose only the cells without NAN
        CC.PixelIdxList = PX(~makeIMG);
        CC.NumObjects = length(PX(~makeIMG));

        segmentimgL = labelmatrix(CC);
        segmentimgz = false(size(segmentimgL));
        segmentimgz(segmentimgL>0)=1;
        If = segmentimgz;
    end
    


if displayTrackingToggle==1
    if centroidDisagreement      
        refineTrackingToggle = 1;
    end
end
% disp(refineTrackingToggle)
if refineTrackingToggle == 1
Tracked = trackingCosmetics(Tracked);
end
refineTrackingToggle=0;


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
     channelimg(~backgroundimg) = channelimg(~backgroundimg)+1000;
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
if runIterateToggle ==0
displayImageFunct(If,channelimg,bkgmedian);
end

end


function contrast_Callback(~,~)
global tcontrast lcontrast updateContrastToggle
prompt = {'High Contrast Percent Limit','Low Contrast Percent Limit'};
dlg_title = 'Contrast limits from 0 to 100';
num_lines = 1;
defaultans = {num2str(tcontrast),num2str(lcontrast)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
tcontrast = str2double(answer{1});
lcontrast = str2double(answer{2});

updateContrastToggle =1;
setSceneAndTime
updateContrastToggle =0;
end

function displayImageFunct(If,channelimg,bkgmedian)
global expDateStr psettings plotSettingsToggle trunccmaplz timeFrames tcontrast lcontrast MainAxes displayTrackingToggle ImageDetails frameToLoad prcntl lprcntl D cmap cmaplz updateContrastToggle


%determine current time Frame
    t = frameToLoad;

%delete old images to keep memory usage low
    axes(MainAxes);
    children = findobj(MainAxes,'Type','image');
    delete(children);

%constrain image axis to be square initially
    MainAxes.Units = 'pixels';
    pos = MainAxes.Position;
    pos(3:4) = [min([pos(3) pos(4)]) min([pos(3) pos(4)])];
    MainAxes.Position = pos;
    MainAxes.Units = 'normalized'; 
    pos = MainAxes.Position;
    pos(2)= 0.5 - pos(4)/2;
    MainAxes.Position = pos;

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
    if updateContrastToggle ==1
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
%         If = bwperim(If);
        dispimg(If>0)=255;
    end

%title the displayed image
    himg = imagesc(uint8(dispimg));
    himgax = get(himg,'Parent');
    himgax.CLim = [0 256];
    ttl = get(himgax,'Title');
    t = ImageDetails.Frame;
    expDateTitleStr = expDateStr;
    [a,~] = regexp(expDateTitleStr,'_');expDateTitleStr(a) = '-';
    set(ttl,'String',[expDateTitleStr ' ' ImageDetails.Scene ' frame ' num2str(t) ' out of ' num2str(timeFrames)]);
    set(ttl,'FontSize',12);
    
    if plotSettingsToggle == 0
        psettings = PlotSettings_callback([],[]);
        plotSettingsToggle=1;
    end
    framesThatMustBeTracked = psettings.framesThatMustBeTracked;

    if ~(t==1)
        if displayTrackingToggle==1
            trajectForPlot = trackingTrajectories(timeFrames);
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
                if ~isempty(h)
                set(h, {'color'}, num2cell(plotcmap,2));
                colormap(cmap);%return colormap so images display properly
                hax = h.Parent;
                hax.Color = 'none';
                himgax.CLim = [0 256];
                himgax.NextPlot = 'replace';
                end
        end
    end
    himgax.YTick = [];
    himgax.XTick = [];
end


%% functions for determining variables

function ImageDetails = InitializeImageDetails

ImageDetails.Scene=[];
ImageDetails.Channel=[];
ImageDetails.Frame=[];

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




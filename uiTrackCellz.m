function uiTrackCellz(FileDate,AutoTrackStr)
global Tracked bkgmedianmat pStruct timeVec expDetailsStruct dirStruct timeSteps DivisionStruct segStruct foStruct xAxisLimits DICimgstack pathStruct segmentimgstack channelimgstack ExportNameKey ExportName plottingTotalOrMedian channelinputs tcontrast lcontrast ThirdPlotAxes SecondPlotAxes togStruct PlotAxes cmap  timeFrames frameToLoad ImageDetails MainAxes SceneList 

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
    togStruct.updateContrast=0;
    togStruct.runIterateToggle =0;
    togStruct.refineTrackingToggle = 1;
    togStruct.plotSettingsToggle =0;
    togStruct.displayTrackingToggle = 0;

    xAxisLimits = [0 50];
    plottingTotalOrMedian = 'median';
    tcontrast = 99;
    lcontrast = 1;
    
    
%determine export details
    ExportNameKey = 'tsi_';
    if ~strcmp(ExportNameKey,'final')
        disp(strcat('Export name key is "',ExportNameKey,'" not FINAL'))
    end
    ExportName = 'fricktrack';
   
    
%initialize global variables  
    channelimgstack =[];
    segmentimgstack =[];
    DICimgstack=[];
    foStruct.cfoName = [];
    foStruct.sfoName =[];
    foStruct.bfoName = [];
    foStruct.nfoName = [];
    foStruct.dfoName=[];
    




%set colormap
    cd(dirStruct.parentdir)
    addpath('Colormaps')
    cmap = gray(256);
%         cmap = colormap(viridis(255));
    % cmap = colormap(magma(255));
    % cmap = colormap(inferno(255));
    % cmap = colormap(plasma(255));
    cmap(256,:)=[1 0 0];

        
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
    segStruct.cell_seg = segInstruct.cell;
    segStruct.background_seg = segInstruct.background;
    channelstoinput = dosestructstruct.channelNames;
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

    timeVec = dosestructstruct.timeVec;
    % timeMatrix(i,:) = timeVec(idxtwo,1:finalFrame)-timeVec(1,stimulationFrame+1); %subtract by time closest to Tgfbeta addition 


%determine how many scenes are present
    dirlist = dir(pathStruct.mstackPath);
    [~,~,~,d] = regexp({dirlist.name},'s[0-9]+');
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
    timeFrames = dim(3);
    bkgmedianmat = zeros(1,timeFrames);
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
    uicontrol('Style', 'pushbutton', 'String', 'Remove area [b]',...
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
        'String','togStruct.displayTrackingToggle [m]',...
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
    
    
    pStruct = loadSegmentParameters([],FileName,dirStruct.exportdir); %loads saved value of pStruct
    timeMat = loadMetadata(FileName,dirStruct.exportdir);
    rtimeMat = round(timeMat);
    stimulationFrame = dosestruct(1).tgfFrame + 1;
    timeVec = rtimeMat(1,:)-rtimeMat(1,stimulationFrame);
    timeSteps = diff(round(timeVec),[],2);
    
    if strcmpi(AutoTrackStr,'AutoTrackCells')
        f.Visible = 'off';
        trackSaveIterateChosen_callback([],[])
        % exportNuclei([],[])
        close(f)
    elseif strcmpi(AutoTrackStr,'AutoExportNuclei')
        exportSegmentedCells([],[])
        close(f)
    elseif strcmpi(AutoTrackStr,'AutoExportTracks')
        f.Visible = 'off';
        exportTrackedCells([],[])
        close(f)
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
global  ImageDetails displaycomments segStruct segStruct
key = get(fig_obj,'CurrentKey');

switch key
    case '1'
        ImageDetails.Channel = segStruct.cell_seg;
        setSceneAndTime
    case '2'
        ImageDetails.Channel = segStruct.nucleus_seg;
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
    case 'b'
        removeArea_Callback([],[]);
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

%% division functions
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

%add cells and link cells
function addareabutton_Callback(~,~) 
 global  ImageDetails Tracked togStruct
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
            zeroImage = zeros(ImageDetails.ImgSize);
            zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
            zerogray = rgb2gray(zeroImage);

            if ~isempty(Tracked{1}.Cellz)
                %if there exists segmenttracking already...then load that. 
                imagio = zeros(ImageDetails.ImgSize);
                imagio(zerogray>0)=1;
                cc = bwconncomp(imagio);
                px = cc.PixelIdxList;
                
                t = ImageDetails.Frame;
                CC = Tracked{t}.Cellz;
                PX = CC.PixelIdxList;        

                idxs = cellfun(@(x) sum(ismember(x,px{1})),PX,'UniformOutput',1);
                index = find(idxs>1);
                if ~isempty(index)
                    newMass = vertcat(PX{index});
                    PX{min(index)} = unique(vertcat(newMass,px{1}));
                    if length(index)>1
                        index(index == min(index))=[];
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
                Tracked{t}.If = IfPerimFunction(CC);
            end
            togStruct.refineTrackingToggle = 1;
            nextbutton_callback([],[]);
        end 
    end
end
%removeArea
function removeArea_Callback(~,~)
global  ImageDetails Tracked togStruct
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
            zeroImage = zeros(ImageDetails.ImgSize);
            zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
            zerogray = rgb2gray(zeroImage);

            if ~isempty(Tracked{1}.Cellz)
                imagio = zeros(ImageDetails.ImgSize);
                imagio(zerogray>0)=1;
                cc = bwconncomp(imagio);
                px = cc.PixelIdxList;

                t = ImageDetails.Frame;
                CC = Tracked{t}.Cellz;
                PX = CC.PixelIdxList;        

                idxs = cellfun(@(x) sum(ismember(x,px{1})),PX,'UniformOutput',1);
                index = find(idxs>1);
                if ~isempty(index)
                    for abc = index
                        oldMass = PX{abc};
                        overlap = ismember(oldMass,px{1});
                        oldMass(overlap)=[];
                        imagio = zeros(ImageDetails.ImgSize);
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
                    CC.PixelIdxList = PX;
                end
                CC.NumObjects = length(CC.PixelIdxList);
                S = regionprops(CC,'Centroid');
                Smat = vertcat(S.Centroid);
                CC.Centroid = Smat;
                Tracked{t}.Cellz = CC;
                Tracked{t}.If = IfPerimFunction(CC);
            end


            togStruct.refineTrackingToggle = 1;
            nextbutton_callback([],[]);
        end 
    end

end


function [cellxx,cellyy,timingkeeper,button] = identifyCellbyClick(t,ginputnum)
            cellxx  =    []; 
            cellyy  =    [];
            timingkeeper = [];
        [xx,yy,button] = ginputmod(ginputnum); %record each click
        if ~isempty(button)
            cellxx  =    xx; 
            cellyy  =    yy;
            timingkeeper = t;
        else
            return
        end
end
function [idx,idxlog,notACell] = findCellByCoordinate(Tracked,ImgSize,cxx,cyy,tkeeper)
%convert coordinates to indices
    cellx = round(cxx);
    celly = round(cyy);
    cellind = sub2ind(ImgSize,celly,cellx);
    
%find the cell
    CC = Tracked{tkeeper}.Cellz;
    PX = CC.PixelIdxList;  
    idxlog = ~cellfun(@(x) sum(x==cellind)>0,PX,'UniformOutput',1);
    idx = find(idxlog==0);
    notACell = isempty(idx);
end
function linkCells_Callback(~,~)
%fast link!
global ImageDetails Tracked togStruct timeFrames

%set accumulating variables equal to empty
cellxx =[];
cellyy =[];
timingkeeper=[];
%choose the initial cell
    togStruct.refineTrackingToggle=0;
    t = ImageDetails.Frame;
    
    notACell = true;
    while notACell %keep running until a cell is clicked
        ginputnum=1;
        [cxx,cyy,tkeeper,button] = identifyCellbyClick(t,ginputnum);
        if isempty(button)
            return
        end
        %find cell
        [idx,idxlog,notACell] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,tkeeper);
    end
    cellxx =[cellxx(:)' cxx];
    cellyy =[cellyy(:)' cyy];
    timingkeeper =[timingkeeper(:)' tkeeper];
    
    a=1;
    while a ==1
        broken = false;
        if button == 1
            for t = tkeeper:timeFrames
            ImageDetails.Frame = t;
            setSceneAndTime
            pause(0.05)
            CC = Tracked{t}.Cellz;
            PX = CC.PixelIdxList;  
            nanidx = isnan(PX{idx});
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
                [idx,idxlog,notACell] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,tkeeper);
            end
            if isempty(button)
                break
            end
            cellxx =[cellxx(:)' cxx];
            cellyy =[cellyy(:)' cyy];
            timingkeeper =[timingkeeper(:)' tkeeper];
            [idx,idxlog] = findCellByCoordinate(Tracked,ImageDetails.ImgSize,cxx,cyy,tkeeper);
        end
        
        if ~broken
            break
        end
    end
        
   
%find the last frame of the intial cell
    %display movie to track to its end?
    togStruct.refineTrackingToggle=1;
    cellx = round(cellxx);
    celly = round(cellyy);
    cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);

    idx = zeros(size(timingkeeper));
    for i = 1:length(timingkeeper)
        t = timingkeeper(i);
        CC = Tracked{t}.Cellz;
        PX = CC.PixelIdxList;  
        idxlog = cellfun(@(x) isempty(find(x==cellind(i),1)),PX,'UniformOutput',1);
        idx(i) = find(idxlog==0);
    end

    %set all preceeding cellIDs to be equal to first cell clicked
    %set all subsequent cellIDs to be equal to last cell clicked
    cellindx = 1:timeFrames;
    cellindx(1:timingkeeper(1))=idx(1);
    cellindx(timingkeeper(end):end) = idx(end);

    %fill in the intermediate cellIDs
    for i=1:length(timingkeeper)-1
        trange=timingkeeper(i):timingkeeper(i+1)-1;
        cellindx(trange) = ones(size(trange)).*idx(i);  
    end


    for i=1:timeFrames
        CC = Tracked{i}.Cellz;
        PX = CC.PixelIdxList;
%         PX{cellindx(1)} = PX{cellindx(i)};
        if cellindx(1)~=cellindx(i)
            PX{cellindx(1)} = PX{cellindx(i)};
            PX{cellindx(i)} = NaN;
            CC.PixelIdxList = PX;
            S = regionprops(CC,'Centroid');
            Smat = vertcat(S.Centroid);
            CC.Centroid = Smat;
            Tracked{i}.Cellz = CC;
            Tracked{i}.If = IfPerimFunction(CC);
        else
            %nothing needed
        end
    end

setSceneAndTime

end


%delete cells
function deletebutton_Callback(~,~) 
    global ImageDetails  Tracked togStruct
    

t = ImageDetails.Frame;
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
     [cellxx,cellyy] = ginput();
       
    cellx = round(cellxx);
    celly = round(cellyy);
    
    cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);
    
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
Tracked{t}.If = IfPerimFunction(CC);
    
togStruct.refineTrackingToggle = 1; 
setSceneAndTime;
   
end
function eliminatebutton_Callback(~,~)

global ImageDetails Tracked togStruct

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
    
    cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);
    
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
Tracked{t}.If = IfPerimFunction(CC);

if button==1
    nextbutton_callback([],[])
elseif button == 3
    prevbutton_callback([],[])
else 
   togStruct.refineTrackingToggle = 1; 
    setSceneAndTime;
end

end
end

%destroy and chosen ones
function destroyAllFramePrevious_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked togStruct

%   determine the frame to load
    t = ImageDetails.Frame;
    CC = Tracked{t}.Cellz;
    PX = CC.PixelIdxList;   

    idxs = false(size(PX));   %choose all cells on frame
    Trackedz = crushThem(Tracked,idxs,[],t); %delete from t=t backwar
    Tracked = Trackedz;
    togStruct.refineTrackingToggle = 1; 
    setSceneAndTime;

end
function destroyAllFrameSubsequent_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked  togStruct

%   determine the frame to load
    t = ImageDetails.Frame;
    CC = Tracked{t}.Cellz;
    PX = CC.PixelIdxList;   

    idxs = false(size(PX));   %choose all cells on frame
    Trackedz = crushThem(Tracked,idxs,t,[]); %delete from t=t onward
    Tracked = Trackedz;
    togStruct.refineTrackingToggle = 1; 
    setSceneAndTime;
end
function destroybuttonAllPrevious_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked togStruct

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
    
      cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
Trackedz = crushThem(Tracked,idxs,[],t);      
Tracked = Trackedz;


    togStruct.refineTrackingToggle = 1; 
    setSceneAndTime;



end
function destroybuttonAllSubsequent_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked  togStruct

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
    
      cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
Trackedz = crushThem(Tracked,idxs,t,[]);      
Tracked = Trackedz;


    togStruct.refineTrackingToggle = 1; 
    setSceneAndTime;



end
function deleteAllonFrame_Callback(~,~)
%delete all cells from one frame
  global  ImageDetails Tracked togStruct
    
    
    
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
Tracked{t}.If = IfPerimFunction(CC);
    
    togStruct.refineTrackingToggle = 1; 
    setSceneAndTime;

end
function destroybutton_Callback(~,~)
%delete a cell from all frames
global ImageDetails Tracked togStruct

t = ImageDetails.Frame;
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
Trackedz = crushThem(Tracked,idxs,[],[]);      
Tracked = Trackedz;


    togStruct.refineTrackingToggle = 1; 
    setSceneAndTime;


end
%choose the cells you want
function chosenOnes_Callback(~,~)
%choose the cells you want
global ImageDetails Tracked togStruct

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
    
      cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);
      
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



    togStruct.refineTrackingToggle = 1; 
    setSceneAndTime;


end
function chosenOnesAllOnFrame_Callback(~,~)
%choose the cells you want
global ImageDetails Tracked togStruct

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


    togStruct.refineTrackingToggle = 1; 
    setSceneAndTime;


end
function erodeOnes_Callback(~,~)
%choose the cells you want
global ImageDetails Tracked togStruct

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
    
      cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);
      
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
        imdub = zeros(ImageDetails.ImgSize);
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
Stacked{i}.If = IfPerimFunction(CC);

end
Trackedz=Stacked;








Tracked = Trackedz;

    togStruct.refineTrackingToggle = 1; 
    setSceneAndTime;


end
function dilateOnes_Callback(~,~)
%choose the cells you want
global ImageDetails Tracked togStruct

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
    
      cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);
      
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
        imdub = zeros(ImageDetails.ImgSize);
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
Stacked{i}.If = IfPerimFunction(CC);

end
Trackedz=Stacked;








Tracked = Trackedz;


    togStruct.refineTrackingToggle = 1; 
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
Stacked{i}.If = IfPerimFunction(CC);
end
Trackedz=Stacked;
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
function PlotCFPnotnorm_callback(~,~)
global togStruct
% tcontrast =99;
% lcontrast =1;
togStruct.cfpNorm = 0;
Plot_callback([],[]);
end
function PlotCFPnorm_callback(~,~)
global togStruct
% tcontrast =99;
% lcontrast =1;
togStruct.cfpNorm = 1;
Plot_callback([],[]);
end
function Plot_callback(~,~)
global segStruct xAxisLimits cmaplz SecondPlotAxes Tracked ImageDetails pathStruct timeFrames frameToLoad PlotAxes togStruct psettings cmap

    if togStruct.plotSettingsToggle == 0
        psettings = PlotSettings_callback([],[]);
        togStruct.plotSettingsToggle=1;
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
    [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,pathStruct.expDirPath,ImageDetails,pathStruct.mstackPath,timeFrames,frameToLoad,PlotAxes,ImageDetails.ImgSize,togStruct.plotSettingsToggle,psettings,makeIMG,makeIMGidx,smooththat);

    % plotStructUI.EGFP = Smad;
    % plotStructUI.mKate = mkate;
    % plotStructUI.Cfp = Cfp;
    % plotStructUI.CfpFC = CfpFC;
    % plotStructUI.SmadFC = SmadFC;
    % plotStructUI.mkateFC = mkateFC;
    % plotStructUI.Smadbkg = Smadbkg;
    % plotStructUI.Cfpbkg = Cfpbkg;
    % plotStructUI.mkatebkg = mkatebkg;

    if strcmp(ImageDetails.Channel,segStruct.cell_seg)
        plotMat = plotStructUI.Smad;
        plotMatFC = plotStructUI.SmadFC;
        ylimit =[0 6];
    elseif strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
%         plotMat = plotStructUI.mkatetotal;
%         plotMatFC = plotStructUI.mkateFCtotal;
        plotMat = plotStructUI.mkate;
        plotMatFC = plotStructUI.mkateFC;
        ylimit =[0 3];
    else
        plotMat = plotStructUI.Smad;
        plotMatFC = plotStructUI.SmadFC;
        ylimit = [0 10];
    end
    
    xmin = 0;
    toplot = plotMatFC;
    idx = true(size(toplot,1),1);
    cmapl = cmaplz;
    idxa = find(idx==1);
    h = plot(SecondPlotAxes,toplot(idx,:)','LineWidth',2);
            
    if togStruct.displayTrackingToggle ==1
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
    
    if togStruct.displayTrackingToggle ==1
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
function plotAxis_callback(~,~)
global xAxisLimits 


prompt = {'xmin','xmax'};
dlg_title = 'set x axis limits...';
inputdlgOutput = inputdlg(prompt,dlg_title);
xAxisLimits = cellfun(@str2num,inputdlgOutput,'UniformOutput',1);
Plot_callback([],[])
end
function Plot_SpecificCell_callback(~,~)
global togStruct cmaplz xAxisLimits ThirdPlotAxes Tracked ImageDetails pathStruct timeFrames frameToLoad PlotAxes psettings

if togStruct.plotSettingsToggle == 0
psettings = PlotSettings_callback([],[]);
togStruct.plotSettingsToggle=1;
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
    cellind = sub2ind(ImageDetails.ImgSize,celly,cellx);
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
[plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,pathStruct.expDirPath,ImageDetails,pathStruct.mstackPath,timeFrames,frameToLoad,PlotAxes,ImageDetails.ImgSize,togStruct.plotSettingsToggle,psettings,makeIMG,makeIMGidx,smooththat);
    
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


    cmapl = cmaplz;
            
    if togStruct.displayTrackingToggle ==1
        for i=1:length(h)
            h(i).Color = cmapl(makeIMGidx(i),:);
        end 
    end

end
function [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,expDirPath,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,ImgSize,plotSettingsToggle,psettings,makeIMG,makeIMGidx,smooththat)
global segStruct pathStruct

PXX = Tracked{1}.Cellz.PixelIdxList;
plotTracesCell = cell(length(PXX),length(Tracked));
    for i = 1:length(Tracked)
        PXX = Tracked{i}.Cellz.PixelIdxList;
        plotTracesCell(:,i) = PXX;
    end
    

cd(pathStruct.mstackPath)
%no bleach correction option yet
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%   open the image files   %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %  open smad img  %
            cd(pathStruct.mstackPath)         
            ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.cell_seg,'*'));
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            cellQ_imgstack = channelfileObject.flatstack;
            
                   %    open cfp img  %
            cd(pathStruct.mstackPath)         
            ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.nucleus_seg,'*'));
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            nuc_imgstack = channelfileObject.flatstack;
            
                    % open background Logical img  %
            cd(pathStruct.segmentPath)         
            ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.background_seg,'*'));
            if length(ff)>1
                ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.background_seg,'*background*'));
            end
    %         ff = dir(strcat(ImageDetails.Channel,'*'));
            filename = char(ff.name);
            channelfileObject = matfile(filename);
            bkglogimgstack = channelfileObject.IfFinal;
            
                    % open nuclear Logical img  %
            cd(pathStruct.segmentPath)         
            ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.nucleus_seg,'*'));
            if length(ff)>1
                ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.nucleus_seg,'*nucleus*'));
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


basalSUB = framesThatMustBeTracked(1)-3;
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
function plotStruct = plotthemfunctionToStructure(Tracked,idScene,pathStruct,timeFrames,makeIMG,makeIMGidx,segStruct) 

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
    

cd(pathStruct.mstackPath)
%no bleach correction option yet
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%   open the image files   %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %  open smad img  %
            cd(pathStruct.mstackPath)         
            ff = dir(strcat('*',idScene,'*',segStruct.cell_seg,'*'));
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
global ExportNameKey ExportName  displaycomments SceneList Tracked pathStruct togStruct psettings

    if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
    end
    
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(pathStruct.trackingPath)
cd ..
    for scenenumber = 1:length(SceneList)
        cd(pathStruct.trackingPath)
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
function xy = labelCells(~,~)


global ExportNameKey ExportName displaycomments Tracked ImageDetails pathStruct frameToLoad togStruct psettings

    if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
    end
    
t = ImageDetails.Frame;
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(pathStruct.trackingPath)
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



function exportTrackedCells(~,~)
global segStruct pathStruct ExportNameKey ExportName dirStruct expDetailsStruct SceneList timeFrames togStruct psettings
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
                    makeIMG = true(size(makeIMG));
                    makeIMGidx = find(makeIMG==1);

                    plotStruct = plotthemfunctionToStructure(trackedArray,idScene,paStruct,tFrames,makeIMG,makeIMGidx,sStruct);
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
        tPath = pathStruct.trackingPath;
        mPath = pathStruct.mstackPath;
        sList = SceneList;
        tFrames = timeFrames;
        cSeg = segStruct.cell_seg;
        nSeg = segStruct.nucleus_seg;
        bSeg = segStruct.background_seg;
        sPath = pathStruct.segmentPath;
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
            cd(dirStruct.exportdir)
            filename = strcat(expDetailsStruct.expDateStr,'_nuclei_export.mat'); 
            save(filename,'exportNucleiStruct');
            

end
function exportNucleiStruct = nucleiToStructure(idScene,mstackPath,timeFrames,cell_seg,nucleus_seg,background_seg,segmentPath)
cd(pathStruct.mstackPath)
%no bleach correction option yet
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%   open the image files   %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %  open smad img  %
        cd(pathStruct.mstackPath)         
        ff = dir(strcat('*',idScene,'*',segStruct.cell_seg,'*'));
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
global ExportNameKey ExportName dirStruct segStruct  pathStruct expDetailsStruct SceneList timeFrames togStruct psettings
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
    cSeg = segStruct.cell_seg;
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
                plotStruct = plotthemfunctionToStructure(Tracked,idScene,pathStruct,timeFrames,makeIMG,makeIMGidx,segStruct); 
                
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
        cd(dirStruct.exportdir)
        filename = strcat(expDetailsStruct.expDateStr,'_tracking_allcells_export.mat'); 
        save(filename,'exportStruct');
end


function exportFrames(~,~)
global ExportNameKey ExportName togStruct ImageDetails pathStruct frameToLoad psettings

    if togStruct.plotSettingsToggle == 0
    psettings = PlotSettings_callback([],[]);
    togStruct.plotSettingsToggle=1;
    end
    
    
%   determine the frame to load
t = ImageDetails.Frame;


framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(pathStruct.trackingPath)
cd ..
CENTROID = struct();
    
        cd(pathStruct.trackingPath)
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
            togStruct.displayTrackingToggle =1;
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
            cd(pathStruct.trackingPath)
            cd ..
            cd ..
            cd('ANNOTATIONS')
            imwrite(X,strcat('s',sceneN,'-t',num2str(i),'.jpg'),'JPEG'); %save image
            cd .. 
            cd('CENTROIDS')
            cd ..
            
%             CENTROIDS.ImageDetails.ImgSize = ImageDetails.ImgSize;
%             CENTROIDS.centroids = xy;
%             CENTROIDS.scene = sceneN;
%             save(strcat('CENTROIDS-',sceneN,'.mat'),'CENTROIDS'); %save CENTROIDS
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        

        
        
        
 
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

Tracked = cell(1,length(timeFrames));
Framez = struct();
Framez.Cellz =[];
Framez.ImageDetails.ImgSize=[];
Framez.filename=[];


for i = 1:length(Tracked)
    Tracked{i} = Framez;
end

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

            togStruct.refineTrackingToggle =1;
            close(h);
setSceneAndTime;
    toc
end

function Tracked = loadTrackedStructure
global pathStruct timeFrames togStruct ImageDetails
if togStruct.runIterateToggle ==0
    cd(pathStruct.trackingPath)
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
        togStruct.refineTrackingToggle=0;
        else
        togStruct.refineTrackingToggle =1;
        end

    else
        Tracked = makeTrackingFile(timeFrames);
    end
else
    Tracked = makeTrackingFile(timeFrames);
end


end
function loadTrackingFile_callback(~,~)
global  Tracked togStruct

Tracked = loadTrackedStructure;
togStruct.refineTrackingToggle =1;
end

%make trajectories for overlay of tracking
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
function Trackedz = trackingCosmetics(Stacked,togStruct)
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
    ImageDetails.ImgSize = Stacked{1}.Cellz.ImageSize;
        for j=1:length(Stacked)
        S = Stacked{j};
        CC = S.Cellz;
        PX = CC.PixelIdxList;
%         PX = Stacked{j}.Cellz.PixelIdxList;
        px = PX(pxidx);
%         CC = Stacked{j}.Cellz;
        CC.Connectivity = connectivity;
        CC.ImageSize = ImageDetails.ImgSize;
        CC.NumObjects = length(px);
        CC.PixelIdxList = px;
        
        if togStruct.displayTrackingToggle==1
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
        S.If = IfPerimFunction(CC);
        Stacked{j} = S;
        end
    Trackedz=Stacked;


end
function trackSaveIterate_callback(~,~)
global togStruct pStruct timeSteps SceneList ImageDetails Tracked pathStruct ExportName timeFrames segStruct

    nucleiDist = pStruct.nucleus.nucDiameter;
    tsteps = timeSteps;
    
togStruct.runIterateToggle =1;
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
            togStruct.refineTrackingToggle =1;
            setSceneAndTime;


        %save
            % saveTrackingFileAs_callback([],[])
            cd(pathStruct.trackingPath)
            filename = 'initial';
            save(strcat(filename,'_',ImageDetails.Scene,'_',ExportName,'.mat'),'Tracked')
%             save(strcat(filename,ExportName,'.mat'),'Tracked')
            
    end
togStruct.runIterateToggle =0;
togStruct.parallel = 0;

end
function trackSaveIterateChosen_callback(~,~)
global SceneList pathStruct ExportName segStruct timeFrames pStruct timeSteps togStruct
    
    %initialize this structure outside of the parfor loop
    psettings = PlotSettings_callback([],[]);
    framesThatMustBeTracked = psettings.framesThatMustBeTracked;
    framesThatMustBeTracked(2) =  min([timeFrames framesThatMustBeTracked(2)+10]);
    nucleiDist = pStruct.nucleus.nucDiameter;
    tsteps = timeSteps;
    
    %remove all global variables before parfor loop
        paStruct = pathStruct;
        sStruct = segStruct;
        tPath = pathStruct.trackingPath;
        sList = SceneList;
        eName = ExportName;
        toggleS = togStruct;
        
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
            Tracked1 = FrickTrackCellsYeah(paStruct,scenestr,sStruct,[],nucleiDist,tsteps);

            %refine tracking if necessary
            Tracked = trackingCosmetics(Tracked1,toggleS);
            
            %display time it takes to track and refine
            disp(scenestr)
            
            cd(tPath)
            filename = 'tsi';
            saveTrackingIteratively(filename,scenestr,eName,Tracked)
            

            %run chosen at two specific frames
            for ftmbt = framesThatMustBeTracked
                frame = ftmbt;
                Tracked1 = specialchosenOnesAllOnFrame_Callback(frame,Tracked);
                Tracked = trackingCosmetics(Tracked1,toggleS);
            end

            cd(tPath)
            filename = 'tsichosen';
            saveTrackingIteratively(filename,scenestr,eName,Tracked)
        end
    end
end
function Tracked = specialchosenOnesAllOnFrame_Callback(frame,Tracked)
    t = frame;
    CC = Tracked{t}.Cellz;
    PX = CC.PixelIdxList;   

    idxs = ~cellfun(@(x) length(x)>1,PX,'UniformOutput',1); %choose all cells on frame

    Trackedz = crushThem(Tracked,~idxs,1,length(Tracked)); 
    Tracked = Trackedz;
end
function saveTrackingIteratively(filename,scenestr,eName,Tracked)
save(strcat(filename,'_',scenestr,'_',eName,'.mat'),'Tracked');
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
save(strcat('commentfinalfricktrack.mat'),'Tracked')
end

%% functions for tracking cells
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
    imgstruct(4).seginstruct = segStruct.cell_seg;
    
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
            segment_Cellnum_array{i} = (1:length(PX))';
            
        end
        
        
%         ArrayStruct.area     =   segment_Area_array;
%         ArrayStruct.ellipt     =   segment_Ellipt_array;
%         ArrayStruct.nucFluor     =   segment_nucFluor_array;
%         ArrayStruct.cellFluor     =   segment_cellFluor_array;
%         ArrayStruct.stdev     =   segment_Stdev_array;
%         ArrayStruct.centroid     =   segment_Centroid_array;
%         ArrayStruct.pixels     =   segment_Pixels_array;
%         ArrayStruct.cellnum     =   segment_Cellnum_array;
%         Frame.Arrayz = ArrayStruct;
%         Tracked{i} = Frame;
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
            centroids = segment_Centroid_array{b};
            centroidsPrev = segment_Centroid_array{a};
            area = segment_Area_array{b};
            areaPrev = segment_Area_array{a};
            nucfluor = segment_nucFluor_array{b};
            nucfluorPrev = segment_nucFluor_array{a};
            cellfluor = segment_cellFluor_array{b};
            cellfluorPrev = segment_cellFluor_array{a};            
            pixels = segment_Pixels_array{b};
            pixelsPrev = segment_Pixels_array{a};

            displacementCutoff = (nucleiDist)*(tsteps(i-1)./10);

            if ~(isempty(pixels) || isempty(pixelsPrev))
               %nearest neighbor distances to determine probability that
               %two cells are the same
                knnnum = 5;
               [distProb,~,distcut] = probSpitter(centroids,centroidsPrev,knnnum,displacementCutoff,[]);
               [areaProbette,~,areacut] = probSpitter(area,areaPrev,size(area,1),[],1.2);
               [nucFluorProbette,~,nuccut] = probSpitter(nucfluor,nucfluorPrev,size(nucfluor,1),[],1.2);
               [cellFluorProbette,~,cellcut] = probSpitter(cellfluor,cellfluorPrev,size(cellfluor,1),[],1.2);
               newProb = distProb.*areaProbette.*nucFluorProbette.*cellFluorProbette;     
               %distProb dimesionas are length(centroidsPrev) x length(centroids)
               

               %centroids [37x2] centroidPrev[34x2] idx[34x3] distProb[34x37];
               trackProb = distProb;
               trackProb(distcut|nuccut)=0;
               [maxvals,idx] = max(trackProb,[],2); %idx is index of input that matches to inputPrev such that input(idx) = inputPrev;
               %idx has length(idx) = length(inputPrev);
               %idx has max(idx) = length(input);
               
               
               %division tracker
               %YFP INCREASES
               %MKATE DECREASES
               %so division is occuring when YFP increases AND MKATE decreases
%                nucdivcut = probSpitterDiv(nucfluor,nucfluorPrev,idx,'lessthan',0.9);
%                celldivcut = probSpitterDiv(cellfluor,cellfluorPrev,idx,'greaterthan',1.1);
%                divcut = nucdivcut&celldivcut;
%                divnumPrev = find(divcut);
%                divnum = idx(divnumPrev);
%                if sum(divcut)>0
%                    stophere=1;
%                end

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


    %now update all the fields
    newpx = cell(size(trackmatrix));
    newpx(:) = {NaN};
    segI = segmentimgstack(:,:,1); %choose the first frame of segmentimgstack to start with
    CC = bwconncomp(segI);
    cnnec = CC.Connectivity;
    imgsize = CC.ImageSize;
    Tracked = cell(1,size(segmentimgstack,3));
    %parfor i = 1:size(segmentimgstack,3)
    for i = 1:size(segmentimgstack,3)
        a=i;
        trackvals=trackmatrix(i,:);
        trackidx = ~isnan(trackvals);
        tracknums = trackvals(trackidx);

        px = segment_Pixels_array{a};
        centroids = segment_Centroid_array{a};
        centnew = nan(length(trackidx),size(centroids,2));
        if ~isempty(px)
            %first pixels
            pxpx = newpx(i,:);
            pxpx(trackidx) = px(tracknums);
            newpx(i,:) = pxpx;

            %then centroids
            centnew(trackidx,:) = centroids(tracknums,:);
        end
        

        

        

        AllCellsPX = newpx(i,:);
        CC = struct();
        Frame = struct();
        CC.ImageSize = imgsize;
        CC.Connectivity = cnnec;
        CC.PixelIdxList = AllCellsPX;
        CC.NumObjects = numel(AllCellsPX);
        CC.Centroid = centnew;
        CC.Division = [];
        Frame.filename = filename;
        Frame.Cellz = CC;
        Frame.If = IfPerimFunction(CC);
%         Frame.Arrayz = ArrayStruct{i};
        Tracked{i} = Frame;
    end
    
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
function cutoffidx = probSpitterDiv(input,inputPrev,idx,greaterORless,relativeCutoff)
valPrev = input(idx);
epsMat = valPrev./inputPrev;
   if strcmpi(greaterORless,'greaterthan')
    cutoffidx = epsMat>relativeCutoff;
   elseif strcmpi(greaterORless,'lessthan')
    cutoffidx = epsMat<relativeCutoff;
   end
end
function [distProb,idx,cutoffidx] = probSpitter(input,inputPrev,knnnum,displacementCutoff,relativeCutoff)
%make a matrix that tells you the probability a value in one vector is the same as anothe

    
    [idx,eps] = knnsearch(input,inputPrev,'K',knnnum); 
    if ~isempty(displacementCutoff)
        eps_scaled = false(size(eps));
    elseif ~isempty(relativeCutoff)
        valPrev = input(idx);
        epsmin = min(inputPrev,valPrev);
        epsmax = max(inputPrev,valPrev);
        eps_scaled = epsmax./epsmin;
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
    
%     distProb(cutoffidx)=0;
    
    
end
function saveTrackingFileAs_callback(~,~)
global  pathStruct Tracked ExportName ImageDetails
cd(pathStruct.trackingPath)
prompt = 'filename of tracking structure to be saved?';
dlg_title = 'save tracking structure as...specific filename';
filename = char(inputdlg(prompt,dlg_title));
save(strcat(filename,'_',ImageDetails.Scene,'_',ExportName,'.mat'),'Tracked')
end

%% Image Display functions
function setSceneAndTime
global bkgmedianmat togStruct DICimgstack foStruct backgroundimgstack segStruct nucleusimgstack segmentimgstack channelimgstack pathStruct frameToLoad ImageDetails Tracked SceneList imgfile 

    %check for empty variables
    cd(pathStruct.mstackPath)
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


    %choose the channel image
        cd(pathStruct.mstackPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.cell_seg,'*'));
        filename = char(ff.name);
        if ~isempty(foStruct.cfoName) %if channelfileObject has been made, check to see if the scene has changed. 
            [a,~] = regexp(foStruct.cfoName,ImageDetails.Scene);
            if isempty(a) %if the scene has changed load the new channelimgstack
                 channelfileObject = matfile(filename);
                 channelimgstack = channelfileObject.flatstack;
                 foStruct.cfoName = char(channelfileObject.Properties.Source);%update foStruct.cfoName
            elseif ~isempty(a) && isempty(channelimgstack)  %if the scene is same but unloaded
                 channelfileObject = matfile(filename);
                 channelimgstack = channelfileObject.flatstack;
                 foStruct.cfoName = char(channelfileObject.Properties.Source);%update foStruct.cfoName
            else
            end
        else %if no foStruct.cfoName, then 
                 channelfileObject = matfile(filename);
                 channelimgstack = channelfileObject.flatstack;
                 foStruct.cfoName = char(channelfileObject.Properties.Source);%update foStruct.cfoName
        end
        cellImg = channelimgstack(:,:,t);
        
    %choose the nucleus image
        cd(pathStruct.mstackPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.nucleus_seg,'*'));
        filename = char(ff.name);
        if ~isempty(foStruct.nfoName) %if channelfileObject has been made, check to see if the scene has changed. 
            [a,~] = regexp(foStruct.nfoName,ImageDetails.Scene);
            if isempty(a) %if the scene has changed load the new channelimgstack
                 nucleusfileObject = matfile(filename);
                 nucleusimgstack = nucleusfileObject.flatstack;
                 foStruct.nfoName = char(nucleusfileObject.Properties.Source);%update foStruct.cfoName
            elseif ~isempty(a) && isempty(nucleusimgstack)  %if the scene is same but unloaded
                 nucleusfileObject = matfile(filename);
                 nucleusimgstack = nucleusfileObject.flatstack;
                 foStruct.nfoName = char(nucleusfileObject.Properties.Source);%update foStruct.cfoName
            else
                %dont do anything
            end
        else %if no foStruct.cfoName, then 
                 nucleusfileObject = matfile(filename);
                 nucleusimgstack = nucleusfileObject.flatstack;
                 foStruct.nfoName = char(nucleusfileObject.Properties.Source);%update foStruct.cfoName
        end
        nucleusImg = nucleusimgstack(:,:,t);
        
        
    %choose the DIC image
        cd(pathStruct.mstackPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*','DIC','*'));
        filename = char(ff.name);
        if ~isempty(foStruct.dfoName) %if channelfileObject has been made, check to see if the scene has changed. 
            [a,~] = regexp(foStruct.dfoName,ImageDetails.Scene);
            if isempty(a) %if the scene has changed load the new channelimgstack
                 DICfileObject = matfile(filename);
                 DICimgstack = DICfileObject.flatstack;
                 foStruct.dfoName = char(DICfileObject.Properties.Source);%update foStruct.cfoName
            elseif ~isempty(a) && isempty(DICimgstack)  %if the scene is same but unloaded
                 DICfileObject = matfile(filename);
                 DICimgstack = DICfileObject.flatstack;
                 foStruct.dfoName = char(DICfileObject.Properties.Source);%update foStruct.cfoName
            else
                %dont do anything
            end
        else %if no foStruct.cfoName, then 
            flist = dir(filename);
            if ~isempty(flist)
                 DICfileObject = matfile(filename);
                 DICimgstack = DICfileObject.flatstack;
                 foStruct.dfoName = char(DICfileObject.Properties.Source);%update foStruct.cfoName
            else
                DICimgstack = false(size(nucleusimgstack));
                foStruct.dfoName = filename;
            end
        end
        DICImg = DICimgstack(:,:,t);
        
         
    %load nucleus segmented image
        cd(pathStruct.segmentPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.nucleus_seg,'*'));
        if length(ff)>1
           ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.nucleus_seg,'*nucleus*'));
        end
        filename = char(ff.name);
        if ~isempty(foStruct.sfoName) %if channelfileObject has been made, check to see if the scene has changed. 
            [a,~] = regexp(foStruct.sfoName,ImageDetails.Scene);
            if isempty(a) %if the scene has changed load the new channelimgstack
                 segmentfileObject = matfile(filename);
                 segmentimgstack = segmentfileObject.IfFinal;
                 foStruct.sfoName = char(segmentfileObject.Properties.Source);%update foStruct.cfoName
            elseif ~isempty(a) && isempty(segmentimgstack)  %if the scene is same but unloaded
                 segmentfileObject = matfile(filename);
                 segmentimgstack = segmentfileObject.IfFinal;
                 foStruct.sfoName = char(segmentfileObject.Properties.Source);%update foStruct.cfoName
            else
                %dont do anything
            end
        else %if no foStruct.cfoName, then 
                 segmentfileObject = matfile(filename);
                 segmentimgstack = segmentfileObject.IfFinal;
                 foStruct.sfoName = char(segmentfileObject.Properties.Source);%update foStruct.cfoName
        end
        segmentimg = segmentimgstack(:,:,t);
        
        
        
    %load background segmented image
        cd(pathStruct.segmentPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.background_seg,'*'));
        if length(ff)>1
            ff = dir(strcat('*',ImageDetails.Scene,'*',segStruct.background_seg,'*background*'));
        end
        filename = char(ff.name);
        if ~isempty(foStruct.bfoName) %if channelfileObject has been made, check to see if the scene has changed. 
            [a,~] = regexp(foStruct.bfoName,ImageDetails.Scene);
            if isempty(a) %if the scene has changed load the new channelimgstack
                 backgroundfileObject = matfile(filename);
                 backgroundimgstack = backgroundfileObject.IfFinal;
                 foStruct.bfoName = char(backgroundfileObject.Properties.Source);%update foStruct.cfoName
            elseif ~isempty(a) && isempty(backgroundimgstack)  %if the scene is same but unloaded
                 backgroundfileObject = matfile(filename);
                 backgroundimgstack = backgroundfileObject.IfFinal;
                 foStruct.bfoName = char(backgroundfileObject.Properties.Source);%update foStruct.cfoName
            else
                %dont do anything
            end
        else %if no foStruct.cfoName, then 
                 backgroundfileObject = matfile(filename);
                 backgroundimgstack = backgroundfileObject.IfFinal;
                 foStruct.bfoName = char(backgroundfileObject.Properties.Source);%update foStruct.cfoName
                bkgmedianmat = zeros(1,size(backgroundimgstack,3));
                for k = 1:size(backgroundimgstack,3)
                    backgroundimg = backgroundimgstack(:,:,k);
                    segimg = segmentimgstack(:,:,k);
                    chanimg = channelimgstack(:,:,k);
                    backgroundimg(segimg) = true;
                    bkgpixels = chanimg(backgroundimg);
                    bkgmedianmat(k) = nanmedian(bkgpixels);
                end
                clear segimg
                clear chanimg
        end
        backgroundimg = backgroundimgstack(:,:,t);
         
        
%choose the segmentation image
    if isempty(Tracked{1}.Cellz) %Tracked is empty try to load tracking or make new tracking structure
        togStruct.refineTrackingToggle = 0;
        cd(pathStruct.segmentPath)
        imgfile = dir(strcat('*',ImageDetails.Scene,'*',segStruct.nucleus_seg,'*'));
        If = segmentimg;
        Tracked = loadTrackedStructure;
        centroidDisagreement=false;
        Ts = Tracked{t};
        CC = Ts.Cellz;
    else  %if there exists segmenttracking already...then use that. 
        Ts = Tracked{t};
        CC = Ts.Cellz;
        PX = CC.PixelIdxList;
        Centroids  = CC.Centroid;
        centroidDisagreement = ~(size(Centroids,2)==CC.NumObjects);
    end
    
    fnames = fieldnames(Ts);
    if sum(strcmp(fnames,'If')>0)
        If = Tracked{t}.If;
    else
        If =  IfPerimFunction(CC);
        Tracked{t}.If = If;
    end



    if togStruct.displayTrackingToggle==1
        if centroidDisagreement      
            togStruct.refineTrackingToggle = 1;
        end
    end
    % disp(togStruct.refineTrackingToggle)
    if togStruct.refineTrackingToggle == 1
        Tracked = trackingCosmetics(Tracked,togStruct);
    end
    togStruct.refineTrackingToggle=0;


    if strcmp(ImageDetails.Channel,segStruct.nucleus_seg)
        channelimg = nucleusImg;
    elseif strcmp(ImageDetails.Channel,segStruct.cell_seg)
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


    bkgmedian = bkgmedianmat(t);
    if togStruct.runIterateToggle ==0
        displayImageFunct(If,channelimg,bkgmedian);
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


if isempty(Tracked{1}.Cellz)
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
function If = IfPerimFunction(CC)
        PX = CC.PixelIdxList;
        makeIMG = cellfun(@(x) length(x)<2,PX,'UniformOutput',1); %choose only the cells without NAN
        CC.PixelIdxList = PX(~makeIMG);
        CC.NumObjects = length(PX(~makeIMG));

        segmentimgL = labelmatrix(CC);
        segmentimgz = false(size(segmentimgL));
        segmentimgz(segmentimgL>0)=1;
        Ifz = segmentimgz;
        If = bwperim(Ifz) | bwperim(imdilate(Ifz,strel('disk',1)));
end
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
        lxy(i) = size(xy{i},1);
    end
    
    traject = nan(max(lxy),2,t);

    for i = 1:t
        traject(1:lxy(i),1:2,i) = xy{i};
    end
end
function displayImageFunct(If,channelimg,bkgmedian)
global oldscenename expDetailsStruct togStruct timeFrames tcontrast lcontrast MainAxes ImageDetails prcntl lprcntl cmap cmaplz

scenestr = ImageDetails.Scene;
tnum = ImageDetails.Frame;

%delete old images to keep memory usage low
    children1 = findobj(MainAxes,'Type','image');
    children2 = findobj(MainAxes,'Type','Line');
    delete(children1);
    delete(children2);
    set(MainAxes,'NextPlot','replace')
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

% %update contrast if contrast values are updated
%     if togStruct.contrastAdjust ==1
%         togStruct.changeSceneOrChannel = 1;
%         oldscenename = ImageDetails.Channel;
%     end

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
        end
        
        dispimg = disprgb;
        
    else  %under normal circumstances
        if togStruct.changeSceneOrChannel ==1
            prcntl = prctile(channelimg(:),tcontrast);
            togStruct.changeSceneOrChannel =0;
        end
        lprcntl = bkgmedian.*0.90;
        scaleFactor = 255./(prcntl - lprcntl);
        dispimg = channelimg.*scaleFactor;
        dispimg = dispimg-(lprcntl.*scaleFactor);
        dispimg(dispimg> 255) =254;
        dispimg(If>0)=255;
    end

%title the displayed image
    image(MainAxes,uint8(dispimg));
    colormap(MainAxes,cmap);

    expDetailsStruct.expDateTitleStr = expDetailsStruct.expDateStr;
    [a,~] = regexp(expDetailsStruct.expDateTitleStr,'_');expDetailsStruct.expDateTitleStr(a) = '-';

    ttl = get(MainAxes,'Title');
    set(ttl,'String',[expDetailsStruct.expDateTitleStr ' ' scenestr ' frame ' num2str(tnum) ' out of ' num2str(timeFrames)]);
    set(ttl,'FontSize',12);
    set(MainAxes,'NextPlot','add')

    %plot cell tracking "tails"
    if ~(tnum==1)
        if togStruct.displayTrackingToggle==1
            trajectForPlot = trackingTrajectories(timeFrames);
            set(MainAxes,'NextPlot','add');
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
                cmaplz = colorcubemodified(length(idx),'colorcube');
                cmapl = cmaplz;
                plotcmap = zeros(length(idxa),3);
                for i=1:length(idxa)
                    plotcmap(i,:) =  cmapl(idxa(i),:);
                end

                x = mainplotX(idx,:);
                y = mainplotY(idx,:);
                t0 =max([1 tnum-20]);
                h = plot(MainAxes,x(:,t0:tnum)',y(:,t0:tnum)','LineWidth',2);
                set(h, {'color'}, num2cell(plotcmap,2));
            end
        end
    end
    set(MainAxes,'Color','none');
    set(MainAxes,'CLim',[0 256]);
    set(MainAxes,'YTick',[]);
    set(MainAxes,'XTick',[]);

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

    cmaplz = ccnew;
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
updateCrossHair(fig,initialState.CrossHair)
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


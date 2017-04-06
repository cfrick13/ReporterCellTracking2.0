function uiTrackCellz
global DICimgstack dfoName cfoName trackingPath background_seg bfoName backgroundimgstack sfoName cell_seg nucleus_seg segmentimgstack channelimgstack segmentPath mstackPath runIterate ExportNameKey ExportName exportdir plottingTotalOrMedian channelinputs adjuster cmapper tcontrast lcontrast ThirdPlotAxes SecondPlotAxes OGExpDate plottingON PlotAxes cmap TC A AA timeFrames frameToLoad ImageDetails MainAxes SceneList displaytracking imgsize ExpDate
adjuster=0;
plottingTotalOrMedian = 'median';
tcontrast = 99;
lcontrast = 1;
% exportdir = 'C:\Users\Kibeom\Desktop\Tracking\Export\';
ExportNameKey = 'final';
% ExportNameKey = 'tsichosen';
if strcmp(ExportNameKey,'final')
else
disp(strcat('Export name key is "',ExportNameKey,'" not FINAL'))
end
ExportName = 'fricktrack';
cfoName = [];
channelimgstack =[];
segmentimgstack =[];
sfoName =[];
bfoName = [];

DICimgstack=[];
dfoName=[];


runIterate =0;
TC = 1;
ImageDetails = InitializeImageDetails;
displaytracking = 0;

%%% set colormap for the images %%%
cmap = colormap(gray(255));
% cmap = colormap(magma(255));
% cmap = colormap(inferno(255));
% cmap = colormap(plasma(255));
cmap(255,:)=[1 0 0];
cmapper = cmap;
close all


plottingON =0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose directory of experiment to track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
        if isempty(b)
            [~,b] = regexp(mdir,'\');
        end
parentdir = mdir(1:b(end));
exportdir = strcat(parentdir,'Export/');

cd(parentdir)
cd ..
A = uigetdir;
AA = 'D:\Users\zeiss\Documents\MATLAB';
cd(A)
mstackName = 'flat mstack';
experimentdir = A;
mstackPath = strcat(experimentdir,'/',mstackName);
segmentName = 'segment mstack';
segmentPath = strcat(experimentdir,'/',segmentName);
trackName = 'tracking files';
trackingPath = strcat(experimentdir,'/',trackName);

cd(experimentdir)
dirlist = dir(trackName);
if isempty(dirlist)
    mkdir(trackName);
end

%determine date of experiment
[a,b] = regexp(A,'201[0-9]');
[c,d] = regexp(A,'exp[0-9]');
ExpDate = A(a:b+6);OGExpDate = A(a:d); [a,~] = regexp(ExpDate,'_');ExpDate(a) = '-';


%subdirectories should include
%> [ flat mstack ]
    %> [ mstack images ]
cd(mstackPath)


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% renamemCherrytoMkate(A,B)
% renamemWRONGtoRIGHT(A,B)
% cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelstoinput = dosestructstruct.channelNames;
channelinputs =channelregexpmaker(channelstoinput);
bkg = dosestructstruct.BACKGROUND;
imgsize = dosestructstruct.dimensions;
backgroundScenes = dosestructstruct.indicesChoice;

BACKGROUND = bkg{1};
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
% f.Visible ='off';
f.Units = 'pixels';
f.Position =[10,10,1800,1200];


buttonwidth = 80;
buttonheight = 50;

ypositions = sort([100:20:1000],'descend');
xpositions = ones(1,length(ypositions)).*1600;

        mmm=1;
htexttwo = uicontrol('Style','text','String','To choose channel push 1, 2, or 3',...
          'Position',[xpositions(mmm)-buttonwidth,ypositions(mmm),buttonwidth+buttonwidth,buttonheight]);
        mmm=mmm+1;
        mmm=mmm+1;
        mmm=mmm+1;

hNextFrame = uicontrol('Style','pushbutton',...
    'String','NextFrame [f]',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@nextbutton_callback);
hPreviousFrame = uicontrol('Style','pushbutton',...
    'String','Previous frame [a]',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@prevbutton_callback);
        mmm=mmm+1;
        mmm=mmm+1;
        
hGoToFrame = uicontrol('Style','pushbutton',...
    'String','Go to Frame',...
    'Position',[xpositions(mmm)-120,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@gotobutton_callback);

        
hFinalFrame = uicontrol('Style','pushbutton',...
    'String','FinalFrame [g]',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@finalbutton_callback);
hFirstFrame = uicontrol('Style','pushbutton',...
    'String','First frame [z]',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@firstbutton_callback);
        mmm=mmm+1;
        mmm=mmm+1;

htextone = uicontrol('Style','text','String','Choose Scene',...
    'Position',[xpositions(mmm)-buttonwidth,ypositions(mmm)-buttonheight./2,buttonwidth+buttonwidth,buttonheight]);
        mmm=mmm+1;
hpopup = uicontrol('Style','popupmenu',...
    'String',SceneList',...
    'Position',[xpositions(mmm)-buttonwidth./2,ypositions(mmm),buttonwidth.*1.5,buttonheight./2],...
    'Callback',@popup_menu_Callback);
        mmm=mmm+1;
        mmm=mmm+1;
        mmm=mmm+1;
hAddArea = uicontrol('Style','pushbutton','String','AddArea [v]',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@addareabutton_Callback);
hLinkCells = uicontrol('Style','pushbutton','String','LinkCells [r]',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@linkCells_Callback);
        mmm=mmm+1;
        mmm=mmm+1;
hDelete = uicontrol('Style','pushbutton','String','Delete [d]',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@deletebutton_Callback);
hDelete = uicontrol('Style','pushbutton','String','DeleteAllOnFrame',...
    'Position',[xpositions(mmm)-120,ypositions(mmm)+10,buttonwidth,buttonheight],...
    'Callback',@deleteAllonFrame_Callback);


hEliminate = uicontrol('Style','pushbutton','String','Eliminate [e]',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@eliminatebutton_Callback);
        mmm=mmm+1;
        mmm=mmm+1;
        
        
hDestroy = uicontrol('Style','pushbutton','String','Destroy',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@destroybutton_Callback);
hDestroy = uicontrol('Style','pushbutton','String','DestroyPrevious',...
    'Position',[xpositions(mmm)-120,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@destroybuttonAllPrevious_Callback);
hDestroy = uicontrol('Style','pushbutton','String','DestroySubsequent',...
    'Position',[xpositions(mmm)-120,ypositions(mmm+2),buttonwidth,buttonheight],...
    'Callback',@destroybuttonAllSubsequent_Callback);
hchosenOnes = uicontrol('Style','pushbutton','String','Chosen Ones',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@chosenOnes_Callback);
hchosenOnes = uicontrol('Style','pushbutton','String','Chosen OnesAllOnFrame',...
    'Position',[xpositions(mmm)+120,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@chosenOnesAllOnFrame_Callback);
       mmm=mmm+1;
       mmm=mmm+1;
       
       
hRemoveArea = uicontrol('Style', 'pushbutton', 'String', 'Remove area',...
    'Position',[xpositions(mmm),ypositions(mmm),buttonwidth,buttonheight/1.5],...
    'Callback',@removeArea_Callback);
        mmm= mmm+1;
hchosenOnesEnd = uicontrol('Style','pushbutton','String','Chosen Ones EndOnly',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*2,buttonheight./2],...
    'Callback',@chosenOnesEnd_Callback);
       mmm=mmm+1;

        
hErode = uicontrol('Style','pushbutton','String','Erode Chosen Only',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*2,buttonheight./2],...
    'Callback',@erodeOnes_Callback);
       mmm=mmm+1;


hDilate = uicontrol('Style','pushbutton','String','Dilate Chosen Only',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*2,buttonheight./2],...
    'Callback',@dilateOnes_Callback);
       mmm=mmm+1;
        mmm=mmm+1;
        mmm=mmm+1;

        
        
hDisplayTracking = uicontrol('Style','pushbutton',...
    'String','DisplayTracking [m]',...
    'Position',[xpositions(mmm)-(buttonwidth./2),ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@displayTrackingButton_Callback);
hTrack = uicontrol('Style','pushbutton',...
    'String','Run Tracking [t]',...
    'Position',[xpositions(mmm)+(buttonwidth./2),ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@trackbutton_Callback);
        mmm=mmm+1;
        mmm=mmm+1;
%         mmm=mmm+1;
        mmm=mmm+1;
 hContrast = uicontrol('Style','pushbutton',...
    'String','contrast user',...
    'Position',[xpositions(mmm),ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@contrast_Callback);       
        
        mmm=mmm+1;
        mmm=mmm+1;      
        
        
        %%%%
        %%%%
hSaveTrackingAs = uicontrol('Style','pushbutton',...
    'String','SaveTrackingAs',...
    'Position',[xpositions(mmm)-buttonwidth./2,ypositions(mmm)-buttonheight,buttonwidth.*2,buttonheight.*2],...
    'Callback',@saveTrackingFileAs_callback);
hTrackSaveIterate = uicontrol('Style','pushbutton',...
    'String','trackSaveIterate',...
    'Position',[xpositions(mmm)-buttonwidth./2 - 100,ypositions(mmm-1)-buttonheight,buttonwidth,buttonheight],...
    'Callback',@trackSaveIterate_callback);
hTrackSaveIterate = uicontrol('Style','pushbutton',...
    'String','TSIchosen',...
    'Position',[xpositions(mmm)-buttonwidth./2 - 100,ypositions(mmm+2)-buttonheight,buttonwidth,buttonheight],...
    'Callback',@trackSaveIterateChosen_callback);
        mmm=mmm+1;
        mmm=mmm+1;
        mmm=mmm+1;
        mmm=mmm+1;
        mmm=mmm+1;
hLoadTracking = uicontrol('Style','pushbutton',...
    'String','LoadTracking',...
    'Position',[xpositions(mmm),ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@loadTrackingFile_callback);
        mmm=mmm+1;
        mmm=mmm+1;
        mmm=mmm+1;
        
hPlot = uicontrol('Style','pushbutton',...
    'String','PLOT!',...
    'Position',[xpositions(mmm)-0,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@Plot_callback);

hPlotCFPnorm = uicontrol('Style','pushbutton',...
    'String','plotCFPnorm?',...
    'Position',[xpositions(mmm)-90,ypositions(mmm)+20,buttonwidth,buttonheight./2],...
    'Callback',@PlotCFPnorm_callback);

hPlotCFPnotnorm = uicontrol('Style','pushbutton',...
    'String','plotCFPnotnorm?',...
    'Position',[xpositions(mmm)-90,ypositions(mmm),buttonwidth,buttonheight./2],...
    'Callback',@PlotCFPnotnorm_callback);


hPlotSpecificCell = uicontrol('Style','pushbutton',...
    'String','Plot Specific Cell!',...
    'Position',[xpositions(mmm)+80,ypositions(mmm)+buttonheight./2,buttonwidth,buttonheight./2],...
    'Callback',@Plot_SpecificCell_callback);
hPlotSettings = uicontrol('Style','pushbutton',...
    'String','Plot Settings!',...
    'Position',[xpositions(mmm)+80,ypositions(mmm),buttonwidth,buttonheight./2],...
    'Callback',@PlotSettings_callback);
      mmm=mmm+1;                      
      mmm=mmm+1;  
      
hExportCells = uicontrol('Style','pushbutton',...
    'String','ExportTrackedCells',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*1.5,buttonheight./2],...
    'Callback',@exportTrackedCells);
hExportCells = uicontrol('Style','pushbutton',...
    'String','ExportAllCells',...
    'Position',[xpositions(mmm)+80,ypositions(mmm),buttonwidth.*1.2,buttonheight./2],...
    'Callback',@exportAllCells);
        mmm=mmm+1; 
        mmm=mmm+1; 
        
hLabelCells = uicontrol('Style','pushbutton',...
    'String','Label Cells',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*1.5,buttonheight./2],...
    'Callback',@labelCells);
        mmm=mmm+1; 
%         mmm=mmm+1; 
        
        
hcomment = uicontrol('Style','pushbutton','String','Comments',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*2,buttonheight./2],...
    'Callback',@comment_Callback);
       mmm=mmm+1;
        mmm=mmm+1;
        mmm=mmm+1;
        
hExportLabelsCells = uicontrol('Style','pushbutton',...
    'String','ExportLabels',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@ exportLabels);
%     'Callback',@ exportFrames);
        mmm=mmm+1; 
        
       
sldc1 = uicontrol('Style', 'slider',...
        'String','channel1',...
        'Min',0,'Max',255,'Value',200,...
        'Position', [xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight./2],...
        'Callback', @PMthreshslider); mmm=mmm+2;
        
       
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
imgdim = 512.*1.8;
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
f.Position = [400 150 1150 800];
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
        deletebutton_Callback([],[]);
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
    case '0'
        displaycomments=1;
    xy = getxy([],[]);
    [~,comments,commentpos,cellidx]=updatecomments(xy);
    setcommentsTracking(comments,commentpos)
    dispxy(xy)
end

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
global   ImageDetails SceneList A trackPath




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
global   ImageDetails SceneList A trackPath

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
global toggleCFPnorm SecondPlotAxes Tracked ImageDetails A trackingPath timeFrames mstackPath frameToLoad PlotAxes imgsize plottingON psettings cmaplz displaytracking cmap

if plottingON == 0
psettings = PlotSettings_callback([],[]);
plottingON=1;
end
framesThatMustBeTracked = psettings.framesThatMustBeTracked;

% for jy = 1:length(framesThatMustBeTracked)
% PX = Tracked{framesThatMustBeTracked(jy)}.Cellz.PixelIdxList;
% %     makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)==1,PX,'UniformOutput',1)); %choose only the cells without NAN
%     makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
% end
% [comments,commentpos,cellidx,plotidx]=commentsforplot(Tracked);
% makeIMG = makeIMG(1,:)& makeIMG(2,:); 

for jy = 1
PX = Tracked{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
%     makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)==1,PX,'UniformOutput',1)); %choose only the cells without NAN
    makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
end
[comments,commentpos,cellidx,plotidx]=commentsforplot(Tracked);
makeIMGidx = find(makeIMG==1);
if ~isempty(plotidx)
iidd = find(~ismember(makeIMGidx,plotidx));
else
    iidd=[];
end
% makeIMG(makeIMGidx(iidd))=0;
% makeIMGidx = find(makeIMG==1);


smooththat=0;
[plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,A,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
% plotStruct = plotthemfunctionToStructure(framesThatMustBeTracked,Tracked,A,ImageDetails,trackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
smooththat=toggleCFPnorm;


% plotStructUI.EGFP = Smad;
% plotStructUI.mKate = mkate;
% plotStructUI.Cfp = Cfp;
% plotStructUI.CfpFC = CfpFC;
% plotStructUI.SmadFC = SmadFC;
% plotStructUI.mkateFC = mkateFC;
% plotStructUI.Smadbkg = Smadbkg;
% plotStructUI.Cfpbkg = Cfpbkg;
% plotStructUI.mkatebkg = mkatebkg;

if strcmp(ImageDetails.Channel,'EGFP')
    plotMat = plotStructUI.Smad;
    plotMatFC = plotStructUI.SmadFC;
elseif strcmp(ImageDetails.Channel,'mKate')
    plotMat = plotStructUI.mkate;
    plotMatFC = plotStructUI.mkateFC;
else
    plotMat = plotStructUI.Smad;
    plotMatFC = plotStructUI.SmadFC;
end
    

xmin = 0;
toplot = plotMatFC;
            idx = true(size(toplot,1),1);
            cmapl = cmaplz;
            idxa = find(idx==1);
h = plot(SecondPlotAxes,toplot(idx,:)');
            
    if displaytracking ==1
        for i=1:length(h)
            h(i).Color = cmapl(idxa(i),:);
        end
    colormap(cmap);    
    end
        for po = iidd  %fade out the commented cells
        h(po).LineStyle = ':';
        end

SecondPlotAxes.XLim = ([xmin 40]);
SecondPlotAxes.YLim = ([0 6]);


toplot = plotMat;
            idx = true(size(toplot,1),1);
            cmapl = cmaplz;
h = plot(PlotAxes,toplot(idx,:)');
    if displaytracking ==1
        for i=1:length(h)
            h(i).Color = cmapl(i,:);
        end
    colormap(cmap);    
    end
        for po = iidd  %fade out the commented cells
        h(po).LineStyle = ':';
        end
        

PlotAxes.XLim = ([xmin size(toplot,2)]);
PlotAxes.XLim = ([xmin 40]);
PlotAxes.YLim = ([0 prctile(toplot(:),98)]);
PlotAxes.YLim = ([prctile(toplot(:),1)./1.2 prctile(toplot(:),98)]);
% PlotAxes.YLim = ([0 1.4e06]);
% PlotAxes.YScale = 'log';


end


function Plot_SpecificCell_callback(~,~)
global plottingTotalOrMedian ThirdPlotAxes toggleCFPnorm Tracked ImageDetails A mstackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings

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
[plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,A,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
smooththat=toggleCFPnorm;
SmadFC = plotStructUI.SmadFC;
if smooththat==1
% toplot = SmadFC./CfpFC;
toplot = SmadFC;
else
toplot = SmadFC;    
end
h = plot(ThirdPlotAxes,toplot','LineWidth',3);
ThirdPlotAxes.XLim = ([0 40]);
ThirdPlotAxes.YLim = ([0 6]);

end

function [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,A,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat)
global plottingTotalOrMedian cell_seg nucleus_seg background_seg segmentPath


plotTracesCell = cell(length(makeIMGidx),length(Tracked));
    for i = 1:length(Tracked)
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

if strcmpi(plottingTotalOrMedian,'total')
    Smad = cellfun(@nansum,cellQ_pxls,'UniformOutput',1);
    % Smad = cellfun(@nanmean,cellQ_pxls,'UniformOutput',1);
    Smad(Smad==single(13579)) = NaN;
    Cfp = cellfun(@nansum,nuc_pxls,'UniformOutput',1);
    Cfp(Cfp==single(13579)) = NaN;
    mkate = cellfun(@nansum,mkatepxls,'UniformOutput',1);
    mkate(mkate==single(13579)) = NaN;
elseif strcmpi(plottingTotalOrMedian,'median')
    Smad = cellfun(@nanmedian,cellQ_pxls,'UniformOutput',1);
    % Smad = cellfun(@nanmean,cellQ_pxls,'UniformOutput',1);
    Smad(Smad==single(13579)) = NaN;
    Cfp = cellfun(@nanmedian,nuc_pxls,'UniformOutput',1);
    Cfp(Cfp==single(13579)) = NaN;
    mkate = cellfun(@nanmedian,mkatepxls,'UniformOutput',1);
    mkate(mkate==single(13579)) = NaN;
end

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
    
    
%     Smad,Cfp,mkate,CfpFC,SmadFC,mkateFC,Smadbkg,Cfpbkg,mkatebkg
plotStructUI.Smad = Smad;
plotStructUI.mkate = mkate;
plotStructUI.Cfp = Cfp;
plotStructUI.CfpFC = CfpFC;
plotStructUI.SmadFC = SmadFC;
plotStructUI.mkateFC = mkateFC;
plotStructUI.Smadbkg = Smadbkg;
plotStructUI.Cfpbkg = Cfpbkg;
plotStructUI.mkatebkg = mkatebkg;
    
end
function plotStruct = plotthemfunctionToStructure(framesThatMustBeTracked,Tracked,A,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat)
global plottingTotalOrMedian cell_seg nucleus_seg background_seg segmentPath
plotStruct = struct();



plotTracesCell = cell(length(makeIMGidx),length(Tracked));
    for i = 1:length(Tracked)
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

%determine median pxl intensities
    Smad = cellfun(@nanmedian,cellQ_pxls,'UniformOutput',1);
        Smad(Smad==single(13579)) = NaN;
    Cfp = cellfun(@nanmedian,nuc_pxls,'UniformOutput',1);
        Cfp(Cfp==single(13579)) = NaN;
    mkate = cellfun(@nanmedian,mkatepxls,'UniformOutput',1);
        mkate(mkate==single(13579)) = NaN;
for i = 1:size(cellQ_pxls,1)
    plotStruct(i).medianNucEGFP = Smad(i,:);
    plotStruct(i).medianNucRFP = Cfp(i,:);
end

%determine total pxl intensities
    Smad = cellfun(@nansum,cellQ_pxls,'UniformOutput',1);
        Smad(Smad==single(13579)) = NaN;
    Cfp = cellfun(@nansum,nuc_pxls,'UniformOutput',1);
        Cfp(Cfp==single(13579)) = NaN;
    mkate = cellfun(@nansum,mkatepxls,'UniformOutput',1);
        mkate(mkate==single(13579)) = NaN;
for i = 1:size(cellQ_pxls,1)
    plotStruct(i).totalNucEGFP = Smad(i,:);
    plotStruct(i).totalNucRFP = Cfp(i,:);
end

%determine mean pxl intensities
    Smad = cellfun(@nanmean,cellQ_pxls,'UniformOutput',1);
    Smad(Smad==single(13579)) = NaN;
    Cfp = cellfun(@nanmean,nuc_pxls,'UniformOutput',1);
    Cfp(Cfp==single(13579)) = NaN;
    mkate = cellfun(@nanmean,mkatepxls,'UniformOutput',1);
    mkate(mkate==single(13579)) = NaN;
for i = 1:size(cellQ_pxls,1)
    plotStruct(i).meanNucEGFP = Smad(i,:);
    plotStruct(i).meanNucRFP = Cfp(i,:);
    plotStruct(i).medianCfpbkg = Cfpbkg;
    plotStruct(i).medianSmadbkg = Smadbkg;
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
global ExportNameKey ExportName OGExpDate displaycomments SceneList Tracked ImageDetails A trackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings

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
global ExportNameKey ExportName exportdir mstackPath OGExpDate SceneList Tracked ImageDetails A trackingPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings
exportStruct = struct();

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(trackingPath)
cd ..
    for scenenumber = 1:length(SceneList)
        cd(trackingPath)
        sceneN = SceneList{scenenumber};
        disp(sceneN)
        ImageDetails.Scene = sceneN;
      

        trackfile = dir(strcat(ExportNameKey,'*',sceneN,'*',ExportName,'.mat'));
%         trackfile = dir('finalfricktrack.mat');
        trackfilename = char({trackfile.name});
        
            if ~isempty(trackfilename)
                load(trackfilename)
                PX = Tracked{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
                makeIMG = false(length(framesThatMustBeTracked),length(PX));
                   
                    for jy = 1:length(framesThatMustBeTracked)
                    PX = Tracked{framesThatMustBeTracked(jy)}.Cellz.PixelIdxList;
%                     makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)==1,PX,'UniformOutput',1)); %choose only the cells without NAN
                    makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
                    end
                    
                makeIMG = makeIMG(1,:)&makeIMG(2,:);
                makeIMGidx = find(makeIMG==1);
%                 makeIMG = true(size(makeIMG));
                smooththat=0;
%                 [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,A,ImageDetails,SceneDirPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
                plotStruct = plotthemfunctionToStructure(framesThatMustBeTracked,Tracked,A,ImageDetails,mstackPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
                
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
                
                
%                 exportStruct(:).scene = sceneN;
        
stophere=1;
%         xlswrite(filename,eggcell,sceneN);
                end
    end
        cd(exportdir)
        filename = strcat(OGExpDate,'_tracking_export.mat'); 
        save(filename,'exportStruct');
end

function exportAllCells(~,~)
global ExportNameKey ExportName exportdir OGExpDate SceneList Tracked ImageDetails A trackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings
exportStruct = struct();

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
                PX = Tracked{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
                makeIMG = zeros(length(framesThatMustBeTracked),length(PX));
%                    
%                     for jy = 1:length(framesThatMustBeTracked)
%                     PX = Tracked{framesThatMustBeTracked(jy)}.Cellz.PixelIdxList;
% %                     makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)==1,PX,'UniformOutput',1)); %choose only the cells without NAN
%                     makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
%                     end
%                     
                makeIMG = makeIMG(1,:)& makeIMG(2,:);
                makeIMG = true(size(makeIMG));
                makeIMGidx = find(makeIMG==1);
                smooththat=0;
%                 [plotStructUI] = plotthemfunction(framesThatMustBeTracked,Tracked,A,ImageDetails,SceneDirPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
                plotStruct = plotthemfunctionToStructure(framesThatMustBeTracked,Tracked,A,ImageDetails,SceneDirPath,timeFrames,frameToLoad,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
                
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
                
                
%                 exportStruct(:).scene = sceneN;
        
stophere=1;
%         xlswrite(filename,eggcell,sceneN);
                end
    end
        cd(exportdir)
        filename = strcat(OGExpDate,'_tracking_export.mat'); 
        save(filename,'exportStruct');
end

function xy = labelCells(~,~)


global ExportNameKey ExportName displaycomments  SceneList Tracked ImageDetails A trackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings

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
global ExportNameKey ExportName  imgsize displaytracking SceneList Tracked ImageDetails A trackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings

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
global ExportNameKey ExportName  imgsize displaytracking SceneList Tracked ImageDetails A trackPath timeFrames frameToLoad PlotAxes imgsize plottingON psettings adjuster

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
function traject = trackingTrajectories(frameToLoad,ImageDetails)
global Tracked imgsize

%   determine the frame to load
t = ImageDetails.Frame;

xy = cell(1,t);
lxy = zeros(1,t);
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;
% makeIMG = cellfun(@(x) length(x)==1,PX,'UniformOutput',1); %choose only the cells without NAN
makeIMG = cellfun(@(x) length(x)<2,PX,'UniformOutput',1); %choose only the cells without NAN

    for i = 1:t %determine centroids from 1:t for plotting a tracking tail 
        CC = Tracked{i}.Cellz;
        %%%%%%%%%%%%%%%%%%%%%
        % S = regionprops(CC,'Centroid');
        % xy{i} = vertcat(S(:).Centroid);
        % lxy(i) = length(xy{i});

        PXX =  CC.PixelIdxList;
        makeCentroids = find((~makeIMG)==1); %you don't want to do this calculation through all the NaN, so index for non NAN
        PX = PXX(~makeIMG);
        mx = nan(1,length(PXX));
        my = nan(1,length(PXX));
            for j = 1:length(PX)
            px = PX{j};
        %     [y,x] = ind2sub(CC.ImageSize,px); %x and y come out reverse of S.Centroid
            y = rem(px-1,imgsize(1))+1; %these two lines replace ind2sub
            x = (px-y)/imgsize(2) + 1;  %these two lines replace ind2sub
    %         mx(j) = sum(x)./numel(x);
    %         my(j) = sum(y)./numel(y);

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
        lxy(i) = length(xy{i});
        %%%%%%%%%%%%%%%%%%%%%
    end
    
    traject = nan(max(lxy),2,t);

    for i = 1:t
        traject(1:lxy(i),1:2,i) = xy{i};
    end
end
function trt = calculateTrackingLogical(Stacked)

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the channel directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(mstackPath)

    if isempty(ImageDetails.Scene)
        ImageDetails.Scene = SceneList{1};
    end

    if isempty(ImageDetails.Channel)
        ImageDetails.Channel = nucleus_seg;
    end
    
ChannelDirectory = dir(strcat('*',ImageDetails.Channel,'_*'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ImageDetails.Frame)
       ImageDetails.Frame = frameToLoad;
       t = ImageDetails.Frame;
    else
       t = ImageDetails.Frame;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   choose the channel image
%options are overlay of background
%overlay of fluorescent channels
%normal image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         imgfile = dir(strcat('*',ImageDetails.Frame,'*.tif'));
        cd(mstackPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',cell_seg,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
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
%                  disp('if -> elseif')
            else
%                 disp('if -> else')
                %dont do anything
            end
        else %if no cfoName, then 
                 channelfileObject = matfile(filename);
                 channelimgstack = channelfileObject.flatstack;
                 cfoName = char(channelfileObject.Properties.Source);%update cfoName
%                  disp('else')
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
    channelimg;
end

displayImageFunct(If,channelimg);
end


function contrast_Callback(~,~)
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

function contrast_CallbackNew(~,~)
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

function displayImageFunctNew(If,channelimg)
global displaycomments timeFrames lprcntlt prcntlt tcontrast lcontrast MainAxes displaytracking ImageDetails frameToLoad prcntlz lprcntlz prcntlk lprcntlk prcntl lprcntl D ExpDate cmap cmaplz adjuster


%determine current time Frame
    t = strcmp(frameToLoad,ImageDetails.Frame);
    t = find(t==1);

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

        %nothing for overlay

    else  %under normal circumstances
            if ifCHANGEofCHANNELorSCENE==1
                dispimg = channelimg-lcontrast;
                dispimg = dispimg./(tcontrast-lcontrast);
                dispimg = dispimg.*255;
                ifCHANGEofCHANNELorSCENE=0;
            end

        dispimg = channelimg-lcontrast;
        dispimg = dispimg./(tcontrast-lcontrast);
        dispimg = dispimg.*255;
        dispimg(dispimg == 255) =254;
        colormap(cmap);
        If = bwperim(If) | bwperim(imdilate(If,strel('disk',1)));
        dispimg(If>0)=255;
    end


himg = imagesc(dispimg);
himgax = get(himg,'Parent');
himgax.CLim = [0 256];
ttl = get(himgax,'Title');
t = ImageDetails.Frame;
set(ttl,'String',[ExpDate ' ' ImageDetails.Scene ' frame ' num2str(t) ' out of ' num2str(timeFrames)]);
set(ttl,'FontSize',12);


    if ~(t==1)
        if displaytracking==1
            traject = trackingTrajectories(frameToLoad,ImageDetails);
            
            himgax.NextPlot = 'add';
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

% saveChannelFiveImages
end

function displayImageFunct(If,channelimg)
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
global runIterate SceneList ImageDetails TC A frameToLoad Tracked trackingPath ExportName timeFrames segmentPath nucleus_seg

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
%             Tracked = FrickTrackCellsYeah(A,frameToLoad,pvalue,[]);
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
global psettings runIterate SceneList ImageDetails TC A frameToLoad Tracked trackingPath ExportName timeFrames segmentPath nucleus_seg

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
%             Tracked = FrickTrackCellsYeah(A,frameToLoad,pvalue,[]);
            TC =1;
            setSceneAndTime;
            
            
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




function  uiSegmentTimeLapseImages

global mstackPath subaxestwo  pStruct subaxes   exportdir  channelinputs adjuster cmapper tcontrast lcontrast   OGExpDate   cmap  A AA timeFrames  ImageDetails  SceneList  imgsize ExpDate
adjuster=0;
tcontrast = 99;
lcontrast = 1;



clearvars -global SceneDirectoryPath

ImageDetails = InitializeImageDetails;

%%% set colormap for the images %%%
cmap = colormap(gray(255));
% cmap = colormap(magma(255));
% cmap = colormap(inferno(255));
% cmap = colormap(plasma(255));
cmap(255,:)=[1 0 0];
cmapper = cmap;
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose directory of experiment to track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set directory to location of code being used (generally external harddrive
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

    cd(parentdir)
    cd ..

%set parent directory from user input
    A = uigetdir;
    AA = 'D:\Users\zeiss\Documents\MATLAB';
    cd(A)
    mstackName = 'flat mstack';
    experimentdir = A;
    mstackPath = strcat(experimentdir,'/',mstackName);
    %subdirectories should include
    %> [ flatfield_corrected ]
        %> [ ####date## smad3g smFISH_scene_s## ]
            %> [ c#_flat ]     [ tiffs ]
                %need to load up the NucleusBinary_flat images
            
           

%determine date of experiment
    cd (mstackPath)
    [a,b] = regexp(A,'201[0-9]');
    [c,d] = regexp(A,'exp[0-9]');
    ExpDate = A(a:b+6);OGExpDate = A(a:d); [a,~] = regexp(ExpDate,'_');ExpDate(a) = '-';
    disp(A)

%load associated metadata
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

    channelstoinput = dosestructstruct.channelNameSwapArray;
    channelinputs =channelregexpmaker(channelstoinput);
    bkg = dosestructstruct.BACKGROUND;
    imgsize = dosestructstruct.dimensions;


%set up regexp parameters
    BACKGROUND = bkg{1};
    bkarray = bkarraymaker(BACKGROUND); %'(s01|s02|s03)'
    bkinputs =channelregexpmaker(bkarray); %'(chanstrA|chanstrB|chanstrC)'


%determine how many scenes are present
    dirlist = dir(mstackPath);
    [~,~,~,d] = regexp({dirlist.name},'s[0-9]+');
    dlog = ~cellfun(@isempty,d,'UniformOutput',1); 
    dcell = d(dlog);
    SceneList = unique(cellfun(@(x) x{1},dcell,'UniformOutput',0));

%remove background scenes from list
    [~,~,~,d] = regexp(SceneList,bkinputs);
    bkgscenelog = cellfun(@isempty,d,'UniformOutput',1);
    SceneList = SceneList(bkgscenelog);





%determine the number of time frames per scene
    cd(A)
    cd(mstackPath)
    filelist = dir('*.mat');
    fnames = {filelist.name};
    fileName = fnames{1};
    fileObject = matfile(fileName);
    dim = size(fileObject,'flatstack');
    timeFrames = dim(3);




%determine the number of channels
folderlist = dir(strcat('*','*'));
channelinputsUnderscore =channelregexpmaker(channelstoinput);

    [~,~,~,channelsListed] = regexp([folderlist.name],channelinputsUnderscore);
    channelList = unique(channelsListed);
            for i=1:length(channelList)
                chan = channelsListed{i};
                [a,~] = regexp(chan,'_');
                chan(a) = [];
                channelList{i} = chan;
            end
            
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Set up  user interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;
% f.Visible ='off';
f.Units = 'pixels';
f.Position =[10,10,1800,1000];


buttonwidth = 80;
buttonheight = 50;

ypositions = sort([100:20:1000],'descend');
xpositions = ones(1,length(ypositions)).*1600;

        mmm=1;
% htexttwo 
uicontrol('Style','text','String','To choose channel push 1, 2, or 3',...
          'Position',[xpositions(mmm),ypositions(mmm),buttonwidth,buttonheight]);
        mmm=mmm+1;mmm=mmm+1;mmm=mmm+1;

% hNextFrame
uicontrol('Style','pushbutton',...
    'String','NextFrame [f]',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@nextbutton_callback);
% hPreviousFrame
uicontrol('Style','pushbutton',...
    'String','Previous frame [a]',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@prevbutton_callback);
        mmm=mmm+1;mmm=mmm+1;
        
% hGoToFrame
uicontrol('Style','pushbutton',...
    'String','Go to Frame',...
    'Position',[xpositions(mmm)-120,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@gotobutton_callback);

        
% hFinalFrame
uicontrol('Style','pushbutton',...
    'String','FinalFrame [g]',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@finalbutton_callback);
% hFirstFrame
uicontrol('Style','pushbutton',...
    'String','First frame [z]',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@firstbutton_callback);
% hFirstFrame
uicontrol('Style','pushbutton',...
    'String','saveSomethingCallback',...
    'Position',[xpositions(mmm)-300,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@saveSomethingCallback);
        mmm=mmm+1; mmm=mmm+1;

% htextone
uicontrol('Style','text','String','Choose Scene',...
    'Position',[xpositions(mmm),ypositions(mmm)-buttonheight./2,buttonwidth,buttonheight]);

% htextone
uicontrol('Style','text','String','Choose Channel',...
    'Position',[xpositions(mmm)-200,ypositions(mmm)-buttonheight./2,buttonwidth,buttonheight]);
        mmm=mmm+1;
% hpopup
uicontrol('Style','popupmenu',...
    'String',SceneList',...
    'Position',[xpositions(mmm)-buttonwidth./2,ypositions(mmm),buttonwidth.*1.5,buttonheight./2],...
    'Callback',@popup_menu_Callback);
        
        
% hpopup
uicontrol('Style','popupmenu',...
    'String',channelList',...
    'Position',[xpositions(mmm)-200-buttonwidth./2,ypositions(mmm),buttonwidth.*1.5,buttonheight./2],...
    'Callback',@popup_menu_Callback_channels);
        mmm=mmm+1;mmm=mmm+1; mmm=mmm+1;
        
  
        %%%%
        %%%%
% hSaveTrackingAs
uicontrol('Style','pushbutton',...
    'String','SaveTrackingAs',...
    'Position',[xpositions(mmm)-buttonwidth./2,ypositions(mmm)-buttonheight,buttonwidth.*2,buttonheight.*2],...
    'Callback',@saveTrackingFileAs_callback);
        mmm=mmm+1; mmm=mmm+1; mmm=mmm+1; mmm=mmm+1; mmm=mmm+1;

% hLoadTracking
uicontrol('Style','pushbutton',...
    'String','LoadTracking',...
    'Position',[xpositions(mmm),ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@loadTrackingFile_callback);
        mmm=mmm+1; mmm=mmm+1; mmm=mmm+1;
       

        
       
f.Visible = 'on'   ;
f.Units = 'normalized';
for i = 1:length(f.Children)
   hhh = f.Children(i);
   hhh.Units = 'normalized';
end
% 
% MainAxes = axes;
% MainAxes.Units = 'pixels';
% MainAxes.XTick=[];
% MainAxes.YTick = [];
% imgdim = 512.*1.8;
% Position = [25 25 imgdim imgdim];
% % Position = [0.1 0.3 0.65 0.65];
% MainAxes.Position = Position;


channelimglength = 9;
xinit = 0.02;
yinit = 0.025;
w = 0.15;
h = 0.15;
xspacefactor = 0.3;
yspacefactor = 1;
x = [xinit xinit+w+(xspacefactor*w) xinit+(w+(xspacefactor*w)).*2 xinit xinit+w+(xspacefactor*w) xinit+(w+(xspacefactor*w)).*2 xinit xinit+w+(xspacefactor*w) xinit+(w+(xspacefactor*w)).*2];
y = fliplr([yinit yinit yinit yinit+h+(yspacefactor*w) yinit+h+(yspacefactor*w) yinit+h+(yspacefactor*w)   yinit+(w+(yspacefactor*w)).*2 yinit+(w+(yspacefactor*w)).*2 yinit+(w+(yspacefactor*w)).*2]);

for i=1:channelimglength
    ax= axes();
    ax.Position = [x(i) y(i) w h];
    ax.Units = 'inches';
    pos = ax.Position;
    pos(4) = pos(3);
    ax.Position = pos;
    ax.Units = 'normalized';
    ax.XTick = [];
    ax.YTick = [];
    subaxes(i) = ax;
end

w=0.15;
h=0.15;
x = [0.62 0.8 0.62 0.8];
y = [0.32 0.32 0.02 0.02];
for i=1:4
    ax= axes();
    ax.Position = [x(i) y(i) w h];
    ax.Units = 'inches';
    pos = ax.Position;
    pos(4) = pos(3);
    ax.Position = pos;
    ax.Units = 'normalized';
    ax.XTick = [];
    ax.YTick = [];
    subaxestwo(i) = ax;
end



%define parameter structure and default parameter values
pStruct = struct();
parameterDefaults.EGFP = [106 1 15 1];
parameterDefaults.CFP = [20 1 2 1];
parameterDefaults.mKate = [40 1 2 1];
parameterDefaults.Hoechst = [20 1 2 1];
parameterDefaults.DIC = [20 1 2 1];
parameterStrings = {'nucDiameter','threshFactor','sigmaScaledToParticle','noparametercurrently'};
for p = 1:length(parameterStrings)
    pString = char(parameterStrings{p});
    for c = 1:length(channelList)
        cstr = char(channelList{c});
        cString = alterChanName(cstr);
        pd = parameterDefaults.(cString);
        pStruct.(cString).(pString) = pd(p); 
    end
end

pStruct = loadSegmentParameters(pStruct,FileName,exportdir); %loads saved value of pStruct

set(f,'KeyPressFcn',@keypress);
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
    disp('no saved parameters')
end

    

end

function saveSomethingCallback(~,~)
global pStruct exportdir OGExpDate

cd(exportdir)
filename = strcat('*',OGExpDate,'*metaData*')
filelist = dir(filename)
savenamebase = char((filelist.name));
savename = strcat(OGExpDate,'-segmentParameters.mat');
save(savename,'pStruct');


end




function updateSliders
global pStruct ImageDetails sliderOne sliderOneTxt

%slider1
sliderx = 0.72;
% slidery = 0.6;
sliderw = 0.1;
sliderh = 0.02;
slidertextw = 0.1;
sliderspace = 0.1;





channel = ImageDetails.Channel;

%nucDiameter
str = 'nucDiameter';
val.(str) = pStruct.(channel).(str);
slidery = 0.75;
minz = 1;
maxz = 400;
sliderOne.(str) = uicontrol('Style', 'slider','String',str,'Min',minz,'Max',maxz,'SliderStep',[1 1]./((maxz-minz)./5),'Value',val.(str),'Position', [1 1 1 1],...
        'Callback', @sliderOneAdjust); 
    sliderOne.(str).Units='normalized';
    sliderOne.(str).Position = [sliderx slidery sliderw sliderh];
sliderOneTxt.(str) = uicontrol('Style','text','Units','Normalized','Position',[1 1 1 1],'String',strcat(str,'=',num2str(uint8(val.(str)))));
    sliderOneTxt.(str).Units= 'Normalized';
    sliderOneTxt.(str).Position = [sliderx-sliderspace slidery slidertextw sliderh];  

    
    
    
% threshFactor
str = 'threshFactor';
val.(str) = pStruct.(channel).(str);
slidery = 0.725;
minz = 0.4;
maxz = 3;
sliderOne.(str) = uicontrol('Style', 'slider','String',str,'Min',minz,'Max',maxz,'SliderStep',[1 1]./(20*(maxz-minz)),'Value',val.(str),'Position', [1 1 1 1],...
        'Callback', @sliderOneAdjust); 
    sliderOne.(str).Units='normalized';
    sliderOne.(str).Position = [sliderx slidery sliderw sliderh];
sliderOneTxt.(str) = uicontrol('Style','text','Units','Normalized','Position',[1 1 1 1],'String',strcat(str,'=',num2str((val.(str)))));
    sliderOneTxt.(str).Units= 'Normalized';
    sliderOneTxt.(str).Position = [sliderx-sliderspace slidery slidertextw sliderh];  

    
   
    
    
% sigmaScaledToParticle (the divide diameter by this factor to get sigma for
% gaussian smoothing)
str = 'sigmaScaledToParticle';
val.(str) = pStruct.(channel).(str);
slidery = 0.7;
minz = 1;
maxz = 40;
sliderOne.(str) = uicontrol('Style', 'slider','String',str,'Min',minz,'Max',maxz,'SliderStep',[1 1]./(maxz-minz),'Value',val.(str),'Position', [1 1 1 1],...
        'Callback', @sliderOneAdjust); 
    sliderOne.(str).Units='normalized';
    sliderOne.(str).Position = [sliderx slidery sliderw sliderh];
sliderOneTxt.(str) = uicontrol('Style','text','Units','Normalized','Position',[1 1 1 1],'String',strcat(str,'=',num2str(uint8(val.(str)))));
    sliderOneTxt.(str).Units= 'Normalized';
    sliderOneTxt.(str).Position = [sliderx-sliderspace slidery slidertextw sliderh];   
  
    
    
% noparemtercurrently
str = 'noparametercurrently';
val.(str) = pStruct.(channel).(str);
slidery = 0.675;
minz = 1;
maxz = 400;
sliderOne.(str) = uicontrol('Style', 'slider','String',str,'Min',minz,'Max',maxz,'SliderStep',[1 1]./(maxz-minz),'Value',val.(str),'Position', [1 1 1 1],...
        'Callback', @sliderOneAdjust); 
    sliderOne.(str).Units='normalized';
    sliderOne.(str).Position = [sliderx slidery sliderw sliderh];
sliderOneTxt.(str) = uicontrol('Style','text','Units','Normalized','Position',[1 1 1 1],'String',strcat(str,'=',num2str(uint8(val.(str)))));
    sliderOneTxt.(str).Units= 'Normalized';
    sliderOneTxt.(str).Position = [sliderx-sliderspace slidery slidertextw sliderh];   

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

function channelinputs =channelregexpmakerUnderscore(channelstoinput)
    channelinputs = '(';
    for i=1:length(channelstoinput) % creates a string of from '(c1|c2|c3|c4)' for regexp functions
        if i ==1
        channelinputs = strcat(channelinputs,channelstoinput{i},'_');
        elseif i < length(channelstoinput)
            channelinputs = strcat(channelinputs,'|',channelstoinput{i},'_');
        else
            channelinputs = strcat(channelinputs,'|',channelstoinput{i},'_)');
        end
    end
end


function sliderOneAdjust(source,~)
    global pStruct  ImageDetails sliderOneTxt
    channel = ImageDetails.Channel;
    str = source.String;
%     threshinput.(str) =source.Value;
%     zerostrel = round(source.Value);

    if strcmpi(str,'threshFactor')
        valupdate = source.Value;  
    else
        valupdate = round(source.Value);
    end
    source.Value = valupdate;
    pStruct.(channel).(str) = valupdate; 
    disp(valupdate)
   source.Visible = 'off';
   sliderOneTxt.(str).String = 'waiting...';
   pause(0.001);
    setSceneAndTime
    disp('done')
    source.Visible = 'on';
    
    sliderOneTxt.(str).String = strcat(str,'=',num2str(pStruct.(channel).(str)));
    
end


function plotTestOut(testOut,channel)
    global subaxestwo

    if strcmp(channel,'mKate')
    stringsToTest = {'rawMinusLPScaled','Ihcf','gradmag2','Ieg'};
    elseif strcmp(channel,'Hoechst')
    stringsToTest = {'rawMinusLPScaled','Ihcf','gradmag2','Ieg'};
    else
%     stringsToTest = {'rawMinusLPScaled','Ih','Ihcd','Shapes'};
    stringsToTest = {'initialIh','imgRawDenoised','Im','imgWWW'};
    end
    for i = 1:length(subaxestwo)
    axes(subaxestwo(i))
    str = stringsToTest{i};
    img = testOut.(str);
    imagesc(img);title(str);
    end
    
%         testOut.img = img;
%         testOut.imgRawDenoised = imgRawDenoised;
%         testOut.imgLowPass = imgLowPass;
%         testOut.rawMinusLP = rawMinusLP;
%         testOut.rawMinusLPScaled = rawMinusLPScaled;
%         testOut.Ih = Ih;
%         testOut.L = zeros[512 512];
%    

end
  


function [IfFinal,testOut] = segmentationNucleus(FinalImage,segmentPath,nucleus_seg,nucleusFileName,pStruct)
    testOut = struct();
    frames = 1;
    img = FinalImage(:,:,frames); 
    if strcmp(nucleus_seg,'Hoechst')
        [~,testOut] = segmentHoechstNuclei(img,nucleus_seg,pStruct,frames);  
    else
        [~,testOut] = segmentNuclei(img,nucleus_seg,pStruct,frames);
    end
    
    IfFinal = false(size(FinalImage));
    for frames = 1:size(FinalImage,3)
        img = FinalImage(:,:,frames); 
        if strcmp(nucleus_seg,'Hoechst')
            [If,~] = segmentHoechstNuclei(img,nucleus_seg,pStruct,frames);  
        else
            [If,~] = segmentNuclei(img,nucleus_seg,pStruct,frames);
        end
        IfFinal(:,:,frames)=If;
    end
                
    %save here if running actual segmentation
end
function [IfFinal,testOut] = segmentationImageBackground(FinalImage,segmentPath,background_seg,backgroundFileName,pStruct)
testOut = struct();
    frames = 1;
    img = FinalImage(:,:,frames); 
    [~,testOut] = segmentCellBackground(img,background_seg,pStruct,frames);
    
    IfFinal = false(size(FinalImage));
    for frames = 1:size(FinalImage,3)
        img = FinalImage(:,:,frames); 
        [If,~] = segmentCellBackground(img,background_seg,pStruct,frames);
        IfFinal(:,:,frames)=If;
    end
                
    %save here if running actual segmentation


% parameters

end
function [IfFinal,testOut] = segmentationDIC(FinalImage,subdirname,scenename,filename,channel)
global  pStruct foldernameglobal
cd(subdirname)


foldername = foldernameglobal;


% parameters
nucDiameter = pStruct.(channel).nucDiameter;
threshFactor = pStruct.(channel).threshFactor;
sigmaScaledToParticle = pStruct.(channel).sigmaScaledToParticle;
kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize

finalerode=2;



% prepareCcodeForAnisotropicDiffusionDenoising(denoisepath)

%start
for frames = 1:size(FinalImage,3)
%Smooth Image using Anisotropic Diffusion
% Options.Scheme :  The numerical diffusion scheme used
%                     'R', Rotation Invariant, Standard Discretization 
%                          (implicit) 5x5 kernel (Default)
%                     'O', Optimized Derivative Kernels
%                     'I', Implicit Discretization (only works in 2D)
%                     'S', Standard Discretization
%                     'N', Non-negativity Discretization
%   Options.T  :      The total diffusion time (default 5)
%   Options.dt :      Diffusion time stepsize, in case of scheme H,R or I
%                     defaults to 1, in case of scheme S or N defaults to
%                     0.15. 
%   Options.sigma :   Sigma of gaussian smoothing before calculation of the
%                     image Hessian, default 1.                   
%   Options.rho :     Rho gives the sigma of the Gaussian smoothing of the 
%                     Hessian, default 1.
%   Options.verbose : Show information about the filtering, values :
%                     'none', 'iter' (default) , 'full'
%   Options.eigenmode : There are many different equations to make an diffusion tensor,
%						this value (only 3D) selects one.
%					    0 (default) : Weickerts equation, line like kernel
%						1 : Weickerts equation, plane like kernel
%						2 : Edge enhancing diffusion (EED)
%						3 : Coherence-enhancing diffusion (CED)
%						4 : Hybrid Diffusion With Continuous Switch (HDCS)
    img = FinalImage(:,:,frames); 
    imgRaw = gaussianBlurz(single(img),ceil(sigma./10),ceil(kernelgsize./10));

    imgW = wiener2(img,[1 20]);
    imgWW = wiener2(imgW,[20 1]);
    imgWWW = wiener2(imgWW,[5 5]);
    imgRawDenoised = imgWWW;
    denoiseVec = single(reshape(imgRawDenoised,size(imgRawDenoised,1)^2,1));
    highpoints = prctile(denoiseVec,95);
    imgRawDenoised(imgRawDenoised>highpoints) = highpoints;
    
 
%         Options.T = 5;
%         Options.dt = 1;
%         Options.Scheme = 'R';
%         Options.rho = 20;
%         Options.sigma = 20;
%         Options.verbose = 'none';ii
% %     imgRawDenoised = CoherenceFilter(imgRaw, Options);
% % imgRawDenoised = imgRaw;

    
    %Based on algorithm of Fast and accurate automated cell boundary determination for fluorescence microscopy by Arce et al (2013)   
    %LOW PASS FILTER THE IMAGE (scale the gaussian filter to diameter of
    %nuclei -- diameter of nuclei is about 50 to 60))
    
    imgLowPass = gaussianBlurz(single(imgRawDenoised),sigma,kernelgsize);
    rawMinusLP = single(imgRawDenoised) -single(imgLowPass);%%%%%%% key step!
    rawMinusLPvec = reshape(rawMinusLP,size(rawMinusLP,1)^2,1);
    globalMinimaValues = prctile(rawMinusLPvec,0.01);
    globalMinimaIndices = find(rawMinusLP < globalMinimaValues);
    LPscalingFactor = imgRawDenoised(globalMinimaIndices)./imgLowPass(globalMinimaIndices);
    imgLPScaled = imgLowPass.*nanmedian(LPscalingFactor);
    rawMinusLPScaled = single(imgRawDenoised) - single(imgLPScaled);


    %determine the threshold by looking for minima in log-scaled histogram
    %of pixels from rawMinusLPScaled
    rawMinusLPScaledContrasted = imadjust(uint16(rawMinusLPScaled));
    vecOG = single(reshape(rawMinusLPScaledContrasted,size(rawMinusLPScaledContrasted,1)^2,1));
    logvecpre = vecOG; logvecpre(logvecpre==0)=[];
    logvec = log10(logvecpre);
    vec = logvec;
    [numbers,bincenters] = hist(vec,prctile(vec,1):(prctile(vec,99)-prctile(vec,1))/1000:max(vec));
    numbersone = medfilt1(numbers, 10); %smooths curve
    numberstwo = medfilt1(numbersone, 100); %smooths curve
    fraction = numberstwo./sum(numberstwo);
    mf = max(fraction);
        %%%%%%%%%%%%%%%%%%%% Important parameters for finding minima of
        %%%%%%%%%%%%%%%%%%%% histogram
        left=0.5*mf;
        slopedown=0.4*mf;
        %%%%%%%%%%%%%%%%%%%%%
    leftedge = find(fraction > left,1,'first');
    insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
    threshLocation = bincenters(leftedge+insideslopedown-1);
    subtractionThreshold = threshLocation;

    if size(subtractionThreshold,1)==size(subtractionThreshold,2)
        else
         subtractionThreshold = mean(threshLocation);
    end


    subtractionThresholdScaled = (10.^subtractionThreshold).*threshFactor;
    subtracted = single(rawMinusLPScaledContrasted)-subtractionThresholdScaled;
    subzero = (subtracted<0);
    Ih = ~subzero;
    Im = Ih;
    If =Im;

    
    
    
    
I = -1.*rawMinusLPScaled;
waterBoundary = Ih;

%gradmag
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(single(I), hy, 'replicate');
Ix = imfilter(single(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

%Smoothing
% I = Ih;
width = round(nucDiameter./4);
se = strel('disk', width);
Io = imopen(I, se);
Ie = imerode(Io, se);
Ieg = gaussianBlurz(Ie,sigma./2,kernelgsize);
%     width = round(nucDiameter./10);
%     Ime = imerode(Ihcf,strel('disk',width));
%     Imeo = imopen(Ime,strel('disk',width));
%     Ieg(~Imeo)=0;
fgm = imregionalmax(Ieg);
width = round(nucDiameter./20);
fgm4 = imdilate(fgm,strel('disk',width));
% fgm4 =fgm;
bw = Im;
D = bwdist(bw);
DL = watershed(D,4);
bgm = DL == 0;
gradmag2 = imimposemin(gradmag, bgm | fgm4);

L = watershed(gradmag2,8);
L(waterBoundary<1) = 0;
% If = L>1;



    
    time = tsn{frames};
    tim = time(2:end);
    IfFinal(:,:,frames)=If;
    
       if frames==1
        testOut.img = img;
        testOut.I = -1.*rawMinusLPScaled;
        testOut.imgRawDenoised = imgRawDenoised;
        testOut.imgLowPass = imgLowPass;
        testOut.rawMinusLP = rawMinusLP;
        testOut.rawMinusLPScaled = rawMinusLPScaled;
        testOut.Ih = Ih;
%         testOut.Ihc = Ihc;
        testOut.Im = Im;
%         testOut.Ihcd = Ihcd;
        testOut.L = zeros([512 512]);
%         testOut.gradmag = gradmag;
        testOut.gradmag = zeros(size(img));
%         testOut.gradmag2 =  gradmag2;
        testOut.gradmag2 = zeros(size(img));
%         testOut.Ie = Ie;
        testOut.Ie = zeros(size(img));
%         testOut.fgm4 = fgm4;
        testOut.fgm4 = zeros(size(img));
%         testOut.Ieg = Ieg;
        testOut.Ieg = zeros(size(img));
        testOut.Shapes = zeros(size(img));
%         testOut.waterBoundary = waterBoundary;

       end
    
end


stophere=1;
end




function bw = gaussianBlurz(im,sigma,kernelgsize,varargin)

filtersize = [kernelgsize kernelgsize];
kernelg = fspecial('gaussian',filtersize,sigma);


gFrame = imfilter(im,kernelg,'repl');

if ~isempty(varargin)
    bw=gFrame.*uint16(varargin{1}>0);
else
    bw=gFrame;
end
end

function keypress(fig_obj,~)
global  ImageDetails displaycomments
key = get(fig_obj,'CurrentKey');

switch key
    case '1'
        ImageDetails.Channel = 'EGFP';
        setSceneAndTime
    case '2'
        ImageDetails.Channel = '_Hoechst';
        setSceneAndTime    
    case '3'
        ImageDetails.Channel = 'mKate';
        setSceneAndTime
    case '4'
        ImageDetails.Channel = 'DIC';
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
global framesForDir ImageDetails 

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = framesForDir{1};
end

Idx = strcmp(ImageDetails.Frame,framesForDir);
idx = find(Idx == 1);
if idx == length(framesForDir)
else
idx = idx + 1;
end
ImageDetails.Frame = framesForDir{idx};


setSceneAndTime
end
function prevbutton_callback(~,~) 
global framesForDir ImageDetails 

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = framesForDir{1};
end

Idx = strcmp(ImageDetails.Frame,framesForDir);
idx = find(Idx == 1);
if idx == 1
else
idx = idx - 1;
end
ImageDetails.Frame = framesForDir{idx};


setSceneAndTime
end
function finalbutton_callback(~,~)
global framesForDir ImageDetails 

idx = length(framesForDir);
ImageDetails.Frame = framesForDir{idx};

setSceneAndTime
end
function firstbutton_callback(~,~)
global framesForDir ImageDetails 
idx = 1;
ImageDetails.Frame = framesForDir{idx};

setSceneAndTime
end
function gotobutton_callback(~,~)
global framesForDir ImageDetails 

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = framesForDir{1};
end

prompt = {'Go to which frame'};
dlg_title = 'Go to frame...';
idx = str2num(cell2mat(inputdlg(prompt,dlg_title)));

ImageDetails.Frame = framesForDir{idx};


setSceneAndTime
end


%choose scenes
function nextscenebutton_Callback(~,~) 
global   ImageDetails SceneList A SceneDirectoryPath




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



cd(A)
cd('flatfield_corrected')
SceneDirectory = dir (strcat('*',ImageDetails.Scene,'*'));
cd(SceneDirectory.name)
SceneDirectoryPath = pwd;
setSceneAndTime
end
function prevscenebutton_Callback(~,~) 
global   ImageDetails SceneList A SceneDirectoryPath

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

cd(A)
cd('flatfield_corrected')
SceneDirectory = dir (strcat('*',ImageDetails.Scene,'*'));
cd(SceneDirectory.name)
SceneDirectoryPath = pwd;
loadTrackingFile_callback([],[])
setSceneAndTime
end


function popup_menu_Callback_channels(source,~)
global ImageDetails Tracked 

Tracked=[];

% Determine the selected data set.
 str = source.String;
 val = source.Value;
 channel = char(str{val});

ImageDetails.Channel = channel;
setSceneAndTime
end

function popup_menu_Callback(source,~) 
global ImageDetails Tracked 

Tracked=[];

% Determine the selected data set.
 str = source.String;
 val = source.Value;
 pvalue = char(str{val});

ImageDetails.Scene = pvalue;
setSceneAndTime

end


%% Image Display functions
function setSceneAndTime
global pStruct timeFrames mstackPath TC A  framesForDir ImageDetails  Tracked SceneList  SceneDirectoryPath imgfile imgsize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the channel directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(mstackPath)

    if isempty(ImageDetails.Scene)
        ImageDetails.Scene = SceneList{1};
    end

SceneDirectoryPath = mstackPath;

    if isempty(ImageDetails.Channel)
        ImageDetails.Channel = 'EGFP';
    end
    
ChannelDirectory = dir(strcat('*',ImageDetails.Channel,'*'));
    if isempty(ChannelDirectory) && ~strcmp(ImageDetails.Channel,'overlay')
        ImageDetails.Channel = 'mKate';
        ChannelDirectory = dir(strcat('*',ImageDetails.Channel,'_*'));
    elseif isempty(ChannelDirectory)
        ChannelDirectory = dir(strcat('*','mKate','_*'));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ImageDetails.Frame)
       ImageDetails.Frame = 1;
    end
t=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose the channel images
%options are overlay of background
%overlay of fluorescent channels
%normal image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgfile = dir(strcat('*',ImageDetails.Frame,'*.tif'));
    if strcmp(ImageDetails.Channel,'BKGbinary')                    %overlay background
        bkgimg = false(imgsize(1),imgsize(2),1);
        bkgimg = imread(char(imgfile.name));
        bkgimg(bkgimg>0) = 1;
        channelimg = ~logical(bkgimg);


%         ChannelDirectory = dir(strcat('*','EGFP_','*'));
%         cd(ChannelDirectory.name)
        cd ..        
        cd('tiffs')
%         ff = dir(strcat('EGFP','*'));
        ff = dir(strcat('*','EGFP','*'));
            if isempty(ff)
                ff = dir(strcat('*','mKate','*'));
            end
        channelimg = single(loadUpTiffStackFrame(char(ff.name),t));
% 
%         cimgfile = dir(strcat('*',ImageDetails.Frame,'*.tif'));
%         cimg = imread(char(cimgfile.name));
%         channelimg = double(cimg);
        prim = imdilate(bwperim(~logical(bkgimg)),strel('square',2));
        channelimg(prim) = max(max(channelimg));
    elseif strcmp(ImageDetails.Channel,'overlay')
           
        
         cd ..        
        cd('tiffs')
        ff = dir(strcat('EGFP','*'));
            if isempty(ff)
               channelimg = zeros(512,512);
            else
            channelimg = single(loadUpTiffStackFrame(char(ff.name),t));
            end
        imgone = channelimg;
        cd ..
        
        cd('tiffs')
        ff = dir(strcat('*','Hoechst','*'));
            if isempty(ff)
               channelimg = zeros(512,512);
            else
            channelimg = single(loadUpTiffStackFrame(char(ff.name),t));
            end
        imgtwo = channelimg;
        cd .. 
        
        
        cd('tiffs')
        ff = dir(strcat('*','mKate_','*'));
        channelimg = single(loadUpTiffStackFrame(char(ff.name),t));
        imgthree = channelimg;
    
        ff = dir(strcat('*','DIC','*'));
        channelimg = single(loadUpTiffStackFrame(char(ff.name),t));
        imgfour = channelimg;
        
        
        channelimg = zeros(size(imgone,1),size(imgone,2),3);
%         imgone(imgone<imgfour) = imgfour(imgone<imgfour);
%         imgtwo(imgtwo<imgfour) = imgfour(imgtwo<imgfour);
%         imgthree(imgthree<imgfour) = imgfour(imgthree<imgfour);
        channelimg(:,:,1) = imgone;
        channelimg(:,:,2) = imgtwo;
        channelimg(:,:,3) = imgthree;
        channelimg(:,:,4) = imgfour;
        
%         
%         channelimg = uint8(channelimg);
%         %make uint8?
        
        

    else
        cd(mstackPath)
        ff = dir(strcat('*',ImageDetails.Scene,'*',ImageDetails.Channel,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
        channelspacing = round(linspace(16,timeFrames,9));
%         channelspacing = [1 4 9 12 15 18 20 22 26];
        if length(channelspacing)>timeFrames
            channelspacing = 1:timeFrames;
        end

        
        fname = char(ff.name);
        fileObject = matfile(fname);
        imgstack = single(fileObject.flatstack);
        channelimgstack = imgstack(:,:,channelspacing);
        img = imgstack(:,:,t);
        channelimg = img;
        
%         channelimg = double(loadUpTiffStackFrame(char(ff.name),t));
    %     channelimg = double(imread(char(imgfile.name)));    %load normal image
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


If = zeros(size(channelimg),'single');

%     subdirname = char(subdir);
%         [sceneinfo,b] = regexp(subdirname,'s[0-9]+');
%         scenename = subdirname(sceneinfo:b);
%         cd(subdirname)

FinalImage = channelimgstack;
subdirname = SceneDirectoryPath;
scenename = ImageDetails.Scene;
filename = imgfile;
filename = char(filename.name);
channel = ImageDetails.Channel;

nucleus_seg = ImageDetails.Channel;
background_seg = ImageDetails.Channel;
nucleusFileName = filename;
backgroundFileName = filename;
segmentPath = mstackPath;

if strcmp(ImageDetails.Channel,'EGFP')
[IfStack,testOut] = segmentationImageBackground(FinalImage,segmentPath,background_seg,backgroundFileName,pStruct);  
elseif strcmp(ImageDetails.Channel,'mKate')
[IfStack,testOut] = segmentationNucleus(FinalImage,segmentPath,nucleus_seg,nucleusFileName,pStruct);
elseif strcmp(ImageDetails.Channel,'Hoechst')
[IfStack,testOut] = segmentationNucleus(FinalImage,segmentPath,nucleus_seg,nucleusFileName,pStruct);
elseif strcmp(ImageDetails.Channel,'CFP')
[IfStack,testOut] = segmentationNucleus(FinalImage,segmentPath,nucleus_seg,nucleusFileName,pStruct);
elseif strcmp(ImageDetails.Channel,'DIC')
[IfStack,testOut] = segmentationDIC(FinalImage,subdirname,scenename,filename,channel);
end


displayImageFunct(IfStack,channelimgstack,channelspacing);
updateSliders

plotTestOut(testOut,channel)




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

function displayImageFunct(IfStack,channelimgstack,channelspacing)
global timeFrames subaxes displaycomments lprcntlt prcntlt tcontrast lcontrast MainAxes displaytracking ImageDetails framesForDir prcntlz lprcntlz prcntlk lprcntlk prcntl lprcntl D ExpDate cmap cmaplz adjuster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   display the images overlayed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axes(MainAxes);
children = findobj(MainAxes,'Type','image');
delete(children);



for i=1:size(channelimgstack,3)
    If = IfStack(:,:,i);
    channelimg = channelimgstack(:,:,i);
    ax = subaxes(i);
    axes(ax);
    
%scale the image brightness
if strcmp(ImageDetails.Channel,'overlay') %when overlay display is desired
    imgone = channelimg(:,:,1);
    imgtwo = channelimg(:,:,2);
    imgthree = channelimg(:,:,3);
    imgfour = channelimg(:,:,4);
    channelimg = zeros(size(channelimg,1),size(channelimg,2),3,'single');
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
        cimgline = reshape(channelimg,[1 size(channelimg,1).*size(channelimg,2)]);
        lprcntl = prctile(cimgline,lcontrast);
%         prcntl = prctile(cimgline,tcontrast)-lprcntl;
        prcntl = prctile(cimgline-lprcntl,tcontrast);


    channelimg = uint8(((channelimg-lprcntl)./prcntl).*255);
    channelimg(channelimg == 255) =254;
    colormap(cmap);
    If = bwperim(If);
    If = imdilate(If,strel('disk',1));
    channelimg(If>0)=255;
end



t = channelspacing(i);
himg = imagesc(channelimg);
himgax = get(himg,'Parent');
himgax.CLim = [0 256];
ttl = get(himgax,'Title');
set(ttl,'String',strcat(ExpDate,'...',ImageDetails.Scene,'...frame ',num2str(t),' out of', num2str(timeFrames)));
set(ttl,'FontSize',12);
himgax.YTick = [];
himgax.XTick = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saveChannelFiveImages
end

%% saving functions
function saveTrackingFileAs_callback(~,~)
global  SceneDirectoryPath Tracked ExportName
cd(SceneDirectoryPath)
prompt = 'filename of tracking structure to be saved?';
dlg_title = 'save tracking structure as...specific filename';
filename = char(inputdlg(prompt,dlg_title));
% save(strcat(filename,ExportName,'.mat'),'Tracked')
end

%% functions for determining variables

function ImageDetails = InitializeImageDetails

ImageDetails.Scene=[];
ImageDetails.Channel=[];
ImageDetails.Frame=[];

end


%% load up images
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

function chanstruct = alterChanName(chan)
[a,~] = regexp(chan,'(_|\W|\s|[0-9])'); %remove underscore, dashes, or whitespace
chanstruct = chan;
chanstruct(a) = [];

% if isempty(chanstruct)
%     chanstruct = 'CFP';
% end
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
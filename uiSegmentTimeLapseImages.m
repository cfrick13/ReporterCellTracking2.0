function  uiSegmentTimeLapseImages

global mstackPath subaxestwo channelList  segInstructList segmentList pStruct subaxes   exportdir  seginputs channelinputs adjuster cmapper tcontrast lcontrast   OGExpDate   cmap  A AA timeFrames  ImageDetails  SceneList  imgsize ExpDate
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
    
    segInstruct = dosestructstruct.segInstruct;
    fnames = fieldnames(segInstruct);
    segmentList = cell(1,length(fnames));
    segmentListDisp = cell(1,length(fnames));
    segInstructList = cell(1,length(fnames));
    for i = 1:length(segmentList)
        str = fnames{i};
        segmentListDisp{i} = [str '-' segInstruct.(str)];
        segmentList{i} = segInstruct.(str);
        segInstructList{i} = str;
    end
    seginputs = channelregexpmaker(fnames);
    
%set up regexp parameters
    BACKGROUND = bkg{1};
    bkarray = bkarraymaker(BACKGROUND); %'(s01|s02|s03)'
    bkinputs =channelregexpmaker(bkarray); %'(chanstrA|chanstrB|chanstrC)'


%determine how many scenes are present
    dirlist = dir(mstackPath);
    [~,~,~,d] = regexp({dirlist.name},'s[0-9]++');
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
    if max(size(dim))>2
        timeFrames = dim(3);
    else
        timeFrames = 1;
    end




%determine the number of channels
folderlist = dir(strcat('*','*'));
channelinputs =channelregexpmaker(channelstoinput);

    [~,~,~,channelsListed] = regexp([folderlist.name],channelinputs);
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
f.Position =[70 90 1800 900];


bW = 80; %buttonWidth
sW = 10; %spacerWidth
bH = 50; %buttonHeigth

%row1
    % hNextFrame
    xPi = 1500;
    yPi = 850;
    xP = xPi;
    yP = yPi;
    uicontrol('Style','pushbutton','String','NextFrame [f]',...
        'Position',[xP,yP,bW,bH],...
        'Callback',@nextbutton_callback);
    % hPreviousFrame
    xP = xP+sW+bW;
    uicontrol('Style','pushbutton','String','Previous frame [a]',...
        'Position',[xP,yP,bW,bH],...
        'Callback',@prevbutton_callback);

    % hGoToFrame
    xP = xP+sW+bW;
    uicontrol('Style','pushbutton','String','Go to Frame',...
        'Position',[xP,yP,bW,bH],...
        'Callback',@gotobutton_callback);
%row2
    % hFirstFrame
    xP = xPi;
    yP = yP-bH-sW;
    uicontrol('Style','pushbutton','String','First frame [z]',...
        'Position',[xP,yP,bW,bH],...
        'Callback',@firstbutton_callback);

    % hFinalFrame
    xP = xP+sW+bW;
    uicontrol('Style','pushbutton','String','FinalFrame [g]',...
        'Position',[xP,yP,bW,bH],...
        'Callback',@finalbutton_callback);

%section3
bW=bW*1.5;
    % htextone
    xP = xPi-(sW*3)-bW;
    yP = yPi;
    uicontrol('Style','text','String','Choose Scene',...
        'Position',[xP,yP,bW,bH./2]);

    % htextone
    xP = xP-(sW)-bW;
    uicontrol('Style','text','String','Channel View',...
        'Position',[xP,yP,bW,bH./2]);
    
    % choose segment channel
    xP = xP-(sW)-bW-bW./2;
    uicontrol('Style','text','String','Segmentation View',...
        'Position',[xP,yP,bW+bW./2,bH./2]);

%section4 (dropdown "popupmenus")
    % hpopupSceneList
    xP = xPi-(sW*3)-bW;
    yP = yP-(sW*2);
    uicontrol('Style','popupmenu',...
        'String',SceneList',...
        'Position',[xP,yP,bW,bH./2],...
        'Callback',@popup_menu_Callback);

    % hpopupChannelList
    xP = xP-(sW)-bW;
    uicontrol('Style','popupmenu',...
        'String',channelList',...
        'Position',[xP,yP,bW,bH./2],...
        'Callback',@popup_menu_Callback_channels);
    
    % hpopupChannelList
    xP = xP-(sW)-bW-bW./2;
    uicontrol('Style','popupmenu',...
        'String',segmentListDisp',...
        'Position',[xP,yP,bW+bW./2,bH./2],...
        'Callback',@popup_menu_Callback_segment);
    
%final Save button
    % hFirstFrame
    xP = xPi;
    yP = yPi-bH-bH-bH-bH-(sW*3);
    uicontrol('Style','pushbutton',...
        'String','saveSomethingCallback',...
        'Position',[xP,yP,bW.*2,bH.*2],...
        'Callback',@saveSomethingCallback);

        
  
        
       
f.Visible = 'on'   ;
% f.Units = 'normalized';
for i = 1:length(f.Children)
   hhh = f.Children(i);
   hhh.Units = 'normalized';
   hhh.FontSize = 8;
end






channelimglength = 9;
xi = 0.02;
yi = 0.025;
w = 0.15;
h = 0.28;
xspf = 0.1; %xspacefactor
yspf = 0.1;

dimm = [3 3];
xvc=[];
yvc=[];
for i = 1:dimm(1)
    xvi = [];
    yvi = [];
    for j = 1:dimm(2)
        xv = xi + (xspf*w)*(j-1) + (w*(j-1));
        xvi = horzcat(xvi,xv);
        
        yv = yi + (yspf*h)*(i-1) + (h*(i-1));
        yvi = horzcat(yvi,yv);
    end
    xvc = horzcat(xvc,xvi);
    yvc = horzcat(yvi,yvc);
end
% x = [xi xi+w+(xspf*w) xi+(w+(xspf*w)).*2 xi xi+w+(xspf*w) xi+(w+(xspf*w)).*2 xi xi+w+(xspf*w) xi+(w+(xspf*w)).*2];
% y = fliplr([yi yi yi yi+h+(yspf*w) yi+h+(yspf*w) yi+h+(yspf*w)   yi+(w+(yspf*w)).*2 yi+(w+(yspf*w)).*2 yi+(w+(yspf*w)).*2]);
x=xvc;
y=yvc;
for i=1:dimm(1)*dimm(2)
    ax= axes();
    ax.Position = [x(i) y(i) w h];
    ax.Units = 'inches';
    pos = ax.Position;
    %make them square
    if pos(4)>pos(3)
        pos(4) = pos(3);
    else
        pos(3) = pos(4);
    end
    ax.Position = pos;
    ax.Units = 'normalized';
    ax.XTick = [];
    ax.YTick = [];
    subaxes(i) = ax;
end

xi = 0.65;
yi = 0.025;
w = 0.15;
h = 0.28;
xspf = 0.1; %xspacefactor
yspf = 0.1;

dimm = [2 2];
xvc=[];
yvc=[];
for i = 1:dimm(1)
    xvi = [];
    yvi = [];
    for j = 1:dimm(2)
        xv = xi + (xspf*w)*(j-1) + (w*(j-1));
        xvi = horzcat(xvi,xv);
        
        yv = yi + (yspf*h)*(i-1) + (h*(i-1));
        yvi = horzcat(yvi,yv);
    end
    xvc = horzcat(xvc,xvi);
    yvc = horzcat(yvc,yvi);
end
x = xvc;
y = yvc;
for i=1:dimm(1)*dimm(2)
    ax= axes();
    ax.Position = [x(i) y(i) w h];
    ax.Units = 'inches';
    pos = ax.Position;
    %make them square
    if pos(4)>pos(3)
        pos(4) = pos(3);
    else
        pos(3) = pos(4);
    end
    ax.Position = pos;
    ax.Units = 'normalized';
    ax.XTick = [];
    ax.YTick = [];
    subaxestwo(i) = ax;
end



%define parameter structure and default parameter values
pStruct = defaultpStructFunc(segInstructList);

pStruct = loadSegmentParameters(pStruct,FileName,exportdir); %loads saved value of pStruct
f.Units='normalized';
f.Position =[0.1,0.2,0.8,0.7];
set(f,'KeyPressFcn',@keypress);
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
    disp('no saved parameters')
end

    

end

function saveSomethingCallback(~,~)
global  exportdir OGExpDate pStruct

cd(exportdir)
filename = strcat('*',OGExpDate,'*metaData*')
filelist = dir(filename)
savenamebase = char((filelist.name));
savename = strcat(OGExpDate,'-segmentParameters.mat');
save(savename,'pStruct');


end




function updateSliders
global pStruct ImageDetails sliderOne sliderOneTxt 

sliderx = 0.72;
sliderw = 0.1;
sliderh = 0.02;
slidertextw = 0.1;
sliderspace = 0.1;


channel = ImageDetails.segInstruct;

%nucDiameter
fnames = fieldnames(pStruct.(channel));
sliderspacing = linspace(0.8,0.7,length(fnames));
    for cyc = 1:length(fnames)
        str = fnames{cyc};
        if strcmp(str,'nucDiameter')
                minz = 1;
                maxz = 400;
                ssa = 1/5;%sliderStepAdjust
        elseif strcmp(str,'threshFactor')
            minz = 0.4;
            maxz = 3;
            ssa = 20;
        elseif strcmp(str,'sigmaScaledToParticle')
            minz = 1;
            maxz = 40;
            ssa = 1;
        elseif strcmp(str,'metthresh')
            minz = 0;
            maxz = 1;
            ssa = 20;
        elseif strcmp(str,'percentSmoothed')
            minz = 1;
            maxz = 100;
            ssa = 1/2.5;
        elseif strcmp(str,'denoise')
            minz = 5;
            maxz = 40;
            ssa = 1/2.5;
        else
            disp('parameter NOT CURRENTLY DEFINED') 
            if sum(strcmp(fnames,'metthresh'))<1
                str = 'metthresh';
                minz = 0;
                maxz = 1;
                ssa = 20;
                pStruct.(channel).(str) = 0.1;
            end
        end
        val.(str) = pStruct.(channel).(str);
        slidery = sliderspacing(cyc);

        
        sliderOne.(str) = uicontrol('Style', 'slider','String',str,'Min',minz,'Max',maxz,'SliderStep',[1 1]./((maxz-minz).*ssa),'Value',val.(str),'Position', [1 1 1 1],...
            'Callback', @sliderOneAdjust); 
        sliderOne.(str).Units='normalized';
        sliderOne.(str).Position = [sliderx slidery sliderw sliderh];
        sliderOneTxt.(str) = uicontrol('Style','text','Units','Normalized','Position',[1 1 1 1],'String',strcat(str,'=',num2str((val.(str)))));
        sliderOneTxt.(str).Units= 'Normalized';
        sliderOneTxt.(str).Position = [sliderx-sliderspace slidery slidertextw sliderh];  

    end

    
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
    channel = ImageDetails.segInstruct;
    str = source.String;
%     threshinput.(str) =source.Value;
%     zerostrel = round(source.Value);

    if strcmpi(str,'threshFactor')
        valupdate = source.Value;  
    elseif strcmpi(str,'metthresh')
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
    stringsToTest = {'rawMinusLPScaled','Inew','gradmag2','Ieg'};
    else
%     stringsToTest = {'rawMinusLPScaled','Ih','Ihcd','Shapes'};
    stringsToTest = {'imgRawDenoised','Inew','imgLowPass','Ih'};
    end
    for i = 1:length(subaxestwo)
    axes(subaxestwo(i))
    str = stringsToTest{i};
    img = testOut.(str);
    imagesc(img);t=title(str);
    t.FontSize=8;
    h=gca;
    h.XTick=[];
    h.YTick=[];
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
%     frames=1;
    [~,testOut] = segmentNuclei(img,nucleus_seg,pStruct,frames);
    
    IfFinal = false(size(FinalImage));
    for frames = 1:size(FinalImage,3)
        img = FinalImage(:,:,frames); 
        [If,~] = segmentNuclei(img,nucleus_seg,pStruct,frames);
        IfFinal(:,:,frames)=If;
    end
                
    %save here if running actual segmentation
end
function [IfFinal,testOut] = segmentationImageBackground(FinalImage,segmentPath,background_seg,backgroundFileName,pStruct)
testOut = struct();
    frames = 1;
    img = FinalImage(:,:,frames); 
    tic
    [~,testOut] = segmentCellBackgroundOLDGREAT(img,background_seg,pStruct,frames);
    toc
    tic
    [~,testOut] = segmentCellBackground(img,background_seg,pStruct,frames);
    toc
    
    IfFinal = false(size(FinalImage));
    for frames = 1:size(FinalImage,3)
        img = FinalImage(:,:,frames); 
        [If,~] = segmentCellBackground(img,background_seg,pStruct,frames);
        IfFinal(:,:,frames)=If;
    end
                
    %save here if running actual segmentation

disp('done')
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
global   ImageDetails Tracked SceneList




Tracked=[];

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


setSceneAndTime
end
function prevscenebutton_Callback(~,~) 
global   ImageDetails Tracked SceneList




Tracked=[];

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

function popup_menu_Callback_segment(source,~)
global ImageDetails Tracked channelinputs seginputs

Tracked=[];

% Determine the selected data set.
 str = source.String;
 val = source.Value;
 channel = char(str{val});

 [~,~,~,d] = regexpi(channel,channelinputs);
ImageDetails.Segment = d{1};

[~,~,~,d] = regexp(channel,seginputs);
ImageDetails.segInstruct = d{1};
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
global pStruct timeFrames mstackPath segInstructList channelList segmentList TC A  framesForDir ImageDetails  Tracked SceneList  SceneDirectoryPath imgfile imgsize

%   determine the channel directory
    cd(mstackPath)
    SceneDirectoryPath = mstackPath;
        if isempty(ImageDetails.Scene)
            ImageDetails.Scene = SceneList{1};
        end
        if isempty(ImageDetails.Channel)
            ImageDetails.Channel = channelList{1};
        end
        if isempty(ImageDetails.Segment)
            ImageDetails.Segment = segmentList{1};
        end
        if isempty(ImageDetails.segInstruct)
            ImageDetails.segInstruct = segInstructList{1};
        end

        
% determine the frame to load
    if isempty(ImageDetails.Frame)
       ImageDetails.Frame = 1;
    end
    t=1;




% define channel spacing

    channelspacing = round(linspace(1,timeFrames,9));
%     channelspacing(2) =69;
    if length(channelspacing)>timeFrames
        channelspacing = 1:min([timeFrames 9]);
    elseif max(channelspacing)>timeFrames 
        channelspacing = 1:min([timeFrames 9]);
    end
%     channelspacing = 1:timeFrames;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose the channel images
%[options are overlay of background // overlay of fluorescent channels // normal image]
        if strcmp(ImageDetails.Channel,'BKGbinary')%overlay background
            cd ..        
            cd('tiffs')
            ff = dir(strcat('*','EGFP','*'));
                if isempty(ff)
                    ff = dir(strcat('*','mKate','*'));
                end
            channelimg = single(loadUpTiffStackFrame(char(ff.name),t));
            prim = imdilate(bwperim(~logical(bkgimg)),strel('square',2));
            channelimg(prim) = max(max(channelimg));
        elseif strcmp(ImageDetails.Channel,'overlay')
            %need to write script for overlay

        else
            %load up channel image to view (not to be segemented)
            cd(mstackPath)
            ff = dir(strcat('*',ImageDetails.Scene,'*',ImageDetails.Channel,'*'));
            fname = char(ff.name);
            fileObject = matfile(fname);
            imgstack = single(fileObject.flatstack);
            channelimgstack = imgstack(:,:,channelspacing);
            img = imgstack(:,:,t);
            channelimg = img;

            %load up  channel image to segment
            cd(mstackPath)
            ff = dir(strcat('*',ImageDetails.Scene,'*',ImageDetails.Segment,'*'));


            fname = char(ff.name);
            fileObject = matfile(fname);
            imgstack = single(fileObject.flatstack);
            segmentimgstack = imgstack(:,:,channelspacing);
           
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


If = zeros(size(channelimg),'single');
FinalImage = segmentimgstack;
subdirname = SceneDirectoryPath;
scenename = ImageDetails.Scene;
filename = fname;
channel = ImageDetails.Channel;

nucleus_seg = ImageDetails.segInstruct;
background_seg = ImageDetails.segInstruct;
nucleusFileName = filename;
backgroundFileName = filename;
segmentPath = mstackPath;

if strcmp(ImageDetails.segInstruct,'background')
[IfStack,testOut] = segmentationImageBackground(FinalImage,segmentPath,background_seg,backgroundFileName,pStruct);  
% elseif strcmpi(ImageDetails.Segment,'Hoechst') && strcmp(ImageDetails.segInstruct,'nucleus')
% [IfStack,testOut] = segmentationNucleusHoechst(FinalImage,segmentPath,nucleus_seg,nucleusFileName,pStruct);
elseif strcmp(ImageDetails.segInstruct,'nucleus')
[IfStack,testOut] = segmentationNucleus(FinalImage,segmentPath,nucleus_seg,nucleusFileName,pStruct);
elseif strcmp(ImageDetails.segInstruct,'cell')
%need a segment cell script  
else
    error('no criteria met')
end


displayImageFunct(IfStack,channelimgstack,channelspacing);
updateSliders

plotTestOut(testOut,ImageDetails.segInstruct)




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
set(ttl,'FontSize',8);
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
ImageDetails.Segment=[];
ImageDetails.segInstruct = [];

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
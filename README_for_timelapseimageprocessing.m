%% set parent directory
mdir = mfilename('fullpath');
[~,b] = regexp(mdir,'Tracking\w*/');
if isempty(b)
    [~,b] = regexp(mdir,'Tracking\w*\');
end
parentdir = mdir(1:b);
exportdir = strcat(parentdir,'Export/');
cd(parentdir)

fullT = tic;

%% set experiments
dateArray = {'2018_08_26 plate c2c12 116 movie_immuno_4i immunoREST_Stain01_44bin_exp04',...
'2018_08_26 plate c2c12 116 movie_immuno_4i immunoREST_Stain02_44bin_exp07',...
'2018_08_26 plate c2c12 116 movie_immuno_4i immunoREST_Stain03_44bin_exp10',...
'2018_08_26 plate c2c12 116 movie_immuno_4i immunoREST_Stain04_44bin_exp15',...
'2018_08_26 plate c2c12 116 movie_immuno_4i immunoREST_Stain05_44bin_exp20',...
'2018_08_26 plate c2c12 116 movie_immuno_4i immunoREST_Stain06_44bin_exp24',...
'2018_08_26 plate c2c12 116 movie_immuno_4i immunoREST_Stain07_44bin_exp27'};

dateArray = {'2018_08_26 plate c2c12 116 movie_immuno_4i immunoREST_Stain03_44bin_focus_exp12',...
    '2018_08_26 plate c2c12 116 movie_immuno_4i immunoREST_Stain04_44bin_focus_exp16',...
    '2018_08_26 plate c2c12 116 movie_immuno_4i immunoREST_Stain05_44bin_focus_exp21'};

dateArray = {'2018_10_02 plate c2c12 116 magicRed rapamycin sapanisertib exp1'};

dateArray = {'2018_10_04 plate nmumg piggybac clones tgf long sub_clone240only exp2'};
dateArray = {'2018_10_13 plate c2c12 116 ampk activator rapamycin tgfbeta exp1'}
dateArray = {'2018_10_22 plate nmmg 235 5LH sort screen exp2'};
dateArray = {'2018_10_24 plate nmmg 235 5LH sort screen MOVIE exp2'}
dateArray = {'2018_10_30 plate c2c12 86 DALgreen autophagy fed starved rapamycin chloroquine exp1'}
dateArray = {'2018_10_30 plate c2c12 116 DALgreen rapamycin chloroquine tgfbeta movie exp1'};
dateArray = {'2018_11_01 plate c2c12 116 bix retinoic acid tgfbeta movie exp1'}

dateArray = {'2018_11_03 plate 2 c2c12 116 DALgreen rapamycin chloroquine tgf movie exp1'};
dateArray = {'2018_11_03 plate 1 c2c12 116 retinoic acid bix tgfbeta movie exp1','2018_11_03 plate 2 c2c12 116 DALgreen rapamycin chloroquine tgf movie exp1'};
dateArray = {'2018_11_03 plate 2 c2c12 116 retinoic acid tgfbeta movie exp2'};
dateArray = {'2018_11_05 plate 1 c2c12 116 retinoic acid tgfbeta movie exp1','2018_11_05 plate 2 c2c12 86 ctgf retinoic acid tgfbeta movie exp1'};
% dateArray = {'2018_11_05 plate 2 c2c12 86 ctgf retinoic acid tgfbeta movie exp1'};
dateArray = {'2018_10_30 plate c2c12 116 DALgreen rapamycin chloroquine tgfbeta movie exp1','2018_11_03 plate 2 c2c12 116 DALgreen rapamycin chloroquine tgf movie exp1'};
dateArray = {'2018_11_07 plate 1 c2c12 116 retinoic acid tgfbeta movie exp1'};
% dateArray = {'2018_11_09 plate 1 nmumg 235+157 ctgf tgfbeta movie exp1'};
dateArray = {'2018_11_07 plate 2 c2c12 116 snail 86 ctgf retinoic acid tgfbeta movie exp1'};
dateArray = {'2018_11_10 plate 1 c2c12 116 retinoic acid tgfbeta movie exp3'};
dateArray = {'2018_11_14 plate c2c12 133 caga12 retinoic acid tgfbeta movie exp1'};
dateArray = {'2018_08_18 plate c2c12 116 24 well movieimmumo immuno exp3'}
dateArray = {'2018_11_30 c2c12 new snail clones round1 plate exp1','2018_11_30 c2c12 new snail clones round1 plate exp2','2018_12_01 c2c12 new snail clones round2and3 plate exp1','2018_12_01 c2c12 new snail clones round2and3 plate exp2'};
dateArray = {'2017_11_04 plate exp1','2017_11_04 plate exp2'}
dateArray = {'2017_04_17 plate exp3','2017_04_17 plate exp4'}
%% specify which functions you which to run
extract_true    = 0;
excel_true      = 1;
flat_true       = 0;
seg_true        = 0;
autotrack_true  = 1;
bleach_true     = 0;

flatINDimg_true = 0;
extractTIFFflat_true    = 0; %extract tiff images for making imageJ movies
TIFFmovie_true    = 0; %extract tiff images for making imageJ movies
MATmovie_true    = 0; %extract tiff images for making imageJ movies

exportFrames = 0;



%% Extraction Loop .mat
if extract_true
    for BB = dateArray
        
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        %extract images from .czi and save image stacks as .mat files
        ExtractMetadataAndImages(B);
    end
end

%% Extraction Loop TIFFS
if extractTIFFflat_true
    for BB = dateArray
        
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        %extract images from .czi and save image stacks as .mat files
        ExtractMetadataAndImagesTIFF(B);
        A=[];
        BACKGROUND = [];
        channelstoinput = {'EGFP','DIC','mKate'};
        FlatfieldCorrectionOfTimeLapseImagesTIFFs(A,B,channelstoinput,BACKGROUND);
        
    end
end

% movies TIFFs
if TIFFmovie_true
    for BB = dateArray
        
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        %extract images from .czi and save image stacks as .mat files
        A=[];
        BACKGROUND = [];
        channelstoinput = {'EGFP','mKate'};
        stimulationFrame = 3;
        MovieMakerBaby([],B,channelstoinput,stimulationFrame)
        
    end
end

% movies TIFFs
if MATmovie_true
    for BB = dateArray
        
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        %extract images from .czi and save image stacks as .mat files
        A=[];
        BACKGROUND = [];
        channelstoinput = {'EGFP','mKate'};
        stimulationFrame = 2;
        MovieMakerBabyMAT([],B,channelstoinput,stimulationFrame)
        
    end
end




%% UserInput (Excel) Details Loop
if excel_true
    for BB = dateArray
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        %extract data in excel file
        cd(parentdir)
        exName = strcat(B,'-metaData.mat');
        makeDoseStructFromXLS(exName);
    end
end

%% Flatfield correction  Loop
if flat_true
    for BB = dateArray
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        % load excel-extracted data
        datequery = strcat(FileName,'*DoseAndScene*');
        cd(exportdir)
        filelist = dir(datequery);
        if isempty(filelist)
            error(strcat('need to run ExtractMetadata for-',FileName));
            %dosestruct = makeDoseStruct; %run function to make doseStruct
        else
            dosestructstruct = load(char(filelist.name));
            dosestruct = dosestructstruct.dosestruct;
            segInstruct = dosestructstruct.segInstruct;
        end
        A=[];
        channelstoinput = dosestructstruct.channelNameSwapArray;
        bkg = dosestructstruct.BACKGROUND;
        BACKGROUND = bkg{1};
        dontsegment = BACKGROUND;
        
        % run Flatflield correction
        BackgroundAndFlatfieldCorrectionOfTimeLapseImages(A,B,channelstoinput,BACKGROUND);
        
    end
end

%% BleachCorrection correction  Loop
if bleach_true
    for BB = dateArray
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        % load excel-extracted data
        datequery = strcat(FileName,'*DoseAndScene*');
        cd(exportdir)
        filelist = dir(datequery);
        if isempty(filelist)
            error(strcat('need to run ExtractMetadata for-',FileName));
            %dosestruct = makeDoseStruct; %run function to make doseStruct
        else
            dosestructstruct = load(char(filelist.name));
            dosestruct = dosestructstruct.dosestruct;
            segInstruct = dosestructstruct.segInstruct;
        end
        A=[];
        channelstoinput = dosestructstruct.channelNameSwapArray;
        bkg = dosestructstruct.BACKGROUND;
        BACKGROUND = bkg{1};
        dontsegment = BACKGROUND;
        
        % run Flatflield correction
        BleachCorrectionOfTimeLapseImages(A,B,channelstoinput,BACKGROUND);
        
    end
end




%% Flatfield correction  Loop
if flatINDimg_true
    for BB = dateArray
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        % load excel-extracted data
        datequery = strcat(FileName,'*DoseAndScene*');
        cd(exportdir)
        filelist = dir(datequery);
        if isempty(filelist)
            error(strcat('need to run ExtractMetadata for-',FileName));
            %dosestruct = makeDoseStruct; %run function to make doseStruct
        else
            dosestructstruct = load(char(filelist.name));
            dosestruct = dosestructstruct.dosestruct;
            segInstruct = dosestructstruct.segInstruct;
        end
        A=[];
        channelstoinput = dosestructstruct.channelNameSwapArray;
        bkg = dosestructstruct.BACKGROUND;
        BACKGROUND = bkg{1};
        dontsegment = BACKGROUND;
        
        % run Flatflield correction
%         BackgroundAndFlatfieldCorrectionOfTimeLapseImages(A,B,channelstoinput,BACKGROUND);
        Seg_Flat_and_Sub_TimeLapseImages(A,B,channelstoinput,BACKGROUND);


        
    end
end


%% Segmentation Loop

if seg_true
    for BB = dateArray
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        % load excel-extracted data
        datequery = strcat(FileName,'*DoseAndScene*');
        cd(exportdir)
        filelist = dir(datequery);
        if isempty(filelist)
            error(strcat('need to run ExtractMetadata for-',FileName));
            %dosestruct = makeDoseStruct; %run function to make doseStruct
        else
            dosestructstruct = load(char(filelist.name));
            dosestruct = dosestructstruct.dosestruct;
            segInstruct = dosestructstruct.segInstruct;
        end
        A=[];
        channelstoinput = dosestructstruct.channelNameSwapArray;
        bkg = dosestructstruct.BACKGROUND;
        BACKGROUND = bkg{1};
        dontsegment = BACKGROUND;
        
        % run segmentation
        %     uiSegmentTimeLapseImages
        SegmentationOfTimeLapseImages(A,B,dontsegment,segInstruct);
    end
end
%     donemail('cfrick@caltech.edu','segmentation complete','segmentation complete')

%% Autotracking and Export Loop
if autotrack_true
    for BB = dateArray
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        uiTrackCellzNEW(B,'AutoExportNuclei')
%         uiTrackCellzNEW(B,'AutoTrackCells')
%         uiTrackCellzNEW(B,'AutoExportTracks')
    end
end

%% export franems
if exportFrames
    for BB = dateArray
        
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        %extract images from .czi and save image stacks as .mat files
        A=[];
        BACKGROUND = [];
        channelstoinput = {'AlexaFluor647','Hoechst'};
        channelstoinput = {'mKate'};
        stimulationFrame = 2;
        frameExporter([],B,channelstoinput,stimulationFrame)
        
    end
end
fullTime = toc(fullT);
disp(['total time for all is = ', num2str(round(fullTime./60,0,'decimals')) ' minutes']);

% donemail('cfrick@caltech.edu','export complete','export complete')




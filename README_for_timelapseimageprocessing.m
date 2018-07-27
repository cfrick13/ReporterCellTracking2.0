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
dateArray = {'2018_02_05 plate exp1','2018_02_05 plate exp2'};
dateArray = {'2018_02_03 plate exp3'};
% dateArray = {'2016_08_01 plate crispr exp1'};
% dateArray = {'2018_03_12 plate exp2'};
% dateArray = {'2018_03_12 plate exp3'};
dateArray = {'2018_03_22 plate james exp1'};
dateArray = {'2018_03_16 plate exp1'};
dateArray = {'2018_03_30 plate james exp1'};
% dateArray = {'2018_03_31 plate exp1'};
% dateArray = {'2018_04_06 plate day3 c2c12 nmumg crispr tfect exp3'};
dateArray = {'2018_04_11 plate james mrna bmp tgf smad135 exp2'};
dateArray = {'2018_03_10 plate exp2'};
dateArray = {'2018_04_25 plate citrine scarlett james mrna exp1'};
% dateArray = {'2018_04_11 plate james mrna bmp tgf smad135 exp2'};

dateArray = {'2018_05_08 plate 116 1050fbs exp1','2018_05_08 plate 116 1050fbs exp2','2018_05_08 plate 116 1050fbs exp3'};
% dateArray = {'2018_05_19 plate c2c12 116 105250fbs exp1'};1
dateArray = {'2018_07_02 plate c2c12 116 inhibitors exp1'};
dateArray = {'2018_07_04 plate c2c12 116 inhibitors exp1'};
dateArray = {'2018_07_11 plate c2c12 116 inhibitors exp1'};
dateArray = {'2018_07_15 plate c2c12 116 bmp exp1'};
dateArray = {'2018_07_15 plate c2c12 116 bmp exp2'};
% dateArray = {'2018_07_17 plate c2c12 116 inhibitors titration exp1'};
% dateArray = {'2018_07_19 plate c2c12 116 post hoc synch exp1'};
dateArray = {'2018_07_25 plate c2c12 116 inhibitors exp1'};


%% specify which functions you which to run
extract_true    = 0;
excel_true      = 1;
flat_true       = 0;
seg_true        = 0;
autotrack_true  = 0;

flatINDimg_true = 0;
extractTIFFflat_true    = 0; %extract tiff images for making imageJ movies
TIFFmovie_true    = 0; %extract tiff images for making imageJ movies
MATmovie_true    = 0; %extract tiff images for making imageJ movies

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
        uiTrackCellzNEW(B,'AutoTrackCells')
        uiTrackCellzNEW(B,'AutoExportTracks')
    end
end

fullTime = toc(fullT);
disp(['total time for all is = ', num2str(round(fullTime./60,0,'decimals')) ' minutes']);

% donemail('cfrick@caltech.edu','export complete','export complete')




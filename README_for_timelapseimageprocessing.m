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


%% specify which functions you which to run
extract_true    = 0;
excel_true      = 1;
flat_true       = 0;
seg_true        = 0;
autotrack_true  = 0;


%% Extraction Loop
if extract_true
    for BB = dateArray
        
        %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)
        
        %time the whole process
        supertic = tic;
        
        %extract images from .czi and save image stacks as .mat files
        ExtractMetadataAndImages(B);
%         ExtractMetadataAndImagesTIFF(B);
%         A=[];
%         BACKGROUND = [];
%         channelstoinput = {'EGFP','DIC'};
%         FlatfieldCorrectionOfTimeLapseImagesTIFFs(A,B,channelstoinput,BACKGROUND);
        
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
        uiTrackCellzNEW(B,'AutoTrackCells')
        uiTrackCellzNEW(B,'AutoExportTracks')
        uiTrackCellzNEW(B,'AutoExportNuclei')
    end
end

fullTime = toc(fullT);
disp(['total time for all is = ', num2str(round(fullTime./60,0,'decimals')) ' minutes']);

% donemail('cfrick@caltech.edu','export complete','export complete')




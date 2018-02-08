%set parent directory
        mdir = mfilename('fullpath');
            [~,b] = regexp(mdir,'Tracking\w*/');
                if isempty(b)
                    [~,b] = regexp(mdir,'Tracking\w*\');
                end
        parentdir = mdir(1:b);
        exportdir = strcat(parentdir,'Export/');
    cd(parentdir)



fullT = tic;

%         dateArray = {'2017_06_22 plate exp1','2017_06_24 plate exp1','2017_06_24 plate exp2',...
%             '2017_06_24 plate exp3','2017_06_24 plate exp4','2017_06_26 plate exp2','2017_07_01 plate exp1',...
%             '2017_07_01 plate exp2','2017_07_03 plate exp1','2017_07_03 plate exp2','2017_07_05 plate exp1'};

        dateArray = {'2017_04_14 plate exp2','2017_04_12 plate exp8','2017_04_12 plate exp5',...
            '2017_04_17 plate exp2','2017_04_17 plate exp3','2017_04_17 plate exp4',...
            '2017_04_29 plate exp1','2017_05_01 plate exp1','2017_05_03 plate exp1','2017_05_03 plate exp2',...
            '2017_06_22 plate exp1','2017_06_24 plate exp1','2017_06_24 plate exp2',...
            '2017_06_24 plate exp3','2017_06_24 plate exp4','2017_06_26 plate exp2','2017_07_01 plate exp1',...
            '2017_07_01 plate exp2','2017_07_03 plate exp1','2017_07_03 plate exp2','2017_07_05 plate exp1',...
            '2017_07_07 plate exp1','2017_07_08 plate exp1','2017_07_14 plate exp1'};
        
                dateArray = {'2017_04_14 plate exp2','2017_04_12 plate exp8','2017_04_12 plate exp5',...
            '2017_04_17 plate exp2','2017_04_17 plate exp3','2017_04_17 plate exp4',...
            '2017_04_29 plate exp1','2017_05_01 plate exp1','2017_05_03 plate exp1','2017_05_03 plate exp2',...
            '2017_06_22 plate exp1','2017_06_24 plate exp1','2017_06_24 plate exp2',...
            '2017_06_24 plate exp3','2017_06_24 plate exp4','2017_06_26 plate exp2','2017_07_01 plate exp1',...
            '2017_07_01 plate exp2'};
        
        
        
        dateArray = {'2017_04_17 plate exp4'};
        dateArray = {'2014_10_01 plate exp1', '2014_09_30 plate exp1'};
        dateArray = {'2017_03_02 plate exp1'};
        dateArray = {'2016_10_08 plate exp3'};
        dateArray = {'2016_10_10 plate exp4','2016_10_10 plate exp5','2016_10_10 plate exp6','2016_10_10 plate exp7','2016_10_10 plate exp8'};
        dateArray = {'2016_10_12 plate exp4'};
        dateArray = {'2017_10_24 plate exp2'};
        dateArray = {'2017_10_25 plate exp1'};
        dateArray = {'2017_10_25 plate exp2'};
        dateArray = {'2017_10_27 plate exp1'};
        dateArray = {'2017_10_30 plate exp1'};
        dateArray = {'2017_10_30 plate exp2'};
        dateArray = {'2017_10_30 plate exp3'};
        dateArray = {'2017_11_01 plate exp1'};
        dateArray = {'2017_11_04 plate exp1'};
        dateArray = {'2017_11_04 plate exp2'};
        dateArray = {'2017_12_05 plate exp1'};
        dateArray = {'2017_12_05 plate exp2','2017_12_05 plate exp3'};
        dateArray = {'2017_12_09 plate exp1'};
        dateArray = {'2017_07_07 plate exp1','2017_07_08 plate exp1'};
        dateArray = {'2017_12_09 plate exp2'};
        dateArray = {'2018_01_02 plate exp1'};
        dateArray = {'2018_01_02 plate exp2'};
        dateArray = {'2018_01_08 plate exp1'};
        dateArray = {'2018_01_08 plate exp2'};
        dateArray = {'2018_01_06 plate exp1'};
        dateArray = {'2018_01_10 plate exp2'};
        dateArray = {'2018_01_10 plate exp1'};
        dateArray = {'2018_01_16 plate james exp1','2018_01_16 plate james exp2','2018_01_16 plate james exp3'};
        dateArray = {'2018_01_16 plate james exp3'};
        dateArray = {'2017_06_24 plate exp3'};
        dateArray = {'2018_01_16 plate exp1'};
        dateArray = {'2018_01_18 plate exp1'};
        dateArray = {'2018_02_03 plate exp1','2018_02_03 plate exp2'};
        dateArray = {'2018_02_03 plate exp3','2018_02_03 plate exp4'};
        dateArray = {'2018_02_03 plate exp4'};
%         dateArray = {'2018_01_18 plate exp2'};
%         dateArray = {'2018_01_20 plate exp1'};
%         dateArray = {'2018_01_20 plate exp2'};
%         dateArray = {'2017_07_01 plate exp1',...
%             '2017_07_01 plate exp2','2017_07_03 plate exp1','2017_07_03 plate exp2','2017_07_05 plate exp1'};
        
%         dateArray = {'2017_04_12 plate exp5','2017_04_29 plate exp1','2017_05_01 plate exp1','2017_05_03 plate exp1','2017_05_03 plate exp2'};
%         dateArray = {'2017_05_01 plate exp1','2017_05_03 plate exp1','2017_05_03 plate exp2'};
        
%         dateArray = {'2017_02_08 plate exp1','2017_02_08 plate exp2'};
%         dateArray = {'2017_04_17 plate exp2','2017_04_17 plate exp3','2017_04_17 plate exp4'};
%         dateArray = {'2017_04_17 plate exp4'};
%         dateArray = {'2017_07_01 plate exp2','2017_07_07 plate exp1','2017_07_08 plate exp1','2017_07_14 plate exp1'};
%         dateArray = {'2017_06_26 plate exp2'};
    for BB = dateArray

    %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)

    %time the whole process
        supertic = tic; 

    %% extract images from .czi and save image stacks as .mat files
%         ExtractMetadataAndImages(B);  
%         ExtractMetadataAndImagesTIFF(B);  

    %% extract data in excel file
        cd(parentdir)
        exName = strcat(B,'-metaData.mat');
        makeDoseStructFromXLS(exName);

    %% load excel-extracted data
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
%             FlatfieldCorrectionOfTimeLapseImagesTIFFs(A,B,channelstoinput,BACKGROUND);
    %% run segmentation
%         uiSegmentTimeLapseImages
        SegmentationOfTimeLapseImages(A,B,dontsegment,segInstruct);
%         donemail('cfrick@caltech.edu','segmentation complete','segmentation complete')

    %% run autotracking algorithms
    %
    uiTrackCellz(B,'AutoTrackCells')
    uiTrackCellz(B,'AutoExportTracks')
    uiTrackCellz(B,'AutoExportNuclei')
        
    end
fullTime = toc(fullT);
disp(['total time for all is = ', num2str(round(fullTime./60,0,'decimals')) ' minutes']);

donemail('cfrick@caltech.edu','export complete','export complete')




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

sbAndWashList = {'2017_04_12 plate exp5 mKate Air objective';'2017_04_12 plate exp8 mKate Air objective';'2017_04_14 plate exp2 mKate Air objective';'2017_04_29 plate exp1 mKate Air objective';'2017_05_01 plate exp1 mKate Air objective';'2017_05_03 plate exp2 mKate Air objective';'2017_06_24 plate exp1 mKate Air objective';'2017_07_01 plate exp2 40x Oil objective';'2017_07_03 plate exp1 40x Oil objective';'2017_07_03 plate exp2 40x Oil objective';'2017_07_05 plate exp1 40x Oil objective';'2017_07_07 plate exp1 40x Oil objective';'2017_07_08 plate exp1 40x Oil objective';'2017_10_23 plate exp1 40x Oil objective';'2017_10_24 plate exp1 40x Oil objective';'2017_10_24 plate exp2 40x Oil objective';'2017_10_25 plate exp1 40x Oil objective';'2017_10_27 plate exp1 40x Oil objective';'2017_10_30 plate exp1 40x Oil objective';'2017_10_30 plate exp2 40x Oil objective';'2017_10_30 plate exp3 40x Oil objective';'2017_11_01 plate exp1 40x Oil objective'};
sbAndWashList = {'2017_04_17 plate exp2 40x Oil objective';'2017_04_17 plate exp3 40x Oil objective';'2017_04_17 plate exp4 40x Oil objective';'2017_05_03 plate exp1 40x Oil objective';'2017_06_22 plate exp1 20x Air objective';'2017_06_24 plate exp2 40x Oil objective';'2017_06_24 plate exp3 40x Oil objective';'2017_06_24 plate exp4 20x Air objective';'2017_07_01 plate exp1 40x Oil objective';'2017_10_25 plate exp2 40x Oil objective';'2017_11_04 plate exp1 40x Oil objective';'2017_11_04 plate exp2 40x Oil objective'};
sbAndWashList = {'2017_07_08 plate exp1'};
sbAndWashList = {'2018_01_29 plate james exp2'};
for j = 1:length(sbAndWashList)
    close all
%     j=7;
    wstr = sbAndWashList{j};
    [~,b] = regexp(wstr,'exp');
    BB = wstr(1:b+1);
%     disp(FileName)
    
    

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
%         makeDoseStructFromXLS(exName);

    %% load excel-extracted data
%         datequery = strcat(FileName,'*DoseAndScene*');
%         cd(exportdir)
%         filelist = dir(datequery);
%             if isempty(filelist)
%                 error(strcat('need to run ExtractMetadata for-',FileName));
%                 %dosestruct = makeDoseStruct; %run function to make doseStruct 
%             else
%                 dosestructstruct = load(char(filelist.name));
%                 dosestruct = dosestructstruct.dosestruct;
%                 segInstruct = dosestructstruct.segInstruct;
%             end
%         A=[];
%         channelstoinput = dosestructstruct.channelNameSwapArray;
%         bkg = dosestructstruct.BACKGROUND;
%         BACKGROUND = bkg{1};
%         dontsegment = BACKGROUND;

%% run Flatflield correction
%         BackgroundAndFlatfieldCorrectionOfTimeLapseImages(A,B,channelstoinput,BACKGROUND);
            channelstoinput = {'DIC','mKate','EGFP','Hoechst','CFP','AlexaFluor647'};
            BACKGROUND = 1:5;
            A=[];
            FlatfieldCorrectionOfTimeLapseImagesTIFFs(A,B,channelstoinput,BACKGROUND);
    %% run segmentation
%         uiSegmentTimeLapseImages
%         SegmentationOfTimeLapseImages(A,B,dontsegment,segInstruct);
%         donemail('cfrick@caltech.edu','segmentation complete','segmentation complete')

    %% run autotracking algorithms
    %
%     uiTrackCellz(B,'AutoTrackCells')
%     uiTrackCellz(B,'AutoExportTracks')
%     uiTrackCellz(B,'AutoExportNuclei')
        
    end
fullTime = toc(fullT);
disp(['total time for all is = ', num2str(round(fullTime./60,0,'decimals')) ' minutes']);

% donemail('cfrick@caltech.edu','export complete','export complete')




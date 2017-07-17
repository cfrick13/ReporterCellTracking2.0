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
            '2017_04_29 plate exp1','2017_05_01 plate exp1','2017_05_03 plate exp1','2017_05_03 plate exp2',...
            '2017_06_22 plate exp1','2017_06_24 plate exp1','2017_06_24 plate exp2',...
            '2017_06_24 plate exp3','2017_06_24 plate exp4','2017_06_26 plate exp2','2017_07_01 plate exp1',...
            '2017_07_01 plate exp2','2017_07_03 plate exp1','2017_07_03 plate exp2','2017_07_05 plate exp1'};
        
        dateArray = {'2017_07_01 plate exp1',...
            '2017_07_01 plate exp2','2017_07_03 plate exp1','2017_07_03 plate exp2','2017_07_05 plate exp1'};
        
%         dateArray = {'2017_04_12 plate exp5','2017_04_29 plate exp1','2017_05_01 plate exp1','2017_05_03 plate exp1','2017_05_03 plate exp2'};
%         dateArray = {'2017_05_01 plate exp1','2017_05_03 plate exp1','2017_05_03 plate exp2'};
        
%         dateArray = {'2017_02_08 plate exp1','2017_02_08 plate exp2'};
        dateArray = {'2017_06_26 plate exp2'};
    for BB = dateArray

    %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)

    %time the whole process
        supertic = tic; 

    %% extract images from .czi and save image stacks as .mat files
%         ExtractMetadataAndImages(B);  

    %% extract data in excel file
        cd(parentdir)
        exName = strcat(B,'-metaData.mat');
%         makeDoseStructFromXLS(exName);

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
    
        
        
    %% run Flatflield correction
%         BackgroundAndFlatfieldCorrectionOfTimeLapseImages(A,B,channelstoinput,BACKGROUND);

    %display time
        totalTime = toc(supertic); %report the timing
%         disp(strcat('total time from extract to flat is=', num2str(totalTime./60),' minutes'));


    
    %% run segmentation
%         uiSegmentTimeLapseImages
%         SegmentationOfTimeLapseImages(A,B,dontsegment,segInstruct);
%         donemail('cfrick@caltech.edu','segmentation complete','segmentation complete')

    %% run autotracking algorithms
%         AutoExportNuclei(B)

        uiTrackCellz(B,'AutoExportNuclei')
        AutoTrackCellz(B)
        uiTrackCellz(B,'AutoTrackCellz')
        
    end
fullTime = toc(fullT);
disp(['total time for all is = ', num2str(fullTime./3600),' hours']);

donemail('cfrick@caltech.edu','export complete','export complete')



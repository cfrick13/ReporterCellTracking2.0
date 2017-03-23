
%*files must be exported as TIFFs from Zen using channel names
%*the export folder must have the same name as the parent folder containing the zen experiment files
% such as parent folder = 'D:/kkim/2016_10_25 bcat cells' then export folder should be 'D:/kkim/2016_10_25 bcat cells/2016_10_25 bcat cells'

    %set parent directory
        mdir = mfilename('fullpath');
            [~,b] = regexp(mdir,'Tracking\w*/');
                if isempty(b)
                    [~,b] = regexp(mdir,'Tracking\w*\');
                end
        parentdir = mdir(1:b);
        exportdir = strcat(parentdir,'Export/');
    cd(parentdir)

    % dateArray = {'2017_03_02 plate exp1'};
%     dateArray = {'2017_03_13 plate exp5','2014_09_27 plate exp1','2017_09_30 plate exp1','2014_10_01 plate exp1'};
%     dateArray = {'2014_09_27 plate exp1','2014_09_30 plate exp1','2014_10_01 plate exp1','2014_10_04 plate exp1','2014_11_21 plate exp1'};
%     dateArray = {'2017_03_20 plate exp1'};
    dateArray = {'2017_01_30 plate exp1','2017_01_30 plate exp2','2017_03_13 plate exp3'};
    for BB = dateArray

    %get directory name from date Array
        B = char(BB);
        FileName = B;
        disp(FileName)

    %time the whole process
        supertic = tic; 

    %extract images from .czi and save image stacks as .mat files
%         ExtractMetadataAndImages(B);    



    %extract data in excel file
        exName = strcat(B,'-metaData.mat');
%         makeDoseStructFromXLS(exName);


    %load excel-extracted data
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
    

    %run Flatflield correction
%         BackgroundAndFlatfieldCorrectionOfTimeLapseImages(A,B,channelstoinput,BACKGROUND);

    %display time
        totalTime = toc(supertic); %report the timing
        disp(strcat('total time from extract to flat is=', num2str(totalTime./60),' minutes'));




        SegmentationOfTimeLapseImages(A,B,dontsegment,segInstruct);
    end





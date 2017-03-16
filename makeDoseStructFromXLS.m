function dosestruct = makeDoseStructFromXLS(savename)
close all
mdir = mfilename('fullpath');
    [~,b] = regexp(mdir,'Tracking\w*/');
        if isempty(b)
            [~,b] = regexp(mdir,'Tracking\w*\');
        end
parentdir = mdir(1:b);
exportdirz = strcat(parentdir,'Export/');
cd(exportdirz);



% FileName = uigetfile('*metaData.mat');%choose file to load
load(savename)
% load(FileName)
FileName = savename;
[~,b] = regexp(FileName,'exp[0-9]');
queryname = strcat(FileName(1:b),'*experimentInfo.xlsx');
filelist = dir(queryname);
if length({filelist.name}) ==1
    [num,txt] = xlsread(char(filelist.name));
else
    error('you need to make a spreadsheet with the following information')
    
% disp(raw)
%     'NumberOfDoses'    'Doses'     'Scenes'    'NumOfTgfAdditions'    'TgfFrames'    'Scenes'    'FlatfieldReference'
%     [            5]    [2.4000]    '1:6'       [                1]    [        8]    '1:31'      '17:21'             
%     [          NaN]    [0.1400]    '7:11'      [              NaN]    [      NaN]    [   NaN]    [               NaN]
%     [          NaN]    [0.0700]    '27:31'     [              NaN]    [      NaN]    [   NaN]    [               NaN]
%     [          NaN]    [0.0400]    '12:16'     [              NaN]    [      NaN]    [   NaN]    [               NaN]
%     [          NaN]    [     0]    '17:26'     [              NaN]    [      NaN]    [   NaN]    [               NaN]
    
end


dosestruct = struct();
numberOfScenes = datastruct.sceneCount;
for i=1:numberOfScenes
   sceneStr = 's00';
   scene = num2str(i);
   sceneStr(end-(length(scene)-1):end)=scene;
   dosestruct(i).scene = sceneStr;
end


% % % % % % % % % % % % %input the number of doses
numberOfDoses = num(1,1);

% % % % % % % %input the Tgfbeta concentrations of each dose
doses = num(1:numberOfDoses,2);


%input numbers for the scene numbers corresponding to each dose
doseToScenez = txt(2:numberOfDoses+1,3);
doseToScene = cellfun(@str2num,doseToScenez,'UniformOutput',0);
%output is a 1,n cell array of matrices corresponding to the scenes for the


%input dose information into the structure
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored
coloringArray = vertcat({dosestruct.(coloringChoice)});

for j=1:length(doseToScene)
    doseToScenemat = doseToScene{j};
    doseToSceneArray=cell(1,length(doseToScenemat));
    for i = 1:length(doseToScenemat)
        dosestr = num2str(doseToScenemat(i)); 
        if length(dosestr)>1
            doseToSceneArray{i} = strcat('s',dosestr);
        else
            doseToSceneArray{i} = strcat('s0',dosestr); 
        end
    end
    
    

indicesChoice =channelregexpmaker(doseToSceneArray);
[~,~,~,d] =  regexp(coloringArray,indicesChoice);
dmat = cellfun(@isempty,d);
indices = ~dmat;
[dosestruct(indices).dose] = deal(doses(j));
[dosestruct(indices).dosestr] = deal(num2str(doses(j)));
end



%%%%%%%%%%%%%%%%%%%%%%  number of Tgfbeta additions %%%%%%%%%%%%%%%%%%
tgfadditions = num(1,4);

%input numbers for the frame(s) after which tgfbeta was added
tgfFrames = num(1:tgfadditions,5);

%input numbers for the scenes for each tgfbeta addition
tgfScenez = txt(2:tgfadditions+1,6);
tgfScenes = cellfun(@str2num,tgfScenez,'UniformOutput',0);
%output is a 1,n cell array of matrices corresponding to the scenes for the
%n doses

%input dose information into the structure
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored
coloringArray = vertcat({dosestruct.(coloringChoice)});

for j=1:length(tgfScenes)
    tgfToScenemat = tgfScenes{j};
    tgfToSceneArray=cell(1,length(tgfToScenemat));
    for i = 1:length(tgfToScenemat)
        dosestr = num2str(tgfToScenemat(i)); 
        if length(dosestr)>1
            tgfToSceneArray{i} = strcat('s',dosestr);
        else
            tgfToSceneArray{i} = strcat('s0',dosestr); 
        end
    end
    
    

indicesChoice =channelregexpmaker(tgfToSceneArray);
[~,~,~,d] =  regexp(coloringArray,indicesChoice);
dmat = cellfun(@isempty,d);
indices = ~dmat;
[dosestruct(indices).tgfFrame] = deal(tgfFrames(j));
[dosestruct(indices).tgfFramestr] = deal(num2str(tgfFrames(j)));
end



%%%%%%%%%%%%%%%%%%%%%%  reference frames %%%%%%%%%%%%%%%%%%
% BACKGROUND = cellfun(@str2num,inputdlgOutput,'UniformOutput',0);
BACKGROUNDz = txt(2,7);
BACKGROUND = cellfun(@str2num,BACKGROUNDz,'UniformOutput',0);

for j=1:length(BACKGROUND)
    tgfToScenemat = BACKGROUND{j};
    tgfToSceneArray=cell(1,length(tgfToScenemat));
    for i = 1:length(tgfToScenemat)
        dosestr = num2str(tgfToScenemat(i)); 
        if length(dosestr)>1
            tgfToSceneArray{i} = strcat('s',dosestr);
        else
            tgfToSceneArray{i} = strcat('s0',dosestr); 
        end
    end
    
    

indicesChoice =channelregexpmaker(tgfToSceneArray);
[~,~,~,d] =  regexp(coloringArray,indicesChoice);
dmat = cellfun(@isempty,d);
indices = ~dmat;
[dosestruct(indices).flatfield] = deal('BACKGROUND');
[dosestruct(dmat).flatfield] = deal('experiment');
end


if size(txt,2)>8
segInstruct.nucleus = txt{2,8};
segInstruct.cell = txt{2,9};
segInstruct.background = txt{2,10};
end

exportdir = exportdirz;
[~,b] = regexp(FileName,'exp[0-9]');
savename = FileName(1:b);
cd(exportdir)
savename = strcat(savename,'-DoseAndSceneData.mat');
save(savename)


 

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
    
    if strcmp(')',channelinputs(end))
    else
        channelinputs(end+1)=')';
    end
end


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


FileName = savename;
% FileName = uigetfile('*metaData.mat');%choose file to load
load(savename)
% load(FileName)
% FileName = savename;
[~,b] = regexp(FileName,'exp[0-9]');
queryname = strcat(FileName(1:b),'*experimentInfo.xlsx');
filelist = dir(queryname);
if length({filelist.name}) ==1
    [num,txt,raw] = xlsread(char(filelist.name));
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
    if numberOfScenes>99
        sceneStr = 's000';
    else
        sceneStr = 's00';
    end
   scene = num2str(i);
   sceneStr(end-(length(scene)-1):end)=scene;
   dosestruct(i).scene = sceneStr;
end


% % % % % % % % % % % % %input the number of doses
numberOfDoses = num(1,7);

% % % % % % % %input the Tgfbeta concentrations of each dose
doses = raw(2:numberOfDoses+1,8);

%input numbers for the scene numbers corresponding to each dose
doseToScenez = raw(2:numberOfDoses+1,9);
doseToScene = cellfun(@str2num,doseToScenez,'UniformOutput',0);
%output is a 1,n cell array of matrices corresponding to the scenes for the


% % % % % % % %input the number of wells imaged
numberOfWells = num(1,1);
wells = raw(2:numberOfWells+1,2);

%input numbers for the scene numbers corresponding to each well
wellsToScenez = raw(2:numberOfWells+1,3);
wellsToScene = cellfun(@str2num,wellsToScenez,'UniformOutput',0);
%output is a 1,n cell array of matrices corresponding to the scenes for the

% % % % % % % %input the number of wells imaged
numberOfConditions = num(1,4);
conditions = raw(2:numberOfConditions+1,5);

%input numbers for the scene numbers corresponding to each well
conditionsToScenez = raw(2:numberOfConditions+1,6);
conditionsToScene = cellfun(@str2num,conditionsToScenez,'UniformOutput',0);
%output is a 1,n cell array of matrices corresponding to the scenes for the


%%%%%%%%%%%%%%%%%%%%%%  number of Tgfbeta additions %%%%%%%%%%%%%%%%%%
tgfadditions = num(1,10);

%input numbers for the frame(s) after which tgfbeta was added
tgfFrames = num(1:tgfadditions,11);

%input numbers for the scenes for each tgfbeta addition
tgfScenez = raw(2:tgfadditions+1,12);
tgfToScene = cellfun(@str2num,tgfScenez,'UniformOutput',0);
%output is a 1,n cell array of matrices corresponding to the scenes for the
%n doses

%input dose information into the structure
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored
sceneList = vertcat({dosestruct.(coloringChoice)});

%input dose information into the structure
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored
sceneList = vertcat({dosestruct.(coloringChoice)});



%deal doses
for j=1:length(doseToScene)
    doseToScenemat = doseToScene{j};
    doseToSceneArray=cell(1,length(doseToScenemat));
    for i = 1:length(doseToScenemat)
        dosestr = num2str(doseToScenemat(i));
        
        %
        
        msv = numberOfScenes;
        if msv>99
            backstr = 's000';
        else
            backstr = 's00';
        end
        backstr(end-(length(dosestr)-1):end) = dosestr;
        doseToSceneArray{i} = backstr;
    end

    
    

indicesChoice =channelregexpmaker(doseToSceneArray);
[~,~,~,d] =  regexp(sceneList,indicesChoice);
dmat = cellfun(@isempty,d);
indices = ~dmat;

if iscellstr(doses(j))
   [dosestruct(indices).dose] = deal(doses{j});
   [dosestruct(indices).dosestr] = deal(doses{j});
elseif iscell(doses(j))
   [dosestruct(indices).dose] = deal(doses{j});
   [dosestruct(indices).dosestr] = deal(num2str(doses{j}));
else
    [dosestruct(indices).dose] = deal(doses(j));
    [dosestruct(indices).dosestr] = deal(num2str(doses(j)));
end

end


%deal Tgf
var = tgfFrames;
varToScene = tgfToScene;
str = 'tgfFrame';
dosestruct = dealData(var,varToScene,dosestruct,str,sceneList);

%deal Wells
var = wells;
varToScene = wellsToScene;
str = 'wells';
dosestruct = dealData(var,varToScene,dosestruct,str,sceneList);

%deal Conditions
var = conditions;
varToScene = conditionsToScene;
str = 'conditions';
dosestruct = dealData(var,varToScene,dosestruct,str,sceneList);


%%%%%%%%%%%%%%%%%%%%%%  reference frames %%%%%%%%%%%%%%%%%%
% BACKGROUND = cellfun(@str2num,inputdlgOutput,'UniformOutput',0);
BACKGROUNDz = raw(2,13);
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
[~,~,~,d] =  regexp(sceneList,indicesChoice);
dmat = cellfun(@isempty,d);
indices = ~dmat;
[dosestruct(indices).flatfield] = deal('BACKGROUND');
[dosestruct(dmat).flatfield] = deal('experiment');
end


if size(txt,2)>13
segInstruct.nucleus = txt{2,14};
segInstruct.cell = txt{2,15};
segInstruct.background = txt{2,16};
end

exportdir = exportdirz;
[~,b] = regexp(FileName,'exp[0-9]');
savename = FileName(1:b);
cd(exportdir)
savename = strcat(savename,'-DoseAndSceneData.mat');
save(savename)


 

end



function dosestruct = dealData(var,varToScene,dosestruct,str,sceneList)
    for j=1:length(varToScene)
        varToScenemat = varToScene{j};
        varToSceneArray=cell(1,length(varToScenemat));
        msv = length(sceneList);
        for i = 1:length(varToScenemat)
            if msv>99
                backstr = 's000';
            else
                backstr = 's00';
            end
            dosestr = num2str(varToScenemat(i));
            backstr(end-(length(dosestr)-1):end) = dosestr;
            varToSceneArray{i} = backstr;
        end
        
    indicesChoice =channelregexpmaker(varToSceneArray);
    [~,~,~,d] =  regexp(sceneList,indicesChoice);
    dmat = cellfun(@isempty,d);
    indices = ~dmat;
    [dosestruct(indices).(str)] = deal(var(j));
        if isnumeric(var(j))
            [dosestruct(indices).([str 'str'])] = deal(num2str(var(j)));
        end
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
    
    if strcmp(')',channelinputs(end))
    else
        channelinputs(end+1)=')';
    end
end


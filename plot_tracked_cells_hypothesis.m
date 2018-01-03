function plot_tracked_cells_hypothesis(FileName,timeCutoff,sortTime,cmapchoice,groupBasedOn,colorBasedOn,subplotGroupBasedOn,labelBasedOn,basalLength,smadnormalizestr,reporternormalizestr,includeStr,excludeStr,SmadTotalOrMedian,ReporterTotalOrMedian,celltracestr)
% close all
%use grabcode to retreive code

%important parameters
%   medfiltnum
%   jittercut
medfiltnum = 3;

%determine the location of the matlab function and establish export
%directory in relation to that filepath
mdir = mfilename('fullpath');
[~,b] = regexp(mdir,'Tracking\w*/');
if isempty(b)
    [~,b] = regexp(mdir,'Tracking\w*\');
end
parentdir = mdir(1:b); %specifies folder in which all analysis is being done


addpath([parentdir 'Colormaps'])
loaddir = strcat(parentdir,'Export'); %specifies where data is exported
cd(loaddir);

[~,b] = regexp(mdir,'/');
if isempty(b)
    [~,b] = regexp(mdir,'\');
end

mfiledir =mdir(1:b(end)); %specifies location of matlab function file

%load the exported tracking structure
cd(loaddir)
load(FileName)

%load metadata associated with the experiment (requires manual input if there is ambiguity)
[a,~] = regexp(FileName,'_tracking');
fileDateStr = FileName(1:a-1);
datequery = strcat(fileDateStr,'*metaData.mat');
cd(loaddir)
filelist = dir(datequery);
if length({filelist.name}) ==1
    metaData = load(char(filelist.name));
else
    error('need to extract metadata')
    filename = uigetfile();
    metaData = load(filename);
end
%determine the timeVector from metaData [dim 1 is scene#, dim 2 is each time frame]
timeVec = metaData.timeVec;


%load information regarding doses and scenes and tgfbeta addition
datequery = strcat(fileDateStr,'*DoseAndScene*');
cd(loaddir)
filelist = dir(datequery);
if isempty(filelist)
    error('run makeDoseStruct')
    dosestruct = makeDoseStruct; %run function to make doseStruct
else
    dosestructstruct = load(char(filelist.name));
    dosestruct = dosestructstruct.dosestruct;
end


%determine the scenes present in the experiment in order to combine metadata
scenestr = 'scene';
sceneListArray = vertcat({exportStruct.(scenestr)});
sceneList = unique(sceneListArray);
sceneListArrayTwo = vertcat({dosestruct.(scenestr)});

%combine the exportStruct information with dosestruct information
for i=1:length(sceneList)
    sceneChoice=sceneList{i};
    indices = strcmp(sceneListArray,sceneChoice);
    indicestwo = strcmp(sceneListArrayTwo,sceneChoice);
    
    fnamesz = fieldnames(dosestruct);
    for p = 1:length(fnamesz)
        fnstr = fnamesz{p};
        var = dosestruct(indicestwo).(fnstr);
        if isnumeric(var)
            [exportStruct(indices).(fnstr)] = deal(var);
        else
            [exportStruct(indices).(fnstr)] = deal(char(var));
        end
    end
end

%combine the exportStruct information with dosestruct information
for i=1:length(sceneList)
    sceneChoice=sceneList{i};
    indices = strcmp(sceneListArray,sceneChoice);
    indicestwo = strcmp(sceneListArrayTwo,sceneChoice);
    
    fnamesz = fieldnames(dosestruct);
    for p = 1:length(fnamesz)
        fnstr1 = 'dosestr';
        fnstr2 = 'conditions';
        fnstr3 = 'expdate';
        var1 = dosestruct(indicestwo).(fnstr1);
        var2 = dosestruct(indicestwo).(fnstr2);
        [exportStruct(indices).doseAndCondition] = deal([char(var1) ' ' char(var2)]);
        [exportStruct(indices).doseconddate] = deal([char(var1) ' ' fileDateStr ' ' char(var2)]);
        [exportStruct(indices).conddate] = deal([fileDateStr ' ' char(var2)]);
        [exportStruct(indices).(fnstr3)] = deal(fileDateStr);
    end
end


%determine details needed for plotting such as when Tgfbeta is added, etc
stimulationFrame = exportStruct(1).tgfFrame;
numberOfFrames = size(timeVec,2);
finalFrame = numberOfFrames;


%need to determine the number of scenes present and choose the time vector
%depending on the scene from which it was imaged
numberOfCells = length(indices);
timeMatrix = zeros(numberOfCells,finalFrame);

for i=1:numberOfCells
    sceneChoice=exportStruct(i).scene;
    idxtwo = strcmp(sceneListArrayTwo,sceneChoice);% sceneListArrayTwo = vertcat({dosestruct.(scenestr)});
    stimulationFrame = dosestruct(idxtwo).tgfFrame+1;
    timeMatrix(i,:) = timeVec(idxtwo,1:finalFrame)-timeVec(1,stimulationFrame); %subtract by time closest to Tgfbeta addition
    exportStruct(i).timeMatrix = timeMatrix(i,:);
    exportStruct(i).cellID = num2str(exportStruct(i).cellID);
    exportStruct(i).cellNum = num2str(i);
end





%% group based on desired criteria, set colormap, etc
%choose all cells that match specified chosen clone input strings
allList = {exportStruct.('doseconddate')};
[~,~,~,d] = regexp(allList,includeStr);
subidx = ~cellfun(@isempty,d,'UniformOutput',1);
exportStruct = exportStruct(subidx);

%exclude all cells that match specified exludeNames input strings
allList = {exportStruct.('doseconddate')};
[~,~,~,d] = regexp(allList,excludeStr);
subidx = ~cellfun(@isempty,d,'UniformOutput',1);
exportStruct = exportStruct(~subidx);




%% extract cell traces and plot
%function to exract the cell traces, normalized and not
indices = true(1,length(exportStruct)); %choose all cells initially
cTracks = struct();
cTracks.Smad = extractTraces(exportStruct,indices,'NucEGFP',finalFrame,stimulationFrame,basalLength,medfiltnum);
cTracks.Reporter = extractTraces(exportStruct,indices,'NucmKate',finalFrame,stimulationFrame,basalLength,medfiltnum);



if strcmp(SmadTotalOrMedian,'nuccytoXabundance')
    smadtraces1 = cTracks.('Smad').('nuccyto').('none');
    smadtraces2 = cTracks.('Smad').('median').('abundance');
    smadtraces = smadtraces1.*smadtraces2;
else
    smadtraces = cTracks.('Smad').(SmadTotalOrMedian).('none');
end

smadtraces1 = cTracks.('Smad').('nuccyto').('none');
smadtraces2 = cTracks.('Smad').('median').('abundance');
smadtracesNC = smadtraces1.*smadtraces2;

reportertraces = cTracks.('Reporter').(ReporterTotalOrMedian).('none');


%% division identifying and interpolating function
cellMat = struct();
[tracemat,divtracemat,allmat,nanidx] = divisionFunction(cTracks,smadtraces,stimulationFrame,timeMatrix,SmadTotalOrMedian);
cellMat.Smad.division = divtracemat;
cellMat.Smad.traces = tracemat;
cellMat.Smad.all = allmat;
nc = cTracks.('Smad').('nuccyto').('none');
cellMat.Smad.NC = nc(nanidx,:);
[tracemat,divtracemat,allmat,~] = divisionFunction(cTracks,reportertraces,stimulationFrame,timeMatrix,ReporterTotalOrMedian);
cellMat.Reporter.division = divtracemat;
cellMat.Reporter.traces = tracemat;
cellMat.Reporter.all = allmat;
nc = cTracks.('Reporter').('nuccyto').('none');
cellMat.Reporter.NC = nc(nanidx,:);


%% update the color values based on normalized or not
exportSubStruct = exportStruct(nanidx);
timeMatrix = timeMatrix(nanidx,:);

%% set cutoff based on time
tvec = timeMatrix(1,:);
timeidx = find(tvec>timeCutoff);
celltraces = cellMat.Smad.all;
celltime = celltraces(:,timeidx(1));
timenan = ~isnan(celltime);

exportStruct = exportSubStruct;
exportSubStruct = exportStruct(timenan);
newCellMat = updatealltraces(cellMat,timenan);
cellMat = newCellMat;
timeMatrix = timeMatrix(timenan,:);


%% sort traces
sortstr = 'Smad';
[sortedCellMat,sortidx] = sorttraces(cellMat,sortstr,timeMatrix,sortTime);
cellMat = sortedCellMat;

%% update the color values based on sorted or not
exportStruct = exportSubStruct;
timeMatrixOld = timeMatrix;
for i = 1:length(sortidx)
    exportSubStruct(i) = exportStruct(sortidx(i));
    timeMatrix(i,:) = timeMatrixOld(sortidx(i),:);
end

%% establish lists and listArrays for coloring and labeling
%establish the color map for plotting
%choose which field based upon which each cell trace will get colored
if strcmp(colorBasedOn,'individually')
    coloringChoice = 'scene';
else
    coloringChoice = colorBasedOn; %color traces based on what
end
[coloringList,coloringArrayList] = listSpitter(exportSubStruct,coloringChoice);
cmap = determineCMAP(cmapchoice,coloringList);

%assign a color array using the created colormap based on the choices above
colormapindividualMatrix = zeros(length(coloringArrayList),size(cmap,2));
for i=1:length(coloringArrayList)
    cA = coloringArrayList{i};
    idx = strcmp(coloringList,cA);
    colormapindividualMatrix(i,:) = cmap(idx,:); %this is used to color each individual trace according to the colorBasedOn criteria selected
end
% randmat = normrnd(0,0.1,size(colormapindividualMatrix,1),3);
% colormapindividualMatrix = colormapindividualMatrix + randmat;
% colormapindividualMatrix(colormapindividualMatrix>1) = 1;
% colormapindividualMatrix(colormapindividualMatrix<0) = 0;

%GROUP BASED ON
[groupList,groupArrayList] = listSpitter(exportSubStruct,groupBasedOn);

%assign a color for grouped median plots
groupColorMap = determineCMAP(cmapchoice,groupList);
groupColorMatrix = zeros(length(groupArrayList),size(groupColorMap,2));
for i=1:length(groupArrayList)
    cA = groupArrayList{i};
    idx = strcmp(groupList,cA);
    groupColorMatrix(i,:) = groupColorMap(idx,:); %this is used to color each individual trace according to the colorBasedOn criteria selected
end

%establish subplot groups
[subplotList,subplotArrayList] = listSpitter(exportSubStruct,subplotGroupBasedOn);

%establish the string array to use for labeling traces
[labelList,labelArrayList] = listSpitter(exportSubStruct,labelBasedOn);




%% how to make correlations

corrMat.Smad.none = cellMat.Smad.all;
corrMat.Reporter.none = cellMat.Reporter.all;
smadtracesNC = cellMat.Smad.NC;
trueCorr = struct();
%1. abundance, foldchange, difference
%2. rate, integral
fnames = fieldnames(corrMat);
for fnum = 1:length(fnames)
    fnstr = fnames{fnum};
    abundance = corrMat.(fnstr).none;
    
    valout = fcanddiff(abundance,stimulationFrame,basalLength,timeMatrix);
    
    corrMat.(fnstr).abundance = valout.abundance;
    corrMat.(fnstr).foldchange = valout.foldchange;
    corrMat.(fnstr).difference = valout.difference;
    corrMat.(fnstr).nuccytoXabundance = (valout.abundance).*smadtracesNC;
    corrMat.(fnstr).zerotoone = valout.zerotoone;
    
    valout = fcanddiff(corrMat.(fnstr).nuccytoXabundance,stimulationFrame,basalLength,timeMatrix);
    
    corrMat.(fnstr).nuccytoXabundanceFC = valout.foldchange;
    corrMat.(fnstr).nuccytoXfoldchange = (valout.abundance).*smadtracesNC;
    corrMat.(fnstr).nuccyto = smadtracesNC;
    
    subfnames = fieldnames(corrMat.(fnstr));
    %abundance, foldchange, difference
    for subfnum = 1:length(subfnames)
        %scalar, rate, integral
        subfnstr = subfnames{subfnum};
        vals = corrMat.(fnstr).(subfnstr);
        [scalar,rate,integral] = rateintegralfunc(vals,stimulationFrame,timeMatrix);
        trueCorr.(fnstr).([subfnstr 'Xscalar']) = scalar;
        trueCorr.(fnstr).([subfnstr 'Xrate']) = rate;
        trueCorr.(fnstr).([subfnstr 'Xintegral']) = integral;
        
    end
end




%% find peaks
sort1 = 'Smad';
sort2 = smadnormalizestr;
sort3 = ReporterTotalOrMedian;
alltracesP = corrMat.(sort1).(sort2);
reptracesP = trueCorr.(sort1).([sort2 'Xscalar']);
timeMatrixP = timeMatrix;

[doseList,doseListArray] = listSpitter(exportSubStruct,'dosestr');

divnum = 3;

savedir = [mfiledir '/' fileDateStr ' - plots'];
savenamestr = ['SORT ' sort2 ' of' ' ' sort1 ' TEST ' 'smad-' smadnormalizestr ' rep-' reporternormalizestr '' sort3 ' ' num2str(divnum) 'classes'];


[sceneList,sceneListArray] = listSpitter(exportSubStruct,'scene');

cellNumMat = horzcat(exportSubStruct.CellNum);


cellNumMat = cellNumMat';
cellNumVec = cellNumMat(:,stimulationFrame);
[sceneList,sceneArrayList] = listSpitter(exportSubStruct,'scene');
stimFrameCellID = cell(size(cellNumVec));
for i = 1:length(cellNumVec)
   cellnumstr = num2str(cellNumVec(i)); 
   scenestr = sceneArrayList{i};
   stimFrameCellID{i} = [scenestr '-' cellnumstr];
end

if strcmp(fileDateStr,'2017_10_30 plate exp2')
    plotTheseCells = '(s01-27|s01-28|s01-29|s01-34|s01-37|s01-41)';
elseif strcmp(fileDateStr,'2017_10_30 plate exp3')
    plotTheseCells = '(s04-14|s04-17|s04-19|s04-22|s04-23|s04-29|s04_31)';
    plotTheseCells = 's04';
end

[~,~,~,d] = regexp(stimFrameCellID,plotTheseCells);    
didx = ~cellfun(@isempty,d,'UniformOutput',1);

    
    
    
    tmat = timeMatrixP(didx,:);
    toteOrMed = {SmadTotalOrMedian,SmadTotalOrMedian};
    traceArray = {'Smad','Smad'};
    tracenormArray = {'abundance','foldchange'};
    for tracenum = 1:length(traceArray)
        
        tracestr = traceArray{tracenum};
        toteOrMedstr = toteOrMed{tracenum};
        tracenormstr = tracenormArray{tracenum};
        alltraces = corrMat.(tracestr).(tracenormstr);
        
        
        divtracesP = cellMat.(tracestr).division;
        tracesforSortP = cellMat.('Smad').all;
        
        
        celltracemat = alltraces(didx,:);
        divtracemat = divtracesP(didx,:);
        tracesforSort = tracesforSortP(didx,:);
        
%         [top1,bot1,~] = booterfunc(alltraces);
%         ylimits = [min(bot1(:)) max(top1(:))];
        ylimits1 = [min(alltraces(:)) prctile(alltraces(:),99.9)];
        ylimits2 = [min(alltraces(:)) prctile(alltraces(:),99.9)];
        xlimits = [min(tmat(:)) max(tmat(:))];
        
        cmap = parula(sum(didx));
        cmap = feval('lines',sum(didx));
        
        D1 = ones(length(cmap),1);
        D2 = 3;
        cmaparray = mat2cell(cmap,D1,D2);
        
        f22 = figure(22);
        f22.Position = [702 520 300*length(doseList) 420];
        subplot(length(traceArray),1,(tracenum-1)+ 1)
            idx = true(size(celltracemat,1),1);
            cellmat = celltracemat(idx,:);
            timemat = tmat(idx,:);
            p = plot(timemat',cellmat'); hold on
            [p.LineWidth] = deal(0.5);
            ylim(ylimits1)
            xlim(xlimits)
            xlabel('minutes')
            ylabel([toteOrMedstr ' ' tracestr ' '  tracenormstr])
            
            [p.Color] = deal(cmaparray{:});
        
        
        f25 = figure(25);
        f25.Position = [702 100 300*length(doseList) 420];
        subplot(length(traceArray),1,(tracenum-1)+ 1)
            cellmat = celltracemat(idx,:);
            timemat = tmat(idx,:);
            
            [top1,bot1,med1] = booterfunc(cellmat);
            nval1 = num2str(length(idx));
            nnan = ~isnan(med1);
            tm = timemat(1,nnan);
            mm1 = med1(nnan);
            mt1 = top1(nnan);
            mb1 = bot1(nnan);
            timevec = horzcat(tm,fliplr(tm)); 
            cellfillvec = horzcat(mm1+mt1,fliplr(mm1-mb1));
            if ~isempty(cellfillvec)
                p1 = plot(tm,mm1);hold on
                p1.Color = cmap(tracenum,:);
                p1.LineWidth = 2;
                p1.DisplayName = ['class-' num2str(i) ' n=' nval1];
                p1.DisplayName = ['n=' nval1];
                fi = fill(timevec,cellfillvec,'g');
                fi.FaceAlpha = 0.2;
                fi.FaceColor = cmap(tracenum,:);
                fi.EdgeColor = fi.FaceColor;
            end
            ylim(ylimits2)
            xlim(xlimits)
            xlabel('minutes')
            ylabel([toteOrMedstr ' ' tracestr ' '  tracenormstr])

        if strcmp(tracestr,'Smad')
            o = findobj(gca,'Type','Line');
            l = legend(gca,o);
            l.FontSize = 7;
        end
        
        
        
    end
    
    f = findobj('Type','Figure');
[f.Color] = deal('w');
h = findobj(f,'Type','Axes');
[h.Color] =  deal([0.95 0.95 0.95]);
[h.GridColor] =  deal([1 1 1]);
[h.XGrid] =  deal('on');
[h.XGrid] =  deal('on');
[h.YGrid] =  deal('on');
[h.GridAlpha] =  deal(1);
[h.LineWidth] = deal(1);
[h.XLim] = deal([-120 360]);
o = findobj(h,'Type','Line');
[o.LineWidth] = deal(2);
f22.Position = [1335 195 258 444];
f25.Position = [500 195 258 444];
f22.InvertHardcopy = 'off';
f25.InvertHardcopy = 'off';
    
    
    
    

olddir = pwd;
cd(mfiledir)
if ~isdir(savedir)
    mkdir(savedir)
end
cd(savedir)
% if ~isdir(savenamestr)
%     mkdir(savenamestr)
% end
% cd(savenamestr)
saveas(f22,[savenamestr 'errorbar.fig'],'fig');
saveas(f22,[savenamestr 'errobar.png'],'png');
saveas(f25,[savenamestr ' traces.fig'],'fig');
saveas(f25,[savenamestr ' traces.png'],'png');
copyfile([mdir '.m'],[savedir '/' savenamestr '.m'])
cd(olddir)


stophere=1;
end




function cmaplz = colorcubemodified(cmaplength,cmapstr)
%default darkmin = 0.2
%default brightmax = 2
darkmin = 0.5;
colorationmin = 0.4;
brightmax = 2;

%generate colormap based on number of cells tracked
cnew=[];
%                 ccc = vertcat(colormap('summer'),colormap('autumn'),colormap('winter'),colormap('spring'));
%                 ccc = vertcat(colormap('hsv'),colormap('hot'));
ccc = flipud(colormap(cmapstr));
cmapccc = zeros(1000,3);
for j = 1:size(cmapccc,2)
    x = linspace(0,1,size(ccc,1));
    v = ccc(:,j);
    xq = linspace(0,1,1000);
    cmapccc(1:length(xq),j) = interp1(x,v,xq);
end

%                 ccc = colormap(jet(1000));
cccyc = 0;
for k = 1:size(cmapccc,1)
    cvec = cmapccc(k,:);
    if sum(cvec)>darkmin && sum(abs(diff(cvec)))>colorationmin && sum(cvec)<brightmax
        cccyc = cccyc+1;
        cnew(cccyc,:) = cvec;
    end
end

ccnew = zeros(cmaplength,3);
for j = 1:size(ccnew,2)
    x = linspace(0,1,size(cnew,1));
    v = cnew(:,j);
    xq = linspace(0,1,cmaplength);
    ccnew(1:length(xq),j) = interp1(x,v,xq);
end

cmaplz = ccnew;
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

function cT = extractTraces(exportStruct,indices,speciesstr,finalFrame,stimulationFrame,basalLength,medfiltnum)
cT = struct();
xtstrArray = {'nuccyto','median','mean','total'};
for xtnum = 1:length(xtstrArray)
    %define how (median, mean, total) and what to quantify (EGFP, RFP)
    xtstr = xtstrArray{xtnum};
    xTracesString = [xtstr speciesstr];
    
    tMatrix = vertcat(exportStruct(indices).timeMatrix);
    tvec = tMatrix(1,:);
    stimulatedFrame = find(tvec>60,1,'first');
    
    % extract the cell traces for the desired number of frames
    cellTracesFull = vertcat(exportStruct(indices).(xTracesString));
    cellTracesUnfilt = cellTracesFull(:,1:finalFrame); %88x50 [needs to be 50x88]
    % cellTraces = medfilt1(cellTracesUnfilt,medfiltnum,[],2,'omitnan','truncate'); %use median filtering
    % cellTraces = smooth(cellTracesUnfilt,0.1,'sgolay'); %use median filtering
    cellTraces = movmean(cellTracesUnfilt,medfiltnum,2,'Endpoints','shrink');
    cellTraces = cellTracesUnfilt;
    
    % normalize by basal values
    if (stimulationFrame-basalLength)<1
        basalLength = stimulationFrame-1;
    end
    %determine cell traces abundance (normalize to basal population median)
    basalValues = cellTraces(:,stimulationFrame-basalLength:stimulationFrame);
    cellTracesAbundance = cellTraces./nanmedian(basalValues(:));
    
    %determine cell traces difference (subtract basal values of individual cells)
    basalValues = cellTracesAbundance(:,stimulationFrame-basalLength:stimulationFrame);
    basalVector = nanmedian(basalValues,2);
    invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
    invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
    basalMatrix = ones(size(cellTraces,2),1)*basalVector';
    cellTracesDifference = cellTracesAbundance - basalMatrix';
    
    %determine cell traces fold-change
    cellTracesFoldChange =cellTracesAbundance.*(invBasalMatrix');
    
    %determine zero to one
    basalValues = cellTraces(:,stimulationFrame-basalLength:stimulationFrame);
    vmin = nanmean(basalValues,2);
    
    stimulatedValues = cellTraces(:,stimulationFrame:stimulatedFrame);
    [vmax,vmaxidx] = max(stimulatedValues,[],2);
    stimulatedValues = nan(3,length(vmaxidx));
    for i = 1:length(vmaxidx)
        vidx = vmaxidx(i)+stimulationFrame-1;
        stimulatedValues(1:3,i) = [cellTraces(i,vidx-1) cellTraces(i,vidx) cellTraces(i,vidx+1)];
    end
    vmax = nanmean(stimulatedValues,1)';
    
    
    scaleFactor = 1./(vmax - vmin);
    scaleMatrix = ones(size(cellTraces,2),1)*scaleFactor';
    scaleFactorMin = vmin.*scaleFactor;
    scaleMatrixMin = ones(size(cellTraces,2),1)*scaleFactorMin';
    cellTracesOne = cellTraces.*(scaleMatrix');
    cellTracesZero = cellTracesOne - scaleMatrixMin';
    
    %assign cellTraces to structure for easy access to data.
    cT.(xtstr).none = cellTraces;
    cT.(xtstr).abundance = cellTracesAbundance;
    cT.(xtstr).difference = cellTracesDifference;
    cT.(xtstr).foldchange = cellTracesFoldChange;
    cT.(xtstr).zerotoone = cellTracesZero;
end
end


function allidx = jittercorrfunc(pmat,jittercut)
%jitter correction in reporter
jittercorrectMat = pmat;
diffrmatz = diff(jittercorrectMat,1,2);
diffrmat = [diffrmatz(:,1) diffrmatz];
upidx = diffrmat>jittercut;
downidx = diffrmat<(-1.*jittercut);
allidx = upidx | downidx;
end

function cmap = determineCMAP(colormapChoice,uniqueColoring)
cmapBig = feval(colormapChoice,length(uniqueColoring));
cmapidx = floor((linspace(1,size(cmapBig,1),length(uniqueColoring))));
cmap = cmapBig(cmapidx,:);
end

function [valList,valArrayList] = listSpitter(exportSubStruct,valBasedOn)
indices = true(1,length(exportSubStruct));
if ~iscell(valBasedOn)
    valArrayList = vertcat({exportSubStruct(indices).(valBasedOn)});
    valList = unique(valArrayList);
else
    for i = 1:length(indices)
        concatstr = [];
        for j = 1:length(valBasedOn)
            cstr = valBasedOn{j};
            newstr = exportSubStruct(i).(cstr);
            concatstr = [concatstr ' - ' newstr];
        end
        valArrayList{i} = concatstr;
    end
    valList = unique(valArrayList);
end
end

function [tracemat,divtracemat,allmat,nanidx] = divisionFunction(cTracks,smadtraces,stimulationFrame,timeMatrix,normalizestr)

if strcmpi(normalizestr,'total')
    appendlog = true;
else
    appendlog = false;
end


celldivtraces = cTracks.Smad.total.foldchange;
repdivtraces = cTracks.Reporter.total.foldchange;

celldivtracesMed = cTracks.Smad.median.foldchange;
repdivtracesMed = cTracks.Reporter.median.foldchange;


nanidx = ~isnan(celldivtraces(:,stimulationFrame));

celltraces = smadtraces(nanidx,:);
tmat = timeMatrix(nanidx,:);
rmat = repdivtraces(nanidx,:);
pmat = celldivtraces(nanidx,:);
rmatmed = repdivtracesMed(nanidx,:);
pmatmed = celldivtracesMed(nanidx,:);
divtracemat = nan(size(rmat));
tracemat = nan(size(rmat));
allmat = nan(size(rmat));

for i = 1:size(rmat,1)
    tvec = round(tmat(i,:));
    pvec = pmat(i,:);
    rvec = rmat(i,:);
    cellvec = celltraces(i,:);
    cellnan = isnan(cellvec);
    ps = smooth(pvec,0.05,'sgolay');
    rs = smooth(rvec,0.05,'sgolay');
    pvec = ps;
    rvec = rs;
    pvec(cellnan) = NaN;
    rvec(cellnan) = NaN;
    
    rdiff = zeros(1,size(rmat,2));
    pdiff = zeros(1,size(rmat,2));
    for j = 1:size(rmat,2)
        afound = find(tvec > tvec(j)-20,1,'first');
        bfound = find(tvec < tvec(j)+20,1,'last');
        a = max([1 afound]);
        b = min([bfound size(rmat,2)]);
        rdiff(1,j) = rvec(b)./rvec(a);
        pdiff(1,j) = pvec(b)./pvec(a);
    end
    pdiff(cellnan) = NaN;
    rdiff(cellnan) = NaN;
    
   
    
    se = strel('disk',2);
    pdiffless = imdilate(pdiff<0.8,se);
    rdiffless = imdilate(rdiff<0.8,se);
    divpoint = pdiffless & rdiffless;
    CC = bwconncomp(divpoint);
    PX = CC.PixelIdxList;
    imglog = false(size(rvec));
    for k = 1:length(PX)
        px = PX{k};
        pxlog = false(size(rvec));
        pxlog(px) = true;
        se = strel('disk',9);
        pxdilate = find(imdilate(pxlog,se));
        rtest = rdiff(pxdilate);
        idx = rtest>1.1 | rtest<0.9;
        pxnew = pxdilate(idx);
        imglog(pxnew) = true;
    end
    se =strel('disk',1);
    imglogc = imclose(imglog,se);
    CC = bwconncomp(imglogc);
    PX = CC.PixelIdxList;
    pxlog = cellfun(@(x) length(x)>1,PX,'UniformOutput',1);
    PXnew = PX(pxlog);
    CC.PixelIdxList = PXnew;
    CC.NumObjects = length(PXnew);
    
    
    %now interpolate through the dividing region
    cellvec = celltraces(i,:);
    cellvectest = cellvec;
    cellvectest(isnan(pdiff))=NaN;
    tracevec = nan(size(rvec));
    divvec = nan(size(rvec));
    pxstartmat = [1 cellfun(@(x) min([x(end)+1 length(rvec)]),PXnew,'UniformOutput',1)];
    pxendmat = [cellfun(@(x) max([1 x(1)-1]),PXnew,'UniformOutput',1) length(rvec)];
    for k = 1:length(pxstartmat)
        pxstart = pxstartmat(k);
        pxend = pxendmat(k);
        if k == 1
            adjustval =0;
        else
            if appendlog
                adjustval = tracevec(pxendmat(k-1))-cellvec(pxstart);
            else
                adjustval = 0;
            end
        end
        tracevec(pxstart:pxend) = cellvec(pxstart:pxend)+adjustval;
        if k>1
            c = min([max([2 pxendmat(k-1)]) length(rvec)-1]);
            d = min([max([2 pxstartmat(k)]) length(rvec)-1]);
            x = [c-1 c d d+1];
            v = [cellvectest(c-1) cellvectest(c) cellvectest(d)+adjustval cellvectest(d+1)+adjustval];
            xq = c:1:d;
            interpp = interp1(x,v,xq);
            divvec(c:d) = interpp;
        end
    end
    allvec = tracevec;
    allvec(~isnan(divvec)) = divvec(~isnan(divvec));
    tracemat(i,:) = tracevec;
    divtracemat(i,:) = divvec;
    allmat(i,:) = allvec;
    
%     figure(2)
%     subplot(1,3,1);
%     plot(tvec,cellvec,'LineWidth',1);hold on
%     plot(tvec,tracevec,'LineWidth',2);hold on
%     plot(tvec,divvec,'k:','LineWidth',2);hold off
%     xlim([-100 1400])
%     subplot(1,3,2);
%     plot(tvec,pvec,'LineWidth',1);hold on
%     plot(tvec,rvec,'LineWidth',1);hold off
%     xlim([-100 1400])
%     ylim([0 4])
%     subplot(1,3,3);
%     plot(tvec,pdiff,'LineWidth',1);hold on
%     plot(tvec,rdiff,'LineWidth',1);hold off
%     xlim([-100 1400])
%     ylim([0 2])
%     if i>130
%         stophere=1;
%     end
end
end

function [smoothtraces,smoothdivtraces,smoothtracesall] = smoothforplot(celltraces,divtraces,tmat)

forsort1 = celltraces;
% forsort1(~isnan(divtraces)) = divtraces(~isnan(divtraces));
tracesforSort = forsort1;

sm = nan(size(tracesforSort));
for i = 1:size(sm,1)
    tvec = tmat(i,:);
    fsvec = tracesforSort(i,:);
    smvec = smooth(tvec,fsvec,0.1,'loess');
    sm(i,:) = smvec;
end
sm(isnan(tracesforSort))=NaN;

smoothdivtraces = nan(size(divtraces));
smoothdivtraces(~isnan(divtraces)) = sm(~isnan(divtraces));

smoothtraces = sm;
smoothtracesall = smoothtraces;
smoothtraces(~isnan(divtraces)) = NaN;

% figure(9)
% subplot(1,2,1);
% plot(tracesforSort');ylim([0 4])
% subplot(1,2,2);
% plot(sm'); ylim([0 4])
% sss=1;
end

function [sortedCellMat,sortidx] = sorttraces(cellMat,sortstr,tmat,sortTime)
% celltraces = cellMat.(sortstr).traces;
% divtraces = cellMat.(sortstr).division;
% forsort1 = celltraces;
% forsort1(~isnan(divtraces)) = divtraces(~isnan(divtraces));
alltraces = cellMat.(sortstr).all;
tracesforSort = alltraces;
tvec = tmat(1,:);


%determine sortidx
sortFrame = find(tvec(1,:)>sortTime);
if isempty(sortFrame)
    disp('sort frame error')
    sortFrame = stimulationFrame+1;
end
spvecForSorting = tracesforSort(:,sortFrame(1));
[~,sortidx] = sort(spvecForSorting);

sortedCellMat = struct();
fnames = fieldnames(cellMat);
for fnum = 1:length(fnames)
    fnstr = fnames{fnum};
    subMat = cellMat.(fnstr);
    subfnames = fieldnames(subMat);
    for subfnum = 1:length(subfnames)
        subfnstr = subfnames{subfnum};
        subtraces = subMat.(subfnstr);
        
        sortedsubcelltraces = zeros(size(tracesforSort));
        for k = 1:length(sortidx)
            sortedsubcelltraces(k,:) = subtraces(sortidx(k),:);
        end
        subtraces = sortedsubcelltraces;
        subMat.(subfnstr) = subtraces;
    end
    sortedCellMat.(fnstr) = subMat;
    
end
end


function newCellMat = updatealltraces(cellMat,timenan)

newCellMat = cellMat;
fnames = fieldnames(cellMat);
for fnum = 1:length(fnames)
    fnstr = fnames{fnum};
    cellfn = cellMat.(fnstr);
    
    fnames2 = fieldnames(cellfn);
    for fnum2 = 1:length(fnames2)
        fnstr2 = fnames2{fnum2};
        
        cellfn2 = cellfn.(fnstr2);
        cellnew = cellfn2(timenan,:);
        newCellMat.(fnstr).(fnstr2) = cellnew;
    end
end
end

%% new functions
function cT = fcanddiff(abundance,stimulationFrame,basalLength,timeMatrix)
cT = struct();
cellTraces = abundance;

% normalize by basal values
if (stimulationFrame-basalLength)<1
    basalLength = stimulationFrame-1;
end
%determine cell traces abundance (normalize to basal population median)
basalValues = cellTraces(:,stimulationFrame-basalLength:stimulationFrame);
cellTracesAbundance = cellTraces./nanmedian(basalValues(:));

%determine cell traces difference (subtract basal values of individual cells)
basalValues = cellTracesAbundance(:,stimulationFrame-basalLength:stimulationFrame);
basalVector = nanmedian(basalValues,2);
invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
basalMatrix = ones(size(cellTraces,2),1)*basalVector';
cellTracesDifference = cellTracesAbundance - basalMatrix';

%determine cell traces fold-change
cellTracesFoldChange =cellTracesAbundance.*(invBasalMatrix');


tvec = timeMatrix(1,:);
stimulatedFrame = find(tvec>60,1,'first');

%determine zero to one
basalValues = cellTraces(:,stimulationFrame-basalLength:stimulationFrame);
vmin = nanmean(basalValues,2);

stimulatedValues = cellTraces(:,stimulationFrame:stimulatedFrame);
[vmax,vmaxidx] = max(stimulatedValues,[],2);
stimulatedValues = nan(3,length(vmaxidx));
for i = 1:length(vmaxidx)
    vidx = vmaxidx(i)+stimulationFrame-1;
    stimulatedValues(1:3,i) = [cellTraces(i,vidx-1) cellTraces(i,vidx) cellTraces(i,vidx+1)];
end
vmax = nanmean(stimulatedValues,1)';


scaleFactor = 1./(vmax - vmin);
scaleMatrix = ones(size(cellTraces,2),1)*scaleFactor';
scaleFactorMin = vmin.*scaleFactor;
scaleMatrixMin = ones(size(cellTraces,2),1)*scaleFactorMin';
cellTracesOne = cellTraces.*(scaleMatrix');
cellTracesZero = cellTracesOne - scaleMatrixMin';

%assign cellTraces to structure for easy access to data.
cT.none = cellTraces;
cT.abundance = cellTracesAbundance;
cT.difference = cellTracesDifference;
cT.foldchange = cellTracesFoldChange;
cT.zerotoone = cellTracesZero;
end

function [scalar,rate,integral] = rateintegralfunc(vals,stimulationFrame,timeMatrix)
%scalar
scalar = vals;

%rate
[rate,~] = gradient(vals,timeMatrix(1,:),timeMatrix(:,1));

%integral
preintegral = cumsum(vals,2,'omitnan');
preintegral(isnan(vals)) = NaN;
basalValues = preintegral(:,stimulationFrame);
basalVector = nanmedian(basalValues,2);
basalMatrix = ones(size(preintegral,2),1)*basalVector';
integral = preintegral - basalMatrix';
end

function [rho,nvalue] = determineCorrelation(pMat,rMat,corrtype,stimulationFrame)
rho = zeros(size(pMat,2),size(rMat,2));
nvalue = zeros(size(pMat,2),size(rMat,2));
for st = stimulationFrame:size(pMat,2) %number of frames
    x = pMat(:,st);
    for rt = st:size(rMat,2)
        y = rMat(:,rt);
        xn = ~isnan(x);
        yn = ~isnan(y);
        cidx = xn&yn;
        xi = x(cidx);
        yi = y(cidx);
        
        %             xtf = ~isoutlierz(xi,'quartiles');
        %             ytf = ~isoutlierz(yi,'quartiles');
        %             oidx = xtf&ytf;
        oidx = true(size(xi));
        xu = xi(oidx);
        yu = yi(oidx);
        
        if isempty(yu) || isempty(xu)
            r =0;
        else
            r = corr(xu,yu,'type',corrtype,'rows','all');
        end
        rho(st,rt) = r;
        nvalue(st,rt) = length(xu);
    end
end
%minimum number of data points for considering correlation
nidx  = nvalue<30;
rho(nidx) = 0;
end

function tf = isoutlierz(xi,funame)
q1 = prctile(xi,25);
q2 = prctile(xi,50);
q3 = prctile(xi,75);
iqr = q3-q1;
outlierdist = iqr*1.5;
lowcut = q1-outlierdist;
highcut = q3+outlierdist;

tf1 = xi>lowcut;
tf2 = xi<highcut;
tf = ~(tf1&tf2);
end

function [top1,bot1,med1] = booterfunc(cellmat1)
nanmat = ~isnan(cellmat1);
nval = sum(nanmat,1);
nvalidx = nval>4;
cellmat1(:,~nvalidx) = NaN;

iter=20;
% botvals = bootstrp(100,@(x) prctile(x,10,1),cellmat1);
% topvals = bootstrp(100,@(x) prctile(x,90,1),cellmat1);
% medvals = bootstrp(100,@(x) prctile(x,50,1),cellmat1);
% 
% med1 = nanmean(medvals,1);
% top1 = nanmean(topvals,1)-med1;
% bot1 = med1 - nanmean(botvals,1);


bt = nan(iter,size(cellmat1,2));
tp = nan(iter,size(cellmat1,2));
md = nan(iter,size(cellmat1,2));
for i = 1:size(cellmat1,2)
    cm = cellmat1(:,i);
    bt(:,i) = bootstrp(iter,@(x) prctile(x,10,1),cm);
    tp(:,i) = bootstrp(iter,@(x) prctile(x,90,1),cm);
    md(:,i) = bootstrp(iter,@(x) prctile(x,50,1),cm);
end


med1 = nanmean(md,1);
top1 = nanmean(tp,1)-med1;
bot1 = med1 - nanmean(bt,1);
end
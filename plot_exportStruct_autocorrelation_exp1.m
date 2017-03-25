function plot_exportStruct_autocorrelation_exp1

close all
%determine the location of the matlab function and establish export
%directory in relation to that filepath
mdir = mfilename('fullpath');
%     [~,b ] = regexp(mdir,'/');
    [~,b] = regexp(mdir,'Tracking\w*/');
        if isempty(b)
%             [~,b] = regexp(mdir,'\');
            [~,b] = regexp(mdir,'Tracking\w*\');
        end
    parentdir = mdir(1:b);
    loaddir = strcat(parentdir,'Export');
    exportdir = strcat(parentdir,'LookingAtData');
cd(loaddir);

[~,b ] = regexp(mdir,'/');
    if isempty(b)
        [~,b] = regexp(mdir,'\');
    end
mfiledir =mdir(1:b(end));

exportdirz = exportdir;
%load the exported tracking structure
% FileName = uigetfile('*export.mat');%choose file to load
FileName = '2014_09_30 plate exp1_tracking_export.mat';
cd(loaddir)
load(FileName)


%load metadata associated with the experiment (requires manual input if
%there is ambiguity
[a,~] = regexp(FileName,'_tracking');
datequery = strcat(FileName(1:a-1),'*metaData.mat');
cd(loaddir)
filelist = dir(datequery);
if length({filelist.name}) ==1
    metaData = load(char(filelist.name));
else
    filename = uigetfile();
    metaData = load(filename);
end

timeVec = metaData.timeVec;

%load information regarding doses and scenes and tgfbeta addition
[a,~] = regexp(FileName,'_tracking');
datequery = strcat(FileName(1:a-1),'*DoseAndScene*');
cd(loaddir)
filelist = dir(datequery);
    if isempty(filelist)
       dosestruct = makeDoseStruct; %run function to make doseStruct 
    else
        dosestructstruct = load(char(filelist.name));
        dosestruct = dosestructstruct.dosestruct;
    end
    
    %exportstruct
    %datastruct
    %dosestruct
    
    
    
    
coloringChoice = 'scene'; %choose which field based upon which each cell trace will get colored 
colormapChoice = 'lines';
darkenFactor = 1.5;
    
%determine the scenes present in the experiment   
scenestr = 'scene';
sceneListArray = vertcat({exportStruct.(scenestr)});
sceneList = unique(sceneListArray);
sceneListArrayTwo = vertcat({dosestruct.(scenestr)});

%combine the exportStruct information with dosesstruct information
for i=1:length(sceneList)
    sceneChoice=sceneList{i};
    indices = strcmp(sceneListArray,sceneChoice);
    indicestwo = strcmp(sceneListArrayTwo,sceneChoice);


    dose = dosestruct(indicestwo).dose;
    frame = dosestruct(indicestwo).tgfFrame;
    
    dosestr = dosestruct(indicestwo).dosestr;
    framestr = dosestruct(indicestwo).tgfFramestr;
    
    
    [exportStruct(indices).dose] = deal(dose);
    [exportStruct(indices).frame] = deal(frame);
    [exportStruct(indices).dosestr] = deal(dosestr);
    [exportStruct(indices).framestr] = deal(framestr);
end

doseListArray = vertcat({exportStruct.dosestr});
doseList = unique(doseListArray);
    

    
%determine details needed for plotting such as when Tgfbeta is added, etc
%medianSmadbkg
stimulationFrame = exportStruct(1).frame;
smadTracesString = 'medianNucEGFP'; %value to plot
% smadTracesString = 'medianSmadbkg';
reporterTracesString = 'medianNucRFP';
numberOfFrames = size(timeVec,2);
finalFrame = numberOfFrames;




%establish the color map for plotting
coloringArray = vertcat({exportStruct.(coloringChoice)});
coloringList = unique(coloringArray);
indices = true(1,length(exportStruct));
    coloringArrayTrunc = vertcat({exportStruct(indices).(coloringChoice)});
    uniqueColoring = unique(coloringArrayTrunc);
    figure(1)
    cmap = colormap(parula(length(coloringList).*2));
    cmap = colormap(colormapChoice)./darkenFactor;
    close 1
    cmap = cmap; %darken the cmap



%assign a color array using the created colormap based on the choices above
colormapMatrix = zeros(length(coloringArrayTrunc),size(cmap,2));
for i=1:length(coloringArrayTrunc)
   cA = coloringArrayTrunc{i};
   idx = strcmp(uniqueColoring,cA);
   colormapMatrix(i,:) = cmap(idx,:);
%    colormapArray{i} = colorNames{idx};
   colormapArray{i} = cmap(idx,:);
end

%need to determine the number of scenes present and choose the time vector
%depending on the scene from which it was imaged
%THIS WORKS FOR NOW BUT NEEDS TO BE CHANGED
numberOfCells = length(indices);
timeMatrix = zeros(numberOfCells,finalFrame);
    coloringArray = vertcat({exportStruct.(coloringChoice)});
    coloringList = unique(coloringArray);
    coloringArrayTrunc = vertcat({exportStruct(indices).(coloringChoice)});
for i=1:numberOfCells
    sceneChoice=exportStruct(i).scene;
    idx = strcmp(sceneListArray,sceneChoice);
    idxtwo = strcmp(sceneListArrayTwo,sceneChoice);
    
   stimulationFrame = dosestruct(idxtwo).tgfFrame;
   timeMatrix(i,:) = timeVec(idxtwo,1:finalFrame)-timeVec(1,stimulationFrame);  
end
% setTequalZeroToStimulation = timeVector(stimulationFrame);
% xtickTimeVector = timeVector - setTequalZeroToStimulation;



indices = true(1,length(exportStruct));
%function to exract the cell traces, normalized and not
[smadCellTracesNorm,smadCellTraces] = extractTraces(exportStruct,indices,smadTracesString,finalFrame,stimulationFrame);
[reporterCellTracesNorm,reporterCellTraces] = extractTraces(exportStruct,indices,reporterTracesString,finalFrame,stimulationFrame);

smadCellTraces = medfilt1(smadCellTraces,3,[],2,'omitnan','truncate');
smadCellTracesNorm = medfilt1(smadCellTracesNorm,3,[],2,'omitnan','truncate');

plottingMat = smadCellTraces;
plottingMatNorm = smadCellTracesNorm;

        f = figure(213);
        
        f.Position = [100 100 400 800];
        f.Color = 'w';
    
        doseList = {'2.4'}; %tgf stimulated cells exp4
    for condidx  = 1:length(doseList)
    
    setvec = -20:4:100;
    axone = subplot(2,1,1);
    axtwo = subplot(2,1,2);
    h(1) = axone;
    h(2) = axtwo;
    
        %general edits to axes
        for i = 1:length(h)
            %specify specific edits to axes
            if i == 1
                ylimmax = 10;
                ylabelstr = 'nuclear Smad3 fluorescence';
                titlestr = 'Abundance of Smad3';
                Position = [0.1900 0.6412 0.7150 0.2838];
 
            elseif i==2
                ylimmax = 6;
                ylabelstr = 'Fold-change in Smad3';
                titlestr = 'Fold-Change of Smad3';
                Position = [0.1900 0.1412 0.7150 0.2838];
            end
            h(i).Position = Position;
            h(i).NextPlot = 'add';
            h(i).TickLength = [0.03 0.03];
            h(i).Box = 'off';
            h(i).TickDir = 'out';
            h(i).XTick = [-20:20:100];
            h(i).LineWidth = 2;
            h(i).XColor = 'k';
            h(i).YColor = 'k';
            h(i).XLim = [-30 100];
            h(i).YLabel.String =ylabelstr;
            h(i).XLabel.String = 'Time (minutes)';
            h(i).FontSize = 14;
            h(i).FontName = 'helvetica';
            h(i).YLim = [0 ylimmax];
            h(i).Title.String = titlestr;
                stem(h(i),0,10,'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',2,'Marker','none');
                t=text(h(i),5,ylimmax*0.9,'+TGFbeta');
                t.FontSize = 12;

        end
            
        for jimmy = 2:length(setvec)       
            drawnow
            
            color1 = [1 0.6 0]./1.1; %cheese
            color2 = [0 1 0]./1.5; %green
            color3 = [1 0 1]./1.5; %purple
            color4 = [0 0.5 1]./1.5; %blue
            colormapMatrix = [color1;color2;color3;color4]; 
            

            plottingMat = smadCellTraces;
            [~,plottingMatforSort,timeMatrixForPlot] = determineSortedpmat(condidx,doseList,stimulationFrame,exportStruct,plottingMat,timeMatrix);        
                sortedpmat = plottingMatforSort./min(plottingMatforSort(:,stimulationFrame));
                tmat = timeMatrixForPlot;
                    p=plot(h(1),tmat(:,1:jimmy)',sortedpmat(:,1:jimmy)','LineWidth',2);hold on
                        set(p, {'color'}, num2cell(colormapMatrix,2));
                        set(p, 'LineWidth',3);
                
                
                

                
        plottingMat = smadCellTracesNorm;
        [~,plottingMatforSort,timeMatrixForPlot] = determineSortedpmat(condidx,doseList,stimulationFrame,exportStruct,plottingMat,timeMatrix);        
            sortedpmat =plottingMatforSort;
            tmat = timeMatrixForPlot;
                    p=plot(h(2),tmat(:,1:jimmy)',sortedpmat(:,1:jimmy)','LineWidth',2);hold on
                        set(p, {'color'}, num2cell(colormapMatrix,2));
                        set(p, 'LineWidth',3);


    olddir = pwd;
    specialdir = 'F:\Frick\moviePlot';
    cd(specialdir)
    tstr = num2str(jimmy);
    tstrz = 't00';
    tstrz(end-length(tstr)+1:end) = tstr;
    filename = strcat('imgplot-',tstrz,'.jpg');
    saveas(gcf,filename,'jpg');
    cd(olddir);






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
end


function [cellTracesNorm,cellTraces] = extractTraces(exportStruct,indices,xTracesString,finalFrame,stimulationFrame)
% extract the cell traces for the desired number of frames
cellTracesFull = vertcat(exportStruct(indices).(xTracesString));
cellTraces = cellTracesFull(:,1:finalFrame); %88x50 [needs to be 50x88]

% normalize by basal values
basalLength = 3;
if (stimulationFrame-basalLength)<1
    basalLength=0;
end

basalVector = nanmedian(cellTraces(:,stimulationFrame-basalLength:stimulationFrame),2);
invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
cellTracesNorm = cellTraces.*(invBasalMatrix');
end


function [timeVec,alphaAuto,alphaAutoErrUp,alphaAutoErrDown,Xi,Ri,Ti,Frankedi] = calcAutoCorrAustin(Ti,Fi,sortFrame)
%determine auto correlation based on method in Sigal et al 2006 paper
%method to calculate auto correlation is from austin et al 2006 Gene
%network shaping of inherent noise spectra (Science). Gives identical
%results to Sigal et al
alphaAutoErrUp=[];
alphaAutoErrDown=[];
numberOfCells = size(Fi,1);
traceLength = size(Fi,2);
[Xi,Ri,Frankedi] = makeRankedMatrix(Fi,sortFrame);



N = traceLength;
M = numberOfCells;
alphaAuto = nan(1,N);
for tau = 1:N
%adjust tau
    j=tau-1;
%determine denominator of autocorrelation
    numerato = nan(1,M);
    denominato = nan(1,M);
    for m = 1:M
        
        %determine numerator;
            numer = nan(1,N-j);
            for n = 1:N-j
                numer(1,n) = Xi(m,n).*Xi(m,n+j);
            end
            numerato(m) = nanmean(numer);
        
        %determine denominator
            denom = nan(1,N);
            for n = 1:N
                denom(1,n) = Xi(m,n).^2;
            end
            denominato(m) = nanmean(denom);
        
    end
    %compute numerator and denominator
    numerator = nanmean(numerato);
    denominator = nanmean(denominato);
    
%determine autocorrelation
    alphaAuto(tau) = numerator./denominator;
end

    timeVec = Ti(1,:);
    
end

function [timeVec,alphaAuto,alphaAutoErrUp,alphaAutoErrDown,Xi,Ri,Ti,Frankedi] = calcAutoCorrAustinTZero(Ti,Fi,sortFrame)
%determine auto correlation based on method in milo et al 2006 paper
alphaAutoErrUp=[];
alphaAutoErrDown=[];
numberOfCells = size(Fi,1);
traceLength = size(Fi,2);
[Xi,Ri,Frankedi] = makeRankedMatrix(Fi,sortFrame);



N = traceLength;
M = numberOfCells;
alphaAuto = nan(1,N);
for tau = 1:N
%adjust tau
    j=tau-1;
%determine denominator of autocorrelation
    numerato = nan(1,M);
    denominato = nan(1,M);
    for m = 1:M
        
        %determine numerator;
            numer = nan(1,N-j);
            for n = 1:1
                numer(1,n) = Xi(m,n).*Xi(m,n+j);
            end
            numerato(m) = nanmean(numer);
        
        %determine denominator
            denom = nan(1,N);
            for n = 1:1
                denom(1,n) = Xi(m,n).^2;
            end
            denominato(m) = nanmean(denom);
        
    end
    %compute numerator and denominator
    numerator = nanmean(numerato);
    denominator = nanmean(denominato);
    
%determine autocorrelation
    alphaAuto(tau) = numerator./denominator;
end

    timeVec = Ti(1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%plotting
    figure(995)
    colormapMatrix = colormap(parula(size(Xi,1)));
    p = plot(Ti'./60,Xi','LineWidth',1.5,'Color',[0.5 0.1 0.1]);hold on
    set(p, {'color'}, num2cell(colormapMatrix,2));
    xlabel('Time (minutes)');
            ylabel('total nuclear fluorescence (au)')
            title('Level of endogenous nuclear NG-Smad3');   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [Xi,Ri,Frankedi] = makeRankedMatrix(Fi,sortFrame)
removeNaN=1;
numberOfCells = size(Fi,1);
traceLength = size(Fi,2);
sortDim = 1; %or 1?

%first sort your poulation based on TgfFrame
    [~,I] = sort(Fi(:,sortFrame),sortDim);

%next build your sorted fluorescence matrix
    Frankedi = zeros(size(Fi));
    for i = 1:numberOfCells
        Frankedi(i,:) = Fi(I(i),:);
    end

%next determine the rank of each trace at a given point
    Ri = zeros(size(Frankedi));
    for j = 1:traceLength
        FrankediFrame = Frankedi(:,j);
        [~,I] = sort(Frankedi(:,j),sortDim); %I(1) is the index of the highest value (not the rank)
        rankI = zeros(size(FrankediFrame));
            rankI(I) = 1:length(I);
            if removeNaN==1
            rankI(I(isnan(FrankediFrame(I)))) = NaN; %remove all NaN values;
            end
            Ri(:,j) = rankI;
    end


%next determine the relative ranks of cells
    Xi = zeros(size(Ri));
    for i = 1:numberOfCells
        for j = 1:traceLength
            Xi(:,j) = Ri(:,j) - nanmean(Ri(:,j));
        end
    end

end


function [sortedpmat,plottingMatforSort,timeMatrixForPlot] = determineSortedpmat(condidx,doseList,stimulationFrame,exportStruct,plottingMat,timeMatrix)
        doseArray = [exportStruct.dose];
        dlog = (doseArray == str2double(doseList{condidx}));
        idx = dlog;
        plottingMatforSort = plottingMat(idx,:);
        plotmatframesorting = plottingMatforSort(:,stimulationFrame);
        [~,indsort] = sort(plotmatframesorting);
        sortedpmat = zeros(size(plottingMatforSort));
        for inin = 1:length(indsort)
            sortedpmat(inin,:) = plottingMatforSort(indsort(inin),:);
        end
        timeMatrixForPlot = timeMatrix(idx,:);
end
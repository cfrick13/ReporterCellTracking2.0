function script_plot_tracked_cells_hypothesis
close all


%set directory to location of code being used (generally external harddrive
%%%%
%determine path to .m file being executed
    mdir = mfilename('fullpath');
        [~,b] = regexp(mdir,'Tracking\w*/');
            if isempty(b)
                [~,b] = regexp(mdir,'Tracking\w*\');
            end
    parentdir = mdir(1:b);
    exportdir = strcat(parentdir,'Export/');

%determine path to gparent folder
    [~,b ] = regexp(parentdir,'/');
        if isempty(b)
            [~,b] = regexp(parentdir,'\');
        end
        gparentdir = parentdir(1:b(end-1));

    cd(parentdir)
    cd ..


%%

% includeArray = {'(clone69|clone80)'};
    includeStr = '2';
%     excludeStr = '(cond|media|->)';
    excludeStr = 'terrible';
%     excludeStr = '(sb|SB|exclude|excluder|wash|0.00|0.02|0.04|0.07|0.06|0.18)';
%     excludeStr = '(sb|SB|exclude|excluder|wash|0.00|0.02|0.04|0.07|0.18|0.54)';
    
%     datestr= '2017_07_08 plate exp1';
%     datestr= '2017_07_07 plate exp1';
%     datestr= '2017_07_01 plate exp2';
%     datestr= '2017_06_26 plate exp2';
    datestr= '2017_04_17 plate exp4';
%     datestr= '2017_07_01 plate exp1';
%     datestr= '2017_07_01 plate exp2';
%     datestr= '2017_06_24 plate exp3';
%     datestr= '2016_10_12 plate exp4';
%     datestr= '2017_10_23 plate exp1';
%     datestr= '2017_10_24 plate exp1';
%     datestr= '2017_10_24 plate exp2';
%     datestr= '2017_10_27 plate exp1';
%     datestr= '2017_10_30 plate exp1';
    datestr= '2017_10_30 plate exp2';
    datestr= '2017_10_30 plate exp3';
%     datestr= '2017_10_25 plate exp1';
%     datestr= '2017_10_25 plate exp2';
%     datestr= '2017_11_01 plate exp1';
%     datestr= '2017_11_04 plate exp1';
%     datestr= '2017_11_04 plate exp2';
%     datestr= '2017_12_09 plate exp1';
    
    FileName = [datestr '_tracking_export.mat'];
    
    %choose from the following (how to quantify nuclei)
    %total, median
    SmadTotalOrMedian = 'median'; % median, mean, total
    ReporterTotalOrMedian = 'total';
    
    %choose from the following
    %doseconddate, conddate, doseAndCondition, conditions, expdate,
    %dosestr, scene, wells, cellid
    groupBasedOn = {'doseAndCondition'};
    labelBasedOn = {'scene','cellID'};
    subplotGroupBasedOn = {'doseAndCondition'};
    colorBasedOn = {'scene','cellID'}; % can be 'individually' or {'cellID','scene'}
    colorBasedOn = {'individually'};
    
    cmapchoice = 'lines';

    celltracestr = 'Reporter'; %Smad, Reporter
    smadnormalizestr = 'foldchange'; %none, abundance, difference, or foldchange, nuccytoXabundance, nuccyto, nuccytoXabundanceFC
    reporternormalizestr = 'difference';
    sortTime = 16; %minutes
    timeCutoff = 100;
    basalLength = 4;
%     determineExportNucleiStructCompiled(FileName,colorBasedOn)
    plot_tracked_cells_hypothesis(FileName,timeCutoff,sortTime,cmapchoice,groupBasedOn,colorBasedOn,subplotGroupBasedOn,labelBasedOn,basalLength,smadnormalizestr,reporternormalizestr,includeStr,excludeStr,SmadTotalOrMedian,ReporterTotalOrMedian,celltracestr)
    stopre=1;
end
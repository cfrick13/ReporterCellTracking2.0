function [TrackedOut,AllImg,segdataout]=DoTrackSegMStack(directory,channels,crop,timecrop)
%
% [Tracked,AllImg] = DoTrackSeg(dirname,channel, crop, timecrop)
%   Track automatically without opening any gui.
%       dirname = name of directory containing all movie files. all files
%           should be the same stage position.
%           e.g.'\Users\ayaron\Movies\2012-01-27\stage1'
%       channel = cell array of strings containing the segmentation color
%           name. e.g. {'YFP'}
%       crop = part of the image to load. [{ystart,yend},{xstart,xend}]
%       timecrop = [startframe, endframe]
%
% DoTrackSeg
% DoTrackSeg(Tracked)
% DoTrackSeg(Tracked,AllImg)
% DoTrackSeg(dirname,channel, crop, timecrop)
%   Open the GUI to edit the tracking

%definitions for bacterial movies
FileType='';%'','bacteria';
BkgType='';%'','template','bacteria';
BkgPrm=[];%the template if it is a template
SegType='';%'','ear','bacteria','bactear';
global segdata
global gcounter
segdata=[];gcounter=0;
segdataout=segdata;
FromTranalyze=-1;

%global variables
Tracked=[];
AllImg=[];
%

%%%uncomment to make background subtraction occur by including external
%%%files
% SegType='';
% BkgType='template';
% BkgPrm={'C:\Users\zeiss\Documents\MATLAB\bkg\','c1','bkg.tif'};

   
%determine matfile directory
    mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
    if isempty(b)
        [~,b] = regexp(mdir,'\');
    end
    parentdir = mdir(1:b(end)); %folder in which matfile exists
    exportdir = strcat(parentdir,'Export/');

    [~,b ] = regexp(parentdir,'/');
    if isempty(b)
        [~,b] = regexp(parentdir,'\');
    end
    gparentdir = parentdir(1:b(end-1)); %folder in which matfile parentdir exists
    
    
%specify important directory names
    mstackName = 'flat mstack';
    trackName = 'tracking files';
    segmentName = 'segment mstack';
  
        
%initialize global variables    

%set colormap

        
% user selects experiment directory to be analyzed

%make tracking file folder if it does not exist

%determine date of experiment

%subdirectories should include
    %> [ flat mstack ]
    %> [ mstack images ]
    %> [ segment mstack ]
    %> [ tracking files ]



ParseArguments(nargin,nargout)

    function ParseArguments(pnumin,pnumout)
        if pnumin >4
            error('Wrong number of arguments')
        end
        if pnumout==3
            RunManuall;%run a different function that i can change manually
        elseif pnumout>0 %automatic Load, Segment, Track, Calculate
            if ~exist('directory','var')
                %open a directory
                directory=uigetdir(pwd,'Select a directory to load files');
            elseif ~ischar(directory) || ~isdir(directory)
                error('First argument is not a valid directory name')
            end
            [~,MovieName]=fileparts(directory);
            if ~exist('channels','var')
                %try to find the possible colors
                olddir=pwd;
                cd(directory);
                f=[dir('*TIF'),dir('*tif')];
                [~,~,~,chlist]=regexp([f.name],'c[0-9]*_');
                chlist=unique(chlist);
                %chlistname=cellfun(@(x) x(3:end), unique(chlist),'unif',0);
                chlistname=chlist;
                chsel=listdlg('ListString', chlistname, 'SelectionMode','single', 'ListSize',[160,80]);
                channels=chlist(chsel);
                cd(olddir);
            elseif ~iscellstr(channels)
                error('Second argument is not a valid cell array of strings')
            end
            if ~exist('crop','var')
                crop=[];
            end
            if ~exist('timecrop','var')
                AutoTrackSeg(directory,channels,crop);
            else
                AutoTrackSeg(directory,channels,crop,timecrop(1),timecrop(2));
            end
            TrackedOut=Tracked;
        else %Open the GUI
            if pnumin==0
                %open a directory
                directory=uigetdir(pwd,'Select a directory to load files');
            end
            if ~exist('crop','var')
                crop=[];
            end
            if ishandle(directory) %called by tranalyze
                FromTranalyze=directory;
                m=guidata(directory);
                directory=m.Tracked;
                if ~isempty(m.AllImg)
                    channels=m.AllImg;
                    pnumin=2;
                end
            end
            if ischar(directory) %should start from files
                if ~ischar(directory) || ~isdir(directory)
                    error('First argument is not a valid directory name or a Tracked structure')
                end
                if ~exist('channels','var')
                    %try to find the possible colors
                    olddir=pwd;
                    cd(directory);
                    
                    
                    % load experiment
                    cd(directory)
                        experimentdir = directory;
                        mstackPath = strcat(experimentdir,'/',mstackName);
                        segmentPath = strcat(experimentdir,'/',segmentName);
                        trackingPath = strcat(experimentdir,'/',trackName);

                    %make tracking file folder if it does not exist
                        cd(experimentdir)

                    %determine date of experiment
                        [a,b] = regexp(experimentdir,'201[0-9]');
                        [~,d] = regexp(experimentdir,'exp[0-9]');
                        ExpDate = experimentdir(a:b+6);expDateStr = experimentdir(a:d); 
                        [a,~] = regexp(ExpDate,'_');ExpDate(a) = '-';
                        
                    %%load helpful metadata
                    cd(exportdir)
                    FileName = expDateStr;
                    datequery = strcat(FileName,'*DoseAndScene*');
                    cd(exportdir)
                    filelist = dir(datequery);
                        if isempty(filelist)
                            error(strcat('need to run ExtractMetadata for-',FileName));
                    %        dosestruct = makeDoseStruct; %run function to make doseStruct 
                        else
                            dosestructstruct = load(char(filelist.name));
                            dosestruct = dosestructstruct.dosestruct;
                        end
                        segInstruct = dosestructstruct.segInstruct;

                    nucleus_seg = segInstruct.nucleus;
                    cell_seg = segInstruct.cell;
                    background_seg = segInstruct.background;
                    channelstoinput = dosestructstruct.channelNames;
                    channelinputs =channelregexpmaker(channelstoinput);
                    bkg = dosestructstruct.BACKGROUND;
                    imgsize = dosestructstruct.dimensions;
                    
                    
                    cd(directory)
%                     dirname = 'flat mstack';
%                     dirlog = isdir(dirname);
%                     if dirlog == false
%                         error('directory does not contain correct file type')
%                     end
%                        
%                     cd(dirname)
                    f=[dir('*mat'),dir('*mat')];
                    [~,~,~,chlist]=regexp([f.name],channelinputs);
                    %chlist=[chlist{:}];
                    chlist=unique(chlist);
                    %chlistname=cellfun(@(x) x(3:end), chlist,'unif',0);
                    chlistname=chlist;
                    chsel=listdlg('ListString', chlistname, 'SelectionMode','single', 'ListSize',[160,80]);
                    channels=chlist(chsel);
                    cd(olddir);
                elseif ~iscellstr(channels)
                    error('Second argument is not a valid cell array of strings')
                end
                %do the segmentation
                if exist('timecrop','var')
                    ParSegment(directory,channels,crop,timecrop(1),timecrop(2));
                else
                    ParSegment(directory,channels,crop);
                end
            elseif iscell(directory) && (pnumin==1 || iscell(channels))%the arguments are the tracked structure and images
                Tracked=directory;
                if exist('channels','var') && ~isempty(channels)
                    AllImg=channels;
                else %open the images
                    AllImg=cell(1,length(Tracked));
                    olddir=pwd;
                    if ~isdir(Tracked{1}.dirname)
                        directory=uigetdir(olddir,'Select a new directory with image files');
                        for i=1:length(Tracked)
                            Tracked{i}.dirname=[directory filesep];
                        end
                    end
                    cd(Tracked{1}.dirname);
                    
                    hh=waitbar(0,'Loading Images...');
                    for cframe=1:length(Tracked)
                        waitbar(cframe/length(Tracked),hh)
                        if isfield(Tracked{cframe},'crop')
                            crop=Tracked{cframe}.crop;
                        else
                            Tracked{cframe}.crop=crop;
                        end
                        %im_bs=loadImage(Tracked{cframe}.filename,crop,Tracked{cframe}.fit,BkgType);
                        im_bs=loadImage(Tracked{cframe}.filename,crop,BkgPrm,BkgType);
                        AllImg{cframe}=im_bs;
                    end
                    delete(hh)
                end
            else
                error('Wrong arguments')
            end
            if ishandle(FromTranalyze)
                DisplayTrack(m.cframe,m.ccell);
                uiwait
                m=[];
                m.Tracked=Tracked;
                m.AllImg=AllImg;
                guidata(FromTranalyze,m)
            else
                DisplayTrack;
                
            end
        end
    end
    function AutoTrackSeg(directory,channels,crop,framestart,framestop)
        olddir=pwd;
        cd(directory);
        directory=[pwd filesep];
        
        imagefiles={};
        for cchannel=channels
            filelist=dir(['*' char(cchannel) '*']);
            switch FileType
                case 'bacteria'
                    fieldpos=strfind(filelist(1).name,'-')+1;
                    fieldpos=fieldpos(end);
                    [~,I]=sort(arrayfun(@(x) str2num(x.name(fieldpos:end-4)), filelist));
                otherwise
                    %directly from metamorph
                    %used to sort by name, now try to sort by date. this
                    %allows for movies that were interupted and then
                    %continued
                    %%[~,~,~,~,fieldpos]=regexp({filelist.name},'_t([^_]+)[\._]');
                    %%fieldpos=[fieldpos{:}];fieldpos=[fieldpos{:}];
                    %%[~,I]=sort(cellfun(@(x) str2num(x),fieldpos));
                    [~,I]=sort(arrayfun(@(x) datenum(x.date), filelist));
            end
            imagefiles(end+1,:)={filelist(I).name};
        end
        
        if ~exist('framestart') || ~exist('framestop') || framestop>length(imagefiles)
            framestart=1;
            framestop=length(imagefiles);
        end

        isWorker=~isempty(getCurrentWorker);
        if ~isWorker
            h=waitbar(0,'Starting AutoTrackSeg...');
        else
            disp(sprintf('%s: Starting AutoTrackSeg...',directory))
        end
        Tracked=cell(1,framestop-framestart+1);
        AllImg=cell(1,framestop-framestart+1);
        for imagenb=framestart:framestop
            Tracked{imagenb}.filename=imagefiles(:,imagenb);
            Tracked{imagenb}.crop=crop;
            Tracked{imagenb}.timecrop=[framestart,framestop];
            %Load and Segment the image
            if ~isWorker
                waitbar((imagenb-1)/length(imagefiles),h,sprintf('Segmenting Image %d',imagenb));
            else
                disp(sprintf('%s: Segmenting Image %d',directory,imagenb))
            end
            [ccell,cimg,cfit]=SegmentFile(imagefiles(:,imagenb),crop,BkgPrm,BkgType,SegType);
            AllImg{imagenb}=cimg;
            Tracked{imagenb}.cells=num2cell(ccell);
            Tracked{imagenb}.filename=imagefiles(:,imagenb);
            Tracked{imagenb}.dirname=directory;
            Tracked{imagenb}.Locked=0;
            Tracked{imagenb}.predictMito=0;
            Tracked{imagenb}.fit=cfit;
            %Track the image
            if ~isWorker
                waitbar((imagenb-2/3)/length(imagefiles),h,sprintf('Tracking Image %d',imagenb));
            else
                disp(sprintf('%s: Tracking Image %d',directory,imagenb))
            end
            if imagenb==framestart
                cells1=[Tracked{framestart}.cells{:}];
                cells1=CalcCellProperties(cells1,AllImg{imagenb});
                Tracked{imagenb}.Locked=0;
                Tracked{imagenb}.cells=num2cell(cells1);
            else
                cells0=[Tracked{imagenb-1}.cells{:}];
                cells1=[Tracked{imagenb}.cells{:}];
                imn0=AllImg{imagenb-1};
                imn1=AllImg{imagenb};
                predictMito=Tracked{imagenb}.predictMito;

                [cells0,cells1,predictMito]=TrackFrame(imn0,cells0,imn1,cells1,predictMito);
                
                Tracked{imagenb-1}.cells=num2cell(cells0);
                Tracked{imagenb-1}.Locked=1;
                Tracked{imagenb}.Locked=0;
                Tracked{imagenb}.predictMito=predictMito;
                Tracked{imagenb}.cells=num2cell(cells1);
            end
            %calulate extra fluorescence channels
            if ~isWorker
                waitbar((imagenb-1/3)/length(imagefiles),h,sprintf('Calculate Fluorescence for Image %d',imagenb));
            else
                disp(sprintf('%s: Calculate Fluorescence for Image %d',directory,imagenb))
            end
            CalcFluorescence(imagenb,imagenb,1)
        end
        if ~isWorker
            delete(h)
        else
            disp(sprintf('%s: Finished Tracking',directory))
        end
        cd(olddir)
    end

    function RunManuall
        
    end

    function ParSegment(directory,channels,crop,framestart,framestop)
        
        olddir=pwd;
        cd(directory);
        directory=[pwd filesep];
        
        imagefiles={};
        for cchannel=channels
            filelist=[dir(['*' char(cchannel) '*.mat']) dir(['*' char(cchannel) '*.MAT'])];
%             filelist=filelist(:);
            filelist=filelist(1);
            [~,ia]=unique({filelist.name});
            filelist=filelist(ia);
            switch FileType
                case 'bacteria'
                    fieldpos=strfind(filelist(1).name,'-')+1;
                    fieldpos=fieldpos(end);
                    [~,I]=sort(arrayfun(@(x) str2num(x.name(fieldpos:end-4)), filelist));
                otherwise %directly from metamorph
                    %used to sort by name, now try to sort by date. this
                    %allows for movies that were interupted and then
                    %continued
                    %%[~,~,~,~,fieldpos]=regexp({filelist.name},'_t([^_]+)[\._]');
                    %%fieldpos=[fieldpos{:}];fieldpos=[fieldpos{:}];
                    %%[~,I]=sort(cellfun(@(x) str2num(x),fieldpos));
                    [~,I]=sort(arrayfun(@(x) datenum(x.date), filelist));
            end
            imagefiles(end+1,:)={filelist(I).name};
        end
        
        
                    cfilename = imagefiles{1};
                    mfile = matfile(char(cfilename));
                    imstack = mfile.flatstack;
        
        
        Tracked=cell(1,length(imagefiles));
        for imagenb=1:size(imstack,3)
            Tracked{imagenb}.filename=imagefiles(:,1);
        end
        
        if ~exist('framestart') || ~exist('framestop') || framestop>size(imstack,3)
            framestart=1;
            framestop=size(imstack,3);
        end
        
        ppool = gcp;
        parallel=ppool.NumWorkers;
        AllImg=cell(1,length(Tracked));
        if parallel==0
            h=waitbar(0,'Segmenting Images...');
        else
            h=[];
            disp('Segmenting Images...')
        end
        
        parfor imagenb=framestart:framestop
%         for imagenb=framestart:framestop
            if parallel==0
                waitbar((imagenb-framestart)/(framestop-framestart),h,sprintf('Segmenting Image %d',imagenb));
            else
                disp(sprintf('Segmenting Image %d',imagenb));
            end
            Tracked{imagenb}.crop=crop;
            Tracked{imagenb}.timecrop=[framestart,framestop];

            [ccell,cimg,cfit]=SegmentFile(imstack(:,:,imagenb),crop,BkgPrm,BkgType,SegType);
            AllImg{imagenb}=cimg;
            Tracked{imagenb}.cells=num2cell(ccell);
            Tracked{imagenb}.filename=cfilename;
            Tracked{imagenb}.dirname=directory;
            Tracked{imagenb}.Locked=0;
            Tracked{imagenb}.predictMito=0;
            Tracked{imagenb}.fit=cfit;

            %display(['Image ' imagefiles{imagenb} ', '  num2str(length(Tracked{imagenb}.cells)) ' cells.']);
        end
        if parallel==0
            delete(h)
        end
        cd(olddir)
    end
    function TrackAll(framestart,framestop)
        
        cells1=[Tracked{framestart}.cells{:}];
        predictMito=Tracked{framestart}.predictMito;
        
        %tic
        h=waitbar(0,'Tracking...');
        for cframe=framestart:framestop
            waitbar((cframe-framestart)/(framestop-framestart),h,['frame ' num2str(cframe) ', cells: ' num2str(length(Tracked{cframe}.cells))]);
            %display(['frame ' num2str(cframe) ', cells: ' num2str(length(Tracked{cframe}.cells))])
            
            cells0=cells1;
            cells1=[Tracked{cframe+1}.cells{:}];
            if ~(Tracked{cframe}.Locked && Tracked{cframe+1}.Locked)
                imn0=AllImg{cframe};
                imn1=AllImg{cframe+1};
                predictMito=Tracked{cframe}.predictMito;
                [cells0,cells1,predictMito]=TrackFrame(imn0,cells0,imn1,cells1,predictMito);
                Tracked{cframe}.cells=num2cell(cells0);
                Tracked{cframe}.Locked=1;
                Tracked{cframe+1}.Locked=0;
                Tracked{cframe+1}.predictMito=predictMito;
            end
            Tracked{cframe+1}.cells=num2cell(cells1);
            if ~ishandle(h)
                return
            end
        end
        delete(h);
        %display(sprintf('finished %.3f', toc))
        
    end
    function [cells0,cells1,predictMito]=TrackFrame(imn0,cells0,imn1,cells1,predictMito)
        recalculate=0;
        [imy,imx]=size(imn1);
        % find registration shift on the mask image
        %pixel [1,1] of the second image is pixel [1,1]+regshift of the
        %first. generally move the first image coords to the second image
        %frame by subtracting regshift
        xcorrim0=min(1,SegRenderNum(cells0,imy,imx));
        xcorrim1=min(1,SegRenderNum(cells1,imy,imx));
        nres=5;
        if max(max(xcorrim0))==0 || max(max(xcorrim1))==0
            %no cells were segmented
            regshift=[0,0];
        else
            xc=xcorr2(xcorrim0(1:nres:end,1:nres:end),xcorrim1(1:nres:end,1:nres:end)); %compute low resolution cross correlation
            [I,J]=find(xc==max(xc(:)),1);
            regshift=nres*[I,J]-[imy,imx]; %determine how far the maxima in the image moves, i think by computing cross correlation      
        end
        % if regshift is small just ignore it
        if max(abs(regshift))<10
            regshift=[0,0];
        end
        
        % Calculate cell properties
        % assume existing: pos, size, mask
        % calculate additional: area, Acom (area-center-of-mass), Ftotal, Fmean, Fmax, Fpixels,  slice
        nbcell0=length(cells0); %number of cells in cells0
        nbcell1=length(cells1);
        cells0=CalcCellProperties(cells0,imn0);
        cells1=CalcCellProperties(cells1,imn1);
        
        % First step: match cells in frame 0 to segments in frame 1.
        % do it by 1. split mitoting cells
        %          2. map cells into segments according to distance.
        %             allow multiple cells to be merged if there is enough energy.
        %          3. check the asignment
        
        %0. if no cells try to look for them
        if nbcell1==0
            for dcell=1:nbcell0 %iterate through the number of cells in image0
                position=min(max(cells0(dcell).pos+cells0(dcell).Acom-[51,51]-regshift,1),[imy,imx]);
                eposition=min(position+100,[imy,imx]);
                imreg=imn1(position(1):eposition(1)-1,position(2):eposition(2)-1);
                if isempty(imreg)
                    continue
                end
                
                ll=localsegment(medfilt2(imreg,[5,5],'symmetric'),SegType);
                [ll2,nn]=bwlabel(ll,8);
                stats = regionprops(ll2, imreg, 'meanIntensity','area','BoundingBox');
                %check that it is similar size and level
                similaritydist=([stats.MeanIntensity]/cells0(dcell).Fmean-1).^2/0.04+([stats.Area]/cells0(dcell).area-1).^2/0.04;
                
                overlapdist=zeros(size(similaritydist));
                origcellreg=SegRenderNum(cells0(dcell),imy,imx);
                origcellreg=origcellreg(position(1):eposition(1)-1,position(2):eposition(2)-1);
                overlapIdx=unique(origcellreg.*ll2);
                
                overlapdist(overlapIdx(2:end))=1;
                
                candict=struct('mask',{},'pos',{},'size',{},'progenitor',{},'descendants',{});
                %for cand=find(similaritydist<4/pi | overlapdist)
                for cand=find(similaritydist<3*4/pi)
                    cellslice=[ceil(stats(cand).BoundingBox(2:-1:1)) ceil(stats(cand).BoundingBox(2:-1:1))+stats(cand).BoundingBox(4:-1:3)-1];
                    cellmask=ll2(cellslice(1):cellslice(3),cellslice(2):cellslice(4))==cand;
                    cellpos=position+cellslice(1:2);
                    candict(end+1)=struct('mask',cellmask,'pos',cellpos,'size',size(cellmask),'progenitor',[],'descendants',[]);
                end
                if size(candict)==0
                    continue
                end
                candict=CalcCellProperties(candict,imn1);
                %should also check that it is not an existing segment (does not overlap with anything)
                overlapcand = any(CalcCellPairProperties(candict,cells1,length(candict),nbcell1)>0,2);
                candict=candict(~overlapcand);
                
                %imshow(SegRender(candict,imx,imx)+2*SegRender(cells1,imx,imx));show()
                if isempty(candict)
                    continue
                end
                % and if more than one candidate, take the closest one
                % maybe later just recompute from scratch
                [~,canddist]=CalcCellPairProperties(candict,cells0(dcell),length(candict),1);
                [~,candselect]=min(canddist);
                if isempty(cells1)
                    cells1=candict(candselect);
                else
                    cells1(end+1)=candict(candselect);
                end
                nbcell1=nbcell1+1;
            end
        end

        if nbcell1==0
            predictMito=[];
            return
        end
        if nbcell0==0
            predictMito=zeros(1,nbcell1);
            return
        end
        
        %1. split mitoting cells. used to predict mitosis. not done
        %anymore.
        % keep original cell for mitoting cells
        OrigCell=[];
        predictMito=[];
        %if any(predictMito)
        %    for j=find(predictMito)
        %       cells0(end+1)=cells0(j); %LB: this is failing sometimes (j>length(cells0))
        %        cells0(end).Fmask=0.5*cells0(end).Fmask;
        %        cells0(j).Fmask=0.5*cells0(j).Fmask;
        %        OrigCell(end+1)=j;
        %    end
        %    nbcell0=length(cells0);
        %    cells0=CalcCellProperties(cells0,imn0);
        %end
        
        %2. map cells
        
        %Calculate properties of cell pairs: overlap, and distance
        [overlap,celldistance]=CalcCellPairProperties(cells0,cells1,nbcell0,nbcell1,regshift);
        %Calculate the cost functions
        overlapnorm0=bsxfun(@rdivide,overlap,[cells0.area]');
        fmapping=bsxfun(@ldivide,[cells0.Ftotal]',0.9*[cells1.Ftotal]);
        cost=overlapnorm0+1./sqrt(512/imx*celldistance+1);
        fused=zeros(nbcell0,nbcell1);
        %cost=bsxfun(@rdivide,overlapnorm0,sum(overlapnorm0))+1./sqrt(512/imx*celldistance+1);
        
        %find the mapping
        cellmapping=zeros(nbcell0,nbcell1);
        cellscore=zeros(nbcell0,nbcell1);
        mappedcells=zeros(1,nbcell0);
        NNmapping=1;
        HAmapping=1;
        
        %if have some links use them
        for ccell=1:nbcell0
            if isempty(cells0(ccell).descendants)
            elseif isscalar(cells0(ccell).descendants)
                cellmapping(ccell,cells0(ccell).descendants)=1;
            else %assume 2 descendants
                cellmapping(ccell,cells0(ccell).descendants(1))=1;
                cellmapping(end+1,:)=0;
                cellmapping(end,cells0(ccell).descendants(2))=1;
                
                cells0(end+1)=cells0(ccell);
                usedpercent=[cells1(cells0(ccell).descendants).Ftotal]./cells0(ccell).Ftotal;
                cells0(ccell).Fmask=usedpercent(1)*cells0(ccell).Fmask;
                cells0(end).Fmask=usedpercent(2)*cells0(end).Fmask;
                cellscore(end+1,:)=0;
                mappedcells(end+1)=0;
                fmapping(end+1,:)=fmapping(ccell,:)/usedpercent(2);
                fmapping(ccell,:)=fmapping(ccell,:)/usedpercent(1);
                cost(end+1,:)=cost(ccell,:);
                OrigCell(end+1)=ccell;
                recalculate=1;
            end
        end
        nbcell0=length(cells0);
        fused=cellmapping.*fmapping;
        cellscore=cellmapping.*cost.*fmapping;
        mappedcells=sum(cellmapping,2);
        
        %complete the mapping
        while sum(NNmapping(:).*HAmapping(:))>0
            fresidue=bsxfun(@times,fmapping,(1-sum(cellmapping./fmapping)));
            fresidue=min(max(fresidue,0.001),1);%cut it at 1 and 0.001
            fcost=cost.*fresidue.*fresidue;
            %prefer orphan segments
            %fcost=bsxfun(@plus,fcost,0.1*(1-sum(cellmapping)));
            fcost=fcost+bsxfun(@times,fcost,0.1*(1-sum(cellmapping)));
            %maybe square it maybe not. for cells that have divided once definitly should be hard to divide again.
            %        fcost=cost*((fmapping*(1-sum(cellmapping/fmapping,axis=0)))**2).clip(max=1,min=0.001)
            fcost(mappedcells==1,:)=0.0001;
            NNmapping=zeros(nbcell0,nbcell1);
            [~,tmp]=max(fcost,[],2);
            NNmapping(sub2ind(size(fcost),1:nbcell0,tmp'))=1;
            %HAmapping=CalcHungarian(1./fcost);
            HAmapping=CalcHungarian(1./max(fcost,0.0001),1/0.01);
            
            cellmapping=cellmapping+NNmapping.*HAmapping;
            cellscore=cellscore+NNmapping.*HAmapping.*fcost;
            mappedcells=sum(cellmapping,2);
            
            %how much of the energy of the original cell is mapped to the
            %new cell
            fused=fused+NNmapping.*HAmapping.*fresidue;
            %if a cell mapped to a low energy segment it probably mitosed so remove .5 of the energy and keep it unmapped
            %TBD: maybe it was mapped to an undersegmented cell. check segmentation.
            for mitcell=find(sum((fused<0.66).*cellmapping,2))'
%                 display('here, please check this code again')
%                continue
                usedpercent=fused(mitcell,cellmapping(mitcell,:)==1);
                cells0(end+1)=cells0(mitcell);
                %maybe 0.9-used and 0.1+used
                cells0(end).Fmask=(1-usedpercent)*cells0(end).Fmask;
                cells0(mitcell).Fmask=(usedpercent)*cells0(mitcell).Fmask;
                nbcell0=nbcell0+1;
                cellmapping(end+1,:)=0;
                cellscore(end+1,:)=0;
                mappedcells(end+1)=0;
                fmapping(end+1,:)=fmapping(mitcell,:)/(1-usedpercent);
                fmapping(mitcell,:)=fmapping(mitcell,:)/usedpercent;
                fused(end+1,:)=0;
                fused(mitcell,:)=fused(mitcell,:)/usedpercent;
                cost(end+1,:)=cost(mitcell,:);
                OrigCell(end+1)=mitcell;
                recalculate=1;
            end
        end
        
        % check if i mapped something into the noise. that is, if a dim object mapped into a bright one.
        % check both Ftotal ratio <5 and Fmean ratio <3
        FtotalR=cellmapping*[cells1.Ftotal]'./[cells0.Ftotal]';
        FmeanR=cellmapping*[cells1.Fmean]'./[cells0.Fmean]';
        tmpcellmapping=cellmapping;
        tmpcellmapping((FtotalR>5) & (FmeanR>3),:)=0;
        areaR=cellmapping*([cells0.area]*tmpcellmapping./[cells1.area])';
        % TBD: it is really a noise if its level is only a fraction of the pixels it cover
        cellmapping((FtotalR>5) & (FmeanR>3) & (areaR>0.85),:)=0;
        
        % if it disappears try to look for it
        for dcell=find(sum(cellmapping,2)==0)'
            position=min(max(cells0(dcell).pos+cells0(dcell).Acom-[51,51]-regshift,1),[imy,imx]);
            eposition=min(position+100,[imy,imx]);
            imreg=imn1(position(1):eposition(1)-1,position(2):eposition(2)-1);
            if isempty(imreg)
                continue
            end
            
            ll=localsegment(medfilt2(imreg,[5,5],'symmetric'),SegType);
            [ll2,nn]=bwlabel(ll,8);
            stats = regionprops(ll2, imreg, 'meanIntensity','area','BoundingBox');
            %check that it is similar size and level
            similaritydist=([stats.MeanIntensity]/cells0(dcell).Fmean-1).^2/0.04+([stats.Area]/cells0(dcell).area-1).^2/0.04;
            
            overlapdist=zeros(size(similaritydist));
            origcellreg=SegRenderNum(cells0(dcell),imy,imx);
            origcellreg=origcellreg(position(1):eposition(1)-1,position(2):eposition(2)-1);
            overlapIdx=unique(origcellreg.*ll2);
            
            overlapdist(overlapIdx(2:end))=1;
            
            if ~isempty(cells1)
                emptycell=struct(cells1(end));
                for cfield=fieldnames(emptycell)'
                    emptycell.(cfield{1})=[];
                end
            else
                emptycell=struct('mask',[],'pos',[],'size',[],'progenitor',[],'descendants',[]);
            end
            candict=emptycell;
            
            %for cand=find(similaritydist<4/pi | overlapdist)
            for cand=find(similaritydist<3*4/pi)
                cellslice=[ceil(stats(cand).BoundingBox(2:-1:1)) ceil(stats(cand).BoundingBox(2:-1:1))+stats(cand).BoundingBox(4:-1:3)-1];
                cellmask=ll2(cellslice(1):cellslice(3),cellslice(2):cellslice(4))==cand;
                cellpos=position+cellslice(1:2);
                candict(end+1)=emptycell;
                candict(end).mask=cellmask;
                candict(end).pos=cellpos;
                candict(end).size=size(cellmask);
            end
            candict(1)=[];
            if isempty(candict)
                continue
            end
            candict=CalcCellProperties(candict,imn1);
            %should also check that it is not an existing segment (does not overlap with anything)
            overlapcand = any(CalcCellPairProperties(candict,cells1,length(candict),nbcell1)>0,2);
            candict=candict(~overlapcand);
            
            %imshow(SegRender(candict,imx,imx)+2*SegRender(cells1,imx,imx));show()
            if isempty(candict)
                continue
            end
            % and if more than one candidate, take the closest one
            % maybe later just recompute from scratch
            [~,canddist]=CalcCellPairProperties(candict,cells0(dcell),length(candict),1);
            [~,candselect]=min(canddist);
            cells1(end+1)=candict(candselect);
            nbcell1=nbcell1+1;
            cellmapping(:,end+1)=0;
            cellmapping(dcell,end)=1;
            recalculate=1;
            
        end
        
        % if two predicted mitotic cells are mapped to the same segment, don't mitose
        addedCells=length(OrigCell);
        for mitocellnum=1:addedCells
            mitocell=OrigCell(addedCells+1-mitocellnum);
            mitocellsister=nbcell0-mitocellnum+1;
            if all(cellmapping(mitocell,:)==cellmapping(mitocellsister,:))
                OrigCell(addedCells+1-mitocellnum)=[];
                cells0(mitocell).Fmask=cells0(mitocell).Fmask+cells0(mitocellsister).Fmask;
                cells0(mitocellsister)=[];
                %renumber the OrigCell numbers since we removed a cell
                OrigCell(OrigCell==mitocellsister)=mitocell;
                OrigCell(OrigCell>mitocellsister)=OrigCell(OrigCell>mitocellsister)-1;
                cellmapping(mitocellsister,:)=[];
                recalculate=1;
            end
        end
        if recalculate==1
            nbcell0=length(cells0);
            cells0=CalcCellProperties(cells0,imn0);
            [overlap,celldistance]=CalcCellPairProperties(cells0,cells1,nbcell0,nbcell1,regshift);
            overlapnorm0=bsxfun(@rdivide,overlap,[cells0.area]');
            recalculate=0;
        end
        
        % 3. split cells
        
        % if two cells map to the same segment, split it
        % three options. touching, partiall overlapp or one withing the other.
        SplitOrigCell=1:nbcell1;
        AllOrigSeed=[];
        for j=find(sum(cellmapping)>1)
            seeds={};
            abscell=find(cellmapping(:,j));
            shifts=bsxfun(@plus,cell2mat({cells0(abscell).pos}'),-cells1(j).pos+[50,50]-regshift);
            change=[0,0;0,-1;1,0;-1,0;0,1];
            Target=zeros([100,100]+size(cells1(j).Fpixels));
            Target=AddCell(Target,[0,0],cells1(j).Fpixels,[50,50]);
            OMtmp=zeros(size(Target));
            OM=zeros(size(Target));
            %start with overlapping cells. maybe take only >0.2 overlapnorm0
            for cell=find(overlap(abscell,j))'
                OM=AddCell(OM,[0,0],cells0(abscell(cell)).Fpixels,shifts(cell,:));
            end
            
            %move them around until they fit best
            origshifts=zeros(size(shifts));
            while any(any(origshifts~=shifts))
                origshifts=shifts;
                for cell=find(overlap(abscell,j))'
                    score=[];
                    OM=AddCell(OM,[0,0],-cells0(abscell(cell)).Fpixels,shifts(cell,:));
                    for shft=change'
                        OM=AddCell(OM,[0,0],cells0(abscell(cell)).Fpixels,shifts(cell,:)+shft');
                        
                        score(end+1)=sum(sum(abs((OM-Target)./(OM+Target+1))));
                        OM=AddCell(OM,[0,0],-cells0(abscell(cell)).Fpixels,shifts(cell,:)+shft');
                    end
                    shifts(cell,:)=shifts(cell,:)+change(find(score==min(score),1),:);
                    OM=AddCell(OM,[0,0],cells0(abscell(cell)).Fpixels,shifts(cell,:));
                end
            end
            %add them to a seed list
            OrigSeed=[];
            for cell=find(overlap(abscell,j))'
                seeds{end+1}=AddCell(zeros(cells1(j).size),[0,0],cells0(abscell(cell)).mask,shifts(cell,:)-[50,50]);
                seeds{end}=seeds{end}.*cells1(j).mask;
                OrigSeed(end+1)=abscell(cell);
            end
            % now add the non overlapping cells
            for cell=find(overlap(abscell,j)==0)'
                % find the minima of missing energy
                [I,J]=find((imerode(OM-Target,cells0(abscell(cell)).mask)==OM-Target).*(Target>0));
                candidatePos=[I,J];
                %if it is empty take the minimal one.
                if isempty(candidatePos)
                    [~,tmp]=min(OM(:)-Target(:));
                    [I,J]=ind2sub(size(Target),tmp);
                    candidatePos=[I,J];
                end
                %find the closest one
                [~,candidate]=min(sum(bsxfun(@minus,candidatePos,shifts(cell,:)+cells0(abscell(cell)).Acom).^2,2));
                OMtmp(candidatePos(candidate,1),candidatePos(candidate,2))=1;
                OMtmp=bwlabel((OM-Target)<min(min(OM-Target))/2);
                blobn=OMtmp(candidatePos(candidate,1),candidatePos(candidate,2));
                OMtmp=(OMtmp==blobn);
                
                while sum(sum(OMtmp.*Target)) < cells0(abscell(cell)).Ftotal
                    prevsize=sum(sum(OMtmp));
                    OMtmp=imdilate(OMtmp,strel('diamond',1)).*(Target>0);
                    if prevsize==sum(sum(OMtmp))
                        break
                    end
                end
                OM=AddCell(OM,[0,0],OMtmp.*Target,[0,0]);
                seeds{end+1}=AddCell(zeros(cells1(j).size),[0,0],OMtmp(51:end,51:end),[0,0]).*cells1(j).mask;
                OrigSeed(end+1)=abscell(cell);
            end
            %           take the cells and expand, adding unasociated pixels
            targetarea=cells1(j).area;
            allseeds=sum(cat(3, seeds{:}),3)>0;
            seedsarea=sum(sum( allseeds>0 ));
            while seedsarea<targetarea
                for n=1:length(seeds)
                    seeds{n}=imdilate(seeds{n},ones(3,3)).*cells1(j).mask.*(1-allseeds)+seeds{n};
                end
                allseeds=sum(cat(3, seeds{:}),3)>0;
                if seedsarea==sum(sum(allseeds));
                    disp('may have some broken cell. check err20120604.')
                    break;%no area change, infinite loop;
                end
                seedsarea=sum(sum( allseeds ));
            end
            %           add the new cells to the dictionary
            %           if n cells are overlapping, split the F between all
            fmeanseeds=[];
            fmeanseeds(1,1,:)=[cells0(OrigSeed).Fmean];
            fmeanseeds=sum(bsxfun(@times,cat(3, seeds{:}),fmeanseeds),3);
            numseeds=sum(cat(3, seeds{:}),3);
            for seednum=1:length(seeds)
                cell=seeds(seednum);
                if ~any(any(cell{1}))
                    continue
                end
                tmp=zeros(size(cell{1})+2);
                tmp(2:end-1,2:end-1)=cell{1};
                tmp=imopen(tmp,strel('diamond', 4));
                if any(any(tmp(2:end-1,2:end-1)))
                    cellOpen=tmp(2:end-1,2:end-1);
                else
                    cellOpen=cell{1};
                end
                bbox=regionprops(cellOpen,'BoundingBox');
                cellslice={floor(bbox.BoundingBox(2))+1:floor(bbox.BoundingBox(2)+bbox.BoundingBox(4)),...
                    floor(bbox.BoundingBox(1))+1:floor(bbox.BoundingBox(1)+bbox.BoundingBox(3))};
                cells1(end+1).mask=cellOpen(cellslice{1},cellslice{2});
                cells1(end).pos=cells1(j).pos+floor(bbox.BoundingBox(2:-1:1));
                cells1(end).Fmask=cells1(end).mask.*cells0(OrigSeed(seednum)).Fmean./fmeanseeds(cellslice{1},cellslice{2});
                %NaN could appear where numseeds is zero. also cell mask is
                %zero there and it should just be zero
                cells1(end).Fmask(isnan(cells1(end).Fmask))=0;
                cells1(end).size=size(cells1(end).mask);
                SplitOrigCell(end+1)=j;
                AllOrigSeed(end+1)=OrigSeed(seednum);
            end
        end
        %   remove the old ones, and recalculate everthing
        if any(sum(cellmapping)>1)
            todelete=sum(cellmapping)>1;
            cells1(todelete)=[];
            for i=1:length(AllOrigSeed)
                cellmapping(:,end+1)=0;
                cellmapping(AllOrigSeed(i),end)=1;
            end
            SplitOrigCell(todelete)=[];
            cellmapping(:,todelete)=[];
            
            nbcell1=length(cells1);
            cells1=CalcCellProperties(cells1,imn1);
            [overlap,celldistance]=CalcCellPairProperties(cells0,cells1,nbcell0,nbcell1,regshift);
            overlapnorm0=bsxfun(@rdivide,overlap,[cells0.area]');
            recalculate=0;
        end
        
        
        
        % Build the mapping:
        %   now we should have nbcell0=nbcell1
        %   we might want to allow cells to apear or disappear but we'll see
        %   if i have a ghost from the previous frame i would try to map all the others and then use it. so give it a cost of 0
        %   if nbcell1 is bigger, there are new cells. i will pad with zero rows. these assignments won't give me any gain so i will try to maximize assigmnment to 'real' cells and only the others will be assigned here.
        %   if nbcell0 is bigger, cells turns into ghosts. again, i will pad with 0.
        cost=(overlapnorm0+1./sqrt(512/imx*celldistance+1));
        fmapping=bsxfun(@ldivide,[cells0.Ftotal]',0.9*[cells1.Ftotal]);
        newcellmapping=CalcHungarian(-cost./exp(abs(log(fmapping))),-0.02);
        %%try to rebuild the matrix for the original unsplitted cells and
        %%see it it is the same as before
        %unsplitcellmapping=zeros(nbcell0,max(SplitOrigCell));
        %for j=1:nbcell1
        %    unsplitcellmapping(:,SplitOrigCell(j))=unsplitcellmapping(:,SplitOrigCell(j))+newcellmapping(:,j);
        %end
        %if ~all(all(newcellmapping==cellmapping))
        %    disp('There might be a problem in the tracking')
        %end
        %cellmapping=newcellmapping;
        
        %cellmapping=CalcHungarian(1./(cost.*fmapping.*fmapping));
        % Predict Mito in the next frame add 0.01 so as not to divide by zero
        FmeanRatio=[cells0.Fmean]*cellmapping./[cells1.Fmean];
        FtRatio=[cells0.Ftotal]*cellmapping./[cells1.Ftotal];
        areaRatio=[cells0.area]*cellmapping./[cells1.area];
        % fix cell mapping for the mitotic cells
        % merge back into the original cell
        for mitocell=OrigCell(end:-1:1)
            cellmapping(mitocell,:)=0.5*(cellmapping(mitocell,:)+cellmapping(end,:));
            cellmapping=cellmapping(1:end-1,:);
            cells0(mitocell).Fmask=cells0(mitocell).Fmask+cells0(end).Fmask;
            cells0(end)=[];
            nbcell0=nbcell0-1;
            recalculate=1;
        end
        if recalculate==1
            nbcell0=length(cells0);
            cells0=CalcCellProperties(cells0,imn0);
            recalculate=0;
        end
        predictMito=areaRatio>1.43;
        % if a mito wasn't predicted but a cell have much less size and F, plus there was an appearing cell, there was a mito
        % #    for cell in predictMito.nonzero()[0]:
        for cell=find(FtRatio>1.43)
            origcell= find(cellmapping(:,cell));
            %is there a relatively close appearing cell
            for newcell=find((cost(origcell,:)>0.2) & (sum(cellmapping)==0))
                %associate it as a descendant
                cellmapping(origcell,newcell)=0.5;
                cellmapping(origcell,cell)=0.5;
            end
        end
        predictMito = predictMito & (sum(cellmapping)==1);
        % if a cell mitosed into two daughters that touch
        % remerge the daughters
        %if mitosed into three, unlink random daugter
        cellstoremove=[];
        for cell=1:nbcell0
            descend=find(cellmapping(cell,:));
            if length(descend)>1
                if length(descend)>2
                    cellmapping(cell,descend)=0.5;
                    for cdes=3:length(descend)
                        cellmapping(cell,descend(cdes))=0;                        
                    end
                    descend=find(cellmapping(cell,:));
                end
                descend_overlap=sum(sum(SegRenderCnt(cells1(descend))==2));
                if descend_overlap>0 && 1<0 %never do that
                    %merge the daughters
                    newcell=struct(cells1(descend(1)));
                    for cfield=fieldnames(newcell)'
                        newcell.(cfield{1})=[];
                    end
                    newcell.mask=min(SegRenderNum(cells1(descend)),1);
                    newcell.pos=min(cell2mat({cells1(descend).pos}'));
                    newcell.size=size(newcell.mask);
                    cells1(descend(1))=newcell;
                    cellmapping(:,descend(1))=sum(cellmapping(:,descend),2);
                    cellstoremove=[cellstoremove descend(2:end)];
                    recalculate=1;
                end
            end
        end
        cellmapping(:,cellstoremove)=[];
        cells1(:,cellstoremove)=[];
        if recalculate==1
            nbcell1=length(cells1);
            cells1=CalcCellProperties(cells1,imn1);
            recalculate=0;
        end
        
        %%%     if toplot==1:
        %%%         imshow(SegRenderNum(cells0));colorbar();figure();imshow(SegRenderNum(cells1));colorbar();show()
        % add final information to the Dictionary and save
        for cell=1:nbcell0
            cells0(cell).descendants=find(cellmapping(cell,:));
        end
        for cell=1:nbcell1
            cells1(cell).progenitor=find(cellmapping(:,cell));
        end
    end

    function CalcFluorescence(framestart,framestop,nowaitbar)
        % Calculate Fluorescence from given segmentation data
        % assume four '_' seperated fields
        % first is name, second is color, third is position, last is time
        % the color name is everything but the first two characters
        remdir=pwd; %remember current directory
        cd(Tracked{framestart}.dirname); %change directory to folder containing images
        
        p = gcp('nocreate');
        parallel = p.NumWorkers;
%         parallel=parpool('size');
        %AllImg=cell(1,length(Tracked));
        h=[];
        if (parallel==0 && ~exist('nowaitbar','var'))
            h=waitbar(0,'Calculating Fluorescence...');
        elseif parallel>0
            disp('Calculating Fluorescence...')
        end
        %parfor
        for cframe=framestart:framestop
            if parallel==0
                if ishandle(h)
                    waitbar((cframe-framestart)/(framestop-framestart),h,sprintf('Calc Fluorescence For Image %d',cframe));
                end
            else
                disp(sprintf('Calc Fluorescence for Image %d',cframe));
            end
            cellsDictNew=[Tracked{cframe}.cells{:}];
            dirname=Tracked{cframe}.dirname;
            filename=Tracked{cframe}.filename;
            crop=Tracked{cframe}.crop;
            fileparts=[];
            switch FileType
                case 'bacteria'
                    fileparts=regexp(filename,'(.+-)(.+)(-.+)','tokens');
                otherwise
                    fileparts=regexp(filename,'(.+_s.+)(c[^_])+(.+)','tokens');
            end                    
            %new version
            newfile = fileparts{1,1}{1,1};
            segmentationcolor = newfile{1,2};
            newfile{1,2}='*';
            newfile=[newfile{:}];
            
            tmp=dir(newfile); %this assumes matlab working directory=directory with the files 
            switch FileType
                case 'bacteria'
                    tmp=regexp({tmp.name},'.+-(.+)-.+','tokens');
                otherwise
                    tmp=regexp({tmp.name},'_s.+(c[^_])+','tokens');
            end

            tmp=vertcat(tmp{:});
            colornames=vertcat(tmp{:});
            %ignore DIC during fluorescence calculation
            %tmp=strfind(colornames, 'DIC');
            %nonDIC=cellfun(@isempty, tmp);
            %colornames=colornames(nonDIC,1);
         
            %segmentation color already exists, so don't recalculate
            %tmp=strfind(colornames, segmentationcolor);
            %nonSGM=cellfun(@isempty, tmp);
            %colornames=colornames(nonSGM,1);
            
            nbcell=length(cellsDictNew);
            for imnum=1:length(colornames)
                %generate the image name by changing the second field
                newfile=fileparts{1,1}{1,1};
                newfile{1,2}=colornames{imnum, 1};
                newfile=[newfile{:}];
                %load the image
                if isempty(crop)
                    im=double(imread([dirname,newfile]));
                else
                    im=double(imread([dirname,newfile] ,'PixelRegion', crop));
                end
                imbs=RemoveBackground(im,BkgPrm,BkgType,colornames{imnum,1});
                imn=medfilt2(imbs,[5,5],'symmetric');
                
                for curcell=1:nbcell
                    cellsDictNew(curcell).Fdata.Fpixels{imnum}=double(imn(cellsDictNew(curcell).pos(1):cellsDictNew(curcell).pos(1)+cellsDictNew(curcell).size(1)-1,...
                        cellsDictNew(curcell).pos(2):cellsDictNew(curcell).pos(2)+cellsDictNew(curcell).size(2)-1)).*cellsDictNew(curcell).Fmask;
                    cellsDictNew(curcell).Fdata.Ftotal(imnum)=sum(cellsDictNew(curcell).Fdata.Fpixels{imnum}(:));
                    cellsDictNew(curcell).Fdata.Fmax(imnum)=max(cellsDictNew(curcell).Fdata.Fpixels{imnum}(:));
                    cellsDictNew(curcell).Fdata.Fmean(imnum)=cellsDictNew(curcell).Fdata.Ftotal(imnum)/cellsDictNew(curcell).area;
                    cellsDictNew(curcell).Fdata.Fname=colornames;
                    cellsDictNew(curcell).Fname=segmentationcolor;
                end
            end
     
            Tracked{cframe}.cells=num2cell(cellsDictNew);
        end
        if parallel==0
            if ishandle(h)
                delete(h)
            end
        end

        cd(remdir); %move back to original Matlab directory
    end

    function [overlap,celldistance]=CalcCellPairProperties(cells0,cells1,nbcell0,nbcell1,regshift)
        if ~exist('regshift','var')
            regshift=[0,0];
        end
        % Calculate properties of cell pairs
        % The overlap matrix
        overlap=zeros(nbcell0,nbcell1);
        celldistance=zeros(nbcell0,nbcell1);
        
        for n1=1:length(cells1)
            [a,b]=ind2sub(size(cells1(n1).mask),find(cells1(n1).mask));
            cell1xy=[a(:),b(:)];
            comdist=sqrt(sum(bsxfun(@minus,cells1(n1).pos+cells1(n1).Acom+regshift,cell2mat({cells0.pos}')+cell2mat({cells0.Acom}')).^2,2));
            for n0=1:length(cells0)
                if comdist(n0)>10*max(cells1(n1).size)
                    celldistance(n0,n1)=comdist(n0)^2;
                else
                    celldistance(n0,n1)=min(sum((bsxfun(@plus,cell1xy,cells1(n1).pos+regshift-cells0(n0).pos-cells0(n0).Acom)).^2,2));
                    if celldistance(n0,n1)< (max(cells1(n1).size)+max(cells0(n0).size))^2
                        overlap_slice=[max(cells0(n0).pos,cells1(n1).pos+regshift)+[1,1] min(cells0(n0).pos+cells0(n0).size,cells1(n1).pos+regshift+cells1(n1).size)];
                        overlap(n0,n1)=sum(sum(cells0(n0).mask((overlap_slice(1):overlap_slice(3))-cells0(n0).pos(1),(overlap_slice(2):overlap_slice(4))-cells0(n0).pos(2)).*...
                            cells1(n1).mask((overlap_slice(1):overlap_slice(3))-cells1(n1).pos(1)-regshift(1),(overlap_slice(2):overlap_slice(4))-cells1(n1).pos(2)-regshift(2))));
                    end
                end
                
            end
        end
        
    end
    function [assignment,cost] = munkres(costMat)
        % MUNKRES   Munkres Assign Algorithm
        %
        % [ASSIGN,COST] = munkres(COSTMAT) returns the optimal assignment in ASSIGN
        % with the minimum COST based on the assignment problem represented by the
        % COSTMAT, where the (i,j)th element represents the cost to assign the jth
        % job to the ith worker.
        %
        
        % This is vectorized implementation of the algorithm. It is the fastest
        % among all Matlab implementations of the algorithm.
        
        % Examples
        % Example 1: a 5 x 5 example
        %{
[assignment,cost] = munkres(magic(5));
[assignedrows,dum]=find(assignment);
disp(assignedrows'); % 3 2 1 5 4
disp(cost); %15
        %}
        % Example 2: 400 x 400 random data
        %{
n=400;
A=rand(n);
tic
[a,b]=munkres(A);
toc                 % about 6 seconds
        %}
        
        % Reference:
        % "Munkres' Assignment Algorithm, Modified for Rectangular Matrices",
        % http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
        
        % version 1.0 by Yi Cao at Cranfield University on 17th June 2008
        
        assignment = false(size(costMat));
        cost = 0;
        
        costMat(costMat~=costMat)=Inf;
        validMat = costMat<Inf;
        validCol = any(validMat);
        validRow = any(validMat,2);
        
        nRows = sum(validRow);
        nCols = sum(validCol);
        n = max(nRows,nCols);
        if ~n
            return
        end
        
        dMat = zeros(n);
        dMat(1:nRows,1:nCols) = costMat(validRow,validCol);
        
        %*************************************************
        % Munkres' Assignment Algorithm starts here
        %*************************************************
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   STEP 1: Subtract the row minimum from each row.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dMat = bsxfun(@minus, dMat, min(dMat,[],2));
        
        %**************************************************************************
        %   STEP 2: Find a zero of dMat. If there are no starred zeros in its
        %           column or row start the zero. Repeat for each zero
        %**************************************************************************
        zP = ~dMat;
        starZ = false(n);
        while any(zP(:))
            [r,c]=find(zP,1);
            starZ(r,c)=true;
            zP(r,:)=false;
            zP(:,c)=false;
        end
        
        while 1
            %**************************************************************************
            %   STEP 3: Cover each column with a starred zero. If all the columns are
            %           covered then the matching is maximum
            %**************************************************************************
            primeZ = false(n);
            coverColumn = any(starZ);
            if ~any(~coverColumn)
                break
            end
            coverRow = false(n,1);
            while 1
                %**************************************************************************
                %   STEP 4: Find a noncovered zero and prime it.  If there is no starred
                %           zero in the row containing this primed zero, Go to Step 5.
                %           Otherwise, cover this row and uncover the column containing
                %           the starred zero. Continue in this manner until there are no
                %           uncovered zeros left. Save the smallest uncovered value and
                %           Go to Step 6.
                %**************************************************************************
                zP(:) = false;
                zP(~coverRow,~coverColumn) = ~dMat(~coverRow,~coverColumn);
                Step = 6;
                while any(any(zP(~coverRow,~coverColumn)))
                    [uZr,uZc] = find(zP,1);
                    primeZ(uZr,uZc) = true;
                    stz = starZ(uZr,:);
                    if ~any(stz)
                        Step = 5;
                        break;
                    end
                    coverRow(uZr) = true;
                    coverColumn(stz) = false;
                    zP(uZr,:) = false;
                    zP(~coverRow,stz) = ~dMat(~coverRow,stz);
                end
                if Step == 6
                    % *************************************************************************
                    % STEP 6: Add the minimum uncovered value to every element of each covered
                    %         row, and subtract it from every element of each uncovered column.
                    %         Return to Step 4 without altering any stars, primes, or covered lines.
                    %**************************************************************************
                    M=dMat(~coverRow,~coverColumn);
                    minval=min(min(M));
                    if minval==inf
                        return
                    end
                    dMat(coverRow,coverColumn)=dMat(coverRow,coverColumn)+minval;
                    dMat(~coverRow,~coverColumn)=M-minval;
                else
                    break
                end
            end
            %**************************************************************************
            % STEP 5:
            %  Construct a series of alternating primed and starred zeros as
            %  follows:
            %  Let Z0 represent the uncovered primed zero found in Step 4.
            %  Let Z1 denote the starred zero in the column of Z0 (if any).
            %  Let Z2 denote the primed zero in the row of Z1 (there will always
            %  be one).  Continue until the series terminates at a primed zero
            %  that has no starred zero in its column.  Unstar each starred
            %  zero of the series, star each primed zero of the series, erase
            %  all primes and uncover every line in the matrix.  Return to Step 3.
            %**************************************************************************
            rowZ1 = starZ(:,uZc);
            starZ(uZr,uZc)=true;
            while any(rowZ1)
                starZ(rowZ1,uZc)=false;
                uZc = primeZ(rowZ1,:);
                uZr = rowZ1;
                rowZ1 = starZ(:,uZc);
                starZ(uZr,uZc)=true;
            end
        end
        
        % Cost of assignment
        assignment(validRow,validCol) = starZ(1:nRows,1:nCols);
        cost = sum(costMat(assignment));
    end
    function cellmapping=CalcHungarian(cost,appear_cost)
        if ~exist('appear_cost','var')
            appear_cost=1/0.02;
        end
        [nbcell0,nbcell1]=size(cost);
        if nbcell0>nbcell1
            cost=[cost zeros(nbcell0,nbcell0-nbcell1)]';
        elseif nbcell0<nbcell1
            cost=[cost;zeros(nbcell1-nbcell0,nbcell1)];
        end
        
        cost=[cost appear_cost*ones(size(cost)); appear_cost*ones(size(cost).*[1,2])];
        cellmapping=double(munkres(cost));
        if nbcell0>nbcell1
            cellmapping=cellmapping';
        end
        cellmapping=cellmapping(1:nbcell0,1:nbcell1);
    end
    function newtarget=AddCell(target, targetpos, cellvalue, cellpos)
        relativepos=cellpos-targetpos;
        %find the coordinates of the region inside the target and inside the cellvalue
        trgtXslice=max(1,relativepos(1)+1):min(size(target,1),relativepos(1)+size(cellvalue,1));
        trgtYslice=max(1,relativepos(2)+1):min(size(target,2),relativepos(2)+size(cellvalue,2));
        cellXslice=max(1,-relativepos(1)+1):min(size(cellvalue,1), size(target,1)-relativepos(1));
        cellYslice=max(1,-relativepos(2)+1):min(size(cellvalue,2), size(target,2)-relativepos(2));
        %put the cell values inside target
        newtarget=target;
        newtarget(trgtXslice,trgtYslice)=target(trgtXslice,trgtYslice)+cellvalue(cellXslice,cellYslice);
    end
    function image=SegRenderNum(cells,imy,imx)
        if isempty(cells)
            if ~exist('imy','var') || ~exist('imx','var')
                image=[];
            else
                image=zeros(imy,imx);
            end
            return
        end
        cellpos=cell2mat({cells.pos}');
        cellsize=cell2mat({cells.size}');
        if ~exist('imy','var')
            ystart=min(cellpos(:,1));
            ystop=max(cellpos(:,1)+cellsize(:,1)-1);
        else
            ystart=1;
            ystop=imy;
        end
        if ~exist('imx','var')
            xstart=min(cellpos(:,2));
            xstop=max(cellpos(:,2)+cellsize(:,2)-1);
        else
            xstart=1;
            xstop=imx;
        end
        
        image=zeros(ystop-ystart+1,xstop-xstart+1);
        for j=1:length(cells)
            image=AddCell(image,[ystart,xstart],(j)*cells(j).mask,cells(j).pos);
        end
    end
    function image=SegRenderCnt(cells,imy,imx)
        if isempty(cells)
            image=[];
            return
        end
        cellpos=cell2mat({cells.pos}');
        cellsize=cell2mat({cells.size}');
        if ~exist('imy','var')
            ystart=min(cellpos(:,1));
            ystop=max(cellpos(:,1)+cellsize(:,1)-1);
        else
            ystart=1;
            ystop=imy;
        end
        if ~exist('imx','var')
            xstart=min(cellpos(:,2));
            xstop=max(cellpos(:,2)+cellsize(:,2)-1);
        else
            xstart=1;
            xstop=imx;
        end
        
        image=zeros(ystop-ystart+1,xstop-xstart+1);
        for j=1:length(cells)
            image=AddCell(image,[ystart,xstart],cells(j).mask,cells(j).pos);
        end
    end

    function DisplayTrack(frame,selectedcell)
        
        %  Construct the main GUI figure
        scrsz=get(0,'ScreenSize'); %detect screen size
        guiW=scrsz(3)*0.8; %width of the GUI
        guiH=scrsz(4)*0.8; %height of the GUI
        MainPos=[guiW/10,guiH/10 ,guiW,guiH]; %position of the GUI
        
        fh=figure('Position',MainPos,...
            'MenuBar','none',...
            'Name','Tranalyze- Movie Analysis Tool',...
            'NumberTitle','off');
        set(fh,'KeyPressFcn',@KeyPressCB,...
            'ResizeFcn',@fhResizeFcn) %make the clicking above figure work
        KPF=uicontrol(fh,'KeyPressFcn',@KeyPressCB,'Position',[1,1,1,1]); %workaround to get the focus for KeyPressFcn

        % attach the global data
        guidata(fh,segdata);
        
        % Construct buttons and menus
        %         fhctl=figure('MenuBar','None'); %deleteLB
        butSz=0.1; %size of buttons
        butLS=0.05+0.35-2.15*butSz; %distance between left of figure and menu buttons (as fraction of GUI fig)
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','<<','Position',...
            [butLS,0.8,butSz,butSz],'Callback',@runbackward)
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','<-','Position',...
            [butLS+1.1*butSz,0.8,butSz,butSz],'Callback',{@iPadKeyPress,'leftarrow'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','->','Position',...
            [butLS+2.2*butSz,0.8,butSz,butSz],'Callback',{@iPadKeyPress,'rightarrow'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','>>','Position',...
            [butLS+3.3*butSz,0.8,butSz,butSz],'Callback',@runforward)
        
        butLS=0.8; %redefine distance between left of figure and menu buttons (as fraction of GUI fig)
        butSz=0.05;
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Color','Position',...
            [butLS,0.75-butSz,butSz,butSz],'Callback',{@iPadKeyPress,'c'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Boundary','Position',...
            [butLS+1.1*butSz,0.75-butSz,butSz,butSz],'Callback',{@iPadKeyPress,'b'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Select','Position',...
            [butLS+2.2*butSz,0.75-butSz,butSz,butSz],'Callback',{@iPadKeyPress,'s'})
        
        
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Delete','Position',...
            [butLS,0.8-3.5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'d'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Add','Position',...
            [butLS+1.1*butSz,0.8-3.5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'a'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','^Autoadd','Position',...
            [butLS+2.2*butSz,0.8-3.5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'a','control'})
        
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','segmeNt','Position',...
            [butLS,0.8-5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'n'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Merge','Position',...
            [butLS+1.1*butSz,0.8-5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'m'})
%        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','meRge all','Position',...
%            [butLS+2.2*butSz,0.8-5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'r'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Editseed','Position',...
            [butLS+2.2*butSz,0.8-5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'e'})
        
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','linK','Position',...
            [butLS,0.8-6.5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'k'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Unlink','Position',...
            [butLS+1.1*butSz,0.8-6.5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'u'})
        
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Track','Position',...
            [butLS,0.8-8.5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'t'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','backYrack','Position',...
            [butLS+1.1*butSz,0.8-8.5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'y'})
%          uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Plot','Position',...
%            [butLS+2.2*butSz,0.8-8.5*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'p'})
        
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Goto','Position',...
            [butLS,0.8-10*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'g'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Lock','Position',...
            [butLS+1.1*butSz,0.8-10*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'l'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','Flour','Position',...
            [butLS+2.2*butSz,0.8-10*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'f'})

        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','autoTrack','Position',...
            [butLS+1.1*butSz,0.8-12*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'t','control'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','autosegN','Position',...
            [butLS,0.8-12*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'n','control'})
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','autoFlour','Position',...
            [butLS+2.2*butSz,0.8-12*butSz,butSz,butSz],'Callback',{@iPadKeyPress,'f','control'})
        
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','ipad/main','Position',...
            [butLS,0.8-14*butSz,butSz,butSz],'Callback',@iPadResize)
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','zoom','Position',...
            [butLS+1.1*butSz,0.8-14*butSz,butSz,butSz],'Callback',@iPadZoom)
        uicontrol( 'Parent', fh, 'Units','normalized','Style','pushbutton','string','save','Position',...
            [butLS+2.2*butSz,0.8-14*butSz,butSz,butSz],'Callback',@CloseAll)
        
        
        % define the variables
        ImgSize=max(cell2mat(cellfun(@size,AllImg,'Uniform',0)'));
        imy=ImgSize(1);
        imx=ImgSize(2);
        AR=imx/imy; %movie aspect ratio (keep constant during figure resizing)
        %         set(fh,'MenuBar','None','KeyPressFcn',@KeyPressCB) %deleteLB
        %colormap(gray)
        iPad=0; %bydefault, the program runs on the main monitor, not on the iPad
        if ~exist('frame','var')
            frame=1;
        end
        State='view';%editing state: view, delete, add,
        StateText='';
        ColorCode=1;
        ColorByNum=0;
        Boundary=1;
        if ~exist('selectedcell','var')
            selectedcell=[];
        end
        
        % plot the movie frames
        hPlotAxes = axes(...    % Axes for plotting the selected plot
            'Parent', fh, ...
            'Units', 'normalized', ...
            'Position',[0.05, 0.05, 0.7*[guiH/guiW*AR, 1]/max([guiH/guiW*AR, 1])]);
        axis([0,imx,0,imy])
        plotframe;
        
        function plotframe(overlay)
            set(fh,'CurrentAxes',hPlotAxes ) %move on the plotting axes
            zoomax=axis;
            Nborder=60;
            Nprog=60;
            Ndesc=40;
            Napp=30;
            Ndisapp=50;
            Nsel=45;
            
            IMborder=zeros(imy,imx);
            IMprog=zeros(imy,imx);
            IMdesc=zeros(imy,imx);
            IMapp=zeros(imy,imx);
            IMdisapp=zeros(imy,imx);
            IMselcell=zeros(imy,imx);
            IMborderSel=zeros(imy,imx); 
            
            im=AllImg{frame};
            if isempty(im)
                [im,cfit]=loadImage(Tracked{frame}.filename,Tracked{frame}.crop,BkgPrm,BkgType);
                AllImg{frame}=im;
                Tracked{frame}.cells={};
                Tracked{frame}.Locked=0;
                Tracked{frame}.fit=cfit;
                
            end
            
            for cell=1:length(Tracked{frame}.cells)
                curcell=Tracked{frame}.cells{cell};
                cmask=zeros(double(curcell.size)+[2,2]);
                cmask(2:end-1,2:end-1)=curcell.mask;
                dpos=double(curcell.pos);
                dsize=double(curcell.size);
                border=cmask-imerode(cmask,strel('disk',2));
                IMborder(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)=...
                    IMborder(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)+...
                    border(2:end-1,2:end-1);
                if ~isfield(curcell,'progenitor')
                    curcell.progenitor=[];
                end
                if isempty(curcell.descendants)
                    IMdisapp(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)=...
                        IMdisapp(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)+curcell.mask;
                elseif isempty(curcell.progenitor)
                    IMapp(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)=...
                        IMapp(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)+curcell.mask;
                else
                    prevcell=Tracked{frame-1}.cells{curcell.progenitor};
                    if ~isscalar(curcell.descendants)
                        IMprog(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)=...
                            IMprog(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)+curcell.mask;
                    end
                    if ~isscalar(prevcell.descendants)
                        IMdesc(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)=...
                            IMdesc(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)+curcell.mask;
                    end
                end
                if any(selectedcell==cell)
                    IMselcell(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)=...
                        IMselcell(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)+curcell.mask;
                    IMborderSel(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)=...
                        IMborderSel(dpos(1):dpos(1)+dsize(1)-1,dpos(2):dpos(2)+dsize(2)-1)+...
                        border(2:end-1,2:end-1); 
                end
                
            end
            
            %im=log(double(im))*50/3-57*5/3;
            %im=(3+log10(max(0.001,double(im))));
            
            im(im<0)=0;
            im=im./max(im(:));
            im=log(im);
            
            im=im-0.4*min(im(~isinf(im)));
            im=20*im./max(im(:));
            %im=20*im+100;
            if ColorCode
                im(IMprog>0)=Nprog;
                im(IMdesc>0)=Ndesc;
                im(IMapp>0)=Napp;
                im(IMdisapp>0)=Ndisapp;
                %im(IMselcell>0)=Nsel; %this plots color filling selected cells
            elseif ColorByNum
                cellsnum=(SegRenderNum([Tracked{frame}.cells{:}],imy,imx)*64/length(Tracked{frame}.cells));
                im(cellsnum>0)=cellsnum(cellsnum>0);
            end
            if Boundary
                im(IMborder>0)=Nborder;
            end
            
            
            im(IMborderSel>0)=Nsel; %this plots only color on the boundary of selected cell
            if Tracked{frame}.Locked
                im(1:16,1:16)=100;
            end
            if exist('overlay','var')
                image(im+overlay)
                axis(zoomax)
            else
                image(im)
                axis(zoomax)
            end
            selectedcellstr=num2str(selectedcell,'%d,');
            selectedcellstr=selectedcellstr(1:end-1);
%             set(gcf,'name',['frame:',num2str(frame),' ', Tracked{frame}.filename{1},' #cell: (' ,selectedcellstr, ')/',num2str(length(Tracked{frame}.cells))])
            %new option
            set(gcf,'name',['frame:',num2str(frame),' ', ' #cell: (' ,selectedcellstr, ')/',num2str(length(Tracked{frame}.cells))])
            title(StateText)
            %give the focus to the keypressfcn uicontrol
            uicontrol(KPF)
            
        end
        
        function fhResizeFcn(src,evnt)
            cPos=get(fh,'Position');
            guiW=cPos(3);
            guiH=cPos(4);
            axes_size=0.7*[guiH/guiW*AR, 1]/max([guiH/guiW*AR, 1]);
            
            set(hPlotAxes,'Parent', fh, ...
                'Units', 'normalized', ...
                'Position',[0.05+0.35-axes_size/2, axes_size]);
            axis([0,imx,0,imy])
        end
        function iPadKeyPress(src,evnt,key,mod)
            evnt2.Key=key;
            evnt2.Modifier='';
            if exist('mod','var')
                evnt2.Modifier={mod};
            end
            figure(fh)
            zoom off
            KeyPressCB(src,evnt2);
        end
        function CloseAll(src,evnt)
            export2wsdlg({'Tracked Data','Images'},{'Tracked','AllImg'},{Tracked,AllImg});
            %            close(fhctl);
            %            close(fh);
        end
        
        function iPadResize(src,evnt)
            if iPad==0
                iPad=1;
                set(fh,'Position',MainPos-[1024,0,0,0]);
            else
                iPad=0;
                set(fh,'Position',MainPos);
            end
        end
        function runforward(src,evnt)
            if get(src,'UserData')
                %stop running
                %curzoom=axis(get(fh,'CurrentAxes'));
                %figure(fh)
                %axis(curzoom*5);
                set(src,'UserData',0)
                set(src,'BackgroundColor',[0.7,0.7,0.7])
                %plotframe
            else
                %start running
                %curzoom=axis(get(fh,'CurrentAxes'));
                %figure(fh)
                %axis(curzoom/5);
                set(src,'UserData',1)
                set(src,'BackgroundColor',[1,0,0])
            end
            while get(src,'UserData')
                ccell=[Tracked{frame}.cells{selectedcell}];
                lastframe=(frame==length(Tracked));
                frame=min(frame+1,length(Tracked));
 
                if ~isempty(ccell)
                    if ~lastframe
                        selectedcell=[ccell.descendants];
                        selectedcell=unique(selectedcell);
                    end
                else
                    selectedcell=[];
                end
                plotframe
                drawnow
            end
        end
        function runbackward(src,evnt)
            if get(src,'UserData')
                set(src,'UserData',0)
                set(src,'BackgroundColor',[0.7,0.7,0.7])
            else
                set(src,'UserData',1)
                set(src,'BackgroundColor',[1,0,0])
            end
            while get(src,'UserData')
                ccell=[Tracked{frame}.cells{selectedcell}];
                firstframe=(frame==1);
                frame=max(1,frame-1);
 
                if ~isempty(ccell)
                    if ~firstframe
                        selectedcell=[ccell.progenitor];
                        selectedcell=unique(selectedcell);
                    end
                else
                    selectedcell=[];
                end
                plotframe
                drawnow
            end
        end
        function iPadZoom(src,evnt)
            if strcmp(get(zoom(gcf),'Enable'),'off')
                set(gcf,'WindowButtonDownFcn',[])
                set(gcf,'WindowButtonUpFcn',[])
            end
            State='view';
            StateText='';
            plotframe
            zoom
        end
        
        function KeyPressCB(src,evnt)
            switch evnt.Key
                case 'leftarrow'
                    ccell=[Tracked{frame}.cells{selectedcell}];
                    firstframe=(frame==1);
                    frame=max(1,frame-1);
 
                    if ~isempty(ccell)
                        if ~firstframe
                            selectedcell=[ccell.progenitor];
                            selectedcell=unique(selectedcell);
                        end
                    else
                        selectedcell=[];
                    end
                    plotframe
                case 'rightarrow'
                    ccell=[Tracked{frame}.cells{selectedcell}];
                    lastframe=(frame==length(Tracked));
                    frame=min(frame+1,length(Tracked));
 
                    if ~isempty(ccell)
                        if ~lastframe
                            selectedcell=[ccell.descendants];
                            selectedcell=unique(selectedcell);
                        end
                    else
                        selectedcell=[];
                    end
                    plotframe
                case 'c' %show/hide cell colorcode
                    if strcmp(evnt.Modifier,'control')
                        ColorByNum=~ColorByNum;
                        ColorCode=0;
                    else
                        ColorCode=~ColorCode;
                        ColorByNum=0;
                    end
                    plotframe
                case 'b' %show/hide cell boundary
                    Boundary=~Boundary;
                    plotframe
                case 'l' %toggle lock state
                    Tracked{frame}.Locked=1-Tracked{frame}.Locked;
                    plotframe
                case 'd' %delete segment
                    DeleteSegment;
                    plotframe
                case 'a' %add segment, control=autoadd blob
                    if strcmp(evnt.Modifier,'control')
                        AutoAddSegment;
                    else
                        AddSegment;
                    end
                    title(StateText);
                case 't' %recalulate tracking
                    if strcmp(evnt.Modifier,'control')
                        TrackAll(frame,length(Tracked)-1);
                        plotframe;
                    else
                        RecalcTrack;
                        iPadKeyPress(src,evnt,'rightarrow')
                    end                    
                case 'y' %back tracking
                    BackTrack;
                case 'r' %reset segmentation
                    ResetSeg;
                case 'm' %merge segments
                    MergeSeg;
                    plotframe;
                case 's' %select cell to follow
                    SelectCell;
                    plotframe;
                case 'g' %goto frame
                    fnum=inputdlg('Goto frame:','Goto frame:');
                    if ~isempty(fnum)
                        frame=str2num(char(fnum));
                    end
                    plotframe;
                case 'x' %fix links
                    if strcmp(evnt.Modifier,'control')
                        FixLinks(1);
                    else
                        FixLinks(0);
                    end
                case 'k' %link cells
                    LinkCells;
                    plotframe;
                case 'n' %resegment
                    if strcmp(evnt.Modifier,'control')
                        AutoReSegment;
                    else
                        ReSegment;
                    end
                case 'p' %plot fluorescence of lineage 
                    PlotLineage;
                case 'f' %recalculate fluorescence in the other channels
                    if strcmp(evnt.Modifier,'control')
                        AutoReFlourescence;
                    else
                        ReFlourescence;
                    end
                case 'e' %edit seeds
                    EditSeeds;
                case 'u' %unlink cells (select parent and child to be separated)
                    UnlinkCells;
                     plotframe;
            end
            
            
        end
        function DeleteSegment
            if strcmp(State,'delete')
                set(gcf,'WindowButtonDownFcn',[])
                set(gcf,'WindowButtonUpFcn',[])
                State='view';
                StateText='';
            else
                set(gcf,'WindowButtonDownFcn',@btndown)
                set(gcf,'WindowButtonUpFcn',[])
                State='delete';
                StateText='Select Cell To Delete';
            end
            
            function btndown(src,evnt)
                p=get(gca,'CurrentPoint');
                coord=p(1,[2,1]);
                allpos=cell2mat(cellfun(@(x) double(x.pos),Tracked{frame}.cells,'Uniform',0)');
                allsize=cell2mat(cellfun(@(x) double(x.size),Tracked{frame}.cells,'Uniform',0)');
                isincell=all(bsxfun(@minus,allsize+allpos,coord)>0 & bsxfun(@minus,allpos,coord)<0,2);
                clicked=find(isincell,1);
                if isempty(clicked)
                    return
                end
                pro=Tracked{frame}.cells{clicked}.progenitor;
                if ~isempty(pro)
                    prodes=Tracked{frame-1}.cells{pro}.descendants;
                    %delete this cell from its progenitor
                    if isscalar(prodes)
                        Tracked{frame-1}.cells{pro}.descendants=[];
                    else
                        Tracked{frame-1}.cells{pro}.descendants(prodes==clicked)=[];
                    end
                end
                if frame>1
                    %renumber the other 'cousins' in their parents
                    for ccell=1:length(Tracked{frame-1}.cells)
                        cousins=Tracked{frame-1}.cells{ccell}.descendants;
                        Tracked{frame-1}.cells{ccell}.descendants=cousins-(cousins>clicked);
                    end
                end
                for des=Tracked{frame}.cells{clicked}.descendants
                    Tracked{frame+1}.cells{des}.progenitor=[];
                end
                if frame<length(Tracked)
                    %renumber the other 'cousins' in their children
                    for ccell=1:length(Tracked{frame+1}.cells)
                        cousins=Tracked{frame+1}.cells{ccell}.progenitor;
                        Tracked{frame+1}.cells{ccell}.progenitor=cousins-(cousins>clicked);
                    end
                end
                Tracked{frame}.cells(clicked)=[];
                Tracked{frame}.Locked=0;
                plotframe;
            end
        end
        function AutoAddSegment
            if strcmp(State,'autoadd')
                set(gcf,'WindowButtonDownFcn',[])
                set(gcf,'WindowButtonUpFcn',[])
                State='view';
                StateText='';
            else
                set(gcf,'WindowButtonDownFcn',@btndown)
                set(gcf,'WindowButtonUpFcn',[])
                State='autoadd';
                StateText='Select Cell To AutoAdd';
            end
            
            function btndown(src,evnt)
                p=get(gca,'CurrentPoint');
                coord=p(1,[2,1]);
                gaussian_size=33;
                im=AllImg{frame};
                imseg=localsegment(medfilt2(im,[5,5],'symmetric'),SegType);
                imsegnum=bwlabel(imseg,4);
                blobnum=imsegnum(round(coord(1)),round(coord(2)));
                if blobnum==0
                    return
                end
                newsegmask=(imsegnum==blobnum);
                bbox=regionprops(newsegmask,'BoundingBox');
                cellslice={floor(bbox.BoundingBox(2))+1:floor(bbox.BoundingBox(2)+bbox.BoundingBox(4)),...
                    floor(bbox.BoundingBox(1))+1:floor(bbox.BoundingBox(1)+bbox.BoundingBox(3))};
                if ~isempty(Tracked{frame}.cells)
                    Tracked{frame}.cells{end+1}=struct(Tracked{frame}.cells{end});
                    for cfield=fieldnames(Tracked{frame}.cells{end})'
                        Tracked{frame}.cells{end}.(cfield{1})=[];
                    end
                else
                    Tracked{frame}.cells{1}=struct('mask',[],'pos',[],'size',[],'progenitor',[],'descendants',[]);
                end
                Tracked{frame}.cells{end}.mask=newsegmask(cellslice{1},cellslice{2});
                Tracked{frame}.cells{end}.pos=floor(bbox.BoundingBox(2:-1:1))+1;
                Tracked{frame}.cells{end}.size=size(Tracked{frame}.cells{end}.mask);
                Tracked{frame}.cells{end}=CalcCellProperties(Tracked{frame}.cells{end},AllImg{frame});
                %once i added a cell, maybe should merge overlapping
                %existing cell
                plotframe;
            end
        end
        function AddSegment
            if strcmp(State,'add')
                set(gcf,'WindowButtonDownFcn',[])
                set(gcf,'WindowButtonUpFcn',[])
                State='view';
                StateText='';
            else
                set(gcf,'WindowButtonDownFcn',@btndown)
                set(gcf,'WindowButtonUpFcn',@btnup)
                State='add';
                StateText='Mark Region To Add';
                lh=[];
                newseg=[];
                
            end
            
            function btndown(src,evnt)
                lh=line('Visible','off');
                newseg=[];
                p=get(gca,'CurrentPoint');
                coord=p(1,[2,1]);
                newseg(end+1,:)=coord;
                set(lh,'XData',newseg(:,2),'YData',newseg(:,1),'Marker','none','Color','r','Visible','on');
                set(src,'WindowButtonMotionFcn',@move)
            end
            function move(src,evnt)
                p=get(gca,'CurrentPoint');
                coord=p(1,[2,1]);
                steps=1+ceil(max(abs(newseg(end,:)-coord)));
                newseg=[newseg;[linspace(newseg(end,1),coord(1),steps); linspace(newseg(end,2),coord(2),steps)]'];
                %newseg(end+1,:)=coord;
                set(lh,'XData',newseg(:,2),'YData',newseg(:,1),'Marker','none','Color','r','Visible','on');
            end
            function btnup(src,evnt)
                coord=newseg(1,:);
                steps=1+ceil(max(abs(newseg(end,:)-coord)));
                newseg=[newseg;[linspace(newseg(end,1),coord(1),steps); linspace(newseg(end,2),coord(2),steps)]'];
                set(src,'WindowButtonMotionFcn',[])
                newsegmask=zeros(imy,imx);
                newsegind=sub2ind(size(newsegmask),round(newseg(:,1)),round(newseg(:,2)));
                newsegmask(newsegind)=1;
                newsegmask=imerode(imfill(imdilate(newsegmask,strel('diamond',1))),strel('diamond',1));
                bbox=regionprops(newsegmask,'BoundingBox');
                cellslice={floor(bbox.BoundingBox(2))+1:floor(bbox.BoundingBox(2)+bbox.BoundingBox(4)),...
                    floor(bbox.BoundingBox(1))+1:floor(bbox.BoundingBox(1)+bbox.BoundingBox(3))};
                if ~isempty(Tracked{frame}.cells)
                    Tracked{frame}.cells{end+1}=struct(Tracked{frame}.cells{end});
                    for cfield=fieldnames(Tracked{frame}.cells{end})'
                        Tracked{frame}.cells{end}.(cfield{1})=[];
                    end
                else
                    Tracked{frame}.cells{1}=struct('mask',[],'pos',[],'size',[],'progenitor',[],'descendants',[]);
                end
                Tracked{frame}.cells{end}.mask=newsegmask(cellslice{1},cellslice{2});
                Tracked{frame}.cells{end}.pos=floor(bbox.BoundingBox(2:-1:1))+1;
                Tracked{frame}.cells{end}.size=size(Tracked{frame}.cells{end}.mask);
                Tracked{frame}.cells{end}=CalcCellProperties(Tracked{frame}.cells{end},AllImg{frame});
                plotframe;
            end 
        end
        function RecalcTrack
            imn0=AllImg{frame};
            imn1=AllImg{frame+1};
            predictMito=Tracked{frame}.predictMito;
            cells0=[Tracked{frame}.cells{:}];
            cells1=[Tracked{frame+1}.cells{:}];
            [cells0,cells1,predictMito]=TrackFrame(imn0,cells0,imn1,cells1,predictMito);
            Tracked{frame}.cells=num2cell(cells0);
            Tracked{frame+1}.cells=num2cell(cells1);
            Tracked{frame}.Locked=1;
            Tracked{frame+1}.Locked=0;
            Tracked{frame+1}.predictMito=predictMito;
            plotframe
        end
        function BackTrack
            imn0=AllImg{frame};
            imn1=AllImg{frame+1};
            predictMito=Tracked{frame}.predictMito;
            cells0=[Tracked{frame}.cells{:}];
            cells1=[Tracked{frame+1}.cells{:}];
            %use track frame from 1 to 0. need to flip the desc and pro.
            %flip prog and desc
            prog={cells0.progenitor};
            [cells0.progenitor]=deal(cells0.descendants);
            [cells0.descendants]=deal(prog{:});
            prog={cells1.progenitor};
            [cells1.progenitor]=deal(cells1.descendants);
            [cells1.descendants]=deal(prog{:});
            [cells1,cells0,predictMito]=TrackFrame(imn1,cells1,imn0,cells0,0);
            %flip back
            prog={cells0.progenitor};
            [cells0.progenitor]=deal(cells0.descendants);
            [cells0.descendants]=deal(prog{:});
            prog={cells1.progenitor};
            [cells1.progenitor]=deal(cells1.descendants);
            [cells1.descendants]=deal(prog{:});
            Tracked{frame}.cells=num2cell(cells0);
            Tracked{frame+1}.cells=num2cell(cells1);
            Tracked{frame}.Locked=1;
            Tracked{frame+1}.Locked=0;
            Tracked{frame+1}.predictMito=predictMito;
            plotframe
        end
        function ResetSeg
            newsegmask=SegRenderNum([Tracked{frame}.cells{:}],imy,imx)>0;
            [newsegmaskl,num]=bwlabel((newsegmask)>0);
            bbox=regionprops(newsegmaskl,'BoundingBox');
            cells=struct('mask',{},'pos',{},'size',{},'progenitor',{},'descendants',{});
            for cseg=1:num
                cellslice={floor(bbox(cseg).BoundingBox(2))+1:floor(bbox(cseg).BoundingBox(2)+bbox(cseg).BoundingBox(4)),...
                    floor(bbox(cseg).BoundingBox(1))+1:floor(bbox(cseg).BoundingBox(1)+bbox(cseg).BoundingBox(3))};
                cells(end+1).mask=double(newsegmaskl(cellslice{1},cellslice{2})==cseg);
                cells(end).pos=floor(bbox(cseg).BoundingBox(2:-1:1))+1;
                cells(end).size=size(cells(end).mask);
            end
            Tracked{frame}.cells=num2cell(cells);
            Tracked{frame}.predictMito=0;
            %remove the links from progenitors and descendants
            if frame<length(Tracked)
                for cell=1:length(Tracked{frame+1}.cells)
                    Tracked{frame+1}.cells{cell}.progenitor=[];
                end
            end
            if frame>1
                for cell=1:length(Tracked{frame-1}.cells)
                    Tracked{frame-1}.cells{cell}.descendants=[];
                end
            end
            plotframe
        end
        function ReSegment
            cells=SegmentImage(AllImg{frame},SegType);
            Tracked{frame}.cells=num2cell(cells);
            Tracked{frame}.predictMito=0;
            %remove the links from progenitors and descendants
            if frame<length(Tracked)
                for cell=1:length(Tracked{frame+1}.cells)
                    Tracked{frame+1}.cells{cell}.progenitor=[];
                end
            end
            if frame>1
                for cell=1:length(Tracked{frame-1}.cells)
                    Tracked{frame-1}.cells{cell}.descendants=[];
                end
            end
            plotframe
        end
        function AutoReSegment
            h=waitbar(0,'Segmenting...');
            for cframe=frame:length(Tracked)
                waitbar((cframe-frame)/(length(Tracked)-frame),h,['Segmenting frame ' num2str(cframe)]);
                cells=SegmentImage(AllImg{cframe},SegType);
                Tracked{cframe}.cells=num2cell(cells);
                Tracked{cframe}.predictMito=0;
                if ~ishandle(h)
                    break
                end
            end
            if ishandle(h)
                delete(h)
            end
            %remove the links from progenitors and descendants
            if cframe<length(Tracked)
                for cell=1:length(Tracked{cframe+1}.cells)
                    Tracked{cframe+1}.cells{cell}.progenitor=[];
                end
            end
            if frame>1
                for cell=1:length(Tracked{frame-1}.cells)
                    Tracked{frame-1}.cells{cell}.descendants=[];
                end
            end
            plotframe
        end
        function ReFlourescence
            CalcFluorescence(frame,frame); 
        end
        function AutoReFlourescence
            CalcFluorescence(frame,length(Tracked)); 
        end

        function EditSeeds
            if strcmp(State,'seed')
                set(gcf,'WindowButtonDownFcn',[])
                set(gcf,'WindowButtonUpFcn',[])
                State='view';
                StateText='';
                %recalcseg;
                plotframe;
            else
                set(gcf,'WindowButtonDownFcn',@btndown)
                set(gcf,'WindowButtonUpFcn',[])
                State='seed';
                StateText='Select to add or remove';
                seedmtx=zeros(size(AllImg{frame}));
%                 [Tracked{frame}.cells{:}];
%                 cellsAcom=cell2mat({ans.pos}')+cell2mat({ans.Acom}');
%                seedmtx(sub2ind(size(seedmtx),cellsAcom(:,1),cellsAcom(:,2)))=1;
                maxpos=cellfun(@(x) find(x.Fpixels(:)==max(x.Fpixels(:)),1,'first'),Tracked{frame}.cells);
                [Tracked{frame}.cells{:}];
                hlen=cell2mat({ans.size}');
                cellsmax=  -1 + cell2mat({ans.pos}') + [mod(maxpos',hlen(:,1)),ceil(maxpos'./hlen(:,1))];
                seedmtx(sub2ind(size(seedmtx),cellsmax(:,1),cellsmax(:,2)))=1;
                plotframe(50*imdilate(seedmtx,[1,1,1;1,1,1;1,1,1]))
            end
            
            function recalcseg(src,evnt)
                im=medfilt2(AllImg{frame},[5,5],'symmetric');
                ll=MarkerControlledWatershedSegmentation(im,seedmtx);
                ll=(ll.*SegRenderNum([Tracked{frame}.cells{:}],512,512))>0;
                %[ll2,nn]=bwlabel(ll,4);
                ll2=ll;
                cells=mask2cells(ll2);
                cells=CalcCellProperties(cells,AllImg{frame});
                Tracked{frame}.cells=num2cell(cells);
                Tracked{frame}.predictMito=0;
                %remove the links from progenitors and descendants
                if frame<length(Tracked)
                    for cell=1:length(Tracked{frame+1}.cells)
                        Tracked{frame+1}.cells{cell}.progenitor=[];
                    end
                end
                if frame>1
                    for cell=1:length(Tracked{frame-1}.cells)
                        Tracked{frame-1}.cells{cell}.descendants=[];
                    end
                end
            end
            function btndown(src,evnt)
                p=get(gca,'CurrentPoint');
                coord=p(1,[2,1]);
                indx=sub2ind(size(seedmtx),round(coord(1)),round(coord(2)));
                tmpseedmtx=zeros(size(seedmtx));
                tmpseedmtx(indx)=1;
                nearindx=find(seedmtx & imdilate(tmpseedmtx,[1,1,1;1,1,1;1,1,1]));
                if ~isempty(nearindx)
                    indx=nearindx;
                end
                
                seedmtx(indx)=~seedmtx(indx);
                recalcseg
                plotframe(50*imdilate(seedmtx,[1,1,1;1,1,1;1,1,1]))
            end
        end
        
        function MergeSeg
            if strcmp(State,'merge')
                set(gcf,'WindowButtonDownFcn',[])
                set(gcf,'WindowButtonUpFcn',[])
                State='view';
                StateText='';
            else
                set(gcf,'WindowButtonDownFcn',@btndown)
                set(gcf,'WindowButtonUpFcn',[])
                State='merge';
                StateText='Select Overlapping Cells To Merge';
            end
            
            function btndown(src,evnt)
                p=get(gca,'CurrentPoint');
                coord=p(1,[2,1]);
                allpos=cell2mat(cellfun(@(x) double(x.pos),Tracked{frame}.cells,'Uniform',0)');
                allsize=cell2mat(cellfun(@(x) double(x.size),Tracked{frame}.cells,'Uniform',0)');
                isincell=all(bsxfun(@minus,allsize+allpos,coord)>0 & bsxfun(@minus,allpos,coord)<0,2);
                clicked=find(isincell);
                if sum(isincell)<2
                    return
                end
                newcell=struct(Tracked{frame}.cells{end});
                for cfield=fieldnames(newcell)'
                    newcell.(cfield{1})=[];
                end
                newcell.mask=min(SegRenderNum([Tracked{frame}.cells{clicked}]),1);
                [Tracked{frame}.cells{clicked}];
                newcell.pos=min(cell2mat({ans.pos}'));
                newcell.size=size(newcell.mask);
                %delete tracking information
                %as i renumber stuff, go from end to beginning as to not
                %change stuf i would refer to later
                for cclicked=clicked(end:-1:1)'
                    pro=Tracked{frame}.cells{cclicked}.progenitor;
                    if ~isempty(pro)
                        prodes=Tracked{frame-1}.cells{pro}.descendants;
                        %delete this cell from its progenitor
                        if isscalar(prodes)
                            Tracked{frame-1}.cells{pro}.descendants=[];
                        else
                            Tracked{frame-1}.cells{pro}.descendants(prodes==cclicked)=[];
                        end
                    end
                    if frame>1
                        %renumber the other 'cousins' in their parents
                        for ccell=1:length(Tracked{frame-1}.cells)
                            cousins=Tracked{frame-1}.cells{ccell}.descendants;
                            Tracked{frame-1}.cells{ccell}.descendants=cousins-(cousins>cclicked);
                        end
                    end
                    for des=Tracked{frame}.cells{cclicked}.descendants
                        Tracked{frame+1}.cells{des}.progenitor=[];
                    end
                    if frame<length(Tracked)
                        %renumber the other 'cousins' in their children
                        for ccell=1:length(Tracked{frame+1}.cells)
                            cousins=Tracked{frame+1}.cells{ccell}.progenitor;
                            Tracked{frame+1}.cells{ccell}.progenitor=cousins-(cousins>cclicked);
                        end
                    end
                end
                Tracked{frame}.cells(clicked)=[];
                Tracked{frame}.cells{end+1}=newcell;
                Tracked{frame}.cells{end}=CalcCellProperties(Tracked{frame}.cells{end},AllImg{frame});
                Tracked{frame}.Locked=0;
                
                plotframe;
            end
        end
        function SelectCell
            if strcmp(State,'select')
                set(gcf,'WindowButtonDownFcn',[])
                set(gcf,'WindowButtonUpFcn',[])
                State='view';
                StateText='';
            else
                set(gcf,'WindowButtonDownFcn',@btndown)
                set(gcf,'WindowButtonUpFcn',[])
                State='select';
                StateText='Select Cell To Follow';
            end
            
            function btndown(src,evnt)
                p=get(gca,'CurrentPoint');
                coord=p(1,[2,1]);
                allpos=cell2mat(cellfun(@(x) double(x.pos),Tracked{frame}.cells,'Uniform',0)');
                allsize=cell2mat(cellfun(@(x) double(x.size),Tracked{frame}.cells,'Uniform',0)');
                isincell=all(bsxfun(@minus,allsize+allpos,coord)>0 & bsxfun(@minus,allpos,coord)<0,2);
                clicked=find(isincell,1);
                selectedcell=clicked;
                plotframe;
            end
        end
        
        function PlotLineage

            if isempty(selectedcell)
                StateText='Please Celect Cell First';
                plotframe
            else
                cells0=selectedcell;
%                 Nframes=length(Tracked);
                Nframes=200;
                
                Rt{1}=zeros(Nframes, 1);
                Rm{1}=zeros(Nframes, 1);
                Rt{2}=zeros(Nframes, 1);
                Rm{2}=zeros(Nframes, 1);
                
                cellnum=[];


                fr=1:Nframes;
                newcells=cells0;
                cells=cells0;
                %get progenitor all the way back
                for k=1:frame
                    frame1=frame-k+1; %counting backwards for progenitors
                    for cell1=cells
                        if cell1==-1
                            continue
                        end
                        Rt{1}(frame1,cells==cell1)=Tracked{frame1}.cells{cell1}.Ftotal;
                        Rm{1}(frame1,cells==cell1)=Tracked{frame1}.cells{cell1}.Fmean;
                        if isfield(Tracked{frame1}.cells{cell1}, 'Fdata') %fluorescence data (not used for segmentation)
                            color{1,1}=Tracked{frame1}.cells{cell1}.Fname;
                            for n=2:length(Tracked{frame1}.cells{cell1}.Fdata.Ftotal)+1
                            Rt{n}(frame1,cells==cell1)=Tracked{frame1}.cells{cell1}.Fdata.Ftotal;
                            Rm{n}(frame1,cells==cell1)=Tracked{frame1}.cells{cell1}.Fdata.Fmean;
                            color{1,n}=Tracked{frame1}.cells{cell1}.Fdata.Fname;
                            end
                        else
                            color{1,1}=' tracking color';
                        end
                        cellnum(frame1,cells==cell1)=cell1;
                        if isempty(Tracked{frame1}.cells{cell1}.progenitor)
                            newcells(cells==cell1)=-1;
                        else
                            newcells(cells==cell1)=Tracked{frame1}.cells{cell1}.progenitor;
                        end
                    end
                    cells=newcells;
                    newcells=cells;
                end
                %get descendants
                newcells=cells0;
                cells=cells0;
                for frame1=frame:Nframes
                    for cell1=cells
                        if cell1==-1
                            continue
                        end
                        Rt{1}(frame1,cells==cell1)=Tracked{frame1}.cells{cell1}.Ftotal;
                        Rm{1}(frame1,cells==cell1)=Tracked{frame1}.cells{cell1}.Fmean;
                        cellnum(frame1,cells==cell1)=cell1;
                        if isfield(Tracked{frame1}.cells{cell1}, 'Fdata') %fluorescence data (not used for segmentation)
                            for n=2:length(Tracked{frame1}.cells{cell1}.Fdata.Ftotal)+1
                            Rt{n}(frame1,cells==cell1)=Tracked{frame1}.cells{cell1}.Fdata.Ftotal;
                            Rm{n}(frame1,cells==cell1)=Tracked{frame1}.cells{cell1}.Fdata.Fmean;
                            end
                        end
                        if isempty(Tracked{frame1}.cells{cell1}.descendants)
                            newcells(cells==cell1)=-1;
                        elseif isscalar(Tracked{frame1}.cells{cell1}.descendants)
                            newcells(cells==cell1)=Tracked{frame1}.cells{cell1}.descendants;
                        else
                            newcells(cells==cell1)=Tracked{frame1}.cells{cell1}.descendants(1);
                            newcells(end+1)=Tracked{frame1}.cells{cell1}.descendants(2);
                            Rt{1}(:,end+1)=Rt{1}(:,cells==cell1);
                            Rm{1}(:,end+1)=Rm{1}(:,cells==cell1);
                            if isfield(Tracked{frame1}.cells{cell1}, 'Fdata') %fluorescence data (not used for segmentation)
                                for n=2:length(Tracked{frame1}.cells{cell1}.Fdata.Ftotal)+1
                                Rt{n}(:,end+1)=Rt{n}(:,cells==cell1);
                                Rm{n}(:,end+1)=Rm{n}(:,cells==cell1);
                                end
                            end
                            cellnum(:,end+1)=cellnum(:,cells==cell1);
                        end
                    end
                    cells=newcells;
                    newcells=cells;
                end
               
                fp = figure(); hold on; set(fp, 'Name', ['Cell ' num2str(cellnum(1,1)) ' in frame 1'])
                subplot(2,2,1); plot(1:length(Rt{1}),Rt{1},'.-'); 
                celltags=regexp(num2str(cellnum(end,:)), '\s*', 'split'); xlabel('time (frames)'); 
                ylabel(['total ' color{1,1}(2:end) ' (a.u.)']);
                legend(celltags); legend('off'); %this labels the cells
                subplot(2,2,2); plot(1:length(Rm{1}),Rm{1},'.-'); 
                celltags=regexp(num2str(cellnum(end,:)), '\s*', 'split'); xlabel('time (frames)'); 
                ylabel(['mean ' color{1,1}(2:end) ' (a.u.)']);
                legend(celltags); legend('off'); %this labels the cells
                
                %plot other colors, if fluorescence information exists
                subplot(2,2,3); plot(1:length(Rt{2}),Rt{2},'.-'); 
                celltags=regexp(num2str(cellnum(end,:)), '\s*', 'split');
                xlabel('time (frames)'); ylabel(['total ' color{1,2}{1,1}(2:end) ' (a.u.)']);
                legend(celltags); legend('off'); %this labels the cells
                subplot(2,2,4); plot(1:length(Rm{2}),Rm{2} ,'.-'); 
                celltags=regexp(num2str(cellnum(end,:)), '\s*', 'split'); xlabel('time (frames)'); ylabel(['mean ' color{1,2}{1,1}(2:end) ' (a.u.)']);
                legend(celltags); legend('off'); %this labels the cells   
            end
         end
         
        function XpandSegs
            segsize=0;
            newsegmask=SegRenderNum([Tracked{frame}.cells{:}],imy,imx)>0;
            im=AllImg{frame};
            while segsize<sum(sum(newsegmask))
                segsize=sum(sum(newsegmask));
                tavg=imfilter(newsegmask.*im,fspecial('average',6),'symmetric');
                navg=imfilter(newsegmask,fspecial('average',6),'symmetric');
                bratio=imdilate(newsegmask,strel('diamond',1)).*(1-newsegmask).*im./tavg.*navg;
                newsegmask=newsegmask+(bratio>0.9);
            end
            
            [newsegmaskl,num]=bwlabel((newsegmask)>0);
            bbox=regionprops(newsegmaskl,'BoundingBox');
            cells=struct('mask',{},'pos',{},'size',{},'progenitor',{},'descendants',{});
            for cseg=1:num
                cellslice={floor(bbox(cseg).BoundingBox(2))+1:floor(bbox(cseg).BoundingBox(2)+bbox(cseg).BoundingBox(4)),...
                    floor(bbox(cseg).BoundingBox(1))+1:floor(bbox(cseg).BoundingBox(1)+bbox(cseg).BoundingBox(3))};
                cells(end+1).mask=double(newsegmaskl(cellslice{1},cellslice{2})==cseg);
                cells(end).pos=floor(bbox(cseg).BoundingBox(2:-1:1))+1;
                cells(end).size=size(cells(end).mask);
            end
            Tracked{frame}.cells=num2cell(cells);
            Tracked{frame}.predictMito=0;
            %remove the links from progenitors and descendants
            if frame<length(Tracked)
                for cell=1:length(Tracked{frame+1}.cells)
                    Tracked{frame+1}.cells{cell}.progenitor=[];
                end
            end
            if frame>1
                for cell=1:length(Tracked{frame-1}.cells)
                    Tracked{frame-1}.cells{cell}.descendants=[];
                end
            end
            plotframe
        end
        function FixLinks(dir)
            cframe=frame;
            cellnb0=length(Tracked{cframe}.cells);
            cellnb1=length(Tracked{cframe+1}.cells);
            transition01=zeros(cellnb0,cellnb1);
            transition10=zeros(cellnb0,cellnb1);
            
            for ccell=1:cellnb0
                desc=Tracked{cframe}.cells{ccell}.descendants;
                transition01(ccell,desc)=1;
            end
            for ccell=1:cellnb1
                prev=Tracked{cframe+1}.cells{ccell}.progenitor;
                transition10(prev,ccell)=1;
            end
            if dir==0
                %use frame 0 to correct frame 1
                for ccell=1:cellnb1
                    prev=find(transition01(:,ccell));
                    Tracked{cframe+1}.cells{ccell}.progenitor=prev;
                end
            else
                %use frame 1 to correct frame 0
                for ccell=1:cellnb0
                    desc=find(transition10(ccell,:));
                    Tracked{cframe}.cells{ccell}.descendants=desc;
                end
            end

            plotframe
        end
        function LinkCells
            if strcmp(State,'link')
                set(gcf,'WindowButtonDownFcn',[])
                set(gcf,'WindowButtonUpFcn',[])
                State='view';
                StateText='';
            else
                set(gcf,'WindowButtonDownFcn',@btndown_prog)
                set(gcf,'WindowButtonUpFcn',[])
                State='link';
                StateText='Select Parent';
            end
            clickedP=[];
            function btndown_prog(src,evnt)
                p=get(gca,'CurrentPoint');
                coord=p(1,[2,1]);
                allpos=cell2mat(cellfun(@(x) double(x.pos),Tracked{frame}.cells,'Uniform',0)');
                allsize=cell2mat(cellfun(@(x) double(x.size),Tracked{frame}.cells,'Uniform',0)');
                allcom=allpos+cell2mat(cellfun(@(x) double(x.Acom),Tracked{frame}.cells,'Uniform',0)');
                [~,clickedP]=min(sum((bsxfun(@minus,allcom,coord)).^2,2));
                if isempty(clickedP)
                    return
                end
                %keep track of selected cell before switching frames
                tmp=[Tracked{frame}.cells{selectedcell}];
                if ~isempty(tmp)
                    selectedcell=[tmp.descendants];
                else
                    selectedcell=[];
                end
                %move to the next frame
                frame=frame+1;
                set(gcf,'WindowButtonDownFcn',@btndown_des)
                set(gcf,'WindowButtonUpFcn',[])
                State='link';
                StateText='Select Descendant';
                plotframe;
            end
            function btndown_des(src,evnt)
                p=get(gca,'CurrentPoint'); 
                coord=p(1,[2,1]);
                allpos=cell2mat(cellfun(@(x) double(x.pos),Tracked{frame}.cells,'Uniform',0)');
                allsize=cell2mat(cellfun(@(x) double(x.size),Tracked{frame}.cells,'Uniform',0)');
                allcom=allpos+cell2mat(cellfun(@(x) double(x.Acom),Tracked{frame}.cells,'Uniform',0)');
                [~,clickedD]=min(sum((bsxfun(@minus,allcom,coord)).^2,2));
                if isempty(clickedD)
                    return
                end
                %keep track properly of selected cell
                tmp=[Tracked{frame}.cells{selectedcell}];
                if ~isempty(tmp)
                    at=tmp.progenitor;
                    selectedcell=unique(at);
                else
                    selectedcell=[];
                end
                frame=frame-1;
                set(gcf,'WindowButtonDownFcn',@btndown_prog)
                set(gcf,'WindowButtonUpFcn',[])
                State='link';
                StateText='Select Parent';

                %adding the correct descendant
                Tracked{frame}.cells{clickedP}.descendants(end+1)=clickedD;
                %make sure not to duplicated descendants
                Tracked{frame}.cells{clickedP}.descendants=unique(Tracked{frame}.cells{clickedP}.descendants);
                %adding the correct progenitor
                Tracked{frame+1}.cells{clickedD}.progenitor=clickedP;
                plotframe;
            end 
        end
        function UnlinkCells
            if strcmp(State,'unlink')
                set(gcf,'WindowButtonDownFcn',[])
                set(gcf,'WindowButtonUpFcn',[])
                State='view';
                StateText='';
            else
                set(gcf,'WindowButtonDownFcn',@btndown_prog)
                set(gcf,'WindowButtonUpFcn',[])
                State='unlink';
                StateText='Select Parent';
            end
            clickedP=[];
            function btndown_prog(src,evnt)
                p=get(gca,'CurrentPoint');
                coord=p(1,[2,1]);
                allpos=cell2mat(cellfun(@(x) double(x.pos),Tracked{frame}.cells,'Uniform',0)');
                allsize=cell2mat(cellfun(@(x) double(x.size),Tracked{frame}.cells,'Uniform',0)');
                allcom=allpos+cell2mat(cellfun(@(x) double(x.Acom),Tracked{frame}.cells,'Uniform',0)');
                [~,clickedP]=min(sum((bsxfun(@minus,allcom,coord)).^2,2));
                if isempty(clickedP)
                    return
                end
                
                %deleting current cell from its descendants
                for desc=Tracked{frame}.cells{clickedP}.descendants
                    if Tracked{frame+1}.cells{desc}.progenitor==clickedP;
                        Tracked{frame+1}.cells{desc}.progenitor=[];
                    end
                end
                
                %clear descendants of current cell
                Tracked{frame}.cells{clickedP}.descendants=[];

                plotframe;
            end 
        end
        
    end

end

function cellsDictNew=CalcCellProperties(cellsDict,imn)
% Calculate properties of cells
% assume existing: pos, size, mask
% area, Acom (area-center-of-mass), Ftotal, Fmean, Fmax, Fpixels
cellsDictNew=cellsDict;
nbcell=length(cellsDict);
for cell=1:nbcell
    %if there is an 'Fmask' then it is the percent of F in every pixel that belongs to that cell. if not it is ones.
    if ~isfield(cellsDictNew(cell),'Fmask') || isempty(cellsDictNew(cell).Fmask)
        cellsDictNew(cell).Fmask=double(cellsDictNew(cell).mask);
    end
    cellsDictNew(cell).area=sum(cellsDictNew(cell).mask(:));
    cellsDictNew(cell).Acom=...
        floor([sum(sum(cellsDictNew(cell).mask,2).*(1:cellsDictNew(cell).size(1))')/cellsDictNew(cell).area,...
        sum(sum(cellsDictNew(cell).mask,1).*(1:cellsDictNew(cell).size(2)))/cellsDictNew(cell).area]);
    cellsDictNew(cell).Fpixels=double(imn(cellsDictNew(cell).pos(1):cellsDictNew(cell).pos(1)+cellsDictNew(cell).size(1)-1,...
        cellsDictNew(cell).pos(2):cellsDictNew(cell).pos(2)+cellsDictNew(cell).size(2)-1)).*cellsDictNew(cell).Fmask;
    cellsDictNew(cell).Ftotal=sum(cellsDictNew(cell).Fpixels(:));
    cellsDictNew(cell).Fmax=max(cellsDictNew(cell).Fpixels(:));
    cellsDictNew(cell).Fmean=cellsDictNew(cell).Ftotal/cellsDictNew(cell).area;
    
    %this might cause errors
    cellsDictNew(cell).Fmedian=median(cellsDictNew(cell).Fpixels(:));
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

function [cells,im_bs,cfit]=SegmentFile(img,crop,BkgPrm,BkgType,SegType)
[im_bs,cfit]=loadImage(img,crop,BkgPrm,BkgType);
cells=SegmentImage(im_bs,SegType);
end
function [im_bs,cfit]=loadImage(img,crop,cfit,BkgType)
if ~exist('cfit','var')
    cfit=[];
end
im={};

    if isempty(crop)
        im = img;
    else
        im = img(crop);
    end

im = double(im);
[im_bs,cfit]=RemoveBackground(im,cfit,BkgType);
end
function [imbs,fitresult]=RemoveBackground(im,cfit,BkgType,other)
if ~exist('cfit','var')
    cfit=[];
end
if ~exist('other','var')
    other=[];
end
switch BkgType
    case 'bacteria'
        [imbs,fitresult]=BactRemoveBackground(im,cfit);
    case 'template'
        [imbs,fitresult]=TemplateRemoveBackground(im,cfit,other);
    otherwise
        [imbs,fitresult]=MamRemoveBackground(im,cfit);
end
end
function [imbs,fitresult]=MamRemoveBackground(im,cfit)
[imy,imx]=size(im);
if isempty(cfit)
    %fit background
    sx=1:5:imx;
    sy=1:5:imy;
    ss=im(1:5:imy,1:5:imx);
    [xInput, yInput, zOutput] = prepareSurfaceData( sx, sy, ss );
    bkgI=zOutput<prctile(zOutput,10);
    xInput=xInput(bkgI);
    yInput=yInput(bkgI);
    zOutput=zOutput(bkgI);
    % Set up fittype and options.
    ft = fittype( 'p00+p10*x+p01*y+p20*x^2+p11*x*y+p02*y^2', 'indep', {'x', 'y'}, 'depend', 'z' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf];
    opts.StartPoint = [300 0 0 0 0 0];
    opts.Upper = [Inf Inf Inf Inf Inf Inf];
    % Fit model to data.
    [fitresult, gof] = fit( [xInput, yInput], zOutput, ft, opts );
else
    fitresult=cfit;
end
[x,y]=meshgrid(1:imx,1:imy);
%malmmalian
L=fitresult(x,y)-300;
L=max(L,max(max(L))/2)+300;
%185 is about the dark value. maybe 200?
imbs=im./L -185./L+185./min(L(:)) -1;


end
function [imbs,fitresult]=BactRemoveBackground(im,cfit)
[imy,imx]=size(im);
if isempty(cfit)
    %fit background
    sx=1:5:imx;
    sy=1:5:imy;
    ss=im(1:5:imy,1:5:imx);
    [xInput, yInput, zOutput] = prepareSurfaceData( sx, sy, ss );
    % Set up fittype and options.
    ft = fittype( 'p00+p10*x+p01*y+p20*x^2+p11*x*y+p02*y^2', 'indep', {'x', 'y'}, 'depend', 'z' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf];
    opts.StartPoint = [300 0 0 0 0 0];
    opts.Upper = [Inf Inf Inf Inf Inf Inf];
    % Fit model to data.
    [fitresult, gof] = fit( [xInput, yInput], zOutput, ft, opts );
else
    fitresult=cfit;
end
[x,y]=meshgrid(1:imx,1:imy);
%L=fitresult(x,y)-200;
%L=max(L,max(max(L))/2)+200;
imbs=(im)./200-1;


end
function [imbs,fitresult]=TemplateRemoveBackground(im,cfit,other)
if exist('other','var') && ~isempty(other)
    template=double(imread([cfit{1},other,cfit{3}]));
else
    template=double(imread([cfit{1},cfit{2},cfit{3}]));
end

imbs=im./template;
imbs=imbs./median(reshape(imbs(1:50,1:50),1,[]));
fitresult=imbs./im;
imbs=imbs-1;

end

function cells=SegmentImage(im_bs,SegType)
switch SegType
    case 'bacteria'
        cells=BactSegmentImage(im_bs,SegType);
    case 'bactear'
        cells=BactEarSegmentImage(im_bs,SegType);
    case 'ear'
        cells=EarSegmentImage(im_bs,SegType);
    otherwise
        cells=MamSegmentImage(im_bs,SegType);
end
cells=CalcCellProperties(cells,im_bs);
end
function cells=MamSegmentImage(im_bs,SegType)
im=medfilt2(im_bs,[5,5],'symmetric');

[imy,imx]=size(im);

noise=2*iqr(im(:)-im_bs(:));
segim=zeros(size(im));
ll=localsegment(im,SegType);% | find_blob(im,33);
[ll2,nn]=bwlabel(ll,4);

stats = regionprops(ll2, im, 'MeanIntensity','area');
mmean=[stats.MeanIntensity];
area=[stats.Area];
%sstd=arrayfun(@(x) std(x.PixelValues), stats)';
%mmin=arrayfun(@(x) min(x.PixelValues), stats)';

%calculate statistic about the neighborhood of the cells
stats2 = regionprops(ll2, imdilate(im.*(1-ll),strel('diamond',2)), 'MeanIntensity');
stats3 = regionprops(ll2-imerode(ll2,strel('diamond',1)),'area');

Ssignal=(mmean).*sqrt(area)/noise;
Sedge=mmean./[stats2.MeanIntensity];
Scircular=[stats3.Area]./sqrt(area);

signify=Scircular<6;
signify=signify &((Ssignal>6 & Sedge>3) | (Ssignal>100 & Sedge>2));

% for cseg=find(signify)
%     %segim=segim+imclose((ll2==cseg),strel('disk',3));
%     segim=segim+(ll2==cseg);
% end

segim=ismember(ll2,find(signify));

%cells=mask2cells(imclose(segim,strel('disk',5)));
%cells=mask2cells(imopen(imclose(segim,strel('disk',10)),strel('disk',10)));
cells=mask2cells(segim);
end
function cells=EarSegmentImage(im_bs,SegType)
im=medfilt2(im_bs,[5,5],'symmetric');

[imy,imx]=size(im);

noise=2*iqr(im(:)-im_bs(:));
segim=zeros(size(im));
ll=localsegment(im,SegType);% | find_blob(im,33);
[ll2,nn]=bwlabel(ll,4);

stats = regionprops(ll2, im, 'MeanIntensity','area');
mmean=[stats.MeanIntensity];
area=[stats.Area];
%sstd=arrayfun(@(x) std(x.PixelValues), stats)';
%mmin=arrayfun(@(x) min(x.PixelValues), stats)';

%calculate statistic about the neighborhood of the cells
stats2 = regionprops(ll2, imdilate(im.*(1-ll),strel('diamond',2)), 'MeanIntensity');
stats3 = regionprops(ll2-imerode(ll2,strel('diamond',1)),'area');

Ssignal=(mmean).*sqrt(area)/noise;
Sedge=mmean./[stats2.MeanIntensity];
Scircular=[stats3.Area]./sqrt(area);

signify=Scircular<6;
signify=signify &((Ssignal>6 & Sedge>3) | (Ssignal>100 & Sedge>2));
signify=Ssignal>100;
signify=Ssignal>100 & area<300 & mmean>10;

% for cseg=find(signify)
%     %segim=segim+imclose((ll2==cseg),strel('disk',3));
%     segim=segim+(ll2==cseg);
% end

segim=ismember(ll2,find(signify));

%cells=mask2cells(imclose(segim,strel('disk',5)));
%cells=mask2cells(imopen(imclose(segim,strel('disk',10)),strel('disk',10)));
cells=mask2cells(segim);
end
function cells=BactSegmentImage(im_bs,SegType)
im=medfilt2(im_bs,[5,5],'symmetric');

[imy,imx]=size(im);

noise=2*iqr(im(:)-im_bs(:));
segim=zeros(size(im));
ll=localsegment(im,SegType);% | find_blob(im,33);
[ll2,nn]=bwlabel(ll,4);

stats = regionprops(ll2, im, 'MeanIntensity','area');
mmean=[stats.MeanIntensity];
area=[stats.Area];
%sstd=arrayfun(@(x) std(x.PixelValues), stats)';
%mmin=arrayfun(@(x) min(x.PixelValues), stats)';

%calculate statistic about the neighborhood of the cells
stats2 = regionprops(ll2, imdilate(im.*(1-ll),strel('diamond',2)), 'MeanIntensity');
stats3 = regionprops(ll2-imerode(ll2,strel('diamond',1)),'area');

Ssignal=(mmean).*sqrt(area)/noise;
Sedge=mmean./[stats2.MeanIntensity];
Scircular=[stats3.Area]./sqrt(area);

%signify=Scircular<6;%maybe also Scircular>3;
%signify=signify &((Ssignal>6 & Sedge>3) | (Ssignal>200 & Sedge>2));

signify=((Ssignal/330+Sedge/2.5)>2) & ((Ssignal/300-Scircular/4.3)>0);
%signify=((Ssignal/165+Sedge/2.5)>2) & ((Ssignal/150-Scircular/4.3)>0);

%usefull to find the good cutoffs
%fnum=1;
%if ~isempty(segdata)
%    fnum=segdata(end)+1;
%end
%segdata=[segdata; Sedge',Ssignal',Scircular',signify',fnum*ones(size(signify))'];
%scatter3(ttt(:,1),ttt(:,2),ttt(:,3),1+100*ttt(:,4),ttt(:,4),'.');
for cseg=find(signify)
    %segim=segim+imclose((ll2==cseg),strel('disk',3));
    segim=segim+(ll2==cseg);
end

overlaps=imclose(segim,strel('diamond',1))-segim;
overlapsCC=bwconncomp(overlaps-(imfilter(overlaps,[0,1,0;1,1,1;0,1,0])>3),4);
stat=regionprops(overlapsCC, im, 'MaxIntensity','MinIntensity','PixelValues');
idx=arrayfun(@(x) prctile(x.PixelValues,90), stat)>2;
overlapsCC.NumObjects=sum(idx);
overlapsCC.PixelIdxList=overlapsCC.PixelIdxList(idx);
segim=segim+(labelmatrix(overlapsCC)>0);

%cells=mask2cells(imclose(segim,strel('disk',5)));
%cells=mask2cells(imopen(imclose(segim,strel('disk',10)),strel('disk',10)));
cells=mask2cells(segim);
end
function cells=BactEarSegmentImage(im_bs,SegType)
im=medfilt2(im_bs,[5,5],'symmetric');

[imy,imx]=size(im);

noise=2*iqr(im(:)-im_bs(:));
segim=zeros(size(im));
ll=localsegment(im,SegType);% | find_blob(im,33);
%[ll2,nn]=bwlabel(ll,4);
ll2=ll;

stats = regionprops(ll2, im, 'MeanIntensity','area','Perimeter');
mmean=[stats.MeanIntensity];
area=[stats.Area];
%sstd=arrayfun(@(x) std(x.PixelValues), stats)';
%mmin=arrayfun(@(x) min(x.PixelValues), stats)';

%calculate statistic about the neighborhood of the cells
%stats2 = regionprops(ll2, imdilate(im.*(1-ll),strel('diamond',2)), 'MeanIntensity');
%stats3 = regionprops(ll2-imerode(ll2,strel('diamond',1)),'area');
stats2 = regionprops(edge(ll2,0,'nothinning').*ll2,im,'MeanIntensity');

Ssignal=(mmean).*sqrt(area)/noise;
%Sedge=mmean./[stats2.MeanIntensity];
Sedge=mmean./[stats2.MeanIntensity];
%Scircular=[stats3.Area]./sqrt(area);
Scircular=[stats.Perimeter]./sqrt(area);
global segdata
global gcounter
gcounter=gcounter+1;
segdata=[segdata; [Ssignal;Sedge;Scircular;mmean;area;stats2.MeanIntensity;stats.Perimeter;gcounter*ones(size(Ssignal))]'];
%signify=Scircular<6;%maybe also Scircular>3;
%signify=signify &((Ssignal>6 & Sedge>3) | (Ssignal>200 & Sedge>2));

signify=((Ssignal/330+Sedge/1)>2) & ((Ssignal/300-Scircular/4.3)>0);
%signify=((Ssignal/165+Sedge/2.5)>2) & ((Ssignal/150-Scircular/4.3)>0);

%usefull to find the good cutoffs
%fnum=1;
%if ~isempty(segdata)
%    fnum=segdata(end)+1;
%end
%segdata=[segdata; Sedge',Ssignal',Scircular',signify',fnum*ones(size(signify))'];
%scatter3(ttt(:,1),ttt(:,2),ttt(:,3),1+100*ttt(:,4),ttt(:,4),'.');

segim=ismember(ll2,find(signify));

%dont know what the next segment does.
%overlaps=imclose(segim,strel('diamond',1))-segim;
%overlapsCC=bwconncomp(overlaps-(imfilter(overlaps,[0,1,0;1,1,1;0,1,0])>3),4);
%stat=regionprops(overlapsCC, im, 'MaxIntensity','MinIntensity','PixelValues');
%idx=arrayfun(@(x) prctile(x.PixelValues,90), stat)>2;
%overlapsCC.NumObjects=sum(idx);
%overlapsCC.PixelIdxList=overlapsCC.PixelIdxList(idx);
%segim=segim+(labelmatrix(overlapsCC)>0);%

%cells=mask2cells(imclose(segim,strel('disk',5)));
%cells=mask2cells(imopen(imclose(segim,strel('disk',10)),strel('disk',10)));
cells=mask2cells(segim.*ll2);
end
function imseg=localsegment(im,SegType)
switch SegType
    case 'bacteria'
        imseg=bactlocalsegment(im);
    case 'bactear'
        imseg=bactearlocalsegment(im);
    case 'ear'
        imseg=earlocalsegment(im);
    otherwise
        imseg=mamlocalsegment(im);
        %imseg=MarkerControlledWatershedSegmentation(im);
end
end

function cells=mask2cells(newsegmask)
[imy,imx]=size(newsegmask);
num=max(newsegmask(:));
if num==1
    %input is just a mask
    [newsegmaskl,num]=bwlabel((newsegmask)>0);
else
    %input is a labeled mask
    newsegmaskl=newsegmask;
end
newsegclose=imclose(newsegmask,strel('diamond',1));
cells=struct('mask',{},'pos',{},'size',{},'progenitor',{},'descendants',{});
for cseg=1:num
    csegclose=imdilate(newsegmaskl==cseg,strel('diamond',1)).*newsegclose;
    %flatten csegclose
    csegclose=csegclose>0;
    bbox=regionprops(2*csegclose,{'Area','BoundingBox'});
    if isempty(bbox)
        continue
    end
    bbox=bbox(end);
    if ~isscalar(bbox)
        continue
    end
    cellslice={floor(bbox.BoundingBox(2))+1:floor(bbox.BoundingBox(2)+bbox.BoundingBox(4)),...
        floor(bbox.BoundingBox(1))+1:floor(bbox.BoundingBox(1)+bbox.BoundingBox(3))};
    cells(end+1).mask=double(csegclose(cellslice{1},cellslice{2}));
    cells(end).pos=floor(bbox.BoundingBox(2:-1:1))+1;
    cells(end).size=size(cells(end).mask);
%    if any(cells(end).pos==1) || any(cells(end).pos-[1,1]+cells(end).size==[imy,imx])
%        cells(end)=[];
%        continue
%    end
end

end
function Iseg = waterSegmentImage(I,  Ibin , minima_depth_thresh )

I = double(I);

%Get an image with basins corresponding to white areas in the bin image
Iseg = -I + max(max(I));

%Suppress all minima below a specified value 'minima_depth_thresh'
Iseg  = imhminlocal(Iseg, minima_depth_thresh);

%Each basin is charachterised by a unique label (i.e. number)
Iseg = watershed(Iseg);

%changing all the labels that were 0 in the original image back to 0.
Iseg(Ibin == 0) = 0;
end
function I2 = imhminlocal(varargin)
% Adaptive h-mnima transform, based on imhmin
% instead of substracting a constant h from I, use h=h(I) a monotone
% non-increasing non-negative function of the grey level in I
% see: ftp://doc.nit.ac.ir/cee/y.baleghi/Advanced%20Image%20Processing/Reference%20Books/3D%20Images%20of%20Materials%20Structures%20Processing%20and%20Analysis.pdf
%
% In practice it substract it h-I/3
% so if h is smaller than max(max(I)) it is set to that value
%
%IMHMIN H-minima transform.
%   I2 = IMHMIN(I,H) suppresses all minima in I whose depth is less than
%   H.  I is an intensity image and H is a nonnegative scalar.
%
%   Regional minima are connected components of pixels with the same
%   intensity value, t, whose external boundary pixels all have a value
%   greater than t.
%
%   By default, IMHMIN uses 8-connected neighborhoods for 2-D images and
%   26-connected neighborhoods for 3-D images.  For higher dimensions,
%   IMHMIN uses CONNDEF(NDIMS(I),'maximal').
%
%   I2 = IMHMIN(I,H,CONN) computes the H-minima transform, where CONN
%   specifies the connectivity.  CONN may have the following scalar
%   values:
%
%       4     two-dimensional four-connected neighborhood
%       8     two-dimensional eight-connected neighborhood
%       6     three-dimensional six-connected neighborhood
%       18    three-dimensional 18-connected neighborhood
%       26    three-dimensional 26-connected neighborhood
%
%   Connectivity may be defined in a more general way for any dimension by
%   using for CONN a 3-by-3-by- ... -by-3 matrix of 0s and 1s.  The 1-valued
%   elements define neighborhood locations relative to the center element of
%   CONN.  CONN must be symmetric about its center element.
%
%   Class support
%   -------------
%   I can be of any nonsparse numeric class and any dimension.  I2 has
%   the same size and class as I.
%
%   Example
%   -------
%       a = 10*ones(10,10);
%       a(2:4,2:4) = 7;  % minima 3 lower than surround
%       a(6:8,6:8) = 2;  % minima 8 lower than surround
%       b = imhmin(a,4); % only the deeper minima survive
%
%   See also CONNDEF, IMEXTENDEDMIN, IMHMAX, IMRECONSTRUCT,
%   IMREGIONALMIN.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.6.4.3 $  $Date: 2005/03/31 16:31:28 $

% Testing notes
% -------------
% I       - N-D, real, full
%         - empty ok
%         - Inf ok
%         - NaNs not allowed
%         - logical flag ignored if present
%
% h       - Numeric scalar; nonnegative; real
%         - Inf ok (doesn't make much sense, though)
%         - NaNs not allowed
%
% conn    - valid connectivity specifier
%
% I2      - same class and size as I

[I,H,conn] = ParseInputs(varargin{:});

%H =  0.1*(max_I+(filter2(filt_gauss,max_I-I)));
%figure;imagesc(H);figure;
%%%%

I = imcomplement(I);

%I2 = imreconstruct(imsubtract(I,h), I, conn);
I2 = imreconstruct(imsubtract(I,H), I, conn);

I2 = imcomplement(I2);

%%%
%%% ParseInputs
%%%
    function [I,H,conn] = ParseInputs(varargin)
        
        if verLessThan('matlab', '7.12.1')
            iptchecknargin(2,3,nargin,mfilename);
        else
            narginchk(2,3)
        end
        
        I = varargin{1};
        H = varargin{2};
        
        if verLessThan('matlab', '7.12.1')
            iptcheckinput(I, {'numeric'}, {'real' 'nonsparse'}, mfilename, 'I', 1);
            iptcheckinput(H, {'numeric'}, {'real'}, mfilename, 'H', 2);
        else
            validateattributes(I, {'numeric'}, {'real' 'nonsparse'}, mfilename, 'I', 1);
            validateattributes(H, {'numeric'}, {'real'}, mfilename, 'H', 2);
        end
        
        H = double(H);
        
        if nargin < 3
            conn = conndef(ndims(I),'maximal');
        else
            conn = varargin{3};
            iptcheckconn(conn, mfilename, 'CONN', 3);
        end
        
    end
end
function imseg=mamlocalsegment(im)

%average_size, about the cell size, reducing noise
average_size=30;%15;

%a more uniform image, keep 1e-5 resolution
Th=imtophat(im,strel('disk',average_size));
Th=round(Th*1e5)/1e5;

%where to put the threshold?
%find where there are the most objects (noise segments)
%and take twice that value
localthreshold=2*(fminsearch(@(t) -max(max(bwlabel((Th+1)>t))),1+mode(Th(Th>0))))-2;
%for every cell take it until it is .5 of its maximum
imseg=Th>localthreshold & (Th./imfilter(Th,fspecial('gauss',50,50)))>0.5;

%now use watershed to split neighbouring cells
localmax=imdilate(Th,strel('diamond',10));
localmax=max(min(localmax(localmax>0)),localmax);
imseg=imseg>0;
imseg=waterSegmentImage(Th./localmax,imseg,0.4)>0;
end
function imseg=bactlocalsegment(im)

%average_size, about the cell size, reducing noise
average_size=100;%15;

%a more uniform image, keep 1e-5 resolution
Th=imtophat(im,strel('disk',average_size));
Th=round(Th*1e5)/1e5;

%where to put the threshold?
%find where there are the most objects (noise segments)
%and take twice that value
localthreshold=2*(fminsearch(@(t) -max(max(bwlabel((Th+1)>t))),1+mode(Th(Th>0))))-2;
%for every cell take it until it is .5 of its maximum
imseg=Th>localthreshold & (Th./imfilter(Th,fspecial('gauss',50,50)))>0.5;

%now use watershed to split neighbouring cells
%localmax=imdilate(Th,strel('diamond',4));
localmax=imdilate(Th,strel('diamond',4));
localmax=max(min(localmax(localmax>0)),localmax);
%imseg=waterSegmentImage(Th./localmax,imseg,0.1)>0;
imseg=waterSegmentImage(Th./localmax,imseg,0.1)>0;

end
function imseg=earlocalsegment(im,seed)
disk1=strel('disk',1);

if exist('seed','var')
    imex=seed;
else    
    imex=imextendedmax(im,2*median(double(im(:))));
end
mask=(watershed(bwdist(imex)))>0;
L=imex.*im;
segarea=0;
for i=1:50
    L=(imdilate(L,disk1).*mask);
end
Lf=(imdilate(L,disk1));
Lf(Lf==0)=max(max(Lf));
imn=imtophat(im./imfilter(Lf,fspecial('gauss',5,5)),strel('disk',50));

%imex=imextendedmax(imfilter(im,fspecial('gauss',3,3)),2*median(double(im(:))));
imex=imdilate(imex,strel('disk',2));
for th=linspace(1,0.3,10)
%maxheight=max(im(imex>0));
%minheight=min(im(imex>0));
%for th=linspace(maxheight,3,10)
    imexarea=0;
    while imexarea<sum(sum(imex))
        imexarea=sum(sum(imex));
        imex=bwmorph(imex,'thicken').*bwmorph(imex,'dilate').*(imn>th) | imex;
    end
end

%imseg=imex;
%thicken leaves 1 pixel chanels and holes. so fix it.
imseg=(imclose(bwlabel(imex),strel('diamond',2))).*(imn>th);
%imseg=imdilate(bwlabel(imex),disk1);
%imsegperim=(imseg-imdilate(imseg,disk1))<0;

end
function imseg=bactearlocalsegment(im,seed)
disk1=strel('disk',1);

if exist('seed','var')
    imex=seed;
else    
    imex=imextendedmax(im,8*median(double(im(:))));
end
mask=(watershed(bwdist(imex)))>0;
L=imex.*im;
segarea=0;
for i=1:50
    L=(imdilate(L,disk1).*mask);
end
Lf=(imdilate(L,disk1));
Lf(Lf==0)=max(max(Lf));
imn=imtophat(im./imfilter(Lf,fspecial('gauss',5,5)),strel('disk',50));

%imex=imextendedmax(imfilter(im,fspecial('gauss',3,3)),2*median(double(im(:))));
imex=imdilate(imex,strel('disk',2));
thmin=fminbnd(@(x) sum(sum(bwperim(imn>x))),0,0.5);
for th=linspace(1,thmin,10)
%maxheight=max(im(imex>0));
%minheight=min(im(imex>0));
%for th=linspace(maxheight,3,10)
    imexarea=0;
    while imexarea<sum(sum(imex))
        imexarea=sum(sum(imex));
        imex=bwmorph(imex,'thicken').*bwmorph(imex,'dilate').*(imn>th) | imex;
    end
end

%imseg=imex;
%thicken leaves 1 pixel chanels and holes. so fix it.
imseg=(imclose(bwlabel(imex),strel('diamond',2))).*(imn>th);
%imseg=imdilate(bwlabel(imex),disk1);
%imsegperim=(imseg-imdilate(imseg,disk1))<0;

end
function imseg=MarkerControlledWatershedSegmentation(im,seed)
if exist('seed','var')
    imex=seed;
else    
    immed=medfilt2(im,[10,10],'symmetric');
    imex=imextendedmax(immed,0);
end

imimp=imimposemin(-im,imex);
imshed=double(watershed(imimp));

%get local information. currently min and max.
CC=bwconncomp(imshed>0);
imsegvalmx=zeros(size(imshed));
imsegvalmn=zeros(size(imshed));

for i=1:length(CC.PixelIdxList)
   imsegvalmx(CC.PixelIdxList{i})=max(im(CC.PixelIdxList{i}));
   imsegvalmn(CC.PixelIdxList{i})=min(im(CC.PixelIdxList{i}));
end

%normalize range and cut at 30%
%imnorm=(im-imsegvalmn)./(imsegvalmx-imsegvalmn);
%imseg=(imnorm>0.3).*(imseg>0);

%find the median level at the perimeter
imnorm=im-imsegvalmn;
threshold=2.6*median(imnorm(bwperim(imshed)));
imseg=(imnorm>threshold).*(imshed>0);

imseg=bwlabel(imseg);
imseg=(ismember(imseg,imseg(imex==1)).*imseg)>0;

end

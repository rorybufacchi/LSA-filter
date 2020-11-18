global gEEG
global gOps
global gH
global gV


% Define necessary EEG variables
gEEG                    = EEG;
gOps                    = ops;
gV.eegTmp               = eegTmp;
gV.filtEEGadapt         = filtEEGadapt;
gV.lbd                  = lbd;
gV.filtEEGstandard      = filtEEGstandard;
gV.nPlChan              = length([gEEG.chanlocs.theta]);
gV.inclEls              = find(cellfun(@(x) ~isempty(x),{gEEG.chanlocs.theta})); % electrodes to include

% Ensure the dimensions are correctly ordered
gV.pltEEG           = pltEEG;
if ismember(1,avDims)
    gV.pltEEG       = permute(gV.pltEEG,[4 1 2 3]);
end

gH.fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
disp('-------------------------------------')
disp('BUSY PLOTTING FIGURE, PLEASE WAIT')
disp('-------------------------------------')

set(gH.fig1,'color','white')


% Make boxes for specifying timpoints
tmpP = [0.91 .7 .05 .05];
gH.txtBox(1) = uicontrol('Style','text',...
    'String','Start (s)',...
    'units','normalized',...
    'BackgroundColor', [1 1 1], ...
    'position',tmpP + [0 .16 0 -.03]);
gH.txtBox(2) = uicontrol('Style','text',...
    'String','End (s)',...
    'units','normalized',...
    'BackgroundColor', [1 1 1], ...
    'position',tmpP + [0 .06 0 -.03]);

gH.tpBox(1) = uicontrol('Style','edit',...
    'String',num2str(gEEG.times(gOps.scalpTP),4),...
    'units','normalized',...
    'position',tmpP + [0 .1 0 0]);
gH.tpBox(2) = uicontrol('Style','edit',...
    'String',num2str(gEEG.times(gOps.scalpTPend),4),...
    'units','normalized',...
    'position',tmpP + [0 0 0 0]);

% Make the 'update timepoints' button
gH.tpChButton = uicontrol('style','push',...
                 'units','normalized',...
                 'position',tmpP - [0 0.07 0 0],...
                 'fontsize',12,...
                 'string','Update',...
                 'callback',{@ClickTpChButton,gH.fig1});
             
% Plot the timecourse data
PlotTimeCourse;
% Set the functions for clicking on the lineplot's backround
set(gH.aH,'ButtonDownFcn',@OnTimecourseClick);
set(gcf, 'WindowButtonUpFcn', @OnMouseRelease);

% Plot the starting intersection line
tmpEv.IntersectionPoint(1) = gEEG.times(gOps.scalpTP);
% % % gOps.scalpTPend = gOps.scalpTP;
tmpEv.IntersectionPoint(2) = 0;
subplot(2,3,2:3);
OnTimecourseClick([],tmpEv);

gV.LastPlotCLicked = '';

% Set the end timepoint and make topoplots - initialise with particular values
OnMouseRelease



% Make the 'run filter' button
EndButton = uicontrol('style','push',...
                 'units','normalized',...
                 'position',[0.8 .02 .1 .05],...
                 'fontsize',14,...
                 'string','Run Filter',...
                 'callback',{@ClickEndButton,gH.fig1});
             
% Function to update start and end times in the box
function UpdateTextBoxes
    global gEEG
    global gOps
    global gH
    
    gH.tpBox(1).String = num2str(gEEG.times(gOps.scalpTP),4);
    gH.tpBox(2).String = num2str(gEEG.times(gOps.scalpTPend),4);
end

% Callback function for the 'update timepoints' button
function ClickTpChButton(varargin)
    global gEEG
    global gOps
    global gH
    global gV
    
    gV.LastPlotCLicked      = 'UpdateButton';
    
    tmpStrt                 = str2num(gH.tpBox(1).String);
    tmpStp                  = str2num(gH.tpBox(2).String);
    
    [~, gOps.scalpTP]       = min(abs(gEEG.times-tmpStrt));
    [~, gOps.scalpTPend]    = min(abs(gEEG.times-tmpStp));
    
    
    % Plot the starting intersection line
    tmpEv.IntersectionPoint(1) = gEEG.times(gOps.scalpTP);
    tmpEv.IntersectionPoint(2) = 0;
    subplot(2,3,2:3);
    OnTimecourseClick([],tmpEv);
    
    gV.LastPlotCLicked      = 'UpdateButton';
    OnMouseRelease;
    
end
             
% Callback function for the 'run filter' button  
function ClickEndButton(varargin)
    global gEEG
    global gOps
    global gH
    global gV
    
    cFig  = varargin{3};
    close(cFig)

end


function OnTimecourseClick(dmy,evnt)
    global gEEG
    global gOps
    global gV   % global variables
    global gH	% global handles
    
    xLims       = xlim;
    yLims       = ylim;
    
    % Update intersection positions
    if exist('evnt','var')
        x           = evnt.IntersectionPoint(1);
        y           = evnt.IntersectionPoint(2);
    else
        x = gEEG.times(gOps.scalpTP);
    end
    
    % if there is already a cross/line plotted, delete it
    children    = get(gca, 'children');
    % Plus one for the reference channel
    linesDiff   =  size(children,1) - size(gV.pltEEG,1) - 1;
    if linesDiff > 0
        delete(children(1:linesDiff));
    end
    
    % plot a cross/line with handles
    hold on;
    % plot(xLims, [y y],'k')
    gH.chTC(2) = plot([x x], yLims,'k');
    set(get(get(gH.chTC(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold off
    
    % Find the index of the nearest timepoint to the click
    [~, gOps.scalpTP] = min(abs(gEEG.times-x));
    
    gV.LastPlotCLicked = 'TimeCourse';

end

function OnChanClick(chanH,evnt)
    global gV
    global gOps
    
    x       = evnt.IntersectionPoint(1);
    y       = evnt.IntersectionPoint(2);
    
    % Find the index of the channel to the click (within some tolerance
    % because otherwise the mirror channels are sometimes the only ones
    % found)
    [chXI]  = find(abs(gV.chX-x) < 0.001);
    [chYI]  = find(abs(gV.chY-y) < 0.001);

    gOps.ref = intersect(chXI,chYI);

    
    % Update the LSA filter
    [lsaEEG lbd]        = LSAfilt(gV.eegTmp,gOps.ref,'adaptive');
    gV.filtEEGadapt     = permute(lsaEEG,[3 1 2]);
    gV.lbd              = permute(lbd,[3 1 2]);
    gV.filtEEGstandard  = permute(LSAfilt(gV.eegTmp,gOps.ref,'standard'),[3 1 2]);    
    % Update the topoplots
    PlotScalpMaps
    
    PlotTimeCourse
    OnTimecourseClick;
    gV.LastPlotCLicked = 'Released';
    OnMouseRelease;
    
    gV.LastPlotCLicked = 'Scalpmap';
    
end

function OnScalpMapClick(tppl,evnt)
    % Set the functions for clicking on a specific channel
    OnChanClick([],evnt)
end

function OnMouseRelease(varargin)
    global gEEG
    global gOps
    global gH
    global gV

    % Only perform dragging if the timecourse was clicked on
    if strcmp(gV.LastPlotCLicked,'TimeCourse')    
        currentPoint = get(gH.aH, 'CurrentPoint');
        x = currentPoint(1);
        
        % Find the index of the nearest timepoint to the mouse release
        [~, gOps.scalpTPend] = min(abs(gEEG.times-x));
        
        UpdateTextBoxes
        
        % Make sure another release can't happen until timecourse is
        % clicked again
        gV.LastPlotCLicked = 'Released';
        
    else
        
        x = gEEG.times(gOps.scalpTPend);
    end

    % Update the assessment of all electrodes as references for the given
    % timepoint
    sn = sign(gOps.scalpTPend - gOps.scalpTP);
    if sn == 0
        sn = 1;
    end
    for iCh =1:gV.nPlChan
        % First dimension is signal channels, second is reference channels
        [tmpEEGadapt(:,:,:,iCh), tmp(:,iCh)] = LSAfilt(nanmean(...
            gV.eegTmp(gOps.scalpTP:sn:gOps.scalpTPend,:,:),1),iCh,'adaptive');
    end
    % First dimension is trials, second is signal channels, third is reference channels
    gV.filtEEGadaptTP = squeeze(tmpEEGadapt);
    % First dimension is signal channels, second is reference channels
    gV.lbdTP = squeeze(tmp);
    % find suggested electrodes
    % number of electrodes for which LSA is applicable
    nApl = sum(gV.lbdTP(gV.inclEls,:) > gOps.lbdThr(1) & gV.lbdTP(gV.inclEls,:) < gOps.lbdThr(2));
    % % %     nApl = sum(gV.lbdTP > gOps.lbdThr(1) & gV.lbdTP < gOps.lbdThr(2),2);
    [nAplS, sInds] = sort(nApl,'descend');
    allEl = gV.inclEls;
    gV.bestLambd = ismember(allEl,sInds(1:ceil(end/4)));
    % Suggested electrodes have many with neat Lambda, and have are
    % 'extreme' in terms of their voltage
    [vS vSInds] = sort(abs(nanmean(gEEG.data(gV.inclEls,gOps.scalpTP:sn:gOps.scalpTPend,:),[2 3])));
    gV.extremeV = ismember(allEl,vSInds((end-ceil(end/4)):end));
    
    % Suggested electrodes should have high total lambda, given lambda < 1
    for iCh = 1:gV.nPlChan
        lSum(iCh) = sum(  find(gV.lbdTP(gV.lbdTP(gV.inclEls,iCh) < gOps.lbdThr(2)))  );
    end
    [sLSum, sLSumInds] = sort(lSum,'descend');
    gV.bestLambdSum = ismember(allEl,sLSumInds(1:ceil(end/4)));
    
    gV.sugEl = gV.bestLambd & gV.extremeV & gV.bestLambdSum;
    
    % Update the topoplots
    PlotScalpMaps;
        
    set(gcf,'CurrentAxes',gH.aH)
    hold on;
    gH.chTC(3) = plot(gH.aH,[x x], ylim(gH.aH),'k');
    set(get(get(gH.chTC(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold off
  
end

function PlotTimeCourse
    global gEEG
    global gOps
    global gV   % global variables
    global gH	% global handles
    
    try
        delete(gH.lgd2);
    end
    
    gH.aH   = subplot(2,3,2:3);
    delete(gH.aH.Children)
    disableDefaultInteractivity(gH.aH);
    % Plot all electrodes
    hold on;
    % Plot the rest of the  electrodes
    plot(gEEG.times,squeeze(gV.pltEEG),'LineWidth',1,'ButtonDownFcn',@OnTimecourseClick);
    % Plot the chosen electrode extra thick
    gH.chTC(1) = plot(gEEG.times,squeeze(gV.pltEEG(gOps.ref,:)),'k','LineWidth',2,'ButtonDownFcn',@OnTimecourseClick);


    xlabel('time (s)');
    ylabel('Potential (\muV)');
    gH.lgd2 = legend(gH.chTC(1),'Current Reference Channel');
    
    title('Click and hold, or fill in values, to select averaging time-window of interest')
    
    hold off;
end

function PlotScalpMaps
    global gEEG
    global gOps
    global gH
    global gV
    
    sn = sign(gOps.scalpTPend - gOps.scalpTP);
    if sn == 0
        sn = 1;
    end
        
    % Plot the topoplot
    gH.tpPl(1) = subplot(2,3,1);
    delete(gH.tpPl(1).Children);
    topoplot(nanmean(gEEG.data (:,gOps.scalpTP:sn:gOps.scalpTPend,gOps.trl),[2 3]),gEEG.chanlocs,'electrodes','ptslabels'); hold off;
    % Move the channel descriptions a little to the right
    for iCh = 1:gV.nPlChan
        tmpS = gH.tpPl(1).Children(iCh);
        tmpF = fields(tmpS);
        if ~isempty(strmatch('Position',tmpF))
            gH.tpPl(1).Children(iCh).Position = gH.tpPl(1).Children(iCh).Position + [0.05 0 0];
        end
    end
    % Extract channel locations
    gV.chX = gH.tpPl(1).Children(gV.nPlChan+1).XData;
    gV.chY = gH.tpPl(1).Children(gV.nPlChan+1).YData;
    % remove small channel markers
    delete(gH.tpPl(1).Children(gV.nPlChan+1))
    % Shift the white ring to see the channels better
    gH.tpPl(1).Children(gV.nPlChan+5).ZData = gH.tpPl(1).Children(gV.nPlChan+5).ZData-0.05;
    % Delete the data
    delete(gH.tpPl(1).Children(gV.nPlChan+[6 7]));
    %----------------------------------------
    % Plot separate channels so that they can be clicked on
    hold on
    PlotChanLocs;
    gH.tpCL(1) = scatter(gV.chX,gV.chY,100,'k','filled'); hold on;
    %----------------------------------------
    title('Channel locations (click to change rererence channel)')
    
    % Plot the marker legend
     lgnNms = {'Current Reference','Suggested reference','LSA is interpretable'};
     lgd = legend(   gH.circs(1:length(gH.circs)),lgnNms(1:length(gH.circs))   );
     basePos = gH.tpPl(1).Position;
     lgd.Position = [basePos(1)*0.5 basePos(2) lgd.Position(3) basePos(4)];
     title(lgd,'Electrode Legend')
    
     hold off;
    
    % Plot the topoplot
    gH.tpPl(2) = subplot(2,3,4);
    topoplot(nanmean(gEEG.data (:,gOps.scalpTP:sn:gOps.scalpTPend,gOps.trl),[2 3]),gEEG.chanlocs,'electrodes','off'); hold off; colorbar
    % Shift the white ring to see the channels better
    gH.tpPl(2).Children(5).ZData = gH.tpPl(2).Children(5).ZData-0.05;
    %----------------------------------------
    % Plot separate channels so that they can be clicked on
    hold on
    PlotChanLocs
    gH.tpCL(2) = scatter(gV.chX,gV.chY,100,'k','filled'); hold on;
    hold off;
    %----------------------------------------
    title('Raw EEG')
    
    gH.tpPl(3) = subplot(2,3,5);
    topoplot(nanmean(gV.filtEEGstandard (:,gOps.scalpTP:sn:gOps.scalpTPend,gOps.trl),[2 3]),gEEG.chanlocs,'electrodes','off'); hold on; colorbar
    gH.tpPl(3).Children(5).ZData = gH.tpPl(3).Children(5).ZData-0.05;
    %----------------------------------------
    hold on
    PlotChanLocs
    gH.tpCL(3) = scatter(gV.chX,gV.chY,100,'k','filled'); hold on;
    hold off;
    %----------------------------------------
    title('Standard Rereferencing')
    gH.tpPl(4) = subplot(2,3,6);
    % Ensure contours are also plotted
    tmpPlt = topoplot(nanmean(gV.filtEEGadaptTP(gOps.trl,:,gOps.ref),[1])',gEEG.chanlocs,'style','contour');
    topoplot(nanmean(gV.filtEEGadaptTP(gOps.trl,:,gOps.ref),[1])',gEEG.chanlocs,'pmask',gV.intp,'style','both'); hold on; colorbar
    % Delete non-useful plots
    delete(gH.tpPl(4).Children(9:14))
    gH.tpPl(4).Children(6).ZData = gH.tpPl(4).Children(6).ZData-0.05;
    % Cange colour of contours so they are less disctracting
    gH.tpPl(4).Children(end).LineColor = [0.8 0.8 0.8];
    %----------------------------------------
    hold on
    PlotChanLocs
    gH.tpCL(4) = scatter(gV.chX,gV.chY,100,'k','filled'); hold on;
    hold off;
    %----------------------------------------
    title('LSA (displaying interpretable data only )')
    
    % Stop the default clicking behaviour so that channels can be clicked
    disableDefaultInteractivity(gH.tpPl(1));
    disableDefaultInteractivity(gH.tpPl(2));
    disableDefaultInteractivity(gH.tpPl(3));
    
    set(gH.tpCL,'ButtonDownFcn',@OnScalpMapClick)
    
%     colormap(BlueWhiteRedDavide)
    colormap redbluecmap
end

function PlotChanLocs
    global gEEG
    global gOps
    global gH
    global gV

    if ~isfield(gV,'chX')
        [eloc, labels, theta, radius, indices] = readlocs(gEEG.chanlocs);
        [gV.chX, gV.chY]                       = pol2cart(deg2rad(theta + 90),radius);
        end
    gH.circs(1) = plot(gV.chX(gOps.ref),gV.chY(gOps.ref),'ok','LineWidth',3,'MarkerSize',20);
    
    % Show electrodes for which LSA is interpretable
    % Lambda between two values and for which the variance is smaller than
    % the reference
    sn = sign(gOps.scalpTPend - gOps.scalpTP);
    if sn == 0
        sn = 1;
    end
    tmpVar = nanvar(nanmean(gEEG.data(gV.inclEls,gOps.scalpTP:sn:gOps.scalpTPend,:),2),[],3);
    gV.intp = gV.lbdTP(gV.inclEls,gOps.ref) > gOps.lbdThr(1) & gV.lbdTP(gV.inclEls,gOps.ref) < gOps.lbdThr(2) & ...
        tmpVar < tmpVar(gOps.ref);
    
    if sum(gV.intp)>0
        plot(gV.chX(gV.intp),gV.chY(gV.intp),'ow','LineWidth',5,'MarkerSize',20);
        gH.circs(3) = plot(gV.chX(gV.intp),gV.chY(gV.intp),'o','LineWidth',.5,'MarkerSize',21,'Color',[1 0 1]);
    end
    
    % Plot suggested Electrode
    if sum(gV.sugEl)>0
        gH.circs(2) = plot(gV.chX(gV.sugEl),gV.chY(gV.sugEl),'o','LineWidth',2,'MarkerSize',25,'Color',[.5 .5 .5]);
    else
        try
            gH.circs(3) = [];
        end
    end
end

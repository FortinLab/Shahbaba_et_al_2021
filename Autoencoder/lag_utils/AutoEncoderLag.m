close all
clear all
clc

%% 
% exportYN = 'N';
exportYN = 'Y';
plotYN = 'Y';
% plotYN = 'N';
% chanceNorm = 'Y';
chanceNorm = 'N';
numReps = 10;
distCalc = 'Mahal';
% distCalc = 'Euclid';
% piWindow = -0.25:0.25:0.5;
% poWindow = -0.5:0.25:0.5;

piWindow = 0:0.25:0.25;
poWindow = -0.25:0.25:0;

odorColors = [44/255, 168/255, 224/255;...
    154/255, 133/255, 122/255;...
    9/255, 161/255, 74/255;...
    128/255, 66/255, 151/255;...
    241/255, 103/255, 36/255];
%%
forestDir = 'D:\Dropbox\02. Work Crap\01. Fortin Lab - Post Doc\00. FL Manuscripts\Shababa et al - Decoding\NPG\Nat Comms\Forest Analysis';
cd(forestDir);
load('Rosetta.mat');
filenames = dir(forestDir);
filenames = {filenames.name};

% Process Poke In Data
piFiles = filenames(cell2mat(cellfun(@(a)~isempty(a), regexp(filenames, '[A-Z|a-z]*_PokeIn'), 'uniformoutput', 0)));
piEncodes = cell(length(piWindow), length(piFiles));
piAnis = cell(1,length(piFiles));
for a = 1:length(piFiles)
    flParts = strsplit(piFiles{a}, '_');
    piAnis{a} = flParts{1};
    % {1:StartTime} {2:CenterTime} {3:TrialID} {4:Trial} {5:Odor} {6:IsInSeq} {7:IsCorrect} {8:Dimenstion1} {9:Dimension2}
    tempEncodes = table2array(readtable(piFiles{a}));
    % %%Comment In To Remove Odor E &/or Pos 5%% %
%     oEorPos5Log = tempEncodes(:,4)==4 | tempEncodes(:,5)==4;
%     tempEncodes(oEorPos5Log,:) = [];
%     %%Comment In To Remove Odor A &/or Pos 1%% %
%     oAorPosAlog = tempEncodes(:,4)==0 | tempEncodes(:,5)==0;
%     tempEncodes(oAorPosAlog,:) = [];
    % %%Comment In To Z-normalize the embeddings
    tempEncodes(:,8) = zscore(tempEncodes(:,8));
    tempEncodes(:,9) = zscore(tempEncodes(:,9));
    for t = 1:length(piWindow)
        piEncodes{t,a} = tempEncodes(tempEncodes(:,1)==piWindow(t),:);
    end
end
[piDistsTrial, piLagDistsTrial, piDists, piLagDists, piRawDists, piRawLagDists] = CalculateInSeqLagDist(piEncodes, piWindow, piAnis, 'PokeIn', distCalc);
[piDistsTrialLagPOS, piDistsTrialLagODR, piDistsPOS, piDistsODR, piLagDistsPOS, piLagDistsODR, piDistsRawLagPOS, piDistsRawLagODR] = CalculateOutSeqLagDist(piEncodes, distCalc);
[piAniOSC, piAniOSI,...
    piAniOSCtrialDistPOS, piAniOSItrialDistPOS,...
    piAniOSCtrialDistODR, piAniOSItrialDistODR,...
    piAniOSCtrialLagDistPOS, piAniOSItrialLagDistPOS,...
    piAniOSCtrialLagDistODR, piAniOSItrialLagDistODR,...
    piAniOSCrawLagDistPOS, piAniOSIrawLagDistPOS,...
    piAniOSCrawLagDistODR, piAniOSIrawLagDistODR] = CalculateOutSeqTrialDist(piEncodes, distCalc);
if strcmp(plotYN, 'N')
    close all
end

% Process Poke Out Data
poFiles = filenames(cell2mat(cellfun(@(a)~isempty(a), regexp(filenames, '[A-z|a-z]*_PokeOut'), 'uniformoutput', 0)));
poEncodes = cell(length(poWindow), length(poFiles));
poAnis = cell(1,length(poFiles));
for a = 1:length(poFiles)
    flParts = strsplit(poFiles{a}, '_');
    poAnis{a} = flParts{1};
    tempEncodes = table2array(readtable(poFiles{a}));
    % %%Comment In To Remove Odor E &/or Pos 5%% %
%     oEorPos5Log = tempEncodes(:,4)==4 | tempEncodes(:,5)==4;
%     tempEncodes(oEorPos5Log,:) = [];
    % %%Comment In To Remove Odor A &/or Pos 1%% %
    oAorPos1log = tempEncodes(:,4)==0 | tempEncodes(:,5)==0;
    tempEncodes(oAorPos1log,:) = [];
    % %%Comment In To Z-normalize the embeddings
    tempEncodes(:,8) = zscore(tempEncodes(:,8));
    tempEncodes(:,9) = zscore(tempEncodes(:,9));
    for t = 1:length(poWindow)
        poEncodes{t,a} = tempEncodes(tempEncodes(:,1)==poWindow(t),:);
    end
end
[poDistsTrial, poLagDistsTrial, poDists, poLagDists, poRawDists, poRawLagDists] = CalculateInSeqLagDist(poEncodes, poWindow, poAnis, 'PokeOut', distCalc);
[poDistsTrialLagPOS, poDistsTrialLagODR, poDistsPOS, poDistsODR, poLagDistsPOS, poLagDistsODR, poDistsRawLagPOS, poDistsRawLagODR] = CalculateOutSeqLagDist(poEncodes, distCalc);
[poAniOSC, poAniOSI,...
    poAniOSCtrialDistPOS, poAniOSItrialDistPOS,...
    poAniOSCtrialDistODR, poAniOSItrialDistODR,...
    poAniOSCtrialLagDistPOS, poAniOSItrialLagDistPOS,...
    poAniOSCtrialLagDistODR, poAniOSItrialLagDistODR,...
    poAniOSCrawLagDistPOS, poAniOSIrawLagDistPOS,...
    poAniOSCrawLagDistODR, poAniOSIrawLagDistODR] = CalculateOutSeqTrialDist(poEncodes, distCalc);
if strcmp(plotYN, 'N')
    close all
end

%% Normalize to Chance
if strcmp(chanceNorm, 'Y')
    chancePIlagDistsTrial = repmat({nan(1,numReps)}, size(piLagDistsTrial));
    chancePIdistsTrialLagPOS = repmat({nan(1,numReps)}, size(piDistsTrialLagPOS));
    chancePIdistsTrialLagODR = repmat({nan(1,numReps)}, size(piDistsTrialLagODR));
    chancePIaniOSCtrialLagDistPOS = repmat({nan(1,numReps)}, size(piAniOSCtrialLagDistPOS));
    chancePIaniOSCtrialLagDistODR = repmat({nan(1,numReps)}, size(piAniOSCtrialLagDistODR));
    chancePIaniOSItrialLagDistPOS = repmat({nan(1,numReps)}, size(piAniOSItrialLagDistPOS));
    chancePIaniOSItrialLagDistODR = repmat({nan(1,numReps)}, size(piAniOSItrialLagDistODR));
    
    chancePOlagDistsTrial = repmat({nan(1,numReps)}, size(poLagDistsTrial));
    chancePOdistsTrialLagPOS = repmat({nan(1,numReps)}, size(poDistsTrialLagPOS));
    chancePOdistsTrialLagODR = repmat({nan(1,numReps)}, size(poDistsTrialLagODR));
    chancePOaniOSCtrialLagDistPOS = repmat({nan(1,numReps)}, size(poAniOSCtrialLagDistPOS));
    chancePOaniOSCtrialLagDistODR = repmat({nan(1,numReps)}, size(poAniOSCtrialLagDistODR));
    chancePOaniOSItrialLagDistPOS = repmat({nan(1,numReps)}, size(poAniOSItrialLagDistPOS));
    chancePOaniOSItrialLagDistODR = repmat({nan(1,numReps)}, size(poAniOSItrialLagDistODR));
    analyStart = tic;
    for rep = 1:numReps 
        assignin('base', 'rep', rep);
        repStart = tic;
        tempPIencodes = piEncodes;
        for pi = 1:length(tempPIencodes(:))
            trialPermCols = sortrows([randperm(size(piEncodes{pi},1))', piEncodes{pi}(:,1:7)]);
            tempPIencodes{pi}(:,1:7) = trialPermCols(:,2:end);
        end
        [~, tempPiLagDistsTrial, ~, ~, ~, ~] = CalculateInSeqLagDist(tempPIencodes, piWindow, piAnis, 'PokeIn', distCalc); close all;
        [tempPiDistsTrialLagPOS, tempPiDistsTrialLagODR, ~, ~, ~, ~, ~, ~] = CalculateOutSeqLagDist(tempPIencodes, distCalc);
        [~, ~,...
            ~, ~,...
            ~, ~,...
            tempPiAniOSCtrialLagDistPOS, tempPiAniOSItrialLagDistPOS,...
            tempPiAniOSCtrialLagDistODR, tempPiAniOSItrialLagDistODR] = CalculateOutSeqTrialDist(tempPIencodes, distCalc);
        for pi = 1:length(chancePIlagDistsTrial(:))
            chancePIlagDistsTrial{pi}(rep) = nanmean(tempPiLagDistsTrial{pi});
            chancePIdistsTrialLagPOS{pi}(rep) = nanmean(tempPiDistsTrialLagPOS{pi});
            chancePIdistsTrialLagODR{pi}(rep) = nanmean(tempPiDistsTrialLagODR{pi});
            chancePIaniOSCtrialLagDistPOS{pi}(rep) = nanmean(tempPiAniOSCtrialLagDistPOS{pi});
            chancePIaniOSCtrialLagDistODR{pi}(rep) = nanmean(tempPiAniOSCtrialLagDistODR{pi});
            chancePIaniOSItrialLagDistPOS{pi}(rep) = nanmean(tempPiAniOSItrialLagDistPOS{pi});
            chancePIaniOSItrialLagDistODR{pi}(rep) = nanmean(tempPiAniOSItrialLagDistODR{pi});
        end
        
        tempPOencodes = poEncodes;
        for po = 1:length(tempPOencodes(:))
            trialPermCols = sortrows([randperm(size(poEncodes{po},1))', poEncodes{po}(:,1:7)]);
            tempPOencodes{po}(:,1:7) = trialPermCols(:,2:end);
        end
        
        [~, tempPoLagDistsTrial, ~, ~, ~, ~] = CalculateInSeqLagDist(poEncodes, poWindow, poAnis, 'PokeOut', distCalc); close all;
        [tempPoDistsTrialLagPOS, tempPoDistsTrialLagODR, ~, ~, ~, ~, ~, ~] = CalculateOutSeqLagDist(poEncodes, distCalc);
        [~, ~,...
            ~, ~,...
            ~, ~,...
            tempPoAniOSCtrialLagDistPOS, tempPoAniOSItrialLagDistPOS,...
            tempPoAniOSCtrialLagDistODR, tempPoAniOSItrialLagDistODR] = CalculateOutSeqTrialDist(poEncodes, distCalc);
        for po = 1:length(poLagDistsTrial(:))
            chancePOlagDistsTrial{po}(rep) = nanmean(tempPoLagDistsTrial{po});
            chancePOdistsTrialLagPOS{po}(rep) = nanmean(tempPoDistsTrialLagPOS{po});
            chancePOdistsTrialLagODR{po}(rep) = nanmean(tempPoDistsTrialLagODR{po});
            chancePOaniOSCtrialLagDistPOS{po}(rep) = nanmean(tempPoAniOSCtrialLagDistPOS{po});
            chancePOaniOSCtrialLagDistODR{po}(rep) = nanmean(tempPoAniOSCtrialLagDistODR{po});
            chancePOaniOSItrialLagDistPOS{po}(rep) = nanmean(tempPoAniOSItrialLagDistPOS{po});
            chancePOaniOSItrialLagDistODR{po}(rep) = nanmean(tempPoAniOSItrialLagDistODR{po});
        end
        if mod(rep,10)==0
            fprintf('Rep %i/%i Done in %.02fs, %.02fmin elapsed\n', rep, numReps, toc(repStart), toc(analyStart)/60);
        end
    end
    fprintf('Permutations complete, it took %.02fmin\n', toc(analyStart)/60);
    
    for r = 1:size(chancePIlagDistsTrial,1)
        for c = 1:size(chancePIlagDistsTrial,2)            
            piLagDistsTrial{r,c} = (piLagDistsTrial{r,c}-mean(chancePIlagDistsTrial{r,c}, 'omitnan'))./std(chancePIlagDistsTrial{r,c}, 'omitnan');
            piDistsTrialLagPOS{r,c} = (piDistsTrialLagPOS{r,c}-mean(chancePIdistsTrialLagPOS{r,c}, 'omitnan'))./std(chancePIdistsTrialLagPOS{r,c}, 'omitnan');
            piDistsTrialLagODR{r,c} = (piDistsTrialLagODR{r,c}-mean(chancePIdistsTrialLagODR{r,c}, 'omitnan'))./std(chancePIdistsTrialLagODR{r,c}, 'omitnan');
            piAniOSCtrialLagDistPOS{r,c} = (piAniOSCtrialLagDistPOS{r,c}-mean(chancePIaniOSCtrialLagDistPOS{r,c}, 'omitnan'))./std(chancePIaniOSCtrialLagDistPOS{r,c}, 'omitnan');
            piAniOSCtrialLagDistODR{r,c} = (piAniOSCtrialLagDistODR{r,c}-mean(chancePIaniOSCtrialLagDistODR{r,c}, 'omitnan'))./std(chancePIaniOSCtrialLagDistODR{r,c}, 'omitnan');
            piAniOSItrialLagDistPOS{r,c} = (piAniOSItrialLagDistPOS{r,c}-mean(chancePIaniOSItrialLagDistPOS{r,c}, 'omitnan'))./std(chancePIaniOSItrialLagDistPOS{r,c}, 'omitnan');
            piAniOSItrialLagDistODR{r,c} = (piAniOSItrialLagDistODR{r,c}-mean(chancePIaniOSItrialLagDistODR{r,c}, 'omitnan'))./std(chancePIaniOSItrialLagDistODR{r,c}, 'omitnan');
        end
    end
    for r = 1:size(poLagDistsTrial,1)
        for c = 1:size(poLagDistsTrial,2)
            poLagDistsTrial{r,c} = (poLagDistsTrial{r,c}-mean(chancePOlagDistsTrial{r,c}, 'omitnan'))./std(chancePOlagDistsTrial{r,c}, 'omitnan');
            poDistsTrialLagPOS{r,c} = (poDistsTrialLagPOS{r,c}-mean(chancePOdistsTrialLagPOS{r,c}, 'omitnan'))./std(chancePOdistsTrialLagPOS{r,c}, 'omitnan');
            poDistsTrialLagODR{r,c} = (poDistsTrialLagODR{r,c}-mean(chancePOdistsTrialLagODR{r,c}, 'omitnan'))./std(chancePOdistsTrialLagODR{r,c}, 'omitnan');
            poAniOSCtrialLagDistPOS{r,c} = (poAniOSCtrialLagDistPOS{r,c}-mean(chancePOaniOSCtrialLagDistPOS{r,c}, 'omitnan'))./std(chancePOaniOSCtrialLagDistPOS{r,c}, 'omitnan');
            poAniOSCtrialLagDistODR{r,c} = (poAniOSCtrialLagDistODR{r,c}-mean(chancePOaniOSCtrialLagDistODR{r,c}, 'omitnan'))./std(chancePOaniOSCtrialLagDistODR{r,c}, 'omitnan');
            poAniOSItrialLagDistPOS{r,c} = (poAniOSItrialLagDistPOS{r,c}-mean(chancePOaniOSItrialLagDistPOS{r,c}, 'omitnan'))./std(chancePOaniOSItrialLagDistPOS{r,c}, 'omitnan');
            poAniOSItrialLagDistODR{r,c} = (poAniOSItrialLagDistODR{r,c}-mean(chancePOaniOSItrialLagDistODR{r,c}, 'omitnan'))./std(chancePOaniOSItrialLagDistODR{r,c}, 'omitnan');
        end
    end
end
%% Export Data
if strcmp(exportYN,'Y')
    if strcmp(chanceNorm, 'N')
% *********************************************Export INSEQ Lag
        ExportMatCSVlag(piLagDistsTrial, piWindow, sprintf('InSeq_LagDistsTrial_PI_%s', distCalc));
        ExportMatCSVlag(poLagDistsTrial, poWindow, sprintf('InSeq_LagDistsTrial_PO_%s', distCalc));

        ExportMatCSVlag(piRawLagDists, piWindow, sprintf('InSeq_LagDistsRaw_PI_%s', distCalc));
        ExportMatCSVlag(poRawLagDists, poWindow, sprintf('InSeq_LagDistsRaw_PO_%s', distCalc));
% *********************************************Export OUTSEQ vs INSEQ Lag
        ExportMatCSVlag(piDistsTrialLagPOS, piWindow, sprintf('OutSeqVsInSeq_LagDistsTrial_POS_PI_%s', distCalc));
        ExportMatCSVlag(piDistsTrialLagODR, piWindow, sprintf('OutSeqVsInSeq_LagDistsTrial_ODR_PI_%s', distCalc));
        
        ExportMatCSVlag(poDistsTrialLagPOS, poWindow, sprintf('OutSeqVsInSeq_LagDistsTrial_POS_PO_%s', distCalc));
        ExportMatCSVlag(poDistsTrialLagODR, poWindow, sprintf('OutSeqVsInSeq_LagDistsTrial_ODR_PO_%s', distCalc));
        
        
        ExportMatCSVlag(piDistsRawLagPOS, piWindow, sprintf('OutSeqVsInSeq_LagDistsRaw_POS_PI_%s', distCalc));
        ExportMatCSVlag(piDistsRawLagODR, piWindow, sprintf('OutSeqVsInSeq_LagDistsRaw_ODR_PI_%s', distCalc));
        
        ExportMatCSVlag(poDistsRawLagPOS, poWindow, sprintf('OutSeqVsInSeq_LagDistsRaw_POS_PO_%s', distCalc));
        ExportMatCSVlag(poDistsRawLagODR, poWindow, sprintf('OutSeqVsInSeq_LagDistsRaw_ODR_PO_%s', distCalc));
% *********************************************Export OUTSEQ vs OUTSEQ Lag
        ExportMatCSVlag(piAniOSCtrialLagDistPOS, piWindow, sprintf('OSC_TrialLagDist_POS_PI_%s', distCalc));
        ExportMatCSVlag(piAniOSCtrialLagDistODR, piWindow, sprintf('OSC_TrialLagDist_ODR_PI_%s', distCalc));
        ExportMatCSVlag(piAniOSItrialLagDistPOS, piWindow, sprintf('OSI_TrialLagDist_POS_PI_%s', distCalc));
        ExportMatCSVlag(piAniOSItrialLagDistODR, piWindow, sprintf('OSI_TrialLagDist_ODR_PI_%s', distCalc));
        
        ExportMatCSVlag(poAniOSCtrialLagDistPOS, poWindow, sprintf('OSC_TrialLagDist_POS_PO_%s', distCalc));
        ExportMatCSVlag(poAniOSCtrialLagDistODR, poWindow, sprintf('OSC_TrialLagDist_ODR_PO_%s', distCalc));
        ExportMatCSVlag(poAniOSItrialLagDistPOS, poWindow, sprintf('OSI_TrialLagDist_POS_PO_%s', distCalc));
        ExportMatCSVlag(poAniOSItrialLagDistODR, poWindow, sprintf('OSI_TrialLagDist_ODR_PO_%s', distCalc));
        
        
        ExportMatCSVlag(piAniOSCrawLagDistPOS, piWindow, sprintf('OSC_RawLagDist_POS_PI_%s', distCalc));
        ExportMatCSVlag(piAniOSCrawLagDistODR, piWindow, sprintf('OSC_RawLagDist_ODR_PI_%s', distCalc));
        ExportMatCSVlag(piAniOSIrawLagDistPOS, piWindow, sprintf('OSI_RawLagDist_POS_PI_%s', distCalc));
        ExportMatCSVlag(piAniOSIrawLagDistODR, piWindow, sprintf('OSI_RawLagDist_ODR_PI_%s', distCalc));
        
        ExportMatCSVlag(poAniOSCrawLagDistPOS, poWindow, sprintf('OSC_RawLagDist_POS_PO_%s', distCalc));
        ExportMatCSVlag(poAniOSCrawLagDistODR, poWindow, sprintf('OSC_RawLagDist_ODR_PO_%s', distCalc));
        ExportMatCSVlag(poAniOSIrawLagDistPOS, poWindow, sprintf('OSI_RawLagDist_POS_PO_%s', distCalc));
        ExportMatCSVlag(poAniOSIrawLagDistODR, poWindow, sprintf('OSI_RawLagDist_ODR_PO_%s', distCalc));
    else        
% *********************************************Export INSEQ Lag
        ExportMatCSVlag(piLagDistsTrial, piWindow, sprintf('InSeq_LagDistsTrial_PI_%s_chanceNorm%i', distCalc, numReps));
        ExportMatCSVlag(poLagDistsTrial, poWindow, sprintf('InSeq_LagDistsTrial_PO_%s_chanceNorm%i', distCalc, numReps));
% *********************************************Export OUTSEQ vs INSEQ Lag
        ExportMatCSVlag(piDistsTrialLagPOS, piWindow, sprintf('OutSeqVsInSeq_LagDistsTrial_POS_PI_%s_chanceNorm%i', distCalc, numReps));
        ExportMatCSVlag(piDistsTrialLagODR, piWindow, sprintf('OutSeqVsInSeq_LagDistsTrial_ODR_PI_%s_chanceNorm%i', distCalc, numReps));
        
        ExportMatCSVlag(poDistsTrialLagPOS, poWindow, sprintf('OutSeqVsInSeq_LagDistsTrial_POS_PO_%s_chanceNorm%i', distCalc, numReps));
        ExportMatCSVlag(poDistsTrialLagODR, poWindow, sprintf('OutSeqVsInSeq_LagDistsTrial_ODR_PO_%s_chanceNorm%i', distCalc, numReps));
% *********************************************Export OUTSEQ vs OUTSEQ Lag
        ExportMatCSVlag(piAniOSCtrialLagDistPOS, piWindow, sprintf('OSC_TrialLagDist_POS_PI_%s_chanceNorm%i', distCalc, numReps));
        ExportMatCSVlag(piAniOSCtrialLagDistODR, piWindow, sprintf('OSC_TrialLagDist_ODR_PI_%s_chanceNorm%i', distCalc, numReps));
        ExportMatCSVlag(piAniOSItrialLagDistPOS, piWindow, sprintf('OSI_TrialLagDist_POS_PI_%s_chanceNorm%i', distCalc, numReps));
        ExportMatCSVlag(piAniOSItrialLagDistODR, piWindow, sprintf('OSI_TrialLagDist_ODR_PI_%s_chanceNorm%i', distCalc, numReps));
        
        ExportMatCSVlag(poAniOSCtrialLagDistPOS, poWindow, sprintf('OSC_TrialLagDist_POS_PO_%s_chanceNorm%i', distCalc, numReps));
        ExportMatCSVlag(poAniOSCtrialLagDistODR, poWindow, sprintf('OSC_TrialLagDist_ODR_PO_%s_chanceNorm%i', distCalc, numReps));
        ExportMatCSVlag(poAniOSItrialLagDistPOS, poWindow, sprintf('OSI_TrialLagDist_POS_PO_%s_chanceNorm%i', distCalc, numReps));
        ExportMatCSVlag(poAniOSItrialLagDistODR, poWindow, sprintf('OSI_TrialLagDist_ODR_PO_%s_chanceNorm%i', distCalc, numReps));
    end
end

%% Plot Summary Figures
if strcmp(plotYN, 'Y')
    % PokeIn lag distance
    figure;
    annotation(gcf,'textbox', [0 0.93 1 0.05],'String', 'ISC PokeIn Aligned Mean Lag Distance',...
        'FontWeight','bold', 'FontSize',12, 'edgecolor', 'none', 'horizontalalignment', 'left');
    for t = 1:size(piLagDists,3)
        curMtx = mean(cell2mat(reshape(cellfun(@(a)a(:,:,t), piDists, 'uniformoutput', 0), [1,1,length(piDists)])),3);
        subplot(size(piLagDists,3),2, sub2ind([2,size(piLagDists,3)], 1,t));
        imagesc(curMtx);
        [r,c] = find(isnan(curMtx));
        for p = 1:length(r)
            patch([c(p)-0.5 c(p)+0.5 c(p)+0.5 c(p)-0.5], [r(p)+0.5 r(p)+0.5 r(p)-0.5 r(p)-0.5],'white', 'EdgeColor', 'white');
        end
        set(gca, 'xtick', 1:5, 'xticklabel', Rosetta(1:5), 'ytick', 1:5, 'yticklabel', Rosetta(1:5));
        subplot(size(piLagDists,3),2, sub2ind([2,size(piLagDists,3)], 2,t));
        BarPlotErrorbars(mean(piLagDists(:,:,t)), std(piLagDists(:,:,t)), {'XTick'}, -4:4);
        title(piWindow(t));
    end
    
    % PokeIn OutSeq Lag Distance
    figure;
    annotation(gcf,'textbox', [0 0.93 1 0.05],'String', 'OSC PokeIn Aligned Mean Lag Distance',...
        'FontWeight','bold', 'FontSize',12, 'edgecolor', 'none', 'horizontalalignment', 'left');
    for t = 1:size(piLagDistsPOS,3)
        curMtx = nanmean(cell2mat(reshape(cellfun(@(a)a(:,:,t), piDistsPOS, 'uniformoutput', 0), [1,1,length(piDistsPOS)])),3);
        subplot(size(piLagDistsPOS,3),4, sub2ind([4,size(piLagDistsPOS,3)], 1,t));
        imagesc(curMtx);
        [r,c] = find(isnan(curMtx));
        for p = 1:length(r)
            patch([c(p)-0.5 c(p)+0.5 c(p)+0.5 c(p)-0.5], [r(p)+0.5 r(p)+0.5 r(p)-0.5 r(p)-0.5],'white', 'EdgeColor', 'white');
        end
        set(gca, 'xtick', 1:5, 'xticklabel', 1:5, 'ytick', 1:5, 'yticklabel', Rosetta(1:5));
        ylabel('Odor');
        xlabel('Position');
        subplot(size(piLagDistsPOS,3),4, sub2ind([4,size(piLagDistsPOS,3)], 2,t));
        BarPlotErrorbars(mean(piLagDistsPOS(:,:,t)), std(piLagDistsPOS(:,:,t)), {'XTick'}, -4:4);
        title(piWindow(t));
    end
    for t = 1:size(piLagDistsODR,3)
        curMtx = nanmean(cell2mat(reshape(cellfun(@(a)a(:,:,t), piDistsODR, 'uniformoutput', 0), [1,1,length(piDistsODR)])),3);
        subplot(size(piLagDistsODR,3),4, sub2ind([4,size(piLagDistsODR,3)], 3,t));
        imagesc(curMtx);
        [r,c] = find(isnan(curMtx));
        for p = 1:length(r)
            patch([c(p)-0.5 c(p)+0.5 c(p)+0.5 c(p)-0.5], [r(p)+0.5 r(p)+0.5 r(p)-0.5 r(p)-0.5],'white', 'EdgeColor', 'white');
        end
        set(gca, 'xtick', 1:5, 'xticklabel', 1:5, 'ytick', 1:5, 'yticklabel', Rosetta(1:5));
        ylabel('Odor');
        xlabel('Position');
        subplot(size(piLagDistsODR,3),4, sub2ind([4,size(piLagDistsODR,3)], 4,t));
        BarPlotErrorbars(mean(piLagDistsODR(:,:,t)), std(piLagDistsODR(:,:,t)), {'XTick'}, -4:4);
        title(piWindow(t));
    end
    
    % PokeOut lag distance
    figure;
    annotation(gcf,'textbox', [0 0.93 1 0.05],'String', 'ISC PokeOut Aligned Mean Lag Distance',...
        'FontWeight','bold', 'FontSize',12, 'edgecolor', 'none', 'horizontalalignment', 'left');
    for t = 1:size(poLagDists,3)
        curMtx = mean(cell2mat(reshape(cellfun(@(a)a(:,:,t), poDists, 'uniformoutput', 0), [1,1,length(poDists)])),3);
        subplot(size(poLagDists,3),2, sub2ind([2,size(poLagDists,3)], 1,t));
        imagesc(curMtx);
        set(gca, 'xtick', 1:5, 'xticklabel', Rosetta(1:5), 'ytick', 1:5, 'yticklabel', Rosetta(1:5));
        subplot(size(poLagDists,3),2, sub2ind([2,size(poLagDists,3)], 2,t));
        BarPlotErrorbars(mean(poLagDists(:,:,t)), std(poLagDists(:,:,t)), {'XTick'}, -4:4);
        title(poWindow(t));
    end
    
    % PokeOut OutSeq Lag Distance
    figure;
    annotation(gcf,'textbox', [0 0.93 1 0.05],'String', 'OSC PokeOut Aligned Mean Lag Distance',...
        'FontWeight','bold', 'FontSize',12, 'edgecolor', 'none', 'horizontalalignment', 'left');
    for t = 1:size(poLagDistsPOS,3)
        curMtx = nanmean(cell2mat(reshape(cellfun(@(a)a(:,:,t), poDistsPOS, 'uniformoutput', 0), [1,1,length(poDistsPOS)])),3);
        subplot(size(poLagDistsPOS,3),4, sub2ind([4,size(poLagDistsPOS,3)], 1,t));
        imagesc(curMtx);
        [r,c] = find(isnan(curMtx));
        for p = 1:length(r)
            patch([c(p)-0.5 c(p)+0.5 c(p)+0.5 c(p)-0.5], [r(p)+0.5 r(p)+0.5 r(p)-0.5 r(p)-0.5],'white', 'EdgeColor', 'white');
        end
        set(gca, 'xtick', 1:5, 'xticklabel', 1:5, 'ytick', 1:5, 'yticklabel', Rosetta(1:5));
        ylabel('Odor');
        xlabel('Position');
        subplot(size(poLagDistsPOS,3),4, sub2ind([4,size(poLagDistsPOS,3)], 2,t));
        BarPlotErrorbars(mean(poLagDistsPOS(:,:,t)), std(poLagDistsPOS(:,:,t)), {'XTick'}, -4:4);
        title(poWindow(t));
    end
    for t = 1:size(poLagDistsODR,3)
        curMtx = nanmean(cell2mat(reshape(cellfun(@(a)a(:,:,t), poDistsODR, 'uniformoutput', 0), [1,1,length(poDistsODR)])),3);
        subplot(size(poLagDistsODR,3),4, sub2ind([4,size(poLagDistsODR,3)], 3,t));
        imagesc(curMtx);
        [r,c] = find(isnan(curMtx));
        for p = 1:length(r)
            patch([c(p)-0.5 c(p)+0.5 c(p)+0.5 c(p)-0.5], [r(p)+0.5 r(p)+0.5 r(p)-0.5 r(p)-0.5],'white', 'EdgeColor', 'white');
        end
        set(gca, 'xtick', 1:5, 'xticklabel', 1:5, 'ytick', 1:5, 'yticklabel', Rosetta(1:5));
        ylabel('Odor');
        xlabel('Position');
        subplot(size(poLagDistsODR,3),4, sub2ind([4,size(poLagDistsODR,3)], 4,t));
        BarPlotErrorbars(mean(poLagDistsODR(:,:,t)), std(poLagDistsODR(:,:,t)), {'XTick'}, -4:4);
        title(poWindow(t));
    end
    
    %% Plot ani outSeq Correct
    figure;
    for t = 1:size(piAniOSC,1)
        tempOSC = piAniOSC{t};
        for a = 1:size(tempOSC,3)
            subplot(size(piAniOSC,1), size(tempOSC,3), (t-1)*size(tempOSC,3)+a);
            for p = 1:size(tempOSC,2)
                for o = 1:size(tempOSC,1)
                    if ~isempty(tempOSC{o,p,a})
                        scatter(tempOSC{o,p,a}(:,1), tempOSC{o,p,a}(:,2), 15,...
                            'markerfacecolor', odorColors(o,:), 'markeredgecolor', odorColors(p,:),...
                            'marker', 'o', 'linewidth', 1);
                        hold on;
                        %                     for tp = 1:16:size(tempOSC{o,p,a},1)
                        %                         plot(tempOSC{o,p,a}(tp:tp+15,1), tempOSC{o,p,a}(tp:tp+15,2), '-k');
                        %                     end
                    end
                end
            end
            if a == 1
                ylabel(piWindow(t), 'fontweight', 'bold');
            end
            if t == 1
                title(piAnis(a));
            end
            drawnow
        end
    end
    
    PlotOSscatters(piAnis, piEncodes, piAniOSC, piWindow, odorColors, 'PokeIn OutSeqCorrect')
    PlotOSscatters(poAnis, poEncodes, poAniOSC, poWindow, odorColors, 'PokeOut OutSeqCorrect');
end
%%



%% CalculateInSeqLagDist
function [trialDists, trialLagDists, aniDists, aniLagDists, rawDists, rawLagDists] = CalculateInSeqLagDist(encodes, window, aniID, alignment, distCalc)
odorColors = [44/255, 168/255, 224/255;...
    154/255, 133/255, 122/255;...
    9/255, 161/255, 74/255;...
    128/255, 66/255, 151/255;...
    241/255, 103/255, 36/255];
lagVals = -4:4;
aniDistsTrial = cell(size(encodes));
aniDistsRaw = cell(size(encodes));
aniDists = cell(1,size(encodes,2));
aniLagDists = nan(size(encodes,2), length(lagVals), size(encodes,1));
for a = 1:size(encodes,2)
    figure;
    annotation(gcf,'textbox', [0 0.93 1 0.05],'String', sprintf('%s InSeq Lag (%s)', aniID{a}, alignment),...
        'FontWeight','bold', 'FontSize',12, 'edgecolor', 'none', 'horizontalalignment', 'left');
    sp = nan(size(encodes,1),3);
    tempAniMeanDist = nan(max(encodes{1}(:,5))+1,max(encodes{1}(:,5))+1,size(encodes,1));
    for t = 1:size(encodes,1)
        sp(1,t) = subplot(3, size(encodes,1), sub2ind([size(encodes,1),3], t,1));
        hold(sp(1,t), 'on');
        title(window(t));
        tempEncodes = encodes{t,a};
        posEncodes = cell(1,max(tempEncodes(:,5))+1);
        trialIDs = cell(max(tempEncodes(:,5))+1,1);
        for pos = min(tempEncodes(:,5))+1:max(tempEncodes(:,5))+1
            posLog = tempEncodes(:,4)==(pos-1);
            isLog = tempEncodes(:,6)==1;
            perfLog = tempEncodes(:,7)==1;
            posEncodes{pos} = tempEncodes(posLog & isLog & perfLog, 8:9);
            trialIDs{pos} = tempEncodes(posLog & isLog & perfLog, 3);
            scatter(sp(1,t), posEncodes{pos}(:,1), posEncodes{pos}(:,2), 5,'marker', 'o', 'markerfacecolor', odorColors(pos,:), 'markeredgecolor', 'none');
        end
        tempRawDists = cell(max(tempEncodes(:,5))+1);
        meanDists = nan(max(tempEncodes(:,5))+1);
        for pos1 = min(tempEncodes(:,5))+1:max(tempEncodes(:,5))+1
            for pos2 = min(tempEncodes(:,5))+1:max(tempEncodes(:,5))+1
                if strcmp(distCalc, 'Euclid')
                    tempDist = squareform(pdist([posEncodes{pos1};posEncodes{pos2}]));
                    distLog = [false(size(posEncodes{pos1},1),size(posEncodes{pos1},1)), true(size(posEncodes{pos1},1),size(posEncodes{pos2},1));...
                        false(size(posEncodes{pos2},1), size(posEncodes{pos1},1)), false(size(posEncodes{pos2},1), size(posEncodes{pos2},1))];
                    tempRawDists{pos1,pos2} = nanmean(reshape(tempDist(distLog), [size(posEncodes{pos1},1), size(posEncodes{pos2},1)]),2);
                elseif strcmp(distCalc, 'Mahal')
                    tempRawDists{pos1,pos2} = mahal(posEncodes{pos1}, posEncodes{pos2});
                end
                meanDists(pos1,pos2) = nanmean(tempRawDists{pos1,pos2});
            end
        end
        aniDistsRaw{t,a} = tempRawDists;
        tempTrialDists = cell(size(tempRawDists));
        for pos2 = min(tempEncodes(:,5))+1:max(tempEncodes(:,5))+1
            for pos1 = min(tempEncodes(:,5))+1:max(tempEncodes(:,5))+1
                uniqIDs = unique(trialIDs{pos1});
                tempTrlMean = nan(size(uniqIDs));
                currRaw = tempRawDists{pos1,pos2};
                for trl = 1:length(uniqIDs)
                    tempTrlMean(trl) = mean(currRaw(trialIDs{pos1}==uniqIDs(trl)));
                end
                tempTrialDists{pos1, pos2} = tempTrlMean;
            end
        end     
        aniDistsTrial{t,a} = tempTrialDists;
        tempAniMeanDist(:,:,t) = meanDists;
        sp(2,t) = subplot(3, size(encodes,1), sub2ind([size(encodes,1),3],t,2));
        imagesc(meanDists);
%         set(sp(2,t), 'ydir', 'normal');
        
        sp(3,t) = subplot(3, size(encodes,1), sub2ind([size(encodes,1),3],t,3));
        lagVect = nan(1,length(lagVals));
        for lag = 1:length(lagVals)
            lagVect(lag) = nanmean(meanDists(logical(triu(ones(max(tempEncodes(:,5))+1), lagVals(lag))) & logical(tril(ones(max(tempEncodes(:,5))+1), lagVals(lag)))));
        end
        lagVect = (lagVect./lagVect(5))-1;
        aniLagDists(a,:,t) = lagVect;
        bar(lagVals, lagVect);
    end
    aniDists{a} = tempAniMeanDist;
    drawnow;
end

tempRawData = reshape(aniDistsRaw, [size(aniDistsRaw,1),1,size(aniDistsRaw,2)]);
tempTrialData = reshape(aniDistsTrial, [size(aniDistsTrial,1),1,size(aniDistsTrial,2)]);
rawDists = cell(size(encodes,1),1);     
trialDists = cell(size(encodes,1),1);
for t = 1:size(encodes,1)
    tempRaw = reshape([tempRawData{t,:,:}],[max(tempEncodes(:,5))+1,max(tempEncodes(:,5))+1,size(encodes,2)]);
    tempTrial = reshape([tempTrialData{t,:,:}],[max(tempEncodes(:,5))+1,max(tempEncodes(:,5))+1,size(encodes,2)]);
    tempRawDists = cell(max(tempEncodes(:,5))+1);
    tempTrialDists = cell(max(tempEncodes(:,5))+1);
    for r = min(tempEncodes(:,5)+1):max(tempEncodes(:,5))+1
        for c = min(tempEncodes(:,5))+1:max(tempEncodes(:,5))+1
            temptempRaw = cellfun(@(a)a', tempRaw(r,c,:), 'uniformoutput', 0);
            temptempTrial = cellfun(@(a)a', tempTrial(r,c,:), 'uniformoutput', 0);
            tempRawDists{r,c} = [temptempRaw{:}];
            tempTrialDists{r,c} = [temptempTrial{:}];
        end
    end
    rawDists{t} = tempRawDists;
    trialDists{t} = tempTrialDists;
end
        
trialLagDists = cell(length(trialDists), length(lagVect));
rawLagDists = cell(length(rawDists), length(lagVect));
for t = 1:length(trialDists)
    for lag = 1:length(lagVals)
        tempRawDists = rawDists{t};
        tempTrialDists = trialDists{t};        
        rawLagDists{t,lag} = [tempRawDists{logical(triu(ones(max(tempEncodes(:,5))+1), lagVals(lag))) & logical(tril(ones(max(tempEncodes(:,5))+1), lagVals(lag)))}];
        trialLagDists{t,lag} = [tempTrialDists{logical(triu(ones(max(tempEncodes(:,5))+1), lagVals(lag))) & logical(tril(ones(max(tempEncodes(:,5))+1), lagVals(lag)))}];
    end
end
end

%% CalculateOutSeqLagDist
function [posAniDistsTrialLag, odrAniDistsTrialLag, aniDistsPOS, aniDistsODR, aniLagDistsPOS, aniLagDistsODR, posAniDistsRawLag, odrAniDistsRawLag] = CalculateOutSeqLagDist(encodes, distCalc)
lagVals = -4:4;
aniDistsTrialPOS = cell(size(encodes));
aniDistsTrialODR = cell(size(encodes));
aniDistsRawPOS = cell(size(encodes));
aniDistsRawODR = cell(size(encodes));
aniDistsPOS = cell(1,size(encodes,2));
aniDistsODR = cell(1,size(encodes,2));
aniLagDistsPOS = nan(size(encodes,2), length(lagVals), size(encodes,1));
aniLagDistsODR = nan(size(encodes,2), length(lagVals), size(encodes,1));

for a = 1:size(encodes,2)
    tempAniDistPOS = nan(max(encodes{1}(:,5))+1,max(encodes{1}(:,5))+1,size(encodes,1));
    tempAniDistODR = nan(max(encodes{1}(:,5))+1,max(encodes{1}(:,5))+1,size(encodes,1));
    for t = 1:size(encodes,1)
        tempEncodes = encodes{t,a};
        posEncodes = cell(1,max(tempEncodes(:,5))+1);
        for pos = 1:5
            posLog = tempEncodes(:,4)==(pos-1);
            isLog = tempEncodes(:,6)==1;
            perfLog = tempEncodes(:,7)==1;
            posEncodes{pos} = tempEncodes(posLog & isLog & perfLog,8:9);
        end
        outSeqEncodes = tempEncodes(tempEncodes(:,6)==0 & tempEncodes(:,7)==1,:);
        osTrialIDs = cell(max(tempEncodes(:,5))+1);
        tempDistsTrialPOS = cell(max(tempEncodes(:,5))+1);
        tempDistsTrialODR = cell(max(tempEncodes(:,5))+1);
        rawDistsPOS = cell(max(tempEncodes(:,5))+1);
        rawDistsODR = cell(max(tempEncodes(:,5))+1);
        meanDistsPOS = nan(max(tempEncodes(:,5))+1);
        meanDistsODR = nan(max(tempEncodes(:,5))+1);
        for posX = min(tempEncodes(:,5))+1:max(tempEncodes(:,5))+1
            for odrY = min(tempEncodes(:,5))+1:max(tempEncodes(:,5))+1
                osVals = outSeqEncodes(outSeqEncodes(:,4)+1==posX & outSeqEncodes(:,5)+1==odrY,8:9);
                osTrialIDs{odrY, posX} = outSeqEncodes(outSeqEncodes(:,4)+1==posX & outSeqEncodes(:,5)+1==odrY,3);
                if ~isempty(osVals)
                    if posX~=odrY
                        pX = posEncodes{posX};
                        oY = posEncodes{odrY};
                    else
                        pX = osVals;
                        oY = osVals;
                    end
                    if strcmp(distCalc, 'Euclid')
                        tempDistPOS = squareform(pdist([osVals; pX]));
                        tempDistODR = squareform(pdist([osVals; oY]));
                        distLogPOS = [false(size(osVals,1),size(osVals,1)), true(size(osVals,1), size(pX,1));...
                            false(size(pX,1), size(osVals,1)), false(size(pX,1), size(pX,1))];
                        distLogODR = [false(size(osVals,1),size(osVals,1)), true(size(osVals,1), size(oY,1));...
                            false(size(oY,1), size(osVals,1)), false(size(oY,1), size(oY,1))];
                        rawDistsPOS{odrY,posX} = nanmean(reshape(tempDistPOS(distLogPOS), [size(osVals,1), size(pX,1)]),2)';
                        rawDistsODR{odrY,posX} = nanmean(reshape(tempDistODR(distLogODR), [size(osVals,1), size(oY,1)]),2)';
                    elseif strcmp(distCalc, 'Mahal')
                        rawDistsPOS{odrY,posX} = mahal(osVals, pX)';
                        rawDistsODR{odrY,posX} = mahal(osVals, oY)';
                    end
                    % Calculate Mean
                    meanDistsPOS(odrY,posX) = nanmean(rawDistsPOS{odrY,posX});
                    meanDistsODR(odrY,posX) = nanmean(rawDistsODR{odrY,posX});
                end
            end
        end
        aniDistsRawPOS{t,a} = rawDistsPOS;
        aniDistsRawODR{t,a} = rawDistsODR;
        for y = 1:size(tempDistsTrialODR,1)
            for x = 1:size(tempDistsTrialPOS,1)
                uniqIDs = unique(osTrialIDs{y,x});
                tempTrlMeanPOS = nan(size(uniqIDs));
                tempTrlMeanODR = nan(size(uniqIDs));
                currRawPOS = rawDistsPOS{y,x};
                currRawODR = rawDistsODR{y,x};
                for trl = 1:length(uniqIDs)
                    tempTrlMeanPOS(trl) = mean(currRawPOS(osTrialIDs{y,x}==uniqIDs(trl)));
                    tempTrlMeanODR(trl) = mean(currRawODR(osTrialIDs{y,x}==uniqIDs(trl)));
                end
                tempDistsTrialPOS{y, x} = tempTrlMeanPOS';
                tempDistsTrialODR{y, x} = tempTrlMeanODR';
            end
        end     
        aniDistsTrialPOS{t,a} = tempDistsTrialPOS;
        aniDistsTrialODR{t,a} = tempDistsTrialODR;
        tempAniDistPOS(:,:,t) = meanDistsPOS;
        tempAniDistODR(:,:,t) = meanDistsODR;
        
        lagVectPOS = nan(1,length(lagVals));
        lagVectODR = nan(1,length(lagVals));
        for lag = 1:length(lagVals)
            lagVectPOS(lag) = nanmean(meanDistsPOS(logical(triu(ones(max(tempEncodes(:,5))+1), lagVals(lag))) & logical(tril(ones(max(tempEncodes(:,5))+1), lagVals(lag)))));
            lagVectODR(lag) = nanmean(meanDistsODR(logical(triu(ones(max(tempEncodes(:,5))+1), lagVals(lag))) & logical(tril(ones(max(tempEncodes(:,5))+1), lagVals(lag)))));
        end
        aniLagDistsPOS(a,:,t) = lagVectPOS;
        aniLagDistsODR(a,:,t) = lagVectODR;
    end
    aniDistsPOS{a} = tempAniDistPOS;
    aniDistsODR{a} = tempAniDistODR;
end
posAniDistsTrial = cell(size(aniDistsTrialPOS,1),1);
odrAniDistsTrial = cell(size(aniDistsTrialODR,1),1);
posAniDistsRaw = cell(size(aniDistsRawPOS,1),1);
odrAniDistsRaw = cell(size(aniDistsRawODR,1),1);
for t = 1:size(aniDistsTrialPOS,1)
    tempTrialDistsPOS = reshape([aniDistsTrialPOS{t,:}], [size(aniDistsTrialPOS{1},1),size(aniDistsTrialPOS{1},2),size(aniDistsTrialPOS,2)]);
    tempTrialDistsODR = reshape([aniDistsTrialODR{t,:}], [size(aniDistsTrialODR{1},1),size(aniDistsTrialODR{1},2),size(aniDistsTrialODR,2)]);
    tempRawDistsPOS = reshape([aniDistsRawPOS{t,:}], [size(aniDistsRawPOS{1},1), size(aniDistsRawPOS{1},2), size(aniDistsRawPOS,2)]);
    tempRawDistsODR = reshape([aniDistsRawODR{t,:}], [size(aniDistsRawODR{1},1), size(aniDistsRawODR{1},2), size(aniDistsRawODR,2)]);
    tempTempTrialDistsPOS = cell(5);
    tempTempTrialDistsODR = cell(5);
    tempTempRawDistsPOS = cell(5);
    tempTempRawDistsODR = cell(5);
    for r = 1:size(aniDistsTrialPOS{1},1)
        for c = 1:size(aniDistsTrialPOS{1},2)
            tempTempTrialDistsPOS{r,c} = [tempTrialDistsPOS{r,c,:}];
            tempTempTrialDistsODR{r,c} = [tempTrialDistsODR{r,c,:}];
            tempTempRawDistsPOS{r,c} = [tempRawDistsPOS{r,c,:}];
            tempTempRawDistsODR{r,c} = [tempRawDistsODR{r,c,:}];
        end
    end
    posAniDistsTrial{t} = tempTempTrialDistsPOS;
    odrAniDistsTrial{t} = tempTempTrialDistsODR;
    posAniDistsRaw{t} = tempTempRawDistsPOS;
    odrAniDistsRaw{t} = tempTempRawDistsODR;
end
posAniDistsTrialLag = cell(length(posAniDistsTrial), length(lagVals));
odrAniDistsTrialLag = cell(length(odrAniDistsTrial), length(lagVals));
posAniDistsRawLag = cell(length(posAniDistsRaw), length(lagVals));
odrAniDistsRawLag = cell(length(odrAniDistsRaw), length(lagVals));
for t = 1:length(posAniDistsTrial)
    for lag = 1:length(lagVals)
        posAniDistsTrialLag{t,lag} = [posAniDistsTrial{t}{logical(triu(ones(5), lagVals(lag))) & logical(tril(ones(5), lagVals(lag)))}];
        odrAniDistsTrialLag{t,lag} = [odrAniDistsTrial{t}{logical(triu(ones(5), lagVals(lag))) & logical(tril(ones(5), lagVals(lag)))}];
        posAniDistsRawLag{t,lag} = [posAniDistsRaw{t}{logical(triu(ones(5), lagVals(lag))) & logical(tril(ones(5), lagVals(lag)))}];
        odrAniDistsRawLag{t,lag} = [odrAniDistsRaw{t}{logical(triu(ones(5), lagVals(lag))) & logical(tril(ones(5), lagVals(lag)))}];
    end
end
end

%% CalculateOutSeqTrialDist
function [aniOSC, aniOSI,...
    aniOSCtrialDistPOS, aniOSItrialDistPOS,...
    aniOSCtrialDistODR, aniOSItrialDistODR,...
    aniOSCtrialLagDistPOS, aniOSItrialLagDistPOS,...
    aniOSCtrialLagDistODR, aniOSItrialLagDistODR,...
    aniOSCrawLagDistPOS, aniOSIrawLagDistPOS,...
    aniOSCrawLagDistODR, aniOSIrawLagDistODR] = CalculateOutSeqTrialDist(encodes, distCalc)
lagVect = -4:4;
aniOSC = cell(size(encodes,1),1);
aniOSI = cell(size(encodes,1),1);
aniOSCtrialDistPOS = cell(size(encodes,1),1);
aniOSItrialDistPOS = cell(size(encodes,1),1);
aniOSCtrialDistODR = cell(size(encodes,1),1);
aniOSItrialDistODR = cell(size(encodes,1),1);
aniOSCtrialLagDistPOS = cell(size(encodes,1),length(lagVect));  
aniOSItrialLagDistPOS = cell(size(encodes,1),length(lagVect));
aniOSCtrialLagDistODR = cell(size(encodes,1),length(lagVect));
aniOSItrialLagDistODR = cell(size(encodes,1),length(lagVect));
    aniOSCrawLagDistPOS = cell(size(encodes,1),length(lagVect));
    aniOSIrawLagDistPOS = cell(size(encodes,1),length(lagVect));
    aniOSCrawLagDistODR = cell(size(encodes,1),length(lagVect));
    aniOSIrawLagDistODR = cell(size(encodes,1),length(lagVect));
for t = 1:size(encodes,1)
    tempOSCraw = cell(5,5,size(encodes,2));
    tempOSCtrialDistsPOS = cell(5,5,size(encodes,2));
    tempOSCtrialDistsODR = cell(5,5,size(encodes,2));
    tempOSCtrialLagDistPOS = cell(length(lagVect),size(encodes,2));
    tempOSCrawLagDistPOS = cell(length(lagVect),size(encodes,2));
    tempOSCtrialLagDistODR = cell(size(encodes,2),length(lagVect));
    tempOSCrawLagDistODR = cell(size(encodes,2),length(lagVect));

    tempOSIraw = cell(5,5,size(encodes,2));
    tempOSItrialDistsPOS = cell(5,5,size(encodes,2));
    tempOSItrialDistsODR = cell(5,5,size(encodes,2));
    tempOSItrialLagDistPOS = cell(length(lagVect),size(encodes,2));
    tempOSIrawLagDistPOS = cell(length(lagVect),size(encodes,2));
    tempOSItrialLagDistODR = cell(size(encodes,2),length(lagVect));
    tempOSIrawLagDistODR = cell(size(encodes,2),length(lagVect));
    for a = 1:size(encodes,2)
        tempEncodes = encodes{t,a};
        posEncodes = cell(1,5);
        for pos = 1:5
            posLog = tempEncodes(:,4)==(pos-1);
            isLog = tempEncodes(:,6)==1;
            perfLog = tempEncodes(:,7)==1;
            posEncodes{pos} = tempEncodes(posLog & isLog & perfLog,8:9);
        end
        outSeqEncodes = tempEncodes(tempEncodes(:,6)==0,:);
        unqTrials = unique(outSeqEncodes(:,3));
        for trl = 1:length(unqTrials)
            curTrlEncodes = outSeqEncodes(outSeqEncodes(:,3)==unqTrials(trl),:);
            tempPerf = curTrlEncodes(1,7)==1;
            tempPos = curTrlEncodes(1,4)+1;
            tempOdr = curTrlEncodes(1,5)+1;
            if strcmp(distCalc, 'Euclid')
                tempDistPOS = squareform(pdist([curTrlEncodes(:,8:9); posEncodes{tempPos}]));
                tempDistODR = squareform(pdist([curTrlEncodes(:,8:9); posEncodes{tempOdr}]));
                distLogPOS = [false(size(curTrlEncodes,1),size(curTrlEncodes,1)), true(size(curTrlEncodes,1), size(posEncodes{tempPos},1));...
                    false(size(posEncodes{tempPos},1), size(curTrlEncodes,1)), false(size(posEncodes{tempPos},1), size(posEncodes{tempPos},1))];
                distLogODR = [false(size(curTrlEncodes,1),size(curTrlEncodes,1)), true(size(curTrlEncodes,1), size(posEncodes{tempOdr},1));...
                    false(size(posEncodes{tempOdr},1), size(curTrlEncodes,1)), false(size(posEncodes{tempOdr},1), size(posEncodes{tempOdr},1))];
                tempPosDist = nanmean(reshape(tempDistPOS(distLogPOS), [size(curTrlEncodes,1), size(posEncodes{tempPos},1)]),2);
                tempOdrDist = nanmean(reshape(tempDistODR(distLogODR), [size(curTrlEncodes,1), size(posEncodes{tempOdr},1)]),2);
            elseif strcmp(distCalc, 'Mahal')
                tempPosDist = mahal(curTrlEncodes(:,8:9), posEncodes{tempPos});
                tempOdrDist = mahal(curTrlEncodes(:,8:9), posEncodes{tempOdr});
            end
            if tempPerf==1
                tempOSCraw{tempOdr, tempPos, a} = [tempOSCraw{tempOdr, tempPos, a}; curTrlEncodes(:,8:9)];
                tempOSCtrialDistsPOS{tempOdr, tempPos, a} = [tempOSCtrialDistsPOS{tempOdr, tempPos, a}, mean(tempPosDist)];
                tempOSCtrialDistsODR{tempOdr, tempPos, a} = [tempOSCtrialDistsODR{tempOdr, tempPos, a}, mean(tempOdrDist)];
            elseif tempPerf==0
                tempOSIraw{tempOdr, tempPos, a} = [tempOSIraw{tempOdr, tempPos, a}; curTrlEncodes(:,8:9)];
                tempOSItrialDistsPOS{tempOdr, tempPos, a} = [tempOSItrialDistsPOS{tempOdr, tempPos, a}, mean(tempPosDist)];
                tempOSItrialDistsODR{tempOdr, tempPos, 2} = [tempOSItrialDistsODR{tempOdr, tempPos, 2}, mean(tempOdrDist)];
            end
        end
        [tempOSCtrialLagDistPOS(:,a), tempOSCrawLagDistPOS(:,a)] = TabluateOutSeqLagPOS(tempOSCraw(:,:,a), lagVect, distCalc);
        [tempOSCtrialLagDistODR(a,:), tempOSCrawLagDistODR(a,:)] = TabluateOutSeqLagODR(tempOSCraw(:,:,a), lagVect, distCalc);
        
        [tempOSItrialLagDistPOS(:,a), tempOSIrawLagDistPOS(:,a)] = TabluateOutSeqLagPOS(tempOSIraw(:,:,a), lagVect, distCalc);
        [tempOSItrialLagDistODR(a,:), tempOSIrawLagDistODR(a,:)] = TabluateOutSeqLagODR(tempOSIraw(:,:,a), lagVect, distCalc);
    end
    aniOSC{t} = tempOSCraw;
    aniOSI{t} = tempOSIraw;
    aniOSCtrialDistPOS{t} = tempOSCtrialDistsPOS;
    aniOSItrialDistPOS{t} = tempOSItrialDistsPOS;
    aniOSCtrialDistODR{t} = tempOSCtrialDistsODR;
    aniOSItrialDistODR{t} = tempOSItrialDistsODR;
    
    for lag = 1:length(lagVect)
        aniOSCtrialLagDistPOS{t,lag} = [tempOSCtrialLagDistPOS{lag,:}];
        aniOSItrialLagDistPOS{t,lag} = [tempOSItrialLagDistPOS{lag,:}];  
        aniOSCtrialLagDistODR{t,lag} = [tempOSCtrialLagDistODR{:,lag}];
        aniOSItrialLagDistODR{t,lag} = [tempOSItrialLagDistODR{:,lag}];
        
        aniOSCrawLagDistPOS{t,lag} = [tempOSCrawLagDistPOS{lag,:}];
        aniOSIrawLagDistPOS{t,lag} = [tempOSIrawLagDistPOS{lag,:}];  
        aniOSCrawLagDistODR{t,lag} = [tempOSCrawLagDistODR{:,lag}];
        aniOSIrawLagDistODR{t,lag} = [tempOSIrawLagDistODR{:,lag}];
    end
end
end

%% Misc Functions
function [trialLagOS, rawLagOS] = TabluateOutSeqLagPOS(osLagMtx, lagVect, distCalc)
% trialLagOS = cell(length(lagVect), 5);
trialLagOS = cell(length(lagVect), 1);
rawLagOS = cell(length(lagVect), 1);
for pos = 1:5
    for o1 = 1:5
        if ~isempty(osLagMtx{o1,pos})
            tempRef = osLagMtx{o1,pos};
            for o2 = 1:5
                if ~isempty(osLagMtx{o2,pos})
                    if strcmp(distCalc, 'Euclid')
                        temp = squareform(pdist([osLagMtx{o2,pos}; tempRef]));
                        distLog = [false(size(osLagMtx{o2,pos},1),size(osLagMtx{o2,pos},1)), true(size(osLagMtx{o2,pos},1), size(tempRef,1));...
                            false(size(tempRef,1), size(osLagMtx{o2,pos},1)), false(size(tempRef,1), size(tempRef,1))];
                        tempDist = nanmean(reshape(temp(distLog), [size(osLagMtx{o2,pos},1), size(tempRef,1)]),2);
                    elseif strcmp(distCalc, 'Mahal')
                        tempDist = mahal(osLagMtx{o2,pos}, tempRef);
                    end
%                     trialLagOS{lagVect==(o2-o1),pos} = [trialLagOS{lagVect==(o2-o1),pos}, mean(tempDist)];
                    trialLagOS{lagVect==(o2-o1)} = [trialLagOS{lagVect==(o2-o1)}, mean(tempDist)];
                    rawLagOS{lagVect==(o2-o1)} = [rawLagOS{lagVect==(o2-o1)}, tempDist'];
                end
            end
        end
    end
end
end

function [trialLagOS, rawLagOS] = TabluateOutSeqLagODR(osLagMtx, lagVect, distCalc)
% trialLagOS = cell(5,length(lagVect));
trialLagOS = cell(1,length(lagVect));
rawLagOS = cell(1,length(lagVect));
for odr = 1:5
    for p1 = 1:5
        if ~isempty(osLagMtx{odr,p1})
            tempRef = osLagMtx{odr,p1};
            for p2 = 1:5
                if ~isempty(osLagMtx{odr,p2})
                    if strcmp(distCalc, 'Euclid')
                        temp = squareform(pdist([osLagMtx{odr,p2}; tempRef]));
                        distLog = [false(size(osLagMtx{odr,p2},1),size(osLagMtx{odr,p2},1)), true(size(osLagMtx{odr,p2},1), size(tempRef,1));...
                            false(size(tempRef,1), size(osLagMtx{odr,p2},1)), false(size(tempRef,1), size(tempRef,1))];
                        tempDist = nanmean(reshape(temp(distLog), [size(osLagMtx{odr,p2},1), size(tempRef,1)]),2);
                    elseif strcmp(distCalc, 'Mahal')
                        tempDist = mahal(osLagMtx{odr,p2}, tempRef);
                    end
%                     trialLagOS{odr,lagVect==(p2-p1)} = [trialLagOS{odr,lagVect==(p2-p1)}; mean(tempDist)];
                    trialLagOS{lagVect==(p2-p1)} = [trialLagOS{lagVect==(p2-p1)}, mean(tempDist)];
                    rawLagOS{lagVect==(p2-p1)} = [rawLagOS{lagVect==(p2-p1)}, tempDist'];
                end
            end
        end
    end
end
end

%%
function PlotOSscatters(anis, encodes, outSeqMat, window, odorColors, tit)
for a = 1:length(anis)
    figure;
    annotation(gcf,'textbox', [0 0.93 1 0.05],'String', sprintf('%s %s', anis{a}, tit),...
        'FontWeight','bold', 'FontSize',12, 'edgecolor', 'none', 'horizontalalignment', 'left');
    for t = 1:size(outSeqMat,1)
        spA(t) = subplot(2, size(outSeqMat,1), t);
        tempEncodes = encodes{t,a};
%         scatter(tempEncodes(tempEncodes(:,6)==1 & tempEncodes(:,7)==1,8), tempEncodes(tempEncodes(:,6)==1 & tempEncodes(:,7)==1,9), 15,...
%             'markerfacecolor', 'k', 'markerfacealpha', 0.2, 'marker', 'o', 'markeredgecolor', 'none');
%         hold on;
        tempOSC = outSeqMat{t}(:,:,a);
        for o = 1:size(tempOSC,1)
            for p = 1:size(tempOSC,2)
                if ~isempty(tempOSC{o,p})
                    scatter(tempOSC{o,p}(:,1), tempOSC{o,p}(:,2), 25,...
                        'markerfacecolor', odorColors(o,:), 'markeredgecolor', odorColors(p,:),...
                        'marker', 'o', 'markerfacealpha', 0.5,'markeredgealpha', 0.25, 'linewidth', 2);
                    hold on;
                end
            end
        end
        for o = 1:size(tempOSC,1)
            for p = 1:size(tempOSC,2)
                if ~isempty(tempOSC{o,p})
                    scatter(mean(tempOSC{o,p}(:,1)), mean(tempOSC{o,p}(:,2)), 150,...
                        'markerfacecolor', odorColors(o,:), 'markeredgecolor', odorColors(p,:),...
                        'marker', 'o', 'linewidth', 3);
                end
            end
        end
        %         for o = 1:size(tempOSC,1)
        %             tempOSodor = cell2mat(tempOSC(o,:)');
        %             if ~isempty(tempOSodor)
        %                 scatter(mean(tempOSodor(:,1)), mean(tempOSodor(:,2)), 100,...
        %                     'markerfacecolor', odorColors(o,:), 'markeredgecolor', 'k',...
        %                     'marker', 'o');
        %             end
        %         end
        if t == 1
            ylabel('Odor Labels', 'fontweight', 'bold');
        end
        title(window(t));
        spB(t) = subplot(2, size(outSeqMat,1), t+size(outSeqMat,1));
%         scatter(tempEncodes(tempEncodes(:,6)==1 & tempEncodes(:,7)==1,8), tempEncodes(tempEncodes(:,6)==1 & tempEncodes(:,7)==1,9), 15,...
%             'markerfacecolor', 'k', 'markerfacealpha', 0.2, 'marker', 'o', 'markeredgecolor', 'none');
%         hold on;
        for p = 1:size(tempOSC,2)
            for o = 1:size(tempOSC,1)
                if ~isempty(tempOSC{o,p})
                    scatter(tempOSC{o,p}(:,1), tempOSC{o,p}(:,2), 15,...
                        'markerfacecolor', odorColors(p,:), 'markeredgecolor', 'none',...
                        'marker', 'o', 'markerfacealpha', 0.75);
                    hold on;
                end
            end
        end
%         for p = 1:size(tempOSC,2)
%             for o = 1:size(tempOSC,1)
%                 if ~isempty(tempOSC{o,p})
%                     scatter(mean(tempOSC{o,p}(:,1)), mean(tempOSC{o,p}(:,2)), 100,...
%                         'markerfacecolor', odorColors(o,:), 'markeredgecolor', odorColors(p,:),...
%                         'marker', 'o', 'linewidth', 3);
%                 end
%             end
%         end
        %         for p = 1:size(tempOSC,2)
        %             tempOSpos = cell2mat(tempOSC(:,p));
        %             if ~isempty(tempOSpos)
        %                 scatter(mean(tempOSpos(:,1)), mean(tempOSpos(:,2)), 100,...
        %                     'markerfacecolor', odorColors(p,:), 'markeredgecolor', 'k',...
        %                     'marker', 'o');
        %             end
        %         end
        lagDist = (size(tempOSC,1)-1)*-1:size(tempOSC,1);
        for lag = 1:length(lagDist)
            tempOSlag = cell2mat(tempOSC(logical(tril(ones(size(tempOSC,1)), lagDist(lag))) & logical(triu(ones(size(tempOSC,1)), lagDist(lag)))));
            if ~isempty(tempOSlag)
                scatter(mean(tempOSlag(:,1)), mean(tempOSlag(:,2)), 100,...
                    'markerfacecolor', 'k', 'markeredgecolor', 'k',...
                    'marker', 'o', 'markerfacealpha', 1/lag);
            end
        end
        if t == 1
            ylabel('Position Labels', 'fontweight', 'bold');
        end
    end
%     linkaxes([spA, spB], 'xy');
end
end

%% 

function ExportMatCSVlag(mat2xport, timeBin, fName)
for t = 1:size(mat2xport,1)
    xportMat = nan(max(cellfun(@(a)length(a),mat2xport(t,:))),size(mat2xport,2));
    for lag = 1:size(mat2xport,2)
        xportMat(1:length(mat2xport{t,lag}),lag) = mat2xport{t,lag};
    end
    writematrix(xportMat, sprintf('%s\\Dists\\%s_%.02f.csv',cd, fName, timeBin(t)));
%     fprintf('%s_%.02f.csv saved\n', fName, timeBin(t));
end
end

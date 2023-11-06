clear
close all

%% load data

%CA1 data
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA1_SlopeDistrib.mat')
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA1perAnimal_slopes_props.mat')

% CA1N
group{1} = ones(size(CA1slopes.N));  
slopes{1} = CA1slopes.N; % all PFs pooled
CA1N.Bck = 0.4850; CA1N.Fwd = 0.0454; CA1N.ns = 0.4696; % based on all PFs pooled (computed in "CanDataAnalysis")
meanSlopeCA1N_animals = mean(meanSlope_CA1(1,:),2);
meanASlopeCA1N_animals = mean(meanASlope_CA1(1,:),2);
PrBck{1} = mean(PropBck_CA1(1,:)); % average across animals
PrFwd{1} = mean(PropFwd_CA1(1,:));
PrNon{1} = mean(PropNonSig_CA1(1,:));


% BTSP CA1N-like
modelN = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_Speed15_Pamp20pA_pCS0_005_500repeats.mat');
% modelN = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_005.mat');
group{2} = 2*ones(length(modelN.COMslope),1);
slopes{2} = modelN.COMslope';
PrBck{2} = modelN.PropBck;
PrFwd{2} = modelN.PropFwd;
PrNon{2} = modelN.PropNonSig;
% minSlope(1) = 

% CA1F
group{3} = 3*ones(size(CA1slopes.F));  
slopes{3} = CA1slopes.F;
CA1F.Bck = 0.3315; CA1F.Fwd = 0.0880; CA1F.ns = 0.5805; % based on all PFs pooled
meanSlopeCA1F_animals = mean(meanSlope_CA1(2,:),2); % average across animals
meanASlopeCA1F_animals = mean(meanASlope_CA1(2,:),2);
PrBck{3} = mean(PropBck_CA1(2,:));
PrFwd{3} = mean(PropFwd_CA1(2,:));
PrNon{3} = mean(PropNonSig_CA1(2,:));

% BTSP CA1F-like
modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_Speed15_Pamp20pA_pCS0_002_500repeats.mat');
% modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_0025.mat');
group{4} = 4*ones(length(modelF.COMslope),1);
slopes{4} = modelF.COMslope';
PrBck{4} = modelF.PropBck;
PrFwd{4} = modelF.PropFwd;
PrNon{4} = modelF.PropNonSig;

%% stats
for i = 1:length(group)
    meanSlope(i) = mean(slopes{i});
    SlopeCI = bootci(1000, @mean, slopes{i});
    SlopeCIlo(i) = SlopeCI(1,1);
    SlopeCIup(i) = SlopeCI(2,1);

    medianASlope(i) = median(abs(slopes{i}));
    ASlopeCI = bootci(1000, @median, abs(slopes{i}));
    ASlopeCIlo(i) = ASlopeCI(1,1);
    ASlopeCIup(i) = ASlopeCI(2,1);
end

%% Load all models for comparison

% BTSP Pamp 15pA
model{1}{1} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp15pA\Workspace_BTSP100repeats_Speed15_Pamp15pA_pCS0_001.mat');
model{1}{2} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp15pA\Workspace_BTSP100repeats_Speed15_Pamp15pA_pCS0_002.mat');
model{1}{3} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp15pA\Workspace_BTSP100repeats_Speed15_Pamp15pA_pCS0_003.mat');
model{1}{4} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp15pA\Workspace_BTSP100repeats_Speed15_Pamp15pA_pCS0_004.mat');
model{1}{5} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp15pA\Workspace_BTSP100repeats_Speed15_Pamp15pA_pCS0_005.mat');
model{1}{6} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp15pA\Workspace_BTSP100repeats_Speed15_Pamp15pA_pCS0_01.mat');
model{1}{7} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp15pA\Workspace_BTSP100repeats_Speed15_Pamp15pA_pCS0_02.mat');

% BTSP Pamp 20pA
model{2}{1} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_BTSP100repeats_Speed15_Pamp20pA_pCS0_001.mat');
model{2}{2} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_BTSP100repeats_Speed15_Pamp20pA_pCS0_002.mat');
model{2}{3} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_BTSP100repeats_Speed15_Pamp20pA_pCS0_003.mat');
model{2}{4} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_BTSP100repeats_Speed15_Pamp20pA_pCS0_004.mat');
model{2}{5} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_BTSP100repeats_Speed15_Pamp20pA_pCS0_005.mat');
model{2}{6} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_BTSP100repeats_Speed15_Pamp20pA_pCS0_01.mat');
model{2}{7} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_BTSP100repeats_Speed15_Pamp20pA_pCS0_02.mat');

% BTSP Pamp 25.5pA
model{3}{1} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_001.mat');
model{3}{2} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_002.mat');
model{3}{3} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_0025.mat');
model{3}{4} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_003.mat');
model{3}{5} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_004.mat');
model{3}{6} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_0045.mat');
model{3}{7} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_005.mat');
model{3}{8} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_01.mat');
model{3}{9} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_02.mat');

% BTSP Pamp 30pA
model{4}{1} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp30pA\Workspace_BTSP100repeats_Speed15_Pamp30pA_pCS0_001.mat');
model{4}{2} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp30pA\Workspace_BTSP100repeats_Speed15_Pamp30pA_pCS0_002.mat');
model{4}{3} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp30pA\Workspace_BTSP100repeats_Speed15_Pamp30pA_pCS0_003.mat');
model{4}{4} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp30pA\Workspace_BTSP100repeats_Speed15_Pamp30pA_pCS0_004.mat');
model{4}{5} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp30pA\Workspace_BTSP100repeats_Speed15_Pamp30pA_pCS0_005.mat');
model{4}{6} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp30pA\Workspace_BTSP100repeats_Speed15_Pamp30pA_pCS0_01.mat');
model{4}{7} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp30pA\Workspace_BTSP100repeats_Speed15_Pamp30pA_pCS0_02.mat');

% BTSP Pamp 80pA
model{5}{1} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp80pA\Workspace_BTSP100repeats_Speed15_Pamp80pA_pCS0_001.mat');
model{5}{2} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp80pA\Workspace_BTSP100repeats_Speed15_Pamp80pA_pCS0_002.mat');
model{5}{3} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp80pA\Workspace_BTSP100repeats_Speed15_Pamp80pA_pCS0_003.mat');
model{5}{4} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp80pA\Workspace_BTSP100repeats_Speed15_Pamp80pA_pCS0_004.mat');
model{5}{5} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp80pA\Workspace_BTSP100repeats_Speed15_Pamp80pA_pCS0_005.mat');
model{5}{6} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp80pA\Workspace_BTSP100repeats_Speed15_Pamp80pA_pCS0_01.mat');
model{5}{7} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp80pA\Workspace_BTSP100repeats_Speed15_Pamp80pA_pCS0_02.mat');

pCS{1} = [0.001 0.002, 0.003, 0.004, 0.005, 0.01, 0.02];
pCS{2} = pCS{1};
pCS{3} = [0.001 0.002, 0.0025, 0.003, 0.004, 0.0045, 0.005, 0.01, 0.02];
pCS{4} = pCS{1};
pCS{5} = pCS{1};

for m = 1:length(model)
    for n = 1:length(model{m})
mdlMaxSlope{m}(n) = max(model{m}{n}.COMslope);
mdlMinSlope{m}(n) = min(model{m}{n}.COMslope);
mdlMeanSlope{m}(n) = mean(model{m}{n}.COMslope);
mdlSlopeCI = bootci(1000, @mean, model{m}{n}.COMslope);
mdlSlopeCIlo{m}(n) = mdlSlopeCI(1,1);
mdlSlopeCIup{m}(n) = mdlSlopeCI(2,1);
mdlPrBck{m}(n) = model{m}{n}.PropBck;
mdlPrFwd{m}(n) = model{m}{n}.PropFwd;
mdlPrNon{m}(n) = model{m}{n}.PropNonSig;
    end
end

%% summary plots for all models

% fwd PFs should be red, Bck PFs blue, NonSig grey or black
blue = [0 0.45 0.74]; % backward
red = [0.85 0.32 0.01]; % forward
grey = [0.5 0.5 0.5]; %nonsig

% specific Pamp

mc = 2;

figure % Proportions as a function of p(CS)
plot(pCS{mc}, mdlPrBck{mc}, 'o-', 'Color', blue, 'LineWidth', 1.5); hold on
plot(pCS{mc}, mdlPrFwd{mc}, 'o-', 'Color', red, 'LineWidth', 1.5); hold on 
plot(pCS{mc}, mdlPrNon{mc}, 'o-', 'Color', grey, 'LineWidth', 1.5); hold on 
yline(PrFwd{3}, '--', 'Color', red, 'LineWidth', 1);
yline(PrBck{3}, '--', 'Color', blue, 'LineWidth', 1);
yline(PrNon{3}, '--', 'Color', grey, 'LineWidth', 1);
legend('B', 'F', 'N', 'CA1F', 'CA1F', 'CA1F', 'Location', 'BestOutside')
ylim([0 1]); xlim([0 max(pCS{mc})])
xlabel('p(CS)')
ylabel('proportion of PFs')
title('speed: 15cm/s, Pamp: 20pA')
box off
axis square

figure % max, min and mean slopes
plot(pCS{mc}, mdlMinSlope{mc}, 'o-', 'Color', blue, 'LineWidth', 1.5); hold on
plot(pCS{mc}, mdlMaxSlope{mc}, 'o-', 'Color', red, 'LineWidth', 1.5); hold on
y = mdlMeanSlope{mc}; yneg = y-mdlSlopeCIlo{mc}; ypos = mdlSlopeCIup{mc}-y;
errorbar(pCS{mc}, y,yneg,ypos,'o-', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3, 'LineWidth', 1.5); hold on
legend('min', 'max', 'mean')
xlim([0 max(pCS{mc})])
xlabel('p(CS)')
ylabel('Shifting slopes (cm/lap)')
title('speed: 15cm/s, Pamp: 20pA')
box off
axis square

clear y yneg ypos

% All Pamps
cmap_Bck = brewermap(length(model)+1,'Blues'); % [0 0.7 1; 0 0.45 0.74; 0 0.2 0.49; 0 0 0.3];  
cmap_Fwd = brewermap(length(model)+1,'Reds');
cmap_Non = flipud(gray(length(model)+1));

figure % Proportions as a function of p(CS) and Pamp
for m = 1:length(model)
plot(pCS{m}, mdlPrBck{m}, '-', 'Color', cmap_Bck(m+1,:), 'LineWidth', 1); hold on
plot(pCS{m}, mdlPrFwd{m}, '-', 'Color', cmap_Fwd(m+1,:), 'LineWidth', 1); hold on 
plot(pCS{m}, mdlPrNon{m}, '-', 'Color', cmap_Non(m+1,:), 'LineWidth', 1); hold on 
end
% yline(PrFwd{3}, '--', 'Color', red, 'LineWidth', 1);
% yline(PrBck{3}, '--', 'Color', blue, 'LineWidth', 1);
% yline(PrNon{3}, '--', 'Color', grey, 'LineWidth', 1);
legend('B 15pA', 'F 15pA', 'N 15pA', 'B 20pA', 'F 20pA', 'N 20pA', 'B 25.5pA', 'F 25.5pA', 'N 25.5pA', 'B 30pA', 'F 30pA', 'N 30pA', 'B 80pA', 'F 80pA', 'N 80pA', 'Location', 'BestOutside')
ylim([0 1]); 
xlabel('p(CS)')
ylabel('proportion of PFs')
title('speed: 15cm/s, Pamp: 15 to 80pA')
box off
axis square

figure % max, min and mean slopes
for m = 1:length(model)
plot(pCS{m}, mdlMinSlope{m}, '-', 'Color', cmap_Bck(m+1,:), 'LineWidth', 1); hold on
plot(pCS{m}, mdlMaxSlope{m}, '-', 'Color', cmap_Fwd(m+1,:), 'LineWidth', 1); hold on
plot(pCS{m}, mdlMeanSlope{m}, '-', 'Color', cmap_Non(m+1,:), 'LineWidth', 1);
% y = mdlMeanSlope; yneg = y-mdlSlopeCIlo; ypos = mdlSlopeCIup-y;
% errorbar(pCS, y,yneg,ypos,'o-', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3, 'LineWidth', 1.5); hold on
end
legend('min', 'max', 'mean')
% xlim([0 max(pCS)])
xlabel('p(CS)')
ylabel('Shifting slopes (cm/lap)')
title('speed: 15cm/s, Pamp: 15 to 80pA')
box off
axis square

figure
nc = 4;
for m = 1:length(model)
y005(m) = mdlMeanSlope{m}(nc); 
y005neg(m) = y005(m)-mdlSlopeCIlo{m}(nc); y005pos(m) = mdlSlopeCIup{m}(nc)-y005(m);
end
errorbar([15, 20, 25.5, 30, 80], y005,y005neg,y005pos,'o', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3, 'LineWidth', 1.5); hold on
xlabel('Params.Pamp (pA)')
ylabel('Mean shifting slope (cm/lap)')
title('speed: 15cm/s, pCS=0.005')
box off
axis square
%% Comparison CA1 vs best models

gp = cat(1, group{:});
COMslopeAll = cat(1, slopes{:});
PropBck = cat(1, PrBck{:});
PropFwd = cat(1, PrFwd{:});
PropNon = cat(1, PrNon{:});
Props = [PropBck PropFwd PropNon];
PropsAll = Props; % PFs pooled from all animals
PropsAll(1,:) = [CA1N.Bck CA1N.Fwd CA1N.ns];
PropsAll(3,:) = [CA1F.Bck CA1F.Fwd CA1F.ns];

%colorcode CA1 vs CA3, N (dark) vs F (light); e.g. CA1F -> mapCA1(1,:)
map = brewermap(8,'Paired'); %cool(4); %brewermap(4,'PiYG')
mapCA1 = map(3:4,:); %brewermap(2,'GnBu');
mapCA3 = map(7:8,:); %brewermap(2,'OrRd');
mapCA = [mapCA1;mapCA3];
colors = [mapCA1(2,:); mapCA1(1,:); mapCA3(2,:); mapCA3(1,:)];
% figure; scatter(1:4,1:4,200,mapCA, 'filled')

figure 
subplot(1,3,1) % violin plots. 
    v = violinplot(COMslopeAll, gp, 'ViolinColor', [0 0 0], 'ShowData', false, 'ShowMean', true);
    v(1).ViolinColor = mapCA1(2,:);
    v(3).ViolinColor = mapCA1(1,:);
    yline(0, 'k--')
    ylabel('Shifting slopes (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA1_N'; '0.005'; 'CA1_F'; '0.002']);
    box off;
    axis square
subplot(1,3,2) % bootstrapped means compared to CA1 mean (averaged across animals)
    y = meanSlope; yneg = y-SlopeCIlo; ypos = SlopeCIup-y;
    % errorbar(1, y(1),yneg(1),ypos(1),'o', 'Color', mapCA1(2,:), 'MarkerFaceColor', mapCA1(2,:), 'MarkerSize', 3); hold on
    errorbar(1:length(group), y,yneg,ypos,'o', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3); hold on
    yline(0, 'k-');
    yline(meanSlopeCA1N_animals, '--', 'Color', mapCA1(2,:), 'LineWidth', 1.5);
    yline(meanSlopeCA1F_animals, '--', 'Color', mapCA1(1,:), 'LineWidth', 1.5);
    xlim([0.5 length(group)+0.5]);
    ylim([-0.5 0.5]); 
    xlabel('models'); ylabel('bootstrapped mean slope (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA1_N'; '0.005'; 'CA1_F'; '0.002']);
    box off;
    axis square
subplot(1,3,3) % props averaged across animals
    b = bar(Props, 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylabel('proportion of PFs'); xlabel('models')
    title('averaged across animals')
    xlim([0.5 length(group)+0.5]);
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA1_N'; '0.005'; 'CA1_F'; '0.002']);
    box off;
    axis square

    figure 
subplot(2,2,1) % violin plots. 
    v = violinplot(COMslopeAll, gp, 'ViolinColor', [0 0 0], 'ShowData', false, 'ShowMean', true);
    v(1).ViolinColor = mapCA1(2,:);
    v(3).ViolinColor = mapCA1(1,:);
    yline(0, 'k--')
    ylabel('Shifting slopes (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA1_N'; '0.005'; 'CA1_F'; '0.002']);
    box off;
    axis square
subplot(2,2,2) % bootstrapped means compared to CA1 mean (pooled PFs)
    y = meanSlope; yneg = y-SlopeCIlo; ypos = SlopeCIup-y;
    errorbar(1:length(group), y,yneg,ypos,'o', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3); hold on
    yline(0, 'k-');
    yline(meanSlope(1), '--', 'Color', mapCA1(2,:), 'LineWidth', 1.5);
    yline(meanSlope(3), '--', 'Color', mapCA1(1,:), 'LineWidth', 1.5);
    xlim([0.5 length(group)+0.5]);
    ylim([-0.5 0.5]); 
    xlabel('models'); ylabel('mean slope (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA1_N'; '0.005'; 'CA1_F'; '0.002']);
    box off;
    axis square
subplot(2,2,3) % median abs(Slopes) compared to CA1 median (pooled PFs)
    y = medianASlope; yneg = y-ASlopeCIlo; ypos = ASlopeCIup-y;
    errorbar(1:length(group), y,yneg,ypos,'o', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3); hold on
    yline(0, 'k-');
    yline(medianASlope(1), '--', 'Color', mapCA1(2,:), 'LineWidth', 1.5);
    yline(medianASlope(3), '--', 'Color', mapCA1(1,:), 'LineWidth', 1.5);
    xlim([0.5 length(group)+0.5]);
    ylim([0 0.5]); 
    xlabel('models'); ylabel('median |slope| (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA1_N'; '0.005'; 'CA1_F'; '0.002']);
    box off;
    axis square
subplot(2,2,4) % props 
    b2 = bar(PropsAll, 'stacked', 'FaceColor','flat');
    b2(1).CData = blue; b2(2).CData = red; b2(3).CData = grey;
    ylabel('proportion of PFs'); xlabel('models')
    xlim([0.5 length(group)+0.5]);
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA1_N'; '0.005'; 'CA1_F'; '0.002']);
    title('all PFs pooled')
    box off;
    axis square
clear
close all

%% load data

%CA1 data
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA1_SlopeDistrib.mat')
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA1perAnimal_slopes_props.mat')
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\\CanData_CA1_MSD.mat')

% CA1N
group{1} = ones(size(CA1slopes.N));  
slopes{1} = CA1slopes.N; % all PFs pooled
CA1N.Bck = 0.4850; CA1N.Fwd = 0.0454; CA1N.ns = 0.4696; % based on all PFs pooled 15 laps minimum (computed in "CanDataAnalysis")
% CA1N.Bck = 0.49; CA1N.Fwd = 0.05; CA1N.ns = 0.46; % based on all PFs pooled 30 laps minimum (computed in "CanDataAnalysis")
meanSlopeCA1N_animals = mean(meanSlope_CA1(1,:),2);
meanASlopeCA1N_animals = mean(meanASlope_CA1(1,:),2);
PrBck{1} = mean(PropBck_CA1(1,:)); % average across animals
PrFwd{1} = mean(PropFwd_CA1(1,:));
PrNon{1} = mean(PropNonSig_CA1(1,:));


% BTSP CA1N-like
modelN = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPhomeo_pCSdynamic\N-like tries\Workspace_instantBTSPhomeo_Pamp80pA_pCSdyn_N-like_Apcs0_65_Bpcs_0_025_Taupcs1_3_SDcs29_Bound60_WUdyntau15s_Ddecay0_2_Pdecay1_5_500samples.mat');
% modelN = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPhomeo_pCSdynamic\N-like tries\Workspace_instantBTSPhomeo_Pamp80pA_pCSdyn_N-like_Apcs0_8_Bpcs_0_03_Taupcs1_1_SDcs25_Bound60_WUdyntau15s_Ddecay0_2_Pdecay1_5_500samples.mat');
% modelN = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\MilsteinModel\Wmax125pA\workspace_pCS0_3%_repeat1rule_500samples_Wmax125pA.mat');
% modelN = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_Speed15_Pamp20pA_pCS0_005_500repeats.mat');
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
CA1F.Bck = 0.3315; CA1F.Fwd = 0.0880; CA1F.ns = 0.5805; % based on all PFs pooled 15 minlaps
% CA1F.Bck = 0.3375; CA1F.Fwd = 0.0955; CA1F.ns = 0.5673; % based on all PFs pooled 30 minlaps
meanSlopeCA1F_animals = mean(meanSlope_CA1(2,:),2); % average across animals
meanASlopeCA1F_animals = mean(meanASlope_CA1(2,:),2);
PrBck{3} = mean(PropBck_CA1(2,:));
PrFwd{3} = mean(PropFwd_CA1(2,:));
PrNon{3} = mean(PropNonSig_CA1(2,:));

% BTSP CA1F-like
modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPhomeo_pCSdynamic\Workspace_instantBTSPhomeo_Pamp80pA_pCSdyn_F-like_Apcs0_6_Bpcs_0_0125_Taupcs1_SDcs21_Bound60_WUdyntau15s_Ddecay0_2_Pdecay1_5_500samples.mat');
%modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPhomeo_pCSdynamic\Workspace_instantBTSPhomeo_Pamp80pA_pCSdynamic_F-like_Apcs0_53_Bpcs_0_025_Taupcs1.35_SDcs20_Bound60_WUdyntau15s_500samples.mat');
% modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\MilsteinModel\Wmax125pA\workspace_pCS0_2%_repeat1rule_500samples_Wmax125pA.mat');
% modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_Speed15_Pamp20pA_pCS0_002_500repeats.mat');
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

%%
% Apcs = [0.8, 0.6, 0.7, 0.5]; % frequency of CSs on second lap. Cf Supp Fig associated to Fig 7. Same max amplitude for Familiar (F) and Novel (N) conditions.
% Tau_pcs = [1.5, 1, 1.25, 0.75]; % in laps (cf Fig 6 and 7 of our manuscript.). Same time constants for F and N (but exp fits on instantaneous MSD suggest 1.15 for N and 0.85 for F)
% Bpcs = [0.05, 0.03, 0.04, 0.02]; % 0 for F, > 0 for N (cf Fig 7: a static pCS is not enough to lead to high Diffusion coeff in later laps)
% % SDcs = [29, 21];
% 
% laps = 0:0.1:28;
% % x=1:300;
% for n = 1:4
% pCSlap(n,:) = Apcs(n).*exp(-laps./Tau_pcs(n)) + Bpcs(n); % lapwise p(CS) not tied to current PF
% % pCSloc(n,:) = (1/(SDcs(n)*sqrt(2*pi))).*exp(-0.5*(x - 300/2).^2/SDcs(n)^2);
% end
% 
% figure
% for n = 1:4
% plot(laps+1, pCSlap(n,:), '-', 'LineWidth', 1,'Color', colors(n,:)); hold on
% end
% yline(0,'k--')
% ylim([-0.03 1])
% xlim([0 15])
% xlabel('post-onset laps')
% ylabel('lap-wise p(CS)')
% box off
% axis square

%%
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
subplot(1,3,2) % bootstrapped means compared to CA1 means averaged across animals)
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


figure 
subplot(1,3,1) % violin plots. 
    v = violinplot(COMslopeAll, gp, 'ViolinColor', [0 0 0], 'ShowData', false, 'ShowMean', true);
    v(1).ViolinColor = mapCA1(2,:);
    v(3).ViolinColor = mapCA1(1,:);
    yline(0, 'k--')
    ylabel('Shifting slopes (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA1N'; '0.3%'; 'CA1F'; '0.2%']);
    box off;
    axis square
subplot(1,3,2) % bootstrapped means compared to CA1 mean (averaged across animals)
    y = medianASlope; yneg = y-ASlopeCIlo; ypos = ASlopeCIup-y;
    errorbar(1:length(group), y,yneg,ypos,'o', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3); hold on
    yline(0, 'k-');
    yline(medianASlope(1), '--', 'Color', mapCA1(2,:), 'LineWidth', 1.5);
    yline(medianASlope(3), '--', 'Color', mapCA1(1,:), 'LineWidth', 1.5);
    xlim([0.5 length(group)+0.5]);
    ylim([0 0.5]); 
    xlabel('models'); ylabel('median |slope| (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA1N'; '0.3%'; 'CA1F'; '0.2%']);
    box off;
    axis square
subplot(1,3,3) % props 
    b2 = bar(PropsAll, 'stacked', 'FaceColor','flat');
    b2(1).CData = blue; b2(2).CData = red; b2(3).CData = grey;
    ylabel('proportion of PFs'); xlabel('models')
    xlim([0.5 length(group)+0.5]);
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA1N'; '0.3%'; 'CA1F'; '0.2%']);
    title('all PFs pooled')
    box off;
    axis square

%% P(CS) dynamic models

TrajsN = cat(1,modelN.COMbin{:}).*6; % in cm
Disp = TrajsN - TrajsN(:,1); % trajectories centered on initial location, with forward direction as positive values
DispSq = Disp.^2;
MSD = mean(DispSq,1,"omitmissing"); % average all neurons for each lap
for lap = 1:size(TrajsN,2)
msdCI(:,lap) = bootci(1000, @nanmean, DispSq(:,lap));
end

TrajsF = cat(1,modelF.COMbin{:}).*6; % in cm
DispF = TrajsF - TrajsF(:,1); % trajectories centered on initial location, with forward direction as positive values
DispFSq = DispF.^2;
MSDf = mean(DispFSq,1,"omitmissing"); % average all neurons for each lap
for lap = 1:size(TrajsF,2)
msdFCI(:,lap) = bootci(1000, @nanmean, DispFSq(:,lap));
end

RegStart = 4;
Xlm = [ [RegStart:length(MSD)]', ones(size([RegStart:length(MSD)]'))];
[Bmsd,BINTmsd,Rmsd,RINTmsd,STATSmsd] = regress(MSD(RegStart:end)', Xlm);
lmMSD = Xlm*Bmsd;
[BmsdF,BINTmsdF,RmsdF,RINTmsdF,STATSmsdF] = regress(MSDf(RegStart:end)', Xlm);
lmMSDf = Xlm*BmsdF;

% BTSP models color code
% cline = lines(7);
% mapBTSP(1,:) = cline(1,:); % N
% mapBTSP(2,:) = cline(1,:)+0.25;% F
Blues = brewermap(9,'*Blues');
mapBTSP(1,:) = Blues(1,:); % N
mapBTSP(2,:) = Blues(4,:); % F

% parameters for [N,F]
Apcs = [0.65, 0.6]; % frequency of CSs on second lap. Cf Supp Fig associated to Fig 7. Same max amplitude for Familiar (F) and Novel (N) conditions.
Tau_pcs = [1.35, 1]; % in laps (cf Fig 6 and 7 of our manuscript.). Same time constants for F and N (but exp fits on instantaneous MSD suggest 1.15 for N and 0.85 for F)
Bpcs = [0.025, 0.0125]; % 0 for F, > 0 for N (cf Fig 7: a static pCS is not enough to lead to high Diffusion coeff in later laps)
SDcs = [29, 21];

laps = 0:0.1:28;
x=1:300;
for n = 1:2
pCSlap(n,:) = Apcs(n).*exp(-laps./Tau_pcs(n)) + Bpcs(n); % lapwise p(CS) not tied to current PF
pCSloc(n,:) = (1/(SDcs(n)*sqrt(2*pi))).*exp(-0.5*(x - 300/2).^2/SDcs(n)^2);
end

figure
tiledlayout('flow', 'TileSpacing','Compact','Padding','Compact')

nexttile % lapwise p(CS)
for n = 1:2
plot(laps+1, pCSlap(n,:), '-', 'LineWidth', 1,'Color', mapBTSP(n,:)); hold on
end
legend('BTSPn', 'BTSPf')
ylim([0 0.7])
xlabel('post-onset laps')
ylabel('lap-wise p(CS)')
box off
axis square

nexttile % position-wise p(CS)
% PF = (1/(SDcs(2)*sqrt(2*pi))).*exp(-0.5*(x - 300/2).^2/13^2); hold on
colororder({'b','k'})
PF = 10.*exp(-0.5*(x - 300/2).^2/13^2); hold on
yyaxis right
plot(x, PF, 'Color', [0.5 0.5 0.5],'LineWidth', 1); hold on
ylabel('Hz')
set(gca, 'YTick', [0, 10], 'YTickLabel', {'0', '10'})
yyaxis left
for n = 1:2
    plot(x, pCSloc(n,:), '-', 'LineWidth', 1,'Color', mapBTSP(n,:)); hold on
end
% xline(150,'k--')
xline([150-60, 150+60], 'k--')
set(gca, 'XTick', [0, 90, 150, 210, 300], 'XTickLabel', {'-150', '-60', '0', '+60', '+150'})
set(gca, 'YTick', [0, round(max(pCSloc,[],"all"),1, 'significant')])
xlabel('COM-centered position on track')
ylabel('position-wise p(CS)')
box off; axis square

nexttile % MSD
plot_ci([1:length(MSD)]-1, [MSD' msdCI(1,:)' msdCI(2,:)'], 'MainLineColor', mapBTSP(1,:), 'MainLineWidth', 1.5, 'PatchColor', mapBTSP(1,:), 'LineColor', mapBTSP(1,:), 'PatchAlpha', 0.6); hold on
plot_ci([1:length(MSDf)]-1, [MSDf' msdFCI(1,:)' msdFCI(2,:)'], 'MainLineColor', mapBTSP(2,:), 'MainLineWidth', 1.5, 'PatchColor', mapBTSP(2,:), 'LineColor', mapBTSP(2,:), 'PatchAlpha', 0.7); hold on
for n = 1:2
    plot([1:length(MSDc(:,n))]-1, MSDc(:,n),'-', 'LineWidth',1,'Color', colors(n,:)); hold on 
end
plot([1:length(MSD)]-1, MSD','-', 'LineWidth',1,'Color', mapBTSP(1,:)); hold on 
plot([1:length(MSDf)]-1, MSDf','-', 'LineWidth',1,'Color', mapBTSP(2,:)); hold on
text(length(MSD), max(lmMSD), ['D = ' num2str( round(Bmsd(1)/2, 2, 'significant') )], 'Color', mapBTSP(1,:)); hold on
text(length(MSD), max(lmMSDf), ['D = ' num2str( round(BmsdF(1)/2, 2, 'significant') )], 'Color', mapBTSP(2,:))
xlabel('Laps after onset')
ylabel('mean squared displacement (cm^2)')
box off
axis square

nexttile % violin plots. 
    v = violinplot(COMslopeAll, gp, 'ViolinColor', [0 0 0], 'ShowData', false, 'ShowMean', true, 'ViolinAlpha', 0.5);
    v(1).ViolinColor = mapCA1(2,:);
    v(2).ViolinColor = mapBTSP(1,:);
    v(3).ViolinColor = mapCA1(1,:);
    v(4).ViolinColor = mapBTSP(2,:);
    yline(0, 'k--')
    ylabel('Shifting slopes (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', {'CA1N'; 'BTSPn'; 'CA1F'; 'BTSPf'});
    box off;
    axis square

nexttile
    y = meanSlope; yneg = y-SlopeCIlo; ypos = SlopeCIup-y;
    errorbar(1:length(group), y,yneg,ypos,'o', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3); hold on
    yline(0, 'k-');
    yline(meanSlope(1), '--', 'Color', mapCA1(2,:), 'LineWidth', 1.5);
    yline(meanSlope(3), '--', 'Color', mapCA1(1,:), 'LineWidth', 1.5);
    xlim([0.5 length(group)+0.5]);
    ylim([-0.5 0.5]); 
    ylabel('mean slope (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', {'CA1N'; 'BTSPn'; 'CA1F'; 'BTSPf'});
    box off;
    axis square

nexttile % props 
    b2 = bar(PropsAll, 'stacked', 'FaceColor','flat');
    b2(1).CData = blue; b2(2).CData = red; b2(3).CData = grey;
    ylabel('fraction of PFs');
    xlim([0.5 length(group)+0.5]);
    set(gca, 'XTick', 1:length(group), 'XTickLabel', {'CA1N'; 'BTSPn'; 'CA1F'; 'BTSPf'});
    % title('all PFs pooled')
    box off;
    axis square


 figure
for n = 1:2
plot(laps+1, pCSlap(n,:), '-', 'LineWidth', 1,'Color', mapBTSP(n,:)); hold on
yline(Bpcs(n),'--', 'LineWidth', 1, 'Color', mapBTSP(n,:))
end
legend('BTSPn', 'BTSPf')
ylim([0 0.1])
xlim([0 15])
xlabel('post-onset laps')
ylabel('lap-wise p(CS)')
box off
axis square
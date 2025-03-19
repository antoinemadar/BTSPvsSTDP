    
clear
close all

%% load data

%CA1 data
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA1_SlopeDistrib.mat')
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA1perAnimal_slopes_props.mat')
group{1} = ones(size(CA1slopes.N));  
slopes{1} = CA1slopes.N;
meanSlopeCA1N_animals = mean(meanSlope_CA1(1,:),2);
meanASlopeCA1N_animals = mean(meanASlope_CA1(1,:),2);
meanSlopeCA1F_animals = mean(meanSlope_CA1(2,:),2);
meanASlopeCA1F_animals = mean(meanASlope_CA1(2,:),2);
% PrBck{1} = mean(PropBck_CA1(2,:),2);
% PrFwd{1} = mean(PropFwd_CA1(2,:),2);
% PrNon{1} = mean(PropNonSig_CA1(2,:),2);
PrBck{1} = 0.4850; PrFwd{1} = 0.0454;  PrNon{1} = 0.4696; % based on all PFs pooled (computed in "CanDataAnalysis")

% classic STDP - baseline input params (average realistic output): PFsd18cm_PFamp10Hz_Wsd10cm
modelA = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\STDP\100samples_RealisticLowFRmap_SongAbbottParams\W_capped\ClassicSTDP20msAmp0.5%\25HzOutputPeak\PFsd18cm_PFamp10Hz_Wsd10cm\Simul3\Workspace.mat');
group{2} = 2*ones(length(modelA.COMslope),1);
slopes{2} = modelA.COMslope';
PrBck{2} = modelA.PropBck;
PrFwd{2} = modelA.PropFwd;
PrNon{2} = modelA.PropNonSig;

% % classic STDP - realistic output, different input params: PFsd10cm_PFamp15Hz_Wsd13cm
% modelB= load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\STDP\100samples_RealisticLowFRmap_SongAbbottParams\W_capped\ClassicSTDP20msAmp0.5%\25HzOutputPeak\PFsd10cm_PFamp15Hz_Wsd13cm\workspace.mat');
% group{3} = 3*ones(length(modelB.COMslope),1);
% slopes{3} = modelB.COMslope';
% PrBck{3} = modelB.PropBck;
% PrFwd{3} = modelB.PropFwd;
% PrNon{3} = modelB.PropNonSig;

% classic STDP - borderline realistic output, with high output FR: PFsd18cm_PFamp15Hz_Wsd10cm
modelB= load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\STDP\100samples_RealisticLowFRmap_SongAbbottParams\W_capped\ClassicSTDP20msAmp0.5%\Input15Hz_45HzMeanPFOut\workspace.mat');
group{3} = 3*ones(length(modelB.COMslope),1);
slopes{3} = modelB.COMslope';
PrBck{3} = modelB.PropBck;
PrFwd{3} = modelB.PropFwd;
PrNon{3} = modelB.PropNonSig;

% model A + dynamic CA3N-like inputs
% modelC = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\STDP\100samples_RealisticLowFRmap_SongAbbottParams\Wnobound_WUdyn_Adapt\HighOutputFR_PFamp15Hz_NarrowWsd\workspace.mat');
modelC = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\STDP\100samples_RealisticLowFRmap_SongAbbottParams\W_capped\ClassicSTDP20msAmp0.5%\25HzOutputPeak\DynamicInputs\CA3Ninputs\workspace.mat');
group{4} = 4*ones(length(modelC.COMslope),1);
slopes{4} = modelC.COMslope';
PrBck{4} = modelC.PropBck;
PrFwd{4} = modelC.PropFwd;
PrNon{4} = modelC.PropNonSig;

% Model B + dynamic CA3_N-like input
 % CA3Fin: 'G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\STDP\100samples_RealisticLowFRmap_SongAbbottParams\W_capped\ClassicSTDP20msAmp0.5%\25HzOutputPeak\DynamicInputs\CA3Finputs\workspace.mat');
modelD = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\STDP\100samples_RealisticLowFRmap_SongAbbottParams\W_capped\ClassicSTDP20msAmp0.5%\Input15Hz_45HzMeanPFOut\withCA3NlikeDynamicInputs\workspace.mat');
group{5} = 5*ones(length(modelD.COMslope),1);
slopes{5} = modelD.COMslope';
PrBck{5} = modelD.PropBck;
PrFwd{5} = modelD.PropFwd;
PrNon{5} = modelD.PropNonSig;

% classic STDP, static inputs, high output FR (


% STDP capped, decay 20ms, high amplitude (5%)

% STDP capped, decay 20ms, unrealistically high amplitude (30%)

% uncapped very high output FR


% 
%% stats
for i = 1:length(group)
    meanSlope(i) = mean(slopes{i});
    SlopeCI = bootci(1000, @mean, slopes{i});
    SlopeCIlo(i) = SlopeCI(1,1);
    SlopeCIup(i) = SlopeCI(2,1);
end


%% plots

gp = cat(1, group{:});
COMslopeAll = cat(1, slopes{:});
PropBck = cat(1, PrBck{:});
PropFwd = cat(1, PrFwd{:});
PropNon = cat(1, PrNon{:});

%colorcode CA1 vs CA3, N (dark) vs F (light); e.g. CA1F -> mapCA1(1,:)
map = brewermap(8,'Paired'); %cool(4); %brewermap(4,'PiYG')
mapCA1 = map(3:4,:); %brewermap(2,'GnBu');
mapCA3 = map(7:8,:); %brewermap(2,'OrRd');
mapCA = [mapCA1;mapCA3];
colors = [mapCA1(2,:); mapCA1(1,:); mapCA3(2,:); mapCA3(1,:)];
% figure; scatter(1:4,1:4,200,mapCA, 'filled')

% fwd PFs should be red, Bck PFs blue, NonSig grey or black
blue = [0 0.45 0.74]; % backward
red = [0.85 0.32 0.01]; % forward
grey = [0.5 0.5 0.5]; %nonsig


figure % violin plots. 

subplot(3,1,1)
v = violinplot(COMslopeAll, gp, 'ViolinColor', [0 0 0], 'ShowData', false, 'ShowMean', true);
v(1).ViolinColor = mapCA1(2,:);
yline(0, 'k--')
ylabel({'COM slopes'; '(cm/lap)'});
set(gca, 'XTick', [], 'XTickLabel', []);
box off;

% figure % bootstrapped means compared to CA1 mean
subplot(3,1,2)
y = meanSlope; yneg = y-SlopeCIlo; ypos = SlopeCIup-y;
errorbar(1, y(1),yneg(1),ypos(1),'o', 'Color', mapCA1(2,:), 'MarkerFaceColor', mapCA1(2,:), 'MarkerSize', 3); hold on
errorbar(2:length(group), y(2:end),yneg(2:end),ypos(2:end),'o', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3); hold on
yline(0, 'k-');
% yline(meanSlopeCA1N_animals, '--', 'Color', mapCA1(2,:));
yline(y(1), '--', 'Color', mapCA1(2,:))
xlim([0.5 length(group)+0.5]);
ylim([-0.5 0.15]); 
ylabel({'bootstrapped'; 'mean slope'; '(cm/lap)'});
set(gca, 'XTick', [], 'XTickLabel', []);
box off;

% props 
% figure
subplot(3,1,3)
Props = [PropBck PropFwd PropNon];
b = bar(Props, 'stacked', 'FaceColor','flat');
b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
ylabel('proportion of PFs'); %xlabel('models')
xlim([0.5 length(group)+0.5]);
set(gca, 'XTick', 1:length(group), 'XTickLabel', {'CA1N data'; 'model A'; 'model B'; 'mdlA + CA3N'; 'mdlB + CA3N'});
box off;
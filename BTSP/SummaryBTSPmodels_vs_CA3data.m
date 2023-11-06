clear
close all

%% load data

%CA3 data
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA3_SlopeDistrib_All.mat')

% CA3N
group{1} = ones(size(CA3slopes.N));  
slopes{1} = CA3slopes.N; % all PFs pooled
CA3N.Bck = 0.3423; CA3N.Fwd = 0.0225; CA3N.ns = 0.6351; % based on all PFs pooled (computed in "CanDataAnalysis")
PrBck{1} = CA3N.Bck; PrFwd{1} = CA3N.Fwd; PrNon{1} = CA3N.ns;
% meanSlopeCA1N_animals = mean(meanSlope_CA1(1,:),2);
% meanASlopeCA1N_animals = mean(meanASlope_CA1(1,:),2);
% PrBck{1} = mean(PropBck_CA1(1,:)); % average across animals
% PrFwd{1} = mean(PropFwd_CA1(1,:));
% PrNon{1} = mean(PropNonSig_CA1(1,:));


% BTSP CA3N-like
modelN = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0017_1sprepost_0.02spostpre.mat');
% modelN = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0015_1sprepost_1.3spostpre.mat');
group{2} = 2*ones(length(modelN.COMslope),1);
slopes{2} = modelN.COMslope';
PrBck{2} = modelN.PropBck;
PrFwd{2} = modelN.PropFwd;
PrNon{2} = modelN.PropNonSig;
% minSlope(1) = 

% CA3F
group{3} = 3*ones(size(CA3slopes.F));  
slopes{3} = CA3slopes.F;
CA3F.Bck = 0.09; CA3F.Fwd = 0.19; CA3F.ns = 0.72; % based on all PFs pooled
PrBck{3} = CA3F.Bck; PrFwd{3} = CA3F.Fwd; PrNon{3} = CA3F.ns;
% meanSlopeCA1F_animals = mean(meanSlope_CA1(2,:),2); % average across animals
% meanASlopeCA1F_animals = mean(meanASlope_CA1(2,:),2);
% PrBck{3} = mean(PropBck_CA1(2,:));
% PrFwd{3} = mean(PropFwd_CA1(2,:));
% PrNon{3} = mean(PropNonSig_CA1(2,:));

% BTSP CA3F-like
modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0015_1sprepost_1.3spostpre.mat');
% modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0015_0.2sprepost_0.5spostpre.mat');
% modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0015_SymmetricBTSP.mat');
% modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.001_SymmetricBTSP.mat');
% modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp25_5pA\Workspace_BTSP100repeats_Speed15_Pamp25_5pA_pCS0_0025.mat');
% modelF = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_Speed15_Pamp20pA_pCS0_002_500repeats.mat');
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

% BTSP pCS = 0.0015, different ta_prepost and tau_postpre combinations
model{1}{1} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0015_1sprepost_0.02spostpre.mat');
model{1}{2} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0015_1sprepost_0.5spostpre.mat');
model{1}{3} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0015_1sprepost_1spostpre.mat');
model{1}{4} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0015_1sprepost_1.3spostpre.mat');

% pCS 0.15% tau_prepost = 1.31s
model{2}{1} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0015_1.31sprepost_0.69spostpre.mat');
model{2}{2} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.0015_SymmetricBTSP.mat');

% BTSP pCS = 0.002, different ta_prepost and tau_postpre combinations
model{3}{1} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.002_1sprepost_0.02spostpre.mat');
model{3}{2} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.002_1.3sprepost_0.4spostpre.mat');
model{3}{3} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\CA3rule\BTSP_SheffieldParams_pCS0.002_1sprepost_0.2spostpre.mat');
model{3}{4} = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_Speed15_Pamp20pA_pCS0_002_500repeats.mat');

tauP{1} = [1, 1, 1, 1];
tauD{1} = [0.02, 0.5, 1, 1.3];
tauP{2} = [1.31, 1.31];
tauD{2} = [0.69, 1.31];
tauP{3} = [1, 1.3, 1, 1.31];
tauD{3} = [0.02, 0.4, 0.2, 0.69];


for m = 1:length(model)
Asym{m} = tauD{m} - tauP{m};
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

%% color code

%colorcode CA1 vs CA3, N (dark) vs F (light); e.g. CA1F -> mapCA1(1,:)
% fwd PFs should be red, Bck PFs blue, NonSig grey or black
blue = [0 0.45 0.74]; % backward
red = [0.85 0.32 0.01]; % forward
grey = [0.5 0.5 0.5]; %nonsig
% 
map = brewermap(8,'Paired'); %cool(4); %brewermap(4,'PiYG')
mapCA1 = map(3:4,:); %brewermap(2,'GnBu');
mapCA3 = map(7:8,:); %brewermap(2,'OrRd');
mapCA = [mapCA1;mapCA3];
colors = [mapCA1(2,:); mapCA1(1,:); mapCA3(2,:); mapCA3(1,:)];
% figure; scatter(1:4,1:4,200,mapCA, 'filled')

% All Pamps
cmap_Bck = brewermap(length(model)+1,'Blues'); % [0 0.7 1; 0 0.45 0.74; 0 0.2 0.49; 0 0 0.3];  
cmap_Fwd = brewermap(length(model)+1,'Reds');
cmap_Non = flipud(gray(length(model)+1));
%% All models, asymmetry variation for different pCS

figure % Proportions as a function of p(CS) and Tau Asymetry
for m = 1:length(model)
plot(Asym{m}, mdlPrBck{m}, '-', 'Color', cmap_Bck(m+1,:), 'LineWidth', 1); hold on
plot(Asym{m}, mdlPrFwd{m}, '-', 'Color', cmap_Fwd(m+1,:), 'LineWidth', 1); hold on 
plot(Asym{m}, mdlPrNon{m}, '-', 'Color', cmap_Non(m+1,:), 'LineWidth', 1); hold on 
end
% yline(PrFwd{3}, '--', 'Color', red, 'LineWidth', 1);
% yline(PrBck{3}, '--', 'Color', blue, 'LineWidth', 1);
% yline(PrNon{3}, '--', 'Color', grey, 'LineWidth', 1);
legend('B p(CS) = 0.15%, tauP = 1s', 'F p(CS) = 0.15% tauP = 1s', 'N p(CS) = 0.15% tauP = 1s','B p(CS) = 0.15%, tauP = 1.31s', 'F p(CS) = 0.15% tauP = 1.31s', 'N p(CS) = 0.15% tauP = 1.31s', 'B p(CS) = 0.2%', 'F p(CS) = 0.2%', 'N p(CS) = 0.2%', 'Location', 'BestOutside')
ylim([0 1]); 
xlabel('Asymmetry')
ylabel('fraction of PFs')
title('speed: 15cm/s, Pamp: 20 pA')
box off
axis square

mc = 1;
figure % Proportions as a function of Tau Asymetry for p(CS) = 0.15%
plot(Asym{mc}, mdlPrBck{mc}, 'o-', 'Color', blue, 'LineWidth', 1.5); hold on
plot(Asym{mc}, mdlPrFwd{mc}, 'o-', 'Color', red, 'LineWidth', 1.5); hold on 
plot(Asym{mc}, mdlPrNon{mc}, 'o-', 'Color', grey, 'LineWidth', 1.5); hold on 
% yline(PrFwd{3}, '--', 'Color', red, 'LineWidth', 1);
% yline(PrBck{3}, '--', 'Color', blue, 'LineWidth', 1);
% yline(PrNon{3}, '--', 'Color', grey, 'LineWidth', 1);
% legend('B', 'F', 'N', 'CA1F', 'CA1F', 'CA1F', 'Location', 'BestOutside')
ylim([0 1]);
xlabel('BTSP Asymmetry')
ylabel('proportion of PFs')
title('speed: 15cm/s, Pamp: 20pA')
box off
axis square

figure % rule   
Dmap = brewermap(length(tauD{1})+4,'Oranges'); %hot(length(tauD{1})+2);
    delay1 = -5000:1:0; % ms
    delay2 = 0:1:5000; %ms
    Wp = 20*exp(delay1/(tauP{1}(1)*1000));
    plot(delay1./1000, Wp, 'k', 'LineWidth', 1.5); hold on
    for n = 1:length(tauD{1})
    Wd(n,:) = 20*exp(-delay2/(tauD{1}(n)*1000));
    plot(delay2./1000, Wd(n,:), 'Color', Dmap(n+4,:), 'LineWidth', 1.5); hold on
    end
    legend('20 ms', '0.5 s', '1 s', '1.3 s')
    xline(0, 'k--'); hold on
    yline(0, 'k--'); 
%     xlim([-4.2 4.2])
%     ylim([0 22])
    xlabel('delay (s)'); ylabel('Potentiation (pA)');
    title(' tau_P = 1 s')
    box off; axis square;

%% plots

gp = cat(1, group{:});
COMslopeAll = cat(1, slopes{:});
PropBck = cat(1, PrBck{:});
PropFwd = cat(1, PrFwd{:});
PropNon = cat(1, PrNon{:});
Props = [PropBck PropFwd PropNon];
PropsAll = Props;
% PFs pooled from all animals
% PropsAll(1,:) = [CA1N.Bck CA1N.Fwd CA1N.ns];
% PropsAll(3,:) = [CA1F.Bck CA1F.Fwd CA1F.ns];

figure 
subplot(2,2,1) % violin plots. 
    v = violinplot(COMslopeAll, gp, 'ViolinColor', [0 0 0], 'ShowData', false, 'ShowMean', true);
    v(1).ViolinColor = mapCA1(2,:);
    v(3).ViolinColor = mapCA1(1,:);
    yline(0, 'k--')
    ylabel('Shifting slopes (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA3_N'; '0.005'; 'CA3_F'; '0.001']);
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
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA3_N'; '0.005'; 'CA3_F'; '0.001']);
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
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA3_N'; '0.005'; 'CA3_F'; '0.001']);
    box off;
    axis square
subplot(2,2,4) % props 
    b2 = bar(PropsAll, 'stacked', 'FaceColor','flat');
    b2(1).CData = blue; b2(2).CData = red; b2(3).CData = grey;
    ylabel('proportion of PFs'); xlabel('models')
    xlim([0.5 length(group)+0.5]);
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA3_N'; '0.005'; 'CA3_F'; '0.001']);
    title('all PFs pooled')
    box off;
    axis square

t = tiledlayout(2,4)
nexttile % rule   
    Dmap = brewermap(length(tauD{1})+4,'RdPu'); %hot(length(tauD{1})+2);
        delay1 = -5000:1:0; % ms
        delay2 = 0:1:5000; %ms
        Wp = 20*exp(delay1/(tauP{1}(1)*1000));
        plot(delay1./1000, Wp, 'k', 'LineWidth', 1.5); hold on
        for n = 1:length(tauD{1})
        Wd(n,:) = 20*exp(-delay2/(tauD{1}(n)*1000));
        plot(delay2./1000, Wd(n,:), 'Color', Dmap(n+4,:), 'LineWidth', 1.5); hold on
        end
%         legend('20 ms', '0.5 s', '1 s', '1.3 s', 'location', 'best outside')
        xline(0, 'k--'); hold on
        yline(0, 'k--'); 
    %     xlim([-4.2 4.2])
    %     ylim([0 22])
        xlabel('delay (s)'); ylabel('Potentiation (pA)');
        title(' tau_P = 1 s')
        box off; axis square;
nexttile(5) % Proportions as a function of Tau Asymetry for p(CS) = 0.15%
    plot(Asym{mc}, mdlPrBck{mc}, 'o-', 'Color', blue, 'LineWidth', 1.5); hold on
    plot(Asym{mc}, mdlPrFwd{mc}, 'o-', 'Color', red, 'LineWidth', 1.5); hold on 
    plot(Asym{mc}, mdlPrNon{mc}, 'o-', 'Color', grey, 'LineWidth', 1.5); hold on 
    ylim([0 1]);
    xlim([-1 0.35])
    xlabel('BTSP Asymmetry')
    ylabel('fraction of PFs')
    title('p(CS) = 0.15%')
    box off
    axis square
nexttile(2, [2,1]) % shift distribs violins
    v = violinplot(COMslopeAll, gp, 'ViolinColor', [0 0 0], 'ShowData', false, 'ShowMean', true);
    v(1).ViolinColor = mapCA3(2,:);
%     v(2).ViolinColor = Dmap(5,:);
    v(3).ViolinColor = mapCA3(1,:);
%     v(4).ViolinColor = Dmap(8,:);
    yline(0, 'k--')
    ylabel('Shifting slopes (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA3_N'; num2str(Asym{1}(1)); 'CA3_F'; '  ' num2str(Asym{1}(end))]);
    box off;
%     axis square
 nexttile(3, [2,1])   
    y = meanSlope; yneg = y-SlopeCIlo; ypos = SlopeCIup-y;
    errorbar(1:length(group), y,yneg,ypos,'o', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3); hold on
    yline(0, 'k-');
    yline(meanSlope(1), '--', 'Color', mapCA3(2,:), 'LineWidth', 1.5);
    yline(meanSlope(3), '--', 'Color', mapCA3(1,:), 'LineWidth', 1.5);
    xlim([0.5 length(group)+0.5]);
    ylim([-0.3 0.15]); 
%     xlabel('models'); 
    ylabel('mean slope (cm/lap)');
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA3_N'; num2str(Asym{1}(1)); 'CA3_F'; '  ' num2str(Asym{1}(end))]);
    box off;
%     axis square
 nexttile(4, [2,1]) 
    b2 = bar(PropsAll, 'stacked', 'FaceColor','flat');
    b2(1).CData = blue; b2(2).CData = red; b2(3).CData = grey;
    ylabel('proportion of PFs'); 
%     xlabel('models')
    xlim([0.5 length(group)+0.5]);
    set(gca, 'XTick', 1:length(group), 'XTickLabel', ['CA3_N'; num2str(Asym{1}(1)); 'CA3_F'; '  ' num2str(Asym{1}(end))]);
%     title('all PFs pooled')
    box off;
%     axis square
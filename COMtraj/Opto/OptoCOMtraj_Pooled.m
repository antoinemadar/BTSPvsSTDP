clear
close all

Minlaps = 15;
Mfirst = load("COMtrajMat10laps_optofirst.mat");
Mlater = load("COMtrajMat10laps_optolater.mat");
COMtrajMat2 = [Mfirst.COMtrajMat2; Mlater.COMtrajMat2];
Tok_first = readtable('Toptofirst.csv');
Tok_later = readtable('Toptolater.csv');
Tok = [Tok_first; Tok_later];

% filename = 'AnqiSourceDataFile.xlsx';
% writetable(Tok,filename);
%% color codes

%colorcode CA1 vs CA3, N (dark) vs F (light); e.g. CA1F -> mapCA1(1,:)
map = brewermap(10,'Reds'); %cool(4); %brewermap(4,'PiYG')
colors = [map(6,:); map(7,:); map(8,:); map(9,:)];
figure; scatter(1:4,1:4,200,colors, 'filled')

% Opto (blue), Ctrl black
mapOpto = map(8,:); %[1 0.2 0.3];
mapCtrl = [0 0 0];
colors = [mapOpto; mapCtrl];

% fwd PFs should be red, Bck PFs blue, NonSig grey or black
blue = [0 0.45 0.74]; % backward
red = [0.85 0.32 0.01]; % forward
grey = [0.5 0.5 0.5]; %nonsig
%% L and R combined, compare Opto_first_day1 with Control_first_day1

OptoIdx = find(Tok.Opto_ON==1);
CtrlIdx = find(Tok.Opto_ON==0);

% Diffusion analysis

TrajsOpto = COMtrajMat2(OptoIdx,:); 
TrajsCtrl = COMtrajMat2(CtrlIdx,:);

DispSqOpto = ( TrajsOpto-TrajsOpto(:,1) ).^2;
DispSqCtrl = ( TrajsCtrl-TrajsCtrl(:,1) ).^2;

msdCtrl = mean( DispSqCtrl , 1 );
msdOpto = mean( DispSqOpto, 1 );

DispSqCAall = {DispSqOpto; DispSqCtrl};
MSDc = [msdOpto', msdCtrl'];

figure(1)
    plot(0:Minlaps-1, msdOpto,'-', 'LineWidth',1,'Color', mapOpto); hold on 
    plot(0:Minlaps-1, msdCtrl,'-', 'LineWidth',1,'Color', mapCtrl); hold on 
    for n = 1:2 %size(MSDc, 2)
        for lap = 1:Minlaps
            msdCA_CI{n}(:,lap) = bootci(1000, @mean, DispSqCAall{n}(:,lap));
        end
    plot_ci(0:Minlaps-1, [MSDc(:,n) msdCA_CI{n}(1,:)' msdCA_CI{n}(2,:)'], 'MainLineColor', colors(n,:), 'PatchColor', colors(n,:), 'LineColor', colors(n,:), 'PatchAlpha', 0.5); hold on
    end
    xlabel('Laps after onset')
    ylabel('Mean Squared Displacement (cm^2)')
    box off
    axis square
    legend('Opto', 'Ctrl', 'Location', 'SouthEast')

% linear regression

figure(2) % Slope distributions
histogram(Tok.COM_slope(OptoIdx), 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'pdf');hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(Tok.COM_slope(CtrlIdx),'EdgeColor', mapCtrl, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'pdf'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
xline(mean(Tok.COM_slope(OptoIdx)),'Color', mapOpto, 'LineWidth', 1.5 );hold on
xline(mean(Tok.COM_slope(CtrlIdx)),'Color', mapCtrl, 'LineWidth', 1.5 );hold on
% xlim([-4, 3])
box off; axis square;
xlabel('Slope (cm/lap)'); ylabel('pdf');

figure % Slope vs R2
subplot(1,2,1)
gscatter(Tok.COM_slope(CtrlIdx),Tok.COM_R2(CtrlIdx), Tok.shiftDir(CtrlIdx), [blue; grey; red], 'o', 5, 'off')
xlabel('lin reg slope (cm/lap)'); ylabel('lin reg R2');
xlim([-5, 5])
% xlim([-4, 3])
title('Ctrl')
box off; axis square;
subplot(1,2,2)
gscatter(Tok.COM_slope(OptoIdx),Tok.COM_R2(OptoIdx), Tok.shiftDir(OptoIdx), [blue; grey; red], 'o', 5, 'off')
xlim([-5, 5])
% xlim([-4, 3])
xlabel('lin reg slope (cm/lap)'); ylabel('lin reg R2');
title('Opto')
box off; axis square;

% proportions of F, B and N
Opto_Fwd_idx = find(Tok.COM_slope(OptoIdx) > 0 & Tok.COM_pval(OptoIdx)<=0.05);
Opto_Bck_idx = find(Tok.COM_slope(OptoIdx) < 0 & Tok.COM_pval(OptoIdx)<=0.05);
Opto_nonsig_idx = find(Tok.COM_pval(OptoIdx)>0.05);

Prop_Bck_OptoAll = length(Opto_Bck_idx)./length(Tok.COM_slope(OptoIdx))
Prop_Fwd_OptoAll = length(Opto_Fwd_idx)./length(Tok.COM_slope(OptoIdx))
Prop_NonSig_OptoAll = length(Opto_nonsig_idx)./length(Tok.COM_slope(OptoIdx))

Ctrl_Fwd_idx = find(Tok.COM_slope(CtrlIdx) > 0 & Tok.COM_pval(CtrlIdx)<=0.05);
Ctrl_Bck_idx = find(Tok.COM_slope(CtrlIdx) < 0 & Tok.COM_pval(CtrlIdx)<=0.05);
Ctrl_nonsig_idx = find(Tok.COM_pval(CtrlIdx)>0.05);

Prop_Bck_CTRall = length(Ctrl_Bck_idx)./length(Tok.COM_slope(CtrlIdx))
Prop_Fwd_CTRall = length(Ctrl_Fwd_idx)./length(Tok.COM_slope(CtrlIdx))
Prop_NonSig_CTRall = length(Ctrl_nonsig_idx)./length(Tok.COM_slope(CtrlIdx))

% bootstrap proportion of shifting PFs
propshiftOpto = propshift(Tok.COM_pval(OptoIdx));
propshiftCtrl = propshift(Tok.COM_pval(CtrlIdx));
[propshiftCI_opto, propshiftboot_opto] = bootci(10000,@propshift, Tok.COM_pval(OptoIdx));
[propshiftCI_ctrl, propshiftboot_ctrl] = bootci(10000,@propshift, Tok.COM_pval(CtrlIdx));
meanPropOpto = mean(propshiftboot_opto);
meanPropCtrl = mean(propshiftboot_ctrl);

DeltaProp = propshiftboot_opto - propshiftboot_ctrl; 
DeltaPropCI = prctile(DeltaProp,[2.5, 97.5]);
pvalProp = length(find(DeltaProp>=0))/length(DeltaProp);

% figure(3)
% subplot(1,3,1)% proportion of forward, backward and stable/non-signif PFs (mean across animals)
%     b = bar([Prop_Bck_CTRall Prop_Fwd_CTRall Prop_NonSig_CTRall; Prop_Bck_OptoAll Prop_Fwd_OptoAll Prop_NonSig_OptoAll], 'stacked', 'FaceColor','flat');
%     b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
%     ylabel('ratio of PFs');
%     set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
%     box off; axis square;
% subplot(1,3,2) 
%     errorbar(1 , meanPropCtrl, meanPropCtrl-propshiftCI_ctrl(1), meanPropCtrl-propshiftCI_ctrl(2), 'o', 'Color', mapCtrl, 'MarkerFaceColor', mapCtrl); hold on
%     errorbar(2 , meanPropOpto, meanPropOpto-propshiftCI_opto(1), meanPropOpto-propshiftCI_opto(2), 'o', 'Color', mapOpto, 'MarkerFaceColor', mapOpto);
%     ylabel('ratio of shifting PFs')
%     xlim([0 3]); ylim([0 0.5])
%     set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
%     box off; axis square;
%     title('bootstrapped 95% CIs')
% subplot(1,3,3)
%     histogram(DeltaProp, 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
%     xline(mean(DeltaProp),'Color', mapOpto, 'LineWidth', 1.5 );hold on
%     xline(0,'k', 'LineWidth', 1.5 );hold on
%     xline(DeltaPropCI, '--', 'Color', mapOpto, 'LineWidth', 1.5 );hold on
%     box off; axis square;
%     xlabel('\Delta ratio (Opto - Ctrl)'); ylabel('probability');
%     title('bootstrapped distribution of Opto - Ctrl')
% 
% figure(4)
%     h = histogram(DeltaProp, 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'cdf');hold on
%     xline(0,'k', 'LineWidth', 1 );hold on
%     box off; axis square;
%     xlabel('\Delta ratio (Opto - Ctrl)'); ylabel('cdf');
%     title('bootstrapped distribution of Opto - Ctrl')
%     BinCenters = h.BinEdges(1:end-1)+h.BinEdges(2:end)./2;
%     PropDecreased = interp1(BinCenters,h.Values,0).*100; % in percent
% 
%     Pvalbootdiff = 1-PropDecreased./100;

% bootstrap proportion of shifting PFs
figure
subplot(1,2,1)% proportion of forward, backward and stable/non-signif PFs (mean across animals)
    b = bar([Prop_Bck_CTRall Prop_Fwd_CTRall Prop_NonSig_CTRall; Prop_Bck_OptoAll Prop_Fwd_OptoAll Prop_NonSig_OptoAll], 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylabel('ratio of PFs');
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
    box off; axis square;
% subplot(1,3,2) 
%     errorbar(1 , meanPropCtrl, meanPropCtrl-propshiftCI_ctrl(1), meanPropCtrl-propshiftCI_ctrl(2), 'o', 'Color', mapCtrl, 'MarkerFaceColor', mapCtrl); hold on
%     errorbar(2 , meanPropOpto, meanPropOpto-propshiftCI_opto(1), meanPropOpto-propshiftCI_opto(2), 'o', 'Color', mapOpto, 'MarkerFaceColor', mapOpto);
%     ylabel('ratio of shifting PFs')
%     xlim([0 3]); ylim([0 0.4])
%     set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
%     box off; axis square;
%     title('bootstrapped 95% CIs')
subplot(1,2,2)
    histogram(DeltaProp, 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaProp), 0.01, mean(DeltaProp)-DeltaPropCI(1),mean(DeltaProp)-DeltaPropCI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-meanPropCtrl 1-meanPropCtrl]);
    xlabel('\Delta fraction of shifting PFs (Opto - Ctrl)'); ylabel('proba');
    title(['bootstrapped exact test: p = ' num2str(pvalProp)])
    set(gca,'view',[90 -90])
    box off; axis square;

%% Lin reg: delta-delta bootstrapped stats for median abs(slope)
[slopeBackCI_CA1N, slopeBackboot_CA1N] = bootci(10000,@median, abs(Tok.COM_slope(CtrlIdx)));
[slopeBackCI_CA1F, slopeBackboot_CA1F] = bootci(10000,@median, abs(Tok.COM_slope(OptoIdx)));
meanSlopeBCA1N = mean(slopeBackboot_CA1N);
meanSlopeBCA1F = mean(slopeBackboot_CA1F);

    % pairwise comparisons N vs F
    DeltaSlopeBCA1_NvsF = slopeBackboot_CA1F - slopeBackboot_CA1N; 
    DeltaSlopeBCA1_NvsF_CI = prctile(DeltaSlopeBCA1_NvsF,[2.5, 97.5]);
    pvalSlopeCA1 = length(find(DeltaSlopeBCA1_NvsF>0))/length(DeltaSlopeBCA1_NvsF);

figure % all PFs pooled. Bar plots of props for all Ca and conditions. Delta F vs N for CA1 and CA3, oriented like in Dabest. 
subplot(1,2,1) % abs(slopes) (all PFs pooled)
    v3 = violinplot(abs(Tok.COM_slope), Tok.Opto_ON);
    v3(1).ViolinColor = mapCtrl;
    v3(2).ViolinColor = mapOpto;
    ylabel('|slope|');
    set(gca, 'XTick', [1,2], 'XTickLabel', ['Ctrl';'Opto']);
    box off; %axis square;
subplot(1,2,2) % CA1 N vs F bootstrapped delta |slope|
    histogram(DeltaSlopeBCA1_NvsF, 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaSlopeBCA1_NvsF), 0.01, mean(DeltaSlopeBCA1_NvsF)-DeltaSlopeBCA1_NvsF_CI(1),mean(DeltaSlopeBCA1_NvsF)-DeltaSlopeBCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    % xlim([-meanSlopeBCA1F 5-meanSlopeBCA1F]);
    xlim([-0.5 0.5]);
    xlabel('\Delta median |slope| (Opto - Ctrl)'); ylabel('proba');
    title('bootstrapped distrib')
    set(gca,'view',[90 -90])
    box off; %axis square;

figure
subplot(1,4,1) % abs(slopes) (all PFs pooled)
    v3 = violinplot(abs(Tok.COM_slope), Tok.Opto_ON);
    v3(1).ViolinColor = mapCtrl;
    v3(2).ViolinColor = mapOpto;
    ylabel('|slope|');
    set(gca, 'XTick', [1,2], 'XTickLabel', ['Ctrl';'Opto']);
    box off; %axis square;
subplot(1,4,2) % CA1 N vs F bootstrapped delta |slope|
    histogram(DeltaSlopeBCA1_NvsF, 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaSlopeBCA1_NvsF), 0.01, mean(DeltaSlopeBCA1_NvsF)-DeltaSlopeBCA1_NvsF_CI(1),mean(DeltaSlopeBCA1_NvsF)-DeltaSlopeBCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 5); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-meanSlopeBCA1F 5-meanSlopeBCA1F]);
    % xlim([-0.3 0.3]);
    xlabel('\Delta median |slope| (Opto - Ctrl)'); ylabel('proba');
    title('bootstrapped distrib')
    set(gca,'view',[90 -90])
    box off;
subplot(1,4,3)% proportion of forward, backward and stable/non-signif PFs (mean across animals)
    b = bar([Prop_Bck_CTRall Prop_Fwd_CTRall Prop_NonSig_CTRall; Prop_Bck_OptoAll Prop_Fwd_OptoAll Prop_NonSig_OptoAll], 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylabel('ratio of PFs');
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
    box off; %axis square;
subplot(1,4,4)
    histogram(DeltaProp, 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaProp), 0.01, mean(DeltaProp)-DeltaPropCI(1),mean(DeltaProp)-DeltaPropCI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 5); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-meanPropCtrl 1-meanPropCtrl]);
    xlabel('\Delta fraction of shifting PFs (Opto - Ctrl)'); ylabel('proba');
    title(['bootstrapped exact test: p = ' num2str(pvalProp)])
    set(gca,'view',[90 -90])
    box off; %axis square;
%% NL reg

[propshiftCI_CA1N, propshiftboot_CA1N] = bootci(10000,@propshiftNL, Tok.RCreg_R2(CtrlIdx));
[propshiftCI_CA1F, propshiftboot_CA1F] = bootci(10000,@propshiftNL, Tok.RCreg_R2(OptoIdx));
meanPropCA1N = mean(propshiftboot_CA1N);
meanPropCA1F = mean(propshiftboot_CA1F);

    % pairwise comparisons Ctrl vs Opto
    DeltaPropCA1_NvsF = propshiftboot_CA1F - propshiftboot_CA1N; 
    DeltaPropCA1_NvsF_CI = prctile(DeltaPropCA1_NvsF,[2.5, 97.5]);
    pvalCA1_NvsF = length(find(DeltaPropCA1_NvsF>0))/length(DeltaPropCA1_NvsF);

% delta bootstrapped stats for Backward shifting PFs
[propBackCI_CA1N, propBackboot_CA1N] = bootci(10000,@propBackNL, Tok.RCreg_R2(CtrlIdx),Tok.RCreg_p1(CtrlIdx));
[propBackCI_CA1F, propBackboot_CA1F] = bootci(10000,@propBackNL, Tok.RCreg_R2(OptoIdx),Tok.RCreg_p1(OptoIdx));
meanPropBCA1N = mean(propBackboot_CA1N);
meanPropBCA1F = mean(propBackboot_CA1F);

    % pairwise comparisons Ctrl vs Opto
    DeltaPropBCA1_NvsF = propBackboot_CA1F - propBackboot_CA1N; 
    DeltaPropBCA1_NvsF_CI = prctile(DeltaPropBCA1_NvsF,[2.5, 97.5]);
    pvalBnl = length(find(DeltaPropBCA1_NvsF>0))/length(DeltaPropBCA1_NvsF);
    
% delta bootstrapped stats for Forward shifting PFs
[propFwdCI_CA1N, propFwdboot_CA1N] = bootci(10000,@propFwdNL, Tok.RCreg_R2(CtrlIdx),Tok.RCreg_p1(CtrlIdx));
[propFwdCI_CA1F, propFwdboot_CA1F] = bootci(10000,@propFwdNL, Tok.RCreg_R2(OptoIdx),Tok.RCreg_p1(OptoIdx));
meanPropFCA1N = mean(propFwdboot_CA1N);
meanPropFCA1F = mean(propFwdboot_CA1F);

    % pairwise comparisons Ctrl vs Opto
    DeltaPropFCA1_NvsF = propFwdboot_CA1F - propFwdboot_CA1N; 
    DeltaPropFCA1_NvsF_CI = prctile(DeltaPropFCA1_NvsF,[2.5, 97.5]);
    pvalFnl = length(find(DeltaPropFCA1_NvsF>0))/length(DeltaPropFCA1_NvsF);

figure % all PFs pooled. 
subplot(1,3,1) % proportion of forward, backward and stable/non-signif PFs (all PFs pooled)
    b = bar([meanPropBCA1N meanPropFCA1N 1-meanPropCA1N; meanPropBCA1F meanPropFCA1F 1-meanPropCA1F], 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylim([0 1]);
    ylabel('proportion of NL shifting PFs');
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
    box off; %axis square;
subplot(1,3,2) % CA1 N vs F bootstrapped delta for Backward
    histogram(DeltaPropBCA1_NvsF, 'EdgeColor', blue, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaPropBCA1_NvsF), 0.01, mean(DeltaPropBCA1_NvsF)-DeltaPropBCA1_NvsF_CI(1),mean(DeltaPropBCA1_NvsF)-DeltaPropBCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-0.2 0.2]);
    xlabel('\Delta ratio (Opto - Ctrl)'); ylabel('proba');
    title(['bootstrapped exact test: p = ' num2str(pvalBnl)])
    set(gca,'view',[90 -90])
    box off; %axis square;
subplot(1,3,3) % forward
    histogram(DeltaPropFCA1_NvsF, 'EdgeColor', red, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaPropFCA1_NvsF), 0.01, mean(DeltaPropFCA1_NvsF)-DeltaPropFCA1_NvsF_CI(1),mean(DeltaPropFCA1_NvsF)-DeltaPropFCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', red, 'MarkerFaceColor', red, 'LineWidth', 1, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-0.2 0.2]);
    xlabel('\Delta ratio (Opto - Ctrl)'); ylabel('proba');
    title(['bootstrapped exact test: p = ' num2str(pvalFnl)])
    set(gca,'view',[90 -90])
    box off;

figure % all PFs pooled. Shifting PFs
subplot(1,2,1) % proportion of forward, backward and stable/non-signif PFs (all PFs pooled)
    b = bar([meanPropBCA1N meanPropFCA1N 1-meanPropCA1N; meanPropBCA1F meanPropFCA1F 1-meanPropCA1F], 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylim([0 1]);
    ylabel('fraction of PFs');
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
    box off; %axis square;
subplot(1,2,2) % CA1 N vs F bootstrapped delta for Backward
    histogram(DeltaPropCA1_NvsF, 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaPropCA1_NvsF), 0.01, mean(DeltaPropCA1_NvsF)-DeltaPropCA1_NvsF_CI(1),mean(DeltaPropCA1_NvsF)-DeltaPropCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-meanPropCA1N 1-meanPropCA1N]);
    xlabel('\Delta fraction of NL shifting PFs (Opto - Ctrl)'); ylabel('proba');
    title(['bootstrapped exact test: p = ' num2str(pvalCA1_NvsF)])
    set(gca,'view',[90 -90])
    box off; %axis square;


%% Difference in NL Amp distrib?
% delta-delta bootstrapped stats for median abs(AMp)
[AmpCI_CA1N, ampBackboot_CA1N] = bootci(10000,@median, abs(Tok.RCreg_p1(CtrlIdx)));
[AmpCI_CA1F, ampBackboot_CA1F] = bootci(10000,@median, abs(Tok.RCreg_p1(OptoIdx)));
meanAmpCA1N = mean(ampBackboot_CA1N);
meanAmpCA1F = mean(ampBackboot_CA1F);

    % pairwise comparisons N vs F
    DeltaAmpCA1_NvsF = ampBackboot_CA1F - ampBackboot_CA1N; 
    DeltaAmpCA1_NvsF_CI = prctile(DeltaAmpCA1_NvsF,[2.5, 97.5]);
    pvalAmp = length(find(DeltaAmpCA1_NvsF>0))/length(DeltaAmpCA1_NvsF);

figure % all PFs pooled. median |Amp|
subplot(1,3,1) % abs(Amp) (all PFs pooled)
    v3 = violinplot(abs(Tok.RCreg_p1), Tok.Opto_ON);
    v3(1).ViolinColor = mapCtrl;
    v3(2).ViolinColor = mapOpto;
    ylabel('|Amp|, cm');
    set(gca, 'XTick', [1,2], 'XTickLabel', ['Ctrl';'Opto']);
    box off; %axis square;
subplot(1,3,2)
    plot([1,2], [meanAmpCA1N meanAmpCA1F], 'k-'); hold on
    errorbar(2 , meanAmpCA1F, meanAmpCA1F-AmpCI_CA1F(1), meanAmpCA1F-AmpCI_CA1F(2), 'o', 'Color', mapOpto, 'MarkerFaceColor', mapOpto); hold on
    errorbar(1 , meanAmpCA1N, meanAmpCA1N-AmpCI_CA1N(1), meanAmpCA1N-AmpCI_CA1N(2), 'o', 'Color', mapCtrl, 'MarkerFaceColor', mapCtrl); hold on
    ylabel('median |Amp|, cm')
    xlim([0 3]); ylim([0 32])
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
    box off; %axis square;
    title('bootstrapped 95% CIs')
subplot(1,3,3) % CA1 N vs F bootstrapped delta |slope|
    histogram(DeltaAmpCA1_NvsF, 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaAmpCA1_NvsF), 0.01, mean(DeltaAmpCA1_NvsF)-DeltaAmpCA1_NvsF_CI(1),mean(DeltaAmpCA1_NvsF)-DeltaAmpCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-20 20]);
    xlabel('\Delta median |Amp| (F - N)'); ylabel('proba');
    title(['bootstrapped exact test: p = ' num2str(pvalAmp)])
    set(gca,'view',[90 -90])
    box off; %axis square;

%% per animals

mice = unique(Tok.mouse);

% Diffusion analysis

% resampling paired test of difference between Opto and Ctrl, taking in account repeated measures in animals: 
% for each animal, take a random sample (with replacement) of COM
% trajectories with the original sample size for that animal; then take the
% difference between MSD in ctrl vs opto. Repeat 1000x, for each animal.
% Combine all the the animal-wise bootstrapped distributions (1000 x number of animals) to obtain the full
% bootstrap distribution. 
% Check whether the mean of the bootstrap distribution is similar to the
% mean of the animal-wise average of MSD -> yes
% Compute 95% CI of the animal-wise mean by 1000x resampling 1 DiffMSDboot per animal and computing the mean across animals. 
%note that this does not take care of the multiple comparison issue... but not sure there is one. 

for m = 1:length(mice)
    
    % animal-wise MSD difference between Opto and Ctrl
    OptoIdx_mouse{m} = find(ismember(Tok.mouse, mice{m}) & Tok.Opto_ON==1);
    CtrlIdx_mouse{m} = find(ismember(Tok.mouse, mice{m}) & Tok.Opto_ON==0);
    TrajsOpto_mouse{m} = COMtrajMat2(OptoIdx_mouse{m},:); 
    TrajsCtrl_mouse{m} = COMtrajMat2(CtrlIdx_mouse{m},:);
    msdCtrl_mouse{m} = mean( TrajsCtrl_mouse{m}.^2 , 1 );
    msdOpto_mouse{m} = mean( TrajsOpto_mouse{m}.^2, 1 );
    DiffMSD(m,:) = msdCtrl_mouse{m} - msdOpto_mouse{m}; 

    if ~isempty(OptoIdx_mouse{m})
    % resampled distribution of MSD
%     OptoSampsize{m} = length(OptoIdx_mouse{m});
%     CtrlSampsize{m} = length(CtrlIdx_mouse{m});
    bootstat_Opto{m} = bootstrp(1000,@mean,TrajsOpto_mouse{m}.^2); % 1000 resamples with replacement (each of original sample size). Should give a 1000 rows x length(Trajectories) matrix, i.e. MSD for a given sample x lap)
    bootstat_Ctrl{m} = bootstrp(1000,@mean,TrajsCtrl_mouse{m}.^2); % 1000 resamples with replacement
    DiffMSDboot{m} = bootstat_Ctrl{m} - bootstat_Opto{m}; % 1000 DiffMSD x number of laps in trajectories 
        for l = 1:size(DiffMSDboot{m},2)
        DiffMSDboot2{l}(:,m) = randsample(DiffMSDboot{m}(:,l),10000,true); % random sampling with replacement
        end
    else
    DiffMSDboot{m} = ones(1000,Minlaps).*NaN;  
        for l = 1:size(DiffMSDboot{m},2)
        DiffMSDboot2{l}(:,m) = ones(10000,1).*NaN; 
        end
    end
end

meanDiffMSD = mean(DiffMSD, 1, 'omitnan'); % raw lapwise MSD averaged across animals 
DiffMSDpermouse = DiffMSD(:);
LapsPerMouse = repmat(0:Minlaps-1, size(DiffMSD,1),1);

DiffMSDboot_all = cat(1,DiffMSDboot{:}); % Full bootstrapped distribution per lap: combine animal-wise bootstrapped distributions of DiffMSD
meanDiffMSDboot = mean(DiffMSDboot_all,1, 'omitnan');
DiffMSDboot_pctl = prctile(DiffMSDboot_all,[2.5, 97.5],1); % 2.5pctl on top row, 97.5pctl on bottom row

% compute bootstrap of the mean and its 95% CI
for l = 1:Minlaps
meanDiffMSDboot2(:,l) = mean(DiffMSDboot2{l},2, 'omitnan');  % 10,000 column vector
end
meanDiffMSDboot2_estim = mean(meanDiffMSDboot2,1); % row vector of Minlaps values
meanDiffMSDboot2_CI = prctile(meanDiffMSDboot2,[2.5, 97.5],1);
meanDiffMSDboot2_onesided = prctile(meanDiffMSDboot2,5,1);

figure % double bootstrapped paired test that MSD is different between Opto and Ctrl
    scatter(LapsPerMouse(:), DiffMSDpermouse, 'go'); hold on
    plot(0:Minlaps-1, meanDiffMSD,'-', 'LineWidth',1,'Color', 'g'); hold on 
    plot_ci(0:Minlaps-1, [meanDiffMSDboot2_estim' meanDiffMSDboot2_CI(1,:)' meanDiffMSDboot2_CI(2,:)'], 'MainLineColor', 'k', 'PatchColor', 'k', 'LineColor', 'k', 'PatchAlpha', 0.5); hold on
    plot(0:Minlaps-1, meanDiffMSDboot2_onesided,'-', 'LineWidth',1,'Color', 'b'); hold on
    plot([0 length(meanDiffMSDboot')], [0 0], 'r--');
    xlabel('Laps after onset')
    ylabel('\DeltaMSD (cm^2)')
    box off
    axis square
%     legend('average across animals', 'animal-wise bootstrapped 95% CI', 'Location', 'BestOutside')
    legend('for each animal', 'average across animals', 'animal-wise paired bootstrapped 95% CI of the mean', 'one-sided a=0.05 test',  'Location', 'BestOutside')

figure % bootstrapped distribution of the DiffMSD
subplot(1,3,1)
    scatter(LapsPerMouse(:), DiffMSDpermouse, 'ko'); hold on
    plot_ci(0:Minlaps-1, [meanDiffMSDboot' DiffMSDboot_pctl(1,:)' DiffMSDboot_pctl(2,:)'], 'MainLineColor', 'g', 'PatchColor', 'k', 'LineColor', 'k', 'PatchAlpha', 0.5); hold on
%     plot(0:Minlaps-1, meanDiffMSD,'-', 'LineWidth',1,'Color', 'g'); hold on 
    plot([0 length(meanDiffMSDboot')], [0 0], 'r--');
    xlabel('Laps after onset')
    ylabel('\DeltaMSD (cm^2)')
    box off
    axis square
%     legend('average across animals', 'animal-wise bootstrapped 95% CI', 'Location', 'BestOutside')
    legend('for each animal', 'animal-wise paired bootstrapped 95% pctiles', 'average across animals', 'Location', 'SouthEast')
subplot(1,3,2)
    c{2} = 1:Minlaps-1; c{1}= -250:10:250;
    dummy = repmat(c{2},size(DiffMSDboot_all,1),1);
    dummy2 = DiffMSDboot_all(:,2:end);
    D(:,2) = dummy(:); D(:,1) = dummy2(:); 
    [N,cc] = hist3(D,'Ctrs',c);
    [X,Y] = meshgrid(c{2},c{1});
    contour(X,Y,N);hold on
    plot([0 length(meanDiffMSDboot')], [0 0], 'k--');
    xlabel('Laps after onset'); ylabel('\DeltaMSD (cm^2)');
    box off; axis square;
    cb = colorbar;
    cb.Label.String = 'PFs';
subplot(1,3,3)
    histogram2(D(:,2), D(:,1), 'DisplayStyle','tile','ShowEmptyBins','off')
    xlabel('Laps after onset'); ylabel('\DeltaMSD (cm^2)');
    box off; axis square;
    cb = colorbar;
    cb.Label.String = 'PFs'


% lin reg proportions per animal and conditions

Ctrl_Fwd_idx = find(Tok.COM_slope(CtrlIdx)>0 & Tok.COM_pval(CtrlIdx)<=0.05);

for m = 1:length(mice)
    IdxFwd_Opto{m} = find(Tok.COM_slope>0 & Tok.COM_pval<=0.05 & ismember(Tok.mouse, mice{m}) & Tok.Opto_ON==1); % & ismember(Tok.Region, 'R')
    IdxBck_Opto{m} = find(Tok.COM_slope<0 & Tok.COM_pval<=0.05 & ismember(Tok.mouse, mice{m}) & Tok.Opto_ON==1);
    IdxNonSig_Opto{m} = find(Tok.COM_pval>0.05 & ismember(Tok.mouse, mice{m}) & Tok.Opto_ON==1);
    IdxO{m} = [IdxFwd_Opto{m}; IdxBck_Opto{m}; IdxNonSig_Opto{m}];
    length_Opto(m) = length(IdxFwd_Opto{m}) + length(IdxBck_Opto{m}) + length(IdxNonSig_Opto{m});
    length_Opto2(m) = length(find(ismember(Tok.mouse, mice{m}) & Tok.Opto_ON==1));

    IdxFwd_Ctrl{m} = find(Tok.COM_slope>0 & Tok.COM_pval<=0.05 & ismember(Tok.mouse, mice{m}) & Tok.Opto_ON==0);
    IdxBck_Ctrl{m} = find(Tok.COM_slope<0 & Tok.COM_pval<=0.05 & ismember(Tok.mouse, mice{m}) & Tok.Opto_ON==0);
    IdxNonSig_Ctrl{m} = find(Tok.COM_pval>0.05 & ismember(Tok.mouse, mice{m}) & Tok.Opto_ON==0);
    IdxC{m} = [IdxFwd_Ctrl{m}; IdxBck_Ctrl{m}; IdxNonSig_Ctrl{m}];
    length_Ctrl(m) = length(IdxFwd_Ctrl{m}) + length(IdxBck_Ctrl{m}) + length(IdxNonSig_Ctrl{m});
    length_Ctrl2(m) = length(find(ismember(Tok.mouse, mice{m}) & Tok.Opto_ON==0));
    
    PropFwd_Opto(m) = length(IdxFwd_Opto{m})/length_Opto(m);
    PropBck_Opto(m) = length(IdxBck_Opto{m})/length_Opto(m);
    PropNonSig_Opto(m) = length(IdxNonSig_Opto{m})/length_Opto(m);
    PropShift_Opto(m) = PropFwd_Opto(m) + PropBck_Opto(m);
    PropTot_Opto(m) = PropFwd_Opto(m) + PropBck_Opto(m) + PropNonSig_Opto(m);

    PropFwd_Ctrl(m) = length(IdxFwd_Ctrl{m})/length_Ctrl(m);
    PropBck_Ctrl(m) = length(IdxBck_Ctrl{m})/length_Ctrl(m);
    PropNonSig_Ctrl(m) = length(IdxNonSig_Ctrl{m})/length_Ctrl(m);
    PropShift_Ctrl(m) = PropFwd_Ctrl(m) + PropBck_Ctrl(m);
    PropTot_Ctrl(m) = PropFwd_Ctrl(m) + PropBck_Ctrl(m) + PropNonSig_Ctrl(m);

    meanSlopeCA1F(m) = mean(Tok.COM_slope(IdxO{m}));
    meanSlopeCA1N(m) = mean(Tok.COM_slope(IdxC{m}));
    medianASlopeCA1F(m) = median(abs(Tok.COM_slope(IdxO{m})));
    medianASlopeCA1N(m) = median(abs(Tok.COM_slope(IdxC{m})));

    region(m) = unique(Tok.Region(ismember(Tok.mouse, mice{m})));
end
PropFwd_CA1 = [PropFwd_Ctrl;PropFwd_Opto];
PropBck_CA1 = [PropBck_Ctrl;PropBck_Opto];
PropNonSig_CA1 = [PropNonSig_Ctrl;PropNonSig_Opto];

% dumF = cat(1,IdxFwd_Opto{:});
% dumB = cat(1,IdxBck_Opto{:});
% dumN = cat(1,IdxNonSig_Opto{:});
% propB = length(dumB)/sum(length_Opto)
% propF = length(dumF)/sum(length_Opto)
% propN = length(dumN)/sum(length_Opto)
% propB+propF+propN
% 
% Prop_Bck_OptoAll = length(Opto_Bck_idx)./length(Tok.COM_slope(OptoIdx))
% Prop_Fwd_OptoAll = length(Opto_Fwd_idx)./length(Tok.COM_slope(OptoIdx))
% Prop_NonSig_OptoAll = length(Opto_nonsig_idx)./length(Tok.COM_slope(OptoIdx))
% Prop_Bck_OptoAll+Prop_Fwd_OptoAll+Prop_NonSig_OptoAll
% 
% dumFc = cat(1,IdxFwd_Ctrl{:});
% dumBc = cat(1,IdxBck_Ctrl{:});
% dumNc = cat(1,IdxNonSig_Ctrl{:});
% propBc = length(dumBc)/sum(length_Ctrl)
% propFc = length(dumFc)/sum(length_Ctrl)
% propNc = length(dumNc)/sum(length_Ctrl)
% propBc+propFc+propNc
% 
% Prop_Bck_CTRall = length(Ctrl_Bck_idx)./length(Tok.COM_slope(CtrlIdx))
% Prop_Fwd_CTRall = length(Ctrl_Fwd_idx)./length(Tok.COM_slope(CtrlIdx))
% Prop_NonSig_CTRall = length(Ctrl_nonsig_idx)./length(Tok.COM_slope(CtrlIdx))
% Prop_Bck_CTRall+Prop_Fwd_CTRall+Prop_NonSig_CTRall

[H_fwdCA1,P_fwdCA1,CI_fwdCA1,STATS_fwdCA1] = ttest(diff(PropFwd_CA1)); % p = 
[H_bckCA1,P_bckCA1,CI_bckCA1,STATS_bckCA1] = ttest(diff(PropBck_CA1)); % p = 
[H_nsCA1,P_nsCA1,CI_nsCA1,STATS_nsCA1] = ttest(diff(PropNonSig_CA1)); % p = 
[p,h,stats] = signrank(diff(PropNonSig_CA1));

figure
subplot(1,2,1)
    plot(PropFwd_CA1,'-o', 'Color', red); hold on
    plot(PropBck_CA1,'-o', 'Color', blue); hold on
    plot(PropNonSig_CA1,'-o', 'Color', grey); hold on
    plot(mean(PropFwd_CA1,2, 'omitnan'),'-', 'Color', red, 'LineWidth', 2); hold on
    plot(mean(PropBck_CA1,2, 'omitnan'),'-', 'Color', blue, 'LineWidth', 2); hold on
    plot(mean(PropNonSig_CA1,2, 'omitnan'),'-', 'Color', grey, 'LineWidth', 2); hold on
    ylim([0 1]); xlim([0 3]);
    ylabel('Prop of PFs'); 
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
    title('Shifting proportions per animal');
    box off; axis square
subplot(1,2,2) % proportion of forward, backward and stable/non-signif PFs (mean across animals)
    b = bar([mean(PropBck_Ctrl, 'omitnan') mean(PropFwd_Ctrl, 'omitnan') mean(PropNonSig_Ctrl, 'omitnan'); mean(PropBck_Opto, 'omitnan') mean(PropFwd_Opto, 'omitnan') mean(PropNonSig_Opto, 'omitnan')], 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylabel('proportion of PFs');
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
    box off; axis square;

%% local functions

function ratio = propshift(pvalArray)
% computes the proportion of significantly shifting PFs
ratio = length(find(pvalArray<=0.05))./length(pvalArray);
end

function ratio = propshiftNL(R2Array)
% computes the proportion of significantly shifting PFs
ratio = length(find(R2Array>=0.4))./length(R2Array);
end

function ratio = propBackNL(R2Array, Amp)
% computes the proportion of significantly backward shifting PFs
%R2array and Amp vectors must be the same length
ratio = length(R2Array(R2Array>=0.4 & Amp<0))./length(R2Array);
end

function ratio = propFwdNL(R2Array, Amp)
% computes the proportion of significantly backward shifting PFs
%pval and slope must be the same length
ratio = length(R2Array(R2Array>=0.4 & Amp>0))./length(R2Array);
end
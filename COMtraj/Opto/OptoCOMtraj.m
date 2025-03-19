clear 
close all

Minlaps = 15; % min number of laps where PF is defined (after interpolation, if interp is used)
restrict = 1; % 0 if analysis on full trajectories for PFs that appeared during Opto ON laps. Full trajectory, not restricted to laps during Opto ON. To compare to PFs appeared in Control during same laps.
% 1 if COMtrajs restricted to laps during Opto ON (which may make too few laps, and too few PFs _ though 8 animals were at 15+ Opto laps)  
binsize = 5; % 5cm. The track was 2m divided in 40 bins and Anqi's table expresses COM position in bins

T_opto = readtable('opto_length_animals.csv', ReadVariableNames=true);

% T_L = readtable('left_df.csv', ReadVariableNames=true);
% T_R = readtable('right_df.csv', ReadVariableNames=true);
% T_L.Region = repmat('L',height(T_L),1);
% T_R.Region = repmat('R',height(T_R),1);
% 
% T_all = [T_L; T_R]; 

T_all = readtable('all_PFs_lap45_NewCriteriaNominDF.csv', ReadVariableNames=true);

%% process data table

COMtrajMat = [];
COMtrajMatNaN = [];
SmallTraj = 0;
Opto = [];
% T.COM_slope = [];                
% T.COM_intcp = [];
% T.COM_R2 = [];
% T.COM_pval = [];  
% T_all.shiftDir = {};
% T.COMinterp_slope = [];
% T.COMinterp_intcp = [];                
% T.COMinterp_R2 = [];                
% T.COMinterp_pval = [];    
% T.RCreg_p1 = [];
% T.RCreg_p2 = [];
% T.RCreg_p3 = [];
% T.RCreg_R2 = [];
% T.RCreg_R2adj = [];

for i = 1:height(T_all)
    midx = find(strcmp(T_all.mouse(i),T_opto.mouse));
    if strcmp(T_all.env(i),'opto_first_day1')||strcmp(T_all.env(i),'control_day1')
        if T_all.emergeLap(i) <= T_opto.opto_first_day1(midx) % If emergeLap <= opto_first_day1 lap (which is the last lap when opto was on), add it to matrix of COMtrajs and compute things.
                
            T_all.PFopto(i) = 1; % add column PF_opto: 1 for PFs that appeared during laps when opto was on. (or equivalent laps for Control) and 0 if not.)
            
            % extract COM trajectory, starting from onset lap, and centered to initial COM location
            OnsetLap = T_all.emergeLap(i);
            OnsetLoc = (T_all{i, 21+OnsetLap}+1).*binsize; % +1 because COM bin is python indexed, with first bin = 0

            if restrict == 0 % analysis on all laps after emergence
            COMtraj = (T_all{i,21+OnsetLap:21+44}+1).*binsize - OnsetLoc; 
            elseif restrict == 1 % analysis on trajectories restricted to laps under opto
            COMtraj = (T_all{i,21+OnsetLap:21+T_opto.opto_first_day1(midx)}+1).*binsize - OnsetLoc; 
            end
            TrajActiveLaps = find(COMtraj);

            % Interpolate NaNs
            x = COMtraj; nanx = isnan(x); laps = 1:numel(x);
            if length(x(~nanx))>1 %exclude PFs with only 1 point
            x(nanx) = interp1(laps(~nanx), x(~nanx), laps(nanx)); % Linear extrapolation. No extrapolation for NaNs at the end
%                 x(nanx) = interp1(laps(~nanx), x(~nanx), laps(nanx), 'previous'); 
            COMtraj2 = x(~isnan(x)); % get rid of end NaNs

            TrajLength = length(COMtraj2);
                if TrajLength >= Minlaps && length(TrajActiveLaps)>=0.3*TrajLength
                    %regression to get COM shifting slope on raw and interp COMtraj
                    [b,~,~,~,stats] = regress(COMtraj2', [ones(TrajLength,1), [1:TrajLength]']); % on interp
        %                 COMtraj_intLR{k} = b(1) + b(2).*[1:TrajLength];
        
                    [b2,~,~,~,stats2] = regress(COMtraj', [ones(numel(x),1), [1:numel(x)]']); % on raw, full length
        %               COMtraj_intLR{k} = b2(1) + b2(2).*[1:numel(x)];
        %               scatter(laps,COMtraj); hold on; plot(1:numel(x),b2(1) + b2(2).*[1:length(COMtraj)])
        
                    % nonlinear regression: 1 - exponential
                    Ynl = COMtraj'; Ynl = Ynl(~isnan(Ynl));
                    Xnl = [0:length(Ynl)-1]'; Xnl = Xnl(~isnan(Ynl));
                    
                      modelfun = fittype(@(p1,p2,p3,x) p1*(1-exp(-x/p2)) + p3);
                      options = fitoptions('Method','NonlinearLeastSquares', 'Algorithm', 'Trust-Region');
                    if b2(2)>0 % positive slope from linear regression
                    params0 = [14 2 0];
                    options.Lower = [0,0.1,-25];
                    options.Upper = [200,100,25];
                    else
                    params0 = [-15 2 0];
                    options.Lower = [-200,0.1,-25];
                    options.Upper = [0,100,25];
                    end
                    options.Startpoint = params0;
                    [mdl,gof,foutput] = fit(Xnl,Ynl,modelfun,options);

                    %store in table
                    T_all.COMinterp_slope(i) = b(2);
                    T_all.COMinterp_intcp(i) = b(1);
                    T_all.COMinterp_R2(i) = stats(1);
                    T_all.COMinterp_pval(i) = stats(3);
    
                    T_all.COM_slope(i) = b2(2);
                    T_all.COM_intcp(i) = b2(1);
                    T_all.COM_R2(i) = stats2(1);
                    T_all.COM_pval(i) = stats2(3);
    
                    if b2(2) > 0 & stats2(3) < 0.05
                    T_all.shiftDir(i) = 1;
                    elseif b2(2) < 0 & stats2(3) < 0.05
                    T_all.shiftDir(i) = -1;    
                    else
                    T_all.shiftDir(i) = 0;    
                    end
    
                    T_all.RCreg_p1(i) = mdl.p1;
                    T_all.RCreg_p2(i) = mdl.p2;
                    T_all.RCreg_p3(i) = mdl.p3;
                    T_all.RCreg_R2(i) = gof.rsquare;
                    T_all.RCreg_R2adj(i) = gof.adjrsquare;
    
                % truncate trajectories to all be the same length and store them as a matrix
                    COMtraj3 = COMtraj2(1:Minlaps); %interpolated
                    COMtraj4 = COMtraj(1:Minlaps); %raw
                    COMtrajMat = [COMtrajMat ; COMtraj3]; % append trajectory to matrix (rows = PFs = PCA observations; columns = laps = variables)
                    COMtrajMatNaN = [COMtrajMatNaN ; COMtraj4];
                    
                    % track control vs opto in a vector of same length as the number of trajectories belonging to first_day1, i.e. in COMtrajMat
                    if strcmp(T_all.env(i),'opto_first_day1')
                    Opto = [Opto; 1]; 
                    else
                    Opto = [Opto; 0]; %Ctrl condition
                    end

                elseif TrajLength < Minlaps
                    SmallTraj = SmallTraj + 1; %count the number of excluded PFs
                    COMtraj3 = zeros(1,Minlaps)*NaN; % row of NaNs to exclude the trajectories with not enough laps
                    COMtraj4 = zeros(1,Minlaps)*NaN;
                    COMtrajMat = [COMtrajMat ; COMtraj3];
                    COMtrajMatNaN = [COMtrajMatNaN ; COMtraj4];

                    T_all.COM_slope(i) = NaN;                
                    T_all.COM_intcp(i) = NaN;
                    T_all.COM_R2(i) = NaN;
                    T_all.COM_pval(i) = NaN;  
                    T_all.shiftDir(i) = NaN;
                    T_all.COMinterp_slope(i) = NaN;
                    T_all.COMinterp_intcp(i) = NaN;                
                    T_all.COMinterp_R2(i) = NaN;                
                    T_all.COMinterp_pval(i) = NaN;    
                    T_all.RCreg_p1(i) = NaN;
                    T_all.RCreg_p2(i) = NaN;
                    T_all.RCreg_p3(i) = NaN;
                    T_all.RCreg_R2(i) = NaN;
                    T_all.RCreg_R2adj(i) = NaN;
                end
            else
                SmallTraj = SmallTraj + 1; %count the number of excluded PFs
                COMtraj3 = zeros(1,Minlaps)*NaN; % row of NaNs to exclude the trajectories with not enough laps
                COMtraj4 = zeros(1,Minlaps)*NaN;
                COMtrajMat = [COMtrajMat ; COMtraj3];
                COMtrajMatNaN = [COMtrajMatNaN ; COMtraj4]; 

                T_all.COM_slope(i) = NaN;                
                T_all.COM_intcp(i) = NaN;
                T_all.COM_R2(i) = NaN;
                T_all.COM_pval(i) = NaN;  
                T_all.shiftDir(i) = NaN;
                T_all.COMinterp_slope(i) = NaN;
                T_all.COMinterp_intcp(i) = NaN;                
                T_all.COMinterp_R2(i) = NaN;                
                T_all.COMinterp_pval(i) = NaN;    
                T_all.RCreg_p1(i) = NaN;
                T_all.RCreg_p2(i) = NaN;
                T_all.RCreg_p3(i) = NaN;
                T_all.RCreg_R2(i) = NaN;
                T_all.RCreg_R2adj(i) = NaN;
            end

        else
            T_all.PFopto(i) = 0;
            T_all.COM_slope(i) = NaN;                
            T_all.COM_intcp(i) = NaN;
            T_all.COM_R2(i) = NaN;
            T_all.COM_pval(i) = NaN;  
            T_all.shiftDir(i) = NaN;
            T_all.COMinterp_slope(i) = NaN;
            T_all.COMinterp_intcp(i) = NaN;                
            T_all.COMinterp_R2(i) = NaN;                
            T_all.COMinterp_pval(i) = NaN;    
            T_all.RCreg_p1(i) = NaN;
            T_all.RCreg_p2(i) = NaN;
            T_all.RCreg_p3(i) = NaN;
            T_all.RCreg_R2(i) = NaN;
            T_all.RCreg_R2adj(i) = NaN;
        end
    else
        T_all.PFopto(i) = 0;
        T_all.COM_slope(i) = NaN;                
        T_all.COM_intcp(i) = NaN;
        T_all.COM_R2(i) = NaN;
        T_all.COM_pval(i) = NaN;  
        T_all.shiftDir(i) = NaN;
        T_all.COMinterp_slope(i) = NaN;
        T_all.COMinterp_intcp(i) = NaN;                
        T_all.COMinterp_R2(i) = NaN;                
        T_all.COMinterp_pval(i) = NaN;    
        T_all.RCreg_p1(i) = NaN;
        T_all.RCreg_p2(i) = NaN;
        T_all.RCreg_p3(i) = NaN;
        T_all.RCreg_R2(i) = NaN;
        T_all.RCreg_R2adj(i) = NaN;
    end
end

COMcol1 = COMtrajMat(:,1);
OkTraj = length(COMcol1(~isnan(COMcol1)))
TotTraj = OkTraj + SmallTraj

COMtrajMat2 = COMtrajMat(~isnan(COMcol1),:); % matrix without NaNs, i.e. without the PFs with not enough laps
COMtrajMatNaN2 = COMtrajMatNaN(~isnan(COMcol1),:);
NaNtraj = isnan(COMcol1);
OKtrajIdx = ~isnan(COMcol1);

Tall = T_all(T_all.PFopto==1,:); %excludes rows not analyzed above
% Tall.COM_slope = T.COM_slope;                
% Tall.COM_intcp = T.COM_intcp;
% Tall.COM_R2 = T.COM_R2;
% Tall.COM_pval = T.COM_pval;  
% Tall.shiftDir = T.shiftDir;
% Tall.COMinterp_slope = T.COMinterp_slope;
% Tall.COMinterp_intcp = T.COMinterp_intcp;                
% Tall.COMinterp_R2 = T.COMinterp_R2;                
% Tall.COMinterp_pval = T.COMinterp_pval ;    
% Tall.RCreg_p1 = T.RCreg_p1;
% Tall.RCreg_p2 = T.RCreg_p2;
% Tall.RCreg_p3 = T.RCreg_p3;
% Tall.RCreg_R2 = T.RCreg_R2;
% Tall.RCreg_R2adj = T.RCreg_R2adj;

%writetable(T_all,'Anqi_all.csv')

Tok = Tall(~isnan(COMcol1),:); % Tok is same height as COMTrajMat2
Tok.Opto_ON = Opto;

writetable(Tok,'Toptofirst.csv')
save('COMtrajMat10laps_optofirst', 'COMtrajMat2')
%% color codes

% Opto (blue), Ctrl black
mapOpto = [0 0.5 1];
mapCtrl = [0 0 0];
colors = [mapOpto; mapCtrl];

% fwd PFs should be red, Bck PFs blue, NonSig grey or black
blue = [0 0.45 0.74]; % backward
red = [0.85 0.32 0.01]; % forward
grey = [0.5 0.5 0.5]; %nonsig

%% L and R combined, compare Opto_first_day1 with Control_first_day1

OptoIdx = find(Opto);
CtrlIdx = find(Opto==0);

% Diffusion analysis

TrajsOpto = COMtrajMat2(OptoIdx,:); 
TrajsCtrl = COMtrajMat2(CtrlIdx,:);

DispSqOpto = ( TrajsOpto-TrajsOpto(:,1) ).^2;
DispSqCtrl = ( TrajsCtrl-TrajsCtrl(:,1) ).^2;

msdCtrl = mean( DispSqCtrl , 1 );
msdOpto = mean( DispSqOpto, 1 );

DispSqCAall = {DispSqOpto; DispSqCtrl};
MSDc = [msdOpto', msdCtrl'];

figure
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

figure % Slope distributions
histogram(Tok.COM_slope(OptoIdx), 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'pdf');hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(Tok.COM_slope(CtrlIdx),'EdgeColor', mapCtrl, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'pdf'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
xline(mean(Tok.COM_slope(OptoIdx)),'Color', mapOpto, 'LineWidth', 1.5 );hold on
xline(mean(Tok.COM_slope(CtrlIdx)),'Color', mapCtrl, 'LineWidth', 1.5 );hold on
% xlim([-4, 3])
box off; axis square;
xlabel('Slope (cm/lap)'); ylabel('pdf');

figure % R2 distributions
histogram(Tok.COM_R2(OptoIdx), 'EdgeColor', mapOpto, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(Tok.COM_R2(CtrlIdx),'EdgeColor', mapCtrl, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
xline(median(Tok.COM_R2(OptoIdx)),'Color', mapOpto, 'LineWidth', 1.5 );hold on
xline(median(Tok.COM_R2(CtrlIdx)),'Color', mapCtrl, 'LineWidth', 1.5 );hold on
% xlim([-4, 3])
box off; axis square;
xlabel('lin reg R2'); ylabel('proba');

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

figure
subplot(1,2,1)% proportion of forward, backward and stable/non-signif PFs (mean across animals)
    b = bar([Prop_Bck_CTRall Prop_Fwd_CTRall Prop_NonSig_CTRall; Prop_Bck_OptoAll Prop_Fwd_OptoAll Prop_NonSig_OptoAll], 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylabel('ratio of PFs');
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['Ctrl';'Opto']);
    box off; %axis square;
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
    box off; %axis square;

%% Lin reg: delta-delta bootstrapped stats for median abs(slope)
[slopeBackCI_CA1N, slopeBackboot_CA1N] = bootci(10000,@median, abs(Tok.COM_slope(CtrlIdx)));
[slopeBackCI_CA1F, slopeBackboot_CA1F] = bootci(10000,@median, abs(Tok.COM_slope(OptoIdx)));
meanSlopeBCA1N = mean(slopeBackboot_CA1N);
meanSlopeBCA1F = mean(slopeBackboot_CA1F);

    % pairwise comparisons N vs F
    DeltaSlopeBCA1_NvsF = slopeBackboot_CA1F - slopeBackboot_CA1N; 
    DeltaSlopeBCA1_NvsF_CI = prctile(DeltaSlopeBCA1_NvsF,[2.5, 97.5]);
%     pvalBca1 = length(find(DeltaSlopeBCA1_NvsF>=0))/length(DeltaSlopeBCA1_NvsF);

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
%     xlim([-meanPropBCA1F 1-meanPropBCA1F]);
    xlim([-0.3 0.3]);
    xlabel('\Delta median |slope| (Opto - Ctrl)'); ylabel('proba');
    title('bootstrapped distrib')
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


%% Difference in Amp distrib?
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
% Combine the animal-wise bootstrapped distribution to obtain the full
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

figure % double bootstrapped paired test that MSD is different between Opto and Ctrl
    scatter(LapsPerMouse(:), DiffMSDpermouse, 'go'); hold on
    plot(0:Minlaps-1, meanDiffMSD,'-', 'LineWidth',1,'Color', 'g'); hold on 
    plot_ci(0:Minlaps-1, [meanDiffMSDboot2_estim' meanDiffMSDboot2_CI(1,:)' meanDiffMSDboot2_CI(2,:)'], 'MainLineColor', 'k', 'PatchColor', 'k', 'LineColor', 'k', 'PatchAlpha', 0.5); hold on
    plot([0 length(meanDiffMSDboot')], [0 0], 'r--');
    xlabel('Laps after onset')
    ylabel('\DeltaMSD (cm^2)')
    box off
    axis square
%     legend('average across animals', 'animal-wise bootstrapped 95% CI', 'Location', 'BestOutside')
    legend('for each animal', 'average across animals', 'animal-wise paired bootstrapped 95% CI of the mean',  'Location', 'BestOutside')

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

figure % Slope vs R2 for each file
for m = 1:length(mice)
subplot(2,4,m)
gscatter(Tok.COM_slope(CtrlIdx_mouse{m}),Tok.COM_R2(CtrlIdx_mouse{m}), Tok.shiftDir(CtrlIdx_mouse{m}), [grey; blue; red], 'o', 5, 'off')
% xlabel('lin reg slope (cm/lap)'); ylabel('lin reg R2');
% xlim([-8, 8])
% xlim([-4, 3])
title(mice{m})
box off; axis square;
end

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

%     region(m) = unique(Tok.Region(ismember(Tok.mouse, mice{m})));
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

clear 
close all

Minlaps = 15; % min number of laps where PF is defined (after interpolation, if interp is used)
Minlaps2 = 30; % min number of laps runned by animal 

% load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA1_COM.mat')
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA1_COM_withPFs.mat')
COMdata.CA1 = data; clear data
COMdata.CA1.animalID = {'1-3'; '1-3'; '1-3'; '1-3'; 'cdc'; 'cfc'; 'cfc'; 'wt1'; 'wt1'};

% load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA3_COM_new_namecorrected.mat')
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA3_COM_withPFs.mat')
COMdata.CA3 = data; clear data
COMdata.CA3.animalID = {'4-2'; '4-2'; '4-2'; '4-2'; '4-2'; '4-2'; '4-1'; '4-1'; '4-1'; '4-1'; '4-1'; '4-1'; '5-3'; '5-3'; '5-3'; '5-1'; '5-1'; '5-1'; '5-3'; '5-3'; '5-3';'5-4';'5-4';'5-4'; '5-1'; '5-1'; '5-1';'5-4';'5-4';'5-4';'7-2';'7-2';'7-2';'7-2';'9-1';'9-1';'9-1' };

fn=fieldnames(COMdata);

SpeedT = readtable('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\IndividualBehavior\All\Can1_SpeedWithFreeze.xlsx');
% note: VR: 1 is f1, 2 is f2, 3 is n1, 4 is n2

COMtrajMat = [];
COMtrajMatNaN = [];
T.COMonset = [];
Region = {};
VR = {};
TotNumPFs = 0;
SmallTraj = 0;
Allshifts = [];

T.animalID = {};
T.Region = {}; 
T.PlaneID = {};
T.f_n = {};
T.VRstr = {};
T.VR = [];
T.PF = [];
T.animalSpeed = [];  
T.meanShift =  [];
T.maxShift = [];
T.maxShift_Lap = [];
T.maxS_LapSpeed = [];
T.maxS_PrevLapSpeed = [];
T.COM_slope = [];                
T.COM_intcp = [];
T.COM_R2 = [];
T.COM_pval = [];  
T.shiftDir = {};
T.COMinterp_slope = [];
T.COMinterp_intcp = [];                
T.COMinterp_R2 = [];                
T.COMinterp_pval = [];    
T.RCreg_p1 = [];
T.RCreg_p2 = [];
T.RCreg_p3 = [];
T.RCreg_R2 = [];
T.RCreg_R2adj = [];
% T.RCreg_pval = [];

%loop through the fields
for a = 1:2 %CA1 or CA3
    fn1 = fieldnames(COMdata.(fn{a})); 
    for i = 1:2 % VR: f or n
        fn2 = fieldnames(COMdata.(fn{a}).(fn1{i}));
        DcomAll{a,i} = [];
        for j = 1:length(fn2) % imaging planes
            
%             SpatialFRout(lap,:) = SpikeCountOut(lap, :)./TimeInBins(lap, :); % FR in each spatial bin, in Hz 
%             %Compute lap-wise COM, SD, skewness, mean and max firing rate
%             COMbin(lap) = sum(SpatialFRout(lap,:).*[1:Nbin])/sum(SpatialFRout(lap,:));
%             COMloc(lap) = sum(SpatialFRout(lap,:).*BinCenters)/sum(SpatialFRout(lap,:));
%             PFsdOut(lap) = sqrt( sum( (BinCenters - COMloc(lap)).^2.*SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ) ); %lap-wise SD of the PF

            COMloc = COMdata.(fn{a}).(fn1{i}).(fn2{j}).COM .*6; % multiply by 6 to get location in cm (this is inherited from Can's code)
            Onset = COMdata.(fn{a}).(fn1{i}).(fn2{j}).onset_lap;
            NumPFs = length(Onset);
            TotNumPFs = TotNumPFs + NumPFs; % track total number of PFs;
            NumLaps = size(COMloc,2); % number of laps during imaging session
            Dcom = diff(COMloc,1,2); % delta COM (in cm) from one lap to the next, for all PFs of the imaging plane
            DcomLaps = size(Dcom,2);

            if DcomLaps >= Minlaps2-1
               DcomAll{a,i} = [DcomAll{a,i}; Dcom(:,1:Minlaps2-1)]; %truncate Dcom and concatenate with previous planes
               Allshifts = [Allshifts; Dcom(:,1:Minlaps2-1)]; % concatenate for all conditions
            end
           
            for k = 1:NumPFs % PFs in the plane           

                % extract COM trajectory, starting from onset lap, and centered to initial COM location
                COMtraj = COMloc(k,Onset(k):end) - COMloc(k,Onset(k)); 
%                 COMtraj = COMloc(k,Onset(k):Onset(k)+Minlaps-1) - COMloc(k,Onset(k));

                % Interpolate NaNs
                x = COMtraj; nanx = isnan(x); laps = 1:numel(x);
                x(nanx) = interp1(laps(~nanx), x(~nanx), laps(nanx)); % Linear extrapolation. No extrapolation for NaNs at the end
%                 x(nanx) = interp1(laps(~nanx), x(~nanx), laps(nanx), 'previous'); 
                COMtraj2 = x(~isnan(x)); % get rid of end NaNs

                TrajLength = length(COMtraj2);

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

%                 options = statset('FunValCheck','off');
%                 modelfun = 'Ynl ~ b1*(1-exp(-Xnl/b2)) + b3';
%                 mdl = fitnlm(Xnl,Ynl,modelfun,params0,'Options', options);
               
%                 figure
%                 plot(Xnl, feval(mdl,Xnl), 'k'); hold on
% %                 plot(Xnl, mdl.Fitted, 'r'); hold on
%                 scatter(Xnl, Ynl, 'ko')

                % truncate trajectories to all be the same length and store them as a matrix
                if TrajLength >= Minlaps
                    COMtraj3 = COMtraj2(1:Minlaps); %interpolated
                    COMtraj4 = COMtraj(1:Minlaps); %raw
                    COMtrajMat = [COMtrajMat ; COMtraj3]; % append trajectory to matrix (rows = PFs = PCA observations; columns = laps = variables)
                    COMtrajMatNaN = [COMtrajMatNaN ; COMtraj4];
                    % COMonset = [COMonset; Onset(k)]; %
                    % Region = [Region; fn{a}]; % CA1 or CA3
                    % VR = [VR; fn1{i}]; % F or N
                elseif TrajLength < Minlaps
                    SmallTraj = SmallTraj + 1; %count the number of excluded PFs
                    COMtraj3 = zeros(1,Minlaps)*NaN; % row of NaNs to exclude the trajectories with not enough laps
                    COMtraj4 = zeros(1,Minlaps)*NaN;
                    COMtrajMat = [COMtrajMat ; COMtraj3];
                    COMtrajMatNaN = [COMtrajMatNaN ; COMtraj4];
                end

                % prepare a table to store different values for each PF 
                T.Region{end+1} = fn{a};
                T.f_n{end+1} = fn1{i};
                T.animalID{end+1} = COMdata.(fn{a}).animalID{j};
                T.PF(end+1) = k;
                T.COMonset(end+1) = Onset(k);

                % parse info from namefile
                S = strfind(fn2(j), 'plain');
                T.PlaneID{end+1} = fn2{j}(S{1}+5);
                VRnum = extract(fn2(j),'00'+digitsPattern); 
                T.VRstr{end+1}= [T.f_n{end} VRnum{end}]; 
                
                % scalar code as in the Speed table
                if contains(T.VRstr{end}, 'f001') | contains(T.VRstr{end}, 'f000')
                  T.VR(end+1) = 1;
                elseif contains(T.VRstr{end}, 'f002')
                    T.VR(end+1) = 2;
                elseif contains(T.VRstr{end}, 'n001') | contains(T.VRstr{end}, 'n000')
                    T.VR(end+1) = 3;
                elseif contains(T.VRstr{end},'n002')
                    T.VR(end+1) = 4;
                else
%                     disp('error')
%                     return
                    T.VR(end+1) = NaN;
                end

                % speed
                T.animalSpeed(end+1) = mean(SpeedT.LapSpeed(strcmp(T.animalID{end}, SpeedT.animalID) & SpeedT.VR == T.VR(end))); %average animal speed over the laps considered
                
                T.meanShift =  [T.meanShift; mean(Dcom(k,:), 'omitnan')];
                [maxShift, maxShift_Lap] = max(abs(Dcom(k,:)));
                T.maxShift = [T.maxShift; maxShift];
                T.maxShift_Lap = [T.maxShift_Lap; maxShift_Lap];
                maxS_LapSpeed = SpeedT.LapSpeed(strcmp(T.animalID{end}, SpeedT.animalID) & SpeedT.VR == T.VR(end) & SpeedT.LapNum == T.maxShift_Lap(end));
                maxS_PrevLapSpeed = SpeedT.LapSpeed(strcmp(T.animalID{end}, SpeedT.animalID) & SpeedT.VR == T.VR(end) & SpeedT.LapNum == T.maxShift_Lap(end)-1);
                if ~isempty(maxS_LapSpeed)
                T.maxS_LapSpeed(end+1) = maxS_LapSpeed;
                else
                T.maxS_LapSpeed(end+1) = NaN;    
                end
                if ~isempty(maxS_PrevLapSpeed)
                T.maxS_PrevLapSpeed(end+1) = maxS_PrevLapSpeed;
                else
                T.maxS_PrevLapSpeed(end+1) = NaN;    
                end
                clear maxS_LapSpeed maxS_PrevLapSpeed
        
                T.COMinterp_slope = [T.COMinterp_slope; b(2)];
                T.COMinterp_intcp = [T.COMinterp_intcp; b(1)];
                T.COMinterp_R2 = [T.COMinterp_R2; stats(1)];
                T.COMinterp_pval = [T.COMinterp_pval; stats(3)];

                T.COM_slope = [T.COM_slope; b2(2)];
                T.COM_intcp = [T.COM_intcp; b2(1)];
                T.COM_R2 = [T.COM_R2; stats2(1)];
                T.COM_pval = [T.COM_pval; stats2(3)];

                if b2(2) > 0 & stats2(3) < 0.05
                T.shiftDir{end+1} = 'F';
                elseif b2(2) < 0 & stats2(3) < 0.05
                T.shiftDir{end+1} = 'B';    
                else
                T.shiftDir{end+1} = 'N';    
                end

                T.RCreg_p1 = [T.RCreg_p1; mdl.p1];
                T.RCreg_p2 = [T.RCreg_p2; mdl.p2];
                T.RCreg_p3 = [T.RCreg_p3; mdl.p3];
                T.RCreg_R2 = [T.RCreg_R2; gof.rsquare];
                T.RCreg_R2adj = [T.RCreg_R2adj; gof.adjrsquare];

%                 T.RCreg_p1 = [T.RCreg_p1; mdl.Coefficients.Estimate(1)];
%                 T.RCreg_p2 = [T.RCreg_p2; mdl.Coefficients.Estimate(2)];
%                 T.RCreg_R2 = [T.RCreg_R2; mdl.Rsquared.Ordinary];
%                 T.RCreg_R2adj = [T.RCreg_R2adj; mdl.Rsquared.Adjusted];
% %                 T.RCreg_pval = [T.RCreg_pval; mdl.ModelFitVsNullModel.Fstats]; %field doesn't seem to exist despite what says documentation...

% 
%                 I = find(strcmp(BehavM.uniqueID, x1.AnimalID)>0);
%                           if ~isempty(I)

%                 figure % COM trajectory (raw + interpolated) and regression
%                 % plot(1:Nlaps, COMloc, 'k-'); hold on
%                 plot(1:length(COMtraj2), COMtraj2, 'b-'); hold on
%                 scatter(1:length(COMtraj), COMtraj, 'r'); hold on
%                 yline(0,'k--');
%                 ylim([-150 150]);
%                 xlabel('lap'); ylabel('COM position (cm)');
%                 % title({'Place Field COM trajectory, slope =' num2str(b(2)) 'pval =' num2str(stats(3))})
            
            end
        end
    end
end
COMcol1 = COMtrajMat(:,1);
OkTraj = length(COMcol1(~isnan(COMcol1)))
TotTraj = OkTraj + SmallTraj
TotNumPFs

COMtrajMat2 = COMtrajMat(~isnan(COMcol1),:); % matrix without NaNs, i.e. without the PFs with not enough laps
COMtrajMatNaN2 = COMtrajMatNaN(~isnan(COMcol1),:);
NaNtraj = isnan(COMcol1);
OKtrajIdx = ~isnan(COMcol1);

% make table
T.animalID = T.animalID';
T.Region = T.Region'; 
T.PlaneID = T.PlaneID';
T.f_n = T.f_n';
T.VRstr = T.VRstr';
T.VR = T.VR';
T.PF = T.PF';
T.COMonset = T.COMonset';
T.animalSpeed = T.animalSpeed';
T.maxS_LapSpeed = T.maxS_LapSpeed';
T.maxS_PrevLapSpeed = T.maxS_PrevLapSpeed';
T.shiftDir = T.shiftDir';

% T = rmfield(T,'RCreg_pval');
Tall = struct2table(T);
% filename = 'Can1All.xlsx';
% writetable(Can1All,filename);

% Groups per conditions, excluding NaN trajectories, i.e. PFs with fewer laps than Minlaps
Tok = Tall(~isnan(COMcol1),:);
CA1_F_idx = find(ismember(Tok.Region,'CA1') & ismember(Tok.f_n,'f'));
CA3_F_idx = find(ismember(Tok.Region,'CA3') & ismember(Tok.f_n,'f'));
CA1_N_idx = find(ismember(Tok.Region,'CA1') & ismember(Tok.f_n,'n'));
CA3_N_idx = find(ismember(Tok.Region,'CA3') & ismember(Tok.f_n,'n'));
CA_VR_idx = {CA1_N_idx; CA1_F_idx; CA3_N_idx; CA3_F_idx};
GroupCA_VR(CA1_N_idx) = 1;
GroupCA_VR(CA1_F_idx) = 2;
GroupCA_VR(CA3_N_idx) = 3;
GroupCA_VR(CA3_F_idx) = 4;

% Groups per conditions, All PFs
CA1_F_idxALL = find(ismember(T.Region,'CA1') & ismember(T.f_n,'f'));
CA3_F_idxALL = find(ismember(T.Region,'CA3') & ismember(T.f_n,'f'));
CA1_N_idxALL = find(ismember(T.Region,'CA1') & ismember(T.f_n,'n'));
CA3_N_idxALL = find(ismember(T.Region,'CA3') & ismember(T.f_n,'n'));
CA_VR_idxALL = {CA1_N_idxALL; CA1_F_idxALL; CA3_N_idxALL; CA3_F_idxALL};
GroupCA_VRall(CA1_N_idxALL) = 1;
GroupCA_VRall(CA1_F_idxALL) = 2;
GroupCA_VRall(CA3_N_idxALL) = 3;
GroupCA_VRall(CA3_F_idxALL) = 4;

% mean trajectories
COMmean = mean(COMtrajMat2,2); %size(COMmean);
BackwardTraj = mean(COMtrajMat2(COMmean<0,:),1); 
ForwardTraj = mean(COMtrajMat2(COMmean>0,:),1);

%colorcode CA1 vs CA3, N (dark) vs F (light); e.g. CA1F -> mapCA1(1,:)
map = brewermap(10,'Paired'); %cool(4); %brewermap(4,'PiYG')
mapCA1 = map(3:4,:); %brewermap(2,'GnBu');
mapCA3 = map(7:8,:); %brewermap(2,'OrRd');
mapCA = [mapCA1;mapCA3];
colors = [mapCA1(2,:); mapCA1(1,:); mapCA3(2,:); mapCA3(1,:)];
figure; scatter(1:4,1:4,200,mapCA, 'filled')

% fwd PFs should be red, Bck PFs blue, NonSig grey or black
blue = [0 0.45 0.74]; % backward
red = [0.85 0.32 0.01]; % forward
grey = [0.5 0.5 0.5]; %nonsig

%% linear regression

Fwd_idx = find(Tok.COM_slope > 0 & Tok.COM_pval<0.05);
Bwd_idx = find(Tok.COM_slope < 0 & Tok.COM_pval<0.05);
Nsig_idx = find(Tok.COM_pval>=0.05);

%% nonlinear regression
% delta-delta bootstrapped stats for median abs(AMp)
[AmpCI_CA1N, ampBackboot_CA1N] = bootci(10000,@median, abs(Tok.RCreg_p1(CA1_N_idx)));
[AmpCI_CA1F, ampBackboot_CA1F] = bootci(10000,@median, abs(Tok.RCreg_p1(CA1_F_idx)));
[AmpCI_CA3N, ampBackboot_CA3N] = bootci(10000,@median, abs(Tok.RCreg_p1(CA3_N_idx)));
[AmpCI_CA3F, ampBackboot_CA3F] = bootci(10000,@median, abs(Tok.RCreg_p1(CA3_F_idx)));
meanAmpCA1N = mean(ampBackboot_CA1N);
meanAmpCA1F = mean(ampBackboot_CA1F);
meanAmpCA3N = mean(ampBackboot_CA3N);
meanAmpCA3F = mean(ampBackboot_CA3F);

    % pairwise comparisons N vs F
    DeltaAmpCA1_NvsF = ampBackboot_CA1F - ampBackboot_CA1N; 
    DeltaAmpCA1_NvsF_CI = prctile(DeltaAmpCA1_NvsF,[2.5, 97.5]);
%     pvalBca1 = length(find(DeltaSlopeBCA1_NvsF>=0))/length(DeltaSlopeBCA1_NvsF);
    DeltaAmpCA3_NvsF = ampBackboot_CA3F - ampBackboot_CA3N; 
    DeltaAmpCA3_NvsF_CI = prctile(DeltaAmpCA3_NvsF,[2.5, 97.5]);
    
    % interaction between CA and N_F factors
    DeltaDeltaAmp = DeltaAmpCA1_NvsF - DeltaAmpCA3_NvsF;
    DeltaDeltaAmp_CI = prctile(DeltaDeltaAmp,[2.5, 97.5]);
    pvalAmp = length(find(DeltaDeltaAmp>=0))/length(DeltaDeltaAmp);

figure % all PFs pooled. median |Amp|
subplot(2,2,1) % abs(Amp) (all PFs pooled)
    v3 = violinplot(abs(Tok.RCreg_p1), GroupCA_VR);
    v3(1).ViolinColor = mapCA1(2,:);
    v3(2).ViolinColor = mapCA1(1,:);
    v3(3).ViolinColor = mapCA3(2,:);
    v3(4).ViolinColor = mapCA3(1,:);
    ylabel('|Amp|');
    set(gca, 'XTick', [1,2,3,4], 'XTickLabel', ['CA1N';'CA1F'; 'CA3N'; 'CA3F']);
    box off; axis square;
subplot(2,2,2)
    plot([1,2], [meanAmpCA1N meanAmpCA1F], 'k-'); hold on
    plot([1,2], [meanAmpCA3N meanAmpCA3F], 'k-'); hold on
    errorbar(2 , meanAmpCA1F, meanAmpCA1F-AmpCI_CA1F(1), meanAmpCA1F-AmpCI_CA1F(2), 'o', 'Color', mapCA1(1,:), 'MarkerFaceColor', mapCA1(1,:), 'LineWidth', 1); hold on
    errorbar(1 , meanAmpCA1N, meanAmpCA1N-AmpCI_CA1N(1), meanAmpCA1N-AmpCI_CA1N(2), 'o', 'Color', mapCA1(2,:), 'MarkerFaceColor', mapCA1(2,:), 'LineWidth', 1); hold on
    errorbar(2 , meanAmpCA3F, meanAmpCA3F-AmpCI_CA3F(1), meanAmpCA3F-AmpCI_CA3F(2), 'o', 'Color', mapCA3(1,:), 'MarkerFaceColor', mapCA3(1,:), 'LineWidth', 1); hold on
    errorbar(1 , meanAmpCA3N, meanAmpCA3N-AmpCI_CA3N(1), meanAmpCA3N-AmpCI_CA3N(2), 'o', 'Color', mapCA3(2,:), 'MarkerFaceColor', mapCA3(2,:), 'LineWidth', 1); hold on
    ylabel('median |Amp|')
    xlim([0 3]); ylim([5 27])
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
    box off; axis square;
    % title('bootstrapped 95% CIs')
subplot(2,2,3)% CA1 N vs F bootstrapped delta |slope|
    histogram(DeltaAmpCA1_NvsF, 'EdgeColor', mapCA1(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaAmpCA1_NvsF), 0.01, mean(DeltaAmpCA1_NvsF)-DeltaAmpCA1_NvsF_CI(1),mean(DeltaAmpCA1_NvsF)-DeltaAmpCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 6); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-17 0.5]);
    xlabel('\Delta median |Amp| (F - N)'); ylabel('proba');
    % title('bootstrapped distrib')
    set(gca,'view',[90 -90])
    box off; axis square;
subplot(2,2,4) % CA3 N vs F bootstrapped delta |slope|
    histogram(DeltaAmpCA3_NvsF, 'EdgeColor', mapCA3(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaAmpCA3_NvsF), 0.01, mean(DeltaAmpCA3_NvsF)-DeltaAmpCA3_NvsF_CI(1), mean(DeltaAmpCA3_NvsF)-DeltaAmpCA3_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 6); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-17 0.5]);
    xlabel('\Delta median |Amp| (F - N)'); ylabel('proba');
    % title('bootstrapped distrib')
    set(gca,'view',[90 -90])
    box off; axis square;

% delta-delta bootstrapped stats for median Tau)
[TauCI_CA1N, Tauboot_CA1N] = bootci(10000,@median, abs(Tok.RCreg_p2(CA1_N_idx)));
[TauCI_CA1F, Tauboot_CA1F] = bootci(10000,@median, abs(Tok.RCreg_p2(CA1_F_idx)));
[TauCI_CA3N, Tauboot_CA3N] = bootci(10000,@median, abs(Tok.RCreg_p2(CA3_N_idx)));
[TauCI_CA3F, Tauboot_CA3F] = bootci(10000,@median, abs(Tok.RCreg_p2(CA3_F_idx)));
meanTauCA1N = mean(Tauboot_CA1N);
meanTauCA1F = mean(Tauboot_CA1F);
meanTauCA3N = mean(Tauboot_CA3N);
meanTauCA3F = mean(Tauboot_CA3F);

    % pairwise comparisons N vs F
    DeltaTauCA1_NvsF = Tauboot_CA1F - Tauboot_CA1N; 
    DeltaTauCA1_NvsF_CI = prctile(DeltaTauCA1_NvsF,[2.5, 97.5]);
%     pvalBca1 = length(find(DeltaSlopeBCA1_NvsF>=0))/length(DeltaSlopeBCA1_NvsF);
    DeltaTauCA3_NvsF = Tauboot_CA3F - Tauboot_CA3N; 
    DeltaTauCA3_NvsF_CI = prctile(DeltaTauCA3_NvsF,[2.5, 97.5]);
    
    % interaction between CA and N_F factors
    DeltaDeltaTau = DeltaTauCA1_NvsF - DeltaTauCA3_NvsF;
    DeltaDeltaTau_CI = prctile(DeltaDeltaTau,[2.5, 97.5]);
    pvalTau = length(find(DeltaDeltaTau>=0))/length(DeltaDeltaTau);

figure % all PFs pooled. Bar plots of props for all Ca and conditions. Delta F vs N for CA1 and CA3, oriented like in Dabest. 
subplot(2,2,1) % Tau (all PFs pooled)
    v3 = violinplot(abs(Tok.RCreg_p2), GroupCA_VR);
    v3(1).ViolinColor = mapCA1(2,:);
    v3(2).ViolinColor = mapCA1(1,:);
    v3(3).ViolinColor = mapCA3(2,:);
    v3(4).ViolinColor = mapCA3(1,:);
    ylabel('Tau');
    set(gca, 'XTick', [1,2,3,4], 'XTickLabel', ['CA1N';'CA1F'; 'CA3N'; 'CA3F']);
    box off; %axis square;
subplot(2,2,2)
    plot([1,2], [meanTauCA1N meanTauCA1F], 'k-'); hold on
    plot([1,2], [meanTauCA3N meanTauCA3F], 'k-'); hold on
    errorbar(2 , meanTauCA1F, meanTauCA1F-TauCI_CA1F(1), meanTauCA1F-TauCI_CA1F(2), 'o', 'Color', mapCA1(1,:), 'MarkerFaceColor', mapCA1(1,:)); hold on
    errorbar(1 , meanTauCA1N, meanTauCA1N-TauCI_CA1N(1), meanTauCA1N-TauCI_CA1N(2), 'o', 'Color', mapCA1(2,:), 'MarkerFaceColor', mapCA1(2,:)); hold on
    errorbar(2 , meanTauCA3F, meanTauCA3F-TauCI_CA3F(1), meanTauCA3F-TauCI_CA3F(2), 'o', 'Color', mapCA3(1,:), 'MarkerFaceColor', mapCA3(1,:)); hold on
    errorbar(1 , meanTauCA3N, meanTauCA3N-TauCI_CA3N(1), meanTauCA3N-TauCI_CA3N(2), 'o', 'Color', mapCA3(2,:), 'MarkerFaceColor', mapCA3(2,:)); hold on
    ylabel('median Tau')
    xlim([0 3]); ylim([0 15])
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
    box off; %axis square;
    title('bootstrapped 95% CIs')
subplot(2,2,3) % CA1 N vs F bootstrapped delta |slope|
    histogram(DeltaTauCA1_NvsF, 'EdgeColor', mapCA1(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaTauCA1_NvsF), 0.01, mean(DeltaTauCA1_NvsF)-DeltaTauCA1_NvsF_CI(1),mean(DeltaTauCA1_NvsF)-DeltaTauCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-6 6]);
    xlabel('\Delta median Tau (F - N)'); ylabel('proba');
    title('bootstrapped distrib')
    set(gca,'view',[90 -90])
    box off; %axis square;
subplot(2,2,4) % CA3 N vs F bootstrapped delta |slope|
    histogram(DeltaTauCA3_NvsF, 'EdgeColor', mapCA3(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaTauCA3_NvsF), 0.01, mean(DeltaTauCA3_NvsF)-DeltaTauCA3_NvsF_CI(1), mean(DeltaTauCA3_NvsF)-DeltaTauCA3_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-6 6]);
    xlabel('\Delta median Tau (F - N)'); ylabel('proba');
    title('bootstrapped distrib')
    set(gca,'view',[90 -90])
    box off;

    % figure % all PFs pooled. median |Amp|
tl =  tiledlayout(1,5, 'TileSpacing','Compact','Padding','Compact');
% nexttile([1,2]) % abs(Amp) (all PFs pooled)
%     v3 = violinplot(abs(Tok.RCreg_p1), GroupCA_VR);
%     v3(1).ViolinColor = mapCA1(2,:);
%     v3(2).ViolinColor = mapCA1(1,:);
%     v3(3).ViolinColor = mapCA3(2,:);
%     v3(4).ViolinColor = mapCA3(1,:);
%     ylabel('|Amp|');
%     set(gca, 'XTick', [1,2,3,4], 'XTickLabel', ['CA1N';'CA1F'; 'CA3N'; 'CA3F']);
%     box off; %axis square;
% nexttile
%     plot([1,2], [meanAmpCA1N meanAmpCA1F], 'k-'); hold on
%     plot([1,2], [meanAmpCA3N meanAmpCA3F], 'k-'); hold on
%     errorbar(2 , meanAmpCA1F, meanAmpCA1F-AmpCI_CA1F(1), meanAmpCA1F-AmpCI_CA1F(2), 'o', 'Color', mapCA1(1,:), 'MarkerFaceColor', mapCA1(1,:), 'LineWidth', 1); hold on
%     errorbar(1 , meanAmpCA1N, meanAmpCA1N-AmpCI_CA1N(1), meanAmpCA1N-AmpCI_CA1N(2), 'o', 'Color', mapCA1(2,:), 'MarkerFaceColor', mapCA1(2,:), 'LineWidth', 1); hold on
%     errorbar(2 , meanAmpCA3F, meanAmpCA3F-AmpCI_CA3F(1), meanAmpCA3F-AmpCI_CA3F(2), 'o', 'Color', mapCA3(1,:), 'MarkerFaceColor', mapCA3(1,:), 'LineWidth', 1); hold on
%     errorbar(1 , meanAmpCA3N, meanAmpCA3N-AmpCI_CA3N(1), meanAmpCA3N-AmpCI_CA3N(2), 'o', 'Color', mapCA3(2,:), 'MarkerFaceColor', mapCA3(2,:), 'LineWidth', 1); hold on
%     ylabel('median |Amp|')
%     xlim([0 3]); ylim([5 27])
%     set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
%     box off; %axis square;
%     % title('bootstrapped 95% CIs')
% nexttile % CA1 N vs F bootstrapped delta |slope|
%     histogram(DeltaAmpCA1_NvsF, 'EdgeColor', mapCA1(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
%     errorbar(mean(DeltaAmpCA1_NvsF), 0.01, mean(DeltaAmpCA1_NvsF)-DeltaAmpCA1_NvsF_CI(1),mean(DeltaAmpCA1_NvsF)-DeltaAmpCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 6); hold on
%     xline(0,'k--', 'LineWidth', 1.5 );hold on
%     xlim([-17 0.5]);
%     xlabel('\Delta median |Amp| (F - N)'); ylabel('proba');
%     % title('bootstrapped distrib')
%     set(gca,'view',[90 -90])
%     box off; %axis square;
% nexttile % CA3 N vs F bootstrapped delta |slope|
%     histogram(DeltaAmpCA3_NvsF, 'EdgeColor', mapCA3(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
%     errorbar(mean(DeltaAmpCA3_NvsF), 0.01, mean(DeltaAmpCA3_NvsF)-DeltaAmpCA3_NvsF_CI(1), mean(DeltaAmpCA3_NvsF)-DeltaAmpCA3_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 6); hold on
%     xline(0,'k--', 'LineWidth', 1.5 );hold on
%     xlim([-17 0.5]);
%     xlabel('\Delta median |Amp| (F - N)'); ylabel('proba');
%     % title('bootstrapped distrib')
%     set(gca,'view',[90 -90])
%     box off; %axis square;

nexttile([1,2]) % Tau (all PFs pooled)
    v3 = violinplot(abs(Tok.RCreg_p2), GroupCA_VR);
    v3(1).ViolinColor = mapCA1(2,:);
    v3(2).ViolinColor = mapCA1(1,:);
    v3(3).ViolinColor = mapCA3(2,:);
    v3(4).ViolinColor = mapCA3(1,:);
    ylabel('Tau');
    set(gca, 'XTick', [1,2,3,4], 'XTickLabel', ['CA1N';'CA1F'; 'CA3N'; 'CA3F']);
    box off; %axis square;
nexttile
    plot([1,2], [meanTauCA1N meanTauCA1F], 'k-'); hold on
    plot([1,2], [meanTauCA3N meanTauCA3F], 'k-'); hold on
    errorbar(2 , meanTauCA1F, meanTauCA1F-TauCI_CA1F(1), meanTauCA1F-TauCI_CA1F(2), 'o', 'Color', mapCA1(1,:), 'MarkerFaceColor', mapCA1(1,:), 'LineWidth', 1); hold on
    errorbar(1 , meanTauCA1N, meanTauCA1N-TauCI_CA1N(1), meanTauCA1N-TauCI_CA1N(2), 'o', 'Color', mapCA1(2,:), 'MarkerFaceColor', mapCA1(2,:), 'LineWidth', 1); hold on
    errorbar(2 , meanTauCA3F, meanTauCA3F-TauCI_CA3F(1), meanTauCA3F-TauCI_CA3F(2), 'o', 'Color', mapCA3(1,:), 'MarkerFaceColor', mapCA3(1,:), 'LineWidth', 1); hold on
    errorbar(1 , meanTauCA3N, meanTauCA3N-TauCI_CA3N(1), meanTauCA3N-TauCI_CA3N(2), 'o', 'Color', mapCA3(2,:), 'MarkerFaceColor', mapCA3(2,:), 'LineWidth', 1); hold on
    ylabel('median Tau')
    xlim([0 3]); ylim([0 12])
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
    box off; %axis square;
    % title('bootstrapped 95% CIs')
nexttile% CA1 N vs F bootstrapped delta |slope|
    histogram(DeltaTauCA1_NvsF, 'EdgeColor', mapCA1(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaTauCA1_NvsF), 0.01, mean(DeltaTauCA1_NvsF)-DeltaTauCA1_NvsF_CI(1),mean(DeltaTauCA1_NvsF)-DeltaTauCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 6); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-6 6]);
    xlabel('\Delta median Tau (F - N)'); ylabel('proba');
    % title('bootstrapped distrib')
    set(gca,'view',[90 -90])
    box off; %axis square;
nexttile % CA3 N vs F bootstrapped delta |slope|
    histogram(DeltaTauCA3_NvsF, 'EdgeColor', mapCA3(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaTauCA3_NvsF), 0.01, mean(DeltaTauCA3_NvsF)-DeltaTauCA3_NvsF_CI(1), mean(DeltaTauCA3_NvsF)-DeltaTauCA3_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 6); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-6 6]);
    xlabel('\Delta median Tau (F - N)'); ylabel('proba');
    % title('bootstrapped distrib')
    set(gca,'view',[90 -90])
    box off;

%% Proportion of B, F and N shifting PFs considering nonlinear regression
% with N defined as Rsquare<threshold

% unpaired delta-delta plot comparing CA fields and familiarity levels

[propshiftCI_CA1N, propshiftboot_CA1N] = bootci(10000,@propshiftNL, Tok.RCreg_R2(CA1_N_idx));
[propshiftCI_CA1F, propshiftboot_CA1F] = bootci(10000,@propshiftNL, Tok.RCreg_R2(CA1_F_idx));
[propshiftCI_CA3N, propshiftboot_CA3N] = bootci(10000,@propshiftNL, Tok.RCreg_R2(CA3_N_idx));
[propshiftCI_CA3F, propshiftboot_CA3F] = bootci(10000,@propshiftNL, Tok.RCreg_R2(CA3_F_idx));
meanPropCA1N = mean(propshiftboot_CA1N);
meanPropCA1F = mean(propshiftboot_CA1F);
meanPropCA3N = mean(propshiftboot_CA3N);
meanPropCA3F = mean(propshiftboot_CA3F);

    % pairwise comparisons N vs F
    DeltaPropCA1_NvsF = propshiftboot_CA1F - propshiftboot_CA1N; 
    DeltaPropCA1_NvsF_CI = prctile(DeltaPropCA1_NvsF,[2.5, 97.5]);
    pvalCA1_NvsF = length(find(DeltaPropCA1_NvsF<=0))/length(DeltaPropCA1_NvsF);
    DeltaPropCA3_NvsF = propshiftboot_CA3F - propshiftboot_CA3N; 
    DeltaPropCA3_NvsF_CI = prctile(DeltaPropCA3_NvsF,[2.5, 97.5]);
    pvalCA3_NvsF = length(find(DeltaPropCA3_NvsF<=0))/length(DeltaPropCA3_NvsF);

    % interaction between CA and N_F factors
    DeltaDelta = DeltaPropCA1_NvsF - DeltaPropCA3_NvsF;
    DeltaDelta_CI = prctile(DeltaDelta,[2.5, 97.5]);

    % CA1 vs CA3 main effect -> different approach than the subsampling.
    CA1idx = [CA1_F_idx; CA1_N_idx ];
    CA3idx = [CA3_F_idx; CA3_N_idx ];
    [propshiftCI_CA1, propshiftboot_CA1] = bootci(10000,@propshiftNL, Tok.COM_pval(CA1idx));
    [propshiftCI_CA3, propshiftboot_CA3] = bootci(10000,@propshiftNL, Tok.COM_pval(CA3idx));
    Delta_Prop_CA1vsCA3 = propshiftboot_CA1 - propshiftboot_CA3;

    % N vs F main effect

% delta-delta bootstrapped stats for Backward shifting PFs
[propBackCI_CA1N, propBackboot_CA1N] = bootci(10000,@propBackNL, Tok.RCreg_R2(CA1_N_idx),Tok.RCreg_p1(CA1_N_idx));
[propBackCI_CA1F, propBackboot_CA1F] = bootci(10000,@propBackNL, Tok.RCreg_R2(CA1_F_idx),Tok.RCreg_p1(CA1_F_idx));
[propBackCI_CA3N, propBackboot_CA3N] = bootci(10000,@propBackNL, Tok.RCreg_R2(CA3_N_idx),Tok.RCreg_p1(CA3_N_idx));
[propBackCI_CA3F, propBackboot_CA3F] = bootci(10000,@propBackNL, Tok.RCreg_R2(CA3_F_idx),Tok.RCreg_p1(CA3_F_idx));
meanPropBCA1N = mean(propBackboot_CA1N);
meanPropBCA1F = mean(propBackboot_CA1F);
meanPropBCA3N = mean(propBackboot_CA3N);
meanPropBCA3F = mean(propBackboot_CA3F);

    % pairwise comparisons N vs F
    DeltaPropBCA1_NvsF = propBackboot_CA1F - propBackboot_CA1N; 
    DeltaPropBCA1_NvsF_CI = prctile(DeltaPropBCA1_NvsF,[2.5, 97.5]);
    pvalBca1 = length(find(DeltaPropBCA1_NvsF>=0))/length(DeltaPropBCA1_NvsF);
    DeltaPropBCA3_NvsF = propBackboot_CA3F - propBackboot_CA3N; 
    DeltaPropBCA3_NvsF_CI = prctile(DeltaPropBCA3_NvsF,[2.5, 97.5]);
    
    % interaction between CA and N_F factors
    DeltaDeltaB = DeltaPropBCA1_NvsF - DeltaPropBCA3_NvsF;
    DeltaDeltaB_CI = prctile(DeltaDeltaB,[2.5, 97.5]);
    pvalB = length(find(DeltaDeltaB>=0))/length(DeltaDeltaB);

% delta-delta bootstrapped stats for Forward shifting PFs
[propFwdCI_CA1N, propFwdboot_CA1N] = bootci(10000,@propFwdNL, Tok.RCreg_R2(CA1_N_idx),Tok.RCreg_p1(CA1_N_idx));
[propFwdCI_CA1F, propFwdboot_CA1F] = bootci(10000,@propFwdNL, Tok.RCreg_R2(CA1_F_idx),Tok.RCreg_p1(CA1_F_idx));
[propFwdCI_CA3N, propFwdboot_CA3N] = bootci(10000,@propFwdNL, Tok.RCreg_R2(CA3_N_idx),Tok.RCreg_p1(CA3_N_idx));
[propFwdCI_CA3F, propFwdboot_CA3F] = bootci(10000,@propFwdNL, Tok.RCreg_R2(CA3_F_idx),Tok.RCreg_p1(CA3_F_idx));
meanPropFCA1N = mean(propFwdboot_CA1N);
meanPropFCA1F = mean(propFwdboot_CA1F);
meanPropFCA3N = mean(propFwdboot_CA3N);
meanPropFCA3F = mean(propFwdboot_CA3F);

    % pairwise comparisons N vs F
    DeltaPropFCA1_NvsF = propFwdboot_CA1F - propFwdboot_CA1N; 
    DeltaPropFCA1_NvsF_CI = prctile(DeltaPropFCA1_NvsF,[2.5, 97.5]);
    DeltaPropFCA3_NvsF = propFwdboot_CA3F - propFwdboot_CA3N; 
    DeltaPropFCA3_NvsF_CI = prctile(DeltaPropFCA3_NvsF,[2.5, 97.5]);

    % interaction between CA and N_F factors
    DeltaDeltaF = DeltaPropFCA1_NvsF - DeltaPropFCA3_NvsF;
    DeltaDeltaF_CI = prctile(DeltaDeltaF,[2.5, 97.5]);

figure % all PFs pooled. Bar plots of props for all Ca and conditions. Delta F vs N for CA1 and CA3, oriented like in Dabest. 
subplot(1,4,1) % proportion of forward, backward and stable/non-signif PFs (all PFs pooled)
    b = bar([meanPropBCA1N meanPropFCA1N 1-meanPropCA1N; meanPropBCA1F meanPropFCA1F 1-meanPropCA1F], 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylim([0 1]);
    ylabel('proportion of PFs');
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['CA1-N';'CA1-F']);
    box off; %axis square;
subplot(1,4,2) % CA1 N vs F bootstrapped delta for Backward and forward shifting
    histogram(DeltaPropBCA1_NvsF, 'EdgeColor', blue, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaPropBCA1_NvsF), 0.01, mean(DeltaPropBCA1_NvsF)-DeltaPropBCA1_NvsF_CI(1),mean(DeltaPropBCA1_NvsF)-DeltaPropBCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1, 'MarkerSize', 8); hold on
    histogram(DeltaPropFCA1_NvsF, 'EdgeColor', red, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaPropFCA1_NvsF), 0.01, mean(DeltaPropFCA1_NvsF)-DeltaPropFCA1_NvsF_CI(1),mean(DeltaPropFCA1_NvsF)-DeltaPropFCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', red, 'MarkerFaceColor', red, 'LineWidth', 1, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
%     xlim([-meanPropBCA1F 1-meanPropBCA1F]);
    xlim([-0.4 0.4]);
    xlabel('\Delta ratio (F - N)'); ylabel('proba');
    title('bootstrapped distrib of diff ratio')
    set(gca,'view',[90 -90])
    box off; %axis square;
subplot(1,4,3) % proportion of forward, backward and stable/non-signif PFs (all PFs pooled)
    b2 = bar([meanPropBCA3N meanPropFCA3N 1-meanPropCA3N; meanPropBCA3F meanPropFCA3F 1-meanPropCA3F], 'stacked', 'FaceColor','flat');
    b2(1).CData = blue; b2(2).CData = red; b2(3).CData = grey;
    ylim([0 1]);
    ylabel('proportion of PFs');
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['CA3-N';'CA3-F']);
    box off; %axis square;
subplot(1,4,4) % CA1 N vs F bootstrapped delta for Backward and forward shifting
    histogram(DeltaPropBCA3_NvsF, 'EdgeColor', blue, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaPropBCA3_NvsF), 0.01, mean(DeltaPropBCA3_NvsF)-DeltaPropBCA3_NvsF_CI(1),mean(DeltaPropBCA3_NvsF)-DeltaPropBCA3_NvsF_CI(2), 'horizontal', 'o', 'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1, 'MarkerSize', 8); hold on
    histogram(DeltaPropFCA3_NvsF, 'EdgeColor', red, 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaPropFCA3_NvsF), 0.01, mean(DeltaPropFCA3_NvsF)-DeltaPropFCA3_NvsF_CI(1),mean(DeltaPropFCA3_NvsF)-DeltaPropFCA3_NvsF_CI(2), 'horizontal', 'o', 'Color', red, 'MarkerFaceColor', red, 'LineWidth', 1, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-0.4 0.4]);
    xlabel('\Delta ratio (F - N)'); ylabel('proba');
    title('bootstrapped distrib of diff ratio')
    set(gca,'view',[90 -90])
    box off; %axis square;

%% local functions

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


%for n = 1:size(CA_VR_idx)
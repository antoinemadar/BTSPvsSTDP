clear 
close all

Minlaps = 15;
Minlaps2 = 30;

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

L = 300; % track length, in cm
Nbin = 50; % number of spatial bins
BinCenters = (1:Nbin).*6; %L/(2*Nbin):L/Nbin:L; %in cm

COMtrajMat = [];
COMonset = [];
Region = {};
VR = {};
TotNumPFs = 0;
SmallTraj = 0;
Allshifts = [];

T.File = {};
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
T.COM_meanPF = [];
T.SD_meanPF = [];
T.PFsdDiff = [];
T.NumNaN = [];

%loop through the fields
for a = 1:2 %CA1 or CA3
    fn1 = fieldnames(COMdata.(fn{a})); 
    for i = 1:2 % VR: f or n
        fn2 = fieldnames(COMdata.(fn{a}).(fn1{i}));
        DcomAll{a,i} = [];
        for j = 1:length(fn2) % imaging planes
            
            COMloc = (COMdata.(fn{a}).(fn1{i}).(fn2{j}).COM) .*6; % multiply by 6 to get location in cm (this is inherited from Can's code)
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
           
            for k = 1:NumPFs % for all PFs in the plane           
            
                % extract COM trajectory, starting from onset lap, and centered to initial COM location
                COMtraj = COMloc(k,Onset(k):end) - COMloc(k,Onset(k)); 

                % Interpolate NaNs
                x = COMtraj; nanx = isnan(x); laps = 1:numel(x);
                x(nanx) = interp1(laps(~nanx), x(~nanx), laps(nanx)); % no extrapolation for NaNs at the end
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
                
                  modelfun = fittype( @(p1,p2,p3,x) p1*(1-exp(-x/p2))+p3 );
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
                    COMtraj3 = COMtraj2(1:Minlaps);
                    COMtrajMat = [COMtrajMat ; COMtraj3]; % append trajectory to matrix (rows = PFs = PCA observations; columns = laps = variables)
                    COMonset = [COMonset; Onset(k)];
                    Region = [Region; fn{a}]; % CA1 or CA3
                    VR = [VR; fn1{i}]; % F or N
                elseif TrajLength < Minlaps
                    SmallTraj = SmallTraj + 1; %count the number of excluded PFs
                    COMtraj3 = zeros(1,Minlaps)*NaN; % row of NaNs to exclude the trajectories with not enough laps
                    COMtrajMat = [COMtrajMat ; COMtraj3];
                end

                % prepare a table to store different values for each PF 
                T.File{end+1} = fn2{j};
                T.PF(end+1) = k;
                T.Region{end+1} = fn{a};

                % parse info from namefile
                S = strfind(fn2(j), 'plain');
                T.PlaneID{end+1} = fn2{j}(S{1}+5);
                T.animalID{end+1} = COMdata.(fn{a}).animalID{j};
                T.f_n{end+1} = fn1{i};
                VRnum = extract(fn2(j),'00'+digitsPattern); 
                T.VRstr{end+1}= [T.f_n{end} VRnum{end}]; % or use an if loop to replace with scalar code as in the Speed table
                
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

                T.RCreg_p1 = [T.RCreg_p1; mdl.p1];
                T.RCreg_p2 = [T.RCreg_p2; mdl.p2];
                T.RCreg_p3 = [T.RCreg_p3; mdl.p3];
                T.RCreg_R2 = [T.RCreg_R2; gof.rsquare];
                T.RCreg_R2adj = [T.RCreg_R2adj; gof.adjrsquare];
                
                % number of laps without PF after onset
                NumNaN = length(find(isnan(COMloc(k,Onset(k):end))));
                T.NumNaN = [T.NumNaN; NumNaN];

                % compute lapwise and average PF s.d.
                PF = COMdata.(fn{a}).(fn1{i}).(fn2{j}).PFs{k}'; % lapwise normalized activity map of a given PF
                
                Nlaps = size(PF,1);
                for lap = 1:Nlaps
                PFsd(lap) = sqrt( sum( (BinCenters - COMloc(k,lap)).^2.*PF(lap,:)/sum(PF(lap,:)) ) ); %lap-wise SD of the PF
                end
                PFsdOn = PFsd(~isnan(PFsd)); % keep only active laps
                NlapOn = length(PFsdOn);
                PFsdDiff = mean(PFsdOn(NlapOn-2:end)) - mean(PFsdOn(1:3));

                meanPF = mean(PF,1);
                COM_meanPF = sum(meanPF.*BinCenters)/sum(meanPF);
                SD_meanPF = sqrt( sum( (BinCenters - COM_meanPF).^2.*meanPF/sum(meanPF) ) );
                
%                 T.lapPFsd{end+1} = PFsd;
%                 T.meanPF{end+1} = meanPF;
                T.COM_meanPF = [T.COM_meanPF; COM_meanPF];
                T.SD_meanPF = [T.SD_meanPF; SD_meanPF];
                T.PFsdDiff = [T.PFsdDiff; PFsdDiff];
                clear PF Nlaps PFsd meanPF COM_meanPF SD_meanPF
            
            end
        end
    end
end

COMcol1 = COMtrajMat(:,1);
OkTraj = length(COMcol1(~isnan(COMcol1)))
TotTraj = OkTraj + SmallTraj
TotNumPFs

COMtrajMat2 = COMtrajMat(~isnan(COMcol1),:); % matrix without NaNs, i.e. without the PFs with not enough laps
NaNtraj = isnan(COMcol1);
OKtrajIdx = ~isnan(COMcol1);

% make table
T.File = T.File';
T.animalID = T.animalID';
T.Region = T.Region'; 
T.PlaneID = T.PlaneID';
T.f_n = T.f_n';
T.VRstr = T.VRstr';
T.VR = T.VR';
T.PF = T.PF';
T.animalSpeed = T.animalSpeed';
T.maxS_LapSpeed = T.maxS_LapSpeed';
T.maxS_PrevLapSpeed = T.maxS_PrevLapSpeed';

% T = rmfield(T,'RCreg_pval');
Can1All = struct2table(T);
% filename = 'Can1All.xlsx';
% writetable(Can1All,filename);

% PCA, tsne, clustering

[coeff,score,latent,tsquared,explained,mu] = pca(COMtrajMat2, 'centered', false);

for nk = 1:20 
[Kidx(:,nk),K_C{nk},sumd{nk},K_D{nk}] = kmeans(score, nk, 'Display', 'final');
totsumd(nk) = sum(sumd{nk},'all'); % within-cluster sums of point-to-centroid distances in the k-by-1 vector sumd
% totKD(nk) = sum(K_D{nk},'all'); %distances from each point to every centroid in the n-by-k matrix D
end

[TSNEpca, loss_pca] = tsne(score);
[TSNEraw, loss_raw] = tsne(COMtrajMat2);

% Groups per conditions, for figures
CA1idx = find(ismember(Region,'CA1'));
CA3idx = find(ismember(Region,'CA3'));
CA1_F_idx = find(ismember(Region,'CA1') & ismember(VR,'f'));
CA3_F_idx = find(ismember(Region,'CA3') & ismember(VR,'f'));
CA1_N_idx = find(ismember(Region,'CA1') & ismember(VR,'n'));
CA3_N_idx = find(ismember(Region,'CA3') & ismember(VR,'n'));
GroupCA_VR(CA1_N_idx) = 1;
GroupCA_VR(CA1_F_idx) = 2;
GroupCA_VR(CA3_N_idx) = 3;
GroupCA_VR(CA3_F_idx) = 4;

%colorcode CA1 vs CA3, N (dark) vs F (light); e.g. CA1F -> mapCA1(1,:)
map = brewermap(8,'Paired'); %cool(4); %brewermap(4,'PiYG')
mapCA1 = map(3:4,:); %brewermap(2,'GnBu');
mapCA3 = map(7:8,:); %brewermap(2,'OrRd');
mapCA = [mapCA1;mapCA3];
colors = [mapCA1(2,:); mapCA1(1,:); mapCA3(2,:); mapCA3(1,:)];

% fwd PFs should be red, Bck PFs blue, NonSig grey or black
blue = [0 0.45 0.74]; % backward
red = [0.85 0.32 0.01]; % forward
grey = [0.5 0.5 0.5]; %nonsig

T2 = Can1All;
Tok = T2(OKtrajIdx,:);

filename = 'Tok.xlsx';
writetable(Tok,filename);
Tpca = array2table([COMtrajMat2, score(:,1), score(:,2), score(:,3)], 'VariableNames',{'COM lap1','COM lap2', 'COM lap3', 'COM lap4', 'COM lap5', 'COM lap6', 'COM lap7', 'COM lap8', 'COM lap9', 'COMlap10', 'COMlap11', 'COMlap12', 'COMlap13', 'COMlap14', 'COMlap15', 'PC1', 'PC2', 'PC3'});
writetable(Tpca,'Tpca.xlsx');

%%
Fwd_idx = find(Tok.COM_slope > 0 & Tok.COM_pval<0.05);
Bwd_idx = find(Tok.COM_slope < 0 & Tok.COM_pval<0.05);
Nsig_idx = find(Tok.COM_pval>=0.05);

for k = 1:OkTraj % for each PFs
    if Tok.COM_slope(k) > 0 & Tok.COM_pval(k) < 0.05
    Tok.shiftDir{k,:} = 'F';
    elseif Tok.COM_slope(k) < 0 & Tok.COM_pval(k) < 0.05
    Tok.shiftDir{k,:} = 'B';    
    else
    Tok.shiftDir{k,:} = 'N';    
    end
end

figure % lin Shift vs Rsquare
gscatter(Tok.COM_slope,Tok.COM_R2, Tok.shiftDir, [blue; grey; red])
xlabel('lin reg slope (cm/lap)'); ylabel('lin reg R2');
title('All conditions')
box off; axis square;

figure % PC1 vs lin R2
gscatter(score(:,1),Tok.COM_R2, Tok.shiftDir, [blue; grey; red])
xlabel('PC1 score'); ylabel('lin reg R2');
box off; axis square;

% figure
% scatter3(Tok.COM_slope(Fwd_idx),Tok.COM_R2(Fwd_idx), score(Fwd_idx,1), 2, 'r'); hold on
% scatter3(Tok.COM_slope(Bwd_idx),Tok.COM_R2(Bwd_idx), score(Bwd_idx,1), 2, 'b'); hold on
% scatter3(Tok.COM_slope(Nsig_idx),Tok.COM_R2(Nsig_idx), score(Nsig_idx,1), 2, 'k'); hold on
% xlabel('lin reg slope (cm/lap)'); ylabel('lin reg R2'); zlabel('PC1 score');
% box off; axis square;

figure % lin R2 vs RC-like R2 and R2 adjusted
% subplot(1,2,1)
% gscatter(Tok.RCreg_R2adj,Tok.COM_R2, Tok.shiftDir, [blue; grey; red], 'o', 5); hold on
% plot([0 1], [0 1], 'k--', 'LineWidth', 2);
% xlabel('RC-like R2'); ylabel('lin reg R2');
% xlim([0 1])
% box off; axis square;
% subplot(1,2,2)
gscatter(Tok.RCreg_R2adj,Tok.COM_R2, Tok.shiftDir, [blue; grey; red], 'o', 5); hold on
plot([0 1], [0 1], 'k--', 'LineWidth', 2);
xlabel('RC-like adjusted R2'); ylabel('lin reg R2');
xlim([0 1])
box off; axis square;

% figure % with R-square, not adjusted
% subplot(1,3,1)
% gscatter(Tok.RCreg_p1,Tok.RCreg_p2, Tok.shiftDir, [blue; grey; red])
% xlabel('RCreg amplitude, in cm'); ylabel('RCreg time constant, in laps');
% % box off; axis square;
% subplot(1,3,2)
% gscatter(abs(Tok.RCreg_p1),Tok.RCreg_R2, Tok.shiftDir, [blue; grey; red])
% ylim([0 1])
% xlabel('RCreg |Amp|, in cm'); ylabel('RCreg R2');
% % box off; axis square;
% subplot(1,3,3)
% gscatter(Tok.RCreg_p2,Tok.RCreg_R2, Tok.shiftDir, [blue; grey; red])
% ylim([0 1])
% xlabel('RCreg time constant, in laps'); ylabel('RCreg R2');

figure % same as above with R2adj
subplot(1,2,1)
gscatter(abs(Tok.RCreg_p1),Tok.RCreg_R2adj, Tok.shiftDir, [blue; grey; red], 'o', 5)
ylim([0 1])
xlabel('RCreg |Amp|, in cm'); ylabel('RCreg R2');
% box off; axis square;
subplot(1,2,2)
gscatter(Tok.RCreg_p2,Tok.RCreg_R2adj, Tok.shiftDir, [blue; grey; red], 'o', 5)
ylim([0 1])
xlabel('RCreg time constant, in laps'); ylabel('RCreg R2');

%% width and width increase and relationship with linear fit
% CA1idx = find(ismember(T2.Region,'CA1'));
% CA3idx = find(ismember(T2.Region,'CA3'));

CA1idx = find(ismember(Tok.Region,'CA1'));
CA3idx = find(ismember(Tok.Region,'CA3'));
CA1sig_idx = find(ismember(Tok.Region,'CA1') & Tok.COM_pval<0.05);
CA3sig_idx = find(ismember(Tok.Region,'CA3') & Tok.COM_pval<0.05);

% meanPF width
figure
histogram(Tok.SD_meanPF(CA1idx), 'Normalization', 'pdf', 'LineWidth', 2, 'DisplayStyle','stairs', 'EdgeColor', mapCA1(2,:)); hold on
xline(median(Tok.SD_meanPF(CA1idx)),'-', 'Color', mapCA1(2,:), 'LineWidth', 1);
histogram(Tok.SD_meanPF(CA3idx), 'Normalization', 'pdf', 'LineWidth', 2, 'DisplayStyle','stairs', 'EdgeColor', mapCA3(2,:)); hold on
xline(median(Tok.SD_meanPF(CA3idx)),'-', 'Color', mapCA3(2,:), 'LineWidth', 1);
xlabel('PF sd (cm)'); ylabel('proba. density');
set(gca, 'XTick', 0:10:80)
title({['CA1F and N combined, median = ' num2str(median(Tok.SD_meanPF(CA1idx))) ' cm']; ['CA3 F and N combined, median = ' num2str(median(Tok.SD_meanPF(CA3idx))) ' cm' ]})
box off; axis square;
mean(Tok.SD_meanPF(CA3idx));

figure
subplot(1,2,1)
histogram(Tok.SD_meanPF(CA1_N_idx), 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA1(2,:)); hold on
xline(median(Tok.SD_meanPF(CA1_N_idx)),'g-');hold on
xlabel('PF s.d. (cm)'); ylabel('probability');
title(['CA1 N, median = ' num2str(median(Tok.SD_meanPF(CA1_N_idx))) ' cm' ])
box off; axis square;
subplot(1,2,2)
histogram(Tok.SD_meanPF(CA3_N_idx), 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA3(2,:)); hold on
xline(median(Tok.SD_meanPF(CA3_N_idx)),'r-');hold on
xlabel('PF s.d. (cm)'); ylabel('probability');
title(['CA3 N, median = ' num2str(median(Tok.SD_meanPF(CA3_N_idx))) ' cm' ])
box off; axis square;

% figure
% subplot(1,2,1)
% histogram(Tok.SD_meanPF(CA1idxOk), 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA1(2,:)); hold on
% xline(median(T2.SD_meanPF(CA1idxOk)),'g-');hold on
% xlabel('PF s.d. (cm)'); ylabel('probability');
% title(['CA1 F and N combined, median = ' num2str(median(T2.SD_meanPF(CA1idxOk))) ' cm' ])
% box off; axis square;
% subplot(1,2,2)
% histogram(Tok.SD_meanPF(CA3idxOk), 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA3(2,:)); hold on
% xline(median(T2.SD_meanPF(CA3idxOk)),'r-');hold on
% xlabel('PF s.d. (cm)'); ylabel('probability');
% title(['CA3 F and N combined, median = ' num2str(median(T2.SD_meanPF(CA3idxOk))) ' cm' ])
% box off; axis square;

figure % PFsd vs abs lin shift
subplot(2,1,1)
gscatter(Tok.SD_meanPF(CA1idx), abs(Tok.COM_slope(CA1idx)), Tok.shiftDir(CA1idx), [blue; grey; red], '.')
ylabel('COM slope (cm/lap)'); xlabel('mean PF s.d. (cm)');
% ylim([-4 3])
title('CA1')
axis square; box off;
subplot(2,1,2)
gscatter(Tok.SD_meanPF(CA3idx), abs(Tok.COM_slope(CA3idx)), Tok.shiftDir(CA3idx), [grey; red; blue], '.')
ylabel('COM slope (cm/lap)'); xlabel('mean PF s.d. (cm)');
% ylim([-4 3])
title('CA3')
box off; axis square;

figure % PFsd vs abs lin shift for sig shifts
subplot(2,1,1)
gscatter(Tok.SD_meanPF(CA1sig_idx), abs(Tok.COM_slope(CA1sig_idx)), Tok.shiftDir(CA1sig_idx), [blue; red], '.')
ylabel('COM slope (cm/lap)'); xlabel('mean PF s.d. (cm)');
% ylim([-4 3])
title('CA1')
box off; axis square
subplot(2,1,2)
gscatter(Tok.SD_meanPF(CA3sig_idx), abs(Tok.COM_slope(CA3sig_idx)), Tok.shiftDir(CA3sig_idx), [red; blue], '.')
ylabel('COM slope (cm/lap)'); xlabel('mean PF s.d. (cm)');
% ylim([-4 3])
title('CA3')
box off; axis square
% there is a clear relationship, which is expected because shifting will
% necessarily induce a larger average, artefactually. The PF may actually
% stay the same size, but the average will increase linearly the bigger the
% shift. 

% change in PF width
figure
histogram(Tok.PFsdDiff(CA1idx), 'Normalization', 'pdf', 'LineWidth', 2, 'DisplayStyle','stairs', 'EdgeColor', mapCA1(2,:)); hold on
xline(median(Tok.PFsdDiff(CA1idx)),'-', 'Color', mapCA1(2,:), 'LineWidth', 1);
histogram(Tok.PFsdDiff(CA3idx), 'Normalization', 'pdf', 'LineWidth', 2, 'DisplayStyle','stairs', 'EdgeColor', mapCA3(2,:)); hold on
xline(median(Tok.PFsdDiff(CA3idx)),'-', 'Color', mapCA3(2,:), 'LineWidth', 1);
xlabel('\Delta PF s.d. (cm)'); ylabel('proba. density');
xlim([-35 35]);
% set(gca, 'XTick', 0:10:80)
title({['CA1F and N combined, median = ' num2str(median(Tok.PFsdDiff(CA1idx))) ' cm']; ['CA3 F and N combined, median = ' num2str(median(Tok.PFsdDiff(CA3idx))) ' cm' ]})
box off; axis square;
mean(Tok.PFsdDiff(CA3idx));

figure
histogram(Tok.PFsdDiff(CA1sig_idx), 'Normalization', 'pdf', 'LineWidth', 2, 'DisplayStyle','stairs', 'EdgeColor', mapCA1(2,:)); hold on
xline(median(Tok.PFsdDiff(CA1sig_idx)),'-', 'Color', mapCA1(2,:), 'LineWidth', 1);
histogram(Tok.PFsdDiff(CA3sig_idx), 'Normalization', 'pdf', 'LineWidth', 2, 'DisplayStyle','stairs', 'EdgeColor', mapCA3(2,:)); hold on
xline(median(Tok.PFsdDiff(CA3sig_idx)),'-', 'Color', mapCA3(2,:), 'LineWidth', 1);
xlabel('\Delta PF s.d. (cm)'); ylabel('proba. density');
xlim([-35 35]);
% set(gca, 'XTick', 0:10:80)
title({['CA1F and N combined, median = ' num2str(median(Tok.PFsdDiff(CA1sig_idx))) ' cm']; ['CA3 F and N combined, median = ' num2str(median(Tok.PFsdDiff(CA3sig_idx))) ' cm' ]})
box off; axis square;

figure
subplot(1,2,1)
histogram(Tok.PFsdDiff(CA1idx), 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA1(2,:)); hold on
xline(mean(Tok.PFsdDiff(CA1idx)),'g-');hold on
xlabel('\Delta PF s.d. (cm)'); ylabel('probability');
title(['CA1 F and N combined, median = ' num2str(median(Tok.PFsdDiff(CA1idx))) ' cm' ])
box off; axis square;
subplot(1,2,2)
histogram(Tok.PFsdDiff(CA3idx), 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA3(2,:)); hold on
xline(mean(Tok.PFsdDiff(CA3idx)),'r-');hold on
xlabel('\Delta PF s.d. (cm)'); ylabel('probability');
title(['CA3 F and N combined, mean = ' num2str(mean(Tok.PFsdDiff(CA3idx))) ' cm' ])
box off; axis square;

figure % PFsdDiff vs lin shift
subplot(1,2,1)
gscatter(Tok.PFsdDiff(CA1idx), Tok.COM_slope(CA1idx), Tok.shiftDir(CA1idx), [blue; grey; red], 'o')
ylim([-4 3])
ylabel('lin reg slope (cm/lap)'); xlabel('\Delta PF s.d. (cm)');
title('CA1')
box off; axis square
subplot(1,2,2)
gscatter(Tok.PFsdDiff(CA3idx), Tok.COM_slope(CA3idx), Tok.shiftDir(CA3idx), [grey; red; blue], 'o')
ylim([-4 3])
ylabel('lin reg slope (cm/lap)'); xlabel('\Delta PF s.d. (cm)');
title('CA3')
box off; axis square

figure % PFsdDiff vs abs lin shift
subplot(1,2,1)
gscatter(Tok.PFsdDiff(CA1idx), abs(Tok.COM_slope(CA1idx)), Tok.shiftDir(CA1idx), [blue; grey; red], 'o')
ylabel('abs slope (cm/lap)'); xlabel('\Delta PF s.d. (cm)');
title('CA1')
box off; axis square
subplot(1,2,2)
gscatter(Tok.PFsdDiff(CA3idx), abs(Tok.COM_slope(CA3idx)), Tok.shiftDir(CA3idx), [grey; red; blue], 'o')
ylabel('abs slope (cm/lap)'); xlabel('\Delta PF s.d. (cm)');
title('CA3')
box off; axis square

[Bca1,~,~,~,statsBca1] = regress(abs(Tok.COM_slope(CA1sig_idx)), [ones(size(CA1sig_idx)), Tok.PFsdDiff(CA1sig_idx)]);
[Bca3,~,~,~,statsBca3] = regress(abs(Tok.COM_slope(CA3sig_idx)), [ones(size(CA3sig_idx)), Tok.PFsdDiff(CA3sig_idx)]);

figure % PFsdDiff vs lin shift for significant shifts only, with regression
% T.COMinterp_slope = [T.COMinterp_slope; b(2)];
%                 T.COMinterp_intcp = [T.COMinterp_intcp; b(1)];
%                 T.COMinterp_R2 = [T.COMinterp_R2; stats(1)];
%                 T.COMinterp_pval = [T.COMinterp_pval; stats(3)];

subplot(1,2,1)
gscatter(Tok.PFsdDiff(CA1sig_idx), abs(Tok.COM_slope(CA1sig_idx)), Tok.shiftDir(CA1sig_idx), [blue; red], 'o')
ylabel('lin reg slope (cm/lap)'); xlabel('\Delta PF s.d. (cm)');
title('CA1')
box off; axis square
subplot(1,2,2)
gscatter(Tok.PFsdDiff(CA3sig_idx), abs(Tok.COM_slope(CA3sig_idx)), Tok.shiftDir(CA3sig_idx), [red; blue], 'o')
ylabel('abs slope (cm/lap)'); xlabel('\Delta PF s.d. (cm)');
title('CA3')
box off; axis square

%% RCreg param distributions

figure
subplot(1,2,1)
histogram(Tok.RCreg_p1, 'Normalization', 'probability')
xlabel('RCreg amplitude, in cm'); ylabel('PF frequency')
box off; axis square;
subplot(1,2,2)
histogram(Tok.RCreg_p2, 40, 'Normalization', 'probability')
xlabel('RCreg time constant, in laps'); ylabel('PF frequency')
box off; axis square;

figure
subplot(1,3,1)
gscatter(Tok.RCreg_p1,Tok.RCreg_p2, Tok.shiftDir, [blue; grey; red])
xlabel('RCreg amplitude, in cm'); ylabel('RCreg time constant, in laps');
% box off; axis square;
subplot(1,3,2)
gscatter(abs(Tok.RCreg_p1),Tok.RCreg_R2, Tok.shiftDir, [blue; grey; red])
ylim([0 1])
xlabel('RCreg |Amp|, in cm'); ylabel('RCreg R2');
% box off; axis square;
subplot(1,3,3)
gscatter(Tok.RCreg_p2,Tok.RCreg_R2, Tok.shiftDir, [blue; grey; red])
ylim([0 1])
xlabel('RCreg time constant, in laps'); ylabel('RCreg R2');
% box off; axis square;

%% examples

% ex stable PF / long Tok.RCreg_p2 (tau, in laps) 
stablePFs = find(Tok.RCreg_p2 == max(Tok.RCreg_p2));
stablePFs1 = find(Tok.RCreg_p2 == max(Tok.RCreg_p2) & Tok.COM_R2<0.05 & Tok.COM_pval>0.05);
stablePFs2 = find(Tok.RCreg_p2 >99.99 & Tok.COM_pval>0.05);
stablePFs3 = find(Tok.RCreg_p2 >50 & Tok.COM_pval>0.05 & Tok.NumNaN ==0);
% stablePFs2 = find(Tok.RCreg_p2 > 50 & Tok.COM_R2<0.05);
stablePFsex = [2286, 1326, 2347, 30, 61, 75]; 

% Ex linear backward shifting 
exPF(1) = find(Tok.COM_slope < -3.6); % largest backward shift
[maxLinR2, exPF(2)] = max(Tok.COM_R2); % best linear fit
exPF(3) = find(Tok.COM_slope < -3.4 & Tok.COM_R2>0.8); % optimal linear fit (largest slope/Rsquare combination)
exPF(14) = find(Tok.COM_R2>0.60 & Tok.RCreg_R2adj < 0.17);
goodlinB = find(Tok.COM_R2>Tok.RCreg_R2adj & Tok.COM_R2>0.8 & Tok.COM_slope < 0);
goodlinB2 = find(Tok.COM_R2>Tok.RCreg_R2adj & Tok.COM_R2>0.5 & Tok.COM_slope < 0 & Tok.NumNaN == 0);
goodlinB3 = find(Tok.COM_R2>Tok.RCreg_R2adj & Tok.COM_R2>0.6 & Tok.COM_slope < 0 & Tok.NumNaN == 4);

% Ex linear forward shifting 
exPF(4) = find(Tok.COM_slope > 2); %largest forward shift
exPF(5) = find(Tok.COM_slope > 1.6 & Tok.COM_R2>0.6); % optimum backward shift --> good late BTSP potential ex
exPF(6) = find(Tok.COM_R2>0.743 & Tok.COM_slope > 0); % best linear fit --? another good late BTSP example
% exPF(15) = find(Tok.COM_R2>0.74 & Tok.RCreg_R2adj < 0.40); % it's actually the same as exPF(6) i.e. ForwardBest
goodlinF = find(Tok.COM_R2>Tok.RCreg_R2adj & Tok.COM_R2>0.6 & Tok.COM_slope > 0);
goodlinF2 = find(Tok.COM_R2>Tok.RCreg_R2adj & Tok.COM_R2>0.4 & Tok.COM_slope > 0 & Tok.NumNaN == 0);
goodlinF3 = find(Tok.COM_R2>0.5 & Tok.COM_slope > 0 & Tok.NumNaN <3);
goodlinF5 = find(Tok.COM_R2>0.4 & Tok.COM_slope > 0)

% ex high absolute PC1 score
[minPC1, exPF(7)] = min(score(:,1)); % min PC1 score
[maxPC1, exPF(8)] = max(score(:,1)); % max PC1 score
exPF(9) = find(Tok.COM_R2>0.7 & score(:,1)<-167 ); % good linear fit and PC1 score

% ex high absolute scores for other PCs
[~, exPC2] = maxk(abs(score(:,2)), 10);
[~, exPC3] = maxk(abs(score(:,3)), 10);
[~, exPC4] = maxk(abs(score(:,4)), 10);
Weird2ndaryPCs = [2428, 2425, 2300, 2180, 2158, 1289, 1242, 1019, 318];

% Ex good RC-like
[maxNLR2, exPF(10)] = max(Tok.RCreg_R2); % 
exPF(11) = find(Tok.RCreg_p2 < 2.3 & Tok.RCreg_R2>0.9); % good example of nice PC1 looking trajectory
betterRCs = find(Tok.RCreg_R2adj > 0.75 & Tok.COM_R2<0.5);

% Ex good RC-like but not signif linear shift
exPF(12) = find(Tok.RCreg_R2 > 0.75 & Tok.COM_R2<0.03);
betterRCs2 = find(Tok.RCreg_R2adj > 0.6 & Tok.COM_pval>0.05);

% ex good RC-like in CA3
goodRC_CA3 = find(Tok.RCreg_R2 > 0.6 & ismember(Tok.Region,'CA3'));

% nonsig linear with poor RC-like
badLinRCmeh = find(Tok.RCreg_R2adj > 0.3 & Tok.RCreg_R2adj <0.35 & Tok.COM_pval>0.05);
badLinRCmeh2 = find(Tok.RCreg_R2adj > 0.1 & Tok.RCreg_R2adj <0.11 & Tok.COM_pval>0.05); % mostly flat, or sputtery with varying width but looking flat
[minR2adj, exPF(13)] = min(Tok.RCreg_R2adj); % also corresponds to highest PC1 

% ex increase in PFsd 
BackSDup = find(Tok.COM_R2>Tok.RCreg_R2adj & Tok.COM_R2>0.5 & Tok.COM_slope < 0 & Tok.PFsdDiff > 10);
FwdSDup = find(Tok.COM_R2>Tok.RCreg_R2adj & Tok.COM_R2>0.5 & Tok.COM_slope > 0 & Tok.PFsdDiff > 4);

% ex decrease in PFsd
BackSDdown = find(Tok.COM_R2>Tok.RCreg_R2adj & Tok.COM_R2>0.5 & Tok.COM_slope < 0 & Tok.PFsdDiff < -8);
FwdSDdown = find(Tok.COM_R2>Tok.RCreg_R2adj & Tok.COM_R2>0.5 & Tok.COM_slope > 0 & Tok.PFsdDiff < -4);

% Ex zigzag (late BTSP) (show where on linR2 vs nlR2 plot it is, or other
% plot used to find the PF)
zigzag(1) = 943;
zigzag(2:4) = exPF(1:3);
zigzag(5) = 512;
zigzag(6) = exPF(5); % optimum backward shift --> good late BTSP potential ex
zigzag(7) = exPF(6); % best linear fit --? another good late BTSP example
zigzag(8) = 1411; % almost midsplit, BTSP event leads to sudden shift
zigzag(9) = 2504; % CA3 best RC
zigzag(10) = 2540; % CA3 good RC and good linear

% Ex linear & Stable PFs
BestLinBck = [511, 1149, 682];
BestLinFwd = [2261, 2290, 549, 2280];
BestStable = [1326, 2286, 2347, 30, 61, 75];
BestLin = [BestLinBck, BestLinFwd, BestStable];



%% Plot ratemaps and COM regressions of chosen example
ex = 986;

Bmap = brewermap(256,'*Spectral');

figure % lapwise rate map
    COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
    RateMap = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).PFs{Tok.PF(ex)};
    % imagesc(rot90(RateMap')); hold on
    imagesc(flipud(RateMap)); hold on
    scatter(1:length(COMbin), 51-COMbin, 10, [1 0.5 0], 'filled');
    set(gca, 'YTick', [1 50], 'YTickLabel', ['  end'; 'start']);
    colormap(jet(256)); %colormap(Bmap); %colormap(flipud(hot(256)));% colormap(flipud(gray(256)));
    c = colorbar; c.Label.String = 'DeltaF/F';
    xlabel('running laps'); ylabel('position on track')
    axis square; box off;
    title({['PF' num2str(ex) ' ' Tok.Region{ex}] ;[' lin R2 = ' num2str(Tok.COM_R2(ex)) ', nl R2 = ' num2str(Tok.RCreg_R2(ex)) ' , PC1 = ' num2str(score(ex,1))]})
%     title(['PF' num2str(ex) ' ' Tok.Region{ex} ' lin R2 = ' num2str(Tok.COM_R2(ex)) ', nl adjR2 = ' num2str(Tok.RCreg_R2adj(ex)) ' , PC1 = ' num2str(score(ex,1))])

figure
    COMa = Tok.COM_slope(ex);
    if COMa > 0
        colorp = red;
    elseif COMa < 0
        colorp = blue;
    else
        colorp = grey;
    end
    COMb= Tok.COM_intcp(ex);
    nlAmp = Tok.RCreg_p1(ex);
    nlTau = Tok.RCreg_p2(ex);
    nlInt = Tok.RCreg_p3(ex);
    onset = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).onset_lap(Tok.PF(ex));
    COMtraj = (COMbin(onset:end) - COMbin(onset)).*6;
    COMlinreg = COMb + COMa.*[1:length(COMtraj)];
    COMnlreg = nlAmp*(1-exp(-[0:0.1:length(COMtraj)-1]./nlTau)) + nlInt;
    scatter(1:length(COMtraj),COMtraj, 50, colorp, 'filled', 'MarkerFaceAlpha', 0.5); hold on 
    plot(1:length(COMtraj),COMlinreg, '-', 'Color', colorp, 'LineWidth', 2); hold on
    plot(1:0.1:length(COMtraj),COMnlreg, '-', 'Color', [0 1 0], 'LineWidth', 2)
    xlabel('laps after onset'); 
    ylabel('onset-centered COM position (cm)');
    ylim([-200 200]); xlim([1 length(COMtraj)])
    axis square; box off;
    title(['linreg Slope = ' num2str(Tok.COM_slope(ex)) ', nlreg Amp (cm) = ' num2str(Tok.RCreg_p1(ex)) ' , nlreg tau (laps) = ' num2str(Tok.RCreg_p2(ex))])

%% plot all PFS of a given group (e.g. betterRCs or betterRCs2)
close all
groupPFs = Weird2ndaryPCs; %[511, 1326, 2261, 2290, 549, 2280] %stablePFsex %2(1:20); %[167, 618, 741, 1503, 1674, 2324, 2327, 1854, 1471, 808, 405, 37];
% f = figure
for n = 1:length(groupPFs)
    figure(n)
    ex = groupPFs(n); 
    subplot(1,2,1)
%     subplot(length(groupPFs)/2,2*2, 1 + 2*(n-1), 'align')
        COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
        RateMap = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).PFs{Tok.PF(ex)};
        imagesc(flipud(RateMap)); hold on
%         scatter(1:length(COMbin), 51-COMbin, 10, [1 1 1], 'filled');
        set(gca, 'YTick', [1 50], 'YTickLabel', ['  end'; 'start']);
        colormap(jet(256)); %colormap colormap(flipud(hot(256)));% colormap(flipud(gray(256)));
        xlabel('running laps'); ylabel('position on track')
        axis square; box off;
%         title(['PF' num2str(ex) ' ' Tok.Region{ex}]);
        title({['PF' num2str(ex) ' ' Tok.Region{ex}] ;['lin R2 = ' num2str(Tok.COM_R2(ex))]; ['nl R2 = ' num2str(Tok.RCreg_R2(ex))]; ['PC1 = ' num2str(score(ex,1))]})
      subplot(1,2,2)
%     subplot(length(groupPFs)/2,2*2,2 + 2*(n-1), 'align')
        COMa = Tok.COM_slope(ex);
        Pa = Tok.COM_pval(ex);
        if Pa < 0.05 && COMa > 0
            colorp = red;
        elseif Pa < 0.05 && COMa < 0
            colorp = blue;
        else
            colorp = grey;
        end
        COMb= Tok.COM_intcp(ex);
        nlAmp = Tok.RCreg_p1(ex);
        nlTau = Tok.RCreg_p2(ex);
        nlInt = Tok.RCreg_p3(ex);
        onset = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).onset_lap(Tok.PF(ex));
        COMtraj = (COMbin(onset:end) - COMbin(onset)).*6;
        COMlinreg = COMb + COMa.*[1:length(COMtraj)];
        COMnlreg = nlAmp*(1-exp(-[0:0.1:length(COMtraj)-1]./nlTau)) + nlInt;
        scatter(1:length(COMtraj),COMtraj, 50, colorp, 'filled', 'MarkerFaceAlpha', 0.5); hold on 
        plot(1:length(COMtraj),COMlinreg, '-', 'Color', colorp, 'LineWidth', 1); hold on
        plot(1:0.1:length(COMtraj),COMnlreg, '-', 'Color', [0 1 0], 'LineWidth', 1)
        xlabel('laps after onset'); 
        ylabel('onset-centered COM position (cm)');
        ylim([-150 150]); xlim([1 length(COMtraj)])
        axis square; box off;
        title({['linreg Slope = ' num2str(Tok.COM_slope(ex))]; ['nlreg Amp (cm) = ' num2str(Tok.RCreg_p1(ex))]; ['nlreg tau (laps) = ' num2str(Tok.RCreg_p2(ex))]})

end
% f.Renderer = 'painters';
% % saveas(f, 'zigzagPFs', 'epsc')
% print(f, '-vector','-dpdf','zigzagPFs.pdf')

%% figure nlreg
clear groupPFs
close all

FigNLreg = [1854, 1674, 1471, 618, 2504, 2327];

groupPFs = FigNLreg;

fNL = figure
for n = 1:length(groupPFs)
    ex = groupPFs(n);
    subplot(2,length(groupPFs), n, 'align')
        COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
        RateMap = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).PFs{Tok.PF(ex)};
        imagesc(RateMap'); hold on
        set(gca, 'XTick', [1 50], 'XTickLabel', ['start'; '  end'], 'YTick', [1 size(RateMap,2)], 'YTickLabel', [' 1'; num2str(size(RateMap,2))]);
        colormap(jet(256)); %colormap colormap(flipud(hot(256)));% colormap(flipud(gray(256)));
%         ylabel('running laps')
%         xlabel('position on track')
        axis square; box off;
        title({['PF' num2str(ex) ' ' Tok.Region{ex}] ;['lin R2 = ' num2str(Tok.COM_R2(ex))]; ['nl R2 = ' num2str(Tok.RCreg_R2(ex))]; ['PC1 = ' num2str(score(ex,1))]})
    subplot(2,length(groupPFs),length(groupPFs)+n, 'align')
            COMa = Tok.COM_slope(ex);
        if COMa > 0
            colorp = red;
        elseif COMa < 0
            colorp = blue;
        else
            colorp = grey;
        end
        COMb= Tok.COM_intcp(ex);
        nlAmp = Tok.RCreg_p1(ex);
        nlTau = Tok.RCreg_p2(ex);
        nlInt = Tok.RCreg_p3(ex);
        onset = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).onset_lap(Tok.PF(ex));
        COMtraj = (COMbin(onset:end) - COMbin(onset)).*6;
        COMlinreg = COMb + COMa.*[1:length(COMtraj)];
        COMnlreg = nlAmp*(1-exp(-[0:0.1:length(COMtraj)-1]./nlTau)) + nlInt;
        scatter(1:length(COMtraj),COMtraj, 30, colorp, 'filled', 'MarkerFaceAlpha', 0.5); hold on 
        plot(1:length(COMtraj),COMlinreg, '-', 'Color', colorp, 'LineWidth', 1); hold on
        plot(1:0.1:length(COMtraj),COMnlreg, '-', 'Color', [0 1 0], 'LineWidth', 1)
        xlabel('laps after onset'); 
%         ylabel('onset-centered COM position (cm)');
        ylim([-80 80]); xlim([1 length(COMtraj)])
        axis square; box off;
        title({['linreg Slope = ' num2str(Tok.COM_slope(ex))]; ['nlreg Amp (cm) = ' num2str(Tok.RCreg_p1(ex))]; ['nlreg tau (laps) = ' num2str(Tok.RCreg_p2(ex))]})

end

% fNL.Renderer = 'painters';
% print(fNL, '-vector','-dpdf','exPFsNL.pdf')

% fNL2 = figure;
% p = panel(fNL2);
% p.pack(2,length(groupPFs));
% % p.margin = 2;
% for n = 1:length(groupPFs)
%     ex = groupPFs(n);
% %     subplot(2,length(groupPFs), n, 'align')
%     p(1, n).select();
%         COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
%         RateMap = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).PFs{Tok.PF(ex)};
%         imagesc(RateMap'); hold on
%         set(gca, 'XTick', [1 50], 'XTickLabel', ['  end'; 'start'], 'YTick', [1 size(RateMap,2)], 'YTickLabel', [' 1'; num2str(size(RateMap,2))]);
%         colormap(jet(256)); %colormap colormap(flipud(hot(256)));% colormap(flipud(gray(256)));
% %         ylabel('running laps')
%         xlabel('position on track')
% %         axis square; box off;
% %         title({['PF' num2str(ex) ' ' Tok.Region{ex}] ;['lin R2 = ' num2str(Tok.COM_R2(ex))]; ['nl R2 = ' num2str(Tok.RCreg_R2(ex))]; ['PC1 = ' num2str(score(ex,1))]})
% %     subplot(2,length(groupPFs),length(groupPFs)+n, 'align')
%     p(2, n).select();
%             COMa = Tok.COM_slope(ex);
%         if COMa > 0
%             colorp = red;
%         elseif COMa < 0
%             colorp = blue;
%         else
%             colorp = grey;
%         end
%         COMb= Tok.COM_intcp(ex);
%         nlAmp = Tok.RCreg_p1(ex);
%         nlTau = Tok.RCreg_p2(ex);
%         nlInt = Tok.RCreg_p3(ex);
%         onset = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).onset_lap(Tok.PF(ex));
%         COMtraj = (COMbin(onset:end) - COMbin(onset)).*6;
%         COMlinreg = COMb + COMa.*[1:length(COMtraj)];
%         COMnlreg = nlAmp*(1-exp(-[0:0.1:length(COMtraj)-1]./nlTau)) + nlInt;
%         scatter(1:length(COMtraj),COMtraj, 30, colorp, 'filled', 'MarkerFaceAlpha', 0.5); hold on 
%         plot(1:length(COMtraj),COMlinreg, '-', 'Color', colorp, 'LineWidth', 1); hold on
%         plot(1:0.1:length(COMtraj),COMnlreg, '-', 'Color', [0 1 0], 'LineWidth', 1)
%         xlabel('laps after onset'); 
% %         ylabel('onset-centered COM position (cm)');
%         ylim([-80 80]); xlim([1 length(COMtraj)])
% %         axis square; box off;
% %         title({['linreg Slope = ' num2str(Tok.COM_slope(ex))]; ['nlreg Amp (cm) = ' num2str(Tok.RCreg_p1(ex))]; ['nlreg tau (laps) = ' num2str(Tok.RCreg_p2(ex))]})
% 
% end


% figure
% ha = tight_subplot(2,length(groupPFs),.03,.1,.1);
% for n = 1:length(groupPFs)
%     ex = groupPFs(n);
%         
%         COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
%         RateMap = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).PFs{Tok.PF(ex)};
%         imagesc(ha(n),RateMap'); hold on
%         set(ha(n), 'XTick', [1 50], 'XTickLabel', ['  end'; 'start'], 'YTick', [1 size(RateMap,2)], 'YTickLabel', [' 1'; num2str(size(RateMap,2))]);
%         colormap(jet(256)); %colormap colormap(flipud(hot(256)));% colormap(flipud(gray(256)));
% %         title({['PF' num2str(ex) ' ' Tok.Region{ex}] ;['lin R2 = ' num2str(Tok.COM_R2(ex))]; ['nl R2 = ' num2str(Tok.RCreg_R2(ex))]; ['PC1 = ' num2str(score(ex,1))]})
% 
%         COMa = Tok.COM_slope(ex);
%         if COMa > 0
%             colorp = red;
%         elseif COMa < 0
%             colorp = blue;
%         else
%             colorp = grey;
%         end
%         COMb= Tok.COM_intcp(ex);
%         nlAmp = Tok.RCreg_p1(ex);
%         nlTau = Tok.RCreg_p2(ex);
%         nlInt = Tok.RCreg_p3(ex);
%         onset = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).onset_lap(Tok.PF(ex));
%         COMtraj = (COMbin(onset:end) - COMbin(onset)).*6;
%         COMlinreg = COMb + COMa.*[1:length(COMtraj)];
%         COMnlreg = nlAmp*(1-exp(-[0:0.1:length(COMtraj)-1]./nlTau)) + nlInt;
%         
%         scatter(ha(5+n),1:length(COMtraj),COMtraj, 30, colorp, 'filled', 'MarkerFaceAlpha', 0.5); hold on 
%         plot(ha(5+n),1:length(COMtraj),COMlinreg, '-', 'Color', colorp, 'LineWidth', 1); hold on
%         plot(ha(5+n),1:0.1:length(COMtraj),COMnlreg, '-', 'Color', [0 1 0], 'LineWidth', 1); hold on
% %         ylim([-80 80]); xlim([1 length(COMtraj)])
% end
% axis(ha, 'square'); 
% axis(ha(6:end), [1 length(COMtraj) -80 80]); % [xmin xmax ymin ymax]
% ha(1).YLabel.String = 'running laps'; ha(1).XLabel.String = 'position on track';
% ha(6).YLabel.String = 'onset-centered COM'; ha(1).XLabel.String = 'laps after onset';

%% fig zigzag (just ratemaps)
clear groupPFs
close all
zigzag1 = [943, 512, 1411, 232, 986, 377]; % main fig 7
zigzag2 = [2180, 1471, 2504, 2428]; % the 2 last ones are CA3. PF2504 and 1471 are plotted in main figure 6
% groupPFs = [943, 956, 512, 1411, 232]; % ranked by PC1 score
groupPFs = zigzag1;

fzigzag2 = figure
for n = 1:length(groupPFs)
    ex = groupPFs(n);
    subplot(2,length(groupPFs), n, 'align')
        COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
        RateMap = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).PFs{Tok.PF(ex)};
        imagesc(RateMap'); hold on
        set(gca, 'XTick', [1 50], 'XTickLabel', ['start'; '  end'], 'YTick', [1 size(RateMap,2)], 'YTickLabel', [' 1'; num2str(size(RateMap,2))]);
        colormap(jet(256)); %colormap colormap(flipud(hot(256)));% colormap(flipud(gray(256)));
%         ylabel('running laps')
%         xlabel('position on track')
        axis square; box off;
        title({['PF' num2str(ex) ' ' Tok.Region{ex}] ;['lin R2 = ' num2str(Tok.COM_R2(ex))]; ['nl R2 = ' num2str(Tok.RCreg_R2(ex))]; ['PC1 = ' num2str(score(ex,1))]})
    subplot(2,length(groupPFs),length(groupPFs)+n, 'align')
        onset = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).onset_lap(Tok.PF(ex));
        COMtraj = (COMbin(onset:end) - COMbin(onset)).*6;
        LapT = num2str(length(COMtraj)-1);
        scatter(COMtraj,0:length(COMtraj)-1, 15, 'k', 'filled'); hold on 
        xlim([-100 100]); 
        ylim([0 length(COMtraj)])
        set(gca, 'YDir','reverse', 'YTick', [0 length(COMtraj)], 'YTickLabel', [' 0'; LapT], 'XTick', [-50 0 50])
%         ylabel('laps after onset'); 
%         ylabel('onset-centered COM position (cm)');
        axis square; box off; 
        clear LapT
%         title({['linreg Slope = ' num2str(Tok.COM_slope(ex))]; ['nlreg Amp (cm) = ' num2str(Tok.RCreg_p1(ex))]; ['nlreg tau (laps) = ' num2str(Tok.RCreg_p2(ex))]})
end
% fzigzag2.Renderer = 'painters';
% print(fzigzag2, '-vector','-dpdf','exPFsZigZag6PFsLine.pdf')

figure % zigzag trajectories for Graphical Summary
tiledlayout('flow', 'TileSpacing','Compact','Padding','Compact')
for n = [1,3,5]
    nexttile
    ex = groupPFs(n);
    COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
        onset = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).onset_lap(Tok.PF(ex));
        COMtraj = (COMbin(onset:end) - COMbin(onset)).*6;
        % Interpolate NaNs
        x = COMtraj; nanx = isnan(x); laps = 1:numel(x);
        x(nanx) = interp1(laps(~nanx), x(~nanx), laps(nanx));
        COMtrajSmooth = smooth(x,5);
        plot(0:length(COMtraj)-1, COMtrajSmooth, 'k-'); hold on 
        ylim([-100 100]); 
        xlim([0 length(COMtraj)])
        axis square; box off; 
        clear LapT x
end


% figure
% tiledlayout('flow', 'TileSpacing','Compact','Padding','Compact')
% for n = 1:length(groupPFs)
%     ex = groupPFs(n);
%     nexttile
%         COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
%         RateMap = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).PFs{Tok.PF(ex)};
%         imagesc(RateMap'); hold on
%         set(gca, 'XTick', [1 50], 'XTickLabel', ['start'; '  end'], 'YTick', [1 size(RateMap,2)], 'YTickLabel', [' 1'; num2str(size(RateMap,2))]);
%         colormap(jet(256)); %colormap colormap(flipud(hot(256)));% colormap(flipud(gray(256)));
%         ylabel('running lap #')
%         xlabel('Track Position')
%         axis square; box off;
% %         title({['PF' num2str(ex) ' ' Tok.Region{ex}] ;['lin R2 = ' num2str(Tok.COM_R2(ex))]; ['nl R2 = ' num2str(Tok.RCreg_R2(ex))]; ['PC1 = ' num2str(score(ex,1))]})
%     nexttile
%         onset = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).onset_lap(Tok.PF(ex));
%         COMtraj = (COMbin(onset:end) - COMbin(onset)).*6;
%         LapT{n} = num2str(length(COMtraj)-1);
%         scatter(COMtraj,0:length(COMtraj)-1, 15, 'k', 'filled'); hold on 
%         xlim([-100 100]); 
%         ylim([1 length(COMtraj)])
%         set(gca, 'YDir','reverse', 'YTick', [0 length(COMtraj)], 'YTickLabel', [' 0'; LapT{n}], 'XTick', [-50 0 50])
%         ylabel('Laps after onset'); 
%         xlabel('COM (cm)');
%         axis square; box off; 
% %         title({['linreg Slope = ' num2str(Tok.COM_slope(ex))]; ['nlreg Amp (cm) = ' num2str(Tok.RCreg_p1(ex))]; ['nlreg tau (laps) = ' num2str(Tok.RCreg_p2(ex))]})
% 
% end

% figure % just rate maps, with COM overlayed to show how it zigzags
% for n = 1:length(groupPFs)
%     ex = groupPFs(n);
%     subplot(2,3, n, 'align')
%         COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
%         RateMap = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).PFs{Tok.PF(ex)};
%         imagesc(RateMap'); hold on
%         set(gca, 'XTick', [1 50], 'XTickLabel', ['start'; '  end'], 'YTick', [1 size(RateMap,2)], 'YTickLabel', [' 1'; num2str(size(RateMap,2))]);
%         colormap(jet(256)); %colormap colormap(flipud(hot(256)));% colormap(flipud(gray(256)));
%         scatter(COMbin,1:length(COMbin), 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 0.2 1]);
% %         plot(COMbin,1:length(COMbin),'w.-', 'Linewidth', 2);
% %         ylabel('running laps')
% %         xlabel('position on track')
%         axis square; box off;
%         title({['PF' num2str(ex) ' ' Tok.Region{ex}] ;['lin R2 = ' num2str(Tok.COM_R2(ex))]; ['nl R2 = ' num2str(Tok.RCreg_R2(ex))]; ['PC1 = ' num2str(score(ex,1))]})
% end


%% fig linear and stable PFs
clear groupPFs
% close all

groupPFs = [511, 2261, 1326];

flin = figure
for n = 1:length(groupPFs)
    ex = groupPFs(n);
    subplot(2,length(groupPFs), n, 'align')
        COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
        RateMap = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).PFs{Tok.PF(ex)};
        imagesc(RateMap'); hold on
        set(gca, 'XTick', [1 50], 'XTickLabel', ['start'; '  end'], 'YTick', [1 size(RateMap,2)], 'YTickLabel', [' 1'; num2str(size(RateMap,2))]);
        colormap(jet(256)); %colormap colormap(flipud(hot(256)));% colormap(flipud(gray(256)));
%         ylabel('running laps')
%         xlabel('position on track')
        axis square; box off;
        title({['PF' num2str(ex) ' ' Tok.Region{ex}] ; ['Slope = ' num2str(Tok.COM_slope(ex))]; ['lin R2 = ' num2str(Tok.COM_R2(ex))]; ['p = ' num2str(Tok.COM_pval(ex))]})

    subplot(2,length(groupPFs),length(groupPFs)+n, 'align')
        COMa = Tok.COM_slope(ex);
        Pa = Tok.COM_pval(ex);
        if Pa < 0.05 && COMa > 0
            colorp = red;
        elseif Pa < 0.05 && COMa < 0
            colorp = blue;
        else
            colorp = grey;
        end
        COMb= Tok.COM_intcp(ex);
        onset = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).onset_lap(Tok.PF(ex));
        COMtraj = (COMbin(onset:end) - COMbin(onset)).*6;
        COMlinreg = COMb + COMa.*[1:length(COMtraj)];
        scatter(1:length(COMtraj),COMtraj, 30, colorp, 'filled', 'MarkerFaceAlpha', 0.5); hold on 
        plot(1:length(COMtraj),COMlinreg, '-', 'Color', colorp, 'LineWidth', 1); hold on
        xlabel('laps after onset'); 
%         ylabel('onset-centered COM position (cm)');
        ylim([-80 80]); xlim([1 length(COMtraj)])
        axis square; box off;
end

flin.Renderer = 'painters';
print(flin, '-vector','-dpdf','exPFsLinear.pdf')

%% Splitting PFs
clear groupPFs
close all

Splitting = [1242, 429, 608,1450, 797, 712, 1327, 1019, 318, 2300]; % ex splitting PFs (discovered as PFs with increasing PFsd + PFs with high values of PC2,3 or 4)
groupPFs = Splitting;


fsplitting = figure
for n = 1:length(groupPFs)
    ex = groupPFs(n);
    subplot(2,5, n, 'align')
        COMbin = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).COM(Tok.PF(ex),:);
        RateMap = COMdata.(Tok.Region{ex}).(Tok.f_n{ex}).(Tok.File{ex}).PFs{Tok.PF(ex)};
        imagesc(RateMap'); hold on
        set(gca, 'XTick', [1 50], 'XTickLabel', ['start'; '  end'], 'YTick', [1 size(RateMap,2)], 'YTickLabel', [' 1'; num2str(size(RateMap,2))]);
        colormap(jet(256)); %colormap colormap(flipud(hot(256)));% colormap(flipud(gray(256)));
%         ylabel('running laps')
%         xlabel('position on track')
        axis square; box off;
        title({['PF' num2str(ex) ' ' Tok.Region{ex}] ;['lin R2 = ' num2str(Tok.COM_R2(ex))]; ['nl R2 = ' num2str(Tok.RCreg_R2(ex))]; ['PC1 = ' num2str(score(ex,1))]})
end

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
                
                  % modelfun = fittype(@(p1,p2,p3,x) p1*(1-exp(-x/p2)) + p3);
                  modelfun = fittype(@(p1,p2,x) p1*(1-exp(-x/p2)));
                  options = fitoptions('Method','NonlinearLeastSquares', 'Algorithm', 'Trust-Region');
                if b2(2)>0 % positive slope from linear regression
                % params0 = [14 2 0];
                % options.Lower = [0,0.1,-25];
                % options.Upper = [200,100,25];
                params0 = [14 2];
                options.Lower = [0,0.1];
                options.Upper = [200,100];
                else
                % params0 = [-15 2 0];
                % options.Lower = [-200,0.1,-25];
                % options.Upper = [0,100,25];
                params0 = [-15 2];
                options.Lower = [-200,0.1];
                options.Upper = [0,100];
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
                T.RCreg_p3 = [T.RCreg_p3; 0]; %[T.RCreg_p3; mdl.p3];
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
% writetable(Tall,filename);

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
%figure; scatter(1:4,1:4,200,mapCA, 'filled')

% fwd PFs should be red, Bck PFs blue, NonSig grey or black
blue = [0 0.45 0.74]; % backward
red = [0.85 0.32 0.01]; % forward
grey = [0.5 0.5 0.5]; %nonsig

% BTSP models color code
% mapBTSP = map(9:10,:);
cline = lines(7);
mapBTSP(2,:) = cline(4,:);
mapBTSP(1,:) = cline(4,:)+0.25;

%% All Trajectories and Mean Squared Displacement analysis

figure % plot all trajectories
xlaps = repmat(1:Minlaps, OkTraj, 1);
plot(xlaps',COMtrajMat2', 'k-'); hold on
plot(1:Minlaps, BackwardTraj, 'Color', blue, 'linewidth', 2);
plot(1:Minlaps, ForwardTraj, 'Color', red, 'linewidth', 2);
yline(0,'k--');
ylim([-200 200]);
xlabel('lap'); ylabel('COM position (cm)');
title('PFs trajectories')
box off; %axis square;
% print('-vector','-dpdf','PFtrajsAll')


% mean squared displacement of CA1F and N, CA3N and F
BTSPmdl1 = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_Speed15_Pamp20pA_pCS0_005_500repeats.mat');
Trajs = cat(1,BTSPmdl1.COMbin{:}).*6; % in cm
Disp = Trajs - Trajs(:,1); % trajectories centered on initial location, with forward direction as positive values
DispSq = Disp.^2;
MSD = mean(DispSq,1); % average all neurons for each lap
rootMSD = sqrt(MSD);
for lap = 1:size(Trajs,2)
msdCI(:,lap) = bootci(1000, @mean, DispSq(:,lap));
end

BTSPmdl2 = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_Speed15_Pamp20pA_noPlasticity_100repeats.mat');
Trajs2 = cat(1,BTSPmdl2.COMbin{:}).*6; % in cm
Disp2 = Trajs2 - Trajs2(:,1); % trajectories centered on initial location, with forward direction as positive values
DispSq2 = Disp2.^2;
MSD2 = mean(DispSq2,1); % average all neurons for each lap
rootMSD2 = sqrt(MSD2);
for lap = 1:size(Trajs2,2)
msd2CI(:,lap) = bootci(1000, @mean, DispSq2(:,lap));
end

BTSPmdl3 = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_BTSP100repeats_Speed15_Pamp20pA_pCS0_02.mat');
Trajs3 = cat(1,BTSPmdl3.COMbin{:}).*6; % in cm
Disp3 = Trajs3 - Trajs3(:,1); % trajectories centered on initial location, with forward direction as positive values
DispSq3 = Disp3.^2;
MSD3 = mean(DispSq3,1); % average all neurons for each lap
rootMSD3 = sqrt(MSD3);
% for lap = 1:size(Trajs3,2)
% msd3CI(:,lap) = bootci(1000, @mean, DispSq3(:,lap));
% end

BTSPmdl4 = load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\BTSP\100samples\instantBTSPdecoupledHomeonorm\Speed15cmpersec_Pamp20pA\Workspace_Speed15_Pamp20pA_pCS0_002_500repeats.mat');
Trajs4 = cat(1,BTSPmdl4.COMbin{:}).*6; % in cm
Disp4 = Trajs4 - Trajs4(:,1); % trajectories centered on initial location, with forward direction as positive values
DispSq4 = Disp4.^2;
MSD4 = mean(DispSq4,1); % average all neurons for each lap
rootMSD4 = sqrt(MSD4);
for lap = 1:size(Trajs4,2)
msd4CI(:,lap) = bootci(1000, @mean, DispSq4(:,lap));
end

% TrajsCA1F = COMtrajMatNaN2(CA1_F_idx,:); % no interpolation (i.e. raw data). Use COMtrajMat2 for interpolated trajs
% TrajsCA1N = COMtrajMatNaN2(CA1_N_idx,:);
% TrajsCA3F = COMtrajMatNaN2(CA3_F_idx,:);
% TrajsCA3N = COMtrajMatNaN2(CA3_N_idx,:);

TrajsCA1F = COMtrajMat2(CA1_F_idx,:); % with interpolation
TrajsCA1N = COMtrajMat2(CA1_N_idx,:);
TrajsCA3F = COMtrajMat2(CA3_F_idx,:);
TrajsCA3N = COMtrajMat2(CA3_N_idx,:);

DispSqCA1N = ( TrajsCA1N-TrajsCA1N(:,1) ).^2;
DispSqCA1F = ( TrajsCA1F-TrajsCA1F(:,1) ).^2;
DispSqCA3N = ( TrajsCA3N-TrajsCA3N(:,1) ).^2;
DispSqCA3F = ( TrajsCA3F-TrajsCA3F(:,1) ).^2;

msdCA1F = mean( DispSqCA1F , 1 );
msdCA1N = mean( DispSqCA1N, 1 );
msdCA3F = mean( DispSqCA3F, 1 );
msdCA3N = mean( DispSqCA3N, 1 );

DispSqCAall = {DispSqCA1N; DispSqCA1F; DispSqCA3N; DispSqCA3F};
MSDc = [msdCA1N', msdCA1F', msdCA3N', msdCA3F'];

figure
subplot(1,4,1)
    plot(0:Minlaps-1, msdCA3F,'-', 'LineWidth',1,'Color', mapCA3(1,:)); hold on 
    plot(0:Minlaps-1, msdCA3N,'-', 'LineWidth',1,'Color', mapCA3(2,:)); hold on 
    for n = 1:2 %size(MSDc, 2)
        for lap = 1:Minlaps
            msdCA_CI{n}(:,lap) = bootci(1000, @mean, DispSqCAall{n}(:,lap));
        end
        plot_ci(0:Minlaps-1, [MSDc(:,n) msdCA_CI{n}(1,:)' msdCA_CI{n}(2,:)'], 'MainLineColor', colors(n,:), 'PatchColor', colors(n,:), 'LineColor', colors(n,:), 'PatchAlpha', 0.5); hold on
    end
    plot_ci([1:length(MSD)]-1, [MSD' msdCI(1,:)' msdCI(2,:)'], 'MainLineColor', mapBTSP(2,:), 'PatchColor', mapBTSP(2,:), 'LineColor', mapBTSP(2,:), 'PatchAlpha', 0.5); hold on
    plot_ci([1:length(MSD2)]-1, [MSD2' msd2CI(1,:)' msd2CI(2,:)'], 'MainLineColor', 'k', 'PatchColor', 'k', 'LineColor', 'k', 'PatchAlpha', 0.5); hold on
    plot_ci([1:length(MSD4)]-1, [MSD4' msd4CI(1,:)' msd4CI(2,:)'], 'MainLineColor', mapBTSP(1,:), 'PatchColor', mapBTSP(1,:), 'LineColor', mapBTSP(1,:), 'PatchAlpha', 0.5); hold on
    xlabel('Laps after onset')
    ylabel('Mean Squared Displacement (cm^2)')
    box off
    axis square
subplot(1,4,2)
    RegStart = 4;
    Xlm = [ [RegStart:Minlaps]', ones(size([RegStart:Minlaps]'))];
    for n = 1:size(MSDc,2)
    Ylm{n} = MSDc(RegStart:end,n);
    [B{n},BINT{n},R{n},RINT{n},STATS{n}] = regress(Ylm{n}, Xlm);
    lm{n} = Xlm*B{n};
    plot([RegStart:Minlaps]'-1, lm{n}, '-', 'LineWidth', 1.5,'Color', colors(n,:)); hold on 
    scatter([1:Minlaps]-1, MSDc(:,n),'o', 'MarkerFaceColor', colors(n,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', colors(n,:), 'MarkerEdgeAlpha', 0.5); hold on
    text(Minlaps, max(lm{n}), ['D = ' num2str( round(B{n}(1)/2, 2, 'significant') )], 'Color', colors(n,:))
    % text(Minlaps+1, max(lm{n})-30, {['slope = ' num2str( round(B{n}(1), 2, 'significant') )]; ['Rsq = ' num2str( round(STATS{n}(1),3, 'significant')*100) '%' ]; ['p = ' num2str( round(STATS{n}(3),3, 'significant')) ]}, 'Color', colors(n,:))
    end
    [Bmsd,BINTmsd,Rmsd,RINTmsd,STATSmsd] = regress(MSD(RegStart:Minlaps)', Xlm);
    lmMSD = Xlm*Bmsd;
    plot([RegStart:Minlaps]'-1, lmMSD, '-', 'LineWidth', 1.5,'Color', mapBTSP(2,:)); hold on 
    scatter([1:length(MSD)]-1, MSD, 'o', 'MarkerEdgeColor', mapBTSP(2,:), 'MarkerEdgeAlpha', 0.5, 'MarkerFaceColor', mapBTSP(2,:), 'MarkerFaceAlpha', 0.5); hold on
    xlim([0 40])
    text(Minlaps, max(lmMSD), ['D = ' num2str( round(Bmsd(1)/2, 2, 'significant') )], 'Color', mapBTSP(2,:))
    % text(Minlaps+1, max(lmMSD)-30, {['slope = ' num2str( round(Bmsd(1), 2, 'significant') )]; ['Rsq = ' num2str( round(STATSmsd(1),3, 'significant')*100) '%' ]; ['p = ' num2str( round(STATSmsd(3),3, 'significant')) ]}, 'Color', cline(4,:))
    [Bmsd4,BINTmsd4,Rmsd4,RINTmsd4,STATSmsd4] = regress(MSD4(RegStart:Minlaps)', Xlm);
    lmMSD4 = Xlm*Bmsd4;
    plot([RegStart:Minlaps]'-1, lmMSD4, '-', 'LineWidth', 1.5,'Color', mapBTSP(1,:)); hold on 
    scatter([1:length(MSD4)]-1, MSD4, 'o', 'MarkerEdgeColor', mapBTSP(1,:), 'MarkerEdgeAlpha', 0.5, 'MarkerFaceColor', mapBTSP(1,:), 'MarkerFaceAlpha', 0.5); hold on
    xlim([0 40])
    text(Minlaps, max(lmMSD4), ['D = ' num2str( round(Bmsd4(1)/2, 2, 'significant') )], 'Color', mapBTSP(1,:))
    xlabel('Laps after onset')
    ylabel('MSD (cm^2)')
    box off
    axis square
subplot(1,4,3)
    plot(1:length(MSD)-1, diff(MSD)./2,'-', 'LineWidth',1, 'Color', cline(4,:)); hold on 
    plot(1:Minlaps-1, diff(msdCA1F)./2,'-', 'LineWidth',1,'Color', mapCA1(1,:)); hold on 
    plot(1:Minlaps-1, diff(msdCA1N)./2,'-', 'LineWidth',1,'Color', mapCA1(2,:)); hold on 
    plot(1:Minlaps-1, diff(msdCA3F)./2,'-', 'LineWidth',1,'Color', mapCA3(1,:)); hold on 
    plot(1:Minlaps-1, diff(msdCA3N)./2,'-', 'LineWidth',1,'Color', mapCA3(2,:)); hold on 
    xlabel('Laps after onset')
    ylabel('Instant. Diffusion Coeff. (cm^2/lap)')
    box off
    axis square
subplot(1,4,4)
    Xeval = 1:Minlaps;
    Xdiff = 2:Minlaps;
    Xdiff2 = 2:0.1:Minlaps;
    modelfun3 = fittype(@(p1,p2,p3, x) p1*exp(-(x-1)/p2)+p3);
    paramsMSD = [100 2 0];
    options3 = fitoptions('Method','NonlinearLeastSquares', 'Algorithm', 'Trust-Region');
    options3.Lower = [0,0.01,0];
    options3.Upper = [1000,100,20];
    options3.Startpoint = paramsMSD;
    for n = 1:size(MSDc,2)
        Dcoef(:,n) = diff(MSDc(:,n))./2;
        [Dmdl{n},gofD{n},outD{n}] = fit(Xdiff', Dcoef(:,n), modelfun3, options3);
        Deval{n} = feval(Dmdl{n}, Xdiff2);
        scatter(Xdiff-1, Dcoef(:,n),'o','MarkerFaceColor', colors(n,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', colors(n,:), 'MarkerEdgeAlpha', 0.5); hold on
        plot(Xdiff2-1, Deval{n}, '-', 'Color', colors(n,:), 'LineWidth', 2); hold on
    end
    DcoefMSD = diff(MSD(1:Minlaps)')./2;
    [Dmsd,gofDmsd,outDmsd] = fit(Xdiff', DcoefMSD, modelfun3, options3);
    DevalMSD = feval(Dmsd, Xdiff2);
    scatter(Xdiff-1, DcoefMSD,'o','MarkerFaceColor', mapBTSP(2,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', mapBTSP(2,:), 'MarkerEdgeAlpha', 0.5); hold on
    plot(Xdiff2-1, DevalMSD, '-', 'Color', mapBTSP(2,:), 'LineWidth', 2); hold on
    xlabel('Laps after onset')
    ylabel('Instant. Diffusion Coeff. (cm^2/lap)')
    box off
    axis square

figure
    for n = 1:size(MSDc, 2)
%         for lap = 1:Minlaps
%             msdCA_CI{n}(:,lap) = bootci(1000, @mean, DispSqCAall{n}(:,lap));
%         end
%         plot_ci(1:Minlaps, [MSDc(:,n) msdCA_CI{n}(1,:)' msdCA_CI{n}(2,:)'], 'MainLineColor', colors(n,:), 'PatchColor', colors(n,:), 'LineColor', colors(n,:), 'PatchAlpha', 0.5); hold on
        plot(0:Minlaps-1, MSDc(:,n),'-', 'LineWidth',1,'Color', colors(n,:)); hold on 
    end
    plot([1:length(MSD2)]-1, MSD2', 'Color', 'k'); hold on   
    plot([1:length(MSD4)]-1, MSD4', 'Color', mapBTSP(1,:)); hold on
    plot([1:length(MSD)]-1, MSD', 'Color', mapBTSP(2,:)); hold on
    legend('CA1N','CA1F','CA3N','CA3F', 'mdl: no plasticity', 'BTSP mdl: p(CS)=0.02','BTSP mdl: p(CS)=0.05','Location', 'BestOutside')
    xlabel('Laps after onset')
    ylabel('Mean Squared Displacement (cm^2)')
    box off
    axis square

figure
RegStart = 4;
Xlm = [ [RegStart:Minlaps]', ones(size([RegStart:Minlaps]'))];
for n = 1:size(MSDc,2)
Ylm{n} = MSDc(RegStart:end,n);
[B{n},BINT{n},R{n},RINT{n},STATS{n}] = regress(Ylm{n}, Xlm);
lm{n} = Xlm*B{n};
plot([RegStart:Minlaps]', lm{n}, '-', 'LineWidth', 1.5,'Color', colors(n,:)); hold on 
scatter(1:Minlaps, MSDc(:,n),'o', 'MarkerFaceColor', colors(n,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', colors(n,:), 'MarkerEdgeAlpha', 0.5); hold on
text(Minlaps+1, max(lm{n}), ['D = ' num2str( round(B{n}(1)/2, 2, 'significant') )], 'Color', colors(n,:))
% text(Minlaps+1, max(lm{n})-30, {['slope = ' num2str( round(B{n}(1), 2, 'significant') )]; ['Rsq = ' num2str( round(STATS{n}(1),3, 'significant')*100) '%' ]; ['p = ' num2str( round(STATS{n}(3),3, 'significant')) ]}, 'Color', colors(n,:))
end
[Bmsd,BINTmsd,Rmsd,RINTmsd,STATSmsd] = regress(MSD(RegStart:Minlaps)', Xlm);
lmMSD = Xlm*Bmsd;
plot([RegStart:Minlaps]', lmMSD, '-', 'LineWidth', 1.5,'Color', cline(4,:)); hold on 
scatter(1:length(MSD), MSD, 'o', 'MarkerEdgeColor', cline(4,:), 'MarkerEdgeAlpha', 0.5, 'MarkerFaceColor', cline(4,:), 'MarkerFaceAlpha', 0.5); hold on
xlim([0 40])
text(Minlaps+1, max(lmMSD), ['D = ' num2str( round(Bmsd(1)/2, 2, 'significant') )], 'Color', cline(4,:))
% text(Minlaps+1, max(lmMSD)-30, {['slope = ' num2str( round(Bmsd(1), 2, 'significant') )]; ['Rsq = ' num2str( round(STATSmsd(1),3, 'significant')*100) '%' ]; ['p = ' num2str( round(STATSmsd(3),3, 'significant')) ]}, 'Color', cline(4,:))
[Bmsd4,BINTmsd4,Rmsd4,RINTmsd4,STATSmsd4] = regress(MSD4(RegStart:Minlaps)', Xlm);
lmMSD4 = Xlm*Bmsd4;
plot([RegStart:Minlaps]', lmMSD4, '-', 'LineWidth', 1.5,'Color', cline(7,:)); hold on 
scatter(1:length(MSD4), MSD4, 'o', 'MarkerEdgeColor', cline(7,:), 'MarkerEdgeAlpha', 0.5, 'MarkerFaceColor', cline(7,:), 'MarkerFaceAlpha', 0.5); hold on
xlim([0 40])
text(Minlaps+1, max(lmMSD4), ['D = ' num2str( round(Bmsd4(1)/2, 2, 'significant') )], 'Color', cline(7,:))
legend(['CA1N: Rsq = ' num2str( round(STATS{1}(1),3, 'significant')*100) '%, p = ' num2str( round(STATS{1}(3),2, 'significant')) ],'', ...
    ['CA1F: Rsq = ' num2str( round(STATS{2}(1),3, 'significant')*100) '%, p = ' num2str( round(STATS{2}(3),2, 'significant')) ],'', ...
    ['CA3N: Rsq = ' num2str( round(STATS{3}(1),3, 'significant')*100) '%, p = ' num2str( round(STATS{3}(3),2, 'significant')) ],'', ...
    ['CA3F: Rsq = ' num2str( round(STATS{4}(1),3, 'significant')*100) '%, p = ' num2str( round(STATS{4}(3),2, 'significant')) ],'', ...
    ['BTSP p(CS)=0.005: Rsq = ' num2str( round(STATSmsd(1),3, 'significant')*100) '%, p = ' num2str( round(STATSmsd(3),2, 'significant')) ],'',...
    ['BTSP p(CS)=0.002: Rsq = ' num2str( round(STATSmsd4(1),3, 'significant')*100) '%, p = ' num2str( round(STATSmsd4(3),2, 'significant')) ],'','Location', 'BestOutside' )
xlabel('laps')
ylabel('mean squared displacement (cm^2)')
box off
axis square

% bootstrapped regression to check if late portion of MSD really has a slope different from 0
% -> resample the neurons and compute the lapwise MSD for each bootstrap, with associated regression
DispSqCAall_late = {DispSqCA1N(:,RegStart:end); DispSqCA1F(:,RegStart:end); DispSqCA3N(:,RegStart:end); DispSqCA3F(:,RegStart:end)};

nboot = 10000;
replacement = true(1);
for n = 1:4%size(MSDc,2) % for each CA group
    for i = 1:nboot
    ResIdx = randsample(1:size(DispSqCAall{n},1), size(DispSqCAall{n},1), replacement);
    MSDcRes = mean( DispSqCAall{n}(ResIdx,:), 1 );
    YlmRes = MSDcRes(RegStart:end)';
    Bres = regress(YlmRes, Xlm);
    SlopeRes(i,n) = Bres(1);
    InterceptRes(i,n) = Bres(2);
    % lmRes = Xlm*Bres;
    % plot([RegStart:Minlaps]', lmRes, '-', 'LineWidth', 1.5,'Color', colors(n,:)); hold on 
    % scatter(1:Minlaps, MSDcRes,'o', 'MarkerFaceColor', colors(n,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', colors(n,:), 'MarkerEdgeAlpha', 0.5); hold on
    % text(Minlaps+1, max(lmRes), ['D = ' num2str( round(Bres(1)/2, 2, 'significant') )], 'Color', colors(n,:))
    end
end




figure
for n = 1:4
Dall = SlopeRes(:,n)./2;
pD{n} = length(Dall(Dall<0))./length(Dall);
    subplot(1,4,n)
    histogram(SlopeRes(:,n)./2, 'EdgeColor', colors(n,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    DiffCoeff_CI{n} = prctile(SlopeRes(:,n)./2,[2.5, 97.5]);
    errorbar(mean(SlopeRes(:,n)./2), 0.01, mean(SlopeRes(:,n)./2)-DiffCoeff_CI{n}(1),mean(SlopeRes(:,n)./2)-DiffCoeff_CI{n}(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 6); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-4 8]);
    xlabel('D (Diff Coef)'); ylabel('proba');
    title(['pval = ' num2str(pD{n})])
    set(gca,'view',[90 -90])
    box off; %axis square;
end

% figure
% for n = 1:2%size(MSDc,2) % for each CA group
% Yall = DispSqCAall_late{n}(:);
% Xall = [ repmat( [RegStart:Minlaps]', size(DispSqCAall_late{n},1) , 1 ), ones( size(Yall) ) ];
% [Ball{n},BINTall{n},Rall{n},RINTall{n},STATSall{n}] = regress(Yall, Xall);
% lmAll{n} = Xlm*Ball{n};
% plot([RegStart:Minlaps]', lmAll{n}, '-', 'LineWidth', 2,'Color', 'k'); hold on 
% plot_ci(0:Minlaps-1, [MSDc(:,n) msdCA_CI{n}(1,:)' msdCA_CI{n}(2,:)'], 'MainLineColor', colors(n,:), 'PatchColor', colors(n,:), 'LineColor', colors(n,:), 'PatchAlpha', 0.5); hold on
% end
% xlabel('laps')
% ylabel('mean squared displacement (cm^2)')
% box off
% axis square

figure % instant diffusion coef D, estimated from derivative of msd and Einstein's equation msd = 2D for brownian mvmt (the 2 comes from the 2 possible directions for random walk in 1 dimension)
plot(2:length(MSD), diff(MSD)./2,'-', 'LineWidth',1, 'Color', cline(4,:)); hold on 
plot(2:Minlaps, diff(msdCA1F)./2,'-', 'LineWidth',1,'Color', mapCA1(1,:)); hold on 
% plot(Xdiff, 150.*exp(-(Xdiff-1)/2))
plot(2:Minlaps, diff(msdCA1N)./2,'-', 'LineWidth',1,'Color', mapCA1(2,:)); hold on 
plot(2:Minlaps, diff(msdCA3F)./2,'-', 'LineWidth',1,'Color', mapCA3(1,:)); hold on 
plot(2:Minlaps, diff(msdCA3N)./2,'-', 'LineWidth',1,'Color', mapCA3(2,:)); hold on 
% xlim([0 15])
xlabel('laps')
ylabel('Instantaneous Diffusion coeff (cm^2/lap)')
box off
axis square

figure % Dcoef fitted with exponential, to smooth and figure out at what lap it stabilizes
% subplot(1,2,1)
    Xeval = 1:Minlaps;
    Xdiff = [2:Minlaps]-1;
    Xdiff2 = [2:0.1:Minlaps]-1;
    modelfun3 = fittype(@(p1,p2,p3, x) p1*exp(-(x-1)/p2)+p3);
    paramsMSD = [100 2 0];
    options3 = fitoptions('Method','NonlinearLeastSquares', 'Algorithm', 'Trust-Region');
    options3.Lower = [0,0.01,0];
    options3.Upper = [1000,100,20];
    options3.Startpoint = paramsMSD;
    for n = 1:size(MSDc,2)
        Dcoef(:,n) = diff(MSDc(:,n))./2;
        [Dmdl{n},gofD{n},outD{n}] = fit(Xdiff', Dcoef(:,n), modelfun3, options3);
        Deval{n} = feval(Dmdl{n}, Xdiff2);
        scatter(Xdiff, Dcoef(:,n),'o','MarkerFaceColor', colors(n,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', colors(n,:), 'MarkerEdgeAlpha', 0.5); hold on
        plot(Xdiff2, Deval{n}, '-', 'Color', colors(n,:), 'LineWidth', 2); hold on
    end
    DcoefMSD = diff(MSD(1:Minlaps)')./2;
    [Dmsd,gofDmsd,outDmsd] = fit(Xdiff', DcoefMSD, modelfun3, options3);
    DevalMSD = feval(Dmsd, Xdiff2);
    scatter(Xdiff, DcoefMSD,'o','MarkerFaceColor', cline(4,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', cline(4,:), 'MarkerEdgeAlpha', 0.5); hold on
    plot(Xdiff2, DevalMSD, '-', 'Color', cline(4,:), 'LineWidth', 2); hold on
    xlabel('laps')
    ylabel('Instantaneous Diffusion coeff (cm^2/lap)')
    box off
    axis square
% subplot(1,2,2) % zoom
figure
    for n = 1:size(MSDc,2)
        scatter(Xdiff, Dcoef(:,n),'o','MarkerFaceColor', colors(n,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', colors(n,:), 'MarkerEdgeAlpha', 0.5); hold on
        plot(Xdiff2, Deval{n}, '-', 'Color', colors(n,:), 'LineWidth', 2); hold on
        yline(Dmdl{n}.p3, '--', 'Color', colors(n,:), 'LineWidth', 1); hold on
    end
    scatter(Xdiff, DcoefMSD,'o','MarkerFaceColor', cline(4,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', cline(4,:), 'MarkerEdgeAlpha', 0.5); hold on
    plot(Xdiff2, DevalMSD, '-', 'Color', cline(4,:), 'LineWidth', 2); hold on
    yline(Dmsd.p3, '--', 'Color', mapBTSP(2,:), 'LineWidth', 1)
    legend('',['CA1N: D = ' num2str( round(Dmdl{1}.p3,2, 'significant' ) ) ', Rsq = ' num2str( round(gofD{1}.rsquare,3, 'significant')*100) '%' ],'', ...
    '',['CA1F: D = ' num2str( round(Dmdl{2}.p3,2, 'significant' ) ) ', Rsq = ' num2str( round(gofD{2}.rsquare,3, 'significant')*100) '%' ],'', ...
    '',['CA3N: D = ' num2str( round(Dmdl{3}.p3,2, 'significant' ) ) ', Rsq = ' num2str( round(gofD{3}.rsquare,3, 'significant')*100) '%' ],'', ...
    '',['CA3F: D = ' num2str( round(Dmdl{4}.p3,2, 'significant' ) ) ', Rsq = ' num2str( round(gofD{4}.rsquare,3, 'significant')*100) '%' ],'', ...
    '',['BTSP mdl: D = ' num2str( round(Dmsd.p3,2, 'significant' ) ) ', Rsq = ' num2str( round(gofDmsd.rsquare,3, 'significant')*100) '%' ],'', 'Location', 'BestOutside' )
    ylim([0 13])
    xlim([2 15])
    xlabel('Laps after onset')
    ylabel('Instant. Diffusion Coeff (cm^2/lap)')
    box off 
    grid on
    axis square

% % fit nonlinear curve to msd and get derivative? Need to find the correct
% % nl function...
% 
% % exponential saturating as for individual trajectories
% paramsMSD = [50 1 0];
% options.Lower = [0,0.1,-25];
% options.Upper = [1000,100,25];
% options.Startpoint = params0;
% [mdlCA3N,gofCA3N,foutputCA3N] = fit([1:Minlaps]',msdCA3N',modelfun,options);
% [mdlCA3F,gofCA3F,foutputCA3F] = fit([1:Minlaps]',msdCA3F',modelfun,options);
% [mdlCA1N,gofCA1N,foutputCA1N] = fit([1:Minlaps]',msdCA1N',modelfun,options);
% [mdlCA1F,gofCA1F,foutputCA1F] = fit([1:Minlaps]',msdCA1F',modelfun,options);
% 
% figure
% plot(mdlCA3N, 'k'); hold on
% plot(1:Minlaps, msdCA3N,'o','Color', mapCA3(2,:)); hold on
% plot(mdlCA3F, 'k'); hold on
% plot(1:Minlaps, msdCA3F,'o','Color', mapCA3(1,:)); hold on
% plot(mdlCA1N, 'k'); hold on
% plot(1:Minlaps, msdCA1N,'o','Color', mapCA1(2,:)); hold on
% plot(mdlCA1F, 'k'); hold on
% plot(1:Minlaps, msdCA1F,'o','Color', mapCA1(1,:)); hold on
% xlabel('laps')
% ylabel('mean squared displacement (cm^2)')
% box off
% axis square
% 
% % smoothdata
% figure 
% plot(1:size(Trajs,2), smoothdata(MSD,2),'b-', 'LineWidth', 1); hold on % BTSP 
% plot(1:Minlaps, smoothdata(msdCA1F,2),'-', 'LineWidth',1,'Color', mapCA1(1,:)); hold on 
% plot(1:Minlaps, smoothdata(msdCA1N,2),'-', 'LineWidth',1,'Color', mapCA1(2,:)); hold on 
% plot(1:Minlaps, smoothdata(msdCA3F,2),'-', 'LineWidth',1,'Color', mapCA3(1,:)); hold on 
% plot(1:Minlaps, smoothdata(msdCA3N,2),'-', 'LineWidth',1,'Color', mapCA3(2,:)); hold on 
% % xlim([1 size(Trajs,2)+1])
% xlabel('laps')
% ylabel('mean squared displacement (cm^2)')
% box off
% axis square
% 
% figure % diffusion coef estimated from derivative of smoothed msd
% plot(2:Minlaps, diff(smoothdata(MSD,2))./2,'-', 'LineWidth',1, 'Color', cline(4,:)); hold on 
% plot(2:Minlaps, diff(smoothdata(msdCA1F))./2,'-', 'LineWidth',1,'Color', mapCA1(1,:)); hold on 
% plot(2:Minlaps, diff(smoothdata(msdCA1N))./2,'-', 'LineWidth',1,'Color', mapCA1(2,:)); hold on 
% plot(2:Minlaps, diff(smoothdata(msdCA3F))./2,'-', 'LineWidth',1,'Color', mapCA3(1,:)); hold on 
% plot(2:Minlaps, diff(smoothdata(msdCA3N))./2,'-', 'LineWidth',1,'Color', mapCA3(2,:)); hold on 
% xlabel('laps')
% ylabel('Instantaneous Diffusion coeff (cm^2/lap)')
% box off
% axis square
% 
% % smoothing splines
% sp = 0.8;
% [MDSsmooth, gofMDS, outMDS]  = fit([1:size(Trajs,2)]', MSD', 'smoothingspline', 'SmoothingParam', sp);
% [mdsCA1Fsmooth, gofCA1Fs, outCA1Fs] = fit([1:Minlaps]', msdCA1F', 'smoothingspline', 'SmoothingParam', sp);
% [mdsCA1Nsmooth, gofCA1Ns, outCA1Ns] = fit([1:Minlaps]', msdCA1N', 'smoothingspline', 'SmoothingParam', sp);
% [mdsCA3Fsmooth, gofCA3Fs, outCA3Fs] = fit([1:Minlaps]', msdCA3F', 'smoothingspline', 'SmoothingParam', sp);
% [mdsCA3Nsmooth, gofCA3Ns, outCA3Ns] = fit([1:Minlaps]', msdCA3N', 'smoothingspline', 'SmoothingParam', sp);
% 
% figure 
% plot(1:size(Trajs,2), MSD,'ob', 'LineWidth', 1); hold on % BTSP 
% plot(1:Minlaps, msdCA1F,'o', 'LineWidth',1,'Color', mapCA1(1,:)); hold on 
% plot(1:Minlaps, msdCA1N,'o', 'LineWidth',1,'Color', mapCA1(2,:)); hold on 
% plot(1:Minlaps, msdCA3F,'o', 'LineWidth',1,'Color', mapCA3(1,:)); hold on 
% plot(1:Minlaps, msdCA3N,'o', 'LineWidth',1,'Color', mapCA3(2,:)); hold on 
% plot(MDSsmooth, 'b'); hold on 
% plot(mdsCA1Fsmooth); hold on 
% plot(mdsCA1Nsmooth); hold on 
% plot(mdsCA3Fsmooth); hold on 
% plot(mdsCA3Nsmooth); hold on 
% % xlim([1 size(Trajs,2)+1])
% xlabel('laps')
% ylabel('mean squared displacement (cm^2)')
% box off
% axis square
% 
% figure % instant diffusion coef D, estimated from derivative of msd spline fit
% % Xeval = [1.9:0.1:Minlaps]';
% % Xdiff = 2:0.1:Minlaps;
% plot(Xdiff, diff(feval(MDSsmooth,Xeval))./2,'-', 'LineWidth',1, 'Color', cline(4,:)); hold on 
% plot(Xdiff, diff(feval(mdsCA1Fsmooth,Xeval))./2,'-', 'LineWidth',1,'Color', mapCA1(1,:)); hold on 
% plot(Xdiff, diff(feval(mdsCA1Nsmooth,Xeval))./2,'-', 'LineWidth',1,'Color', mapCA1(2,:)); hold on 
% plot(Xdiff, diff(feval(mdsCA3Fsmooth,Xeval))./2,'-', 'LineWidth',1,'Color', mapCA3(1,:)); hold on 
% plot(Xdiff, diff(feval(mdsCA3Nsmooth,Xeval))./2,'-', 'LineWidth',1,'Color', mapCA3(2,:)); hold on 
% xlim([0 30])
% xlabel('laps')
% ylabel('Instantaneous Diffusion coeff (cm^2/lap)')
% box off
% axis square

%% PCA + kmeans clustering, tSNE, Averages all converge to same nonlinear template trajectory
% PCA, tsne, clustering

[coeff,score,latent,tsquared,explained,mu] = pca(COMtrajMat2, 'centered', false); % PCA on interpolated data
% [coeff,score,latent,tsquared,explained,mu] = pca(COMtrajMatNaN2, 'centered', false, 'Algorithm', 'als'); % PCA on raw data without interpolation, without tossing trajectories with a few missing values

Tpca = [COMtrajMat2, score(:,1), score(:,2), score(:,3)];

for nk = 1:20 
[Kidx(:,nk),K_C{nk},sumd{nk},K_D{nk}] = kmeans(score, nk, 'Display', 'final');
totsumd(nk) = sum(sumd{nk},'all'); % within-cluster sums of point-to-centroid distances in the k-by-1 vector sumd
% totKD(nk) = sum(K_D{nk},'all'); %distances from each point to every centroid in the n-by-k matrix D
end

[TSNEpca, loss_pca] = tsne(score);
[TSNEraw, loss_raw] = tsne(COMtrajMat2);

figure
clear Xnl
Xnl = [1:Minlaps]'-1;
YnlF = ForwardTraj';
YnlB = BackwardTraj';
modelfun = fittype(@(p1,p2,p3,x) p1*(1-exp(-x/p2)) + p3);
params0F = [14 2 0];
optionsFlow = [0,0.1,-25];
optionsFUp = [200,100,25];
params0B = [-15 2 0];
optionsBlow = [-200,0.1,-25];
optionsBUp = [0,100,25];
% modelfun = fittype(@(p1,p2,x) p1*(1-exp(-x/p2)));
% params0F = [14 2];
% optionsFlow = [0,0.1]
% optionsFUp = [200,100];
% params0B = [-15 2];
% optionsBlow = [-200,0.1];
% optionsBUp = [0,100];
% title('positive and negative averages VS A*(1-exp(x/tau))')
clear options %F
options = fitoptions('Method','NonlinearLeastSquares', 'Algorithm', 'Trust-Region');
options.Lower = optionsFlow;
options.Upper = optionsFUp;
options.Startpoint = params0F;
[mdlF,gofF,foutputF] = fit(Xnl,YnlF,modelfun,options);
clear options %B
options = fitoptions('Method','NonlinearLeastSquares', 'Algorithm', 'Trust-Region');
options.Lower = optionsBlow;
options.Upper = optionsBUp;
options.Startpoint = params0B;
[mdlB,gofB,foutputB] = fit(Xnl,YnlB,modelfun,options);

% RC = 14*(1-exp(-Xnl./2));
% RC2 = -15*(1-exp(-Xnl./2));
plot([1:0.1:Minlaps]-1, feval(mdlF,[1:0.1:Minlaps]-1), 'Color', red, 'linewidth', 2); hold on
plot([1:0.1:Minlaps]-1, feval(mdlB,[1:0.1:Minlaps]-1), 'Color', blue, 'linewidth', 2); hold on
% plot(Xnl, mdl.Fitted, 'r'); hold on
% plot([1:0.1:Minlaps], RC, red); hold on
plot(0:Minlaps-1, ForwardTraj,'o', 'Color', red, 'linewidth', 2);
% plot([1:0.1:Minlaps], RC2, blue); hold on
plot(0:Minlaps-1, BackwardTraj, 'o', 'Color', blue, 'linewidth', 2);
xlabel('laps from PF onset')
ylabel('COM shift (cm)')
axis square; box off;
title('positive and negative averages VS A*(1-exp(x/tau))+p3')
clear Xnl

figure % plot trajectories as matrix image
[sortedCOMmean, COMsortedIndex] = sort(COMmean, 'descend');
COMtrajMatOrd = COMtrajMat2(COMsortedIndex,:);
imagesc(COMtrajMatOrd);
xlabel('Laps'); ylabel('Place Fields');
colormap(brewermap(256,'*RdBu')); 
caxis([-50, 50]);
xticks([1 15])
xticklabels({'0','14'})
c = colorbar; c.Label.String = 'COM location from onset position (cm)';
box off; 
axis square;
title('PFs trajectories')

figure % variance explained by PCs
bar(explained);
xlabel('PC#'); ylabel('Variance explained');
box off; axis square;

figure % plot 3 first PCs
PCsToPlot = 3;
xlapsPC = repmat(1:Minlaps, PCsToPlot,1);
plot(xlapsPC',coeff(:,1:PCsToPlot))
yline(0,'k--');
% ylim([-150 150]);
xlabel('lap'); ylabel('COM position (cm)');
legend('PC1', 'PC2', 'PC3');
title('Principal Components')
box off; axis square;

figure % plot all PCs
for p = 1:Minlaps
subplot(5,Minlaps/5,p)
plot(0:Minlaps-1,coeff(:,p), 'k')
yline(0,'k--');
% ylim([-150 150]);
% xlabel('lap'); ylabel('COM position (cm)');
title(['PC' num2str(p)])
box off;
end

figure % project trajectories in 3 PC space
scatter3(score(:,1),score(:,2),score(:,3));
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
box off; axis square;

% distribution of all PFs projection on each PC (to see if some are multimodal)
figure % plot all PCs
for p = 1:Minlaps
subplot(5,3,p)
histogram(score(:,p), 'normalization', 'pdf')
xlabel(['PC' num2str(p) ' proj']); ylabel('pdf');
box off;
end

% figure % plot all trajectories reconstructed from pca, to check it works
% xlaps = repmat(1:Minlaps, OkTraj, 1);
% plot(xlaps',[score*coeff']', 'k-')
% yline(0,'k--');
% ylim([-300 300]);
% xlabel('lap'); ylabel('COM position (cm)');
% title('PFs trajectories')
% box off; axis square;

% % plot average Traj for PFs with negative or positive PC1 score
% figure
% scoreOK = score(score(:,1)<0,:);
% reconstruct = [scoreOK*coeff'];
% xlaps = repmat(1:Minlaps, size(scoreOK,1), 1);
% plot(xlaps',reconstruct', 'k-'); hold on
% plot(1:Minlaps, mean(reconstruct,1), 'r');
% yline(0,'k--');
% ylim([-300 300]);
% xlabel('lap'); ylabel('COM position (cm)');
% title('PFs trajectories for criteria')
% box off; axis square;

figure % kmeans clustering on pca scores: finding number of clusters
plot(1:20,totsumd,'b'); hold on
xlabel('number of clusters')
ylabel('total sum of within-cluster distances')
title('kmeans clustering goodness-of-fit')
axis square; box off;

cmax = 6; % number of clusters
cmapk = lines(cmax);

figure % how do the clusters look like
gscatter(score(:,1),score(:,2), Kidx(:,cmax), cmapk)
xlabel('PC1'); ylabel('PC2');
box off; axis square;

% plot average traj for each Kmeans cluster
figure
xlaps = repmat(1:Minlaps, cmax, 1);
for c = 1:cmax
scoreOK = score(Kidx(:,cmax)==c,:);
reconstruct = [scoreOK*coeff'];
TrajC(:,c) = mean(reconstruct,1);
% plot(1:Minlaps,TrajC(:,c)); hold on
clear scoreOK reconstruct
end
plot(xlaps',TrajC, 'LineWidth', 2); hold on
yline(0,'k--');
% ylim([-300 300]);
xlabel('lap'); ylabel('COM position (cm)');
title('average PFs trajectories by clusters')
box off; axis square;

% PC1 vs PC2 for all CA_VR groups
figure
gscatter(score(:,1),score(:,2), GroupCA_VR, colors)
legend('CA1_N','CA1_F','CA3_N', 'CA3_F');
xlabel('PC1'); ylabel('PC2');
box off; axis square;

figure % PC1 violin plots. 
v2 = violinplot(score(:,1), GroupCA_VR);
v2(1).ViolinColor = mapCA1(2,:);
v2(2).ViolinColor = mapCA1(1,:);
v2(3).ViolinColor = mapCA3(2,:);
v2(4).ViolinColor = mapCA3(1,:);
xlabel('CA-VR groups'); ylabel('PC1');
set(gca, 'XTick', [1:4], 'XTickLabel', ['CA1N';'CA1F'; 'CA3N'; 'CA3F']);
box off; axis square;

figure % distrib of 
% histogram(score(CA1_N_idx,1), 'FaceColor', mapCA1(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'Normalization', 'cdf'); hold on
% histogram(score(CA3_N_idx,1), 'FaceColor', mapCA3(2,:), 'EdgeColor', 'none','FaceAlpha', 0.5, 'Normalization', 'cdf'); hold on
histogram(score(CA1_N_idx,1),30, 'EdgeColor', mapCA1(2,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf');hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(score(CA3_N_idx,1),30, 'EdgeColor', mapCA3(2,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf'); hold on
histogram(score(CA1_F_idx,1),30, 'EdgeColor', mapCA1(1,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(score(CA3_F_idx,1),30, 'EdgeColor', mapCA3(1,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
% legend('CA1_N','CA3_N','CA1_F', 'CA3_F');
xlabel('PC1'); ylabel('cdf');
box off; axis square;

figure % distrib of 
histogram(abs(score(CA1_N_idx,1)),30, 'EdgeColor', mapCA1(2,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf');hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(abs(score(CA3_N_idx,1)),30, 'EdgeColor', mapCA3(2,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf'); hold on
histogram(abs(score(CA1_F_idx,1)),30, 'EdgeColor', mapCA1(1,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(abs(score(CA3_F_idx,1)),30, 'EdgeColor', mapCA3(1,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
% legend('CA1_N','CA3_N','CA1_F', 'CA3_F');
xlabel('abs(PC1)'); ylabel('cdf');
box off; axis square;

figure % distrib of 
histogram(score(CA1_N_idx,2), 'FaceColor', mapCA1(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'Normalization', 'pdf'); hold on
histogram(score(CA3_N_idx,2), 'FaceColor', mapCA3(2,:), 'EdgeColor', 'none','FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(score(CA1_F_idx,2), 'EdgeColor', mapCA1(1,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'pdf') %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(score(CA3_F_idx,2), 'EdgeColor', mapCA3(1,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'pdf') %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
legend('CA1_N','CA3_N','CA1_F', 'CA3_F');
xlabel('PC2'); ylabel('pdf');
box off; axis square;

% try different ways to fill the data in (e.g. not linear interpolation but using preceding value), or even do PCA on raw data with
% missing values (which pca can handle in different ways), to check if it
% changes things radically. 

% figure 
% gscatter(TSNEpca(:,1),TSNEpca(:,2), Kidx(:,cmax))
% xlabel('tsne1'); ylabel('tsne2');
% title('tsne on PCA scores, K-means groups')
% box off; axis square;
% 
% figure % identify 
% gscatter(TSNEpca(:,1),TSNEpca(:,2), GroupCA_VR, colors)
% legend('CA1_N','CA1_F','CA3_N', 'CA3_F');
% xlabel('tsne1'); ylabel('tsne2');
% title('tsne on PCA scores, CA_VR groups')
% box off; axis square;

figure
gscatter(TSNEraw(:,1),TSNEraw(:,2), GroupCA_VR, colors)
legend('CA1_N','CA1_F','CA3_N', 'CA3_F');
xlabel('tsne1'); ylabel('tsne2');
title('tsne on raw scores, CA_VR groups')
box off; axis square;

figure % distrib of 
histogram(TSNEraw(CA1_N_idx,1), 'FaceColor', mapCA1(2,:),'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(TSNEraw(CA1_F_idx,1), 'FaceColor', mapCA1(1,:), 'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(TSNEraw(CA3_N_idx,1), 'FaceColor', mapCA3(2,:), 'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(TSNEraw(CA3_F_idx,1), 'FaceColor', mapCA3(1,:), 'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
legend('CA1_N','CA1_F','CA3_N', 'CA3_F');
xlabel('tsne1 raw'); ylabel('pdf');
box off; axis square;

figure % same fig as above but in violin plots. Need to add bootstrap to compare means and show difference for CA3 F with other conditions
v1 = violinplot(TSNEraw(:,1),GroupCA_VR);
v1(1).ViolinColor = mapCA1(2,:);
v1(2).ViolinColor = mapCA1(1,:);
v1(3).ViolinColor = mapCA3(2,:);
v1(4).ViolinColor = mapCA3(1,:);
xlabel('CA-VR groups'); ylabel('tsne1 raw');
box off; axis square;

% figure % 
% violinplot(TSNEpca(:,1),GroupCA_VR)
% xlabel('CA-VR groups'); ylabel('tsne1 pca');
% box off; axis square;

% % plot average Traj for all PFs with tSNE(1) between 15 and 30 (difference zone for CA3 F)
% figure
% OKidx = find(TSNEraw(:,1) > 15 & TSNEraw(:,1) < 30 );
% scoreOK = score(OKidx,:);
% reconstruct = [scoreOK*coeff'];
% xlaps = repmat(1:Minlaps, size(scoreOK,1), 1);
% plot(xlaps',reconstruct', 'k-'); hold on
% plot(1:Minlaps, mean(reconstruct,1), 'r');
% yline(0,'k--');
% ylim([-300 300]);
% xlabel('lap'); ylabel('COM position (cm)');
% title('PFs trajectories for criteria')
% box off; axis square;

% plot average Traj for CA3_F PFs with tSNE(1) between 15 and 30 (high density zone for CA3 F)
% figure
% OKidx2 = find(TSNEraw(CA3_F_idx,1) < -10 );
% scoreOK = score(OKidx2,:);
% reconstruct = [scoreOK*coeff'];
% xlaps = repmat(1:Minlaps, size(scoreOK,1), 1);
% plot(xlaps',reconstruct', 'k-'); hold on
% plot(1:Minlaps, mean(reconstruct,1), 'r');
% yline(0,'k--');
% ylim([-300 300]);
% xlabel('lap'); ylabel('COM position (cm)');
% title('PFs trajectories for criteria')
% box off; axis square;
% 
% % plot average Traj for all PFs with tSNE(1) between -10 and -30 (difference zone for CA3 F)
% figure
% OKidx = find(TSNEraw(:,1) < -10 & TSNEraw(:,1) > -30 );
% scoreOK = score(OKidx,:);
% reconstruct = [scoreOK*coeff'];
% xlaps = repmat(1:Minlaps, size(scoreOK,1), 1);
% plot(xlaps',reconstruct', 'k-'); hold on
% plot(1:Minlaps, mean(reconstruct,1), 'r');
% yline(0,'k--');
% ylim([-300 300]);
% xlabel('lap'); ylabel('COM position (cm)');
% title('PFs trajectories for criteria')
% box off; axis square;

%% linear regressions

% Tok = Tall(OKtrajIdx,:); % need to exclude the PFs that were excluded for the PCA analysis to compare to PC1. 
% It also makes sense to excludes PFs that are not defined on a min number of laps when studying lap-wise dynamics. 

Fwd_idx = find(Tok.COM_slope > 0 & Tok.COM_pval<0.05);
Bwd_idx = find(Tok.COM_slope < 0 & Tok.COM_pval<0.05);
Nsig_idx = find(Tok.COM_pval>=0.05);

Tok.SuccessSig = ones(height(Tok), 1); Tok.SuccessSig(Nsig_idx) = 0;
Tok.SuccessB = zeros(height(Tok)); Tok.SuccessB(Bwd_idx) = 1;
Tok.SuccessF = zeros(height(Tok)); Tok.SuccessF(Bwd_idx) = 1;

% TdabestCA1propSig.F = table2array(Tok(CA1_F_idx,"SuccessSig"));
% TdabestCA1propSig.N = table2array(Tok(CA1_N_idx,"SuccessSig"));
% NanpadCA1 = ones(length(CA1_N_idx)-length(CA1_F_idx),1); NanpadCA1(:) = NaN;
% TdabestCA1propSig.F(end+1:end+length(CA1_N_idx)-length(CA1_F_idx),1) = NanpadCA1;
% TdabestCA1propSig = struct2table(TdabestCA1propSig);
% 
% TdabestCA1slope.N = table2array(Tok(CA1_N_idx,"COM_slope"));
% TdabestCA1slope.F = table2array(Tok(CA1_F_idx,"COM_slope"));
% TdabestCA1slope.F(end+1:end+length(CA1_N_idx)-length(CA1_F_idx),1) = NanpadCA1;
% TdabestCA1slope = struct2table(TdabestCA1slope);
% 
% TdabestCA3propSig.F = table2array(Tok(CA3_F_idx,"SuccessSig"));
% TdabestCA3propSig.N = table2array(Tok(CA3_N_idx,"SuccessSig"));
% Nanpad = ones(length(CA3_N_idx)-length(CA3_F_idx),1); Nanpad(:) = NaN;
% TdabestCA3propSig.F(end+1:end+length(CA3_N_idx)-length(CA3_F_idx),1) = Nanpad;
% TdabestCA3propSig = struct2table(TdabestCA3propSig);

% Tok.shiftDir = {};
% Tok.shiftDir{Fwd_idx} = 'F';
% Tok.shiftDir{Bwd_idx} = 'B';
% Tok.shiftDir{Nsig_idx} = 'N';

% for k = 1:OkTraj % for each PFs
%     if Tok.COM_slope(k) > 0 & Tok.COM_pval(k) < 0.05
%     Tok.shiftDir{k,:} = 'F';
%     elseif Tok.COM_slope(k) < 0 & Tok.COM_pval(k) < 0.05
%     Tok.shiftDir{k,:} = 'B';    
%     else
%     Tok.shiftDir{k,:} = 'N';    
%     end
% end

% lin reg results (like in Can's paper) for all, CA1 and CA3
figure
gscatter(Tok.COM_slope,Tok.COM_R2, Tok.shiftDir, [blue; grey; red])
xlabel('lin reg slope (cm/lap)'); ylabel('lin reg R2');
title('All conditions')
box off; axis square;

figure
CA1idx = [CA1_F_idx; CA1_N_idx ];
gscatter(Tok.COM_slope(CA1idx),Tok.COM_R2(CA1idx), Tok.shiftDir(CA1idx), [blue; grey; red], 'o', 5)
xlabel('lin reg slope (cm/lap)'); ylabel('lin reg R2');
title('CA1')
box off; axis square;

figure
gscatter(Tok.COM_slope(CA1_N_idx),Tok.COM_R2(CA1_N_idx), Tok.shiftDir(CA1_N_idx), [grey; blue; red], 'o', 5)
xlabel('lin reg slope (cm/lap)'); ylabel('lin reg R2');
title('CA1N')
box off; axis square;

figure
CA3idx = [CA3_F_idx; CA3_N_idx ];
gscatter(Tok.COM_slope(CA3idx),Tok.COM_R2(CA3idx), Tok.shiftDir(CA3idx), [grey; red; blue ], 'o', 5)
xlabel('lin reg slope (cm/lap)'); ylabel('lin reg R2');
xlim([-4 3])
title('CA3')
box off; axis square;

% histograms of CA1 and CA3 slopes in different conditions
figure
histogram(Tok.COM_slope(CA1_N_idx),30, 'EdgeColor', mapCA1(2,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf');hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(Tok.COM_slope(CA3_N_idx),30, 'EdgeColor', mapCA3(2,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf'); hold on
histogram(Tok.COM_slope(CA1_F_idx),30, 'EdgeColor', mapCA1(1,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(Tok.COM_slope(CA3_F_idx),30, 'EdgeColor', mapCA3(1,:), 'LineWidth', 2, 'DisplayStyle','stairs', 'Normalization', 'cdf'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
% legend('CA1_N','CA3_N','CA1_F', 'CA3_F');
xlabel('Slope (cm/lap)'); ylabel('cdf');
box off; axis square;

figure
subplot(1,2,1)
histogram(Tok.COM_slope(CA1_N_idx), 'EdgeColor', mapCA1(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
histogram(Tok.COM_slope(CA1_F_idx),'EdgeColor', 'none', 'FaceColor', mapCA1(1,:), 'LineWidth', 1, 'DisplayStyle','bar', 'Normalization', 'probability'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
xline(mean(Tok.COM_slope(CA1_N_idx)),'Color', mapCA1(2,:), 'LineWidth', 1.5 );hold on
xline(mean(Tok.COM_slope(CA1_F_idx)),'Color', 'g', 'LineWidth', 1.5 );hold on
xlim([-4, 3])
box off; axis square;
xlabel('Slope (cm/lap)'); ylabel('pdf');
subplot(1,2,2)
histogram(Tok.COM_slope(CA3_N_idx), 'EdgeColor', mapCA3(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability'); hold on
histogram(Tok.COM_slope(CA3_F_idx), 'EdgeColor', 'none', 'FaceColor', mapCA3(1,:), 'LineWidth', 1, 'DisplayStyle','bar', 'Normalization', 'probability'); hold on %'FaceAlpha', 0.5, 'Normalization', 'pdf'); hold on
xlim([-4, 3])
xline(mean(Tok.COM_slope(CA3_N_idx)),'Color', mapCA3(2,:), 'LineWidth', 1.5 );hold on
xline(mean(Tok.COM_slope(CA3_F_idx)),'Color', 'r', 'LineWidth', 1.5 );hold on
xlabel('Slope (cm/lap)'); ylabel('pdf');
box off; axis square;

figure
hCA1 = histogram(Tok.COM_slope(CA1idx), 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA1(2,:)); hold on
% hCA1 = histfit(Tok.COM_slope(CA1idx), 10, 'normal'); hold on
xline(mean(Tok.COM_slope(CA1idx)),'g-');hold on
% Ygauss = 0.135*exp(-0.5*([-2:0.05:2]-mean(Tok.COM_slope(CA1idx))).^2/std(Tok.COM_slope(CA1idx))^2);
% plot([-2:0.05:2], Ygauss, 'k')
xlabel('lin reg slope (cm/lap)'); ylabel('proba');
title(['CA1 F and N combined, mean = ' num2str(mean(Tok.COM_slope(CA1idx))) ' cm/lap' ])
box off; axis square;

figure
hCA1N = histogram(Tok.COM_slope(CA1_N_idx), 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA1(2,:)); hold on
xline(mean(Tok.COM_slope(CA1_N_idx)),'g-');hold on
Ygauss = 0.22*exp(-0.5*([-2:0.05:2]-mean(Tok.COM_slope(CA1_N_idx))).^2/std(Tok.COM_slope(CA1_N_idx))^2);
plot([-2:0.05:2], Ygauss, 'k')
xlabel('lin reg slope (cm/lap)'); ylabel('proba');
title(['CA1 N, mean = ' num2str(mean(Tok.COM_slope(CA1_N_idx))) ' cm/lap' ])
box off; axis square;

figure
hCA1F = histogram(Tok.COM_slope(CA1_F_idx), 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA1(2,:)); hold on
xline(mean(Tok.COM_slope(CA1_F_idx)),'g-');hold on
Ygauss = 0.18*exp(-0.5*([-2:0.05:2]-mean(Tok.COM_slope(CA1_F_idx))).^2/std(Tok.COM_slope(CA1_F_idx))^2);
plot([-2:0.05:2], Ygauss, 'k')
xlabel('lin reg slope (cm/lap)'); ylabel('proba');
title(['CA1 F, mean = ' num2str(mean(Tok.COM_slope(CA1_F_idx))) ' cm/lap' ])
box off; axis square;

figure
h = histogram(Tok.COM_slope(CA3idx), [-1.7 -1.5 -1.3 -1.1 -0.9 -0.7 -0.5 -0.3 -0.1 0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7], 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA3(2,:)); hold on
xline(mean(Tok.COM_slope(CA3idx)),'r-');hold on
Ygauss = 0.3*exp(-0.5*([-2:0.05:2]-mean(Tok.COM_slope(CA3idx))).^2/std(Tok.COM_slope(CA3idx))^2);
plot([-2:0.05:2], Ygauss, 'k')
xlabel('lin reg slope (cm/lap)'); ylabel('proba');
title('CA3')
box off; axis square;

figure
h2 = histogram(Tok.COM_slope(CA3_N_idx), [-1.7 -1.5 -1.3 -1.1 -0.9 -0.7 -0.5 -0.3 -0.1 0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7], 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA3(2,:)); hold on
xline(mean(Tok.COM_slope(CA3_N_idx)),'r-');hold on
xlabel('lin reg slope (cm/lap)'); ylabel('proba');
title('CA3_N')
box off; axis square;

figure
h3 = histogram(Tok.COM_slope(CA3_F_idx), [-1.7 -1.5 -1.3 -1.1 -0.9 -0.7 -0.5 -0.3 -0.1 0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7], 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', mapCA3(2,:)); hold on
xline(mean(Tok.COM_slope(CA3_F_idx)),'r-');hold on
xlabel('lin reg slope (cm/lap)'); ylabel('proba');
title('CA3_F')
box off; axis square;

% CA1slopes.SlopeHistFandN = hCA1; CA1slopes.FandN = Tok.COM_slope(CA1idx); CA1slopes.F = Tok.COM_slope(CA1_F_idx); CA1slopes.N = Tok.COM_slope(CA1_N_idx);
CA3slopes.SlopeHistFandN = hCA1; CA3slopes.FandN = Tok.COM_slope(CA3idx); CA3slopes.F = Tok.COM_slope(CA3_F_idx); CA3slopes.N = Tok.COM_slope(CA3_N_idx);
% save('CA3_SlopeDistrib_All', 'CA3slopes')

% CA3_FandN_SlopeHist = h; CA3_N_SlopeHist = h2; CA3_F_SlopeHist = h3;
% save('CA3_SlopeDistrib', 'CA3_FandN_SlopeHist', 'CA3_N_SlopeHist', 'CA3_F_SlopeHist')
% save('CA1_SlopeDistrib', 'CA1slopes')

% delta-delta bootstrapped stats for median abs(slope)
[slopeBackCI_CA1N, slopeBackboot_CA1N] = bootci(10000,@median, abs(Tok.COM_slope(CA1_N_idx)));
[slopeBackCI_CA1F, slopeBackboot_CA1F] = bootci(10000,@median, abs(Tok.COM_slope(CA1_F_idx)));
[slopeBackCI_CA3N, slopeBackboot_CA3N] = bootci(10000,@median, abs(Tok.COM_slope(CA3_N_idx)));
[slopeBackCI_CA3F, slopeBackboot_CA3F] = bootci(10000,@median, abs(Tok.COM_slope(CA3_F_idx)));
meanSlopeBCA1N = mean(slopeBackboot_CA1N);
meanSlopeBCA1F = mean(slopeBackboot_CA1F);
meanSlopeBCA3N = mean(slopeBackboot_CA3N);
meanSlopeBCA3F = mean(slopeBackboot_CA3F);

    % pairwise comparisons N vs F
    DeltaSlopeBCA1_NvsF = slopeBackboot_CA1F - slopeBackboot_CA1N; 
    DeltaSlopeBCA1_NvsF_CI = prctile(DeltaSlopeBCA1_NvsF,[2.5, 97.5]);
%     pvalBca1 = length(find(DeltaSlopeBCA1_NvsF>=0))/length(DeltaSlopeBCA1_NvsF);
    DeltaSlopeBCA3_NvsF = slopeBackboot_CA3F - slopeBackboot_CA3N; 
    DeltaSlopeBCA3_NvsF_CI = prctile(DeltaSlopeBCA3_NvsF,[2.5, 97.5]);
    
    % interaction between CA and N_F factors
    DeltaDeltaSlopeB = DeltaSlopeBCA1_NvsF - DeltaSlopeBCA3_NvsF;
    DeltaDeltaSlopeB_CI = prctile(DeltaDeltaSlopeB,[2.5, 97.5]);
    pvalSlopeB = length(find(DeltaDeltaSlopeB>=0))/length(DeltaDeltaSlopeB);

figure % all PFs pooled. Bar plots of props for all Ca and conditions. Delta F vs N for CA1 and CA3, oriented like in Dabest. 
subplot(1,3,1) % abs(slopes) (all PFs pooled)
    v3 = violinplot(abs(Tok.COM_slope), GroupCA_VR);
    v3(1).ViolinColor = mapCA1(2,:);
    v3(2).ViolinColor = mapCA1(1,:);
    v3(3).ViolinColor = mapCA3(2,:);
    v3(4).ViolinColor = mapCA3(1,:);
    ylabel('|slope|');
    set(gca, 'XTick', [1,2,3,4], 'XTickLabel', ['CA1N';'CA1F'; 'CA3N'; 'CA3F']);
    box off; %axis square;
subplot(1,3,2) % CA1 N vs F bootstrapped delta |slope|
    histogram(DeltaSlopeBCA1_NvsF, 'EdgeColor', mapCA1(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaSlopeBCA1_NvsF), 0.01, mean(DeltaSlopeBCA1_NvsF)-DeltaSlopeBCA1_NvsF_CI(1),mean(DeltaSlopeBCA1_NvsF)-DeltaSlopeBCA1_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
%     xlim([-meanPropBCA1F 1-meanPropBCA1F]);
    xlim([-0.3 0.2]);
    xlabel('\Delta median |slope| (F - N)'); ylabel('proba');
    title('bootstrapped distrib')
    set(gca,'view',[90 -90])
    box off; %axis square;
subplot(1,3,3) % CA3 N vs F bootstrapped delta |slope|
    histogram(DeltaSlopeBCA3_NvsF, 'EdgeColor', mapCA3(2,:), 'LineWidth', 1, 'DisplayStyle','stairs', 'Normalization', 'probability');hold on
    errorbar(mean(DeltaSlopeBCA3_NvsF), 0.01, mean(DeltaSlopeBCA3_NvsF)-DeltaSlopeBCA3_NvsF_CI(1), mean(DeltaSlopeBCA3_NvsF)-DeltaSlopeBCA3_NvsF_CI(2), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 8); hold on
    xline(0,'k--', 'LineWidth', 1.5 );hold on
    xlim([-0.3 0.2]);
    xlabel('\Delta median |slope| (F - N)'); ylabel('proba');
    title('bootstrapped distrib')
    set(gca,'view',[90 -90])

% proportions of Back and Forward shifting PFs, all pooled

% in CA1
CA1N_Fwd_idx = find(Tok.COM_slope(CA1_N_idx) > 0 & Tok.COM_pval(CA1_N_idx)<0.05);
CA1N_Bck_idx = find(Tok.COM_slope(CA1_N_idx) < 0 & Tok.COM_pval(CA1_N_idx)<0.05);
CA1N_nonsig_idx = find(Tok.COM_pval(CA1_N_idx)>=0.05);
CA1F_Fwd_idx = find(Tok.COM_slope(CA1_F_idx) > 0 & Tok.COM_pval(CA1_F_idx)<0.05);
CA1F_Bck_idx = find(Tok.COM_slope(CA1_F_idx) < 0 & Tok.COM_pval(CA1_F_idx)<0.05);
CA1F_nonsig_idx = find(Tok.COM_pval(CA1_F_idx)>=0.05);

Prop_Bck_CA1Nall = length(CA1N_Bck_idx)./length(Tok.COM_slope(CA1_N_idx))
Prop_Fwd_CA1Nall = length(CA1N_Fwd_idx)./length(Tok.COM_slope(CA1_N_idx))
Prop_NonSig_CA1Nall = length(CA1N_nonsig_idx)./length(Tok.COM_slope(CA1_N_idx))

Prop_Bck_CA1Fall = length(CA1F_Bck_idx)./length(Tok.COM_slope(CA1_F_idx))
Prop_Fwd_CA1Fall = length(CA1F_Fwd_idx)./length(Tok.COM_slope(CA1_F_idx))
Prop_NonSig_CA1Fall = length(CA1F_nonsig_idx)./length(Tok.COM_slope(CA1_F_idx))

%in CA3
CA3N_Fwd_idx = find(Tok.COM_slope(CA3_N_idx) > 0 & Tok.COM_pval(CA3_N_idx)<0.05);
CA3N_Bck_idx = find(Tok.COM_slope(CA3_N_idx) < 0 & Tok.COM_pval(CA3_N_idx)<0.05);
CA3N_nonsig_idx = find(Tok.COM_pval(CA3_N_idx)>=0.05);
CA3F_Fwd_idx = find(Tok.COM_slope(CA3_F_idx) > 0 & Tok.COM_pval(CA3_F_idx)<0.05);
CA3F_Bck_idx = find(Tok.COM_slope(CA3_F_idx) < 0 & Tok.COM_pval(CA3_F_idx)<0.05);
CA3F_nonsig_idx = find(Tok.COM_pval(CA3_F_idx)>=0.05);

Prop_Bck_CA3Nall = length(CA3N_Bck_idx)./length(Tok.COM_slope(CA3_N_idx))
Prop_Fwd_CA3Nall = length(CA3N_Fwd_idx)./length(Tok.COM_slope(CA3_N_idx))
Prop_NonSig_CA3Nall = length(CA3N_nonsig_idx)./length(Tok.COM_slope(CA3_N_idx))

Prop_Bck_CA3Fall = length(CA3F_Bck_idx)./length(Tok.COM_slope(CA3_F_idx))
Prop_Fwd_CA3Fall = length(CA3F_Fwd_idx)./length(Tok.COM_slope(CA3_F_idx))
Prop_NonSig_CA3Fall = length(CA3F_nonsig_idx)./length(Tok.COM_slope(CA3_F_idx))

% unpaired delta-delta plot comparing CA fields and familiarity levels

% propshiftCA1F = 1 - Prop_NonSig_CA1Fall;
% propshiftCA1N = 1 - Prop_NonSig_CA1Nall;
% propshiftCA3F = 1 - Prop_NonSig_CA3Fall;
% propshiftCA3N = 1 - Prop_NonSig_CA3Nall;
[propshiftCI_CA1N, propshiftboot_CA1N] = bootci(10000,@propshift, Tok.COM_pval(CA1_N_idx));
[propshiftCI_CA1F, propshiftboot_CA1F] = bootci(10000,@propshift, Tok.COM_pval(CA1_F_idx));
[propshiftCI_CA3N, propshiftboot_CA3N] = bootci(10000,@propshift, Tok.COM_pval(CA3_N_idx));
[propshiftCI_CA3F, propshiftboot_CA3F] = bootci(10000,@propshift, Tok.COM_pval(CA3_F_idx));
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
    [propshiftCI_CA1, propshiftboot_CA1] = bootci(10000,@propshift, Tok.COM_pval(CA1idx));
    [propshiftCI_CA3, propshiftboot_CA3] = bootci(10000,@propshift, Tok.COM_pval(CA3idx));
    Delta_Prop_CA1vsCA3 = propshiftboot_CA1 - propshiftboot_CA3;

    % N vs F main effect

% delta-delta bootstrapped stats for Backward shifting PFs
[propBackCI_CA1N, propBackboot_CA1N] = bootci(10000,@propBack, Tok.COM_pval(CA1_N_idx),Tok.COM_slope(CA1_N_idx));
[propBackCI_CA1F, propBackboot_CA1F] = bootci(10000,@propBack, Tok.COM_pval(CA1_F_idx),Tok.COM_slope(CA1_F_idx));
[propBackCI_CA3N, propBackboot_CA3N] = bootci(10000,@propBack, Tok.COM_pval(CA3_N_idx),Tok.COM_slope(CA3_N_idx));
[propBackCI_CA3F, propBackboot_CA3F] = bootci(10000,@propBack, Tok.COM_pval(CA3_F_idx),Tok.COM_slope(CA3_F_idx));
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
[propFwdCI_CA1N, propFwdboot_CA1N] = bootci(10000,@propFwd, Tok.COM_pval(CA1_N_idx),Tok.COM_slope(CA1_N_idx));
[propFwdCI_CA1F, propFwdboot_CA1F] = bootci(10000,@propFwd, Tok.COM_pval(CA1_F_idx),Tok.COM_slope(CA1_F_idx));
[propFwdCI_CA3N, propFwdboot_CA3N] = bootci(10000,@propFwd, Tok.COM_pval(CA3_N_idx),Tok.COM_slope(CA3_N_idx));
[propFwdCI_CA3F, propFwdboot_CA3F] = bootci(10000,@propFwd, Tok.COM_pval(CA3_F_idx),Tok.COM_slope(CA3_F_idx));
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


% resampling stats to check whether props and mean slope are different 
% btw CA1 vs CA3 regardless of sample size
nboot = 1000;
replacement = true(1);
for n = 1:nboot
    Res = randsample(CA1idx, length(CA3idx), replacement);
    ResMslope(n) = mean(Tok.COM_slope(Res));
    ResMAslope(n) = median(abs(Tok.COM_slope(Res)));

    IdxFwd_CA1res = find(Tok.COM_slope(Res)>0 & Tok.COM_pval(Res)<=0.05);
    IdxBck_CA1res = find(Tok.COM_slope(Res) < 0 & Tok.COM_pval(Res)<=0.05);
    IdxNonSig_CA1res = find(Tok.COM_pval(Res)>0.05);
    length_CA1res = length(IdxFwd_CA1res) + length(IdxBck_CA1res) + length(IdxNonSig_CA1res);
    PropFwd_CA1res(n) = length(IdxFwd_CA1res)/length_CA1res;
    PropBck_CA1res(n) = length(IdxBck_CA1res)/length_CA1res;
    PropNonSig_CA1res(n) = length(IdxNonSig_CA1res)/length_CA1res;
    PropShift_CA1res(n) = PropFwd_CA1res(n) + PropBck_CA1res(n);
end
    IdxFwd_CA3res = find(Tok.COM_slope(CA3idx)>0 & Tok.COM_pval(CA3idx)<=0.05);
    IdxBck_CA3res = find(Tok.COM_slope(CA3idx) < 0 & Tok.COM_pval(CA3idx)<=0.05);
    IdxNonSig_CA3res = find(Tok.COM_pval(CA3idx)>0.05);
    length_CA3res = length(IdxFwd_CA3res) + length(IdxBck_CA3res) + length(IdxNonSig_CA3res);
    PropFwd_CA3res = length(IdxFwd_CA3res)/length_CA3res;
    PropBck_CA3res = length(IdxBck_CA3res)/length_CA3res;
    PropNonSig_CA3res = length(IdxNonSig_CA3res)/length_CA3res;
    PropShift_CA3res = PropFwd_CA3res + PropBck_CA3res;

figure
subplot(2,2,1)
    histogram(ResMslope, 'EdgeColor', 'none', 'FaceColor', mapCA1(2,:), 'Normalization', 'probability'); hold on
    xline(mean(ResMslope),'g-', 'LineWidth', 2);hold on
    xline(mean(Tok.COM_slope(CA3idx)),'-', 'Color', mapCA3(2,:), 'LineWidth', 2);hold on
    xline(0, 'k--');
    xlim([-0.4 0.1])
    xlabel('mean slope (cm/lap)')
    ylabel('probability')
%     title('resampled CA1 vs CA3')
    axis square; box off;
subplot(2,2,2)
    histogram(ResMAslope, 'EdgeColor', 'none', 'FaceColor', mapCA1(2,:), 'Normalization', 'probability'); hold on
    xline(mean(ResMAslope),'g-', 'LineWidth', 2);hold on
    xline(median(abs(Tok.COM_slope(CA3idx))),'-', 'Color', mapCA3(2,:), 'LineWidth', 2);hold on
    % xline(0, 'k-');
    % xlim([-0.4 0.1])
    xlabel('median absolute slope (cm/lap)')
    ylabel('probability')
%     title('resampled CA1 vs CA3')
    axis square; box off;
subplot(2,2,3)
    histogram(PropShift_CA1res, 'EdgeColor', 'none', 'FaceColor', mapCA1(2,:), 'Normalization', 'probability'); hold on
    xline(mean(PropShift_CA1res),'g-', 'LineWidth', 2); hold on
    xline(PropShift_CA3res,'-', 'Color', mapCA3(2,:), 'LineWidth', 2); 
    xlabel('proportion of lin shifting PFs')
    ylabel('probability')
%     title('resampled CA1 vs CA3')
    axis square; box off;
subplot(2,2,4)
    histogram(PropFwd_CA1res, 'EdgeColor', 'none', 'FaceColor', mapCA1(2,:), 'Normalization', 'probability'); hold on
    xline(mean(PropFwd_CA1res),'g-', 'LineWidth', 2); hold on
    xline(PropFwd_CA3res,'-', 'Color', mapCA3(2,:), 'LineWidth', 2); 
    xlabel('proportion of forward shifting PFs')
    ylabel('probability')
%     title('resampled CA1 vs CA3')
    axis square; box off;

% lin reg proportions per animal and conditions
CA1mice = unique(Tok.animalID(CA1_N_idx));
CA3mice = unique(Tok.animalID(CA3_N_idx));

for m = 1:length(CA1mice)
    IdxFwd_F{m} = find(Tok.COM_slope>0 & Tok.COM_pval<=0.05 & ismember(Tok.animalID, CA1mice{m}) & ismember(Tok.f_n, 'f'));
    IdxBck_F{m} = find(Tok.COM_slope < 0 & Tok.COM_pval<=0.05 & ismember(Tok.animalID, CA1mice{m}) & ismember(Tok.f_n, 'f'));
    IdxNonSig_F{m} = find(Tok.COM_pval>0.05 & ismember(Tok.animalID, CA1mice{m}) & ismember(Tok.f_n, 'f'));
    IdxCA1F{m} = [IdxFwd_F{m}; IdxBck_F{m}; IdxNonSig_F{m}];
    length_CA1F(m) = length(IdxFwd_F{m}) + length(IdxBck_F{m}) + length(IdxNonSig_F{m});
    length_CA1F2(m) = length(find(ismember(Tok.animalID, CA1mice{1}) & ismember(Tok.f_n, 'f')));

    IdxFwd_N{m} = find(Tok.COM_slope>0 & Tok.COM_pval<=0.05 & ismember(Tok.animalID, CA1mice{m}) & ismember(Tok.f_n, 'n'));
    IdxBck_N{m} = find(Tok.COM_slope < 0 & Tok.COM_pval<=0.05 & ismember(Tok.animalID, CA1mice{m}) & ismember(Tok.f_n, 'n'));
    IdxNonSig_N{m} = find(Tok.COM_pval>0.05 & ismember(Tok.animalID, CA1mice{m}) & ismember(Tok.f_n, 'n'));
    IdxCA1N{m} = [IdxFwd_N{m}; IdxBck_N{m}; IdxNonSig_N{m}];
    length_CA1N(m) = length(IdxFwd_N{m}) + length(IdxBck_N{m}) + length(IdxNonSig_N{m});
    length_CA1N2(m) = length(find(ismember(Tok.animalID, CA1mice{m}) & ismember(Tok.f_n, 'n')));
    
    PropFwd_CA1F(m) = length(IdxFwd_F{m})/length_CA1F(m);
    PropBck_CA1F(m) = length(IdxBck_F{m})/length_CA1F(m);
    PropNonSig_CA1F(m) = length(IdxNonSig_F{m})/length_CA1F(m);
    PropShift_CA1F(m) = PropFwd_CA1F(m) + PropBck_CA1F(m);
    PropTot_CA1F(m) = PropFwd_CA1F(m) + PropBck_CA1F(m) + PropNonSig_CA1F(m);

    PropFwd_CA1N(m) = length(IdxFwd_N{m})/length_CA1N(m);
    PropBck_CA1N(m) = length(IdxBck_N{m})/length_CA1N(m);
    PropNonSig_CA1N(m) = length(IdxNonSig_N{m})/length_CA1N(m);
    PropShift_CA1N(m) = PropFwd_CA1N(m) + PropBck_CA1N(m);
    PropTot_CA1N(m) = PropFwd_CA1N(m) + PropBck_CA1N(m) + PropNonSig_CA1N(m);

    meanSlopeCA1F(m) = mean(Tok.COM_slope(IdxCA1F{m}));
    meanSlopeCA1N(m) = mean(Tok.COM_slope(IdxCA1N{m}));
    medianASlopeCA1F(m) = median(abs(Tok.COM_slope(IdxCA1F{m})));
    medianASlopeCA1N(m) = median(abs(Tok.COM_slope(IdxCA1N{m})));
end
PropFwd_CA1 = [PropFwd_CA1N;PropFwd_CA1F];
PropBck_CA1 = [PropBck_CA1N;PropBck_CA1F];
PropNonSig_CA1 = [PropNonSig_CA1N;PropNonSig_CA1F];

[H_fwdCA1,P_fwdCA1,CI_fwdCA1,STATS_fwdCA1] = ttest(diff(PropFwd_CA1)); % p = 0.057
[H_bckCA1,P_bckCA1,CI_bckCA1,STATS_bckCA1] = ttest(diff(PropBck_CA1)); % p = 0.0203
[H_nsCA1,P_nsCA1,CI_nsCA1,STATS_nsCA1] = ttest(diff(PropNonSig_CA1)); % p = 0.0346

figure
plot(PropFwd_CA1,'-o', 'Color', red); hold on
plot(PropBck_CA1,'-o', 'Color', blue); hold on
plot(PropNonSig_CA1,'-o', 'Color', grey); hold on
plot(mean(PropFwd_CA1,2),'-', 'Color', red, 'LineWidth', 2); hold on
plot(mean(PropBck_CA1,2),'-', 'Color', blue, 'LineWidth', 2); hold on
plot(mean(PropNonSig_CA1,2),'-', 'Color', grey, 'LineWidth', 2); hold on
ylim([0 1]); xlim([0 3]);
ylabel('Prop of PFs'); 
set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
title('CA1');
box off; axis square

clear IdxFwd_F IdxBck_F IdxNonSig_F IdxFwd_N IdxBck_N IdxNonSig_N

for m = 1:length(CA3mice)
    IdxFwd_F{m} = find(Tok.COM_slope>0 & Tok.COM_pval<=0.05 & ismember(Tok.animalID, CA3mice{m}) & ismember(Tok.f_n, 'f'));
    IdxBck_F{m} = find(Tok.COM_slope < 0 & Tok.COM_pval<=0.05 & ismember(Tok.animalID, CA3mice{m}) & ismember(Tok.f_n, 'f'));
    IdxNonSig_F{m} = find(Tok.COM_pval>0.05 & ismember(Tok.animalID, CA3mice{m}) & ismember(Tok.f_n, 'f'));
    IdxCA3F{m} = [IdxFwd_F{m}; IdxBck_F{m}; IdxNonSig_F{m}];
    length_CA3F(m) = length(IdxFwd_F{m}) + length(IdxBck_F{m}) + length(IdxNonSig_F{m});
    length_CA3F2(m) = length(find(ismember(Tok.animalID, CA3mice{1}) & ismember(Tok.f_n, 'f')));

    IdxFwd_N{m} = find(Tok.COM_slope>0 & Tok.COM_pval<=0.05 & ismember(Tok.animalID, CA3mice{m}) & ismember(Tok.f_n, 'n'));
    IdxBck_N{m} = find(Tok.COM_slope < 0 & Tok.COM_pval<=0.05 & ismember(Tok.animalID, CA3mice{m}) & ismember(Tok.f_n, 'n'));
    IdxNonSig_N{m} = find(Tok.COM_pval>0.05 & ismember(Tok.animalID, CA3mice{m}) & ismember(Tok.f_n, 'n'));
    IdxCA3N{m} = [IdxFwd_N{m}; IdxBck_N{m}; IdxNonSig_N{m}];
    length_CA3N(m) = length(IdxFwd_N{m}) + length(IdxBck_N{m}) + length(IdxNonSig_N{m});
    length_CA3N2(m) = length(find(ismember(Tok.animalID, CA3mice{m}) & ismember(Tok.f_n, 'n')));

    PropFwd_CA3F(m) = length(IdxFwd_F{m})/length_CA3F(m);
    PropBck_CA3F(m) = length(IdxBck_F{m})/length_CA3F(m);
    PropNonSig_CA3F(m) = length(IdxNonSig_F{m})/length_CA3F(m);
    PropShift_CA3F(m) = PropFwd_CA3F(m) + PropBck_CA3F(m);
    PropTot_CA3F(m) = PropFwd_CA3F(m) + PropBck_CA3F(m) + PropNonSig_CA3F(m);

    PropFwd_CA3N(m) = length(IdxFwd_N{m})/length_CA3N(m);
    PropBck_CA3N(m) = length(IdxBck_N{m})/length_CA3N(m);
    PropNonSig_CA3N(m) = length(IdxNonSig_N{m})/length_CA3N(m);
    PropShift_CA3N(m) = PropFwd_CA3N(m) + PropBck_CA3N(m);
    PropTot_CA3N(m) = PropFwd_CA3N(m) + PropBck_CA3N(m) + PropNonSig_CA3N(m);

    meanSlopeCA3F(m) = mean(Tok.COM_slope(IdxCA3F{m}));
    meanSlopeCA3N(m) = mean(Tok.COM_slope(IdxCA3N{m}));
    medianASlopeCA3F(m) = median(abs(Tok.COM_slope(IdxCA3F{m})));
    medianASlopeCA3N(m) = median(abs(Tok.COM_slope(IdxCA3N{m})));
end
PropFwd_CA3 = [PropFwd_CA3N;PropFwd_CA3F];
PropBck_CA3 = [PropBck_CA3N;PropBck_CA3F];
PropNonSig_CA3 = [PropNonSig_CA3N;PropNonSig_CA3F];

[H_fwdCA3,P_fwdCA3,CI_fwdCA3,STATS_fwdCA3] = ttest(diff(PropFwd_CA3)) % p = 0.0514
[H_bckCA3,P_bckCA3,CI_bckCA3,STATS_bckCA3] = ttest(diff(PropBck_CA3)) % p = 0.0076
[H_nsCA3,P_nsCA3,CI_nsCA3,STATS_nsCA3] = ttest(diff(PropNonSig_CA3)) % p = 0.0528

figure
plot(PropFwd_CA3,'-o', 'Color', red); hold on
plot(PropBck_CA3,'-o', 'Color', blue); hold on
plot(PropNonSig_CA3,'-o', 'Color', grey); hold on
plot(mean(PropFwd_CA3,2),'-', 'Color', red, 'LineWidth', 2); hold on
plot(mean(PropBck_CA3,2),'-', 'LineWidth', 2, 'Color', blue); hold on
plot(mean(PropNonSig_CA3,2),'-',  'LineWidth', 2, 'Color', grey); hold on
ylim([0 1]); xlim([0 3]);
ylabel('Prop of PFs'); 
set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
title('CA3');
box off; axis square

%stacked bars
figure
% proportion of forward, backward and stable/non-signif PFs (mean across animals)
b = bar([mean(PropBck_CA1N) mean(PropFwd_CA1N) mean(PropNonSig_CA1N); mean(PropBck_CA1F) mean(PropFwd_CA1F) mean(PropNonSig_CA1F); mean(PropBck_CA3N) mean(PropFwd_CA3N) mean(PropNonSig_CA3N); mean(PropBck_CA3F) mean(PropFwd_CA3F) mean(PropNonSig_CA3F)], 'stacked', 'FaceColor','flat');
b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
ylabel('proportion of PFs');
set(gca, 'XTick', [1, 2, 3, 4], 'XTickLabel', ['CA1-N';'CA1-F'; 'CA3-N'; 'CA3-F']);
box off; axis square;

% CA1 vs CA3 prop of shifting cells / mean Slope; 
PropShift_CA1 = [PropShift_CA1N;PropShift_CA1F];
PropShift_CA3 = [PropShift_CA3N;PropShift_CA3F];
meanSlope_CA1 = [meanSlopeCA1N;meanSlopeCA1F];
meanSlope_CA3 = [meanSlopeCA3N;meanSlopeCA3F];
medianASlope_CA1 = [medianASlopeCA1N;medianASlopeCA1F];
medianASlope_CA3 = [medianASlopeCA3N;medianASlopeCA3F];

wmeanSlopeCA1N = length_CA1N.*meanSlopeCA1N./sum(length_CA1N);
wmeanSlopeCA1F = length_CA1F.*meanSlopeCA1F./sum(length_CA1F);
wmeanSlopeCA3N = length_CA3N.*meanSlopeCA3N./sum(length_CA3N);
wmeanSlopeCA3F = length_CA3F.*meanSlopeCA3F./sum(length_CA3F);
wmeanSlope_CA1 = [length_CA1N.*meanSlopeCA1N./size(Tok,1);length_CA1F.*meanSlopeCA1F./size(Tok,1)];
wmeanSlope_CA3 = [length_CA3N.*meanSlopeCA3N./size(Tok,1);length_CA3F.*meanSlopeCA3F./size(Tok,1)];

[H_psCA3,P_psCA3,CI_psCA3,STATS_psCA3] = ttest(diff(PropShift_CA3)) % p = 0.0528
[H_psCA1,P_psCA1,CI_psCA1,STATS_psCA1] = ttest(diff(PropShift_CA1)) % p = 0.0346
[H_CA1vs3N,P_CA1vs3N,CI_CA1vs3N,STATS_CA1vs3N] = ttest2(PropShift_CA1N, PropShift_CA3N) % p = 0.0568
[H_CA1vs3F,P_CA1vs3F,CI_CA1vs3F,STATS_CA1vs3F] = ttest2(PropShift_CA1F, PropShift_CA3F) % p = 0.1595

[H_msCA3,P_msCA3,CI_msCA3,STATS_msCA3] = ttest(diff(meanSlope_CA3)) % p = 0.24
[H_msCA1,P_msCA1,CI_msCA1,STATS_msCA1] = ttest(diff(meanSlope_CA1)) % p = 0.0219
[H_msCA1vs3N,P_msCA1vs3N,CI_msCA1vs3N,STATS_msCA1vs3N] = ttest2(meanSlopeCA1N, meanSlopeCA3N) % p =0.09
[H_msCA1vs3F,P_msCA1vs3F,CI_msCA1vs3F,STATS_msCA1vs3F] = ttest2(meanSlopeCA1F, meanSlopeCA3F) % p=0.2 

[H_wmsCA3,P_wmsCA3,CI_wmsCA3,STATS_wmsCA3] = ttest(diff(wmeanSlope_CA3)) % p = 0.1240
[H_wmsCA1,P_wmsCA1,CI_wmsCA1,STATS_wmsCA1] = ttest(diff(wmeanSlope_CA1)) % p = 0.3214
[H_wmsCA1vs3N,P_wmsCA1vs3N,CI_wmsCA1vs3N,STATS_wmsCA1vs3N] = ttest2(wmeanSlopeCA1N, wmeanSlopeCA3N) % p =0.2390
[H_wmsCA1vs3F,P_wmsCA1vs3F,CI_wmsCA1vs3F,STATS_wmsCA1vs3F] = ttest2(wmeanSlopeCA1F, wmeanSlopeCA3F) % p=0.1377 

[Hslopes,Pslopes(1),CIslopes,STATSslopes] = ttest(medianASlopeCA1N, medianASlopeCA1F);
[HslopesCA3,Pslopes(2),CIslopesCA3,STATSslopesCA3] = ttest(medianASlopeCA3N, medianASlopeCA3F);
[HslopesCA1vs3N,Pslopes(3),CIslopesCA1vs3N,STATSslopesCA1vs3N] = ttest2(medianASlopeCA1N, medianASlopeCA3N);
[HslopesCA1vs3F,Pslopes(4),CIslopesCA1vs3F,STATSslopesCA1vs3F] = ttest2(medianASlopeCA1F, medianASlopeCA3F);
pBonf = Pslopes*4; %highly conservative bonferroni correction (multiplication with number of comparisons being implicitely made) 

%linear mixed effect model 
mice = [CA1mice;CA1mice;CA3mice;CA3mice];
F_N = [repmat({'N'}, length(CA1mice), 1);repmat({'F'}, length(CA1mice), 1); repmat({'N'}, length(CA3mice), 1); repmat({'F'}, length(CA3mice), 1) ];
CA = [repmat({'CA1'}, 2*length(CA1mice), 1);repmat({'CA3'}, 2*length(CA3mice), 1)];
PropShift = [PropShift_CA1N'; PropShift_CA1F'; PropShift_CA3N'; PropShift_CA3F'];
% Slopes = [meanSlopeCA1N, meanSlopeCA1F, meanSlopeCA3N, meanSlopeCA3F]';
% Slopes = [wmeanSlopeCA1N, wmeanSlopeCA1F, wmeanSlopeCA3N, wmeanSlopeCA3F]';
Slopes = [medianASlopeCA1N, medianASlopeCA1F, medianASlopeCA3N, medianASlopeCA3F]';

tbl1 = table(mice, CA, F_N, PropShift, Slopes, 'VariableNames',{'mice', 'CA', 'F_N', 'PropShift', 'Slopes'});
tbl1.CA = categorical(tbl1.CA);
tbl1.mice = categorical(tbl1.mice);
tbl1.F_N = categorical(tbl1.F_N);

lme = fitlme(tbl1,'PropShift ~ 1 + CA + F_N + (1 + F_N | mice)'); % interaction CA * F_N excluded because not signif. 
anova(lme) % Fixed factors CA signif (p = 0.02889) and F_N signif (p = 0.00061)

lme_slope = fitlme(tbl1,'Slopes ~ 1 + CA * F_N + (1 + F_N | mice)'); % interaction CA * F_N excluded because not signif. 
anova(lme_slope) % 
% CA signif for abs(mean(slope)) (0.00036) but not F_N (0.9). For slope, Fixed factors CA signif (p = 0.0071) and F_N almost signif (p = 0.051) when number of Pfs per animal not taken into account. When weighted by num of PFs per animal and conditions/total recorded PFs across animals and conditions, CA not signif (0.0886) and F_N signif (0.0438)

% dabest.ID = [CA1mice; CA1mice ; CA3mice; CA3mice];
% dabest.CA = CA; %[repmat('CA1',2*length(CA1mice),1); repmat('CA3',2*length(CA3mice),1)];
% dabest.NF = [repmat('N',length(CA1mice),1); repmat('F',length(CA1mice),1); repmat('N',length(CA3mice),1); repmat('F',length(CA3mice),1)];
% dabest.propB = [PropBck_CA1N'; PropBck_CA1F'; PropBck_CA3N'; PropBck_CA3F'];
% dabest.propF = [PropFwd_CA1N'; PropFwd_CA1F'; PropFwd_CA3N'; PropFwd_CA3F'];
% dabest.propshift = PropShift;
% dabest.absSlopes = Slopes;
% dabest = struct2table(dabest);
% writetable(dabest,'dabest.csv')

figure % stats per animals and conditions
subplot(1,2,1) %mean abs slope per animal and conditions
    plot([0.8 2.2], mean(medianASlope_CA1,2),'ok-', 'LineWidth', 0.75, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', mapCA1(2,:), 'MarkerSize', 8); hold on
    plot([0.8 2.2],mean(medianASlope_CA3,2),'dk-', 'LineWidth', 0.75, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', mapCA3(2,:), 'MarkerSize', 8); hold on
    plot(medianASlope_CA1,'-', 'Color', mapCA1(2,:), 'LineWidth', 1); hold on
    plot(medianASlope_CA3,'-', 'Color', mapCA3(2,:), 'LineWidth', 1); hold on
    xlim([0 3]);
    ylim([0 0.7]);
    % ylim([-0.5 0.5]);
    ylabel('median |slope| (cm/lap)'); 
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
    box off; axis square
subplot(1,2,2) %proportions of shifting cells per animal and conditions
    plot([0.8 2.2], mean(PropShift_CA1,2),'ok-', 'LineWidth', 0.75, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', mapCA1(2,:), 'MarkerSize', 8); hold on
    plot([0.8 2.2], mean(PropShift_CA3,2),'dk-', 'LineWidth', 0.75, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', mapCA3(2,:), 'MarkerSize', 8); hold on
    plot(PropShift_CA1,'-', 'Color', mapCA1(2,:), 'LineWidth', 1); hold on
    plot(PropShift_CA3,'-', 'Color', mapCA3(2,:), 'LineWidth', 1); hold on
    ylim([0 1]); xlim([0 3]);
    ylabel('Prop of linearly shifting PFs'); 
    set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
    box off; axis square

% save('CA1perAnimal_slopes_props', 'meanSlope_CA1', 'meanASlope_CA1', 'PropShift_CA1', 'PropFwd_CA1', 'PropBck_CA1', 'PropNonSig_CA1');

% figure % weighted (by number of PFs) mean slope per animal and conditions
% plot(mean(wmeanSlope_CA1,2),'-', 'LineWidth', 2, 'Color', mapCA1(2,:)); hold on
% plot(mean(wmeanSlope_CA3,2),'-', 'LineWidth', 2, 'Color', mapCA3(2,:)); hold on
% legend('CA1', 'CA3')
% plot(wmeanSlope_CA1,'o-', 'Color', mapCA1(2,:)); hold on
% plot(wmeanSlope_CA3,'d-', 'Color', mapCA3(2,:)); hold on
% xlim([0 3]);
% % ylim([0 1]);
% ylabel('weighted mean slope (cm/lap)'); 
% set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
% box off; axis square

%% are the PFs with good linear fits and the PFs with large PC1 score the same? 
% Or some are linear and some just have a strong early shift?

% PC1 score vs linear fit
figure
gscatter(score(:,1),Tok.COM_R2, Tok.shiftDir, [blue; grey; red], 'o', 5)
xlabel('PC1 score'); ylabel('lin reg R2');
box off; axis square;

figure
scatter3(Tok.COM_slope(Fwd_idx),Tok.COM_R2(Fwd_idx), score(Fwd_idx,1), 2, 'r'); hold on
scatter3(Tok.COM_slope(Bwd_idx),Tok.COM_R2(Bwd_idx), score(Bwd_idx,1), 2, 'b'); hold on
scatter3(Tok.COM_slope(Nsig_idx),Tok.COM_R2(Nsig_idx), score(Nsig_idx,1), 2, 'k'); hold on
xlabel('lin reg slope (cm/lap)'); ylabel('lin reg R2'); zlabel('PC1 score');
box off; axis square;

figure
scatter(Tok.COM_slope(Fwd_idx), score(Fwd_idx,1), 2, red); hold on
scatter(Tok.COM_slope(Bwd_idx), score(Bwd_idx,1), 2, blue); hold on
scatter(Tok.COM_slope(Nsig_idx), score(Nsig_idx,1), 2, grey); hold on
xlabel('lin reg slope (cm/lap)'); ylabel('PC1 score');
box off; axis square;

%% linear Rsquare vs nonlin RC-like Rsquare
RCregR2adjAll = Tok.RCreg_R2adj;
RCregR2adjAll(RCregR2adjAll<0) = 0;
RCregR2 = Tok.RCreg_R2;
RCregR2(RCregR2<0) = 0;

figure
gscatter(Tok.RCreg_R2,Tok.COM_R2, Tok.shiftDir, [blue; grey; red], 'o', 5)
xlabel('nonlin RC-like R2'); ylabel('lin reg R2');
xlim([0 1])
box off; axis square;

figure
g=gramm('x',RCregR2,'y',Tok.COM_R2,'color',Tok.shiftDir);
g.set_color_options('map', [blue; red; grey])
g.geom_point() 
g.stat_cornerhist()
g.set_names('x','nl R2','y','linear R2','color','direction (lin reg)')
g.draw();
plot(g.facet_axes_handles, [0 max(RCregR2)], [0 max(RCregR2)], 'k--', 'LineWidth', 1.5);
set([g.results.geom_point_handle],'MarkerSize',3);
% set(g.facet_axes_handles,'XLim', [-0.05 1.5], 'YLim', [-0.05 1.5])

figure
[~, pval] = ttest(RCregR2-Tok.COM_R2)
[pu, hu, statsu] = signrank(RCregR2-Tok.COM_R2)
scatterCornerHist(RCregR2, Tok.COM_R2); hold on
xlabel('nonlin RC-like R2'); ylabel('lin reg R2');
title({['ttest p = ' num2str(round(pval, 2, 'significant')) ', mean = ' num2str(mean(round(RCregR2 - Tok.COM_R2, 2, 'significant')))]...
    ; ['utest p = ' num2str(round(pu, 2, 'significant')), ' median = ' num2str(median(round(RCregR2 - Tok.COM_R2, 2, 'significant')))]})

% figure % same as below but highlighting the nonlinear shift direction
% colormap(brewermap(2,'*RdBu'));
% scatter(RCregR2,Tok.COM_R2, 20, Tok.RCreg_p1)
% xlabel('nonlin R2'); ylabel('lin reg R2');
% xlim([0 1])
% % caxis([-200 200])
% box off; axis square;

figure % same as below but highlighting the nonlinear shift direction
AmpMap = colormap(brewermap(256,'*RdBu'));
scatter(RCregR2,Tok.COM_R2, 10, Tok.RCreg_p1, 'filled', 'MarkerFaceAlpha', 0.7); hold on
plot([0 1], [0 1], '--', 'LineWidth', 1.5, 'Color', [0.8 0.8 0.8]);
xlabel('nonlin R2'); ylabel('lin reg R2');
% xlim([-0.05 1])
% ylim([-0.05 1])
xlim([-0.0719 1.5106])
ylim([-0.0628 1.3182])
set(gca, 'Position', [0.1110 0.1384 0.5731 0.8216])
caxis([-50 50])
box off; 
axis square;
clb2 = colorbar; clb2.Label.String = 'nl Amp (cm)'; 

% figure
% g=gramm('x',RCregR2,'y',Tok.COM_R2,'color',Tok.RCreg_p1);
% % g.set_color_options('map', AmpMap)
% g.set_continuous_color('colormap', 'cool', 'CLim', [-50 50])
% g.geom_point() 
% g.stat_cornerhist()
% g.set_names('x','nl R2','y','linear R2','color','nl Amp (cm)')
% g.draw();
% plot(g.facet_axes_handles, [0 max(RCregR2)], [0 max(RCregR2)], 'k--', 'LineWidth', 1.5);
% set([g.results.geom_point_handle],'MarkerSize',3);

% figure
% colormap(brewermap(2,'*RdBu'));
% scatter(Tok.RCreg_p1,Tok.RCreg_p2, 20, Tok.RCreg_p1)
% xlabel('nl Amp (cm)'); ylabel('nl Tau (laps)');
% % caxis([-200 200])
% box off; axis square;

%focus on nonsig linear shifts and plot Amp vs Tau vs R2 (color coded) 
figure
colormap(flipud(copper))
% colormap(flipud(parula))
[RCregR2Sorted,nlR2SortIdx] = sort(RCregR2(Nsig_idx));
nlAmpNS = Tok.RCreg_p1(Nsig_idx); nlTauNS = Tok.RCreg_p2(Nsig_idx);
scatter(nlAmpNS(nlR2SortIdx), nlTauNS(nlR2SortIdx), 30, RCregR2Sorted, 'filled', 'MarkerFaceAlpha', 0.8)
xlabel('nl Amp (cm)')
ylabel('nl Tau (laps)')
clb = colorbar;
clb.Label.String = 'nl R-square'; 
% clb.Limits = [0 max(RCregR2adjAll)];
box off; axis square
title('nl params and fits of non linearly shifting PFs (p>0.05)')

% distribution of nonlin RCreg parameters (Amplitude and time constant)
figure
subplot(1,2,1)
histogram(Tok.RCreg_p1, 'Normalization', 'probability')
xlabel('RCreg amplitude, in cm'); ylabel('PF frequency')
box off; axis square;
subplot(1,2,2)
histogram(Tok.RCreg_p2, 40, 'Normalization', 'probability')
xlabel('RCreg time constant, in laps'); ylabel('PF frequency')
box off; axis square;

figure % nonlin RCreg params as a function of Rsquare
[RampR2,PampR2]=corrcoef(abs(Tok.RCreg_p1), Tok.RCreg_R2);
subplot(1,2,1)
gscatter(abs(Tok.RCreg_p1),Tok.RCreg_R2, Tok.shiftDir, [blue; grey; red])
ylim([0 1])
xlabel('absolute Amp (cm)'); ylabel('RCreg R2');
axis square
% box off; axis square;
subplot(1,2,2)
gscatter(Tok.RCreg_p2,Tok.RCreg_R2, Tok.shiftDir, [blue; grey; red])
ylim([0 1])
xlabel('Tau (laps)'); ylabel('RCreg R2');
axis square
% box off; axis square;

% Are abs(Amp) and Tau correlated? Can they be explained by a single parameter like p(CS)?
[b_AmpTau,~,~,~,stats_AmpTau] = regress(Tok.RCreg_p2, [ones(size(abs(Tok.RCreg_p1))), abs(Tok.RCreg_p1)]);
[Ramptau,Pamptau]=corrcoef(abs(Tok.RCreg_p1), Tok.RCreg_p2);
fAT = figure
subplot(2,1,1)
h = histogram2(abs(Tok.RCreg_p1),Tok.RCreg_p2,[10 10],'DisplayStyle','tile','ShowEmptyBins','off', 'Normalization', 'count');
% xlabel('absolute Amp (cm)'); ylabel('Tau (laps)');
box off; axis square;
cb = colorbar;
cb.Label.String = 'number of PFs'; cb.Limits = [0 max(h.Values, [], 'all')];
subplot(2,1,2)
[RCregR2SortedAll,nlR2SortIdxAll] = sort(RCregR2);
scatter(abs(Tok.RCreg_p1), Tok.RCreg_p2, 20, RCregR2SortedAll, 'filled', 'MarkerFaceAlpha', 0.8); hold on
colormap(gca, flipud(copper))
% xlabel('absolute Amp (cm)'); ylabel('Tau (laps)');
clb = colorbar;
clb.Label.String = 'nl R-square'; 
axis square
box off
% fAT.Units = 'Normalized'
[ax1, h1] = suplabel('absolute Amp (cm)')
[ax2, h2] = suplabel('Tau (laps)','y')
h1.Position = [0.5 0 0];
h2.Position = [0.25 0.5 0];
% fAT.Renderer = 'painters';
% print(fAT, '-vector','-dpdf','Amp_vs_Tau_DistribandCorr.pdf')

figure % absolute Amp and tau as a function of each other, for different conditions. 4 panels. 
subplot(2,2,1)
gscatter(abs(Tok.RCreg_p1(CA1_N_idx)),Tok.RCreg_p2(CA1_N_idx), Tok.shiftDir(CA1_N_idx), [grey; blue; red], '.', 10, 'off')
% xlabel('RCreg amplitude (cm)'); ylabel('RCreg tau (laps)');
xlim([0 210])
title('CA1N')
axis square, box off
subplot(2,2,2)
gscatter(abs(Tok.RCreg_p1(CA1_F_idx)),Tok.RCreg_p2(CA1_F_idx), Tok.shiftDir(CA1_F_idx), [blue; grey; red], '.', 10, 'off')
% xlabel('RCreg amplitude (cm)'); ylabel('RCreg tau (laps)');
xlim([0 210])
title('CA1F')
axis square, box off
subplot(2,2,3)
gscatter(abs(Tok.RCreg_p1(CA3_N_idx)),Tok.RCreg_p2(CA3_N_idx), Tok.shiftDir(CA3_N_idx), [blue; grey; red], '.', 10, 'off')
% xlabel('RCreg amplitude (cm)'); ylabel('RCreg tau (laps)');
xlim([0 210])
title('CA3N')
axis square, box off
subplot(2,2,4)
gscatter(abs(Tok.RCreg_p1(CA3_F_idx)),Tok.RCreg_p2(CA3_F_idx), Tok.shiftDir(CA3_F_idx), [grey; red; blue], '.', 10, 'off')
% xlabel('RCreg amplitude (cm)'); ylabel('RCreg tau (laps)');
xlim([0 210])
title('CA3F')
axis square, box off
suplabel('|Amp| (cm)');
suplabel('Tau (laps)','y');

figure % what different Taus correspond to
tau_range = [0.1, 1, 5, 10, 20, 50, 100];
Xeval2 = [1:0.1:Minlaps]-1;
for n = 1: length(tau_range)
    curv_ex{n} = 1*(1-exp(-Xeval2/tau_range(n)));
    plot(Xeval2, curv_ex{n}, 'Color', cline(n,:) ); hold on 
end
legend('1', '3', '5', '10', '20', '50', '100')
xlabel('laps from PF onset')
ylabel('COM shift (cm)')
axis square; box off;
% large tau = flat, either linear shifting or stable
%large tau & large Amp quadrant => linear shifting
% large tau & small Amp => stable PF
% small tau and large Amp => early large shift (if tau<1 => abrupt shift on first lap)
% small tau and small Amp => early small shift

figure
% distributions of Amp for each condition (pooled across animals)
subplot(2,2,1)
for n = 1:size(CA_VR_idx)
hCAp1{n} = histogram(Tok.RCreg_p1(CA_VR_idx{n}), 40, 'Normalization', 'cdf', 'DisplayStyle', 'stairs', 'EdgeColor', colors(n,:), 'FaceColor', 'none', 'LineWidth', 1.5); hold on
text(0, 0.005+n*0.05, ['mean = ' num2str(round( mean(Tok.RCreg_p1(CA_VR_idx{n})),2,'significant') ) ' cm' ], 'Color', colors(n,:) );
end
xlabel('nlfit Amp (cm)'); ylabel('cdf');
box off; axis square;

% distributions of abs(Amp) for each condition (pooled across animals)
subplot(2,2,2)
for n = 1:size(CA_VR_idx)
hCAp1{n} = histogram(abs(Tok.RCreg_p1(CA_VR_idx{n})), 40, 'Normalization', 'cdf', 'DisplayStyle', 'stairs', 'EdgeColor', colors(n,:), 'FaceColor', 'none', 'LineWidth', 1.5); hold on
text(0, 0.005+n*0.05, ['median = ' num2str(round( median(abs(Tok.RCreg_p1(CA_VR_idx{n}))),2,'significant') ) ' cm' ], 'Color', colors(n,:) );
end
xlabel('nlfit abs(Amp) (cm)'); ylabel('cdf');
box off; axis square;

% distributions of time constant for each conditions (pooled across animals)
subplot(2,2,3)
for n = 1:size(CA_VR_idx)
hCAp2{n} = histogram(Tok.RCreg_p2(CA_VR_idx{n}),50, 'Normalization', 'cdf', 'DisplayStyle', 'stairs', 'EdgeColor', colors(n,:), 'FaceColor', 'none', 'LineWidth', 1.5); hold on
text(50, 0.003+n*0.05, ['median = ' num2str(round( median(Tok.RCreg_p2(CA_VR_idx{n})),2,'significant') ) ' laps' ], 'Color', colors(n,:) );
end
ylim([0 1])
xlabel('nlfit Tau (laps)'); ylabel('cdf');
box off; axis square;

% distributions of Rsquare for each conditions (pooled across animals)
subplot(2,2,4)
for n = 1:size(CA_VR_idx)
RCregR2adj{n} = Tok.RCreg_R2adj(CA_VR_idx{n});
RCregR2adj{n}(RCregR2adj{n}<0) = 0;
hCAR2{n} = histogram(RCregR2adj{n}, 40, 'Normalization', 'cdf', 'DisplayStyle', 'stairs', 'EdgeColor', colors(n,:), 'FaceColor', 'none', 'LineWidth', 1.5); hold on
text(50, 0.003+n*0.05, ['median = ' num2str(round( median(Tok.RCreg_R2adj(CA_VR_idx{n})),2,'significant') ) ' laps' ], 'Color', colors(n,:) );
end
% xlim([0 1])
xlabel('nlfit Rsquare'); ylabel('cdf');
box off; axis square;

% animal-wise stats for tau, Amp or abs(Amp), Rsq, prop of Fwd, prop Bck
% for tau: a) prop (or mean?) of PFs < 15 (or 20?). Threshold defined as
% the elbow of the distribution of taus. 
tau_thresh = 10;
pctl = 66;
% CA1mice = unique(Tok.animalID(CA1_N_idx));
% CA3mice = unique(Tok.animalID(CA3_N_idx));
% CAmice = {CA1mice; CA3mice};
% for c = 1:length(CAmice)
for m = 1:length(CA1mice)
    nl_IdxFwd_CA1F{m} = find(Tok.RCreg_p1>0 & ismember(Tok.animalID, CA1mice{m}) & ismember(Tok.f_n, 'f'));
    nl_PropFwd_CA1F(m) = length(nl_IdxFwd_CA1F{m})/length_CA1F(m);
    nl_PropBck_CA1F(m) = 1-nl_PropFwd_CA1F(m);
    propTauThresh_CA1F(m) = length(find(Tok.RCreg_p2(IdxCA1F{m}) < tau_thresh))/length(IdxCA1F{m});
    medianTau_CA1F(m) = median(Tok.RCreg_p2(IdxCA1F{m}));
    pctTau_CA1F(m) = prctile(Tok.RCreg_p2(IdxCA1F{m}),pctl);
    medianAbsAmp_CA1F(m) = median(abs(Tok.RCreg_p1(IdxCA1F{m})));
    meanAmp_CA1F(m) = mean(Tok.RCreg_p1(IdxCA1F{m}));
    RCregR2adj_CA1F{m} = Tok.RCreg_R2adj(IdxCA1F{m});
    RCregR2adj_CA1F{m}(RCregR2adj_CA1F{m}<0) = 0;
    medianR2_CA1F(m) = median(RCregR2adj_CA1F{m});
    meanPC1_CA1F(m) = mean(score(IdxCA1F{m}));
    medianAbsPC1_CA1F(m) = median(abs(score(IdxCA1F{m})));
    score_CA1F = score(IdxCA1F{m});
    medianPC1fwd_CA1F(m) = median(score_CA1F(score_CA1F>=0));
    medianPC1bck_CA1F(m) = median(score_CA1F(score_CA1F<0));
    propPC1fwd_CA1F(m) = length(score_CA1F(score_CA1F>=0))./length(score_CA1F);
    
    nl_IdxFwd_CA1N{m} = find(Tok.RCreg_p1>0 & ismember(Tok.animalID, CA1mice{m}) & ismember(Tok.f_n, 'n'));
    nl_PropFwd_CA1N(m) = length(nl_IdxFwd_CA1N{m})/length_CA1N(m);
    nl_PropBck_CA1N(m) = 1-nl_PropFwd_CA1N(m);
    propTauThresh_CA1N(m) = length(find(Tok.RCreg_p2(IdxCA1N{m}) < tau_thresh))/length(IdxCA1N{m});
    medianTau_CA1N(m) = median(Tok.RCreg_p2(IdxCA1N{m}));
    pctTau_CA1N(m) = prctile(Tok.RCreg_p2(IdxCA1N{m}),pctl);
    medianAbsAmp_CA1N(m) = median(abs(Tok.RCreg_p1(IdxCA1N{m})));
    meanAmp_CA1N(m) = mean(Tok.RCreg_p1(IdxCA1N{m}));
    RCregR2adj_CA1N{m} = Tok.RCreg_R2adj(IdxCA1N{m});
    RCregR2adj_CA1N{m}(RCregR2adj_CA1N{m}<0) = 0;
    medianR2_CA1N(m) = median(RCregR2adj_CA1N{m});
    meanPC1_CA1N(m) = mean(score(IdxCA1N{m}));
    medianAbsPC1_CA1N(m) = median(abs(score(IdxCA1N{m})));
    score_CA1N = score(IdxCA1N{m});
    medianPC1fwd_CA1N(m) = median(score_CA1N(score_CA1N>=0));
    medianPC1bck_CA1N(m) = median(score_CA1N(score_CA1N<0));
    propPC1fwd_CA1N(m) = length(score_CA1N(score_CA1N>=0))./length(score_CA1N);

    clear scoreCA1F scoreCA1N
end

for m = 1:length(CA3mice)
    nl_IdxFwd_CA3F{m} = find(Tok.RCreg_p1>0 & ismember(Tok.animalID, CA3mice{m}) & ismember(Tok.f_n, 'f'));
    nl_PropFwd_CA3F(m) = length(nl_IdxFwd_CA3F{m})/length_CA3F(m);
    nl_PropBck_CA3F(m) = 1-nl_PropFwd_CA3F(m);
    propTauThresh_CA3F(m) = length(find(Tok.RCreg_p2(IdxCA3F{m}) < tau_thresh))/length(IdxCA3F{m});
    medianTau_CA3F(m) = median(Tok.RCreg_p2(IdxCA3F{m}));
    pctTau_CA3F(m) = prctile(Tok.RCreg_p2(IdxCA3F{m}),pctl);
    medianAbsAmp_CA3F(m) = median(abs(Tok.RCreg_p1(IdxCA3F{m})));
    meanAmp_CA3F(m) = mean(Tok.RCreg_p1(IdxCA3F{m}));
    RCregR2adj_CA3F{m} = Tok.RCreg_R2adj(IdxCA3F{m});
    RCregR2adj_CA3F{m}(RCregR2adj_CA3F{m}<0) = 0;
    medianR2_CA3F(m) = median(RCregR2adj_CA3F{m});
    meanPC1_CA3F(m) = mean(score(IdxCA3F{m}));
    medianAbsPC1_CA3F(m) = median(abs(score(IdxCA3F{m})));
    
    score_CA3F = score(IdxCA3F{m});
    medianPC1fwd_CA3F(m) = median(score_CA3F(score_CA3F>=0));
    medianPC1bck_CA3F(m) = median(score_CA3F(score_CA3F<0));
    propPC1fwd_CA3F(m) = length(score_CA3F(score_CA3F>=0))./length(score_CA3F);
    
    nl_IdxFwd_CA3N{m} = find(Tok.RCreg_p1>0 & ismember(Tok.animalID, CA3mice{m}) & ismember(Tok.f_n, 'n'));
    nl_PropFwd_CA3N(m) = length(nl_IdxFwd_CA3N{m})/length_CA3N(m);
    nl_PropBck_CA3N(m) = 1-nl_PropFwd_CA3N(m);
    propTauThresh_CA3N(m) = length(find(Tok.RCreg_p2(IdxCA3N{m}) < tau_thresh))/length(IdxCA3N{m});
    medianTau_CA3N(m) = median(Tok.RCreg_p2(IdxCA3N{m}));
    pctTau_CA3N(m) = prctile(Tok.RCreg_p2(IdxCA3N{m}),pctl);
    medianAbsAmp_CA3N(m) = median(abs(Tok.RCreg_p1(IdxCA3N{m})));
    meanAmp_CA3N(m) = mean(Tok.RCreg_p1(IdxCA3N{m}));
    RCregR2adj_CA3N{m} = Tok.RCreg_R2adj(IdxCA3N{m});
    RCregR2adj_CA3N{m}(RCregR2adj_CA3N{m}<0) = 0;
    medianR2_CA3N(m) = median(RCregR2adj_CA3N{m});
    meanPC1_CA3N(m) = mean(score(IdxCA3N{m}));
    medianAbsPC1_CA3N(m) = median(abs(score(IdxCA3N{m})));
    
    score_CA3N = score(IdxCA3N{m});
    medianPC1fwd_CA3N(m) = median(score_CA3N(score_CA3N>=0));
    medianPC1bck_CA3N(m) = median(score_CA3N(score_CA3N<0));
    propPC1fwd_CA3N(m) = length(score_CA3N(score_CA3N>=0))./length(score_CA3N);
    
    clear scoreCA3F scoreCA3N
end
% end

propTau_CA1 = [propTauThresh_CA1N; propTauThresh_CA1F];
propTau_CA3 = [propTauThresh_CA3N; propTauThresh_CA3F];
medianTau_CA1 = [medianTau_CA1N; medianTau_CA1F];
medianTau_CA3 = [medianTau_CA3N; medianTau_CA3F];
pctTau_CA1 = [pctTau_CA1N; pctTau_CA1F];
pctTau_CA3 = [pctTau_CA3N; pctTau_CA3F];
medianAbsAmp_CA1 = [medianAbsAmp_CA1N;medianAbsAmp_CA1F];
medianAbsAmp_CA3 = [medianAbsAmp_CA3N;medianAbsAmp_CA3F];
meanAmp_CA1 = [meanAmp_CA1N; meanAmp_CA1F];
meanAmp_CA3 = [meanAmp_CA3N; meanAmp_CA3F];
medianR2_CA1 = [medianR2_CA1N; medianR2_CA1F];
medianR2_CA3 = [medianR2_CA3N; medianR2_CA3F];
nlPropBck_CA1 = [nl_PropBck_CA1N; nl_PropBck_CA1F];
nlPropBck_CA3 = [nl_PropBck_CA3N; nl_PropBck_CA3F];
meanPC1_CA1 = [meanPC1_CA1N; meanPC1_CA1F];
meanPC1_CA3 = [meanPC1_CA3N; meanPC1_CA3F];
medianAbsPC1_CA1 = [medianAbsPC1_CA1N; medianAbsPC1_CA1F];
medianAbsPC1_CA3 = [medianAbsPC1_CA3N; medianAbsPC1_CA3F];
medianPC1fwd_CA1 = [medianPC1fwd_CA1N; medianPC1fwd_CA1F];
medianPC1fwd_CA3 = [medianPC1fwd_CA3N; medianPC1fwd_CA3F];
medianPC1bck_CA1 = [medianPC1bck_CA1N; medianPC1bck_CA1F];
medianPC1bck_CA3 = [medianPC1bck_CA3N; medianPC1bck_CA3F];
propPC1fwd_CA1 = [propPC1fwd_CA1N; propPC1fwd_CA1F];
propPC1fwd_CA3 = [propPC1fwd_CA3N; propPC1fwd_CA3F];

%figures
RCregStats_CA1 = {propTau_CA1; medianTau_CA1; pctTau_CA1; medianAbsAmp_CA1; meanAmp_CA1; medianR2_CA1; nlPropBck_CA1; meanPC1_CA1; medianAbsPC1_CA1; medianPC1fwd_CA1; medianPC1bck_CA1; propPC1fwd_CA1};
RCregStats_CA3 = {propTau_CA3; medianTau_CA3; pctTau_CA3; medianAbsAmp_CA3; meanAmp_CA3; medianR2_CA3; nlPropBck_CA3; meanPC1_CA3; medianAbsPC1_CA3; medianPC1fwd_CA3; medianPC1bck_CA3; propPC1fwd_CA3};
RCregStats = { ['ratio of PFs with tau < ' num2str(tau_thresh) 'laps']; 'median tau, laps'; [num2str(pctl) 'th pctile tau, laps']; 'median abs(Amp), cm'; 'mean Amp, cm'; 'nl adj R-square'; 'ratio of PFs w/ Amp < 0'; 'mean PC1'; 'median Abs(PC1)'; 'median PC1>0'; 'median PC1<0'; 'proportion PC1>0' };

figure
RCvar1 = [1 2 4 5];
for p = 1:length(RCvar1) %1:length(RCregStats)-1
subplot(2,2,p)
n = RCvar1(p);
plot(RCregStats_CA1{n},'-', 'Color', mapCA1(2,:)); hold on
plot(RCregStats_CA3{n},'-', 'Color', mapCA3(2,:)); hold on
plot(mean(RCregStats_CA1{n},2), 'ok-', 'LineWidth', 0.75, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', mapCA1(2,:), 'MarkerSize', 8); hold on
plot(mean(RCregStats_CA3{n},2),'dk-', 'LineWidth', 0.75, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', mapCA3(2,:), 'MarkerSize', 8); hold on
% legend('CA1 means', 'CA3 means', 'CA1 individuals', 'CA3 individuals')
xlim([0 3]);
if min([RCregStats_CA1{n}, RCregStats_CA3{n}], [], 'all')>0
    ylim([0 max([RCregStats_CA1{n}, RCregStats_CA3{n}], [], 'all')+0.1*max([RCregStats_CA1{n}, RCregStats_CA3{n}], [], 'all')])
end
ylabel(RCregStats{n}); 
set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
box off; 
axis square
end

figure
RCvar2 = [3 6:12];
for p = 1:length(RCvar2)
subplot(3,3,p)
n = RCvar2(p);
plot(nanmean(RCregStats_CA1{n},2),'o-', 'LineWidth', 2, 'Color', mapCA1(2,:)); hold on
plot(nanmean(RCregStats_CA3{n},2),'d-', 'LineWidth', 2, 'Color', mapCA3(2,:)); hold on
plot(RCregStats_CA1{n},'-', 'Color', mapCA1(2,:)); hold on
plot(RCregStats_CA3{n},'-', 'Color', mapCA3(2,:)); hold on
% legend('CA1 means', 'CA3 means', 'CA1 individuals', 'CA3 individuals')
xlim([0 3]);
if min([RCregStats_CA1{n}, RCregStats_CA3{n}], [], 'all')>0
    ylim([0 max([RCregStats_CA1{n}, RCregStats_CA3{n}], [], 'all')+0.1*max([RCregStats_CA1{n}, RCregStats_CA3{n}], [], 'all')])
end
ylabel(RCregStats{n}); 
set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
box off; 
axis square
end

figure
RCvar3 = [9];
for p = 1:length(RCvar3)
% subplot(1,2,p)
n = RCvar3(p);
plot(RCregStats_CA1{n},'-', 'Color', mapCA1(2,:), 'LineWidth', 0.5); hold on
plot(RCregStats_CA3{n},'-', 'Color', mapCA3(2,:), 'LineWidth', 0.5); hold on
plot(mean(RCregStats_CA1{n},2),'ok-', 'LineWidth', 0.75, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', mapCA1(2,:), 'MarkerSize', 15); hold on
plot(mean(RCregStats_CA3{n},2),'dk-', 'LineWidth', 0.75, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', mapCA3(2,:), 'MarkerSize', 15); hold on
% legend('CA1 means', 'CA3 means', 'CA1 individuals', 'CA3 individuals')
xlim([0 3]);
if min([RCregStats_CA1{n}, RCregStats_CA3{n}], [], 'all')>=0
    ylim([0 max([RCregStats_CA1{n}, RCregStats_CA3{n}], [], 'all')+0.1*max([RCregStats_CA1{n}, RCregStats_CA3{n}], [], 'all')])
end
ylabel(RCregStats{n}); 
set(gca, 'XTick', [1, 2], 'XTickLabel', ['N';'F']);
box off; 
axis square
end

%linear mixed effect model 
% mice = [CA1mice;CA1mice;CA3mice;CA3mice];
% F_N = [repmat({'N'}, length(CA1mice), 1);repmat({'F'}, length(CA1mice), 1); repmat({'N'}, length(CA3mice), 1); repmat({'F'}, length(CA3mice), 1) ];
% CA = [repmat({'CA1'}, 2*length(CA1mice), 1);repmat({'CA3'}, 2*length(CA3mice), 1)];
medianAbsAmp = [medianAbsAmp_CA1N, medianAbsAmp_CA1F, medianAbsAmp_CA3N, medianAbsAmp_CA3F]';
medianTau = [medianTau_CA1N, medianTau_CA1F, medianTau_CA3N, medianTau_CA3F]';
propTau = [propTauThresh_CA1N, propTauThresh_CA1F, propTauThresh_CA3N, propTauThresh_CA3F]';
propPC1fwd = [propPC1fwd_CA1N, propPC1fwd_CA1F, propPC1fwd_CA3N, propPC1fwd_CA3F]';
medianAbsPC1 = [medianAbsPC1_CA1N, medianAbsPC1_CA1F, medianAbsPC1_CA3N, medianAbsPC1_CA3F]';

% tbl1 = table(mice, CA, F_N, PropShift, Slopes, 'VariableNames',{'mice', 'CA', 'F_N', 'PropShift', 'Slopes'});
tbl1.AbsAmp = medianAbsAmp;
tbl1.Tau = medianTau;
tbl1.propTau = propTau;
tbl1.propPC1fwd = propPC1fwd;
tbl1.medianAbsPC1 = medianAbsPC1;

lme_Tau = fitlme(tbl1,'Tau ~ 1 + CA + F_N + (1 + F_N | mice)'); % interaction CA * F_N excluded because not signif. 
anova(lme_Tau) % Fixed factors CA ns (p = 0.865) and F_N ns (p = 0.098)

lme_propTau = fitlme(tbl1,'propTau ~ 1 + CA + F_N + (1 + F_N | mice)'); % interaction CA * F_N excluded because not signif. 
anova(lme_propTau) % Fixed factors CA ns (p = 0.349) and F_N sig (p = 0.038)

lme_Amp = fitlme(tbl1,'AbsAmp ~ 1 + CA + F_N + (1 + F_N | mice)'); % interaction CA * F_N excluded because not signif. 
anova(lme_Amp) % Fixed factors CA sig (p = 0.0053) and F_N sig (p < 0.0001)

lme_AbsPC1 = fitlme(tbl1,'medianAbsPC1 ~ 1 + CA + F_N + (1 + F_N | mice)'); % interaction CA * F_N excluded because not signif. 
anova(lme_AbsPC1) % Fixed factors CA sig (p = 0.0095) and F_N sig (p = 0.00024)

lme_propPC1 = fitlme(tbl1,'propPC1fwd ~ 1 + CA + F_N + (1 + F_N | mice)'); % interaction CA * F_N excluded because not signif. (although it looks like it could, if more n)
anova(lme_propPC1) % fixed factors CA and F_N not significant

% ? need to write a function to perform resampling test (generalize to more
% than 2 distrib) and plot in dabest mode, with possibility to choose
% statistic (mean, median or other percentile). 

% how many PFs are more linear, more nonlin, or well explained by both, or
% by neither?

% plot all trajectories with nonlinear Rsquare > 0.5 and better than linear
% R2

% plot all trajectories with nonlinear Rsquare better than linear R2 and
% linear pval > 0.05 (i.e, those that were not explained by  

% scatter plot of those PFs with good nonlinear fits: linear slope vs p(1) and p(2)

% plot proportion of: lin forward (best fit), lin backward (best fit), nlinforward, nlin backward, and nonexplained (would need to access the Ftest pval of fitnlm for that...) 

%% lap-by-lap COM shift on onset-centered trajectories (interpolated): do large shift occur after 1st five laps post-onset? their fqcy?
% this will get at whether there are potentially BTSP events occuring late
% after onset

DcomAll2 = diff(COMtrajMat2,1,2); % on interpolated data
% DcomAll2 = diff(COMtrajMatNaN2,1,2); % on raw data

xPOL = 1:Minlaps-1;
Xdcom2 = repmat(1:Minlaps-1, size(DcomAll2,1),1);

% check for all different conditions
Dcom2_CA1F = diff(COMtrajMat2(CA1_F_idx,:),1,2); % on interpolated data
Dcom2_CA1N = diff(COMtrajMat2(CA1_N_idx,:),1,2);
Dcom2_CA3F = diff(COMtrajMat2(CA3_F_idx,:),1,2);
Dcom2_CA3N = diff(COMtrajMat2(CA3_N_idx,:),1,2);

% Dcom2_CA1F = diff(COMtrajMatNaN2(CA1_F_idx,:),1,2); % on raw data
% Dcom2_CA1N = diff(COMtrajMatNaN2(CA1_N_idx,:),1,2);
% Dcom2_CA3F = diff(COMtrajMatNaN2(CA3_F_idx,:),1,2);
% Dcom2_CA3N = diff(COMtrajMatNaN2(CA3_N_idx,:),1,2);

% number of PFs and frequency, per lap, with shift higher (or lower) than threshold
SThresh = [3, 6, 10, 20]; 
SThresh_string = string(SThresh);

for s = 1:length(SThresh)
for l = 1:size(DcomAll2,2)
nPFsup(s,l) = length(find(abs(DcomAll2(:,l)) > SThresh(s))); 
nPFsub(s,l) = length(find(abs(DcomAll2(:,l)) <= SThresh(s)));

nPFsup_CA1F(s,l) = length(find(abs(Dcom2_CA1F(:,l)) > SThresh(s)));
nPFsup_CA1N(s,l) = length(find(abs(Dcom2_CA1N(:,l)) > SThresh(s)));
nPFsup_CA3F(s,l) = length(find(abs(Dcom2_CA3F(:,l)) > SThresh(s)));
nPFsup_CA3N(s,l) = length(find(abs(Dcom2_CA3N(:,l)) > SThresh(s))); 
end
end


t7 = tiledlayout('flow', 'TileSpacing','Compact','Padding','Compact');

ax1 = nexttile;
colormap(ax1, parula(256))
h = histogram2(Xdcom2,DcomAll2,'DisplayStyle','tile','ShowEmptyBins','off', 'Normalization', 'count');
xlabel('post-onset lap'); ylabel('lap-to-lap COM shift (cm)');
ylim([-50, 50])
% title('lapwise Delta COM')
box off; axis square;
cb = colorbar;
cb.Label.String = 'number of PFs'; cb.Limits = [0 max(h.Values, [], 'all')];

for s = 1:length(SThresh)
nexttile
plot(xPOL, nPFsup_CA1F(s,:)./length(Dcom2_CA1F(:,1)), 'LineWidth', 1,'Color', mapCA1(1,:));hold on
plot(xPOL, nPFsup_CA1N(s,:)./length(Dcom2_CA1N(:,1)), 'LineWidth', 1,'Color', mapCA1(2,:));hold on
plot(xPOL, nPFsup_CA3F(s,:)./length(Dcom2_CA3F(:,1)), 'LineWidth', 1,'Color', mapCA3(1,:));hold on
plot(xPOL, nPFsup_CA3N(s,:)./length(Dcom2_CA3N(:,1)), 'LineWidth', 1,'Color', mapCA3(2,:));hold on
plot(xPOL, nPFsup(s,:)./length(DcomAll2(:,1)), 'k-', 'LineWidth',2);
yline(0.5, 'k--')
xlim([0 Minlaps]); ylim([0 1]);
xlabel('post-onset lap'); ylabel('frequency of PFs');
title(['Shifts > ' num2str(SThresh(s)) ' cm'])
% legend(['CA1F >' num2str(SThresh(s)) ' cm'], 'CA1N >', 'CA3F >', 'CA3N >', 'All >', 'Location', 'Best' )
box off; axis square;
end

ax2 = nexttile;
for s = 1:length(SThresh)
plot(xPOL, nPFsup(s,:)./length(DcomAll2(:,1))); hold on
end
% plot(xPOL, nPFsub./length(DcomAll2(:,1)), 'r-');
xlim([0 Minlaps]); ylim([0 1]);
xlabel('post-onset lap'); ylabel('frequency of PFs');
legend(SThresh_string, 'Location', 'Best')
% legend(['< ' num2str(SThresh) ' cm'], ['> ' num2str(SThresh) ' cm'])
box off; axis square;

% exportgraphics(t7,'Lap2LapShifts_multThresh.pdf','BackgroundColor','none', 'ContentType','vector')


% f7.Renderer = 'painters';
% print(f7, '-vector','-dpdf','Lap2LapShifts_multThresh.pdf')

% CCL: larger shifts can occur any time after onset, but are more frequent in the beginning. it seems there is a trade off in early laps between the 2 first
% thresholds: a progressive increase for the smallest shifts vs larger shift in 1st lap and a
% progressive decrease from then on that stabilizes somewhat after 5 laps. => plot as lines the % of shifts under and over that threshold.

% for l = 1:Minlaps-1
% nPFsup_CA1F(l) = length(find(abs(Dcom2_CA1F(:,l)) > SThresh)); 
% nPFsub_CA1F(l) = length(find(abs(Dcom2_CA1F(:,l)) <= SThresh));
% 
% nPFsup_CA1N(l) = length(find(abs(Dcom2_CA1N(:,l)) > SThresh)); 
% nPFsub_CA1N(l) = length(find(abs(Dcom2_CA1N(:,l)) <= SThresh));
% 
% nPFsup_CA3F(l) = length(find(abs(Dcom2_CA3F(:,l)) > SThresh)); 
% nPFsub_CA3F(l) = length(find(abs(Dcom2_CA3F(:,l)) <= SThresh));
% 
% nPFsup_CA3N(l) = length(find(abs(Dcom2_CA3N(:,l)) > SThresh)); 
% nPFsub_CA3N(l) = length(find(abs(Dcom2_CA3N(:,l)) <= SThresh));
% end
% figure
% plot(2:Minlaps, nPFsub_CA1F./length(Dcom2_CA1F(:,1)),'--', 'LineWidth',1,'Color', mapCA1(1,:)); hold on
% plot(2:Minlaps, nPFsup_CA1F./length(Dcom2_CA1F(:,1)), 'LineWidth', 1,'Color', mapCA1(1,:));hold on
% plot(2:Minlaps, nPFsub_CA1N./length(Dcom2_CA1N(:,1)),'--','LineWidth',1,'Color', mapCA1(2,:)); hold on
% plot(2:Minlaps, nPFsup_CA1N./length(Dcom2_CA1N(:,1)), 'LineWidth', 1,'Color', mapCA1(2,:));hold on
% plot(2:Minlaps, nPFsub_CA3F./length(Dcom2_CA3F(:,1)),'--', 'LineWidth',1,'Color', mapCA3(1,:)); hold on
% plot(2:Minlaps, nPFsup_CA3F./length(Dcom2_CA3F(:,1)), 'LineWidth', 1,'Color', mapCA3(1,:));hold on
% plot(2:Minlaps, nPFsub_CA3N./length(Dcom2_CA3N(:,1)),'--', 'LineWidth',1,'Color', mapCA3(2,:)); hold on
% plot(2:Minlaps, nPFsup_CA3N./length(Dcom2_CA3N(:,1)), 'LineWidth', 1,'Color', mapCA3(2,:));hold on
% plot(2:Minlaps, nPFsub./length(DcomAll2(:,1)), 'k--', 'LineWidth',2); hold on
% plot(2:Minlaps, nPFsup./length(DcomAll2(:,1)), 'k-', 'LineWidth',2);
% xlim([1 Minlaps]); ylim([0 1]);
% xlabel('post-onset lap'); ylabel('frequency of PFs');
% legend(['CA1F <' num2str(SThresh) ' cm'], 'CA1F >', 'CA1N <', 'CA1N >', 'CA3F <', 'CA3F >', 'CA3N <', 'CA3N >', 'All <', 'All >', 'Location', 'BestOutside' )
% box off; axis square;

%% absolute lap-wise pop analysis: same as above but on running laps (i.e. centered on VR start rather than PF onset, and not interpolated)

% plot PC1 score against onset lap: when do PFs with large shift on first 5 laps after onset occur? 
figure
% subplot(1,2,1)
% scatter(COMonset, score(:,1))
% xlabel('COM onset lap'); ylabel('PF PC1 score');
% box off; axis square;
% subplot(1,2,2)
histogram2(T.COMonset(~isnan(COMcol1)), score(:,1), 'DisplayStyle','tile','ShowEmptyBins','off')
xlabel('COM onset lap'); ylabel('PF PC1 score');
box off; axis square;
cb = colorbar;
cb.Label.String = 'PFs'
%normalize with number of PFs per onset lap? -> the normalization I think
%would show that the distribution of PC1 scores is independent of onset
%laps. There's just less of PFs appearing later on, so less opportunities to get a
%large PC1 score. 

% same as above but comparing CA_VR groups

% Are there sudden large shifts and how frequent are they? 

CA1F = 1*ones(1,numel(DcomAll{1,1})); %group1
CA1N = 2*ones(1, numel(DcomAll{1,2})); %group2
CA3F = 3*ones(1, numel(DcomAll{2,1})); %group3
CA3N = 4*ones(1, numel(DcomAll{2,2})); %group4
CAVRall = [CA1F, CA1N, CA3F, CA3N];

% what is a large shift?
% distribution of lap-by-lap shifts across all conditions, Minlaps2 # of laps and PFs 
% (Can did that analysis for 25 laps and significantly linearly shifting PFs only (probably just CA1 F or N) 
% => she found a similar gaussian, mostly between -50 and 50cm centered on something close to 0) 
figure 
subplot(1,2,1)
histogram(Allshifts, 'Normalization', 'count');
xlabel('lap-wise COM shift (cm)'); ylabel('PFs');
box off; axis square;
subplot(1,2,2)
violinplot(Allshifts, CAVRall)
xlabel('CA-VR groups'); ylabel('COM shifts');
box off; axis square;


% lap-wise COM shift distrib as a function of lap
GpToPlot = DcomAll{1,2}; % change the dataset to plot here to affect all following plots
Xdcom = repmat(2:Minlaps2, size(GpToPlot,1),1);
Xedges = [1.5:Minlaps2+0.5];

% figure
% subplot(1,2,1)
% scatter(Xdcom, GpToPlot)
% xlabel('running lap'); ylabel('COM shift (cm)');
% box off; axis square;
% subplot(1,2,2)
% scatter(Xdcom, GpToPlot)
% xlabel('running lap'); ylabel('COM shift');
% ylim([-20, 20])
% box off; axis square;

% figure
%     Yedges = [-250, -100, -50, -25, -10, -5, -1, 1, 5, 10, 25, 50, 100, 250];
%     % Yedges1 = logspace(0,2,15);
%     % Yedges = sort(2.5*[-Yedges1, Yedges1]);
%     h1 = histogram2(Xdcom,GpToPlot,Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','off', 'Normalization', 'countdensity');
%     xlabel('running lap'); ylabel('COM shift');
% %     ylim([-25, 25])
%     box off; axis square;
%     cb1 = colorbar;
%     cb1.Label.String = 'density of PFs'; cb1.Limits = [0 max(h1.Values, [], 'all')];

figure
    Yedges = [-250:5:250];
    h2 = histogram2(Xdcom,GpToPlot,Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on', 'Normalization', 'count');
    xlabel('running lap'); ylabel('COM shift');
    ylim([-50, 50])
    box off; axis square;
    cb2 = colorbar;
    cb2.Label.String = 'number of PFs'; cb2.Limits = [0 max(h2.Values, [], 'all')];

figure
    PFsPerLap = sum(h2.Values, 2);
    PFsPLmat = repmat(PFsPerLap,1,size(h2.Values,2));
    NormHisto = h2.Values./PFsPLmat;
    imagesc(rot90(NormHisto));
    set(gca, 'YTick', [40.5,50.5,60.5], 'YTickLabel', [' 50'; '  0'; '-50']);
    xlabel('running lap'); ylabel('COM shift (cm)');
    ylim([40.5, 60.5]) % corresponds to shifts of -50 and 50 in Yedges -> replace tick labels later
    box off; axis square;
    cb3 = colorbar;
    cb3.Label.String = 'PFs count / total PFs per lap'; cb3.Limits = [0 max(NormHisto, [], 'all')];

figure % mean lap-wise COM shift (should lead to same ccl than in Can's paper) -> it does, although it's more noisy because she plotted a moving average to smooth COM location
DcomMean_CA1F = mean(DcomAll{1,1}, 1, 'omitnan');
DcomMean_CA1N = mean(DcomAll{1,2}, 1, 'omitnan');
DcomMean_CA3F = mean(DcomAll{2,1}, 1, 'omitnan');
DcomMean_CA3N = mean(DcomAll{2,2}, 1, 'omitnan');
plot(2:Minlaps2, DcomMean_CA1N, 'b-', 2:Minlaps2, DcomMean_CA1F, 'c-', 2:Minlaps2, DcomMean_CA3N, 'r-', 2:Minlaps2, DcomMean_CA3F, 'm-');
legend('CA1F','CA1N','CA3F', 'CA3N')
xlabel('running lap'); ylabel('mean COM shift (cm)');

figure
subplot(1,2,1)
scatter(T.maxShift_Lap, T.maxShift)
xlabel('running lap'); ylabel('max COM shift per PF (cm/lap)');
box off; axis square;
subplot(1,2,2)
histogram2(T.maxShift_Lap, T.maxShift,'DisplayStyle','tile','ShowEmptyBins','off', 'Normalization', 'count');
xlabel('running lap'); ylabel('max COM shift per PF (cm/lap)');
box off; axis square;

%% Speed vs Shift

figure
scatter(T.maxS_PrevLapSpeed, T.maxS_LapSpeed)
xlabel('Previous Lap Speed (cm/s)'); ylabel('Lap Speed');
box off; axis square;
%CCL: speed on lap and previous lap are correlated

figure
scatter(T.maxS_PrevLapSpeed, T.maxShift)
xlabel('Previous Lap Speed (cm/s)'); ylabel('max absolute COM shift per PF (cm/lap)');
box off; axis square;

figure
histogram2(T.maxS_PrevLapSpeed, T.maxShift, 'DisplayStyle','tile','ShowEmptyBins','off')
xlabel('Previous Lap Speed (cm/s)'); ylabel('max absolute COM shift per PF (cm/lap)');
box off; axis square;
cb = colorbar;
cb.Label.String = 'PFs';

figure
subplot(1,2,1)
scatter(T.maxS_LapSpeed, T.maxShift)
xlabel('Lap Speed (cm/s)'); ylabel('max absolute COM shift per PF (cm/lap)');
box off; axis square;
subplot(1,2,2)
histogram2(T.maxS_LapSpeed, T.maxShift, 'DisplayStyle','tile','ShowEmptyBins','off')
xlabel('Lap Speed (cm/s)'); ylabel('max absolute COM shift per PF (cm/lap)');
xlim([0 max(T.maxS_LapSpeed)])
box off; axis square;
% cb = colorbar;
% cb.Label.String = 'PFs';

figure
boxchart(T.animalSpeed,abs(T.meanShift))
xlabel('mean animal speed (cm/s)'); ylabel('absolute mean COM shift per PF (cm/lap)');
box off; axis square;

figure
boxchart(T.animalSpeed,abs(T.COM_slope)); hold on
xlabel('mean animal speed (cm/s)'); ylabel('absolute COM slope per PF (cm/lap)');
box off; axis square;

figure
boxchart(T.animalSpeed(T.COM_pval<0.05),abs(T.COM_slope(T.COM_pval<0.05))); hold on
[B,~,~,~,StatsB] = regress(abs(T.COM_slope(T.COM_pval<0.05)), [ones(size(T.animalSpeed(T.COM_pval<0.05))), T.animalSpeed(T.COM_pval<0.05)]); % regression
Yfit = B(1) + B(2)*[min(T.animalSpeed(T.COM_pval<0.05)):max(T.animalSpeed(T.COM_pval<0.05))];
plot([min(T.animalSpeed(T.COM_pval<0.05)):max(T.animalSpeed(T.COM_pval<0.05))],Yfit, 'g-', 'markersize', 2);
xlabel('mean animal speed (cm/s)'); ylabel('absolute COM slope per PF (cm/lap)');
title(['sig PFs: R-sqare = ' num2str(StatsB(1)) ' p = ' num2str(StatsB(3))])
box off; axis square;

% same plots but for proportion of sig shifting PFs rather than actual shift speed

% do linear model to check influence of various variables (CA3 vs CA1, lap speed, f vs n) on absolute shift

%% local functions

function ratio = propshift(pvalArray)
% computes the proportion of significantly shifting PFs
ratio = length(find(pvalArray<=0.05))./length(pvalArray);
end

function ratio = propBack(pval, slope)
% computes the proportion of significantly backward shifting PFs
%pval and slope must be the same length
ratio = length(pval(pval<=0.05 & slope<0))./length(pval);
end

function ratio = propFwd(pval, slope)
% computes the proportion of significantly backward shifting PFs
%pval and slope must be the same length
ratio = length(pval(pval<=0.05 & slope>0))./length(pval);
end
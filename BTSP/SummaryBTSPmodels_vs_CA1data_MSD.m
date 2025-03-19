clear 
close all

Minlaps = 30; % min number of laps where PF is defined (after interpolation, if interp is used)
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

DispCA1N = ( TrajsCA1N-TrajsCA1N(:,end) );
DispCA1F = ( TrajsCA1F-TrajsCA1F(:,end) );
DispCA3N = ( TrajsCA3N-TrajsCA3N(:,end) );
DispCA3F = ( TrajsCA3F-TrajsCA3F(:,end) );
mdCA1F = mean( DispCA1F , 1 );
mdCA1N = mean( DispCA1N, 1 );
mdCA3F = mean( DispCA3F, 1 );
mdCA3N = mean( DispCA3N, 1 );
madCA1F = mean( abs(DispCA1F) , 1 );
madCA1N = mean( abs(DispCA1N), 1 );
madCA3F = mean( abs(DispCA3F), 1 );
madCA3N = mean( abs(DispCA3N), 1 );
DispCAall = {DispCA1N; DispCA1F; DispCA3N; DispCA3F};
MDc = [mdCA1N', mdCA1F', mdCA3N', mdCA3F'];
MADc = [madCA1N', madCA1F', madCA3N', madCA3F'];

figure 
subplot(1,2,1)
    for n = 1:2 %size(MSDc, 2)
        for lap = 1:Minlaps
            mdCA_CI{n}(:,lap) = bootci(1000, @mean, DispCAall{n}(:,lap));
        end
        plot_ci(0:Minlaps-1, [MDc(:,n) mdCA_CI{n}(1,:)' mdCA_CI{n}(2,:)'], 'MainLineColor', colors(n,:), 'PatchColor', colors(n,:), 'LineColor', colors(n,:), 'PatchAlpha', 0.5); hold on
    end
    xlabel('Laps after onset')
    ylabel('Mean Displacement (cm)')
    box off
    axis square
subplot(1,2,2)
    for n = 1:2 %size(MSDc, 2)
        for lap = 1:Minlaps
            madCA_CI{n}(:,lap) = bootci(1000, @mean, abs(DispCAall{n}(:,lap)));
        end
        plot_ci(0:Minlaps-1, [MADc(:,n) madCA_CI{n}(1,:)' madCA_CI{n}(2,:)'], 'MainLineColor', colors(n,:), 'PatchColor', colors(n,:), 'LineColor', colors(n,:), 'PatchAlpha', 0.5); hold on
    end
    xlabel('Laps after onset')
    ylabel('Mean Absolute Displacement (cm)')
    box off
    axis square

% figure 
% errorbar(1:size(Trajs,2), MSD,MSD-msdCI(1,:),msdCI(2,:)-MSD,'-', 'LineWidth', 1, 'Color', cline(4,:)); hold on % BTSP 
% errorbar(1:size(Trajs2,2), MSD2,MSD2-msd2CI(1,:),msd2CI(2,:)-MSD2,'k-', 'LineWidth', 1); hold on % no plasticity
% % errorbar(1:size(Trajs3,2), MSD3,MSD3-msd3CI(1,:),msd3CI(2,:)-MSD3,'r-', 'LineWidth', 1); hold on % no plasticity
% plot(1:Minlaps, msdCA1F,'-', 'LineWidth',1,'Color', mapCA1(1,:)); hold on 
% plot(1:Minlaps, msdCA1N,'-', 'LineWidth',1,'Color', mapCA1(2,:)); hold on 
% plot(1:Minlaps, msdCA3F,'-', 'LineWidth',1,'Color', mapCA3(1,:)); hold on 
% plot(1:Minlaps, msdCA3N,'-', 'LineWidth',1,'Color', mapCA3(2,:)); hold on 
% % xlim([1 size(Trajs,2)+1])
% xlabel('laps')
% ylabel('mean squared displacement (cm^2)')
% box off
% axis square

% figure
% subplot(1,4,1)
%     plot(0:Minlaps-1, msdCA3F,'-', 'LineWidth',1,'Color', mapCA3(1,:)); hold on 
%     plot(0:Minlaps-1, msdCA3N,'-', 'LineWidth',1,'Color', mapCA3(2,:)); hold on 
%     for n = 1:2 %size(MSDc, 2)
%         for lap = 1:Minlaps
%             msdCA_CI{n}(:,lap) = bootci(1000, @mean, DispSqCAall{n}(:,lap));
%         end
%         plot_ci(0:Minlaps-1, [MSDc(:,n) msdCA_CI{n}(1,:)' msdCA_CI{n}(2,:)'], 'MainLineColor', colors(n,:), 'PatchColor', colors(n,:), 'LineColor', colors(n,:), 'PatchAlpha', 0.5); hold on
%     end
%     plot_ci([1:length(MSD)]-1, [MSD' msdCI(1,:)' msdCI(2,:)'], 'MainLineColor', mapBTSP(2,:), 'PatchColor', mapBTSP(2,:), 'LineColor', mapBTSP(2,:), 'PatchAlpha', 0.5); hold on
%     plot_ci([1:length(MSD2)]-1, [MSD2' msd2CI(1,:)' msd2CI(2,:)'], 'MainLineColor', 'k', 'PatchColor', 'k', 'LineColor', 'k', 'PatchAlpha', 0.5); hold on
%     plot_ci([1:length(MSD4)]-1, [MSD4' msd4CI(1,:)' msd4CI(2,:)'], 'MainLineColor', mapBTSP(1,:), 'PatchColor', mapBTSP(1,:), 'LineColor', mapBTSP(1,:), 'PatchAlpha', 0.5); hold on
%     xlabel('Laps after onset')
%     ylabel('Mean Squared Displacement (cm^2)')
%     box off
%     axis square
% subplot(1,4,2)
%     RegStart = 4;
%     Xlm = [ [RegStart:Minlaps]', ones(size([RegStart:Minlaps]'))];
%     for n = 1:size(MSDc,2)
%     Ylm{n} = MSDc(RegStart:end,n);
%     [B{n},BINT{n},R{n},RINT{n},STATS{n}] = regress(Ylm{n}, Xlm);
%     lm{n} = Xlm*B{n};
%     plot([RegStart:Minlaps]'-1, lm{n}, '-', 'LineWidth', 1.5,'Color', colors(n,:)); hold on 
%     scatter([1:Minlaps]-1, MSDc(:,n),'o', 'MarkerFaceColor', colors(n,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', colors(n,:), 'MarkerEdgeAlpha', 0.5); hold on
%     text(Minlaps, max(lm{n}), ['D = ' num2str( round(B{n}(1)/2, 2, 'significant') )], 'Color', colors(n,:))
%     % text(Minlaps+1, max(lm{n})-30, {['slope = ' num2str( round(B{n}(1), 2, 'significant') )]; ['Rsq = ' num2str( round(STATS{n}(1),3, 'significant')*100) '%' ]; ['p = ' num2str( round(STATS{n}(3),3, 'significant')) ]}, 'Color', colors(n,:))
%     end
%     [Bmsd,BINTmsd,Rmsd,RINTmsd,STATSmsd] = regress(MSD(RegStart:Minlaps)', Xlm);
%     lmMSD = Xlm*Bmsd;
%     plot([RegStart:Minlaps]'-1, lmMSD, '-', 'LineWidth', 1.5,'Color', mapBTSP(2,:)); hold on 
%     scatter([1:length(MSD)]-1, MSD, 'o', 'MarkerEdgeColor', mapBTSP(2,:), 'MarkerEdgeAlpha', 0.5, 'MarkerFaceColor', mapBTSP(2,:), 'MarkerFaceAlpha', 0.5); hold on
%     xlim([0 40])
%     text(Minlaps, max(lmMSD), ['D = ' num2str( round(Bmsd(1)/2, 2, 'significant') )], 'Color', mapBTSP(2,:))
%     % text(Minlaps+1, max(lmMSD)-30, {['slope = ' num2str( round(Bmsd(1), 2, 'significant') )]; ['Rsq = ' num2str( round(STATSmsd(1),3, 'significant')*100) '%' ]; ['p = ' num2str( round(STATSmsd(3),3, 'significant')) ]}, 'Color', cline(4,:))
%     [Bmsd4,BINTmsd4,Rmsd4,RINTmsd4,STATSmsd4] = regress(MSD4(RegStart:Minlaps)', Xlm);
%     lmMSD4 = Xlm*Bmsd4;
%     plot([RegStart:Minlaps]'-1, lmMSD4, '-', 'LineWidth', 1.5,'Color', mapBTSP(1,:)); hold on 
%     scatter([1:length(MSD4)]-1, MSD4, 'o', 'MarkerEdgeColor', mapBTSP(1,:), 'MarkerEdgeAlpha', 0.5, 'MarkerFaceColor', mapBTSP(1,:), 'MarkerFaceAlpha', 0.5); hold on
%     xlim([0 40])
%     text(Minlaps, max(lmMSD4), ['D = ' num2str( round(Bmsd4(1)/2, 2, 'significant') )], 'Color', mapBTSP(1,:))
%     xlabel('Laps after onset')
%     ylabel('MSD (cm^2)')
%     box off
%     axis square
% subplot(1,4,3)
%     plot(1:length(MSD)-1, diff(MSD)./2,'-', 'LineWidth',1, 'Color', cline(4,:)); hold on 
%     plot(1:Minlaps-1, diff(msdCA1F)./2,'-', 'LineWidth',1,'Color', mapCA1(1,:)); hold on 
%     plot(1:Minlaps-1, diff(msdCA1N)./2,'-', 'LineWidth',1,'Color', mapCA1(2,:)); hold on 
%     plot(1:Minlaps-1, diff(msdCA3F)./2,'-', 'LineWidth',1,'Color', mapCA3(1,:)); hold on 
%     plot(1:Minlaps-1, diff(msdCA3N)./2,'-', 'LineWidth',1,'Color', mapCA3(2,:)); hold on 
%     xlabel('Laps after onset')
%     ylabel('Instant. Diffusion Coeff. (cm^2/lap)')
%     box off
%     axis square
% subplot(1,4,4)
%     Xeval = 1:Minlaps;
%     Xdiff = 2:Minlaps;
%     Xdiff2 = 2:0.1:Minlaps;
%     modelfun3 = fittype(@(p1,p2,p3, x) p1*exp(-(x-1)/p2)+p3);
%     paramsMSD = [100 2 0];
%     options3 = fitoptions('Method','NonlinearLeastSquares', 'Algorithm', 'Trust-Region');
%     options3.Lower = [0,0.01,0];
%     options3.Upper = [1000,100,20];
%     options3.Startpoint = paramsMSD;
%     for n = 1:size(MSDc,2)
%         Dcoef(:,n) = diff(MSDc(:,n))./2;
%         [Dmdl{n},gofD{n},outD{n}] = fit(Xdiff', Dcoef(:,n), modelfun3, options3);
%         Deval{n} = feval(Dmdl{n}, Xdiff2);
%         scatter(Xdiff-1, Dcoef(:,n),'o','MarkerFaceColor', colors(n,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', colors(n,:), 'MarkerEdgeAlpha', 0.5); hold on
%         plot(Xdiff2-1, Deval{n}, '-', 'Color', colors(n,:), 'LineWidth', 2); hold on
%     end
%     DcoefMSD = diff(MSD(1:Minlaps)')./2;
%     [Dmsd,gofDmsd,outDmsd] = fit(Xdiff', DcoefMSD, modelfun3, options3);
%     DevalMSD = feval(Dmsd, Xdiff2);
%     scatter(Xdiff-1, DcoefMSD,'o','MarkerFaceColor', mapBTSP(2,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', mapBTSP(2,:), 'MarkerEdgeAlpha', 0.5); hold on
%     plot(Xdiff2-1, DevalMSD, '-', 'Color', mapBTSP(2,:), 'LineWidth', 2); hold on
%     xlabel('Laps after onset')
%     ylabel('Instant. Diffusion Coeff. (cm^2/lap)')
%     box off
%     axis square

figure % MSD as lines without CI
    for n = 1:size(MSDc, 2)
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

save('CanData_CA1_MSD', "MSDc", '-mat')

figure % MSD as points, no CI, + lin reg fit to estimate Diff Coef
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

% bootstrapped linear regression on MSD to check if late portion of MSD really has a slope different from 0
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
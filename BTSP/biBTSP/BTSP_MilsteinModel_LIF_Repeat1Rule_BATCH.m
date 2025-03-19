%BTSPmodelLIF_Repeat1Rule_BATCH
clear 
close all
tic
%% Parameters
% Inputs params
Params.N = 100; % number of input neurons
Params.L = 300; % length of track in cm
Params.PFsd = 18; % Standard deviation of mean PF, in cm
Params.PFamp = 10; %peak FR in Hz
Params.Nlaps = 30; % number of laps
Params.period = 20; % lap duration, in sec
Params.dt = 0.001; % time resolution in sec

% Synapses params
Params.Imax = 85e-12; % max synaptic weight, a current in Amps. % NolanWiles2011 uses 10-226pA. Mehta uses 1.5-4nA. SongAbbot use conductances with G between 0 and 150pS, but not sure about the initial weight matrix. 
Params.Idecay = 10e-3; % EPSC time constant, in sec
Params.Wsd = 10; %1.6*Params.PFsd/(Params.L/Params.N); %in cm, standard deviation of initial synaptic weight vector (for gaussian connectivity)
Params.maxW = 1; % from 0 to 1, proportion of Imax defining the initial maximum weight

% Output neuron LIF params
Params.Rm = 100e6; % membrane resistance, in Ohm
Params.tau_m = 20e-3 ; %membrane time constant, in sec
Params.Vrest = -70e-3 ; % resting membrane potential, in Volt
Params.Vthr = -54e-3; % spike threshold, in Volt
Params.Vreset = -60e-3; % reset potential after spike, in Volt

Params.Adapt = 0; % 1 or 0, to implement spike-rate adaptation or not
Params.Ek = -70e-3; % K+ equilibrium potential, in Volt
Params.dSRA = 0.06; % increment value of spike-rate adaptation variable modeling adapting K+ leak current (cf DayanAbbott fig 5.6 and eq 5.13)
Params.tau_sra = 100e-3; % time constant for exponential decay of SRA variable

% Plasticity params
Params.CSproba = 0.2/100; % proba of an output spike to be a plateau potential/complex spike (1.8% on average in Bittner et al. 2017)
    %Milstein-model biBTSTP rule params (average of fits across all Bittner+Milstein recordings)
% Params.Pdecay = 1.66; % Synaptic Eligibility Trace decay time constant, in seconds (average: 0.9; Milstein Fig 5: 2.5s; for Network: 1.66s) 
% Params.Ddecay = 0.44; % Global Instructive Signal (from CS) decay time constant, in seconds (average: 0.5; Milstein Fig 5: 1.5s; for Network: 0.44s) 
% Params.k_pot = 1.1; % in s-1 (average: 2.27; Milstein Fig 5: 1.7; for Network: 1.1) 
% Params.k_dep = 0.7; % in s-1 (average: 0.33; Milstein Fig 5: 0.204; for Network: 0.424) 
% Params.a_pot = 0.415; % sigmoid threshold for potentiation gain (average: 0.24; Milstein Fig 5: 0.5; for Network: 0.415) 
% Params.b_pot = 4.4; % sigmoid slope for potentiation gain (average: 30.32; Milstein Fig 5: 4; for Network: 4.40) 
% Params.a_dep = 0.026; % thresh for depression gain (average: 0.09; Milstein Fig 5: 0.01; for Network: 0.026) 
% Params.b_dep = 20.04; % slope for depression gain (average: 2260; Milstein Fig 5: 44.44; for Network: 20.04) 
% Params.Wmax = 1.5*Imax; % 
Params.Pdecay = 2; % Synaptic Eligibility Trace decay time constant, in seconds (average: 0.9; Milstein Fig 5: 2.5s; for Network: 1.66s) 
Params.Ddecay = 1.5; % Global Instructive Signal (from CS) decay time constant, in seconds (average: 0.5; Milstein Fig 5: 1.5s; for Network: 0.44s) 
Params.k_pot = 1.2; % in s-1 (average: 2.27; Milstein Fig 5: 1.7; for Network: 1.1) 
Params.k_dep = 0.2; % in s-1 (average: 0.33; Milstein Fig 5: 0.204; for Network: 0.424) 
Params.a_pot = 0.5; % sigmoid threshold for potentiation gain (average: 0.24; Milstein Fig 5: 0.5; for Network: 0.415) 
Params.b_pot = 4; % sigmoid slope for potentiation gain (average: 30.32; Milstein Fig 5: 4; for Network: 4.40) 
Params.a_dep = 0.01; % thresh for depression gain (average: 0.09; Milstein Fig 5: 0.01; for Network: 0.026) 
Params.b_dep = 44; % slope for depression gain (average: 2260; Milstein Fig 5: 44.44; for Network: 20.04) 
Params.Wmax = 125e-12; % in Amps. e.g. = 1.4*Params.Imax

% output PF 
Params.Nbin = 50; % Number of bins in which the length of the track is divided

% number of simulations
Nsim = 500;

%% Simulations
for n = 1:Nsim
    output(n) = BTSP_MilsteinModel_LIFadapt(Params);

% % Warning: loading all that in the workspace demands a lot of RAM when Nsim is large 
    % I{n} = output(n).I;
    % V{n} = output(n).V;
    % D{n} = output(n).D;
    % P{n} = output(n).P;
    % W{n} = output(n).W;
    % plateau{n} = output(n).plateau;
%     Spiketimes_in{n} = output.Spiketimes_in;
%     Spiketimes_out{n} = output(n).Spiketimes_out;

    COMslope(n) = output(n).COMslope;
    shiftPval(n) = output(n).shiftPval;
    shiftR2(n) = output(n).shiftR2;
    SD_meanPF(n) = output(n).SD_meanPF;
%     PFsdDiff(n) = output(n).PFsdOut(end) - output(n).PFsdOut(1);
    PFsdDiff(n) = mean(output(n).PFsdOut(Params.Nlaps-2:end)) - mean(output(n).PFsdOut(1:3));
    SpatialFRout{n} = output(n).PF;
    COMbin{n} = output(n).COMbin;
    CSbin{n} = output(n).CSbin;
end
    
% group as forward, backward and nonsignificant shifting
IdxFwd = find(shiftPval<=0.05 & COMslope>0);
IdxBck = find(shiftPval<=0.05 & COMslope<0);
IdxNonSig = find(shiftPval>0.05);

PropFwd = length(IdxFwd)/Nsim;
PropBck = length(IdxBck)/Nsim;
PropNonSig = length(IdxNonSig)/Nsim;

toc
%% figures
blue = [0 0.45 0.74]; % backward
red = [0.85 0.32 0.01]; % forward
grey = [0.5 0.5 0.5]; %nonsig

% figure(1) % rule   
%     delay1 = -5000:1:0; % ms
%     delay2 = 0:1:5000; %ms
%     Wd = (Params.Pamp)*exp(-delay2/(Params.Ddecay*1000));
%     Wp = (Params.Pamp)*exp(delay1/(Params.Pdecay*1000));
%     plot([delay1 delay2]./1000, [Wp*10^12 Wd*10^12], 'k'); hold on
%     xline(0, 'k--'); hold on
%     yline(0, 'k--'); 
%     xlim([-4.2 4.2])
%     ylim([0 22])
%     xlabel('delay (s)'); ylabel('Weight change (pA)');
%     title(['Adapt: ' num2str(Params.Adapt), ', Shunting: ' num2str(Params.Shunting), ', capWeights: ' num2str(Params.capWeights) ', Homeo: ' num2str(Params.HomeoNorm) ', WUdyn: ' num2str(Params.WUdyn)])
%     box off; axis square;

figure(2) % R2 vs slope
    scatter(COMslope(IdxFwd), shiftR2(IdxFwd), [], red); hold on
    scatter(COMslope(IdxBck), shiftR2(IdxBck), [], blue); hold on
    scatter(COMslope(IdxNonSig), shiftR2(IdxNonSig), [], grey);
    % xlim([-2, 1.5])
    xlabel('Shift (COM slope), cm/lap');
    ylabel('R-square');
    box off; axis square;
% add bootstrap of mean or median to fig 2

figure(3) % violin plot of slopes, with mean and median + line at 0 + 95% bootstrap CI
    plot([0.5 1.5], [0 0], 'k--') 
    %yline(0,'k--'); hold on
    violinplot(COMslope, [], 'Width', 0.3, 'ViolinColor', [0, 0, 0], 'ShowMean', true(1));
    ylabel('Shift (COM slope), cm/lap');
    box off; axis square;

figure(4) % slopes, zoom
    %yline(0,'k--'); hold on
    plot([0.5 1.5], [0 0], 'k--') 
    violinplot(COMslope, [], 'Width', 0.3, 'ViolinColor', [0, 0, 0], 'ShowMean', true(1));
    ylabel('Shift (COM slope), cm/lap');
    ylim([-1, 1])
    box off; axis square;

figure(5) % proportion of forward, backward and stable/non-signif PFs
    b = bar(1, [PropBck PropFwd PropNonSig], 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylabel('proportion');
    box off; axis square;

figure(6)
    plot([0.5 1.5], [0 0], 'k--') 
    %yline(0,'k--'); hold on
    violinplot(PFsdDiff, [], 'Width', 0.3, 'ViolinColor', [0, 0, 0], 'ShowMean', true(1));
%     ylim([-1 11])
    ylabel('PF \DeltaS.D., cm (last 3 laps - 1st 3 laps)');
    box off; axis square;

figure(7) % Width increase vs slope
    % PFsdDiff_mean = mean(PFsdDiff);
    % PFsdDiff_CI = bootci(1000, @mean, PFsdDiff);
    % PFsdDiff_CIlo = PFsdDiff_CI(1,1);
    % PFsdDiff_CIup = PFsdDiff_CI(2,1);
    % y = PFsdDiff_mean; yneg = y-PFsdDiff_CIlo; ypos = PFsdDiff_CIup-y;
    plot([-3 1.5], [0 0], 'k--') ; hold on
    scatter(COMslope(IdxFwd), PFsdDiff(IdxFwd), [], red, 'filled', 'MarkerFaceAlpha', 0.5); hold on
    scatter(COMslope(IdxBck), PFsdDiff(IdxBck), [], blue, 'filled', 'MarkerFaceAlpha', 0.5); hold on
    scatter(COMslope(IdxNonSig), PFsdDiff(IdxNonSig), [], grey, 'filled', 'MarkerFaceAlpha', 0.4); hold on
    % errorbar(0, y,yneg,ypos,'o', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 5, 'LineWidth', 2); hold on
    xlabel('Shift (COM slope), cm/lap');
    ylabel('PF \DeltaS.D., cm (last 3 laps - 1st 3 laps)');
    box off; axis square;

%% examples

% find extrema examples
[minShift, idxMin] = min(COMslope);
[maxShift, idxMax] = max(COMslope);
[minPFsD, idxMin2] = min(PFsdDiff);
[maxPFsdD, idxMax2] = max(PFsdDiff);

% choose Neuron to plot
Ni = idxMin;

% recompute unsaved time and spatial axes for some figures
speed = Params.L/Params.period; % cm/sec
EndTime = Params.Nlaps*Params.L/speed; %in sec
Trun = 0:Params.dt:EndTime; Trun = Trun(1:end-1); % ms resolution
Tlap1 = 0:Params.dt:Params.period; Tlap1 = Tlap1(1:end-1);
Lap1 = linspace(0,Params.L,Params.period/Params.dt); % positions during 1 lap
Run = repmat(Lap1,[1 Params.Nlaps]);%positions during the whole run
NewLapIdx = find(Run==0);
NewLapTimes = Trun(NewLapIdx);

figure %FR map of extreme backward and forward shifts
    [row_lap,column_bin] = find(CSbin{Ni});
    imagesc(SpatialFRout{Ni}); hold on
    scatter(COMbin{Ni},1:Params.Nlaps, 30, [1 0 0], 'filled'); hold on
    scatter(column_bin, row_lap, 10, [0 1 1], 'filled');
    xlabel('spatial bins'); ylabel('lap');
    legend('Complex spike','COM', 'Location', 'BestOutside');
    colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
    c = colorbar; c.Label.String = 'Firing rate (Hz)';
    box off; axis square;

figure %spatial FR lap 1 vs last lap
    subplot(1,2,1)
        plot(SpatialFRout{Ni}(1,:), 'k'); hold on
        plot(SpatialFRout{Ni}(end,:), 'r'); hold on
        xlabel('spatial bins'); ylabel('Firing rate (Hz)');
        box off; axis square;
    subplot(1,2,2)
        plot(SpatialFRout{Ni}(1,:)./max(SpatialFRout{Ni}(1,:)), 'k'); hold on
        plot(SpatialFRout{Ni}(end,:)./max(SpatialFRout{Ni}(end,:)), 'r'); hold on
        xlabel('spatial bins'); ylabel('norm FR relative to peak');
        box off; axis square;

% [~, idxR2] = maxk(shiftR2, 20);
% [~, idxPval] = maxk(shiftPval, 20);
% gpidx = find(shiftPval<0.05 & shiftR2>0.9);
% 
% Ni = gpidx; 
% for n = 1:length(Ni)
%     m = Ni(n);
% figure %FR map of extreme backward and forward shifts
%     [row_lap,column_bin] = find(CSbin{m});
%     imagesc(SpatialFRout{m}); hold on
%     scatter(COMbin{m},1:Params.Nlaps, 30, [1 0 0], 'filled'); hold on
%     scatter(column_bin, row_lap, 10, [0 1 1], 'filled');
%     xlabel('spatial bins'); ylabel('lap');
%     legend('Complex spike','COM', 'Location', 'BestOutside');
%     colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
%     c = colorbar; c.Label.String = 'Firing rate (Hz)';
%     box off; axis square;
% end

% Ni = 1;
% figure %BTSP variables on given lap for initial peak input
% LapNum = row_lap;
% x1 = 0; % in s
% x2 = Params.period; % in s 
% % CSidx = find(CS_laps{LapNum});
% InputNum = Params.N/2; %dWmaxIn; %Nbittner
% maxD = max(D{Ni}(NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1));
% maxP = max(P{Ni}(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1));
% normD = D{Ni}(NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./maxD;
% normP = P{Ni}(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./maxP;
% plot(Tlap1, normD, 'r'); hold on
% plot(Tlap1, normP, 'k'); hold on
% % scatter(Tlap1(CSidx), ones(size(CSidx)).*0.5, 'oc', 'filled')
% legend('normIS', 'normET')
% xlim([x1 x2]);
% xlabel('time (s)');
% title(['input ' num2str(InputNum)])
% box off;
% 
% figure % Weight changes over time, full run
% imagesc([0 Params.Nlaps*Params.period],[1 Params.N], W{Ni}.*10^12); hold on
% scatter(NewLapTimes, (Params.N/2)*ones(1,Params.Nlaps), 10, [1 0 0], '+');
% xlabel('time (s)'); ylabel('input neuron');
% colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
% c = colorbar; c.Label.String = 'Synaptic weight (pA)';
% box off;
%% mean square displacement (random walk analysis)

Trajs = cat(1,COMbin{:}).*6; % in cm
Disp = Trajs - Trajs(:,1); % trajectories centered on initial location, with forward direction as positive values
DispSq = Disp.^2;
MSD = mean(DispSq,1,"omitmissing"); % average all neurons for each lap
rootMSD = sqrt(MSD);
for lap = 1:size(Trajs,2)
msdCI(:,lap) = bootci(1000, @nanmean, DispSq(:,lap));
end
[B,BINT,R,Rint,STATS] = regress(MSD', [ones(size(Trajs,2),1), [1:size(Trajs,2)]'] );
linMDL = B(1)+B(2)*[1:size(Trajs,2)];

figure
errorbar(1:size(Trajs,2), MSD,MSD-msdCI(1,:),msdCI(2,:)-MSD,'o-', 'Color', 'k', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 3, 'LineWidth', 1.5); hold on
plot([1 size(Trajs,2)], [0 MSD(end)], 'g'); hold on
plot(1:size(Trajs,2), linMDL, 'r'); hold on
xlim([1 size(Trajs,2)+1])
xlabel('laps after emergence')
ylabel('mean squared displacement (cm^2)')
box off
axis square

figure
plot(1:size(Trajs,2), rootMSD,'k-')
xlim([1 size(Trajs,2)+1])
ylim([1 size(Trajs,2)+1])
xlabel('laps after emergence')
ylabel('root mean squared displacement (cm)')
box off
axis square

figure
plot(1:size(Trajs,2), MSD,'k-')
xlim([1 size(Trajs,2)+1])
ylim([1 size(Trajs,2)+1])
xlabel('laps after emergence')
ylabel('root mean squared displacement (cm)')
box off
axis square
%% figure to adapt to present code
% figure
%     subplot(2,2,1)
%         plot(1:Nlaps, meanFRout_lap, 'r+'); hold on
%         plot(1:Nlaps, maxFRout_lap, 'k-');
%         xlabel('lap'); ylabel('Firing Rate (Hz)')
%         legend('mean','max', 'Location', 'BestOutside');
%         box off; axis square;
%     subplot(2,2,2)
%         plot(1:Nlaps, PFsd, 'k-');
%         xlabel('lap'); ylabel('PF sd (cm)');
%         box off; axis square;
%     subplot(2,2,3)
%         plot(1:Nlaps, PFskew, 'k-');
%         xlabel('lap'); ylabel('PF skewness');
%         box off; axis square;
%     subplot(2,2,4)
%         plot(BinCenters, meanFRmap, 'k-'); hold on
%         xline(COM_meanPF, 'r'); hold on
%         xline(L/2, 'k--'); 
%         xlabel('position'); ylabel('mean Firing Rate (Hz) ');
%         legend('mean PF', 'COM', 'track center', 'Location', 'BestOutside');
%         title({'SD = ' num2str(SD_meanPF) ', Skew = ' num2str(Skew_meanPF)})
%         box off; axis square;
% 
% figure % Shift: COM regression
% plot(1:Nlaps, COMloc, 'k-'); hold on
% plot(1:Nlaps, COMtraj, 'b-'); hold on
% yline(COMloc(1),'r');
% legend('COM','lin reg','COM #1')
% ylim([0 300]);
% xlabel('lap'); ylabel('COM position (cm)');
% title({'Place Field COM trajectory, slope =' num2str(b(2)) 'pval =' num2str(stats(3))})
% box off; axis square;

% %% figure to plot when only 1 simul. Comment out otherwise.
% figure % Plasticity variables for a given lap and input neuron
% 
% LapNum = 1;
% InputNum = 66;
% x1 = 0; % in s
% x2 = 20; % in s 
% 
% Spiketimes_in_OK = Spiketimes_in{Ni}{InputNum}(Spiketimes_in{Ni}{InputNum} < NewLapTimes(LapNum+1) & Spiketimes_in{Ni}{InputNum} >= NewLapTimes(LapNum));
% 
%     subplot(3,1,1)
%         plot(Tlap1, D{Ni}(NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'r'); hold on
%         plot(repmat(Spiketimes_out{Ni}{LapNum},2,1), repmat([0; 1], 1, length(Spiketimes_out{Ni}{LapNum})), 'k');
%         xlim([x1 x2]);
%         xlabel('time (s)'); 
%         title('output raster + Post-before-Pre variable')
%         box off;
%     subplot(3,1,2)
%         plot(Tlap1, P{Ni}(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'r'); hold on
%         plot( repmat( Spiketimes_in_OK, 2, 1), repmat( [0; 1], 1, length(Spiketimes_in_OK) ), 'k' ) ;
%         xlim([x1 x2]);
%         xlabel('time (s)'); 
%         title(['input raster #' num2str(InputNum) '+ Pre-before-Post variable'])
%         box off;
%     subplot(3,1,3)
%         plot(Tlap1, W{Ni}(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./Params.Imax, 'r'); hold on
%         if Params.WUdyn == 1
%         plot(Tlap1, Wtarget{Ni}(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./Params.Imax, 'k'); hold on
%         end
%         xlim([x1 x2]);
%         xlabel('time (s)'); 
%         ylabel('Norm. synaptic strength')
%         title(['input #' num2str(InputNum)] )
%         box off;
% 
% figure % Weight changes over time, full run
% % subplot(2,1,1)
% imagesc([0 Params.Nlaps*Params.period],[1 Params.N], W{Ni}.*10^12); hold on
% scatter(NewLapTimes, (Params.N/2)*ones(1,Params.Nlaps), 10, [1 0 0], '+');
% xlabel('time (s)'); ylabel('input neuron');
% colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
% c = colorbar; c.Label.String = 'Synaptic weight (pA)';
% box off; 
% axis square;
% 
% figure % initial weights vs end weights
%     subplot(1,2,1)
%         plot(1:Params.N, W{Ni}(:,1).*10^12, 'k'); hold on
%         plot(1:Params.N, W{Ni}(:,end).*10^12, 'r'); hold on
%         legend('start', 'end', 'Location', 'Best');
%         xlabel('input neuron'); ylabel('synaptic weight (pA)')
%         box off; 
%         axis square;
%     subplot(1,2,2)
%         plot(1:Params.N, W{Ni}(:,1)./max(W{Ni}(:,1)), 'k'); hold on
%         plot(1:Params.N, W{Ni}(:,end)./max(W{Ni}(:,end)), 'r'); hold on
%         legend('start', 'end', 'Location', 'Best');
%         xlabel('input neuron'); ylabel('synaptic weight relative to peak')
%         box off; 
%         axis square;
% 
% 
% 

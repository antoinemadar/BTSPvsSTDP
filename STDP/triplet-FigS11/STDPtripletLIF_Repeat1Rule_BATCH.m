%STDPmodelLIF_Repeat1Rule_BATCH
clear 
close all

% for dynamic inputs: slope distribution
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA3_SlopeDistrib.mat')

tic
%% Parameters
% Inputs params
Params.N = 100; % number of input neurons
Params.L = 300; % length of track in cm
Params.PFsd = 18; % Standard deviation of mean PF, in cm
Params.PFamp = 15; %peak FR in Hz
Params.Nlaps = 30; % number of laps
Params.period = 20; % lap duration, in sec
Params.dt = 0.001; % time resolution in sec
Params.InShift = 0; % 0 if static input, 1 if shift
Params.PslopeBinCenters = [-1.6 -1.4, -1.2 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6]; %in cm/lap
Context = 'N'; % N for novel, F for Familiar, NandF for all

if ismember('F', Context)
    Params.Pslope = CA3_F_SlopeHist.Values;
elseif ismember('N', Context)
    Params.Pslope = CA3_N_SlopeHist.Values;
elseif ismember('FandN', Context)
    Params.Pslope = CA3_FandN_SlopeHist.Values;
end

% Synapses params
Params.Imax = 85e-12; % max synaptic weight, a current in Amps. 
Params.Idecay = 10e-3; % EPSC time constant, in sec
Params.Wsd = 10; %1.6*Params.PFsd/(Params.L/Params.N); %in neurons, standard deviation of initial synaptic weight vector (for gaussian connectivity)
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

Params.Shunting = 0; % implement shunting inhib (capping I) if = 1, 0 otherwise. 
Params.Icap = 350e-12; % in Amps

% Plasticity params

%ZenkeGerstner2015
Params.Pdecay = 20e-3; % Pre-before-Post time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Params.Ddecay = 20e-3; % Post-before-Pre time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Params.Pdecay2 = 100e-3; % tau_x from All-to-All model in Table 4 of Pfistergerstner2006
Params.Ddecay2 = 100e-3; % tau_y
Params.Pamp = 0.5/100; % A2+ in PfisterGerstner2006. Peak Pre-before-Post weight change for paired protocol (0.5% in SongAbbot2000, 0.6pA with no Imax in MehtaWilson2000). 
Params.Damp = 0.5/100; %0.16/100; % A2- in PfisterGerstner2006. Peak Post-before-Pre weight change for paired protocol (-0.525% in SongAbbott2000, -90% of Pamp in MehtaWilson2000, which didn't have a maximum weight)
Params.Pamp2 = 0; % A3- in PfisterGerstner2006
Params.Damp2 = 0.5/100; % A3+ in PfisterGerstner2006

%WangBi2005 and PfisterGerstner2006
% Params.Pdecay = 16.8e-3; % Pre-before-Post time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
% Params.Ddecay = 33.7e-3; % Post-before-Pre time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
% Params.Pdecay2 = 946e-3; % tau_x from All-to-All model in Table 4 of Pfistergerstner2006
% Params.Ddecay2 = 27e-3; % tau_y
% Params.Pamp = 0.61/100; % A2+ in PfisterGerstner2006. Peak Pre-before-Post weight change for paired protocol (0.5% in SongAbbot2000, 0.6pA with no Imax in MehtaWilson2000). 
% Params.Damp = 0.16/100; % A2- in PfisterGerstner2006. Peak Post-before-Pre weight change for paired protocol (-0.525% in SongAbbott2000, -90% of Pamp in MehtaWilson2000, which didn't have a maximum weight)
% Params.Pamp2 = 0.14/100; % A3- in PfisterGerstner2006
% Params.Damp2 = 0.67/100; % A3+

Params.capWeights = 1; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 

Params.HomeoNorm = 0; % 0 if no homeostatic normalization rule, 1 if yes.

Params.WUdyn = 0; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Params.Wrest = 0; % Synaptic weight resting state, in Amps (towards which any synapse would relax if Wtarget = 0) 
Params.tau_wu = 5; % weight update time constant, in sec

% output PF 
Params.Nbin = 50; % Number of bins in which the length of the track is divided. 50 for CanSheffield2021 (6 cm), 120 for MouJi2018 (2.5 cm)

% number of simulations
Nsim = 100;

%% Simulations

for n = 1:Nsim
    output(n) = STDPtriplet(Params);
    % output(n) = STDPplus_LIFadapt(Params);
    
% % Warning: loading all that in the workspace demands a lot of RAM. 
%     I{n} = output(n).I;
%     V{n} = output(n).V;
%     D{n} = output(n).D;
%     P{n} = output(n).P;
%     W{n} = output(n).W;
%     Wtarget{n} = output(n).Wtarget;
% 
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
    meanFRmap(n,:) = output(n).meanFRmap;
    PFoutPeak(n) = max(meanFRmap(n,:));
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

figure(1) %STDP rule   
    delay1 = -100:1:0; % ms
    delay2 = 0:1:100; %ms
    Wd = (Params.Damp*100)*exp(delay1/(Params.Ddecay*1000));
    Wp = (Params.Pamp*100)*exp(-delay2/(Params.Pdecay*1000));
    plot([delay1 delay2], [Wd Wp], 'k'); hold on
    xline(0, 'k--'); hold on
    yline(0, 'k--'); 
    xlabel('delay (ms)'); ylabel('Weight change (% of EPSC_{max})');
    title(['Adapt: ' num2str(Params.Adapt), ', Shunting: ' num2str(Params.Shunting), ', capWeights: ' num2str(Params.capWeights) ', Homeo: ' num2str(Params.HomeoNorm) ', WUdyn: ' num2str(Params.WUdyn)])
    box off; axis square;

figure(2) % R2 vs slope
    scatter(COMslope(IdxFwd), shiftR2(IdxFwd), 50, red); hold on
    scatter(COMslope(IdxBck), shiftR2(IdxBck), 50, blue); hold on
    scatter(COMslope(IdxNonSig), shiftR2(IdxNonSig), 50, grey);
    xlabel('Shift (COM slope), cm/lap');
    xlim([-4 3]); ylim([0 1]);
    ylabel('R-square');
    box off; axis square;
% add bootstrap of mean or median to fig 2

    figure(3) % R2 vs slope zoom
    scatter(COMslope(IdxFwd), shiftR2(IdxFwd), 100, red, 'filled', 'MarkerFaceAlpha',0.6); hold on
    scatter(COMslope(IdxBck), shiftR2(IdxBck), 100, blue, 'filled', 'MarkerFaceAlpha',0.6); hold on
    scatter(COMslope(IdxNonSig), shiftR2(IdxNonSig), 100, grey, 'filled', 'MarkerFaceAlpha',0.6);
    xlabel('Shift (COM slope), cm/lap');
    ylabel('R-square');
    box off; axis square;

% figure % violin plot of slopes, with mean and median + line at 0 + 95% bootstrap CI
%     plot([0.5 1.5], [0 0], 'k--') 
%     %yline(0,'k--'); hold on
%     violinplot(COMslope, [], 'Width', 0.3, 'ViolinColor', [0, 0, 0], 'ShowMean', true(1));
%     ylabel('Shift (COM slope), cm/lap');
%     box off; axis square;

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
clear PFsdDiff
for n = 1:Nsim
    PFsdDiff(n) = mean(output(n).PFsdOut(Params.Nlaps-2:end)) - mean(output(n).PFsdOut(1:3));
end
    plot([0.5 1.5], [0 0], 'k--') 
    %yline(0,'k--'); hold on
    violinplot(PFsdDiff, [], 'Width', 0.3, 'ViolinColor', [0, 0, 0], 'ShowMean', true(1));
    ylim([-3 14])
    ylabel('PF \DeltaS.D., cm (last 3 laps - 1st 3 laps)');
    box off; axis square;

figure(7) % Width increase vs slope
    scatter(abs(COMslope(IdxFwd)), PFsdDiff(IdxFwd), 'r'); hold on
    scatter(abs(COMslope(IdxBck)), PFsdDiff(IdxBck), 'b'); hold on
    scatter(abs(COMslope(IdxNonSig)), PFsdDiff(IdxNonSig), 'k');
%     ylim([0 11])
    xlabel('Absolute Shift, cm/lap');
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
    imagesc(SpatialFRout{Ni}); hold on
    scatter(COMbin{Ni},1:Params.Nlaps, 30, [1 0 0], 'filled');
    xlabel('spatial bins'); ylabel('lap');
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

tl = tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
title(tl, ['STDP dWmax = ' num2str(Params.Pamp*100) '%, STDP tau = ' num2str(Params.Pdecay.*10^3) ' ms, PFin peak = ' num2str(Params.PFamp) ' Hz, Speed = ' num2str(Params.L/Params.period) ' cm/s'])

nexttile
    scatter(COMslope(IdxFwd), shiftR2(IdxFwd), 30, red); hold on
    scatter(COMslope(IdxBck), shiftR2(IdxBck), 30, blue); hold on
    scatter(COMslope(IdxNonSig), shiftR2(IdxNonSig), 50, grey);
    xlabel('Shift (COM slope), cm/lap');
    xlim([-4 3]); ylim([0 1]);
    ylabel('R-square');
    box off; axis square;
nexttile
    clear PFsdDiff
    for n = 1:Nsim
        PFsdDiff(n) = mean(output(n).PFsdOut(Params.Nlaps-2:end)) - mean(output(n).PFsdOut(1:3));
    end
    plot([0.5 1.5], [0 0], 'k--') 
    %yline(0,'k--'); hold on
    violinplot(PFsdDiff, [], 'Width', 0.3, 'ViolinColor', [0, 0, 0], 'ShowMean', true(1));
    % ylim([-30 1])
    ylabel('PF \DeltaSD (cm)'); %(last 3 laps - first 3 laps)
    box off; axis square;
nexttile
    imagesc(SpatialFRout{Ni}); hold on
    scatter(COMbin{Ni},1:Params.Nlaps, 20, [1 0 0], 'filled');
    xlabel('spatial bins'); ylabel('lap');
    colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
    c = colorbar; c.Label.String = 'Firing rate (Hz)';
    box off; axis square;
    title('Largest backward slope')
nexttile
    plot(mean(SpatialFRout{Ni}(1:3,:)), 'k'); hold on
    plot(mean(SpatialFRout{Ni}(end-2:end,:)), 'r'); hold on
    legend('laps 1-3', 'last 3 laps', 'Location','best')
    xlabel('spatial bins'); ylabel('Firing rate (Hz)');
    box off; axis square;

tl2 = tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
title(tl2, ['STDP dWmax = ' num2str(Params.Pamp*100) '%, PFin peak = ' num2str(Params.PFamp) ' Hz'])

% nexttile % ex of lapwise PF
%     imagesc(SpatialFRout{Ni}); hold on
%     scatter(COMbin{Ni},1:Params.Nlaps, 30, [1 0 0], 'filled');
%     ylabel('Lap #');
%     xlabel('Position (cm)')
%     set(gca, 'XTick', [0 Params.Nbin], 'XTickLabel', ['0' '300'])
%     colormap(flipud(gray(256)));
%     c = colorbar; c.Label.String = 'Firing rate (Hz)';
%     box off; axis square;
nexttile % ex average rate map
    plot(mean(SpatialFRout{Ni},1), 'k'); hold on
    ylabel('Firing rate (Hz)');
    xlabel('Position (cm)')
    xlim([0 Params.Nbin])
    set(gca, 'XTick', [0, Params.Nbin], 'XTickLabel', {'0', '300'})
    box off; axis square;
    title('example PF (30 lap average)')
nexttile % average rate map across all simulations
    plot(mean(meanFRmap,1), 'k'); hold on
    ylabel('Firing rate (Hz)');
    xlabel('Position (cm)')
    xlim([0 Params.Nbin])
    set(gca, 'XTick', [0, Params.Nbin], 'XTickLabel', {'0', '300'})
    box off; axis square;    
    title('average rate map across all PFs')
nexttile % PF sd for all simulations
    histogram(SD_meanPF, 'Normalization', 'probability', 'LineWidth', 2, 'DisplayStyle','stairs'); hold on
    xline(mean(SD_meanPF),'k-'); hold on
    xlabel('PF sd (cm)'); ylabel('probability');
    box off; axis square;
    title(['mean = ' num2str(mean(SD_meanPF)) 'cm'])
nexttile % PFout peak FR (of the average PF) for all simulations
    histogram(PFoutPeak, 'Normalization', 'probability', 'LineWidth', 2, 'DisplayStyle','stairs'); hold on
    xline(mean(PFoutPeak),'k-'); hold on
    xlabel('Peak PFout (Hz)'); ylabel('probability');
    title(['mean = ' num2str(mean(PFoutPeak)) 'Hz'])
    box off; axis square;


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


% figure %STDP variables for a given lap and input neuron
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

% figure % Weight changes over time, full run
% % subplot(2,1,1)
% imagesc([0 Params.Nlaps*Params.period],[1 Params.N], W{Ni}.*10^12); hold on
% scatter(NewLapTimes, (Params.N/2)*ones(1,Params.Nlaps), 10, [1 0 0], '+');
% xlabel('time (s)'); ylabel('input neuron');
% colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
% c = colorbar; c.Label.String = 'Synaptic weight (pA)';
% box off; 
% axis square;

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




clear; close all

%% Parameters
%Params for inputs
N = 100; % number of input neurons
L = 185 % length of track in cm (300 for Can, 185 for Bittner and Milstein). For Milstein experiments, I use 2*185 to make sure Vm ramps are not bleeding into next lap + see what happens far from CS, but still use same smoothing parameters
PFsd = 18; % cm (~median in Can data) 
PFamp = 10; %peak FR in Hz
Nlaps = 23 % 30; % number of laps
speed = 25 %15; % cm/sec, (15cm/s average in DongSheffield2021, 25cm/s in Milstein2021 network model)
dt = 0.001; % time resolution in sec

% Synapses params
Imax = 85e-12 % max synaptic weight, a current in Amps. % NolanWiles2011 uses 10-226pA. Mehta uses 1.5-4nA. SongAbbot use dimensionless conductances with G between 0 and Gmax, but not sure about the initial weight matrix. 
Params.Idecay = 10e-3; % EPSC time constant, in sec
Wsd = 10; %1.6*PFsd/(L/N); %standard deviation of initial synaptic weight vector (for gaussian connectivity)
maxW = 1; % from 0 to 1, proportion of Imax defining the initial maximum weight

%Integrate and Fire params for output neuron
Params.Rm = 100e6; % membrane resistance, in Ohm (100MOhm in SongAbbot2000)
Params.tau_m = 20e-3 ; %membrane time constant, in sec (25ms in YuShouval2006, 20ms in SongAbbot2000)
Params.Vrest = -70e-3 ; % resting membrane potential, in Volt (-70mV in SongAbbot2000, -60mV in YuShouval2006)
Params.Vthr = -54e-3; % spike threshold, in Volt (-54mV in SongAbbot2000)
Params.Vreset = -60e-3; % reset potential after spike, in Volt (-60mV in SongAbbot2000) 

Params.Adapt = 0; % 1 or 0, to implement spike-rate adaptation or not
Params.Ek = -70e-3; % K+ equilibrium potential, in Volt
Params.dSRA = 0.06; % increment value of spike-rate adaptation variable modeling adapting K+ leak current (cf DayanAbbott fig 5.6 and eq 5.13)
Params.tau_sra = 100e-3; % time constant for exponential decay of SRA variable

% Plasticity params
Params.Pdecay = 1.31; % Pre-before-Post time constant, in sec (1.31s on average in Bittner et al. 2017)
Params.Ddecay = 0.69; % Post-before-Pre time constant, in sec (0.69s on average in Bittner et al. 2017)
Params.Pamp = 20e-12; % peak weight change, in Amps (3*baseline EPSP for pairings with 5*10 EPSPs at 10Hz, with initial EPSP at 2mV, in Bittner et al. 2017)
Params.PostPreBoost = 1.1; % scaling factor, as a function of Pamp, for the postpre variable (triggered on CS, not on spikes, so with little to no opportunities for temporal summation). Optimized on Bittner2017 invitro data
Params.capWeights = 0; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 
Params.HomeoNorm = 1; % 0 if no homeostatic normalization rule, 1 if yes.
Params.WUdyn = 0; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Params.Wrest = 0; % Synaptic weight resting state, in Amps (towards which any synapse would relax if Wtarget = 0) 
Params.tau_wu = 5; % weight update time constant, in sec

Params.CSproba = []; % proba of an output spike to be a complex spike (make it an empty array when you want CSs in specific locations)
% params for adding CSs at specific times:
addCS = 1; % 0 if no predefined CS, 1 if you want to simulate CS induction
CSlap = 11:13; % lap on which CS will occur
CStime = 3; %-L./(2.*speed); % time before COM of initial PF, in sec. If defined as empty array, it puts CS at a time point corresponding to input 46 COM (closest to 77pA)

% output PF 
Nbin = 50; % Number of bins in which the length of the track is divided (50 bins of 6cm in DongSheffield2021, 100 bins of 1.85cm in Milstein2021)
Nbin2 = round(L/1.85); % for spatial smoothing like in Milstein2021

%% Inputs
%Place fields (FR = gaussian function of position x)
x = 1:L; % track positions binned by cm
PFcom = 1:L/N:L; % PF centers spread uniformly along the track, with L/N spacing in cm
FRin = zeros(N,length(x));
for i = 1:N
    FRinNorm(i,x) = exp(-0.5*(x-PFcom(i)).^2/PFsd^2); %normalized from 0 to 1
    FRin(i,x) = PFamp*exp(-0.5*(x-PFcom(i)).^2/PFsd^2); % real FR 
end

%Running model with constant speed
period = L./speed; % in sec
EndTime = Nlaps*period; %in sec
Trun = 0:dt:EndTime; Trun = Trun(1:end-1); % ms resolution
Tlap1 = 0:dt:period; Tlap1 = Tlap1(1:end-1);
Lap1 = linspace(0,L,period/dt); % positions during 1 lap
Run = repmat(Lap1,[1 Nlaps]);%positions during the whole run

%convert trajectory into firing rate
for i = 1:length(Run)
    [minval, idx] = min(abs(x-Run(i))); %find closest spatial bin in x for each position in trajectory
    FRrun(:,i) = FRin(:,idx); % Firing rate at that position for each neuron
end

% nonhomogenous Poisson spike generator along trajectory
InRasters = poissrnd(FRrun.*dt); %matrix of size (N,length(Tlap1)) containing, for each input neuron, the number of spikes in each 1ms bin drawn from Poisson distribution
for i = 1:N
    Spiketimes_in{i} = Trun( find(InRasters(i,:)) ); %find positive values = when a spike occured
end

% gaussian connectivity with max in center of track, akin to Mehta and Wilson 2000
NeuronID = 1:N;
W = zeros(N,length(Trun));
Wmax = maxW*Imax;
W(:,1) = Wmax*exp(-0.5*(NeuronID - N/2).^2/Wsd^2); % vector of N rows following a gaussian distribution of Wmax amplitude 

%% input-output dynamics

% CS trace: 
CS = zeros(1,length(Trun)); 
if addCS == 1
%add 1 CS 3s before initial COM of current PF (which is in the middle)
    COMidx = length(Tlap1)./2;
    % [Dmin, COMidx] = min(abs(Lap1-(L./2))); % find middle of trajectory for 1st lap
    COM1idx_Laps = find(Run==Run(COMidx)); % initial COM idx on all laps (assuming constant speed)
    CS(COM1idx_Laps(CSlap)-(CStime./dt)) = 1;
end

% run LIF + plasticity model
Params.dt = dt; Params.Imax = Imax;
[I, V, V2, CS, OutRaster, D, P, Dact, Pact, SRA, W, W2, Wtarget]  = BTSPplus_LIFadapt_ODEs(Params, W, InRasters, CS);

Spiketimes_outAll = Trun(logical(OutRaster));
SpikeLoc_outAll = Run(logical(OutRaster)); %spikes location on track
CSloc_All = Run(logical(CS)); % CSs location on track 

% segment in laps
NewLapIdx = find(Run==0);
NewLapTimes = Trun(NewLapIdx);

% smooth Vm traces by using a 3Hz low-pass filter like in Bittner2017
% wrap-around padding of Hamming window size and zero-phase filtering
Fcut = 3; %cut-off frequency, in Hz
Fs = 1./dt; %sampling frequency in Hz
NyqF = Fs/2; % Nyquist frequency
HWin = 0.2; % size of hamming window, in seconds
fir_n = HWin.*Fs-1; % filter order
nFcut = Fcut/NyqF; % normalized cutoff frequency (proprtion of the Nyquist frequency)
blo = fir1(fir_n,nFcut); % design a FIR filter of fir_n'th order, using hamming window, with cut-off frequency at 3Hz
% freqz(blo,1,2^16,Fs)

Vpad = [V2(length(V2)-HWin*Fs:end), V2, V2(1:HWin*Fs)]; % pad V2 on the left to avoid weird thing in beginning
Vlo = filtfilt(blo, 1, Vpad);
Vlo = Vlo(HWin*Fs+1:end-HWin*Fs-1); % remove pad sections =>should be same length as V2

% compute spatial firing rate and related PF properties
Track_dbin = L./Nbin; % bin size, in cm
Track_dbin2 = L./Nbin2; %for spatial smoothing like Milstein2021
TrackBins = Lap1(1):Track_dbin:Lap1(end); % Spatial bin edges
TrackBins2 = Lap1(1):Track_dbin2:Lap1(end); % for spatial smoothing like Milstein2021
BinCenters = L/(2*Nbin):L/Nbin:L; % in cm
BinCenters2 = L/(2*Nbin2):L/Nbin2:L; % in cm
Run_bin = discretize(Run, TrackBins2); %binned position
for lap = 1:Nlaps
    if lap == Nlaps
    OutRaster_laps{lap} = OutRaster(NewLapIdx(lap):end);
    LapN = Run(NewLapIdx(lap):end);
    LapNbin = Run_bin(NewLapIdx(lap):end);
    CS_laps{lap} = CS(NewLapIdx(lap):end);
    Vlo_laps{lap} = Vlo(NewLapIdx(lap):end);
    else
    OutRaster_laps{lap} = OutRaster(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    LapN = Run(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    LapNbin = Run_bin(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    CS_laps{lap} = CS(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    Vlo_laps{lap} = Vlo(NewLapIdx(lap):NewLapIdx(lap+1)-1);   
    end
    Spiketimes_out{lap} = Tlap1(find(OutRaster_laps{lap}));
    SpikeLoc_out{lap} = Lap1(find(OutRaster_laps{lap}));
    CStimes{lap} = Tlap1(find(CS_laps{lap}));
    CSloc{lap} = Lap1(find(CS_laps{lap}));
    % compute place field
    SpikeCountOut(lap, :) = histcounts(SpikeLoc_out{lap}, TrackBins); % number of spikes in each spatial bin, for given lap
    CSbin(lap,:) = histcounts(CSloc{lap}, TrackBins); % number of CSs pin each spatial bins, for given lap
    TimeInBins(lap,:) = dt*histcounts(LapN, TrackBins); % compute time spent (in sec) in each spatial bin (no assumption on trajectories)
    SpatialFRout(lap,:) = SpikeCountOut(lap, :)./TimeInBins(lap, :); % FR in each spatial bin, in Hz 
    %Compute lap-wise COM, SD, skewness, mean and max firing rate
    COMbin(lap) = sum(SpatialFRout(lap,:).*[1:Nbin])/sum(SpatialFRout(lap,:));
    COMloc(lap) = sum(SpatialFRout(lap,:).*BinCenters)/sum(SpatialFRout(lap,:));
    PFsd(lap) = sqrt( sum( (BinCenters - COMloc(lap)).^2.*SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ) ); %lap-wise SD of the PF
    PFskew(lap) = sum( ((BinCenters-COMloc(lap))./PFsd(lap) ).^3.* SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ); %lap-wise skewness of the PF 
    meanFRout_lap(lap) = sum(OutRaster_laps{lap})/(dt*length(LapN)); % average FR over time, in Hz
    maxFRout_lap(lap) = max(SpatialFRout(lap,:)); % spatial bin with max FR, in Hz
   
    for b = 1:Nbin2
        Vbin(lap,b) = mean(Vlo_laps{lap}(LapNbin==b));  % average Vm for each spatial bin (using Milstein's binsize)
    end
end
meanFRmap = mean(SpatialFRout,1); %average FR map (i.e. average PF)
COM_meanPF = sum(meanFRmap.*BinCenters)/sum(meanFRmap);
SD_meanPF = sqrt( sum( (BinCenters - COM_meanPF).^2.*meanFRmap/sum(meanFRmap) ) );
Skew_meanPF = sum( ((BinCenters-COM_meanPF)/SD_meanPF ).^3.* meanFRmap/sum(meanFRmap) );

[b,~,~,~,stats] = regress(COMloc', [ones(Nlaps,1), [1:Nlaps]']);
COMtraj = b(1) + b(2)*[1:Nlaps];

if addCS == 1
    addCStime = CStimes{CSlap(1)}; % time of CS during induction lap, where 0 is start of CSlap
    addCSloc = CSloc{CSlap(1)}; % CS location on track, in cm (not binned)
    addCSbin = find(CSbin(CSlap(1),:)); % spatial bin of the added CS

% average of spatially binned Vm traces before or after induction +
% savitzky-golay smoothing (as in Milstein 2021)
sgOrder = 3; sgWinLen = 21;

muVbin.pre = mean(Vbin(1:CSlap(1)-1,:),1);
padStart.pre = muVbin.pre(end-sgWinLen:end); padEnd.pre = muVbin.pre(1:sgWinLen);
muVbin.pre_smooth = sgolayfilt([padStart.pre, muVbin.pre, padEnd.pre], sgOrder,sgWinLen);
muVbin.pre_smooth = muVbin.pre_smooth(sgWinLen+1:end-sgWinLen-1);

muVbin.post = mean(Vbin(CSlap(end)+1:end,:),1);
padStart.post = muVbin.post(end-sgWinLen:end); padEnd.post = muVbin.post(1:sgWinLen);
muVbin.post_smooth = sgolayfilt([padStart.post, muVbin.post, padEnd.post], sgOrder,sgWinLen);
muVbin.post_smooth = muVbin.post_smooth(sgWinLen+1:end-sgWinLen-1);

% peak
RelmuVbin.pre_smooth = muVbin.pre_smooth - Params.Vrest; % relative to Vrest. Baseline is now zero
RelmuVbin.post_smooth = muVbin.post_smooth - Params.Vrest; % use pre-induction baseline like in Milstein2021

    [Peak.pre, PeakIdx.pre] = max(RelmuVbin.pre_smooth);
    [Peak.post, PeakIdx.post] = max(RelmuVbin.post_smooth);
    
    PeakLoc.pre = BinCenters2(PeakIdx.pre); % location of the peak on the track
    PeakLoc.post = BinCenters2(PeakIdx.post);

% % ramp starts and ends at 15% of peak (cf Bittner2017), or location of lowest Vm, whichever is lower (min Vm may never be as low as 15% of the peak)
NormVbin.pre_smooth = RelmuVbin.pre_smooth./Peak.pre; % average Vm normalized to peak
NormVbin.post_smooth = RelmuVbin.post_smooth./Peak.post;

    nVmStart.pre = min( [min(NormVbin.pre_smooth), NormVbin.pre_smooth == 0.15] ); % norm Vm value at start of ramp
    StartIdx.pre = find(NormVbin.pre_smooth == nVmStart.pre); % , corresponding start idx
    rampStartLoc.pre = BinCenters2(StartIdx.pre);
    
    nVmStart.post = min( [min(NormVbin.post_smooth), NormVbin.post_smooth == 0.15] ); % norm Vm value at start of ramp, and corresponding start idx
    StartIdx.post = find(NormVbin.post_smooth == nVmStart.post);
    rampStartLoc.post = BinCenters2(StartIdx.post);

    % ramp end, rise length, decay length? 
    % Need to compute only if initial connectivity exactly like in Milstein2021 after 1st induction, 
    % to make sure ramp is not too big and long for the track, as it currently is with
    % params inherited from my STDP models (Imax and Wsd)

% spatial deltaVm (difference after induction)
DeltaVm_space = muVbin.post_smooth-muVbin.pre_smooth; % plot against BinCenter

    dVmPeak1 = DeltaVm_space(PeakIdx.pre); % DeltaVm at peak 1
    dVmPeak2 = DeltaVm_space(PeakIdx.post); % deltaVm at the new peak
    dVmCS = DeltaVm_space(addCSbin); % deltaVm at location of CS

% temporal pre, post and deltaVm
muVtime.pre = mean(cat(1, Vlo_laps{1:CSlap(1)-1}), 1);
muVtime.post = mean(cat(1,Vlo_laps{CSlap(end)+1:end}), 1);
DeltaVm_time = muVtime.post - muVtime.pre; % plot against Tlap1
    % peaks times and deltaVm at CStime
    [PeakT.pre, PeakTidx.pre] = max(smoothdata(muVtime.pre, 'sgolay'));
    [PeakT.post, PeakTidx.post] = max(smoothdata(muVtime.post, 'sgolay'));
    
    PeakTime.pre = Tlap1(PeakTidx.pre); % time of the peak during lap
    PeakTime.post = Tlap1(PeakTidx.post);
    
    dVmPeak1T = DeltaVm_time(PeakTidx.pre); % DeltaVm at peak 1
    dVmPeak2T = DeltaVm_time(PeakTidx.post); % deltaVm at the new peak
    dVmCStime = DeltaVm_time(Tlap1==addCStime); % deltaVm at time of CS

% center the spatial and temporal Vm traces on CS location and time: 
% (like in Milstein for easier comparison, with negative values corresponding to before the CS, i.e. not like in my learning rule display)
CScenteredLoc = BinCenters2 - addCSloc;
CScenteredTime = Tlap1 - addCStime;

% effective Pamp ~ infering dWmax
Wchange_prepost = W2 - W(:,1); % weight change after CS induction
dWeff = diff(W2, 1, 2); % instantaneous weight change at every time point (columns) for all synapses (rows)
[effdWmaxAll, dWmaxAllIdx] = max(dWeff, [], 2); 
[effdWmax, dWmaxIn] = max(effdWmaxAll);
fdWmax_time = Trun(dWmaxAllIdx(dWmaxIn)); % time point where the max weight change occurs
effPamp = effdWmax/P(dWmaxIn,dWmaxAllIdx(dWmaxIn)) % deduce effective Pamp
[effdWcsbin, fdWcsbin_idx] = max(dWeff(addCSbin,:)); % focus on input synapse active when CS occured. Should be where dW is highest (regardless of whether dW is dependent on synaptic weight or not)

end
%% figures
close all

figure % input current and Vm on last lap
    subplot(2,1,1)
        plot(Tlap1, I(NewLapIdx(end):end).*10^12) %plot(Tlap1, I(1:NewLapIdx(2)-1).*10^12, 'k'); 
        xlabel('time (s)'); ylabel('input current (pA)')
        box off;
        title('input current on last lap');
    subplot(2,1,2)
        plot(Tlap1, V(NewLapIdx(end):end).*10^3) %plot(Tlap1, V(1:NewLapIdx(2)-1).*10^3, 'k'); 
        xlabel('time (s)'); ylabel('membrane potential (mV)')
        box off;
        title('Vm on last lap (LIF model)');

figure % lapwise PF
[row_lap,column_bin] = find(CSbin);
imagesc(SpatialFRout); hold on
scatter(column_bin, row_lap, 10, [0 1 1], 'filled'); hold on
scatter(COMbin,1:Nlaps, 10, [1 0 0], 'filled');
legend('Complex spike','COM', 'Location', 'BestOutside');
xlabel('spatial bins'); ylabel('lap');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Firing rate (Hz)';
box off; axis square

figure %BTSP variables on given lap
Wsyn =77e-12; % in Amps
[val, Nbittner] = min(abs(W(:,1)-Wsyn));

LapNum = CSlap(2);
InputNum = dWmaxIn; %N/2; %Nbittner
x1 = 0; % in s
x2 = period; % in s 
CSidx = find(CS_laps{LapNum});

subplot(3,1,1)
    plot(Tlap1, D(NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'r'); hold on
    plot(repmat(Spiketimes_out{LapNum},2,1), repmat([0; 1], 1, length(Spiketimes_out{LapNum})), 'k'); hold on
%     stem(Tlap1, CS_laps{LapNum}, 'c');
    scatter(Tlap1(CSidx), ones(size(CSidx)).*0.5, 'oc', 'filled')
    xlim([x1 x2]);
    xlabel('time (s)'); 
    title('output raster, CS raster, Post-before-Pre variable')
    box off;
subplot(3,1,2)
    plot(Tlap1, max(P(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)).*InRasters(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'k'); hold on
    plot(Tlap1, P(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'r'); hold on
    xlim([x1 x2]);
    xlabel('time (s)'); 
    title(['input raster #' num2str(InputNum) '+ Pre-before-Post variable'])
    box off;
subplot(3,1,3)
    plot(Tlap1, W(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'b--'); hold on
    plot(Tlap1, W2(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'r'); hold on
    if Params.WUdyn == 1
    plot(Tlap1, Wtarget(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'k'); hold on
    end
    xlim([x1 x2]);
    xlabel('time (s)'); 
    ylabel('Norm. synaptic strength')
    title(['input #' num2str(InputNum)] )
    box off;

figure %BTSP variables on end lap
subplot(2,1,1)
    plot(Tlap1, D(NewLapIdx(end):end), 'r'); hold on %plot(Tlap1, D(1:NewLapIdx(2)-1), 'r'); hold on
    plot(Tlap1, OutRaster_laps{Nlaps}, 'k'); hold on
    stem(Tlap1, CS_laps{Nlaps}, 'r');
    xlabel('time (s)'); 
    title('output raster + Post-before-Pre weight update')
    box off;
subplot(2,1,2)
    plot(Tlap1, P(1,NewLapIdx(end):end), 'r'); hold on
    plot(Tlap1, InRasters(1,NewLapIdx(end):end), 'k');
    xlabel('time (s)'); 
    title('input raster #1 + Pre-before-Post weight update')
    box off;

figure % Weight changes over time, full run
subplot(2,1,1)
imagesc([0 Nlaps*period],[1 N], W.*10^12); hold on
scatter(NewLapTimes, (N/2)*ones(1,Nlaps), 10, [1 0 0], '+');
xlabel('time (s)'); ylabel('input neuron');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Synaptic weight (pA)';
title('W')
box off;
subplot(2,1,2)
imagesc([0 Nlaps*period],[1 N], W2.*10^12); hold on
scatter(NewLapTimes, (N/2)*ones(1,Nlaps), 10, [1 0 0], '+');
xlabel('time (s)'); ylabel('input neuron');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Synaptic weight (pA)';
title('W2')
box off;

% figure
% plot(Trun, W(25,:).*10^12, 'r'); hold on
% plot(Trun, W(15,:).*10^12, 'k');
% scatter(NewLapTimes, W(1,NewLapIdx).*10^12, 10, [0 0 1], '+'); 
% legend('input neuron 25', 'input neuron 15', 'input 1 (new lap start)')
% xlabel('time (s)'); ylabel('input current (pA)')
% title('Synaptic weight of input neuron #1');
% box off;

figure % PF properties
    subplot(2,2,1)
        plot(1:Nlaps, meanFRout_lap, 'r+'); hold on
        plot(1:Nlaps, maxFRout_lap, 'k-');
        xlabel('lap'); ylabel('Firing Rate (Hz)')
        legend('mean','max', 'Location', 'BestOutside');
        box off; axis square;
    subplot(2,2,2)
        plot(1:Nlaps, PFsd, 'k-');
        xlabel('lap'); ylabel('PF sd (cm)');
        box off; axis square;
    subplot(2,2,3)
        plot(1:Nlaps, PFskew, 'k-');
        xlabel('lap'); ylabel('PF skewness');
        box off; axis square;
    subplot(2,2,4)
        plot(BinCenters, meanFRmap, 'k-'); hold on
        xline(COM_meanPF, 'r'); hold on
        xline(L/2, 'k--'); 
        xlabel('position'); ylabel('mean Firing Rate (Hz) ');
        legend('mean PF', 'COM', 'track center', 'Location', 'BestOutside');
        title({'SD = ' num2str(SD_meanPF) ', Skew = ' num2str(Skew_meanPF)})
        box off; axis square;

figure % PF trajectory and linear regression
plot(1:Nlaps, COMloc, 'k-'); hold on
plot(1:Nlaps, COMtraj, 'b-'); hold on
yline(COMloc(1),'r');
legend('COM','lin reg','COM #1')
ylim([0 300]);
xlabel('lap'); ylabel('COM position (cm)');
title({'Place Field COM trajectory, slope =' num2str(b(2)) 'pval =' num2str(stats(3))})

if addCS == 1

    figure % Vm (raw and low-pass filtered) on lap 1 and end lap
% LapNum = CSlap;
    subplot(2,2,1)
        plot(Tlap1, V(NewLapIdx(1):NewLapIdx(2)-1).*10^3) %plot(Tlap1, V(1:NewLapIdx(2)-1).*10^3, 'k'); 
        xlabel('time (s)'); ylabel('membrane potential (mV)')
        box off;
        title('first lap');
    subplot(2,2,2)
        plot(Tlap1, V(NewLapIdx(end):end).*10^3)
        xlabel('time (s)'); ylabel('membrane potential (mV)')
        box off;
        title('last lap');
%         title(['Lap ' num2str(LapNum)] )
    subplot(2,2,3)
        plot(Tlap1, Vlo(NewLapIdx(1):NewLapIdx(2)-1).*10^3) %plot(Tlap1, V(1:NewLapIdx(2)-1).*10^3, 'k'); 
        xlabel('time (s)'); ylabel('low-passed Vm (mV)')
        box off;
        title('first lap');
    subplot(2,2,4)
        plot(Tlap1, Vlo(NewLapIdx(end):end).*10^3)
        xlabel('time (s)'); ylabel('low-passed Vm (mV)')
        box off;
        title('last lap');
%         title(['Lap ' num2str(LapNum)] )

t2 = tiledlayout('flow', 'TileSpacing','Compact','Padding','Compact')
title(t2,['Pamp = ' num2str(Params.Pamp.*10^12) 'pA'])
% figure % lap-wise Vm spatially binned (like PF but subthreshold Vm instead of FR)
% [row_lap,column_bin] = find(CSbin);

ax2 = nexttile
    imagesc(Vbin.*10^3); hold on
    scatter(column_bin.*2, row_lap, 10, [0 1 1], 'filled'); hold on
%     legend('CS','COM', 'Location', 'BestOutside');
    xlabel('spatial bins'); ylabel('lap');
    % colormap(hot(256)); 
    colormap(ax2, brewermap(256,'*YlOrRd'));
    c2 = colorbar; c2.Label.String = 'membrane potential (mV)'; c2.Ticks = [min(Vbin.*10^3,[], 'all') max(Vbin.*10^3,[], 'all')]; c2.TickLabels = [min(Vbin.*10^3,[], 'all') max(Vbin.*10^3,[], 'all')];
    % caxis([-75 -55])
    box off; axis square
nexttile
    plot(BinCenters2, muVbin.pre_smooth.*10^3, 'k-'); hold on
    plot(BinCenters2, muVbin.post_smooth.*10^3, 'b-');
    xline(addCSloc, 'r')
%     legend('pre', 'post', 'CS', 'Location', 'BestOutside')
    xlabel('position'); ylabel('membrane potential (mV)')
    axis square
    box off;
    title('pre and post spatial Vm');
ax1 = nexttile
    imagesc([0 Nlaps*period],[1 N], W2.*10^12); hold on
    scatter(NewLapTimes, (N/2)*ones(1,Nlaps), 10, [1 0 0], '+');
%     legend('new lap', 'location', 'bestoutside')
    xlabel('time (s)'); ylabel('input neuron');
    colormap(ax1, flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
    c = colorbar; c.Label.String = 'Synaptic weight (pA)'; c.Ticks = [0 max(W2.*10^12,[], 'all')]; c.TickLabels = [0 max(W2.*10^12,[], 'all')];
    title('W2')
    box off; axis square
nexttile
    plot(CScenteredLoc, DeltaVm_space.*10^3, 'k-'); hold on
    scatter(CScenteredLoc(PeakIdx.pre), dVmPeak1.*10^3, 'k^', 'filled');hold on
    scatter(CScenteredLoc(PeakIdx.post), dVmPeak2.*10^3, 'b^', 'filled');hold on
    xline(0, 'r')
%     legend('\DeltaVm', 'Peak 1', 'Peak 2', 'CS', 'Location', 'BestOutside')
    xlim([floor(min(CScenteredLoc)) ceil(max(CScenteredLoc))]);
%     ylim([-9 9])
    xlabel('distance from CS (cm)'); ylabel('membrane potential (mV)')
    box off; axis square;
    title('spatial delta Vm (post-pre induction) ');

% exportgraphics(t2,'BTSP_MilsteinExp1CS_MilsteinParams_80pA_CS3sPostCOM1.pdf','BackgroundColor','none', 'ContentType','vector')

% average spatial Vm ("raw" and smoothed) before and after induction + show CS time
figure
    subplot(2,2,1)
        plot(BinCenters2, muVbin.pre.*10^3, 'k-')
        xlabel('position'); ylabel('membrane potential (mV)')
        box off;
        title('laps before CS');
    subplot(2,2,2)
        plot(BinCenters2, muVbin.post.*10^3, 'k-')
        xline(addCSloc, 'c')
        xlabel('position'); ylabel('membrane potential (mV)')
        box off;
        title('laps after CS');
%         title(['Lap ' num2str(LapNum)] )
    subplot(2,2,3)
        plot(BinCenters2, muVbin.pre_smooth.*10^3, 'k-')
        xlabel('position'); ylabel('membrane potential (mV)')
        box off;
        title('laps before CS (smoothed)');
    subplot(2,2,4)
        plot(BinCenters2, muVbin.post_smooth.*10^3, 'k-')
        xline(addCSloc, 'c')
        xlabel('position'); ylabel('membrane potential (mV)')
        box off;
        title('last after CS (smoothed)');

figure % average temporal Vm before and after induction   
    subplot(1,2,1)
        plot(Tlap1, muVtime.pre.*10^3, 'k-'); hold on
        xlabel('time (s)'); ylabel('membrane potential (mV)')
        box off;
        title('laps before CS');
    subplot(1,2,2)
        plot(Tlap1, muVtime.post.*10^3, 'k-')
        xline(addCStime, 'c')
        xlabel('time (s)'); ylabel('membrane potential (mV)')
        box off;
        title('laps after CS');   

figure % initial Vm vs delta Vm
    subplot(1,2,1) % using spatial smoothed values
        scatter(muVbin.pre_smooth.*10^3, DeltaVm_space.*10^3, 'ok')
        xlabel('pre-induction Vm (mV)'); ylabel('deltaVm (mV)')
        box off;
        title('Spatial (smoothed)'); 
    subplot(1,2,2)% using temporal
        scatter(muVtime.pre.*10^3, DeltaVm_time.*10^3, 'ok')
        xlabel('pre-induction Vm (mV)'); ylabel('deltaVm (mV)')
        box off;
        title('Temporal');

figure % final Vm vs delta Vm
    subplot(1,2,1) % using spatial smoothed values
        scatter(muVbin.pre_smooth.*10^3, muVbin.post_smooth.*10^3, 'ok')
        xlabel('pre-induction Vm (mV)'); ylabel('post-induction Vm (mV)')
        box off;
        title('spatial (smoothed)'); 
    subplot(1,2,2)% using temporal
        scatter(muVtime.pre.*10^3, muVtime.post.*10^3, 'ok')
        xlabel('pre-induction Vm (mV)'); ylabel('post-induction Vm (mV)')
        box off;
        title('Temporal');

% figure % spatial delta Vm (post-pre induction) + show CS location
% plot(BinCenters, (muVbin.pre_smooth-muVbin.post_smooth).*10^3, 'k-')
% xline(addCSloc, 'c')
% xlabel('position'); ylabel('membrane potential (mV)')
% box off;
% title('post - pre spatial Vm ');

figure % spatial delta Vm (post-pre induction) centered on CS position
plot(CScenteredLoc, DeltaVm_space.*10^3, 'k-'); hold on
scatter(CScenteredLoc(PeakIdx.pre), dVmPeak1.*10^3, 'k^');hold on
scatter(CScenteredLoc(PeakIdx.post), dVmPeak2.*10^3, 'b^');hold on
xline(0, 'c')
xlabel('distance from CS (cm)'); ylabel('membrane potential (mV)')
box off;
title('spatial delta Vm (post-pre induction) '); 

figure % temporal delta Vm (post-pre induction) centered on CS time
plot(CScenteredTime, DeltaVm_time.*10^3, 'k-'); hold on
scatter(CScenteredTime(PeakTidx.pre), dVmPeak1T.*10^3, 'k^');hold on
scatter(CScenteredTime(PeakTidx.post), dVmPeak2T.*10^3, 'b^');hold on
xline(0, 'c')
xlabel('time from CS'); ylabel('membrane potential (mV)')
box off;
title('temporal delta Vm (post-pre induction) '); 

figure % ramp relative to baseline and normalized to peak
    subplot(2,2,1)
        plot(BinCenters2, RelmuVbin.pre_smooth.*10^3, 'k-')
        xlabel('position'); ylabel('membrane potential relative to baseline (mV)')
        box off;
        title('laps before CS');
    subplot(2,2,2)
        plot(BinCenters2, RelmuVbin.post_smooth.*10^3, 'k-')
        xlabel('position'); ylabel('membrane potential relative to baseline (mV)')
        box off;
        title('laps after CS');
    subplot(2,2,3)
        plot(BinCenters2, NormVbin.pre_smooth, 'k-')
        xlabel('position'); ylabel('norm Vm')
        box off;
        title('laps before CS');
    subplot(2,2,4)
        plot(BinCenters2, NormVbin.post_smooth, 'k-')
        xlabel('position'); ylabel('norm Vm')
        box off;
        title('laps after CS');

figure    % zoom on BTSP variables + effective weight update
x1 = addCStime-0.3; % in s
x2 = addCStime+0.3; % in s 

subplot(4,1,1)
    plot(Tlap1, D(NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'r'); hold on
    plot(repmat(Spiketimes_out{LapNum},2,1), repmat([0; 1], 1, length(Spiketimes_out{LapNum})), 'k'); hold on
%     stem(Tlap1, CS_laps{LapNum}, 'c');
    scatter(Tlap1(CSidx), ones(size(CSidx)).*0.5, 'oc', 'filled')
    xlim([x1 x2]);
    xlabel('time (s)'); 
    title('output raster, CS raster, Post-before-Pre variable')
    box off;
subplot(4,1,2)
    plot(Tlap1, P(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'r'); hold on
    plot(Tlap1, max(P(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)).*InRasters(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'k');
    xlim([x1 x2]);
    xlabel('time (s)'); 
    title(['input raster #' num2str(InputNum) '+ Pre-before-Post variable'])
    box off;
subplot(4,1,3)
    plot(Tlap1, W(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'b--'); hold on
    plot(Tlap1, W2(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'r'); hold on
    if Params.WUdyn == 1
    plot(Tlap1, Wtarget(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'k'); hold on
    end
    xlim([x1 x2]);
    xlabel('time (s)'); 
    ylabel('Norm. synaptic strength')
    title(['input #' num2str(InputNum)] )
    box off;
subplot(4,1,4)    
%     plot(Tlap1, Wchange_prepost(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1).*10^12, 'r'); hold on
    plot(Tlap1, dWeff(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1).*10^12, 'g'); hold on
    xlim([x1 x2]);
    xlabel('time (s)'); 
    ylabel('dWeff (pA)')
    title(['input #' num2str(InputNum)] )
    box off;

t1 = tiledlayout(5, 1, 'TileSpacing','Compact','Padding','Compact')
%     figure    % zoom on BTSP variables + effective weight update
% subplot(5,1,1)
nexttile
    plot(Tlap1, D(NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'Color', [0.5 0.2 0.5]); hold on
    plot(repmat(Spiketimes_out{LapNum},2,1), repmat([0; 1], 1, length(Spiketimes_out{LapNum})), 'k'); hold on
%     stem(Tlap1, CS_laps{LapNum}, 'c');
    scatter(Tlap1(CSidx), ones(size(CSidx)).*0.5, 'oc', 'filled')
    xlim([x1 x2]);
    xlabel('time (s)'); 
    legend('Post-before-Pre var.', 'CS', 'output spikes', 'Location', 'bestoutside')
    box off; 
% subplot(5,1,2)
nexttile
    plot(Tlap1, P(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'Color', [1    0.5    0]); hold on
    plot(Tlap1, max(P(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)).*InRasters(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'k');
    legend('Pre-before-Post var.', 'input spikes', 'Location', 'bestoutside')
    xlim([x1 x2]);
    xlabel('time (s)'); 
    title(['input #' num2str(InputNum)])
    box off; 
% subplot(5,1,3)
nexttile
    plot(Tlap1, W(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'b--'); hold on
    plot(Tlap1, W2(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'r'); hold on
    if Params.WUdyn == 1
    plot(Tlap1, Wtarget(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'k'); hold on
    end
    xlim([x1 x2]);
    xlabel('time (s)'); 
    ylabel('Norm. synaptic strength')
    legend('BTSP alone', 'BTSP + syn. norm.', 'Location', 'bestoutside')
    title(['input #' num2str(InputNum)] )
    box off; 
% subplot(5,1,4)
nexttile
    plot(Tlap1, P(50, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'Color', [1    0.5    0]); hold on
    plot(Tlap1, max(P(50, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)).*InRasters(50,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'k');
    legend('Pre-before-Post var.', 'input spikes', 'Location', 'bestoutside')
    xlim([x1 x2]);
    xlabel('time (s)'); 
    title('input #50')
    box off; 
% subplot(5,1,5)
nexttile
    plot(Tlap1, W(50,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'b--'); hold on
    plot(Tlap1, W2(50,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'r'); hold on
    if Params.WUdyn == 1
    plot(Tlap1, Wtarget(50,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2(:)), 'k'); hold on
    end
    xlim([x1 x2]);
    xlabel('time (s)'); 
    ylabel('Norm. synaptic strength')
    legend('BTSP alone', 'BTSP + syn. norm.', 'Location', 'bestoutside')
    title('input #50' )
    box off; 

%     exportgraphics(t1,'BTSPvariables_MilsteinParams_80pA_CS1sPostCOM.pdf','BackgroundColor','none', 'ContentType','vector')
end
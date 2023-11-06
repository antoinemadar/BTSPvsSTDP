%InVivoMilstein2021_1CSbatch_PampSpace
clear
close all

%% Parameters
%Params for inputs
N = 100; % number of input neurons
L = 300 % length of track in cm (300 for Can, 185 for Bittner and Milstein). For Milstein experiments, I use 2*185 to make sure Vm ramps are not bleeding into next lap + see what happens far from CS, but still use same smoothing parameters
PFsd = 18; % cm (~median in Can data) 
PFamp = 10; %peak FR in Hz
Nlaps = 21 % 30; % number of laps
speed = 15; % cm/sec, (15cm/s average in DongSheffield2021, 25cm/s in Milstein2021 network model)
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

Params.PostPreBoost = 1.1; % scaling factor, as a function of Pamp, for the postpre variable (triggered on CS, not on spikes, so with little to no opportunities for temporal summation). Optimized on Bittner2017 invitro data
Params.capWeights = 0; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 
Params.HomeoNorm = 1; % 0 if no homeostatic normalization rule, 1 if yes.
Params.WUdyn = 0; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Params.Wrest = 0; % Synaptic weight resting state, in Amps (towards which any synapse would relax if Wtarget = 0) 
Params.tau_wu = 5; % weight update time constant, in sec

Params.CSproba = []; % proba of an output spike to be a complex spike (make it an empty array when you want CSs in specific locations)
CSlap = 11; % lap on which CS will occur 
repeats = 1; % number of repeats with same CS time

% output PF 
Nbin = 50; % Number of bins in which the length of the track is divided (50 bins of 6cm in DongSheffield2021, 100 bins of 1.85cm in Milstein2021)
Nbin2 = round(L/1.85); % for spatial smoothing like in Milstein2021

Params.dt = dt; Params.Imax = Imax;

%% Param to vary: Pamp
% PampRange = [4,4,4,6,6,6,8,8,8,10,10,10,15,15,15,20,20,20,30,30,30,40,40,40,60,60,60,80,80,80].*10^-12; %in Amp
PampRange = [4,4,4,6,6,6,8,8,8,10,10,10,11,11,11,12,12,12,13,13,13,15,15,15,20,20,20,30,30,30,60,60,60,80,80,80].*10^-12; %in Amp

%% design FIR for smoothing Vm: with wrap-around padding of Hamming window size and zero-phase filtering
Fcut = 3; %cut-off frequency, in Hz
Fs = 1./dt; %sampling frequency in Hz
NyqF = Fs/2; % Nyquist frequency
HWin = 2; % size of hamming window, in seconds
fir_n = HWin.*Fs-1; % filter order
nFcut = Fcut/NyqF; % normalized cutoff frequency (proprtion of the Nyquist frequency)
blo = fir1(fir_n,nFcut); % design a FIR filter of fir_n'th order, using hamming window, with cut-off frequency at 3Hz
%         freqz(blo,1,2^16,Fs)

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

% segment in laps and spatial bins
NewLapIdx = find(Run==0);
NewLapTimes = Trun(NewLapIdx);
Track_dbin = L./Nbin; % bin size, in cm
Track_dbin2 = L./Nbin2; %for spatial smoothing like Milstein2021
TrackBins = Lap1(1):Track_dbin:Lap1(end); % Spatial bin edges
TrackBins2 = Lap1(1):Track_dbin2:Lap1(end); % for spatial smoothing like Milstein2021
BinCenters = L/(2*Nbin):L/Nbin:L; % in cm
BinCenters2 = L/(2*Nbin2):L/Nbin2:L; % in cm
Run_bin = discretize(Run, TrackBins2); %binned position

% Time bins for downsampling
TimeBins = linspace(Tlap1(1),Tlap1(end),Nbin2+1); % bin edges with Nbin2 time bins
TimeDS = linspace(TimeBins(2)/2, TimeBins(end)-TimeBins(2)/2, Nbin2); % bin centers
Tlap1_bin = discretize(Tlap1, TimeBins);

%convert trajectory into firing rate
for i = 1:length(Run)
    [minval, idx] = min(abs(x-Run(i))); %find closest spatial bin in x for each position in trajectory
    FRrun(:,i) = FRin(:,idx); % Firing rate at that position for each neuron
end

% gaussian connectivity with max in center of track, akin to Mehta and Wilson 2000
NeuronID = 1:N;
W = zeros(N,length(Trun));
Wmax = maxW*Imax;
W(:,1) = Wmax*exp(-0.5*(NeuronID - N/2).^2/Wsd^2); % vector of N rows following a gaussian distribution of Wmax amplitude 

%% loop through PampRange

for u = 1:length(PampRange)
  
    Params.Pamp = PampRange(u); % peak weight change, in Amps (3*baseline EPSP for pairings with 5*10 EPSPs at 10Hz, with initial EPSP at 2mV, in Bittner et al. 2017)

    
    %% input-output dynamics

for n = 1:4:N % vary location of CS along COM of input neurons. Use range that matches Milstein 2021 (0 to 140cm between CS and PF. Actually, the 140cm extrema correspond to initial peak near edge of track. I cannot go that far. My range is 0 to 185/2=92.5cm if I use a 185cm track like them. But I could use a longer track... 
    
    % make CS trace with 1 CS at time corresponding to COM of input n
    CS{n} = zeros(1,length(Trun)); 
    [minvalN(n), COMidxN(n)] = min(abs(Lap1-PFcom(n)));
    COMidxN_Laps = find(Run==Run(COMidxN(n))); % COM idx of input cell n on all laps (assuming constant speed)
    CS{n}(COMidxN_Laps(CSlap)) = 1;
    
    for r = 1:repeats
        
        % nonhomogenous Poisson spike generator along trajectory
        InRasters{n,r} = poissrnd(FRrun.*dt); %matrix of size (N,length(Tlap1)) containing, for each input neuron, the number of spikes in each 1ms bin drawn from Poisson distribution
        
        % run LIF + plasticity model
        [ ~, ~, V2{n,r}, CS{n,r}, OutRaster{n,r}, D{n,r}, P{n,r}, ~, ~, ~, ~, W2{n,r}, ~ ]  = BTSPplus_LIFadapt_ODEs(Params, W, InRasters{n,r}, CS{n});
        
%         Spiketimes_outAll = Trun(logical(OutRaster{n,r}));
%         SpikeLoc_outAll = Run(logical(OutRaster{n,r})); %spikes location on track
%         CSloc_All = Run(logical(CS{n})); % CSs location on track 
        
        % smooth Vm traces by using a 3Hz low-pass filter like in Bittner2017     
        Vpad = [V2{n,r}(length(V2{n,r})-HWin*Fs:end), V2{n,r}, V2{n,r}(1:HWin*Fs)]; % zero-pad V2 on the left to avoid weird thing in beginning
        Vlo{n,r} = filtfilt(blo, 1, Vpad);
        Vlo{n,r} = Vlo{n,r}(HWin*Fs+1:end-HWin*Fs-1); % remove pad sections =>should be same length as V2
        clear Vpad

        % compute spatial firing rate and related PF properties
        for lap = 1:Nlaps
            if lap == Nlaps
            OutRaster_laps{lap} = OutRaster{n,r}(NewLapIdx(lap):end);
            LapN = Run(NewLapIdx(lap):end);
            LapNbin = Run_bin(NewLapIdx(lap):end);
            CS_laps{n}{lap} = CS{n}(NewLapIdx(lap):end);
            Vlo_laps{n,r}{lap} = Vlo{n,r}(NewLapIdx(lap):end);
            else
            OutRaster_laps{lap} = OutRaster{n,r}(NewLapIdx(lap):NewLapIdx(lap+1)-1);
            LapN = Run(NewLapIdx(lap):NewLapIdx(lap+1)-1);
            LapNbin = Run_bin(NewLapIdx(lap):NewLapIdx(lap+1)-1);
            CS_laps{n}{lap} = CS{n}(NewLapIdx(lap):NewLapIdx(lap+1)-1);
            Vlo_laps{n,r}{lap} = Vlo{n,r}(NewLapIdx(lap):NewLapIdx(lap+1)-1);   
            end
            Spiketimes_out{n,r}{lap} = Tlap1(find(OutRaster_laps{lap}));
%             SpikeLoc_out{n,r}{lap} = Lap1(find(OutRaster_laps{n,r}{lap}));
            CStimes{n}{lap} = Tlap1(find(CS_laps{n}{lap}));
            CSloc{n}{lap} = Lap1(find(CS_laps{n}{lap}));
            % compute place field
%             SpikeCountOut{n,r}(lap, :) = histcounts(SpikeLoc_out{n,r}{lap}, TrackBins); % number of spikes in each spatial bin, for given lap
            CSbin{n}(lap,:) = histcounts(CSloc{n}{lap}, TrackBins); % number of CSs in each spatial bins, for given lap
            CSbin2{n}(lap,:) = histcounts(CSloc{n}{lap}, TrackBins2); % number of CSs in each spatial bins, for given lap
%             TimeInBins(lap,:) = dt*histcounts(LapN, TrackBins); % compute time spent (in sec) in each spatial bin (no assumption on trajectories)
%             SpatialFRout{n,r}(lap,:) = SpikeCountOut{n,r}(lap, :)./TimeInBins(lap, :); % FR in each spatial bin, in Hz 
            %Compute lap-wise COM, SD, skewness, mean and max firing rate
%             COMbin{n,r}(lap) = sum(SpatialFRout{n,r}(lap,:).*[1:Nbin])/sum(SpatialFRout{n,r}(lap,:));
%             COMloc{n,r}(lap) = sum(SpatialFRout{n,r}(lap,:).*BinCenters)/sum(SpatialFRout{n,r}(lap,:));
%             PFsd(lap) = sqrt( sum( (BinCenters - COMloc(lap)).^2.*SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ) ); %lap-wise SD of the PF
%             PFskew(lap) = sum( ((BinCenters-COMloc(lap))./PFsd(lap) ).^3.* SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ); %lap-wise skewness of the PF 
%             meanFRout_lap(lap) = sum(OutRaster_laps{lap})/(dt*length(LapN)); % average FR over time, in Hz
%             maxFRout_lap(lap) = max(SpatialFRout(lap,:)); % spatial bin with max FR, in Hz
           
            for b = 1:Nbin2
                Vbin{n,r}(lap,b) = mean(Vlo_laps{n,r}{lap}(LapNbin==b));  % average Vm for each spatial bin (using Milstein's binsize)
                Vtime{n,r}(lap,b) = mean(Vlo_laps{n,r}{lap}(Tlap1_bin==b)); %average Vm for each temporal bin (as many as spatial bins)
            end

        end
%         meanFRmap{n,r} = mean(SpatialFRout{n,r},1); %average FR map (i.e. average PF)
%         COM_meanPF{n,r} = sum(meanFRmap{n,r}.*BinCenters)/sum(meanFRmap{n,r});
%         SD_meanPF = sqrt( sum( (BinCenters - COM_meanPF).^2.*meanFRmap/sum(meanFRmap) ) );
%         Skew_meanPF = sum( ((BinCenters-COM_meanPF)/SD_meanPF ).^3.* meanFRmap/sum(meanFRmap) );
        
        addCStime{n} = CStimes{n}{CSlap}; % time of CS during induction lap, where 0 is start of CSlap
        [~, addCStimebin{n}] = min(abs(TimeDS-addCStime{n}));
        addCSidx{n} = find(CS_laps{n}{CSlap}); % equivalent to: find(Tlap1==addCStime{n});
        addCSloc{n} = CSloc{n}{CSlap}; % CS location on track, in cm (not binned)
        addCSlocR{n,r} = CSloc{n}{CSlap}; % tracked for all repetitions r
        addCSbin2{n} = find(CSbin2{n}(CSlap,:)); % spatial bin of the added CS (based on nbin2, i.e. number of bins as in Milstein2021)
        
        % average of spatially binned Vm traces before or after induction +
        % savitzky-golay smoothing (as in Milstein 2021)
        sgOrder = 3; sgWinLen = 21;
        
        muVbin{n,r}.pre = mean(Vbin{n,r}(1:CSlap-1,:),1);
        padStart.pre = muVbin{n,r}.pre(end-sgWinLen:end); padEnd.pre = muVbin{n,r}.pre(1:sgWinLen);
        muVbin{n,r}.pre_smooth = sgolayfilt([padStart.pre, muVbin{n,r}.pre, padEnd.pre], sgOrder,sgWinLen);
        muVbin{n,r}.pre_smooth = muVbin{n,r}.pre_smooth(sgWinLen+1:end-sgWinLen-1);
        
        muVbin{n,r}.post = mean(Vbin{n,r}(CSlap+1:end,:),1);
        padStart.post = muVbin{n,r}.post(end-sgWinLen:end); padEnd.post = muVbin{n,r}.post(1:sgWinLen);
        muVbin{n,r}.post_smooth = sgolayfilt([padStart.post, muVbin{n,r}.post, padEnd.post], sgOrder,sgWinLen);
        muVbin{n,r}.post_smooth = muVbin{n,r}.post_smooth(sgWinLen+1:end-sgWinLen-1);
        clear padStart padEnd
        
        % peak
        RelmuVbin{n,r}.pre_smooth = muVbin{n,r}.pre_smooth - Params.Vrest; % relative to Vrest. Baseline is now zero
        RelmuVbin{n,r}.post_smooth = muVbin{n,r}.post_smooth - Params.Vrest; % use pre-induction baseline like in Milstein2021
        
            [Peak{n,r}.pre, PeakIdx{n,r}.pre] = max(RelmuVbin{n,r}.pre_smooth);
            [Peak{n,r}.post, PeakIdx{n,r}.post] = max(RelmuVbin{n,r}.post_smooth);
            
            PeakLoc{n,r}.pre = BinCenters2(PeakIdx{n,r}.pre); % location of the peak on the track
            PeakLoc{n,r}.post = BinCenters2(PeakIdx{n,r}.post);

       % Distance of Peak 1 from peak 2 (new PF = y-axis) Peak1 from CS (distance from plateau = x-axis)
       DistPeaks{n,r} = abs(PeakLoc{n,r}.post - PeakLoc{n,r}.pre);
       DistPeak1toCS{n,r} = abs(PeakLoc{n,r}.pre - addCSloc{n});
       DistPeak2toCS{n,r} = abs(PeakLoc{n,r}.post - addCSloc{n});
        
%         % % ramp starts and ends at 15% of peak (cf Bittner2017), or location of lowest Vm, whichever is lower (min Vm may never be as low as 15% of the peak)
%         NormVbin.pre_smooth = RelmuVbin.pre_smooth./Peak.pre; % average Vm normalized to peak
%         NormVbin.post_smooth = RelmuVbin.post_smooth./Peak.post;
%         
%             nVmStart.pre = min( [min(NormVbin.pre_smooth), NormVbin.pre_smooth == 0.15] ); % norm Vm value at start of ramp
%             StartIdx.pre = find(NormVbin.pre_smooth == nVmStart.pre); % , corresponding start idx
%             rampStartLoc.pre = BinCenters2(StartIdx.pre);
%             
%             nVmStart.post = min( [min(NormVbin.post_smooth), NormVbin.post_smooth == 0.15] ); % norm Vm value at start of ramp, and corresponding start idx
%             StartIdx.post = find(NormVbin.post_smooth == nVmStart.post);
%             rampStartLoc.post = BinCenters2(StartIdx.post);
        
            % ramp end, rise length, decay length? 
            % Need to compute only if initial connectivity exactly like in Milstein2021 after 1st induction, 
            % to make sure ramp is not too big and long for the track, as it currently is with
            % params inherited from my STDP models (Imax and Wsd)
        
%         % spatial deltaVm (difference after induction)
        DeltaVm_space{n,r} = muVbin{n,r}.post_smooth-muVbin{n,r}.pre_smooth; % plot against BinCenters2
%         
%             dVmPeak1{n,r} = DeltaVm_space{n,r}(PeakIdx{n,r}.pre); % DeltaVm at peak 1 position
%             dVmPeak2{n,r} = DeltaVm_space{n,r}(PeakIdx{n,r}.post); % deltaVm at the new peak position
%             dVmCS{n,r} = DeltaVm_space{n,r}(addCSbin2{n}); % deltaVm at location of CS
% 
%             % pad with NaNs in order to center on CS and, later, compute average centered on CS. 
%             NaNpadL = NaN(1,Nbin2-addCSbin2{n});
%             NaNpadR = NaN(1,addCSbin2{n}-1);
%             DeltaVm_spacePad{n,r} = [NaNpadL, DeltaVm_space{n,r}, NaNpadR]; % each trace is thus 2*size(BinCenters2)-1= (2*Nbin2)-1
%             clear NaNpadL NaNpadR
        
        % temporal pre, post and deltaVm (smooth and downsample)
%         muVtime.pre = mean(cat(1, Vlo_laps{n,r}{1:CSlap-1}), 1);
%         muVtime.post = mean(cat(1,Vlo_laps{n,r}{CSlap+1:end}), 1);

        muVtimeD{n,r}.pre = mean(Vtime{n,r}(1:CSlap-1,:),1); % downsample to have Nbin2 temporal bins
        muVtimeD{n,r}.post = mean(Vtime{n,r}(CSlap+1:end,:),1);


%         padStartT = muVtime.pre(end-sgWinLen:end); padEndT = muVtime.pre(1:sgWinLen);
%         DummyTpre = sgolayfilt([padStartT, muVtime.pre, padEndT], sgOrder,sgWinLen); % smooth using the same sgolay filter as for the spatial traces
%         DummyTpre = DummyTpre(sgWinLen+1:end-sgWinLen-1); %remove pads
%         muVtimeD{n,r}.pre = DummyTpre(1:dt2/dt:end); % downsample the smoothed voltage trace to store and plot less datapoints
%         clear padStartT padEndT DummyTpre
%         muVtime.post = mean(cat(1,Vlo_laps{n,r}{CSlap+1:end}), 1);
%         padStartT = muVtime.post(end-sgWinLen:end); padEndT = muVtime.post(1:sgWinLen);
%         DummyTpost = sgolayfilt([padStartT, muVtime.pre, padEndT], sgOrder,sgWinLen); % smooth using the same sgolay filter as for the spatial traces
%         DummyTpost = DummyTpost(sgWinLen+1:end-sgWinLen-1); %remove pads
%         muVtimeD{n,r}.post = DummyTpost(1:dt2/dt:end); % downsample the smoothed voltage trace to store and plot less datapoints
%         clear padStartT padEndT DummyTpost

%         DeltaVm_time{n,r} = muVtime{n,r}.post - muVtime{n,r}.pre; % plot against Tlap1
          DeltaVm_timeD{n,r} = muVtimeD{n,r}.post - muVtimeD{n,r}.pre; % plot against downsample time vector Tlap1(1:dt2/dt:end)
        

%             % peaks times and deltaVm at CStime
%             [PeakT{n,r}.pre, PeakTidx{n,r}.pre] = max(smoothdata(muVtime{n,r}.pre, 'sgolay'));
%             [PeakT{n,r}.post, PeakTidx{n,r}.post] = max(smoothdata(muVtime{n,r}.post, 'sgolay'));
%             
%             PeakTime{n,r}.pre = Tlap1(PeakTidx{n,r}.pre); % time of the peak during lap
%             PeakTime{n,r}.post = Tlap1(PeakTidx{n,r}.post);
%             
%             dVmPeak1T{n,r} = DeltaVm_time{n,r}(PeakTidx{n,r}.pre); % DeltaVm at peak 1
%             dVmPeak2T{n,r} = DeltaVm_time{n,r}(PeakTidx{n,r}.post); % deltaVm at the new peak
%             dVmCStime{n,r} = DeltaVm_time{n,r}(Tlap1==addCStime{n}); % deltaVm at time of CS

            % pad with NaNs in order to center on CS and, later, compute average centered on CS. 
%             NaNpadL = NaN( 1, cast((Tlap1(end)-addCStime{n})./dt, "uint16") );
%             NaNpadR = NaN(1, cast(addCStime{n}/dt, "uint16") );
%             DeltaVm_timePad{n,r} = [NaNpadL, DeltaVm_time{n,r}, NaNpadR]; % each trace is thus 2*length(Tlap1)
%             clear NaNpadL NaNpadR

%             NaNpadL = NaN( 1, cast((Tlap1(end)-addCStime{n})./dt2, "uint16") ); %divide by time step size
%             NaNpadR = NaN(1, cast(addCStime{n}/dt2, "uint16") );
%             DeltaVm_timePadD{n,r} = [NaNpadL, DeltaVm_timeD{n,r}, NaNpadR]; % each trace is thus 2*length(Tlap1)
%             clear NaNpadL NaNpadR

%             NaNpadL = NaN(1,Nbin2-addCStimebin{n});
%             NaNpadR = NaN(1,addCStimebin{n}-1);
%             DeltaVm_timePadD{n,r} = [NaNpadL, DeltaVm_timeD{n,r}, NaNpadR]; % each trace is thus 2*size(BinCenters2)-1= (2*Nbin2)-1
%             clear NaNpadL NaNpadR
%         
        % center the spatial and temporal Vm traces on CS location and time: 
        % (like in Milstein for easier comparison, with negative values corresponding to before the CS, i.e. not like in my learning rule display)
%         CScenteredLoc{n} = BinCenters2 - addCSloc{n};
        CScenteredLocR{n,r} = BinCenters2 - addCSloc{n}; % track for repetitions
        CScenteredTime{n} = Tlap1 - addCStime{n};
%         CScenteredTimeR{n,r} = Tlap1 - addCStime{n}; % track for repetitions
        CScenteredTimeD{n,r} = TimeDS - addCStime{n};

%         % prepost dW = f(time or distance from CS, initial W), i.e. equivalent to deltaVm 3D plots 
%         Wpre{n,r} = W2{n,r}(:,1);
%         Wpost{n,r} = W2{n,r}(:,end);
%         dW_prepost{n,r} = Wpost{n,r} - Wpre{n,r}; % weight change after CS induction, to plot against CScenteredTimeD or CScenteredLocR (assuming N = Nbin2)

        % effective weight changes, corresponding weights before plasticity, corresponding CS-centered times
        dWeff = diff(W2{n,r}, 1, 2); % instantaneous weight change at every time point (columns) for all synapses (rows): potentiation when there is a input spike or a CS, depression when dW changes on other inputs
%         dWeffCSlap{n,r} = dWeff(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2); % dWeff on CSlap
%         [dWeffnz_R, dWeffnz_C, dWeffnz{n,r}] = find(dWeffCSlap{n,r}); % non-zero (+ and -) values on dWeffCSlap
%         dWeffnz_Idx = find(dWeffCSlap{n,r}); % same but with linear idx
%         
%         Pcslap = P{n,r}(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2);
%         Dcslap = D{n,r}(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2);
%         InRastersCSlap = InRasters{n,r}(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2);
%         CScslap = CS{n}(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2);
% %         dWmaxSpk{n,r} = dWeffnz{n,r}./ ( Pcslap(dWeffnz_Idx)- Dcslap(dWeffnz_C)' ); % max and min inst
%         dWmaxSpk{n,r} = dWeffnz{n,r}./ ( Pcslap(dWeffnz_Idx).*InRastersCSlap(dWeffnz_Idx)+[Dcslap(dWeffnz_C).*CScslap(dWeffnz_C)*Params.PostPreBoost]' ); % estimate instant weight change per spike, CS time excluded
%         dWmaxSpk_forCSt = dWeffnz{n,r}()./ Pcslap(dWeffnz_Idx); % keeps values at CS time
% %         maxDW = max( dWeffnz{n,r}./ ( Pcslap(dWeffnz_Idx) );
%         %min(dWmaxSpk{n,r}(~isinf(dWmaxSpk{n,r})).*10^12)
% 
% 
%         W2init = W2{n,r}(:,NewLapIdx(1):NewLapIdx(2)-2); % W2 on lap1 (i.e. pre-induction plottable against Tlap1-1idx);
%         W2CSlap  = W2{n,r}(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2); 
%         W2init_nz{n,r} = W2init(dWeffnz_R,1); % initial W for each non-zero dWeff
%         W2bef{n,r} = W2CSlap(dWeffnz_Idx-1); % weight just before change, for each non-zero dWeff
% 
%         CScenteredTimeRW{n,r} = CScenteredTime{n}(dWeffnz_C); % same size as dWeffnz{n,r}       
%         
%         dWmaxSpk_CSt{n,r} = dWmaxSpk_forCSt(CScenteredTimeRW{n,r}==0); % estimate instant weight change per spike, at CS time
%         W2init_CSt{n,r} = W2init_nz{n,r}(CScenteredTimeRW{n,r}==0);

        %Effective dWmax
        [effdWmaxAll, dWmaxAllIdx] = max(dWeff, [], 2); 
        [effdWmax, dWmaxIn] = max(effdWmaxAll);
        fdWmax_time(n,r) = Trun(dWmaxAllIdx(dWmaxIn))-(Tlap1(end).*(CSlap-1)+addCStime{n}); % time from CS when the max weight change occurs
        effPamp(n,r) = effdWmax/P{n,r}(dWmaxIn,dWmaxAllIdx(dWmaxIn)); % deduce effective Pamp

%         W2bef_F(n,r) = W2{n,r}(dWmaxIn, dWmaxAllIdx(dWmaxIn)-1); % Weight just before max change for the input where max change has been detected

        clear dWeff W2init CScenteredTimeW W2CSlap effdWmaxAll dWmaxAllIdx effdWmax dWmaxIn

    end

end

% AllDeltaVm_space = cat(1, DeltaVm_space{:});
% AllCScenteredLoc = cat(1,CScenteredLocR{:});
% 
% AllDeltaVm_spacePad = cat(1, DeltaVm_spacePad{:});
% MeanDeltaVm_spacePad = mean(AllDeltaVm_spacePad,1, 'omitnan');
% AllDeltaVm_timePadD = cat(1, DeltaVm_timePadD{:});
% MeanDeltaVm_timePadD = mean(AllDeltaVm_timePadD,1, 'omitnan');

%% Correlations
clear T1 T2 T3

% distance peak 1 to CS vs distance between peaks

T1.DistPeak1toCS = [DistPeak1toCS{:}]';
T1.DistPeaks = [DistPeaks{:}]';
T1 = struct2table(T1);
LM{u} = fitlm(T1, 'DistPeaks ~ DistPeak1toCS');
% LMnoint{u} = fitlm(T1, 'DistPeaks ~ DistPeak1toCS - 1'); %no intercept
[Rdist1 Pdist1] = corrcoef(T1.DistPeak1toCS, T1.DistPeaks);

T2.DistPeak1toCS = T1.DistPeak1toCS( T1.DistPeak1toCS < 50); % keeping only CSs happening ~close to peak1 (50cm for Milstein params)
T2.DistPeaks = T1.DistPeaks(T1.DistPeak1toCS < 50);
T2 = struct2table(T2);
LMclose{u} = fitlm(T2, 'DistPeaks ~ DistPeak1toCS');
[Rdist2 Pdist2] = corrcoef(T2.DistPeak1toCS, T2.DistPeaks);
% LMclose_noint{u} = fitlm(T2, 'DistPeaks ~ DistPeak1toCS-1');

beta(:,u) = LM{u}.Coefficients.Estimate;
Rsquared(u) = LM{u}.Rsquared.Ordinary;
[pval(u),F(u),r(u)] = coefTest(LM{u}); % p-val, F value, numerator degrees of freedom r

beta2(:,u) = LMclose{u}.Coefficients.Estimate;
Rsquared2(u) = LMclose{u}.Rsquared.Ordinary;
[pval2(u),F2(u),r2(u)] = coefTest(LMclose{u});

R2dist1(u) = Rdist1(1,2)^2;
R2dist2(u) = Rdist2(1,2)^2;
Pvaldist1(u) = Pdist1(1,2);
Pvaldist2(u) = Pdist2(1,2);

% figure % distance peak 1 to CS vs distance between peaks
% scatter([DistPeak1toCS{:}], [DistPeaks{:}], 'k', 'filled'); hold on
% plot([0 L/2], [0 L/2], 'r--'); % identity line
% xlabel('Distance between Peak1 and CS (cm)');
% ylabel('Distance between Peaks (cm)')
% title(['Pamp:' num2str(PampRange(u)) ', R2 = ' num2str(LM.Rsquared.Ordinary*100) '%, R2noint = ' num2str(LMnoint.Rsquared.Ordinary*100) '%'] )
% box off; axis square
% 
% figure
% plot(LM)
% 
% figure
% plot(LMnoint)

% initial Vm vs delta Vm or FinalVm (using spatially binned delta Vm, to have less points)
AllmuVbin = [muVbin{:}]; % [AllmuVbin.pre_smooth] is a row vector
T3.AllDeltaVspace = [DeltaVm_space{:}].*10^3;
T3.initialVm = [AllmuVbin.pre_smooth].*10^3;
T3.finalVm = [AllmuVbin.post_smooth].*10^3;
T3 = struct2table(T3);
[Rid,Pid] = corrcoef(T3.initialVm, T3.AllDeltaVspace);
[Rid2,Pid2] = corrcoef(T3.initialVm, T3.finalVm);
RsquaredT3delta(u) = Rid(1,2)^2;
RsquaredT3final(u) = Rid2(1,2)^2;
PvalT3delta(u) = Pid(1,2);
PvalT3final(u) = Pid2(1,2);

% figure
% scatter(T3.initialVm,T3.AllDeltaVspace,'.k')
% 
% figure
% scatter(T3.initialVm,T3.finalVm,'.k')

% Time from CS vs pre-induction Vm vs deltaVm (Downsampled)
AllmuVtime = [muVtimeD{:}];
SurfX2 = [CScenteredTimeD{:}]; SurfY2 = ([AllmuVtime.pre]-Params.Vrest).*10^3; SurfZ2 = [DeltaVm_timeD{:}].*10^3;
[surface2, Sgoodness2, Soutput2] = fit([SurfX2' SurfY2'], SurfZ2','linearinterp'); % note that 'lowess' method would give a smooth surface + go outside defined bounds, but takes too much time to run. Need to do that on downsampled data.
f2 = plot(surface2);
MaxInterp_dVm(u) = max(f2.ZData(:));
MaxDeltaVm(u) = max(SurfZ2);
MaxEffPamp(u) = max(effPamp(:)).*10^12;

% clearvars -except MaxEffPamp MaxDeltaVm MaxInterp_dVm LM LMnoint LMclose LMclose_noint mdl mdl2

end
%% Summary Figures across multiple simulations of the above experiment, using different Params.Pamp
PampRangepA = PampRange.*10^12;
PampRangeU = unique(PampRangepA);

PampRangeIdx = zeros(size(PampRangepA));
maxFampMeans = zeros(size(PampRangeU));
MeanR2 = zeros(size(PampRangeU));
MeanSlope = zeros(size(PampRangeU));
MeanSlopeClose = zeros(size(PampRangeU));

for n = 1:length(PampRangeU)
    idx = find(PampRangepA == PampRangeU(n));
    maxFampMeans(n) = mean(MaxEffPamp(idx));
    MeanR2(n) = mean(Rsquared(idx));
    MeanSlope(n) = mean(beta(2,idx));
    MeanSlopeClose(n) = mean(beta2(2,idx));
    MeanR2T3final(n) = mean(RsquaredT3final(idx));
    MeanPvalT3final(n) = mean(PvalT3final(idx));
    PampRangeIdx(idx) = n;
end

% dWmax and dVm_max as a function of Params.Pamp

Ynl1 = MaxEffPamp';
Xnl1 = PampRange'.*10^12;                
% modelfun = fittype( @(p1,p2,x) p1*(1-exp(-x/p2)) );
modelfun = fittype(@(a,b,x) a+b.*log(x) );
options = fitoptions('Method','NonlinearLeastSquares', 'Algorithm', 'Trust-Region');
options.Startpoint = [1 2];
[mdl1,gof,foutput] = fit(Xnl1,Ynl1,modelfun,options);
x1 = 1:0.01:max(PampRange.*10^12,[], 'all')+1;
% Pcoef = polyfit(Xpoly1,Ypoly1,6);
% Poly1 = polyval(Pcoef,x1);

% figure
toptim = tiledlayout(1,4, 'TileSpacing','Compact','Padding','Compact');
title(toptim,'Milstein Params: 185cm track, animal speed: 25 cm/s')

nexttile
% subplot(1,4,1)
    % plot(x1, feval(mdl1,x1), 'r--');hold on
    plot(PampRangeU, maxFampMeans, 'r--'); hold on
    scatter(PampRangepA, MaxEffPamp, 'ok', 'filled', 'MarkerFaceAlpha', 0.4); hold on
    xlim([0 max(PampRange.*10^12)+1])
    ylim([0 max(MaxEffPamp,[], 'all')+1])
    xlabel('Params.Pamp (pA)')
    ylabel('effective dWmax estimate (pA)')
%     title('185cm track, animal speed: 25 cm/s')
    axis square
    box off

nexttile
% subplot(1,4,2) % R2 for distances corr
    Ynl2 = Rsquared'*100;               
    % modelfun2 = fittype( @(k,mp,s,x) k.*( 1./( 1+( ( x.*(1-mp) )./( mp.*(1-x) ) ).^s ) ) ) ;
    modelfun2 = fittype( @(k,mp,s,x) k.*(tanh((x-mp)/s)/2+0.5) )  ;
    options2 = fitoptions('Method','NonlinearLeastSquares', 'Algorithm', 'Trust-Region');
    % options2.Startpoint = [90 15 3];
    options2.Startpoint = [90 15 3];
    [mdl2,gof,foutput] = fit(Xnl1,Ynl2,modelfun2,options2);
    fv2 = feval(mdl2,x1);
    [~, R50] = min(abs(fv2-50));
    
    plot(x1, fv2, 'r-'); hold on
    scatter(PampRangepA, Rsquared*100, 'ok', 'filled', 'MarkerFaceAlpha', 0.4); hold on
    plot([30 30], [0 MeanR2(PampRangeU == 30)*100], 'k--');
    plot([0 30], [MeanR2(PampRangeU == 30)*100 MeanR2(PampRangeU == 30)*100], 'k--')
    % plot([mdl2.mp mdl2.mp], [0 feval(mdl2,mdl2.mp)], 'k--');
    % plot([0 mdl2.mp], [feval(mdl2,mdl2.mp) feval(mdl2,mdl2.mp)], 'k--')
%     plot([x1(R50) x1(R50)], [0 feval(mdl2,x1(R50))], 'k--');
%     plot([0 x1(R50)], [feval(mdl2,x1(R50)) feval(mdl2,x1(R50))], 'k--')
    xlim([0 max(PampRange.*10^12)+1])
    % ylim([0 100])
    xlabel('Params.Pamp (pA)')
    ylabel('Distance correlations, R2 (%)')
%     title('300cm track, animal speed: 15 cm/s')
    axis square
    box off

nexttile
% subplot(1,4,3) % slope for distance corr, model with and without intercept
    FO = fit(PampRange'.*10^12, beta(2,:)', 'smoothingspline');
    scatter(PampRange.*10^12, beta(2,:), 'ok', 'filled', 'MarkerFaceAlpha', 0.4); hold on % scatter(PampRange.*10^12, beta(2,:), 'o', 'filled', 'MarkerFaceColor', [0.5 0.2 0.5], 'MarkerFaceAlpha', 0.5); hold on
    plot(FO); legend off
    % plot(PampRangeU,MeanSlope, '--', 'Color', [0.5 0.2 0.5]); hold on
    % % scatter(PampRange.*10^12, beta2(2,:), 'o', 'filled', 'MarkerFaceColor', [1 0.7 0.1], 'MarkerFaceAlpha', 0.5); hold on
    % % plot(PampRangeU,MeanSlopeClose, '--', 'Color', [1 0.7 0.1]); hold on
    % % legend('all Peak1-CS distances', 'mean', 'DistPeak1toCS < 50cm', 'mean',  'Location', 'Best')
    xlim([min(PampRange.*10^12)-0.1 max(PampRange.*10^12)+1])
    xlabel('Params.Pamp (pA)')
    ylabel('Distance correlations: slope')
%     title('300cm track, animal speed: 15 cm/s')
    axis square
    box off

nexttile
% subplot(1,4,4) % correlation between initial and final Vm 
    scatter(PampRangepA, RsquaredT3final*100, 'ok', 'filled', 'MarkerFaceAlpha', 0.4); hold on
    plot(PampRangeU, MeanR2T3final*100, 'r--'); hold on
    xlim([0 max(PampRange.*10^12)+1])
    ylim([0 100])
    xlabel('Params.Pamp (pA)')
    ylabel('Pre vs Post-induction Vm corr., R2 (%)')
%     title('300cm track, animal speed: 15 cm/s')
    axis square
    box off

figure % correlation between initial and final Vm 
scatter(PampRangepA, PvalT3final, 'ok', 'filled', 'MarkerFaceAlpha', 0.4); hold on
plot(PampRangeU, MeanPvalT3final, 'r--'); hold on
% yline(0.05,'k--')
xlim([0 max(PampRange.*10^12)+1])
xlabel('Params.Pamp (pA)')
ylabel('Pre vs Post-induction Vm correlation, P-value')
% title('300cm track, animal speed: 15 cm/s')
axis square
box off
%InVivoMilstein2021_1CSbatch
clear
close all

%% Parameters
% Inputs params
N = 65 %100; % number of input neurons
L = 185 %300 % length of track in cm (300 for Can, 185 for Bittner and Milstein). 
PFsd = 18; % cm (~median in Can data) 
PFamp = 10; %peak FR in Hz
Nlaps = 23; % 30; % number of laps
speed = 25; %15; % cm/sec, (15cm/s average in DongSheffield2021, 25cm/s in Milstein2021 network model)
dt = 0.001; % time resolution in sec

% Synapses params
Imax = 85e-12; % max synaptic weight, a current in Amps. % NolanWiles2011 uses 10-226pA. Mehta uses 1.5-4nA. SongAbbot use dimensionless conductances with G between 0 and Gmax, but not sure about the initial weight matrix. 
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

%Milstein-model biBTSTP rule params (average of fits across all Bittner+Milstein recordings)
% Params.Pdecay = 1.66; % Synaptic Eligibility Trace decay time constant, in seconds (average: 0.9; Milstein Fig 5: 2.5s; for Network: 1.66s) 
% Params.Ddecay = 0.44; % Global Instructive Signal (from CS) decay time constant, in seconds (average: 0.5; Milstein Fig 5: 1.5s; for Network: 0.44s) 
% Params.k_pot = 1.1; % in s-1 (average: 2.27; Milstein Fig 5: 1.7; for Network: 1.1) 
% Params.k_dep = 0.7; % in s-1 (average: 0.33; Milstein Fig 5: 0.204; for Network: 0.424) 
% Params.a_pot = 0.415; % sigmoid threshold for potentiation gain (average: 0.24; Milstein Fig 5: 0.5; for Network: 0.415) 
% Params.b_pot = 4.4; % sigmoid slope for potentiation gain (average: 30.32; Milstein Fig 5: 4; for Network: 4.40) 
% Params.a_dep = 0.026; % thresh for depression gain (average: 0.09; Milstein Fig 5: 0.01; for Network: 0.026) 
% Params.b_dep = 20.04; % slope for depression gain (average: 2260; Milstein Fig 5: 44.44; for Network: 20.04) 
% Params.Wmax = 1.5*Imax; % 57e-12; 
Params.Pdecay = 2; % Synaptic Eligibility Trace decay time constant, in seconds (average: 0.9; Milstein Fig 5: 2.5s; for Network: 1.66s) 
Params.Ddecay = 1.5; % Global Instructive Signal (from CS) decay time constant, in seconds (average: 0.5; Milstein Fig 5: 1.5s; for Network: 0.44s) 
Params.k_pot = 1.2; % in s-1 (average: 2.27; Milstein Fig 5: 1.7; for Network: 1.1) 
Params.k_dep = 0.2; % in s-1 (average: 0.33; Milstein Fig 5: 0.204; for Network: 0.424) 
Params.a_pot = 0.5; % sigmoid threshold for potentiation gain (average: 0.24; Milstein Fig 5: 0.5; for Network: 0.415) 
Params.b_pot = 4; % sigmoid slope for potentiation gain (average: 30.32; Milstein Fig 5: 4; for Network: 4.40) 
Params.a_dep = 0.01; % thresh for depression gain (average: 0.09; Milstein Fig 5: 0.01; for Network: 0.026) 
Params.b_dep = 44; % slope for depression gain (average: 2260; Milstein Fig 5: 44.44; for Network: 20.04) 
Params.Wmax = Imax; % 57e-12;

% BTSP-induction params
Params.CSproba = []; % proba of an output spike to be a plateau potential/complex spike (1.8% on average in Bittner et al. 2017)
CSlap = 11%:13; % lap on which CS will occur if addCS=1
repeats = 1; % number of repeats with same CS time

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

% segment in laps and spatial bins
NewLapIdx = find(Run==0);
NewLapTimes = Trun(NewLapIdx);
EndLapIdx = find(Run==L);
laps_end = zeros(size(Run));
laps_end(EndLapIdx) = 1;

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

% design FIR for smoothing Vm: with wrap-around padding of Hamming window size and zero-phase filtering
Fcut = 3; %cut-off frequency, in Hz
Fs = 1./dt; %sampling frequency in Hz
NyqF = Fs/2; % Nyquist frequency
HWin = 2; % size of hamming window, in seconds
fir_n = HWin.*Fs-1; % filter order
nFcut = Fcut/NyqF; % normalized cutoff frequency (proprtion of the Nyquist frequency)
blo = fir1(fir_n,nFcut); % design a FIR filter of fir_n'th order, using hamming window, with cut-off frequency at 3Hz
%         freqz(blo,1,2^16,Fs)

%% input-output dynamics

for n = 1:3:N % vary location of CS along COM of input neurons. Use range that matches Milstein 2021 (0 to 140cm between CS and PF. Actually, the 140cm extrema correspond to initial peak near edge of track. I cannot go that far. My range is 0 to 185/2=92.5cm if I use a 185cm track like them. But I could use a longer track... 
    
    % make CS trace with 1 CS at time corresponding to COM of input n
    CS{n} = zeros(1,length(Trun)); 
    plateau = CS;
    [minvalN(n), COMidxN(n)] = min(abs(Lap1-PFcom(n)));
    COMidxN_Laps = find(Run==Run(COMidxN(n))); % COM idx of input cell n on all laps (assuming constant speed)
    CS{n}(COMidxN_Laps(CSlap)) = 1;
    for l = 1:length(CSlap)
    plateau{n}(COMidxN_Laps(CSlap(l)):COMidxN_Laps(CSlap(l))+300) = 1; % 300ms plateau
    end
    
    for r = 1:repeats
        
        % nonhomogenous Poisson spike generator along trajectory
        InRasters{n,r} = poissrnd(FRrun.*dt); %matrix of size (N,length(Tlap1)) containing, for each input neuron, the number of spikes in each 1ms bin drawn from Poisson distribution
        
        % run LIF + plasticity model
        Params.dt = dt; Params.Imax = Imax;
        [ ~, ~, V2{n,r}, ~, OutRaster{n,r}, D{n,r}, P{n,r}, ~, W2{n,r}, ~]  = MilsteinBTSP_LIFadapt_ODEs(Params, W, InRasters{n,r}, plateau{n}, laps_end);

%         Spiketimes_outAll = Trun(logical(OutRaster{n,r}));
%         SpikeLoc_outAll = Run(logical(OutRaster{n,r})); %spikes location on track
%         CSloc_All = Run(logical(CS{n})); % CSs location on track 
        
        % smooth Vm traces by using a 3Hz low-pass filter like in Bittner2017     
        Vpad = [V2{n,r}(length(V2{n,r})-HWin*Fs:end), V2{n,r}, V2{n,r}(1:HWin*Fs)]; % wrap-around pad
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
        
        % spatial deltaVm (difference after induction)
        DeltaVm_space{n,r} = muVbin{n,r}.post_smooth-muVbin{n,r}.pre_smooth; % plot against BinCenters2
        
            dVmPeak1{n,r} = DeltaVm_space{n,r}(PeakIdx{n,r}.pre); % DeltaVm at peak 1 position
            dVmPeak2{n,r} = DeltaVm_space{n,r}(PeakIdx{n,r}.post); % deltaVm at the new peak position
            dVmCS{n,r} = DeltaVm_space{n,r}(addCSbin2{n}); % deltaVm at location of CS

            % pad with NaNs in order to center on CS and, later, compute average centered on CS. 
            NaNpadL = NaN(1,Nbin2-addCSbin2{n});
            NaNpadR = NaN(1,addCSbin2{n}-1);
            DeltaVm_spacePad{n,r} = [NaNpadL, DeltaVm_space{n,r}, NaNpadR]; % each trace is thus 2*size(BinCenters2)-1= (2*Nbin2)-1
            clear NaNpadL NaNpadR
        
        % temporal pre, post and deltaVm (smooth and downsample)
        muVtime.pre = mean(cat(1, Vlo_laps{n,r}{1:CSlap-1}), 1);
        muVtime.post = mean(cat(1,Vlo_laps{n,r}{CSlap+1:end}), 1);

        muVtimeD{n,r}.pre = mean(Vtime{n,r}(1:CSlap-1,:),1); % using Vtime, a downsampled version to have Nbin2 temporal bins
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

            NaNpadL = NaN(1,Nbin2-addCStimebin{n});
            NaNpadR = NaN(1,addCStimebin{n}-1);
            DeltaVm_timePadD{n,r} = [NaNpadL, DeltaVm_timeD{n,r}, NaNpadR]; % each trace is thus 2*size(BinCenters2)-1= (2*Nbin2)-1
            clear NaNpadL NaNpadR
        
        % center the spatial and temporal Vm traces on CS location and time: 
        % (like in Milstein for easier comparison, with negative values corresponding to before the CS, i.e. not like in my learning rule display)
%         CScenteredLoc{n} = BinCenters2 - addCSloc{n};
        CScenteredLocR{n,r} = BinCenters2 - addCSloc{n}; % track for repetitions
        CScenteredTime{n} = Tlap1 - addCStime{n};
%         CScenteredTimeR{n,r} = Tlap1 - addCStime{n}; % track for repetitions
        CScenteredTimeD{n,r} = TimeDS - addCStime{n};

        % prepost dW = f(time or distance from CS, initial W), i.e. equivalent to deltaVm 3D plots 
        Wpre{n,r} = W2{n,r}(:,1);
        Wpost{n,r} = W2{n,r}(:,end);
        dW_prepost{n,r} = Wpost{n,r} - Wpre{n,r}; % weight change after CS induction, to plot against CScenteredTimeD or CScenteredLocR (assuming N = Nbin2)

%         % effective weight changes, corresponding weights before plasticity, corresponding CS-centered times
%         dWeff = diff(W2{n,r}, 1, 2); % instantaneous weight change at every time point (columns) for all synapses (rows): potentiation when there is a input spike or a CS, depression when dW changes on other inputs
%         dWeffCSlap{n,r} = dWeff(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2); % dWeff on CSlap
%         [dWeffnz_R, dWeffnz_C, dWeffnz{n,r}] = find(dWeffCSlap{n,r}); % non-zero (+ and -) values on dWeffCSlap
%         dWeffnz_Idx = find(dWeffCSlap{n,r}); % same but with linear idx
% 
%         Pcslap = P{n,r}(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2);
%         Dcslap = D{n,r}(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2);
%         InRastersCSlap = InRasters{n,r}(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2);
%         CScslap = CS{n}(:,NewLapIdx(CSlap):NewLapIdx(CSlap+1)-2);
% %         dWmaxSpk{n,r} = dWeffnz{n,r}./ ( Pcslap(dWeffnz_Idx)- Dcslap(dWeffnz_C)' ); % max and min inst
%         % dWmaxSpk{n,r} = dWeffnz{n,r}./ ( Pcslap(dWeffnz_Idx).*InRastersCSlap(dWeffnz_Idx)+[Dcslap(dWeffnz_C).*CScslap(dWeffnz_C)*Params.PostPreBoost]' ); % estimate instant weight change per spike, CS time excluded
%         % dWmaxSpk_forCSt = dWeffnz{n,r}()./ Pcslap(dWeffnz_Idx); % keeps values at CS time
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
%         % dWmaxSpk_CSt{n,r} = dWmaxSpk_forCSt(CScenteredTimeRW{n,r}==0); % estimate instant weight change per spike, at CS time
%         W2init_CSt{n,r} = W2init_nz{n,r}(CScenteredTimeRW{n,r}==0);
% 
%         %Effective dWmax
%         [effdWmaxAll, dWmaxAllIdx] = max(dWeff, [], 2); 
%         [effdWmax, dWmaxIn] = max(effdWmaxAll);
%         fdWmax_time(n,r) = Trun(dWmaxAllIdx(dWmaxIn))-(Tlap1(end).*(CSlap-1)+addCStime{n}); % time from CS when the max weight change occurs
%         effPamp(n,r) = effdWmax/P{n,r}(dWmaxIn,dWmaxAllIdx(dWmaxIn)); % deduce effective Pamp
% 
%         W2bef_F(n,r) = W2{n,r}(dWmaxIn, dWmaxAllIdx(dWmaxIn)-1); % Weight just before max change for the input where max change has been detected

        clear dWeff W2init CScenteredTimeW W2CSlap effdWmaxAll dWmaxAllIdx effdWmax dWmaxIn

    end

end

AllDeltaVm_space = cat(1, DeltaVm_space{:});
AllCScenteredLoc = cat(1,CScenteredLocR{:});

AllDeltaVm_spacePad = cat(1, DeltaVm_spacePad{:});
MeanDeltaVm_spacePad = mean(AllDeltaVm_spacePad,1, 'omitnan');
% AllDeltaVm_timePad = cat(1, DeltaVm_timePad{:});
% MeanDeltaVm_timePad = mean(AllDeltaVm_timePad,1, 'omitnan');
AllDeltaVm_timePadD = cat(1, DeltaVm_timePadD{:});
MeanDeltaVm_timePadD = mean(AllDeltaVm_timePadD,1, 'omitnan');

%% figures

% figure %BTSP variables on given lap
LapNum = CSlap;
InputNum = 50 %dWmaxIn{n,r}; %N/2; %Nbittner
x1 = 0; % in s
x2 = period; % in s 

% subplot(4,1,1)
%     plot(Tlap1, D{n,r}(NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'r'); hold on
%     plot(repmat(Spiketimes_out{n,r}{LapNum},2,1), repmat([0; 1], 1, length(Spiketimes_out{n,r}{LapNum})), 'k'); hold on
%     scatter(Tlap1(addCSidx{n}), ones(size(addCSidx{n})).*0.5, 'oc', 'filled')
%     xlim([x1 x2]);
%     xlabel('time (s)'); 
%     title('output raster, CS, Post-before-Pre variable')
%     box off;
% subplot(4,1,2)
%     plot(Tlap1, max(P{n,r}(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)).*InRasters{n,r}(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'k'); hold on % input spiketrain
%     plot(Tlap1, P{n,r}(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'r'); hold on
%     xlim([x1 x2]);
%     xlabel('time (s)'); 
%     title(['input raster #' num2str(InputNum) '+ Pre-before-Post variable'])
%     box off;
% subplot(4,1,3)
% %     plot(Tlap1, Wraw{n,r}(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2{n,r}(:)), 'b--'); hold on
%     plot(Tlap1, W2{n,r}(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./max(W2{n,r}(:)), 'r'); hold on
%     xlim([x1 x2]);
%     xlabel('time (s)'); 
%     ylabel('Norm. synaptic strength')
%     title(['input #' num2str(InputNum)] )
%     box off;
% % subplot(4,1,4)    
% %     plot(Tlap1(1:end-1), dWeffCSlap{n,r}(InputNum,:).*10^12, 'g'); hold on
% %     xlim([x1 x2]);
% %     xlabel('time (s)'); 
% %     ylabel('dWeff (pA)')
% %     title(['input #' num2str(InputNum)] )
% %     box off;

% example
dummy = 1:2:N;
n = dummy(5);
r = 1;

clear T1 T3
figure
t1 = tiledlayout('flow', 'TileSpacing','Compact','Padding','Compact')
% nexttile
%     % f = figure; % pre and post spatial Vm for given {n,r}
%     plot(BinCenters2, muVbin{n,r}.pre_smooth.*10^3, 'k-'); hold on
%     plot(BinCenters2, muVbin{n,r}.post_smooth.*10^3, 'b-');
%     legend('pre', 'post')
%     xlabel('position'); ylabel('membrane potential (mV)')
%     xline(addCSloc{n}, 'r')
%     axis square
%     box off;
%     title('pre and post spatial Vm');
%     % f.Renderer = 'painters'   
%     % print('BTSP_PotentiationRules.pdf','-dpdf')
% 
% nexttile % CS-centered spatial deltaVm for given {n,r}
%     plot(CScenteredLocR{n,r}, DeltaVm_space{n,r}.*10^3, 'k-'); hold on
%     scatter(CScenteredLocR{n,r}(PeakIdx{n,r}.pre), dVmPeak1{n,r}.*10^3, 'k^', 'filled');hold on
%     scatter(CScenteredLocR{n,r}(PeakIdx{n,r}.post), dVmPeak2{n,r}.*10^3, 'b^', 'filled');hold on
%     xline(0, 'r')
%     xlabel('distance from CS (cm)'); ylabel('membrane potential (mV)')
%     xlim([floor(min(CScenteredLocR{n,r})) ceil(max(CScenteredLocR{n,r}))]);
%     axis square
%     box off;
%     title('spatial delta Vm (post-pre induction) '); 

% summary figures

nexttile % position of peak 1, CS and peak 2 for all n and r (like fig1D in Milstein2021)
    % add xlines to show the range of the PF before induction
    PeakLocAll = [PeakLoc{:}]; PeaksLocPre = [PeakLocAll.pre]; PeaksLocPost = [PeakLocAll.post]; 
    [addCSlocRall, sortIdx] = sort([addCSlocR{:}]);
    scatter(PeaksLocPre(sortIdx), 1:size([PeakLocAll.pre],2), 'k|', 'LineWidth', 2);hold on
    scatter(addCSlocRall, 1:size([PeakLocAll.pre],2), 'vr', 'filled', 'MarkerFaceAlpha', 0.5);hold on
    scatter(PeaksLocPost(sortIdx), 1:size([PeakLocAll.pre],2), 'b^', 'filled','MarkerFaceAlpha', 0.5);hold on
    legend('Peak1', 'CS', 'Peak2', 'Location', 'Northwest')
    xlim([0 L])
    xlabel('position on track (cm)')
    ylabel('simulated cell index')
    axis square; box off

nexttile % CS-centered DeltaVm_time profiles + mean
    CScenteredTds = linspace(-TimeDS(end),TimeDS(end), 2*Nbin2-1);
    plot(CScenteredTds,AllDeltaVm_timePadD.*10^3,'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5); hold on
    plot(CScenteredTds,MeanDeltaVm_timePadD.*10^3, 'k', 'LineWidth', 2); hold on
    xline(0, 'r')
    xlim([-Tlap1(end) Tlap1(end)])
    xlabel('time from CS (s)'); ylabel('membrane potential (mV)')
    box off; axis square
    title('temporal delta Vm (post-pre induction) ');

nexttile % CS-centered DeltaVm_space profiles + mean
    CScenteredPos = -(Nbin2-1)*L/Nbin2:L/Nbin2:(Nbin2-1)*L/Nbin2;
    plot(CScenteredPos,AllDeltaVm_spacePad.*10^3,'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5); hold on
    plot(CScenteredPos,MeanDeltaVm_spacePad.*10^3, 'k', 'LineWidth', 2); hold on
    xline(0, 'r')
    xlim([-L L])
    xlabel('distance from CS (cm)'); ylabel('membrane potential (mV)')
    box off; axis square
    title('spatial delta Vm (post-pre induction) ');

nexttile % distance peak 1 to CS vs distance between peaks
    T1.DistPeak1toCS = [DistPeak1toCS{:}]';
    T1.DistPeaks = [DistPeaks{:}]';
    T1 = struct2table(T1);
    LM = fitlm(T1, 'DistPeaks ~ DistPeak1toCS');
    Rsquared = LM.Rsquared.Ordinary;
    [Rdist1 Pdist1] = corrcoef(T1.DistPeak1toCS, T1.DistPeaks);
    scatter([DistPeak1toCS{:}], [DistPeaks{:}], 'k', 'filled'); hold on
    plot([0 L/2], [0 L/2], 'r--'); % identity line, to regress to
    xlabel('Distance between Peak1 and CS (cm)');
    ylabel('Distance between Peaks (cm)')
    title(['R2 = ' num2str(Rsquared)])
    box off; axis square

nexttile % scatter plot: initial Vm vs delta Vm (using spatially binned delta Vm, to have less points)
    AllmuVbin = [muVbin{:}]; % [AllmuVbin.pre_smooth] is a row vector
    AllDeltaVspace = [DeltaVm_space{:}];
    T3.AllDeltaVspace = [DeltaVm_space{:}].*10^3;
    T3.initialVm = [AllmuVbin.pre_smooth].*10^3;
    T3.finalVm = [AllmuVbin.post_smooth].*10^3;
    T3 = struct2table(T3);
    [Rid,Pid] = corrcoef(T3.initialVm, T3.AllDeltaVspace);
    [Rid2,Pid2] = corrcoef(T3.initialVm, T3.finalVm);
    RsquaredT3delta = Rid(1,2)^2;
    RsquaredT3final = Rid2(1,2)^2;
    scatter([AllmuVbin.pre_smooth].*10^3, AllDeltaVspace.*10^3, '.k')
    xlabel('pre-induction Vm (mV)'); ylabel('deltaVm (mV)')
    xlim([Params.Vrest.*10^3, Params.Vthr.*10^3])
    ylim([-10 10])
    box off; axis square
%     title('Spatial (smoothed)');
    title(['R2 = ' num2str(RsquaredT3delta)])

% scatter plot: initial Vm vs final Vm (using spatially binned delta Vm, to have less points)
nexttile 
        scatter([AllmuVbin.pre_smooth].*10^3, [AllmuVbin.post_smooth].*10^3, '.k')
        xlabel('pre-induction Vm (mV)'); ylabel('post-induction Vm (mV)')
        xlim([Params.Vrest.*10^3, Params.Vthr.*10^3])
        ylim([Params.Vrest.*10^3, Params.Vthr.*10^3])
        box off; axis square
%         title('Spatial (smoothed)');
        title(['R2 = ' num2str(RsquaredT3final)])

exportgraphics(t1,'BTSP_MilsteinModel_1CSBatch_SheffieldParams_A.pdf','BackgroundColor','none', 'ContentType','vector')

figure
t2 = tiledlayout('flow', 'TileSpacing','Compact','Padding','Compact');
nexttile % 3D plot: Time from CS vs pre-induction Vm vs deltaVm (Downsampled)
    AllmuVtime = [muVtimeD{:}];
    scatter3([CScenteredTimeD{:}], ([AllmuVtime.pre]-Params.Vrest).*10^3, [DeltaVm_timeD{:}].*10^3,  '.k')
    xlabel('Time from CS (cm)'); 
    ylabel('pre-induction Vm (mV)')
    zlabel('deltaVm (mV)')
    axis square; box off
    view(-150, 14)
nexttile
    SurfX2 = [CScenteredTimeD{:}]; SurfY2 = ([AllmuVtime.pre]-Params.Vrest).*10^3; SurfZ2 = [DeltaVm_timeD{:}].*10^3;
    [surface2, Sgoodness2, Soutput2] = fit([SurfX2' SurfY2'], SurfZ2','linearinterp'); % note that 'lowess' method would give a smooth surface + go outside defined bounds, but takes too much time to run. Need to do that on downsampled data.
    f2 = plot(surface2); f2.EdgeColor = 'none';
    xlabel('Time from CS (s)'); 
    ylabel('initial Vm ramp (mV)')
    zlabel('deltaVm (mV)')
    axis square; box off
    view(-150, 14)
    colormap(brewermap(256,'*RdBu'));
nexttile % same as above but flat map, like in Milstein2021: Distance from CS vs pre-induction Vm vs deltaVm
    colormap(brewermap(256,'*RdBu'));
    f2 = plot(surface2); f2.EdgeColor = 'none'; hold on
    xline(0,'k--', 'LineWidth', 2)
    xlabel('Time from CS (s)'); 
    ylabel('initial Vm ramp (mV)')
    zlabel('deltaVm (mV)')
    axis square; box off; grid off
    view(0, 90) % flat matrix view as in Milstein 2021
    caxis([-max(f2.ZData(:)), max(f2.ZData(:))]); xlim([-4.9, 4.9])
    c = colorbar; c.Label.String = '\Delta Vm (mV)';
    title([ 'max interpolated \DeltaVm =' num2str(max(f2.ZData(:))) ' mV'])
% c = colorbar; c.Label.String = '\Delta Vm';

exportgraphics(t2,'BTSP_MilsteinModel__1CSBatch_SheffieldParams.pdf','BackgroundColor','none', 'ContentType','vector')


% figure % scatter plot: initial Vm vs delta Vm (using spatially binned delta Vm, to have less points)
%     AllmuVbin = [muVbin{:}]; % [AllmuVbin.pre_smooth] is a row vector
%     AllDeltaVspace = [DeltaVm_space{:}];
%     subplot(2,1,1)
%     scatter([AllmuVbin.pre_smooth].*10^3, AllDeltaVspace.*10^3, '.k')
%     xlabel('pre-induction Vm (mV)'); ylabel('deltaVm (mV)')
%     xlim([Params.Vrest.*10^3, Params.Vthr.*10^3])
%     ylim([-10 10])
%     box off; axis square
%     title('Spatial (smoothed)');
%     subplot(2,1,2) % as a 2D distribution
%     histogram2([AllmuVbin.pre_smooth].*10^3, AllDeltaVspace.*10^3, 'DisplayStyle','tile', 'Normalization', 'probability');
%     xlabel('pre-induction Vm (mV)'); ylabel('deltaVm (mV)')
%     xlim([Params.Vrest.*10^3, Params.Vthr.*10^3])
%     ylim([-10 10])
%     box off; axis square;  grid off
%     cb = colorbar;
%     cb.Label.String = 'probability';
% 
% figure % scatter plot: initial Vm vs final Vm (using spatially binned delta Vm, to have less points)
%     subplot(2,1,1)
%         scatter([AllmuVbin.pre_smooth].*10^3, [AllmuVbin.post_smooth].*10^3, '.k')
%         xlabel('pre-induction Vm (mV)'); ylabel('post-induction Vm (mV)')
%         xlim([Params.Vrest.*10^3, Params.Vthr.*10^3])
%         ylim([Params.Vrest.*10^3, Params.Vthr.*10^3])
%         box off; axis square
%         title('Spatial (smoothed)');
%     subplot(2,1,2) % as a 2D distribution
%         histogram2([AllmuVbin.pre_smooth].*10^3, [AllmuVbin.post_smooth].*10^3, 'DisplayStyle','tile', 'Normalization', 'probability');
%         xlabel('pre-induction Vm (mV)'); ylabel('post-induction Vm (mV)')
%         xlim([Params.Vrest.*10^3, Params.Vthr.*10^3])
%         ylim([Params.Vrest.*10^3, Params.Vthr.*10^3])
%         box off; axis square; grid off
%         cb = colorbar;
%         cb.Label.String = 'probability';

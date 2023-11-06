function output = STDP_LIF(params)
%
%params: struct with the following fields
% 
%% Parameters
%Params for inputs
N = params.N; % number of input neurons
L = params.L; % length of track in cm
PFsd = params.PFsd; % cm (~median in Can data) 
PFamp = params.PFamp; %peak FR in Hz
Nlaps = params.Nlaps; % number of laps
period = params.period; % lap duration, in sec
dt = params.dt; % time resolution in sec

% Synapses params
Imax = params.Imax; % max synaptic weight, a current in Amps. % NolanWiles2011 uses 10-226pA. Mehta uses 1.5-4nA. SongAbbot use dimensionless conductances with G between 0 and Gmax, but not sure about the initial weight matrix. 
Idecay = params.Idecay; % EPSC time constant, in sec
Wsd = params.Wsd; %standard deviation of initial synaptic weight vector (for gaussian connectivity)
maxW = params.maxW; % from 0 to 1, proportion of Imax defining the initial maximum weight

%Integrate and Fire params for output neuron
Rm = params.Rm; % membrane resistance, in Ohm
tau_m = params.tau_m ; %membrane time constant, in sec (25ms in YuShouval2006)
Vrest = params.Vrest ; % resting membrane potential, in Volt
Vthr = params.Vthr; % spike threshold, in Volt
Vreset = params.Vreset; % reset potential after spike, in Volt

% STDP params
Pdecay = params.Pdecay; % Pre-before-Post time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Ddecay = params.Ddecay; % Post-before-Pre time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Pamp = params.Pamp*Imax; % peak Pre-before-Post weight change, in percent of Imax (0.5% in SongAbbot2000, 0.6% in MehtaWilson2000)
Damp = params.Damp*Imax; % peak Post-before-Pre weight change, in percent of Imax (0.525% in SongAbbott2000, 90% in MehtaWilson2000, which didn't have a maximum weight)

% output PF 
Nbin = params.Nbin; % Number of bins in which the length of the track is divided

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
speed = L/period; % cm/sec
EndTime = Nlaps*L/speed; %in sec
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

%% input-output dynamics using Euler's forward method to solve ODEs
I = zeros(1,length(Trun)); 
V = I;
D = I;
OutRaster = I;
P = zeros(N,length(Trun));
Pact = P;
Dact = P;
V(1) = Vrest;

for i = 1:length(Trun)-1
    % input current to output neuron, EPSCs modelled as pulse with exponential decay
    I(i+1)  = I(i) - dt*I(i)./Idecay + sum(InRasters(:,i).*W(:,i));
    % I&F output neuron
    V(i+1) = V(i) + dt * (Vrest - V(i) + I(i) * Rm) / tau_m; % subthreshold membrane potential
    if V(i+1) >= Vthr
       V(i) = 0; %to visualize spikes with peak at 0 mV
       V(i+1) = Vreset;
       OutRaster(i) = 1;
    end
    % STDP variables P and D + update Synaptic weights W
    D(i+1) = D(i) - dt*D(i)/Ddecay + OutRaster(i); % output spiketrain convolved with Post-before-Pre update rule, in %
    for n = 1:N % for each presynaptic neuron
        P(n,i+1) = P(n,i) - dt*P(n,i)/Pdecay + InRasters(n,i); % input spiketrains convolved with Pre-before-Post update rule, in %
        Pact(n,i+1) = OutRaster(i)*P(n,i)*Pamp;
        Dact(n,i+1) = InRasters(n,i)*D(i)*Damp;
        W(n,i+1) = W(n,i) + Pact(n,i+1) + Dact(n,i+1); % weight update
        if W(n,i+1) < 0
            W(n,i+1) = 0;
        elseif W(n,i+1) > Imax
            W(n,i+1) = Imax;
        end           
    end    
end
Spiketimes_outAll = Trun(find(OutRaster));
SpikeLoc_outAll = Run(find(OutRaster)); %spikes location on track

% segment output raster by lap + compute spatial firing rate
NewLapIdx = find(Run==0);
NewLapTimes = Trun(NewLapIdx);
Track_dbin = L./Nbin; % bin size, in cm
TrackBins = Lap1(1):Track_dbin:Lap1(end); % Spatial bin edges
BinCenters = L/(2*Nbin):L/Nbin:L;
for lap = 1:Nlaps
    if lap == Nlaps
    OutRaster_laps{lap} = OutRaster(NewLapIdx(lap):end);
    LapN = Run(NewLapIdx(lap):end);
    else
    OutRaster_laps{lap} = OutRaster(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    LapN = Run(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    end
    Spiketimes_out{lap} = Tlap1(find(OutRaster_laps{lap}));
    SpikeLoc_out{lap} = Lap1(find(OutRaster_laps{lap}));
    % compute lap-wise place field
    SpikeCountOut(lap, :) = histcounts(SpikeLoc_out{lap}, TrackBins); % number of spikes in each spatial bin
    TimeInBins(lap,:) = dt*histcounts(LapN, TrackBins); % time spent (in sec) in each spatial bin (no assumption on trajectories)
    SpatialFRout(lap,:) = SpikeCountOut(lap, :)./TimeInBins(lap, :); % FR in each spatial bin, in Hz 
    %Compute lap-wise COM, SD, skewness, mean and max firing rate
    COMbin(lap) = sum(SpatialFRout(lap,:).*[1:Nbin])/sum(SpatialFRout(lap,:));
    COMloc(lap) = sum(SpatialFRout(lap,:).*BinCenters)/sum(SpatialFRout(lap,:));
    PFsd(lap) = sqrt( sum( (BinCenters - COMloc(lap)).^2.*SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ) ); %lap-wise SD of the PF
    PFskew(lap) = sum( ((BinCenters-COMloc(lap))./PFsd(lap) ).^3.* SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ); %lap-wise skewness of the PF 
    meanFRout_lap(lap) = sum(OutRaster_laps{lap})/(dt*length(LapN)); % average FR over time, in Hz
    maxFRout_lap(lap) = max(SpatialFRout(lap,:)); % spatial bin with max FR, in Hz
end
meanFRmap = mean(SpatialFRout,1); %average FR map (i.e. average PF)
COM_meanPF = sum(meanFRmap.*BinCenters)/sum(meanFRmap);
SD_meanPF = sqrt( sum( (BinCenters - COM_meanPF).^2.*meanFRmap/sum(meanFRmap) ) );
Skew_meanPF = sum( ((BinCenters-COM_meanPF)/SD_meanPF ).^3.* meanFRmap/sum(meanFRmap) );

[b,~,~,~,stats] = regress(COMloc', [ones(Nlaps,1), [1:Nlaps]']);
COMtraj = b(1) + b(2)*[1:Nlaps];

%% output
output.Spiketimes_out = Spiketimes_out;
output.SpikeLoc_out = SpikeLoc_out;
output.PF = SpatialFRout; % matrix of FR per bin, Nlaps rows and Nbin columns
output.COMloc = COMloc; %lap-wise COM
output.COMbin = COMbin;
output.PFsd = PFsd; %lap-wise SD (PF width metric)
output.PFskew = PFskew; % lap-wise PF skewness
output.meanFRout_lap = meanFRout_lap; %average FR for each lap
output.maxFRout_lap = maxFRout_lap; %bin with max FR, for each lap
output.COMslope = b(2); % slope of the lap-wise COM regression (cm/lap) 
output.shiftR2 = stats(1); % Rsquare of the lap-wise COM regression.
output.shiftPval = stats(3); %pvalue of the lap-wise COM regression.
output.meanFRmap = meanFRmap; % average PF over all laps
output.COM_meanPF = COM_meanPF; % COM of average PF
output.SD_meanPF = SD_meanPF; % SD of average PF
output.Skew_meanPF = Skew_meanPF; % skewness of average PF

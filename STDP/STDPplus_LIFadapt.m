function output = STDPplus_LIFadapt(params)
%
%params: struct with the following fields:
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
InShift = params.InShift; % 0 if static input, 1 if shift
PslopeBinCenters = params.PslopeBinCenters;
Pslope = params.Pslope;

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
Adaptation = params.Adapt; % 1 or 0, to implement spike-rate adaptation or not
Ek = params.Ek; % K+ equilibrium potential, in Volt
dSRA = params.dSRA; % increment value of spike-rate adaptation variable modeling adapting K+ leak current (cf DayanAbbott fig 5.6 and eq 5.13)
tau_sra = params.tau_sra; % time constant for exponential decay of SRA variable
Shunting = params.Shunting; % implement shunting inhib (capping I) if = 1, 0 otherwise. 
Icap = params.Icap;

% Plasticity params
Pdecay = params.Pdecay; % Pre-before-Post time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Ddecay = params.Ddecay; % Post-before-Pre time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Pamp = params.Pamp*Imax; % peak Pre-before-Post weight change, in percent of Imax (0.5% in SongAbbot2000, 0.6% in MehtaWilson2000)
Damp = params.Damp*Imax; % peak Post-before-Pre weight change, in percent of Imax (0.525% in SongAbbott2000, 90% in MehtaWilson2000, which didn't have a maximum weight)
capWeights = params.capWeights; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 
HomeoNorm = params.HomeoNorm; % 0 if no homeostatic normalization rule, 1 if yes.
WUdyn = params.WUdyn; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Wrest = params.Wrest;
tau_wu = params.tau_wu; % weight update time constant, in sec

% output PF 
Nbin = params.Nbin; % Number of bins in which the length of the track is divided

%% Inputs

% randomly select slope for each input neuron according to bin proba
InSlopes = randsample(PslopeBinCenters,N,true,Pslope); 

%Place fields (FR = gaussian function of position x)
x = 1:L; % track positions binned by cm
PFcom = linspace(1,L,N) %1:L/N:L; % PF centers spread uniformly along the track, with L/N spacing in cm

if InShift == 1 % dynamic inputs
   PFcomD = zeros(Nlaps,N);
   for l = 1:Nlaps
   PFcomD(l,:) = InSlopes.*(l-1) + PFcom;
   end
   PFcomD(PFcomD < 1)=1;
   PFcomD(PFcomD > L)=L;
end

for i = 1:N
    if InShift == 0
    FRinNorm(i,x) = exp(-0.5*(x-PFcom(i)).^2/PFsd^2); %normalized from 0 to 1
    FRin(i,x) = PFamp*exp(-0.5*(x-PFcom(i)).^2/PFsd^2); % real FR 
    elseif InShift == 1 % dynamic inputs
        for l = 1:Nlaps
            FRin(l,i,x) = PFamp*exp(-0.5*(x-PFcomD(l,i)).^2/PFsd^2); % real FR 
        end
    end
end

%Running model with constant speed
speed = L/period; % cm/sec
EndTime = Nlaps*L/speed; %in sec
Trun = 0:dt:EndTime; Trun = Trun(1:end-1); % ms resolution
Tlap1 = 0:dt:period; Tlap1 = Tlap1(1:end-1);
Lap1 = linspace(0,L,period/dt); % positions during 1 lap
Run = repmat(Lap1,[1 Nlaps]);%positions during the whole run

% Lap indices and times for given trajectory
NewLapIdx = find(Run==0);
% NewLapTimes = Trun(NewLapIdx);

%convert trajectory into firing rate
   for i = 1:length(Run)
   [~, idx] = min(abs(x-Run(i))); %find closest spatial bin in x for each position in trajectory
       if InShift == 0
        FRrun(:,i) = FRin(:,idx); % Firing rate at that position for each neuron
       elseif InShift == 1  
        [~, l] = min(abs(NewLapIdx-i)); % what lap is it?
        FRrun(:,i) = FRin(l,:,idx);
       end
   end

% figure
% imagesc(FRrun);
% xlabel('time (s)'); ylabel('input neuron')
% % xlim([500000 600000])
% colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
% c = colorbar; c.Label.String = 'Firing rate (Hz)';
% box off; axis square

% nonhomogenous Poisson spike generator along trajectory
InRasters = poissrnd(FRrun.*dt); %matrix of size (N,length(Trun)) containing, for each input neuron, the number of spikes in each 1ms bin drawn from Poisson distribution
InRasters(InRasters>=1) = 1; 
InRasters(InRasters < 1) = 0; 
InRastersTF = logical(InRasters);
Spiketimes_in = cell(N,1);
for i = 1:N
    Spiketimes_in{i} = Trun( InRastersTF(i,:) ); %find positive values = when a spike occured
end

% for i = 1:N
%     Spiketimes_in{i} = Trun( find(InRasters(i,:)) ); %find positive values = when a spike occured
% end

% gaussian connectivity with max in center of track, akin to Mehta and Wilson 2000
NeuronID = 1:N;
W = zeros(N,length(Trun));
Wmax = maxW*Imax;
W(:,1) = Wmax*exp(-0.5*(NeuronID - N/2).^2/Wsd^2); % vector of N rows following a gaussian distribution of Wmax amplitude 
Wtarget = W;

%% input-output dynamics using Euler's forward method to solve ODEs
I = zeros(1,length(Trun)); 
V = I;
SRA = I;
D = I;
OutRaster = I;
P = zeros(N,length(Trun));
Pact = P;
Dact = P;
V(1) = Vrest;

for i = 1:length(Trun)-1
    
    % input current to output neuron, EPSCs modelled as pulse with exponential decay
    I(i+1)  = I(i) - dt*I(i)./Idecay + sum(InRasters(:,i).*W(:,i));
    if Shunting == 1 
        if I(i+1) > Icap % implements shunting inhib rather than bounding weights
        I(i+1) = Icap;
        end
    end

    % output neuron
    if Adaptation == 0
        V(i+1) = V(i) + dt * (Vrest - V(i) + I(i) * Rm) / tau_m; % subthreshold membrane potential
    elseif Adaptation == 1
        V(i+1) = V(i) + dt .* (Vrest - V(i) -SRA(i).*(V(i)-Ek) + I(i) .* Rm) ./ tau_m;
        SRA(i+1) = SRA(i) - dt .* SRA(i)./tau_sra; % spike-rate-adaptation variable with exponential decay to 0
    end
    
    if V(i+1) >= Vthr
       SRA(i+1) = SRA(i+1) + dSRA;
       V(i) = 0; % spikes with peak at 0 mV
       V(i+1) = Vreset;
       OutRaster(i) = 1;
    end

    % STDP variables
    D(i+1) = D(i) - dt*D(i)/Ddecay + OutRaster(i); % output spiketrain convolved with Post-before-Pre update rule, in %
    for n = 1:N % for each presynaptic neuron
        P(n,i+1) = P(n,i) - dt*P(n,i)/Pdecay + InRasters(n,i); % input spiketrains convolved with Pre-before-Post update rule, in %
        Pact(n,i+1) = OutRaster(i)*P(n,i)*Pamp;
        Dact(n,i+1) = InRasters(n,i)*D(i)*Damp;
        
        % update Synaptic weights W, with or without dynamics
        if WUdyn == 0 % instantaneous update
            W(n,i+1) = W(n,i) + Pact(n,i+1) + Dact(n,i+1); % weight update
        elseif WUdyn == 1 % dynamic update with time constant tau_wu
            Wtarget(n,i+1) = Wtarget(n,i) + Pact(n,i+1) + Dact(n,i+1);
            W(n,i+1) = W(n,i) + dt .* (Wrest - W(n,i) + Wtarget(n,i+1))./ tau_wu;    
        end
        
        % Hard bounds on weights
        if W(n,i+1) < 0
            W(n,i+1) = 0;
        elseif W(n,i+1) > Imax && capWeights == 1 % as in SongAbbott2000, weights are capped to prevent runaway plasticity and FR
            W(n,i+1) = Imax;
        end
    end
    
    if HomeoNorm == 1 % implement heterosynaptic competition or not
       W(:,i+1) = W(:,i+1).*sum(W(:,1))./sum(W(:,i+1)); % Homeostatic synaptic normalization, multiplicative method
    end
end
% Spiketimes_outAll = Trun(find(OutRaster));
% SpikeLoc_outAll = Run(find(OutRaster)); %spikes location on track

%% compute spatial properties by lap
Track_dbin = L./Nbin; % bin size, in cm
TrackBins = Lap1(1):Track_dbin:Lap1(end); % Spatial bin edges
BinCenters = L/(2*Nbin):L/Nbin:L;

Track_dbinMou = 2.5; % 2.5cm bins like in MouJi2018, for comparison
TrackBinsMou = Lap1(1):Track_dbinMou:Lap1(end);

OutRaster_laps = cell(Nlaps,1);
Spiketimes_out = cell(Nlaps,1);
SpikeLoc_out = cell(Nlaps,1);
SpikeCountOut = zeros(Nlaps,Nbin);
TimeInBins = zeros(Nlaps,Nbin);
% TimeInBinsMou = zeros(Nlaps,120);
SpatialFRout = zeros(Nlaps,Nbin);
% SpatialFRoutMou = zeros(Nlaps,120);
COMbin = zeros(1, Nlaps);
COMloc = zeros(1,Nlaps);
PFsdOut = zeros(1,Nlaps);
PFskew = zeros(1,Nlaps);
meanFRout_lap = zeros(1,Nlaps);
maxFRout_lap = zeros(1,Nlaps);

for lap = 1:Nlaps
    if lap == Nlaps
    OutRaster_laps{lap} = OutRaster(NewLapIdx(lap):end);
    LapN = Run(NewLapIdx(lap):end);
    else
    OutRaster_laps{lap} = OutRaster(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    LapN = Run(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    end
    Spiketimes_out{lap} = Tlap1(logical(OutRaster_laps{lap}));
    SpikeLoc_out{lap} = Lap1(logical(OutRaster_laps{lap}));
    % compute lap-wise place field
    SpikeCountOut(lap, :) = histcounts(SpikeLoc_out{lap}, TrackBins); % number of spikes in each spatial bin
    SpikeCountOutMou(lap, :) = histcounts(SpikeLoc_out{lap}, TrackBinsMou);
    TimeInBins(lap,:) = dt*histcounts(LapN, TrackBins); % time spent (in sec) in each spatial bin (no assumption on trajectories)
    TimeInBinsMou(lap,:) = dt*histcounts(LapN, TrackBinsMou);
    SpatialFRout(lap,:) = SpikeCountOut(lap, :)./TimeInBins(lap, :); % FR in each spatial bin, in Hz 
    SpatialFRoutMou(lap,:) = SpikeCountOutMou(lap, :)./TimeInBinsMou(lap, :);
    %Compute lap-wise COM, SD, skewness, mean and max firing rate
    COMbin(lap) = sum(SpatialFRout(lap,:).*[1:Nbin])/sum(SpatialFRout(lap,:));
    COMloc(lap) = sum(SpatialFRout(lap,:).*BinCenters)/sum(SpatialFRout(lap,:));
    PFsdOut(lap) = sqrt( sum( (BinCenters - COMloc(lap)).^2.*SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ) ); %lap-wise SD of the PF
    PFskew(lap) = sum( ((BinCenters-COMloc(lap))./PFsdOut(lap) ).^3.* SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ); %lap-wise skewness of the PF 
    meanFRout_lap(lap) = sum(OutRaster_laps{lap})/(dt*length(LapN)); % average FR over time, in Hz
    maxFRout_lap(lap) = max(SpatialFRout(lap,:)); % spatial bin with max FR, in Hz
end
meanFRout = sum(OutRaster)./EndTime; %in Hz
meanFRmap = mean(SpatialFRout,1); %average FR map (i.e. average PF)
meanFRmapMou = mean(SpatialFRoutMou,1);
COM_meanPF = sum(meanFRmap.*BinCenters)/sum(meanFRmap);
SD_meanPF = sqrt( sum( (BinCenters - COM_meanPF).^2.*meanFRmap/sum(meanFRmap) ) );
Skew_meanPF = sum( ((BinCenters-COM_meanPF)/SD_meanPF ).^3.* meanFRmap/sum(meanFRmap) );

[b,~,~,~,stats] = regress(COMloc', [ones(Nlaps,1), [1:Nlaps]']);
COMtraj = b(1) + b(2)*[1:Nlaps];

%% output

output.PF = SpatialFRout; % matrix of FR per bin, Nlaps rows and Nbin columns
output.COMloc = COMloc; %lap-wise COM
output.COMbin = COMbin;
output.PFsdOut = PFsdOut; %lap-wise SD (PF width metric)
output.PFskew = PFskew; % lap-wise PF skewness
output.meanFRout = meanFRout; %average FR of the output neuron
output.meanFRout_lap = meanFRout_lap; %average FR for each lap
output.maxFRout_lap = maxFRout_lap; %bin with max FR, for each lap
output.COMslope = b(2); % slope of the lap-wise COM regression (cm/lap) 
output.COMtraj = COMtraj; % PF COM trajectory inferred from regression
output.shiftR2 = stats(1); % Rsquare of the lap-wise COM regression.
output.shiftPval = stats(3); %pvalue of the lap-wise COM regression.
output.meanFRmap = meanFRmap; % average PF over all laps
output.meanFRmapMou = meanFRmapMou; % average PF over all laps computed with 2.5cm bins like in Mou and Ji 2018
output.COM_meanPF = COM_meanPF; % COM of average PF
output.SD_meanPF = SD_meanPF; % SD of average PF
output.Skew_meanPF = Skew_meanPF; % skewness of average PF

% % Warning: loading all that in the workspace demands a lot of RAM.
% output.I = I;
% output.V = V;
% output.D = D;
% output.P = P;
% output.W = W;
% output.Wtarget = Wtarget;
% output.Spiketimes_in = Spiketimes_in;
% output.Spiketimes_out = Spiketimes_out;
% output.SpikeLoc_out = SpikeLoc_out;

clearvars -EXCEPT output
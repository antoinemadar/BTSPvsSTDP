function output = BTSPplus_LIFadapt(params)
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
Adaptation = params.Adapt; % 1 or 0, to implement spike-rate adaptation or not
Ek = params.Ek; % K+ equilibrium potential, in Volt
dSRA = params.dSRA; % increment value of spike-rate adaptation variable modeling adapting K+ leak current (cf DayanAbbott fig 5.6 and eq 5.13)
tau_sra = params.tau_sra; % time constant for exponential decay of SRA variable
% Shunting = params.Shunting; % implement shunting inhib (capping I) if = 1, 0 otherwise. 
% Icap = params.Icap;

% Plasticity params
CSproba = params.CSproba; % proba of an output spike to be a plateau potential/complex spike (1.8% on average in Bittner et al. 2017)
pCSdynamic = params.pCSdynamic;
Apcs = params.Apcs; % max frequency of CSs. Cf Supp Fig associated to Fig 7. Same max amplitude for Familiar (F) and Novel (N) conditions.
Tau_pcs = params.Tau_pcs; % in laps (cf Fig 6 and 7 of our manuscript.). Same time constants for F and N (but exp fits on instantaneous MSD suggest 1.15 for N and 0.85 for F)
Bpcs = params.Bpcs; % N>F>0. 0 results in no shift induced in later laps (cf Fig 7: a static pCS is not enough to lead to high Diffusion coeff in later laps)
SDcs = params.SDcs; % in cm (same as Params.L); Standard deviation of the gaussian defining where a CS can occur around the current COM of the weights
Bound = params.Bound;
LowBound = round(L/2)-Bound;
HighBound = round(L/2)+Bound;

Pdecay = params.Pdecay; % Pre-before-Post time constant, in sec (1.31s on average in Bittner et al. 2017)
Ddecay = params.Ddecay; % Post-before-Pre time constant, in sec (0.69s on average in Bittner et al. 2017)
Pamp = params.Pamp; % peak weight change, in Amps (3*baseline EPSP for pairings with 10 EPSPs at 10Hz, in Bittner et al. 2017)
PostPreBoost = params.PostPreBoost; % amplitude of weight change for postpre rule which is triggered on CS, not on spikes (so, to be somewhat balanced, it should = 10 so that amplitude is the Bittner2017 amplitude/1CS and not /10spikes like the prepost rule)

capWeights = params.capWeights; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 

HomeoNorm = params.HomeoNorm; % 0 if no homeostatic normalization rule, 1 if yes.

WUdyn = params.WUdyn; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Wrest = params.Wrest;
tau_wu = params.tau_wu; % weight update time constant, in sec

% output PF 
Nbin = params.Nbin; % Number of bins in which the length of the track is divided

%% Inputs
%Place fields (FR = gaussian function of position x)
x = 1:L; % track positions binned by cm
PFcom = 1:L/N:L; % PF centers spread uniformly along the track, with L/N spacing in cm
FRin = zeros(N,length(x));
for i = 1:N
%     FRinNorm(i,x) = exp(-0.5*(x-PFcom(i)).^2/PFsd^2); %normalized from 0 to 1
    FRin(i,x) = PFamp*exp(-0.5*(x-PFcom(i)).^2/PFsd^2); % real FR in Hz
end

%Running model with constant speed
speed = L/period; % cm/sec
EndTime = Nlaps*L/speed; %in sec
Trun = 0:dt:EndTime; Trun = Trun(1:end-1); % ms resolution
Tlap1 = 0:dt:period; Tlap1 = Tlap1(1:end-1);
Lap1 = linspace(0,L,period/dt); % positions during 1 lap
Run = repmat(Lap1,[1 Nlaps]);%positions during the whole run

% segment by lap
NewLapIdx = find(Run==0);
    % NewLapTimes = Trun(NewLapIdx);
laps_start = zeros(size(Run));
% laps_start(Run==0) = 1; % 1 at beginning of a new lap
for nl = 1:Nlaps
    laps_start(NewLapIdx(nl)) = nl; % tracks the lap number
end
laps_end = zeros(size(Run));
laps_end(Run==L) = 1; % 1 on end of laps
lapL = find(laps_end,1); % number of indices in a lap

%convert trajectory into firing rate
FRrun = zeros(N,length(Run));
for i = 1:length(Run)
    [~, idx] = min(abs(x-Run(i))); %find closest spatial bin in x for each position in trajectory
    FRrun(:,i) = FRin(:,idx); % Firing rate at that position for each neuron
end

% nonhomogenous Poisson spike generator along trajectory
InRasters = poissrnd(FRrun.*dt); %matrix of size (N,length(Tlap1)) containing, for each input neuron, the number of spikes in each 1ms bin drawn from Poisson distribution
InRasters(InRasters>=1) = 1;  
InRasters(InRasters < 1) = 0; 
InRastersTF = logical(InRasters);
Spiketimes_in = cell(N,1);
for i = 1:N
    Spiketimes_in{i} = Trun( InRastersTF(i,:) ); %find positive values = when a spike occured
end

% gaussian connectivity with max in center of track, akin to Mehta and Wilson 2000
NeuronID = 1:N;
W = zeros(N,length(Trun));
Wmax = maxW*Imax;
W(:,1) = Wmax*exp(-0.5*(NeuronID - N/2).^2/Wsd^2); % vector of N rows following a gaussian distribution of Wmax amplitude 

%% for dynamic p(CS)
if pCSdynamic == 1
laps = 0:28;
pCSdyn = Apcs.*exp(-laps./Tau_pcs) + Bpcs; % lapwise p(CS) for lap 1 to 29 
CSok = binornd(1,pCSdyn); %determine if there is a CS in each lap, given lapwise proba
end

%% input-output dynamics using Euler's forward method to solve ODEs
I = zeros(1,length(Trun)); 
V = I;
D = I;
SRA = I;
OutRaster = I;
CS = I;
P = zeros(N,length(Trun));
Pact = P;
Dact = P;
V(1) = Vrest;
Wtarget = W;
W2 = W;

for i = 1:length(Trun)-1
    % input current to output neuron, EPSCs modelled as pulse with exponential decay
    if HomeoNorm == 1
    I(i+1)  = I(i) - dt*I(i)./Idecay + sum(InRasters(:,i).*W2(:,i));
    else
    I(i+1)  = I(i) - dt*I(i)./Idecay + sum(InRasters(:,i).*W(:,i)); 
    end
    
    % I&F output neuron
    
    if Adaptation == 0
        V(i+1) = V(i) + dt * (Vrest - V(i) + I(i) * Rm) / tau_m; % subthreshold membrane potential
    elseif Adaptation == 1
        V(i+1) = V(i) + dt .* (Vrest - V(i) -SRA(i).*(V(i)-Ek) + I(i) .* Rm) ./ tau_m;
        SRA(i+1) = SRA(i) - dt .* SRA(i)./tau_sra; % spike-rate-adaptation variable with exponential decay to 0
    end

    if V(i+1) >= Vthr
       V(i) = 0; %to visualize spikes with peak at 0 mV
       V(i+1) = Vreset;
       SRA(i+1) = SRA(i+1) + dSRA;
       OutRaster(i) = 1;
       if pCSdynamic == 0
          CS(i) = binornd(1, CSproba);
       end
    end
    
    if pCSdynamic==1 && laps_start(i)>0 && laps_start(i)<Nlaps && CSok(laps_start(i))>0 % at the beginning of a new lap where a CS will occur, select CS location and index
        if HomeoNorm == 1
            Wcom = sum(W2(:,i)'.*PFcom)/sum(W2(:,i)); % find COM of the current weights
        else
            Wcom = sum(W(:,i)'.*PFcom)/sum(W(:,i)); % find COM of the current weights
        end
        CSx = SDcs.*randn(1) + Wcom; % position of the CS on the track is randomly sampled from a normal distribution centered on current weight COM and Params.SDcs standard deviation
        CSx = max(LowBound,min(CSx,HighBound)); % restrict CSx between Low and High track boundaries (e.g. 75 and 225cm, for a +/-75cm around middle of the 300cm track)
        [~, CSi] = min(abs(Lap1-CSx));% find the corresponding index of CSx
        CS(CSi+i-1) = 1;
        clear CSx CSi Wcom
    end

    % BTSP variables P and D + update Synaptic weights W
    D(i+1) = D(i) - dt*D(i)/Ddecay + CS(i); % output train of Complex Spikes convolved with Post-before-Pre update rule, in %
    for n = 1:N % for each presynaptic neuron
        P(n,i+1) = P(n,i) - dt*P(n,i)/Pdecay + InRasters(n,i); % input spiketrains convolved with Pre-before-Post update rule, in %
        Pact(n,i+1) = CS(i)*P(n,i)*Pamp; %evaluate pre-before-post variable at time of CS
        Dact(n,i+1) = InRasters(n,i)*D(i)*Pamp*PostPreBoost; %evaluate post-before-pre variable at time of an input spike
        
        % update Synaptic weights W, with or without dynamics
        if WUdyn == 1 %&& HomeoNorm == 0 % dynamic update with time constant tau_wu
            Wtarget(n,i+1) = Wtarget(n,i) + Pact(n,i+1) + Dact(n,i+1);
            W(n,i+1) = W(n,i) + dt .* (Wrest - W(n,i) + Wtarget(n,i+1))./ tau_wu;
        else
            W(n,i+1) = W(n,i) + Pact(n,i+1) + Dact(n,i+1);
        end

        % Hard bounds on weights
        if W(n,i+1) < 0
            W(n,i+1) = 0;
        elseif W(n,i+1) > Imax && capWeights == 1 % as in SongAbbott2000, weights are capped to prevent runaway plasticity and FR
            W(n,i+1) = Imax;
        end          
    end
    
    if HomeoNorm == 1 % implement heterosynaptic competition or not  
        % if WUdyn == 0
            W2(:,i+1) = W(:,i+1).*sum(W(:,1))./sum(W(:,i+1)); % Homeostatic synaptic normalization, multiplicative method
        % else % WUdyn == 1
        %     Wtarget(:,i+1) = W(:,i+1).*sum(W(:,1))./sum(W(:,i+1));
        % %     W2(n,i+1) = W(n,i) + dt .* (Wrest - W(n,i) + Wtarget(n,i+1))./ tau_wu;      % this cannot work here because it has to be in the for loop through the n inputs          
        %     end
        % end
    end
end
% Spiketimes_outAll = Trun(logical(OutRaster));
% SpikeLoc_outAll = Run(logical(OutRaster)); %spikes location on track
% CSloc_All = Run(logical(CS)); % CSs location on track

%% segment output raster by lap + compute spatial firing rate

% bin the track
Track_dbin = L./Nbin; % bin size, in cm
TrackBins = Lap1(1):Track_dbin:Lap1(end); % Spatial bin edges
BinCenters = L/(2*Nbin):L/Nbin:L;

%pre-allocate variables
OutRaster_laps = cell(Nlaps,1);
CS_laps = cell(Nlaps,1);
CStimes = cell(Nlaps,1);
CSloc = cell(Nlaps,1);
Spiketimes_out = cell(Nlaps,1);
SpikeLoc_out = cell(Nlaps,1);
SpikeCountOut = zeros(Nlaps,Nbin);
CSbin = zeros(Nlaps,Nbin);
TimeInBins = zeros(Nlaps,Nbin);
SpatialFRout = zeros(Nlaps,Nbin);
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
    CS_laps{lap} = CS(NewLapIdx(lap):end);
    else
    OutRaster_laps{lap} = OutRaster(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    CS_laps{lap} = CS(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    LapN = Run(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    end
    Spiketimes_out{lap} = Tlap1(logical(OutRaster_laps{lap}));
    SpikeLoc_out{lap} = Lap1(logical(OutRaster_laps{lap}));
    CStimes{lap} = Tlap1(logical(CS_laps{lap}));
    CSloc{lap} = Lap1(logical(CS_laps{lap}));
    
    % compute lap-wise place field
    SpikeCountOut(lap, :) = histcounts(SpikeLoc_out{lap}, TrackBins); % number of spikes in each spatial bin
    CSbin(lap,:) = histcounts(CSloc{lap}, TrackBins); 
    TimeInBins(lap,:) = dt*histcounts(LapN, TrackBins); % time spent (in sec) in each spatial bin (no assumption on trajectories)
    SpatialFRout(lap,:) = SpikeCountOut(lap, :)./TimeInBins(lap, :); % FR in each spatial bin, in Hz 
    
    %Compute lap-wise COM, SD, skewness, mean and max firing rate
    COMbin(lap) = sum(SpatialFRout(lap,:).*[1:Nbin])/sum(SpatialFRout(lap,:));
    COMloc(lap) = sum(SpatialFRout(lap,:).*BinCenters)/sum(SpatialFRout(lap,:));
    PFsdOut(lap) = sqrt( sum( (BinCenters - COMloc(lap)).^2.*SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ) ); %lap-wise SD of the PF
    PFskew(lap) = sum( ((BinCenters-COMloc(lap))./PFsdOut(lap) ).^3.* SpatialFRout(lap,:)/sum(SpatialFRout(lap,:)) ); %lap-wise skewness of the PF 
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
output.CStimes = CStimes;
output.CSloc = CSloc;
output.CSbin = CSbin;
output.PF = SpatialFRout; % matrix of FR per bin, Nlaps rows and Nbin columns
output.PFsdOut = PFsdOut; %lap-wise SD (PF width metric)
output.PFskew = PFskew; % lap-wise PF skewness
output.meanFRout_lap = meanFRout_lap; %average FR for each lap
output.maxFRout_lap = maxFRout_lap; %lap-wise PF peak FR (Hz) (defined as bin with max FR)
output.COMloc = COMloc; %lap-wise COM
output.COMbin = COMbin;
output.COMtraj = COMtraj; % PF COM trajectory inferred from regression
output.COMslope = b(2); % PF shift: slope of the lap-wise COM regression (cm/lap) 
output.shiftR2 = stats(1); % Rsquare of the lap-wise COM regression.
output.shiftPval = stats(3); %pvalue of the lap-wise COM regression.
output.meanFRmap = meanFRmap; % average PF over all laps
output.COM_meanPF = COM_meanPF; % COM of average PF
output.SD_meanPF = SD_meanPF; % SD of average PF
output.Skew_meanPF = Skew_meanPF; % skewness of average PF

% % Warning: loading all that in the workspace demands a lot of RAM when Nsim >> 1
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
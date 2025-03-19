function output = BTSP_MilsteinModel_LIFadapt(params)
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

% BTSP-induction params
CSproba = params.CSproba; % proba of an output spike to be a plateau potential/complex spike (1.8% on average in Bittner et al. 2017)

%Milstein biBTSTP rule params
ETdecay = params.Pdecay; % 
CSdecay = params.Ddecay; %
k_pot = params.k_pot;
k_dep = params.k_dep;
a_pot = params.a_pot;
b_pot = params.b_pot;
a_dep = params.a_dep;
b_dep = params.b_dep;
Wmax = params.Wmax;

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

laps_end = zeros(size(Run));
laps_end(Run==L) = 1; % 1 on end of laps

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
W(:,1) = maxW*Imax*exp(-0.5*(NeuronID - N/2).^2/Wsd^2); % vector of N rows following a gaussian distribution of Wmax amplitude 

% segment by lap
NewLapIdx = find(Run==0);
% NewLapTimes = Trun(NewLapIdx);

%% input-output dynamics using Euler's forward method to solve ODEs
lapL = find(laps_end,1); % number of indices in a lap

%initial values
I = zeros(1,size(W,2)); 
CS = zeros(1,length(Trun));
plateau = CS;
V = I;
D = I;
SRA = I;
OutRaster = I;
P = zeros(size(W));
V(1) = Vrest;
V2 = V;
W2 = W;

for i = 1:size(W,2)-1 % time points

    % input current to output neuron, EPSCs modelled as pulse with exponential decay
    I(i+1)  = I(i) - dt*I(i)./Idecay + sum(InRasters(:,i).*W(:,i)); 
    
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
       CS(i) = binornd(1, CSproba);
           if CS(i) == 1
           plateau(i:i+299) = 1;
           end
    end
    V2(i+1) = V(i+1); % to keep a version of V without spikes

    % BTSP variables P and D + update Synaptic weights W
    D(i+1) = D(i) - dt*D(i)/CSdecay + plateau(i); % Global Instructive Signal from the CS

    for n = 1:size(W,1) % for each presynaptic neuron
        P(n,i+1) = P(n,i) - dt*P(n,i)/ETdecay + InRasters(n,i); % Synaptic Eligibility Trace
        
        % update Synaptic weights W at the end of every lap during which a CS occured
        if laps_end(i)==1 && sum(CS(i-lapL+1:i))>0 % at the end of laps, when the lap contained at least 1 CS
            peakIS = max(D(i-lapL+1:i));
            normIS = D(i-lapL+1:i)./peakIS;% normalize IS trace of the whole lap
            normET = P(n,i-lapL+1:i)./max(P(n,i-lapL+1:i));% normalize ET traces of the whole lap (normET is a matrix, with rows being synapses and columns time points)
            pot_rate = trapz(dt, q_sigmo(normET.*normIS,a_pot,b_pot));
            dep_rate = trapz(dt, q_sigmo(normET.*normIS,a_dep,b_dep));
            update = W(n,i) + (Wmax - W(n,i))*k_pot*pot_rate - W(n,i)*k_dep*dep_rate; % I think we don't divide the update term by the duration of a lap (by dt if we had updated at all timepoints) because trapz above already integrates over dt
            W(n,i+1) = update;
            W(n,i+1) = max(0, min(Wmax,update)); % hard bounds, the way implemented in Milstein et al 2021 (although the formulation above should already have soft bounds, if 0 < k constants * plasticity rates < 1...)
            clear normIS normET pot_rate dep_rate update
        elseif laps_end(i)==1 && sum(CS(i-lapL+1:i))==0 && max(D(i-lapL+1:i))>0 % for laps with no CS but residual IS trace from the decay
            normIS = D(i-lapL+1:i)./peakIS;% normalize IS trace to the peak of the last lap with a CS
            normET = P(n,i-lapL+1:i)./max(P(n,i-lapL+1:i));% normalize ET traces of the whole lap (normET is a matrix, with rows being synapses and columns time points)
            pot_rate = trapz(dt, q_sigmo(normET.*normIS,a_pot,b_pot));
            dep_rate = trapz(dt, q_sigmo(normET.*normIS,a_dep,b_dep));
            update = W(n,i) + (Wmax - W(n,i))*k_pot*pot_rate - W(n,i)*k_dep*dep_rate; % I think we don't divide the update term by the duration of a lap (by dt if we had updated at all timepoints) because trapz above already integrates over dt
            W(n,i+1) = update;
            W(n,i+1) = max(0, min(Wmax,update));
        else
            W(n,i+1) = W(n,i);
        end         
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
% output.plateau = plateau;
% output.Spiketimes_in = Spiketimes_in;
% output.Spiketimes_out = Spiketimes_out;
% output.SpikeLoc_out = SpikeLoc_out;

clearvars -EXCEPT output

%% local function
    function q = q_sigmo(x,a,b) % sigmoid defined on [0,1], being 0 at x=0 and 1 at x=1.
s = 1./(1+exp(-b.*(x-a)));
s0 = 1./(1+exp(-b.*(0-a)));
s1 = 1./(1+exp(-b.*(1-a)));
q = (s-s0)./(s1-s0);
    end
end
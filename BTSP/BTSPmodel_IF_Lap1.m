clear; close all

%% Parameters
%Params for inputs
N = 50; % number of input neurons
L = 300; % length of track in cm
PFsd = 18; %cm (from Can's CA3 data)
PFamp = 10; %peak FR in Hz
Nlaps = 10; % number of laps
period = 20; % lap duration, in sec
dt = 0.001; % time resolution in sec

% Synapses params
Imax = 145e-12; % max synaptic weight, a current in Amps. % NolanWiles2011 uses 10-226pA. Mehta uses 1.5-4nA. SongAbbot use dimensionless conductances with G between 0 and Gmax, but not sure about the initial weight matrix. 
Idecay = 10e-3; % EPSC time constant, in sec
Wsd = 1.6*PFsd/(L/N); %standard deviation of initial synaptic weight vector (for gaussian connectivity) 
maxW = 1; % from 0 to 1, proportion of Imax defining the initial maximum weight

%Integrate and Fire params for output neuron (values from SongAbbot2000 if not specified otherwise)
Rm = 100e6; % membrane resistance, in Ohm
tau_m = 20e-3 ; %membrane time constant, in sec
Vrest = -70e-3 ; % resting membrane potential, in Volt
Vthr = -54e-3; % spike threshold, in Volt
Vreset = -60e-3; % reset potential after spike, in Volt

% BTSP params
CSproba = 1/100; % proba of an output spike to be a plateau potential/complex spike (1.8% on average in Bittner et al. 2017)
Pdecay_PrePost = 1.31; % Pre-before-Post time constant, in sec (1.31s on average in Bittner et al. 2017)
Pdecay_PostPre = 0.69; % Post-before-Pre time constant, in sec (0.69s on average in Bittner et al. 2017)

% output PF 
Nbin = 50; % Number of bins in which the length of the track is divided

%% Inputs
%Place fields (FR = gaussian function of position x)
x = 1:L; % track positions binned by cm
PFcom = 1:L/N:L; % PF centers spread uniformly along the track, with L/N spacing in cm
FRin = zeros(N,length(x));
for i = 1:N
    FRinNorm(i,x) = exp(-0.5*(x-PFcom(i)).^2/PFsd^2); %normalized from 0 to 1
    FRin(i,x) = PFamp*exp(-0.5*(x-PFcom(i)).^2/PFsd^2); % real FR 
end

figure
imagesc(FRin); 
xlabel('position (cm)'); ylabel('input neurons');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Firing rate (Hz)';
box off; axis square 

%Running model with constant speed
speed = L/period; % cm/sec
EndTime = Nlaps*L/speed; %in sec
Trun = 0:dt:EndTime; Trun = Trun(1:end-1); % ms resolution
Tlap1 = 0:dt:period; Tlap1 = Tlap1(1:end-1);
Lap1 = linspace(0,L,period/dt); % positions during 1 lap
Run = repmat(Lap1,[1 Nlaps]);%positions during the whole run

figure
plot(Trun(1:end), Run, 'k');
xlabel('time (s)'); ylabel('position (cm)')
box off; axis square

%convert trajectory into firing rate
for i = 1:length(Lap1)
    [minval, idx] = min(abs(x-Lap1(i))); %find closest spatial bin in x for each position in trajectory
    FRlap1(:,i) = FRin(:,idx); % Firing rate at that position for each neuron
end

figure
plot(Tlap1, FRlap1, 'k');
xlabel('time (s)'); ylabel('Firing Rate (Hz)')
box off; axis square

% nonhomogenous Poisson spike generator along trajectory
InRasters = poissrnd(FRlap1.*dt); %matrix of size (N,length(Tlap1)) containing, for each input neuron, the number of spikes in each 1ms bin drawn from Poisson distribution
for i = 1:N
    Spiketimes_in{i} = Tlap1( find(InRasters(i,:)) ); %find positive values = when a spike occured
end

figure
[spikex, spikey] = makedisplayrasters(Spiketimes_in, 0);
displaydisplayrasters(spikex, spikey); hold on
grid off
axis square
xlabel('time, (s)')
ylabel('input neurons')

% convolve input spiketrains with exponential decay to get Post-before-Pre
% weight rule (Depression in classic STDP curves)
% Dshape = exp(-[0::endpoint]./Decay); % 
% D = conv2(InRasters, Dshape(:)); 

% Synaptic weights at time 0 between all input neurons and 1 output neuron (vector of length N) 
% MehtaWilson2000 started with a symmetric connectivity matrix. I could do something
% similar, assuming a CA1 PC already has a PF with COM in the center of
% track and strength of connection following a gaussian. 
% or 
% I could randomly assign weigths to start with. from a distribution that
% favors clusters, e.g. Poisson or exponential, or just from uniform
% distrib, based on BittnerMagee2015 finding that spatial inputs are
% about equal.
% or 
% For BTSP model of emergence: equal weights for all synapses. PF could appear with BTSP or
% inhomogeneities in input neurons spiking. 

% gaussian connectivity with max in center of track, aking to Mehta and Wilson 2000
NeuronID = 1:N;
W = zeros(N,length(Tlap1));
Wmax = maxW*Imax;
W(:,1) = Wmax*exp(-0.5*(NeuronID - N/2).^2/Wsd^2); % vector of N rows following a gaussian distribution of Wmax amplitude 

figure
stem(1:N, W(:,1)*10^12, 'k');
xlabel('Input Neurons'); ylabel('Synaptic Weights (pA) ')
box off; axis square

%% input-output dynamics using Euler's forward method to solve ODEs
I = zeros(1,length(Tlap1)); 
V = I;
OutRaster = I;
CS = I;
D = I;
P = zeros(N,length(Tlap1));
V(1) = Vrest;

for i = 1:length(Tlap1)-1
    % input current to output neuron, EPSCs modelled as pulse with exponential decay
    I(i+1)  = I(i) - dt*I(i)./Idecay + sum(InRasters(:,i).*W(:,i));
    % I&F output neuron
    V(i+1) = V(i) + dt * (Vrest - V(i) + I(i) * Rm) / tau_m; % subthreshold membrane potential
    if V(i+1) >= Vthr
       V(i) = 0; %to visualize spikes with peak at 0 mV
       V(i+1) = Vreset;
       OutRaster(i) = 1;
       CS(i) = binornd(1, CSproba); %Bernoulli trial (0 or 1) with p(success) = CSproba, to model plateau ptential (aka Complex Spike)occurrence
    end
    % STDP variables + update Synaptic weights
    D(i+1) = D(i) - dt*D(i)/Pdecay_PostPre + CS(i); % output spiketrain convolved with Post-before-Pre update rule, in %
    for n = 1:N % for each presynaptic neuron
        P(n,i+1) = P(n,i) - dt*P(n,i)/Pdecay_PrePost + InRasters(n,i); % input spiketrains convolved with Pre-before-Post update rule, in %
        W(n,i+1) = W(n,i) + CS(i)*P(n,i)*Pamp + InRasters(n,i)*D(i)*Pamp; % weight update
        if W(n,i+1) < 0
            W(n,i+1) = 0;
        elseif W(n,i+1) > Imax
            W(n,i+1) = Imax;
        end           
    end    
end
Spiketimes_out{1} = Tlap1(find(OutRaster));
SpikeLoc_out = Lap1(find(OutRaster)); %spikes location on track
CSloc_All = Run(find(CS)); % CSs location on track

% compute place field
Track_dbin = L./Nbin; % bin size, in cm
TrackBins = Lap1(1):Track_dbin:Lap1(end);
SpikeCountOut = histcounts(SpikeLoc_out, TrackBins);
TimeInBins = dt*histcounts(Lap1, TrackBins); % compute time spent (in sec) in each spatial bin (not assuming a linear running model but real trajectories)
SpatialFRout = SpikeCountOut./TimeInBins; % in Hz   

figure 
subplot(2,2,1)
plot(Tlap1, I.*10^12)
xlabel('time (s)'); ylabel('input current (pA)')
box off;
subplot(2,2,2)
plot(Tlap1, V.*10^3)
xlabel('time (s)'); ylabel('membrane potential (mV)')
box off;
subplot(2,2,3)
[spikex, spikey] = makedisplayrasters(Spiketimes_out, 0);
displaydisplayrasters(spikex, spikey); hold on
grid off
xlim([0,Tlap1(end)]); xlabel('time, (s)')
ylabel('lap')
title('output neuron')
subplot(2,2,4)
imagesc(SpatialFRout); 
xlabel('spatial bins'); ylabel('lap');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Firing rate (Hz)';
box off;

figure %STDP variables
subplot(2,1,1)
plot(Tlap1, D, 'r'); hold on
plot(Tlap1, OutRaster, 'k');
xlabel('time (s)'); 
title('output raster + Post-before-Pre weight update')
box off;
subplot(2,1,2)
plot(Tlap1, P(1,:), 'r'); hold on
plot(Tlap1, InRasters(1,:), 'k');
xlabel('time (s)'); 
title('input raster #1 + Pre-before-Post weight update')
box off;

figure % Weight changes over time
subplot(2,1,1)
imagesc(W.*10^12);
xlabel('time'); ylabel('input neuron');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Synaptic weight (pA)';
box off;
subplot(2,1,2)
plot(Tlap1, W(1,:).*10^12);
xlabel('time (s)'); ylabel('input current (pA)')
box off;







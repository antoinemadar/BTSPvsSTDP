clear; close all

%% Parameters
%Params for inputs
N = 100; % number of input neurons
L = 300; % length of track in cm
PFsd = 18; %cm (10 mode 16-18 mean from Can's CA3 data)
PFamp = 20; %peak FR in Hz
Nlaps = 1; % number of laps
period = 20; % lap duration, in sec
dt = 0.001; % time resolution in sec

% Synapses params
Imax = 85e-12; % max synaptic weight, a current in Amps. % NolanWiles2011 uses 10-226pA. Mehta uses 1.5-4nA. SongAbbot use dimensionless conductances with G between 0 and Gmax, but not sure about the initial weight matrix. 
Idecay = 10e-3; % EPSC time constant, in sec
Wsd = 10; %1.6*PFsd/(L/N); %standard deviation of initial synaptic weight vector (for gaussian connectivity) 
maxW = 1; % from 0 to 1, proportion of Imax defining the initial maximum weight

%Integrate and Fire params for output neuron (values from SongAbbot2000 if not specified otherwise)
Rm = 100e6; % membrane resistance, in Ohm
tau_m = 20e-3 ; %membrane time constant, in sec
Vrest = -70e-3 ; % resting membrane potential, in Volt
Vthr = -54e-3; % spike threshold, in Volt
Vreset = -60e-3; % reset potential after spike, in Volt
% K = 0.7; % divisive inhibition constant as in Yu Shouva2006 (value unknown)
Adaptation = 0; % 1 or 0, to implement spike-rate adaptation or not
Ek = -70e-3; % K+ equilibrium potential, in Volt
dSRA = 0.06; % increment value of spike-rate adaptation variable modeling adapting K+ leak current (cf DayanAbbott fig 5.6 and eq 5.13)
tau_sra = 100e-3; % time constant for exponential decay of SRA variable

% STDP params
Pdecay = 20e-3; % Pre-before-Post time constant, in sec (20ms in SongAbbott2000, 10ms in MehtaWilson2000)
Ddecay = 20e-3; % Post-before-Pre time constant, in sec (20ms in SongAbbott2000, 10ms in MehtaWilson2000)
Pmax = 0.5/100*Imax; % max Pre-before-Post weight change, in percent of Imax (0.5% in SongAbbot2000, 0.6% in MehtaWilson2000)
Dmax = -0.5/100*Imax; % max Post-before-Pre weight change, in percent of Imax (0.525% in SongAbbott2000, 90% of Pamp in MehtaWilson2000)

% output PF 
Nbin = 50; % Number of bins in which the length of the track is divided

%% Inputs
%Place fields (FR = gaussian function of position x)
x = 1:L; % track positions binned by cm
PFcom = linspace(1, L, N); %1:L/N:L; gives floating point issues for large N % PF centers spread uniformly along the track, with L/N spacing in cm
% PFsd = PFwidth/6; % for gaussian, 4sd = 95%, 6sd = 99.7% of data
FRin = zeros(N,length(x));
for i = 1:N
    FRinNorm(i,x) = exp(-0.5*(x-PFcom(i)).^2/PFsd^2); %normalized from 0 to 1
    FRin(i,x) = PFamp*exp(-0.5*(x-PFcom(i)).^2/PFsd^2); % real FR 
end
meanFRin = max(mean(FRin,2)); % mean FR of the middle input neuron, in Hz

figure
imagesc(FRin(round(N/2),:)); 
xlabel('position (cm)');
colormap(jet(256)); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Firing rate (Hz)';
box off; 
% title(['max(meanFRin) = ' num2str(meanFRin) ' Hz'])

figure
imagesc(FRin); 
xlabel('position (cm)'); ylabel('input neurons');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Firing rate (Hz)';
box off; axis square 
title(['max(meanFRin) = ' num2str(meanFRin) ' Hz'])

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
InRasters(InRasters>=1) = 1; 
InRasters(InRasters == 0) = 0; 
InRastersTF = logical(InRasters);
for i = 1:N
    Spiketimes_in{i} = Trun( InRastersTF(i,:) ); %find positive values = when a spike occured
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

%define intial values and length
I = zeros(1,length(Tlap1)); 
V = I;
SRA = I;
OutRaster = I;
D = I;
P = zeros(N,length(Tlap1));
V(1) = Vrest;


for i = 1:length(Tlap1)-1
    % input current to output neuron, EPSCs modelled as pulse with exponential decay
    I(i+1)  = I(i) - dt*I(i)./Idecay + sum(InRasters(:,i).*W(:,i));
    % I&F output neuron
    if Adaptation == 0
        V(i+1) = V(i) + dt * (Vrest - V(i) + I(i) * Rm) / tau_m; % subthreshold membrane potential
    elseif Adaptation == 1
        V(i+1) = V(i) + dt * (Vrest - V(i) -SRA(i)*(V(i)-Ek) + I(i) * Rm) / tau_m; %LIF with spike rate adaptation
        SRA(i+1) = SRA(i) - dt * SRA(i)/tau_sra; % spike-rate-adaptation variable with exponential decay to 0
    end
    
    if V(i+1) >= Vthr
       SRA(n,i+1) = SRA(i+1) + dSRA;
       V(i) = 0; %to visualize spikes with peak at 0 mV
       V(i+1) = Vreset;
       OutRaster(i) = 1;
    end
        
    % STDP variables + update Synaptic weights
    D(i+1) = D(i) - dt*D(i)/Ddecay + OutRaster(i); % output spiketrain convolved with Post-before-Pre update rule, in %
    for n = 1:N % for each presynaptic neuron
        P(n,i+1) = P(n,i) - dt*P(n,i)/Pdecay + InRasters(n,i); % input spiketrains convolved with Pre-before-Post update rule, in %
        W(n,i+1) = W(n,i) + OutRaster(i)*P(n,i)*Pmax + InRasters(n,i)*D(i)*Dmax; % weight update
        if W(n,i+1) < 0
            W(n,i+1) = 0;
        elseif W(n,i+1) > Imax
            W(n,i+1) = Imax;
        end           
    end    
end
Spiketimes_out{1} = Tlap1(find(OutRaster));
SpikeLoc_out = Lap1(find(OutRaster)); %spikes location on track
meanFR_Lap1 = sum(OutRaster)./EndTime;

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
title(['meanFR: ' num2str(meanFR_Lap1) ' Hz'])

figure %STDP variables
InputNum = 66;
x1 = 11.5;
x2 = 12;

subplot(3,1,1)
plot(Tlap1, D, 'r'); hold on
plot(Tlap1, OutRaster, 'k');
xlabel('time (s)'); 
xlim([x1 x2]);
title('output raster + Post-before-Pre variable')
box off;
subplot(3,1,2)
plot(Tlap1, P(InputNum,:), 'r'); hold on
plot(Tlap1, InRasters(InputNum,:), 'k');
xlim([x1 x2]);
title(['input raster #' num2str(InputNum) '  + Pre-before-Post variable'])
box off;
subplot(3,1,3)
plot(Tlap1, W(InputNum,:)./Imax, 'k'); hold on
xlim([x1 x2]);
xlabel('time (s)'); 
ylabel('Norm. synaptic strength')
title(['input #' num2str(InputNum)] )
box off;

figure % Weight changes over time
subplot(2,1,1)
imagesc(W.*10^12);
xlabel('time'); ylabel('input neuron');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Synaptic weight (pA)';
box off;
subplot(2,1,2)
plot(Trun, W(60,:).*10^12, 'r'); hold on
plot(Trun, W(69,:).*10^12, 'k'); hold on
xlabel('time (s)'); ylabel('input current (pA)')
box off; 
axis square;







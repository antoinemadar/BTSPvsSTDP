clear; close all

%% Parameters
%Params for inputs
N = 100; % number of input neurons
L = 300; % length of track in cm
PFsd = 18; % cm (~median in Can data) 
PFamp = 10; %peak FR in Hz
Nlaps = 30; % number of laps
period = 20; % lap duration, in sec
dt = 0.001; % time resolution in sec

% Synapses params
Imax = 85e-12; % max synaptic weight, a current in Amps. % NolanWiles2011 uses 10-226pA. Mehta uses 1.5-4nA. SongAbbot use dimensionless conductances with G between 0 and Gmax, but not sure about the initial weight matrix. 
Idecay = 10e-3; % EPSC time constant, in sec
Wsd = 10; %standard deviation of initial synaptic weight vector (for gaussian connectivity)
maxW = 1; % from 0 to 1, proportion of Imax defining the initial maximum weight

%Integrate and Fire params for output neuron
Rm = 100e6; % membrane resistance, in Ohm (100MOhm in SongAbbot2000)
tau_m = 20e-3 ; %membrane time constant, in sec (25ms in YuShouval2006, 20ms in SongAbbot2000)
Vrest = -70e-3 ; % resting membrane potential, in Volt (-70mV in SongAbbot2000, -60mV in YuShouval2006)
Vthr = -54e-3; % spike threshold, in Volt (-54mV in SongAbbot2000)
Vreset = -60e-3; % reset potential after spike, in Volt (-60mV in SongAbbot2000) 

% Plasticity params
CSproba = 0.018; % proba of an output spike to be a plateau potential/complex spike (1.8% on average in Bittner et al. 2017)
Pdecay = 1.31; % Pre-before-Post time constant, in sec (1.31s on average in Bittner et al. 2017)
Ddecay = 0.69; % Post-before-Pre time constant, in sec (0.69s on average in Bittner et al. 2017)
Pamp = 0.045*Imax; % peak weight change, in Amps (3*baseline EPSP for pairings with 10 EPSPs at 10Hz, in Bittner et al. 2017)
HomeoNorm = 1; % 0 if no homeostatic normalization rule, 1 if yes.

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

% segment in laps
NewLapIdx = find(Run==0);

%% input-output dynamics using Euler's forward method to solve ODEs
I = zeros(1,length(Trun)); 
V = I;
D = I;
OutRaster = I;
CS = I;
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
       CS(i) = binornd(1, CSproba); %Bernoulli trial (0 or 1) with p(success) = CSproba, to model plateau ptential (aka Complex Spike)occurrence
    end
    % BTSP variables + update Synaptic weights W
    D(i+1) = D(i) - dt*D(i)/Ddecay + CS(i); % train of output complex spikes convolved with Post-before-Pre update rule, in %
    for n = 1:N % for each presynaptic neuron
        P(n,i+1) = P(n,i) - dt*P(n,i)/Pdecay + InRasters(n,i); % input spiketrains convolved with Pre-before-Post update rule, in %
        Pact(n,i+1) = CS(i)*P(n,i)*Pamp;
        Dact(n,i+1) = InRasters(n,i)*D(i)*Pamp;
        W(n,i+1) = W(n,i) + Pact(n,i+1) + Dact(n,i+1); % weight update
%         if W(n,i+1) < 0
%             W(n,i+1) = 0;
%         elseif W(n,i+1) > Imax % this upper bound can be deleted if implementing homeostatic rule below. But that will allow higher max FR and thus increase likelihood of complex spikes.
%             W(n,i+1) = Imax;
%         end           
    end 

    if HomeoNorm == 1 % implement heterosynaptic competition or not
       W2(:,i+1) = W(:,i+1);
       W(:,i+1) = W(:,i+1).*sum(W(:,1))./sum(W(:,i+1)); % Homeostatic synaptic normalization, multiplicative method
        
    end

end
Spiketimes_outAll = Trun(find(OutRaster));
SpikeLoc_outAll = Run(find(OutRaster)); %spikes location on track
CSloc_All = Run(find(CS)); % CSs location on track

% segment output raster by lap + compute spatial firing rate

NewLapTimes = Trun(NewLapIdx);
Track_dbin = L./Nbin; % bin size, in cm
TrackBins = Lap1(1):Track_dbin:Lap1(end); % Spatial bin edges
BinCenters = L/(2*Nbin):L/Nbin:L; % in cm
for lap = 1:Nlaps
    if lap == Nlaps
    OutRaster_laps{lap} = OutRaster(NewLapIdx(lap):end);
    LapN = Run(NewLapIdx(lap):end);
    CS_laps{lap} = CS(NewLapIdx(lap):end);
    else
    OutRaster_laps{lap} = OutRaster(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    LapN = Run(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    CS_laps{lap} = CS(NewLapIdx(lap):NewLapIdx(lap+1)-1);
    end
    Spiketimes_out{lap} = Tlap1(find(OutRaster_laps{lap}));
    SpikeLoc_out{lap} = Lap1(find(OutRaster_laps{lap}));
    CStimes{lap} = Tlap1(find(CS_laps{lap}));
    CSloc{lap} = Lap1(find(CS_laps{lap}));
    % compute place field
    SpikeCountOut(lap, :) = histcounts(SpikeLoc_out{lap}, TrackBins); % number of spikes in each spatial bin
    CSbin(lap,:) = histcounts(CSloc{lap}, TrackBins); 
    TimeInBins(lap,:) = dt*histcounts(LapN, TrackBins); % compute time spent (in sec) in each spatial bin (no assumption on trajectories)
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

%% figures

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

% figure % output spiketrains + spatial FR by lap
% subplot(2,1,1)
% [spikex, spikey] = makedisplayrasters(Spiketimes_out, 0);
% displaydisplayrasters(spikex, spikey); hold on
% grid off;
% xlim([0,Tlap1(end)]); xlabel('time, (s)')
% ylabel('lap')
% title('output neuron')
% subplot(2,1,2)
% [spikex, spikey] = makedisplayrasters(CStimes, 0);
% displaydisplayrasters(spikex, spikey); hold on
% grid off;
% xlim([0,Tlap1(end)]); xlabel('time, (s)')
% ylabel('lap')
% title('Complex Spikes')

figure
[row_lap,column_bin] = find(CSbin);
imagesc(SpatialFRout); hold on
scatter(column_bin, row_lap, 10, [0 1 1], 'filled'); hold on
scatter(COMbin,1:Nlaps, 10, [1 0 0], 'filled');
legend('Complex spike','COM', 'Location', 'BestOutside');
xlabel('spatial bins'); ylabel('lap');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Firing rate (Hz)';
box off; axis square

figure %STDP variables on last lap
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
box off;
subplot(2,1,2)
imagesc([0 Nlaps*period],[1 N], W2.*10^12); hold on
scatter(NewLapTimes, (N/2)*ones(1,Nlaps), 10, [1 0 0], '+');
xlabel('time (s)'); ylabel('input neuron');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Synaptic weight (pA)';
box off;

% plot(Trun, W(25,:).*10^12, 'r'); hold on
% plot(Trun, W(15,:).*10^12, 'k');
% scatter(NewLapTimes, W(1,NewLapIdx).*10^12, 10, [0 0 1], '+'); 
% legend('input neuron 25', 'input neuron 15', 'input 1 (new lap start)')
% xlabel('time (s)'); ylabel('input current (pA)')
% title('Synaptic weight of input neuron #1');
% box off;

figure
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

figure
plot(1:Nlaps, COMloc, 'k-'); hold on
plot(1:Nlaps, COMtraj, 'b-'); hold on
yline(COMloc(1),'r');
legend('COM','lin reg','COM #1')
ylim([0 300]);
xlabel('lap'); ylabel('COM position (cm)');
title({'Place Field COM trajectory, slope =' num2str(b(2)) 'pval =' num2str(stats(3))})

% plot(Tlap1, P(1,NewLapIdx(1):NewLapIdx(2)-1), 'r'); hold on
% for i = 1:length(Spiketimes_out{1})
%     xline(Spiketimes_out{1}(i),'k');
% end
% xlim([10, 12]);
% subplot(2,1,2)
% plot(Tlap1, Pact(1,NewLapIdx(1):NewLapIdx(2)-1), 'r'); hold on
% plot(Tlap1, Dact(1,NewLapIdx(1):NewLapIdx(2)-1), 'b'); hold on
% % plot(Tlap1, InRasters(1,NewLapIdx(1):NewLapIdx(2)-1), 'k'); hold on
% % plot(Tlap1, D(NewLapIdx(1):NewLapIdx(2)-1), 'r'); hold on
% for i = 1:length(Spiketimes_out{1})
%     xline(Spiketimes_out{1}(i),'k'); hold on
% end
% xlim([10, 12]);
% xlabel('time (s)');





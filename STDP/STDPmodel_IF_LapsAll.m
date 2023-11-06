clear; close all

load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA3_SlopeDistrib.mat')

tic 

%% Parameters
%Params for inputs
N = 1000; % number of input neurons
L = 300; % length of track in cm
PFsd = 18; % cm (10 ~ median in our CA3 data, mean = 18) 
PFamp = 10; %peak FR in Hz
Nlaps = 30; % number of laps
period = 20; % lap duration, in sec
dt = 0.001; % time resolution in sec
InShift = 0; % 0 if static input, 1 if shift
Context = 'N'; % N for novel, F for Familiar, NandF for all

% Synapses params
Imax = 12e-12; % max synaptic weight, a current in Amps. % NolanWiles2011 uses 10-226pA. Mehta uses 1.5-4nA. SongAbbot use dimensionless conductances with G between 0 and Gmax, but not sure about the initial weight matrix. 
Idecay = 10e-3; % EPSC time constant, in sec
Wsd = 100; % in neurons. Standard deviation of initial synaptic weight vector (for gaussian connectivity)
maxW = 1; % from 0 to 1, proportion of Imax defining the initial maximum weight

%Integrate and Fire params for output neuron
Rm = 100e6; % membrane resistance, in Ohm (100MOhm in SongAbbot2000)
tau_m = 20e-3 ; %membrane time constant, in sec (25ms in YuShouval2006, 20ms in SongAbbot2000)
Vrest = -70e-3 ; % resting membrane potential, in Volt (-70mV in SongAbbot2000, -60mV in YuShouval2006)
Vthr = -54e-3; % spike threshold, in Volt (-54mV in SongAbbot2000)
Vreset = -60e-3; % reset potential after spike, in Volt (-60mV in SongAbbot2000) 
% Icap = 350e-12; % in Amp (simple way to enforce shunting inhib)
% K = 0.7; % divisive inhibition constant as in Yu Shouva2006 (value unknown)

Adaptation = 0; % 1 or 0, to implement spike-rate adaptation or not
Ek = -70e-3; % K+ equilibrium potential, in Volt
dSRA = 0.06; % increment value of spike-rate adaptation variable modeling adapting K+ leak current (cf DayanAbbott fig 5.6 and eq 5.13)
tau_sra = 100e-3; % time constant for exponential decay of SRA variable

% Plasticity params
Pdecay = 20e-3; % Pre-before-Post time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Ddecay = 20e-3; % Post-before-Pre time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Pamp = 0.005*Imax; % peak Pre-before-Post weight change (0.5% in SongAbbot2000, 0.6pA with no Imax in MehtaWilson2000)
Damp = -0.005*Imax; % peak Post-before-Pre weight change (0.525% in SongAbbott2000, 90% of Pamp in MehtaWilson2000, which didn't have a maximum weight)

capWeights = 1; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 

HomeoNorm = 0; % 0 if no homeostatic normalization rule, 1 if yes.

WUdyn = 0; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Wrest = 0;
tau_wu = 5; % weight update time constant, in sec

% output PF 
Nbin = 50; % Number of bins in which the length of the track is divided

%% Inputs

% for dynamic inputs: slope distributions
PslopeBinCenters = [-1.6 -1.4, -1.2 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6]; %in cm/lap
if ismember('F', Context)
    Pslope = CA3_F_SlopeHist.Values;
elseif ismember('N', Context)
    Pslope = CA3_N_SlopeHist.Values;
elseif ismember('FandN', Context)
    Pslope = CA3_FandN_SlopeHist.Values;
end

% randomly select slope for each input neuron according to bin proba
InSlopes = randsample(PslopeBinCenters,N,true,Pslope); 

%Place fields (FR = gaussian function of position x)
x = 1:L; % track positions binned by cm
PFcom = linspace(1,L,N); %1:L/N:L; % PF centers spread uniformly along the track, with L/N spacing in cm
% PFsd = PFwidth/6; % for gaussian, 4sd = 95%, 6sd = 99.7% of data
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

% figure % for static inputs
% imagesc(FRin); 
% xlabel('position (cm)'); ylabel('input neurons');
% colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
% c = colorbar; c.Label.String = 'Firing rate (Hz)';
% box off; axis square 

% figure % example input neuron with dynamic field
% [Fval, Nex1] = max(InSlopes);
% % [Bval, Nex2] = min(InSlopes);
% Nex = 50;
% FRbyLap = squeeze(FRin(:,Nex,:));
% imagesc(FRbyLap)
% xlabel('position (cm)'); ylabel('lap');
% colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
% c = colorbar; c.Label.String = 'Firing rate (Hz)';
% box off; axis square

%Running model with constant speed
speed = L/period; % cm/sec
EndTime = Nlaps*L/speed; %in sec
Trun = 0:dt:EndTime; Trun = Trun(1:end-1); % ms resolution
Tlap1 = 0:dt:period; Tlap1 = Tlap1(1:end-1);
Lap1 = linspace(0,L,period/dt); % positions during 1 lap
Run = repmat(Lap1,[1 Nlaps]);%positions during the whole run

% figure
% plot(Trun(1:end), Run, 'k');
% xlabel('time (s)'); ylabel('position (cm)')
% box off; axis square

% Lap indices and times for given trajectory
NewLapIdx = find(Run==0);
NewLapTimes = Trun(NewLapIdx);

%convert trajectory into firing rate
   for i = 1:length(Run)
   [minval, idx] = min(abs(x-Run(i))); %find closest spatial bin in x for each position in trajectory
       if InShift == 0
        FRrun(:,i) = FRin(:,idx); % Firing rate at that position for each neuron
       elseif InShift == 1  
        [val, l] = min(abs(NewLapIdx-i)); % what lap is it?
        FRrun(:,i) = FRin(l,:,idx);
       end
   end

figure
imagesc(FRrun);
xlabel('time (s)'); ylabel('input neuron')
% xlim([500000 600000])
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Firing rate (Hz)';
box off; axis square

% nonhomogenous Poisson spike generator along trajectory
InRasters = poissrnd(FRrun.*dt); %matrix of size (N,length(Run)) containing, for each input neuron, the number of spikes in each 1ms bin drawn from Poisson distribution
InRasters(InRasters>=1) = 1; 
InRasters(InRasters == 0) = 0; 
InRastersTF = logical(InRasters);
for i = 1:N
    Spiketimes_in{i} = Trun( InRastersTF(i,:) ); %find positive values = when a spike occured
end

% figure
% [spikex, spikey] = makedisplayrasters(Spiketimes_in, 0);
% displaydisplayrasters(spikex, spikey); hold on
% grid off
% axis square
% xlabel('time, (s)')
% ylabel('input neurons')

% gaussian connectivity with max in center of track, akin to Mehta and Wilson 2000
% Wsd = N/6; %standard deviation of synaptic weight distribution
NeuronID = 1:N;
W = zeros(N,length(Trun));
Wmax = maxW*Imax;
W(:,1) = Wmax*exp(-0.5*(NeuronID - N/2).^2/Wsd^2); % vector of N rows following a gaussian distribution of Wmax amplitude 
W2 = W;
Wtarget = W;

figure
stem(1:N, W(:,1)*10^12, 'k');
xlabel('Input Neurons'); ylabel('Synaptic Weights (pA) ')
box off; axis square

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

    if HomeoNorm == 1
    I(i+1)  = I(i) - dt*I(i)./Idecay + sum(InRasters(:,i).*W2(:,i));
    else
    I(i+1)  = I(i) - dt*I(i)./Idecay + sum(InRasters(:,i).*W(:,i)); 
    end
%             if exist(Icap, 'var') 
%                 if I(i+1) > Icap % implements shunting inhib rather than bounding weights
%                 I(i+1) = Icap;
%                 end
%             end
    % I&F output neuron
    if Adaptation == 0
        V(i+1) = V(i) + dt * (Vrest - V(i) + I(i) * Rm) / tau_m; % subthreshold membrane potential
    elseif Adaptation == 1
        V(i+1) = V(i) + dt * (Vrest - V(i) -SRA(i)*(V(i)-Ek) + I(i) * Rm) / tau_m;
        SRA(i+1) = SRA(i) - dt * SRA(i)/tau_sra; % spike-rate-adaptation variable with exponential decay to 0
    end
    
    if V(i+1) >= Vthr
       SRA(i+1) = SRA(i+1) + dSRA;
       V(i) = 0; %to visualize spikes with peak at 0 mV
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
            W(n,i+1) = W(n,i) + dt * (Wrest - W(n,i) + Wtarget(n,i+1)) / tau_wu;    
        end
        
        % Hard bounds on weights
        if W(n,i+1) < 0
            W(n,i+1) = 0;
        elseif W(n,i+1) > Imax && capWeights == 1 % as in SongAbbott2000, weights are capped to prevent runaway plasticity and FR
            W(n,i+1) = Imax;
        end

    end
    
    if HomeoNorm == 1 % implement heterosynaptic competition or not
       W2(:,i+1) = W(:,i+1).*sum(W(:,1))./sum(W(:,i+1)); % Homeostatic synaptic normalization, multiplicative method
    end


end
Spiketimes_outAll = Trun(find(OutRaster));
SpikeLoc_outAll = Run(find(OutRaster)); %spikes location on track

% compute output spatial firing rate by lap
Track_dbin = L./Nbin; % bin size, in cm
TrackBins = Lap1(1):Track_dbin:Lap1(end); % Spatial bin edges
BinCenters = L/(2*Nbin):L/Nbin:L; % in cm
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
    % compute place field
    SpikeCountOut(lap, :) = histcounts(SpikeLoc_out{lap}, TrackBins); % number of spikes in each spatial bin
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

toc

%% figures

figure % input current and Vm on last lap
subplot(3,1,1)
plot(Tlap1, I(NewLapIdx(end):end).*10^12, 'k')
xlabel('time (s)'); ylabel('input current (pA)')
box off;
title('input current on last lap');
subplot(3,1,2)
plot(Tlap1, V(NewLapIdx(end):end).*10^3, 'k')
xlabel('time (s)'); ylabel('membrane potential (mV)')
box off;
title('Vm on last lap (LIF model)');
subplot(3,1,3)
plot(Tlap1, SRA(NewLapIdx(end):end), 'k')
xlabel('time (s)'); ylabel('SRA variable (no unit)')
box off;
title('SRA on last lap (LIF model)');

% figure % output spiketrains 
% [spikex, spikey] = makedisplayrasters(Spiketimes_out, 0);
% displaydisplayrasters(spikex, spikey); hold on
% grid off; axis square;
% xlim([0,Tlap1(end)]); xlabel('time, (s)')
% ylabel('lap')
% title('output neuron')

figure %spatial FR by lap
imagesc(SpatialFRout); hold on
scatter(COMbin,1:Nlaps, 30, [1 0 0], 'filled');
xlabel('spatial bins'); ylabel('lap');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Firing rate (Hz)';
box off; axis square;

figure %spatial FR lap 1 vs last lap
subplot(1,2,1)
plot(SpatialFRout(1,:), 'k'); hold on
plot(SpatialFRout(end,:), 'r'); hold on
xlabel('spatial bins'); ylabel('Firing rate (Hz)');
box off; axis square;
subplot(1,2,2)
plot(SpatialFRout(1,:)./max(SpatialFRout(1,:)), 'k'); hold on
plot(SpatialFRout(end,:)./max(SpatialFRout(end,:)), 'r'); hold on
xlabel('spatial bins'); ylabel('norm FR relative to peak');
box off; axis square;

figure %STDP variables for a given lap and input neuron
LapNum = 1;
InputNum = 66;
x1 = 0; % in s
x2 = 20; % in s 

subplot(3,1,1)
plot(Tlap1, D(NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'r'); hold on
plot(repmat(Spiketimes_out{LapNum},2,1), repmat([0; 1], 1, length(Spiketimes_out{LapNum})), 'k');
% plot(Tlap1, OutRaster_laps{LapNum}, 'k');
xlim([x1 x2]);
xlabel('time (s)'); 
title('output raster + Post-before-Pre variable')
box off;
subplot(3,1,2)
plot(Tlap1, P(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'r'); hold on
plot(Tlap1, InRasters(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'k');
xlim([x1 x2]);
xlabel('time (s)'); 
title(['input raster #' num2str(InputNum) '+ Pre-before-Post variable'])
box off;
subplot(3,1,3)
plot(Tlap1, W(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./Imax, 'r'); hold on
if WUdyn == 1
plot(Tlap1, Wtarget(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./Imax, 'k'); hold on
end
xlim([x1 x2]);
xlabel('time (s)'); 
ylabel('Norm. synaptic strength')
title(['input #' num2str(InputNum)] )
box off;

figure %STDP variables on last lap
subplot(3,1,1)
plot(Tlap1, D(NewLapIdx(end):end), 'r'); hold on
plot(Tlap1, OutRaster_laps{Nlaps}, 'k');
xlim([x1 x2]);
xlabel('time (s)'); 
title('output raster + Post-before-Pre weight update')
box off;
subplot(3,1,2)
plot(Tlap1, P(InputNum,NewLapIdx(end):end), 'r'); hold on
plot(Tlap1, InRasters(InputNum,NewLapIdx(end):end), 'k');
xlim([x1 x2]);
xlabel('time (s)'); 
title(['input raster #' num2str(InputNum) '+ Pre-before-Post weight update'])
box off;
subplot(3,1,3)
plot(Tlap1, W(InputNum,NewLapIdx(end):end)./Imax, 'r'); hold on
if WUdyn == 1
    plot(Tlap1, Wtarget(InputNum,NewLapIdx(end):end)./Imax, 'k'); hold on
end
xlim([x1 x2]);
xlabel('time (s)'); 
ylabel('Norm. synaptic strength')
title(['input #' num2str(InputNum)] )
box off;
% print('-vector','-dpdf','WUdynamics.pdf')

figure %STDP variables for a given lap and input neuron, organized differently, for Figure 2 illustration
LapNum = 29;
InputNum = 66;
x1 = 12.45; % in s
x2 = 13; % in s 
subplot(5,1,1) % Pre before Post variable (Potentiation)
    plot(Tlap1, P(InputNum, NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'Color', [1    0.5    0]); hold on
    xlim([x1 x2]);
    set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
    box off;
subplot(5,1,2) % input spiketrain
    SpkTimesIN = Spiketimes_in{InputNum}(Spiketimes_in{InputNum}>Trun(NewLapIdx(LapNum)) & Spiketimes_in{InputNum}<Trun(NewLapIdx(LapNum+1)-1));
    SpkTin = SpkTimesIN-Trun(NewLapIdx(LapNum)); % center time on the current lap
    plot(repmat(SpkTin,2,1), repmat([0; 1], 1, length(SpkTin)), 'k');hold on
    xlim([x1 x2]);
    set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
    box off
subplot(5,1,3) % output spiketrain
    plot(repmat(Spiketimes_out{LapNum},2,1), repmat([0; 1], 1, length(Spiketimes_out{LapNum})), 'k');
    % plot(Tlap1, OutRaster_laps{LapNum}, 'k');
    xlim([x1 x2]);
    set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
    box off;
subplot(5,1,4) % depression (Post-before-pre) variable
    plot(Tlap1, D(NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1), 'Color', [0.5 0.2 0.5]); hold on
%     title(['input #' num2str(InputNum)] )
    xlim([x1 x2]);
    set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
    box off;
subplot(5,1,5)
    plot(Tlap1, W(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./Imax, 'r'); hold on
    if WUdyn == 1
    plot(Tlap1, Wtarget(InputNum,NewLapIdx(LapNum):NewLapIdx(LapNum+1)-1)./Imax, 'k'); hold on
    end
    xlim([x1 x2]);
    xlabel('time (s)');
    set(gca, 'YTick', [], 'YTickLabel', [])
    ylabel({'Norm.', 'synaptic', 'strength'})
    box off;

figure % Weight changes over time, full run
% subplot(2,1,1)
imagesc([0 Nlaps*period],[1 N], W.*10^12); hold on
scatter(NewLapTimes, (N/2)*ones(1,Nlaps), 10, [1 0 0], '+');
xlabel('time (s)'); ylabel('input neuron');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Synaptic weight (pA)';
box off; 
axis square;
% subplot(2,1,2)
% plot(Trun, W(25,:).*10^12, 'r'); hold on
% plot(Trun, W(20,:).*10^12, 'k'); hold on
% scatter(NewLapTimes, W(1,NewLapIdx).*10^12, 10, [0 0 1], '+'); 
% legend('input 25', 'input 20', 'input 1: new laps', 'Location', 'BestOutside');
% xlabel('time (s)'); ylabel('input current (pA)')
% box off; 
% axis square;

figure % Weight changes over time, full run, zoomed in
imagesc([0 Nlaps*period],[30 70], W(30:70,:).*10^12); hold on
scatter(NewLapTimes, (N/2)*ones(1,Nlaps), 10, [1 0 0], '+');
xlabel('time (s)'); ylabel('input neuron');
colormap(flipud(gray(256))); %colormap(brewermap(256,'*YlOrRd'));
c = colorbar; c.Label.String = 'Synaptic weight (pA)';
box off; 
axis square;

% current just before 1st spike of InputNum for each lap + net weight change for InputNum on each lap
PtsNum = 50;  % number of datapoints over which the input current is averaged    
InitSpkIdx = find(InRasters(InputNum,NewLapIdx(1):NewLapIdx(2)-1),1);
for l = 1:Nlaps-1
    Ilap = I(NewLapIdx(l):NewLapIdx(l+1)-1);
    Ipre(l) = mean(Ilap(InitSpkIdx-PtsNum:InitSpkIdx-1)); % average current over DataPtsNum data points before initial spike in input raster
    NetWchange(:,l) = W(:,NewLapIdx(l+1)-1) - W(:,NewLapIdx(l)); % net weight change for all inputs on each lap
    WendLap(:,l) = W(:,NewLapIdx(l+1));
end  
figure
subplot(3,1,1)
    plot(1:Nlaps-1, Ipre.*10^12, 'k-'); hold on
    xlabel('Lap number'); ylabel('input current (pA)')
    box off; 
    axis square;
subplot(3,1,2)
    plot( SpatialFRout(1:Nlaps-1,round(InputNum/(N/Nbin))-1), 'b-' ); hold on
    plot( SpatialFRout(1:Nlaps-1,round(InputNum/(N/Nbin))), 'k-' )
    xlabel('Lap number'); ylabel('Output FR (Hz)')
    title('blue: bin before, black: bin of Input Neuron')
    box off; 
    axis square;
subplot(3,1,3)
    plot(1:Nlaps-1, WendLap(InputNum,:).*10^12, 'r-')
    % plot(1:Nlaps-1, NetWchange(InputNum,:).*10^12, 'r-')
    ylim([0 max(WendLap(InputNum,:).*10^12+1)])
    xlabel('Lap number'); ylabel('Weight at end of lap (pA)')
    title(['Input #' num2str(InputNum)])
    box off; 
    axis square;


figure
plot(1:Nlaps-1, Ipre./max(Ipre), 'k-'); hold on
plot( SpatialFRout(1:Nlaps-1,round(InputNum/(N/Nbin)))/max(SpatialFRout(1:Nlaps-1,round(InputNum/(N/Nbin)))), 'g-' ); hold on
plot( SpatialFRout(1:Nlaps-1,round(InputNum/(N/Nbin))-1)/max(SpatialFRout(1:Nlaps-1,round(InputNum/(N/Nbin)-1))), 'b-' ); hold on
plot(1:Nlaps-1, WendLap(InputNum,:)./max(WendLap(InputNum,:)), 'r-')
xlabel('Lap number'); title(['Input #' num2str(InputNum) ' Ipre (Black), FR (green), FRpre (Blue), norm. weight at end of lap (Red), normalized'])
box off; 
axis square;

figure % initial weights vs end weights
subplot(1,2,1)
plot(1:N, W(:,1).*10^12, 'k'); hold on
plot(1:N, W(:,end).*10^12, 'r'); hold on
legend('start', 'end', 'Location', 'Best');
xlabel('input neuron'); ylabel('synaptic weight (pA)')
box off; 
axis square;
subplot(1,2,2)
plot(1:N, W(:,1)./max(W(:,1)), 'k'); hold on
plot(1:N, W(:,end)./max(W(:,end)), 'r'); hold on
legend('start', 'end', 'Location', 'Best');
xlabel('input neuron'); ylabel('synaptic weight relative to peak')
box off; 
axis square;

figure
    subplot(2,2,1)
        plot(1:Nlaps, meanFRout_lap, 'r+'); hold on
        plot(1:Nlaps, maxFRout_lap, 'k-');
        xlabel('lap'); ylabel('Firing Rate (Hz)')
%         legend('mean','max', 'Location', 'BestOutside');
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
%         legend('mean PF', 'COM', 'track center', 'Location', 'BestOutside');
%         title({'SD = ' num2str(SD_meanPF) ', Skew = ' num2str(Skew_meanPF)})
        box off; axis square;

figure % Shift: COM regression
plot(1:Nlaps, COMloc, 'k-'); hold on
plot(1:Nlaps, COMtraj, 'b-'); hold on
yline(COMloc(1),'r');
legend('COM','lin reg','COM #1')
ylim([0 300]);
xlabel('lap'); ylabel('COM position (cm)');
title({'Place Field COM trajectory, slope =' num2str(b(2)) 'pval =' num2str(stats(3))})
box off; axis square;

figure % Shift Zoom
plot(1:Nlaps, COMloc, 'k-'); hold on
plot(1:Nlaps, COMtraj, 'b-'); hold on
yline(COMloc(1),'r');
legend('COM','lin reg','COM #1')
ylim([125 175]);
xlabel('lap'); ylabel('COM position (cm)');
title({'Place Field COM trajectory, slope =' num2str(b(2)) 'pval =' num2str(stats(3))})
box off; axis square;






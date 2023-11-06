%simulation of in vitro demonstration of BTSP rule as in Bittner et al. 2017
clear 
close all
%% Parameters

% Synapses params
Params.Idecay = 10e-3; % EPSC time constant, in sec
Wsyn = 77e-12 % 63.75e-12; % synaptic weight of the stimulated input fiber (or ensemble of fibers), in Amps. This stimulation-evoked current was adjusted to get 2mV amplitude for initial EPSP 
Imax = 85e-12; % max synaptic weight, in Amps (used in Sheffield or Milstein's style simulations). Not relevant in this experiment, but Pamp can be expressed in terms of it, and the parameter is required in BTSPsimple_LIFadapt_ODEs in case weight capping is implemented
% Vsyn = 2e-3; % in V. Stim Amplitude before plasticity

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
Params.Pamp = 3.8e-12 %0.06*Wsyn; %0.045*Imax; % peak weight change, in Amps (3*baseline EPSP for pairings with 10 EPSPs repeated 5 times, in Bittner et al. 2017 => 3/50=0.06)
Params.PostPreBoost = 1.15; % amplitude of weight change for postpre rule which is triggered on CS, not on spikes (so, to be somewhat balanced, it should = 10 so that amplitude is the Bittner2017 amplitude/1CS and not /10spikes like the prepost rule)
Params.capWeights = 0; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 
Params.WUdyn = 0; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Params.Wrest = 0; % Synaptic weight resting state, in Amps (towards which any synapse would relax if Wtarget = 0) 
Params.tau_wu = 5; % weight update time constant, in sec

% Protocol params
pairings = 5; %number of pairings
NumSpikes = 10; % number of spikes
StimFR = 20; % Hz
stimD = 15; % Duration of stim protocol (including input spikes and CS induction). CS will be in the middle of that stim duration. Choose duration > 10 to add single input at beginning without affecting pairings 
ITI = 15; % Interval between stim protocols, in sec
dt = 0.001; % time resolution in sec
Delay = [-4000, -3000, -2000, -1000, -750, -500, -250, -100, -10, -1, 1, 10, 100, 250, 500, 750, 1000, 2000, 3000, 4000]; % delay between last stim and CS, in ms. Negative values correspond to post-before-pre condition.

%% Protocol

StimLen = stimD/dt;
Tsweep = dt:dt:stimD+ITI; % ms resolution
Trun = dt:dt:pairings*(stimD+ITI);
TswLen = length(Tsweep);

% CS trace
CS = zeros(1,StimLen); 
CSidx = round(StimLen/2); % middle of the stim protocol
CS(CSidx) = 1; % put CS in the middle of the stim protocol
CS = [CS, zeros(1,ITI/dt)]; % pad with ITI
CSall = repmat(CS, 1, pairings); % collate all pairings+ITI

W = zeros(length(Delay),length(CSall));

InRasters = W;
I = W; V = W; V2 = W; 
OutRaster= W; 
D = W; P = W; Pact = W; Dact = W; SRA = W; Wtarget = W;
Params.dt = dt; Params.Imax = Imax;

for d = 1:length(Delay)

    % make inputs
    InRaster = zeros(size(Tsweep));
    SpiketrainD = NumSpikes/StimFR; %in sec
    ISI = 1/StimFR; % = SpiketrainD/NumSpikes, in sec
    SpikeIdx = 1:ISI/dt:SpiketrainD/dt;
    
    closestSpkIdx = CSidx - Delay(d); % if delay is < 0 => stim spikes will occur after CS, i.e. testing post-pre rule
    
        if closestSpkIdx < CSidx % pre before post
            InRaster(closestSpkIdx-SpikeIdx+1) = 1;
        elseif closestSpkIdx > CSidx % post before pre
            InRaster(SpikeIdx+closestSpkIdx-1) = 1;
        end
    
    InRasters(d,:) = repmat(InRaster, 1, pairings);

        % add 1 stim spike at beginning and at the end (1s before, to let EPSP unravel) to measure pre and post-induction EPSP amplitude
    InRasters(d, 1) = 1; 
    InRasters(d, end-1/dt) = 1;
    clear InRaster

    % Synaptic weight (in pA)
    W(d,1) = Wsyn; %Vsyn/Params.Rm; % in Amp. Only 1 synapse. 
    
    % run LIF + plasticity model
    [I(d,:), V(d,:), V2(d,:), CSall, OutRaster(d,:), D(d,:), P(d,:), Dact(d,:), Pact(d,:), SRA(d,:), W(d,:), Wtarget(d,:)]  = BTSPsimple_LIFadapt_ODEs(Params, W(d,:), InRasters(d,:), CSall) ;
    
    % get post-induction EPSP amplitude and normalize to pre-induction amplitude (i.e. Vsyn)
    EPSPpost(d) = max(V2(d, end-1/dt:end)); %in Volts
    EPSPpre(d) = max(V2(d, 1:100)); % in Volts
%     EPSPmax(d) = max(V2(d,:)); %in Volts
    
    EPSPampl(d).pre = 10^3*(EPSPpre(d)-Params.Vrest); %in mv
    EPSPampl(d).post = 10^3*(EPSPpost(d)-Params.Vrest); %in mV
    normEPSPampl(d) = EPSPampl(d).post/EPSPampl(d).pre; % amplitude of post-induction EPSP normalized to pre-induction EPSP

end

%% figure

Didx = 10; % Delay index to plot. Corresponds to Delay(D) in ms

% InRasters and CS raster
figure
    plot(repmat(Trun(logical(InRasters(Didx,:))),2,1), repmat([0; 1], 1, length(find(InRasters(Didx,:)))), 'k'); hold on
    scatter(Trun(logical(CSall)), ones(size(pairings)).*0.5, 'oc', 'filled')
    xlabel('time (s)'); 
    title(['input raster: delay = ' num2str(Delay(Didx)) ' ms'])
    box off;

figure % zoom on 1st pairing + align time on 1st CS
TrunCScentered = Trun - Trun(CSidx); 
    plot(repmat(TrunCScentered(logical(InRasters(Didx,:))),2,1), repmat([0; 1], 1, length(find(InRasters(Didx,:)))), 'k'); hold on
    scatter(TrunCScentered(logical(CSall)), ones(size(pairings)).*0.5, 'oc', 'filled')
    xlim([-1 1])
    xlabel('time (s)');
    title(['input raster, pairing 1: delay = ' num2str(Delay(Didx)) ' ms'])
    box off;

figure % Vm on first and last pairing
CStimes = Trun(logical(CSall));
TrunCScentered2 = Trun - CStimes(end); %time centered on last CS
    plot(TrunCScentered, V(Didx,:).*10^3, 'k'); hold on
    plot(TrunCScentered2, V2(Didx,:).*10^3, 'r'); hold on
    xlim([-1 1])
    xlabel('time (s)'); ylabel('membrane potential (mV)')
    title(['1st (black) and last (red) pairing: delay = ' num2str(Delay(Didx)) ' ms'])
    box off;

figure %post-induction EPSP at the end
lastEPSPtime = Trun(end-1/dt);
    plot(Trun-lastEPSPtime, V2(Didx,:).*10^3, 'k'); hold on
    xlim([-1 1])
    xlabel('time (s)'); ylabel('membrane potential (mV)')
    title(['post-induction EPSP: delay = ' num2str(Delay(Didx)) ' ms'])
    box off;

% post before Pre variable on first pairing, to check how it sums and
% whether decay appropriately falls close to 0 after 5 or 6 s. 
figure
plot(TrunCScentered, P(Didx,:), 'r'); hold on
plot(TrunCScentered, InRasters(Didx,:), 'k');
xlim([-1 15])
xlabel('time (s)'); 
title('input raster #1 + Pre-before-Post weight update')
box off;

% synaptic weight vs time
figure
plot(Trun, W(Didx,:)*10^12, 'k'); hold on
plot(repmat(Trun(logical(CSall)),2,1), repmat([0; max(W(Didx,:))*10^12], 1, length(find(CSall))), 'c');
xlabel('time (s)');
ylabel('synaptic weight (pA)')
title(['delay = ' num2str(Delay(Didx)) ' ms'])

% plot normalized post-induction EPSP as a function of delay (with the rule overlaid)

delay1 = -5:dt:0; % in sec
delay2 = 0:dt:5; % in sec
Wd = (max(normEPSPampl(Delay<0))-1)*exp(-delay2/Params.Ddecay)+1;
Wp = (max(normEPSPampl(Delay>0))-1)*exp(delay1/Params.Pdecay)+1;

figure
plot(delay2, Wd, 'b--' ); hold on
plot(delay1, Wp, 'r--' ); hold on
scatter(-Delay.*dt, normEPSPampl, 'ok'); hold on
% scatter(Delay.*dt, normEPSPampl, 'ok', 'filled', 'MarkerFaceAlpha', 0.6); hold on
% scatter(Delay(Delay<0).*dt, normEPSPampl(Delay<0), 'ob', 'filled'); hold on
% scatter(Delay(Delay>0).*dt, normEPSPampl(Delay>0), 'or', 'filled')
ylabel('norm. EPSP amplitude')
xlabel('CS (post) - Stim (pre) delay (s)')
% ylim([1 max(normEPSPampl)+0.5])

tiledlayout(1,4, 'TileSpacing','Compact','Padding','Compact')

nexttile
    plot(repmat(Trun(logical(InRasters(Didx,:))),2,1), repmat([0; 1], 1, length(find(InRasters(Didx,:)))), 'k'); hold on
    scatter(Trun(logical(CSall)), ones(size(pairings)).*0.5, 'oc', 'filled')
    xlabel('time (s)'); 
    xlim([-5 150])
    title(['input raster: delay = ' num2str(Delay(Didx)) ' ms'])
    box off; axis square;

nexttile % zoom on 1st pairing + align time on 1st CS
    plot(repmat(TrunCScentered(logical(InRasters(Didx,:))),2,1), repmat([0; 1], 1, length(find(InRasters(Didx,:)))), 'k'); hold on
    scatter(TrunCScentered(logical(CSall)), ones(size(pairings)).*0.5, 'oc', 'filled')
    xlim([-1 1])
    xlabel('time (s)');
    title(['input raster, pairing 1: delay = ' num2str(Delay(Didx)) ' ms'])
    box off; axis square;

nexttile % Vm on first and last pairing
    plot(TrunCScentered, V(Didx,:).*10^3, 'k'); hold on
    plot(TrunCScentered2, V2(Didx,:).*10^3, 'r'); hold on
    xlim([-1 1])
    xlabel('time (s)'); ylabel('membrane potential (mV)')
    title(['1st (black) and last (red) pairing: delay = ' num2str(Delay(Didx)) ' ms'])
    box off; axis square;

nexttile
    plot(delay2, Wd, '--', 'Color', [0.5 0.2 0.5] ); hold on
    plot(delay1, Wp, '--', 'Color', [1    0.5    0] ); hold on
%     scatter(-Delay.*dt, normEPSPampl, 'ok'); hold on
    scatter(-Delay(Delay<0).*dt, normEPSPampl(Delay<0), 'o', 'MarkerEdgeColor', [0.5 0.2 0.5], 'MarkerFaceColor', [0.5 0.2 0.5], 'MarkerFaceAlpha', 0.6); hold on
    scatter(-Delay(Delay>0).*dt, normEPSPampl(Delay>0), 'or', 'MarkerEdgeColor', [1    0.5    0], 'MarkerFaceColor', [1    0.5    0], 'MarkerFaceAlpha', 0.6)
    ylabel('norm. EPSP amplitude')
    xlabel('Stim (pre)- CS (post) delay (s)')
    box off; axis square;













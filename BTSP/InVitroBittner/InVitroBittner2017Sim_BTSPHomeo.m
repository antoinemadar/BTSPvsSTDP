%simulation of in vitro demonstration of BTSP rule as in Bittner et al. 2017
clear 
close all
%% Parameters

% Synapses params
Params.Idecay = 10e-3; % EPSC time constant, in sec
Wsyn = 77e-12 % 63.75e-12; % synaptic weight of the stimulated input fiber (or ensemble of fibers), in Amps. This stimulation-evoked current was adjusted to get 2mV amplitude for initial EPSP 
Imax = 85e-12; % max synaptic weight, in Amps (used in Sheffield or Milstein's style simulations). Not relevant in this experiment, but Pamp can be expressed in terms of it, and the parameter is required in BTSPsimple_LIFadapt_ODEs in case weight capping is implemented
N = 100; % number of inputs
Wsd = 10; %1.6*PFsd/(L/N); %standard deviation of initial synaptic weight vector (for gaussian connectivity)

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
Params.Pamp = 4.2e-12; 
Params.PostPreBoost = 1.1; % amplitude of weight change for postpre rule which is triggered on CS, not on spikes (so, to be somewhat balanced, it should = 10 so that amplitude is the Bittner2017 amplitude/1CS and not /10spikes like the prepost rule)
Params.CSproba = [];
Params.capWeights = 0; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 
Params.HomeoNorm = 1; % 0 if no homeostatic normalization rule, 1 if yes.
Params.WUdyn = 0; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Params.Wrest = 0; % Synaptic weight resting state, in Amps (towards which any synapse would relax if Wtarget = 0) 
Params.tau_wu = 5; % weight update time constant, in sec
Params.WUdyn = 0; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Params.Wrest = 0; % Synaptic weight resting state, in Amps (towards which any synapse would relax if Wtarget = 0) 
Params.tau_wu = 5; % weight update time constant, in sec

% Protocol params
pairings = 5; %number of pairings
NumSpikes = 10; % number of spikes
StimFR = 20; % Hz
stimD = 15; % Duration of stim protocol (including input spikes and CS induction), in sec: up to 4 sec before or 4 after 
ITI = 15; % Interval between stim protocols, in sec
dt = 0.001; % time resolution in sec
Delay = [-4000, -3000, -2000, -1000, -750, -500, -250, -100, -10, -1, 1, 10, 100, 250, 500, 750, 1000, 2000, 3000, 4000]; % delay between last stim and CS, in ms. Negative values correspond to post-before-pre condition.

%% Set Up

StimLen = stimD/dt;
Tsweep = dt:dt:stimD+ITI; % ms resolution
Trun = dt:dt:pairings*(stimD+ITI);
TswLen = length(Tsweep);

NeuronID = 1:N;
W = zeros(N,length(Trun));
Wmax = Imax;
W(:,1) = Wmax*exp(-0.5*(NeuronID - N/2).^2/Wsd^2); % vector of N rows following a gaussian distribution of Wmax amplitude 
[val, Nbittner] = min(abs(W(:,1)-Wsyn));

% stimulated input neuron
Nidx = Nbittner;

%% Protocol 

% CS trace
CS = zeros(1,StimLen); 
CSidx = round(StimLen/2);
CS(CSidx) = 1;
CS = [CS, zeros(1,ITI/dt)];
CSall = repmat(CS, 1, pairings);

InRasters = {};
I = {}; V = {}; V2 = {}; 
OutRaster= {}; 
D = {}; P = {}; Pact = {}; Dact = {}; SRA = {}; Wraw = {}; W2={}; Wtarget = {}; CSbis = {};
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
        InRasters{d} = zeros(N,length(Trun));
        InRasters{d}(Nidx,:) = repmat(InRaster, 1, pairings);
        
        % add 1 stim spike at beginning and before the end to measure pre and post-induction EPSP amplitude
        InRasters{d}(Nidx, 1) = 1;
        InRasters{d}(Nidx, end-1/dt) = 1;
        clear InRaster
        
        % run LIF + plasticity model
        [I{d}, V{d}, V2{d}, CSbis{d}, OutRaster{d}, D{d}, P{d}, Dact{d}, Pact{d}, SRA{d}, Wraw{d}, W2{d}, Wtarget{d}]  = BTSPplus_LIFadapt_ODEs(Params, W, InRasters{d}, CSall) ;
        
        % get post-induction EPSP amplitude and normalize to pre-induction amplitude (i.e. Vsyn)
        EPSPpost(d) = max(V2{d}(end-1/dt:end)); %in Volts
        EPSPpre(d) = max(V2{d}(1:100)); % in Volts
        EPSPampl(d).pre = 10^3*(EPSPpre(d)-Params.Vrest); %in mv
        EPSPampl(d).post = 10^3*(EPSPpost(d)-Params.Vrest); %in mV
        normEPSPampl(d) = EPSPampl(d).post/EPSPampl(d).pre; % amplitude of post-induction EPSP normalized to pre-induction EPSP
    
    end

%% figure
Didx = 10;

figure % Vm on first and last pairing
CStimes = Trun(logical(CSall));
TrunCScentered = Trun - Trun(CSidx); 
TrunCScentered2 = Trun - CStimes(end); %time centered on last CS
    plot(TrunCScentered, V{Didx}.*10^3, 'k'); hold on
    plot(TrunCScentered2, V2{Didx}.*10^3, 'r'); hold on
    xlim([-1 1])
    xlabel('time (s)'); ylabel('membrane potential (mV)')
    title(['1st (black) and last (red) pairing: delay = ' num2str(Delay(Didx)) ' ms'])
    box off;

% for given Nidx, plot normalized post-induction EPSP as a function of delay (with the rule overlaid)

delay1 = -5:dt:0; % in sec
delay2 = 0:dt:5; % in sec
Wd = (max(normEPSPampl(Delay<0))-1)*exp(delay1/Params.Ddecay)+1;
Wp = (max(normEPSPampl(Delay>0))-1)*exp(-delay2/Params.Pdecay)+1;

figure
plot(delay1, Wd, 'b--' ); hold on
plot(delay2, Wp, 'r--' ); hold on
scatter(Delay.*dt, normEPSPampl, 'ok'); hold on
ylabel('norm. EPSP amplitude')
xlabel('CS (post) - Stim (pre) delay (s)')
% ylim([1 max(normEPSPampl)+0.5])

% 
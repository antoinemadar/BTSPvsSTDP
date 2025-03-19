function [I, V, V2, CS, OutRaster, D, P, Dact, Pact, SRA, W, W2, Wtarget]  = BTSPplus_LIFadapt_ODEs(params, W, InRasters, CS) 
% 
% [I, V, V2, CS, OutRaster, D, P, Dact, Pact, SRA, W, W2, Wtarget]  = BTSPplus_LIFadapt_ODEs(params, W, InRasters, CS) 
% 
% W and InRasters are N x M matrices, with N = number of input neurons and M = number of time points
% CS is a binary row vector of length M, taking value 1 when a CS occurs

%% params

dt = params.dt;

% synaptic params
Idecay = params.Idecay;
Imax = params.Imax;

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
Pdecay = params.Pdecay; % Pre-before-Post time constant, in sec (1.31s on average in Bittner et al. 2017)
Ddecay = params.Ddecay; % Post-before-Pre time constant, in sec (0.69s on average in Bittner et al. 2017)
Pamp = params.Pamp; % *WARNING: not like in BTSPplus_LIFadapt* the provided value needs to already be a current in Amps
PostPreBoost = params.PostPreBoost; % amplitude of weight change for postpre rule which is triggered on CS, not on spikes

capWeights = params.capWeights; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 

HomeoNorm = params.HomeoNorm; % 0 if no homeostatic normalization rule, 1 if yes.

WUdyn = params.WUdyn; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Wrest = params.Wrest;
tau_wu = params.tau_wu; % weight update time constant, in sec

%% input-output dynamics using Euler's forward method to solve ODEs

%initial values
I = zeros(1,size(W,2)); 
V = I;
D = I;
SRA = I;
OutRaster = I;
P = zeros(size(W));
Pact = P;
Dact = P;
V(1) = Vrest;
V2 = V;
W2 = W;
Wtarget = W;

for i = 1:size(W,2)-1

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
           if ~isempty(CSproba)
           CS(i) = binornd(1, CSproba);
           end
    end
    V2(i+1) = V(i+1); % to keep a version of V without spikes

    % BTSP variables P and D + update Synaptic weights W
    D(i+1) = D(i) - dt*D(i)/Ddecay + CS(i); % output train of Complex Spikes convolved with Post-before-Pre update rule, in %
    for n = 1:size(W,1) % for each presynaptic neuron
        P(n,i+1) = P(n,i) - dt*P(n,i)/Pdecay + InRasters(n,i); % input spiketrains convolved with Pre-before-Post update rule, in %
        Pact(n,i+1) = CS(i)*P(n,i)*Pamp; %evaluate pre-before-post variable at time of CS
        Dact(n,i+1) = InRasters(n,i)*D(i)*Pamp*PostPreBoost; %evaluate post-before-pre variable at time of an input spike
        
        % update Synaptic weights W, with or without dynamics
        if WUdyn == 1 && HomeoNorm == 0 % dynamic update with time constant tau_wu
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
        if WUdyn == 0
            W2(:,i+1) = W(:,i+1).*sum(W(:,1))./sum(W(:,i+1)); % Homeostatic synaptic normalization, multiplicative method
        else % WUdyn == 1
            Wtarget(:,i+1) = W(:,i+1).*sum(W(:,1))./sum(W(:,i+1));
            W2(n,i+1) = W(n,i) + dt .* (Wrest - W(n,i) + Wtarget(n,i+1))./ tau_wu;
        end
    end
end


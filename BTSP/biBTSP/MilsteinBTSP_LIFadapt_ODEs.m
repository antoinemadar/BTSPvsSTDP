function [I, V,V2, CS, OutRaster, D, P, SRA, W, W2]  = MilsteinBTSP_LIFadapt_ODEs(params, W, InRasters, CS, laps_end) 

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

% capWeights = params.capWeights; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 

%% input-output dynamics using Euler's forward method to solve ODEs
lapL = find(laps_end,1); % number of indices in a lap

%initial values
I = zeros(1,size(W,2)); 
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
           if ~isempty(CSproba)
           CS(i) = binornd(1, CSproba);
               if CS(i) == 1
               CS(i:i+299) = 1;
               end
           end
    end
    V2(i+1) = V(i+1); % to keep a version of V without spikes

    % BTSP variables P and D + update Synaptic weights W
    D(i+1) = D(i) - dt*D(i)/CSdecay + CS(i); % Global Instructive Signal from the CS

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

%% local function
    function q = q_sigmo(x,a,b) % sigmoid defined on [0,1], being 0 at x=0 and 1 at x=1.
s = 1./(1+exp(-b.*(x-a)));
s0 = 1./(1+exp(-b.*(0-a)));
s1 = 1./(1+exp(-b.*(1-a)));
q = (s-s0)./(s1-s0);
    end

end
%simulation of in vitro demonstration of BTSP rule as in Bittner et al. 2017
clear 
close all
%% Parameters

% Synapses params
Idecay = 10e-3; % EPSC time constant, in sec
% Wsyn = 100e-12; % synaptic weight of the stimulated input fiber (or ensemble of fibers), in Amps. Based on WangBi2005

% Plasticity params
Pdecay = 16.8e-3; % Pre-before-Post time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Ddecay = 33.7e-3; % Post-before-Pre time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Pdecay2 = 946e-3; % tau_x from All-to-All model in Table 4 of Pfistergerstner2006
Ddecay2 = 27e-3; % tau_y

Pamp = 0.61/100; % A2+ in PfisterGerstner2006. Peak Pre-before-Post weight change for paired protocol (0.5% in SongAbbot2000, 0.6pA with no Imax in MehtaWilson2000). 
Damp = 0.16/100; % A2- in PfisterGerstner2006. Peak Post-before-Pre weight change for paired protocol (-0.525% in SongAbbott2000, -90% of Pamp in MehtaWilson2000, which didn't have a maximum weight)
Pamp2 = 0.14/100; % A3- in PfisterGerstner2006
Damp2 = 0.67/100; % A3+

% Protocol params
pairings = 60; %number of pairings
delta = 10e-3; % delta time (in sec) between the input spike and output spike of a triplet protocol. Assumes a symetric triplet protocol
ITI = 1/0.05; % Interval between stim protocols, in sec (same as BiPoo2005)
dt = 1e-3; % time resolution in sec

%% Design Protocol
InSweep = zeros(2,11);
OutSweep = zeros(2,11);

% pre-post-pre
InSweep(1,1) = 1;
OutSweep(1,1+delta/dt) = 1;
InSweep(1,1+2*delta/dt) = 1;

% post-pre-post
OutSweep(2,1) = 1;
InSweep(2,1+delta/dt) = 1;
OutSweep(2,1+2*delta/dt) = 1;

% put rasters together
spacer = zeros(2,ITI/dt);
SweepIn = [InSweep spacer];
SweepOut = [OutSweep spacer];

InRaster = repmat(SweepIn, 1, pairings);
OutRaster = repmat(SweepOut, 1, pairings);

% add 1 input spike at beginning and before the end to measure pre and post-induction EPSP amplitude
InRaster = [spacer InRaster spacer];
% InRaster(:,1) = 1;
% InRaster(:,end-1/dt) = 1;

OutRaster = [spacer OutRaster spacer];

% time vector
Trun = dt:dt:size(InRaster,2)*dt;

%% Run Model

% W = Wsyn*ones(2,length(Trun));
W = ones(2,length(Trun)); % weights in arbitrary unit
Wstdp = W;
I = zeros(2,length(Trun)); % currents in arbitrary unit, same as weight
D = I; D2 = I; P = I; P2 = I; Pact = I; Dact = I;

for n = 1:size(InRaster,1) % for each protocol, pre-post-pre or post-pre-post
for i = 1:size(W,2)-1

I(n,i+1)  = I(n,i) - dt*I(n,i)./Idecay + InRaster(n,i).*W(n,i);

D(n,i+1) = D(n,i) - dt*D(n,i)/Ddecay + OutRaster(n,i); % output spiketrain convolved with Post-before-Pre update rule, in %
D2(n,i+1) = D2(n,i) - dt*D2(n,i)/Ddecay2 + OutRaster(n,i); % postsynaptic triplet variable
P(n,i+1) = P(n,i) - dt*P(n,i)/Pdecay + InRaster(n,i); % input spiketrains convolved with Pre-before-Post update rule, in %
P2(n,i+1) = P2(n,i) - dt*P2(n,i)/Pdecay2 + InRaster(n,i); % presynaptic triplet variable
Pact(n,i+1) = OutRaster(n,i)*P(n,i+1)*(Pamp + Damp2*D2(n,i)); % total potentiation term
Dact(n,i+1) = InRaster(n,i)*D(n,i+1)*(Damp + Pamp2*P2(n,i)); % total depression term

W(n,i+1) = W(n,i) + Pact(n,i+1) - Dact(n,i+1); % weight update
Wstdp(n,i+1) = Wstdp(n,i) + OutRaster(n,i)*P(n,i+1)*Pamp - InRaster(n,i)*D(n,i+1)*Damp; % weight update without the triplet terms. i.e. pure STDP

% for i = 2:size(W,2)-1
% 
% I(n,i+1)  = I(n,i) - dt*I(n,i)./Idecay + InRaster(n,i).*W(n,i);
% 
% D(n,i+1) = D(n,i) - dt*D(n,i)/Ddecay + OutRaster(n,i); % output spiketrain convolved with Post-before-Pre update rule, in %
% D2(n,i+1) = D2(n,i) - dt*D2(n,i)/Ddecay2 + OutRaster(n,i); % postsynaptic triplet variable
% P(n,i+1) = P(n,i) - dt*P(n,i)/Pdecay + InRaster(n,i); % input spiketrains convolved with Pre-before-Post update rule, in %
% P2(n,i+1) = P2(n,i) - dt*P2(n,i)/Pdecay2 + InRaster(n,i); % presynaptic triplet variable
% Pact(n,i) = OutRaster(n,i)*P(n,i)*(Pamp + Damp2*D2(n,i-1)); % total potentiation term
% Dact(n,i) = InRaster(n,i)*D(n,i)*(Damp + Pamp2*P2(n,i-1)); % total depression term
% 
% W(n,i+1) = W(n,i) + Pact(n,i) - Dact(n,i); % weight update
end

STDPratio(n) = W(n,end)/W(n,1);
STDPdelta(n) = W(n,end)-W(n,1);

figure % pre-post-pre
plot(Trun,W(n,:),'r-'); hold on
plot(Trun,Wstdp(n,:),'k-')
xlabel('time (sec)')
ylabel('synaptic weight')
xlim([0 Trun(end)]);
ylim([1 1.25]);
legend('Triplet rule', 'Doublet rule', 'Location', 'northwest')
axis square
box off
end
        
% % get post-induction EPSC amplitude and normalize to pre-induction amplitude
% EPSCpost(n) = max(I(n,end-1/dt:end));
% EPSCpre(n) = max(I(n,1:100)); 
% STDPratio(n) = EPSCpost/EPSCpre; % amplitude of post-induction EPSC normalized to pre-induction EPSP
% STDPdelta(n) = EPSCpost-EPSCpre;    

%% figure

figure % 1st pre-post-pre
subplot(3,1,1) % input raster and presynaptic STDP variables
    plot(repmat(Trun(InRaster(1,:)>0),2,1), repmat([0; 1], 1, length(Trun(InRaster(1,:)>0))), 'k', 'LineWidth', 1.5); hold on
    plot(Trun, P(1,:), 'b'); hold on
    plot(Trun, P2(1,:), 'r'); hold on
    % legend('input spikes', 'presynaptic doublet variable', 'presynaptic triplet variable', 'Location', 'bestoutside')
    title('presynaptic spikes and variables')
    xlim([ITI-0.01 ITI+0.15]);
    set(gca, 'XTickLabel', 'off', 'XTick', [])
    box off
subplot(3,1,2) % outputput raster and postsynaptic STDP variables
    plot(repmat(Trun(OutRaster(1,:)>0),2,1), repmat([0; 1], 1, length(Trun(OutRaster(1,:)>0))), 'k', 'LineWidth', 1.5); hold on
    plot(Trun, D(1,:), 'b'); hold on
    plot(Trun, D2(1,:), 'r'); hold on
    % legend('output spikes', 'postsynaptic doublet variable', 'postsynaptic triplet variable', 'Location', 'bestoutside')
    xlim([ITI-0.01 ITI+0.15]);
    title('postsynaptic spikes and variables')
    set(gca, 'XTickLabel', 'off', 'XTick', [])
    box off
subplot(3,1,3) % weights
    plot(Trun,W(1,:),'r-'); hold on
    plot(Trun,Wstdp(1,:),'k--')
    % legend('Triplet rule', 'Doublet rule', 'presynaptic triplet variable', 'Location', 'bestoutside')
    xlabel('time (sec)')
    title('synaptic weight')
    xlim([ITI-0.01 ITI+0.15]);
    legend('Triplet rule', 'Doublet rule', 'Location', 'northeast')
    box off

figure % 1st post-pre-post
subplot(3,1,1) % input raster and presynaptic STDP variables
    plot(repmat(Trun(InRaster(2,:)>0),2,1), repmat([0; 1], 1, length(Trun(InRaster(2,:)>0))), 'k', 'LineWidth', 1.5); hold on
    plot(Trun, P(2,:), 'b'); hold on
    plot(Trun, P2(2,:), 'r'); hold on
    title('presynaptic spikes and variables')
    % legend('input spikes', 'presynaptic doublet variable', 'presynaptic triplet variable', 'Location', 'bestoutside')
    xlim([ITI-0.01 ITI+0.15]);
    box off
    % axis square
subplot(3,1,2) % input raster and presynaptic STDP variables
    plot(repmat(Trun(OutRaster(2,:)>0),2,1), repmat([0; 1], 1, length(Trun(OutRaster(2,:)>0))), 'k', 'LineWidth', 1.5); hold on
    plot(Trun, D(2,:), 'b'); hold on
    plot(Trun, D2(2,:), 'r'); hold on
    title('postsynaptic spikes and variables')
    % legend('output spikes', 'postsynaptic doublet variable', 'postsynaptic triplet variable', 'Location', 'bestoutside')
    xlim([ITI-0.01 ITI+0.15]);
    box off
    % axis square
subplot(3,1,3) % weights
    plot(Trun,W(2,:),'r-'); hold on
    plot(Trun,Wstdp(2,:),'k--')
    % legend('Triplet rule', 'Doublet rule', 'presynaptic triplet variable', 'Location', 'bestoutside')
    xlabel('time (sec)')
    title('synaptic weight')
    xlim([ITI-0.01 ITI+0.15]);
    legend('Triplet rule', 'Doublet rule', 'Location', 'southeast')
    box off
    % axis square

% figure % 2nd pre-post-pre
% subplot(3,1,1) % input raster and presynaptic STDP variables
%     plot(repmat(Trun(InRaster(1,:)>0),2,1), repmat([0; 1], 1, length(Trun(InRaster(1,:)>0))), 'k'); hold on
%     plot(Trun, P(1,:), 'b'); hold on
%     plot(Trun, P2(1,:), 'r'); hold on
%     xlim([2*ITI-0.01 2*ITI+0.15]);
%     box off
% subplot(3,1,2) % input raster and presynaptic STDP variables
%     plot(repmat(Trun(OutRaster(1,:)>0),2,1), repmat([0; 1], 1, length(Trun(OutRaster(1,:)>0))), 'k'); hold on
%     plot(Trun, D(1,:), 'b'); hold on
%     plot(Trun, D2(1,:), 'r'); hold on
%     xlim([2*ITI-0.01 2*ITI+0.15]);
%     box off
% subplot(3,1,3) % weights
%     plot(Trun,W(1,:),'k-'); hold on
%     plot(Trun,Wstdp(1,:),'r--')
%     xlabel('time (sec)')
%     ylabel('synaptic weight')
%     xlim([2*ITI-0.01 2*ITI+0.15]);
%     box off
% 
% figure % 
% subplot(3,1,1) % input raster and presynaptic STDP variables
%     plot(repmat(Trun(InRaster(1,:)>0),2,1), repmat([0; 1], 1, length(Trun(InRaster(1,:)>0))), 'k'); hold on
%     plot(Trun, P(1,:), 'b'); hold on
%     plot(Trun, P2(1,:), 'r'); hold on
%     xlim([0 2*ITI+1]);
%     box off
% subplot(3,1,2) % input raster and presynaptic STDP variables
%     plot(repmat(Trun(OutRaster(1,:)>0),2,1), repmat([0; 1], 1, length(Trun(OutRaster(1,:)>0))), 'k'); hold on
%     plot(Trun, D(1,:), 'b'); hold on
%     plot(Trun, D2(1,:), 'r'); hold on
%     xlim([0 2*ITI+1]);
%     box off
% subplot(3,1,3) % weights
%     plot(Trun,W(1,:),'k-')
%     xlabel('time (sec)')
%     ylabel('synaptic weight')
%     xlim([0 2*ITI+1]);
%     box off
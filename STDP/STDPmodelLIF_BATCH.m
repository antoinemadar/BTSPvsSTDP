%STDPmodelLIF_BATCH
clear 
close all
tic
%% Fixed Parameters
%Params for inputs
Params.N = 50; % number of input neurons
Params.L = 300; % length of track in cm
Params.PFsd = 18; % Standard deviation of mean PF, in cm
Params.PFamp = 10; %peak FR in Hz
Params.Nlaps = 30; % number of laps
Params.period = 20; % lap duration, in sec
Params.dt = 0.001; % time resolution in sec

% Synapses params
Params.Imax = 145e-12; % max synaptic weight, a current in Amps. % NolanWiles2011 uses 10-226pA. Mehta uses 1.5-4nA. SongAbbot use conductances with G between 0 and 150pS, but not sure about the initial weight matrix. 
Params.Idecay = 10e-3; % EPSC time constant, in sec
Params.Wsd = 1.6*Params.PFsd/(Params.L/Params.N); %standard deviation of initial synaptic weight vector (for gaussian connectivity)
Params.maxW = 1; % from 0 to 1, proportion of Imax defining the initial maximum weight

%Integrate and Fire params for output neuron
Params.Rm = 100e6; % membrane resistance, in Ohm
Params.tau_m = 20e-3 ; %membrane time constant, in sec
Params.Vrest = -70e-3 ; % resting membrane potential, in Volt
Params.Vthr = -54e-3; % spike threshold, in Volt
Params.Vreset = -60e-3; % reset potential after spike, in Volt

% output PF 
Params.Nbin = 50; % Number of bins in which the length of the track is divided

%% STDP params to vary
decay = 10:10:100; %ms
Amp = 0.5:0.5:5; %in %

delay1 = -500:1:0; % ms
delay2 = 0:1:500; %ms
for n = 1:length(decay)
    for i = 1:length(Amp)
        
    %STDP rules    
    Params.Pdecay = decay(n)/1000; % Pre-before-Post time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
    Params.Ddecay = decay(n)/1000; % Post-before-Pre time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
    Params.Pamp = Amp(i)/100; % peak Pre-before-Post weight change, in percent of Imax (0.5% in SongAbbot2000, 0.6% in MehtaWilson2000)
    Params.Damp = -Amp(i)/100; % peak Post-before-Pre weight change, in percent of Imax (0.525% in SongAbbott2000, 90% in MehtaWilson2000, which didn't have a maximum weight)
    
    Wd{n,i} = -Amp(i)*exp(delay1/decay(n));
    Wp{n,i} = Amp(i)*exp(-delay2/decay(n));

    % model
    output(n,i) = STDP_LIF(Params);
    
    % results
    COMslope(n,i) = output(n,i).COMslope;
    shiftPval(n,i) = output(n,i).shiftPval;
    shiftR2(n,i) = output(n,i).shiftR2;
    SD_meanPF(n,i) = output(n,i).SD_meanPF;
    PFpeakFR(n,i) = max(output(n,i).meanFRmap);

    end
end

[row_p,col_p] = find(shiftPval<0.05);
[row_p_sd,col_p_sd] = find(shiftPval<0.05 & SD_meanPF<20);

%% figures

figure(1) %STDP rules
subplot(1,2,1)
    for n = 1:length(decay)
    plot([delay1 delay2], [Wd{n,1} Wp{n,1}] );hold on
    end
    legend(split(num2str(decay)), 'Location', 'BestOutside')
    xline(0); hold on
    yline(0); 
    xlabel('delay (ms)'); ylabel('Weight change (% of EPSC_{max})');
    title('max Weight change = 0.5%')
    box off; axis square;
subplot(1,2,2)
    for i = 1:length(Amp)
    plot([delay1 delay2], [Wd{2,i} Wp{2,i}] );hold on
    end
    legend(split(num2str(Amp)), 'Location', 'BestOutside')
    xline(0); hold on
    yline(0); 
    xlim([-100 100]);
    xlabel('delay (ms)'); ylabel('Weight change (% of EPSC_{max})');
    title('Tau_{STDP} = 20 ms')
    box off; axis square;

figure(2)
subplot(1,2,1)
    plot(decay, COMslope(:,1), 'k-');
    xlabel('Tau_{STDP} (ms)');
    ylabel('COM shifting slope (cm/lap)');
    title('max Weight change = 0.5%')
    box off; axis square;
subplot(1,2,2)
    plot(Amp, COMslope(2,:), 'k-');
    xlabel('max Weight change (%)');
    ylabel('COM shifting slope (cm/lap)');
    title('Tau_{STDP} = 20 ms')
    box off; axis square;

figure(3)
imagesc(Amp,decay, COMslope); hold on
scatter(Amp(col_p), decay(row_p), 'k*'); hold on
scatter(Amp(col_p_sd), decay(row_p_sd), 'ko');
xlabel('max Weight change (%)'); ylabel('Tau_{STDP} (ms)');
colormap(brewermap(256,'*RdBu')); 
caxis([-1.5 1.5]);
c = colorbar; c.Label.String = 'COM linear regression: shifting slope (cm/lap)';
box off; axis square;
title('sample size per square = 1')

figure(4)
% for n = 1:length(decay)
%     for i = 1:length(Amp)
%         SD_meanPF(n,i) = output(n,i).SD_meanPF;
%     end
% end
imagesc(Amp,decay, SD_meanPF); hold on
xlabel('max Weight change (%)'); ylabel('Tau_{STDP} (ms)');
colormap(brewermap(256,'oranges'));
c = colorbar; c.Label.String = 'SD of mean PF (cm)';
box off; axis square;
title('sample size per square = 1')

figure(5)
xval = split(num2str(Amp));
yval = split(num2str(decay));
% for n = 1:length(decay)
%     for i = 1:length(Amp)
%         PFpeakFR(n,i) = max(output(n,i).meanFRmap);
%     end
% end
imagesc(Amp, decay, PFpeakFR); hold on
xlabel('max Weight change (% of Imax)'); ylabel('Tau_{STDP} (ms)');
% set(gca, 'XTick', 1:length(Amp), 'XTickLabel', xval, 'YTick', 1:length(decay), 'YTickLabel', yval);
colormap(flipud(gray(256)));
c = colorbar; c.Label.String = 'mean PF peak FR (Hz)';
box off; axis square;
title('sample size per square = 1')

figure(6)
% for n = 1:length(decay)
%     for i = 1:length(Amp)
%         shiftR2(n,i) = output(n,i).shiftR2;
%     end
% end
imagesc(Amp, decay, shiftR2); hold on
xlabel('max Weight change (% of Imax)'); ylabel('Tau_{STDP} (ms)');
% set(gca, 'XTick', 1:length(Amp), 'XTickLabel', xval, 'YTick', 1:length(decay), 'YTickLabel', yval);
colormap(brewermap(256,'purples'));
caxis([0 1]);
c = colorbar; c.Label.String = 'COM linear regression: R-squared';
box off; axis square;
title('sample size per square = 1')

figure
% for n = 1:length(decay)
%     for i = 1:length(Amp)
%         shiftPval(n,i) = output(n,i).shiftPval;
%     end
% end
imagesc(Amp, decay, shiftPval); hold on
xlabel('max Weight change (% of Imax)'); ylabel('Tau_{STDP} (ms)');
% set(gca, 'XTick', 1:length(Amp), 'XTickLabel', xval, 'YTick', 1:length(decay), 'YTickLabel', yval);
colormap(brewermap(256,'*greens'));
caxis([0 0.1]);
c = colorbar; c.Label.String = 'COM linear regression: p-value';
box off; axis square;
title('sample size per square = 1')

toc
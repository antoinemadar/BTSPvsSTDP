%STDPmodelLIF_BATCH_ParamSpaceRepeat
clear 
close all

% for dynamic inputs: slope distribution
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA3_SlopeDistrib.mat')
load('G:\My Drive\LabReasearch\Projects\Sheffield\Can1-Shifting\CanData\CA1_SlopeDistrib.mat')

tic
%% Parameters
% Inputs params
Params.N = 100; % number of input neurons
Params.L = 300; % length of track in cm
Params.PFsd = 18; % Standard deviation of mean PF, in cm
Params.PFamp = 10; %peak FR in Hz
Params.Nlaps = 30; % number of laps
Params.period = 6; % lap duration, in sec
Params.dt = 0.001; % time resolution in sec
Params.InShift = 0; % 0 if static input, 1 if shift
Params.PslopeBinCenters = [-1.6 -1.4, -1.2 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6]; %in cm/lap
Context = 'N'; % N for novel, F for Familiar, NandF for all

if ismember('F', Context)
    Params.Pslope = CA3_F_SlopeHist.Values;
elseif ismember('N', Context)
    Params.Pslope = CA3_N_SlopeHist.Values;
elseif ismember('FandN', Context)
    Params.Pslope = CA3_FandN_SlopeHist.Values;
end

% Synapses params
Params.Imax = 85e-12; % max synaptic weight, a current in Amps. % NolanWiles2011 uses 10-226pA. Mehta uses 1.5-4nA. SongAbbot use conductances with G between 0 and 150pS, but not sure about the initial weight matrix. 
Params.Idecay = 10e-3; % EPSC time constant, in sec
Params.Wsd = 10; %1.6*Params.PFsd/(Params.L/Params.N); %in cm, standard deviation of initial synaptic weight vector (for gaussian connectivity)
Params.maxW = 1; % from 0 to 1, proportion of Imax defining the initial maximum weight

% Output neuron LIF params
Params.Rm = 100e6; % membrane resistance, in Ohm
Params.tau_m = 20e-3 ; %membrane time constant, in sec
Params.Vrest = -70e-3 ; % resting membrane potential, in Volt
Params.Vthr = -54e-3; % spike threshold, in Volt
Params.Vreset = -60e-3; % reset potential after spike, in Volt

Params.Adapt = 0; % 1 or 0, to implement spike-rate adaptation or not
Params.Ek = -70e-3; % K+ equilibrium potential, in Volt
Params.dSRA = 0.06; % increment value of spike-rate adaptation variable modeling adapting K+ leak current (cf DayanAbbott fig 5.6 and eq 5.13)
Params.tau_sra = 100e-3; % time constant for exponential decay of SRA variable

Params.Shunting = 0; % implement shunting inhib (capping I) if = 1, 0 otherwise. 
Params.Icap = 350e-12; % in Amps

% Plasticity params
Params.Pdecay = 20e-3; % Pre-before-Post time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Params.Ddecay = 20e-3; % Post-before-Pre time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
Params.Pamp = 0.5/100; % in terms of Imax. Peak Pre-before-Post weight change (0.5% in SongAbbot2000, 0.6pA with no Imax in MehtaWilson2000)
Params.Damp = -0.5/100; % in terms of Imax. Peak Post-before-Pre weight change (-0.525% in SongAbbott2000, -90% of Pamp in MehtaWilson2000, which didn't have a maximum weight)

Params.capWeights = 1; %cap weights at Imax if = 1. Runnaway weights are allowed if set to 0. 

Params.HomeoNorm = 0; % 0 if no homeostatic normalization rule, 1 if yes.

Params.WUdyn = 0; %  Weight Update dynamics: 0 if instantaneous update, 1 if dynamic update
Params.Wrest = 0; % Synaptic weight resting state, in Amps (towards which any synapse would relax if Wtarget = 0) 
Params.tau_wu = 5; % weight update time constant, in sec

% output PF 
Params.Nbin = 50; % Number of bins in which the length of the track is divided

% number of simulations
Nsim = 20;

%% STDP params to vary
STDPrule = 0 ; % 0 or 1. 1 if varying plasticity params, 0 if varying input params controling FRout range

if STDPrule == 0
% I choose not to vary Imax, because I want to distinguish between the effect of output FR and max weight change, which would be confounded here. 
% I designed pairs of PFamp and PFsd values to cover a realistic range of meanFRin, meanFRout and PeakFRout
PFsd = [10 10 18 18 18 26 26 26] ; % Standard deviation of mean input PF, in cm (cf range in Can CA3 data + YuShouval2006). Lee2015 suggests it increases along proximodistal axis
PFamp = [15 20 10 15 20 10 12 15]; % in Hz (cf Mou et al. 2018, LeeKnierim2004b, and Lee Knierim 2015). Lee2004b and 2015 suggest mean is around 9Hz in CA3, range [2 to 50Hz]
FRin = [1.2, 1.7, 1.5, 2.2, 3, 2.2, 2.6, 3.3]; % approx. mean FR of each input neurons for the corresponding PFsd and PFamp stated above. This is not a parameter, but used for plots
Wsd = [7 10 13]; % SD of the gaussian connectivity vector (unit = number of input neurons)
param1 = FRin; 
param2 = Wsd;
else
decay = [10 20 30 50 100]; %[10 20 30] % in ms
Amp = [0.5 1 2 4 10]; %[1 5 10 15 30] %in %
param1 = decay;
param2 = Amp;
end
delay1 = -500:1:0; % ms
delay2 = 0:1:500; %ms

%% Simulations
dummy = 0;
    for p1 = 1:length(param1)
    for p2 = 1:length(param2)
        
        dummy = dummy + 1;
        modelgp(p1,p2) = dummy;

%         modelname{p1,p2} = strcat( num2str(PFsd(p1)), 'cm, ', num2str(PFamp(p1)), 'Hz, ', num2str(param2(p2)), ' i.n.' );
%         for n = 1:Nsim
%                     % params matrices
%             param1Mat(p1,p2,n) = param1(p1);
%             param2Mat(p1,p2,n) = param2(p2);
%             modelgpMat(p1,p2,n) = modelgp(p1,p2);
%             modelnameMat{p1,p2,n} = modelname(p1,p2);
%            
%         end
%     end
%     end

        if STDPrule == 1
        modelname{p1,p2} = strcat(num2str(param1(p1)), ' ms - ', num2str(param2(p2)), ' %');
        Wd{p1,p2} = -Amp(p2)*exp(delay1/decay(p1));
        Wp{p1,p2} = Amp(p2)*exp(-delay2/decay(p1));
        else
        modelname{p1,p2} = strcat( num2str(PFsd(p1)), 'cm, ', num2str(PFamp(p1)), 'Hz, ', num2str(param2(p2)), ' i.n.' );
%         modelname{p1,p2} = strcat('FRin ', num2str(param1(p1)), ' Hz, ', 'Wsd ', num2str(param2(p2)));
        end

        for n = 1:Nsim

            if STDPrule == 0 % params controlling FRout
                Params.PFamp = PFamp(p1); %peak FR in Hz
                Params.PFsd = PFsd(p1); % Standard deviation of mean PF, in cmC
                Params.Wsd = Wsd(p2); % 
            else %STDP rules
                Params.Pdecay = decay(p1)/1000; % Pre-before-Post time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
                Params.Ddecay = decay(p1)/1000; % Post-before-Pre time constant, in sec (20ms in SongAbbott2000 and YuShouval2006, 10ms in MehtaWilson2000)
                Params.Pamp = Amp(p2)/100; % peak Pre-before-Post weight change, in percent of Imax (0.5% in SongAbbot2000, 0.6% in MehtaWilson2000)
                Params.Damp = -Amp(p2)/100; % peak Post-before-Pre weight change, in percent of Imax (0.525% in SongAbbott2000, 90% in MehtaWilson2000, which didn't have a maximum weight)
                
            end
    
            % model
            output = STDPplus_LIFadapt(Params);
            
            % package results
            COMslope(p1, p2, n) = output.COMslope;
            shiftPval(p1, p2, n) = output.shiftPval;
            shiftR2(p1, p2, n) = output.shiftR2;
            PFpeak(p1,p2,n) = max(output.meanFRmap);
            PFpeakMou(p1,p2,n) = max(output.meanFRmapMou);
            FRout(p1, p2, n) = output.meanFRout;

            SD_meanPF(p1, p2, n) = output.SD_meanPF;
%             PFsdDiff(p1, p2, n) = output.PFsdOut(end) - output(n).PFsdOut(1);
            PFsdDiff(p1, p2, n) = mean(output.PFsdOut(Params.Nlaps-2:end)) - mean(output.PFsdOut(1:3));
%             SpatialFRout{p1, p2, n} = output.PF;
%             COMbin{p1, p2, n} = output.COMbin;

            clear output

            % params matrices
            param1Mat(p1,p2,n) = param1(p1);
            param2Mat(p1,p2,n) = param2(p2);
            modelgpMat(p1,p2,n) = modelgp(p1,p2);
            modelnameMat{p1,p2,n} = modelname(p1,p2);
           
        end
        
        meanSlope(p1, p2) = mean(COMslope(p1, p2, :));
        meanAbsSlope(p1, p2) = mean(abs(COMslope(p1, p2, :)));
        meanFRout(p1, p2) = mean(FRout(p1,p2,:));
        meanPFpeak(p1,p2) = mean(PFpeak(p1, p2, :));
        meanPFpeakMou(p1,p2) = mean(PFpeakMou(p1, p2, :));
        meanPFsd(p1,p2) = mean(SD_meanPF(p1, p2, :));
        if STDPrule == 0
            meanFRin(p1,p2) = FRin(p1);
        end

        % bootstrap of mean slope and FRout
        Slopes = squeeze(COMslope(p1, p2, :));
        FRouts = squeeze(FRout(p1, p2, :));
        PFpeaks = squeeze(PFpeak(p1,p2,:));
        PFpeaksMou = squeeze(PFpeakMou(p1,p2,:));
        PFoutSD = squeeze(SD_meanPF(p1, p2, :));

        SlopeCI = bootci(1000, @mean, Slopes);
        AbsSlopeCI = bootci(1000, @mean, abs(Slopes));
        FRoutCI = bootci(1000, @mean, FRouts);
        PFpeakCI = bootci(1000, @mean, PFpeaks);
        PFpeakMouCI = bootci(1000, @mean, PFpeaksMou);
        PFoutSDCI = bootci(1000, @mean, PFoutSD);

        SlopeCIlo(p1, p2) = SlopeCI(1,1);
        SlopeCIup(p1, p2) = SlopeCI(2,1);
        AbsSlopeCIup(p1, p2) = AbsSlopeCI(2,1);
        AbsSlopeCIlo(p1, p2) = AbsSlopeCI(1,1);
        PFoutSDCIlo(p1, p2) = PFoutSDCI(1,1);
        PFoutSDCIup(p1, p2) = PFoutSDCI(2,1);
        SlopeCIup(p1, p2) = SlopeCI(2,1);
        FRoutCIlo(p1, p2) = FRoutCI(1,1);
        FRoutCIup(p1, p2) = FRoutCI(2,1);
        PFpeakCIlo(p1, p2) = PFpeakCI(1,1);
        PFpeakCIup(p1, p2) = PFpeakCI(2,1);
        PFpeakMouCIlo(p1, p2) = PFpeakMouCI(1,1);
        PFpeakMouCIup(p1, p2) = PFpeakMouCI(2,1);
        clear Slopes FRouts PFpeaks SlopeCI AbsSlopeCI FRoutCI PFpeakCI PFpeakMouCI

        % compute shifting PFs proportions
        Idx.Fwd = find(shiftPval(p1, p2, :)<=0.05 & COMslope(p1, p2, :)>0);
        Idx.Bck = find(shiftPval(p1, p2, :)<=0.05 & COMslope(p1, p2, :)<0);
        Idx.NonSig = find(shiftPval(p1, p2, :)>0.05);
        PropFwd(p1,p2) = length(Idx.Fwd)/Nsim;
        PropBck(p1,p2) = length(Idx.Bck)/Nsim;
        PropNonSig(p1,p2) = length(Idx.NonSig)/Nsim;
        clear Idx

    end
    end

% Overall prop of shifting PFs
Idx.Fwd = find(shiftPval<0.05 & COMslope>0);
Idx.Bck = find(shiftPval<0.05 & COMslope<0);
Idx.Non = find(shiftPval>0.05); 
PropTotFwd = length(Idx.Fwd)/(Nsim*p1*p2);
PropTotBck = length(Idx.Bck)/(Nsim*p1*p2);
PropTotNon = length(Idx.Non)/(Nsim*p1*p2);

toc

%% plots

% fwd PFs should be red, Bck PFs blue, NonSig grey or black
blue = [0 0.45 0.74]; % backward
red = [0.85 0.32 0.01]; % forward
grey = [0.5 0.5 0.5]; %nonsig

if STDPrule == 1 % exploring STDP rules

    yval = split(num2str(param1));
    xval = split(num2str(param2));
    
    figure(1) % stdp rules
        for p1 = 1:length(decay)
            for p2 = 1:length(Amp)
                leg(p1,p2) = strcat(yval(p1), ' ms - ', xval(p2), ' %');
                plot([delay1 delay2], [Wd{p1,p2} Wp{p1,p2}], 'k');hold on
%                 plot([delay1 delay2], [Wd{2,1} Wp{2,1}], 'g', 'LineWidth', 2);hold on
            end
        end
%         legend(leg(:), 'Location', 'BestOutside')
        xline(0);
        yline(0); 
        xlabel('pre-post delay (ms)'); ylabel('Weight change (% of EPSC_{max})');
        box off; axis square;
%         subplot(1,2,2)
%             for i = 1:length(Amp)
%             plot([delay1 delay2], [Wd{2,i} Wp{2,i}] );hold on
%             end
%             legend(split(num2str(Amp)), 'Location', 'BestOutside')
%             xline(0); hold on
%             yline(0); 
%             xlim([-100 100]);
%             xlabel('delay (ms)'); ylabel('Weight change (% of EPSC_{max})');
%             title('Tau_{STDP} = 20 ms')
%             box off; axis square;

    figure(2) % shifting slopes 
    % imagesc([Amp(1) Amp(end)], [decay(1) decay(end)], meanSlope); hold on
    imagesc(meanSlope); hold on
    set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
    xlabel('max Weight change (%)'); ylabel('Tau_{STDP} (ms)');
    colormap(brewermap(256,'*RdBu')); 
%     caxis([-ceil(max(meanSlope(:))), ceil(max(meanSlope(:)))]);
    caxis([-(max(abs(meanSlope(:)))), max(abs(meanSlope(:)))]);
    c = colorbar; c.Label.String = 'Shifting slope (cm/lap)';
    for p1 = 1:length(param1)-1
        yline(p1+0.5, 'k-', 'LineWidth', 0.5);
    end
    for p2 = 1:length(param2)-1
        xline(p2+0.5, 'k-', 'LineWidth', 0.5);
    end
    axis square
    
    figure(3) % FRout 
    imagesc(meanFRout); hold on
    set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
    xlabel('max Weight change (%)'); ylabel('Tau_{STDP} (ms)');
    colormap(brewermap(256,'oranges')) % colormap(parula(256));
    % caxis([0 ceil(max(meanFRout, [], 'all'))]);
    c = colorbar; c.Label.String = 'FR out (Hz)';
    for p1 = 1:length(param1)-1
        yline(p1+0.5, 'k-', 'LineWidth', 0.5);
    end
    for p2 = 1:length(param2)-1
        xline(p2+0.5, 'k-', 'LineWidth', 0.5);
    end
    axis square

    figure(4) % peak PF 
    imagesc(meanPFpeak); hold on
    set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
    xlabel('max Weight change (%)'); ylabel('Tau_{STDP} (ms)');
    colormap(brewermap(256,'reds')); %colormap(hot(256)); 
    % caxis([0 ceil(max(meanPFpeak, [], 'all'))]);
    c = colorbar; c.Label.String = 'PF peak (Hz)';
    for p1 = 1:length(param1)-1
        yline(p1+0.5, 'k-', 'LineWidth', 0.5);
    end
    for p2 = 1:length(param2)-1
        xline(p2+0.5, 'k-', 'LineWidth', 0.5);
    end
    axis square

    figure(5) % output PF width 
    imagesc(meanPFsd); hold on
    set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
    xlabel('max Weight change (%)'); ylabel('Tau_{STDP} (ms)');
    colormap(brewermap(256,'purples')); 
    % caxis([0 ceil(max(meanPFsd, [], 'all'))]);
    c = colorbar; c.Label.String = 'output PF s.d. (cm)';
    for p1 = 1:length(param1)-1
        yline(p1+0.5, 'k-', 'LineWidth', 0.5);
    end
    for p2 = 1:length(param2)-1
        xline(p2+0.5, 'k-', 'LineWidth', 0.5);
    end
    axis square

    % props (bar plots) for each param set
    figure(6)
    PrFwd = PropFwd(:); %reshape into a column vector, where each column was concatenated
    PrBck = PropBck(:); 
    PrNon = PropNonSig(:);
%     PrFwd = reshape(PropFwd', [], 1); %reshape into a column vector, where each column was concatenated
%     PrBck = reshape(PropBck', [], 1); 
%     PrNon = reshape(PropNonSig', [], 1);
    Props = [PrBck PrFwd PrNon];
    b = bar(Props, 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylabel('proportion of PFs'); xlabel('models (Tau _ max \DeltaW)')
    set(gca, 'XTick', 1:length(modelgp(:)), 'XTickLabel', modelname(:));
    box off;

    
    COMslopeAll = COMslope(:); FRoutAll = FRout(:);

    % Frout (or PFpeak) vs mean Slope with bootstrapped CI
    % MarkerEdgeColor of individual datapoints in black for significant shifts
    figure(7)
    x = meanFRout(:); y = meanSlope(:); sig = shiftPval(:);
    xneg = x-FRoutCIlo(:); xpos = FRoutCIup(:)-x;
    yneg = y-SlopeCIlo(:); ypos = SlopeCIup(:)-y;
    scatter(FRoutAll, COMslopeAll, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    scatter(FRoutAll(sig<0.05), COMslopeAll(sig<0.05), 'o', 'MarkerEdgeColor', 'k'); hold on
    errorbar(x,y,yneg,ypos,xneg,xpos,'o', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'r'); hold on
    yline(0,'--'); hold on
    % yline(mean(CA1slopes.N),'g');
%     ylim([-0.4 0.4])
    xlabel('FR out (Hz)'); ylabel('Shifting slope (cm/lap)'); 
    box off; axis square;   

    figure(8)
    clear x xneg xpos
    x = meanPFpeak(:); PFpeakAll = PFpeak(:);
    xneg = x-PFpeakCIlo(:); xpos = PFpeakCIup(:)-x;
    scatter(PFpeakAll, COMslopeAll, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    scatter(PFpeakAll(sig<0.05), COMslopeAll(sig<0.05), 'o', 'MarkerEdgeColor', 'k'); hold on
    errorbar(x,y,yneg,ypos,xneg,xpos,'o', 'Color', 'r', 'MarkerFaceColor', 'r'); hold on
    yline(0,'--');
%     ylim([-0.4 0.4])
    xlabel('peak PFout (Hz)'); ylabel('Shifting slope (cm/lap)'); 
    box off; axis square;
    
    % do the same 2 figures with absolute slope
    figure(9)
    clear x xneg xpos y yneg ypos
    x = meanFRout(:); y = meanAbsSlope(:); sig = shiftPval(:);
    AbsCOMslopeAll = abs(COMslopeAll); FRoutAll = FRout(:);
    xneg = x-FRoutCIlo(:); xpos = FRoutCIup(:)-x;
    yneg = y-AbsSlopeCIlo(:); ypos = AbsSlopeCIup(:)-y;
    scatter(FRoutAll, AbsCOMslopeAll, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    scatter(FRoutAll(sig<0.05), AbsCOMslopeAll(sig<0.05), 'o', 'MarkerEdgeColor', 'k'); hold on
    errorbar(x,y,yneg,ypos,xneg,xpos,'o', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'r'); hold on
    yline(mean(abs(CA1slopes.N)),'g');
%     ylim([0 0.5])
    xlabel('FR out (Hz)'); ylabel('Absolute shifting slope (cm/lap)'); 
    box off; axis square;
    
    figure(10)
    clear x xneg xpos
    x = meanPFpeak(:); PFpeakAll = PFpeak(:);
    xneg = x-PFpeakCIlo(:); xpos = PFpeakCIup(:)-x;
    scatter(PFpeakAll, AbsCOMslopeAll, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    scatter(PFpeakAll(sig<0.05), AbsCOMslopeAll(sig<0.05), 'o', 'MarkerEdgeColor', 'k'); hold on
    errorbar(x,y,yneg,ypos,xneg,xpos,'o', 'Color', 'r', 'MarkerFaceColor', 'r'); hold on
    yline(mean(abs(CA1slopes.N)),'g');
%     ylim([0 0.4])
    xlabel('peak PFout (Hz)'); ylabel('Absolute shifting slope (cm/lap)'); 
    box off; axis square;

    % props (bar plots) for all
    figure(11)
    PropsAll = [PropTotBck PropTotFwd PropTotNon];
    b2 = bar(1, PropsAll, 'stacked', 'FaceColor','flat');
    b2(1).CData = blue; b2(2).CData = red; b2(3).CData = grey;
    ylabel('proportion of PFs');
    % set(gca, 'XTick', [1:length(PrFwd)], 'XTickLabel', ['CA1-N';'CA1-F'; 'CA3-N'; 'CA3-F']);
    box off;
    
    % slopes vs R2, colored by signif 
    figure(12)
    scatter(COMslope(Idx.Bck), shiftR2(Idx.Bck), 100, blue, 'filled', 'MarkerFaceAlpha',0.6); hold on
    scatter(COMslope(Idx.Fwd), shiftR2(Idx.Fwd), 100, red, 'filled', 'MarkerFaceAlpha',0.6); hold on
    scatter(COMslope(Idx.Non), shiftR2(Idx.Non), 100, grey, 'filled', 'MarkerFaceAlpha',0.6); hold on
    xlabel('Shifting slope (cm/lap)'); ylabel('R-square');
    title('All params sets')
    box off; axis square;
    
    figure(13)
    scatter(COMslope(Idx.Bck), shiftR2(Idx.Bck), 'MarkerEdgeColor', blue); hold on
    scatter(COMslope(Idx.Fwd), shiftR2(Idx.Fwd), 'MarkerEdgeColor', red); hold on
    scatter(COMslope(Idx.Non), shiftR2(Idx.Non), 'MarkerEdgeColor', grey); hold on
    xlabel('Shifting slope (cm/lap)'); ylabel('R-square');
    xlim([-1.5 1.5])
    title('All params sets')
    box off; axis square;

    figure(14) % violin plots of each models aligned
%     violinplot(COMslopeAll, modelgpMat(:), 'ShowMean', true); hold on
%     violinplot(COMslopeAll, modelnameMat(:), 'ShowMean', true); hold on
    violinplot(COMslopeAll, cat(1,modelnameMat{:}), 'GroupOrder', reshape(modelname', [], 1), 'ShowMean', true); hold on
%     y = meanSlope(:); yneg = y-SlopeCIlo(:); ypos = SlopeCIup(:)-y;
    y = reshape(meanSlope', [], 1); yneg = y - reshape(SlopeCIlo', [], 1); ypos = reshape(SlopeCIup', [], 1)-y;
    errorbar(1:length(modelgp(:)), y,yneg,ypos,'.', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'r'); hold on
    yline(0, 'k--', 'linewidth', 2)
    xlabel('models (Tau _ max \DeltaW)'); ylabel('Shifting Slope (cm/lap)');
    xlim([0 length(modelgp(:))+1])
%     set(gca, 'XTick', 1:length(modelgp(:)), 'XTickLabel', modelname(:));
    box off;

%     figure
%     tbl.modelnameMat = cat(1,modelnameMat{:});
%     tbl.mdlnameMatcat = categorical(tbl.modelnameMat,reshape(modelname', [], 1));
%     tbl.COMslopeAll = COMslopeAll;
% %     T = struct2table(tbl);
%     
%     boxchart(tbl.modelgpMat(:),COMslopeAll); hold on
%     y = meanSlope(:); yneg = y-SlopeCIlo(:); ypos = SlopeCIup(:)-y;
%     errorbar(1:length(modelgp(:)), y,yneg,ypos,'.', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
%     xlabel('models (Tau _ max \DeltaW)'); ylabel('Shifting Slope (cm/lap)');
%     xlim([0 length(modelgp(:))+1])
%     set(gca, 'XTick', 1:length(modelgp(:)), 'XTickLabel', modelname(:));
%     box off;

    % slopes vs output PF sd (and delta sd)
    figure(15)
    clear x xneg xpos y yneg ypos
    x = meanPFsd(:); y = meanSlope(:);
    COMslopeAll = COMslope(:); SD_meanPFAll = SD_meanPF(:);
    xneg = x-PFoutSDCIlo(:); xpos = PFoutSDCIup(:)-x;
    yneg = y-SlopeCIlo(:); ypos = SlopeCIup(:)-y;
    scatter(SD_meanPFAll, COMslopeAll, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    scatter(SD_meanPFAll(sig<0.05), COMslopeAll(sig<0.05), 'o', 'MarkerEdgeColor', 'k'); hold on
    errorbar(x,y,yneg,ypos,xneg,xpos,'o', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'r'); hold on
    yline(0,'--'); hold on
    % yline(mean(CA1slopes.N),'g');
%     ylim([-0.4 0.4])
    xlabel('PFout s.d. (cm)'); ylabel('Shifting slope (cm/lap)'); 
    box off; axis square

    figure % best to plot when param 1 is maintained constant
    x = repmat(param2', 1, length(param1)); x = x(:); y = meanSlope(:);
    param2Matall = param2Mat(:);
    scatter(param2Matall, COMslopeAll, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    scatter(param2Matall(sig<0.05), COMslopeAll(sig<0.05), 'o', 'MarkerEdgeColor', 'k'); hold on
    errorbar(x,y,yneg,ypos,'o', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'r'); hold on
    yline(0,'--'); hold on
    xlabel('max weight change (% EPSCmax)'); ylabel('Shifting slope (cm/lap)'); 
    box off; axis square;

tl =  tiledlayout(2,4, 'TileSpacing','Compact','Padding','Compact');
title(tl,['PFin peak = ' num2str(Params.PFamp) ' Hz, Speed = ' num2str(Params.L/Params.period) ' cm/s'])
    nexttile
        for p1 = 1:length(decay)
            for p2 = 1:length(Amp)
                plot([delay1 delay2], [Wd{p1,p2} Wp{p1,p2}], 'k');hold on
            end
        end
        xline(0);
        yline(0); 
        xlabel('pre-post delay (ms)'); ylabel('Weight change (% of EPSC_{max})');
        box off; axis square;
    nexttile
        imagesc(meanSlope); hold on
        set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
        xlabel('max Weight change (%)'); ylabel('Tau_{STDP} (ms)');
        colormap(brewermap(256,'*RdBu')); 
    %     caxis([-ceil(max(meanSlope(:))), ceil(max(meanSlope(:)))]);
        caxis([-(max(abs(meanSlope(:)))), max(abs(meanSlope(:)))]);
        c = colorbar; c.Label.String = 'Shifting slope (cm/lap)';
        for p1 = 1:length(param1)-1
            yline(p1+0.5, 'k-', 'LineWidth', 0.5);
        end
        for p2 = 1:length(param2)-1
            xline(p2+0.5, 'k-', 'LineWidth', 0.5);
        end
        axis square
     nexttile([1 2])
        PrFwd2 = reshape(PropFwd', [], 1); %reshape into a column vector, where each column was concatenated
        PrBck2 = reshape(PropBck', [], 1); 
        PrNon2 = reshape(PropNonSig', [], 1);
        Props2 = [PrBck2 PrFwd2 PrNon2];
        b2 = bar(Props2, 'stacked', 'FaceColor','flat');
        b2(1).CData = blue; b2(2).CData = red; b2(3).CData = grey;
        ylabel('proportion of PFs'); 
%         xlabel('models (Tau - max \DeltaW)')
%         set(gca, 'XTick', 1:length(modelgp(:)), 'XTickLabel', reshape(modelname', [], 1));
        set(gca, 'XTick', [], 'XTickLabel', []);
        box off;
    nexttile
        clear x xneg xpos
        x = meanPFpeak(:); PFpeakAll = PFpeak(:);
        xneg = x-PFpeakCIlo(:); xpos = PFpeakCIup(:)-x;
        scatter(PFpeakAll, COMslopeAll, 15, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
        scatter(PFpeakAll(sig<0.05), COMslopeAll(sig<0.05), 15, 'o', 'MarkerEdgeColor', 'k'); hold on
        errorbar(x,y,yneg,ypos,xneg,xpos,'.', 'MarkerSize', 8, 'LineWidth', 0.75, 'Color', 'r', 'MarkerFaceColor', 'r'); hold on
        yline(0,'k--', 'linewidth', 1);
        xlabel('peak PFout (Hz)'); ylabel('Shifting slope (cm/lap)'); 
        box off; axis square; 
     nexttile
        clear x xneg xpos y yneg ypos
        x = meanPFsd(:); y = meanSlope(:);
        COMslopeAll = COMslope(:); SD_meanPFAll = SD_meanPF(:);
        xneg = x-PFoutSDCIlo(:); xpos = PFoutSDCIup(:)-x;
        yneg = y-SlopeCIlo(:); ypos = SlopeCIup(:)-y;
        scatter(SD_meanPFAll, COMslopeAll, 15, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
        scatter(SD_meanPFAll(sig<0.05), COMslopeAll(sig<0.05),15, 'o', 'MarkerEdgeColor', 'k'); hold on
        errorbar(x,y,yneg,ypos,xneg,xpos,'.', 'Color', 'r', 'MarkerSize', 8,'LineWidth', 0.75, 'MarkerFaceColor', 'r'); hold on
        yline(0,'k--', 'linewidth', 1); hold on
        % yline(mean(CA1slopes.N),'g');
    %     ylim([-0.4 0.4])
        xlabel('PFout s.d. (cm)'); ylabel('Shifting slope (cm/lap)'); 
        box off; axis square
     nexttile([1,2])
        yline(0, 'k--', 'linewidth', 1)
        map = brewermap(8,'Paired'); %cool(4); %brewermap(4,'PiYG')
        mapCA1 = map(3:4,:);
        yline(mean(CA1slopes.N), '--', 'Color', mapCA1(2,:), 'LineWidth', 1);
        yline(mean(CA1slopes.F), '--', 'Color', mapCA1(1,:), 'LineWidth', 1);
        vp = violinplot(COMslopeAll, cat(1,modelnameMat{:}), 'GroupOrder', reshape(modelname', [], 1), 'ShowMean', true); hold on
        vcolors = lines(length(param2));
%         vcolors = vcolors([1:length(param2)-1, end], :);
        ax = gca;
        for v1 = 1:length(param2)
            for v2 = 1:length(param1)
                item = v2+(v1-1)*length(param2);
                vp(item).ViolinColor = vcolors(v2,:);
                ax.XTickLabel{item} = sprintf('\\color[rgb]{%f,%f,%f}%s', vcolors(v2,:), ax.XTickLabel{item});
            end
        end
        y = reshape(meanSlope', [], 1); yneg = y - reshape(SlopeCIlo', [], 1); ypos = reshape(SlopeCIup', [], 1)-y;
        errorbar(1:length(modelgp(:)), y,yneg,ypos,'.', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 0.75, 'MarkerFaceColor', 'r'); hold on
        xlabel('models (Tau - max \DeltaW)'); ylabel('Shifting Slope (cm/lap)');
        xlim([0 length(modelgp(:))+1])
        box off;

        figure
        yline(0, 'k--', 'linewidth', 1)
        map = brewermap(8,'Paired'); %cool(4); %brewermap(4,'PiYG')
        mapCA1 = map(3:4,:);
        yline(mean(CA1slopes.N), '--', 'Color', mapCA1(2,:), 'LineWidth', 1);
        yline(mean(CA1slopes.F), '--', 'Color', mapCA1(1,:), 'LineWidth', 1);
        vp = violinplot(COMslopeAll, cat(1,modelnameMat{:}), 'GroupOrder', reshape(modelname', [], 1), 'ShowMean', true); hold on
        vcolors = lines(length(param2)+1);
        vcolors = vcolors([1:length(param2)-1, end], :);
        ax = gca;
        for v1 = 1:length(param2)
            for v2 = 1:length(param1)
                item = v2+(v1-1)*length(param2);
                vp(item).ViolinColor = vcolors(v2,:);
                ax.XTickLabel{item} = sprintf('\\color[rgb]{%f,%f,%f}%s', vcolors(v2,:), ax.XTickLabel{item});
            end
        end
        y = reshape(meanSlope', [], 1); yneg = y - reshape(SlopeCIlo', [], 1); ypos = reshape(SlopeCIup', [], 1)-y;
        errorbar(1:length(modelgp(:)), y,yneg,ypos,'.', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 0.75, 'MarkerFaceColor', 'r'); hold on
        xlabel('models (Tau - max \DeltaW)'); ylabel('Shifting Slope (cm/lap)');
        xlim([0 length(modelgp(:))+1])
        ylim([-1 0.5])
        box off;
%         exportgraphics(tl,'ParamSpace_AmpDecay_PFin10HzSpeedMehta.pdf','BackgroundColor','none', 'ContentType','vector')

else % STDPrule == 0 -> exploring FRinput and output variations

    yval1 = split(num2str(PFsd));
    yval2 = repmat('cm-', length(param1),1);
    yval3 = split(num2str(PFamp));
    yval4 = repmat('Hz', length(param1),1);
    yval = strcat(yval1, yval2, yval3, yval4);
    xval = split(num2str(Wsd));
    FRinVal = split(num2str(FRin));
    
    figure(1)
    % imagesc([Wsd(1) Wsd(end)], [FRin(1) FRin(end)], meanFRin); hold on
    imagesc(meanFRin); hold on
    set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
    ylabel('input [PF s.d., PF peak]'); xlabel('Connectivity s.d. (neurons)');
    colormap(brewermap(256,'greys')); 
    c = colorbar; c.Label.String = 'FRin (Hz)';
    for p1 = 1:length(param1)-1
        yline(p1+0.5, 'k-', 'LineWidth', 0.5);
    end
    for p2 = 1:length(param2)-1
        xline(p2+0.5, 'k-', 'LineWidth', 0.5);
    end
    
    figure(2) % slopes
    % imagesc([Wsd(1) Wsd(end)], [FRin(1) FRin(end)], meanSlope); hold on
    % set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', FRinVal);
    % ylabel('FRin (Hz)'); 
    imagesc(meanSlope); hold on
    set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
    ylabel('input [PF s.d., PF peak]')
    xlabel('Connectivity s.d. (cm)');
    colormap(brewermap(256,'*RdBu')); 
    caxis([-0.4 0.4]);
    c = colorbar; c.Label.String = 'Shifting slope (cm/lap)';
    for p1 = 1:length(param1)-1
        yline(p1+0.5, 'k-', 'LineWidth', 0.5);
    end
    for p2 = 1:length(param2)-1
        xline(p2+0.5, 'k-', 'LineWidth', 0.5);
    end
    
    figure(3) % FRout 
    % imagesc([Wsd(1) Wsd(end)], [FRin(1) FRin(end)], meanFRout); hold on
    %ylabel('FR in (Hz)');
    imagesc(meanFRout); hold on
    set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
    ylabel('input [PF s.d., PF peak]')
    xlabel('Connectivity s.d.');
    colormap(brewermap(256,'oranges')) % colormap(parula(256));
    caxis([0 ceil(max(meanFRout, [], 'all'))]);
    % caxis([0 10]);
    c = colorbar; c.Label.String = 'FR out (Hz)';
    for p1 = 1:length(param1)-1
        yline(p1+0.5, 'k-', 'LineWidth', 0.5);
    end
    for p2 = 1:length(param2)-1
        xline(p2+0.5, 'k-', 'LineWidth', 0.5);
    end
    
    figure(4) % peak PF 
    % imagesc([Wsd(1) Wsd(end)], [FRin(1) FRin(end)], meanPFpeak); hold on
    % ylabel('FRin (Hz)');
    imagesc(meanPFpeak); hold on
    set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
    ylabel('input [PF s.d., PF peak]')
    xlabel('Connectivity s.d.');
    colormap(brewermap(256,'reds')); %colormap(hot(256)); 
    % caxis([0 ceil(max(meanPFpeak, [], 'all'))]);
    c = colorbar; c.Label.String = 'PF peak (Hz)';
    for p1 = 1:length(param1)-1
        yline(p1+0.5, 'k-', 'LineWidth', 0.5);
    end
    for p2 = 1:length(param2)-1
        xline(p2+0.5, 'k-', 'LineWidth', 0.5);
    end
    
    figure(5) % output PF width 
    % imagesc([Wsd(1) Wsd(end)], [FRin(1) FRin(end)], meanPFsd); hold on
    % ylabel('FRin (Hz)');
    imagesc(meanPFsd); hold on
    set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
    ylabel('input [PF s.d., PF peak]')
    xlabel('Connectivity s.d.');
    colormap(brewermap(256,'purples')); 
    % caxis([0 ceil(max(meanPFsd, [], 'all'))]);
    c = colorbar; c.Label.String = 'output PF s.d. (cm)';
    for p1 = 1:length(param1)-1
        yline(p1+0.5, 'k-', 'LineWidth', 0.5);
    end
    for p2 = 1:length(param2)-1
        xline(p2+0.5, 'k-', 'LineWidth', 0.5);
    end

    % props (bar plots) for each param set
    figure(6)
    PrFwd = PropFwd(:); %reshape into a column vector, where each column was concatenated
    PrBck = PropBck(:); 
    PrNon = PropNonSig(:);
    Props = [PrBck PrFwd PrNon];
    b = bar(Props, 'stacked', 'FaceColor','flat');
    b(1).CData = blue; b(2).CData = red; b(3).CData = grey;
    ylabel('proportion of PFs');
    % set(gca, 'XTick', [1:length(PrFwd)], 'XTickLabel', ['CA1-N';'CA1-F'; 'CA3-N'; 'CA3-F']);
    box off;
    
    % Frout (or PFpeak) vs mean Slope with bootstrapped CI
    % MarkerEdgeColor of individual datapoints in black for significant shifts
    figure(7)
    x = meanFRout(:); y = meanSlope(:); sig = shiftPval(:);
    COMslopeAll = COMslope(:); FRoutAll = FRout(:);
    xneg = x-FRoutCIlo(:); xpos = FRoutCIup(:)-x;
    yneg = y-SlopeCIlo(:); ypos = SlopeCIup(:)-y;
    scatter(FRoutAll, COMslopeAll, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    scatter(FRoutAll(sig<0.05), COMslopeAll(sig<0.05), 'o', 'MarkerEdgeColor', 'k'); hold on
    errorbar(x,y,yneg,ypos,xneg,xpos,'o', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'r'); hold on
    yline(0,'--'); hold on
    % yline(mean(CA1slopes.N),'g');
    ylim([-0.4 0.4])
    xlabel('FR out (Hz)'); ylabel('Shifting slope (cm/lap)'); 
    box off; axis square;
    
    figure(8)
    clear x xneg xpos
    x = meanPFpeak(:); PFpeakAll = PFpeak(:);
    xneg = x-PFpeakCIlo(:); xpos = PFpeakCIup(:)-x;
    scatter(PFpeakAll, COMslopeAll, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    scatter(PFpeakAll(sig<0.05), COMslopeAll(sig<0.05), 'o', 'MarkerEdgeColor', 'k'); hold on
    errorbar(x,y,yneg,ypos,xneg,xpos,'o', 'Color', 'r', 'MarkerFaceColor', 'r'); hold on
    yline(0,'--');
    ylim([-0.4 0.4])
    xlabel('peak PFout (Hz)'); ylabel('Shifting slope (cm/lap)'); 
    box off; axis square;
    
    % do the same 2 figures with absolute slope
    figure(9)
    clear x xneg xpos y yneg ypos
    x = meanFRout(:); y = meanAbsSlope(:); sig = shiftPval(:);
    AbsCOMslopeAll = abs(COMslopeAll); FRoutAll = FRout(:);
    xneg = x-FRoutCIlo(:); xpos = FRoutCIup(:)-x;
    yneg = y-AbsSlopeCIlo(:); ypos = AbsSlopeCIup(:)-y;
    scatter(FRoutAll, AbsCOMslopeAll, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    scatter(FRoutAll(sig<0.05), AbsCOMslopeAll(sig<0.05), 'o', 'MarkerEdgeColor', 'k'); hold on
    errorbar(x,y,yneg,ypos,xneg,xpos,'o', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'r'); hold on
    yline(mean(abs(CA1slopes.F)),'g');
    ylim([0 0.5])
    xlabel('FR out (Hz)'); ylabel('Absolute shifting slope (cm/lap)'); 
    box off; axis square;
    
    figure(10)
    clear x xneg xpos
    x = meanPFpeak(:); PFpeakAll = PFpeak(:);
    xneg = x-PFpeakCIlo(:); xpos = PFpeakCIup(:)-x;
    scatter(PFpeakAll, AbsCOMslopeAll, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    scatter(PFpeakAll(sig<0.05), AbsCOMslopeAll(sig<0.05), 'o', 'MarkerEdgeColor', 'k'); hold on
    errorbar(x,y,yneg,ypos,xneg,xpos,'o', 'Color', 'r', 'MarkerFaceColor', 'r'); hold on
    yline(mean(abs(CA1slopes.F)),'g');
    ylim([0 0.4])
    xlabel('peak PFout (Hz)'); ylabel('Absolute shifting slope (cm/lap)'); 
    box off; axis square;
    
    % the slopes are larger than I expected for low FRout... 
    % Relationship between with Wsd, FRin, PFsd and PFamp:
    
    for p1 = 1:length(param1)
        for p2 = 1:length(param2)
            for n = 1:Nsim
                WsdMat(p1,p2,n) = Wsd(p2);
                PFampMat(p1,p2,n) = PFamp(p1);
                PFsdMat(p1,p2,n) = PFsd(p1);
                FRinMat(p1,p2,n) = FRin(p1);
            end
        end
    end
    
    figure(11)
    scatter3(FRoutAll, WsdMat(:), COMslopeAll,  'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    xlabel('FR out (Hz)'); zlabel('Shifting slope (cm/lap)'); ylabel('Connectivity S.D. (cm)')
    box off; axis square;
    
    figure(12)
    scatter3(FRinMat(:), WsdMat(:), COMslopeAll,  'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    xlabel('FR in (Hz)'); zlabel('Shifting slope (cm/lap)'); ylabel('Connectivity S.D. (cm)')
    box off; axis square;
    
    figure(13)
    scatter3(PFampMat(:), PFsdMat(:), COMslopeAll,  'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    xlabel('input PF peak (Hz)'); zlabel('Shifting slope (cm/lap)'); ylabel('input PF width (cm)')
    box off; axis square;
    
    % props (bar plots) for all
    figure(14)
    PropsAll = [PropTotBck PropTotFwd PropTotNon];
    b2 = bar(1, PropsAll, 'stacked', 'FaceColor','flat');
    b2(1).CData = blue; b2(2).CData = red; b2(3).CData = grey;
    ylabel('proportion of PFs');
    % set(gca, 'XTick', [1:length(PrFwd)], 'XTickLabel', ['CA1-N';'CA1-F'; 'CA3-N'; 'CA3-F']);
    box off;
    
    % slopes vs R2, colored by signif 
    figure(15)
    scatter(COMslope(Idx.Bck), shiftR2(Idx.Bck), 100, blue, 'filled', 'MarkerFaceAlpha',0.6); hold on
    scatter(COMslope(Idx.Fwd), shiftR2(Idx.Fwd), 100, red, 'filled', 'MarkerFaceAlpha',0.6); hold on
    scatter(COMslope(Idx.Non), shiftR2(Idx.Non), 100, grey, 'filled', 'MarkerFaceAlpha',0.6); hold on
    xlabel('Shifting slope (cm/lap)'); ylabel('R-square');
    title('All params sets')
    box off; axis square;
    
    figure(16)
    scatter(COMslope(Idx.Bck), shiftR2(Idx.Bck), 'MarkerEdgeColor', blue); hold on
    scatter(COMslope(Idx.Fwd), shiftR2(Idx.Fwd), 'MarkerEdgeColor', red); hold on
    scatter(COMslope(Idx.Non), shiftR2(Idx.Non), 'MarkerEdgeColor', grey); hold on
    xlabel('Shifting slope (cm/lap)'); ylabel('R-square');
    xlim([-1.5 1.5])
    title('All params sets')
    box off; axis square;

    figure(17)
    clear x xneg xpos y ypos yneg
    x = meanPFpeak(:); y = meanSlope(:);
    xMou = meanPFpeakMou(:);
    PFpeakAll = PFpeak(:);
    PFpeakMouAll = PFpeakMou(:);
    xneg = x-PFpeakCIlo(:); xpos = PFpeakCIup(:)-x;
    xMouneg = xMou-PFpeakMouCIlo(:); xMoupos = PFpeakMouCIup(:)-xMou;
    Xreg = [PFpeakAll, ones(size(PFpeakAll))];
    B = regress(PFpeakMouAll,Xreg);
    Xplot = [0, 1; max(PFpeakAll)+1, 1];
    Yplot = Xplot*B;
    Yplot2 = B(1)*Xplot(:,1)+B(2);
    Yreal = [4 35];
    Xreal = (Yreal - B(2))/B(1);
    plot(Xplot(:,1), Yplot2, 'r--', 'LineWidth', 1); hold on
    plot([0 max([x ; xMou])], [0 max([x ; xMou])], 'k--', 'LineWidth', 1); hold on
    for r = 1:2
    plot([0 Xreal(r)], [Yreal(r) Yreal(r)], 'g-', 'LineWidth', 1)
    plot([Xreal(r) Xreal(r)], [0 Yreal(r)], 'g-', 'LineWidth', 1)
    end
%     yline(4,'g--', 'LineWidth', 1)
%     yline(35,'g--', 'LineWidth', 1)
    scatter(PFpeakAll, PFpeakMouAll, 15, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
    errorbar(x,xMou,xMouneg,xMoupos,xneg,xpos,'.', 'MarkerSize', 8, 'LineWidth', 0.75, 'Color', 'r', 'MarkerFaceColor', 'r'); hold on
    xlabel('peak PFout with 6cm bins (Hz)'); ylabel('peak PFout with 2.5cm bins (Hz)'); 
    xlim([0 110])
    title(['Y = ' num2str(B(1)) 'X + ' num2str(B(2))])
    box off; 
    axis square;

    tl =  tiledlayout(2,4, 'TileSpacing','Compact','Padding','Compact');
    if Params.InShift == 1
        DynIn = ['CA3' Context];
    else
        DynIn = 'No';
    end
    title(tl,['Dynamic Inputs: ' DynIn ', STDP tau = ' num2str(Params.Pdecay*1000) 'ms, dWmax = ' num2str(Params.Pamp*100) '%, EPSCmax = ' num2str(Params.Imax*10^12) 'pA, Speed = ' num2str(Params.L/Params.period) ' cm/s'])
%     nexttile
%         delay1 = -100:1:0; % ms
%         delay2 = 0:1:100; %ms
%         Wd = (Params.Damp*100)*exp(delay1/(Params.Ddecay*1000));
%         Wp = (Params.Pamp*100)*exp(-delay2/(Params.Pdecay*1000));
%         plot([delay1 delay2], [Wd Wp], 'k'); hold on
%         xline(0);
%         yline(0); 
%         xlabel('pre-post delay (ms)'); ylabel('Weight change (% of EPSC_{max})');
%         box off; axis square;        
%     ax1 = nexttile;
%         imagesc(meanFRin); hold on
%         set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
%         ylabel('input [PF s.d., PF peak]'); xlabel('Connectivity s.d. (neurons)');
%         colormap(ax1, brewermap(256,'greys')); 
%         c = colorbar; c.Label.String = 'FRin (Hz)';
%         for p1 = 1:length(param1)-1
%             yline(p1+0.5, 'k-', 'LineWidth', 0.5);
%         end
%         for p2 = 1:length(param2)-1
%             xline(p2+0.5, 'k-', 'LineWidth', 0.5);
%         end
%         axis square
    ax2 = nexttile;
        imagesc(meanSlope); hold on
%         set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', 1:length(param1), 'YTickLabel', yval);
        set(gca, 'XTick', 1:length(param2), 'XTickLabel', xval, 'YTick', [], 'YTickLabel', []);
        ylabel('input [PF s.d., PF peak]'); xlabel('Connectivity s.d. (neurons)');
        colormap(ax2, brewermap(256,'*RdBu')); 
    %     caxis([-ceil(max(meanSlope(:))), ceil(max(meanSlope(:)))]);
        caxis([-(max(abs(meanSlope(:)))), max(abs(meanSlope(:)))]);
        c = colorbar; c.Label.String = 'COM slope (cm/lap)';
        for p1 = 1:length(param1)-1
            yline(p1+0.5, 'k-', 'LineWidth', 0.5);
        end
        for p2 = 1:length(param2)-1
            xline(p2+0.5, 'k-', 'LineWidth', 0.5);
        end
        axis square
    nexttile
        clear x xneg xpos y yneg ypos
        x = meanPFsd(:); y = meanSlope(:);
        COMslopeAll = COMslope(:); SD_meanPFAll = SD_meanPF(:);
%         xneg = x-PFoutSDCIlo(:); xpos = PFoutSDCIup(:)-x;
%         yneg = y-SlopeCIlo(:); ypos = SlopeCIup(:)-y;
        scatter(SD_meanPFAll, COMslopeAll, 15, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
        scatter(SD_meanPFAll(sig<0.05), COMslopeAll(sig<0.05),15, 'o', 'MarkerEdgeColor', 'k'); hold on
%         errorbar(x,y,yneg,ypos,xneg,xpos,'.', 'Color', 'r', 'MarkerSize', 8,'LineWidth', 0.75, 'MarkerFaceColor', 'r'); hold on
        yline(0,'k--', 'linewidth', 1); hold on
        % yline(mean(CA1slopes.N),'g');
    %     ylim([-0.4 0.4])
        xlabel('PFout s.d. (cm)'); ylabel('COM slope (cm/lap)'); 
        box off; axis square
     nexttile([1 2])
        PrFwd2 = reshape(PropFwd', [], 1); %reshape into a column vector, where each column was concatenated
        PrBck2 = reshape(PropBck', [], 1); 
        PrNon2 = reshape(PropNonSig', [], 1);
        Props2 = [PrBck2 PrFwd2 PrNon2];
        b2 = bar(Props2, 'stacked', 'FaceColor','flat');
        b2(1).CData = blue; b2(2).CData = red; b2(3).CData = grey;
        ylabel('proportion of PFs'); 
%         xlabel('models (Tau - max \DeltaW)')
%         set(gca, 'XTick', 1:length(modelgp(:)), 'XTickLabel', reshape(modelname', [], 1));
        set(gca, 'XTick', [], 'XTickLabel', []);
        box off;
%      nexttile([1,2])
%         clear x xneg xpos y ypos yneg
%         x = meanFRout(:); y = meanSlope(:); sig = shiftPval(:);
%         COMslopeAll = COMslope(:); FRoutAll = FRout(:);
%         xneg = x-FRoutCIlo(:); xpos = FRoutCIup(:)-x;
%         yneg = y-SlopeCIlo(:); ypos = SlopeCIup(:)-y;
%         scatter(FRoutAll, COMslopeAll, 15, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
%         scatter(FRoutAll(sig<0.05), COMslopeAll(sig<0.05), 15, 'o', 'MarkerEdgeColor', 'k'); hold on
%         errorbar(x,y,yneg,ypos,xneg,xpos,'.', 'MarkerSize', 8, 'LineWidth', 0.75, 'Color', 'r', 'MarkerFaceColor', 'r'); hold on
%         yline(0,'--'); hold on
%         % yline(mean(CA1slopes.N),'g');
%         xlabel('FR out (Hz)'); ylabel('COM slope (cm/lap)'); 
%         box off;
     nexttile([1,2])
        clear x xneg xpos y ypos yneg
        x = meanPFpeak(:); y = meanSlope(:);
        PFpeakAll = PFpeak(:);
        xneg = x-PFpeakCIlo(:); xpos = PFpeakCIup(:)-x;
        yneg = y-SlopeCIlo(:); ypos = SlopeCIup(:)-y;
        scatter(PFpeakAll, COMslopeAll, 15, 'o', 'filled', 'MarkerFaceColor', grey, 'MarkerFaceAlpha', 0.5); hold on
        scatter(PFpeakAll(sig<0.05), COMslopeAll(sig<0.05), 15, 'o', 'MarkerEdgeColor', 'k'); hold on
        errorbar(x,y,yneg,ypos,xneg,xpos,'.', 'MarkerSize', 8, 'LineWidth', 0.75, 'Color', 'r', 'MarkerFaceColor', 'r'); hold on
        set(gca, 'XTick', [0 5 10 20:20:120])
        yline(0,'k--', 'linewidth', 1);
        xlabel('peak PFout (Hz)'); ylabel('COM slope (cm/lap)'); 
        box off; 
        %axis square; 
     nexttile([1,2])
        yline(0, 'k--', 'linewidth', 1)
        map = brewermap(8,'Paired'); %cool(4); %brewermap(4,'PiYG')
        mapCA1 = map(3:4,:);
        yline(mean(CA1slopes.N), '--', 'Color', mapCA1(2,:), 'LineWidth', 1);
%         yline(mean(CA1slopes.F), '--', 'Color', mapCA1(1,:), 'LineWidth', 1);
        vp = violinplot(COMslopeAll, cat(1,modelnameMat{:}), 'GroupOrder', reshape(modelname', [], 1), 'ShowMean', true); hold on
        y = reshape(meanSlope', [], 1); yneg = y - reshape(SlopeCIlo', [], 1); ypos = reshape(SlopeCIup', [], 1)-y;
        errorbar(1:length(modelgp(:)), y,yneg,ypos,'.', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 0.75, 'MarkerFaceColor', 'r'); hold on
        xlabel('models (PFin sd, PFin peak, Connectivity sd)'); ylabel('COM Slope (cm/lap)');
        xlim([0 length(modelgp(:))+1])
        box off;

end
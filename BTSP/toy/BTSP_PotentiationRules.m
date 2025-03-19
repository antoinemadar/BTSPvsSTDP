clear
close all

PampRange = [4,6,8,10,11,12,13,15,20, 25.5, 30,60,80]; %in Amp
Params.Pdecay = 1.31; % Pre-before-Post time constant, in sec (1.31s on average in Bittner et al. 2017)
Params.Ddecay = 0.69; % Post-before-Pre time constant, in sec (0.69s on average in Bittner et al. 2017)

f = figure; %STDP rule   
    delay1 = -5000:1:0; % ms
    delay2 = 0:1:5000; %ms
    for n = 1:length(PampRange)
    Wd = (PampRange(n))*exp(-delay2/(Params.Ddecay*1000));
    Wp = (PampRange(n))*exp(delay1/(Params.Pdecay*1000));
    plot([delay1 delay2]./1000, [Wp Wd], 'k'); hold on
%     clear Wd Wp
    end
%     xline(0, 'k--'); hold on
%     yline(0, 'k--'); 
    xlim([-4.2 4.2])
%     ylim([0 22])
    xlabel('delay (s)'); ylabel('Weight change (pA)');
    box off; axis square;

f.Renderer = 'painters'   
print('BTSP_PotentiationRules.pdf','-dpdf')
%% script to replot data from runTRESimulations.

for testcase = 1:8
    switch(testcase)
        case 1
            casename = 'Des_d_Iso_85mmTip';
            probeDesign = 'd';
            %Sigma = SigmaIso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
        case 2
            casename = 'Des_d_Iso_170mmTip';
            probeDesign = 'd';
            %Sigma = SigmaIso;
            A = 71;
            B = 54;
            rho = 170;
            % define West tool design parameters.
            r = 32;
        case 3
            casename = 'Des_e_Iso_85mmTip';
            probeDesign = 'e';
            %Sigma = SigmaIso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
        case 4
            casename = 'Des_e_Iso_170mmTip';
            probeDesign = 'e';
            %Sigma = SigmaIso;
            A = 71;
            B = 54;
            rho = 170;
            % define West tool design parameters.
            r = 32;
        case 5
            casename = 'Des_d_Aniso_85mmTip';
            probeDesign = 'd';
            %Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
        case 6
            casename = 'Des_d_Aniso_170mmTip';
            probeDesign = 'd';
            %Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 170;
            % define West tool design parameters.
            r = 32;
        case 7
            casename = 'Des_e_Aniso_85mmTip';
            probeDesign = 'e';
            %Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 85;
            % define West tool design parameters.
            r = 32;
        case 8
            casename = 'Des_e_Aniso_170mmTip';
            probeDesign = 'e';
            %Sigma = SigmaAniso;
            A = 71;
            B = 54;
            rho = 170;
            % define West tool design parameters.
            r = 32;
        otherwise
            error('Ummm... that case does not exist');
    end

    cd(casename)
    filename = sprintf('%s.csv',casename);
    DataIn = csvread(filename, 0, 2)
    
    PredictedRMS = DataIn(:,2);
    MeasuredRMS = DataIn(:,3);
    RMSPercentDiff = DataIn(:,4);
    
    figFontSize = 18;
    
    figure(1);
    hist(RMSPercentDiff);
    axis([-3,3,0,300]);
    titlestring = sprintf('Case %d: Histogram of the %% Differences for RMS', testcase);
    title(titlestring, 'fontsize',figFontSize);
    set(gca,'fontsize',figFontSize);
    figurefilename = sprintf('Histogram_%s', casename);
    print('-depsc', '-tiff', '-r300', figurefilename);
    
    figure(2);
    plot(PredictedRMS, MeasuredRMS, '.');
    hold on;
    meanRMSPercentDiff = mean(RMSPercentDiff);
    stdRMSPercentDiff = std(RMSPercentDiff);
    maxValue = max([PredictedRMS; MeasuredRMS]);
    minValue = min([PredictedRMS; MeasuredRMS]);
    [R,P] = corrcoef([PredictedRMS, MeasuredRMS]);
    rho = R(1,2);
    plot([minValue maxValue], [minValue, maxValue], 'k');
    xlabel('Predicted TRE RMS (mm)', 'fontsize',figFontSize);
    ylabel('Simulated TRE RMS (mm)', 'fontsize',figFontSize);
    set(gca,'fontsize',figFontSize);
    titlestring = sprintf('Case %d: Comparison of Predicted and Simulated TRE RMS', testcase);
    subtitlestring = sprintf('\\rho = %2.3f', rho);
    title({titlestring, subtitlestring}, 'fontsize',figFontSize);
    hold off;
    axis([minValue, maxValue, minValue, maxValue]);
    figurefilename = sprintf('PredvMeas_%s', casename);
    print('-depsc', '-tiff', '-r300', figurefilename);
    
    copyfile('*.eps', 'E:\docs\research\phd\publications\journal_papers\IEEE\anisotropic_TRE' );
    
    cd ..
end
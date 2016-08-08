%% Plot modeled CSTD with comparision to measured CSTD

% Transformation of Conc_monthly in Conc_monthly_cstd

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cstd_age = [25]+1;               %%  % +1 because age 0 is the first columne.             
        cstd_age_max = 60+1;             %%  % +1 because age 0 is included
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % automatic plot code 
        cstd_tm = (emp_year_all - Y_begin)*M; 

        Conc_monthly_cstd  = zeros(cstd_age_max, length(emp_year_all));

        for i = 1:length(emp_year_all)
            Conc_monthly_cstd(:,i)  = flipud(Conc_monthly_all_0to80(emp_year_all(i)-cstd_age_max-Y_begin+2+people/2:emp_year_all(i)-Y_begin+1+people/2, cstd_tm(i)));
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(chooseOpt,'intake') 

    % Choose congener
    if chooseCongener == 1     
        OPT_chem_PCB28;
        congener = 'PCB28'; 
    elseif chooseCongener == 2 
        OPT_chem_PCB52;
        congener = 'PCB52'; 
    elseif chooseCongener == 3 
        OPT_chem_PCB118;
        congener = 'PCB118'; 
    elseif chooseCongener == 4 
        OPT_chem_PCB138;
        congener = 'PCB138'; 
    elseif chooseCongener == 5 
        OPT_chem_PCB153;
        congener = 'PCB153'; 
    elseif chooseCongener == 6 
        OPT_chem_PCB180;
        congener = 'PCB180'; 
    elseif chooseCongener == 7
        OPT_chem_PCB170;
        congener = 'PCB170'; 
    elseif chooseCongener == 8 
        OPT_chem_HCB;
        congener = 'HCB'; 
    elseif chooseCongener == 9 
        OPT_chem_DDE;
        congener = 'DDE';     
    end
      
    figure1 = figure(chooseCongener+40); %clf
%     figure_title = [congener ': 5-year-mean CSD from 1994 until 2009'];        
%     annotation(figure1,'textbox', [0.35 0.97 0.82 0.03], 'String',{figure_title}, 'FontWeight','bold', ...
%         'FontSize',14, 'FitBoxToText','off', 'LineStyle','none');

    set(gcf,'Color',[1,1,1]); 
    set(gcf, 'PaperPositionMode','auto')     %important! don't know why... :-)
    set(0,'DefaultAxesFontName', 'Times New Roman'); set(0,'DefaultTextFontname', 'Times New Roman')
    set(0,'DefaultTextFontsize', 16); set(0,'DefaultAxesFontsize', 16)
    set(gcf, 'Units','centimeters','position', [-50 0 23 17])

    plot(emp_year_all, Conc_monthly_cstd(cstd_age,:),'*-r'); % fitted modeled CSTD
    hold on
    plot(CSTD_intake(:,1), CSTD_intake(:,chooseCongener+1), '*k'); % measured CSTD, +1 because first row is emp_year_all
    
    % calculate Tdec from modeled CSTD
    ppp = polyfit(emp_year_all,log(Conc_monthly_cstd(cstd_age,:)),1);
    Tdec_modeled_cstd = log(2)/abs(ppp(1));

    title([congener ' with CSTD-based-Tdec = ' num2str(Tdec_modeled_cstd) ' years and input Tdec = ' num2str(Th_dec_down) ' years'] )         
    ylabel('concentration in ng/g lipid')
    xlabel('time')
    
elseif strcmp(chooseOpt,'original')
    
      % Choose congener
            if chooseCongener == 1     
                OPT_chem_PCB28;
                congener = 'PCB28'; 
            elseif chooseCongener == 2 
                OPT_chem_PCB52;
                congener = 'PCB52'; 
            elseif chooseCongener == 3 
                OPT_chem_PCB101;
                congener = 'PCB101'; 
            elseif chooseCongener == 4 
                OPT_chem_PCB118;
                congener = 'PCB118'; 
            elseif chooseCongener == 5 
                OPT_chem_PCB138;
                congener = 'PCB138'; 
            elseif chooseCongener == 6 
                OPT_chem_PCB153;
                congener = 'PCB153'; 
            elseif chooseCongener == 7 
                OPT_chem_PCB180;
                congener = 'PCB180'; 
            elseif chooseCongener == 8
                OPT_chem_PCB170;
                congener = 'PCB170'; 
            elseif chooseCongener == 9 
                OPT_chem_HCB;
                congener = 'HCB'; 
            elseif chooseCongener == 10 
                OPT_chem_bHCH;
                congener = 'bHCH';      
            elseif chooseCongener == 11 
                OPT_chem_gHCH;
                congener = 'gHCH'; 
            elseif chooseCongener == 12 
                OPT_chem_DDE;
                congener = 'DDE'; 
            end
            
    figure1 = figure(chooseCongener+40); %clf
%     figure_title = [congener ': 5-year-mean CSD from 1994 until 2009'];        
%     annotation(figure1,'textbox', [0.35 0.97 0.82 0.03], 'String',{figure_title}, 'FontWeight','bold', ...
%         'FontSize',14, 'FitBoxToText','off', 'LineStyle','none');

    set(gcf,'Color',[1,1,1]); 
    set(gcf, 'PaperPositionMode','auto')     %important! don't know why... :-)
    set(0,'DefaultAxesFontName', 'Times New Roman'); set(0,'DefaultTextFontname', 'Times New Roman')
    set(0,'DefaultTextFontsize', 16); set(0,'DefaultAxesFontsize', 16)
    set(gcf, 'Units','centimeters','position', [-50 0 23 17])

    plot(emp_year_all, Conc_monthly_cstd(cstd_age,:),'*-r'); % fitted modeled CSTD
    hold on
    plot(CSTD_original(:,1), CSTD_original(:,chooseCongener+1), '*-k');     
    
    % calculate Tdec from modeled CSTD
    ppp = polyfit(emp_year_all,log(Conc_monthly_cstd(cstd_age,:)),1);
    Tdec_modeled_cstd = log(2)/abs(ppp(1));
    
    title([congener ' with CSTD-based-Tdec = ' num2str(Tdec_modeled_cstd) ' years and fitted Tdec = ' num2str(Th_dec_down) ' years'] )         
    ylabel('concentration in ng/g lipid')
    xlabel('time')
end

    



% 
%     
%     for i = 1:length(emp_year)
%         subplot(4,4,i)
%         semilogy(age_group',emp_data_noNANs(:,i),'o'); hold on 
%         semilogy(age_group', mod_data_nonNaNs(:,i),'or');
%         semilogy((1/12:1:79+1/12)', mod_data_nonNaNs_all(:,i),'r');
%         xlim([0 80]); 
%         
%          if chooseCongener == 1     
%                 ylim([0.1 100]);  
%             elseif chooseCongener == 2 
%                 ylim([0.1 100]);  
%             elseif chooseCongener == 3 
%                 ylim([0.1 100]);  
%             elseif chooseCongener == 4 
%                 ylim([0.1 1000]); 
%             elseif chooseCongener == 5 
%                 ylim([0.1 1000]); 
%             elseif chooseCongener == 6 
%                 ylim([0.1 1000]);  
%             elseif chooseCongener == 7
%                 ylim([0.1 1000]);  
%             elseif chooseCongener == 8 
%                 ylim([0.1 1000]);
%             elseif chooseCongener == 9 
%                 ylim([0.1 1000]);
%             elseif chooseCongener == 10 
%                ylim([0.1 1000]);   
%         end       

     
%     end   
   
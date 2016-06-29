%% Plot cross-sectional data - only in combination with backbone_op_matrix_euler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  
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
            elseif chooseCongener == 10 
                OPT_chem_DDT;
                congener = 'DDT';      
            end
      
    figure1 = figure(100); %clf
    figure_title = [congener ': 5-year-mean CSD from 1994 until 2009'];        
    annotation(figure1,'textbox', [0.35 0.97 0.82 0.03], 'String',{figure_title}, 'FontWeight','bold', ...
        'FontSize',14, 'FitBoxToText','off', 'LineStyle','none');

    set(gcf,'Color',[1,1,1]); 
    set(gcf, 'PaperPositionMode','auto')     %important! don't know why... :-)
    set(0,'DefaultAxesFontName', 'Times New Roman'); set(0,'DefaultTextFontname', 'Times New Roman')
    set(0,'DefaultTextFontsize', 10); set(0,'DefaultAxesFontsize', 10)
    set(gcf, 'Units','centimeters','position', [-45 0 23 17])

    for i = 1:length(emp_year)
        subplot(4,4,i)
        %semilogy(age_group',emp_data_noNANs(:,i),'o'); hold on     %% empirical DATA
        semilogy(age_group', mod_data_nonNaNs(:,i),'or'); hold on   %% modeled DATA
        semilogy((1/12:1:79+1/12)', mod_data_nonNaNs_all(:,i),'r');
        xlim([0 80]); 
        
         if chooseCongener == 1     
                ylim([0.1 100]);  
            elseif chooseCongener == 2 
                ylim([0.1 100]);  
            elseif chooseCongener == 3 
                ylim([0.1 100]);  
            elseif chooseCongener == 4 
                ylim([0.1 1000]); 
            elseif chooseCongener == 5 
                ylim([0.1 1000]); 
            elseif chooseCongener == 6 
                ylim([0.1 1000]);  
            elseif chooseCongener == 7
                ylim([0.1 1000]);  
            elseif chooseCongener == 8 
                ylim([0.1 1000]);
            elseif chooseCongener == 9 
                ylim([0.1 1000]);
            elseif chooseCongener == 10 
               ylim([0.1 1000]);   
        end       

        title(num2str(emp_year(i)))
    end   
   
%% figure of choice of year
year_choice = 1978; 
year_choice_m = (year_choice-Y_begin)*12;

    figure2 = figure(2); %clf
    set(gcf,'Color',[1,1,1]); 
    set(gcf, 'PaperPositionMode','auto')     %important! don't know why... :-)
    set(0,'DefaultAxesFontName', 'Times New Roman'); set(0,'DefaultTextFontname', 'Times New Roman')
    set(0,'DefaultTextFontsize', 10); set(0,'DefaultAxesFontsize', 10)
    set(gcf, 'Units','centimeters','position', [-50 0 23 17])

    plot((1/12:1:57+1/12)', flipud(Conc_monthly_all_0to80(101:158, year_choice_m)), 'r')
    title(['CSD in year ' num2str(year_choice)])

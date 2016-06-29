%% plot intake

fabs = 0.9;
intake = Iref*1000/D/fabs; % ng/g/month --> ng/kg/d, this is true intake, not uptake anymore

%% plots 

    figure(chooseCongener+20)
        set(gcf,'Color',[1,1,1]); hold on 
        set(gcf, 'PaperPositionMode','auto')     %important! don't know why... :-)
        set(0,'DefaultAxesFontName', 'Arial'); set(0,'DefaultTextFontname', 'Arial')
        set(0,'DefaultTextFontsize', 10); set(0,'DefaultAxesFontsize', 10)
        set(gcf, 'Units','centimeters','position', [-17 5 15 10])

        plot(time, intake, 'r', 'LineWidth', 2)
        % plot(time, intake./max(intake(:)), 'Color', col)
        ylabel('intake in ng/kg/d')
        xlabel('time')

        disp(['Intake in peak year ' num2str(Y_peak) ' is ' num2str(max(intake)) ' ng/kg/d' ]) 
        
        Y_chosen = 2000; 
        intake_y_chosen = max(intake)*exp(kdec_down*(Y_peak-Y_chosen)); % ng/kg/d
        
        disp(['Intake in year ' num2str(Y_chosen) ' is ' num2str(intake_y_chosen) ' ng/kg/d' ]) 
   
        if strcmp(chooseOpt,'intake')
            figure(chooseCongener+20); hold on 
            plot(int_emp_year, int_emp*1000/D, '*k') % ng/g/months --> ng/kg/d
        
        elseif strcmp(chooseOpt,'original')
            figure(chooseCongener+20); hold on 
            plot(int_emp_year, int_emp*1000/D, '*b') % ng/g/months --> ng/kg/d
            
        end 

        
 
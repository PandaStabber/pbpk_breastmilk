%% Plot cross-sectional data - in combination with PK_model_basis.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

figure1 = figure(1); clf
figure_title = ['modeled age groups from 1994 until 2009'];        
annotation(figure1,'textbox', [0.35 0.97 0.82 0.03], 'String',{figure_title}, 'FontWeight','bold', ...
    'FontSize',14, 'FitBoxToText','off', 'LineStyle','none');

set(gcf,'Color',[1,1,1]); 
set(gcf, 'PaperPositionMode','auto')     %important! don't know why... :-)
set(0,'DefaultAxesFontName', 'Times New Roman'); set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontsize', 10); set(0,'DefaultAxesFontsize', 10)
set(gcf, 'Units','centimeters','position', [-50 0 23 17])

for i = 1994:2009
    subplot(4,4,i-1994+1)
    semilogy(age_group', mod_data_nonNaNs(:,i-1994+1),'or');
    xlim([0 80]); ylim([0.1 100]); 
    title(num2str(i))
end   

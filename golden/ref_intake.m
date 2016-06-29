function [Iref] = ref_intake(Ipeak, M, Y_begin, Y_peak, year, kdec_up, kdec_down)
% REF_INTAKE Calculation of symmetric exponential increase and decrease of adult reference intake [ng/gbw/month]
%   First-order decrease and increase rate constants are fixed.  .
%   kdec_up, kdec_dwon are in units of 1/yrs

% Iref = zeros(length(year),1);      % [ng/g/month] adult refernce intake for all year
% Iref(1:(Y_peak-Y_begin)*M,1)       = Fit_Var(1)*exp(kdec.*(year(1:(Y_peak-Y_begin)*M)-Y_peak));               
% Iref((1+(Y_peak-Y_begin)*M):end,1) = Fit_Var(1)*exp(-kdec.*(year((Y_peak-Y_begin)*M+1:length(year))-Y_peak)); 

Iref = zeros(1, length(year));      % [ng/g/month] adult refernce intake for all year
Iref(1,1:(Y_peak-Y_begin)*M)        = Ipeak.*exp(kdec_up.*(year(1:(Y_peak-Y_begin)*M)-Y_peak));                 
Iref(1,(1+(Y_peak-Y_begin)*M):end)  = Ipeak.*exp(-kdec_down.*(year((Y_peak-Y_begin)*M+1:length(year))-Y_peak)); 
 
% Iref = ones(1,length(year)) * 5;   % [ng/gbw/month] adult refernce intake for all year
%                                               % constant Iref as test/comparison  

end


%% adult reference intake with plateau
% exponential increase and decrease of adult refernce intake [ng/gbw/week] with peak intake in 2004
% function [Iref] = ref_intake(Fit_Var, M, Y_begin, Y_peak, year, kdec)
% yr_pl = 5; % number of year during which plateau intake
% 
% Iref = zeros(1,length(year)); % [ng/gbw/week] adult refernce intake for all year
% Iref(1,1:(Y_peak-Y_begin)*M)                                  = Fit_Var*exp(kdec.*(year(1:(Y_peak-Y_begin)*M)-Y_peak));      
% Iref(1,(1+(Y_peak-Y_begin)*M):(1+(Y_peak-Y_begin)*M)+yr_pl*M) = Fit_Var;
% Iref(1,(1+(Y_peak-Y_begin)*M)+yr_pl*M:end)                    = Fit_Var*exp(-kdec.*(year((1+(Y_peak-Y_begin)*M)+yr_pl*M:end)-Y_peak-yr_pl));  
% 
% end

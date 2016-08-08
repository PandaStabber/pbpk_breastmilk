%% Chemical specific parameters -- HCB

% time points
    Y_begin = 1921;                         % [year] first year considered in model
    Y_end   = Y_begin+people/2+age_max_y-1; % [year] last year considered in model
    Y_peak  = 1977;                         % [year] reference year, i.e. peak exposure in this year according to 03_exposureCZ(premonitoring).xlsx
    time    = (Y_begin+1/M:1/M:Y_end)';     % [year] times between Y_begin and Y_end in monthly steps

% Different fitting options
if strcmp(chooseOpt,'intake')
    Th_elim = 10*M;            % [months] intrinsic elimination half-life in months, Ritter et al. 2011
    kelim   = log(2)/Th_elim;   % [1/months] intrinsic elimination rate constamt in months
    % from CSTD_input_data.xlsx --> Ruprich et al. (2009)
    Th_dec_down  = 5;                     % [years] doubling time and half-life of intake exposure in YEARS (unit is correct) from Cupr air monitoring
    Th_dec_up    = Th_dec_down;             % [years]
    kdec_up      = log(2)/Th_dec_up;        % [1/yr]
    kdec_down    = log(2)/Th_dec_down;      % [1/yr]
    % known intake
    int_data      = xlsread('CSTD_input_data.xlsx', 'Mikes et al 2012', 'B25:C35'); % 

    int_emp_year  = int_data(:,1); 
    int_emp       = int_data(:,2)*1000;  % ug/kg/d --> ng/kg/d
    int_emp       = int_emp/1000*D;      % ng/kg/d --> ng/g/months
    
elseif strcmp(chooseOpt,'int_kelim')        
    % from CSTD_input_data.xlsx --> Ruprich et al. (2009)
    Th_dec_down  = 4.3;                     % [years] doubling time and half-life of intake exposure in years from Cupr air monitoring
    Th_dec_up    = Th_dec_down;             % [years]
    kdec_up      = log(2)/Th_dec_up;        % [1/yr]
    kdec_down    = log(2)/Th_dec_down;      % [1/yr]

elseif strcmp(chooseOpt,'original')        
    % known intake
    int_data      = xlsread('CSTD_input_data.xlsx', 'Mikes et al 2012', 'B25:C35'); % 
    int_emp_year  = int_data(:,1); 
    int_emp       = int_data(:,2)*1000;  % ug/kg/d --> ng/kg/d
    int_emp       = int_emp/1000*D;      % ng/kg/d --> ng/g/months

elseif strcmp(chooseOpt,'telim')        
    % from CSTD_input_data.xlsx --> Ruprich et al. (2009)
    Th_dec_down  = 4.3;                     % [years] doubling time and half-life of intake exposure in years from Cupr air monitoring
    Th_dec_up    = Th_dec_down;             % [years]
    kdec_up      = log(2)/Th_dec_up;        % [1/yr]
    kdec_down    = log(2)/Th_dec_down;      % [1/yr]
    % known intake
    int_data      = xlsread('CSTD_input_data.xlsx', 'Mikes et al 2012', 'B25:C35'); % 
    int_emp_year  = int_data(:,1); 
    int_emp       = int_data(:,2)*1000;  % ug/kg/d --> ng/kg/d
    int_emp       = int_emp/1000*D;      % ng/g/months
    
end

%% Alternative kdec
    % exponetial increase and decrease of adult reference intake, according to
    % Breivik_2007 max emission scenario: http://www.nilu.no/projects/globalpcb/globalpcb2.htm
    %Th_dec_up    = 4.3;                     % [years] doubling time and half-life of intake exposure in years
    %Th_dec_down  = 4.3;                     % [years] from 01_Dietary consumption.xlsx
   
%% Empirical data - data_all_mean generated in load_emp_data_primi.m
    % concentration 
    emp_data      = data_all_mean(8*no_ag+1:9*no_ag,2:end);   % extration of emp data from full matrix
    emp_year      = [1994:1:2003,2005:1:2009];              % years with data 

    emp_data(emp_data == 0) = NaN; % convertion of LoD to NANs --> columns with ONLY NANs will be deleted in control_O.m

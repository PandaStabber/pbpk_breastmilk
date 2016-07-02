%% Forward model for interpretation of POPs in breast milk for Czech mothers from 1994-2009

clc; clear; close all; format long; format compact; tic

% profile on 

% Time variables (single-value parameters --> input_single matrix)
    M   = 12;                        % [months/year] unit conversion factor
    D   = 30;                        % [days/month]  unit conversion factor
    people    = 2*100;               % total number of individuals, females only, ***230***should be the max number
    age_max_y = 80;                  % [years] maximum age of a person in years
    age_max_m = age_max_y * M;       % [months] maximum age of a person in months
    age_mom_y = 25;                  % [years] age of mother at birth in years 
    age_mom_m = age_mom_y * M;       % [months] age of mother at birth in months
    age_y     = (1:age_max_y);       % [years] age in years with interval of 1 year
    age_m     = (1:age_max_m);       % [months] age in weeks with interval of 1 month
    age_ym    = (1/M:1/M:age_max_y); % [years] age in years with interval of 1 month
    age_group     = 15:5:45;

% empirical data 
    load data_all_mean.mat
    emp_year_all = 1994:2009;    
   
   
%% change variables accordingly
    fabs    = 0.9; % [-], absorption factor, generic value
    I_peak  = 80*fabs/1000*D; % ng/kg/d --> ng/g/month
    no_ag   = 7;    % number of age groups
                            % 6:PCB180, 7:PCB170, 8:HCB, 9:DDE, 10:DDT
    % Choose optimization 
        chooseOpt = 'intake';
        %chooseOpt = 'int_kelim';
        %chooseOpt = 'original';     
    OPT_chem_HCB
    chooseCongener = 8;     % 1:PCB28, 2:PCB52, 3:PCB118, 4:PCB138, 5:PCB153, 
      
%% Number of runs
    nRuns = 1;    
    
% Number of congeners    
%     nCongeners = 1; 

% Choose body weight 
    bodyweight_cz()
    
% P value
    load P_ritter; P = P_ritter(1:960); Psource = 'Pritter'; % is my default choice, see page 13 in book *cstd*, also CSTD_Pritter_Plorber_check.fig
    
% breastfeeding 
    t_lac   = 6;                % [months] 6 months exclusively --> generic value
    bm_flip = 0.0034 * log(age_ym(1:t_lac)) + 0.0414; % Verner (2013), EHP
    bm_m(1:t_lac)           = (-0.0024 .* age_ym(1,1:t_lac) + 0.0063) .* Mbw(1,1:t_lac)*24*D; 
%     disp(bm_m)
    % bm_m(2:t_lac)         = 0; 
    klac_in(1:t_lac,1)      = bm_m.*bm_flip./Mlip(1,1:t_lac); 
    klac_out = zeros(1,age_max_m); 
    klac_out(1,age_mom_m+1:age_mom_m+t_lac) = bm_m.*bm_flip./Mlip(1,age_mom_m+1:age_mom_m+t_lac); 
       
%% acutal calculation and optimization 

%for chooseCongener = 1:nCongeners
    
    %% prepartatin of input parameters 
    Uptake_monthly_all  = zeros(people, numel(time), nRuns);
    Conc_monthly_all_0to80  = zeros(people, numel(time), nRuns);
    Conc_monthly      = zeros(people, numel(time));  % [ng/glip] concentration
    Mbw_monthly       = Conc_monthly;  % [g] body weight 
    Mlip_monthly      = Conc_monthly;  % [glip] lipid weight 
    P_monthly         = Conc_monthly;  % [-] proportionality factor
    Flux_monthly      = Conc_monthly;  % [ng] monthly transfer of mass  
    Flux_elim_monthly = Conc_monthly;  % [ng/glip] monthly transfer of conc  
    Flux_lac_monthly  = Conc_monthly;  % [ng/glip] monthly transfer of conc  
    Flux_gro_monthly  = Conc_monthly;  % [ng/glip] monthly transfer of conc Kelim_monthly  = zeros(people, people, nume(time));  % [1/month] elimination rate constant
    Kelim_monthly     = zeros(people, people, numel(time));  % [1/month] elimination rate constant
    Klac_monthly      = Kelim_monthly;  % [1/month] lactation rate constant
    Kgro_monthly      = Kelim_monthly;  % [1/month] growth dilution rate constant

    for a = 1:people/2 % FEMALES who give birth
        % females 1
        Mbw_monthly  (a,(a-1)*M+1:age_max_y*M+(a-1)*M)   = Mbw;
        Mlip_monthly (a,(a-1)*M+1:age_max_y*M+(a-1)*M)   = Mlip; 
        P_monthly    (a,(a-1)*M+1:age_max_y*M+(a-1)*M)   = P; 
        Kgro_monthly (a,a,(a-1)*M+1:age_max_y*M+(a-1)*M) = -kgrowth;  % !! negativ sign !! 
        Klac_monthly (a,a,(a-1)*M+1:age_max_y*M+(a-1)*M) = -klac_out; % !! negativ sign !!
       
        if a > age_mom_y
            Klac_monthly(a,a-age_mom_y,age_mom_m+1+(a-age_mom_y-1)*M:age_mom_m+(a-age_mom_y-1)*M+t_lac) = klac_in; 
        end       
    end
    
    for a = people/2+1:people % FEMALES who did not give birth and are compared to emp. data
        % females 2
        Mbw_monthly  (a,(a-people/2-1)*M+1:age_max_y*M+(a-people/2-1)*M)   = Mbw;
        Mlip_monthly (a,(a-people/2-1)*M+1:age_max_y*M+(a-people/2-1)*M)   = Mlip;
        P_monthly    (a,(a-people/2-1)*M+1:age_max_y*M+(a-people/2-1)*M)   = P;
        Kgro_monthly (a,a-people/2,(a-people/2-1)*M+1:age_max_y*M+(a-people/2-1)*M) = -kgrowth;  % !! negativ sign !! 
        
        if a > people/2 + age_mom_y
            Klac_monthly(a,a-people/2-age_mom_y,age_mom_m+1+(a-people/2-age_mom_y-1)*M:age_mom_m+(a-people/2-age_mom_y-1)*M+t_lac) = klac_in;
        end
    end
    
    %% runs    
    for run = 1:nRuns
                             
            %age-dependent kelim
            %t_met = Th_elim .*(Mlip./Mlip_ref).*(Mliv_ref./Mliv).^(2/3); 
            %kelim = log(2)./(t_met*M); % 1/months 
            
       % prepartatin of input parameters - ONLY KELIM BECAUSE THIS IS CONGENER-SPECIFIC
            for a = 1:people/2 
                Kelim_monthly(a,a,(a-1)*M+1:age_max_y*M+(a-1)*M)          = -kelim';    % !! negativ sign !! females 1
                Kelim_monthly(a+people/2,a,(a-1)*M+1:age_max_y*M+(a-1)*M) = -kelim';    % !! negativ sign !! females 2
            end
      
            Ktot_monthly = Kelim_monthly + Klac_monthly + Kgro_monthly; 
        
            % ref_intake        
            [Iref] = ref_intake(I_peak, M, Y_begin, Y_peak, time, kdec_up, kdec_down);
            
                % Iref(1440:end) = 0; % postban no intake anymore
            
            Iref_monthly   = zeros(people, numel(time));  % [ng/g/month] monthly adult reference intake 

            % multi-individual intake
            % females 1
            for a = 1:people/2
                Iref_monthly (a,(a-1)*M+1:age_max_y*M+(a-1)*M)   = Iref(1,(a-1)*M+1:age_max_y*M+(a-1)*M);
            end

            % females 2
            Iref_monthly(people/2+1:end,:) = Iref_monthly(1:people/2,:); % females 2 copy of females 1

            Intake_monthly = Iref_monthly.*P_monthly.*Mbw_monthly./Mlip_monthly; % [ng/glip/month]
            Intake_monthly(isnan(Intake_monthly)) = 0;
            
            
            %% actual calculation
            % euler - matrix form
            no_timestep   = 10;                             % number of timesteps per month
            dt            = 1/no_timestep;                  % fraction of timestep per month
            conc          = zeros(people, no_timestep); 
            flux          = conc; 
            flux_elim     = conc;
            flux_lac      = conc;
            flux_gro      = conc;
            conc_initials = zeros(people, numel(time));
            conc_ini      = zeros(people,1);

            count = age_mom_y+1; 
            birthlist = (age_mom_m:M:age_mom_m+(people/2-count)*M)'; % months at which the conc of mother is taken as initial conc for baby

            for month = 1:numel(time); 
                I     = Intake_monthly(:,month); 
                K     = Ktot_monthly(:,:,month);
                Kelim = Kelim_monthly(:,:,month);
                Klac  = Klac_monthly(:,:,month); 
                Kgro  = Kgro_monthly(:,:,month);
                
                conc_initials(:,month) = conc_ini; 

                for timestep = 1:no_timestep

                    if timestep == 1      
                        conc(:,timestep) = conc_ini + I * dt + K * conc_ini * dt; 
                        flux(:,timestep) = -(K * conc_ini * dt); 
                        flux_elim(:,timestep) = -(Kelim * conc_ini * dt); 
                        flux_lac(:,timestep)  = -(Klac * conc_ini * dt); 
                        flux_gro(:,timestep)  = -(Kgro * conc_ini * dt);  
                    else
                        conc(:,timestep) = conc(:,timestep-1) + I * dt + K * conc(:,timestep-1) * dt;
                        flux(:,timestep) = flux(:,timestep-1) + -(K * conc(:,timestep-1) * dt);
                        flux_elim(:,timestep) = flux_elim(:,timestep-1) + -(Kelim * conc(:,timestep-1) * dt);
                        flux_lac(:,timestep) = flux_lac(:,timestep-1) + -(Klac * conc(:,timestep-1) * dt);
                        flux_gro(:,timestep) = flux_gro(:,timestep-1) + -(Kgro * conc(:,timestep-1) * dt);
                    end

                end

                Conc_monthly(:,month) = conc(:,no_timestep); 
                Flux_monthly(:,month) = flux(:,no_timestep); 
                Flux_elim_monthly(:,month) = flux_elim(:,no_timestep);
                Flux_lac_monthly(:,month) = flux_lac(:,no_timestep);
                Flux_gro_monthly(:,month) = flux_gro(:,no_timestep);

                conc_ini = conc(:,no_timestep);  

                if ~isempty(find(birthlist(:,1) == month, 1)) == 1 && people > age_mom_y
                    conc_ini(count) = Conc_monthly(count-age_mom_y,age_mom_m+(find(birthlist(:,1) == month) - 1)*12); 
                    conc_ini(count+people/2) = Conc_monthly(count-age_mom_y,age_mom_m+(find(birthlist(:,1) == month) - 1)*12); 
%                     conc_ini(count) = 0; 
                    % explanation: (find(birthlist(:,1) == month) - 1) gives me the row number of the desired month in birthlist
                    count = count+1;
                end

            end
                     
            % Conc_monthly 
                Conc_monthly(Conc_monthly == 0) = NaN; % vor der Geburt

                for a = 1:people/2 % Females 1
                    Conc_monthly(a, age_max_m+1 + (a-1)*M : end) = NaN; % nach dem Tod 
                end 
         
                for a = people/2+1:people % Females 2 
                    Conc_monthly(a, age_max_m+1 + (a-people/2-1)*M : end) = NaN; % nach dem Tod 
                end 
                
                Conc_monthly_rerun = Conc_monthly; % insert the initial concentration for plotting age 0-80 and not 1/12 to 80 1/12
                
                for b = 1:people/2-age_mom_y
                    Conc_monthly_rerun(age_mom_y+b, age_mom_m + (b-1)*M)          = conc_initials(age_mom_y+b, age_mom_m + (b-1)*M+1); 
                    Conc_monthly_rerun(age_mom_y+b+people/2, age_mom_m + (b-1)*M) = conc_initials(age_mom_y+b+people/2, age_mom_m+ (b-1)*M+1); 
                end
                Conc_monthly_all_0to80(:,:,run) = Conc_monthly_rerun;  
                
  
               
                
        %% determine fluxes - check mass balance 
            % version 1
            Kdiag_loss_elim = zeros(people,numel(time));  % loss via elimination (minus rate constant)
            Kdiag_loss_lac  = Kdiag_loss_elim;            % loss via lactation (minus rate constant)
            Kdiag_in_lac_f  = zeros(people-age_mom_y,numel(time));  % intake via lactation (plus rate constant)
            Kdiag_loss_gro  = Kdiag_loss_elim;            % loss via growth (minus rate constant)

            Flux_loss_elim  = Kdiag_loss_elim;            % mass lost via elimination [ng/month]
            Flux_loss_lac   = Kdiag_loss_elim;            % mass lost via lactation [ng/month]
            Flux_in_lac     = Kdiag_loss_elim;            % mass intake via lactation [ng/month]
            Flux_loss_gro   = Kdiag_loss_elim;            % conc loss via growth [ng/glip/month]

            for ii = 1:numel(time)
                Kdiag_loss_elim(:,ii)  = diag(Kelim_monthly(:,:,ii));
                Kdiag_loss_lac(:,ii)   = diag(Klac_monthly(:,:,ii)); %Kdiag_loss_lac = 0 because klac_out = 0
                Kdiag_in_lac_f(:,ii)   = diag(Klac_monthly(:,:,ii), -(age_mom_y)); %-age_mom_y steht für Verschiebung in Diagonale
                Kdiag_loss_gro(:,ii)   = diag(Kgro_monthly(:,:,ii)); 
            end

            Kdiag_in_lac_fem  = [zeros(age_mom_y, numel(time)); Kdiag_in_lac_f];
            Kdiag_sum         = Kdiag_loss_elim + Kdiag_loss_lac + Kdiag_loss_gro + Kdiag_in_lac_fem;

            for ii = 1:people
                Flux_loss_elim(ii,:)    = Flux_monthly(ii,:).*Kdiag_loss_elim(ii,:)./Kdiag_sum(ii,:);   
                Flux_loss_lac(ii,:)     = Flux_monthly(ii,:).*Kdiag_loss_lac(ii,:)./Kdiag_sum(ii,:);     
                Flux_in_lac_fem (ii,:)  = Flux_monthly(ii,:).*Kdiag_in_lac_fem(ii,:)./Kdiag_sum(ii,:);    
                Flux_loss_gro(ii,:)     = Flux_monthly(ii,:).*Kdiag_loss_gro(ii,:)./Kdiag_sum(ii,:);    
            end

            Ball_fem  = sum(Intake_monthly(1,1:age_max_m)) - sum(Flux_monthly(1, 1:age_max_m)) - Conc_monthly(1,age_max_m)

            % version 2
            Bal4_fems  = zeros(people/2,1);

            for i = 1:people/2  % FEMALES ONLY 
                Bal4_fems(i,1) = sum(Intake_monthly   (i,M*(i-1)+1:age_max_m+M*(i-1))) ...    % in  via food, dust, etc
                               - sum(Flux_elim_monthly(i,M*(i-1)+1:age_max_m+M*(i-1))) ...    % minus flux_tot = 1. minus elimination
                               - sum(Flux_gro_monthly (i,M*(i-1)+1:age_max_m+M*(i-1))) ...    %                  2. minus growth 
                               - sum(Flux_lac_monthly (i,M*(i-1)+1:age_max_m+M*(i-1))) ...    %                  3. minus lactatin
                               - Conc_monthly(i, age_max_m+M*(i-1)) ...                       % left at the end 
                               + conc_initials(i,age_mom_m+1+M*(i-age_mom_y-1));              % in  at birth
            end
               
            Bal4_fems_max = max(Bal4_fems)
            disp('mass balance error')
            

        %% extraction of CSD sets from Conc_monthly for FEMALES 2
        CSD_all = zeros(length(emp_year_all),age_max_y); 
        mo_year_all = (emp_year_all-Y_begin)*M;

        for i = 1:length(emp_year_all)
            CSD_all(i,:) = Conc_monthly_all_0to80(emp_year_all(i)-age_max_y-Y_begin+2+people/2:emp_year_all(i)-Y_begin+1+people/2, mo_year_all(i))';  % CSD concentration of age 0 - 80 years
        end 
        CSD_all = fliplr(CSD_all)'; 
        mod_data_nonNaNs_all = CSD_all; 
        mod_data_nonNaNs = mod_data_nonNaNs_all((age_group+1)', :);
        
    end % end run

    Mass_monthly_all_0to80 = Conc_monthly_all_0to80.*Mlip_monthly; 
    
%     figure_name = ['141001_oi_' scenario '_' num2str(Y_peak) '_' num2str(round(Th_dec_up)) 'yr_' Psource '_' congener]; 
%     saveas(figure(chooseCongener), figure_name); 
    
%end % end chooseCongener

%% plots

%
% % graphs_CSD
% graphs_CSTD
%  graphs_intake
%  graphs_CSD_congener

% graphs_LD


% profile viewer
toc




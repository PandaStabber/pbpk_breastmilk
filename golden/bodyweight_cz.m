%% Bodyweight 


%% body weight .mat files
load bw_cz.mat      % [kg]
load flip_cz.mat    % [%]
flip = (flip/100)'; % [-]

Mbw = (Mbw.*1000)';           % [g] 
Mbw_ini  = 3.1744*1000;       % [g]
flip_ini = 0.215;             % [-]
Mlip_ini = Mbw_ini*flip_ini;  % [glip]
Mlip = Mbw.*flip;             % [glip]

kgrowth = zeros(1,age_max_m); % rate constant of growth, [1/month]
kgrowth(1,1) = (Mlip(1,1)-Mlip_ini)/Mlip(1,1);

for a = 2:age_max_m
    kgrowth (1,a) = (Mlip(1,a)-Mlip(1,a-1))/Mlip(1,a); % [1/month] 
	
end


Mbw_av  = 69; % [kg]



%% Calculation of body weight from 03_Measures_CZ.xlsx obtained from Ondrej Mikes
% 
% bw_all = xlsread('03_Measures_CZ.xlsx','final_upload','A2:B67');
% flip_all1 = xlsread('03_Measures_CZ.xlsx','final_upload','C2:E48');
% flip_all2 = xlsread('03_Measures_CZ.xlsx','final_upload','C60:E66');
% 
% flip_all = [flip_all1(:,1) flip_all1(:,3); 
%             flip_all2(:,1) flip_all2(:,3)];
% figure(1)
% plot(bw_all(:,1), bw_all(:,2),'*')
% 
% figure(2)
% plot(flip_all(:,1), flip_all(:,2),'*')



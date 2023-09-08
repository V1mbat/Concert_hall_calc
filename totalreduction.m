clear 
close all
clc
%% new window values [35 53 67 60 52 53 66]
%%
L_facade = [71 73 75 80 66 54 50]; % SPL at facade
NC_15 = [47 36 29 22 17 14 12];
NC_20 = [51 40	33 26 22 19	17];   
NC_40 = [64 56 50 45 41 39 38];

f_c = 10^3 * (2.^((-4:2)));

c_0 = 340;                      % speed of sound
rho_0 = 1.2;                    % density air 
%% double wall
% you can adjust all those parameters and put in material data you found
% the format of the data should be like below

% wall 1
wall1.rho = 900;                % density wall 1
wall1.t = 0.15;                 % thickness wall 1
wall1.m = wall1.rho * wall1.t;  % mass per unit area wall 1
wall1.eta = 0.025;              % loss factor wall 1
wall1.E = 10E9;                 % young's modulus wall 1
wall1.nu = 0.3;                 % poisson ratio wall 1

% wall 2
wall2.rho = 900;                % density wall 2
wall2.t = 0.1;                  % thickness wall 2
wall2.m = wall2.rho * wall2.t;  % mass per unit area wall 2
wall2.eta = 0.015;              % loss factor wall 2
wall2.E = 10E9;                 % young's modulus wall 2
wall2.nu = 0.3;                 % poisson ratio wall 2

% glass 1
glass1.rho = 2500;               % density wall 3
glass1.t = 0.0042;               % thickness wall 3
glass1.m = glass1.rho * glass1.t;% mass per unit area wall 3
glass1.eta = 1E-4;               % loss factor wall 3
glass1.E = 80E9;                 % young's modulus wall 3
glass1.nu = 0.3;                 % poisson ratio wall 3

% glass 2
glass2.rho = 2500;               % density wall 1
glass2.t = 0.015;                % thickness wall 1
glass2.m = glass2.rho * glass2.t;% mass per unit area wall 1
glass2.eta = 1E-4;               % loss factor wall 3
glass2.E = 80E9;                 % young's modulus wall 3
glass2.nu = 0.3;                 % poisson ratio wall 3

% glass 3
glass3.rho = 2500;               % density wall 1
glass3.t = 0.0035;               % thickness wall 1
glass3.m = glass3.rho * glass3.t;% mass per unit area wall 1
glass3.eta = 1E-4;               % loss factor wall 3
glass3.E = 80E9;                 % young's modulus wall 3
glass3.nu = 0.3;                 % poisson ratio wall 3

% glass 4
glass4.rho = 2500;               % density wall 1
glass4.t = 0.012;                % thickness wall 1
glass4.m = glass4.rho * glass4.t;% mass per unit area wall 1
glass4.eta = 1E-4;               % loss factor wall 3
glass4.E = 80E9;                 % young's modulus wall 3
glass4.nu = 0.3;                 % poisson ratio wall 3


% distance between the walls
d_1 = 0.25;                     % distance between wall 1 and wall 2
d1_glass = 0.004;               % distance between the glass panels
d2_glass = 0.006;
d_34 = 0.25;                    % distance between wall 3 and wall 4

f = 1:10000;                    % frequency domain

doubleglass1.rho = 2500*(glass1.t + glass2.t)/(glass1.t + glass2.t + d1_glass);
doubleglass1.m = glass1.m + glass2.m;
doubleglass1.t = (glass1.t + glass2.t + d1_glass);  % mass per unit area wall 1
doubleglass1.eta = 1E-4;               % loss factor wall 3
doubleglass1.E = 80E9;                 % young's modulus wall 3
doubleglass1.nu = 0.3;                 % poisson ratio wall 3

doubleglass2.rho = 2500*(glass3.t + glass4.t)/(glass3.t + glass4.t + d2_glass);
doubleglass2.m = glass3.m + glass4.m;
doubleglass2.t = (glass3.t + glass4.t + d2_glass);  % mass per unit area wall 1
doubleglass2.eta = 1E-4;               % loss factor wall 3
doubleglass2.E = 80E9;                 % young's modulus wall 3
doubleglass2.nu = 0.3;                 % poisson ratio wall 3

%% reduction indexi of doublewalls and each wall 
r12_double = mean_oct(f_c, doublewall(f,wall1,wall2,d_1,false));  % reduction index for doublewall (wall 1 and wall 2)
r_glass1 = doublewall(f,glass1,glass2,d1_glass,false);  % reduction index for doublewall (wall 3 and wall 4) 
r_glass2 = doublewall(f,glass3,glass4,d2_glass,false);

doubleglass1.r = r_glass1;
doubleglass2.r = r_glass2;

r_doubleglass = mean_oct(f_c, doubleglass(f,doubleglass1,doubleglass2,d_34));

r1_single = mean_oct(f_c, singlewall(f,wall1));
r2_single = mean_oct(f_c, singlewall(f,wall2));
r3_single = singlewall(f,glass1);
r4_single = singlewall(f,glass2);

%% for having multiple wall sections:
% this calculates the overall reduction index of the wall, when part of the
% wall is a double wall consisting of wall 1 and wall 2 and part of it is a
% double wall consisting of wall 3 and wall 4
% adjust S_12 and S_34 for the areas of the partitions

S_12 = 5850;        % area of the double wall made out of wall 1 and wall 2
S_34 = 117.6;         % area of the double wall made out of wall 3 and wall 4

r_glass = r_doubleglass;
r_comb = 10*log10((S_12+S_34)./(S_12*10.^(-r12_double/10) + S_34*10.^(-r_glass/10)));
r_test = r1_single + r2_single;
%% room parameters

S_hall = S_12 + S_34;           % area of the outer wall of the hall in m^2
V_hall = 21727;                 % volume of hall
T_hall = 2.0;                   % reverb time hall
A_hall = 0.163*V_hall/T_hall;   % equivilant absorption area in hall

L_hall = L_facade + 3 - r_comb + 10*log10(S_hall/A_hall);

S_lobby = S_12 + S_34;           % area of the outer wall of the hall in m^2
V_lobby = 18070;                 % volume of hall
T_lobby = 2.0;                   % reverb time hall
A_lobby = 0.163*V_lobby/T_lobby; % equivilant absorption area in hall

L_lobby = L_facade + 3 - r1_single + 10*log10(S_lobby/A_lobby);

%% plotting
figure
semilogx(f_c, L_hall)
hold on
grid on 
semilogx(f_c, NC_15)
semilogx(f, r3_single)
semilogx(f, r4_single)
thickenall_big
xlabel('f in Hz')
ylabel('L_{SPL} in dB')
legend('L_{hall}','NC_{20}', 'r glass 1', 'r glass 2','NumColumns',2)
xlim([50 5000])
ylim([-20 50])
legend(Location='southoutside')

figure
semilogx(f_c, r_comb)
legend('combined reduction index')
grid on
thickenall_big

figure
semilogx(f_c, r_doubleglass)
legend('Reduction Index Window')
ylim([20 100])
grid on
thickenall_big
clear 
close all
clc

%%
L_facade = [73 75 80 66 54 50]; % SPL at facade
NC_20 = [40	33 26 22 19	17];    % noise requirement in concert hall

f_c = 10^3 * (2.^((-3:2)));
%% double wall

wall1.rho = 700;                % density wall 1
wall1.t = 0.3;                  % thickness wall 1
wall1.m = wall1.rho * wall1.t;  % mass per unit area wall 1
wall1.eta = 0.015;              % loss factor wall 1
wall1.E = 10E9;                 % young's modulus wall 1
wall1.nu = 0.3;                 % poisson ratio wall 1

wall2.rho = 900;                % density wall 2
wall2.t = 0.3;                  % thickness wall 2
wall2.m = wall2.rho * wall2.t;  % mass per unit area wall 2
wall2.eta = 0.015;              % loss factor wall 2
wall2.E = 10E9;                 % young's modulus wall 2
wall2.nu = 0.3;                 % poisson ratio wall 2

d_z = 0.05;                     % distance between walls
c_0 = 340;                      % speed of sound
rho_0 = 1.2;                    % density air 

f = 1:10000;                    % frequency domain

r_doublewall = doublewall(f,wall1,wall2,d_z);
% 
reduction_index2 = singlewall(f,wall1);
% 
reduction_index3 = singlewall(f,wall2);

%% lobby
red_oct = mean_oct(f_c, r_doublewall);
red2_oct = mean_oct(f_c, reduction_index2);

S_lobby = 1100;                 % area of the outer wall of the lobby in m^2
V_lobby = 12070;                % volume of lobby
T_lobby = 1.0;                  % reverb time lobby
A_lobby = 0.163*V_lobby/T_lobby;% equivilant absorption area in lobby

S_hall = 1770;                  % area of the outer wall of the lobby in m^2
V_hall = 18070;                 % volume of lobby
T_hall = 2.0;                   % reverb time lobby
A_hall = 0.163*V_hall/T_hall;   % equivilant absorption area in lobby

L_lobby = L_facade + 3 - red_oct + 10*log10(S_lobby/A_lobby);
L_hall = L_facade + 3 - red2_oct + 10*log10(S_hall/A_hall);

%% plotting
figure
semilogx(f_c, L_lobby)
hold on
grid on 
semilogx(f_c, NC_20)
thickenall_big
xlabel('f in Hz')
ylabel('L_{SPL} in dB')
legend('L_{Lobby}','NC_{20}')

figure
semilogx(f, r_doublewall)
hold on
semilogx(f, reduction_index2)
semilogx(f, reduction_index3)
legend('double wall', 'wall 1', 'wall 2')
thickenall_big

figure
semilogx(f_c, L_hall)
hold on
grid on 
semilogx(f_c, NC_20)
thickenall_big
xlabel('f in Hz')
ylabel('L_{SPL} in dB')
legend('L_{Hall}','NC_{20}')


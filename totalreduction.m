clear 
close all
clc
%% double wall

wall1.rho = 700;                % density wall 1
wall1.t = 0.3;                  % thickness wall 1
wall1.m = wall1.rho * wall1.t;  % mass per unit area wall 1
wall1.eta = 0.025;             % loss factor wall 1
wall1.E = 10E9;                 % young's modulus wall 1
wall1.nu = 0.3;                 % poisson ratio wall 1

wall2.rho = 800;                % density wall 2
wall2.t = 0.2;                  % thickness wall 2
wall2.m = wall2.rho * wall2.t;  % mass per unit area wall 2
wall2.eta = 0.025;             % loss factor wall 2
wall2.E = 10E9;                 % young's modulus wall 2
wall2.nu = 0.3;                 % poisson ratio wall 2

d_z = 0.05;             % distance between walls
c_0 = 340;              % speed of sound
rho_0 = 1.2;            % density air 

f = 1:10000;            % frequency domain

% para.K = para.rho_air*para.c^2/para.d_z;
% para.m_eff = para.m1*para.m2/(para.m1+para.m2);
% 
% para.f_d = 1/(2*pi)*sqrt(para.K/para.m_eff);
% 
% lamda = para.c./f;      % wavelength
%%

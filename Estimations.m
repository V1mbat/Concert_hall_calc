clc
close all
clear

save_T = false;
%% reverb time and absorption of concert hall
% Mepfit Room  
% All tables are for the octave bands from 63 Hz to 4 kHz
f_oct = [63 125 250 500 1000 2000 4000];
c = 340;                                            % speed of sound

if save_T == false
    load T_max.mat;
end
%%

abs_occ = [0.6 0.62 0.72 0.78 0.81 0.84 0.85];        % absorption coeff seated
abs_empty = [0.54 0.58 0.68 0.74 0.77 0.78 0.80];      % absorption coeff empty
att_coeff = [0 0 0 0.001 0.002 0.004 0.0086];         % attanuation coeff

abs_absorber = [0.7 0.7	0.95 0.95 0.95 0.9 0.9];      % absorption coeff of absorber (example)
abs_absorber_low = [0.7 0.75 0.85 0.4 0.2 0.1 0.1];    % absorption coeff of low freq absorber
abs_wood = [0.15 0.15 0.11 0.10 0.07 0.06 0.07];
abs_glass = [0.2 0.18 0.06 0.04 0.05 0.02 0.02];
abs_brick = [0.03 0.03 0.03 0.03 0.04 0.05 0.07];
abs_stone = [0.01 0.01 0.01 0.015 0.02 0.02 0.02];

%%
% things you should adjust
V = 21000;                                          % Volume of concert hall
S_peraudiencemember = 0.5*0.55;                     % area per seat
S_absorber = 200;                                   % area of normal absorbers
S_absorber_low = 130;                               % area of low frequency absorbers

% 0/100 for high T_60
% 200/130 for low T_60

S_audience = 2300 * S_peraudiencemember;            % area of whole audience
S_concerthall_wall = 5850;
S_window = 117.6;

A_concerthall_wall = abs_wood*S_concerthall_wall;
A_absorber = - S_absorber * log(1-abs_absorber);    % eq. absorption area of absorbers
A_absorber_low = - S_absorber_low * log(1-abs_absorber_low);
A_glass = abs_glass*S_window;
A_occ = - S_audience*log(1-abs_occ);                % eq. absorption area of audience
A_empty = - S_audience*log(1-abs_empty);            % eq. absorption area of seats

A_total_occ = A_concerthall_wall + A_absorber + A_absorber_low + A_glass + A_occ; 
A_total_empty = A_concerthall_wall + A_absorber + A_absorber_low + A_glass + A_empty; 

% reverb time when seated and empty
T_60_calc_occ = 24*log(10)*V./(c*(4*att_coeff*V + A_total_occ));
T_60_calc_empty = 24*log(10)*V./(c*(4*att_coeff*V + A_total_empty));

%% clarity
r = 5:1:50;       % distance

E_direct = zeros(1,length(r));
E_direct(:) = 100./r.^2;

E_early_empty = zeros(length(T_60_calc_empty), length(r));
E_late_empty = zeros(length(T_60_calc_empty), length(r));
C_80_empty = zeros(length(T_60_calc_empty), length(r));     % clarity when emtpy

E_early_occ = zeros(length(T_60_calc_occ), length(r));
E_late_occ = zeros(length(T_60_calc_occ), length(r));
C_80_occ = zeros(length(T_60_calc_occ), length(r));         % clarity when fully seated

for i = 1:length(T_60_calc_empty)

    E_early_empty(i,:) = (31200*T_60_calc_empty(i)/V)*exp(-0.04.*r./T_60_calc_empty(i))*(1-exp(-1.11/T_60_calc_empty(i)));
    E_late_empty(i,:) = (31200*T_60_calc_empty(i)/V)*exp(-0.04.*r./T_60_calc_empty(i))*exp(-1.11/T_60_calc_empty(i));
    C_80_empty(i,:) = 10*log10((E_direct+E_early_empty(i,:))./E_late_empty(i,:));

    E_early_occ(i,:) = (31200*T_60_calc_occ(i)/V)*exp(-0.04.*r./T_60_calc_occ(i))*(1-exp(-1.11/T_60_calc_occ(i)));
    E_late_occ(i,:) = (31200*T_60_calc_occ(i)/V)*exp(-0.04.*r./T_60_calc_occ(i))*exp(-1.11/T_60_calc_occ(i));
    C_80_occ(i,:) = 10*log10((E_direct+E_early_empty(i,:))./E_late_empty(i,:));
end

%% strength

A = 0.163*V./T_60_calc_occ;             % equivilant absorption area
Q = 1;                                  % directivity
delta = 6.91./T_60_calc_empty;          % decay
ref = Q/(4*pi*10^2);                    % reference
dir = Q./(4*pi*r.^2);                   % direct component

diff = zeros(length(A), length(r));     % diffuse component
G = zeros(length(A), length(r));

for i = 1:length(A)
    diff(i,:) = 4./A(i) * exp(-2*delta(i)*r./c);
    G(i,:) = 10*log10((dir + diff(i,:))./ref);  % Strength for all six octave bands
end

%% rehersal room

% things you should adjust 
S_absorber_rehersal = 70;
S_absorber_low_rehersal = 20; 

n_person = 1:50;                             % number of people in rehersal room
V_rehersal = 1190;                       % Volume of rehersal room
S_wall = 538;         % area of the walls, floor and ceiling in rehersal room
% end

A_absorber_rehersal = - S_absorber_rehersal * log(1-abs_absorber);
A_absorber_rehersal_low = - S_absorber_low_rehersal * log(1-abs_absorber_low);
A_person = [0.02 0.05 0.16 0.25 0.58 0.86 1.03]; % eq. absorption area per person
A_wall = S_wall*abs_wood;                   % eq. absorption area of walls


% reverberation time of rehersal room
T_60_rehersal = zeros(length(n_person), length(A_wall));
for i = n_person
    T_60_rehersal(i,:) = 0.163*V_rehersal./(n_person(i)*A_person + A_wall + A_absorber_rehersal_low + A_absorber_rehersal);
end

%% lobby

V_lobby1 = 6275;
S_ceiling_L1 = 880;
S_walls_L1 = 780;
S_floor_L1 = 1255;
S_window_L1 = 125;
S_absorber_L1 = 50;
S_absorber_low_L1 = 250;

V_lobby2 = 4500;
S_ceiling_L2 = 1090;
S_walls_L2 = 680;
S_floor_L2 = 900;
S_window_L2 = 125;
S_absorber_L2 = 50;
S_absorber_low_L2 = 150;


A_ceiling_L1 = S_ceiling_L1*abs_wood;
A_walls_L1 = S_walls_L1*abs_brick;
A_floor_L1 = S_floor_L1*abs_stone;
A_window_L1 = S_window_L1*abs_glass;
A_absorber_L1 = - S_absorber_L1 * log(1-abs_absorber);
A_absorber_low_L1 = - S_absorber_low_L1 * log(1-abs_absorber);
A_total_L1 = A_ceiling_L1 + A_walls_L1 + A_floor_L1 + A_window_L1 + A_absorber_L1 + A_absorber_low_L1;

A_ceiling_L2 = S_ceiling_L2*abs_wood;
A_walls_L2 = S_walls_L2*abs_brick;
A_floor_L2 = S_floor_L2*abs_stone;
A_window_L2 = S_window_L2*abs_glass;
A_absorber_L2 = - S_absorber_L2 * log(1-abs_absorber);
A_absorber_low_L2 = - S_absorber_low_L2 * log(1-abs_absorber);
A_total_L2 = A_ceiling_L2 + A_walls_L2 + A_floor_L2 + A_window_L2 + A_absorber_L2 + A_absorber_low_L2;

T_60_Lobby1 = 24*log(10)*V_lobby1./(c*(4*att_coeff*V_lobby1 + A_total_L1));
T_60_Lobby2 = 24*log(10)*V_lobby2./(c*(4*att_coeff*V_lobby2 + A_total_L2));

%% MEPFIT
V_mepfit = 2.7*25*5 + 11*23*5 + 28*16*5;
S_mepfit = 2*(2.7*25 + 11*23 + 28*16) + 2*(25+23+28)*5;

A_mepfit = S_mepfit*abs_wood;

%%
if save_T
    T_max = [T_60_calc_empty; T_60_calc_occ];
    save('T_max.mat','T_max');
    return;
end
%% Plots 
% adjust colors / add more colors
col1 = '[0, 0.4470, 0.7410]';
col2 = '[0.8500 0.3250 0.0980]';

%%
% T_60 in concert hall
fig_name = 'T_60_group2';
f = figure;
semilogx(f_oct, T_60_calc_empty,'--', color = col1)
hold on
semilogx(f_oct, T_60_calc_occ, color = col1)
semilogx(f_oct, T_max(1,:),'--', color = col2)
semilogx(f_oct, T_max(2,:), color = col2)
thickenall_big;
xticks(f_oct);
xticklabels({'63', '125', '250', '500', '1k', '2k', '4k'})
xlim([0 5000])
ylim([0 2.5])
xlabel('f in Hz')
ylabel('T_{60} in s')

%set(gcf, 'color', 'none');

exportPlot(f,'Plots/', fig_name, true)

%% 

% clarity and strength
fig_name = 'C_80_G_group2';
f = figure;

plot(r, C_80_occ(4,:), color = col1)
hold on
plot(r, G(4,:), color = col2)
thickenall_big;
grid on
xlabel('Distance in m')
ylabel('C_{80}/G in dB')
%legend('Clarity C_{80}','Strength G')
xlim([5 35])
ylim([-2 10])
exportPlot(f,'Plots/', fig_name, true)

return;
%%

% lobby
figure
semilogx(f_oct, T_60_Lobby1)
hold on
grid on
semilogx(f_oct, T_60_Lobby2)
thickenall_big;
xticks(f_oct);
xticklabels({'63','125', '250', '500', '1k', '2k', '4k'})
xlim([0 5000])
ylim([0 3])
xlabel('f in Hz')
ylabel('T_{60} in s')
legend('Lobby1','Lobby2')

% T_60 in rehersal room
figure
semilogx(f_oct, T_60_rehersal(1,:))
hold on
grid on
semilogx(f_oct, T_60_rehersal(25,:))
semilogx(f_oct, T_60_rehersal(50,:))
thickenall_big;
xticks(f_oct);
xticklabels({'63', '125', '250', '500', '1k', '2k', '4k'})
%xlim([0 5000])
ylim([0 3])
xlabel('f in Hz')
ylabel('T_{60} in s')
legend('1 Musician','25 Musicians', '50 Musicians')

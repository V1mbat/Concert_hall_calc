clc
close all
clear
% to do:
% add absorption to change T_60
% done: T_60 of rehersal hall
% done: C_80 and G dependant on T_60_calc
% desired T_60 for rehersal room?

%% reverb time and absorption of audience and seats
% All tables are for the octave bands from 125 Hz to 4 kHz

V = 18000;          %Volume
c = 340;            % speed of sound

f_oct = [125 250 500 1000 2000 4000];

area_peraudiencemember = 0.5*0.55;                  % area per seat

area_audience = 2300 * area_peraudiencemember;      % area of whole audience

abs_seated = [0.62 0.72 0.8 0.83 0.84 0.85];        % absorption coeff seated
abs_empty = [0.54 0.62 0.68 0.7 0.68 0.66];         % absorption coeff empty
att_coeff = [0 0 0.0024 0.0042 0.0089 0.0262];      % attanuation coeff

% reverb time when seated and empty
T_60_calc_seated = 24*log(10)*V./(c*(4*att_coeff*V - area_audience*log(1-abs_seated)));
T_60_calc_empty = 24*log(10)*V./(c*(4*att_coeff*V - area_audience*log(1-abs_empty)));

%% clarity
r = 1:0.1:50;

E_direct = zeros(1,length(r));
E_direct(:) = 100./r.^2;

E_early_empty = zeros(length(T_60_calc_empty), length(r));
E_late_empty = zeros(length(T_60_calc_empty), length(r));
C_80_empty = zeros(length(T_60_calc_empty), length(r));

E_early_occ = zeros(length(T_60_calc_seated), length(r));
E_late_occ = zeros(length(T_60_calc_seated), length(r));
C_80_occ = zeros(length(T_60_calc_seated), length(r));

for i = 1:length(T_60_calc_empty)

    E_early_empty(i,:) = (31200*T_60_calc_empty(i)/V)*exp(-0.04.*r./T_60_calc_empty(i))*(1-exp(-1.11/T_60_calc_empty(i)));
    E_late_empty(i,:) = (31200*T_60_calc_empty(i)/V)*exp(-0.04.*r./T_60_calc_empty(i))*exp(-1.11/T_60_calc_empty(i));
    C_80_empty(i,:) = 10*log10((E_direct+E_early_empty(i,:))./E_late_empty(i,:));

    E_early_occ(i,:) = (31200*T_60_calc_seated(i)/V)*exp(-0.04.*r./T_60_calc_seated(i))*(1-exp(-1.11/T_60_calc_seated(i)));
    E_late_occ(i,:) = (31200*T_60_calc_seated(i)/V)*exp(-0.04.*r./T_60_calc_seated(i))*exp(-1.11/T_60_calc_seated(i));
    C_80_occ(i,:) = 10*log10((E_direct+E_early_empty(i,:))./E_late_empty(i,:));
end

%% strength

A = 0.163*V./T_60_calc_seated;          % equivilant absorption area
Q = 1;                                  % directivity
delta = 6.91./T_60_calc_empty;          % decay
ref = Q/(4*pi*10^2);
dir = Q./(4*pi*r.^2);

diff = zeros(length(A), length(r));
G = zeros(length(A), length(r));

for i = 1:length(A)
    diff(i,:) = 4./A(i) * exp(-2*delta(i)*r./c);
    G(i,:) = 10*log10((dir + diff(i,:))./ref);
end

%% rehersal room

abs_area_person = [0.05 0.16 0.25 0.58 0.86 1.03];
n_person = 100;                             % number of people in rehersal room

V_rehersal = 20*10*5;                          % Volume of rehersal room
area_wall = 2 * (20*10 + 20*5 + 10*5);      % area of the walls, floor and ceiling in rehersal room

abs_wall = [0.2 0.18 0.15 0.12 0.1 0.1];    % absorption coefficient of wall
abs_area_wall = area_wall*abs_wall;         % eq. absorption area of walls

% reverberation time of rehersal room
T_60_rehersal = 0.163*V_rehersal./(n_person*abs_area_person + abs_area_wall);

%% plots

% clarity
figure
plot(r, C_80_occ)
hold on
thickenall_big;
grid on
xlabel('Distance r in m')
ylabel('Clarity C_{80} in dB')
legend('125','250','500','1k','2k','4k')
%ylim([-5 20])

% strength
figure
plot(r, G)
thickenall_big;
grid on 
legend('125','250','500','1k','2k','4k')
xlabel('Distance r in m')
ylabel('Strength G in dB')
xlim([0 40])
ylim([-5 20])

% T_60
figure
semilogx(f_oct, T_60_calc_empty)
hold on
grid on
semilogx(f_oct, T_60_calc_seated)
thickenall_big;
xticks(f_oct);
xticklabels({'125', '250', '500', '1k', '2k', '4k'})
xlim([0 5000])
ylim([0 4])
xlabel('f in Hz')
ylabel('T_{60} in s')
legend('empty','seated')


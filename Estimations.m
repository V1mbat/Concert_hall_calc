clc
close all
clear
%% clarity
V = 18000;
T = 2;

r = 1:0.1:50;

E_direct = 100./r.^2;

E_early = (31200*T/V)*exp(-0.04.*r./T)*(1-exp(-1.11/T));

E_late = (31200*T/V)*exp(-0.04.*r./T)*exp(-1.11/T);

C_80 = 10*log10((E_direct+E_early)./E_late);

%% strength

A = 0.163*V/T;      % equivilant absorption area
Q = 1;              % directivity
delta = 6.91/T;     % decay
c = 340;            % speed of sound
dir = Q./(4*pi*r.^2);
diff = 4/A * exp(-2*delta*r./c);
ref = Q/(4*pi*10^2);

G = 10*log10((dir + diff)./ref);

%% absorption layer
f = 1:500;
c = 340;
lambda = c./f;

abs = lambda./4;

%% reverb time and absorption of audience and seats
% All tables are for the octave bands from 125 Hz to 4 kHz

f_oct = [125 250 500 1000 2000 4000];

area_peraudiencemember = 0.5*0.55;                  % area per seat

area_audience = 2300 * area_peraudiencemember;      % area of whole audience

L_facade = [73 75 80 66 54 50];                     % SPL at facade
NC_20 = [40	33	26	22	19	17];                    % noise requirement in concert hall
abs_seated = [0.62 0.72 0.8 0.83 0.84 0.85];        % absorption coeff seated
abs_empty = [0.54 0.62 0.68 0.7 0.68 0.66];         % absorption coeff empty
att_coeff = [0 0 0.0024 0.0042 0.0089 0.0262];      % attanuation coeff

T_60_calc_seated = 24*log(10)*V./(c*(4*att_coeff*V - area_audience*log(1-abs_seated)));
T_60_calc_empty = 24*log(10)*V./(c*(4*att_coeff*V - area_audience*log(1-abs_empty)));

%% plots

figure
plot(r, C_80)
thickenall_big;
grid on
xlabel('Distance r in m')
ylabel('Clarity C_{80} in dB')
ylim([-5 20])

figure
plot(r, G)
thickenall_big;
grid on 
xlabel('Distance r in m')
ylabel('Strength G in dB')
xlim([0 40])
ylim([-5 20])

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


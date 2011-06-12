%T = 2210;  % kelvin
%T = 2503.7
T = 2503.7
%T = 3000
%T = 3500
%T = 1200
%T = 800
rs = 5.90; % bohr
%rho = 219  % kg/m^3
%rho = 400  % kg/m^3
%rho = 310  % kg/m^3
%rho = 470  % kg/m^3
%rho = 740  % kg/m^3
%rho = 930  % kg/m^3
rho = 380
%rho = 300

% Na 
 M = 22.989770 % Na
% rho = 219
 Tc = 2503

% Cs
% M = 132.90545 % Cs
% rho = 380
% Tc = 1651


% Rb
% M = 85.4678
% rho = 290
% Tc = 1744

% K
% M = 39.0983
% rho = 180
% Tc = 1905

% Li
% M = 6.941
% rho = 110
% Tc = 3000

n = rho*1000*5.2917721e-11^3/M*6.022145e23;

N = 16;    % particles
kB = 1.3807e-23; % Joules/Kelvin
kB = kB * 6.24150974e18;  % eV/Kelvin
kB = kB / 27.211383;       % hartree/Kelvin
kBT = T * kB;

r_s = (3/(4*pi*n))^(1/3)
Ef = (3*n/pi)^(2/3)*pi^2/2
Tf = Ef / kB
ratio = Tc/Tf
beta = 1.0/kBT;

V1 = 4*pi/3 * rs^3 * N;
V2 = N/n;
L1 = V1^(1/3);
L2 = V2^(1/3)
tau = 0.125 % inverse hartrees
M = round(beta/tau)
Tsim = 1/(kB*tau*M)

%T = 2210;  % kelvin
%T = 2503.7
T = 2503.7
rs = 5.90; % bohr
%rho = 219  % kg/m^3
rho = 240  % kg/m^3

n = rho*1000*5.2917721e-11^3/22.98977*6.022145e23;

N = 16;    % particles
kB = 1.3807e-23; % Joules/Kelvin
kB = kB * 6.24150974e18;  % eV/Kelvin
kB = kB / 27.211383       % hartree/Kelvin
kBT = T * kB;

beta = 1.0/kBT


V1 = 4*pi/3 * rs^3 * N;
V2 = N/n;
L1 = V1^(1/3);
L2 = V2^(1/3)
tau = 0.125 % inverse hartrees
M = round(beta/tau)
Tsim = 1/(kB*tau*M)

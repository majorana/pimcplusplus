T = 2210;  % kelvin
rs = 5.90; % bohr
N = 16;    % particles
kB = 1.3807e-23 % Joules/Kelvin
kB = kB * 6.24150974e18;  % eV/Kelvin
kB = kB / 27.211383;       % hartree/Kelvin
kBT = T * kB;

beta = 1.0/kBT

V = 4*pi/3 * rs^3 * N
L = V^(1/3)
tau = 0.0625 % inverse hartrees
M = beta/tau

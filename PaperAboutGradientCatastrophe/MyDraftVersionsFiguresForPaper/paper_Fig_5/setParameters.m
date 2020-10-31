%% --------------------------------------------------------------------- %%
%% Physical parameters.
sigmaX = 3.2e-6;       % Elliptical waveguide size along the x-axis. (in metres)
sigmaY = 4.9e-6;       % Elliptical waveguide size along the y-axis. (in metres)
a0 = 19e-6;            % Distance between waveguide A and waveguide B. (in metres)
a = sqrt(3) * a0;      % Waveguide separation. (in metres)
n0 = 1.47;             % Background refractive index.
nA = 2.6e-3;           % Refractive index modulation sublattice A.
nB = 2.8e-3;           % Refractive index modulation sublattice B.
n2 = 3e-20;            % Fused silica Kerr coefficient determined
                       % the strength of self-focusing nonlinearity. (in m^2/W)
lambda = 1.65e-6;                % Wavelength. (in metres)
k0 = 2 * pi * n0 / lambda;       % Wavenumber.
%% --------------------------------------------------------------------- %%
%% Dimensionless variables.
w0 = 10e-6;            % Transverse length scale in the xy-plane. (in metres)
z0 = 2 * k0 * w0^2;    % Characteristic length along the z-axis.
I0 = 1e+16;            % Characteristic irradiance scale. (in W/m^2)
sigmaX_d = sigmaX / w0;          % Dimensionless waveguide x-size.
sigmaY_d = sigmaY / w0;          % Dimensionless waveguide y-size.
a_d = a / w0;                    % Dimensionless waveguide separation.
a0_d = a0 / w0;
nA_d = k0 * z0 * nA / n0;        % Normalised potential depth.
nB_d = k0 * z0 * nB / n0;        % Normalised potential depth.
gNL = I0 * k0 * z0 * n2 / n0;    % Normalised nonlinearity strength.
%% --------------------------------------------------------------------- %%
%% Parameters of lattice geometry.
numCellsX = 2^0;      % Number of unit cells in x-direction.
numCellsY = 2^3;      % Number of unit cells in y-direction.
Lx = numCellsX * a_d;
Ly = numCellsY * 3 * a0_d;
%% --------------------------------------------------------------------- %%
%% Parameters of numerical scheme.
Nx = numCellsX * 2^5;
Ny = numCellsY * 2^5;
%% --------------------------------------------------------------------- %%
%% Parameters of TBM model.
t_d=0.382;
M_d=0.191;
VD=sqrt(3)/2*t_d*a_d;
eta_d=t_d*a_d^2/8;
%% Parameters of I_0(x).
widthWP = 5 * a_d / sqrt(2);


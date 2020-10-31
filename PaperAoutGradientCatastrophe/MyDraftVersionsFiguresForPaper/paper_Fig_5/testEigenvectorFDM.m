%% --------------------------------------------------------------------- %%
% clear;
% close all;
setParameters;
%% --------------------------------------------------------------------- %%
x = makeSpatialGrid( Nx, Lx );
y = makeSpatialGrid( Ny, Ly );
latticePotential = makeLatticePotential( x, y, numCellsX, numCellsY, ...
                                         a0_d, nA_d, nB_d, sigmaX_d, sigmaY_d, 1 );
%% --------------------------------------------------------------------- %%
tic;
kappaFDM = 0.668241871033659;

% omegaFDM = 0.972383523828923;
% omegaFDM = 1.026264930037803;
omegaFDM = 0.901707145624187;
% omegaFDM = 0.974293003055888;

% omegaFDM = 1.838742458254975;
% omegaFDM = 1.868913966297406;

kappaFDM = ( 2 * pi / 3 ) / a_d;
omegaFDM = 0.929758407192387;
% omegaFDM = 0.947879101849542;
% omegaFDM = 0.889337910549110;
% omegaFDM = 1.005277478352531;


psiFDM = linearMode_FDM( latticePotential, Lx, Ly, kappaFDM, omegaFDM );
toc;
%% --------------------------------------------------------------------- %%
figure;
subplot( 1, 3, 1 );
imagesc( x( 1, 1 : 1 : end ), y( 1, 1 : 1 : end ), abs( psiFDM ).^2 );
box on;
xlim( [ x( 1, 1 ), x( 1, end ) ] );
ylim( [ y( 1, 1 ), y( 1, end ) ] );
subplot( 1, 3, 2 );
imagesc( x( 1, 1 : 1 : end ), y( 1, 1 : 1 : end ), real( psiFDM ) );
box on;
xlim( [ x( 1, 1 ), x( 1, end ) ] );
ylim( [ y( 1, 1 ), y( 1, end ) ] );
subplot( 1, 3, 3 );
imagesc( x( 1, 1 : 1 : end ), y( 1, 1 : 1 : end ), imag( psiFDM ) );
box on;
xlim( [ x( 1, 1 ), x( 1, end ) ] );
ylim( [ y( 1, 1 ), y( 1, end ) ] );
%% --------------------------------------------------------------------- %%

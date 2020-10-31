%% --------------------------------------------------------------------- %%
addpath('C:\Users\smoli\Dropbox\Topological Nonlinear 2D\Paraxial_ver_2')
addpath('C:\Users\smoli\Dropbox\Topological Nonlinear 2D\Paraxial_ver_2\Tests')
close all;
% format long;
clear max_diff
strDataSave = 'on';
strVisualization = 'on';
%% --------------------------------------------------------------------- %%
%% Parameters of the problem of edge pulses propagation along topological domain walls.
setParameters;
numCellsX = 2^6;
Nx = numCellsX * Nx;
Lx = numCellsX * a_d;
nonlinearCoefficient = 0 * gNL;
%-------------------------------------------------------------------------%
z = 0;
dz = 0.01;
Num_Save_Zsteps = 2;
Num_Extra_Zsteps = 1000;
if( nonlinearCoefficient == 0 )
    problemCase = 'LinearCase';
else
    problemCase = 'NonlinearCase';
end
%% --------------------------------------------------------------------- %%
%% Building of grids in direct and inverse spaces, a lattice potential, and a dispersion operator.
[ x, kx ] = makeSpatialGrid( Nx, Lx );
[ y, ky ] = makeSpatialGrid( Ny, Ly );
latticePotential = makeLatticePotential( x, y, numCellsX, numCellsY, ...
                                         a0_d, nA_d, nB_d, sigmaX_d, sigmaY_d, 1 );
kDispersion = zeros( Ny, Nx );
for iy = 1 : Ny
    for ix = 1 : Nx
        kDispersion( iy, ix ) = - ( kx( 1, ix )^2 + ky( 1, iy )^2 );
    end
end
%% --------------------------------------------------------------------- %%
if( strcmp( strDataSave, 'on' ) )
    foldername = 'Results/';
    if( ~exist( foldername, 'dir' ) )
       mkdir( foldername );
    end
    foldername = strcat( foldername, 'Dynamics/' );
    if( ~exist( foldername, 'dir' ) )
       mkdir( foldername );
    end
    foldername = strcat( foldername, 'Data/' );
    if( ~exist( foldername, 'dir' ) )
       mkdir( foldername );
    end
    foldername = strcat( foldername, problemCase, '/' );
    if( ~exist( foldername, 'dir' ) )
       mkdir( foldername );
    end
    filename = strcat( foldername, 'x.dat' );
    save( filename, 'x', '-ascii', '-double' );
    filename = strcat( foldername, 'y.dat' );
    save( filename, 'y', '-ascii', '-double' );
end
%% --------------------------------------------------------------------- %%
%% Constraction of an initial wave function.
amplitudeWP = 0.5;
xpositionWP = ( numCellsX / 4 ) * a_d;
psi = zeros( Ny, Nx );
for ixCell = 1 : numCellsX
    psi( :, ( Nx / numCellsX ) * ( ixCell - 1 ) + 1 : ( Nx / numCellsX ) * ixCell ) = psiFDM( :, : );
end
for iy = 1 : Ny
    psi( iy, : ) = amplitudeWP.* exp( - ( ( x - xpositionWP )./ widthWP ).^2./2 ).* ...
                                 exp( - 1i.* kappaFDM.* ( x - xpositionWP ) ).* psi( iy, : );
end
if( strcmp( strVisualization, 'on' ) )
    fig = figure;
    maxSqrAbsPsi = zeros( Num_Save_Zsteps + 1, 4 );
    maxSqrAbsPsi( 1, : ) = dynamicsVisualization( psi, x, y, z );
end
if( strcmp( strDataSave, 'on' ) )
    realPsi = real( psi );
    imagPsi = imag( psi );
    filename = strcat( foldername, 'realPsi_z=', num2str( z ), '.dat' );
    save( filename, 'realPsi', '-ascii', '-double' );
    filename = strcat( foldername, 'imagPsi_z=', num2str( z ), '.dat' );
    save( filename, 'imagPsi', '-ascii', '-double' );
end
%% --------------------------------------------------------------------- %%
%% Numerical propagation of the paraxial equation with the split-step Fourier method.
tic;
% figure;
for jz = 1 : Num_Save_Zsteps
    for iz = 1 : Num_Extra_Zsteps
        psi = exp( ( 1i * dz / 2 ).* ( latticePotential + nonlinearCoefficient.* abs( psi ).^2 ) ).* psi;
        psi = fft2( psi );
        psi = exp( ( 1i * dz ).* kDispersion ).* psi;
        psi = ifft2( psi );
        psi = exp( ( 1i * dz / 2 ).* ( latticePotential + nonlinearCoefficient.* abs( psi ).^2 ) ).* psi;
        z = z + dz;
    end    
    if( strcmp( strVisualization, 'on' ) )
        maxSqrAbsPsi( jz + 1, : ) = dynamicsVisualization( psi, x, y, z );
    end
    if( strcmp( strDataSave, 'on' ) )
        realPsi = real( psi );
        imagPsi = imag( psi );
        filename = strcat( foldername, 'realPsi_z=', num2str( z ), '.dat' );
        save( filename, 'realPsi', '-ascii', '-double' );
        filename = strcat( foldername, 'imagPsi_z=', num2str( z ), '.dat' );
        save( filename, 'imagPsi', '-ascii', '-double' );
    end
end
toc;
%% --------------------------------------------------------------------- %%
if( strcmp( strVisualization, 'on' ) && strcmp( strDataSave, 'on' ) )
    filename = strcat( foldername, 'maxSqrAbsPsi.dat' );
    save( filename, 'maxSqrAbsPsi', '-ascii', '-double' );
end
%% --------------------------------------------------------------------- %%

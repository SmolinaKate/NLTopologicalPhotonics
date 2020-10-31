%% --------------------------------------------------------------------- %%
close all;
% format long;
strDataSave = 'on';
strVisualization = 'on';
%% --------------------------------------------------------------------- %%
%% Parameters of the problem of edge pulses propagation along topological domain walls.
setParameters;
numCellsX = 2^6;
Nx = numCellsX * Nx;
Lx = numCellsX * a_d;
nonlinearCoefficient = 1 * gNL;
%-------------------------------------------------------------------------%
z = 0;
dz = 0.01;
Num_Save_Zsteps = 7000;
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
%% Range of amplitudes.
amplitudeWP_array=0.4:0.025:0.55;
%% --------------------------------------------------------------------- %%
MaxAbsPsi=zeros(Num_Save_Zsteps,length(amplitudeWP_array));
max_diff=zeros(Num_Save_Zsteps,length(amplitudeWP_array));
max_diff_norm=zeros(Num_Save_Zsteps,length(amplitudeWP_array));
for amplitudeWPind = 1:1:length(amplitudeWP_array)
        tic;
%% --------------------------------------------------------------------- %%
%% Initial distribution.
        amplitudeWP = amplitudeWP_array(amplitudeWPind); 
        xpositionWP = ( numCellsX / 4 ) * a_d;
        psi = zeros( Ny, Nx );
        for ixCell = 1 : numCellsX
            psi( :, ( Nx / numCellsX ) * ( ixCell - 1 ) + 1 : ( Nx / numCellsX ) * ixCell ) = psiFDM( :, : );
        end
        for iy = 1 : Ny
            psi( iy, : ) = amplitudeWP.* exp( - ( ( x - xpositionWP )./ widthWP ).^2./2 ).* ...
                                         exp( - 1i.* kappaFDM.* ( x - xpositionWP ) ).* psi( iy, : );
        end
%       Finding the position of the domain wall.
        [max_val, numxpositionWP] = (max( abs( psi( Ny/2, : ) ).^2 ));
        [max_val_1, y_domain_wall ] = (max( abs( psi( :, numxpositionWP ) ).^2 ));
        line_domain_wall(1:length(x))=y(y_domain_wall);
%       Plot fig to check the position of the domain wall.
%       figure; imagesc(x,y,abs(psi)); hold on;
%       plot(x,line_domain_wall,'r','LineWidth',2); axis equal;
%       xlabel('$x$'); ylabel('$y$')
%       MakeFigNice
%% Numerical propagation of the paraxial equation with the split-step Fourier method.
        for jz = 1 : Num_Save_Zsteps
                psi = exp( ( 1i * dz / 2 ).* ( latticePotential + nonlinearCoefficient.* abs( psi ).^2 ) ).* psi;
                psi = fft2( psi );
                psi = exp( ( 1i * dz ).* kDispersion ).* psi;
                psi = ifft2( psi );
                psi = exp( ( 1i * dz / 2 ).* ( latticePotential + nonlinearCoefficient.* abs( psi ).^2 ) ).* psi;
                z = z + dz;
%% Calculating maps.
%% A map of maxs:
%        Finding the position of waveguides' centers.
                c=0;
                for ind=1:1:numCellsX
                    number_of_array_x_where_cen_of_waveg=ceil((sqrt(3)*a0_d*(1-0.25)+sqrt(3)*a0_d*c)/abs(x(1)-x(2)));
                    if(number_of_array_x_where_cen_of_waveg<length(x))
                        x_where_cen_of_waveg(ind)=x(number_of_array_x_where_cen_of_waveg);
                        Dist_where_cen_of_waveg(ind)=abs( psi( y_domain_wall, number_of_array_x_where_cen_of_waveg ) ).^2;
                        Dist_where_cen_of_waveg_norm(ind)=abs( psi( y_domain_wall, number_of_array_x_where_cen_of_waveg ) ).^2./max(abs( psi( y_domain_wall, : ) ).^2);
                    end
                    c=c+1;
                end
%% Envelope (z,x_(positions of waveguides)):
                Envelope(jz,:,amplitudeWPind)=Dist_where_cen_of_waveg(:);
      %Plot fig to check the position of the waveguides' centers.
%       figure; plot(x, abs(psi( y_domain_wall,:).^2)); 
%       hold on; plot(x_where_cen_of_waveg,Dist_where_cen_of_waveg)
%       xlabel('$x$')
%       ylabel('$|\psi|^2$')
%       MakeFigNice
                Dist_diff=diff(Dist_where_cen_of_waveg);
                Dist_diff_norm=diff(Dist_where_cen_of_waveg_norm);
%% A map of differences (derivatives):
                max_diff(jz,amplitudeWPind)=max(abs(Dist_diff));
%% A map of normalized differences (derivatives):
                max_diff_norm(jz,amplitudeWPind)=max(abs(Dist_diff_norm));
                clear Dist_where_cen_of_waveg Dist_diff x_where_cen_of_waveg
        end
        toc;
end

z_array=dz:dz:dz*Num_Save_Zsteps;

filename = strcat( foldername, 'max_diff.dat' );
save( filename, 'max_diff', '-ascii', '-double' );

filename = strcat( foldername, 'max_diff_norm.dat' );
save( filename, 'max_diff_norm', '-ascii', '-double' );

filename = strcat( foldername, 'z_array.dat' );
save( filename, 'z_array', '-ascii', '-double' );

filename = strcat( foldername, 'amplitudeWP_array.dat' );
save( filename, 'amplitudeWP_array', '-ascii', '-double' );

filename = strcat( foldername, 'Envelope.mat' );
save(filename,'Envelope','Envelope'); 




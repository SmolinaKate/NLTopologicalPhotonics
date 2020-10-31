%% --------------------------------------------------------------------- %%
clear;  
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster\profile_evolution_with_d2_data';
%% --------------------------------------------------------------------- %%
link=1;
a=abs(2/link/sqrt(3));
% a=0.1;
Kp=4*pi/3/a;
Kp1=4*pi/3/a;
a0=sqrt(3)*a;
dt = 2;

Num_Save_Tsteps = 3;
Num_Extra_Tsteps = 2000;
l_x = 1/link;
N_x = 2^8;
L_x = l_x * N_x;
x = zeros( 1, N_x );
for ix = 1 : N_x
    x( 1, ix ) = ( ix - ( 1 + N_x / 2 ) ) * l_x;
end
k_x = zeros( 1, N_x );
for ix = 1 : N_x
    k_x( 1, ix ) = ( ix - ( 1 + N_x / 2 ) ) * ( 2 * pi / L_x );
end
k_x = fftshift( k_x );
l_y = 1/link;
N_y = 2^8;
L_y = l_y * N_y;

y = zeros( 1, N_y );
for iy = 1 : N_y
    y( 1, iy ) = ( iy - ( 1 + N_y / 2 ) ) * l_y;
end

k_y = zeros( 1, N_y );
for iy = 1 : N_y
    k_y( 1, iy ) = ( iy - ( 1 + N_y / 2 ) ) * ( 2 * pi / L_y );
end
k_y = fftshift(k_y);

max(k_y)
min(k_y)
Kp

M=0.07;
g=0.01;
B=0.001;
A=1;


tbr=4*M^2*sqrt(exp(1))/g^2/sqrt(B)/A^2;
[ X, Y ] = meshgrid( x, y );
psiA = zeros( N_y, N_x );
psiB = zeros( N_y, N_x );
epsilon = zeros(N_y,N_x);
y0=0; 
x0=0;

for ii=1:1:N_x
    if abs(x(ii)-x0)<l_x
       Nx0 = ii; 
    end
end

Nep=find(y==y0);
Nem=find(y==y0);

for ix = 1 : N_x
    for iy = 1 : N_y        
     epsilon( iy, ix ) = M* ( 2 * (heaviside( y( 1, iy ))-0.5));    
     I1(ix)= A*exp ( -B*(x(ix)-x0)^2);     
     k (ix)=  - g*I1(ix)/4;
     [psiA(iy,ix), psiB(iy,ix)] = YSolitonDistribution ( y(iy), - y0, M, k(ix), g, I1(ix));
         if ix>1
             xp = x(1):l_x:x(ix);
             kp=-g/4.*(A*(exp ( -B*(xp-x0).^2)));
             phi = trapz(xp,kp);
             psiA( iy, ix ) = psiA( iy, ix ) * exp ( 1i*phi);
             psiB( iy, ix ) = psiB( iy, ix ) * exp ( 1i*phi);
         else
             psiA( iy, ix ) = psiA( iy, ix ) * exp ( 1i*x(ix)*k(1));
             psiB( iy, ix ) = psiB( iy, ix ) * exp ( 1i*x(ix)*k(1));
         end  
    end       
end


file_name = sprintf('initial_parameters');
ToSaveFileName = [file_name '.mat'];
ToSaveFileName1 = fullfile(path, ToSaveFileName);
save(ToSaveFileName1,'psiA','psiB','x','y','dt','Num_Save_Tsteps','Num_Extra_Tsteps','Nem','M','g','A','B');                    

figure;
Visualization( psiA, psiB, x, y, Nem, Nx0, [ 0, x0, y0 ] );

%% --------------------------------------------------------------------- %%
tic;
strVisualization = 'on';
t = 0;    xSoliton = x0;    ySoliton = y0;
SolitonTrajectory = zeros( Num_Save_Tsteps + 1, 3 );
SolitonTrajectory( 1, : ) = [ t, xSoliton, ySoliton ];
der_x = max(diff(abs( psiA(Nem,:)).^2)/l_x);
der_y = max(diff(abs( psiA(:,Nx0)).^2)/l_y);

der=[der_y, der_x];

for it = 1 : Num_Save_Tsteps
     file_name = sprintf('pr_ev#%it',it);
     for jt = 1 : Num_Extra_Tsteps
                [ psiA, psiB ] = NonlinearStep( dt/2, psiA, psiB, epsilon,g);
                [ psiA, psiB ] = LinearStepD2_1( dt, psiA, psiB, k_x, k_y, link, a);
                [ psiA, psiB ] = NonlinearStep( dt/2, psiA, psiB, epsilon,g);
                max_val = max(diff(abs( psiA(Nem,:)).^2)/l_x);
                der=[der,max_val];
                t = t + dt;
      end     
                     [ xSoliton, ySoliton ] = SolitonPosition( psiA, psiB, x, y );
                    SolitonTrajectory( it + 1, : ) = [ t, xSoliton, ySoliton ];
                    if( strcmp( strVisualization, 'on' ) )
                        Visualization(psiA, psiB, x, y, Nem, Nx0, SolitonTrajectory( 1 : it + 1, : ) );
                    else
                        disp( strcat( num2str( it ), ' / ', num2str( Num_Save_Tsteps ) ) );
                    end
        ToSaveFileName = [file_name '.mat'];
        ToSaveFileName1 = fullfile(path, ToSaveFileName);
        save(ToSaveFileName1,'psiA','psiB');                    
                    
end
toc;
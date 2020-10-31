%% --------------------------------------------------------------------- %%
clear;    close all;
%% --------------------------------------------------------------------- %%
dt = 0.05;
Num_Save_Tsteps = 30;
Num_Extra_Tsteps = 30;
l_x = 0.01;
N_x = 2^10;
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
l_y = 0.01;
N_y = 2^10;
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
%% --------------------------------------------------------------------- %%
epsilon = zeros(N_y,N_x);
y0=0; 
x0=0;
Nx0=find(x==0);
Nep=find(y==y0);
Nem=find(y==y0);
M=2;
I1=1; 
g=1;
[ X, Y ] = meshgrid( x, y );
psiA = zeros( N_y, N_x );
psiB = zeros( N_y, N_x );
for ix = 1 : N_x
    for iy = 1 : N_y    
          epsilon( iy, ix ) = M*( 2 * (heaviside( y( 1, iy ) - y0)-0.5));
          k = - g*I1/4;
          [psiA(iy,ix), psiB(iy,ix)] = YSolitonDistribution ( y(iy), y0, M, k, g, I1);
          psiA(iy,ix) = psiA(iy,ix).* exp (1i*k.*x(ix));
          psiB(iy,ix) = psiB(iy,ix).* exp (1i*k.*x(ix));
    end       
end


figure; 


subplot( 2, 2, 4 );
hold on;
plot(abs(psiB(:,Nx0)),y, 'r', 'Linewidth',1);
ylabel( '$|\psi_2|$');
xlabel( '$y$');

        

grid on; 
box on;

         
subplot( 2, 2, 3 );
hold on;
imagesc( x, y, abs( psiB )); axis equal; colormap hot;
hold off;
box on;
daspect( [ 1, 1, 1 ] );
xlabel( '$x$' );
ylabel( '$y$' );
       
         
subplot( 2, 2, 2 );
hold on;
plot(abs(psiA(:,Nx0)),y, 'r', 'Linewidth',1);
ylabel( '$|\psi_1|$' );
xlabel( '$y$');
        
grid on; 
box on;

         
subplot( 2, 2, 1 );
hold on;
imagesc( x, y, abs( psiA )); axis equal; colormap hot;
hold off;
box on;
daspect( [ 1, 1, 1 ] );
xlabel( '$x$');
ylabel( '$y$');
       

        
        
a = findobj(gcf); 
size = 24;

allaxes  = findall(a,'Type','axes');
% alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

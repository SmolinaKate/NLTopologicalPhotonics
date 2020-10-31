%% --------------------------------------------------------------------- %%
clear;    close all;
%% --------------------------------------------------------------------- %%
dt = 0.05;
Num_Save_Tsteps = 30;
Num_Extra_Tsteps = 30;
l_x = 0.01;
N_x = 2^9;
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
N_y = 2^9;
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


r=250;
figure('Renderer', 'painters', 'Position', [10 10 2*r r])
clf

% create panel
p = panel();
border = 0; 
size = 20;


p.pack(1, 3);


p(1, 2).marginright = border; 
p(1, 2).marginleft = border; 
p(1, 1).marginright = border; 
p(1, 1).marginleft = border; 
p(1,3).marginright = border; 
p(1, 3).marginleft = border; 


p(1, 1).select();
plot(abs(psiA(:,Nx0)),y, 'r', 'Linewidth',1);
xlabel( '$|\psi_2|$');
ylabel( '$y$');
set(gca, 'xtick', [0:0.5:1], 'ytick', [-5:2.5:5],'TickDir','out'); 
xlim( [0,1] );
ylim( [-2.5,2.5] );
grid on;
box on; 


p(1, 2).select();
hold on;
imagesc( x, y, abs( psiA )); axis equal; colormap hot;
yticklabels([])
hold off;
box on;
%daspect( [ 1, 1, 1 ] );
xlabel( '$x$');
set(gca, 'xtick', [-2:2:2], 'ytick', [-2:2:2],'TickDir','out'); 
xlim( [-2.5,2.5] );
ylim( [-2.5,2.5] );
grid on;
box on; 
ax = gca;
ax.YAxisLocation = 'left';

p(1, 3).select();
plot(epsilon(:,1),y, 'r', 'Linewidth',1);
xlabel( '$M$');
% ylabel( '$y$');
yticklabels([])
set(gca, 'xtick', [-2:2:2], 'ytick', [-5:2.5:5],'TickDir','out'); 
xlim( [-2,2] );
ylim( [-2.5,2.5] );
grid on;
box on; 



a = findobj(gcf); 

allaxes  = findall(a,'Type','axes');
% alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

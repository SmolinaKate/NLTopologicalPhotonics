%% --------------------------------------------------------------------- %%
clear all;    close all;
addpath('C:\Users\smoli\Dropbox\Topological Nonlinear 2D\split operator method')
map1 = imread('C:\Users\smoli\Dropbox\Topological Nonlinear 2D\Figures_for_paper\paper_Fig_1\Lattice_ver_2.png', 'png','BackgroundColor',[1 1 1]); 
addpath('C:\Users\smoli\Dropbox\Topological Nonlinear 2D\Paraxial_ver_2\Postprocessing\Colormaps\cbrewer')
CT=cbrewer('seq', 'GnBu' , 50, 'pchip');

%% --------------------------------------------------------------------- %%
Num_Save_Tsteps = 30;
Num_Extra_Tsteps = 30;
l_x = 0.01;
N_x = 2^9;
L_x = l_x * N_x;
x = zeros( 1, N_x );
for ix = 1 : N_x
    x( 1, ix ) = ( ix - ( 1 + N_x / 2 ) ) * l_x;
end

l_y = 0.01;
N_y = 2^9;
L_y = l_y * N_y;
y = zeros( 1, N_y );
for iy = 1 : N_y
    y( 1, iy ) = ( iy - ( 1 + N_y / 2 ) ) * l_y;
end
epsilon = zeros(N_y,N_x);
y0=0; 
x0=0;
Nx0=find(x==0);
Nep=find(y==y0);
Nem=find(y==y0);
M=1;
I1=0.5; 
g=0.75;
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
figure('Renderer', 'painters', 'Position', [10 10 600 700])
clf

% create panel
p = panel();
border =1; 
p.pack(2, 5);
p(2,1).repack(0.35)
p(2,5).repack(0.0001);
%p(2,1).select()
%imshow(map1);

p(2, 1).marginright = border; 
p(2, 1).marginleft = border; 
p(2, 2).marginright = border; 
p(2, 2).marginleft = border; 
p(2, 3).marginright = border; 
p(2, 3).marginleft = border; 
p(2, 4).marginright = border; 
p(2, 4).marginleft = border; 


p(2, 2).select();
color_line(squeeze(abs(psiA(:,Nx0))).^2+squeeze(abs(psiA(:,Nx0))).^2,y,squeeze(abs(psiA(:,Nx0)))+squeeze(abs(psiA(:,Nx0))).^2,'LineWidth',2)
colormap(CT);
%set(gca, 'xtick', [0:0.5:1], 'ytick', ...
%    [-2.5:2.5:2.5],'TickDir','in','XAxisLocation',...
%    'bottom','YAxisLocation','right','TickLength',[0.05, 0.01]);
set(gca, 'xtick', [0:0.5:1], 'ytick', ...
    [min(y),0,max(y)],'TickLength',[0.05, 0.01]);
yticklabels([])
xlim( [0,1] );
ylim( [min(x),max(x)] );
ax1 = gca; % current axes
ax1.XColor = [2,62,130]./255;
ax1.YColor = [2,62,130]./255;


ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
hold on; plot(epsilon(:,1),y, '--black', 'Linewidth',1,'Parent',ax2);
set(gca, 'xtick', [min(epsilon(:,1)):1:max(epsilon(:,1))], 'ytick', [],...
    'TickDir','in','XAxisLocation','top','YAxisLocation','right','TickLength',[0.05, 0.01]); 
xlim( [min(epsilon(:,1)),max(epsilon(:,1))]);
ylim( [min(y),max(y)] );




p(2, 3).select();
imagesc( x, y, abs( psiA ).^2+abs( psiB ).^2); 
text(0.7,2.25,'$\mathrm{(c)}$');
colormap(CT);
yticklabels([])
hold off;
box on;
set(gca, 'XAxisLocation','top','TickLength',[0.05, 0.01]);
set(gca, 'XTick', [-2:2:2], 'XTickLabel', {'-2', '0', '2'}, 'fontname','symbol');
set(gca, 'YTick', [-2.5:2.5:2.5],'YTickLabel', {'-2.5', '0', '2.5'},'TickDir','in',...
    'XAxisLocation','top', 'YAxisLocation','right','TickLength',[0.05, 0.01]);
xlim( [-2.5,2.5] );
ylim( [-2.5,2.5] );
a = gca;
b = copyobj(a, gcf);
set(b, 'Layer','top');
box on; 
ax = gca;



text(-4.5,2.8,'$M$')
text(-4.15,-2.9,'$\mathcal{I}$','Color',[2,62,130]./255)
text(-14.3,-0.6,'$a_{0}$')
text(-17.35+0.8,2.1,'$g$')
text(-16.7+0.8,1.22,'$\kappa$')
text(-18.4+0.8,-2.1,'$-M$')
text(-17.6+0.8,-1.7,'$M$')
text(-8.7,0.2,'$x$')
text(0.75,2.85,'$x$')
text(-17.3+0.8,0,'$0$')
text(-17.3+0.8,0.8,'$y$')
text(3,0.8,'$y$')
text(-9.3,2.25,'$\mathrm{(a)}$')
text(-4.6,2.25,'$\mathrm{(b)}$')


a = findobj(gcf); 

allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
size=18;
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(alllines,'LineWidth',2);
str=['Fig_0_ver_2','.png'];
print('-dpng','-r600',str)


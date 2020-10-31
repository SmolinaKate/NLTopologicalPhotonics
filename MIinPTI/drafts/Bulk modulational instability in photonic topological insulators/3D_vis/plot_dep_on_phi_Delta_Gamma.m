
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\Bulk modulational instability in photonic topological insulators\data';
FileName1=fullfile(path, '3D_map_cross.mat');
load(FileName1)
D = incr_cross_band;
h = vol3d('cdata',D,'texture','2D','XData',G,'YData',Delta,'ZData',phi_ar);
view(3);  
%axis tight;  daspect([1 1 .4])
alphamap('rampup');
alphamap(.5.* alphamap);
colormap('hot')
grid on;
box on;



xlabel('$\Delta$')
ylabel('$\Gamma$')
zlabel('$\varphi$')


size=10;
a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(alllines,'LineWidth',1)

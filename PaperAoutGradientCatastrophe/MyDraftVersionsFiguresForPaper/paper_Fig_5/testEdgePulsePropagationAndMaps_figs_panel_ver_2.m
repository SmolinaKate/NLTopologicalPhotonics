clear all;
close all;
setParameters
FontSize = 11;
border=4;
LineWidth=1;
addpath( 'Postprocessing/FigureSettings/' );
addpath( 'Postprocessing/Colormaps/cbrewer/' );
%% --------------------------------------------------------------------- %%
foldername = 'Results/Dynamics/Figs/';
if( ~exist( foldername, 'dir' ) )
   mkdir( foldername );
end
foldernameLC = 'Results/Dynamics/Data/LinearCase/';
foldernameNC = 'Results/Dynamics/Data/NonlinearCase/';
%% --------------------------------------------------------------------- %%
xLC = load( strcat( foldernameLC, 'x.dat' ) );
yLC = load( strcat( foldernameLC, 'y.dat' ) );
maxSqrAbsPsiLC = load( strcat( foldernameLC, 'maxSqrAbsPsi.dat' ) );
Num_Save_Zsteps_LC = length( maxSqrAbsPsiLC( :, 1 ) );
zLC = maxSqrAbsPsiLC( :, 1 );
xNC = load( strcat( foldernameNC, 'x.dat' ) );
yNC = load( strcat( foldernameNC, 'y.dat' ) );
maxSqrAbsPsiNC = load( strcat( foldernameNC, 'maxSqrAbsPsi.dat' ) );
Num_Save_Zsteps_NC = length( maxSqrAbsPsiNC( :, 1 ) );
zNC = maxSqrAbsPsiNC( :, 1 );
%% --------------------------------------------------------------------- %%
%-------------------------------------------------------------------------%
xmin = min( [ xLC( 1, 1 ), xNC( 1, 1 ) ] );
xmax = min( [ xLC( 1, end ), xNC( 1, end ) ] );
ymin = min( [ yLC( 1, 1 ), yNC( 1, 1 ) ] );
ymax = min( [ yLC( 1, end ), yNC( 1, end ) ] );
maxCdata = max( max( [ maxSqrAbsPsiLC( :, 4 ), maxSqrAbsPsiNC( :, 4 ) ] ) );
%-------------------------------------------------------------------------%
addpath( 'Postprocessing/Colormaps/cbrewer/' );
cmap_1 = cbrewer( 'seq', 'Blues', 50, 'pchip' );
cmap_1 = [ [ 1, 1, 1 ]; cmap_1( 1 : end, : ) ];
cmap_1 = brighten( cmap_1, -0.5 );
cmap_2 = cbrewer( 'seq', 'Reds', 50, 'pchip' );
cmap_2 = [ [ 1, 1, 1 ]; cmap_2( 1 : end, : ) ];      
cmap_2 = brighten( cmap_2, -0.5 );
cmap = [ cmap_2; cmap_1 ];
rmpath( 'Postprocessing/Colormaps/cbrewer/' );
%-------------------------------------------------------------------------%
array_z=[1 2];
array_z_dim=[10 20];
%% --------------------------------------------------------------------- %%
%% Figs of evolution. 
%% --------------------------------------------------------------------- %%

figure;
set(gcf, 'Units','centimeters', 'Position',[0 0 8.6*1.5 8.6*1.5])
p = panel();
p.pack(2, 1);
p(1,1).pack(1,2);
p(1,1,1,1).pack(3,1);
p(1,1,1,2).pack(3,1);
p(2,1).pack(1,3);
p(2).repack(0.35);

p(1,1).marginright = border; 
p(1,1).marginleft = border; 
p(1,1).margintop = border; 
p(1,1).marginbottom = border; 
p(1).marginbottom = 2*border; 

for ind_array_z = 1 : 1 : 2
    n=array_z(ind_array_z);
    iz=n;
    
    realPsi = load( strcat( foldernameLC, 'realPsi_z=', num2str( array_z_dim(ind_array_z) ), '.dat' ) );
    imagPsi = load( strcat( foldernameLC, 'imagPsi_z=', num2str( array_z_dim(ind_array_z) ), '.dat' ) );
    psiLC = realPsi + 1i * imagPsi;
    realPsi = load( strcat( foldernameNC, 'realPsi_z=', num2str( array_z_dim(ind_array_z) ), '.dat' ) );
    imagPsi = load( strcat( foldernameNC, 'imagPsi_z=', num2str( array_z_dim(ind_array_z) ), '.dat' ) );
    psiNC = realPsi + 1i * imagPsi;  

    p(1,1,1,ind_array_z,2,1).select();
    hold on;
    [ delta_y, iy ] = min( abs( yLC - maxSqrAbsPsiLC( iz, 3 ) ) );
     plot( xLC, abs( psiLC( iy, : ) ).^2, ...
     'Color', [ 0.75, 0, 0 ] , 'LineWidth', 1.5*LineWidth, 'LineStyle', '-' );
     [ delta_y, iy ] = min( abs( yNC - maxSqrAbsPsiNC( iz, 3 ) ) );
     plot( xNC, abs( psiNC( iy, : ) ).^2, ...
     'Color', [ 51, 171, 249 ] / 255, 'LineWidth', LineWidth, 'LineStyle', '-' );
    hold off;
    box on;
    xlim( [ xmin, xmax ] );
    ylim( [ 0, maxCdata ] );
    if ( ind_array_z ==1) 
        ylabel( '$\mathcal{I}_{wd}$');
    else
        yticklabels([])		
    end
    xticklabels([])	
     
    p_NC=p(1,1,1,ind_array_z,1,1).select();
    strTitle = sprintf( '$z_d= %g$', floor( array_z_dim(ind_array_z)));
    title( strTitle );
    plot( 0, 0, 'k' );
    imagesc( xNC, yNC, maxCdata + abs( psiNC ).^2 );
    
    cmap = [ cmap_2; cmap_1 ];
    colormap(p_NC, cmap)
    
    box on;
    xlim( [ xmin, xmax ] );
    ylim( [ ymin, ymax ] );
    caxis( [ 0, 2 * maxCdata ] );
    xticklabels([])	
    if ( ind_array_z ==1) 
        h=text(-140,0,'$y_d$');
        set(h,'Rotation',90);
    else
        yticklabels([])		
    end
    a = gca;
    b = copyobj(a, gcf);
    set(b, 'Xcolor', [0 0 0], 'YColor', [0 0 0], 'XTickLabel', [], 'YTickLabel', [],'Layer','top');

    p_LC=p(1,1,1,ind_array_z,3,1).select();
    imagesc( xLC, yLC, maxCdata + abs( psiLC ).^2 );
    text(65,-35,'$x_d$')

    cmap = [ cmap_1; cmap_2 ];
    colormap(p_LC, cmap)
    
    box on;
    xlim( [ xmin, xmax ] );
    ylim( [ ymin, ymax ] );
    caxis( [ 0, 2 * maxCdata ] );
    if ( ind_array_z ==1) 
        h=text(-140,0,'$y_d$');
        set(h,'Rotation',90);
    else
        yticklabels([])		
    end
    a = gca;
    b = copyobj(a, gcf);
    set(b, 'Xcolor', [0 0 0], 'YColor', [0 0 0], 'XTickLabel', [], 'YTickLabel', [],'Layer','top');    
end
%% --------------------------------------------------------------------- %%
rmpath( 'Postprocessing/FigureSettings/' );
%% --------------------------------------------------------------------- %%
addpath( 'Dynamics/Data/NonlinearCase' );
max_diff = load( strcat( foldernameNC, 'max_diff', '.dat' ) );
max_diff_norm = load( strcat( foldernameNC, 'max_diff_norm', '.dat' ) );
z_array = load( strcat( foldernameNC, 'z_array', '.dat' ) );
load( strcat( foldernameNC, 'Envelope', '.mat' ) );
amplitudeWP_array = load( strcat( foldernameNC, 'amplitudeWP_array', '.dat' ) );
[z_array_mg,amplitudeWP_mg]=meshgrid(z_array,amplitudeWP_array);
%% --------------------------------------------------------------------- %%
%% Theoretical curve:
z_br=2*M_d^2*sqrt(2.71828182846)/gNL^2./amplitudeWP_mg.^4.*widthWP./VD;
%% --------------------------------------------------------------------- %%
for inda=1:1:length(amplitudeWP_array)
    indz=2;
    while (max_diff_norm(indz,inda)<1.3*max_diff_norm(1,inda)) && (indz<length(z_array))
        indz=indz+1;
    end
    z_br_num(inda)=z_array(indz);
end
p(2,1,1,3).select();
plot(amplitudeWP_array(2:end).^2,z_br_num(2:end),'x','Color',[ 0.75, 0, 0 ]);
hold on; plot(amplitudeWP_array.^2,z_br,'--','Color',[ 0.75, 0, 0 ]);
text(0.3,0,'$F_{0d}$')
text(0.16,32,'$z_{*d}$')
box on;
MakeFigNice
print( '-dpng', 'Fig_5_d' , '-r600', '-opengl' );

function Visualization(psiA, psiB, x, y, Ne, Nx0, SolitonTrajectory)
%-------------------------------------------------------------------------%
xmin = x( 1, 1 );    xmax = x( 1, end );
ymin = y( 1, 1 );    ymax = y( 1, end );
%-------------------------------------------------------------------------%

strTitle = sprintf( 't= %g', SolitonTrajectory( end, 1 ) );

subplot( 1, 3, 3 );
hold on;
plot(x, abs( psiA(Ne,:)).^2,'DisplayName',[strTitle]);
hold off;
box on;

xlim( [ xmin, xmax ] ); 
xlabel( 'x', 'FontSize', 11, 'FontName', 'Times New Roman', ...
             'FontWeight', 'normal', 'FontAngle', 'normal' );
ylabel( '$|\psi_1(y=y_0)|^2$', 'FontSize', 11, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex' );

legend(gca,'show')
box on;


subplot( 1, 3, 2 );
hold on;
imagesc( x, y, abs( psiA ).^2 + abs( psiB ).^2 );
hold off;
box on;
daspect( [ 1, 1, 1 ] );
xlim( [ xmin, xmax ] );    ylim( [ ymin, ymax ] );
xlabel( 'x', 'FontSize', 11, 'FontName', 'Times New Roman', ...
             'FontWeight', 'normal', 'FontAngle', 'normal' );
ylabel( 'y', 'FontSize', 11, 'FontName', 'Times New Roman', ...
             'FontWeight', 'normal', 'FontAngle', 'normal' );
set( gca, 'FontSize', 11, 'FontName', 'Times New Roman', ...
          'FontWeight', 'normal', 'FontAngle', 'normal', ...
          'LineWidth', 1, 'XDir', 'normal', 'YDir', 'normal', ...
          'XMinorTick', 'on', 'YMinorTick', 'on' );
subplot( 1, 3, 1 );
hold on;
plot(y, abs( psiA(:,Nx0)).^2);
hold off;
xlim( [ ymin, ymax ] ); 
xlabel( 'y', 'FontSize', 11, 'FontName', 'Times New Roman', ...
             'FontWeight', 'normal', 'FontAngle', 'normal' );
% ylabel( '$S_x$', 'FontSize', 11, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex' );
ylabel( '$|\psi_1(x=x_0)|^2$', 'FontSize', 11, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex' );
box on;


% title( strTitle, 'FontSize', 11, 'FontName', 'Times New Roman', ...
%                  'FontWeight', 'normal', 'FontAngle', 'normal' );
% legend (strTitle)


% %-------------------------------------------------------------------------%
 drawnow;
%  print( '-dpng', length( SolitonTrajectory( :, 1 ) ) - 1 );
%-------------------------------------------------------------------------%
end
clear all;
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster\profile_evolution_data';
file_name = sprintf('initial_parameters');
ToSaveFileName = [file_name '.mat'];
ToSaveFileName1 = fullfile(path, ToSaveFileName);
load(ToSaveFileName1);                    

Num_of_fig = Num_Save_Tsteps + 1;
figure;
subplot (2,Num_of_fig,Num_of_fig + 1)
hold on;
plot(x,abs( psiA(Nem,:)).^2,'r', 'Linewidth',1);
ylabel( '$|\psi_1(y=y_0)|^2$');
xlabel( '$x$');
        
grid on; 
box on;
         
subplot (2,Num_of_fig,1)
hold on;
imagesc( x, y, abs( psiA).^2); colormap hot;
hold off;
xlabel( '$x$');
ylabel( '$y$' );
        
box on;        
        
        
t=4*M^2*sqrt(2.71828)/g^2./sqrt(B)./A^2        
t=4*M^2*sqrt(2.71828)/g^2./sqrt(B)./A^2 *0.7 

for it = 1 : Num_Save_Tsteps
     file_name = sprintf('pr_ev#%it',it);
     ToSaveFileName = [file_name '.mat'];
     ToSaveFileName1 = fullfile(path, ToSaveFileName);
     load(ToSaveFileName1,'psiA','psiB');                    
      
        
        
        subplot (2,Num_of_fig,Num_of_fig + it + 1)
        hold on;
        plot(x,abs( psiA(Nem,:)).^2,'r', 'Linewidth',1);
        xlabel( '$x$');

        grid on; 
        box on;

        subplot (2,Num_of_fig,1 + it)
        hold on;
        imagesc( x, y, abs( psiA).^2); colormap hot;
        hold off;
        xlabel( '$x$' );

        box on;        

end


a = findobj(gcf); 
size = 12;

allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alllines,'LineWidth',1);
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

        
        
        
        
M=0.1;
g=M/3;
A=1;
B=[0.0006 0.0007 0.0008 0.0009 0.001];

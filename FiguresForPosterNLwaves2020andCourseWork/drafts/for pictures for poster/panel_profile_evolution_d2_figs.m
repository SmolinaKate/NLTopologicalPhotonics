clear all;
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster\profile_evolution_with_d2_data';
file_name = sprintf('initial_parameters');
ToSaveFileName = [file_name '.mat'];
ToSaveFileName1 = fullfile(path, ToSaveFileName);
load(ToSaveFileName1);                    

Num_of_fig = Num_Save_Tsteps + 1;
%%параметры границ и осей 
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);
ymin1 = 0;
ymax1 = 2;

border = 1;       
ar_xticks = [floor(xmin/100)*100:100:floor(xmax/100)*100];
ar_yticks = [-50:50:50];
ar_yticks2 = [ymin1:0.5:ymax1];

%%сдвиг для непрерывной картинки по опрокидыванию
step_max=40;


% clf
r=250;
figure('Renderer', 'painters', 'Position', [10 10 Num_of_fig*r r])
clf

% create panel
p = panel();

Num_of_fig = Num_Save_Tsteps + 1;

p.pack(2, Num_of_fig);

p(2, 1).select();
plot(x,abs( psiA(Nem,:)).^2,'r', 'Linewidth',1);
xlim( [ xmin, xmax ] );
ylim( [ ymin1, ymax1 ] );
p(2, 1).marginright = border; 
p(2, 1).marginleft = border; 
set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks2,'TickDir','out'); 
xticklabels([-100:100:100])		
grid on; 
box on;


p(1).marginbottom = border; 
p(2).margintop = border; 


p(1, 1).select();
hold on;
imagesc( x, y, abs( psiA).^2); colormap hot;
xlim( [ xmin, xmax ] );
ylim( [ ymin, ymax ] );
hold off;
p(1, 1).marginright = border; 
p(1, 1).marginleft = border; 
box on;
set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks,'TickDir','out');
xticklabels([])		


% p(1, Num_Save_Tsteps+1).select();
p(1, 1).select();
p(1, 1).ylabel( '$y$');

% p(2, Num_Save_Tsteps+1).select();
p(2, 1).select();
p(2, 1).ylabel( '$|\psi_1(y=y_0)|^2$');

t=4*M^2*sqrt(2.71828)/g^2./sqrt(B)./A^2        
t=4*M^2*sqrt(2.71828)/g^2./sqrt(B)./A^2 *0.7 

for it = 1 : Num_Save_Tsteps
     file_name = sprintf('pr_ev#%it',it);
     ToSaveFileName = [file_name '.mat'];
     ToSaveFileName1 = fullfile(path, ToSaveFileName);
     load(ToSaveFileName1,'psiA','psiB');                    
      
       [M,N] = max(abs(psiA(Nem,:)));
       [Ny,Nx] = size(psiA);
        for ix = 1:1:Ny 
             psiA_new(ix,:)= circshift(psiA(ix,:), - N + ceil(Nx/2) + step_max*it);
        end

      
        p(1, 1 + it).select();
        imagesc( x, y, abs( psiA_new).^2); colormap hot;
        xlim( [ xmin, xmax ] );
        ylim( [ ymin, ymax ] );
        set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks,'TickDir','out');
        xticklabels([])	
        yticklabels([])	
        p(1, 1 + it).marginright = border; 
        p(1, 1 + it).marginleft = border; 
        hold off;
        box on;   
        
        p(2, it + 1).select();
        hold on;
        plot(x,abs( psiA_new(Nem,:)).^2,'r', 'Linewidth',1);
        xlim( [ xmin, xmax ] );
        ylim( [ ymin1, ymax1 ] );
        set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks2,'TickDir','out');
        xticklabels([])	
        yticklabels([])	
        p(2, it + 1).marginright = border; 
        p(2, it + 1).marginleft = border; 
        grid on; 
        box on;
end


xlabel( '$x$');


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



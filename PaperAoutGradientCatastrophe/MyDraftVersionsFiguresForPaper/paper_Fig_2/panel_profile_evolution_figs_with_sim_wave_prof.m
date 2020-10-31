clear all; close all;

path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\Figures_for_paper\paper_Fig_2\profile_evolution_data';
file_name = sprintf('initial_parameters');
ToSaveFileName = [file_name '.mat'];
ToSaveFileName1 = fullfile(path, ToSaveFileName);
load(ToSaveFileName1);                    

%%параметры границ и осей 
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);
ymin1 = 0;
ymax1 = 2;

y1=40;
y2=88;

border = 1;       
ar_xticks = [floor(xmin/100)*100:100:floor(xmax/100)*100];
ar_yticks = [-50:50:50];
ar_yticks2 = [ymin1:1:ymax1];

%%сдвиг для непрерывной картинки по опрокидыванию
step_max=1;
Num_of_fig = Num_Save_Tsteps + 1;


% clf
r=250;
figure('Renderer', 'painters', 'Position', [10 10 Num_of_fig*r/1.5/2.7 r*1.5])
clf

% create panel
p = panel();


p.pack(4, 4);

p(2, 1).select();
plot(x,abs( psiA(Nem,:)).^2,'r', 'Linewidth',1);
xlim( [ xmin, xmax ] );
ylim( [ ymin1, ymax1 ] );
p(2, 1).marginright = border; 
p(2, 1).marginleft = border; 
set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks2,'TickDir','out'); 
xticklabels([0:1:3])		
grid on; 
box on;


p(1).marginbottom = border; 
p(2).margintop = border; 


p(3).marginbottom = border; 
p(4).margintop = border; 

p(1, 1).select();
hold on;
imagesc( x, y(y1:y2), abs( psiA(y1:y2,:)).^2); colormap hot; axis equal;  
xlim( [ xmin, xmax ] );
ylim( [ y(y1), y(y2) ] );
hold off;
p(1, 1).marginright = border; 
p(1, 1).marginleft = border; 
box on;
set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks,'TickDir','out');
xticklabels([])		

p(1, 1).select();
p(1, 1).ylabel( '$y$');

p(2, 1).select();
p(2, 1).ylabel( '$|\psi_1(y=0)|^2$');

p(3, 1).select();
p(3, 1).ylabel( '$y$');

p(4, 1).select();
p(4, 1).ylabel( '$|\psi_1(y=0)|^2$');

[Ip,treshhold,x_sim_wave]=SimWaveChar(dt,Num_Extra_Tsteps,Num_Save_Tsteps,l_x,A,B,g,M);


for it = 1 : Num_Save_Tsteps
    if it<=3 %%строю первую группу из четырех картинок
       file_name = sprintf('pr_ev#%it',it);
       ToSaveFileName = [file_name '.mat'];
       ToSaveFileName1 = fullfile(path, ToSaveFileName);
       load(ToSaveFileName1,'psiA','psiB'); 
       
       [M,N] = max(abs(psiA(Nem,:)));
       [Ny,Nx] = size(psiA);
        for ix = 1:1:Ny 
             %psiA_1(ix,:)= circshift(psiA(ix,:), - N + ceil(Nx/2) + step_max*it);
             psiA_1(ix,:)= circshift(psiA(ix,:), - N + ceil(Nx/2));
             psiA_new(ix,:)=flip(psiA_1(ix,:));
        end
        
        p(1, 1 + it).select();
        imagesc( x, y(y1:y2), abs( psiA_new(y1:y2,:)).^2);colormap hot; axis equal;  
        xlim( [ xmin, xmax ] );
        ylim( [ y(y1), y(y2)] );
        set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks,'TickDir','out');
        xticklabels([])	
        yticklabels([])	
        p(1, 1 + it).marginright = border; 
        p(1, 1 + it).marginleft = border; 
        hold off;
        box on;   
        
        p(2, it + 1).select();
        hold on;
        plot(x,abs( psiA_new(Nem,:)).^2./max(abs( psiA_new(Nem,:)).^2),'r', 'Linewidth',1);
        xlim( [ xmin, xmax ] );
        ylim( [ ymin1, ymax1 ] );
        set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks2,'TickDir','out');
        xticklabels([(it+it*2):1:(it+it*2+3)])	
        yticklabels([])	
        p(2, it + 1).marginright = border; 
        p(2, it + 1).marginleft = border; 
        grid on; 
        box on;
        hold on;
        
        s_new = SimWaveProf(Ip,Num_Extra_Tsteps,treshhold,it,x_sim_wave,Nx);
        plot(x,(s_new)./max(s_new),'--');
        hold off;
    end
    
    
    
    
    if (it>3)&&(it<=5) %%строю вторую группу из четырех картинок -- которые снизу и для которых надо строить простую волну
       file_name = sprintf('pr_ev#%it',it);
       ToSaveFileName = [file_name '.mat'];
       ToSaveFileName1 = fullfile(path, ToSaveFileName);
       load(ToSaveFileName1,'psiA','psiB'); 
       
       [M,N] = max(abs(psiA(Nem,:)));
       [Ny,Nx] = size(psiA);
        for ix = 1:1:Ny 
             %psiA_1(ix,:)= circshift(psiA(ix,:), - N + ceil(Nx/2) + step_max*it);
             psiA_1(ix,:)= circshift(psiA(ix,:), - N + ceil(Nx/2));
             psiA_new(ix,:)=flip(psiA_1(ix,:));
        end
        
        p(3, -3 + it).select();
        imagesc( x, y(y1:y2), abs( psiA_new(y1:y2,:)).^2);colormap hot; axis equal;  
        xlim( [ xmin, xmax ] );
        ylim( [ y(y1), y(y2)] );
        set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks,'TickDir','out');
        xticklabels([])	
        yticklabels([])	
        p(3, -3 + it).marginright = border; 
        p(3, -3 + it).marginleft = border; 
        hold off;
        box on;   
        
        p(4, it - 3).select();
        hold on;
        plot(x,abs( psiA_new(Nem,:)).^2./max(abs( psiA_new(Nem,:)).^2),'r', 'Linewidth',1);
        xlim( [ xmin, xmax ] );
        ylim( [ ymin1, ymax1 ] );
        set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks2,'TickDir','out');
        xticklabels([(it+it*2):1:(it+it*2+3)])	
        yticklabels([])	
        p(4, -3 + it).marginright = border; 
        p(4, -3 + it).marginleft = border; 
        grid on; 
        box on;
        hold on;
        
        s_new = SimWaveProf(Ip,Num_Extra_Tsteps,treshhold,it,x_sim_wave,Nx);
        plot(x,(s_new)./max(s_new),'--');
        hold off;
    end
    
    
    if (it>6)%%строю картинки без простой волны
       file_name = sprintf('pr_ev#%it',it);
       ToSaveFileName = [file_name '.mat'];
       ToSaveFileName1 = fullfile(path, ToSaveFileName);
       load(ToSaveFileName1,'psiA','psiB'); 
       
       [M,N] = max(abs(psiA(Nem,:)));
       [Ny,Nx] = size(psiA);
        for ix = 1:1:Ny 
             %psiA_1(ix,:)= circshift(psiA(ix,:), - N + ceil(Nx/2) + step_max*it);
             psiA_1(ix,:)= circshift(psiA(ix,:), - N + ceil(Nx/2));
             psiA_new(ix,:)=flip(psiA_1(ix,:));
        end
        
        p(3, -4 + it).select();
        imagesc( x, y(y1:y2), abs( psiA_new(y1:y2,:)).^2);colormap hot; axis equal;  
        xlim( [ xmin, xmax ] );
        ylim( [ y(y1), y(y2)] );
        set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks,'TickDir','out');
        xticklabels([])	
        yticklabels([])	
        p(3, -4 + it).marginright = border; 
        p(3, -4 + it).marginleft = border; 
        hold off;
        box on;   
        
        p(4, it - 4).select();
        hold on;
        plot(x,abs( psiA_new(Nem,:)).^2,'r', 'Linewidth',1);
        xlim( [ xmin, xmax ] );
        ylim( [ ymin1, ymax1 ] );
        set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks2,'TickDir','out');
        xticklabels([(it+it*2-3):1:(it+it*2+3-3)])	
        yticklabels([])	
        p(4, -4 + it).marginright = border; 
        p(4, -4 + it).marginleft = border; 
        grid on; 
        box on;
    end
    
end


text(115,-0.1,'$10^2 x$')
text(-1000,7.5,'$t=0$')
text(-700,7.5,'$t=240$')
text(-400,7.5,'$t=480$')
text(-100,7.5,'$t=720$')
text(-1000,1.4,'$t=960$')
text(-700,1.4,'$t=1200$')
text(-400,1.4,'$t=1680$')
text(-100,1.4,'$t=1920$')


%%параметры текста и толщин линий
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




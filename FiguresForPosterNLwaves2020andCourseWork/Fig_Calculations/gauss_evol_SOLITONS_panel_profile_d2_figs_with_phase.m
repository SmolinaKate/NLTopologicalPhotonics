%%%здесь строятся картинки с второй производной -- где массив из 6
%%%картинок, показывающий формирование солитона, и его фазу



clear all;
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster\profile_evolution_with_d2_data';
file_name = sprintf('initial_parameters');
ToSaveFileName = [file_name '.mat'];
ToSaveFileName1 = fullfile(path, ToSaveFileName);
load(ToSaveFileName1);                    


%%сдвиг для непрерывной картинки по опрокидыванию
step_max=40;
Num_of_fig = 6;


% clf
r=250;
figure('Renderer', 'painters', 'Position', [0 0 Num_of_fig*r r])
clf

% create panel
p = panel();


p.pack(4, Num_of_fig);

p(2, 1).select();
plot(x,abs( psiA(Nem,:)).^2,'r', 'Linewidth',1);
grid on; 
box on;


%txt='$5\cdot10^{-2}$';

%text(-175,5.5*10^(-2),txt)


p(1, 1).select();
hold on;
imagesc( x, y, abs( psiA(:,:)).^2); colormap hot; 
hold off;
box on;



p(1, 1).select();
p(1, 1).ylabel( '$y$');

% p(2, Num_Save_Tsteps+1).select();
p(2, 1).select();
p(2, 1).ylabel( '$|\psi_1(y=y_0)|^2$');


ar_for_fig = [4,5,6,14,47];

for it = 1 : (Num_of_fig - 1)
       j = ar_for_fig(it);
       file_name = sprintf('pr_ev#%it',j);
       ToSaveFileName = [file_name '.mat'];
       ToSaveFileName1 = fullfile(path, ToSaveFileName);
       load(ToSaveFileName1,'psiA','psiB'); 
       
       [M,N] = max(abs(psiA(Nem,:)));
       [Ny,Nx] = size(psiA);
        for ix = 1:1:Ny 
             psiA_new(ix,:)= circshift(psiA(ix,:), - N + ceil(Nx/2) + step_max*it);
        end

      
        p(1, 1 + it).select();
        imagesc( x, y, abs( psiA_new(:,:)).^2);colormap hot;  
    
        hold off;
        box on;   
        
        p(2, it + 1).select();
        hold on;
        plot(x,abs( psiA_new(Nem,:)).^2,'r', 'Linewidth',1);

        grid on; 
        box on;

end

Ylm=ylim;                         
Xlm=xlim;                          
Xlb=max(Xlm);                    
Ylb=- Ylm(1); 
%hXLbl=xlabel('$10^2 x$','Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center'); 
txt='$10^2 x$';

text(125,-0.52*10^(-2)*2,txt)
%%%создание подписи справа























%%%%%%для фазлвый картинок

file_name = sprintf('initial_parameters');
ToSaveFileName = [file_name '.mat'];
ToSaveFileName1 = fullfile(path, ToSaveFileName);
load(ToSaveFileName1);                    
p(4, 1).select();
plot(x,angle( psiA(Nem,:)).^2,'r', 'Linewidth',1);

grid on; 
box on;


%txt='$5\cdot10^{-2}$';

%text(-175,5.5*10^(-2),txt)


p(3, 1).select();
hold on;
imagesc( x, y, angle( psiA(:,:)).^2); colormap hot; 
hold off;

box on;



p(3, 1).select();
p(3, 1).ylabel( '$y$');

% p(2, Num_Save_Tsteps+1).select();
p(4, 1).select();
p(4, 1).ylabel( '$|\psi_1(y=y_0)|^2$');


ar_for_fig = [4,5,6,14,47];

for it = 1 : (Num_of_fig - 1)
       j = ar_for_fig(it);
       file_name = sprintf('pr_ev#%it',j);
       ToSaveFileName = [file_name '.mat'];
       ToSaveFileName1 = fullfile(path, ToSaveFileName);
       load(ToSaveFileName1,'psiA','psiB'); 
       
       [M,N] = max(abs(psiA(Nem,:)));
       [Ny,Nx] = size(psiA);
        for ix = 1:1:Ny 
             psiA_new(ix,:)= circshift(psiA(ix,:), - N + ceil(Nx/2) + step_max*it);
        end

      
        p(3, 1 + it).select();
        imagesc( x, y, angle( psiA_new(:,:)).^2);colormap hot; 

        hold off;
        box on;   
        
        p(4, it + 1).select();
        hold on;
        plot(x,angle( psiA_new(Nem,:)).^2,'r', 'Linewidth',1);

        grid on; 
        box on;

end








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




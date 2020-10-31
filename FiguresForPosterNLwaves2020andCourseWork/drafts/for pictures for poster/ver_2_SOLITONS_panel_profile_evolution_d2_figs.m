%%%здесь строятся картинки с второй производной -- где массив из 6 картинок, показывающий формирование солитона, которые потом
%%%будут использованы в курсовой! (и в постере на нелинейные волны тоже они были, только нужно внимательно смотреть на параметры)



clear all;
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster\profile_evolution_with_d2_data';
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
ymax1 = 50*10^(-3);

y1=40;
y2=88;

border = 1;       
ar_xticks = [floor(xmin/100)*100:100:floor(xmax/100)*100];
ar_yticks = [-50:50:50];
ar_yticks2 = [ymin1:(ymax1 - ymin1)/1:ymax1];

%%сдвиг для непрерывной картинки по опрокидыванию
step_max=40;
Num_of_fig = 6;


% clf
r=250;
figure('Renderer', 'painters', 'Position', [0 0 Num_of_fig*r r])
clf

% create panel
p = panel();


p.pack(2, Num_of_fig);

p(2, 1).select();
plot(x,abs( psiA(Nem,:)).^2,'r', 'Linewidth',1);
xlim( [ xmin, xmax ] );
ylim( [ ymin1, ymax1 ] );
p(2, 1).marginright = border; 
p(2, 1).marginleft = border; 
set(gca, 'xtick', ar_xticks, 'ytick', [0:0.05:0.05],'TickDir','out'); 
xticklabels([0:1:3])		
grid on; 
box on;


%txt='$5\cdot10^{-2}$';

%text(-175,5.5*10^(-2),txt)

p(1).marginbottom = border; 
p(2).margintop = border; 


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
        plot(x,abs( psiA_new(Nem,:)).^2,'r', 'Linewidth',1);
        xlim( [ xmin, xmax ] );
        ylim( [ ymin1, ymax1 ] );
        set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks2,'TickDir','out');
        xticklabels([(it+it*2):1:(it+it*2+3)])	
        yticklabels([])	
        p(2, it + 1).marginright = border; 
        p(2, it + 1).marginleft = border; 
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

% p(1, Num_of_fig).select();
% set(gca, 'Box', 'on');  % Turn off the box surrounding the whole axes
% axesPosition = get(gca, 'Position');           % Get the current axes position
% hNewAxes = axes('Position', axesPosition, ...  % Place a new axes on top...
%                 'Color', 'none', ...           %   ... with no background color
%                 'YLim', [ar_yticks(1),ar_yticks(end)], ...            %   ... and a different scale
%                 'YAxisLocation', 'right', ...  %   ... located on the right
%                 'XTick', [], ...               %   ... with no x tick marks
%                 'Box', 'off');                 %   ... and no surrounding box


% 
% p(2, Num_of_fig).select();
% set(gca, 'Box', 'on');  % Turn off the box surrounding the whole axes
% axesPosition = get(gca, 'Position');           % Get the current axes position
% hNewAxes = axes('Position', axesPosition, ...  % Place a new axes on top...
%                 'Color', 'none', ...           %   ... with no background color
%                 'YLim', [ar_yticks2(1),ar_yticks2(end)], ...            %   ... and a different scale
%                 'YAxisLocation', 'right', ...  %   ... located on the right
%                 'XTick', [], ...               %   ... with no x tick marks
%                 'Box', 'off');                 %   ... and no surrounding box
% 




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




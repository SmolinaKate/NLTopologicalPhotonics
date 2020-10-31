clear all;
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster\profile_evolution_data';
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
figure('Renderer', 'painters', 'Position', [10 10 Num_of_fig*r r])
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
set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks2,'TickDir','out'); 
xticklabels([0:1:3])		
grid on; 
box on;


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




%%%ПРОСТАЯ ВОЛНА


%%для построения нужна больщая сетка по икс, чтобы характеристики не вышли
%%за ее пределы
t=0:dt:Num_Extra_Tsteps*Num_Save_Tsteps*dt;
dx=l_x;
N_x_sim_wave  = 2^14;%2^12 -- работает
L_x = l_x * N_x_sim_wave;
x_sim_wave = zeros( 1, N_x_sim_wave );
for ix = 1 : N_x_sim_wave
    x_sim_wave ( 1, ix ) = ( ix - ( 1 + N_x_sim_wave  / 2 ) ) * l_x;
end

%%построение характеристик для простой волны 
treshhold=10^(-10);
I_0 = (A*exp ( - B*(x_sim_wave-0).^2) + treshhold);
v_g = 1 - g^2*I_0.^2/4/M^2;
for i=1:1:length(v_g)
xP(i,:) = x_sim_wave(i) + v_g (i).*t;
end

% t=2*M^2*sqrt(2.71828)/g^2./sqrt(B)./A^2;


Ip = zeros(length(I_0), length(t));
for it=2:1:length(t)
    for j=1:1:length(I_0)
         if round(xP(j,it)/dx) > 0 && round((xP(j,it)-x_sim_wave(1))/dx) <= length(I_0)
            Ip(round((xP(j,it)-x_sim_wave(1))/dx),it)=I_0(j);
         end
    end
end


% Ip = zeros(length(I_0), length(t));
% Ip(:,1)=I_0;
% 
% for it=2:1:length(t)
%      for j=1:1:length(I_0)
%         Ip(j,it) = Ip(,it-1);
%       end
% end


for it = 1 : Num_Save_Tsteps
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
        
        
        %if(it<3)
        %%%простая волна
        x_approx = [];
        r=1;
        for j=1:1:length(Ip(:,1))
                if (abs(Ip(j,it*Num_Extra_Tsteps)))>=treshhold
                   x_approx(r)=x_sim_wave(j);
                   Prof_br(r)=Ip(j,it*Num_Extra_Tsteps);
                   r=r+1;
                end
        end         
        
        vq = interp1(x_approx,Prof_br,x_sim_wave);
        [a,b]=max(vq);
        s=vq(b-floor(Nx/2)+1:b+floor(Nx/2));
        [M,N] = max(abs(s));
        [Ny,Nx] = size(s);
        %if(it==2)
         %  s(262)=1.99;
         %  s(265)=1.98;
        %end
        s_new(:)= circshift(s(:), - N + ceil(Nx/2));
        plot(x,(s_new)./max(s_new),'--');
        hold off;
        %end
        clear x_approx r Prof_br
end


Ylm=ylim;                         
Xlm=xlim;                          
Xlb=max(Xlm);                    
Ylb=- Ylm(1); 
%hXLbl=xlabel('$10^2 x$','Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center'); 

txt = {'$10^2 x$'};
text(125,-0.43,txt)



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




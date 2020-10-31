%%%%%%%%%%
%% 1) PERF_PAR
%% 2) B=B/4
%% 3) B=B/4 and M = M/2
%% 4) A=A/10
%% 5) A=5A
%% 6) B = B/30
%% 7) A = 10A
%% 8) g = g/2

%%%%%%%%%%

clear all; close all

path='D:\Users\smoli\Desktop\test\sol_NSE';
file_name = sprintf('initial_parameters');
ToSaveFileName = [file_name '.mat'];
ToSaveFileName1 = fullfile(path, ToSaveFileName);
load(ToSaveFileName1); 







M0=M;
% figure;
% plot(y,abs( psiA(:,512/2)).^2,'r');
% hold on; plot(y,exp(- 2.*M0.*abs(y)).*A,'--b');
% xlabel('$y$')
% ylabel('$f$')
% legend('|\psi_1(x=x_0)|^2','A e^{-2 M_0 |y|}')
% box on;

% a = findobj(gcf); 
% size1 = 12;
% allaxes  = findall(a,'Type','axes');
% alllines = findall(a,'Type','line');
% alltext  = findall(a,'Type','text');
% set(allaxes,'FontSize', size1, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal');
% set(alllines,'LineWidth',1);
% set(alltext,'FontSize', size1, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        




Num_of_fig = Num_Save_Tsteps + 1;
%%параметры границ и осей 
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

border = 10;       
ar_xticks = [floor(xmin/100)*100:100:floor(xmax/100)*100];
ar_yticks = [-50:50:50];

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
%xlim( [ xmin, xmax ] );
%ylim( [ ymin1, ymax1 ] );
p(2, 1).marginright = border; 
p(2, 1).marginleft = border; 
% set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks2,'TickDir','out'); 
xticklabels([])		
grid on; 
box on;


p(1).marginbottom = border; 
p(2).margintop = border; 


p(1, 1).select();
hold on;
imagesc( x, y, abs( psiA).^2); colormap hot;
daspect( [ 1, 1, 1 ] );
xlim( [ xmin, xmax ] );
ylim( [ ymin, ymax ] );
hold off;
p(1, 1).marginright = border; 
p(1, 1).marginleft = border; 
box on;
set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks,'TickDir','out');
xticklabels([])	

hold off;
figure(2)
a=abs( psiA(:,:)).^2;
contour(a);
daspect( [ 1, 1, 1 ] );
hold off;



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
        daspect( [ 1, 1, 1 ] );
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
        plot(x,abs( psiA_new(Nem,:)).^2,'r');
    %    xlim( [ xmin, xmax ] );
     %   ylim( [ ymin1, ymax1 ] );
%         set(gca, 'xtick', ar_xticks, 'ytick', ar_yticks2,'TickDir','out');
        xticklabels([])	
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

 hold off;
 a=abs( psiA_new(:,:)).^2;
 figure(2); 
 hold on;
 contour(a)
 

 dot = zeros (128,512);
 dot (128/2,3*512/4) = abs(max(max(psiA_new(:,:))));
 
 hold on;
 contour(dot) 
 
 
% link=1;
% a=abs(2/link/sqrt(3));
% eta = link*a^2/8;
%  
% [m,n] = max(abs( psiA(128/2,:)).^2); 
% figure;
% plot(y,abs( psiA(:,n)).^2,'r');
% hold on; plot(y,exp(- 2.*M0.*abs(y)).*A,'--b');
% xlabel('$y$')
% ylabel('$f$')
% legend('|\psi_1(x=x_0)|^2','A e^{-2 M_0 |y|}')
% box on;
% 
% a = findobj(gcf); 
% size1 = 12;
% allaxes  = findall(a,'Type','axes');
% alllines = findall(a,'Type','line');
% alltext  = findall(a,'Type','text');
% set(allaxes,'FontSize', size1, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal');
% set(alllines,'LineWidth',1);
% set(alltext,'FontSize', size1, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
% 
% 
%  
%  
 
 
 
 
 
 
 
 
%  
% link=1;
% a=abs(2/link/sqrt(3));
% % a=0.1;
% Kp=4*pi/3/a;
% Kp1=4*pi/3/a;
% a0=sqrt(3)*a;
% dt =2;
% eta = link*a^2/8;
%  
% g=eta;
% B=eta^2;
% A=sqrt(B);
% M = link/10;
% v=1;
% omega = -0.0004;
% beta = 2*g*eta*M^2;
% % omega = - beta*A^2*eta/2 - 1/eta*(v/2 - 2*M*eta^2)^2 - 4 *M*eta^2*(v/2 - 2*M*eta^2)^2; 
% omegao= - (omega + 1/eta*(v/2 - 2*M*eta^2)^2 + 4 *M*eta^2*(v/2 - 2*M*eta^2)^2)/eta;
% 
% %f = sqrt(omegao*2/beta)./cosh(sqrt(omegao).*x);
% f = A./cosh(sqrt(B).*x);
% 
% figure(4);  
% plot(x,abs( psiA_new(Nem,:)).^2,'r')
% hold on; plot(x,f);

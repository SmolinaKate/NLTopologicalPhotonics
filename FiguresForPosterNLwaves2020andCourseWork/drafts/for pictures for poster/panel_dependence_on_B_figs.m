clear all;
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster\dependence_on_B_data';

% шаг из моделирования
dt = 1;
%пороговое значение производной (в это число раз превышает начальное значение)
cr_val_der = 3;

M=0.1;
g=M/3;
A=1;
B=[0.0006 0.0007 0.0008 0.0009 0.001];


t=2*M^2*sqrt(2.71828)/g^2./sqrt(B)./A^2;
t1=ceil(t(1))-1;
t2=ceil(t(2));
t3=ceil(t(3))-1;
t4=ceil(t(4));
t5=ceil(t(5))-1;



r=250;
figure('Renderer', 'painters', 'Position', 1.5*[10 10 2*r r])
clf

% create panel
p = panel();
border = 20; 

p.pack(1, 2);


p(1, 2).marginright = border; 
p(1, 2).marginleft = border; 
p(1, 1).marginright = border; 
p(1, 1).marginleft = border; 


p(1, 1).select();


i=0;
ee=1400;
for B=0.0006:0.0001:0.001 
%загрузка распределения производной для каждого из значений параметра B             
            i = i+1;
            file_name = sprintf('values_der_B_#%d',B);
            ToSaveFileName = [file_name '.mat'];
            ToSaveFileName1 = fullfile(path, ToSaveFileName);
            load(ToSaveFileName1);
            
            
        der1 = der;
        tau=dt:dt:dt*(length(der1)-1);
%макс значение производной по у        
        level(1:length(der1)-1)=der1(1);

        [m,j] = min(der1);
        to(i) = tau(j);
        
%время, определенное из cr_val_der        
        c=3;
        while der1(c)<= cr_val_der*der1(2)
            c = c+1;
        end
        tcr(i) = tau(c);
        
        
        if (i==1)
%           hold on;  plot(tau,level,'r')
%           str{i} = ['\Delta_y(t=0)'];
            hold on; plot(tau(2:ee),der1(2:ee),'b')
            hold on; plot([t1 t1],[0 0.14],'--b','HandleVisibility','off')
            str{i} = ['B = 6 \cdot 10^{-4}']; 
            
        end
        if (i==2)
            hold on; plot(tau(2:ee),der1(2:ee),'c')
            hold on; plot([t2 t2],[0 0.14],'--c','HandleVisibility','off')
            str{i} = ['B = 7 \cdot 10^{-4}']; 
        end
        if (i==3)
            hold on; plot(tau(2:ee),der1(2:ee),'m')
            hold on; plot([t3 t3],[0 0.14],'--m','HandleVisibility','off')
            str{i} = ['B = 8 \cdot 10^{-4}']; 
        end
        if (i==4)
            hold on; plot(tau(2:ee),der1(2:ee),'black')
            hold on; plot([t4 t4],[0 0.14],'--black','HandleVisibility','off')
            str{i} = ['B = 9 \cdot 10^{-4}']; 
        end
        if (i==5)
            hold on; plot(tau(2:ee),der1(2:ee),'green')
            hold on; plot([t5 t5],[0 0.14],'--green','HandleVisibility','off')
            str{i} = ['B = 10 \cdot 10^{-4}']; 
        end
end
legend(str,'Location','northwest')
grid on;
ylabel( '$\cdot 10^{-2} \Delta$');
xlabel( '$\cdot 10^2, t$');
set(gca, 'xtick', [0:200:1400], 'ytick', [0:0.03:0.14],'TickDir','out'); 
xticklabels([0:2:14])
yticklabels([0:3:14])
xlim( [ 0, 1400 ] );
ylim( [ 0, 0.14] );
box on


            


B=[0.0006 0.0007 0.0008 0.0009 0.001];



p(1, 2).select();

plot(B,t,'b')
hold on; plot(B,tcr,'or','LineWidth',1)
hold on; plot(B,to,'oc','LineWidth',1)
hold on; plot(B, 0.76.*t,'--b','LineWidth',1)
hold on; plot(B, 1.34.*t,'--b','LineWidth',1)

box on;
grid on;
ylabel('$\cdot 10^{2}, t^*$')

xlabel('$ \cdot 10^{-4}, B$')
% lgd = legend('t_0','t_1','t_2','t_0 \cdot 0.76','t_0 \cdot 1.55');
lgd = legend('t_0','t_1','t_2');
lgd.NumColumns = 3;
set(gca, 'xtick', [6:1:10].*10^(-4),'ytick', [600:300:2000],'TickDir','out'); 
xticklabels([6:1:10])
yticklabels([6:3:20])
xlim( [ 6, 10 ].*10^(-4) );


%%%очь игрек (можно справа)
ax = gca;
ax.YAxisLocation = 'left';



a = findobj(gcf); 
size = 16;
%20

allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alllines,'LineWidth',1);
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

        
        
        
        
figure;

plot(B,t,'r','color',[1 0.5 0])
hold on; plot(B,tcr,'o','color',[1 0 0])
%hold on; plot(B,to,'oc','LineWidth',1)
% hold on; plot(B, 0.76.*t,'--b','LineWidth',1)
hold on; plot(B, 1.35.*t,'--','color',[1 0.5 0])

box on;
grid on;
ylabel('$\cdot 10^{2}, t^*$')

xlabel('$ \cdot 10^{-4}, B$')
% lgd = legend('t_0','t_1','t_2','t_0 \cdot 0.76','t_0 \cdot 1.55');
lgd = legend('t*','t_1','1.35 \cdot t_0');
lgd.NumColumns = 3;
set(gca, 'xtick', [6:1:10].*10^(-4),'ytick', [900:200:2000],'TickDir','out'); 
xticklabels([6:1:10])
yticklabels([9:2:20])
xlim( [ 6, 10 ].*10^(-4) );


%%%очь игрек (можно справа)
ax = gca;
ax.YAxisLocation = 'left';



a = findobj(gcf); 
size = 16;
%20

allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alllines,'LineWidth',1);
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');



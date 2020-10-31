clear all;
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster\dependence_on_B_data';

% шаг из моделирования
fc=1;
dt = 1.*fc;
%пороговое значение производной (в это число раз превышает начальное значение)

M=0.1;
g=M/10;
A=1;
B=[0.0006 0.0007 0.0008 0.0009 0.001];


t=2*M^2*sqrt(2.71828)/g^2./sqrt(B)./A^2;
t1=ceil(t(1));
t2=ceil(t(2));
t3=ceil(t(3));
t4=ceil(t(4));
t5=ceil(t(5));


w=0.03

figure; 
       
subplot( 1, 2, 1 );        

i=0;

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

%время, которое определяет min роивзодной     
        [m,n] = min (der1);
        to(i) = tau(n);
%         
%время, определенное из min*2        
        c=n;
        while der1(c)<= 0.005
            c = c+1;
        end
        tcr(i) = tau(c);
        
        
        if (i==1)
%             hold on;  plot(tau,level,'r')
            hold on; plot(tau(2:t1/fc),der1(2:t1/fc),'b')
            hold on; plot([t1 t1],[0 w],'--b')
        end
        if (i==2)
            hold on; plot(tau(2:t2/fc),der1(2:t2/fc),'c')
            hold on; plot([t2 t2],[0 w],'--c')
        end
        if (i==3)
            hold on; plot(tau(2:t3/fc),der1(2:t3/fc),'m')
            hold on; plot([t3 t3],[0 w],'--m')
        end
        if (i==4)
            hold on; plot(tau(2:t4/fc),der1(2:t4/fc),'black')
            hold on; plot([t4 t4],[0 w],'--black')
        end
        if (i==5)
            hold on; plot(tau(2:t5/fc),der1(2:t5/fc),'green')
            hold on; plot([t5 t5],[0 w],'--green')
        end
end

legend('\Delta_y(t=0)','0.0006','0.0006','0.0007','0.0007','0.0008','0.0008','0.0009','0.0009','0.001','0.001')
grid on;
ylabel( '$\Delta$');
xlabel( '$t$');
box on



B=[0.0006 0.0007 0.0008 0.0009 0.001];



subplot( 1, 2, 2 ); 

plot(B,t,'b')
hold on; plot(B,tcr,'or','LineWidth',1)
% hold on; plot(B,to,'oc','LineWidth',1)
hold on; plot(B, 0.87.*t,'--b','LineWidth',1)
box on;
grid on;
ylabel('$t*$')
xlabel('$B$')
legend('t_0','t_1','t_0 \cdot 0.7')



a = findobj(gcf); 
size = 11;
%20

allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alllines,'LineWidth',1);
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');


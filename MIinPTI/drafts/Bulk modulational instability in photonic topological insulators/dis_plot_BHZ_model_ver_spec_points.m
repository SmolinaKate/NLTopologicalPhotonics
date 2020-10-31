
clear all; close all;


J_2=1/sqrt(2);
J_1=1;
G=1;
g=-G;
M=-1/2;
Delta = M+ 4*J_2;
beta=J_2;

start_E=-10;
end_E=10;

% figure('Renderer', 'painters', 'Position', [10 10 1000 250])
% clf
p = panel();
p.pack(1, 2);


p(1, 1).select();
I_0=2/G.*(-(Delta - 4*J_2) + sqrt(2*J_1^2.*(4*J_2-Delta)./J_2));
I_0=3;
[w,k]=meshgrid([start_E:0.005:end_E],[-3:0.005:3]);
z= (w+g*I_0/2).^2 - 2*J_1^2.*k.^2 - ((M+beta*k.^2)-g*I_0*(M+beta*k.^2)./2./(w+g*I_0)).^2;
v = [0,0];
contour(k,w,z,v,'b');
grid on; box on;
ylabel('$E$')
xlabel('$p_y$')
title(['(a)~~$I_0 = $',num2str(I_0)])

%характерные энергии
ad_loop = - sqrt((4*J_2-Delta)/J_2)-1.5;


E_1=Delta-4*J_2+G*I_0;
E_2=-Delta+4*J_2+G*I_0;
E_3=G*I_0/2;
E_1_b=G*I_0/2-sqrt(2*J_1^2*(4*J_2-Delta)/J_2);
E_2_b=G*I_0/2+sqrt(2*J_1^2*(4*J_2-Delta)/J_2);
E_3_b=G*I_0;
text(0,E_1,'$\leftarrow E_1^0$')
text(0,E_2,'$\leftarrow E_2^0$')
text(0,E_3,'$\leftarrow E_3^0$')
text(ad_loop,E_1_b,'$E_{1 BHZ}^0 \rightarrow$')
text(ad_loop,E_2_b,'$E_{2 BHZ}^0 \rightarrow$')
text(ad_loop,E_3_b,'$E_{3 BHZ}^0 \rightarrow$')


% p(1,2).select();
% I_0=2/G.*(-(Delta - 4*J_2) - sqrt(2*J_1^2.*(4*J_2-Delta)./J_2));
% [w,k]=meshgrid([start_E:0.05:end_E],[-3:0.05:3]);
% z= (w+g*I_0/2).^2 - 2*J_1^2.*k.^2 - ((M+beta*k.^2)-g*I_0*(M+beta*k.^2)./2./(w+g*I_0)).^2;
% v = [0,0];
% contour(k,w,z,v,'b');
% grid on; box on;
% ylabel('$E$')
% xlabel('$p_y$')
% title(['(b)~~$I_0 = $',num2str(I_0)])
% 
% 
% %характерные энергии
% 
% E_1=Delta-4*J_2+G*I_0;
% E_2=-Delta+4*J_2+G*I_0;
% E_3=G*I_0/2;
% E_1_b=G*I_0/2-sqrt(2*J_1^2*(4*J_2-Delta)/J_2);
% E_2_b=G*I_0/2+sqrt(2*J_1^2*(4*J_2-Delta)/J_2);
% E_3_b=G*I_0;
% text(0,E_1,'$\leftarrow E_1^0$')
% text(0,E_2,'$\leftarrow E_2^0$')
% text(0,E_3,'$\leftarrow E_3^0$')
% text(ad_loop,E_1_b,'$E_{1 BHZ}^0 \rightarrow$')
% text(ad_loop,E_2_b,'$E_{2 BHZ}^0 \rightarrow$')
% text(ad_loop,E_3_b,'$E_{3 BHZ}^0 \rightarrow$')
% 
% 
% 
% 
% 
size=10;
a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

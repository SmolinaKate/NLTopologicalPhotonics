clear all; 
close all;

J_2=1/sqrt(2);
J_1=1;
G=1;
Delta =2;

%determine the plot area of dispersion curves
start_E=-10;
end_E=10;
step_E=0.05;
start_k=-3;
end_k=3;
step_k=0.05;


p = panel();
p.pack(1, 3);


p(1, 1).select();
I_0=1;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
contour(k,E,z,v,'b');
grid on; box on;
ylabel('$E$')
xlabel('$p_y$')
title(['(a)~~$I_0 = $',num2str(I_0)])


E_1=Delta-4*J_2+G*I_0;
E_2=-Delta+4*J_2+G*I_0;
E_3=G*I_0/2;
%BHZ's energies
E_1_b=G*I_0/2-sqrt(2*J_1^2*(4*J_2-Delta)/J_2);
E_2_b=G*I_0/2+sqrt(2*J_1^2*(4*J_2-Delta)/J_2);
E_3_b=G*I_0;
text(0,E_1,'$\leftarrow E_1^0$')
text(0,E_2,'$\leftarrow E_2^0$')
text(0,E_3,'$\leftarrow E_3^0$')
ad_loop = - sqrt((4*J_2-Delta)/J_2)-1.5;
text(ad_loop,E_1_b,'$E_{1 BHZ}^0 \rightarrow$')
text(ad_loop,E_2_b,'$E_{2 BHZ}^0 \rightarrow$')
text(ad_loop,E_3_b,'$E_{3 BHZ}^0 \rightarrow$')


p(1,2).select();
I_0=3;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
contour(k,E,z,v,'b');
grid on; box on;
ylabel('$E$')
xlabel('$p_y$')
title(['(b)~~$I_0 = $',num2str(I_0)])


%BHZ's energies
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




p(1,3).select();
I_0=5;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
contour(k,E,z,v,'b');
grid on; box on;
ylabel('$E$')
xlabel('$p_y$')
title(['(b)~~$I_0 = $',num2str(I_0)])


%BHZ's energies
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




size=10;
a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

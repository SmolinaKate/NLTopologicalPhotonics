clear all; 
close all;
LineSize=1;%%LineWidth for all figs
border=4;%%margins
sizeFont=12;
figure('Position', [10 50 400 550]);%%Fig's size
%determine the plot area of dispersion curves (a-b)
E1=-2;
E2=2;
step_E=0.01;
start_k=-2;
end_k=2;
step_k=0.01;
%panel packs
p = panel();
p.pack(2, 1);
p(1,1).pack(1,3);
p(2,1).pack(1,2);
p(2,1,1,1).pack(3,1);
%free space to make 3pi/2 and Gamma I0 = - m_eff at the same line
p(2,1,1,1,3).repack(0.1);
%
p(2,1,1,2).repack(0.6);
p(1,1,1,1).marginright = border; 
p(1,1,1,2).marginleft = border; 
p(1,1,1,2).marginright = border; 
p(1,1,1,3).marginleft = border; 
p(1,1,1,3).marginright = border; 
p(2,1,1,2).marginleft = border; 
p(2,1,1,2).marginright = border; 
p(2).margintop = border; 
p(2).marginbottom = border; 
p(1).margintop = border; 
p(1).marginbottom = border+10; 
p(2,1,1,1).margintop=border+2;
p(2,1,1,1).marginbottom=border+2;
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% (a,b,c)
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%
J_2=1/sqrt(2);
J_1=1;
G=1;
m_eff=1/2;
Delta = 4*J_2+m_eff;
p(1,1,1,1).select();
I_0=1;
start_E=E1+G*I_0;
end_E=E2+G*I_0;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
hold on; contour(k,E-G*I_0,z,v,'r','LineWidth',LineSize); box on;
ylabel('$E-\Gamma I_0$')
title(['(a)~~$\Gamma I_0 = $',num2str(I_0)])
set(gca, 'xtick', [-2 0 2], 'ytick', [-2 -1 0 1 2],'TickDir','in');

p(1,1,1,2).select();
I_0=2;
start_E=E1+G*I_0;
end_E=E2+G*I_0;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
hold on; contour(k,E-G*I_0,z,v,'r','LineWidth',LineSize);
box on;
title(['(b)~~$\Gamma I_0 = $',num2str(I_0)])
yticklabels([])		
set(gca, 'xtick', [-2 0 2], 'ytick', [-2 -1 0 1 2],'TickDir','in');

p(1,1,1,3).select();
I_0=3;
start_E=E1+G*I_0;
end_E=E2+G*I_0;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
contour(k,E-G*I_0,z,v,'r','LineWidth',LineSize);
hold on; box on;
title(['(c)~~$\Gamma I_0 = $',num2str(I_0)])
yticklabels([])		
set(gca, 'xtick', [-2 0 2], 'ytick', [-2 -1 0 1 2],'TickDir','in');
G=1;
m_eff=-1/2;
Delta = 4*J_2+m_eff;

p(1,1,1,1).select();
I_0=1;
start_E=E1+G*I_0;
end_E=E2+G*I_0;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
hold on; contour(k,E-G*I_0,z,v,'b','LineWidth',LineSize); box on;
set(gca, 'xtick', [-2 0 2], 'ytick', [-2 -1 0 1 2],'TickDir','in');

p(1,1,1,2).select();
I_0=2;
start_E=E1+G*I_0;
end_E=E2+G*I_0;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
hold on; contour(k,E-G*I_0,z,v,'b','LineWidth',LineSize);
box on;
set(gca, 'xtick', [-2 0 2], 'ytick', [-2 -1 0 1 2],'TickDir','in');

p(1,1,1,3).select();
I_0=3;
start_E=E1+G*I_0;
end_E=E2+G*I_0;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
contour(k,E-G*I_0,z,v,'b','LineWidth',LineSize);
hold on; box on;
set(gca, 'xtick', [-2 0 2], 'ytick', [-2 -1 0 1 2],'TickDir','in');
%all text
text(0.8,-2.3,'$p_y$')
text(-3.9,-2.3,'$p_y$')
text(-8.5,-2.3,'$p_y$')
text(-10.8,-2.8,'(d)')
text(-10.8,-4.65,'(e)')
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% (d,e)
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%
clear k
step_G=0.001;
G=-4:step_G:4;
J_2=1/sqrt(2);
J_1=1;
I_0 = 1;
inc(1:length(G))=0;
k(1:length(G))=0;
for q=1:1:length(G)
  if (-(G(q)*I_0+J_1^2/J_2)>0)
     inc(q) = abs(G(q)*I_0+J_1^2/J_2);
  end
end
for q=1:1:length(G)
   if (-(G(q)*I_0+J_1^2/J_2)>0)
     k(q)=(sqrt(abs(G(q).*I_0*J_2+J_1^2)/J_2^2));
   end 
end

p(2,1,1,1,1,1).select();
hold on; plot(G,inc,'r')
hold on; plot([-sqrt(2) -sqrt(2)], [0 max(inc)],'--b')
set(gca, 'xtick', [], 'ytick', [-2 -1 0 1 2],'TickDir','in');
ylabel('$\mathrm{max}_{p_y}(\mathrm{Im}[\lambda({\bf{p}})])$')
box on; 
xlim( [ G(1), G(end)] );
ylim( [ inc(end), inc(1)] );

p(2,1,1,1,2,1).select();
hold on; plot(G,k,'r')
hold on; plot([-sqrt(2) -sqrt(2)], [0 max(k)],'--b')
xlabel('$\Gamma I_0 = - 2 m_{\mathrm{eff}}$')
ylabel('$|{\bf{p}}_c|$')
set(gca, 'xtick', [-4 -2 0 2 4], 'ytick', [-2 -1 0 1 2],'TickDir','in');
box on; 
xlim( [ G(1), G(end)] );
ylim( [ k(end), k(1)] );
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% (f)
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%
J1=1;
J2=J1/sqrt(2);
I_0 = 1;
m_eff=[1/2,-1/2];
Delta = 4*J2+m_eff;
G=3;
phi_ar = [0:0.01:2*pi];
for ind_phi=1:1:length(phi_ar)
for ind_Delta=1:1:length(Delta)
    phi=phi_ar(ind_phi);               
    g=G*I_0;
    c=-2*(Delta(ind_Delta)-4*J2)-G*I_0;
        k= -2:0.005:2;
        x1=-(1/sqrt(2))*sqrt(4*J1^2.*k.^2 + 2*c*J2.*k.^2 + 2*g*J2.*k.^2 + 2*J2^2.*k.^4 -... 
             sqrt(2)*sqrt((exp(-2*1i*phi).*k.^2.*(c^2*(-1 + exp(2*1i*phi))^2*J1^2+... 
             2*c*(-1 + exp(2*1i*phi))^2*g*J1^2 + 2*exp(2*1i*phi)*g^2*J2^2.*k.^2))));
        x2=-(1/sqrt(2))*sqrt(4*J1^2.*k.^2 + 2*c*J2.*k.^2 + 2*g*J2.*k.^2 + 2*J2^2.*k.^4 +... 
             sqrt(2)*sqrt((exp(-2*1i*phi).*k.^2.*(c^2*(-1 + exp(2*1i*phi))^2*J1^2+... 
             2*c*(-1 + exp(2*1i*phi))^2*g*J1^2 + 2*exp(2*1i*phi)*g^2*J2^2.*k.^2))));
        x3=(1/sqrt(2))*sqrt(4*J1^2.*k.^2 + 2*c*J2.*k.^2 + 2*g*J2.*k.^2 + 2*J2^2.*k.^4 -... 
             sqrt(2)*sqrt((exp(-2*1i*phi).*k.^2.*(c^2*(-1 + exp(2*1i*phi))^2*J1^2+... 
             2*c*(-1 + exp(2*1i*phi))^2*g*J1^2 + 2*exp(2*1i*phi)*g^2*J2^2.*k.^2))));
        x4=(1/sqrt(2))*sqrt(4*J1^2.*k.^2 + 2*c*J2.*k.^2 + 2*g*J2.*k.^2 + 2*J2^2.*k.^4 +... 
            sqrt(2)*sqrt((exp(-2*1i*phi).*k.^2.*(c^2*(-1 + exp(2*1i*phi))^2*J1^2+... 
            2*c*(-1 + exp(2*1i*phi))^2*g*J1^2 + 2*exp(2*1i*phi)*g^2*J2^2.*k.^2))));
       incr(1) = (max(imag(x1)));  
       incr(2) = (max(imag(x2)));
       incr(3) = (max(imag(x3)));
       incr(4) = (max(imag(x4)));
       incr_cross_band(ind_phi,ind_Delta) = (double(max(incr)));
end
end

p(2,1,1,2).select(); polarplot(phi_ar,(squeeze(incr_cross_band(:,1)))','r')
title(['(f)~$\mathrm{max}_{p_y}(\mathrm{Im}[\lambda({\bf{p}}]),~\Gamma I_0 =$',num2str(I_0*G)])
hold on;  polarplot(phi_ar,(squeeze(incr_cross_band(:,2)))','b')
box on;
pax = gca;
rlim([0 max(max(incr_cross_band))])
rticks(pax,[0 0.5 1 1.5])
rticklabels(pax,{'','0.5','1','1.5'})
thetaticks(pax,[0 45 90 135 180 225 270 315 360])
thetaticklabels(pax,{'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})

a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(alllines,'LineWidth',LineSize)

str=['Fig_2','.png'];
print('-dpng','-r500',str)%%resol 500

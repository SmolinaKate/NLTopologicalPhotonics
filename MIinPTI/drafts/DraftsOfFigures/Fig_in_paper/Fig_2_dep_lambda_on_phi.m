clear all;
close all;
J1=1;
J2=J1/sqrt(2);
I_0 = 1;
m_eff=[1/2,-1/2];
Delta = 4*J2+m_eff;
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
G=2;
phi_ar = [0:0.005:2*pi];
for ind_phi=1:1:length(phi_ar)
for ind_Delta=1:1:length(Delta)
    phi=phi_ar(ind_phi);               
    g=G*I_0;
    c=-2*(Delta(ind_Delta)-4*J2)-G*I_0;
        k= -2:0.001:2;
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
%figure('Position', [10 50 10*10^2 3*10^2]);
ed(1:length(phi_ar))=1;
figure;
p = panel();
p.pack(1, 2);
p(1, 1).select(); polarplot(phi_ar,(squeeze(incr_cross_band(:,1)))','r')
title(['$\Im_\lambda(\varphi),~\Gamma I_0 =$',num2str(I_0*G)])
hold on;  polarplot(phi_ar,(squeeze(incr_cross_band(:,2)))','b')
hold on;  polarplot(phi_ar,ed.*min(incr_cross_band(:,1)),'--r')
box on;
pax = gca;
rlim([0 max(max(incr_cross_band))])
rticks(pax,[0 0.5 1])
rticklabels(pax,{'','0.5','1'})
thetaticks(pax,[0 45 90 135 180 225 270 315 360])
thetaticklabels(pax,{'0', '45', '90', '135', '180', '225', '270', '315', '360'})

size=12;
a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(alllines,'LineWidth',1)
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%                       
G=3;
for ind_phi=1:1:length(phi_ar)
for ind_Delta=1:1:length(Delta)
    phi=phi_ar(ind_phi);               
    g=G*I_0;
    c=-2*(Delta(ind_Delta)-4*J2)-G*I_0;
        k= -2:0.001:2;
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

p(1, 2).select(); polarplot(phi_ar,squeeze(incr_cross_band(:,1))','r')
title(['$\Im_\lambda(\varphi),~\Gamma I_0 =$',num2str(I_0*G)])
hold on;  polarplot(phi_ar,squeeze(incr_cross_band(:,2))','b')
hold on;  polarplot(phi_ar,ed.*min(incr_cross_band(:,1)),'--r')
box on;
rlim([0 max(max(incr_cross_band))])
pax = gca;
rlim([0 max(max(incr_cross_band))])
rticks(pax,[0 0.5 1 1.5])
rticklabels(pax,{'','0.5','1','1.5'})
thetaticks(pax,[0 45 90 135 180 225 270 315 360])
thetaticklabels(pax,{'0', '45', '90', '135', '180', '225', '270', '315', '360'})
%thetaticklabels(pax,{'0','\pi/4','\pi/2','3 \pi/4','\pi','5 \pi/4','3 \pi/2','7 \pi/4','2 \pi'})

a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(alllines,'LineWidth',1)
clear all;
close all;

sizeFont=12;
LineSize=1;
figure('Position', [10 50 700 300]);
p = panel();
p.pack(1, 2);
border=4;
p(1,2).marginleft = border; 
p(1,1).marginright = border; 
p(1,1).repack(0.545);


J2=1/sqrt(2);
J1=1;
I_0 = 1;
step=0.05;
G=-4:step:4;
Delta = -5:step:5;
phi_ar = [0,pi/2];

incr_cross_band=zeros(length(G),length(Delta),length(phi_ar));


for ind_phi=1:1:length(phi_ar)
    phi=phi_ar(ind_phi);
for j = 1:1:length(G)
for d=1:1:length(Delta)                
    g=G(j)*I_0;
    c=-2*(Delta(d)-4*J2)-G(j)*I_0;
    if (I_0/2>=(Delta(d) - 4*J2)/G(j)) && (I_0/2>=-(Delta(d) - 4*J2)/G(j))
        k= -2:0.03:2;
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
       incr_cross_band(j,d,ind_phi) = (double(max(incr)));
     else
       incr_cross_band(j,d,ind_phi) = -1;
     end
end
end
end


delta_c=J1^2/2/J2+4*J2;
bound_stab=2/I_0.*(Delta-4*J2-J1^2/J2);
ed(1:length(Delta))=1;


p(1, 1).select();
title('(a)~~$\varphi=0$')
imagesc(Delta,G,squeeze(incr_cross_band(:,:,1))); box on;
hold on; plot(Delta,-2.*Delta./I_0+8*J2/I_0,'--r')
hold on; plot(Delta,2.*Delta./I_0-8*J2/I_0,'--r')
hold on; plot(Delta,bound_stab,'black')
hold on; plot(Delta,-bound_stab,'black')
hold on; plot([2*sqrt(2) 2*sqrt(2)],[-4 4],'--green')
hold on; plot([delta_c delta_c],[-5 5],'green')
hold on; plot(Delta,ed.*(-2*J1^2/J2/I_0),'b');
hold on; plot(Delta,ed.*(2*J1^2/J2/I_0),'b');
set(gca, 'xtick',[ -5 -2.5 0 2.5  5], 'ytick',[ -4 -2 0 2 4],'TickDir','in','Layer','top');
xlim( [Delta(1) Delta(end)] );
ylim( [G(1) G(end)] );
xlabel('Detuning $\Delta/ J_1$')
ylabel('Nonlinearity $\Gamma I_0 / J_1$')
colorbar('westoutside')


p(1, 2).select();
title('(b)~~$\varphi=\pi/2$')
imagesc(Delta,G,squeeze(incr_cross_band(:,:,2))); box on;
hold on; plot(Delta,-2.*Delta./I_0+8*J2/I_0,'--r')
hold on; plot(Delta,2.*Delta./I_0-8*J2/I_0,'--r')
hold on; plot(Delta,bound_stab,'black')
hold on; plot(Delta,-bound_stab,'black')
hold on; plot([2*sqrt(2) 2*sqrt(2)],[-4 4],'--green')
hold on; plot([delta_c delta_c],[-5 5],'green')
hold on; plot(Delta,ed.*(-2*J1^2/J2/I_0),'b');
hold on; plot(Delta,ed.*(2*J1^2/J2/I_0),'b');
hold on; plot(Delta,-2.*Delta./I_0+8*J2/I_0,'--r')
set(gca, 'xtick',[ -5 -2.5 0 2.5  5], 'ytick',[ -4 -2 0 2 4],'TickDir','in','Layer','top');
xlim( [Delta(1) Delta(end)] );
ylim( [G(1) G(end)] );
yticklabels([])
xlim( [Delta(1) Delta(end)] );
ylim( [G(1) G(end)] );
xlabel('Detuning $\Delta/ J_1$')
colorbar('eastoutside','YTickLabel',[])

text(4.2,4.35,'$\mathrm{Im}[\lambda({\bf{p}}_c)] \geq 2$')
text(6.8,-3.7,'$0$')

a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(alllines,'LineWidth',LineSize)



clear all;
close all;


J2=1/sqrt(2);
J1=1;
I_0 = 1;
step=0.05;
G=-4:step:4;
Delta = -4:step:4;
phi_ar = 0:0.01:pi/2;

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


h = vol3d('cdata',incr_cross_band,'texture','2D','XData',G,'YData',Delta,'ZData',phi_ar);
view(3);  
%axis tight;  daspect([1 1 .4])
alphamap('rampup');
alphamap(.5.* alphamap);
colormap('hot')
grid on;
box on;



xlabel('$\Delta$')
ylabel('$\Gamma$')
zlabel('$\varphi$')


size=10;
a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(alllines,'LineWidth',1)

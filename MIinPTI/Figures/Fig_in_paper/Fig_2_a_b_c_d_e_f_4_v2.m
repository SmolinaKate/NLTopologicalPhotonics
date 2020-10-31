clear all; 
close all;
LineSize=1;%%LineWidth for all figs
border=4;%%margins
sizeFont=12;
figure('Position', [10 50 400 450]);%%Fig's size
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
p(2,1,1,1,3).repack(0.01);
%
p(2,1,1,2).repack(0.56);
p(1,1,1,1).marginright = border; 
p(1,1,1,2).marginleft = border; 
p(1,1,1,2).marginright = border; 
p(1,1,1,3).marginleft = border; 
p(1,1,1,3).marginright = border; 
p(2,1,1,2).marginleft = border; 
p(2,1,1,2).marginright = border; 
p(2).margintop = border*1+0; 
p(2).marginbottom = border; 
p(1).margintop = border; 
p(1).marginbottom = border+10.8; 
p(2,1,1,1).margintop=border; 
p(2,1,1,1).marginbottom=1*border+2*0; % a vertical gap between the left bottom panels
p(2,1,1,2).marginbottom=10; 

%% Identify
% p.select('all');
% % display() the panel object at the prompt
% p.identify();
% p.show();
% p.fontname = 'Times New Roman';
% p.fontsize= sizeFont;

% p.de.margin = 2.5;
% p(1,1).marginbottom = 2.5;
% p(2).marginleft = 3;
p.margin = [12 8 6 6]; % r b l t
% return
%% 
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% (a,b,c)
%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%
J_2=1/sqrt(2);
J_1=1;

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

p(1,1,1,2).select();
I_0=2.4; % 2*2^(1/4);
start_E=E1+G*I_0;
end_E=E2+G*I_0;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
hold on; contour(k,E-G*I_0,z,v,'b','LineWidth',LineSize);
box on;

p(1,1,1,3).select();
I_0=3;
start_E=E1+G*I_0;
end_E=E2+G*I_0;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
contour(k,E-G*I_0,z,v,'b','LineWidth',LineSize);
hold on; box on;

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
hold on; contour(k,E-G*I_0,z,v,'--r','LineWidth',LineSize); box on;
ylabel('$E-\Gamma I_0$')
title(['(a)~~$\Gamma I_0 = $',num2str(I_0)])
set(gca, 'xtick', [-2 -1 0 1 2],'XTickLabel', {'-2',' ','0',' ', '2'}, 'ytick', [-2 -1 0 1 2],'TickDir','in','TickLength',[0.015, 0.015]);
xlim([-2 2])
ylim([-2 2])

p(1,1,1,2).select();
I_0=2.4;
start_E=E1+G*I_0;
end_E=E2+G*I_0;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
hold on; contour(k,E-G*I_0,z,v,'--r','LineWidth',LineSize);
box on;
title(['(b)~~$\Gamma I_0 = $',num2str(I_0)])
yticklabels([])		
set(gca, 'xtick', [-2 -1 0 1 2],'XTickLabel', {'-2',' ','0',' ', '2'}, 'ytick', [-2 -1 0 1 2],'TickDir','in','TickLength',[0.015, 0.015]);
xlim([-2 2])
ylim([-2 2])

p(1,1,1,3).select();
I_0=3;
start_E=E1+G*I_0;
end_E=E2+G*I_0;
[E,k]=meshgrid([start_E:step_E:end_E],[start_k:step_k:end_k]);
z= ((E-G*I_0/2).^2 - 2*J_1^2.*k.^2).*(E-G*I_0).^2 - ((Delta - 4*J_2)+J_2*k.^2).^2.*(E-G*I_0/2).^2;
v = [0,0];
contour(k,E-G*I_0,z,v,'--r','LineWidth',LineSize);
hold on; box on;
title(['(c)~~$\Gamma I_0 = $',num2str(I_0)])
yticklabels([])		
set(gca, 'xtick', [-2 -1 0 1 2],'XTickLabel', {'-2',' ','0',' ', '2'}, 'ytick', [-2 -1 0 1 2],'TickDir','in','TickLength',[0.015, 0.015]);
xlim([-2 2])
ylim([-2 2])

leg = legend({'nontrivial','trivial'},'Location','southeast') ;
legend('boxoff')
rect = [0.12, 0.46, .25, .05];
set(leg, 'Position', rect)

%all text
text(0.8,-2.3,'$p_y$')
text(-3.9,-2.3,'$p_y$')
text(-8.5,-2.3,'$p_y$')
% text(-10.8,-2.8,'(d)')
% text(-10.8,-4.65,'(e)')
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
hold on; plot(G,inc,'m')
set(gca, 'xtick', [-4 -2 0 2 4], ...
  'XTickLabel', {' ',' ',' ',' ', ' '},  ...
    'ytick', [-2 -1 0 1 2 3],'TickDir','in','TickLength',[0.015, 0.015]);
% ylabel('$\mathrm{max}_{p_y}(\mathrm{Im}[\lambda({\bf{p}})])$')
ylabel('$\mathrm{Im}[\lambda(${\boldmath${{p}}$}$_{\mathrm{c}})]$')
box on; 
xlim( [ G(1), G(end)] );
% ylim( [ inc(end), inc(1)] );
ylim( [ 0, 3.001] );
% hold on; plot([-sqrt(2) -sqrt(2)], [0 max(inc)],'--k')
hold on; plot([-sqrt(2) -sqrt(2)], [0 3],'--k','Linewidth',0.5)
text(2.42,2.56,'(d)')

p(2,1,1,1,2,1).select();
hold on; plot(G,k,'m')
hold on; 
xlabel('$\Gamma I_0 = - 2 m_{\mathrm{eff}}$')
ylabel('$|${\boldmath${{p}}$}$_{\mathrm{c}}|$')
set(gca, 'xtick', [-4 -2 0 2 4], 'ytick', [-2 -1 0 1 2],'TickDir','in','TickLength',[0.015, 0.015]);
box on; 
xlim( [ G(1), G(end)] );
% ylim( [ k(end), k(1)] );
ylim( [ 0, 2] );
% plot([-sqrt(2) -sqrt(2)], [0 max(k)],'--k')
plot([-sqrt(2) -sqrt(2)], [0 2],'--k','Linewidth',0.5)
% text(2,4,'(d)')
text(2.42,1.68,'(e)')


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

p(2,1,1,2).select(); 
polarplot(phi_ar,(squeeze(incr_cross_band(:,1)))','--r')
% set(gca, 'FontName','Time New Roman')
% title(['(f)~$\mathrm{max}_{p_y}(\mathrm{Im}[\lambda({\bf{p}}]),~\Gamma I_0 =$',num2str(I_0*G)])
title(['~(f)~$~~\mathrm{Im}[\lambda(${\boldmath${{p}}$}$_{\mathrm{c}})],~\Gamma I_0 =$',num2str(I_0*G)])
title2=get(gca,'Title'); % Get handle to the title object
currentTitlePosition=get(title2,'Position');    % Get the title's current position
set(title2,'Position',currentTitlePosition-[0 0.35 0]); % Shift the position 
% title(['(f)~$\mathrm{max}_{p_y}(\mathrm{Im}[\lambda({\bf{p}}]),~\Gamma I_0 =$',num2str(I_0*G)],'position',[0 -1])
% set(get(gca,'title'),'Position',[5.5 0.0 1.00011]) 
%  text(4, 3, ['(f)~$\mathrm{max}_{p_y}(\mathrm{Im}[\lambda({\bf{p}}]),~\Gamma I_0 =$',num2str(I_0*G)])
hold on;  polarplot(phi_ar,(squeeze(incr_cross_band(:,2)))','b')
box on;
pax = gca;
rlim([0 max(max(incr_cross_band))])
rticks(pax,[0 0.5 1 1.5])
rticklabels(pax,{'','0.5','1','1.5'})
thetaticks(pax,[0 45 90 135 180 225 270 315 360])
thetaticklabels(pax,{'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
pax.FontName = 'Times New Roman';
%%
p(1,1,1,1).select();
text(1.4,-1.8,'I')
p(1,1,1,2).select();
text(1.4,-1.8,'II')

p(1,1,1,1).select();
plot(0,-0.5,'redo','MarkerSize',2.4,'MarkerFaceColor', 'r');
% text(-0.1,-1.7,'$p_I$','color','red')

p(1,1,1,2).select();
plot(0.01,-1.2,'redo','MarkerSize',2.4,'MarkerFaceColor', 'r');
% text(-0.8,-1,'$p_0$')
pB=sqrt(4-Delta/J2);
plot(pB(2),0,'blueo','MarkerSize',2.4,'MarkerFaceColor', 'b');
plot(-pB(2),0,'blueo','MarkerSize',2.4,'MarkerFaceColor', 'b');

p(1,1,1,3).select();
pB=sqrt(4-Delta/J2);
plot(pB(2),0,'blueo','MarkerSize',2.4,'MarkerFaceColor', 'b');
plot(-pB(2),0,'blueo','MarkerSize',2.4,'MarkerFaceColor', 'b');
plot(0.01,-1.5,'redo','MarkerSize',2.4,'MarkerFaceColor', 'r');
%%

a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(alllines,'LineWidth',LineSize)
%% Printing
str=['Fig_2_revised','.png'];
print('-dpng','-r500',str)%%resol 500
set(gcf, 'PaperPositionMode', 'auto')
% print -depsc2 Switching.eps
% print( '-depsc2', 'Fig3', '-painters' );
print( '-depsc2', 'figure2_v5', '-painters' ); 


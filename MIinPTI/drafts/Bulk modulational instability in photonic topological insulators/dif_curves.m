clear all;
close all;


step_p=0.01;
G=-4:step_p:4;
J_2=1/sqrt(2);
J_1=1;
I_0 = 1;
Delta = -4:step_p:4;

background(1:length(Delta),1:length(G))=0;

delta_zv=J_1^2/2/J_2+4*J_2;

bound_stab=2/I_0.*(Delta-4*J_2-J_1^2/J_2);

bound_stab_a=2/I_0.*(-Delta+4*J_2-J_1^2/J_2);

global_bif_p=2/I_0.*(-(Delta - 4*J_2) + sqrt(2*J_1^2.*(4*J_2-Delta)./J_2));
global_bif_m=2/I_0.*((Delta - 4*J_2) - sqrt(2*J_1^2.*(4*J_2-Delta)./J_2));

extra_bif_m=2/I_0.*((Delta - 4*J_2) + sqrt(2*J_1^2.*(4*J_2-Delta)./J_2));
extra_bif_p=2/I_0.*(-(Delta - 4*J_2) - sqrt(2*J_1^2.*(4*J_2-Delta)./J_2));

cross_ad_p = -2*sqrt(2*J_1^2.*(4*J_2-Delta)./J_2/I_0^2);
cross_ad_m = 2*sqrt(2*J_1^2.*(4*J_2-Delta)./J_2/I_0^2);

Gamma_from_inc_ad_loop = sqrt(8*J_1^2.*(-Delta - J_1/2/J_2 + 4*J_2)./J_2./I_0^2);

% test
% [y,x]=meshgrid(G,Delta); %% x is G y is Delta
% v = [0,0];
% z=x.^2/4+y.^2-1;

figure;
imagesc(Delta,G,background); colormap hot; box on;
hold on; plot(Delta,Gamma_from_inc_ad_loop,'--r')
hold on; plot(Delta,-Gamma_from_inc_ad_loop,'--r')
hold on; plot(Delta,-2.*Delta./I_0+8*J_2/I_0,'r')
hold on; plot(Delta,2.*Delta./I_0-8*J_2/I_0,'r')
hold on; plot([2*sqrt(2) 2*sqrt(2)],[-4 4],'--m')
hold on; plot([delta_zv delta_zv],[-4 4],'m')
hold on; plot(Delta,bound_stab,'black')
hold on; plot(Delta,-bound_stab,'black')
hold on; plot(Delta,global_bif_p,'--c')
hold on; plot(Delta,global_bif_m,'--c')
hold on; plot(Delta,extra_bif_m,'--green')
hold on; plot(Delta,extra_bif_p,'--green')
hold on; plot(Delta,cross_ad_p,'*green')
hold on; plot(Delta,cross_ad_m,'*green')
xlabel('$\Delta$')
ylabel('$\Gamma$')


size=10;
a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(alllines,'LineWidth',2)


figure;
imagesc(Delta,G,background); colormap hot; box on;
hold on; plot(Delta,-2.*Delta./I_0+8*J_2/I_0,'r')
hold on; plot(Delta,2.*Delta./I_0-8*J_2/I_0,'r')
hold on; plot(Delta,2.*Delta./I_0,'c')
hold on; plot(Delta,-2.*Delta./I_0,'c')
hold on; plot(Delta,bound_stab,'black')
hold on; plot(Delta,-bound_stab,'black')
hold on; plot(Delta,bound_stab_a,'green')
hold on; plot(Delta,-bound_stab_a,'green')

xlabel('$\Delta$')
ylabel('$\Gamma$')


size=10;
a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(alllines,'LineWidth',2)


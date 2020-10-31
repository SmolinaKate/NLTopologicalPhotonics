clear all;
close all;

% 
% for j = 1:1:length(G)
% for d=1:1:length(Delta)
%     map_to_det_area_of_t(d,j)=4*(Delta(d)-4*J_2)*J_2/(2*J_1^2);
%     if(map_to_det_area_of_t(d,j)<=1) &&(map_to_det_area_of_t(d,j)>0)
%         sgn(d,j)=1;
%     else
%         sgn(d,j)=-1;
%     end
% end
% end

%%%%TEST IMAGESC
% close all
% clear map
% G=-1:1:1;
% J_2=1/sqrt(2);
% J_1=1;
% I_0 = 1;
% Delta = -4:step_p:4;
% 
% for j = 1:1:length(G)
% for d=1:1:length(Delta)
%     map(j,d)=1;
% end
% end
% 
% figure; imagesc(Delta,G,map); box on; 
% hold on;plot(Delta, 1./4.*Delta)
% xlabel('$\Delta$')
% ylabel('$G$')



step_p=0.5;
G=-4:step_p:4;
J_2=1/sqrt(2);
J_1=1;
I_0 = 1;
Delta = -4:step_p:4;


for j = 1:1:length(G)
for d=1:1:length(Delta)
    
    g=-G(j);
    M = Delta(d) - 4*J_2;
    beta=J_2;
 
    %верхн€€ ветка
    psi_2=0;
    psi_1=sqrt(I_0);
    omega = M-g*I_0;
    k_x=0;
    syms Omega x
    kappa_x=[-10:0.5:10];
    parfor ix = 1:1:length(kappa_x) %дл€ каждого каппа составл€ем матрицу
                Mat=[-Omega-omega+(M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*(psi_1)^2,- g*psi_1^2,k_x- sqrt(2).*J_1.*kappa_x(ix),0;
                k_x- sqrt(2).*J_1.*kappa_x(ix),0,-Omega - omega  - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2, - g*psi_2^2; 
                - g*psi_1^2,Omega - omega + (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_1^2,0,k_x + sqrt(2).*J_1.*kappa_x(ix); 
                0,k_x+ sqrt(2).*J_1.*kappa_x(ix),- g*psi_2^2,Omega - omega - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2]; 
                eqn = det(Mat) == 0; 
                freq(ix,:) = solve(eqn,Omega); 
                incr(ix) = (max(imag(freq(ix,:))));                   
    end
    incr_up_band(j,d) = (double(max(incr)));
%     
%     
%     %%нижн€€ ветка
%     omega = -M-g*I_0;
%     psi_1=0;
%     psi_2=sqrt(I_0);
%     k_x=0;
%     syms Omega x
%     kappa_x= [-10:0.5:10];
%     parfor ix = 1:1:length(kappa_x) %дл€ каждого каппа составл€ем матрицу
%                 Mat=[-Omega-omega+(M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*(psi_1)^2,- g*psi_1^2,k_x- sqrt(2).*J_1.*kappa_x(ix),0;
%                 k_x- sqrt(2).*J_1.*kappa_x(ix),0,-Omega - omega  - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2, - g*psi_2^2; 
%                 - g*psi_1^2,Omega - omega + (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_1^2,0,k_x + sqrt(2).*J_1.*kappa_x(ix); 
%                 0,k_x+ sqrt(2).*J_1.*kappa_x(ix),- g*psi_2^2,Omega - omega - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2]; 
%                 eqn = det(Mat) == 0; 
%                 freq(ix,:) = solve(eqn,Omega); 
%                 incr(ix) = (max(imag(freq(ix,:))));                   
%     end
%     incr_l_band(j,d) = (double(max(incr)));
%     
%     
    %%крест
    omega =-g*I_0/2;
    k_x=0;
    if (I_0/2>=-M/g) && (I_0/2>=M/g)
        psi_1=sqrt(I_0/2+M/g);
        psi_2=sqrt(I_0/2-M/g);
        syms Omega x
    kappa_x= [-10:0.5:10];
    parfor ix = 1:1:length(kappa_x) %дл€ каждого каппа составл€ем матрицу
                Mat=[-Omega-omega+(M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*(psi_1)^2,- g*psi_1^2,k_x- sqrt(2).*J_1.*kappa_x(ix),0;
                k_x- sqrt(2).*J_1.*kappa_x(ix),0,-Omega - omega  - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2, - g*psi_2^2; 
                - g*psi_1^2,Omega - omega + (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_1^2,0,k_x + sqrt(2).*J_1.*kappa_x(ix); 
                0,k_x+ sqrt(2).*J_1.*kappa_x(ix),- g*psi_2^2,Omega - omega - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2]; 
                eqn = det(Mat) == 0; 
                freq(ix,:) = solve(eqn,Omega); 
                incr(ix) = (max(imag(freq(ix,:))));                   
    end
        incr_cross_band(j,d) = (double(max(incr)));
    else
        incr_cross_band(j,d) = -1;
    end

end
end

delta_zv=J_1^2/2/J_2+4*J_2;

bound_stab=2/I_0.*(Delta-4*J_2-J_1^2/J_2);

global_bif_p=2/I_0.*(-(Delta - 4*J_2) + sqrt(2*J_1^2.*(4*J_2-Delta)./J_2));
global_bif_m=2/I_0.*((Delta - 4*J_2) - sqrt(2*J_1^2.*(4*J_2-Delta)./J_2));

extra_bif_m=2/I_0.*((Delta - 4*J_2) + sqrt(2*J_1^2.*(4*J_2-Delta)./J_2));
extra_bif_p=2/I_0.*(-(Delta - 4*J_2) - sqrt(2*J_1^2.*(4*J_2-Delta)./J_2));

cross_ad_p = -2*sqrt(2*J_1.*(4*J_2-Delta)./J_2/I_0^2);
cross_ad_m = 2*sqrt(2*J_1.*(4*J_2-Delta)./J_2/I_0^2);



% test
% [y,x]=meshgrid(G,Delta); %% x is G y is Delta
% v = [0,0];
% z=x.^2/4+y.^2-1;
[y,x]=meshgrid(G,Delta);
v = [0,0];
z_1= I_0/2+sqrt(I_0^2./4-2*J_1.^2.*sqrt((4*J_2-y)./J_2).^2./x.^2);
z_2= I_0/2-sqrt(I_0^2./4-2*J_1.^2.*sqrt((4*J_2-y)./J_2).^2./x.^2);

p = panel();
p.pack(1, 2);


p(1, 1).select();
title('(a)~~The instability of the upper band')
imagesc(Delta,G,abs(incr_up_band)); box on;
hold on; plot(Delta,-2.*Delta./I_0+8*J_2/I_0,'r')
hold on; plot(Delta,2.*Delta./I_0-8*J_2/I_0,'r')
hold on; plot([2*sqrt(2) 2*sqrt(2)],[-4 4],'--m')
hold on; plot([delta_zv delta_zv],[-4 4],'m')
hold on; plot(Delta,bound_stab,'+black')
hold on; plot(Delta,global_bif_p,'--c')
hold on; plot(Delta,global_bif_m,'--c')
hold on; plot(Delta,extra_bif_m,'--green')
hold on; plot(Delta,extra_bif_p,'--green')
hold on; plot(Delta,cross_ad_p,'*green')
hold on; plot(Delta,cross_ad_m,'*green')
hold on; contour(Delta,G,real(z_1),v,'*black');
hold on; contour(Delta,G,real(z_2),v,'*black');
xlabel('$\Delta$')
ylabel('$\Gamma$')
% 
% 
% 
% 
% p(2, 1).select();
% title('(b)~~The instability of the lower band')
% imagesc(Delta,G,abs(incr_l_band)); box on;
% hold on; plot(Delta,-2.*Delta./I_0+8*J_2/I_0,'r')
% hold on; plot(Delta,2.*Delta./I_0-8*J_2/I_0,'r')
% hold on; plot([2*sqrt(2) 2*sqrt(2)],[-4 4],'--m')
% hold on; plot(Delta,bif_of_loops_p,'--black')
% hold on; plot(Delta,bif_of_loops_m,'--black')
% 
% xlabel('$\Delta$')
% ylabel('$\Gamma$')






p(1, 2).select();
title('(c)~~The instability of the cross point')
imagesc(Delta,G,incr_cross_band); box on;
hold on; plot(Delta,-2.*Delta./I_0+8*J_2/I_0,'r')
hold on; plot(Delta,2.*Delta./I_0-8*J_2/I_0,'r')
hold on; plot([2*sqrt(2) 2*sqrt(2)],[-4 4],'--m')
hold on; plot(Delta,bound_stab,'+black')
hold on; plot([delta_zv delta_zv],[-4 4],'m')
hold on; plot(Delta,global_bif_p,'--c')
hold on; plot(Delta,global_bif_m,'--c')
hold on; plot(Delta,extra_bif_m,'--green')
hold on; plot(Delta,extra_bif_p,'--green')
hold on; plot(Delta,cross_ad_p,'*green')
hold on; plot(Delta,cross_ad_m,'*green')
hold on; contour(Delta,G,real(z_1),v,'*black');
hold on; contour(Delta,G,real(z_2),v,'*black');
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

clear all;
close all;
 
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\Bulk modulational instability in photonic topological insulators\data';


J_2=1/sqrt(2);
J_1=1;
I_0 = 1;
step_p=0.5;
G=-4:step_p:4;
Delta = -4:step_p:4;


for j = 1:1:length(G)
for d=1:1:length(Delta)
    
    g=-G(j);
    M = Delta(d) - 4*J_2;
    beta=J_2;
 
    omega =-g*I_0/2;
    if (I_0/2>=-M/g) && (I_0/2>=M/g)
        psi_1=sqrt(I_0/2+M/g);
        psi_2=sqrt(I_0/2-M/g);
        syms Omega x
        kappa_x= [-1:0.1:1];
    parfor ix = 1:1:length(kappa_x) %для каждого каппа составляем матрицу
                Mat=[-Omega-omega+(M+beta*kappa_x(ix)^2) - 2*g*(psi_1)^2,-g*psi_1^2,- sqrt(2).*J_1.*kappa_x(ix),0;
                - sqrt(2).*J_1.*kappa_x(ix),0,-Omega - omega  - (M+beta*kappa_x(ix)^2) - 2*g*psi_2^2, - g*psi_2^2; 
                -g*psi_1^2,Omega - omega + (M+beta*kappa_x(ix)^2) - 2*g*psi_1^2,0,sqrt(2).*J_1.*kappa_x(ix); 
                0,sqrt(2).*J_1.*kappa_x(ix),-g*psi_2^2,Omega - omega - (M+beta*kappa_x(ix)^2) - 2*g*psi_2^2]; 
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


figure; 
p = panel();
p.pack(1, 1);
p(1, 1).select();
title('The instability of the cross point')
imagesc(Delta,G,incr_cross_band); box on;
hold on; plot(Delta,-2.*Delta./I_0+8*J_2/I_0,'r')
hold on; plot(Delta,2.*Delta./I_0-8*J_2/I_0,'r')
hold on; plot([2*sqrt(2) 2*sqrt(2)],[-4 4],'--m')
hold on; plot(Delta,bound_stab,'black')
hold on; plot(Delta,-bound_stab,'black')
hold on; plot([delta_zv delta_zv],[-4 4],'m')
xlabel('$\Delta$')
ylabel('$\Gamma$')


ed(1:length(Delta))=1;
hold on; plot(Delta,ed.*(-2*J_1^2/J_2/I_0),'b');
hold on; plot(Delta,ed.*(2*J_1^2/J_2/I_0),'b');


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
        
        
%file_name = sprintf('map_cross');
%ToSaveFileName = [file_name '.mat'];
%ToSaveFileName1 = fullfile(path, ToSaveFileName);
%save(ToSaveFileName1,'G','Delta','incr_cross_band');                    


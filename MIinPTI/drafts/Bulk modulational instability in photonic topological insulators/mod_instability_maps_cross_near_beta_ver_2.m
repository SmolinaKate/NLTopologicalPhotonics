clear all
close all;

path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\Bulk modulational instability in photonic topological insulators\data';

step_p=0.25;
J_2=1/sqrt(2);
J_1=1;
I_0 = 1;
G=-4:step_p:4;
Delta = -4:step_p:4;


for j = 1:1:length(G)
for d=1:1:length(Delta)
    p_b= - sqrt((4*J_2-Delta(d))./J_2);
    root=I_0^2/4-2*J_1^2.*p_b.^2./G(j)^2;
    mod_psi_1_in2=I_0/2+sqrt(root);
    mod_psi_2_in2=I_0/2-sqrt(root);

    if (4*J_2-Delta(d)>=0) && (mod_psi_1_in2>=0)&& (mod_psi_2_in2>=0)&& root>=0
        
            g=-G(j);
            M = Delta(d) - 4*J_2;
            beta=J_2;
            %CROSS
            psi_1=sqrt(mod_psi_1_in2);
            psi_2=sqrt(mod_psi_2_in2);
            omega = -g*I_0;
            k_x=p_b;
            syms Omega 
            kappa_x= [-10:0.5:10];
            parfor ix = 1:1:length(kappa_x) %для каждого каппа составляем матрицу
                        Mat=[-Omega-omega+(M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*(psi_1)^2,- g*psi_1^2,- sqrt(2).*J_1.*k_x- sqrt(2).*J_1.*kappa_x(ix),0;
                        - sqrt(2).*J_1.*k_x- sqrt(2).*J_1.*kappa_x(ix),0,-Omega - omega  - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2, - g*psi_2^2; 
                        - g*psi_1^2,Omega - omega + (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_1^2,0,-sqrt(2).*J_1.*k_x + sqrt(2).*J_1.*kappa_x(ix); 
                        0,-sqrt(2).*J_1.*k_x+ sqrt(2).*J_1.*kappa_x(ix),- g*psi_2^2,Omega - omega - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2]; 
                        eqn = det(Mat) == 0; 
                        freq(ix,:) = solve(eqn,Omega); 
                        incr(ix) = (max(imag(freq(ix,:))));                   
            end
            incr_cross_band(j,d) = (double(max(incr)));
    else
            incr_cross_band(j,d) = 0;
    end

end
end



figure; 
p = panel();
p.pack(1, 1);
p(1, 1).select();
imagesc(Delta,G,abs(incr_cross_band)); box on;
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
        
file_name = sprintf('map_ad_cross');
ToSaveFileName = [file_name '.mat'];
ToSaveFileName1 = fullfile(path, ToSaveFileName);
save(ToSaveFileName1,'G','Delta','incr_cross_band');                    

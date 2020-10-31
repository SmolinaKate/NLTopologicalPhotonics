clear all;
close all;
 
path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\Bulk modulational instability in photonic topological insulators\data';


J2=1/sqrt(2);
J1=1;
I_0 = 1;
step_p=0.05;
G=-4:step_p:4;
Delta = -4:step_p:4;


phi_ar = 0:0.05:pi/2;



for ind_phi=1:1:length(phi_ar)
    phi=phi_ar(ind_phi);
        for j = 1:1:length(G)
        for d=1:1:length(Delta)                
                  g=G(j)*I_0;
                  c=-2*(Delta(d)-4*J2)-G(j)*I_0;
                   if (I_0/2>=(Delta(d) - 4*J2)/G(j)) && (I_0/2>=-(Delta(d) - 4*J2)/G(j))
                                k= [-2:0.1:2];
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

%         delta_zv=J1^2/2/J2+4*J2;
% 
%         bound_stab=2/I_0.*(Delta-4*J2-J1^2/J2);
% 
%         p(1,1).select();
%         title(['The instability of the cross point for $\varphi$ = ',num2str(phi)])
%         imagesc(Delta,G,incr_cross_band); box on;
%         hold on; plot(Delta,-2.*Delta./I_0+8*J2/I_0,'r')
%         hold on; plot(Delta,2.*Delta./I_0-8*J2/I_0,'r')
%         hold on; plot([2*sqrt(2) 2*sqrt(2)],[-4 4],'--m')
%         hold on; plot(Delta,bound_stab,'black')
%         hold on; plot(Delta,-bound_stab,'black')
%         hold on; plot([delta_zv delta_zv],[-4 4],'m')
%         xlabel('$\Delta$')
%         ylabel('$\Gamma$')
% 
% 
%         ed(1:length(Delta))=1;
%         hold on; plot(Delta,ed.*(-2*J1^2/J2/I_0),'b');
%         hold on; plot(Delta,ed.*(2*J1^2/J2/I_0),'b');
% 
       


% size=10;
% a = findobj(gcf); 
% allaxes  = findall(a,'Type','axes');
% alllines = findall(a,'Type','line');
% alltext  = findall(a,'Type','text');
% set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal');
% set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
% set(alllines,'LineWidth',1)
        
file_name = sprintf('3D_map_cross');
ToSaveFileName = [file_name '.mat'];
ToSaveFileName1 = fullfile(path, ToSaveFileName);
save(ToSaveFileName1,'G','Delta','phi_ar','incr_cross_band');   

% dlmwrite('C:\Users\smoli\Dropbox\Topological Nonlinear 2D\Bulk modulational instability in photonic topological insulators\data\3d_map.dat',incr_cross_band)


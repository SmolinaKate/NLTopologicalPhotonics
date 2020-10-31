clear all
close all

%%%%%%%%BULK STATE (NonDis)
k_min=-2;
delta_k=0.06;
delta_w=delta_k;
k_max=2;
I_0=3;
I_1=I_0/2;
M=1;
g=1;
eta=1;

[w,k]=meshgrid([k_min:delta_k:k_max],[k_min:delta_k:k_max]);
z= (w+g*I_0/2).^2 - (k-eta.*k^2).^2 - (M-g*I_0*M./2./(w+g*I_0)).^2;
v = [0,0];
figure(2);
contour(k,w,z,v,'b','LineWidth',2);
hold on;
axis equal
grid on;
ylabel('$\omega$','Interpreter', 'latex','FontSize', 15)
xlabel('$k$','Interpreter', 'latex','FontSize', 15)
title(['I_0 = ',num2str(I_0),';I_1 = ',num2str(I_1),';M = ',num2str(M),';g = ',num2str(g)]);

kappa_x=[-2:0.1:2];

for ind_k=1:1:length(k)
        ind_k
        k_x=k(ind_k);
        omega = w(ind_k);
        psi_2=I_0/2-M*I_0/(2*(g*I_0+ omega));
        psi_1=I_0/2+M*I_0/(2*(g*I_0+ omega));
        syms Omega
        %for iy = 1:1:length(kappa_y)
            %iy
        for ix = 1:1:length(kappa_x)
            Mat=[-Omega-omega+M - 2*g*(psi_1)^2,- g*psi_1^2,k_x+kappa_x(ix)- eta*k_x^2-eta*kappa_x(ix)^2,0;
            k_x+kappa_x(ix)- eta*k_x^2-eta*kappa_x(ix)^2,0,-Omega - omega  - M - 2*g*psi_2^2, - g*psi_2^2; 
            - g*psi_1^2,Omega - omega + M - 2*g*psi_1^2,k_x-kappa_x(ix)- eta*k_x^2-eta*kappa_x(ix)^2,0; 
            0,k_x-kappa_x(ix)- eta*k_x^2-eta*kappa_x(ix)^2,- g*psi_2^2,Omega - omega - M - 2*g*psi_2^2]; 
            eqn = det(Mat) == 0;
            freq(ix,:) = solve(eqn,Omega);
            frq1(ix,:) = double(freq(ix,:));
            frq(ix,:) = sort(frq1(ix,:));
            frq(ix,:) = imag(frq1(ix,:));
            [max_val(ix),max_n(ix)] = max(frq(ix,:));
        end
        incr(ind_k) = max(max_val);
end
figure(2); hold on;
plot(k,incr,'.r')















%%%%%%%%BULK STATE (MI)
% M=1;
% g=1;
% eta=1;
clear all; close all;
link=1;%link=1
a=abs(2/link/sqrt(3));
eta = link*a^2/8;

M = link/10;
g=sqrt(eta);


I0=1;
omega=-2;
k_x=0;

kappa_x=[0:0.1:10];
psi_2=I0/2-M*I0/(2*(g*I0+ omega));
psi_1=I0/2+M*I0/(2*(g*I0+ omega));

kappa_y=0;
k_y=0;
syms Omega

%for iy = 1:1:length(kappa_y)
    %iy
for ix = 1:1:length(kappa_x)
            Mat=[-Omega-omega+M - 2*g*(psi_1)^2,- g*psi_1^2,k_x+kappa_x(ix)- eta*k_x.^2-eta*kappa_x(ix)^2,0;
            k_x+kappa_x(ix)- eta*k_x.^2-eta*kappa_x(ix)^2,0,-Omega - omega  - M - 2*g*psi_2^2, - g*psi_2^2; 
            - g*psi_1^2,Omega - omega + M - 2*g*psi_1^2,k_x-kappa_x(ix)- eta*k_x.^2-eta*kappa_x(ix)^2,0; 
            0,k_x-kappa_x(ix)- eta*k_x.^2-eta*kappa_x(ix)^2,- g*psi_2^2,Omega - omega - M - 2*g*psi_2^2]; 
            eqn = det(Mat) == 0;
            freq(ix,:) = solve(eqn,Omega);
            frq1(ix,:) = double(freq(ix,:));
            frq(ix,:) = sort(frq1(ix,:));
            frq(ix,:) = imag(frq1(ix,:));
end


figure;
for n=1:1:4
    if (n==1)
            plot(kappa_x,real(frq(:,n)),'.r');
            hold on; plot(kappa_x,imag(frq(:,n)),'b')
            legend('Re(\Omega)','Im(\Omega)') 
    else
            plot(kappa_x,real(frq(:,n)),'.r','HandleVisibility','off');
            hold on; plot(kappa_x,imag(frq(:,n)),'b','HandleVisibility','off')
    end
end

xlabel('$\kappa_x$')
ylabel('$\Omega$')

grid on; box on;
a = findobj(gcf); 
size = 12;

allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alllines,'LineWidth',1);
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
        
        

       
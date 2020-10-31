clear all
close all;
%%%%%%%%BULK STATE (NonDis)

%параметры
% k_min=-1;
% delta_k=0.15;
% delta_w=delta_k;
% k_max=1;
k_min=-0.24;
delta_k=0.04;
delta_w=delta_k;
k_max=0.24;
% I_0=1;
I_0=3;
M=1;
g=1;


%построение нелин дисперсии
%[w,k]=meshgrid([k_min:delta_w:k_max],[k_min:delta_k:k_max]);
%z= (w+g*I_0/2).^2 - k.^2 - (M-g*I_0*M./2./(w+g*I_0)).^2;
%v = [0,0];
%figure(2);
%contour(k,w,z,v,'b','LineWidth',2);
%hold on;
%axis equal
%grid on;
%ylabel('$\omega$','Interpreter', 'latex','FontSize', 15)
%xlabel('$k$','Interpreter', 'latex','FontSize', 15)
%title(['I_0 = ',num2str(I_0),';I_1 = ',num2str(I_1),';M = ',num2str(M),';g = ',num2str(g)]);

%составление матрицы w_solved(размер по к,потенциальное число корней дисперс уравнения (это же многочлен))
k=[k_min:delta_k:k_max];
syms w
for i=1:1:length(k)
    i
    eq=(w+g*I_0/2).^2 - k(i).^2 - (M-g*I_0*M./2./(w+g*I_0)).^2==0;
    w_solved_str=solve(eq,w);
    w_solved(i,1:size(w_solved_str)) = sort(double(w_solved_str));%все корни, включая комплексные -- последние надо исключить
end

cr_val=1000000;%будем присваивать частоте это значение, если у нее есть мнимая часть -- такие корни нам сейчас не надо учитывать
for i=1:1:length(k)
    i
    for j=1:1:size(w_solved_str)
        if imag(w_solved(i,j))~=0||(w_solved(i,j)==0)
            w_solved(i,j)=cr_val;
        end
    end
end

figure(2);
for i=1:1:length(k)
    i
          j=1;
      if w_solved(i,j)~=cr_val %если у частоты нет мнимой части, то построим точку
           hold on; plot(k(i),w_solved(i,j),'.c','LineWidth',2);
      end
          j=2;
      if w_solved(i,j)~=cr_val %если у частоты нет мнимой части, то построим точку
           hold on; plot(k(i),w_solved(i,j),'.r','LineWidth',2);
      end
          j=3;
      if w_solved(i,j)~=cr_val %если у частоты нет мнимой части, то построим точку
           hold on; plot(k(i),w_solved(i,j),'.b','LineWidth',2);
      end
          j=4;
      if w_solved(i,j)~=cr_val %если у частоты нет мнимой части, то построим точку
           hold on; plot(k(i),w_solved(i,j),'.m','LineWidth',2);
      end
    
end
ylabel('$\omega$')
xlabel('$k$')

% kappa_x=[-4:2*delta_k:4];
kappa_x=[-150:20*delta_k:150];


for ind_k=1:1:length(k)
for ind_w=1:1:size(w_solved_str) 
        ind_k
        k_x=k(ind_k);
        omega = w_solved(ind_k,ind_w); %зафиксировали точку по оси к и номер решения уравнения для частоты
        if(w_solved(ind_k,ind_w)~=cr_val) %рассматриваем только действительные \omega
            psi_2=sqrt(I_0/2-M*I_0/(2*(g*I_0+ omega)));
            psi_1=sqrt(I_0/2+M*I_0/(2*(g*I_0+ omega)));
            syms Omega
            parfor ix = 1:1:length(kappa_x) %для каждого каппа составляем матрицу
                Mat=[-Omega-omega+M - 2*g*(psi_1)^2,- g*psi_1^2,k_x+kappa_x(ix),0;
                k_x+kappa_x(ix),0,-Omega - omega  - M - 2*g*psi_2^2, - g*psi_2^2; 
                - g*psi_1^2,Omega - omega + M - 2*g*psi_1^2,0,k_x-kappa_x(ix); 
                0,k_x-kappa_x(ix),- g*psi_2^2,Omega - omega - M - 2*g*psi_2^2]; 
                eqn = det(Mat) == 0; %приравн детерминант к нулю
                freq(ix,:) = solve(eqn,Omega); %получившиеся частоты 
                frq1(ix,:) = double(freq(ix,:));
                frq(ix,:) = sort(frq1(ix,:));
                frq(ix,:) = imag(frq1(ix,:));
                [max_val(ix),max_n(ix)] = max(frq(ix,:)); %икремент для фиксированного каппа (разные \omega и k)
            end
            incr_2(ind_k,ind_w) = max(max_val); %макс инкремент по всем каппа
        end
end
end

a = findobj(gcf); 
size = 12;


figure(2); hold on; plot(k,incr_2(:,1),'+c')
hold on; plot(k,incr_2(:,2),'<r')

allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alllines,'LineWidth',1);
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
grid on; box on; hold on;

hold on; plot(k,incr_2(:,3),'*b')
hold on; plot(k,incr_2(:,4),'om')






% figure;
% for n=1:1:4
%     if (n==1)
%             plot(kappa_x,real(frq(:,n)),'.r');
%             hold on; plot(kappa_x,imag(frq1(:,n)),'.b')
%             legend('Re(\Omega)','Im(\Omega)') 
%     else
%             plot(kappa_x,real(frq(:,n)),'.r','HandleVisibility','off');
%             hold on; plot(kappa_x,imag(frq1(:,n)),'.b','HandleVisibility','off')
%     end
% end












%%%%%%%
%%%%%%%
%%%%%%%

%%%%%%%%BULK STATE (MI)
% M=0.1;
% g=M/3;
% I0=1;
% omega=-g*I0/4;
% k_x=-g*I0/4;
% 
% kappa_x=[-2:0.1:2];
% psi_2=I0/2-M*I0/(2*(g*I0+ omega));
% psi_1=I0/2+M*I0/(2*(g*I0+ omega));
% 
% kappa_y=0;
% k_y=0;
% syms Omega
% 
% %for iy = 1:1:length(kappa_y)
%     %iy
% for ix = 1:1:length(kappa_x)
%     Mat=[-Omega-omega+M - 2*g*(psi_1)^2,- g*psi_1^2,k_x+kappa_x(ix)-1i*kappa_y- 1i*k_y,0;
%     k_x+kappa_x(ix)+1i*kappa_y + 1i*k_y,0,-Omega - omega  - M - 2*g*psi_2^2, - g*psi_2^2; 
%     - g*psi_1^2,Omega - omega + M - 2*g*psi_1^2,k_x-kappa_x(ix)-1i*kappa_y + 1i*k_y,0; 
%     0,k_x-kappa_x(ix)+1i*kappa_y - 1i*k_y,- g*psi_2^2,Omega - omega - M - 2*g*psi_2^2]; 
% 
%     eqn = det(Mat) == 0;
%     freq(ix,:) = solve(eqn,Omega);
%     frq1(ix,:) = double(freq(ix,:));
%    %frq(ix,:) = sort(frq1(ix,:));
%     frq(ix,:) = sort(frq1(ix,:));
% end
% 
% 
% figure;
% for n=1:1:4
%     if (n==1)
%             plot(kappa_x,real(frq(:,n)),'.r');
%             hold on; plot(kappa_x,imag(frq(:,n)),'b')
%             legend('Re(\Omega)','Im(\Omega)') 
%     else
%             plot(kappa_x,real(frq(:,n)),'.r','HandleVisibility','off');
%             hold on; plot(kappa_x,imag(frq(:,n)),'b','HandleVisibility','off')
%     end
% end
% 
% xlabel('$\kappa_x$')
% ylabel('$\Omega$')
% 
% grid on; box on;
% a = findobj(gcf); 
% size = 12;
% 
% allaxes  = findall(a,'Type','axes');
% alllines = findall(a,'Type','line');
% alltext  = findall(a,'Type','text');
% 
% set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal');
% set(alllines,'LineWidth',1);
% set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

%%%%%%%%bulk STATE
% clear all; 
% 
% A_0=1;
% B_0=1;
% 
% 
% link=1;%link=1
% a=abs(2/link/sqrt(3));
% eta = link*a^2/8;
% M=0.1;
% g=M/3;
% 
% 
% omega = -g/2;
% 
% 
% f = (omega*M)/((omega + g*( A_0^2 + B_0^2))^2);
% mu =eta;
% 
% 
% A_1 = B_0*f/2/A_0;
% A_2 = 1 + f/2;
% B_1 = - A_0*f/2/B_0;
% B_2 = 1 - f/2;
% 
% 
% 
% kappa_x=([-1:0.04:1]);
% kappa_y=0;
% syms Omega;
% 
% 
% for ix = 1:1:length(kappa_x)
%     Mat=[-Omega*A_1 + (kappa_x(ix) -  1i*kappa_y - mu*(kappa_x(ix) +  1i*kappa_y)^2)*B_2 + M*A_1 - 2*g*A_0^2*A_1, -Omega*A_2 + (kappa_x(ix) -  1i*kappa_y - mu*(kappa_x(ix) +  1i*kappa_y)^2)*B_1 + M*A_2 - 2*g*A_0^2*A_2, - g*A_0^2*A_1,  - g*A_0^2*A_2;
%     -Omega*B_2 + (kappa_x(ix)+1i*kappa_y - mu*(kappa_x(ix) - 1i*kappa_y)^2)*A_1 - M*B_2 - 2*g*B_0^2*B_2, -Omega*B_1 + (kappa_x(ix) + 1i*kappa_y-mu*(kappa_x(ix) -  1i*kappa_y)^2)*A_2 - M*B_1 + 2*g*B_0^2*B_1, - g*B_0^2*B_2, - g*B_0^2*B_1;   
%      - g*A_0^2*A_1, - g*A_0^2*A_2,Omega*A_1 + ( -kappa_x(ix) -  1i *kappa_y - mu*(-kappa_x(ix)+1i*kappa_y)^2)*B_2 + M*A_1 - 2*g*A_0^2*A_1, Omega*A_2 + ( -kappa_x(ix) -  1i*kappa_y - mu*(-kappa_x(ix) +  1i*kappa_y)^2)*B_1 + M*A_2 - 2*g*A_0^2*A_2; 
%     - g*B_0^2*B_2, - g*B_0^2*B_1,Omega*B_2 + (-kappa_x(ix)+1i*kappa_y-mu*(-kappa_x(ix) -  1i*kappa_y)^2)*A_1 - M*B_2 - 2*g*B_0^2*B_2, Omega*B_1 + (-kappa_x(ix) + 1i*kappa_y - mu*(-kappa_x(ix) -  1i*kappa_y)^2)*A_2 - M*B_1 + 2*g*B_0^2*B_1];
%     
%     eqn = det(Mat) == 0;
%     freq(ix,:) = solve(eqn,Omega);
%     frq1(ix,:) = double(freq(ix,:));
%     frq(ix,:) = sort(frq1(ix,:));
% end
% 
% 
% 
% figure; imagesc([1:1:4],kappa_x,imag(frq))
% xlabel('$\Omega_n$')
% ylabel('$\kappa_x$')
% 
% figure; box on; grid on; 
% for i = 1:1:4
%     hold on; plot(kappa_x,imag(frq(:,i)))
% end
% 
% ylabel('$\Omega_n$')
% xlabel('$\kappa_x$')
% 
% a = findobj(gcf); 
% size = 12;
% 
% allaxes  = findall(a,'Type','axes');
% alllines = findall(a,'Type','line');
% alltext  = findall(a,'Type','text');
% 
% set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal');
% set(alllines,'LineWidth',1);
% set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
%             'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
%         
        
        
        
        

% A=1;
% B1=[0.0006 0.0007 0.0008 0.0009 0.001];


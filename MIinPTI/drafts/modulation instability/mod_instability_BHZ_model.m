clear all
close all;
%%%%%%%%BULK STATE (NonDis)

%параметры
% k_min=-1;
% delta_k=0.15;
% delta_w=delta_k;
% k_max=1;
k_min=-2;
delta_k=0.04;
delta_w=delta_k;
k_max=2;
beta= - 1;
% I_0=1;
I_0=1;
M=1;
g=1;


%составление матрицы w_solved(размер по к,потенциальное число корней дисперс уравнения (это же многочлен))
k=[k_min:delta_k:k_max];
syms w
for i=1:1:length(k)
    i
    eq=(w+g*I_0/2).^2 - k(i).^2 - ((M+beta*k(i).^2)-g*I_0*(M+beta*k(i).^2)./2./(w+g*I_0)).^2==0;
    w_solved_str=solve(eq,w);
    %%%%%
    %%%%%%
    %%%%%%% для петли сделать тут сорт
    w_solved(i,1:size(w_solved_str)) =sort(double(w_solved_str));%все корни, включая комплексные -- последние надо исключить
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

kappa_x=[-4:10*delta_k:4];
%kappa_x=[-150:20*delta_k:150];


for ind_k=1:1:length(k)
for ind_w=1:1:size(w_solved_str) 
        ind_k
        k_x=k(ind_k);
        omega = w_solved(ind_k,ind_w); %зафиксировали точку по оси к и номер решения уравнения для частоты
        if(w_solved(ind_k,ind_w)~=cr_val) %рассматриваем только действительные \omega
            psi_2=sqrt(I_0/2-(M+beta*k_x.^2)*I_0/(2*(g*I_0+ omega)));
            psi_1=sqrt(I_0/2+(M+beta*k_x.^2)*I_0/(2*(g*I_0+ omega)));
            syms Omega
            parfor ix = 1:1:length(kappa_x) %для каждого каппа составляем матрицу
                Mat=[-Omega-omega+(M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*(psi_1)^2,- g*psi_1^2,k_x+kappa_x(ix),0;
                k_x+kappa_x(ix),0,-Omega - omega  - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2, - g*psi_2^2; 
                - g*psi_1^2,Omega - omega + (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_1^2,0,k_x-kappa_x(ix); 
                0,k_x-kappa_x(ix),- g*psi_2^2,Omega - omega - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2]; 
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


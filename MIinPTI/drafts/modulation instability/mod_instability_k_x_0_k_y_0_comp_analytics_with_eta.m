clear all
close all;

%%1)
I_0=1;
M=1;
g=1;
eta=1;
omega = M-g*I_0;
psi_2=0;
psi_1=sqrt(I_0);
syms Omega x
kappa_x=[-2:0.1:2];
 for ix = 1:1:length(kappa_x) 
                Mat=[-Omega-omega+M - 2*g*(psi_1)^2,- g*psi_1^2,kappa_x(ix)+eta*kappa_x(ix)^2,0;
                kappa_x(ix)+eta*kappa_x(ix)^2,0,-Omega - omega  - M - 2*g*psi_2^2, - g*psi_2^2; 
                - g*psi_1^2,Omega - omega + M - 2*g*psi_1^2,0,-kappa_x(ix)+eta*kappa_x(ix)^2; 
                0,-kappa_x(ix)+eta*kappa_x(ix)^2,- g*psi_2^2,Omega - omega - M - 2*g*psi_2^2]; 
                eqn = det(Mat) == 0; 
                freq(ix,:) = solve(eqn,Omega); 
                frq1(ix,:) = double(freq(ix,:));
                frq(ix,:) = sort(frq1(ix,:));
end

figure(2); hold on;
for n=1:1:4
    plot(kappa_x,imag(frq(:,n)),'om')
end
grid on; box on;

figure(2); hold on;
for n=1:1:4
    plot(kappa_x,real(frq(:,n)),'+r')
end

%%%для сравнения со случаем ета = 0 -- аналитические формулы, полученные
%%%без дисперсии
x1 = sqrt((g*I_0)^2 + sqrt(g*I_0 - 2*M)*sqrt((g*I_0)^3 - 6*(g*I_0)^2*M - 4*(g*I_0).*kappa_x.^2 + 12*(g*I_0)*M^2 ...
- 8.*kappa_x.^2*M - 8*M^3) - 4*(g*I_0)*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
x2 = -sqrt((g*I_0)^2 + sqrt(g*I_0 - 2*M)*sqrt((g*I_0)^3 - 6*(g*I_0)^2*M - 4*(g*I_0).*kappa_x.^2 + 12*(g*I_0)*M^2 ...
- 8.*kappa_x.^2*M - 8*M^3) - 4*(g*I_0)*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
x3 = sqrt((g*I_0)^2 - sqrt(g*I_0 - 2*M)*sqrt((g*I_0)^3 - 6*(g*I_0)^2*M - 4*(g*I_0).*kappa_x.^2 + 12*(g*I_0)*M^2 ...
- 8.*kappa_x.^2*M - 8*M^3) - 4*(g*I_0)*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
x4 = -sqrt((g*I_0)^2 - sqrt(g*I_0 - 2*M)*sqrt((g*I_0)^3 - 6*(g*I_0)^2*M - 4*(g*I_0).*kappa_x.^2 + 12*(g*I_0)*M^2 ...
- 8.*kappa_x.^2*M - 8*M^3) - 4*(g*I_0)*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);


figure(2); hold on; plot(kappa_x,imag(x1),'--m')
plot(kappa_x,imag(x2),'--m')
plot(kappa_x,imag(x3),'--m')
plot(kappa_x,imag(x4),'--m')


figure(2); hold on; plot(kappa_x,real(x1),'--r')
plot(kappa_x,real(x2),'--r')
plot(kappa_x,real(x3),'--r')
plot(kappa_x,real(x4),'--r')


ylabel('$\Omega$')
xlabel('$\kappa_x$')

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
        
%%%в случае ета = 0 зависимость от интенсивности просто построить -- нужно использовать полученные ранее аналитические формулы        
kappa_x=-10:0.5:10;
I_0=0:0.25:3;
for i=1:1:length(I_0)       
    x1 = sqrt((g*I_0(i))^2 + sqrt(g*I_0(i) - 2*M)*sqrt((g*I_0(i))^3 - 6*(g*I_0(i))^2*M - 4*(g*I_0(i)).*kappa_x.^2 + 12*(g*I_0(i))*M^2 ...
    - 8.*kappa_x.^2*M - 8*M^3) - 4*(g*I_0(i))*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
    x2 = -sqrt((g*I_0(i))^2 + sqrt(g*I_0(i) - 2*M)*sqrt((g*I_0(i))^3 - 6*(g*I_0(i))^2*M - 4*(g*I_0(i)).*kappa_x.^2 + 12*(g*I_0(i))*M^2 ...
    - 8.*kappa_x.^2*M - 8*M^3) - 4*(g*I_0(i))*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
    x3 = sqrt((g*I_0(i))^2 - sqrt(g*I_0(i) - 2*M)*sqrt((g*I_0(i))^3 - 6*(g*I_0(i))^2*M - 4*(g*I_0(i)).*kappa_x.^2 + 12*(g*I_0(i))*M^2 ...
    - 8.*kappa_x.^2*M - 8*M^3) - 4*(g*I_0(i))*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
    x4 = -sqrt((g*I_0(i))^2 - sqrt(g*I_0(i) - 2*M)*sqrt((g*I_0(i))^3 - 6*(g*I_0(i))^2*M - 4*(g*I_0(i)).*kappa_x.^2 + 12*(g*I_0(i))*M^2 ...
    - 8.*kappa_x.^2*M - 8*M^3) - 4*(g*I_0(i))*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
    [incr(1),n]=max(imag(x1));
    [incr(2),n]=max(imag(x2));
    [incr(3),n]=max(imag(x3));
    [incr(4),n]=max(imag(x4));
    inc(i)=max(incr);
end

figure(10); plot(I_0,inc,'*')

clear inc incr
%%%в случае ета != 0 зависимость от интенсивности не построить с помощью аналитической формулы -- надо считать для каждой итенсивности в лоб       
for i=1:1:length(I_0)
    i
    omega = M-g*I_0(i);
    psi_2=0;
    psi_1=sqrt(I_0(i));

        for ix = 1:1:length(kappa_x) 
                    Mat=[-Omega-omega+M - 2*g*(psi_1)^2,- g*psi_1^2,kappa_x(ix)+eta*kappa_x(ix)^2,0;
                    kappa_x(ix)+eta*kappa_x(ix)^2,0,-Omega - omega  - M - 2*g*psi_2^2, - g*psi_2^2; 
                    - g*psi_1^2,Omega - omega + M - 2*g*psi_1^2,0,-kappa_x(ix)+eta*kappa_x(ix)^2; 
                    0,-kappa_x(ix)+eta*kappa_x(ix)^2,- g*psi_2^2,Omega - omega - M - 2*g*psi_2^2]; 
                    eqn = det(Mat) == 0; 
                    freq(ix,:) = solve(eqn,Omega); 
                    frq1(ix,:) = double(freq(ix,:));
                    frq(ix,:) = sort(frq1(ix,:));
        end
        %%%для каждого из массива частот (их 4, т к уравнение относительно
        %%%частоты 4 степени) ищем макс инкремент и потом находим максимум
        %%%из них 
                    [incr(1),n]=max(imag(frq(:,1)));
                    [incr(2),n]=max(imag(frq(:,2)));
                    [incr(3),n]=max(imag(frq(:,3)));
                    [incr(4),n]=max(imag(frq(:,4)));
                    inc(i)=max(incr);

end

figure(10); hold on; plot(I_0,inc,'o')
legend ('\eta = 0','\eta =1')
ylabel('$max(Im\Omega)$')
xlabel('$I_0$')

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
        


        
%%%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%
clear all
%%
%%2)
I_0=1;
M=1;
g=1;
omega = -M-g*I_0;
psi_1=0;
psi_2=sqrt(I_0);
k_x=0;
syms Omega x
kappa_x=[-2:0.1:2];
for ix = 1:1:length(kappa_x) 
                Mat=[-Omega-omega+M - 2*g*(psi_1)^2,- g*psi_1^2,k_x+kappa_x(ix),0;
                k_x+kappa_x(ix),0,-Omega - omega  - M - 2*g*psi_2^2, - g*psi_2^2; 
                - g*psi_1^2,Omega - omega + M - 2*g*psi_1^2,0,k_x-kappa_x(ix); 
                0,k_x-kappa_x(ix),- g*psi_2^2,Omega - omega - M - 2*g*psi_2^2]; 
                eqn = det(Mat) == 0; 
                freq(ix,:) = solve(eqn,Omega); 
                frq1(ix,:) = double(freq(ix,:));
                frq(ix,:) = sort(frq1(ix,:));
end

figure(4); hold on;
for n=1:1:4
    plot(kappa_x,imag(frq(:,n)),'om')
end
grid on; box on;

figure(4); hold on;
for n=1:1:4
    plot(kappa_x,real(frq(:,n)),'+r')
end

x1 = sqrt((g*I_0)^2 - sqrt((g*I_0) + 2*M)*sqrt((g*I_0)^3 + 6*(g*I_0)^2*M - 4*(g*I_0).*kappa_x.^2 ...
    + 12*(g*I_0)*M^2 + 8.*kappa_x.^2*M + 8*M^3) + 4*(g*I_0)*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
x2 = -sqrt((g*I_0)^2 - sqrt((g*I_0) + 2*M)*sqrt((g*I_0)^3 + 6*(g*I_0)^2*M - 4*(g*I_0).*kappa_x.^2 ...
    + 12*(g*I_0)*M^2 + 8.*kappa_x.^2*M + 8*M^3) + 4*(g*I_0)*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
x3 = sqrt((g*I_0)^2 + sqrt((g*I_0) + 2*M)*sqrt((g*I_0)^3 + 6*(g*I_0)^2*M - 4*(g*I_0).*kappa_x.^2 ...
    + 12*(g*I_0)*M^2 + 8.*kappa_x.^2*M + 8*M^3) + 4*(g*I_0)*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
x4 = -sqrt((g*I_0)^2 + sqrt((g*I_0) + 2*M)*sqrt((g*I_0)^3 + 6*(g*I_0)^2*M - 4*(g*I_0).*kappa_x.^2 ...
    + 12*(g*I_0)*M^2 + 8.*kappa_x.^2*M + 8*M^3) + 4*(g*I_0)*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);



figure(4); hold on; plot(kappa_x,imag(x1),'--m')
plot(kappa_x,imag(x2),'--m')
plot(kappa_x,imag(x3),'--m')
plot(kappa_x,imag(x4),'--m')


figure(4); hold on; plot(kappa_x,real(x1),'--r')
plot(kappa_x,real(x2),'--r')
plot(kappa_x,real(x3),'--r')
plot(kappa_x,real(x4),'--r')


ylabel('$\Omega$')
xlabel('$\kappa_x$')

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
        
        
kappa_x=-5:0.01:50;
I_0=0:0.1:130;
for i=1:1:length(I_0)       
    x1 = sqrt((g*I_0(i))^2 - sqrt((g*I_0(i)) + 2*M)*sqrt((g*I_0(i))^3 + 6*(g*I_0(i))^2*M - 4*(g*I_0(i)).*kappa_x.^2 ...
    + 12*(g*I_0(i))*M^2 + 8.*kappa_x.^2*M + 8*M^3) + 4*(g*I_0(i))*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
    x2 = -sqrt((g*I_0(i))^2 - sqrt((g*I_0(i)) + 2*M)*sqrt((g*I_0(i))^3 + 6*(g*I_0(i))^2*M - 4*(g*I_0(i)).*kappa_x.^2 ...
    + 12*(g*I_0(i))*M^2 + 8.*kappa_x.^2*M + 8*M^3) + 4*(g*I_0(i))*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
    x3 = sqrt((g*I_0(i))^2 + sqrt((g*I_0(i)) + 2*M)*sqrt((g*I_0(i))^3 + 6*(g*I_0(i))^2*M - 4*(g*I_0(i)).*kappa_x.^2 ...
    + 12*(g*I_0(i))*M^2 + 8.*kappa_x.^2*M + 8*M^3) + 4*(g*I_0(i))*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
    x4 = -sqrt((g*I_0(i))^2 + sqrt((g*I_0(i)) + 2*M)*sqrt((g*I_0(i))^3 + 6*(g*I_0(i))^2*M - 4*(g*I_0(i)).*kappa_x.^2 ...
    + 12*(g*I_0(i))*M^2 + 8.*kappa_x.^2*M + 8*M^3) + 4*(g*I_0(i))*M + 2.*kappa_x.^2 + 4*M^2)/sqrt(2);
    [incr(1),n]=max(imag(x1));
    [incr(2),n]=max(imag(x2));
    [incr(3),n]=max(imag(x3));
    [incr(4),n]=max(imag(x4));
    inc(i)=max(incr);
end

figure; plot(I_0,inc)

ylabel('$max(Im\Omega)$')
xlabel('$I_0$')

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
        
%%%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%       
clear all
%%
%%3)
I_0=50;
M=1;
g=1;
omega =-g*I_0/2;
psi_1=sqrt(I_0/2+M/g);
psi_2=sqrt(I_0/2-M/g);
k_x=0;
syms Omega x
kappa_x=[-20:1:20];
for ix = 1:1:length(kappa_x) 
                Mat=[-Omega-omega+M - 2*g*(psi_1)^2,- g*psi_1^2,k_x+kappa_x(ix),0;
                k_x+kappa_x(ix),0,-Omega - omega  - M - 2*g*psi_2^2, - g*psi_2^2; 
                - g*psi_1^2,Omega - omega + M - 2*g*psi_1^2,0,k_x-kappa_x(ix); 
                0,k_x-kappa_x(ix),- g*psi_2^2,Omega - omega - M - 2*g*psi_2^2]; 
                eqn = det(Mat) == 0; 
                freq(ix,:) = solve(eqn,Omega); 
                frq1(ix,:) = double(freq(ix,:));
                frq(ix,:) = sort(frq1(ix,:));
end

figure(4); hold on;
for n=1:1:4
    plot(kappa_x,imag(frq(:,n)),'om')
end
grid on; box on;

figure(4); hold on;
for n=1:1:4
    plot(kappa_x,real(frq(:,n)),'+r')
end
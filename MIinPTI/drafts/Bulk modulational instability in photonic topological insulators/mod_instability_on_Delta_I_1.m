clear all
close all;


% G=-1;
% J_2= -1:0.25:1;
% J_1=1;
% I_0_mat = 10:0.5:15;

G=-1;
J_2= -1:0.5:1;
J_1=1;
I_0_mat = 5:1:10;

Delta = -5;

for j = 1:1:length(I_0_mat)
for d=1:1:length(J_2)
    map_to_det_area_of_cross(d,j)=I_0_mat(j)-2*(Delta-4*J_2(d))/(-G);
    map_to_det_area_of_cross_1(d,j)=I_0_mat(j)+2*(Delta-4*J_2(d))/(-G);
    if(map_to_det_area_of_cross(d,j))>0 && (map_to_det_area_of_cross_1(d,j))>0
        sgn(d,j)=1;
    else
        sgn(d,j)=-1;
    end
end
end



p = panel();
p.pack(2, 2);

p(1, 1).select();
title('(a)~~The existence of the point of the cross')
imagesc(I_0_mat,J_2,sgn); box on;
xlabel('$I_0$')
ylabel('$J_2$')

for j = 1:1:length(I_0_mat)
for d=1:1:length(J_2)
    
    I_0 = I_0_mat(j);
    g=-G;
    M = Delta - 4*J_2(d);
    beta=J_2(d);
 
    %%%верхн€€ ветка
    psi_2=0;
    psi_1=sqrt(I_0);
    omega = M-g*I_0;
    psi_2=0;
    psi_1=sqrt(I_0);
    k_x=0;
    syms Omega x
    kappa_x= - sqrt(2).*J_1.*[-10:0.5:10];
    parfor ix = 1:1:length(kappa_x) %дл€ каждого каппа составл€ем матрицу
                Mat=[-Omega-omega+(M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*(psi_1)^2,- g*psi_1^2,k_x+kappa_x(ix),0;
                k_x+kappa_x(ix),0,-Omega - omega  - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2, - g*psi_2^2; 
                - g*psi_1^2,Omega - omega + (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_1^2,0,k_x-kappa_x(ix); 
                0,k_x-kappa_x(ix),- g*psi_2^2,Omega - omega - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2]; 
                eqn = det(Mat) == 0; 
                freq(ix,:) = solve(eqn,Omega); 
                incr(ix) = (max(imag(freq(ix,:))));                   
    end
    incr_up_band(d,j) = (double(max(incr)));
    
    
    %%%нижн€€ ветка
    omega = -M-g*I_0;
    psi_1=0;
    psi_2=sqrt(I_0);
    k_x=0;
    syms Omega x
    kappa_x= - sqrt(2).*J_1.*[-10:0.5:10];
    parfor ix = 1:1:length(kappa_x) %дл€ каждого каппа составл€ем матрицу
                Mat=[-Omega-omega+(M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*(psi_1)^2,- g*psi_1^2,k_x+kappa_x(ix),0;
                k_x+kappa_x(ix),0,-Omega - omega  - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2, - g*psi_2^2; 
                - g*psi_1^2,Omega - omega + (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_1^2,0,k_x-kappa_x(ix); 
                0,k_x-kappa_x(ix),- g*psi_2^2,Omega - omega - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2]; 
                eqn = det(Mat) == 0; 
                freq(ix,:) = solve(eqn,Omega); 
                incr(ix) = (max(imag(freq(ix,:))));                   
    end
    incr_l_band(d,j) = (double(max(incr)));
    
    
    %%%крест
    omega =-g*I_0/2;
    if (I_0/2+M/g)>0 && (I_0/2-M/g)>0
        psi_1=sqrt(I_0/2+M/g);
        psi_2=sqrt(I_0/2-M/g);
        syms Omega x
        kappa_x= - sqrt(2).*J_1.*[-10:0.5:10];
        parfor ix = 1:1:length(kappa_x) %дл€ каждого каппа составл€ем матрицу
                    Mat=[-Omega-omega+(M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*(psi_1)^2,- g*psi_1^2,k_x+kappa_x(ix),0;
                    k_x+kappa_x(ix),0,-Omega - omega  - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2, - g*psi_2^2; 
                    - g*psi_1^2,Omega - omega + (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_1^2,0,k_x-kappa_x(ix); 
                    0,k_x-kappa_x(ix),- g*psi_2^2,Omega - omega - (M+beta*k_x^2+beta*kappa_x(ix)^2) - 2*g*psi_2^2]; 
                    eqn = det(Mat) == 0; 
                    freq(ix,:) = solve(eqn,Omega); 
                    incr(ix) = (max(imag(freq(ix,:))));                   
        end
        incr_cross_band(d,j) = (double(max(incr)));
    else
        incr_cross_band(d,j) = -1;
    end

end
end


p(1, 2).select();
title('(b)~~The instability of the upper band')
imagesc(I_0_mat,J_2,incr_up_band); box on;
xlabel('$I_0$')
ylabel('$J_2$')

p(2, 1).select();
title('(c)~~The instability of the lower band')
imagesc(I_0_mat,J_2,incr_l_band); box on;
xlabel('$I_0$')
ylabel('$J_2$')

p(2, 2).select();
title('(d)~~The instability of the cross point')
imagesc(I_0_mat,J_2,incr_cross_band); box on; 
xlabel('$I_0$')
ylabel('$J_2$')


size=10;
a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

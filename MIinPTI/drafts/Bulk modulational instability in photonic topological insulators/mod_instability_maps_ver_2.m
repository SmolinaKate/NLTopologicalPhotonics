clear all
close all;


G=0:0.25:4;
J_2= 1;
J_1=sqrt(2);
I_0 = 1;
Delta = 0:0.25:5;

for j = 1:1:length(G)
for d=1:1:length(Delta)
    map_to_det_area_of_cross(d,j)=I_0-2*(Delta(d)-4*J_2)/(-G(j));
    map_to_det_area_of_cross_1(d,j)=I_0+2*(Delta(d)-4*J_2)/(-G(j));
    if(map_to_det_area_of_cross(d,j))>=0 && (map_to_det_area_of_cross_1(d,j))>=0
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
imagesc(G,Delta,sgn); box on;
ylabel('$\Delta$')
xlabel('$G$')

for j = 1:1:length(G)
for d=1:1:length(Delta)
    
    g=-G(j);
    M = Delta(d) - 4*J_2;
    beta=J_2;
 
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
    k_x=0;
    if (I_0/2+M/g)>=0 && (I_0/2-M/g)>=0
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
imagesc(G,Delta,incr_up_band); box on;
ylabel('$\Delta$')
xlabel('$G$')


p(2, 1).select();
title('(c)~~The instability of the lower band')
imagesc(G,Delta,incr_l_band); box on;
ylabel('$\Delta$')
xlabel('$G$')


p(2, 2).select();
title('(d)~~The instability of the cross point')
imagesc(G,Delta,incr_cross_band); box on;
ylabel('$\Delta$')
xlabel('$G$')



size=10;
a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

clear all;
close all;

N=4;
p = panel();
p.pack(2, N);
%solve (determinant[{{g-x+b,g,-k,0};{-k,0,-x-b+C,0};{g,x+g+b,0,k};{0,k,0,x-b+C}}]=0) for x
step_p=0.06;
G=-4:step_p:4;
J_2=1/sqrt(2);
J_1=1;
I_0 = 1;

kappa=[-5:0.0001:5];
g=G.*I_0;
b=kappa.^2.*J_2;
k=J_1*sqrt(2).*kappa;


C=0;
for i=1:1:length(g)
    x1 =sqrt(-sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i).^2 - 4.*b.*C.^3 ...
        - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2) + 2.*b.^2 - 2.*b.*C + 2.*b.*g(i) + C.^2 + 2.*k.^2)./sqrt(2);
    x2 = sqrt(1/2.*sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i)^2 - 4.*b.*C.^3 ...
        - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2) + b.^2 - b.*C + b.*g(i) + C.^2/2 + k.^2);
    [x1_im,n1]=(max(imag(x1)));
    [x2_im,n2]=(max(imag(x2)));
    if x1_im >= x2_im
        ar_to_plor_kappa(i)=1;
        inc(i)=x1_im;
        kap_inst(i)=kappa(n1);
    else
        ar_to_plor_kappa(i)=1;
        inc(i)=x2_im;
        kap_inst(i)=kappa(n2);
    end
    
    if(x1_im==0)&&(x2_im==0)
        ar_to_plor_kappa(i)=0;
    end
end
p(1, 1).select();
title(['(a)~~$C=0$'])
hold on; plot(G,inc,'r');box on;
hold on; plot([-sqrt(2) -sqrt(2)], [0 max(inc)],'--b')
% for q=1:1:length(g)
%  if (-(G(q)*I_0+J_1^2/J_2)>0)
%    hold on; plot(G(q),-(G(q)*I_0+J_1^2/J_2),'.black')
%  end
% end
xlabel('$\Gamma$')
ylabel('$\Im(\lambda)$')

p(2, 1).select();
for h=1:1:length (G)
   if ar_to_plor_kappa(h)>0
 %       hold on; plot(G(h),-(sqrt(abs(G(h).*I_0*J_2+J_1^2)/J_2^2)),'oc')
        hold on; plot(G(h),kap_inst(h),'.r');box on;
   else
        hold on; plot(G(h),0,'.b');box on;
   end    
end
xlabel('$\Gamma$')
ylabel('$\kappa_y^{max}$')




















clear inc
C=8*J_2;
for i=1:1:length(g)
    x1 =sqrt(-sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i).^2 - 4.*b.*C.^3 ...
        - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2) + 2.*b.^2 - 2.*b.*C + 2.*b.*g(i) + C.^2 + 2.*k.^2)./sqrt(2);
    x2 = sqrt(1/2.*sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i)^2 - 4.*b.*C.^3 ...
        - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2) + b.^2 - b.*C + b.*g(i) + C.^2/2 + k.^2);
    [x1_im,n1]=(max(imag(x1)));
    [x2_im,n2]=(max(imag(x2)));
    if x1_im >= x2_im
        ar_to_plor_kappa(i)=1;
        inc(i)=x1_im;
        kap_inst(i)=kappa(n1);
    else
        ar_to_plor_kappa(i)=1;
        inc(i)=x2_im;
        kap_inst(i)=kappa(n2);
    end
    
    if(x1_im==0)&&(x2_im==0)
        ar_to_plor_kappa(i)=0;
    end
end
p(1, 2).select();
title(['(b)~~$C=8 J_2$'])
hold on; plot(G,inc,'r')
box on;
xlabel('$\Gamma$')
ylabel('$\Im(\lambda)$')
p(2, 2).select();
for h=1:1:length (G)
   if ar_to_plor_kappa(h)>0
        hold on; plot(G(h),kap_inst(h),'.r');box on;
   else
        hold on; plot(G(h),0,'.b');box on;
   end    
end
xlabel('$\Gamma$')
ylabel('$\kappa_y^{max}$')


















clear inc
C=2*J_1^2/J_2;
for i=1:1:length(g)
    x1 =sqrt(-sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i).^2 - 4.*b.*C.^3 ...
        - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2) + 2.*b.^2 - 2.*b.*C + 2.*b.*g(i) + C.^2 + 2.*k.^2)./sqrt(2);
    x2 = sqrt(1/2.*sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i)^2 - 4.*b.*C.^3 ...
        - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2) + b.^2 - b.*C + b.*g(i) + C.^2/2 + k.^2);
    [x1_im,n1]=(max(imag(x1)));
    [x2_im,n2]=(max(imag(x2)));
    if x1_im >= x2_im
        ar_to_plor_kappa(i)=1;
        inc(i)=x1_im;
        kap_inst(i)=kappa(n1);
    else
        ar_to_plor_kappa(i)=1;
        inc(i)=x2_im;
        kap_inst(i)=kappa(n2);
    end
    
    if(x1_im==0)&&(x2_im==0)
        ar_to_plor_kappa(i)=0;
    end
end
p(1, 3).select();
title(['(c)~~$C=2 J_1^2/J_2$'])
hold on; plot(G,inc,'r')
box on;
xlabel('$\Gamma$')
ylabel('$\Im(\lambda)$')
p(2, 3).select();
for h=1:1:length (G)
   if ar_to_plor_kappa(h)>0
        hold on; plot(G(h),kap_inst(h),'.r');box on;
   else
        hold on; plot(G(h),0,'.b');box on;
   end    
end
xlabel('$\Gamma$')
ylabel('$\kappa_y^{max}$')
















clear inc
C=-2*J_1^2/J_2;
for i=1:1:length(g)
    x1 =sqrt(-sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i).^2 - 4.*b.*C.^3 ...
        - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2) + 2.*b.^2 - 2.*b.*C + 2.*b.*g(i) + C.^2 + 2.*k.^2)./sqrt(2);
    x2 = sqrt(1/2.*sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i)^2 - 4.*b.*C.^3 ...
        - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2) + b.^2 - b.*C + b.*g(i) + C.^2/2 + k.^2);
    [x1_im,n1]=(max(imag(x1)));
    [x2_im,n2]=(max(imag(x2)));
    if x1_im >= x2_im
        ar_to_plor_kappa(i)=1;
        inc(i)=x1_im;
        kap_inst(i)=kappa(n1);
    else
        ar_to_plor_kappa(i)=1;
        inc(i)=x2_im;
        kap_inst(i)=kappa(n2);
    end
    
    if(x1_im==0)&&(x2_im==0)
        ar_to_plor_kappa(i)=0;
    end
end
p(1, 4).select();
title(['(d)~~$C=-2 J_1^2/J_2$'])
hold on; plot(G,inc,'r')
box on;
xlabel('$\Gamma$')
ylabel('$\Im(\lambda)$')
p(2, 4).select();
for h=1:1:length (G)
   if ar_to_plor_kappa(h)>0
        hold on; plot(G(h),kap_inst(h),'.r');box on;
   else
        hold on; plot(G(h),0,'.b');box on;
   end    
end
xlabel('$\Gamma$')
ylabel('$\kappa_y^{max}$')



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


% checking formulas
% g(i)=1;
% G(i)=1;
% C=3;
%  x1 =sqrt(-sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i).^2 - 4.*b.*C.^3 ...
%         - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2) + 2.*b.^2 - 2.*b.*C + 2.*b.*g(i) + C.^2 + 2.*k.^2)./sqrt(2);
%     x2 = sqrt(1/2.*sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i)^2 - 4.*b.*C.^3 ...
%         - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2) + b.^2 - b.*C + b.*g(i) + C.^2/2 + k.^2);
%         sq=sqrt(4.*b.^2.*C^2 + 8.*b.^2.*C.*g(i) + 4.*b.^2.*g(i).^2 - 4.*b.*C.^3 ...
%         - 4.*b.*C.^2.*g(i) + C.^4 + 4.*C.^2.*k.^2 + 8.*C.*g(i).*k.^2);
% 
%     
%  F=2*sqrt(C^4/4+(C+G(i)*I_0)^2*J_2^2.*kappa.^4 + kappa.^2.*(-J_2*C^2*(C+G(i)*I_0)+4*J_1*C*G(i)*I_0+2*C^2*J_1));
%     
%  D=J_2^2.*kappa.^4+kappa.^2.*(J_2*I_0*G(i)+2*J_1^2-C*J_2)+C^2/2;
%  
%  l_1=sqrt(-F+2.*D)./sqrt(2);
%  
%  l_2=sqrt(F+2.*D)./sqrt(2);
%  close all
% 
%  figure; plot(F,'b'); hold on; plot(F,'--r');  
%   figure; plot(imag(l_1),'b'); hold on; plot(imag(x1),'--r');  
%    figure; plot(l_2,'b'); hold on; plot(x2,'--r');  
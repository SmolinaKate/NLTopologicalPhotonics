clear 
close all

J1=1;
J2=J1/sqrt(2);
Delta=J1*2*sqrt(2);

kx=-2*pi:0.01:2*pi;
ky= kx;
omega_an=zeros(length(ky),length(ky));
omega_an1=zeros(length(ky),length(ky));
for j=1:1:length(ky)
omega_an_plus(j,:)= +sqrt(4*J1^2*(1+cos(kx)*cos(ky(j)))+(Delta+2*J2.*(cos(kx)-cos(ky(j)))).^2);
omega_an_minus(j,:)=-sqrt(4*J1^2*(1+cos(kx)*cos(ky(j)))+(Delta+2*J2.*(cos(kx)-cos(ky(j)))).^2);
end

figure(1);
hold on
surfc(kx,ky,omega_an_plus,'FaceAlpha',0.6);shading interp;
surf(kx,ky,omega_an_minus);shading interp;
%plot3([x1,x2,x3,x4,x1],[y1,y2,y3,y4,y1],[z1,z2,z3,z4,z1])
%plot3(pi*[0.5,0.5,-0.5,-0.5],pi*[0.5,-0.5,-0.5,0.5],pi*[0,0,0,0],'-k','LineWidth',0.5)
plot3(pi*[1,1,-1,-1,1],pi*[1,-1,-1,1,1],pi*[0,0,0,0,0],'-k','LineWidth',0.5)
scatter3(pi,0,0, 'MarkerFaceColor',[0 0 0])
scatter3(pi,pi,0, 'MarkerFaceColor',[0 0 0])
scatter3(0,pi,0, 'MarkerFaceColor',[0 0 0])
scatter3(0,0,0, 'MarkerFaceColor',[0 0 0])

view(80,25)

zlabel('$E_L(k_x,k_y,\Delta)$');
xlabel('$k_x$')
ylabel('$k_y$')


text(0,0,0.5,'$\Gamma$')
text(pi,0,0.5,'$X$')
text(0,pi,0.5,'$Y$')
text(pi,pi,0.5,'$M$')
           
sizeFont=15;
a = findobj(gcf); 
allaxes  = findall(a,'Type','axes');
alltext  = findall(a,'Type','text');
set(allaxes,'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');




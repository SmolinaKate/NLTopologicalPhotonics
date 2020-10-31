clear all; close all; 
addpath('C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster')

% eta=0.1;
mu=[0:0.001:1.5];
Lambda=[0:0.01:10];
for im=1:1:length(mu)
    for il=1:1:length(Lambda)
 %       T(im,il)=2*sqrt(2.718281828459)*real(sqrt(eta./mu(im).^3 ...
  %      -eta^2./Lambda(il).^2./mu(im).^4-1/64));
       T(il,im)=sqrt(2.718281828459)*real(sqrt(1./mu(im).^3 ...
       -1./Lambda(il).^2./mu(im).^4-1))/4;
    end
end

[Lambdag,mug]=meshgrid(Lambda,mu);
% z= eta./mug.^3 ...
%         -eta^2./Lambdag.^2./mug.^4-1/64-1/4/2.718281828459;
z=(2.718281828459)*((1./mug.^3 ...
       -1./Lambdag.^2./mug.^4-1))/4^2 - 1;

v = [0,0];

figure; 
imagesc( Lambda,mu, T'); colorbar
hold on; contour(Lambda, mu, z,v,'--r'); box on;
% imagesc(T); colorbar
ylabel('$\tilde{\mu}$')
xlabel('$\tilde{\Lambda}$')

text(10.5,1.55,'$t^*/ T_{MI}$');

a = findobj(gcf); 
size = 16;
allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alllines,'LineWidth',2);
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(gca,'YDir','normal')
str=['Fig_S2_ver_2','.png'];
%print('-dpng','-r500',str)
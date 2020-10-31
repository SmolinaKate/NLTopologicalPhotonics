clear all; close all; 
addpath('C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster')
link=1;
a=abs(2/link/sqrt(3));
dt =1;
eta = link*a^2/8;
M = link/3;
mu_in_2=eta*M;
g=mu_in_2*3*M;
Omega =M^2*eta-0.01;
t_end=300;
f = @(t,Y) [Y(2); -(Omega./eta-M^2).*Y(1)+ g^4/96^2/M^4.* Y(1).^9./eta^2- g/4/eta.*Y(1).^3];
y1 = [0:0.001:0.4];
y2 = [0:0.001:0.05];
syms x
eqn =-(Omega./eta-M^2).*x+ g^4/96^2/M^4.* x.^9./eta^2- g/4/eta.*x.^3==0;
V = vpasolve(eqn,x,[-1000 1000]);

[x,y] = meshgrid(y1,y2);
size(x);
size(y);
u = zeros(size(x));
v1 = zeros(size(x));

t=0;
for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v1(i) = Yprime(2);
end


figure('Position', [10 50 800 300]);
p = panel();
border =15; 
p.pack(2, 4);
p(1).repack(0.001);
p(2,1).repack(0.01);

p(2, 2).marginright = border; 
p(2, 2).marginleft = border; 
p(2, 3).marginright = border; 
p(2, 3).marginleft = border; 
p(2, 4).marginright = border; 
p(2, 4).marginleft = border; 


p(2, 3).select();
hold on;plot([V(3),V(3)],[0,0],'oc')
hold on;plot([V(4),V(4)],[0,0],'oc')
% hold on;plot([V(5),V(5)],[0,0],'oc')
%h1=quiver(x,y,u,v1,'black'); figure(gcf)
% set(h1,'AutoScale','on', 'AutoScaleFactor', 0.5)
box on; 
xlabel('$\mathcal{A}$')
% ylabel('$d \mathcal{A}/d \zeta$')
hold on


y20 =  10^(-18);
y21 = sqrt(-Omega/eta+M^2)*y20;
[ts,ys] = ode113(f,[0,t_end],[y20,y21]);
max_f=max(ys(:,1));
plot(ys(:,1),ys(:,2))
hold off

 
 p(2, 2).select();
 clear f U
 f=[-1.5*max_f:max_f/100:1.5*max_f];       
 U=-( -(Omega./eta-M^2).*f.^2./2+ ( g^4/96^2/M^4.* f.^10./eta^2)./10- g/4/eta.*f.^4./4);
 plot(f,U,'c');
%  hold on; plot([max_f,max_f],[min(U),max(U)],'--c')
%  hold on; plot([-max_f,-max_f],[min(U),max(U)],'--c')
 xlim( [ min(f), max(f)] );
 ylim([ min(U), max(U)] );
 box on; 
 xlabel('$\mathcal{A}$')
 ylabel('$U$')

 

p(2, 4).select();
% plot(ts,ys)
cutn=13;
zeta=ts(cutn:end)-ts(cutn);
Sol_num=ys(cutn:end,1);
[m,zeta_max]=max(Sol_num);
Sol_S=sqrt(8*(M^2*eta-Omega)/g)./cosh((zeta-zeta(zeta_max)).*sqrt(M^2-Omega/eta));
scale=sqrt(M^2-Omega/eta);
plot(Sol_num,(zeta-zeta(zeta_max)).*scale);
hold on; plot(Sol_S,(zeta-zeta(zeta_max)).*scale,'--black')
box on; 
xlim( [ min(Sol_num), max(Sol_num)] );
ylim([ min((zeta-zeta(zeta_max)).*scale), max((zeta-zeta(zeta_max)).*scale)] );
%ylabel('$\zeta \sqrt{M_0^2 - \omega_s/\eta}$')
xlabel('$\mathcal{A}$')

h=text(-1.87,3.5-zeta(zeta_max).*scale,'$d \mathcal{A}/d \zeta$');
set(h,'Rotation',90);
h=text(-0.25,1-zeta(zeta_max).*scale,'$\zeta \sqrt{M_0^2 - \omega_s/\eta}$');
set(h,'Rotation',90);
text(-2.2,1.5+10-zeta(zeta_max).*scale,'$(a)$');
text(1,1.5+10-zeta(zeta_max).*scale,'$(c)$');
text(-0.6,1.5+10-zeta(zeta_max).*scale,'$(b)$');
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

str=['Fig_S_2','.png'];

g*max_f^2/M
eta*M
%print('-dpng','-r500',str)
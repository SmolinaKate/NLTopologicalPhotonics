clear all; close all; 
addpath('C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster')
link=1;
a=abs(2/link/sqrt(3));
dt =1;
eta = link*a^2/8;
M = link/3;
g=eta*3;
Omega =M^2*eta-0.01;
t_end=340;
f = @(t,Y) [Y(2); -(Omega./eta-M^2).*Y(1)+ ( g^4/96^2/M^4.* Y(1).^4./eta^2).^2.*Y(1)- g/4/eta.*Y(1).^3];
y1 = [0:0.01:0.4];
y2 = [-0.05:0.001:0.05];
syms x
eqn =-(Omega./eta-M^2).*x+ ( g^4/96^2/M^4.* x.^4./eta^2).^2.*x- g/4/eta.*x.^3==0;
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
p.pack(1, 5);
p(1,1).repack(0.01);

p(1, 2).marginright = border; 
p(1, 2).marginleft = border; 
p(1, 3).marginright = border; 
p(1, 3).marginleft = border; 
p(1, 4).marginright = border; 
p(1, 4).marginleft = border; 
p(1, 5).marginright = border; 
p(1, 5).marginleft = border; 


p(1, 2).select();
hold on;plot([V(3),V(3)],[0,0],'oc')
hold on;plot([V(4),V(4)],[0,0],'oc')
% hold on;plot([V(5),V(5)],[0,0],'oc')
%h1=quiver(x,y,u,v1,'black'); figure(gcf)
% set(h1,'AutoScale','on', 'AutoScaleFactor', 0.5)
box on; 
xlabel('$\mathcal{A}$')
ylabel('$d \mathcal{A}/d \zeta$')
hold on


y20 =  10^(-30);
y21 = sqrt(-Omega/eta+M^2)*y20;
[ts,ys] = ode45(f,[0,t_end],[y20,y21]);
max_f=max(ys(:,1));
plot(ys(:,1),ys(:,2))
hold off


 p(1, 4).select();
 clear f U
 f=[0:0.01:V(5)+1];       
 U=-( -(Omega./eta-M^2).*f.^2./2+ ( g^4/96^2/M^4.* f.^4./eta^2).^2.*f.^2./10- g/4/eta.*f.^4./4);
 plot(f,U,'c');
 box on; 
 xlabel('$\mathcal{A}$')
 ylabel('$U$')


 
 p(1, 3).select();
 clear f U
 f=[-1.5*max_f:max_f/100:1.5*max_f];       
 U=-( -(Omega./eta-M^2).*f.^2./2+ ( g^4/96^2/M^4.* f.^4./eta^2).^2.*f.^2./10- g/4/eta.*f.^4./4);
 plot(f,U,'c');
%  hold on; plot([max_f,max_f],[min(U),max(U)],'--c')
%  hold on; plot([-max_f,-max_f],[min(U),max(U)],'--c')
 xlim( [ min(f), max(f)] );
 ylim([ min(U), max(U)] );
 box on; 
 xlabel('$\mathcal{A}$')
 ylabel('$U$')

 

p(1, 5).select();
plot( ys(50:end,1),ts(50:end)-ts(50))
box on; 
ylabel('$\zeta$')
xlabel('$\mathcal{A}$')




a = findobj(gcf); 
size = 12;

allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alllines,'LineWidth',2);
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

str=['phase_por_1','.png'];
print('-dpng','-r500',str)
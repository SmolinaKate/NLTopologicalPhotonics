clear all; close all; 

eta=0.1667;
link=1;

g=eta*10;
M = link/10;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % Omega = - 2;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % v=1.158;


Omega = -2;
v=1.1558;

%%%%%%%%%%%%% числ решение уравнения
syms x
eqn =Omega/eta.*x+( v/2/eta -  g^2/16/eta/M^2/6.*x.^4).^2.*x  -g/4/eta.*x.^3==0;
V = vpasolve(eqn,x,[-1000 1000])
%%%%%%%%%%%%%


f = @(t,Y) [Y(2); Omega./eta.*Y(1)+ ( v./2./eta -  g^2/16/M^2.* Y(1).^4/6./eta).^2.*Y(1)- g/4/eta.*Y(1).^3];


y1 = [-0.9:0.01:0.9];
y2 = [-0.9:0.01:0.9];

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


figure;
h1=quiver(x,y,u,v1,'black'); figure(gcf)
xlabel('f')
ylabel('df')
set(h1,'AutoScale','on', 'AutoScaleFactor', 20)
hold on


y20 =  10^(-20);
[ts,ys] = ode45(f,[0,355],[y20,0]);
plot(ys(:,1),ys(:,2))
hold off

              

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

        
        
figure; plot(ts, ys(:,1))
xlabel('$\xi$')
ylabel('f')


% % % Ed(1:length(squeeze(ys(:,1)))) = 1;
% % % sf=v./2./eta.*Ed - g.^2/16/eta/M^2/6.*squeeze(ys(:,1)').^2; 
% % % dt = abs(ts(1) - ts(2));
% % % int=cumsum(sf)*dt; 
% % % 
% % % 
% % % 
% % % 
% % % ys_new = ys(:,1).*exp(1i.*int);

path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster\1_soliton';
file_name = sprintf('soliton');
ToSaveFileName = [file_name '.mat'];
ToSaveFileName1 = fullfile(path, ToSaveFileName);
save(ToSaveFileName1,'ts','ys');  

a = findobj(gcf); 

allaxes  = findall(a,'Type','axes');
alllines = findall(a,'Type','line');
alltext  = findall(a,'Type','text');

set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alllines,'LineWidth',1);
set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');

        
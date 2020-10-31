
% linear case
M=1;
k=[-5:0.01:5];
omega_1=-k;
omega_2=k;
omega_3=sqrt(M^2+k.^2);
omega_4=-sqrt(M^2+k.^2);
figure;
plot(k,omega_1,'cyan','LineWidth',2);
hold on; plot(k,omega_2,'cyan','LineWidth',2); 
hold on; plot(k,omega_3,'b','LineWidth',2); 
hold on; plot(k,omega_4,'b','LineWidth',2); 
grid on;
ylabel('$\omega$','Interpreter', 'latex','FontSize', 15)
xlabel('$k$','Interpreter', 'latex','FontSize', 15)
xticks([-5:1:5])

% nonlinear case
I_0=3;
I_1=I_0/2;
M=1;
g=1;

omega_1=-g*I_0/2;
omega_2=-M - g*I_0;
omega_3=M - g*I_0;
[w,k]=meshgrid([-10:0.05:10],[-10:0.05:10]);
%z= 4*w.^4 + w.^3*(12*g*I_0) + w.^2.*(13*g^2*I_0^2 - 4*M^2 -4.*k.^2)+w.*(6*g^3*I_0^3 - 4*M^2*g*I_0 - 8*k^2*g*I_0)+g^4*I_0^4-M^2*g^2*I_0^2-k.^2*4*g^2*I_0^2;
z= (w+g*I_0/2).^2 - k.^2 - (M-g*I_0*M./2./(w+g*I_0)).^2;
z_1=w+g/2*I_1+k;
z_2=w+g/2*I_1-k;
v = [0,0];
figure(2);
subplot(2,2,1)
contour(k,w,z,v,'b','LineWidth',2);
% hold on; contour(k,w,z_2,v,'cyan','LineWidth',2);
% hold on; plot (0,omega_2,'ko');
% hold on; plot (0,omega_3,'ko');
axis equal
grid on;
xticks([-10:1:10])
yticks([-10:1:10])
ylabel('$\omega$','Interpreter', 'latex','FontSize', 15)
xlabel('$k$','Interpreter', 'latex','FontSize', 15)
title(['I_0 = ',num2str(I_0),';I_1 = ',num2str(I_1),';M = ',num2str(M),';g = ',num2str(g)]);



I_0=1;
I_1=I_0/2;
M=1;
g=1;

omega_1=-g*I_0/2;
omega_2=-M - g*I_0;
omega_3=M - g*I_0;
[w,k]=meshgrid([-10:0.05:10],[-10:0.05:10]);
z= (w+g*I_0/2).^2 - k.^2 - (M-g*I_0*M./2./(w+g*I_0)).^2;
z_1=w+g/2*I_1+k;
z_2=w+g/2*I_1-k;
v = [0,0];
figure(2);
subplot(2,2,2)
contour(k,w,z,v,'b','LineWidth',2);
% hold on; contour(k,w,z_1,v,'cyan','LineWidth',2);
% hold on; contour(k,w,z_2,v,'cyan','LineWidth',2);
% hold on; plot (0,omega_2,'ko');
% hold on; plot (0,omega_3,'ko');
axis equal
grid on;
xticks([-10:1:10])
yticks([-10:1:10])
ylabel('$\omega$','Interpreter', 'latex','FontSize', 15)
xlabel('$k$','Interpreter', 'latex','FontSize', 15)
title(['I_0 = ',num2str(I_0),';I_1 = ',num2str(I_1),';M = ',num2str(M),';g = ',num2str(g)]);



I_0=3;
I_1=I_0/2;
M=1;
g=1;

omega_1=-g*I_0/2;
omega_2=-M - g*I_0;
omega_3=M - g*I_0;
[w,k]=meshgrid([-10:0.05:10],[-10:0.05:10]);
z= (w+g*I_0/2).^2 - k.^2 - (M-g*I_0*M./2./(w+g*I_0)).^2;
z_1=w+g/2*I_1+k;
z_2=w+g/2*I_1-k;
v = [0,0];
figure(2);
subplot(2,2,3)
contour(k,w,z,v,'b','LineWidth',2);
hold on; contour(k,w,z_1,v,'cyan','LineWidth',2);
hold on; contour(k,w,z_2,v,'cyan','LineWidth',2);
% hold on; plot (0,omega_2,'ko');
% hold on; plot (0,omega_3,'ko');
axis equal
grid on;
xticks([-10:1:10])
yticks([-10:1:10])
ylabel('$\omega$','Interpreter', 'latex','FontSize', 15)
xlabel('$k$','Interpreter', 'latex','FontSize', 15)
title(['I_0 = ',num2str(I_0),';I_1 = ',num2str(I_1),';M = ',num2str(M),';g = ',num2str(g)]);


I_0=1;
I_1=I_0/2;
M=1;
g=1;

omega_1=-g*I_0/2;
omega_2=-M - g*I_0;
omega_3=M - g*I_0;
[w,k]=meshgrid([-10:0.05:10],[-10:0.05:10]);
z= (w+g*I_0/2).^2 - k.^2 - (M-g*I_0*M./2./(w+g*I_0)).^2;
z_1=w+g/2*I_1+k;
z_2=w+g/2*I_1-k;
v = [0,0];
figure(2);
subplot(2,2,4)
contour(k,w,z,v,'b','LineWidth',2);
hold on; contour(k,w,z_1,v,'cyan','LineWidth',2);
hold on; contour(k,w,z_2,v,'cyan','LineWidth',2);
% hold on; plot (0,omega_2,'ko');
% hold on; plot (0,omega_3,'ko');
axis equal
grid on;
xticks([-10:1:10])
yticks([-10:1:10])
ylabel('$\omega$','Interpreter', 'latex','FontSize', 15)
xlabel('$k$','Interpreter', 'latex','FontSize', 15)
title(['I_0 = ',num2str(I_0),';I_1 = ',num2str(I_1),';M = ',num2str(M),';g = ',num2str(g)]);





% 4*w.^4 + w.^3*(12*g*I_0) + w.^2.*(13*g^2*I_0^2 - 4*M^2 -2.*k.^2)+w.*(6*g^3*I_0^3 - 4*M^2*g*I_0 - 4*k^2*g*I_0)+g^4*I_0^4-M^2*g^2*I_0^2-k.^2*2*g^2*I_0^2 =0

% ezplot( 0 == 4*w.^4 + w.^3*(12*g*I_0) + w.^2.*(13*g^2*I_0^2 - 4*M^2 -2.*k.^2)+w.*(6*g^3*I_0^3 - 4*M^2*g*I_0 - 4*k^2*g*I_0)+g^4*I_0^4-M^2*g^2*I_0^2-k.^2*2*g^2*I_0^2,[-0.3 0.3 1 7])

% figure;
% f = @(w,k) 4*w.^4 + w.^3*(12*g*I_0) + w.^2.*(13*g^2*I_0^2 - 4*M^2 -2.*k.^2)+w.*(6*g^3*I_0^3 - 4*M^2*g*I_0 - 4*k^2*g*I_0)+g^4*I_0^4-M^2*g^2*I_0^2-k.^2*2*g^2*I_0^2;
% fimplicit(f,[-0.3 0.3 2.6 5.3])


% g=100;
% M_0=2;
% ky=0;
% omega_5= - g/2 + (sqrt(k.^2+ky^2))./(sqrt(1 - (4*M_0^2)/(g^2)));
% omega_6= - g/2 - (sqrt(k.^2+ky^2))./(sqrt(1 - (4*M_0^2)/(g^2)));
% 
% figure; 
% plot(k,omega_5,'cyan','LineWidth',2);
% hold on; plot(k,omega_6,'cyan','LineWidth',2);
% ylabel('$\omega$','Interpreter', 'latex','FontSize', 15)
% xlabel('$k$','Interpreter', 'latex','FontSize', 15)
% grid on;




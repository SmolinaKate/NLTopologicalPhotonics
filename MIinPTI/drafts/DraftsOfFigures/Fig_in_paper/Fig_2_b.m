clear all;
close all;

step_G=0.001;
G=-4:step_G:4;
J_2=1/sqrt(2);
J_1=1;
I_0 = 1;


inc(1:length(G))=0;
k(1:length(G))=0;

for q=1:1:length(G)
  if (-(G(q)*I_0+J_1^2/J_2)>0)
     inc(q) = abs(G(q)*I_0+J_1^2/J_2);
  end
end
for q=1:1:length(G)
   if (-(G(q)*I_0+J_1^2/J_2)>0)
     k(q)=-(sqrt(abs(G(q).*I_0*J_2+J_1^2)/J_2^2));
   end 
end

p = panel();
p.pack(2, 1);
p(1, 1).select();
hold on; plot(G,inc,'r')
hold on; plot([-sqrt(2) -sqrt(2)], [0 2],'--b')
xlabel('$\Gamma$')
ylabel('$\Im(\lambda)$')
box on; grid on;

p(2, 1).select();
hold on; plot(G,k,'r')
hold on; plot([-sqrt(2) -sqrt(2)], [-2 0],'--b')
xlabel('$\Gamma$')
ylabel('$\kappa_y^{max}$')
box on; grid on;


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
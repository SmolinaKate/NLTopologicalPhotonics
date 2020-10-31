
clear all
addpath('C:\Users\smoli\Dropbox\Topological Nonlinear 2D\for pictures for poster')

path='C:\Users\smoli\Dropbox\Topological Nonlinear 2D\modulation instability\data';



st=0.03;
A=0.9600;
for A=0.9600:st:A+st*5
    
            figure('Renderer', 'painters', 'Position', [20 20 4*250 250])
            p = panel();
            border = 20; 

            p.pack(1, 3);


            for ind=1:1:3
                p(1, ind).marginright = border; 
                p(1, ind).marginleft = border; 
            end


        M=0.1;
        g=M/3;
        I1=A^2;
        omega=-g*I1/4;
        k_x=-g*I1/4;

        kappa_x=[-0.6:0.01:0.6];
        kappa_y =0;
        k_y=0;
        psi_2= - A;
        psi_1=A;

        syms Omega


        for ix = 1:1:length(kappa_x)
            Mat=[-Omega-omega+M - 2*g*(psi_1)^2, - g*psi_1^2,k_x+kappa_x(ix),0;
            k_x+kappa_x(ix),0,-Omega - omega  - M - 2*g*psi_2^2, - g*psi_2^2; 
            - g*psi_1^2,Omega - omega + M - 2*g*psi_1^2,k_x-kappa_x(ix),0; 
            0,k_x-kappa_x(ix),- g*psi_2^2,Omega - omega - M - 2*g*psi_2^2]; 

            eqn = det(Mat) == 0;
            freq(ix,:) = solve(eqn,Omega);
            frq1(ix,:) = double(freq(ix,:));
            frq2(ix,:) = sort(frq1(ix,:));
            frq(ix,:) =frq2(ix,:);
        end


        p(1, 1).select();
        box on; grid on; 
        for i = 1:1:4
            hold on; plot(kappa_x,imag(frq(:,i)))
        end

        ylabel('$Im(\Omega_n)$')
        xlabel('$\kappa_x$')

      
        
        
       file_name = sprintf('values_der_A_#%d',A);
       ToSaveFileName = [file_name '.mat'];
       ToSaveFileName1 = fullfile(path, ToSaveFileName);
       load(ToSaveFileName1);
        
        f_signal= fft(abs( psiA(Nem,:)).^2);
        
        
        p(1, 2).select(); box on; grid on;
        plot( x, abs( psiA(Nem,:)).^2);
        xlabel('$x$')
        ylabel('$|\psi_1(y=0)|^2$')
        
        p(1, 3).select(); box on; grid on;
        %plot( k_x(ceil(length(k_x)/2):ceil(length(k_x)/2)+100), abs(f_signal(ceil(length(k_x)/2):ceil(length(k_x)/2)+100))./max(abs(f_signal))); 
        
        plot( k_x(1:70), abs(f_signal(1:70))./max(abs(f_signal))); 
        xlabel('$k_x$')
        ylabel('$|f(|\psi_1(y=0)|^2)_\omega|$')
        
        
        a = findobj(gcf); 
        size = 11;

        allaxes  = findall(a,'Type','axes');
        alllines = findall(a,'Type','line');
        alltext  = findall(a,'Type','text');

        set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
                    'FontWeight', 'normal', 'FontAngle', 'normal');
        set(alllines,'LineWidth',1);
        set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
                    'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');


end        


A=[0.9800:st:A+st*5];
k_num=[0.31 0.2551 0.2763 0.318 0.3188 0.3188];
k_ins=[0.25 0.27 0.28 0.3 0.31 0.33];

figure; box on; grid on; 
plot(A,k_ins,'o');
hold on; plot(A,k_num,'+');
legend('(a)','(c)')
xlabel('A')
ylabel('k')
 a = findobj(gcf); 
        size = 11;

        allaxes  = findall(a,'Type','axes');
        alllines = findall(a,'Type','line');
        alltext  = findall(a,'Type','text');

        set(allaxes,'FontSize', size, 'FontName', 'Times New Roman', ...
                    'FontWeight', 'normal', 'FontAngle', 'normal');
        set(alllines,'LineWidth',1);
        set(alltext,'FontSize', size, 'FontName', 'Times New Roman', ...
                    'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');


        

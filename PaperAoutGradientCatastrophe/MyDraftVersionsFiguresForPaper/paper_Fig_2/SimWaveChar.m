function [Ip,treshhold,x_sim_wave] = SimWaveChar(dt,Num_Extra_Tsteps,Num_Save_Tsteps,l_x,A,B,g,M)

t=0:dt:Num_Extra_Tsteps*Num_Save_Tsteps*dt;
dx=l_x;
N_x_sim_wave  = 2^13;
L_x = l_x * N_x_sim_wave;%%для построения нужна больщая сетка по икс, чтобы характеристики не вышли
x_sim_wave = zeros( 1, N_x_sim_wave );
for ix = 1 : N_x_sim_wave
    x_sim_wave ( 1, ix ) = ( ix - ( 1 + N_x_sim_wave  / 2 ) ) * l_x;
end

%%построение характеристик для простой волны 
treshhold=10^(-10);
I_0 = (A*exp ( - B*(x_sim_wave-0).^2) + treshhold);
v_g = 1 - g^2*I_0.^2/4/M^2;
for i=1:1:length(v_g)
xP(i,:) = x_sim_wave(i) + v_g (i).*t;%%характеритистики
end

% t=2*M^2*sqrt(2.71828)/g^2./sqrt(B)./A^2;

Ip = zeros(length(I_0), length(t));
for it=2:1:length(t)
    for j=1:1:length(I_0)
         if round(xP(j,it)/dx) > 0 && round((xP(j,it)-x_sim_wave(1))/dx) <= length(I_0)
            Ip(round((xP(j,it)-x_sim_wave(1))/dx),it)=I_0(j);%%амплитуда сохраняется на этих характеристиках
         end
    end
end


end
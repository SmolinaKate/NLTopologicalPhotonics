function s_new = SimWaveProf(Ip,Num_Extra_Tsteps,treshhold,it,x_sim_wave,Nx)
        x_approx = [];
        r=1;
        for j=1:1:length(Ip(:,1))
                if (abs(Ip(j,it*Num_Extra_Tsteps)))>=treshhold
                   x_approx(r)=x_sim_wave(j);
                   Prof_br(r)=Ip(j,it*Num_Extra_Tsteps);
                   r=r+1;
                end
        end         
        
        vq = interp1(x_approx,Prof_br,x_sim_wave);
        [a,b]=max(vq);
        s=vq(b-floor(Nx/2)+1:b+floor(Nx/2));
        [M,N] = max(abs(s));
        [Ny,Nx] = size(s);
        s_new(:)= circshift(s(:), - N + ceil(Nx/2));
end
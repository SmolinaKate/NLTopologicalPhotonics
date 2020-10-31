function [ psiA, psiB ] = YSolitonDistribution (y, y0, M, k, g, I1)
        omega = - k - g*I1/2;
        Omega = sqrt (M^2 + k^2);
        delta = atan (k/M);
    if y>=y0
            y00 = 1/sqrt(Omega^2 - omega^2)*atanh(sqrt(Omega^2-omega^2)/(Omega-omega)*tan(pi/4+delta/2));
            alphan = atan ((Omega - omega)./sqrt(Omega^2 - omega^2).*tanh (sqrt(Omega^2 - omega^2).*(y-y0+y00)));
            alphas = alphan - delta/2;
            rho = 2*(Omega.*cos(2.*alphan)-omega)./g./(1 + cos(2*alphas).^2);
            psiA_s = real(sqrt(2.*rho).*cos(alphas));
            psiB_s =  - real(sqrt(2.*rho).*sin(alphas));
            psiA = psiA_s;
            psiB = psiB_s;
    else
            y00 = 1/sqrt(Omega^2 - omega^2)*atanh (sqrt(Omega+omega)/sqrt(Omega-omega)*tan(pi/4+delta/2));
            alphan = atan ((Omega - omega)./sqrt(Omega^2 - omega^2).*tanh (sqrt(Omega^2 - omega^2).*(-y+y0+y00)));
            alphas = alphan - delta/2;
            rho = 2*(Omega.*cos(2.*alphan)-omega)./g./(1 + cos(2*alphas).^2);
            psiA_s = real(sqrt(2.*rho).*cos(alphas));
            psiB_s =  - real(sqrt(2.*rho).*sin(alphas));
            psiA = - psiB_s;
            psiB = - psiA_s;
    end

end


function [ psiA, psiB ] = NonlinearStep( dt, psiA, psiB, epsilon,g )
%-------------------------------------------------------------------------%
psiA = exp( ( 1i * dt ).* ( g*abs( psiA ).^2 - epsilon ) ).* psiA;
psiB = exp( ( 1i * dt ).* ( g*abs( psiB ).^2 + epsilon ) ).* psiB;
%-------------------------------------------------------------------------%
end

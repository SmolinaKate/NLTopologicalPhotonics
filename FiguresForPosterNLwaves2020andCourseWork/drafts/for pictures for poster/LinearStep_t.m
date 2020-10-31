function [ psiA, psiB ] = LinearStep_t( dt, psiA, psiB, k_x, k_y, link, a)
%-------------------------------------------------------------------------%
N_x = length( k_x );
N_y = length( k_y );
psiA = fft2( psiA );
psiB = fft2( psiB );
%-------------------------------------------------------------------------%
for iy = 1 : N_y
    for ix = 1 : N_x
        if( k_x( 1, ix ) == 0 && k_y( 1, iy ) == 0 )
            psiA( iy, ix ) = psiA( iy, ix );
            psiB( iy, ix ) = psiB( iy, ix );
        else
            h=(sqrt(3)*link*a/2*(k_x( 1, ix) -1i*k_y( 1, iy)));
            lambda=sqrt(h*conj(h));
            a = 1;
            b = 1;
            c = ( conj(h) ) / ( + lambda );
            d = ( conj(h) ) / ( - lambda );
            detA = a * d - b * c;
            a_inv = d / detA;
            b_inv = - b / detA;
            c_inv = - c / detA;
            d_inv = a / detA;
            s1 = a_inv * psiA( iy, ix ) + b_inv * psiB( iy, ix );
            s2 = c_inv * psiA( iy, ix ) + d_inv * psiB( iy, ix );
            s1 = exp( ( + 1i * dt ) * ( - lambda ) ) * s1;
            s2 = exp( ( + 1i * dt ) * ( + lambda ) ) * s2;
            psiA( iy, ix ) = a * s1 + b * s2;
            psiB( iy, ix ) = c * s1 + d * s2;
        end
    end
end
%-------------------------------------------------------------------------%
psiA = ifft2( psiA );
psiB = ifft2( psiB );
%-------------------------------------------------------------------------%
end

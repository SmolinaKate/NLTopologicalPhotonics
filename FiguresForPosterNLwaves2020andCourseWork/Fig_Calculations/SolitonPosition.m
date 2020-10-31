function [ xMax, yMax ] = SolitonPosition( psiA, psiB, x, y )
%-------------------------------------------------------------------------%
Intensity = abs( psiA ).^2 + abs( psiB ).^2;
[ maxValue, maxIndex ] = max( Intensity( : ) );
[ iyMax, ixMax ] = ind2sub( size( Intensity ), maxIndex );
xMax = x( 1, ixMax );
yMax = y( 1, iyMax );
%-------------------------------------------------------------------------%
end

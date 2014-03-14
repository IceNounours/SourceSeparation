function [ z ] = Whiten( x )
%WHITEN Whiten the signal x (i.e othronormalize the signals)
%   Perform a Gram-Schmidt orthonormalisation
    numLines = size(x,1);
    numRows = size(x,2);
    fNumRows = single(numRows);
    
    z = x;

    C = z*z' ./ fNumRows;
    
    % orthogonalize
    for i=1:numLines
        for j=1:(i-1)
            z(i,:) = z(i,:) - (C(i,j)/C(j,j)) * z(j,:);
        end
    end
    
    % scale to unit variance
    z = sqrt(fNumRows) .* normr(z);
end


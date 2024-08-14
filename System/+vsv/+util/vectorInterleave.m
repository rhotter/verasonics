function vInterleave = vectorInterleave(v1,v2)
%VECTORINTERLEAVE Interleaves vector 1 and vector 2 
%
%   vInterleave(i) = v1(round( i/2 )); % if i is odd
%   vInterleave(i) = v1( i/2 ); % if i is even
%
%   vInterleave is a column vector
% Version 1.0 | 2020-05-01 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 


% a = {'Pram1', 'Param2', 'Pram3'};
% b = {1 2 3};
% col_interleave = reshape([a(:) b(:)]',2*size(a,1), [])' 
% col_interleave =
% 
%   3×2 cell array
% 
%     {'Pram1' }    {[1]}
%     {'Param2'}    {[2]}
%     {'Pram3' }    {[3]}
%     

% interleaveM  = reshape([v1(:) v2(:)]',2*size(v1,1), []);
% vInterleave  = interleaveM(:);

vInterleave(1:2:2*numel(v1)) = v1;
vInterleave(2:2:2*numel(v1)) = v2;

end


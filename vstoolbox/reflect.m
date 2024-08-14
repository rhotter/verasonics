%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = reflect(x);

flipped = size(x,1)==1 ;
   if flipped
       x = x(:);
   end

dx = x - x(1);

y = x(1) - dx(end:-1:2,1) ;

if flipped
    y=y(:).';
end

end

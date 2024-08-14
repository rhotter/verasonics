function w=window(fhand,varargin);


w = feval(fhand,varargin{:});
end

function mustBeScalarOrEmpty( obj )
%MUSTBESCALAROREMPTY Summary of this function goes here
%   Detailed explanation goes here

    if ~vsv.util.isScalarOrEmpty(obj)
        error('vsv:util:mustBeScalarOrEmpty', ....
            'Object must be a scalar or empty');
    end

end


classdef (Hidden) AbstractInputSetting < vsv.extio.mixin.ParamImportable
%AbstractInputSetting basic class defining extIOSettings
%  
%
% Version 1.0 | 2020-06-16 
% $Author: Dr. Jack Potter, Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
    properties
        
        % buffer settings are part of all ExtIO settings 
        % @type vsv.extio.BufferSetting
        BufferSetting vsv.extio.BufferSetting ...
                        {vsv.util.mustBeScalarOrEmpty(BufferSetting)} ...
                            = vsv.extio.BufferSetting.empty();
    end
    
    methods
        
        function this = AbstractInputSetting(varargin)
            
            this@vsv.extio.mixin.ParamImportable(varargin{:});
            if isempty(this.BufferSetting)
                this.BufferSetting = vsv.extio.BufferSetting();
            end
            
        end
    end
    
end


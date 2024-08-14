function [ fInd ] = cellfind( ce, ex, exact, usei )
%CELLFIND uses regular expression to find a string in a cellstring array
%  
%  If exact is false, the function will return a logical array with true
%  everywhere were the regular expression would find at least one match. 
% 
%  The function can handle different input types for the cell array ce. It
%  can be mixed of any type with the following rules of conversions:
%   - numeric values are converted to chars using num2string
%   - strings ("A String") will be converted to chars
%   - All other types will be converted to ''
%
%  NOTE: This function is undocumented and may change in future releases
%
% Usage:
%   ce = {'test.mat', 'ha.llo', 'file.mat' 'Test.mat'};
%   ex = '[a-z][\.]mat';
%   [ fInd ] = vsv.util.cellfind( ce, ex );
%
%   usei = true;
%   [ fInd ] = vsv.util.cellfind( ce, ex , [], true);
%
%   ex = 'test';
%   exact = true;
%   [ fInd ] = vsv.util.cellfind( ce, ex , exact);
%
% Parameters:
%   @param ce - @type string cell array, that contains the names to search 
%               for, if is char ce will be converted to a cell
%   @param ex - @type char the expression to search for
%   @param exact- @type logical, true if match needs to be exact, false
%                 other, @optional, @default false
%   @param usei - @type logical, true if regexpi should be used, false if regexp
%                 should be used, (@optional), default = false
%
% Return values:
%   @return fInd - @type logical array, fInd(i) is true if ce{i} matches expression
%                  ex
%
% $Author: Dr. Daniel Rohrbach,
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%

invalidSearchStringErrorID  = 'cellfind:invalidSearchString';
invaldiSearchStringErrorMsg = 'ERROR cannot parse cell. Cell must be a cell string';

if ischar(ex) || ( isstring(ex) && length(ex) == 1 )
    if nargin > 2 && ~isempty(exact) && exact
        if ischar(ce) || iscellstr(ce) || isstring(ce)
            fInd  = strcmp( ex,  ce );
        else
            if iscell(ce) 
                ce   = makeCellCompatible(ce);
                fInd = vsv.util.cellfind(ce, ex, true);
            else
                error( invalidSearchStringErrorID, ...
                       invaldiSearchStringErrorMsg);
            end
        end
    else
        if ischar(ce) || iscellstr(ce) || isstring(ce)

            % this would result in inconsistent results with exact which
            % returns true for empy matches 
            if ( isempty(ex) && ischar(ex) ) || ( isstring(ex) && isempty( convertStringsToChars(ex) )) 
                exact = true;
                fInd  = vsv.util.cellfind( ce, ex, exact);
            else
                if ischar(ce)
                    ce = {ce};
                elseif isstring(ce) && length(ce) == 1
                    ce = { convertStringsToChars(ce) };
                end

                if nargin > 3 && usei
                    fInd  = cellfun(@(x)~isempty(x), regexpi( ce, ex, 'once') );
                else
                    fInd  = cellfun(@(x)~isempty(x), regexp( ce, ex, 'once' ) );
                end
            end
        else
            % last try check whether we can convert it to a cell
            if iscell(ce)
                ce = makeCellCompatible(ce);
                if nargin > 3 && usei
                    [ fInd ] = vsv.util.cellfind( ce, ex, [], usei );
                else
                    [ fInd ] = vsv.util.cellfind( ce, ex );
                end

            else
                error( invalidSearchStringErrorID, ...
                       invaldiSearchStringErrorMsg);
            end

        end
    end
else
    error('cellfind:invalidPattern', ...
          [ 'ERROR cannot parse pattern. Pattern must be a char row vector,' ...
            'a cell array of char row vectors, or a string array.']);
end

end

function ce = makeCellCompatible(ce)
    % convert numebers to chars
    ce = cellfun( @(x)convertNum2String(x), ce, 'UniformOutput', false );
    % convert non chars to empty
    ind     = cellfun( @(x)~ischar(x), ce );
    ce(ind) = {''};
end

function output = convertNum2String(input)
    
    output = '';
    
    if ischar(input)
        output = input;
    elseif isnumeric(input)
        output = num2str(input);
    elseif isstring(input)
        output = convertStringsToChars(input);
    end
    
    
end


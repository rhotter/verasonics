function [ C ] = text2cell(varargin)
%text2cell(<filename>,delimiter)
%  Converts the lines of text between the first two instances of the delimiter
%  in the file specified. (If filename is missing, the calling file is used.)
%  A line which contains the delimiter in quotes is not recognized as an instance
%  of the delimiter.

% March 2020 VTS-1480 Remove '_KeyWord' suffix, for 2019b compatibility
% (for callbacks or external functions embedded in user scripts, with
% persistent variable declaration)

switch nargin
    case 0
        error('text2cell: Inputs - (<filename>,delimiter)  missing filename will use calling file.');
    case 1
        % Assume only delimiter specified.
        % Get the name of the file that called this function.
        ST = dbstack(1,'-completenames');
        if isempty(ST)
           error('text2cell: no caller file. Please run entire script or assign filename.');
        end
        filename = ST.file;
        if ~ischar(varargin{1})
            error('text2cell: delimiter must be a string in quotes or a character array.');
        end
        delimiter = varargin{1};
    case 2
        if ~ischar(varargin{1})
            error('text2cell: filename must be a string in quotes or a character array.');
        end
        filename = varargin{1};
        if ~ischar(varargin{2})
            error('text2cell: delimiter must be a string in quotes or a character array.');
        end
        delimiter = varargin{2};
    otherwise
        error('text2cell: too many inputs. Inputs - (<filename>,delimiter)');
end

fid = fopen(filename);
if fid<0, error('text2cell: Couldn''t open %s for reading.',filename); end

n = 0;
tline = fgetl(fid);
while ischar(tline)
   m = strfind(tline, delimiter);
   if isempty(m)
       tline = fgetl(fid);
       continue
   else
       m = strfind(tline, ['''' delimiter '''']);
       if ~isempty(m)  % skip any line with the delimiter in quotes
           tline = fgetl(fid);
           continue
       end
   end
   n = n + 1;
   tline = fgetl(fid);
   j = 1;
   while ischar(tline)
       m = strfind(tline, delimiter);
       if ~isempty(m)
           fclose(fid);
           return
       end
       tline = erase(tline, '_KeyWord'); % VTS-1480 eliminate characters added to bypass avoid errors
       C{j} = tline;
       j = j + 1;
       tline = fgetl(fid);
   end
end
if n==0, error('text2cell: Delimiter ''%s'' not found in %s.', delimiter, filename); end
fclose(fid);

end

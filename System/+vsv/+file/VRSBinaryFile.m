classdef VRSBinaryFile < handle
    %VSXBINFILE is the main class representing the VSX binary file format
    %
    % Usage:
    %
    %   % create an object for reading vsx binary files
    %   bin = vsv.file.VRSBinaryFile;
    %   % use readFile function for reading the binary file
    %   bin.readFile;
    %   data = bin.formatData;
    %
    % Examples
    %   see example.file.vsxBinFileExample for more details
    %
    % $Author: Dr. Daniel Rohrbach,
    % Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
    % and Trademark Office.
    %
    properties(SetAccess = private)

        % @type vsv.file.VRSBinaryHeader the file header
        header;

        % @type numeric the data of the file
        data;

    end

    properties(Constant)
        % java class to help encoding the data if compressed
        coder = vsv.file.DictCoder(16);
    end

    properties(Dependent=true)

        % @type logical indicating whether a time tag was enabled
        isTimeTag;

        % @type char a user defined comment
        comment;

        % @type char a user defined study ID
        studyID;

        % @type char a user defined sample ID
        sampleID;

        % @type numeric time tag month
     	filePath;
    end

    properties(Dependent=true, SetAccess=private)

        % @type logical indicating whether a time tag was enabled
        version;

        % @type numeric time tag second
     	dim;

        % @type numeric time tag hour
     	numDataPoints;

        % @type numeric time tag day
     	dataType;

    end


    methods

        function data = formatData(this)
        % Returns the data as a 3D volume object with three dimensions
        %
        %   The data is storred as a linear byte stream and need to be
        %   formatted in order to represent the original data set
        %   structure.
        %
        % Usage:
        %   data = bin.formatData();
        % Parameters:
        %
        % Return values:
        %   data - the formatted data
        %

            dimension  = this.header.dim;

            % tolerate one mismatch
            if( abs( prod( dimension ) - numel(this.data) ) == 0 )
                data = reshape(this.data, dimension(:)' );
            else
                N = prod(dimension);
                % try to fix
                data = zeros( N, 1, class(this.data) );
                if N < numel( this.data )
                    data(:) = this.data(1:N);
                else
                    data(1:numel(this.data)) = this.data(:);
                end

                data = reshape( data, dimension(:)' );
            end

        end

        function clear(this)
        % clears the header and data of this file object.
        %   this is called when reading the file

            this.header.clearHeader();
            this.data  = [];
        end


        function this = VRSBinaryFile(file)
        % Constructor for creating the vsx binary file format
        %
        % Usage:
        %   bin = vsv.file.VRSBinaryFile(); - will construct an empty storage
        %   object that can be used for loading vrs data
        %   bin = vsv.file.VRSBinaryFile(file); - will try to read given file
        %
        % Parameters:
        %   @param file - @type char a file name of the file to load, (optional)
        % Return values:
        %   the object

            this.header = vsv.file.VRSBinaryHeader;
            if nargin >= 1
                isOK = this.readFile(file);
                if ~isOK
                    error('ERROR: Unable to read file. Cannot create object');
                end
            end

        end

        function isOK = saveAsMatFile(this, file, saveRaw, objname)
        % Saves the binary file object as a mat file
        %
        %   This will create a mat file and saves the object to the file.
        %
        % Usage:
        %   bin.saveAsMatFile(); - opens a dialog for selecting the file
        %   for storage
        %   bin.saveAsMatFile(file); - will try to save the data to file
        %   bin.saveAsMatFile(file, saveRaw) - if saveRaw is true the data
        %   will be converted to matlab numeric matrix
        %   bin.saveAsMatFile(file, [], objname) - specify the object name
        %   inside the mat file
        %   bin.saveAsMatFile(file, saveRaw, objname)
        %   bin.saveAsMatFile([], saveRaw, objname)
        %   bin.saveAsMatFile([], [], objname)
        %
        % Parameters:
        %   @param file - @type char the file name (optional)
        %   @param saveRaw - @type logical indicating the whether to save raw data
        %   @param objname - @type char the name of the object inside the function
        % Return values:
        %   @return isOK - @type logical true if saving was successful, false
        %   otherwise
        %

            if nargin < 3 || vsv.util.isNULL(saveRaw)
                saveRaw = false;
            end
            if nargin < 4 || vsv.util.isNULL(objname) || isempty(objname)
                objname = 'vrsData';
            end

            if nargin < 2 || isempty(file)
                % the file path based from the header
                fip = this.header.filePath;
                if vsv.util.isNULL(fip)
                    file = '';
                    fn   = 'VRS_binaryFile';
                else
                    % remove extension
                    file = fip;
                    [fp, fn] = fileparts(file);
                    file = fullfile( fp, [fn '.mat'] );
                end
                file = vsv.file.FileTools.obtainFilenamePut(file, '.mat', 'Get Mat file to save vrs file.', fn);
            end

            %
            if ~isempty(file)

                if( saveRaw )
                    data = this.formatData(); %#ok
                    info = this.header.toStruct();%#ok
                    eval( [ objname ' = data;'] );
                    eval( [ objname 'Info = info;'] );
                    save(file, objname, [objname 'Info'] );
                    isOK =true;
                else
                    eval( [ objname ' = this;'] );
                    save(file, objname);
                    isOK =true;
                end
            else
                isOK = false;
            end


        end

        function [isOK] = convertToMatFile(this, file, newMatFile, varargin)
        % reads a vsx file and converts it to a mat file
        %
        %   This will first read in the binary file specified by file and
        %   then saves the file as a mat file
        %
        %   this function uses readFile and saveAsMatFile, see these
        %   functions for documentation.
        %
        %   If no newMatFile is given the function will create a file name
        %   for the mat file based on the vrs file name but with .mat as
        %   the file tag.
        %
        %   The user can also specify two parameters (i.e., saveRaw,
        %   objname) that allow to specify the type of data that is saved
        %   to the mat file
        %
        %   Usage:
        %     [isOK] = bin.convertToMatFile(file) - reads file and stores
        %     as mat file with same file name but with .mat filetag
        %     [isOK] = bin.convertToMatFile(file, newMatFile) - reads file
        %     and stores mat file in newMatFile
        %     [isOK] = bin.convertToMatFile(file, newMatFile, saveRaw, objname) -
        %     see saveAsMatFile for information
        %   Parameters:
        %     this -
        %     @param file - the file name (optional)
        %     @param newMatFile - @type char a new mat file
        %   Return values:
        %     @return isOK @type logical true if successful, false otherwise

            if nargin < 1
                file = '';
            end

            isOK = this.readFile(file);

            if(isOK)
                if nargin > 2
                    this.saveAsMatFile(newMatFile, varargin{:});
                else
                    this.saveAsMatFile([], varargin{:});
                end
            end

        end

        function [isOK, msg] = readFile(this, file)
        % Function to read a vsx file
        %  this reads the header first then the data
        %
        % Usage:
        %   bin.readFile( "file.vrs" ); - read file.vrs
        %   bin.readFile( ); - opens a dialog for reading a file
        % Parameters:
        %   this the object reference
        %   @param file (optional) @type char a string specifying a file to load
        % Return values:
        %   @return isOK @type logical which is TRUE if file operation was
        %   successful false otherwise
        %   @return msg @type char a message string indicating the error, or empty
        %   if isOK is TRUE
        %

            isOK = false;
            msg  = '';

            if nargin < 2
                file = [];
            end

            this.clear();
            fid = this.header.readHeaeder( file );

            % something went wrong or the user aborted
            if isempty(fid) || fid < 0
                isOK = false;
                msg  = 'File not found!';
            else
                if ~isempty( this.header.numDataPoints ) && ~isempty( this.header.dataType)
                    this.data = fread(fid, this.header.numDataPoints, this.header.getFormatString());
                    this.data = cast(this.data, this.header.getFormatString());
                    fclose(fid);
                    if this.header.compression
                        if this.header.isSignedCompress
                            this.data = this.coder.decodeSigned( this.data );
                            this.data = this.data(1:end-1);
                        else
                            this.data = this.coder.decode( this.data );
                            this.data = this.data(1:end-1);
                            if ischar(this.data)
                                this.data = int16(this.data);
                            end
                        end
                        %this.data(end) = [];
                    end
                    isOK = true;
                end
            end

        end

        function setData( this, data )
        % allows to set some user defined data for storage and import into
        % VSX
        %
        % The function will modify the header information to fit the given
        % data and prepare the object for file read access.
        %
        % Parameters :
        %   @param data - the data which is to be added to this object
        %

            info = whos('data');
            switch info.class
                case 'double'
                    this.header.dataType = vsv.common.VCDataType.Vc_DOUBLE;
                case 'single'
                    this.header.dataType = vsv.common.VCDataType.Vc_FLOAT;
                case 'int8'
                    this.header.dataType = vsv.common.VCDataType.Vc_S8;
                case 'int16'
                    this.header.dataType = vsv.common.VCDataType.Vc_S16;
                case 'int32'
                    this.header.dataType = vsv.common.VCDataType.Vc_S32;
                case 'int64'
                    this.header.dataType = vsv.common.VCDataType.Vc_S64;
                case 'uint8'
                    this.header.dataType = vsv.common.VCDataType.Vc_U8;
                case 'uint16'
                    this.header.dataType = vsv.common.VCDataType.Vc_U16;
                case 'uint32'
                    this.header.dataType = vsv.common.VCDataType.Vc_U32;
                case 'uint64'
                    this.header.dataType = vsv.common.VCDataType.Vc_U64;
                otherwise
                    error('ERROR: data Type is not supported');
            end

            % get the size of the data
            si = size(data);
            this.header.dim = ones( 1, 4);
            leng = min( [ length(si), 4 ]);
            for i = 1 : leng
                this.header.dim(i) = si(i);
            end

            % compression must be zero this function does not compress data
            this.header.compression = 0;
            this.header.isSignedCompress = 0;
            % assign the number of data Points
            this.header.numDataPoints = numel(data);
            this.data = data;

        end

        function [isOK, msg] = writeFile(this, file)
        % Function to read a vsx file
        %  this reads the header first then the data
        %
        % Usage:
        %   bin.readFile( "file.vrs" ); - read file.vrs
        %   bin.readFile( ); - opens a dialog for reading a file
        % Parameters:
        %   @param this the object reference
        %   @param file (optional) @type char a string specifying a file to load
        % Return values:
        %   @return isOK @type logical which is TRUE if file operation was
        %   successful false otherwise
        %   @return msg @type char a message string indicating the error, or empty
        %   if isOK is TRUE
        %

            isOK = false;
            msg  = '';

            if nargin < 2
                file = [];
            end

            dd = this.formatData;
            switch class(this.data)
                case 'double'
                    dd = double( dd );
                case 'single'
                    dd = single( dd );
                case 'int8'
                    dd = int8( dd );
                case 'int16'
                    dd = int16( dd );
                case 'int32'
                    dd = int32( dd );
                case 'int64'
                    dd = int64( dd );
                case 'uint8'
                    dd = uint8( dd );
                case 'uint16'
                    dd = uint16( dd );
                case 'uint32'
                    dd = uint32( dd );
                case 'uint64'
                    dd = uint64( dd );
                otherwise
                    error('ERROR: data Type is not supported');
            end

            % make sure the data matches for writing
            this.setData(dd);

            fid = this.header.writeHeaeder( file );

            % something went wrong or the user aborted
            if isempty(fid) || fid < 0
                isOK = false;
                msg  = 'Cannot write header!';
            else
                try
                    fwrite(fid, this.data(:), class(this.data) );
                catch
                    try
                        fclose(fid);
                    catch
                    end
                    error(lasterror);%#ok
                end
                fclose(fid);
            end

        end

    end


    %% getter and setter
    methods

        function headerStruct = getHeaderInfo(this)
        % converts the header to a struct and returns the information
        %
        % Usage:
        %   info = bin.getHeaderInfo( );
        %
        % Return values:
        %   @return headerStruct   - @type struct the header information
        %

            headerStruct = this.header.toStruct();
        end

        function [sec, min, hour, day, month, year] = getTimeTag(this)
        % get the time tag information
        %
        % Usage:
        %   [sec, min, hour, day, month, year] = bin.setTimeTag( );
        %
        % Return values:
        %   @return sec   - @type numeric the seconds
        %   @return min   - @type numeric the minutes
        %   @return hour  - @type numeric the hours
        %   @return day   - @type numeric the day
        %   @return month - @type numeric the month
        %   @return year  - @type numeric the year
        %
            sec   = this.header.sec;
            min   = this.header.min;
            hour  = this.header.hour;
            day   = this.header.day;
            month = this.header.month;
            year  = this.header.year;
        end

        function setTimeTag(this, sec, min, hour, day, month, year)
        % set the time tag information
        %
        % Usage:
        %   bin.setTimeTag( sec, min, hour, day, month, year);
        %
        % Parameters:
        %   @param sec   - @type numeric the seconds
        %   @param min   - @type numeric the minutes
        %   @param hour  - @type numeric the hours
        %   @param day   - @type numeric the day
        %   @param month - @type numeric the month
        %   @param year  - @type numeric the year
        %

            vsv.util.Errors.checkErr('sec', sec, 'numeric');
            vsv.util.Errors.checkErr('sec', sec, 'size', 1);

            vsv.util.Errors.checkErr('min', min, 'numeric');
            vsv.util.Errors.checkErr('min', min, 'size', 1);

            vsv.util.Errors.checkErr('hour', hour, 'numeric');
            vsv.util.Errors.checkErr('hour', hour, 'size', 1);

            vsv.util.Errors.checkErr('day', day, 'numeric');
            vsv.util.Errors.checkErr('day', day, 'size', 1);

            vsv.util.Errors.checkErr('month', month, 'numeric');
            vsv.util.Errors.checkErr('month', month, 'size', 1);

            vsv.util.Errors.checkErr('year', year, 'numeric');
            vsv.util.Errors.checkErr('year', year, 'size', 1);

            this.header.sec   = sec;
            this.header.min   = min;
            this.header.hour  = hour;
            this.header.day   = day;
            this.header.month = month;
            this.header.year  = year;

        end

        function set.isTimeTag( this, istime )
            if isnumeric(istime)
                istime = logical(istime);
            end
            vsv.util.Errors.checkErr('isTimeTag', istime, 'logical');
            vsv.util.Errors.checkErr('isTimeTag', istime, 'size', 1);
            this.header.istime = istime;
        end

        function set.comment( this, comment )
            if vsv.util.isNULL(comment)
                comment = '';
            end
            vsv.util.Errors.checkErr('comment', comment, 'char');
            this.header.comment = comment;
        end

        function set.studyID(this, studyID )
            if vsv.util.isNULL(studyID)
                studyID = '';
            end
            vsv.util.Errors.checkErr('studyID', studyID, 'char');
            this.header.studyID = studyID;
        end

        % @type char a user defined sample ID
        function set.sampleID(this, sampleID )
            if vsv.util.isNULL(sampleID)
                sampleID = '';
            end
            vsv.util.Errors.checkErr('sampleID', sampleID, 'char');
            this.header.sampleID = sampleID;
        end

        % @type numeric time tag month
     	function set.filePath(this, filePath )
            if vsv.util.isNULL(filePath)
                filePath = '';
            end
            vsv.util.Errors.checkErr('filePath', filePath, 'char');
            this.header.filePath = filePath;
        end


        function val = get.isTimeTag(this)
            val = this.header.isTimeTag;
        end

        function val = get.comment(this)
            val = this.header.comment;
        end

        function val = get.studyID(this)
            val = this.header.studyID;
        end

        function val = get.sampleID(this)
            val = this.header.sampleID;
        end

     	function val = get.filePath(this)
            val = this.header.filePath;
        end

        function val = get.version(this)
            val = this.header.version;
        end

     	function val = get.dim(this)
            val = this.header.dim;
        end

     	function val = get.numDataPoints(this)
            val = this.header.numDataPoints;
        end

     	function val = get.dataType(this)
            val = this.header.dataType;
        end

    end
end


classdef VRSBinaryHeader < handle
    %VSXBINHEADER Class representing the header of a vsx binary file
    %   The vsx binary file format representaiton
    %
    %   File header specification:
    %
    %   Name Tag:   Bytes:      comment:
    %   version     1xuint16    a version number in case file specification
    %                           cahtes in future releases
    %
    %   compression 1xuint8     0 if no compression 1 BL compression,
    %   isTimeTag   1xtimeTag   0 if no time info was written determines the
    %                           following bytes (0/1), 1 if time tag was
    %                           written
    %   sec         1/0xuint8
    %   min         1/0xuint8
    %   hour        1/0xuint8
    %   day         1/0xuint8
    %   month       1/0xuint8
    %   year        1/0xuint8
    %
    %   lSID        1xuint64    length of the following string studyID
    %   studyID     lSIDxuint64 study ID string
    %   lAID        1xuint64    length of the following string sampleID
    %   sampleID    lAIDxuint64 sample ID string
    %   lC          1xuint64    length of the following string sampleID
    %   comment     lCxuint64   comment string
    %
    %   dim         4xuint64    dimension of the array data
    %   numDataPoints
    %               1xunit64    number of data points
    %   dataType    1xuint8     data type
    %
    % If compression == 0 prod( dim(:)) should be equal numDataPoints. If
    % compression was enables. The dimension refers to the dimension of the
    % uncompressed data
    %
    % $Author: Dr. Daniel Rohrbach,
    % Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
    % and Trademark Office.
    %

    properties

        % @type numeric version of the file header.
        version;

        % @type logical indicating whether a time tag was enabled
        isTimeTag;

        % @type numeric indicating the compression to use
        compression;

        % @type numeric indicating the compression to use
        isSignedCompress;

        % @type numeric time tag second
     	sec;

        % @type numeric time tag minute
     	min;

        % @type numeric time tag hour
     	hour;

        % @type numeric time tag day
     	day;

        % @type numeric time tag month
     	month;

        % @type numeric time tag year
     	year;

        % @type char a user defined comment
        comment;

        % @type char a user defined study ID
        studyID;

        % @type char a user defined sample ID
        sampleID;

        % @type numeric dimension of stored data
        dim;

        % @type number of Datapoints
        numDataPoints;

        % @type number of Datapoints
        dataType;

        % @type char absolute file path
        filePath;

    end



    methods(Static = true, Access=private)

        function str = readStudyString(fid)

            str = '';
            leng = fread(fid, 1, 'uint64');
            if(leng>0)
                str = char( fread(fid, leng, 'char')');
            end

        end

        function writeStudyString(fid, str)

            if isempty(str)
                str = '';
            end

            leng = length(str);
            fwrite(fid, leng, 'uint64');

            if(leng>0)
                fwrite( fid, str, 'char');
            end

        end


    end

    methods(Access=private)

        function setDataType(this, num)
            switch num
                case 0
                    this.dataType = vsv.common.VCDataType.Vc_U8;
                case 1
                    this.dataType = vsv.common.VCDataType.Vc_U16;
                case 2
                    this.dataType = vsv.common.VCDataType.Vc_U32;
                case 3
                    this.dataType = vsv.common.VCDataType.Vc_U64;
                case 4
                    this.dataType = vsv.common.VCDataType.Vc_S8;
                case 5
                    this.dataType = vsv.common.VCDataType.Vc_S16;
                case 6
                    this.dataType = vsv.common.VCDataType.Vc_S32;
                case 7
                    this.dataType = vsv.common.VCDataType.Vc_S64;
                case 8
                    this.dataType = vsv.common.VCDataType.Vc_DOUBLE;
                case 9
                    this.dataType = vsv.common.VCDataType.Vc_FLOAT;
            end

        end

        function ndt = getNumDataType(this)

            if isempty(this.dataType)
                ndt = vsv.common.VCDataType.Vc_S16;
            else
                switch this.dataType
                    case vsv.common.VCDataType.Vc_U8
                        ndt = 0;
                    case vsv.common.VCDataType.Vc_U16
                        ndt = 1;
                    case vsv.common.VCDataType.Vc_U32
                        ndt = 2;
                    case vsv.common.VCDataType.Vc_U64
                        ndt = 3;
                    case vsv.common.VCDataType.Vc_S8
                        ndt = 4;
                    case vsv.common.VCDataType.Vc_S16
                        ndt = 5;
                    case vsv.common.VCDataType.Vc_S32
                        ndt = 6;
                    case vsv.common.VCDataType.Vc_S64
                        ndt = 7;
                    case vsv.common.VCDataType.Vc_DOUBLE
                        ndt = 8;
                    case vsv.common.VCDataType.Vc_FLOAT
                        ndt = 9;
                    otherwise
                        error('ERROR: Unexpected Error, data type is unknown');
                end
            end
        end


    end

    methods(Access=protected)
        function initHeader(this)
        % initialize the header with empty strings for char parameters
        %
        %  Default init value is [] which is numeric. This makes sure that
        %  the char parameters of the class are initialized with the
        %  correct data type


            this.comment  = '';
            this.studyID  = '';
            this.sampleID = '';
            this.filePath = '';
        end

    end

    methods

        function st = toStruct(this)
        % converts this object to a struct
        %
        % Usage:
        %   st = he.toStruct();
        % Return values:
        %   @return st - the struct representation of this object

            st.version     = this.version     ;
            st.isTimeTag   = this.isTimeTag   ;
            st.compression = this.compression ;
            st.sec   = this.sec         ;
            st.min   = this.min         ;
            st.hour  = this.hour        ;
            st.day   = this.day         ;
            st.month = this.month       ;
            st.year  = this.year        ;
            st.dim   = this.dim         ;
            st.numDataPoints = this.numDataPoints ;
            st.dataType      = this.dataType      ;
            st.comment       = this.comment;
            st.studyID       = this.studyID;
            st.sampleID      = this.sampleID;
            st.filePath      = this.filePath;
        end


        function format = getFormatString(this)
        % returns the char representation of the data format specified by
        % dataType
            format = vsv.common.VCDataType.getDataString( this.dataType );
        end

        function clearHeader(this)
        % clears the header parameter (i.e., sets them to [] )
        %
        % Usage:
        %   this.clearHeader()
        % Parameters:
        %   this the object reference
        %

            this.version     = [] ;
            this.isTimeTag   = [] ;
            this.compression = [];
            this.sec         = [];
            this.min         = [];
            this.hour        = [];
            this.day         = [];
            this.month       = [];
            this.year        = [];
            this.dim         = [];
            this.numDataPoints = [];
            this.dataType      = [];
            this.initHeader();
            this.isSignedCompress = [];

        end



        function this = VRSBinaryHeader(file)
            %VSXBINHEADER Construct an instance of this class
            %   Detailed explanation goes here
            this.initHeader();
            if nargin >= 1
                isOK = this.readHeaeder(file);
                if ~isOK
                    error('ERROR: Unable to read file. Cannot create object');
                end
            end

        end

        function fid = readHeaeder(this, file )
        % Reads the header from the file
        %
        % the function will read the header from a vrs binary file
        %
        % if file is a char it will open the file. The function will not
        % close the file if the user specifies a return parameter.
        %
        % Usage:
        %   fid = bin.readHeader();
        %   fid = bin.readHeader(file);
        %   bin.readHeader(...); - this will close the file after reading
        %
        % Parameters:
        %   @param file - @type the file either a char giving the absolute path
        %   or a file pointer (optional)
        % Return values:
        %   @return fid - the file pointer to the file set to the data block ready
        %   to read in the file

            fid = [];
            if nargin < 2 || isempty(file)
                [file ] = vsv.file.FileTools.obtainFilename( [], {'*.vrs', 'VSX Files (*.vrs)'; ...
                                                                  '*.vsx', 'VSX Files (*.vsx)'; ...
                                                                  '*.*',   'All Files (*.*)'},   ...
                                                                  'Pick a vsxHeader file' );
            end

            if isempty(file)
                return;
            end

            if ischar(file)
                fid = fopen(file);
            end

            % put reading in try and catch to make sure that we release the
            % file pointer in case something goes wrong

            if isempty(fid) || fid <=0
                return;
            end
            try
                [fp, ~, ~] = fileparts(file);

                % if filepath is empty file was obtained from current folder
                if isempty(fp)
                    fullfile(cd, file);
                end

                this.filePath = file;

                frewind(fid);
                this.version          = fread(fid, 1, 'uint16');
                this.compression      = fread(fid, 1, 'uint8');
                this.isSignedCompress = fread(fid, 1, 'uint8');

                this.isTimeTag   = fread(fid, 1, 'uint8');
                if( this.isTimeTag > 0)
                    this.sec    = fread(fid, 1, 'uint8');
                    this.min    = fread(fid, 1, 'uint8');
                    this.hour   = fread(fid, 1, 'uint8');
                    this.day    = fread(fid, 1, 'uint8');
                    this.month  = fread(fid, 1, 'uint8');
                    this.year   = fread(fid, 1, 'uint8') + 1900;
                end
                this.studyID  = this.readStudyString(fid);
                this.sampleID = this.readStudyString(fid);
                this.comment  = this.readStudyString(fid);


                this.dim           = fread(fid, 4, 'uint64');
                this.numDataPoints = fread(fid, 1, 'uint64');
                this.setDataType( fread(fid, 1, 'uint8' ) );
            catch
                % something went wrong make sure the file gets closed
                try
                    fclose(fid);
                catch
                end
                fid = -1;
                % rethrow the error
                error(lasterror);%#ok thats OK
            end

            % make sure to close the file if not used
            if nargout < 1
                fclose(fid);
            end

        end


        function fid = writeHeaeder(this, file )
        % Writes the header to the file
        %
        % the function will write the header from a vrs binary file object
        % header
        %
        % if file is a char it will open the file and DISCARD! all content
        % in the file. The function will not close the file if the user
        % specifies a return parameter.
        %
        % Usage:
        %   fid = bin.readHeader();
        %   fid = bin.readHeader(file);
        %   bin.readHeader(...); - this will close the file after writing
        %
        % Parameters:
        %   @param file - @type the file either a char giving the absolute path
        %   or a file pointer (optional)
        % Return values:
        %   @return fid - the file pointer to the file set to the data block ready
        %   to read in the file

            fid = [];
            if nargin < 2 || isempty(file)
                [file ] = vsv.file.FileTools.obtainFilenamePut( [], {'*.vrs', 'VSX Files (*.vrs)'; ...
                                                                  '*.vsx', 'VSX Files (*.vsx)'; ...
                                                                  '*.*',   'All Files (*.*)'},   ...
                                                                  'Pick a vsxHeader file' );
            end

            if isempty( this.version )
                this.version = 1;
            end

            if isempty( this.compression )
                this.compression = 0;
            end

            if isempty( this.isSignedCompress )
                this.isSignedCompress = 0;
            end

            if isempty( this.isTimeTag )
                this.isTimeTag = false;
            end

            if isempty(file)
                return;
            end

            if ischar(file)
                fid = fopen(file, 'w');
            end

            if isempty(fid) || fid <=0
                return;
            end

            try
                [fp, ~, ~] = fileparts(file);

                % if filepath is empty file was obtained from current folder
                if isempty(fp)
                    fullfile(cd, file);
                end

                this.filePath = file;

                frewind(fid);

                fwrite(fid, this.version,           'uint16'); % version
                fwrite(fid, this.compression,       'uint8');  % compression
                fwrite(fid, this.isSignedCompress,  'uint8');  % isSignedCompression

                fwrite(fid, this.isTimeTag, 'uint8'); % time tag
                if( this.isTimeTag > 0)
                    fwrite(fid, this.sec,       'uint8');
                    fwrite(fid, this.min,       'uint8');
                    fwrite(fid, this.hour,      'uint8');
                    fwrite(fid, this.day,       'uint8');
                    fwrite(fid, this.month,     'uint8');
                    fwrite(fid, this.year-1900, 'uint8');
                end

                this.writeStudyString(fid, this.studyID);
                this.writeStudyString(fid, this.sampleID);
                this.writeStudyString(fid, this.comment);

                % make sure tmp is not 0
                if( length(this.dim) ~= 4 )
                    tmp  = ones(1, 4);
                    leng = min( [ length(this.dim), 4]);%#ok not so good, but OK
                    for i = 1:leng
                        tmp(i) = this.dim(i);
                    end
                    this.dim = tmp(1:4);
                end

                fwrite(fid, this.dim, 'uint64');
                if isempty(this.numDataPoints)
                    this.numDataPoints = [];
                end
                fwrite(fid, this.numDataPoints, 'uint64');
                if isempty(this.numDataPoints)
                    this.numDataPoints = [];
                end

                fwrite(fid, this.getNumDataType(), 'uint8' );
            catch
                % something went wrong make sure the file gets closed
                try
                    fclose(fid);
                catch
                end
                fid = -1;
                % rethrow the error
                error(lasterror);%#ok thats OK
            end

            % make sure to close the file if not used
            if nargout < 1
                fclose(fid);
            end

        end
    end


end


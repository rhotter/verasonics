%% VRS binary file import and export example
% The following tutorial demonstrates file import and export of
% storage-class binary files into and from the Matlab work space.
%
% NOTE: This file requires a test vrs file in the current folder. Use the
% storage class example script
% SetUpL11_5vFlashAngles_Export_SequenceStorage_frames to create a vrs
% binary file and place a file in the folder where you want to execut this
% example. Make sure that folder contains a file TestFile030.vrs
%
% $Author: Dr. Daniel Rohrbach,
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%

if ~exist( 'TestFile030.vrs', 'file' )
    error('TestFile030.vrs does not exist');
end

%% Simple file import
% The following demonstrates of how to import binary Storage-class files.
% The verasonics toolbox provides an easy method for loading data into the
% Matlab work space

[data, binaryFile] = vsv.file.readVRSFile('TestFile030.vrs' );
whos data
figure; imagesc( log10( abs( double(data) ) ) ); colormap('gray');

%% Advanced file import
% The binaryFile is a handle to the file and can be used to read files and
% retreive information about the stored data. The user can create their own
% binary file handle objects for loading and exporting data

% generate the file object
binaryFile = vsv.file.VRSBinaryFile;

% read the file
binaryFile.readFile( 'TestFile030.vrs' );

% get header information
headerInfo = binaryFile.getHeaderInfo( )

% get the formatted data. This will format the data to have the right
% dimensions. If the data was compressed the formatted data will be
% automatically uncompressed
data = binaryFile.formatData( );
whos data

figure; imagesc( log10( abs( double(data) ) ) ); colormap('gray');

%% Alternative file import
% The user does not need to specify the file name. If the file name is not
% specified ui-dialog for selecting a binary file will open.
binaryFile.readFile( );

%%
% Alternatively, the binary file object can also be constructed directly
% with a file name. This will automatically load the specific file.
binaryFile = vsv.file.VRSBinaryFile('TestFile030.vrs');
binaryFile.getHeaderInfo()

%%
% If read file is called multiply times using the same file object the
% content in the file object is first cleared and then the new data is
% loaded into the object. Object generation is slow in Matlab. Reusing the
% same file object for several files can increase performance when reading
% several files.
%
% The utility functions also provide tools for efficiently listing files in
% a folder.

% define a folder
folder = cd;
% use regular expressions to list a set of files
% the following will list files that contain the TestFile as a base file
% name followed by exact 3 digits and a file tag .vrs. See regular
% expressions for more details
ls = vsv.file.FileTools.listFilesRegexp(folder, 'TestFile[0-9]{3}\.vrs');

% create one file object
binaryFile = vsv.file.VRSBinaryFile;

% plot all files in a folder to create a short movie zene
figure;
for i = 1:length(ls)
    % read in the next file while reusing the storage object
    binaryFile.readFile( fullfile( folder, ls{i} ) );
    data = binaryFile.formatData();
    imagesc( 10*log10(abs(double(data) ) ) );
    colormap(gray);
    drawnow
end

%% Read only header information
% It is possible to only retrieve  header information from the files without
% reading the data.

header = vsv.file.VRSBinaryHeader();
header.readHeaeder( 'TestFile030.vrs' );

%%
% Reading header information follows the same rules as reading data. The
% following calls are valid as well:
%

% read header with UI support
header.readHeaeder( );

% read header directly while constructing the header object
header = vsv.file.VRSBinaryHeader( 'TestFile030.vrs' );

%%
% Information about the header and file can be retreived directly from the
% header object using dot-notation similar to a struct:

% read the user comment
comment = header.comment;
% read whether compression was used to save the file
dim = header.dim;
% read the day written if isTimeTag was turned on
if header.isTimeTag
    day = header.day;
else
    day = 0;
end

disp( [ 'File comment = ' comment '; Dimension = [' num2str(dim') ']'] );


%%
% The header information can be converted directly to a struct for
% convenience

str = header.toStruct;
disp(str)

%% File Export
%
% The utility functions also provide tools for exporting data from the
% Matlab workspace to a binary file. Currently, only uncompressed data can
% be written using the binary file support in Matlab.
%
% The easiest and safest way of importing data is to use a previously
% exported file and modify its data and save it back to file using the
% Matlab utility functions. A script for importing data can be easily
% generated from a script that exports data.

% import a binary file
binaryFile = vsv.file.VRSBinaryFile('TestFile030.vrs');

% get formatted data
data   = binaryFile.formatData;

% do some data modification, e.g., filter RF-data using an average filter
kernel = [1 1 1 1]';
fdata = filter(kernel, 1, data);

% use set function to set data, convert back to int16 otherwise data will
% be save as double
binaryFile.setData( int16(fdata) );
binaryFile.comment = 'Filtered Data';

% write the file to disk
binaryFile.writeFile( 'FilteredTestFile030.vrs' );




function Trans = computeHvMuxSHIAperture(Trans)
% Create Trans.HVMux.SHIAperture array based on existing
% Trans.HVMux.ApertureES and the HVMuxMap for the probe identified by the
% Trans structure input argument.
%
% This function assumes the input argument Trans structure includes an
% HVMux.ApertureES array that has already been generated and checked for
% validity and compatibility with the probe constraints.  No additional
% validity checking is done here; for each column in HVMux.ApertureES a
% corresponding column in HVMux.SHIAperture will be created, to select all
% elements identified by non-zero entries in the Aperture array.
%
% If HVMux.SHIAperture is empty or does not exist, it will be created with
% a new column for every column in HVMux.ApertureES.  If SHIAperture does
% exist, it will be assumed that all existing entries actually match the
% corresponding entries in HVMux.ApertureES and new columns to SHIAperture
% will be created only for Aperture columns beyond the end of the current
% SHIAperture array.
%
% Note that for probes using predefined SHIAperture tables such as the
% GEL3-12D the number of entries in SHIAperture and Aperture will be equal
% so this function will return without doing anything.
%
% The HVMuxMap table identifies the byte number and bit position within
% that byte that needs to be set to select the element identified by the
% row index of the table.

% revision history:
% May 2020 VTS-343 new "SHIAperture" function replaces "VDAS_Aperture", with
% changes to allow use of HVMux programming tables 256 bytes long instead
% of the original 255 byte limit.
% Note that SHIAperture is a Matlab double; the old VDASAperture was a
% Matlab uint8!!
%
% April 2019 VTS-1201 Check and report error condition if SHIAperture
% table will not fit in available memory on SHI Module
%
% July 2018 VTS-852: add ability for mux UTA's to support probes
% that do not have 1:1 mapping from elements to connector signals (i.e. the
% inputs to the mux switches)

if ~isfield(Trans, 'HVMux')
    % this is not a mux probe so return without doing anything
    return
elseif ~isfield(Trans.HVMux, 'ApertureES') || isempty(Trans.HVMux.ApertureES)
    error('computeHvMuxSHIAperture: Trans.HVMux.ApertureES field missing or empty');
end

numAper = size(Trans.HVMux.ApertureES, 2); % number of Apertures that have been defined
if ~isfield(Trans.HVMux, 'SHIAperture') || isempty(Trans.HVMux.SHIAperture)
    % need to create the entire table
    numSHIAper = 0;
    % create field in case it didn't exist
    Trans.HVMux.SHIAperture = [];
else
    numSHIAper = size(Trans.HVMux.SHIAperture, 2); % number of SHIApertures that already exist
end
if numSHIAper == numAper
    % return since no new SHIApertures are needed
    return
elseif numSHIAper > numAper
    error('computeHvMuxSHIAperture: Trans.HVMux has SHIAperture entries with no corresponding ApertureES');
end

% There are some new SHI apertures to create, so get the mux table
if isfield(Trans.HVMux, 'utaSF') && Trans.HVMux.utaSF > 0
    % This is a UTA HVMux array, not an HVMux probe
    switch Trans.HVMux.utaSF
        case 2
            % This is the UTA 260-Mux so use MuxMap for it, not the probe
            [MuxMap, SHIAperLgth] = getMuxMap('UTA260Mux');
        case 3
            % This is the UTA 1024-Mux so use MuxMap for it, not the probe
            [MuxMap, SHIAperLgth] = getMuxMap('UTA1024Mux');
        case 6
            % This is the UTA 1024-Mux so use MuxMap for it, not the probe
            [MuxMap, SHIAperLgth] = getMuxMap('SHI2048Mux');
        otherwise
            error('computeHvMuxSHIAperture: Unrecognized UTA special feature index in Trans.HVMux.utaSF.');
    end
else
    % Not a UTA with HVMux, so this is an HVMux probe; get the following
    % two fields from the probe's Trans structure
    MuxMap = Trans.HVMux.MuxMap;
    SHIAperLgth = Trans.HVMux.SHIAperLgth;
end
if numSHIAper > 0 && SHIAperLgth ~= size(Trans.HVMux.SHIAperture, 1)
    % required number of entries in each column of SHIAperture does not
    % match existng table
    error('computeHvMuxSHIAperture: Trans.HVMux has SHIAperture entries with incorrect length');
end

% VTS-1201 determine size of SHIAperture table and compare to available
% memory; don't waste time querying the SHI if memory required is less than
% half of what's in all early SHI or UTA baseboard modules
% Note tables stored in SHI add one more byte at beginning to specify the
% length, so add one to SHIAperLgth to account for that.
memRqd = (SHIAperLgth+1) * numAper;
if memRqd > 2^14
    % Import Java library packages so we can make Hal query
    import com.verasonics.hal.hardware.*
    import com.verasonics.hal.shi.*
    % attempt to open the HW to query for memory size
    hwOpenResult = Hardware.openHardware(false);
    if(hwOpenResult == HwOpenResult.success)
        result = Hardware.getInt(IntAttr.shiHvMuxRamInstalledInKb);
        memSize = result*1024; % convert to bytes from Kib
    else
        % use default value
        memSize = 2^15; % memory available on all early SHI and UTA baseboard modules
    end
    if memRqd > memSize
        fprintf(2, '\n ERROR!\ncomputeHvMuxSHIAperture: Insufficient memory space for HvMuxSHIaperture Tables.\n');
        fprintf(2, '%d Apertures of length %d bytes exceeds available memory space of %d KiBytes.\n', numAper, (SHIAperLgth + 1), memSize/1024);
        error(' ');
    end
end

% create the new SHIAperture entries
for apnum = (numSHIAper+1):numAper
    [activeEL, ~, activeES] = find(Trans.HVMux.ApertureES(:, apnum));
    if Trans.HVMux.utaSF == 0
        % This is an HVMux probe, where the inputs to the HVMux switches
        % are mapped directly to elements. In this case activeEL, the
        % indices of the non-zero values in the HVMux.Aperture array, tell
        % us which HVMux switches to turn on.  The signals at the connector
        % pins are the outputs of the HVMux switch array.
        % UTA.TransConnector will map those signals to hardware channels.
        MuxIn = activeEL;
    else
        % This is a UTA module with it's own HVMux switching, so the inputs
        % to the HVMux switches are the signals at the connector pins which
        % may not have a one-to-one mapping to the elements.  That mapping
        % is specified by Trans.Connector, which has been used to set the
        % values in HVMux.Aperture. In this case activeCS, the actual
        % non-zero values in the HVMux.Aperture array, tell us which HVMux
        % switches to turn on. UTA.TransConnector defines the mapping
        % through the UTA module mux switches from connector signals to
        % hardware channels at the output of the HVMux switch array.
        MuxIn = activeES;
    end
    NewSHIAperCol = zeros(SHIAperLgth,1); % create the new column filled with zeros
    for i = 1:length(MuxIn)
        elnum = MuxIn(i);
        if MuxMap(elnum, 2) > 0
            NewSHIAperCol(MuxMap(elnum, 1)) = uint8(bitset(NewSHIAperCol(MuxMap(elnum, 1)), MuxMap(elnum, 2)));
        end
    end
    Trans.HVMux.SHIAperture = double([Trans.HVMux.SHIAperture, NewSHIAperCol]);
end

end

function [ UTA ] = computeUTA( UTAtype, Connector )
% computeUTA: define the connector pin assignments, channel to pin mapping,
% and other characteristics of a UTA adapter module
% Copyright 2001-2023 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
%
% USAGE:    UTA = computeUTA( UTAtype, Connector )
%
% Input arguments
%   'UTAtype'   is the 1 X 4 UTA ID code read from the HW system
% 	'Connector' identifies which connector(s) are to be used if the
%               UTA module supports more than one.  This argument is optional and
%               will default to 1 if not provided.  When provided, it is
%               expected to be in the format defined for
%               Resurce.Parameters.Connector in a user SetUp script
%
% Output arguments
%   'UTA'       structure with fields as defined below:
%       UTA.UTAname: a string with the Marketing name of the UTA Module,
%               for example 'UTA 260-MUX'
%       UTA.UTAdesc: a string that can be appended after the UTAname string,
%               to provide more details of connector type(s), number of
%               channels, etc. (used in Version output listing)
%       UTA.mfrPartNum: a string with the manufacturer's name and part
%               number for the connector being used.
%       UTA.numCh: total number of channels required in the HW system.
%               The field Resource.Parameters.numTransmit will be set equal
%               to UTA.numCh. NOTE: the actual number of channels available
%               at connector pins may be less than numCh, since numCh must
%               include all channels from each active Channel Group even if
%               some of them are unused. (Due to the constraints of the Hal
%               DMA generation SW)
%       UTA.activeCG: a 1 X 16 array treated as logicals, with each entry
%               identifying whether the associated channel group contains any
%               active channels.  The first entry applies to the first CG in the
%               first acquisition board slot, and the last entry applies to the
%               second CG on the fourth acq board slot. VTS-1600 8 board
%               slots
%       UTA.ChPinMap: an R X C array of doubles, where R and C are the number
%               of rows and columns in the connector manufacturer's
%               representation of the physical arrangement of pins in the
%               connector.  Each entry in the array identifies the signal
%               connected to that connector pin in the UTA module:
%                  = -1 identifies any pin that is not a ground or
%                   an element signal (could be a no-connect, control signal,
%                   power supply, etc.).  Refer to the UTA module
%                   specification for names and definitions of these signals.
%                  = 0 identifies a ground pin
%                  = 1 to numCh, identifies the
%                   element signal number connected to that pin at the
%                   connector, based on the connector pinout
%       UTA.ChPinNames: a numCh X 2 cell array of string values, derived from
%               UTA.ChPinMap.  The first column contains the active system
%               channel numbers present at the connector, in ascending numeric
%               order. The second column identifies the connector pin name/number
%               to which that channel signal is connected, using the pin naming
%               convention defined by the connector manufacturer.
%       UTA.TransConnector: The mapping through the UTA module from element
%               signals at the connector to the active channels in the
%               hardware system.  To create the overall mapping from
%               transducer elements to system channels, VSX initialization
%               will combine Trans.ConnectorES with UTA.TransConnector to
%               create Trans.Connector as follows:
%                  Trans.Connector = UTA.TransConnector(Trans.ConnectorES);

% Revision History
% Jan. 2020: VTS-1600 V-512 system support (8 board slots)
% Nov. 2019: VTS-1513 Hal query didn't work if no HW system present
% Oct. 2018: VTS-960 improved support for UTA configuration checks and
%            reporting
% Aug. 2018: VTS-872 support alternate use deriving numCh from boards that
%            are present when Trans structure not defined
% May  2018: VTS-808 eliminating use of Resource.Parameters.numTransmit and
%            deriving channel usage from Trans.Connector
%            VTS-770 added recognition of six new UTA modules
% Dec. 2017: Added FUS Elite 3000 3D use of UTA 260-D with Vantage 64, 64 LE, and 128 systems.
% Aug. 2017: Add new field UTA.elBiasEna (VTS-598)

import com.verasonics.hal.hardware.*

if nargin == 1
    Connector = 1; % default to first connector if not specified
elseif max(Connector) > UTAtype(3)
    fprintf(2, 'ERROR computeUTA: requested connector number exceeds connector count for specified UTA module.\n');
    error('computeUTA exiting');
end

if evalin('base', 'exist(''Trans'', ''var'') && isfield(Trans, ''numelements'') && ~isempty(Trans.numelements)')
    % Trans structure with numelements exists so use it to determine number
    % of active element signals at the connector interface
    Trans = evalin('base', 'Trans');
    if isfield(Trans, 'ConnectorES') && ~isempty(Trans.ConnectorES)
        % if ConnectorES exists, use it to get a more accurate indication of
        % number of active element signals at connector
        maxEL = max(Trans.ConnectorES(:, 1)); % use first column, in case parallel channels are being used
    else % assign default using numelements if Trans.Connector does not exist
        maxEL = Trans.numelements;
    end
else
    % Trans does not exist, so this call to computeUTA may be from a test
    % script trying to determine how many channels are available.  Set
    % maxEL to zero as a dummy flag, and query hardware for number of
    % acquisition boards present; each UTA will then fill in a default
    % maxEL value based on numBoards and requested Connector.
    maxEL = 0;
    hwOpenResult = Hardware.openHardware(); % Does nothing if hardware is already open.
    if(hwOpenResult == HwOpenResult.success)
        % hw system is present so query for number of boards
        numBoards  = Hardware.getInt(IntAttr.numAcquisitionBoards);
        if numBoards == 0
            % If no boards are present we can only run in simulate mode;
            % for simulation set numBoards to four so the maximum possible
            % number of channels can be used.
            numBoards = 4;
        end
    else
        % Can't access HW system so again set numBoards to four so the maximum possible
        % number of channels can be used for simulation
        numBoards = 4;
    end

end
% VTS-1600 Initialize UTA.activeCG to all zeros, then set entries for
% active CG's to 1
UTA.activeCG = zeros(1, 16); % 16 CG's on 8 boards

if isequal(UTAtype, [1, 1, 1, 0])
    %% UTA 260-S
    % Single HDI connector (or SHI with one connector)
    UTA.UTAname = 'UTA 260-S';
    UTA.UTAdesc = '  single Verasonics 260 pin, 128 ch. Connector';
    UTA.mfrPartNum = 'Cannon DL5-260  260 pin ZIF connector';

    HDIpins = mapHDIpins(Connector); % subfunction defined after the end of this function
    UTA.ChPinMap = HDIpins.ChPinMap;
    % Determine default maxEL when Trans not defined:
    if maxEL == 0 && numBoards > 1
        % at least two boards are available so we can use all 128 channels
        % at connector
        maxEL = 128;
    end

    if maxEL > 64
        UTA.numCh = 128;
        UTA.activeCG([1, 2, 7, 8]) = 1; % acquisition modules in slots 1 and 4
        UTA.ChPinNames = HDIpins.ChPinNames;
    else
        % using only first 64 elements at connector: Vantage 64 compatible script
        UTA.numCh = 64;
        UTA.activeCG(1:2) = 1; % acquisition module in slot 1 only
        UTA.ChPinNames = cell(64, 2);
        UTA.ChPinNames = HDIpins.ChPinNames(1:64, :);
        ChPins = UTA.ChPinMap;
        ChPins(UTA.ChPinMap > 64) = -1; % force unused channels 65:192 to -1
        UTA.ChPinMap = ChPins;
    end
    UTA.TransConnector = (1:UTA.numCh)'; % 1-to-1 mapping for UTA 260-S and single SHI
    UTA.elBiasEna = 1;

elseif isequal(UTAtype, [1, 1, 1, 2])
    %% UTA 260-MUX
    % Single HDI connector with 64:128 HVMux for use with Vantage 64
    % system (one acquisition board)
    UTA.UTAname = 'UTA 260-MUX';
    UTA.UTAdesc = '  Verasonics 260 pin, 128 ch. Connector with 128:64 HVMux';
    UTA.mfrPartNum = 'Cannon DL5-260  260 pin ZIF connector';
    UTA.numCh = 64;
    UTA.activeCG(1:2) = 1; % single acquisition module in slot 1
    HDIpins = mapHDIpins(Connector); % subfunction defined after the end of this function
    UTA.ChPinMap = HDIpins.ChPinMap;
    UTA.ChPinNames = HDIpins.ChPinNames;
    UTA.TransConnector = ([1:UTA.numCh, 1:UTA.numCh])'; % 2-to-1 HVMux mapping for UTA 260-MUX
    UTA.elBiasEna = 0;


elseif isequal(UTAtype, [1, 1, 2, 0])
    %% UTA 260-D
    % Dual HDI connector (or SHI with dual connectors)
    UTA.UTAname = 'UTA 260-D';
    UTA.UTAdesc = '  dual Verasonics 260 pin, 128 ch. Connectors';
    UTA.mfrPartNum = 'Cannon DL5-260  260 pin ZIF connector';

    HDIpins = mapHDIpins(Connector); % subfunction defined after the end of this function
    UTA.ChPinMap = HDIpins.ChPinMap;

    if isequal(Connector, [1 2])
        % Determine default maxEL when Trans not defined:
        if maxEL == 0 && numBoards == 4
            % all four boards are available so we can use all 128 channels
            % at connector
            maxEL = 256;
        end
        if maxEL > 128
            % using both connectors for 256 channels; mapping is 1-1 to
            % first connector and 1-(1+128) for second connector
            UTA.numCh = 256;
            UTA.activeCG(1:8) = 1; % all four boards in slots 1, 2, 3, 4
            UTA.ChPinNames = HDIpins.ChPinNames;
        else
            % using 260-D with both connectors but only 128 or fewer
            % elements, so assume this is the "FUS Elite-like" utilization
            % of only elements 1:64 on connector 1 and 65:128 on second
            % connector, for compatibility with two-board systems (V-64LE
            % or V-128).
            UTA.numCh = 128;
            UTA.activeCG([1, 2, 7, 8]) = 1; % use only the two boards in slots 1, 4
            UTA.ChPinNames = cell(128, 2);
            UTA.ChPinNames(:, 1) = HDIpins.ChPinNames(1:128, 1);
            UTA.ChPinNames(1:64, 2) = HDIpins.ChPinNames(1:64, 2);
            UTA.ChPinNames(65:128, 2) = HDIpins.ChPinNames(193:256, 2);
        end
    elseif isequal(Connector, 1)
        % left connector by itself
        % Determine default maxEL when Trans not defined:
        if maxEL == 0 && numBoards == 4
            % all four boards are available so we can use all 128 channels
            % at connector
            maxEL = 128;
        end
        if maxEL > 64
            UTA.numCh = 128;
            UTA.activeCG(1:4) = 1; % two boards in slots 1, 2
            UTA.ChPinNames = HDIpins.ChPinNames;
        else
            % using only first 64 elements at connector 1: Vantage 64 compatible script
            UTA.numCh = 64;
            UTA.activeCG(1:2) = 1; % acquisition module in slot 1 only
            UTA.ChPinNames = cell(64, 2);
            UTA.ChPinNames = HDIpins.ChPinNames(1:64, :);
            ChPins = UTA.ChPinMap;
            ChPins(UTA.ChPinMap > 64) = -1; % force unused channels 65:192 to -1
            UTA.ChPinMap = ChPins;
        end
    elseif isequal(Connector, 2)
        % right connector
        UTA.numCh = 128;
        UTA.activeCG(5:8) = 1; % two boards in slots 3, 4
        UTA.ChPinNames = HDIpins.ChPinNames;
    end
    UTA.TransConnector = (1:UTA.numCh)'; % 1-to-1 mapping for UTA 260-D and dual SHI, for all three connector selections
    UTA.elBiasEna = 1;

elseif isequal(UTAtype, [1, 2, 1, 0])
    %% breakout board or custom adapter; define 1:1 mappings for number of
    % boards actually present.
    UTA.UTAname = 'Custom or breakout board';
    UTA.UTAdesc = [];
    UTA.mfrPartNum = [];
    UTA.numCh = [];
    UTA.activeCG = [];
    UTA.ChPinMap = [];
    UTA.ChPinNames = [];
    UTA.TransConnector = [];
    UTA.elBiasEna = [];


elseif isequal(UTAtype, [1, 3, 1, 0])
    %% UTA 360 MSProbes
    % Cannon 360 pin ZIF Connector
    UTA.UTAname = 'UTA 360';
    UTA.UTAdesc = '  MS series 360 pin, 256 ch. ZIF Connector';
    UTA.mfrPartNum = 'ITT Cannon DLM6-360  360 pin ZIF connector';
    % This UTA module cannot be used on any system configuration other than
    % Vantage 256, so always enable all 8 CG's to enforce that.
    UTA.numCh = 256;
    UTA.activeCG(1:8) = 1; % all four boards in slots 1, 2, 3, 4
    UTA.elBiasEna = 0;

    UTA.TransConnector  =  [253 063 236 059 189 062 161 036 175 116 162 098 243 125 247 099 ...
                            238 039 184 048 188 106 245 058 166 117 235 111 165 040 167 110 ...
                            221 031 204 027 157 030 129 004 143 084 130 066 211 093 215 067 ...
                            206 007 152 016 156 074 213 026 134 085 203 079 133 008 135 078 ...
                            191 061 173 121 192 124 225 101 244 115 164 100 179 064 252 037 ...
                            176 041 251 113 190 042 250 057 168 053 174 047 169 105 163 046 ...
                            159 029 141 089 160 092 193 069 212 083 132 068 147 032 220 005 ...
                            144 009 219 081 158 010 218 025 136 021 142 015 137 073 131 014 ...
                            255 127 240 120 256 126 228 097 180 118 231 102 239 128 249 103 ...
                            242 043 246 045 254 104 248 119 170 050 171 051 232 044 227 114 ...
                            223 095 208 088 224 094 196 065 148 086 199 070 207 096 217 071 ...
                            210 011 214 013 222 072 216 087 138 018 139 019 200 012 195 082 ...
                            186 123 172 056 187 122 230 033 177 054 226 034 241 060 183 035 ...
                            178 107 182 109 185 038 181 055 234 052 237 049 233 108 229 112 ...
                            154 091 140 024 155 090 198 001 145 022 194 002 209 028 151 003 ...
                            146 075 150 077 153 006 149 023 202 020 205 017 201 076 197 080 ]';


    % Connector pin mapping from UTA module spec

    UTA.ChPinMap = [67     3    66     2    65     1   193   129   194   130   195   131
                     5     0    68     0     4     0   196     0   132     0   197     0
                     8    71     7    70     6    69   133   198   134   199     0   135
                     0    73     0     9     0    72     0   200     0   136     0     0
                    76    12    75    11    74    10   201   137   202   138     0     0
                    79    15    78    14    77    13   203   139     0     0   204   140
                    17     0    80     0    16     0   205     0     0     0   141     0
                    83    19     0    82    18    81   206   142   207     0   143   208
                     0    84     0     0     0    20     0   144     0   209     0   145
                    86    22     0     0    85    21   210   146   211   147   212   148
                     0     0    88    24    87    23   213   149   214   150   215   151
                     0     0    89     0    25     0   216     0   152     0   217     0
                    28     0    91    27    90    26   153   218   154   219   155   220
                     0    93     0    29     0    92     0   156     0   221     0   157
                    96    32    95    31    94    30   222   158   223   159   224   160
                    99    35    98    34    97    33   225   161   226   162   227   163
                    37     0   100     0    36     0   228     0   164     0   229     0
                    40   103    39   102    38   101   165   230   166   231     0   167
                     0   105     0    41     0   104     0   232     0   168     0     0
                   108    44   107    43   106    42   233   169   234   170     0     0
                   111    47   110    46   109    45   235   171     0     0   236   172
                    49     0   112     0    48     0   237     0     0     0   173     0
                   115    51     0   114    50   113   238   174   239     0   175   240
                     0   116     0     0     0    52     0   176     0   241     0   177
                   118    54     0     0   117    53   242   178   243   179   244   180
                     0     0   120    56   119    55   245   181   246   182   247   183
                     0     0   121     0    57     0   248     0   184     0   249     0
                    60     0   123    59   122    58   185   250   186   251   187   252
                     0   125     0    61     0   124     0   188     0   253     0   189
                   128    64   127    63   126    62   254   190   255   191   256   192];

    ChPins = UTA.ChPinMap .* (UTA.ChPinMap > 0);
    [row, col, chnum] = find(ChPins);
    ChRow = zeros(256, 1);
    ChRow(chnum(1:256)) = row(1:256);
    ChCol = zeros(256, 1);
    ChCol(chnum(1:256)) = col(1:256);

    Rowname = 'ABCDEFGHJKLMNPRSTUVWXYZabcdefg';
    UTA.ChPinNames = cell(256, 2);
    for i=1:256
        UTA.ChPinNames{i, 1} = ['CH ', num2str(i, '%03d')];
        UTA.ChPinNames{i, 2} = [Rowname(ChRow(i)), num2str(ChCol(i), '%03d')];
    end

elseif isequal(UTAtype, [1, 4, 1, 0])
    %% UTA 408  Verasonics
    % Verasonics 408 pin connector
    UTA.UTAname = 'UTA 408';
    UTA.UTAdesc = '  Verasonics 408 pin, 256 ch. Connector';
    UTA.mfrPartNum = 'ITT Cannon DLP408R 408 pin ZIF connector';

    % create an array listing channel assignments to connector pins as
    % given in Marc's UTA module specifiction
    % a value of 0 is a ground pin
    % a value of -1 is any other signal pin or NC
    UTA.ChPinMap = [ -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1;...
                     -1  -1 143 212 154 223 165 234 176 245 187 256;...
                     -1  -1 138 207 149 218 160 229 171 240 182 251;...
                     -1  -1   0   0   0   0   0   0   0   0   0   0;...
                     -1  -1 202 144 213 155 224 166 235 177 246 188;...
                     -1  -1 206 148 217 159 228 170 239 181 250 192;...
                     -1   0   0   0   0   0   0   0   0   0   0   0;...
                    133 201 142 211 153 222 164 233 175 244 186 255;...
                    129 197 139 208 150 219 161 230 172 241 183 252;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                    193 134 203 145 214 156 225 167 236 178 247 189;...
                    196 137 205 147 216 158 227 169 238 180 249 191;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                    132 200 141 210 152 221 163 232 174 243 185 254;...
                    130 198 140 209 151 220 162 231 173 242 184 253;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                    194 135 204 146 215 157 226 168 237 179 248 190;...
                    195 136  12  82  23  93  34 104  45 115  56 126;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                    131 199  77  18  88  29  99  40 110  51 121  62;...
                     67   7  76  17  87  28  98  39 109  50 120  61;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                      2  71  11  81  22  92  33 103  44 114  55 125;...
                      3  72  13  83  24  94  35 105  46 116  57 127;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                     68   8  78  19  89  30 100  41 111  52 122  63;...
                     66   6  75  16  86  27  97  38 108  49 119  60;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                      1  70  10  80  21  91  32 102  43 113  54 124;...
                      4  73  14  84  25  95  36 106  47 117  58 128;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                     69   9  79  20  90  31 101  42 112  53 123  64;...
                     65   5  74  15  85  26  96  37 107  48 118  59;...
                      0   0   0   0   0   0   0   0   0   0   0   0];

    % Determine default maxEL when Trans not defined:
    if maxEL == 0 && numBoards == 4
        % all four boards are available so we can use all 256 channels
        % at connector
        maxEL = 256;
    end
    if maxEL > 128
        % using all 256 element pins on a four board system
        ChPins = UTA.ChPinMap .* (UTA.ChPinMap > 0);
        UTA.numCh = 256;
        UTA.activeCG(1:8) = 1; % all four boards in slots 1, 2, 3, 4
    else
        % Only using the first 128 connector elements; this allows use with
        % a V-128 or V-64LE two-board configuration as well as V-256 with
        % only the CG's in board slots 1 and 4 enabled
        UTA.numCh = 128;
        UTA.activeCG([1, 2, 7, 8]) = 1; % two boards in slots 1, 4
        ChPins = UTA.ChPinMap;
        ChPins(UTA.ChPinMap > 64 & UTA.ChPinMap < 193) = -1; % force unused channels 65:192 to -1
        ChPins(UTA.ChPinMap > 192) = ChPins(UTA.ChPinMap > 192) - 128; % change slot 4 channel numbers from 193:256 to 65:128
        UTA.ChPinMap = ChPins;
        ChPins = UTA.ChPinMap .* (UTA.ChPinMap > 0); % force -1 entries to zero
    end
    UTA.TransConnector = (1:UTA.numCh)'; % 1-to-1 mapping for UTA 408 Verasonics, for both 128 and 256 element configurations
    UTA.elBiasEna = 0; % will change to 2 in 3.4 release, when SW support for the baseboard element bias supply is added.

    nch = UTA.numCh;
    [row, col, chnum] = find(ChPins);
    ChRow = zeros(nch, 1);
    ChRow(chnum(1:nch)) = row(1:nch);
    ChCol = zeros(nch, 1);
    ChCol(chnum(1:nch)) = col(1:nch);

    Colname = 'ABCDEFGHJKLM';
    UTA.ChPinNames = cell(nch, 2);
    for i=1:nch
        UTA.ChPinNames{i, 1} = ['CH ', num2str(i, '%03d')];
        UTA.ChPinNames{i, 2} = [Colname(ChCol(i)), num2str(ChRow(i), '%03d')];
    end

elseif isequal(UTAtype, [1, 4, 4, 0])
    %% UTA Adapter Test STE
    % Verasonics Adapter STE "faked" on-board transducer.
    % provides controllable open/termination for all elements
    % and selectability of any connector.  Scales for any qty of
    % system channels based on what's present.
    UTA.UTAname = 'Verasonics Adapter STE';
    UTA.UTAdesc = ' with Faked Connectors';
    UTA.mfrPartNum = 'P01567-01';
    UTA.numCh = 0;
    % note UTA.activeCG already initialized to all zeros; set desired
    % entries below (VTS-1600)
    for i = 1:4
      if com.verasonics.hal.hardware.Hardware.isAcqBoardDetected(i - 1)  % if board in slot, add 64ch and turn on the CGs
        UTA.numCh = UTA.numCh + 64;
        UTA.activeCG(2*i-1 : 2*i) = [1 1];
      end
    end
    UTA.ChPinMap = [ -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1;...
                     -1  -1 143 212 154 223 165 234 176 245 187 256;...
                     -1  -1 138 207 149 218 160 229 171 240 182 251;...
                     -1  -1   0   0   0   0   0   0   0   0   0   0;...
                     -1  -1 202 144 213 155 224 166 235 177 246 188;...
                     -1  -1 206 148 217 159 228 170 239 181 250 192;...
                     -1   0   0   0   0   0   0   0   0   0   0   0;...
                    133 201 142 211 153 222 164 233 175 244 186 255;...
                    129 197 139 208 150 219 161 230 172 241 183 252;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                    193 134 203 145 214 156 225 167 236 178 247 189;...
                    196 137 205 147 216 158 227 169 238 180 249 191;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                    132 200 141 210 152 221 163 232 174 243 185 254;...
                    130 198 140 209 151 220 162 231 173 242 184 253;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                    194 135 204 146 215 157 226 168 237 179 248 190;...
                    195 136  12  82  23  93  34 104  45 115  56 126;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                    131 199  77  18  88  29  99  40 110  51 121  62;...
                     67   7  76  17  87  28  98  39 109  50 120  61;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                      2  71  11  81  22  92  33 103  44 114  55 125;...
                      3  72  13  83  24  94  35 105  46 116  57 127;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                     68   8  78  19  89  30 100  41 111  52 122  63;...
                     66   6  75  16  86  27  97  38 108  49 119  60;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                      1  70  10  80  21  91  32 102  43 113  54 124;...
                      4  73  14  84  25  95  36 106  47 117  58 128;...
                      0   0   0   0   0   0   0   0   0   0   0   0;...
                     69   9  79  20  90  31 101  42 112  53 123  64;...
                     65   5  74  15  85  26  96  37 107  48 118  59;...
                      0   0   0   0   0   0   0   0   0   0   0   0];

    ChPins = UTA.ChPinMap .* (UTA.ChPinMap > 0);
    [row, col, chnum] = find(ChPins);
    ChRow = zeros(256, 1);
    ChRow(chnum(1:256)) = row(1:256);
    ChCol = zeros(256, 1);
    ChCol(chnum(1:256)) = col(1:256);

    Colname = 'ABCDEFGHJKLM';
    UTA.ChPinNames = cell(256, 2);
    for i=1:256
        UTA.ChPinNames{i, 1} = ['CH ', num2str(i, '%03d')];
        UTA.ChPinNames{i, 2} = [Colname(ChCol(i)), num2str(ChRow(i), '%03d')];
    end
    UTA.TransConnector = (1:UTA.numCh)'; % 1-to-1 mapping all test configurations (elements don't really exist)
    UTA.elBiasEna = 0;

elseif isequal(UTAtype, [1, 4, 1, 6])
    %% HIFU Adapter STE
    % Verasonics HIFU Adapter STE with "faked" on-board transducer.
    % Provides transmit load resistors coupled to large Heat sinks with
    % built-in temperature sensors, to allow testing HIFU systems at full
    % HIFU transmit power levels.  Scales for any qty of system channels
    % based on what's present.  Special feature value of 6 identifies the
    % HIFU adapter STE to differentiate it from ordinary UTA adapter modules
    UTA.UTAname = 'HIFU Adapter STE';
    UTA.UTAdesc = ' with Faked Connectors, 1:1 channel mapping for all board counts';
    UTA.mfrPartNum = 'TBD';
    UTA.numCh = 0;
    UTA.elBiasEna = 0;
    % note UTA.activeCG already initialized to all zeros; set desired
    % entries below
    if maxEL == 0
        maxEL = 64 * numBoards;
    end
    if maxEL > 128
        % enable 256 channels; requires V-256
        UTA.numCh = 256;
        UTA.activeCG(1:8) = 1; % all four boards in slots 1, 2, 3, 4
        UTA.TransConnector = (1:256)';
    elseif maxEL > 64
        % enable 128 channels; requires V-128 or 64 LE
        UTA.numCh = 128;
        UTA.activeCG([1, 2, 7, 8]) = 1; % two boards in slots 1, 4
        UTA.TransConnector = (1:128)';
    else
        % enable 64 channels; requires V-64 or 32 LE
        UTA.numCh = 64;
        UTA.activeCG(1:2) = 1; % one board in slot 1
        UTA.TransConnector = (1:64)';
    end


elseif isequal(UTAtype, [1, 6, 3, 1]) || isequal(UTAtype, [1, 12, 3, 1])
    %% UTA 160-DH/32 Lemo (released version)
    % Dual Hypertac connectors with relays for selecting 32
    % single-element Lemo connectors.  Connector type 6 supports disconnect
    % sensing, but with problems at high transmit voltages due to poor
    % ground connections in the probe (see VTS-792).  Connector type 12 is
    % the same module, but with the disconnect sense pin shorted to ground
    % so it can no longer sense connect/ disconnect
    UTA.UTAname = 'UTA 160-DH/32 LEMO';
    UTA.UTAdesc = '  dual Hypertac 160 pin, 128 ch. plus 32 Lemo coax connectors';
    UTA.mfrPartNum = 'Hypertac HJ501133019  160 pin ZIF connector';
    UTA.elBiasEna = 0;

    %  Hypertac Connector pin labelling as viewed from in front of the
    %  system; this table shows the element signal number at each signal
    %  pin and zero for the ground pins in rows E and F.  The three columns
    %  9, 10, 11 do not exist- they represent the keyway in the middle of
    %  the connector.
% column # 19  18  17  16  15  14  13  12 11 10  9   8   7   6   5   4   3   2   1      row #
    EsPinMap = [ ...
            4   8  12  16  20  24  28  32  0  0  0  36  40  44  48  52  56  60  64; ... % K
            3   7  11  15  19  23  27  31  0  0  0  35  39  43  47  51  55  59  63; ... % J
            2   6  10  14  18  22  26  30  0  0  0  34  38  42  46  50  54  58  62; ... % H
            1   5   9  13  17  21  25  29  0  0  0  33  37  41  45  49  53  57  61; ... % G
            0   0   0   0   0   0   0   0  0  0  0   0   0   0   0   0   0   0   0; ... % F
            0   0   0   0   0   0   0   0  0  0  0   0   0   0   0   0   0   0   0; ... % E
           65  69  73  77  81  85  89  93  0  0  0  97 101 105 109 113 117 121 125; ... % D
           66  70  74  78  82  86  90  94  0  0  0  98 102 106 110 114 118 122 126; ... % C
           67  71  75  79  83  87  91  95  0  0  0  99 103 107 111 115 119 123 127; ... % B
           68  72  76  80  84  88  92  96  0  0  0 100 104 108 112 116 120 124 128];    % A
% column # 19  18  17  16  15  14  13  12 11 10  9   8   7   6   5   4   3   2   1      row #

    % Assign pin name for each element signal
    numES = max(EsPinMap(:)); % number of elements signals in the map
    [row, col, esnum] = find(EsPinMap);
    EsRow = zeros(numES, 1);
    EsRow(esnum(1:numES)) = row(1:numES);
    EsCol = zeros(numES, 1);
    EsCol(esnum(1:numES)) = col(1:numES);
    rowname = 'KJHGFEDCBA';
    EsPinNameConn = cell(numES, 1);
    for i=1:numES
        EsPinNameConn{i} = [num2str(20-EsCol(i), '%d'), rowname(EsRow(i))];
    end


    % Element signal pin map for the array of 32 Lemo connectors as viewed
    % from the front of the system.  This is connector number 3
    EsPinMapLemo = [ ...
            0 0 0 0 0 4   8  12  16  20  24  28  32 0 0 0 0 0 0; ...
            0 0 0 0 0 3   7  11  15  19  23  27  31 0 0 0 0 0 0; ...
            0 0 0 0 0 2   6  10  14  18  22  26  30 0 0 0 0 0 0; ...
            0 0 0 0 0 1   5   9  13  17  21  25  29 0 0 0 0 0 0];

    if isequal(Connector, [1 2]) || isequal(Connector, [1 2 3])
        % using both connectors with or without Lemo's;
        % Determine default maxEL when Trans not defined:
        if maxEL == 0 && numBoards == 4
            % all four boards are available so we can use all 128 channels
            % at each connector
            maxEL = 256;
        elseif maxEL == 0 && numBoards == 2
            % two boards are available so we can use 64 channels
            % at each connector
            maxEL = 128;
        end
        if maxEL > 128
            % when using more than 128, the presumption is more than 64 are
            % being enabled on at least one connector so we enable all 256
            % channels, requiring a V-256 system
            UTA.numCh = 256;
            UTA.activeCG(1:8) = 1; % both CG's of all four boards
            numbds = 4;
            UTA.TransConnector = [(1:32),(193:224),(65:96),(129:160), (33:64),(225:256),(97:128),(161:192)]';
        elseif maxEL > 64
            % more than 64 but not more than 128, so only the first 64
            % elements will be mapped at each connector, allowing use on
            % two-board systems V-128, V-64LE as well as V-256
            UTA.numCh = 128;
            UTA.activeCG([1, 2, 7, 8]) = 1; % both CG's of boards 1, 4
            EsPinMap = EsPinMap .* (EsPinMap < 65);
            numbds = 2;
            UTA.TransConnector = [(1:32), (65:96), (33:64), (97:128)]';
        else
            % 64 or fewer elements but both connectors selected, so we
            % assume script is using only the first 32 on each connector,
            % allowing use on V-64 and V-32LE as well as all larger systems
            % (using only the board in slot 1)
            UTA.numCh = 64;
            UTA.activeCG(1:2) = 1; % both CG's of board 1
            EsPinMap = EsPinMap .* (EsPinMap < 33);
            numbds = 1;
            UTA.TransConnector = (1:64)';
        end
        EsPinMap2 = EsPinMap + 32*numbds*(EsPinMap > 0);
        EsPinMap = [EsPinMap; -ones(1, 19); EsPinMap2];
        % Create pin name for all active element signals on both connectors
        EsPinName = cell(UTA.numCh, 1);
        numPerCon = UTA.numCh / 2;
        for esnum = 1:numPerCon
            EsPinName{esnum, 1} = ['J1-', EsPinNameConn{esnum, 1}];
            EsPinName{numPerCon + esnum, 1} = ['J2-', EsPinNameConn{esnum, 1}];
        end
    elseif isequal(Connector, 1) || isequal(Connector, [1 3])
        % upper connector only
        % Determine default maxEL when Trans not defined:
        if maxEL == 0 && numBoards == 4
            % all four boards are available so we can use all 128 channels
            % at connector 1
            maxEL = 128;
        elseif maxEL == 0 && numBoards == 2
            % two boards are available so we can use 64 channels
            % at connector 1
            maxEL = 64;
        end
        if maxEL > 64
            % using all 128 element signals at one connector requires a
            % V-256 with all four boards present.
            UTA.numCh = 128;
            UTA.activeCG([1, 3, 5, 7]) = 1; % first CG of all four boards
            UTA.TransConnector = [(1:32), (97:128), (33:96)]';
        elseif maxEL > 32
            % only enable the first 64 element signals at connector,
            % allowing use on V-128 and V-64LE as well as V-256
            UTA.numCh = 64;
            UTA.activeCG([1, 7]) = 1; % first CG of boards 1, 4
            EsPinMap = EsPinMap .* (EsPinMap < 65);
            UTA.TransConnector = (1:64)';
        else
            % only enable the first 32 element signals at connector,
            % allowing use on all systems including V-64 and V-32LE
            UTA.numCh = 32;
            UTA.activeCG(1) = 1; % first CG of board 1
            EsPinMap = EsPinMap .* (EsPinMap < 33);
            UTA.TransConnector = (1:32)';
        end
        % Create pin name for all active element signals
        EsPinName = cell(UTA.numCh, 1);
        for esnum = 1:UTA.numCh
            EsPinName{esnum, 1} = ['J1-', EsPinNameConn{esnum, 1}];
        end
    elseif isequal(Connector, 2)
        % lower connector only
        % Determine default maxEL when Trans not defined:
        if maxEL == 0 && numBoards == 4
            % all four boards are available so we can use all 128 channels
            % at connector 2
            maxEL = 128;
        elseif maxEL == 0 && numBoards == 2
            % two boards are available so we can use 64 channels
            % at connector 2
            maxEL = 64;
        end
        if maxEL > 64
            % using all 128 element signals at one connector requires a
            % V-256 with all four boards present.
            UTA.numCh = 128;
            UTA.activeCG([2, 4, 6, 8]) = 1; % second CG of all four boards
            UTA.TransConnector = [(1:32), (97:128), (33:96)]';
        elseif maxEL > 32
            % only enable the first 64 element signals at connector,
            % allowing use on V-128 and V-64LE as well as V-256
            UTA.numCh = 64;
            UTA.activeCG([2, 8]) = 1; % second CG of boards 1, 4
            EsPinMap = EsPinMap .* (EsPinMap < 65);
            UTA.TransConnector = (1:64)';
        else
            % only enable the first 32 element signals at connector,
            % allowing use on all systems including V-64 and V-32LE
            UTA.numCh = 32;
            UTA.activeCG(2) = 1; % second CG of board 1
            EsPinMap = EsPinMap .* (EsPinMap < 33);
            UTA.TransConnector = (1:32)';
        end
        % Create pin name for all active element signals
        EsPinName = cell(UTA.numCh, 1);
        for esnum = 1:UTA.numCh
            EsPinName{esnum, 1} = ['J2-', EsPinNameConn{esnum, 1}];
        end
    elseif isequal(Connector, [2 3])
        % lower connector plus Lemo's
        % Determine default maxEL when Trans not defined:
        if maxEL == 0 && numBoards == 4
            % all four boards are available so we can use all 128 channels
            % at connector 2 plus 32 Lemo
            maxEL = 160;
        elseif maxEL == 0 && numBoards == 2
            % two boards are available so we can use 64 channels
            % at connector 2 plus 32 Lemo
            maxEL = 96;
        end
        if maxEL > 96
            % enabling all 128 element signals at connector 2, requiring a
            % V-256 with all four boards present.
            UTA.numCh = 160;
            UTA.activeCG([1, 2, 4, 6, 8]) = 1; % second CG of all four boards plus Lemo
            UTA.TransConnector = [(1:64), (129:160), (65:128)]';
        elseif maxEL > 64
            % only enable the first 64 element signals at connector 2,
            % allowing use on V-128 and V-64LE as well as V-256
            UTA.numCh = 96;
            UTA.activeCG([1, 2, 8]) = 1; % second CG of boards 1, 4 plus Lemo
            EsPinMap = EsPinMap .* (EsPinMap < 65);
            UTA.TransConnector = (1:96)';
        else
            % only enable the first 32 element signals at connector 2,
            % allowing use on all systems including V-64 and V-32LE
            UTA.numCh = 64;
            UTA.activeCG(1:2) = 1; % both CG's of board 1
            EsPinMap = EsPinMap .* (EsPinMap < 33);
            UTA.TransConnector = (1:64)';
        end
        EsPinMap = EsPinMap + 32*(EsPinMap > 0); % offset of 32 for Lemos
        EsPinMap = [EsPinMapLemo; -ones(1, 19); EsPinMap];
        % Create pin name for all active element signals at Lemos and J2
        EsPinName = cell(UTA.numCh, 1);
        for esnum = 1:32
            EsPinName{esnum, 1} = ['LEMO #', num2str(esnum, '%d')];
        end
        for esnum = 1:(UTA.numCh - 32)
            EsPinName{esnum+32, 1} = ['J2-', EsPinNameConn{esnum, 1}];
        end
    elseif isequal(Connector, 3)
        % Lemo connectors by themselves; allows use on all systems
        UTA.numCh = 32;
        UTA.activeCG(1) = 1; % first CG on board 1
        % only the first CG board 1 for Lemo's
        EsPinMap = EsPinMapLemo;
        UTA.TransConnector = (1:32)';
        EsPinName = cell(UTA.numCh, 1);
        for esnum = 1:32
            EsPinName{esnum, 1} = ['LEMO #', num2str(esnum, '%d')];
        end
    end

    UTA.EsPinMap = EsPinMap;
    % pin mapping to channels for ChPinMap and ChPinNames
    TC = [-1; 0; UTA.TransConnector];
    % map element signals through TransConnector to get ChPinMap, but
    % still map all zero entries (not element signals) to zero
    UTA.ChPinMap = TC(UTA.EsPinMap + 2);
    % create the ChPinNames array
    UTA.ChPinNames = cell(UTA.numCh, 3);
    for esnum=1:UTA.numCh
        chnum = UTA.TransConnector(esnum);
        UTA.ChPinNames{chnum, 3} = ['ES ', num2str(esnum, '%03d')];
        UTA.ChPinNames{chnum, 2} = EsPinName{esnum, 1};
        UTA.ChPinNames{chnum, 1} = ['CH ', num2str(chnum, '%03d')];
    end


elseif isequal(UTAtype, [1, 12, 5, 4])
    %% UTA 160-SH/8 LEMO
    % Single Hypertac connector (with no disconnect sensing) and four pairs
    % of two relay-switched Lemo connectors.  Relay control is special
    % feature 4
    UTA.UTAname = 'UTA 160-SH/8 LEMO';
    UTA.UTAdesc = '  single Hypertac 160 pin, 128 ch. Connector plus 8 Lemo';
    UTA.mfrPartNum = 'Hypertac HJ501133019  160 pin ZIF connector';
    UTA.elBiasEna = 0;

    %  Hypertac Connector pin labelling as viewed from in front of the
    %  system; this table shows the element signal number at each signal
    %  pin and zero for the ground pins in rows E and F.  The three columns
    %  9, 10, 11 do not exist- they represent the keyway in the middle of
    %  the connector.
% column # 19  18  17  16  15  14  13  12 11 10  9   8   7   6   5   4   3   2   1      row #
    UTA.EsPinMap = [ ...
            4   8  12  16  20  24  28  32  0  0  0  36  40  44  48  52  56  60  64; ... % K
            3   7  11  15  19  23  27  31  0  0  0  35  39  43  47  51  55  59  63; ... % J
            2   6  10  14  18  22  26  30  0  0  0  34  38  42  46  50  54  58  62; ... % H
            1   5   9  13  17  21  25  29  0  0  0  33  37  41  45  49  53  57  61; ... % G
            0   0   0   0   0   0   0   0  0  0  0   0   0   0   0   0   0   0   0; ... % F
            0   0   0   0   0   0   0   0  0  0  0   0   0   0   0   0   0   0   0; ... % E
           65  69  73  77  81  85  89  93  0  0  0  97 101 105 109 113 117 121 125; ... % D
           66  70  74  78  82  86  90  94  0  0  0  98 102 106 110 114 118 122 126; ... % C
           67  71  75  79  83  87  91  95  0  0  0  99 103 107 111 115 119 123 127; ... % B
           68  72  76  80  84  88  92  96  0  0  0 100 104 108 112 116 120 124 128];    % A
% column # 19  18  17  16  15  14  13  12 11 10  9   8   7   6   5   4   3   2   1      row #

    % Determine default maxEL when Trans not defined:
    if maxEL == 0 && numBoards > 1
        % at least two boards are available so we can use all 128 channels
        % at connector
        maxEL = 128;
    end
    if maxEL > 64
        % enable all 128 element signals at Hypertac connector, requiring
        % two-board or four-board system
        UTA.numCh = 128;
        UTA.TransConnector = [(64:-1:1), (128:-1:65)]';
        % 1-to-N reversed mapping separately for each board
        UTA.activeCG([1, 2, 7, 8]) = 1; % acquisition modules in slots 1 and 4
    else
        % using only 64 or fewer connector elements: Vantage 64 and 32LE
        % compatible as well as all larger systems
        UTA.numCh = 64;
        UTA.TransConnector = (UTA.numCh:-1:1)'; % 1-to-N reversed mapping
        UTA.activeCG(1:2) = 1; % acquisition module in slot 1 only
        UTA.EsPinMap = UTA.EsPinMap .* (UTA.EsPinMap < 65);
    end

    % pin mapping to channels for ChPinMap and ChPinNames
    TC = [0; UTA.TransConnector];
    % map element signals through TransConnector to get ChPinMap, but
    % still map all zero entries (not element signals) to zero
    UTA.ChPinMap = TC(UTA.EsPinMap + 1);
    conCH = UTA.numCh;
    [row, col, chnum] = find(UTA.EsPinMap);
    ChRow = zeros(conCH, 1);
    ChRow(chnum(1:conCH)) = row(1:conCH);
    ChCol = zeros(conCH, 1);
    ChCol(chnum(1:conCH)) = col(1:conCH);
    rowname = 'KJHGFEDCBA';
    UTA.ChPinNames = cell(conCH, 3);
    for i=1:conCH
        chnum = UTA.TransConnector(i);
        UTA.ChPinNames{chnum, 3} = ['ES ', num2str(i, '%03d')];
        UTA.ChPinNames{chnum, 2} = ['J1-', num2str(20-ChCol(i), '%d'), rowname(ChRow(i))];
        UTA.ChPinNames{chnum, 1} = ['CH ', num2str(chnum, '%03d')];
    end
elseif isequal(UTAtype, [1, 7, 1, 0])
    %% UTA 408-GE
    % GE 408 pin connector
    UTA.UTAname = 'UTA 408-GE';
    UTA.UTAdesc = '  408 pin, 256 channel Connector (GE D-series probes)';
    UTA.mfrPartNum = 'ITT Cannon DLP408R 408 pin ZIF connector';
    UTA.elBiasEna = 1;

    % create an array listing channel assignments to connector pins as
    % given in Marc's UTA module specifiction
    % a value of 0 is a ground pin
    % a value of -1 is any other signal pin or NC
    UTA.ChPinMap =   [   0     0     0     0     0     0     0     0     0     0     0    -1; ...
                       212     0   132     0   144   173   172   225   232   253    -1    -1; ...
                       208   201   220   133   140   161     0   185     0   249    -1    -1; ...
                       211     0   131     0   143   174   167   226   231   254    -1    -1; ...
                       207   202   219   134   139   162     0   186     0   250    -1    -1; ...
                       210     0   130     0   142   175   166   227   230   255    -1    -1; ...
                       196   203   218   135   138   163     0   187     0   251    -1    -1; ...
                       209     0   129     0   141   176   165   228   229   256    -1    -1; ...
                       195   204   217   136   137   164     0   188     0   252    -1    -1; ...
                       200     0   148     0   160   177   184   244   248   237    -1    -1; ...
                       194   213   221   149   156   168     0   189     0   233    -1    -1; ...
                       199     0   147     0   159   178   183   243   247   238    -1    -1; ...
                       193   214   222   150   155   169     0   190     0   234    -1    -1; ...
                       198     0   146     0   158   179   182   242   246   239    -1    -1; ...
                       206   215   223   151   154   170     0   191     0   235    -1    -1; ...
                       197     0   145     0   157   180   181   241   245   240    -1     0; ...
                       205   216   224   152   153   171     0   192     0   236    -1     0; ...
                        13     0    65     0    76   112   117    40    45    64    -1    -1; ...
                         1    24    29    72    73   107     0   128     0    56    -1    -1; ...
                        14     0    66     0    77   110   118    39    46    62    -1    -1; ...
                         2    23    30    71    74   109     0   127     0    54    -1    -1; ...
                        15     0    67     0    78   108   119    38    47    60     0    -1; ...
                         3    22    31    70    75   111     0   126     0    52     0    -1; ...
                        16     0    68     0    79   106   113    37    48    63    -1    -1; ...
                         4    21    32    69    88   105     0   125     0    50    -1    -1; ...
                         9     0    80     0    93   104   114    36    41    61    -1    -1; ...
                         5    20    25    87    89   103     0   124     0    55    -1    -1; ...
                        10     0    81     0    94   102   115    35    42    59    -1    -1; ...
                         6    19    26    86    90   101     0   123     0    53    -1    -1; ...
                        11     0    82     0    95   100   116    34    43    58    -1    -1; ...
                         7    18    27    85    91    99     0   122     0    51    -1    -1; ...
                        12     0    83     0    96    98   120    33    44    57    -1    -1; ...
                         8    17    28    84    92    97     0   121     0    49    -1    -1; ...
                         0    -1     0     0     0     0     0     0     0     0     0    -1];


% The lines below show GE's assignment of element numbers to connector
% pins, as given in Marc's GE UTA module specification.  Note that GE
% numbers the elements from 0 to 255;  The entry with value 256 is actually
% element number 0, to differentiate it from all the other zero entries
% that actually represent ground pins.  The entries numbered 1 through 255
% represent GE's element numbers.
%    UTA.EL_PinMap = [  0   0   0   0   0   0   0   0   0   0   0  -1;...
%                     159   0 222   0 255  31  62  94  95 127  -1  -1;...
%                     158 190 191 223 254  30   0  63   0 126  -1  -1;...
%                     157   0 220   0 253  29  60  92  93 125  -1  -1;...
%                     156 188 189 221 252  28   0  61   0 124  -1  -1;...
%                     155   0 218   0 251  27  58  90  91 123  -1  -1;...
%                     154 186 187 219 250  26   0  59   0 122  -1  -1;...
%                     153   0 216   0 249  25  56  88  89 121  -1  -1;...
%                     152 184 185 217 248  24   0  57   0 120  -1  -1;...
%                     151   0 214   0 247  23  54  86  87 119  -1  -1;...
%                     150 182 183 215 246  22   0  55   0 118  -1  -1;...
%                     149   0 212   0 245  21  52  84  85 117  -1  -1;...
%                     148 180 181 213 244  20   0  53   0 116  -1  -1;...
%                     147   0 210   0 243  19  50  82  83 115  -1  -1;...
%                     146 178 179 211 242  18   0  51   0 114  -1  -1;...
%                     145   0 208   0 241  17  48  80  81 113  -1   0;...
%                     144 176 177 209 240  16   0  49   0 112  -1   0;...
%                     143   0 206   0 239  15  46  78  79 111  -1  -1;...
%                     142 174 175 207 238  14   0  47   0 110  -1  -1;...
%                     141   0 204   0 237  13  44  76  77 109  -1  -1;...
%                     140 172 173 205 236  12   0  45   0 108  -1  -1;...
%                     139   0 202   0 235  11  42  74  75 107   0  -1;...
%                     138 170 171 203 234  10   0  43   0 106   0  -1;...
%                     137   0 200   0 233   9  40  72  73 105  -1  -1;...
%                     136 168 169 201 232   8   0  41   0 104  -1  -1;...
%                     135   0 198   0 231   7  38  70  71 103  -1  -1;...
%                     134 166 167 199 230   6   0  39   0 102  -1  -1;...
%                     133   0 196   0 229   5  36  68  69 101  -1  -1;...
%                     132 164 165 197 228   4   0  37   0 100  -1  -1;...
%                     131   0 194   0 227   3  34  66  67  99  -1  -1;...
%                     130 162 163 195 226   2   0  35   0  98  -1  -1;...
%                     129   0 192   0 225   1  32  64  65  97  -1  -1;...
%                     128 160 161 193 224   256 0  33   0  96  -1  -1;...
%                       0  -1   0   0   0   0   0   0   0   0   0  -1];

% The lines below can be used to define a generic Trans.Connector array for
% the GE connector, showing the mapping from transducer elements to
% connector channels as defined by the UTA.  Note that the indexing into
% the Trans.Connector array represents channels numbered from 1 to 256, so
% the element index values used in our SW will always be one greater than
% the corresponding element number as defined by GE.  For 128 element (or
% smaller) probes GE uses the element positions 64:191 from the connector
% (elements 65:192 for our counting from 1 instead of zero).  So for a
% two-board system supporting only 128 channels, the Trans.Connector array
% would be the middle 128 entries of the array shown here, with 128
% % subtracted from all entries greater than 192.

    % Determine default maxEL when Trans not defined:
    if maxEL == 0 && numBoards == 4
        % all four boards are available so we can use all 256 channels
        % at connector
        maxEL = 256;
    end

    if maxEL > 128
        % enable all 256 element pins; requires a V-256 four board system
        ChPins = UTA.ChPinMap .* (UTA.ChPinMap > 0);
        UTA.numCh = 256;
        UTA.activeCG(1:8) = 1; % all four boards in slots 1, 2, 3, 4
        UTA.TransConnector = [  97    98    99   100   101   102   103   104   105   106   111   108   109   110   107   112, ...
                               171   180   170   179   169   178   168   177   164   176   163   175   162   174   161   173, ...
                               120   121   116   122   115   123   114   124   113   125   119   126   118   127   117   128, ...
                               181   192   182   191   183   190   184   189   165   188   166   187   167   186   172   185, ...
                                33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                               241   245   242   246   243   247   244   248   228   229   227   230   226   231   225   232, ...
                                49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                               236   240   235   239   234   238   233   237   252   256   251   255   250   254   249   253, ...
                                 8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                               205   197   206   198   193   199   194   200   195   209   196   210   207   211   208   212, ...
                                17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                               216   224   215   223   214   222   213   221   204   217   203   218   202   219   201   220, ...
                                83    84    82    85    81    86    80    87    68    69    67    70    66    71    65    72, ...
                               145   152   146   151   147   150   148   149   129   136   130   135   131   134   132   133, ...
                                92    96    91    95    90    94    89    93    88    79    75    78    74    77    73    76, ...
                               153   157   154   158   155   159   156   160   137   141   138   142   139   143   140   144 ]';
    else
        % define the requirements for a V-128 or V-64LE two-board
        % configuration; note this configuration will also be used on a
        % four-board system when 128 or fewer connector elements are active
        UTA.numCh = 128;
        UTA.activeCG([1, 2, 7, 8]) = 1; % two boards in slots 1, 4
        ChPins = UTA.ChPinMap;
        ChPins(UTA.ChPinMap > 64 & UTA.ChPinMap < 193) = -1; % force unused channels 65:192 to -1
        ChPins(UTA.ChPinMap > 192) = ChPins(UTA.ChPinMap > 192) - 128; % change slot 4 channel numbers from 193:256 to 65:128
        UTA.ChPinMap = ChPins;
        ChPins = UTA.ChPinMap .* (UTA.ChPinMap > 0); % force -1 entries to zero
        UTA.TransConnector = [   33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                                113   117   114   118   115   119   116   120   100   101    99   102    98   103    97   104, ...
                                 49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                                108   112   107   111   106   110   105   109   124   128   123   127   122   126   121   125, ...
                                  8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                                 77    69    78    70    65    71    66    72    67    81    68    82    79    83    80    84, ...
                                 17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                                 88    96    87    95    86    94    85    93    76    89    75    90    74    91    73    92 ]';
    end


    nch = UTA.numCh;
    [row, col, chnum] = find(ChPins);
    ChRow = zeros(nch, 1);
    ChRow(chnum(1:nch)) = row(1:nch);
    ChCol = zeros(nch, 1);
    ChCol(chnum(1:nch)) = col(1:nch);

    Colname = 'ABCDEFGHJKLM';
    UTA.ChPinNames = cell(nch, 2);
    for i=1:nch
        UTA.ChPinNames{i, 1} = ['CH ', num2str(i, '%03d')];
        UTA.ChPinNames{i, 2} = [Colname(ChCol(i)), num2str(ChRow(i), '%03d')];
    end

elseif isequal(UTAtype, [ 1  8  1  3 ])
    %% UTA 1024-Mux
    % note, all aperture and channel mapping is done in the example script
    % Trans.HVMux structure must be added to Trans structure in user script
    % by calling "computeUTAMux1024".
    UTA.UTAname = 'UTA 1024-MUX';
    UTA.UTAdesc = '  1024:256 mux of direct-connect cable';
    UTA.mfrPartNum = 'P01668-01';
    UTA.numCh = 256;
    UTA.activeCG(1:8) = 1; % all four boards in slots 1, 2, 3, 4
    UTA.elBiasEna = 0;
    UTA.TransConnector = repmat((1:256)', 4, 1);

elseif isequal(UTAtype, [ 1  8  1  0 ])
    %% UTA 256-Direct
    % Similar to UTA 1024 MUX but without the muxes; direct connection from
    % probe cable to the 256 system channels.  Can be used with all system
    % configurations, so adjust numCh to match

    if maxEL == 0
        maxEL = 64 * numBoards;
    end
    if maxEL > 128
        % enable 256 channels; requires V-256
        UTA.numCh = 256;
        UTA.activeCG(1:8) = 1; % all four boards in slots 1, 2, 3, 4
        UTA.TransConnector = (1:256)';
    elseif maxEL > 64
        % enable 128 channels; requires V-128 or 64 LE
        UTA.numCh = 128;
        UTA.activeCG([1, 2, 7, 8]) = 1; % two boards in slots 1, 4
        UTA.TransConnector = (1:128)';
    else
        % enable 64 channels; requires V-64 or 32 LE
        UTA.numCh = 64;
        UTA.activeCG(1:2) = 1; % one board in slot 1
        UTA.TransConnector = (1:64)';
    end
    UTA.UTAname = 'UTA 256-DIRECT';
    UTA.UTAdesc = '  256 channel direct connections';
    UTA.mfrPartNum = []; % TBD
    UTA.elBiasEna = 0;

elseif isequal(UTAtype, [ 1  13  1  0 ])
    %% SHI 512-Direct
    % VTS-1600 Preliminary version of V-512 with 512 CH direct connect UTA
    % Direct connection from probe cable to the 512 system channels.  Can
    % be used only with V-512 system configuration

    % enable 512 channels; requires V-512
    UTA.numCh = 512;
    UTA.activeCG = ones(1, 16); % all 8 boards/ 16 CG's must be present
    UTA.TransConnector = (1:512)';

    UTA.UTAname = 'V-512 SHI 512-DIRECT';
    UTA.UTAdesc = '  512 channel direct connections';
    UTA.mfrPartNum = []; % TBD
    UTA.elBiasEna = 0;

elseif isequal(UTAtype, [ 1  13  1  6 ])
    %% SHI 2048-MUX
    % Note: This is a captive SHI module for the V-512 system
    % configuration, mapping 2048 transducer element signals to the
    % system's 512 channels through an array of 4:1 HVMux switches within
    % the SHI module. For scripts using this SHI the Trans.HVMux structure
    % must be added to Trans structure in user script by calling
    % "computeSHIMux2048".
    UTA.UTAname = 'V-512 SHI 2048-MUX';
    UTA.UTAdesc = '  2048:512 mux of direct-connect cable';
    UTA.mfrPartNum = 'TBD';
    UTA.numCh = 512;
    UTA.activeCG(1:16) = 1; % all eight boards in slots 1-8 of V-512 chassis
    UTA.elBiasEna = 0;
    UTA.TransConnector = repmat((1:512)', 4, 1);
    
    
elseif isequal(UTAtype, [1, 9, 1, 0])
    %% UTA 156-U
    % Single Ultrasonix format connector 128 channels
    UTA.UTAname = 'UTA 156-U';
    UTA.UTAdesc = '  single Ultrasonix format 156 pin, 128 ch. Connector';
    UTA.mfrPartNum = 'Cannon DL5-156  156 pin ZIF connector';

    % Note element-to-channel mapping for this module does not allow using
    % only the first 64 elements on a single-board system; this adapter
    % will only be compatible with Vantage 64 LE and higher systems.

    % Ground pin for disconnect sensing is pin # F1 on Cannon DL-156

    % a six-bit jumper ID code is provided through pins A1:A6 of the
    % connector.

    % Ultrasonix numbers the element signals from 0 to 127 on the connector
    % pinout; add one to get the Vantage element numbering used in
    % Trans.Connector etc.

    UTA.numCh = 128;
    UTA.activeCG([1, 2, 7, 8]) = 1; % acquisition modules in slots 1 and 4

    UTA.TransConnector = [ 1   3  67   8   9  11  75  16  17  19  83  24  25  27  91  32 ...
                          33  35  99  40  41  43 107  48  49  51 115  56  57  59 123  64 ...
                           2   4  68  72  10  12  76  80  18  20  84  88  26  28  92  96 ...
                          34  36 100 104  42  44 108 112  50  52 116 120  58  60 124 128 ...
                          66   5   6  70  74  13  14  78  82  21  22  86  90  29  30  94 ...
                          98  37  38 102 106  45  46 110 114  53  54 118 122  61  62 126 ...
                          65  69   7  71  73  77  15  79  81  85  23  87  89  93  31  95 ...
                          97 101  39 103 105 109  47 111 113 117  55 119 121 125  63 127]';
    UTA.elBiasEna = 0;


    UTA.EsPinMap = [  -1  -1  -1  -1  -1  -1; ...
                      -1  -1  -1  -1  -1  -1; ...
                       1  33   0  -1  65  97; ...
                       2  34  66  98  35   3; ...
                      67  99   4  36 100  68; ...
                       0   5  37  69 101   0; ...
                       6  38  70 102  39   7; ...
                      71 103   8  40 104  72; ...
                       9  41   0   0  73 105; ...
                      10  42  74 106  43  11; ...
                      75 107  12  44 108  76; ...
                       0  13  45  77 109   0; ...
                      14  46  78 110  47  15; ...
                      79 111  16  48 112  80; ...
                      17  49   0   0  81 113; ...
                      18  50  82 114  51  19; ...
                      83 115  20  52 116  84; ...
                       0  21  53  85 117   0; ...
                      22  54  86 118  55  23; ...
                      87 119  24  56 120  88; ...
                      25  57   0   0  89 121; ...
                      26  58  90 122  59  27; ...
                      91 123  28  60 124  92; ...
                       0  29  61  93 125   0; ...
                      30  62  94 126  63  31; ...
                      95 127  32  64 128  96];

    EsPinMap = UTA.EsPinMap .* (UTA.EsPinMap > 0);% force -1 entries to zero

    % pin mapping to channels for ChPinMap and ChPinNames
    TC = [0; UTA.TransConnector];
    % map connector signals through TransConnector to get ChPinMap, but
    % still map all zero entries (not element signals) to zero
    UTA.ChPinMap = TC(EsPinMap + 1);
    conCH = UTA.numCh;
    [row, col, chnum] = find(UTA.ChPinMap);
    ChRow = zeros(conCH, 1);
    ChRow(chnum(1:conCH)) = row(1:conCH);
    ChCol = zeros(conCH, 1);
    ChCol(chnum(1:conCH)) = col(1:conCH);
    Rowname = 'ABCDEFGHJKLMNPRSTUVWXYZabc';
    UTA.ChPinNames = cell(conCH, 3);
    for esnum=1:conCH
        chnum = UTA.TransConnector(esnum);
        UTA.ChPinNames{chnum, 1} = ['CH ', num2str(chnum, '%03d')];
        UTA.ChPinNames{chnum, 2} = ['J1-', Rowname(ChRow(chnum)), num2str(ChCol(chnum), '%d')];
        UTA.ChPinNames{chnum, 3} = ['ES ', num2str(esnum, '%03d')];
    end


elseif isequal(UTAtype, [1, 10, 5, 4])
    %% UTA 160-SI/8 LEMO
    % Single I-Pex connector and four pairs of two relay-switched Lemo
    % connectors.  Relay control is special feature 4
    UTA.UTAname = 'UTA 160-SI/8 LEMO';
    UTA.UTAdesc = '  single Ipex 160 pin, 128 ch. Connector plus 8 Lemo';
    UTA.mfrPartNum = 'I-pex part # 30046-160T-03F';

    % The connector signal mapping listed below is as shown in the
    % connector pinout in the module spec: two columns of 80 pins.  The
    % first column is even pins 2 thorugh 160 (top to bottom) and the
    % second column is odd pins 1 through 159.
    UTA.EsPinMap = [  0   0; 112 128;  80  96;  48  64;  16  32;   0   0; 111 127;  79  95;  47  63; 15  31; ...
                      0   0; 110 126;  78  94;  46  62;  14  30;   0   0; 109 125;  77  93;  45  61; 13  29; ...
                      0   0; 108 124;  76  92;  44  60;  12  28;   0   0; 107 123;  75  91;  43  59; 11  27; ...
                      0   0; 106 122;  74  90;  42  58;  10  26; 105 121;  73  89;  41  57;   9  25;  0   0; ...
                      0   0; 104 120;  72  88;  40  56;   8  24; 103 119;  71  87;  39  55;   7  23;  0   0; ...
                    102 118;  70  86;  38  54;   6  22;   0   0; 101 117;  69  85;  37  53;   5  21;  0   0; ...
                    100 116;  68  84;  36  52;   4  20;   0   0;  99 115;  67  83;  35  51;   3  19;  0   0; ...
                     98 114;  66  82;  34  50;   2  18;   0   0;  97 113;  65  81;  33  49;   1  17;  0   0];

    % Determine default maxEL when Trans not defined:
    if maxEL == 0 && numBoards > 1
        % at least two boards are available so we can use all 128 channels
        % at connector
        maxEL = 128;
    end

    if maxEL > 64
        UTA.numCh = 128;
        % connector signal to system channel mapping:
        UTA.TransConnector = [(32:-1:1), (33:96), (128:-1:97)]';
        UTA.activeCG([1, 2, 7, 8]) = 1; % acquisition modules in slots 1 and 4
    else
        % using only 64 or fewer channels: Vantage 64 compatible script
        UTA.numCh = 64;
        UTA.TransConnector = [(32:-1:1), (33:64)]'; % connector signal to system channel mapping
        UTA.activeCG(1:2) = 1; % acquisition module in slot 1 only
        UTA.EsPinMap = UTA.EsPinMap .* (UTA.EsPinMap < 65);
    end
    UTA.elBiasEna = 0;
    % pin mapping to channels for ChPinMap and ChPinNames
    TC = [0; UTA.TransConnector];
    % map connector signals through TransConnector to get ChPinMap, but
    % still map all zero entries (not element signals) to zero
    UTA.ChPinMap = TC(UTA.EsPinMap + 1);
    conCH = UTA.numCh;
    [row, col, chnum] = find(UTA.ChPinMap);
    ChRow = zeros(conCH, 1);
    ChRow(chnum(1:conCH)) = row(1:conCH);
    ChCol = zeros(conCH, 1);
    ChCol(chnum(1:conCH)) = col(1:conCH);
    UTA.ChPinNames = cell(conCH, 3);
    for esnum=1:conCH
        chnum = UTA.TransConnector(esnum);
        UTA.ChPinNames{chnum, 1} = ['CH ', num2str(chnum, '%03d')];
        UTA.ChPinNames{chnum, 2} = ['J1 pin ', num2str((2*ChRow(chnum)+1-ChCol(chnum)), '%d')];
        UTA.ChPinNames{chnum, 3} = ['ES ', num2str(esnum, '%03d')];
    end


elseif isequal(UTAtype, [1, 11, 1, 5])
    %% UTA 64-LEMO
    % 64 individual single-element Lemo coax connectors  Note there is no
    % disconnect sensing and all 64 connectors are treated as a single
    % "connector" that cannot be deselected.  Special feature 5 simply
    % identifies this module as a depopulated version of the UTA 128 LEMO
    UTA.UTAname = 'UTA 64-LEMO';
    UTA.UTAdesc = '  single-element coax connectors';
    UTA.mfrPartNum = 'Lemo EPB.00.250.NTN (Lemo type 00 coax connector)';
    UTA.ChPinMap = [1:8; 9:16; 17:24; 25:32; 33:40; 41:48; 49:56; 57:64];

    UTA.numCh = 64;
    UTA.activeCG(1:2) = 1; % acquisition module in slot 1 only

    UTA.TransConnector = (1:UTA.numCh)'; % 1-to-1 mapping
    UTA.elBiasEna = 0;

    UTA.ChPinNames = cell(UTA.numCh, 2);
    for i=1:UTA.numCh
        UTA.ChPinNames{i, 1} = ['CH ', num2str(i, '%d')];
        UTA.ChPinNames{i, 2} = ['Lemo #', num2str(i, '%d')];
    end


elseif isequal(UTAtype, [1, 11, 1, 0])
    %% UTA 128-LEMO
    % 128 individual single-element Lemo coax connectors  Note there is no
    % disconnect sensing and all 128 connectors are treated as a single
    % "connector" that cannot be deselected
    UTA.UTAname = 'UTA 128-LEMO';
    UTA.UTAdesc = '  single-element coax connectors';
    UTA.mfrPartNum =  'Lemo EPB.00.250.NTN (Lemo type 00 coax connector)';
    % define ChPinMap to represent how they are arranged on front of UTA
    UTA.ChPinMap = [65:72; 73:80; 81:88; 89:96; 1:8; 9:16; 17:24; 25:32; ...
            33:40; 41:48; 49:56; 57:64; 97:104; 105:112; 113:120; 121:128];

    % Determine default maxEL when Trans not defined:
    if maxEL == 0 && numBoards > 1
        % at least two boards are available so we can use all 128 channels
        % at connector
        maxEL = 128;
    end

    if maxEL > 64
        UTA.numCh = 128;
        UTA.activeCG([1, 2, 7, 8]) = 1; % acquisition module in slots 1 and 4
    else
        UTA.numCh = 64;
        UTA.activeCG(1:2) = 1; % acquisition module in slot 1 only
    end

    UTA.TransConnector = (1:UTA.numCh)'; % 1-to-1 mapping
    UTA.elBiasEna = 0;

    UTA.ChPinNames = cell(UTA.numCh, 2);
    for i=1:UTA.numCh
        UTA.ChPinNames{i, 1} = ['CH ', num2str(i, '%d')];
        UTA.ChPinNames{i, 2} = ['Lemo #', num2str(i, '%d')];
    end

else
    % unrecognized UTA configuration
    error(['computeUTA: Unrecognized UTAtype value of [ ', num2str(UTAtype, '%3d'), ' ].']);
end

end

function HDIpins = mapHDIpins(Connector)
% function to create HDI connector pin maps for both single and dual HDI
% adapters

HDIpins.ChPinMap = [  0  -1  -1  -1  -1  -1   0   0;...
                     -1  -1  -1  -1  -1  -1  -1  -1;...
                     -1  -1  -1  -1  -1  -1  -1  -1;...
                     -1  -1  -1  -1  -1  -1  -1  -1;...
                     -1  -1  -1  -1  -1  -1  -1   0;...
                     -1  -1  -1  -1  -1  -1  -1  -1;...
                     -1  -1  -1  -1  -1  -1  -1  -1;...
                     -1  -1  -1  -1  -1  -1  -1  -1;...
                     -1  -1  -1  -1  -1  -1  -1  -1;...
                     93  94  95  96  97  98  99 100;...
                     92  91  90  89 104 103 102 101;...
                     85  86  87  88 105 106 107 108;...
                     84  83  82  81 112 111 110 109;...
                     77  78  79  80 113 114 115 116;...
                     76  75  74  73 120 119 118 117;...
                     69  70  71  72 121 122 123 124;...
                     68  67  66  65 128 127 126 125;...
                      4   3   2   1  64  63  62  61;...
                      5   6   7   8  57  58  59  60;...
                     12  11  10   9  56  55  54  53;...
                     13  14  15  16  49  50  51  52;...
                     20  19  18  17  48  47  46  45;...
                     21  22  23  24  41  42  43  44;...
                     28  27  26  25  40  39  38  37;...
                     29  30  31  32  33  34  35  36;...
                      0  -1  -1  -1  -1  -1  -1   0];

HDIpins.ChPinMap = [zeros(26, 1), HDIpins.ChPinMap, zeros(26, 1)];
ChPins = HDIpins.ChPinMap .* (HDIpins.ChPinMap > 0);
[row, col, chnum] = find(ChPins);
ChRow = zeros(128, 1);
ChRow(chnum(1:128)) = row(1:128);
ChCol = zeros(128, 1);
ChCol(chnum(1:128)) = col(1:128);

Rowname = 'ABCDEFGHJKLMNPRSTUVWXYZabc';
if isequal(Connector, [1 2])
    % using 256 channels through both connectors
    HDIpins.ChPinNames = cell(256, 2);
    Cstr = 'J1-';
    for i=1:128
        HDIpins.ChPinNames{i, 1} = ['CH ', num2str(i, '%03d')];
        HDIpins.ChPinNames{i, 2} = [Cstr, Rowname(ChRow(i)), num2str(ChCol(i), '%d')];
    end
    Cstr = 'J2-';
    for i=1:128
        HDIpins.ChPinNames{i+128, 1} = ['CH ', num2str(i+128, '%03d')];
        HDIpins.ChPinNames{i+128, 2} = [Cstr, Rowname(ChRow(i)), num2str(ChCol(i), '%d')];
    end
elseif Connector > 0
    % using a single connector
    HDIpins.ChPinNames = cell(128, 2);
    Cstr = ['J', num2str(Connector), '-'];
    for i=1:128
        HDIpins.ChPinNames{i, 1} = ['CH ', num2str(i, '%03d')];
        HDIpins.ChPinNames{i, 2} = [Cstr, Rowname(ChRow(i)), num2str(ChCol(i), '%d')];
    end
end

end


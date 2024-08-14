% Copyright 2001-2019 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.
%
% Inputs:
%    SensorList: 1xN list of sensor indexes.
% Outputs:
%    SensorStrings: Nx1 cell array indicating the sensor and the result in text.
%    SensorData: Nx1 double array of readings.
%
% If no input argument is given, the output SensorStrings will be sensor names and output SensorData as index values.
%
function [SensorStrings, SensorData] = GetTemps(SensorList)

import com.verasonics.hal.hardware.*

% Three on TPC.
TpcSensors = { ...
    TpcReg.txBusSepicTemp, NaN, TpcReg.txBusSepicCoolingWindowTemp;
    TpcReg.txGdSensorTemp, NaN, TpcReg.txGdSensorCoolingWindowTemp;
    TpcReg.rxhSensorTemp,  NaN, TpcReg.rxhSensorCoolingWindowTemp;
};

% Six on Backplane.
BkpSensors = { ...
    BkpReg.ambientTemp, BkpReg.configAmbientTempErrorLevel, BkpReg.ambientCoolingWindowTemp;
    BkpReg.exhaustTemp, BkpReg.configExhaustTempErrorLevel, BkpReg.exhaustCoolingWindowTemp;
    BkpReg.pexSmpsTemp, BkpReg.configPexSmpsTempErrorLevel, BkpReg.pexSmpsCoolingWindowTemp;
    BkpReg.pex8748Temp, BkpReg.configPex8748TempErrorLevel, BkpReg.pex8748CoolingWindowTemp;
    BkpReg.tmp431aTemp, BkpReg.configTmp431aTempErrorLevel, BkpReg.tmp431aCoolingWindowTemp;
    BkpReg.pexDieTemp,  BkpReg.configPexDieTempErrorLevel,  BkpReg.pexDieCoolingWindowTemp;
};

% Per board. Mainly eight readings with sixteen NTC individual sensors.
AcqBrdSensors = { ...
    AcqReg.cgd1Temp, AcqReg.configCgd1TempErrorLevel, AcqReg.cgd1CoolingWindowTemp;
    AcqReg.cgd2Temp, AcqReg.configCgd2TempErrorLevel, AcqReg.cgd2CoolingWindowTemp;
    AcqReg.topSmpsTemp, AcqReg.configTopSmpsTempErrorLevel, AcqReg.topSmpsCoolingWindowTemp;
    AcqReg.bottomSmpsTemp, AcqReg.configBottomSmpsTempErrorLevel, AcqReg.bottomSmpsCoolingWindowTemp;
    AcqReg.topShieldTemp, AcqReg.configTopShieldTempErrorLevel, AcqReg.topShieldCoolingWindowTemp;
    AcqReg.bottomShieldTemp, AcqReg.configBottomShieldTempErrorLevel, AcqReg.bottomShieldCoolingWindowTemp;
    AcqReg.topNtcMaxTemp, AcqReg.configTopNtcTempErrorLevel, AcqReg.topNtcMaxCoolingWindowTemp;
    AcqReg.bottomNtcMaxTemp, AcqReg.configBottomNtcTempErrorLevel, AcqReg.bottomNtcMaxCoolingWindowTemp;
    AcqReg.topNtc1Temp, NaN, NaN;
    AcqReg.topNtc2Temp, NaN, NaN;
    AcqReg.topNtc3Temp, NaN, NaN;
    AcqReg.topNtc4Temp, NaN, NaN;
    AcqReg.topNtc5Temp, NaN, NaN;
    AcqReg.topNtc6Temp, NaN, NaN;
    AcqReg.topNtc7Temp, NaN, NaN;
    AcqReg.topNtc8Temp, NaN, NaN;
    AcqReg.bottomNtc1Temp, NaN, NaN;
    AcqReg.bottomNtc2Temp, NaN, NaN;
    AcqReg.bottomNtc3Temp, NaN, NaN;
    AcqReg.bottomNtc4Temp, NaN, NaN;
    AcqReg.bottomNtc5Temp, NaN, NaN;
    AcqReg.bottomNtc6Temp, NaN, NaN;
    AcqReg.bottomNtc7Temp, NaN, NaN;
    AcqReg.bottomNtc8Temp, NaN, NaN;
};

% Get the descriptive string for the enum.
AcqBrdSensors = appendDescription(AcqBrdSensors, 'Acq Board %%d %s');
TpcSensors = appendDescription(TpcSensors, 'TPC Board %s');
BkpSensors = appendDescription(BkpSensors, 'Backplane Board %s');
TpcAndBkpSensors = [TpcSensors; BkpSensors];
AllSensors = [AcqBrdSensors; TpcAndBkpSensors];

acqBrdSensorCnt = size(AcqBrdSensors, 1);
tpcSensorCnt = size(TpcSensors, 1);
bkpSensorCnt = size(BkpSensors, 1);
maxAcqSensorCnt = 4*acqBrdSensorCnt;
maxSensorCnt = maxAcqSensorCnt + tpcSensorCnt + bkpSensorCnt;

if nargin == 0
    SensorList = 1:maxSensorCnt;
    SensorHelp = true;
else
    SensorHelp = false;
end

if Hardware.openHardware(false) ~= HwOpenResult.success
    fprintf(2, 'Could not open Hardware.');
    SensorHelp = true;
end

SensorStrings = cell(numel(SensorList), 1);
if ~SensorHelp
    SensorData = zeros(numel(SensorList), 1);
else
    SensorData = (1:numel(SensorList))';
end

maxStringLen = max(cellfun(@length, AllSensors(:,4)));
padString = @(s, m) sprintf(sprintf('%%-%ds : ', m), s);

rowIdx = 0;
for n=1:maxSensorCnt
    if ~any(SensorList == n) % not listed to be checked.
        continue
    end
    rowIdx = rowIdx+1;
    errorLevel = NaN;
    coolingWindow = [NaN NaN];
    if n > maxAcqSensorCnt
        % TPC or Backplane sensor.
        m = SensorList(rowIdx) - maxAcqSensorCnt;
        SensorStrings{rowIdx} = padString(TpcAndBkpSensors{m,4}, maxStringLen);
        if ~SensorHelp
            tempValue = peekTemperatureValue(TpcAndBkpSensors{m,1});
            if isa(TpcAndBkpSensors{m,2}, 'com.verasonics.hal.hardware.FpgaRegister')
                errorLevel = peekErrorLevel(TpcAndBkpSensors{m,2});
            end
            if isa(TpcAndBkpSensors{m,3}, 'com.verasonics.hal.hardware.FpgaRegister')
                coolingWindow = peekCoolingWindow(TpcAndBkpSensors{m,3});
            end
            SensorData(rowIdx) = tempValue;
            SensorStrings{rowIdx} = [SensorStrings{rowIdx} ...
                formatValueString(tempValue, errorLevel, coolingWindow)];
        end
    else
        % It's on an acquisition board.  Figure out board index 0-3 used by HAL.
        brdIdx = ceil(n/acqBrdSensorCnt) - 1;
        m = SensorList(rowIdx) - acqBrdSensorCnt*brdIdx;
        SensorStrings{rowIdx} =  padString(sprintf(AcqBrdSensors{m,4}, brdIdx+1), maxStringLen);
        if ~SensorHelp
            if ~Hardware.isAcqBoardDetected(brdIdx)
                SensorData(rowIdx) = 0;
                SensorStrings{rowIdx} = sprintf('Acq Board %d is not populated. :  N/A', brdIdx+1);
            else
                tempValue = peekTemperatureValue(brdIdx, AcqBrdSensors{m,1});
                if isa(AcqBrdSensors{m,2}, 'com.verasonics.hal.hardware.FpgaRegister')
                    errorLevel = peekErrorLevel(brdIdx, AcqBrdSensors{m,2});
                end
                if isa(AcqBrdSensors{m,3}, 'com.verasonics.hal.hardware.FpgaRegister')
                    coolingWindow = peekCoolingWindow(brdIdx, AcqBrdSensors{m,3});
                end
                SensorData(rowIdx) = tempValue;
                SensorStrings{rowIdx} = [SensorStrings{rowIdx} ...
                    formatValueString(tempValue, errorLevel, coolingWindow)];
            end
        end
    end
    if SensorHelp
        SensorData(rowIdx) = n;
        SensorStrings{rowIdx} = [SensorStrings{rowIdx} sprintf('Index %d', n)];
    end
end
end

function v = peekTemperatureValue(varargin)
    bitMask = hex2dec('0FFF');
    rv = peekRegister(varargin{:});
    v = bitand(rv, bitMask) / 16;
end

function v = peekErrorLevel(varargin)
    bitMask = hex2dec('0FFF');
    rv = peekRegister(varargin{:});
    v = bitand(rv, bitMask) / 16;
end

function v = peekCoolingWindow(varargin)
    bitMask1 = hex2dec('00FF');
    bitMask2 = hex2dec('FF00');
    rv = peekRegister(varargin{:});
    v = zeros(1, 2);
    v(1) = bitand(rv, bitMask1);
    v(2) = bitshift(bitand(rv, bitMask2), -8);
end

function v = peekRegister(varargin)
    import com.verasonics.hal.hardware.*
    register = varargin{end};
    switch class(register)
        case 'com.verasonics.hal.hardware.AcqReg'
            boardIndex = varargin{1};
            v = Memory.peekAcqReg(boardIndex, register);
        case 'com.verasonics.hal.hardware.BkpReg'
            v = Memory.peekBkpReg(register);
        case 'com.verasonics.hal.hardware.TpcReg'
            v = Memory.peekTpcReg(register);
        otherwise
            error('Unknown register type');
    end
end

function Sensors = appendDescription(Sensors, formatSpec)
    Sensors = [Sensors cell(size(Sensors, 1), 1)];
    for i=1:size(Sensors,1)
        Sensors{i,end} = sprintf(formatSpec, Sensors{i,1});
    end
end

function s = formatValueString(t, e, c)
    valueFormatString = '%3.1f C (Error Level %s C, Cooling Window [%s %s] C)';
    function fv = fval(v)
        if isnan(v)
            fv = '----';
        else
            fv = num2str(v, '%3.1f');
        end
    end
    s = sprintf(valueFormatString, t, fval(e), fval(c(1)), fval(c(2)));
end

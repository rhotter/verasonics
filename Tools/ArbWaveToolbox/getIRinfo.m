function [IRstruct,clip] = getIRinfo(cirFileCode,IRevalElemOffsetInd);
%getIRinfo:  get cir information from cir File Code
%[IRstruct,clip] = getIRinfo(cirFileCode,IRevalElemOffsetInd);
% Signature:
%  cirFileCode:
%       'L7-4' ;
%       'L22-14vX' ;
%       'L22-14vX:3005c' ;

%john flynn 9/14/2018

% IR should be sampled at 250MHz
switch cirFileCode
    case {1,'L7-4'} %cirFileCode --------------
        cir_file1 =  'cirest0701_code_acoustic-L74-23.mat';

        try
            cirS = load(cir_file1,'cirS');cirS = cirS.cirS;

        catch
            cirS = load(fullfile(datadir,cir_file1),'cirS');cirS = cirS.cirS;
        end
        clip=struct;  %
        clip.startInd = 1;
        clip.endInd  = length(cirS.CIR.gateSignalModel);
        clipOffset = clip;  %temp
        [cirS_NoOffset] = parseCIRstruct(cirS,clip,0);
        cirS_ElemOffset = cirS_NoOffset;

    case {21,...
            'L11-5v' , ...
            }
        clip=struct;
        IR_version = 'published';
        switch IR_version
            case {1,'deprecated'}
                error('deprecated')
                cir_file1 =  'ir_L115v_180006_c.mat';
                %some hardcoded h-file -specific parameters:
                %clip.startInd = 25; clip.endInd  = 110;
            case 2
                cir_file1 =  'ir_L115v_180008_stability_d.mat';
                %some hardcoded h-file -specific parameters:
                clip.startInd = 25; clip.endInd  = 110;
            case {'published'}
                cir_file1= 'L11-5v_info';

                %    clip.startInd = 6; clip.endInd  = 80;%legacy info file
                clip.startInd = 1; clip.endInd  = []; %published

            otherwise
                error('bad switch')
        end
        clipOffset = clip;  %temp

        %get the IR file ..............
        if isempty(strfind(cir_file1,'.mat'))
            cir_file1=[cir_file1,'.mat'];
        end
        try
            cirS = load(cir_file1);
        catch
            lasterr
            cirS = load(fullfile(datadir,cir_file1));
        end
        %parse the IR file ..................
        try
            error('debug - old format, needs updating.')
            %legacy format noted here:
            cirS = cirS.cirS;
        catch
            [cirS_NoOffset] = parseCIRstruct(cirS,clip,0);
            [cirS_ElemOffset] = parseCIRstruct(cirS,clipOffset,IRevalElemOffsetInd);
        end
        %overwrite
        cirS = cirS_NoOffset;

    case {31  ...
         ... PROMOTED TO METHOD 32:   ,'L22-14vX', ...
            ...'L22-14vX:3005c' ...
            'L22-14vX:3005d' ...
            ...
            } %cirFileCode --------------
warning([mfilename,': Using legacy Impulse Response file - performance may be compromised.'])
warning([mfilename,': To use released IR file, specify "L22-14vX" for "xdcr" in encode_dp2.m .'])
        clip=struct;
        %    cir_file1 =  'ir_L2214vx_3005_c.mat';
        cir_file1 =  'ir_L2214vx_3005_d.mat';
        %some hardcoded h-file -specific parameters:
        clip.startInd = 25; clip.endInd  = 110;
        clipOffset = clip;  %temp

        %get the IR file ..............
        try
            cirS = load(cir_file1);
        catch
            lasterr
            cirS = load(fullfile(datadir,cir_file1));
        end
        %parse the IR file ..................
        try
            error('debug - old format, needs updating.')
            cirS = cirS.cirS;
        catch
            try
                [cirS_NoOffset] = parseCIRstruct(cirS,clip,0);
                [cirS_ElemOffset] = parseCIRstruct(cirS,clipOffset,IRevalElemOffsetInd);
            catch
                error('could not parse lR structure 5788593894587')
            end
        end
        %overwrite
        cirS = cirS_NoOffset;
   case {32  ...
            ,'L22-14vX' ...
            } %cirFileCode --------------
          clip=struct;
        IR_version = 'published';
        switch IR_version
                case {'published'}
                cir_file1= 'L22-14vX_info';

                 clip.startInd = 1; clip.endInd  = []; %published

            otherwise
                error('bad switch')
        end
        clipOffset = clip;  %temp

        %get the IR file ..............
        try
            cirS = load(cir_file1);
        catch
            lasterr
            cirS = load(fullfile(datadir,cir_file1));
        end
        %parse the IR file ..................
        try
            error('debug - old format, needs updating.')
            cirS = cirS.cirS;
        catch
            try
                [cirS_NoOffset] = parseCIRstruct(cirS,clip,0);
                [cirS_ElemOffset] = parseCIRstruct(cirS,clipOffset,IRevalElemOffsetInd);
            catch
                error('could not parse lR structure 5788593894587')
            end
        end
        %overwrite
        cirS = cirS_NoOffset;

    otherwise %cirFileCode --------------
        error('bad switch')
end %cirFileCode --------------


IRstruct=struct;
IRstruct.cir_file = cir_file1;
IRstruct.clip = clip;
IRstruct.clipOffset = clipOffset;
IRstruct.IR = cirS;
IRstruct.cirFileCode = cirFileCode;
IRstruct.IR_ElemOffset = cirS_ElemOffset;

% IRstruct. = ;

end   %getIRinfo %%%%%%%%%%%%%%%%%%%%%%%%%%
%MAIN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cirS]=parseCIRstruct(cirS,clip, ElemOffset)

if isfield(cirS,'info');
    cirFileFormat = '2019';
else
    cirFileFormat = 'legacy';
end

try  %level1
    switch lower(cirFileFormat)
        case {'2019','2019.a'}%cirFileFormat switch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            TXclock= ...
                cirS.info.impulseResponse.txClock;
            FsRX = TXclock;
            % txChannInfoString = 'error(''TX channel not defined'')';
            txChannInfoString = 'TX channel not defined';
            txChann = [];

            if ElemOffset ~= 0,
                error( ' non-zero ElemOffset not supported ')
            end

            cirEstModel = 'gateSignalModel';

            h250w = ...
                cirS.info.impulseResponse.impulseResponseTXRate0dBGain;
            CIR = struct(cirEstModel,[]);
            CIR.(cirEstModel)=h250w;
            cirS.CIR =CIR ;

        case {'legacy'} %cirFileFormat switch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            TXclock = 250;

            try
                txChannInfoString= 'cirS.Params.txEl';

                txChann = eval(txChannInfoString);
            catch
                try
                    txChannInfoString= 'cirS.Params.TX(2).Apod';
                    disp([mfilename,':info[5472487]:txChannInfoString=',txChannInfoString])
                    txChann= find(eval(txChannInfoString));
                catch
                    disp([mfilename,'info(578475307): "txChann" not defined for this CIR file.'])
                    txChann = [];
                end
            end
            try
                FsRX = ...
                    cirS.Params.Receive(2).demodFrequency*...
                    cirS.Params.Receive(2).samplesPerWave;
            catch
                FsRX =[];
            end

        otherwise   %cirFileFormat switch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            lasterr
            error('unknown cir file format 5875874783975 ')

    end %cirFileFormat switch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if isempty(FsRX)
        disp([mfilename,'info(45233535): "FsRX" not defined for this CIR file: using TxClock.'])
        FsRX = TXclock;
    end

    timing2 = struct('RXfreq',FsRX,'TXclock',TXclock,'Flpf',[]);
    cirS.timing = timing2;
    [cirS_cond] = getConditionedIR(cirS,txChann+ElemOffset,clip);
    cirS = cirS_cond;
    cirS.ElemOffset = ElemOffset;


catch %level1
    lasterr
    warning([mfilename,': failed CIR parse.'])
    whos cirS
    cirS = [];

end %level1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [condH] = getConditionedIR(cirS,txChann,clip);

        CIREstModel = 'gateSignalModel';

        if isfield(cirS,'Hall')
            hRX =double(cirS.Hall(:,txChann));
        else
            hRX =cirS.CIR.(CIREstModel);
        end
        if ~isempty(clip)
            startInd=clip.startInd;
            endInd=clip.endInd;
        end
        if isempty(startInd)
            startInd = 1;
        end
        if isempty(endInd)
            endInd=length(hRX);
        end
        hRX = hRX(startInd:endInd);

        %ensure h sample rate is correct:
        timing = cirS.timing;
        FsRX = timing.RXfreq;
        TXclock = timing.TXclock;

        M= makeIntFromFraction(FsRX);
        %should be really close to ints:
        P = round(TXclock*M);Q = round(FsRX*M);
        h250 = resample(hRX,P , Q);
        cirS_temp = struct;
        cirS_temp.CIR.(CIREstModel)=h250;
        cirS_temp.timing = timing;
        cirS_temp.creator = mfilename;

        condH =cirS_temp;
    end %getConditionedIR

end %parseCIRstruct %%%%%%%%%%%%%%%%%


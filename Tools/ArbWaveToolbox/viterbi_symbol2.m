%viterbi_symbol.m
%Define symbol period processed by Viterbi algorithm.
% Called by viterbi_factor.m

 disp([mfilename,': symbolRateSetup: ',num2str(symbolRateSetup)])

 switch symbolRateSetup ,
    case {1 } %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=4;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        FcFactored = 9;Wbeta=3;
        shortPulseFunc=@fir1;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta)};
     case {2} %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=4;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        FcFactored = 7;Wbeta=2; %-30.3
        shortPulseFunc=@fir1;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta)};
    case {3} %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=4;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        FcFactored = 7;Wbeta=2.1; %
        shortPulseFunc=@fir1;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta)};
    case {4} %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=4;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        FcFactored = 7;Wbeta=2.5; %
        shortPulseFunc=@fir1;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta)};
    case {5} %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=3;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        FcFactored = 7;Wbeta=3; %
        shortPulseFunc=@fir1;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta)};
    case {6} %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=3;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        firRegFact = 1.05;
        FcFactored = 7;Wbeta=3; %
        shortPulseFunc=@fir1reg;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};
    case {7} %symbolRateSetup %** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=3;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        firRegFact = 1.05;
        FcFactored = 3;Wbeta=3; %
        shortPulseFunc=@fir1reg;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};
    case {8} %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=3;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        firRegFact = 1.05;
        FcFactored = 5;Wbeta=3; %
        shortPulseFunc=@fir1reg;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};
    case {9} %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=3;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        firRegFact = 1.05;
        FcFactored = 12;Wbeta=3; %
        shortPulseFunc=@fir1reg;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};
      case {10} %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=4;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        FcFactored = 10.1;Wbeta=3; %
        shortPulseFunc=@fir1;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta)};
    case {31} %symbolRateSetup %** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=3;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        firRegFact = 1.05;
     %   FcFactored = [17 ];Wbeta=3; %
        FcFactored = [23 ];Wbeta=3; %
        shortPulseFunc=@fir1reg;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};
  case {32} %symbolRateSetup %** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=3;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        firRegFact = 1.05;
        FcFactored = [9 18 ];Wbeta=3; %
        shortPulseFunc=@fir1reg;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};
     case {33} %symbolRateSetup %** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=3;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        firRegFact = 1.05;
       FcFactored = [15.625 ];Wbeta=3; %
         shortPulseFunc=@fir1reg;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};
    case {34} %symbolRateSetup %** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpp=3;  %clocks per pulse
        Fv=250/cpp;  %viterbi symbol rate/ frequency
        firRegFact = 1.05;
        FcFactored = [26 ];Wbeta=3; %
        shortPulseFunc=@fir1reg;
        shortPulseFuncParamList=...
            {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};
     case {35} %symbolRateSetup %** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         cpp=3;  %clocks per pulse
         Fv=250/cpp;  %viterbi symbol rate/ frequency
         firRegFact = 1.05;
         FcFactored = [20 ];Wbeta=3; %
         shortPulseFunc=@fir1reg;
         shortPulseFuncParamList=...
             {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};

     case {41,'L11-5v.std'} %symbolRateSetup %** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         cpp=3;  %clocks per pulse
         Fv=250/cpp;  %viterbi symbol rate/ frequency
         firRegFact = 1.05;
         FcFactored = 7.8125*[ 0.5 1.0 ];Wbeta=3; %
    %     FcFactored = 7.8125*[ 0.5 1.0 ];Wbeta=3; %
         shortPulseFunc=@fir1reg;
         shortPulseFuncParamList=...
             {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};

     case {42,'L11-5v.hi','L11=5v.hirate'} %symbolRateSetup %** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         cpp=2;  %clocks per pulse
         Fv=250/cpp;  %viterbi symbol rate/ frequency
         firRegFact = 1.05;
         FcFactored = 7.8125*[ 0.5 1.0 ];Wbeta=3; %
    %     FcFactored = 7.8125*[ 0.5 1.0 ];Wbeta=3; %
         shortPulseFunc=@fir1reg;
         shortPulseFuncParamList=...
             {Nfir,[FcFactored]*2/Fv,kaiser(Nfir+1,Wbeta),firRegFact};

     otherwise %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('bad switch')
end %symbolRateSetup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFF = {'cpp' ...
    'Fv' ...
    'FcFactored' ...
    'shortPulseFunc' ...
    ... removing, should be Nfir (set by Nalphabet)  'NfactPulse' ...
    };
paramStruct.cpp=[];
NNN = length(FFF);
for kp=1:NNN,
    commstr=[   'paramStruct.(FFF{',num2str(kp),'})=',  FFF{kp} ,';'];
    disp(commstr)
    eval(commstr)
end

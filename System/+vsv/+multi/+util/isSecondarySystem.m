function isSecondary = isSecondarySystem()
% The multi-system secondary nodes are headless.
% See if this host is headless
%
% Copyright (C) 2001-2023, Verasonics, Inc.  All worldwide rights and
% remedies under all intellectual property laws and industrial property
% laws are reserved.

headlessStr = java.lang.System.getProperty('java.awt.headless');
isHeadless = ~isempty(headlessStr) && headlessStr.equals("true");
isSecondaryIP = false;

try
    % Get all of the IP address for this node.
    myIP = java.net.InetAddress.getLocalHost();
    allIPs = java.net.InetAddress.getAllByName(myIP.getCanonicalHostName());

    numElements = allIPs.length;
    for i = 1 : numElements
        hostAddress = allIPs(i).getHostAddress();

        % The multi-system secondary nodes have an IP address pattern of: 10.10.*.2
        isSecStartIP = hostAddress.startsWith("10.10.");
        isSecEndIP = hostAddress.endsWith("2");

        isSecondaryIP = isSecStartIP && isSecEndIP;
        if isSecondaryIP
            break;
        end
    end
catch
    warning('vsv:activate:networkProblem', ...
        ['Cannot determine if this is a secondary system because ' ...
        'of a failure to get host IP address.']);
end

isSecondary = isHeadless && isSecondaryIP;
end

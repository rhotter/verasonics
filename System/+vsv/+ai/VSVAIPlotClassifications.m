classdef (Sealed) VSVAIPlotClassifications < handle
% VSVAIPLOTCLASSIFICATIONS Plots users defined or classifier defined
% classifications
%
% @type vsv.ai.VSVAIPlotClassification
%
% Version 1.0 | 2019-10-24
% $Author: Ellyn Cashdollar & Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%

    methods(Static)

        function plotHandle = plotBoundaries(ax, bn, color)
        % plotBoundaries Plots classification boundaries on image
        %
        % @param ax @type Axes
        % @param bn @type cell
        % @param color @type char


            hold(ax, 'on');

            nBN = length(bn);
            plotHandle = cell(1, nBN );

            for k = 1:nBN
                bb            = bn{k};
                plotHandle{k} = plot(ax, bb(:, 2), bb(:, 1), color, 'LineWidth', 2 );
            end
            hold( ax, 'off');
            drawnow;
        end

    end
end


classdef Logger
%LOGGER Enable/disable detailed event logging from shared libraries.
%
% @untested

    methods (Static)

        function start(logSeverity, logCategories, logFilename, logFolderPath)
            % Turn on detailed log reporting.
            %
            % @param logSeverity (optional) @type char.
            %        'i' Informational (default)
            %        'w' Warnings
            %        'e' Errors
            % @param logCategories (optional) @type char.
            %        'all' Show all log entries. (default)
            %        'dma' Show DMA log entries.
            %        'all_except_dma' Show all except DMA. (default)
            %        'hal' Show HAL log entries.
            %        'functions' Show function entry/exit log entries.
            % @param logFilename (optional) @type char.
            % @param logFolderPath (optional) @type char.

            import com.verasonics.common.logger.*
            import com.verasonics.hal.logger.*
            import com.verasonics.vantage.logger.*

            % Check input parameters.
            if nargin > 4
                error('Logger:start:invalidArgs', 'Expected 0-4 input arguments');
            end

            % Set the user specified logging level.
            level = LoggingLevel.info;
            if(nargin > 0)
                if strcmp(logSeverity, 'i')
                    level = LoggingLevel.info;
                elseif strcmp(logSeverity,'w')
                    level = LoggingLevel.warn;
                elseif strcmp(logSeverity,'e')
                    level = LoggingLevel.error;
                else
                    error('Logger:start:invalidArgs', 'Expected values are "i", "w", "e"');
                end
            end

            % Set the user specified category.
            categories = (Logger.CATEGORY_EVERYTHING - HalLogger.CATEGORY_HAL_DMA);
            if (nargin > 1)
                if strcmp(logCategories,'all')
                    categories = Logger.CATEGORY_EVERYTHING;
                elseif strcmp(logCategories,'dma')
                    categories = (HalLogger.CATEGORY_HAL_DMA + HalLogger.CATEGORY_HAL_INTERRUPTS);
                elseif strcmp(logCategories, 'all_except_dma')
                    % leave as default value
                elseif strcmp(logCategories,'hal')
                    categories = (Logger.CATEGORY_COMMON_ALL + HalLogger.CATEGORY_HAL_ALL);
                elseif strcmp(logCategories,'functions')
                    categories = Logger.CATEGORY_FUNCTIONS;
                elseif strcmp(logCategories,'rdma')
                    categories = (VantageLogger.CATEGORY_VANTAGE_RDMA);
                else
                    error('Logger:start:invalidArgs', 'Unsupported category specified.');
                end
            end

            % Set the user specified log filename.
            logName = sprintf('VSV_%i%02i%02i_%02i%02i%02i', fix(clock));
            if(nargin > 2)
                logName = logFilename;
            end

            logDir = '.';
            if nargin > 3
                logDir = logFolderPath;
            end

            % Configure and enable logging.
            if ~Logger.configureFileLogging(logDir, logName, 'log', 6) || ...
                    ~Logger.setFileLogging(true) || ...
                    ~Logger.setLogging(true) || ...
                    ~Logger.setLoggingLevel(level) || ...
                    ~Logger.setLoggerShowingSourceCodeInfo(true) || ...
                    ~Logger.setLoggerShowingTimestamp(true) || ...
                    ~Logger.setLoggingCategories(categories)
                error('Failed to configure logger.')
            end
            Logger.logToConsoleForThisProcess(false);
        end

        function stop()
            % Turn off logging and close the log file.
            import com.verasonics.common.logger.*

            % Turn off logging and close the logger.
            Logger.setLogging(false);
            Logger.setFileLogging(false);
            Logger.closeLogger();
        end

        function result = getNumWarnings()
            % Get the number of warnings and/or errors that were logged.
            import com.verasonics.common.logger.*
            result = Logger.getNumWarnings();
        end

        function printWarnings()
            % Print the warnings and/or errors that were logged.
            import com.verasonics.common.logger.*
            Logger.printCachedWarnings();
        end

    end
end


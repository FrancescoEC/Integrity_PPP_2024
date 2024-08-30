function [gps_alm_file, glo_alm_file] = find_alm(GPS_week_start)

% [gps_alm_file, glo_alm_file] = find_alm(GPS_week_start);
%
% Search the Matlab path for the most recent GPS (and GLONASS) almanacs.
% The seach starts with the GPS_week_start and proceeds backward in time.
%
% Inputs:
%   GPS_week_start - GPS week number to begin almanac seach.  (optional)
%                     If no start week is given, the computer clock time will be used.
%                     The files that are searched for have the following
%                     naming convention
%                           gps###.alm or yuma###.txt for GPS Almanacs
%                           glo###.alm for GLONASS almanacs.
%                     Almanac names are stored as mod(1024) weeks to conform to YUMA standards.
%                     You may input the GPS_week with or without the 1024
%                     rollover weeks. FIND_ALM will automatically mod(1024)
%                     when searching for matching file names.
%                     This function will search backward across the 1024 rollover
%                     to find the most recent almanac.    See UTC2GPS for more
%                     information on handeling the 1024 week rollover times.  This
%                     function will search back for 250 weeks.  If there are no
%                     almanacs within the 250 weeks previous to the start date,
%                     blank output file name(s) are returned.
% Outputs:
%   gps_alm_file   - Most recent GPS almanac file name (string).
%   glo_alm_file   - Most recent GLONASS almanac file name (string).
%
% See also READYUMA

% Written by: Jimmy LaMance 8/20/97
% Modified by: Maria Eagen for new GPS week rollover handling 10/19/2002
% Copyright (c) 1998 by Constell, Inc.

% functions called: ERR_CHK

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Check inputs to determine if a valid week number was provided, if not
% get the time from the computer clock
if nargin < 1
    now_utc = clock;
    now_gps = utc2gps(now_utc);
    GPS_week_start = now_gps(1);
end % if nargin == 1

% Take care of week 1024 rollovers in the GPS start week
max_lookback = 250; 			% maximum number of weeks to look backward for an almanac

GPS_week_start = rem(GPS_week_start,1024);

% Get the current Matlab version
% Get the current Matlab version
matlab_version = version;
matlab_version = str2num(matlab_version(1));

% If the Matlab version is 5.x and the DEBUG_MODE flag is not set
% then set up the error checking structure and call the error routine.
if matlab_version >= 5.0
    estruct.func_name = 'FIND_ALM';

    % Develop the error checking structure with required dimension, matching
    % dimension flags, and input dimensions.
    estruct.variable(1).name = 'GPS_week_start';
    estruct.variable(1).req_dim = [1 1];
    estruct.variable(1).var = GPS_week_start;

    % Call the error checking function
    stop_flag = err_chk(estruct);

    if stop_flag == 1
        fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
            estruct.func_name);
        return
    end % if stop_flag == 1
end % if matlab_version >= 5.0 & isempty(DEBUG_MODE)

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

GPS_file_found = 0;
GLONASS_file_found = 0;

% set the output files to the null string
gps_alm_file = blanks(0);
glo_alm_file = blanks(0);

% set the starting week
this_week = GPS_week_start;

% start looping through GPS weeks for the most recent GPS almanac
num_looks = 0;
while GPS_file_found == 0 & this_week >= 0 & num_looks < max_lookback
    % set up the file name for the GPS almanac using the user provided week
    %Fixed handeling of input start GPS weeks < 1024.  Modified file search string to force a 1, 2, 3, and 4 digit search.
    for i=1:4
        format_string = sprintf('gps%%0%dd.alm',i);
        testfile = sprintf(format_string,this_week);

        % see if the GPS almanac exists
        if exist(testfile) == 2    % the file exists
            GPS_file_found = 1;
            gps_alm_file = testfile;
            break;
        end % if exist(testfile) == 2
    end % for

    for i=1:4
        format_string = sprintf('yuma%%0%dd.txt',i);
        testfile = sprintf(format_string,this_week);

        % see if the GPS almanac exists
        if exist(testfile) == 2    % the file exists
            GPS_file_found = 1;
            gps_alm_file = testfile;
            break;
        end % if exist(testfile) == 2
    end % for

    % decrease the week number to check
    if this_week == 0 & rollover_flag == 1
        this_week = 1024;
    end % if this_week == 0 & rollover_flag == 1

    this_week = this_week - 1;
    num_looks = num_looks + 1;
end % while GPS_file_found == 0 & this_week >= 0

% repeat the process for the GLONASS almanac
% set the starting week
this_week = GPS_week_start;

% start looping through GPS weeks for the most recent GLONASS almanac
num_looks = 0;
while GLONASS_file_found == 0 & this_week >= 0 & num_looks < max_lookback
    for i=1:4
        format_string = sprintf('glo%%0%dd.alm',i);
        testfile = sprintf(format_string,this_week);

        % see if the GLONASS almanac exists
        if exist(testfile) == 2    % the file exists
            GLONASS_file_found = 1;
            glo_alm_file = testfile;
            break;
        end % if exist(testfile) == 2
    end % for

    % decrease the week number to check
    if this_week == 0 & rollover_flag == 1
        this_week = 1024;
    end % if this_week == 0 & rollover_flag == 1

    this_week = this_week - 1;
    num_looks = num_looks + 1;
end % while GPS_file_found == 0 & this_week >= 0

% verify that files were found, if not report a warning to the user
if GLONASS_file_found == 0 & GPS_file_found == 0
    fprintf('No GPS or GLONASS almanacs were found in the Matlab path.\n');
    fprintf('File names will return as a null string.\n');
    return
end % if GLONASS_file_found == 0 & GPS_file_found == 0

%if GLONASS_file_found == 0
%    fprintf('No GLONASS almanacs were found in the Matlab path.\n');
%    fprintf('The GLONASS file name will return as a null string.\n');
%    return
%end % if GLONASS_file_found == 0

if GPS_file_found == 0
    fprintf('No GPS almanacs were found in the Matlab path.\n');
    fprintf('The GPS file name will return as a null string.\n');
    return
end % if GPS_file_found == 0

%%%%% END ALGORITHM CODE %%%%%

% end of FIND_ALM

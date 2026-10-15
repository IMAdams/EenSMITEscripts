function T = extract_oex_cycle_time(rootDir, outCsv)
% EXTRACT_OEX_CYCLE_TIME  Batch-extract "cycle time" from cellSens .oex files.
%
%   T = extract_oex_cycle_time()              % prompts for a folder
%   T = extract_oex_cycle_time(rootDir)       % searches rootDir recursively
%   T = extract_oex_cycle_time(rootDir, csv)  % saves to a custom csv path
%   T = extract_oex_cycle_time(rootDir, false)% skip saving
%
% By default the table is saved as <rootDir>/cycle_times.csv.
%
% Looks for:
%   <attribute name="cycle time">
%       <double val="49.954..."></double>
%
% Output table columns:
%   File        - file name
%   Folder      - containing folder
%   NMatches    - number of cycle time entries found (check files where > 1)
%   CycleTime   - first match (assumed ms; verify against .vsi timestamps)
%   FPS         - 1000 / CycleTime
%   AllValues   - all matches, semicolon-separated

    if nargin < 1 || isempty(rootDir)
        rootDir = uigetdir(pwd, 'Select folder containing .oex files');
        if isequal(rootDir, 0), T = table(); return; end
    end
    if nargin < 2 || isempty(outCsv)
        outCsv = fullfile(rootDir, 'cycle_times.csv');
    elseif islogical(outCsv) || isnumeric(outCsv)
        if outCsv
            outCsv = fullfile(rootDir, 'cycle_times.csv');
        else
            outCsv = '';   % false/0 = don't save
        end
    end

    files = dir(fullfile(rootDir, '**', '*.oex'));
    nFiles = numel(files);
    if nFiles == 0
        warning('No .oex files found under %s', rootDir);
        T = table();
        return;
    end

    % Tolerant of whitespace/case differences and attribute spacing
    pat = '<attribute\s+name="cycle time"\s*>\s*<double\s+val="([^"]+)"';

    File      = strings(nFiles, 1);
    Folder    = strings(nFiles, 1);
    NMatches  = zeros(nFiles, 1);
    CycleTime = nan(nFiles, 1);
    FPS       = nan(nFiles, 1);
    AllValues = strings(nFiles, 1);

    for i = 1:nFiles
        fpath = fullfile(files(i).folder, files(i).name);
        File(i)   = string(files(i).name);
        Folder(i) = string(files(i).folder);

        txt = read_text_any_encoding(fpath);
        tok = regexp(txt, pat, 'tokens', 'ignorecase');
        vals = cellfun(@(c) str2double(c{1}), tok);

        NMatches(i) = numel(vals);
        if ~isempty(vals)
            CycleTime(i) = vals(1);
            FPS(i)       = 1000 / vals(1);   % assumes cycle time is in ms
            AllValues(i) = strjoin(string(vals), ';');
        end
    end

    T = table(File, Folder, NMatches, CycleTime, FPS, AllValues);

    % Flag anything unexpected
    nNone  = sum(NMatches == 0);
    nMulti = sum(NMatches > 1);
    if nNone > 0
        warning('%d file(s) had no cycle time match (different structure or encoding).', nNone);
    end
    if nMulti > 0
        warning('%d file(s) had multiple cycle time entries; check AllValues.', nMulti);
    end

    if ~isempty(outCsv)
        writetable(T, outCsv);
        fprintf('Wrote %s\n', outCsv);
    end
end


function txt = read_text_any_encoding(fpath)
% Read a file as text, handling UTF-8 and UTF-16 (with or without BOM).
    fid = fopen(fpath, 'r');
    if fid < 0
        error('Could not open %s', fpath);
    end
    cleanup = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, '*uint8')';

    if numel(bytes) >= 2 && isequal(bytes(1:2), uint8([255 254]))
        txt = native2unicode(bytes(3:end), 'UTF-16LE');
    elseif numel(bytes) >= 2 && isequal(bytes(1:2), uint8([254 255]))
        txt = native2unicode(bytes(3:end), 'UTF-16BE');
    elseif ~isempty(bytes) && mean(bytes == 0) > 0.3
        % Looks like BOM-less UTF-16LE (lots of null bytes)
        txt = native2unicode(bytes, 'UTF-16LE');
    else
        txt = native2unicode(bytes, 'UTF-8');
    end
end
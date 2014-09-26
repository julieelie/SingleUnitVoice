function [calfilename] = AddVoicenames( MatfilePath )
%Open the voice random file to add the names of individuals stim files of
%the original matrix (needed for MI calculations)
%   Detailed explanation goes here
load(MatfilePath,'VoiceTypeSel');
load(MatfilePath,'subject');
[Path, Matfile] = fileparts(MatfilePath);

if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                calfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',subject,['RandPVoi_' Matfile(9:end) '.mat']);
            end
        elseif strcmp(strtok(username), 'elie')
            calfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',['RandPVoi_' Matfile(9:end) '.mat']);
        end
else
    calfilename=fullfile('/auto','k6','julie','matfile',subject,['RandPVoi_' Matfile(9:end) '.mat']);
end

save(calfilename, 'VoiceTypeSel', '-append');
fprintf(1,'VoiceTypeSel appended to %s\n', calfilename)
end


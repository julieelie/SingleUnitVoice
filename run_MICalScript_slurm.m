cd /auto/k6/julie/matfile
resultsDirectory='/auto/k6/julie/matfile';
addpath('/auto/k1/queued');

input_dir=pwd;
Subjects = dir(input_dir);
%cmd = 'Cal_MIperarea_ConfVoice(''%s'');';
cmd = 'Cal_MIperarea_RandVoice(''%s'');';

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
%     if strcmp(Indiv, 'GreBlu9508M')
        Idir=fullfile(input_dir, Indiv);
        %matfiles=dir(fullfile(Idir,'ConfVoi_*.mat'));
        matfiles=dir(fullfile(Idir,'RandPVoi*.mat'));
        lm=length(matfiles);
        SS_Ind=zeros(lm,1);
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
            %if ~isempty(strfind(matfiles(ff).name, 'e12'))
%                 SS_Ind(ff)=1;
%             end
%             if ~isempty(strfind(matfiles(ff).name, 'e21'))
%                 SS_Ind(ff)=1;
%             end
        end
        Indices=find(SS_Ind);
        LM=length(Indices);
        %Randlist=randperm(LM);
        for hh=1:LM
            MATName=matfiles(Indices(hh)).name;
            jobParams = struct;
            jobParams.partition = 'all';
            jobParams.cpus = 2;
            jobParams.memory = 500;
            jobParams.out = fullfile(resultsDirectory, Indiv, sprintf('slurmout_MICal_%s.txt', MATName));
            jobParams.err = jobParams.out;
            MatfilePath=fullfile(Idir, MATName);
            icmd = sprintf(cmd, MatfilePath); 
            fprintf('%s. Calling slurm_sbatch with command %s\n',Indiv, icmd);
            slurm_sbatch(icmd,jobParams);
            
        end
    end
end
exit

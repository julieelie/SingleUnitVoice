function [calfilename]=Voice_selectivity_ConfusionMatrix(MatfilePath)
Plotfig=0; %set this to 1 if you want to see confusion matrices 
Cal = struct();
Res=load(MatfilePath, 'VocType', 'Trials', 'Trials_BG','TDT_wavfiles', 'SectionLength', 'Site', 'subject','Original_wavfiles', 'ESex', 'Eage');

%% Extract Emitter ID
Ename = cell(length(Res.VocType),1);
Esex = cell(length(Res.VocType),1);
Eage = cell(length(Res.VocType),1);
Erelated = cell(length(Res.VocType),1);
for ii=1:length(Res.VocType)
    [Path,File,Ext] = fileparts(Res.Original_wavfiles{ii});
    Ename{ii} = File(1:10);
    Esex{ii} = File(12);
    Eage{ii} = File(13);
    Erelated{ii} = File(14);
end

%% Get rid of unwanted stims
StimsRestricted = 1:length(Res.VocType);
StimsRestricted = setdiff(StimsRestricted, find(strcmp(Res.VocType, 'mlnoise')));
StimsRestricted = setdiff(StimsRestricted, find(strcmp(Res.VocType, 'Wh')));
StimsRestricted = setdiff(StimsRestricted, find(strcmp(Ename, 'STRFxx0000')));
StimTypeR=unique(Res.VocType(StimsRestricted));
NamesR = unique(Ename(StimsRestricted));
SexR = unique(Esex(StimsRestricted));
AgeR = unique(Eage(StimsRestricted));% Here I'm not using the one in h5 files because of NaN for age and related and sex for some of the song stimuli. This is a default of GenStimProtocol.m 
RelatedR = unique(Erelated(StimsRestricted));
NAgeR = length(AgeR);
NSexR = length(SexR);
NstimTypeR = length(StimTypeR);

%% Reorganize the data according to ID, sex, age and relatedness of the emitter in order to obtain the right confusion matrix
Original_wavfiles=Res.Original_wavfiles;
SpikeTrains_temp=Res.Trials;
SpikeTrains=cell([size(SpikeTrains_temp,1) + size(Res.Trials_BG,1) 1]);
TDTWavfiles=Res.TDT_wavfiles;
Sections_len_temp = Res.SectionLength;
Sections_len = zeros([size(Sections_len_temp,1) + size(Res.Trials_BG,1) 1]);
VocTypeList=Res.VocType;
VoiceTypeSel=cell([size(VocTypeList,1)+ size(Res.Trials_BG,1) 1]);
ORIwav=cell([size(VocTypeList,1)+ size(Res.Trials_BG,1) 1]);
TDTwav=ORIwav;
VoiceIndices = Sections_len;

cc=0;

for aa=1:NAgeR
    ag = AgeR(aa);
    selectorInd=find(strcmp(Eage, ag));%select one age
    for ss = 1:NSexR
        sex = SexR(ss);%select one sex
        selectorIndsex = intersect(selectorInd, find(strcmp(Esex, sex)));
        Nnames = unique(Ename(selectorIndsex));
        for nn=1:length(Nnames)
            name = Nnames(nn);%select one name
            selectorIndname = intersect(selectorIndsex, find(strcmp(Ename,name)));
            NselectorInd=length(selectorIndname);
            for cc_temp=1:NselectorInd
                cc=cc+1;
                VoiceIndices(cc) = selectorIndname(cc_temp);
                SpikeTrains{cc}=SpikeTrains_temp{selectorIndname(cc_temp)};
                Sections_len(cc)=Sections_len_temp(selectorIndname(cc_temp));
                VoiceTypeSel{cc} = cell2mat(strcat(ag,sex,Erelated{selectorIndname(cc_temp)},'_',name));
                TDTwav{cc}=TDTWavfiles{selectorInd(cc_temp)};
                ORIwav{cc}=Original_wavfiles{selectorInd(cc_temp)};
            end
        end
    end
end


SpikeTrains_temp=Res.Trials_BG;
NBG=size(Res.Trials_BG,1);
Rand_BG=randperm(NBG);
for cc_temp=1:20
    cc=cc+1;
    SpikeTrains{cc}=SpikeTrains_temp{Rand_BG(cc_temp)};
    Sections_len(cc)=600;
    VoiceTypeSel{cc}='BG';
    TDTwav{cc}=TDTWavfiles{Rand_BG(cc_temp)};
end
SpikeTrains=SpikeTrains(1:cc);
Sections_len=Sections_len(1:cc);
VoiceTypeSel=VoiceTypeSel(1:cc);
TDTwav=TDTwav(1:cc);
nfiles=cc;
VoiceIndices=VoiceIndices(1:cc);
ORIwav=ORIwav(1:cc);


%% Set up parameters for all the Info calculations    
% Window sizes for calculations that depend on window length
winSize = [2 5 10 20 30 50 100 300 600];    % Window size for static Info calculation and for ideal observer
% winSize = [50];    % Window size for static Info calculation and for ideal observer
ns = length(winSize);

% Initialize all return values to zero
%***Ideal Observer Calculations (by Trial Template)***
percorrectB = zeros(1,ns);
mi_confusionB = zeros(1,ns);

percorrectVT = zeros(1,ns);
mi_confusionVT = zeros(1,ns);

percorrectB_BGwo = zeros(1,ns);
mi_confusionB_BGwo = zeros(1,ns);

percorrectVT_BGwo = zeros(1,ns);
mi_confusionVT_BGwo = zeros(1,ns);

zdistTB = zeros(1,ns);              % Uses std for other for both self and other in KL Divergence calc
pzdistTB = zeros(1,ns);
mizdistTB = zeros(1,ns);
mizdistncTB = zeros(1,ns);

zdistSB = zeros(1,ns);              % Recenters "other" distrubtions to zero to minimize additive variance across songs
pzdistSB = zeros(1,ns);             % and collapses distances across all trials of a single stimulus
mizdistSB = zeros(1,ns);
mizdistncSB = zeros(1,ns);
confusionMatrix = cell(1,ns);
confusionMatrixVT = cell(1,ns);
neventMatrix = zeros(1,ns);


%***Gamma Model Calculations***
% gamma_mutual_info = 0;
% gamma_noise_entropy = 0;
% gamma_spike_rate = 0;
% gamma_const = 0;
% rate_bandwidth = 0;
% rate_gamma = 0;
% fano_factor = 0;

%***Rate Information Calculations***
% rate_info_biased = zeros(1,ns);
% rate_info_bcorr = zeros(1,ns);
% rate_info_stderr = zeros(ns,2);

%% Calculations!!
if nfiles
    % With VR_distanceB
    for is=1:ns
        fprintf(1, 'Calculating confusion matrix with winSize = %d\n', winSize(is));
        % Reinsert this if using the confusion matrix in version B
        [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, confusion_matrix, nevent_matrix] = info_distanceB(nfiles, SpikeTrains , Sections_len, @VR_distanceB, winSize(is));

        %Calculating miconfusion without the BG category.
        selectorBG=find(strcmp(VoiceTypeSel, 'BG'));
        confusion_matrix_BGwo=confusion_matrix;
        confusion_matrix_BGwo(:,selectorBG)=[];
        confusion_matrix_BGwo(selectorBG,:)=[];
        confusion_matrix_BGwo = confusion_matrix_BGwo ./sum(sum(confusion_matrix_BGwo));
        mi_conf_BGwo=info_matrix(confusion_matrix_BGwo);
        pc_BGwo=sum(diag(confusion_matrix_BGwo));
        
        % creating the compiled confusion matrix per voice category
        VoiceTypeCM=unique(VoiceTypeSel);
        IndBG = find(strcmp(VoiceTypeCM, 'BG'));
        VoiceTypeCM = [VoiceTypeCM(1:(IndBG-1)); VoiceTypeCM((IndBG+1):end); VoiceTypeCM(IndBG)];%place BG at the end
        NVoiceTypeCM=length(VoiceTypeCM);
        confusion_matrix_VoiceType = zeros(NVoiceTypeCM, NVoiceTypeCM);
        for vtR=1:NVoiceTypeCM
            stR=VoiceTypeCM(vtR);
            selectorR=strcmp(VoiceTypeSel, stR);
            selectorIndR=find(selectorR);
            for vtC = 1:NVoiceTypeCM
                stC=VoiceTypeCM(vtC);
                selectorC=strcmp(VoiceTypeSel, stC);
                selectorIndC=find(selectorC);
                confusion_matrix_VoiceType(vtR,vtC)=sum(sum(confusion_matrix(selectorIndR, selectorIndC)));
            end
        end
        
        mi_confVT=info_matrix(confusion_matrix_VoiceType);
        pcVT=sum(diag(confusion_matrix_VoiceType));
        
        % compiled confusion matrix without BG category
        selectorBG=find(strcmp(VoiceTypeCM, 'BG'));
        confusion_matrix_VoiceType_BGwo=confusion_matrix_VoiceType;
        confusion_matrix_VoiceType_BGwo(selectorBG,:)=[];
        confusion_matrix_VoiceType_BGwo(:,selectorBG)=[];
        confusion_matrix_VoiceType_BGwo =  confusion_matrix_VoiceType_BGwo ./ sum(sum(confusion_matrix_VoiceType_BGwo));
        mi_confVT_BGwo = info_matrix(confusion_matrix_VoiceType_BGwo);
        pcVT_BGwo = sum(diag(confusion_matrix_VoiceType_BGwo));
        
        if Plotfig==1
            figure(is*10);
            imagesc(confusion_matrix);
            colorbar;
            xlabel('Model Vocalization');
            ylabel('Actual Vocalization');
            title(sprintf('All vocalizations\n Winsize=%d miconf=%0.2f PCC=%0.2f',winSize(is), mi_conf, pc));
            fprintf(1, 'Pcorrect = %f, MIconfusion = %f\n Song Calculation zd = %f pzd = %f MIzd = %f MIzd(nc) = %f\n Trial Calculation zd = %f pzd = %f MIzd = %f MIzd(nc) = %f\n', ...
            pc, mi_conf, zdS, pzdS, mi_zdS, mi_zdS_nc, zdT, pzdT, mi_zdT, mi_zdT_nc);
    
            figure(is*10+1);
            imagesc(confusion_matrix_VoiceType);
            colorbar;
            xlabel('Model vocalization');
            ylabel('Actual vocalization');
            title(sprintf('Vocalization Categories p(x,y)\n Winsize=%d miconf=%0.2f PCC=%0.2f', winSize(is), mi_confVT, pcVT));
            set(gca(), 'Ytick', 1:NVoiceTypeCM);
            set(gca(), 'YTickLabel', VoiceTypeCM);
            set(gca(), 'Xtick', 1:NVoiceTypeCM);
            set(gca(), 'XTickLabel', VoiceTypeCM);
        
             % calculating p(y/x) pour les articles!
            confusion_matrix_VT_p1=zeros(size(confusion_matrix_VoiceType));
            px=sum(confusion_matrix_VoiceType, 2);
            for cc=1:size(confusion_matrix_VoiceType,2)
                confusion_matrix_VT_p1(:,cc)=confusion_matrix_VoiceType(:,cc)./px;
            end
        
            figure(is*10+2);
            imagesc(confusion_matrix_VT_p1);
            colorbar;
            xlabel('Model vocalization');
            ylabel('Actual vocalization');
            title(sprintf('Vocalization Categories p(x/y)\n Winsize=%d', winSize(is)));
            set(gca(), 'Ytick', 1:NVoiceTypeCM);
            set(gca(), 'YTickLabel', VoiceTypeCM);
            set(gca(), 'Xtick', 1:NVoiceTypeCM);
            set(gca(), 'XTickLabel', VoiceTypeCM);
        end
        
        percorrectB(is) = pc;
        mi_confusionB(is) = mi_conf;
        percorrectVT(is) = pcVT;
        mi_confusionVT(is) = mi_confVT;
        percorrectB_BGwo(is) = pc_BGwo;
        mi_confusionB_BGwo(is) = mi_conf_BGwo;
        percorrectVT_BGwo(is) = pcVT_BGwo;
        mi_confusionVT_BGwo(is) = mi_confVT_BGwo;
        
        zdistTB(is) = zdT;
        pzdistTB(is) = pzdT;
        mizdistTB(is) = mi_zdT;
        mizdistncTB(is) = mi_zdT_nc;
        zdistSB(is) = zdS;
        pzdistSB(is) = pzdS;
        mizdistSB(is) = mi_zdS;
        mizdistncSB(is) = mi_zdS_nc;
        confusionMatrix{is} = confusion_matrix;
        confusionMatrixVT{is} = confusion_matrix_VoiceType;
        neventMatrix(is) = nevent_matrix;
    end
    

    % Gamma model information calculation
    %onset_time = 200;  % the first 200 ms of the response are removed from the data set
%    fprintf(1, 'Calculating Gamma model information\n');
%    onset_time = 0;
%    spiketrain = spike_times_to_train(NstimType, SpikeTrains, Sections_len, onset_time);
    %[info, noiseentropy, totalentropy, gamma_const, rate_bandwidth, rate_gamma, fano_factor]= gamma_info(spiketrain);

    %gamma_mutual_info = info(1)*1000;                    % Info rates in bits/second
    %gamma_noise_entropy = noiseentropy(1)*1000;       % Noiseentropy in bits/second
%    gamma_spike_rate = mean(mean(spiketrain))*1000;  % Spiking rate is spikes/second

    % Rate Information
    %fprintf(1, 'Calculating Rate Information\n');
    %[rate_info_biased, rate_info_bcorr, rate_info_stderr] = findInfoSR(spiketrain,winSize);

end

%% Store Values
fprintf(1, 'Storing values\n');
Cal.subject=Res.subject;
Cal.StimType=VoiceTypeCM;
Cal.winSize=winSize;
Cal.percorrectB = percorrectB;
Cal.mi_confusionB = mi_confusionB;
Cal.percorrectCT = percorrectVT;
Cal.mi_confusionCT = mi_confusionVT;
Cal.percorrectB_BGwo = percorrectB_BGwo;
Cal.mi_confusionB_BGwo = mi_confusionB_BGwo;
Cal.percorrectCT_BGwo = percorrectVT_BGwo;
Cal.mi_confusionCT_BGwo = mi_confusionVT_BGwo;
Cal.zdistSB = zdistSB;              % Recenters "other" distrubtions to zero to minimize additive variance across songs
Cal.pzdistSB = pzdistSB;             % and collapses distances across all trials of a single stimulus
Cal.mizdistSB = mizdistSB;
Cal.mizdistncSB = mizdistncSB;
%Cal.gamma_mutual_info = gamma_mutual_info;
%Cal.gamma_noise_entropy = gamma_noise_entropy;
%Cal.gamma_spike_rate = gamma_spike_rate;
%Cal.totalentropy=totalentropy;
%Cal.gamma_const=gamma_const;
%Cal.rate_bandwidth=rate_bandwidth;
%Cal.rate_gamma=rate_gamma;
%Cal.fano_factor=fano_factor;
%Cal.rate_info_biased=rate_info_biased;
%Cal.rate_info_bcorr=rate_info_bcorr;
%Cal.rate_info_stderr=rate_info_stderr;
% Cal.Best_nbPC=Best_nbPC;
% Cal.R2A=R2A;
% Cal.ModelPredict = ModelPredict;
% Cal.LogLikelihood = LL;
% Cal.Pvalue=Pvalue; %result of the anova on each model
% Cal.PValLRatio = PValLRatio;
% Cal.SignifModelCompare = h;
% Cal.NeuralResponse = NeuralResponse;
% Cal.VocType=voc;
Cal.confusionMatrix = confusionMatrix;
Cal.confusionMatrixCT = confusionMatrixVT;
Cal.neventMatrix=neventMatrix;
Cal.VoiceTypeSel=VoiceTypeSel;
Cal.TDTwav=TDTwav;
Cal.ORIwav=ORIwav;
Cal.VoiceIndices=VoiceIndices;




if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                calfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',Res.subject,['ConfVoi_' Res.Site '.mat']);
            end
        elseif strcmp(strtok(username), 'elie')
            calfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',['ConfVoi_' Res.Site '.mat']);
        end
else
    calfilename=fullfile('/auto','k6','julie','matfile',Res.subject,['ConfVoi_' Res.Site '.mat']);
end

save(calfilename, '-struct', 'Cal');
fprintf(1,'done making calculus on %s\nData save under %s\n',MatfilePath, calfilename);
clear Cal Res winSize percorrectB mi_confusionB zdistSB pzdistSB mizdistSB mizdistncSB confusionMatrix confusionMatrixCT
%clear gamma_mutual_info gamma_noise_entropy totalentropy gamma_const
%rate_bandwidth rate_gamma fano_factor rate_info_biased rate_info_bcorr rate_info_stderr
end 

 
%%%% LO SCRIPT CONVERTE E MERGIA TUTTI I FILE REGISTRATI NELLA STESSA
%%%% SESSIONE DI REGISTRAZIONE
%%%% VA INSERITO IL PERCORSO ALL'INTERNO DEL QUALE CI SONO I VARI SCRIPTS UTILI ALLA CONVERSIONE DAL FILE OPEN EPHYS


clear
close all
clc


%%%%%%%%%%  INPUT  %%%%%%%%%%%%

TASKS= {
    'Wifi_20240430_Clinico_TestaFissa';...
    'Wifi_20240430_Clinico_TestaLibera';...
    'Wifi_20240430_Gabbione'};

raw_files_directory='\\fs01.hpc.unipr.it\bonini01\WIRELESS_Wifi_Router_Impianti2021\Acquisizioni Deuteron\Wifi_20240430\';
Deuteron_file_type='DT6'; % pu� essere 'DT2' o 'DT4' o 'DT6' a seconda che siano 
                          % stati registrati 32 o 64 o 128 canali, vedere nella cartella "name_mrg"
start_and_stop_deuteron{1}=[6.8060 1117.1148]; % in secondi (presi dall'Excel "Start_and_Stop..."
start_and_stop_deuteron{2}=[3.6482 937.6556]; % in secondi
start_and_stop_deuteron{3}=[4.5617 1984.9776]; % in secondi

freq=32000;  % frequenza dei file per il merge, 32000Hz per il Deuteron
nCH=128;
% stepCH={'1:4', '5:8', '9:12', '13:16', '17:20', '21:24', '25:28', '29:32'}; 
stepCH={'1:32','33:64','65:96','97:128'};
gap=10; % 10s (ma non viene applicato nello script, serve solo per creare gli Start_and_Stop)
name_mrg=char('Wifi_20240430');   % Mettere il nome del file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% INIZIA LA CONVERSIONE DI OGNI SINGOLO TASK


% Creo e salvo gli Start and Stop per poi fare il reshift di 10s degli spikes 
% (dato che non posso farlo prima del sorting)
mkdir([raw_files_directory,'\File Merged no gap\Start and Stop'])
Start_and_Stop=zeros(size(TASKS,1),2);
for task=1:size(TASKS,1)-1
    Start_and_Stop(task+1,2)=start_and_stop_deuteron{task}(2)-start_and_stop_deuteron{task}(1) + Start_and_Stop(task,1);
    Start_and_Stop(task+1,1)=Start_and_Stop(task+1,2) + gap;
end
save([raw_files_directory,'\File Merged no gap\Start and Stop\Start_and_Stop_',name_mrg],'Start_and_Stop'); 



% Converto i task in .mat
mkdir([raw_files_directory,'\File singoli task'])
for task=1:size(TASKS,1)

    cd([raw_files_directory '\' TASKS{task}]); 
    files=dir(['*',Deuteron_file_type]);
    
    %%% A QUESTO PUNTO CONVERTIAMO IMPORTIAMO IN MATLAB I SINGOLI CANALI DEUTERON
    %%% OSS. CHE I SINGOLI FILE "NEURO000.DTx" CONTENGONO I DATI DI TUTTI I
    %%% CANALI IN UN CERTO BREVE FINESTRA TEMPORALE (es. 4s), NON � COME CON I
    %%% .continuous DELL'OPEN EPHYS IN CUI OGNI FILE � UN CANALE PER TUTTA
    %%% LA DURATA DELLA REGISTRAZIONE
    metaData=GetMetaData(Deuteron_file_type); % Estrae i parametri del logger a partire dall'estensione dei raw file
    switch Deuteron_file_type
        case 'DT2'
            time_duration_file=8.192; % seconds
        case 'DT4'
            time_duration_file=4.096; % seconds
        case 'DT6'
            time_duration_file=2.048; % seconds
    end
    
%     data_cell=cell(length(files),1);
    data=zeros(length(files)*freq*time_duration_file,nCH); 
    contatore=0;
    for t=1:length(files) % estraggo tutti i frammenti di ~4sec l'uno
        file_ID=fopen(files(t).name);
        all_channels_t=fread(file_ID,'uint16'); % each data point of neural data is a 16 bit word
        fclose(file_ID);
        all_channels_t=reshape(all_channels_t,metaData.numChannels,[]); % data are now in form of channels x samples
        all_channels_t_converted=metaData.voltageRes*(all_channels_t - 2^(metaData.numADCBits - 1)); % conversione in Volt, vedi "RatLog_instructions.pdf" 
        data((1:size(all_channels_t_converted,2))+contatore,:)=transpose(all_channels_t_converted);
        contatore=contatore+size(all_channels_t_converted,2);
    end
%     data=cell2mat(data_cell);    
    
    % Metto a zero i frammenti della traccia con dropped blocks, che
    % vengono messi dal Deuteron a -6.5mV di default
    data(data(:,1)<-6E-3,:)=0;  
    
    %%% Estraggo gli start e gli stop x eliminare tutto ci� che viene prima
    %%% dello Start e dopo lo Stop
    onset=round(freq*start_and_stop_deuteron{task}(1));
    offset=round(freq*start_and_stop_deuteron{task}(2));

    %%%%%%%%%%%%% Creo le variabili canali e taglio in corrispondenza di start e stop %%%%%%%%%%%%%%%%
    var=cell(1,nCH);
    for n=1:nCH
        if n<=10
            assignin('base',['CH','00',num2str(n-1)],data(onset:offset,n)); % "n-1" perch� il deuteron parte da 0 col conteggio canali
            var{1,n}=['CH','00',num2str(n-1)];
        elseif n>=11 && n<=100
            assignin('base',['CH','0',num2str(n-1)],data(onset:offset,n));
            var{1,n}=['CH','0',num2str(n-1)];
        elseif n>=101
            assignin('base',['CH',num2str(n-1)],data(onset:offset,n));
            var{1,n}=['CH',num2str(n-1)];
        end
    end
    cd(raw_files_directory);
    save(['File singoli task\',TASKS{task}],var{:},'onset','offset'); 
    clearvars -except Deuteron_file_type task stepCH name_mrg freq TASKS task nCH raw_files_directory start_and_stop_deuteron;
end

for blocco_canali=1:size(stepCH,2)
    nCH=str2num(stepCH{blocco_canali});            % RANGE DI CANALI CHE SI VOGLIONO MERGIARE
    savename_p=[name_mrg '_' stepCH{blocco_canali}];  % NOME DEL FILE DA SALVARE
    savename= strrep(savename_p, ':', '-');  % FREQUENZA DI CAMPIONAMENTO

    %%% rinominare i canali di ogni task t1 t2 t3...tn
    for n=1:size(TASKS,1)

        % Carico i canali specifici del blocco canali (in modo da
        % risparmiare tempo)
        cd('File singoli task\')
        list_of_file_variables=whos('-file',TASKS{n});      
        load(TASKS{n},list_of_file_variables(nCH).name);
        cd ..
        
        for i=nCH
            if i<=10
                assignin('base',['CH','00',num2str(i-1),'t',num2str(n)],eval(['CH','00',num2str(i-1)])); % "i-1" perch� il deuteron parte da 0 col conteggio canali
            elseif i>=11 && i<=100
                assignin('base',['CH','0',num2str(i-1),'t',num2str(n)],eval(['CH','0',num2str(i-1)]));
            elseif i>=101
                assignin('base',['CH',num2str(i-1),'t',num2str(n)],eval(['CH',num2str(i-1)]));
            end
        end
    end
    
    
    %---------- MERGIAMO I FILES-----------
    var=[];
    for i=nCH
        if i<=10
            assignin('base',['CH','00',num2str(i-1)],[]); % replace the CH*** variable with [] before assign the merged tasks (it may not be mandatory)
            var{end+1}=['CH','00',num2str(i-1)];
        elseif i>=11 && i<=100
            assignin('base',['CH','0',num2str(i-1)],[]);
            var{end+1}=['CH','0',num2str(i-1)];
        elseif i>=101
            assignin('base',['CH',num2str(i-1)],[]);
            var{end+1}=['CH',num2str(i-1)];
        end
        for n=1:size(TASKS,1)
            if i<=10
                assignin('base',['CH','00',num2str(i-1)],[eval(['CH','00',num2str(i-1)]);eval(['CH','00',num2str(i-1),'t',num2str(n)])]);
            elseif i>=11 && i<=100
                assignin('base',['CH','0',num2str(i-1)],[eval(['CH','0',num2str(i-1)]);eval(['CH','0',num2str(i-1),'t',num2str(n)])]);
            elseif i>=101
                assignin('base',['CH',num2str(i-1)],[eval(['CH',num2str(i-1)]);eval(['CH',num2str(i-1),'t',num2str(n)])]);
            end
        end
    end
    %%%%%% definiamo le variabili che verranno salvate nel file finale
    save(['File Merged no gap\',savename],var{:},'-v7.3');
    clearvars -except task TASKS stepCH name_mrg freq raw_files_directory;
end


clear
close all
clc
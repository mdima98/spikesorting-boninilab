%% ====================================== Example Script ====================================== %%
%    This script takes a single file of logger data and extracts it and converts it to neural data

%% ====================================== Initialization ====================================== %%
filePath = ''; % Enter the path for the location of your file
fileName = 'NEUR1490.DT6'; % Enter your filename


%% ========================================= Read data ========================================= %%

myFile = fullfile(filePath, fileName);
fid = fopen(myFile);
data = fread(fid, 'uint16'); % each data point of neural data is a 16 bit word
fclose(fid);
ext = fileName(end-2:end); % We know what type of file this is ba1.sed on the file extension ('DT2, 'DT4', DT8' or 'DAT')


%% ================================Prepare data for processing ================================ %%

metaData = GetMetaData(ext); % checks what kind of logger it is and sets parameters. This must be the original extension of the file.
dataMatrix = reshape(data', metaData.numChannels, []); % data are now in form of channels x samples

%% =============================== Conversion to neural data ==============================================%%

neuralData = metaData.voltageRes*(dataMatrix - 2^(metaData.numADCBits - 1)); % conversion of data to neural data

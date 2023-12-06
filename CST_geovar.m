
%Set Parameters that is going to be changed
params_str = {'Dy','h'};% param names
params = [600,1000];    % param nominal value
param_per = 5/100;      % percentage change

%Load LHS 
savePath ='C:\Users\cmpieterse\Documents\1_PhD_Work\Geo_var\';
filename ='Xs1.mat';
Xs = load([savePath filename]);
Xs = Xs.Xs;

%Get number of simulations
runs =size(Xs);
runs = runs(1);

%Calculate parameter values for LHS
new_params = (params-params*(param_per))+(params*(param_per)).*Xs(:,1);
new_params = round(new_params,2);

%Open CST and get simulation info
filePathName = 'C:\Users\cmpieterse\Documents\1_PhD_Work\Antenna_geovar_test.cst';
C = CSTsimulation(filePathName,false);
C = C.init;
C = C.getSimInfo;


% 
% C = C.updateParams(params_str,params);
% C = C.runSim;
% [S,Sn,freqS] = C.readSparResults;
% FF = C.readFarFieldResults;


%%
run_count =1;
for kk = 1:1:runs
    C = C.updateParams(params_str,new_params(kk,:));
    C = C.runSim;
    %Get results
    [S,Sn,freqS] = C.readSparResults;
    FF = C.readFarFieldResults;
    
    parseobj = inputParser;
    parseobj.FunctionName = 'readCSTffs';

    typeValidator_pathName = @(x) isa(x,'char');
    parseobj.addRequired('pathname',typeValidator_pathName);

    for ii = 1:C.Nports
        pathName = [C.filePath,'\',C.fileName,['\Result\Farfield Source [',num2str(ii),']']];
    end

    %Open the data file
    if ~strcmp(pathName(end-3:end),'.ffs')
        pathName = [pathName,'.ffs'];
    end

    fid = fopen(pathName);
    % Read the main header info
    freqMarker = '// #Frequencies';
    powerFreqMarker = '// Radiated/Accepted/Stimulated Power , Frequency';
    NphNthMarker = '// >> Total #phi samples, total #theta samples';
    fieldMarker = '// >> Phi, Theta, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi):';

    fCount = 0;
    read = 1;
    while read
        a = fgetl(fid);
        if strcmp(a,freqMarker) % Read the number of frequencies
            Nf = fscanf(fid,'%i',1);
        end

        if strncmp(a,powerFreqMarker,length(powerFreqMarker))
            PF = fscanf(fid,'%f',[Nf*4,1]);
            [Prad,Pacc,Pstim,freq] = deal(zeros(1,Nf));
            for ii = 1:Nf
                Prad(ii) = PF((ii-1)*4+1);
                Pacc(ii) = PF((ii-1)*4+2);
                Pstim(ii) = PF((ii-1)*4+3);
                freq(ii) = PF((ii-1)*4+4);
            end
        end
        if strncmp(a,NphNthMarker,length(NphNthMarker))
            NphNth = fscanf(fid,'%i %i',[2,1]);
            Nph = NphNth(1);
            Nth = NphNth(2);
        end

        if strncmp(a,fieldMarker,length(fieldMarker))
            fCount = fCount + 1;
            fDataForm = '%f %f %f %f %f %f';
            fData = fscanf(fid,fDataForm,[6,Nth*Nph])';
            Eth(:,fCount) = fData(:,3) + 1i.*fData(:,4);
            Eph(:,fCount) = fData(:,5) + 1i.*fData(:,6);
            if fCount == Nf
                read = 0;
            end
        end
    end

    fclose(fid);

    radEff = Prad./Pacc;

    %n =1;
    filename = ['GeoVar' '_' num2str(kk) '_' params_str{1} '_' num2str(new_params(kk,1)) '_' params_str{2} '_' num2str(new_params(kk,2)) '.txt'];
    savePath = 'C:\Users\cmpieterse\Documents\1_PhD_Work\Geo_var\';
    writecell(params_str,[savePath filename]);
    writematrix(new_params(kk,:),[savePath filename],'WriteMode','append');
    writematrix('E1',[savePath filename],'WriteMode','append');
    writematrix(Eth,[savePath filename],'WriteMode','append');
    writematrix('E2',[savePath filename],'WriteMode','append');
    writematrix(Eph,[savePath filename],'WriteMode','append');
    writematrix('Prad',[savePath filename],'WriteMode','append');
    writematrix(Prad,[savePath filename],'WriteMode','append');
    writematrix('radEff',[savePath filename],'WriteMode','append');
    writematrix(radEff,[savePath filename],'WriteMode','append');
    writematrix('S11',[savePath filename],'WriteMode','append');
    writematrix(S,[savePath filename],'WriteMode','append');
    writematrix('freqS',[savePath filename],'WriteMode','append');
    writematrix(freqS,[savePath filename],'WriteMode','append');
end

function [FFref] = getReferenceFarField(filePathName, params_str,params,savePath)
C = CSTsimulation(filePathName,false);
C = C.init;
C = C.getSimInfo;
C = C.updateParams(params_str,params);
C = C.runSim;
[S,Sn,freqS] = C.readSparResults;
FFref = C.readFarFieldResults;
save([savePath 'FF.mat'],'FF');
end


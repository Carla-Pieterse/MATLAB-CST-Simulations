%% Get Directivity Error for changing amount of CBFP coefficients
y=1;
D =(4.*pi.*FF.getU)./FF.pradInt;
index_abvHor = FF.th <= pi/2;
nmax =8; %23;
nmin =1;%16;
save = 0; %1=print to txt
for n =24%nmin:1:nmax %1:1:12%FF.Nf %Max number of basis functons = Nf
    %tic
    NumBasis = n;
    objCBFP = CBFP(FF,[],[],NumBasis);
    W = objCBFP.farField2Coeffs(FF,[],NumBasis);

    F_F= objCBFP.coeffs2FarField(W,FF.freq);
    % E-Field Pattern Error
    E_fine = [FF.E1;FF.E2];
    E_expn = [F_F.E1;F_F.E2];
    e =20.*log10(sqrt(1./FF.Nang*2.*sum(abs(E_fine-E_expn).^2)));
    e_Err(y,:) = e;

    % D-Field Pattern Error
    D_r =(4.*pi.*F_F.getU)./F_F.pradInt;

    d =10.*log10(sqrt(1./length(index_abvHor).*sum(abs(D(index_abvHor,:)-D_r(index_abvHor,:)).^2)));
    d_Err(y,:) =d;

    %Save out some directivities to txt
    if save==1 && n>=nmin &&n<=nmax
        filename = ['D_CBFP' '_' num2str(n) '.txt'];
        savePath = 'C:\Users\cmpieterse\Documents\1_PhD_Work\Reports\Progress Reports\2023_05_10 Threshold paper\'; % Put the string of the path you want to save your file in
        writematrix(D_r,[savePath filename],'WriteMode','overwrite')
    end

    y=y+1;
end

%% Write angle pairs to Matlab
% 
% angles = [FF.x FF.y];
% filename = ['phi_theta' '.txt'];
% savePath = 'C:\Users\cmpieterse\Documents\1_PhD_Work\Reports\Progress Reports\2023_05_10 Threshold paper\'; % Put the string of the path you want to save your file in
% writematrix(angles,[savePath filename],'WriteMode','overwrite')
% 
% 
%% Read Tant 
% 
% filename = ['Gal_down' '.txt'];
% savePath = 'C:\Users\cmpieterse\Documents\1_PhD_Work\Reports\Progress Reports\2023_05_10 Threshold paper\';
% Gal_down = dlmread([savePath filename]);
% 
% filename = ['Gal_up' '.txt'];
% Gal_up = dlmread([savePath filename]);
%% Error Plot
% 
% figure
% plot(1:1:NumBasis, max(d_Err'),'-b',LineWidth=2)
% hold on
% plot(1:1:NumBasis, mean(d_Err'),'-r',LineWidth=2)
% plot(1:1:NumBasis, median(d_Err'),'-k',LineWidth=2)
% plot(1:1:NumBasis, min(d_Err'),'-g',LineWidth=2)
% xlabel('Number of coefficients',FontSize=14)
% ylabel('Directivity Error')
% %ylabel('$10\log{\frac{|Tant_{o}-Tant_{r}|}{Tant_{o}}}$','interpreter','latex',FontSize=16)
% ax = gca;
% ax.FontSize = 16;
% legend('max','mean','median','min')

%% Observer
%
% Sky = Haslam('MHz',-2.5);
% Sky = Sky.setLocation([-30.71131, 21.4476236,1151]);
% EOR = [];
% antCoor = CoordinateSystem;
% includeGround =true;
% includeAtm = false;
% obj = Observer(FF,Sky,EOR,antCoor,includeGround,includeAtm);
% time = datetime(2019,7,1,0,0,0); %Gal Up
% underSampleFactor =1;%3;
% Tant = obj.getTant(time,underSampleFactor);
% figure
% Tant.plotInputMap
% figure
% Tant.plotTantFreq
% Tant_Gal_up = Tant.Tant;
% 
% time = datetime(2019,10,1,0,0,0); %Gal Down
% Tant = obj.getTant(time,underSampleFactor);
% figure
% Tant.plotInputMap
% figure
% Tant.plotTantFreq
% Tant_Gal_down = Tant.Tant;

%% Different type of errors
%
% e =20.*log10(mse(abs(E_fine(:,x)-E_expn(:,x)).^2));
%e =20.*log10(sqrt(mean(abs(E_fine-E_expn)).^2));
%e =20.*log10(sqrt(2./FF.Nang.*sum(abs(E_fine-E_expn).^2)));
%d =10.*log10(mse(D-D_r).^2);
%d= 10.*log10(sqrt(mean(D-D_r).^2));
%         d_mse(x)=10.*log10(sum(D-D_r)./sum(D.^2));
%         d_mse_Err(x,y)=d_mse(x);
%         d_dbi = 10.*log10(sum(FF.Directivity_dBi-F_F.Directivity_dBi)./FF.Nang);
%         d_dbi_error(x)=d_dbi;
% end
%d_dbi = 10.*log10(sum(FF.Directivity_dBi-F_F.Directivity_dBi)./FF.Nang);
%Error(x,y) = e(x);
%d(x)= 10.*log10(mean(abs(D_r(x,:)-D(x,:))));
%d_dbi = 10.*log10(sqrt(sum((FF.Directivity_dBi-F_F.Directivity_dBi).^2)./FF.Nang));

%% Coeff change for two different CST simulations
% 
% NumBasis=9; %14; 
% objCBFP = CBFP(FF,[],[],NumBasis);
% W = objCBFP.farField2Coeffs(FF,[],NumBasis);
% objCBFP_J = CBFP(FF_John,[],[],NumBasis);
% W_J = objCBFP.farField2Coeffs(FF_John,[],NumBasis);  
% figure, plot(FF.freq, real(W_J{1,1}(:,:)-W{1,1}(:,:)))
% figure, plot(FF.freq, real(W{1,1}(:,:)),'--')
% hold on
% plot(FF.freq, real(W_J{1,1}(:,:)))
% 
% obj = FarFieldReconstructor(FF,1,[],1,[]);
% obj_J = FarFieldReconstructor(FF_John,1,[],1,[]);
% 
% figure
% hold on
% for ii = 1:1:obj.NumBasis
%         plot(abs(obj.beta{1,ii}-obj_J.beta{1,ii}),'-o')
% end

%% 'CBFP' on directivity
% D = FF.getDirectivity;
% FM =D;
% [U,S,V] = svd(FM,'econ');
% Snorm = S(1,1);
% sigmaMat = diag(S)./Snorm;
% threshold = 1e-3;
% index = sigmaMat>1e-3;
% iB = sum(index);
% range_B = 1:1:iB;
% Umat = U(:,range_B);
% Smat = S(range_B,range_B);
% Vmat = V(:,range_B);
%                     
% % Calculate weights for all FarFields in the FM matrix
% Wmat = pinv(Umat)*FM;
% W = mat2cell(Wmat,NR,Nfs);
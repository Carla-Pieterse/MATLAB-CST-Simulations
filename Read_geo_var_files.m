savePath ='C:\Users\cmpieterse\Documents\1_PhD_Work\Geo_var\';
filename ='GeoVar_1_Dy_600_h_1000.txt';
FFref = load([savePath 'FF.mat']);

fid = fopen([savePath filename]);
Nf = FFref.FF.Nf;
Nang = FFref.FF.Nang;


E1= readmatrix([savePath filename], 'FileType', 'text', 'Range', [4, 1, Nang+3, Nf]);
E2= readmatrix([savePath filename], 'FileType', 'text', 'Range', [Nang+5, 1, Nang*2+4, Nf]);
Prad = readmatrix([savePath filename], 'FileType', 'text', 'Range', [Nang*2+6, 1, Nang*2+6, Nf]);
radEff = readmatrix([savePath filename], 'FileType', 'text', 'Range', [Nang*2+8, 1, Nang*2+8, Nf]);

 FF1 = FarField(FFref.FF.x,FFref.FF.y,E1,E2,FFref.FF.freq,Prad,radEff,...
    'coorType',FFref.FF.coorType,'polType',FFref.FF.polType,'gridType',FFref.FF.gridType,'freqUnit',FFref.FF.freqUnit,...
    'r',FFref.FF.r,'orientation',FFref.FF.orientation,'earthLocation',FFref.FF.earthLocation,'time',FFref.FF.time);


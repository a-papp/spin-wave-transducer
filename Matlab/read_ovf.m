function data = read_ovf(filename)


%% read header
file=fopen(filename);
text=textscan(file,'%s',500);
fclose(file);

idx=find(ismember(text{1},'xnodes:'));

data.I = str2double(text{1}{idx+1});
data.J = str2double(text{1}{idx+4});
data.K = str2double(text{1}{idx+7});

idx=find(ismember(text{1},'xstepsize:'));

data.dx = str2double(text{1}{idx+1});
data.dy = str2double(text{1}{idx+4});
data.dz = str2double(text{1}{idx+7});

idx=find(ismember(text{1},'time:'));
data.time = str2double(text{1}{idx+1});

%% read data
file=fopen(filename);
datacell=textscan(file,'%f %f %f',data.I*data.J*data.K,'commentstyle','#');
fclose(file);

data.X = permute(reshape(datacell{1}, [data.I data.J data.K]),[2 1 3]);
data.Y = permute(reshape(datacell{2}, [data.I data.J data.K]),[2 1 3]);
data.Z = permute(reshape(datacell{3}, [data.I data.J data.K]),[2 1 3]);


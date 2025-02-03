%% generate parameter sweep for mumax files
keywords = {'@f_exc','@file_re','@file_im'};    % keywords to be replaced

dr = 'mumax';       % folder for new files
[~,~] = mkdir(dr);  % create folder if needed
[~,~] = mkdir(dr+"/images");  % create folder if needed
f = [0:0.02:0.4 0.41:0.01:0.6 0.62:0.02:1] + 13.0;  % frequencies in GHz

save("parameters.mat","f","dr")

lines = readlines('template.mx3');  % load template file
for k = 1:length(keywords)
    line_idx(k) = find(contains(lines,keywords{k}),1); % find keyword line number
end

for i = 1:length(f)
    assignment = {['f_exc := ',num2str(f(i)),'e9'],...  % Lines to replace with
                  ['file_re := "',dr,'/CPW_DL_Real_',num2str(1),'.ohf"'],...
                  ['file_im := "',dr,'/CPW_DL_Imag_',num2str(1),'.ohf"']};  
    newfile = ['mumax_f',num2str(i),'.mx3'];    % new filename with index

    lines(line_idx) = assignment;   % replace full line with assignment
    
    writelines(lines,fullfile(dr,newfile)); % write file
end

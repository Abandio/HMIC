

function run_HMIC()
clear;
close all;
clc;

folderId = 1

% 1 x v ,n x dv ;n +x 1
dataNames = {'MSRCv1'}

%%  调试1
alpha =  [1e-7];
beta =  [1e-3];


rho = [1e-3];
ps = [0.06]

for idx = 1:length(dataNames)

    dataName = dataNames{idx};
    res = [];
    ores = [];
    sum_Zs={};
    preYs = {};
    [X,gt,kernel_type] = dataLoader(dataName);
    
    if(strcmp(kernel_type, 'None'))
        % kernel_type = 'Linear'
        kernel_type = 'Gaussian'
        % kernel_type = 'PolyPlus'
    end
    norm_type = 'SP'
    % norm_type = 'Nuclear'

    
    n = numel(gt);
    c = numel(unique(gt));
    m = c
    
    while exist(['Res/'  dataNames{idx} '-' kernel_type '-' norm_type '-' num2str(folderId)], 'dir')
        folderId = folderId + 1;
    end
    mkdir(['Res/' dataNames{idx} '-' kernel_type '-' norm_type '-' num2str(folderId)])
    
    for k = 1

    for seed = 100
        rand('seed',seed)
        for i = 1: length(alpha)
            for j = 1 :length(beta)
                for r = 1 : length(rho)
                    for p = 1:length(ps)
                    [Label,oresults,sum_Z] = HMIC(X,c,m,alpha(i),beta(j),rho(r),ps(p),seed,kernel_type,norm_type); 
                       
                    disp(['alpha=' num2str(alpha(i)) '-beta='  num2str(beta(j)) '-rho=' num2str(rho(r)) '-p=' num2str(ps(p))]);
                    results = ClusteringMeasure1(gt,Label)

                    %[ACC MIhat Purity  P R F RI];
                    res = [res;[seed m alpha(i) beta(j) rho(r) ps(p) results]];
                    ores = [ores;[seed m alpha(i) beta(j) rho(r) ps(p)  oresults']];
                    sum_Zs{end+1} = sum_Z;
                    preYs{end+1} = Label;
            end
        %% 保存
        % 循环直到找到不存在的文件夹

        save(['Res/'  dataNames{idx} '-'  kernel_type '-' norm_type '-' num2str(folderId) '/Results-' dataName '.mat'],'res')
        save(['Res/'  dataNames{idx} '-' kernel_type '-' norm_type '-' num2str(folderId) '/ObjResults-' dataName '.mat'],'ores')
        save(['Res/'  dataNames{idx} '-'  kernel_type '-' norm_type '-' num2str(folderId) '/X-' dataName '.mat'],'sum_Zs')
        save(['Res/'  dataNames{idx} '-'  kernel_type '-' norm_type '-' num2str(folderId) '/Y-' dataName '.mat'],'preYs')
                end
            end
            end
        end
    end
    end


end


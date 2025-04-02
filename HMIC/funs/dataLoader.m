function [X,gt,kernel_type] = dataLoader(dataName)
    kernerl_type='None'
    dataset = load(['T:\Files\projs\datasets\转换后的数据\' dataName '.mat']);
    if strcmp(dataName,'3sources')
        X = dataset.X';
        gt = dataset.Y;
    elseif strcmp(dataName, 'bbcsport_2view')
        X = dataset.data;
        gt = dataset.truth;
     elseif strcmp(dataName,  'CESC')
        X = dataset.fea;
        gt = dataset.gt;
    elseif strcmp(dataName, 'NGs')
        X = dataset.X;
        for cellIdx = 1:length(dataset.X)
            X{cellIdx} = X{cellIdx}';
        end
        gt = dataset.gt;
    elseif strcmp(dataName,  'prokaryotic')
        X = dataset.X;
        for cellIdx = 1:length(dataset.X)
            X{cellIdx} = X{cellIdx}';
        end
        gt = dataset.gt;
    elseif strcmp(dataName,'MSRCv1')
        X = dataset.X;
        gt = dataset.Y;
        kernel_type = 'Gaussian'
    elseif strcmp(dataName,'data_OutScene')
        X = dataset.X;
        gt = dataset.Y;
    elseif strcmp(dataName,'NH_p4660')
        X = dataset.X_train;
        for cellIdx = 1:length(X)
            X{cellIdx} = X{cellIdx}';
        end
        gt = dataset.truth;
        kernel_type = 'Gaussian'
    elseif strcmp(dataName,'Cora')
        X{1} = dataset.coracites;
        X{2} = dataset.coracontent;
        X{3} = dataset.corainbound;
        X{4} = dataset.coraoutbound';
        gt = dataset.y;
    elseif strcmp(dataName,'Citeseer')
        X = dataset.X;
        gt = dataset.y;
    elseif strcmp(dataName,'YoutubeFace_sel_fea')
        X = dataset.X';
        gt = dataset.Y;
    elseif strcmp(dataName,'AwA_fea')
        X = dataset.X';
        gt = dataset.Y;
    elseif strcmp(dataName,'Caltech101-7')
        X = dataset.X';
        gt = dataset.Y;
    elseif strcmp(dataName,'Caltech101-20')
        X = dataset.X';
        gt = dataset.Y;
    elseif strcmp(dataName,'Cifar10_test_4deep')
        X = dataset.fea;
        gt = dataset.gt;
    elseif strcmp(dataName,'Handwritten4')
        X = dataset.X;
        gt = dataset.Y;
    elseif strcmp(dataName,'BRCA')
        X = dataset.fea;
        gt = dataset.gt;
    elseif strcmp(dataName,'BRCA-2V')
        X = dataset.fea;
        gt = dataset.gt;
    elseif strcmp(dataName,'BRCA-3V')
        X = dataset.fea;
        gt = dataset.gt;
    elseif strcmp(dataName,'CESC')
        X = dataset.fea;
        gt = dataset.gt;
    elseif strcmp(dataName,'CESC-2V')
        X = dataset.fea;
        gt = dataset.gt;
    elseif strcmp(dataName,'CESC-3V')
        X = dataset.fea;
        gt = dataset.gt;
        
    elseif strcmp(dataName,'TCGA-100')
        X = dataset.fea;
        gt = dataset.gt;        
    elseif strcmp(dataName,'TCGA-200')
        X = dataset.fea;
        gt = dataset.gt;        
    elseif strcmp(dataName,'TCGA-500')
        X = dataset.fea;
        gt = dataset.gt;        
    elseif strcmp(dataName,'TCGA-1000')
        X = dataset.fea;
        gt = dataset.gt;        
    elseif strcmp(dataName,'TCGA-2000')
        X = dataset.fea;
        gt = dataset.gt;
	 elseif strcmp(dataName,'Caltech101')
        X = dataset.X;
        gt = dataset.Y';
     elseif strcmp(dataName,'SUNRGBD_fea')
        X = dataset.X';
        gt = dataset.Y;
    elseif strcmp(dataName,'STL10')
        X = dataset.X';
        gt = dataset.Y';
    elseif strcmp(dataName,'SUNRGBD_fea')
        X = dataset.X';
        gt = dataset.Y;
    elseif strcmp(dataName,'STL10')
        X = dataset.X';
        gt = dataset.Y';
    elseif strcmp(dataName,'Handwritten4')
        X = dataset.X;
        gt = dataset.Y;
    elseif strcmp(dataName,'Caltech101-all_fea')
        X = dataset.X';
        gt = dataset.Y;
    elseif strcmp(dataName,'scene15')
        X = dataset.data;
        gt = dataset.gt;
    elseif strcmp(dataName,'COIL20')
         X = dataset.data;
        for cellIdx = 1:length(dataset.data)
            X{cellIdx} = X{cellIdx}';
        end
        gt = dataset.truth;
    elseif strcmp(dataName,'ALOI-100')
        X = dataset.fea;
        gt = dataset.gt;
    elseif strcmp(dataName,'cifar100')
         X = dataset.data;
        for cellIdx = 1:length(dataset.data)
            X{cellIdx} = X{cellIdx}';
        end
        gt = dataset.truelabel{1};
    elseif strcmp(dataName,'BDGP_fea')
        X = dataset.X';
        gt = dataset.Y;

    elseif strcmp(dataName,'CCV') 
        X = dataset.X';
        gt = dataset.Y;
        kernel_type = 'Gaussian'
    elseif strcmp(dataName,'ACM') 
        X = dataset.X';
        gt = dataset.Y;
        kernel_type = 'Gaussian'

    end
    kernel_type = 'Gaussian'
end
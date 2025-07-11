% Load data
tic
outdir = '/media/proaction/data1/andre/fmroi/schaefer100';
ftpath = '/media/proaction/data1/andre/fmroi/schaefer100/zfeatmat.mat';
partinfopath = '/media/proaction/data1/andre/fmroi/participant_info.csv';

% Check if oudir exists, otherwise it creates outdir
if ~isfolder(outdir)
    mkdir(outdir)
end

auxft = load(ftpath);
varname = fieldnames(auxft);
ft_full = auxft.(varname{1});
clear auxft

partinfo = readmatrix(partinfopath,'FileType','delimitedtext');
y_full = partinfo(:,2);


% Class info
class_labels = [1, 2, 3, 4];
class_names = {'CONTROL', 'SCHZ', 'BIPOLAR', 'ADHD'};

% parameters
nreps = 100;
nperm = 1000;
pcavarexp = 50;


% Initialize results storage
pair_names = {};
cmat_all = {};
f1_all = [];
p_f1_all = [];
acc_all = [];
p_acc_all = [];
prec_all = [];
p_prec_all = [];
rec_all = [];
p_rec_all = [];
acc_real_std  = [];
prec_real_std = [];
rec_real_std  = [];
f1_real_std   = [];
acc_real_reps = [];
prec_real_reps = [];
rec_real_reps = [];
f1_real_reps = [];

fprintf('Running binary classification with LPOCV + permutation test (parallelized)...\n');

for i = 1:length(class_labels)-1
    for j = i+1:length(class_labels)
        label_i = class_labels(i);
        label_j = class_labels(j);
        name_i = class_names{i};
        name_j = class_names{j};
        pair_name = sprintf('%s_vs_%s', name_i, name_j);
        fprintf('\n=== %s ===\n', pair_name);

        % Select data for the current pair
        idx = (y_full == label_i) | (y_full == label_j);
        xpair = ft_full(idx,:);
        ypair = y_full(idx);

        % Relabel classes as 1 and 2
        ypair_bin = zeros(size(ypair));
        ypair_bin(ypair == label_i) = 1;
        ypair_bin(ypair == label_j) = 2;

        n1 = sum(ypair_bin == 1);
        n2 = sum(ypair_bin == 2);
        n = min(n1, n2);

        % Initialize metric storage
        cmat_real_all = [];
        acc_real_all = [];
        f1_real_all = [];
        prec_real_all = [];
        rec_real_all = [];
        acc_perm_all = zeros(nperm, nreps);
        f1_perm_all = zeros(nperm, nreps);
        prec_perm_all = zeros(nperm, nreps);
        rec_perm_all = zeros(nperm, nreps);

        for rep = 1:nreps
            % Downsample to the smallest class size
            idx1 = find(ypair_bin == 1);
            idx2 = find(ypair_bin == 2);
            idx1 = idx1(randperm(length(idx1), n));
            idx2 = idx2(randperm(length(idx2), n));
            idx_sampled = [idx1; idx2];

            x = xpair(idx_sampled,:);
            y = ypair_bin(idx_sampled);

            % PCA reduction
            [coeff, score, ~, ~, explained] = pca(x);
            ncomp = find(cumsum(explained) >= pcavarexp, 1);
            x_pca = score(:, 1:ncomp);
            x_pca = zscore(x_pca);

            x1 = x_pca(y == 1, :);
            x2 = x_pca(y == 2, :);
            y1 = y(y == 1);
            y2 = y(y == 2);

            % 1. Bloco real (sem permutação)
            total_pairs = n * n;
            ytrue_real = zeros(2*total_pairs,1);
            ypred_real = zeros(2*total_pairs,1);

            for a = 1:n
                for b = 1:n
                    count = (a - 1) * n + b;
                    train_idx1 = setdiff(1:n, a);
                    train_idx2 = setdiff(1:n, b);

                    xtest = [x1(a,:); x2(b,:)];
                    xtrain = [x1(train_idx1,:); x2(train_idx2,:)];

                    ytest = [y1(a,:); y2(b,:)];
                    ytrain = [y1(train_idx1,:); y2(train_idx2,:)];

                    % Real model
                    model_real = svmtrain(ytrain, xtrain, '-s 0 -t 0 -q');
                    yhat_real = svmpredict(ytest, xtest, model_real, '-q');
                    
                    ytrue_real((2*count)-1:2*count) = ytest;
                    ypred_real((2*count)-1:2*count) = yhat_real;
                end
            end

            % Compute real metrics
            acc = mean(ytrue_real == ypred_real);
            cmat = confusionmat(ytrue_real, ypred_real, 'Order', [1 2]);
            prec = diag(cmat) ./ sum(cmat, 1)';
            rec = diag(cmat) ./ sum(cmat, 2);
            f1 = 2 * (prec .* rec) ./ (prec + rec);
            prec(isnan(prec)) = 0;
            rec(isnan(rec)) = 0;
            f1(isnan(f1)) = 0;
            
            cmat_real_all = cat(3,cmat_real_all,cmat);
            acc_real_all(end+1) = acc;
            prec_real_all(end+1) = mean(prec);
            rec_real_all(end+1) = mean(rec);
            f1_real_all(end+1) = mean(f1);
            
            
            % 2. Bloco permutação (paralelo)
            ytrue_perm_all = zeros(2*total_pairs,nperm);
            ypred_perm_all = zeros(2*total_pairs,nperm);

            parfor p = 1:nperm %parfor
                ytrue_permfold = zeros(2*total_pairs, 1);
                ypred_permfold = zeros(2*total_pairs, 1);

                y_perm = y(randperm(length(y)));
                x1 = x_pca(y_perm == 1, :);
                x2 = x_pca(y_perm == 2, :);
                y1 = y_perm(y_perm == 1);
                y2 = y_perm(y_perm == 2);

                for a = 1:n
                    for b = 1:n
                        count = (a - 1) * n + b;
                        train_idx1 = setdiff(1:n, a);
                        train_idx2 = setdiff(1:n, b);

                        xtest = [x1(a,:); x2(b,:)];
                        xtrain = [x1(train_idx1,:); x2(train_idx2,:)];

                        ytest = [y1(a,:); y2(b,:)];
                        ytrain = [y1(train_idx1,:); y2(train_idx2,:)];

                        model_perm = svmtrain(ytrain, xtrain, '-s 0 -t 0 -q');
                        yhat_perm = svmpredict(ytest, xtest, model_perm, '-q');

                        ytrue_permfold((2*count)-1:2*count) = ytest;
                        ypred_permfold((2*count)-1:2*count) = yhat_perm;
                    end
                end

                ytrue_perm_all(:,p) = ytrue_permfold;
                ypred_perm_all(:,p) = ypred_permfold;
            end
            
            % Compute permuted metrics
            for p = 1:nperm
                accp = mean(ytrue_perm_all(:,p) == ypred_perm_all(:,p));
                cmatp = confusionmat(ytrue_perm_all(:,p), ypred_perm_all(:,p), 'Order', [1 2]);
                precp = diag(cmatp) ./ sum(cmatp, 1)';
                recp = diag(cmatp) ./ sum(cmatp, 2);
                f1p = 2 * (precp .* recp) ./ (precp + recp);
                precp(isnan(precp)) = 0;
                recp(isnan(recp)) = 0;
                f1p(isnan(f1p)) = 0;

                acc_perm_all(p, rep) = accp;
                prec_perm_all(p, rep) = mean(precp);
                rec_perm_all(p, rep) = mean(recp);
                f1_perm_all(p, rep) = mean(f1p);
            end
        end

        % Average real metrics
        acc_real = mean(acc_real_all);
        prec_real = mean(prec_real_all);
        rec_real = mean(rec_real_all);
        f1_real = mean(f1_real_all);

        acc_real_std(end+1)  = std(acc_real_all);
        prec_real_std(end+1) = std(prec_real_all);
        rec_real_std(end+1)  = std(rec_real_all);
        f1_real_std(end+1)   = std(f1_real_all);

        acc_real_reps(end+1,:) = acc_real_all;
        prec_real_reps(end+1,:) = prec_real_all;
        rec_real_reps(end+1,:) = rec_real_all;
        f1_real_reps(end+1,:) = f1_real_all;


        % Compute permutation distributions and p-values
        acc_perm_dist = mean(acc_perm_all, 2);
        prec_perm_dist = mean(prec_perm_all, 2);
        rec_perm_dist = mean(rec_perm_all, 2);
        f1_perm_dist = mean(f1_perm_all, 2);

        p_acc = mean(acc_perm_dist >= acc_real);
        p_prec = mean(prec_perm_dist >= prec_real);
        p_rec = mean(rec_perm_dist >= rec_real);
        p_f1 = mean(f1_perm_dist >= f1_real);

        % Store metrics
        pair_names{end+1} = pair_name;
        cmat_all{end+1} = cmat_real_all;
        acc_all(end+1) = acc_real;
        p_acc_all(end+1) = p_acc;
        prec_all(end+1) = prec_real;
        p_prec_all(end+1) = p_prec;
        rec_all(end+1) = rec_real;
        p_rec_all(end+1) = p_rec;
        f1_all(end+1) = f1_real;
        p_f1_all(end+1) = p_f1;

        % Display metrics
        fprintf('Accuracy: %.2f%% (p = %.4f)\n', 100*acc_real, p_acc);
        fprintf('Precision: %.2f%% (p = %.4f)\n', 100*prec_real, p_prec);
        fprintf('Recall: %.2f%% (p = %.4f)\n', 100*rec_real, p_rec);
        fprintf('F1-score: %.2f%% (p = %.4f)\n', 100*f1_real, p_f1);
    end
end

% Build and display result table
results_table = table(pair_names', f1_all', p_f1_all', acc_all', p_acc_all', ...
                      prec_all', p_prec_all', rec_all', p_rec_all', ...
                      'VariableNames', {'pair_name', 'f1', 'p_f1', 'acc', 'p_acc', ...
                                        'precision', 'p_precision', 'recall', 'p_recall'});

fprintf('\n=== SUMMARY TABLE ===\n');
disp(results_table);

t = toc;

[~, cpu] = system('cat /proc/cpuinfo | grep "model name" | head -1');
cpool = parcluster('local');

exp_info = struct();
exp_info.date = datetime('now');
exp_info.script_name = mfilename;
exp_info.data_path = ftpath;
exp_info.ft = ft_full;
exp_info.participant_info_path = partinfopath;
exp_info.participant_info = partinfo;
exp_info.class_labels = class_labels;
exp_info.class_names = class_names;
exp_info.nreps = nreps;
exp_info.nperm = nperm;
exp_info.ncomp = ncomp;
exp_info.varexp = pcavarexp;
exp_info.pca_criterion = sprintf('%.1f%% of explained variance', pcavarexp);
exp_info.classification_method = 'Linear SVM (LIBSVM), Leave-Pair-Out CV, downsampling to equal groups';
exp_info.matlab_version = version;
exp_info.cpu = strtrim(cpu);
exp_info.workers = cpool.NumWorkers;
exp_info.processing_duration = sprintf('Total processing time: %.2f seconds', t);
exp_info.os = ['Operating system: ', system_dependent('getos')];
exp_info.hardware = 'Machine with 2x RTX 6000 (used for parallel CPU only)';
exp_info.comments = 'Permutation test for statistical significance. All metrics are averaged across repetitions.';

% Results
exp_info.results_table = results_table;
exp_info.confmat_all = cmat_all;
exp_info.f1_all = f1_all;
exp_info.acc_allreps = f1_real_reps;
exp_info.p_f1_all = p_f1_all;
exp_info.acc_all = acc_all;
exp_info.acc_allreps = acc_real_reps;
exp_info.p_acc_all = p_acc_all;
exp_info.prec_all = prec_all;
exp_info.prec_allreps = prec_real_reps;
exp_info.p_prec_all = p_prec_all;
exp_info.rec_all = rec_all;
exp_info.rec_allreps = rec_real_reps;
exp_info.p_rec_all = p_rec_all;
exp_info.acc_real_std  = acc_real_std;
exp_info.prec_real_std = prec_real_std;
exp_info.rec_real_std  = rec_real_std;
exp_info.f1_real_std   = f1_real_std;
exp_info.pair_names = pair_names;

outfn = ['svm_results_pcavar',num2str(pcavarexp),'_nreps',num2str(nreps),'_nperm',num2str(nperm),'.mat'];
tablefn = ['svm_resultstable_pcavar',num2str(pcavarexp),'_nreps',num2str(nreps),'_nperm',num2str(nperm),'.csv'];
save(fullfile(outdir,outfn),'results_table', 'exp_info');
writetable(results_table,fullfile(outdir,tablefn));
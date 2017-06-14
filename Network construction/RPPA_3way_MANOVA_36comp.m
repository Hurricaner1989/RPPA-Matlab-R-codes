clear all; clc;

[num1, txt1, raw1] = xlsread('RPPAdata.xlsx');

LCC1 = raw1(1:46,:)';
LCC9 = raw1([1 47:91],:)';
MCF7 = raw1([1 92:136],:)';

%% 3- way mANOVA
data1 = cell2mat(LCC1(2:end,2:end));
data9 = cell2mat(LCC9(2:end,2:end));
data7 = cell2mat(MCF7(2:end,2:end));


data = [data1 data9 data7];

for j = 1:size(data,1)
    sdata = sort(data(j,:));
    for i = 1:size(data,2)
        
    end
end

C = combnk(1:size(data1,1),2);
M = [];
for i = 1:3
    for j = 1:5
        for l= 1:3
            for k = 1:3
                M = [M; [i j l]];
            end
        end
    end
end

MA3={};
for i = 1:size(C,1)
    i
    v = C(i,:);
    X = [M data(v,:)'];
    [MAO3, E, dfE] = MAOV3(X);
    MA3 = [MA3; {MAO3 E dfE}];
end

save ('RPPA_3way_MANOVA.mat')

%% Compute q-value and then proceed with the rest of the analysis
clear all; clc
load ('RPPA_3way_MANOVA.mat')

pv_MANOVA = [];
for i = 1:size(MA3,1)
    pv_MANOVA(i,:) = cell2mat(MA3{i,1}(2:end,6));
end

qv_MANOVA = [];
for i = 1:size(pv_MANOVA,2)
    qv_MANOVA(:,i) = mafdr(pv_MANOVA(:,i),'BHFDR','true');
end

% Based on the stats decide the cutoff to use to select the significant
% protein pairs.
median(qv_MANOVA,1)
mean(qv_MANOVA)
std(qv_MANOVA)
max(qv_MANOVA)
quantile(qv_MANOVA(:,1,0.6,1))

c = 3e-15 %1406 % This cut-off changes based on your data.
c1=c
c2=c;
c3=c;

in1 = find(qv_MANOVA(:,1)<c1 & qv_MANOVA(:,2)<c1 & qv_MANOVA(:,3)<c1...
    & qv_MANOVA(:,4)<c2 & qv_MANOVA(:,5)<c2 & qv_MANOVA(:,6)<c2...
    & qv_MANOVA(:,7)<c3);



%% Pairwise comparison (36 comparison)

k = numel(in1);
m = size(C,1);
cutoff = (c1*k/m)/8
sig_pair = C(in1,:);               % significant pairs from MANOVA

idx_NEG = {1:3 4:6 7:9}; %% be careful here to check when you have a new dataset..
idx_KDG = {10:12 13:15 16:18;...
    19:21 22:24 25:27;...
    28:30 31:33 34:36;...
    37:39 40:42 43:45};

T_LCC1 = {};
T_LCC9 = {};
T_MCF7 = {};

group = [ones(3,1); ones(3,1)*2];
for i = 1:size(sig_pair,1)
    pair = sig_pair(i,:);
    E = MA3{in1(i),2};
    dfE = MA3{in1(i),3};
    X1 = [M(1:45,:) data1(pair,:)'];
    X9 = [M(46:90,:) data9(pair,:)'];
    X7 = [M(91:135,:) data7(pair,:)'];
    for j = 1:4
        for k = 1:3
            % LCC1
            s1 = ['T_LCC1.G' num2str(j) '{' num2str(i) ',' num2str(k) '}'...
                '=MAOV_ttest([group X1([idx_NEG{' num2str(k)...
                '} idx_KDG{' num2str(j) ',' num2str(k) '}], [4:5])], E, dfE);'];
            eval(s1)
            % LCC9
            s9 = ['T_LCC9.G' num2str(j) '{' num2str(i) ',' num2str(k) '}'...
                '=MAOV_ttest([group X9([idx_NEG{' num2str(k)...
                '} idx_KDG{' num2str(j) ',' num2str(k) '}], [4:5])], E, dfE);'];
            eval(s9)
            % MCF7
            s7 = ['T_MCF7.G' num2str(j) '{' num2str(i) ',' num2str(k) '}'...
                '=MAOV_ttest([group X7([idx_NEG{' num2str(k)...
                '} idx_KDG{' num2str(j) ',' num2str(k) '}], [4:5])], E, dfE);'];
            eval(s7)
        end
    end
end

%% Score protein-pairs

num_sig = [];
uni_protein = unique(sig_pair(:));
for i = 1:size(data,1)
    num_sig(i,:) = [i numel(find(sig_pair(:)==i))];
end

T1  = {T_LCC1.G1; T_LCC1.G2; T_LCC1.G3; T_LCC1.G4};
T9  = {T_LCC9.G1; T_LCC9.G2; T_LCC9.G3; T_LCC9.G4};
T7  = {T_MCF7.G1; T_MCF7.G2; T_MCF7.G3; T_MCF7.G4};
pv = [];
for j = 1:4
    for i = 1:size(sig_pair,1)
        for k = 1:3
            pv(i,j,k) = cell2mat(T1{j}{i,k}(2,5));
        end
    end
end
for j = 5:8
    for i = 1:size(sig_pair,1)
        for k = 1:3
            pv(i,j,k) = cell2mat(T9{j-4}{i,k}(2,5));
        end
    end
end
for j = 9:12
    for i = 1:size(sig_pair,1)
        for k = 1:3
            pv(i,j,k) = cell2mat(T7{j-8}{i,k}(2,5));
        end
    end
end

score = [];
corr = [];
FC = [];
for j = 1:4
    for i = 1:size(sig_pair,1)
        pin = sig_pair(i,:);
        for l = 1:3
            mpv = min(pv(find(pv(:,j,l)>0),j,l));
            y1 = [(data1(pin(1),idx_NEG{l})-mean(data1(pin(1),idx_NEG{l})))'...
                (data1(pin(2),idx_NEG{l})-mean(data1(pin(2),idx_NEG{l})))'];
            y2 = [(data1(pin(1),idx_KDG{j,l})-mean(data1(pin(1),idx_KDG{j,l})))'...
                (data1(pin(2),idx_KDG{j,l})-mean(data1(pin(2),idx_KDG{j,l})))'];
            for k = 1:2
                FC(pin(k),j,l) = mean(data1(pin(k),idx_KDG{j,l}))/mean(data1(pin(k),idx_NEG{l}));
                if FC(pin(k),j,l) < 1
                    FC(pin(k),j,l) = -1/FC(pin(k),j,l);
                end
            end
            Y = [y1; y2];
            R = corrcoef(Y);
            num = num_sig(sig_pair(i,1),2)*num_sig(sig_pair(i,2),2);
            if cell2mat(T1{j}{i,l}(2,5)) == 0
                score(i,j,l) = abs(R(1,2))*log10(num/(mpv/10));
            else
                score(i,j,l) = abs(R(1,2))*log10(num/cell2mat(T1{j}{i,l}(2,5)));
            end
            corr(i,j,l) = R(1,2);
        end
    end
end
for j = 5:8
    for i = 1:size(sig_pair,1)
        pin = sig_pair(i,:);
        for l = 1:3
            mpv = min(pv(find(pv(:,j,l)>0),j,l));
            y1 = [(data9(pin(1),idx_NEG{l})-mean(data9(pin(1),idx_NEG{l})))'...
                (data9(pin(2),idx_NEG{l})-mean(data9(pin(2),idx_NEG{l})))'];
            y2 = [(data9(pin(1),idx_KDG{j-4,l})-mean(data9(pin(1),idx_KDG{j-4,l})))'...
                (data9(pin(2),idx_KDG{j-4,l})-mean(data9(pin(2),idx_KDG{j-4,l})))'];
            for k = 1:2
                FC(pin(k),j,l) = mean(data9(pin(k),idx_KDG{j-4,l}))/mean(data9(pin(k),idx_NEG{l}));
                if FC(pin(k),j,l) < 1
                    FC(pin(k),j,l) = -1/FC(pin(k),j,l);
                end
            end
            Y = [y1; y2];
            R = corrcoef(Y);
            num = num_sig(sig_pair(i,1),2)*num_sig(sig_pair(i,2),2);
            if cell2mat(T9{j-4}{i,l}(2,5)) == 0
                score(i,j,l) = abs(R(1,2))*log10(num/(mpv/10));
            else
                score(i,j,l) = abs(R(1,2))*log10(num/cell2mat(T9{j-4}{i,l}(2,5)));
            end
            corr(i,j,l) = R(1,2);
        end
    end
end
for j = 9:12
    for i = 1:size(sig_pair,1)
        pin = sig_pair(i,:);
        for l = 1:3
            mpv = min(pv(find(pv(:,j,l)>0),j,l));
            y1 = [(data7(pin(1),idx_NEG{l})-mean(data7(pin(1),idx_NEG{l})))'...
                (data7(pin(2),idx_NEG{l})-mean(data7(pin(2),idx_NEG{l})))'];
            y2 = [(data7(pin(1),idx_KDG{j-8,l})-mean(data7(pin(1),idx_KDG{j-8,l})))'...
                (data7(pin(2),idx_KDG{j-8,l})-mean(data7(pin(2),idx_KDG{j-8,l})))'];
            for k = 1:2
                FC(pin(k),j,l) = mean(data7(pin(k),idx_KDG{j-8,l}))/mean(data7(pin(k),idx_NEG{l}));
                if FC(pin(k),j,l) < 1
                    FC(pin(k),j,l) = -1/FC(pin(k),j,l);
                end
            end
            Y = [y1; y2];
            R = corrcoef(Y);
            num = num_sig(sig_pair(i,1),2)*num_sig(sig_pair(i,2),2);
            if cell2mat(T7{j-8}{i,l}(2,5)) == 0
                score(i,j,l) = abs(R(1,2))*log10(num/(mpv/10));
            else
                score(i,j,l) = abs(R(1,2))*log10(num/cell2mat(T7{j-8}{i,l}(2,5)));
            end
            corr(i,j,l) = R(1,2);
        end
    end
end
% hist(score(:,1))

F = {};
highScore_pairs = {};

% sc = 1.64; % 90% quantile cutoff for a standard normal distribution
sc = 2; % 95% quantile cutoff for a standard normal distribution
for j = 1:12
    for k = 1:3
        F{j,k} = find(score(:,j,k)>nanmean(score(:,j,k))+sc*nanstd(score(:,j,k)));
        highScore_pairs{j,k} = sig_pair(F{j,k},:);
    end
end

%%

[~,~,names] = xlsread('Protein_Symbols.xlsx'); % The protein names are used in Cystoscape.

symbol = names(2:end, [2 4 5 6]);

%% Network files needed that can be used to construct network in Cystoscape and also for topology analysis
pname = raw1(1,2:end);
filename=[{'LCC1_CYR61_48h.xlsx'} {'LCC1_CYR61_96h.xlsx'} {'LCC1_CYR61_144h.xlsx'};...
    {'LCC1_POLR2B_48h.xlsx'} {'LCC1_POLR2B_96h.xlsx'} {'LCC1_POLR2B_144h.xlsx'};...
    {'LCC1_PSMC5_48h.xlsx'} {'LCC1_PSMC5_96h.xlsx'} {'LCC1_PSMC5_144h.xlsx'};...
    {'LCC1_TOB1_48h.xlsx'} {'LCC1_TOB1_96h.xlsx'} {'LCC1_TOB1_144h.xlsx'};...
    {'LCC9_CYR61_48h.xlsx'} {'LCC9_CYR61_96h.xlsx'} {'LCC9_CYR61_144h.xlsx'};...
    {'LCC9_POLR2B_48h.xlsx'} {'LCC9_POLR2B_96h.xlsx'} {'LCC9_POLR2B_144h.xlsx'};...
    {'LCC9_PSMC5_48h.xlsx'} {'LCC9_PSMC5_96h.xlsx'} {'LCC9_PSMC5_144h.xlsx'};...
    {'LCC9_TOB1_48h.xlsx'} {'LCC9_TOB1_96h.xlsx'} {'LCC9_TOB1_144h.xlsx'}
    {'MCF7_CYR61_48h.xlsx'} {'MCF7_CYR61_96h.xlsx'} {'MCF7_CYR61_144h.xlsx'};...
    {'MCF7_POLR2B_48h.xlsx'} {'MCF7_POLR2B_96h.xlsx'} {'MCF7_POLR2B_144h.xlsx'};...
    {'MCF7_PSMC5_48h.xlsx'} {'MCF7_PSMC5_96h.xlsx'} {'MCF7_PSMC5_144h.xlsx'};...
    {'MCF7_TOB1_48h.xlsx'} {'MCF7_TOB1_96h.xlsx'} {'MCF7_TOB1_144h.xlsx'}];

for j = 4:4:12
    for k = 1:3
        M1 = [symbol(highScore_pairs{j,k}(:,1),[1 2 4]) num2cell(FC(highScore_pairs{j,k}(:,1),j,k))...
            symbol(highScore_pairs{j,k}(:,2),[1 2 4]) num2cell(FC(highScore_pairs{j,k}(:,2),j,k))...
            num2cell(corr(F{j,k},j,k))];
        xlswrite(filename{j,k},M1)
    end
end

%% Node file

% symbol(highScore_pairs1{1,k}(:,1),3)

filename2=[{'uLCC1_CYR61_48h.xlsx'} {'uLCC1_CYR61_96h.xlsx'} {'uLCC1_CYR61_144h.xlsx'};...
    {'uLCC1_POLR2B_48h.xlsx'} {'uLCC1_POLR2B_96h.xlsx'} {'uLCC1_POLR2B_144h.xlsx'};...
    {'uLCC1_PSMC5_48h.xlsx'} {'uLCC1_PSMC5_96h.xlsx'} {'uLCC1_PSMC5_144h.xlsx'};...
    {'uLCC1_TOB1_48h.xlsx'} {'uLCC1_TOB1_96h.xlsx'} {'uLCC1_TOB1_144h.xlsx'};...
    {'uLCC9_CYR61_48h.xlsx'} {'uLCC9_CYR61_96h.xlsx'} {'uLCC9_CYR61_144h.xlsx'};...
    {'uLCC9_POLR2B_48h.xlsx'} {'uLCC9_POLR2B_96h.xlsx'} {'uLCC9_POLR2B_144h.xlsx'};...
    {'uLCC9_PSMC5_48h.xlsx'} {'uLCC9_PSMC5_96h.xlsx'} {'uLCC9_PSMC5_144h.xlsx'};...
    {'uLCC9_TOB1_48h.xlsx'} {'uLCC9_TOB1_96h.xlsx'} {'uLCC9_TOB1_144h.xlsx'}
    {'uMCF7_CYR61_48h.xlsx'} {'uMCF7_CYR61_96h.xlsx'} {'uMCF7_CYR61_144h.xlsx'};...
    {'uMCF7_POLR2B_48h.xlsx'} {'uMCF7_POLR2B_96h.xlsx'} {'uMCF7_POLR2B_144h.xlsx'};...
    {'uMCF7_PSMC5_48h.xlsx'} {'uMCF7_PSMC5_96h.xlsx'} {'uMCF7_PSMC5_144h.xlsx'};...
    {'uMCF7_TOB1_48h.xlsx'} {'uMCF7_TOB1_96h.xlsx'} {'uMCF7_TOB1_144h.xlsx'}];

for j = 4:4:12
    for k = 1:3
        u = unique(highScore_pairs{j,k});
        M2 = [symbol(u,[2 3 4]) num2cell(FC(u,j,k)) num2cell(num_sig(u,2))];
        numel(u)
        xlswrite(filename2{j,k},M2)
    end
end

%% Edge file
filename3=[{'eLCC1_CYR61_48h.xlsx'} {'eLCC1_CYR61_96h.xlsx'} {'eLCC1_CYR61_144h.xlsx'};...
    {'eLCC1_POLR2B_48h.xlsx'} {'eLCC1_POLR2B_96h.xlsx'} {'eLCC1_POLR2B_144h.xlsx'};...
    {'eLCC1_PSMC5_48h.xlsx'} {'eLCC1_PSMC5_96h.xlsx'} {'eLCC1_PSMC5_144h.xlsx'};...
    {'eLCC1_TOB1_48h.xlsx'} {'eLCC1_TOB1_96h.xlsx'} {'eLCC1_TOB1_144h.xlsx'};...
    {'eLCC9_CYR61_48h.xlsx'} {'eLCC9_CYR61_96h.xlsx'} {'eLCC9_CYR61_144h.xlsx'};...
    {'eLCC9_POLR2B_48h.xlsx'} {'eLCC9_POLR2B_96h.xlsx'} {'eLCC9_POLR2B_144h.xlsx'};...
    {'eLCC9_PSMC5_48h.xlsx'} {'eLCC9_PSMC5_96h.xlsx'} {'eLCC9_PSMC5_144h.xlsx'};...
    {'eLCC9_TOB1_48h.xlsx'} {'eLCC9_TOB1_96h.xlsx'} {'eLCC9_TOB1_144h.xlsx'}
    {'eMCF7_CYR61_48h.xlsx'} {'eMCF7_CYR61_96h.xlsx'} {'eMCF7_CYR61_144h.xlsx'};...
    {'eMCF7_POLR2B_48h.xlsx'} {'eMCF7_POLR2B_96h.xlsx'} {'eMCF7_POLR2B_144h.xlsx'};...
    {'eMCF7_PSMC5_48h.xlsx'} {'eMCF7_PSMC5_96h.xlsx'} {'eMCF7_PSMC5_144h.xlsx'};...
    {'eMCF7_TOB1_48h.xlsx'} {'eMCF7_TOB1_96h.xlsx'} {'eMCF7_TOB1_144h.xlsx'}];

for j = 4:4:12
    for k = 1:3
        M3 = [];
        for l = 1:size(highScore_pairs{j,k},1)
            s3 = [char(symbol(highScore_pairs{j,k}(l,1),2)) ' (pp) ' char(symbol(highScore_pairs{j,k}(l,2),2))];
            M3 = [M3; {s3}  num2cell(corr(F{j,k}(l),j,k)) num2cell(score(F{j,k}(l),j,k))];
        end
        xlswrite(filename3{j,k},M3)
    end
end

%%
save ('RPPA_3way_MANOVA_36Comp.mat')
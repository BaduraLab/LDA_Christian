%the first 130 lines are to format the data, after that the calculations and
%preprocessing begin. If you have a question about what the inputs and
%outputs are for some functions, look at the document, 
%'function descriptions'

%% LOGS
%on 3/10/20 the calculation for the inverse is changed to
%pinv(Sw).*(MC1-MC2), this changed the eigenvalues to be more of the same
%order
%instead of eigenvector used, the argument itself is used as vectors 
%good source:
%https://cvhci.anthropomatik.kit.edu/download/publications/IEE_SIU_2LDA.pdf
%the eigenvalues are now displayed at the axis of the barplots
%variables dimensions need to be n-c=30-2=28

% 7/10/20
%inv is used because pinv is pseudo inv.
% the LDA is done on d-(n-c) variables, so bad variables are filtered
% before LDA, preferable the discrete variables 

% 13/10/20
% loading plot to see the grouping of the variables over LD1 and 2
%heatmap of mean differences std1,std2 and a score, based on that filtered
%
%%
%% loading the data
clc;
tic
% clear all;
% close all
toc
distinput=input('Do you want to see all the distributions? (yes/no) ','s');
shuffleinput=input('Do you want to wait a lot of time to shuffle the data? (yes/no) ','s');
%%
NumGenotype1 = 15; %Class size

% colors=[1 0 0 
%     1 0.5 0.5
%     1 1 1
%     0.5 1 0.5
%     0 1 0];
% extracting data of first dataset
% [numdata, textdata] = xlsread('PI3Kanalysis_transformeddata.xlsx');
% txtdata=textdata(1,2:21);
% nmdata = numdata(:, 2:end);
% genotypelabels=numdata(:,1);
%% extracting single values dataset
[dat, textdata] = xlsread('Copy of Immunology Project_ LDA table_inesbackup.xlsx','One-point variables');
txtdata=textdata(2,3:end);
% sorting the data in 2 genotype groups
[genotypelabels, sorting_index]=sort(dat(:,1),'ascend');
numdata=dat(sorting_index,:);
% getting the IDs (and order) for the single point values
IDs_single=textdata(3:end,1);
IDs_single=IDs_single(sorting_index);
%% extracting multiple values dataset
[mdat, multtextdata] = xlsread('Copy of Immunology Project_ LDA table_inesbackup.xlsx','multi-point variables (for slop');
% sorting the data in 2 genotype groups
[mgenotypelabels, sorting_index]=sort(mdat(:,1),'ascend');
multdata=mdat(sorting_index,:);
% getting the IDs (and order) for the multiple points values
IDs_multi=multtextdata(3:end,1);
IDs_multi=IDs_multi(sorting_index);

% disp(char(IDs_multi)==char(IDs_single))  %checking if the labels are the
% same for both single and multiple points

%isolating the data for the slopes
Rdata = multdata(:, 2:5);
YMaqdata = multdata(:, 6:8);
YMrev1data = multdata(: , 9:11);
YMrev2data = multdata(:,12: 14);
ELtrialdurdata=multdata(:, 15:19);
ELtimedata =multdata(:, 20:24);
ELleavingbeforedata=multdata(:, 25:29);
EL_light_air_ratiodata=multdata(:, 30:34);
ELshortstepsdata=multdata(:, 35:39);
ELjumpsdata=multdata(:, 40:44);
ELmisstepsdata=multdata(:, 45:49);

%calculating the slopes and putting them in a matrix
multfeatdata=[scatterfit(Rdata)',scatterfit(YMaqdata)',scatterfit(YMrev1data)',scatterfit(YMrev2data)',scatterfit(ELtrialdurdata)',scatterfit(ELtimedata)',scatterfit(ELleavingbeforedata)',scatterfit(EL_light_air_ratiodata)',scatterfit(ELshortstepsdata)',scatterfit(ELjumpsdata)',scatterfit(ELmisstepsdata)'];
%setting text for the variables
multfeattxt={'R','YM: aq','YM: rev1','YM: rev2','EL: trial duration','EL: time on ladder','EL: leaving before','EL: light/air ratio','EL: short steps','EL: jumps','EL: missteps'};


%% visualizing the datapoints  
nm=numdata(:, 2:end);
figure("Name","single point values")
for num=1:length(txtdata)
    subplot(4, 7,num)
    scatter(nm(1:15,num),1:15,'MarkerFaceColor',[0,0.5,0.5])
    hold on
    scatter(nm(16:end,num),1:15,'MarkerFaceColor',[1,0,0.5])
    title(char(txtdata(num)))
    hold off
end
legend('Class 1','Class 2')


figure("Name","multipoint slope values")
for num=1:length(multfeattxt)
    subplot(3, 4,num)
    scatter(multfeatdata(1:15,num),1:15,'MarkerFaceColor',[0,0.5,0.5])
    hold on
    scatter(multfeatdata(16:end,num),1:15,'MarkerFaceColor',[1,0,0.5])
    title(char(multfeattxt(num)))
    hold off
end
%% concating the single points with the slopes
txtdata=[txtdata,multfeattxt];
numdata=[numdata,multfeatdata];
genotypelabels=numdata(:, 1); %the genotype label (1 or 2)
nmdata = numdata(:, 2:end); %the actual data
%% deleting faulty variables
% delete:
%'YMrev2data'
nmdata(:,30)=[];
txtdata(:,30)=[];
% %to delete all the speed related
nmdata(:,23)=[];
nmdata(:,15)=[];
nmdata(:,7)=[];
nmdata(:,6)=[];
nmdata(:,2)=[];
txtdata(:,23)=[];
txtdata(:,15)=[];
txtdata(:,7)=[];
txtdata(:,6)=[];
txtdata(:,2)=[];
disp(txtdata)

delete=false;

if delete
    % % %SI transition test
    idelete=find(strcmp(txtdata, 'SI: total distance test'))
    % disp(txtdata(:,idelete))
    nmdata(:,idelete)=[];
    txtdata(:,idelete)=[]

    % % %{'SI: total distance BL'}
    idelete=find(strcmp(txtdata, 'YM: aq'))
    % disp(txtdata(:,idelete))
    nmdata(:,idelete)=[];
    txtdata(:,idelete)=[]

    % % {'SI: ratio transitions cups test'}
    idelete=find(strcmp(txtdata, 'SI: transitions test'))
    % disp(txtdata(:,idelete))
    nmdata(:,idelete)=[];
    txtdata(:,idelete)=[]
    
    % % {'SI: ratio transitions cups test'}
    idelete=find(strcmp(txtdata, 'EL: short steps'))
    % disp(txtdata(:,idelete))
    nmdata(:,idelete)=[];
    txtdata(:,idelete)=[]
end

% %'YMaqdata'
% disp(txtdata(:,28))
% nmdata(:,28)=[];
% txtdata(:,28)=[];

% nmdata(:,26)=[];
% txtdata(:,26)=[];
%'YM: total platform not achieved'
% disp(txtdata(:,4))
% nmdata(:,4)=[];
% txtdata(:,4)=[];

%%
nmdata_copy=nmdata; %made a copy so, I can keep the original data and make a dataset without outliers
% the 2 classes
Class_1=nmdata(1:NumGenotype1,:);
Class_2=nmdata(NumGenotype1+1:end,:);
% scatter(nmdata(1:15,4),1:15,'b')
% hold on
% scatter(nmdata(16:end,4),1:15,'r')
% hold off
% toc
[subrows, subcollumns]=subplotratio(length(txtdata)); % setting the ratio for subplots
[length_c,length_r]=size(nmdata);
%distributions before outlier filtering
yes='yes';

distbool=convertCharsToStrings(distinput);
if distbool==yes
    classdistplot(zscore_omit_nan(nmdata_copy(1:NumGenotype1,:)),txtdata,1,'',subrows, subcollumns)
    classdistplot(zscore_omit_nan(nmdata_copy(NumGenotype1+1:end,:)),txtdata,2,'',subrows, subcollumns)
    classdistplot(zscore_omit_nan(nmdata_copy),txtdata,12,'',subrows, subcollumns);
end
%% Outlier plot

figure('Name','outlier plot','Position', get(0, 'Screensize')); 

%making 2 outlier matrices (a 1 if its an outlier)

outliermatrix1=isoutlier(zscore_omit_nan(Class_1));
outliermatrix2=isoutlier(zscore_omit_nan(Class_2));

outliers=[outliermatrix1;outliermatrix2]; %overall outlier matrix
%plotting the outliers for every class
%class 1
subplot(2,1,1)
colormap([1 1 1
    0 0.5 0.5]);
imagesc(outliermatrix1);
set(gca,'XTick', 1:length(txtdata), 'xticklabel', txtdata,'YTick',[1:length(genotypelabels)]);
xtickangle(45);
ylabel ('Mice ID');
grid on;
title 'Outliers genotype 1';
subplot(2,1,2)

%class 2
colormap([1 1 1
    1 0 0.5]);
imagesc(outliermatrix2);
set(gca,'XTick', 1:length(txtdata), 'xticklabel', txtdata,'YTick',[1:length(genotypelabels)]);
xtickangle(45);
ylabel ('Mice ID');
grid on;
title 'Outliers genotype 2';
%% setting the outliers to NaN

% if the value is an outlier, put it to NaN
for i=1:size(Class_1,1)
    for j=1:length_r
        if outliermatrix1(i,j)==1
           Class_1(i,j)=NaN;
        end
    end
end

% same but for class 2
for i=1:size(Class_2,1)
    for j=1:length_r
        if outliermatrix2(i,j)==1
           Class_2(i,j)=NaN;
        end
    end
end


%% calculating z-score
% zscore and interpolation with mean, to be used for outlier visualisation
% and distribution
zscoredClass1=interpolateNaN(zscore_omit_nan(Class_1));
zscoredClass2=interpolateNaN(zscore_omit_nan(Class_2));
zscored_outliers=[zscoredClass1;zscoredClass2]; %concating them

%zscore over the whole variable
zscoreddata=zscore_omit_nan([Class_1;Class_2]);
%interpolation of the mean of the class, this is the data you will use
zscored_interp_data=[interpolateNaN(zscoreddata(1:NumGenotype1,:));interpolateNaN(zscoreddata(NumGenotype1+1:end,:))];

%zscore of the data with outliers, just in case you want to use the ouliers
for i=1:length_r
    zscore_rawdata(:,i)=zscore_omit_nan(nmdata_copy(:,i));
end

%t-score of the data, if one chooses to use it
tscored_data_filtered=(zscored_interp_data*10)+50;
tscored_data=(zscore_rawdata*10)+50;

%% normplots of the data without the outliers
%%distributions after outlier filtering
if distbool==yes
    classdistplot(zscoredClass1,txtdata,1,' without outliers',subrows, subcollumns)
    classdistplot(zscoredClass2,txtdata,2,' without outliers',subrows, subcollumns)
end
%% choosing the data to do LDA on

% this is an example dataset
% load fisheriris.mat
% NumGenotype1=50;
% Class1=meas(1:NumGenotype1,:)
% Class2=meas(NumGenotype1+1:100,:)
% txtdata=["sepal lenght","sepal width","pedal length","pedal width"];


% random data
% randomdata=randn (size(nmdata));
% % uncomment this to make a really biased random set with 2 very distinctive classes
% % randomdata=sort(randomdata,1); 
% Class1=randomdata(1:NumGenotype1,:);
% Class2=randomdata(NumGenotype1+1:end,:);

% raw data
% Class1=nmdata_copy(1:NumGenotype1,:);
% Class2=nmdata_copy(NumGenotype1+1:end,:);

%% 
% only outlier filtering
% Class1=interpolateNaN(Class_1,outliermatrix1);
% Class2=interpolateNaN(Class_2,outliermatrix2);

% only z score
% Class1=zscore_rawdata(1:NumGenotype1,:);
% Class2=zscore_rawdata(NumGenotype1+1:end,:);

% outlier filtering and z score
Class1=zscored_interp_data(1:NumGenotype1,:);
Class2=zscored_interp_data(NumGenotype1+1:end,:);
%%

[txtsorted,txtI]=sort(txtdata);
figure('Name','correlation matrix of class 1');
heatmap(txtsorted',txtsorted',corrcoef(Class1(:,txtI)));
colormap(jet);
figure('Name','correlation matrix of class 2');
heatmap(txtsorted',txtsorted',corrcoef(Class2(:,txtI)));
colormap(jet);
%%
figure('Name','correlation matrix of the whole variable');
heatmap(txtsorted',txtsorted',abs(corrcoef(zscored_interp_data(:,txtI))));
colormap(autumn);
%%
figure('Name','correlation matrix of class 2-Class1');
heatmap(txtsorted',txtsorted',corrcoef(Class2(:,txtI))-corrcoef(Class1(:,txtI)));
colormap(jet);
%%

[~, within]=scattermat([Class1;Class2],genotypelabels)
score=(abs(mean(Class1)-mean(Class2))')./(std(Class1)'+std(Class2)')
meandiff=abs(mean(Class1)-mean(Class2))'
results=[meandiff,std(Class1)',std(Class2)',score]
figure("Name","scores of the variables")
heatmap({'meandifference', 'std1','std2','score: (m1-m2)/(s1+s2)'},txtdata',results)
[scoresorted,scoreI]=sort(score,'descend')
disp('The best scores')
scorechosen=txtdata(scoreI(1:10));

%%
% amount=3;
% CorrMatrix=abs(corrcoef(zscored_interp_data));
% % disp(size(CorrMatrix,1));
% for i=1:size(CorrMatrix,1)
%     CorrMatrix(i,1:i)=0;
% end
% 
% sortedevery=sort(sort(CorrMatrix),2,'descend');
% dependend=zeros(amount,2);
% 
% for i=1:amount
%     for j=1:size(CorrMatrix,1)
%         [row,col]=find(CorrMatrix==sortedevery(end,i));
%         dependend(i,:)=[row,col];
%     end
% end
% 
% disp(sort(dependend(:,1),'descend'));
% disp(txtdata(dependend(:,2)));
% deleting=sort(dependend(:,1),'descend');
% for i=1:amount
%     disp(deleting(i))
%     Class1(:,deleting(i))=[];
%     Class2(:,deleting(i))=[];
%     txtdata(deleting(i))=[];
% end

%%
% figure('Name','zscored data')
% for num=1:length(txtdata)
%     subplot(4, 8,num)
%     scatter(Class1(:,num),1:15,'MarkerFaceColor',[0,0.5,0.5])
%     hold on
%     scatter(Class2(:,num),1:15,'MarkerFaceColor',[1,0,0.5])
%     title(char(txtdata(num)))
%     hold off
% end
% hold off

%t-score
% Class1=(Class1*10)+50;
% Class2=(Class2*10)+50;
toc

%% LDA
%LDA on the data

disp('if "badly" scaled warning: just means, that there are more variables than samples')
[LD1_proj,LD2_proj,Proj_Vector,LDs,sorted_eigenvalues,bartext_LDA,bardata_LDA,LDA_summary,indexorder]=LDA(Class1,Class2,txtdata);
%plotting the LDA
LDAplot(LD1_proj,LD2_proj,bartext_LDA,bardata_LDA,size(Class1,1),sorted_eigenvalues);
title('LDA plot on the whole data without outliers')
toc

%% choosing the best components by iteratively removing the last 1
% loopindexorder=indexorder
% for i=1:length(txtdata)-3
%     LDAorder=flip(loopindexorder); %need to flip it, because in the barplots its ascending
%     n_LDA=length(txtdata)-i; % iD LDA
%     LDA28order=LDAorder(1:n_LDA); % 10 best variables
%     C1_LD28=Class1(:,LDA28order); %new class 1 with the 10 best variables
%     C2_LD28=Class2(:,LDA28order); %new class 2 with the 10 best variables
%     LD28_text=txtdata(LDA28order);  %text data of thr 10 best variables
% 
%     [LD1_proj28,LD2_proj28,Proj_Vector28,D28,sorted_D28,bartext_LDA28,bardata_LDA28,LDA_sumary28,loopindexorder]=LDA(C1_LD28,C2_LD28,LD28_text);
%     disp(txtdata(loopindexorder(1:3)))
%     if n_LDA==10 | n_LDA==3
%         LDAplot(LD1_proj28,LD2_proj28,bartext_LDA28,bardata_LDA28,size(C1_LD28,1),sorted_D28);
%         title('LDA plot best components')
%     end
% end
% LDAorder=flip(loopindexorder); %need to flip it, because in the barplots its ascending
% n_LDA=length(txtdata)-i; % iD LDA
% LDA28order=LDAorder(1:n_LDA); % 10 best variables
% data123=[C1_LD28(:,LDA28order);C2_LD28(:,LDA28order)];
% LDvsLDplot(data123(:,1),data123(:,2),data123(:,3),NumGenotype1,txtdata(LDA28order));
% title('The 3 best components plotted vs eachother')
%% checking first 5 LD's vs the rest
% for j=1:5
%     figure('Name',sprintf("LD%d vs the rest",j),'Position', get(0, 'Screensize'))
%     for i=1:size(Class1,2)
%         subplot(subrows, subcollumns, i) ;
%         LD1new=j; % LD 1-5 to be plotted against
%         LD2new=i; % all the other LDs
%         LD1new_proj=LDA_summary(:,LD1new); % LD 1-5 to be plotted against
%         LD2new_proj=LDA_summary(:,LD2new); % all the other LDs
%         LDAsubplot(LD1new_proj,LD2new_proj,size(Class1,1),LD1new,LD2new);
%         title(sprintf('LDA plot on LD%d vs LD%d',LD1new,LD2new))
%         grid
%     end
% end

%% doing lDA only with the most weighting components
% c=2;
% n_LDA=length(txtdata)-((size(Class1,1)+size(Class2,1))-c); % 28D LDA
% LDA28order=indexorder(n_LDA+1:end); % 10 best variables
% C1_LD28=Class1(:,LDA28order); %new class 1 with the 28 best variables
% C2_LD28=Class2(:,LDA28order); %new class 2 with the 28 best variables
% txtdata=txtdata(LDA28order);  %text data of thr 28 best variables
% 
% [LD1_proj28,LD2_proj28,Proj_Vector28,D28,sorted_D28,bartext_LDA28,bardata_LDA28,LDA_sumary28,indexorder]=LDA(C1_LD28,C2_LD28,txtdata);
% LDAplot(LD1_proj28,LD2_proj28,bartext_LDA28,bardata_LDA28,size(C1_LD28,1),sorted_D28);
% title('LDA plot on 28 best components')
% toc
%%
LDAorder=flip(indexorder); %need to flip it, because in the barplots its ascending

%% LDA on first 10 features

n_LDA=10; % 10D LDA
LDA10order=LDAorder(1:n_LDA); % 10 best variables
C1_LD10=Class1(:,LDA10order); %new class 1 with the 10 best variables
C2_LD10=Class2(:,LDA10order); %new class 2 with the 10 best variables
LD10_text=txtdata(LDA10order);  %text data of thr 10 best variables

[LD1_proj10,LD2_proj10,Proj_Vector10,D10,sorted_D10,bartext_LDA10,bardata_LDA10,LDA_sumary10,indexorder10]=LDA(C1_LD10,C2_LD10,LD10_text);

LDAplot(LD1_proj10,LD2_proj10,bartext_LDA10,bardata_LDA10,size(C1_LD10,1),sorted_D10);
title('LDA plot on 10 best components')
toc


%%

LDAorder=flip(indexorder);
n_LDA=3; %3 dimensional LDA
LDA3order=LDAorder(1:n_LDA); % 3 best variables
C1_LD123=Class1(:,LDA3order); %new class 1 with the 3 best variables
C2_LD123=Class2(:,LDA3order); %new class 2 with the 3 best variables
LD123_text=txtdata(LDA3order); %text data of thr 3 best variables

[LD1_proj3,LD2_proj3,Proj_Vector3,D3,sorted_D3,bartext_LDA3,bardata_LDA3,LDA_sumary3,~]=LDA(C1_LD123,C2_LD123,LD123_text);

LDAplot(LD1_proj3,LD2_proj3,bartext_LDA3,bardata_LDA3,size(C1_LD123,1),sorted_D3);
title('LDA plot on only the 3 best components')
%% plotting the 3 most weighting components vs eachother
data123=[Class1(:,LDA3order);Class2(:,LDA3order)];
LDvsLDplot(data123(:,1),data123(:,2),data123(:,3),NumGenotype1,LD123_text);
title('The 3 best components plotted vs eachother')

%% random data with 2 dominant features
%LDA on random data with 2 predefined dominant features, to see how the
%validation responds to biased data

shufflebool=convertCharsToStrings(shuffleinput);
yes='yes';
if shufflebool==yes
    feat1=1; %index of first dominant variable
    feat2=2; %index of second dominant variable
    amount=28
    %random data with the same dimensions
    randomdata=zscore(randn(size(nmdata(:,1:amount))));
    %making it biased, by ordening the data from min to max
    randomdata(:,feat1)=sort(randomdata(:,feat1)); 
    randomdata(:,feat2)=sort(randomdata(:,feat2));

    Class1r=randomdata(1:NumGenotype1,:);
    Class2r=randomdata(NumGenotype1+1:end,:);

    disp("the chosen dominant features are:")
    disp(txtdata(feat1));
    disp(txtdata(feat2));
%     figure('Name','zscored data')
%     for num=1:length(txtdata(1:amount))
%         subplot(5, 8,num)
%         scatter(Class1r(:,num),1:15,'MarkerFaceColor',[0,0.5,0.5])
%         hold on
%         scatter(Class2r(:,num),1:15,'MarkerFaceColor',[1,0,0.5])
%         title(char(txtdata(num)))
%         hold off
%     end
    
    %doing the LDA
    [LD1_projr,LD2_projr,Proj_Vectorr,LDsr,sorted_eigenvaluesr,bartext_LDAr,bardata_LDAr,LDA_summaryr,indexorder_r]=LDA(Class1r,Class2r,txtdata(1:amount));
    LDAplot(LD1_projr,LD2_projr,bartext_LDAr,bardata_LDAr,size(Class1r,1),sorted_eigenvaluesr);
    
    title(strcat('LDA plot on random data with dominant features ',txtdata(feat1),' and  ',txtdata(feat2)))
%     %% plotting the 3 most weighting components vs eachother
    LDAorder_r=flip(indexorder_r)
    n_LDA=3; %3 dimensional LDA
    LDA3order=LDAorder_r(1:n_LDA); % 3 best variables
    LD123_text=txtdata(LDA3order); %text data of thr 3 best variables

    data123=[Class1r(:,LDA3order);Class2r(:,LDA3order)];
    LDvsLDplot(data123(:,1),data123(:,2),data123(:,3),NumGenotype1,LD123_text);
    title('The 3 best components plotted vs eachother')
    
    %% shuffling the individuals for verification
    tic
    % one can choose or to do label shuffling on the real dataset or the random
    % dominant dataset

    Shufflematrix=zscored_interp_data; % to shuffle the labels of the real data
    disp('randomizing is done on the real data')

    % Shufflematrix=randomdata;   % to validate random data with the random dominant features
    % bartext_LDA=bartext_LDAr; %need to uncomment this if random data is chosen
    % disp('randomizing is done on the dominant random data')

    counter12=0; %to count if the 2nd and 1st features of LD1 appear together while
    %shuffling
    counter123=0; %to count if the 3rd, 2nd and 1st features of LD1 appear together while
    %shuffling

    % matrix for every variable how many times it appears on every place
    count_all=zeros(length(txtdata),length(txtdata));
    shuffle_amount=10000; %amount of shuffling times
    disp('shuffling starts')
    tic
    warning('off','all')

    for i=1:shuffle_amount
        % if you want shuffle the collumns, but thats just a transformation of the feature space
    %     f_shuffle=randperm(20); 
        f_shuffle=linspace(1,length(txtdata),length(txtdata)); %length of amount of variables
        % shuflling the matrix
        shuffled=Shufflematrix(randperm(size(Shufflematrix,1)),f_shuffle);
        % new classes of shuffled individuals
        Class1_s=shuffled(1:NumGenotype1,:);
        Class2_s=shuffled(NumGenotype1+1:end,:);
        % doing the LDA
        [LD1_proj_s,LD2_proj_s,Proj_Vector_s,D_s,sorted_D_s,bartext_LDA_s,bardata_LDA_s,LDA_summary_s,~]=LDA(Class1_s,Class2_s,txtdata(f_shuffle));

    %     counting if both on 1st and 2nd
        if isequal(bartext_LDA_s(1,end),bartext_LDA(1,end)) && isequal(bartext_LDA_s(1,end-1),bartext_LDA(1,end-1))
            counter12=counter12+1;
        %     counting if on 1st 2nd and 3rd place
            if isequal(bartext_LDA_s(1,end-2),bartext_LDA(1,end-2))
                counter123=counter123+1;
            end
        end

    % counting for every possible place for every variable
    %first row is counting for every variable if its on the first place etc
    % the diagonal (from upper right to bottom left ) is counting if on their own place in respect to LD1
        for j=1:length(txtdata)
            if find(strcmp(bartext_LDA_s(1,:), txtdata(j))) ==length(txtdata)
                count_all(1,j)=count_all(1,j)+1;
            end
            for ii=1:length(txtdata)-1
                if find(strcmp(bartext_LDA_s(1,:), txtdata(j))) ==length(txtdata)-ii
                    count_all(ii+1,j)=count_all(ii+1,j)+1;
                end
            end
        end

        if mod(i,1000)==0
            disp(sprintf("shuflled %d times",i));
        end
    end
    warning('on','all')
    clc
    toc
    disp('shuffling is done')
    % %%

    % setting the indexing, so the order is the same as LD1
    for i=1:length(txtdata)
        indx(i)=find(strcmp(txtdata,bartext_LDA(1,i)));
    end
    % ordering the the numbers in the order of LD1
    ld1count=count_all(:,indx);
    flipindx=flip(indx);
    for i=1:length(txtdata)
        counter_all(i)=count_all(i,flipindx(i));
    end
    counter_all=flip(counter_all);
    disp(strcat('percentage of the 1st feature= ',num2str(ld1count(1,end)/shuffle_amount*100),'%'))
    disp(strcat('percentage of the 2nd feature= ',num2str(ld1count(2,end-1)/shuffle_amount*100),'%'))
    disp(strcat('percentage of the 3rd feature= ',num2str(ld1count(3,end-2)/shuffle_amount*100),'%'))
    disp(strcat('percentage of the 2 best features together= ',num2str(counter12/shuffle_amount*100),'%'))
    disp(strcat('percentage of the 3 best features together= ',num2str(counter123/shuffle_amount*100),'%'))

    % formatting the text for the bargraph
    counttxt=categorical(bartext_LDA(1,:));
    counttxt=reordercats(counttxt,bartext_LDA(1,:));

    % bar of amount on their own place in respect to LD1
    figure('Name', 'bar of times on their respective places after shuffling')
    b6=bar(counttxt,counter_all,'FaceColor',[0,0.5,0.5]);
    labels1 = string(b6.YData/shuffle_amount*100);
    text(b6.XData,b6.YData,strcat(labels1,' %'),'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    title(sprintf('amount of times the features of LD1 are on their own place while shuffling %d times', shuffle_amount));


    % bar of amount on 1st place after shuffling
    figure('Name', 'bar of amount on 1st place after shuffling')
    b1=bar(counttxt,ld1count(1,:),'FaceColor',[0,0.5,0.5]);
    labels1 = string(b1.YData/shuffle_amount*100);
    text(b1.XData,b1.YData,strcat(labels1,' %'),'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    title(sprintf('amount of times the features of LD1 are on the 1st place while shuffling %d times', shuffle_amount));

    % bar of amount on 1st place after shuffling of the first 3 features of LD1
    figure('Name', '3 bar of amount on 1st place after shuffling')
    count123=categorical(bartext_LDA(1,end-2:end));
    count123=reordercats(count123,bartext_LDA(1,end-2:end));
    data123=ld1count(1,end-2:end);
    b2=bar(count123,data123,'FaceColor',[0,0.5,0.5]);
    labels1 = string(b2.YData/shuffle_amount*100);
    text(b2.XData,b2.YData,strcat(labels1,' %'),'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    title(sprintf('amount of times the first 3 features of LD1 are on the 1st place while shuffling %d times', shuffle_amount));

    % bar of amount on 2nd place after shuffling
    figure('Name', 'bar of amount on 2nd place after shuffling')
    b3=bar(counttxt,ld1count(2,:),'FaceColor',[0,0.5,0.5]);
    labels1 = string(b3.YData/shuffle_amount*100);
    text(b3.XData,b3.YData,strcat(labels1,' %'),'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    title(sprintf('amount of times the features of LD1 are on the 2nd place while shuffling %d times', shuffle_amount));

    % bar of amount on 3rd place after shuffling
    figure('Name', 'bar of amount on 3rd place after shuffling')
    b4=bar(counttxt,ld1count(3,:),'FaceColor',[0,0.5,0.5]);
    labels1 = string(b4.YData/shuffle_amount*100);
    text(b4.XData,b4.YData,strcat(labels1,' %'),'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    title(sprintf('amount of times the features of LD1 are on the 3rd place while shuffling %d times', shuffle_amount));

    % bar of amount on first 3 places after shuffling
    figure('Name', 'bar of amount on first 3 places after shuffling')
    b5=bar(counttxt,ld1count','FaceColor',[0,0.5,0.5]);

    title(sprintf('amount of times the features of LD1 are on place 1, 2 and 3 while shuffling %d times', shuffle_amount));
    legend('1st place','2nd place','3rd place');
    %%

%%% histogram of the 1st place while the shuffling
figure("Name","Distribution of the first place shuffling")
h=histogram(ld1count/shuffle_amount*100);

title("Distribution of the first place shuffling")
else
    clc
    disp('wise choice')
    disp('Continuing to the heatmaps...')

end
% hold off
%% heatmaps contributions of features per LD
%deleting the conjugate equivalents
% LDs(:,10)=[];
% LDs(:,7)=[];
% LDs(:,4)=[];
Contributions=real(LDs)./sum(abs(real(LDs)))*100;
for i=1:length(txtdata)
    LDtxt{i}=strcat('LD',num2str(i));
end

%%
% heatmap of contribution percentage of every feature per LD
figure('Name', 'heatmap features per LD')
heatmap(LDtxt(1:10),txtdata(txtI),abs(Contributions(txtI,1:10)));
colormap(hot);
title('absolute contribution percentage of every feature per LD')
%%
[orderedtxt,txtorder]=sort(txtdata(LDA10order));
txtordering=LDA10order(txtorder);
figure('Name', 'heatmap features per LD')
heatmap(LDtxt(1:10),txtdata(txtordering),abs(Contributions(txtordering,1:10)));
colormap(hot);
title('absolute contribution percentage of the 10 best features per LD')

%%
%LD1 vs LD2
contr_eigevalues=sorted_eigenvalues/sum(abs(sorted_eigenvalues))*100;
figure ('Position', get(0, 'Screensize'), 'color','w');
biplot(Contributions(:,1:2), 'Color', 'r','VarLabel',txtdata);
grid on;
% xlim ([-0.45 0.45]); ylim ([-0.45 0.45]);
%if you want to see the outliers, you can use these limits
%xlim ([min(coeff(:,1)) max(coeff(:,1))]); ylim ([min(coeff(:,2)) max(coeff(:,2))]);
xlabel(sprintf('LD1 (%1.1f%%)',contr_eigevalues(1)));
ylabel(sprintf('LD2 (%1.1f%%)',contr_eigevalues(2)));
title({'LDA Loading Plot for all variables'});
return
%% LD1 vs LD2 for 10
Contributions10_1=bardata_LDA10(1,:)./sum(abs(bardata_LDA10(1,:)))*100;
Contributions10_2=bardata_LDA10(2,:)./sum(abs(bardata_LDA10(2,:)))*100;
Contributions10=[Contributions10_1;Contributions10_2]
contr_eigevalues=sorted_D10/sum(abs(sorted_D10))*100;
figure ('Position', get(0, 'Screensize'), 'color','w');
biplot(Contributions10', 'Color', 'r','VarLabel',bartext_LDA10);
grid on;
% xlim ([-0.45 0.45]); ylim ([-0.45 0.45]);
%if you want to see the outliers, you can use these limits
%xlim ([min(coeff(:,1)) max(coeff(:,1))]); ylim ([min(coeff(:,2)) max(coeff(:,2))]);
xlabel(sprintf('LD1 (%1.1f%%)',contr_eigevalues(1)));
ylabel(sprintf('LD2 (%1.1f%%)',contr_eigevalues(2)));
title({'LDA Loading Plot for the 10 best variables'});
 %%
% heatmap of sorted features by all the contribution over the range of LD_thresh
% LD_thresh=10; % the range of LD's that contribute to the order of features
% [maxord, indmax]=sort(sum(abs(Contributions(:,1:LD_thresh)),2));
% contfeatures=Contributions(indmax,1:10);
% 
% figure('Name', 'heatmap sorted for LD')
% heatmap(LDtxt(1:10),txtdata(indmax),abs(contfeatures));
% colormap(jet);
% title('sorted contribution percentage of every feature per LD')
% 
% % heatmap of sorted features by the eigenvalues of the LD's
% [maxord, indmax2]=sort(sum(sorted_eigenvalues'.*abs(Contributions)));
% contfeatures=Contributions(indmax2,:);
% figure('Name', 'heatmap weighted per LD')
% heatmap(LDtxt(1:length(txtdata)),txtdata(indmax2),contfeatures);
% colormap(jet);
% title('sorted contribution weighted by their eigenvalue of every feature per LD')
% 
% % heatmap of sorted features by the eigenvalues of the first 'LD_amount' LD's 
% LD_amount=10;
% figure('Name', 'heatmap weighted per LD for the first 10')
% heatmap(LDtxt(1:LD_amount),txtdata(indmax2),abs(contfeatures(:,1:LD_amount)));
% colormap(jet);
% title('sorted contribution weighted by their eigenvalue of every feature per LD')
% 


%% PCA
[coeff_f,score_f,latent_f] = pca(zscored_interp_data); %ALS uses least-square approx. to fill the missing data, apparently with a lot of randomness
%zscored data and rawdata
% [coeff,score,latent,tsquared,explained] = pca(zscore_rawdata,'algorithm','svd'); % ALS is good with missing data, but in the raw data there isnt any missing data
% [coeff_r,score_r,latent_r,tsquared_r,explained_r] = pca(nmdata_copy,'algorithm','als'); %ALS uses least-square approx. to fill the missing data, apparently with a lot of randomness
%%
%PCA classes
%for first 3 PCA's
C1_PC123 = score_f(1:NumGenotype1,1:3);
C2_PC123 = score_f(NumGenotype1+1:end,1:3);
txt_data123=txtdata(1:3);
% 
% %for first 10 PCA's
C1_PC10 = score_f(1:NumGenotype1,1:10);
C2_PC10 = score_f((NumGenotype1+1):end,1:10);
txt_data10=txtdata(1:10);

%PCA plot for PC1
PC1=coeff_f(:,1);
[sorted_PC1,PC1_I]=sort(abs(PC1),'ascend');
PC_txt=txtdata(PC1_I);
 figure ('Position', get(0, 'Screensize'), 'color','w');
    bar(abs(100*sorted_PC1));
    for i=1:length(sorted_PC1)
        h=bar(i,sorted_PC1(i));
        if  PC1(i)<0
            set(h,'FaceColor','k');
        else
            set(h,'FaceColor','w');
        end
        hold on
    end
    hold on
    set(gca,'xticklabel',PC_txt);
    xtickangle(20)
    xticks(1:length(PC_txt));
    grid on;
    ylabel('Contribution to PC1 (%%)');
    title('PCA contribution plot for PC1');
    hold off
%% LDA on PCA first 3 PCA
[LD1_proj_PCA,LD2_proj_PCA,Proj_VectorPCA,D,sorted_D,bartext_LDA_PCA,bardata_LDA_PCA,LDA_sumaryPCA,~]=LDA(C1_PC123,C2_PC123,txt_data123);
[bartext_LDA_PCA,bardata_LDA_PCA]=LDAPCAconv(Proj_VectorPCA,coeff_f,txtdata); %converting the weights
LDAplot(LD1_proj_PCA,LD2_proj_PCA,bartext_LDA_PCA,bardata_LDA_PCA,size(C1_PC123,1),sorted_D);
title('LDA plot on first 3 PCA')
%% LDA on first 10 PCA
[LD1_proj_PCA10,LD2_proj_PCA10,Proj_VectorPCA10,D,sorted_D,bartext_LDA_PCA10,bardata_LDA_PCA10,LDA_sumaryPCA10,~]=LDA(C1_PC10,C2_PC10,txt_data10);
[bartext_LDA_PCA10,bardata_LDA_PCA10]=LDAPCAconv(Proj_VectorPCA10,coeff_f,txtdata); %converting the weights
LDAplot(LD1_proj_PCA10,LD2_proj_PCA10,bartext_LDA_PCA10,bardata_LDA_PCA10,size(C1_PC10,1),sorted_D);
title('LDA plot on first 10 PCA')
%% PC1vsPC2vsPC3 plot
PCA_txt={'PC1','PC2','PC3'};
LDvsLDplot(score_f(:,1),score_f(:,2),score_f(:,3),NumGenotype1,PCA_txt);
title('The LDA of the 3 PCA plotted vs eachother')
toc
% %% SVM
% X = [LD1_proj;LD2_proj]';
% % X=zscoreddata;
% y = genotypelabels;
% % SVMModel= fitclinear(X,y);
% 
% figure("Name","SVM on LDA")
% g=gscatter(X(:,1),X(:,2),y);
% hold on
% 
% constrain=1/(abs((max(g(1).XData)-min(g(2).XData)))*0.05);
% SVMModel = fitcsvm(X,y,'BoxConstraint',constrain);
% classOrder = SVMModel.ClassNames;
% sv = SVMModel.SupportVectors;
% 
% plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
% legend('Class 1','Class 2','out or on the boundary')
% hold off
% %The support vectors are observations that occur on or beyond their estimated class boundaries.
% %You can adjust the boundaries (and, therefore, the number of support vectors) 
% %by setting a box constraint during training using the 'BoxConstraint' name-value pair argument.
% %% k crossvalidate
% kfold=10;
% indices = crossvalind('Kfold',genotypelabels,kfold);
% cp = classperf(genotypelabels);
% LDAspace = [LD1_proj;LD2_proj]';
% for i = 1:kfold
%     test = (indices == i);
%     train = ~test;
%     class = classify(LDAspace(test,:),LDAspace(train,:),genotypelabels(train,:));
%     classperf(cp,class,test);
% end
% fprintf('cp errorrate after k-crossvalidate on LDA= %d \n', cp.ErrorRate);
%% deleting individuals one per one

% kfolddata=zscoreddata;
% newamount=NumGenotype1;
% count_all_k=zeros(length(txtdata),length(txtdata));
% for i=1:size(zscoreddata,1)-4
%     delete_row=randperm(size(kfolddata,1),1);
%     if delete_row<=newamount
%         if newamount>1
%             newamount=newamount-1;
%         else
%             newamount=1;
%         end
%     end
%     kfolddata(delete_row,:)=[];
%     Class1_k=kfolddata(1:newamount,:);
%     Class2_k=kfolddata(newamount+1:end,:);
    
%     [LD1_proj_k,LD2_proj_k,Proj_Vector_k,D_k,sorted_k,bartext_LDA_k,bardata_LDA_k,LDA_summary_k]=LDA(Class1_k,Class2_k,txtdata);
%     
%     % counting if on their own place in respect to LD1
%     for j=1:length(txtdata)
%         if find(strcmp(bartext_LDA_k(1,:), txtdata(j))) ==length(txtdata)
%             count_all_k(1,j)=count_all_k(1,j)+1;
%         end
%         for ii=1:length(txtdata)-1
%             if find(strcmp(bartext_LDA_k(1,:), txtdata(j))) ==length(txtdata)-ii
%                 count_all_k(ii+1,j)=count_all_k(ii+1,j)+1;
%             end
%         end
%     end
% end
%% pre-processing functions
%function to interpolate all NaN with the mean of the column


%calculating the z-score on data with NaN values
function zscored_outliers=zscore_omit_nan(Class)
    %z-score in matlab gives NaN for the whole column if there is one NaN
    %present, so made loop to only zscore over the non-NaN values 

    %mean and std omitting the NaN's per column
    M1 = mean(Class, 'omitnan');
    S1 = std(Class, 'omitnan');

    %own z-score omitting the NaN's
    [length_c,length_r]=size(Class);
    zscored_outliers=zeros(length_c,length_r);


    for j=1:length_r
    n=size(Class,1)-length(find(isnan(Class(:,j))));
    %         n=size(Class,1)
        for i=1:length_c
            if isnan(Class(i,j))
                zscored_outliers(i,j)=Class(i,j); % the NaN stay in the z-score matrix, because pca can handle NaN
            else
                zscored_outliers(i,j)=sqrt(n)*(Class(i,j)-M1(j))/S1(j);  %only values get z scored
            end
        end
    end


end

function Class=interpolateNaN(Class)
    [length_c,length_r]=size(Class);
    %interpolate the NaN's
    for i=1:length_c
        for j=1:length_r
            if isnan(Class(i,j))
                % maybe interpolating with the mean of the z-score?
                Class(i,j)=mean(Class(:,j), 'omitnan');
            end
        end
    end
end

% plot function to show the distribution and the outliers
function classdistplot(Class,txtdata,n,outlier,subrows, subcollumns)
% figure('Name',strcat(sprintf('normplot for Class %d',n),outlier),'Position', get(0, 'Screensize'))
% for i=1:size(Class,2)
%     subplot(subrows, subcollumns, i) ;
%     x=Class(:,i);
%     normplot(x)
%     title(char(txtdata(i)));
%     grid
% end

figure('Name',strcat(sprintf('distribution for Class %d',n),outlier),'Position', get(0, 'Screensize'))

for i=1:size(Class,2)
    subplot(subrows, subcollumns, i) ;
    x=Class(:,i);
    [xsort,I]=sort(x);
    outliers=isoutlier(xsort);
    y=linspace(-max(x),max(x),length(x));
    cm=zeros(size(Class,1),3);
    for j=1:size(Class,1)
        cm(j,:)=[0,0.5,0.5];
        if outliers(j)
            cm(j,:)=[1,0,0.5];
        else
            cm(j,:)=[0,0.5,0.5];
        end
    end
    for j=1:size(Class,1)
        scatter(xsort(j),y(j),'filled','MarkerFaceColor',cm(j,:),'MarkerEdgeColor','k')
%         text(xsort(j)-0.5,y(j)+0.2,num2str(I(j)));
        hold on
    end
    % making the normal distribution line
    r = normrnd(0,1,[1,1000]);
    plot(sort(r),linspace(-max(x),max(x),1000),'--k')
    title(char(txtdata(i)));
    grid
    hold off

end
end

% calculating the row collumn ratio for the subplots
function [subrows, subcollumns]=subplotratio(classsize)
    mvec=zeros(classsize,1);
    for i=1:classsize
        m=mod(classsize,i);
        if m==0 
            mvec(i)=i;
        end
    end
    div=nonzeros(mvec);
    dom1=div(ceil(length(div)/2));
    dom2=classsize/dom1;
    subrows=max(dom1,dom2);
    subcollumns=min(dom1,dom2);
    if subrows/subcollumns<=2
        subrows=max(dom1,dom2);
        subcollumns=min(dom1,dom2);
    else
        [subrows, subcollumns]=subplotratio(classsize+1);
    end
end

% fit the data with a slope
function slopes=scatterfit(data)
    amount=size(data,2);
    for i=1:size(data,1)
        x=linspace(0,amount,amount);
        P=polyfit(x,data(i,:),1);
        slope = P(1);
        slopes(i)=slope;
    end
end
%% main LDA functions
% main LDA function

function [LD1_proj,LD2_proj,Proj_Vector,V_sorted,sorted_D,bartext_LDA,bardata_LDA,LDA_sumary,text_I_1]=LDA(Class1,Class2,txtdata)

    %mean of every features over all samples
    MC1=mean(Class1);
    MC2=mean(Class2);
    O_M=((size(Class1,1)*MC1)+(size(Class2,1)*MC2))/(size(Class1,1)+size(Class2,1)); %overall mean

    %between class scatter matrix
    cov1=size(Class1,1)*(MC1-O_M)'*(MC1-O_M); 
    cov2=size(Class2,1)*(MC2-O_M)'*(MC2-O_M);
    Sb=cov1+cov2;

    %within scatter matrix Class 1
    Class1_mu=Class1-repmat(MC1,size(Class1,1),1); %data-mean
    C1_mat=Class1_mu'*Class1_mu;     %Sw=(x-mean)*(x-mean)'
    
    %within scatter matrix Class 2
    Class2_mu=Class2-repmat(MC2,size(Class2,1),1); %data-mean
    C2_mat=Class2_mu'*Class2_mu;     %Sw=(x-mean)*(x-mean)'
    
    %total within scatter matrix 
    Sw=C1_mat+C2_mat;
   
    [B, W]=scattermat([Class1;Class2],[zeros(size(Class1,1),1);ones(size(Class1,1),1)]);
    %projection_vectors=eig(Sw^(-1)*Sb)
%     invSw_Sb=pinv(Sw)*Sb;
%       invSw_Sb=pinv(Sw).*(MC1-MC2);
%     invSw_Sb=Sb/Sw;
%     invSw_Sb=Sw\Sb;

    %https://nl.mathworks.com/help/matlab/ref/balance.html
    %balance before
%     argument=balance(invSw_Sb)
    LD1=pinv(W).*(MC1-MC2);
%     LD1=inv(W)*B;
    argument=LD1;
    [V,D]=eig(argument);
    %projection should be  made on eigenvector with the highest (real) eigenvalue
    %this means, the vectors that seperates the classes the best
    [sorted_D,I]=sort(abs(diag(D)),'descend');
    V_sorted=V(:,I);
    
    %sort also the arguments by their eigenvalue right?

    LD1=LD1(:,I);
    for i=1:size(D,1)
        cumulative_value(i)=sum(sorted_D(1:i))/sum(abs(sorted_D))*100; %cumulative value of the eigenvalues

    end

    % Get the 2 highest contributing vectors to plot (LD1 and LD2)
%     Proj_Vector=[LD1(:,1),LD1(:,2)]';
    Proj_Vector=[real(V_sorted(:,1)),real(V_sorted(:,2))]';

    %getting the sorting indices for the bar plots
    [LD1_sort,text_I_1]=sort(abs(Proj_Vector(1,:)),'ascend');
    [LD2_sort,text_I_2]=sort(abs(Proj_Vector(2,:)),'ascend');
    LD1_sorted=Proj_Vector(1,text_I_1);
    LD2_sorted=Proj_Vector(2,text_I_2);


    %getting the contribution percentage per feature and
    %sorting the text of the features for LD1
    contribution_value_LD1=LD1_sorted/sum(abs(LD1_sort))*100;
    sorted_features_LD1=txtdata(text_I_1);
    %getting the contribution percentage per feature and
    %sorting the text of the features for LD2
    contribution_value_LD2=LD2_sorted/sum(abs(LD2_sort))*100;
    sorted_features_LD2=txtdata(text_I_2);

    
    %formatting the text and data, to reduce the outputs
    bartext_LDA=[sorted_features_LD1;sorted_features_LD2];
    bardata_LDA=[contribution_value_LD1;contribution_value_LD2];
    

    % dot product of the projection vector and the data to map the data
    data=[Class1;Class2]; % to get the size of both classes combined
    for i=1:length(text_I_1)
        LDA_sumary(:,i)=data*LD1(:,i);        
    end
    
    % dot product of the projection vector and the data to map the data but
    % now with the sorted vectors
    LD1_proj=(data*Proj_Vector(1,:)')';
    LD2_proj=(data*Proj_Vector(2,:)')';

end

% main plot function for LDA
function LDAplot(LD1_proj,LD2_proj,bartext_LDA,bardata_LDA,Class_size,sorted_eigenvalues)
    

    eigenvalues=sorted_eigenvalues/sum(abs(sorted_eigenvalues))*100;


    %unpacking the text and data for the bar plots
    
    sorted_features_LD1=bartext_LDA(1,:);
    sorted_features_LD2=bartext_LDA(2,:);
    
    %non-absolute contributions
    magn_contr_LD1=bardata_LDA(1,:);
    magn_contr_LD2=bardata_LDA(2,:);
    
    %absolute contributions
    contribution_value_LD1=abs(magn_contr_LD1);
    contribution_value_LD2=abs(magn_contr_LD2);

    %contribution bar plot for LD2
    %setting the locus of 10 variable line
       n_10_LDA=10;
    if n_10_LDA>length(contribution_value_LD1)
    n_10_LDA=length(contribution_value_LD1);
    end
    loc_10_line=max(0,length(contribution_value_LD1)-n_10_LDA+0.5);

    %setting the locus of 3 variable line
    n_3_LDA=3;
    if n_3_LDA>length(contribution_value_LD1)
    n_3_LDA=length(contribution_value_LD1);
    end
    loc_3_line=max(0,length(contribution_value_LD1)-n_3_LDA+0.5);

    
    figure ('Name','LDA plot','Position', get(0, 'Screensize'), 'color','w')
    set(gca, 'FontName', 'Myriad Pro','Fontsize',8)
    subplot(4,4,[1,5,9])
    hold on
   
    for i=1:size(contribution_value_LD2,2)
        h=barh(i,abs(magn_contr_LD2(i)));
        if  magn_contr_LD2(i)<0
            set(h,'FaceColor','k');
        else
            set(h,'FaceColor','w');
        end
    end
    
    %plotting the lines and text
    plot([0 max(contribution_value_LD2)],[loc_10_line loc_10_line],'--k')
    text(max(contribution_value_LD2)-3,loc_10_line+0.4, sprintf('%1.1f%% cumulative contribution', sum(contribution_value_LD2(end-n_10_LDA+1:end))/sum(abs(contribution_value_LD2))*100),'Color','k');
    hold on
    plot([0 max(contribution_value_LD2)],[loc_3_line loc_3_line],'--k')
    text(max(contribution_value_LD2)-3,loc_3_line+0.4, sprintf('%1.1f%% cumulative contibution', sum(contribution_value_LD2(end-n_3_LDA+1:end))/sum(abs(contribution_value_LD2))*100),'Color','k');
  
    %contribution label of the whole bargraph
    hold on
    ytickangle(45)
    yticks(1:length(bartext_LDA(1,:)));
    set(gca,'yticklabel',sorted_features_LD2(:));
    grid on;
    xlabel(strcat(sprintf('Contribution to LD%s (%%);\n_ ', num2str(2)),sprintf(' Eigenvalue contributes %1.1f%%', eigenvalues(2))));
    hold off



    %contribution bar plot for LD1
    %setting the locus of 10 variable line
    n_10_LDA=10;
    if n_10_LDA>length(contribution_value_LD1)
    n_10_LDA=length(contribution_value_LD1);
    end
    loc_10_line=max(0,length(contribution_value_LD1)-n_10_LDA+0.5);

    %setting the locus of 3 variable line
    n_3_LDA=3;
    if n_3_LDA>length(contribution_value_LD1)
    n_3_LDA=length(contribution_value_LD1);
    end
    loc_3_line=max(0,length(contribution_value_LD1)-n_3_LDA+0.5);

    subplot(4,4,[14,16])
    hold on
    for i=1:size(contribution_value_LD1,2)
    h=bar(i,abs(magn_contr_LD1(i)));
    if  magn_contr_LD1(i)<0
        set(h,'FaceColor','k');
    else
        set(h,'FaceColor','w');
    end
    end
    hold on
    
    %plotting the lines and text
    plot([loc_10_line loc_10_line],[0 max(contribution_value_LD1)],'--k')
    text(loc_10_line-0.1, max(contribution_value_LD1)+2, sprintf('%1.1f%% cumulative contribution', sum(contribution_value_LD1(end-n_10_LDA+1:end))/sum(abs(contribution_value_LD1))*100),'Color','black');
    hold on
    plot([loc_3_line loc_3_line],[0 max(contribution_value_LD1)],'--k')
    text(loc_3_line-0.1, max(contribution_value_LD1)+2, sprintf('%1.1f%% cumulative contribution', sum(abs(contribution_value_LD1(end-n_3_LDA+1:end)))/sum(contribution_value_LD1)*100),'Color','black');
    hold on
    set(gca,'xticklabel',sorted_features_LD1);
    xtickangle(20)
    xticks(1:length(bartext_LDA(1,:)));
    grid on;
    ylabel(strcat(sprintf('Contribution to LD%s (%%);\n_ ', num2str(2)),sprintf(' Eigenvalue contributes %1.1f%%', eigenvalues(1))));

    %LDA plot 
    subplot(4,4,[2,4,10])
    for i=1:Class_size
        scatter(LD1_proj(i), LD2_proj(i),'filled','MarkerFaceColor',[0,0.5,0.5])
        ax = gca;
        %putting the origin in the middle
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        
        %these 2 only work in matlab 2020
%         xline(0,'--k');
%         yline(0,'--k');
        
        hold on
        scatter(LD1_proj(Class_size+i), LD2_proj(Class_size+i),'filled','MarkerFaceColor',[1,0,0.5])
%%uncomment for the number labels
        dx=abs(mean(LD1_proj(Class_size))/10);
        dy=abs(mean(LD2_proj(Class_size))/10);
        text(LD1_proj(i)-dx,LD2_proj(i)+dy,num2str(i));

        dx=abs(mean(LD1_proj(Class_size+i))/10);
        dy=abs(mean(LD2_proj(Class_size+i))/10);
        text(LD1_proj(Class_size+i)-dx,LD2_proj(Class_size+i)+dy,num2str(i));
        hold on
    end
    slope=eigenvalues(1)/eigenvalues(2);

    ylimits=ylim;
    line(ylimits/slope,ylimits,'LineStyle','--')
    line(-ylimits/slope,ylimits,'LineStyle','--')
    xlabel('LD1');
    ylabel('LD2');
    hold off
    
    
end

% plotting the scores of the first 3 LD's vs eachother
function LDvsLDplot(LD1,LD2,LD3,Class_size,txtdata)
    LD1_text=char(txtdata(1));
    LD2_text=char(txtdata(2));
    LD3_text=char(txtdata(3));
    figure ('Name','LDvsLD plot','Position', get(0, 'Screensize'), 'color','w')
    title('LD1 vs LD2 vs LD3')
    subplot (3,3,1);
    scatter(LD1(1:Class_size), LD2(1:Class_size),'MarkerFaceColor',[0,0.5,0.5])
    hold on
    scatter(LD1(Class_size+1:end), LD2(Class_size+1:end),'MarkerFaceColor',[1,0,0.5])
    hold on
    xlabel(LD1_text);
    ylabel(LD2_text);
    hold off
    
    subplot (3,3,2);
    scatter(LD1(1:Class_size), LD3(1:Class_size),'MarkerFaceColor',[0,0.5,0.5])
    hold on
    scatter(LD1(Class_size+1:end), LD3(Class_size+1:end),'MarkerFaceColor',[1,0,0.5])
    hold on
    xlabel(LD1_text);
    ylabel(LD3_text);
    hold off
    
    subplot (3,3,3);
    scatter(LD2(1:Class_size), LD3(1:Class_size),'MarkerFaceColor',[0,0.5,0.5])
    hold on
    scatter(LD2(Class_size+1:end), LD3(Class_size+1:end),'MarkerFaceColor',[1,0,0.5])
    hold on
    xlabel(LD2_text);
    ylabel(LD3_text);
    hold off
    
    subplot(3,3,[4,5,6,7,8,9]); 
    scatter3(LD1(1:Class_size),LD2(1:Class_size), LD3(1:Class_size),'MarkerFaceColor',[0,0.5,0.5]);
    hold on
    scatter3(LD1(Class_size+1:end),LD2(Class_size+1:end), LD3(Class_size+1:end),'MarkerFaceColor',[1,0,0.5])
    xlabel('LD1');
    ylabel('LD2');
    zlabel('LD3')
    hold off
end

%% 
% subplot function for all the different LD vs LD plots
function LDAsubplot(LD1_proj,LD2_proj,Class_size,number1,number2)
    for i=1:Class_size
        scatter(LD1_proj(i), LD2_proj(i),'MarkerFaceColor',[0,0.5,0.5])
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        dx=abs(mean(LD1_proj(Class_size))/10);
        dy=abs(mean(LD2_proj(Class_size))/10);
        text(LD1_proj(i)-dx,LD2_proj(i)+dy,num2str(i));
        hold on
        scatter(LD1_proj(Class_size+i), LD2_proj(Class_size+i),'MarkerFaceColor',[1,0,0.5])
        dx=abs(mean(LD1_proj(Class_size+i))/10);
        dy=abs(mean(LD2_proj(Class_size+i))/10);
        text(LD1_proj(Class_size+i)-dx,LD2_proj(Class_size+i)+dy,num2str(i));
        hold on
    end
    ylabel(sprintf('LD%d',number1));
    ylabel(sprintf('LD%d',number2));
    hold off
    
    
end

%% LDA on PCA function (converting)
% converting the weights for doing LDA on PCA
function [bartext_LDA,bardata_LDA]=LDAPCAconv(Proj_Vector,coef,txtdata)
    %setting how many PC's
    coeff=coef(:,1:size(Proj_Vector,2));
    %multiply in every PC dimension the weights of that PC
    for i=1:size(coeff,2)
        weigth_matrix1(:,i)=Proj_Vector(1,i)*coeff(:,i);
        weigth_matrix2(:,i)=Proj_Vector(2,i)*coeff(:,i);
    end
    %summing up all the PC weights per feature 
    weigth_vector1=sum(weigth_matrix1,2);
    weigth_vector2=sum(weigth_matrix2,2);
    %getting the sorting indices for the bar plots
    [LD1_sorted,text_I_1]=sort(abs(weigth_vector1),'ascend');
    [LD2_sorted,text_I_2]=sort(abs(weigth_vector2),'ascend');

    for i=1:size(weigth_vector1,1)
        %getting the contribution percentage per feature and
        %sorting the text of the features for LD1
        contribution_value_LD1(i)=LD1_sorted(i)/sum(LD1_sorted)*100;
        sorted_features_LD1(i)=txtdata(text_I_1(i));
        %getting the contribution percentage per feature and
        %sorting the text of the features for LD2
        contribution_value_LD2(i)=LD2_sorted(i)/sum(LD2_sorted)*100;
        sorted_features_LD2(i)=txtdata(text_I_2(i));
    end
    
    %formatting the text and data, to reduce the outputs
    bartext_LDA=[sorted_features_LD1;sorted_features_LD2];
    bardata_LDA=[contribution_value_LD1;contribution_value_LD2];
end

%% calculating the scores on different of already known LD's and the plots
% doing LDA on different LD's if the projector vectors are already known,
% this function is obsolete now and isnt used
function [LD1_proj,LD2_proj,bartext_LDA,bardata_LDA]=LDA_known(Class1,Class2,Proj_Vectors,txtdata)
    
    
    [LD1_sort,text_I_1]=sort(abs(Proj_Vectors(1,:)),'ascend');
    [LD2_sort,text_I_2]=sort(abs(Proj_Vectors(2,:)),'ascend');
    for i=1:size(Proj_Vectors,2)
        LD1_sorted(i)=Proj_Vectors(1,text_I_1(i));
        LD2_sorted(i)=Proj_Vectors(2,text_I_2(i));
    end


    for i=1:size(Proj_Vectors,2)
        %getting the contribution percentage per feature and
        %sorting the text of the features for LD1
        contribution_value_LD1(i)=LD1_sorted(i)/sum(LD1_sort)*100;
        sorted_features_LD1(i)=txtdata(text_I_1(i));
        %getting the contribution percentage per feature and
        %sorting the text of the features for LD2
        contribution_value_LD2(i)=LD2_sorted(i)/sum(LD2_sort)*100;
        sorted_features_LD2(i)=txtdata(text_I_2(i));
    end
    
    %formatting the text and data, to reduce the outputs
    bartext_LDA=[sorted_features_LD1;sorted_features_LD2];
    bardata_LDA=[contribution_value_LD1;contribution_value_LD2];
    

    % dot product of the projection vector and the data to map the data
    data=[Class1;Class2]; % to get the size of both classes combined
    for i=1:length(text_I_1)
        LDA_sumary(:,i)=data(:,text_I_1(length(text_I_1)-i+1));        
    end
    
    for i=1:size(data,1)
        LD1_proj(i)=dot(Proj_Vectors(1,:),data(i,:));
        LD2_proj(i)=dot(Proj_Vectors(2,:),data(i,:));
    end
end

function [B W]=scattermat(data,Y)
    %FUNCTION THAT CALCULATES SCATTER MATRIX:
    % B:BETWEEN CLASS SCATTER MATRIX
    % W:WITHIN CLASS SCATTER MATRIX
    %

    [~, l]=size(data); %CALCULATE SIZE OF DATA
    clases=unique(Y); %GET VECTOR OF CLASSES
    tot_clases=length(clases); %HOW MANY CLASSES
    B=zeros(l,l); %INIT B AND W
    W=zeros(l,l);
    overallmean=mean(data); %MEAN OVER ALL DATA
    for i=1:tot_clases
        clasei = find(Y==clases(i)); %GET DATA FOR EACH CLASS
        xi=data(clasei,:);
        total_length=size(Y,1);
        mci=mean(xi); %MEAN PER CLASS
        xi=xi-repmat(mci,length(clasei),1); %Xi-MeanXi
        % W=W+xi'*xi; %CALCULATE W
        W=W+((length(clasei)./total_length)*(xi'*xi));
        % B=B+length(clasei)*(mci-overallmean)'*(mci-overallmean); %CALCULATE B
        B=B+((length(clasei)./total_length)*(overallmean-mci)'*(overallmean-mci));
    end
end


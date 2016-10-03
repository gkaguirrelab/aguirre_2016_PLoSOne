% Figures and Analyses for BlindnessNOS structure paper

clear;
close all;

%  CAN SWITCH BETWEEN PROCESSING THE AUDITORY AND VISUAL DATA AT THIS STAGE

data_type = 'visual';
%data_type = 'auditory';

% Place to save figures
OutFileStem=['Figures_'];
FigIndex=1;

% Create a header to mark out the text output

TAB=char(9);  % Set up a tab to be used in text output
fprintf('\n\n\n***************************************************\n\n')

% Load a .csv data file that has all the integrated subject data and
% information. The measures in this file are the product of scripts
% produced by Rito Datta that act upon the anatomical and fMRI data.

PrimaryDataFileName=['Aguirre_2016_PLoSOne_DataMatrix.csv'];
T = readtable(PrimaryDataFileName);
clear filename;

% Extract a vector that indicates if the subject was studied under
% Protocol:
% 1 - prior to June 23, 2008, or for the "weather" or "migraine" subjects
% 2 - July 2008 or later
% 3 - MPRAGE file from Oxford

protocolNumber = T.ProtocolNumber;

% Extract a set of vectors from this table to idenitify subjects by group,
% gender, age, etc.

gender = double(strcmp(T.Gender,'Male')); % Female = 0, Male = 1
ages = T.Age;
acuity = T.Acuity;
blindnessOnset=T.VisionLossOnset;
SubjectID=T.ExternalSubjectID;
NumSubjects = length(SubjectID);

% This set of variables contain indicies that identify sub-populations of
% the total set of subjects studied

indexsight=find(strcmp(T.Group,'Sighted')==1)';
indexanophthalmic=find(strcmp(T.Group,'Anophthalmic')==1)';
indexcongenital=find(strcmp(T.Group,'Congenital')==1)';
indexpostnatal=find(strcmp(T.Group,'Postnatal')==1)';

indexCongenitalNLP=find(strcmp(T.Group,'Congenital').*strcmp(acuity,'NLP')==1)';
indexCongenitalLP=find(strcmp(T.Group,'Congenital').*strcmp(acuity,'LP')==1)';

indexrpe65=find(strcmp(T.Group,'rpe65')==1)';
indexcrb1=find(strcmp(T.Group,'crb1')==1)';
indexcep290=find(strcmp(T.Group,'CEP290')==1)';

indexlca=[indexrpe65 indexcrb1 indexcep290];
indexblind=[indexanophthalmic indexcongenital indexpostnatal indexlca];

indexpostnatal_beforepuberty = find( (strcmp(T.Group,'Postnatal').*(T.VisionLossOnset<=14)) ==1 )';
indexpostnatal_afterpuberty = find( (strcmp(T.Group,'Postnatal').*(T.VisionLossOnset>=15)) ==1 )';

% Identify the sighted subjects studied under protocol 1 and 2
indexsightProtocol1=find((strcmp(T.Group,'Sighted')==1).*(protocolNumber==1))';
indexsightProtocol2=find((strcmp(T.Group,'Sighted')==1).*(protocolNumber==2))';



% This is a matrix of the non-MPRAGE derived data

othermeasures=NaN([NumSubjects,3]);

othermeasures(:,1)=T.V1meanCBF;       % V1 CBF values already scaled by global CBF
othermeasures(:,2)=(T.lhV1SenNoiseBeta + T.rhV1SenNoiseBeta + T.lhV1RevNoiseBeta + T.rhV1RevNoiseBeta);
othermeasures(:,3)=(T.FASplenium+T.FAOpticRadiation)/2;

% Now assemble a matrix that contains the core anatomical values,
%  prior to adjustments for age, gender, overall size

if strcmp(data_type,'visual')
    raw_data=zeros(NumSubjects,8);
    measures={'lhV1thick','rhV1thick','lhV1SurfArea','rhV1SurfArea','LHpericalVol','RHpericalVol','ChiasmVol','LGNjacobian'};%,'SpleniumVol'};
    measure_type={'thickness','thickness','area','area','volume','volume','volume','jacobian'};%,'volume'};
    raw_data(:,1)= T.lhV1thick;
    raw_data(:,2)= T.rhV1thick;
    raw_data(:,3)= T.lhV1SurfArea;
    raw_data(:,4)= T.rhV1SurfArea;
    raw_data(:,5)= T.LHpericalVol;
    raw_data(:,6)= T.RHpericalVol;
    raw_data(:,7)= T.ChiasmVol;
    raw_data(:,8)= T.LGNjacobian;
    %raw_data(:,9)= T.SpleniumVol;
    
    fprintf('\n\nVISUAL DATA ANALYSIS\n\n');
    OutFileStem=[OutFileStem 'Visual/'];
end

if strcmp(data_type,'auditory')
    raw_data=zeros(NumSubjects,7);
    measures={'lhA1thick','rhA1thick','lhA1SurfArea','rhA1SurfArea','LHtranstemporalVol','RHtranstemporalVol','MGNjacobian'};
    measure_type={'thickness' 'thickness' 'area' 'area' 'volume' 'volume' 'jacobian'};
    raw_data(:,1)= T.lhA1thick;
    raw_data(:,2)= T.rhA1thick;
    raw_data(:,3)= T.lhA1SurfArea;
    raw_data(:,4)= T.rhA1SurfArea;
    raw_data(:,5)= T.LHtranstemporalVol;
    raw_data(:,6)= T.RHtranstemporalVol;
    raw_data(:,7)= T.MGNjacobian;
    
    fprintf('\n\nAUDITORY DATA ANALYSIS\n\n');
    OutFileStem=[OutFileStem 'Auditory/'];
end

sizer=size(raw_data);
NumMeasures=sizer(2);
clear sizer;


% Obtain some whole brain structural measures that will be used to adjust
% the data for individual differences in size

scalerNames={'global thickness','global area','supratentorial volume','intracranial volume'};
thicknessScaler=(T.lhCortexthick+T.rhCortexthick)/2;  % Scale by whole brain average cortical thickness
areaScaler=T.lhCortexSurfArea+T.rhCortexSurfArea; % Scale by whole brian surface area
volumeScalerA=T.SupraTenVol;
volumeScalerB=T.EstimatedTotalIntraCranialVol;

acuity=T.Acuity;

clear T;  % Done with the table


%% Create a little subject demographics table

fprintf('\n\nSubject Demographics:\n\n');
fprintf(['Sighted: ' num2str(sum(gender(indexsight))) 'm / ' num2str(length(indexsight)-sum(gender(indexsight))) 'f, age= ' num2str(mean(ages(indexsight)),'%1.0f') ' ± ' num2str(std(ages(indexsight)),'%1.0f') '\n']);
fprintf(['All blind: ' num2str(sum(gender(indexblind))) 'm / ' num2str(length(indexblind)-sum(gender(indexblind))) 'f, age= ' num2str(mean(ages(indexblind)),'%1.0f') ' ± ' num2str(std(ages(indexblind)),'%1.0f') '\n']);
fprintf(['  [Anophthalmic, congenital, LCA]: ' num2str(sum(gender([indexanophthalmic indexcongenital indexlca]))) 'm / ' num2str(length([indexanophthalmic indexcongenital indexlca])-sum(gender([indexanophthalmic indexcongenital indexlca]))) 'f, age= ' num2str(mean(ages([indexanophthalmic indexcongenital indexlca])),'%1.0f') ' ± ' num2str(std(ages([indexanophthalmic indexcongenital indexlca])),'%1.0f') '\n']);
fprintf(['  Postnatal: ' num2str(sum(gender(indexpostnatal))) 'm / ' num2str(length(indexpostnatal)-sum(gender(indexpostnatal))) 'f, age= ' num2str(mean(ages(indexpostnatal)),'%1.0f') ' ± ' num2str(std(ages(indexpostnatal)),'%1.0f') '\n']);
fprintf('\n\n');



%% Adjust the data to account for age, gender, and whole brain scaling effects
% Store and report the beta values associated with these adjustments.

% Create a vector that models the difference between protocol 1 and 2 at
% the Penn site

protocolCovariate=zeros(NumSubjects,1);
protocolCovariate(indexsightProtocol1)=1/length(indexsightProtocol1);
protocolCovariate(indexsightProtocol2)=(-1)/length(indexsightProtocol2);


adjusted_data=zeros(NumSubjects,NumMeasures);

indexAdjustGroup=[indexblind indexsight];
groupVector=zeros(NumSubjects,1);
groupVector(indexblind)=-1;
groupVector(indexsight)=1;
fprintf('_Effects of age and gender in the scaled BLIND and SIGHTED data_\n');
fprintf('Note that the effect of protcol was only modeled in the sighted control subjects\n');
fprintf('Beta (p), for each covariate [age, age^2, age^3, protocol1or2, gender, sizescalerA, sizescalerA^2, sizescalerA^3, sizerscalerB, sizescalerB^2, sizerscalerB^3]\n\n');

for i=1:NumMeasures
    if strcmp(measure_type(i),'thickness')
        [regressionMatrix,~]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup), gender(indexAdjustGroup), thicknessScaler(indexAdjustGroup));
    end
    if strcmp(measure_type(i),'area')
        [regressionMatrix,~]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup),gender(indexAdjustGroup), areaScaler(indexAdjustGroup));
    end
    if strcmp(measure_type(i),'volume')
        [regressionMatrix,~]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup),gender(indexAdjustGroup), volumeScalerA(indexAdjustGroup), volumeScalerB(indexAdjustGroup));
    end
    if strcmp(measure_type(i),'jacobian')
        [regressionMatrix,~]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup),gender(indexAdjustGroup), volumeScalerA(indexAdjustGroup), volumeScalerB(indexAdjustGroup));
    end
    regressionMatrix=[regressionMatrix;ones(1,length(indexAdjustGroup))];
    y=raw_data(indexAdjustGroup,i);
    weights=y*0+1;
    [b,~,stats] = glmfit(regressionMatrix',y,'normal','constant','off','weights',weights);
    adjusted_data((indexAdjustGroup),i)=stats.resid+b(end);
    Outline=[char(measures(i)) TAB];
    for betas=1:length(b)-1
        Outline=[Outline num2str(b(betas),'%1.2e') ' (' num2str(stats.p(betas),'%1.3f') ')' TAB];
    end
    fprintf([Outline  '\n']);
    clear b;
    clear stats;
    clear y;
    clear regressionMatrix;
end
fprintf('\n\n');


%% Test for an interaction effect of group (blind or sighted) with the other
% adjustments
fprintf('_Test for an interaction of group (blind or sighted) with the other adjustment covariates_\n');

for i=1:NumMeasures
    if strcmp(measure_type(i),'thickness')
        [~,interactionMatrix]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup), gender(indexAdjustGroup), thicknessScaler(indexAdjustGroup));
    end
    if strcmp(measure_type(i),'area')
        [~,interactionMatrix]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup),gender(indexAdjustGroup), areaScaler(indexAdjustGroup));
    end
    if strcmp(measure_type(i),'volume')
        [~,interactionMatrix]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup),gender(indexAdjustGroup), volumeScalerA(indexAdjustGroup), volumeScalerB(indexAdjustGroup));
    end
    if strcmp(measure_type(i),'jacobian')
        [~,interactionMatrix]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup),gender(indexAdjustGroup), volumeScalerA(indexAdjustGroup), volumeScalerB(indexAdjustGroup));
    end
    
    y=adjusted_data(indexAdjustGroup,i);
    interactionModel=fitglm(interactionMatrix',y,'weights',weights);
    [interaction_p,interaction_F,interaction_d] = coefTest(interactionModel);
    Outline=[char(measures(i)) TAB];
    Outline=[Outline 'F(' num2str(interaction_d) ')=' num2str(interaction_F,'%1.3f') ', p=' num2str(interaction_p,'%1.3f') ];
    fprintf([Outline  '\n']);
    
    clear b;
    clear stats;
    clear y;
    clear interactionMatrix;
end % measures
fprintf('\n\n');


clear raw_data;  % done with the raw measures
clear indexgroup;
clear groups;
clear scaled_data;


%% Report the mean and SEM of the measures in the blind and sighted groups,

fprintf('Mean ± SD for each adjusted measure for each group:\n');
fprintf('        Sighted     Blind:\n');
for i=1:NumMeasures
    Outline=['Measure ' measures{i} '  ' ];
    Outline=[Outline num2str(mean(adjusted_data(indexsight,i)),'%1.2f') '±' num2str(std(adjusted_data(indexsight,i)),'%1.2f') '  '];
    Outline=[Outline num2str(mean(adjusted_data(indexblind,i)),'%1.2f') '±' num2str(std(adjusted_data(indexblind,i)),'%1.2f') '  '];
    fprintf([Outline '\n']);
end
% now report the size scaling measures
fprintf('\n');
Outline=['Measure ' scalerNames{1} '  ' ];
Outline=[Outline num2str(mean(thicknessScaler(indexsight)),'%1.2f') '±' num2str(std(thicknessScaler(indexsight)),'%1.2f') '  '];
Outline=[Outline num2str(mean(thicknessScaler(indexblind)),'%1.2f') '±' num2str(std(thicknessScaler(indexblind)),'%1.2f') '  '];
fprintf([Outline '\n']);
Outline=['Measure ' scalerNames{2} '  ' ];
Outline=[Outline num2str(mean(areaScaler(indexsight)),'%1.2f') '±' num2str(std(areaScaler(indexsight)),'%1.2f') '  '];
Outline=[Outline num2str(mean(areaScaler(indexblind)),'%1.2f') '±' num2str(std(areaScaler(indexblind)),'%1.2f') '  '];
fprintf([Outline '\n']);
Outline=['Measure ' scalerNames{3} '  ' ];
Outline=[Outline num2str(mean(volumeScalerA(indexsight)),'%1.2f') '±' num2str(std(volumeScalerA(indexsight)),'%1.2f') '  '];
Outline=[Outline num2str(mean(volumeScalerA(indexblind)),'%1.2f') '±' num2str(std(volumeScalerA(indexblind)),'%1.2f') '  '];
fprintf([Outline '\n']);
Outline=['Measure ' scalerNames{4} '  ' ];
Outline=[Outline num2str(mean(volumeScalerB(indexsight)),'%1.2f') '±' num2str(std(volumeScalerB(indexsight)),'%1.2f') '  '];
Outline=[Outline num2str(mean(volumeScalerB(indexblind)),'%1.2f') '±' num2str(std(volumeScalerB(indexblind)),'%1.2f') '  '];
fprintf([Outline '\n']);
fprintf('\n\n');


fprintf('Cohen-s d effect sizes for blind vs. sighted for each adjusted measure:\n');
for i=1:NumMeasures
    d=(mean(adjusted_data(indexsight,i))-mean(adjusted_data(indexblind,i)))/ ( ( std(adjusted_data(indexsight,i)) + std(adjusted_data(indexblind,i)) ) /2);
    fprintf(['measure ' measures{i} ': ' num2str(d,'%1.1f') '\n']);
    clear d;
end
% now report the size scaling measures
fprintf('\n');
Outline=['Measure ' scalerNames{1} '  ' ];
Outline=[Outline num2str( ( mean(thicknessScaler(indexsight)) - mean(thicknessScaler(indexblind)) ) / ( std(thicknessScaler(indexsight)) + std(thicknessScaler(indexblind)) ) /2,'%1.2f')];
fprintf([Outline '\n']);
Outline=['Measure ' scalerNames{2} '  ' ];
Outline=[Outline num2str( ( mean(areaScaler(indexsight)) - mean(areaScaler(indexblind)) ) / ( std(areaScaler(indexsight)) + std(areaScaler(indexblind)) ) /2,'%1.2f')];
fprintf([Outline '\n']);
Outline=['Measure ' scalerNames{3} '  ' ];
Outline=[Outline num2str( ( mean(volumeScalerA(indexsight)) - mean(volumeScalerA(indexblind)) ) / ( std(volumeScalerA(indexsight)) + std(volumeScalerA(indexblind)) ) /2,'%1.2f')];
fprintf([Outline '\n']);
Outline=['Measure ' scalerNames{4} '  ' ];
Outline=[Outline num2str( ( mean(volumeScalerB(indexsight)) - mean(volumeScalerB(indexblind)) ) / ( std(volumeScalerB(indexsight)) + std(volumeScalerB(indexblind)) ) /2,'%1.2f')];
fprintf([Outline '\n']);

fprintf('\n\n');



fprintf('p value of a t-test between the blind and sighted for each adjusted measure:\n');
for i=1:NumMeasures
    [~,pVal]=ttest2(adjusted_data(indexsight,i),adjusted_data(indexblind,i));
    fprintf(['measure ' measures{i} ': ' num2str(pVal,'%1.5f') '\n']);
    clear d;
end
% now report the size scaling measures
fprintf('\n');
Outline=['Measure ' scalerNames{1} '  ' ];
Outline=[Outline num2str( ( mean(thicknessScaler(indexsight)) - mean(thicknessScaler(indexblind)) ) / ( std(thicknessScaler(indexsight)) + std(thicknessScaler(indexblind)) ) /2,'%1.2f')];
fprintf([Outline '\n']);
Outline=['Measure ' scalerNames{2} '  ' ];
Outline=[Outline num2str( ( mean(areaScaler(indexsight)) - mean(areaScaler(indexblind)) ) / ( std(areaScaler(indexsight)) + std(areaScaler(indexblind)) ) /2,'%1.2f')];
fprintf([Outline '\n']);
Outline=['Measure ' scalerNames{3} '  ' ];
Outline=[Outline num2str( ( mean(volumeScalerA(indexsight)) - mean(volumeScalerA(indexblind)) ) / ( std(volumeScalerA(indexsight)) + std(volumeScalerA(indexblind)) ) /2,'%1.2f')];
fprintf([Outline '\n']);
Outline=['Measure ' scalerNames{4} '  ' ];
Outline=[Outline num2str( ( mean(volumeScalerB(indexsight)) - mean(volumeScalerB(indexblind)) ) / ( std(volumeScalerB(indexsight)) + std(volumeScalerB(indexblind)) ) /2,'%1.2f')];
fprintf([Outline '\n']);
fprintf('\n\n');


fprintf('p value of a t-test between the congenitally and postnatally blind for each adjusted measure:\n');
for i=1:NumMeasures
    [~,pVal]=ttest2(adjusted_data(indexpostnatal,i),adjusted_data([indexanophthalmic indexcongenital indexlca],i));
    fprintf(['measure ' measures{i} ': ' num2str(pVal,'%1.5f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('p value of a t-test between the sighted and postnatally blind for each adjusted measure:\n');
for i=1:NumMeasures
    [~,pVal]=ttest2(adjusted_data(indexpostnatal,i),adjusted_data([indexsight],i));
    fprintf(['measure ' measures{i} ': ' num2str(pVal,'%1.5f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('p value of a t-test between congenitally blind with light perception and without for each adjusted measure:\n');
for i=1:NumMeasures
    [~,pVal]=ttest2(adjusted_data([indexCongenitalLP],i),adjusted_data([indexCongenitalNLP],i));
    fprintf(['measure ' measures{i} ': ' num2str(pVal,'%1.5f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('Spearmans rho correlation between age at blindness onset and measure in the postnatally blind:\n');
fprintf('Num postnatal subjects with a defined date of vision loss = ');
fprintf(num2str(length(indexpostnatal)));
fprintf('\n');

for i=1:NumMeasures
    [rho,pVal]=corr(adjusted_data(indexpostnatal,i),blindnessOnset(indexpostnatal),'type','Spearman');
    fprintf(['measure ' measures{i} '- rho = ' num2str(rho,'%1.5f') ', p = ' num2str(pVal,'%1.5f') '\n']);
end
fprintf('\n\n');




%% Create bar plots of the group means and SEM for measure

figtmp=figure('name','Group Means ± SEM each measure');

groupLabels={'[anop + congenital + lca' 'postnatal' 'sighted'};

for i=1:NumMeasures
    subplot(3,3,i);
    tmpGroupDimScoreMeans=[mean(adjusted_data([indexanophthalmic indexcongenital indexlca],i)) mean(adjusted_data(indexpostnatal,i)) mean(adjusted_data(indexsight,i))];
    tmpGroupDimScoreSEM=[std(adjusted_data([indexanophthalmic indexcongenital indexlca],i))/sqrt(length([indexanophthalmic indexcongenital indexlca])) std(adjusted_data(indexpostnatal,i))/sqrt(length(indexpostnatal)) std(adjusted_data(indexsight,i))/sqrt(length(indexsight))];
    bar([1 2 3],tmpGroupDimScoreMeans);
    hold on
    errorbar([1 2 3],tmpGroupDimScoreMeans,tmpGroupDimScoreSEM,'o');
    hold on
    xlim([0 4]);
    set(gca,'Xtick',1:3,'XTickLabel',groupLabels);
    pbaspect([2 1 1])
    box off;
    title([measures{i}]);
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;






%% Prepare for the cluster analysis

% Z-transform the data within each measure
%  The adjustment is performed on all subjects (including LCA1), but the
%  parameters of the adjustment are calculated excluding the LCA1 subjects.

z_data=NaN([NumSubjects NumMeasures]);

for i=1:NumMeasures
    z_data(:,i)=(adjusted_data(:,i) - mean(adjusted_data([indexblind indexsight],i))) / std(adjusted_data([indexblind indexsight],i));
end

% run through the measures and multiply any "thickness" value by -1
% to render it a measure of "thinness", which is more easily interpreted as
% altered in blindness

for i=1:NumMeasures
    if strcmp(measure_type(i),'thickness')
        z_data(:,i)=z_data(:,i)*(-1);
    end
end

data=z_data;
clear z_data;


% Conduct a cluster analysis to demonstrate that the grouping of the
%   measures makes sense

clusterDistance='euclidean';
clusterLinkage='ward';
clusterMax=3;


Y_blind = pdist(data(indexblind,:)',clusterDistance);
Z_blind = linkage(Y_blind,clusterLinkage);
T_blind = cluster(Z_blind,'maxclust',clusterMax);
leafOrder_blind = optimalleaforder(Z_blind,Y_blind);

Y_sight = pdist(data(indexsight,:)',clusterDistance);
Z_sight = linkage(Y_sight,clusterLinkage);
T_sight = cluster(Z_sight,'maxclust',clusterMax);
leafOrder_sight = optimalleaforder(Z_sight,Y_sight);

% Bootstrap analysis to determine how often this precise tree results from
% the data
numBoots=100;

clusterReplicationBlind=[0,0,0];
for boots=1:numBoots
    indexResample=randsample(indexblind,length(indexblind),true);
    Y_boot = pdist(data(indexResample,:)',clusterDistance);
    Z_boot = linkage(Y_boot,clusterLinkage);
    T_boot = cluster(Z_boot,'maxclust',clusterMax);
    if T_boot(1)==T_boot(2)
        clusterReplicationBlind(1)=clusterReplicationBlind(1)+1;
    end
    if (T_boot(3)==T_boot(4) || T_boot(4)==T_boot(5) || T_boot(5)==T_boot(6))
        clusterReplicationBlind(2)=clusterReplicationBlind(2)+1;
    end
    if T_boot(7)==T_boot(8)
        clusterReplicationBlind(3)=clusterReplicationBlind(3)+1;
    end
end

fprintf('Replication of cluster assignment over 1000 bootstraps,\n');
fprintf('Blind subject data: \n');
for i=1:3
    fprintf(['Cluster ' num2str(i) ':  ' num2str(clusterReplicationBlind(i)/numBoots,'%1.3f')  '\n']);
end
fprintf('\n\n');



clusterReplicationSight=[0,0,0];
for boots=1:numBoots
    indexResample=randsample(indexsight,length(indexsight),true);
    Y_boot = pdist(data(indexResample,:)',clusterDistance);
    Z_boot = linkage(Y_boot,clusterLinkage);
    T_boot = cluster(Z_boot,'maxclust',clusterMax);
    if T_boot(1)==T_boot(2)
        clusterReplicationSight(1)=clusterReplicationSight(1)+1;
    end
    if (T_boot(3)==T_boot(4) || T_boot(4)==T_boot(5) || T_boot(5)==T_boot(6))
        clusterReplicationSight(2)=clusterReplicationSight(2)+1;
    end
    if T_boot(7)==T_boot(8)
        clusterReplicationSight(3)=clusterReplicationSight(3)+1;
    end
end

fprintf('Replication of cluster assignment over 1000 bootstraps,\n');
fprintf('Sighted subject data: \n');
for i=1:3
    fprintf(['Cluster ' num2str(i) ':  ' num2str(clusterReplicationSight(i)/numBoots,'%1.3f')  '\n']);
end
fprintf('\n\n');




% Plot the dendrograms

leafOrder=[1 2 3 5 4 6 7 8];

figtmp=figure('name','Dendrogram Blind');
H = dendrogram(Z_blind,'Orientation','left','ColorThreshold','default','Reorder',leafOrder);
set(H,'LineWidth',2)
title('Dendrogram Blind');
xlim([4 15]);
pbaspect([2 2 1])
box off;
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;

leafOrder=[1 2 3 4 5 6 7 8];

figtmp=figure('name','Dendrogram Sight');
H = dendrogram(Z_sight,'Orientation','left','ColorThreshold','default','Reorder',leafOrder);
set(H,'LineWidth',2)
title('Dendrogram Sight');
xlim([4 15]);
pbaspect([2 2 1])
box off;
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


fprintf('\n\nThe cophenetic correlation coefficient for the hierarchical cluster trees:\n\n');
fprintf(['Sight: ' num2str(cophenet(Z_sight,Y_sight),'%1.2f') '\n']);
fprintf(['Blind: ' num2str(cophenet(Z_blind,Y_blind),'%1.2f') '\n']);
fprintf('\n\n');




% Create the distances matrix implied by the dendrograms, first for the
% blind

% Create the distance matrices
[~,dendrogramDistances]=cophenet(Z_blind,Y_blind);
dendrogramDistancesMatrix=squareform(dendrogramDistances);

% Plot the observed distance matrix
figtmp=figure('name','Blind subjects measure distances');
m=squareform(Y_blind);
for i=1:NumMeasures
    m(i,i)=NaN;
end
h = imagesc(m,[4 12]);
axis square;
cmap='gray';
colormap(cmap);
colorbar;
set(gca, 'visible', 'off') ;
title('Blind subjects measure distances');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;

% Plot the modeled distance matrix
figtmp=figure('name','Distances matrix modeled by dendrogram for blind');
m=dendrogramDistancesMatrix;
for i=1:NumMeasures
    m(i,i)=NaN;
end
h = imagesc(m,[4 12]);
axis square;
cmap='gray';
colormap(cmap);
colorbar;
set(gca, 'visible', 'off') ;
title('Distances matrix modeled by dendrogram for blind');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;

% Now create the distances matrix implied by the dendrograms for the
% sighted

% Create the distance matrices
[~,dendrogramDistances]=cophenet(Z_sight,Y_sight);
dendrogramDistancesMatrix=squareform(dendrogramDistances);

% Plot the observed distance matrix
figtmp=figure('name','Sighted subjects measure distances');
m=squareform(Y_sight);
for i=1:NumMeasures
    m(i,i)=NaN;
end
h = imagesc(m,[4 12]);
axis square;
cmap='gray';
colormap(cmap);
colorbar;
set(gca, 'visible', 'off') ;
title('Sighted subjects measure distances');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;

% Plot the modeled distance matrix
figtmp=figure('name','Distances matrix modeled by dendrogram for sighted');
m=dendrogramDistancesMatrix;
for i=1:NumMeasures
    m(i,i)=NaN;
end
h = imagesc(m,[4 12]);
axis square;
cmap='gray';
colormap(cmap);
colorbar;
set(gca, 'visible', 'off') ;
title('Distances matrix modeled by dendrogram for sighted');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;


%% Average the data within clusters

all_score=NaN([NumSubjects NumMeasures]);
all_score(:,1)=mean(data(:,[1 2]),2);
all_score(:,2)=mean(data(:,[3 4 5 6]),2);
all_score(:,3)=mean(data(:,[7 8]),2);


%% Report means and t-tests of the groups in clusters

fprintf('p value of a t-test between the blind and sighted for the clusters:\n');
for i=1:3
    [~,pVal]=ttest2(all_score(indexsight,i),all_score([indexblind],i));
    fprintf(['cluster ' num2str(i) ': ' num2str(pVal,'%1.5f') '\n']);
end
fprintf('\n\n');

fprintf('p value of a t-test between the congenitally and postnatally blind for the clusters:\n');
for i=1:3
    [~,pVal]=ttest2(all_score(indexpostnatal,i),all_score([indexanophthalmic indexcongenital indexlca],i));
    fprintf(['cluster ' num2str(i) ': ' num2str(pVal,'%1.5f') '\n']);
end
fprintf('\n\n');

fprintf('p value of a t-test between the sighted and postnatally blind for the clusters:\n');
for i=1:3
    [~,pVal]=ttest2(all_score(indexpostnatal,i),all_score([indexsight],i));
    fprintf(['cluster ' num2str(i) ': ' num2str(pVal,'%1.5f') '\n']);
end
fprintf('\n\n');


% Create bar plots of the group means and SEM for each dimension

figtmp=figure('name','Group Means ± SEM each dimension');
titles={'clust1' 'clust2' 'clust3'};

PCLabels={'sighted' 'postnatal' 'congen'};

for i=1:3
    subplot(3,1,i);
    tmpGroupDimScoreMeans=[mean(all_score(indexsight,i)) mean(all_score(indexpostnatal,i)) mean(all_score([indexanophthalmic indexcongenital indexlca],i))];
    tmpGroupDimScoreSEM=[std(all_score(indexsight,i))/sqrt(length(indexsight)) std(all_score(indexpostnatal,i))/sqrt(length(indexpostnatal)) std(all_score([indexanophthalmic indexcongenital indexlca],i))/sqrt(length([indexanophthalmic indexcongenital indexlca])) ];
    bar([1 2 3],tmpGroupDimScoreMeans);
    hold on
    errorbar([1 2 3],tmpGroupDimScoreMeans,tmpGroupDimScoreSEM,'o');
    hold on
    xlim([0 4]);
    ylim([-1 1]);
    set(gca,'Xtick',1:3,'XTickLabel',PCLabels);
    pbaspect([2 1 1])
    box off;
    title([titles{i}]);
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;






%%%%%%%%%%%%%%%%%%%%%%%
%% Time to work on the "other measures"
%%%%%%%%%%%%%%%%%%%%%%%


%% Report the number of subjects who have each kind of "other" measure

fprintf('\n\nNumber of subjects with additional brain measures:\n\n');
fprintf( '         CBF  Sentences FA\n');
fprintf(['Sighted: ' num2str(length(indexsight)-sum(isnan(othermeasures(indexsight,:)))) '\n']);
fprintf(['Blind:   ' num2str(length(indexblind)-sum(isnan(othermeasures(indexblind,:)))) '\n']);
fprintf('\n\n');


LabelsPC={'Clust1' 'Clust2' 'Clust3'};
LabelsMeasures={'CBF' 'Sentences' 'FA'};


%% Test if there are differences between the control groups in protocol 1
% vs. protocol 2 measures


[~,p_vals,~,stats] = ttest2(othermeasures(indexsightProtocol1,:),othermeasures(indexsightProtocol2,:));
fprintf('Unpaired t-test comparing the other measure in the control group between protocol 1 and 2.\n');
fprintf('Measure, t-test.\n');
for i=1:3
    Outline=[LabelsMeasures{i}];
    Outline=[Outline ', t(' num2str(stats.df(i)) ') = ' num2str(stats.tstat(i),'%1.2f') ', p = ' num2str(p_vals(i),'%1.2g') '\n'];
    fprintf(Outline);
end
fprintf('\n\n');

%% Report the means, and t-test of the sighted vs. the blind

% Obtain the mean and SEM of the non-MPRAGE measures, broken down
%  by sub-group. Create a plot of this.

LabelsMeasures={'CBF' 'Sentences' 'FA'};
LabelsGroups={'Sight' 'Postnatal' 'Congen'};

OtherMeasureMean=zeros(3,3);
OtherMeasureSEM=zeros(3,3);

% Use a common ranges to plot these results
 
YRanges=[0 3; -5 10; 0.3 0.6];
XRanges=YRanges;

for i=1:3 % Loop across measures
    OtherMeasureMean(i,1)=nanmean(othermeasures(indexsight,i));
    OtherMeasureSEM(i,1)=nanstd(othermeasures(indexsight,i)) / sqrt(length(indexsight));
    
    OtherMeasureMean(i,2)=nanmean(othermeasures(indexpostnatal,i));
    OtherMeasureSEM(i,2)=nanstd(othermeasures(indexpostnatal,i)) / sqrt(length(indexpostnatal));
    
    OtherMeasureMean(i,3)=nanmean(othermeasures(indexcongenital,i));
    OtherMeasureSEM(i,3)=nanstd(othermeasures(indexcongenital,i)) / sqrt(length(indexcongenital));
end

figtmp=figure('name','Mean of other measures by group ±SEM');
titles={'CBF' 'Sentences' 'FA'};

for i=1:3
    subplot(3,1,i);
    errorbar([1 2 3],OtherMeasureMean(i,:),OtherMeasureSEM(i,:),'o');
    hold on
    ylim(YRanges(i,:));
    xlim([1 3]);
    set(gca,'Xtick',1:3,'XTickLabel',LabelsGroups);
    pbaspect([2 1 1])
    box off;
    title([titles{i}]);
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


[~,p_vals,~,stats] = ttest2(othermeasures(indexsight,:),othermeasures(indexblind,:));
fprintf('Compare blind and all sighted on the CBF, Sentences, and FA measures.\n');
fprintf('Measure: [mean sighted , mean blind], t-test.\n');
for i=1:3
    Outline=[LabelsMeasures{i} ': [' num2str(nanmean(othermeasures(indexsight,i)),'%1.2f') ', ' num2str(nanmean(othermeasures(indexblind,i)),'%1.2f') ']'];
    Outline=[Outline ', t(' num2str(stats.df(i)) ') = ' num2str(stats.tstat(i),'%1.2f') ', p = ' num2str(p_vals(i),'%1.2g') '\n'];
    fprintf(Outline);
end
fprintf('\n\n');

[~,p_vals,~,stats] = ttest2(othermeasures(indexpostnatal,:),othermeasures(indexcongenital,:));
fprintf('Compare early blind and late blindon the CBF, Sentences, and FA measures.\n');
fprintf('Measure: [mean early , mean late blind], t-test.\n');
for i=1:3
    Outline=[LabelsMeasures{i} ': [' num2str(nanmean(othermeasures(indexpostnatal,i)),'%1.2f') ', ' num2str(nanmean(othermeasures(indexcongenital,i)),'%1.2f') ']'];
    Outline=[Outline ', t(' num2str(stats.df(i)) ') = ' num2str(stats.tstat(i),'%1.2f') ', p = ' num2str(p_vals(i),'%1.2g') '\n'];
    fprintf(Outline);
end
fprintf('\n\n');







% Calculate F tests associated with the three anatomical clusters modeling the other
% measures

d=3; % three data variables to model (CSF, sentences, FA)
p=3; % the three anatomical variation dimensions
n=NumSubjects; % NumSubjects rows, although most are NaNs

% Get the other measures, with the mean across all subjects removed
Y=[othermeasures([indexsight indexblind],1) othermeasures([indexsight indexblind],2) othermeasures([indexsight indexblind],3)];
for i=1:3
    Y(:,i)=Y(:,i)-nanmean(Y(:,i));
end

% Create a matrix of regressors which has just each anatomical dimension.
X=[all_score([indexsight indexblind],1) all_score([indexsight indexblind],2) all_score([indexsight indexblind],3)];
for i=1:p
    X(:,i)=X(:,i)-nanmean(X(:,i));
end

clear B; clear Bint; clear modelSEM;
fprintf('F test of the 3 PC model of the measures for all subjects:\n');
for i=1:3
    [B(i,:),Bint,~,~,tmp_stats]=regress(Y(:,i),X);
    modelSEM(i,:)=(Bint(:,2)-Bint(:,1))/3.92; % convert 95% CI to SEM
    denom_f_df=length([indexsight indexblind])-sum(isnan(Y(:,i)))-3;
    Outline=[titles{i} '- F (3,' num2str(denom_f_df) ') = ' num2str(tmp_stats(2),'%1.2f') ', p = ' num2str(tmp_stats(3),'%1.2g') '\n'];
    fprintf(Outline);
end
fprintf('\n\n');


%% Report the beta weights on the cluster components and SEM of these

% Create a plot of this

figtmp=figure('name','Weights ±SEM on the components');
titles={'CBF' 'Sentences' 'FA'};
PCLabels={'Clust1' 'Clust2' 'Clust3'};
PCPlotYRages=[-.25 .25; -1.5 1.5; -.025 .025];


for i=1:3
    subplot(3,1,i);
    bar([1 2 3],B(i,1:3));
    hold on
    errorbar([1 2 3],B(i,1:3),modelSEM(i,1:3),'o');
    hold on
    xlim([0 5]);
    ylim(PCPlotYRages(i,:));
    set(gca,'Xtick',1:3,'XTickLabel',PCLabels);
    pbaspect([2 1 1])
    box off;
    title([titles{i}]);
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;

%%Calculate what the p-value is for each anatomical predictor once
% the group effect is removed.

% Get a matrix of the other measures with no mean centering
Y=[othermeasures(:,1) othermeasures(:,2) othermeasures(:,3)];

% Create a matrix of regressors which has the intercept, group effect, and
% then the mean-centered effect of each anatomical dimension.
groupRegressor=zeros(1,NumSubjects);
groupRegressor(indexblind)=-1;
groupRegressor(indexsight)=1;
X=[all_score(:,1) all_score(:,2) all_score(:,3) groupRegressor'];
for i=1:p
    X(indexsight,i)=X(indexsight,i)-nanmean(X(indexsight,i));
    X(indexblind,i)=X(indexblind,i)-nanmean(X(indexblind,i));
    indexnan=find(isnan(othermeasures(:,i)));
    X(indexnan,i)=nan;
    indexpos=find(X(:,i)==1);
    X(indexpos,i)=1/indexpos;
    indexneg=find(X(:,i)==(-1));
    X(indexneg,i)=(-1)/indexneg;
end

fprintf('Here is the p-value for each anatomical dimension predictor, after removing group effect:\n');
clear B;
for i=1:3
    [B(i,:),~,tmp_stats] = glmfit([X ones(1,NumSubjects)'],Y(:,i),'normal','constant','off');
    tmp_stats.p(1:3)
end
fprintf('\n\n');


%% Create three scatter plots of model fit vs. measure

figtmp=figure('name','Correlation with other measures');

titles={'Fit vs CBF All' 'Fit vs Sentences All' 'Fit vs FA All'};
 
YRanges=[0 3; -5 10; 0.3 0.6];

for i=1:3
    subplot(2,2,i);
    hold on
    plot((B(i,1:5)*[X ones(1,NumSubjects)']'),Y(:,i),'.b');
    r=nancorr((B(i,1:5)*[X ones(1,NumSubjects)']'),Y(:,i));
    lsline
    ylim(YRanges(i,:));
    xlim(YRanges(i,:));
    plot((B(i,1:5)*[X(indexsight,:) ones(1,length(indexsight))']'),Y(indexsight,i),'.k');
    plot((B(i,1:5)*[X([indexpostnatal indexcongenital],:) ones(1,length([indexpostnatal indexcongenital]))']'),Y([indexpostnatal indexcongenital],i),'.r');
    pbaspect([1 1 1])
    axis equal;
    box off;
    title([titles{i} '  r=' num2str(r)]);
hold off
end
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;














% Post-hoc: Check if there is a difference in correlation of the early
% blind compared to the late blind / sighted of anatomical dimension 3 and
% the optic radiation / splenium FA

fprintf(['Corr of blind FA with anatomy: ' num2str(nancorr((B(3,1:3)*X([indexcongenital indexpostnatal],1:3)'),othermeasures([indexcongenital indexpostnatal],3)),'%1.2f') '\n' ]);
fprintf(['Corr of sighted FA with anatomy: ' num2str(nancorr((B(3,1:3)*X([ indexsight],1:3)'),othermeasures([ indexsight],3)),'%1.2f') '\n' ]);
fprintf('\n\n');

fprintf(['Corr of blind BOLD fMRI with anatomy: ' num2str(nancorr((B(2,1:3)*X([indexcongenital indexpostnatal],1:3)'),othermeasures([indexcongenital indexpostnatal],2)),'%1.2f') '\n' ]);
fprintf(['Corr of sighted BOLD fMRI with anatomy: ' num2str(nancorr((B(2,1:3)*X([indexsight],1:3)'),othermeasures([indexsight],2)),'%1.2f') '\n' ]);
fprintf('\n\n');

% Report here the correlation within group between V1 THICKNESS and BOLD
% fMRI cross modal response. To do so, we take the "all_score" covariate,
% which has the data after being grouped into a cluster. The first cluster
% had earlier been multiplied by (-1) to render it V1 THINNESS. To make the
% reported value easy to interpret, we multiply again by -1 here to return
% the value to cortical thickness.
fprintf(['Corr of blind BOLD fMRI with V1 thickness: ' num2str((-1)*nancorr(all_score([indexcongenital indexpostnatal],1),othermeasures([indexcongenital indexpostnatal],2)),'%1.2f') '\n' ]);
fprintf(['Corr of sighted BOLD fMRI with V1 thickness: ' num2str((-1)*nancorr(all_score([indexsight],1),othermeasures([indexsight],2)),'%1.2f') '\n' ]);
fprintf('\n\n');

% Drop a footer to mark out the end of BlindnessNOS text output

fprintf('\n\n\n***************************************************\n\n')



close all;

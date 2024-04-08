%% Jackknife

%% Get Data
clc
clear

totalsubj=13;
grouplevel_res.grouplevel_alldata=[];

% Combine data from all participants
for subj=1:totalsubj

    % Load
    filename=sprintf('ResultsSubjectOfficial%i',subj);
    load(filename,'alldataclean')

    % Combine
    grouplevel_res.grouplevel_alldata=[grouplevel_res.grouplevel_alldata; alldataclean(:,1:16)];
end

save('grouplevel_res','grouplevel_res')

%% Calculate Normal Mean at Each Contrast
% Calculate Objective Performance for each Contrast and Condition
fullresultsobj=[];

for idxcond=1:3 % for each condition
    conditiondata=grouplevel_res.grouplevel_alldata(find(grouplevel_res.grouplevel_alldata{:,"Condition"}==idxcond),:); % get trials
    results=[];
    for idxcont=1:10 % find each contrast
        contrastindex=find(conditiondata{:,"Contrast Level"}==idxcont);
        results(:,idxcont)=conditiondata{contrastindex,"Correct/Incorrect"};
    end
    fullresultsobj{idxcond}=results; % calculate mean across trials for each contrast
end

% Prepare Data for Psignifit (levels x correct x total)
% Prepare data structs
gl_dataforpsignifit.irrobj=(1:10)'; % Stimulus Levels
gl_dataforpsignifit.intobj=gl_dataforpsignifit.irrobj; % copy for other conditions
gl_dataforpsignifit.rhyobj=gl_dataforpsignifit.irrobj;

gl_dataforpsignifit.irrsubj=(1:10)'; % Stimulus Levels
gl_dataforpsignifit.intsubj=gl_dataforpsignifit.irrsubj; % copy for other conditions
gl_dataforpsignifit.rhysubj=gl_dataforpsignifit.irrsubj;

% Get Correct/incorrect and binary subj response for each condition
indxirr=find(grouplevel_res.grouplevel_alldata{:,"Condition"} == 3);
indxint=find(grouplevel_res.grouplevel_alldata{:,"Condition"} == 2);
indxrhy=find(grouplevel_res.grouplevel_alldata{:,"Condition"} == 1);

irreval=grouplevel_res.grouplevel_alldata{indxirr,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
inteval=grouplevel_res.grouplevel_alldata{indxint,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
rhyeval=grouplevel_res.grouplevel_alldata{indxrhy,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};

% Calculate nCorrect per condition and Stimulus Level

% Get Correct/incorrect and binary subj response for each stimulus Level
for stimlev=1:10
    indxtempirr=find(irreval(:,1) == stimlev);
    indxtempint=find(inteval(:,1) == stimlev);
    indxtemprhy=find(rhyeval(:,1) == stimlev);

    irreval_currlevobj=irreval(indxtempirr,2);
    inteval_currlevobj=inteval(indxtempint,2);
    rhyeval_currlevobj=rhyeval(indxtemprhy,2);

    irreval_currlevsubj=irreval(indxtempirr,3);
    inteval_currlevsubj=inteval(indxtempint,3);
    rhyeval_currlevsubj=rhyeval(indxtemprhy,3);

    gl_dataforpsignifit.irrobj(stimlev,2)=sum(irreval_currlevobj);
    gl_dataforpsignifit.intobj(stimlev,2)=sum(inteval_currlevobj);
    gl_dataforpsignifit.rhyobj(stimlev,2)=sum(rhyeval_currlevobj);

    gl_dataforpsignifit.irrsubj(stimlev,2)=sum(irreval_currlevsubj);
    gl_dataforpsignifit.intsubj(stimlev,2)=sum(inteval_currlevsubj);
    gl_dataforpsignifit.rhysubj(stimlev,2)=sum(rhyeval_currlevsubj);

    % Get the total n of trials at this level for each condition
    gl_dataforpsignifit.irrobj(stimlev,3)=length(find(grouplevel_res.grouplevel_alldata{indxirr,"Contrast Level"} == stimlev));
    gl_dataforpsignifit.intobj(stimlev,3)=length(find(grouplevel_res.grouplevel_alldata{indxint,"Contrast Level"} == stimlev));
    gl_dataforpsignifit.rhyobj(stimlev,3)=length(find(grouplevel_res.grouplevel_alldata{indxrhy,"Contrast Level"} == stimlev));

    gl_dataforpsignifit.irrsubj(stimlev,3)=gl_dataforpsignifit.irrobj(stimlev,3);
    gl_dataforpsignifit.intsubj(stimlev,3)=gl_dataforpsignifit.intobj(stimlev,3);
    gl_dataforpsignifit.rhysubj(stimlev,3)=gl_dataforpsignifit.rhyobj(stimlev,3);
end

% Mean Objective

for idxcond=1:3
    means.objmeancondcontrast(idxcond,:)=mean(fullresultsobj{idxcond},1,'omitnan');
end

%%%% Calculate Subjective Performance for each Contrast and Condition

fullresultssubj=[];

% Make table samples x contrasts (subjective)
for idxcond=1:3 % for each condition
    conditiondata=grouplevel_res.grouplevel_alldata(find(grouplevel_res.grouplevel_alldata{:,"Condition"}==idxcond),:); % get trials
    results=[];
    for idxcont=1:10 % find each contrast
        contrastindex=find(conditiondata{:,"Contrast Level"}==idxcont);
        results(:,idxcont)=conditiondata{contrastindex,"Binary Visibility"};
    end
    fullresultssubj{idxcond}=results; % calculate mean across trials for each contrast
end

% Mean Subjective
means.subjmeancondcontrast=[];
for idxcond=1:3
    means.subjmeancondcontrast(idxcond,:)=mean(fullresultssubj{idxcond},1,'omitnan');
end

% Save

grouplevel_res.means.obj=means.objmeancondcontrast;
grouplevel_res.means.subj=means.subjmeancondcontrast;
grouplevel_res.dataforpsignifit=gl_dataforpsignifit;

save('grouplevel_res','grouplevel_res')

%% Fit and Plot Curve to normal average
% Prepare struct with fitting options
fitting_options_obj=struct;
fitting_options_subj=struct;

fitting_options_obj.expType = 'nAFC';
fitting_options_obj.expN = 2;
fitting_options_subj.expType = 'YesNo';

fitting_options_obj.sigmoidName = 'logistic';
fitting_options_subj.sigmoidName = 'logistic';

% Fit Objective Curves

grouplevel_res.avgcurves.irrobj = psignifit(grouplevel_res.dataforpsignifit.irrobj,fitting_options_obj);
grouplevel_res.avgcurves.intobj = psignifit(grouplevel_res.dataforpsignifit.intobj,fitting_options_obj);
grouplevel_res.avgcurves.rhyobj = psignifit(grouplevel_res.dataforpsignifit.rhyobj,fitting_options_obj);

% Fit Subjective Curves
grouplevel_res.avgcurves.irrsubj = psignifit(grouplevel_res.dataforpsignifit.irrsubj,fitting_options_subj);
grouplevel_res.avgcurves.intsubj = psignifit(grouplevel_res.dataforpsignifit.intsubj,fitting_options_subj);
grouplevel_res.avgcurves.rhysubj = psignifit(grouplevel_res.dataforpsignifit.rhysubj,fitting_options_subj);


% Plot
plotOptions1.lineColor = [0,0,0];
plotOptions1.dataColor = [0,0,0];
plotOptions1.CIthresh = true;
plotOptions1.dataSize=2;
plotOptions1.lineWidth = 1.5;
plotOptions2.lineColor = [1,0,0];
plotOptions2.dataColor = [1,0,0];
plotOptions2.lineWidth = 1.5;
plotOptions2.dataSize=2;
plotOptions2.CIthresh = true;
plotOptions3.lineColor = [0,0.7,0.5];
plotOptions3.dataColor = [0,0.7,0.5];
plotOptions3.dataSize=2;
plotOptions3.lineWidth = 1.5;
plotOptions3.CIthresh = true;

figure(1);
[hline]=plotPsych(grouplevel_res.avgcurves.irrobj,plotOptions1);
hold on
[hline2]=plotPsych(grouplevel_res.avgcurves.intobj,plotOptions2);
[hline3]=plotPsych(grouplevel_res.avgcurves.rhyobj,plotOptions3);
title('Objective Average Curves')
legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')
hold off

figure(2);
[hline]=plotPsych(grouplevel_res.avgcurves.irrsubj,plotOptions1);
hold on
[hline2]=plotPsych(grouplevel_res.avgcurves.intsubj,plotOptions2);
[hline3]=plotPsych(grouplevel_res.avgcurves.rhysubj,plotOptions3);
title('Subjective Average Curves')
legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')
hold off

save('grouplevel_res','grouplevel_res')
%% Calculate Average Leaving one out
for partcount=1:totalsubj % thirteen times calculate the means

    reduced_data=[];

    % Leave one participant out
    partsample=1:totalsubj;
    partsample(partcount)=[];
    for partidx=partsample
        idx=find(grouplevel_res.grouplevel_alldata{:,"Subject"}==partidx);
        partdata=grouplevel_res.grouplevel_alldata(idx,:);
        reduced_data=[reduced_data; partdata];
    end
    
    % Calculate Objective Performance for each Contrast and Condition
    fullresultsobj=[];

    for idxcond=1:3 % for each condition
        conditiondata=reduced_data(find(reduced_data{:,"Condition"}==idxcond),:); % get trials
        results=[];
        for idxcont=1:10 % find each contrast
            contrastindex=find(conditiondata{:,"Contrast Level"}==idxcont);
            results(:,idxcont)=conditiondata{contrastindex,"Correct/Incorrect"};
        end
        fullresultsobj{idxcond}=results; % calculate mean across trials for each contrast
    end

    % Prepare Data for Psignifit (levels x correct x total)
    % Prepare data structs
    gl_dataforpsignifit_jk{partcount}.irrobj=(1:10)'; % Stimulus Levels
    gl_dataforpsignifit_jk{partcount}.intobj=gl_dataforpsignifit_jk{partcount}.irrobj; % copy for other conditions
    gl_dataforpsignifit_jk{partcount}.rhyobj=gl_dataforpsignifit_jk{partcount}.irrobj;

    gl_dataforpsignifit_jk{partcount}.irrsubj=(1:10)'; % Stimulus Levels
    gl_dataforpsignifit_jk{partcount}.intsubj=gl_dataforpsignifit_jk{partcount}.irrsubj; % copy for other conditions
    gl_dataforpsignifit_jk{partcount}.rhysubj=gl_dataforpsignifit_jk{partcount}.irrsubj;

    % Get Correct/incorrect and binary subj response for each condition
    indxirr=find(reduced_data{:,"Condition"} == 3);
    indxint=find(reduced_data{:,"Condition"} == 2);
    indxrhy=find(reduced_data{:,"Condition"} == 1);

    irreval=reduced_data{indxirr,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
    inteval=reduced_data{indxint,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
    rhyeval=reduced_data{indxrhy,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};

    % Calculate nCorrect per condition and Stimulus Level

    % Get Correct/incorrect and binary subj response for each stimulus Level
    for stimlev=1:10
        indxtempirr=find(irreval(:,1) == stimlev);
        indxtempint=find(inteval(:,1) == stimlev);
        indxtemprhy=find(rhyeval(:,1) == stimlev);

        irreval_currlevobj=irreval(indxtempirr,2);
        inteval_currlevobj=inteval(indxtempint,2);
        rhyeval_currlevobj=rhyeval(indxtemprhy,2);

        irreval_currlevsubj=irreval(indxtempirr,3);
        inteval_currlevsubj=inteval(indxtempint,3);
        rhyeval_currlevsubj=rhyeval(indxtemprhy,3);

        gl_dataforpsignifit_jk{partcount}.irrobj(stimlev,2)=sum(irreval_currlevobj);
        gl_dataforpsignifit_jk{partcount}.intobj(stimlev,2)=sum(inteval_currlevobj);
        gl_dataforpsignifit_jk{partcount}.rhyobj(stimlev,2)=sum(rhyeval_currlevobj);

        gl_dataforpsignifit_jk{partcount}.irrsubj(stimlev,2)=sum(irreval_currlevsubj);
        gl_dataforpsignifit_jk{partcount}.intsubj(stimlev,2)=sum(inteval_currlevsubj);
        gl_dataforpsignifit_jk{partcount}.rhysubj(stimlev,2)=sum(rhyeval_currlevsubj);

        % Get the total n of trials at this level for each condition
        gl_dataforpsignifit_jk{partcount}.irrobj(stimlev,3)=length(find(reduced_data{indxirr,"Contrast Level"} == stimlev));
        gl_dataforpsignifit_jk{partcount}.intobj(stimlev,3)=length(find(reduced_data{indxint,"Contrast Level"} == stimlev));
        gl_dataforpsignifit_jk{partcount}.rhyobj(stimlev,3)=length(find(reduced_data{indxrhy,"Contrast Level"} == stimlev));

        gl_dataforpsignifit_jk{partcount}.irrsubj(stimlev,3)=gl_dataforpsignifit_jk{partcount}.irrobj(stimlev,3);
        gl_dataforpsignifit_jk{partcount}.intsubj(stimlev,3)=gl_dataforpsignifit_jk{partcount}.intobj(stimlev,3);
        gl_dataforpsignifit_jk{partcount}.rhysubj(stimlev,3)=gl_dataforpsignifit_jk{partcount}.rhyobj(stimlev,3);
    end

    % Mean Objective

    for idxcond=1:3
        means.objmeancondcontrast_jk{partcount}(idxcond,:)=mean(fullresultsobj{idxcond},1,'omitnan');
    end

    %%%% Calculate Subjective Performance for each Contrast and Condition

    fullresultssubj=[];

    % Make table samples x contrasts (subjective)
    for idxcond=1:3 % for each condition
        conditiondata=reduced_data(find(reduced_data{:,"Condition"}==idxcond),:); % get trials
        results=[];
        for idxcont=1:10 % find each contrast
            contrastindex=find(conditiondata{:,"Contrast Level"}==idxcont);
            results(:,idxcont)=conditiondata{contrastindex,"Binary Visibility"};
        end
        fullresultssubj{idxcond}=results; % calculate mean across trials for each contrast
    end

    % Mean Subjective
    means.subjmeancondcontrast=[];
    for idxcond=1:3
        means.subjmeancondcontrast_jk{partcount}(idxcond,:)=mean(fullresultssubj{idxcond},1,'omitnan');
    end

    % Save

    grouplevel_res_jk.obj{partcount}=means.objmeancondcontrast_jk;
    grouplevel_res_jk.subj{partcount}=means.subjmeancondcontrast_jk;
end

grouplevel_res_jk.means.obj=grouplevel_res_jk.obj;
grouplevel_res_jk.means.subj=grouplevel_res_jk.subj;
save('gl_dataforpsignifit_jk','gl_dataforpsignifit_jk')


%% Fit Curves

% Bring data into correct structure: levelxcorrectxtotal (done above)

%load('grouplevel_res','gl_dataforpsignifit_jk')

% Prepare struct with fitting options
fitting_options_obj=struct;
fitting_options_subj=struct;

fitting_options_obj.expType = 'nAFC';
fitting_options_obj.expN = 2;
fitting_options_subj.expType = 'YesNo';

fitting_options_obj.sigmoidName = 'logistic';
fitting_options_subj.sigmoidName = 'logistic';

for partcount=1:totalsubj

    % Fit Objective Curves

    grouplevel_res_jk_curve.irrobj{partcount} = psignifit(gl_dataforpsignifit_jk{partcount}.irrobj,fitting_options_obj);
    grouplevel_res_jk_curve.intobj{partcount} = psignifit(gl_dataforpsignifit_jk{partcount}.intobj,fitting_options_obj);
    grouplevel_res_jk_curve.rhyobj{partcount} = psignifit(gl_dataforpsignifit_jk{partcount}.rhyobj,fitting_options_obj);

    % Fit Subjective Curves
    grouplevel_res_jk_curve.irrsubj{partcount} = psignifit(gl_dataforpsignifit_jk{partcount}.irrsubj,fitting_options_subj);
    grouplevel_res_jk_curve.intsubj{partcount} = psignifit(gl_dataforpsignifit_jk{partcount}.intsubj,fitting_options_subj);
    grouplevel_res_jk_curve.rhysubj{partcount} = psignifit(gl_dataforpsignifit_jk{partcount}.rhysubj,fitting_options_subj);


    % Plot
    plotOptions1.lineColor = [0,0,0];
    plotOptions1.dataColor = [0,0,0];
    plotOptions1.CIthresh = true;
    plotOptions1.dataSize=2;
    plotOptions1.lineWidth = 1.5;
    plotOptions2.lineColor = [1,0,0];
    plotOptions2.dataColor = [1,0,0];
    plotOptions2.lineWidth = 1.5;
    plotOptions2.dataSize=2;
    plotOptions2.CIthresh = true;
    plotOptions3.lineColor = [0,0.7,0.5];
    plotOptions3.dataColor = [0,0.7,0.5];
    plotOptions3.dataSize=2;
    plotOptions3.lineWidth = 1.5;
    plotOptions3.CIthresh = true;

    figure(1);
    subplot(5,3,partcount); [hline]=plotPsych(grouplevel_res_jk_curve.irrobj{partcount},plotOptions1);
    hold on
    subplot(5,3,partcount); [hline2]=plotPsych(grouplevel_res_jk_curve.intobj{partcount},plotOptions2);
    subplot(5,3,partcount); [hline3]=plotPsych(grouplevel_res_jk_curve.rhyobj{partcount},plotOptions3);
    title1=sprintf('Objective Without Subj%i',partcount);
    title(title1)
    if partcount==totalsubj
    legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')
    end
    hold off

    figure(2);
    subplot(5,3,partcount); [hline]=plotPsych(grouplevel_res_jk_curve.irrsubj{partcount},plotOptions1);
    hold on
    subplot(5,3,partcount); [hline2]=plotPsych(grouplevel_res_jk_curve.intsubj{partcount},plotOptions2);
    subplot(5,3,partcount); [hline3]=plotPsych(grouplevel_res_jk_curve.rhysubj{partcount},plotOptions3);
    title2=sprintf('Subjective Without Subj%i',partcount);
    title(title2)
    if partcount==totalsubj
    legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')
    end
    hold off
end
save('grouplevel_res_jk','grouplevel_res_jk')

% Extract new midpoints
for partcount=1:totalsubj
    grouplevel_res_jk.thresholdsubj{partcount,3}=grouplevel_res_jk_curve.irrsubj{1,partcount}.Fit(1);
    grouplevel_res_jk.thresholdsubj{partcount,2}=grouplevel_res_jk_curve.intsubj{1,partcount}.Fit(1);
    grouplevel_res_jk.thresholdsubj{partcount,1}=grouplevel_res_jk_curve.rhysubj{1,partcount}.Fit(1);

    grouplevel_res_jk.thresholdobj{partcount,3}=grouplevel_res_jk_curve.irrobj{1,partcount}.Fit(1);
    grouplevel_res_jk.thresholdobj{partcount,2}=grouplevel_res_jk_curve.intobj{1,partcount}.Fit(1);
    grouplevel_res_jk.thresholdobj{partcount,1}=grouplevel_res_jk_curve.rhyobj{1,partcount}.Fit(1);
end

save('grouplevel_res_jk','grouplevel_res_jk')


%% Jackknife Stats: Repeated measure one way ANOVA on threshold values
load('grouplevel_res_jk.mat')

% Define Within-Subj Variable
Condition = [1 2 3]'; % 1=Rhythm, 2=Interval, 3=Irregular

obj_mat=cell2mat(grouplevel_res_jk.thresholdobj);
subj_mat=cell2mat(grouplevel_res_jk.thresholdsubj);

% Convert into table
results_jk_table_obj = table(obj_mat(:,1),obj_mat(:,2),obj_mat(:,3), ...
    'VariableNames',{'t1','t2','t3'});
results_jk_table_subj = table(subj_mat(:,1),subj_mat(:,2),subj_mat(:,3), ...
    'VariableNames',{'t1','t2','t3'});

rm_obj = fitrm(results_jk_table_obj,'t1-t3 ~ 1','WithinDesign',Condition);
ranovatbl_obj = ranova(rm_obj);
rm_subj = fitrm(results_jk_table_subj,'t1-t3 ~ 1','WithinDesign',Condition);
ranovatbl_subj = ranova(rm_subj);

save('jk_stats','ranovatbl_subj',"ranovatbl_obj","rm_subj","rm_obj",'obj_mat','subj_mat')

% Correct F values
F_corr_obj=ranovatbl_obj.F(1)/((totalsubj-1)^2); % fcorr=F/(n-1)^2 with n= number of subsamples
F_corr_subj=ranovatbl_subj.F(1)/(totalsubj-1)^2;

% % Post-Hoc
% multcompare(rm_obj,'t1-t3')

%% Plot Jackknifing Means and Data Points
figure;
plot(mean(obj_mat(:,1)),5,'.r','MarkerSize',30)
hold on
plot(mean(obj_mat(:,2)),5,'.g','MarkerSize',30)
plot(mean(obj_mat(:,3)),5,'.b','MarkerSize',30)

plot(obj_mat(:,1),5,'.r','MarkerSize',10)
plot(obj_mat(:,2),5,'.g','MarkerSize',10)
plot(obj_mat(:,3),5,'.b','MarkerSize',10)
title('Objective Jackknifing Thresholds')
legend('Rhythm','Interval','Irregular')

figure;
plot(mean(subj_mat(:,1)),5,'.r','MarkerSize',30)
hold on
plot(mean(subj_mat(:,2)),5,'.g','MarkerSize',30)
plot(mean(subj_mat(:,3)),5,'.b','MarkerSize',30)

plot(subj_mat(:,1),5,'.r','MarkerSize',10)
plot(subj_mat(:,2),5,'.g','MarkerSize',10)
plot(subj_mat(:,3),5,'.b','MarkerSize',10)

title('Subjective Jackknifing Thresholds')
legend('Rhythm','Interval','Irregular')


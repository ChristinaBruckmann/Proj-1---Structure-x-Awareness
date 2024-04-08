% f(x) = L / 1 + e ^ -k(x-x0)
%
% x0, the x value of the sigmoid's midpoint;
% L, the curve's maximum value;
% k, the logistic growth rate or steepness of the curve

%% Test Staircase
clearvars

numberofparticipants=25;
numberoftrialsperpart=100;
x0=rand(numberofparticipants,1); %simulated PSE
x=0:.001:1; % levels of x
npart = 1:numberofparticipants;
L = 1;
kpart = randi(30,numberofparticipants,1); %simulated slope

%startvalue = x0 * 1.5;
%startvalue(startvalue>1)=1;
startvalue = 1;

for nfactor=1

    stepsize=0.1; %constant
    % startvalue=ones(10,1);
    % xvalues=startvalue;
    % UpDownChanges=1;
    Response=[];
    %trialsinarow=15;

    for n = npart
        xvalues=startvalue(n);
        ReductionFactor = stepsize;
        Response=[];
        similarxvalues=0;

        for nn=1:length(x)
            y(n,nn) = 1 / ( 1 + exp( -kpart(n) * (x(nn)-x0(n))));
        end

        for nn = 1 : numberoftrialsperpart

            pResponse=interp1(x, ...
                y(n,:),xvalues(nn),'nearest');

            %         end

            temp=rand(1)+ 0.1 * wgn( 1,1, 0 );

            if pResponse > temp
                Response(nn) = 1; %response towards the Up
                xvalues(nn+1)= xvalues(nn) - ReductionFactor;
                if xvalues(nn+1)<0
                    xvalues(nn+1)=0;
                end

            elseif pResponse <= temp
                Response(nn) = 0;
                xvalues(nn+1)= xvalues(nn) + 3*ReductionFactor;
                if xvalues(nn+1)>1
                    xvalues(nn+1)=1;
                end
            end

            %             if nn >= 2 && Response(nn) ~= Response(nn-1)
            %                 UpDownChanges = UpDownChanges + 2;
            %             end

        end


        [~, indx(n,1)] = min(abs(0.75-y(n,:)));
        upConv(n,1)=x(indx(n));


        ResultsStairCaseTest{n,1}=xvalues;
        ResultsStairCaseTest{n,2}=upConv(n);
        ResultsStairCaseTest{n,3}=y(n,:);

    end
end

%%
figure('Name','Simulated psychometric functions');
for n=npart
    plot(x,y(n,:));
    hold on
end



figure

if numberofparticipants<=100
    for n=1:numberofparticipants

        subplot(round(numberofparticipants/5),round(numberofparticipants/5),n), plot(ResultsStairCaseTest{n,1},'x-')
        axis([0 100 0 1])
        hold on
        line([0 100],[x0(n) x0(n)],'Color','red')
        title(num2str(x0(n)))

    end
else
    for n=1:100

        subplot(round(100/5),round(100/20),n), plot(ResultsStairCaseTest{n,1},'x-')
        axis([0 100 0 1])
        hold on
        line([0 100],[x0(n) x0(n)],'Color','red')
        title(num2str(x0(n)))

    end
end

figure

if numberofparticipants<=100
    for n=1:size(y,1)

        subplot(round(numberofparticipants/5),round(numberofparticipants/20),n), plot(x,y(n,:),'-')
        axis([0 1 0 1])
        hold on
        line([x0(n) x0(n)],[0 1],'Color','red')

    end
else
    for n=1:100

        subplot(round(100/5),round(100/20),n), plot(x,y(n,:),'-')
        axis([0 1 0 1])
        hold on
        line([x0(n) x0(n)],[0 1],'Color','red')

    end

end
for n=1:numberofparticipants

    temp=ResultsStairCaseTest{n,1}<=x0(n)*1.05 & ResultsStairCaseTest{n,1}>=x0(n)*.95;

    for nn=1:numberoftrialsperpart

        if nn<numberoftrialsperpart-trialsinarow

            if sum(temp(nn:nn+trialsinarow))==trialsinarow
                ResultsStairCaseTest{n,4}=nn+trialsinarow;
                break
            end
        else
            ResultsStairCaseTest{n,4}=nan;
        end
    end

end

Results(nfactor,:)=[mean([ResultsStairCaseTest{:,4}],'omitnan'),std([ResultsStairCaseTest{:,4}],'omitnan'),...
    sum(isnan([ResultsStairCaseTest{:,4}]))/10];
%             similarxvalues=similarxvalues+1;
%         elseif xvalues(nn)>=x0(n)*.95
%             similarxvalues=similarxvalues+1;

end

figure, plot(Results(:,1),'x-')
hold on
plot(Results(:,3),'rx-')

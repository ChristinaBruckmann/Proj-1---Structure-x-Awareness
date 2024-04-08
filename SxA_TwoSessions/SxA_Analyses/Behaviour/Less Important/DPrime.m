 % calculate dprime
                hit=sum(condData(:,1)==1 & condData(:,2)==1)/sum(condData(:,1)==1);
                FA=sum(condData(:,1)==0 & condData(:,2)==1)/sum(condData(:,1)==0);
 
                if hit==1 % correction for extreme values
                    hit=(sum(condData(:,1)==1)-0.5)/sum(condData(:,1)==1);
                elseif hit==0
                    hit=0.5/sum(condData(:,1)==1);
                end
                
                if FA==1
                    FA=(sum(condData(:,1)==0)-0.5)/sum(condData(:,1)==0);
                elseif FA==0
                    FA=0.5/sum(condData(:,1)==0);
                end
                
                dP=icdf('Normal', hit, 0,1)-icdf('Normal', FA, 0,1); % dprime
                subDataDP=[subDataDP, dP];


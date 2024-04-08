clear
clc
% Tet
trials=48; % Amount of trials per block
blocks=10; % total blocks
participants=10;

for part=1:participants
    for block=1:blocks
        % 75% correct on average:

        % 50 percent are chance
        for t=1:trials/2
            data(t)=randi(2,1)-1;
        end

        % the other 50 are always correct
        data(trials/2+1:trials)=1;

        % calculate percent correct
        pc(block)=mean(data);
    end
    pcmean(part)=mean(pc);
    subplot(2,5,part); plot(1:blocks, pc)
end
% Plot  block performance
figure;plot(1:participants, pcmean)
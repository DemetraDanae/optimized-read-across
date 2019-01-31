%% Tropsha's acceptability criteria

function [verdict, R2] = tropsha(final)

    pred_av = mean(final(:,2));
    exp_av = mean(final(:,1));
    
    %R^2 calculation
    cor = corrcoef(final(:,2),final(:,1));
    R2 = cor(1,2)^2;
    %R2 = sum((final(:,1)-exp_av).*(final(:,2)-pred_av))/sqrt(sum((final(:,1)-exp_av).^2)*sum((final(:,2)-pred_av).^2));

    %k calculation
    k_ = sum((final(:,2)-exp_av).*(final(:,1)-pred_av))/sum(final(:,2).^2);
    k = sum(final(:,2).*final(:,1))/sum(final(:,1).^2);

    % R0^2 calculation
    yro_ = k*final(:,1);
    yro = k_*final(:,2);

    R02 = 1-(sum((final(:,1)-yro(:)).^2)/sum((final(:,1)-exp_av).^2));
    R02_ = 1-(sum((final(:,2)-yro_(:)).^2)/sum((final(:,2)-pred_av).^2));

    %q2
    q2 = 1-sum((final(:,1)-final(:,2)).^2)/sum((final(:,1)-exp_av).^2);
   
    %Accurate predictions?
    if (R2>0.6) && (q2>0.5) && (abs(R2-R02)<0.3) && ((((R2-R02)/R2<0.1) && (0.85<=k<=1.15))|(((R2-R02_)/R2<0.1) && (0.85<=k_<=1.15)))
        verdict = 'pass';
    else
        verdict = 'fail';
    end

end
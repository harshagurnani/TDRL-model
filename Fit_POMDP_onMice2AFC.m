
function [Val_2_min] = Fit_POMDP_onMice2AFC(x,Data)

data=Data.data;
AnimalID=Data.ID;

alphaSame =x(1);  % learning rate
alphaDiff =x(2);  % learning rate
da_val    =x(3);  % value of DA
noiseSTD  =x(4);  % sensory noise
beta      =x(5);  % softmax param

alpha_L = alphaSame - alphaDiff;   % if alphaDiff = 0, then alpha_L an R are the same. useful for later striatal data
alpha_R = alphaSame + alphaDiff;


% outcomes of correct trials:
reward = [1+da_val, 1, 1+da_val, 1;
    1, 1+da_val, 1+da_val, 1] ;

contrast =unique(data(:,2));

contrastNum = length(unique(data(:,2)));

trialN  =   length(data);

% initialise
QLL    = [1];
QRR    = [1];
QLR    = [1];
QRL    = [1];
QL     = nan (trialN,1);
QR     = nan (trialN,1);
action = nan (trialN,1);
Phat   = nan (trialN,2);
P      = nan (trialN,1);

for trials = 1 : trialN
    
    if data (trials, 11) ==1 && data (trials, 1) ==1   % at the onset of the day initialise Q values..but don't initilaise if sessions ran in one day
        
        QLL = [1];
        QRR = [1];
        QLR = [1];
        QRL = [1];
        
    end
    
    currentContrast = data (trials, 2);
    
    % sensory noise
    contrast_withnoise        = currentContrast + noiseSTD * randn;
    
    Belief = normpdf(contrast,contrast_withnoise,noiseSTD);
    Belief = Belief ./sum(Belief);            % normalise belief
    Belief_L=sum(Belief(1:floor(contrastNum/2)));
    Belief_R=sum(Belief(floor(contrastNum/2)+1:end));
    
    QL(trials) = Belief_L*QLL + Belief_R*QRL;
    QR(trials) = Belief_L*QLR + Belief_R*QRR;
    
    % generate prob (softmax)
    
    Phat(trials,1) = exp(beta * QL(trials)) ./ (exp(beta * QL(trials)) + exp(beta * QR(trials)) );
    Phat(trials,2)= exp(beta * QR(trials)) ./ (exp(beta * QL(trials)) + exp(beta * QR(trials)) );
    
   
    % action based on animal (when we fit, it is animal's choice to be considered)
    action(trials) = data(trials,3);
    
    % trial reward
    
    if currentContrast<0 && action(trials)==-1
        Reward = reward(1,data(trials,8)) ;
        
    elseif currentContrast>0 && action(trials)==1
        Reward = reward(2,data(trials,8)) ;
    elseif currentContrast==0
        if action(trials)==-1
            
            Reward = reward(1,data(trials,8)) ;
        elseif action(trials)==1
            
            Reward = reward(2,data(trials,8)) ;
        end
        
    else
        Reward = 0 ;
    end
    
    if action(trials) == -1  % left
        
        delta        = Reward - QL(trials);
        
        QLL = QLL + alpha_L*delta*Belief_L ;
        QRL = QRL + alpha_R*delta*Belief_R ;
        
    else   % right
        
        
        delta        = Reward - QR(trials);
        
        QLR = QLR + alpha_L*delta*Belief_L ;
        QRR = QRR + alpha_R*delta*Belief_R;
        
    end
    
end

action(action==1)=2;
action(action==-1)=1;

for i=1:trials
    P(i,:)=Phat(i,action(i));
end


Val_2_min= -sum(log2(P));


function [ TrialStimuli, TrialBlocks, action, correct, QL, QR] = DummyRun_IdenticalBelief(x, data)


TrialStimuli = data(:,2);
TrialBlocks = data(:,8);

alpha     =x(1);  % learning rate
da_val    =x(2);  % value of DA
noiseSTD  =x(3);  % sensory noise
beta      =x(4);  % softmax param

alpha_L = alpha;   % if alphaDiff = 0, then alpha_L an R are the same. useful for later striatal data
alpha_R = alpha;


% outcomes of correct trials:
reward = [1+da_val, 1, 1+da_val, 1;
          1, 1+da_val, 1+da_val, 1] ;

contrast =unique(data(:,2));
Belief_PC = erf(contrast/noiseSTD) ;

% contrastNum = length(unique(data(:,2)));

trialN  =   length(data);
%% initialise
% set run numbers (odd)
iterN = 21;                 
% if nargin>3
%     iterN = 2*floor(varargin{1}/2)+1;
% end
% trialN = length(contrastTrials);

% initialise variables, for speed
action  = nan(trialN,iterN);
correct = nan(trialN, iterN);
% correct = nan(trialN,iterN);
QL      = nan(trialN,iterN);
QR      = nan(trialN,iterN);




Phat   = nan (trialN,2);
P      = nan (trialN,1);

for iter = 1:iterN
    % initialise
    QLL    = 1;
    QRR    = 1;
    QLR    = 1;
    QRL    = 1;
for trials = 1 : trialN
    
    if data (trials, 11) ==1 && data (trials, 1) ==1   % at the onset of the day initialise Q values..but don't initilaise if sessions ran in one day
        
        QLL = 1;
        QRR = 1;
        QLR = 1;
        QRL = 1;
        
    end
    
    currentContrast = TrialStimuli(trials);
    
    % sensory noise
%     contrast_withnoise        = currentContrast + noiseSTD * randn;
    
%     Belief = normpdf(contrast,contrast_withnoise,noiseSTD);
%     Belief = Belief ./sum(Belief);            % normalise belief
%     Belief_L=sum(Belief(1:floor(contrastNum/2)));
%     Belief_R=sum(Belief(floor(contrastNum/2)+1:end));
    Belief_L = Belief_PC( contrast==currentContrast );
    Belief_R = 1-Belief_L;
    QL(trials, iter) = Belief_L*QLL + Belief_R*QRL;
    QR(trials, iter) = Belief_L*QLR + Belief_R*QRR;
    
    % generate prob (softmax)
    
    Phat(trials,1) = exp(beta * QL(trials, iter)) ./ (exp(beta * QL(trials, iter)) + exp(beta * QR(trials, iter)) );
    Phat(trials,2)= exp(beta * QR(trials, iter)) ./ (exp(beta * QL(trials, iter)) + exp(beta * QR(trials, iter)) );
    
   
    % action based on animal (when we fit, it is animal's choice to be considered)
    if rand<Phat(trials,1)
        action(trials, iter) = -1;
    else
        action(trials, iter) = 1;
    end
    
       
%    action(trials) = data(trials,3);
    
    % trial reward
    
    if currentContrast<0 && action(trials, iter)==-1
        Reward = reward(1,TrialBlocks(trials)) ;
        correct(trials, iter) = 1;
    elseif currentContrast>0 && action(trials, iter)==1
        Reward = reward(2,TrialBlocks(trials)) ;
        correct(trials, iter) = 1;
    elseif currentContrast==0
        if action(trials, iter)==-1
            
            Reward = reward(1,TrialBlocks(trials)) ;
        elseif action(trials, iter)==1
            
            Reward = reward(2,TrialBlocks(trials)) ;
        end
        correct(trials, iter) = 1;
    else
        Reward = 0 ;
        correct(trials, iter) = 0;
    end
    
    if action(trials, iter) == -1  % left
        
        delta        = Reward - QL(trials, iter);
        
        QLL = QLL + alpha_L*delta*Belief_L ;
        QRL = QRL + alpha_R*delta*Belief_R ;
        
    else   % right
        
        
        delta        = Reward - QR(trials, iter);
        
        QLR = QLR + alpha_L*delta*Belief_L ;
        QRR = QRR + alpha_R*delta*Belief_R;
        
    end
    
end
end
% 
% action(action==1)=2;
% action(action==-1)=1;

% for i=1:trials
%     P(i,:)=Phat(i,action(i));
% end
% 
% 
% NLL= -sum(log2(P));


function [NLL, NLL_L, NLL_R] = Get_NLL_ModelWithBias_TempTesting(params, Data, includeContrast, varargin)
%% 4 Parameter Model - Learning rate, Dopamine Value, Sensory Noise, Action Bias

%%% NEW CHANGES:
%%% changed calculation of PC - use mean across iterations rather than mode
%%% response.
%%% belief is calculated at trial stimuli, but prob density at zero
%%% contrast is shared 

data = Data.data;

alpha     = params(1);  % learning rate
da_val    = params(2);  % value of DA
noiseSTD  = params(3);  % sensory noise
Qbias     = params(4);  % bias in decision/value

% set alpha values
alpha_L = alpha;
alpha_R = alpha;


% outcomes of correct trials:
reward = [  1+da_val , 1         , 1+da_val, 1;
            1        , 1+da_val  , 1+da_val, 1] ;

% set contrast
contrastTrials = data(:,2);
contrast = unique(data(:,2))';
contrastNum = length(contrast);




%% initialise
% set run numbers (odd)
iterN = 1;                 
if nargin>3
    iterN = 2*floor(varargin{1}/2)+1;
end
trialN = length(contrastTrials);


% initialise variables, for speed
action  = nan(trialN,iterN);
correct = nan(trialN,iterN);
QL      = nan(trialN,iterN);
QR      = nan(trialN,iterN);

parfor iter = 1:iterN
   
   % initalise Q values for each iteration
   QLL(1,:) = 1;
   QRR(1,:) = 1;
   QLR(1,:) = 1;
   QRL(1,:) = 1;
   
   % start model
   for trials = 1 : trialN
      
      % at the onset of the day initialise Q values..but don't initialise
      % if sessions ran in one day
      if data(trials, 11) ==1 && data(trials, 1) ==1
         
         QLL = 1;
         QRR = 1;
         QLR = 1;
         QRL = 1;
         
      end
      
      % set contrast
      currentContrast = contrastTrials(trials);
      
      % add sensory noise
      contrast_withnoise  = currentContrast + noiseSTD * randn;
      
      % define Belief for every possible contrast value
      Belief = normpdf(contrast,contrast_withnoise,noiseSTD);
      Belief = Belief ./sum(Belief);            % normalise belief
      
      % define Belief_L and Belief_R
      % corresponds to the model's belief about side of stimulus
      % Belief_L + Belief_R = 1
      Belief_L=sum(Belief(contrast < 0))+Belief(contrast==0)/2;
      Belief_R=sum(Belief(contrast>0))  +Belief(contrast==0)/2;
      
      %initialise Q values for this iteration
      QL(trials,iter) = Belief_L*QLL + Belief_R*QRL;
      QR(trials,iter) = Belief_L*QLR + Belief_R*QRR;    
      
      
      % action based on model (action <-- max(QL,QR))
      if QL(trials,iter) > QR(trials,iter)+ Qbias           % Add fixed bias throughout - without changing prediction errors
         
         action(trials,iter) = -1;
         
      elseif QL(trials,iter) < QR(trials,iter)+ Qbias
         
         action(trials,iter) = 1;
         
      else
         if rand > 0.5
            action(trials,iter) = 1;
         else
            action(trials,iter) = -1;
         end
      end
      
      
      % trial reward based on model
      if currentContrast<0 && action(trials,iter)==-1
         
         correct(trials,iter) = 1;
         Reward = reward(1,data(trials,8));
         
      elseif currentContrast>0 && action(trials,iter)==1
         
         correct(trials,iter) = 1;
         Reward = reward(2,data(trials,8));
         
      elseif currentContrast==0
         
         if rand > 0.5
            correct(trials,iter) = 1; 
            if action(trials,iter)==-1
%                
%                correct(trials,iter) = 1; %-------------------- correct = -1???
               Reward = reward(1,data(trials,8));
               
            elseif action(trials,iter)==1
               
%                correct(trials,iter) = 1;
               Reward = reward(2,data(trials,8));
               
            end
            
         else
            
            correct(trials,iter) = 0;
            Reward = 0 ;
            
         end
         
      else
         
         correct(trials,iter) = 0;
         Reward = 0 ;
         
      end
      
      
      % calculate delta, and update Q values for the model
      if action(trials,iter) == -1  % left
         
         delta   = Reward - QL(trials,iter);
         
         QLL     = QLL + alpha_L*delta*Belief_L ;
         QRL     = QRL + alpha_R*delta*Belief_R ;
         
      else   % right
         
         delta   = Reward - QR(trials,iter);
         
         QLR     = QLR + alpha_L*delta*Belief_L;
         QRR     = QRR + alpha_R*delta*Belief_R;
         
      end
      
   end
   
end

modeAction = mode(action,2);

xx = contrast;
nn = nan(length(xx)*length(unique(data(:,8))),1);
pp = nn;
probs = nn;

nn_L=nn;        nn_R=nn;
pp_L=nn;        pp_R=nn;
probs_L=nn;     probs_R=nn;

c=1;
for j=unique(data(:,8))'
   for i=contrast
      
      
      nn(c) = length(data(data(:,2)==i & data(:,8)==j,2));                          % Total number of trials
      pp(c) = length(data(data(:,2)==i & data(:,8)==j & data(:,3)==1,2))/nn(c);     % Observed P(R)
      probs(c) = 0;
      for iter = 1:iterN
        probs(c) = probs(c) + length(data(data(:,2)==i & data(:,8)==j & action(:, iter)==1,2))/nn(c); % Estimated P(R) for each stimulus
      end
      probs(c) = probs(c)/iterN;
      
      %Trials after correct left choices
      id = [false; data(2:end,2)==i & data(2:end,8)==j                    & data(1:end-1,3)==-1  &  data(1:end-1,10)==1];
      nn_L(c) = sum(id); 
      id = [false;data(2:end,2)==i & data(2:end,8)==j & data(2:end,3)==1  & data(1:end-1,3)==-1  & data(1:end-1,10)==1];
      pp_L(c) = sum(id)/nn_L(c);            % Observed P(R)
      probs_L(c)=0;
      for iter = 1:iterN
        id = [false;data(2:end,2)==i & data(2:end,8)==j                           & action(1:end-1,iter)==-1  & correct(1:end-1,iter)==1];
        model_nn_L= sum(id);
        id = [false;data(2:end,2)==i & data(2:end,8)==j & action(2:end, iter)==1  & action(1:end-1,iter)==-1  & correct(1:end-1,iter)==1];
        probs_L(c) = probs_L(c) + sum(id)/model_nn_L; % Estimated P(R) 
      end
      probs_L(c) = probs_L(c)/iterN;        % Average across iterations
      
      
      %Trials after correct right choices
      id = [false;data(2:end,2)==i & data(2:end,8)==j                     & data(1:end-1,3)==1  & data(1:end-1,10)==1];
      nn_R(c) = sum(id);
      id = [false;data(2:end,2)==i & data(2:end,8)==j & data(2:end,3)==1  & data(1:end-1,3)==1  & data(1:end-1,10)==1];
      pp_R(c) = sum(id)/nn_R(c);            % Observed P(R)
      probs_R(c)=0;
      for iter = 1:iterN
        id = [false;data(2:end,2)==i & data(2:end,8)==j                           & action(1:end-1,iter)==1  & correct(1:end-1,iter)==1];
        model_nn_R= sum(id);
        id = [false;data(2:end,2)==i & data(2:end,8)==j & action(2:end, iter)==1  & action(1:end-1,iter)==1  & correct(1:end-1,iter)==1];
        probs_R(c) = probs_R(c) +sum(id)/model_nn_R; % Estimated P(R) 
      end
      probs_R(c) = probs_R(c)/iterN;        % Average across iterations
      
      c=c+1;

   end
   
   % if probs is 1 or 0, NLL will output infinity
   probs(probs==0)=eps;
   probs(probs==1)=1-eps;
   
   probs_L(probs_L==0)=eps;
   probs_L(probs_L==1)=1-eps;
   
   probs_R(probs_R==0)=eps;
   probs_R(probs_R==1)=1-eps;

end
contrastToUse = (includeContrast==1);
% Note: HG. Assuming Binomial RV of Rightward(R) Actions at each
% contrast(per block) 'c', Probability is given as:
% P(c) = choose(nn,k) * probs^(k) * (1-probs)^(nn-k)       
%           where  pp(c) = k(c)/nn(c)   or k is the number of R actions for contrast c  
% Total probability is  product of these probabilities for all contrasts.
% NLL = - Sum(log(P(c))
% - log(P(c)) = - log(choose(nn,k)) -  nn*[ k/nn log(probs) + (nn-k)/nn log(1-probs) ]
%           = CONSTANT K + nn*[ pp*log(probs) + (1-pp)*log(1-probs) ]
% Constant K is constant for all simulations so only need to maximise second term 
NLL   = - sum(nn(contrastToUse).*(pp(contrastToUse).*log(probs(contrastToUse))+(1-pp(contrastToUse)).*log(1-probs(contrastToUse))));
NLL_L = - sum(nn_L(contrastToUse).*(pp_L(contrastToUse).*log(probs_L(contrastToUse))+(1-pp_L(contrastToUse)).*log(1-probs_L(contrastToUse))));
NLL_R = - sum(nn_R(contrastToUse).*(pp_R(contrastToUse).*log(probs_R(contrastToUse))+(1-pp_R(contrastToUse)).*log(1-probs_R(contrastToUse))));

end


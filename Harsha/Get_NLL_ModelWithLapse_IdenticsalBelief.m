function [NLL] = Get_NLL_ModelWithLapse_IdenticsalBelief(params, Data, includeContrast)
%% 4 Parameter Model - Learning rate, Dopamine Value, Sensory Noise, Action Bias

%%% NEW CHANGES:
%%% changed calculation of PC - use mean across iterations rather than mode
%%% response.
%%% belief is calculated at a grid of stimuli, not necessarily the ones seen by the animal, 
%%% but prob density at zero contrast is shared between Left and Right

data = Data.data;

alpha     = params(1);  % learning rate
da_val    = params(2);  % value of DA
noiseSTD  = params(3);  % sensory noise
Qbias     = params(4);  % bias in decision/value
lambdaL   = params(5);  % probability of lapse to left
lambdaR   = params(6);  % probability of lapse to right
beta      = params(7);

% set alpha values
alpha_L = alpha;
alpha_R = alpha;
netlapse = lambdaL+lambdaR;

% outcomes of correct trials:
reward = [  1+da_val , 1         , 1+da_val, 1;
            1        , 1+da_val  , 1+da_val, 1] ;

% set contrast
contrastTrials = data(:,2);
contrast = unique(data(:,2))';
contrastNum = length(contrast);

% define Belief for every possible contrast value
Belief_PC = erf(contrast/noiseSTD) ;





%% initialise
% set run numbers (odd)
iterN = 1;                 
% if nargin>3
%     iterN = 2*floor(varargin{1}/2)+1;
% end
trialN = length(contrastTrials);
Phat   = nan (trialN,2);
P      = nan (trialN,1);


% initialise variables, for speed
action  = nan(trialN,iterN);
correct = nan(trialN,iterN);
QL      = nan(trialN,iterN);
QR      = nan(trialN,iterN);
%onluy 1 iteration
% for iter = 1:iterN
   iter = 1;
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
      currentContrast = data(trials,2);

      % Belief at contrast
      Belief_L = Belief_PC( contrast==currentContrast );
      Belief_R = 1-Belief_L;
      
      %initialise Q values for this iteration
      QL(trials, iter) = Belief_L*QLL + Belief_R*QRL;
      QR(trials, iter) = Belief_L*QLR + Belief_R*QRR;
    
      
      % generate prob (softmax)    
      Phat(trials,1) = lambdaL + (1-netlapse)* exp(beta * QL(trials, iter)) ./(exp(beta * QL(trials, iter)) + exp(beta * QR(trials, iter)) );
      Phat(trials,2) = lambdaR + (1-netlapse)* exp(beta * QR(trials, iter)) ./(exp(beta * QL(trials, iter)) + exp(beta * QR(trials, iter)) );
    
      action(trials, iter) = data(trials, 3);
      
      % trial reward based on animal's choice
      if currentContrast<0 && action(trials,iter)==-1
         
         correct(trials,iter) = 1;
         Reward = reward(1,data(trials,8));
         
      elseif currentContrast>0 && action(trials,iter)==1
         
         correct(trials,iter) = 1;
         Reward = reward(2,data(trials,8));
         
      elseif currentContrast==0         
         if data(trials,10)==1
            correct(trials,iter) = 1;             
            if action(trials,iter)==-1                
               Reward = reward(1,data(trials,8));              
            elseif action(trials,iter)==1               
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
   
% end

action(action==1)=2;
action(action==-1)=1;

for i=1:trials
    P(i,:)=Phat(i,action(i));
end


NLL= -sum(log2(P));

end


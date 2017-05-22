
function [data_model, action, correct] = RunPOMDP_GS_NLL(Data,params)

data = Data.data;

alpha = params(1);
DA_val = params(2);
noiseSTD = params(3);

contrastTrials = data(:,2);

% outcomes of correct trials:
reward = [1+DA_val, 1, 1+DA_val, 1;
          1, 1+DA_val, 1+DA_val, 1] ;


% set run numbers
iterN = 99;                    % model values are averaged over iterations
trialN = length(contrastTrials);


% set alpha values
alpha_L = alpha;
alpha_R = alpha;


% set contrast
contrast = unique(data(:,2))';
contrastNum = length(contrast);


% initialise variables, for speed
action = nan(trialN,iterN);
correct = zeros(trialN,iterN);
QL = nan (trialN,iterN);
QR = nan (trialN,iterN);

BeepValue = nan (trialN,iterN);
StimPE    = nan (trialN,iterN);
delta    = nan (trialN,iterN);
 
parfor iter = 1:iterN
   
   % initalise Q values for each iteration
   QLL(1,:) = 1;
   QRR(1,:) = 1;
   QLR(1,:) = 1;
   QRL(1,:) = 1;
   
   % start model
   for trials = 1 : trialN
      
      % at the onset of the day initialise Q values..but don't initilaise
      % if sessions ran in one day
      if data (trials, 11) ==1 && data (trials, 1) ==1
         
         QLL = 1;
         QRR = 1;
         QLR = 1;
         QRL = 1;
         
      end
      
      BeepValue(trials,iter) = (QLL + QLR + QRL + QRR ) ./ 4;

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
      Belief_L=sum(Belief(1:floor(contrastNum/2)));
      Belief_R=sum(Belief(floor(contrastNum/2)+1:end));
      
      %initialise Q values for this iteration
      QL(trials,iter) = Belief_L*QLL + Belief_R*QRL  ;
      QR(trials,iter) = Belief_L*QLR + Belief_R*QRR  ;
      
      
      % action based on model (action <-- max(QL,QR))
      if QL(trials,iter) > QR(trials,iter)
         
         action(trials,iter) = -1;
         StimPE(trials,iter) = QL(trials,iter) - BeepValue(trials,iter); 

         
      elseif QL(trials,iter) < QR(trials,iter)
         
         action(trials,iter) = 1;
         StimPE(trials,iter) = QR(trials,iter) - BeepValue(trials,iter); 

         
      else
         if rand > 0.5
            action(trials,iter) = 1;
            StimPE(trials,iter) = QR(trials,iter) - BeepValue(trials,iter); 

         else
            action(trials,iter) = -1;
             StimPE(trials,iter) = QL(trials,iter) - BeepValue(trials,iter); 

         end
      end
      

      
      
      % trial reward based on model
      if currentContrast<0 && action(trials,iter)==-1
         
         Reward = reward(1,data(trials,8));
         correct(trials,iter)=1;
      elseif currentContrast>0 && action(trials,iter)==1
         
         Reward = reward(2,data(trials,8));
         correct(trials,iter)=1;
      elseif currentContrast==0
         
         if rand > 0.5
            correct(trials,iter)=1;
            if action(trials,iter)==-1
               
               Reward = reward(1,data(trials,8));
               
            elseif action(trials,iter)==1
               
               Reward = reward(2,data(trials,8));
               
            end
            
         else
            
            Reward = 0 ;
            
         end
         
      else
         
         Reward = 0 ;
         
      end
      
      
      % calculate delta, and update Q values for the model
      if action(trials,iter) == -1  % left
         
         delta(trials, iter)   = Reward - QL(trials,iter);
         
         QLL     = QLL + alpha_L*delta(trials, iter)*Belief_L ;
         QRL     = QRL + alpha_R*delta(trials, iter)*Belief_R ;
         
      else   % right
         
         delta(trials, iter)   = Reward - QR(trials,iter);
         
         QLR     = QLR + alpha_L*delta(trials, iter)*Belief_L;
         QRR     = QRR + alpha_R*delta(trials, iter)*Belief_R;
         
      end
      
   end
   
end

% set output
% Note: HG. Maybe turn this inyo a structure to allow for wasier parsing of
% useful data
data_model = [data,  mean(action,2), mean(QL,2), mean(QR,2),...
   mean(BeepValue,2), mean(StimPE,2), mean(delta,2), mode(action,2)];

end

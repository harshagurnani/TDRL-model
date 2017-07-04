function [ TrialStimuli, TrialBlocks, action, correct, QL, QR ] = DummyRun_WithLapse( params )
% Run model with parameters for fake stimuli
% shifts after diff/easy stimuli

% PARAMS
alpha       = params(1);
DA_val      = params(2);
noiseSTD    = params(3);
Qbias       = params(4);
lambda      = params(5);
lapseAct    = params(6);       %-1 or 1
if (lapseAct ~= -1 && lapseAct ~= 1)
    lapseAct = 1;
    warning('Parameter 6: Lapse Action must be 1 or -1. Unaccepted value provided - defaulting to 1 (Right action)')
end

stimuli = [-0.5 -0.25 -0.12 -0.05 0 0.05 0.12 0.25 0.5]';
nStimuli = length(stimuli);
nTrials = 30000;
nIters = 99;

%Generate random stims
TrialStimuli = stimuli( unidrnd( nStimuli, [nTrials,1]) );

blocks = [1,2, 3, 4];
nBlocks = length(blocks);

%Create blocks
TrialBlocks = ones(nTrials,1);
TrialsPerSess = 500;
nSess = floor(nTrials/(TrialsPerSess * nBlocks));
for jj=1:nSess
   blockID = blocks(1+mod(jj,nBlocks));
   TrialBlocks(1+(jj-1)*TrialsPerSess:jj*TrialsPerSess) = blockID;
end

% reward 
BlockReward = [1+DA_val,     1   , 1+DA_val, 1;
                   1   , 1+DA_val, 1+DA_val, 1 ];
        
% New days
action = nan( nTrials, nIters);
correct = zeros(nTrials, nIters);
newDay = zeros( nTrials,1);
for jj=501:nTrials
   if sum(newDay(jj-500:jj-1)) == 0
       if rand<0.0025
          newDay(jj) = 1;  
       end
   end
end

Bstimuli = [min(stimuli):0.02:max(stimuli)];

Bstimuli = [-1:0.02:1];
Bstimuli  = stimuli;
QL = nan(nTrials, nIters); QR = QL;

% Run model
for iter = 1:nIters
   QLL = 1; QLR = 1; QRL = 1; QRR = 1;
   for trial = 1:nTrials
      if newDay(trial) == 1
          QLL = 1; QLR = 1; QRL = 1; QRR = 1;
      end
      BlockID = TrialBlocks(trial);
      lapse_rate = lambda;
      if BlockID == 3
          lapse_rate = 1.5*lambda;
      elseif BlockID == 4
          lapse_rate = 0;
      end
      reward = BlockReward(:, BlockID);
      
      % default - no reward unless correct action
      rewardGiven = 0;
      
      %Internal noisy representation of stimulus
      stimResp = randn*noiseSTD + TrialStimuli(trial);
      
      %Belief state
      Belief = normpdf(Bstimuli,stimResp,noiseSTD);
      Belief = Belief./sum(Belief);
      
      Belief_L = sum(Belief(Bstimuli<0)) + 0.5*Belief(Bstimuli==0);
      Belief_R = sum(Belief(Bstimuli>0)) + 0.5*Belief(Bstimuli==0);
      
      % Q-values "Expected reward"
      QL(trial, iter) = Belief_L * QLL + Belief_R * QRL;
      QR(trial, iter) = Belief_L * QLR + Belief_R * QRR;
      
      % Action
      if rand < lapse_rate
          % Lapse 
          action(trial, iter) = lapseAct;
      else
          
          if QL(trial, iter) > QR(trial, iter) + Qbias
              % left more valuable
              action(trial, iter) = -1;

          elseif QL(trial, iter) < QR(trial, iter) + Qbias
              % right more valuable
              action(trial, iter) = 1;

          else
            % Random action
            if rand< 0.5
                action(trial, iter) = -1;
            else 
                action(trial, iter) = 1;
            end
          end

      end
      
      % correct or error
      if TrialStimuli(trial) < 0  &&  action(trial, iter) == -1
        %Correct left action
        correct( trial, iter) = 1;
        rewardGiven = reward(1);
            
         
      elseif TrialStimuli(trial) > 0 && action(trial, iter) == 1
        % Correct Right action
        correct(trial, iter) = 1;
        rewardGiven = reward(2);
        
      elseif TrialStimuli(trial) == 0
         if rand() < 0.5
             % Correct randomly
             rewardGiven =  reward( 1+ (action(trial,iter)==1) );
             correct( trial,iter) = 1;
         end
      end
    
      if action(trial, iter) == -1
        PredError = rewardGiven - QL(trial, iter);
        QLL = QLL + alpha*PredError*Belief_L;
        QRL = QRL + alpha*PredError*Belief_R;
      else
        PredError = rewardGiven - QR(trial, iter);
        QLR = QLR + alpha*PredError*Belief_L;
        QRR = QRR + alpha*PredError*Belief_R;

      end
    
    
   end

end
end
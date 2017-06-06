function [ TrialStimuli, TrialBlocks, action, correct ] = DummyRun( params )
% Run model with parameters for fake stimuli
% shifts after diff/easy stimuli


alpha = params(1);
DA_val = params(2);
noiseSTD = params(3);
Qbias = params(4);


stimuli = [-0.5 -0.25 -0.12 -0.05 0 0.05 0.12 0.25 0.5];
nStimuli = length(stimuli);
nTrials = 10000;
nIters = 99;

%Generate random stims
TrialStimuli = stimuli( unidrnd( nStimuli, [nTrials,nIters]) );

blocks = [1,2];
nBlocks = length(blocks);

%Create blocks
TrialBlocks = ones(nTrials,1);
TrialsPerSess = 700;
nSess = floor(nTrials/(TrialsPerSess * nBlocks));
for jj=1:nSess
   blockID = blocks(1+mod(jj,nBlocks));
   TrialBlocks(1+(jj-1)*TrialsPerSess:jj*TrialsPerSess) = blockID;
end

% reward 
BlockReward = [1+DA_val,     1   , 1+DA_val, 1;
                   1   , 1+DA_val, 1+DA_val, 1 ];
        

action = nan( nTrials, nIters);
correct = zeros(nTrials, nIters);

QLL = 1; QLR = 1; QRL = 1; QRR = 1;
QL = nan(nTrials, nIters); QR = QL;

% Run model
for iter = 1:nIters
   for trial = 1:nTrials
      BlockID = TrialBlocks(trial);
      reward = BlockReward(:, BlockID);
      
      % default - no reward unless correct action
      rewardGiven = 0;
      
      %Internal noisy representation of stimulus
      stimResp = randn*noiseSTD + TrialStimuli(trial, iter);
      
      %Belief state
      Belief = normpdf(stimuli,stimResp,noiseSTD);
      Belief = Belief./sum(Belief);
      
      Belief_L = sum(Belief(stimuli<0));
      Belief_R = sum(Belief(stimuli>=0));
      
      % Q-values "Expected reward"
      QL(trial, iter) = Belief_L * QLL + Belief_R * QRL;
      QR(trial, iter) = Belief_L * QLR + Belief_R * QRR     + Qbias;
      
      % Action
      
      if QL(trial, iter) > QR(trial, iter)
          % left more valuable
          action(trial, iter) = -1;
          
      elseif QL(trial, iter) < QR(trial, iter)
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
      
      
      % correct or error
      if TrialStimuli(trial, iter) < 0  &&  action(trial, iter) == -1
        %Correct left action
        correct( trial, iter) = 1;
        rewardGiven = reward(1);
            
         
      elseif TrialStimuli(trial, iter) > 0 && action(trial, iter) == 1
        % Correct Right action
        correct(trial, iter) = 1;
        rewardGiven = reward(2);
        
      elseif TrialStimuli(trial, iter) == 0
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
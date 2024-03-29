function [ TrialStimuli, TrialBlocks, action, correct, QL, QR ] = DummyRun_noisyBelief( params )
% Run model with parameters for fake stimuli
% Belief becomes noisy on one side stimuli


alpha       = params(1);
DA_val      = params(2);
noiseSTD    = params(3);
Qbias       = params(4);

StimSide    = sign(params(5));
sigmaChange = abs(params(5));
if StimSide == 1
    noiseBL = noiseSTD;
    noiseBR = sigmaChange*noiseSTD;
else
    noiseBL = sigmaChange*noiseSTD;
    noiseBR = noiseSTD;
end
noise0 = sqrt( (noiseBL^2 + noiseBR^2)/2 );

stimuli     = [-0.5 -0.25 -0.12 -0.05 0 0.05 0.12 0.25 0.5]';
nStimuli    = length(stimuli);
nTrials     = 10000;
nIters      = 10;

%Generate random stims
TrialStimuli = stimuli( unidrnd( nStimuli, [nTrials,1]) );

blocks = [1,2];
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
        

action = nan( nTrials, nIters);
correct = zeros(nTrials, nIters);
newDay = zeros( nTrials,1);
for jj=2001:nTrials
   if sum(newDay(jj-1000:jj-1)) == 0
       if rand<0.003
          newDay(jj) = 1;  
       end
   end
end

Bstimuli = [min(stimuli):0.02:max(stimuli)];
Belief = nan(size(Bstimuli));
% Bstimuli = [-1:0.02:1];
% Bstimuli  = stimuli;
QL = nan(nTrials, nIters); QR = QL;
% Run model

nSteps = 1;     %no. of trials that stim lasts
for iter = 1:nIters
   QLL = 1; QLR = 1; QRL = 1; QRR = 1;
   noiseBL_trial = noiseSTD; noiseBR_trial = noiseSTD; noise0_trial = noiseSTD;
   laser = nSteps;
   reset = 1;
   for trial = 1:nTrials
      if newDay(trial) == 1
          QLL = 1; QLR = 1; QRL = 1; QRR = 1;
      end
      BlockID = TrialBlocks(trial);
      reward = BlockReward(:, BlockID);
%       laser = laser+1;
      if laser >= nSteps && reset==0
             noiseBL_trial = noiseSTD; noiseBR_trial = noiseSTD; noise0_trial = noiseSTD;
             reset = 1;
      end
      % default - no reward unless correct action
      rewardGiven = 0;
      
      %Internal noisy representation of stimulus
      stimResp = randn*noiseSTD + TrialStimuli(trial);
      
      %Belief state
      Belief(Bstimuli <0) = normpdf(Bstimuli(Bstimuli <0),stimResp,noiseBL_trial);
      Belief(Bstimuli >0) = normpdf(Bstimuli(Bstimuli >0),stimResp,noiseBR_trial);
      Belief(Bstimuli==0) = normpdf(0,stimResp,noise0_trial);
      Belief = Belief./sum(Belief);
      
      Belief_L = sum(Belief(Bstimuli<0)) + 0.5*Belief(Bstimuli==0);
      Belief_R = sum(Belief(Bstimuli>0)) + 0.5*Belief(Bstimuli==0);
      
      % Q-values "Expected reward"
      QL(trial, iter) = Belief_L * QLL + Belief_R * QRL;
      QR(trial, iter) = Belief_L * QLR + Belief_R * QRR;
      
      % Action
      
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
      
      
      % correct or error
      if TrialStimuli(trial) < 0  &&  action(trial, iter) == -1
        %Correct left action
        correct( trial, iter) = 1;
        rewardGiven = reward(1);
        if TrialBlocks(trial) == 1
            noiseBL_trial = noiseBL;
            noiseBR_trial = noiseBR;
            noise0_trial = noise0;
            laser = 0;    
            reset = 0;
        else
            laser = laser+1;
        end
      elseif TrialStimuli(trial) > 0 && action(trial, iter) == 1
        % Correct Right action
        correct(trial, iter) = 1;
        rewardGiven = reward(2);
        if TrialBlocks(trial) == 2
            noiseBL_trial = noiseBL;
            noiseBR_trial = noiseBR;
            noise0_trial = noise0;
            laser = 0;    
            reset = 0;
        else
            laser = laser+1;
        end
      elseif TrialStimuli(trial) == 0
         if rand() < 0.5
             % Correct randomly
             rewardGiven =  reward( 1+ (action(trial,iter)==1) );
             correct( trial,iter) = 1;
            if 2*(TrialBlocks(trial))-3 == action(trial,iter)
                noiseBL_trial = noiseBL;
                noiseBR_trial = noiseBR;
                noise0_trial = noise0;
                laser = 0;
                reset = 0;
            else
                laser = laser+1;
            end
         end
      else
          laser = laser+1;
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
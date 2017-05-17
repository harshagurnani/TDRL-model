function  data = LoadSavedDataForBehExp (AnimalID, ExpID)

if strcmp(ExpID,'1')   % loading Exp 1 and other useful sessions for that exp
    if strcmp(AnimalID,'SS030_DA')
        data = load('SS030_DA_Exp1_GoodSession1_Bait[]');
        data=data.data;
        
    elseif strcmp(AnimalID,'SS031_DA')
        data = load('SS031_DA_Exp1_GoodSession1_Bait[]');
        tested_stims=[-.5 -0.2 0 0.2 0.5];
        data=data.data;
        data(~ismember(data(:,2),tested_stims),:)=[];
         
    elseif strcmp(AnimalID,'SS037')
        data = load('SS037_Exp1_GoodSession1_Bait[]');
        %tested_stims=[-.5 -0.2 -0.05  0.05 0.2 0.5];
        data=data.data;     
        %data(~ismember(data(:,2),tested_stims),:)=[];
        
    elseif strcmp(AnimalID,'SS040')
        data = load('SS040_Exp1_GoodSession1_Bait[]');
        %tested_stims = [ -0.2  -0.1   -0.05 0.05  0.1   0.2 ];
        data=data.data;    
        %data(~ismember(data(:,2),tested_stims),:)=[];
        
    elseif strcmp(AnimalID,'EJ008DA')
        data = load('EJ008DA_Exp1_GoodSession1_Bait[]');
        %tested_stims = [ -0.5  -0.25   -0.12  0  0.12    0.25    0.5 ];
        data=data.data;      
        %data(~ismember(data(:,2),tested_stims),:)=[];
        
    elseif strcmp(AnimalID,'ALK05')
        data = load('ALK05_Exp1_GoodSession1_Bait[]');
        data=data.data;
        
    elseif strcmp(AnimalID,'ALK028')
        data = load('ALK028_Exp1_GoodSession[]_Bait[]');
        data=data.data;
        
    elseif strcmp(AnimalID,'ALK011')
        data = load('ALK011_Exp1_GoodSession[]_Bait[]');
        data=data.data;
        
    elseif strcmp(AnimalID,'ALK017')
        data = load('ALK017_Exp1_GoodSession[]_Bait[]');
        data=data.data;
        
    elseif strcmp(AnimalID,'ALK045')
        data = load('ALK045_Exp1_GoodSession[]_Bait[]');
        data=data.data;
    end
    
    
    
elseif strcmp(ExpID,'5') % loading Exp 5 and other useful sessions for that exp
    
    if strcmp(AnimalID,'SS031_DA')
        data = load('SS031_DA_Exp5_GoodSession[]_Bait[]');
        %data = load('SS031_DA_Exp5_GoodSession[]_Bait0');
        data=data.data;
    elseif strcmp(AnimalID,'SS037')
        data = load('SS037_Exp5_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'SS040')
        data = load('SS040_Exp5_GoodSession[]_Bait[]');
        %data = load('SS040_Exp5_GoodSession[]_Bait0');
        data=data.data;
    elseif strcmp(AnimalID,'ALK01')
        data = load('ALK01_Exp5_GoodSession[]_Bait[]');
        %data = load('ALK01_Exp5_GoodSession[]_Bait0');
        data=data.data;
    elseif strcmp(AnimalID,'ALK03')
        data = load('ALK03_Exp5_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK05')
        data = load('ALK05_Exp5_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK008')
        data = load('ALK008_Exp5_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK011')
        data = load('ALK011_Exp5_GoodSession[]_Bait[]');
        %data = load('ALK011_Exp5_GoodSession[]_Bait0');
        data=data.data;
    elseif strcmp(AnimalID,'EJ008DA')
        data = load('EJ008DA_Exp5_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK017')
        data = load('ALK017_Exp5_GoodSession[]_Bait[]');
        data=data.data;
    end
    
elseif strcmp(ExpID,'12') % loading Exp 12 and other useful sessions for that exp
    if strcmp(AnimalID,'ALK024')
        data = load('ALK024_Exp12_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK05')
        data = load('ALK05_Exp12_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK009')
        data = load('ALK009_Exp12_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK014')
        data = load('ALK014_Exp12_GoodSession[]_Bait[]');
        data=data.data;
    end
    
elseif strcmp(ExpID,'22') % loading Exp 22 and other useful sessions for that exp
    if strcmp(AnimalID,'ALK022')
        %        data = load('ALK022_Exp22_GoodSession1_BaitEmpty');
        data = load('ALK022_Exp22_GoodSession1_Bait[]');    
        data=data.data;
    elseif strcmp(AnimalID,'ALK025')
        data = load('ALK025_Exp22_GoodSession1_Bait[]');      
        data=data.data;
    elseif strcmp(AnimalID,'ALK027')
        data = load('ALK027_Exp22_GoodSession1_Bait[]');    
        data=data.data;
    end
    
elseif strcmp(ExpID,'23') % loading Exp 23 and other useful sessions for that exp
    if strcmp(AnimalID,'ALK013')
        %        data = load('ALK022_Exp22_GoodSession1_BaitEmpty');
        data = load('ALK013_Exp23_GoodSession[]_Bait[]');   
        data=data.data;
    elseif strcmp(AnimalID,'ALK017')
        data = load('ALK017_Exp23_GoodSession[]_Bait[]'); 
        data=data.data;
    elseif strcmp(AnimalID,'ALK028')
        data = load('ALK028_Exp23_GoodSession[]_Bait[]');  
        data=data.data;
    end
    
    
elseif strcmp(ExpID,'25') % loading Exp 25 and other useful sessions for that exp
    
    if strcmp(AnimalID,'ALK022')
        data = load('ALK022_Exp25_GoodSession[]_Bait0');
        data=data.data;
    elseif strcmp(AnimalID,'ALK025')
        data = load('ALK025_Exp25_GoodSession[]_Bait0');
        data=data.data;
    elseif strcmp(AnimalID,'ALK027')
        data = load('ALK027_Exp25_GoodSession[]_Bait0');
        data=data.data;
    end
    
elseif strcmp(ExpID,'6_DrD') % loading Exp 6_DrD and other useful sessions for that exp
    if strcmp(AnimalID,'ALK034')
        data = load('ALK034_Exp6DrD_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK038')
        data = load('ALK038_Exp6DrD_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK039')
        data = load('ALK039_Exp6DrD_GoodSession[]_Bait[]');
        data=data.data;
    end
elseif strcmp(ExpID,'12_DrD') % loading Exp 12_DrD and other useful sessions for that exp
    if strcmp(AnimalID,'ALK034')
        data = load('ALK034_Exp12DrD_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK038')
        data = load('ALK038_Exp12DrD_GoodSession[]_Bait[]');
        data=data.data;
    elseif strcmp(AnimalID,'ALK039')
        data = load('ALK039_Exp12DrD_GoodSession[]_Bait[]');
        data=data.data;
    end
end

function [all_AIC]= get_AIC( SampledParams, Top5Params, NLL )

    nP = length( SampledParams );
    all_AIC = nan(5,1);
    sampleNum = nan(nP, 1);
    for Ptype = 1:nP
        sampleNum( Ptype ) = length( SampledParams{Ptype} );
    end
    
    indxP = cell( nP, 1);
    for model = 1:5
        modelP = Top5Params( model, :);
        
        %first params is alpha
        alpha_diff = abs( SampledParams{1} - Top5Params( model, 1) );
        indxP{1} = find( alpha_diff == min(alpha_diff),1); 
        for Ptype = 2:nP
           indxP{Ptype} = find( modelP( Ptype ) == SampledParams{Ptype},1 ); 
        end
        nll_indx = sub2ind( sampleNum, indxP{:} );
        all_AIC( model ) = NLL(nll_indx) + 2*nP;
    end
end
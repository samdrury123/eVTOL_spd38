function Cut_Data = ts_structured_cut(Data,bid,Ist,Ien,Jst,Jen,Kst,Ken)
% Function to take a structured cut and return all secondary variables

    % Take out data of interest
    Data_Temp = Data{bid+1};
    
    % Check if ends have been specified
    ni = Data_Temp.attribute.ni;
    nj = Data_Temp.attribute.nj;
    nk = Data_Temp.attribute.nk;
    
    if strcmp('en',Ien) == 1
        Ien = ni;
    end
    if strcmp('en',Jen) == 1
        Jen = nj;
    end    
    if strcmp('en',Ken) == 1
        Ken = nk;
    end       
    
    if strcmp('en',Ist) == 1
        Ist = ni;
    end
    if strcmp('en',Jst) == 1
        Jst = nj;
    end    
    if strcmp('en',Kst) == 1
        Kst = nk;
    end     
    
    % Get names of all items in block, cut arrays by ien,ist etc and keep
    % the rest as is.
    names = fieldnames(Data_Temp);
    for i = 1:length(names)
        name = names{i};
        if ndims(Data_Temp.(name)) == 3
            Cut_Data.(name) = squeeze(Data_Temp.(name)(Ist:Ien,Jst:Jen,Kst:Ken));
        else
            Cut_Data.(name) = Data_Temp.(name);
        end
    end
       
end


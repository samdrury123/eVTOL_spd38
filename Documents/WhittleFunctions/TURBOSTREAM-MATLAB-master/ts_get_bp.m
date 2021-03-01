function [bp] = ts_get_bp(Data, bp_id)
% Read blockproperty data into a cell array of same format as coords and
% connectivity data

% Initialise bp
bp = cell(length(Data),1);

if length(Data) > 1
    for n = 1:length(Data)

        % Move bp data into cell array
        if isfield(Data{n},bp_id) == 1
            bp{n} = Data{n}.(bp_id);
        else
            bp{n} = [];
        end

    end
else
    bp = Data{1}.(bp_id);
end

end
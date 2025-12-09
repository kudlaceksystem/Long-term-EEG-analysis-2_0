function h = prepareFigures(stg, h, figDesc)
    for kfig = 1 : numel(figDesc.Name)
        fd = figDesc.(figDesc.Name(kfig));
        nm = fd.Name;
        if isfield(fd, "PositionCm")
            h.f.(nm) = figure("Units", "centimeters", "Position", fd.PositionCm);
        elseif isfield(fd, "subplotHeCm")
            if stg.sbNCol == 1
                wiCm = stg.figWidth1Cm;
            elseif stg.sbNCol == 2
                wiCm = (stg.figWidth1Cm + stg.figWidth2)/2;
            else
                wiCm = stg.figWidth2Cm;
            end
            heCm = min(25, fd.subplotHeCm*stg.sbNRow);
            positionCm = [20, 2, wiCm, heCm];
            h.f.(nm) = figure("Units", "centimeters", "Position", positionCm);
        end
        h.f.(nm).Units = "pixels";
        h.f.(nm).Name = nm;
        h.f.(nm).Color = [1 1 1];
    end
end

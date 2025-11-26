function getSubjAndSubjClr(subjToPlot, subjList, colorfulSubjects)
    arguments (Input)
    subjToPlot
    subjList
    colorfulSubjects
    end

    arguments (Output)
    end

    global stg
    if colorfulSubjects
        hf = figure;
        axes;
        subjClr = get(gca, 'ColorOrder');
        close(hf);
        delete(hf); % Dummy axes to get Matlab default color order
        subjClr = [max(1 - (1 - subjClr)*0.75, 0); max(1 - (1 - subjClr)*1.2, 0)]; % Each subject has different color
    else
        subjClr = ones(numel(subjList), 1)*stg.uniformSubjectColor; % All subjects have red
    end
    subjInd = ismember(subjList, subjToPlot);
    subjClr = subjClr(subjInd, :);
    stg.subjColor = subjClr;
    stg.subjNumber = find(subjInd);
end
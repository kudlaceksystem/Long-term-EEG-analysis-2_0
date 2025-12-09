function [spx, spy, spWi, spHe, numr, numc] = getSubplotXYWH(stg, h, d, margGlob, marg)
    % Note that the term subplot means a room for all possible plots of a given subject, not necessarily a single subplot
    % plotName .. name of the plot derived from the name of the calling functions
    % margGlob .. margins around the whole page, left, bottom, right, top
    % marg ...... margins around each subject's plots, left, bottom, right, top
    % spx ....... x-coordinate of lower left corner of the subjects' subplots in the figure, it is a vector, ksubj-th subject will use spx(mod(ksubj - 1, numc) + 1)
    % spy ....... y-coordinate of lower left corner of the subjects' subplots in the figure, it is a vector, ksubj-th subject will use spy(ceil(ksubj/numc))
    % spWi ...... width of subjects subplot
    % spHe ...... height of subjects subplot
    % numr ...... in how many rows the subjects will be plotted
    % numc ...... in how many columns the subjects will be plotted
    nm = d.Name;
    h.f.(nm).Units = "centimeters";
    figPos = h.f.(nm).Position;
    % % % % if contains(plotName, 'All') || contains(plotName, 'Pop')
    % % % %     numc = 1;
    % % % %     numr = 1;
    % % % % else
        numc = stg.sbNCol;
        numr = stg.sbNRow;
    % % % % end
    splWi = (figPos(3) - margGlob(1) - margGlob(3))/numc; % Subplot including labels width
    spWi = splWi - marg(1) - marg(3); % Subplot width
    splHe = (figPos(4) - margGlob(2) - margGlob(4))/numr; % Subplot including labels height
    spHe = splHe - marg(2) - marg(4); % Subplot height
    spx = (0 : numc-1)*splWi + margGlob(1) + marg(1); % Subplot - x-coordinate of its lower left corner
    spy = (numr-1 : -1 : 0)*splHe + margGlob(2) + marg(2); % Subplot - y-coordinate of its lower left corner
end
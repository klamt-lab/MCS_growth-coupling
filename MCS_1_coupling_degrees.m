% This script reproduces the results from Figure 4:
% It computes all MCS for the weakly and strongly growth-coupled and 
% substrate uptake production of ethanol with E. coli.
% pGCP: production potential at max growth rate
% wGCP: production at max growth rate
% sGCP: prodution at all positive growth rates
% SUCP: prodcution in all flux states
%
% MCS computation 
%
% % Content:
%   0) Start parallel pool (if available) to speed up computation
%   1) Setup model
%   2) Define Target and Desired regions for MCS computation
%   3) Run MCS computation
%   4) Validate MCS
%   5) Plot results
%
% Correspondence: cellnetanalyzer@mpi-magdeburg.mpg.de
% -Mar 2021
%

atpm = true; % atpm = true  - forces a minimum ATP maintenance rate (as in the original iML1515 model)
             % atpm = false - lifts the ATP maintenance constraint
maxSolutions = inf;
coupling = {'potential growth-coupling' 'weak growth-coupling' 'strong growth-coupling' 'substrate uptake coupling'};
maxCost = 3;
options.mcs_search_mode = 2;
verbose = 1;

%% 0) Starting CNA and Parallel pool (for faster FVA), defining computation settings
addpath(fullfile(fileparts(mfilename('fullpath')),'functions'));
if ~exist('cnan','var')
    startcna(1)
end

% if availabe, use parpool
if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
       (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && isempty(gcp('nocreate')) %#ok<DCRENAME>
    parpool();
    wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
    numworkers = getfield(gcp('nocreate'),'NumWorkers');
else
    numworkers = 0;
end

%% 1) Model setup
% load model from file
load(which('iML1515.mat'));
load(which('iML1515geneNames.mat'));
cnap = CNAcobra2cna(iML1515,0);
cnap = block_non_standard_products(cnap);
cnap.reacMin(ismember(cnap.reacID,{'EX_glc__D_e'})) = -10;
cnap.path = 'iMLcore';
load(which('core.mat'));
% Reduce model to core network
cnap = CNAdeleteReaction(cnap,find(~ismember(cellstr(cnap.reacID),core_reacs)));
cnap = CNAdeleteSpecies(cnap,find(~ismember(cellstr(cnap.specID),core_specs)),0);
% replace gene numbers with names
for i = 1:length(ecoliGeneNames)
    cnap.reacNotes = strrep(cnap.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
end
% Remove minimum ATP Maintenance if indicated. Substrate uptake can then be fixed to 10 as it remains the
% only heterologous bound of the model.
if ~atpm
	cnap.reacMin(~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'ATPM')))) = 0;
    cnap.reacMax(~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'EX_glc__D_e')))) = -10;
end
% replace bounds with inf
cnap.reacMin(cnap.reacMin==-1000) = -inf;
cnap.reacMax(cnap.reacMax== 1000) =  inf;
% lump several gene subunits
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'ATPS4rpp')),'geneProductAssociation','atpS*');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'NADH16pp')),'geneProductAssociation','nuo*');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'NADH17pp')),'geneProductAssociation','nuo*');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'NADH18pp')),'geneProductAssociation','nuo*');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'FRD2')),'geneProductAssociation','frd*');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'FRD3')),'geneProductAssociation','frd*');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'CYTBO3_4pp')),'geneProductAssociation','cyo*');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'THD2pp')),'geneProductAssociation','pnt*');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'PDH')),'geneProductAssociation','ace* and lpd');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'AKGDH')),'geneProductAssociation','sucAB and lpd');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'SUCOAS')),'geneProductAssociation','sucCD');
cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'SUCDi')),'geneProductAssociation','sdh*'); % sdhA,B,C,D
[~,~,genes,gpr_rules] = CNAgenerateGPRrules(cnap);

%% 2) Define MCS setup
% All genes are knockable apart from pseudo-gene that marks spontanous reactions
gkoCost = ones(length(genes),1);
gkoCost(ismember(genes,'spontanous')) = nan;
% Reactions are not knockable apart from O2 supply
rkoCost = nan(cnap.numr,1);
rkoCost(strcmp(cellstr(cnap.reacID),'EX_o2_e')) = 1;
[full_cnap, rmap] = CNAintegrateGPRrules(cnap);
product_rID   = 'EX_etoh_e';
substrate_rID = 'EX_glc__D_e';
biomass_rID   = 'BIOMASS_Ec_iML1515_core_75p37M';
for coupling_i = coupling
    clear('modules');
    % Desired fluxes for all coupling cases: growth >= 0.05/h
    modules{1}.sense = 'desired';
    modules{1}.type  = 'lin_constraints';
    [modules{1}.V(1,:),modules{1}.v(1,:)] = genV([{biomass_rID} {'>='} 0.05 ],cnap);

    disp(' ');
    disp('=============');
    disp(char(coupling_i));
    disp('=============');
    switch char(coupling_i)
        case 'potential growth-coupling'
            modules{1}.type  = 'bilev_w_constr';
            % Desired flux states with maximal growth and positive ethanol production
            [modules{1}.V(2,:),modules{1}.v(2,:)] = genV([{product_rID} {'>='} 0.01 ],cnap);
            [modules{1}.c(1,:),~]                 = genV([{biomass_rID} {'>='} 0    ],cnap); % maximize biomass
        case 'weak growth-coupling'
            modules{2}.sense = 'target';
            modules{2}.type  = 'bilev_w_constr';
            % Target flux states with maximal growth and no ethanol production
            [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{product_rID} {'<='} 0    ],cnap);
            [modules{2}.V(2,:),modules{2}.v(2,:)] = genV([{biomass_rID} {'>='} 0.05 ],cnap);
            [modules{2}.c(1,:),~]                 = genV([{biomass_rID} {'>='} 0 ],cnap); % maximize biomass
        case 'strong growth-coupling'
            modules{2}.sense = 'target';
            modules{2}.type = 'lin_constraints';
            [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{product_rID}   {'<='}  0     ],cnap);
            [modules{2}.V(2,:),modules{2}.v(2,:)] = genV([{biomass_rID}   {'>='}  0.001 ],cnap);
        case 'substrate uptake coupling'
            modules{2}.sense = 'target';
            modules{2}.type = 'lin_constraints';
            [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{product_rID}   {'<='}  0   ],cnap);
            modules_su = modules;
        otherwise
            error(['Coupling mode ''' char(coupling_i) ''' not found.']);
    end
    %% 3) Run MCS computation
    tic;
    [rmcs, mcs, full_cnap] = ...
        CNAgeneMCSEnumerator3(cnap,modules,...
            rkoCost,[],...
            maxSolutions,maxCost,...
            gkoCost,[],[],...
            options,verbose);
    disp(['Computation time: ' num2str(toc) ' s']);
    %% 4) Verify MCS
    disp('Verifying MCS');
    if ~isempty(mcs) % if mcs have been found
        if isfield(modules{1},'c') % pGCP
            [valid_T, valid_D] = verify_mcs(full_cnap,mcs,[],[],modules{1}.c*rmap,modules{1}.V*rmap,modules{1}.v);
        elseif isfield(modules{2},'c') % wGCP
            [valid_T, valid_D] = verify_mcs(full_cnap,mcs,modules{2}.V*rmap,modules{2}.v,modules{2}.c*rmap,modules{1}.V*rmap,modules{1}.v);
        else % sGCP & SUCP
            [valid_T, valid_D] = verify_mcs(full_cnap,mcs,modules{2}.V*rmap,modules{2}.v,[],modules{1}.V*rmap,modules{1}.v);
        end
        valid = ~valid_T & valid_D;
        if all(valid)
            disp('MCS VALID');
        else
            disp('MCS INVALID');
        end
    else
        valid = [];
        valid_T = [];
        valid_D = [];
    end
    if isempty(mcs)
        mcs = nan;
        text = 'no MCS found';
    else
        switch char(coupling_i)
        case 'potential growth-coupling'
            comp_time_pGCP = toc;
            mcs_pGCP   = mcs;
            disp([num2str(size(mcs_pGCP,2))   ' MCS for potential growth-coupling (after '      num2str(comp_time_pGCP) ' s)']);
        case 'weak growth-coupling'
            comp_time_wGCP = toc;
            mcs_wGCP   = mcs;
            disp([num2str(size(mcs_wGCP,2))   ' MCS for weak growth-coupling (after '      num2str(comp_time_wGCP) ' s)']);
        case 'strong growth-coupling'
            comp_time_sGCP = toc;
            mcs_sGCP = mcs;
            disp([num2str(size(mcs_sGCP,2)) ' MCS for strong growth-coupling (after '    num2str(comp_time_sGCP) ' s)']);
        case 'substrate uptake coupling'
            comp_time_SUCP = toc;
            mcs_SUCP  = mcs;
            disp([num2str(size(mcs_SUCP,2))  ' MCS for substrate uptake coupling (after ' num2str(comp_time_SUCP) ' s)']);
        end
    end
end

%% 5) Plot Results
disp('========')
disp('Results:')
reac_names = full_cnap.reacID;
reac_names(full_cnap.rType ~= 'r',:) = reac_names(full_cnap.rType ~= 'r',[4:end '   ']);
warning('off','MATLAB:hg:AutoSoftwareOpenGL');
f1 = figure;
p = plot_mcs_relationships(reac_names,mcs_wGCP,mcs_sGCP,mcs_SUCP);
f2 = figure;
p = plot_mcs_relationships(reac_names,mcs_pGCP,mcs_wGCP,mcs_sGCP);
warning('on','MATLAB:hg:AutoSoftwareOpenGL');
disp('Finished.');

%% Supplementary functions. 
%Generate Target or desired region from text.
function [V,v] = genV(constraints,cnap)
% Generate Vectors V and v so that V*r <= v
% input: some constraint seperated 3 three cells. e.g.:
%           r_1 + r_4 / r_3 - r_2    |    >=    |    a
    V = [];
    v = [];
    rMin = nan(cnap.numr,1);
    rMax = nan(cnap.numr,1);
    reacID = cellstr(cnap.reacID);

    for j = 1:size(constraints,1)
        % get right hand side
        a = constraints{j,3};
        if ischar(a)
            a = str2double(a);
            if isnan(a)
                error('error in ''t'' of target region definition');
            end
        end
        % get direction of inequality
        switch constraints{j,2}
            case '<='
                eqop = 1;
            case '>='
                eqop = -1;
            otherwise
                error('please define inequality either as ''<='' or ''>=''');
        end
        % split into numerator an divisor and get coefficients
        numDiv = strtrim(split(constraints{j,1},'/'));
        % get variables and coefficients for vector
        [num, cnum] = findReacAndCoeff(numDiv(1),reacID);
            % fractional constraint ispreprocessed
        if length(numDiv) == 2
            [div, cdiv] = findReacAndCoeff(numDiv(2),reacID);
            [rMin(div), rMax(div)] = CNAfluxVariability(cnap,[],[],-1,div,[],[],0);
            % check if all reactions take identical signs (+ or -)
            % this is needed to do the equation rearrangement. If the signs
            % are ambigous, a lot of case differentiations would be needed,
            % what is not done here.
            if any(rMin(div).*rMax(div) < 0)
                error(['reactions that are part of the divisor inside the '... 
                        'target constraints must not span a positive AND negative range']);
            end
            rDir = sign(rMin(div)+rMax(div));
            if any(rDir.*cdiv' > 0) &&  any(rDir.*cdiv' < 0)
                error(['reactions that are part of the divisor inside the '...
                       'target constraints must all have the same direction '...
                       '(positive or negative). Please check if coefficients and '...
                       'reaction ranges lead to all positive or all negative variables.']);
            end
            if sign(sum(rDir.*cdiv')) == -1 % if the divisor is all negative
                % change direction of inequality
                eqop = -eqop;
            end
            cdiv = -a*cdiv; % transformation to following form:
            a = 0;
            % (cnum num) - a*cdiv1*div1 - a*cdiv2*div2 <= 0
            num  = [num,   div];
            cnum = [cnum, cdiv];
        end
        % constraint is generated
        if eqop == -1  % if inequality is >= the whole term is multiplied with -1
            cnum = -cnum;
            a    = -a;
        end
        v(j,1)   = a;
        V(j,:) = full(sparse(num,1,cnum,length(reacID),1));
    end
end

function [ridx,coeff] = findReacAndCoeff(eq,reacID)
    coeff = [];
    r = cell.empty(1,0);
    ridx = [];
    for strPart = strsplit(char(eq))
        str = char(regexprep(strPart,'^(\s|-|\.|\()*|(\s|-|\.|\))*$','')); % remove leading and tailing special characters
        if ~isempty(str)
            v = regexp(reacID, ['^' str '$'], 'match');
            if any(~cellfun(@isempty,v))
                r(end+1) = {str};
                ridx = [ridx find(~cellfun(@isempty,v))'];
            end
        end
    end
    for k = 1:length(r)
        c = regexp(eq, ['(\s|\d|-|\.)*?(?=' r{k} '(\s|$))'], 'match');
        c = regexprep(char(c{:}),'\s','');
        switch c
            case ''
                coeff(k) = 1;
            case '+'
                coeff(k) = 1;
            case '-'
                coeff(k) = -1;
            otherwise
                coeff(k) = str2double(c);
        end
    end
end

function [valid_T, valid_D] = verify_mcs(cnap,mcs,T,t,c,D,d)
    if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
           (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && ~isempty(gcp('nocreate')) %#ok<DCRENAME>
        numworkers = getfield(gcp('nocreate'),'NumWorkers');
    else
        numworkers = 0;
    end
    mcs = mcs<0;
    valid_T = nan(size(mcs,2),1);
    valid_D = nan(size(mcs,2),1);
    parfor (i = 1:size(mcs,2),numworkers)
        cnap_valid = cnap;
        cnap_valid.reacMin(mcs(:,i)) = 0;
        cnap_valid.reacMax(mcs(:,i)) = 0;
        cnap_valid_opt = cnap_valid;
        if ~isempty(c)
            cnap_valid.objFunc = c';
            fv = CNAoptimizeFlux(cnap_valid,[], [], 2, -1);
            cnap_valid_opt.reacMin(logical(c)) = fv(logical(c));
            cnap_valid_opt.reacMax(logical(c)) = fv(logical(c));
        end
        if isempty(T)
            valid_T(i) = 0;
            valid_D(i) = testRegionFeas(cnap_valid_opt,D,d,2);
        else
            valid_T(i) = testRegionFeas(cnap_valid_opt,T,t,2);
            valid_D(i) = testRegionFeas(cnap_valid,D,d,2);
        end
    end
end

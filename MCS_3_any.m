% This script reproduces the results from Table 3:
% It computes a random MCS for the weakly and strongly growth-coupled
% and substrate uptake production of 10 different products with E. coli.
% pGCP: production potential at max growth rate
% wGCP: production at max growth rate
% sGCP: prodution at all positive growth rates
% SUCP: prodcution in all flux states
% Computations are repeated 12 times and have a time limit of 2 hours.
%
% MCS computation 
%
% % Content:
%   0) Start parallel pool (if available) to speed up computation
%   1) Setup model (add pathways if necessary)
%   2) Define Target and Desired regions for MCS computation
%   3) Run MCS computation
%   4) Validate MCS
%   5) Plot results
%
% Correspondence: cellnetanalyzer@mpi-magdeburg.mpg.de
% -Mar 2021
%

%% User settings

% select coupling type:
% potential growth-coupling
% weak growth-coupling
% strong growth-coupling
% substrate uptake coupling
coupling = 'potential growth-coupling';

% select product
%  1: ethanol
%  2: lysine
%  3: glutamate
%  4: isobutanol
%  5: 1,4-BDO
%  6: 2,3-BDO
%  7: itaconic acid
%  8: isoprene
%  9: butane
% 10: methacrylic acid
% 11: resveratrol
% 12: bisabolene
product = 1;

% select number and time limit for computations
num_iter = 12; % 12 computations
options.milp_time_limit = 7200; % 2h

% choose whether a new instance of MATLAB should be started for each MCS computation.
% This measure ensures memory cleanup between the MCS runs. (usually not needed with Windows)
solve_in_new_process = 1;

options.mcs_search_mode = 1; % find a random MCS
maxSolutions = 1;
maxCost = 60;
verbose = 1;

%% 0) Starting CNA and Parallel pool (for faster FVA), defining computation settings
addpath(fullfile(fileparts(mfilename('fullpath')),'functions'));
if ~exist('cnan','var')
    try
        startcna(1)
    catch
        error('CellNetAnalyzer could not be started. Make sure that it is installed and added to the MATLAB path.');
    end
end
if ismember(3,find(LP_solver_availability(true)))
    options.milp_solver = 'cplex';
elseif ismember(4,find(LP_solver_availability(true)))
    options.milp_solver = 'gurobi';
else
    error('This script requires either CPLEX or Gurobi. Make sure that at least one of these two is set up with CNA.');
end

% if availabe, use parpool
if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
       (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && isempty(gcp('nocreate')) && ~solve_in_new_process %#ok<DCRENAME>
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

% % Uncomment this for computation in core network
% load(which('core.mat'));
% % Reduce model to core network
% cnap = CNAdeleteReaction(cnap,find(~ismember(cellstr(cnap.reacID),core_reacs)));
% cnap = CNAdeleteSpecies(cnap,find(~ismember(cellstr(cnap.specID),core_specs)),0);

% replace gene numbers with names
for i = 1:length(ecoliGeneNames)
    cnap.reacNotes = strrep(cnap.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
end
% Load heterologous pathways if necessary
[product_rID,species,reactions] = load_pathway(product);
for spec = species
    cnap = CNAaddSpeciesMFN(cnap,spec.spec_id,0,spec.spec_name);
    cnap = CNAsetGenericSpeciesData(cnap,cnap.nums,'fbc_chemicalFormula',char(spec.fbc_chemicalFormula),'fbc_charge',double(spec.fbc_charge));
end
for reac = reactions
    cnap = CNAaddReactionMFN(cnap, reac.reac_id, reac.equation, reac.lb, reac.ub,0,nan,nan,'',0,0,0,0);
    cnap = CNAsetGenericReactionData(cnap,cnap.numr,'geneProductAssociation',char(reac.fbc_geneProductAssociation));
end

% Checking mass and and charge balances after pathway additions
check_mass_balance(cnap);

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
substrate_rID = 'EX_glc__D_e';
biomass_rID   = 'BIOMASS_Ec_iML1515_core_75p37M';

% Desired fluxes for all coupling cases: growth >= 0.05/h
modules{1}.sense = 'desired';
modules{1}.type  = 'lin_constraints';
[modules{1}.V(1,:),modules{1}.v(1,:)] = genV([{biomass_rID} {'>='} 0.05 ],cnap);

%% compute thresholds:
% wGCP: At 20% of maximal growth rate, maximize production rate. 
%       20% of this rate should be ensured at maximal growth after 
%       the interventions.
% sGCP: At 20% of maximal growth rate, maximize production rate. 
%       20% of this rate devided by the 20% of the maximal growth rate
%       is the ratio between production and growth that should be attained
%       in all flux states.
% SUCP: At 20% of maximal growth rate, maximize production rate. 
%       20% of this rate, devided by the maximal substrate uptake rate
%       should be attained in all flux states.

disp('Compute production (yield) thresholds.');
% 1. compute max growth
cnap.objFunc(:) = 0;
cnap.objFunc(ismember(cnap.reacID,{biomass_rID})) = -1;
fv = CNAoptimizeFlux(cnap,[],[],2,-1);
r_bm_max20 = 0.2*fv(ismember(cnap.reacID,{biomass_rID}));
% 2. compute max production at 20% growth
cnap.objFunc(:) = 0;
cnap.objFunc(ismember(cnap.reacID,{product_rID})) = -1;
fv_fix = nan(cnap.numr,1);
fv_fix(ismember(cnap.reacID,{biomass_rID})) = r_bm_max20;
fv = CNAoptimizeFlux(cnap,fv_fix,[],2,-1);
r_p_20 = 0.2*fv(ismember(cnap.reacID,{product_rID}));
% 3. Thresholds for sGCP and SUCP
Y_PBM = r_p_20/r_bm_max20;
Y_PS  = r_p_20/-fv(ismember(cnap.reacID,{substrate_rID}));

clear('modules');
% Desired fluxes for all coupling cases: growth >= 0.05/h
modules{1}.sense = 'desired';
modules{1}.type  = 'lin_constraints';
[modules{1}.V(1,:),modules{1}.v(1,:)] = genV([{biomass_rID} {'>='} 0.05 ],cnap);

disp(' ');
disp('=============');
disp(coupling);
disp('=============');
switch coupling
    case 'potential growth-coupling'
        modules{1}.type  = 'bilev_w_constr';
        % Desired flux states with maximal growth and positive ethanol production
        [modules{1}.V(2,:),modules{1}.v(2,:)] = genV([{product_rID} {'>='} r_p_20 ],cnap);
        [modules{1}.c(1,:),~]                 = genV([{biomass_rID} {'>='} 0    ],cnap); % maximize biomass
    case 'weak growth-coupling'
        modules{2}.sense = 'target';
        modules{2}.type  = 'bilev_w_constr';
        [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{product_rID} {'<='} r_p_20 ],cnap);
        [modules{2}.V(2,:),modules{2}.v(2,:)] = genV([{biomass_rID} {'>='} 0.05   ],cnap);
        [modules{2}.c(1,:),~]                 = genV([{biomass_rID} {'>='} 0      ],cnap); % maximize biomass
    case 'strong growth-coupling'
        modules{2}.sense = 'target';
        modules{2}.type = 'lin_constraints';
        [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{[product_rID ' / ' biomass_rID]}   {'<='}  Y_PBM ],cnap);
        [modules{2}.V(2,:),modules{2}.v(2,:)] = genV([{biomass_rID}   {'>='}  0.01 ],cnap);
    case 'substrate uptake coupling'
        modules{2}.sense = 'target';
        modules{2}.type = 'lin_constraints';
        [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{[product_rID ' / -' substrate_rID]}    {'<='}  Y_PS],cnap);
    otherwise
        error(['Coupling mode ''' coupling ''' not found.']);
end

comptime = nan(num_iter,1);
mcs = nan(full_cnap.numr,num_iter);
for i = 1:num_iter
    if solve_in_new_process    
       displ(['Running MCS computation in new process: ' num2str(i) '/' num2str(num_iter)],verbose);
       [mcs_i, comptime(i)] = MCS_enum_thread(1,cnap,modules,...
                                        rkoCost,...
                                        maxSolutions,maxCost,...
                                        gkoCost,...
                                        options,verbose);
    else
        displ(['Running MCS computation: ' num2str(i) '/' num2str(num_iter)],verbose);
        tic;
        [rmcs, mcs_i, full_cnap] = ...
            CNAgeneMCSEnumerator3(cnap,modules,...
                rkoCost,[],...
                maxSolutions,maxCost,...
                gkoCost,[],[],...
                options,verbose);
        comptime(i) = toc;
    end
    if ~isempty(mcs_i)
        mcs(:,i) = mcs_i(:,find(sum(abs(mcs_i),1) == min(sum(abs(mcs_i),1)),1));
    end
end
mcs_size = sum(abs(mcs),1);
succ = ~isnan(mcs_size);
size_smallest = min(mcs_size(succ));
size_average = mean(mcs_size(succ));
comptime_succ = mean(comptime(mcs_size == min(mcs_size)));
comptime_shortest = min(comptime);
comptime_mean = mean(comptime);

comptime_succ = num2str(comptime_succ);

%% 4) Verify MCS
disp('Verifying MCS');
if ~isempty(mcs(:,succ)) % if mcs have been found
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
if isempty(succ)
    mcs = nan;
    text = 'no MCS found';
else
    text = ['av. mcs size: ' num2str(mean(mcs_size(succ))) ' | valid/found/num_runs: ' num2str(sum(valid)) '/' num2str(sum(succ)) '/' num2str(num_iter) ...
            ' | ' strjoin(cellstr(num2str(mcs_size')),',') ' | times: ' strjoin(cellstr(num2str(comptime)),',')];
end

disp('========')
disp(['Results - any MCS - ' coupling ': ' product_rID])
disp(['Shortest computation time: ' num2str(comptime_shortest) ' s']);
disp(['Mean computation time: ' num2str(comptime_mean) ' s']);
disp(['Successful computations: ' num2str(sum(succ))]);
disp(['Smallest MCS size: ' num2str(size_smallest)]);
disp(['Average MCS size: ' num2str(size_average)]);
disp(['Mean computation time of shortest MCS: ' num2str(comptime_succ) ' s']);
disp(text);
disp('Finished.');


%% Supplementary function. Generate Target or desired region from text.

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
            fv = CNAoptimizeFlux(cnap_valid,[], [], 0, -1);
            cnap_valid_opt.reacMin(logical(c)) = fv(logical(c));
            cnap_valid_opt.reacMax(logical(c)) = fv(logical(c));
        end
        if isempty(T)
            valid_T(i) = 0;
            valid_D(i) = testRegionFeas(cnap_valid_opt,D,d);
        else
            valid_T(i) = testRegionFeas(cnap_valid_opt,T,t);
            valid_D(i) = testRegionFeas(cnap_valid,D,d);
        end
    end
end

function [mcs, comptime] = MCS_enum_thread(gene_mcs,cnap,modules,...
                                            rkoCost,...
                                            maxSolutions,maxCost,...
                                            gkoCost,...
                                            options,verbose)
        % start MCS Enumeration in separate MATLAB instance
        if ~isempty(getenv('SLURM_TEMP'))
            tempdir = getenv('SLURM_TEMP');
        else
            tempdir = eval('tempdir');
        end
        rng('shuffle');
        wdir = [tempdir num2str(randi([0,1e6-1])) filesep];
        mkdir(wdir);
        filename = [wdir 'ws_' num2str(randi([0,1e6-1])) '.mat'];
        save(filename,'gene_mcs','cnap','modules','rkoCost','maxSolutions','maxCost','gkoCost','options','verbose');
        wd = pwd;
        cd(evalin('base','cnan.cnapath'));
        system(['matlab -nosplash -nodesktop -r "addpath(''' genpath(fileparts(mfilename('fullpath'))) ''');' ... % add project path
               'addpath(''' evalin('base','cnan.cnapath') ''');startcna(1);' ... % add CNA path and start CNA
               'compute_mcs_ext(''' wdir ''',''' filename ''');exit()"']); % compute
        cd(wd);
        load(filename,'mcs','comptime','status');
        if status ~= 0
            comptime = nan;
        end
        rmdir(wdir,'s');
end
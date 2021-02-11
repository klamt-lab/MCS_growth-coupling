% This script reproduces the results from Table 3:
% It computes the smallest MCS for the weakly and strongly growth-coupled
% and substrate uptake production of 10 different products with E. coli.
% sGCP: prodution at all positive growth rates
% ACP: prodution at all flux states with ATP maintenance
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
% weak growth-coupling
% strong growth-coupling
% substrate uptake coupling
coupling = 'weak growth-coupling';

% select product
%  1: ethanol
%  2: isobutanol
%  3: 1,4-BDO
%  4: 2,3-BDO
%  5: itaconic acid
%  6: tryptophane
%  7: lysine
%  8: glutamate
%  9: isoprene
% 10: octyl acetate
% 11: butane
% 12: methacrylic acid
% 13: resveratrol
% 14: bisabolene
% 15: 4-hydroxycoumarin
product = 1;

% choose whether a new instance of MATLAB should be started for each MCS computation.
% This measure ensures memory cleanup between the MCS runs. (usually not needed with Windows)
solve_in_separ_thread = 1;

options.mcs_search_mode = 3;
maxSolutions = 1;
verbose = 1;

%% 0) Starting CNA and Parallel pool (for faster FVA), defining computation settings
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
% Remove minimum ATP Maintenance. Substrate uptake can then be fixed to 10 as it remains the
% only heterologous bound of the model.
cnap.reacMin(~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'ATPM')))) = 0;
cnap.reacMax(~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'EX_glc__D_e')))) = -10;

[specs, reacs] = load_pathway(product);

for spec = specs
    cnap = CNAaddSpeciesMFN(cnap,spec);
    cnap = CNAsetGenericSpeciesData(cnap,cnap.nums,'fbc_chemicalFormula',char(spec.fbc_chemicalFormula),'fbc_charge',double(spec.fbc_charge));
end
for reac = reacs
    cnap = CNAaddReactionMFN(cnap,reac);
    cnap = CNAsetGenericReactionData(cnap,cnap.numr,'geneProductAssociation',char(reac.fbc_geneProductAssociation));
end

productID = num2str(productID,'%02i');
filesinPath = dir('StrainBooster/my_mcs_setups');
prod = regexp({filesinPath.name}, ['^P',productID,'_.*.xls.*'], 'match'); % PXX_ABC
prod = prod{~cellfun(@isempty,prod)};
[~,prod_name] = fileparts(prod{:});
[~,model_name] = fileparts(cnap.path);

disp(['Loading reactions and species from file: ' strjoin(prod,', ')]);
try %% Add new reactions to model from product-xls (if any were defined)
    cnap = CNAaddSpecsAndReacsFromFile(cnap,prod{:});
catch err
    warning('off','backtrace')
    warning( getReport( err, 'extended', 'hyperlinks', 'on' ) );
    warning([char(prod) ': no reactions were added to model']);
    warning('on','backtrace')
end
%% 3. Checking mass and and charge balances
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
[~,~,genes,gpr_rules] = CNAgenerateGPRrules(cnap);


%% 2) Define MCS setup
% reaction indices used in Target and Desired regions
productID = num2str(productID,'%02i');
filesinPath = dir('StrainBooster/my_mcs_setups');
prod = regexp({filesinPath.name}, ['^P',productID,'_.*.xls.*'], 'match'); % PXX_ABC
prod = prod{~cellfun(@isempty,prod)};
[~,prod_name] = fileparts(prod{:});
[~,model_name] = fileparts(cnap.path);

disp(['Loading reactions and species from file: ' strjoin(prod,', ')]);
try %% Add new reactions to model from product-xls (if any were defined)
    cnap = CNAaddSpecsAndReacsFromFile(cnap,prod{:});
catch err
    warning('off','backtrace')
    warning( getReport( err, 'extended', 'hyperlinks', 'on' ) );
    warning([char(prod) ': no reactions were added to model']);
    warning('on','backtrace')
end
%% 3. Checking mass and and charge balances
check_mass_balance(cnap);

%% 4. Get GPR-associations and set knockability
if gene_mcs
    % For gene MCS
    % lump gene subunits
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'ATPS4rpp')),'geneProductAssociation','atpS*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'NADH16pp')),'geneProductAssociation','nuo*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'NADH17pp')),'geneProductAssociation','nuo*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'NADH18pp')),'geneProductAssociation','nuo*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'FRD2')),'geneProductAssociation','frd*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'FRD3')),'geneProductAssociation','frd*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'CYTBO3_4pp')),'geneProductAssociation','cyo*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'THD2pp')),'geneProductAssociation','pnt*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'PDH')),'geneProductAssociation','ace* and lpd');
    [~,~,genes,gpr_rules] = CNAgenerateGPRrules(cnap);
    gkoCost = ones(length(genes),1);
    % pseudo-gene that marks spontanous reactions is not knockable
%     gkoCost(ismember(genes,{'spontanous' 'sufA' 'sufB' 'sufC' 'sufD' 'sufE' 'sufS'})) = nan; % maybe also: 'iscA' 'iscS' 'iscU'
    gkoCost(ismember(genes,'spontanous')) = nan;
    rkoCost = nan(cnap.numr,1);
    rkoCost(strcmp(cellstr(cnap.reacID),'EX_o2_e')) = 1;
    [full_cnap, rmap] = CNAintegrateGPRrules(cnap);
else
    % For reaction MCS
    % knockables: All reactions with gene rules + O2 exchange as a potential knockout
    gpr = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
    koCost = double(cellfun(@(x) ~isempty(x),gpr));
    notknock_gene = ~cellfun(@isempty,regexp(gpr,'(spontanous|phoE|ompF|ompN|ompC)','match'));
%     notknock_gene = ~cellfun(@isempty,regexp(gpr,'(spontanous|phoE|ompF|ompN|ompC|sufA|sufB|sufC|sufD|sufE|sufS)','match')); % maybe also iscA|iscS|iscU
    koCost(koCost==0 | notknock_gene) = nan;
%     koCost(koCost==0) = nan;
    koCost(strcmp(cellstr(cnap.reacID),'EX_o2_e')) = 1;
    full_cnap = cnap;
    rmap = eye(cnap.numr);
end
%% 5. Construct MCS problem from excel file
table_prod = loadSpecReacXLStoStrArray(char(prod));
product_rID   = strsplit(char(readCells(table_prod,'product',1)));
product_rID   = char(product_rID(end));
substrate_rID = 'EX_glc__D_e';
biomass_rID   = 'BIOMASS_Ec_iML1515_core_75p37M';

modules{1}.sense = 'desired';
modules{1}.type  = 'lin_constraints';
[modules{1}.V(1,:),modules{1}.v(1,:)] = genV([{biomass_rID} {'>='} 0.05 ],cnap);

%% compute thresholds:
% wGCP: At 20% of maximal growth rate, maximize production rate. 
%       30% of this rate should be ensured at maximal growth after 
%       the interventions.
% sGCP: At 20% of maximal growth rate, maximize production rate. 
%       30% of this rate devided by the 20% of the maximal growth rate
%       is the ratio between production and growth that should be attained
%       in all flux states.
% SUCP: At 20% of maximal growth rate, maximize production rate. 
%       30% of this rate, devided by the maximal substrate uptake rate
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

for coupling_i = coupling
    modules{2} = struct;
    clear('T','t');
    disp(' ');
    disp('=============');
    disp(char(coupling_i));
    disp('=============');
    switch char(coupling_i)
        case 'weak growth-coupling'
            modules{2}.sense = 'target';
            modules{2}.type  = 'bilev_w_constr';
            [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{product_rID} {'<='} r_p_20 ],cnap);
            [modules{2}.V(2,:),modules{2}.v(2,:)] = genV([{biomass_rID} {'>='} 0.05   ],cnap);
            [modules{2}.c(1,:),~]                 = genV([{biomass_rID} {'>='} 0      ],cnap); % maximize biomass
            modules_wg = modules;
        case 'strong growth-coupling'
            [T_sg(1,:),t_sg(1,:)] = genV([{[product_rID ' / ' biomass_rID]}   {'<='}  Y_PBM ],cnap);
            [T_sg(2,:),t_sg(2,:)] = genV([{biomass_rID}   {'>='}  0.01 ],cnap);
            modules{2}.sense = 'target';
            modules{2}.type = 'lin_constraints';
            [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{[product_rID ' / ' biomass_rID]}   {'<='}  Y_PBM ],cnap);
            [modules{2}.V(2,:),modules{2}.v(2,:)] = genV([{biomass_rID}   {'>='}  0.01 ],cnap);
            modules_sg = modules;
        case 'substrate uptake coupling'
            [T_su(1,:),t_su(1,:)] = genV([{[product_rID ' / -' substrate_rID]}    {'<='}  Y_PS], cnap);
%             [T_su(2,:),t_su(2,:)] = genV([{substrate_rID} {'<='}  -1    ],cnap);
            modules{2}.sense = 'target';
            modules{2}.type = 'lin_constraints';
            [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{[product_rID ' / -' substrate_rID]}    {'<='}  Y_PS],cnap);
            modules_su = modules;
        otherwise
            error(['Coupling mode ''' char(coupling_i) ''' not found.']);
    end
    % 12 computations with 2h each
    options.milp_time_limit = 7200; % return best solution after 24h
    num_iter = 12;
    comptime = nan(num_iter,1);
    mcs = nan(full_cnap.numr,num_iter);
    for i = 1:num_iter
        if solve_in_separ_thread    
            displ(['Running MCS computation in separate thread: ' num2str(i) '/' num2str(num_iter)],verbose);
            if gene_mcs  %#ok<*UNRCH>
               [mcs_i, comptime(i)] = MCS_enum_thread(gene_mcs,cnap,modules,...
                                                rkoCost,...
                                                maxSolutions,maxCost,...
                                                gkoCost,...
                                                options,verbose);
            else
                [mcs_i, comptime(i)] = MCS_enum_thread(gene_mcs,cnap,  modules,...
                                                koCost,...
                                                maxSolutions,maxCost,...
                                                [],...
                                                options,verbose);
            end
        else
            displ(['Running MCS computation: ' num2str(i) '/' num2str(num_iter)],verbose);
            tic;
            if gene_mcs  %#ok<*UNRCH>
                [rmcs, mcs_i, full_cnap] = ...
                    CNAgeneMCSEnumerator3(cnap,modules,...
                        rkoCost,[],...
                        maxSolutions,maxCost,...
                        gkoCost,[],[],...
                        options,verbose);
            else
                 [mcs_i, status] = CNAMCSEnumerator3(cnap,modules,...
                        koCost,[],...
                        maxSolutions,maxCost,...
                        options,verbose);
            end
            comptime(i) = toc;
        end
        if ~isempty(mcs_i)
            mcs(:,i) = mcs_i(:,find(sum(abs(mcs_i),1) == min(sum(abs(mcs_i),1)),1));
        end
    end
    mcs_size = sum(abs(mcs),1);
    succ = ~isnan(mcs_size);
    size_succ = min(mcs_size(succ));
    comptime_succ = mean(comptime(mcs_size == min(mcs_size)));
    comptime_mean = mean(comptime);

    disp(['Mean computation time: ' num2str(comptime_mean) ' s']);
    disp(['Successful computations: ' num2str(sum(succ))]);
    disp(['Smallest MCS size: ' num2str(size_succ)]);
    disp(['Mean computation time of shortest MCS: ' num2str(comptime_succ) ' s']);
    
    comptime_succ = num2str(comptime_succ);

    %     if toc > options.milp_time_limit
    %         comptime = ['(timeout) ' comptime];
    %     end

    disp('========')
    disp('Results:')
    datetime.setDefaultFormats('default','yyyy-MM-dd hh:mm:ss');
    if ~isempty(mcs(:,succ)) % if mcs have been found
        if isfield(modules{2},'c')
            [valid_T, valid_D] = verify_mcs(full_cnap,mcs(:,succ),modules{2}.V*rmap,modules{2}.v,modules{2}.c*rmap,modules{1}.V*rmap,modules{1}.v);
        else
            [valid_T, valid_D] = verify_mcs(full_cnap,mcs(:,succ),modules{2}.V*rmap,modules{2}.v,[],modules{1}.V*rmap,modules{1}.v);
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
        if all(valid)
            text = ['av. mcs size: ' num2str(mean(mcs_size(succ))) ' | valid/found/num_runs: ' num2str(sum(valid)) '/' num2str(sum(succ)) '/' num2str(num_iter) ...
                    ' | ' strjoin(cellstr(num2str(mcs_size')),',') ' | times: ' strjoin(cellstr(num2str(comptime)),',')];
        end
    end
    disp(['Smallest MCS size for weak growth-coupling: ' num2str(size_succ) ' (after ' num2str(comptime_succ) ' s)']);

    switch char(coupling_i)
        case 'weak growth-coupling'
            comp_time_wg = comptime_succ;
            min_MCS_wg   = num2str(size_succ);
            weakGC_text  = text;
        case 'strong growth-coupling'
            comp_time_sg = comptime_succ;
            min_MCS_sg   = num2str(size_succ);
            strongGC_text = text;
        case 'substrate uptake coupling'
            comp_time_su = comptime_succ;
            min_MCS_su   = num2str(size_succ);
            substUp_text = text;
    end
end

if gene_mcs
    reac_names = full_cnap.reacID;
    reac_names(full_cnap.rType ~= 'r',:) = reac_names(full_cnap.rType ~= 'r',[4:end '   ']);
else
    reac_names = cnap.reacID;
end

if exist('setups','var')
    xls_filename = which('2_smallest_MCS.xls');
    disp('Reading file..');
    tab = readXls(xls_filename);
    [row,col,sheet]   = ind2sub(size(tab),find(strcmp(strtrim(tab),num2str(setups(str2double(getenv('SLURM_ARRAY_TASK_ID')))))));

    if isempty(col)
        disp(['array ID ' num2str(setups(str2double(getenv('SLURM_ARRAY_TASK_ID')))) ' not listed.']);
    else
        if ismember('weak growth-coupling',coupling)
            tab{row+1,col} = min_MCS_wg;
            tab{row+2,col} = weakGC_text;
            tab{row+3,col} = [comp_time_wg ' (' char(datetime) ')'];
        end
        if ismember('strong growth-coupling',coupling)
            tab{row+1,col} = min_MCS_sg;
            tab{row+2,col} = strongGC_text;
            tab{row+3,col} = [comp_time_sg ' (' char(datetime) ')'];
        end
        if ismember('substrate uptake coupling',coupling)
            tab{row+1,col} = min_MCS_su;
            tab{row+2,col} = substUp_text;
            tab{row+3,col} = [comp_time_su ' (' char(datetime) ')'];
        end
        disp('Writing file..');
        writecell(tab,xls_filename);
    end
end
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

% read certain cells
function C = readCells(table,keyword,cols) % returns all cells below the given keyword until the first empty line
    [row,col,sheet]   = ind2sub(size(table),find(strcmp(strtrim(table),keyword)));
    for i = 1:length(sheet)
        lastrow = row(i)+find(strcmp(strtrim(table(row(i):end,col(i),sheet)),''),1,'first')-2;
        if isempty(lastrow)
            lastrow = size(table,1);
        end
        C{i} = table((row(i)+1):lastrow,col(i):(col(i)+cols-1),sheet(i));
    end
    if ~exist('C','var')
        error(['ERROR: Keyword ''' keyword ''' was not found in any sheet.'])
    elseif length(C) == 1
        C = C{:};
    else
        disp(['WARNING: Keyword ''' keyword ''' was found multiple times in the sheets.']);
    end
end

function tab = readXls(filename)
    sheets = sheetnames(filename);
    for i=1:length(sheets)
        workSheet = readcell(filename,'sheet',i);
        if ~isempty(workSheet)
            tab(1:size(workSheet,1),1:size(workSheet,2),i) = workSheet;
        end
    end
    tab = strTable(tab);
    
    function rawmat = strTable(rawmat)
        [rows, cols, shts] = size(rawmat);
        for row = 1:rows
            for col = 1:cols
                for sh = 1:shts
                    if isnumeric(rawmat{row,col,sh})
                        if isnan(rawmat{row,col,sh}) || isempty(rawmat{row,col,sh})
                            rawmat(row,col,sh) = {''};
                        else
                            rawmat{row,col,sh} = num2str(rawmat{row,col,sh});
                        end
                    elseif isa(rawmat{row,col,sh},'missing')
                        rawmat(row,col,sh) = {''};
                    end
                end
            end
        end
    end
end

function [valid_T, valid_D] = verify_mcs_strong_substrate(cnap,mcs,T,t,D,d)
    mcs = mcs<0;
    valid_T = nan(size(mcs,2),1);
    valid_D = nan(size(mcs,2),1);
    parfor i = 1:size(mcs,2)
        cnap_valid = cnap;
        cnap_valid.reacMin(mcs(:,i)) = 0;
        cnap_valid.reacMax(mcs(:,i)) = 0;
		valid_T(i) = testRegionFeas(cnap_valid,T,t,2);
        valid_D(i) = testRegionFeas(cnap_valid,D,d,2);
    end
end

function [valid_T, valid_D] = verify_mcs_weak(cnap,mcs,T,t,c,D,d)
    mcs = mcs<0;
    valid_T = nan(size(mcs,2),1);
    valid_D = nan(size(mcs,2),1);
    parfor i = 1:size(mcs,2)
        cnap_valid = cnap;
        cnap_valid.reacMin(mcs(:,i)) = 0;
        cnap_valid.reacMax(mcs(:,i)) = 0;
        cnap_valid_opt = cnap_valid;
        cnap_valid.objFunc = c';
        fv = CNAoptimizeFlux(cnap_valid,[], [], 2, -1);
        cnap_valid_opt.reacMin(logical(c)) = fv(logical(c));
        cnap_valid_opt.reacMax(logical(c)) = fv(logical(c));
		valid_T(i) = testRegionFeas(cnap_valid_opt,T,t,2);
        valid_D(i) = testRegionFeas(cnap_valid,D,d,2);
    end
end

function [valid_T, valid_D] = verify_mcs(cnap,mcs,T,t,c,D,d)
    mcs = mcs<0;
    valid_T = nan(size(mcs,2),1);
    valid_D = nan(size(mcs,2),1);
    parfor i = 1:size(mcs,2)
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
		valid_T(i) = testRegionFeas(cnap_valid_opt,T,t,2);
        valid_D(i) = testRegionFeas(cnap_valid,D,d,2);
    end
end

function [specs, reacs] = load_pathway(product)
    switch product
        case  1 % ethanol
        case  2 % isobutanol
        case  3 % 1,4-BDO
        case  4 % 2,3-BDO
        case  5 % itaconic acid
        case  6 % tryptophane
        case  7 % lysine
        case  8 % glutamate
        case  9 % isoprene
        case 10 % octyl acetate
        case 11 % butane
        case 12 % methacrylic acid
        case 13 % resveratrol
        case 14 % bisabolene
        case 15 % 4-hydroxycoumarin
    end
end

function ko_ki_text = mcs2text(mcs,reacID)
    kos = find(mcs<0);
    kis = find(mcs>0);
    if ~isempty(kis) 
        kis = [strjoin(strcat('+',cellstr(reacID(kis,:))),', ') ', '];
    else
        kis = [];
    end
    if ~isempty(kos) 
        kos = strjoin(strcat('-',cellstr(reacID(kos,:))),', ');
    else
        kos = [];
    end
    ko_ki_text = [kis kos];
end

function [mcs, comptime] = MCS_enum_thread(gene_mcs,cnap,modules,...
                                            rkoCost,...
                                            maxSolutions,maxCost,...
                                            gkoCost,...
                                            options,verbose)
        if ~isempty(getenv('SLURM_TEMP'))
            tempdir = getenv('SLURM_TEMP');
        else
            tempdir = eval('tempdir');
        end                              
        wdir = [tempdir num2str(randi([0,1e6-1])) filesep];
        mkdir(wdir);
        filename = [wdir 'ws_' num2str(randi([0,1e6-1])) '.mat'];
        save(filename,'gene_mcs','cnap','modules','rkoCost','maxSolutions','maxCost','gkoCost','options','verbose');
        wd = pwd;
        cd(evalin('base','cnan.cnapath'));
        system(['matlab -nosplash -nodesktop -r "compute_MCS_ext(''' wdir ''',''' filename ''')"']);
        cd(wd);
        load(filename,'mcs','comptime');
        rmdir(wdir,'s');
end
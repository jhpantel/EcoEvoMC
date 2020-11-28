
clear

%% set parameters and initialize variables
% set parameter values to loop through
niches = [6.8 8.5 1.5]; 
olap = [0.5];    % difference in xmaxs between patches; level of overlap between resource distributions
%conx = [0];      % rate of exponential decay in dispersal prob (not implemented in result runs)
drate = [0 0.01 0.1];     % dispersal rates
env_ch = [0 0.00001 0.0004];  % rate of linear environmental change
ICs = [1];    % which set of initial conditions from seed runs to use?

for icz = 1:length(ICs)
for evar = 1:length(olap)
for d = 1:length(drate)
for v = 1:length(env_ch)
for gp = 1:length(niches)
    % load in initial patch communities from separate "seed" files (the results of
    % single patch runs of this code, with seeding=1)
    ikk = num2str(ICs(icz));
    nn = num2str(niches(gp) *10);
    eval(['load N_ICs_nch' nn '_k3_ic' ikk '.mat'])
    NIC = N_ic;
    XIC = X_ic;
    clearvars -except NIC XIC gp niches olap evar drate d env_ch v icz ikk ICs
    
for run = 1:1
    
tic

rng('shuffle');

interval = 200;   % the number of generations that we run the simulation before saving all records, dropping dead weight and continuing

if env_ch(1,v) == 1*10^-5
    gennum = 200002;  % if env change is slow, we want more gens to see full picture of dynamics
else
    gennum=50002;
end

patchnum = 3; 

morphcount = zeros(1,gennum);   % record for number of morphs generated up to generation t

% seeding with newly generated pops?
seeding = 1;
% NB: seeding and filein variables should always have opposite values
% seeding with results from IC runs?
filein = 0;

%parameters for dispersal
mig_thresh = drate(1,d); %dispersal rate 
introduce = 0;  % if introducing dispersal midway through

% parameters for env change
weather = 0;
stormy = env_ch(1,v);

if seeding ==1
    sn = patchnum;
    morphcount(1,1) = patchnum;
end

if filein == 1
    sn = size(XIC,2);   % sn = indicator for how many pops there are in system
    morphcount(1,1) = sn;
end

% create records for the first time interval
traitx = single(zeros(patchnum,sn,(interval+2))); 
snpop = single(zeros(patchnum,sn,(interval+2))); 
dispersed_to = single(zeros(patchnum,sn,(interval+2)));
dispersed_from = single(zeros(patchnum,sn,(interval+2)));
mut_totals = single(zeros(patchnum,sn,(interval+2)));
births = single(zeros(patchnum,sn,(interval+2)));
deathtoll = single(zeros(patchnum,sn,(interval+2)));
totalcomp = single(zeros(patchnum,sn,(interval+2)));

morphs = single(zeros(1,sn));   
species_id = single(zeros(3,sn)); % column 1: trait value; column 2: species id; column 3: patch id

% generate xmax values:

% generating XMAX with olap parameter (above)
lap = olap(1,evar);
if patchnum == 1
    lap = 0;
end
if lap == 0
xmax = zeros(1,patchnum);
else
winglen = (patchnum-1)/2;
leftwing = -winglen*lap;
xmax = [leftwing:lap:-leftwing];
initland = xmax;
end

% since IC runs all had Xmax=0, adjust species traits for new environmental
% (trait) optima, and load into species_id array
if filein == 1
    ini = XIC;
    num_sn_patch = zeros(1,patchnum);
for kc = 1:patchnum
    ini(kc,ini(kc,:)~=0) = ini(kc,ini(kc,:)~=0) + xmax(kc);
    num_sn_patch(1,kc) = length(ini(kc,ini(kc,:)~=0));
end
    idx = find(ini);
    morphs(1,:) = ini(idx);
    init_x = morphs(1,:);     % for easy access later
    species_id(1,1:sn) = morphs;
    species_id(2,1:sn) = 1:sn;
    
    orig_id = zeros(2,length(init_x));
    orig_id(1,:) = init_x;
    tick=1;
    for k = 1:patchnum
        kpops = length(find(ini(k,:)));
        orig_id(2,tick:tick+kpops-1) = k;
        tick=tick+kpops;
    end
    species_id(3,1:sn) = orig_id(2,:);
    
    
    % initialize main records
    traitx(:,:,1) = ini;
    snpop(:,:,1) = NIC;
    init_n = NIC;  % for easy access later 
    
end
    
% ignore section below, for seeding
if seeding ==1     
 % generating xmaxs with olap parameter (above)
 lap = olap(1,evar);
 if patchnum == 1
     lap = 0;
 end
 if lap == 0
 xmax = zeros(1,patchnum);
 else
 winglen = (patchnum-1)/2;
 leftwing = -winglen*lap;
 xmax = [leftwing:lap:-leftwing];
 
 end
 morphs(1,:) = ((rand(1,patchnum)*1)-.5) + xmax;
 init_x = morphs(1,:);
 
 species_id = zeros(1,1);  % species id not applied when seeding==1
 traitx(:,1:patchnum,1) = diag(morphs);
 
 snpop(:,1:patchnum,1) = 500;
 init_n = snpop(:,1:patchnum,1);
end


% creating a matrix of randomized coordinates for each patch on a 10x10
% grid (ignore, if we don't have nariation in dispersal connectivity)
location = zeros(patchnum,2);
location(:,1) = rand(patchnum,1);
location(:,2) = rand(patchnum,1);

Kmax = zeros(1,patchnum);
Kwidth = zeros(1,patchnum);

% store Kmax and Kwidth values
space = zeros(1,patchnum);
for k = 1:patchnum
    Kwidth(1,k) = 1;
    Kmax(1,k) = 10000;
    load = @(x) Kmax(1,k)*exp((-(x-xmax(1,k)).^2/(Kwidth(1,k)).^2));
    space(1,k) = integral(load,-100,100);
end

%calc K(x) for initial pops
Ks = zeros(patchnum,sn);
Ks(:,1:sn) = Kx(Kmax, Kwidth, xmax, patchnum, morphs);


% since lambda = 0, connectivity is divided evenly among all patches (ie,
% no decay with distance)

%calculate probabilistic connectivity between all pairs of patches (i,j) as
%a function of distance between patches -> k by k matrix. Prob distribution
%for patch i of dispersing to all other patches is described by all (i,j)
%cell is in row i of this matrix. It's a discrete conditional distribution
%(conditional on distance between patches). I also save idxpatches, a logical matrix
%indicating all cells {(i,j), i~= j}.

if patchnum > 1

connectivity = zeros(patchnum,patchnum);
distances = zeros(patchnum,patchnum);
expo = zeros(patchnum,patchnum);
idxpatches = zeros(patchnum-1,patchnum-1);

lambda = 0;   % NB: for all sim runs, the pdf of dispersal was equiprobable among patches
if lambda == 0
    binner = 1/(patchnum-1);
    pdfk = 0:binner:1;
    pdf = repelem(pdfk,patchnum,1);

else
for ki = 1:patchnum
    for kj = 1:patchnum
        connectivity(ki,kj) = sqrt((location(ki,1) - location(kj,1))^2 + (location(ki,2) - location(kj,2))^2);%/2;
        distances(ki,kj) = connectivity(ki,kj);
        connectivity(ki,kj) = exp(lambda*(-(connectivity(ki,kj))));
        if connectivity(ki,kj) == 1
            connectivity(ki,kj) = 0;
        end
        expo(ki,kj) = connectivity(ki,kj);
    end
    connectivity(ki,connectivity(ki,:)~=0) = connectivity(ki,connectivity(ki,:)~=0)/sum(connectivity(ki,connectivity(ki,:)~=0));
    idxpatches(ki,:) = find(connectivity(ki,:));
end

pdf = zeros(patchnum,patchnum); % will contain bin "cutoffs". Width of each "bin" is probability of mutant landing on that patch
obn = zeros(patchnum,patchnum-1);  %temp
for k = 1:patchnum
    tempalter = idxpatches(k,:);   % grabs the indexes of alternative patches from k from the earlier initialized matrix
    for bins = 1:patchnum-1
        obn(k,bins) = pdf(k,bins) + connectivity(k,tempalter(bins));  % add connectivity of alternative matrices up (they sum to one) to obtain bin cutoffs
        pdf(k,bins+1) = obn(k,bins);
    end
end

end
end



%parameters for mutation
mutsd = .05;  %mutant trait standard deviation around parent trait
mutthresh = 1*10^-5;  

%parameter for growth rate
r = 1.9;

%parameters for competition function (alpha)
p = 2;  
niche = niches(1,gp)*.1;

if niches(1,gp) == 1.5   % avoiding bug with niche params above
    niche = 1.5;
end

%hard or soft extinction?
hard = 0;

extn = 2; %stochastic extinction threshold
risk = .025; %stochastic extinction risk

%% define useful functions
% Gaussian comp = @(alpha,niche,p) single(exp(-abs(alpha/niche)^p));

% Johansson's alpha function:
johcomp = @(alpha,niche) single(1/(1+((alpha^2)/(2*(niche^2)))));

mut_count = @(n,m_prob) sum(n <= m_prob);

mutation = @(x,mutsd) x*mutsd;

parmut = @(a,b) (a+b);

pars = @(c) cell2mat(c);

mig_count = @(n,d_prob) sum(n <= d_prob);

snnew = single(1);

t = 2; %starting on gen 2
%% large gen loop

for gen = 2:gennum %each gen
    
    
    %% Introduce dispersal, enviro change, or other changes in a particular gen
    
    % introduce dispersal (normally dispersal is present from beginning)
    if gen == 200000
        if introduce == 1
            mig_thresh = drate(1,d);
        end
    end
    
    % if you want to introduce moving trait maxima ("environmental change")
    if gen == 2
        weather = 1; % currently weather is turned on in t=2
        %scape = zeros(gennum-20000,patchnum);
    end
    
    if weather == 1
        % commented out section below creates brownian environmental noise
%     flip = zeros(1,patchnum);   
%     for k = 1:patchnum
%         god = rand;
%         if god > 0.5
%             flip(1,k) = 1;
%         else
%             flip(1,k) = -1;
%         end
%     end
%     delta = stormy * rand(1,patchnum);
%     delta = flip .* delta;
    %scape(gen-+1,:) = xmax + delta;
    %xmax = scape(gen-+1,:);
    
    xmax = xmax + stormy;   % increase all environmental optima by 'stormy'
    Ks(:,1:sn) = Kx(Kmax, Kwidth, xmax, patchnum, morphs);    % update all existing pops' K(x)
    
    end
    
    %% save time record, delete useless info, begin new record
    
    if mod((gen-2),interval) == 0 && gen > 2 %cut off bottom 200 gens
        if patchnum > 1
        idxalltraits = single(find(sum(snpop(:,:,interval+1))));   
        else
        idxalltraits = single(find(snpop(1,:,interval+1)));
        end
        % idxalltraits = index of existing populations in current records
        
        % create temporary records, for after dropping/resaving
        t = 2;
        seedpops = snpop(:,idxalltraits,(interval+1));
        seedtraits = traitx(:,idxalltraits,(interval+1));
        tempKs = Ks(:,idxalltraits);
        tempmorf = morphs(1,idxalltraits);
        seedimig = dispersed_to(:,idxalltraits,(interval+1));
        seedemig = dispersed_from(:,idxalltraits,(interval+1));
        seedmuts = mut_totals(:,idxalltraits,(interval+1));
        seedbrth = births(:,idxalltraits,(interval+1));
        seeddeth = deathtoll(:,idxalltraits,(interval+1));
        seedcomp = totalcomp(:,idxalltraits,(interval+1));
        
        time = num2str(gen-2);
        
        %create var names, with flag for interval, and save as current
        %records (not the temp ones)
        pop = matlab.lang.makeValidName(['oldpops_',time]);
        tra = matlab.lang.makeValidName(['oldtrat_',time]);
        immig = matlab.lang.makeValidName(['oldimig_',time]);
        emmig = matlab.lang.makeValidName(['oldemig_',time]);
        muts = matlab.lang.makeValidName(['oldmuts_',time]);
        birt = matlab.lang.makeValidName(['oldbrth_',time]);
        deth = matlab.lang.makeValidName(['olddead_',time]);
        acomp = matlab.lang.makeValidName(['oldcomp_',time]);
        
        eval([pop '= snpop(:,:,1:interval);']);
        eval([tra '= traitx(:,:,1:interval);']);
        eval([immig '= dispersed_to(:,:,1:interval);']);
        eval([emmig '= dispersed_from(:,:,1:interval);']);
        eval([muts '= mut_totals(:,:,1:interval);']);
        eval([birt '= births(:,:,1:interval);']);
        eval([deth '= deathtoll(:,:,1:interval);']);
        eval([acomp '= totalcomp(:,:,1:interval);']);
        
        ic0 = num2str(icz);
        thru = num2str(run);
        ni = num2str(niches(1,gp));
        ol = num2str(olap(1,evar));
        dispe = num2str(drate(1,d));
        env = num2str(env_ch(1,v));
        pnum = num2str(patchnum);
        filename = matlab.lang.makeValidName(['final_524_nch' ni '_ol' ol '_d' dispe '_env' env]);
        
        % save the records from this interval to output file
        if gen == interval+2
            save(filename);
            save(filename,pop,tra,immig,emmig,muts,birt,deth,acomp,'-append');
        else
            save(filename,pop,tra,immig,emmig,muts,birt,deth,acomp,'-append');
        end

        clearvars snpop traitx morphs 
        clearvars -regexp ^old
        
        
        %create new records, as the temp ones from above
        traitx = single(zeros(patchnum,length(idxalltraits),(interval+2))); 
        snpop = single(zeros(patchnum,length(idxalltraits),(interval+2))); 
        dispersed_to = single(zeros(patchnum,length(idxalltraits),(interval+2)));
        dispersed_from = single(zeros(patchnum,length(idxalltraits),(interval+2)));
        mut_totals = single(zeros(patchnum,length(idxalltraits),(interval+2)));
        births = single(zeros(patchnum,length(idxalltraits),(interval+2)));
        deathtoll = single(zeros(patchnum,length(idxalltraits),(interval+2)));
        totalcomp = single(zeros(patchnum,length(idxalltraits),(interval+2)));
        
        snpop(:,:,1) = seedpops;
        traitx(:,:,1) = seedtraits;
        dispersed_to(:,:,1) = seedimig;
        dispersed_from(:,:,1) = seedemig;
        mut_totals(:,:,1) = seedmuts;
        births(:,:,1) = seedbrth;
        deathtoll(:,:,1) = seeddeth;
        totalcomp(:,:,1) = seedcomp;

        morphs = tempmorf;
        Ks = tempKs;
        
        idxalltraits = single(1:length(idxalltraits));
        
        sn = length(idxalltraits);
        
    end
    
    % idxalltraits = index of existing populations in current records
    if patchnum > 1
        idxalltraits = single(find(sum(snpop(:,:,t-1))));    
        else
        idxalltraits = single(find(snpop(1,:,t-1)));
    end
    
    if snnew == 1
        alph = repmat(morphs(1,idxalltraits),length(morphs(1,idxalltraits)),1)';
        alph = bsxfun(@minus,alph,morphs(1,idxalltraits));
        %alph = arrayfun(@(a) expcomp(a,niche,p), alph);  %GAUSSIAN
        alph = arrayfun(@(a) johcomp(a,niche), alph);
    end
    
    %% Calculate births
    %temporary list of net births for each morph, to use for calculating mutants later
    snbirths = r*snpop(:,idxalltraits,t-1);
    dice = floor(snbirths); 
    snbirths = real(snbirths);
    dice = real(dice);
    log = logical(snbirths);
    
    %% update traits
    % existing traits from last gen are input into current gen
    traitx(:,idxalltraits,t) = traitx(:,idxalltraits,t-1);  
    parrec = traitx(:,idxalltraits,t-1);   % record of 'parents', for calculating species_id
    
    %% Determine if mutants were born
    if isempty(idxalltraits) == 0
    rcount = arrayfun(@rand, dice, log, 'UniformOutput' , false);   % list of rand numbers for each morph-pop, euqal to number of births
    mutnum = cellfun(@(x) mut_count(x,mutthresh), rcount);    % how many mutations are there?
    else
    mutnum = [];
    end
    
    %% If mutants were born, generate mutant lists for each patch-morph-pop
    
    all_traits_added = [];
    tempidx = find(mutnum);
    mutlog = logical(mutnum);

    if isempty(tempidx) == 0
         
    partraits = arrayfun(@repelem, traitx(:,idxalltraits,t), mutnum, 'UniformOutput' , false);  % lists of parent traits (with length=number of mutants per morph-pop)

    mutations = arrayfun(@randn,mutlog,mutnum,'UniformOutput',false);   % generate randomized mutations
    
    mutations = cellfun(@(x) mutation(x,mutsd), mutations, 'UniformOutput', false);   % mutliply mutations by mutsd, scaling them
    
    newtraits = cellfun(@(a,b) parmut(a,b), partraits, mutations,'UniformOutput',false);    % add mutations to parent traits to get mutant traits
    if seeding == 0
    newtraits = cellfun(@(x) no_clones(x,species_id(1,:)), newtraits, 'UniformOutput', false);  % apply no_clones function to make sure no 2 morphs ever have the same trait (doesn't change results)
    end
    new = sum(sum(mutnum));    %number of mutants slash number of new pops
    
    %update big records
    sn = sn + new;   % update sn for new pops
    morphcount(1,gen) = morphcount(1,gen-1) + new;

    % update records for births and mutant numbers
    mut_totals(:,idxalltraits,t) = mutnum(:,:);
    births(:,idxalltraits,t) = snbirths(:,:);
     
%     %update major records
     temp=dispersed_to;
     dispersed_to = single(zeros(patchnum,sn,(interval+2)));
     dispersed_to(:,1:sn-new,:)=temp;
%     
     temp=dispersed_from;
     dispersed_from = single(zeros(patchnum,sn,(interval+2)));
     dispersed_from(:,1:sn-new,:)=temp;
%     
     temp=mut_totals;
     mut_totals = single(zeros(patchnum,sn,(interval+2)));
     mut_totals(:,1:sn-new,:)=temp;
     
     temp=births;
     births = single(zeros(patchnum,sn,(interval+2)));
     births(:,1:sn-new,:)=temp;
     
     temp=totalcomp;
     totalcomp = single(zeros(patchnum,sn,(interval+2)));
     totalcomp(:,1:sn-new,:)=temp;
    
     temp=deathtoll;
     deathtoll = single(zeros(patchnum,sn,(interval+2)));
     deathtoll(:,1:sn-new,:)=temp;
    if seeding == 0
    temp = species_id;
    history = size(species_id,2);
    species_id = zeros(3,history+new);
    species_id(:,1:history)=temp;
    newtree = size(species_id,2);
    end
% add mutant pops to snpop and traitx, as well as storing in species_id
    ticker = 1;
    id = 1;
    for k = 1:patchnum
    
    muts_byk = cell2mat(newtraits(k,:));
    newlength = length(muts_byk);

    if isempty(muts_byk) == 0
    
    numss = mutnum(k,find(mutnum(k,:)));
    parss = parrec(k,find(mutlog(k,:)));
    if seeding == 0
    idk = 1;
    for litter = 1:length(numss)
        litt = numss(litter);
        parid = species_id(2,species_id(1,:)==parss(litter));

        for tree = 1:litt
            species_id(1,newtree-new+id) = muts_byk(idk);
            species_id(2,newtree-new+id) = parid;
            species_id(3,newtree-new+id) = k;
            id = id+1;
            idk = idk+1;
        end
    end
    end
    end
    
    
    snpop(k,(sn-new)+ticker:sn-new+ticker+newlength-1,t) = single(ones(1,newlength)); %add mutant individual pops to snpop

    traitx(k,(sn-new)+ticker:sn-new+ticker+newlength-1,t) = single(muts_byk); 
    
    morphs(1,(sn-new)+ticker:sn-new+ticker+newlength-1) = single(muts_byk);
    
    %where_from(1,(max(morphcount)-new)+ticker:max(morphcount)-new+ticker+newlength-1) = k;
    
    ticker = ticker+newlength;
    
    all_traits_added = horzcat(all_traits_added,muts_byk);
    end
    
    % calculate K(x) for all new mutants
    Ks(:,(sn-new)+1:sn) = Kx(Kmax,Kwidth,xmax,patchnum,all_traits_added);
    
    %subtract mutnum from children (rN) populations
    snbirths = snbirths - single(mutnum);
    
    % mutant population chunk    
    muttz = snpop(:,(sn-new)+1:sn,t);
    
    children = horzcat(snbirths, muttz);
    snbirths = children;
    % snbirths now includes mutant children in separate columns 

        % calculate new idxAll for end of gen t (only populations yet
        % recorded in snpop in gen t are the mutants, so append that index
        % to previous, parent index (idxalltraits)
        
        if patchnum > 1
            shiut = find(sum(snpop(:,:,t)));
            idxalltraits_post = horzcat(idxalltraits,shiut);           
        else
            shiut = find(snpop(1,:,t));
            idxalltraits_post = horzcat(idxalltraits,shiut);           
        end
        
    else
        idxalltraits_post = idxalltraits;
    end
    
    % NB: at this point, idxalltraits_post indexes the children+mutant
    % populations, while idxalltraits indexes the parent populations 
    
    %% compute alpha matrix for existing pops
    if isempty(tempidx) == 0
        alph = repmat(morphs(1,idxalltraits_post),length(morphs(1,idxalltraits_post)),1)';
        alph = bsxfun(@minus,alph,morphs(1,idxalltraits_post));
        %alph = arrayfun(@(a) expcomp(a,niche,p), alph);  %GAUSSIAN
        alph = arrayfun(@(a) johcomp(a,niche), alph);
    end
            
    
    %% calculate competition/density terms    
    sigcomp = alph * snpop(:,idxalltraits_post,t-1)'    ;   % total competition for each existing morph-pop
    dierate = sigcomp ./ Ks(:,idxalltraits_post)'  ;%divide each sigcomp term by corresponding K (denominator in logistic term)
    totalcomp(:,idxalltraits_post,t) = sigcomp'  ;
       
    %% calculate (1-aN/K) values
    perf = ones(length(idxalltraits_post),patchnum);
    survivors = perf - dierate  ;% 1 minus logistic term is the proportion of each existing pop that survives

    
    %% ecology for parent populations    
    if isempty(tempidx) == 0
    snpop(:,idxalltraits,t) = snpop(:,idxalltraits,t-1) + (snbirths(:,1:length(idxalltraits)) .* survivors(1:length(idxalltraits),:)') ;
    
    deathtoll(:,idxalltraits,t) = (snbirths(:,1:length(idxalltraits)) .* dierate(1:length(idxalltraits),:)');   % for death record
    
    else
    % ecology for parent populations, no mutants
    snpop(:,idxalltraits,t) = snpop(:,idxalltraits,t-1) + (snbirths .* survivors');
    morphcount(1,gen) = morphcount(1,gen-1);
    
    deathtoll(:,idxalltraits,t) = (snbirths .* dierate');   % update deaths and births
    births(:,idxalltraits,t) = snbirths(:,:);
    end
    
% delete NaN values (created during the record updates, yet to be filled)
    wick = find(isnan(snpop(:,:,t)));
    if isempty(wick) == 0
        temp = snpop(:,:,t);
        temp(wick) = 0;
        snpop(:,:,t) = temp;
    end
    
    %% Extinctions
    %version 1: growing populations don't die
    if hard == 1
    dead = find(snpop(:,:,t)< extn & snpop(:,:,t) < snpop(:,:,t-1));
    
    temp1 = snpop(:,:,t);
    temp2 = traitx(:,:,t);
    
    temp1(dead) = single(0);
    temp2(dead) = single(0);
    
    snpop(:,:,t) = temp1;
    traitx(:,:,t) = temp2;
    
    else
    
    %version 2: growing populations subject to stochastic extinction
    %w/ a hard threshold at 1 (currently in use)
    
    toast = find(snpop(:,:,t)< extn);
    gone = find(snpop(:,:,t)<1);
    
    fate = rand(1,length(toast));
    
    snip = find(fate < risk);
    
    kill = toast(snip);
    
    temp1 = snpop(:,:,t);
    temp2 = traitx(:,:,t);
    
    temp1(kill) = single(0);
    temp2(kill) = single(0);
    
    temp1(gone) = single(0);
    temp2(gone) = single(0);
    
    snpop(:,:,t) = temp1;
    traitx(:,:,t) = temp2;
    
    end  
    
    
%% Dispersal

if mig_thresh > 0   % if dispersal is happening
 mignum = single(zeros(patchnum,sn)); % each cell (k,i) contains the number of migrants OUT OF morph population i on patch k in this gen
 added = single(zeros(patchnum,sn));  % each cell (k,i) contains the number of migrants ADDED TO morph population i on patch k from ALL other patches.

 % index of existing pops, for dispersal
 if patchnum > 1
 idxallpops = single(find(sum(snpop(:,:,t)))); 
 else
 idxallpops = single(find(snpop(1,:,t)));
 end
 
     pops = floor(snpop(:,idxallpops,t)); %this is a (1,sn) list of all morph populations (incl 0) on patch k
     pops = real(pops);
     miglog = logical(pops);
     npops = length(idxallpops);
     
     rlist = arrayfun(@rand, pops, miglog, 'UniformOutput' , false); %cell array correpsonding to existing morph pops, from left to right in "pops" list. A cell corresponding to a population of size n contains a list of n random (0,1) numbers
     migs = cellfun(@(x) mig_count(x,mig_thresh), rlist); % this is a (1,sn) array containing the total number of migrants from each existing morph pop on patch k in this generation
     
     mignum(:,idxallpops) = migs;
     
     mignumlog = logical(migs);
     
     miglist = arrayfun(@rand, migs, mignumlog, 'UniformOutput' , false);  %cell array corresponding to existing pops, left to right. A cell corresponding to a migrant pool of size m contains a list of random (0,1)s of length m.
   
     % which alt patch will each migrant from each morph-pop go to? (if
     % there is a connectivity distribution)
     pathlist_patch = cellfun(@(x) histcounts(x,pdf(k,:))', miglist','UniformOutput',false);
     
     disper = single(zeros(patchnum,patchnum,npops));
     
     for migpop = 1:npops
         temp1 = cell2mat(pathlist_patch(migpop,:));
         temp2 = zeros(patchnum,patchnum);
         temp2(2:patchnum,:) = temp2(2:patchnum,:) + tril(temp1);
         temp2(1:patchnum-1,:) = temp2(1:patchnum-1,:) + triu(temp1,1);
         disper(:,:,migpop) = single(temp2);
     end
     temp3 = sum(disper,2);
     temp3 = single(reshape(temp3,[patchnum,npops]));
     
     added(:,idxallpops) = temp3;   % this is an array of the number of indivs of species i added to morph-pop i on patch k
     
     dispersed_to(:,:,t) = single(added(:,:));    % update immigration and emigration records
     dispersed_from(:,:,t) = single(mignum(:,:));
    
 trait_filler = repelem(morphs,patchnum,1); % temporary; imagine all traits present in every patch in this gen

 snpop(:,:,t) = snpop(:,:,t) - single(mignum(:,:)); % subtract all emmigrants from all patches
 snpop(:,:,t) = snpop(:,:,t) + single(added(:,:));  % add all immigrants from all patches to other patches

 pl = logical(snpop(:,:,t));  %which populations now exist in which patches?
 traitx(:,:,t) = single(trait_filler .* pl);   %update trait record according to which populations now exist in which patches.

end
    
    % calculate new idxAll for end of gen t
    if patchnum > 1
        idxalltraits_end = find(sum(snpop(:,:,t)));
    else
        idxalltraits_end = find(snpop(1,:,t)); % ???
    end


%% if mutants added, no new idxAll or alpha. If mutants added, new idxAll and alpha
if length(idxalltraits_post) == length(idxalltraits_end) && sum(idxalltraits_post == idxalltraits_end) == length(idxalltraits_post)
    snnew = 0;
    clearvars added all_traits_added dead dierate disper idxallpops ...
        miglist miglog mignum mignumlog migs mutations mutlog ...
        mutnum newtraits partraits pathlist_patch perf please pops present2 ...
        present_post rcount rlist sigcomp snbirths survivors temp1 temp2 temp3...
        trait_filler
        
else %clear everything including alpha idxalltraits 
    snnew = 1;
    clearvars added all_traits_added alpha dead dierate disper idxallpops ...
        idxalltraits miglist miglog mignum mignumlog migs mutations mutlog ...
        mutnum newtraits partraits pathlist_patch perf please pops present2 ...
        present_post rcount rlist sigcomp snbirths survivors temp1 temp2 temp3
end

t = t+1;

if mod(gen,2000) == 0
gen
end

end

telap = toc;
save(filename,'morphcount','telap','init_x','init_n','species_id','-append')

filename

clearvars -except run gp niches olap drate evar d env_ch v NIC XIC icz ikk ICs

end
end
end
end
end
end
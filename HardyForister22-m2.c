//a SLiM model of the joint evolution of quantitative habitat-specific performance and preference. In this version pleiotropies affect habitat preference as well as performance


initialize() {
	//##Set Some Run Parameters ##################################################
	defineConstant("Selection_strength", ss);
	defineConstant("RunId", rn);
	prefRat = (pp + 1) * pr;
	defineConstant("Mutation_mix", c(pp,1,1,prefRat));
	defineConstant("Mutation_var",c(mv));
	defineConstant("delto", delo); //rate at which performance optima change
	defineConstant("K1", 500); // host 1 carrying capacity
	defineConstant("K2", 500); //host 2 carrying capacity

	// QTL-related constants used below
     defineConstant("QTL_mu", c(0, 0, 0)); // pleio effect means
    defineConstant("QTL_cov", c(cv)); // |max out cov |started with 0.25 pleio effect covariance
    defineConstant("QTL_sigma", matrix(c(mv,QTL_cov,QTL_cov,QTL_cov,mv,QTL_cov,QTL_cov,QTL_cov,mv), nrow=3));
    //Fixed optima for each host 
    defineConstant("o1", 0.0);
    defineConstant("o2", 0.0);
	//############################################################################

	initializeSLiMModelType("nonWF");
	initializeMutationRate(1e-7); //7
	initializeMutationType("m1", 0.5, "f", 0.0);   // neutral mutations
	initializeMutationType("m2", 0.5, "n", 0.0, 0.0);   // Pleiotropic QTLs
	initializeMutationType("m3", 0.5, "n", 0.0, 0.0);  // trait 1 CN QTLs
	initializeMutationType("m4", 0.5, "n", 0.0, 0.0);  // trait 2 CN QTLs
	initializeMutationType("m5", 0.5, "n", 0.0, Mutation_var); // QTLs for preference
	m1.convertToSubstitution = T;
	m2.convertToSubstitution = F;
	m3.convertToSubstitution = F;
	m4.convertToSubstitution = F;
	m5.convertToSubstitution = F;
	m2.mutationStackPolicy = "l";
	m3.mutationStackPolicy = "l";
	m4.mutationStackPolicy = "l";
	m5.mutationStackPolicy = "l";

	//##Set the mix of Pleiotropic and CN QTLs####################################
	// g1 is for neutral sites. g2 is for performance and preference QTLs
	initializeGenomicElementType("g1", c(m1), c(1.0));
	initializeGenomicElementType("g2", c(m2,m3,m4,m5), Mutation_mix);
	//############################################################################

	// chromosome of length 100 kb with 10 QTL regions interspered with neutral regions
    for (index in 1:10){
        initializeGenomicElement(g1, index*10000, index*10000 + 4999);
		initializeGenomicElement(g2, index*10000 - 5000, index*10000 - 1);
	}
	// Set recombination rate
	initializeRecombinationRate(1e-4); //this had been set to 1e-8
}

// Kill the non-neutrallness of the preference QTLs
// We don't need to do this for performance QTLs, because we code our own
// mutation generators for those. (Lets us log all mutations at those sites.)
fitness(m5) { return 1.0; } 

//In non-WF models we need to explicitly model reproduction
reproduction(){
    litter_lambda = 1.5; //lots of fun stuff with 2.0
    litterSize = rpois(1, litter_lambda);
    partner = subpop.sampleIndividuals(1);
    for (j in seqLen(litterSize)){
        subpop.addCrossed(individual, partner);
        }
}

// Specifify performance QTL mutation effect generators
mutation(m2){
	// draw mutation effect for new m2 pleiotropic mutation
	Peffect = rmvnorm(1, QTL_mu, QTL_sigma);// * Mutation_var; this to scale 
    mut.setValue("ep0", Peffect[0]);
    mut.setValue("ep1", Peffect[1]);
    mut.setValue("ep2", Peffect[2]);
	Peffects = sim.getValue("all_P_effects");
	sim.setValue("all_P_effects", rbind(Peffects, Peffect));
	return T;
}

mutation(m3){
	// draw mutation effect for new m3 conditionally neutral mutation
	CNeffect = rnorm(1, 0, Mutation_var);
	mut.setValue("ecn", CNeffect);
	CN3effects = sim.getValue("all_CN3_effects");
	sim.setValue("all_CN3_effects", c(CN3effects, CNeffect));
	return T;
}

mutation(m4){
   // draw mutation effect for new m4 conditionally neutral mutation
	CNeffect = rnorm(1, 0, Mutation_var);
	mut.setValue("ecn", CNeffect);
	CN4effects = sim.getValue("all_CN4_effects");
	sim.setValue("all_CN4_effects", c(CN4effects, CNeffect));
	return T;
}

// Set up two subpopulations, specify their starting sizes
1 early() {
	sim.addSubpop("p1", 500);
	sim.addSubpop("p2", 500);
	sim.setValue("p1_opt", o1);
	sim.setValue("p2_opt", o2);
}

//Do habitat choice as a function of m5 and m2 mutations
//To keep things speedy all of the calculations are vectorized
early(){
    inds = sim.subpopulations.individuals;
    inds_m5 = inds.sumOfMutationsOfType(m5);
    muts2 = inds.genomes.mutationsOfType(m2);
    pleioPref = size(muts2) ? sum(muts2.getValue("ep2")) else 0.0;
    preffy = inds_m5 + pleioPref;
    
    //Here's a new implementation that's not as prone to shutting down migration
    //we'll skip the inertia thing
    pref_p1 = 1/(1+exp(-1*preffy));
    choice = ifelse(runif(inds.size()) < pref_p1, 1, 2);
    //preference only gets turned in to migration if an individual is not where they'd
    //prefer to be
    moving = inds[choice != inds.subpopulation.id];
    from_p1 = moving[moving.subpopulation == p1];
    from_p2 = moving[moving.subpopulation == p2];
    p2.takeMigrants(from_p1);
    p1.takeMigrants(from_p2);
}


//Scale fitness by performance phenotypes.
//In non-WF models fitness affects survival.
1000: early() {
    p1_opt = sim.getValue("p1_opt");
    p2_opt = sim.getValue("p2_opt");
    Do = rexp(1,delo);
    p1_new_opt = p1_opt + Do; //0.1;//rnorm(1,0,0.1); //0.05 will sink em
    p2_new_opt = p2_opt + Do; //0.1;//rnorm(1,0,0.5);
    sim.setValue("p1_opt", p1_new_opt);
    sim.setValue("p2_opt", p2_new_opt);
    
    inds = sim.subpopulations.individuals;
    inds[inds.age > 0].fitnessScaling = 0.0; // non-overlapping generations (buggier)
	for (ind in inds){
		// construct phenotypes from additive effects of QTL mutations
		muts2 = ind.genomes.mutationsOfType(m2);
		muts3 = ind.genomes.mutationsOfType(m3);
		muts4 = ind.genomes.mutationsOfType(m4);
		baseP1 = size(muts2) ? sum(muts2.getValue("ep0")) else 0.0;
		baseP2 = size(muts2) ? sum(muts2.getValue("ep1")) else 0.0;
		cnP1 = size(muts3) ? sum(muts3.getValue("ecn")) else 0.0;
		cnP2 = size(muts4) ? sum(muts4.getValue("ecn")) else 0.0;
		phenotype1 = baseP1 + cnP1; //the sum of P and CN mutations
		phenotype2 = baseP2 + cnP2;
		ind.setValue("phenotype1", phenotype1);
		ind.setValue("phenotype2", phenotype2);

		// calculate fitness effects
	    scale = dnorm(p1_opt, p1_opt, Selection_strength[0]); //quantiles, mean, sd
		effect0 = dnorm(p1_opt, phenotype1, Selection_strength[0]) / scale;
		effect1 = dnorm(p2_opt, phenotype2, Selection_strength[0]) / scale;
		
		if (ind.subpopulation.id==1){
		    ind.fitnessScaling = effect0;
		    ind.tagF = effect0;
		    ind.setValue("GoI", effect1);
		} else {
		    ind.fitnessScaling = effect1;
		    ind.setValue("GoI", effect0);
		    ind.tagF = effect1;
		}
	}
}


1: early() {
    p1.fitnessScaling = min((K1 / p1.individualCount), 1.0);
    p2.fitnessScaling = min((K2 / p2.individualCount), 1.0);
}

//Start printing outputs after 200 generations
1000:2000 late() {
    inds = sim.subpopulations.individuals;
	if (!size(inds)){
        sim.simulationFinished();
    } else {
		mp1P1 = size(p1.individuals) ? mean(p1.individuals.getValue("phenotype1")) else 0.0;
		mp2P2 = size(p2.individuals) ? mean(p2.individuals.getValue("phenotype2")) else 0.0;
		
		mp1 = mean(inds.getValue("phenotype1"));
		mp2 = mean(inds.getValue("phenotype2"));
		
		inds_m5_p1 = p1.individuals.sumOfMutationsOfType(m5);
        muts2_p1 = p1.genomes.mutationsOfType(m2);
        pleioPref_p1 = size(muts2_p1) ? sum(muts2_p1.getValue("ep2")) else 0.0;
        preffy_p1 = inds_m5_p1 + pleioPref_p1;
        inds_m5_p2 = p2.individuals.sumOfMutationsOfType(m5);
        muts2_p2 = p2.genomes.mutationsOfType(m2);
        pleioPref_p2 = size(muts2_p2) ? sum(muts2_p2.getValue("ep2")) else 0.0;
        preffy_p2 = inds_m5_p2 + pleioPref_p2;
		mPrefP1 = size(inds_m5_p1) ? mean(preffy_p1) else 0.0;
		mPrefP2 = size(inds_m5_p2) ? mean(preffy_p2) else 0.0;
		pref_p1 = 1/(1+exp(-1*mPrefP1));
		pref_p2 = 1/(1+exp(-1*mPrefP2));
		Np1 = size(p1.individuals);
		Np2 = size(p2.individuals);
			
		cat("-------------------------------\n");
      	cat("Output for end of generation " + (sim.generation - 1) + ":\n\n");
       		
       	// Output population additive QTL-based phenotypes
       	catn("p1 mean phenotype = " + mp1);
       	catn("p2 mean phenotype = " + mp2);
       	catn("p1 mean preference = " + mPrefP1 + " | " + pref_p1);
       	catn("p2 mean preference = " + mPrefP2 + " | " + pref_p2);
       	catn("p1 mean fitness = " + mean(p1.individuals.tagF)); 
       	catn("p2 mean fitness = " + mean(p2.individuals.tagF));
       	catn("p1 size = " + Np1);
       	catn("p2 size = " + Np2);
	    
	    //bWd is an index of specialization:
	    bWd = mean(inds.tagF) - mean(inds.getValue("GoI"));
	    catn("bWd = " + bWd);
	    H = calcHeterozygosity(inds.genomes,muts=c(inds.genomes.mutationsOfType(m2),inds.genomes.mutationsOfType(m3),inds.genomes.mutationsOfType(m4)));
	    //calc Shanon's
	    //Np1 = size(p1.individuals);
        //Np2 = size(p2.individuals);
        Nall = size(inds);
        pi = (Np1+0.0001)/Nall;
        pj = (Np2+0.0001)/Nall;
        Hp = -1 * (pi*log(pi) + pj*log(pj));
        Shan = Hp; 
        mW = mean(inds.tagF);
        //fst = calcFST(p1.genomes, p2.genomes);
        writeFile("Highway.csv", paste(c(mp1,mp2,mPrefP1,mPrefP2,Np1,Np2,bWd,(sim.generation-1),ss,mv,pp,pr,cv,H,Shan,delto,mW,RunId),sep=','), append=T);
	}
}



/* epistatic.c */

/* ***************************************************** */

#include "libhdr"
#include "ranlib.h"
#define NN 5001  /* max number of NIND and gen */
#define MM 2001  /* max number of NCRO */
#define maxmpt 5001
#define normalthreshold 30

int muts, mutsL, lastinmutantspoissontable, lastinmutantspoissontableL, lastinrecombinantpoissontable;
int NMAX, NIND[NN], NCRO, NLOCI, TOTLOCI, reps_alive_gen[101], NSEGLOCNP;
int gen, tgen, generations, i, j, k, l, m, rep, replicates, born, dead;
int mutants_ocurred, RM[NN], initialgenNP[MM][31], initialgen[MM][31], p1, p2;
int ran_i, ran_j, ran_k, ran_l, ran_h, countNoSS, NoSS_k[60001], NoSS_l[60001];
int gm[NN][MM][2], gmo[NN][MM][2], gmb[NN][MM][2], sm[NN][MM][2];
int NINDNP, gmNP[NN][MM][2], dom, C, scaling;
int a, b, n, nn[NN], father[NN], mother[NN], assign;
int family[NN], maxfec[NN], count, neutral;
int tagPS[NN], tagFS[NN], tagHS[NN], failed[NN], mark[NN];
int fathero[NN], mothero[NN], tagPSo[NN], tagFSo[NN], tagHSo[NN];
int episK, episType, num_epis_loci, num_seg_loci;
int n1, n2, nsup, min, max;
int num_homo_whole, homo[MM][31], num_homo[MM];


double genval_sf[NN], genval_sv[NN], pm_sf[NN], pm_sv[NN], pm_sfb[NN], pm_svb[NN];
double Lambda_s, Lambda_L, MAXFEC, DS_MAXFEC;
double L, w, AA, Aa, aa, q[MM][31], qNP[MM][31], LEQ_NP, LEQ_NPL, LEQ_NPf, LEQ_NPv, TLEQ_NP, D_NP;
double ave_s, addedfixed_sNP, addedfixed_sf, addedfixed_sv, alpha_s, beta_s, k_s, ave_hs, zs;
double mutantspoissontable[maxmpt], mutantspoissontableL[maxmpt], recombinantpoissontable[maxmpt];
double sNP[MM][31], hsNP[MM][31], s[MM][31], hs[MM][31];
double H, PS, leqLv[NN], leqNLv[NN], leqLf[NN], leqNLf[NN], leqT[NN];
double smat[NN][NN], mat[NN][NN];
double sk2, ne, fvk, relaxedsel;
double cum_pm_sf[NN], Wfo, Wvo;
double numPSborn, numFSborn, numHSborn;
double sup[NN][100], g_sf_sup[NN][100], g_sv_sup[NN][100], W[NN][100];
double sEPISmin;

int fecunditysearch();
void disorder_NoSS();

struct acc gmean_sf_zygotes[NN], gmean_sv_zygotes[NN];
struct acc mean_N[NN], gmean_sf[NN], gmean_sv[NN], SK2[NN], Hw[NN], fp[NN], Fp[NN];
struct acc LEQ_Lf[NN], LEQ_Lv[NN], LEQ_NLf[NN], LEQ_NLv[NN], LEQ_T[NN], TIME_EXT, P_PSBORN[NN], P_FSBORN[NN], P_HSBORN[NN];
struct acc NeVk[NN], Fvk[NN], num_EPIS[NN], num_SEG[NN];
struct acc mean_N[NN], gmean_sf[NN], gmean_sv[NN], SK2[NN], Hw[NN], fp[NN], Fp[NN];
struct acc gmean_sf_sup[NN], gmean_sv_sup[NN], PROD[NN];
struct acc AVE_num_homo[NN], AVE_num_homo_whole[NN];
struct acc FREQ_HOMO[NN];

/* ********************** Distribution of genes  *********************** */

struct acc q_ng_00[101], q_ng_00_01[101], q_ng_01_02[101], q_ng_02_03[101];
struct acc q_ng_03_04[101], q_ng_04_05[101], q_ng_05_06[101], q_ng_06_07[101];
struct acc q_ng_07_08[101], q_ng_08_09[101], q_ng_09_10[101], q_ng_10[101];

struct acc s_ng_0[101], s_ng_0_106[101], s_ng_106_104[101], s_ng_104_102[101], s_ng_102_101[101];
struct acc s_ng_01_02[101], s_ng_02_04[101], s_ng_04_06[101], s_ng_06_10[101], s_ng_10[101];

struct acc s_ng_0_fix[101], s_ng_0_106_fix[101], s_ng_106_104_fix[101], s_ng_104_102_fix[101], s_ng_102_101_fix[101];
struct acc s_ng_01_02_fix[101], s_ng_02_04_fix[101], s_ng_04_06_fix[101], s_ng_06_10_fix[101], s_ng_10_fix[101];

struct acc s_ng_0_lost[101], s_ng_0_106_lost[101], s_ng_106_104_lost[101], s_ng_104_102_lost[101], s_ng_102_101_lost[101];
struct acc s_ng_01_02_lost[101], s_ng_02_04_lost[101], s_ng_04_06_lost[101], s_ng_06_10_lost[101], s_ng_10_lost[101];

struct acc hs_ng_00_02[101], hs_ng_02_04[101], hs_ng_04_06[101], hs_ng_06_08[101], hs_ng_08_10[101];

struct acc hs_ng_00_02_fix[101], hs_ng_02_04_fix[101], hs_ng_04_06_fix[101], hs_ng_06_08_fix[101], hs_ng_08_10_fix[101];

struct acc hs_ng_00_02_lost[101], hs_ng_02_04_lost[101], hs_ng_04_06_lost[101], hs_ng_06_08_lost[101], hs_ng_08_10_lost[101];

FILE *fptr, *fgen, *frep, *fdat, *fpop, *fout, *fped, *fWrep, *fSAL;

/* ***************************************************** */

main()
{
	fgen = fopen ("genfile.dat","w");
	fSAL = fopen ("salida.dat","w");
//	fout = fopen ("outfile.dat","w");
//	fped = fopen ("pedfile.dat","w");
	fWrep = fopen ("Wreplicates","w");


	getinputs();
	recombination_masks();
	natural_population();

	if (tracelevel!=0) 	fptr = fopen ("dfilename.dat","w");

	for (rep=1; rep<=replicates; rep++)
	{
		// counter for pedfile
		n = 0;

		frep = fopen ("repfile.dat","a");
		fprintf (frep,"replicate %d\n", rep);
		fclose(frep);

		if (tracelevel!=0) fprintf (fptr,"\n***** Replicate %d *****\n", rep);

		natural_data();
		sample();
		NIND[0] = NMAX;

		for (gen=0; gen<=generations; gen++)
		{
			if (NIND[gen] == 1)	break;

			if (tracelevel!=0) fprintf (fptr,"\n***** Generation %d *****\n\n", gen);
			if (tracelevel!=0)	fprintf (fptr,"NIND = %d  MAXFEC = %f  DS_MAXFEC = %f\n", NIND[gen], MAXFEC, DS_MAXFEC);

			if (NIND[gen] == 0)
			{
//				fprintf(fgen,"%d ", gen);
				accum(&TIME_EXT, gen);
				break;
			}

			poisson_tables();
			accum(&mean_N[gen], (double)NIND[gen]);

//			mutation_neutral ();
			mutation_selected ();
//			if (tracelevel!=0) dumpoffspringaftermutation();
			neutral_genes ();
			selected_genes ();
			genotypic_values();
			if ((episK > 1)&&(episType==1))		genotypic_values_epistasis_ncro();
			if ((episK > 1)&&(episType==2))		genotypic_values_epistasis_whole();
			if (gen%tgen == 0) 	if (NIND[gen] != 0)
			{
				reps_alive_gen[gen/tgen] ++;
				distribution_qsh();
			}
			phenotypeB();
			if (tracelevel!=0) dumpphenotypes();

			disorder_parents();

			coancestry_matrix();
			dumpparents();

			fecundities();
			mating_RC();

			calculations();

//			if (tracelevel!=0) dumpoffspring();
		}
  	}

	for (rep=1; rep<=replicates; rep++)
	if (sup[rep][generations]==1.0) 
	{
		nsup++;

		for (gen=0; gen<=generations; gen++)
		{
			accum (&gmean_sf_sup[gen], g_sf_sup[rep][gen]);
			accum (&gmean_sv_sup[gen], g_sv_sup[rep][gen]);
			if (gen==0)	W[nsup][gen] = (g_sv_sup[rep][gen]);
			else		W[nsup][gen] = (g_sf_sup[rep][gen-1] * g_sv_sup[rep][gen]);

//			printf("nsup=%d rep=%d gen=%d fec=%f via%f W=%f\n", nsup, rep, gen, g_sf_sup[rep][gen], g_sv_sup[rep][gen], W[nsup][gen]);
		}
	}

//	ORDER FITNESSES
	for (gen=0; gen<=generations; gen++)
	for (n1=1; n1<=nsup; n1++)
	for (n2=n1; n2<=nsup; n2++)
	{
		if (W[n1][gen] > W[n2][gen])
		{
			w=W[n1][gen];
			W[n1][gen]=W[n2][gen];
			W[n2][gen]=w;
		} 
	}

//	printf("\n");
//	for (gen=0; gen<=generations; gen++)
//	for (n1=1; n1<=nsup; n1++)
//	{
//		printf("gen=%d nsup=%d W=%f\n", gen, n1, W[n1][gen]);
//	}

//	for (gen=0; gen<=generations; gen++)
//	printf("\ngen=%d  min=%d max=%d Wmin=%f  Wmax=%f\n", gen, (int)(nsup*0.025), (int)(nsup*0.975), W[(int)(nsup*0.025)][gen], W[(int)(nsup*0.975)][gen]);


	printout();
	distribution_out();

	writeseed();
}

/* ***************************************************** */

getinputs()
{
	tracestart();
	getseed();
	getintandskip("NMAX (max 1000):",&NMAX,1,1000);
	getrealandskip("MAXFEC per couple:",&MAXFEC,1.0,100.0);
	getrealandskip("Stand. Dev. MAXFEC:",&DS_MAXFEC,0.0,100.0);
	getrealandskip("Probability of partial selfing (99:Random): ",&PS,0.0,99.0);
	getrealandskip("Length of genome in Morgans (99:FreeRecom) :",&L,0.0,99.0);
	getintandskip("NCRO (min 1, max 2000):",&NCRO,1,2000);
	getintandskip("NLOCI (first is neutral)(min 2, max 30):",&NLOCI,2,30);
	TOTLOCI = NCRO * NLOCI;
	getrealandskip("Lambda_s :",&Lambda_s,0.0,(double)infinity);
	getrealandskip("Lambda_L (s=1,h=0.02):",&Lambda_L,0.0,(double)infinity);

	getrealandskip("Beta_s :",&beta_s,0.0,(double)infinity);

	getrealandskip("Average |s| :",&ave_s,0.0,1.0);
	alpha_s = beta_s / ave_s;

	getintandskip("dom (constant 0, variable 1):",&dom,0,1);
	getrealandskip("Average h:",&ave_hs,0.0,(double)infinity);
	k_s = alpha_s * (pow((2.0*ave_hs),((-1.0)/beta_s))-1.0);

	getrealandskip("Relaxation factor (0: no relaxation; 1: full relaxaton) :",&relaxedsel,0.0,1.0);
	getintandskip("Neutral model :",&neutral,0,1);
	getintandskip("Scaling (0: no; 1: yes) :",&scaling,0,1);
	getintandskip("Epistatic effect (episK=1: no epistasis; 2: effect squared; 3: effect^Nhomo) :",&episK,1,3);
	getintandskip("Epistatic type (episType=1: ncro; 2: whole) :",&episType,1,2);
     getrealandskip("s EPIS min (0.0-1.0) :",&sEPISmin,0.0,1.0);

	getintandskip("Number of generations :",&generations,1,1000);
	getintandskip("Number of replicates :",&replicates,1,infinity);

	tgen = 1;
}

/* **************************************************** */

recombination_masks ()
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
}

/* ***************************************************** */

natural_population ()
{
	int g0, g1, dinitialgen, dN;
	
	double dadds, dadda, ds, da, dhs, dha;

	/* ***** take genotypic values of natural population ***** */

	fpop=fopen("popfile","r");

	fscanf(fpop,"%d", &dN);
	NINDNP = dN;

	for (i=0; i<dN; i++)
	for (k=0; k<NCRO; k++)
	{
		fscanf(fpop,"%d%d", &g0, &g1);
		gmNP[i][k][0] = g0;
		gmNP[i][k][1] = g1;
	}

	fclose(fpop);

	/* ***** take effects of genes ***** */

	fdat=fopen("datafile","r");

	fscanf(fdat,"%lf%lf", &dadds, &dadda);
	addedfixed_sNP = dadds;
	addedfixed_sNP = addedfixed_sNP + ((1.0 - addedfixed_sNP) * relaxedsel);

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		fscanf(fdat,"%lf%lf%lf%lf%d", &ds, &da, &dhs, &dha, &dinitialgen);
		sNP[k][l] = ds;
		if (sNP[k][l] == (-1))	sNP[k][l] = sNP[k][l];
		else				sNP[k][l] = sNP[k][l] - (sNP[k][l] * relaxedsel);
		hsNP[k][l] = dhs;
		
		initialgenNP[k][l] = dinitialgen;
	}

	/* ***** estimate LEQ in the natural population ***** */

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NINDNP; i++)
		{
			if (((gmNP[i][k][0] & RM[l])==RM[l])&&((gmNP[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gmNP[i][k][0] & RM[l])!=RM[l])&&((gmNP[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNP[k][l] = (aa/(double)NINDNP)+(Aa/(2.0*(double)NINDNP));
		LEQ_NP += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
		if (sNP[k][l] == (-1.0))	LEQ_NPL += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);

			if (l%2==0)	LEQ_NPf += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
			else		LEQ_NPv += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);

		TLEQ_NP += -(sNP[k][l] * qNP[k][l]);
		if (qNP[k][l] > 0.0)
		{
			D_NP += (sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0);
			NSEGLOCNP ++;
		}
	} 

//	COMMENT
//	printf("\n LEQ_NP = %f", LEQ_NP);
//	for (k=0; k<NCRO; k++)	for (l=1; l<NLOCI; l++)
//	printf("\n k=%d l=%d sNP=%f  hNP=%f inigen=%d", k, l, sNP[k][l], hsNP[k][l], initialgenNP[k][l]);

	fclose(fdat);
}

/* ***************************************************** */

natural_data()
{
	addedfixed_sf = sqrt(addedfixed_sNP);
	addedfixed_sv = sqrt(addedfixed_sNP);

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		s[k][l] = sNP[k][l];
		hs[k][l] = hsNP[k][l];
		
		initialgen[k][l] = initialgenNP[k][l];
//COMMENT
		if (tracelevel!=0)    if ((k==0)&&(l<=10))    fprintf(fptr,"\n k=%d l=%d s=%f  h=%f  inigen=%d", k, l, s[k][l], hs[k][l], initialgen[k][l]);
	}
}

/* ***************************************************** */

sample ()
{
	int g;
	
	if (NMAX == 2)
	{
		ran_i = (int)(uniform() * NINDNP);
		do
		{
			ran_j = (int)(uniform() * NINDNP);
		}
		while (ran_j == ran_i);

		for (k=0; k<NCRO; k++)
		{
			g=gmNP[0][k][0]; gmNP[0][k][0]=gmNP[ran_i][k][0]; gmNP[ran_i][k][0]=g;
			g=gmNP[0][k][1]; gmNP[0][k][1]=gmNP[ran_i][k][1]; gmNP[ran_i][k][1]=g;

			g=gmNP[1][k][0]; gmNP[1][k][0]=gmNP[ran_j][k][0]; gmNP[ran_j][k][0]=g;
			g=gmNP[1][k][1]; gmNP[1][k][1]=gmNP[ran_j][k][1]; gmNP[ran_j][k][1]=g;
		}
	}
	else
	{
		for (i=0; i<NMAX; i++)
		{
			ran_i = (int)(uniform() * NINDNP);

			for (k=0; k<NCRO; k++)
			{
				g=gmNP[i][k][0]; gmNP[i][k][0]=gmNP[ran_i][k][0]; gmNP[ran_i][k][0]=g;
				g=gmNP[i][k][1]; gmNP[i][k][1]=gmNP[ran_i][k][1]; gmNP[ran_i][k][1]=g;
			}
		}
	}
	for (i=0; i<NMAX; i++)
	for (k=0; k<NCRO; k++)
	{
		gm[i][k][0]=gmNP[i][k][0];
		gm[i][k][1]=gmNP[i][k][1];
	}
}

/* ***************************************************** */

poisson_tables ()
{   
    /* SELECTED LOCI WITH POISSON (2NL) NEW MUTATIONS */

    if ( (exp(-2.0*(double)NIND[gen]*Lambda_s) != 0.0)&&(2.0*(double)NIND[gen]*Lambda_s < normalthreshold) )
    generatepoissontable(2.0*(double)NIND[gen]*Lambda_s, &lastinmutantspoissontable, mutantspoissontable, maxmpt-1);

    if ( (exp(-2.0*(double)NIND[gen]*Lambda_L) != 0.0)&&(2.0*(double)NIND[gen]*Lambda_L < normalthreshold) )
    generatepoissontable(2.0*(double)NIND[gen]*Lambda_L, &lastinmutantspoissontableL, mutantspoissontableL, maxmpt-1);

  /* NUMERO DE RECOMBINACIONES POISSON CON MEDIA L*/

    if ( (exp(-L) != 0.0) && (L < normalthreshold) )
    generatepoissontable(L, &lastinrecombinantpoissontable, recombinantpoissontable, maxmpt-1);
}

/* ***************************************************** */

mutation_neutral ()
{
    /* NEUTRAL GENES: (POISSON) 2N(Lambda_a) NEW MUTATIONS (TWO DIRECTIONS) */

    muts = mutationnumber();

    if (tracelevel!=0)    fprintf(fptr,"\n New neutral mutants = %d\n", muts);

    for (m=0; m<muts; m++)
    {
		ran_i = (int)(uniform()*NIND[gen]);
		ran_k = (int)(uniform()*NCRO);
		ran_h = (int)(uniform()*2.0);
		if ( (gm[ran_i][ran_k][ran_h] & RM[0])==RM[0] )
		    gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] & (~RM[0]));
		else	gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] | RM[0]);
    }
}

/* ***************************************************** */

void mutation_selected ()
{
	/* SELECTED GENES: (POISSON) 2N(Lambda_s + Lambda_L) NEW MUTATIONS */

	muts = mutationnumber();
	mutsL = mutationnumberL();

	if ( (muts + mutsL) == 0 ) goto label;

	if (tracelevel!=0)    fprintf(fptr,"\nNew selected mutants  muts=%d  mutsL=%d\n", muts, mutsL);

	countNoSS = 0;

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		if (q[k][l]==0.0)
		{
			countNoSS += 1;
			NoSS_k[countNoSS-1] = k;
			NoSS_l[countNoSS-1] = l;
		}
	}
   
//	if (tracelevel!=0)    fprintf(fptr,"\n countNoSS=%d\n", countNoSS);

	mutants_ocurred = 0;

	if (countNoSS != 0)
	{
		disorder_NoSS (NoSS_k,NoSS_l);

		for (m=0; m<countNoSS; m++)
		{
			if (mutants_ocurred==(muts+mutsL))    goto label;

			ran_i = (int)(uniform()*NIND[gen]);
   			ran_h = (int)(uniform()*2.0);
    			gm[ran_i][NoSS_k[m]][ran_h]=(gm[ran_i][NoSS_k[m]][ran_h] | RM[NoSS_l[m]]);
		   	mutants_ocurred += 1;
			initialgen[NoSS_k[m]][NoSS_l[m]] = gen;

			/* ****** Lethal mutations ****** */

			if(mutants_ocurred <= mutsL)
			{
				s[NoSS_k[m]][NoSS_l[m]] = (-1.0);
				hs[NoSS_k[m]][NoSS_l[m]] = 0.02;
			}
			else
			{
	     			/* ****** values of s, hs ****** */

	    			zs = gengam (alpha_s, beta_s);
				if (zs > 1.0)	zs=1.0;
				s[NoSS_k[m]][NoSS_l[m]] = (-zs) - ((-zs) * relaxedsel);
				if (dom==0)	hs[NoSS_k[m]][NoSS_l[m]] = ave_hs;
				else		hs[NoSS_k[m]][NoSS_l[m]] = uniform() * exp(k_s*s[NoSS_k[m]][NoSS_l[m]]);
			}
//			if (tracelevel!=0)    fprintf(fptr,"alpha_s=%f  beta_s=%f\n", alpha_s, beta_s);
//			if (tracelevel!=0)    fprintf(fptr,"s=%f\n", s[NoSS_k[m]][NoSS_l[m]]);
			if (tracelevel!=0)    fprintf(fptr,"ran_i=%d  k=%d  l=%d   ran_h=%d\n", ran_i, NoSS_k[m],  NoSS_l[m], ran_h);
		}
	}

	for (m=mutants_ocurred; m<muts; m++)
	{
		ran_i = (int)(uniform()*NIND[gen]);
		ran_k = (int)(uniform()*NCRO);
		do {ran_l = (int)(uniform()*NLOCI);}   while (ran_l==0);
		ran_h = (int)(uniform()*2.0);

//		fprintf(frep,"recurrent mutation\n");
		
		gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] | RM[ran_l]);

//		if (tracelevel!=0)    fprintf(fptr,"ran_i=%d  ran_k=%d  ran_l=%d  ran_h=%d\n", ran_i, ran_k, ran_l, ran_h);
	}
	label: /* end of mutations */;
}

/* ***************************************************** */

int mutationnumber()
{
	int r;
	if ((2.0*(double)NIND[gen]*Lambda_s < normalthreshold) && (exp(-2.0*(double)NIND[gen]*Lambda_s) != 0.0) )
	{
		r = poisson(lastinmutantspoissontable, mutantspoissontable);
	}
	else r = (int)( normal(2.0*(double)NIND[gen]*Lambda_s, sqrt(2.0*(double)NIND[gen]*Lambda_s)) );
	return(r);
}

/* ***************************************************** */

int mutationnumberL()
{
	int r;
	if ((2.0*(double)NIND[gen]*Lambda_L < normalthreshold) && (exp(-2.0*(double)NIND[gen]*Lambda_L) != 0.0) )
	{
		r = poisson(lastinmutantspoissontableL, mutantspoissontableL);
	}
	else r = (int)( normal(2.0*(double)NIND[gen]*Lambda_L, sqrt(2.0*(double)NIND[gen]*Lambda_L)) );
	return(r);
}

/* ***************************************************** */

void disorder_NoSS (NoSS_k,NoSS_l)
int NoSS_k[], NoSS_l[];
{
	int a, b, rnd;
	
	for (i=0; i<countNoSS-1; i++)
	{
	   rnd=(int)(uniform()*(countNoSS-i));
	   a=NoSS_k[countNoSS-1-i]; NoSS_k[countNoSS-1-i]=NoSS_k[rnd]; NoSS_k[rnd]=a;
	   b=NoSS_l[countNoSS-1-i]; NoSS_l[countNoSS-1-i]=NoSS_l[rnd]; NoSS_l[rnd]=b;
	}
}

/* ***************************************************** */

dumpoffspringaftermutation()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\nOffspring after mutation (gm0 gm1)\n");	
	for (i=0; i<NIND[gen]; i++)   if (i<=10)    fprintf(fptr,"%d  %d\n",gm[i][0][0],gm[i][0][1]);
}

/* ***************************************************** */

neutral_genes ()
{
	H = 0.0;

	for (k=0; k<NCRO; k++)
	{
	    AA=0.0; Aa=0.0; aa=0.0;

	    for (i=0; i<NIND[gen]; i++)
	    {
	        if (((gm[i][k][0] & RM[0])==RM[0])&&((gm[i][k][1] & RM[0])==RM[0]))	aa+=1.0;
	        else    if (((gm[i][k][0] & RM[0])!=RM[0])&&((gm[i][k][1] & RM[0])!=RM[0])) AA+=1.0;
	        else		Aa+=1.0;
	    }

	    q[k][0] = (aa/(double)NIND[gen])+(Aa/(2.0*(double)NIND[gen]));
	    H += 2.0 * q[k][0] * (1.0 - q[k][0]);

//	    if (tracelevel!=0)    fprintf(fptr,"%d\t0\t%1.0f\t%1.0f\t%1.0f\t%f\t%d\n",k,AA,Aa,aa,q[k][0],initialgen[k][0]);	
	}
	accum(&Hw[gen], H/(double)NCRO);
}

/* ***************************************************** */

selected_genes()
{
	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NIND[gen]; i++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		q[k][l] = (aa/(double)NIND[gen])+(Aa/(2.0*(double)NIND[gen]));
		if ((gen == 0) && (q[k][l] == 0.0))    initialgen[k][l]=(-99);

		accum(&FREQ_HOMO[gen], aa/(double)NIND[gen]);

//COMMENT
//		if (tracelevel!=0)    fprintf(fptr,"\n k=%d l=%d q=%f  s=%f  h=%f", k, l, q[k][l], s[k][l], hs[k][l]);

//		if (q[k][l]==1.0)
//		{
//			/* make it wild type */
//			for (i=0; i<NIND[gen]; i++)
//			{
//				gm[i][k][0] = ( gm[i][k][0] & (~RM[l]) );
//				gm[i][k][1] = ( gm[i][k][1] & (~RM[l]) );
//			}
//			initialgen[k][l]=(-99);
//		if (fec_via == 0)
//		{
//			if (l%2==0)	addedfixed_sf *= (1.0 + s[k][l]);
//			else	addedfixed_sv *= (1.0 + s[k][l]);
//		}
//		else
//		{
//			if (l%3==0)	addedfixed_sf *= (1.0 + s[k][l]);
//			else	addedfixed_sv *= (1.0 + s[k][l]);
//		}
//			q[k][l]=0.0;
//			s[k][l]=0.0; hs[k][l]=0.0;
//		}
	} 
}

/* ***************************************************** */

genotypic_values ()
{
//COMMENT
	if (tracelevel!=0)	fprintf(fptr,"\n Genotypic values \n");

	for (i=0; i<NIND[gen]; i++)
	{
		leqLf[i] = 0.0;
		leqNLf[i] = 0.0;
		leqLv[i] = 0.0;
		leqNLv[i] = 0.0;
		leqT[i] = 0.0;

		genval_sf[i] = addedfixed_sf;
		genval_sv[i] = addedfixed_sv;

		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if (initialgen[k][l] != (-99))
		{
			if (tracelevel!=0)
			{
				if ( ((gm[i][k][0] & RM[l])==RM[l]) || ((gm[i][k][1] & RM[l])==RM[l]) )	if ((i<=5)&&(k<=2))    fprintf(fptr,"gen=%d i=%d k=%d l=%d s=%f h=%f q=%f ini=%d\n", gen, i, k, l, s[k][l], hs[k][l], q[k][l], initialgen[k][l]);
			}

			if (l%2==0)
			{
				// fecundity
		    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	genval_sf[i] *= (1.0 + s[k][l]);
				else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
				else	genval_sf[i] *= (1.0 + (s[k][l]*hs[k][l]));
			}
			else
			{
				// viability
		    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	genval_sv[i] *= (1.0 + s[k][l]);
				else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
				else	genval_sv[i] *= (1.0 + (s[k][l]*hs[k][l]));
			}

			if (s[k][l]==(-1.0))
			{
				leqT[i] += -(s[k][l] * q[k][l]);

				if (l%2==0)	leqLf[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
				else		leqLv[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
			}
			else
			{
				leqT[i] += -(s[k][l] * q[k][l]);

				if (l%2==0)	leqNLf[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
				else		leqNLv[i] += 2.0 * ((s[k][l]*hs[k][l]) - (s[k][l]/2.0)) * q[k][l] * (1.0 - q[k][l]);
			}
		}
//COMMENT
		if (tracelevel!=0)	fprintf(fptr," %d    genval_sf = %f    genval_sv = %f\n", i, genval_sf[i], genval_sv[i]);

		if (genval_sf[i]<0)	genval_sf[i]=0;
		if (genval_sv[i]<0)	genval_sv[i]=0;
	}
}

/* ***************************************************** */

genotypic_values_epistasis_ncro ()
{
	if (tracelevel!=0)	fprintf(fptr,"\n Values for epistasis \n");

	num_seg_loci=0;
	num_epis_loci=0;

	for (i=0; i<NIND[gen]; i++)
	{
		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)	homo[k][l] = 0;

		// Fecundity and viability

		for (k=0; k<NCRO; k++)
		{
			num_homo[k] = 0;

			for (l=1; l<NLOCI; l++)
			if (initialgen[k][l] != (-99))
			{
				num_seg_loci ++;
		    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l])) 
				{
					homo[k][l] = 1;
					num_homo[k] ++;

				}
			}
		}
		accum (&AVE_num_homo[gen+1], num_homo[k]);

		for (k=0; k<NCRO; k++)
		{
			if (num_homo[k] >= 2)
			for (l=1; l<NLOCI; l++)
			{
				if (homo[k][l] == 1)
				{
					// Fecundity
					if (episK==2)		genval_sf[i] *= pow( (1.0 + s[k][l]) , (double)episK-1.0);
					else if (episK==3)	genval_sf[i] *= pow( (1.0 + s[k][l]) , (double)num_homo[k]-1.0);

					// Viability
					if (episK==2)		genval_sv[i] *= pow( (1.0 + s[k][l]) , (double)episK-1.0);
					else if (episK==3)	genval_sv[i] *= pow( (1.0 + s[k][l]) , (double)num_homo[k]-1.0);

					num_epis_loci ++;
					if (tracelevel!=0)	if ((i<=10)&&(k==10)&&(l<=10))    fprintf(fptr,"fec i=%d k=%d l=%d s[k][l]=%f numepis=%d num_homo=%d\n",
								i, k, l, s[k][l], num_epis_loci, num_homo[k]);
				}
			}
		}

		accum(&num_EPIS[gen], (double)num_epis_loci);
		accum(&num_SEG[gen], (double)num_seg_loci);

		if (tracelevel!=0)	if (i<=10)    fprintf(fptr," %d    genval_sf = %f    genval_sv = %f\n", i, genval_sf[i], genval_sv[i]);

		if (genval_sf[i]<0)	genval_sf[i]=0;
		if (genval_sv[i]<0)	genval_sv[i]=0;
	}
}

/* ***************************************************** */

genotypic_values_epistasis_whole ()
{
//	if (tracelevel!=0)
//	{
//		fprintf(fptr,"\n\nGENERATION %d\n\n", gen);
//		for (i=0; i<NIND[gen]; i++)
//		{
//			fprintf(fptr,"IND %d  ", i);
//			for (k=0; k<NCRO; k++)
//			for (l=1; l<NLOCI; l++)
//			{
//				if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
//				{
//						fprintf(fptr,"%d_%d aa ", k, l);
//				}
//		    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /**/;
//				else
//				{
//					fprintf(fptr,"%d_%d Aa ",k, l);
//				}
//			}
//			fprintf(fptr,"\n");
//		}
//		
//	}

	if (tracelevel!=0)	fprintf(fptr,"\n Values for epistasis whole genome \n");

	num_seg_loci=0;
	num_epis_loci=0;

	for (i=0; i<NIND[gen]; i++)
	{
		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)	homo[k][l] = 0;

		// FECUNDITY and VIABILITY

		num_homo_whole = 0;

		for (k=0; k<NCRO; k++)
		{
			for (l=1; l<NLOCI; l++)
			{
				num_seg_loci ++;

		    		if (-s[k][l] > sEPISmin)
				if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l])) 
				{
					homo[k][l] = 1;
					num_homo_whole ++;
					if (tracelevel!=0)	fprintf(fptr,"epis i=%d k=%d l=%d -s=%f sEPISmin=%f ini=%d num_homo_whole=%d\n",
								i, k, l, -s[k][l], sEPISmin, initialgen[k][l], num_homo_whole);
				}
			}
		}
		accum (&AVE_num_homo_whole[gen+1], num_homo_whole);

		if (num_homo_whole >= 2)
		{
			// Fecundity
			for (k=0; k<NCRO; k++)
			for (l=2; l<(NLOCI-2); l+=2)
			{
				if (homo[k][l] == 1)
				{
					if (episK==2)		genval_sf[i] *= pow( (1.0 + s[k][l]) , (double)episK-1.0);
					else if (episK==3)	genval_sf[i] *= pow( (1.0 + s[k][l]) , (double)num_homo_whole-1.0);
					num_epis_loci ++;
				}
			}

			// Viability
			for (k=0; k<NCRO; k++)
			for (l=1; l<(NLOCI-2); l+=2)
			{
				if (homo[k][l] == 1)
				{
					if (episK==2)		genval_sv[i] *= pow( (1.0 + s[k][l]) , (double)episK-1.0);
					else if (episK==3)	genval_sv[i] *= pow( (1.0 + s[k][l]) , (double)num_homo_whole-1.0);
					num_epis_loci ++;
				}
			}
		}

		accum(&num_EPIS[gen], (double)num_epis_loci);
		accum(&num_SEG[gen], (double)num_seg_loci);

		if (genval_sf[i]<0)	genval_sf[i]=0;
		if (genval_sv[i]<0)	genval_sv[i]=0;

		if (tracelevel!=0)	fprintf(fptr," %d    genval_sf = %f    genval_sv = %f\n", i, genval_sf[i], genval_sv[i]);
	}
}

/* ***************************************************** */

distribution_qsh()
{
	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	if (initialgen[k][l] != (-99))
	{
		if (q[k][l] == 0.0)					accum (&q_ng_00[gen/tgen],  1.0);
		else if ( (q[k][l] > 0.0) && (q[k][l] < 0.1) )	accum (&q_ng_00_01[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.1) && (q[k][l] < 0.2) )	accum (&q_ng_01_02[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.2) && (q[k][l] < 0.3) )	accum (&q_ng_02_03[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.3) && (q[k][l] < 0.4) )	accum (&q_ng_03_04[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.4) && (q[k][l] < 0.5) )	accum (&q_ng_04_05[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.5) && (q[k][l] < 0.6) )	accum (&q_ng_05_06[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.6) && (q[k][l] < 0.7) )	accum (&q_ng_06_07[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.7) && (q[k][l] < 0.8) )	accum (&q_ng_07_08[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.8) && (q[k][l] < 0.9) )	accum (&q_ng_08_09[gen/tgen],  1.0);
		else if ( (q[k][l] >= 0.9) && (q[k][l] < 1.0) )	accum (&q_ng_09_10[gen/tgen],  1.0);
		else if (q[k][l] == 1.0)	accum (&q_ng_10[gen/tgen],  1.0);

		if (q[k][l] == 0.0)
		{
			if (s[k][l] == 0.0)							accum (&s_ng_0_lost[gen/tgen],  1.0);		
			else if ( (s[k][l] < 0.0) && (s[k][l] > (-0.000001)) )		accum (&s_ng_0_106_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.000001)) && (s[k][l] > (-0.0001)) )	accum (&s_ng_106_104_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.0001)) && (s[k][l] > (-0.01)) )	accum (&s_ng_104_102_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.01)) && (s[k][l] > (-0.1)) )		accum (&s_ng_102_101_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.1)) && (s[k][l] > (-0.2)) )		accum (&s_ng_01_02_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.2)) && (s[k][l] > (-0.4)) )		accum (&s_ng_02_04_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.4)) && (s[k][l] > (-0.6)) )		accum (&s_ng_04_06_lost[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.6)) && (s[k][l] > (-1.0)) )		accum (&s_ng_06_10_lost[gen/tgen],  1.0);
			else if (s[k][l] == (-1.0))						accum (&s_ng_10_lost[gen/tgen],  1.0);
		}
		else if (q[k][l] == 1.0)
		{
			if (s[k][l] == 0.0)							accum (&s_ng_0_fix[gen/tgen],  1.0);		
			else if ( (s[k][l] < 0.0) && (s[k][l] > (-0.000001)) )		accum (&s_ng_0_106_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.000001)) && (s[k][l] > (-0.0001)) )	accum (&s_ng_106_104_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.0001)) && (s[k][l] > (-0.01)) )	accum (&s_ng_104_102_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.01)) && (s[k][l] > (-0.1)) )		accum (&s_ng_102_101_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.1)) && (s[k][l] > (-0.2)) )		accum (&s_ng_01_02_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.2)) && (s[k][l] > (-0.4)) )		accum (&s_ng_02_04_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.4)) && (s[k][l] > (-0.6)) )		accum (&s_ng_04_06_fix[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.6)) && (s[k][l] > (-1.0)) )		accum (&s_ng_06_10_fix[gen/tgen],  1.0);
			else if (s[k][l] == (-1.0))						accum (&s_ng_10_fix[gen/tgen],  1.0);
		}
		else
		{
			if (s[k][l] == 0.0)		accum (&s_ng_0[gen/tgen],  1.0);		
			else if ( (s[k][l] < 0.0) && (s[k][l] > (-0.000001)) )		accum (&s_ng_0_106[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.000001)) && (s[k][l] > (-0.0001)) )	accum (&s_ng_106_104[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.0001)) && (s[k][l] > (-0.01)) )	accum (&s_ng_104_102[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.01)) && (s[k][l] > (-0.1)) )		accum (&s_ng_102_101[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.1)) && (s[k][l] > (-0.2)) )		accum (&s_ng_01_02[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.2)) && (s[k][l] > (-0.4)) )		accum (&s_ng_02_04[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.4)) && (s[k][l] > (-0.6)) )		accum (&s_ng_04_06[gen/tgen],  1.0);
			else if ( (s[k][l] <= (-0.6)) && (s[k][l] > (-1.0)) )		accum (&s_ng_06_10[gen/tgen],  1.0);
			else if (s[k][l] == (-1.0))		accum (&s_ng_10[gen/tgen],  1.0);
		}

		if (q[k][l] == 0.0)
		{
			if ( (hs[k][l] > 0.0) && (hs[k][l] <= 0.2) )		accum (&hs_ng_00_02_lost[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.2) && (hs[k][l] <= 0.4) )	accum (&hs_ng_02_04_lost[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.4) && (hs[k][l] <= 0.6) )	accum (&hs_ng_04_06_lost[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.6) && (hs[k][l] <= 0.8) )	accum (&hs_ng_06_08_lost[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.8) && (hs[k][l] <= 1.0) )	accum (&hs_ng_08_10_lost[gen/tgen],  1.0);
		}
		else if (q[k][l] == 1.0)
		{
			if ( (hs[k][l] > 0.0) && (hs[k][l] <= 0.2) )		accum (&hs_ng_00_02_fix[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.2) && (hs[k][l] <= 0.4) )	accum (&hs_ng_02_04_fix[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.4) && (hs[k][l] <= 0.6) )	accum (&hs_ng_04_06_fix[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.6) && (hs[k][l] <= 0.8) )	accum (&hs_ng_06_08_fix[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.8) && (hs[k][l] <= 1.0) )	accum (&hs_ng_08_10_fix[gen/tgen],  1.0);
		}
		else
		{
			if ( (hs[k][l] > 0.0) && (hs[k][l] <= 0.2) )		accum (&hs_ng_00_02[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.2) && (hs[k][l] <= 0.4) )	accum (&hs_ng_02_04[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.4) && (hs[k][l] <= 0.6) )	accum (&hs_ng_04_06[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.6) && (hs[k][l] <= 0.8) )	accum (&hs_ng_06_08[gen/tgen],  1.0);
			else if ( (hs[k][l] > 0.8) && (hs[k][l] <= 1.0) )	accum (&hs_ng_08_10[gen/tgen],  1.0);
		}
//COMMENT
//		if (tracelevel!=0)	fprintf(fptr,"k=%d l=%d s=%f hs=%f q=%f ini=%d\n", k, l, s[k][l], hs[k][l], q[k][l], initialgen[k][l]);
	}
}

/* **************************************************************** */

phenotypeB()
{
	int ii, it;
	double gsum_sf=0.0, gsum_sv=0.0, sum_leqNLf=0.0, sum_leqLf=0.0, sum_leqNLv=0.0, sum_leqLv=0.0, sum_leqT=0.0;

	// scaling to generation 0 fitness

	if (scaling == 1)	if (gen == 0)
	{
		Wfo=0.0; Wvo=0.0;
		for (i=0; i<NIND[gen]; i++)
		{
			Wfo += genval_sf[i]/(double)NIND[gen];
			Wvo += genval_sv[i]/(double)NIND[gen];
		}
	}

	if (tracelevel!=0)	fprintf(fptr,"\nWfo=%f Wvo=%f \n", Wfo, Wvo);

	for (i=0; i<NIND[gen]; i++)
	{
		if (scaling == 1)
		{
			pm_sf[i] = genval_sf[i] / Wfo;
			pm_sv[i] = genval_sv[i] / Wvo;
		}
		else
		{
			pm_sf[i] = genval_sf[i];
			pm_sv[i] = genval_sv[i];
		}

		sum_leqNLf += leqNLf[i];
		sum_leqLf += leqLf[i];
		sum_leqNLv += leqNLv[i];
		sum_leqLv += leqLv[i];
		sum_leqT += leqT[i];

		gsum_sf += pm_sf[i];
		gsum_sv += pm_sv[i];
	}

	accum (&LEQ_NLf[gen], sum_leqNLf/(double)NIND[gen]);
	accum (&LEQ_Lf[gen], sum_leqLf/(double)NIND[gen]);
	accum (&LEQ_NLv[gen], sum_leqNLv/(double)NIND[gen]);
	accum (&LEQ_Lv[gen], sum_leqLv/(double)NIND[gen]);
	accum (&LEQ_T[gen], sum_leqT/(double)NIND[gen]);

	accum (&gmean_sf[gen], gsum_sf/(double)NIND[gen]);
	accum (&gmean_sv[gen], gsum_sv/(double)NIND[gen]);

	sup[rep][gen] = 1.0;
	g_sf_sup[rep][gen] = gsum_sf/(double)NIND[gen];
	g_sv_sup[rep][gen] = gsum_sv/(double)NIND[gen];

//	printf("rep=%d gen=%d sup=%f gsf=%f gsv=%f\n", rep, gen, sup[rep][gen], g_sf_sup[rep][gen], g_sv_sup[rep][gen]);

	fprintf(fWrep,"%d %d %f\n", rep, gen, (gsum_sf/(double)NIND[gen])*(gsum_sv/(double)NIND[gen]));

	if (tracelevel!=0)   fprintf(fptr,"\nLEQ_NLf = %f  LEQ_Lf = %f  LEQ_NLv = %f  LEQ_Lv = %f  LEQ_T = %f  mean fec = %f  mean via = %f\n",
				sum_leqNLf/(double)NIND[gen], sum_leqLf/(double)NIND[gen], 
				sum_leqNLv/(double)NIND[gen], sum_leqLv/(double)NIND[gen], sum_leqT/(double)NIND[gen],
				gsum_sf/(double)NIND[gen], gsum_sv/(double)NIND[gen]);

	// Neutral case
	if (neutral == 1)
	for (i=0; i<NIND[gen]; i++)
	{
		pm_sf[i] = 1.0;
		pm_sv[i] = 1.0;
	}
}

/* ***************************************************** */

dumpphenotypes()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\nFecundity and viability of parents\n");
	for (i=0; i<NIND[gen]; i++)		fprintf(fptr,"%d  %f %f\n", i, pm_sf[i], pm_sv[i]);
}

/* ***************************************************** */

disorder_parents()
{
	for (i=0; i<NIND[gen]; i++)
	{
		ran_i=(int)(uniform()*NIND[gen]);
		for (k=0; k<NCRO; k++)
		{
			a=gm[i][k][0]; gm[i][k][0]=gm[ran_i][k][0]; gm[ran_i][k][0]=a;
			a=gm[i][k][1]; gm[i][k][1]=gm[ran_i][k][1]; gm[ran_i][k][1]=a;
		}
		w=pm_sf[i]; pm_sf[i]=pm_sf[ran_i]; pm_sf[ran_i]=w;
		w=pm_sv[i]; pm_sv[i]=pm_sv[ran_i]; pm_sv[ran_i]=w;

		a=father[i]; father[i]=father[ran_i]; father[ran_i]=a;
		a=mother[i]; mother[i]=mother[ran_i]; mother[ran_i]=a;
	}
}

/* ***************************************************** */

coancestry_matrix()
{
	if (gen==0)
	{
		for (i=0; i<NIND[gen]; i++)
		for (j=0; j<NIND[gen]; j++)
		{
			if (i == j)	mat[i][j] = 0.5;
			else		mat[i][j] = 0.0;
		}
	}
	else
	{
		for (i=0; i<NIND[gen]; i++)
		for (j=0; j<NIND[gen]; j++)
		{
			if (i == j)	mat[i][j] = 0.5 * (1.0 + smat[father[i]][mother[i]]);
			else		mat[i][j] = 0.25 * (smat[father[i]][father[j]] + smat[father[i]][mother[j]] + smat[mother[i]][father[j]] + smat[mother[i]][mother[j]]);
		}
	}

	for (i=0; i<NIND[gen]; i++)
	for (j=0; j<NIND[gen]; j++)
	{
		smat[i][j] = mat[i][j];
		accum(&fp[gen], mat[i][j]); 
		accum(&Fp[gen], (2.0*mat[i][i]-1.0)); 
	}

//	if (tracelevel!=0)    
//	for (i=0; i<NIND[gen]; i++)
//	{
//		fprintf(fptr, "gen=%d   i=%d   fat=%d   mot=%d   Fp=%f\n", gen, i, father[i], mother[i], 2.0*mat[i][i]-1.0);
//	}

	/* **************** PEDIGREE FILE ***************** */

//	if (gen==0)	fprintf(fped, "ind   fath   moth   fit   Fp   gen\n");
//	for (i=0; i<NIND[gen]; i++)
//	{
//		n ++;
//		if (gen == 0)	fprintf(fped, "%d   %5d   %5d   %f   %f   %d\n", n-1, 0, 0, pm_sf[i]*pm_sf[i], 2.0*mat[i][i]-1.0, gen);
//		else if (gen == 1)	fprintf(fped, "%d   %5d   %5d   %f   %f   %d\n", n-1, father[i], mother[i], pm_sf[i]*pm_sf[i], 2.0*mat[i][i]-1.0, gen);
//		else			fprintf(fped, "%d   %5d   %5d   %f   %f   %d\n", n-1, father[i]+nn[gen-2], mother[i]+nn[gen-2], pm_sf[i]*pm_sf[i], 2.0*mat[i][i]-1.0, gen);
//		nn[gen] = n;
//	}
}

/* ***************************************************** */

dumpparents()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\nParents\n");
	for (i=0; i<NIND[gen]; i++)   fprintf(fptr,"i=%d  father=%d  mother=%d  pm_sf=%f pm_sv=%f  Fp=%f\n", i, father[i], mother[i], pm_sf[i], pm_sv[i], 2.0*mat[i][i]-1.0);
}

/* ***************************************************** */

fecundities()
{
	double cum=0.0;

	/********* keep genotypes of parents ***********/ 

	for (i=0; i<NIND[gen]; i++)
	for (k=0; k<NCRO; k++)
	{
		sm[i][k][0]=gm[i][k][0];
		sm[i][k][1]=gm[i][k][1];
	}

	/********* Accumulate psf fecundity values ***********/

	for (i=0; i<NIND[gen]; i++)
	{
		cum += pm_sf[i];
		cum_pm_sf[i] = cum;
	}
	for (i=0; i<NIND[gen]; i++)	cum_pm_sf[i] = cum_pm_sf[i] / cum;

//	COMMENT
	if (tracelevel!=0)	for (i=0; i<NIND[gen]; i++)	fprintf (fptr," %d    cum_pm_sf = %f\n", i, cum_pm_sf[i]);
}

/* ***************************************************** */

int fecunditysearch(double *array, int size)
{
	int i, p;
	double r;

	r = uniform();
	if ( (r >= 0.0) && (r <= array[0]) )	p = 0; 
	else for (i=1; i<size; i++)
		if ( (r > array[i-1]) && (r <= array[i]) )	p = i; 

	return(p);
}

/* ***************************************************** */

mating_RC()
{
	int rnd, a;
	double K;
	for (i=0; i<NIND[gen]; i++)	family[i] = 0;
	for (i=0; i<NMAX; i++)	failed[i] = 1;
	
	/******************************************/

	if (tracelevel!=0)	fprintf(fptr,"\nMatings MAXFEC %f\n", MAXFEC);

	/******************* Round of random contributions mating ***********************/

	if (tracelevel!=0)	fprintf (fptr,"\nRound of random contributions mating\n\n");

	born = 0;
	dead = 0;

	if (NIND[gen]%2 != 0)	NIND[gen] = NIND[gen] - 1;

	for (i=0; i<NIND[gen]; i+=2)
	for (m=1; m<=2; m++) 
	{
		p1 = i;
		p2 = p1 + 1;

		K = normal(MAXFEC, DS_MAXFEC) * sqrt(pm_sf[p1] * pm_sf[p2]);
		maxfec[i] = (int)(K);
//		maxfec[i] = (int)(K + 0.5);
//		if (tracelevel!=0)	fprintf(fptr,"\nK  %f    macfec %d\n", K, maxfec[i]);

		for (j=0; j<maxfec[i]; j++)
		{
			if (tracelevel!=0)	fprintf (fptr,"Progeny %d - %d  maxfec %d  p1 = %d    p2 = %d\n", i, j, maxfec[i], p1, p2);
			gametes_offspring();
			bornoffspring();
		}
	}

	// First assignation of progeny to population

	if (tracelevel!=0)	fprintf (fptr,"First assignation of progeny to population\n\n");

	for (b=1; b<=born; b++)	mark[b] = 0;
	assign = 0;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	fprintf (fptr,"Available born before first assignment b %d  fa %d  mo %d  mark %d\n", b, fathero[b], mothero[b], mark[b]);
	if (tracelevel!=0)	fprintf (fptr,"\n");

	// Randomise born progeny

	for (b=1; b<=born; b++)
	{
		rnd = (int)(uniform()*born) + 1;

		for (k=0; k<NCRO; k++)
		{
		   	a = gmo[b][k][0]; gmo[b][k][0] = gmo[rnd][k][0]; gmo[rnd][k][0] = a;
		   	a = gmo[b][k][1]; gmo[b][k][1] = gmo[rnd][k][1]; gmo[rnd][k][1] = a;
		}
	   	a = fathero[b]; fathero[b] = fathero[rnd]; fathero[rnd] = a;
	   	a = mothero[b]; mothero[b] = mothero[rnd]; mothero[rnd] = a;
	   	a = tagPSo[b]; tagPSo[b] = tagPSo[rnd]; tagPSo[rnd] = a;
	   	a = tagFSo[b]; tagFSo[b] = tagFSo[rnd]; tagFSo[rnd] = a;
	   	a = tagHSo[b]; tagHSo[b] = tagHSo[rnd]; tagHSo[rnd] = a;
	   	a = mark[b]; mark[b] = mark[rnd]; mark[rnd] = a;
	}
	if (tracelevel!=0)	for (b=1; b<=born; b++)
	fprintf (fptr,"Randomised born b %d  mark_b %d  fa %d  mo %d\n", b, mark[b], fathero[b], mothero[b]);

	// Assignation

	if (assign == NMAX) goto endofmating_RC;

	if (tracelevel!=0)	for (b=1; b<=born; b++)	if (mark[b] == 0)
	fprintf (fptr,"\nAvailable born before final assignment b %d  fa %d  mo %d  mark %d", b, fathero[b], mothero[b], mark[b]);
	if (tracelevel!=0)	fprintf (fptr,"\n\n");

	for (i=0; i<NMAX; i++)  if (failed[i] == 1)
	{
		if (assign == NMAX) goto endofmating_RC;

		for (b=1; b<=born; b++)  if (mark[b] == 0)
		{
			assign ++;
			for (k=0; k<NCRO; k++)
			{
			   	gm[i][k][0] = gmo[b][k][0];
			   	gm[i][k][1] = gmo[b][k][1];
			}
			father[i] = fathero[b];
			mother[i] = mothero[b];
			tagPS[i] = tagPSo[b];
			tagFS[i] = tagFSo[b];
			tagHS[i] = tagHSo[b];
			family[father[i]] += 1;
			family[mother[i]] += 1;
			mark[b] = 1;
			failed[i] = 0;
			if (tracelevel!=0)	fprintf (fptr,"After final assignment i %d  b %d  mark_b %d  fa %d  mo %d failed_i %d\n", i, b, mark[b], fathero[b], mothero[b], failed[i]);
			break;
		}
	}

	/******************* End of rounds of mating ***********************/

	endofmating_RC: /******/;

	NIND[gen+1] = assign;
}

/* ***************************************************** */

printtriedoffspring()
{
	if (tracelevel!=0)
	{
		fprintf (fptr,"try %d  Tried parents  p1 = %d    p2 = %d   (p1f=%f p1v=%f p2f=%f p2v=%f)\n",
						j, p1, p2, pm_sf[p1], pm_sv[p1], pm_sf[p2], pm_sv[p2]);
	}
}

/* ***************************************************** */

gametes_offspring()
{
	int EE[MM], FF[MM], marker;
	int numberrecs, nr, pointrec[MM][31], ncrorec[MM], rndk, rndl;
	double rnd;

	/* ******************* Gametes with free recombination ******************* */

	    if(L==99.0)
	    {	    /* ******************* Free recombination ******************* */
		for (k=0; k<NCRO; k++)
		{
		   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
			FF[k] = ~EE[k];
		   	gmb[j][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
		}
// COMMENT
//		if (tracelevel!=0)   fprintf (fptr,"j=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n",
//		j, EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1],
//		sm[p1][2][0], sm[p1][2][1], gmb[j][0][0], gmb[j][1][0], gmb[j][2][0]);

		for (k=0; k<NCRO; k++)
		{
		   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
		   	FF[k] = ~EE[k];
		   	gmb[j][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
		}
//		if (tracelevel!=0)   fprintf (fptr,"j=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", j,
//		EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 
//		sm[p1][2][0], sm[p1][2][1], gmb[j][0][0], gmb[j][1][0], gmb[j][2][0]);
	   }

	    else
	    {	    /* ************** Restricted recombination ***************** */

			/* ****** Chromosome from father ****** */

			numberrecs = recombinationnumber();

// ************ NO RECOMBINATION IN MALES ******************************

			numberrecs = 0;

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}
			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      		{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)
				{
					EE[k] = ~EE[k];
				}
				FF[k] = ~EE[k];
				gmb[j][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
			}

// COMMENT
/*			if (tracelevel!=0)   fprintf (fptr,"j=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n",
			j, EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1],
			sm[p1][2][0], sm[p1][2][1], gmb[j][0][0], gmb[j][1][0], gmb[j][2][0]);

			if (tracelevel!=0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
			}
*/
			/* ****** Chromosome from mother ****** */

			numberrecs = recombinationnumber();
			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}
			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
			}

			marker = 1;


			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      		{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)
				{
					EE[k] = ~EE[k];
				}
				FF[k] = ~EE[k];
				gmb[j][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
			}
// COMMENT
/*			if (tracelevel!=0)   fprintf (fptr,"j=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", j,
			EE[0], EE[1], EE[2], sm[p2][0][0], sm[p2][0][1], sm[p2][1][0], sm[p2][1][1], 
			sm[p2][2][0], sm[p2][2][1], gmb[j][0][1], gmb[j][1][1], gmb[j][2][1]);

			if (tracelevel!=0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<50; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
			}
*/
	    }

	/* ****************** Genotypic value of offspring ****************** */

	pm_sfb[j] = addedfixed_sf;
	pm_svb[j] = addedfixed_sv;

	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	{
		// fecundity
		if ( (l%2 == 0) && (initialgen[k][l] != (-99)) ) 
		{
			if (((gmb[j][k][0] & RM[l])==RM[l])&&((gmb[j][k][1] & RM[l])==RM[l]))  	pm_sfb[j] *= (1.0 + s[k][l]);
			else    if (((gmb[j][k][0] & RM[l])!=RM[l])&&((gmb[j][k][1] & RM[l])!=RM[l])) 		/* AA */;
			else	pm_sfb[j] *= (1.0 + (s[k][l]*hs[k][l]));
		}
		// viability
		if ( (l%2 != 0) && (initialgen[k][l] != (-99)) ) 
		if (initialgen[k][l] != (-99))
		{
			if (((gmb[j][k][0] & RM[l])==RM[l])&&((gmb[j][k][1] & RM[l])==RM[l]))  	pm_svb[j] *= (1.0 + s[k][l]);
			else    if (((gmb[j][k][0] & RM[l])!=RM[l])&&((gmb[j][k][1] & RM[l])!=RM[l])) 		/* AA */;
			else	pm_svb[j] *= (1.0 + (s[k][l]*hs[k][l]));
		}
	}

	// scaling by generation 0 fitness

	if (scaling == 1)
	{
		pm_sfb[j] = pm_sfb[j] / Wfo;
		pm_svb[j] = pm_svb[j] / Wvo;
	}
	else
	{
		pm_sfb[j] = pm_sfb[j];
		pm_svb[j] = pm_svb[j];
	}

	if (tracelevel!=0)	fprintf (fptr,"Offspring   fec = %f   via = %f\n", pm_sfb[j], pm_svb[j]);

	if (tracelevel!=0)	fprintf (fptr,"Offspring   born = %d   dead = %d\n", born, dead);

	if ((born+dead) > 0)  accum (&gmean_sf_zygotes[gen+1], pm_sfb[j]);
	if ((born+dead) > 0)  accum (&gmean_sv_zygotes[gen+1], pm_svb[j]);

	if (tracelevel!=0)	fprintf (fptr,"gen = %d   mean_pm_sfb = %f   mean_pm_svb = %f\n", gen, accmean(&gmean_sf_zygotes[gen+1]), accmean(&gmean_sv_zygotes[gen+1]));

	if (neutral == 1) pm_svb[j] = 1.0;
}

/* ***************************************************** */

int recombinationnumber ()
{
	int r;
	if ((L < normalthreshold) && (exp(-L) != 0.0) )
	{
		r = poisson(lastinrecombinantpoissontable, recombinantpoissontable);
	}
	else r = (int)normal(L, sqrt(L));
	return(r);
}

/* ***************************************************** */

bornoffspring()
{
	/********* Born or dead individuals *********/

	if (uniform() <= pm_svb[j])
	{
		born ++;

		for (k=0; k<NCRO; k++)
		{
		   	gmo[born][k][0] = gmb[j][k][0];
		   	gmo[born][k][1] = gmb[j][k][1];
		}
		fathero[born] = p1;
		mothero[born] = p2;
		if (p1 == p2)	tagPSo[born] = 1;
		if (gen > 0)
		{	
			if ( ((father[p1] == father[p2]) && (mother[p1] == mother[p2])) ||
				((father[p1] == mother[p2]) && (mother[p1] == father[p2])) )	tagFSo[born] = 1;
			else if ((father[p1] == father[p2]) || (father[p1] == mother[p2]) ||
				(mother[p1] == father[p2]) || (mother[p1] == mother[p2]))	tagHSo[born] = 1;
		}
		if (tracelevel!=0)
		{
			fprintf (fptr,"tagPSo=%d tagFSo=%d tagHSo=%d\n", tagPSo[born], tagFSo[born], tagHSo[born]);
			fprintf (fptr,"Born %d\n\n", born);
		}
	}
	else
	{
		if (tracelevel!=0)	fprintf (fptr,"Dead\n\n");
	}
}

/* ***************************************************** */

calculations()
{
	double family_sum=0.0, family_sum2=0.0;

	numPSborn=0.0;
	numFSborn=0.0;
	numHSborn=0.0;

	for (i=0; i<assign; i++)
	{
		numPSborn += (double)tagPS[i] / assign;
		numFSborn += (double)tagFS[i] / assign;
		numHSborn += (double)tagHS[i] / assign;
	}

	/************ P_SELF *************/

	if (assign > 0)		accum(&P_PSBORN[gen], numPSborn);
	if (assign > 0)		accum(&P_FSBORN[gen], numFSborn);
	if (assign > 0)		accum(&P_HSBORN[gen], numHSborn);

	/************ SK2, Ne(Sk2) and FVk *************/

	for (i=0; i<assign; i++)
	{
		family_sum += family[i];
		family_sum2 += (family[i] * family[i]);
	}
	sk2 = ((family_sum2-(family_sum*family_sum/assign))/assign);
	if (assign != 0)	accum(&SK2[gen], sk2);

	ne = (4.0 * assign) / (2.0 + sk2);
	if (assign != 0)	accum(&NeVk[gen], ne);

	fvk = 1.0 - pow(1.0 - (1/(2.0*ne)), gen);
	if (assign != 0)	accum(&Fvk[gen], fvk);

	if (tracelevel!=0)  fprintf (fptr,"\nContributions from parents\n");
	if (tracelevel!=0)  for (i=0; i<NIND[gen]; i++)  fprintf (fptr,"family(%d) = %d\n", i, family[i]);
	if (tracelevel!=0)	fprintf (fptr,"SK2 = %5.3f\n", sk2);

	for (i=0; i<assign; i++) { tagPS[i]=0; tagFS[i]=0; tagHS[i]=0; }
	for (i=1; i<=born; i++) { tagPSo[i]=0; tagFSo[i]=0; tagHSo[i]=0; }

	/************ MEAN PRODUCTIVITY *************/

	accum (&PROD[gen], born);

}

/* ***************************************************** */

dumpoffspring()
{
	if (tracelevel==0)   return (0);

	fprintf(fptr,"\n Offspring before mutation (gm0 gm1)\n");	
	for (i=0; i<NIND[gen]; i++)   fprintf(fptr,"%d  %d\n",gm[i][0][0],gm[i][0][1]);
}

/* ***************************************************** */

printout()
{
	fprintf(fgen,"\n\nNMAX=%d   L=%f  PS=%4.2f	N.S-LOCI=%d    N.N-LOCI=%d    reps=%d   gens=%d\n    MAXFEC=%f    DS_MAXFEC=%f    Lambda_s=%f   Lambda_L=%f\n    ave_s=%f   relax=%f    neutral=%d   scaling=%d   episK=%d   episType=%d\n    beta_s=%f   alpha_s=%f  dom=%d   k_s=%f   ave_hs=%f    addedfixed_sNP=%f\n",
		NMAX, L, PS, TOTLOCI-NCRO, NCRO, replicates, generations, MAXFEC, DS_MAXFEC, Lambda_s, Lambda_L, ave_s, relaxedsel, neutral, scaling, episK, episType, beta_s, alpha_s, dom, k_s, ave_hs, addedfixed_sNP);

	fprintf(fgen,"LEQ(2dpq) = %f    LEQL = %f    d = %f    LEQf = %f    LEQv = %f    TLEQ(sq) = %f\n", LEQ_NP, LEQ_NPL, D_NP/NSEGLOCNP, LEQ_NPf, LEQ_NPv, TLEQ_NP);
	fprintf(fgen,"TIME TO EXTINCTION = %f   var = %f    NSEGLOCNP = %d    num_epis_loci = %d\n", accmean(&TIME_EXT), variance(&TIME_EXT), NSEGLOCNP, num_epis_loci);
	
	fprintf(fgen,"\ngen   PROD      Wfzyg      Wvzyg     N          fp        Fp        Fvk       NeVk      Wf        Wv     Wf(sup)  Wv(sup)   LEQ_NLf   LEQ_Lf    LEQ_NLv   LEQ_Lv    LEQ_T     SK2   FSborn    HSborn    Hw    EPIS    num_loci    F    W    lnW    W(sup)    min   max   lnW(sup)   min   max   FREQ_HOMO\n");

	min=(int)(nsup*0.025);
	if (min==0)	min=1;
	max=(int)(nsup*0.975);

	for (gen=0; gen<=generations; gen++)
	{
		if (accmean(&mean_N[gen]) > 0)
		{
			fprintf(fgen, "%d   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",
			gen, accmean(&PROD[gen]), accmean(&gmean_sf_zygotes[gen]), accmean(&gmean_sv_zygotes[gen]), accmean(&mean_N[gen]), accmean(&fp[gen]), accmean(&Fp[gen]), accmean(&Fvk[gen]), accmean(&NeVk[gen]), accmean(&gmean_sf[gen]), accmean(&gmean_sv[gen]), 
			accmean(&gmean_sf_sup[gen]), accmean(&gmean_sv_sup[gen]), accmean(&LEQ_NLf[gen]),
			accmean(&LEQ_Lf[gen]), accmean(&LEQ_NLv[gen]), accmean(&LEQ_Lv[gen]), accmean(&LEQ_T[gen]), accmean(&SK2[gen]),
			accmean(&P_FSBORN[gen]), accmean(&P_HSBORN[gen]), accmean(&Hw[gen]), accmean(&num_EPIS[gen]), accmean(&num_SEG[gen]),
			(accmean(&Fp[gen-1])+accmean(&Fp[gen]))/2.0, accmean(&gmean_sf[gen-1])*accmean(&gmean_sv[gen]), log(accmean(&gmean_sf[gen-1])*accmean(&gmean_sv[gen])),
			accmean(&gmean_sf_sup[gen-1])*accmean(&gmean_sv_sup[gen]), W[min][gen], W[max][gen],
			log(accmean(&gmean_sf_sup[gen-1])*accmean(&gmean_sv_sup[gen])), log(W[max][gen]), log(W[min][gen]), accmean(&FREQ_HOMO[gen]));
		}
		else
			fprintf(fgen, "%d\n", gen);
	}

//	SALIDA GLOBAL

	fprintf(fSAL,"\ngen PROD  Wzy  logWzy  N    LOCI  K   L     s     h     B     e  t  snp    LEQnp LEQnl  Hw    PrE     LEQnl LEQl   F      W       lnW       W(sup)  min   max   lnW(sup)    min       max       PrSurv   homo  homoW  FREQ_HOMO\n");

	for (gen=0; gen<=generations; gen++)
	{
		if (accmean(&mean_N[gen]) > 0)
		{
		   if ( (gen==0) || (gen==1) )
		   {
			fprintf(fSAL, "%d  %4.2f  1.00  0.00  %4.2f  %4.0d  %d  %4.2f  %4.2f  %4.2f  %4.2f  %d  %d  %d   %4.2f  %4.2f  %5.3f  %6.3f  %4.2f  %4.2f  %5.3f  %5.3f  %9.6f  %5.3f  %5.3f  %5.3f  %9.6f  %9.6f  %9.6f  %4.2f  %4.2f  %4.2f  %f\n",
			gen, accmean(&PROD[gen]), accmean(&mean_N[gen]), TOTLOCI, (int)MAXFEC, Lambda_s, ave_s, ave_hs, beta_s, episK, episType, NSEGLOCNP, LEQ_NP, LEQ_NPL, accmean(&Hw[gen]), accmean(&num_EPIS[gen])*100/accmean(&num_SEG[gen]),
			accmean(&LEQ_NLf[gen])+accmean(&LEQ_NLv[gen]), accmean(&LEQ_Lf[gen])+accmean(&LEQ_Lv[gen]),
			(accmean(&Fp[gen-1])+accmean(&Fp[gen]))/2.0, accmean(&gmean_sf[gen-1])*accmean(&gmean_sv[gen]), log(accmean(&gmean_sf[gen-1])*accmean(&gmean_sv[gen])),
			accmean(&gmean_sf_sup[gen-1])*accmean(&gmean_sv_sup[gen]), W[min][gen], W[max][gen],
			log(accmean(&gmean_sf_sup[gen-1])*accmean(&gmean_sv_sup[gen])), log(W[max][gen]), log(W[min][gen]), reps_alive_gen[gen]*100.0/replicates, accmean(&AVE_num_homo[gen]), accmean(&AVE_num_homo_whole[gen]), accmean(&FREQ_HOMO[gen]) );
		   }
		   else
		   {
			fprintf(fSAL, "%d  %4.2f  %4.2f  %4.2f  %4.2f  %4.0d  %d  %4.2f  %4.2f  %4.2f  %4.2f  %d  %d  %d   %4.2f  %4.2f  %5.3f  %6.3f  %4.2f  %4.2f  %5.3f  %5.3f  %9.6f  %5.3f  %5.3f  %5.3f  %9.6f  %9.6f  %9.6f  %4.2f  %4.2f  %4.2f  %f\n",
			gen, accmean(&PROD[gen]), accmean(&gmean_sf_zygotes[gen-1])*accmean(&gmean_sv_zygotes[gen]), log(accmean(&gmean_sf_zygotes[gen-1])*accmean(&gmean_sv_zygotes[gen])), accmean(&mean_N[gen]), TOTLOCI, (int)MAXFEC, Lambda_s, ave_s, ave_hs, beta_s, episK, episType, NSEGLOCNP, LEQ_NP, LEQ_NPL, accmean(&Hw[gen]), accmean(&num_EPIS[gen])*100/accmean(&num_SEG[gen]),
			accmean(&LEQ_NLf[gen])+accmean(&LEQ_NLv[gen]), accmean(&LEQ_Lf[gen])+accmean(&LEQ_Lv[gen]),
			(accmean(&Fp[gen-1])+accmean(&Fp[gen]))/2.0, accmean(&gmean_sf[gen-1])*accmean(&gmean_sv[gen]), log(accmean(&gmean_sf[gen-1])*accmean(&gmean_sv[gen])),
			accmean(&gmean_sf_sup[gen-1])*accmean(&gmean_sv_sup[gen]), W[min][gen], W[max][gen],
			log(accmean(&gmean_sf_sup[gen-1])*accmean(&gmean_sv_sup[gen])), log(W[max][gen]), log(W[min][gen]), reps_alive_gen[gen]*100.0/replicates, accmean(&AVE_num_homo[gen]), accmean(&AVE_num_homo_whole[gen]) , accmean(&FREQ_HOMO[gen]));
		   }
		}
		else
			fprintf(fSAL, "%d\n", gen);
	}


	//	outfile summary  Ls - Ll - s - hs - NMAX - MAXFEC - PS - gen - TimeEXT - N - Mf - Mv - LEQ_NLf - LEQ_Lf - LEQ_NLv - LE_Lv - LE_T - SK2 - selfB - Hw  

	// gen 0
//	fprintf(fout, "%f  %f  %f  %f  %d  %d  %f  0  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", Lambda_s, Lambda_L, ave_s, ave_hs, NMAX, MAXFEC, PS,
//		accmean(&TIME_EXT), accmean(&mean_N[0]), accmean(&gmean_sf[0]), accmean(&gmean_sv[0]), accmean(&LEQ_NLf[0]),
//		accmean(&LEQ_Lf[0]), accmean(&LEQ_NLv[0]), accmean(&LEQ_Lv[0]), accmean(&LEQ_T[0]), accmean(&SK2[0]), 
//		accmean(&P_FSBORN[0]), accmean(&P_HSBORN[0]), accmean(&Hw[0]));
	// gen 10
//	fprintf(fout, "%f  %f  %f  %f  %d  %d  %f  10  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", Lambda_s, Lambda_L, ave_s, ave_hs, NMAX, MAXFEC, PS,
//		accmean(&TIME_EXT), accmean(&mean_N[10]), accmean(&gmean_sf[10]), accmean(&gmean_sv[10]), accmean(&LEQ_NLf[10]),
//		accmean(&LEQ_Lf[10]), accmean(&LEQ_NLv[10]), accmean(&LEQ_Lv[10]), accmean(&LEQ_T[10]), accmean(&SK2[10]),
//		accmean(&P_FSBORN[10]), accmean(&P_HSBORN[10]), accmean(&Hw[10]));
	// gen generations
//	fprintf(fout, "%f  %f  %f  %f  %d  %d  %f  %d  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", Lambda_s, Lambda_L, ave_s, ave_hs, NMAX, MAXFEC, PS,
//		generations, accmean(&TIME_EXT), accmean(&mean_N[generations]), accmean(&gmean_sf[generations]), accmean(&gmean_sv[generations]), accmean(&LEQ_NLf[generations]),
//		accmean(&LEQ_Lf[generations]), accmean(&LEQ_NLv[generations]), accmean(&LEQ_T[generations]), accmean(&LEQ_Lv[generations]),
//		accmean(&P_FSBORN[generations]), accmean(&P_HSBORN[generations]), accmean(&Hw[generations]));
//	fprintf(fout, "\n\n\n");
}

/* ***************************************************** */

distribution_out()
{
	fprintf (fgen, "\n");
	for (i=0; i<=(generations/tgen); i++) 	fprintf (fgen, "Extinction risk gen %d = %f\n", i*tgen, (replicates-reps_alive_gen[i])*100.0/replicates);

	fprintf (fgen, "\nGen		");
	for (i=0; i<=(generations/tgen); i++) 	fprintf (fgen, "%d     ", i*tgen);
	fprintf (fgen, "\n");

	// distribution of q values in the population

		fprintf (fgen, "q=0.0          ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_00[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.0-0.1      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_00_01[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.1-0.2      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_01_02[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.2-0.3      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_02_03[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.3-0.4      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_03_04[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.4-0.5      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_04_05[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.5-0.6      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_05_06[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.6-0.7      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_06_07[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.7-0.8      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_07_08[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.8-0.9      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_08_09[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=0.9-1.0      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_09_10[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "q=1.0          ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&q_ng_10[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

	// distribution of s values in the population

		fprintf (fgen, "LOST\n");

		fprintf (fgen, "s=0.0         ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=0.0-10-6    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0_106_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-6-10-4  ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_106_104_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-4-0.01 ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_104_102_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.01-0.1  ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_102_101_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.1-0.2    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_01_02_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.2-0.4    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_02_04_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.4-0.6    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_04_06_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.6-1.0    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_06_10_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-1.0        ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_10_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "SEGREGATING\n");

		fprintf (fgen, "s=0.0         ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=0.0-10-6    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0_106[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-6-10-4  ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_106_104[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-4-0.001 ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_104_102[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.001-0.1  ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_102_101[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.1-0.2    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_01_02[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.2-0.4    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_02_04[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.4-0.6    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_04_06[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.6-1.0    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_06_10[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-1.0        ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_10[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "FIXED\n");

		fprintf (fgen, "s=0.0         ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=0.0-10-6    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_0_106_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-6-10-4  ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_106_104_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-10-4-0.001 ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_104_102_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.001-0.1  ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_102_101_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.1-0.2    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_01_02_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.2-0.4    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_02_04_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.4-0.6    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_04_06_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-0.6-1.0    ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_06_10_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "s=-1.0        ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&s_ng_10_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n\n");

	// distribution of hs values in the population

		fprintf (fgen, "LOST\n");

		fprintf (fgen, "hs=0.0-0.2      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_00_02_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.2-0.4      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_02_04_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.4-0.6      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_04_06_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.6-0.8      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_06_08_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.8-1.0      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_08_10_lost[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "SEGREGATING\n");

		fprintf (fgen, "hs=0.0-0.2      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_00_02[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.2-0.4      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_02_04[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.4-0.6      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_04_06[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.6-0.8      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_06_08[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.8-1.0      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_08_10[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "FIXED\n");

		fprintf (fgen, "hs=0.0-0.2      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_00_02_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.2-0.4      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_02_04_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.4-0.6      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_04_06_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.6-0.8      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_06_08_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");

		fprintf (fgen, "hs=0.8-1.0      ");
		for (i=0; i<=10; i++) if (NIND[i*tgen] != 0) fprintf (fgen, "%8.2f  ", accsum(&hs_ng_08_10_fix[i])/reps_alive_gen[i]);
		fprintf (fgen, "\n");
}

/* ***************************************************** */



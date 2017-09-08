/* */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define NSTACK 50
#define pi 3.1416
#define round 0.01
#define ndd 20

/******************* global variables *************************************/
 unsigned int	i ,j , k, kk, ii, jj, nt, np, jt, nav, *npix, *npiy,  nhist, clmi,
		 					i0, i0max, ntot, navi, clmerr, ihisterr, *nltoti, i0i, i0ii,
							nrtt, *nprt, np1, n00, wl, sl;
 int ihist, icl, ihisto, i00, nth, nthmin, nthmax, nc, nthp, ncp, ihistdd;
 float		dum, x, y, *kl, *knt;
 double   xl0, yl0, *xi, *yi, clm, amin, amax, aav, xm, ym, xmi, ymi,
  				clij, *ac, *acq, *acm, *dd,  *actr, *acqtr, *acmtr, *ddtr, ddi, gewm,
			 			*histc, *histdd, *nhistc, *nhistdd, cltot, tres, xli, yli, cosij,
						*noro, aavmin, aavmax, ddmax, clmm, *cltoti, *msdcltr, *msdcl,
						*msdcl2tr, *msdcl2, *histddtr, *nhistddtr, *nnormt,
						cltotj, ntotj, clmmj, rtt, *xrt, *yrt, *lstep, dt, xxm,
						cstep, astep, *histctr, *nhistctr, clmp, aminp, amaxp, xx, *xiyi,
					*ddttr, *rmsdclttr, *msdclttr, *ddt, *rmsdclt, *msdclt, rmsdcl;

 char		buf[70], bufe[70], trackchar[4];
 char *bufein, *fstnumb;
 char *del=" ,;:";
 FILE*		eing;
 FILE*		ausg;

/***************************************************************************/

/******************* global functions ************************************/
void ainit1();
void ainit2();
void pardet();

/***************************************************************************/

int main(int argc, char *argv[])
{
/*-------------------------------------------------------------------------*/
 //unsigned long ill;

  if (argc != 27 )
  {
    fprintf(stderr,"\n\n  command line parameters wrong , correct : \n\n"
			"./persist_2d_26  input_file   LaminA_111109_510.dat "
			" nt 510  np 37  used_points 37 nav 4 start_position 0 "
			" cstep 10 astep 10  clmp 100 aminp -20 amaxp 80 cutoff 0.1 dt 10 \n\n\n");
    exit(2);
  }

  nt  = ((unsigned int) atol(argv[4]));
  np  = ((unsigned int) atol(argv[6]));
	np1  = ((unsigned int) atol(argv[8]));
	nav = ((unsigned int) atol(argv[10]));
	i00 = ((unsigned int) atol(argv[12]));
  cstep = ((double) atof(argv[14]));
	astep = ((double) atof(argv[16]));
  clmp  = ((double) atof(argv[18]));
  aminp = ((double) atof(argv[20]));
	amaxp = ((double) atof(argv[22]));
	tres = ((double) atof(argv[24]));
	dt = ((double) atof(argv[26]));

	ainit1();

	if ((eing = fopen(argv[2], "r")) == NULL)
  {
      fprintf(stderr,"File \"%s\" not found !\n",argv[2]);
      exit(1);
  }
	fprintf(stderr,"\n Experiment : \"%s\"\n\n",argv[2]);

	pardet();
	n00=0;

	fclose(eing);
	nc = floor(clmm/cstep)+1;
	clm = nc*cstep;
	nthmin = floor(fabs(aavmin/astep))+1;
	amin = -astep*nthmin;
	nthmax = floor(fabs(aavmax/astep))+1;
	amax = astep*nthmax;
	nth = nthmin+nthmax;


	fprintf(stderr,"\n clm : %3.3f, amin : %3.3f, amax: %3.3f "
			//" nc : %d, nthmin : %d, nthmax : %d, nth : %d
			"\n",
					clm, amin, amax
							//, nc, nthmin, nthmax, nth
				 );
	//scanf("%f", &dum);

	ainit2();

	clmm=aavmin=aavmax=0.; clmerr=ihisterr=0;	nrtt=0; wl=sl=0;

	if ((eing = fopen(argv[2], "r")) == NULL)
  {
      fprintf(stderr,"File \"%s\" not found !\n",argv[2]);
      exit(1);
  }
	//fprintf(stderr,"\n Experiment : \"%s\"\n\n",argv[2]);

	for (k=0; k< nt; k++) /*Track number*/
  {
		fgets(bufein,24*(np+4),eing);
		ii=0;
		//fprintf(stderr," Track : %u \n %s \n", k+1, bufein);
		for(fstnumb=strtok(bufein,del);fstnumb;fstnumb=strtok(NULL,del))
		{
			xm=xx;
			xx=atof(fstnumb);
			xiyi[ii]=xx;
			ii++;
			//fprintf(stderr," Track : %u, Point : %u, %3.6f  \n", k+1, ii, xx);
		}
		if(((ii/2)*2!=ii && xx!=0)|| ((ii/2)*2==ii && xx==0))
		{
			fprintf(stderr," different numbers of x and y values in track "
					" %u  :\n", k+1);
			wl ++;
			//scanf("%f", &dum);
			continue;
		}
		if((ii/2)*2!=ii && xx==0)
		{
			//fprintf(stderr," space at end :\n");
			//fprintf(stderr," Track : %u, Point : %u, %3.6f %3.6f \n", k+1, ii, xm, xx);
			ii--;
			xx=xm;
			//fprintf(stderr," Track : %u, Point : %u, %3.6f %3.6f \n", k+1, ii, xm, xx);
			//scanf("%f", &dum);
		}
		if(xx==-1)
		{
			npix[k]=npiy[k]=(ii-2)/2;
		}
		else
		{
			npix[k]=npiy[k]=(ii)/2;
		}
		if(npix[k]<np) sl++;
		//scanf("%f", &dum);
		for (i=0; i<npix[k]; i++)
		{
			xi[i] = xiyi[i];
			if(xx==-1) yi[i] = xiyi[i+npix[k]+1];
			else yi[i] = xiyi[i+npix[k]];
			/*fprintf(stderr," Track : %u, Point : %u, x : %3.3f, y: %3.3f\n"
												, k+1, i+1, xi[i], yi[i]);*/
		}
		//scanf("%f", &dum);
		if(npix[k] != npiy[k])
		{
			/*fprintf(stderr," Track : %u, npx : %u, npy : %u\n",k+1,
							npix[k], npiy[k]);*/
			continue;
		}

		if(np1<npix[k])npix[k]=np1;
		// np1 used points
		fscanf(eing, "\n");
		/* start direction */

		if(npix[k] < nav)
		{
			fprintf(stderr," Track %u too short, npx : %u\n",k+1, npix[k]);
		}
		else
		{
			rtt=0.;

	  	for (i0=nav; i0<(npix[k]); i0++)
				/*i0=nav lower limit of the window for determination of the persistence length*/
			{ //i0i: averaging interval for the starting directon
				i0i=i0-nav;
				xm=ym=0.;
				for (j=i0i+1; j<=i0i+nav-1; j++)
				{
					xm +=xi[j]-xi[i0i];
					ym +=yi[j]-yi[i0i];
				} // starting directon
				if(sqrt(xm*xm+ym*ym)>0.)
				{
					xl0 = xm/sqrt(xm*xm+ym*ym);
					yl0 = ym/sqrt(xm*xm+ym*ym);
					if((fabs(sqrt(xl0*xl0+yl0*yl0))-1.)>round)
					{
						fprintf(stderr," Normalization !!\n");
						//scanf("%f", &dum);
						continue;
					}

				//persistence length histograms
				clij=0.;
				for (i=i0-nav+1; i<=i0; i++)
				clij += sqrt((xi[i]-xi[i-1])*(xi[i]-xi[i-1]) +
						          					(yi[i]-yi[i-1])*(yi[i]-yi[i-1]));
				// determine contour lengths for initial averaging interval
				for (i=i0+1; i< npix[k]; i++)
				{
					clij += sqrt((xi[i]-xi[i-1])*(xi[i]-xi[i-1]) +
						          											 	 (yi[i]-yi[i-1])*(yi[i]-yi[i-1]));
					icl =	floor(nc*(clij)/(clm));
					// index of contour lengths, bin centres
					if(clmm< clij) clmm=clij;

					aav = (xi[i]-xi[i0-nav])*xl0 + (yi[i]-yi[i0-nav])*yl0;
					ddi = (xi[i]-xi[i0-nav])*(xi[i]-xi[i0-nav])+
							(yi[i]-yi[i0-nav])*(yi[i]-yi[i0-nav]);
					actr[icl]+=aav;
					acqtr[icl]+=aav*aav;
					ddtr[icl]+=ddi;
					msdcltr[icl]+=sqrt(ddi)/clij;
					msdcl2tr[icl]+=ddi/(clij);

					if(aav<aavmin)aavmin=aav;
					if(aavmax<aav)aavmax=aav;
					if(aav>amin)ihist = floor(nth*(aav-amin)/(amax-amin));
					else ihist=0;
					// Histogram index
					/*fprintf(stderr,"Track : %u, i0 : %u, , i : %u , aav : %3.3f, "
								"  ihist : %u, xl0 :%3.3f, yl0 :%3.3f\n",
								k+1, i0, i, aav, ihist, xl0, yl0 );*/
								//scanf("%f", &dum);
					if(ihist<0 ||ihist>=nth || i>npix[k] || i<(i0+1))
					{
						fprintf(stderr,"Track : %u, i0 : %u, , i : %u , aav : %3.3f, "
								"  ihist : %u, xl0 :%3.3f, yl0 :%3.3f\n",
						k+1, i0, i, aav, ihist, xl0, yl0 );
						ihisterr ++;
						//scanf("%f", &dum);
					}
					else
					{
						if(icl<0||icl>(nc-1))
						{
							fprintf(stderr,"clm too small!!, Track:%u, i0:%u, i:%u, cl:%3.3f\n"
									, k+1, i0, i, clij);
							clmerr ++;
							//scanf("%f", &dum);
						}
						else
						{
							histctr[ihist*nc+icl]++;
							// index of contour length assigned to lines of the histogram matrix
							nhistctr[icl]++;
							//if((k/100)*100==k)
							/*fprintf(stderr," k : %u, i0 : %u, i : %u, icl : %u, ihist : %u "
										" aav : %3.3f, nhistctr: %3.3f, histctr : %3.3f\n",
							k, i0, i, icl, ihist, aav, nhistctr[icl], histctr[ihist*nc+icl]);*/
							//scanf("%f", &dum);
							if(aav>amax)
							{
								fprintf(stderr," aav>amax! aav : %3.3f, i0 %u, i %u, icl %u\n",
												aav, i0, i, icl);
								//scanf("%f", &dum);
							}

							if(kl[icl]<k) // Track k contributes to bin icl
							{
								kl[icl] = k;
								knt[icl]++;
							/*fprintf(stderr," k : %u, icl : %u, kl[icl] : %f, knt[icl] : %f\n",
												k, icl, kl[icl], knt[icl]);
							scanf("%f", &dum);*/
							}
						}
						//fprintf(stderr," i0 : %u, i : %u ihist : %u, icl : %u\n",
						//				i0, i, ihist, icl);
					}
				} /* i*/
				//if(i0==(i00+nav) && (k/10)*10==k) scanf("%f", &dum);
				} // valid start direction
				else
				{
					fprintf(stderr," Start direction: Track: %u, i0: %u x, y: %3.3f %3.3f\n"
					,k+1, i0, xm, ym);
					n00++;
					//scanf("%f", &dum);
				}
			} /* i0 */

	// rms only
	  	for (i0=0; i0<npix[k]; i0++)
			{
				clij=0.;
				for (i=i0+1; i<=i0; i++)
				clij += sqrt((xi[i]-xi[i-1])*(xi[i]-xi[i-1]) +
						          											 	 (yi[i]-yi[i-1])*(yi[i]-yi[i-1]));
				for (i=i0+1; i< npix[k]; i++)
				{
					clij += sqrt((xi[i]-xi[i-1])*(xi[i]-xi[i-1]) +
						          											 	 (yi[i]-yi[i-1])*(yi[i]-yi[i-1]));

					ddi = (xi[i]-xi[i0])*(xi[i]-xi[i0])+
							(yi[i]-yi[i0])*(yi[i]-yi[i0]);
					ddttr[i-i0]+=ddi;
					if(clij>0.)
					{
						rmsdclttr[i-i0]+=sqrt(ddi)/clij;
						msdclttr[i-i0]+=ddi/(clij);
						ihistdd = floor(ndd*sqrt(ddi)/clij);
					}
					else
					{
						rmsdclttr[i-i0]+=1.;
						msdclttr[i-i0]+=0.;
						ihistdd = ndd-1;
					}
					/*fprintf(stderr," Track: %u, i0: %u, i: %u, ihistdd, %u "
							" ddi: %3.3f, rmsdcl: %3.3f\n"
					,k+1, i0, i, ihistdd, ddi, sqrt(ddi)/clij);*/
					histddtr[ihistdd*np+i-i0]++;
					nhistddtr[i-i0]++;
				}
				//scanf("%f", &dum);
			} /* i0 */

		/* Determination of histogram values for the longest persistent path */
		} /*valid track*/

		for (kk=0; kk< nth; kk++)
		{
			for (ii=0; ii< nc; ii++)
			{
				if(nhistctr[ii]>0.)
				{
					histc[kk*nc+ii] += histctr[kk*nc+ii]/nhistctr[ii];
				}
			}
		}
		for (ii=0; ii< nc; ii++)
		{
			if(nhistctr[ii]>0.)
			{
				ac[ii]	+= actr[ii]/nhistctr[ii];
				acq[ii] += acqtr[ii]/nhistctr[ii];
				dd[ii]	+= ddtr[ii]/nhistctr[ii];
				msdcl[ii]+=msdcltr[ii]/nhistctr[ii];
				msdcl2[ii]+=msdcl2tr[ii]/nhistctr[ii];
			}
		}
		for (kk=0; kk< nth; kk++)
		{
			for (ii=0; ii< nc; ii++)
			{
				histctr[kk*nc+ii] = 0.;
				nhistctr[ii]= actr[ii]=acqtr[ii]=ddtr[ii]=msdcltr[ii]=msdcl2tr[ii]=0.;
			}
		}

		for (kk=0; kk< ndd; kk++)
		{
			for (ii=0; ii< np; ii++)
			{
				if(nhistddtr[ii]>0.)
				{
					histdd[kk*np+ii] += histddtr[kk*np+ii]/nhistddtr[ii];
				}
			}
		}
		for (ii=0; ii< np; ii++)
		{
			if(nhistddtr[ii]>0.)
			{
				ddt[ii]	+= ddttr[ii]/nhistddtr[ii];
				rmsdclt[ii] += rmsdclttr[ii]/nhistddtr[ii];
				msdclt[ii]	+= msdclttr[ii]/nhistddtr[ii];
				nnormt[ii] +=1.;
				/*fprintf(stderr," Track: %u, i-i0: %u, dd: %3.3f, rmsdcl: %3.3f, msdcl: %3.3f\n"
				,k+1, ii, ddttr[ii]/nhistddtr[ii], rmsdclttr[ii]/nhistddtr[ii],
	  		msdclttr[ii]/nhistddtr[ii]);*/
			}
		}
		//scanf("%f", &dum);
		for (kk=0; kk< ndd; kk++)
		{
			for (ii=0; ii< np; ii++)
			{
				histddtr[kk*nc+ii] = 0.;
				nhistddtr[ii]=ddttr[ii]=rmsdclttr[ii]=msdclttr[ii]=0.;
			}
		}


		//fprintf(stderr,"Track : %u, lstep : %f\n", k+1, lstep[k]);
	} /* k  Tracks */
  fclose(eing);

	for (ii=0; ii< nc; ii++) nhistc[ii] = 0.;
	for (kk=0; kk< nth; kk++)
	{
		for (ii=0; ii< nc; ii++)
		{
			nhistc[ii] += histc[kk*nc+ii];
		}
	}

	for (ii=0; ii< np; ii++) nhistdd[ii] = 0.;
	for (kk=0; kk< ndd; kk++)
	{
		for (ii=0; ii< np; ii++)
		{
			nhistdd[ii] += histdd[kk*np+ii];
		}
	}


	//Determination of maximum values
	for (i=0; i< nc; i++)
	{
		gewm = 0.;
		for (k=0; k< nth; k++)
		{
			if(histc[k*nc+i]> gewm)
			{
				acm[i]=k*((amax-amin)/nth)+amin;
				gewm = histc[k*nc+i];
			}
		}
	}

	fprintf(stderr,"\n");
	strcpy(bufe,argv[2]);
	strncat(bufe,"_a(cl)",6);
	ausg  = fopen(bufe, "w");
	for (k=0; k< nc; k++)
	{
		if(knt[k]>(tres*nt))
		{
			fprintf(stderr," cl: %3.3f, a(cl)): %3.3f, rms: %3.3f, rms/cl: %3.3f, "
					" ms/cl: %3.3f, fraction tracks: %3.3f %3.3f\n",
		 		k*(clm/nc)+0.5*cstep, ac[k]/nhistc[k], sqrt(dd[k]/nhistc[k]),
					 msdcl[k]/nhistc[k], msdcl2[k]/nhistc[k], knt[k]/nt, nhistc[k]/nt);
			fprintf(ausg," %e %e %e %e %e \n",
				k*(clm/nc)+0.5*cstep, ac[k]/nhistc[k], sqrt(dd[k]/nhistc[k]),
				msdcl[k]/nhistc[k], msdcl2[k]/nhistc[k]);
		}
	}
	fclose(ausg);
	fprintf(stderr,"\n");

	//scanf("%f", &dum);
	strcpy(bufe,argv[2]);
	strncat(bufe,"_speed",6);
	ausg  = fopen(bufe, "w");
	for (k=0; k< nt; k++)
	{
		fprintf(ausg," %u  %e\n", k, lstep[k]/dt);
	}
	fclose(ausg);

	strcpy(bufe,argv[2]);
	strncat(bufe,"_histc",6);
	ausg  = fopen(bufe, "w");
	fprintf(ausg," 0 ");
	//fprintf(stderr," 0 ");
	ncp = floor(clmp/cstep);
  for (i=0; i< ncp; i++)
	{
		fprintf(ausg," %e", (i+0.5)*cstep);
		//fprintf(stderr," %e", (i+0.5)*cstep);
	}

	fprintf(ausg," \n ");
	//fprintf(stderr," \n ");

  for (k=0; k< nth; k++)
	{
		aav = (k+0.5)*((amax-amin)/nth)+amin;
		//fprintf(stderr," %e,  %e,  %e", aav, aminp, amaxp);
		if(aav>aminp && aav<amaxp)
		{
			fprintf(ausg," %e", aav);
			//fprintf(stderr," %e ", aav);
			for (i=0; i< ncp; i++)
			{
				//if(knt[i]>(tres*nt))
				{
					fprintf(ausg," %e", histc[k*nc+i]/nhistc[i]);
					//fprintf(stderr," %e", histc[k*nc+i]/nhistc[i]);
				}
				/*else
				{
					fprintf(ausg," 0.");
					//fprintf(stderr," 0.");
				}*/
			}

		fprintf(ausg," \n ");
		//fprintf(stderr," \n ");
		}
		//fprintf(ausg," nhist : %u \n ", nhist);
	}
	fclose(ausg);


	strcpy(bufe,argv[2]);
	strncat(bufe,"_histt",6);
	ausg  = fopen(bufe, "w");
	fprintf(ausg," 0 ");
	//fprintf(stderr," 0 ");

  for (i=0; i< np; i++)
	{
		fprintf(ausg," %e", (i+1.)*dt);
		//fprintf(stderr," %e", (i+0.5)*cstep);
	}

	fprintf(ausg," \n ");
	//fprintf(stderr," \n ");

  for (k=0; k< ndd; k++)
	{
		rmsdcl = (k*1.)/ndd;
		//fprintf(stderr," %e,  %e,  %e", aav, aminp, amaxp);
			fprintf(ausg," %e", rmsdcl);
			//fprintf(stderr," %e ", aav);
			for (i=1; i< np; i++)
			{
					fprintf(ausg," %e", histdd[k*np+i]/nhistdd[i]);
					//fprintf(stderr," %e", histc[k*nc+i]/nhistc[i]);
			}

		fprintf(ausg," \n ");
	}
	fclose(ausg);



	cltotj=0.;
	ntotj=0;
	clmmj=0.;
	strcpy(bufe,argv[2]);
	strncat(bufe,"_(t)",4);
	ausg  = fopen(bufe, "w");
	for (i=i00+1; i< np1; i++)
	{
  	fprintf(stderr," step: %2u, time: %3.3f min, mean contour length: %3.3f,"
				" mean speed %3.3f um/min, rms/cl: %3.3f, ms/cl: %3.3f\n",
		i, i*dt, i*cltot/(nt*np1), cltoti[i]/(nltoti[i]*dt), rmsdclt[i]/nnormt[i],
					 msdclt[i]/nnormt[i]);
		fprintf(ausg," %e %e %e %e %e \n",
			i*dt, i*cltot/(nt*np1), cltoti[i]/(nltoti[i]*dt), rmsdclt[i]/nnormt[i],
					 msdclt[i]/nnormt[i]);
		clmmj+=cltoti[i]/nltoti[i];
		cltotj+=cltoti[i];
		ntotj+=nltoti[i];
	}
	fclose(ausg);

	fprintf(stderr,"\n mean step length : %3.3f, mean total length : %3.3f\n\n",
					 cltotj/ntotj, clmmj);
  //fprintf(stderr," mittlere Schrittlänge : %f \n", cltot/ntot);
  if(clmerr>0. ||ihisterr>0.)
	fprintf(stderr," Histogramm error: clmerr : %u , ihisterr %u \n\n",
					clmerr, ihisterr);
	if(n00>0) fprintf(stderr," %u undetermined start directions \n", n00);
	if(wl>0) fprintf(stderr," %u wrong input lines \n", wl);
	if(sl>0) fprintf(stderr," %u short input lines \n", sl);

	fprintf(stderr,"\n Experiment : \"%s\"\n\n",argv[2]);

	return(0);
}

void ainit1()
{
  if ((npix = (unsigned int*)calloc(nt, sizeof(unsigned int))) == NULL)
  {
    fprintf(stderr,"not enough memory for *npix !\n");
    exit(1);
  }
  if ((npiy = (unsigned int*)calloc(nt, sizeof(unsigned int))) == NULL)
  {
    fprintf(stderr,"not enough memory for *npiy !\n");
    exit(1);
  }
  if ((lstep = (double*)calloc(nt, sizeof(double))) == NULL)
  {
    fprintf(stderr,"not enough memory for *lstep !\n");
    exit(1);
  }
	for (k=0; k< nt; k++) lstep[k]=0.;
  if ((xi = (double*)calloc(np, sizeof(double))) == NULL)
  {
    fprintf(stderr,"not enough memory for *xi !\n");
    exit(1);
  }
  if ((yi = (double*)calloc(np, sizeof(double))) == NULL)
  {
    fprintf(stderr,"not enough memory for *yi !\n");
    exit(1);
  }
	if ((cltoti = (double*)calloc(np, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *cltoti !\n");
    exit(1);
  }
	if ((nltoti = (unsigned int*)calloc(np, sizeof(unsigned int*))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *cltoti !\n");
    exit(1);
  }
	for (k=0; k< np; k++)
	{
		cltoti[k] = 0.;
		nltoti[k] = 0;
	}
	if ((nprt = (unsigned int*)calloc(nt, sizeof(unsigned int*))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *nprt !\n");
    exit(1);
  }
	if ((bufein = (char*)calloc(24*(np+4), sizeof(char*))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *bufein !\n");
    exit(1);
  }
	if ((xiyi = (double*)calloc(2*(np+2), sizeof(double))) == NULL)
  {
    fprintf(stderr,"not enough memory for *xiyi !\n");
    exit(1);
  }
}

void ainit2()
{
	if ((histc = (double*)calloc((nc+1)*(nth+1), sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *histc !\n");
    exit(1);
  }
	if ((histctr = (double*)calloc((nc+1)*(nth+1), sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *histctr !\n");
    exit(1);
  }
	if ((histdd = (double*)calloc((np+1)*(nth+1), sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *histdd !\n");
    exit(1);
  }
	if ((histddtr = (double*)calloc((np+1)*(nth+1), sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *histddtr !\n");
    exit(1);
  }
	for (k=0; k< (nc+1)*(nth+1); k++) histc[k] = histctr[k] = 0.;
	for (k=0; k< (np+1)*(nth+1); k++) histdd[k] = histddtr[k] = 0.;
	if ((nhistdd = (double*)calloc(np+1, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *nhistdd !\n");
    exit(1);
  }
	if ((nnormt = (double*)calloc(np+1, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *nnormt !\n");
    exit(1);
  }
	if ((nhistddtr = (double*)calloc(np+1, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *nhistddtr !\n");
    exit(1);
  }
	if ((nhistc = (double*)calloc(nc+1, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *nhist !\n");
    exit(1);
  }
	if ((nhistctr = (double*)calloc(nc+1, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *nhisttr !\n");
    exit(1);
  }

	if ((noro = (double*)calloc(nc+1, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *noro !\n");
    exit(1);
  }
	for (k=0; k< nc+1; k++) nhistc[k] = nhistctr[k] = noro[k] = 0.;
	for (k=0; k< np+1; k++) nhistdd[k] = nhistddtr[k] = nnormt[k] = 0.;

	if ((ac = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *ac !\n");
    exit(1);
  }
	if ((acq = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *acq !\n");
    exit(1);
  }
	if ((acm = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *acm !\n");
    exit(1);
  }
	if ((dd = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *dd !\n");
    exit(1);
  }
	if ((actr = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *actr !\n");
    exit(1);
  }
	if ((acqtr = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *acqtr !\n");
    exit(1);
  }
	if ((acmtr = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *acmtr !\n");
    exit(1);
  }
	if ((ddtr = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *ddtr !\n");
    exit(1);
  }
		if ((msdcltr = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *mscdcltr !\n");
    exit(1);
  }
		if ((msdcl = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *mscdcl !\n");
    exit(1);
  }
		if ((msdcl2tr = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *mscdcl2tr !\n");
    exit(1);
  }
		if ((msdcl2 = (double*)calloc(nc, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *mscdcl2 !\n");
    exit(1);
  }

	if ((ddttr = (double*)calloc(np, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *ddttr !\n");
    exit(1);
  }
	if ((ddt = (double*)calloc(np, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *ddt !\n");
    exit(1);
  }
	if ((rmsdclttr = (double*)calloc(np, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *rmsdclttr !\n");
    exit(1);
  }
	if ((msdclttr = (double*)calloc(np, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *msdclttr !\n");
    exit(1);
  }
	if ((rmsdclt = (double*)calloc(np, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *rmsdclt !\n");
    exit(1);
  }
	if ((msdclt = (double*)calloc(np, sizeof(double))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *msdclt !\n");
    exit(1);
  }

	for (k=0; k< nc; k++) ac[k] = dd[k] = acq[k] = actr[k] = ddtr[k] = acqtr[k]=
				msdcl[k]=msdcltr[k]=msdcl2[k]=msdcl2tr[k]=0.;
	for (k=0; k< np; k++) ddttr[k] = rmsdclttr[k] = msdclttr[k] = ddt[k] =
				 rmsdclt[k] = msdclt[k] =0.;

	if ((kl = (float*)calloc(nc, sizeof(float))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *kl !\n");
    exit(1);
  }
	if ((knt = (float*)calloc(nc, sizeof(float))) == NULL)
  {
  	fprintf(stderr,"not enough memory for *knt !\n");
    exit(1);
  }
	for (k=0; k< nc; k++)
	{
		kl[k] = -1.;
		knt[k] = 0.;
	}
}

#include "pardet_26mm.c"

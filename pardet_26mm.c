void pardet()
{
	cltot=clmm=aavmin=aavmax=ddmax=0.; ntot=clmerr=ihisterr=0;

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
			/*fprintf(stderr," different numbers of x and y values in track "
					" %u  :\n", k+1);
			wl ++;*/
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
		//if(npix[k]<np) sl++;
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
			//fprintf(stderr," Track : %u, npx : %u, npy : %u\n",k+1, npix[k], npiy[k]);
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
			for (i=1; i<npix[k]; i++)
			{
				cltot+=sqrt((xi[i]-xi[i-1])*(xi[i]-xi[i-1])+
										(yi[i]-yi[i-1])*(yi[i]-yi[i-1]));
				ntot++;
				cltoti[i]+=sqrt((xi[i]-xi[i-1])*(xi[i]-xi[i-1])+
												(yi[i]-yi[i-1])*(yi[i]-yi[i-1]));
				lstep[k]+=sqrt((xi[i]-xi[i-1])*(xi[i]-xi[i-1])+
												(yi[i]-yi[i-1])*(yi[i]-yi[i-1]));
				nltoti[i]++;
				//fprintf(stderr," Track : %u, Point : %u\n", k+1, i);
			}
			lstep[k]/=npix[k];
			//scanf("%f", &dum);
			//rtt=0.;

	  	for (i0=nav; i0<(npix[k]); i0++)
				/*i0=nav lower limit of the window for determination of the persistence length*/
			{ //i0i: averaging interval for the starting directon
				i0i=i0-nav;
				xm=ym=0.;
				for (j=i0i+1; j<=i0i+nav-1; j++)
				{
					xm +=xi[j]-xi[i0i];
					ym +=yi[j]-yi[i0i];
				} // starting direction
				if(sqrt(xm*xm+ym*ym)>0.)
				{
					xl0 = xm/sqrt(xm*xm+ym*ym);
					yl0 = ym/sqrt(xm*xm+ym*ym);
					if((fabs(sqrt(xl0*xl0+yl0*yl0))-1.)>round)
					{
						fprintf(stderr," Normalization !!\n");
						scanf("%f", &dum);
					}

				// persistence length histograms
				clij=0.;
				for (i=i0-nav+1; i<=i0; i++)
				clij += sqrt((xi[i]-xi[i-1])*(xi[i]-xi[i-1]) +
						          											 	 (yi[i]-yi[i-1])*(yi[i]-yi[i-1]));
				// determine accumulated distance for initial averaging interval
				for (i=i0+1; i< npix[k]; i++) //
				{
					clij += sqrt((xi[i]-xi[i-1])*(xi[i]-xi[i-1]) +
						          											 	 (yi[i]-yi[i-1])*(yi[i]-yi[i-1]));
					icl =	floor(nc*(clij)/(clm));
					/* index of accumulated distance (cl), bin centres */
					if(clmm< clij) clmm=clij;

					aav = (xi[i]-xi[i0-nav])*xl0 + (yi[i]-yi[i0-nav])*yl0;
					ddi = (xi[i]-xi[i0-nav])*(xi[i]-xi[i0-nav])+
							(yi[i]-yi[i0-nav])*(yi[i]-yi[i0-nav]);
					if(aav<aavmin)aavmin=aav;
					if(aavmax<aav)aavmax=aav;
					if(ddmax<ddi)ddmax=ddi;
				}// i
				}
				else
				{
					fprintf(stderr," start direction! track: %u, i0: %u x, y: %3.3f %3.3f\n"
					,k+1, i0, xm, ym);
					//scanf("%f", &dum);
				}
			} // i0
		} //valid track
	} //k  Tracks

	/*fprintf(stderr,"\n clmm : %f3.3, aavmin : %f3.3, aavmax: %f3.3\n",
					clmm, aavmin, aavmax);*/
	//scanf("%f", &dum);
}

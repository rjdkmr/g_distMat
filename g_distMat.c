/*
 * This file is part of g_distMat
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2014  Rajendra Kumar
 *
 * g_distMat is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * g_distMat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with g_distMat.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#include <math.h>
#include <string.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>


#include "vec.h"
#include "typedefs.h"
#include "filenm.h"
#include "statutil.h"
#include "futil.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "index.h"
#include "tpxio.h"
#include "rmpbc.h"
#include "pbc.h"

void *status;
pthread_t *thread;
pthread_attr_t attr;
pthread_mutex_t dist_mutex;
int ithread;

typedef struct {
	int nA, nB, nframes;
	int *resndxA, *resndxB, *natmresA, *natmresB;
	real cutoff, **mean, **var, **cmap, **sumdist, **sumsqdist;
	rvec *coord;
	int nthreads;
} DIST_MAT;

DIST_MAT distance_matrix;

void CopyRightMsg() {

    const char *copyright[] = {
            "                                                                        ",
            "                           :-)  g_distMat (-:                           ",
            "                                                                        ",
            "                         Author: Rajendra Kumar                         ",
            "                                                                        ",
            "                    Copyright (C) 2014  Rajendra Kumar                  ",
            "                                                                        ",
            "                                                                        ",
            "g_distMat is a free software: you can redistribute it and/or modify     ",
            "it under the terms of the GNU General Public License as published by    ",
            "the Free Software Foundation, either version 3 of the License, or       ",
            "(at your option) any later version.                                     ",
            "                                                                        ",
            "g_distMat is distributed in the hope that it will be useful,             ",
            "but WITHOUT ANY WARRANTY; without even the implied warranty of          ",
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ",
            "GNU General Public License for more details.                            ",
            "                                                                        ",
            "You should have received a copy of the GNU General Public License       ",
            "along with g_distMat.  If not, see <http://www.gnu.org/licenses/>.       ",
            "                                                                        ",
            "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     ",
            "\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     ",
            "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   ",
            "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT    ",
            "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   ",
            "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED",
            "TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  ",
            "PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  ",
            "LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    ",
            "NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      ",
            "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            ",
            "                                                                        ",
            "                           :-)  g_distMat (-:                       ",
            "                                                                        ",
            "                                                                        "
    };
    int i = 0;
    char *str;
    for(i=0; i<35; i++) {
        str = strdup(copyright[i]);
        fprintf(stderr,"%s\n", str);
    }
}



void set_num_threads(int nt)	{
	distance_matrix.nthreads = nt;
}


void make_thread_index()	{
	long i = 0, start, end;

	int NTHREADS = distance_matrix.nthreads;

	for(i=0; i<NTHREADS; i++){
		if ((distance_matrix.nA%NTHREADS) < 1)	{
			start = ceil(i*((distance_matrix.nA/NTHREADS)));
			end = start + ceil(distance_matrix.nA/NTHREADS);
		}
		else	{
			start = ceil(i*((distance_matrix.nA/NTHREADS)));
			end = start + ceil(distance_matrix.nA/NTHREADS);
			if (i==NTHREADS-1)
				end = start + ceil(distance_matrix.nA/NTHREADS)+distance_matrix.nA%NTHREADS;
		}
		//end = ceil((i+1)*(distance_matrix.nA/NTHREADS))+NTHREADS;
		fprintf(stdout,"\n%ld %ld\n",start, end);
	}
}

// Function to calculate minimum distance between two residues
real calc_min_dist (int resi, int *resndxA, int *natmresA, int resj, int *resndxB, int *natmresB, rvec *x)
{
	real dist=999.0, r;
	int i, j,ii,jj;

	for(i=0;i<natmresA[resi];i++)
	{
		ii = resndxA[resi] + i;
		for(j=0;j<natmresB[resj];j++)
		{
			jj = resndxB[resj] + j;
			r = ((x[ii][XX]-x[jj][XX])*(x[ii][XX]-x[jj][XX]))+
			    ((x[ii][YY]-x[jj][YY])*(x[ii][YY]-x[jj][YY])) +
			    ((x[ii][ZZ]-x[jj][ZZ])*(x[ii][ZZ]-x[jj][ZZ]));

			r = sqrt(r);

			if (r<dist)
				dist = r;
		}
	}
	return dist;
}

// Initialization of residue indexing for group A
void init_ndx_distmat_grp_A(t_topology top, int *index, int isize)	{
	int i=0, nres=0, prevres;

	prevres = top.atoms.atom[index[0]].resind;
	nres  = 1;

	snew(distance_matrix.natmresA, nres);
	distance_matrix.natmresA[nres-1] = 0;

	snew(distance_matrix.resndxA, nres);

	for(i=0;i<isize;i++)
	{
		if(top.atoms.atom[index[i]].resind == prevres)	{
			distance_matrix.natmresA[nres-1]++;
		}

		if(top.atoms.atom[index[i]].resind != prevres)
		{
			nres++;

			srenew(distance_matrix.natmresA, nres);
			srenew(distance_matrix.resndxA, nres);

			distance_matrix.natmresA[nres-1] = 0;
			distance_matrix.natmresA[nres-1]++;
			distance_matrix.resndxA[nres-1] = index[i];

			prevres = top.atoms.atom[index[i]].resind;
		}
	}
	distance_matrix.nA = nres;
	fprintf(stderr,"There are %d residues with %d atoms in first group\n",nres,isize);
}

// Initialization of residue indexing for group B
void init_ndx_distmat_grp_B(t_topology top, int *index, int isize)	{
	int i=0, nres=0, prevres;

	prevres = top.atoms.atom[index[0]].resind;
	nres  = 1;

	snew(distance_matrix.natmresB, nres);
	distance_matrix.natmresB[nres-1] = 0;

	snew(distance_matrix.resndxB, nres);

	for(i=0;i<isize;i++)
	{
		if(top.atoms.atom[index[i]].resind == prevres)	{
			distance_matrix.natmresB[nres-1]++;
		}

		if(top.atoms.atom[index[i]].resind != prevres)
		{
			nres++;

			srenew(distance_matrix.natmresB, nres);
			srenew(distance_matrix.resndxB, nres);

			distance_matrix.natmresB[nres-1] = 0;
			distance_matrix.natmresB[nres-1]++;
			distance_matrix.resndxB[nres-1] = index[i];

			prevres = top.atoms.atom[index[i]].resind;
		}
	}
	distance_matrix.nB = nres;
	fprintf(stderr,"There are %d residues with %d atoms in second group\n",nres,isize);
}

// Initialization of output data
void init_outdata_distmat()	{
	int i=0, j=0;

	snew(distance_matrix.sumdist, distance_matrix.nA);
	snew(distance_matrix.sumsqdist, distance_matrix.nA);
	snew(distance_matrix.mean, distance_matrix.nA);
	snew(distance_matrix.var, distance_matrix.nA);
	snew(distance_matrix.cmap, distance_matrix.nA);


	for(i=0; i<distance_matrix.nA; i++){
		snew(distance_matrix.sumdist[i], distance_matrix.nB);
		snew(distance_matrix.sumsqdist[i], distance_matrix.nB);
		snew(distance_matrix.mean[i], distance_matrix.nB);
		snew(distance_matrix.var[i], distance_matrix.nB);
		snew(distance_matrix.cmap[i], distance_matrix.nB);


		for(j=0; j<distance_matrix.nB; j++){
			distance_matrix.sumdist[i][j] = 0;
			distance_matrix.sumsqdist[i][j] = 0;
			distance_matrix.mean[i][j] = 0;
			distance_matrix.var[i][j] = 0;
			distance_matrix.cmap[i][j] = 0;
		}

	}
	distance_matrix.nframes = 0;
}

void calculate_dist_mat() {
	int i=0, j=0;
	real dist=0;

	for(i=0; i<distance_matrix.nA; i++){
		for(j=0; j<distance_matrix.nB; j++){
			 dist = calc_min_dist(i, distance_matrix.resndxA, distance_matrix.natmresA, j, distance_matrix.resndxB, distance_matrix.natmresB, distance_matrix.coord);
			 distance_matrix.sumdist[i][j] += dist;
			 distance_matrix.sumsqdist[i][j] += (dist*dist);
		}
	}
}


void *calculate_dist_mat_pthread (void *arg) {
	int NTHREADS = distance_matrix.nthreads;
	long t;
	t = (long)arg;
	int start, end;

	if ((distance_matrix.nA%NTHREADS) < 1)	{
		start = ceil(t*((distance_matrix.nA/NTHREADS)));
		end = start + ceil(distance_matrix.nA/NTHREADS);
	}
	else	{
		start = ceil(t*((distance_matrix.nA/NTHREADS)));
		end = start + ceil(distance_matrix.nA/NTHREADS);
		if (t==NTHREADS-1)
			end = start + ceil(distance_matrix.nA/NTHREADS)+distance_matrix.nA%NTHREADS;
	}

	int i=start, j=0;
	real dist=0;

	if(NTHREADS==1){
		start = 0;
		end = distance_matrix.nA;
	}
	if(end>distance_matrix.nA)
		end=distance_matrix.nA;


	for(i=start; i<end; i++){
		for(j=0; j<distance_matrix.nB; j++){
			 dist = calc_min_dist(i, distance_matrix.resndxA, distance_matrix.natmresA, j, distance_matrix.resndxB, distance_matrix.natmresB, distance_matrix.coord);
			 pthread_mutex_lock(&dist_mutex);

			 distance_matrix.sumdist[i][j] += dist;

			 distance_matrix.sumsqdist[i][j] += (dist*dist);

			 if (dist <= distance_matrix.cutoff)
				 distance_matrix.cmap[i][j]++;

			 pthread_mutex_unlock(&dist_mutex);
		}
	}
	if(NTHREADS>1)
		pthread_exit(NULL);
	return (void *)0;
}

void calc_mean_var()	{
	int i=0, j=0, n = distance_matrix.nframes;
	real sum=0, sumsq=0;
	for(i=0; i<distance_matrix.nA; i++){
		for(j=0; j<distance_matrix.nB; j++){
			sum = distance_matrix.sumdist[i][j];
			sumsq = distance_matrix.sumsqdist[i][j];

			distance_matrix.mean[i][j] = sum/n ;
			distance_matrix.var[i][j] = (sumsq - ( (sum*sum) /n)) / (n-1);
			distance_matrix.cmap[i][j] /= n;

		}

	}
}


int gmx_distMat(int argc,char *argv[])
{
	const char *desc[] = {
			"It calculates average minimum-distance matrix of residues between two atom-groups.",
			"Simultaneously, it calculates variance  and standard-deviation matrices.",
			"Also, it calculates fraction of contact-map over entire trajectory for the ",
			"residues that are within a minimum distance of \"-ct\" option value.\n\n",
			"To speed up the calculation, it uses all available cores of the CPU using",
			"multi-threading. Number of threads/cores could be change by \"-nt\" option."
	};
	static real cutoff = 1.0;
	int NTHREADS = sysconf( _SC_NPROCESSORS_ONLN );
	t_pargs pa[] = {
			{ "-ct",   FALSE, etREAL, {&cutoff}, "cut-off distance (nm) for contact map" },
			{ "-nt",  FALSE, etINT, {&NTHREADS}, "number of threads for multi-threading" }
	};

	t_filenm   fnm[] = {
			{ efTRX, "-f",  NULL, ffREAD },
			{ efTPS, NULL,  NULL, ffREAD },
			{ efNDX, NULL,  NULL, ffOPTRD },
			{ efDAT, "-mean", "average", ffWRITE },
			{ efDAT, "-var", "variance", ffOPTWR },
			{ efDAT, "-std", "stdeviation", ffOPTWR },
			{ efDAT, "-cmap", "contact_map", ffOPTWR }
	};
#define NFILE asize(fnm)

	FILE *fMean, *fVar, *fStd, *fCmap;
	t_topology top;
	int ePBC;

	int isizeA, isizeB;
	atom_id *indexA, *indexB;
	char *grpnameA, *grpnameB;

	int i, j, trxnat;
	long nt;
	t_trxstatus *trjstatus;
	rvec *x;
	char title[256];

	real time;
	matrix box;
	output_env_t oenv;
	gmx_rmpbc_t gpbc = NULL;

	// Copyright message
	CopyRightMsg();

	// Parse command line argument and print all options
	parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);

	// Set number of threads
	set_num_threads(NTHREADS);

	// Set cut-off distance for contact map
	distance_matrix.cutoff = cutoff;

	// Reading tpr filr
	read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&x,NULL,box,FALSE);

	// Selection of first index group
	fprintf(stderr,"Select first group:\n");
	get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&isizeA,&indexA,&grpnameA);

	// Selection of second index group
	fprintf(stderr,"Select second group:\n");
	get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&isizeB,&indexB,&grpnameB);

	// Initialization for first group
	init_ndx_distmat_grp_A(top, indexA, isizeA);

	// Initialization for second group
	init_ndx_distmat_grp_B(top, indexB, isizeB);

	// Initialization for output arrays
	init_outdata_distmat();

	trxnat=read_first_x(oenv, &trjstatus,ftp2fn(efTRX,NFILE,fnm), &time, &distance_matrix.coord, box);
	gpbc = gmx_rmpbc_init(&top.idef,ePBC,trxnat,box);

	do {
		distance_matrix.nframes++;

		if(NTHREADS>1)	{

			//make_thread_index();
			//exit(1);

			thread = malloc(sizeof(pthread_t) * NTHREADS);
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			pthread_mutex_init(&dist_mutex, NULL);

			for(nt=0;nt<NTHREADS;nt++)
				ithread = pthread_create(&thread[nt],&attr, calculate_dist_mat_pthread, (void *) nt);
				//mi_calculate((void *)&i);

			 //Free attribute and wait for the other threads
			pthread_attr_destroy(&attr);
			for(nt=0; nt<NTHREADS; nt++) {
				ithread = pthread_join(thread[nt], &status);
				if (ithread) {
					printf("ERROR; return code from pthread_join() is %d\n", ithread);
					exit(-1);
				 }
			 }
		}

		else
			calculate_dist_mat();

	}while (read_next_x(oenv,trjstatus, &time, trxnat, distance_matrix.coord, box));

	// Final calculation of average, var, std-deviation and contact-map
	calc_mean_var();

	//Writing output file
	fMean = opt2FILE("-mean", NFILE, fnm, "w");
	fVar = opt2FILE("-var", NFILE, fnm, "w");
	fStd = opt2FILE("-std", NFILE, fnm, "w");
	fCmap = opt2FILE("-cmap", NFILE, fnm, "w");

	for(i=0; i<distance_matrix.nA; i++){
		for(j=0; j<distance_matrix.nB; j++){
			fprintf(fMean, "%5.3f ", distance_matrix.mean[i][j]);
			fprintf(fVar, "%5.3f ", distance_matrix.var[i][j]);
			fprintf(fStd, "%5.3f ", sqrt(distance_matrix.var[i][j]));
			fprintf(fCmap, "%5.3f ", distance_matrix.cmap[i][j]);
		}
		fprintf(fMean,"\n");
		fprintf(fVar,"\n");
		fprintf(fStd,"\n");
		fprintf(fCmap,"\n");
	}

	fprintf(stdout, "\n\nThanks for using g_distMat!!!\n");
	return 0;
}

int main (int argc,char *argv[])	{
	gmx_distMat(argc,argv);
	return 0;
}

//#define mu(f)				(1./((f) / mu1 + (1. - (f)) / mu2))
//#define rho(f)			(1./((f) / rho1 + (1. - (f)) / rho2))
//#define FILTERED			1

#include "constants.h"

#if dimension == 3
#include "lambda2.h"
#endif

double drop_time_1file, bubble_time_1file, writefile_time_1file, simulation_time_1file;
double drop_time_total, bubble_time_total, writefile_time_total, simulation_time_total;
clock_t simulation_str_time, simulation_end_time;

static void remove_droplets(scalar c, bool droplet, double dropsize);
int refinebetweentwofluid(const double xmin, const double xmax, const double ymin, const double ymax, scalar f1, scalar f2, const int lvl);

struct CFDValues cfdbv;

// runnning example: ./file n6 x12 r3000 w500 h0.5 ts0.01 te0.250
// this means level min = 6, level max = 12, Reynolds = 3000, Weber = 500, H (pool depth) = 0.5*(DROPDIAMETER)

int main(int argc, char **argv)
{
	simulation_str_time = clock();
	simulation_time_1file = 0.0;
	writefile_time_1file = 0.0;
	drop_time_1file = 0.0;
	bubble_time_1file = 0.0;
	numericalmainvalues(argv, argc, &cfdbv);
	;
	size(cfdbv.domainsize);
#if AXI
	;
#else
	origin(0, -cfdbv.domainsize / 2., -cfdbv.domainsize / 2.);
#endif
	int initialgrid = pow(2, LEVELmin);
	init_grid(initialgrid);
	;
	rho1 = cfdbv.rhoL;
	rho2 = cfdbv.rhoG;
	mu1 = cfdbv.muL;
	mu2 = cfdbv.muG;
	f.sigma = cfdbv.Sigma;
	;
	TOLERANCE = 1e-6;
	run();
	;
	return 1;
}

event defaults(i = 0)
{
	interfaces = list_add(NULL, f);
	interfaces = list_add(interfaces, fdrop);
}

event init(i = 0)
{
	if (restore(file = FILENAME_LASTFILE))
	{
#if AXI
		boundary((scalar *){fm});
		//boundary({p});
/*
		scalar * dlist = dump_list (all);
		int varNO = list_len(dlist);
		printf("%d", varNO);
		for (scalar s in dlist)
			printf(", %s", s.name);
		printf("\r\n");
		free(dlist);
*/
#endif
	}
	else
	{
		double x0 = cfdbv.pooldepth + cfdbv.initialdis + cfdbv.diameter * 0.50;
		;
		refine(sq(x - x0) + sq(y) + sq(z) < sq(0.50 * cfdbv.diameter + cfdbv.refinegap) && sq(x - x0) + sq(y) + sq(z) > sq(0.50 * cfdbv.diameter - cfdbv.refinegap) && level < LEVELmax);
		refine((x < (cfdbv.pooldepth + 2.0 * cfdbv.refinegap) && x > (cfdbv.pooldepth - 2.0 * cfdbv.refinegap)) && level < LEVELmax);
		;
		foreach ()
		{
			fdrop[] = 0.0;
			f[] = 0.0;
			if (sq(x - x0) + sq(y) + sq(z) < sq(0.50 * cfdbv.diameter))
			{
				fdrop[] = 1.0;
				f[] = 1.0;
				u.x[] = -cfdbv.vel;
			}
			else if (x <= cfdbv.pooldepth)
				f[] = 1.0;
		};
		clock_t timestr, timeend;
		timestr = clock();
		FILE *fp;
		char name[100], tmp[50];
		sprintf(name, FILENAME_DURATION);
		sprintf(tmp, "-CPU%02d.plt", pid());
		strcat(name, tmp);
		fp = fopen(name, "w");
		fprintf(fp, "Variables = Iteration DeltaTime CriticalTime PhysicalTime LastDuration DropDuration BubbleDuration FileDuration CellNumber TotalLastDuration TotalDropDuration TotalBubbleDuration TotalFileDuration\r\nzone\r\n");
		fclose(fp);
		timeend = clock();
		writefile_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
	}
}
/*
event acceleration(i++)
{
//	FroudeNumber = VELOCITY / sqrt(GRAVITY * DROP_DIAMETER);
	if (fabs(GRAVITY) > 1.0e-9 && fabs(VELOCITY) > 1.0e-9)
	{
		face vector av = a;
		foreach_face(x)
			av.x[] -= GRAVITY * DROP_DIAMETER / (VELOCITY * VELOCITY);
	}
}
*/
#if REMOVE_DROP_YESNO == 'y'
event drop_remove(i += REMOVE_DROP_PERIOD)
{
	clock_t timestr, timeend;
	timestr = clock();
	if (t > cfdbv.timecontact)
	{
		remove_droplets(f, 1, REMOVE_DROP_SIZE);
		remove_droplets(fdrop, 1, REMOVE_DROP_SIZE);
	}
	timeend = clock();
	drop_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
}
#endif

#if REMOVE_BUBBLE_YESNO == 'y'
event drop_remove(i += REMOVE_BUBBLE_PERIOD)
{
	clock_t timestr, timeend;
	timestr = clock();
	if (t > cfdbv.timecontact)
	{
		remove_droplets(f, 0, REMOVE_BUBBLE_SIZE);
		remove_droplets(fdrop, 0, REMOVE_BUBBLE_SIZE);
	}
	timeend = clock();
	bubble_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
}
#endif

event adapt(i++)
{
	double refine[4];
	refine[0] = pow(10.0, REFINE_VALUE_0);
	refine[1] = pow(10.0, REFINE_VALUE_1);
	refine[2] = pow(10.0, REFINE_VALUE_2);
	refine[3] = pow(10.0, REFINE_VALUE_3);
	adapt_wavelet(REFINE_VAR, (double[]){refine[0], refine[1], refine[2], refine[3]}, maxlevel = LEVELmax, minlevel = LEVELmin);
	;
	//refinebetweentwofluid(0.95 * cfdbv.pooldepth, 1.05 * cfdbv.pooldepth, 0.0, 0.50 * cfdbv.diameter, f, fdrop, 13);
}

event showiteration(i++)
{
	switch (pid())
	{
	case 0:
	{
		char name[500], tmp[100];
		if (t - cfdbv.timecontact < 0.0)
			sprintf(name, "i%05d_dt%.2e_tb%.3f_P%02d", i, dt, t - cfdbv.timecontact, (int)(100.0 * t / MAX_TIME));
		else
			sprintf(name, "i%05d_dt%.2e_ta%.3f_P%02d", i, dt, t - cfdbv.timecontact, (int)(100.0 * t / MAX_TIME));
		sprintf(tmp, "_Re%d_We%d", (int)cfdbv.Reynolds, (int)cfdbv.Weber);
		strcat(name, tmp);
#if AXI
		sprintf(tmp, "_AXI");
		strcat(name, tmp);
#else
#if dimension == 3
		sprintf(tmp, "_3D");
		strcat(name, tmp);
#else
		sprintf(tmp, "_2D");
		strcat(name, tmp);
#endif
#endif
		sprintf(tmp, "_L%02d%02d", LEVELmin, LEVELmax);
		strcat(name, tmp);
#if REMOVE_DROP_YESNO == 'y'
		sprintf(tmp, "_RD%.1fP%02d", REMOVE_DROP_SIZE, REMOVE_DROP_PERIOD);
		strcat(name, tmp);
#endif
#if REMOVE_BUBBLE_YESNO == 'y'
		sprintf(tmp, "_RB%.1fP%02d", REMOVE_BUBBLE_SIZE, REMOVE_BUBBLE_PERIOD);
		strcat(name, tmp);
#endif
		printf("%s\r\n", name);
	}
	}
}

event end(t = cfdbv.timecontact + cfdbv.timeend)
{
	FILE *fp;
	char name[500], tmp[100];
	sprintf(name, FILENAME_ENDOFRUN);
	sprintf(tmp, "-CPU%02d.txt", pid());
	strcat(name, tmp);
	fp = fopen(name, "w");
	fprintf(fp, "SimulationTime %e\r\nRemoveDropTime %e\r\nRemoveBubbleTime %e\r\nWriteFileTime %e\r\n", simulation_time_total, drop_time_total, bubble_time_total, writefile_time_total);
	fclose(fp);
}

event outputfiles(t += SAVE_FILE_EVERY)
{
	clock_t timestr, timeend;
	timestr = clock();
	static FILE *fp;
	char name[500], tmp[100];
	;
	foreach ()
	{
		if (fdrop[] < 0.0)
			fdrop[] = 0.0;
		else if (fdrop[] > 1.0)
			fdrop[] = 1.0;
	};
	// write data file
	;
	foreach ()
		pressure[] = p[];
//	p.nodump = false;
	sprintf(name, "DataALL-%.4f", t);
	dump(file = name);
	;
	dump(file = FILENAME_LASTFILE);
	/*	;
	scalar * dlist = dump_list (all);
	int varNO = list_len(dlist);
	printf("%d", varNO);
	for (scalar s in dlist)
		printf(", %s", s.name);
	printf("\r\n");
	free(dlist);
*/;
	// duration file
	;
	timeend = clock();
	writefile_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
	int cellnumber = 0;
	foreach ()
		cellnumber++;
	simulation_end_time = clock();
	double estimatetimeleft;
	char LDc[100], TDc[100], ETLc[100];
	simulation_time_1file = (double)(simulation_end_time - simulation_str_time) / CLOCKS_PER_SEC;
	;
	simulation_time_total += simulation_time_1file;
	bubble_time_total += bubble_time_1file;
	drop_time_total += drop_time_1file;
	writefile_time_total += writefile_time_1file;
	if (t == 0.0)
		estimatetimeleft = 0.0;
	else
		estimatetimeleft = simulation_time_total * (cfdbv.timecontact + cfdbv.timeend) / t - simulation_time_total;
	timecalculation(simulation_time_1file, LDc);
	timecalculation(simulation_time_total, TDc);
	timecalculation(estimatetimeleft, ETLc);
	sprintf(name, FILENAME_DURATION);
	sprintf(tmp, "-CPU%02d.plt", pid());
	strcat(name, tmp);
	fp = fopen(name, "a");
	//fprintf(fp, "Variables = Iteration DeltaTime CriticalTime PhysicalTime LastDuration DropDuration BubbleDuration FileDuration CellNumber TotalLastDuration TotalDropDuration TotalBubbleDuration TotalFileDuration\r\nzone\r\n");
	fprintf(fp, "%d %e %e %e %e %e %e %e %d %e %e %e %e\r\n", i, dt, t - cfdbv.timecontact, t, simulation_time_1file, drop_time_1file, bubble_time_1file, writefile_time_1file, cellnumber, simulation_time_total, drop_time_total, bubble_time_total, writefile_time_total);
	fclose(fp);
	simulation_str_time = clock();
	simulation_time_1file = 0.0;
	writefile_time_1file = 0.0;
	drop_time_1file = 0.0;
	bubble_time_1file = 0.0;
	;
	switch (pid())
	{
	case 0:
	{
		printf("\r\nData Files are Written!\r\nDuration Last: %s\r\nTotal: %s, Time Left: %s\r\n\r\n", LDc, TDc, ETLc);
		break;
	}
	}
}

static void remove_droplets(scalar c, bool droplet, double dropsize)
{
	int j, n;
	scalar m[];
	const double THR = R_VOFLIMIT; //THRESHOLD
	const double delta = cfdbv.domainsize / pow(2.0, LEVELmax);
	double realD2;
	;
	foreach ()
		m[] = (droplet ? (c[] > THR) : (c[] < (1. - THR)));
	n = tag(m);
	double v[n];
	for (j = 0; j < n; j++)
		v[j] = 0.0;
	foreach_leaf()
	{
		if (m[] > 0)
		{
			j = m[] - 1;
#if dimension == 3
			v[j] += Delta * Delta * Delta * (droplet ? c[] : 1. - c[]); // c[]
#else
			v[j] += Delta * Delta * (droplet ? c[] : 1. - c[]); // c[]
#endif
			//v[j] += dv()*c[];
		}
	}
#if _MPI
	MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	/*	FILE *fp;
	char name[500];
	sprintf(name, "vj-CPU[%02d].txt", pid());
	fp = fopen (name, "a");
	fprintf(fp, "%d, %d, %e\r\n", droplet, n, minvol);
	for (j = 0; j < n; j++)
		fprintf(fp, "%e, ", v[j]);
	fclose(fp);
*/;
	foreach ()
	{
		if (m[] > 0)
		{
			j = m[] - 1;
			realD2 = 4.0 * v[j] / R_PI;
			if (realD2 < dropsize * dropsize * delta * delta)
				c[] = droplet ? 0. : 1.;
		}
	}
	boundary({c});
}

int refinebetweentwofluid(const double xmin, const double xmax, const double ymin, const double ymax, scalar f1, scalar f2, const int lvl)
{
	scalar f1tmp[], f2tmp[];
	const double threshold = 1.0e-3;
	foreach ()
	{
		f1tmp[] = 0.0;
		f2tmp[] = 0.0;
		foreach_dimension()
		{
			if (fabs(f1[2] - f1[-2]) > threshold)
			{
				f1tmp[] = 1.0;
			}
			if (fabs(f2[2] - f2[-2]) > threshold)
			{
				f2tmp[] = 1.0;
			}
		}
	}
	refine(fabs(f1tmp[] - 0.0) < threshold && fabs(f2tmp[] - 1.0) < threshold && level < lvl);
	//	refine( fabs(gtmp[] - 1.0) < 1.0e-3 && x > xmin && x < xmax && y > ymin && y < ymax && level < lvl );
	return 1;
}
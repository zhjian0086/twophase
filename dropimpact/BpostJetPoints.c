#include "constants.h"

#include "jetrecognition.h"
double jetroot_time_1file, jetroot_time_total;

// ./DLpost tb0.05 ts0.01 te0.06
// means:
// time-begin: 0.05
// time-step: 0.01
// time-end: 0.06

#if dimension == 3
#include "lambda2.h"
#endif

double total_time_1file, total_time_total;
clock_t simulation_str_time, simulation_end_time;

double Time_BGN = 0.0, Time_STP = 0.1, Time_END = 0.2;

void readfromargPOST(int argc, char **argv);

int main(int argc, char **argv)
{
	simulation_str_time = clock();
	jetroot_time_1file = 0.0;
	jetroot_bin_time_total = 0.0;
	total_time_1file = 0.0;
	total_time_total = 0.0;
	readfromargPOST(argc, argv);
	if (argc < 4)
	{
		printf("Run like this:\r\n./DLpost 0.0 0.001 1.0\r\nStart Time, Step Size, End Time\r\n");
		return 1;
	}
	printf("Time_BGN-%f__Time_STP-%f__Time_END-%f\r\n", Time_BGN, Time_STP, Time_END);
	;
	size(DOMAIN_WIDTH);
#if AXI
	;
#else
	origin(0, -DOMAIN_WIDTH / 2., -DOMAIN_WIDTH / 2.);
#endif
	run();
	;
	return 1;
}

void readfromargPOST(int argc, char **argv)
{
	int i, j;
	char tmp[100];
	if (argc < 2)
		return;
	for (i = 1; i < argc; i++)
	{
		switch (argv[i][0])
		{
		case 't':
		case 'T':
		{
			switch (argv[i][1])
			{
			case 'b':
			case 'B':
			{
				for (j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				Time_BGN = atof(tmp);
				break;
			}
			case 's':
			case 'S':
			{
				for (j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				Time_STP = atof(tmp);
				break;
			}
			case 'e':
			case 'E':
			{
				for (j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				Time_END = atof(tmp);
				break;
			}
			}
			break;
		}
		}
	}
}

event defaults(i = 0)
{
	interfaces = list_add(NULL, f);
	interfaces = list_add(interfaces, fdrop);
}

event loadfiles(i = 0)
{
	int iloop;
	const double iloopmax = (Time_END - Time_BGN) / Time_STP + 1.0001;
	char nameloadfile[500];
	char ETL[500], TMThis1[500], TMTotal[500];
	double tc, estimatetimeleft;
	int cellnumber;
	scalar varLeft[], varRight[];
	scalar alpha[], kappa[], omega[];
	vector vn[];
	clock_t timestr, timeend, timestrtotal, timeendtotal;
	;
	char JRThis1[500], JRTotal[500];
	double **JetRoots;
	JetRoots = (int **)calloc(5, sizeof(int *));
	for (iloop = 0; iloop < 5; iloop++)
		JetRoots[iloop] = (int *)calloc(2, sizeof(int));
	;
	timestrtotal = clock();
	for (tc = Time_BGN, iloop = 0; tc <= Time_END + Time_STP * 0.01; tc += Time_STP, iloop++)
	{
		// loading
		;
		sprintf(nameloadfile, "DataALL-%.4f", tc);
		restore(file = nameloadfile);
		;
		// calculating
		;
		vorticity(u, omega);
		curvature(f, kappa);
		reconstruction(f, vn, alpha);
		cellnumber = 0;
		foreach ()
			cellnumber++;
		printf("==========----------==========----------==========\r\n");
		printf("time: %.4f ****** cell number: %d\r\n", tc, cellnumber);
		printf("==========----------==========----------==========\r\n");
		;
		// Jet recognition
		;
		timestr = clock();
		char nameTecplotJT[500];
		sprintf(nameTecplotJT, "%s-%.4f.plt", NAMETECPLOTJT, tc);
		FindJetRoots(iloop, tc, kappa, omega, 'n', JetRoots);
		output_tecplot2D_jetcoordinates(nameTecplotJT, tc, JetRoots, angle);
		timeend = clock();
		jetroot_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
		;
		// calculating total time
		;
		jetroot_time_total += jetroot_time_1file;
		timecalculation(jetroot_time_1file, JRThis1);
		timecalculation(jetroot_time_total, JRTotal);
		;
		timeendtotal = clock();
		total_time_1file += (double)(timeendtotal - timestrtotal) / CLOCKS_PER_SEC;
		total_time_total += total_time_1file;
		;
		estimatetimeleft = total_time_total * iloopmax / (iloop + 1.0) - total_time_total;
		timecalculation(estimatetimeleft, ETL);
		timecalculation(total_time_1file, TMThis1);
		timecalculation(total_time_total, TMTotal);
		;
		// finish post-processing
		;
		printf("==========----------==========----------==========\r\n");
		printf("LAST FILE DURATIONS:\r\n");
		printf("jet root duration: %s\r\n", JRThis1);
		printf("total duration: %s\r\n", TMThis1);
		printf("ALL FILES DURATIONS UNTIL NOW:\r\n");
		printf("jet root total duration: %s\r\n", JRThis1);
		jetroot_time_1file = 0.0;
		printf("total duration: %s\r\n", TMTotal);
		total_time_1file = 0.0;
		printf("==========----------==========----------==========\r\n");
		printf("Estimated time left: %s\r\n", ETL);
		printf("done with t = %.4f\r\n", tc);
		printf("==========----------==========----------==========\r\n");
	}
	for (iloop = 0; iloop < 5; iloop++)
		free(JetRoots[iloop]);
	free(JetRoots);
}

event end(i = 0)
{
}

#include "constants.h"

#define TECPLOTFOLDER		"tecplotfiles"
#define NAMEMACRO       	"macro.mcr"
double tecplot_time_1file, tecplot_time_total;

#include "tecplot.h"

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
	tecplot_time_1file = 0.0;
	tecplot_time_total = 0.0;
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
	char TPThis1[500], TPTotal[500];
	const double angle = 90.0;
	char onlyf[500], nameTecplotND[500], nameTecplotCC[500], nameTecplotF1[500], nameTecplotF2[500];
	char namesTecplot[10][500], nameTecplotBinND[500], nameTecplotBinCC[500];
	strcpy(TPThis1, "mkdir ");
	strcat(TPThis1, TECPLOTFOLDER);
	system(TPThis1);
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
		// tecplot output
		;
		timestr = clock();
		printf("==========----------==========----------==========\r\n");
		printf("time: %.4f ****** tecplot interface (f1, f2).\r\n");
		sprintf(nameTecplotF1, "%s/%s-%.4f.plt", TECPLOTFOLDER, NAMETECPLOTF1, tc);
		sprintf(nameTecplotF2, "%s/%s-%.4f.plt", TECPLOTFOLDER, NAMETECPLOTF2, tc);
		output_tecplot2D_Intrfc(nameTecplotF1, f, tc, angle);
		output_tecplot2D_Intrfc(nameTecplotF2, fdrop, tc, angle);
		printf("done!\r\n");
		printf("==========----------==========----------==========\r\n");
		;
		printf("==========----------==========----------==========\r\n");
		printf("time: %.4f ****** tecplot nodal (iso-lines-enable).\r\n");
		sprintf(nameTecplotND, "%s/%s-%.4f.plt", TECPLOTFOLDER, NAMETECPLOTND, tc);
		sprintf(nameTecplotBinND, "%s/%s-%.4f.plt", TECPLOTFOLDER, NAMETECPLOTBINND, tc);
		onlyf[0] = 'y';
		onlyf[1] = 'y';
		onlyf[2] = 'y';
		onlyf[3] = 'y';
		onlyf[4] = 'y';
		output_tecplot2D_nodal(nameTecplotND, tc, {fdrop, u.x, u.y, omega, pressure}, onlyf, angle);
		printf("done!\r\n");
		printf("==========----------==========----------==========\r\n");
		;
		printf("==========----------==========----------==========\r\n");
		printf("time: %.4f ****** tecplot cell center (numerical values).\r\n");
		sprintf(nameTecplotCC, "%s/%s-%.4f.plt", TECPLOTFOLDER, NAMETECPLOTCC, tc);
		sprintf(nameTecplotBinCC, "%s/%s-%.4f.plt", TECPLOTFOLDER, NAMETECPLOTBINCC, tc);
		onlyf[0] = 'y';
		onlyf[1] = 'y';
		onlyf[2] = 'y';
		onlyf[3] = 'y';
		onlyf[4] = 'y';
		output_tecplot2D_cellcenter(nameTecplotCC, tc, {fdrop, u.x, u.y, omega, pressure}, onlyf, angle);
		printf("done!\r\n");
		printf("==========----------==========----------==========\r\n");
		;
		printf("==========----------==========----------==========\r\n");
		printf("time: %.4f ****** tecplot ND, ASCII to BINARY, one file -> one file\r\n");
		strcpy(namesTecplot[0], nameTecplotND);
		strcpy(namesTecplot[1], nameTecplotF1);
		strcpy(namesTecplot[2], nameTecplotF2);
		macrotecplotonefile(namesTecplot, 3, nameTecplotBinND);
		printf("done!\r\n");
		printf("==========----------==========----------==========\r\n");
		;
		printf("==========----------==========----------==========\r\n");
		printf("time: %.4f ****** tecplot CC, ASCII to BINARY, one file -> one file\r\n");
		strcpy(namesTecplot[0], nameTecplotCC);
		strcpy(namesTecplot[1], nameTecplotF1);
		strcpy(namesTecplot[2], nameTecplotF2);
		macrotecplotonefile(namesTecplot, 3, nameTecplotBinCC);
		printf("done!\r\n");
		printf("==========----------==========----------==========\r\n");
		;
		timeend = clock();
		tecplot_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
		;
		// calculating total time
		;
		tecplot_time_total += tecplot_time_1file;
		timecalculation(tecplot_time_1file, TPThis1);
		timecalculation(tecplot_time_total, TPTotal);
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
		printf("tecplot duration: %s\r\n", TPThis1);
		printf("total duration: %s\r\n", TMThis1);
		printf("ALL FILES DURATIONS UNTIL NOW:\r\n");
		printf("tecplot total duration: %s\r\n", TPThis1);
		tecplot_time_1file = 0.0;
		printf("total duration: %s\r\n", TMTotal);
		total_time_1file = 0.0;
		printf("==========----------==========----------==========\r\n");
		printf("Estimated time left: %s\r\n", ETL);
		printf("done with t = %.4f\r\n", tc);
		printf("==========----------==========----------==========\r\n");
	}
    printf("==========----------==========----------==========\r\n");
    printf("remove all the ascii\r\n");
	strcpy(TPThis1, "rm ");
	strcat(TPThis1, TECPLOTFOLDER);
	strcat(TPThis1, "/out-tecplot*");
	system(TPThis1);
    printf("done.\r\n");
    printf("==========----------==========----------==========\r\n");
// 	printf("==========----------==========----------==========\r\n");
//     printf("convert all binary files to one binary file.\r\nTime_BGN-%f__Time_STP-%f__Time_END-%f\r\n", Time_BGN, Time_STP, Time_END);
// 	macrotecplotallfiles(NAMETECPLOTBINND, Time_BGN, Time_STP, Time_END);
// 	macrotecplotallfiles(NAMETECPLOTBINCC, Time_BGN, Time_STP, Time_END);
// 	printf("done!\r\n");
// 	printf("==========----------==========----------==========\r\n");
}

event end(i = 0)
{
}

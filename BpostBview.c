#include "constants.h"

#define BVIEW_YESNO		'y'
#define TECPLOT_YESNO	'n'
#define JETROOT_YESNO	'n'

#define BVIEWFOLDER		"bviewfiles"
double bview_time_1file, bview_time_total;
double SCALE = 1.0, YPOSITION = 0.50, NOPIXELS = 1000;

#include "view.h"

// install bview: http://basilisk.fr/src/gl/INSTALL
// qcc -Wall -O2 DLpost.c -o DLpost -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
// qcc -Wall -O2 DLpost.c -o DLpost -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11

// ./DLpost tb0.05 ts0.01 te0.06 s20.0 y0.6 p1000
// means:
// time-begin: 0.05
// time-step: 0.01
// time-end: 0.06
// scale: 20.0
// y position from the bottom 0.6 (60 %)
// number of pixels in y-direction: 1000 (in the x direction is twice)

// ffmpeg -framerate 1 -i out-bview-omega-%06d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out-bview-omega-video.mp4

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
	bview_time_1file = 0.0;
	bview_time_total = 0.0;
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
		case 's':
		case 'S':
		{
			for (j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			SCALE = atof(tmp);
			break;
		}
		case 'y':
		case 'Y':
		{
			for (j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			YPOSITION = atof(tmp);
			break;
		}
		case 'p':
		case 'P':
		{
			for (j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			NOPIXELS = (int)atof(tmp);
			break;
		}
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
	char BVThis1[500], BVTotal[500];
	strcpy(BVThis1, "mkdir ");
	strcat(BVThis1, BVIEWFOLDER);
	system(BVThis1);
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
		// bview output
		;
		timestr = clock();
		printf("==========----------==========----------==========\r\n");
		printf("bview image producing.\r\n");
		char nameBview[500], textBview[500];
		;
		foreach ()
		{
			varLeft[] = level;
			varRight[] = omega[] * f[];
		}
		boundary({varLeft});
		boundary({varRight});
		sprintf(nameBview, "%s/out-bview-level-vorticity-%09d.png", BVIEWFOLDER, (int)(round(fabs(tc) * 1000000)));
		view(width = 2 * NOPIXELS, height = NOPIXELS, quat = {0, 0, -0.707, 0.707}, sx = SCALE, sy = SCALE, ty = -1.0 * SCALE * YPOSITION);
		clear();
		draw_vof("f", lw = 5);
		draw_vof("fdrop", lw = 5);
		squares("varLeft");
		box(notics = true);
		cells(lw = 1); // showing the grid
		mirror({-1})
		{
			draw_vof("f", lw = 5);
			draw_vof("fdrop", lw = 5);
			squares("varRight", linear = true, min = -10.0, max = 10.0); //
			box(notics = true);
		}
		sprintf(textBview, "L: grid, R: vorticity, time: %.4f", tc);
		draw_string(textBview, lw = 5, size = 100);
		save(nameBview);
		;
		foreach ()
		{
			varLeft[] = 1.0 - f[] + 2.0 * fdrop[] < 0.0 ? 0.0 : (1.0 - f[] + 2.0 * fdrop[] > 2.0 ? 2.0 : 1.0 - f[] + 2.0 * fdrop[]);
			//varLeft[] = p[] * f[];
			varRight[] = pressure[] * f[];
		}
		boundary({varLeft});
		boundary({varRight});
		sprintf(nameBview, "%s/out-bview-fraction-pressure-%09d.png", BVIEWFOLDER, (int)(round(fabs(tc) * 1000000)));
		view(width = 2 * NOPIXELS, height = NOPIXELS, quat = {0, 0, -0.707, 0.707}, sx = SCALE, sy = SCALE, ty = -1.0 * SCALE * YPOSITION);
		clear();
		draw_vof("f", lw = 5);
		draw_vof("fdrop", lw = 5);
		squares("varLeft", linear = true, min = 0.0, max = 2.0);
		box(notics = true);
//		cells(lw = 1); // showing the grid
		mirror({-1})
		{
			draw_vof("f", lw = 5);
			draw_vof("fdrop", lw = 5);
			squares("varRight", linear = true, min = 0.0, max = 5.0); //
			box(notics = true);
		}
		sprintf(textBview, "L: fraction, R: pressure, time: %.4f", tc);
		draw_string(textBview, lw = 5, size = 100);
		save(nameBview);
		;
		printf("done!\r\n");
		printf("==========----------==========----------==========\r\n");
		timeend = clock();
		bview_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
		;
		// calculating total time
		;
		bview_time_total += bview_time_1file;
		timecalculation(bview_time_1file, BVThis1);
		timecalculation(bview_time_total, BVTotal);
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
		printf("bview duration: %s\r\n", BVThis1);
		printf("total duration: %s\r\n", TMThis1);
		printf("ALL FILES DURATIONS UNTIL NOW:\r\n");
		printf("bview total duration: %s\r\n", BVThis1);
		bview_time_1file = 0.0;
		printf("total duration: %s\r\n", TMTotal);
		total_time_1file = 0.0;
		printf("==========----------==========----------==========\r\n");
		printf("Estimated time left: %s\r\n", ETL);
		printf("done with t = %.4f\r\n", tc);
		printf("==========----------==========----------==========\r\n");
	}
}

event end(i = 0)
{
}

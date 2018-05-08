//#include "axi.h"
//#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "curvature.h"

#define DIM_NONDIM_EXP			'n' // d: dimension; n: nondimension; e: experimentalization

#if DIM_NONDIM_EXP == 'd' || DIM_NONDIM_EXP == 'D'

#define BUBBLE_DIAMETER			2.0e-3
#define RHO_L				851.0
#define RHO_G				1.2
#define MU_L				1.3e-3
#define MU_G				1.94e-5
#define SIGMA				17.6e-3
#define GRAVITY				9.806
//
#define RHO_GL				0.0
#define MU_GL				0.0
#define GALILEI				0.0
#define EOTVOS				0.0

#elif DIM_NONDIM_EXP == 'n' || DIM_NONDIM_EXP == 'N'

#define EOTVOS				10.0
#define GALILEI				10.0
#define RHO_GL				0.0010
#define MU_GL				0.010
#define RVERTICAL_RHORIZONTAL_RATIO	0.70
//
#define BUBBLE_DIAMETER			0.0
#define RHO_L				0.0
#define MU_L				0.0
#define SIGMA				0.0
#define RHO_G				0.0
#define MU_G				0.0
#define GRAVITY				0.0

#elif DIM_NONDIM_EXP == 'e' || DIM_NONDIM_EXP == 'E'

#define EOTVOS				300.0
#define GALILEI				1000.0
#define BUBBLE_DIAMETER			2.0e-3
#define SIGMA				17.6e-3
#define RHO_L				816.0
#define RHO_G				1.2041
#define MU_G				1.94e-5
//
#define VELOCITY			0.0
#define RHO_GL				0.0
#define MU_GL				0.0
#define MU_L				0.0

#endif

#define INITAL_GRID_LEVEL		6
#define MAX_GRID_LEVEL			9
#define DOMAIN_WIDTH			20.00
#define INITIAL_DISTANCE		0.05
#define REFINE_GAP			0.02
#define MAX_TIME			30.00

#define REFINE_VAR			{f, u.x, u.y, omega}
#define REFINE_VAR_TEXT			"f, u.x, u.y, omega"
#define REFINE_VALUE_0			-6
#define REFINE_VALUE_1			-1
#define REFINE_VALUE_2			-1
#define REFINE_VALUE_3			-4

#define R_VOFLIMIT			1.0e-6

#define REMOVE_DROP_YESNO		'n'
#define REMOVE_DROP_SIZE		10.0 // equivalent radius base on the maximum refinement
#define REMOVE_DROP_PERIOD		10
#define REMOVE_BUBBLE_YESNO		'n'
#define REMOVE_BUBBLE_SIZE		5.0 // equivalent radius base on the maximum refinement
#define REMOVE_BUBBLE_PERIOD		10

#define SAVE_FILE_EVERY			0.1 // critical time
#define SHOWITERATION			1

#define FILENAME_DATA			"data"
#define FILENAME_DURATION		"duration"
#define FILENAME_PARAMETERS		"parameters.txt"
#define FILENAME_ENDOFRUN		"endofrun"
#define FILENAME_LASTFILE		"lastfile"

#define R_PI					3.1415926535897932384626433832795

int LEVELmin = INITAL_GRID_LEVEL, LEVELmax = MAX_GRID_LEVEL;
double maxruntime = HUGE;
scalar fdrop[], pressure[];

struct CFDValues {
	double rhoL, rhoG, muL, muG, Sigma, gravity;
	double Galilei, Eotvos;
	double radius, domainsize, refinegap, initialdis;
};

void readfromarg(char **argv, int argc, struct CFDValues *bvalues);

int numericalmainvalues(char **argv, int argc, struct CFDValues *bvalues)
{
	bvalues->rhoL = 1.0;
	bvalues->radius = 1.0;
	bvalues->gravity = 1.0;
	;
	bvalues->Galilei = -1.0;
	bvalues->Eotvos = -1.0;
	readfromarg(argv, argc, bvalues);
	switch(DIM_NONDIM_EXP)
	{
	case 'd':
	case 'D':
	{
		break;
	}
	case 'n':
	case 'N':
	{
		if (bvalues->Galilei < 0.0)
			bvalues->Galilei = GALILEI;
		if (bvalues->Eotvos < 0.0)
			bvalues->Eotvos = EOTVOS;
		bvalues->muL = (bvalues->rhoL * sqrt(bvalues->gravity * bvalues->radius) * bvalues->radius / bvalues->Galilei);
		bvalues->Sigma = (bvalues->rhoL * bvalues->gravity * bvalues->radius * bvalues->radius / bvalues->Eotvos);
		bvalues->rhoG = RHO_GL * bvalues->rhoL;
		bvalues->muG = MU_GL * bvalues->muL;
		break;
	}
	case 'e':
	case 'E':
	{
		break;
	}
	}
	bvalues->domainsize = DOMAIN_WIDTH * bvalues->radius;
	bvalues->initialdis = INITIAL_DISTANCE * bvalues->radius;
	bvalues->refinegap = REFINE_GAP * bvalues->radius;
	;
	switch (pid())
	{
	case 0:
	{
		printf("Galilei: %f --- Eotvos: %f\r\n", bvalues->Galilei, bvalues->Eotvos);
		FILE *fp;
		fp = fopen (FILENAME_PARAMETERS, "w");
		fprintf (fp, "Physical / Numerical Values\r\n");
		fprintf (fp, "Radius: %.3e / %.3e\r\n", BUBBLE_DIAMETER, bvalues->radius);
		fprintf (fp, "Rho(L): %.3e / %.3e\r\n", RHO_L, bvalues->rhoL);
		fprintf (fp, "Rho(G): %.3e / %.3e\r\n", RHO_G, bvalues->rhoG);
		fprintf (fp, "Mu(L): %.3e / %.3e\r\n", MU_L, bvalues->muL);
		fprintf (fp, "MU(G): %.3e / %.3e\r\n", MU_G, bvalues->muG);
		fprintf (fp, "Sigma: %.3e / %.3e\r\n", SIGMA, bvalues->Sigma);
		fprintf (fp, "Galilei: %.2f\r\n", bvalues->Galilei);
		fprintf (fp, "Eotvos: %.2f\r\n", bvalues->Eotvos);
		fprintf (fp, "Radius Horizontal: %.3e\r\n", bvalues->radius / sqrt (RVERTICAL_RHORIZONTAL_RATIO));
		fprintf (fp, "Radius Vertical: %.3e\r\n", bvalues->radius * sqrt (RVERTICAL_RHORIZONTAL_RATIO));
		fprintf (fp, "\r\n");
		fprintf (fp, "Level Max: %d\r\n", LEVELmax);
		fprintf (fp, "Level Min: %d\r\n", LEVELmin);
		fprintf (fp, "Domain Size: %.2f\r\n", bvalues->domainsize);
		fprintf (fp, "Initial Distance: %.2f\r\n", bvalues->initialdis);
		fprintf (fp, "Refine Gap: %.2f\r\n", bvalues->refinegap);
		fclose (fp);
		break;
	}
	}
	return 1;
}

void readfromarg(char **argv, int argc, struct CFDValues *bvalues)
{
	int i, j;
	char tmp[100];
	if (argc < 2)
		return;
	for (i = 1; i < argc; i++)
	{
		switch(argv[i][0])
		{
		case 'g':
		case 'G':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Galilei = atof(tmp);
			break;
		}
		case 'e':
		case 'E':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Eotvos = atof(tmp);
			break;
		}
		case 'x':
		case 'X':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			LEVELmax = atoi(tmp);
			break;
		}
		case 'n':
		case 'N':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			LEVELmin = atoi(tmp);
			break;
		}
		}
	}
}

int timecalculation(double t, char *chartime)
{
	int d, h, m, s;
	if(t < 60.0)
	{
		d = 0;
		h = 0;
		m = 0;
		s = (int) t;
	}
	else if(t < 3600.0)
	{
		d = 0;
		h = 0;
		m = (int) (t / 60.0);
		s = (int) (t - m*60.0);
	}
	else if(t < 3600.0*24.0)
	{
		d = 0;
		h = (int) (t / 3600.0);
		m = (int) ((t - h*3600.0) / 60.0);
		s = (int) (t - h*3600.0 - m*60.0);
	}
	else
	{
		d = (int) (t / 3600.0 / 24.0);
		h = (int) ((t - d*3600.0*24.0) / 3600.0);
		m = (int) ((t - d*3600.0*24.0 - h*3600.0) / 60.0);
		s = (int) (t - d*3600.0*24.0 - h*3600.0 - m*60.0);
	}
	sprintf(chartime, "%d:%02d:%02d:%02d", d, h, m, s);
	return 1;
}


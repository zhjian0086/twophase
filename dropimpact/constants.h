#include "axi.h"
//#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "curvature.h"

#define DIM_NONDIM_EXP			'n' // d: dimension; n: nondimension; e: experimentalization

#if DIM_NONDIM_EXP == 'd' || DIM_NONDIM_EXP == 'D'

#define VELOCITY			1.56
#define DROP_DIAMETER		2.0e-3
#define RHO_L				851.0
#define RHO_G				1.2
#define MU_L				1.3e-3
#define MU_G				1.94e-5
#define SIGMA				17.6e-3
#define GRAVITY				9.81
//
#define RHO_GL				0.0
#define MU_GL				0.0
#define REYNOLDS			0.0
#define WEBER				0.0

#elif DIM_NONDIM_EXP == 'n' || DIM_NONDIM_EXP == 'N'

#define WEBER				100.0
#define REYNOLDS			100.0
#define RHO_GL				(0.0012) // air-water at 25C: 0.001187503
#define MU_GL				(0.0210) // air-water at 25C: 0.020898876
//
#define VELOCITY			0.0
#define DROP_DIAMETER		0.0
#define RHO_L				0.0
#define MU_L				0.0
#define SIGMA				0.0
#define RHO_G				0.0
#define MU_G				0.0
#define GRAVITY				0.0

#elif DIM_NONDIM_EXP == 'e' || DIM_NONDIM_EXP == 'E'

#define WEBER				300.0
#define REYNOLDS			1000.0
#define DROP_DIAMETER		2.0e-3
#define SIGMA				17.6e-3
#define RHO_L				816.0
#define RHO_G				1.2041
#define MU_G				1.94e-5
//
#define VELOCITY			0.0
#define RHO_GL				0.0
#define MU_GL				0.0
#define MU_L				0.0
#define GRAVITY				0.0

#endif

#define INITAL_GRID_LEVEL		6
#define MAX_GRID_LEVEL			9
#define DOMAIN_WIDTH			5.00
#define POOL_DEPTH				0.10
#define INITIAL_DISTANCE		0.05
#define REFINE_GAP				0.01
#define MAX_TIME				0.20
#define SAVE_FILE_EVERY			0.01 // critical time

#define REFINE_VAR				{f, u.x, u.y, fdrop}
#define REFINE_VAR_TEXT			"f, u.x, u.y, fdrop"
#define REFINE_VALUE_0			-6
#define REFINE_VALUE_1			-3
#define REFINE_VALUE_2			-3
#define REFINE_VALUE_3			-6

#define REMOVE_DROP_YESNO		'y'
#define REMOVE_DROP_SIZE		10.0 // equivalent diameter base on the maximum refinement
#define REMOVE_DROP_PERIOD		10
#define REMOVE_BUBBLE_YESNO		'y'
#define REMOVE_BUBBLE_SIZE		10.0 // equivalent diameter base on the maximum refinement
#define REMOVE_BUBBLE_PERIOD	10


#define FILENAME_DATA			"data"
#define FILENAME_DURATION		"duration"
#define FILENAME_PARAMETERS		"parameters.txt"
#define FILENAME_ENDOFRUN		"endofrun"
#define FILENAME_LASTFILE		"lastfile"

#define R_VOFLIMIT				1.0e-9
#define R_PI					3.1415926535897932384626433832795

int LEVELmin = INITAL_GRID_LEVEL, LEVELmax = MAX_GRID_LEVEL;
double maxruntime = HUGE;
scalar fdrop[], pressure[];

struct CFDValues {
	double rhoL, rhoG, muL, muG, Sigma;
	double vel, Reynolds, Weber;
	double diameter, domainsize, refinegap, pooldepth, initialdis;
	double timecontact, timeend, timestep;
};

void readfromarg(char **argv, int argc, struct CFDValues *bvalues);

int numericalmainvalues(char **argv, int argc, struct CFDValues *bvalues)
{
	double velocity = VELOCITY, mu_l = MU_L;
	bvalues->rhoL = 1.0;
	bvalues->vel = 1.0;
	bvalues->diameter = 1.0;
	;
	bvalues->Reynolds = -1.0;
	bvalues->Weber = -1.0;
	bvalues->pooldepth = -1.0;
	bvalues->timeend = -1.0;
	bvalues->timestep = -1.0;
	;
	readfromarg(argv, argc, bvalues);
	switch(DIM_NONDIM_EXP)
	{
	case 'd':
	case 'D':
	{
		bvalues->Reynolds = (RHO_L * VELOCITY * DROP_DIAMETER / mu_l);
		bvalues->Weber = (RHO_L * VELOCITY * VELOCITY * DROP_DIAMETER / SIGMA);
		bvalues->muL = (bvalues->rhoL * bvalues->vel * bvalues->diameter / bvalues->Reynolds);
		bvalues->Sigma = (bvalues->rhoL * bvalues->vel * bvalues->vel * bvalues->diameter / bvalues->Weber);
		bvalues->rhoG = (RHO_G / RHO_L) * bvalues->rhoL;
		bvalues->muG = (MU_G / mu_l) * bvalues->muL;
		break;
	}
	case 'n':
	case 'N':
	{
		if (bvalues->Reynolds < 0.0)
			bvalues->Reynolds = REYNOLDS;
		if (bvalues->Weber < 0.0)
			bvalues->Weber = WEBER;
		bvalues->muL = (bvalues->rhoL * bvalues->vel * bvalues->diameter / bvalues->Reynolds);
		bvalues->Sigma = (bvalues->rhoL * bvalues->vel * bvalues->vel * bvalues->diameter / bvalues->Weber);
		bvalues->rhoG = RHO_GL * bvalues->rhoL;
		bvalues->muG = MU_GL * bvalues->muL;
		break;
	}
	case 'e':
	case 'E':
	{
		if (bvalues->Reynolds < 0.0)
			bvalues->Reynolds = REYNOLDS;
		if (bvalues->Weber < 0.0)
			bvalues->Weber = WEBER;
		bvalues->muL = (bvalues->rhoL * bvalues->vel * bvalues->diameter / bvalues->Reynolds);
		bvalues->Sigma = (bvalues->rhoL * bvalues->vel * bvalues->vel * bvalues->diameter / bvalues->Weber);
		velocity = sqrt (bvalues->Weber * SIGMA / (DROP_DIAMETER * RHO_L));
		mu_l = (RHO_L * velocity * DROP_DIAMETER / bvalues->Reynolds);
		bvalues->muG = (MU_G / mu_l) * bvalues->muL;
		bvalues->rhoG = (RHO_G / RHO_L) * bvalues->rhoL;
		break;
	}
	}
	bvalues->domainsize = DOMAIN_WIDTH * bvalues->diameter;
	if (bvalues->pooldepth < 0.0)
		bvalues->pooldepth = POOL_DEPTH * bvalues->diameter;
	bvalues->initialdis = INITIAL_DISTANCE * bvalues->diameter;
	bvalues->refinegap = REFINE_GAP * bvalues->diameter;
	;
	bvalues->timecontact = bvalues->initialdis / bvalues->vel;
	if (bvalues->timeend < 0.0)
		bvalues->timeend = MAX_TIME;
	if (bvalues->timestep < 0.0)
		bvalues->timestep = SAVE_FILE_EVERY;
	;
	switch (pid())
	{
	case 0:
	{
		printf("R: %f --- W: %f --- H: %f\r\n", bvalues->Reynolds, bvalues->Weber, bvalues->pooldepth);
		FILE *fp;
		fp = fopen (FILENAME_PARAMETERS, "w");
		fprintf (fp, "Physical / Numerical Values\r\n");
		fprintf (fp, "Diameter: %.3e / %.3e\r\n", DROP_DIAMETER, bvalues->diameter);
		fprintf (fp, "Velocity: %.3e / %.3e\r\n", velocity, bvalues->vel);
		fprintf (fp, "Rho(L): %.3e / %.3e\r\n", RHO_L, bvalues->rhoL);
		fprintf (fp, "Rho(G): %.3e / %.3e\r\n", RHO_G, bvalues->rhoG);
		fprintf (fp, "Mu(L): %.3e / %.3e\r\n", mu_l, bvalues->muL);
		fprintf (fp, "MU(G): %.3e / %.3e\r\n", MU_G, bvalues->muG);
		fprintf (fp, "Sigma: %.3e / %.3e\r\n", SIGMA, bvalues->Sigma);
		fprintf (fp, "Reynolds: %.2f\r\n", bvalues->Reynolds);
		fprintf (fp, "Weber: %.2f\r\n", bvalues->Weber);
		fprintf (fp, "\r\n");
		fprintf (fp, "Level Max: %d\r\n", LEVELmax);
		fprintf (fp, "Level Min: %d\r\n", LEVELmin);
		fprintf (fp, "Domain Size: %.2f\r\n", bvalues->domainsize);
		fprintf (fp, "Pool Depth: %.2f\r\n", bvalues->pooldepth);
		fprintf (fp, "Initial Distance: %.2f\r\n", bvalues->initialdis);
		fprintf (fp, "Refine Gap: %.2f\r\n", bvalues->refinegap);
		fprintf (fp, "Contact Time: %.2f\r\n", bvalues->timecontact);
		fprintf (fp, "Domain Size: %.2f\r\n", bvalues->domainsize);
		fprintf (fp, "\r\n");
		fprintf (fp, "Refine Variables: %s\r\n", REFINE_VAR_TEXT);
		fprintf (fp, "Refine Variables powers: %d, %d, %d, %d\r\n", REFINE_VALUE_0, REFINE_VALUE_1, REFINE_VALUE_2, REFINE_VALUE_3);
		fprintf (fp, "\r\n");
		fprintf (fp, "Remove Drop YesNo: %c\r\n", REMOVE_DROP_YESNO);
		fprintf (fp, "Remove Drop Size: %f\r\n", REMOVE_DROP_SIZE);
		fprintf (fp, "Remove Drop Period: %d\r\n", REMOVE_DROP_PERIOD);
		fprintf (fp, "\r\n");
		fprintf (fp, "Remove Bubble YesNo: %c\r\n", REMOVE_BUBBLE_YESNO);
		fprintf (fp, "Remove Bubble Size: %f\r\n", REMOVE_BUBBLE_SIZE);
		fprintf (fp, "Remove Bubble Period: %d\r\n", REMOVE_BUBBLE_PERIOD);
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
		case 'r':
		case 'R':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Reynolds = atof(tmp);
			break;
		}
		case 'w':
		case 'W':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Weber = atof(tmp);
			break;
		}
		case 'h':
		case 'H':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->pooldepth = atof(tmp);
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
		case 't':
		case 'T':
		{
			switch(argv[i][1])
			{
			case 'e':
			case 'E':
			{
				for(j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				bvalues->timeend = atof(tmp);
				break;
			}
			case 's':
			case 'S':
			{
				for(j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				bvalues->timestep = atof(tmp);
				break;
			}
			}
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


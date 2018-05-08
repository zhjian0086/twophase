
#define NAMETECPLOTJT			"out-tecplot-JT"
#define NAMETECPLOTND			"out-tecplot-ND"
#define NAMETECPLOTCC			"out-tecplot-CC"
#define NAMETECPLOTF1			"out-tecplot-F1"
#define NAMETECPLOTF2			"out-tecplot-F2"
#define NAMETECPLOTBINND		"out-bin-tecplot-ND"
#define NAMETECPLOTBINCC		"out-bin-tecplot-CC"

int system(const char *command);

scalar cal_intersectionarea();
double calculate_output_value(char type, double xxx, double yyy, scalar fract, scalar myvar);
void cellcornerinterpolate(scalar cellvar, double *pointvar);
int findintersectionpoints(double nx, double ny, double alpha, double xc, double yc, double delta, double xi[2], double yi[2], char typei[2]);

int macrotecplotonefile(char names[10][500], int NO, char namebin[500])
{
    int i;
    char command[500];
    FILE *fp;
	;
	fp = fopen(NAMEMACRO, "w");
	fprintf(fp, "#!MC 1410");
	fprintf(fp, "\r\n$!READDATASET ");
	fprintf(fp, "'");
	for (i = 0; i < NO; i++)
		fprintf(fp, " \"%s\" ", names[i]);
	fprintf(fp, "'");
	fprintf(fp, "\r\nREADDATAOPTION = NEW");
	fprintf(fp, "\r\nRESETSTYLE = YES");
	fprintf(fp, "\r\nVARLOADMODE = BYNAME");
	fprintf(fp, "\r\nASSIGNSTRANDIDS = YES");
	fprintf(fp, "\r\nVARNAMELIST = '\"x\" \"y\" \"f\" \"fdrop\" \"u.x\" \"u.y\" \"omega\" \"pressure\"'");
	fprintf(fp, "\r\n$!WRITEDATASET  \"%s\"", namebin);
	fprintf(fp, "\r\nINCLUDETEXT = NO");
	fprintf(fp, "\r\nINCLUDEGEOM = NO");
	fprintf(fp, "\r\nINCLUDEDATASHARELINKAGE = YES");
	fprintf(fp, "\r\nBINARY = YES");
	fprintf(fp, "\r\nUSEPOINTFORMAT = NO");
	fprintf(fp, "\r\nPRECISION = 9");
	fprintf(fp, "\r\nTECPLOTVERSIONTOWRITE = TECPLOTCURRENT");
	fclose(fp);
	;
	strcpy(command, "/usr/local/tecplot360ex/bin/tec360 -b -p ");
	strcat(command, NAMEMACRO);
	system(command);
	;
    return 0;
}

int macrotecplotallfiles(const char name[500], double Time_BGN, double Time_STP, double Time_END)
{
    double tc;
    int iloop, NO;
    char command[500], *filenames;
    FILE *fp;
    ;
    NO = (int) (1.1 * (Time_END - Time_BGN) / Time_STP);
    filenames = (char *)calloc(500 * NO, sizeof(char));
	;
    fp = fopen(NAMEMACRO, "w");
    fprintf(fp, "#!MC 1410\r\n$!READDATASET ");
    fprintf(fp, "'");
    for (tc = Time_BGN, iloop = 0; tc <= Time_END + Time_STP * 0.01; tc += Time_STP, iloop++)
    {
        fprintf(fp, " \"%s-%.4f.plt\" ", name, tc);
    }
    fprintf(fp, "'");
    fprintf(fp, "\r\nREADDATAOPTION = NEW");
    fprintf(fp, "\r\nRESETSTYLE = YES");
    fprintf(fp, "\r\nVARLOADMODE = BYNAME");
    fprintf(fp, "\r\nASSIGNSTRANDIDS = YES");
    fprintf(fp, "\r\nVARNAMELIST = '\"x\" \"y\" \"f\" \"fdrop\" \"u.x\" \"u.y\" \"omega\" \"pressure\"'");
    fprintf(fp, "\r\n$!WRITEDATASET  \"%s.plt\"", name);
    fprintf(fp, "\r\nINCLUDETEXT = NO");
    fprintf(fp, "\r\nINCLUDEGEOM = NO");
    fprintf(fp, "\r\nINCLUDEDATASHARELINKAGE = YES");
    fprintf(fp, "\r\nBINARY = YES");
    fprintf(fp, "\r\nUSEPOINTFORMAT = NO");
    fprintf(fp, "\r\nPRECISION = 9");
    fprintf(fp, "\r\nTECPLOTVERSIONTOWRITE = TECPLOTCURRENT");
    fclose(fp);
    ;
    strcpy(command, "/usr/local/tecplot360ex/bin/tec360 -b -p ");
    strcat(command, NAMEMACRO);
    system(command);
    ;
    free(filenames);
    return 0;
}

void output_tecplot2D_cellcenter(char *name, double time, scalar *outlist, char *onlyf, double rotationangle)
{
	scalar *list = dump_list(outlist);
	;
	int i, j, cellnumber;
	//	int varNO = list_len(list);
	double xxx[4], yyy[4];
	const double angle = rotationangle * R_PI / 180.0;
	FILE *fp;
	scalar vartmp[];
	;
	cellnumber = 0;
	foreach ()
		cellnumber++;
	;
	fp = fopen(name, "w");
	fprintf(fp, "variables = x, y, f");
	i = 0;
	for (scalar s in list)
	{
		if (strcmp(s.name, "cm"))
			fprintf(fp, ", %s", s.name);
	}
	fprintf(fp, "\r\nZONE T=CellCenterData"
				"\r\nN = %d, E = %d, DATAPACKING = BLOCK, ZONETYPE = FEQUADRILATERAL",
			cellnumber * 4, cellnumber);
	fprintf(fp, "\r\nVARLOCATION = (NODAL, NODAL, CELLCENTERED");
	for (scalar s in list)
	{
		if (strcmp(s.name, "cm"))
			fprintf(fp, ", CELLCENTERED");
	}
	fprintf(fp, ")");
	fprintf(fp, "\r\nSOLUTIONTIME = %.3e", time);
	;
	foreach_leaf()
	{
		xxx[0] = x - 0.50 * Delta;
		xxx[1] = x + 0.50 * Delta;
		xxx[2] = x + 0.50 * Delta;
		xxx[3] = x - 0.50 * Delta;
		yyy[0] = y - 0.50 * Delta;
		yyy[1] = y - 0.50 * Delta;
		yyy[2] = y + 0.50 * Delta;
		yyy[3] = y + 0.50 * Delta;
		fprintf(fp, "\r\n%e", xxx[0] * cos(angle) - yyy[0] * sin(angle));
		fprintf(fp, "\r\n%e", xxx[1] * cos(angle) - yyy[1] * sin(angle));
		fprintf(fp, "\r\n%e", xxx[2] * cos(angle) - yyy[2] * sin(angle));
		fprintf(fp, "\r\n%e", xxx[3] * cos(angle) - yyy[3] * sin(angle));
	}
	foreach_leaf()
	{
		xxx[0] = x - 0.50 * Delta;
		xxx[1] = x + 0.50 * Delta;
		xxx[2] = x + 0.50 * Delta;
		xxx[3] = x - 0.50 * Delta;
		yyy[0] = y - 0.50 * Delta;
		yyy[1] = y - 0.50 * Delta;
		yyy[2] = y + 0.50 * Delta;
		yyy[3] = y + 0.50 * Delta;
		fprintf(fp, "\r\n%e", yyy[0] * cos(angle) + xxx[0] * sin(angle));
		fprintf(fp, "\r\n%e", yyy[1] * cos(angle) + xxx[1] * sin(angle));
		fprintf(fp, "\r\n%e", yyy[2] * cos(angle) + xxx[2] * sin(angle));
		fprintf(fp, "\r\n%e", yyy[3] * cos(angle) + xxx[3] * sin(angle));
	}
	foreach_leaf()
		fprintf(fp, "\r\n%e", f[]);
	j = 0;
	for (scalar s in list)
	{
		if (strcmp(s.name, "cm"))
		{
			if (!strcmp(s.name, "fdrop"))
			{
				if (onlyf[j] == 'y')
				{
					foreach ()
					{
						vartmp[] = 1.0 - f[] + 2.0 * fdrop[];
						if (vartmp[] < 0.0)
							vartmp[] = 0.0;
						else if (vartmp[] > 2.0)
							vartmp[] = 2.0;
					}
				}
				else
				{
					foreach ()
						vartmp[] = fdrop[];
				}
			}
			else
			{
				if (!strcmp(s.name, "u.x"))
				{
					foreach ()
						vartmp[] = u.x[] * cos(angle) - u.y[] * sin(angle);
				}
				else if (!strcmp(s.name, "u.y"))
				{
					foreach ()
						vartmp[] = u.y[] * cos(angle) + u.x[] * sin(angle);
				}
				else
				{
					foreach ()
						vartmp[] = s[];
				}
				if (onlyf[j] == 'y')
				{
					foreach ()
						vartmp[] *= f[];
				}
			}
			foreach_leaf()
				fprintf(fp, "\r\n%e", vartmp[]);
			j++;
		}
	}
	i = 0;
	foreach_leaf()
	{
		fprintf(fp, "\r\n%d %d %d %d", i * 4 + 1, i * 4 + 2, i * 4 + 3, i * 4 + 4);
		i++;
	}
	fclose(fp);
	return;
}

void output_tecplot2D_jetcoordinates(char *name, double tc, int **JetRoots, double rotationangle)
{
	int i, k;
	const double angle = rotationangle * R_PI / 180.0;
	double **xjet, **yjet;
	FILE *fp;
	;
	xjet = (double **)calloc(JetRoots[0][0] + 1, sizeof(double *));
	yjet = (double **)calloc(JetRoots[0][0] + 1, sizeof(double *));
	for (i = 0; i <= JetRoots[0][0]; i++)
	{
		xjet[i] = (double *)calloc(2, sizeof(double));
		yjet[i] = (double *)calloc(2, sizeof(double));
	}
	k = 0;
	foreach_leaf()
	{
		for (i = 1; i <= JetRoots[0][0]; i++)
		{
			if (k == JetRoots[i][0])
			{
				xjet[i][0] = x * cos(angle) - y * sin(angle);
				yjet[i][0] = y * cos(angle) + x * sin(angle);
			}
			if (k == JetRoots[i][1])
			{
				xjet[i][1] = x * cos(angle) - y * sin(angle);
				yjet[i][1] = y * cos(angle) + x * sin(angle);
			}
		}
		k++;
	}
	fp = fopen(name, "w");
	for (k = 1; k <= JetRoots[0][0]; k++)
	{
		fprintf(fp,
				"\r\nvariables = x, y"
				"\r\nZONE T=\"JetData%d\""
				"\r\nN = %d, E = %d, DATAPACKING = BLOCK, ZONETYPE = FETRIANGLE"
				"\r\nSOLUTIONTIME = %.4e",
				k, 3, 1, tc);
		fprintf(fp, "\r\n%e %e %e", xjet[k][0], xjet[k][1], xjet[k][0]);
		fprintf(fp, "\r\n%e %e %e", yjet[k][0], yjet[k][1], yjet[k][0]);
		fprintf(fp, "\r\n%d %d %d", 1, 2, 3);
	}
	fclose(fp);
	;
	// free memory
	;
	for (i = 0; i <= JetRoots[0][0]; i++)
	{
		free(xjet[i]);
		free(yjet[i]);
	}
	free(xjet);
	free(yjet);
	;
	return;
}

void output_tecplot2D_nodal(char *name, double time, scalar *outlist, char *onlyf, double rotationangle)
{
	scalar *list = dump_list(outlist);
	;
	int i, j, p, cellnumber;
	//	int varNO = list_len(list);
	double *f1inter, *varinter, *varinter2, xxx[4], yyy[4];
	const double angle = rotationangle * R_PI / 180.0;
	FILE *fp;
	;
	cellnumber = 0;
	foreach ()
		cellnumber++;
	;
	f1inter = (double *)calloc(cellnumber * 4, sizeof(double));
	varinter = (double *)calloc(cellnumber * 4, sizeof(double));
	varinter2 = (double *)calloc(cellnumber * 4, sizeof(double));
	;
	cellcornerinterpolate(f, f1inter);
	;
	fp = fopen(name, "w");
	fprintf(fp, "variables = x, y, f");
	i = 0;
	for (scalar s in list)
	{
		if (strcmp(s.name, "cm"))
			fprintf(fp, ", %s", s.name);
	}
	fprintf(fp, "\r\nZONE T=NodalData"
				"\r\nN = %d, E = %d, DATAPACKING = BLOCK, ZONETYPE = FEQUADRILATERAL",
			cellnumber * 4, cellnumber);
	fprintf(fp, "\r\nSOLUTIONTIME = %.3e", time);
	;
	foreach_leaf()
	{
		xxx[0] = x - 0.50 * Delta;
		xxx[1] = x + 0.50 * Delta;
		xxx[2] = x + 0.50 * Delta;
		xxx[3] = x - 0.50 * Delta;
		yyy[0] = y - 0.50 * Delta;
		yyy[1] = y - 0.50 * Delta;
		yyy[2] = y + 0.50 * Delta;
		yyy[3] = y + 0.50 * Delta;
		fprintf(fp, "\r\n%e", xxx[0] * cos(angle) - yyy[0] * sin(angle));
		fprintf(fp, "\r\n%e", xxx[1] * cos(angle) - yyy[1] * sin(angle));
		fprintf(fp, "\r\n%e", xxx[2] * cos(angle) - yyy[2] * sin(angle));
		fprintf(fp, "\r\n%e", xxx[3] * cos(angle) - yyy[3] * sin(angle));
	}
	foreach_leaf()
	{
		xxx[0] = x - 0.50 * Delta;
		xxx[1] = x + 0.50 * Delta;
		xxx[2] = x + 0.50 * Delta;
		xxx[3] = x - 0.50 * Delta;
		yyy[0] = y - 0.50 * Delta;
		yyy[1] = y - 0.50 * Delta;
		yyy[2] = y + 0.50 * Delta;
		yyy[3] = y + 0.50 * Delta;
		fprintf(fp, "\r\n%e", yyy[0] * cos(angle) + xxx[0] * sin(angle));
		fprintf(fp, "\r\n%e", yyy[1] * cos(angle) + xxx[1] * sin(angle));
		fprintf(fp, "\r\n%e", yyy[2] * cos(angle) + xxx[2] * sin(angle));
		fprintf(fp, "\r\n%e", yyy[3] * cos(angle) + xxx[3] * sin(angle));
	}
	for (i = 0; i < cellnumber; i++)
	{
		fprintf(fp, "\r\n%e", f1inter[0 + i * 4]);
		fprintf(fp, "\r\n%e", f1inter[1 + i * 4]);
		fprintf(fp, "\r\n%e", f1inter[2 + i * 4]);
		fprintf(fp, "\r\n%e", f1inter[3 + i * 4]);
	};
	j = 0;
	for (scalar s in list)
	{
		if (strcmp(s.name, "cm"))
		{
			cellcornerinterpolate(s, varinter);
			if (!strcmp(s.name, "fdrop"))
			{
				if (onlyf[j] == 'y')
				{
					for (i = 0; i < cellnumber; i++)
					{
						for (p = 0; p < 4; p++)
						{
							varinter[p + i * 4] = 1.0 - f1inter[p + i * 4] + 2.0 * varinter[p + i * 4];
							if (varinter[p + i * 4] < 0.0)
								varinter[p + i * 4] = 0.0;
							else if (varinter[p + i * 4] > 2.0)
								varinter[p + i * 4] = 2.0;
						}
					}
				}
			}
			else
			{
				if (!strcmp(s.name, "u.x"))
				{
					cellcornerinterpolate(u.y, varinter2);
					for (i = 0; i < cellnumber; i++)
					{
						varinter[0 + i * 4] = varinter[0 + i * 4] * cos(angle) - varinter2[0 + i * 4] * sin(angle);
						varinter[1 + i * 4] = varinter[1 + i * 4] * cos(angle) - varinter2[1 + i * 4] * sin(angle);
						varinter[2 + i * 4] = varinter[2 + i * 4] * cos(angle) - varinter2[2 + i * 4] * sin(angle);
						varinter[3 + i * 4] = varinter[3 + i * 4] * cos(angle) - varinter2[3 + i * 4] * sin(angle);
					}
				}
				if (!strcmp(s.name, "u.y"))
				{
					cellcornerinterpolate(u.x, varinter2);
					for (i = 0; i < cellnumber; i++)
					{
						varinter[0 + i * 4] = varinter[0 + i * 4] * cos(angle) + varinter2[0 + i * 4] * sin(angle);
						varinter[1 + i * 4] = varinter[1 + i * 4] * cos(angle) + varinter2[1 + i * 4] * sin(angle);
						varinter[2 + i * 4] = varinter[2 + i * 4] * cos(angle) + varinter2[2 + i * 4] * sin(angle);
						varinter[3 + i * 4] = varinter[3 + i * 4] * cos(angle) + varinter2[3 + i * 4] * sin(angle);
					}
				}
				if (onlyf[j] == 'y')
				{
					for (i = 0; i < cellnumber; i++)
					{
						varinter[0 + i * 4] *= f1inter[0 + i * 4];
						varinter[1 + i * 4] *= f1inter[1 + i * 4];
						varinter[2 + i * 4] *= f1inter[2 + i * 4];
						varinter[3 + i * 4] *= f1inter[3 + i * 4];
					}
				}
			}
			i = 0;
			foreach_leaf()
			{
				fprintf(fp, "\r\n%e", varinter[0 + i * 4]);
				fprintf(fp, "\r\n%e", varinter[1 + i * 4]);
				fprintf(fp, "\r\n%e", varinter[2 + i * 4]);
				fprintf(fp, "\r\n%e", varinter[3 + i * 4]);
				i++;
			}
			j++;
		}
	};
	i = 0;
	foreach_leaf()
	{
		fprintf(fp, "\r\n%d %d %d %d", i * 4 + 1, i * 4 + 2, i * 4 + 3, i * 4 + 4);
		i++;
	};
	fclose(fp);
	free(f1inter);
	free(varinter);
	free(varinter2);
	return;
}

void cellcornerinterpolate(scalar cellvar, double *pointvar)
{
	int i;
	double xtmp, ytmp, vtmp;
	;
	i = 0;
	foreach_leaf()
	{
		xtmp = x - 0.50 * Delta;
		ytmp = y - 0.50 * Delta;
		vtmp = interpolate(cellvar, xtmp, ytmp);
		if (vtmp == nodata)
			vtmp = 0.0;
		pointvar[0 + i * 4] = vtmp;
		;
		xtmp = x + 0.50 * Delta;
		ytmp = y - 0.50 * Delta;
		vtmp = interpolate(cellvar, xtmp, ytmp);
		if (vtmp == nodata)
			vtmp = 0.0;
		pointvar[1 + i * 4] = vtmp;
		;
		xtmp = x + 0.50 * Delta;
		ytmp = y + 0.50 * Delta;
		vtmp = interpolate(cellvar, xtmp, ytmp);
		if (vtmp == nodata)
			vtmp = 0.0;
		pointvar[2 + i * 4] = vtmp;
		;
		xtmp = x - 0.50 * Delta;
		ytmp = y + 0.50 * Delta;
		vtmp = interpolate(cellvar, xtmp, ytmp);
		if (vtmp == nodata)
			vtmp = 0.0;
		pointvar[3 + i * 4] = vtmp;
		;
		i++;
	}
	return;
}

void output_tecplot2D_Intrfc(char *name, scalar intrfc, double time, double rotationangle)
{
	int interfacepoints, cellnumber, iii;
	char typei[2], *typeintersect;
	double xi[2], yi[2], *xintersect, *yintersect;
	FILE *fp;
	const double angle = rotationangle * R_PI / 180.0;
	scalar alpha[];
	vector n[];
	;
	cellnumber = 0;
	foreach ()
		cellnumber++;
	xintersect = (double *)calloc(sizeof(double), cellnumber);
	yintersect = (double *)calloc(sizeof(double), cellnumber);
	typeintersect = (char *)calloc(sizeof(char), cellnumber);
	;
	reconstruction(intrfc, n, alpha);
	interfacepoints = 0;
	foreach_leaf()
	{
		if (intrfc[] > R_VOFLIMIT && intrfc[] < 1.0 - R_VOFLIMIT)
		{
			findintersectionpoints(n.x[], n.y[], alpha[], x, y, Delta, xi, yi, typei);
			xintersect[interfacepoints] = xi[0] * cos(angle) - yi[0] * sin(angle);
			yintersect[interfacepoints] = yi[0] * cos(angle) + xi[0] * sin(angle);
			typeintersect[interfacepoints] = typei[0];
			interfacepoints++;
			xintersect[interfacepoints] = xi[1] * cos(angle) - yi[1] * sin(angle);
			yintersect[interfacepoints] = yi[1] * cos(angle) + xi[1] * sin(angle);
			typeintersect[interfacepoints] = typei[1];
			interfacepoints++;
		}
	}
	fp = fopen(name, "w");
	switch (interfacepoints)
	{
	case 0:
	{
		fprintf(fp,
				"\r\nvariables = x, y"
				"\r\nZONE T=InterfaceData"
				"\r\nN = 3, E = 1, DATAPACKING = BLOCK, ZONETYPE = FETRIANGLE"
				"\r\nSOLUTIONTIME = %.4e",
				time);
		foreach ()
		{
			fprintf(fp, "\r\n%e %e %e", x, x, x);
			fprintf(fp, "\r\n%e %e %e", y, y, y);
			fprintf(fp, "\r\n1 2 3");
			break;
		}
		break;
	}
	default:
	{
		fprintf(fp,
				"\r\nvariables = x, y"
				"\r\nZONE T=InterfaceData"
				"\r\nN = %d, E = %d, DATAPACKING = BLOCK, ZONETYPE = FETRIANGLE"
				"\r\nSOLUTIONTIME = %.4e",
				3 * interfacepoints / 2, interfacepoints / 2, time);
		for (iii = 0; iii < interfacepoints - 1; iii += 2)
			fprintf(fp, "\r\n%e %e %e", xintersect[iii], xintersect[iii + 1], xintersect[iii]);
		for (iii = 0; iii < interfacepoints - 1; iii += 2)
			fprintf(fp, "\r\n%e %e %e", yintersect[iii], yintersect[iii + 1], yintersect[iii]);
		for (iii = 1; iii < 3 * interfacepoints / 2; iii += 3)
			fprintf(fp, "\r\n%d %d %d", iii, iii + 1, iii + 2);
		break;
	}
	}
	fclose(fp);
	;
	free(xintersect);
	free(yintersect);
}

scalar cal_intersectionarea()
{
	scalar alpha[];
	vector n[];
	scalar inarea[];
	char typei[2];
	double xi[2], yi[2], nx, ny, alp;
	;
	reconstruction(f, n, alpha);
	foreach ()
	{
		if (f[] > R_VOFLIMIT && f[] < 1.0 - R_VOFLIMIT)
		{
			nx = n.x[];
			ny = n.y[];
			alp = alpha[];
			findintersectionpoints(nx, ny, alp, x, y, Delta, xi, yi, typei);
#if AXI
			inarea[] = sqrt((xi[1] - xi[0]) * (xi[1] - xi[0]) + (yi[1] - yi[0]) * (yi[1] - yi[0])) * 2.0 * 3.1415926536 * 0.50 * (yi[1] + yi[0]);
#else
			inarea[] = sqrt((xi[1] - xi[0]) * (xi[1] - xi[0]) + (yi[1] - yi[0]) * (yi[1] - yi[0]));
#endif
		}
		else
			inarea[] = 0.0;
	}
	return inarea;
}

int findintersectionpoints(double nx, double ny, double alpha, double xc, double yc, double delta, double xi[2], double yi[2], char typei[2])
{
	int ppp = 0;
	double xtmp[2], ytmp[2], underflow = 1.0e-6;
	if (fabs(nx) < underflow)
	{
		ytmp[0] = (alpha - nx * (-0.50)) / ny;
		ytmp[1] = (alpha - nx * (+0.50)) / ny;
		xi[ppp] = xc + (-0.50) * delta;
		yi[ppp] = yc + (ytmp[0]) * delta;
		typei[ppp] = 'l';
		(ppp)++;
		xi[ppp] = xc + (+0.50) * delta;
		yi[ppp] = yc + (ytmp[1]) * delta;
		typei[ppp] = 'r';
		(ppp)++;
	}
	else if (fabs(ny) < underflow)
	{
		xtmp[0] = (alpha - ny * (-0.50)) / nx;
		xtmp[1] = (alpha - ny * (+0.50)) / nx;
		xi[ppp] = xc + (xtmp[0]) * delta;
		yi[ppp] = yc + (-0.50) * delta;
		typei[ppp] = 'b';
		(ppp)++;
		xi[ppp] = xc + (xtmp[1]) * delta;
		yi[ppp] = yc + (+0.50) * delta;
		typei[ppp] = 't';
		(ppp)++;
	}
	else
	{
		xtmp[0] = (alpha - ny * (-0.50)) / nx;
		xtmp[1] = (alpha - ny * (+0.50)) / nx;
		ytmp[0] = (alpha - nx * (-0.50)) / ny;
		ytmp[1] = (alpha - nx * (+0.50)) / ny;

		if (-0.50 <= ytmp[0] && ytmp[0] <= +0.50)
		{
			xi[ppp] = xc + (-0.50) * delta;
			yi[ppp] = yc + (ytmp[0]) * delta;
			typei[ppp] = 'l';
			(ppp)++;
		}
		if (-0.50 <= xtmp[0] && xtmp[0] <= +0.50)
		{
			xi[ppp] = xc + (xtmp[0]) * delta;
			yi[ppp] = yc + (-0.50) * delta;
			typei[ppp] = 'b';
			(ppp)++;
		}
		if (-0.50 <= ytmp[1] && ytmp[1] <= +0.50)
		{
			xi[ppp] = xc + (+0.50) * delta;
			yi[ppp] = yc + (ytmp[1]) * delta;
			typei[ppp] = 'r';
			(ppp)++;
		}
		if (-0.50 <= xtmp[1] && xtmp[1] <= +0.50)
		{
			xi[ppp] = xc + (xtmp[1]) * delta;
			yi[ppp] = yc + (+0.50) * delta;
			typei[ppp] = 't';
			(ppp)++;
		}
	}
	return true;
}

void output_tecplot3D(char *name, int iteration, double time, bool swap_xy)
{
	FILE *fp;
	int iii, cellnumber, *lcell;																			//, interfacepoints
	double *xcell, *ycell, *zcell, *fcell, *Fcell, *ocell, *ucell, *wcell, *vcell, *pcell, *nablaV, *dcell; // ocell: omega *** dcell: Delta of the cell *** lcell: level of the cell
	//double *xintersect, *yintersect, *nx, *ny, *aintersect; // aintersect: alpha intersect
	double divtmp; //, xi[2], yi[2];
	;
	scalar omega[]; //, alpha[];
	//vector n[];
	;
	vorticity(u, omega);
	//reconstruction(f, n, alpha);
	;
	cellnumber = 0;
	foreach ()
		cellnumber++;
	xcell = (double *)calloc(sizeof(double), cellnumber);
	ycell = (double *)calloc(sizeof(double), cellnumber);
	zcell = (double *)calloc(sizeof(double), cellnumber);
	dcell = (double *)calloc(sizeof(double), cellnumber);
	fcell = (double *)calloc(sizeof(double), cellnumber);
	Fcell = (double *)calloc(sizeof(double), cellnumber);
	ocell = (double *)calloc(sizeof(double), cellnumber);
	ucell = (double *)calloc(sizeof(double), cellnumber);
	vcell = (double *)calloc(sizeof(double), cellnumber);
	wcell = (double *)calloc(sizeof(double), cellnumber);
	pcell = (double *)calloc(sizeof(double), cellnumber);
	nablaV = (double *)calloc(sizeof(double), cellnumber);
	lcell = (int *)calloc(sizeof(int), cellnumber);
	;
	//xintersect = (double *)calloc(sizeof(double), cellnumber);
	//yintersect = (double *)calloc(sizeof(double), cellnumber);
	//aintersect = (double *)calloc(sizeof(double), cellnumber);
	//nx = (double *)calloc(sizeof(double), cellnumber);
	//ny = (double *)calloc(sizeof(double), cellnumber);
	;
	iii = 0;
	foreach ()
	{
		if (swap_xy)
		{
			xcell[iii] = z;
			ycell[iii] = y;
			zcell[iii] = x;
		}
		else
		{
			xcell[iii] = x;
			ycell[iii] = y;
			zcell[iii] = z;
		}
		dcell[iii] = Delta;
		fcell[iii] = f[];
		Fcell[iii] = f[] * (1.0 + fdrop[]);
		ocell[iii] = omega[];
		if (swap_xy)
		{
			ucell[iii] = u.z[];
			vcell[iii] = u.y[];
			wcell[iii] = u.x[];
		}
		else
		{
			ucell[iii] = u.x[];
			vcell[iii] = u.y[];
			wcell[iii] = u.z[];
		}
		pcell[iii] = p[];
		//nx[iii] = n.x[];
		//ny[iii] = n.y[];
		//aintersect[iii] = alpha[];
		lcell[iii] = level;
		divtmp = 0.;
		foreach_dimension()
			divtmp += u.x[1] - u.x[-1];
		divtmp /= 2.0 * Delta;
		nablaV[iii] = divtmp;
		iii++;
	}
	fp = fopen(name, "w");
	fprintf(fp, "variables = x, y, z, vof, u, v, w, p, omega, NablaV, CPU, level");
	fprintf(fp, "\r\nZONE T=\"tc[%f]\""
				"\r\nN = %d, E = %d, DATAPACKING = BLOCK, ZONETYPE = FEBRICK"
				"\r\nVARLOCATION = (NODAL, NODAL, NODAL"
				", CELLCENTERED, CELLCENTERED, CELLCENTERED"
				", CELLCENTERED, CELLCENTERED, CELLCENTERED"
				", CELLCENTERED, CELLCENTERED, CELLCENTERED)"
				"\r\nSOLUTIONTIME = %e",
			time, cellnumber * 8, cellnumber, time);
	for (iii = 0; iii < cellnumber; iii++)
	{
		fprintf(fp, "\r\n%e", xcell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", xcell[iii] + 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", xcell[iii] + 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", xcell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", xcell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", xcell[iii] + 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", xcell[iii] + 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", xcell[iii] - 0.50 * dcell[iii]);
	}
	for (iii = 0; iii < cellnumber; iii++)
	{
		fprintf(fp, "\r\n%e", ycell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", ycell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", ycell[iii] + 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", ycell[iii] + 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", ycell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", ycell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", ycell[iii] + 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", ycell[iii] + 0.50 * dcell[iii]);
	}
	for (iii = 0; iii < cellnumber; iii++)
	{
		fprintf(fp, "\r\n%e", zcell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", zcell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", zcell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", zcell[iii] - 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", zcell[iii] + 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", zcell[iii] + 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", zcell[iii] + 0.50 * dcell[iii]);
		fprintf(fp, "\r\n%e", zcell[iii] + 0.50 * dcell[iii]);
	}
	for (iii = 0; iii < cellnumber; iii++)
		fprintf(fp, "\r\n%e", Fcell[iii]);
	for (iii = 0; iii < cellnumber; iii++)
		fprintf(fp, "\r\n%e", ucell[iii]);
	for (iii = 0; iii < cellnumber; iii++)
		fprintf(fp, "\r\n%e", vcell[iii]);
	for (iii = 0; iii < cellnumber; iii++)
		fprintf(fp, "\r\n%e", wcell[iii]);
	for (iii = 0; iii < cellnumber; iii++)
		fprintf(fp, "\r\n%e", pcell[iii]);
	for (iii = 0; iii < cellnumber; iii++)
		fprintf(fp, "\r\n%e", ocell[iii]);
	for (iii = 0; iii < cellnumber; iii++)
		fprintf(fp, "\r\n%e", nablaV[iii]);
	for (iii = 0; iii < cellnumber; iii++)
		fprintf(fp, "\r\n%d", pid());
	for (iii = 0; iii < cellnumber; iii++)
		fprintf(fp, "\r\n%d", lcell[iii]);
	;
	for (iii = 0; iii < cellnumber; iii++)
		fprintf(fp, "\r\n%d %d %d %d %d %d %d %d", 1 + iii * 8, 2 + iii * 8, 3 + iii * 8, 4 + iii * 8, 5 + iii * 8, 6 + iii * 8, 7 + iii * 8, 8 + iii * 8);
	;
	fclose(fp);
	;
	free(nablaV);
	free(dcell);
	free(xcell);
	free(ycell);
	free(zcell);
	free(fcell);
	free(Fcell);
	free(ocell);
	free(ucell);
	free(vcell);
	free(wcell);
	free(pcell);
	free(lcell);
}

// type-> r: reqular, l: in liquid only, s: second interface
void output_tecplot2D_AllTime(char *name, double time, int iloop, char type, int NO, scalar fract, scalar myvar1, scalar myvar2, int swap_xy)
{
	FILE *fp;
	int cellnumber = 0;
	double xtmp, ytmp, vtmp[4];
	;
	foreach ()
		cellnumber++;
	switch (iloop)
	{
	case 0:
	{
		fp = fopen(name, "w");
		break;
	}
	default:
	{
		fp = fopen(name, "a");
		break;
	}
	}
	switch (NO)
	{
	case 2:
	{
		fprintf(fp, "variables = x, y, v1, v2");
		break;
	}
	default:
	{
		fprintf(fp, "variables = x, y, v");
		break;
	}
	}
	fprintf(fp, "\r\nZONE T=\"tc[%.3f]\""
				"\r\nN = %d, E = %d, DATAPACKING=BLOCK, ZONETYPE = FEQUADRILATERAL"
				"\r\nSOLUTIONTIME = %.3e",
			time, cellnumber * 4, cellnumber, time);
	switch (swap_xy)
	{
	case 0: // 2D and 3D case
	{
		foreach_leaf()
		{
			fprintf(fp, "\r\n%e", x - 0.50 * Delta);
			fprintf(fp, "\r\n%e", x + 0.50 * Delta);
			fprintf(fp, "\r\n%e", x + 0.50 * Delta);
			fprintf(fp, "\r\n%e", x - 0.50 * Delta);
		}
		foreach_leaf()
		{
			fprintf(fp, "\r\n%e", y - 0.50 * Delta);
			fprintf(fp, "\r\n%e", y - 0.50 * Delta);
			fprintf(fp, "\r\n%e", y + 0.50 * Delta);
			fprintf(fp, "\r\n%e", y + 0.50 * Delta);
		}
		break;
	}
	default: // axi-symmetric case
	{
		foreach_leaf()
		{
			fprintf(fp, "\r\n%e", y - 0.50 * Delta);
			fprintf(fp, "\r\n%e", y - 0.50 * Delta);
			fprintf(fp, "\r\n%e", y + 0.50 * Delta);
			fprintf(fp, "\r\n%e", y + 0.50 * Delta);
		}
		foreach_leaf()
		{
			fprintf(fp, "\r\n%e", x - 0.50 * Delta);
			fprintf(fp, "\r\n%e", x + 0.50 * Delta);
			fprintf(fp, "\r\n%e", x + 0.50 * Delta);
			fprintf(fp, "\r\n%e", x - 0.50 * Delta);
		}
		break;
	}
	}
	foreach_leaf()
	{
		xtmp = x - 0.50 * Delta;
		ytmp = y - 0.50 * Delta;
		vtmp[0] = calculate_output_value(type, xtmp, ytmp, fract, myvar1);
		;
		xtmp = x + 0.50 * Delta;
		ytmp = y - 0.50 * Delta;
		vtmp[1] = calculate_output_value(type, xtmp, ytmp, fract, myvar1);
		;
		xtmp = x + 0.50 * Delta;
		ytmp = y + 0.50 * Delta;
		vtmp[2] = calculate_output_value(type, xtmp, ytmp, fract, myvar1);
		;
		xtmp = x - 0.50 * Delta;
		ytmp = y + 0.50 * Delta;
		vtmp[3] = calculate_output_value(type, xtmp, ytmp, fract, myvar1);
		;
		fprintf(fp, "\r\n%e", vtmp[0]);
		fprintf(fp, "\r\n%e", vtmp[1]);
		fprintf(fp, "\r\n%e", vtmp[2]);
		fprintf(fp, "\r\n%e", vtmp[3]);
	}
	switch (NO)
	{
	case 2:
	{
		foreach_leaf()
		{
			xtmp = x - 0.50 * Delta;
			ytmp = y - 0.50 * Delta;
			vtmp[0] = calculate_output_value(type, xtmp, ytmp, fract, myvar2);
			;
			xtmp = x + 0.50 * Delta;
			ytmp = y - 0.50 * Delta;
			vtmp[1] = calculate_output_value(type, xtmp, ytmp, fract, myvar2);
			;
			xtmp = x + 0.50 * Delta;
			ytmp = y + 0.50 * Delta;
			vtmp[2] = calculate_output_value(type, xtmp, ytmp, fract, myvar2);
			;
			xtmp = x - 0.50 * Delta;
			ytmp = y + 0.50 * Delta;
			vtmp[3] = calculate_output_value(type, xtmp, ytmp, fract, myvar2);
			;
			fprintf(fp, "\r\n%e", vtmp[0]);
			fprintf(fp, "\r\n%e", vtmp[1]);
			fprintf(fp, "\r\n%e", vtmp[2]);
			fprintf(fp, "\r\n%e", vtmp[3]);
		}
		break;
	}
	default:
	{
		break;
	}
	};
	cellnumber = 0;
	foreach_leaf()
	{
		fprintf(fp, "\r\n%d %d %d %d", 0 + cellnumber * 4 + 1, 1 + cellnumber * 4 + 1, 2 + cellnumber * 4 + 1, 3 + cellnumber * 4 + 1);
		cellnumber++;
	};
	fclose(fp);
}

// type-> r: reqular, l: in liquid only, s: second interface
double calculate_output_value(char type, double xxx, double yyy, scalar fract, scalar myvar)
{
	double vtmp, ftmp, rtmp = 0.0;
	vtmp = interpolate(myvar, xxx, yyy);
	if (vtmp == nodata)
		vtmp = 0.0;
	switch (type)
	{
	case 'r':
	case 'R':
	{
		rtmp = vtmp;
		break;
	}
	case 'l':
	case 'L':
	{
		ftmp = interpolate(fract, xxx, yyy);
		if (ftmp == nodata)
			ftmp = 0.0;
		rtmp = vtmp * ftmp;
		break;
	}
	case 's':
	case 'S':
	{
		ftmp = interpolate(fract, xxx, yyy);
		if (ftmp == nodata)
			ftmp = 0.0;
		rtmp = 1.0 - ftmp + 2.0 * vtmp;
		break;
	}
	default:
	{
		printf("Flag error in calculate_output_value.\r\n");
		break;
	}
	}
	return rtmp;
}

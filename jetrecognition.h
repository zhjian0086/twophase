int WriteTmprlData(char *name, int Iteration, double tc, double **xyDKO, int *JetRoots);
int WriteExtremums(char *name, int Iteration, double tc, double **xyDKO, double *kppF, double *omgF, int CellNO, int *Kindex, int KindexNO, int *Oindex, int OindexNO, int **JetRoots);
int Write2DJetRoot(char *name, int Iteration, double tc, double **xyDKO, double *kppF, double *omgF, int CellNO, int *Kindex, int KindexNO, int *Oindex, int OindexNO, int **JetRoots);

int OrderInterface(double tc, scalar ux, scalar uy, scalar kappa, scalar omega);
int FindAllChainLoops(double tc, int interfaceNO, const double deltaMIN, double **xyDelta, int **ChainLoopIndex, char writeNBinfo);
int FindJetRoots(int iteration, double tc, scalar kappa, scalar omega, char writeNBinfo, int **JetRoots);
int JetRecognition(int iteration, double tc, double **xyDKOjet, int CellNO, int **JetRoots, char writeNBinfo);

int findlocalmaxminIndependent(double *var, int NO, int *index, char *type, int *noindex);
int filtervariables(double *Var, int CellNO, double *VarF, char smoothing_lowpass, int filterlevels, int *filterrepeats, int *filterpoints);
int averagesmooth(double *varin, int NO, int order, double *varout);
int lowpassfilter(double *varin, int NO, int order, double *varout);
int mergeloops(int L1, int L2, int P1, int P2, int **ChainLoopIndex);
int checkneighbcorner(int C1, int C2, int **cellNC);
int removecell(int i, int **cellneighb, int **cellcorner, char *flag);

#define MARGINPOINTS 			10
#define MARGINPOINTSSKIP		5
#define THRESHOLDMAXVALUE 		5	// percent
#define THRESHOLDDROPK 			5	// percent
#define THRESHOLDINITIALGAP 	50	// points

int FindJetRoots(int iteration, double tc, scalar kappa, scalar omega, char writeNBinfo, int **JetRoots)
{
	int i, j, k;
	int interfaceNO, MaxChain, MaxChainCellNO, itmp[2];
	const int VarNO = 5; // x, y, Delta, Kappa, Omega
	int *IndexTotalLeaf, **ChainLoopIndex;
	double **xyDKO, **xyDKOjet; // xyDKO: x, y, Delta, Kappa, Omega
	double *kpp, *omg;
	char **VarFileName;
	double deltaMIN;
	;
	// count interface cells
	;
	deltaMIN = 1.0e12;
	interfaceNO = 0;
	foreach ()
	{
		if (f[] > R_VOFLIMIT && f[] < 1.0 - R_VOFLIMIT)
		{
			if (deltaMIN > 0.50 * Delta)
				deltaMIN = 0.50 * Delta;
			interfaceNO++;
		}
	};
	// memory allocation
	;
	kpp = (double *)calloc(interfaceNO, sizeof(double));
	omg = (double *)calloc(interfaceNO, sizeof(double));
	xyDKOjet = (double **)calloc(VarNO, sizeof(double *));
	xyDKO = (double **)calloc(VarNO, sizeof(double *));
	for (i = 0; i < VarNO; i++)
	{
		xyDKOjet[i] = (double *)calloc(interfaceNO, sizeof(double));
		xyDKO[i] = (double *)calloc(interfaceNO, sizeof(double));
	}
	IndexTotalLeaf = (int *)calloc(interfaceNO, sizeof(int));
	ChainLoopIndex = (int **)calloc(interfaceNO, sizeof(int *));
	for (i = 0; i < interfaceNO; i++)
	{
		ChainLoopIndex[i] = (int *)calloc(interfaceNO, sizeof(int));
	}
	VarFileName = (char **)calloc(100, sizeof(char *));
	for (i = 0; i < 100; i++)
		VarFileName[i] = (char *)calloc(100, sizeof(char));
	;
	// transfer the cell data to array
	;
	k = 0;
	i = 0;
	foreach_leaf()
	{
		if (f[] > R_VOFLIMIT && f[] < 1.0 - R_VOFLIMIT)
		{
			xyDKO[0][i] = x;
			xyDKO[1][i] = y;
			xyDKO[2][i] = 0.50 * Delta;
			xyDKO[3][i] = kappa[];
			xyDKO[4][i] = omega[];
			IndexTotalLeaf[i] = k;
			i++;
		}
		k++;
	};
	// find the chain loops of interface cells
	;
	FindAllChainLoops(tc, interfaceNO, deltaMIN, xyDKO, ChainLoopIndex, writeNBinfo);
	;
	// find the largest chain
	;
	MaxChain = 0;
	MaxChainCellNO = 0;
	for (k = 1; k <= ChainLoopIndex[0][0]; k++)
	{
		if (MaxChainCellNO < ChainLoopIndex[k][0])
		{
			MaxChainCellNO = ChainLoopIndex[k][0];
			MaxChain = k;
		}
	}
	for (k = 1; k <= ChainLoopIndex[MaxChain][0]; k++)
	{
		j = ChainLoopIndex[MaxChain][k];
		for (i = 0; i < VarNO; i++)
			xyDKOjet[i][k - 1] = xyDKO[i][j];
	};
	// find jet roots
	;
	JetRecognition(iteration, tc, xyDKOjet, MaxChainCellNO, JetRoots, writeNBinfo);
	for (k = 1; k <= JetRoots[0][0]; k++)
	{
		itmp[0] = JetRoots[k][0];
		itmp[1] = JetRoots[k][1];
		JetRoots[k][0] = ChainLoopIndex[MaxChain][itmp[0] + 1];
		JetRoots[k][1] = ChainLoopIndex[MaxChain][itmp[1] + 1];
	};
	// export the data
	;
	switch (writeNBinfo)
	{
	case 'y':
	case 'Y':
	{
		WriteTmprlData("out-jet-temporal-1.plt", iteration, tc, xyDKO, JetRoots[1]);
		WriteTmprlData("out-jet-temporal-2.plt", iteration, tc, xyDKO, JetRoots[2]);
		WriteTmprlData("out-jet-temporal-3.plt", iteration, tc, xyDKO, JetRoots[3]);
		break;
	}
	}
	// leaf-index of the jet-roots
	;
	for (k = 1; k <= JetRoots[0][0]; k++)
	{
		itmp[0] = JetRoots[k][0];
		itmp[1] = JetRoots[k][1];
		JetRoots[k][0] = IndexTotalLeaf[itmp[0]];
		JetRoots[k][1] = IndexTotalLeaf[itmp[1]];
	};
	// free the allocated memory
	;
	free(IndexTotalLeaf);
	for (i = 0; i < interfaceNO; i++)
	{
		free(ChainLoopIndex[i]);
	}
	free(ChainLoopIndex);
	for (i = 0; i < VarNO; i++)
	{
		free(xyDKO[i]);
		free(xyDKOjet[i]);
	}
	free(xyDKO);
	free(kpp);
	free(omg);
	free(xyDKOjet);
	for (i = 0; i < 100; i++)
		free(VarFileName[i]);
	free(VarFileName);
	return true;
}

int JetRecognition(int iteration, double tc, double **xyDKOjet, int CellNO, int **JetRoots, char writeNBinfo)
{
	int i, k, p;
	int JetRootMethod, jetp;
	int filterrepeats[100], filterpoints[100];
	double *kppF, *omgF, *VarI, *VarO;
	int *Kindex, *Oindex, KindexNO, OindexNO;
	char flag, *Ktype, *Otype;
	;
	// allocating memory
	;
	kppF = (double *)calloc(CellNO, sizeof(double));
	omgF = (double *)calloc(CellNO, sizeof(double));
	VarI = (double *)calloc(CellNO, sizeof(double));
	VarO = (double *)calloc(CellNO, sizeof(double));
	Kindex = (int *)calloc(1000, sizeof(int));
	Oindex = (int *)calloc(1000, sizeof(int));
	Ktype = (char *)calloc(1000, sizeof(char));
	Otype = (char *)calloc(1000, sizeof(char));
	;
	// filtering Kappa and Omega
	;
	filterrepeats[0] = 8;
	filterrepeats[1] = 4;
	filterrepeats[2] = 2;
	filterrepeats[3] = 1;
	filterrepeats[4] = 2;
	filterrepeats[5] = 4;
	filterrepeats[6] = 8;
	filterpoints[0] = 1;
	filterpoints[1] = 2;
	filterpoints[2] = 4;
	filterpoints[3] = 8;
	filterpoints[4] = 4;
	filterpoints[5] = 2;
	filterpoints[6] = 1;
	filtervariables(xyDKOjet[3], CellNO, kppF, 's', 7, filterrepeats, filterpoints);
	filtervariables(xyDKOjet[4], CellNO, omgF, 's', 7, filterrepeats, filterpoints);
	;
	// find local maximums
	;
	findlocalmaxminIndependent(kppF, CellNO, Kindex, Ktype, &KindexNO);
	findlocalmaxminIndependent(omgF, CellNO, Oindex, Otype, &OindexNO);
	;
	// using 3 methods two find the jet roots
	;
	JetRoots[0][0] = 3;
	;
	// analyze: Kmin start and Kmin end
	;
	JetRootMethod = 1;
	flag = 'n';
	jetp = 0;
	JetRoots[JetRootMethod][jetp] = -1;
	for (p = 0; p < KindexNO; p++)
	{
		if (Ktype[p] == 'n')
		{
			JetRoots[JetRootMethod][jetp] = Kindex[p];
			flag = 'y';
		}
		if (flag == 'y')
			break;
	}
	flag = 'n';
	jetp = 1;
	JetRoots[JetRootMethod][jetp] = -1;
	for (p = KindexNO - 1; p >= 0; p--)
	{
		if (Ktype[p] == 'n')
		{
			JetRoots[JetRootMethod][jetp] = Kindex[p];
			flag = 'y';
		}
		if (flag == 'y')
			break;
	}
	if (JetRoots[JetRootMethod][0] == -1 || JetRoots[JetRootMethod][1] == -1)
	{
		JetRoots[JetRootMethod][0] = -1;
		JetRoots[JetRootMethod][1] = -1;
	};
	// analyze: Kmin-Omin start and Kmin-Omax end
	;
	JetRootMethod = 2;
	flag = 'n';
	jetp = 0;
	JetRoots[JetRootMethod][jetp] = -1;
	for (p = 0; p < KindexNO; p++)
	{
		for (i = 0; i < OindexNO; i++)
		{
			for (k = Oindex[i] - MARGINPOINTS; k <= Oindex[i] + MARGINPOINTS; k++)
			{
				if (Kindex[p] == k && Ktype[p] == 'n' && Otype[i] == 'n')
				{
					JetRoots[JetRootMethod][jetp] = Kindex[p];
					flag = 'y';
				}
			}
		}
		if (flag == 'y')
			break;
	}
	flag = 'n';
	jetp = 1;
	JetRoots[JetRootMethod][jetp] = -1;
	for (p = KindexNO - 1; p >= 0; p--)
	{
		for (i = 0; i < OindexNO; i++)
		{
			for (k = Oindex[i] - MARGINPOINTS; k <= Oindex[i] + MARGINPOINTS; k++)
			{
				if (Kindex[p] == k && Ktype[p] == 'n' && Otype[i] == 'x')
				{
					JetRoots[JetRootMethod][jetp] = Kindex[p];
					flag = 'y';
				}
			}
		}
		if (flag == 'y')
			break;
	}
	if (JetRoots[JetRootMethod][0] == -1 || JetRoots[JetRootMethod][1] == -1)
	{
		JetRoots[JetRootMethod][0] = -1;
		JetRoots[JetRootMethod][1] = -1;
	};
	// analyze: Omin start and Omax end
	;
	JetRootMethod = 3;
	flag = 'n';
	jetp = 0;
	JetRoots[JetRootMethod][jetp] = -1;
	for (p = 0; p < OindexNO; p++)
	{
		if (Otype[p] == 'n')
		{
			JetRoots[JetRootMethod][jetp] = Oindex[p];
			flag = 'y';
		}
		if (flag == 'y')
			break;
	}
	flag = 'n';
	jetp = 1;
	JetRoots[JetRootMethod][jetp] = -1;
	for (p = OindexNO - 1; p >= 0; p--)
	{
		if (Otype[p] == 'x')
		{
			JetRoots[JetRootMethod][jetp] = Oindex[p];
			flag = 'y';
		}
		if (flag == 'y')
			break;
	}
	if (JetRoots[JetRootMethod][0] == -1 || JetRoots[JetRootMethod][1] == -1)
	{
		JetRoots[JetRootMethod][0] = -1;
		JetRoots[JetRootMethod][1] = -1;
	};
	// writing files
	;
	switch (writeNBinfo)
	{
	case 'y':
	case 'Y':
	{
		WriteExtremums("out-jet-extremum.plt", iteration, tc, xyDKOjet, kppF, omgF, CellNO, Kindex, KindexNO, Oindex, OindexNO, JetRoots);
		Write2DJetRoot("out-jet-2D-roots.plt", iteration, tc, xyDKOjet, kppF, omgF, CellNO, Kindex, KindexNO, Oindex, OindexNO, JetRoots);
		break;
	}
	}
	// free memory
	;
	free(kppF);
	free(omgF);
	free(VarI);
	free(VarO);
	free(Kindex);
	free(Oindex);
	free(Ktype);
	free(Otype);
	;
	return true;
}

int filtervariables(double *Var, int CellNO, double *VarF, char smoothing_lowpass, int filterlevels, int *filterrepeats, int *filterpoints)
{
	int i, l, k;
	double *VarI, *VarO;
	;
	// allocating memory
	;
	VarI = (double *)calloc(CellNO, sizeof(double));
	VarO = (double *)calloc(CellNO, sizeof(double));
	;
	// filtering
	;
	for (i = 0; i < CellNO; i++)
		VarI[i] = Var[i];
	switch (smoothing_lowpass)
	{
	case 's':
	case 'S':
	{
		for (l = 0; l < filterlevels; l++)
		{
			for (k = 0; k < filterrepeats[l]; k++)
			{
				averagesmooth(VarI, CellNO, filterpoints[l], VarO);
				averagesmooth(VarO, CellNO, filterpoints[l], VarI);
			}
		}
		break;
	}
	case 'l':
	case 'L':
	{
		for (l = 0; l < filterlevels; l++)
		{
			for (k = 0; k < filterrepeats[l]; k++)
			{
				lowpassfilter(VarI, CellNO, filterpoints[l], VarO);
				lowpassfilter(VarO, CellNO, filterpoints[l], VarI);
			}
		}
		break;
	}
	default:
	{
		printf("Wrong parameter in filtering!!\r\n");
		break;
	}
	}
	for (i = 0; i < CellNO; i++)
		VarF[i] = VarI[i];
	;
	// free the memory
	;
	free(VarI);
	free(VarO);
	return true;
}

// ChainLoopIndex[0][0]: number of the loops
// ChainLoopIndex[i][0]: size of the loop (i)
// ChainLoopIndex[i][1-j]: indeces of loop i start from 1-j
int FindAllChainLoops(double tc, int interfaceNO, const double deltaMIN, double **xyDelta, int **ChainLoopIndex, char writeNBinfo)
{
	int i, j, k, pj, ki, m, n, p, ibfr, icrn, looptotalNO, loopindexNO;
	int **cellneighb, **cellcorner, *cellloopno, *cellloopindex, *CellNB;
	char *flag, iflag, mainflag;
	double xtmp, ytmp, vtmp;
	double **f1;
	const double vofmid = 0.50;
	;
	// memory allocation
	;
	cellloopno = (int *)calloc(interfaceNO, sizeof(int));
	CellNB = (int *)calloc(interfaceNO, sizeof(int));
	cellloopindex = (int *)calloc(interfaceNO, sizeof(int));
	cellneighb = (int **)calloc(interfaceNO, sizeof(int *));
	cellcorner = (int **)calloc(interfaceNO, sizeof(int *));
	f1 = (double **)calloc(interfaceNO, sizeof(double *));
	for (j = 0; j < interfaceNO; j++)
	{
		cellneighb[j] = (int *)calloc(10, sizeof(int));
		cellcorner[j] = (int *)calloc(10, sizeof(int));
		f1[j] = (double *)calloc(10, sizeof(double));
	}
	flag = (char *)calloc(interfaceNO, sizeof(char));
	;
	// initializing
	;
	for (i = 0; i < interfaceNO; i++)
		flag[i] = 'n';
	;
	// fraction in corner points of a cell
	;
	for (i = 0; i < interfaceNO; i++)
	{
		xtmp = xyDelta[0][i];
		ytmp = xyDelta[1][i];
		vtmp = interpolate(f, xtmp, ytmp);
		if (vtmp == nodata)
			vtmp = 0.0;
		f1[i][0] = vtmp; // center
		;
		xtmp = xyDelta[0][i] - xyDelta[2][i];
		ytmp = xyDelta[1][i] - xyDelta[2][i];
		vtmp = interpolate(f, xtmp, ytmp);
		if (vtmp == nodata)
			vtmp = 0.0;
		f1[i][1] = vtmp; // left-bottom
		;
		xtmp = xyDelta[0][i] + xyDelta[2][i];
		ytmp = xyDelta[1][i] - xyDelta[2][i];
		vtmp = interpolate(f, xtmp, ytmp);
		if (vtmp == nodata)
			vtmp = 0.0;
		f1[i][2] = vtmp; // right-bottom
		;
		xtmp = xyDelta[0][i] + xyDelta[2][i];
		ytmp = xyDelta[1][i] + xyDelta[2][i];
		vtmp = interpolate(f, xtmp, ytmp);
		if (vtmp == nodata)
			vtmp = 0.0;
		f1[i][3] = vtmp; // right-top
		;
		xtmp = xyDelta[0][i] - xyDelta[2][i];
		ytmp = xyDelta[1][i] + xyDelta[2][i];
		vtmp = interpolate(f, xtmp, ytmp);
		if (vtmp == nodata)
			vtmp = 0.0;
		f1[i][4] = vtmp; // left-top
	};
	// find neighbors / corners for each cell
	;
	for (j = 0; j < interfaceNO; j++)
	{
		for (i = j + 1; i < interfaceNO; i++)
		{
			if ((fabs(fabs(xyDelta[0][i] - xyDelta[0][j]) - fabs(xyDelta[2][i] + xyDelta[2][j])) < 1.0e-6 &&
				 fabs(fabs(xyDelta[1][i] - xyDelta[1][j]) - fabs(xyDelta[2][i] - xyDelta[2][j])) < 1.0e-6) ||
				(fabs(fabs(xyDelta[0][i] - xyDelta[0][j]) - fabs(xyDelta[2][i] - xyDelta[2][j])) < 1.0e-6 &&
				 fabs(fabs(xyDelta[1][i] - xyDelta[1][j]) - fabs(xyDelta[2][i] + xyDelta[2][j])) < 1.0e-6))
			{
				cellneighb[j][++cellneighb[j][0]] = i;
				cellneighb[i][++cellneighb[i][0]] = j;
			}
			if (fabs(fabs(xyDelta[0][i] - xyDelta[0][j]) - fabs(xyDelta[2][i] + xyDelta[2][j])) < 1.0e-6 &&
				fabs(fabs(xyDelta[1][i] - xyDelta[1][j]) - fabs(xyDelta[2][i] + xyDelta[2][j])) < 1.0e-6)
			{
				cellcorner[j][++cellcorner[j][0]] = i;
				cellcorner[i][++cellcorner[i][0]] = j;
			}
		}
	};
	// filter neighbors
	;
	// // big 3 neighbor cells (modify from 2-3-3-2 to 2-2-2-2)
	for (i = 0; i < interfaceNO; i++)
	{
		if (fabs(xyDelta[2][i] - 2.0 * deltaMIN) < 1.0e-6 && cellneighb[i][0] == 3)
		{
			// neighbors
			n = 0;
			m = 0;
			p = -1;
			for (j = 1; j <= cellneighb[i][0]; j++)
			{
				if (cellneighb[cellneighb[i][j]][0] == 3)
				{
					n++;
					p = cellneighb[i][j];
				}
				if (cellneighb[cellneighb[i][j]][0] == 2 || cellneighb[cellneighb[i][j]][0] == 1)
					m++;
			}
			if (n == 1 && m == 2)
			{
				ibfr = 0;
				for (j = 1; j <= cellneighb[p][0]; j++)
					if (i != cellneighb[p][j])
						cellneighb[p][++ibfr] = cellneighb[p][j];
				cellneighb[p][0]--;
				ibfr = 0;
				for (j = 1; j <= cellneighb[i][0]; j++)
					if (p != cellneighb[i][j])
						cellneighb[i][++ibfr] = cellneighb[i][j];
				cellneighb[i][0]--;
			}
		}
	}
	// remove 4 (or more) neighbor cells
	p = 0;
	for (i = 0; i < interfaceNO; i++)
		if (cellneighb[i][0] >= 4)
			CellNB[p++] = i;
	for (m = 0; m < p; m++)
	{
		i = CellNB[m];
		removecell(i, cellneighb, cellcorner, flag);
	}
	// remove unnecessary cells with 3 neighbors
	for (p = 0; p < interfaceNO; p++)
	{
		if (fabs(xyDelta[2][p] - 2.0 * deltaMIN) < 1.0e-6 && cellneighb[p][0] == 3)
		{
			i = cellneighb[p][1];
			j = cellneighb[p][2];
			k = cellneighb[p][3];
			n = 0;
			if (checkneighbcorner(i, j, cellneighb) || checkneighbcorner(i, j, cellcorner))
				n++;
			if (checkneighbcorner(j, k, cellneighb) || checkneighbcorner(j, k, cellcorner))
				n++;
			if (checkneighbcorner(k, i, cellneighb) || checkneighbcorner(k, i, cellcorner))
				n++;
			if (n == 2)
			{
				// remove cell[i]
				removecell(p, cellneighb, cellcorner, flag);
			}
		}
	};
	// find the chain of cell-loops
	;
	looptotalNO = 0;
	do
	{
		// find the first cell
		mainflag = 'n';
		loopindexNO = 0;
		for (i = 0; i < interfaceNO; i++)
		{
			if (cellneighb[i][0] == 1 && cellcorner[i][0] == 0 && flag[i] == 'n')
			{
				ChainLoopIndex[++looptotalNO][++loopindexNO] = i;
				flag[i] = 'y';
				mainflag = 'y';
				i = interfaceNO + 1;
			}
		}
		// find the chain of points
		if (mainflag == 'y')
		{
			ibfr = ChainLoopIndex[looptotalNO][loopindexNO];
			ChainLoopIndex[looptotalNO][++loopindexNO] = cellneighb[ibfr][1];
			icrn = ChainLoopIndex[looptotalNO][loopindexNO];
			flag[icrn] = 'y';
			do
			{
				iflag = 'n';
				switch (cellneighb[icrn][0])
				{
				case 0:
				{
					for (i = 1; i <= cellcorner[icrn][0]; i++)
					{
						if (flag[cellcorner[icrn][i]] == 'n')
						{
							ibfr = icrn;
							icrn = cellcorner[icrn][i];
							ChainLoopIndex[looptotalNO][++loopindexNO] = icrn;
							flag[ChainLoopIndex[looptotalNO][loopindexNO]] = 'y';
							iflag = 'y';
							i = cellcorner[icrn][0] + 1;
						}
					}
					break;
				}
				case 1:
				{
					if (flag[cellneighb[icrn][1]] == 'y')
					{
						for (i = 1; i <= cellcorner[icrn][0]; i++)
						{
							if (flag[cellcorner[icrn][i]] == 'n')
							{
								ibfr = icrn;
								icrn = cellcorner[icrn][i];
								ChainLoopIndex[looptotalNO][++loopindexNO] = icrn;
								flag[ChainLoopIndex[looptotalNO][loopindexNO]] = 'y';
								iflag = 'y';
								i = cellcorner[icrn][0] + 1;
							}
						}
					}
					else
					{
						ibfr = icrn;
						icrn = cellneighb[icrn][1];
						ChainLoopIndex[looptotalNO][++loopindexNO] = icrn;
						flag[ChainLoopIndex[looptotalNO][loopindexNO]] = 'y';
						iflag = 'y';
					}
					break;
				}
				case 2:
				{
					if (flag[cellneighb[icrn][1]] == 'y' && flag[cellneighb[icrn][2]] == 'n')
					{
						ibfr = icrn;
						icrn = cellneighb[icrn][2];
						ChainLoopIndex[looptotalNO][++loopindexNO] = icrn;
						flag[ChainLoopIndex[looptotalNO][loopindexNO]] = 'y';
						iflag = 'y';
					}
					else if (flag[cellneighb[icrn][1]] == 'n' && flag[cellneighb[icrn][2]] == 'y')
					{
						ibfr = icrn;
						icrn = cellneighb[icrn][1];
						ChainLoopIndex[looptotalNO][++loopindexNO] = icrn;
						flag[ChainLoopIndex[looptotalNO][loopindexNO]] = 'y';
						iflag = 'y';
					}
					break;
				}
				case 3:
				{
					vtmp = 1.0e12;
					k = -1;
					for (j = 1; j <= cellneighb[icrn][0]; j++)
					{
						i = cellneighb[icrn][j];
						if (flag[i] == 'y')
							;
						else
						{
							// if (vtmp > fabs(f1[i][0] - vofmid))
							// {
							// 	vtmp = fabs(f1[i][0] - vofmid);
							// 	k = i;
							// }
							if (vtmp > fabs(0.250 * (f1[i][1] + f1[i][2] + f1[i][3] + f1[i][4]) - vofmid))
							{
								vtmp = fabs(0.250 * (f1[i][1] + f1[i][2] + f1[i][3] + f1[i][4]) - vofmid);
								k = i;
							}
						}
					}
					if (k != -1)
					{
						ibfr = icrn;
						icrn = k;
						ChainLoopIndex[looptotalNO][++loopindexNO] = icrn;
						flag[ChainLoopIndex[looptotalNO][loopindexNO]] = 'y';
						iflag = 'y';
					}
					break;
				}
				}
			} while (iflag == 'y');
			ChainLoopIndex[looptotalNO][0] = loopindexNO;
		}
	} while (mainflag == 'y');
	ChainLoopIndex[0][0] = looptotalNO;
	;
	// check the connection of chains
	;
	for (k = 1; k < ChainLoopIndex[0][0]; k++)
	{
		for (p = k + 1; p <= ChainLoopIndex[0][0]; p++)
		{
			for (i = 1; i <= ChainLoopIndex[k][0]; i++)
			{
				ki = ChainLoopIndex[k][i];
				for (j = 1; j <= ChainLoopIndex[p][0]; j++)
				{
					pj = ChainLoopIndex[p][j];
					iflag = 'n';
					for (n = 1; n <= cellneighb[ki][0]; n++)
					{
						if (cellneighb[ki][n] == pj)
						{
							// connection found
							mergeloops(k, p, i, j, ChainLoopIndex);
							printf("Merge loops %d and %d from neighbor-point, from %d\r\n", k, p, pj);
							iflag = 'y';
						}
					}
					if (iflag == 'n')
					{
						for (n = 1; n <= cellcorner[ki][0]; n++)
						{
							if (cellcorner[ki][n] == pj)
							{
								// connection found
								mergeloops(k, p, i, j, ChainLoopIndex);
								printf("Merge loops %d and %d from corner-points, from %d\r\n", k, p, pj);
							}
						}
					}
				}
			}
		}
	}
	for (k = 1; k < ChainLoopIndex[0][0]; k++)
	{
		for (p = k + 1; p <= ChainLoopIndex[0][0]; p++)
		{
			for (i = 1; i <= ChainLoopIndex[k][0]; i++)
			{
				ki = ChainLoopIndex[k][i];
				for (j = 1; j <= ChainLoopIndex[p][0]; j++)
				{
					pj = ChainLoopIndex[p][j];
					iflag = 'n';
					for (n = 1; n <= cellneighb[ki][0]; n++)
					{
						if (cellneighb[ki][n] == pj)
						{
							// connection found
							mergeloops(k, p, i, j, ChainLoopIndex);
							printf("Merge loops %d and %d from neighbor-point, from %d\r\n", k, p, pj);
							iflag = 'y';
						}
					}
					if (iflag == 'n')
					{
						for (n = 1; n <= cellcorner[ki][0]; n++)
						{
							if (cellcorner[ki][n] == pj)
							{
								// connection found
								mergeloops(k, p, i, j, ChainLoopIndex);
								printf("Merge loops %d and %d from corner-points, from %d\r\n", k, p, pj);
							}
						}
					}
				}
			}
		}
	}
	for (k = 1; k < ChainLoopIndex[0][0]; k++)
	{
		for (p = k + 1; p <= ChainLoopIndex[0][0]; p++)
		{
			for (i = 1; i <= ChainLoopIndex[k][0]; i++)
			{
				ki = ChainLoopIndex[k][i];
				for (j = 1; j <= ChainLoopIndex[p][0]; j++)
				{
					pj = ChainLoopIndex[p][j];
					iflag = 'n';
					for (n = 1; n <= cellneighb[ki][0]; n++)
					{
						if (cellneighb[ki][n] == pj)
						{
							// connection found
							mergeloops(k, p, i, j, ChainLoopIndex);
							printf("Merge loops %d and %d from neighbor-point, from %d\r\n", k, p, pj);
							iflag = 'y';
						}
					}
					if (iflag == 'n')
					{
						for (n = 1; n <= cellcorner[ki][0]; n++)
						{
							if (cellcorner[ki][n] == pj)
							{
								// connection found
								mergeloops(k, p, i, j, ChainLoopIndex);
								printf("Merge loops %d and %d from corner-points, from %d\r\n", k, p, pj);
							}
						}
					}
				}
			}
		}
	};
	// check the direction of chains
	;
	for (k = 1; k <= ChainLoopIndex[0][0]; k++)
	{
		if (xyDelta[1][ChainLoopIndex[k][1]] > xyDelta[1][ChainLoopIndex[k][ChainLoopIndex[k][0]]])
		{
			for (j = 1; j <= ChainLoopIndex[k][0] / 2; j++)
			{
				i = ChainLoopIndex[k][0] + 1 - j;
				icrn = ChainLoopIndex[k][j];
				ChainLoopIndex[k][j] = ChainLoopIndex[k][i];
				ChainLoopIndex[k][i] = icrn;
			}
		}
	};
	// update indeces
	;
	for (i = 0; i < interfaceNO; i++)
	{
		cellloopno[i] = -1;
		cellloopindex[i] = -1;
	}
	for (k = 1; k <= ChainLoopIndex[0][0]; k++)
	{
		for (j = 1; j <= ChainLoopIndex[k][0]; j++)
		{
			i = ChainLoopIndex[k][j];
			cellloopno[i] = k;
			cellloopindex[i] = j;
		}
	};
	// write the output files
	;
	//	printf("Total Interfaces: %d\r\n", interfaceNO);
	//	printf("Total Chains: %d\r\n", looptotalNO);
	//	for (i = 1; i <= ChainLoopIndex[0][0]; i++)
	//		printf("Chain Loop Number: %d\r\n", ChainLoopIndex[i][0]);
	switch (writeNBinfo)
	{
	case 'y':
	case 'Y':
	{
		// data all neighbors and corners
		FILE *fp;
		char name[500], ctmp[100];
		strcpy(name, "out-jet-NB-all-");
		sprintf(ctmp, "%.4f.plt", tc);
		strcat(name, ctmp);
		fp = fopen(name, "w");
		fprintf(fp, "variables = x, y, cell, chainloopNO, chainloopINDEX, neighb, neighb1, neighb2, neighb3, neighb4, corner, corner1, corner2, corner3, corner4");
		fprintf(fp, "\r\nZONE T=\"tc[%.3f]\"\r\nN = %d, E = %d, DATAPACKING = BLOCK, ZONETYPE = FEQUADRILATERAL", tc, interfaceNO * 4, interfaceNO);
		fprintf(fp, "\r\nVARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED)");
		fprintf(fp, "\r\nSOLUTIONTIME = %.3e", tc);
		for (i = 0; i < interfaceNO; i++)
		{
			fprintf(fp, "\r\n%e", xyDelta[0][i] - xyDelta[2][i]);
			fprintf(fp, "\r\n%e", xyDelta[0][i] + xyDelta[2][i]);
			fprintf(fp, "\r\n%e", xyDelta[0][i] + xyDelta[2][i]);
			fprintf(fp, "\r\n%e", xyDelta[0][i] - xyDelta[2][i]);
		}
		for (i = 0; i < interfaceNO; i++)
		{
			fprintf(fp, "\r\n%e", xyDelta[1][i] - xyDelta[2][i]);
			fprintf(fp, "\r\n%e", xyDelta[1][i] - xyDelta[2][i]);
			fprintf(fp, "\r\n%e", xyDelta[1][i] + xyDelta[2][i]);
			fprintf(fp, "\r\n%e", xyDelta[1][i] + xyDelta[2][i]);
		}
		for (i = 0; i < interfaceNO; i++)
			fprintf(fp, "\r\n%d", i);
		for (i = 0; i < interfaceNO; i++)
			fprintf(fp, "\r\n%d", cellloopno[i]);
		for (i = 0; i < interfaceNO; i++)
			fprintf(fp, "\r\n%d", cellloopindex[i]);
		for (j = 0; j < 5; j++)
			for (i = 0; i < interfaceNO; i++)
				fprintf(fp, "\r\n%d", cellneighb[i][j]);
		for (j = 0; j < 5; j++)
			for (i = 0; i < interfaceNO; i++)
				fprintf(fp, "\r\n%d", cellcorner[i][j]);
		for (i = 0; i < interfaceNO; i++)
			fprintf(fp, "\r\n%d %d %d %d", i * 4 + 1, i * 4 + 2, i * 4 + 3, i * 4 + 4);
		fclose(fp);
		break;
	}
	};
	// free the memory
	;
	for (j = 0; j < interfaceNO; j++)
	{
		free(cellneighb[j]);
		free(cellcorner[j]);
		free(f1[j]);
	}
	free(cellneighb);
	free(cellcorner);
	free(f1);
	;
	free(cellloopindex);
	free(cellloopno);
	free(flag);
	free(CellNB);
	return true;
}

int removecell(int i, int **cellneighb, int **cellcorner, char *flag)
{
	int j, k, ij, nb;
	// neighbors
	for (j = 1; j <= cellneighb[i][0]; j++)
	{
		ij = cellneighb[i][j];
		nb = 0;
		for (k = 1; k <= cellneighb[ij][0]; k++)
		{
			if (i != cellneighb[ij][k])
				cellneighb[ij][++nb] = cellneighb[ij][k];
		}
		cellneighb[ij][0] = nb;
	}
	// corners
	for (j = 1; j <= cellcorner[i][0]; j++)
	{
		ij = cellcorner[i][j];
		nb = 0;
		for (k = 1; k <= cellcorner[ij][0]; k++)
		{
			if (i != cellcorner[ij][k])
				cellcorner[ij][++nb] = cellcorner[ij][k];
		}
		cellcorner[ij][0] = nb;
	}
	// excluding
	cellneighb[i][0] = 0;
	cellcorner[i][0] = 0;
	flag[i] = 'y';
	;
	return true;
}

int checkneighbcorner(int C1, int C2, int **cellNC)
{
	int i;
	for (i = 1; i <= cellNC[C1][0]; i++)
	{
		if (C2 == cellNC[C1][i])
		{
			return true;
		}
	}
	return false;
}

int mergeloops(int L1, int L2, int P1, int P2, int **ChainLoopIndex)
// L: loop index
// P: point index
// C: connection index
{
	int i, *sizenewloop, sizetmp;
	sizenewloop = (int *)calloc(ChainLoopIndex[L1][0] + ChainLoopIndex[L2][0] + 1, sizeof(int));
	sizetmp = 0;
	if (P1 > ChainLoopIndex[L1][0] / 2)
	{
		for (i = 1; i <= P1; i++)
			sizenewloop[++sizetmp] = ChainLoopIndex[L1][i];
	}
	else
	{
		for (i = ChainLoopIndex[L1][0]; i >= P1; i--)
			sizenewloop[++sizetmp] = ChainLoopIndex[L1][i];
	}
	if (P2 > ChainLoopIndex[L2][0] / 2)
	{
		for (i = P2; i >= 1; i--)
			sizenewloop[++sizetmp] = ChainLoopIndex[L2][i];
	}
	else
	{
		for (i = P2; i <= ChainLoopIndex[L2][0]; i++)
			sizenewloop[++sizetmp] = ChainLoopIndex[L2][i];
	}
	// reset loops L1 and L2
	for (i = 1; i <= ChainLoopIndex[L1][0]; i++)
		ChainLoopIndex[L1][i] = -1;
	ChainLoopIndex[L1][0] = 0;
	for (i = 1; i <= ChainLoopIndex[L2][0]; i++)
		ChainLoopIndex[L2][i] = -1;
	ChainLoopIndex[L2][0] = 0;
	// fill L1
	ChainLoopIndex[L1][0] = sizetmp;
	for (i = 1; i <= ChainLoopIndex[L1][0]; i++)
		ChainLoopIndex[L1][i] = sizenewloop[i];
	// free memory
	free(sizenewloop);
	return true;
}

int FindCellIsoLines(int cellNO, double **Var, const double VarRef, double **xyDelta, double **xp, double **yp, int *pno)
// Var[i][j]: "i" is cell-index, "j" is 0-1-2-3-4 which are the center, left-bottom, right-bottom, right-top, left-top
{
	int i, j, k;
	;
	// find iso-line points for a reference value
	;
	for (i = 0; i < cellNO; i++)
	{
		j = 1;
		k = 2;
		{
			if ((Var[i][j] > 0.50 && Var[i][k] <= 0.50) || (Var[i][j] < 0.50 && Var[i][k] >= 0.50))
			{
				xp[i][pno[i]] = xyDelta[0][i] + (2.0 * xyDelta[2][i] / (Var[i][k] - Var[i][j]) * (VarRef - Var[i][j]) - xyDelta[2][i]);
				yp[i][pno[i]] = xyDelta[1][i] - xyDelta[2][i];
				pno[i]++;
			}
		}
		j = 2;
		k = 3;
		{
			if ((Var[i][j] > 0.50 && Var[i][k] <= 0.50) || (Var[i][j] < 0.50 && Var[i][k] >= 0.50))
			{
				xp[i][pno[i]] = xyDelta[0][i] + xyDelta[2][i];
				yp[i][pno[i]] = xyDelta[1][i] - (2.0 * xyDelta[2][i] / (Var[i][k] - Var[i][j]) * (VarRef - Var[i][j]) - xyDelta[2][i]);
				pno[i]++;
			}
		}
		j = 3;
		k = 4;
		{
			if ((Var[i][j] > 0.50 && Var[i][k] <= 0.50) || (Var[i][j] < 0.50 && Var[i][k] >= 0.50))
			{
				xp[i][pno[i]] = xyDelta[0][i] - (2.0 * xyDelta[2][i] / (Var[i][k] - Var[i][j]) * (VarRef - Var[i][j]) - xyDelta[2][i]);
				yp[i][pno[i]] = xyDelta[1][i] + xyDelta[2][i];
				pno[i]++;
			}
		}
		j = 4;
		k = 1;
		{
			if ((Var[i][j] > 0.50 && Var[i][k] <= 0.50) || (Var[i][j] < 0.50 && Var[i][k] >= 0.50))
			{
				xp[i][pno[i]] = xyDelta[0][i] - xyDelta[2][i];
				yp[i][pno[i]] = xyDelta[1][i] + (2.0 * xyDelta[2][i] / (Var[i][k] - Var[i][j]) * (VarRef - Var[i][j]) - xyDelta[2][i]);
				pno[i]++;
			}
		}
	}
	return true;
}

int findlocalmaxminIndependent(double *var, int NO, int *index, char *type, int *noindex)
{
	int i, j, k, index2ndno, index3rdno, total, p, str;
	int *index1st, *index2nd, *index3rd;
	double minmaxtmp, maxall, threshold, *befr, *aftr, *both;
	char flag, *flag1st;
	;
	flag1st = (char *)calloc(NO, sizeof(char));
	index1st = (int *)calloc(NO, sizeof(int));
	index2nd = (int *)calloc(NO, sizeof(int));
	index3rd = (int *)calloc(NO, sizeof(int));
	befr = (double *)calloc(NO, sizeof(double));
	aftr = (double *)calloc(NO, sizeof(double));
	both = (double *)calloc(NO, sizeof(double));
	;
	// start threshold
	;
	threshold = -1.0;
	minmaxtmp = 0.0;
	for (i = 0; i < THRESHOLDINITIALGAP; i++)
		minmaxtmp += fabs(var[i]);
	if (threshold < (1.0 + THRESHOLDDROPK * 0.01) * minmaxtmp / (THRESHOLDINITIALGAP * 1.0))
		threshold = (1.0 + THRESHOLDDROPK * 0.01) * minmaxtmp / (THRESHOLDINITIALGAP * 1.0);
	;
	// check derivative signs
	total = 0;
	for (i = MARGINPOINTSSKIP; i < NO - MARGINPOINTSSKIP; i++)
	{
		if (fabs(var[i]) > threshold)
			break;
	}
	str = i;
	for (i = str; i < NO - MARGINPOINTSSKIP; i++)
	{
		flag1st[i] = 'z';
		k = 0;
		for (j = 1; j <= MARGINPOINTSSKIP; j++)
			if (var[i - (j - 1)] - var[i - j] > 0.0 && var[i + (j - 1)] - var[i + j] > 0.0)
				k++;
		j = MARGINPOINTSSKIP;
		if (k == MARGINPOINTSSKIP)
		{
			flag1st[i] = 'x';
			befr[i] = fabs(var[i] - var[i - j]);
			aftr[i] = fabs(var[i] - var[i + j]);
			index1st[total++] = i;
		}
		k = 0;
		for (j = 1; j <= MARGINPOINTSSKIP; j++)
			if (var[i - (j - 1)] - var[i - j] < 0.0 && var[i + (j - 1)] - var[i + j] < 0.0)
				k++;
		j = MARGINPOINTSSKIP;
		if (k == MARGINPOINTSSKIP)
		{
			flag1st[i] = 'n';
			befr[i] = fabs(var[i] - var[i - j]);
			aftr[i] = fabs(var[i] - var[i + j]);
			index1st[total++] = i;
		}
		for (j = 1; j <= NO; j++)
		{
			switch (flag1st[i])
			{
			case 'x':
			{
				if (var[i - (j - 1)] - var[i - j] >= 0.0)
				{
					if (befr[i] < fabs(var[i] - var[i - j]))
						befr[i] = fabs(var[i] - var[i - j]);
				}
				else
					j = NO + 1;
				break;
			}
			case 'n':
			{
				if (var[i - (j - 1)] - var[i - j] <= 0.0)
				{
					if (befr[i] < fabs(var[i] - var[i - j]))
						befr[i] = fabs(var[i] - var[i - j]);
				}
				else
					j = NO + 1;
				break;
			}
			default:
			{
				j = NO + 1;
				break;
			}
			}
		}
		for (j = 1; j <= NO; j++)
		{
			switch (flag1st[i])
			{
			case 'x':
			{
				if (var[i + (j - 1)] - var[i + j] >= 0.0)
				{
					if (aftr[i] < fabs(var[i] - var[i + j]))
						aftr[i] = fabs(var[i] - var[i + j]);
				}
				else
					j = NO + 1;
				break;
			}
			case 'n':
			{
				if (var[i + (j - 1)] - var[i + j] <= 0.0)
				{
					if (aftr[i] < fabs(var[i] - var[i + j]))
						aftr[i] = fabs(var[i] - var[i + j]);
				}
				else
					j = NO + 1;
				break;
			}
			default:
			{
				j = NO + 1;
				break;
			}
			}
		}
		//both[i] = befr[i];
		//if (both[i] < aftr[i])
		//	both[i] = aftr[i];
		//both[i] = 0.50 * (befr[i] + aftr[i]);
		both[i] = befr[i];
		if (both[i] > aftr[i])
			both[i] = aftr[i];
	};
	// filter the extremums to the second level (remove close extremums)
	;
	index2ndno = 0;
	index2nd[0] = index1st[0];
	for (i = 0; i < total; i++)
	{
		flag = 'n';
		p = index1st[i];
		//minmaxtmp = both[index1st[i]];
		minmaxtmp = fabs(var[index1st[i]]);
		for (k = index1st[i] + 1; k <= index1st[i] + MARGINPOINTS; k++)
		{
			for (j = 0; j < total; j++)
			{
				if (k == index1st[j])
				{
					//if (minmaxtmp < both[k])
					if (minmaxtmp < fabs(var[k]))
					{
						index2nd[index2ndno] = k;
						//minmaxtmp = both[k];
						minmaxtmp = fabs(var[k]);
						flag = 'y';
					}
					i = j;
				}
			}
		}
		if (flag == 'n')
			index2nd[index2ndno++] = p;
		else
			index2ndno++;
	};
	// global threshold
	;
	maxall = -1.0;
	for (i = 0; i < NO; i++)
	{
		if (maxall < fabs(var[i]))
			maxall = fabs(var[i]);
	}
	threshold = THRESHOLDMAXVALUE / 100.0 * maxall;
	;
	// threshold checking
	;
	index3rdno = 0;
	for (i = 0; i < index2ndno; i++)
	{
		if (fabs(var[index2nd[i]]) > threshold || fabs(both[index2nd[i]]) > threshold)
			index3rd[index3rdno++] = index2nd[i];
	};
	// sort by amplitude
	;
	//for (i = 0; i < index2ndno + 1; i++)
	//{
	//	for (j = 1; j < index2ndno; j++)
	//	{
	//		if (fabs(var[index2nd[j]]) > fabs(var[index2nd[j - 1]]))
	//		{
	//			k = index2nd[j];
	//			index2nd[j] = index2nd[j - 1];
	//			index2nd[j - 1] = k;
	//		}
	//	}
	//}
	;
	// sort by one/two side(s) (before-after) derivatives
	;
	//for (i = 0; i < index2ndno + 1; i++)
	//{
	//	for (j = 1; j < index2ndno; j++)
	//	{
	//		if (both[index2nd[j]] > both[index2nd[j - 1]])
	//		{
	//			k = index2nd[j];
	//			index2nd[j] = index2nd[j - 1];
	//			index2nd[j - 1] = k;
	//		}
	//	}
	//}
	;
	// sort by index
	;
	for (i = 0; i < index3rdno + 1; i++)
	{
		for (j = 1; j < index3rdno; j++)
		{
			if (index3rd[j] < index3rd[j - 1])
			{
				k = index3rd[j];
				index3rd[j] = index3rd[j - 1];
				index3rd[j - 1] = k;
			}
		}
	}
	// finishing
	for (i = 0; i < index3rdno; i++)
	{
		index[i] = index3rd[i];
		type[i] = flag1st[index3rd[i]];
	}
	*noindex = index3rdno;
	;
	free(flag1st);
	free(index1st);
	free(index2nd);
	free(index3rd);
	free(befr);
	free(aftr);
	free(both);
	return 1;
}

int lowpassfilter(double *varin, int NO, int order, double *varout)
{
	int i;
	double dr2nd[9] = {0.0, 0.0, 0.0, 0.50, 1.0, 0.50, 0.0, 0.0, 0.0};
	double dr4th[9] = {0.0, 0.0, -1.0 / 8.0, 1.0 / 2.0, 5.0 / 4.0, 1.0 / 2.0, -1.0 / 8.0, 0.0, 0.0};
	double dr6th[9] = {0.0, 1.0 / 32.0, -3.0 / 16.0, 15.0 / 32.0, 11.0 / 8.0, 15.0 / 32.0, -3.0 / 16.0, 1.0 / 32.0, 0.0};
	double dr8th[9] = {-1.0 / 128.0, 1.0 / 16.0, -7.0 / 32.0, 7.0 / 16.0, 93.0 / 64.0, 7.0 / 16.0, -7.0 / 32.0, 1.0 / 16.0, -1.0 / 128.0};
	for (i = 0; i < 9; i++)
	{
		dr2nd[i] *= 0.50;
		dr4th[i] *= 0.50;
		dr6th[i] *= 0.50;
		dr8th[i] *= 0.50;
	}
	switch (order)
	{
	case 2:
	{
		for (i = 1; i < NO - 2; i++)
			varout[i] = dr2nd[3] * varin[i - 1] + dr2nd[4] * varin[i] + dr2nd[5] * varin[i + 1];
		break;
	}
	case 4:
	{
		i = 1;
		varout[i] = dr2nd[3] * varin[i - 1] + dr2nd[4] * varin[i] + dr2nd[5] * varin[i + 1];
		i = NO - 2;
		varout[i] = dr2nd[3] * varin[i - 1] + dr2nd[4] * varin[i] + dr2nd[5] * varin[i + 1];
		;
		for (i = 2; i < NO - 3; i++)
			varout[i] = dr4th[2] * varin[i - 2] + dr4th[3] * varin[i - 1] + dr4th[4] * varin[i] + dr4th[5] * varin[i + 1] + dr4th[6] * varin[i + 2];
		break;
	}
	case 6:
	{
		i = 1;
		varout[i] = dr2nd[3] * varin[i - 1] + dr2nd[4] * varin[i] + dr2nd[5] * varin[i + 1];
		i = NO - 2;
		varout[i] = dr2nd[3] * varin[i - 1] + dr2nd[4] * varin[i] + dr2nd[5] * varin[i + 1];
		;
		i = 2;
		varout[i] = dr4th[2] * varin[i - 2] + dr4th[3] * varin[i - 1] + dr4th[4] * varin[i] + dr4th[5] * varin[i + 1] + dr4th[6] * varin[i + 2];
		i = NO - 3;
		varout[i] = dr4th[2] * varin[i - 2] + dr4th[3] * varin[i - 1] + dr4th[4] * varin[i] + dr4th[5] * varin[i + 1] + dr4th[6] * varin[i + 2];
		;
		for (i = 3; i < NO - 4; i++)
			varout[i] = dr6th[1] * varin[i - 3] + dr6th[2] * varin[i - 2] + dr6th[3] * varin[i - 1] + dr6th[4] * varin[i] + dr6th[5] * varin[i + 1] + dr6th[6] * varin[i + 2] + dr6th[7] * varin[i + 3];
		break;
	}
	case 8:
	{
		i = 1;
		varout[i] = dr2nd[3] * varin[i - 1] + dr2nd[4] * varin[i] + dr2nd[5] * varin[i + 1];
		i = NO - 2;
		varout[i] = dr2nd[3] * varin[i - 1] + dr2nd[4] * varin[i] + dr2nd[5] * varin[i + 1];
		;
		i = 2;
		varout[i] = dr4th[2] * varin[i - 2] + dr4th[3] * varin[i - 1] + dr4th[4] * varin[i] + dr4th[5] * varin[i + 1] + dr4th[6] * varin[i + 2];
		i = NO - 3;
		varout[i] = dr4th[2] * varin[i - 2] + dr4th[3] * varin[i - 1] + dr4th[4] * varin[i] + dr4th[5] * varin[i + 1] + dr4th[6] * varin[i + 2];
		;
		i = 3;
		varout[i] = dr6th[1] * varin[i - 3] + dr6th[2] * varin[i - 2] + dr6th[3] * varin[i - 1] + dr6th[4] * varin[i] + dr6th[5] * varin[i + 1] + dr6th[6] * varin[i + 2] + dr6th[7] * varin[i + 3];
		i = NO - 4;
		varout[i] = dr6th[1] * varin[i - 3] + dr6th[2] * varin[i - 2] + dr6th[3] * varin[i - 1] + dr6th[4] * varin[i] + dr6th[5] * varin[i + 1] + dr6th[6] * varin[i + 2] + dr6th[7] * varin[i + 3];
		;
		for (i = 4; i < NO - 5; i++)
			varout[i] = dr8th[0] * varin[i - 4] + dr8th[1] * varin[i - 3] + dr8th[2] * varin[i - 2] + dr8th[3] * varin[i - 1] + dr8th[4] * varin[i] + dr8th[5] * varin[i + 1] + dr8th[6] * varin[i + 2] + dr8th[7] * varin[i + 3] + dr8th[8] * varin[i + 4];
		break;
	}
	}
	return 1;
}

int averagesmooth(double *varin, int NO, int order, double *varout)
{
	int i, j;
	for (i = order; i < NO - order; i++)
	{
		varout[i] = varin[i];
		for (j = 0; j < order; j++)
			varout[i] += varin[i - j] + varin[i + j];
		varout[i] /= (2.0 * order + 1.0);
	}
	for (i = 0; i < order; i++)
		varout[i] = varout[order];
	for (i = NO - order; i < NO; i++)
		varout[i] = varout[NO - order - 1];
	return 1;
}

int WriteTmprlData(char *name, int Iteration, double tc, double **xyDKO, int *JetRoots)
{
	double angle, dx, dy;
	FILE *fpout;
	switch (Iteration)
	{
	case 0:
	{
		fpout = fopen(name, "w");
		fprintf(fpout, "variables = time, angle, xT, yT, KappaT, OmegaT, xB, yB, KappaB, OmegaB");
		fprintf(fpout, "\r\nzone T = \"Jet Data\"\r\n");
		break;
	}
	default:
	{
		fpout = fopen(name, "a");
		break;
	}
	}
	if (JetRoots[0] == -1 || JetRoots[1] == -1)
	{
		fclose(fpout);
		return 1;
	}
	dx = xyDKO[0][JetRoots[0]] - xyDKO[0][JetRoots[1]];
	dy = xyDKO[1][JetRoots[1]] - xyDKO[1][JetRoots[0]];
	angle = atan2(dy, dx) * 180.0 / 3.1415926535897932384626433832795;
	if (angle < 0.0)
		angle += 360.0;
	fprintf(fpout, "%e %e ", tc, angle);
	fprintf(fpout, "%e %e ", xyDKO[0][JetRoots[0]], xyDKO[1][JetRoots[0]]);
	fprintf(fpout, "%e %e ", xyDKO[3][JetRoots[0]], xyDKO[4][JetRoots[0]]);
	fprintf(fpout, "%e %e ", xyDKO[0][JetRoots[1]], xyDKO[1][JetRoots[1]]);
	fprintf(fpout, "%e %e\r\n", xyDKO[3][JetRoots[1]], xyDKO[4][JetRoots[1]]);
	fclose(fpout);
	;
	return 1;
}

int Write2DJetRoot(char *name, int Iteration, double tc, double **xyDKO, double *kppF, double *omgF, int CellNO, int *Kindex, int KindexNO, int *Oindex, int OindexNO, int **JetRoots)
{
	int i, k, f;
	FILE *fp;
	switch (Iteration)
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
	// original values: all cells
	fprintf(fp, "variables = x, y, Delta, Kappa, Omega, KappaFiltered, OmegaFiltered");
	fprintf(fp, "\r\nZONE T=\"AllData[%.4f]\"\r\n"
				"N = %d, E = %d, DATAPACKING = BLOCK, ZONETYPE = FEQUADRILATERAL",
			tc, CellNO * 4, CellNO);
	fprintf(fp, "\r\nVARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED)");
	fprintf(fp, "\r\nSOLUTIONTIME = %.4e", tc);
	for (i = 0; i < CellNO; i++)
	{
		fprintf(fp, "\r\n%e", xyDKO[0][i] - xyDKO[2][i]);
		fprintf(fp, "\r\n%e", xyDKO[0][i] + xyDKO[2][i]);
		fprintf(fp, "\r\n%e", xyDKO[0][i] + xyDKO[2][i]);
		fprintf(fp, "\r\n%e", xyDKO[0][i] - xyDKO[2][i]);
	}
	for (i = 0; i < CellNO; i++)
	{
		fprintf(fp, "\r\n%e", xyDKO[1][i] - xyDKO[2][i]);
		fprintf(fp, "\r\n%e", xyDKO[1][i] - xyDKO[2][i]);
		fprintf(fp, "\r\n%e", xyDKO[1][i] + xyDKO[2][i]);
		fprintf(fp, "\r\n%e", xyDKO[1][i] + xyDKO[2][i]);
	}
	for (i = 0; i < CellNO; i++)
		fprintf(fp, "\r\n%e", xyDKO[2][i]);
	for (i = 0; i < CellNO; i++)
		fprintf(fp, "\r\n%e", xyDKO[3][i]);
	for (i = 0; i < CellNO; i++)
		fprintf(fp, "\r\n%e", xyDKO[4][i]);
	for (i = 0; i < CellNO; i++)
		fprintf(fp, "\r\n%e", kppF[i]);
	for (i = 0; i < CellNO; i++)
		fprintf(fp, "\r\n%e", omgF[i]);
	for (i = 0; i < CellNO; i++)
		fprintf(fp, "\r\n%d %d %d %d", i * 4 + 1, i * 4 + 2, i * 4 + 3, i * 4 + 4);
	for (f = 1; f <= JetRoots[0][0]; f++)
	{
		// jet roots
		fprintf(fp, "variables = x, y, Delta, Kappa, Omega, KappaFiltered, OmegaFiltered");
		fprintf(fp, "\r\nZONE T=\"JetRoot[%d][%.4f]\"\r\n"
					"N = %d, E = %d, DATAPACKING = BLOCK, ZONETYPE = FEQUADRILATERAL",
				f, tc, 2 * 4, 2);
		fprintf(fp, "\r\nVARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED)");
		fprintf(fp, "\r\nSOLUTIONTIME = %.4e", tc);
		if (JetRoots[f][0] == -1 || JetRoots[f][1] == -1)
		{
			for (i = 0; i < 2; i++)
			{
				fprintf(fp, "\r\n%e", xyDKO[0][i] - xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[0][i] + xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[0][i] + xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[0][i] - xyDKO[2][i]);
			}
			for (i = 0; i < 2; i++)
			{
				fprintf(fp, "\r\n%e", xyDKO[1][i] - xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[1][i] - xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[1][i] + xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[1][i] + xyDKO[2][i]);
			}
			for (i = 0; i < 2; i++)
				fprintf(fp, "\r\n%e", xyDKO[2][i]);
			for (i = 0; i < 2; i++)
				fprintf(fp, "\r\n%e", xyDKO[3][i]);
			for (i = 0; i < 2; i++)
				fprintf(fp, "\r\n%e", xyDKO[4][i]);
			for (i = 0; i < 2; i++)
				fprintf(fp, "\r\n%e", kppF[i]);
			for (i = 0; i < 2; i++)
				fprintf(fp, "\r\n%e", omgF[i]);
		}
		else
		{
			for (k = 0; k < 2; k++)
			{
				i = JetRoots[f][k];
				fprintf(fp, "\r\n%e", xyDKO[0][i] - xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[0][i] + xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[0][i] + xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[0][i] - xyDKO[2][i]);
			}
			for (k = 0; k < 2; k++)
			{
				i = JetRoots[f][k];
				fprintf(fp, "\r\n%e", xyDKO[1][i] - xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[1][i] - xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[1][i] + xyDKO[2][i]);
				fprintf(fp, "\r\n%e", xyDKO[1][i] + xyDKO[2][i]);
			}
			for (k = 0; k < 2; k++)
			{
				i = JetRoots[f][k];
				fprintf(fp, "\r\n%e", xyDKO[2][i]);
			}
			for (k = 0; k < 2; k++)
			{
				i = JetRoots[f][k];
				fprintf(fp, "\r\n%e", xyDKO[3][i]);
			}
			for (k = 0; k < 2; k++)
			{
				i = JetRoots[f][k];
				fprintf(fp, "\r\n%e", xyDKO[4][i]);
			}
			for (k = 0; k < 2; k++)
			{
				i = JetRoots[f][k];
				fprintf(fp, "\r\n%e", kppF[i]);
			}
			for (k = 0; k < 2; k++)
			{
				i = JetRoots[f][k];
				fprintf(fp, "\r\n%e", omgF[i]);
			}
		}
		for (i = 0; i < 2; i++)
			fprintf(fp, "\r\n%d %d %d %d", i * 4 + 1, i * 4 + 2, i * 4 + 3, i * 4 + 4);
	}
	fclose(fp);
	return true;
}

int WriteExtremums(char *name, int Iteration, double SolutionTime, double **xyDKO, double *kppF, double *omgF, int CellNO, int *Kindex, int KindexNO, int *Oindex, int OindexNO, int **JetRoots)
{
	int i, k, f;
	FILE *fpout;
	switch (Iteration)
	{
	case 0:
	{
		fpout = fopen(name, "w");
		break;
	}
	default:
	{
		fpout = fopen(name, "a");
		break;
	}
	}
	// original values: kappa
	fprintf(fpout, "variables = index, kappa");
	fprintf(fpout, "\r\nzone T = \"Original Kappa\"\r\nSOLUTIONTIME = %.4f\r\n", SolutionTime);
	for (i = 0; i < CellNO; i++)
	{
		fprintf(fpout, "%d ", i);
		fprintf(fpout, "%e ", xyDKO[3][i]);
		fprintf(fpout, "\r\n");
	}
	// filtered values: kappa
	fprintf(fpout, "variables = index, kappa");
	fprintf(fpout, "\r\nzone T = \"Filtered Kappa\"\r\nSOLUTIONTIME = %.4f\r\n", SolutionTime);
	for (i = 0; i < CellNO; i++)
	{
		fprintf(fpout, "%d ", i);
		fprintf(fpout, "%e ", kppF[i]);
		fprintf(fpout, "\r\n");
	}
	// original values: omega
	fprintf(fpout, "variables = index, omega");
	fprintf(fpout, "\r\nzone T = \"Original Omega\"\r\nSOLUTIONTIME = %.4f\r\n", SolutionTime);
	for (i = 0; i < CellNO; i++)
	{
		fprintf(fpout, "%d ", i);
		fprintf(fpout, "%e ", xyDKO[4][i]);
		fprintf(fpout, "\r\n");
	}
	// filtered values: omega
	fprintf(fpout, "variables = index, omega");
	fprintf(fpout, "\r\nzone T = \"Filtered Omega\"\r\nSOLUTIONTIME = %.4f\r\n", SolutionTime);
	for (i = 0; i < CellNO; i++)
	{
		fprintf(fpout, "%d ", i);
		fprintf(fpout, "%e ", omgF[i]);
		fprintf(fpout, "\r\n");
	}
	// extremums: kappa
	fprintf(fpout, "variables = index, kappa");
	fprintf(fpout, "\r\nzone T = \"Extremums Kappa\"\r\nSOLUTIONTIME = %.4f\r\n", SolutionTime);
	switch (KindexNO)
	{
	case 0:
	{
		i = 0;
		fprintf(fpout, "%d ", i);
		fprintf(fpout, "%e ", kppF[i]);
		fprintf(fpout, "\r\n");
		break;
	}
	default:
	{
		for (k = 0; k < KindexNO; k++)
		{
			i = Kindex[k];
			fprintf(fpout, "%d ", i);
			fprintf(fpout, "%e ", kppF[i]);
			fprintf(fpout, "\r\n");
		}
		break;
	}
	}
	// extremums: omega
	fprintf(fpout, "variables = index, omega");
	fprintf(fpout, "\r\nzone T = \"Extremums Omega\"\r\nSOLUTIONTIME = %.4f\r\n", SolutionTime);
	switch (OindexNO)
	{
	case 0:
	{
		i = 0;
		fprintf(fpout, "%d ", i);
		fprintf(fpout, "%e ", omgF[i]);
		fprintf(fpout, "\r\n");
		break;
	}
	default:
	{
		for (k = 0; k < OindexNO; k++)
		{
			i = Oindex[k];
			fprintf(fpout, "%d ", i);
			fprintf(fpout, "%e ", omgF[i]);
			fprintf(fpout, "\r\n");
		}
		break;
	}
	}
	for (f = 1; f <= JetRoots[0][0]; f++)
	{
		// jet values: kappa
		fprintf(fpout, "variables = index, kappa");
		fprintf(fpout, "\r\nzone T = \"Jet Kappa [%d]\"\r\nSOLUTIONTIME = %.4f\r\n", f, SolutionTime);
		if (JetRoots[f][0] == -1 || JetRoots[f][1] == -1)
		{
			i = 0;
			{
				fprintf(fpout, "%d ", i);
				fprintf(fpout, "%e ", kppF[i]);
				fprintf(fpout, "\r\n");
			}
			i = 0;
			{
				fprintf(fpout, "%d ", i);
				fprintf(fpout, "%e ", kppF[i]);
				fprintf(fpout, "\r\n");
			}
		}
		else
		{
			i = JetRoots[f][0];
			{
				fprintf(fpout, "%d ", i);
				fprintf(fpout, "%e ", kppF[i]);
				fprintf(fpout, "\r\n");
			}
			i = JetRoots[f][1];
			{
				fprintf(fpout, "%d ", i);
				fprintf(fpout, "%e ", kppF[i]);
				fprintf(fpout, "\r\n");
			}
		}
		// jet values: omega
		fprintf(fpout, "variables = index, omega");
		fprintf(fpout, "\r\nzone T = \"Jet Omega [%d]\"\r\nSOLUTIONTIME = %.4f\r\n", f, SolutionTime);
		if (JetRoots[f][0] == -1 || JetRoots[f][1] == -1)
		{
			i = 0;
			{
				fprintf(fpout, "%d ", i);
				fprintf(fpout, "%e ", omgF[i]);
				fprintf(fpout, "\r\n");
			}
			i = 0;
			{
				fprintf(fpout, "%d ", i);
				fprintf(fpout, "%e ", omgF[i]);
				fprintf(fpout, "\r\n");
			}
		}
		else
		{
			i = JetRoots[f][0];
			{
				fprintf(fpout, "%d ", i);
				fprintf(fpout, "%e ", omgF[i]);
				fprintf(fpout, "\r\n");
			}
			i = JetRoots[f][1];
			{
				fprintf(fpout, "%d ", i);
				fprintf(fpout, "%e ", omgF[i]);
				fprintf(fpout, "\r\n");
			}
		}
	}
	fclose(fpout);
	return 1;
}

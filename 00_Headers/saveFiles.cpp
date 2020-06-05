//#ifndef SAVEFILES_H
//#define SAVEFILES_H

using namespace std;

//******************************************************************//
//Mindenfele mellekinfo kimento fv-ek
//******************************************************************//
//list<string> kimentese, mely vegyes informaciokat tartalmaz, 2 soronkent: egy sor leiras, egy sor adat (utana egy ures sor kerul kiiratasra)
//Használat pl: list<string> infos; infos.push_back("Energia minimum:"); infos.push_back("22"); saveLS_infos(infos,"Egyebek");
void saveLS_infos(list<string> infos, string fileName)
{
	fileName += ".inf";
	int num = infos.size();

	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	
	for(int i=0; i<num; i=i+2)
	{
			save << infos.front() << "\t";
			infos.pop_front();
			save << infos.front() << endl;
			infos.pop_front();
	}

	save.close();	
}

//******************************************************************//
// Az alabbi a fuggvenyek adatokat mentenek ki ugy, hogy azok alkalmasak legyenek GNUplot-ban valo abrazolasra.
// Az egyszeru adatok [.dat] kiterjesztest kapnak.
// Meg kell adni a a file nevet [fileName], tovabba az elso sorba plussz informacio adhato, mely #-tel lesz gnuplot komment. 
//******************************************************************//
// double-tipusu vector kimentese
void saveVecDoub(vector<double> V, string fileName, string VecInfo)
{
	VecInfo = (VecInfo == "") ? "# No additional information about data" : "#" + VecInfo;
	
	fileName += ".dat";
	int num = V.size();

	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << VecInfo << endl;

	for(int i=0; i<num; i++)
			save << i << '\t' << V[i] << endl;

	save.close();	
}

//int-tipusu vector kimentese
void saveVecInt(vector<int> V, string fileName, string VecInfo)
{
	VecInfo = (VecInfo == "") ? "# No additional information about data" : "#" + VecInfo;
	
	fileName += ".dat";
	int num = V.size();

	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << VecInfo << endl;

	for(int i=0; i<num; i++)
			save << i << '\t' << V[i] << endl;

	save.close();	
}

//bool-tipusu vector kimentese
void saveVecBool(vector<bool> V, string fileName, string VecInfo)
{
	VecInfo = (VecInfo == "") ? "# No additional information about data" : "#i	" + VecInfo;
	
	fileName += ".dat";
	int num = V.size();

	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << VecInfo << endl;

	for(int i=0; i<num; i++)
			save << i << '\t' << V[i] << endl;

	save.close();	
}

// double-tipusu matrix kimentese
void saveMatDoub(vector<vector<double> >  M, string fileName, string firstLine)
{
	firstLine = (firstLine == "") ? "# No additional information about data" : "#" + firstLine;
	fileName += ".dat";
	int row = M[0].size(); int col = M.size();
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
			save << M[j][i] << '\t';

		save << endl;
	}
	save.close();	
}

// double-tipusu NEGYZETES matrix kimentese, grid-elt 3D-s GNUplot abrahoz.
void saveMatDoub3d(vector<vector<double> >  M, string fileName, string firstLine)
{
	firstLine = (firstLine == "") ? "# No additional information about data" : "#" + firstLine;
	fileName += ".dat";
	int N = M.size();
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			save << i << '\t' << j << '\t' << M[i][j] << endl;
		}
		save << endl;
	}	
	save.close();
}

// double-tipusu matrix kimentese, grid-elt 3D-s GNUplot abrahoz, allithato tengelyekkel
void saveMatDoub3d_wAxes(vector<vector<double> >  M, double x_min, double x_max, double dx, double y_min, double y_max, double dy, string fileName, string firstLine)
{
	firstLine = (firstLine == "") ? "# No additional information about data" : "#" + firstLine;
	fileName += ".dat";
	//int N = M.size();
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;
	int i = -1;
	int j = -1;

	for(double x=x_min; x<=x_max; x+=dx)
	{
		i++; j=-1;
		for(double y=y_min; y<=y_max; y+=dy)
		{
			j++;
			save << x << '\t' << y << '\t' << M[i][j] << endl;
		}
		save << endl;
	}	
	save.close();
}

// double-tipusu matrixot vektorkent ment ki fileba
void saveMatDoubToVec(vector<vector<double> >  M, string fileName, string firstLine)
{
	firstLine = (firstLine == "") ? "# No additional information about data" : "#" + firstLine;
	fileName += ".dat";
	int row = M[0].size(); int col = M.size();
	int count = 0;
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
		{
			save << count << '\t' << M[j][i] << endl;
			count++;
		}
	}
	save.close();	
}

// int-tipusu matrix kimentese
void saveMatInt(vector<vector<int> >  M, string fileName, string firstLine)
{
	firstLine = (firstLine == "") ? "# No additional information about data" : "#" + firstLine;
	fileName += ".dat";
	int row = M[0].size(); int col = M.size();
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
			save << M[j][i] << '\t';

		save << endl;
	}
	save.close();	
}

// int-tipusu NEGYZETES matrix kimentese, grid-elt 3D-s GNUplot abrahoz.
void saveMatInt3d(vector<vector<int> >  M, string fileName, string firstLine)
{
	firstLine = (firstLine == "") ? "# No additional information about data" : "#" + firstLine;
	fileName += ".dat";
	int N = M.size();
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			save << i << '\t' << j << '\t' << M[i][j] << endl;
		}
		save << endl;
	}	
	save.close();
}

// bool-tipusu matrix kimentese
void saveMatBool(vector<vector<bool> >  M, string fileName, string firstLine)
{
	firstLine = (firstLine == "") ? "# True = 1, False = 0" : "# True = 1, False = 0" + firstLine;
	fileName += ".dat";
	int row = M[0].size(); int col = M.size();
	int tmp;
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
		{
			tmp = (M[j][i] == true) ? 1 : 0;
			save << tmp << '\t';
		}

		save << endl;
	}
	save.close();	
}

// bool-tipusu NEGYZETES matrix kimentese, grid-elt 3D-s GNUplot abrahoz.
void saveMatBool3d(vector<vector<bool> >  M, string fileName, string firstLine)
{
	firstLine = (firstLine == "") ? "# True = 1, False = 0" : "# True = 1, False = 0" + firstLine;
	fileName += ".dat";
	int N = M.size();
	int tmp;
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			tmp = (M[i][j] == true) ? 1 : 0;
			save << i << '\t' << j << '\t' << tmp << endl;
		}
		save << endl;
	}	
	save.close();
}

//******************************************************************//
//Az alabbi fuggvenyek PBM/PGM file-okat mentenek ki, ahol +1/true = feher (0), -1/false = fekete (1/255)
//******************************************************************//
// bool-tipusu vector-t X-Y oldalu racsba ment
void saveVecBool_meshPBM(vector<bool> V, int X, int Y, int pixelSize, string fileName)
{
	fileName += ".pbm";
	int count = 0;
	int countTmp = 0;
	int tmp;
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << "P1" << endl;
	save << "# " << fileName << endl;
	save << X*pixelSize << "\t" << Y*pixelSize << endl;

	for(int i=0; i<Y; i++)
	{
		for (int a=0; a<pixelSize; a++)
		{
			countTmp = count;
			for(int j=0; j<X; j++)
			{
				tmp = (V[countTmp] == true) ? 0 : 1;
				
				for (int b=0; b<pixelSize; b++)
					save << tmp << '\t';
				
				countTmp++;
			}
			save << endl;
		}
		count+=X;
	}
	save.close();	
}

// bool-tipusu matrixot ment ki PBM-be
void saveMatBool_meshPBM(vector<vector<bool> > M, int pixelSize, string fileName)
{
	fileName += ".pbm";
	int col = M[0].size();
	int row = M.size();
	int tmp;
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << "P1" << endl;
	save << "# " << fileName << endl;
	save << col*pixelSize << "\t" << row*pixelSize << endl;

	for (int i=0; i<row; i++)
	{
		for (int a=0; a<pixelSize; a++)
		{
			for (int j=0; j<col; j++)
			{	
				tmp = (M[i][j] == true) ? 0 : 1;

				for (int b=0; b<pixelSize; b++)
					save << tmp << '\t';			
			}
			save << endl;
		}
	}
	save.close();	
}

// int-tipusu matrixot ment ki PGM-be; max 255 lehet az ertek! (Pl: compareMat matrix kimentese)
void saveMatInt_meshPGM(vector<vector<int> > M, string fileName, int pixelSize)
{
	fileName += ".pgm";
	int col = M[0].size();
	int row = M.size();

	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << "P2" << endl << col * pixelSize << " " << row * pixelSize << endl << "255" << endl;

	for (int i=0; i<row; i++)
	{
		for (int a=0; a<pixelSize; a++)
		{
			for (int j=0; j<col; j++)
			{	
				for (int b=0; b<pixelSize; b++)
				{
					save << M[i][j] << '\t';			
				}
			}
			save << endl;
		}
	}
	save.close();	
}

void saveMatInt_meshPGM_forJ(vector<vector<int> > J, string fileName, int pixelSize)
{
	fileName += ".pgm";
	int col = J[0].size();
	int row = J.size();

	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << "P2" << endl << col * pixelSize << " " << row * pixelSize << endl << "255" << endl;

	for (int i=0; i<row; i++)
	{
		for (int a=0; a<pixelSize; a++)
		{
			for (int j=0; j<col; j++)
			{	
				for (int b=0; b<pixelSize; b++)
				{
					if (J[i][j] == -1)
					{
						save << "0" << '\t';
					}
					else if (J[i][j] == 1)
					{
						save << "255" << '\t';
					}
					else
					{
						save << "127" << '\t';
					}
				}
			}
			save << endl;
		}
	}
	save.close();	
}

//******************************************************************//
//Az alabbi fuggvenyek Valoszinuseg-suruseg eloszlas-fv.t (PDDF) készítenek egy adatsorból
//******************************************************************//
// double-tipusu vectorbol keszit Valoszinuseg-suruseg eloszlas-fv-t; 0-n van bin
void savePDDF_vec(vector<double> datas, double binSize, double minBorder, double maxBorder, string fileName, string firstLine)
{
	int dataNum = datas.size();
	double min = datas[0]; double max = datas[0];
	double tmp;
		
	if (minBorder == 0 && maxBorder == 0)		//be lehet állítani egy min és max értéket. Ha 0-0-t adunk meg, akkor ez megkeresi a globális min és max-ot határoknak.
	{
		//min es max keresese
		for (int i=1; i<dataNum; i++)
		{
			if (datas[i] > max) max = datas[i];
			if (datas[i] < min) min = datas[i];
		}
	}
	else
	{
		min = minBorder;
		max = maxBorder;
	}

	int distSize = int(floor(((max-min)/binSize)+1));		//ennyi elemu lesz az eloszlas
	vector<double> distrPD(distSize,0.0);

	for (int i=0; i<dataNum; i++)
	{
		tmp = min-(binSize/2);
		for (int j=0; j<distSize; j++)
		{
			if (datas[i] >= tmp && datas[i] <= tmp+binSize)
			{
				distrPD[j]++;
				break;
			}
			else
			{
				tmp += binSize;
			}
		}
	}

	//Normalas
	for (int i=0; i<distSize; i++)
		distrPD[i] /= dataNum;			//igy az integral 1 lesz...

	//A szeleken fel-binek vannak
	distrPD[0] *= 2;
	distrPD[distSize-1] *= 2;

	//kiiratas:
	if (firstLine == "") firstLine = "# No additional information about data";
	
	fileName += ".dat";
	firstLine = "#"+firstLine;
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<distSize; i++)
			save << min+(i*binSize) << '\t' << distrPD[i] << endl;

	save.close();	
	
}

// double-tipusu vectorbol keszit Valoszinuseg-suruseg eloszlas-fv-t; 0-mellett vannak a bin-ek
void savePDDF_vec_alter(vector<double> datas, double binSize, double minBorder, double maxBorder, string fileName, string firstLine)		
{
	int dataNum = datas.size();
	double min = datas[0]; double max = datas[0];
	double tmp;
		
	if (minBorder == 0 && maxBorder == 0)		//be lehet állítani egy min és max értéket. Ha 0-0-t adunk meg, akkor ez megkeresi a globális min és max-ot határoknak.
	{
		//min es max keresese
		for (int i=1; i<dataNum; i++)
		{
			if (datas[i] > max) max = datas[i];
			if (datas[i] < min) min = datas[i];
		}
	}
	else
	{
		min = minBorder;
		max = maxBorder;
	}

	int distSize = int(floor((max-min)/binSize));		//ennyi elemu lesz az eloszlas
	vector<double> distrPD(distSize,0.0);

	for (int i=0; i<dataNum; i++)
	{
		tmp = min;
		for (int j=0; j<distSize; j++)
		{
			if (datas[i] >= tmp && datas[i] <= tmp+binSize)
			{
				distrPD[j]++;
				break;
			}
			else
			{
				tmp += binSize;
			}
		}
	}

	//Normalas
	for (int i=0; i<distSize; i++)
		distrPD[i] /= dataNum;			//igy az integral 1 lesz...

	//kiiratas:
	if (firstLine == "") firstLine = "# No additional information about data";
	
	fileName += ".dat";
	firstLine = "#"+firstLine;
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<distSize; i++)
			save << min+(binSize/2.0)+(i*binSize) << '\t' << distrPD[i] << endl;

	save.close();	
	
}

// double-tipusu matrixbol keszit Valoszinuseg-suruseg eloszlas-fv-t; 0-n van bin
void savePDDF_mat(vector<vector<double> >  datas, int colNum, double binSize, double minBorder, double maxBorder, string fileName, string firstLine)		
	//ha nem allitunk be also és felso hatart, akkor megkeresi a min es max ertekeket hataroknak
	//colNum: ebben az oszlopban levo adatokbol szamolja a PDDF-et (mint a savePDDF_vec!)
{
	int dataNum = datas[0].size();
	double min = datas[colNum][0]; double max = datas[colNum][0];
	double tmp;
		
	if (minBorder == 0 && maxBorder == 0)
	{
		//min es max keresese
		for (int i=1; i<dataNum; i++)
		{
			if (datas[colNum][i] > max) max = datas[colNum][i];
			if (datas[colNum][i] < min) min = datas[colNum][i];
		}
	}
	else
	{
		min = minBorder;
		max = maxBorder;
	}

	int distSize = int(floor(((max-min)/binSize)+1));
	vector<double> distrPD(distSize,0.0);

	for (int i=0; i<dataNum; i++)
	{
		tmp = min-(binSize/2);
		for (int j=0; j<distSize; j++)
		{
			if (datas[colNum][i] >= tmp && datas[colNum][i] <= tmp+binSize)
			{
				distrPD[j]++;
				break;
			}
			else
			{
				tmp += binSize;
			}
		}
	}

	//Normalas
	for (int i=0; i<distSize; i++)
		distrPD[i] /= double(dataNum);			//igy az integral 1 lesz...

	//A szeleken fel-binek vannak
	distrPD[0] *= 2;
	distrPD[distSize-1] *= 2;

	if (firstLine == "") firstLine = "# No additional information about data";
	
	fileName += ".dat";
	firstLine = "#"+firstLine;
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<distSize; i++)
			save << min+(i*binSize) << '\t' << distrPD[i] << endl;

	save.close();	
	
}

// double-tipusu matrixbol keszit Valoszinuseg-suruseg eloszlas-fv-t; 0 mellett vannak a binek
void savePDDF_mat_alter(vector<vector<double> >  datas, int colNum, double binSize, double minBorder, double maxBorder, string fileName, string firstLine)		
	//ha nem allitunk be also és felso hatart, akkor megkeresi a min es max ertekeket hataroknak
	//colNum: ebben az oszlopban levo adatokbol szamolja a PDDF-et (mint a savePDDF_vec!)
{
	int dataNum = datas.size();
	double min = datas[colNum][0]; double max = datas[colNum][0];
	double tmp;
		
	if (minBorder == 0 && maxBorder == 0)
	{
		//min es max keresese
		for (int i=1; i<dataNum; i++)
		{
			if (datas[colNum][i] > max) max = datas[colNum][i];
			if (datas[colNum][i] < min) min = datas[colNum][i];
		}
	}
	else
	{
		min = minBorder;
		max = maxBorder;
	}

	int distSize = int(floor((max-min)/binSize));		//ennyi elemu lesz az eloszlas
	vector<double> distrPD(distSize,0.0);

	for (int i=0; i<dataNum; i++)
	{
		tmp = min;
		for (int j=0; j<distSize; j++)
		{
			if (datas[colNum][i] >= tmp && datas[colNum][i] <= tmp+binSize)
			{
				distrPD[j]++;
				break;
			}
			else
			{
				tmp += binSize;
			}
		}
	}

	//Normalas
	for (int i=0; i<distSize; i++)
		distrPD[i] /= double(dataNum);			//igy az integral 1 lesz...

	if (firstLine == "") firstLine = "# No additional information about data";
	
	fileName += ".dat";
	firstLine = "#"+firstLine;
	
	ofstream save;
	const char *file_char;
	file_char = fileName.c_str();
	save.open(file_char);
	save << firstLine << endl;

	for(int i=0; i<distSize; i++)
			save << min+(i*binSize) << '\t' << distrPD[i] << endl;

	save.close();	
	
}

//#endif
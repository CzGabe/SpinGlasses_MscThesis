//#ifndef SAVEFILES_H
//#define SAVEFILES_H

using namespace std;

//int tipusu rnd CompleteGraph J matrixot letrehozo fuggveny
vector<vector<int> >  create_rndJint_CG (int N, MTRand& MTrnd)
{
	vector<int> VecZeros(N,0);
	vector<vector<int> >  matrix(N,VecZeros);		//A fo-atloban vegig "0" van...
			
	for (int i=0; i<N; i++)
	{
		for (int j=i; j<N; j++)
		{
			if (i == j)
				matrix[i][j] = 0;
			else
				matrix[i][j] = matrix [j][i] = (MTrnd() < 0.5) ? +1 : -1;
		}
	}
	
	return matrix;
}

///bool tipusu rnd spin-vektort letrehozo fuggveny
vector<bool>  create_rndSpinBool_CG (int N, MTRand& MTrnd)
{
	vector<bool>  vecSpin(N,false);
		
	for (int i=0; i<N; i++)
		vecSpin[i] = (MTrnd() < 0.5) ? true : false;

	return vecSpin;
}

//******************************************************************//
//string konverterek
//******************************************************************//
string convert_DtoS(double number)
{
   stringstream ss;		
   ss << number;		
   return ss.str();		
}

string convert_ItoS(int number)
{
   stringstream ss;		
   ss << number;		
   return ss.str();		
}

string convert_DtoFix2S(double number)
{
   stringstream ss;		
   ss << fixed << setprecision (2) << number;		
   return ss.str();		
}

//3 -> "003"; 45 -> "045"; 123 -> "123"
string convert_Ito3str(int number)
{
	string strNum = "";

	if(number<10)
	{
		strNum = "00" + convert_ItoS(number);
	}
	else if(number<100)
	{
		strNum = "0" + convert_ItoS(number);
	}
	else
	{
		strNum = convert_ItoS(number);
	}

	return strNum;
}

//3 -> "00003"; 45 -> "00045"; 123 -> "00123"
string convert_Ito5str(int number)
{
	string strNum = "";

	if(number<10)
	{
		strNum = "0000" + convert_ItoS(number);
	}
	else if(number<100)
	{
		strNum = "000" + convert_ItoS(number);
	}
	else if(number<1000)
	{
		strNum = "00" + convert_ItoS(number);
	}
	else if(number<10000)
	{
		strNum = "0" + convert_ItoS(number);
	}
	else
	{
		strNum = convert_ItoS(number);
	}

	return strNum;
}

//visszaadja egy vektor átlagát
double getVectorsAVG (vector<double> V)
{
	int N = V.size();
	double avgV = 0.0;
	for (int i=0; i<N; i++)
		avgV += V[i];

	avgV /= N;
	return avgV;
}

//******************************************************************//
//Matrix szorzasok
//******************************************************************//
vector<vector<double> >  mulSqMatrixes(vector<vector<double> >  A, vector<vector<double> >  B)
{
	int N = A.size();
	vector<double> vecZeros(N,0.0);
	vector<vector<double> >  C(N,vecZeros);

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			for (int k=0; k<N; k++)
				C[i][j] += A[i][k] * B[k][j];			
			
	return C;
}

vector<vector<double> >  mulSqMatrixes_inInt2(vector<vector<int> >  A, vector<vector<int> >  B)
{
	int N = A.size();
	vector<double> vecZeros(N,0.0);
	vector<vector<double> >  C(N,vecZeros);

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			for (int k=0; k<N; k++)
				C[i][j] += A[i][k] * B[k][j];			
			
	return C;
}

vector<vector<double> >  mulSqMatrixes_inInt1(vector<vector<int> >  A, vector<vector<double> >  B)
{
	int N = A.size();
	vector<double> vecZeros(N,0.0);
	vector<vector<double> >  C(N,vecZeros);

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			for (int k=0; k<N; k++)
				C[i][j] += A[i][k] * B[k][j];			
			
	return C;
}

//Visszaadja egy vector elemeinek a szorasat
double get_SqrtVar_ofVec(vector<double> V)
{
	double varV = 0.0;
	int N = V.size();
	double avgV = 0.0;
	double tmp;

	for (int i=0; i<N; i++)
		avgV += V[i];

	avgV /= double(N);

	for (int i=0; i<N; i++)
	{
		tmp = V[i]-avgV;
		varV += tmp * tmp;
	}

	varV /= N;

	return sqrt(varV);
}

//Visszaadja egy matrix (összes!) elemeinek a szorasat
double get_SqrtVar_ofMat(vector<vector<double> > M)
{
	double varM = 0.0;
	int row = M.size(); int col = M[0].size();
	int N = row*col;
	double avgM = 0.0;
	double tmp;

	for (int i=0; i<row; i++)
		for (int j=0; j<col; j++)
			avgM += M[i][j];

	avgM /= double(N);

	for (int i=0; i<row; i++)
	{
		for (int j=0; j<col; j++)
		{
			tmp = M[i][j]-avgM;
			varM += tmp * tmp;
		}
	}

	varM /= N;

	return sqrt(varM);
}

/////----- Program-specifikus fuggvenyek -----/////
//J matrix felbontasa (J_pos[i][j]: 1, ha van kapcsolat; 0, ha nincs | J_val[i][j]: 1, ha +1-es a kapcsolat; 0, ha -1-es)
void giveVal_and_Pos_matrixes(vector<vector<int> >  J, vector<vector<bool> > & J_pos, vector<vector<bool> > & J_val)
{
	int N=J.size();

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (J[i][j] == 0)
			{
				J_pos[i][j] = false;
				J_val[i][j] = false;		//ez vegulis mind1...
			}
			else if (J[i][j] == 1)
			{
				J_pos[i][j] = true;
				J_val[i][j] = true;
			}
			else if (J[i][j] == -1)
			{
				J_pos[i][j] = true;
				J_val[i][j] = false;
			}
		}
	}
}


//Visszaadja egy vektorba egy (double) Matrix felso haromszoget (atlo nelkul) - CDF abrakhoz szukseges
vector<double> getUpperTri_fromSqMatrix_doub (vector<vector<double> >  M)
{
	int size = M.size();
	vector<double> V(int((size*(size-1))/2));
	int count = 0;

	for (int i=0; i<size; i++)
	{
		for (int j=i+1; j<size; j++)
		{
			V[count] = M[i][j];
			count++;
		}
	}

	return V;
}

//Visszaadja egy matrixnak azon sorok szamait, ahol a szum(abs(xi)) max ill min
void giveExtrLines_of_matrix (vector<vector<double> >  M, int& lineMin, int& lineMax, double& minVal, double& maxVal)
{
	int N = M.size();
	double summ = 0.0;
	lineMin = lineMax = 0;

	for (int j=1; j<N; j++)		//j=1tol, h i!=j teljesuljon
		summ += abs(M[0][j]);

	minVal = summ; lineMax = 0;	//elso sor lesz a min
	maxVal = summ; lineMin = 0;	//elso sor lesz a max
	
	for (int i=1; i<N; i++)
	{
		summ = 0.0;
		for (int j=0; j<N; j++)
			if (i!=j) summ += abs(M[i][j]);
		
		if (summ > maxVal)
		{
			maxVal = summ;
			lineMax = i;
		}

		if (summ < minVal)
		{
			minVal = summ;
			lineMin = i;
		}
	}
}

//visszaadja egy CG J matrix +/- frusztralt haromszogeinek szamat stringben....
string fr_arany_J(vector<vector<int> > J)
{
	string arany = "";
	int N = J.size();
	int posFr = 0; int negFr = 0;
	int tmp = 0;

	for(int i=0; i<N; i++)
	{
		for(int j=i+1; j<N; j++)
		{
			for(int k=j+1; k<N; k++)
			{
				tmp = J[i][j] * J[i][k] * J[j][k];
				if (tmp == 1)	posFr++;
				if (tmp == -1)	negFr++;
			}
		}
	}

	arany = " posFR="+convert_ItoS(posFr)+",negFr="+convert_ItoS(negFr)+",negArany="+convert_DtoS(double(negFr)/double(posFr+negFr));

	return arany;
}

//a matrix 2. oszlopa adatot tartalmaz, mig ez a fuggveny beallitja hozza az x-tengelyt az elso [0] oszlopba 0-tol kezdve freq lepesekkel
void setXaxis_forVecInMat_D(vector<vector<double> >& M, double minVal, double freq)
{
	int size = M[0].size();
	for (int i=0; i<size; i++)
		M[0][i] = minVal + double(i)*freq; 
}

void setXaxis_forVecInMat_I(vector<vector<int> >& M, int minVal, int freq)
{
	int size = M[0].size();
	for (int i=0; i<size; i++)
		M[0][i] = minVal + i*freq; 
}

void givePosNegNum(vector<vector<int> > J, int& numPos, int& numNeg)
{
	int N = J.size();
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (J[i][j] == 1) numPos++;
			if (J[i][j] == -1) numNeg++;
		}
	}
}
//#endif
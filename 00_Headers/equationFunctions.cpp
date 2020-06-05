//#ifndef SAVEFILES_H
//#define SAVEFILES_H

//*************************************************//
//MMC-cWPS fuggetlen, egyenleteket leiro fuggvenyek//
//*************************************************//

using namespace std;

//[1] keplet: Adott J-re es s-re az energia szamitasa
double getEnergy (vector<vector<bool> >  J_pos, vector<vector<bool> >  J_val, vector<bool> s, vector<double> h)
{
	int N = J_pos.size();
	double energy = 0.0;

	for (int i=0; i<N; i++)
	{
		for (int j=i+1; j<N; j++)
		{
			if (J_pos[i][j] == true)		
			{
				if ((J_val[i][j] ^ (s[i] ^ s[j])) == true)		//ekvivalens ezzel: energy -= J[i][j] * s[i] * s[j];
					energy--;	
				else
					energy++;	
			}
		}
	}

	energy /= sqrt(double(N));		//normalas
	//energy /= double(N);

	for (int i=0; i<N; i++)
	{
		if (s[i] == true)										//ekvivalens ezzel: energy -= h * s[i];
			energy -= h[i];
		else
			energy += h[i];
	}

	return energy;
}

double getEnergy_justPairing (vector<vector<bool> >  J_pos, vector<vector<bool> >  J_val, vector<bool> s)
{
	int N = J_pos.size();
	double energy = 0.0;

	for (int i=0; i<N; i++)
	{
		for (int j=i+1; j<N; j++)
		{
			if (J_pos[i][j] == true)		
			{
				if ((J_val[i][j] ^ (s[i] ^ s[j])) == true)		//ekvivalens ezzel: energy -= J[i][j] * s[i] * s[j];
					energy--;	
				else
					energy++;	
			}
		}
	}

	energy /= sqrt(double(N));		//normalas

	return energy;
}

//[4]. keplet: A magnesseg kiszamitasa
int getMagnetism(vector<bool> spin)
{
	int m = 0;
	int N = spin.size();
	
	for (int i=0; i<N; i++)
		m += (spin[i] == true) ? +1 : -1;
	
	//m /= double(N);

	return m;
}

//m_pddf_sum 2.oszlopaba szamolja az eloszlast m alapjan
void getPDDF_forM(vector<vector<int> > m, vector<vector<int> > & m_pddf_sum)
{
	int dataNum = m[0].size();
	int distSize = m_pddf_sum[0].size();
	
	for (int i=0; i<dataNum; i++)
	{
		for (int j=0; j<distSize; j++)
		{
			if (m[1][i] == m_pddf_sum[0][j])
			{
				m_pddf_sum[1][j]++;
				break;
			}
		}
	}
}

//[9] keplet: Kapcsolt-Korrelacios Matrix kiszamitasa (Connected Correlation Matrix)
vector<vector<double> >  getConCorrMatrix(vector<double> s_avg, vector<vector<double> >  s2_avg)
{
	int N = s_avg.size();
	vector<double> vecZeros(N,0.0);
	vector<vector<double> >  C(N,vecZeros);
	
	for (int i=0; i<N; i++)
		for (int j=i+1; j<N; j++)
			C[i][j] = C[j][i] = s2_avg[i][j] - (s_avg[i]*s_avg[j]);		//onmagaval nem korrelalt (0) 

	return C;
}

//[10] keplet: Korrelacios Matrix kiszamitasa
vector<vector<double> >  getCorrMatrix(vector<vector<double> >  s2_avg)
{
	int N = s2_avg.size();
	vector<double> vecZeros(N,0.0);
	vector<vector<double> >  C(N,vecZeros);
	
	for (int i=0; i<N; i++)
		for (int j=i+1; j<N; j++)
			C[i][j] = C[j][i] = s2_avg[i][j];		//onmagaval nem korrelalt (0) 

	return C;
}

//visszaadja <si><sj> matrixot
vector<vector<double> > getSiSj(vector<double> s_avg)
{
	int N = s_avg.size();
	vector<double> vecZeros(N,0.0);
	vector<vector<double> >  sisj(N,vecZeros);

	for (int i=0; i<N; i++)
		for (int j=i; j<N; j++)
			sisj[i][j] = sisj[j][i] = s_avg[i] * s_avg[j];

	return sisj;
}

//[13] keplet: int tipusu J matrix frusztraltsagi fokat visszaado fuggvenyy
double getJfrustrate_int(vector<vector<int> >  J)
{
	int N = J.size();
	double fi_J = 0.0;
	
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			for (int k=0; k<N; k++)
				fi_J += J[i][j] * J[j][k] * J[i][k];

	fi_J /= N*(N-1)*(N-2);
	
	return fi_J;
}

//[14] keplet: Null (ill legkisebb H) teru es adott k teru kovariancia matrixok "tavolsaga"
double get_dC0Ck(vector<vector<double> > C0, vector<vector<double> > Ck)
{
	double d = 0.0;
	int N = C0.size();

	for (int i=0; i<N; i++)
		for (int j=i+1; j<N; j++)
			d += abs(C0[i][j] - Ck[i][j]);

	d /= double(N)*double(N-1);

	return d;
}

//[15] keplet: Adott k teres kovariancia es J matrix "tavolsaga"
double get_dCkJ(vector<vector<double> > Ck, vector<vector<int> > J)
{
	double d = 0.0;
	int N = Ck.size();

	////signum
	//for (int i=0; i<N; i++)
	//	for (int j=i+1; j<N; j++)
	//		Ck[i][j] = (Ck[i][j] > 0) ? 1.0 : ( (Ck[i][j] < 0) ? -1.0 : 0.0);

	for (int i=0; i<N; i++)
		for (int j=i+1; j<N; j++)
			d += Ck[i][j] * double(J[i][j]);

	d /= double(N)*double(N-1);

	return d;
}



//#endif
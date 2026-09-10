#include <ilcplex/ilocplex.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <iostream>
ILOSTLBEGIN      

using namespace std;

typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix; 
typedef IloArray<IloNumArray>    NumMatrix;
typedef IloArray<IloIntVarArray> FloatMatrix;

const int Mnum = 13;   //M 
const int Nnum = 14;   // N 
const double eps = 0.30;
const int service_hours = 6;
const int P = 1000000; 
const int cycle = 20;  const double Acycle = 20.0;  
const int realization_num = 1000; const double Arealization_num = 1000.0;  

const double PI = 3.141593;


inline double compare(double x, double y)
{
	double min;
	min = (x < y) ? x : y;
	return min;
}
double findavg(double* a, int length) {
	double sum = 0;
	for (int i = 0; i < length; i++)
		sum += a[i];
	return sum / length;
}

double findmax(double* a, int length) {
	int i; double maxvalue = 0.0;
	for (i = 0; i < length; i++) {
		if (i == 0) {
			maxvalue = a[i];
		}
		else {
			if (maxvalue < a[i]) {
				maxvalue = a[i];
			}
		}
	}
	return maxvalue;
}

double findmin(double* a, int length) {
	int i; double minvalue = 0.0;
	for (i = 0; i < length; i++) {
		if (i == 0) {
			minvalue = a[i];
		}
		else {
			if (minvalue > a[i]) {
				minvalue = a[i];
			}
		}
	}
	return minvalue;
}

double random_uniform_distribution(double min, double max) {
	return min + (max - min) * rand() / (RAND_MAX + 1.0);   
}


double random_normal_distribution_BoxMuller(double mean, double stdv) {
	double u_1;
	double u_2;
	double R;
	double theta;
	u_1 = random_uniform_distribution(0, 1);
	u_2 = random_uniform_distribution(0.01, 1);
	R = sqrt((-2) * log(u_2));
	theta = 2 * PI * u_1;
	return mean + (R * sin(theta)) * stdv;
}

double random_triangle_distribution(double low, double miu, double upper) {

	double u = (double)rand() / RAND_MAX;
	double tri;
	if (u <= (miu - low) / (upper - low))

		tri = low + sqrt(u * (upper - low) * (miu - low));
	else
		tri = upper - sqrt((1 - u) * (upper - low) * (upper - miu));
	return tri;
}

IloEnv env;
//******************************Multi *************************
void Multi_solution(double solu[10], double f[Mnum], double h[Mnum], int c[Mnum][Nnum], double Q[Mnum], double D_mean[Nnum], double DL[Nnum], double DU[Nnum], double sumDL, double TI)
{
	IloNumArray YYY(env, Mnum), II(env, Mnum), III(env, Mnum); 
	IloNumArray uu(env, Nnum), vv(env, Mnum);
	ofstream eout("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/eoutoutput.txt");
	
	eout << "Starting" << endl;
	try
	{
		IloEnv env;
		
		clock_t multi_start, multi_finish;
		double multi_duration = 0;
		multi_start = clock();
		clock_t multi_MPstart, multi_MPfinish;
		double multi_MPduration = 0.0;
		multi_MPstart = clock();

		IloModel model(env);                         
		IloCplex cplex(model);
		int i, j;
		double  I_final[Mnum];
		int cs = 0;         
		NumVarMatrix x(env, Mnum);                  
		for (i = 0; i < Mnum; i++)
		{
			x[i] = IloNumVarArray(env, Nnum, 0, 1, ILOINT);   
		}
		IloNumVarArray I(env, Mnum, 0, P, ILOFLOAT);  
		IloNumVarArray y(env, Mnum, 0.0, 1.0, ILOINT);     

		NumVarMatrix q(env, Mnum);                  
		for (i = 0; i < Mnum; i++)
		{
			q[i] = IloNumVarArray(env);
		}
		

		IloExpr obj1(env);
		IloExpr LC(env);
		for (int i = 0; i < Mnum; i++)
		{
			for (int j = 0; j < Nnum; j++)
			{
				LC += c[i][j] * x[i][j];        
			}
		}
		

		for (int i = 0; i < Mnum; i++)
		{
			obj1 += f[i] * y[i] + h[i] * I[i];
		}
		model.add(IloMinimize(env, LC + obj1));
		obj1.clear();
		LC.clear();
		LC.end();
		obj1.end();
		eout << "Objective function successful" << endl;
		
		for (int i = 0; i < Mnum; i++)
		{
			for (int j = 0; j < Nnum; j++)
			{
				model.add(x[i][j] <= y[i]);
				model.add(x[i][j] >= 0);
			}
		}
		for (int i = 0; i < Mnum; i++)
		{
			model.add(Q[i] * y[i] - I[i] >= 0); 
		}
		for (int i = 0; i < Mnum; i++)
		{
			model.add(I[i] >= 0);
		}
		eout << "Constraint3 successful" << endl;
		
		
		//////**************************** Single Set*****************
		
		for (int j = 0; j < Nnum; j++) 
		{
			
			for (int i = 0; i < Mnum; i++) {
				
				q[i].add(IloNumVar(env, 0.0, P, ILOFLOAT)); 
				
				model.add(P * x[i][j] - q[i][j] >= 0);
				model.add(q[i][j] <= I[i]);
			}

			
			IloExpr expr2(env);
			for (int i = 0; i < Mnum; i++) {
				expr2 += q[i][j];
				
			}

			double TTi = 0.0;
			double S_no_j = 0.0;
			for (int k = 0; k < Nnum; k++) {
				if (k != j) {
					S_no_j += DL[k];
				}
			}


			TTi = TI - S_no_j;
			double minbd = compare(TTi, DU[j]);

			model.add(minbd <= expr2);
			expr2.clear();
			expr2.end();
		}
		

		if (!cplex.solve())
		{
			env.error() << "Failed to optimize LP." << endl;
		}
		else
		{
			cplex.setOut(env.getNullStream());
			cplex.setWarning(env.getNullStream());
			cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
			
			cplex.setParam(IloCplex::EpGap, 0);	//!!!relative gap
			cplex.setParam(IloCplex::EpAGap, 0);
			
			cplex.extract(model);               
			cplex.exportModel("transport.lp");
			////////////////////////////////////////////////////////////
			cplex.solve();
			if (cplex.getStatus() == IloAlgorithm::Infeasible)
				cout << "No Solution" << endl;
			eout << "Solution status: " << cplex.getStatus() << endl;
			
		}
		
		IloNum obj = 0;
		int XXXX[Mnum][Nnum];
		double  Q[Mnum][Nnum];
		double Real_Sum_I;
		Real_Sum_I = 0;

		for (i = 0; i < Mnum; i++)
		{
			for (j = 0; j < Nnum; j++)
			{
				
				if (cplex.getValue(x[i][j]) > 0.5) {
					XXXX[i][j] = 1;
				}
				else {
					XXXX[i][j] = 0;
				}
				
				eout << "x[" << i << "][" << j << "]= " << XXXX[i][j] << "   ";
				solu[1] += XXXX[i][j];
			}
			eout << endl;
		}
		eout << "----------------------------" << endl;
		
		for (int i = 0; i < Mnum; i++)
		{
			II[i] = cplex.getValue(I[i]);


			if (cplex.getValue(y[i]) > 0.5) {
				YYY[i] = 1;
			}
			else {
				YYY[i] = 0;
			}

			eout << "y[" << i << "]=" << YYY[i] << "          ";
			
			eout << "I[" << i << "]=" << II[i] << "          ";
			
			solu[2] += YYY[i];
			Real_Sum_I += II[i];
		}
		eout << "----------------------------" << endl;
		eout << "Total Arc" << solu[1] << endl;
		eout << "Real_Sum_I=" << Real_Sum_I << endl;
		eout << "Total open facility" << solu[2] << endl;

		////////////////////////////////////////////////////////////////////////////////////
		eout << "----------------" << endl;
		eout << "Best objective value：";
		obj = cplex.getObjValue();                   
		eout << obj << endl;

		multi_MPfinish = clock();
		multi_MPduration += ((double)multi_MPfinish - multi_MPstart) / ((double)CLK_TCK);
		/////////////////////////////////////////////////////////////////////////////**********************
		clock_t multi_SPstart, multi_SPfinish;
		double multi_SPduration = 0.0;

		double cut = 0;
		do
		{
			multi_SPstart = clock();
			
			eout << "--------------------------------------------SP starting-----------------------------------------------" << endl;
			eout << "---------------------------------------------SP starting----------------------------------------------" << endl;
			///////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			IloModel submodel(env);  
			IloNumVarArray v(env, Mnum, 0, 1, ILOINT);
			IloNumVarArray u(env, Nnum, 0, 1, ILOINT);
			IloNumVar z(env, 0.0, IloInfinity, ILOFLOAT);

			submodel.add(v);
			submodel.add(u);
			submodel.add(z);
														
			IloExpr exp0(env);
										
			for (int i = 0; i < Mnum; i++)
			{
				exp0 += II[i] * v[i];
			}

			IloExpr exp1(env);  
			for (int j = 0; j < Nnum; j++)
			{
				exp1 += DU[j] * u[j];
			}
			eout << "-----------submodel obejctive function-----------------" << endl;
			submodel.add(IloMinimize(env, exp0 - exp1 + z));
			eout << "obejctive function successful" << endl;
			eout << "Starting input constraints..." << endl;
			exp0.clear(); exp0.end();
			

			eout << "Two：" << endl;
			
			IloExpr exp2(env);
			for (int j = 0; j < Nnum; j++)
			{
				exp2 += DL[j] * (1 - u[j]);
			}
			submodel.add(z - exp1 + TI - exp2 >= 0);
			exp1.clear(); exp1.end();
			exp2.clear(); exp2.end();
			eout << "Constraint successful" << endl;
			eout << endl;
																
			////////////////////////////////////////////////////////////////////////////////////
			eout << "Most have" << Mnum * Nnum << "constraints：" << endl;
			for (int j = 0; j < Nnum; j++)
			{
				for (int i = 0; i < Mnum; i++)
				{
					if (XXXX[i][j] > 0.5)
					{
						submodel.add(v[i] >= u[j]);
						
					}
					else continue;
				}
			}
			eout << "Constraint successful" << endl;
			eout << "Starting..." << endl << endl;
			////////////////////////////////////////////////////////////////////////////////////
			double zout = 0.0;

			eout << "-----------submodel result：-----------------" << endl;
			
			IloCplex subsolver(submodel);
			
			subsolver.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
			
			subsolver.setParam(IloCplex::EpGap, 0);	//!!!relative gap
			subsolver.setParam(IloCplex::EpAGap, 0);

			cplex.setOut(env.getNullStream());
			cplex.setWarning(env.getNullStream());
			
			subsolver.solve();

			if (cplex.getStatus() == IloAlgorithm::Infeasible)
				eout << "No Solution" << endl;
			eout << "Solution status: " << cplex.getStatus() << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////

			double objyan = 0.0; 
			zout = subsolver.getValue(z);
			eout << "z=" << zout << endl;
			
			vector<int> uu(Nnum);
			vector<int> vv(Mnum);

			for (int i = 0; i < Mnum; i++)
			{
				if (subsolver.getValue(v[i]) > 0.5) {
					vv[i] = 1;
				}
				else {
					vv[i] = 0;
				}
				
				eout << "v[" << i << "]=" << vv[i] << " ";
				cout << "v[" << i << "]=" << vv[i] << " ";
			}
			eout << "----------------------------" << endl;

			
			for (int j = 0; j < Nnum; j++)
			{
				if (subsolver.getValue(u[j]) > 0.5) {
					uu[j] = 1;
				}
				else {
					uu[j] = 0;
				}

				eout << "u[" << j << "]=" << uu[j] << " ";
			}
			eout << endl << "---------------------------------------" << endl;
			eout << "Best objective value：";
			objyan = subsolver.getObjValue();
			eout << " objyan=" << objyan << "   if bigger than 0" << endl;
			
			zout = subsolver.getValue(z);
			eout << "z=" << zout << endl;
			subsolver.clearModel();
			subsolver.clear();
			subsolver.end();
			///////////////////////////////////////////////////////////////////////////////
			eout << "if" << endl;
			u.clear(); v.clear();
			u.end(); v.end(); z.end();//uu.end();vv.end();
			submodel.end(); //
			

			multi_SPfinish = clock();
			multi_SPduration += ((double)multi_SPfinish - multi_SPstart) / ((double)CLK_TCK);
			if (objyan < -0.10)
			{

				multi_MPstart = clock();
				cut++;
				eout << "--------------Ite continue--------------" << endl;
				
				///////////////////////////////////////////////////////////////////////////////
				IloNum ccx = 0, ccd = 0; double ccy = 0.0, mincc = 0.0;
				for (int j = 0; j < Nnum; j++)
				{
					ccx += DU[j] * uu[j];
					ccd += DL[j] * (1 - uu[j]);
				}
				eout << "ccx=" << ccx << endl;                     
				eout << "ccd=" << ccd << endl;
				ccy = TI - ccd;
				mincc = compare(ccx, ccy);
				eout << "mincc=" << mincc << endl;
				///////////////////////////////////////////////////////////////////////////////

				for (int i = 0; i < Mnum; i++) {
					
					q[i].add(IloNumVar(env, 0.0, P, ILOFLOAT)); 
						
					model.add(q[i][q[i].getSize() - 1] <= I[i]);
				}
				
				IloExpr expr3(env);
				for (int i = 0; i < Mnum; i++)
				{
					expr3 += q[i][q[i].getSize() - 1];
				}

				model.add(expr3 - mincc >= 0.0);               
				expr3.clear();
				expr3.end();

				
				for (int i = 0; i < Mnum; i++)
				{
					IloExpr X_is(env);

					for (j = 0; j < Nnum; j++)
					{
						X_is += P * x[i][j] * uu[j];
					}
					model.add(q[i][q[i].getSize() - 1] <= X_is);
						
					X_is.clear();
					X_is.end();
				}
				
				
				
				cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
				
				cplex.setOut(env.getNullStream());
				cplex.setWarning(env.getNullStream());
				cplex.solve();

				if (cplex.getStatus() == IloAlgorithm::Infeasible)
					eout << "No Solution" << endl;

				eout << "Solution status: " << cplex.getStatus() << endl;
				////////////////////////////////////////////////////////////////////////////////////
				eout << "-----------------Ite" << cut << "solution：-------------------" << endl;
				
				Real_Sum_I = 0.0;
				solu[1] = 0.0; solu[2] = 0.0;
				
				for (i = 0; i < Mnum; i++)
				{
					for (j = 0; j < Nnum; j++)
					{
						if (cplex.getValue(x[i][j]) > 0.5) {
							XXXX[i][j] = 1;
						}
						else {
							XXXX[i][j] = 0;
						}
						
						eout << "x[" << i << "][" << j << "]= " << XXXX[i][j] << "   ";
						solu[1] += XXXX[i][j];
						
					}
					eout << endl;
				}
				eout << "----------------------------" << endl;
				
				for (int i = 0; i < Mnum; i++)
				{
					II[i] = cplex.getValue(I[i]);
					YYY[i] = cplex.getValue(y[i]);
					
					eout << "y[" << i << "]=" << YYY[i] << "          ";
					
					eout << "I[" << i << "]=" << II[i] << " " << endl;
					solu[2] += YYY[i];
					I_final[i] = II[i]; 
					Real_Sum_I += II[i];
				}
				eout << "Total Arc" << solu[1] << endl;
				eout << "Real_Sum_I=" << Real_Sum_I << endl;
				eout << "Total open facility" << solu[2] << endl;
				////////////////////////////////////////////////////////////////////////////////////
				obj = cplex.getObjValue();
				eout << "New objective value" << obj << endl;
				vv.clear(); uu.clear();
				vv.end(); uu.end();
				multi_MPfinish = clock();
				multi_MPduration += ((double)multi_MPfinish - multi_MPstart) / ((double)CLK_TCK);
			}

			else break;
		} while (1);

		multi_finish = clock();
		multi_duration = ((double)multi_finish - multi_start) / ((double)CLK_TCK);
		

		ofstream xout("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/XIJ.txt");
		for (i = 0; i < Mnum; i++)
		{
			for (j = 0; j < Nnum; j++)
			{

				xout << "   " << XXXX[i][j] << "   ";
				
			}
			xout << endl;
		}
		

		ofstream yout("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/YI.txt");
		for (i = 0; i < Mnum; i++)
		{
			yout << "   " << YYY[i] << "   ";
		}
		ofstream Iout("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/II.txt");
		for (i = 0; i < Mnum; i++)
		{
			Iout << "   " << II[i] << "   ";
		}
		YYY.end();  I.end(); y.end(); II.end();    
		cplex.clearModel();
		cplex.end();
		model.end();
		eout << "------------------------End------------------------" << endl;
		eout.end;
		eout << "Total Arc" << solu[1] << endl;
		eout << "Real_Sum_I=" << Real_Sum_I << endl;
		eout << "Total open facility" << solu[2] << endl;
		////////////////////////////////////////////////////////////////////////////////////

		eout << "New objective value" << obj << endl;

		/////////////////////////////////////////////////////////////////////////////////////
		eout << "_______Real_Sum_I   =   " << Real_Sum_I << endl;

		env.end();

		solu[0] = multi_duration;
		solu[3] = Real_Sum_I;
		solu[5] = obj;
		eout << " CPU time: " << solu[0] << endl;
		solu[8] = multi_MPduration;
		solu[9] = multi_SPduration;
		solu[10] = cut;

		eout << "MP Time: " << multi_MPduration << endl;
		eout << "SP Time: " << multi_SPduration << endl;
		eout << " cut number: " << cut << endl << endl;
		
	}
	catch (IloAlgorithm::CannotExtractException& e)
	{
		std::cerr <<
			"CannoExtractException: " << e << std::endl; IloExtractableArray failed = e.getExtractables();

		for (IloInt i = 0; i < failed.getSize(); ++i) std::cerr <<
			"\t" << failed[i] << std::endl;
		// Handle exception ... 
	}
	catch (IloException& ex)
	{
		cerr << "Error: " << ex << endl;
		cerr << ex.getMessage();
		ex.end();
	}
	catch (...)
	{
		cerr << "Error" << endl;
	}
}

//******************************Single *************************

void Single_solution(double single_solu[10], double f[Mnum], double h[Mnum], int c[Mnum][Nnum], double Q[Mnum], double D_mean[Nnum], double DL[Nnum], double DU[Nnum], double sumDL, double TI) //函数声明，就是具体的求解方法  
{
	IloNumArray YYY(env, Mnum), II(env, Mnum), III(env, Mnum); 
	IloNumArray uu(env, Nnum), vv(env, Mnum);
	ofstream gout("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/goutoutput.txt");
	
	gout << "Starting" << endl;
	try
	{
		int LimtTime = 0;
		LimtTime = 60;

		IloEnv env;
		clock_t single_start, single_finish, single_terminate;
		double solve = 1.0;
		double single_duration = 0.0;
		single_start = clock();
		clock_t single_MPstart, single_MPfinish;
		double single_MPduration = 0.0;
		single_MPstart = clock();

		IloModel model(env);                         
		IloCplex cplex(model);
		int i, j;
		double  I_final[Mnum];
		int cs = 0;         
		NumVarMatrix x(env, Mnum);                  
		for (i = 0; i < Mnum; i++)
		{
			x[i] = IloNumVarArray(env, Nnum, 0, 1, ILOBOOL);   
		}
		IloNumVarArray I(env, Mnum, 0, P, ILOFLOAT);   //Ii 
		IloNumVarArray y(env, Mnum, 0.0, 1.0, ILOBOOL);      // YYYi 

		NumVarMatrix q(env, Mnum);                    //XXXij 
		for (i = 0; i < Mnum; i++)
		{
			q[i] = IloNumVarArray(env);
		}

		IloExpr obj1(env);
		IloExpr LC(env);
		for (int i = 0; i < Mnum; i++)
		{
			for (int j = 0; j < Nnum; j++)
			{
				LC += c[i][j] * x[i][j];         
			}
		}
		

		for (int i = 0; i < Mnum; i++)
		{
			obj1 += f[i] * y[i] + h[i] * I[i];
		}
		model.add(IloMinimize(env, LC + obj1));
		obj1.clear();
		LC.clear();
		LC.end();
		obj1.end();
		gout << "Objective function successful" << endl;
		

		for (int i = 0; i < Mnum; i++)
		{
			for (int j = 0; j < Nnum; j++)
			{
				model.add(x[i][j] <= y[i]);
				model.add(x[i][j] >= 0);
			}
		}
		for (int i = 0; i < Mnum; i++)
		{
			model.add(Q[i] * y[i] - I[i] >= 0);  
		}
		
		for (int i = 0; i < Mnum; i++)
		{
			model.add(I[i] >= 0);
		}
		gout << "Constraint 3 successful" << endl;
		//******************************Single constraint *************************

		for (j = 0; j < Nnum; j++)
		{
			IloExpr single(env);
			single.clear();
			for (int i = 0; i < Mnum; i++)

			{
				single += x[i][j];  
			}

			model.add(single == 1);
			
			single.end();
		}
		

		//*********************************************

		for (int j = 0; j < Nnum; j++)  
		{

			for (int i = 0; i < Mnum; i++) {
				
				q[i].add(IloNumVar(env, 0.0, P, ILOFLOAT)); 
				model.add(P * x[i][j] - q[i][j] >= 0);
				model.add(q[i][j] <= I[i]);
			}


			IloExpr expr2(env);
			for (int i = 0; i < Mnum; i++) {
				expr2 += q[i][j];

			}

			double TTi = 0.0;
			double S_no_j = 0.0;
			for (int k = 0; k < Nnum; k++) {
				if (k != j) {
					S_no_j += DL[k];
				}
			}


			TTi = TI - S_no_j;
			double minbd = compare(TTi, DU[j]);

			model.add(minbd <= expr2);
			expr2.clear();
			expr2.end();
		}

		
		gout << "............11....." << endl;


		if (!cplex.solve())
		{
			env.error() << "Failed to optimize LP." << endl;
		}
		else
		{
			cplex.setOut(env.getNullStream());
			cplex.setWarning(env.getNullStream());
			cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
			
			cplex.setParam(IloCplex::EpGap, 0);	//!!!relative gap
			cplex.setParam(IloCplex::EpAGap, 0);
			
			cplex.extract(model);            
			cplex.exportModel("transport.lp");
			////////////////////////////////////////////////////////////
			cplex.solve();
			if (cplex.getStatus() == IloAlgorithm::Infeasible)
				cout << "No Solution" << endl;
			gout << "Solution status: " << cplex.getStatus() << endl;
			
		}
		
		IloNum obj = 0;

		int  XXXX[Mnum][Nnum];
		double Real_Sum_I;
		Real_Sum_I = 0;

		for (i = 0; i < Mnum; i++)
		{
			for (j = 0; j < Nnum; j++)
			{
				if (cplex.getValue(x[i][j]) > 0.5) {
					XXXX[i][j] = 1;
				}
				else {
					XXXX[i][j] = 0;
				}
				
				gout << "x[" << i << "][" << j << "]= " << XXXX[i][j] << "   ";
				single_solu[1] += XXXX[i][j];
			}
			gout << endl;
		}
		gout << "----------------------------" << endl;
		
		for (int i = 0; i < Mnum; i++)
		{
			II[i] = cplex.getValue(I[i]);

			if (cplex.getValue(y[i]) > 0.5) {
				YYY[i] = 1;
			}
			else {
				YYY[i] = 0;
			}
			
			gout << "y[" << i << "]=" << YYY[i] << "          ";
			
			gout << "I[" << i << "]=" << II[i] << "          ";
			
			single_solu[2] += YYY[i];
			Real_Sum_I += II[i];
		}
		gout << "----------------------------" << endl;
		gout << "Total Arc" << single_solu[1] << endl;
		gout << "Real_Sum_I=" << Real_Sum_I << endl;
		gout << "Total open facility" << single_solu[2] << endl;

		////////////////////////////////////////////////////////////////////////////////////
		gout << "----------------" << endl;
		gout << "Best objective value：";
		obj = cplex.getObjValue();                   
		gout << obj << endl;

		single_MPfinish = clock();
		single_MPduration += ((double)single_MPfinish - single_MPstart) / ((double)CLK_TCK);
		/////////////////////////////////////////////////////////////////////////////**********************
		clock_t single_SPstart, single_SPfinish;
		double single_SPduration = 0.0;


		double cut = 0;
		do
		{
			single_terminate = clock();
			if (((double)single_terminate - single_start) / ((double)CLK_TCK) > LimtTime)
			{
				break;
			}                                                              
			single_SPstart = clock();
			
			gout << "--------------------------------------------SP starting-----------------------------------------------" << endl;
			gout << "---------------------------------------------SP starting----------------------------------------------" << endl;
			///////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			IloModel submodel(env);  
			IloNumVarArray v(env, Mnum, 0, 1, ILOINT);
			IloNumVarArray u(env, Nnum, 0, 1, ILOINT);
			IloNumVar z(env, 0.0, IloInfinity, ILOFLOAT);
			
			submodel.add(v);
			submodel.add(u);
			submodel.add(z);

			IloExpr exp0(env); 

			for (int i = 0; i < Mnum; i++)
			{
				exp0 += II[i] * v[i];
			}
			

			IloExpr exp1(env);  
			for (int j = 0; j < Nnum; j++)
			{
				exp1 += DU[j] * u[j];
			}
			gout << "-----------submodel objective function-----------------" << endl;
			submodel.add(IloMinimize(env, exp0 - exp1 + z));
			gout << "objective function successful" << endl;
			gout << "Constraint..." << endl;
			exp0.clear(); exp0.end();
			
			gout << "Two：" << endl;
			
			IloExpr exp2(env);
			for (int j = 0; j < Nnum; j++)
			{
				exp2 += DL[j] * (1 - u[j]);
			}
			submodel.add(z - exp1 + TI - exp2 >= 0);
			exp1.clear(); exp1.end();
			exp2.clear(); exp2.end();
			gout << "Constraint successful" << endl;
			gout << endl;
			
			////////////////////////////////////////////////////////////////////////////////////
			gout << "Most have" << Mnum * Nnum << "constraints：" << endl;
			for (int j = 0; j < Nnum; j++)
			{
				for (int i = 0; i < Mnum; i++)
				{
					if (XXXX[i][j] > 0.5)
					{
						submodel.add(v[i] >= u[j]);
					}
					else continue;
				}
			}

			gout << "Constraint successful" << endl;
			gout << "Starting..." << endl << endl;
			////////////////////////////////////////////////////////////////////////////////////
			double zout = 0.0;

			gout << "-----------submodel result：-----------------" << endl;
			IloCplex subsolver(submodel);
			subsolver.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
			
			subsolver.setParam(IloCplex::EpGap, 0);	//!!!relative gap
			subsolver.setParam(IloCplex::EpAGap, 0);
			
			subsolver.setParam(IloCplex::TiLim, 3600); 
														////////////////////////////////////////////////////////////////////////////////////
			cplex.setOut(env.getNullStream());
			cplex.setWarning(env.getNullStream());
			subsolver.extract(submodel);
			subsolver.solve();

			if (subsolver.getStatus() == IloAlgorithm::Infeasible)
				gout << "No Solution" << endl;
			gout << "Solution status: " << subsolver.getStatus() << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////
			double objyan = 0.0; 
			zout = subsolver.getValue(z);
			gout << "z=" << zout << endl;
			
			vector<int> uu(Nnum);
			vector<int> vv(Mnum);

			for (int i = 0; i < Mnum; i++)
			{
				if (subsolver.getValue(v[i]) > 0.5) {
					vv[i] = 1;
				}
				else {
					vv[i] = 0;
				}
				
				gout << "v[" << i << "]=" << vv[i] << " ";
			}
			gout << "----------------------------" << endl;
			for (int j = 0; j < Nnum; j++)
			{
				if (subsolver.getValue(u[j]) > 0.5) {
					uu[j] = 1;
				}
				else {
					uu[j] = 0;
				}
				
				gout << "u[" << j << "]=" << uu[j] << " ";
			}
			gout << endl << "---------------------------------------" << endl;
			gout << "Best objective value：";
			objyan = subsolver.getObjValue();
			gout << " objyan=" << objyan << "   if bigger than 0" << endl;
			zout = subsolver.getValue(z);
			gout << "z=" << zout << endl;
			subsolver.clearModel();
			subsolver.clear();
			subsolver.end();
			///////////////////////////////////////////////////////////////////////////////
			gout << "if" << endl;
			u.clear(); v.clear();
			u.end(); v.end(); z.end();//uu.end();vv.end();
			submodel.end(); //
			single_SPfinish = clock();
			single_SPduration += ((double)single_SPfinish - single_SPstart) / ((double)CLK_TCK);

			

			if (objyan < -0.1) 
			{

				single_MPstart = clock();
				cut++;
				gout << "--------------Ite continue--------------" << endl;
				///////////////////////////////////////////////////////////////////////////////
				IloNum ccx = 0, ccd = 0; double ccy = 0.0, mincc = 0.0;
				for (int j = 0; j < Nnum; j++)
				{
					ccx += DU[j] * uu[j];
					ccd += DL[j] * (1 - uu[j]);
				}
				gout << "ccx=" << ccx << endl;                    
				gout << "ccd=" << ccd << endl;
				ccy = TI - ccd;
				gout << "ccy=" << ccy << endl;
				mincc = compare(ccx, ccy);
				gout << "mincc=" << mincc << endl;
				///////////////////////////////////////////////////////////////////////////////

				for (int i = 0; i < Mnum; i++) {
					q[i].add(IloNumVar(env, 0.0, P, ILOFLOAT)); 
					model.add(q[i][q[i].getSize() - 1] <= I[i]);
					
				}

				IloExpr expr3(env);
				for (int i = 0; i < Mnum; i++)
				{
					expr3 += q[i][q[i].getSize() - 1];
				}

				model.add(expr3 - mincc >= 0.0);                
				expr3.clear();
				expr3.end();


				for (int i = 0; i < Mnum; i++)
				{
					IloExpr X_is(env);

					for (j = 0; j < Nnum; j++)
					{
						if (uu[j] > 0.5) {
							X_is += Q[i] * x[i][j]; //X_is += P * x[i][j] * uu[j];
						}
						
					}
					model.add(q[i][q[i].getSize() - 1] <= X_is);
					
					X_is.clear();
					X_is.end();
				}

				cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
				
				cplex.setParam(IloCplex::TiLim, LimtTime);
				cplex.setOut(env.getNullStream());
				cplex.setWarning(env.getNullStream());
				cplex.extract(model);
				cplex.solve();


				////////////////////////////////////////////////////////////////////////////////////
				gout << "-----------------Ite" << cut << "solution：-------------------" << endl;
				
				Real_Sum_I = 0.0;
				single_solu[1] = 0.0; single_solu[2] = 0.0;
				

				for (i = 0; i < Mnum; i++)
				{
					int Size_q = 0;
					Size_q = q[i].getSize();
					gout << "q[" << i << "]： " << Size_q << "; value：" << cplex.getValue(q[i][Size_q - 1]) << "   ";
					gout << endl;
				}
				gout << "----------------------------" << endl;


				for (i = 0; i < Mnum; i++)
				{
					for (j = 0; j < Nnum; j++)
					{

						if (cplex.getValue(x[i][j]) > 0.5) {
							XXXX[i][j] = 1;
						}
						else {
							XXXX[i][j] = 0;
						}

						
						gout << "x[" << i << "][" << j << "]= " << XXXX[i][j] << "   ";
						single_solu[1] += XXXX[i][j];
						
					}
					gout << endl;
				}
				gout << "----------------------------" << endl;
				
				for (int i = 0; i < Mnum; i++)
				{
					II[i] = cplex.getValue(I[i]);
					YYY[i] = cplex.getValue(y[i]);
					
					gout << "y[" << i << "]=" << YYY[i] << "          ";
					
					gout << "I[" << i << "]=" << II[i] << " " << endl;
					single_solu[2] += YYY[i];
					I_final[i] = II[i]; 
					Real_Sum_I += II[i];
				}
				gout << "Total Arc" << single_solu[1] << endl;
				gout << "Real_Sum_I=" << Real_Sum_I << endl;
				gout << "Total open facility" << single_solu[2] << endl;
				////////////////////////////////////////////////////////////////////////////////////
				obj = cplex.getObjValue();

				gout << "===================================" << endl;


				gout << "New objective value" << obj << endl;
				vv.clear(); uu.clear();
				vv.end(); uu.end();
				single_MPfinish = clock();
				single_MPduration += ((double)single_MPfinish - single_MPstart) / ((double)CLK_TCK);
			}
			else break;
		} while (1);

		single_finish = clock();
		single_duration = ((double)single_finish - single_start) / ((double)CLK_TCK);
		////////////////////////////////////////////////////////////////////////////////////

		ofstream xout("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/single_XIJ.txt");
		for (i = 0; i < Mnum; i++)
		{
			for (j = 0; j < Nnum; j++)
			{

				xout << "   " << XXXX[i][j] << "   ";
				
			}
			xout << endl;
		}


		ofstream yout("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/single_YI.txt");
		for (i = 0; i < Mnum; i++)
		{
			yout << "   " << YYY[i] << "   ";
		}
		ofstream Iout("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/single_II.txt");
		for (i = 0; i < Mnum; i++)
		{
			Iout << "   " << II[i] << "   ";
		}
		YYY.end();  I.end(); y.end(); II.end();   
		cplex.clearModel();
		cplex.end();
		model.end();
		gout << "------------------------End------------------------" << endl;
		gout.end;
		gout << "Total Arc" << single_solu[1] << endl;
		gout << "Real_Sum_I=" << Real_Sum_I << endl;
		gout << "Total open facility" << single_solu[2] << endl;
		////////////////////////////////////////////////////////////////////////////////////

		gout << "New objective value" << obj << endl;

		/////////////////////////////////////////////////////////////////////////////////////
		gout << "Real_Sum_I   =   " << Real_Sum_I << endl;

		env.end();

		single_solu[0] = single_duration;
		single_solu[3] = Real_Sum_I;
		single_solu[5] = obj;
		gout << " CPU time: " << single_solu[0] << endl;
		single_solu[8] = single_MPduration;
		single_solu[9] = single_SPduration;
		single_solu[10] = cut;

		gout << "MP Time: " << single_MPduration << endl;
		gout << "SP Time: " << single_SPduration << endl;
		gout << " cut number: " << cut << endl << endl;
		
	}
	catch (IloAlgorithm::CannotExtractException& e)
	{
		std::cerr <<
			"CannoExtractException: " << e << std::endl; IloExtractableArray failed = e.getExtractables();

		for (IloInt i = 0; i < failed.getSize(); ++i) std::cerr <<
			"\t" << failed[i] << std::endl;
		// Handle exception ... 
	}
	catch (IloException& ex)
	{
		cerr << "Error: " << ex << endl;
		cerr << ex.getMessage();
		ex.end();
	}
	catch (...)
	{
		cerr << "Error" << endl;
	}
}


//FillRate&Chance

IloModel ScLP(env);	
IloCplex ScLPSolver(ScLP);
IloObjective ScLPObj = IloAdd(ScLP, IloMaximize(env));

IloRangeArray ScLPCon1(env);
IloRangeArray ScLPCon2(env);
IloRangeArray ScLPCon3(env);
NumVarMatrix tau(env, Mnum);
void FormulateScLP(IloModel ScLP, IloObjective ScLPObj, IloRangeArray ScLPCon1, IloRangeArray ScLPCon3, int x[Mnum][Nnum], double I[Mnum])
{
	//----- Objective function.-----//
	IloExpr Obj_temp(env);
	char tauName[20];
	for (IloInt i = 0; i < Mnum; i++)
	{
		tau[i] = IloNumVarArray(env, Nnum);
		for (IloInt j = 0; j < Nnum; j++)
		{
			sprintf_s(tauName, "tau(%d,%d)", i, j);
			tau[i][j] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT, tauName);
		}
	}
	for (IloInt i = 0; i < Mnum; i++)
	{
		for (IloInt j = 0; j < Nnum; j++)
		{
			Obj_temp += tau[i][j];
		}
	}
	ScLPObj.setExpr(Obj_temp);	//obj=max sum_ij tau
	Obj_temp.end();

	//----- Constrain No.1.-----// sum_j tau[i][j] <= I_i \forall i
	ScLP.remove(ScLPCon1);
	ScLPCon1.clear();
	for (IloInt i = 0; i < Mnum; i++)
	{
		IloExpr Con1Temp(env);
		for (IloInt j = 0; j < Nnum; j++)
		{
			Con1Temp += tau[i][j];
		}
		ScLPCon1.add(Con1Temp - I[i] <= 0);
		Con1Temp.end();
	}
	ScLP.add(ScLPCon1);

	//----- Constrain No.3.-----// tau[i][j] =0 \forall X[i][j]=0
	ScLP.remove(ScLPCon3);
	ScLPCon3.clear();
	for (IloInt i = 0; i < Mnum; i++)
	{
		for (IloInt j = 0; j < Nnum; j++)
		{
			if (abs(x[i][j]) < 0.1) //if X[i][j]=0
			{
				ScLPCon3.add(tau[i][j] <= 0);
			}
		}
	}
	ScLP.add(ScLPCon3);
}
void FillRateLoop(double fill[2], IloModel ScLP, IloRangeArray ScLPCon2, IloCplex& ScLPSolver, double D[realization_num][Nnum])
{
	//Button_JOC: 0->JOC model; 1->Convex model; 2->Chebyshev modell 3->SAA
	double AveFillRate = 0;
	double AveChance = 0;
	double TotalDemand = 0;
	double FillRate = 0;
	int Chance = 0;

	for (IloInt n = 0; n < realization_num; n++)
	{
		//cout<<endl<<"---------- Scenario: "<<n+1<<"/"<<nbScenarioFR<<" ----------"<<endl;

		//realize demand
		TotalDemand = 0;
		for (IloInt j = 0; j < Nnum; j++)
		{
			TotalDemand += D[n][j];
		}

		//----- Constrain No.2.-----// sum_i tau[i][j] <= D_j \forall j
		ScLP.remove(ScLPCon2);
		ScLPCon2.clear();	//clear ScLPCon2 in ScLP
		for (IloInt j = 0; j < Nnum; j++)
		{
			IloExpr Con2Temp(env);
			for (IloInt i = 0; i < Mnum; i++)
			{
				Con2Temp += tau[i][j];
			}
			ScLPCon2.add(Con2Temp - D[n][j] <= 0);
			Con2Temp.end();
		}
		ScLP.add(ScLPCon2);

		////----- Print the formulation: ScLP -----//
		//cout<<endl<<"*** the ScLP formulation"<<endl;
		//cout<<endl<<ScLP<<endl;

		//----- Solution -----//
		//cout<<endl<<"----- Solving ScLP -----"<<endl;
		ScLPSolver.solve();

		//cout<<endl<<"----- *** Solution of ScLP -----"<<endl;
		//cout<<"Status:"<<'\t'<<ScLPSolver.getStatus()<<endl;
		//cout<<"SubObj:"<<'\t'<<ScLPSolver.getObjValue()<<endl;
		//cout<<endl<<"Flow from i to j, tau[i][j]"<<endl;
		//PrintMatrix(ScLPSolver,tau,"S","D");

		//Calculate FillRate & Chance
		FillRate = (double)ScLPSolver.getObjValue() / TotalDemand;
		//	if (abs(FillRate - 1)<0.001)
		if (TotalDemand - (double)ScLPSolver.getObjValue() < 0.0001)
			Chance = 1;
		else
			Chance = 0;

		AveFillRate += FillRate;
		AveChance += Chance;

		//cout<<endl<<"*** TotalDemand = "<<TotalDemand;
		//cout<<endl<<"*** Fill rate = "<<FillRate <<endl;	//<<setprecision(5)<<fixed

	}	//END: scenario loop

	AveFillRate = (double)AveFillRate / Arealization_num;
	AveChance = (double)AveChance / Arealization_num;
	cout << "Average fill rate = " << AveFillRate << endl;
	cout << "Average chance = " << AveChance << endl;
	fill[0] = AveFillRate;
	fill[1] = AveChance;
}


int main(int argc, char** argv)
{
	IloInt i, j;
	ofstream fout("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/foutPoutput.txt");
	ofstream result("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/resultoutput.txt");

	fout << "-------------------------P List-----------------------" << endl;
	fout << "Q:" << endl;

	double Q[Mnum], Q_yuan[Mnum];
	
	ifstream Qfile("C:/Users/20148/Desktop/small case/CPLEX_programming1/input/Qi.txt");
	for (int i = 0; i < 13; i++)
	{
		Qfile >> Q_yuan[i];
	}
	for (int i = 0; i < Mnum; i++)
	{
		
		Q[i] = Q_yuan[(i % 13)];
		
		fout << "Q[" << i << "]=" << Q[i] << ";   ";
		
	}
	fout << endl;
	for (int i = 0; i < Mnum; i++)
	{
		fout << Q[i] << "  ";
	}
	Qfile.close();
	fout << endl << endl;

	int  c[Mnum][Nnum], C[13][14];
	double D[13][14]; 
	double f[Mnum], F[13], h[Mnum], H[13];
	fout << "Facility cost:" << endl;
	ifstream Ffile("C:/Users/20148/Desktop/small case/CPLEX_programming1/input/Fi.txt");
	for (int i = 0; i < 13; i++)
	{
		Ffile >> F[i];
	}
	for (int i = 0; i < Mnum; i++)
	{
		f[i] = F[(i % 13)];
		
		fout << "f[" << i << "]=" << f[i] << ";";
	}
	Ffile.close();
	fout << endl;
	for (int i = 0; i < Mnum; i++)
	{
		fout << f[i] << "  ";
	}
	fout << endl << endl;

	fout << "H cost:" << endl;
	ifstream Hfile("C:/Users/20148/Desktop/small case/CPLEX_programming1/input/Hi.txt");
	//Hi
	for (int i = 0; i < 13; i++)
	{
		Hfile >> H[i];
	}
	for (int i = 0; i < Mnum; i++)
	{
		h[i] = H[(i % 13)];
		
		fout << "h[" << i << "]=" << h[i] << ";";
		
	}fout << endl << endl;
	for (int i = 0; i < Mnum; i++)
	{
		fout << h[i] << "  ";
	}fout << endl;
	Hfile.close();


	fout << "D cost:" << endl;      
	ifstream Cfile("C:/Users/20148/Desktop/small case/CPLEX_programming1/input/Dij.txt");
	for (int k = 0; k < 13; k++) {
		for (int j = 0; j < 14; j++) {
			Cfile >> D[k][j];
		}
	}

	double Service_distance = 0.0;
	Service_distance = 60 * service_hours;

	for (int k = 0; k < 13; k++) {
		for (int j = 0; j < 14; j++) {
			if (D[k][j] <= Service_distance) {
				C[k][j] = 0;
			}
			else {
				C[k][j] = 10000000;
			}
		}
	}

	for (int i = 0; i < Mnum; i++)
	{
		for (int j = 0; j < Nnum; j++)
		{
			c[i][j] = C[(i % 13)][(j % 14)];
			
			fout << '\t' << "c[" << i << "][" << j << "]= " << c[i][j] << " ";
		}
		fout << endl;         
	}
	Cfile.close();


	fout << endl << "Dem inf" << endl;
	

	double DFF[14], DF[Nnum], miu[Nnum], miu2[Nnum], V[Nnum];     
	double cputime[cycle][2], Inventory[cycle][2], obj[cycle][2],
		FillrateN[cycle][2], chanceN[cycle][2], FillrateU[cycle][2], chanceU[cycle][2], FillrateT[cycle][2], chanceT[cycle][2];
	double improvement[cycle], improvement_total = 0.0;
	int OpenF[cycle][2], link[cycle][2];
	double multi_time = 0.0, multi_I = 0.0, multi_obj = 0.0,
		multi_fillU = 0.0, multi_chanceU = 0.0, multi_fillN = 0.0, multi_chanceN = 0.0, multi_fillT = 0.0, multi_chanceT = 0.0;
	double single_time = 0.0, single_facility = 0.0, single_links = 0.0, single_I = 0.0, single_obj = 0.0,
		single_fillrateU = 0.0, single_chanceU = 0.0, single_fillrateN = 0.0, single_chanceN = 0.0, single_fillrateT = 0.0, single_chanceT = 0.0;
	int multi_facility = 0, multi_links = 0;

	double multi_MPtime[cycle][2], MMP_time = 0.0;
	double multi_SPtime[cycle][2], MSP_time = 0.0;
	int multi_cut[cycle][2], Multi_cut = 0.0;

	double single_MPtime[cycle][2], SMP_time = 0.0;
	double single_SPtime[cycle][2], SSP_time = 0.0;
	int single_cut[cycle][2], Single_cut = 0.0;

	for (int k = 0; k < cycle; k++)
	{
		
		ifstream DMfile("C:/Users/20148/Desktop/small case/CPLEX_programming1/input/MeanDj.txt");
		fout << "**********************" << k << "times Dem inf*************************" << endl;
		
		for (int j = 0; j < 14; j++)
		{
			DMfile >> DFF[j];
		}
		for (int j = 0; j < Nnum; j++)
		{
			DF[j] = DFF[(j % 14)];    
			miu[j] = random_uniform_distribution(0.9 * DF[j], 1.1 * DF[j]);
			
			V[j] = random_uniform_distribution(50, 150);
			fout << "DF[" << j << "]=" << DF[j] << ";   ";
			fout << "miu[" << j << "]=" << miu[j] << ";   ";
			fout << "V[" << j << "]=" << V[j] << ";   ";
			fout << endl;
		}
		DMfile.close();
		double std[Nnum];
		for (j = 0; j < Nnum; j++)
		{
			std[j] = random_uniform_distribution(10, 100);
		}
		double sumDL = 0, sumDU = 0, sumDM = 0;
		double D[Nnum][100];
		double DL[Nnum];
		double DU[Nnum];
		double DM[Nnum], djw[Nnum];
		fout << endl;

		for (int j = 0; j < Nnum; j++)
		{
			for (int cpj = 0; cpj < 100; cpj++)
			{
				D[j][cpj] = random_normal_distribution_BoxMuller(miu[j], V[j]);
				
				if (D[j][cpj] <= 0) D[j][cpj] = 0.0;
				else continue;
			}
		}
		for (int j = 0; j < Nnum; j++)
		{
			DM[j] = findavg(D[j], 100);
			DL[j] = findmin(D[j], 100);
			DU[j] = findmax(D[j], 100);
			fout << "Node" << j << "average dem：" << DM[j] << endl;
			fout << "Node" << j << "interval：" << DL[j] << " to " << DU[j] << endl;
		}
		
		fout << "Lower：" << endl;
		for (int j = 0; j < Nnum; j++)
		{
			fout << DL[j] << "    ";
		}
		fout << endl << "Upper：" << endl;
		for (int j = 0; j < Nnum; j++)
		{
			fout << DU[j] << "    ";
		}
		fout << endl;
		for (int j = 0; j < Nnum; j++)
		{
			
			sumDM += DM[j];
			sumDU += DU[j];
			sumDL += DL[j];
		}
		fout << "Sum upper：" << sumDU << endl;
		fout << "Sum lower：" << sumDL << endl;
		fout << "Sum ave：" << sumDM << endl << endl;

		double he = 0.0, TI = 0.0;
		for (j = 0; j < Nnum; j++)
		{
			he += (DU[j] - DL[j]) * (DU[j] - DL[j]);
		}
		TI = sumDM + sqrt(-log(eps) * he / 2);
		
		fout << "TI = " << TI << endl << endl;
		

		ofstream dddout("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/doutDemoutput.txt");
		//////////////////////////////// normal //////////////////////////////////////////////////
		dddout << " normal " << endl;
		static double dn[realization_num][Nnum], du[realization_num][Nnum], dt[realization_num][Nnum];
		for (int k = 0; k < realization_num; k++)
		{
			for (j = 0; j < Nnum; j++)
			{
				dn[k][j] = double(random_normal_distribution_BoxMuller(miu[j], V[j]));

				if (dn[k][j] <= 0) dn[k][j] = 0.0;
				else continue;
				dddout << "dn[" << j << "][" << k << "]= " << dn[j][k] << " ";
			}
			dddout << endl;
		}

		result << "======================== " << k << " tiems result： ========================" << endl;

		//========================================Multi result======================================		
		double* solu = new double[10];
		double  multi_fh = 0.0, multi_lc = 0.0;

		Multi_solution(solu, f, h, c, Q, DM, DL, DU, sumDL, TI);
		result << "*****Multi*****" << endl;
		cputime[k][0] = solu[0];     multi_time += cputime[k][0];
		multi_MPtime[k][0] = solu[8]; MMP_time += multi_MPtime[k][0];
		multi_SPtime[k][0] = solu[9]; MSP_time += multi_SPtime[k][0];
		multi_cut[k][0] = solu[10];
		if (multi_cut[k][0] < 0.5) multi_cut[k][0] = 0;
		Multi_cut += multi_cut[k][0];
		

		ifstream OPin("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/open.txt");
		
		ifstream Yin("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/YI.txt");
		int Multi_Yi[Mnum];
		for (i = 0; i < Mnum; i++)
		{
			Yin >> Multi_Yi[i];
		}
		OpenF[k][0] = 0; Inventory[k][0] = 0.0;
		for (i = 0; i < Mnum; i++)
		{
			OpenF[k][0] += Multi_Yi[i];
		}
		multi_facility += OpenF[k][0];

		ifstream Iin("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/II.txt");
		double Multi_II[Mnum] = { 0.0 };
		for (i = 0; i < Mnum; i++)
		{
			Iin >> Multi_II[i];
			Inventory[k][0] += Multi_II[i];
		}
		multi_I += Inventory[k][0];
		for (i = 0; i < Mnum; i++)
		{
			multi_fh += f[i] * Multi_Yi[i] + h[i] * Multi_II[i];
			result << "  Y[" << i << "] = " << Multi_Yi[i];
			result << "  I[" << i << "] = " << Multi_II[i] << endl;
		}

		ifstream xin("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/XIJ.txt");
		int Multi_Xij[Mnum][Nnum];
		link[k][0] = 0;
		for (i = 0; i < Mnum; i++)
		{
			for (j = 0; j < Nnum; j++) {
				xin >> Multi_Xij[i][j];
			}
		}
		for (i = 0; i < Mnum; i++)
		{
			for (j = 0; j < Nnum; j++) {
				result << "  X[" << i << "][" << j << "]=" << Multi_Xij[i][j];
				
				multi_lc += c[i][j] * Multi_Xij[i][j];
			}
			result << endl;
		}result << endl << endl;
		
		for (i = 0; i < Mnum; i++)
		{
			for (j = 0; j < Nnum; j++)
			{
				
				if (Multi_Xij[i][j] > 0.1)
					link[k][0] += 1;
				else link[k][0] += 0;
				
			}
			
		}
		
		multi_links += link[k][0];
		
		obj[k][0] = multi_fh + multi_lc;         multi_obj += obj[k][0];
		Yin.close();
		Iin.close();
		xin.close();
		
		double* fill = new double[2];

		result << "Time :" << cputime[k][0] << endl;
		result << "MP time :" << multi_MPtime[k][0] << endl;
		result << "SP time :" << multi_SPtime[k][0] << endl;
		result << "cut number :" << multi_cut[k][0] << endl;
		result << "open facility :" << OpenF[k][0] << endl;
		result << "Arc :" << link[k][0] << endl;
		result << "Inv :" << Inventory[k][0] << endl;
		result << "Objective value :" << obj[k][0] << endl;
		
		

		FormulateScLP(ScLP, ScLPObj, ScLPCon1, ScLPCon3, Multi_Xij, Multi_II);
		
		FillRateLoop(fill, ScLP, ScLPCon2, ScLPSolver, dn);
		FillrateN[k][0] = fill[0]; multi_fillN += FillrateN[k][0];
		chanceN[k][0] = fill[1]; multi_chanceN += chanceN[k][0];
		result << "FillRateN :" << FillrateN[k][0] << endl;
		result << "ChanceN :" << chanceN[k][0] << endl;

		result << ' ' << '\t' << "CPUtime" << '\t' << "MPtime" << '\t' << "SPtime" << '\t' << "setupF" << '\t' << "Arcs" << '\t' << "Inv" << '\t' << "Obj" << '\t' << "cut"
			<< '\t' << "FN" << '\t' << "CN" << endl;
		result << "Ins" << k << '\t' << cputime[k][0] << '\t' << multi_MPtime[k][0] << '\t' << multi_SPtime[k][0] << '\t' << OpenF[k][0] << '\t' << link[k][0] << '\t'
			<< Inventory[k][0] << '\t' << obj[k][0] << '\t' << multi_cut[k][0]
			<< '\t' << FillrateN[k][0] << '\t' << chanceN[k][0] << endl;
		
		delete[] solu;
		
		//========================================Single result======================================
		double* single_solu = new double[10];
		double  single_fh = 0.0, single_lc = 0.0;
		Single_solution(single_solu, f, h, c, Q, DM, DL, DU, sumDL, TI);
		result << "*****Single*****" << endl;
		cputime[k][1] = single_solu[0];     single_time += cputime[k][1];

		single_MPtime[k][0] = single_solu[8]; SMP_time += single_MPtime[k][0];
		single_SPtime[k][0] = single_solu[9]; SSP_time += single_SPtime[k][0];
		single_cut[k][0] = single_solu[10];
		if (single_cut[k][0] < 0.5) single_cut[k][0] = 0;
		Single_cut += single_cut[k][0];
		

		ifstream single_OPin("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/single_open.txt");
		single_OPin >> OpenF[k][1];

		ifstream single_Yin("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/single_YI.txt");
		double single_Yi[Mnum] = { 0.0 };
		for (i = 0; i < Mnum; i++)
		{
			single_Yin >> single_Yi[i];
		}
		ifstream single_Iin("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/single_II.txt");
		double single_II[Mnum] = { 0.0 };
		Inventory[k][1] = 0.0;
		for (i = 0; i < Mnum; i++)
		{
			single_Iin >> single_II[i];
			Inventory[k][1] += single_II[i];
		}
		single_I += Inventory[k][1];
		OpenF[k][1] = 0;
		for (i = 0; i < Mnum; i++)
		{
			result << "  Y[" << i << "] = " << single_Yi[i];
			OpenF[k][1] += single_Yi[i];
			result << "  I[" << i << "] = " << single_II[i] << endl;
			single_fh += f[i] * single_Yi[i] + h[i] * single_II[i];
		}
		single_facility += OpenF[k][1];
		single_Yin.close();
		single_Iin.close();
		ifstream single_xin("C:/Users/20148/Desktop/small case/CPLEX_programming1/output/single_XIJ.txt");
		int single_Xij[Mnum][Nnum];
		link[k][1] = 0;
		for (i = 0; i < Mnum; i++)
		{
			for (j = 0; j < Nnum; j++) {
				single_xin >> single_Xij[i][j];
			}
		}
		for (i = 0; i < Mnum; i++)
		{
			for (j = 0; j < Nnum; j++) {
				result << "  X[" << i << "][" << j << "]=" << single_Xij[i][j];
				link[k][1] += single_Xij[i][j];
				single_lc += c[i][j] * single_Xij[i][j];
			}
			result << endl;
		}
		//		-----------Objective value----------;
		obj[k][1] = single_fh + single_lc;         single_obj += obj[k][1];
		single_links += link[k][1];
		single_xin.close();

		double* single_fill = new double[2];

		result << "Time :" << cputime[k][1] << endl;
		result << "MP time :" << single_MPtime[k][0] << endl;
		result << "SP time :" << single_SPtime[k][0] << endl;
		result << "cut number :" << single_cut[k][0] << endl;
		result << "Open facility :" << OpenF[k][1] << endl;
		result << "Arc :" << link[k][1] << endl;
		result << "Inv :" << Inventory[k][1] << endl;
		result << "Objective value :" << obj[k][1] << endl;
		

		FormulateScLP(ScLP, ScLPObj, ScLPCon1, ScLPCon3, single_Xij, single_II);
		
		FillRateLoop(fill, ScLP, ScLPCon2, ScLPSolver, dn);
		FillrateN[k][1] = fill[0]; single_fillrateN += FillrateN[k][1];
		chanceN[k][1] = fill[1]; single_chanceN += chanceN[k][1];
		result << "FillRateN :" << FillrateN[k][1] << endl;
		result << "chanceN :" << chanceN[k][1] << endl;

		result << ' ' << '\t' << "CPUtime" << '\t' << "MPtime" << '\t' << "SPtime" << '\t' << "setupF" << '\t' << "Arcs" << '\t' << "Inv" << '\t' << "Obj" << '\t' << "cut"
			<< '\t' << "FN" << '\t' << "CN" << endl;

		result << "Ins" << k << '\t' << cputime[k][1] << '\t' << single_MPtime[k][0] << '\t' << single_SPtime[k][0] << '\t' << OpenF[k][1] << '\t' << link[k][1] << '\t'
			<< Inventory[k][1] << '\t' << obj[k][1] << '\t' << single_cut[k][0]
			<< '\t' << FillrateN[k][1] << '\t' << chanceN[k][1] << endl;
		
		delete[] single_solu;
		delete[] single_fill;
		
		delete[] fill;


	}   
	result << endl << endl << endl;
	//======================================Average result================================================
	result << "======================== Average result ========================" << endl;
	result << "*****Multi*****" << endl;
	result << ' ' << '\t' << "CPUtime" << '\t' << "MPtime" << '\t' << "SPtime" << '\t' << "setupF" << '\t' << "Arcs" << '\t' << "Inv" << '\t' << "Obj" << '\t' << "cut"
		<< '\t' << "FN" << '\t' << "CN" << endl;
	
	for (int k = 0; k < cycle; k++)
	{
		result << "Ins" << k << '\t' << cputime[k][0] << '\t' << multi_MPtime[k][0] << '\t' << multi_SPtime[k][0] << '\t' << OpenF[k][0] << '\t' << link[k][0] << '\t'
			<< Inventory[k][0] << '\t' << obj[k][0] << '\t' << multi_cut[k][0]
			<< '\t' << FillrateN[k][0] << '\t' << chanceN[k][0] << endl;
	}

	result << "Ave" << '\t' << multi_time / Acycle << '\t' << MMP_time / Acycle << '\t' << MSP_time / Acycle << '\t' << multi_facility / Acycle << '\t' << multi_links / Acycle << '\t'
		<< multi_I / Acycle << '\t' << multi_obj / Acycle << '\t' << Multi_cut / Acycle
		<< '\t' << multi_fillN / Acycle << '\t'
		<< multi_chanceN / Acycle << endl << endl;


	result << "Ave CPUtime = " << multi_time / Acycle << endl;
	result << "Ave MPtime :" << MMP_time / Acycle << endl;
	result << "Ave SPtime :" << MSP_time / Acycle << endl;
	result << "Ave cut number :" << Multi_cut / Acycle << endl;
	result << "  Total open facility = " << multi_facility << "  Ave open facility = " << multi_facility / Acycle << endl;
	result << "  Total Arc = " << multi_links << "  Ave Arc = " << multi_links / Acycle << endl;
	result << "  Total Inv = " << multi_I << "  Ave Inv = " << multi_I / Acycle << endl;
	result << "  Total objective value :" << multi_obj << "  Ave objective value :" << multi_obj / Acycle << endl;

	result << "Ave FillRateN = " << multi_fillN / Acycle << endl;
	result << "Ave chanceN = " << multi_chanceN / Acycle << endl;
	


	result << "*****Single*****" << endl;
	result << "*****Single*****" << endl;
	result << ' ' << '\t' << "CPUtime" << '\t' << "MPtime" << '\t' << "SPtime" << '\t' << "setupF" << '\t' << "Arcs" << '\t' << "Inv" << '\t' << "Obj" << '\t' << "cut"
		<< '\t' << "FN" << '\t' << "CN" << endl;

	for (int k = 0; k < cycle; k++)
	{
		result << "Ins" << k << '\t' << cputime[k][1] << '\t' << single_MPtime[k][0] << '\t' << single_SPtime[k][0] << '\t' << OpenF[k][1] << '\t' << link[k][1] << '\t'
			<< Inventory[k][1] << '\t' << obj[k][1] << '\t' << single_cut[k][0]
			<< '\t' << FillrateN[k][1] << '\t' << chanceN[k][1] << endl;
	}

	result << "Ave" << '\t' << single_time / Acycle << '\t' << SMP_time / Acycle << '\t' << SSP_time / Acycle << '\t' << single_facility / Acycle << '\t' << single_links / Acycle << '\t'
		<< single_I / Acycle << '\t' << single_obj / Acycle << '\t' << Single_cut / Acycle
		<< '\t' << single_fillrateN / Acycle << '\t'
		<< single_chanceN / Acycle << endl << endl;


	result << "Ave CPUtime = " << single_time / Acycle << endl;
	result << "Ave MPtime :" << SMP_time / Acycle << endl;
	result << "Ave SPtime :" << SSP_time / Acycle << endl;
	result << "Ave cut number :" << Single_cut / Acycle << endl;

	result << "  Total open facility = " << single_facility <<
		"  Ave open facility = " << single_facility / Acycle << endl;
	result << "  Total Arc = " << single_links <<
		"  Ave Arc = " << single_links / Acycle << endl;
	result << "  Total Inv = " << single_I <<
		"  Ave Inv = " << single_I / Acycle << endl;
	result << "  Total objective value :" << single_obj <<
		"  Ave objective value :" << single_obj / Acycle << endl;
	result << "  Total FillRateN = " << single_fillrateN <<
		"  Ave FillRateN = " << single_fillrateN / Acycle << endl;
	result << "  Total chanceN = " << single_chanceN <<
		"；  Ave chanceN = " << single_chanceN / Acycle << endl << endl;

	result << "EPS = " << eps << endl << endl;
	
	return 0;
}


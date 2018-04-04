#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <map>
using namespace std;

double randnu ();
double roundin (double);
int rand_l_lim(int comp_size);
int rand_w_lim(int comp_size);
int rand_tox_lim(int comp_size);
void setzero(double &val);//??
void carry_out_reaction();//??
int kil;
int num_defects;
int end_of_str_traps;
double k31, k13, mul_fac;
double E1, E2, R, SHW;
double compartment_length;
double z_comp_length;
int num_compartment, domains_per_nm;
int j;
int nr1, nr2, nf1, nf2, ndz1, ndz2;
double q1, q2, q3;
int c1, c2, c3, c4;
int k1, k2, k3;
double somevar, somevar2, somevar3;
double alpha0, yemp, tau;
double kc1,kc2;
double temp13, temp11, temp12;
double width, length, height;
double delVot = 0,delVt_prev;
int w_lim, l_lim, z_lim;

#define WIDTH  50
#define LENGTH 20
#define HEIGHT 1

#define fpsi(VGS,VFB,psi,bot) VGS - VFB - psi - bot*sqrt( (2*q*ESi*E0*ND) * (VT*exp(-psi/VT) + psi - VT + exp(-2*phif/VT)*(VT*exp(psi/VT) - psi - VT)))/cox
#define gaussian(x,mu,sigma) exp(- pow((x-mu)/sigma,2)/2)/(sqrt(2*pi)*sigma)

double alpha_k13_total;
double alpha_k31_total;
double temp_1, temp_2;
double err;
vector < pair <double, double> > E1_arr, E2_arr, SHW_arr, R_arr;
vector <double> E1_store, E2_store, SHW_store, R_store, gauss_arr;

double psi_str, FSiO2;
int emi, emiprev, previous_emi, previous_emiH2;
double str_Nit[10000];
double str_time[10000];
double t,t_prev,t_base;
int prev,curr;
int counter, f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16;
double r, r1, r2;

double kk1, kk2, kk3;

struct threeD {
	threeD(int x, int y, int z) : x(x), y(y), z(z) {}
	threeD(const threeD &c) {x = c.x; y = c.y; z = c.z;}
	int x;
	int y;
	int z;

    bool operator==(const threeD &o) const {
     return ((x == o.x) && (z == o.z) && (y == o.y));
     }

	bool operator<(const threeD &o)  const {
		return ((x < o.x) || ((x == o.x) && (y < o.y)) || ((x == o.x) && (y == o.y) && (z < o.z)));
	}
};

struct locate_threeD {
    int a;
    int b;
    int c;

    locate_threeD(int a, int b, int c) : a(a), b(b), c(c) {}

    bool operator () (const threeD &o) const{
        return (o.x==a && o.y==b && o.z==c);
    }
};

vector <threeD> init_sites;
vector < double > k13_arr, k31_arr;
multimap<int, threeD> bulk_defects;
map <threeD, pair<double, double> > defect_rates;

/**********************************************************************************/
/********************************    Device Data     ******************************/
/**********************************************************************************/

const double q = 1.6e-19; // Farad, charge on an electron
const double ESi = 11.8; // Relative permittivity of Silicon
const double ESiO2 = 3.9; // Relative permittivity of Silicon-di-oxide
const double E0 = 8.854e-14; // Absolute permittivity in farad/cm
const double ni0 = 1e10; // Intrensic carrier concentration at 300 K in /cm^3
const double Eg = 1.12; // Bandgap of Si at 300 K in eV
const double Kb= 1.380658e-23; // Boltzman constant in J/K

/**********************************************************************************/

//// Changeable parameters
    double Temp = 130;
    double Vgs_dc = 2; //DC stess voltage in Volts
    double Vgr_dc = 0; //DC recovery voltage in Volts
	double Tox = HEIGHT*1e-7;             //cm
	double dev_area = WIDTH*LENGTH*1E-14; // cm2
	double ND = 1e18; //  Channel Doping in atoms/cm^3
    double N0 = 5e20; //  Density of gate insulator traps in /cm^3
    double vfb = 0.4992; // Flatband voltage in Volts
    double E1m = -0.3;     // eV
    double E1s = 0.2;        // eV
    double E2m = 0;         // eV
    double E2s = 0.05;        // eV
    double EV = 0;         // eV
	double SHWm = 3;      // 5.9eV
	double SHWs = 0.5;       // eV
	double Rm = 0.5;       //0.52
	double Rs = 0;         // eV
    double Gamma = 2e-8; // eV cm/V
/*********************************************************************************/

	double VT1 = (Kb*300)/q;    // Thermal voltage Vt@ 300 K = 0.0259eV in eV
	double VT= (Kb*(273+Temp))/q;   // Thermal voltage Vt@ 273+Temp in eV
	double ni = ni0*(pow(static_cast<double>((273+Temp)/300),1.5))*exp(-Eg/2*(1/VT-1/VT1)); // in cm^-3
	double phif=VT*log(ND/ni);                // Bulk potential in V
	double beta = 1/VT;                       // in 1/eV ;34.7meV for 130 deg C

	double kt=Kb*(273+Temp);// Joules
	double cox=(ESiO2*E0)/Tox; // in F/cm^2
	double Ldp=sqrt((ESi*E0*kt)/(q*q*ND)); // Debye Length in cm

    double pi=22.0/7;
	double mp=0.5*9.1e-31;                   // hole mass in kg (http://ecee.colorado.edu/~bart/book/effmass.htm)
    double vp=sqrt(8*kt/(pi*mp))*100;         // should be thermal velocity (check the formula) in cm/s
    double sp=(1e-15) ;                       // cross-section area in cm^2
    double h=6.62607004e-34;                  // in Joule*seconds
    double NV=2*(pow(static_cast<double>(2*pi*mp*Kb*(273+Temp)/(pow(h,2))),1.5))*1e-6; // in /cm^3
    double p0=(pow(ni,2))/ND;                 // Bulk minority concentration  /cm^3

void generate_rand_params()
{
    double g1;
    if(E1s==0) E1 = E1m;
    else
    {
        g1 = randnu();
        for(int i = 0; i < E1_arr.size(); i++)
        {
            if(g1 <= E1_arr[i].second)
            {
                E1 = abs(E1_arr[i].first) < 1e-8 ? 0 : E1_arr[i].first;
                break;
            }
        }
    }

    if(E2s==0) E2 = E2m;
    else
    {
        g1 = randnu();
        for(int i = 0; i < E2_arr.size(); i++)
        {
            if(g1 <= E2_arr[i].second)
            {
                E2 = abs(E2_arr[i].first) < 1e-8 ? 0 : E2_arr[i].first;
                break;
            }
        }

    }

    if(SHWs==0) SHW = SHWm;
    else
    {
        g1 = randnu();
        for(int i = 0; i < SHW_arr.size(); i++)
        {
            if(g1 <= SHW_arr[i].second)
            {
                SHW = abs(SHW_arr[i].first) < 1e-8 ? 0 : SHW_arr[i].first;
                break;
            }
        }
    }

    if(Rs==0) R = Rm;
    else
    {


        g1 = randnu();
        for(int i = 0; i < R_arr.size(); i++)
        {
            if(g1 <= R_arr[i].second)
            {
                R = abs(R_arr[i].first) < 1e-8 ? 0 : R_arr[i].first;    //Check 1e-8.Might want to adjust for different values of sigma.
                break;
            }
        }
    }
}


void ox_field(double VGb_str)
{

//    if(VGb_str == 1.5) FSiO2 = 10e6;
//    else if(VGb_str == 0) FSiO2 = 1e6;
	//Solving for surface potential
	int cnt=1;

	int maxiter = 2500;
	double lpsi;
	double rpsi;
    int b1;

	if(VGb_str >= vfb)
    {
        lpsi = 0;
        rpsi = 3*phif;
        b1 = 1;
    }

    else
    {
        lpsi = -2*phif;
        rpsi = 0;
        b1 = -1;
    }

	double fvall = fpsi(VGb_str,vfb,lpsi,b1);
    //cout<<"fvall "<<fvall<<endl;

	double fvalu = fpsi(VGb_str,vfb,rpsi,b1);
    //cout<<"fvalu "<<fvalu<<endl;
    //cout<<VGb_str<<endl;
    double newpsi;
	double tmp=0.0;
	double tol=1e-10;
	double fvalleft,fvalright;

	if (fvall*fvalu < 0.0)
	{
	    //cout<<"hello"<<endl;
		newpsi = (lpsi+rpsi)/2;
		do {
			fvalleft  = fpsi(VGb_str,vfb,lpsi,b1);
			fvalright = fpsi(VGb_str,vfb,newpsi,b1);

			if (fvalleft*fvalright < 0)
					rpsi=newpsi;
				else
					lpsi=newpsi;

            tmp = (lpsi+rpsi)/2;
            newpsi=tmp;
            cnt=cnt+1;
		}

		while(cnt < maxiter && fabs(fvalright) > tol );
		psi_str=newpsi;
	}

    //cout<<psi_str<<endl;

	//Surface Electric Field
	double FSi = sqrt(2)*(VT/Ldp)*sqrt((exp(-psi_str/VT)+(psi_str/VT)-1)+((exp(-2*phif/VT))*(exp(psi_str/VT)-(psi_str/VT)-1)));

	//Electric Field in SiO2
	FSiO2 = FSi*ESi/ESiO2;
	//cout<<FSiO2<<endl;
}

void rate(double xvecbytox)
{
    double ps=p0*(exp(beta*psi_str)); // Surface minority carrier concentration near the interface which is inversion carrier density in  /cm^3
    //cout<<psi_str<<endl;
    double pf_ps = ps*vp*sp;
    double pf_nv = NV*vp*sp;
    //cout<<"pf_ps "<<pf_ps<<endl;
    //cout<<"pf_nv "<<pf_nv<<endl;

    //cout<<"NV "<<NV<<", ps "<<ps<<endl;

    double delET = E1-EV;

    // CHECK SIGN!
    double epsT1_str = E2 - delET - Gamma*(1-xvecbytox)*FSiO2;   // Need to check sign

    //cout<<FSiO2<<endl;
    //cout<<"epsT1_str "<<epsT1_str<<endl;
    double eps12d_str = SHW/(pow((1 + R),2)) + R*epsT1_str/(1+R) + R*(pow(static_cast<double>(epsT1_str),2))/(4*SHW);
    //double eps21d_str = SHW/(pow((1 + R),2)) - epsT1_str/(1+R) + R*(pow(static_cast<double>(epsT1_str),2))/(4*SHW);
    double eps21d_str = eps12d_str - epsT1_str;

    //cout<<pf_ps<<","<<pf_nv<<endl;
    //cout<<beta<<endl;

    //cout<<"eps21d_str "<<eps21d_str<<endl;
    //cout<<"eps12d_str "<<eps12d_str<<endl;

    k13 = pf_ps*exp(-beta*eps12d_str);
    //cout<<beta<<endl;
		if(k13 > pf_ps)
			k13 = pf_ps;

    k31 = pf_nv*exp(-beta*eps21d_str);
		if(k31 > pf_nv)
			k31 = pf_nv;


    //k13 = 0.2;
    //k31 = 0.8;
}

int main(){


	srand (99981);
	//cout<<rand()<<endl;
	double sw_time_1 [] = {1e-6, 1e3};

    //****************Device Dimensions*********************//

	num_defects = (int) LENGTH*WIDTH*HEIGHT*1/2;

	compartment_length = 0.5e-9;
	z_comp_length = 0.5e-9;
	mul_fac = 1/(compartment_length*compartment_length*z_comp_length*1e6);

	float modder = (1e13)/(num_defects);

	width = WIDTH*1e-9;
	w_lim = (int) ceil(width/compartment_length);
	length = LENGTH*1e-9;
	l_lim = (int) ceil(length/compartment_length);
	height = HEIGHT*1e-9;
	z_lim =  (int) ceil(height/z_comp_length);


    //******************************************************//


	bool hit1 = 0;
	bool hit2 = 0;
	bool hit3 = 0;
	bool hit4 = 0;
	bool hit5 = 0;

	E1_store.resize(num_defects);
	E2_store.resize(num_defects);
	R_store.resize(num_defects);
    SHW_store.resize(num_defects);

	ofstream myfile1,myfile,myfile2,myfile3;
	// Uncomment for writing in single text file
	//myfile1.open ("Rec.txt", ios::out | ios::trunc);
	//myfile.open ("Str.txt", ios::out | ios::trunc);

	myfile1.open ("Rec.csv", ios::out | ios::trunc);
	myfile.open ("Str.csv", ios::out | ios::trunc);

    myfile2.open ("Str_data.txt", ios::out | ios::trunc);
    myfile3.open ("Rec_data.txt", ios::out | ios::trunc);
    multimap <int, threeD>::iterator itstr,itrec;
    map<threeD, pair <double, double> >::iterator itrate;

    int cn = 0;
	int flagger = 0;
    int d1 = 0;

//**********************************************************************
// PDFs for all 4 Gaussian distributions

// Store probabilities of Gaussian distributed RV from mu-3*sigma to mu+3*sigma
// mu = 1 and sigma=1
somevar = -2;
somevar2 = 0;
for(int i = 0; i < 20; i++)
{
    somevar += 0.3;
    somevar2 += gaussian(somevar,1,1);
    gauss_arr.push_back(somevar2);
}

for(int i = 0; i<20; i++)
{
    gauss_arr[i] = gauss_arr[i]/somevar2;
    //cout<<gauss_arr[i]<<endl;
}

if(E1s != 0)
{
    somevar = E1m-3*E1s;
    somevar3 = 6*E1s/20;

    for(int i = 0; i<20; i++)
    {
        somevar += somevar3;
        E1_arr.push_back(make_pair(somevar,gauss_arr[i]));
    }
}

if(E2s != 0)
{
    somevar = E2m-3*E2s;
    somevar3 = 6*E2s/20;

    for(int i = 0; i<20; i++)
    {
        somevar += somevar3;
        E2_arr.push_back(make_pair(somevar,gauss_arr[i]));
    }
}

if(SHWs != 0)
{
    somevar = SHWm-3*SHWs;
    somevar3 = 6*SHWs/20;

    for(int i = 0; i<20; i++)
    {
        somevar += somevar3;
        SHW_arr.push_back(make_pair(somevar,gauss_arr[i]));
    }
}

if(Rs != 0)
{
    somevar = Rm-3*Rs;
    somevar3 = 6*Rs/20;

    for(int i = 0; i<20; i++)
    {
        somevar += somevar3;
        R_arr.push_back(make_pair(somevar,gauss_arr[i]));
    }
}

t_prev = sw_time_1[0];
while(t_prev <= sw_time_1[1])
{
    myfile<<t_prev<<",";
    myfile1<<t_prev<<",";
    t_prev = t_prev + pow(10,floor( log10(t_prev)));
}

myfile<<t_prev<<endl;
myfile1<<t_prev<<endl;

for (int ii = 0; ii < 2; ii++)
{


//**********************************************************************

    //Rate constants for Stress
    //k13 = 0.001;
	//k31 = 0.00001;

	hit1 = 0;
	hit2 = 0;
	hit3 = 0;
	hit4 = 0;
	hit5 = 0;

	d1 =0;
    alpha_k13_total = 0;                   // Initial total forward propensity: from ground state to transport state
    alpha_k31_total = 0;                  // Initial reverse propensity:

//  Distribute defect sites spatially in the bulk
	for(d1= 0; d1<num_defects; d1++){
		c1 = rand_l_lim(l_lim);
		c2 = rand_w_lim(w_lim);
		c3 = rand_tox_lim(z_lim);

		vector<threeD>::iterator it;
		it = find_if(init_sites.begin(), init_sites.end(), locate_threeD(c1,c2,c3));
		if (it != init_sites.end())
        {
			d1--;
			continue;
		}

		else
        {
			init_sites.push_back(threeD(c1,c2,c3));
			defect_rates.insert(make_pair(threeD(c1,c2,c3), pair<double, double> (0,0)));
		}
	}


//*********************************************************************************
//******************************     STRESS     ***********************************
//*********************************************************************************


	emi = 0;
    emiprev = 0;
    delVot = 0;
    delVt_prev = 0;
	t = sw_time_1[0];

	counter = 0;
//	ofstream myfile1,myfile,myfile2,myfile3;
//	// Uncomment for writing in single text file
//	//myfile1.open ("Rec.txt", ios::out | ios::trunc);
//	//myfile.open ("Str.txt", ios::out | ios::trunc);
//	//
//
//	  myfile1.open ("Rec.csv", ios::out | ios::app);
//	  myfile.open ("Str.csv", ios::out | ios::app);
//
//    myfile2.open ("Str_data.txt", ios::out | ios::trunc);
//    myfile3.open ("Rec_data.txt", ios::out | ios::trunc);
//    multimap <int, threeD>::iterator itstr,itrec;
//    map<threeD, pair <double, double> >::iterator itrate;
	flagger =0;


    ox_field(Vgs_dc);
    //kc1 = 0.001;
    //kc2 = 0.009;
    //cout<<"Field "<<FSiO2<<endl;
    //cout<<"pf_ps "<<pf_ps<<endl;
    //cout<<"pf_nv "<<pf_nv<<endl;

    //cout<<psi_str<<endl;
    cn = 0;
    for(itrate = defect_rates.begin(); itrate!=defect_rates.end(); itrate++)
    {
        generate_rand_params();
        E1_store[cn] = E1;
        E2_store[cn] = E2;
        SHW_store[cn] = SHW;
        R_store[cn] = R;
        cn += 1;
//        cout<<"E1 "<<E1<<endl;
//        cout<<"E2 "<<E2<<endl;
//        cout<<"Shw "<<SHW<<endl;
//        cout<<"R "<<R<<endl;

        rate((double) (itrate->first.z)/z_lim);
//        if(k31>10000)
//        cout<<k13<<" "<<k31<<endl;

        (itrate->second).first = k13;
        //cout<<k13<<" ";
        alpha_k13_total += k13;
        (itrate->second).second = k31;
        //cout<<k31<<endl;
    }
    //cout<<E1_store[0]<<" "<<E1_store[20]<<" "<<E1_store[40]<<endl;
    cout<<alpha_k13_total<<endl;


//    cout<<"E1_str "<<E1<<endl;
//    cout<<"E2_str "<<E2<<endl;
//    cout<<"SHW_str "<<SHW<<endl;
//    cout<<"R_str "<<R<<endl;
//    cout<<"field "<<FSiO2<<endl;
//    cout<<"surface potential "<<psi_str<<endl;

//    cout<<"k31_str "<<k31<<endl;
//    cout<<"k13_str "<<k13<<endl;

//    cout<<"Field "<<FSiO2<<endl;
//    cout<<(int) -0.9<<endl;
	while (t < sw_time_1[1] && flagger == 0){


		r1 = randnu();
		while (r1 == 0){
			r1 = randnu();
		}

		alpha0 =  alpha_k13_total + alpha_k31_total;
//		cout<<"alpha0 "<<alpha0<<endl;
//		cout<<"alpha_k13_total "<<alpha_k13_total<<endl;
//		cout<<"alpha_k31_total "<<alpha_k31_total<<endl;
		temp_1 = alpha_k31_total/alpha0;
		temp_2 = temp_1 + (alpha_k13_total/alpha0);

		tau = (1/alpha0)*log(1/r1);
		if (tau < 0){
			flagger = 1;
			myfile << "Error"<<endl;
		}

        if(t+tau <= sw_time_1[1])
        {
            r2 = randnu();
            while (r2 ==1){
                r2 = randnu();
            }

            carry_out_reaction();
        }

        if (t == sw_time_1[0])
        {
            t_prev = pow(10,floor(log10(t)));
            //cout<<tau;
            //myfile<<t_prev<<' '<<emiprev<<endl;
            //t_prev = t_prev*pow(10,0.1);
        }
            t = t + tau;

            while (t_prev < t && t_prev <= sw_time_1[1])

            {
                // Uncomment when writing in single file
                //myfile<<t_prev<<' '<<emiprev<<endl;

                //myfile<<emiprev<<",";
                myfile<<delVt_prev<<",";
                t_prev = t_prev + pow(10,floor( log10(t_prev))); // linear scale, check at t = 1


                //t_prev = t_prev*pow(10,0.1);                    //log scale
                //if(t_prev==1) cout<<log10(t_prev);
            }

            if(t <= sw_time_1[1])
            {
                //myfile<<t<<' '<<emi<<endl; //Comment this out if you need evenly spaced time vector from e-6 to e3

                emiprev = emi;
                delVt_prev = delVot;
            }


//        else
//        {
//            t = t + tau;
//            //t_base = floor(log10(t)*10)/10;
//            while (t_prev < t && t_prev < sw_time_1[1])
//            {
//                myfile<<t_prev<<' '<<emiprev<<endl;
//                //t_prev = t_prev + pow(10,floor(log10(t_prev)));   //linear scale
//                t_prev = t_prev*pow(10,0.1);                        //log scale
//            }
//
//            //
//            if(t < sw_time_1[1])           // Doesn't work for t = 1000
//            {
//                myfile<<t<<' '<<emi<<endl;
//                emiprev = emi;
//            }
////            else
////            {
////                myfile<<sw_time_1[1]<<' '<<emi<<endl;
////                //emi--;      // t>1000 had already been updated and needs to be reverted
////                cout<<((t_prev)<=sw_time_1[1])<<endl;
////            }


//        }

		counter = counter + 1;
	}

//    while(t_prev < sw_time_1[1])
//    {
//        myfile<<t_prev<<' '<<emi<<endl;
//        //t_prev = t_prev + pow(10,floor(log10(t_prev)));   //linear scale
//        t_prev = t_prev*pow(10,0.1);                        //log scale
//    }

    // Uncomment when writing in single file
    //myfile<<sw_time_1[1]<<' '<<emi<<endl;


      //myfile<<emi<<endl;
      myfile<<delVot<<endl;



	end_of_str_traps = emi;     //CHECK!

	myfile2 << emi <<" out of "<< num_defects<<" defects broken in bulk"<<endl;

	for (itstr = bulk_defects.begin(); itstr != bulk_defects.end(); itstr++) {
		myfile2<<(itstr->second).x<<" "<<(itstr->second).y<<" "<<(itstr->second).z<<endl;
	}
	myfile2<<endl;

//*********************************************************************************
//****************************     RECOVERY     ***********************************
//*********************************************************************************

	t = sw_time_1[0];
//  fill(i,i+sizeof(i),2);

	counter = 0;



    ox_field(Vgr_dc);
//    kc1 = 0.001;
//    kc2 = 0.05;
//    cout<<psi_str<<endl;
//    cout<<FSiO2<<endl;
      cn = 0;
    for(itrate = defect_rates.begin(); itrate != defect_rates.end(); itrate++)
    {
        // This can be done since ordering of index of map won't change
        E1 = E1_store[cn];
        E2 = E2_store[cn];
        SHW = SHW_store[cn];
        R = R_store[cn];
        cn += 1;

        // Uodate forward and backwrd rates at possible defect sites
        rate((double) (itrate->first.z)/z_lim);
        (itrate->second).first = k13;
        (itrate->second).second = k31;
    }

    //cout<<E1_store[0]<<" "<<E1_store[20]<<" "<<E1_store[40]<<endl;
    //cout<<E1_store.size()<<endl;

    alpha_k13_total = 0;
    alpha_k31_total = 0;

// After rates have been updated, recalculation of alphak13_tot and alphak31 is required

// alpha_k31_total is updated from existing bulk defects
    for(map<int,threeD>::iterator o = bulk_defects.begin(); o != bulk_defects.end(); o++)
    {
        alpha_k31_total += defect_rates[o->second].second;
    }
    //cout<<alpha_k31_total<<endl;
    //cout<<"k31_rec "<<k31<<endl;

// alpha_k13_total is updated from uncharged init_sites
    for(vector<threeD>::iterator o = init_sites.begin(); o != init_sites.end(); o++)
    {
        alpha_k13_total += defect_rates[*o].first;
    }
    //cout<<alpha_k13_total<<endl;
    //cout<<"k13_rec "<<k13<<endl;

    while (t < sw_time_1[1] && flagger == 0){

		r1 = randnu();
		while (r1 ==0){
			r1 = randnu();
		}

		alpha0 =  alpha_k13_total+ alpha_k31_total;
        temp_1 = alpha_k31_total/alpha0;
		temp_2 = temp_1 + (alpha_k13_total/alpha0);

		tau = (1/alpha0)*log(1/r1);
		if (tau < 0){
			flagger = 1;
			myfile << "Error"<<endl;
		}

        if(t+tau <= sw_time_1[1])
        {
            r2 = randnu();
            while (r2 ==1){
                r2 = randnu();
            }

            carry_out_reaction();

        }

        if (t == sw_time_1[0])
        {
            t_prev = pow(10,floor(log10(t)));
            //myfile1<<t_prev<<' '<<emiprev<<endl;
            //t_prev = t_prev*pow(10,0.1);
        }

            t = t + tau;

            while(t_prev < t && t_prev < sw_time_1[1])
            {
                // Uncomment for single file
                //myfile1<<t_prev<<' '<<emiprev<<endl;


                //myfile1<<emiprev<<",";
                myfile1<<delVt_prev<<",";
                t_prev = t_prev + pow(10,floor(log10(t_prev))); // linear scale


                //t_prev = t_prev*pow(10,0.1);                      //log scale

            }

            if(t < sw_time_1[1])
            {
                //myfile1<<t<<' '<<emi<<endl; //Comment this out if you need evenly spaced time vector from 1e-6 to 1e3

                delVt_prev = delVot;
                emiprev = emi;
            }


		counter = counter + 1;

		if (t > 1e-5 && hit1 == 0) {
			hit1 = 1;
			myfile3 <<"t = 1e-5"<<endl;
			for (itrec = bulk_defects.begin(); itrec != bulk_defects.end(); itrec++) {
				myfile3<<(itrec->second).x<<" "<<(itrec->second).y<<" "<<(itrec->second).z<<endl;
			}
			myfile3<<endl;
			myfile3<<(end_of_str_traps - emi)<<" defects recovered"<<endl;
		}

		if (t > 1e-3 && hit2 == 0) {
			hit2 = 1;
			myfile3 <<"t = 1e-3"<<endl;
			for (itrec = bulk_defects.begin(); itrec != bulk_defects.end(); itrec++) {
				myfile3<<(itrec->second).x<<" "<<(itrec->second).y<<" "<<(itrec->second).z<<endl;
			}
			myfile3<<endl;
			myfile3<<(end_of_str_traps - emi)<<" defects recovered"<<endl;
		}

		if (t > 1e-1 && hit3 == 0) {
			hit3 = 1;
			myfile3 <<"t = 1e-1"<<endl;
			for (itrec = bulk_defects.begin(); itrec != bulk_defects.end(); itrec++) {
				myfile3<<(itrec->second).x<<" "<<(itrec->second).y<<" "<<(itrec->second).z<<endl;
			}
			myfile3<<endl;
			myfile3<<(end_of_str_traps - emi)<<" defects recovered"<<endl;
		}

		if (t > 0.999 && hit4 == 0) {
			hit4 = 1;
			myfile3 <<"t = 1e0"<<endl;
			for (itrec = bulk_defects.begin(); itrec != bulk_defects.end(); itrec++) {
				myfile3<<(itrec->second).x<<" "<<(itrec->second).y<<" "<<(itrec->second).z<<endl;
			}
			myfile3<<endl;
			myfile3<<(end_of_str_traps - emi)<<" defects recovered"<<endl;
		}

		if (t > 999 && hit5 == 0) {
			hit5 = 1;
			myfile3 <<"t = 1e3"<<endl;
			for (itrec = bulk_defects.begin(); itrec != bulk_defects.end(); itrec++) {
				myfile3<<(itrec->second).x<<" "<<(itrec->second).y<<" "<<(itrec->second).z<<endl;
			}
			myfile3<<endl;
			myfile3<<(end_of_str_traps - emi)<<" defects recovered"<<endl;
		}
    }

    // Uncomment for single file
    //myfile1<<sw_time_1[1]<<' '<<emi<<endl;


    //myfile1<<emi<<endl;
    myfile1<<delVot<<endl;

    init_sites.clear();
    defect_rates.clear();
    bulk_defects.clear();

}
    myfile1.close();
    myfile3.close();
    myfile2.close();
    myfile.close();

	return 0;
}

double randnu(){
    return (rand() * 1.0) / (RAND_MAX);
}

int rand_l_lim(int l_lim){
	int randllim = (int) floor((randnu()*(l_lim)));
	if(randllim==l_lim) return l_lim-1;
	return randllim;
}

int rand_w_lim(int w_lim){
	int randwlim = (int) floor((randnu()*(w_lim)));
	if(randwlim==w_lim) return w_lim-1;
	return randwlim;
}

int rand_tox_lim(int z_lim){
	int randtoxlim = (int) floor((randnu()*(z_lim)));
	if(randtoxlim==z_lim) return z_lim-1;
	return randtoxlim;
}

double roundin(double d) {
	return floor(d + 0.5);
}

void setzero(double &val) {
	if (val > -1e-5 && val < 1e-5) {
		val = 0;
	}
}

void carry_out_reaction() {
	//K31
	if (r2 < temp_1) {
		temp11 = 0;
		multimap <int, threeD>::iterator it;
		for (it = bulk_defects.begin(); it != bulk_defects.end(); it++){

			temp12 = temp11/alpha0;
			temp11 = temp11 + defect_rates[it->second].second;
			temp13 = temp11/alpha0;

			if ((r2 >= temp12) && (r2 < temp13)) {
                init_sites.push_back(it->second);
                alpha_k13_total = alpha_k13_total + defect_rates[it->second].first;
                alpha_k31_total = alpha_k31_total - defect_rates[it->second].second;

                delVot -= 1000*q*Tox*(1 - ((it->second).z)/z_lim)/(dev_area*ESiO2*E0);
                bulk_defects.erase(it);
                //setzero(alpha_k31_total);
                emi--;
                break;
			}
		}
	}

	//K13
	else if (r2 < temp_2) {
		temp11 = temp_1*alpha0;
        vector <threeD>::iterator it;
		for (it = init_sites.begin(); it != init_sites.end(); it++){

			temp12 = temp11/alpha0;
			temp11 = temp11 + defect_rates[*it].first;
			temp13 = temp11/alpha0;
			if ((r2 >= temp12) && (r2 < temp13)) {

				bulk_defects.insert(pair<int, threeD> (0,*it));
				alpha_k13_total = alpha_k13_total - defect_rates[*it].first;
				alpha_k31_total = alpha_k31_total + defect_rates[*it].second;

				delVot += 1000*q*Tox*(1 - (it->z)/z_lim)/(dev_area*ESiO2*E0);
				init_sites.erase(it);

                //setzero(alpha_k13_total);
				emi++;
				break;
			}
		}
	}
}

#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
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
double compartment_length;
double z_comp_length;
int num_compartment, domains_per_nm;
int j;
int nr1, nr2, nf1, nf2, ndz1, ndz2;
double q1, q2, q3;
int c1, c2, c3, c4;
int k1, k2, k3;

double alpha0, yemp, tau;
double temp13, temp11, temp12;
double width, length, height;
int w_lim, l_lim, z_lim;

#define WIDTH  50
#define LENGTH 20
#define HEIGHT 1

#define fpsi(VGS,VFB,psi) VGS-VFB-psi-((ESi*E0*sqrt(2)*(VT/Ldp)*sqrt((exp(-psi/VT)+(psi/VT)-1)+((exp(-2*phif/VT))*(exp(psi/VT)-(psi/VT)-1))))/cox)

double alpha_k13_total;
double alpha_k31_total;


double temp_1, temp_2;
double err;

double psi_str, FSiO2;
int emi, emiH2, previous_emi, previous_emiH2;
double str_Nit[10000];
double str_time[10000];
double t,t_prev;
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
    double Vgs_dc = 1.5; //DC stess voltage in Volts
    double Vgr_dc = 0; //DC recovery voltage in Volts
	double Tox = HEIGHT*1e-7;
	double ND = 5e16; //  Channel Doping in atoms/cm^3
    double N0 = 6.5e12; //  Channel Doping in atoms/cm^3
    double vfb = 0.4992; // Flatband voltage in Volts
    double E1 = 0.65;
    double E2 = 0;
    double EV = 0;
	double SHW = 5.9;
	double R = 0.52;
    double gamma = 3.2e-7;
/*********************************************************************************/

	double VT1 = (Kb*300)/q;// Thermal voltage Vt@ 300 K = 0.0259eV in eV
	double VT= (Kb*(273+Temp))/q;//Thermal voltage Vt@ 273+Temp in eV
	double ni = ni0*(pow(((273+Temp)/300),1.5))*exp(-Eg/2*(1/VT-1/VT1)); // in cm^-3
	double phif=VT*log(ND/ni);                // Bulk potential in V
	double beta = 1/VT;                       // in /V

	double kt=Kb*(273+Temp);//// Joules
	double cox=(ESiO2*E0)/Tox; // in F/cm^2
	double Ldp=sqrt((ESi*E0*kt)/(q*q*ND)); ////// Debye Length in cm

    double pi=22.0/7;
	double mp=0.5*9.1e-31;                   // hole mass in kg (http://ecee.colorado.edu/~bart/book/effmass.htm)
    double vp=sqrt(8*kt/(pi*mp))*100;         // should be thermal velocity (check the formula) in cm/s
    double sp=(1e-15) ;                       // cross-section area in cm^2
    double h=6.62607004e-34;                  // in Joule*seconds
    double NV=2*(pow((2*pi*mp*Kb*(273+Temp)/(pow(h,2))),1.5))*1e-6; // in /cm^3
    double p0=(pow(ni,2))/ND;                 // Bulk minority concentration  /cm^3

void ox_field(double VGb_str)
{
	//Solving for surface potential
	int cnt=1;

	int maxiter = 2500;
	double lpsi=0;
	double rpsi = 3*phif;

	double fvall = fpsi(VGb_str,vfb,lpsi);

	double fvalu = fpsi(VGb_str,vfb,rpsi);

    double newpsi;
	double tmp=0.0;
	double tol=1e-10;
	double fvalleft,fvalright;

	if (fvall*fvalu < 0.0)
	{
		newpsi = (lpsi+rpsi)/2;
		do {
			fvalleft  = fpsi(VGb_str,vfb,lpsi);
			fvalright = fpsi(VGb_str,vfb,newpsi);

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

	//Surface Electric Field
	double FSi=sqrt(2)*(VT/Ldp)*sqrt((exp(-psi_str/VT)+(psi_str/VT)-1)+((exp(-2*phif/VT))*(exp(psi_str/VT)-(psi_str/VT)-1)));

	//Electric Field in SiO2
	FSiO2=FSi*ESi/ESiO2;
}

void rate(double xvecbytox)
{
    double ps=p0*(exp(beta*psi_str)); // Surface minority carrier concentration near the interface which is inversion carrier density in  /cm^3
    double pf_ps = ps*vp*sp;
    double pf_nv = NV*vp*sp;

    double delET = E1-EV;
    double epsT1_str = E2- delET - gamma*(1-xvecbytox)*FSiO2;

    double eps12d_str = SHW/(pow((1 + R),2)) + R*epsT1_str/(1+R) + R*(pow(epsT1_str,2))/(4*SHW);
    double eps21d_str = SHW/(pow((1 + R),2)) - epsT1_str/(1+R) + R*(pow(epsT1_str,2))/(4*SHW);

    k13 = pf_ps*exp(-beta*eps12d_str);
		if(k13 > pf_ps)
			k13 = pf_ps;

    k31 = pf_nv*exp(-beta*eps21d_str);
		if(k31 > pf_nv)
			k31 = pf_nv;

}

int main(){

	srand (99989);
	double sw_time_1 [] = {1e-6, 1e3};

	num_defects = (int) LENGTH*WIDTH*HEIGHT*1/10;

	compartment_length = 0.5e-9;
	z_comp_length = 0.5e-9;




//  Rate constants for Stress
//    k13 = 9;
//	k31 = 1;

	mul_fac = 1/(compartment_length*compartment_length*z_comp_length*1e6);

	float modder = (1e13)/(num_defects);

	width = WIDTH*1e-9;
	w_lim = (int) ceil(width/compartment_length);
	length = LENGTH*1e-9;
	l_lim = (int) ceil(length/compartment_length);
	height = HEIGHT*1e-9;
	z_lim =  (int) ceil(height/z_comp_length);


    ox_field(Vgs_dc);

    for(int o = 0; o <= z_lim-1; o++)
    {
        rate(o/z_lim);
        k13_arr.push_back(k13);
        k31_arr.push_back(k31);
    }

	bool hit1 = 0;
	bool hit2 = 0;
	bool hit3 = 0;
	bool hit4 = 0;
	bool hit5 = 0;

	int d1 =0;
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
            alpha_k13_total += k13_arr[c3];
			init_sites.push_back(threeD(c1,c2,c3));
		}
	}



//*********************************************************************************
//******************************     STRESS     ***********************************
//*********************************************************************************



	emi = 0;

	t = sw_time_1[0];
	t_prev = t;
	float i[] = {2,2,2,2,2,2,2,2,2,2,2};
	counter = 0;
	ofstream myfile1,myfile,myfile2,myfile3;
	myfile1.open ("Rec.txt", ios::out | ios::trunc);
	myfile.open ("Str.txt", ios::out | ios::trunc);
    myfile2.open ("Str_data.txt", ios::out | ios::trunc);
    myfile3.open ("Rec_data.txt", ios::out | ios::trunc);
    multimap <int, threeD>::iterator itstr,itrec;

	kil=1;
	int flagger =0;

	while (t < sw_time_1[1] && flagger == 0){

		kil =0;
		r1 = randnu();
		while (r1 ==0){
			r1 = randnu();
		}

		alpha0 =  alpha_k13_total+ alpha_k31_total;

		tau = (1/alpha0)*log(1/r1);
		if (tau < 0){
			flagger = 1;
			myfile << "Error"<<endl;
		}

		r2 = randnu();
		while (r2 ==1){
			r2 = randnu();
		}

		temp_1 = alpha_k31_total/alpha0;
		temp_2 = temp_1 + (alpha_k13_total/alpha0);

		carry_out_reaction();

        if (t == sw_time_1[0])
        {
            t = t + tau;
            curr = (int) floor(log10(t));
            prev = curr;
        }

        else
        {
            t = t + tau;
            curr = (int) ceil(log10(t));
        }

        if(curr!=prev)
        {
            prev = curr;
            myfile<<t<<' '<<emi<<endl;
            t_prev = t + 0.1*pow(10,curr);
        }

        else
        {
            if(t > t_prev)
            {
                myfile<<t<<' '<<emi<<endl;
                t_prev = t + 0.1*pow(10,curr);
            }
        }

//        if (( t < 1e-6*i[0] ) && ( t >= 1e-6*(i[0]-1) ) && (i[0]<11)){
//			myfile <<t<<' '<< emi <<endl;
//
//			i[0]=i[0]+1;
//		}
//
//		else if(( t < 1e-5*i[1] ) && ( t >= 1e-5*(i[1]-1) ) && (i[1]<11)){
//			myfile <<t<<' '<< emi <<endl;
//			i[1]=i[1]+1;
//		}
//
//		else if(( t < 1e-4*i[2] ) && ( t >= 1e-4*(i[2]-1) ) && (i[2]<11)){
//			myfile <<t<<' '<< emi <<endl;
//			i[2]=i[2]+1;
//		}
//
//		else if(( t < 1e-3*i[3] ) && ( t >= 1e-3*(i[3]-1) ) && (i[3]<11)){
//			myfile <<t<<' '<< emi <<endl;
//			i[3]=i[3]+1;
//		}
//
//		else if(( t < 1e-2*i[4] ) && ( t >= 1e-2*(i[4]-1) ) && (i[4]<11)){
//			myfile <<t<<' '<< emi <<endl;
//			i[4]=i[4]+1;
//		}
//
//		else if(( t < 1e-1*i[5] ) && ( t >= 1e-1*(i[5]-1) ) && (i[5]<11)){
//			myfile <<t<<' '<< emi <<endl;
//			i[5]=i[5]+1;
//		}
//
//		else if(( t < i[6] ) && ( t >= (i[6]-1) ) && (i[6]<11)){
//			myfile <<t<<' '<< emi <<endl;
//			i[6]=i[6]+1;
//		}
//
//		else if(( t < 1e1*i[7] ) && ( t >= 1e1*(i[7]-1) ) && (i[7]<11)){
//			myfile <<t<<' '<< emi <<endl;
//			i[7]=i[7]+1;
//		}
//
//		else if(( t < 1e2*i[8] ) && ( t >= 1e2*(i[8]-1) ) && (i[8]<11)){
//			myfile <<t<<' '<< emi <<endl;
//			i[8]=i[8]+1;
//		}
//
//		else if(( t < 1e3*i[9] ) && ( t >= 1e3*(i[9]-1) ) && (i[9]<11)){
//			myfile <<t<<' '<< emi <<endl;
//			i[9]=i[9]+1;
//		}
//
//		else if(( t < 1e4*i[10] ) && ( t >= 1e4*(i[10]-1) ) && (i[10]<11)){
//			myfile <<t<<' '<< emi <<endl;
//			i[10]=i[10]+1;
//		}

		counter = counter + 1;
	}

	//myfile <<t<<' '<<emi<<endl;

	end_of_str_traps = emi;

	myfile2 << emi <<" out of "<< num_defects<<" defects broken in bulk"<<endl;

	for (itstr = bulk_defects.begin(); itstr != bulk_defects.end(); itstr++) {
		myfile2<<(itstr->second).x<<" "<<(itstr->second).y<<" "<<(itstr->second).z<<endl;
	}
	myfile2<<endl;



//*********************************************************************************
//****************************     RECOVERY     ***********************************
//*********************************************************************************

	t = sw_time_1[0];
    t_prev = t;
//  fill(i,i+sizeof(i),2);
	for (int k=0; k<11; k++) {
    i[k] = 2;
    }
	counter = 0;

//  Rate constants for recovery
//	k13 = 1;
//	k31 = 6;

//  Total Propensities
    ox_field(Vgr_dc);

    for(int o = 0; o < z_lim; o++)
    {
        rate(o/z_lim);
        k13_arr[o] = k13;
        k31_arr[o] = k31;
    }
    alpha_k13_total = 0;
    alpha_k31_total = 0;

    for(map<int,threeD>::iterator o = bulk_defects.begin(); o != bulk_defects.end(); o++)
    {
        alpha_k31_total += k31_arr[o->second.z];
    }

    for(vector<threeD>::iterator o = init_sites.begin(); o != init_sites.end(); o++)
    {
        alpha_k13_total += k13_arr[o->z];
    }


    while (t < sw_time_1[1] && flagger == 0){

		kil = 0;
		r1 = randnu();
		while (r1 ==0){
			r1 = randnu();
		}

		alpha0 =  alpha_k13_total+ alpha_k31_total;

		tau = (1/alpha0)*log(1/r1);
		if (tau < 0){
			flagger = 1;
			myfile << "Error"<<endl;
		}

		r2 = randnu();
		while (r2 ==1){
			r2 = randnu();
		}

		temp_1 = alpha_k31_total/alpha0;
		temp_2 = temp_1 + (alpha_k13_total/alpha0);

		carry_out_reaction();

        if (t == sw_time_1[0])
        {
            t = t + tau;
            curr = (int) floor(log10(t));
            prev = curr;
        }

        else
        {
            t = t + tau;
            curr = (int) ceil(log10(t));
        }

        if(curr!=prev)
        {
            prev = curr;
            myfile1<<t<<' '<<emi<<endl;
            t_prev = t + 0.1*pow(10,curr);
        }

        else
        {
            if(t > t_prev)
            {
                myfile1<<t<<' '<<emi<<endl;
                t_prev = t + 0.1*pow(10,curr);
            }
        }


//		if (( t < 1e-6*i[0] ) && ( t >= 1e-6*(i[0]-1) ) && (i[0]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[0]=i[0]+1;
//		}
//
//		else if(( t < 1e-5*i[1] ) && ( t >= 1e-5*(i[1]-1) ) && (i[1]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[1]=i[1]+1;
//		}
//
//		else if(( t < 1e-4*i[2] ) && ( t >= 1e-4*(i[2]-1) ) && (i[2]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[2]=i[2]+1;
//		}
//
//		else if(( t < 1e-3*i[3] ) && ( t >= 1e-3*(i[3]-1) ) && (i[3]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[3]=i[3]+1;
//		}
//
//		else if(( t < 1e-2*i[4] ) && ( t >= 1e-2*(i[4]-1) ) && (i[4]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[4]=i[4]+1;
//		}
//
//		else if(( t < 1e-1*i[5] ) && ( t >= 1e-1*(i[5]-1) ) && (i[5]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[5]=i[5]+1;
//		}
//
//		else if(( t < i[6] ) && ( t >= (i[6]-1) ) && (i[6]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[6]=i[6]+1;
//		}
//
//		else if(( t < 1e1*i[7] ) && ( t >= 1e1*(i[7]-1) ) && (i[7]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[7]=i[7]+1;
//		}
//
//		else if(( t < 1e2*i[8] ) && ( t >= 1e2*(i[8]-1) ) && (i[8]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[8]=i[8]+1;
//		}
//
//		else if(( t < 1e3*i[9] ) && ( t >= 1e3*(i[9]-1) ) && (i[9]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[9]=i[9]+1;
//		}
//
//		else if(( t < 1e4*i[10] ) && ( t >= 1e4*(i[10]-1) ) && (i[10]<11)){
//			myfile1 <<t<<' '<< emi <<endl;
//			i[10]=i[10]+1;
//		}

//        myfile <<t<< ' '<< emi <<endl;

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

	return 0;
}

double randnu(){
    return (rand() * 1.0) / (RAND_MAX);
}

int rand_l_lim(int l_lim){
	int randllim = (int) roundin((randnu()*(l_lim-1)));
	return randllim;
}

int rand_w_lim(int w_lim){
	int randwlim = (int) roundin((randnu()*(w_lim-1)));
	return randwlim;
}

int rand_tox_lim(int z_lim){
	int randtoxlim = (int) roundin((randnu()*(z_lim-1)));
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
			temp11 = temp11 + k31;
			temp13 = temp11/alpha0;

			if ((r2 >= temp12) && (r2 < temp13)) {
                init_sites.push_back(it->second);
                alpha_k13_total = alpha_k13_total + k13_arr[it->second.z];
                alpha_k31_total = alpha_k31_total - k31_arr[it->second.z];
                bulk_defects.erase(it);
                setzero(alpha_k31_total);
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
			temp11 = temp11 + k13;
			temp13 = temp11/alpha0;
			if ((r2 >= temp12) && (r2 < temp13)) {

				bulk_defects.insert(pair<int, threeD> (0,*it));
				alpha_k13_total = alpha_k13_total - k13_arr[it->z];
				alpha_k31_total = alpha_k31_total + k31_arr[it->z];
				init_sites.erase(it);

                setzero(alpha_k13_total);
				emi++;
				break;
			}
		}
	}
}

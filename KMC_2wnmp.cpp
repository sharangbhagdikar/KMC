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
double ini_trap1_check;
double ini_trap2_check;
double width, length, height;
int w_lim, l_lim, z_lim;

#define WIDTH  50
#define LENGTH 20
#define HEIGHT 1

double alpha_k13_total;
double alpha_k31_total;


double temp_1, temp_2;
double err;

double Y1, Y2;
int emi, emiH2, previous_emi, previous_emiH2;
double str_Nit[10000];
double str_time[10000];
double t;
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
multimap<int, threeD> bulk_defects;

int main(){

	srand (99989);
	double sw_time_1 [] = {1e-6, 1e3};

	num_defects = (int) LENGTH*WIDTH*HEIGHT*1/10;

	compartment_length = 0.5e-9;
	z_comp_length = 0.5e-9;

//  Rate constants for Stress
    k13 = 9e10;
	k31 = 1e10;

	mul_fac = 1/(compartment_length*compartment_length*z_comp_length*1e6);

	float modder = (1e13)/(num_defects);

	width = WIDTH*1e-9;
	w_lim = (int) ceil(width/compartment_length);
	length = LENGTH*1e-9;
	l_lim = (int) ceil(length/compartment_length);
	height = HEIGHT*1e-9;
	z_lim =  (int) ceil(height/z_comp_length);

	bool hit1 = 0;
	bool hit2 = 0;
	bool hit3 = 0;
	bool hit4 = 0;
	bool hit5 = 0;

	int d1 =0;

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
		}
	}

    alpha_k13_total = num_defects*k13;    // Initial total forward propensity: from ground state to transport state
    alpha_k31_total = 0;                  // Initial reverse propensity:




//*********************************************************************************
//******************************     STRESS     ***********************************
//*********************************************************************************



	emi = 0;

	t = sw_time_1[0];
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

		t = t + tau;
        myfile <<t<<' '<<emi<<endl;
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

//  fill(i,i+sizeof(i),2);
	for (int k=0; k<11; k++) {
    i[k] = 2;
    }
	counter = 0;

//  Rate constats for recovery
	k13 = 1e10;
	k31 = 6e10;

//  Total Propensities
    alpha_k31_total = bulk_defects.size()*k31;
    alpha_k13_total = init_sites.size()*k13;


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

		t = t + tau;

		if (( t < 1e-6*i[0] ) && ( t >= 1e-6*(i[0]-1) ) && (i[0]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[0]=i[0]+1;
		}

		else if(( t < 1e-5*i[1] ) && ( t >= 1e-5*(i[1]-1) ) && (i[1]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[1]=i[1]+1;
		}

		else if(( t < 1e-4*i[2] ) && ( t >= 1e-4*(i[2]-1) ) && (i[2]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[2]=i[2]+1;
		}

		else if(( t < 1e-3*i[3] ) && ( t >= 1e-3*(i[3]-1) ) && (i[3]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[3]=i[3]+1;
		}

		else if(( t < 1e-2*i[4] ) && ( t >= 1e-2*(i[4]-1) ) && (i[4]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[4]=i[4]+1;
		}

		else if(( t < 1e-1*i[5] ) && ( t >= 1e-1*(i[5]-1) ) && (i[5]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[5]=i[5]+1;
		}

		else if(( t < i[6] ) && ( t >= (i[6]-1) ) && (i[6]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[6]=i[6]+1;
		}

		else if(( t < 1e1*i[7] ) && ( t >= 1e1*(i[7]-1) ) && (i[7]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[7]=i[7]+1;
		}

		else if(( t < 1e2*i[8] ) && ( t >= 1e2*(i[8]-1) ) && (i[8]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[8]=i[8]+1;
		}

		else if(( t < 1e3*i[9] ) && ( t >= 1e3*(i[9]-1) ) && (i[9]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[9]=i[9]+1;
		}

		else if(( t < 1e4*i[10] ) && ( t >= 1e4*(i[10]-1) ) && (i[10]<11)){
			myfile1 <<t<<' '<< emi <<endl;
			i[10]=i[10]+1;
		}

        //myfile <<t<< ' '<< emi <<endl;

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
                bulk_defects.erase(it);
                alpha_k13_total = alpha_k13_total + k13;
                alpha_k31_total = alpha_k31_total - k31;
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
				init_sites.erase(it);

                alpha_k13_total = alpha_k13_total - k13;
				alpha_k31_total = alpha_k31_total + k31;
				setzero(alpha_k13_total);
				emi++;
				break;
			}
		}
	}
}

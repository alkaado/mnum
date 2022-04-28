#include <inttypes.h>
#include "../lib/cnumlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <portaudio.h>
#ifdef _DEBUG_
#endif

#define _BUF_LEN_ 4096


void sigsin(vector_t* x, double ampl , freq_t freq, double dt){
	for(int i = 0 ; i < vlen(x) ; i++){
		complex tmp = vit(x, i);
		vit(x, i) = 
			cadd(tmp, 
			cplx(ampl * sin(2 * M_PI* freq * i * dt), 0));
	}
}


void out_data(vector_t* vx, vector_t* vy){
	for(int i = 0 ; i < vlen(vx) ; i++){
		printf("%f %f %f\n",
		creal(vit(vx, i)),
		creal(vit(vy, i)),
		cimag(vit(vy, i)));
	} 
}

void out_data1(vector_t* vx){
	for(int i = 0 ; i < vlen(vx) ; i++){
		printf("%f\n",
		creal(vit(vx, i)));
	} 
}

void out_data2(vector_t* vx, vector_t* vy, vector_t* svx, vector_t* svy){
	for(int i = 0 ; i < vlen(vx) ; i++){
		printf("%f %f %f\n",
		creal(vit(vx, i)),
		creal(vit(vy, i))
		);
	} 

	for(int i = 0 ; i < vlen(svx) ; i++){
		printf("  %f %f \n",
		creal(vit(svx, i)),
		creal(vit(svy, i))
		);
	} 
}


void test_filter(){
	double dt = 1.0 / 60000;
	vector(org, _BUF_LEN_);
	vector(transf, _BUF_LEN_);
	vector(squaref, _BUF_LEN_);
	vector(freq_vect, _BUF_LEN_);
	vector(out, _BUF_LEN_);
	vector(ref, _BUF_LEN_);
	vector(triang, _BUF_LEN_);
	linspace(vx, 0, (dt * _BUF_LEN_), _BUF_LEN_);

	sigsin(org, 200.0, 19000, dt);
	sigsin(org, 100.0, 15000, dt);
	sigsin(org, 60.0, 10000, dt);
	sigsin(org, 60.0, 4000, dt);
	sigsin(org, 30.0, 2000, dt);
	sigsin(org, 100.0, 800, dt);
	//FFT_square(squaref, dt, 10000 , 500);
	//FFT_square(squaref, dt, 10000 , 1000);
	FFT_triang(triang, dt, 10000 , 1000);
	FFT_triang(triang, dt, 15000 , 1000);
	//FFT_triang(triang, dt, 19000 , 1000);
	FFT_freq(freq_vect, dt);
	FFT(transf, org);

	vmul(transf, triang);
	IFT(transf, out);
		
	//out_data(freq_vect, squaref);	
	//out_data(freq_vect, varg(transf));	
	vrotl(out, _BUF_LEN_ / 4);
	out_data(vx, out);	
	//out_data(freq_vect, vmuli(vabs(vdivi(transf, _BUF_LEN_)), 2.0));	
	//out_data(freq_vect, vabs(transf));	
	//out_data(freq_vect, vmul(transf, squaref));	
	//out_data(vx, org);	
	//out_data(freq_vect, triang);	
}

void testShift(){
	linspace(vx, 1, 12 , 12);
	linspace(vx1, 13, 24 , 12);
	out_data1(vx);

	vshla(vx, vx1, 3);	
	vrotl(vx1,3);
	printf("----- %d \n", vx->_act_el_);
	out_data1(vx);

	vshla(vx, vx1, 3);	
	vrotl(vx1,3);
	printf("----- %d \n", vx->_act_el_);
	out_data1(vx);
}


void testMatrix(){
	matrix(A, 4, 4);	
	minit(A, 
		0, 1, 0, 2,
		1, 2, 3, 0, 
		-1, 2, 4, 1,
		2, 0, 1, 1,
		);

	printf("--------\n");
	printf("A:\n");
	mdump(A);	

	printf("--------\n");
	printf("A - 1:\n");
	msubi(A, cplx(1, 0));
	mdump(A);	

	printf("--------\n");
	printf("A + 1:\n");
	maddi(A, cplx(1, 0));
	mdump(A);	

	matrix(B, 3, 3);
	minit(B,
		3, 1, 1,
		4, 2, -1,
		-2, -1, 1);

	matrix(C, 4, 4);
	minit(C,
		2, 1, 1, 1,
		-1, 2, 3, 2,
		-1, 0, 0, 0,
		1, 2, 3, 1
		);
	printf("--------\n");
	printf("Det:\n");
	matrix_t* Alg = mcinv(C);
	mdump(Alg);
	printf("\n");
	printf("--------\n");
}


void testEquation(){
	matrix(A, 3, 3);
	minit(A, 
		2, 3, 2,
		1, -2, -3,
		3, 1, 1
		);

	mdump(A);
	printf("---------\n");

	matrix(Y, 1, 3);
	minit(Y,
		3,
		-3,
		4
		);

	mdump(Y);
	printf("---------\n");
	matrix_t* X = msolve(A, Y);
	mdump(X);
	mfree(X);
	printf("---------\n");
}

void testApprox(){
	vector(svx, 7);
	vinit(svx, 
		0,
		0.2,
		0.4,
		0.6,
		0.8,
		1.0,
		1.2
		);

	vector(svy, 7);
	vinit(svy,
		0.08,
		0.2,
		0.25,
		0.27,
		0.22,
		0.15,
		0.01
		);

	linspace(vx, 0, 1.5, 100);
	zeros(vy, 100);
	approx_poly(vx, vy, 5, svx, svy);
	out_data(vx, vy);

}

void testApproxTryg(){
	vector(svx, 7);
	vinit(svx, 
		0,
		0.2,
		0.4,
		0.6,
		0.8,
		1.0,
		1.2,
		1.4
		);

	vector(svy, 7);
	vinit(svy,
		0.08,
		0.2,
		0.25,
		0.27,
		0.22,
		0.15,
		0.1
		);

	linspace(vx, 0, 5, 100);
	zeros(vy, 100);
	approx_tryg(vx, vy, 2, svx, svy);
	out_data(vx, vy);
}

void testInterpPoly(){
	vector(svx, 7);
	vinit(svx, 
		0,
		0.2,
		0.4,
		0.6,
		0.8,
		1.0,
		1.2,
		1.4
		);

	vector(svy, 7);
	vinit(svy,
		0.08,
		0.2,
		0.25,
		0.27,
		0.22,
		0.15,
		0.1
		);

	linspace(vx, 0, 1, 100);
	zeros(vy, 100);
	interp_poly(vx, vy, svx, svy);
	out_data(vx, vy);
}

void testInterpSpline1(){
	vector(svx, 7);
	vinit(svx, 
		0,
		0.2,
		0.4,
		0.6,
		0.8,
		1.0,
		1.2,
		1.4
		);

	vector(svy, 7);
	vinit(svy,
		0.08,
		0.2,
		0.25,
		0.27,
		0.22,
		0.15,
		0.1
		);

	linspace(vx, 0, 1, 100);
	zeros(vy, 100);
	interp_spline1(vx, vy, svx, svy);
	out_data(vx, vy);
}


void testInterpSpline3(){
	double f(double x){
		return x + cos(2 * x);
	};

	linspace(svx, 0, 5, 20);
	zeros(svy, 20);
	for(int i = 0 ; i < vlen(svx) ; i++){
		vit(svy, i) = cplx(f(creal(vit(svx, i))), 0);
	}

	linspace(vx, 0, 5, 500);
	zeros(vy, 500);
	interp_poly(vx, vy, svx, svy);
	//approx_poly(vx, vy, 6, svx, svy);
	out_data(vx, vy);
}

void testIntegr(){

	void fx1(vector_t* vx, vector_t* vy){
		vcopy(vy, vx);	
	}

	linspace(vx, 0, 2, 50);
	zeros(vy, 50);

	fx1(vx, vy);

	complex S = integr(vx, vy, 0, 2, INTEGR_SQ);
	cprintf(S);
	printf("\n");
	S = integr(vx, vy, 0, 2, INTEGR_TR);
	cprintf(S);
	printf("\n");
	S = integr(vx, vy, 0, 2, INTEGR_PA);
	cprintf(S);
	printf("\n");
}

void testMresize(){
	matrix_t* M = mmalloc(2, 2);
	minit(M, 1, 2,
		 3, 4);

	M = mresize(M, 3, 2);
	mdump(M);
	mfree(M);
}

void testFroots(){
	linspace(vx, 0, 7, 8);
	zeros(vy, 8);
	vit(vy, 1) = cplx(1, 0);
	vit(vy, 3) = cplx(-1, 0);
	vit(vy, 4) = cplx(1, 0);
	vit(vy, 6) = cplx(-1, 0);

	double f(double x){
		vector(tvx, 1);
		vector(tvy, 1);
		vit(tvx, 0) = cplx(x, 0);
		interp_poly(tvx, tvy, vx, vy);
		return creal(vit(tvy, 0));
	};


	vector_t* vr = froots(f, 0, 7, 0.001);
	vdump(vr);
	linspace(ix, 0, 7, 100);
	zeros(iy, 100);
	//interp_poly(ix, iy, vx, vy);
	//out_data(ix, iy);
}

int main(){
	//testMatrix();
	//testEquation();
	//testApprox();
	//testApproxTryg();
	//testInterpPoly();
	//testInterpSpline1();
	//testShift();
	//test_filter();
	//testIntegr();
	//testInterpSpline3();
	//testMresize();
	testFroots();
	/*PaError err = Pa_Initialize();

	if( err != paNoError ){

		Pa_Terminate();
	}*/

	return 0;
};

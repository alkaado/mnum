#ifndef _CNUMLIB_H_
#define _CNUMLIB_H_

#include <math.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#define creal(z) z.real
#define cimag(z) z.imag
#define cabs(z) sqrt( pow( creal(z), 2 ) + pow( cimag(z), 2 ))
#define carg(z) atan(cimag(z) / creal(z))
#define cplx(re, im) (complex){re, im}
#define sgn(arg) ( arg < 0.0 ? -1.0 : 1.0 )
#define cexp(base, arg) cplx(base * cos( sgn(arg) * arg ), base * sgn(arg) * sin( sgn(arg) * arg))
#define cprintf(c) if( cimag(c) > 0 ){printf("%.4f + %.4fj", creal(c), cimag(c));}else{printf("%.4f", creal(c));}


#define vindex(v, i) 	( (v->_act_el_ + i) < v->_length_ ? (v->_act_el_ + i) : ( v->_act_el_ + i - v->_length_ ) )
#define vit(v, i) 	(v->_data_[vindex(v, (i))])

#define vrotl(v, i)	v->_act_el_ = ( v->_act_el_ + (i)) % v->_length_
#define vrotr(v, i)	v->_act_el_ = ( v->_act_el_ - (i)) % v->_length_; \
			v->_act_el_ = ( v->_act_el_ < 0) ? v->_act_el_ + v->_length_ : v->_act_el_;

#define vector(v, size) \
	vector_t* v; \
	{ \
		static complex ___##v##tmp___[size]; \
		static vector_t ___##v##___ = (vector_t){size, 0, size - 1, 0, ___##v##tmp___}; \
		v = &___##v##___; \
	}


#define linspace(v, x1, x2, n) \
	vector_t* v; \
	{ \
		static complex ___##v##tmp___[n]; \
		static vector_t ___##v##___ = (vector_t){n, 0, (n) - 1, 0 , ___##v##tmp___}; \
		v = &___##v##___; \
		for( int i = 0 ; i < n ; i++){ \
			double ___##v##xx1___  = x1; \
			double ___##v##xx2___  = x2; \
			vit(v, i) = cplx(___##v##xx1___ + \
			(___##v##xx2___ - ___##v##xx1___ ) * i / (n - 1), \
			0.0); \
		}; \
	}

#define zeros(v, n) \
	linspace(v, 0.0, 0.0, n)

#define ones(v, n) \
	linspace(v, 1.0, 1.0, n)

#define NUMARGS(...) (sizeof( (double[]) {__VA_ARGS__}) / sizeof(double))

#define vinit(v, ...) \
	{ \
		double ___##v##vec___[] = {__VA_ARGS__}; \
		for(int i = 0 ; i < vlen(v) ; i++){ \
			vit(v, i) = cplx(___##v##vec___[i], 0.0); \
		} \
	}


#define minit(m, ...) { \
	double ___##m##vec___[] = {__VA_ARGS__}; \
	for(int i = 0 ; i < msize(m, 0) ; i++){ \
		for(int j = 0 ; j < msize(m,1) ; j++){ \
			mit(m, i, j) = cplx(___##m##vec___[j + i * msize(m, 1)], 0.0); \
		} \
	} \
}

#define mvec(m) (m->_mat_data_)

/*
* Retruns the item of matrix by given row and col number
* Indexes starts from 0.
* m - (matrix_t*) matrix
* r - (uint32_t) row index
* c - (uint32_t) col index
*/
#define mit(m, r, c)	( vit(m->_mat_data_, (c) + mindex(m,(r))))
#define mindex(m, r)	( m->_col_len_ * (r) )

/*
* Creates a matrix with given size
* m - matrix variable name
* rlen - (uint32_t) row length 
* clen - (uint32_t) col length
*/
#define matrix(m, rlen, clen) \
	matrix_t* m; \
	{ \
		vector(___vec##m##___, (rlen) * (clen)); \
		static matrix_t ___mat##m##___; \
		___mat##m##___._row_len_ = clen; \
		___mat##m##___._col_len_ = rlen; \
		___mat##m##___._mat_data_ = ___vec##m##___; \
		m = &___mat##m##___; \
	}

/*
* Returns size of the matrix
* dim - determines the dimension for which the size should be
* returned: 0 - ydim (number of rows), 1 - xdim (number of
* columns)
*/
#define msize(m, dim)	( (dim) == 0 ? m->_row_len_ : m->_col_len_)

typedef double freq_t;

typedef struct {
	double real;
	double imag;
} complex;

typedef struct{
	int32_t _length_;
	int32_t _act_el_;
	uint32_t _last_el_;
	uint32_t _first_el_;
	complex* _data_;	
} vector_t;

typedef struct{
	int32_t _row_len_;
	int32_t _col_len_;
	vector_t* _mat_data_;
} matrix_t;

/******************************************************************* 
* Complex numbers
********************************************************************/
complex cmul(complex z1, complex z2);			// z3 = z1 * z2
complex cdiv(complex z1, complex z2);			// z3 = z1 / z2
complex cadd(complex z1, complex z2);			// z3 = z1 + z2
complex csub(complex z1, complex z2);			// z3 = z1 - z2
complex cpow(complex z1, float p);			// z3 = z1^p

/******************************************************************* 
* Vectors
********************************************************************/
uint32_t vlen( vector_t* v );				//length of v
vector_t* vmalloc(uint32_t len);			
void vfree(vector_t* v);
vector_t* vresize(vector_t* a, int32_t newSize);
vector_t* vcopy( vector_t* a, vector_t* b);		// a = b
vector_t* vshla( vector_t* a, vector_t* b, int32_t s); 	// last a = b << s
vector_t* vmul( vector_t* a, vector_t* b);		// a = a * b
vector_t* vadd( vector_t* a, vector_t* b);		// a = a + b
vector_t* vsub( vector_t* a, vector_t* b);		// a = a - b
vector_t* vpow( vector_t* a, float p);			// a = a^p
vector_t* vabs( vector_t* a);				// a = |a|
vector_t* vdivi( vector_t* a, double vi);		// a = a/vi
vector_t* vmuli( vector_t* a, double vi);		// a = a * vi
vector_t* varg( vector_t* a);				// arg(a)
vector_t* vreal( vector_t* a);				// re(a)
vector_t* vimag( vector_t* a);				// im(a)
void vdump(vector_t* v);


/******************************************************************* 
* Matrixes
********************************************************************/
matrix_t* mmalloc(uint32_t r, uint32_t c);
matrix_t* mresize(matrix_t* m, uint32_t r, uint32_t c);
void mfree(matrix_t* m);
void mcopy(matrix_t* m1, matrix_t* m2);			// m1 = m2
void mdump(matrix_t* m);				// show(m)
matrix_t* madd(matrix_t* m1, matrix_t* m2);		// m1 = m1 + m2
matrix_t* mcadd(matrix_t* m1, matrix_t* m2);		// m3 = m1 + m2
matrix_t* msub(matrix_t* m1, matrix_t* m2);		// m1 = m1 - m2
matrix_t* mcsub(matrix_t* m1, matrix_t* m2);		// m3 = m1 - m2
matrix_t* mcmul(matrix_t* m1, matrix_t* m2);		// m3 = m1 * m2
matrix_t* mctr(matrix_t* m1);				// m2 = m1'
matrix_t* mcdiv(matrix_t* m1, matrix_t* m2);		// m3 = m1 \ m2
matrix_t* mrmcol(matrix_t* m1, uint32_t c);		// m1 = m1[ -c ]
matrix_t* mrmrow(matrix_t* m1, uint32_t r);		// m1 = m1[ -r ]
complex mdet(matrix_t* m1);				// det(m1)
matrix_t* mcalgc(matrix_t* m1);				// m1^D
matrix_t* mcinv(matrix_t* m1);				// m2 = m1^-1
matrix_t* maddi(matrix_t* m1, complex a);		// m1 = m1 + a
matrix_t* msubi(matrix_t* m1, complex a);		// m1 = m1 - a
matrix_t* mmuli(matrix_t* m1, complex a);		// m1 = m1 * a
matrix_t* mdivi(matrix_t* m1, complex a);		// m1 = m1 / a
matrix_t* mcaddi(matrix_t* m1, complex a);		// m2 = m1 + a
matrix_t* mcsubi(matrix_t* m1, complex a);		// m2 = m1 - a
matrix_t* mcmuli(matrix_t* m1, complex a);		// m2 = m1 * a
matrix_t* mcdivi(matrix_t* m1, complex a);		// m2 = m1 / a
vector_t* mdiag(matrix_t* m);				// diag(m)

/******************************************************************* 
* Numerical operation
********************************************************************/
void interp_poly(vector_t* vx, vector_t* vy, vector_t* svx, vector_t* svy);
void interp_spline1(vector_t* vx, vector_t* vy, vector_t* svx, vector_t* svy);
void interp_spline3(vector_t* vx, vector_t* vy, vector_t* svx, vector_t* svy);

matrix_t* approx_poly_ind(vector_t* vx, vector_t* vy, uint16_t n);
void approx_poly(vector_t* vx, vector_t* vy, uint16_t n, vector_t* svx, vector_t* svy);
matrix_t* approx_tryg_ind(vector_t* vx, vector_t* vy, uint16_t N);
void approx_tryg(vector_t* vx, vector_t* vy, uint16_t N, vector_t* svx, vector_t* svy);

#define INTERP_SPLINE1 	0
#define INTERP_SPLINE3	1
#define INTERP_POLY	2
typedef double (*froots_f) (double);
vector_t* froots(froots_f f, double x0, double x1, double prec);

#define INTEGR_SQ	0
#define INTEGR_TR	1
#define INTEGR_PA	2
complex integr(vector_t* vx, vector_t* vy, double x0, double xn, int method);

/******************************************************************* 
* Fourier Transform
********************************************************************/
void DFT(vector_t* X, vector_t* x);
void IDFT(vector_t* X, vector_t* x);
void FFT(vector_t* X, vector_t* x);
void IFT(vector_t* X, vector_t* x);
void FFT_freq(vector_t* vf, double dt);
void FFT_square(vector_t* sq, float dt, freq_t base, freq_t width);
void FFT_triang(vector_t* sq, float dt, freq_t base, freq_t width);


/******************************************************************* 
* IMPLEMENTATION
********************************************************************/
complex cmul(complex z1, complex z2){
	complex result = {creal(z1) * creal(z2) - cimag(z1) * cimag(z2) , creal(z1) * cimag(z2) + cimag(z1) * creal(z2)};
	return result;
}

complex cdiv(complex z1, complex z2){
	creal(z2) = creal(z2) != 0 ? 1.0/creal(z2) : 0;
	cimag(z2) = cimag(z2) != 0 ? 1.0/cimag(z2) : 0;

	return cmul(z1, z2);
}

complex cadd(complex z1, complex z2){
	complex result = {creal(z1) + creal(z2) , cimag(z1) + cimag(z2)};
	return result;
}

complex csub(complex z1, complex z2){
	complex result = {creal(z1) - creal(z2) , cimag(z1) - cimag(z2)};
	return result;
}

complex cpow(complex z1, float p){
	double a = cabs(z1);
	double ar = (creal(z1) == 0.0 && cimag(z1) == 0.0 ) ? 1.0 : carg(z1);

	complex result = cexp(pow(a, p), p * ar);
	return result;
}

vector_t* vmalloc(uint32_t len){
	vector_t* tmp = malloc(sizeof(vector_t));
	tmp->_data_ = malloc(len * sizeof(complex));
	tmp->_length_ = len;
	tmp->_act_el_ = 0;
	tmp->_last_el_ = len - 1;
	tmp->_first_el_ = 0;
	return tmp;
}

vector_t* vresize(vector_t* a, int32_t newSize){
	vector_t* tmp = vmalloc( newSize );

	for(int i = 0 ; i < ( vlen(a) < vlen(tmp) ? vlen(a) : vlen(tmp)  ) ; i++){
		vit(tmp, i) = vit(a, i);
	}
	vfree(a);

	return tmp;
}

void vfree(vector_t* v){
	free(v->_data_);
	free(v);
}

//shift left vector b into vector a
//a is shifted left and b is added at the end
vector_t* vshla( vector_t* a, vector_t* b, int32_t s){
	s = s % vlen(a);

	for(int i = 0 ; i < s ; i++){
		vit(a, i) = vit(b, i);	
	}

	vrotl(a, s);

	return a;
}

uint32_t vlen( vector_t* v ){
	if( v == NULL )
		return 0;

	return v->_length_;
}

vector_t* vcopy( vector_t* a, vector_t* b){
	if(vlen(a) > vlen(b))
		return NULL;
	memcpy(a->_data_, b->_data_, vlen(a) * sizeof(complex));
	return a;
}

vector_t* vmul( vector_t* a, vector_t* b){
	if( vlen(a) != vlen(b) )
		return NULL;

	for( int i = 0 ; i < vlen(a) ; i++){
		vit(a, i) = cmul(vit(a,i), vit(b,i));	
	}

	return a;
}

vector_t* vadd( vector_t* a, vector_t* b){
	if( vlen(a) != vlen(b) )
		return NULL;

	for( int i = 0 ; i < vlen(a) ; i++){
		vit(a, i) = cadd(vit(a,i), vit(b,i));	
	}

	return a;
}

vector_t* vsub( vector_t* a, vector_t* b){
	if( vlen(a) != vlen(b) )
		return NULL;

	for( int i = 0 ; i < vlen(a) ; i++){
		vit(a, i) = csub(vit(a,i), vit(b,i));	
	}

	return a;
}

vector_t* vpow( vector_t* a, float p){
	for( int i = 0 ; i < vlen(a) ; i++){
		vit(a, i) = cpow(vit(a, i), p);	
	}

	return a;
}

vector_t* vabs( vector_t* a){
	for( int i = 0 ; i < vlen(a) ; i++){
		vit(a, i) = cplx(cabs(vit(a, i)), 0.0);	
	}

	return a;
}

vector_t* vdivi( vector_t* a, double vi){
	for( int i = 0 ; i < vlen(a) ; i++){
		vit(a, i) = cplx(creal(vit(a, i)) / vi, 
				cimag(vit(a, i)) / vi);	
	}

	return a;
}

vector_t* vmuli( vector_t* a, double vi){
	for( int i = 0 ; i < vlen(a) ; i++){
		vit(a, i) = cplx(creal(vit(a, i)) * vi, 
				cimag(vit(a, i)) * vi);	
	}

	return a;
}

vector_t* varg( vector_t* a){
	for( int i = 0 ; i < vlen(a) ; i++){
		vit(a, i) = cplx(carg(vit(a, i)), 0.0);	
	}

	return a;
}

vector_t* vreal( vector_t* a){
	for( int i = 0 ; i < vlen(a) ; i++){
		vit(a, i) = cplx(creal(vit(a, i)), 0.0);	
	}

	return a;
}

vector_t* vimag( vector_t* a){
	for( int i = 0 ; i < vlen(a) ; i++){
		vit(a, i) = cplx(cimag(vit(a, i)), 0.0);	
	}

	return a;
}

void vdump(vector_t* v){
	printf("[");
	for(int i = 0 ; i < vlen(v); i++){
		cprintf(vit(v, i));
		if( i < vlen(v) - 1 )
			printf(", ");
	}
	printf("]\n");
}

matrix_t* mmalloc(uint32_t r, uint32_t c){
	matrix_t* tmp = malloc(sizeof(matrix_t));
	tmp->_mat_data_ = vmalloc(r * c);
	tmp->_row_len_ = c;
	tmp->_col_len_ = r;
	return tmp;
}

matrix_t* mresize(matrix_t* m, uint32_t r, uint32_t c){
	matrix_t* tmp = mmalloc(r, c);

	for(int i = 0 ; i < ( (msize(m, 0) < msize(tmp, 0) ) ? msize(m, 0) : msize(tmp, 0)) ; i++){
		for(int j = 0 ; j < ( (msize(m, 1) < msize(tmp, 1) ) ? msize(m, 1) : msize(tmp, 1)) ; j++){
			mit(tmp, i, j) = mit(m, i, j);

		}
	}

	mfree(m);
	return tmp;
}

void mfree(matrix_t* m){
	vfree(m->_mat_data_);
	free(m);
}

void mcopy(matrix_t* m1, matrix_t* m2){
	m1->_row_len_ = m2->_row_len_;
	m1->_col_len_ = m2->_row_len_;
	vcopy(m1->_mat_data_, m2->_mat_data_);
}

void mdump(matrix_t* m){
	printf("[");
	for(int i = 0 ; i < msize(m, 0); i++){
		for( int j = 0 ; j < msize(m, 1) ; j++){
			if( j == 0 && i > 0)
				printf(" ");
				
			cprintf(mit(m, i, j));
			if( j < msize(m, 1) - 1)
				printf(", ");
		}
		if( i < msize(m, 0) - 1)
			printf("\n");
	}
	printf("]\n");
}

matrix_t* madd(matrix_t* m1, matrix_t* m2){
	if( msize(m1, 0) != msize(m2, 0) || msize(m1, 1) != msize(m2, 1))
		return NULL;
	vadd(mvec(m1), mvec(m2));
	return m1;
}

matrix_t* mcadd(matrix_t* m1, matrix_t* m2){
	if( msize(m1, 0) != msize(m2, 0) || msize(m1, 1) != msize(m2, 1))
		return NULL;

	matrix_t* m3 = mmalloc(msize(m1,0), msize(m1, 1));;
	mcopy(m3, m1);
	vadd(mvec(m3), mvec(m2));
	return m3;
}

matrix_t* msub(matrix_t* m1, matrix_t* m2){
	if( msize(m1, 0) != msize(m2, 0) || msize(m1, 1) != msize(m2, 1))
		return NULL;
	vsub(mvec(m1), mvec(m2));
	return m1;
}

matrix_t* mcsub(matrix_t* m1, matrix_t* m2){
	if( msize(m1, 0) != msize(m2, 0) || msize(m1, 1) != msize(m2, 1))
		return NULL;
	matrix_t* m3 = mmalloc(msize(m1,0), msize(m1, 1));;
	mcopy(m3, m1);
	vsub(mvec(m3), mvec(m2));
	return m3;
}

matrix_t* mcmul(matrix_t* m1, matrix_t* m2){
	
	complex mmulr1(matrix_t* A, matrix_t* B, uint32_t n, uint32_t offsetA, uint32_t offsetB){

		if( msize(A, 1) <= n || msize(B, 0) <= n )
			return cplx(0 , 0);

		complex An = mit(A, offsetA, n);
		complex Bn = mit(B, n, offsetB);

		An = cmul(An, Bn);
		An = cadd(An, mmulr1(A, B, n + 1, offsetA, offsetB));
		return An;
	}

	void mmulr2(matrix_t* A, matrix_t* B, matrix_t* C, uint32_t r, uint32_t c){
		if( msize(A, 0) <= r || msize(B, 1) <= c)
			return;

		complex An = mmulr1(A, B, 0, r, c);
		mmulr2(A, B, C, r, c + 1);
		mit(C, r, c) = An;
	}

	if( msize(m1, 1) != msize(m2, 0)){
		printf("Number of columns in A is different that rows in B");
		return NULL;
	}

	matrix_t* m3 = mmalloc(msize(m2, 1), msize(m1, 0));
	for( int i = 0 ; i < msize(m1, 0) ; i++){
		mmulr2(m1, m2, m3, i, 0);
	}

	return m3;
}

matrix_t* mctr(matrix_t* m1){
	matrix_t* tmp = mmalloc(msize(m1, 0), msize(m1,1));	
	for(int i = 0 ; i < msize(m1, 0) ; i++){
		for(int j = 0 ; j < msize(m1, 1) ; j++){
			mit(tmp, j, i) = mit(m1, i, j);
		}
	}

	return tmp;
}

matrix_t* msolve(matrix_t* A, matrix_t* Y){
	matrix_t* Ainv = mcinv(A);
	matrix_t* ans = mcmul(Ainv, Y);
	mfree(Ainv);
	return ans;
}

matrix_t* mrmcol(matrix_t* m1, uint32_t c){
	void mrmcol1(vector_t* V, uint32_t c, uint32_t n){

		if( n == vlen(V))
			return;

		if( n >= c )
			vit(V, n) = vit(V, n + 1);

		mrmcol1(V, c, n + 1);
	}

	if( c >= msize(m1, 1) )
		return NULL;

	for( int i = 0 ; i < msize(m1, 0) ; i++){
		mrmcol1(m1->_mat_data_, c + i * (msize(m1, 1) - 1), 0);
	}

	m1->_col_len_ = m1->_col_len_ - 1 ;
}

matrix_t* mrmrow(matrix_t* m1, uint32_t r){
	void mrmrow1(vector_t* V, uint32_t r, uint32_t n){

		if( n == vlen(V))
			return;

		if( n >= r )
			vit(V, n) = vit(V, n + 1);

		mrmrow1(V, r, n + 1);
	}

	if( r >= msize(m1, 0) )
		return NULL;

	mrmrow1(m1->_mat_data_, r * msize(m1, 1), 0);
	for( int i = 0 ; i < msize(m1, 1) - 1 ; i++){
		mrmrow1(m1->_mat_data_, r * (msize(m1, 1)), 0);
	}

	m1->_row_len_ = m1->_row_len_ - 1 ;
}


complex mdet(matrix_t* m1){
	complex mdet1(matrix_t* M, uint32_t r, uint32_t c){
		if( msize(M, 0) == 2 ){
			complex a00 = mit(M, 0, 0);
			complex a01 = mit(M, 0, 1);
			complex a10 = mit(M, 1, 0);
			complex a11 = mit(M, 1, 1);
			return csub(cmul(a00, a11), cmul(a10, a01));
		}

		complex ANS = cplx(0, 0);
		matrix_t* N = mmalloc(msize(M,0), msize(M, 1));
		for( int i = 0 ; i < msize(M, 0) ; i++){
			mcopy(N, M);
			complex anm = mit(M, r, i);
			mrmcol(N, i);
			mrmrow(N, r);
			complex ans = mdet1(N, r, i);
			ans = cmul(ans, cplx(pow(-1, r + i + 2), 0));
			ans = cmul(ans, anm);
			ANS = cadd(ANS, ans);
		}
		mfree(N);
		return ANS;
	}

	if(msize(m1, 1) == 1){
		return mit(m1, 0, 0);
	}else if(msize(m1, 0) == msize(m1, 1)){
		return mdet1(m1, 0, 0);
	} 	
	return cplx(0, 0);
}

matrix_t* mcalgc(matrix_t* m1){
	if(msize(m1, 0) == msize(m1, 1)){
		matrix_t* A = mmalloc(msize(m1,0), msize(m1, 1));
		matrix_t* B = mmalloc(msize(m1,0), msize(m1, 1));
		for(int i = 0 ; i < msize(m1, 0) ; i++){
			for(int j = 0 ; j < msize(m1, 1); j++){
				mcopy(B, m1);
				mrmrow(B, i);
				mrmcol(B, j);
				complex aij = cplx(pow(-1, j + i + 2), 0);
				complex M = mdet(B);
				mit(A, i, j) = cmul(aij, M);
			}
		}
		mfree(B);
		return A;
	}
	return NULL;
}

matrix_t* mcinv(matrix_t* m1){
	if(msize(m1, 0) == msize(m1, 1)){
		matrix_t* AD = mcalgc(m1);
		matrix_t* ADT = mctr(AD);
		complex detA = mdet(m1);
		matrix_t* A = mmuli(ADT, cdiv(cplx(1, 0), detA));
		mfree(AD);
		return A;
	}

	return NULL;
}

matrix_t* maddi(matrix_t* m1, complex a){
	for( int i = 0 ; i < msize(m1, 0) ; i++){
		for( int j = 0 ; j < msize(m1, 1) ; j++){
			mit(m1, i, j) = cadd(mit(m1, i, j), a);
		}
	}
	return m1;
}

matrix_t* msubi(matrix_t* m1, complex a){
	for( int i = 0 ; i < msize(m1, 0) ; i++){
		for( int j = 0 ; j < msize(m1, 1) ; j++){
			mit(m1, i, j) = csub(mit(m1, i, j), a);
		}
	}
	return m1;
}

matrix_t* mmuli(matrix_t* m1, complex a){
	for( int i = 0 ; i < msize(m1, 0) ; i++){
		for( int j = 0 ; j < msize(m1, 1) ; j++){
			mit(m1, i, j) = cmul(mit(m1, i, j), a);
		}
	}
	return m1;
}

matrix_t* mdivi(matrix_t* m1, complex a){
	for( int i = 0 ; i < msize(m1, 0) ; i++){
		for( int j = 0 ; j < msize(m1, 1) ; j++){
			mit(m1, i, j) = cdiv(mit(m1, i, j), a);
		}
	}
	return m1;
}

matrix_t* mcaddi(matrix_t* m1, complex a){
	matrix_t* ans = mmalloc(msize(m1, 0), msize(m1, 1));
	for( int i = 0 ; i < msize(m1, 0) ; i++){
		for( int j = 0 ; j < msize(m1, 1) ; j++){
			mit(ans, i, j) = cadd(mit(m1, i, j), a);
		}
	}
	return ans;
}

matrix_t* mcsubi(matrix_t* m1, complex a){
	matrix_t* ans = mmalloc(msize(m1, 0), msize(m1, 1));
	for( int i = 0 ; i < msize(m1, 0) ; i++){
		for( int j = 0 ; j < msize(m1, 1) ; j++){
			mit(ans, i, j) = csub(mit(m1, i, j), a);
		}
	}
	return ans;
}

matrix_t* mcmuli(matrix_t* m1, complex a){
	matrix_t* ans = mmalloc(msize(m1, 0), msize(m1, 1));
	for( int i = 0 ; i < msize(m1, 0) ; i++){
		for( int j = 0 ; j < msize(m1, 1) ; j++){
			mit(ans, i, j) = cmul(mit(m1, i, j), a);
		}
	}
	return ans;
}

matrix_t* mcdivi(matrix_t* m1, complex a){
	matrix_t* ans = mmalloc(msize(m1, 0), msize(m1, 1));
	for( int i = 0 ; i < msize(m1, 0) ; i++){
		for( int j = 0 ; j < msize(m1, 1) ; j++){
			mit(ans, i, j) = cdiv(mit(m1, i, j), a);
		}
	}
	return ans;
}

vector_t* mdiag(matrix_t* m){
	vector_t* diag = vmalloc(msize(m, 0));

	for(int i = 0 ; i < msize(m, 0) ; i++){
		vit(diag, i) = mit(m, i, i);
	}

	return diag;
}


void interp_poly(vector_t* vx, vector_t* vy, vector_t* svx, vector_t* svy){

	complex interp_poly1(complex x, vector_t* va, vector_t* vx){
		complex ans = vit(va, 0);

		for(int i = 1 ; i < vlen(va) ; i++){
			complex tmp = vit(va, i);

			for( int j = 0 ; j < i ; j++){
				tmp = cmul(tmp, 
						csub(x,
						     vit(vx, j)
						     ));
			}
			ans = cadd(ans, tmp);
		}
		return ans;
	}

	matrix_t* v = mmalloc(vlen(svx), vlen(svx));
	
	for(int i = 0 ; i < vlen(svy); i++){
		mit(v, i, 0) = vit(svy, i);
	}

	for(int m = 1 ; m < vlen(svx) ; m++){
		for(int k = m ; k < vlen(svx) ; k++){
			complex dy = mit(v, k, m-1);
			dy = csub(dy, mit(v, k-1, m-1));

			complex dx = vit(svx, k);
			dx = csub(dx, vit(svx, k - m));

			mit(v, k, m) = cdiv(dy, dx);
		}
	}

	vector_t* va = mdiag(v);
	mfree(v);

	for(int i = 0 ; i <vlen(vx) ; i++){
		vit(vy, i) = interp_poly1(vit(vx, i), va, svx); 
	}
	
	vfree(va);
}

void interp_spline1(vector_t* vx, vector_t* vy, vector_t* svx, vector_t* svy){

	complex phi1(complex x, vector_t* svx,  vector_t* svy, int i0, int i1){
		if( ( i1 - i0 ) == 1 ){
			matrix(A, 2, 2);
			minit(A,
				1, creal(vit(svx, i0)),
				1, creal(vit(svx, i1))
			);

			matrix(Y, 1, 2);
			minit(Y,
				creal(vit(svy, i0)),
				creal(vit(svy, i1))
			);

			matrix_t* X = msolve(A, Y);
			complex a = mit(X, 0, 0);
			complex b = mit(X, 1, 0);
			complex ans = cmul(x, b);
			ans = cadd(ans, a);
			mfree(X);
			return ans;
		}else{
			int i00 = i0;
			int i01 = i0 + ( (i1 - i0) % 2 ? (i1 - i0 - 1) : (i1 - i0) ) / 2;
			int i10 = i01;
			int i11 = i1;

			if ( creal(x) >= creal(vit(svx, i00)) && creal(x) <= creal(vit(svx, i01)) )
				return phi1(x, svx, svy, i00, i01);
			else if ( creal(x) >= creal(vit(svx, i10)) && creal(x) <= creal(vit(svx, i11)) )
				return phi1(x, svx, svy, i10, i11);
			else
				return cplx(0, 0);
		}
	}

	for(int j = 0 ; j < vlen(vx); j++){
		vit(vy, j) = phi1(vit(vx, j), svx, svy, 0, vlen(svx) - 1);
	}
}

void interp_spline3(vector_t* vx, vector_t* vy, vector_t* svx, vector_t* svy){
	uint32_t n = vlen(svx) - 1;
	double a = creal(vit(svx, 0));
	double b = creal(vit(svx, n));
	double h = (b - a) / n;

	double phi(double x, int32_t i){
		double xi(int32_t i){ return a + i * h; }
		
		if( x >= xi(i-2) && x <= xi(i-1)){
			return (1.0/pow(h, 3)) * pow(x - xi(i-2), 3);
		}else if( x >= xi(i-1)  && x <= xi(i)){
			return (1.0/pow(h, 3)) * (pow(x - xi(i-2), 3) - 4*pow(x - xi(i-1), 3));
		}else if( x >= xi(i)  && x <= xi(i+1)){
			return (1.0/pow(h, 3)) * (pow(xi(i+2) - x, 3) - 4*pow(xi(i+1) - x, 3));
		}else if( x >= xi(i+1)  && x <= xi(i+2)){
			return (1.0/pow(h, 3)) * pow(xi(i+2) - x, 3);
		}else{
			return 0;
		}
	}

	double alpha(){
		double x0 = creal(vit(svx, 0));
		double x1 = creal(vit(svx, 1));
		double y0 = creal(vit(svy, 0));
		double y1 = creal(vit(svy, 1));

		return (y1 - y0)/(x1 - x0);
	}

	double beta(){
		double x0 = creal(vit(svx, n - 1));
		double x1 = creal(vit(svx, n));
		double y0 = creal(vit(svy, n - 1));
		double y1 = creal(vit(svy, n));

		return (y1 - y0)/(x1 - x0);
	}

	double S3(double x, vector_t* vc){
		double S = phi(x, -1) * (creal(vit(vc, 1)) - (h/3)*alpha());
		
		for(int i = 0 ; i < n + 1; i++){
			S += phi(x, i) * creal(vit(vc, i));
		}
		
		S += phi(x, n + 1) * (creal(vit(vc, n - 1)) + (h/3)*beta());
		return S;
	}

	
	//Fullfil M
	matrix_t* M = mmalloc(n + 1, n + 1);
	mit(M, 0, 0) = cplx(4, 0);
	mit(M, 0, 1) = cplx(2, 0);

	for( int i = 1 ; i < n ; i++){
		mit(M, i, i-1) = cplx(1, 0);
		mit(M, i, i  ) = cplx(4, 0);
		mit(M, i, i+1) = cplx(1, 0);
	}
	
	mit(M, n, n-1 ) = cplx(2, 0);
	mit(M, n, n   ) = cplx(4, 0);

	//Fullfil Y
	matrix_t* Y = mmalloc(1,n+1);
	mit(Y, 0, 0) = cadd(vit(svy,0), cplx(h/3 * alpha(), 0));

	for( int i = 1 ; i < n ; i++){
		mit(Y, i, 0) = vit(svy, i);	
	}

	mit(Y, n, 0) = csub(vit(svy, n), cplx(h/3 * beta(), 0));
	matrix_t* C = msolve(M, Y);
	
	for( int i = 0 ; i < vlen(vx) ; i++ ){
		vit(vy, i) = cplx(S3(creal(vit(vx, i)), mvec(C)), 0);	
	}
	
	mfree(M);
	mfree(Y);
	mfree(C);
}

matrix_t* approx_poly_ind(vector_t* vx, vector_t* vy, uint16_t n){
	matrix_t* M = mmalloc(n, vlen(vx));
	matrix_t* Y = mmalloc(1, vlen(vy));

	for(int i = 0 ; i < msize(M, 0) ; i++){
		for( int m = 0 ; m < msize(M, 1); m++){
			mit(M, i, m) = cpow(vit(vx, i), m);
		}
	}

	vcopy(mvec(Y), vy);
	matrix_t* MT = mctr(M);
	matrix_t* MTM = mcmul(MT, M);
	matrix_t* MTY = mcmul(MT, Y);
	matrix_t* MX = msolve(MTM, MTY);

	mfree(M);
	mfree(MT);
	mfree(MTM);
	mfree(MTY);
	return MX;
}

void approx_poly(vector_t* vx, vector_t* vy, uint16_t n, vector_t* svx, vector_t* svy){
	matrix_t* ma = approx_poly_ind(svx, svy, n);
	vector_t* va = mvec(ma);

	for(int i = 0 ; i < vlen(vx) ; i++){
		complex S = cplx(0,0);

		for(int m = 0 ; m < vlen(va) ; m++){
			complex ai = vit(va, m);
			complex xi = vit(vx, i);
			xi = cpow(xi, m);
			xi = cmul(xi, ai);
			S = cadd(S, xi);		
		}
		vit(vy, i) = S;
	}

	mfree(ma);
}


matrix_t* approx_tryg_ind(vector_t* vx, vector_t* vy, uint16_t N){
	matrix_t* M = mmalloc(1 + (N - 1) * 2, vlen(vx));
	matrix_t* Y = mmalloc(1, vlen(vy));
	uint32_t n = vlen(vx) - 1;
	double a = creal(vit(vx, 0));
	double b = creal(vit(vx, n));
	double h = (b - a) / n;
	double l = ((n + 1.0) / 2.0) * h;
	double c = M_PI / l;

	for(int i = 0 ; i < msize(M, 0) ; i++){
		int cs = 0;
		int j = 1;
		for( int m = 0 ; m < msize(M, 1); m++){
			if( m == 0){
				mit(M, i, m) = cplx(1, 0);
			}else if(cs){
				mit(M, i, m) =
				cplx(cos(creal(vit(vx, i)) *
				j * c), 0);
				j++;
			}else{
				mit(M, i, m) =
				cplx(sin(creal(vit(vx, i)) *
				j * c), 0);
			}
			cs = cs ? 0 : 1;
		}
	}


	vcopy(mvec(Y), vy);
	matrix_t* MT = mctr(M);
	matrix_t* MTM = mcmul(MT, M);
	matrix_t* MTY = mcmul(MT, Y);
	matrix_t* MX = msolve(MTM, MTY);

	mfree(M);
	mfree(MT);
	mfree(MTM);
	mfree(MTY);
	return MX;
}

void approx_tryg(vector_t* vx, vector_t* vy, uint16_t N, vector_t* svx, vector_t* svy){
	matrix_t* ma = approx_tryg_ind(svx, svy, N);
	vector_t* va = mvec(ma);

	uint32_t n = vlen(svx) - 1;
	double a = creal(vit(svx, 0));
	double b = creal(vit(svx, n));
	double h = (b - a) / n;
	double l = ((n + 1.0) / 2.0) * h;
	double c = M_PI / l;


	for(int i = 0 ; i < vlen(vx) ; i++){
		complex S = vit(va, 0);
		complex xi = vit(vx, i);
		int j = 1;

		for(int m = 1 ; m < vlen(va) ; m++){
			double ci = cos(j * c * creal(xi));
			double si = sin(j * c * creal(xi));
			ci = ci * creal(vit(va, m));
			si = si * creal(vit(va, m + 1));
			S = cadd(S, cplx(ci, 0) );		
			S = cadd(S, cplx(si, 0) );		
			m++;
			j++;
		}
		vit(vy, i) = S;
	}

	mfree(ma);
}

void FFT_freq(vector_t* vf, double dt){
	uint32_t n = vlen(vf);
	float ndt = n * dt;
	uint32_t N = n/2;
	float dt2 = 1.0 / (2.0 * dt);
	
	for(int32_t i = 0 ; i < N  ; i++){
		creal(vit(vf, i)) = i / ndt;
		creal(vit(vf, (i + N))) = creal(vit(vf, i)) - dt2;
	}
}

void DFT(vector_t* X, vector_t* x){
	uint32_t n = vlen(X);
	complex S = cplx(0,0);

	for(uint32_t k = 0 ; k < n ; k++){
		for(uint32_t i = 0 ; i < n ; i++){
			S = cadd(S, cexp(creal(vit(x, i)), -2 * M_PI/n * k * i));
		}
		vit(X, k) = S;
		S = cplx(0,0);
	}
}

void IDFT(vector_t* X, vector_t* x){
	uint32_t n = vlen(X);
	complex S = cplx(0,0);

	for(uint32_t k = 0 ; k < n ; k++){
		for(uint32_t i = 0 ; i < n ; i++){
			S = cadd(S, cmul(vit(X, i), cexp(1,  2 * M_PI/n * k * i)));
		}
		vit(x, k) = cplx(creal(S)/n, 0);
		S = cplx(0,0);
	}
}

void FFT(vector_t* X, vector_t* x){
	uint32_t n = vlen(X);
	uint32_t N = n / 2;
	complex WN = cexp(1.0, -2.0 * M_PI/n);
	complex WNK = cplx(1,0);
	complex WNK2 = cplx(1, 0);
	complex WNK2I = cplx(1, 0);
	complex S1 = cplx(0,0);
	complex S2 = cplx(0,0);

	for(uint32_t k = 0 ; k < n ; k++){
		for(uint32_t i = 0 ; i < N ; i++){
			S1 = cadd(S1, cmul(vit(x, (i << 1)), WNK2I ));
			S2 = cadd(S2, cmul(vit(x, ((i << 1) + 1)), WNK2I ));
			WNK2I = cmul(WNK2I, WNK2);
		}
		vit(X, k) = cadd(S1, cmul(S2, WNK));
		WNK = cmul(WNK, WN);
		WNK2 = cmul(WNK, WNK);
		S1 = cplx(0,0);
		S2 = cplx(0,0);
	}
}

void IFT(vector_t* X, vector_t* x){
	uint32_t n = vlen(X);
	uint32_t N = n / 2;
	complex WN = cexp(1.0, 2.0 * M_PI/n);
	complex WNK = cplx(1, 0);
	complex WNK2 = cplx(1, 0);
	complex WNK2I = cplx(1, 0);
	complex S1 = cplx(0,0);
	complex S2 = cplx(0,0);

	for(uint32_t k = 0 ; k < n ; k++){
		for(uint32_t i = 0 ; i < N ; i++){
			S1 = cadd(S1, cmul(vit(X, (i << 1)), WNK2I ));
			S2 = cadd(S2, cmul(vit(X, (i << 1) + 1), WNK2I ));
			WNK2I = cmul(WNK2I, WNK2);
		}
		vit(x, k) = cplx(creal(cadd(S1, cmul(S2, WNK)))/n, 0);
		WNK = cmul(WNK, WN);
		WNK2 = cmul(WNK, WNK);
		S1 = cplx(0,0);
		S2 = cplx(0,0);
	}
}

void FFT_square(vector_t* sq, float dt, freq_t base, freq_t width){
	uint32_t n = vlen(sq);	
	float ndt = n * dt;
	uint32_t N = n/2;
	float dt2 = 1.0 / (2.0 * dt);
	
	for(int32_t i = 0 ; i < N  ; i++){
		float tf = i / ndt;
		if (tf >= ( base - width/2 ) && tf <= (base + width / 2)) {
			vit(sq, i) = cplx(1.0 , 0.0);
			vit(sq, n - i) = cplx(1.0 , 0.0);
		}
	}
}

void FFT_triang(vector_t* sq, float dt, freq_t base, freq_t width){
	uint32_t n = vlen(sq);	
	float ndt = n * dt;
	uint32_t N = n/2;
	float dt2 = 1.0 / (2.0 * dt);
	
	for(int32_t i = 0 ; i < N  ; i++){
		float tf = i / ndt;
		
		if (tf >= ( base - width/2 ) && tf <= base) {
			double x = tf - (base - width/2);
			float a = 1.0 / (width/2);
			float y = x * a;
			vit(sq, i) = cplx(y , 0.0);
			vit(sq, (n - i)) = cplx(y , 0.0);
		} else if( tf >= base && tf <= (base + width/2)){
			double x = tf - base;
			float a = 1.0 / (width/2);
			float y = 1.0 - x * a;
			vit(sq, i) = cplx(y , 0.0);
			vit(sq, (n - i)) = cplx(y , 0.0);
		}
	}
}

complex integr(vector_t* vx, vector_t* vy, double x0, double xn, int method){
	complex integr_sq(vector_t* vx, vector_t* vy, double x0, double xn){
		complex S = cplx(0, 0);

		for(int i = 0 ; i < vlen(vx) - 1 ; i++){
			double xi = cabs(vit(vx, i));
			double xi1 = cabs(vit(vx, i + 1));
			double dx = xi1 - xi;
			complex yi = vit(vy, i + 1);

			if( xi >= x0 && xi1 <= xn){
				complex tmp = cmul(yi, cplx(dx, 0));
				S = cadd(S, tmp);
			}
		}
		return S;
	}

	complex integr_tr(vector_t* vx, vector_t* vy, double x0, double xn){
		complex S = cplx(0, 0);

		for(int i = 0 ; i < vlen(vx) - 1 ; i++){
			double xi = cabs(vit(vx, i));
			double xi1 = cabs(vit(vx, i + 1));
			double dx = xi1 - xi;
			complex yi = vit(vy, i);
			complex yi1 = vit(vy, i + 1);

			if( xi >= x0 && xi1 <= xn){
				complex tmp = cadd(yi, yi1);
				tmp = cmul(tmp, cplx(dx / 2 , 0));
				S = cadd(S, tmp);
			}
		}
		return S;
	}

	complex integr_pa(vector_t* vx, vector_t* vy, double X0, double Xn){
		complex S = cplx(0, 0);

		for(int i = 0 ; i < vlen(vx) - 3 ; i+=4){
			double x0 = cabs(vit(vx, i));
			double x1 = cabs(vit(vx, i + 1));
			double x2 = cabs(vit(vx, i + 2));
			double dx = x2 - x0;
			complex y0 = vit(vy, i);
			complex y1 = vit(vy, i + 1);
			complex y2 = vit(vy, i + 2);

			if( x0 >= X0 && x2 <= Xn){
				complex tmp = cmul(y1, cplx(4, 0));
				tmp = cadd(tmp, y0);
				tmp = cadd(tmp, y2);
				tmp = cmul(tmp, cplx(dx/3, 0));
				S = cadd(S, tmp);
			}
		}
		return S;
	}
	
	if(method == INTEGR_SQ){
		return integr_sq(vx, vy, x0, xn);	
	}else if(method == INTEGR_TR){
		return integr_tr(vx, vy, x0, xn);	
	}else if(method == INTEGR_PA){
		return integr_pa(vx, vy, x0, xn);	
	}
	
	return integr_sq(vx, vy, x0, xn);	

}


vector_t* froots(froots_f f, double x0, double x1, double prec){

	vector_t* froots1(double x0, double x1, vector_t* vr){
		double a = f(x0);
		double b = f(x1);
		double x00 = x0;
		double x01 = x0 + (x1 - x0) / 2;
		double x10 = x01;
		double x11 = x1;

		if( (a*b <= 0) && ( fabs(x0-x1) <= prec )){
			if(vr == NULL){
				vr = vmalloc(1);
			}else{
				vr = vresize(vr, vlen(vr) + 1);
			}
			vit(vr, vlen(vr) - 1) = cplx(x0, 0.0);

			return vr;

		}else if(  ( fabs(x1-x0) > prec ) ){
			vr = froots1(x00, x01, vr);
			vr = froots1(x10, x11, vr);
			return vr;
		}

		return vr;
	}


	vector_t* vr = NULL;
	vr = froots1(x0, x1, vr);
	return vr;
}


#endif

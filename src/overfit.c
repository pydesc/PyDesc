/*
 * Copyright 2011 Maciek Dziubinski 
 * Copyright 2003-2011, 2017 Pawel Daniluk
 * Copyright 2005 Bosco K Ho
 * 
 * This file is part of PyDesc.
 * 
 * PyDesc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PyDesc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with PyDesc.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include"overfit.h"

#include"simple_macros.h"

#ifndef SOLARIS
#ifndef PGI

#define _mix0000(v0,v1) __builtin_ia32_shufps(v0, v1, 0x00)  // 0 -> 0000 -> 0000 -> 0x00
#define _mix1111(v0,v1) __builtin_ia32_shufps(v0, v1, 0x55)  // 0 -> 1111 -> 1111 -> 0x55
#define _mix2222(v0,v1) __builtin_ia32_shufps(v0, v1, 0xAA)  // 0 -> 2222 -> 2222 -> 0xAA

#define _1st(v) _mix0000(v,v)
#define _2nd(v) _mix1111(v,v)
#define _3rd(v) _mix2222(v,v)



#define _zero ({float tmp[4]={.0,.0,.0,.0}; LOAD(tmp);})

#define _sum(v) ({float tmp[4]; __builtin_ia32_storeups(tmp, v); tmp[0]+tmp[1]+tmp[2];})


typedef v4sf t_coord;


#define LOAD(X) (__builtin_ia32_loadups(X))
#define STORE(X, Y) (__builtin_ia32_storeups((X), (Y)))

#define ADD3(X,Y,Z) ((X)+(Y)+(Z))
#define ADD(X,Y) ((X)+(Y))
#define SUB(X,Y) ((X)-(Y))
#define MUL(X,Y) ((X)*(Y))
#define DIV(X,Y) ((X)/(Y))


#define _allput(i) ({float _i=(i); float _tmp[4]={_i,_i,_i,_i}; LOAD(_tmp);})
#define ISZERO(X) ({float tmp[4]; STORE(tmp, (X)); (tmp[0]==0.0 && tmp[1]==0.0 && tmp[2]==0.0)?1:0;})

#else /* SOLARIS */

typedef struct {
	float x;
	float y;
	float z;
	float crap;  // it is easier when structure occupies 128 bits
} t_coord;

t_coord _zero={0, 0, 0, 0};


#define LOAD(X) ({ float *p=(X); t_coord res=*(t_coord *)(p); res;})
#define STORE(X, Y) ({ float *XX=(X); t_coord YY=(Y); XX[0]=YY.x; XX[1]=YY.y; XX[2]=YY.z; XX[3]=YY.crap;})

#define ADD3(X,Y,Z) ({ t_coord XX=(X); t_coord YY=(Y);  t_coord ZZ=(Z); t_coord res={XX.x+YY.x+ZZ.x, XX.y+YY.y+ZZ.y, XX.z+YY.z+ZZ.z, XX.crap+YY.crap+ZZ.crap}; res;})
#define ADD(X,Y) ({ t_coord XX=(X); t_coord YY=(Y); t_coord res={XX.x+YY.x, XX.y+YY.y, XX.z+YY.z, XX.crap+YY.crap}; res;})
#define SUB(X,Y) ({ t_coord XX=(X); t_coord YY=(Y); t_coord res={XX.x-YY.x, XX.y-YY.y, XX.z-YY.z, XX.crap-YY.crap}; res;})
#define MUL(X,Y) ({ t_coord XX=(X); t_coord YY=(Y); t_coord res={XX.x*YY.x, XX.y*YY.y, XX.z*YY.z, XX.crap*YY.crap}; res;})
#define DIV(X,Y) ({ t_coord XX=(X); t_coord YY=(Y); t_coord res={XX.x/YY.x, XX.y/YY.y, XX.z/YY.z, XX.crap/YY.crap}; res;})

#define _1st(v) ({ t_coord V=(v); t_coord res={V.x, V.x, V.x, V.x}; res;})
#define _2nd(v) ({ t_coord V=(v); t_coord res={V.y, V.y, V.y, V.y}; res;})
#define _3rd(v) ({ t_coord V=(v); t_coord res={V.z, V.z, V.z, V.z}; res;})
#define _sum(v) ({ t_coord V=(v); V.x+V.y+V.z;})

#define _allput(i) ({float _i=(i); t_coord res={_i,_i,_i,_i}; res;})
#define ISZERO(X) ({t_coord XX=(X); (XX.x==0.0 && XX.y==0.0 && XX.z==0.0)?1:0;})

#endif /* PGI */
#endif /* SOLARIS */

/*
t_coord vAB0;
t_coord vAB1;
t_coord vAB2;

t_coord vU0;
t_coord vU1;
t_coord vU2;
*/

#define MAX_THREADS 128

t_coord vstartA[MAX_THREADS];
t_coord vstartB[MAX_THREADS];

t_coord vSumA[MAX_THREADS];
t_coord vSumB[MAX_THREADS];
t_coord vSumAA[MAX_THREADS];
t_coord vSumBB[MAX_THREADS];
t_coord vSumAB0[MAX_THREADS];
t_coord vSumAB1[MAX_THREADS];
t_coord vSumAB2[MAX_THREADS];

int n_points[MAX_THREADS];

int tokens[MAX_THREADS];


void __attribute__ ((constructor)) overfit_init(void)
{
#pragma omp critical(overfit_token)
	{
//		printf("overfit_init\n");
		memset(tokens, 0, sizeof(tokens));
	}
}


int overfit_get_token()
{
	int res=-1;

#pragma omp critical(overfit_token)
	{
		for(int i=0; i<MAX_THREADS; i++) if(!tokens[i]) {
			tokens[i]=1;
			res=i;
			break;
		}

	}

	if(res<0) {
		printf("out of overfit tokens\n");
		abort();
	}

	return res;
}

void overfit_release_token(int token)
{
#pragma omp critical(overfit_token)
	{
		tokens[token]=0;
	}
}


void overfit_reset(int token)
{
	n_points[token]=0;
	vSumA[token]=_zero;
	vSumB[token]=_zero;
	vSumAA[token]=_zero;
	vSumBB[token]=_zero;
	vSumAB0[token]=_zero;
	vSumAB1[token]=_zero;
	vSumAB2[token]=_zero;
	vstartA[token]=_zero;
	vstartB[token]=_zero;
}

static void modify_sums(t_overfit_sums *sums, t_coord newStartA, t_coord newStartB);

void overfit_setstart(int token, float p1[4], float p2[4])
{
    vstartA[token]=LOAD(p1);
    vstartB[token]=LOAD(p2);

	if(n_points[token]!=0) {
		t_overfit_sums tmp_sums;
		overfit_sumsave_str(token,&tmp_sums);

		modify_sums(&tmp_sums,vstartA[token],vstartB[token]);

		vSumA[token]=LOAD(tmp_sums.A);
		vSumB[token]=LOAD(tmp_sums.B);
		vSumAA[token]=LOAD(tmp_sums.AA);
		vSumBB[token]=LOAD(tmp_sums.BB);
		vSumAB0[token]=LOAD(tmp_sums.AB0);
		vSumAB1[token]=LOAD(tmp_sums.AB1);
		vSumAB2[token]=LOAD(tmp_sums.AB2);
	}

}

void overfit_add(int token, float p1[4], float p2[4])
{
	n_points[token]++;

	t_coord cA=SUB(LOAD(p1),vstartA[token]);
	t_coord cB=SUB(LOAD(p2),vstartB[token]);
	
	vSumA[token]=ADD(vSumA[token], cA);
	vSumB[token]=ADD(vSumB[token], cB);
	vSumAA[token]=ADD(vSumAA[token], MUL(cA, cA));
	vSumBB[token]=ADD(vSumBB[token], MUL(cB, cB));

	vSumAB0[token]=ADD(vSumAB0[token], MUL(_1st(cA), cB));
	vSumAB1[token]=ADD(vSumAB1[token], MUL(_2nd(cA), cB));
	vSumAB2[token]=ADD(vSumAB2[token], MUL(_3rd(cA), cB));
}

float overfit_dist(float p1[4], float p2[4])
{
	t_coord vp1=LOAD(p1);
	t_coord vp2=LOAD(p2);

	t_coord vtmp=SUB(vp2, vp1);
	vtmp=MUL(vtmp, vtmp);

	return sqrt(_sum(vtmp));
}
	
float overfit_trans_dist_sqr(float p1[4], float p2[4], t_transrot TR)
{
	t_coord vp1=LOAD(p1);
	t_coord vp2=LOAD(p2);

	t_coord vU0=LOAD(&(TR.U[0][0]));
	t_coord vU1=LOAD(&(TR.U[1][0]));
	t_coord vU2=LOAD(&(TR.U[2][0]));
	
	t_coord vTr=LOAD(TR.Tr);
	
	t_coord vp3=ADD(ADD3(MUL(vU0, _1st(vp2)), MUL(vU1, _2nd(vp2)), MUL(vU2, _3rd(vp2))), vTr);

	t_coord vtmp=SUB(vp3, vp1);
	vtmp=MUL(vtmp, vtmp);

	return _sum(vtmp);
}
	
void overfit_trans(float p2[4], float p1[4], t_transrot TR)
{
	t_coord vp1=LOAD(p1);

	t_coord vU0=LOAD(&(TR.U[0][0]));
	t_coord vU1=LOAD(&(TR.U[1][0]));
	t_coord vU2=LOAD(&(TR.U[2][0]));
	
	t_coord vTr=LOAD(TR.Tr);
	
	t_coord vp2=ADD(ADD3(MUL(vU0, _1st(vp1)), MUL(vU1, _2nd(vp1)), MUL(vU2, _3rd(vp1))), vTr);

	STORE(p2, vp2);
}

float overfit_trans_dist(float p1[4], float p2[4], t_transrot TR)
{
	return sqrt(overfit_trans_dist_sqr(p1, p2, TR));
}

void overfit_trans_combine(t_transrot *TR, t_transrot TR1, t_transrot TR2)
{
	t_coord vU10=LOAD(&(TR1.U[0][0]));
	t_coord vU11=LOAD(&(TR1.U[1][0]));
	t_coord vU12=LOAD(&(TR1.U[2][0]));
	t_coord vTr1=LOAD(TR1.Tr);
	
	t_coord vU20=LOAD(&(TR2.U[0][0]));
	t_coord vU21=LOAD(&(TR2.U[1][0]));
	t_coord vU22=LOAD(&(TR2.U[2][0]));
	t_coord vTr2=LOAD(TR2.Tr);

	t_coord vU30=ADD3(MUL(_1st(vU10), vU20), MUL(_2nd(vU10), vU21), MUL(_3rd(vU10), vU22));
	t_coord vU31=ADD3(MUL(_1st(vU11), vU20), MUL(_2nd(vU11), vU21), MUL(_3rd(vU11), vU22));
	t_coord vU32=ADD3(MUL(_1st(vU12), vU20), MUL(_2nd(vU12), vU21), MUL(_3rd(vU12), vU22));

	t_coord vTr3=ADD(ADD3(MUL(vU20, _1st(vTr1)), MUL(vU21, _2nd(vTr1)), MUL(vU22, _3rd(vTr1))), vTr2);

	STORE(TR->Tr, vTr3);
	STORE(&(TR->U[0][0]), vU30);
	STORE(&(TR->U[1][0]), vU31);
	STORE(&(TR->U[2][0]), vU32);
}

static void new_starts(t_overfit_sums *sums1, t_overfit_sums *sums2, t_coord *newStartA_p, t_coord *newStartB_p)
{
	int N=sums1->N+sums2->N;

    if(N==0) {
        *newStartA_p=_zero;
        *newStartB_p=_zero;
    } else {
        t_coord actA1=ADD(LOAD(sums1->A), MUL(LOAD(sums1->startA), _allput(sums1->N)));
        t_coord actB1=ADD(LOAD(sums1->B), MUL(LOAD(sums1->startB), _allput(sums1->N)));
        t_coord actA2=ADD(LOAD(sums2->A), MUL(LOAD(sums2->startA), _allput(sums2->N)));
        t_coord actB2=ADD(LOAD(sums2->B), MUL(LOAD(sums2->startB), _allput(sums2->N)));

        *newStartA_p=DIV(ADD(actA1, actA2), _allput(N));
        *newStartB_p=DIV(ADD(actB1, actB2), _allput(N));
    }
}

void modify_sums(t_overfit_sums *sums, t_coord newStartA, t_coord newStartB)
{
	t_coord deltaA=SUB(newStartA, LOAD(sums->startA)); // deltaA = newStartA - oldStartB
	t_coord deltaB=SUB(newStartB, LOAD(sums->startB)); // analogously

    t_coord A=LOAD(sums->A);
    t_coord B=LOAD(sums->B);

    t_coord vN=_allput((float)sums->N);

    t_coord deltaAvN=MUL(deltaA, vN);
    t_coord deltaBvN=MUL(deltaB, vN);

	//saving changes:
	
	// updating starts:
	STORE(sums->startA, newStartA);
	STORE(sums->startB, newStartB);

	// AA' = AA - 2*A*deltaA + N*deltaA^2 = AA - deltaA*(2*A - N*deltaA)
	STORE(sums->AA, SUB(LOAD(sums->AA), MUL(deltaA, SUB(MUL(A, _allput(2.0)), deltaAvN))));

	// BB' = BB - 2*B*deltaB + N*deltaB^2 = BB - deltaB*(2*B - N*deltaB)
	STORE(sums->BB, SUB(LOAD(sums->BB), MUL(deltaB, SUB(MUL(B, _allput(2.0)), deltaBvN))));


	// AB0' = AB0 - deltaA_0*B - deltaB*A_0 + N*deltaA_0*deltaB = AB0 - (deltaB*A_0 + deltaA_0*(B - N*deltaB))
	STORE(sums->AB0, SUB(LOAD(sums->AB0), ADD(MUL(deltaB, _1st(A)), MUL(_1st(deltaA), SUB(B, deltaBvN)))));

	// AB1' = AB1 - deltaA_1*B - deltaB*A_1 + N*deltaA_1*deltaB = AB1 - (deltaB*A_1 + deltaA_1*(B - N*deltaB))
	STORE(sums->AB1, SUB(LOAD(sums->AB1), ADD(MUL(deltaB, _2nd(A)), MUL(_2nd(deltaA), SUB(B, deltaBvN)))));

	// AB2' = AB2 - deltaA_2*B - deltaB*A_2 + N*deltaA_2*deltaB = AB2 - (deltaB*A_2 + deltaA_2*(B - N*deltaB))
	STORE(sums->AB2, SUB(LOAD(sums->AB2), ADD(MUL(deltaB, _3rd(A)), MUL(_3rd(deltaA), SUB(B, deltaBvN)))));


	// these two sums, A' and B', are last because they were needed un-updated above:
	// A' = A - N*deltaA
	STORE(sums->A, SUB(A, deltaAvN));

	// B' = B - N*deltaB
	STORE(sums->B, SUB(B, deltaBvN));

}

static void adapt_sums(t_overfit_sums *sums1, t_overfit_sums *sums2)
{
	t_coord newStartA;
	t_coord newStartB;

	new_starts(sums1, sums2, &newStartA, &newStartB);

	modify_sums(sums1, newStartA, newStartB);
	modify_sums(sums2, newStartA, newStartB);
}

void overfit_sumadd_str(int token, t_overfit_sums *sums)
{
	t_overfit_sums sums1;
	overfit_sumsave_str(token,&sums1); // copying global sums to "sums1"

    t_overfit_sums sums2=*sums;

	if(!(ISZERO(SUB(LOAD(sums1.startA), LOAD(sums2.startA))) && ISZERO(SUB(LOAD(sums1.startB), LOAD(sums2.startB)))))
	{
		//sums, in general, have different "startsX" (primarily held
		//in global variables "vstartX[token]"), where X=A,B

		//calling "adapt_sums" is the default procedure. 
		//the alternative (i.e. "vstartX"s are the same) occurs 
		//in particular when in the function "overfit_setstart" 
		//"vstartX"s are set to "_zero"
		adapt_sums(&sums1, &sums2);
	}

	n_points[token]=sums1.N+sums2.N;
	
	vstartA[token]=LOAD(sums1.startA); // sums1.startX==sums2.startX by now
	vstartB[token]=LOAD(sums1.startB); // so it's immaterial which sum is here

	vSumA[token]=ADD(LOAD(sums1.A), LOAD(sums2.A));
	vSumB[token]=ADD(LOAD(sums1.B), LOAD(sums2.B));
	vSumAA[token]=ADD(LOAD(sums1.AA), LOAD(sums2.AA));
	vSumBB[token]=ADD(LOAD(sums1.BB), LOAD(sums2.BB));
	vSumAB0[token]=ADD(LOAD(sums1.AB0), LOAD(sums2.AB0));
	vSumAB1[token]=ADD(LOAD(sums1.AB1), LOAD(sums2.AB1));
	vSumAB2[token]=ADD(LOAD(sums1.AB2), LOAD(sums2.AB2));
}

void overfit_sumsave_str(int token, t_overfit_sums *sums)
{
	STORE(sums->startA, vstartA[token]);
	STORE(sums->startB, vstartB[token]);
	STORE(sums->A, vSumA[token]);
	STORE(sums->AA, vSumAA[token]);
	STORE(sums->B, vSumB[token]);
	STORE(sums->BB, vSumBB[token]);
	STORE(sums->AB0, vSumAB0[token]);
	STORE(sums->AB1, vSumAB1[token]);
	STORE(sums->AB2, vSumAB2[token]);

	sums->N=n_points[token];
}

float overfit(int N, float *xa, float *xb, t_transrot *TR)
{
	int token=overfit_get_token();

	overfit_reset(token);

	overfit_setstart(token, xa, xb);

	for(int i=0; i<N; i++) overfit_add(token, xa+4*i, xb+4*i);

	float res=fast_overfit(token, TR);

	overfit_release_token(token);

	return res;
}

/* vector functions using c arrays */

static inline void normalize(double a[3])
{
	double  b;

	b = sqrt((double)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
	a[0] /= b;
	a[1] /= b;
	a[2] /= b;
}

static inline double dot(double a[3], double b[3])
{
	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

static inline void cross(double a[3], double b[3], double c[3])
{
	a[0] = b[1]*c[2] - b[2]*c[1];
	a[1] = b[2]*c[0] - b[0]*c[2];
	a[2] = b[0]*c[1] - b[1]*c[0];
}


static void setup_rotation(int token, double R[3][3], double* E0)
{
	t_overfit_sums sums;

	overfit_sumsave_str(token, &sums);

	int n_list=sums.N;
	
	/* initialize */
	for(int i=0; i<3; i++) for(int j=0; j<3; j++) R[i][j] = 0.0;
	*E0 = 0.0;

	/* 
	 * E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n) 
	 */
	for(int i=0; i<3; i++) {
		*E0 += sums.AA[i]-sums.A[i]*sums.A[i]/n_list; 
		*E0 += sums.BB[i]-sums.B[i]*sums.B[i]/n_list; 
	}

	/*
	 * correlation matrix R:   
	 *   R[i,j) = sum(over n): y(n,i) * x(n,j)  
	 *   where x(n) and y(n) are two vector sets   
	 */
	for(int i=0; i<3; i++) {
		R[0][i] = sums.AB0[i]-sums.A[0]*sums.B[i]/n_list; 
		R[1][i] = sums.AB1[i]-sums.A[1]*sums.B[i]/n_list; 
		R[2][i] = sums.AB2[i]-sums.A[2]*sums.B[i]/n_list; 
	}

	*E0 *= 0.5;
}


#define ROTATE(a,i,j,k,l) { g = a[i][j]; h = a[k][l]; a[i][j] = g-s*(h+g*tau); a[k][l] = h+s*(g-h*tau); }
/*   
 * jacobi3
 *
 *    computes eigenval and eigen_vec of a real 3x3
 * symmetric matrix. On output, elements of a that are above 
 * the diagonal are destroyed. d[1..3] returns the 
 * eigenval of a. v[1..3][1..3] is a matrix whose 
 * columns contain, on output, the normalized eigen_vec of
 * a. n_rot returns the number of Jacobi rotations that were required.
 */
static int jacobi3(double a[3][3], double d[3], double v[3][3], int* n_rot)
{
	int count, k, i, j;
	double tresh, theta, tau, t, sum, s, h, g, c, b[3], z[3];

	/*Initialize v to the identity matrix.*/
	for (i=0; i<3; i++) { 
		for (j=0; j<3; j++) v[i][j] = 0.0;
		v[i][i] = 1.0;
	}

	/* Initialize b and d to the diagonal of a */
	for (i=0; i<3; i++) b[i] = d[i] = a[i][i];

	/* z will accumulate terms */
	for (i=0; i<3; i++) z[i] = 0.0; 

	*n_rot = 0;

	/* 50 tries */
	for (count=0; count<50; count++) {

		/* sum off-diagonal elements */
		sum = 0.0;
		for (i=0; i<2; i++) for (j=i+1; j<3; j++) sum += fabs(a[i][j]);

		/* if converged to machine underflow */
		if (sum == 0.0) return(1);

		/* on 1st three sweeps... */
		if (count < 3) tresh = sum * 0.2 / 9.0;    
		else tresh = 0.0;      

		for (i=0; i<2; i++) {
			for (j=i+1; j<3; j++) {
				g = 100.0 * fabs(a[i][j]);

				/*  after four sweeps, skip the rotation if
				 *   the off-diagonal element is small 
				 */
				if(count>3  &&  fabs(d[i])+g==fabs(d[i]) &&  fabs(d[j])+g==fabs(d[j])) {
					a[i][j] = 0.0;
				} else if(fabs(a[i][j]) > tresh) {
					h = d[j] - d[i];

					if (fabs(h)+g == fabs(h)) {
						t = a[i][j] / h;
					} else {
						theta = 0.5 * h / (a[i][j]);
						t = 1.0 / ( fabs(theta) + (double)sqrt(1.0 + theta*theta) );
						if (theta < 0.0) t = -t;
					}

					c = 1.0 / (double) sqrt(1 + t*t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[i][j];

					z[i] -= h;
					z[j] += h;
					d[i] -= h;
					d[j] += h;

					a[i][j] = 0.0;

					for (k=0; k<=i-1; k++) ROTATE(a, k, i, k, j);

					for (k=i+1; k<=j-1; k++) ROTATE(a, i, k, k, j);

					for (k=j+1; k<3; k++) ROTATE(a, i, k, j, k);

					for (k=0; k<3; k++) ROTATE(v, k, i, k, j);

					++(*n_rot);
				}
			}
		}

		for (i=0; i<3; i++) {
			b[i] += z[i];
			d[i] = b[i];
			z[i] = 0.0;
		}
	}

	printf("Too many iterations in jacobi3\n");
	return (0);
}  

/* 
 * diagonalize_symmetric 
 *
 *    Diagonalize a 3x3 matrix & sort eigenval by size
 */
static int diagonalize_symmetric(double matrix[3][3], double eigen_vec[3][3], double eigenval[3])
{
	int n_rot, i, j, k;
	double vec[3][3];
	double val; 

	if(!jacobi3(matrix, eigenval, vec, &n_rot)) {
		printf("convergence failed\n");
		return (0);
	}

	/* sort solutions by eigenval */
	for(i=0; i<3; i++) {
		k = i;
		val = eigenval[i];

		for(j=i+1; j<3; j++) if(eigenval[j] >= val)	{ 
			k = j;
			val = eigenval[k];
		}

		if(k != i) {
			eigenval[k] = eigenval[i];
			eigenval[i] = val;
			for (j=0; j<3; j++) {
				val = vec[j][i];
				vec[j][i] = vec[j][k];
				vec[j][k] = val;
			}
		}
	}

	/* transpose such that first index refers to solution index */
	for(i=0; i<3; i++) for(j=0; j<3; j++) eigen_vec[i][j] = vec[j][i];

	return (1);
}



/*
 * calculate_rotation_matrix() 
 *
 *   calculates the rotation matrix U and the
 * rmsd from the R matrix and E0:
 */
static int calculate_rotation_matrix(double R[3][3], double U[3][3], double E0, double* residual)
{
	int i, j, k;
	double Rt[3][3], RtR[3][3];
	double left_eigenvec[3][3], right_eigenvec[3][3], eigenval[3];
	double v[3];
	double sigma;

	/* build Rt, transpose of R  */
	for(i=0; i<3; i++) for(j=0; j<3; j++) Rt[i][j] = R[j][i];

	/* make symmetric RtR = Rt X R */
	for(i=0; i<3; i++) for(j=0; j<3; j++) {
		RtR[i][j] = 0.0;
		for (k = 0; k<3; k++) RtR[i][j] += Rt[k][i] * R[j][k];
	}

	if(!diagonalize_symmetric(RtR, right_eigenvec, eigenval)) return(0);

	/* right_eigenvec's should be an orthogonal system but could be left
	 * or right-handed. Let's force into right-handed system.
	 */
	cross(&right_eigenvec[2][0], &right_eigenvec[0][0], &right_eigenvec[1][0]);

	/* From the Kabsch algorithm, the eigenvec's of RtR
	 * are identical to the right_eigenvec's of R.
	 * This means that left_eigenvec = R x right_eigenvec 
	 */
	for(i=0; i<3; i++) for(j=0; j<3; j++) left_eigenvec[i][j] = dot(&right_eigenvec[i][0], &Rt[j][0]);

	for(i=0; i<3; i++) normalize(&left_eigenvec[i][0]);

	/* 
	 * Force left_eigenvec[2] to be orthogonal to the other vectors.
	 * First check if the rotational matrices generated from the 
	 * orthogonal eigenvectors are in a right-handed or left-handed
	 * co-ordinate system - given by sigma. Sigma is needed to
	 * resolve this ambiguity in calculating the RMSD.
	 */
	cross(v, &left_eigenvec[0][0], &left_eigenvec[1][0]);
	if (dot(v, &left_eigenvec[2][0]) < 0.0) sigma = -1.0;
	else sigma = 1.0;

	for (i=0; i<3; i++) left_eigenvec[2][i] = v[i]; 

	/* calc optimal rotation matrix U that minimises residual */
	for(i=0;i<3; i++) for(j=0; j<3; j++) {
		U[i][j] = 0.0;
		for(k=0; k<3; k++) U[i][j] += left_eigenvec[k][i] * right_eigenvec[k][j];
	}

	*residual = E0 - (double) sqrt(fabs(eigenval[0])) - (double) sqrt(fabs(eigenval[1])) - sigma * (double) sqrt(fabs(eigenval[2]));

	return (1);
}



float fast_overfit(int token, t_transrot *TR)
{
	double dU[3][3];
	double Eo;
    double residual=.0;
	double R[3][3];

	setup_rotation(token, R, &Eo);
	
	calculate_rotation_matrix(R, dU, Eo, &residual);

	if((TR!=NULL)) {
		t_overfit_sums sums;
		overfit_sumsave_str(token, &sums);

		float A[4];
		float B[4];

		memset(A, 0, sizeof(A));
		memset(A, 0, sizeof(A));

		for(int i=0; i<3; i++) {
			A[i]=sums.startA[i]+sums.A[i]/sums.N;
			B[i]=sums.startB[i]+sums.B[i]/sums.N;
		}

		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) TR->U[j][i]=(float) dU[j][i];
			TR->U[i][3]=0;

			TR->Tr[i]=A[i]-(B[0]*TR->U[0][i]+B[1]*TR->U[1][i]+B[2]*TR->U[2][i]);
		}
		
		TR->Tr[3]=0;
	}

	residual = fabs(residual); /* avoids the awkward case of -0.0 */
	

	if(isinf(residual)) {
		printf("thread %d: res: %lf n_points: %d\n", token, residual, n_points[token]);
		abort();
	}

	float res=(float)sqrt( fabs((double) (residual)*2.0/((double)n_points[token])) );

//	res=rintf(res*10000.0)*0.0001;
//	res=rintf(res*1000.0)*0.001;
//	res=rintf(res*100.0)*0.01;

	return res;
}


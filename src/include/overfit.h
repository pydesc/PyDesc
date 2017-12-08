/*
 * Copyright 2003-2010, 2017 Pawel Daniluk
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
 


#ifndef __OVERFIT_H__
#define __OVERFIT_H__


#ifdef __cplusplus
extern "C" {
#endif

typedef float v4sf __attribute__ ((vector_size(16)));
//typedef float v4sf __attribute__ ((mode(V4SF)));
typedef float float_al __attribute__ ((aligned(16)));

typedef struct {
	int N; 
	float startA[4]; // "cA" in "overfit_add" were taken in a coordinate 
	                //system with a center=="startA"
	float startB[4]; // "cB" in "overfit_add" were taken in a coordinate 
	                //system with a center=="startB"
	float A[4]; 
	float AA[4]; 
	float B[4]; 
	float BB[4]; 
	float AB0[4]; 
	float AB1[4]; 
	float AB2[4];
} t_overfit_sums;

typedef struct {
	float U[3][4];
	float Tr[4];
} t_transrot;

int overfit_get_token(void);
void overfit_release_token(int token);

float overfit_dist(float p1[4], float p2[4]);
void overfit_reset(int token);
void overfit_add(int token, float p1[4], float p2[4]);
void overfit_sumadd(int token, int N, float A[4], float AA[4], float B[4], float BB[4], float AB0[4], float AB1[4], float AB2[4]) __attribute__ ((deprecated));
void overfit_sumsave(int token, int *N, float A[4], float AA[4], float B[4], float BB[4], float AB0[4], float AB1[4], float AB2[4]) __attribute__ ((deprecated));
/*
float overfit_trans_dist(float p1[4], float p2[4], float U[3][4], float Tr[4]);
float overfit_trans_dist_sqr(float p1[4], float p2[4], float U[3][4], float Tr[4]);
void overfit_trans_combine(float U[3][4], float Tr[4], float U1[3][4], float Tr1[4], float U2[3][4], float Tr2[4]);
float fast_overfit(float U[3][4], float Tr[4]);
*/
float overfit_trans_dist(float p1[4], float p2[4], t_transrot TR);
float overfit_trans_dist_sqr(float p1[4], float p2[4], t_transrot TR);
void overfit_trans_combine(t_transrot *TR, t_transrot TR1, t_transrot TR2);
void overfit_trans(float p2[4], float p1[4], t_transrot TR);
float fast_overfit(int token, t_transrot *TR);
//float fast_overfit1(float U[3][4], float Tr[4]);

void overfit_sumadd_str(int token, t_overfit_sums *sums);
void overfit_sumsave_str(int token, t_overfit_sums *sums);

void overfit_setstart(int token, float p1[4], float p2[4]);

float overfit(int N, float *xa, float *xb, t_transrot *TR);
//float overfit(int N, float *xa, float *xb, float U[3][4], float Tr[4]);
//float overfit1(int N, float *xa, float *xb, float U[3][4], float Tr[4]);

void print_t_coord(float printed[4]);
void subtract_t_coord(float p1[4], float center[4], float out[4]);

#ifdef __cplusplus
};
#endif

#endif

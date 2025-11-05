/*
 * Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <assert.h>
#include <inttypes.h>
#include "sonLib.h"
#include "substitutionC.h"
#include "sequenceGraphC.h"
#include "xyzModelC.h"

#define STATE_NO 10

const int64_t D_1I_1M_1 = 0;
const int64_t I_1I_0D_0 = 1;
const int64_t I_1I_0M_0 = 2;
const int64_t I_1S_0S_0 = 3;
const int64_t M_1I_1D_1 = 4;
const int64_t M_1I_1M_1 = 5;
const int64_t W_0I_0I_1 = 6;
const int64_t W_0S_0I_1 = 7;
const int64_t X_0I_0I_1 = 8;
const int64_t D_1I_1D_1 = 9;

inline int64_t isSilent(int64_t state) {
	return FALSE;
}

inline int64_t isXYDelete(int64_t state) {
	if(state == D_1I_1D_1) return TRUE;
	return FALSE;
}

inline int64_t isXInsert(int64_t state) {
	if(state == I_1I_0D_0) return TRUE;
	if(state == I_1I_0M_0) return TRUE;
	if(state == I_1S_0S_0) return TRUE;
	return FALSE;
}

inline int64_t isYInsert(int64_t state) {
	if(state == W_0I_0I_1) return TRUE;
	if(state == W_0S_0I_1) return TRUE;
	if(state == X_0I_0I_1) return TRUE;
	return FALSE;
}

inline int64_t isXDelete(int64_t state) {
	if(state == D_1I_1M_1) return TRUE;
	return FALSE;
}

inline int64_t isYDelete(int64_t state) {
	if(state == M_1I_1D_1) return TRUE;
	return FALSE;
}

inline int64_t isMatch(int64_t state) {
	if(state == M_1I_1M_1) return TRUE;
	return FALSE;
}

struct CombinedTransitionModel *constructCombinedTransitionModel(float DX, float DZ, int64_t includeRoot, struct ParameterStruct *pS) {
	struct CombinedTransitionModel *temp = st_malloc(sizeof(struct CombinedTransitionModel));

	st_logInfo("Building combined transition model, DX: %f, DY %f\n", DX, DZ);
	temp->includeRoot = includeRoot;
	st_logInfo("Is root? " INT_STRING " \n", temp->includeRoot);
	temp->subModelX = constructHKYSubModel(DX,
		0.300000, 0.200000,
		0.200000, 0.300000,
		2.000000);

	temp->subModelY = constructHKYSubModel(DZ,
		0.300000, 0.200000,
		0.200000, 0.300000,
		2.000000);

	temp->ancestorProbs = st_malloc(sizeof(float)*ALPHABET_SIZE);
	temp->ancestorProbs[0] = 0.300000;
	temp->ancestorProbs[1] = 0.200000;
	temp->ancestorProbs[2] = 0.200000;
	temp->ancestorProbs[3] = 0.300000;

	int64_t i;
	int64_t j;
	for(i=0; i<ALPHABET_SIZE; i++) {
		st_logInfo("Stationary frequency %i value %f \n", i, temp->subModelX->stationaryDistribution[i]);
	}
	for(i=0; i<ALPHABET_SIZE; i++) {
		for(j=0; j<ALPHABET_SIZE; j++) {
			st_logInfo("Substitition probability (X Branch) %i %i, forward value %f, backward value %f\n", i, j, temp->subModelX->forward[i*ALPHABET_SIZE + j], temp->subModelX->backward[i*ALPHABET_SIZE + j]);
		}
	}
	for(i=0; i<ALPHABET_SIZE; i++) {
		for(j=0; j<ALPHABET_SIZE; j++) {
			st_logInfo("Substitition probability (Z Branch) %i %i, forward value %f, backward value %f\n", i, j, temp->subModelY->forward[i*ALPHABET_SIZE + j], temp->subModelY->backward[i*ALPHABET_SIZE + j]);
		}
	}
	float gO = pS->gO;
	st_logInfo("Parameter gO, value %f\n", gO);
	float gC = pS->gC;
	st_logInfo("Parameter gC, value %f\n", gC);
	float e = pS->e;
	st_logInfo("Parameter e, value %f\n", e);
	float aX = gO*DX;
	st_logInfo("Parameter aX, value %f\n", aX);
	float bX = gO*DX;
	st_logInfo("Parameter bX, value %f\n", bX);
	float gX = gC;
	st_logInfo("Parameter gX, value %f\n", gX);
	float dX = gC;
	st_logInfo("Parameter dX, value %f\n", dX);
	float maX = 1-aX;
	st_logInfo("Parameter maX, value %f\n", maX);
	float mbX = 1-bX;
	st_logInfo("Parameter mbX, value %f\n", mbX);
	float mgX = 1-gX;
	st_logInfo("Parameter mgX, value %f\n", mgX);
	float mdX = 1-dX;
	st_logInfo("Parameter mdX, value %f\n", mdX);
	float aZ = gO*DZ;
	st_logInfo("Parameter aZ, value %f\n", aZ);
	float bZ = gO*DZ;
	st_logInfo("Parameter bZ, value %f\n", bZ);
	float gZ = gC;
	st_logInfo("Parameter gZ, value %f\n", gZ);
	float dZ = gC;
	st_logInfo("Parameter dZ, value %f\n", dZ);
	float maZ = 1-aZ;
	st_logInfo("Parameter maZ, value %f\n", maZ);
	float mbZ = 1-bZ;
	st_logInfo("Parameter mbZ, value %f\n", mbZ);
	float mgZ = 1-gZ;
	st_logInfo("Parameter mgZ, value %f\n", mgZ);
	float mdZ = 1-dZ;
	st_logInfo("Parameter mdZ, value %f\n", mdZ);
	float me = 1-e;
	st_logInfo("Parameter me, value %f\n", me);
	temp->ftD_1I_1D_1_D_1I_1M_1 = LOG(gX*mgZ*maZ*mbZ*e);
	temp->tbtD_1I_1D_1_D_1I_1M_1 = LOG(gX*mgZ*maZ*mbZ*1.0);
	st_logInfo("From state D_1I_1D_1, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1D_1_D_1I_1M_1, temp->tbtD_1I_1D_1_D_1I_1M_1);
	temp->ftD_1I_1M_1_D_1I_1M_1 = LOG(gX*maZ*mbZ*e);
	temp->tbtD_1I_1M_1_D_1I_1M_1 = LOG(gX*maZ*mbZ*1.0);
	st_logInfo("From state D_1I_1M_1, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1M_1_D_1I_1M_1, temp->tbtD_1I_1M_1_D_1I_1M_1);
	temp->ftI_1I_0D_0_D_1I_1M_1 = LOG(mdX*mgZ*maZ*bX*mbZ*e);
	temp->tbtI_1I_0D_0_D_1I_1M_1 = LOG(mdX*mgZ*maZ*bX*mbZ*1.0);
	st_logInfo("From state I_1I_0D_0, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0D_0_D_1I_1M_1, temp->tbtI_1I_0D_0_D_1I_1M_1);
	temp->ftI_1I_0M_0_D_1I_1M_1 = LOG(mdX*maZ*bX*mbZ*e);
	temp->tbtI_1I_0M_0_D_1I_1M_1 = LOG(mdX*maZ*bX*mbZ*1.0);
	st_logInfo("From state I_1I_0M_0, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0M_0_D_1I_1M_1, temp->tbtI_1I_0M_0_D_1I_1M_1);
	temp->ftI_1S_0S_0_D_1I_1M_1 = LOG(mdX*maZ*bX*mbZ*e);
	temp->tbtI_1S_0S_0_D_1I_1M_1 = LOG(mdX*maZ*bX*mbZ*1.0);
	st_logInfo("From state I_1S_0S_0, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1S_0S_0_D_1I_1M_1, temp->tbtI_1S_0S_0_D_1I_1M_1);
	temp->ftM_1I_1D_1_D_1I_1M_1 = LOG(maX*mgZ*maZ*bX*mbZ*e);
	temp->tbtM_1I_1D_1_D_1I_1M_1 = LOG(maX*mgZ*maZ*bX*mbZ*1.0);
	st_logInfo("From state M_1I_1D_1, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1D_1_D_1I_1M_1, temp->tbtM_1I_1D_1_D_1I_1M_1);
	temp->ftM_1I_1M_1_D_1I_1M_1 = LOG(maX*maZ*bX*mbZ*e);
	temp->tbtM_1I_1M_1_D_1I_1M_1 = LOG(maX*maZ*bX*mbZ*1.0);
	st_logInfo("From state M_1I_1M_1, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1M_1_D_1I_1M_1, temp->tbtM_1I_1M_1_D_1I_1M_1);
	temp->ftS_0S_0S_0_D_1I_1M_1 = LOG(maX*maZ*bX*mbZ*e);
	temp->tbtS_0S_0S_0_D_1I_1M_1 = LOG(maX*maZ*bX*mbZ*1.0);
	st_logInfo("From state S_0S_0S_0, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftS_0S_0S_0_D_1I_1M_1, temp->tbtS_0S_0S_0_D_1I_1M_1);
	temp->ftW_0I_0I_1_D_1I_1M_1 = LOG(mdZ*bX*mbZ*e);
	temp->tbtW_0I_0I_1_D_1I_1M_1 = LOG(mdZ*bX*mbZ*1.0);
	st_logInfo("From state W_0I_0I_1, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0I_0I_1_D_1I_1M_1, temp->tbtW_0I_0I_1_D_1I_1M_1);
	temp->ftW_0S_0I_1_D_1I_1M_1 = LOG(mdZ*bX*mbZ*e);
	temp->tbtW_0S_0I_1_D_1I_1M_1 = LOG(mdZ*bX*mbZ*1.0);
	st_logInfo("From state W_0S_0I_1, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0S_0I_1_D_1I_1M_1, temp->tbtW_0S_0I_1_D_1I_1M_1);
	temp->ftX_0I_0I_1_D_1I_1M_1 = LOG(mdZ*mbZ*e);
	temp->tbtX_0I_0I_1_D_1I_1M_1 = LOG(mdZ*mbZ*1.0);
	st_logInfo("From state X_0I_0I_1, To state D_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftX_0I_0I_1_D_1I_1M_1, temp->tbtX_0I_0I_1_D_1I_1M_1);
	temp->ftD_1I_1D_1_E_0E_0E_0 = LOG(gX*gZ*me);
	temp->tbtD_1I_1D_1_E_0E_0E_0 = LOG(gX*gZ*1.0);
	st_logInfo("From state D_1I_1D_1, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1D_1_E_0E_0E_0, temp->tbtD_1I_1D_1_E_0E_0E_0);
	temp->ftD_1I_1M_1_E_0E_0E_0 = LOG(gX*maZ*me);
	temp->tbtD_1I_1M_1_E_0E_0E_0 = LOG(gX*maZ*1.0);
	st_logInfo("From state D_1I_1M_1, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1M_1_E_0E_0E_0, temp->tbtD_1I_1M_1_E_0E_0E_0);
	temp->ftI_1I_0D_0_E_0E_0E_0 = LOG(mdX*gZ*me);
	temp->tbtI_1I_0D_0_E_0E_0E_0 = LOG(mdX*gZ*1.0);
	st_logInfo("From state I_1I_0D_0, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0D_0_E_0E_0E_0, temp->tbtI_1I_0D_0_E_0E_0E_0);
	temp->ftI_1I_0M_0_E_0E_0E_0 = LOG(mdX*maZ*me);
	temp->tbtI_1I_0M_0_E_0E_0E_0 = LOG(mdX*maZ*1.0);
	st_logInfo("From state I_1I_0M_0, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0M_0_E_0E_0E_0, temp->tbtI_1I_0M_0_E_0E_0E_0);
	temp->ftI_1S_0S_0_E_0E_0E_0 = LOG(mdX*maZ*me);
	temp->tbtI_1S_0S_0_E_0E_0E_0 = LOG(mdX*maZ*1.0);
	st_logInfo("From state I_1S_0S_0, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1S_0S_0_E_0E_0E_0, temp->tbtI_1S_0S_0_E_0E_0E_0);
	temp->ftM_1I_1D_1_E_0E_0E_0 = LOG(maX*gZ*me);
	temp->tbtM_1I_1D_1_E_0E_0E_0 = LOG(maX*gZ*1.0);
	st_logInfo("From state M_1I_1D_1, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1D_1_E_0E_0E_0, temp->tbtM_1I_1D_1_E_0E_0E_0);
	temp->ftM_1I_1M_1_E_0E_0E_0 = LOG(maX*maZ*me);
	temp->tbtM_1I_1M_1_E_0E_0E_0 = LOG(maX*maZ*1.0);
	st_logInfo("From state M_1I_1M_1, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1M_1_E_0E_0E_0, temp->tbtM_1I_1M_1_E_0E_0E_0);
	temp->ftS_0S_0S_0_E_0E_0E_0 = LOG(maX*maZ*me);
	temp->tbtS_0S_0S_0_E_0E_0E_0 = LOG(maX*maZ*1.0);
	st_logInfo("From state S_0S_0S_0, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftS_0S_0S_0_E_0E_0E_0, temp->tbtS_0S_0S_0_E_0E_0E_0);
	temp->ftW_0I_0I_1_E_0E_0E_0 = LOG(mdZ*me);
	temp->tbtW_0I_0I_1_E_0E_0E_0 = LOG(mdZ*1.0);
	st_logInfo("From state W_0I_0I_1, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0I_0I_1_E_0E_0E_0, temp->tbtW_0I_0I_1_E_0E_0E_0);
	temp->ftW_0S_0I_1_E_0E_0E_0 = LOG(mdZ*me);
	temp->tbtW_0S_0I_1_E_0E_0E_0 = LOG(mdZ*1.0);
	st_logInfo("From state W_0S_0I_1, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0S_0I_1_E_0E_0E_0, temp->tbtW_0S_0I_1_E_0E_0E_0);
	temp->ftX_0I_0I_1_E_0E_0E_0 = LOG(mdZ*me);
	temp->tbtX_0I_0I_1_E_0E_0E_0 = LOG(mdZ*1.0);
	st_logInfo("From state X_0I_0I_1, To state E_0E_0E_0, Forward Parameter %f Traceback parameter %f \n", temp->ftX_0I_0I_1_E_0E_0E_0, temp->tbtX_0I_0I_1_E_0E_0E_0);
	temp->ftD_1I_1D_1_I_1I_0D_0 = LOG(mgX*aX);
	temp->tbtD_1I_1D_1_I_1I_0D_0 = LOG(mgX*aX);
	st_logInfo("From state D_1I_1D_1, To state I_1I_0D_0, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1D_1_I_1I_0D_0, temp->tbtD_1I_1D_1_I_1I_0D_0);
	temp->ftI_1I_0D_0_I_1I_0D_0 = LOG(dX);
	temp->tbtI_1I_0D_0_I_1I_0D_0 = LOG(dX);
	st_logInfo("From state I_1I_0D_0, To state I_1I_0D_0, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0D_0_I_1I_0D_0, temp->tbtI_1I_0D_0_I_1I_0D_0);
	temp->ftM_1I_1D_1_I_1I_0D_0 = LOG(aX);
	temp->tbtM_1I_1D_1_I_1I_0D_0 = LOG(aX);
	st_logInfo("From state M_1I_1D_1, To state I_1I_0D_0, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1D_1_I_1I_0D_0, temp->tbtM_1I_1D_1_I_1I_0D_0);
	temp->ftD_1I_1M_1_I_1I_0M_0 = LOG(mgX*aX);
	temp->tbtD_1I_1M_1_I_1I_0M_0 = LOG(mgX*aX);
	st_logInfo("From state D_1I_1M_1, To state I_1I_0M_0, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1M_1_I_1I_0M_0, temp->tbtD_1I_1M_1_I_1I_0M_0);
	temp->ftI_1I_0M_0_I_1I_0M_0 = LOG(dX);
	temp->tbtI_1I_0M_0_I_1I_0M_0 = LOG(dX);
	st_logInfo("From state I_1I_0M_0, To state I_1I_0M_0, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0M_0_I_1I_0M_0, temp->tbtI_1I_0M_0_I_1I_0M_0);
	temp->ftM_1I_1M_1_I_1I_0M_0 = LOG(aX);
	temp->tbtM_1I_1M_1_I_1I_0M_0 = LOG(aX);
	st_logInfo("From state M_1I_1M_1, To state I_1I_0M_0, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1M_1_I_1I_0M_0, temp->tbtM_1I_1M_1_I_1I_0M_0);
	temp->ftI_1S_0S_0_I_1S_0S_0 = LOG(dX);
	temp->tbtI_1S_0S_0_I_1S_0S_0 = LOG(dX);
	st_logInfo("From state I_1S_0S_0, To state I_1S_0S_0, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1S_0S_0_I_1S_0S_0, temp->tbtI_1S_0S_0_I_1S_0S_0);
	temp->ftS_0S_0S_0_I_1S_0S_0 = LOG(aX);
	temp->tbtS_0S_0S_0_I_1S_0S_0 = LOG(aX);
	st_logInfo("From state S_0S_0S_0, To state I_1S_0S_0, Forward Parameter %f Traceback parameter %f \n", temp->ftS_0S_0S_0_I_1S_0S_0, temp->tbtS_0S_0S_0_I_1S_0S_0);
	temp->ftD_1I_1D_1_M_1I_1D_1 = LOG(mgX*maX*gZ*mbX*e);
	temp->tbtD_1I_1D_1_M_1I_1D_1 = LOG(mgX*maX*gZ*mbX*1.0);
	st_logInfo("From state D_1I_1D_1, To state M_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1D_1_M_1I_1D_1, temp->tbtD_1I_1D_1_M_1I_1D_1);
	temp->ftD_1I_1M_1_M_1I_1D_1 = LOG(mgX*maX*maZ*mbX*bZ*e);
	temp->tbtD_1I_1M_1_M_1I_1D_1 = LOG(mgX*maX*maZ*mbX*bZ*1.0);
	st_logInfo("From state D_1I_1M_1, To state M_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1M_1_M_1I_1D_1, temp->tbtD_1I_1M_1_M_1I_1D_1);
	temp->ftI_1I_0D_0_M_1I_1D_1 = LOG(mdX*gZ*mbX*e);
	temp->tbtI_1I_0D_0_M_1I_1D_1 = LOG(mdX*gZ*mbX*1.0);
	st_logInfo("From state I_1I_0D_0, To state M_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0D_0_M_1I_1D_1, temp->tbtI_1I_0D_0_M_1I_1D_1);
	temp->ftI_1I_0M_0_M_1I_1D_1 = LOG(mdX*maZ*mbX*bZ*e);
	temp->tbtI_1I_0M_0_M_1I_1D_1 = LOG(mdX*maZ*mbX*bZ*1.0);
	st_logInfo("From state I_1I_0M_0, To state M_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0M_0_M_1I_1D_1, temp->tbtI_1I_0M_0_M_1I_1D_1);
	temp->ftI_1S_0S_0_M_1I_1D_1 = LOG(mdX*maZ*mbX*bZ*e);
	temp->tbtI_1S_0S_0_M_1I_1D_1 = LOG(mdX*maZ*mbX*bZ*1.0);
	st_logInfo("From state I_1S_0S_0, To state M_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1S_0S_0_M_1I_1D_1, temp->tbtI_1S_0S_0_M_1I_1D_1);
	temp->ftM_1I_1D_1_M_1I_1D_1 = LOG(maX*gZ*mbX*e);
	temp->tbtM_1I_1D_1_M_1I_1D_1 = LOG(maX*gZ*mbX*1.0);
	st_logInfo("From state M_1I_1D_1, To state M_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1D_1_M_1I_1D_1, temp->tbtM_1I_1D_1_M_1I_1D_1);
	temp->ftM_1I_1M_1_M_1I_1D_1 = LOG(maX*maZ*mbX*bZ*e);
	temp->tbtM_1I_1M_1_M_1I_1D_1 = LOG(maX*maZ*mbX*bZ*1.0);
	st_logInfo("From state M_1I_1M_1, To state M_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1M_1_M_1I_1D_1, temp->tbtM_1I_1M_1_M_1I_1D_1);
	temp->ftS_0S_0S_0_M_1I_1D_1 = LOG(maX*maZ*mbX*bZ*e);
	temp->tbtS_0S_0S_0_M_1I_1D_1 = LOG(maX*maZ*mbX*bZ*1.0);
	st_logInfo("From state S_0S_0S_0, To state M_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftS_0S_0S_0_M_1I_1D_1, temp->tbtS_0S_0S_0_M_1I_1D_1);
	temp->ftW_0I_0I_1_M_1I_1D_1 = LOG(mdZ*mbX*bZ*e);
	temp->tbtW_0I_0I_1_M_1I_1D_1 = LOG(mdZ*mbX*bZ*1.0);
	st_logInfo("From state W_0I_0I_1, To state M_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0I_0I_1_M_1I_1D_1, temp->tbtW_0I_0I_1_M_1I_1D_1);
	temp->ftW_0S_0I_1_M_1I_1D_1 = LOG(mdZ*mbX*bZ*e);
	temp->tbtW_0S_0I_1_M_1I_1D_1 = LOG(mdZ*mbX*bZ*1.0);
	st_logInfo("From state W_0S_0I_1, To state M_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0S_0I_1_M_1I_1D_1, temp->tbtW_0S_0I_1_M_1I_1D_1);
	temp->ftD_1I_1D_1_M_1I_1M_1 = LOG(mgX*maX*mgZ*maZ*mbX*mbZ*e);
	temp->tbtD_1I_1D_1_M_1I_1M_1 = LOG(mgX*maX*mgZ*maZ*mbX*mbZ*1.0);
	st_logInfo("From state D_1I_1D_1, To state M_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1D_1_M_1I_1M_1, temp->tbtD_1I_1D_1_M_1I_1M_1);
	temp->ftD_1I_1M_1_M_1I_1M_1 = LOG(mgX*maX*maZ*mbX*mbZ*e);
	temp->tbtD_1I_1M_1_M_1I_1M_1 = LOG(mgX*maX*maZ*mbX*mbZ*1.0);
	st_logInfo("From state D_1I_1M_1, To state M_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1M_1_M_1I_1M_1, temp->tbtD_1I_1M_1_M_1I_1M_1);
	temp->ftI_1I_0D_0_M_1I_1M_1 = LOG(mdX*mgZ*maZ*mbX*mbZ*e);
	temp->tbtI_1I_0D_0_M_1I_1M_1 = LOG(mdX*mgZ*maZ*mbX*mbZ*1.0);
	st_logInfo("From state I_1I_0D_0, To state M_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0D_0_M_1I_1M_1, temp->tbtI_1I_0D_0_M_1I_1M_1);
	temp->ftI_1I_0M_0_M_1I_1M_1 = LOG(mdX*maZ*mbX*mbZ*e);
	temp->tbtI_1I_0M_0_M_1I_1M_1 = LOG(mdX*maZ*mbX*mbZ*1.0);
	st_logInfo("From state I_1I_0M_0, To state M_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0M_0_M_1I_1M_1, temp->tbtI_1I_0M_0_M_1I_1M_1);
	temp->ftI_1S_0S_0_M_1I_1M_1 = LOG(mdX*maZ*mbX*mbZ*e);
	temp->tbtI_1S_0S_0_M_1I_1M_1 = LOG(mdX*maZ*mbX*mbZ*1.0);
	st_logInfo("From state I_1S_0S_0, To state M_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1S_0S_0_M_1I_1M_1, temp->tbtI_1S_0S_0_M_1I_1M_1);
	temp->ftM_1I_1D_1_M_1I_1M_1 = LOG(maX*mgZ*maZ*mbX*mbZ*e);
	temp->tbtM_1I_1D_1_M_1I_1M_1 = LOG(maX*mgZ*maZ*mbX*mbZ*1.0);
	st_logInfo("From state M_1I_1D_1, To state M_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1D_1_M_1I_1M_1, temp->tbtM_1I_1D_1_M_1I_1M_1);
	temp->ftM_1I_1M_1_M_1I_1M_1 = LOG(maX*maZ*mbX*mbZ*e);
	temp->tbtM_1I_1M_1_M_1I_1M_1 = LOG(maX*maZ*mbX*mbZ*1.0);
	st_logInfo("From state M_1I_1M_1, To state M_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1M_1_M_1I_1M_1, temp->tbtM_1I_1M_1_M_1I_1M_1);
	temp->ftS_0S_0S_0_M_1I_1M_1 = LOG(maX*maZ*mbX*mbZ*e);
	temp->tbtS_0S_0S_0_M_1I_1M_1 = LOG(maX*maZ*mbX*mbZ*1.0);
	st_logInfo("From state S_0S_0S_0, To state M_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftS_0S_0S_0_M_1I_1M_1, temp->tbtS_0S_0S_0_M_1I_1M_1);
	temp->ftW_0I_0I_1_M_1I_1M_1 = LOG(mdZ*mbX*mbZ*e);
	temp->tbtW_0I_0I_1_M_1I_1M_1 = LOG(mdZ*mbX*mbZ*1.0);
	st_logInfo("From state W_0I_0I_1, To state M_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0I_0I_1_M_1I_1M_1, temp->tbtW_0I_0I_1_M_1I_1M_1);
	temp->ftW_0S_0I_1_M_1I_1M_1 = LOG(mdZ*mbX*mbZ*e);
	temp->tbtW_0S_0I_1_M_1I_1M_1 = LOG(mdZ*mbX*mbZ*1.0);
	st_logInfo("From state W_0S_0I_1, To state M_1I_1M_1, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0S_0I_1_M_1I_1M_1, temp->tbtW_0S_0I_1_M_1I_1M_1);
	temp->ftD_1I_1D_1_W_0I_0I_1 = LOG(mgX*maX*mgZ*aZ);
	temp->tbtD_1I_1D_1_W_0I_0I_1 = LOG(mgX*maX*mgZ*aZ);
	st_logInfo("From state D_1I_1D_1, To state W_0I_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1D_1_W_0I_0I_1, temp->tbtD_1I_1D_1_W_0I_0I_1);
	temp->ftD_1I_1M_1_W_0I_0I_1 = LOG(mgX*maX*aZ);
	temp->tbtD_1I_1M_1_W_0I_0I_1 = LOG(mgX*maX*aZ);
	st_logInfo("From state D_1I_1M_1, To state W_0I_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1M_1_W_0I_0I_1, temp->tbtD_1I_1M_1_W_0I_0I_1);
	temp->ftI_1I_0D_0_W_0I_0I_1 = LOG(mdX*mgZ*aZ);
	temp->tbtI_1I_0D_0_W_0I_0I_1 = LOG(mdX*mgZ*aZ);
	st_logInfo("From state I_1I_0D_0, To state W_0I_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0D_0_W_0I_0I_1, temp->tbtI_1I_0D_0_W_0I_0I_1);
	temp->ftI_1I_0M_0_W_0I_0I_1 = LOG(mdX*aZ);
	temp->tbtI_1I_0M_0_W_0I_0I_1 = LOG(mdX*aZ);
	st_logInfo("From state I_1I_0M_0, To state W_0I_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0M_0_W_0I_0I_1, temp->tbtI_1I_0M_0_W_0I_0I_1);
	temp->ftM_1I_1D_1_W_0I_0I_1 = LOG(maX*mgZ*aZ);
	temp->tbtM_1I_1D_1_W_0I_0I_1 = LOG(maX*mgZ*aZ);
	st_logInfo("From state M_1I_1D_1, To state W_0I_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1D_1_W_0I_0I_1, temp->tbtM_1I_1D_1_W_0I_0I_1);
	temp->ftM_1I_1M_1_W_0I_0I_1 = LOG(maX*aZ);
	temp->tbtM_1I_1M_1_W_0I_0I_1 = LOG(maX*aZ);
	st_logInfo("From state M_1I_1M_1, To state W_0I_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1M_1_W_0I_0I_1, temp->tbtM_1I_1M_1_W_0I_0I_1);
	temp->ftW_0I_0I_1_W_0I_0I_1 = LOG(dZ);
	temp->tbtW_0I_0I_1_W_0I_0I_1 = LOG(dZ);
	st_logInfo("From state W_0I_0I_1, To state W_0I_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0I_0I_1_W_0I_0I_1, temp->tbtW_0I_0I_1_W_0I_0I_1);
	temp->ftI_1S_0S_0_W_0S_0I_1 = LOG(mdX*aZ);
	temp->tbtI_1S_0S_0_W_0S_0I_1 = LOG(mdX*aZ);
	st_logInfo("From state I_1S_0S_0, To state W_0S_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1S_0S_0_W_0S_0I_1, temp->tbtI_1S_0S_0_W_0S_0I_1);
	temp->ftS_0S_0S_0_W_0S_0I_1 = LOG(maX*aZ);
	temp->tbtS_0S_0S_0_W_0S_0I_1 = LOG(maX*aZ);
	st_logInfo("From state S_0S_0S_0, To state W_0S_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftS_0S_0S_0_W_0S_0I_1, temp->tbtS_0S_0S_0_W_0S_0I_1);
	temp->ftW_0S_0I_1_W_0S_0I_1 = LOG(dZ);
	temp->tbtW_0S_0I_1_W_0S_0I_1 = LOG(dZ);
	st_logInfo("From state W_0S_0I_1, To state W_0S_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0S_0I_1_W_0S_0I_1, temp->tbtW_0S_0I_1_W_0S_0I_1);
	temp->ftD_1I_1D_1_X_0I_0I_1 = LOG(gX*mgZ*aZ);
	temp->tbtD_1I_1D_1_X_0I_0I_1 = LOG(gX*mgZ*aZ);
	st_logInfo("From state D_1I_1D_1, To state X_0I_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1D_1_X_0I_0I_1, temp->tbtD_1I_1D_1_X_0I_0I_1);
	temp->ftD_1I_1M_1_X_0I_0I_1 = LOG(gX*aZ);
	temp->tbtD_1I_1M_1_X_0I_0I_1 = LOG(gX*aZ);
	st_logInfo("From state D_1I_1M_1, To state X_0I_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1M_1_X_0I_0I_1, temp->tbtD_1I_1M_1_X_0I_0I_1);
	temp->ftX_0I_0I_1_X_0I_0I_1 = LOG(dZ);
	temp->tbtX_0I_0I_1_X_0I_0I_1 = LOG(dZ);
	st_logInfo("From state X_0I_0I_1, To state X_0I_0I_1, Forward Parameter %f Traceback parameter %f \n", temp->ftX_0I_0I_1_X_0I_0I_1, temp->tbtX_0I_0I_1_X_0I_0I_1);
	temp->ftD_1I_1D_1_D_1I_1D_1 = LOG(gX*gZ*e);
	temp->tbtD_1I_1D_1_D_1I_1D_1 = LOG(gX*gZ*1.0);
	st_logInfo("From state D_1I_1D_1, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1D_1_D_1I_1D_1, temp->tbtD_1I_1D_1_D_1I_1D_1);
	temp->ftD_1I_1M_1_D_1I_1D_1 = LOG(gX*maZ*bZ*e);
	temp->tbtD_1I_1M_1_D_1I_1D_1 = LOG(gX*maZ*bZ*1.0);
	st_logInfo("From state D_1I_1M_1, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftD_1I_1M_1_D_1I_1D_1, temp->tbtD_1I_1M_1_D_1I_1D_1);
	temp->ftI_1I_0D_0_D_1I_1D_1 = LOG(mdX*gZ*bX*e);
	temp->tbtI_1I_0D_0_D_1I_1D_1 = LOG(mdX*gZ*bX*1.0);
	st_logInfo("From state I_1I_0D_0, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0D_0_D_1I_1D_1, temp->tbtI_1I_0D_0_D_1I_1D_1);
	temp->ftI_1I_0M_0_D_1I_1D_1 = LOG(mdX*maZ*bX*bZ*e);
	temp->tbtI_1I_0M_0_D_1I_1D_1 = LOG(mdX*maZ*bX*bZ*1.0);
	st_logInfo("From state I_1I_0M_0, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1I_0M_0_D_1I_1D_1, temp->tbtI_1I_0M_0_D_1I_1D_1);
	temp->ftI_1S_0S_0_D_1I_1D_1 = LOG(mdX*maZ*bX*bZ*e);
	temp->tbtI_1S_0S_0_D_1I_1D_1 = LOG(mdX*maZ*bX*bZ*1.0);
	st_logInfo("From state I_1S_0S_0, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftI_1S_0S_0_D_1I_1D_1, temp->tbtI_1S_0S_0_D_1I_1D_1);
	temp->ftM_1I_1D_1_D_1I_1D_1 = LOG(maX*gZ*bX*e);
	temp->tbtM_1I_1D_1_D_1I_1D_1 = LOG(maX*gZ*bX*1.0);
	st_logInfo("From state M_1I_1D_1, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1D_1_D_1I_1D_1, temp->tbtM_1I_1D_1_D_1I_1D_1);
	temp->ftM_1I_1M_1_D_1I_1D_1 = LOG(maX*maZ*bX*bZ*e);
	temp->tbtM_1I_1M_1_D_1I_1D_1 = LOG(maX*maZ*bX*bZ*1.0);
	st_logInfo("From state M_1I_1M_1, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftM_1I_1M_1_D_1I_1D_1, temp->tbtM_1I_1M_1_D_1I_1D_1);
	temp->ftS_0S_0S_0_D_1I_1D_1 = LOG(maX*maZ*bX*bZ*e);
	temp->tbtS_0S_0S_0_D_1I_1D_1 = LOG(maX*maZ*bX*bZ*1.0);
	st_logInfo("From state S_0S_0S_0, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftS_0S_0S_0_D_1I_1D_1, temp->tbtS_0S_0S_0_D_1I_1D_1);
	temp->ftW_0I_0I_1_D_1I_1D_1 = LOG(mdZ*bX*bZ*e);
	temp->tbtW_0I_0I_1_D_1I_1D_1 = LOG(mdZ*bX*bZ*1.0);
	st_logInfo("From state W_0I_0I_1, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0I_0I_1_D_1I_1D_1, temp->tbtW_0I_0I_1_D_1I_1D_1);
	temp->ftW_0S_0I_1_D_1I_1D_1 = LOG(mdZ*bX*bZ*e);
	temp->tbtW_0S_0I_1_D_1I_1D_1 = LOG(mdZ*bX*bZ*1.0);
	st_logInfo("From state W_0S_0I_1, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftW_0S_0I_1_D_1I_1D_1, temp->tbtW_0S_0I_1_D_1I_1D_1);
	temp->ftX_0I_0I_1_D_1I_1D_1 = LOG(mdZ*bZ*e);
	temp->tbtX_0I_0I_1_D_1I_1D_1 = LOG(mdZ*bZ*1.0);
	st_logInfo("From state X_0I_0I_1, To state D_1I_1D_1, Forward Parameter %f Traceback parameter %f \n", temp->ftX_0I_0I_1_D_1I_1D_1, temp->tbtX_0I_0I_1_D_1I_1D_1);

	float tltD_1I_1D_1_D_1I_1D_1 = gX*gZ*e;
	 float fA[STATE_NO];
	 float fA1[STATE_NO];
	 float fA2[STATE_NO];
	float *l1;
	float *l2;
	float *l3;
	l1 = fA1;
	l2 = fA2;
	for(i=0; i<STATE_NO; i++) {
		l1[i] = 0.0;
		l2[i] = 0.0;
		fA[i] = 0.0;
	}
	l2[D_1I_1D_1] = 1.0;
	for(i=0; i<1000; i++) {
		fA[D_1I_1D_1] += l1[D_1I_1D_1];
		l2[D_1I_1D_1] += l1[D_1I_1D_1] * tltD_1I_1D_1_D_1I_1D_1;
		l3 = l1;
		l1 = l2;
		l2 = l3;
		for(j=0; j<STATE_NO; j++)
			l2[j] = 0.0;
	}
	temp->ltD_1I_1D_1_D_1I_1D_1 = LOG(fA[D_1I_1D_1]);
	st_logInfo("From state D_1I_1D_1, To state D_1I_1D_1, LOOP PARAMETER %f \n", temp->ltD_1I_1D_1_D_1I_1D_1);
	return temp;
}

inline void insertXFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge,void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float)) {
	assignFn(aDS, D_1I_1D_1, I_1I_0D_0, model->ftD_1I_1D_1_I_1I_0D_0 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1I_0D_0, I_1I_0D_0, model->ftI_1I_0D_0_I_1I_0D_0 + edge->edgeScore + edge->subScore);
	assignFn(aDS, M_1I_1D_1, I_1I_0D_0, model->ftM_1I_1D_1_I_1I_0D_0 + edge->edgeScore + edge->subScore);
	assignFn(aDS, D_1I_1M_1, I_1I_0M_0, model->ftD_1I_1M_1_I_1I_0M_0 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1I_0M_0, I_1I_0M_0, model->ftI_1I_0M_0_I_1I_0M_0 + edge->edgeScore + edge->subScore);
	assignFn(aDS, M_1I_1M_1, I_1I_0M_0, model->ftM_1I_1M_1_I_1I_0M_0 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1S_0S_0, I_1S_0S_0, model->ftI_1S_0S_0_I_1S_0S_0 + edge->edgeScore + edge->subScore);
}

inline void insertXFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge,void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float, float)) {
	if (model->includeRoot) {
		float i;
		i = model->ftD_1I_1D_1_I_1I_0D_0 + edge->edgeScore + edge->subScore;
		assignFn(aDS, D_1I_1D_1, I_1I_0D_0, i, i);
		i = model->ftI_1I_0D_0_I_1I_0D_0 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1I_0D_0, I_1I_0D_0, i, i);
		i = model->ftM_1I_1D_1_I_1I_0D_0 + edge->edgeScore + edge->subScore;
		assignFn(aDS, M_1I_1D_1, I_1I_0D_0, i, i);
		i = model->ftD_1I_1M_1_I_1I_0M_0 + edge->edgeScore + edge->subScore;
		assignFn(aDS, D_1I_1M_1, I_1I_0M_0, i, i);
		i = model->ftI_1I_0M_0_I_1I_0M_0 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1I_0M_0, I_1I_0M_0, i, i);
		i = model->ftM_1I_1M_1_I_1I_0M_0 + edge->edgeScore + edge->subScore;
		assignFn(aDS, M_1I_1M_1, I_1I_0M_0, i, i);
		i = model->ftI_1S_0S_0_I_1S_0S_0 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1S_0S_0, I_1S_0S_0, i, i);
	}
	else {
		assignFn(aDS, D_1I_1D_1, I_1I_0D_0, model->ftD_1I_1D_1_I_1I_0D_0 + edge->edgeScore + edge->subScore, model->tbtD_1I_1D_1_I_1I_0D_0 + edge->edgeScore);
		assignFn(aDS, I_1I_0D_0, I_1I_0D_0, model->ftI_1I_0D_0_I_1I_0D_0 + edge->edgeScore + edge->subScore, model->tbtI_1I_0D_0_I_1I_0D_0 + edge->edgeScore);
		assignFn(aDS, M_1I_1D_1, I_1I_0D_0, model->ftM_1I_1D_1_I_1I_0D_0 + edge->edgeScore + edge->subScore, model->tbtM_1I_1D_1_I_1I_0D_0 + edge->edgeScore);
		assignFn(aDS, D_1I_1M_1, I_1I_0M_0, model->ftD_1I_1M_1_I_1I_0M_0 + edge->edgeScore + edge->subScore, model->tbtD_1I_1M_1_I_1I_0M_0 + edge->edgeScore);
		assignFn(aDS, I_1I_0M_0, I_1I_0M_0, model->ftI_1I_0M_0_I_1I_0M_0 + edge->edgeScore + edge->subScore, model->tbtI_1I_0M_0_I_1I_0M_0 + edge->edgeScore);
		assignFn(aDS, M_1I_1M_1, I_1I_0M_0, model->ftM_1I_1M_1_I_1I_0M_0 + edge->edgeScore + edge->subScore, model->tbtM_1I_1M_1_I_1I_0M_0 + edge->edgeScore);
		assignFn(aDS, I_1S_0S_0, I_1S_0S_0, model->ftI_1S_0S_0_I_1S_0S_0 + edge->edgeScore + edge->subScore, model->tbtI_1S_0S_0_I_1S_0S_0 + edge->edgeScore);
	}
}

inline void insertYFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge,void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float)) {
	assignFn(aDS, D_1I_1D_1, W_0I_0I_1, model->ftD_1I_1D_1_W_0I_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, D_1I_1M_1, W_0I_0I_1, model->ftD_1I_1M_1_W_0I_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1I_0D_0, W_0I_0I_1, model->ftI_1I_0D_0_W_0I_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1I_0M_0, W_0I_0I_1, model->ftI_1I_0M_0_W_0I_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, M_1I_1D_1, W_0I_0I_1, model->ftM_1I_1D_1_W_0I_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, M_1I_1M_1, W_0I_0I_1, model->ftM_1I_1M_1_W_0I_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, W_0I_0I_1, W_0I_0I_1, model->ftW_0I_0I_1_W_0I_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1S_0S_0, W_0S_0I_1, model->ftI_1S_0S_0_W_0S_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, W_0S_0I_1, W_0S_0I_1, model->ftW_0S_0I_1_W_0S_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, D_1I_1D_1, X_0I_0I_1, model->ftD_1I_1D_1_X_0I_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, D_1I_1M_1, X_0I_0I_1, model->ftD_1I_1M_1_X_0I_0I_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, X_0I_0I_1, X_0I_0I_1, model->ftX_0I_0I_1_X_0I_0I_1 + edge->edgeScore + edge->subScore);
}

inline void insertYFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge,void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float, float)) {
	if (model->includeRoot) {
		float i;
		i = model->ftD_1I_1D_1_W_0I_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, D_1I_1D_1, W_0I_0I_1, i, i);
		i = model->ftD_1I_1M_1_W_0I_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, D_1I_1M_1, W_0I_0I_1, i, i);
		i = model->ftI_1I_0D_0_W_0I_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1I_0D_0, W_0I_0I_1, i, i);
		i = model->ftI_1I_0M_0_W_0I_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1I_0M_0, W_0I_0I_1, i, i);
		i = model->ftM_1I_1D_1_W_0I_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, M_1I_1D_1, W_0I_0I_1, i, i);
		i = model->ftM_1I_1M_1_W_0I_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, M_1I_1M_1, W_0I_0I_1, i, i);
		i = model->ftW_0I_0I_1_W_0I_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, W_0I_0I_1, W_0I_0I_1, i, i);
		i = model->ftI_1S_0S_0_W_0S_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1S_0S_0, W_0S_0I_1, i, i);
		i = model->ftW_0S_0I_1_W_0S_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, W_0S_0I_1, W_0S_0I_1, i, i);
		i = model->ftD_1I_1D_1_X_0I_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, D_1I_1D_1, X_0I_0I_1, i, i);
		i = model->ftD_1I_1M_1_X_0I_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, D_1I_1M_1, X_0I_0I_1, i, i);
		i = model->ftX_0I_0I_1_X_0I_0I_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, X_0I_0I_1, X_0I_0I_1, i, i);
	}
	else {
		assignFn(aDS, D_1I_1D_1, W_0I_0I_1, model->ftD_1I_1D_1_W_0I_0I_1 + edge->edgeScore + edge->subScore, model->tbtD_1I_1D_1_W_0I_0I_1 + edge->edgeScore);
		assignFn(aDS, D_1I_1M_1, W_0I_0I_1, model->ftD_1I_1M_1_W_0I_0I_1 + edge->edgeScore + edge->subScore, model->tbtD_1I_1M_1_W_0I_0I_1 + edge->edgeScore);
		assignFn(aDS, I_1I_0D_0, W_0I_0I_1, model->ftI_1I_0D_0_W_0I_0I_1 + edge->edgeScore + edge->subScore, model->tbtI_1I_0D_0_W_0I_0I_1 + edge->edgeScore);
		assignFn(aDS, I_1I_0M_0, W_0I_0I_1, model->ftI_1I_0M_0_W_0I_0I_1 + edge->edgeScore + edge->subScore, model->tbtI_1I_0M_0_W_0I_0I_1 + edge->edgeScore);
		assignFn(aDS, M_1I_1D_1, W_0I_0I_1, model->ftM_1I_1D_1_W_0I_0I_1 + edge->edgeScore + edge->subScore, model->tbtM_1I_1D_1_W_0I_0I_1 + edge->edgeScore);
		assignFn(aDS, M_1I_1M_1, W_0I_0I_1, model->ftM_1I_1M_1_W_0I_0I_1 + edge->edgeScore + edge->subScore, model->tbtM_1I_1M_1_W_0I_0I_1 + edge->edgeScore);
		assignFn(aDS, W_0I_0I_1, W_0I_0I_1, model->ftW_0I_0I_1_W_0I_0I_1 + edge->edgeScore + edge->subScore, model->tbtW_0I_0I_1_W_0I_0I_1 + edge->edgeScore);
		assignFn(aDS, I_1S_0S_0, W_0S_0I_1, model->ftI_1S_0S_0_W_0S_0I_1 + edge->edgeScore + edge->subScore, model->tbtI_1S_0S_0_W_0S_0I_1 + edge->edgeScore);
		assignFn(aDS, W_0S_0I_1, W_0S_0I_1, model->ftW_0S_0I_1_W_0S_0I_1 + edge->edgeScore + edge->subScore, model->tbtW_0S_0I_1_W_0S_0I_1 + edge->edgeScore);
		assignFn(aDS, D_1I_1D_1, X_0I_0I_1, model->ftD_1I_1D_1_X_0I_0I_1 + edge->edgeScore + edge->subScore, model->tbtD_1I_1D_1_X_0I_0I_1 + edge->edgeScore);
		assignFn(aDS, D_1I_1M_1, X_0I_0I_1, model->ftD_1I_1M_1_X_0I_0I_1 + edge->edgeScore + edge->subScore, model->tbtD_1I_1M_1_X_0I_0I_1 + edge->edgeScore);
		assignFn(aDS, X_0I_0I_1, X_0I_0I_1, model->ftX_0I_0I_1_X_0I_0I_1 + edge->edgeScore + edge->subScore, model->tbtX_0I_0I_1_X_0I_0I_1 + edge->edgeScore);
	}
}

inline void deleteXFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge,void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float)) {
	assignFn(aDS, D_1I_1D_1, D_1I_1M_1, model->ftD_1I_1D_1_D_1I_1M_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, D_1I_1M_1, D_1I_1M_1, model->ftD_1I_1M_1_D_1I_1M_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1I_0D_0, D_1I_1M_1, model->ftI_1I_0D_0_D_1I_1M_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1I_0M_0, D_1I_1M_1, model->ftI_1I_0M_0_D_1I_1M_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1S_0S_0, D_1I_1M_1, model->ftI_1S_0S_0_D_1I_1M_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, M_1I_1D_1, D_1I_1M_1, model->ftM_1I_1D_1_D_1I_1M_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, M_1I_1M_1, D_1I_1M_1, model->ftM_1I_1M_1_D_1I_1M_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, W_0I_0I_1, D_1I_1M_1, model->ftW_0I_0I_1_D_1I_1M_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, W_0S_0I_1, D_1I_1M_1, model->ftW_0S_0I_1_D_1I_1M_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, X_0I_0I_1, D_1I_1M_1, model->ftX_0I_0I_1_D_1I_1M_1 + edge->edgeScore + edge->subScore);
}

inline void deleteXFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge,void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float, float)) {
	if (model->includeRoot) {
		float i;
		i = model->ftD_1I_1D_1_D_1I_1M_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, D_1I_1D_1, D_1I_1M_1, i, i);
		i = model->ftD_1I_1M_1_D_1I_1M_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, D_1I_1M_1, D_1I_1M_1, i, i);
		i = model->ftI_1I_0D_0_D_1I_1M_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1I_0D_0, D_1I_1M_1, i, i);
		i = model->ftI_1I_0M_0_D_1I_1M_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1I_0M_0, D_1I_1M_1, i, i);
		i = model->ftI_1S_0S_0_D_1I_1M_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1S_0S_0, D_1I_1M_1, i, i);
		i = model->ftM_1I_1D_1_D_1I_1M_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, M_1I_1D_1, D_1I_1M_1, i, i);
		i = model->ftM_1I_1M_1_D_1I_1M_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, M_1I_1M_1, D_1I_1M_1, i, i);
		i = model->ftW_0I_0I_1_D_1I_1M_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, W_0I_0I_1, D_1I_1M_1, i, i);
		i = model->ftW_0S_0I_1_D_1I_1M_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, W_0S_0I_1, D_1I_1M_1, i, i);
		i = model->ftX_0I_0I_1_D_1I_1M_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, X_0I_0I_1, D_1I_1M_1, i, i);
	}
	else {
		assignFn(aDS, D_1I_1D_1, D_1I_1M_1, model->ftD_1I_1D_1_D_1I_1M_1 + edge->edgeScore + edge->subScore, model->tbtD_1I_1D_1_D_1I_1M_1 + edge->edgeScore);
		assignFn(aDS, D_1I_1M_1, D_1I_1M_1, model->ftD_1I_1M_1_D_1I_1M_1 + edge->edgeScore + edge->subScore, model->tbtD_1I_1M_1_D_1I_1M_1 + edge->edgeScore);
		assignFn(aDS, I_1I_0D_0, D_1I_1M_1, model->ftI_1I_0D_0_D_1I_1M_1 + edge->edgeScore + edge->subScore, model->tbtI_1I_0D_0_D_1I_1M_1 + edge->edgeScore);
		assignFn(aDS, I_1I_0M_0, D_1I_1M_1, model->ftI_1I_0M_0_D_1I_1M_1 + edge->edgeScore + edge->subScore, model->tbtI_1I_0M_0_D_1I_1M_1 + edge->edgeScore);
		assignFn(aDS, I_1S_0S_0, D_1I_1M_1, model->ftI_1S_0S_0_D_1I_1M_1 + edge->edgeScore + edge->subScore, model->tbtI_1S_0S_0_D_1I_1M_1 + edge->edgeScore);
		assignFn(aDS, M_1I_1D_1, D_1I_1M_1, model->ftM_1I_1D_1_D_1I_1M_1 + edge->edgeScore + edge->subScore, model->tbtM_1I_1D_1_D_1I_1M_1 + edge->edgeScore);
		assignFn(aDS, M_1I_1M_1, D_1I_1M_1, model->ftM_1I_1M_1_D_1I_1M_1 + edge->edgeScore + edge->subScore, model->tbtM_1I_1M_1_D_1I_1M_1 + edge->edgeScore);
		assignFn(aDS, W_0I_0I_1, D_1I_1M_1, model->ftW_0I_0I_1_D_1I_1M_1 + edge->edgeScore + edge->subScore, model->tbtW_0I_0I_1_D_1I_1M_1 + edge->edgeScore);
		assignFn(aDS, W_0S_0I_1, D_1I_1M_1, model->ftW_0S_0I_1_D_1I_1M_1 + edge->edgeScore + edge->subScore, model->tbtW_0S_0I_1_D_1I_1M_1 + edge->edgeScore);
		assignFn(aDS, X_0I_0I_1, D_1I_1M_1, model->ftX_0I_0I_1_D_1I_1M_1 + edge->edgeScore + edge->subScore, model->tbtX_0I_0I_1_D_1I_1M_1 + edge->edgeScore);
	}
}

inline void deleteYFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge,void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float)) {
	assignFn(aDS, D_1I_1D_1, M_1I_1D_1, model->ftD_1I_1D_1_M_1I_1D_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, D_1I_1M_1, M_1I_1D_1, model->ftD_1I_1M_1_M_1I_1D_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1I_0D_0, M_1I_1D_1, model->ftI_1I_0D_0_M_1I_1D_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1I_0M_0, M_1I_1D_1, model->ftI_1I_0M_0_M_1I_1D_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, I_1S_0S_0, M_1I_1D_1, model->ftI_1S_0S_0_M_1I_1D_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, M_1I_1D_1, M_1I_1D_1, model->ftM_1I_1D_1_M_1I_1D_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, M_1I_1M_1, M_1I_1D_1, model->ftM_1I_1M_1_M_1I_1D_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, W_0I_0I_1, M_1I_1D_1, model->ftW_0I_0I_1_M_1I_1D_1 + edge->edgeScore + edge->subScore);
	assignFn(aDS, W_0S_0I_1, M_1I_1D_1, model->ftW_0S_0I_1_M_1I_1D_1 + edge->edgeScore + edge->subScore);
}

inline void deleteYFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge,void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float, float)) {
	if (model->includeRoot) {
		float i;
		i = model->ftD_1I_1D_1_M_1I_1D_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, D_1I_1D_1, M_1I_1D_1, i, i);
		i = model->ftD_1I_1M_1_M_1I_1D_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, D_1I_1M_1, M_1I_1D_1, i, i);
		i = model->ftI_1I_0D_0_M_1I_1D_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1I_0D_0, M_1I_1D_1, i, i);
		i = model->ftI_1I_0M_0_M_1I_1D_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1I_0M_0, M_1I_1D_1, i, i);
		i = model->ftI_1S_0S_0_M_1I_1D_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, I_1S_0S_0, M_1I_1D_1, i, i);
		i = model->ftM_1I_1D_1_M_1I_1D_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, M_1I_1D_1, M_1I_1D_1, i, i);
		i = model->ftM_1I_1M_1_M_1I_1D_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, M_1I_1M_1, M_1I_1D_1, i, i);
		i = model->ftW_0I_0I_1_M_1I_1D_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, W_0I_0I_1, M_1I_1D_1, i, i);
		i = model->ftW_0S_0I_1_M_1I_1D_1 + edge->edgeScore + edge->subScore;
		assignFn(aDS, W_0S_0I_1, M_1I_1D_1, i, i);
	}
	else {
		assignFn(aDS, D_1I_1D_1, M_1I_1D_1, model->ftD_1I_1D_1_M_1I_1D_1 + edge->edgeScore + edge->subScore, model->tbtD_1I_1D_1_M_1I_1D_1 + edge->edgeScore);
		assignFn(aDS, D_1I_1M_1, M_1I_1D_1, model->ftD_1I_1M_1_M_1I_1D_1 + edge->edgeScore + edge->subScore, model->tbtD_1I_1M_1_M_1I_1D_1 + edge->edgeScore);
		assignFn(aDS, I_1I_0D_0, M_1I_1D_1, model->ftI_1I_0D_0_M_1I_1D_1 + edge->edgeScore + edge->subScore, model->tbtI_1I_0D_0_M_1I_1D_1 + edge->edgeScore);
		assignFn(aDS, I_1I_0M_0, M_1I_1D_1, model->ftI_1I_0M_0_M_1I_1D_1 + edge->edgeScore + edge->subScore, model->tbtI_1I_0M_0_M_1I_1D_1 + edge->edgeScore);
		assignFn(aDS, I_1S_0S_0, M_1I_1D_1, model->ftI_1S_0S_0_M_1I_1D_1 + edge->edgeScore + edge->subScore, model->tbtI_1S_0S_0_M_1I_1D_1 + edge->edgeScore);
		assignFn(aDS, M_1I_1D_1, M_1I_1D_1, model->ftM_1I_1D_1_M_1I_1D_1 + edge->edgeScore + edge->subScore, model->tbtM_1I_1D_1_M_1I_1D_1 + edge->edgeScore);
		assignFn(aDS, M_1I_1M_1, M_1I_1D_1, model->ftM_1I_1M_1_M_1I_1D_1 + edge->edgeScore + edge->subScore, model->tbtM_1I_1M_1_M_1I_1D_1 + edge->edgeScore);
		assignFn(aDS, W_0I_0I_1, M_1I_1D_1, model->ftW_0I_0I_1_M_1I_1D_1 + edge->edgeScore + edge->subScore, model->tbtW_0I_0I_1_M_1I_1D_1 + edge->edgeScore);
		assignFn(aDS, W_0S_0I_1, M_1I_1D_1, model->ftW_0S_0I_1_M_1I_1D_1 + edge->edgeScore + edge->subScore, model->tbtW_0S_0I_1_M_1I_1D_1 + edge->edgeScore);
	}
}

inline void matchFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model,struct Edge *edgeX, struct Edge *edgeY, void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float)) {
	 float m[ALPHABET_SIZE];
	 float i;
	 float j;
	multiplyWV(edgeX->wV, edgeY->wV, m, ALPHABET_SIZE);
	i = LOG(combineWV(m, model->ancestorProbs, ALPHABET_SIZE));
	j = model->ftD_1I_1D_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
	assignFn(aDS, D_1I_1D_1, M_1I_1M_1, j);
	j = model->ftD_1I_1M_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
	assignFn(aDS, D_1I_1M_1, M_1I_1M_1, j);
	j = model->ftI_1I_0D_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
	assignFn(aDS, I_1I_0D_0, M_1I_1M_1, j);
	j = model->ftI_1I_0M_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
	assignFn(aDS, I_1I_0M_0, M_1I_1M_1, j);
	j = model->ftI_1S_0S_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
	assignFn(aDS, I_1S_0S_0, M_1I_1M_1, j);
	j = model->ftM_1I_1D_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
	assignFn(aDS, M_1I_1D_1, M_1I_1M_1, j);
	j = model->ftM_1I_1M_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
	assignFn(aDS, M_1I_1M_1, M_1I_1M_1, j);
	j = model->ftW_0I_0I_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
	assignFn(aDS, W_0I_0I_1, M_1I_1M_1, j);
	j = model->ftW_0S_0I_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
	assignFn(aDS, W_0S_0I_1, M_1I_1M_1, j);
}

inline void matchFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model,struct Edge *edgeX, struct Edge *edgeY, void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float, float)) {
	 float m[ALPHABET_SIZE];
	 float i;
	 float j;
	 float k;
	 float l;
	multiplyWV(edgeX->wV, edgeY->wV, m, ALPHABET_SIZE);
	i = LOG(sumWV(m, ALPHABET_SIZE));
	j = LOG(combineWV(m, model->ancestorProbs, ALPHABET_SIZE));
	if (model->includeRoot) {
		l = model->ftD_1I_1D_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, D_1I_1D_1, M_1I_1M_1, l, l);
		l = model->ftD_1I_1M_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, D_1I_1M_1, M_1I_1M_1, l, l);
		l = model->ftI_1I_0D_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, I_1I_0D_0, M_1I_1M_1, l, l);
		l = model->ftI_1I_0M_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, I_1I_0M_0, M_1I_1M_1, l, l);
		l = model->ftI_1S_0S_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, I_1S_0S_0, M_1I_1M_1, l, l);
		l = model->ftM_1I_1D_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, M_1I_1D_1, M_1I_1M_1, l, l);
		l = model->ftM_1I_1M_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, M_1I_1M_1, M_1I_1M_1, l, l);
		l = model->ftW_0I_0I_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, W_0I_0I_1, M_1I_1M_1, l, l);
		l = model->ftW_0S_0I_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, W_0S_0I_1, M_1I_1M_1, l, l);
	}
	else {
		k = model->tbtD_1I_1D_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
		l = model->ftD_1I_1D_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, D_1I_1D_1, M_1I_1M_1, l, k);
		k = model->tbtD_1I_1M_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
		l = model->ftD_1I_1M_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, D_1I_1M_1, M_1I_1M_1, l, k);
		k = model->tbtI_1I_0D_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
		l = model->ftI_1I_0D_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, I_1I_0D_0, M_1I_1M_1, l, k);
		k = model->tbtI_1I_0M_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
		l = model->ftI_1I_0M_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, I_1I_0M_0, M_1I_1M_1, l, k);
		k = model->tbtI_1S_0S_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
		l = model->ftI_1S_0S_0_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, I_1S_0S_0, M_1I_1M_1, l, k);
		k = model->tbtM_1I_1D_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
		l = model->ftM_1I_1D_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, M_1I_1D_1, M_1I_1M_1, l, k);
		k = model->tbtM_1I_1M_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
		l = model->ftM_1I_1M_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, M_1I_1M_1, M_1I_1M_1, l, k);
		k = model->tbtW_0I_0I_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
		l = model->ftW_0I_0I_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, W_0I_0I_1, M_1I_1M_1, l, k);
		k = model->tbtW_0S_0I_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + i;
		l = model->ftW_0S_0I_1_M_1I_1M_1 + edgeX->edgeScore + edgeY->edgeScore + j;
		assignFn(aDS, W_0S_0I_1, M_1I_1M_1, l, k);
	}
}

inline void silentFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, float *cell,void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float)) {
	assignFn(aDS, D_1I_1M_1, D_1I_1D_1, model->ftD_1I_1M_1_D_1I_1D_1);
	assignFn(aDS, I_1I_0D_0, D_1I_1D_1, model->ftI_1I_0D_0_D_1I_1D_1);
	assignFn(aDS, I_1I_0M_0, D_1I_1D_1, model->ftI_1I_0M_0_D_1I_1D_1);
	assignFn(aDS, I_1S_0S_0, D_1I_1D_1, model->ftI_1S_0S_0_D_1I_1D_1);
	assignFn(aDS, M_1I_1D_1, D_1I_1D_1, model->ftM_1I_1D_1_D_1I_1D_1);
	assignFn(aDS, M_1I_1M_1, D_1I_1D_1, model->ftM_1I_1M_1_D_1I_1D_1);
	assignFn(aDS, W_0I_0I_1, D_1I_1D_1, model->ftW_0I_0I_1_D_1I_1D_1);
	assignFn(aDS, W_0S_0I_1, D_1I_1D_1, model->ftW_0S_0I_1_D_1I_1D_1);
	assignFn(aDS, X_0I_0I_1, D_1I_1D_1, model->ftX_0I_0I_1_D_1I_1D_1);
	 float fA[STATE_NO];
	float i;
	fA[D_1I_1D_1] = cell[D_1I_1D_1];
	i = LOG_ZERO;
	LOG_PLUS_EQUALS(&i, fA[D_1I_1D_1] + model->ltD_1I_1D_1_D_1I_1D_1);
	cell[D_1I_1D_1] = i;
}

inline void silentFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float, float)) {
	if (model->includeRoot) {
	}
	else {
	}
}

inline void deleteFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, int64_t, int64_t, float, float)) {
	if (model->includeRoot) {
		assignFn(aDS, D_1I_1D_1, D_1I_1D_1, model->ftD_1I_1D_1_D_1I_1D_1, model->ftD_1I_1D_1_D_1I_1D_1);
		assignFn(aDS, D_1I_1M_1, D_1I_1D_1, model->ftD_1I_1M_1_D_1I_1D_1, model->ftD_1I_1M_1_D_1I_1D_1);
		assignFn(aDS, I_1I_0D_0, D_1I_1D_1, model->ftI_1I_0D_0_D_1I_1D_1, model->ftI_1I_0D_0_D_1I_1D_1);
		assignFn(aDS, I_1I_0M_0, D_1I_1D_1, model->ftI_1I_0M_0_D_1I_1D_1, model->ftI_1I_0M_0_D_1I_1D_1);
		assignFn(aDS, I_1S_0S_0, D_1I_1D_1, model->ftI_1S_0S_0_D_1I_1D_1, model->ftI_1S_0S_0_D_1I_1D_1);
		assignFn(aDS, M_1I_1D_1, D_1I_1D_1, model->ftM_1I_1D_1_D_1I_1D_1, model->ftM_1I_1D_1_D_1I_1D_1);
		assignFn(aDS, M_1I_1M_1, D_1I_1D_1, model->ftM_1I_1M_1_D_1I_1D_1, model->ftM_1I_1M_1_D_1I_1D_1);
		assignFn(aDS, W_0I_0I_1, D_1I_1D_1, model->ftW_0I_0I_1_D_1I_1D_1, model->ftW_0I_0I_1_D_1I_1D_1);
		assignFn(aDS, W_0S_0I_1, D_1I_1D_1, model->ftW_0S_0I_1_D_1I_1D_1, model->ftW_0S_0I_1_D_1I_1D_1);
		assignFn(aDS, X_0I_0I_1, D_1I_1D_1, model->ftX_0I_0I_1_D_1I_1D_1, model->ftX_0I_0I_1_D_1I_1D_1);
	}
	else {
		assignFn(aDS, D_1I_1D_1, D_1I_1D_1, model->ftD_1I_1D_1_D_1I_1D_1, model->tbtD_1I_1D_1_D_1I_1D_1);
		assignFn(aDS, D_1I_1M_1, D_1I_1D_1, model->ftD_1I_1M_1_D_1I_1D_1, model->tbtD_1I_1M_1_D_1I_1D_1);
		assignFn(aDS, I_1I_0D_0, D_1I_1D_1, model->ftI_1I_0D_0_D_1I_1D_1, model->tbtI_1I_0D_0_D_1I_1D_1);
		assignFn(aDS, I_1I_0M_0, D_1I_1D_1, model->ftI_1I_0M_0_D_1I_1D_1, model->tbtI_1I_0M_0_D_1I_1D_1);
		assignFn(aDS, I_1S_0S_0, D_1I_1D_1, model->ftI_1S_0S_0_D_1I_1D_1, model->tbtI_1S_0S_0_D_1I_1D_1);
		assignFn(aDS, M_1I_1D_1, D_1I_1D_1, model->ftM_1I_1D_1_D_1I_1D_1, model->tbtM_1I_1D_1_D_1I_1D_1);
		assignFn(aDS, M_1I_1M_1, D_1I_1D_1, model->ftM_1I_1M_1_D_1I_1D_1, model->tbtM_1I_1M_1_D_1I_1D_1);
		assignFn(aDS, W_0I_0I_1, D_1I_1D_1, model->ftW_0I_0I_1_D_1I_1D_1, model->tbtW_0I_0I_1_D_1I_1D_1);
		assignFn(aDS, W_0S_0I_1, D_1I_1D_1, model->ftW_0S_0I_1_D_1I_1D_1, model->tbtW_0S_0I_1_D_1I_1D_1);
		assignFn(aDS, X_0I_0I_1, D_1I_1D_1, model->ftX_0I_0I_1_D_1I_1D_1, model->tbtX_0I_0I_1_D_1I_1D_1);
	}
}

void destructCombinedTransitionModel(struct CombinedTransitionModel *model) { free(model); }

int64_t stateNo() { return STATE_NO; }

float *startStates(struct CombinedTransitionModel *model) {
	float *i; int64_t j; i = st_malloc(sizeof(float)*STATE_NO);
	for (j=0;j<STATE_NO; j++) i[j] = LOG_ZERO;
	i[D_1I_1D_1] = model->ftS_0S_0S_0_D_1I_1D_1;
	i[D_1I_1M_1] = model->ftS_0S_0S_0_D_1I_1M_1;
	i[I_1S_0S_0] = model->ftS_0S_0S_0_I_1S_0S_0;
	i[M_1I_1D_1] = model->ftS_0S_0S_0_M_1I_1D_1;
	i[M_1I_1M_1] = model->ftS_0S_0S_0_M_1I_1M_1;
	i[W_0S_0I_1] = model->ftS_0S_0S_0_W_0S_0I_1;
	return i;
}

float *endStates(struct CombinedTransitionModel *model) {
	float *i; int64_t j; i = st_malloc(sizeof(float)*STATE_NO);
	for (j=0;j<STATE_NO; j++) i[j] = LOG_ZERO;

	i[D_1I_1D_1] = model->ftD_1I_1D_1_E_0E_0E_0;
	i[D_1I_1M_1] = model->ftD_1I_1M_1_E_0E_0E_0;
	i[I_1I_0D_0] = model->ftI_1I_0D_0_E_0E_0E_0;
	i[I_1I_0M_0] = model->ftI_1I_0M_0_E_0E_0E_0;
	i[I_1S_0S_0] = model->ftI_1S_0S_0_E_0E_0E_0;
	i[M_1I_1D_1] = model->ftM_1I_1D_1_E_0E_0E_0;
	i[M_1I_1M_1] = model->ftM_1I_1M_1_E_0E_0E_0;
	i[W_0I_0I_1] = model->ftW_0I_0I_1_E_0E_0E_0;
	i[W_0S_0I_1] = model->ftW_0S_0I_1_E_0E_0E_0;
	i[X_0I_0I_1] = model->ftX_0I_0I_1_E_0E_0E_0;
	return i;
}

struct ParameterStruct *constructParamStruct(int argc, char *argv[]) {
	int64_t i;
	char *mod;
	struct ParameterStruct *pM = st_malloc(sizeof(struct ParameterStruct));
	float floatParser;
	pM->gO = 0.0025;
	pM->gC = 0.7;
	pM->e = 0.999;
	for(i=1; i<argc; i++) {
		mod = argv[i];
		if(mod[0] != '-') {
			continue;
		}
		assert(mod[0] == '-');
		switch(mod[1]) {
			case 'A':
				sscanf(argv[++i], "%f", &floatParser);
				pM->gO = floatParser;
			break;
			case 'B':
				sscanf(argv[++i], "%f", &floatParser);
				pM->gC = floatParser;
			break;
			case 'C':
				sscanf(argv[++i], "%f", &floatParser);
				pM->e = floatParser;
			break;
		}
	}
	return pM;
}

void destructParamStruct(struct ParameterStruct *pM) {
	free(pM);
}

void printParamStruct() {
	fprintf(stderr, "	-A [FLOAT] value of parameter gO , default: 0.0025 \n");
	fprintf(stderr, "	-B [FLOAT] value of parameter gC , default: 0.7 \n");
	fprintf(stderr, "	-C [FLOAT] value of parameter e , default: 0.999 \n");
}


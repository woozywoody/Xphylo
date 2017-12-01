#include "uthash.h"

#ifndef _GTREE_H_
#define _GTREE_H_

#define NAMELEN 100
typedef struct _tag_Family Family;
typedef struct _tag_SnpKey SnpKey;
typedef struct _tag_Pos2Snp Pos2Snp;
typedef struct _tag_Snp2Family Snp2Family;
typedef struct _tag_nodeArray NodeArray;


struct _tag_nodeArray{
	int n;
	int maxn;
	Family* nodeArray;
};

struct _tag_Family{
	char familyId[NAMELEN];
	Family* children;
	Family* brother;
	Family* parent;
	int tmpweight;
	double weight;
	double pathWeight;
	double forwardPathWeight;
	UT_hash_handle hh;
};

struct _tag_SnpKey{
	int pos;
	char alt[NAMELEN];
};

struct _tag_Snp2Family{
	char snp[NAMELEN];
	Family** familyP;
	int ns; //N family share 
	int mns; //At most mns family share
	SnpKey snpkey;
	UT_hash_handle hh;
};

struct _tag_Pos2Snp{
	SnpKey snpkey;
	Snp2Family** snpsP;
	int ns; //ns aliases
	int mns; //At most mns aliases
	UT_hash_handle hh;
};

NodeArray Ytree;

typedef void (Family_Printf)(Family*);

int addFamily( Family* parent, Family* current );

int familyDisplay( Family* node, Family_Printf* pFunc, int gap, char div );

int familyDisplayDot( NodeArray* P, char* filename );

Family* getLikelyFamily( NodeArray* T );

int setSnpsFamily( Snp2Family* snp, Family* cNode );

int setSnpsName( Pos2Snp* psnp, Snp2Family* fsnp );

int snpsFamilyWeight1( Pos2Snp* psnp );

int initFamilyWeight( NodeArray* P );

int setForwardPathWight( Family* node );

int setPathWight( Family* node );

int freePosHash( Pos2Snp* P );

int freeSnpHash( Snp2Family* P );

int freeFamilyHash( Family* P );

int freeYtree( NodeArray* P );

int showPath( Family* P);


extern int Verbose;
#endif 

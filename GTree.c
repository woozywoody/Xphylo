#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include "uthash.h"
#include "GTree.h"

int addFamily( Family* parent, Family* current ){

	Family* tmpF;
	current->parent = parent;
	if( ! parent )//Insert root
		return 0;

	if( ! parent->children )
		parent->children = current; //The first son
	else{ //Other sons
		tmpF = parent->children;
		while( tmpF->brother )
			tmpF = tmpF->brother;
		tmpF->brother = current;
	}	
	return 0;	
}

//Print Tree
int familyDisplay( Family* node, Family_Printf* pFunc, int gap, char div )
{
	int i = 0;
	Family* cnode;

	if( (node != NULL) && (pFunc != NULL) )
	{
		for( i=0; i<gap; i++ )
			fprintf(stderr,"%c ", div);
	
		pFunc(node);
	
		cnode = node->children;
		while( cnode ){
			familyDisplay( cnode, pFunc, gap+1, div);
			cnode = cnode->brother;
		}
	}
	return 0;
}

//Display dot tree
int familyDisplayDotR( Family* node, FILE* f )
{
	Family* cnode;

	if( node && node->children )
	{
		cnode = node->children;
		while( cnode ){
			fprintf( f,"\"%s:%.2lf\" -> \"%s:%.2lf\";\n", node->familyId, node->weight, cnode->familyId, cnode->weight );
			familyDisplayDotR( cnode, f );
			cnode = cnode->brother;
		}
	}
	return 0;
}

//Display dot tree
int familyDisplayDot( NodeArray* P, char* filename )
{
	FILE* f;
	int i;
	Family* node = P->nodeArray;
	
	if( (f=fopen(filename,"w")) == NULL ){
		fprintf( stderr, "Can not open file %s.\n",filename );
		return 1;
	}

	fprintf( f, "digraph G {\nrankdir=LR\n" );
	familyDisplayDotR( node, f );
	
	for( i=0; i<P->n; i++ ){
		if( node[i].weight > 0 )
			fprintf( f, "\"%s:%.2lf\" [color=red, style=filled];\n", node[i].familyId, node[i].weight );
	}
	fprintf( f, "}\n" );

	return 0;
}

//Set path_weight
int setPathWight( Family* node )
{
	Family* cnode;
	if( node )
	{
		if( node->parent )
			node->pathWeight = node->weight + node->parent->pathWeight;
		else
			node->pathWeight = node->weight;

		cnode = node->children;
		while( cnode ){
			setPathWight( cnode );
			cnode = cnode->brother;
		}
	}
	return 0;
}

//Set path_weight
int setForwardPathWight( Family* node )
{
	Family* cnode;
	if( node )
	{
		cnode = node->children;
		while( cnode ){
			setForwardPathWight( cnode );
			node->forwardPathWeight += cnode->forwardPathWeight + cnode->weight;
			cnode = cnode->brother;
		}
	}
	return 0;
}

//Clear flags for next sample
int initFamilyWeight( NodeArray* P )
{
	int i;
	Family* node = P->nodeArray;

	for( i=0; i<P->n; i++,node++ ){
		node->weight = 0;
		node->tmpweight = 0;
		node->pathWeight = 0;
		node->forwardPathWeight = 0;
	}
	return 0;
}

//Compare function for qsort
static int weightCmp(const void *A, const void *B)
{
	Family* a=((Family**)A)[0];
	Family* b=((Family**)B)[0];
	return  b->pathWeight - a->pathWeight > 0 ? 1 : -1;
}

//Clear flag tmpweight on path from root to P
int setTmpweightOnPath( Family* P, int n )
{
	Family* tmpP = P;
	while( tmpP ){
		tmpP->tmpweight = n;
		tmpP = tmpP->parent;
	}
	return 0;
}

//If there are tow paths with the same weight, trace their common ancestor.
Family* traceAncestor(Family* A, Family* B)
{
	Family* ancestor = A;
	setTmpweightOnPath( A, 0 );
	setTmpweightOnPath( B, 1 );
	//use tmpweight as a flag. 11111111111110000000000000000000   A
	//                                      \1111111111111111     B
	while( ancestor->tmpweight == 0 )
		ancestor = ancestor->parent;

	return ancestor;
}

//Find out the family with the maximum weight.
Family* getLikelyFamily( NodeArray* T )
{
	int i = 0;
	Family* first;
	Family* p=T->nodeArray;
	Family** tmpT=(Family**)malloc( sizeof(Family*) * T->n );

	while( i<T->n )
		tmpT[i++] = p++;

	qsort( tmpT, T->n, sizeof(Family*), weightCmp );

	first = tmpT[0];
	for( i=1; i<T->n; i++ ){//Get common ancestor. use first as ancestor.
		if( fabs( tmpT[i-1]->pathWeight - tmpT[i]->pathWeight) > 0.1 )
			break;
		else
			first = traceAncestor( first, tmpT[i] );
	}
	free(tmpT);

	return first;
}

//Set snp's family. A snp may be shared by many families.
int setSnpsFamily( Snp2Family* snp, Family* cNode )
{
	int i;
	
	for( i=0; i < snp->ns; i++ )
		if( snp->familyP[i] == cNode ) //Already exist
			return 0;

	if( snp->ns >= snp->mns ){
		snp->mns += 2;
		snp->familyP = (Family**)realloc( snp->familyP, sizeof(Family *) * snp->mns );
	}

	snp->familyP[snp->ns++] = cNode;
	return 0;
}

//Set snp's name. A snp may have many aliases.
int setSnpsName( Pos2Snp* psnp, Snp2Family* fsnp )
{
	int i;
	
	for( i=0; i < psnp->ns; i++ )
		if( psnp->snpsP[i] == fsnp ) //Already exist
			return 0;

	if( psnp->ns >= psnp->mns ){
		psnp->mns += 2;
		psnp->snpsP = (Snp2Family **)realloc( psnp->snpsP, sizeof(Snp2Family *) * psnp->mns );
	}

	psnp->snpsP[psnp->ns++] = fsnp;
	return 0;
}

//Weight++
int snpsFamilyWeight1( Pos2Snp* psnp )
{
	int i,j;
	int uniqscore = 0;
	double rv_uniqscore = 0;
	Snp2Family* fsnp;

	for( j=0; j<psnp->ns; j++ ){
		fsnp = psnp->snpsP[j];
		for( i=0; i<fsnp->ns; i++ )
			fsnp->familyP[i]->tmpweight++;
	}

	for( j=0; j<psnp->ns; j++ ){ //Won't plus more than 1 on each node.
		fsnp = psnp->snpsP[j];
		for( i=0; i<fsnp->ns; i++ ){
			if( fsnp->familyP[i]->tmpweight >0 ){
				uniqscore++;
			}
		}
	}

	rv_uniqscore = (double)1/uniqscore;
	for( j=0; j<psnp->ns; j++ ){ //Won't plus more than 1 on each node.
		fsnp = psnp->snpsP[j];
		for( i=0; i<fsnp->ns; i++ ){
			if( fsnp->familyP[i]->tmpweight >0 ){
				fsnp->familyP[i]->weight += rv_uniqscore;
				fsnp->familyP[i]->tmpweight = 0;
			}
		}
	}
	return 0;
}

//Show path
int showPath( Family* P )
{
	if( !P )
		return 0;

	showPath( P->parent );
	printf(" -> %s",P->familyId);
	return 0;
}

//Free hash
int freePosHash( Pos2Snp* P )
{
	Pos2Snp *s, *tmp;
	HASH_ITER( hh, P, s, tmp ){
		HASH_DEL( P, s );
		free( s->snpsP );
		free( s );
	}
	return 0;
}

//Free hash
int freeSnpHash( Snp2Family* P )
{
	Snp2Family *s, *tmp;
	HASH_ITER( hh, P, s, tmp ){
		HASH_DEL( P, s );
		free( s->familyP );
		free( s );
	}
	return 0;
}

//Free hash
int freeFamilyHash( Family* P )
{
	Family *s, *tmp;
	HASH_ITER( hh, P, s, tmp ){
		HASH_DEL( P, s );
	}
	return 0;
}

//Free Ytree
int freeYtree( NodeArray* P )
{
	freeFamilyHash( P->nodeArray );
	free( P->nodeArray );
	P->nodeArray = NULL;	
	P->maxn = 0;
	P->n = 0;
	return 0;
}


#include "GTree.h"
#include "uthash.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>

#define BUFFERLEN 100000
#define LINELEN 100000

Family* Root;
Family* nodeHash;
Snp2Family* snpHash;
Pos2Snp* posHash;

int Verbose = 0;
int print_usage( void );
void printf_data( Family* data );
int loadTree( char* f);
int loadSnpIndex( char* f );
int loadSampleVcf( char* f);

int main(int argc, char** argv)
{
	char treeFile[NAMELEN]="\0";
	char indexFile[NAMELEN]="\0";
	int vcfStart=0, vcfEnd=0;
	char next_option;
	const char* short_options = "hf:i:d:v";
	const struct option long_options[] = {
		{"help",0,NULL,'h'},
		{"verbose",0,NULL,'v'},
		{"ref",1,NULL,'f'},
		{"input",1,NULL,'i'},
		{"index",1,NULL,'d'},
		{NULL,0,NULL,0},
	};

	do{
		next_option = getopt_long(argc, argv, short_options, long_options, NULL);
		switch(next_option){
			case 'h':
				print_usage();
				break;
			case 'f':
				sscanf(optarg,"%s", treeFile);
				break;
			case 'd':
				sscanf(optarg,"%s", indexFile);
				break;
			case 'i':
				vcfStart = optind-1;
				while( optind < argc && argv[optind][0] != '-' )
					optind++;
				vcfEnd = optind;
				break;	
			case 'v':
				Verbose = 1;
				break;;
			case -1:
				break;
			default:
				fprintf(stderr,"Unexpected arguments %d, exit to system.\n", next_option );
				print_usage();
				break;
		}
	}while(next_option!=-1);


	if( treeFile[0] == '\0' || indexFile[0] == '\0' || vcfEnd == 0 )
		print_usage();
	
	loadTree( treeFile );
	loadSnpIndex( indexFile );
	
	if( Verbose ){
		fprintf( stderr, "\nThe tree is:\n" );
		familyDisplay( Root, printf_data, 0, '.' );
	}

	int i;
	for( i=vcfStart; i<vcfEnd; i++ )
		loadSampleVcf( argv[i] );

	freePosHash( posHash );
	freeSnpHash( snpHash );
	freeYtree( &Ytree );

	return 0;
}

int print_usage(void)
{
	fprintf(stderr,"Version\n"
			"\tAuthor:   Wang Hailong,wanghailong@firstdimension.net\n"
			"\tCompany:  FirstDimension\n"
			"\tVersion:  0.1\n"
			"\tCreated:  03/21/2017\n\n");
	
	fprintf(stderr,"Description\n"
		"\t-h|--help		Show this help\n"
		"\t-f|--ref	<file>	Tree file\n"
		"\t-d|--index	<file>	SNP index file\n"
		"\t-i|--input	<file1> [file2 file3 ...]\n"
		"\t			Vcfs\n"
		"\t-v|--verbose		Show more information\n\n");
	exit(1);
	return 0;
}

void printf_data(Family* data)
{
	fprintf(stderr,"%s\n", data->familyId);
}

int loadTree( char* filename )
{	
	FILE* f = fopen( filename, "r" );
	if( !f ){
		fprintf(stderr, "Cannot open file %s\n", filename);
		return 1;
	}

	char line[LINELEN];
	char id[NAMELEN];
	char alias[NAMELEN];
	char parent[NAMELEN];
	char snps[LINELEN];
	char *tmpstr;
	int len,i,nrow;
	Snp2Family* snp;
	Family* pNode;
	Family* cNode;
	Family* tmpNode;
	Ytree.nodeArray = (Family*)calloc( BUFFERLEN, sizeof(Family) );
	Ytree.n = 0;
	Ytree.maxn = BUFFERLEN;

	while( NULL != fgets(line, LINELEN-1, f) )
	{//Add all genotypes into hash nodeHash
		tmpNode = NULL;
		if( Ytree.n == Ytree.maxn ){
			fprintf(stderr,"Please set BUFFERLEN in makefile more bigger and make again\n");
			exit(1);
		}
			
		sscanf( line, "%s", id );
		cNode = Ytree.nodeArray + Ytree.n;
		strncpy( cNode->familyId, id, NAMELEN );

		HASH_FIND_STR( nodeHash, cNode->familyId, tmpNode );
		if( tmpNode ){
			if( Verbose )
				fprintf( stderr, "Type %s duplicated\n", cNode->familyId );
			//exit(1);
		}else
			HASH_ADD_STR( nodeHash, familyId, cNode );

		Ytree.n++;
	}

	nrow = 0;
	rewind( f );
	while( NULL != fgets(line, LINELEN-1, f) ){
		pNode = NULL;
		sscanf( line, "%s %s %s %[^\n]", id, alias, parent, snps );

		cNode = Ytree.nodeArray + nrow;
		
		if( strncmp( id, parent, NAMELEN ) == 0 || strncmp( "-", parent, NAMELEN ) == 0 || strncmp( "Root", id, NAMELEN ) == 0 ){ //ROOT
			HASH_FIND_STR( nodeHash, id, pNode );
			if( pNode ){
				Root = pNode;
				addFamily( NULL, pNode );
			}else{
				fprintf(stderr,"Cannot find genotype %s\n", id);
				exit(1);
			}
		}else{ //Children
			HASH_FIND_STR( nodeHash, parent, pNode );
			if( pNode ){
				addFamily( pNode, cNode );
			}else{
				fprintf(stderr, "Cannot find genotype %s\n", parent );
				exit(1);
			}
		}
		nrow++;

		len = strlen(snps);
		for( i=len-1; i>=0; i-- ){
			tmpstr = NULL;
			if( i == 0 && !isspace(snps[i]) )
				tmpstr = snps;
			else if( isspace(snps[i]) || snps[i]=='\0' ){
				if( !isspace(snps[i+1]) && snps[i+1]!='\0' ){//Suffix word array
					tmpstr = snps+i+1;
				}else
					snps[i] = '\0';
			}

			if( tmpstr ){
				HASH_FIND_STR( snpHash, tmpstr, snp );
				if( snp ){
					if( Verbose )
						fprintf( stderr, "Marker %s duplicated\n", snp->snp );
				}else{
					snp = (Snp2Family*)calloc( 1, sizeof(Snp2Family) ); //For struct hash key, the memery should be initialized.
					strncpy( snp->snp, tmpstr, NAMELEN );
					HASH_ADD_STR( snpHash, snp, snp );
				}
				setSnpsFamily( snp, cNode );
				snps[i]='\0';
			}
		}
	}

	fclose( f );
	return 0;
}

int str2upper( char* a ){
	int i=0;
	while( a[i] != '\0' && i<NAMELEN ){
		a[i] = toupper(a[i]);
		i++;
	}
	return 0;
}

int loadSnpIndex( char* filename )
{
	FILE* f = fopen( filename, "r" );
	if( !f ){
		fprintf(stderr, "Cannot open file %s\n", filename);
		return 1;
	}

	char snpName[NAMELEN];
	char line[LINELEN];
	char ref[NAMELEN];
	char alt[NAMELEN];
	Snp2Family* fsnp;
	Pos2Snp* psnp;
	SnpKey psnpkey;
	int pos;
	
	while( NULL != fgets(line, LINELEN-1, f) ){
		fsnp = NULL;
		psnp = NULL;
		sscanf( line, "%s %d %[^-]->%s", snpName, &pos, ref, alt );
		str2upper( alt );
		
		HASH_FIND_STR( snpHash, snpName, fsnp );
		if( fsnp ){
			memset( &psnpkey, '\0', sizeof(SnpKey) );
			psnpkey.pos = pos;
			strncpy( psnpkey.alt, alt, NAMELEN );

			HASH_FIND( hh, posHash, &psnpkey, sizeof(SnpKey), psnp );
			if( psnp ){
				if( Verbose )
					fprintf( stderr, "chr\t%d\t%s duplicated\n", psnp->snpkey.pos, psnp->snpkey.alt );
			}else{
				psnp = (Pos2Snp*)calloc( 1, sizeof(Pos2Snp) );
				psnp->snpkey = psnpkey;
				HASH_ADD( hh, posHash, snpkey, sizeof(SnpKey), psnp );
			}
			setSnpsName( psnp, fsnp );
		}
		//else
		//	fprintf(stderr, "Cannot find snp name %s in function loadSnpIndex\n", snpName);
	}
	
	fclose( f );
	return 0;
}

int loadSampleVcf( char* filename )
{
	FILE* f = fopen( filename, "r" );
	if( !f ){
		fprintf(stderr, "Cannot open file %s\n", filename);
		return 1;
	}

	int Nsnp4phy=0;
	char line[LINELEN];
	char ref[NAMELEN];
	Family* result;
	Pos2Snp* psnp;
	SnpKey snpkey;

	while( NULL != fgets( line, LINELEN-1, f ) ){
		psnp = NULL;
		if( line[0] == '#' )
			continue;
		memset( &snpkey, '\0', sizeof(SnpKey) );
		sscanf( line,"%*s %d %s %s", &(snpkey.pos), ref, snpkey.alt );
		
		//Alter 1
		str2upper( snpkey.alt );
		HASH_FIND( hh, posHash, &snpkey, sizeof(SnpKey), psnp );
		if( psnp ){
			snpsFamilyWeight1( psnp );
			Nsnp4phy++;
		}

		//Alter 2
		str2upper( ref );
		if( strcmp(ref,snpkey.alt) != 0 ){
			psnp = NULL;
			strcpy( snpkey.alt, ref );
			HASH_FIND( hh, posHash, &snpkey, sizeof(SnpKey), psnp );
			if( psnp ){
				snpsFamilyWeight1( psnp );
				Nsnp4phy++;
			}
		}
	}

	setPathWight( Root );
	setForwardPathWight( Root );
	result = getLikelyFamily( &Ytree );

	printf("%s:\t%s\t, %.2lf in %d snps surpported. Path is ", filename, result->familyId, result->pathWeight + result->forwardPathWeight, Nsnp4phy);
	showPath( result );
	printf("\n");

	fclose( f );

	if( Verbose ){
		char dotfile[]="./treeDot.gvz";
		familyDisplayDot( &Ytree, dotfile );
	}
	initFamilyWeight( &Ytree );

	return 0;
}

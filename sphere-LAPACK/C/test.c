/*
 tightly connected sphere elements (C.version)
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/types.h>
#include <sys/resource.h>
#include <sys/time.h>

#include <mpi.h>

/* static functions */
double* vectorR(int);
double** matrixR(int,int);
double dmin1(double,double);

void free_vectorR(double*);
void free_matrixR(double**);

long int* vectorR2(int);
void free_vectorR2(long int*);


/* Main */
int main(){
	double vol,area,qvol,cond0,cond,radimax;
	double hconv,t0,surf,del,delq,coef1,coef2;
	double *xc,*yc,*zc;
	double **amat;
	double *rhs;

        double *xx;
        double *amat2;

	double dx,dy,dz;
	double eps;
	double coef,resid;
        
        double t1,t2,tw;
        double dflops;

	int nx,ny,nz,n;
	int r,z,p,q,dd;
	FILE *fp;
	int i,j,k,icou,iter,ier,ic;
	int n1,n3;

        int *piv;
        int info;
        int nn,inc;
                 


/*
	+------+
	| INIT |
	+------+
*/
	if( (fp=fopen("inp1","r")) == NULL){
		fprintf(stdout,"input file cannot be opened!\n");
		return 1;
	}
		fscanf(fp,"%d %d %d",&nx,&ny,&nz);
		fscanf(fp,"%lf %lf %lf",&dx,&dy,&dz);
		fscanf(fp,"%lf %lf %lf %lf",&vol,&area,&qvol,&cond0);
		fscanf(fp,"%lf %lf %lf",&hconv,&t0,&surf);
		fscanf(fp,"%lf",&radimax);
		fscanf(fp,"%lf %lf",&coef1,&coef2);
	fclose(fp);

	n=nx*ny*nz;

 	xc=vectorR(n);
	yc=vectorR(n);
	zc=vectorR(n);

	icou= 0;
	for(k=0;k<nz;k++){
		for(j=0;j<ny;j++){
			for(i=0;i<nx;i++){
				icou++;
				xc[icou-1]=(double)i * dx;
				yc[icou-1]=(double)j * dy;
				zc[icou-1]=(double)k * dz;
			}
		}
	}

/*
	+--------+
	| MATRIX |
	+--------+
*/
	amat=matrixR(n,n);
	rhs =vectorR(n);
 
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			amat[i][j]=0.0;
		}
	}
	for(i=0;i<n;i++) rhs[i]=0.0;

	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if( j != i ){
				del=sqrt( (xc[i]-xc[j])*(xc[i]-xc[j]) 
					+ (yc[i]-yc[j])*(yc[i]-yc[j])
					+ (zc[i]-zc[j])*(zc[i]-zc[j]));
				if( del <= radimax ){
					cond=cond0/pow(10.0,dmin1(del,20.0));
					coef=cond*area/del;
					amat[i][j] =coef;
					amat[i][i]+=-coef;
				}
			}
		}
	}

	for(i=0;i<n;i++){
		delq=sqrt( xc[i]*xc[i]+yc[i]*yc[i]+zc[i]*zc[i]);
		rhs[i]=-qvol*(coef1*vol+coef2*delq*vol);
	}

	i=0;
	for(k=0;k<nz;k++){
		for(j=0;j<ny;j++){
			ic=k*nx*ny+j*nx+i;
			amat[ic][ic]+=-hconv*surf;
			rhs[ic]+=-hconv*surf*t0;
		}
	}
/*
	+-----------------------+
	| LAPACK Call           |
	+-----------------------+
*/
        amat2 = vectorR(n*n);
        k = 0;
        for(i=0;i<n;i++) {
          for(j=0;j<n;j++) {
            amat2[k] = amat[i][j];
            k++;
	  }
        }
        piv = vectorR2(n*n);

	free_matrixR(amat);
        nn = n;
        inc = 1;      

        t1 = MPI_Wtime();

	// --------------------------------------------------
        // Please write LAPACK call in here.
        // dgesv_(.........);
	dgesv_(&nn, &inc , amat2, &nn, piv, rhs , &nn , &info);

        t2 = MPI_Wtime();
 
        tw = t2-t1;
	printf("time  = %lf [sec.] \n",tw);

        dflops = (2.0/3.0)*(double)n*(double)n*(double)n;
	dflops += 7.0/2.0*(double)n*(double)n;
	dflops += 4.0/3.0*(double)n;

        dflops = dflops / tw /1.0e6;
	printf(" %6.2lf [MFLOPS] \n",dflops);

        free_vectorR(amat2);
        free_vectorR2(piv);

/*
	+--------+
	| OUTPUT |
	+--------+
*/

	n1=1;
	n3=3;
/*	MESH */
	if( (fp=fopen("sphere.fld","w")) == NULL){
		fprintf(stdout,"output file cannot be opened!\n");
		return 1;
	}
		fprintf(fp,"# AVS field file\n");
		fprintf(fp,"ndim=%d\n",n3);
		fprintf(fp,"dim1=%d\n",nx);
		fprintf(fp,"dim2=%d\n",ny);
		fprintf(fp,"dim3=%d\n",nz);
		fprintf(fp,"nspace=%d\n",n3);
		fprintf(fp,"veclen=%d\n",n1);
		fprintf(fp,"data=float\n");
		fprintf(fp,"field=uniform\n");
		fprintf(fp,"label=temperature\n");
		fprintf(fp,"variable 1 file=./spheredata filetype=ascii\n");
	fclose(fp);
/*	RESULTS */
	if( (fp=fopen("spheredata","w")) == NULL){
		fprintf(stdout,"output file cannot be opened!\n");
		return 1;
	}
		for(i=0;i<n;i++)
			fprintf(fp,"%e\n",rhs[i]);
		for(i=n-nx;i<n;i++)
			fprintf(stdout,"%d %e\n",i,rhs[i]);
	fclose(fp);


	if( (fp=fopen("spheredata_ans","r")) == NULL){
		fprintf(stdout,"input file cannot be opened!\n");
		return 1;
	}
	        xx = vectorR(n);

		for(i=0;i<n;i++)
		       fscanf(fp,"%lf\n",&xx[i]);

	fclose(fp);
        
        eps = 0.0;
	for(i=0;i<n;i++) {
          eps = eps + (rhs[i]-xx[i])*(rhs[i]-xx[i]);
        }

        eps = sqrt(eps);
        printf("err = %e \n",eps);

        
        free_vectorR(xx);

     return 0;
}
/**************************
  allocate vector
**************************/
double* vectorR(int m)
{
	double *a;
	int i;
	if ( ( a=(double * )malloc( m * sizeof( double ) ) ) == NULL ) {
		printf("Error:Memory does not enough! a in vector \n");
		exit(1);
	}
/*** zero clear ***/  
	for(i=0;i<m;i++){
		a[i]=0.0;
	}
	return a;
}
long int* vectorR2(int m)
{
	long int *a;
	int i;
	if ( ( a=(long int * )malloc( m * sizeof( long int ) ) ) == NULL ) {
		printf("Error:Memory does not enough! a in vector \n");
		exit(1);
	}
	return a;
}

/**************************
  delete vector
**************************/
void free_vectorR(double *a)
{
	free( (double *) a );
}
void free_vectorR2(long int *a)
{
	free( (long int *) a );
}

/**************************
  allocate matrix 
**************************/
double** matrixR(int m,int n)
{
	double **aa;
	int i,j;
	if ( ( aa=(double ** )malloc( m * sizeof( double * ) ) ) == NULL ) {
		printf("Error:Memory does not enough! aa in matrix \n");
		exit(1);
	}
	if ( ( aa[0]=(double * )malloc( m * n * sizeof( double ) ) ) == NULL ) {
		printf("Error:Memory does not enough! aa in matrix \n");
		exit(1);
	}
	for(i=1;i<m;i++) aa[i]=aa[i-1]+n;
/*** zero clear ***/  
	for(j=0;j<n;j++){
		for(i=0;i<m;i++){
			aa[i][j]=0.0;
		}
	}
	return aa;
}
/**************************
  delete matrix 
**************************/
void free_matrixR(double **aa)
{
	free( (double *) aa[0] );
	free( (double **) aa );
}
/**************************
  dmin1
**************************/
double dmin1(double a,double b)
{
	if( a < b ){
		return a;
	}else{
		return b;
	}
}

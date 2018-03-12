/*
 tightly connected sphere elements (C.version)
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/types.h>
#include <sys/resource.h>
#include <sys/time.h>

/* static functions */
double* vectorR(int);
double** matrixR(int,int);
double dmin1(double,double);

void free_vectorR(double*);
void free_matrixR(double**);

double my_clock()
{
  struct rusage RU;
  getrusage(RUSAGE_SELF, &RU);
  return RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec*1e-6;
}

/* Main */
int main(){
	double vol,area,qvol,cond0,cond,radimax;
	double hconv,t0,surf,del,delq,coef1,coef2;
	double *xc,*yc,*zc;
	double **amat;
	double *rhs;
	double *phi;
	double **w;
        
        double t1,t2,tw;

	int nx,ny,nz,n;
	int r,z,p,q,dd;
	FILE *fp;
	int i,j,k,icou,iter,ier,ic;
	int n1,n3;
	double dx,dy,dz;
	double eps;
	double coef,resid;
	double alpha,beta,bnrm2,dnrm2,c1;
	double rho,rho1;
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
	+---------------+
	| CG iterations |
	+---------------+
*/
	eps=1.0e-8;
	w=matrixR(n,4);
	phi=vectorR(n);

        t1=my_clock();

	for(i=0;i<n;i++){
		for(j=0;j<4;j++){
			w[i][j]=0.0;
		}
	}
	for(i=0;i<n;i++) phi[i]=0.0;

	r=0;
	z=1;
	q=1;
	p=2;
	dd=3;

	for(i=0;i<n;i++){
		w[i][dd]=1.0/amat[i][i];
	}
/*
	 {r0}= {b} - [A]{xini}
*/
	for(i=0;i<n;i++){
		w[i][r]=rhs[i];
		for(j=0;j<n;j++){
			w[i][r]+=-amat[i][j]*phi[j];
		}
	}
	bnrm2=0.0;
	for(i=0;i<n;i++){
		bnrm2+=rhs[i]*rhs[i];
	}
/**********************************************************************/
	for(iter=0;iter<n;iter++){
/* {z}= [Minv]{r} */
		for(i=0;i<n;i++){
			w[i][z]=w[i][dd]*w[i][r];
		}
/*  RHO= {r}{z} */
		rho=0.0;
		for(i=0;i<n;i++){
			rho+=w[i][r]*w[i][z];
		}

/*  {p} = {z} if      ITER=1  
	   BETA= RHO / RHO1  otherwise */

		if( iter == 0 ){
			for(i=0;i<n;i++){
				w[i][p]=w[i][z];
			}
		}else{
			beta=rho/rho1;
			for( i=0;i<n;i++){
				w[i][p]=w[i][z]+beta*w[i][p];
			}
		}
/* {q}= [A]{p} */
		for(i=0;i<n;i++){
			w[i][q]=0.0;
			for(j=0;j<n;j++){
				w[i][q]+=amat[i][j]*w[j][p];
			}
		}
/* ALPHA= RHO / {p}{q} */
		c1=0.0;
		for(i=0;i<n;i++){
			c1+=w[i][p]*w[i][q];
		}
		alpha=rho/c1;
/*	{x} = {x} + ALPHA*{p}
	{r}= {r} - ALPHA*{q} */

		for(i=0;i<n;i++){
			phi[i]+=alpha*w[i][p];
			w[i][r]+=-alpha*w[i][q];
		}

		dnrm2=0.0;
		for(i=0;i<n;i++){
			dnrm2+=w[i][r]*w[i][r];
		}

		resid=sqrt(dnrm2/bnrm2);

		fprintf(stdout,"%d %e\n",iter+1,resid);

		if( resid <= eps ) goto FINISH;

		rho1=rho;
	}
/**********************************************************************/
	ier=1;
	FINISH:

        t2=my_clock();

        tw = t2 - t1;

        printf("time  = %lf [sec.] \n",tw);

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
			fprintf(fp,"%e\n",phi[i]);
		for(i=n-nx;i<n;i++)
			fprintf(stdout,"%d %e\n",i,phi[i]);
	fclose(fp);
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
/**************************
  delete vector
**************************/
void free_vectorR(double *a)
{
	free( (double *) a );
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

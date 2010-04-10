

#include <fftw3.h>
#include <math.h>


//
//f is the NxNxN source term in row major order
//phi is the potential also in row major
//
//
void hockney(const double* f, int N, \
				double xmin, double xmax, \
				double ymin, double ymax, \
				double zmin, double zmax, \
				double* phi){

	int i,j,k,indx,bigindx;

	double h = (xmax - xmin) / N;

	int bigN = 2*N;
	//temp arrays to store the padded data
	fftw_complex* bigf  = (fftw_complex*) fftw_malloc( sizeof(fftw_complex)*bigN*bigN*bigN);

	fftw_plan forwardplan = fftw_plan_dft_3d(bigN, bigN, bigN, bigf, bigf, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan inverseplan = fftw_plan_dft_3d(bigN, bigN, bigN, bigf, bigf, FFTW_BACKWARD, FFTW_ESTIMATE);

	

	//initialize bigf
    for(i=0; i < bigN ; i++){
        for(j=0; j < bigN ; j++){
            for(k=0; k < bigN ; k++){
                bigindx = k + bigN*(j + bigN*i);
				bigf[bigindx][0] = 0.0;
				bigf[bigindx][1] = 0.0;
			}
		}
	}

	//copy f into the front-bottom-left of bigf
    //indx tracks the row major order style
    for(i=0; i < N ; i++){
        for(j=0; j < N ; j++){
            for(k=0; k < N ; k++){
                indx = k + N*(j + N*i);
                bigindx = k + bigN*(j + bigN*i);
                bigf[bigindx][0] = f[indx]; 
            }
        }
    }


	fftw_execute(forwardplan);
	fftw_execute(inverseplan);


	//copy into phi
    for(i=0; i < N ; i++){
        for(j=0; j < N ; j++){
            for(k=0; k < N ; k++){
                indx = k + N*(j + N*i);
                bigindx = k + bigN*(j + bigN*i);
                phi[indx] = (1.0/(bigN*bigN*bigN))*bigf[bigindx][0];
			}
		}
	}


	fftw_destroy_plan(forwardplan);
	fftw_destroy_plan(inverseplan); 

	fftw_free(bigf);

	return;
}


int main(){

	int N = 6;
	double xmin = -5;
	double xmax = 5;
	double ymin = -5;
	double ymax = 5;
	double zmin = -5;
	double zmax = 5;

	double* f = (double*)malloc(sizeof(double)*N*N*N);
	double* phi = (double*)malloc(sizeof(double)*N*N*N);
	
	int i,j,k,indx;


    for(i=0; i < N ; i++){
        for(j=0; j < N ; j++){
            for(k=0; k < N ; k++){
                indx = k + N*(j + N*i);
                f[indx] = indx;
			}
		}
	}

	
	hockney(f, N, xmin, xmax, ymin, ymax, zmin, zmax, phi);

    
	for(i=0; i < N ; i++){
        for(j=0; j < N ; j++){
            for(k=0; k < N ; k++){
                indx = k + N*(j + N*i);
                printf("%f, difference\n", f[indx]-phi[indx]);
			}
		}
	}
	

	free(f);
	free(phi);

	return 0;
}

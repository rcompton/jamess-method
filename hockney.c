#include <fftw3.h>
#include <math.h>

#define PI  3.1415927

//
//spatial -> discrete coords
//
//A B
//C D
//has origin at bottom left
//(like in real life)
//(i=0,j=0) is C's lower left corner
//no funny matrix indexing...
void shift_out_of_space_3d(double* a, int N){

    double* b = (double*)malloc(sizeof(double)*N*N*N);
	
    int i,j,k, indx1, indx2;

	for(i=0;i<N*N*N;i++)
		b[i] = a[i];

    //replace Af with Db
    for(i=0; i < N/2; i++){
        for(j=0; j < N/2; j++){
            for(k=0; k < N/2; k++)
            {
                indx1 = k+ N/2 + N*(j + N*(i+N/2));
                indx2 = k + N*(j+N/2 + N*(i));
                a[indx1] = b[indx2];
            }
        }
    }

    //replace Bf with Cb
    for(i=0; i < N/2; i++){
        for(j=0; j < N/2; j++){
            for(k=0; k < N/2; k++)
            {
                indx1 = k+N/2 + N*(j+N/2 + N*(i+N/2));
                indx2 = k + N*(j + N*(i));
                a[indx1] = b[indx2];
            }
        }
    }

    //replace Cf with Bb
    for(i=0; i < N/2; i++){
        for(j=0; j < N/2; j++){
            for(k=0; k < N/2; k++)
            {
                indx1 = k+ N*(j + N*(i+N/2));
                indx2 = k+N/2 + N*(j+N/2 + N*(i));
                a[indx1] = b[indx2];
            }
        }
    }

    //replace Df with Ab
    for(i=0; i < N/2; i++){
        for(j=0; j < N/2; j++){
            for(k=0; k < N/2; k++)
            {
                indx1 = k+ N*(j+N/2 + N*(i+N/2));
                indx2 = k+N/2 + N*(j + N*(i));
                a[indx1] = b[indx2];
            }
        }
    }

    //replace Ab with Df
    for(i=0; i < N/2; i++){
        for(j=0; j < N/2; j++){
            for(k=0; k < N/2; k++)
            {
                indx1 = k+N/2+ N*( + N*(i));
                indx2 = k + N*(j+N/2 + N*(i+N/2));
                a[indx1] = b[indx2];
            }
        }
    }

    //replace Bb with Cf
    for(i=0; i < N/2; i++){
        for(j=0; j < N/2; j++){
            for(k=0; k < N/2; k++)
            {
                indx1 = k + N/2 + N*(j+N/2 + N*(i));
                indx2 = k + N*(j + N*(i+N/2));
                a[indx1] = b[indx2];
            }
        }
    }

    //replace Cb with Bf
    for(i=0; i < N/2; i++){
        for(j=0; j < N/2; j++){
            for(k=0; k < N/2; k++)
            {
                indx1 = k + N*(j + N*(i));
                indx2 = k+N/2 + N*(j+N/2 + N*(i+N/2));
                a[indx1] = b[indx2];
            }
        }
    }

    //replace Db with Af
    for(i=0; i < N/2; i++){
        for(j=0; j < N/2; j++){
            for(k=0; k < N/2; k++)
            {
                indx1 = k + N*(j+N/2 + N*(i));
                indx2 = k+N/2 + N*(j + N*(i+N/2));
                a[indx1] = b[indx2];
            }
        }
    }

	free(b);

    return;
}



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

	double w1,w2,w3;

	//temp arrays to store the padded data
	fftw_complex* bigf  = (fftw_complex*) fftw_malloc( sizeof(fftw_complex)*bigN*bigN*bigN);
	double* symbol = (double*) malloc(sizeof(double)*bigN*bigN*bigN);


	//make plan
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

	for(i=0; i < bigN ; i++){
		for(j=0; j < bigN ; j++){
			for(k=0; k < bigN ; k++){
				bigindx = k + bigN*(j + bigN*i);
				w1 = 2*PI*(-N + i)/bigN;
				w2 = 2*PI*(-N + j)/bigN;
				w3 = 2*PI*(-N + k)/bigN;

				symbol[bigindx] = -1.0/(w1*w1 + w2*w2 + w3*w3);

			}
		}
	}

	shift_out_of_space_3d(symbol, bigN);

	symbol[0] = 1.0;

	//multiply in freq space
	for(i=0; i < pow(bigN, 3) ; i++){
		bigf[bigindx][0] *= symbol[i]; 
		bigf[bigindx][1] *= symbol[i]; 
	}



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

	int N = 100;
	double xmin = -5;
	double xmax = 5;
	double ymin = -5;
	double ymax = 5;
	double zmin = -5;
	double zmax = 5;

	double x,y,z,r,h;

	x = (xmax-xmin)/N;

	double* f = (double*)malloc(sizeof(double)*N*N*N);
	double* phi = (double*)malloc(sizeof(double)*N*N*N);

	int i,j,k,indx;


	for(i=0; i < N ; i++){
		for(j=0; j < N ; j++){
			for(k=0; k < N ; k++){
				x = xmin + i*h;
				y = ymin + j*h;
				z = zmin + k*h;
				r = sqrt(x*x + y*y + z*z);

				indx = k + N*(j + N*i);
				
				f[indx] = 0;
			}
		}
	}


	hockney(f, N, xmin, xmax, ymin, ymax, zmin, zmax, phi);


	double norm = 0.0;
	for(i=0; i < N ; i++){
		for(j=0; j < N ; j++){
			for(k=0; k < N ; k++){
				indx = k + N*(j + N*i);
				norm += pow(f[indx]-phi[indx],2);
			}
		}
	}
	printf("%f, difference\n", norm);


	free(f);
	free(phi);

	return 0;
}

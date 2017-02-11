#include <cstdlib>
#include <cstdio>
#include <cstring>

#define _MAX_(_X_, _Y_) ((_X_)>(_Y_) ? (_X_) : (_Y_)) 
#define _MIN_(_X_, _Y_) ((_X_)<(_Y_) ? (_X_) : (_Y_)) 

const int		NX			= 100;
const int		NY			= 100;
const double	CFL			= 0.1;
const double	XMIN		= 0.0;
const double	XMAX		= 1.0;
const double	YMIN		= 0.0;
const double	YMAX		= 1.0;
const double	HX			= (XMAX - XMIN) / NX;
const double	HY			= (YMAX - YMIN) / NY;
const double	DT			= CFL*_MIN_(HX*HX, HY*HY);
const double	TMAX		= 1.0;
const int		SAVE_STEP	= 100;


double **u, **u_old;

void init() 
{
	u = new double*[NX];
	u_old = new double*[NX];
	for (int i = 0; i < NX; i++) {
		u[i] = new double[NY];
		u_old[i] = new double[NY];

		memset(u[i], 0, sizeof(double)*NY);
	}
}

void done()
{
	for (int i = 0; i < NX; i++) {
		delete[] u[i];
		delete[] u_old[i];
	}
	delete[] u;
	delete[] u_old;
}

void bnd_cond()
{
	for (int i = 0; i < NX; i++) {
		u[i][0] = 0.0;
		u[i][NY-1] = 0.0;
	}

	for (int i = 0; i < NY; i++) {
		u[0][i] = 1.0;
		u[NX-1][i] = 0.0;
	}
}

void swap()
{
	for (int i = 0 ; i < NX; i++) {
		memcpy(u_old[i], u[i], sizeof(double)*NY);
	}
}

void time_step()
{
	for (int i = 1; i < NX-1; i++) {
		for (int j = 1; j < NY-1; j++) {
			u[i][j] = u_old[i][j] + DT*(
				(u_old[i-1][j] - 2.0*u_old[i][j] + u_old[i+1][j])/(HX*HX) +
				(u_old[i][j-1] - 2.0*u_old[i][j] + u_old[i][j+1])/(HY*HY)
			);
		}
	}
}

void save_vtk(int step)
{
	FILE *fp;
	char f_name[128];

	sprintf(f_name, "result_%010d.vtk", step);
	fp = fopen(f_name, "w");

	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "Temperature\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET STRUCTURED_GRID\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", NX+1, NY+1, 1);
	fprintf(fp, "POINTS %d float\n", (NX+1)*(NY+1));
	for (int j = 0 ; j <= NY; j++) {
		for (int i = 0 ; i <= NX; i++) {
			fprintf(fp, "%f %f 0.0\n", XMIN+HX*i, YMIN+j*HY);
		}
	}
	fprintf(fp, "CELL_DATA %d\n", NX*NY);
	fprintf(fp, "SCALARS Temperature float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			fprintf(fp, "%f\n", u[i][j]);
		}
	}
	
	fclose(fp);

	printf("File '%s' is saved.\n", f_name);
}

int main(int argc, char** argv)
{
	init();
	bnd_cond();
	swap();
	double t = 0.0;
	int step = 0;

	while (t < TMAX) {
		t += DT;
		step++;
		time_step();
		bnd_cond();
		if (step%SAVE_STEP == 0) {
			save_vtk(step);
		}
		swap();
	}
	done();
}

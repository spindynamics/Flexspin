#include "objects.hpp"

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct{
        	double x;
        	double y;
 		    double z;
} point;

typedef struct{
       	double x0;
       	double y0;
		double z0;
		double Lx;
		double Ly;
		double sigma;
} surface_z;

typedef struct{
       	double x0;
       	double y0;
		double z0;
		double Ly;
		double Lz;
		double sigma;
} surface_x;
typedef struct{
       	double x0;
       	double y0;
		double z0;
		double Lz;
		double Lx;
		double sigma;
} surface_y;

double IparaXX(surface_x S1,surface_x S2);
double IparaYY(surface_y S1,surface_y S2);
double IparaZZ(surface_z S1,surface_z S2);

double IperpXY(surface_x S1,surface_y S2);
double IperpZX(surface_z S1,surface_x S2);
double IperpYZ(surface_y S1,surface_z S2);

double F220(double x, double y, double z);
double F022(double x, double y, double z);
double F202(double x, double y, double z);

double F211(double x, double y, double z);
double F121(double x, double y, double z);
double F112(double x, double y, double z);

tensor coeffdemag(const Layer& S1, const Layer& S2)
{
	double     x01, y01, z01, x02, y02, z02;
	double     Lx1, Ly1, Lz1, Lx2, Ly2, Lz2 ;

	surface_x  Sx1, Sx2;
	surface_y  Sy1, Sy2;
	surface_z  Sz1, Sz2;

	double     Ndipxx, Ndipyy,Ndipzz;
	double     Ndipxy, Ndipxz, Ndipyx, Ndipzx, Ndipyz, Ndipzy;

	tensor N;

	Lx1=S1.size.x;   Ly1=S1.size.y;   Lz1=S1.size.z;
	x01=S1.center.x; y01=S1.center.y; z01=S1.center.z;

	Lx2=S2.size.x;   Ly2=S2.size.y;   Lz2=S2.size.z;
	x02=S2.center.x; y02=S2.center.y; z02=S2.center.z;
	
	//std::cout << "Lx1 = " << Lx1 << ", Ly1 = " << Ly1 << ", Lz1 = " << Lz1 << std::endl;
	//std::cout << "x01 = " << x01 << ", y01 = " << y01 << ", z01 = " << z01 << std::endl;
	//std::cout << "Lx2 = " << Lx2 << ", Ly2 = " << Ly2 << ", Lz2 = " << Lz2 << std::endl;
	//std::cout << "x02 = " << x02 << ", y02 = " << y02 << ", z02 = " << z02 << std::endl;



	Ndipxx=0.;
	/*--- interaction surfaces left -left------*/
	Sx1.x0=x01 - Lx1/2.; Sx1.y0=y01; Sx1.z0=z01; Sx1.Lz=Lz1; Sx1.Ly=Ly1; Sx1.sigma=-1;
	Sx2.x0=x02 - Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Lz=Lz2; Sx2.Ly=Ly2; Sx2.sigma=-1;
	Ndipxx += Sx1.sigma*Sx2.sigma*IparaXX(Sx1, Sx2)/(4.*M_PI);
	/*--- interaction surfaces left-right------*/
	Sx1.x0=x01 - Lx1/2.; Sx1.y0=y01; Sx1.z0=z01; Sx1.Lz=Lz1; Sx1.Ly=Ly1; Sx1.sigma= -1;
	Sx2.x0=x02 + Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Lz=Lz2; Sx2.Ly=Ly2; Sx2.sigma=  1;
	Ndipxx += Sx1.sigma*Sx2.sigma*IparaXX(Sx1, Sx2)/(4.*M_PI);
	/*--- interaction surfaces right-left------*/
	Sx1.x0=x01 + Lx1/2.; Sx1.y0=y01; Sx1.z0=z01; Sx1.Lz=Lz1; Sx1.Ly=Ly1; Sx1.sigma=  1;
	Sx2.x0=x02 - Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Lz=Lz2; Sx2.Ly=Ly2; Sx2.sigma= -1;
	Ndipxx += Sx1.sigma*Sx2.sigma*IparaXX(Sx1, Sx2)/(4.*M_PI);
	/*--- interaction surfaces up - up------*/
	Sx1.x0=x01 + Lx1/2.; Sx1.y0=y01; Sx1.z0=z01; Sx1.Lz=Lz1; Sx1.Ly=Ly1; Sx1.sigma= 1;
	Sx2.x0=x02 + Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Lz=Lz2; Sx2.Ly=Ly2; Sx2.sigma= 1;
	Ndipxx += Sx1.sigma*Sx2.sigma*IparaXX(Sx1, Sx2)/(4.*M_PI);


    Ndipyy=0.;
	/*--- interaction surfaces back -back------*/
	Sy1.x0=x01; Sy1.y0=y01- Ly1/2.; Sy1.z0=z01; Sy1.Lz=Lz1; Sy1.Lx=Lx1; Sy1.sigma=-1;
	Sy2.x0=x02; Sy2.y0=y02- Ly2/2.; Sy2.z0=z02; Sy2.Lz=Lz2; Sy2.Lx=Lx2; Sy2.sigma=-1;
	Ndipyy += Sy1.sigma*Sy2.sigma*IparaYY(Sy1, Sy2)/(4.*M_PI);
	/*--- interaction surfaces back-front------*/
	Sy1.x0=x01; Sy1.y0=y01- Ly1/2.; Sy1.z0=z01; Sy1.Lz=Lz1; Sy1.Lx=Lx1; Sy1.sigma=-1;
	Sy2.x0=x02; Sy2.y0=y02+ Ly2/2.; Sy2.z0=z02; Sy2.Lz=Lz2; Sy2.Lx=Lx2; Sy2.sigma= 1;
	Ndipyy += Sy1.sigma*Sy2.sigma*IparaYY(Sy1, Sy2)/(4.*M_PI);
	/*--- interaction surfaces front -back------*/
	Sy1.x0=x01; Sy1.y0=y01+ Ly1/2.; Sy1.z0=z01; Sy1.Lz=Lz1; Sy1.Lx=Lx1; Sy1.sigma= 1;
	Sy2.x0=x02; Sy2.y0=y02- Ly2/2.; Sy2.z0=z02; Sy2.Lz=Lz2; Sy2.Lx=Lx2; Sy2.sigma=-1;
	Ndipyy += Sy1.sigma*Sy2.sigma*IparaYY(Sy1, Sy2)/(4.*M_PI);
	/*--- interaction surfaces front-front------*/
	Sy1.x0=x01; Sy1.y0=y01+ Ly1/2.; Sy1.z0=z01; Sy1.Lz=Lz1; Sy1.Lx=Lx1; Sy1.sigma=1;
	Sy2.x0=x02; Sy2.y0=y02+ Ly2/2.; Sy2.z0=z02; Sy2.Lz=Lz2; Sy2.Lx=Lx2; Sy2.sigma=1;
	Ndipyy += Sy1.sigma*Sy2.sigma*IparaYY(Sy1, Sy2)/(4.*M_PI);


    Ndipzz=0.;
	/*--- interaction surfaces down -down------*/
	Sz1.x0=x01; Sz1.y0=y01; Sz1.z0=z01 - Lz1/2.; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma=-1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 - Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma=-1;
	Ndipzz += Sz1.sigma*Sz2.sigma*IparaZZ(Sz1, Sz2)/(4.*M_PI);
	/*--- interaction surfaces down - up------*/
	Sz1.x0=x01; Sz1.y0=y01; Sz1.z0=z01 - Lz1/2.; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma=-1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 + Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma=+1;
	Ndipzz += Sz1.sigma*Sz2.sigma*IparaZZ(Sz1, Sz2)/(4.*M_PI);
	/*--- interaction surfaces up - down------*/
	Sz1.x0=x01; Sz1.y0=y01; Sz1.z0=z01 + Lz1/2.; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma= 1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 - Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma=-1;
	Ndipzz += Sz1.sigma*Sz2.sigma*IparaZZ(Sz1, Sz2)/(4.*M_PI);
	/*--- interaction surfaces up - up------*/
	Sz1.x0=x01; Sz1.y0=y01; Sz1.z0=z01 + Lz1/2.; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma= 1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 + Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma= 1;
	Ndipzz += Sz1.sigma*Sz2.sigma*IparaZZ(Sz1, Sz2)/(4.*M_PI);
    //printf("\n");

    Ndipxy=0.;
    /*--- interaction surfaces left-back------*/
	Sx1.x0=x01 - Lx1/2.;Sx1.y0=y01; Sx1.z0=z01;  Sx1.Ly=Ly1; Sx1.Lz=Lz1; Sx1.sigma=-1;
	Sy2.x0=x02; Sy2.y0=y02 - Ly2/2.; Sy2.z0=z02; Sy2.Lx=Lx2; Sy2.Lz=Lz2; Sy2.sigma=-1;
	Ndipxy += Sx1.sigma*Sy2.sigma*IperpXY(Sx1, Sy2)/(4.*M_PI);
    /*--- interaction surfaces left front------*/
	Sx1.x0=x01 - Lx1/2.;Sx1.y0=y01; Sx1.z0=z01;  Sx1.Ly=Ly1; Sx1.Lz=Lz1; Sx1.sigma=-1;
	Sy2.x0=x02; Sy2.y0=y02 + Ly2/2.; Sy2.z0=z02; Sy2.Lx=Lx2; Sy2.Lz=Lz2; Sy2.sigma= 1;
	Ndipxy += Sx1.sigma*Sy2.sigma*IperpXY(Sx1, Sy2)/(4.*M_PI);
    /*--- interaction surfaces right back-------*/
	Sx1.x0=x01 + Lx1/2.;Sx1.y0=y01; Sx1.z0=z01;  Sx1.Ly=Ly1; Sx1.Lz=Lz1; Sx1.sigma= 1;
	Sy2.x0=x02; Sy2.y0=y02 - Ly2/2.; Sy2.z0=z02; Sy2.Lx=Lx2; Sy2.Lz=Lz2; Sy2.sigma=-1;
	Ndipxy += Sx1.sigma*Sy2.sigma*IperpXY(Sx1, Sy2)/(4.*M_PI);
    /*--- interaction surfaces right front------*/
	Sx1.x0=x01 + Lx1/2.;Sx1.y0=y01; Sx1.z0=z01;  Sx1.Ly=Ly1; Sx1.Lz=Lz1; Sx1.sigma= 1;
	Sy2.x0=x02; Sy2.y0=y02 + Ly2/2.; Sy2.z0=z02; Sy2.Lx=Lx2; Sy2.Lz=Lz2; Sy2.sigma= 1;
	Ndipxy += Sx1.sigma*Sy2.sigma*IperpXY(Sx1, Sy2)/(4.*M_PI);

	//printf("Ndipxy= %lf\t", Ndipxy/(Lx1*Ly1*Lz1));
	//printf("Ndipxy= %lf\n", Ndipxy/(Lx2*Ly2*Lz2));

    Ndipyx=0.;
    /*--- interaction surfaces back - left------*/
	Sy1.x0=x01; Sy1.y0=y01 - Ly1/2.;Sy1.z0=z01;  Sy1.Lx=Lx1; Sy1.Lz=Lz1; Sy1.sigma=-1;
	Sx2.x0=x02 - Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Ly=Ly2; Sx2.Lz=Lz2; Sx2.sigma=-1;
	Ndipyx += Sy1.sigma*Sx2.sigma*IperpXY(Sx2, Sy1)/(4.*M_PI);
    /*--- interaction surfaces  back - right------*/
	Sy1.x0=x01; Sy1.y0=y01 - Ly1/2.;Sy1.z0=z01;  Sy1.Lx=Lx1; Sy1.Lz=Lz1; Sy1.sigma=-1;
	Sx2.x0=x02 + Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Ly=Ly2; Sx2.Lz=Lz2; Sx2.sigma= 1;
	Ndipyx += Sy1.sigma*Sx2.sigma*IperpXY(Sx2, Sy1)/(4.*M_PI);
    /*--- interaction surfaces front left-------*/
	Sy1.x0=x01; Sy1.y0=y01 + Ly1/2.;Sy1.z0=z01;  Sy1.Lx=Lx1; Sy1.Lz=Lz1; Sy1.sigma= 1;
	Sx2.x0=x02 - Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Ly=Ly2; Sx2.Lz=Lz2; Sx2.sigma=-1;
	Ndipyx += Sy1.sigma*Sx2.sigma*IperpXY(Sx2, Sy1)/(4.*M_PI);
	/*--- interaction surfaces front right------*/
	Sy1.x0=x01; Sy1.y0=y01 + Ly1/2.;Sy1.z0=z01;  Sy1.Lx=Lx1; Sy1.Lz=Lz1; Sy1.sigma= 1;
	Sx2.x0=x02 + Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Ly=Ly2; Sx2.Lz=Lz2; Sx2.sigma= 1;
	Ndipyx += Sy1.sigma*Sx2.sigma*IperpXY(Sx2, Sy1)/(4.*M_PI);

	//printf("Ndipyx= %lf\t", Ndipyx/(Lx1*Ly1*Lz1));
	//printf("Ndipyx= %lf\n", Ndipyx/(Lx2*Ly2*Lz2));

    //printf("\n");
    Ndipxz=0.;
    /*--- interaction surfaces left- bottom------*/
	Sx1.x0=x01 - Lx1/2.; Sx1.y0=y01; Sx1.z0=z01; Sx1.Lz=Lz1; Sx1.Ly=Ly1; Sx1.sigma=-1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 - Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma=-1;
	Ndipxz += Sz2.sigma*Sx1.sigma*IperpZX(Sz2, Sx1)/(4.*M_PI);
    /*--- interaction surfaces right- bottom------*/
	Sx1.x0=x01 + Lx1/2.; Sx1.y0=y01; Sx1.z0=z01; Sx1.Lz=Lz1; Sx1.Ly=Ly1; Sx1.sigma= 1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 - Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma=-1;
	Ndipxz += Sz2.sigma*Sx1.sigma*IperpZX(Sz2, Sx1)/(4.*M_PI);
    /*--- interaction surfaces left- top------*/
	Sx1.x0=x01 - Lx1/2.; Sx1.y0=y01; Sx1.z0=z01; Sx1.Lz=Lz1; Sx1.Ly=Ly1; Sx1.sigma=-1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 + Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma= 1;
	Ndipxz += Sz2.sigma*Sx1.sigma*IperpZX(Sz2, Sx1)/(4.*M_PI);
    /*--- interaction surfaces right -top------*/
	Sx1.x0=x01 + Lx1/2.; Sx1.y0=y01; Sx1.z0=z01; Sx1.Lz=Lz1; Sx1.Ly=Ly1; Sx1.sigma= 1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 + Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma= 1;
	Ndipxz += Sz2.sigma*Sx1.sigma*IperpZX(Sz2, Sx1)/(4.*M_PI);

	//printf("Ndipxz= %lf\t", Ndipxz/(Lx1*Ly1*Lz1));
	//printf("Ndipxz= %lf\n", Ndipxz/(Lx2*Ly2*Lz2));

    Ndipzx=0.;
    /*--- interaction surfaces bottom- left------*/
	Sz1.x0=x01; Sz1.y0=y01; Sz1.z0=z01 - Lz1/2.; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma=-1;
	Sx2.x0=x02 - Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Lz=Lz2; Sx2.Ly=Ly2; Sx2.sigma=-1;
	Ndipzx += Sz1.sigma*Sx2.sigma*IperpZX(Sz1, Sx2)/(4.*M_PI);
    /*--- interaction surfaces bottom- right------*/
	Sz1.x0=x01; Sz1.y0=y01; Sz1.z0=z01 - Lz1/2.; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma=-1;
	Sx2.x0=x02 + Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Lz=Lz2; Sx2.Ly=Ly2; Sx2.sigma=+1;
	Ndipzx += Sz1.sigma*Sx2.sigma*IperpZX(Sz1, Sx2)/(4.*M_PI);
    /*--- interaction surfaces top- left------*/
	Sz1.x0=x01; Sz1.y0=y01; Sz1.z0=z01 + Lz1/2.; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma= 1;
	Sx2.x0=x02 - Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Lz=Lz2; Sx2.Ly=Ly2; Sx2.sigma=-1;
	Ndipzx += Sz1.sigma*Sx2.sigma*IperpZX(Sz1, Sx2)/(4.*M_PI);
    /*--- interaction surfaces top- right------*/
	Sz1.x0=x01; Sz1.y0=y01; Sz1.z0=z01 + Lz1/2.; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma= 1;
	Sx2.x0=x02 + Lx2/2.; Sx2.y0=y02; Sx2.z0=z02; Sx2.Lz=Lz2; Sx2.Ly=Ly2; Sx2.sigma= 1;
	Ndipzx += Sz1.sigma*Sx2.sigma*IperpZX(Sz1, Sx2)/(4.*M_PI);

	//printf("Ndipzx= %lf\t", Ndipzx/(Lx1*Ly1*Lz1));
	//printf("Ndipzx= %lf\n", Ndipzx/(Lx2*Ly2*Lz2));
    //printf("\n");

    Ndipyz=0.;
    /*--- interaction surfaces back - bottom------*/
	Sy1.x0=x01; Sy1.y0=y01 - Ly1/2.;Sy1.z0=z01;  Sy1.Lx=Lx1; Sy1.Lz=Lz1; Sy1.sigma=-1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 - Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma=-1;
	Ndipyz += Sy1.sigma*Sz2.sigma*IperpYZ(Sy1, Sz2)/(4.*M_PI);
    /*--- interaction surfaces  back - top------*/
	Sy1.x0=x01; Sy1.y0=y01 - Ly1/2.;Sy1.z0=z01;  Sy1.Lx=Lx1; Sy1.Lz=Lz1; Sy1.sigma=-1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 + Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma= 1;
	Ndipyz += Sy1.sigma*Sz2.sigma*IperpYZ(Sy1, Sz2)/(4.*M_PI);
	/*--- interaction surfaces front bottom-------*/
	Sy1.x0=x01; Sy1.y0=y01 + Ly1/2.;Sy1.z0=z01;  Sy1.Lx=Lx1; Sy1.Lz=Lz1; Sy1.sigma= 1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 - Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma=-1;
	Ndipyz += Sy1.sigma*Sz2.sigma*IperpYZ(Sy1, Sz2)/(4.*M_PI);
	/*--- interaction surfaces front top------*/
	Sy1.x0=x01; Sy1.y0=y01 + Ly1/2.;Sy1.z0=z01;  Sy1.Lx=Lx1; Sy1.Lz=Lz1; Sy1.sigma= 1;
	Sz2.x0=x02; Sz2.y0=y02; Sz2.z0=z02 + Lz2/2.; Sz2.Lx=Lx2; Sz2.Ly=Ly2; Sz2.sigma=1;
	Ndipyz += Sy1.sigma*Sz2.sigma*IperpYZ(Sy1, Sz2)/(4.*M_PI);

	//printf("Ndipyz= %lf\t", Ndipyz/(Lx1*Ly1*Lz1));
	//printf("Ndipyz= %lf\n", Ndipyz/(Lx2*Ly2*Lz2));


    Ndipzy=0.;
    /*--- interaction surfaces bottom - back------*/
	Sz1.x0=x01; Sz1.y0=y01;         Sz1.z0=z01- Lz1/2; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma=-1;
	Sy2.x0=x02; Sy2.y0=y02 - Ly2/2.;Sy2.z0=z02;        Sy2.Lx=Lx2; Sy2.Lz=Lz2; Sy2.sigma=-1;
	Ndipzy += Sz1.sigma*Sy2.sigma*IperpYZ(Sy2,Sz1)/(4.*M_PI);
    /*--- interaction surfaces  bottom - front------*/
	Sz1.x0=x01; Sz1.y0=y01;         Sz1.z0=z01- Lz1/2; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma=-1;
	Sy2.x0=x02; Sy2.y0=y02 + Ly2/2.;Sy2.z0=z02;        Sy2.Lx=Lx2; Sy2.Lz=Lz2; Sy2.sigma= 1;
	Ndipzy += Sz1.sigma*Sy2.sigma*IperpYZ(Sy2,Sz1)/(4.*M_PI);
	/*--- interaction surfaces top - back-------*/
	Sz1.x0=x01; Sz1.y0=y01;         Sz1.z0=z01+ Lz1/2; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma= 1;
	Sy2.x0=x02; Sy2.y0=y02 - Ly2/2.;Sy2.z0=z02;        Sy2.Lx=Lx2; Sy2.Lz=Lz2; Sy2.sigma=-1;
	Ndipzy += Sz1.sigma*Sy2.sigma*IperpYZ(Sy2,Sz1)/(4.*M_PI);
	/*--- interaction surfaces top -front------*/
	Sz1.x0=x01; Sz1.y0=y01;         Sz1.z0=z01+ Lz1/2; Sz1.Lx=Lx1; Sz1.Ly=Ly1; Sz1.sigma= 1;
	Sy2.x0=x02; Sy2.y0=y02 + Ly2/2.;Sy2.z0=z02;        Sy2.Lx=Lx2; Sy2.Lz=Lz2; Sy2.sigma= 1;
	Ndipzy += Sz1.sigma*Sy2.sigma*IperpYZ(Sy2,Sz1)/(4.*M_PI);

	//printf("Ndipzy= %lf\t", Ndipzy/(Lx1*Ly1*Lz1));
	//printf("Ndipzy= %lf\n", Ndipzy/(Lx2*Ly2*Lz2));
	//printf("\n");

    //printf("Ndipxx= %lf\t", Ndipxx/(Lx1*Ly1*Lz1));
	//printf("Ndipxx= %lf\n", Ndipxx/(Lx2*Ly2*Lz2));
	//printf("Ndipyy= %lf\t", Ndipyy/(Lx1*Ly1*Lz1));
	//printf("Ndipyy= %lf\n", Ndipyy/(Lx2*Ly2*Lz2));
	//printf("Ndipzz= %lf\t", Ndipzz/(Lx1*Ly1*Lz1));
	//printf("Ndipzz= %lf\n", Ndipzz/(Lx2*Ly2*Lz2));


	const double volume = S1.size.x * S1.size.y * S1.size.z;

	N.xx=Ndipxx/volume;    N.xy=Ndipxy/volume;    N.xz=Ndipxz/volume;
	N.yx=Ndipyx/volume;    N.yy=Ndipyy/volume;    N.yz=Ndipyz/volume;
	N.zx=Ndipzx/volume;    N.zy=Ndipzy/volume;    N.zz=Ndipzz/volume;

		
	//std::cout << "Ndipxx = " << Ndipxx << ", Ndipxy = " << Ndipxy << ", Ndipxz = " << Ndipxz << std::endl;
	//std::cout << "Ndipyx = " << Ndipyx << ", Ndipyy = " << Ndipyy << ", Ndipyz = " << Ndipyz << std::endl;
	//std::cout << "Ndipzx = " << Ndipzx << ", Ndipzy = " << Ndipzy << ", Ndipzz = " << Ndipzz << std::endl;
	//std::cout << "Volume = " << S1.volume << std::endl;
	
    return (N);
}

/*******************************************************************/
double F220(double x, double y, double z)
{
	double s, R=sqrt(x*x+y*y+z*z);

 	if ( (x!=0.) && (y!=0.) && (z!=0.) )
 	{
		s=     0.5*x*(y*y-z*z)*0.5*log((R+x)/(R-x));
        	s= s + 0.5*y*(x*x-z*z)*0.5*log((R+y)/(R-y));
        	s= s - x*y*z*atan(x*y/(z*R));
        	s= s + 1/6.0*R*(2.0*z*z-x*x-y*y);
	}
	if ( (x==0.) && (y!=0.) && (z!=0.) )
 	{
        	s=  -0.5*y*z*z*0.5*log((R+y)/(R-y));
        	s= s + 1/6.0*R*(2.0*z*z-y*y);
	}
 	if ( (x!=0.) && (y==0.) && (z!=0.) )
 	{
		s=   -0.5*x*z*z*0.5*log((R+x)/(R-x));
        	s= s + 1/6.0*R*(2.0*z*z-x*x);
	}
 	if ( (x!=0.) && (y!=0.) && (z==0.) )
 	{
		s=     0.5*x*y*y*0.5*log((R+x)/(R-x));
        	s= s + 0.5*y*x*x*0.5*log((R+y)/(R-y));
        	s= s + 1/6.0*R*(-x*x-y*y);
	}
	if ( (x==0.) && (y==0.) && (z!=0.) )    s=  1/6.0*R*(2.0*z*z);

	if ( (x==0.) && (y!=0.) && (z==0.) )    s= - 1/6.0*R*y*y;

	if ( (x!=0.) && (y==0.) && (z==0.) )    s= - 1/6.0*R*x*x;

	if ( (x==0.) && (y==0.) && (z==0.)  )   s=    0.;


	return (s);

}
/*******************************************************************/
double F022(double x, double y, double z)
{
	double s=F220(z,x,y);
	return (s);
}
/*******************************************************************/
double F202(double x, double y, double z)
{
	double s=F220(y,z,x);
	return (s);
}
/*******************************************************************/
double F121(double x, double y, double z)
{
        double s, R=sqrt(x*x+y*y+z*z);

	if ( (x!=0.) && (y!=0.) && (z!=0.) )
 	{
		s=x*y*z*0.5*log((R+y)/(R-y));
        	s=s+0.5*x*(y*y-1/3.0*x*x)*0.5*log((R+z)/(R-z));
        	s=s+0.5*z*(y*y-1/3.0*z*z)*0.5*log((R+x)/(R-x));
        	s=s-1.0/6.0*y*y*y*atan(x*z/(y*R));
		s=s-0.5*y*(x*x*atan(y*z/(x*R))+z*z*atan(x*y/(z*R)));
		s=s-x*z*R/3.0;
	}

	if ( (x!=0.) && (y==0.) && (z!=0.) )
 	{

		s=  0.5*x*(y*y-1/3.0*x*x)*0.5*log((R+z)/(R-z));
		s=s+0.5*z*(y*y-1/3.0*z*z)*0.5*log((R+x)/(R-x));
		s=s-x*z*R/3.0;
	}

        if ( (x!=0.) && (y!=0.) && (z==0.) )	s=0.;
        if ( (x==0.) && (y!=0.) && (z!=0.) )	s=0.;

        if ( (x==0.) && (y==0.) && (z!=0.) )	s=0.;
        if ( (x==0.) && (y!=0.) && (z==0.) )	s=0.;
        if ( (x!=0.) && (y==0.) && (z==0.) )	s=0.;

        if ( (x==0.) && (y==0.) && (z==0.) )	s=0.;


        return (s);
}
/***********************************************************/
double F112(double x, double y, double z)
{
    	double s = F121(y,z,x);

	return (s);
}
/***********************************************************/
double F211(double x, double y, double z)
{
    	double s = F121(z,x,y);
	return (s);
}
/****************************************************************/
double IparaXX(surface_x S1,surface_x S2)
{
	double aux;
	double  z1=S1.z0 - 0.5*S1.Lz, z2=S1.z0 + 0.5*S1.Lz;
	double  y1=S1.y0 - 0.5*S1.Ly, y2=S1.y0 + 0.5*S1.Ly;

	double  zp1=S2.z0-0.5*S2.Lz, zp2=S2.z0+0.5*S2.Lz;
	double  yp1=S2.y0-0.5*S2.Ly, yp2=S2.y0+0.5*S2.Ly;

	double  x=S1.x0, xp=S2.x0;


	aux  =  F220(y2-yp1, z2-zp1, x-xp) -  F220(y2-yp2, z2-zp1, x-xp);
	aux += -F220(y2-yp1, z2-zp2, x-xp) +  F220(y2-yp2, z2-zp2, x-xp);

	aux += -F220(y2-yp1, z1-zp1, x-xp) +  F220(y2-yp2, z1-zp1, x-xp);
	aux +=  F220(y2-yp1, z1-zp2, x-xp) -  F220(y2-yp2, z1-zp2, x-xp);

	aux += -F220(y1-yp1, z2-zp1, x-xp) +  F220(y1-yp2, z2-zp1, x-xp);
	aux +=  F220(y1-yp1, z2-zp2, x-xp) -  F220(y1-yp2, z2-zp2, x-xp);

	aux +=  F220(y1-yp1, z1-zp1, x-xp) -  F220(y1-yp2, z1-zp1, x-xp);
	aux += -F220(y1-yp1, z1-zp2, x-xp) +  F220(y1-yp2, z1-zp2, x-xp);

	return (aux);
}
/****************************************************************/
double IparaYY(surface_y S1,surface_y S2)
{
	double aux;
	double  x1=S1.x0 - 0.5*S1.Lx, x2=S1.x0 + 0.5*S1.Lx;
	double  z1=S1.z0 - 0.5*S1.Lz, z2=S1.z0 + 0.5*S1.Lz;

	double  xp1=S2.x0-0.5*S2.Lx, xp2=S2.x0+0.5*S2.Lx;
	double  zp1=S2.z0-0.5*S2.Lz, zp2=S2.z0+0.5*S2.Lz;

	double  y=S1.y0, yp=S2.y0;


	aux  =  F220(z2-zp1, x2-xp1, y-yp) -  F220(z2-zp2, x2-xp1, y-yp);
	aux += -F220(z2-zp1, x2-xp2, y-yp) +  F220(z2-zp2, x2-xp2, y-yp);

	aux += -F220(z2-zp1, x1-xp1, y-yp) +  F220(z2-zp2, x1-xp1, y-yp);
	aux +=  F220(z2-zp1, x1-xp2, y-yp) -  F220(z2-zp2, x1-xp2, y-yp);

	aux += -F220(z1-zp1, x2-xp1, y-yp) +  F220(z1-zp2, x2-xp1, y-yp);
	aux +=  F220(z1-zp1, x2-xp2, y-yp) -  F220(z1-zp2, x2-xp2, y-yp);

	aux +=  F220(z1-zp1, x1-xp1, y-yp) -  F220(z1-zp2, x1-xp1, y-yp);
	aux += -F220(z1-zp1, x1-xp2, y-yp) +  F220(z1-zp2, x1-xp2, y-yp);

	return (aux);
}
/****************************************************************/
double IparaZZ(surface_z S1,surface_z S2)
{
	double aux;
	double  x1=S1.x0 - 0.5*S1.Lx, x2=S1.x0 + 0.5*S1.Lx;
	double  y1=S1.y0 - 0.5*S1.Ly, y2=S1.y0 + 0.5*S1.Ly;

	double  xp1=S2.x0-0.5*S2.Lx, xp2=S2.x0+0.5*S2.Lx;
	double  yp1=S2.y0-0.5*S2.Ly, yp2=S2.y0+0.5*S2.Ly;

	double  z=S1.z0, zp=S2.z0;


	aux  =  F220(x2-xp1, y2-yp1, z-zp) -  F220(x2-xp2, y2-yp1, z-zp);
	aux += -F220(x2-xp1, y2-yp2, z-zp) +  F220(x2-xp2, y2-yp2, z-zp);

	aux += -F220(x2-xp1, y1-yp1, z-zp) +  F220(x2-xp2, y1-yp1, z-zp);
	aux +=  F220(x2-xp1, y1-yp2, z-zp) -  F220(x2-xp2, y1-yp2, z-zp);

	aux += -F220(x1-xp1, y2-yp1, z-zp) +  F220(x1-xp2, y2-yp1, z-zp);
	aux +=  F220(x1-xp1, y2-yp2, z-zp) -  F220(x1-xp2, y2-yp2, z-zp);

	aux +=  F220(x1-xp1, y1-yp1, z-zp) -  F220(x1-xp2, y1-yp1, z-zp);
	aux += -F220(x1-xp1, y1-yp2, z-zp) +  F220(x1-xp2, y1-yp2, z-zp);

	return (aux);
}

/****************************************************************/
double IperpXY(surface_x S1,surface_y S2)
{
	double aux;

	double  x =S1.x0;

	double  y1=S1.y0 - 0.5*S1.Ly,
	        y2=S1.y0 + 0.5*S1.Ly;

	double  z1=S1.z0 - 0.5*S1.Lz,
	        z2=S1.z0 + 0.5*S1.Lz;

	double  xp1=S2.x0 - 0.5*S2.Lx,
	        xp2=S2.x0 + 0.5*S2.Lx;

	double  yp =S2.y0;

	double  zp1=S2.z0 - 0.5*S2.Lz,
	        zp2=S2.z0 + 0.5*S2.Lz;

	aux  =   F112(x-xp1,y2-yp,z2-zp1) - F112(x-xp1,y2-yp,z2-zp2);
	aux += - F112(x-xp2,y2-yp,z2-zp1) + F112(x-xp2,y2-yp,z2-zp2);

	aux += - F112(x-xp1,y1-yp,z2-zp1) + F112(x-xp1,y1-yp,z2-zp2);
	aux += + F112(x-xp2,y1-yp,z2-zp1) - F112(x-xp2,y1-yp,z2-zp2);

	aux += - F112(x-xp1,y2-yp,z1-zp1) + F112(x-xp1,y2-yp,z1-zp2);
	aux += + F112(x-xp2,y2-yp,z1-zp1) - F112(x-xp2,y2-yp,z1-zp2);

	aux += + F112(x-xp1,y1-yp,z1-zp1) - F112(x-xp1,y1-yp,z1-zp2);
	aux += - F112(x-xp2,y1-yp,z1-zp1) + F112(x-xp2,y1-yp,z1-zp2);

	//printf("XY: %lf\n",aux/1.e-27);

	return (aux);
}
/****************************************************************/
double IperpYZ(surface_y S1,surface_z S2)
{
	double aux;

	double  x1=S1.x0 - 0.5*S1.Lx, x2=S1.x0 + 0.5*S1.Lx;
	double  y =S1.y0;
	double  z1=S1.z0 - 0.5*S1.Lz, z2=S1.z0 + 0.5*S1.Lz;

	double  xp1=S2.x0 - 0.5*S2.Lx, xp2=S2.x0 + 0.5*S2.Lx;
	double  yp1=S2.y0 - 0.5*S2.Ly, yp2=S2.y0 + 0.5*S2.Ly;
	double  zp =S2.z0;

	aux  =   F211(x2-xp1,y-yp1,z2-zp)- F211(x2-xp1,y-yp2,z2-zp);
	aux +=  -F211(x2-xp2,y-yp1,z2-zp)+ F211(x2-xp2,y-yp2,z2-zp);

	aux +=  -F211(x2-xp1,y-yp1,z1-zp)+ F211(x2-xp1,y-yp2,z1-zp);
	aux +=  +F211(x2-xp2,y-yp1,z1-zp)- F211(x2-xp2,y-yp2,z1-zp);

	aux +=  -F211(x1-xp1,y-yp1,z2-zp)+ F211(x1-xp1,y-yp2,z2-zp);
	aux +=  +F211(x1-xp2,y-yp1,z2-zp)- F211(x1-xp2,y-yp2,z2-zp);

	aux +=  +F211(x1-xp1,y-yp1,z1-zp)- F211(x1-xp1,y-yp2,z1-zp);
	aux +=  -F211(x1-xp2,y-yp1,z1-zp)+ F211(x1-xp2,y-yp2,z1-zp);
	//printf("YZ: %lf\n",aux/1.e-27);

	return (aux);
}
/****************************************************************/
double IperpZX(surface_z S1,surface_x S2)
{
	double aux;
	double  x1=S1.x0 - 0.5*S1.Lx, x2=S1.x0 + 0.5*S1.Lx;
	double  y1=S1.y0 - 0.5*S1.Ly, y2=S1.y0 + 0.5*S1.Ly;
	double  z=S1.z0;

	double  xp =S2.x0;
	double  yp1=S2.y0-0.5*S2.Ly, yp2=S2.y0+0.5*S2.Ly;
	double  zp1=S2.z0-0.5*S2.Lz, zp2=S2.z0+0.5*S2.Lz;

	aux  =   F121(x2-xp,y2-yp1,z-zp1) - F121(x2-xp,y2-yp1,z-zp2);
	aux +=  -F121(x2-xp,y2-yp2,z-zp1) + F121(x2-xp,y2-yp2,z-zp2);

	aux +=  -F121(x2-xp,y1-yp1,z-zp1) + F121(x2-xp,y1-yp1,z-zp2);
	aux +=  +F121(x2-xp,y1-yp2,z-zp1) - F121(x2-xp,y1-yp2,z-zp2);

	aux +=  -F121(x1-xp,y2-yp1,z-zp1) + F121(x1-xp,y2-yp1,z-zp2);
	aux +=  +F121(x1-xp,y2-yp2,z-zp1) - F121(x1-xp,y2-yp2,z-zp2);

	aux +=  +F121(x1-xp,y1-yp1,z-zp1) - F121(x1-xp,y1-yp1,z-zp2);
	aux +=  -F121(x1-xp,y1-yp2,z-zp1) + F121(x1-xp,y1-yp2,z-zp2);

	//printf("ZX: %lf\n",aux/1.e-27);

	return (aux);
}

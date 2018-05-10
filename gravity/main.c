#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"

#define REDUCE_YES_NO 'n'

#if REDUCE_YES_NO == 'y'
    #include "reduced.h"
#endif

const int LevelMax = 13;
const double R_G = -0.50;

int main()
{
    size(2.0);
    origin(0.0, 0.0);
    init_grid(64);
    f.sigma = 0.10;
    rho1 = 1000.0;
    rho2 = 1.0;
    mu1 = 1e-1;
    mu2 = 1e-3;
    ;
#if REDUCE_YES_NO == 'y'
    G.x = R_G;
    Z.x = 0.0;
#endif
    ;
    run();
}

#if REDUCE_YES_NO == 'n'
event acceleration(i++)
{
    face vector av = a;
    foreach_face(x)
        av.x[] += R_G;
}
#endif

event init(t = 0)
{
    double x0 = 1.0;
    double r0 = 0.10;
    ;
    refine(sq(x - x0) + sq(y) < sq(1.1 * r0) && sq(x - x0) + sq(y) > sq(0.9 * r0) && level < LevelMax);
    ;
    foreach ()
    {
        f[] = 0.0;
        if (sq(x - x0) + sq(y) < sq(1.0 * r0))
        {
            f[] = 1.0;
        }
    }
}

event adapt(i++)
{
    adapt_wavelet({f, u.x, u.y}, (double[]){1.0e-6, 1.0e-3, 1.0e-3}, maxlevel = LevelMax, minlevel = 6);
}

event writefile(t += 0.1)
{
    char nameFile[500];
    sprintf(nameFile, "file-%.2f", t);
    printf("time: %f      iteration: %d      file is written!\r\n", t, i);
    dump(file = nameFile);
}

event end(t = 1.0)
{
    printf("\r\n-------\r\nEND!\r\n");
}
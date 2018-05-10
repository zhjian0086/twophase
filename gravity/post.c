#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"

#include "view.h"

// install bview: http://basilisk.fr/src/gl/INSTALL
// qcc -Wall -O2 post.c -o post -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

// ./BpostBview tb0.1 ts0.01 te0.2 s10.0 y0.2 p1000
// means:
// time-begin: 0.1
// time-step: 0.01
// time-end: 0.2
// scale: 10.0
// y position from the bottom 0.2 (20 %)
// number of pixels in y-direction: 1000 (in the x direction is twice)

// ffmpeg -framerate 1 -i out-bview-omega-%06d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out-bview-omega-video.mp4

const double R_G = -0.50;

int main()
{
    run();
}

event init (t = 0)
{
    char nameFile[500], nameBview[500];
    double tc, xvelocity, yvelocity, area, xcenter, ycenter;
    FILE *fp;
    fp = fopen ("changemyname.txt", "w");
    fprintf (fp, "time, xcenter(B), xcenter(N), ycenter(B), ycenter(N), xvelocity(B), xvelocity(N), yvelocity(B), yvelocity(N)\r\n");
    for (tc = 0.0; tc <= 1.0; tc += 0.1)
    {
        printf ("load time %.2f\r\n", tc);
        sprintf (nameFile, "file-%.2f", tc);
    	restore (file = nameFile);
        ;
        area = 0.0;
        xvelocity = 0.0;
        yvelocity = 0.0;
        xcenter = 0.0;
        ycenter = 0.0;
        foreach()
        {
            area += Delta * Delta * f[];
            xvelocity += u.x[] * Delta * Delta * f[];
            yvelocity += u.y[] * Delta * Delta * f[];
            xcenter += x * Delta * Delta * f[];
            ycenter += y * Delta * Delta * f[];
        }
        xvelocity /= area;
        yvelocity /= area;
        xcenter /= area;
        ycenter /= area;
        printf ("Bulk mass ___ Numerical\r\n");
        printf ("x-center-pos -> %f ___ %f\r\n", 1.0 + 0.50 * R_G * tc * tc, xcenter);
        printf ("y-center-pos -> %f ___ %f\r\n", 0.0, ycenter);
        printf ("x - velocity -> %f ___ %f\r\n", R_G * tc, xvelocity);
        printf ("y - velocity -> %f ___ %f\r\n", 0.0, yvelocity);
        fprintf (fp, "%.3f, %f, %f, %f, %f, %f, %f, %f, %f\r\n", tc, 1.0 + 0.50 * R_G * tc * tc, xcenter, 0.0, ycenter, R_G * tc, xvelocity, 0.0, yvelocity);
        ;
        sprintf(nameBview, "out-bview-%.2f.png", tc);
        view(width = 1000, height = 500, quat = {0, 0, -0.707, 0.707}, sx = 1.0, sy = 1.0, ty = -0.50);
        clear();
        draw_vof("f", lw = 5);
        squares("u.x");
        box(notics = true);
        //cells(lw = 1); // showing the grid
        mirror({-1})
        {
            draw_vof("f", lw = 5);
            squares("u.y", linear = true); // , min = -10.0, max = 10.0
            box(notics = true);
        }
        draw_string("x-velocity ___ y-velocity", lw = 5, size = 100);
        save(nameBview);
    }
    fclose(fp);
}

event end(t = 0.0)
{
    printf("\r\n-------\r\nEND!\r\n");
}
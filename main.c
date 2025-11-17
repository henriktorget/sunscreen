#include <stdio.h>
#include <math.h>

int dmyjdn(int day, int month, int year) 
// Day, month year to julian date number
{
    int m = (int) (((float) month - 14) / 12);
    int y1 = year + 4800 + m;
    int y2 = (int) (year + 4900 + m) / 100;

    int jdn1 = (int) (1461 * y1) / 4;
    int jdn2 = (367 * (month - 2 - 12 * m)) / 12;
    int jdn3 = ( (int) (3 * y2 ) / 4);

    return (jdn1 + jdn2 - jdn3 + day - 32075);

}

double hmstjdf(int hour, int minute, int second, int timezone)
// Julian date fraction, from hour, minute and second. 
{

    double h = ((double) hour + timezone - 12) / 24;
    double m = (double) minute / 1440;
    double s = (double) second / 86400;

    return h + m + s;

}

double ydt(int year) 
// Delta time, from year. Currently only supports 2005 -> 2050
{

    int t = year - 2000;

    return 62.92 + 0.32217 * t + 0.005589 * (t * t);
}

double hl0(double jme)
// Heliocentric L from julian ephemeris millennium
{
    double l0terms[64][3] =
    {
        {175347046, 0,        0},
        {3341656,   4.6692568, 6283.07585},
        {34894,     4.6261,    12566.1517},
        {3497,      2.7441,    5753.3849},
        {3418,      2.8289,    3.5231},
        {3136,      3.6277,    77713.7715},
        {2676,      4.4181,    7860.4194},
        {2343,      6.1352,    3930.2097},
        {1324,      0.7425,    11506.7698},
        {1273,      2.0371,    529.691},
        {1199,      1.1096,    1577.3435},
        {990,       5.233,     5884.927},
        {902,       2.045,     26.298},
        {857,       3.508,     398.149},
        {780,       1.179,     5223.694},
        {753,       2.533,     5507.553},
        {505,       4.583,     18849.228},
        {492,       4.205,     775.523},
        {357,       2.92,      0.067},
        {317,       5.849,     11790.629},
        {284,       1.899,     796.298},
        {271,       0.315,     10977.079},
        {243,       0.345,     5486.778},
        {206,       4.806,     2544.314},
        {205,       1.869,     5573.143},
        {202,       2.458,     6069.777},
        {156,       0.833,     213.299},
        {132,       3.411,     2942.463},
        {126,       1.083,     20.775},
        {115,       0.645,     0.98},
        {103,       0.636,     4694.003},
        {102,       0.976,     15720.839},
        {102,       4.267,     7.114},
        {99,        6.21,      2146.17},
        {98,        0.68,      155.42},
        {86,        5.98,      161000.69},
        {85,        1.3,       6275.96},
        {85,        3.67,      71430.7},
        {80,        1.81,      17260.15},
        {79,        3.04,      12036.46},
        {75,        1.76,      5088.63},
        {74,        3.5,       3154.69},
        {74,        4.68,      801.82},
        {70,        0.83,      9437.76},
        {62,        3.98,      8827.39},
        {61,        1.82,      7084.9},
        {57,        2.78,      6286.6},
        {56,        4.39,      14143.5},
        {56,        3.47,      6279.55},
        {52,        0.19,      12139.55},
        {52,        1.33,      1748.02},
        {51,        0.28,      5856.48},
        {49,        0.49,      1194.45},
        {41,        5.37,      8429.24},
        {41,        2.4,       19651.05},
        {39,        6.17,      10447.39},
        {37,        6.04,      10213.29},
        {37,        2.57,      1059.38},
        {36,        1.71,      2352.87},
        {36,        1.78,      6812.77},
        {33,        0.59,      17789.85},
        {30,        0.44,      83996.85},
        {30,        2.74,      1349.87},
        {25,        3.16,      4690.48}
    };

    double l0 = 0;

    for (int i = 0; i < 64; i++)
    {
        l0 += l0terms[i][0] * cos(l0terms[i][1] + (l0terms[i][2] * jme));
    }

    return l0;

}

double hl1(double jme)
{
    double l1terms[34][3] =
    {
        {628331966747.0, 0,        0},
        {206059,         2.678235, 6283.07585},
        {4303,           2.6351,   12566.1517},
        {425,            1.59,     3.523},
        {119,            5.796,    26.298},
        {109,            2.966,    1577.344},
        {93,             2.59,     18849.23},
        {72,             1.14,     529.69},
        {68,             1.87,     398.15},
        {67,             4.41,     5507.55},
        {59,             2.89,     5223.69},
        {56,             2.17,     155.42},
        {45,             0.40,     796.30},
        {36,             0.47,     775.52},
        {29,             2.65,     7.11},
        {21,             5.34,     0.98},
        {19,             1.85,     5486.78},
        {19,             4.97,     213.30},
        {17,             2.99,     6275.96},
        {16,             0.03,     2544.31},
        {16,             1.43,     2146.17},
        {15,             1.21,     10977.08},
        {12,             2.83,     1748.02},
        {12,             3.26,     5088.63},
        {12,             5.27,     1194.45},
        {12,             2.08,     4694.00},
        {11,             0.77,     553.57},
        {10,             1.30,     6286.60},
        {10,             4.24,     1349.87},
        { 9,             2.70,     242.73},
        { 9,             5.64,     951.72},
        { 8,             5.30,     2352.87},
        { 6,             2.65,     9437.76},
        { 6,             4.67,     4690.48}
    };

    double l1 = 0;

    for (int i = 0; i < 34; i++)
    {
        l1 += l1terms[i][0] * cos(l1terms[i][1] + (l1terms[i][2] * jme));
    }

    return l1;
}

double hl2(double jme)
{
    double l2terms[20][3] =
    {
        {52919, 0,       0},
        {8720,  1.0721,  6283.0758},
        {309,   0.867,   12566.152},
        {27,    0.05,    3.52},
        {16,    5.19,    26.30},
        {16,    3.68,    155.42},
        {10,    0.76,    18849.23},
        { 9,    2.06,    77713.77},
        { 7,    0.83,    775.52},
        { 5,    4.66,    1577.34},
        { 4,    1.03,    7.11},
        { 4,    3.44,    5573.14},
        { 3,    5.14,    796.30},
        { 3,    6.05,    5507.55},
        { 3,    1.19,    242.73},
        { 3,    6.12,    529.69},
        { 3,    0.31,    398.15},
        { 3,    2.28,    553.57},
        { 2,    4.38,    5223.69},
        { 2,    3.75,    0.98}
    };


    double l2 = 0;

    for (int i = 0; i < 20; i++)
    {
        l2 += l2terms[i][0] * cos(l2terms[i][1] + (l2terms[i][2] * jme));
    }

    return l2;
}

double hl3(double jme)
{
    double l3terms[7][3] =
    {
        {289, 5.844, 6283.076},
        {35,  0,     0},
        {17,  5.49,  12566.15},
        {3,   5.20,  155.42},
        {1,   4.72,  3.52},
        {1,   5.30,  18849.23},
        {1,   5.97,  242.73}
    };

    double l3 = 0;

    for (int i = 0; i < 7; i++)
    {
        l3 += l3terms[i][0] * cos(l3terms[i][1] + (l3terms[i][2] * jme));
    }

    return l3;
}

double hl4(double jme)
{
    double l4terms[3][3] =
    {
        {114, 3.142, 0},
        {  8, 4.13,  6283.08},
        {  1, 3.84,  12566.15}
    };

    double l4 = 0;

    for (int i = 0; i < 3; i++)
    {
        l4 += l4terms[i][0] * cos(l4terms[i][1] + (l4terms[i][2] * jme));
    }

    return l4;
}

double hb0(double jme)
{
    double b0terms[5][3] =
    {
        {280, 3.199, 84334.662},
        {102, 5.422, 5507.553},
        { 80, 3.88,  5223.69},
        { 44, 3.70,  2352.87},
        { 32, 4.00,  1577.34}
    };

    double b0 = 0;

    for (int i = 0; i < 5; i++)
    {
        b0 += b0terms[i][0] * cos(b0terms[i][1] + (b0terms[i][2] * jme));
    } 

    return b0;
}

double hb1(double jme)
{
    double b1terms[2][3] =
    {
        {9, 3.90, 5507.55},
        {6, 1.73, 5223.69}
    };

    double b1 = 0;

    for (int i = 0; i < 2; i++)
    {
        b1 += b1terms[i][0] * cos(b1terms[i][1] + (b1terms[i][2] * jme));
    } 

    return b1;

}

double hr0(double jme)
{
    double r0terms[40][3] =
    {
        {100013989, 0,      0},
        {1670700,   3.0984635, 6283.07585},
        {13956,     3.05525,   12566.1517},
        {3084,      5.1985,    77713.7715},
        {1628,      1.1739,    5753.3849},
        {1576,      2.8469,    7860.4194},
        {925,       5.453,     11506.77},
        {542,       4.564,     3930.21},
        {472,       3.661,     5884.927},
        {346,       0.964,     5507.553},
        {329,       5.90,      5223.694},
        {307,       0.299,     5573.143},
        {243,       4.273,     11790.629},
        {212,       5.847,     1577.344},
        {186,       5.022,     10977.079},
        {175,       3.012,     18849.228},
        {110,       5.055,     5486.778},
        {98,        0.89,      6069.78},
        {86,        5.69,      15720.84},
        {86,        1.27,      161000.69},
        {65,        0.27,      17260.15},
        {63,        0.92,      529.69},
        {57,        2.01,      83996.85},
        {56,        5.24,      71430.7},
        {49,        3.25,      2544.31},
        {47,        2.58,      775.52},
        {45,        5.54,      9437.76},
        {43,        6.01,      6275.96},
        {39,        5.36,      4694.00},
        {38,        2.39,      8827.39},
        {37,        0.83,      19651.05},
        {37,        4.90,      12139.55},
        {36,        1.67,      12036.46},
        {35,        1.84,      2942.46},
        {33,        0.24,      7084.9},
        {32,        0.18,      5088.63},
        {32,        1.78,      398.15},
        {28,        1.21,      6286.6},
        {28,        1.90,      6279.55},
        {26,        4.59,      10447.39}
    };


    double r0 = 0;

    for (int i = 0; i < 40; i++)
    {
        r0 += r0terms[i][0] * cos(r0terms[i][1] + (r0terms[i][2] * jme));
    } 

    return r0;

}

double hr1(double jme)
{
    double r1terms[10][3] =
    {
        {103019, 1.10749, 6283.07585},
        {1721,   1.06440, 12566.1517},
        {702,    3.14200, 0},
        {32,     1.02,    18849.23},
        {31,     2.84,    5507.55},
        {25,     1.32,    5223.69},
        {18,     1.42,    1577.34},
        {10,     5.91,    10977.08},
        { 9,     1.42,    6275.96},
        { 9,     0.27,    5486.78}
    };



    double r1 = 0;

    for (int i = 0; i < 10; i++)
    {
        r1 += r1terms[i][0] * cos(r1terms[i][1] + (r1terms[i][2] * jme));
    } 

    return r1;

}

double hr2(double jme)
{
    double r2terms[6][3] =
    {
        {4359, 5.7846, 6283.0758},
        {124,  5.579,  12566.152},
        {12,   3.14,   0},
        { 9,   3.63,   77713.77},
        { 6,   1.87,   5573.14},
        { 3,   5.47,   18849.23}
    };




    double r2 = 0;

    for (int i = 0; i < 6; i++)
    {
        r2 += r2terms[i][0] * cos(r2terms[i][1] + (r2terms[i][2] * jme));
    } 

    return r2;

}

double hr3(double jme)
{
    double r3terms[2][3] =
    {
        {145, 4.273, 6283.076},
        {  7, 3.92,  12566.15}
    };





    double r3 = 0;

    for (int i = 0; i < 2; i++)
    {
        r3 += r3terms[i][0] * cos(r3terms[i][1] + (r3terms[i][2] * jme));
    } 

    return r3;

}


int main(void){

    // Input
    
    int day = 11;
    int month = 11;
    int year = 2025;
    int hour = 16;
    int minute = 33;
    int second = 14;
    int timezone = 0;

    // Calculate julian date from gregorian

    int jdn = dmyjdn(day, month, year);
    double jdf = hmstjdf(hour, minute, second, timezone);
    double jd = jdf + (double) jdn;
    double dt = ydt(year);
    double jde = jd + dt / 86400;
    double jc = (jd - (double) 2451545) / 36525;
    double jce = (jde - (double) 2451545) / 36525;
    double jme = jce / 10; 

    printf("jdn  - %d\n", jdn);
    printf("jdf  - %f\n", jdf);
    printf("jd   - %f\n", jd);
    printf("dt   - %f\n", dt);
    printf("jde  - %f\n", jde);
    printf("jc   - %f\n", jc);
    printf("jce  - %f\n", jce);
    printf("jme  - %f\n", jme);

    // Calculate earth heliocentric lon, lat and radius vector

    double l0 = hl0(jme);
    double l1 = hl1(jme);
    double l2 = hl2(jme);
    double l3 = hl3(jme);
    double l4 = hl4(jme);
    double l5 = cos(3.14);

    printf("l0  - %f\n", l0);
    printf("l1  - %f\n", l1);
    printf("l2  - %f\n", l2);
    printf("l3  - %f\n", l3);
    printf("l4  - %f\n", l4);
    printf("l5  - %f\n", l5);

    double lrad = ( l0 + l1 * jme + l2 * pow(jme, 2) + l3 * pow(jme, 3) + l4 * 
                  pow(jme, 4) + l5 * pow(jme, 5)) / pow(10, 8);

    printf("lrad  - %f\n", lrad);
    
    double ldeg = (lrad * 180) / 3.1415926535898;
    
    printf("ldeg  - %f\n", ldeg);

    double lf = ldeg / 360;
    double lk = floor(lf);
    double l = ldeg - 360 * lk;

    printf("l - %f\n", l);

    double b0 = hb0(jme);
    double b1 = hb1(jme);

    printf("b0 - %f\n", b0);
    printf("b1 - %f\n", b1);

    double brad = (b0 + b1 * jme) / pow(10, 8);
    double b = (brad * 180) / 3.1415926535898;

    printf("b - %f\n", b);

    double r0 = hr0(jme);
    double r1 = hr1(jme);
    double r2 = hr2(jme);
    double r3 = hr3(jme);
    double r4 = 4 * cos(2.56 + (6283.08 * jme));

    printf("r0 - %f\n", r0);
    printf("r1 - %f\n", r1);
    printf("r2 - %f\n", r2);
    printf("r3 - %f\n", r3);
    printf("r4 - %f\n", r4);

    double r = ( r0 + r1 * jme + r2 * pow(jme, 2) + r3 * pow(jme, 3) + r4 * 
                  pow(jme, 4)) / pow(10, 8);

    printf("r - %f\n", r);

    

    return 0;
}


p=0;
#include "IBC/LBM/HRRp1.h"
p=1;
#include "IBC/LBM/HRRp2_p7.h"
p=2;
#include "IBC/LBM/HRRp2_p7.h"
p=3;
#include "IBC/LBM/HRRp2_p7.h"
p=4;
#include "IBC/LBM/HRRp2_p7.h"
p=5;
#include "IBC/LBM/HRRp2_p7.h"
p=6;
#include "IBC/LBM/HRRp2_p7.h"
p=7;
#include "IBC/LBM/HRRp8_p11.h"
p=8;
#include "IBC/LBM/HRRp8_p11.h"
p=9;
#include "IBC/LBM/HRRp8_p11.h"
p=10;
#include "IBC/LBM/HRRp8_p11.h"
p=11;
#include "IBC/LBM/HRRp12_p15.h"
p=12;
#include "IBC/LBM/HRRp12_p15.h"
p=13;
#include "IBC/LBM/HRRp12_p15.h"
p=14;
#include "IBC/LBM/HRRp12_p15.h"
p=15;
#include "IBC/LBM/HRRp16_p19.h"
p=16;
#include "IBC/LBM/HRRp16_p19.h"
p=17;
#include "IBC/LBM/HRRp16_p19.h"
p=18;
#include "IBC/LBM/HRRp16_p19.h"

a1hrr_1 = a1hrr_1 + sigma*a1pr_1 ;//xx
a1hrr_2 = a1hrr_2 + sigma*a1pr_2 ;//xy
a1hrr_3 = a1hrr_3 + sigma*a1pr_3 ;//xz
a1hrr_4 = a1hrr_4 + sigma*a1pr_4 ;//yx
a1hrr_5 = a1hrr_5 + sigma*a1pr_5 ;//yy
a1hrr_6 = a1hrr_6 + sigma*a1pr_6 ;//yz
a1hrr_7 = a1hrr_7 + sigma*a1pr_7 ;//zx
a1hrr_8 = a1hrr_8 + sigma*a1pr_8 ;//zy
a1hrr_9 = a1hrr_9 + sigma*a1pr_9 ;//zz

a3neq_1 = 2.*Ux*a1hrr_2 +    Uy*a1hrr_1;
a3neq_2 = 2.*Ux*a1hrr_3 +    Uz*a1hrr_1;
a3neq_3 =    Ux*a1hrr_5 + 2.*Uy*a1hrr_2;
a3neq_4 =    Ux*a1hrr_9 + 2.*Uz*a1hrr_3;
a3neq_5 =    Uy*a1hrr_9 + 2.*Uz*a1hrr_6;
a3neq_6 = 2.*Uy*a1hrr_6 +    Uz*a1hrr_5;

coef1 = a3neq_1+a3neq_5;
coef2 = a3neq_4+a3neq_3;
coef3 = a3neq_6+a3neq_2;
coef4 = a3neq_1-a3neq_5;
coef5 = a3neq_4-a3neq_3;
coef6 = a3neq_6-a3neq_2;



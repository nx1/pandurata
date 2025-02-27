Debugging Commands
------------------
b main
run

display yn
display y
display y1
display y2
display y_ck
display del
display del_ck
display v_
display ph_v
display ph_v_p
display p_
display po_
display ph_p
display p_hat
display n_hat
display n_p_hat
display ph_v_hat
display f_
display f0_
display E0
display E_f
display E_i
display steps

E0
--
The energy or norm of a four-velocity for a photon (timelike trajectory) should always be:
$$ E_0 = g_{\mu\nu} v^\mu v^\nu = -1 $$

Similarly Minkowski norm squared of the four-momentum for a photon should be:
$$ >p \cdot p = \eta_{\mu\nu} p^\mu p^\nu = p_\nu p^\nu = -{E^2 \over c^2} + |\mathbf p|^2 = -m^2 c^2 = 0$$

This holds for a while until after a certain amount of scattering steps, this no longer is true
with E0 having extremely large values.

The specific variable E0 has different calculations throughout the code, these are not equivalent.

Line 1629
```
E0 = g_dn_ph[0][0]    * v_[0] * v_[0]
   + 2.*g_dn_ph[0][3] * v_[0] * v_[3]
   + g_dn_ph[1][1]    * v_[1] * v_[1]
   + g_dn_ph[2][2]    * v_[2] * v_[2]
   + g_dn_ph[3][3]    * v_[3] * v_[3];
```

Line 1843 (same as above but uses k_ = ph_v
```
E0 = g_dn_ph[0][0] * k_[0] * k_[0]
+ 2.*g_dn_ph[0][3] * k_[0] * k_[3]
   + g_dn_ph[1][1] * k_[1] * k_[1]
   + g_dn_ph[2][2] * k_[2] * k_[2]
   + g_dn_ph[3][3] * k_[3] * k_[3];
```

$$ p \dot p = p_{\nu} p^{\nu} = -p^0 p^0 + p^1 p^1 + p^2 p^2 + p^3 p^3 = 0 $$ for a photon
Line 1859
```
E0 = - ph_v_hat[0]*ph_v_hat[0]
     + ph_v_hat[1]*ph_v_hat[1]
     + ph_v_hat[2]*ph_v_hat[2]
     + ph_v_hat[3]*ph_v_hat[3];
```

Metric Tensor
-------------
The values for metric tensors for g_dn_ph and g_up_ph numerically seem okay, they only depend
on r and theta, (y[1] and y[2]) which also stay in reasonable values so I think this is fine.
see: calc_g.c


Large Velocities
----------------
Displaying the variables after a few scattering steps shows that the four-velocity/momentums are blowing
up. Which in turn is probably what is causing what should be Lorentz invariant energy quantities to also
explode.

yn     = {16.17, 1.49, 1.36, 3.57, 26108191945.85, -58082609630.26, -2177840683.17, -111443265954.37}
y      = {16.10, 1.49, 1.36, 3.55, 26108196543.17, -50732885498.87, -2352283491.96, -111443265954.37}
y2     = {16.10, 1.49, 1.36, 3.57, 26108196543.17, -58082609630.26, -2177840683.17, -111443265954.37}
y_ck   = {16.10, 1.49, 1.36, 3.57, 26108196543.17, -58082609630.26, -2177840683.17, -111443265954.37}

v_       = {-238.61, -0.25, 0.05, -56.11}
ph_v     = {505457216208.37, -1339019729.72, -962485244.44, 119131863260.28}
ph_v_p   = {495950988197.86, -170254163.20, -1504497254.70, 116245422811.50}
p_       = {26108191945.85, -58082609630.26, -2177840683.17, -111443265954.37}

ph_p     = {0.23, -3.12, -0.53, -1.38}

p_hat    = {-0.61, -0.70, 0.39}
n_p_hat  = {0.97, -0.09, -0.23}
ph_v_hat = {-36999618428.74, -36122144382.08, -7915068973.84, -1230500816.43} 

f_  = {-17.51, 0.03, -0.66, -4.11}
f0_ = {3.77e-15, -0.15, 0.04, 0.33}

E0    = 7994701568
E_f   = 4647491668541440
E_i   = -7120707127.34
steps = 213


One thing that seems weird to me is that ph_v_hat is not a unit vector, this is in contrast to all
the other hat parameters in the code:

f_hat   = {0.72, -0.66, 0.22}
p_hat   = {-0.61, -0.70, 0.39}
n_hat   = {0.58, 0.80, 0.18}
r_hat   = {0.48, 0.42, -0.09}
fp_hat  = {0.69, 0.63, -0.35}
z_hat   = {0.00, 0.00, 1.00}
n_p_hat = {-0.28, -0.74, -0.61}
e_x_hat = {-0.43, -0.81, 0.40}
e_y_hat = {0.90, -0.43, 0.10}
e_z_hat = {0.09, 0.40, 0.91}
f_v_hat = {0.00, -0.77, 0.54, -0.34}

could it be that we just need to call normalize(ph_v_hat) ?
There seems to be no need for it in the thin disc code...


cashkarp.c
----------
A cash-karp method is used to solve a set of ODEs,
is used to update the state vector (y) by a time step (dt).

The time step (dt) is set via

Line 1076:
```
if (erro > 0) {
   dt = 0.9*dt*pow((erro/err),0.2);
} else {
   dt = 0.9*dt*pow((erro/err),0.25);
}
if ((yn[1] < rr[Nr/2])&&(dt > yn[1]/20.)) {
   dt = yn[1]/20.;
}
```

err:

   1e-8
1.25e-11
1.10e-11
2.81e-12

where
erro = 





$$ dt = 0.9 \cdot dt \cdot \left( \frac{\text{erro}}{\text{err}} \right)^{0.2}, \quad \text{if } \text{erro} > 0 $$
$$ dt = 0.9 \cdot dt \cdot \left( \frac{\text{erro}}{\text{err}} \right)^{0.25}, \quad \text{if } \text{erro} \leq 0 $$
$$ \text{If } yn[1] < rr[Nr/2] \text{ and } dt > \frac{yn[1]}{20}, \text{ then } dt = \frac{yn[1]}{20} $$



cashkarp(dt, y, y_ck, del_ck);

it uses 6 function evaluations, each time calculating gradients using the accel()
function in accel.c

accel(y0,k1);
accel(yn,k2);
...
accel(yn,k6);

The accel function only gives derivatives for certain state vectors, this means
that some are not affected by the cash-karp directly, namely:
  
  dy[0] = 0; # time 
  dy[4] = 0; # 0th component of the 4-velocity / momentum
  dy[7] = 0; # 4th component of the 4-velocity / momentum (phi component)
  
yet the y[7] is still blowing up, which means this must be happening outside of the cash-karp step
atleast for this.

accel.c
-------
This function is responsible for calculating the gradients for state vector.
accel(double y[], double dy[])

It updates the dy values in-place based on some GR calculations and existing values of y

dy[0] = 0 (time) (controlled through timestep dt)
dy[1] =   (r     position)
dy[2] =   (theta position)
dy[3] =   (phi   position)
dy[4] = 0 (Energy) Conserved
dy[5] =   (p_r)
dy[6] =   (p_theta)
dy[7] = 0 (p_phi)  Conserved (conservation of angular momentum?)

There is some run-away issue with the four momentum

The gradient values for the p_r and p_theta are calculated through

dy[5] = - domg_dr*y[7]
        + 0.5*al2*o_p_mod*(2.0*o_Sig*(r-M-r*Del*o_Sig)*y[5]*y[5] - 2.0*r*o_Sig2*y[6]*y[6] + dobr_dr*y[7]*y[7])
        + p_mod/alp*dalp_dr;

dy[6] = - domg_dt*y[7]
        + 0.5*al2*o_p_mod*(2.*a2*sthcth*o_Sig2*(Del*y[5]*y[5]+y[6]*y[6]) + dobr_dt*y[7]*y[7])
        + p_mod/alp*dalp_dt;
		
debugging vars: (after many runs)
						ULX Simulation       Thin Disc
display y[5]           -9031759795471.2773   0.6946     !  
display dy[5]          -9057963982161.4414   0.0642     !
display y[6]           12726295505.6018      1.9928     !
display dy[6]          202450956362.4897     -0.5082    !
display y[7]           -9697128920091.7969   1.4087     !
display dy[7]          0.0000                0.0000    
display domg_dr        -0.3580               -0.0341    
display domg_dt        0.0003                -0.0006    
display dobr_dr        -0.1487               -0.1448    
display dobr_dt        -0.1552               0.8032    
display dalp_dr        1.0504                0.1101
display dalp_dt        -0.0084               0.0118
display al2            0.0125                0.4390
display o_p_mod        -1.7448e-12           -1.4518
display o_Sig          0.4423                0.0845
display o_Sig2         0.1957                0.0071
display r              1.4906                3.3571
display M              1.0                   1.0
display alp            0.1119                0.6626
display a2             0.8100                0.8100
display sthcth         0.2132                -0.4619
display o_Sig2         0.1957                0.0071    
display Del            0.0507                5.3661    
display p_mod          -573121192854.2236    -0.6888   !
display alp            0.1119                0.6626    

All three terms that compose dy[5] and dy[6] (by terms i mean parts of the sum on each line)
depend in some way on y[6] or y[7] (p_mod depends on them too)

p_mod   = y[4]+omg*y[7];

But note also, that y[7] which should not be updated here since dy[7] = 0 has also exploded.
This means that this is becoming large elsewhere.


yn[7] is updated on lines 1549 & 2066 using the same command
yn[7] = g_dn_ph[0][3]*ph_v_p[0]+g_dn_ph[3][3]*ph_v_p[3];

This is causing large jumps
Old value = -51.035163494211218
New value = -274.06950080424929

Old value = -274.06950080424929
New value = -40756.620070307777


but we're at the same thing with  ph_v_p having becoming really large for some reason

Tensor is fine I'm pretty sure:
g_dn_ph[0][0] = 0.30535696963551939
g_dn_ph[1][1] = 33.31604469954646
g_dn_ph[2][2] = 2.3125081294250074
g_dn_ph[0][3] = -1.1248576023477985
g_dn_ph[3][3] = 3.9260454910448006

ph_v_p = {158362.77084587753, 119.45157239510925, -642.52607207884296, 34991.684879326378}



















































		

For the thin simulation these values stay relatively stable  < 1 throughout the simulation
however for the ULX simulation, there is a compounding effect that causes them to both blow up over
time.



boost() function
----------------
There is a function called boost defined in tensor_math.c that is called to update ph_v_hat and f_v_hat


Line 1886
```
boost(beta,n_hat,ph_v_hat);
boost(beta,n_hat,f_v_hat);

rdsh = ph_v_hat[0] / E_i;
for (j=0;j<=3;j++) ph_v_hat[j]=ph_v_hat[j]/rdsh;
```


Line 1984
```
//BOOST BACK INTO CORONA FRAME, THEN CONVERT TO COORDINATE BASIS
beta = -beta;
//printf("post-scatter p0  %12.5g\n",ph_v_hat[0]);
E_f = ph_v_hat[0];

E0 = - ph_v_hat[0]*ph_v_hat[0]
	 + ph_v_hat[1]*ph_v_hat[1]
	 + ph_v_hat[2]*ph_v_hat[2]
	 + ph_v_hat[3]*ph_v_hat[3];
// if (fabs(E0) > 1e-4) {
//   printf("neg energy %ld %ld %ld %ld %g %g %g %g %g %g\n",it,ip,ir,iph,
//       E0,dot_g4(g_up_ph,p_,p_),dot_g4(g_dn_ph,v_,v_),yn[1],yn[2],rdsh);
//       sleep(1);
//}

boost(beta,n_hat,ph_v_hat);
boost(beta,n_hat,f_v_hat);

rdsh = ph_v_hat[0] / E_f;
```

The definition of boost in 

void boost(double beta, double n_[], double p_[])
{
  double gamma,Lambda[4][4],p_p[4];
  int i,j;
  gamma = 1./sqrt(1.-beta*beta);
  Lambda[0][0]=gamma;
  for (i=1;i<=3;i++) {
    Lambda[0][i] = -beta*gamma*n_[i-1];
    Lambda[i][0] = Lambda[0][i];
    for (j=1;j<=3;j++) {
      Lambda[i][j] = (gamma-1.)*n_[i-1]*n_[j-1];
    }
    Lambda[i][i] = Lambda[i][i]+1.;
  }
  for (i=0;i<=3;i++) {
    p_p[i]=0;
    for (j=0;j<=3;j++) p_p[i]=p_p[i]+Lambda[i][j]*p_[j];
  }
  for (i=0;i<=3;i++) p_[i]=p_p[i];
}

it simply update p_ in direction n_ by factor beta, I think this is fine.



LU Decomposition
----------------
https://en.wikipedia.org/wiki/LU_decomposition
https://docs.itascacg.com/itasca900/common/fish/doc/fish_manual/fish_fish/matrix_utilities/fish_matrix.ludcmp.html
https://docs.itascacg.com/flac3d700/common/fish/doc/fish_manual/fish_fish/matrix_utilities/fish_matrix.lubksb.html

adata & adata1    A : Matrix   the 1 versions just have 1-indexing since the following functions are 1-indexed
bdata & bdata1    B : Vector   

solves Ax = B

ludcmp_js(adata1,4,indx,&dd);     // LU decomposition of matrix adata1, result is stored in-place
lubksb_js(adata1,4,indx,bdata1);  // Used the decomposed matrix to solve the system of equations through back-substitution
                                  // Result replaced for x is stored in place of B (bdata1)

The vector B holds
bdata1[1] = kap1i;
bdata1[2] = kap2i;
bdata1[3] = 0;
bdata1[4] = 0;

While the matrix A holds:
adata1[1][1] = -aa * cos(t)*k_[1] + r*sin(t)*aa*k_[2];
adata1[1][2] = aa*cos(t)*k_[0] - aa*sin(t)*sin(t)*k_[3];
adata1[1][3] = r*sin(t) * ((r2+a2)*k_[3]-aa*k_[0]);
adata1[1][4] = a2*cos(t)*sin(t)*sin(t)*k_[1] + r*sin(t)*(-(r2+a2)*k_[2]);
adata1[2][1] = r*(-k_[1])+a2*sin(t)*cos(t) * k_[2];
adata1[2][2] = r*(k_[0]-aa*sin(t)*sin(t) * k_[3]);
adata1[2][3] = aa*sin(t)*cos(t)*((r2+a2) * k_[3]-aa*k_[0]);
adata1[2][4] = r*aa*sin(t)*sin(t)*k_[1] - aa*sin(t)*cos(t)*(-(r2+a2)*k_[2]);
for (j=0;j<=3;j++) {
   adata1[3][j+1]=p_[j];
   adata1[4][j+1]=po_[j];
}


ludcmp_js
---------
Remember the arrays are 1 indexed for these functions so the 0th index is just empty

// Matrix LU decomposition
// don't need to show 0th row because it's empty.
// 0th column is also all zeros
break ludcmp_js
display *a[1]@5
display *a[2]@5
display *a[3]@5
display *a[4]@5
display *vv@5

// Use this to break at the end of the ludcmp_js
break 65

This seems okay, matricies after 5000 runs of ludcmp_js look like,

// Thin Simulation
*a[1]@5 = { 0.0000,  2.1301, -1.0609,  1.5349, -0.1273 }
*a[2]@5 = { 0.0000, -0.0766,  0.2452,  0.1176,  0.0000 }
*a[3]@5 = { 0.0000,  3.2105,  0.0000, -1.1726, -2.7756 }
*a[4]@5 = { 0.0000,  0.0172,  0.0743,  0.0299, -0.6102 }

// ULX Simulation
*a[1]@5 = { 0.0000,  1.0237,  0.0384, -0.0983,  0.0120 }
*a[2]@5 = { 0.0000,  0.0360,  0.9815,  0.0035,  0.0000 }
*a[3]@5 = { 0.0000, -0.0029,  0.0000,  0.0306,  0.0000 }
*a[4]@5 = { 0.0000, -0.0002,  0.0000, -0.0007, -0.0177 }


lubksb_js
---------
break ludcmp_js
display *a[1]@5
display *a[2]@5
display *a[3]@5
display *a[4]@5
display *b@5

After 5000 calls runs:

// Thin simulation
*a[1]@5 = [  0.0000   4.9743  -4.1380   1.5686   0.6174 ]
*a[2]@5 = [  0.0000  -0.0540   0.0944   0.0847  -0.0000 ]
*a[3]@5 = [  0.0000  -0.0167  -0.7307   0.0880  -0.6579 ]
*a[4]@5 = [  0.0000   3.1928  -0.0000 -16.1631 -10.6341 ]

*b@5    = [  0.0000   0.0000  -0.0354   0.0551   0.0053 ]


// ULX Simulation
*a[1]@5 = [  0.0000   1.6775   0.0317  -1.0782  -0.0241 ]
*a[2]@5 = [  0.0000   0.0089   0.7916   0.0096   0.0000 ]
*a[3]@5 = [  0.0000  -0.0991   0.0000   0.1728  -0.0000 ]
*a[4]@5 = [  0.0000   0.0017  -0.0001   0.0107  -0.1996 ]

*b@5    = [  0.0000  -0.0000  -0.0092   0.7576   0.0070 ]

there is a step near the end of ludbskb_js that divides by the
diagonal b[i] = sum/a[i][i];
If the diagonal is very small then this can cause the number to explode


kap1i and kap2i
---------------
These two variables are used for the 1st and 2nd entry of bdata
They also depend on ph_v_p so there is a sort of recursive dependence on these.
which causes them to become large.

not sure what they are exactly but presumably some GR quantity related
to the photon's trajectory

kap1i = aa*cos(t)*((ph_v_p[0]*f_[1]-ph_v_p[1]*f_[0])
      + aa*sin(t)*sin(t)*(ph_v_p[1]*f_[3]-ph_v_p[3]*f_[1]))
      + r*sin(t)*((r*r+aa*aa)*(ph_v_p[3]*f_[2]-ph_v_p[2]*f_[3])
      - aa*(ph_v_p[0]*f_[2]-ph_v_p[2]*f_[0]));

kap2i = r*(ph_v_p[0]*f_[1]-ph_v_p[1]*f_[0]
      + aa*sin(t)*sin(t)*(ph_v_p[1]*f_[3]-ph_v_p[3]*f_[1]))
      - aa*sin(t)*cos(t)*((r*r+aa*aa)*(ph_v_p[3]*f_[2]-ph_v_p[2]*f_[3])
      - aa*(ph_v_p[0]*f_[2]-ph_v_p[2]*f_[0]));
	
There is also kap10 and kap20 which are the same but has
ph_v_p --> ph_v
f_     --> f0_


k_
--
bdata1 is also updated to be k_

Line 1852:
for (j=0;j<=3;j++) bdata1[j+1] = k_[j];





for (j=0;j<=3;j++) ph_v_hat[j]=bdata1[j+1];
so ph_v_hat is large:
ph_v_hat = {-116939.55517341282, 10007.029636625028, -12082.723979739048, -115882.38305815132}




ph_v_p[i] += e_lfs[j][i] * ph_v_hat[j];
thus ph_v_p becomes large
ph_v_p = {21406019.784500774, 30641.054405532283, -5127.7881375568077, 4846900.3290679231}







y2 = y_ck
y2 = yn 
y_ck is used for Cash-Karp solving ODE

in cashkarp.c

ph_v_hat is very large: 
{-38487292861.987518, -36689585590.023933, -9274490589.7120113, 7009268534.3561964}
i feel like this should be normalized, isn't _hat supposed to be a unit vector?

ph_v_p:
for (i=0;i<=3;i++) {           (2027)                   
   ph_v_p[i] = 0;                                 
   f_[i] = 0;                                     
   for (j=0;j<=3;j++) {                           
      ph_v_p[i] += e_lfs[j][i]*ph_v_hat[j];       
      f_[i]     += e_lfs[j][i]*f_v_hat[j];        
   }                                              
}                                                 
deg = degp;                                       


yn:

for (j=0;j<=7;j++) {   (1572)
   yn[j]=y2[j];      
}  

yn[4] = g_dn_ph[0][0]*ph_v_p[0] + g_dn_ph[0][3]*ph_v_p[3];  (2048)
yn[5] = g_dn_ph[1][1]*ph_v_p[1];                          
yn[6] = g_dn_ph[2][2]*ph_v_p[2];                          
yn[7] = g_dn_ph[0][3]*ph_v_p[0] + g_dn_ph[3][3]*ph_v_p[3];                  



# This is done at the start and doesn't seem like too much of an issue
ph_v[0] = e_lf[0][0]*(1.)+e_lf[1][0]*p_hat[0] + e_lf[2][0]*p_hat[1]+e_lf[3][0]*p_hat[2];   1014
ph_v[1] = e_lf[0][1]*(1.)+e_lf[1][1]*p_hat[0] + e_lf[2][1]*p_hat[1]+e_lf[3][1]*p_hat[2];
ph_v[2] = e_lf[0][2]*(1.)+e_lf[1][2]*p_hat[0] + e_lf[2][2]*p_hat[1]+e_lf[3][2]*p_hat[2];
ph_v[3] = e_lf[0][3]*(1.)+e_lf[1][3]*p_hat[0] + e_lf[2][3]*p_hat[1]+e_lf[3][3]*p_hat[2];
f0_[0]  = e_lf[0][0]*(0.)+e_lf[1][0]*f_hat[0] + e_lf[2][0]*f_hat[1]+e_lf[3][0]*f_hat[2];
f0_[1]  = e_lf[0][1]*(0.)+e_lf[1][1]*f_hat[0] + e_lf[2][1]*f_hat[1]+e_lf[3][1]*f_hat[2];

yn factors here are large though
ph_v[0] = g_up_ph[0][0]*yn[4] + g_up_ph[0][3]*yn[7];   1052
ph_v[1] = g_up_ph[1][1]*yn[5];                      
ph_v[2] = g_up_ph[2][2]*yn[6];                      
ph_v[3] = g_up_ph[0][3]*yn[4] + g_up_ph[3][3]*yn[7];


ph_v is large
k_[j]=ph_v[j]   (1612)
so k_ becomes large

k_ becomes large 
for (j=0;j<=3;j++) bdata1[j+1]=k_[j];
so bdata1 becomes large



main loop
---------

The condition for continue to scatter a photon is all contained within the loop

while ((yn[1]>1.02*Rhor) && (yn[1]<Rshell*0.99999) && (steps<1000) && (A_fact>0)) {
	
There are various steps inside this loop
1. Updating the photon position & Momentum by time step dt using the Cash Karp function.
2. Disk Interaction
3. Coronal Scattering
4. Polarisation
5. Calculations





//CONVERT BACK FROM DISK FRAME TO COORDINATE BASIS
for (i=0;i<=3;i++) {
  ph_v_p[i]=0;
  f_[i]=0;
  for (j=0;j<=3;j++) {
	 
	 f_[i]+=e_lfs[j][i]*f_v_hat[j];
  }
}

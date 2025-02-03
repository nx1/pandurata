#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define Nr 96
#define Nth 64
#define Nph 48
#define Nstart 624
#define Nstop 624
#define Nstep 1

#define index2(a,b) ((Nfreq+1)*(b)+(a))
#define index3(a,b,c) (8*(N+1)*(c)+8*(b)+(a))
#define index4(a,b,c,d) (4*(Nlat+1)*(N+1)*(d)+4*(Nlat+1)*(c)+4*(b)+(a))
#define indexdex(b,c,d) ((Nlat+1)*(N+1)*(d)+(Nlat+1)*(c)+(b))
#define indexijk(a,b,c) (Nph*Nth*(a)+Nph*(b)+(c))
#define indexthijk(a,b,c) (Nph*(Nth+1)*(a)+Nph*(b)+(c))
#define movdex(a,b) ((N+1)*(b)+(a))

int main(void)
{
  float *rho_ijk,*T_ijk,*bb_ijk,*tau_ijk,*Ut_ijk,*Ur_ijk,*Uz_ijk,*Up_ijk,
    *phtop_ik,*phbot_ik,z1,z2,z3,z4,z5,z6;
  long it,i_start,i_stop,i,j,k,
    ilat,ik,id,ih,iu,ilo,ihi,jlo,jhi,klo,khi,i10,j10,k10;
  FILE *infile,*outfile,
    *rho_file,*T_file,*bb_file,*Ut_file,*Ur_file,*Uz_file,*Up_file,*tau_file;
  char fname[18];

  T_ijk = (float *)malloc(Nr*Nth*Nph*sizeof(float));
  rho_ijk = (float *)malloc(Nr*Nth*Nph*sizeof(float));
  Ut_ijk = (float *)malloc(Nr*Nth*Nph*sizeof(float));
  Ur_ijk = (float *)malloc(Nr*Nth*Nph*sizeof(float));
  Uz_ijk = (float *)malloc(Nr*Nth*Nph*sizeof(float));
  Up_ijk = (float *)malloc(Nr*Nth*Nph*sizeof(float));
  tau_ijk = (float *)malloc(Nr*(Nth+1)*Nph*sizeof(float));
  bb_ijk = (float *)malloc(Nr*Nth*Nph*sizeof(float));

  for (it=Nstart;it<=Nstop;it+=Nstep) {
    //build filename string
    ik = 0;
    id = it;
    ik = (id-fmod(id,1000))/1000;
    id = id-ik*1000;
    ih = (id-fmod(id,100))/100;
    id = id-ih*100;
    iu = fmod(id,10);
    id = (id-fmod(id,10))/10;

    printf("%d %d %d %d %d\n",it,ik,ih,id,iu);

    strcpy(fname,"rh_0000_new.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    rho_file = fopen(fname, "r");
    strcpy(fname,"te_0000_new.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    T_file = fopen(fname, "r");
    strcpy(fname,"u0_0000_new.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    Ut_file = fopen(fname, "r");
    strcpy(fname,"u1_0000_new.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    Ur_file = fopen(fname, "r");
    strcpy(fname,"u2_0000_new.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    Uz_file = fopen(fname, "r");
    strcpy(fname,"u3_0000_new.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    Up_file = fopen(fname, "r");
    strcpy(fname,"bb_0000_new.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    bb_file = fopen(fname, "r");
    /*
    rho_file = fopen("rho_data_new.dat","r");
    T_file = fopen("T_data_new.dat","r");
    Ut_file = fopen("u0_data_new.dat","r");
    Ur_file = fopen("u1_data_new.dat","r");
    Uz_file = fopen("u2_data_new.dat","r");
    Up_file = fopen("u3_data_new.dat","r");
    tau_file = fopen("tau_data_new.dat","r");
    rho_file = fopen("rho_data_fake.dat","r");
    T_file = fopen("T_data_fake.dat","r");
    Ut_file = fopen("u0_data_fake.dat","r");
    Ur_file = fopen("u1_data_fake.dat","r");
    Uz_file = fopen("u2_data_fake.dat","r");
    Up_file = fopen("u3_data_fake.dat","r");
    tau_file = fopen("tau_data_fake.dat","r");
    */
    for (k=0;k<Nph;k++) {
      for (j=0;j<Nth;j++) {
	for (i=0;i<Nr;i++) {
	  fscanf(rho_file,"%f",&z1);
	  rho_ijk[indexijk(i,j,k)]=z1;
	  fscanf(T_file,"%f",&z1);
	  T_ijk[indexijk(i,j,k)]=z1;
	  fscanf(Ut_file,"%f",&z1);
	  Ut_ijk[indexijk(i,j,k)]=z1;
	  fscanf(Ur_file,"%f",&z1);
	  Ur_ijk[indexijk(i,j,k)]=z1;
	  fscanf(Uz_file,"%f",&z1);
	  Uz_ijk[indexijk(i,j,k)]=z1;
	  fscanf(Up_file,"%f",&z1);
	  Up_ijk[indexijk(i,j,k)]=z1;
	  fscanf(bb_file,"%f",&z1);
	  bb_ijk[indexijk(i,j,k)]=z1;
	}
      fscanf(rho_file,"\n");
      fscanf(T_file,"\n");
      fscanf(Ut_file,"\n");
      fscanf(Ur_file,"\n");
      fscanf(Uz_file,"\n");
      fscanf(Up_file,"\n");
      fscanf(bb_file,"\n");
      }
    } 
    fclose(rho_file);
    fclose(T_file);
    fclose(Ut_file);
    fclose(Ur_file);
    fclose(Uz_file);
    fclose(Up_file);
    fclose(bb_file);
     
    strcpy(fname,"ta_0000_new.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    tau_file = fopen(fname, "r");
     for (k=0;k<Nph;k++) {
        for (j=0;j<Nth+1;j++) {
           for (i=0;i<Nr;i++) {
              fscanf(tau_file,"%f",&z1);
              tau_ijk[indexthijk(i,j,k)]=z1;
           }
           fscanf(tau_file,"\n");
        }
     }
    fclose(tau_file);

    strcpy(fname,"rh_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(rho_ijk,sizeof(float),Nr*Nth*Nph,outfile);
    fclose(outfile);

    strcpy(fname,"te_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(T_ijk,sizeof(float),Nr*Nth*Nph,outfile);
    fclose(outfile);
    
    strcpy(fname,"u0_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(Ut_ijk,sizeof(float),Nr*Nth*Nph,outfile);
    fclose(outfile);

    strcpy(fname,"u1_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(Ur_ijk,sizeof(float),Nr*Nth*Nph,outfile);
    fclose(outfile);

    strcpy(fname,"u2_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(Uz_ijk,sizeof(float),Nr*Nth*Nph,outfile);
    fclose(outfile);

    strcpy(fname,"u3_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(Up_ijk,sizeof(float),Nr*Nth*Nph,outfile);
    fclose(outfile);

    strcpy(fname,"bb_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(bb_ijk,sizeof(float),Nr*Nth*Nph,outfile);
    fclose(outfile);

    strcpy(fname,"ta_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(tau_ijk,sizeof(float),Nr*(Nth+1)*Nph,outfile);
    fclose(outfile);
  }
  return(0);
}

// #######################################
// # ZWECK: Konvertieren von binären output-files, 
// # generiert mittels mpiio ins imd chkpt-format
// #############################################


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>

int main(int argc, char** argv)
{

  FILE* infile=fopen(argv[1],"rb");
  char IMD[4];
  long natoms=0;
  short obs=0;
  short disp=0;
  int len;

  double bxx,bxy,bxz;
  double byx,byy,byz;
  double bzx,bzy,bzz;
  len=fread(&IMD, sizeof(char), 4, infile); //IMD
  fseek(infile,3,SEEK_SET);
  len=fread(&disp, sizeof(short),1,infile); //disp
  fseek(infile,5,SEEK_SET);
  len=fread(&natoms, sizeof(long),1,infile); //natoms
  fseek(infile,13,SEEK_SET);
  len=fread(&obs, sizeof(short),1,infile); //obs

  //BOX
  fseek(infile,15,SEEK_SET);
  fread(&bxx, sizeof(double),1,infile);//box_x.x
  fseek(infile,23,SEEK_SET);
  fread(&bxy, sizeof(double),1,infile);//box_x.y
  fseek(infile,31,SEEK_SET);
  fread(&bxz, sizeof(double),1,infile);//box_x.z

  fseek(infile,39,SEEK_SET);
  fread(&byx, sizeof(double),1,infile);//box_y.x
  fseek(infile,47,SEEK_SET);
  fread(&byy, sizeof(double),1,infile);//box_y.y
  fseek(infile,55,SEEK_SET);
  fread(&byz, sizeof(double),1,infile);//box_y.z
  
  fseek(infile,63,SEEK_SET);
  fread(&bzx, sizeof(double),1,infile); //box_z.x
  fseek(infile,71,SEEK_SET);
  fread(&bzy, sizeof(double),1,infile); //box_z.y
  fseek(infile,79,SEEK_SET);
  fread(&bzz, sizeof(double),1,infile); //boxzy.z
  

printf("incar:%s, disp:%d,atoms:%ld,obs:%d\n",IMD,disp,natoms,obs);
printf("bxx:%f,bxy:%f,bxz:%f\n", bxx,byx,byz);
printf("byx:%f,byy:%f,byz:%f\n", byx,byy,byz);
printf("bzx:%f,bzy:%f,bzz:%f\n", bzx,bzy,bzz);

 fseek(infile,disp,SEEK_SET);  
 double nr;
 double mass;
 double x,y,z;
 double vx,vy,vz;
 int count; 

 FILE* outfile=fopen(argv[2],"w");
 //Write header
 fprintf(outfile,"#F A 1 1 1 3 3\n");
 fprintf(outfile,"#C number type mass x y z vx vy vz\n");
 fprintf(outfile,"#X %.15e %.15e %.15e\n", bxx,bxy,bxz);
 fprintf(outfile,"#Y %.15e %.15e %.15e\n", byx,byy,byz);
 fprintf(outfile,"#Z %.15e %.15e %.15e\n", bzx,bzy,bzz);
 fprintf(outfile,"#E\n");
 
 printf("writing atoms..\n");
 while(!feof(infile)) 
 {
    fread(&nr, sizeof(double),1,infile);
    fread(&mass, sizeof(double),1,infile);
    fread(&x, sizeof(double),1,infile);
    fread(&y, sizeof(double),1,infile);
    fread(&z, sizeof(double),1,infile);
    fread(&vx, sizeof(double),1,infile);
    fread(&vy, sizeof(double),1,infile);
    fread(&vz, sizeof(double),1,infile);
    //printf("nr:%f, mass: %f, x:%f, y:%f, z:%f, vx:%f, vy:%f, vz:%f\n",nr,mass,x,y,z,vx,vy,vz);  
    fprintf(outfile,"%d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", (int) nr, 0, mass,x,y,z,vx,vy,vz);
    count++;
 }  
  
  printf("written %d atoms\n",count-1);
  fclose(infile);
  fclose(outfile);
  return 0;
}

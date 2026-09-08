// Triangle traversal, clipping order and Float32 stores match the JS reference.
__attribute__((export_name("blanket_contact"))) void blanket_contact(uint32_t ids_p,uint32_t weights_p,uint32_t indices_p,uint32_t vertices,uint32_t indices,uint32_t projected_p,uint32_t heights_p,int columns,int rows,double dx,double dz,double bedX,double bedZ){
  const uint32_t *ids=u32(ids_p),*ix=u32(indices_p);const double *w=d64(weights_p),*x=d64(x_p);double *p=d64(projected_p),*heights=d64(heights_p);
  for(uint32_t vertex=0;vertex<vertices;vertex++){
    uint32_t j=vertex*3;p[j]=p[j+1]=p[j+2]=0;
    for(int k=0;k<4;k++){uint32_t o=vertex*4+k,id=ids[o]*3;for(int a=0;a<3;a++)p[j+a]+=x[id+a]*w[o];}
    p[j]=(p[j]-bedX+.063)/dx;p[j+2]=(p[j+2]-bedZ+.026)/dz;
  }
  for(uint32_t t=0;t<indices;t+=3){
    uint32_t a=ix[t]*3,b=ix[t+1]*3,c=ix[t+2]*3;
    double ax=p[a],az=p[a+2],bx=p[b]-ax,bz=p[b+2]-az,cx=p[c]-ax,cz=p[c+2]-az,det=bx*cz-bz*cx;
    if(cabs(det)<1e-10)continue;
    double minX=dmax(0,__builtin_ceil(dmin(dmin(ax,p[b]),p[c]))),maxX=dmin(columns-1,__builtin_floor(dmax(dmax(ax,p[b]),p[c])));
    double minZ=dmax(0,__builtin_ceil(dmin(dmin(az,p[b+2]),p[c+2]))),maxZ=dmin(rows-1,__builtin_floor(dmax(dmax(az,p[b+2]),p[c+2])));
    for(int z=(int)minZ;z<=maxZ;z++)for(int x=(int)minX;x<=maxX;x++){
      double u=((x-ax)*cz-(z-az)*cx)/det,v=(bx*(z-az)-bz*(x-ax))/det;
      if(u< -1e-8||v< -1e-8||u+v>1.00000001)continue;
      double y=p[a+1]+u*(p[b+1]-p[a+1])+v*(p[c+1]-p[a+1]);int i=z*columns+x;heights[i]=dmax(heights[i],y+.00065);
    }
  }
}
static int blanket_clip(const double *input,int count,double *output,double nx,double nz,double offset){
  int inside=1,outside=1;double distances[10];
  for(int i=0;i<count;i++){int j=i*3;double d=input[j]*nx+input[j+2]*nz-offset;distances[i]=d;if(d<0)inside=0;else outside=0;}
  if(inside){for(int i=0;i<count*3;i++)output[i]=input[i];return count;}if(outside)return 0;
  int size=0;
  for(int i=0;i<count;i++){
    int j=i*3,previous=(i+count-1)%count,k=previous*3;double d=distances[i],e=distances[previous];
    if((d>=0)!=(e>=0)){double t=e/(e-d);for(int a=0;a<3;a++)output[size*3+a]=input[k+a]+t*(input[j+a]-input[k+a]);size++;}
    if(d>=0){for(int a=0;a<3;a++)output[size*3+a]=input[j+a];size++;}
  }return size;
}
__attribute__((export_name("blanket_clearance"))) int blanket_clearance(uint32_t indices_p,uint32_t indices,uint32_t projected_p,uint32_t cloth_p,int columns,int rows,double dx,double dz,double bedX,double bedZ){
  const float *positions=f32(surface_positions_p);const uint32_t *ix=u32(indices_p);double *p=d64(projected_p);float *cloth=f32(cloth_p);
  for(uint32_t j=0;j<surface_count*3;j+=3){p[j]=(positions[j]-bedX+.063)/dx;p[j+1]=positions[j+1];p[j+2]=(positions[j+2]-bedZ+.026)/dz;}
  int changed=0;double source[9],aPoly[30],bPoly[30];
  for(uint32_t t=0;t<indices;t+=3){
    double minX=1.0/0.0,maxX=-1.0/0.0,minZ=1.0/0.0,maxZ=-1.0/0.0,top=-1.0/0.0;
    for(int k=0;k<3;k++){
      uint32_t j=ix[t+k]*3;int o=k*3;double x=p[j],z=p[j+2];source[o]=x;source[o+1]=p[j+1];source[o+2]=z;
      minX=dmin(minX,x);maxX=dmax(maxX,x);minZ=dmin(minZ,z);maxZ=dmax(maxZ,z);top=dmax(top,p[j+1]);
    }
    double endZ=dmin(rows-2,__builtin_floor(maxZ)),endX=dmin(columns-2,__builtin_floor(maxX));
    for(int z=(int)dmax(0,__builtin_floor(minZ));z<=endZ;z++)for(int x=(int)dmax(0,__builtin_floor(minX));x<=endX;x++){
      int a=z*columns+x,b=a+1,c=a+columns,d=c+1;
      if(top+.00045<=dmin(dmin(dmin(cloth[a*3+1],cloth[b*3+1]),cloth[c*3+1]),cloth[d*3+1]))continue;
      int count=blanket_clip(source,3,bPoly,1,0,x);if(!count)continue;
      count=blanket_clip(bPoly,count,aPoly,-1,0,-x-1);if(!count)continue;
      count=blanket_clip(aPoly,count,bPoly,0,1,z);if(!count)continue;
      count=blanket_clip(bPoly,count,aPoly,0,-1,-z-1);if(!count)continue;
      for(int half=0;half<2;half++){
        int size=blanket_clip(aPoly,count,bPoly,half?1:-1,half?1:-1,(half?1:-1)*(x+z+1));
        int i0=(half?d:a)*3+1,i1=(half?c:b)*3+1,i2=(half?b:c)*3+1;
        for(int k=0;k<size;k++){
          double u=bPoly[k*3]-x,v=bPoly[k*3+2]-z;
          double w0=dmax(0,half?u+v-1:1-u-v),w1=dmax(0,half?1-u:u),w2=dmax(0,half?1-v:v);
          double height=(double)cloth[i0]*w0+(double)cloth[i1]*w1+(double)cloth[i2]*w2,depth=bPoly[k*3+1]+.00045-height;
          if(depth<=0)continue;
          double target=bPoly[k*3+1]+.00045001,d0=dmax(0,target-cloth[i0]),d1=dmax(0,target-cloth[i1]),d2=dmax(0,target-cloth[i2]);
          double amount=dmin(1,(depth+1e-8)/(w0*d0+w1*d1+w2*d2));
          cloth[i0]=(float)(cloth[i0]+amount*d0);cloth[i1]=(float)(cloth[i1]+amount*d1);cloth[i2]=(float)(cloth[i2]+amount*d2);changed=1;
        }
      }
    }
  }return changed;
}

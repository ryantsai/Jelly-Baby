// Included in the existing kernel: all cage arrays remain authoritative.
// Box record: center[3], axes[9], halfSize[3], motion index (-1 = static).
// Motion record: speed, pivotY, pivotZ, inertia. Shared records preserve impulse order.
static inline double cabs(double v){return __builtin_fabs(v);}
static void collision_point(uint32_t vertex,const double *x,double *p){
  const uint32_t *ids=u32(surface_ids_p)+vertex*4;const double *w=d64(surface_weights_p)+vertex*4;
  p[0]=p[1]=p[2]=0;
  for(int k=0;k<4;k++)for(int a=0;a<3;a++)p[a]+=x[ids[k]*3+a]*w[k];
}
static void collision_bounds(double *b,double magnitude,double error){
  const double *x=d64(x_p);
  for(int a=0;a<3;a++){b[a]=1.0/0.0;b[a+3]=-1.0/0.0;}
  for(uint32_t i=0;i<node_count*3;i+=3)for(int a=0;a<3;a++){b[a]=dmin(b[a],x[i+a]);b[a+3]=dmax(b[a+3],x[i+a]);}
  for(int a=0;a<3;a++){
    double center=(b[a]+b[a+3])*.5,radius=(b[a+3]-b[a])*.5*magnitude+cabs(center)*error;
    double pad=1e-10*(1+cabs(center)+radius);b[a]=center-radius-pad;b[a+3]=center+radius+pad;
  }
}
static int collision_overlaps(const double *b,const double *box,double margin){
  double d[3],r[3];for(int a=0;a<3;a++){d[a]=(b[a]+b[a+3])*.5-box[a];r[a]=(b[a+3]-b[a])*.5;}
  for(int a=0;a<3;a++){
    const double *n=box+3+a*3;
    if(cabs(d[0]*n[0]+d[1]*n[1]+d[2]*n[2])>box[12+a]+margin+r[0]*cabs(n[0])+r[1]*cabs(n[1])+r[2]*cabs(n[2]))return 0;
  }return 1;
}
static void collision_contact(uint32_t vertex,double denominator,const double *p,double nx,double ny,double nz,double depth,double *motion){
  if(denominator<1e-15)return;
  const uint32_t *ids=u32(surface_ids_p)+vertex*4;const double *w=d64(surface_weights_p)+vertex*4,*im=d64(inverse_mass_p);
  double *x=d64(x_p),*v=d64(velocity_p),normalVelocity=0;
  for(int k=0;k<4;k++){
    uint32_t id=ids[k],j=id*3;double amount=im[id]*w[k]*depth/denominator;
    x[j]+=amount*nx;x[j+1]+=amount*ny;x[j+2]+=amount*nz;
  }
  for(int k=0;k<4;k++){uint32_t j=ids[k]*3;normalVelocity+=w[k]*(v[j]*nx+v[j+1]*ny+v[j+2]*nz);}
  double mv=0,mi=0;
  if(motion){
    mv=0*nx+(motion[0]*(p[2]-motion[2]))*ny+(-motion[0]*(p[1]-motion[1]))*nz;
    double jacobian=ny*(p[2]-motion[2])-nz*(p[1]-motion[1]);mi=jacobian*jacobian/motion[3];
  }
  double relative=normalVelocity-mv;if(relative>=0)return;
  double impulse=-relative/(denominator+dmax(0,mi));
  for(int k=0;k<4;k++){
    uint32_t id=ids[k],j=id*3;double amount=im[id]*w[k]*impulse;
    v[j]+=amount*nx;v[j+1]+=amount*ny;v[j+2]+=amount*nz;
  }
  if(motion)motion[0]+=((-impulse*ny)*(p[2]-motion[2])+(-impulse*nz)*(-(p[1]-motion[1])))/motion[3];
}
static int collision_throw(const uint32_t *vertices,uint32_t samples,const double *boxes,const uint32_t *candidates,uint32_t count,double margin,double restitution){
  const double *previous=d64(previous_p),*mass=d64(mass_p);double *x=d64(x_p),*v=d64(velocity_p);
  double vx=0,vy=0,vz=0,dx=0,dy=0,dz=0;
  for(uint32_t i=0;i<node_count;i++){
    uint32_t j=i*3;double w=mass[i]/total_mass;
    vx+=w*v[j];vy+=w*v[j+1];vz+=w*v[j+2];dx+=w*(x[j]-previous[j]);dy+=w*(x[j+1]-previous[j+1]);dz+=w*(x[j+2]-previous[j+2]);
  }
  if(vx*vx+vy*vy+vz*vz<.8*.8||dx*dx+dy*dy+dz*dz<1e-10)return 0;
  double first=1,nx=0,ny=0,nz=0;
  for(uint32_t s=0;s<samples;s++){
    double p[3];collision_point(vertices[s],previous,p);
    for(uint32_t i=0;i<count;i++){
      const double *box=boxes+candidates[i]*16;if(box[15]>=0)continue;
      double enter=0,leave=first,ax=0,ay=0,az=0,px=p[0]-box[0],py=p[1]-box[1],pz=p[2]-box[2];
      for(int axis=0;axis<3;axis++){
        const double *n=box+3+axis*3;double extent=box[12+axis]+margin;
        double q=px*n[0]+py*n[1]+pz*n[2],travel=dx*n[0]+dy*n[1]+dz*n[2];
        if(cabs(travel)<1e-15){if(cabs(q)>=extent){leave=-1;break;}continue;}
        double a=(-extent-q)/travel,b=(extent-q)/travel,near=dmin(a,b),far=dmax(a,b);
        if(near>=enter){enter=near;double sign=travel>0?-1:1;ax=sign*n[0];ay=sign*n[1];az=sign*n[2];}
        leave=dmin(leave,far);if(enter>leave)break;
      }
      if(enter<=leave&&enter<first&&(ax!=0||ay!=0||az!=0)){
        double remaining=-(1-dmax(0,enter))*(dx*ax+dy*ay+dz*az);
        if(remaining>margin*.25){first=enter;nx=ax;ny=ay;nz=az;}
      }
    }
  }
  if(first==1)return 0;
  double fraction=dmax(0,first-1e-5),incoming=vx*nx+vy*ny+vz*nz,impulse=incoming<0?-(1+restitution)*incoming:0;
  for(uint32_t j=0;j<node_count*3;j+=3){
    x[j]=previous[j]+fraction*dx;x[j+1]=previous[j+1]+fraction*dy;x[j+2]=previous[j+2]+fraction*dz;
    v[j]+=impulse*nx;v[j+1]+=impulse*ny;v[j+2]+=impulse*nz;
  }return 1;
}
__attribute__((export_name("collision_boxes"))) int collision_boxes(uint32_t vertices_p,uint32_t denominators_p,uint32_t samples,double magnitude,double error,uint32_t boxes_p,uint32_t count,uint32_t candidates_p,uint32_t motions_p,double margin,double restitution,double floor){
  const uint32_t *vertices=u32(vertices_p);const double *denominators=d64(denominators_p),*boxes=d64(boxes_p);uint32_t *candidates=u32(candidates_p);
  double b[6];collision_bounds(b,magnitude,error);uint32_t found=0;
  for(uint32_t i=0;i<count;i++)if(collision_overlaps(b,boxes+i*16,margin))candidates[found++]=i;
  if(!found)return 0;
  int changed=collision_throw(vertices,samples,boxes,candidates,found,margin,restitution),exhaustive=0;
  if(!changed)for(int iteration=0;iteration<2;iteration++){
    int iterationChanged=0;
    for(uint32_t s=0;s<samples;s++){
      double p[3];collision_point(vertices[s],d64(x_p),p);
      for(uint32_t i=0;i<(exhaustive?count:found);i++){
        uint32_t index=exhaustive?i:candidates[i];const double *box=boxes+index*16;
        double dx=p[0]-box[0],dy=p[1]-box[1],dz=p[2]-box[2];
        double qx=dx*box[3]+dy*box[4]+dz*box[5],hx=box[12]+margin,hy=box[13]+margin,hz=box[14]+margin;
        if(cabs(qx)>=hx)continue;
        double qy=dx*box[6]+dy*box[7]+dz*box[8];if(cabs(qy)>=hy)continue;
        double qz=dx*box[9]+dy*box[10]+dz*box[11];if(cabs(qz)>=hz)continue;
        double depth=hx-cabs(qx),side=qx;const double *n=box+3;
        if(hy-cabs(qy)<depth){depth=hy-cabs(qy);side=qy;n=box+6;}
        if(hz-cabs(qz)<depth){depth=hz-cabs(qz);side=qz;n=box+9;}
        double sign=side<0?-1:1;
        collision_contact(vertices[s],denominators[s],p,sign*n[0],sign*n[1],sign*n[2],depth,box[15]<0?0:d64(motions_p)+(uint32_t)box[15]*4);
        changed=iterationChanged=1;collision_point(vertices[s],d64(x_p),p);
        if(!exhaustive){i=index;exhaustive=1;}
      }
    }if(!iterationChanged)break;
  }
  if(changed)stabilize_contacts(floor);return changed;
}
__attribute__((export_name("collision_cylinder"))) int collision_cylinder(uint32_t vertices_p,uint32_t denominators_p,uint32_t samples,double magnitude,double error,double cx,double cz,double radius,double minY,double maxY,double margin,double floor){
  double boundary=radius+margin,squared=boundary*boundary,b[6];collision_bounds(b,magnitude,error);
  if(b[3]<cx-boundary||b[0]>cx+boundary||b[4]<minY||b[1]>maxY||b[5]<cz-boundary||b[2]>cz+boundary)return 0;
  double dx=dmax(dmax(b[0]-cx,0),cx-b[3]),dz=dmax(dmax(b[2]-cz,0),cz-b[5]);if(dx*dx+dz*dz>squared)return 0;
  int changed=0;const uint32_t *vertices=u32(vertices_p);const double *denominators=d64(denominators_p);
  for(int iteration=0;iteration<2;iteration++){
    int iterationChanged=0;
    for(uint32_t s=0;s<samples;s++){
      double p[3];collision_point(vertices[s],d64(x_p),p);if(p[1]<=minY||p[1]>=maxY)continue;
      double x=p[0]-cx,z=p[2]-cz,d=x*x+z*z;if(d>=squared)continue;
      double distance=dsqrt(d),nx=distance<1e-9?1:x/distance,nz=distance<1e-9?0:z/distance;
      collision_contact(vertices[s],denominators[s],p,nx,0,nz,boundary-distance,0);changed=iterationChanged=1;
    }if(!iterationChanged)break;
  }
  if(changed)stabilize_contacts(floor);return changed;
}

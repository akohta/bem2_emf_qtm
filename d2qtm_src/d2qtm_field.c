#include "bem2_emf_qtm.h"

int Hz_s(double complex *Hz,double *r,int type,DOMD *md)
{
  double modx(int *l,double x,double d);
  
  double _Complex CC[7],k,A=0.0,B=0.0,ce;
  double F=0.0,lkxd,rt[2];
  int s,id,i,l;

  if(md->PD==0){
    lkxd=0.0;
    rt[0]=r[0];    rt[1]=r[1];
    ce=1.0;
  }
  else { // md->PD==1
    rt[0]=modx(&l,r[0],md->qd.d);
    rt[1]=r[1];
    lkxd=(double)l*md->qd.kx*md->qd.d;
    ce=cos(lkxd)+I*sin(lkxd);
  }

  id=domain_id(rt,md);
  k=md->wd.k0*md->n[id];

  for(s=1;s<=md->bd.sb[id].Ne;s++){
    if     (type==0){
      if(md->PD==1 && id==0) coefficient_qpgf_GL(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_GL(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else if(type==1){
      if(md->PD==1 && id==0) coefficient_qpgf_GK(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_GK(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else if(type==2){
      if(md->PD==1 && id==0) coefficient_qpgf_HP(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_HP(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else if(type==3){
      if(md->PD==1 && id==0) coefficient_qpgf_DE(CC,rt,md->bd.sb[id].eid[s],md);
      else  coefficient_DE(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else{
      if(md->PD==1 && id==0) coefficient_qpgf_NV(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_NV(CC,rt,md->bd.sb[id].eid[s],k,md);
    }

    for(i=0;i<3;i++){
      A+= md->bd.sb[id].dudn[s][i]*CC[i];
      B+=-md->bd.sb[id].u[s][i]*CC[3+i];
    }
    F+=creal(CC[6]);
  }
  if(id==0) *Hz=ce*(A+B)/(1.0+F);
  else *Hz=ce*(A+B)/F;

  return id;
}

int Hz_t(double _Complex *Hz,double *rt,int type,DOMD *md)
{
  int id;

  id=Hz_s(Hz,rt,type,md);
  if(id==0) *Hz+=infd_Hz(rt,&(md->wd));
  return id;
}

int Hz_i(double _Complex *Hz,double *rt,int type,DOMD *md)
{
  *Hz=infd_Hz(rt,&(md->wd));
  return domain_id(rt,md);
}

int Er_s(double complex *Er,double *r,int type,DOMD *md)
{
  double modx(int *l,double x,double d);
  
  double complex CC[7],k,Ax=0.0,Bx=0.0,Ay=0.0,By=0.0,ce;
  double F=0.0,lkxd,rt[2];
  int s,id,i,l;

  if(md->PD==0){
    lkxd=0.0;
    rt[0]=r[0];    rt[1]=r[1];
    ce=1.0;
  }
  else { // md->PD==1
    rt[0]=modx(&l,r[0],md->qd.d);
    rt[1]=r[1];
    lkxd=(double)l*md->qd.kx*md->qd.d;
    ce=cos(lkxd)+I*sin(lkxd);
  }

  id=domain_id(rt,md);
  k=md->wd.k0*md->n[id];

  for(s=1;s<=md->bd.sb[id].Ne;s++){
    if     (type==0){
      if(md->PD==1 && id==0) coefficient_qpgf_GL(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_GL(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else if(type==1){
      if(md->PD==1 && id==0) coefficient_qpgf_GK(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_GK(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else if(type==2){
      if(md->PD==1 && id==0) coefficient_qpgf_HP(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_HP(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else if(type==3){
      if(md->PD==1 && id==0) coefficient_qpgf_DE(CC,rt,md->bd.sb[id].eid[s],md);
      else  coefficient_DE(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else{
      if(md->PD==1 && id==0) coefficient_qpgf_NV(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_NV(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    
    for(i=0;i<3;i++){
      Ax+= md->bd.sb[id].dvdn[s][i]*CC[i];
      Bx+=-md->bd.sb[id].v[s][i]*CC[3+i];
      Ay+= md->bd.sb[id].dwdn[s][i]*CC[i];
      By+=-md->bd.sb[id].w[s][i]*CC[3+i];
    }
    F+=creal(CC[6]);
  }
  if(id==0) F+=1.0;
  Er[0]=ce*(Ax+Bx)/F;
  Er[1]=ce*(Ay+By)/F;

  return id;
}

int Er_t(double complex *Er,double *rt,int type,DOMD *md)
{
  double complex ei[2];
  int id;

  id=Er_s(Er,rt,type,md);
  if(id==0){
    infd_Er(ei,rt,&(md->wd));
    Er[0]+=ei[0];
    Er[1]+=ei[1];
  }
  return id;
}

int Er_i(double complex *Er,double *rt,int type,DOMD *md)
{
  infd_Er(Er,rt,&(md->wd));
  return domain_id(rt,md);
}

int HE_s(double complex *Hz,double complex *Er,double *r,int type,DOMD *md)
{
  double modx(int *l,double x,double d);
  
  double complex CC[7],k;
  double complex A=0.0,Ax=0.0,Ay=0.0,B=0.0,Bx=0.0,By=0.0,ce;
  double F=0.0,lkxd,rt[2];
  int s,id,i,l;

  if(md->PD==0){
    lkxd=0.0;
    rt[0]=r[0];
    rt[1]=r[1];
    ce=1.0;
  }
  else { // md->PD==1
    rt[0]=modx(&l,r[0],md->qd.d);
    rt[1]=r[1];
    lkxd=(double)l*md->qd.kx*md->qd.d;
    ce=cos(lkxd)+I*sin(lkxd);
  } 

  id=domain_id(rt,md);
  k=md->wd.k0*md->n[id];

  for(s=1;s<=md->bd.sb[id].Ne;s++){
    if     (type==0){
      if(md->PD==1 && id==0) coefficient_qpgf_GL(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_GL(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else if(type==1){
      if(md->PD==1 && id==0) coefficient_qpgf_GK(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_GK(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else if(type==2){
      if(md->PD==1 && id==0) coefficient_qpgf_HP(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_HP(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else if(type==3){
      if(md->PD==1 && id==0) coefficient_qpgf_DE(CC,rt,md->bd.sb[id].eid[s],md);
      else  coefficient_DE(CC,rt,md->bd.sb[id].eid[s],k,md);
    }
    else{
      if(md->PD==1 && id==0) coefficient_qpgf_NV(CC,rt,md->bd.sb[id].eid[s],md);
      else coefficient_NV(CC,rt,md->bd.sb[id].eid[s],k,md);
    }

    for(i=0;i<3;i++){
      A += md->bd.sb[id].dudn[s][i]*CC[i];      B +=-md->bd.sb[id].u[s][i]*CC[3+i];
      Ax+= md->bd.sb[id].dvdn[s][i]*CC[i];      Bx+=-md->bd.sb[id].v[s][i]*CC[3+i];
      Ay+= md->bd.sb[id].dwdn[s][i]*CC[i];      By+=-md->bd.sb[id].w[s][i]*CC[3+i];
    }
    F+=creal(CC[6]);
  }
  if(id==0) F+=1.0;
  *Hz=ce*(A+B)/F;
  Er[0]=ce*(Ax+Bx)/F;
  Er[1]=ce*(Ay+By)/F;

  return id;
}

int HE_t(double complex *Hz,double complex *Er,double *rt,int type,DOMD *md)
{
  double complex th,te[2];
  int id;

  id=HE_s(Hz,Er,rt,type,md);
  if(id==0){
    infd_HE(&th,te,rt,&(md->wd));
    *Hz+=th;
    Er[0]+=te[0];
    Er[1]+=te[1];
  }
  return id;
}

int HE_i(double complex *Hz,double complex *Er,double *rt,int type,DOMD *md)
{
  infd_HE(Hz,Er,rt,&(md->wd));
  return domain_id(rt,md);
}

double complex Hz_bv(int did,double eta_t,int t,DOMD *md)
{
  double complex CC[7],k;
  double complex A=0.0,B=0.0;
  double rt[2],rs[2][3],F=0.0;
  int i,s,tid,atid,sid;

  if(-1.0 < eta_t && eta_t < 1.0){
    k=md->wd.k0*md->n[did];
    tid=md->bd.sb[did].eid[t];
    atid=abs(tid);

    for(i=0;i<3;i++){
      rs[0][i]=md->bd.x[atid][i];
      rs[1][i]=md->bd.y[atid][i];
    }
    if(tid>0) rs_eta(rt,rs, eta_t);
    else      rs_eta(rt,rs,-eta_t);

    for(s=1;s<=md->bd.sb[did].Ne;s++){
      sid=md->bd.sb[did].eid[s];
      
      if(md->PD==1 && did==0){
        if(s!=t)   coefficient_qpgf_NV(CC,rt,sid,md);
        else       coefficient_qpgf_bd_eta(CC,eta_t,sid,md);
      }
      else {
        if(s!=t) coefficient_NV(CC,rt,sid,k,md);
        else coefficient_bd_eta(CC,eta_t,sid,k,md);
      }

      for(i=0;i<3;i++){
        A+= md->bd.sb[did].dudn[s][i]*CC[i];
        B+=-md->bd.sb[did].u[s][i]*CC[3+i];
      }
      F+=creal(CC[6]);
    }
    if(did==0) return (A+B)/(1.0+F);
    else return (A+B)/F;
  }
  else {
    printf("Hz_bv() region error. eta_t=%15.14e. Exit...\n",eta_t);
    exit(1);
  }
}

void HE_bv(double complex *Hz,double complex *Er,int did,double eta_t,int t,DOMD *md)
{
  double complex CC[7],k;
  double complex A=0.0,B=0.0,Ax=0.0,Bx=0.0,Ay=0.0,By=0.0;
  double rt[2],rs[2][3],F=0.0;
  int i,s,tid,atid,sid;

  if(-1.0 < eta_t && eta_t < 1.0){
    k=md->wd.k0*md->n[did];
    tid=md->bd.sb[did].eid[t];
    atid=abs(tid);

    for(i=0;i<3;i++){
      rs[0][i]=md->bd.x[atid][i];
      rs[1][i]=md->bd.y[atid][i];
    }
    if(tid>0) rs_eta(rt,rs, eta_t);
    else      rs_eta(rt,rs,-eta_t);

    for(s=1;s<=md->bd.sb[did].Ne;s++){
      sid=md->bd.sb[did].eid[s];
      
      if(md->PD==1 && did==0){
        if(s!=t)   coefficient_qpgf_NV(CC,rt,sid,md);
        else       coefficient_qpgf_bd_eta(CC,eta_t,sid,md);
      }
      else {
        if(s!=t) coefficient_NV(CC,rt,sid,k,md);
        else coefficient_bd_eta(CC,eta_t,sid,k,md);
      }

      for(i=0;i<3;i++){
        A += md->bd.sb[did].dudn[s][i]*CC[i];   B +=-md->bd.sb[did].u[s][i]*CC[3+i];
        Ax+= md->bd.sb[did].dvdn[s][i]*CC[i];   Bx+=-md->bd.sb[did].v[s][i]*CC[3+i];
        Ay+= md->bd.sb[did].dwdn[s][i]*CC[i];   By+=-md->bd.sb[did].w[s][i]*CC[3+i];
      }
      F+=creal(CC[6]);
    }
    if(did==0){
      *Hz=(A+B)/(1.0+F);
      Er[0]=(Ax+Bx)/(1.0+F);
      Er[1]=(Ay+By)/(1.0+F);
    }
    else{
      *Hz=(A+B)/F;
      Er[0]=(Ax+Bx)/F;
      Er[1]=(Ay+By)/F;
    }
  }
  else {
    printf("d2qtm_field.c, HE_bv() region error. eta_t=%15.14e. Exit...\n",eta_t);
    exit(1);
  }
}

double complex dHzdt_bv_node_ndmtd(int did,int t,int tn,DOMD *md)
{
  double complex FD_1(int order,double h,double x0,int did,int s,DOMD *md);
  
  double dx,dy,i_aJ;
  int tid,atid;

  tid=md->bd.sb[did].eid[t];
  atid=abs(tid);

  if(tid>0){
    dx=md->bd.dx[tid][tn];
    dy=md->bd.dy[tid][tn];
  }
  else {
    if(tn==0){
      dx=-md->bd.dx[atid][1];
      dy=-md->bd.dy[atid][1];
    }
    else if(tn==1){
      dx=-md->bd.dx[atid][0];
      dy=-md->bd.dy[atid][0];
    }
    else {
      dx=-md->bd.dx[atid][2];
      dy=-md->bd.dy[atid][2];
    }
  }
  i_aJ=1.0/sqrt(dx*dx+dy*dy);

  return i_aJ*FD_1(FDO,FDH,md->bd.et[tn],did,t,md);
}

/////////////////////////////////////////////////////////////////////////
double modx(int *l,double x,double d)
{
  double t,dp2,tl,mx;

  dp2=0.5*d;

  if(x==0.0){
    *l=0;
    return x;
  }
  else if(x>0.0){
    t=(x+dp2)/d;
    mx=d*modf(t,&tl);
    *l=(int)tl;
    return mx-dp2;
  }
  else if(x<0.0){
    t=(x-dp2)/d;
    mx=d*modf(t,&tl);
    *l=(int)tl;
    return mx+dp2;
  }
  else {
    printf("Error in modx()! x=%e Exit...\n",x);
    exit(1);
  }
}

double complex FD_1(int order,double h,double x0,int did,int s,DOMD *md)
{
  double complex H;

  if(order==2){
    H=0.5*Hz_bv(did,x0+h,s,md)-0.5*Hz_bv(did,x0-h,s,md);
    return H/h;
  }
  if(order==4){
    H= 2.0/ 3.0*(Hz_bv(did,x0+1.0*h,s,md)-Hz_bv(did,x0-1.0*h,s,md))
                 -1.0/12.0*(Hz_bv(did,x0+2.0*h,s,md)-Hz_bv(did,x0-2.0*h,s,md));
    return H/h;
  }
  if(order==6){
    H= 3.0/ 4.0*(Hz_bv(did,x0+1.0*h,s,md)-Hz_bv(did,x0-1.0*h,s,md))
                 -3.0/20.0*(Hz_bv(did,x0+2.0*h,s,md)-Hz_bv(did,x0-2.0*h,s,md))
                 +1.0/60.0*(Hz_bv(did,x0+3.0*h,s,md)-Hz_bv(did,x0-3.0*h,s,md));
    return H/h;
  }
  else {
    printf("d2qtm_field.c, FD_1(), order number error! order=%d. exit\n",order);
    exit(1);
  }
}



asymptStDev<-function(mdr){
  
  w=mdr$weights/mdr$n
  
  y=all.vars(as.formula(mdr$frm))[1]
  p=mdr$data[[y]]
  
  #calculate q  vector
  q1=p*w
  q0=(1-p)*w
  q=c(q0,q1)
  
  #calculate derivative
  p=q1/(q1+q0)
  f=p-mdr$residuals
  dq0=mdr$residuals*mdr$residuals-2*mdr$residuals*p/(q0+q1)
  dq1=2*p*mdr$residuals+mdr$residuals*mdr$residuals
  
  vol=asymptSDMultinomial(p=c(q0,q1), derivative=c(dq0,dq1))
  return(vol/sqrt(mdr$n))
}

asymptSDMultinomial<-function(p,derivative){
  vec = derivative
  vnsq_1  = sum(p*vec*vec)
  
  k=length(p)
  vnsq_2=0
  
  f<-function(j){
    v=vec[j]*vec
    v=v*p[j]
    v=v*p
    return(sum(v))
  }
  
  vv=sapply(c(1:k),f)
  vnsq_2=sum(vv)
  
  vnsq  = vnsq_1 - vnsq_2
  return (sqrt(vnsq))
}



asymptoticTest<-function(mdr){
  #calculate asymptotic min eps
  vol = asymptStDev(mdr)
  qt=qnorm(1-mdr$alpha,0,1)
  aps = mdr$min.distance^2 + qt*vol
  return(sqrt(aps))
}
